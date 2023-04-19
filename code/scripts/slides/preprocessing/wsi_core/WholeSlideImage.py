import re
import os
import multiprocessing as mp
import cv2
import numpy as np

import openslide
from PIL import Image
import math

from .wsi_utils import save_hdf5, local_average, compute_law_feats, filter_ROI, thresh_ROI, \
    floodfill_ROI, contour_ROI, remove_black_areas, isWhitePatch_S, isBlackPatch_S

from .util_classes import isInContourV1, isInContourV2, isInContourV3_Easy, isInContourV3_Hard, \
    Contour_Checking_fn
from ..utils import load_pkl, save_pkl
from tifffile import imread

Image.MAX_IMAGE_PIXELS = 225705984000


class WholeSlideImage(object):
    def __init__(self, path, ROIsdf, tumorareasfolder):

        """
        Args:
            path (str): fullpath to WSI file
            ROIsdf (DataFrame): df with ROI coordinates for each slide
        """
        self.path = path
        self.name = re.search(r'(?:[^\\\/](?!(\\|\/)))+(?=.mrxs)', path).group(0)
        self.wsi = openslide.open_slide(path)
        self.level_downsamples = self.wsi.level_downsamples
        self.level_dim = self.wsi.level_dimensions
        self.ROIs = self.get_ROIs(ROIsdf)
        self.tumor_areas = self.get_tumor_areas(tumorareasfolder)

        self.contours_tissue = None
        self.contours_tumor = None
        self.hdf5_file = None

    @staticmethod
    def isInContours(cont_check_fn, pt):
        if cont_check_fn(pt):
            return 1
        else:
            return 0

    def getOpenSlide(self):
        return self.wsi

    def get_ROIs(self, df):
        rows = df[df.slideid == self.name]
        return rows.values

    def get_tumor_areas(self, tumor_areas_folder):
        filename = f"{self.name}-labels.ome.tif"
        if filename in os.listdir(tumor_areas_folder):
            filepath = os.path.join(tumor_areas_folder, filename)
        else:
            filepath = ""
        return filepath

    def create_tumor_mask(self, slide, seg_level):
        y_nonzero, x_nonzero, _ = np.nonzero(slide)
        yminp, ymaxp, xminp, xmaxp = np.min(y_nonzero), np.max(y_nonzero), np.min(x_nonzero), np.max(x_nonzero)
        tumor_img = Image.fromarray(imread(self.tumor_areas, series=0, level=0)).resize((xmaxp - xminp, ymaxp - yminp))
        mask = np.zeros((self.level_dim[seg_level][1], self.level_dim[seg_level][0]), dtype=np.uint8)
        mask[yminp:ymaxp, xminp:xmaxp] = np.asarray(tumor_img)
        mask = Image.fromarray(mask)
        return mask

    def load_ROI(self, seg_level, roi_idx):
        check_mask = True
        slide = self.wsi.read_region((0, 0), seg_level, self.level_dim[seg_level]).convert("RGB")
        coeff = self.level_downsamples[seg_level]
        l, t, r, b = self.ROIs[roi_idx][1:]
        slide_roi = slide.crop((l // coeff, t // coeff, r // coeff, b // coeff))
        if self.tumor_areas:
            mask = self.create_tumor_mask(slide, seg_level)
            mask = mask.crop((l // coeff, t // coeff, r // coeff, b // coeff))
            if np.asarray(mask).sum() == 0:
                check_mask = False
            tumor = np.zeros(np.asarray(slide_roi).shape)
            tumor += (np.repeat(np.expand_dims(np.asarray(mask) == 1, -1), 3, -1)) * 255
            slide_roi = cv2.bitwise_and(np.asarray(slide_roi), tumor.astype(np.uint8))
        slide_roi = remove_black_areas(slide_roi)
        return slide_roi, check_mask

    def initSegmentation(self, mask_file):
        # load segmentation results from pickle file
        import pickle
        asset_dict = load_pkl(mask_file)
        self.contours_tissue = asset_dict['tissue']

    def saveSegmentation(self, mask_file):
        # save segmentation results using pickle
        asset_dict = {'tissue': self.contours_tissue}
        save_pkl(mask_file, asset_dict)

    def rebase_contours(self, contours, seg_level, roi_idx, area_min):
        left, top = self.ROIs[roi_idx][1:3]
        new_contours = []
        for contour in contours:
            c = contour.copy()
            area = cv2.contourArea(c)
            if area > area_min:
                c[:, :, 0] = c[:, :, 0] + left / self.level_downsamples[seg_level]
                c[:, :, 1] = c[:, :, 1] + top / self.level_downsamples[seg_level]
                new_contours.append(c)
        return new_contours

    def segmentTissue(self, seg_level=0, window_avg=3, window_eng=3, thresh=30, keep_grayscale=False, inv=False,
                      start=(10, 10), area_min=2e5):
        """
            Remove artefacts and create contours via filtering, thresholding and flooding
        """

        contours_slide = []

        for roi_idx in range(len(self.ROIs)):
            img, check_mask = self.load_ROI(seg_level, roi_idx)
            if not check_mask:
                print(f"Mask empty for ROI n°{roi_idx}, continuing")
                continue
            img_avg = local_average(np.asarray(img), window_avg, keep_grayscale)
            law_feats = compute_law_feats(img_avg, window_eng)
            filterred_roi = filter_ROI(law_feats[:, :, 3].astype(np.uint8()))
            threshed_roi = thresh_ROI(filterred_roi, thresh, inv)
            flooded_roi = floodfill_ROI(threshed_roi, start)
            contours = contour_ROI(flooded_roi)
            contours = self.rebase_contours(contours, seg_level, roi_idx, area_min)
            for contour in contours:
                contours_slide.append(contour)
        scale = self.level_downsamples[seg_level]
        if len(contours_slide) == 0:
            self.contours_tissue = []
            print(f"No contours found for slide {self.name}")
            return
        self.contours_tissue = self.scaleContourDim(contours_slide, scale)
        return

    def visWSI(self, vis_level=0, color=(255, 0, 0), line_thickness=250, max_size=None, top_left=None, bot_right=None,
               custom_downsample=1, view_slide_only=False,
               number_contours=False, seg_display=True, annot_display=True):

        downsample = self.level_downsamples[vis_level]
        scale = [1 / downsample, 1 / downsample]

        if top_left is not None and bot_right is not None:
            top_left = tuple(top_left)
            bot_right = tuple(bot_right)
            w, h = tuple((np.array(bot_right) * scale).astype(int) - (np.array(top_left) * scale).astype(int))
            region_size = (w, h)
        else:
            top_left = (0, 0)
            region_size = self.level_dim[vis_level]

        img = np.array(self.wsi.read_region(top_left, vis_level, region_size).convert("RGB"))

        if not view_slide_only:
            offset = tuple(-(np.array(top_left) * scale).astype(int))
            line_thickness = int(line_thickness * math.sqrt(scale[0] * scale[1]))
            if self.contours_tissue is not None and seg_display:
                if not number_contours:
                    cv2.drawContours(img, self.scaleContourDim(self.contours_tissue, scale),
                                     -1, color, line_thickness, lineType=cv2.LINE_8, offset=offset)

                else:  # add numbering to each contour
                    for idx, cont in enumerate(self.contours_tissue):
                        contour = np.array(self.scaleContourDim(cont, scale))
                        M = cv2.moments(contour)
                        cX = int(M["m10"] / (M["m00"] + 1e-9))
                        cY = int(M["m01"] / (M["m00"] + 1e-9))
                        # draw the contour and put text next to center
                        cv2.drawContours(img, [contour], -1, color, line_thickness, lineType=cv2.LINE_8, offset=offset)
                        cv2.putText(img, "{}".format(idx), (cX, cY),
                                    cv2.FONT_HERSHEY_SIMPLEX, 2, (255, 0, 0), 10)

        img = Image.fromarray(img)

        w, h = img.size
        if custom_downsample > 1:
            img = img.resize((int(w / custom_downsample), int(h / custom_downsample)))

        if max_size is not None and (w > max_size or h > max_size):
            resizeFactor = max_size / w if w > h else max_size / h
            img = img.resize((int(w * resizeFactor), int(h * resizeFactor)))

        return img

    @staticmethod
    def scaleContourDim(contours, scale):
        return [np.array(cont * scale, dtype='int32') for cont in contours]

    def process_contours(self, save_path, patch_level=0, patch_size=256, step_size=256, **kwargs):
        save_path_hdf5 = os.path.join(save_path, str(self.name) + '.h5')
        print("Creating patches for: ", self.name, "...", )
        n_contours = len(self.contours_tissue)
        print("Total number of contours to process: ", n_contours)
        fp_chunk_size = math.ceil(n_contours * 0.05)
        init = True
        for idx, cont in enumerate(self.contours_tissue):
            if (idx + 1) % fp_chunk_size == fp_chunk_size:
                print('Processing contour {}/{}'.format(idx, n_contours))

            asset_dict, attr_dict = self.process_contour(cont, patch_level, save_path, patch_size, step_size, **kwargs)
            if len(asset_dict) > 0:
                if init:
                    save_hdf5(save_path_hdf5, asset_dict, attr_dict, mode='w')
                    init = False
                else:
                    save_hdf5(save_path_hdf5, asset_dict, mode='a')

        return self.hdf5_file

    def process_contour(self, cont, patch_level, save_path, patch_size=256, step_size=256,
                        contour_fn='four_pt', use_padding=True, top_left=None, bot_right=None):
        start_x, start_y, w, h = cv2.boundingRect(cont) if cont is not None else (
            0, 0, self.level_dim[patch_level][0], self.level_dim[patch_level][1])

        patch_downsample = (int(self.level_downsamples[patch_level]), int(self.level_downsamples[patch_level]))
        ref_patch_size = (patch_size * patch_downsample[0], patch_size * patch_downsample[0])

        img_w, img_h = self.level_dim[0]
        if use_padding:
            stop_y = start_y + h
            stop_x = start_x + w
        else:
            stop_y = min(start_y + h, img_h - ref_patch_size[1] + 1)
            stop_x = min(start_x + w, img_w - ref_patch_size[0] + 1)

        print("Bounding Box:", start_x, start_y, w, h)
        print("Contour Area:", cv2.contourArea(cont))

        if bot_right is not None:
            stop_y = min(bot_right[1], stop_y)
            stop_x = min(bot_right[0], stop_x)
        if top_left is not None:
            start_y = max(top_left[1], start_y)
            start_x = max(top_left[0], start_x)

        if bot_right is not None or top_left is not None:
            w, h = stop_x - start_x, stop_y - start_y
            if w <= 0 or h <= 0:
                print("Contour is not in specified ROI, skip")
                return {}, {}
            else:
                print("Adjusted Bounding Box:", start_x, start_y, w, h)

        if isinstance(contour_fn, str):
            if contour_fn == 'four_pt':
                cont_check_fn = isInContourV3_Easy(contour=cont, patch_size=ref_patch_size[0], center_shift=0.5)
            elif contour_fn == 'four_pt_hard':
                cont_check_fn = isInContourV3_Hard(contour=cont, patch_size=ref_patch_size[0], center_shift=0.5)
            elif contour_fn == 'center':
                cont_check_fn = isInContourV2(contour=cont, patch_size=ref_patch_size[0])
            elif contour_fn == 'basic':
                cont_check_fn = isInContourV1(contour=cont)
            else:
                raise NotImplementedError
        else:
            assert isinstance(contour_fn, Contour_Checking_fn)
            cont_check_fn = contour_fn

        step_size_x = step_size * patch_downsample[0]
        step_size_y = step_size * patch_downsample[1]

        x_range = np.arange(start_x, stop_x, step=step_size_x)
        y_range = np.arange(start_y, stop_y, step=step_size_y)
        x_coords, y_coords = np.meshgrid(x_range, y_range, indexing='ij')
        coord_candidates = np.array([x_coords.flatten(), y_coords.flatten()]).transpose()

        num_workers = mp.cpu_count()
        if num_workers > 4:
            num_workers = 4
        pool = mp.Pool(num_workers)

        iterable = [(coord, cont_check_fn, self.path, patch_level, patch_size) for coord in
                    coord_candidates]
        results = pool.starmap(WholeSlideImage.process_coord_candidate, iterable)
        pool.close()
        results = np.array([result for result in results if result is not None])

        print('Extracted {} coordinates'.format(len(results)))

        if len(results) > 1:
            asset_dict = {'coords': results}

            attr = {'patch_size': patch_size,  # To be considered...
                    'patch_level': patch_level,
                    'downsample': self.level_downsamples[patch_level],
                    'downsampled_level_dim': tuple(np.array(self.level_dim[patch_level])),
                    'level_dim': self.level_dim[patch_level],
                    'name': self.name,
                    'save_path': save_path}

            attr_dict = {'coords': attr}
            return asset_dict, attr_dict

        else:
            return {}, {}

    @staticmethod
    def is_black_white(patch):
        # white :0,5 ou 0.8
        return isBlackPatch_S(patch, rgbThresh=2, percentage=0.7) or isWhitePatch_S(patch, rgbThresh=250,
                                                                                    percentage=0.8)

    @staticmethod
    def process_coord_candidate(coord, cont_check_fn, path, patch_level, patch_size):
        if WholeSlideImage.isInContours(cont_check_fn, coord):
            wsi = openslide.open_slide(path)
            patch = wsi.read_region(coord, patch_level, tuple([patch_size, patch_size])).convert('RGB')
            if not WholeSlideImage.is_black_white(patch):
                return coord
            else:
                return None
        else:
            return None
