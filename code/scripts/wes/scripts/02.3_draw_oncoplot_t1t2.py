# -*- coding: utf-8 -*-
"""
@created: May 03 2022
@modified: Jun 06 2023
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Oncoplot-like figure detailing all alterations of paired tumor samples.
"""

import argparse
import numpy as np
import pandas as pd
import re
import sys

# prettypy
import matplotlib.cm as cm
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

# comut
from comut import comut

# utils
sys.path.append("../common/functions")
from utils import combine_all_alterations, combine_tcga_alterations
from utils_comut import add_clinical_benefit_plot, get_ordered_alterations_and_samples, add_vaf_for_mutations
from utils_comut import get_tables, get_mappings, draw_comut_plot

def add_cohort(df, df_sam, col_sam_id, col_sam_wes_id):
    if "Cohort" in df:
        del df["Cohort"]
    df_sam_wes = df_sam[[col_sam_wes_id, "Cohort"]].dropna(subset=[col_sam_wes_id])
    df_sam_wes = df_sam_wes.rename(columns={col_sam_wes_id: col_sam_id})
    return df.merge(df_sam_wes, how="left", on=col_sam_id)


def select_alterations_t2_vs_t1(df_alt, subs_all, col_sub_id, col_sam_id, col_gen, col_alt, col_alt_det, col_t_vaf,
                                col_inc, t_vaf_inc_thresh, cohort_t1="DAISY Pre-treatment",
                                cohort_t2="DAISY Post-treatment"):
    # select alterations in t2 not in t1
    assert df_alt["Cohort"].isnull().sum()==0
    dfs_alt_t2_vs_t1 = []
    for sub in subs_all:
        df_alt_sub = df_alt.loc[df_alt[col_sub_id]==sub].copy()
        df_alt_sub_t1 = df_alt_sub.loc[df_alt_sub["Cohort"]==cohort_t1].copy()
        df_alt_sub_t2 = df_alt_sub.loc[df_alt_sub["Cohort"]==cohort_t2].copy()

        # identify alterations in t2 but not in t1
        alt_t2_vs_t1 = list(set(df_alt_sub_t2[col_alt]).difference(set(df_alt_sub_t1[col_alt])))
        print("-INFO: %d alterations in t2 not in t1" % len(alt_t2_vs_t1))

        # identify alterations in t1 and t2 and with increased VAF in t2
        cols_row = [col_gen, col_alt_det]
        col_row = "_".join(cols_row)
        df_alt_sub[col_row] = df_alt_sub[cols_row].apply("_".join, axis=1)
        df_alt_sub_t1[col_row] = df_alt_sub_t1[cols_row].apply("_".join, axis=1)
        df_alt_sub_t2[col_row] = df_alt_sub_t2[cols_row].apply("_".join, axis=1)
        alt_t2_and_t1 = list(set(df_alt_sub_t2[col_row]).intersection(set(df_alt_sub_t1[col_row])))
        df_alt_sub_t1_t1_and_t2 = df_alt_sub_t1.loc[df_alt_sub_t1[col_row].isin(alt_t2_and_t1)].copy()
        df_alt_sub_t2_t1_and_t2 = df_alt_sub_t2.loc[df_alt_sub_t2[col_row].isin(alt_t2_and_t1)].copy()
        df_alt_sub_t1_t1_and_t2 = df_alt_sub_t1_t1_and_t2.loc[~df_alt_sub_t1_t1_and_t2[col_t_vaf].isnull()].copy()
        df_alt_sub_t2_t1_and_t2 = df_alt_sub_t2_t1_and_t2.loc[~df_alt_sub_t2_t1_and_t2[col_t_vaf].isnull()].copy()
        df_alt_sub_t1_t1_and_t2 = df_alt_sub_t1_t1_and_t2.set_index(col_row)
        df_alt_sub_t2_t1_and_t2 = df_alt_sub_t2_t1_and_t2.set_index(col_row)
        df_tvaf_diff = df_alt_sub_t2_t1_and_t2[col_t_vaf] - df_alt_sub_t1_t1_and_t2[col_t_vaf]
        alt_t1_and_t2_inc_vaf = df_tvaf_diff.loc[df_tvaf_diff > t_vaf_inc_thresh].index.tolist()
        print("-INFO: %d alterations with increased VAF in t2" % len(alt_t1_and_t2_inc_vaf))

        df_alt_sub_t2_vs_t1 = df_alt_sub.loc[df_alt_sub[col_alt].isin(alt_t2_vs_t1)].copy()
        df_alt_sub_t1_and_t2_inc_vaf = df_alt_sub_t2.loc[df_alt_sub_t2[col_row].isin(alt_t1_and_t2_inc_vaf)].copy()
        df_alt_sub_t2_vs_t1[col_inc] = "No"
        df_alt_sub_t1_and_t2_inc_vaf[col_inc] = "Yes"
        dfs_alt_t2_vs_t1.append(pd.concat((df_alt_sub_t2_vs_t1, df_alt_sub_t1_and_t2_inc_vaf)))

    return pd.concat(dfs_alt_t2_vs_t1)


def select_alterations_seen_recurrently(df_alt, df_cln, col_sam_id, col_alt_cat, col_alt, threshold_cna=3, threshold_oth=2):
    # don't base selection on alterations seen in unmatched DNA samples
    mask_matched = ~df_cln["Sample_Id_DNA_N"].isnull()
    sams_matched = df_cln[mask_matched][col_sam_id]
    df_alt_matched = df_alt.loc[df_alt[col_sam_id].isin(sams_matched)]

    df_alt_matched_cna = df_alt_matched.loc[df_alt_matched[col_alt_cat].isin(["Amplification", "Deletion"])]
    df_alt_matched_oth = df_alt_matched.loc[~df_alt_matched[col_alt_cat].isin(["Amplification", "Deletion"])]

    df_alt_cnt_cna = df_alt_matched_cna[[col_sam_id, col_alt]].drop_duplicates()[col_alt].value_counts()
    alts_keep_cna  = df_alt_cnt_cna.loc[df_alt_cnt_cna >= threshold_cna].index.tolist()
    df_alt_cnt_oth = df_alt_matched_oth[[col_sam_id, col_alt]].drop_duplicates()[col_alt].value_counts()
    alts_keep_oth  = df_alt_cnt_oth.loc[df_alt_cnt_oth >= threshold_oth].index.tolist()
    alts_keep = alts_keep_cna + alts_keep_oth

    mask_drop = ~df_alt[col_alt].isin(alts_keep)
    df_alt = df_alt.loc[~mask_drop].copy()
    print("-INFO: dropped %d/%d alterations not seen recurrently" % (sum(mask_drop), len(mask_drop)))

    return df_alt


def select_alterations_and_level(df_alt, cat_ign, alt_lvl, col_sam_id, col_alt_cat, col_alt, col_gen, col_pth, pth=None):
    # remove alterations not in coding regions
    vcs_keep = ["3'UTR", "5'UTR", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site",
                "Translation_Start_Site", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins"]
    mask_mut = df_alt[col_alt_cat].isin(["Mutation", "Indel"])
    mask_vcs = df_alt["Variant_Classification"].isin(vcs_keep)
    mask_drop = mask_mut & (~mask_vcs)
    df_alt = df_alt.loc[~mask_drop].copy()
    print("-INFO: dropped %d/%d alterations silent or outside coding regions" % (sum(mask_drop), len(mask_drop)))

    # remove alterations from selected categories if any
    mask_drop = df_alt[col_alt_cat].isin(cat_ign)
    df_alt = df_alt.loc[~mask_drop].copy()
    print("-INFO: dropped %d/%d alterations from ignore categories" % (sum(mask_drop), len(mask_drop)))

    # choose how alterations are aggregated
    if alt_lvl == "alt":
        pass
    elif alt_lvl == "gen":
        df_alt[col_alt] = df_alt[col_gen]
    elif alt_lvl.startswith("pth"):
        df_pth = pd.read_table(pth)

        # select hallmark pathway
        if alt_lvl == "pth_c6":
            df_pth = df_pth.loc[df_pth["gs_cat"]=="C6"]
            df_pth = df_pth.rename(columns={"gs_name": col_pth, "gene_symbol": col_gen})
            df_pth[col_pth] = df_pth[col_pth].apply(lambda x: " ".join(x.split("_")))
        elif alt_lvl == "pth_h":
            df_pth = df_pth.loc[df_pth["gs_cat"]=="H"]
            df_pth = df_pth.rename(columns={"gs_name": col_pth, "gene_symbol": col_gen})
            df_pth[col_pth] = df_pth[col_pth].apply(lambda x: x.replace("HALLMARK_", ""))
            df_pth[col_pth] = df_pth[col_pth].apply(lambda x: " ".join([e.title() for e in x.split("_")]))

        df_alt = df_alt.merge(df_pth[[col_pth, col_gen]].drop_duplicates(), how="left", on=col_gen)
        df_alt[col_alt] = df_alt[col_pth]
        mask_drop = df_alt[col_pth].isnull()
        df_alt = df_alt.loc[~mask_drop]
        print("-INFO: dropped %d/%d alterations not mapped to a pathway" % (sum(mask_drop), len(mask_drop)))

    return df_alt


def select_repeated_alterations(df_alt, col_sam_id, col_alt_cat, col_alt, col_inc, col_t_vaf, name_t_vaf_inc,
                                col_isin_t1t2=None):
    priority_inc = ["Yes", "No"]
    priority_isin_t2 = ["Yes", "No"]
    priority_clf = ["Mutation", "Indel", "Amplification", "Deletion"]
    priority_rep = ['Additional Mutation', 'Additional Indel', 'Additional Deletion', 'Additional Amplification']

    # if a sample has 2 mutations in gene XXXX but also 1 indel with VAF increased in XXXX, the indel
    # with VAF increased will be prioritized
    if col_inc in df_alt:
        df_alt[col_inc] = pd.Categorical(df_alt[col_inc], priority_inc[::-1])
        df_alt[col_alt_cat] = pd.Categorical(df_alt[col_alt_cat], priority_clf[::-1])
        df_alt = df_alt.sort_values(by=[col_sam_id, col_alt, col_inc, col_alt_cat, col_t_vaf], ascending=False)

    # splits alterations according to absence/presence of duplicates
    mask_dup = df_alt.duplicated(subset=[col_sam_id, col_alt], keep=False)
    df_alt_a = df_alt.loc[~mask_dup].copy()
    df_alt_b = df_alt.loc[mask_dup].copy()

    # for alterations not in duplicates and with col_inc="Yes"
    #  - keep col_alt_cat intact
    #  - duplicate rows and change col_alt_cat to name_t_vaf_inc where applicable
    if col_inc in df_alt_a:
        mask_t_vaf_inc = df_alt_a[col_inc]=="Yes"
        df_alt_a[col_alt_cat] = df_alt_a[col_alt_cat].astype(str)
        df_alt_a_1 = df_alt_a.loc[~mask_t_vaf_inc].copy()
        df_alt_a_2 = df_alt_a.loc[mask_t_vaf_inc].copy()
        df_alt_a_3 = df_alt_a_2.copy()
        df_alt_a_3[col_alt_cat] = name_t_vaf_inc
        df_alt_a = pd.concat((df_alt_a_1, df_alt_a_2, df_alt_a_3))

    # in alterations with duplicates, select the alteration that will be displayed
    if col_isin_t1t2 is not None:
        df_alt_b[col_isin_t1t2] = pd.Categorical(df_alt_b[col_isin_t1t2], priority_isin_t2[::-1])
        cols_sort_b = [col_sam_id, col_alt, col_isin_t1t2, col_alt_cat, col_t_vaf]
    else:
        cols_sort_b = [col_sam_id, col_alt, col_alt_cat, col_t_vaf]
    df_alt_b = df_alt_b.sort_values(by=cols_sort_b, ascending=False)

    mask_keep_first = ~df_alt_b.duplicated(subset=[col_sam_id, col_alt], keep="first")
    df_alt_b_1 = df_alt_b.loc[mask_keep_first].copy()
    df_alt_b_2 = df_alt_b.loc[~mask_keep_first].copy()

    # dont consider additional mutations only through increased VAF
    if col_inc in df_alt_b_2:
        df_alt_b_2[col_alt_cat] = pd.Categorical(df_alt_b_2[col_alt_cat], priority_clf)
        df_alt_b_2[col_alt_cat] = df_alt_b_2[col_alt_cat].astype(str)
        df_alt_b_2[col_alt_cat] = "Additional " + df_alt_b_2[col_alt_cat].astype(str)
        df_alt_b_2[col_alt_cat] = pd.Categorical(df_alt_b_2[col_alt_cat], priority_rep[::-1])
        df_alt_b_2 = df_alt_b_2.sort_values(by=[col_sam_id, col_alt, col_alt_cat, col_t_vaf], ascending=False)
        mask_keep_first = ~df_alt_b_2.duplicated(subset=[col_sam_id, col_alt], keep="first")
        df_alt_b_2 = df_alt_b_2.loc[mask_keep_first]
        print("-INFO: dropping %d alterations with > 2 repeats in the same alteration and sample" % len(~mask_keep_first))

        df_alt_b_1[col_alt_cat] = df_alt_b_1[col_alt_cat].astype(str)
        df_alt_b_2[col_alt_cat] = df_alt_b_2[col_alt_cat].astype(str)
        df_alt_b = pd.concat((df_alt_b_1, df_alt_b_2))
    else:
        df_alt_b = df_alt_b_1

    return pd.concat((df_alt_a, df_alt_b))


def main(args):
    # define util variable names
    col_gen = "Hugo_Symbol"
    col_pth = "Pathway_Name"
    col_alt = "Alterations"
    col_alt_cat = "Alteration_Category"
    col_alt_det = "Alteration_Detail"
    col_sub_id = "Subject_Id"
    col_sam_id = "Sample_Id"
    col_t_vaf = "t_vaf"
    col_inc = "VAF_Increased"
    font_mode = "tex"

    if font_mode=="tex":
        plt.rcParams['text.usetex'] = True
    else:
        plt.rcParams['text.usetex'] = False

    # load data for subjects with both t1 and t2 samples
    df_alt, df_cln = combine_all_alterations(alt=args.alt, cln=args.cln, cna=args.cna, mut=args.mut, col_gen=col_gen,
                                             col_alt=col_alt, col_alt_det=col_alt_det, col_alt_cat=col_alt_cat,
                                             col_sub_id=col_sub_id, col_sam_id=col_sam_id, subs=None,
                                             samples_select="t1t2", keep_alt_det=True, cna_selection="focal_high_level")

    # load data for all samples
    df_alt_all, df_cln_all = combine_all_alterations(alt=args.alt, cln=args.cln, cna=args.cna, mut=args.mut,
                                                     col_gen=col_gen, col_alt=col_alt, col_alt_det=col_alt_det,
                                                     col_alt_cat=col_alt_cat, col_sub_id=col_sub_id,
                                                     col_sam_id=col_sam_id, subs=None, samples_select="all",
                                                     keep_alt_det=True, cna_selection="focal_high_level")

    # add Cohort
    cohort_t1 = "DAISY Pre-treatment"
    cohort_t2 = "DAISY Post-treatment"

    col_sam_wes_id = "Sample_Id_WES"
    df_sam = pd.read_excel(args.sam)
    df_sam = df_sam.loc[df_sam["Used_In_WES_Analyses"]=="YES"]
    df_sam["Cohort"] = np.nan
    df_sam.loc[df_sam[col_sam_wes_id].str.endswith("T1"), "Cohort"] = cohort_t1
    df_sam.loc[df_sam[col_sam_wes_id].str.endswith("T2"), "Cohort"] = cohort_t2

    df_alt = add_cohort(df_alt, df_sam, col_sam_id, col_sam_wes_id)
    df_cln = add_cohort(df_cln, df_sam, col_sam_id, col_sam_wes_id)
    df_alt_all = add_cohort(df_alt_all, df_sam, col_sam_id, col_sam_wes_id)
    df_cln_all = add_cohort(df_cln_all, df_sam, col_sam_id, col_sam_wes_id)

    # load data for tcga BRCA samples
    df_alt_tcga, df_cln_tcga = combine_tcga_alterations(cln=args.cln_tcga, cna=args.cna_tcga, mut=args.mut_tcga,
                                                        col_gen=col_gen, col_alt=col_alt, col_alt_det=col_alt_det,
                                                        col_alt_cat=col_alt_cat, col_sub_id=col_sub_id,
                                                        col_sam_id=col_sam_id, keep_alt_det=True,
                                                        cna_selection="focal_high_level")

    if font_mode=="tex":
        name_t_vaf_inc = 'VAF increase $>$ %.2g' % args.vaf_inc
    else:
        name_t_vaf_inc = 'VAF increase > %.2g' % args.vaf_inc

    borders = ['Additional Mutation', 'Additional Indel', 'Additional Deletion', 'Additional Amplification',
               name_t_vaf_inc]

    # df_alt_top is used to retain the list of alterations per sample. this data is used for drawing the top barplot
    df_alt_top = df_alt.copy()

    # remove alterations from certain categories and decide on the level of aggregation
    df_alt_top = select_alterations_and_level(df_alt=df_alt_top, cat_ign=args.cat_ign, alt_lvl="alt",
                                              col_sam_id=col_sam_id, col_alt_cat=col_alt_cat, col_alt=col_alt,
                                              col_gen=col_gen, col_pth=col_pth, pth=args.pth)

    sams_t1 = df_cln.loc[(df_cln["Cohort"]==cohort_t1)][col_sam_id].tolist()
    sams_t2 = df_cln.loc[(df_cln["Cohort"]==cohort_t2)][col_sam_id].tolist()
    subs_t1t2 = df_cln[col_sub_id].unique().tolist()

    if args.alt_lvl == "alt" or args.agg_mod == "v1":
        # select alterations in t2 not in t1 or in t2 and t1 but with increased vaf in t2
        df_alt_t2 = select_alterations_t2_vs_t1(df_alt=df_alt, subs_all=subs_t1t2, col_sub_id=col_sub_id,
                                                col_sam_id=col_sam_id, col_gen=col_gen, col_alt=col_alt,
                                                col_alt_det=col_alt_det, col_t_vaf=col_t_vaf, col_inc=col_inc,
                                                t_vaf_inc_thresh=args.vaf_inc, cohort_t1=cohort_t1,
                                                cohort_t2=cohort_t2)

        # remove alterations from certain categories and decide on the level of aggregation
        df_alt_t2 = select_alterations_and_level(df_alt=df_alt_t2, cat_ign=args.cat_ign, alt_lvl=args.alt_lvl,
                                                 col_sam_id=col_sam_id, col_alt_cat=col_alt_cat, col_alt=col_alt,
                                                 col_gen=col_gen, col_pth=col_pth, pth=args.pth)
    elif args.agg_mod in ["v2", "v3"]:
        # remove alterations from certain categories and decide on the level of aggregation
        df_alt = select_alterations_and_level(df_alt=df_alt, cat_ign=args.cat_ign, alt_lvl=args.alt_lvl,
                                              col_sam_id=col_sam_id, col_alt_cat=col_alt_cat, col_alt=col_alt,
                                              col_gen=col_gen, col_pth=col_pth, pth=args.pth)

        # select alterations in t2 not in t1 or in t2 and t1 but with increased vaf in t2
        df_alt_t2 = select_alterations_t2_vs_t1(df_alt=df_alt, subs_all=subs_t1t2, col_sub_id=col_sub_id,
                                                col_sam_id=col_sam_id, col_gen=col_gen, col_alt=col_alt,
                                                col_alt_det=col_alt_det, col_t_vaf=col_t_vaf, col_inc=col_inc,
                                                t_vaf_inc_thresh=args.vaf_inc, cohort_t1=cohort_t1,
                                                cohort_t2=cohort_t2)


    df_alt_all = select_alterations_and_level(df_alt=df_alt_all, cat_ign=args.cat_ign, alt_lvl=args.alt_lvl,
                                              col_sam_id=col_sam_id, col_alt_cat=col_alt_cat, col_alt=col_alt,
                                              col_gen=col_gen, col_pth=col_pth, pth=args.pth)

    df_alt_tcga = select_alterations_and_level(df_alt=df_alt_tcga, cat_ign=args.cat_ign, alt_lvl=args.alt_lvl,
                                               col_sam_id=col_sam_id, col_alt_cat=col_alt_cat, col_alt=col_alt,
                                               col_gen=col_gen, col_pth=col_pth, pth=args.pth)

    # select alterations seen in at least 2 samples or 3 samples (if CNA only)
    df_alt_t2 = select_alterations_seen_recurrently(df_alt_t2, df_cln, col_sam_id, col_alt_cat, col_alt,
                                                    threshold_cna=3, threshold_oth=2)

    if args.sel_mod=="A":
        #### 
        #### METHOD A: all hits of a given alteration are shown, even hits seen in both samples
        ####

        # select alterations for the plot from pool of alterations in df_alt_t2
        df_alt_plot = df_alt_all.loc[df_alt_all[col_sam_id].isin(sams_t1+sams_t2)].copy()
        mask_alt_plot = df_alt_plot[col_alt].isin(df_alt_t2[col_alt])
        df_alt_plot = df_alt_plot.loc[mask_alt_plot].copy()

        # add VAF inc column
        col_sam_alt = "Sample_Alt"
        df_alt_plot[col_sam_alt] = df_alt_plot[[col_sam_id, col_alt, col_alt_det]].apply("/".join, axis=1)
        df_alt_t2[col_sam_alt] = df_alt_t2[[col_sam_id, col_alt, col_alt_det]].apply("/".join, axis=1)
        assert df_alt_plot[col_sam_alt].nunique()==df_alt_plot.shape[0]
        assert df_alt_t2[col_sam_alt].nunique()==df_alt_t2.shape[0]
        df_alt_plot = df_alt_plot.merge(df_alt_t2[[col_sam_alt, col_inc]], how="left", on=col_sam_alt)
        del df_alt_plot[col_sam_alt]
        df_alt_plot[col_inc] = df_alt_plot[col_inc].fillna("No")

        # in case an alteration is repeated in a sample,
        #  - select 1 alteration from a priority list, prioritizing alterations seen in both T1/T2
        #  - select repeat and encode using "Additional" tag
        #  - drop other repeats

        col_sub_alt = "Subject_Alt"
        df_alt_plot[col_sub_alt] = df_alt_plot[[col_sub_id, col_alt, col_alt_det]].apply("/".join, axis=1)
        df_alt_plot_t1 = df_alt_plot.loc[df_alt_plot[col_sam_id].isin(sams_t1)].copy()
        df_alt_plot_t2 = df_alt_plot.loc[df_alt_plot[col_sam_id].isin(sams_t2)].copy()
        sub_alt_t1t2 = list(set(df_alt_plot_t1[col_sub_alt]).intersection(set(df_alt_plot_t2[col_sub_alt])))

        col_isin_t1t2 = "Alt_In_T1_T2"
        df_alt_plot_t1[col_isin_t1t2] = "No"
        df_alt_plot_t1.loc[df_alt_plot_t1[col_sub_alt].isin(sub_alt_t1t2), col_isin_t1t2] = "Yes"
        df_alt_plot_t2[col_isin_t1t2] = "No"
        df_alt_plot_t2.loc[df_alt_plot_t2[col_sub_alt].isin(sub_alt_t1t2), col_isin_t1t2] = "Yes"

        df_alt_plot_t1 = select_repeated_alterations(df_alt=df_alt_plot_t1, col_sam_id=col_sam_id, col_alt_cat=col_alt_cat,
                                                     col_alt=col_alt, col_inc=col_inc, col_t_vaf=col_t_vaf,
                                                     name_t_vaf_inc=name_t_vaf_inc, col_isin_t1t2=col_isin_t1t2)

        df_alt_plot_t2 = select_repeated_alterations(df_alt=df_alt_plot_t2, col_sam_id=col_sam_id, col_alt_cat=col_alt_cat,
                                                     col_alt=col_alt, col_inc=col_inc, col_t_vaf=col_t_vaf,
                                                     name_t_vaf_inc=name_t_vaf_inc, col_isin_t1t2=col_isin_t1t2)

    elif args.sel_mod=="B":
        ####
        #### METHOD B: only acquired or lost hits of a given alteration are shown
        ####

        # in case an alteration is repeated in a sample,
        #  - select 1 alteration from a priority list
        #  - select repeat and encode using "Additional" tag
        #  - drop other repeats
        df_alt_plot_t2 = select_repeated_alterations(df_alt=df_alt_t2, col_sam_id=col_sam_id, col_alt_cat=col_alt_cat,
                                                     col_alt=col_alt, col_inc=col_inc, col_t_vaf=col_t_vaf,
                                                     name_t_vaf_inc=name_t_vaf_inc)

        # select alterations in t1 samples from pool of alterations in df_alt_t2
        mask_sam_t1 = df_alt["Sample_Id"].isin(sams_t1)
        mask_alt_t1 = df_alt[col_alt].isin(df_alt_plot_t2[col_alt])
        df_alt_plot_t1 = df_alt.loc[mask_sam_t1 & mask_alt_t1].copy()

        col_sub_alt = "Subject_Alt"
        df_alt_plot_t2[col_sub_alt] = df_alt_plot_t2[[col_sub_id, col_alt, col_alt_det]].apply("/".join, axis=1)
        df_alt_plot_t1[col_sub_alt] = df_alt_plot_t1[[col_sub_id, col_alt, col_alt_det]].apply("/".join, axis=1)
        mask_sub_alt_t1 = df_alt_plot_t1[col_sub_alt].isin(df_alt_plot_t2[col_sub_alt])
        col_isin_t1t2 = "Alt_In_T1_T2"
        df_alt_plot_t1[col_isin_t1t2] = np.where(mask_sub_alt_t1, "Yes", "No")
        df_alt_plot_t1 = select_repeated_alterations(df_alt=df_alt_plot_t1, col_sam_id=col_sam_id,
                                                     col_alt_cat=col_alt_cat, col_alt=col_alt, col_inc=col_inc,
                                                     col_t_vaf=col_t_vaf, name_t_vaf_inc=name_t_vaf_inc,
                                                     col_isin_t1t2=col_isin_t1t2)
    else:
        raise ValueError("Unsupported value '%s' for sel_mod. Choose 'A' or 'B'." % args.sel_mod)

    # if mode v3, remove alterations seen at least once in t1
    if args.agg_mod=="v3":
        mask_alt_t2_rmv = df_alt_plot_t2[col_alt].isin(df_alt_plot_t1[col_alt])
        df_alt_plot_t2 = df_alt_plot_t2.loc[~mask_alt_t2_rmv].copy()

    # combine alterations in T1 and T2 samples
    if args.agg_mod=="v1":
        df_alt_plot = df_alt_plot_t2
    elif args.agg_mod=="v2":
        df_alt_plot = pd.concat((df_alt_plot_t2, df_alt_plot_t1))
    elif args.agg_mod=="v3":
        df_alt_plot = df_alt_plot_t2

    # select alterations in all samples from pool of alterations in df_alt_plot
    mask_alt_all = df_alt_all[col_alt].isin(df_alt_plot[col_alt])
    df_alt_all = df_alt_all.loc[mask_alt_all].copy()
    df_alt_all = select_repeated_alterations(df_alt=df_alt_all, col_sam_id=col_sam_id, col_alt_cat=col_alt_cat,
                                             col_alt=col_alt, col_inc=col_inc, col_t_vaf=col_t_vaf,
                                             name_t_vaf_inc=name_t_vaf_inc)


    # select alterations in tcga samples from pool of alterations in df_alt_plot
    mask_alt_tcga = df_alt_tcga[col_alt].isin(df_alt_plot[col_alt])
    df_alt_tcga = df_alt_tcga.loc[mask_alt_tcga].copy()
    df_alt_tcga = select_repeated_alterations(df_alt=df_alt_tcga, col_sam_id=col_sam_id, col_alt_cat=col_alt_cat,
                                             col_alt=col_alt, col_inc=col_inc, col_t_vaf=col_t_vaf,
                                             name_t_vaf_inc=name_t_vaf_inc)

    # set values for cohort
    df_alt_tcga["Cohort"] = "TCGA-BRCA"
    df_cln_tcga["Cohort"] = "TCGA-BRCA"

    # add rows with VAF category for col_alt_cat in order to represent the VAF in the oncoplot
    df_alt_plot, t_vaf_labs = add_vaf_for_mutations(df_alt=df_alt_plot, df_cln=df_cln, col_sub_id=col_sub_id,
                                                    col_sam_id=col_sam_id, col_alt=col_alt, col_alt_det=col_alt_det,
                                                    col_alt_cat=col_alt_cat, col_t_vaf=col_t_vaf, mode=font_mode)

    # parameters for the plot 
    cols_comut = ["sample", "category", "value"]

    # add Clinical_Benefit_Plot column
    df_cln = add_clinical_benefit_plot(df_cln, mode="simple")

    # order samples based on clinical features
    # categories_benefit = ["PR", "CR", "PD", "SD w/o benefit", "SD w/ benefit", "N/A"]
    categories_benefit = ["Yes", "No", "N/A"]
    df_cln["Clinical_Benefit_Plot"] = pd.Categorical(df_cln["Clinical_Benefit_Plot"], categories_benefit)
    df_cln = df_cln.sort_values(by=["Clinical_Benefit_Plot", "HER2_Cohort"])

    alts_ordered = sorted(df_alt_plot[col_alt].unique().tolist())[::-1]
    sams_t1_ordered = df_cln.loc[df_cln["Cohort"]==cohort_t1][col_sam_id].tolist()
    sams_t2_ordered = df_cln.loc[df_cln["Cohort"]==cohort_t2][col_sam_id].tolist()
    sams_ordered = sams_t2_ordered + sams_t1_ordered

    if args.agg_mod=="v1":
        sams_ordered = sams_t2_ordered
    else:
        sams_ordered = sams_t2_ordered + sams_t1_ordered

    # subset tables
    df_cln = df_cln.loc[df_cln[col_sam_id].isin(sams_ordered)].copy()
    df_alt_plot = df_alt_plot.loc[df_alt_plot[col_sam_id].isin(sams_ordered)].copy()
    df_alt_top = df_alt_top.loc[df_alt_top[col_sam_id].isin(sams_ordered)].copy()

    # create tables for side bar
    df_alt_side = pd.concat((df_alt_all, df_alt_tcga))
    df_cln_side = pd.concat((df_cln_all, df_cln_tcga))

    # data for the plot
    categorical_data = ["bpn", "her2", "ben", "mtc"]
    categorical_names = ["Cohort", "HER2 Cohort", "ORR", "Matched DNA"]
    categorical_pairs = {d: n for d,n in zip(categorical_data, categorical_names)}
    dfs_data = get_tables(df_cln=df_cln, df_alt=df_alt_plot, col_alt=col_alt, col_alt_cat=col_alt_cat, col_sub_id=col_sub_id,
                          col_sam_id=col_sam_id, alts_ordered=alts_ordered, df_alt_top=df_alt_top,
                          df_alt_side=df_alt_side, df_cln_side=df_cln_side, categorical_pairs=categorical_pairs)
    mappings = get_mappings(t_vaf_labs=t_vaf_labs, name_t_vaf_inc=name_t_vaf_inc, borders=borders)


    # compute widths and heights - calibrated on full cohort with 108 samples.
    width_left = 2
    width_middle = 0.25*len(sams_ordered)

    if args.alt_lvl == "alt":
        shadow_width_left = 1.6
    elif args.alt_lvl == "gen":
        shadow_width_left = 1
    elif args.alt_lvl.startswith("pth"):
        shadow_width_left = 2.4

    if args.alt_lvl == "pth_c6" and args.agg_mod=="v1":
        height_middle = 0.1*len(alts_ordered)
        height_top = max(0.3*height_middle, 2)
    else:
        height_middle = 0.15*len(alts_ordered)
        height_top = max(0.6*height_middle, 4)

    ignored_plots=["Alteration Type", "Same patient bis", "Cohort"]
    ignored_values=["Absent"] + t_vaf_labs
    labels_orders = {categorical_pairs["ben"]: categories_benefit}

    # draw figure
    comut_alt = draw_comut_plot(dfs_data=dfs_data, mappings=mappings, sams_ordered=sams_ordered,
                                alts_ordered=alts_ordered, t_vaf_labs=t_vaf_labs, categorical_data=categorical_data,
                                categorical_names=categorical_names, borders=borders,
                                label_bar_top="Nonsyn.\nAlterations", label_bar_side="Fraction of samples",
                                mode_bar_side="Percent", width_left=width_left, width_middle=width_middle,
                                shadow_width_left=shadow_width_left, height_middle=height_middle, height_top=height_top,
                                stacked_bar_side=False, ncol_legend=2, bbox_to_anchor_legend=(1, 1),
                                ignored_plots=ignored_plots, ignored_values=ignored_values, labels_orders=labels_orders)

    # save figure
    fig = comut_alt.figure
    fig.savefig(args.output, dpi=300, bbox_inches='tight')
    print("-plot saved at %s" % args.output)


# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Oncoplot-like figure detailing treatment resistances.")
    parser.add_argument("--cln", type=str, help="Path to input curated clincal table.")
    parser.add_argument('--alt', type=str, help='Path to table of alterations.')
    parser.add_argument('--mut', type=str, help='Path to table of all mutations.')
    parser.add_argument('--cna', type=str, help='Path to table of all CNAs.')
    parser.add_argument("--sam", type=str, help="Paths to table of all samples.")
    parser.add_argument('--cln_tcga', type=str, help='Path to table of tcga clinical table')
    parser.add_argument('--mut_tcga', type=str, help='Path to table of tcga mutations.')
    parser.add_argument('--cna_tcga', type=str, help='Path to table of tcga CNAs.')
    parser.add_argument('--pth', type=str, help='Path to table of all pathways.')
    parser.add_argument('--cat_ign', type=str, nargs="*", help='List of alterations categories to be ignored',
                        default=[])
    parser.add_argument('--alt_lvl', type=str, help='Level to which alterations are aggregated.',
                        default="gen")
    parser.add_argument('--agg_mod', type=str, help='Mode of aggregation if alt_lvl is not "alt"',
                        default="v2")
    parser.add_argument("--sel_mod", type=str,
                        help="Choose A for showing all hits, B for showing only differential hits.",
                        default="B")
    parser.add_argument("--vaf_inc", type=float,
                        help="Mutations with a VAF increase greater than this value will be treated as if acquired.",
                        default=0.2)
    parser.add_argument('--output', type=str,  help='Paths to output oncoplot-like.',
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
