# -*- coding: utf-8 -*-
"""
@created: Jul 07 2022
@modified: Jul 15 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Useful functions for drawing oncoplot-like figures.
"""

import numpy as np
import pandas as pd

# matplotlib
import matplotlib.cm as cm
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

# comut
from comut import comut

# matplotlib
import matplotlib.cm as cm

def drop_duplicates_list(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


def get_ordered_alterations_and_samples(df_alt, sams, col_gen, col_alt, col_sub_id, col_sam_id, threshold_alt,
                                        complete_gene=False, contiguous_sub=False):
    """
    From a table of alteration data per sample, return a list of ordered samples for the
    columns of the oncoplot and a list or ordred alterations for the rows of the oncoplot.

    Parameters
    ----------
    df_alt: dataframe
        Must contain the columns `col_alt`, `col_gen`, `col_sub_id`, and `col_sam_id`.
    sams: list
        List of all samples. This is needed as some samples may be absent from the table `df_alt`.
    col_gen: str
        Name designating genes.
    col_alt: str
        Name designating alterations.
    col_sub_id: str
        Name designating the subject id. Useful if you want to position next to each other samples from the same
        subject. Otherwise may be identical to `col_sam_id`.
    col_sam_id: str
        Name designating the sample id.
    threshold_alt: double
        Value in [0,1]. Only alterations seen more frequently than this threshold will be retained.
    complete_gene: bool, optional
        Set to True to list contiguously all alterations from the same in case the most frequent alteration of this gene
        passes the threshold `threshold_alt`.
    contiguous_sub: bool, optional
        Set to False to prevent reordering of samples so that samples coming from the same subject are contiguous.

    Returns
    -------
    alts_ordered, sams_ordered: tuple
        Tuple of 2 ordered lists.
    """
    # order the alterations from least frequent to most frequent
    # if complete_gene is set to True, the second, third, etc. most frequent alteration of the gene are
    # positioned after the most frequent alteration.
    if not complete_gene:
        df_alt_cnt = df_alt[col_alt].value_counts()
        df_alt_cnt = df_alt_cnt.to_frame("Count").reset_index()
        df_alt_cnt = df_alt_cnt.sort_values(by=["Count", "index"])
        df_alt_cnt = df_alt_cnt.set_index("index")["Count"]
        alts_ordered = df_alt_cnt.loc[df_alt_cnt > len(sams)*threshold_alt].index.tolist()
    else:
        # select only alterations seen more frequently than alt_threshold
        df_alt_cnt = df_alt[col_alt].value_counts()
        df_alt_cnt = df_alt_cnt.to_frame("Count").reset_index()
        df_alt_cnt = df_alt_cnt.sort_values(by=["Count", "index"])
        df_alt_cnt = df_alt_cnt.set_index("index")["Count"]
        alts_ordered = df_alt_cnt.loc[df_alt_cnt > len(sams)*threshold_alt].index.tolist()

        df_alt_sub = df_alt.loc[df_alt[col_alt].isin(alts_ordered)]
        alts_orderings = df_alt_sub[col_gen].value_counts().index[::-1].tolist()

        for alt_gen in alts_orderings:
            df_alt_gen = df_alt.loc[df_alt[col_gen]==alt_gen]
            alts_ordered_gen = df_alt_gen[col_alt].value_counts().index.tolist()[::-1]
            if len(alts_ordered_gen)>1:
                index_highest = alts_ordered.index(alts_ordered_gen[-1]) - len(alts_ordered)
                alts_ordered = [x for x in alts_ordered if x not in alts_ordered_gen]
                for i, alt in enumerate(alts_ordered_gen):
                    index_insert = index_highest + i + 1
                    if index_insert==0:
                        alts_ordered += [alt]
                    else:
                        alts_ordered.insert(index_insert, alt)
                    index_highest -= 1

    # samples are ordered using stratified lexicographical ordering
    cols_sort = [col_alt, col_sam_id, col_sub_id]
    cols_sort = drop_duplicates_list(cols_sort)
    df_alt_sam = df_alt[cols_sort].drop_duplicates()
    df_alt_sam["Mutated"] = 1
    df_alt_binary = df_alt_sam.pivot(index=col_alt, columns=col_sam_id, values="Mutated")
    df_alt_binary = df_alt_binary.loc[alts_ordered[::-1],:].fillna("0")
    df_alt_binary[df_alt_binary==1] = "1"
    df_alt_string = df_alt_binary.apply(lambda x: "".join(x.tolist()), axis=0)
    sams_alts_ordered = df_alt_string.sort_values(ascending=False).index.tolist()

    # if 2 samples come form the same subject, there are place next to each other
    if contiguous_sub:
        sams_alts_orderings = df_alt_sam[col_sub_id].unique().tolist()

        for subject in sams_alts_orderings:
            df_subject = df_alt.loc[df_alt[col_sub_id]==subject]
            sams_alts_ordered_sub = sorted(df_subject[col_sam_id].unique().tolist())
            if len(sams_alts_ordered_sub)>1:
                indexes_all =[sams_alts_ordered.index(sam) for sam in sams_alts_ordered_sub]
                index_highest = min(indexes_all) - len(sams_alts_ordered)
                sams_alts_ordered = [x for x in sams_alts_ordered if x not in sams_alts_ordered_sub]
                index_highest += len(sams_alts_ordered_sub)
                for i, alt in enumerate(sams_alts_ordered_sub):
                    index_insert = index_highest + i
                    if index_insert==0:
                        sams_alts_ordered += [alt]
                    else:
                        sams_alts_ordered.insert(index_insert, alt)
                    index_highest -= 1

    # samples with 0 detected alteration are appended at the end
    sams_no_alt = sorted(list(set(sams).difference(set(sams_alts_ordered))))
    sams_ordered = sams_alts_ordered + sams_no_alt

    print("-%s alerations will be displayed" % len(alts_ordered))
    print("-%s samples will be displayed" % len(sams_ordered))

    return alts_ordered, sams_ordered


def add_clinical_benefit_plot(df_cln, mode="simple"):
    """
    Add a `Clinical_Benefit_Plot` column that will serve to order samples and that will be displayed in the plots.

    Parameters
    ----------
    df_cln: dataframe
        Must contain the columns `Best_Overall_Response`, `Clinical_Benefit`, `Objective_Response_Confirmed`.

    Returns
    -------
    df_cln: identical to input `df_cln` with an extra column.
    """
    mask_cr_a = (df_cln["Best_Overall_Response"]=="COMPLETE") & (df_cln["Objective_Response_Confirmed"]=="Confirmed")
    mask_pr_a = (df_cln["Best_Overall_Response"]=="PARTIAL") & (df_cln["Objective_Response_Confirmed"]=="Confirmed")
    mask_cr_b = (df_cln["Best_Overall_Response"]=="COMPLETE") & (df_cln["Objective_Response_Confirmed"]=="Unconfirmed")
    mask_pr_b = (df_cln["Best_Overall_Response"]=="PARTIAL") & (df_cln["Objective_Response_Confirmed"]=="Unconfirmed")
    mask_pd = (df_cln["Best_Overall_Response"]=="PROGRESSION")
    mask_sd_a = (df_cln["Best_Overall_Response"]=="STABLE") & (df_cln["Clinical_Benefit"]=="YES")
    mask_sd_b = (df_cln["Best_Overall_Response"]=="STABLE") & (df_cln["Clinical_Benefit"]=="NO")
    df_cln["Clinical_Benefit_Plot"] = "N/A"

    if mode=="simple":
        df_cln.loc[mask_cr_a, "Clinical_Benefit_Plot"] = "Yes"
        df_cln.loc[mask_pr_a, "Clinical_Benefit_Plot"] = "Yes"
        df_cln.loc[mask_cr_b, "Clinical_Benefit_Plot"] = "No"
        df_cln.loc[mask_pr_b, "Clinical_Benefit_Plot"] = "No"
        df_cln.loc[mask_pd, "Clinical_Benefit_Plot"] = "No"
        df_cln.loc[mask_sd_a, "Clinical_Benefit_Plot"] = "No"
        df_cln.loc[mask_sd_b, "Clinical_Benefit_Plot"] = "No"
    else:
        df_cln.loc[mask_cr_a, "Clinical_Benefit_Plot"] = "CR"
        df_cln.loc[mask_pr_a, "Clinical_Benefit_Plot"] = "PR"
        df_cln.loc[mask_pd, "Clinical_Benefit_Plot"] = "PD"
        df_cln.loc[mask_sd_a, "Clinical_Benefit_Plot"] = "SD w/ benefit"
        df_cln.loc[mask_sd_b, "Clinical_Benefit_Plot"] = "SD w/o benefit"


    return df_cln


def get_tables(df_cln, df_alt, col_alt, col_alt_cat, col_sub_id, col_sam_id, alts_ordered, df_alt_top=None,
               df_alt_side=None, df_cln_side=None, categorical_pairs={}, col_side_barplot="Cohort"):
    """
    Prepare a dict of tables for the oncoplot.

    Parameters
    ----------
    df_cln: dataframe
        Must contain the columns `col_sub_id` and `col_sam_id`.
    df_alt: dataframe
        Must contain the columns `col_alt`, `col_gen`, `col_sub_id`, and `col_sam_id`.
    col_alt: str
        Name designating alterations.
    col_alt_cat: str
        Name designating alteration class. Used for the coloring of cells in the oncoplot.
    col_sub_id: str
        Name designating the subject id. Useful if you want to position next to each other samples from the same
        subject. Otherwise may be identical to `col_sam_id`.
    col_sam_id: str
        Name designating the sample id.
    alts_ordered: list
        List of ordered alterations.
    df_alt_top: dataframe
        Data for the top barplot. If not specified, default to all alterations in df_alt except for "VAF alterations".
    df_alt_side: dataframe
        Data for the side barplot. If not specified, default to alterations in the plot except for "VAF alterations".
    df_cln_side: dataframe
        Data for total of the side barplot. If not specified, default to `df_cln`.
    categorical_pairs: dict
        Keys are fixed tags defined in the code like 'her2', 'ben', etc. while values are the labels to appear on the
        plot.

    Returns
    -------
    dfs_data: dict
        Dict of tables for the oncoplot.
    """

    dfs_data = {}
    cols_comut = ["sample", "category", "value"]
    df_alt_no_vaf = df_alt.loc[~df_alt[col_alt_cat].apply(lambda x: "VAF" in x)]
    df_alt_no_vaf_inplot = df_alt_no_vaf.loc[df_alt_no_vaf[col_alt].isin(alts_ordered)]

    # dataframe for indicator same patient
    df_cnt = df_cln.groupby(col_sub_id)[col_sam_id].size()
    subs_mult = df_cnt.loc[df_cnt > 1].index.tolist()
    df_ind = df_cln.loc[df_cln[col_sub_id].isin(subs_mult)][[col_sam_id, col_sub_id]]
    df_ind = df_ind.rename(columns={col_sam_id: "sample", col_sub_id: "group"})
    map_grp = {sub: i for i, sub in enumerate(subs_mult)}
    df_ind["group"] = df_ind["group"].map(map_grp)
    dfs_data["ind"] = df_ind

    # dataframe for center comut plot
    df_cen = df_alt.rename(columns={col_sam_id: "sample", col_alt: "category", col_alt_cat: "value"})
    df_cen = df_cen[cols_comut].drop_duplicates()
    dfs_data["cen"] = df_cen

    # dataframe for categorical data biopsy number
    df_bpn = df_cln.rename(columns={col_sam_id: "sample", "Cohort": "value"})
    if "bpn" in categorical_pairs:
        df_bpn["category"] = categorical_pairs["bpn"]
    else:
        df_bpn["category"] = "Cohort"
    df_bpn = df_bpn[cols_comut].drop_duplicates()
    dfs_data["bpn"] = df_bpn

    # dataframe for categorical data her2 cohort
    df_her2 = df_cln.rename(columns={col_sam_id: "sample", "HER2_Cohort": "value"})
    if "her2" in categorical_pairs:
        df_her2["category"] = categorical_pairs["her2"]
    else:
        df_her2["category"] = "HER2 Cohort"
    df_her2.loc[df_her2["value"].isnull(), "value"] = "N/A"
    df_her2 = df_her2[cols_comut].drop_duplicates()
    dfs_data["her2"] = df_her2

    # dataframe for categorical data clinical benefit
    df_ben = df_cln.rename(columns={col_sam_id: "sample", "Clinical_Benefit_Plot": "value"})
    if "ben" in categorical_pairs:
        df_ben["category"] = categorical_pairs["ben"]
    else:
        df_ben["category"] = "Clinical Benefit"
    df_ben = df_ben[cols_comut].drop_duplicates()
    dfs_data["ben"] = df_ben

    # dataframe for categorical data matched normal / unmatched
    mask_matched = ~df_cln["Sample_Id_DNA_N"].isnull()
    df_cln.loc[mask_matched, "Matched DNA"] = "Matched normal"
    df_cln.loc[~mask_matched, "Matched DNA"] = "Unmatched"
    df_mtc = df_cln.rename(columns={col_sam_id: "sample", "Matched DNA": "value"})
    if "mtc" in categorical_pairs:
        df_mtc["category"] = categorical_pairs["mtc"]
    else:
        df_mtc["category"] = "Matched DNA"
    df_mtc = df_mtc[cols_comut].drop_duplicates()
    dfs_data["mtc"] = df_mtc

    # dataframe for top barplot
    if df_alt_top is None:
        df_alt_top = df_alt_no_vaf
    df_bur = df_alt_top.groupby(col_sam_id)[col_alt_cat].value_counts().unstack(level=-1)
    df_bur = df_bur.fillna(0).astype(int).reset_index()
    df_bur = df_bur.rename(columns={col_sam_id: "sample"})
    dfs_data["bur"] = df_bur

    # dataframe for side barplot
    if df_alt_side is None:
        df_alt_side = df_alt_no_vaf_inplot[[col_sam_id, col_alt, col_side_barplot]].drop_duplicates()
    if df_cln_side is None:
        df_cln_side = df_cln
    df_cnt = df_alt_side.groupby([col_alt, col_side_barplot]).size().unstack(level=-1)
    s_tot = df_cln_side[[col_sam_id, col_side_barplot]].drop_duplicates()[col_side_barplot].value_counts()
    df_frq = 100*df_cnt/s_tot
    df_frq = df_frq.fillna(0)
    df_frq.columns = df_frq.columns.tolist()
    df_cnt.columns = df_cnt.columns.tolist()
    df_frq = df_frq.reset_index().rename(columns={col_alt: "category"})
    df_cnt = df_cnt.reset_index().rename(columns={col_alt: "category"})
    dfs_data["frq"] = df_frq
    dfs_data["cnt"] = df_cnt

    # dataframe for stars
    df_sym = df_alt_no_vaf_inplot.rename(columns={col_sam_id: "sample", col_alt: "category", "Annotated": "value"})
    df_sym = df_sym[["sample", "category", "value"]]
    df_sym = df_sym.loc[df_sym["value"]=="Yes"].copy()
    df_sym["value"] = "Annotated OncoKB/CIViC"
    dfs_data["sym"] = df_sym

    return dfs_data


def get_mappings(t_vaf_labs, name_t_vaf_inc=None, borders=[], darkgrey_frq=False, col_side_barplot="Cohort"):
    cmap_tab20 = cm.get_cmap("tab20")
    cmap_tab20b = cm.get_cmap("tab20b")
    cmap_tab10 = cm.get_cmap("tab10")
    cmap_green = cm.get_cmap("Greens")

    mappings = {}
    mappings["cen"] = {"Mutation": "#FFBC42",
                       "Indel": "#8F2D56",
                       "Amplification": "#D81159",
                       "Deletion": "#0496FF",
                       "Multiple": "#BC8034",
                       "CNA_Other": "#EDE6F2",
                       "Additional Mutation": {'facecolor':'none', 'edgecolor':'#BC8034', 'linewidth': 3},
                       "Additional Indel": {'facecolor':'none', 'edgecolor':'#362023', 'linewidth': 3},
                       "Additional Deletion": {'facecolor':'none', 'edgecolor':'#247BA0', 'linewidth': 3},
                       "Additional Amplification": {'facecolor':'none', 'edgecolor': '#B20D30', 'linewidth': 3}}

    if name_t_vaf_inc is not None:
        mappings["cen"][name_t_vaf_inc] = {'facecolor':'none', 'edgecolor':'#04724D', 'linewidth': 3}

    for i, t_vaf_lab in enumerate(t_vaf_labs):
        color_rgba = list(cmap_green(0.1 + 0.899 * i/(len(t_vaf_labs)-1)))
        color_rgba = [round(x*255) for x in color_rgba]
        color_hex = '#{:02x}{:02x}{:02x}'.format(*tuple(color_rgba))
        mappings["cen"][t_vaf_lab] = color_hex

    mappings["bpn"] = {"DAISY-T1": "#FFE381", "DAISY-T2": "#CEC288", "TCGA-BRCA": "#9F9FED",
                       "DAISY Pre-treatment": "#FFE381", "DAISY Post-treatment": "#CEC288"}

    mappings["her2"] = {"COHORT 1": "#00BBF9", "COHORT 2": "#F9D4BB", "COHORT 3": "#F15BB5",
                       "N/A": "lightgrey"}
    mappings["ben"] = {"Yes": "#60D394", "No": "#FF9B85", "CR": "#AAF683", "PR": "#60D394", "SD w/ benefit": "#FFD97D",
                       "SD w/o benefit": "#FF9B85", "PD": "#EE6055", "N/A": "lightgrey"}
    mappings["mtc"] = {"Matched normal": "#D7BCE8", "Unmatched": "#8884FF"}
    mappings["bur"] = {k:v for k,v in mappings["cen"].items() if k not in borders}

    if darkgrey_frq:
        mappings["frq"] = {"Number of samples": "darkgrey"}
    else:
        if col_side_barplot=="Cohort":
            mappings["frq"] = mappings["bpn"]
        elif col_side_barplot=="Clinical_Benefit_Plot":
            mappings["frq"] = mappings["ben"]
        else:
            raise NotImplementedError("Please add colors for side barplot of %s" % col_side_barplot)

    mappings["sym"] = {"Annotated OncoKB/CIViC": ("black", "black")}

    return mappings


def add_vaf_for_mutations(df_alt, df_cln, col_sub_id, col_sam_id, col_alt, col_alt_det, col_alt_cat, col_t_vaf,
                          mode="tex"):
    """
    If the VAF is represented in the oncoplot colors, then only one mutation may be shown for each combination of row
    and column. In this case, let us select the mutation with the highest VAF. Exception to this rule:
      - if there are 2 alterations or more and there is another paired sample,  prefer the mutation that is also seen
      in the paired sample

    Parameters
    ----------
    df_cln: dataframe
        Must contain the columns `col_sub_id` and `col_sam_id`.
    df_alt: dataframe
        Must contain the columns `col_alt`, `col_gen`, `col_sub_id`, and `col_sam_id`.
    col_sub_id: str
        Name designating the subject id. Useful if you want to position next to each other samples from the same
        subject. Otherwise may be identical to `col_sam_id`.
    col_sam_id: str
        Name designating the sample id.
    col_alt: str
        Name designating alterations.
    col_alt_det: str
        Name designating detailed alterations.
    col_alt_cat: str
        Name designating alteration class. Used for the coloring of cells in the oncoplot.
    col_t_vaf: str
        Name designating the VAF. Used for the coloring of cells in the oncoplot.
    mode: str
        Choose 'tex' or 'unicode'.

    Returns
    -------
    df_alt, t_vaf_labs: tuple dataframe, list
        Updated table with 1 additional row for each mutation/indel with the VAF category as a value for `col_alt_cat` and
        list of VAF labels.
    """
    df_cln_gby = df_cln.groupby([col_sub_id])["Biopsy_Number"].nunique()
    subjects_paired = df_cln_gby[df_cln_gby > 1].index.tolist()
    samples_paired = df_cln.loc[df_cln[col_sub_id].isin(subjects_paired)][col_sam_id].unique().tolist()
    mask_double = df_alt[[col_sam_id, col_alt, col_alt_cat]].duplicated(keep=False)
    mask_paired = df_alt[col_sam_id].isin(samples_paired)
    mask_a = mask_double & mask_paired
    mask_b = mask_double & ~mask_paired
    mask_c = ~mask_double
    df_alt_a = df_alt.loc[mask_a]
    df_alt_b = df_alt.loc[mask_b]
    df_alt_c = df_alt.loc[mask_c]

    df_alt_a_cnt = df_alt_a.groupby([col_sub_id, col_alt, col_alt_det]).size().to_frame("Count_Alt").reset_index()
    df_alt_a_uni = df_alt_a.merge(df_alt_a_cnt, how="left", on=[col_sub_id, col_alt, col_alt_det])
    df_alt_a_uni = df_alt_a_uni.sort_values(by=["Count_Alt", col_alt, col_alt_det, "t_vaf"])
    df_alt_a_uni = df_alt_a_uni.drop_duplicates(subset=[col_sam_id, col_alt], keep="first")

    df_alt_b_uni = df_alt_b.sort_values(by=[col_alt, "t_vaf"], ascending=False)
    df_alt_b_uni = df_alt_b_uni.drop_duplicates(subset=[col_sam_id, col_alt], keep="first")

    df_alt = pd.concat((df_alt_a_uni, df_alt_b_uni, df_alt_c))

    # create column col_alt_cat that will be the basis for the colors of the oncoplot
    # in order to allow a double annotation of VAF category and alteration identity for mutations,
    # add artificial rows where col_alt_cat will the rows
    t_vaf_bins = [0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1+ 1e-5]
    t_vaf_labs = []
    for i, t_vaf_bin in enumerate(t_vaf_bins[:-1]):
        if i == 0:
            if mode=="tex":
                t_vaf_labs.append("VAF $<$ %.2g" % (t_vaf_bins[i+1]))
            else:
                t_vaf_labs.append(u"VAF < %.2g" % (t_vaf_bins[i+1]))
        elif i < len(t_vaf_bins)-2:
            if mode=="tex":
                t_vaf_labs.append("%.2g $\leq$ VAF $<$ %.2g" % (t_vaf_bins[i], t_vaf_bins[i+1]))
            else:
                t_vaf_labs.append(u"%.2g ≤ VAF < %.2g" % (t_vaf_bins[i], t_vaf_bins[i+1]))
        else:
            if mode=="tex":
                t_vaf_labs.append("%.2g $\leq$ VAF $\leq$ 1" % (t_vaf_bins[i]))
            else:
                t_vaf_labs.append("%.2g ≤ VAF ≤ 1" % (t_vaf_bins[i]))

    df_alt_mut = df_alt.loc[df_alt[col_alt_cat].isin(["Mutation", "Indel"])].copy()
    df_alt_mut[col_alt_cat] = pd.cut(df_alt_mut[col_t_vaf], bins=t_vaf_bins, labels=t_vaf_labs)
    df_alt = pd.concat((df_alt_mut, df_alt), axis=0)

    return df_alt, t_vaf_labs


def convert_num_to_str(x):
    try:
        if int(x)==x:
            y = "%d" % x
        else:
            raise ValueError
    except:
        try:
            y = "%.1f" % x
            if y=="nan":
                y = x
        except:
            y = x

    return y


def draw_comut_plot(dfs_data, mappings, sams_ordered, alts_ordered, t_vaf_labs, categorical_data=[],
                    categorical_names=[], borders=None, width_left=None, width_middle=None, shadow_width_left=None,
                    height_middle=None, height_top=None, hspace=0.01, wspace=0.05, label_bar_top="Drivers",
                    label_bar_top_fontsize=12, stacked_bar_top=True, stacked_bar_side=True,
                    label_bar_side="Number of samples", label_bar_side_fontsize=12, mode_bar_side="Number",
                    sample_indicators=False, ncol_legend=1, bbox_to_anchor_legend=(1,1),
                    ignored_plots=["Alteration Type", "Same patient bis"], ignored_values=[], gap_between_groups=0.2,
                    gap_within_groups=0.05, labels_orders={}, show_bar_side=True, show_symbols=False,
                    show_legend_cohorts=True, unit_width=0.15, unit_height=0.20):
    """
    Instantiate a `CoMut` object and populate with data. Graphical parameters are adjusted based on user parameters.


    Parameters
    ----------
    dfs_data: dict
        Dict of dataframes with data for the figure.
    mappings: dict
        Dict of colors. Keys should match keys of `dfs_data`
    sams_ordered: list
        Ordered samples.
    alts_ordered: list
        Ordered alterations.
    t_vaf_labs: list
        Names of the VAF classes.
    categorical_data: list
        List of keys of `dfs_data` holding categorical data to be displayed at the top.
    categorical_names: list
        List of the labels that will be displayed for the categorical data.
    ncol_legend: int
        Number of columns for legend.


    Returns
    ----------
    comut_alt: a `CoMut` object
        Populated with data and parameters set for the rendering the figure.
    """

    comut_alt = comut.CoMut()
    comut_alt.samples = sams_ordered

    priority = t_vaf_labs[::-1] + ["Mutation", "Amplification", "Deletion", "Indel"]
    value_order = t_vaf_labs[::-1] + ["Mutation", "Amplification", "Deletion", "Indel"]

    if sample_indicators:
        # add indicators first, since they will be at the bottom
        indicator_kwargs = {'color': 'black', 'marker': 'o', 'linewidth': 1, 'markersize': 5}

        comut_alt.add_sample_indicators(data=dfs_data["ind"], name='Same patient', plot_kwargs=indicator_kwargs,
                                        xtick_show=True, xtick_fontdict={"fontsize": 8})

        comut_alt.add_categorical_data(data=dfs_data["cen"], name='Alterations', category_order=alts_ordered,
                                       borders=borders, priority=priority, value_order=value_order,
                                       mapping=mappings["cen"], xtick_show=False, xtick_fontdict={"fontsize": 8},
                                       ytick_style='italic', ytick_fontdict={"fontsize": 12})

        comut_alt.add_sample_indicators(data=dfs_data["ind"], name='Same patient bis', plot_kwargs=indicator_kwargs,
                                        xtick_show=True, xtick_fontdict={"fontsize": 8})
    else:
        comut_alt.add_categorical_data(data=dfs_data["cen"], name='Alterations', category_order=alts_ordered,
                                       borders=borders, priority=priority, value_order=value_order,
                                       mapping=mappings["cen"], xtick_show=True, xtick_fontdict={"fontsize": 8},
                                       ytick_style='italic', ytick_fontdict={"fontsize": 12})

    for cat_data, cat_name in zip(categorical_data, categorical_names):
        comut_alt.add_categorical_data(data=dfs_data[cat_data], name=cat_name,
                                       mapping=mappings[cat_data], xtick_show=False,
                                       ytick_style='normal', ytick_fontdict={"fontsize": 12})

    comut_alt.add_bar_data(data=dfs_data["bur"], name='Alteration Type',
                           mapping=mappings["bur"], ytick_fontdict={"fontsize": 12}, stacked=stacked_bar_top,
                           ylabel=label_bar_top, ylabel_fontsize=label_bar_top_fontsize, ylabel_rotation=90)

    if show_bar_side:
        if mode_bar_side=="Number":
            data = dfs_data["cnt"]
        elif mode_bar_side=="Percent":
            data = dfs_data["frq"]

        comut_alt.add_side_bar_data(data=data, paired_name='Alterations', name="Number of samples",
                                    position="left", mapping=mappings["frq"], xtick_fontdict={"fontsize": 12},
                                    stacked=stacked_bar_side, xlabel=label_bar_side,
                                    xlabel_fontsize=label_bar_side_fontsize, xlabel_rotation=0,
                                    gap_between_groups=gap_between_groups, gap_within_groups=gap_within_groups)

    if show_symbols:
        comut_alt.add_scatter_data(data=dfs_data["sym"], paired_name='Alterations', name='Annotations',
                                   mapping=mappings["sym"], marker="*", markersize=8)


    # widths
    if show_bar_side:
        if width_left is None:
            width_left = 4
        if width_middle is None:
            width_middle = unit_width*len(sams_ordered)
        if shadow_width_left is None:
            shadow_width_left = 2.5

        total_width = width_left + width_middle + shadow_width_left
        r_width_left = width_left/total_width
        r_width_middle = width_middle/total_width
        r_shadow_width_left = shadow_width_left/total_width
        widths = [r_width_left, r_width_middle]
    else:
        if width_middle is None:
            width_middle = unit_width*len(sams_ordered)
        total_width = width_middle
        r_shadow_width_left=None
        widths = None

    # heights
    if height_middle is None:
        height_middle = unit_height*len(alts_ordered)
    if height_top is None:
        height_top = max(1.75*unit_height*height_middle, 2.5)

    total_height = height_top + height_middle
    heights = {'Alteration Type': height_top}

    # render the plot
    comut_alt.plot_comut(x_padding=0.04,
                         y_padding=0.04,
                         tri_padding=0.03,
                         figsize=(total_width,total_height),
                         hspace=hspace,
                         wspace=wspace,
                         heights=heights,
                         widths=widths,
                         shadow_width_left=r_shadow_width_left)

    # configure legends
    handles_more = []
    labels_more = []
    titles_more = []

    # legend for VAF colors
    labels = t_vaf_labs
    colors = [mappings["cen"][lab] for lab in labels]
    handles = [mpatches.Patch(color=color) for color in colors]

    handles_more.append(handles)
    labels_more.append(labels)
    titles_more.append("Alteration VAF")

    # legend for cohorts
    df_frq = dfs_data["frq"]
    n_cohorts = df_frq.shape[1]-1
    if n_cohorts > 1 and show_legend_cohorts:
        labels = df_frq.columns[1:].tolist()
        colors = [mappings["frq"][lab] for lab in labels]
        handles = [mpatches.Patch(color=color) for color in colors]

        handles_more.append(handles)
        labels_more.append(labels)
        titles_more.append("Cohort")

    if show_symbols:
        handles = []
        labels = []
        for label, color in mappings["sym"].items():
            if len(label.split("&"))==1:
                handles.append(mlines.Line2D([], [], color=color[0], marker='*', linestyle='None', markersize=10))
                labels.append(label)

        handles_more.append(handles)
        labels_more.append(labels)
        titles_more.append("Annotation")

    comut_alt.add_unified_legend(ncol=ncol_legend, axis_name="Alterations", ignored_plots=ignored_plots,
                                 ignored_values=ignored_values, handles_more=handles_more, labels_more=labels_more,
                                 titles_more=titles_more, bbox_to_anchor=bbox_to_anchor_legend,
                                 labels_orders=labels_orders)

    # last adjustments

    # show spines central comut
    axis_name = "Alterations"
    for loc in ['bottom']:
        comut_alt.axes[axis_name].spines[loc].set_visible(True)

    # labels barplot
    axis_name = "Alteration Type"
    text_a = comut_alt.axes[axis_name].get_yticklabels()[0].get_text()
    text_b = comut_alt.axes[axis_name].get_yticklabels()[-1].get_text()
    text_a = "%d" % int(round(float(text_a)))
    text_b = "%d" % int(round(float(text_b)))
    comut_alt.axes[axis_name].set_yticklabels([text_a, text_b])

    # calculate the percentage of samples with that gene mutated, rounding and adding a percent sign
    if show_bar_side:
        axis_name = "Number of samples"
        df_frq = dfs_data["frq"]
        df_frq = df_frq.set_index("category").loc[alts_ordered].reset_index()
        n_bars = df_frq.shape[1]-1
        height = (1-gap_between_groups-(n_bars-1)*gap_within_groups)/n_bars

        yticks = []
        yticklabels = []

        for i,col in enumerate(df_frq.columns[1:]):
            y_range = -0.5 + gap_between_groups/2 + (i+0.5) * height + i * gap_within_groups + \
                    np.arange(0.5, len(alts_ordered))
            yticks += y_range.tolist()
            percentages = df_frq[col].apply(convert_num_to_str)
            yticklabels += percentages.tolist()

        # set location of yticks
        comut_alt.axes[axis_name].set_yticks(yticks)

        # set labels of yticks
        comut_alt.axes[axis_name].set_yticklabels(yticklabels, ha="right")

        # set label size
        if n_bars==1:
            labelsize = 12
        elif n_bars==2:
            labelsize = 8
        elif n_bars >= 3:
            labelsize = 5

        # move the ytick labels inside the bar graph
        comut_alt.axes[axis_name].tick_params(axis='y', pad=0, length=0, labelsize=labelsize)

        # Make y axis visible (by default it is not visible)
        comut_alt.axes[axis_name].get_yaxis().set_visible(True)

        # move y axis ticks to the right
        comut_alt.axes[axis_name].yaxis.tick_right()

        # x axis ticks
        if mode_bar_side=="Number":
            x_max = dfs_data["cnt"].set_index("category").max().max()
        elif mode_bar_side=="Percent":
            x_max = dfs_data["frq"].set_index("category").max().max()
        x_max_str = convert_num_to_str(x_max)
        comut_alt.axes[axis_name].set_xticks([0, x_max])
        comut_alt.axes[axis_name].set_xticklabels([0, x_max_str], fontdict={"fontsize": 12}, style="normal", rotation=0)

    # move vertical adjustment of ytick labels of burden plot
    axis_name = "Alteration Type"

    # set labels of yticks
    yticklabels = comut_alt.axes[axis_name].get_yticklabels()
    comut_alt.axes[axis_name].set_yticklabels(yticklabels, va="bottom")

    return comut_alt
