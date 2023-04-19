# -*- coding: utf-8 -*-
"""
@created: May 02 2022
@modified: Feb 27 2023
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Oncoplot-like figure detailing treatment resistances
"""

import argparse
import numpy as np
import pandas as pd
import re
import sys

import matplotlib.pyplot as plt

# utils
sys.path.append("../common/functions")
from utils import combine_all_alterations
from utils_comut import add_clinical_benefit_plot, get_ordered_alterations_and_samples, add_vaf_for_mutations
from utils_comut import get_tables, get_mappings, draw_comut_plot

from scipy.stats import mannwhitneyu, boschloo_exact
from statsmodels.sandbox.stats.multicomp import multipletests

def main(args):
    # define util variable names
    col_lvl = "Sen_Level_Simple"
    col_gen = "Hugo_Symbol_%s" % col_lvl
    col_alt = "Alteration_%s" % col_lvl
    col_alt_cat = "Alteration_Category"
    col_alt_det = "Alteration"
    col_sub_id = "Subject_Id"
    col_sam_id = "Sample_Id"
    col_t_vaf = "t_vaf"
    font_mode = "tex"

    if font_mode=="tex":
        plt.rcParams['text.usetex'] = True
    else:
        plt.rcParams['text.usetex'] = False

    # load data
    df_alt, df_cln = combine_all_alterations(alt=args.alt, cln=args.cln, cna=None, mut=None, col_gen=col_gen,
                                             col_alt=col_alt, col_alt_det=col_alt_det, col_alt_cat=col_alt_cat,
                                             col_lvl=col_lvl, col_sub_id=col_sub_id, col_sam_id=col_sam_id, subs=None,
                                             samples_select=args.samples_select, keep_alt_det=False)

    # special cases
    # for TP53, use distinction from Res_Level_Simple
    col_lvl_oth = "Res_Level_Simple"
    col_gen_oth = "Hugo_Symbol_%s" % col_lvl_oth
    col_alt_oth = "Alteration_%s" % col_lvl_oth
    mask = df_alt[col_gen] == "TP53"
    df_alt.loc[mask, col_alt] = df_alt.loc[mask, [col_gen_oth, col_alt_oth]].apply(" ".join, axis=1)

    # parameters for the plot 
    sams = df_cln[col_sam_id].tolist()
    cols_comut = ["sample", "category", "value"]
    alts_ordered, sams_ordered = get_ordered_alterations_and_samples(df_alt=df_alt, sams=sams,
                                                                     col_gen=col_gen, col_alt=col_alt,
                                                                     col_sub_id=col_sub_id, col_sam_id=col_sam_id,
                                                                     threshold_alt=args.threshold_alt,
                                                                     complete_gene=True, contiguous_sub=True)

    # add Clinical_Benefit_Plot column
    df_cln = add_clinical_benefit_plot(df_cln, mode="simple")

    # order samples based on clinical features
    categories_benefit = ["Yes", "No"]
    df_cln["Clinical_Benefit_Plot"] = pd.Categorical(df_cln["Clinical_Benefit_Plot"], categories_benefit)
    df_sam_order = pd.DataFrame({col_sam_id: sams_ordered, "Order": np.arange(0, len(sams_ordered))})
    df_cln = df_cln.merge(df_sam_order, how="left", on=col_sam_id)
    df_cln = df_cln.sort_values(by=["Clinical_Benefit_Plot", "HER2_Cohort", "Order"])
    sams_ordered = df_cln[col_sam_id].tolist()
    df_alt = df_alt.merge(df_cln[["Subject_Id", "Clinical_Benefit_Plot"]], how="left")

    # add rows with VAF category for col_alt_cat in order to represent the VAF in the oncoplot
    df_alt, t_vaf_labs = add_vaf_for_mutations(df_alt=df_alt, df_cln=df_cln, col_sub_id=col_sub_id,
                                               col_sam_id=col_sam_id, col_alt=col_alt, col_alt_det=col_alt_det,
                                               col_alt_cat=col_alt_cat, col_t_vaf=col_t_vaf)

    # data for the plot
    if args.samples_select in ["all", "t1t2"]:
        col_side_barplot = "Cohort"
    else:
        col_side_barplot = "Clinical_Benefit_Plot"
    categorical_data = ["her2", "ben", "mtc"]
    categorical_names = ["HER2 Cohort", "ORR", "Matched DNA"]
    categorical_pairs = {d: n for d,n in zip(categorical_data, categorical_names)}
    dfs_data = get_tables(df_cln=df_cln, df_alt=df_alt, col_alt=col_alt, col_alt_cat=col_alt_cat, col_sub_id=col_sub_id,
                          col_sam_id=col_sam_id, alts_ordered=alts_ordered, categorical_pairs=categorical_pairs,
                          col_side_barplot=col_side_barplot)
    mappings = get_mappings(t_vaf_labs=t_vaf_labs, col_side_barplot=col_side_barplot)

    height_middle = 0.15*len(alts_ordered)
    height_top = 8
    ignored_plots = ["Alteration Type", "Same patient bis"]
    ignored_values=["Absent"] + t_vaf_labs
    labels_orders = {categorical_pairs["ben"]: categories_benefit}

    # draw figure
    if args.samples_select in ["all", "t1t2"]:
        sample_indicators = True
        stacked_bar_side = True
        show_legend_cohorts = True
    else:
        sample_indicators = False
        stacked_bar_side = False
        show_legend_cohorts = False

    comut_alt = draw_comut_plot(dfs_data=dfs_data, mappings=mappings, sams_ordered=sams_ordered,
                                alts_ordered=alts_ordered, t_vaf_labs=t_vaf_labs, categorical_data=categorical_data,
                                categorical_names=categorical_names, label_bar_top="Drivers", shadow_width_left=2.2,
                                mode_bar_side="Percent", label_bar_side="Fraction of samples",
                                sample_indicators=sample_indicators, height_middle=height_middle, height_top=height_top,
                                stacked_bar_side=stacked_bar_side, labels_orders=labels_orders,
                                ignored_values=ignored_values, label_bar_top_fontsize=16, label_bar_side_fontsize=16,
                                show_legend_cohorts=show_legend_cohorts, unit_width=0.25, unit_height=0.20)

    # save figure
    fig = comut_alt.figure
    fig.savefig(args.output, dpi=300, bbox_inches='tight')
    print("-plot saved at %s" % args.output)

    # for t1, run fisher-boschloo tests
    if args.samples_select=="t1":

        alts_test, _ = get_ordered_alterations_and_samples(df_alt=df_alt, sams=sams,
                                                           col_gen=col_gen, col_alt=col_alt,
                                                           col_sub_id=col_sub_id, col_sam_id=col_sam_id,
                                                           threshold_alt=0.05,
                                                           complete_gene=False, contiguous_sub=False)

        # some alterations are perfectly correlated, treat them as one event
        alts_groups = {"CCDN1 Amp": ["CCND1 Amp", "FGF4 Amp", "FGF3 Amp", "FGF19 Amp"],
                       "STK11 Del": ["STK11 Del", "TCF3 Del"]}

        df_cnt = dfs_data["cnt"].rename(columns={"category": "Alteration"})
        margins = df_cln["Clinical_Benefit_Plot"].value_counts()
        df_tests = pd.DataFrame(columns=["Alteration", "P_Value"])

        for alt in alts_test[::-1]:
            alts_group = [alt]
            for alts_group_search in alts_groups.values():
                if alt in alts_group_search:
                    alts_group = alts_group_search
                    break
            if alt==alts_group[0]:
                cnt_yes = df_cnt.loc[df_cnt["Alteration"]==alt, "Yes"].iloc[0]
                cnt_no = df_cnt.loc[df_cnt["Alteration"]==alt, "No"].iloc[0]

                table = np.array([[cnt_yes, cnt_no],
                                    [margins["Yes"]-cnt_yes, margins["No"]-cnt_no]])
                out_test = boschloo_exact(table=table, alternative="two-sided", n=100)
                pvalue = out_test.pvalue

                df_tests = df_tests.append({"Alteration": alt, "P_Value": pvalue}, ignore_index=True)

        # adjust pvalues
        df_tests["P_Value_Adj"] = multipletests(pvals=df_tests["P_Value"], alpha=0.1, method="fdr_bh")[1]
        df_tests.to_excel(args.output_tests, index=False)

# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Oncoplot-like figure detailing treatment resistances.")
    parser.add_argument("--cln", type=str, help="Path to input curated clincal table.")
    parser.add_argument('--alt', type=str, help='Path to table of alterations.')
    parser.add_argument("--sam", type=str, help="Paths to table of all samples.")
    parser.add_argument('--threshold_alt', type=float, default=0.03,
                        help='Only alterations with a frequency of at least this value will be displayed')
    parser.add_argument('--samples_select', type=str, default="t1", help='Choose "all", "t1" or "t1t2".')
    parser.add_argument('--output', type=str,  help='Paths to output oncoplot-like.')
    parser.add_argument('--output_tests', type=str,
                        help='Paths to output tables of tests of associations with clinical response.')
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
