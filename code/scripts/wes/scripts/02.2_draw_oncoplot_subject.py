# -*- coding: utf-8 -*-
"""
@created: May 02 2022
@modified: Jan 06 2023
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Oncoplot-like figure detailing all alterations of a pair of samples.
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
import palettable
from comut import comut

# utils
sys.path.append("../common/functions")
from utils import combine_all_alterations
from utils_comut import add_clinical_benefit_plot, add_vaf_for_mutations
from utils_comut import get_tables, get_mappings, draw_comut_plot

def main(args):
    # define util variable names
    col_lvl = "Sen_Level_Simple"
    col_gen = "Hugo_Symbol"
    col_alt = "Alteration"
    col_alt_det = "Alteration_Detail"
    col_alt_cat = "Alteration_Category"
    col_sub_id = "Subject_Id"
    col_sam_id = "Sample_Id"
    col_t_vaf = "t_vaf"
    plt.rcParams['text.usetex'] = False

    df_cln = pd.read_table(args.cln)
    df_cln = df_cln.loc[df_cln["Keep"]==1]
    df_cln_cnt = df_cln.groupby(col_sub_id)["Biopsy_Number"].nunique()

    df_alt, df_cln = combine_all_alterations(alt=args.alt, cna=args.cna, mut=args.mut, cln=args.cln,
                                             col_gen=col_gen, col_alt=col_alt, col_alt_det=col_alt_det,
                                             col_alt_cat=col_alt_cat, col_sub_id=col_sub_id, col_sam_id=col_sam_id,
                                             col_lvl=col_lvl, subs=args.subject, select_keep_cln=False)

    # parameters for the plot 
    cols_comut = ["sample", "category", "value"]
    sams_ordered = sorted(df_cln[col_sam_id].tolist())
    alts_ordered = sorted(df_alt[col_alt].unique().tolist())

    # remove amplifications/deletions if asked by the user
    if args.mode_cna.lower()=="wo_cna":
        mask_drop = df_alt["Alteration_Category"].isin(["Amplification", "Deletion"])
        df_alt = df_alt.loc[~mask_drop].copy()

    # parameters for the plot 
    cols_comut = ["sample", "category", "value"]
    alts_ordered = sorted(df_alt[col_alt].unique().tolist())[::-1]
    sams_ordered = sorted(df_alt[col_sam_id].unique().tolist())

    # add rows with VAF category for col_alt_cat in order to represent the VAF in the oncoplot
    df_alt, t_vaf_labs = add_vaf_for_mutations(df_alt=df_alt, df_cln=df_cln, col_sub_id=col_sub_id,
                                               col_sam_id=col_sam_id, col_alt=col_alt, col_alt_det=col_alt_det,
                                               col_alt_cat=col_alt_cat, col_t_vaf=col_t_vaf)


    # add Clinical_Benefit_Plot column
    df_cln = add_clinical_benefit_plot(df_cln, mode="simple")
    categories_benefit = ["Yes", "No", "N/A"]
    df_cln["Clinical_Benefit_Plot"] = pd.Categorical(df_cln["Clinical_Benefit_Plot"], categories_benefit)

    # data for the plot
    categorical_data = ["ben", "mtc"]
    categorical_names = ["ORR", "Matched DNA"]
    categorical_pairs = {d: n for d,n in zip(categorical_data, categorical_names)}
    dfs_data = get_tables(df_cln=df_cln, df_alt=df_alt, col_alt=col_alt, col_alt_cat=col_alt_cat, col_sub_id=col_sub_id,
                          col_sam_id=col_sam_id, alts_ordered=alts_ordered, categorical_pairs=categorical_pairs)
    mappings = get_mappings(t_vaf_labs=t_vaf_labs)

    # sizes
    width_middle = 0.5*len(sams_ordered)
    if len(alts_ordered) < 100:
        height_middle = 0.15*len(alts_ordered)
    else:
        height_middle = 0.2*len(alts_ordered)
    height_top = max(6, 0.15*height_middle)
    ignored_plots = ["Alteration Type", "Same patient bis"]
    ignored_values=["Absent"] + t_vaf_labs
    labels_orders = {categorical_pairs["ben"]: categories_benefit}

    # draw figure
    comut_alt = draw_comut_plot(dfs_data=dfs_data, mappings=mappings, sams_ordered=sams_ordered,
                                alts_ordered=alts_ordered, t_vaf_labs=t_vaf_labs, categorical_data=categorical_data,
                                categorical_names=categorical_names, label_bar_top="All\nAlterations",
                                height_middle=height_middle, width_middle=width_middle,
                                height_top=height_top, stacked_bar_side=False, labels_orders=labels_orders,
                                ignored_values=ignored_values, label_bar_top_fontsize=16, show_bar_side=False,
                                show_symbols=True)

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
    parser.add_argument('--mode_cna', type=str, help='Choose "w_cna" or "wo_cna".', default="wo_cna")
    parser.add_argument('--subject', type=str, help='Subject id')
    parser.add_argument('--output', type=str,  help='Paths to output oncoplot-like.')
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
