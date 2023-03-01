# -*- coding: utf-8 -*-
"""
@created: May 03 2022
@modified: Jan 06 2022
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
import sys
sys.path.append("../common/functions")

from utils import combine_all_alterations

# functions ============================================================================================================

def main(args):
    # define util variable names
    col_gen = "Hugo_Symbol"
    col_alt = "Alteration"
    col_alt_det = "Alteration_Detail"
    col_sub_id = "Subject_Id"
    col_sam_id = "Sample_Id"
    col_alt_cat = "Alteration_Category"
    col_lvl = "Sen_Level_Simple"

    df_alt, df_cln = combine_all_alterations(alt=args.alt, cna=args.cna, mut=args.mut, cln=args.cln,
                                             col_gen=col_gen, col_alt=col_alt, col_alt_det=col_alt_det,
                                             col_sub_id=col_sub_id, col_sam_id=col_sam_id, col_alt_cat=col_alt_cat,
                                             col_lvl=col_lvl, samples_select="all", keep_alt_det=True,
                                             high_level_cnas_only=False)

    # select data for gen and add rows for samples not altered
    cols_cln = [col_sub_id, col_sam_id, "Cluster_Id_DNA_P"]
    cols_alt = [col_sub_id, col_sam_id, col_gen, col_alt, "t_vaf", "Annotated"]
    df_alt_gen = df_alt.loc[df_alt[col_gen]==args.gen]
    df_fin_gen = df_cln[cols_cln].merge(df_alt_gen[cols_alt], how="left", on=[col_sub_id, col_sam_id])
    df_fin_gen = df_fin_gen.sort_values(by=[col_sub_id, col_sam_id])
    mask_nna = ~df_fin_gen["t_vaf"].isnull()
    df_fin_gen.loc[mask_nna, "t_vaf"] = df_fin_gen.loc[mask_nna, "t_vaf"].apply(lambda x: "%.3g %%" % (x*100))

    # empty col_alt when equals to gene symbol
    mask_gen = df_fin_gen[col_alt]==df_fin_gen[col_gen]
    df_fin_gen.loc[mask_gen, col_alt] = np.nan

    # select patients
    if args.select:
        df_sub = df_fin_gen[[col_sub_id, "t_vaf"]].dropna(how="any")
        df_fin_gen = df_fin_gen.loc[df_fin_gen[col_sub_id].isin(df_sub[col_sub_id])]

    # save
    df_fin_gen.to_excel(args.output, index=False)
    print("-file saved at %s" % args.output)


# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Oncoplot-like figure detailing treatment resistances.")
    parser.add_argument("--cln", type=str, help="Path to input curated clincal table.")
    parser.add_argument('--alt', type=str, help='Path to table of alterations.')
    parser.add_argument('--mut', type=str, help='Path to table of all mutations.')
    parser.add_argument('--cna', type=str, help='Path to table of all CNAs.')
    parser.add_argument('--gen', type=str, help='Gene symbol')
    parser.add_argument('--select', action="store_true", default=False,
                        help='If True, select only samples from patients bearing an alteration in at least 1 sample')
    parser.add_argument('--output', type=str,  help='Paths to output table.')
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
