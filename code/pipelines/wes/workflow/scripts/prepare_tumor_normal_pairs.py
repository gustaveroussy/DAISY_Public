import argparse
import numpy as np
import pandas as pd

def main(args):
    if args.cohort in ["prism", "met500"]:
        filepath = "../../data/%s/clinical/curated/bio_%s_in_design_curated.tsv" % (args.cohort, args.cohort)
        cols_sel = ["Subject_Id", "Sample_Id", "Sample_Type"]
        df_bio = pd.read_table(filepath)[cols_sel]
        df_bio = df_bio.loc[df_bio["Sample_Type"]!="RNA_T"]

        df = pd.DataFrame(columns=["Subject_Id", "DNA_T", "DNA_N"])
        for subject in df_bio["Subject_Id"].unique():
            sids_n = df_bio.loc[(df_bio["Subject_Id"]==subject) & (df_bio["Sample_Type"]=="DNA_N"), "Sample_Id"]
            sids_t = df_bio.loc[(df_bio["Subject_Id"]==subject) & (df_bio["Sample_Type"]=="DNA_T"), "Sample_Id"]
            if len(sids_n)>1:
                print("-subject %s has the following normals %s" % (subject, ";".join(sids_n)))
            if len(sids_n)==0:
                print("-subject %s has no normal sample" % subject)
            if len(sids_t)==0:
                print("-subject %s has no tumor sample" % subject)
            for sid_t in sids_t:
                for sid_n in sids_n:
                    row = {"Subject_Id": subject, "DNA_T": sid_t, "DNA_N": sid_n}
                    df = df.append(row, ignore_index=True)
    elif args.cohort == "tcga":
        filepath = "../../data/tcga/wes/summary/mc3_unfiltered.tsv"
        df = pd.read_table(filepath)
        df = df.rename(columns={"Tumor_Sample_Id": "DNA_T", "Normal_Sample_Id": "DNA_N"})
        df = df[["Subject_Id", "DNA_T", "DNA_N"]]
    elif args.cohort == "tcga_validation":
        filepath = "../../data/tcga_validation/wes/bam_grch37/gdc-manifest.txt"
        df = pd.read_table(filepath)
        df = df.rename(columns={"id": "File_Id_Legacy", "filename": "File_Name_Legacy"})

        filepath = "../../data/tcga/clinical/curated/bio_tcga_in_design_curated.tsv"
        cols_sel = ["Subject_Id", "Aliquot_Id", "Sample_Type", "File_Id_Legacy"]
        df_bio = pd.read_table(filepath)[cols_sel]
        df_bio["File_Id_Legacy"] = df_bio["File_Id_Legacy"].apply(lambda x: x.split("|") if type(x)==str else [])
        df_bio = df_bio.explode("File_Id_Legacy")
        df_bio = df_bio.loc[df_bio["Sample_Type"].isin(["DNA_T", "DNA_N"])]
        df = df.merge(df_bio, how="left", on="File_Id_Legacy")
        df_n = df.loc[df["Sample_Type"]=="DNA_N"]
        df_n = df_n[["Subject_Id", "Aliquot_Id"]].rename(columns={"Aliquot_Id": "DNA_N"})
        df_t = df.loc[df["Sample_Type"]=="DNA_T"]
        df_t = df_t[["Subject_Id", "Aliquot_Id"]].rename(columns={"Aliquot_Id": "DNA_T"})
        df = df_t.merge(df_n, how="left", on="Subject_Id")

    # add DNA_P
    df["DNA_P"] = df[["DNA_T", "DNA_N"]].fillna("NA").apply("_vs_".join, axis=1)

    # add Project_TCGA_More, MSKCC_Oncotree and Civic_Disease
    if args.cohort in ["prism", "met500"]:
        filepath = "../../data/%s/clinical/curated/cln_%s_in_design_curated.tsv" % (args.cohort, args.cohort)
    elif args.cohort=="tcga":
        filepath = "../../data/%s/clinical/curated_other/cln_%s_all_curated.tsv" % (args.cohort, args.cohort)
    elif args.cohort=="tcga_validation":
        filepath = "../../data/tcga/clinical/curated/cln_tcga_in_design_curated.tsv"

    df_cln = pd.read_table(filepath)
    cols_cln = ["Project_TCGA_More", "MSKCC_Oncotree", "Civic_Disease", "Gender"]
    cols_cln = [x for x in cols_cln if x in df_cln.columns]
    df = df.merge(df_cln[cols_cln + ["Subject_Id"]].drop_duplicates(), how="left", on="Subject_Id")

    # save
    df.to_csv("config/tumor_normal_pairs.tsv", index=False, sep="\t")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare samples table.")
    parser.add_argument("--cohort", type=str, help="Choose one of prism, met500, tcga.", default="met500")
    args = parser.parse_args()

    main(args)
