import argparse
import numpy as np
import pandas as pd

def main(args):
    if args.cohort in ["prism", "met500"]:
        filepath = "../../data/%s/clinical/curated/bio_%s_in_design_curated.tsv" % (args.cohort, args.cohort)
    elif args.cohort == "tcga":
        filepath = "../../data/%s/clinical/curated_other/bio_%s_all_curated.tsv" % (args.cohort, args.cohort)
    elif args.cohort == "tcga_validation":
        filepath = "../../data/tcga/clinical/curated/bio_tcga_in_design_curated.tsv"

    # read samples data
    df_bio = pd.read_table(filepath, low_memory=False)

    # keep only DNA
    df_bio = df_bio.loc[df_bio["Sample_Type"].isin(["DNA_N", "DNA_T"])]

    if args.cohort == "tcga":
        df_bio["Sample_Id"] = df_bio["Aliquot_Id"]

    # as for MET500 we do not have the target files used, intersection of all target files of PRISM samples
    # will be used as a proxy
    if args.cohort=="met500":
        df_bio["Capture_Kit"] = "all_targets_intersect"
    elif args.cohort=="tcga_validation":
        df_bio["Capture_Kit"] = "all_targets_intersect"

    cols_ids = [x for x in df_bio.columns if x.endswith("_Id")]
    cols_bio = ["Biopsy_Type", "Sample_Type", "Capture_Kit", "Date_Last_State", "Date_Received",
                "Platform", "Study", "Comment_State_Integragen"]
    cols_bio = [x for x in cols_bio if x in df_bio.columns]

    if args.cohort == "prism":
        # extract fastqs data
        df_bio[["FASTQ_1", "FASTQ_2"]] = df_bio["Datafile_Path"].str.split(",").apply(pd.Series)
        df_bio[["FASTQ_1_Size", "FASTQ_2_Size"]] = df_bio["FASTQ_Sizes"].str.split(",").apply(pd.Series)
        df_bio["Dataset"] = df_bio["FASTQ_1"].apply(lambda x: "/".join(x.split("/")[:5]))
        df_bio["FASTQ_1_Name"] = df_bio["FASTQ_1"].apply(lambda x: x.split("/")[-1] if type(x)==str else np.nan)
        df_bio["FASTQ_2_Name"] = df_bio["FASTQ_2"].apply(lambda x: x.split("/")[-1] if type(x)==str else np.nan)

        # irods sample id extracted from fastq name
        mask_a = df_bio["FASTQ_1"].apply(lambda x: x.endswith("_R1.fastq.gz"))
        mask_b = df_bio["FASTQ_2"].apply(lambda x: x.endswith("_R2.fastq.gz"))
        df_bio.loc[mask_a, "Irods_Id"] = df_bio.loc[mask_a]["FASTQ_1_Name"].apply(lambda x: x.split("_R1.fastq.gz")[0])

        mask_a = df_bio["FASTQ_1"].apply(lambda x: x.endswith("_R1_fastq.gz"))
        mask_b = df_bio["FASTQ_2"].apply(lambda x: x.endswith("_R2_fastq.gz"))
        df_bio.loc[mask_a, "Irods_Id"] = df_bio.loc[mask_a]["FASTQ_1_Name"].apply(lambda x: x.split("_R1_fastq.gz")[0])
        df_bio.loc[mask_b, "Irods_Id"] = df_bio.loc[mask_b]["FASTQ_2_Name"].apply(lambda x: x.split("_R2_fastq.gz")[0])
        df_bio.loc[mask_b, "Irods_Id"] = df_bio.loc[mask_b]["FASTQ_2_Name"].apply(lambda x: x.split("_R2_fastq.gz")[0])

        # add status irods for review
        df_bio["Status"] = "OK IRODS"
        cols_fastqs = ["Dataset", "FASTQ_1", "FASTQ_2", "FASTQ_1_Name", "FASTQ_2_Name", "Irods_Id", "Status"]
    elif args.cohort == "met500":
        # extract fastqs data
        df_bio[["FASTQ_1", "FASTQ_2"]] = df_bio["Datafile_Path"].str.split(",").apply(pd.Series)
        df_bio[["FASTQ_1_Size", "FASTQ_2_Size"]] = df_bio["FASTQ_Sizes"].str.split(",").apply(pd.Series)
        df_bio["FASTQ_1_Name"] = df_bio["FASTQ_1"].apply(lambda x: x.split("/")[-1] if type(x)==str else np.nan)
        df_bio["FASTQ_2_Name"] = df_bio["FASTQ_2"].apply(lambda x: x.split("/")[-1] if type(x)==str else np.nan)

        # crap sample id extracted from fastq name
        mask_a = df_bio["FASTQ_1"].apply(lambda x: x.endswith("_1.fastq.gz") if type(x)==str else False)
        mask_b = df_bio["FASTQ_2"].apply(lambda x: x.endswith("_2.fastq.gz") if type(x)==str else False)

        def irods_id(x, i=1):
            if type(x)==str:
                return x.split("_dbGaP-26018_%d.fastq.gz" % i)[0]
            else:
                return x

        df_bio.loc[mask_a, "Irods_Id"] = df_bio.loc[mask_a]["FASTQ_1_Name"].apply(irods_id, i=1)
        df_bio.loc[mask_b, "Irods_Id"] = df_bio.loc[mask_b]["FASTQ_2_Name"].apply(irods_id, i=2)

        # add status irods for review
        cols_fastqs = ["FASTQ_1", "FASTQ_2", "FASTQ_1_Name", "FASTQ_2_Name", "Irods_Id"]
    elif args.cohort == "tcga_validation":
        filepath = "../../data/tcga_validation/wes/bam_grch37/gdc-manifest.txt"
        df_gdc = pd.read_table(filepath)
        df_gdc = df_gdc.rename(columns={"id": "File_Id_Legacy"})

        df_bio["File_Id_Legacy"] = df_bio["File_Id_Legacy"].apply(lambda x: x.split("|") if type(x)==str else [])
        df_bio = df_bio.explode("File_Id_Legacy")
        df_bio = df_bio.loc[df_bio["Sample_Type"].isin(["DNA_T", "DNA_N"])]
        df_bio = df_bio.merge(df_gdc, how="inner", on="File_Id_Legacy")

        cols = ["Subject_Id", "Aliquot_Id", "Sample_Type", "Capture_Kit", "filename", "File_Id_Legacy"]
        df_bio = df_bio[cols]
        df_bio = df_bio.rename(columns={"Aliquot_Id": "Sample_Id", "File_Id_Legacy": "BAM_Id", "filename": "BAM_Name"})
        df_bio["FASTQ_1_Name"] = df_bio["Sample_Id"] + "_1.fastq.gz"
        df_bio["FASTQ_2_Name"] = df_bio["Sample_Id"] + "_2.fastq.gz"

    # add Project_TCGA_More, MSKCC_Oncotree and Civic_Disease
    if args.cohort in ["prism", "met500"]:
        filepath = "../../data/%s/clinical/curated/cln_%s_in_design_curated.tsv" % (args.cohort, args.cohort)
    elif args.cohort=="tcga":
        filepath = "../../data/%s/clinical/curated_other/cln_%s_all_curated.tsv" % (args.cohort, args.cohort)
    elif args.cohort=="tcga_validation":
        filepath = "../../data/tcga/clinical/curated/cln_tcga_in_design_curated.tsv"

    df_cln = pd.read_table(filepath)
    cols_cln = ["Project_TCGA_More", "MSKCC_Oncotree", "Civic_Disease"]
    cols_cln = [x for x in cols_cln if x in df_cln.columns]
    df_bio = df_bio.merge(df_cln[cols_cln + ["Subject_Id"]], how="left", on="Subject_Id")

    # save
    if args.cohort in ["prism", "met500"]:
        cols = cols_ids+cols_bio+cols_fastqs+cols_cln
    elif args.cohort=="tcga":
        cols = cols_ids+cols_bio+cols_cln
    elif args.cohort=="tcga_validation":
        cols = df_bio.columns.tolist()

    if args.n_splits==1:
        filepath = "config/samples.tsv"
        df_bio[cols].to_csv(filepath, index=False, sep="\t")
    else:
        # split tcga samples in batches because otherwise the snakemake DAG is too big
        # and the pipeline becomes very slow
        subjects_all = df_bio["Subject_Id"].unique().tolist()
        n = len(subjects_all)
        s = np.int(np.ceil(n/args.n_splits))
        for i in range(args.n_splits):
            subjects_i = subjects_all[(i*s):((i+1)*s)]
            filepath_i = "config/samples_%d.tsv" % (i+1)
            df_bio_i = df_bio.loc[df_bio["Subject_Id"].isin(subjects_i)]
            df_bio_i[cols].to_csv(filepath_i, index=False, sep="\t")

    # # save
    # filepath = "config/samples_tumor_type_mskcc_oncotree.tsv"
    # df_bio_mskcc_oncotree = df_bio[["Sample_Id", "MSKCC_Oncotree"]].drop_duplicates()
    # df_bio_mskcc_oncotree = df_bio_mskcc_oncotree.rename(columns={"Sample_Id": "SAMPLE_ID",
    #                                                               "MSKCC_Oncotree": "ONCOTREE_CODE"})
    # df_bio_mskcc_oncotree.to_csv(filepath, index=False, sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare samples table.")
    parser.add_argument("--cohort", type=str, help="Choose one of prism, met500, tcga, tcga_validation.",
                        default="tcga_validation")
    parser.add_argument("--n_splits", type=int, default=1, help="Number of splits.")
    args = parser.parse_args()

    main(args)
