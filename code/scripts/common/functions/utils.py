import gzip
import numpy as np
import os
import pandas as pd
import re

def read_header(path, prefix):
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as file:
            header = [x for x in file.readlines() if x.startswith(prefix)]
    else:
        with open(path, "r") as file:
            header = [x for x in file.readlines() if x.startswith(prefix)]
    return header


def join_unique(x):
    return "|".join(list(set(x)))


def concatenate_cna(df_alt, df_cna, col_sam_id, col_gen, col_alt, col_alt_det, col_alt_cat, col_alt_cat_sim):
    if df_cna.shape[0]==0:
        return df_alt
    else:
        cols_ids = ["Aliquot_Id", "Sample_Id", "Cluster_Id", "Subject_Id"]
        cols_ids = [x for x in cols_ids if x in df_cna]

        cols_det = ["Copy_Number_More", "TCN_EM:LCN_EM"]
        cols = cols_ids + [col_gen, col_alt_cat, col_alt_cat_sim, col_alt, col_alt_det, "Annotated"]
        df_cna["Annotated"] = "No"
        df_cna[col_gen] = df_cna["Hugo_Symbol"]
        mask_del = df_cna["Copy_Number"]==-2
        df_cna.loc[mask_del, col_alt_cat] = "Deletion"
        df_cna.loc[mask_del, col_alt_cat_sim] = "Del"
        df_cna.loc[mask_del, col_alt] = "Del"
        df_cna.loc[mask_del, col_alt_det] = df_cna.loc[mask_del, cols_det].apply(" - ".join, axis=1)

        mask_amp = df_cna["Copy_Number"]==2
        df_cna.loc[mask_amp, col_alt_cat] = "Amplification"
        df_cna.loc[mask_amp, col_alt_cat_sim] = "Amp"
        df_cna.loc[mask_amp, col_alt] = "Amp"
        df_cna.loc[mask_amp, col_alt_det] = df_cna.loc[mask_amp, cols_det].apply(" - ".join, axis=1)

        mask_oth = (~mask_del & ~mask_amp) & (~df_cna["Copy_Number_More"].isnull())
        df_cna.loc[mask_oth, col_alt_det] = df_cna.loc[mask_oth, cols_det].apply(" - ".join, axis=1)
        df_cna.loc[mask_oth, col_alt] = df_cna.loc[mask_oth, "Copy_Number_More"]
        df_cna.loc[mask_oth, col_alt_cat] = "CNA_Other"
        df_cna.loc[mask_oth, col_alt_cat_sim] = "CNA_Other"

        df_alt = pd.concat((df_alt, df_cna[cols]), axis=0)
        df_alt = df_alt.drop_duplicates(subset=[col_sam_id, col_gen, col_alt_det], keep="first")

        return df_alt


def concatenate_mut(df_alt, df_mut, col_sam_id, col_gen, col_alt, col_alt_det, col_alt_cat, col_alt_cat_sim):
    if df_mut.shape[0]==0:
        return df_alt
    else:
        cols_ids = ["Aliquot_Id", "Sample_Id", "Cluster_Id", "Subject_Id"]
        cols_ids = [x for x in cols_ids if x in df_mut]
        cols = cols_ids + [col_gen, "Variant_Classification", col_alt_cat, col_alt_cat_sim, col_alt, col_alt_det,
                           "Annotated", "t_vaf"]

        df_mut["Annotated"] = "No"
        df_mut[col_gen] = df_mut["Hugo_Symbol"]
        get_mut_det = lambda x: re.sub("^p.", "", x["HGVSp_Short"]) if type(x["HGVSp_Short"])==str else x["HGVSc"]
        hgvsp_short_split = lambda x: x.split("p.")[1].replace("%3D", "=") if type(x)==str else x
        exon_split = lambda x: x.split("/")[0] if type(x)==str else ""

        mask_ins = df_mut["Variant_Classification"].isin(["Frame_Shift_Ins", "In_Frame_Ins"])
        df_mut.loc[mask_ins, col_alt_cat] = "Indel"
        df_mut.loc[mask_ins, col_alt_cat_sim] = "Mut"
        df_mut.loc[mask_ins, col_alt] = "Exon " + df_mut.loc[mask_ins,"EXON"].apply(exon_split) + " Ins"
        df_mut.loc[mask_ins, col_alt_det] = df_mut.loc[mask_ins, ["HGVSp_Short", "HGVSc"]].apply(get_mut_det, axis=1)

        mask_del = df_mut["Variant_Classification"].isin(["Frame_Shift_Del", "In_Frame_Del"])
        df_mut.loc[mask_del, col_alt_cat] = "Indel"
        df_mut.loc[mask_del, col_alt_cat_sim] = "Mut"
        df_mut.loc[mask_del, col_alt] = "Exon " + df_mut.loc[mask_del,"EXON"].apply(exon_split) + " Del"
        df_mut.loc[mask_del, col_alt_det] = df_mut.loc[mask_del, ["HGVSp_Short", "HGVSc"]].apply(get_mut_det, axis=1)

        mask_mut = (~mask_ins) & (~mask_del)
        mask_nul_hgvsp = df_mut["HGVSp_Short"].isnull()
        df_mut.loc[mask_mut, col_alt_cat] = "Mutation"
        df_mut.loc[mask_mut, col_alt_cat_sim] = "Mut"
        df_mut.loc[mask_mut, col_alt] = df_mut.loc[mask_mut, "HGVSp_Short"].apply(hgvsp_short_split)
        df_mut.loc[mask_mut & mask_nul_hgvsp, col_alt] = \
                df_mut.loc[mask_mut & mask_nul_hgvsp, "Variant_Classification"]
        df_mut.loc[mask_mut, col_alt_det] = df_mut.loc[mask_mut, ["HGVSp_Short", "HGVSc"]].apply(get_mut_det, axis=1)

        # for splice, set col_alt to Splice_Site
        mask_splice = df_mut["Variant_Classification"]=="Splice_Site"
        df_mut.loc[mask_splice, col_alt] = "Splice_Site"
        mask_splice = df_mut[col_alt_det].apply(lambda x: "splice" in x if type(x)==str else False)
        df_mut.loc[mask_splice, col_alt_det] = "Splice_Site"

        # fill null col_alt_det
        mask_null = df_mut[col_alt_det].isnull()
        df_mut.loc[mask_null, col_alt_det] = df_mut.loc[mask_null, "Variant_Classification"]

        df_alt = pd.concat((df_alt, df_mut[cols]), axis=0)
        df_alt = df_alt.drop_duplicates(subset=[col_sam_id, col_gen, col_alt_det], keep="first")

        return df_alt


def process_alt(df_alt, df_cna, df_mut, sam_list, col_sam_id, col_gen, col_alt, col_alt_det, col_alt_cat,
                col_alt_cat_sim, keep_alt_det=False):

    if df_cna is not None:
        df_cna = df_cna.loc[df_cna[col_sam_id].isin(sam_list)].copy()
    if df_mut is not None:
        df_mut = df_mut.loc[df_mut[col_sam_id].isin(sam_list)].copy()

    # concatenate cna
    if df_cna is not None:
        df_alt = concatenate_cna(df_alt, df_cna, col_sam_id, col_gen, col_alt, col_alt_det, col_alt_cat, col_alt_cat_sim)

    # concatenate mut
    if df_mut is not None:
        df_alt = concatenate_mut(df_alt, df_mut, col_sam_id, col_gen, col_alt, col_alt_det, col_alt_cat, col_alt_cat_sim)

    # there should not be any null col_alt_det
    assert df_alt[col_alt_det].isnull().sum()==0

    # concatenate gene name and alteration
    # for add gene symbol as prefix except for
    #  - alterations that are already equal to gene symbol
    df_gen_alt = df_alt[[col_gen, col_alt]].drop_duplicates().fillna("NA")
    df_gen_alt_gby = df_gen_alt.groupby(col_gen)[col_alt].agg("|".join)
    mask_a = (df_alt[col_alt] != df_alt[col_gen])
    df_alt.loc[mask_a, col_alt] = df_alt.loc[mask_a, [col_gen, col_alt]].apply(lambda x: " ".join(x.dropna()), axis=1)

    # if don't keep alteration detail, only the gene name is kept unless only 1 alt is associated with the gene,
    # in which case replace by alt
    if not keep_alt_det:
        df_gen_alt = df_alt[[col_gen, col_alt]].drop_duplicates()
        df_gen_alt_gby = df_gen_alt.groupby(col_gen).size()
        gen_1_alt_only = df_gen_alt_gby[df_gen_alt_gby==1].index.tolist()
        mask_1_alt_only = df_alt[col_gen].isin(gen_1_alt_only)
        df_alt.loc[~mask_1_alt_only,col_alt] = df_alt.loc[~mask_1_alt_only, col_gen]

    return df_alt


def combine_all_alterations(alt, cln, cna, mut, col_gen, col_alt, col_alt_det, col_alt_cat, col_sub_id, col_sam_id,
                            subs=None, samples_select="all", keep_alt_det=False, cna_selection="focal_all"):

    # get data
    df_alt = pd.read_table(alt)
    df_cln = pd.read_table(cln)

    # create col_alt_cat_sim for naming the alerations and group some categories for coloring
    col_alt_cat_sim = "%s_Simple" % col_alt_cat
    df_alt[col_alt_cat_sim] = df_alt[col_alt_cat].replace({"Del": "Mut", "Ins": "Mut", "Amplification": "Amp",
                                                           "Deletion": "Del"})
    df_alt[col_alt_cat] = df_alt[col_alt_cat].replace({"Mut": "Mutation", "Ins": "Indel", "Del": "Indel"})

    if cna is not None:
        df_cna = pd.read_table(cna)
        if cna_selection=="focal_high_level":
            mask_keep = ~df_cna["Copy_Number"].isnull()
            df_cna = df_cna.loc[mask_keep].copy()
            print("-retained %d/%d CNAs that are focal and HL/ML or HD" % (sum(mask_keep), mask_keep.shape[0]))
        elif cna_selection=="focal_non_neutral":
            mask_keep = ~df_cna["Copy_Number_More"].isnull()
            df_cna = df_cna.loc[mask_keep].copy()
            print("-retained %d/%d CNAs that are focal and non-neutral" % (sum(mask_keep), mask_keep.shape[0]))
        elif cna_selection=="focal_all":
            print("-retained all %d focal CNAs (some may be neutral)" % df_cna.shape[0])
        else:
            choices = ["focal_high_level", "focal_non_neutral", "focal_all"]
            raise ValueError("Unknown value '%s' for cna_selection. Choose one of: %s" % (cna_selection, choices))
    else:
        df_cna = None

    if mut is not None:
        df_mut = pd.read_table(mut, sep="\t", skiprows=len(read_header(path=mut, prefix="##")), low_memory=False)
        df_mut["HGVSp_Short"] = df_mut["HGVSp_Short"].apply(lambda x: x.replace("%3D", "=") if type(x)==str else x)
    else:
        df_mut = None

    # add Sample_Id to all tables and select only samples that pass QC unless specified otherwise by the user
    df_cln = df_cln.loc[df_cln["QC_Final_Decision"]==1]
    col_tsb = "Tumor_Sample_Barcode"
    col_nsb = "Matched_Norm_Sample_Barcode"
    df_cln["Sample_Id"] = df_cln["Sample_Id_DNA_T"]
    df_cln["Cluster_Id"] = df_cln["Cluster_Id_DNA_T"]
    cols_cln = ["Cluster_Id_DNA_P","Sample_Id","Subject_Id", "Cluster_Id"]

    if df_cna is not None:
        df_cna["Cluster_Id_DNA_P"] = df_cna[[col_tsb, col_nsb]].fillna("NA").apply("_vs_".join, axis=1)
        df_cna = df_cna.loc[df_cna["Cluster_Id_DNA_P"].isin(df_cln["Cluster_Id_DNA_P"])].copy()
        df_cna = df_cna.merge(df_cln[cols_cln], how="left", on="Cluster_Id_DNA_P")

    if df_mut is not None:
        df_mut["Cluster_Id_DNA_P"] = df_mut[[col_tsb, col_nsb]].fillna("NA").apply("_vs_".join, axis=1)
        df_mut = df_mut.loc[df_mut["Cluster_Id_DNA_P"].isin(df_cln["Cluster_Id_DNA_P"])].copy()
        df_mut = df_mut.merge(df_cln[cols_cln], how="left", on="Cluster_Id_DNA_P")

        # for mut, add t_vaf
        df_mut["t_vaf"] = df_mut["t_alt_count"]/df_mut["t_depth"]

    # select samples
    mask_t1 = df_cln["Sample_Id_DNA_T"].str.endswith("T1").fillna(False)
    mask_t2 = df_cln["Sample_Id_DNA_T"].str.endswith("T2").fillna(False)
    df_cln.loc[mask_t1, "Cohort"] = "DAISY-T1"
    df_cln.loc[mask_t2, "Cohort"] = "DAISY-T2"

    dfs_cln = []
    subs_t1 = df_cln.loc[df_cln["Cohort"]=="DAISY-T1", col_sub_id].tolist()
    subs_t2 = df_cln.loc[df_cln["Cohort"]=="DAISY-T2", col_sub_id].tolist()
    df_cnt_bpn = df_cln.groupby(col_sub_id)["Biopsy_Number"].nunique()
    subs_t1t2 = df_cnt_bpn[df_cnt_bpn > 1].index.tolist()
    sams_t1 = df_cln.loc[df_cln["Cohort"]=="DAISY-T1", col_sam_id].tolist()
    sams_t2 = df_cln.loc[df_cln["Cohort"]=="DAISY-T2", col_sam_id].tolist()
    if samples_select=="all":
        dfs_cln.append(df_cln)
    if "t1" in samples_select.split("_"):
        dfs_cln.append(df_cln.loc[df_cln[col_sam_id].isin(sams_t1)].copy())
    if "t2" in samples_select.split("_"):
        dfs_cln.append(df_cln.loc[df_cln[col_sam_id].isin(sams_t2)].copy())
    if "t1only" in samples_select.split("_"):
        subs_t1only = list(set(subs_t1).difference(set(subs_t2)))
        dfs_cln.append(df_cln.loc[df_cln[col_sub_id].isin(subs_t1only)].copy())
    if "t2only" in samples_select.split("_"):
        subs_t2only = list(set(subs_t2).difference(set(subs_t1)))
        dfs_cln.append(df_cln.loc[df_cln[col_sub_id].isin(subs_t2only)].copy())
    elif "t1t2" in samples_select.split("_"):
        dfs_cln.append(df_cln.loc[df_cln[col_sub_id].isin(subs_t1t2)].copy())
    df_cln = pd.concat(dfs_cln)

    if subs is not None:
        if type(subs)==str:
            df_cln = df_cln.loc[df_cln[col_sub_id]==subs]
        elif type(subs)==list:
            df_cln = df_cln.loc[df_cln[col_sub_id].isin(subs)]
        else:
            raise ValueError("-unsupported type for subs")

    # samples selection
    sam_list = df_cln[col_sam_id].tolist()

    # combined annotated and not annotated events
    df_alt["Annotated"] = "Yes"
    df_alt = df_alt.loc[df_alt[col_sam_id].isin(sam_list)].copy()
    # if keep_alt_det:
    #     # where col_alt_det is not NA, replace the value of col_alt by col_alt_det
    #     # # exception: for Amplification and Deletion, ignore col_alt_det
    #     # mask_amp_del = df_alt["Alteration_Category"].isin(["Amplification", "Deletion"])
    #     # df_alt.loc[mask_amp_del, col_alt_det] = np.nan
    #     mask_alt_det_nna = ~df_alt[col_alt_det].isnull()
    #     df_alt.loc[mask_alt_det_nna, col_alt] = df_alt.loc[mask_alt_det_nna, col_alt_det]

    # process alterations
    df_alt = process_alt(df_alt, df_cna, df_mut, sam_list, col_sam_id, col_gen, col_alt, col_alt_det, col_alt_cat,
                         col_alt_cat_sim, keep_alt_det)

    return df_alt, df_cln


def combine_tcga_alterations(cln, cna, mut, col_gen, col_alt, col_alt_det, col_alt_cat, col_sub_id,
                             col_sam_id, keep_alt_det=False, cna_selection="focal_all"):

    col_alt_cat_sim = "%s_Simple " % col_alt_cat
    df_alt = pd.DataFrame()

    # load cln
    df_cln = pd.read_table(cln)
    df_cln[col_sam_id] = df_cln["Sample_Id_DNA_T"]

    # load cna
    if cna is not None:
        df_cna = pd.read_table(cna)
        df_cna[col_sub_id] = df_cna["Tumor_Sample_Barcode"].str[:12]
        df_cna[col_sam_id] = df_cna["Tumor_Sample_Barcode"]
        df_ids_cna = pd.read_table(os.path.join(os.path.dirname(cna), "sample_list.tsv"))

        if cna_selection=="focal_high_level":
            mask_keep = ~df_cna["Copy_Number"].isnull()
            df_cna = df_cna.loc[mask_keep].copy()
            print("-retained %d/%d CNAs that are focal and HL/ML or HD" % (sum(mask_keep), mask_keep.shape[0]))
        elif cna_selection=="focal_non_neutral":
            mask_keep = ~df_cna["Copy_Number_More"].isnull()
            df_cna = df_cna.loc[mask_keep].copy()
            print("-retained %d/%d CNAs that are focal and non-neutral" % (sum(mask_keep), mask_keep.shape[0]))
        elif cna_selection=="focal_all":
            print("-retained all %d focal CNAs (some may be neutral)" % df_cna.shape[0])
        else:
            choices = ["focal_high_level", "focal_non_neutral", "focal_all"]
            raise ValueError("Unknown value '%s' for cna_selection. Choose one of: %s" % (cna_selection, choices))

    else:
        df_cna = None

    # load mut
    if mut is not None:
        df_mut = pd.read_table(mut, sep="\t", skiprows=len(read_header(path=mut, prefix="##")))
        df_mut["HGVSp_Short"] = df_mut["HGVSp_Short"].apply(lambda x: x.replace("%3D", "=") if type(x)==str else x)
        df_mut[col_sub_id] = df_mut["Tumor_Sample_Barcode"].str[:12]
        df_mut[col_sam_id] = df_mut["Tumor_Sample_Barcode"]
        df_mut["t_vaf"] = df_mut["t_alt_count"]/df_mut["t_depth"]
        df_ids_mut = pd.read_table(os.path.join(os.path.dirname(mut), "sample_list.tsv"))
    else:
        df_mut = None

    # select ids common to cln, mut, cna
    sam_list = df_cln[col_sam_id].unique().tolist()

    if df_cna is not None:
        df_cna = df_cna.loc[df_cna[col_sam_id].isin(sam_list)].copy()
        sam_list = list(set(sam_list).intersection(set(df_ids_cna["Tumor_%s" % col_sam_id])))
    if df_mut is not None:
        df_mut = df_mut.loc[df_mut[col_sam_id].isin(sam_list)].copy()
        sam_list = list(set(sam_list).intersection(set(df_ids_mut["Tumor_%s" % col_sam_id])))

    # filter on common samples
    mask_keep = df_cln[col_sam_id].isin(sam_list)
    print("-retained %d/%d patients with all data available" % (sum(mask_keep), mask_keep.shape[0]))
    df_cln = df_cln.loc[mask_keep].copy()

    # process alterations
    df_alt = process_alt(df_alt, df_cna, df_mut, sam_list, col_sam_id, col_gen, col_alt, col_alt_det, col_alt_cat,
                         col_alt_cat_sim, keep_alt_det)

    return df_alt, df_cln
