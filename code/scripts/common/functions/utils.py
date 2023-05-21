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


def concatenate_cna(df_alt, df_cna, col_gen, col_alt, col_alt_det, col_alt_cat, col_alt_cat_sim, keep_alt_det):
    if df_cna.shape[0]==0:
        return df_alt
    else:
        cols_ids = ["Aliquot_Id", "Sample_Id", "Cluster_Id", "Subject_Id",]
        cols_ids = [x for x in cols_ids if x in df_cna]

        cols_det = ["Copy_Number_More", "TCN_EM:LCN_EM"]
        cols = cols_ids + [col_gen, col_alt_cat, col_alt_cat_sim, col_alt, col_alt_det, "Annotated"]
        df_cna["Annotated"] = "No"
        df_cna[col_gen] = df_cna["Hugo_Symbol"]
        mask_del = df_cna["Copy_Number"]==-2
        df_cna.loc[mask_del, col_alt_cat] = "Deletion"
        df_cna.loc[mask_del, col_alt_cat_sim] = "Del"
        df_cna.loc[mask_del, col_alt_det] = df_cna.loc[mask_del, cols_det].apply(" - ".join, axis=1)
        if keep_alt_det:
            df_cna.loc[mask_del, col_alt] = df_cna.loc[mask_del, col_alt_det]
        else:
            df_cna.loc[mask_del, col_alt] = "Del"

        mask_amp = df_cna["Copy_Number"]==2
        df_cna.loc[mask_amp, col_alt_cat] = "Amplification"
        df_cna.loc[mask_amp, col_alt_cat_sim] = "Amp"
        df_cna.loc[mask_amp, col_alt_det] = df_cna.loc[mask_amp, cols_det].apply(" - ".join, axis=1)
        if keep_alt_det:
            df_cna.loc[mask_amp, col_alt] = df_cna.loc[mask_amp, col_alt_det]
        else:
            df_cna.loc[mask_amp, col_alt] = "Amp"

        return pd.concat((df_alt, df_cna[cols]), axis=0)


def concatenate_mut(df_alt, df_mut, col_gen, col_alt, col_alt_det, col_alt_cat, col_alt_cat_sim):
    if df_mut.shape[0]==0:
        return df_alt
    else:
        cols_ids = ["Aliquot_Id", "Sample_Id", "Cluster_Id", "Subject_Id",]
        cols_ids = [x for x in cols_ids if x in df_mut]
        cols = cols_ids + [col_gen, "Variant_Classification", col_alt_cat, col_alt_cat_sim, col_alt, col_alt_det,
                           "Annotated", "t_vaf"]

        df_mut["Annotated"] = "No"
        df_mut[col_gen] = df_mut["Hugo_Symbol"]
        hgvsp_short_split = lambda x: x.split("p.")[1].replace("%3D", "=") if type(x)==str else x
        exon_split = lambda x: x.split("/")[0] if type(x)==str else ""
        mask_ins = df_mut["Variant_Classification"].isin(["Frame_Shift_Ins", "In_Frame_Ins"])
        mask_del = df_mut["Variant_Classification"].isin(["Frame_Shift_Del", "In_Frame_Del"])
        df_mut.loc[mask_ins, col_alt_cat] = "Indel"
        df_mut.loc[mask_ins, col_alt_cat_sim] = "Mut"
        df_mut.loc[mask_ins, col_alt] = "Exon " + df_mut.loc[mask_ins,"EXON"].apply(exon_split) + " Ins"
        df_mut.loc[mask_ins, col_alt_det] = df_mut.loc[mask_ins, "HGVSp_Short"].apply(hgvsp_short_split)
        df_mut.loc[mask_del, col_alt_cat] = "Indel"
        df_mut.loc[mask_del, col_alt_cat_sim] = "Mut"
        df_mut.loc[mask_del, col_alt_det] = df_mut.loc[mask_del, "HGVSp_Short"].apply(hgvsp_short_split)
        df_mut.loc[mask_del, col_alt] = "Exon " + df_mut.loc[mask_del,"EXON"].apply(exon_split) + " Del"

        mask_mut = (~mask_ins) & (~mask_del)
        mask_nul_hgvsp = df_mut["HGVSp_Short"].isnull()
        df_mut.loc[mask_mut, col_alt_cat] = "Mutation"
        df_mut.loc[mask_mut, col_alt_cat_sim] = "Mut"
        df_mut.loc[mask_mut, col_alt] = df_mut.loc[mask_mut, "HGVSp_Short"].apply(hgvsp_short_split)
        df_mut.loc[mask_mut & mask_nul_hgvsp, col_alt] = \
                df_mut.loc[mask_mut & mask_nul_hgvsp, "Variant_Classification"]

        # for splice, set Alteration to Splice_Site
        mask_splice = df_mut["Variant_Classification"]=="Splice_Site"
        df_mut.loc[mask_splice, col_alt] = "Splice_Site"

        return pd.concat((df_alt, df_mut[cols]), axis=0)


def process_alt(df_alt, df_cna, df_mut, samples_all, col_sam_id, col_gen, col_alt, col_alt_det, col_alt_cat,
                col_alt_cat_sim, col_lvl=None, keep_alt_det=False):

    if df_cna is not None:
        df_cna = df_cna.loc[df_cna[col_sam_id].isin(samples_all)].copy()
    if df_mut is not None:
        df_mut = df_mut.loc[df_mut[col_sam_id].isin(samples_all)].copy()

    # concatenate cna
    if df_cna is not None:
        df_alt = concatenate_cna(df_alt, df_cna, col_gen, col_alt, col_alt_det, col_alt_cat, col_alt_cat_sim,
                                 keep_alt_det)

    # concatenate mut
    if df_mut is not None:
        df_alt = concatenate_mut(df_alt, df_mut, col_gen, col_alt, col_alt_det, col_alt_cat, col_alt_cat_sim)

    # where col_alt_det is NA, replace the value of col_alt_det by col_alt
    mask_col_alt_det_na = df_alt[col_alt_det].isnull()
    df_alt.loc[mask_col_alt_det_na, col_alt_det] = df_alt.loc[mask_col_alt_det_na, col_alt]

    # concatenate gene name and alteration
    # for add gene symbol as prefix except for
    #  - alterations that are already equal to gene symbol
    #  - genes with only "Oncogenic, no Level"
    df_gen_alt = df_alt[[col_gen, col_alt]].drop_duplicates().fillna("NA")
    df_gen_alt_gby = df_gen_alt.groupby(col_gen)[col_alt].agg("|".join)
    genes_no_lvl_only = df_gen_alt_gby[df_gen_alt_gby=="Oncogenic, no Level"].index.tolist()
    mask_a = df_alt[col_gen].isin(genes_no_lvl_only)
    mask_b = (~mask_a) & (df_alt[col_alt] != df_alt[col_gen])
    df_alt.loc[mask_a,col_alt] = df_alt.loc[mask_a, col_gen]
    df_alt.loc[mask_b,col_alt] = df_alt.loc[mask_b, [col_gen, col_alt]].apply(lambda x: " ".join(x.dropna()), axis=1)

    # for alterations where alterations is gene + Oncogenic, no Level is only
    #   - if only 1 alt, replace by alt
    #   - otherwise, replace by concatenated categories
    if not keep_alt_det:
        mask_c = df_alt[col_alt] == (df_alt[col_gen] + " Oncogenic, no Level")
        df_gen_alt_det = df_alt.loc[(~mask_a) & mask_c][[col_gen, col_alt, col_alt_det]].drop_duplicates()
        df_gen_alt_det_gby = df_gen_alt_det.groupby(col_alt)[col_alt_det].size()
        alt_1_alt_only = df_gen_alt_det_gby[df_gen_alt_det_gby==1].index.tolist()
        mask_d = df_alt[col_alt].isin(alt_1_alt_only)
        df_alt.loc[mask_d,col_alt] = df_alt.loc[mask_d, [col_gen, col_alt_det]].apply(" ".join, axis=1)

        # for alterations that are gene + Oncogenic, no Level or gene,
        # and for which other alterations from the same gene coexist, 
        # replace col_alt by col_gene + col_alt_cat (grouped by col_alt)
        mask_e = (~mask_a) & mask_c & (~mask_d)
        if sum(mask_e) > 0:
            col_alt_new = "%s_New" % col_alt
            df_gen_alt_cat = df_alt.loc[mask_e][[col_gen, col_alt, col_alt_cat_sim]].drop_duplicates()
            df_gen_alt_cat = df_gen_alt_cat.groupby([col_gen, col_alt]).agg({col_alt_cat_sim: join_unique}).reset_index()
            df_gen_alt_cat[col_alt_new] = df_gen_alt_cat[[col_gen, col_alt_cat_sim]].apply(" ".join, axis=1)

            df_alt_mask = df_alt.loc[mask_e]
            df_alt_othr = df_alt.loc[~mask_e]
            df_alt_mask = df_alt_mask.merge(df_gen_alt_cat[[col_gen, col_alt, col_alt_new]], how="left", on=[col_gen,
                                                                                                             col_alt])
            df_alt_mask[col_alt] = df_alt_mask[col_alt_new]
            del df_alt_mask[col_alt_new]
            df_alt = pd.concat((df_alt_mask, df_alt_othr))

        # for genes where col alt is gen and not all alterations have the same level,
        #  - if col alt is only 1 alteration, replace by alteration
        #  - otherwise replace, by categories concatenate
        df_gen_lvl = df_alt[[col_gen, col_lvl]].drop_duplicates()
        df_gen_lvl[col_lvl] = df_gen_lvl[col_lvl].fillna("N/A")
        df_gen_lvl_gby = df_gen_lvl.groupby(col_gen)[col_lvl].size()
        genes_one_lvl_only = df_gen_lvl_gby[df_gen_lvl_gby==1].index.tolist()
        mask_one_lvl = df_alt[col_gen].isin(genes_one_lvl_only)
        mask_gen = (df_alt[col_gen] == df_alt[col_alt]) & ~(mask_one_lvl) & ~(mask_a)
        df_gen_alt_det = df_alt.loc[mask_gen][[col_gen, col_alt, col_alt_det]].drop_duplicates()
        df_gen_alt_det_gby = df_gen_alt_det.groupby(col_alt)[col_alt_det].size()
        alt_1_alt_only = df_gen_alt_det_gby[df_gen_alt_det_gby==1].index.tolist()
        mask_f_1 = df_alt[col_alt].isin(alt_1_alt_only) & ~(mask_one_lvl) & ~(mask_a) & mask_gen
        mask_f_2 = ~(df_alt[col_alt].isin(alt_1_alt_only)) & ~(mask_one_lvl) & ~(mask_a) & mask_gen

        df_alt.loc[mask_f_1,col_alt] = df_alt.loc[mask_f_1, [col_gen, col_alt_det]].apply(" ".join, axis=1)
        df_alt.loc[mask_f_2,col_alt] = df_alt.loc[mask_f_2, [col_gen, col_alt_cat_sim]].apply(" ".join, axis=1)

        # for genes where col alt is not gen,
        #  - if col alt is only 1 alteration, replace by alteration
        df_gen_alt_det = df_alt.loc[~mask_gen][[col_gen, col_alt, col_alt_det]].drop_duplicates()
        df_gen_alt_det_gby = df_gen_alt_det.groupby(col_alt)[col_alt_det].size()
        alt_1_alt_only = df_gen_alt_det_gby[df_gen_alt_det_gby==1].index.tolist()
        mask_g = df_alt[col_alt].isin(alt_1_alt_only)
        df_alt.loc[mask_g,col_alt] = df_alt.loc[mask_g, [col_gen, col_alt_det]].apply(" ".join, axis=1)

    # drop duplicated alterations, priotirizing Annotated=="Yes"
    df_alt = df_alt.drop_duplicates(subset=[col_sam_id, col_alt], keep="first").copy()

    return df_alt


def combine_all_alterations(alt, cln, cna, mut, col_gen, col_alt, col_alt_det, col_alt_cat, col_sub_id, col_sam_id,
                            col_lvl=None, subs=None, samples_select="all", keep_alt_det=False,
                            high_level_cnas_only=True, select_keep_cln=True):

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
        if high_level_cnas_only:
            mask_keep = ~df_cna["Copy_Number"].isnull()
            df_cna = df_cna.loc[mask_keep].copy()
            print("-retained %d/%d CNAs that are -2 or +2" % (sum(mask_keep), mask_keep.shape[0]))
    else:
        df_cna = None

    if mut is not None:
        df_mut = pd.read_table(mut, sep="\t", skiprows=len(read_header(path=mut, prefix="##")), low_memory=False)
        df_mut["HGVSp_Short"] = df_mut["HGVSp_Short"].apply(lambda x: x.replace("%3D", "=") if type(x)==str else x)
    else:
        df_mut = None

    # add Sample_Id to all tables
    col_tsb = "Tumor_Sample_Barcode"
    col_nsb = "Matched_Norm_Sample_Barcode"
    df_cln["Sample_Id"] = df_cln["Sample_Id_DNA_T"]
    df_cln["Cluster_Id"] = df_cln["Cluster_Id_DNA_T"]
    cols_cln = ["Cluster_Id_DNA_P","Sample_Id","Subject_Id", "Cluster_Id"]

    if df_cna is not None:
        df_cna["Cluster_Id_DNA_P"] = df_cna[[col_tsb, col_nsb]].fillna("NA").apply("_vs_".join, axis=1)
        df_cna = df_cna.merge(df_cln[cols_cln], how="left", on="Cluster_Id_DNA_P")

    if df_mut is not None:
        df_mut["Cluster_Id_DNA_P"] = df_mut[[col_tsb, col_nsb]].fillna("NA").apply("_vs_".join, axis=1)
        df_mut = df_mut.merge(df_cln[cols_cln], how="left", on="Cluster_Id_DNA_P")

        # for mut, add t_vaf
        df_mut["t_vaf"] = df_mut["t_alt_count"]/df_mut["t_depth"]


    # select only samples that pass QC unless specified otherwise by the user
    if select_keep_cln:
        df_cln = df_cln.loc[df_cln["QC_Final_Decision"]==1]

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
    samples_all = df_cln[col_sam_id].tolist()

    # combined annotated and not annotated events
    df_alt["Annotated"] = "Yes"
    df_alt = df_alt.loc[df_alt[col_sam_id].isin(samples_all)].copy()
    if keep_alt_det:
        # where col_alt_det is not NA, replace the value of col_alt by col_alt_det
        # # exception: for Amplification and Deletion, ignore col_alt_det
        # mask_amp_del = df_alt["Alteration_Category"].isin(["Amplification", "Deletion"])
        # df_alt.loc[mask_amp_del, col_alt_det] = np.nan
        mask_alt_det_nna = ~df_alt[col_alt_det].isnull()
        df_alt.loc[mask_alt_det_nna, col_alt] = df_alt.loc[mask_alt_det_nna, col_alt_det]

    # process alterations
    df_alt = process_alt(df_alt, df_cna, df_mut, samples_all, col_sam_id, col_gen, col_alt, col_alt_det, col_alt_cat,
                         col_alt_cat_sim, col_lvl, keep_alt_det)

    return df_alt, df_cln


def combine_tcga_alterations(cln, cna, mut, col_gen, col_alt, col_alt_det, col_alt_cat, col_lvl, col_sub_id,
                             col_sam_id, keep_alt_det=False):

    col_alt_cat_sim = "%s_Simple " % col_alt_cat
    df_alt = pd.DataFrame()

    # load cln
    df_cln = pd.read_table(cln)
    df_cln[col_sam_id] = df_cln["Sample_Id_DNA_T"]

    # load cna
    if cna is not None:
        df_cna = pd.read_table(cna)
        mask_keep = ~df_cna["Copy_Number"].isnull()
        df_cna = df_cna.loc[mask_keep].copy()
        print("-retained %d/%d CNAs that are -2 or +2" % (sum(mask_keep), mask_keep.shape[0]))
        df_cna[col_sub_id] = df_cna["Tumor_Sample_Barcode"].str[:12]
        df_cna[col_sam_id] = df_cna["Tumor_Sample_Barcode"]
        df_ids_cna = pd.read_table(os.path.join(os.path.dirname(cna), "sample_list.tsv"))
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
    samples_all = df_cln[col_sam_id].unique().tolist()

    if df_cna is not None:
        df_cna = df_cna.loc[df_cna[col_sam_id].isin(samples_all)].copy()
        samples_all = list(set(samples_all).intersection(set(df_ids_cna["Tumor_%s" % col_sam_id])))
    if df_mut is not None:
        df_mut = df_mut.loc[df_mut[col_sam_id].isin(samples_all)].copy()
        samples_all = list(set(samples_all).intersection(set(df_ids_mut["Tumor_%s" % col_sam_id])))

    # filter on common samples
    mask_keep = df_cln[col_sam_id].isin(samples_all)
    print("-retained %d/%d patients with all data available" % (sum(mask_keep), mask_keep.shape[0]))
    df_cln = df_cln.loc[mask_keep].copy()

    # process alterations
    df_alt = process_alt(df_alt, df_cna, df_mut, samples_all, col_sam_id, col_gen, col_alt, col_alt_det, col_alt_cat,
                         col_alt_cat_sim, col_lvl, keep_alt_det)

    return df_alt, df_cln
