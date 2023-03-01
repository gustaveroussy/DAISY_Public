import argparse
import numpy as np
import pandas as pd
import re
import functools
print = functools.partial(print, flush=True)

def print_status_changes(df_table):
    old_ok = ~df_table["Status_Old"].isnull()
    new_ok = ~df_table["Status"].isnull()
    total = df_table.shape[0]
    changed_to_ok = sum(~old_ok & new_ok)
    changed_to_nok = sum(old_ok & ~new_ok)
    print("%d/%d samples changed from Not OK to OK" % (changed_to_ok, total))
    print("%d/%d samples changed from OK to Not OK" % (changed_to_nok, total))


def main(args):
    df_table = pd.read_table(args.table)
    if "Status" in df_table:
        df_table = df_table.rename(columns={"Status": "Status_Old"})

    df_irods = pd.read_table(args.irods, header=None)
    df_irods.columns = ["Dataset"]
    mask_keep = df_irods["Dataset"].apply(lambda x: x.startswith("collection:"))
    df_irods = df_irods.loc[mask_keep].copy()
    df_irods["Dataset"] = df_irods["Dataset"].apply(lambda x: re.sub("collection: ", "", x))
    df_irods["Status"] = "OK IRODS"

    col_x = "FASTQ_1"
    col_y = "Dataset"
    df_table[col_y] = df_table[col_x].apply(lambda x: "/".join(x.split("/")[:-2]))
    df_table = df_table.merge(df_irods, how="left", on="Dataset")

    if "Status_Old" in df_table:
        print_status_changes(df_table)
        del df_table["Status_Old"]

    df_table.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Create a status indicator indicating whether the dataset is still
                                     available in IRODS""")
    parser.add_argument("--table", type=str, help="Path to table of samples")
    parser.add_argument("--irods", type=str, help="Path to output of imeta ls -C command.")
    parser.add_argument("--output", type=str, help="Path to output table with status.")

    args = parser.parse_args()
    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
