import os
import subprocess

filepath = "imeta_output.tsv"
with open(filepath, "r") as file:
    datasets = [x.strip() for x in file.readlines() if x.strip()!="----"]
    datasets = [x.split("collection:")[1].strip() for x in datasets]

header = ["Dataset", "FASTQ_I_Name", "FASTQ_1_Name", "FASTQ_2_Name"]
print("\t".join(header))

for dataset in datasets:
    output = subprocess.check_output(["ils", "-r", dataset+"/archive"])
    lines = [x.strip().split(":")[0] for x in output.split("\n") if x.strip()!=""]
    # if len(lines)==4:
    #     print("\t".join(lines))
    if len(lines)!=4:
        print("\t".join(lines))
        print(dataset)
