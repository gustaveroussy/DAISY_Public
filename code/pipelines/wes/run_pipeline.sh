#!/bin/bash

n_splits=$1

for ((i=1; i<=${n_splits}; i+=1));
do
    printf "\n\nRunning batch $i/$n_splits ...\n\n"
    cp config/samples_${i}.tsv config/samples.tsv
    snakemake -s workflow/Snakefile --profile slurm --jobs 100 --rerun-incomplete
done
