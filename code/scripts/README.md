# Organisation

The code for each separate part of the analysis is located in the folders listed below along with a brief description of
the analysis performed in each folder.

- `common`
    - Contains scripts and files for setting up conda environments that were used for running the analyses.
    - Hosts the code for useful functions used throughout the code. 
- `slides`
    - Contains the code for the machine-learning pipeline developed to analyze HER2 slides. See `slides/README.md` for
      more details.
- `wes`
    - Contains scripts for preparing table that combine driver point mutations, small indels, and focal CNAs.
    - Hosts the code for drawing Fig. 4A, Fig. 4B.

# Per analysis

## Organisation

The workflow of each analysis (except for `slides`) is divided in rules that are assembled in the snakemake subfiles of
`rules` which are themselves included in the global snakemake file `Snakefile`. The rules call scripts
located in the `scripts` folder.

The parameters of the analysis may be modified in the `config/config.yaml` file. All libraries needed for running the
analysis are installed during the first steps of the pipeline through [conda virtual
environments](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) whose specifications
are given in the folder `common/envs`.

## How to run the analysis?

In order to run the analysis, ensure the `snakemake` and `conda` commands are available (if you are on a HPC using
slurm, you may use `module load` to load [modules](https://curc.readthedocs.io/en/latest/compute/modules.html) providing
you with these commands). You can launch the full pipeline via

```
snakemake --jobs [n_jobs] --profile [your_profile]
```

where `[n_jobs]` is the maximum number of CPUs (or jobs if you are on a cluster) to be used/run in parallel and
`[your_profile]` your snakemake profile (read the [snakemake
documentation](<https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles>)).

In case you are only interested in running part of the analysis so as to reproduce one of the results file, you may do
so with

```
snakemake --jobs [n_jobs] --profile [your_profile] [path/to/file/you/want/to/reproduce]
```
