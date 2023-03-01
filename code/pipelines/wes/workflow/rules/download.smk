####
#### TCGA specific ####
####

rule download_somatic_mc3:
    input:
        client=config["download_gdc"]["client"],
        manifest=config["download_gdc"]["manifests"]["somatic_mc3"],
        token=config["download_gdc"]["token"]
    output:
        expand("%s/data/somatic_mc3/{id}" % R_FOLDER, id=ids_gdc_somatic_mc3)
    benchmark:
        "%s/download/download_somatic_mc3/mc3.tsv" % B_FOLDER
    log:
        "%s/download/download_somatic_mc3/mc3.log" % L_FOLDER
    params:
        dir="%s/data/somatic_mc3" % R_FOLDER
    threads: 16
    resources:
        queue="shortq",
        mem_mb=24000,
        time_min=240
    shell:
        """
        ./{input.client} download \
            --manifest {input.manifest} \
            --token-file {input.token} \
            --dir {params.dir} \
            --n-processes {threads} \
            --log-file {log}
        """
