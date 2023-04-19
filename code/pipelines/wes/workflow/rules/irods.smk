rule irods_check_status:
    log:
        "%s/irods/irods_check_status.log" % L_FOLDER
    conda:
        "../envs/python.yaml"
    input:
        table = config["samples"],
        imeta = config["imeta"]
    output:
        status = "%s/irods/irods_check_status.tsv" % R_FOLDER
    resources:
        queue = "shortq"
    threads: 1
    shell:
        """
        python workflow/scripts/01.1_irods_check_status.py \
            --table {input.table} \
            --irods {input.imeta} \
            --output {output.status} &> {log}
        """


rule irods_get_fastq:
    log:
        "%s/irods/irods_get_fastq_{sample}.log" % L_FOLDER
    input:
        table = config["samples"]
    params:
        dir = "%s/data/fastq" % R_FOLDER,
        irods_1 = lambda w: get_fastqs_irods(w)["r1"],
        irods_2 = lambda w: get_fastqs_irods(w)["r2"],
        local_1 = lambda w: get_fastqs_local(w)["r1"],
        local_2 = lambda w: get_fastqs_local(w)["r2"],
        name_1 = lambda w: get_fastqs_names(w)["r1"],
        name_2 = lambda w: get_fastqs_names(w)["r2"]
    output:
        r1=temp("%s/data/fastq/{sample}_R1.fastq.gz" % R_FOLDER),
        r2=temp("%s/data/fastq/{sample}_R2.fastq.gz" % R_FOLDER)
    shell:
        """
        module unload anaconda3
        bash workflow/scripts/01.2_irods_get_fastq.sh \
            -d {params.dir} \
            -a {params.irods_1} \
            -b {params.irods_2} \
            -x {params.local_1} \
            -y {params.local_2} \
            -m {params.name_1} \
            -n {params.name_2} &> {log}
        """
