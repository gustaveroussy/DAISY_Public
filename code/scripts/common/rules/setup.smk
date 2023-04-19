rule setup_conda:
    conda:
        "../envs/DAISY.yaml"
    params:
        comut = "../common/external/comut",
    output:
        touch("../common/logs/setup_conda.done")
    shell:
        """
        pip install -e {params.comut}
        """
