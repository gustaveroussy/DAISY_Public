rule draw_oncoplot_samples:
    input:
        config["data"]["inputs"]["cln_curated"],
        config["data"]["outputs"]["alterations"],
        config["data"]["inputs"]["sam_used"],
        config["data"]["inputs"]["ids_anonym"],
        config["setup"]["DAISY_done"]
    output:
        "%s/oncoplot/oncoplot_samples_{select}.pdf" % R_FOLDER,
        "%s/oncoplot/tests_response_{select}.xlsx" % R_FOLDER
    benchmark:
        "%s/draw_oncoplot_samples_{select}.tsv" % B_FOLDER
    log:
        "%s/draw_oncoplot_samples_{select}.log" % L_FOLDER
    conda:
        config["setup"]["DAISY"]
    params:
        threshold_alt = 0.03
    threads: 1
    shell:
        """
        python -u scripts/02.1_draw_oncoplot_samples.py \
            --cln {input[0]} \
            --alt {input[1]} \
            --sam {input[2]} \
            --ids_anonym {input[3]} \
            --threshold_alt {params.threshold_alt} \
            --samples_select {wildcards.select} \
            --output {output[0]} \
            --output_tests {output[1]} &> {log}
        """


rule draw_oncoplot_subject:
    wildcard_constraints:
        mode_cna = "wo_cna|w_cna"
    input:
        config["data"]["inputs"]["cln_curated"],
        config["data"]["outputs"]["alterations"],
        config["data"]["inputs"]["mut_pass"],
        config["data"]["inputs"]["cna_pass"],
        config["setup"]["DAISY_done"]
    output:
        "%s/oncoplot/oncoplot_subject_{subject}_{mode_cna}.pdf" % R_FOLDER
    benchmark:
        "%s/draw_oncoplot_subject_{subject}_{mode_cna}.tsv" % B_FOLDER
    log:
        "%s/draw_oncoplot_subject_{subject}_{mode_cna}.log" % L_FOLDER
    conda:
        config["setup"]["DAISY"]
    threads: 1
    shell:
        """
        python -u scripts/02.2_draw_oncoplot_subject.py \
            --cln {input[0]} \
            --alt {input[1]} \
            --mut {input[2]} \
            --cna {input[3]} \
            --mode_cna {wildcards.mode_cna} \
            --subject {wildcards.subject} \
            --output {output} &> {log}
        """


rule draw_oncoplot_t1t2:
    wildcard_constraints:
        alt_lvl = "gen|pth"
    input:
        config["data"]["inputs"]["cln_curated"],
        config["data"]["outputs"]["alterations"],
        config["data"]["inputs"]["mut_pass"],
        config["data"]["inputs"]["cna_pass"],
        config["data"]["inputs"]["cln_tcga_curated"],
        config["data"]["inputs"]["mut_tcga_pass"],
        config["data"]["inputs"]["cna_tcga_pass"],
        config["data"]["inputs"]["pathways"],
        config["data"]["inputs"]["ids_anonym"],
        config["setup"]["DAISY_done"]
    output:
        "%s/oncoplot/oncoplot_t1t2_diff_per_{alt_lvl}_{agg_mod}.pdf" % R_FOLDER
    benchmark:
        "%s/draw_oncoplot_t1t2_{alt_lvl}_{agg_mod}.tsv" % B_FOLDER
    log:
        "%s/draw_oncoplot_t1t2_{alt_lvl}_{agg_mod}.log" % L_FOLDER
    conda:
        config["setup"]["DAISY"]
    threads: 1
    shell:
        """
        python -u scripts/02.3_draw_oncoplot_t1t2.py \
            --cln {input[0]} \
            --alt {input[1]} \
            --mut {input[2]} \
            --cna {input[3]} \
            --cln_tcga {input[4]} \
            --mut_tcga {input[5]} \
            --cna_tcga {input[6]} \
            --pth {input[7]} \
            --ids_anonym {input[8]} \
            --alt_lvl {wildcards.alt_lvl} \
            --agg_mod {wildcards.agg_mod} \
            --output {output} &> {log}
        """
