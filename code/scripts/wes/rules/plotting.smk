rule draw_heatmap_all_alterations:
    input:
        config["data"]["inputs"]["cln_curated"],
        config["data"]["outputs"]["alterations"],
        config["setup"]["DAISY_done"]
    output:
        "%s/heatmap_all/tables_heatmap_{direction}.xlsx" % R_FOLDER,
        "%s/heatmap_all/heatmap_{direction}.pdf" % R_FOLDER
    benchmark:
        "%s/draw_heatmap_all_alterations_{direction}.tsv" % B_FOLDER
    log:
        "%s/draw_heatmap_all_alterations_{direction}.log" % L_FOLDER
    conda:
        config["setup"]["DAISY"]
    params:
        cohorts = ["daisy"],
        barplot_height = lambda w: config["heatmap_all"][w.direction]["barplot_height"],
        width = lambda w: config["heatmap_all"][w.direction]["width"],
        height = lambda w: config["heatmap_all"][w.direction]["height"]
    threads: 4
    shell:
        """
        Rscript scripts/02.1_draw_heatmap_all_alterations.R \
            --cohorts {params.cohorts} \
            --clns {input[0]} \
            --alts {input[1]} \
            --direction {wildcards.direction} \
            --n_cores {threads} \
            --output_tables {output[0]} \
            --output_plot_barplot_height {params.barplot_height} \
            --output_plot_width {params.width} \
            --output_plot_height {params.height} \
            --output_plot {output[1]} \
            --log {log}
        """


rule draw_oncoplot_samples:
    input:
        config["data"]["inputs"]["cln_curated"],
        config["data"]["outputs"]["alterations"],
        config["data"]["inputs"]["sam_used"],
        config["data"]["inputs"]["ids_anonym"],
        config["setup"]["DAISY_done"]
    output:
        "%s/oncoplot/oncoplot_samples_{select}.pdf" % R_FOLDER,
    benchmark:
        "%s/draw_oncoplot_samples_{select}.tsv" % B_FOLDER
    log:
        "%s/draw_oncoplot_samples_{select}.log" % L_FOLDER
    conda:
        config["setup"]["DAISY"]
    params:
        threshold_alt = 0.03,
        output_tests = lambda w: "%s/oncoplot/tests_response_{select}.xlsx" % R_FOLDER if w.select=="t1" else []
    threads: 1
    shell:
        """
        python -u scripts/02.2_draw_oncoplot_samples.py \
            --cln {input[0]} \
            --alt {input[1]} \
            --sam {input[2]} \
            --ids_anonym {input[3]} \
            --threshold_alt {params.threshold_alt} \
            --samples_select {wildcards.select} \
            --output {output} \
            --output_tests {params.output_tests} &> {log}
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
        python -u scripts/02.3_draw_oncoplot_subject.py \
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
        "%s/oncoplot/oncoplot_t1t2_diff_per_{alt_lvl}_{agg_mod}_{sel_mod}_{vaf_inc}.pdf" % R_FOLDER
    benchmark:
        "%s/draw_oncoplot_t1t2_{alt_lvl}_{agg_mod}_{sel_mod}_{vaf_inc}.tsv" % B_FOLDER
    log:
        "%s/draw_oncoplot_t1t2_{alt_lvl}_{agg_mod}_{sel_mod}_{vaf_inc}.log" % L_FOLDER
    conda:
        config["setup"]["DAISY"]
    threads: 1
    shell:
        """
        python -u scripts/02.4_draw_oncoplot_t1t2.py \
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
            --sel_mod {wildcards.sel_mod} \
            --vaf_inc {wildcards.vaf_inc} \
            --output {output} &> {log}
        """
