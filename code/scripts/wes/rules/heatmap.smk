rule draw_heatmap_all_alterations:
    log:
        "%s/draw_heatmap_all_alterations_{direction}.log" % L_FOLDER
    input:
        alterations = "%s/alterations/aggregated_alterations.tsv" % R_FOLDER,
        samples = expand("%s/selection/selection_samples_{cohort}.tsv" % R_FOLDER,
            cohort=config["data"]["cohorts"]),
        counts = "%s/selection/selection_tumor_types.tsv" % R_FOLDER,
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    output:
        tables = "%s/heatmap_all/tables_heatmap_{direction}.xlsx" % R_FOLDER,
        plot = "%s/heatmap_all/heatmap_{direction}.svg" % R_FOLDER,
        plot_paper = "%s/F5_{direction}.svg" % F_FOLDER,
    params:
        width = lambda w: config["heatmap_all"][w.direction]["width"],
        height = lambda w: config["heatmap_all"][w.direction]["height"]
    resources:
        partition = "cpu_med",
        mem_mb = 8000,
        time = "01:00:00"
    threads: 4
    shell:
        """
        Rscript workflow/scripts/02.1_draw_heatmap_all_alterations.R \
            --alterations {input.alterations} \
            --samples {input.samples} \
            --counts {input.counts} \
            --direction {wildcards.direction} \
            --n_cores {threads} \
            --output_tables {output.tables} \
            --output_plot_width {params.width} \
            --output_plot_height {params.height} \
            --output_plot {output.plot} \
            --output_plot_paper {output.plot_paper} \
            --log {log}
        """
