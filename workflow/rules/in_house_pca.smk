rule fair_rnaseq_salmon_quant_bioinfokit_pca:
    input:
        quant="results/{species}.{build}.{release}/Quantification/TPM.transcripts.tsv",
    output:
        loadings_correlation=report(
            "results/{species}.{build}.{release}/Quantification/PCA/Loadings_Correlation_heatmap.png",
            caption="../report/loadings_correlation_heatmap.rst",
            category="Quantification",
            subcategory="PCA",
        ),
        scree_plot=report(
            "results/{species}.{build}.{release}/Quantification/PCA/Screeplot.png",
            caption="../report/screeplot.rst",
            category="Quantification",
            subcategory="PCA",
        ),
        pca_2d=report(
            "results/{species}.{build}.{release}/Quantification/PCA/PCA_2d.png",
            caption="../report/pca_2d.rst",
            category="Quantification",
            subcategory="PCA",
        ),
        pca_3d=report(
            "results/{species}.{build}.{release}/Quantification/PCA/PCA_3d.png",
            caption="../report/pca_3d.rst",
            category="Quantification",
            subcategory="PCA",
        ),
        biplot=report(
            "results/{species}.{build}.{release}/Quantification/PCA/Biplot.png",
            caption="../report/biplot.rst",
            category="Quantification",
            subcategory="PCA",
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_rnaseq_salmon_quant_bioinfokit_pca/{species}.{build}.{release}/pca.log",
    benchmark:
        "benchmark/fair_rnaseq_salmon_quant_bioinfokit_pca/{species}.{build}.{release}/pca.tsv"
    conda:
        "../envs/bioinfokit.yaml"
    script:
        "../scripts/bioinfokit_pca.py"
