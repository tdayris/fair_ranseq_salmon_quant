rule fair_rnaseq_salmon_quant_in_house_pca:
    input:
        quant="results/{species}.{build}.{release}/Quantification/TPM.genes.tsv",
    output:
        pca_yaml=temp(
            "tmp/fair_rnaseq_salmon_quant_in_house_pca/{species}.{build}.{release}/pca.yaml"
        ),
        corr_yaml=temp(
            "tmp/fair_rnaseq_salmon_quant_in_house_pca/{species}.{build}.{release}/correlation.yaml"
        ),
        pca_png=report(
            "results/{species}.{build}.{release}/Quantification/PCA.png",
            caption="../report/in_house_pca.rst",
            category="Quantification",
            labels={
                "species": "{species}.{build}.{release}",
            },
        ),
        scree_png=report(
            "results/{species}.{build}.{release}/Quantification/Scree.png",
            caption="../report/in_house_scree.rst",
            category="Quantification",
            labels={
                "species": "{species}.{build}.{release}",
            },
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_rnaseq_salmon_quant_in_house_pca/{species}.{build}.{release}/pca.log",
    benchmark:
        "benchmark/fair_rnaseq_salmon_quant_in_house_pca/{species}.{build}.{release}/pca.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/fair_rnaseq_salmon_quant_in_house_pca.py"
