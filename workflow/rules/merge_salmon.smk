rule aggregate_salmon_counts:
    input:
        quant=expand(
            "tmp/salmon/quant/{species}.{build}.{release}/{sample}/quant.genes.sf",
            sample=samples.sample_id,
        ),
        tx2gene="resources/{species}.{build}.{release}/tx2gene.tsv",
    output:
        tsv=protected(
            "results/{species}.{build}.{release}/Quantification/{counts}.{targets}.tsv"
        ),
    log:
        "logs/python/aggregate_salmon/{species}.{build}.{release}/{counts}.{targets}.log",
    log:
        "benchmark/python/aggregate_salmon/{species}.{build}.{release}/{counts}.{targets}.tsv",
    threads: 1
    params:
        header=False,
        position=False,
        gencode=True,
        genes=lambda wildcards: str(wildcards.counts).lower().startswith("gene"),
        index_label=True,
        fillna=0,
        column=lambda wildcards: "NumReads" if str(wildcards.counts).lower() == "raw" else "TPM",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/merge_salmon_quant.py"