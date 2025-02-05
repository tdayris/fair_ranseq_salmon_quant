rule fair_rnaseq_salmon_quant_aggregate_salmon_gene_counts:
    input:
        quant=expand(
            "tmp/fair_rnaseq_salmon_quant_salmon_quant_reads/{sample.species}.{sample.build}.{sample.release}/{sample.sample_id}/quant.genes.sf",
            sample=lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}'",
                within=samples,
            ),
        ),
        tx2gene=lambda wildcards: get_id2gene(wildcards),
    output:
        tsv=protected(
            "results/{species}.{build}.{release}/Quantification/{counts}.genes.tsv"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir=tmp,
    log:
        "logs/fair_rnaseq_salmon_quant_aggregate_salmon_gene_counts/{species}.{build}.{release}/{counts}.genes.log",
    benchmark:
        "benchmark/fair_rnaseq_salmon_quant_aggregate_salmon_gene_counts/{species}.{build}.{release}/{counts}.genes.tsv"
    params:
        header=False,
        position=False,
        gencode=True,
        genes=True,
        index_label=True,
        fillna=0,
        column=lambda wildcards: (
            "NumReads" if str(wildcards.counts).lower() == "raw" else "TPM"
        ),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/fair_rnaseq_salmon_quant_merge_salmon_quant.py"


use rule fair_rnaseq_salmon_quant_aggregate_salmon_gene_counts as fair_rnaseq_salmon_quant_aggregate_salmon_transcript_counts with:
    input:
        quant=expand(
            "tmp/fair_rnaseq_salmon_quant_salmon_quant_reads/{sample.species}.{sample.build}.{sample.release}/{sample.sample_id}/quant.sf",
            sample=lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}'",
                within=samples,
            ),
        ),
        tx2gene=lambda wildcards: get_tx2gene(wildcards),
    output:
        tsv=protected(
            "results/{species}.{build}.{release}/Quantification/{counts}.transcripts.tsv"
        ),
    log:
        "logs/fair_rnaseq_salmon_quant_aggregate_salmon_transcript_counts/{species}.{build}.{release}.{counts}.log",
    benchmark:
        "benchmark/fair_rnaseq_salmon_quant_aggregate_salmon_transcript_counts/{species}.{build}.{release}.{counts}.tsv"
    params:
        header=False,
        position=False,
        gencode=True,
        genes=False,
        index_label=True,
        fillna=0,
        column=lambda wildcards: (
            "NumReads" if str(wildcards.counts).lower() == "raw" else "TPM"
        ),
