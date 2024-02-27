rule fair_rnaseq_salmon_quant_aggregate_salmon_gene_counts:
    input:
        quant=expand(
            "tmp/fair_rnaseq_salmon_quant/salmon_quant_reads/{sample.species}.{sample.build}.{sample.release}/{sample.sample_id}/quant.genes.sf",
            sample=lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}'",
                within=samples,
            ),
        ),
        tx2gene=getattr(
            lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=genomes,
            ),
            "id_to_gene",
            "reference/annotation/{species}.{build}.{release}.id_to_gene.tsv",
        ),
    output:
        tsv=protected(
            "results/{species}.{build}.{release}/Quantification/{counts}.genes.tsv"
        ),
    log:
        "logs/fair_rnaseq_salmon_quant/aggregate_salmon_gene_counts/{species}.{build}.{release}/{counts}.genes.log",
    benchmark:
        "benchmark/fair_rnaseq_salmon_quant/aggregate_salmon_gene_counts/{species}.{build}.{release}/{counts}.genes.tsv"
    threads: 1
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
        "../scripts/merge_salmon_quant.py"


use rule fair_rnaseq_salmon_quant_aggregate_salmon_gene_counts as fair_rnaseq_salmon_quant_aggregate_salmon_transcript_counts with:
    input:
        quant=expand(
            "tmp/fair_rnaseq_salmon_quant/salmon_quant_reads/{sample.species}.{sample.build}.{sample.release}/{sample.sample_id}/quant.sf",
            sample=lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}'",
                within=samples,
            ),
        ),
        tx2gene=getattr(
            lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=genomes,
            ),
            "id_to_gene",
            "reference/annotation/{species}.{build}.{release}.t2g.tsv",
        ),
    output:
        tsv=protected(
            "results/{species}.{build}.{release}/Quantification/{counts}.transcripts.tsv"
        ),
    log:
        "logs/fair_rnaseq_salmon_quant/aggregate_salmon_transcript_counts/{species}.{build}.{release}.{counts}.log",
    benchmark:
        "benchmark/fair_rnaseq_salmon_quant/aggregate_salmon_transcript_counts/{species}.{build}.{release}.{counts}.tsv"
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
