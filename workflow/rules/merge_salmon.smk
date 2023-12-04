rule aggregate_raw_gene_counts:
    input:
        quant=expand(
            "tmp/salmon/quant/{species}.{build}.{release}/{sample}/quant.genes.sf",
            sample=samples.sample_id,
        ),
        tx2gene="resources/{species}.{build}.{release}/tx2gene.tsv",
    output:
        tsv=protected("results/{species}.{build}.{release}/Quantification/Raw.genes.tsv"),
    log:
        "logs/python/aggregate_salmon/{species}.{build}.{release}/genes.raw.log",
    log:
        "benchmark/python/aggregate_salmon/{species}.{build}.{release}/genes.raw.tsv",
    threads: 1
    params:
        header=False,
        position=False,
        gencode=True,
        genes=True,
        index_label=True,
        fillna="Unknown",
        column="NumReads",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/merge_salmon_quant.py"


use rule aggregate_raw_gene_counts as aggregate_tpm_gene_counts with:
    input:
        quant=expand(
            "tmp/salmon/quant/{species}.{build}.{release}/{sample}/quant.genes.sf",
            sample=samples.sample_id,
        ),
        tx2gene="resources/{species}.{build}.{release}/tx2gene.tsv",
    output:
        tsv=protected("results/{species}.{build}.{release}/Quantification/TPM.genes.tsv"),
    log:
        "logs/python/aggregate_salmon/{species}.{build}.{release}/genes.tpm.log",
    log:
        "benchmark/python/aggregate_salmon/{species}.{build}.{release}/genes.tpm.tsv",
    threads: 1
    params:
        header=False,
        position=False,
        gencode=True,
        genes=True,
        index_label=True,
        fillna="Unknown",

use rule aggregate_raw_gene_counts as aggregate_tpm_transcripts_counts with:
    input:
        quant=expand(
            "tmp/salmon/quant/{species}.{build}.{release}/{sample}/quant.sf",
            sample=samples.sample_id,
        ),
        tx2gene="resources/{species}.{build}.{release}/tx2gene.tsv",
    output:
        tsv=protected("results/{species}.{build}.{release}/Quantification/TPM.transcripts.tsv"),
    log:
        "logs/python/aggregate_salmon/{species}.{build}.{release}/transcripts.tpm.log",
    log:
        "benchmark/python/aggregate_salmon/{species}.{build}.{release}/transcripts.tpm.tsv",
    threads: 1
    params:
        header=False,
        position=False,
        gencode=True,
        genes=False,
        index_label=True,
        fillna="Unknown",


use rule aggregate_raw_gene_counts as aggregate_raw_transcripts_counts with:
    input:
        quant=expand(
            "tmp/salmon/quant/{species}.{build}.{release}/{sample}/quant.sf",
            sample=samples.sample_id,
        ),
        tx2gene="resources/{species}.{build}.{release}/tx2gene.tsv",
    output:
        tsv=protected("results/{species}.{build}.{release}/Quantification/Raw.transcripts.tsv"),
    log:
        "logs/python/aggregate_salmon/{species}.{build}.{release}/transcripts.raw.log",
    log:
        "benchmark/python/aggregate_salmon/{species}.{build}.{release}/transcripts.raw.tsv",
    threads: 1
    params:
        header=False,
        position=False,
        gencode=True,
        genes=False,
        index_label=True,
        fillna="Unknown",
        column="NumReads",