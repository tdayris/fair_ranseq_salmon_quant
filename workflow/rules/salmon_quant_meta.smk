module salmon_tximport:
    meta_wrapper:
        f"{snakemake_wrappers_prefix}/meta/bio/salmon_tximport"
    config:
        config


use rule salmon_quant_reads from salmon_tximport as fair_rnaseq_salmon_quant_salmon_quant_reads with:
    input:
        unpack(get_salmon_quant_reads_input),
    output:
        quant=temp(
            "tmp/fair_rnaseq_salmon_quant_salmon_quant_reads/{species}.{build}.{release}/{sample}/quant.sf"
        ),
        quant_gene=temp(
            "tmp/fair_rnaseq_salmon_quant_salmon_quant_reads/{species}.{build}.{release}/{sample}/quant.genes.sf"
        ),
        lib=temp(
            "tmp/fair_rnaseq_salmon_quant_salmon_quant_reads/{species}.{build}.{release}/{sample}/lib_format_counts.json"
        ),
        aux_info=temp(
            directory(
                "tmp/fair_rnaseq_salmon_quant_salmon_quant_reads/{species}.{build}.{release}/{sample}/aux_info"
            )
        ),
        cmd_info=temp(
            "tmp/fair_rnaseq_salmon_quant_salmon_quant_reads/{species}.{build}.{release}/{sample}/cmd_info.json"
        ),
        libparams=temp(
            directory(
                "tmp/fair_rnaseq_salmon_quant_salmon_quant_reads/{species}.{build}.{release}/{sample}/libParams"
            )
        ),
        logs=temp(
            directory(
                "tmp/fair_rnaseq_salmon_quant_salmon_quant_reads/{species}.{build}.{release}/{sample}/logs"
            )
        ),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: 20 * 1024 * attempt,
        runtime=lambda wildcards, attempt: 45 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_rnaseq_salmon_quant_salmon_quant_pair_ended_reads/{sample}.{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_rnaseq_salmon_quant_salmon_quant_pair_ended_reads/{sample}.{species}.{build}.{release}.tsv"
    params:
        libtype=lookup_config(dpath="params/salmon/libtype", default="A"),
        extra=lookup_config(
            dpath="params/salmon/quant",
            default="--numBootstraps 100 --gcBias --seqBias --posBias",
        ),


use rule tximport from salmon_tximport as fair_rnaseq_salmon_quant_tximport with:
    input:
        quant=expand(
            "tmp/fair_rnaseq_salmon_quant_salmon_quant_reads/{sample.species}.{sample.build}.{sample.release}/{sample.sample_id}/quant.sf",
            sample=lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}'",
                within=samples,
            ),
        ),
        quant_genes=expand(
            "tmp/fair_rnaseq_salmon_quant_salmon_quant_reads/{sample.species}.{sample.build}.{sample.release}/{sample.sample_id}/quant_genes.sf",
            sample=lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}'",
                within=samples,
            ),
        ),
        lib=expand(
            "tmp/fair_rnaseq_salmon_quant_salmon_quant_reads/{sample.species}.{sample.build}.{sample.release}/{sample.sample_id}/lib_format_counts.json",
            sample=lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}'",
                within=samples,
            ),
        ),
        aux_info=expand(
            "tmp/fair_rnaseq_salmon_quant_salmon_quant_reads/{sample.species}.{sample.build}.{sample.release}/{sample.sample_id}/aux_info",
            sample=lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}'",
                within=samples,
            ),
        ),
        cmd_info=expand(
            "tmp/fair_rnaseq_salmon_quant_salmon_quant_reads/{sample.species}.{sample.build}.{sample.release}/{sample.sample_id}/cmd_info.json",
            sample=lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}'",
                within=samples,
            ),
        ),
        libparams=expand(
            "tmp/fair_rnaseq_salmon_quant_salmon_quant_reads/{sample.species}.{sample.build}.{sample.release}/{sample.sample_id}/libParams",
            sample=lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}'",
                within=samples,
            ),
        ),
        logs=expand(
            "tmp/fair_rnaseq_salmon_quant_salmon_quant_reads/{sample.species}.{sample.build}.{sample.release}/{sample.sample_id}/logs",
            sample=lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}'",
                within=samples,
            ),
        ),
        tx_to_gene=expand(
            "reference/annotation/{genome.species}.{genome.build}.{genome.release}/{genome.species}.{genome.build}.{genome.release}.id_to_gene.tsv",
            genome=lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}'",
                within=genomes,
            ),
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 8,
        runtime=lambda wildcards, attempt: attempt * 45,
        tmpdir=tmp,
    output:
        txi=temp(
            "tmp/fair_rnaseq_salmon_quant_tximport/{species}.{build}.{release}/SummarizedExperimentObject.RDS"
        ),
    params:
        extra=lookup_config(
            dpath="params/fair_rnaseq_salmon_quant_tximport",
            default="type='salmon', ignoreTxVersion=TRUE, ignoreAfterBar=TRUE",
        ),
    log:
        "logs/fair_rnaseq_salmon_quant_tximport/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_rnaseq_salmon_quant_tximport/{species}.{build}.{release}.tsv"
