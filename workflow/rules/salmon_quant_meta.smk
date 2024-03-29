module salmon_tximport:
    meta_wrapper:
        "v3.7.0/meta/bio/salmon_tximport"
    config:
        config


use rule salmon_decoy_sequences from salmon_tximport as fair_rnaseq_salmon_quant_salmon_decoy_sequences with:
    input:
        transcriptome=getattr(
            lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=genomes,
            ),
            "cdna_fasta",
            "reference/sequences/{species}.{build}.{release}.transcripts.fasta",
        ),
        genome=getattr(
            lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=genomes,
            ),
            "dna_fasta",
            "reference/sequences/{species}.{build}.{release}.dna.fasta",
        ),
    output:
        gentrome=temp("reference/sequences/{species}.{build}.{release}.gentrome.fasta"),
        decoys=temp("reference/sequences/{species}.{build}.{release}.decoys.txt"),
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: 512 * attempt,
        runtime=lambda wildcards, attempt: 25 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_rnaseq_salmon_quant/salmon_decoy_sequences/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_rnaseq_salmon_quant/salmon_decoy_sequences/{species}.{build}.{release}.tsv"


use rule salmon_index_gentrome from salmon_tximport as fair_rnaseq_salmon_quant_salmon_index_gentrome with:
    input:
        sequences="reference/sequences/{species}.{build}.{release}.gentrome.fasta",
        decoys="reference/sequences/{species}.{build}.{release}.decoys.txt",
    output:
        temp(
            multiext(
                "reference/salmon_index/{species}.{build}.{release}/{species}.{build}.{release}/",
                "complete_ref_lens.bin",
                "ctable.bin",
                "ctg_offsets.bin",
                "duplicate_clusters.tsv",
                "info.json",
                "mphf.bin",
                "pos.bin",
                "pre_indexing.log",
                "rank.bin",
                "refAccumLengths.bin",
                "ref_indexing.log",
                "reflengths.bin",
                "refseq.bin",
                "seq.bin",
                "versionInfo.json",
            )
        ),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: 20 * 1024 * attempt,
        runtime=lambda wildcards, attempt: 45 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_rnaseq_salmon_quant/salmon_index_gentrome/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_rnaseq_salmon_quant/salmon_index_gentrome/{species}.{build}.{release}.tsv"
    params:
        extra=lookup(dpath="params/salmon/index", within=config),


use rule salmon_quant_reads from salmon_tximport as fair_rnaseq_salmon_quant_salmon_quant_reads with:
    input:
        unpack(get_salmon_quant_reads_input),
    output:
        quant=temp(
            "tmp/fair_rnaseq_salmon_quant/salmon_quant_reads/{species}.{build}.{release}/{sample}/quant.sf"
        ),
        quant_gene=temp(
            "tmp/fair_rnaseq_salmon_quant/salmon_quant_reads/{species}.{build}.{release}/{sample}/quant.genes.sf"
        ),
        lib=temp(
            "tmp/fair_rnaseq_salmon_quant/salmon_quant_reads/{species}.{build}.{release}/{sample}/lib_format_counts.json"
        ),
        aux_info=temp(
            directory(
                "tmp/fair_rnaseq_salmon_quant/salmon_quant_reads/{species}.{build}.{release}/{sample}/aux_info"
            )
        ),
        cmd_info=temp(
            "tmp/fair_rnaseq_salmon_quant/salmon_quant_reads/{species}.{build}.{release}/{sample}/cmd_info.json"
        ),
        libparams=temp(
            directory(
                "tmp/fair_rnaseq_salmon_quant/salmon_quant_reads/{species}.{build}.{release}/{sample}/libParams"
            )
        ),
        logs=temp(
            directory(
                "tmp/fair_rnaseq_salmon_quant/salmon_quant_reads/{species}.{build}.{release}/{sample}/logs"
            )
        ),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: 20 * 1024 * attempt,
        runtime=lambda wildcards, attempt: 45 * attempt,
        tmpdir=tmp,
    log:
        "logs/fair_rnaseq_salmon_quant/salmon_quant_pair_ended_reads/{sample}.{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_rnaseq_salmon_quant/salmon_quant_pair_ended_reads/{sample}.{species}.{build}.{release}.tsv"
    params:
        libtype=lookup(dpath="params/salmon/libtype", within=config),
        extra=lookup(dpath="params/salmon/quant", within=config),


use rule tximport from salmon_tximport as fair_rnaseq_salmon_quant_tximport with:
    input:
        quant=expand(
            "tmp/fair_rnaseq_salmon_quant/salmon_quant_reads/{sample.species}.{sample.build}.{sample.release}/{sample.sample_id}/quant.sf",
            sample=lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}'",
                within=samples,
            ),
        ),
        quant_genes=expand(
            "tmp/fair_rnaseq_salmon_quant/salmon_quant_reads/{sample.species}.{sample.build}.{sample.release}/{sample.sample_id}/quant_genes.sf",
            sample=lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}'",
                within=samples,
            ),
        ),
        lib=expand(
            "tmp/fair_rnaseq_salmon_quant/salmon_quant_reads/{sample.species}.{sample.build}.{sample.release}/{sample.sample_id}/lib_format_counts.json",
            sample=lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}'",
                within=samples,
            ),
        ),
        aux_info=expand(
            "tmp/fair_rnaseq_salmon_quant/salmon_quant_reads/{sample.species}.{sample.build}.{sample.release}/{sample.sample_id}/aux_info",
            sample=lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}'",
                within=samples,
            ),
        ),
        cmd_info=expand(
            "tmp/fair_rnaseq_salmon_quant/salmon_quant_reads/{sample.species}.{sample.build}.{sample.release}/{sample.sample_id}/cmd_info.json",
            sample=lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}'",
                within=samples,
            ),
        ),
        libparams=expand(
            "tmp/fair_rnaseq_salmon_quant/salmon_quant_reads/{sample.species}.{sample.build}.{sample.release}/{sample.sample_id}/libParams",
            sample=lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}'",
                within=samples,
            ),
        ),
        logs=expand(
            "tmp/fair_rnaseq_salmon_quant/salmon_quant_reads/{sample.species}.{sample.build}.{sample.release}/{sample.sample_id}/logs",
            sample=lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}'",
                within=samples,
            ),
        ),
        tx_to_gene=expand(
            "reference/annotation/{genome.species}.{genome.build}.{genome.release}.id_to_gene.tsv",
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
            "tmp/fair_rnaseq_salmon_quant/tximport/{species}.{build}.{release}/SummarizedExperimentObject.RDS"
        ),
    params:
        extra=lookup(dpath="params/tximport", within=config),
    log:
        "logs/fair_rnaseq_salmon_quant/tximport/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_rnaseq_salmon_quant/tximport/{species}.{build}.{release}.tsv"
