module salmon_tximport:
    meta_wrapper:
        "v3.3.3/meta/bio/salmon_tximport"
    config:
        config


use rule salmon_decoy_sequences from salmon_tximport with:
    input:
        transcriptome="reference/{species}.{build}.{release}.cdna.fasta",
        genome="reference/{species}.{build}.{release}.dna.fasta",
    output:
        gentrome=temp("reference/{species}.{build}.{release}.gentrome.fasta"),
        decoys=temp("reference/{species}.{build}.{release}.decoys.txt"),
    log:
        "logs/salmon/decoy_sequence/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/salmon/decoy_sequence/{species}.{build}.{release}.tsv"


use rule salmon_index_gentrome from salmon_tximport with:
    input:
        sequences="reference/{species}.{build}.{release}.gentrome.fasta",
        decoys="reference/{species}.{build}.{release}.decoys.txt",
    output:
        multiext(
            "reference/{species}.{build}.{release}/salmon_index_{species}.{build}.{release}",
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
        ),
    log:
        "logs/salmon/index/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/salmon/index/{species}.{build}.{release}.tsv"
    params:
        extra=config.get("params", {}).get("salmon", {}).get("index", ""),


use rule salmon_quant_reads from salmon_tximport with:
    input:
        unpack(get_salmon_quant_reads_input),
        r="tmp/fastp/trimmed/{sample}.{stream}.fastq",
        index=multiext(
            "reference/{species}.{build}.{release}/salmon_index_{species}.{build}.{release}",
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
        ),
        gtf="resources/{species}.{build}.{release}.gtf",
    output:
        quant=temp("tmp/salmon/quant/{species}.{build}.{release}/{sample}/quant.sf"),
        quant_gene=temp(
            "tmp/salmon/quant/{species}.{build}.{release}/{sample}/quant.genes.sf"
        ),
        lib=temp(
            "tmp/salmon/quant/{species}.{build}.{release}/{sample}/lib_format_counts.json"
        ),
        aux_info=temp(
            directory("tmp/salmon/quant/{species}.{build}.{release}/{sample}/aux_info")
        ),
        cmd_info=temp(
            "tmp/salmon/quant/{species}.{build}.{release}/{sample}/cmd_info.json"
        ),
        libparams=temp(
            directory(
                "tmp/salmon/quant/{species}.{build}.{release}/{sample}/libParams"
            )
        ),
        logs=temp(
            directory("tmp/salmon/quant/{species}.{build}.{release}/{sample}/logs")
        ),
    log:
        "logs/salmon/quant/{sample}.log",
    benchmark:
        "benchmark/salmon/quant/{sample}.tsv"
    params:
        libtype="A",
        extra=config.get("params", {}).get("salmon", {}).get("quant", ""),


rule tximport:
    input:
        quant=expand(
            "tmp/salmon/quant/{species}.{build}.{release}/{sample}/quant.sf",
            sample=samples.sample_id,
        ),
        lib=expand(
            "tmp/salmon/quant/{species}.{build}.{release}/{sample}/lib_format_counts.json",
            sample=samples.sample_id,
        ),
        aux_info=expand(
            "tmp/salmon/quant/{species}.{build}.{release}/{sample}/aux_info",
            sample=samples.sample_id,
        ),
        cmd_info=expand(
            "tmp/salmon/quant/{species}.{build}.{release}/{sample}/cmd_info.json",
            sample=samples.sample_id,
        ),
        libparams=expand(
            "tmp/salmon/quant/{species}.{build}.{release}/{sample}/libParams",
            sample=samples.sample_id,
        ),
        logs=expand(
            "tmp/salmon/quant/{species}.{build}.{release}/{sample}/logs",
            sample=samples.sample_id,
        ),
        tx_to_gene="resources/{species}.{build}.{release}/tx2gene.tsv",
    output:
        txi=temp(
            "tmp/tximport/{species}.{build}.{release}/SummarizedExperimentObject.RDS"
        ),
    params:
        extra="type='salmon'",
    log:
        "logs/tximport/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/tximport/{species}.{build}.{release}.tsv"
