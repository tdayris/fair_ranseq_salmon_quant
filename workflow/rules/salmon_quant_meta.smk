module salmon_tximport:
    meta_wrapper:
        "v3.3.3/meta/bio/salmon_tximport"
    config:
        config


use rule salmon_decoy_sequences from salmon_tximport with:
    input:
        unpack(get_salmon_decoy_sequences_input),
    output:
        gentrome=temp("reference/sequences/{species}.{build}.{release}.gentrome.fasta"),
        decoys=temp("reference/sequences/{species}.{build}.{release}.decoys.txt"),
    log:
        "logs/salmon/decoy_sequence/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/salmon/decoy_sequence/{species}.{build}.{release}.tsv"


use rule salmon_index_gentrome from salmon_tximport with:
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
    log:
        "logs/salmon/index/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/salmon/index/{species}.{build}.{release}.tsv"
    params:
        extra=config.get("params", {}).get("salmon", {}).get("index", ""),


use rule salmon_quant_reads from salmon_tximport with:
    input:
        unpack(get_salmon_quant_reads_input),
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
        "logs/salmon/quant/{sample}.{species}.{build}.{release}.log",
    benchmark:
        "benchmark/salmon/quant/{sample}.{species}.{build}.{release}.tsv"
    params:
        libtype="A",
        extra=config.get("params", {}).get("salmon", {}).get("quant", ""),


rule tximport:
    input:
        unpack(get_tximport_input),
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
