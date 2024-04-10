rule fair_rnaseq_salmon_quant_fastp_trimming_pair_ended:
    input:
        sample=expand(
            "tmp/fair_fastqc_multiqc/link_or_concat_pair_ended_input/{sample}.{stream}.fastq.gz",
            stream=stream_list,
            allow_missing=True,
        ),
    output:
        trimmed=temp(
            expand(
                "tmp/fair_rnaseq_salmon_quant/fastp_trimming_pair_ended/{sample}.{stream}.fastq.gz",
                stream=stream_list,
                allow_missing=True,
            )
        ),
        html=report(
            "results/QC/report_pe/{sample}.html",
            caption="../report/fastp.rst",
            category="Quality Controls",
            subcategory="Trimming",
        ),
        json=temp(
            "tmp/fair_rnaseq_salmon_quant/fastp_trimming_pair_ended/{sample}.json"
        ),
    threads: 5
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 1.5,
        runtime=lambda wildcards, attempt: attempt * 60,
        tmpdir=tmp,
    log:
        "logs/fair_rnaseq_salmon_quant/fastp_trimming_pair_ended/{sample}.log",
    benchmark:
        "benchmark/fair_rnaseq_salmon_quant/fastp_trimming_pair_ended/{sample}.tsv"
    params:
        adapters=lookup(dpath="params/fastp/adapters", within=config),
        extra=lookup(dpath="params/fastp/extra", within=config),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/fastp"


use rule fair_rnaseq_salmon_quant_fastp_trimming_pair_ended as fair_rnaseq_salmon_quant_fastp_trimming_single_ended with:
    input:
        sample=expand(
            "tmp/fair_fastqc_multiqc/link_or_concat_single_ended_input/{sample}.fastq.gz",
            allow_missing=True,
        ),
    output:
        trimmed=temp(
            "tmp/fair_rnaseq_salmon_quant/fastp_trimming_single_ended/{sample}.fastq.gz"
        ),
        html=temp(
            "tmp/fair_rnaseq_salmon_quant/fastp_trimming_single_ended/{sample}.html"
        ),
        json=temp(
            "tmp/fair_rnaseq_salmon_quant/fastp_trimming_single_ended/{sample}.json"
        ),
    log:
        "logs/fair_rnaseq_salmon_quant/fastp_trimming_single_ended/{sample}.log",
    benchmark:
        "benchmark/fair_rnaseq_salmon_quant/fastp_trimming_single_ended/{sample}.tsv"
