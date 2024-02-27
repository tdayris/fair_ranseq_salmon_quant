rule quantification_multiqc_report:
    input:
        fastqc_single_ended=collect(
            "results/QC/report_pe/{sample.sample_id}_fastqc.zip",
            sample=lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}' & downstream_file != downstream_file",
                within=samples,
            ),
        ),
        fastqc_pair_ended=collect(
            "results/QC/report_pe/{sample.sample_id}.{stream}_fastqc.zip",
            sample=lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}' & downstream_file == downstream_file",
                within=samples,
            ),
            stream=stream_list,
        ),
        fastp_pe_json=expand(
            "tmp/fair_rnaseq_salmon_quant/fastp_trimming_pair_ended/{sample.sample_id}.json",
            sample=lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}' & downstream_file == downstream_file",
                within=samples,
            ),
        ),
        fastp_se_json=expand(
            "tmp/fair_rnaseq_salmon_quant/fastp_trimming_single_ended/{sample.sample_id}.json",
            sample=lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}' & downstream_file != downstream_file",
                within=samples,
            ),
        ),
        quant=expand(
            "tmp/fair_rnaseq_salmon_quant/salmon_quant_reads/{sample.species}.{sample.build}.{sample.release}/{sample.sample_id}/quant.sf",
            sample=lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}'",
                within=samples,
            ),
        ),
        quant_genes=expand(
            "tmp/fair_rnaseq_salmon_quant/salmon_quant_reads/{sample.species}.{sample.build}.{sample.release}/{sample.sample_id}/quant.genes.sf",
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
    output:
        report(
            "results/{species}.{build}.{release}/QC/MultiQC_Quantification.html",
            caption="../report/multiqc.rst",
            category="Quality Controls",
            subcategory="General",
            labels={
                "report": "html",
                "organism": "{species}.{build}.{release}",
                "step": "Quantification",
            },
        ),
        "results/{species}.{build}.{release}/QC/MultiQC_Quantification_data.zip",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (2 * 1024) * attempt,
        runtime=lambda wildcards, attempt: int(60 * 0.5) * attempt,
        tmpdir="tmp",
    params:
        extra=lookup(dpath="params/multiqc", within=config),
        use_input_files_only=True,
    log:
        "logs/fair_rnaseq_salmon_quant/multiqc/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_rnaseq_salmon_quant/multiqc/{species}.{build}.{release}.tsv"
    wrapper:
        "v3.3.6/bio/multiqc"
