rule quantification_multiqc_report:
    input:
        unpack(get_rnaseq_salmon_quant_multiqc_report_input),
    output:
        report(
            "results/QC/MultiQC_Quantification.html",
            caption="../report/multiqc.rst",
            category="Quality Controls",
            subcategory="General",
            labels={
                "report": "html",
                "step": "Quantification",
            },
        ),
        "results/QC/MultiQC_Quantification_data.zip",
    threads: 1
    resources:
        # Reserve 2Gb per attempt
        mem_mb=lambda wildcards, attempt: (2 * 1024) * attempt,
        # Reserve 30min per attempt
        runtime=lambda wildcards, attempt: int(60 * 0.5) * attempt,
        tmpdir="tmp",
    params:
        extra=config.get("params", {}).get("multiqc", "--module fastqc --module fastp --module salmon --zip-data-dir --verbose --no-megaqc-upload --no-ansi --force"),
        use_input_files_only=True,
    log:
        "logs/multiqc.log",
    benchmark:
        "benchmark/multiqc.tsv"
    wrapper:
        "v3.3.3/bio/multiqc"