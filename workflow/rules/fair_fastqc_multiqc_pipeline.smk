module fair_fastqc_multiqc:
    snakefile:
        github("tdayris/fair_fastqc_multiqc", path="workflow/Snakefile", tag="1.0.1")
    config:
        {
            "samples": config.get("samples", "config/samples.csv"),
            "params": {
                "fastqc": config.get("params", {}).get("fastqc"),
                "multiqc": config.get("params", {}).get("multiqc"),
            },
            "genomes": config.get("genomes", "config/genomes.csv"),
        }


use rule * from fair_fastqc_multiqc as fair_fastqc_multiqc_*


use rule multiqc_report from fair_fastqc_multiqc as rnaseq_salmon_quant_multiqc_report with:
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
            },
        ),
        "results/QC/MultiQC_Quantification_data.zip",
