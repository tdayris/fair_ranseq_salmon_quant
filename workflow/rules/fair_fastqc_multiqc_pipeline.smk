module fair_fastqc_multiqc:
    snakefile:
        github("tdayris/fair_fastqc_multiqc", path="workflow/Snakefile", tag="2.2.4")
    config:
        config


use rule * from fair_fastqc_multiqc
