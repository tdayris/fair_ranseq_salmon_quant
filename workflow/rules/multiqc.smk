rule fair_rnaseq_salmon_quant_multiqc_config:
    input:
        "tmp/fair_fastqc_multiqc_bigr_logo.png",
    output:
        temp(
            "tmp/fair_rnaseq_salmon_quant/{species}.{build}.{release}/multiqc_config.yaml"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512,
        runtime=lambda wildcards, attempt: attempt * 5,
        tmpdir=tmp,
    log:
        "logs/fair_rnaseq_salmon_quant/multiqc_config/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_rnaseq_salmon_quant/multiqc_config/{species}.{build}.{release}.tsv"
    params:
        extra=lambda wildcards, input: {
            "title": "RNA-Seq gene abundance estimation report",
            "subtitle": "From raw fastq to decoy aware gene abundance estimation with Salmon",
            "intro_text": (
                f"This pipeline building this report has been "
                f"fine tuned for {wildcards.species}.{wildcards.build}"
                f".{wildcards.release}. The sequenced library is expected "
                f"to be bulk RNA-Seq. Wet-lab experimental design was not "
                f"considered."
            ),
            "report_comment": (
                "This report was generated using: "
                "https://github.com/tdayris/fair_rnaseq_salmon_quant"
            ),
            "show_analysis_paths": False,
            "show_analysis_time": False,
            "custom_logo": input[0],
            "custom_logo_url": "https://www.gustaveroussy.fr/en",
            "custom_logo_title": "Bioinformatics Platform @ Gustave Roussy",
            "report_header_info": [
                {"Contact E-mail": "bigr@gustaveroussy.fr"},
                {"Applivation type": "Bulk RNA-Seq"},
                {"Project Type": "Quantification (Abundance estimation)"},
            ],
            "software_versions": {
                "Quality controls": {
                    "fastqc": "1.12.1",
                    "fastq_screen": "0.15.3",
                    "bowtie1": "1.3.1",
                    "bowtie2": "2.5.2",
                    "multiqc": "1.21.0",
                    "deeptools": "3.5.5",
                    "python": "3.12.0",
                    "seaborn": "0.13.2",
                    "pyyaml": "6.0.1",
                    "pandas": "2.1.3",
                    "agat": "1.2.0",
                },
                "Trimming": {
                    "fastp": "0.23.4",
                },
                "Quantification": {
                    "salmon": "1.10.3",
                    "gzip": "1.13",
                    "pzip2": "1.0.8",
                },
                "Pipeline": {
                    "snakemake": "8.5.3",
                    "fair_genome_indexer": "3.5.0",
                    "fair_fastqc_multiqc": "2.2.8",
                    "fair_rnaseq_salmon_quant": "1.1.0",
                },
            },
            "disable_version_detection": True,
            "run_modules": [
                "fastqc",
                "fastq_screen",
                "salmon",
                "fastp",
                "deeptools",
            ],
            "report_section_order": {
                "fastqc": {"order": 1000},
                "fastq_screen": {"before": "fastqc"},
                "fastp": {"before": "fastq_screen"},
                "salmon": {"before": "fastp"},
                "deeptools": {"before": "salmon"},
                "software_versions": {"before": "deeptools"},
            },
        },
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/fair_rnaseq_salmon_quant_multiqc_config.py"


rule fair_rnaseq_salmon_quant_multiqc_report:
    input:
        "tmp/fair_rnaseq_salmon_quant/{species}.{build}.{release}/multiqc_config.yaml",
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
            "tmp/fair_rnaseq_salmon_quant_fastp_trimming_pair_ended/{sample.sample_id}.json",
            sample=lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}' & downstream_file == downstream_file",
                within=samples,
            ),
        ),
        fastp_se_json=expand(
            "tmp/fair_rnaseq_salmon_quant_fastp_trimming_single_ended/{sample.sample_id}.json",
            sample=lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}' & downstream_file != downstream_file",
                within=samples,
            ),
        ),
        quant=expand(
            "tmp/fair_rnaseq_salmon_quant_salmon_quant_reads/{sample.species}.{sample.build}.{sample.release}/{sample.sample_id}/quant.sf",
            sample=lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}'",
                within=samples,
            ),
        ),
        quant_genes=expand(
            "tmp/fair_rnaseq_salmon_quant_salmon_quant_reads/{sample.species}.{sample.build}.{sample.release}/{sample.sample_id}/quant.genes.sf",
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
        tmpdir=tmp,
    params:
        extra=lookup_config(
            dpath="params/fair_rnaseq_salmon_quant/multiqc",
            default="--verbose --no-megaqc-upload --no-ansi --force",
        ),
        use_input_files_only=True,
    log:
        "logs/fair_rnaseq_salmon_quant/multiqc/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_rnaseq_salmon_quant/multiqc/{species}.{build}.{release}.tsv"
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/multiqc"


rule fair_rnaseq_salmon_quant_unzip_multiqc_data:
    input:
        "results/{species}.{build}.{release}/QC/MultiQC_Quantification_data.zip",
    output:
        multiext(
            "tmp/fair_rnaseq_salmon_quant/unzip_multiqc_data/{species}.{build}.{release}/multiqc_",
            "citations.txt",
            "data.json",
            "fastp.txt",
            "fastqc.txt",
            "general_stats.txt",
            "salmon.txt",
            "software_versions.txt",
            "sources.txt",
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512,
        runtime=lambda wildcards, attempt: attempt * 2,
        tmpdir=tmp,
    log:
        "logs/fair_rnaseq_salmon_quant/unzip_mqc_data/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_rnaseq_salmon_quant/unzip_mqc_data/{species}.{build}.{release}.tsv"
    params:
        outdir=lambda wildcards, output: os.path.dirname(output[0]),
        extra=lookup_config(dpath="params/fair_rnaseq_salmon_quant/unzip", default="-o"),
    conda:
        "../envs/bash.yaml"
    shell:
        "unzip -d {params.outdir} {params.extra} {input} > {log} 2>&1"
