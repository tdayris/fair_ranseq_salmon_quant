rule fair_rnaseq_salmon_quant_mito_rrna_percents:
    input:
        table="results/{species}.{build}.{release}/Quantification/TPM.genes.tsv",
    output:
        mito_csv="tmp/fair_rnaseq_salmon_quant/mito_rrna_percents/{species}.{build}.{release}/mito.csv",
        rrna_csv="tmp/fair_rnaseq_salmon_quant/mito_rrna_percents/{species}.{build}.{release}/rrna.csv",
        mito_png="results/{species}.{build}.{release}/QC/Mitochodrial_ratio.png",
        rrna_png="results/{species}.{build}.{release}/QC/rRNA_ratio.png",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 2,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir=tmp,
    log:
        "logs/fair_rnaseq_salmon_quant/mito_rrna_percents/{species}.{build}.{release}/TPM.genes.log",
    benchmark:
        "benchmark/fair_rnaseq_salmon_quant/mito_rrna_percents/{species}.{build}.{release}/TPM.genes.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/fair_rnaseq_salmon_quant_in_house_qc_mt_ribo_percents.py"


rule fair_rnaseq_salmon_quant_qc_table:
    input:
        mito="tmp/fair_rnaseq_salmon_quant/mito_rrna_percents/{species}.{build}.{release}/mito.csv",
        rrna="tmp/fair_rnaseq_salmon_quant/mito_rrna_percents/{species}.{build}.{release}/rrna.csv",
        general="tmp/fair_rnaseq_salmon_quant/unzip_multiqc_data/{species}.{build}.{release}/multiqc_general_stats.txt",
        salmon=collect(
            "tmp/fair_rnaseq_salmon_quant_salmon_quant_reads/{sample.species}.{sample.build}.{sample.release}/{sample.sample_id}/aux_info",
            sample=lookup(
                query="species == '{species}' & build == '{build}' & release == '{release}'",
                within=samples,
            ),
        ),
    output:
        stats="results/{species}.{build}.{release}/QC/Stats.csv.gz",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir=tmp,
    log:
        "logs/fair_rnaseq_salmon_quant/qc_table/{species}.{build}.{release}.log",
    benchmark:
        "benchmark/fair_rnaseq_salmon_quant/qc_table/{species}.{build}.{release}.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/fair_rnaseq_salmon_quant_qc_table.py"
