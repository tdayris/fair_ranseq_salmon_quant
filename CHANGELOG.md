# 1.2.1

## Features:

* Snakemake wrappers up to 3.10.1
* fair_genome_indexer up to 3.4.4
* fair_fastqc_multiqc up to 2.2.7
* log, benchmark, and temp files renamed to follow rule names

# 1.2.0

## Features:

* rRNA ratio and Mitocholdrial ratio added to QC


# 1.1.1

## Features:

* Resources reservation have been set accordingly to user requests
* Documentation update

# 1.1.0

## Features:

* all configuration keys are now optional
* fair_fastqc_multiqc_pipeline update to 2.2.3
* fair_genome_indexer update to 3.4.0
* snakemake-wrappers update to 3.7.0
* MultiQC configuration file

## Fixes:

* Documentation errors

# 1.0.2

## Features:

* Use `lookup` and `branch` instead of hand made functions
* Pipeline can be used through mamba + apptainer
* fair_fastqc_multiqc_pipeline version 2.0.4
* fair_genome_indexer version 3.1.4
* snakemake wrappers v3.3.6

# 1.0.1

## Features

* Salmon quant merge now annotates transcripts as well as genes
* DataVzrd to explore quantification if needed

## Fix

* Missing fastp report page

# 1.0.0

## Features

* Control/Clean fastq files
* Estimate quantification
* Build report (snakemake + datavzrd)
* Build a single table with all counts from all samples
* Snakemake-wrappers v3.3.3
* Snakemake v8+ compatible
* fair_genome_indexer v3.0.0
