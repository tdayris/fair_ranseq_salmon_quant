SHELL=/usr/bin/bash
.ONESHELL:
.SHELLFLAGS := -euic
.DELETE_ON_ERROR: linter_info.txt format.txt pipeline.txt

MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

# List of snakefiles to lint
snakefiles := \
	"../workflow/rules/common.smk" \
	"../workflow/rules/datavzrd.smk" \
	"../workflow/rules/fair_genome_indexer_pipeline.smk" \
	"../workflow/rules/fair_fastqc_multiqc_pipeline.smk" \
	"../workflow/rules/fastp.smk" \
	"../workflow/rules/in_house_pca.smk" \
	"../workflow/rules/in_house_qc.smk" \
	"../workflow/rules/merge_salmon.smk" \
	"../workflow/rules/multiqc.smk" \
	"../workflow/rules/salmon_quant_meta.smk" \
	"../workflow/Snakefile"


python_scripts := \
	"../workflow/scripts/bioinfokit_pca.py" \
	"../workflow/scripts/build_datavzrd_yaml.py" \
	"../workflow/scripts/fair_rnaseq_salmon_quant_multiqc_config.py" \
	"../workflow/scripts/in_house_qc_mt_ribo_percents.py" \
	"../workflow/scripts/merge_salmon_quant.py" \
	"../workflow/scripts/qc_table.py" \
	"../workflow/scripts/regenerate_genomes.py"

.PHONY: all
all: report.txt


dag.png: linter.txt
	snakemake -s ../workflow/Snakefile \
		--rulegraph --forceall | dot -Tpng > dag.png && \
	cp --verbose dag.png ../dag.png


report.txt: pipeline.txt
	snakemake -s ../workflow/Snakefile \
		--report report.zip 2>&1 | tee report.txt


pipeline.txt: linter.txt format.txt black.txt
	snakemake -s ../workflow/Snakefile \
		--cores 2 \
		--software-deployment-method apptainer conda\
		2>&1 | tee pipeline.txt


format.txt:
	snakefmt $(snakefiles) 2>&1 | tee format.txt


linter.txt: format.txt black.txt
	snakemake -s ../workflow/Snakefile \
		--lint 2>&1 | tee linter.txt


black.txt:
	black $(python_scripts) 2>&1 | tee black.txt

.PHONY: clean
clean:
	rm -rvf linter.txt \
		   format.txt \
		   reference \
		   results \
		   tmp \
		   logs \
		   benchmark \
		   report.html \
		   report.txt \
		   pipeline.txt \
		   genomes.csv \
		   dag.png \
		   black.txt \
		   report.zip \
		   report
