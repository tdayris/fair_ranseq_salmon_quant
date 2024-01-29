Snakemake workflow used to estimate transcripts/genes abundance with Salmon.

## Usage

The usage of this workflow is described in the [Snakemake workflow catalog](https://snakemake.github.io/snakemake-workflow-catalog?usage=tdayris/fair_rnaseq_salmon_quant) 
it is also available [locally](https://github.com/tdayris/fair_rnaseq_salmon_quant/blob/main/workflow/report/usage.rst) on a single page.

## Results

A complete description of the results can be found here in [workflow reports](https://github.com/tdayris/fair_rnaseq_salmon_quant/blob/main/workflow/report/results.rst).

## Material and Methods

The tools used in this pipeline are described [here](https://github.com/tdayris/fair_rnaseq_salmon_quant/blob/main/workflow/report/material_methods.rst) textually.

### Index and genome sequences with [`fair_genome_indexer`](https://github.com/tdayris/fair_genome_indexer/tree/main)

#### Get DNA sequences

| Step                             | Commands                                                                                                         |
| -------------------------------- | ---------------------------------------------------------------------------------------------------------------- |
| Download DNA Fasta from Ensembl  | [ensembl-sequence](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/reference/ensembl-sequence.html) |
| Remove non-canonical chromosomes | [pyfaidx](https://github.com/mdshw5/pyfaidx)                                                                     |
| Index DNA sequence               | [samtools](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/samtools/faidx.html)                     |
| Creatse sequence Dictionary      | [picard](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/picard/createsequencedictionary.html)      |

#### Get genome annotation (GTF)

| Step                                                       | Commands                                                                                                             |
| ---------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------- |
| Download GTF annotation                                    | [ensembl-annotation](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/reference/ensembl-annotation.html) |
| Fix format errors                                          | [Agat](https://agat.readthedocs.io/en/latest/tools/agat_convert_sp_gff2gtf.html)                                     |
| Remove non-canonical chromosomes, based on above DNA Fasta | [Agat](https://agat.readthedocs.io/en/latest/tools/agat_sq_filter_feature_from_fasta.html)                           |
| Remove `<NA>` Transcript support levels                    | [Agat](https://agat.readthedocs.io/en/latest/tools/agat_sp_filter_feature_by_attribute_value.html)                   |


#### Quality controls

| Step     | Wrapper                                                                                                                          |
| -------- | -------------------------------------------------------------------------------------------------------------------------------- |
| FastQC   | [fastqc-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/fastqc.html)                                       |
| MultiQC  | [multiqc-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/multiqc.html)                                     |

### Read abundance estimation with [salmon-tximport meta-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/meta-wrappers/salmon_tximport.html)

#### Indexation

| Step                                | Wraper                                                                                                              |
| ----------------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| Create gentrome and decoy sequences | [generate-decoy](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/salmon/decoys.html#bio-salmon-decoys) |
| Index decoy aware gentrome          | [salmon-index](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/salmon/index.html#bio-salmon-index)     |


#### Abundance estimation

| Step                    | Wrapper                                                                                                                 |
| ----------------------- | ----------------------------------------------------------------------------------------------------------------------- |
| Trimm raw reads         | [fastp](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/fastp.html)                                        |
| Estimate abundances     | [salmon-quand](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/salmon/quant.html#bio-salmon-quant)         |
| Aggregate counts in TSV | in-house [script](https://github.com/tdayris/fair_rnaseq_salmon_quant/blob/main/workflow/scripts/regenerate_genomes.py) |
| Aggregate counts in R   | [tximport](https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/tximport.html#bio-tximport)                     |