
Reference
=========

Alongside with this report, you may find a directory called `reference`.
You shall find all requested files in it. By default, the following
files are present:

::

    reference/
    ├── XXX.all.vcf
    ├── XXX.cdna.fasta
    ├── XXX.cdna.fasta.fai
    ├── XXX.dna.dict
    ├── XXX.dna.fasta
    ├── XXX.dna.fasta.fai
    ├── XXX.tsv
    └── XXX.gtf


+---------------+-----------------------------+
| Extension     | Content                     |
+===============+=============================+
| `.gtf`        | Genome annotation           |
+---------------+-----------------------------+
| `.tsv`        | Genome id-to-name           |
+---------------+-----------------------------+
| `.fasta`      | Genome sequences            |
+---------------+-----------------------------+
| `.fasta.fai`  | Genome sequences index      |
+---------------+-----------------------------+
| `.dict`       | Genome sequences dictionary |
+---------------+-----------------------------+
| `.vcf`        | Genome known variations     |
+---------------+-----------------------------+

These files are quite volumous and are not embeded in this HTML page. Please
find them directly on file system.

Results
=======

```
results/
├── QC
│   ├── MultiQC_FastQC_data.zip                 # Raw fastq file quality tables
│   ├── MultiQC_FastQC.html                     # Raw fastq file quality report
│   └── report_pe
│       ├── YYY_fastqc.zip                      # Per-sample quality table
│       └── YYY.html                            # Per-sample quality report
└── XXX
    ├── QC
    │   ├── Mitochodrial_ratio.png              # Mitochodrial ratio across samples
    │   ├── MultiQC_Quantification_data.zip     # Complete pipeline statistics as TSV
    │   ├── MultiQC_Quantification.html         # Quantification quality report
    │   ├── rRNA_ratio.png                      # rRNA ratio across samples
    │   └── Stats.csv.gz                        # (pseudo-)mapping statitstics
    └── Quantification
        ├── html_reports/                       # View counts in your web-browser
        ├── Raw.genes.tsv                       # Raw gene abundance estimation
        ├── Raw.transcripts.tsv                 # Raw transcripts abundance estimation
        ├── TPM.genes.tsv                       # Normalized gene counts
        └── TPM.transcripts.tsv                 # Normalized transcripts counts
```


