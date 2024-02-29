Matierial and methods
=====================

Genome information was download from Ensembl. Samtools_ [#samtoolspaper]_ 
and Picard_ [#gatkpaper]_ were used to index genome sequences.
Agat_ [#agatpaper]_ was used to correct common issues found in Ensembl
genome annotation files. Salmon_ [#salmonpaper]_ was used to preduce and index a 
`decoy-aware gentrome`_ based on both DNA and cDNA sequences.

Raw fastq_ files were trimmed using Fastp_ [#fastppaper]_. Salmon_ [#salmonpaper]
performed the pseudo-mapping and estimation of transcripts abundance. The count
aggregation was performed with tximport_ [#tximportpaper].

The  whole pipeline was powered by Snakemake_ [#snakemakepaper]_.

This pipeline is freely available on Github_, details about installation
usage, and resutls can be found on the `Snakemake workflow`_ page.

.. [#samtoolspaper] Li, Heng, et al. "The sequence alignment/map format and SAMtools." bioinformatics 25.16 (2009): 2078-2079.
.. [#gatkpaper] McKenna, Aaron, et al. "The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data." Genome research 20.9 (2010): 1297-1303.
.. [#agatpaper] Dainat J. AGAT: Another Gff Analysis Toolkit to handle annotations in any GTF/GFF format.  (Version v0.7.0). Zenodo. https://www.doi.org/10.5281/zenodo.3552717
.. [#fastppaper] Chen, Shifu, et al. "fastp: an ultra-fast all-in-one FASTQ preprocessor." Bioinformatics 34.17 (2018): i884-i890.
.. [#salmonpaper] Patro, Rob, et al. "Salmon provides fast and bias-aware quantification of transcript expression." Nature methods 14.4 (2017): 417-419.
.. [#tximportpaper] Love, Michael I., Charlotte Soneson, and Mark D. Robinson. "Importing transcript abundance datasets with tximport." Dim Txi. Inf. Rep. Sample 1.1 (2017): 5.
.. [#snakemakepaper] Köster, Johannes, and Sven Rahmann. "Snakemake—a scalable bioinformatics workflow engine." Bioinformatics 28.19 (2012): 2520-2522.

.. _Snakemake: https://snakemake.readthedocs.io
.. _Github: https://github.com/tdayris/fair_rnaseq_salmon_quant
.. _`Snakemake workflow`: https://snakemake.github.io/snakemake-workflow-catalog?usage=tdayris/fair_rnaseq_salmon_quant
.. _Picard: https://snakemake-wrappers.readthedocs.io/en/v3.3.6/wrappers/picard/createsequencedictionary.html
.. _Samtools: https://snakemake-wrappers.readthedocs.io/en/v3.3.6/wrappers/samtools/faidx.html
.. _Agat: https://agat.readthedocs.io/en/latest/index.html
.. _Salmon: https://snakemake-wrappers.readthedocs.io/en/v3.3.6/wrappers/salmon.html
.. _`decoy-aware gentrome`: https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode
.. _Fastp: https://snakemake-wrappers.readthedocs.io/en/v3.3.6/wrappers/fastp.html
.. _fastq: https://fr.wikipedia.org/wiki/FASTQ
.. _tximport: https://snakemake-wrappers.readthedocs.io/en/v3.3.6/wrappers/tximport.html