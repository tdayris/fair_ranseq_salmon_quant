module fair_genome_indexer:
    snakefile:
        github("tdayris/fair_genome_indexer", path="workflow/Snakefile", tag="2.1.0")
    config:
        config


use rule * from fair_genome_indexer
