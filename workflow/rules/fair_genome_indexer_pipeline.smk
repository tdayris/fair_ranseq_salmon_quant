module fair_genome_indexer:
    snakefile:
        github("tdayris/fair_genome_indexer", path="workflow/Snakefile", tag="2.3.2")
    config:
        {"genomes": config.get("genome.csv", "config/genome.csv")}


use rule * from fair_genome_indexer as fair_genome_indexer_*
