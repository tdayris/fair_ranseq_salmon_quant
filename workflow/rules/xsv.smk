rule xsv_sort:
    input:
        table="tmp/salmon/quant/{species}.{build}.{release}/{sample}/quant{content}.sf"
    output:
        temp("tmp/xsv/sort/{species}.{build}.{release}/{sample}.csv")
    log:
        "logs/xsv/sort/{species}.{build}.{release}/{sample}.log"
    benchmark:
        "benchmark/xsv/sort/{species}.{build}.{release}/{sample}.tsv"
    params:
        subcommand="sort",
        extra="--delimiter $'\t'",
    wrapper:
        f"{snakemake_wrappers_version}/utils/xsv"
