import csv
import os
import pandas
import os
import snakemake
import snakemake.utils

from collections import defaultdict
from typing import Any

snakemake.utils.min_version("8.5.3")


container: "docker://snakemake/snakemake:v8.16.0"


# Load and check configuration file
configfile: "config/config.yaml"


snakemake.utils.validate(config, "../schemas/config.schema.yaml")

# Load and check samples properties table
sample_table_path: str = config.get("samples", "config/samples.csv")
with open(sample_table_path, "r") as sample_table_stream:
    dialect: csv.Dialect = csv.Sniffer().sniff(sample_table_stream.readline())
    sample_table_stream.seek(0)

samples: pandas.DataFrame = pandas.read_csv(
    filepath_or_buffer=sample_table_path,
    sep=dialect.delimiter,
    header=0,
    index_col=None,
    comment="#",
    dtype=str,
)
samples = samples.where(samples.notnull(), None)
snakemake.utils.validate(samples, "../schemas/samples.schema.yaml")

# This is here for compatibility with
genome_table_path: str = config.get("genomes")
if genome_table_path:
    with open(genome_table_path, "r") as genome_table_stream:
        dialect: csv.Dialect = csv.Sniffer().sniff(genome_table_stream.readline())
        genome_table_stream.seek(0)

    genomes: pandas.DataFrame = pandas.read_csv(
        filepath_or_buffer=genome_table_path,
        sep=dialect.delimiter,
        header=0,
        index_col=None,
        comment="#",
        dtype=str,
    )
    genomes = genomes.where(genomes.notnull(), None)
else:
    genomes: pandas.DataFrame = samples[
        ["species", "build", "release"]
    ].drop_duplicates(keep="first", ignore_index=True)
    genomes.to_csv("genomes.csv", sep=",", index=False, header=True)
    config["genomes"] = "genomes.csv"

snakemake.utils.validate(genomes, "../schemas/genomes.schema.yaml")


report: "../report/workflow.rst"


release_list: list[str] = list(set(genomes.release.tolist()))
build_list: list[str] = list(set(genomes.build.tolist()))
species_list: list[str] = list(set(genomes.species.tolist()))
datatype_list: list[str] = ["dna", "cdna", "gentrome"]
stream_list: list[str] = ["1", "2"]
tmp: str = f"{os.getcwd()}/tmp"
snakemake_wrappers_prefix: str = "v3.13.7"


wildcard_constraints:
    sample=r"|".join(samples.sample_id),
    release=r"|".join(release_list),
    build=r"|".join(build_list),
    species=r"|".join(species_list),
    datatype=r"|".join(datatype_list),
    stream=r"|".join(stream_list),


def lookup_config(
    dpath: str, default: str | None = None, config: dict[str, Any] = config
) -> str:
    """
    Run lookup function with default parameters in order to search a key in configuration and return a default value
    """
    value: str | None = default

    try:
        value = lookup(dpath=dpath, within=config)
    except LookupError:
        value = default
    except WorkflowError:
        value = default

    return value


def lookup_genomes(
    wildcards: snakemake.io.Wildcards,
    key: str,
    default: str | list[str] | None = None,
    genomes: pandas.DataFrame = genomes,
) -> str:
    """
    Run lookup function with default parameters in order to search user-provided sequence/annotation files
    """
    query = str(
        "species == '{wildcards.species}' & build == '{wildcards.build}' & release == '{wildcards.release}'".format(
            wildcards=wildcards
        )
    )
    return getattr(lookup(query=query, within=genomes), key, default)


def get_gtf(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final genome annotation
    """
    default = str(
        "reference/annotation/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.gtf".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(wildcards, key="gtf", default=default, genomes=genomes)


def get_tx2gene(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to final tx2gene table
    """
    default = str(
        "reference/annotation/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.t2g.tsv".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(wildcards, key="tx_to_gene", default=default)


def get_id2gene(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to final id2gene table
    """
    default = str(
        "reference/annotation/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.id_to_gene.tsv".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(wildcards, key="id_to_gene", default=default)


def get_sample_information(
    wildcards: snakemake.io.Wildcards, samples: pandas.DataFrame
) -> dict[str, str | None]:
    """
    Return sample information for a given {sample} wildcards

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe samples and their input files

    Return (dict[str, str | None]):
    Sample information
    """
    result: str | None = samples.loc[(samples["sample_id"] == str(wildcards.sample))]
    if len(result) > 0:
        return next(iter(result.to_dict(orient="index").values()))
    return defaultdict(lambda: None)


def get_reference_genome_data(
    wildcards: snakemake.io.Wildcards,
    genomes: pandas.DataFrame = genomes,
) -> dict[str, str | None]:
    """
    Return genome information for a given set of {species, build, release} wildcards

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    genomes   (pandas.DataFrame)      : Describe genomes and reference file(s)

    Return (dict[str, str | None]):
    Genome information
    """
    result: str | None = genomes.loc[
        (genomes["species"] == str(wildcards.species))
        & (genomes["build"] == str(wildcards.build))
        & (genomes["release"] == str(wildcards.release))
    ]
    if len(result) > 0:
        return next(iter(result.to_dict(orient="index").values()))
    return defaultdict(lambda: None)


def get_salmon_quant_reads_input(
    wildcards: snakemake.io.Wildcards,
    samples: pandas.DataFrame = samples,
    genomes: pandas.DataFrame = genomes,
    config: dict[str, Any] = config,
) -> dict[str, list[str]]:
    """
    Return expected input files for Salmon mapping, according to user-input,
    and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome
    genomes   (pandas.DataFrame)      : Describe genome and index files paths
    config    (dict[str, Any])        : Configuration file

    Return (dict[str, list[str]]):
    Dictionnary of all input files as required by Salmon's snakemake-wrapper
    """
    species: str = str(wildcards.species)
    release: str = str(wildcards.release)
    build: str = str(wildcards.build)
    sample_data: dict[str, str | None] = get_sample_information(wildcards, samples)
    downstream_file: str | None = sample_data.get("downstream_file")
    genome_data: dict[str, str | None] = get_reference_genome_data(wildcards, genomes)
    salmon_index: list[str] | None = genome_data.get("salmon_index")
    if salmon_index:
        salmon_index = list(map(str, Path(salmon_index).iterdir()))
    else:
        salmon_index = multiext(
            f"reference/salmon_index/{species}.{build}.{release}/{species}.{build}.{release}/",
            "complete_ref_lens.bin",
            "ctable.bin",
            "ctg_offsets.bin",
            "duplicate_clusters.tsv",
            "info.json",
            "mphf.bin",
            "pos.bin",
            "pre_indexing.log",
            "rank.bin",
            "refAccumLengths.bin",
            "ref_indexing.log",
            "reflengths.bin",
            "refseq.bin",
            "seq.bin",
            "versionInfo.json",
        )

    results: dict[str, str | list[str]] = {
        "index": ancient(salmon_index),
        "gtf": get_gtf(wildcards),
    }

    if downstream_file or not pandas.isna(downstream_file):
        results["r1"] = str(
            f"tmp/fair_rnaseq_salmon_quant_fastp_trimming_pair_ended/{wildcards.sample}.1.fastq.gz"
        )
        results["r2"] = str(
            f"tmp/fair_rnaseq_salmon_quant_fastp_trimming_pair_ended/{wildcards.sample}.2.fastq.gz"
        )
    else:
        results["r"] = str(
            f"tmp/fair_rnaseq_salmon_quant_fastp_trimming_single_ended/{wildcards.sample}.fastq.gz"
        )

    return results


def get_rnaseq_salmon_quant_target(
    wildcards: snakemake.io.Wildcards,
    samples: pandas.DataFrame = samples,
    config: dict[str, Any] = config,
) -> dict[str, list[str]]:
    """
    Return the expected list of output files at the end of the pipeline

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome
    config    (dict[str, Any])        : Configuration file

    Return (dict[str, List(str)]):
    Dictionnary of expected output files
    """
    results: dict[str, list[str]] = {
        "multiqc": [
            "results/QC/MultiQC_FastQC.html",
            "results/QC/MultiQC_FastQC_data.zip",
        ],
        "quant": [],
        "datavzrd": [],
        "stats_csv": [],
        "quant_pca": [],
    }

    sample_iterator = zip(
        samples.species,
        samples.build,
        samples.release,
    )

    genome_data = list(
        set(
            f"{species}.{build}.{release}"
            for species, build, release in sample_iterator
        )
    )

    for genome_version in genome_data:
        for counts in ["Raw", "TPM"]:
            for targets in ["transcripts", "genes"]:
                results["quant"].append(
                    f"results/{genome_version}/Quantification/{counts}.{targets}.tsv"
                )

        results["multiqc"].append(
            f"results/{genome_version}/QC/MultiQC_Quantification.html"
        )
        results["datavzrd"].append(
            f"results/{genome_version}/Quantification/html_reports/TPM.genes"
        )
        results["stats_csv"].append(f"results/{genome_version}/QC/Stats.csv.gz")
        results["quant_pca"].append(f"results/{genome_version}/Quantification/PCA.png")
        results["quant_pca"].append(
            f"results/{genome_version}/Quantification/Scree.png"
        )

    results["multiqc"] = list(set(results["multiqc"]))
    return results
