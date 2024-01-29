import csv
import pandas
import snakemake
import snakemake.utils

from collections import defaultdict
from typing import Any

snakemake.utils.min_version("7.29.0")

# containerized: "docker://snakemake/snakemake:v7.32.4"
# containerized: "docker://mambaorg/micromamba:git-8440cec-jammy-cuda-12.2.0"
# containerized: "docker://condaforge/mambaforge:23.3.1-1"


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


wildcard_constraints:
    sample=r"|".join(samples.sample_id),
    release=r"|".join(release_list),
    build=r"|".join(build_list),
    species=r"|".join(species_list),
    datatype=r"|".join(datatype_list),
    stream=r"|".join(stream_list),


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
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame
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
    salmon_index: list[str] | None = genome_data.get(
        "salmon_index"
    )
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
        "gtf": ancient(genome_data.get("gtf", f"reference/annotation/{species}.{build}.{release}.gtf")),
    }

    if downstream_file or not pandas.isna(downstream_file):
        results["r1"] = f"tmp/fastp/trimmed/{wildcards.sample}.1.fastq"
        results["r2"] = f"tmp/fastp/trimmed/{wildcards.sample}.2.fastq"
    else:
        results["r"] = f"tmp/fastp/trimmed/{wildcards.sample}.fastq"

    return results



def get_salmon_decoy_sequences_input(
    wildcards: snakemake.io.Wildcards,
    genomes: pandas.DataFrame = genomes,
) -> dict[str, list[str]]:
    """
    Return expected input files for Salmon indexation, according to user-input,
    and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe genome files

    Return (dict[str, list[str]]):
    Dictionnary of all input files as required by snakemake-wrapper
    """
    genome_data: dict[str, str | None] = get_reference_genome_data(wildcards, genomes)
    species: str = str(wildcards.species)
    build: str = str(wildcards.build)
    release: str = str(wildcards.release)
    return {
        "transcriptome": genome_data.get("transcripts_fasta", f"reference/sequences/{species}.{build}.{release}.transcripts.fasta"),
        "genome": genome_data.get("dna_fasta", f"reference/sequences/{species}.{build}.{release}.dna.fasta"),
    }


def get_fastp_trimming_input(
    wildcards: snakemake.io.Wildcards,
    samples: pandas.DataFrame = samples,
    config: dict[str, Any] = config,
) -> dict[str, list[str]]:
    """
    Return expected input files for Fastp mapping, according to user-input,
    and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome
    config    (dict[str, Any])        : Configuration file

    Return (dict[str, list[str]]):
    Dictionnary of all input files as required by Fastp's snakemake-wrapper
    """
    sample_data: dict[str, str | None] = get_sample_information(wildcards, samples)
    downstream_file: str | None = sample_data.get("downstream_file")
    if downstream_file and not pandas.isna(downstream_file):
        return {
            "sample": [sample_data["upstream_file"], downstream_file],
        }

    return {
        "sample": [sample_data["upstream_file"]],
    }


def get_aggregate_salmon_counts_input(
    wildcards: snakemake.io.Wildcards,
    samples: pandas.DataFrame = samples,
    config: dict[str, Any] = config,
) -> dict[str, list[str]]:
    """
    Return expected input files for Salmon quant aggregation script,
    according to user-input, and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome
    config    (dict[str, Any])        : Configuration file

    Return (dict[str, list[str]]):
    Dictionnary of all input files as required by Salmon quant aggregation script
    """
    species: str = str(wildcards.species)
    build: str = str(wildcards.build)
    release: str = str(wildcards.release)
    samples_list: list[str] = list(
        samples.loc[
            (samples.species == species)
            & (samples.build == build)
            & (samples.release == release)
        ].sample_id
    )

    aggregate_salmon_counts_input: dict[str, str |list[str]] = {
        "quant": expand(
            "tmp/salmon/quant/{species}.{build}.{release}/{sample}/{quant_file}", # quant.genes.sf",
            quant_file=(["quant.genes.sf"] if str(wildcards.release).lower().startswith("gene") else ["quant.sf"]),
            sample=samples_list,
            species=[species],
            build=[build],
            release=[release],
        ),
    }

    if str(wildcards.counts).lower().startswith("gene"):
        aggregate_salmon_counts_input["tx2gene"] = ancient(f"reference/annotation/{species}.{build}.{release}.id_to_gene.tsv")
    else:
        aggregate_salmon_counts_input["tx2gene"] = ancient(f"reference/annotation/{species}.{build}.{release}.t2g.tsv")

    return aggregate_salmon_counts_input


def get_tximport_input(
    wildcards: snakemake.io.Wildcards,
    samples: pandas.DataFrame = samples,
    config: dict[str, Any] = config,
) -> dict[str, list[str]]:
    """
    Return expected input files for tximport wrapper,
    according to user-input, and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome
    config    (dict[str, Any])        : Configuration file

    Return (dict[str, list[str]]):
    Dictionnary of all input files as required by tximport wrapper
    """
    species: str = str(wildcards.species)
    build: str = str(wildcards.build)
    release: str = str(wildcards.release)
    samples_list: list[str] = list(
        samples.loc[
            (samples.species == species)
            & (samples.build == build)
            & (samples.release == release)
        ].sample_id
    )

    return {
        "quant": expand(
            "tmp/salmon/quant/{species}.{build}.{release}/{sample}/quant.sf",
            sample=samples_list,
            species=[species],
            build=[build],
            release=[release],
        ),
        "quant_genes": expand(
            "tmp/salmon/quant/{species}.{build}.{release}/{sample}/quant.genes.sf",
            sample=samples_list,
            species=[species],
            build=[build],
            release=[release],
        ),
        "lib": expand(
            "tmp/salmon/quant/{species}.{build}.{release}/{sample}/lib_format_counts.json",
            sample=samples_list,
            species=[species],
            build=[build],
            release=[release],
        ),
        "aux_info": expand(
            "tmp/salmon/quant/{species}.{build}.{release}/{sample}/aux_info",
            sample=samples_list,
            species=[species],
            build=[build],
            release=[release],
        ),
        "cmd_info": expand(
            "tmp/salmon/quant/{species}.{build}.{release}/{sample}/cmd_info.json",
            sample=samples_list,
            species=[species],
            build=[build],
            release=[release],
        ),
        "libparams": expand(
            "tmp/salmon/quant/{species}.{build}.{release}/{sample}/libParams",
            sample=samples_list,
            species=[species],
            build=[build],
            release=[release],
        ),
        "logs": expand(
            "tmp/salmon/quant/{species}.{build}.{release}/{sample}/logs",
            sample=samples_list,
            species=[species],
            build=[build],
            release=[release],
        ),
        "tx_to_gene": ancient(f"reference/annotation/{species}.{build}.{release}.id_to_gene.tsv"),
    }


def get_rnaseq_salmon_quant_multiqc_report_input(
    wildcards: snakemake.io.Wildcards, samples: pandas.DataFrame = samples
) -> dict[str, list[str]]:
    """
    Return expected input files for MultiQC report, according to user-input,
    and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome

    Return (dict[str, list[str]]):
    Dictionnary of all input files as required by MultiQC's snakemake-wrapper
    """
    results: dict[str, list[str]] = {"fastqc": [], "salmon": []}
    datatype: str = "dna"
    sample_iterator = zip(
        samples.sample_id,
        samples.species,
        samples.build,
        samples.release,
    )
    for sample, species, build, release in sample_iterator:
        sample_data: dict[str, str | None] = get_sample_information(
            snakemake.io.Wildcards(fromdict={"sample": sample}), samples
        )
        if sample_data.get("downstream_file"):
            results["fastqc"].append(f"results/QC/report_pe/{sample}.1_fastqc.zip")
            results["fastqc"].append(f"results/QC/report_pe/{sample}.2_fastqc.zip")
        else:
            results["fastqc"].append(f"results/QC/report_pe/{sample}_fastqc.zip")

        results["salmon"].append(
            f"tmp/salmon/quant/{species}.{build}.{release}/{sample}/quant.sf"
        )
        results["salmon"].append(
            f"tmp/salmon/quant/{species}.{build}.{release}/{sample}/lib_format_counts.json"
        )
        results["salmon"].append(
            f"tmp/salmon/quant/{species}.{build}.{release}/{sample}/aux_info"
        )
        results["salmon"].append(
            f"tmp/salmon/quant/{species}.{build}.{release}/{sample}/cmd_info.json"
        )
        results["salmon"].append(
            f"tmp/salmon/quant/{species}.{build}.{release}/{sample}/libParams"
        )
        results["salmon"].append(
            f"tmp/salmon/quant/{species}.{build}.{release}/{sample}/logs"
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
            "results/QC/MultiQC_Quantification.html",
            "results/QC/MultiQC_Quantification_data.zip",
        ],
        "quant": [],
    }

    sample_iterator = zip(
        samples.sample_id,
        samples.species,
        samples.build,
        samples.release,
    )

    for sample, species, build, release in sample_iterator:
        for counts in ["Raw", "TPM"]:
            for targets in ["transcripts", "genes"]:
                results["quant"].append(
                    f"results/{species}.{build}.{release}/Quantification/{counts}.{targets}.tsv"
                )

    return results
