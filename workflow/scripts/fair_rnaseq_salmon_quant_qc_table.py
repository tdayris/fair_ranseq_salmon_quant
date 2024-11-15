# coding: utf-8

import json
import logging
import pandas

from pathlib import Path, PurePath
from typing import Any


def read_salmon_qc(path: str) -> pandas.DataFrame:
    """
    !WARNING! Legacy
    MultiQC update does not return this file anymore.

    Return MultiQC stat table as a data frame

    Parameters:
    path (str): Path to multiqc salmon table

    Results:
    (pandas.DataFrame): Salmon table
    """
    logging.info(f"Loading Salmon's MultiQC table at {path=}")
    mqc: pandas.DataFrame = pandas.read_csv(
        path,
        sep="\t",
        header=0,
        index_col=0,
    )
    mqc = mqc[["salmon_version", "frag_length_mean", "frag_length_sd"]]
    mqc.columns = [
        "Salmon_version",
        "Salmon_Mean_fragment_length",
        "Salmon_Fragment_length_sd",
    ]
    logging.debug(mqc.head())
    return mqc


def read_salmon_auxdir(paths: list[str | Path]) -> pandas.DataFrame:
    """
    Search for meta_info.json files and load them
    """
    frag_length = {}
    for path in paths:
        if isinstance(path, str):
            path = Path(path) / "meta_info.json"

        if not path.exists():
            raise FileNotFoundError(f"Could not find Salmon's meta_json file at {path}")

        with path.open() as json_stream:
            meta_json = json.load(json_stream)

        sample_id = PurePath(path).parts[-3]
        frag_length[sample_id] = {
            "Salmon_version": meta_json["salmon_version"],
            "Salmon_Mean_fragment_length": meta_json["frag_length_mean"],
            "Salmon_Fragment_length_sd": meta_json["frag_length_sd"],
        }
    return pandas.DataFrame.from_dict(frag_length, orient="index")


def read_general_stat_table(
    path: str, samples_to_keep: list[str] | pandas.Series
) -> pandas.DataFrame:
    """
    Return MultiQC general statistics as a data frame

    Parameters:
    path (str): Path to the general stats table

    Return:

    """
    logging.info(f"Loading MultiQC's general table at {path=}")
    mqc: pandas.DataFrame = pandas.read_csv(
        path,
        sep="\t",
        header=0,
        index_col=0,
    )
    mqc = mqc[mqc.index.isin(samples_to_keep)]
    mqc = mqc[
        [
            "Salmon_mqc_generalstats_salmon_percent_mapped",
            "Salmon_mqc_generalstats_salmon_num_mapped",
            "Salmon_mqc_generalstats_salmon_library_types",
            "Salmon_mqc_generalstats_salmon_strand_mapping_bias",
            "fastp_mqc_generalstats_fastp_after_filtering_gc_content",
            "fastp_mqc_generalstats_fastp_pct_adapter",
            "fastp_mqc_generalstats_fastp_filtering_result_passed_filter_reads",
        ]
    ]
    mqc.columns = [
        "Salmon_Percent_mapped",
        "Salmon_Num_mapped",
        "Salmon_Library_Type",
        "Salmon_Strand_mapping_bias",
        "Fastp_GC_content",
        "Fastp_Percent_adapters",
        "Fastp_nb_read_passing_trimming",
    ]
    logging.debug(mqc.head())
    return mqc


def read_mito_ratio_table(path: str) -> pandas.DataFrame:
    """
    Return the Mitochondrial counts ratio as a data frame

    Parameters:
    path (str): Path to the table

    Return:
    (pandas.DataFrame): The mito ratio table
    """
    logging.info(f"Loading mito ratio available at {path=}")
    mito: pandas.DataFrame = pandas.read_csv(
        path,
        sep=",",
        header=0,
        index_col=0,
    )[["Ratio", "mito_counts"]]
    mito.columns = ["Mito_Ratio", "mito_counts"]
    logging.debug(mito.head())
    return mito


def read_rrna_ratio_table(path: str) -> pandas.DataFrame:
    """
    Return the rRNA counts ration as a data frame

    Parameters:
    path (str): Path to the table

    Return:
    (pandas.DataFrame): The rrna ratio table
    """
    logging.info(f"Loading rRNA ratio table at {path=}")
    rrna: pandas.DataFrame = pandas.read_csv(
        path,
        sep=",",
        header=0,
        index_col=0,
    )[["Ratio", "rRNA_counts"]]
    rrna.columns = ["rRNA_Ratio", "rRNA_counts"]
    logging.debug(rrna.head())
    return rrna


if __name__ == "__main__":
    logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)

    logging.getLogger("matplotlib.font_manager").disabled = True
    logging.getLogger("PIL.PngImagePlugin").disabled = True

    # Loading data, and building stats table
    # salmon: pandas.DataFrame = read_salmon_qc(path=snakemake.input.salmon)
    salmon: pandas.DataFrame = read_salmon_auxdir(paths=snakemake.input.salmon)
    general: pandas.DataFrame = read_general_stat_table(
        path=snakemake.input.general, samples_to_keep=salmon.index.tolist()
    )

    stats: pandas.DataFrame = pandas.merge(
        left=salmon,
        right=general,
        left_index=True,
        right_index=True,
        how="left",
    )

    mito: pandas.DataFrame = read_mito_ratio_table(path=snakemake.input.mito)
    stats = pandas.merge(
        left=stats,
        right=mito,
        right_index=True,
        left_index=True,
        how="left",
    )

    rrna: pandas.DataFrame = read_rrna_ratio_table(path=snakemake.input.rrna)
    stats = pandas.merge(
        left=stats,
        right=rrna,
        right_index=True,
        left_index=True,
        how="left",
    )

    # Saving stats table on disk
    stats.to_csv(
        snakemake.output.stats,
        sep=",",
        header=True,
        index=True,
    )
