#!/usr/bin/python3.8
# conding: utf-8

"""
Join multiple Salmon files, rename columns and
add genes name if annotation is provided
"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2020, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


import logging
import pandas
import numpy

from os.path import dirname, basename
from snakemake.utils import makedirs


def read_tx2gene(
    path: str, genes: bool = False, header: bool = False
) -> pandas.DataFrame:
    """
    This function reads a TSV file containing the following columns:

    1: ens_gene
    2: ens_transcript
    3: ext_gene

    And returns it as a DataFrame
    """
    if str(snakemake.wildcards.counts).lower().startswith("gene"):
        t2g: pandas.DataFrame = pandas.read_csv(
            path,
            sep="\t",
            index_col=None,
            header=None,
            dtype=str,
        )
        t2g.columns = [
            "Ensembl_Gene_ID", 
            "Gene_Name",
        ]
        t2g.drop_duplicates(inplace=True)
        t2g.set_index("Ensembl_Gene_ID",  inplace=True)
    else:
        t2g: pandas.DataFrame = pandas.read_csv(
            path,
            sep=",",
            index_col=None,
            header=None,
            dtype=str,
        )
        t2g.columns = ["Ensembl_Gene_ID", "Ensembl_Transcript_ID", "Gene_Name"]
        t2g.drop_duplicates(inplace=True)
        t2g.set_index("Ensembl_Transcript_ID",  inplace=True)

    return t2g


def read_salmon(path: str) -> pandas.DataFrame:
    """
    This function reads a single salmon quant.sf or quant.genes.sf
    and returns it as a pandas DataFrame
    """
    df: pandas.DataFrame = pandas.read_csv(
        path,
        sep="\t",
        index_col=0,
        header=0,
        na_values="",
    )

    if snakemake.params.get("drop_patch", False) is True:
        df.index = [i.split(".")[0] for i in df.index.tolist()]

    return df


logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)

merged_frame: pandas.DataFrame | None = None

for quant in snakemake.input["quant"]:
    logging.debug(f"Reading {quant}")
    data: pandas.DataFrame = read_salmon(quant)

    logging.debug("Cleaning dataframe")
    sample_id: str = basename(dirname(quant))
    if len(suffix := snakemake.params.get("suffix", "")) > 0:
        sample_id = sample_id[: -len(suffix)]
    if len(prefix := snakemake.params.get("prefix", "")) > 0:
        sample_id = sample_id[len(prefix) :]

    data = data[[snakemake.params.get("column", "TPM")]]
    data.columns = [sample_id]

    logging.debug("Merging dataframe")
    try:
        merged_frame = pandas.merge(
            merged_frame, data, left_index=True, right_index=True
        )
    except TypeError:
        merged_frame = data

merged_frame.fillna(0)
logging.debug("Merged quant frame:")
logging.debug(merged_frame.head())

if snakemake.params.get("gencode", False) is True:
    logging.debug("Removing gencode patch ids")
    merged_frame = merged_frame.set_index(
        pandas.DataFrame(merged_frame.index.str.split(".").tolist())[0]
    )

if (fillna := snakemake.params.get("fillna", None)) is not None:
    merged_frame.fillna(fillna, inplace=True)

if snakemake.params.get("drop_null", False) is True:
    logging.debug("Removing null values")
    merged_frame = merged_frame.loc[~(merged_frame == 0).all(axis=1)]

if snakemake.params.get("drop_na", False) is True:
    logging.debug("Removing NA values")
    merged_frame.dropna(axis=0, how="all", inplace=True)

if (tr2gene_path := snakemake.input.get("tx2gene", None)) is not None:
    logging.debug("Adding gene names")

    t2g = read_tx2gene(
        tr2gene_path,
        snakemake.params.get("genes", False),
        snakemake.params.get("header", False),
    )
    logging.debug("tx2gene:")
    logging.debug(t2g.head())

    merged_frame = pandas.merge(
        merged_frame, t2g, left_index=True, right_index=True, how="left"
    )

if snakemake.params.get("genes", False) is True:
    merged_frame.index.name = "Ensembl_Gene_ID"
else:
    merged_frame.index.name = "Ensembl_Transcript_ID"


logging.debug("Saving DataFrame to disk")
merged_frame.to_csv(
    snakemake.output["tsv"],
    sep="\t",
    index=True,
    header=True,
    index_label=(
        "target_id" if snakemake.params.get("index_label", False) is True else False
    ),
)
