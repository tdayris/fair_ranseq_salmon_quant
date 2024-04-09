# coding: utf-8

import logging
import matplotlib.pyplot as plt
import pandas
import seaborn


def load_counts(path: str) -> pandas.DataFrame:
    """
    Load Salmon counts in memory

    Parameters:
    path (str): Path to counts

    Return:
    (pandas.DataFrame): Loaded counts
    """
    logging.info(f"Loading quantification from {path=}")
    tpm: pandas.DataFrame = pandas.read_csv(path, sep="\t", header=0, index_col=0)
    logging.debug(tpm.shape)
    return tpm


def get_mito_counts(
    counts: pandas.DataFrame, reverse: bool = False
) -> pandas.DataFrame:
    """
    Extract mitochondrial counts in all quantified samples

    Parameters:
    counts  (pandas.DataFrame)  : Complete count table
    reverse (bool)              : Reverse the search

    Return:
    (pandas.DataFrame): The extracted counts
    """
    if reverse:
        logging.info("Extracting *not* mitochondial counts...")
        moti = counts.loc[~counts.Gene_Name.str.upper().str.startswith("MT-", na=False)]
    else:
        logging.info("Extracting mitochondial counts...")
        moti = counts.loc[counts.Gene_Name.str.upper().str.startswith("MT-", na=False)]

    logging.debug(moti.shape)
    return moti


def get_rrna_counts(
    counts: pandas.DataFrame, reverse: bool = False
) -> pandas.DataFrame:
    """
    Extract mitochondrial counts in all quantified samples

    Parameters:
    counts  (pandas.DataFrame)  : Complete count table
    reverse (bool)              : Reverse the search

    Return:
    (pandas.DataFrame): The extracted counts
    """
    if reverse:
        logging.info("Extracting *not* rRNA counts...")
        rrna = counts.loc[~counts.Gene_Name.str.lower().str.contains("rrna")]
    else:
        logging.info("Extracting rRNA counts...")
        rrna = counts.loc[counts.Gene_Name.str.lower().str.contains("rrna")]

    logging.debug(rrna.shape)
    return rrna


def compute_sum_expression(counts: pandas.DataFrame) -> pandas.DataFrame:
    """
    Compute sum expression, knowing that a NaN sums there is
    no red count at all on the gene.

    Parameters:
    couts (pandas.DataFrame) : The count table to compute per-gene sum expression


    Return:
    (pandas.DataFrame): The per-gene sum expression
    """
    logging.info("Computing sum expression...")
    sums: pandas.DataFrame = counts.fillna(0).sum(numeric_only=True)
    sums.columns = ["Sample_id", "Expression"]
    logging.debug(sums.head())
    return sums


def compute_ratios(
    numerator: pandas.DataFrame, denominator: pandas.DataFrame
) -> pandas.DataFrame:
    """
    Compute expression ratios, knowing that a missing key
    sums there was no read counts

    Parameters:
    numerator       (pandas.DataFrame): The rrna/mito counts
    denominator     (pandas.DataFrame): The non-rrna or non-mito counts

    Return:
    (pandas.DataFrame): The computed ratios
    """
    logging.info("Computing ratios...")
    div: pandas.DataFrame = numerator.div(other=denominator, axis="rows", fill_value=0)
    div = pandas.concat([div, numerator, denominator], axis=1)
    div.sort_index(inplace=True)
    div.reset_index(inplace=True)
    div.columns = ["Sample_id", "Ratio", "numerator", "denominator"]
    logging.debug(div.head())
    return div


def to_png(data: pandas.DataFrame, title: str, output: str) -> None:
    """
    Save input data as a barplot to the given file

    Parameters:
    data    (pandas.DataFrame)  : Data to plot
    title   (str)               : Graph title
    output  (str)               : Path to output PNG
    """
    logging.info(f"Saving graph to {output=}")
    plt.figure(figsize=(1.024, 2.048), dpi=100)
    seaborn.set_theme(style="whitegrid")
    seaborn.set(font_scale=0.55)

    logging.debug(data.head())
    logging.debug(data.columns)
    g: seaborn.FacetGrid = seaborn.catplot(
        data=data,
        kind="bar",
        y="Sample_id",
        x="Ratio",
    )
    g.despine(left=True)
    g.set_axis_labels(title, "Samples")

    plt.savefig(output, bbox_inches="tight", dpi=1000)
    plt.close()


if __name__ == "__main__":
    logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)

    logging.getLogger("matplotlib.font_manager").disabled = True
    logging.getLogger("PIL.PngImagePlugin").disabled = True

    counts: pandas.DataFrame = load_counts(path=snakemake.input.table)
    total_sum: pandas.DataFrame = compute_sum_expression(counts=counts)

    # Working on mito counts
    mito: pandas.DataFrame = get_mito_counts(counts=counts)
    mito = compute_sum_expression(counts=mito)

    mito_ratio: pandas.DataFrame = compute_ratios(numerator=mito, denominator=total_sum)
    mito_ratio.columns = ["Sample_id", "Ratio", "mito_counts", "Total_counts"]
    mito_ratio.to_csv(snakemake.output.mito_csv, sep=",", header=True, index=False)
    to_png(
        data=mito_ratio,
        title="Mitocondrial RNA ratio",
        output=snakemake.output.mito_png,
    )

    # Working on rrna counts
    rrna: pandas.DataFrame = get_rrna_counts(counts=counts)
    rrna = compute_sum_expression(counts=rrna)

    rrna_ratios: pandas.DataFrame = compute_ratios(
        numerator=rrna, denominator=total_sum
    )
    rrna_ratios.columns = ["Sample_id", "Ratio", "rRNA_counts", "Total_counts"]
    rrna_ratios.to_csv(snakemake.output.rrna_csv, sep=",", header=True, index=False)
    to_png(data=rrna_ratios, title="rRNA ratio", output=snakemake.output.rrna_png)
