# coding: utf-8

import logging
import pandas
import numpy
import seaborn
import yaml
import matplotlib.pyplot

from pathlib import Path
from sklearn.decomposition import PCA


def read_counts(path: str) -> pandas.DataFrame:
    """
    Return loaded counts as a pandas DataFrame

    Parameters:
    path    (str): Path to count file

    Return:
    (pandas.DataFrame): Loaded counts
    """
    # Load from file
    counts: pandas.DataFrame = pandas.read_csv(
        filepath_or_buffer=path,
        sep="\t",
        header=0,
        index_col=0,
    )
    logger.debug(counts.shape)

    # Keeping numeric columns
    counts = counts.select_dtypes(include=numpy.number)
    logger.debug(counts.shape)

    # Dropping NA counts
    counts.dropna(axis=1, how="all", inplace=True)
    logger.debug(counts.shape)

    # Drop null rows
    counts = counts.loc[counts.sum(axis=1) != 0]
    logger.debug(counts.shape)

    logger.debug(counts.head())
    return counts


def perform_pca(counts: pandas.DataFrame) -> dict[str, pandas.DataFrame, numpy.ndarray]:
    """
    Return a fitted PCA object

    Parameters:
    counts  (pandas.DataFrame): Loaded and filtered counts

    Return:
    """
    # Prepare PCA
    nbc: int = min(counts.shape)
    skpca = PCA(n_components=nbc)

    # Perform and fit PCA
    sktransform = skpca.fit_transform(counts.T)
    skvar = skpca.explained_variance_ratio_
    results = pandas.DataFrame(
        sktransform,
        columns=[f"PC{i}" for i in range(1, nbc + 1, 1)],
        index=counts.columns.tolist(),
    )

    # Return results
    return {
        "pca": results,
        "variance": skvar,
    }


def perform_correlation(pca: pandas.DataFrame) -> pandas.DataFrame:
    """
    Perform dataframe correlation

    Parameters:
    pca (pandas.DataFrame): PCA results
    """
    correlation: pandas.DataFrame = pca.corr(method="pearson", numeric_only=True)
    correlation.columns = pca.index
    correlation.index = pca.index
    return correlation


def corr_to_yaml(corr: pandas.DataFrame) -> dict[str, str | dict[str, str | float]]:
    """
    Return MultiQC heatmap configuration as a dictionary

    Parameters:
    corr    (padnas.DataFrame): PCA correlation results

    Return:
    (dict[str, str | dict[str, str | float]]): MultiQC config
    """
    logging.info("Building PCA correlation for MultiQC.")
    return {
        "parent_id": "In_house_QC",
        "parent_name": "In-house quality controls",
        "parent_description": "This section contains PCA and heatmaps built through in-house scripts.",
        "id": "Sample_correlation",
        "section_name": "Heatmap of sample correlation",
        "description": "Pearson correlation between samples.",
        "plot_type": "heatmap",
        "pconfig": {
            "id": "sample_heatmap",
            "title": "Heatmap of sample correlation",
        },
        "data": corr.to_dict(orient="index"),
    }


def pca_to_yaml(pca: pandas.DataFrame) -> dict[str, str | dict[str, str | float]]:
    """
    Return MultiQC scatterplot configuration as a dictionary

    Parameters:
    pca (pandas.DataFrame): PCA results

    Return:
    (dict[str, str | dict[str, str | float]]): MultiQC config
    """
    logging.info("Building PCA section for MultiQC.")
    return {
        "parent_id": "In_house_QC",
        "id": "Sample_PCA",
        "section_name": "Sample PCA",
        "description": "Axes 1 & 2 of the PCA among samples.",
        "plot_type": "scatter",
        "pconfig": {
            "id": "pca_scatter_plot",
            "title": "PCA plot",
            "xlab": "PC1",
            "ylab": "PC2",
        },
        "data": pca.to_dict(orient="index"),
    }


def plot_loadings(variance: numpy.ndarray, out_png: str | Path) -> None:
    """
    Save PCA explained variance as a PNG file
    """
    logging.info(f"Building scree plot and saving it to {out_png=}")
    var_df = pandas.DataFrame(variance, columns=["Loadings"])
    var_df["PCA_Axe"] = [f"PC{i}" for i in range(1, len(var_df) + 1)]
    var_df = var_df.head(min(20, len(var_df)))

    seaborn.set()
    seaborn.barplot(var_df, y="Loadings", x="PCA_Axe")

    matplotlib.pyplot.plot(
        range(0, len(variance)),
        numpy.cumsum(variance),
        c="red",
        label="Cumulative Explained Variance",
    )

    matplotlib.pyplot.legend(loc="upper left")
    matplotlib.pyplot.xlabel("Number of components")
    matplotlib.pyplot.ylabel("Explained variance (eignenvalues)")
    matplotlib.pyplot.title("Scree plot")

    matplotlib.pyplot.savefig(
        out_png,
        bbox_inches="tight",
    )

    matplotlib.pyplot.cla()
    matplotlib.pyplot.clf()
    matplotlib.pyplot.close()


def plot_pca(pca: pandas.DataFrame, out_png: str | Path) -> None:
    """
    Save first 2 PCA axes on disk as a PNG
    """
    logging.info(f"Building and saving pca to {out_png=}")
    print(pca)
    seaborn.set()

    seaborn.relplot(
        data=pca,
        x="PC1",
        y="PC2",
    )
    seaborn.rugplot(
        data=pca,
        x="PC1",
        y="PC2",
    )

    matplotlib.pyplot.title("First axes of the PCA")
    matplotlib.pyplot.savefig(
        out_png,
        bbox_inches="tight",
    )
    matplotlib.pyplot.cla()
    matplotlib.pyplot.clf()
    matplotlib.pyplot.close()


def main(
    salmon_raw_counts: str | Path,
    out_loadings_png: str | Path,
    out_pca_png: str | Path,
    out_corr_yaml: str | Path,
    out_pca_yaml: str | Path,
) -> None:
    """
    Perform correlation and PCA from Salmon counts
    """
    salmon: pandas.DataFrame = read_counts(path=salmon_raw_counts)

    pca = perform_pca(counts=salmon)

    correlation: pandas.DataFrame = perform_correlation(pca=pca["pca"])

    corr_yaml = corr_to_yaml(corr=correlation)
    with open(out_corr_yaml, "w") as out_corr_stream:
        yaml.dump(
            corr_yaml,
            out_corr_stream,
            default_flow_style=False,
        )

    pca_yaml = pca_to_yaml(pca=pca["pca"])
    with open(out_pca_yaml, "w") as out_pca_stream:
        yaml.dump(
            pca_yaml,
            out_pca_stream,
            default_flow_style=False,
        )

    plot_loadings(
        variance=pca["variance"],
        out_png=out_loadings_png,
    )
    plot_pca(
        pca=pca["pca"],
        out_png=out_pca_png,
    )


if __name__ == "__main__":
    logger = logging.getLogger(__name__)
    logging.basicConfig(
        level=logging.DEBUG,
        filename=snakemake.log[0],
        filemode="w",
    )
    logging.getLogger("matplotlib.font_manager").disabled = True
    logging.getLogger("PIL.PngImagePlugin").disabled = True

    main(
        salmon_raw_counts=snakemake.input[0],
        out_loadings_png=snakemake.output["scree_png"],
        out_pca_png=snakemake.output["pca_png"],
        out_corr_yaml=snakemake.output["corr_yaml"],
        out_pca_yaml=snakemake.output["pca_yaml"],
    )
