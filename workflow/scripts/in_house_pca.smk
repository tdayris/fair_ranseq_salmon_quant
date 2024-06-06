# coding: utf-8

import pandas
import numpy
import sklearn
import yaml


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
    counts = counts[counts.select_dtypes(['number']) or counts.select_dtypes([numpy.number]]
    logger.debug(counts.shape)

    # Dropping NA counts
    counts.dropna(axis=1, how='all', inplace=True)
    logger.debug(counts.shape)

    # Drop null counts
    counts = counts[~ counts.isull().all(axis=1)]
    logger.debug(counts.shape)

    logger.debug(counts.head())
    return counts


def perform_pca(counts: pandas.DataFrame) -> sklearn.decomposition.PCA:
    """
    Return a fitted PCA object

    Parameters:
    counts  (pandas.DataFrame): Loaded and filtered counts

    Return:
    """
    # Prepare PCA
    nbc: int = min(counts.shape)
    skpca = sklearn.decomposition.PCA(n_components=nbc)

    # Perform and fit PCA
    sktransform = skpca.fit_transform(counts.T)
    skvar = skpca.explained_variance_ratio_
    results = pandas.DataFrame(
        sktransform,
        columns=[f"PC{i}" for i in range(1, nbc+1, 1)],
        index=counts.columns.tolist()
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
    return pca.corr(method="pearson", numeric_only=True)

def corr_to_yaml(corr: pandas.DataFrame) -> dict[str, str | dict[str, str | float]]:
    """
    Return MultiQC heatmap configuration as a dictionary

    Parameters:
    corr    (padnas.DataFrame): PCA correlation results

    Return:
    (dict[str, str | dict[str, str | float]]): MultiQC config
    """
    return {
        "parent_id": "In_house_QC",
        "parent_name": "In-house quality controls",
        "parent_description": "This section contains PCA and heatmaps built through in-house scripts."
        "id": "Sample_correlation",
        "section_name": "Heatmap of sample correlation",
        "description": "Pearson correlation between samples.",
        "plot_type": "heatmap",
        "pconfig": {
            "id": "sample_heatmap",
            "title": "Heatmap of sample correlation",
        },
        "data": corr.to_dict(),
    }

def pca_to_yaml(pca: pandas.DataFrame) -> dict[str, str | dict[str, str | float]]:
    """
    Return MultiQC scatterplot configuration as a dictionary

    Parameters:
    pca (pandas.DataFrame): PCA results

    Return:
    (dict[str, str | dict[str, str | float]]): MultiQC config
    """
    return {
        "parent_id": "In_house_QC",
        "id": "Sample_PCA",
        "section_name": ""
    }

    
    
