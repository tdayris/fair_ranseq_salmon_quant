"""Snakemake-wrapper for building new genome file"""

__author__ = "Thibault Dayris"
copyright__ = "Copyright 2024, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import os
import pandas


def replace_col(df: pandas.DataFrame, column: str, default: str) -> pandas.DataFrame:
    """
    Replace null values and unreachable paths 
    in a given column with provided defaults.

    Parameters:
    df      (pandas.DataFrame)  : Existing genomes information
    column  (str)               : The column to modify
    default (str)               : The default value in case of missing value or column

    Return  (pandas.DataFrame)  : Updated genomes information
    """
    datatype: str = "dna"

    if column in df.columns:
        # Case user did (partially?) filled the column of interest
        new_values: list[str] = []

        df_iterator = zip(df.species, df.build, df.release, df[column], strict=False)

        for species, build, release, value in df_iterator:
            if value and os.path.exists(value):
                # User value was not correct
                new_values.append(os.path.realpath(value))
            else:
                # User value is correct and should be kept
                new_values.append(
                    os.path.realpath(
                        default.format(
                            species=species,
                            build=build,
                            release=release,
                            datatype=datatype,
                        )
                    )
                )

    else:
        # Case user did not filled the column of interest
        df_iterator = zip(df.species, df.build, df.release)

        df[column] = [
            os.path.realpath(
                default.format(
                    species=species, build=build, release=release, datatype=datatype
                )
            )
            for species, build, release in df_iterator
        ]

    return df


# Main programm
genomes: pandas.DataFrame = snakemake.params.genomes

# Bowtie2 index
genomes = replace_col(
    df=genomes,
    column="salmon_index",
    default="reference/salmon_index/{species}.{build}.{release}/{species}.{build}.{release}/",
)
genomes.to_csv(path_or_buf=snakemake.output.genomes, sep=",", header=True, index=False)
