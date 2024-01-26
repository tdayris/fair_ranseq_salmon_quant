# codinf: utf-8

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2024, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import logging
import yaml

from typing import Any


def write_config(data: dict[str, Any], path: str) -> None:
    """
    Save dictionary in a yaml-formatted file

    Parameters:
    data    (dict[str, Any]): datavzrd configuration
    path    (str)           : path to output yaml file

    Return: None
    """
    logging
    with open(path, "w") as yaml_stream:
        yaml.dump(data=data, stream=yaml_stream, default_flow_style=False)


def fill_dict(
    targets: str, counts: str, table_path: str, table_description: str
) -> dict[str, Any]:
    """
    Build datavzrd dictionary

    Parameters:
    targets             (str): Gene or Transcripts
    counts              (str): Raw counts, or TPM
    table_path          (str): Path to the data table
    table_description   (str): Header text in the HTML index table

    Return (dict[str, Any]):
    Datavzrd configuration
    """
    return {
        "name": f"Quantification of {targets}, units: {counts} counts",
        "datasets": {
            "salmon": {
                "path": table_path,
                "separator": "\t",
                "headers": 1,
                "offer-excel": False,
            },
        },
        "views": {
            "salmon": {
                "dataset": "salmon",
                "desc": table_description,
            },
        },
    }


table_description: str = f"""
# Salmon quantification

Please find below the Salmon quantification of {snakemake.wildcards.targets}, using {snakemake.wildcards.counts} counts.
"""

if __name__ == "__main__":
    logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)
    try:
        datavzrd_configuration: dict[str, Any] = fill_dict(
            sample_name=snakemake.wildcards.sample,
            table_path=snakemake.input.table,
            table_description=table_description,
        )
        write_config(data=datavzrd_configuration, path=snakemake.output.yaml)
    except Exception as e:
        logging.error(e)
        raise e
