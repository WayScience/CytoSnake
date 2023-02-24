"""
Module: helper_funcs.py

This module contains helper functions for snakemake workflows.

"""


from pathlib import Path

from cytosnake.utils.config_utils import load_meta_path_configs

# loading in config as global variables
PATHS = load_meta_path_configs()
RESULTS_DIR = "../results"


# ------------------------------
# Inputs helper functions
# ------------------------------
def get_data_folder() -> Path:
    """Returns absolute path that points to the data folder.

    Returns
    -------
    Path
        path to data folder

    Raises
    ------
    KeyError
        Raised if the `data` path is not found in the `_paths.yaml`

    """
    meta_configs = load_meta_path_configs()
    try:
        data_dir_path = Path(meta_configs["project_dir"]["data"])
    except KeyError:
        raise KeyError("Unable to find data folder path in project_dir config")

    return data_dir_path


def get_plate_data() -> str:
    """Returns absolute path were the data are located

    Returns
    -------
    str
        returns a string that contains a snakemake wild card that collects all
        plate data

    """
    # load in the meta_path configs
    data_dir_path = get_data_folder()
    return str(data_dir_path / "{file_name}.sqlite")


def get_barcodes() -> str:
    """Obtains path to the bardcode file within data folder

    Returns
    -------
    str
        return string with barcode name wildcard

    """
    data_dir_path = get_data_folder()
    return str(data_dir_path / "{barcode_name}.txt")


def get_metadata_dir() -> str:
    """Obtains path to the metadata folder with in data folder

    Returns
    -------
    str
        returns string that contains metadata wildcard
    """
    data_dir_path = get_data_folder()
    return str(data_dir_path / "{metadata}")


# ------------------------------
# Output helper functions
# ------------------------------
def cell_count_output() -> str:
    """Generates output path for cell counts;

    Returns
    -------
    Path
        path were the cell count data will be produced.
    """
    results_path = Path(PATHS["projet_dir_path"]) / "results"
    output_name = "cell_counts"
    ext = ".csv"

    # building the string
    return str(results_path / f"{output_name}.{ext}")
