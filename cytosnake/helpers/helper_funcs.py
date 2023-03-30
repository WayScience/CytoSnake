"""
Module: helper_funcs.py

This module contains helper functions for snakemake workflows.

"""


from typing import Optional
from snakemake.io import expand
from pathlib import Path
from cytosnake.utils.config_utils import load_meta_path_configs
from cytosnake.guards.path_guards import is_valid_path

# loading in config as global variables
PATHS = load_meta_path_configs()


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
    # Barcodes are optional. If not added, set to "None"
    try:
        barcode_path = PATHS["project_dir"]["data_directory_contents"]["barcode"]
    except KeyError:
        barcode_path = None

    return barcode_path


def get_metadata_dir() -> str:
    """Obtains path to the metadata folder with in data folder

    Returns
    -------
    str
        returns string that contains metadata wildcard
    """

    try:
        metadata_dir = PATHS["project_dir"]["data_directory_contents"]["metadata"]
    except KeyError:
        raise KeyError("Unable to find Metadata folder")

    # checking if it is a directory
    if not Path(metadata_dir).is_dir():
        raise TypeError("Metadata file must be a directory not a file")

    return metadata_dir


# ------------------------------
# Output helper functions
# ------------------------------
def cell_count_output() -> str:
    """Generates output path for cell counts;

    Returns
    -------
    str
        path were the cell count data will be produced.
    """
    results_path = Path(PATHS["project_dir_path"]) / "results"
    output_name = "{file_name}_cell_counts"
    ext = "csv"

    # building the string
    return str(results_path / f"{output_name}.{ext}")


def aggregate_output() -> str:
    """Generates output path for cell counts.

    Returns
    -------
    list[str], str
        path to aggregate output file.
    """

    # components of output string
    results_path = Path(PATHS["project_dir_path"]) / "results"
    output_name = "{file_name}_aggregate"
    ext = "csv.gz"

    # constructing file output string
    return str(results_path / f"{output_name}.{ext}")


def annotated_output() -> str:
    """Generates output path for annotated dataset

    Returns
    -------
    str
        path to annotated output file
    """
    # components of output string for annotated profile
    results_path = Path(PATHS["project_dir_path"]) / "results"
    output_name = "{file_name}_annotated"
    ext = "csv.gz"

    # constructing file output string
    return str(results_path / f"{output_name}.{ext}")


def normalized_output() -> str:
    """Generates output path for normalized dataset

    Returns
    -------
    str
        path to normalized output file
    """
    # components of output string for annotated profile
    results_path = Path(PATHS["project_dir_path"]) / "results"
    output_name = "{file_name}_normalized"
    ext = "csv.gz"

    # constructing file output string
    return str(results_path / f"{output_name}.{ext}")


def selected_features_output() -> str:
    """Generates output path for selected features dataset

    Returns
    -------
    str
        path to selected features output file
    """
    # components of output string for annotated profile
    results_path = Path(PATHS["project_dir_path"]) / "results"
    output_name = "{file_name}_features"
    ext = "csv.gz"

    # constructing file output string
    return str(results_path / f"{output_name}.{ext}")


def consensus_output() -> str:
    """Generates output path for consensus  dataset

    Returns
    -------
    str
        path to selected features output file
    """
    # components of output string for annotated profile
    results_path = Path(PATHS["project_dir_path"]) / "results"
    output_name = "consensus_profile"
    ext = "csv.gz"

    # constructing file output string
    return str(results_path / f"{output_name}.{ext}")


def converted_output():
    """Generates output path for parquet profiles

    Returns
    -------
    str
        path to generated parquet files

    """
    data_path = Path(PATHS["project_dir_path"]) / "data"
    output_name = "{file_name}"
    ext = "parquet"

    # constructing file output string
    return str(data_path / f"{output_name}.{ext}")


# ------------------------------
# Formatting I/O functions
# ------------------------------
def get_file_basenames(
    dir_path: Path | str, ext_target: Optional[None | str] = None
) -> list[str]:
    """Obtains all base name of files within a given directory.


    If `ext_target` is set to None, all file basenames within the provided
    directory path will be returned.

    Parameters
    ----------
    dir_path : Path | str
        target directory path

    ext_target : Optional[None | str]
        Allows to obtains file basenames

    Returns
    -------
    list[str]
        list of file base names
    """

    # checking if valid path (Type and if it exists)
    if not is_valid_path(dir_path):
        raise TypeError("Must provide a valid path. Path or str types")

    # -- change to Path type if str type
    if isinstance(dir_path, str):
        dir_path = Path(dir_path)

    # getting all file base names
    glob_query = f"*.{ext_target}"
    if ext_target is None:
        return [_file.stem for _file in dir_path.glob("*")]
    else:
        return [_file.stem for _file in dir_path.glob(glob_query)]
