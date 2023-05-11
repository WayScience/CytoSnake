"""
module: datapath.py

This module contains fuctions that build stirng path that are snakemake compatible.
"""

import pathlib
from typing import Optional

import cytosnake.utils.cyto_paths as cyto_paths
from cytosnake.guards.path_guards import is_valid_path
from cytosnake.utils.config_utils import load_general_configs, load_meta_path_configs

# constant wildcard name
WILDCARD_BASENAME = "{basename}"

# loading in config files
META_PATH_CONFIGS = load_meta_path_configs()
GENERAL_CONFIGS = load_general_configs()


# -----------------------
# input paths from configs
# -----------------------
def get_barcodes() -> str:
    """Obtains path to the bardcode file within data folder

    Returns
    -------
    str
        return string with barcode name wildcard

    """
    # Barcodes are optional. If not added, set to "None"
    try:
        barcode_path = META_PATH_CONFIGS["project_dir"]["data_directory_contents"][
            "barcode"
        ]
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
        metadata_dir = META_PATH_CONFIGS["project_dir"]["data_directory_contents"][
            "metadata"
        ]
    except KeyError as e:
        raise KeyError("Unable to find Metadata folder") from e

    # checking if it is a directory
    if not pathlib.Path(metadata_dir).is_dir():
        raise TypeError("Metadata file must be a directory not a file")

    return metadata_dir


# --------------------------
# Main path builder function
# --------------------------
def build_path(data_type: str, use_converted: bool) -> str:
    """Builds a output string path pointing to a specific dataset

    data_type : str
        type of dataset string path to be build.

    use_converted : bool
        indicating to use converted plate dataset. However if this is True and
        data_type != 'plate_data' then the `use_converted` logic will be ignored
        since it is only meant to be used for plate data

    Returns
    -------
    str
        string path pointing to dataset

    Raises
    ------
    ValueError
        Raised if the `data_type` is not found within configs.
    """

    # loading all data configs
    data_configs = GENERAL_CONFIGS["data_configs"]

    # loading in the config
    if data_type not in data_configs.keys():
        raise ValueError(f"`{data_type}` is not a supported dataset")

    # loading in data_type specific configs
    selected_datatype = data_configs[data_type]

    # loading default path components into varaibles
    header = cyto_paths.get_results_dir_path()
    suffix = selected_datatype["suffix"]
    compression = selected_datatype["compression_ext"]
    ext = selected_datatype["ext"]

    # edit ext if compression is used
    if compression is not None:
        ext = f"{ext}.{compression}"

    # plate_data requires some of the defaults to be change since it is in a different
    # folder
    if data_type == "plate_data":
        header = cyto_paths.get_project_root() / "data"

        # this will determine if users want to use their converted dataset
        if use_converted:
            ext = selected_datatype["converted_ext"]

            # check if there are any converted files
            # if there's none, default to inital plate data type extension
            n_converted = list(len(header.glob(f"*.{ext}")))
            if n_converted == 0:
                ext = selected_datatype["ext"]

    # building path
    return f"{header}/{WILDCARD_BASENAME}_{suffix}.{ext}"


# --------------------
# Helper i/o functions
# --------------------
def get_all_basenames(
    dirpath: pathlib.Path, ext_target: Optional[str | None] = None
) -> list[str]:
    """Gets the basenames of all files from a given directory. If there's a need for a
    more targeted search, use the `target_ext` which will get all the file basenames
    that have that extension.

    Parameters
    ----------
    dirpath : str
        path to directory to be searched
    target_ext : Optional[str  |  None], optional
        allows a more targeted search if given an extension, by default None

    Returns
    -------
    list[str]
        list of file basenames
    """

    # check if the path given is a directory path
    if not is_valid_path():
        raise TypeError("`dirpath` must be a valid path")
    if not dirpath.is_dir():
        raise NotADirectoryError("`dirpath` must be a directory, not a file")

    # getting all file base names
    glob_query = f"*.{ext_target}"
    if ext_target is None:
        return [_file.stem for _file in dirpath.glob("*")]
    else:
        return [_file.stem for _file in dirpath.glob(glob_query)]
