"""
module: datapath.py


The datapath module contains functions that construct string paths that are snakemake
compatible.

Datapaths obtain pathing information from the `cyto_paths` modules, which contains
pathing information after executing the `init` mode. With this information, the
datapath module ensures that each path produced are snakemake compatible.
"""

import pathlib
from collections import defaultdict
from typing import Optional

from cytosnake.guards.path_guards import is_valid_path
from cytosnake.utils import cyto_paths
from cytosnake.utils.config_utils import load_general_configs, load_meta_path_configs

# loading in config files
META_PATH_CONFIGS = load_meta_path_configs()
GENERAL_CONFIGS = load_general_configs()


# -----------------------
# input paths from configs
# -----------------------
def get_barcodes() -> str | None:
    """Obtains path to the bardcode file within data folder. Returns None if no barcode
    was provided in the `init` mode.

    Returns
    -------
    str | None
        return string with barcode name wildcard, else None
    """
    # create a default dictionary that contains default value if key is not found
    barcode_path = defaultdict(lambda: None)

    # update default dict with loaded parameters
    barcode_path |= META_PATH_CONFIGS["project_dir"]["data_directory_contents"]

    # returns None if not found, else the absolute barcode path
    return barcode_path["barcode"]


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
def build_path(input_type: str, use_converted: Optional[bool] = False) -> str:
    """Builds an output string path pointing to a specific dataset

    input_type : str
        type of dataset string path to be build.

    use_converted : Optional[bool]
        indicating to use converted plate dataset. However if this is True and
        data_type != 'plate_data' then the `use_converted` logic will be ignored
        since it is only meant to be used for plate data. (default=False)

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
    if input_type not in data_configs["data_types"].keys():
        raise ValueError(f"`{input_type}` is not a supported dataset")

    # loading in data_type specific configs
    selected_datatype = data_configs["data_types"][input_type]

    # loading default path components into variables
    header = cyto_paths.get_results_dir_path()
    suffix = f"_{selected_datatype['suffix']}"
    compression = selected_datatype["compression_ext"]
    ext = selected_datatype["file_ext"]

    # edit ext if compression is used
    if compression is not None:
        ext = f"{ext}.{compression}"

    # plate_data requires some of the defaults to be change since it is in a different
    # folder
    if input_type == "plate_data":
        header = cyto_paths.get_project_root() / "data"
        suffix = ""

        # this will determine if users want to use their converted dataset
        if use_converted:
            ext = selected_datatype["converted_ext"]

    # building path
    return f"{header}/{{basename}}{suffix}.{ext}"


# --------------------
# Helper i/o functions
# --------------------
def get_all_basenames(
    dirpath: str | pathlib.Path, ext_target: Optional[str | None] = None
) -> list[str]:
    """Gets the basenames of all files from a given directory. If there's a need for a
    more targeted search, use the `target_ext` which will get all the file basenames
    that have that extension.

    Parameters
    ----------
    dirpath : str | pathlib.Path
        path to directory to be searched
    target_ext : Optional[str  |  None], optional
        allows a more targeted search if given an extension, by default None

    Returns
    -------
    list[str]
        list of file basenames
    """

    # check if the path given is a directory path
    if isinstance(dirpath, str):
        dirpath = pathlib.Path(dirpath)
    if not is_valid_path(dirpath):
        raise TypeError("`dirpath` must be a valid path")
    if not dirpath.is_dir():
        raise NotADirectoryError("`dirpath` must be a directory, not a file")

    # getting all file base names
    if ext_target is None:
        return [_file.stem for _file in dirpath.glob("*")]

    # if ext_target is provided
    return [_file.stem for _file in dirpath.glob(f"*.{ext_target}")]
