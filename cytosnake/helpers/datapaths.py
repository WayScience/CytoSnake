"""
module: datapath.py

This module contains fuctions that build stirng path that are snakemake compatible.
"""

from cytosnake.utils.config_utils import load_general_configs, load_meta_path_configs

# loading in config files
META_PATH_CONFIGS = load_meta_path_configs()
GENERAL_CONFIGS = load_general_configs()


def build_path(data_type: str) -> str:
    """Builds a string path pointing to a specific dataset

    data_type : str
        type of dataset string path to be build.

    Returns
    -------
    str
        string path pointing to dataset
    """

    # loading all data configs
    data_configs = GENERAL_CONFIGS["data_config"]

    # loading in the config
    if data_type not in data_configs.keys():
        raise ValueError(f"`{data_type}` is not a supported dataset")

    # loading path components into varaibles
    header = None
    basename_wildcard = "{basename}"
    suffix = None
    ext = None
    compression = None

    # edit ext if compression is used
    if compression is not None:
        ext = f"{ext}.{compression}"

    # building path
    return f"{header}/{basename_wildcard}_{suffix}.{ext}"
