"""
Module: helper_funcs.py

This module contains helper functions for snakemake workflows to access with. 

"""


from pathlib import Path
from cytosnake.utils.config_utils import load_meta_path_configs


def get_data() -> str:
    """Returns absolute path were the data is located

    Returns
    -------
    str
        returns a string that contains a snakemake wild card that collects all
        plate data

    Raises
    ------
    FileNotFoundError
        Raised if the data folder is not found in the current project folder
    """

    # load in the meta_path configs
    meta_configs = load_meta_path_configs()
    data_dir_path = Path(meta_configs["project_dir"]["data"])
    print(data_dir_path)

    # check if the path exists
    if not data_dir_path.exists():
        raise FileNotFoundError(
            "Unable to find data path in project directory"
        )

    # return snakemake format string path
    return str(data_dir_path / "{file_name}.sqlite")
