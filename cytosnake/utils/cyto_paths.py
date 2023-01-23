"""
module: cyto_paths.py

This module will contain functions that handles `cytosnake's` pathing
"""

from typing import Optional
from pathlib import Path
import cytosnake

# cytosnake imports
from cytosnake.utils.file_utils import file_search
from cytosnake.guards.path_guards import is_valid_path
from cytosnake.utils.config_utils import load_meta_path_configs


def get_meta_path() -> Path:
    """returns meta path configurational file path.

    Returns
    -------
    Path
        Path object pointing to `_paths.yaml` config file
    """

    # getting project directory
    proj_dir = get_project_root() / ".cytosnake"
    if not proj_dir.exists():
        raise FileNotFoundError("Current directory is not a project folder")

    return proj_dir.absolute()


def is_cytosnake_dir(dir_path: Optional[str | Path] = None) -> bool:
    """Checks if the current directory has been set up for cytosnake. Searches
    for the `.cytosnake` file in current or specified directory.

    Parameters
    ----------
    dir_path : Optional[str | Path]
        Path to directory. If None, current directory will be used.
        [Default=None]

    Returns
    -------
    bool
        True if directory has been initialized for cytosnake, else False
    """

    # setting to current directory if not path is specified
    if dir_path is None:
        dir_path = Path().absolute()

    # type checking
    if not is_valid_path(dir_path):
        _type = type(dir_path)
        raise ValueError(
            f"dir_path must be either str or Path type, not: {_type}"
        )
    elif isinstance(dir_path, str):
        dir_path = Path(str).absolute()

    # get project older
    cyto_proj_dir = dir_path / ".cytosnake"
    return bool(cyto_proj_dir.exists())


def get_cytosnake_package_path() -> Path:
    """Returns paths where the package is installed

    Return
    ------
    Path
        Returns absolute path of where the package is installed

    Raises
    ------
    FileNotFoundError
        Raised if the `.git` folder is not found.
    """

    # get location of this file
    # -- check if he ".git" folder exists
    project_path = Path(cytosnake.__path__[0]).parent
    git_path = project_path / ".git"
    if not git_path.exists:
        raise FileNotFoundError("Unable to find cytosnake package path")

    return project_path


def get_project_root() -> Path:
    """Returns complete path where the analysis is taken place. The function
    will check if `.cytosnake` folder exists, if not an error will be raised.

    Return
    ------
    Path
        Returns absolute path of the project folder

    Raises
    ------
    FileNotFoundError
        Raised with the current directory is not a project folder
    """

    # get current working directory
    project_dir = Path().absolute()

    # check if the `.cytosnake` folder exist
    project_folder = project_dir / ".cytosnake"
    if not project_folder.exists():
        raise FileNotFoundError("Current directory is not a project folder")

    return project_dir


def get_workflow_fpaths() -> dict:
    """Obtains all file paths located in the `workflows` folder as a dictionary.

    Returns
    -------
    dict
        Structured dictioanry containing directory names and paths are key value
        pairs

    Raises
    ------
    FileNotFoundError
        Raised if the file workflows directory is not found
    """

    # find workflow path in project directory
    proj_root_path = get_project_root()
    workflow_path = proj_root_path / "workflows"
    if not is_valid_path(workflow_path):
        raise FileNotFoundError("Unable to find workflow directory")

    return file_search(workflow_path)


def get_config_fpaths() -> dict:
    """Obtains all file paths located in the `configs` folder as a dictionary.

    Returns
    -------
    dict
        structured dictionary directory name and file paths as key value pairs
    """
    proj_root_path = get_project_root()
    config_path = proj_root_path / "configs"
    if not is_valid_path(config_path):
        raise FileNotFoundError("Unable to find workflow directory")

    return file_search(config_path)


def get_project_dirpaths() -> dict:
    """returns a dictionary containing directory name and path as key value
    pairs.

    Returns
    -------
    dict
        directory name and path as key value pairs.
    """

    # getting project root folder
    proj_root_path = get_project_root()

    # finding all files
    list_of_files = proj_root_path.glob("*")

    # creating dict containing {dir_name : dir_path}
    # -- ignore unwanted directories
    ignore_dirs = [".vscode", ".git", "CytoSnake.egg-info"]

    # -- collect dir names and paths
    all_dirs = {}
    for _file in list_of_files:
        if _file.name in ignore_dirs:
            continue
        elif _file.is_dir():
            all_dirs[_file.name] = str(_file.absolute())

    return all_dirs


def get_data_path() -> Path:
    """Returns path pointing where the data folder is

    Returns
    -------
    Path
        Path to data folder in project directory
    """

    # load in `_paths.yaml` meta data
    loaded_meta_paths = load_meta_path_configs()
    return Path(loaded_meta_paths["project_dir"]["data"])

    