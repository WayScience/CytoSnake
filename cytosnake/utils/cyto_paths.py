"""
module: cyto_paths.py

This module will contain functions that handles `cytosnake's` pathing
"""
import pathlib
from typing import Optional, TypeVar

# cytosnake imports
import cytosnake
from cytosnake.guards.path_guards import is_valid_path
from cytosnake.utils.file_utils import file_search, find_project_dir

# user defined types without importing module
Namespace = TypeVar("Namespace")


def get_meta_path() -> pathlib.Path:
    """returns meta path configurational file path.

    Returns
    -------
    Path
        Path object pointing to `_paths.yaml` config file
    """

    # getting project directory
    proj_dir = find_project_dir() / ".cytosnake"
    if proj_dir is None:
        raise FileNotFoundError("Current directory is not a project folder")

    return proj_dir.absolute()


# TODO: this should be placed in the guards module
def is_cytosnake_dir(dir_path: Optional[str | pathlib.Path] = None) -> bool:
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
        dir_path = pathlib.Path().absolute()

    # type checking
    if not is_valid_path(dir_path):
        _type = type(dir_path)
        raise ValueError(f"dir_path must be either str or Path type, not: {_type}")
    if isinstance(dir_path, str):
        dir_path = pathlib.Path(str).absolute()

    # get project older
    cyto_proj_dir = dir_path / ".cytosnake"
    return bool(cyto_proj_dir.exists())


def get_cytosnake_package_path() -> pathlib.Path:
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
    project_path = pathlib.Path(cytosnake.__path__[0]).parent
    git_path = project_path / ".git"
    if not git_path.exists:
        raise FileNotFoundError("Unable to find cytosnake package path")

    return project_path


def get_project_root() -> pathlib.Path:
    """Returns complete path where cytosnake performs the analysis. The function
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
    project_dir = find_project_dir()
    if project_dir is None:
        raise NotADirectoryError("Unable to find project directory")

    return project_dir


def get_results_dir_path() -> pathlib.Path:
    """Returns path to results directory. This is where all the outputs generated
    from the workflow will be placed.

    Returns
    -------
    pathlib.Path
        absolute path to results directory
    """

    # building results path
    # will not be check if because the results directory is not generated
    # until the workflow has been executed
    return get_project_root() / "results"


def get_workflow_fpaths() -> dict:
    """Obtains all file paths located in the `workflows` folder as a dictionary.

    Returns
    -------
    dict
        Structured dictionary containing directory names and paths are key value
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


def get_config_dir_path() -> pathlib.Path:
    """Returns path to configuration folder

    Returns
    -------
    Path
        Path to config directory
    """

    proj_root_path = get_project_root()
    config_path = proj_root_path / "configs"
    if not is_valid_path(config_path):
        raise FileNotFoundError("Unable to find config directory")

    return config_path


def get_cytosnake_config_path() -> pathlib.Path:
    """Returns absolute path to CytoSnake's general config file .

    Returns
    -------
    pathlib.Path
        path to cytosnake general config file
    """
    config_dir_path = get_config_dir_path() / "configuration.yaml"
    if not is_valid_path(config_dir_path):
        raise FileNotFoundError("Unable to find CytoSnake general config file")

    return config_dir_path


def get_config_fpaths() -> dict:
    """Obtains all file paths located in the `configs` folder as a dictionary.

    Returns
    -------
    dict
        structured dictionary with directory name and file paths as key value pairs
    """
    return file_search(get_config_dir_path())


def get_project_dirpaths(args: Namespace) -> dict:
    """returns a dictionary containing directory name and path as key value
    pairs.

    Parameters
    ----------
    args : Namespace
        Uses argparse's  Namespace object to add additional information into the
        data section in `_path.yaml`

    Returns
    -------
    dict
        directory name and path as key value pairs.
    """

    # getting project root folder
    proj_root_path = get_project_root()

    # creating dict containing {dir_name : dir_path}
    # -- ignore unwanted directories
    ignore_dirs = [".vscode", ".git", "CytoSnake.egg-info"]

    # -- collect dir names and paths
    all_dirs = {}
    for _file in proj_root_path.glob("*"):
        if _file.name in ignore_dirs:
            continue

        # if the file is a directory, get paths of individual files
        # -- this is done recursively.
        if _file.is_dir():
            abs_path = str(_file.absolute())
            all_dirs[_file.name] = abs_path

            # use the namespace arguments to write _paths.yaml
            if _file.name.lower() == "data":
                abs_data_path = pathlib.Path(abs_path)

                # creating  dictionary for the plate data
                data_dir_conts = {}
                for platename in args.data:
                    name = platename.split(".")[0]
                    platedata_path = abs_data_path / platename
                    data_dir_conts[name] = str(platedata_path)

                # adding single files paths
                data_dir_conts["data"] = abs_path
                if args.barcode is None:
                    data_dir_conts["barcode"] = None
                else:
                    data_dir_conts["barcode"] = str(abs_data_path / args.barcode)
                data_dir_conts["metadata"] = str(abs_data_path / args.metadata)
                all_dirs["data_directory_contents"] = data_dir_conts

    return all_dirs
