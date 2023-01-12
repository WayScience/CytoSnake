"""
module: cyto_paths.py

This module will contain functions that handles `cytosnake's` pathing
"""

from pathlib import Path
import cytosnake


def get_cytosnake_package_path() -> Path:
    """Returns paths where the package is installed

    Return
    ------
    Path
        Returns absolute path to
    """

    # get location of this file
    # -- check if he ".git" folder exists
    project_path = Path(cytosnake.__file__)
    git_path = project_path / ".git"
    if not git_path.exists:
        raise FileNotFoundError("Unable to find cytosnake package path")

    return project_path


def get_project_root() -> Path:
    """Returns complete path where the analysis is taken place. The function
    will check if `.cytosnake` folder exists, if not an error will be raised.
    """

    # get current working directory
    cwd = Path().absolute()

    # check if the `.cytosnake` folder exist
    project_folder = cwd / ".cytosnake"
    if not project_folder.exist():
        raise FileNotFoundError("Current directory is not ")
