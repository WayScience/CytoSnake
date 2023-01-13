"""
Module: cytosnake_setup.py

This modele is responsible for setting up and recording file paths where analysis
are being conducted. 

It will be in charge of setting up path config files in order for cytosnake 

"""
import json
import shutil

from pathlib import Path
from cytosnake.utils.cyto_paths import (
    get_cytosnake_package_path,
    get_project_root,
    get_config_fpaths,
    get_workflow_fpaths,
    get_project_dirpaths,
)

from cytosnake.common.errors import display_error


def create_cytosnake_dir() -> None:
    """Sets up a `.cytosnake` folder that contains meta information that helps
    cytosnake's cli functionality. This will include configurations that
    contains states, paths and settings.


    Returns
    -------
    None
        Generates a `.cytosnake` file in current directory

    Raises:
    FileExistsError
        Raised by pathlib.Path if `.cytosnake' exists in cur
    ."""

    # create the `.cytosnake` folder in where `cytosnake init` mode was executed
    # - if the current director has a `.cytosnake` FileExistError will be raised
    cwd = Path().absolute()
    project_folder_path = cwd / ".cytosnake"
    try:
        project_folder_path.mkdir(exist_ok=False)
    except FileExistsError as e:
        msg = "'.cytosnake' directory exists in current directory"
        display_error(e, msg=msg)


# TODO: This should be in the `file_utils` module, may trigger circular imports
# -- this is a temporary solution
def transport_project_files() -> None:
    """obtains the necessary files from software package and transport them into
    current working directory"""

    # project directory
    proj_path = get_project_root()

    # find software package directory path
    pkg_path = get_cytosnake_package_path()

    # get all important folders
    # - workflows: contains all workflow files (envs, rules)
    # - configs: contains yaml files providing workflow and cli configurations
    target_dirs = ["workflows", "configs"]
    for target_dir in target_dirs:

        # construct source directory path
        src_path = pkg_path / target_dir
        if not src_path.exists():
            raise FileNotFoundError(f"{src_path} does not exist")

        # move to dest: working project directory
        shutil.copytree(src_path, proj_path)


def generate_meta_path_configs() -> None:
    """constructs a `_paths.yaml` file that contains the necessary path
    information for file handling. this allows `cytosnake` to know which folders
    to look for when performing tasks.

    returns
    -------
    none
        creates a `_paths.yaml` file in the `.cytosnake` project directory
    """

    # get all pathing info
    dir_paths = {"project_dir": get_project_dirpaths()}
    config_fpaths = {"config_dir": get_config_fpaths()}
    workflow_fpath = {"workflow_dir": get_workflow_fpaths()}

    # merge all path dictionaries into one
    all_paths_dict = dir_paths | config_fpaths | workflow_fpath

    # write dict into yaml file
    save_path = get_project_root() / ".cytosnake" / "_paths.yaml"
    with open(save_path, "w") as stream:
        json.dump(all_paths_dict, stream, indent=4)


def setup_cytosnake_env() -> None:
    """main wrapper function that sets up current directory into a `cytosnake`
    project directory. this means that all the analysis being conducted within
    the current directory will allow
    """

    # create '.cytosnake' file
    create_cytosnake_dir()

    # move necessary files to directory
    transport_project_files()

    # create `_paths.yaml` path meta data in `.cytosnake` project dir
    generate_meta_path_configs()
