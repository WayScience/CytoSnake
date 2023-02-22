"""
Module: cytosnake_setup.py

This modules sets up the current directory as a project directory. Project
directories are dictated by the presence of `.cytosnake`, which contains
meta data for cytosnake to use.

"""
import json
import shutil

from pathlib import Path
from cytosnake.utils.cyto_paths import (
    get_cytosnake_package_path,
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
        Raised by pathlib.Path if `.cytosnake' already exists
    """

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
    proj_path = Path().absolute()

    # find software package directory path
    pkg_path = get_cytosnake_package_path()

    # get all important folders
    # - workflows: contains all workflow files (envs, rules)
    # - configs: contains yaml files providing workflow and cli configurations
    target_dirs = ["workflows", "configs"]
    for target_dir in target_dirs:

        # construct source directory path
        src_path = pkg_path / target_dir
        target_dst = proj_path / target_dir
        if not src_path.exists():
            raise FileNotFoundError(f"{src_path} does not exist")

        # move to dest: working project directory
        shutil.copytree(src_path, target_dst)


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
    cwd_path = Path().absolute()
    cwd_str = str(cwd_path)
    dir_paths = {"project_dir": get_project_dirpaths()}
    config_fpaths = {"config_dir": get_config_fpaths()}
    workflow_fpath = {"workflow_dir": get_workflow_fpaths()}
    proj_dir_path = {"project_dir_path": cwd_str}

    # merge all path dictionaries into one
    all_paths_dict = proj_dir_path | dir_paths | config_fpaths | workflow_fpath

    # write dict into yaml file
    # -- this steps generates the `_paths.yaml`
    save_path = cwd_path / ".cytosnake" / "_paths.yaml"
    with open(save_path, mode="w", encoding="utf-8") as stream:
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
