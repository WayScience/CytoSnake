from pathlib import Path

import yaml

from cytosnake.common.errors import WorkflowNotFoundError
from cytosnake.guards.path_guards import is_valid_path
from cytosnake.utils import cyto_paths


def load_configs(config_path: str | Path) -> dict:
    """Returns a dictionary of given configurations

    Parameters
    ----------
    config_path : str | Path
        path to config file

    Returns
    -------
    dict
        configuration dictionary

    Raises
    ------
    FileNotFoundError
        raised if provided config file paths is invalid
    """

    # check if config path is a valid path
    if not is_valid_path(config_path):
        raise FileNotFoundError(f"Invalid config path provided: {config_path}")
    if isinstance(config_path, str):
        config_path = Path(config_path).resolve(strict=True)

    # loading in config_path
    with open(config_path, "r") as yaml_contents:
        loaded_configs = yaml.safe_load(yaml_contents)

    return loaded_configs


def load_general_configs() -> dict:
    """Loads cytosnake's general configurations

    Returns:
    -------
    dict
        dictionary containing the cytosnake general configs
    """
    config_dir_path = cyto_paths.get_config_dir_path() / "configuration.yaml"
    return load_configs(config_dir_path)


def load_meta_path_configs() -> dict:
    """Loads the metadata path from `.cytosnake/_paths.yaml` file

    Returns
    -------
    dict
        meta path contents from the `_paths.yaml` file
    """

    # construct path to `.cytosnake/_paths.yaml` file
    meta_path = cyto_paths.get_meta_path() / "_paths.yaml"
    return load_configs(meta_path)


def load_workflow_path(wf_name: str) -> Path:
    """Loads in configurations and returns path pointing to
    workflow

    Parameters
    ----------
    wf_name : str
        workflow name

    Returns
    -------
    pathlib.PosixPath
        Path to workflow

    Raises
    ------
    WorkFlowNotFound
        Raised if the desired workflow is not found.
    """
    # load configurations from `.cytosnake`
    meta_path = cyto_paths.get_meta_path() / "_paths.yaml"

    # loading loading workflow paths
    general_config = load_configs(meta_path)
    workflows = general_config["workflow_dir"]["workflow"]

    # checking if the workflow exists
    if wf_name not in workflows.keys():
        raise WorkflowNotFoundError(f"Unable to find {wf_name} workflow")

    # returning workflow path
    return Path(workflows[wf_name])


def load_data_path_configs():
    """Returns path pointing where the data folder is

    Returns
    -------
    Path
        Path to data folder in project directory
    """

    # load in `_paths.yaml` meta data
    loaded_meta_paths = load_meta_path_configs()
    return Path(loaded_meta_paths["project_dir"]["data"])


def load_workflow_paths_config() -> dict:
    # load in _path.yaml and select key where all workflow paths are
    loaded_meta_paths = load_meta_path_configs()
    return loaded_meta_paths["workflow_dir"]["workflow"]


def load_cytosnake_configs() -> dict:
    """Loads in CytoSnake's general configuration

    Returns
    -------
    dict
        CytoSnake configs
    """

    # gets absolute path to cytosnake configs and load the configs
    cytosnake_config_path = cyto_paths.get_cytosnake_config_path()
    return load_configs(cytosnake_config_path)
