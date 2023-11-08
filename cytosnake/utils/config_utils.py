import pathlib
from typing import Any

import yaml

from cytosnake.common.errors import WorkflowNotFoundError
from cytosnake.guards.path_guards import is_valid_path
from cytosnake.utils import cyto_paths


def load_configs(config_path: str | pathlib.Path) -> dict:
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
        config_path = pathlib.Path(config_path).resolve(strict=True)

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


def load_workflow_path(wf_name: str) -> pathlib.Path:
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
    return pathlib.Path(workflows[wf_name])


def load_data_path_configs():
    """Returns path pointing where the data folder is

    Returns
    -------
    Path
        Path to data folder in project directory
    """

    # load in `_paths.yaml` meta data
    loaded_meta_paths = load_meta_path_configs()
    return pathlib.Path(loaded_meta_paths["project_dir"]["data"])


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


def update_config(
    config_file_path: str | pathlib.Path,
    new_key: str,
    new_value: str | Any,
    update=False,
) -> None:
    """This updates config level in the upper level.

    Parameters
    ----------
    key : str
        key to add into the config
    value : str | Number
        Value to add to the given key.
    update : bool
        Update value to existing key

    Raises
    ------
    ValueError
        raised if the key value cannot be converted into a string
    TypeError
        raised if there is a type error with any parameter
    FileNotFoundError
        raised if the config file is not found
    KeyError
        raised if given key exists within the config file
    """

    # type checking
    if not isinstance(config_file_path, pathlib.Path):
        if isinstance(config_file_path, str):
            config_file_path = pathlib.Path(config_file_path).resolve(strict=True)
        raise TypeError(
            "'config_file_path' must be str or pathlib.Path which is able to resolve as a Path."
            f"not: {type(config_file_path)}"
        )
    if not isinstance(new_key, str):
        try:
            new_key = str(new_key)
        except Exception as exc:
            raise ValueError("Unable to convert key value into a string") from exc
    if not isinstance(new_value, (str, float, int)):
        raise TypeError("value must either be a string, float, or int")

    # open config file
    loaded_configs = load_configs(config_path=config_file_path)

    # updating configs
    new_key_value_pair = {new_key: new_value}
    if new_key in [keys for keys in loaded_configs.keys()]:
        if update:
            loaded_configs.update(new_key_value_pair)
        raise KeyError(
            f"key: `{new_key}` already exists in the config. "
            "Please use a different key to prevent key clashing"
        )
    loaded_configs.update(new_key_value_pair)

    # capture empty dictionaries before writing
    if not loaded_configs:
        raise RuntimeError(
            "Empty dictionary captured when entering new key value pairs"
        )

    # write our the updated config
    with open(config_file_path, mode="w", encoding="utf-8") as stream:
        yaml.dump(loaded_configs, stream)
