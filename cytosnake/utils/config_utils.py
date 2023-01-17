from pathlib import Path 
import yaml

from cytosnake.utils.cyto_paths import get_meta_path
from cytosnake.guards.path_guards import is_valid_path
from cytosnake.common.errors import WorkflowNotFoundError


def load_configs(config_path: str | Path) -> dict:
    """Returns a dictionary of given configurations

    Parameters
    ----------
    config_path : Union[str, Path]
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
        raise FileNotFoundError("Invalid config path provided")
    if isinstance(config_path, str):
        config_path = Path(config_path).absolute()
    if not config_path.is_absolute():
        config_path = config_path.absolute()

    # loading in config_path
    with open(config_path, "r") as yaml_contents:
        loaded_configs = yaml.safe_load(yaml_contents)

    return loaded_configs


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
    meta_path = get_meta_path() / "_paths.yaml"

    # loading loading workflow paths
    general_config = load_configs(meta_path)
    workflows = general_config["workflow_dir"]["workflow"]

    # checking if the workflow exists
    if wf_name not in workflows.keys():
        raise WorkflowNotFoundError(f"Unable to find {wf_name} workflow")

    # returning workflow path
    return Path(workflows[wf_name])
