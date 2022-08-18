from pathlib import Path
import yaml

from ..common.errors import WorkflowNotFoundError


def load_configs(config: str) -> dict:
    """Returns a dictionary of given configurations

    Parameters
    ----------
    config : str
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
    config_path_obj = Path(config)
    if not config_path_obj.is_file():
        e_msg = "Unable to find Deep Profiler aggregation configuration file"
        raise FileNotFoundError(e_msg)

    config_path = config_path_obj.absolute()
    with open(config_path, "r") as yaml_contents:
        loaded_configs = yaml.safe_load(yaml_contents)

    return loaded_configs


def load_workflow_path(wf_name: str) -> str:
    """Loads in configurations and returns path pointing to
    workflow

    Parameters
    ----------
    wf_name : str
        workflow name

    Returns
    -------
    str
        Path to workflow

    Raises
    ------
    WorkFlowNotFound
        Raised if the desired workflow is not found.
    """
    # loading loading workflow paths
    general_config = load_configs("../../configs/configuration.yaml")
    workflows = general_config["workflow_paths"]

    # checking if the workflow exists
    if wf_name not in workflows.keys():
        raise WorkflowNotFoundError(f"Unable to find {wf_name} workflow")

    # returning workflow path
    return workflows[wf_name]


