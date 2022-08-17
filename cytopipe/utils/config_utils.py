from pathlib import Path
import yaml

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