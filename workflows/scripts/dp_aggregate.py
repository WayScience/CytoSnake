from pathlib import Path

import snakemake
import yaml
from pycytominer.cyto_utils.DeepProfiler_processing import (
    AggregateDeepProfiler,
    DeepProfilerData,
)


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
        e_msg = "Unable to find aggregation configuration file"
        raise FileNotFoundError(e_msg)

    config_path = config_path_obj.absolute()
    with open(config_path, "r") as yaml_contents:
        loaded_configs = yaml.safe_load(yaml_contents)

    return loaded_configs


def aggregate_dp_profiles(
    dp_dir_path: str,
    index_file: str,
    outname: str,
    dp_data_config: str,
    aggregator_config: str,
) -> None:
    """_summary_

    Parameters
    ----------
    dp_dir_path : str
        path to directory that contains deep profiler data
    index_file : str
        path to metadata index file
    outname : str
        name of generated
    aggregator_config: str
        path pointing to deep profiler Aggregator class configurations
    agg_config : str
        path pointing deep profiler aggregation method configurations
    """

    # loading config files
    dp_data_configs = load_configs(dp_data_config)
    aggregator_configs = load_configs(aggregator_config)

    # loading and setting up input files
    dp_data_params = dp_data_configs["deep_profiler_data_configs"]["params"]
    dp_data = DeepProfilerData(
        profile_dir=dp_dir_path,
        index_file=index_file,
        file_extension=dp_data_params["file_extension"],
        filename_delimiter=dp_data_params["filename_delimiter"],
    )

    # placing single cell dp_data into Aggregator Class
    aggregator_params = aggregator_configs["AggregatorDeepProfiler_configs"]["params"]
    aggregator = AggregateDeepProfiler(
        deep_data=dp_data,
        aggregate_operation=aggregator_params["aggregate_operation"],
        aggregate_on=aggregator_params["aggregate_on"],
    )

    # -- aggregate singel cell and save it
    dp_agg_df = aggregator.aggregate_deep()
    dp_agg_df.to_csv(outname)


if __name__ == "__main__":

    # collecting snakemake inputs
    dp_path = str(snakemake.input["dp_features_dir"])
    index_file_path = str(snakemake.input["index_file"])
    out_name = str(snakemake.output)
    dp_data_config_path = str(snakemake.params["dp_data_configs"])
    aggregator_config_path = str(snakemake.params["aggregator_configs"])

    # executing function
    aggregate_dp_profiles(
        dp_dir_path=dp_path,
        index_file=index_file_path,
        dp_data_config=dp_data_config_path,
        aggregator_config=aggregator_config_path,
        outname=out_name,
    )
