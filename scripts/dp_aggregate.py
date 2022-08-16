import logging
from pathlib import Path

import yaml
from pycytominer.cyto_utils.DeepProfiler_processing import (
    DeepProfilerData,
    SingleCellDeepProfiler,
)


def aggregate_dp_profiles(
    dp_dir_path: str,
    index_file: str,
    outname: str,
    config: str,
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
    config : str
        path to config file
    """

    # loading config file
    dp_aggregate_config_path_obj = Path(config)
    if not dp_aggregate_config_path_obj.is_file():
        e_msg = "Unable to find consensus configuration file"
        logging.error(e_msg)
        raise FileNotFoundError(e_msg)

    dp_aggregate_config_path = dp_aggregate_config_path_obj.absolute()
    with open(dp_aggregate_config_path, "r") as yaml_contents:
        dp_aggregate_configs = yaml.safe_load(yaml_contents)
        logging.info("Consensus configuration loaded")

    # loading and setting up input files
    dp_data_configs = dp_aggregate_configs["deep_profiler_data_configs"]["params"]
    dp_data = DeepProfilerData(
        profile_dir=dp_dir_path,
        index_file=index_file,
        file_extension=dp_data_configs["file_extension"],
        filename_delimiter=dp_data_configs["filename_delimiter"],
    )

    # saving aggregate data
    dp_data.to_csv(outname)


if __name__ == "__mani__":

    # collecting snakemake inputs
    dp_path = str(snakemake.input["dp_features_dir"])
    index_file_path = str(snakemake.input["index_file"])
    out_name = str(snakemake.output)
    config_path = str(snakemake.params["consensus_configs"])

    # executing function
    aggregate_dp_profiles(
        dp_dir_path=dp_path,
        index_file=index_file_path,
        outname=out_name,
        config=config_path,
    )
