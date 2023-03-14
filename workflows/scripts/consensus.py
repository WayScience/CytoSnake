import logging
from pathlib import Path

import pandas as pd
import yaml
from pycytominer import consensus
from pycytominer.operations import get_na_columns


def build_consensus(
    profile_list: list, consensus_file_out: str, config: str, log_file: str
) -> None:
    """Concatenates all normalized aggregated features into one
    pandas DataFrame

    Parameters
    ----------
    profile_list : list
        list of paths pointing to normalized aggregated features
    consensus_file_out : str
        path where to write consensus signature file
    config : str
        Path pointing to configuration file
    log_file : str
        Path pointing to log file

    Returns
    -------
    pd.DataFrame
        concatenated normalized aggregated features
    """

    # initiating Logger
    log_path = Path(log_file).absolute()
    logging.basicConfig(
        filename=log_path,
        encoding="utf-8",
        level=logging.DEBUG,
        format="%(asctime)s.%(msecs)03d - %(levelname)s - %(thread)d - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logging.info("Building consensus profile")

    # loading config file
    consensus_path_obj = Path(config)
    if not consensus_path_obj.is_file():
        e_msg = "Unable to find consensus configuration file"
        logging.error(e_msg)
        raise FileNotFoundError(e_msg)

    consensus_config_path = consensus_path_obj.absolute()
    with open(consensus_config_path, "r") as yaml_contents:
        consensus_config = yaml.safe_load(yaml_contents)["consensus_config"]["params"]
        logging.info("Consensus configuration loaded")

    concat_df = pd.concat(
        [pd.read_csv(profile_path) for profile_path in profile_list], sort=True
    ).rename(
        {
            "Image_Metadata_Plate": "Metadata_Plate",
            "Image_Metadata_Well": "Metadata_Well",
        },
        axis="columns",
    )

    # realignment of the meta data column names
    concat_metadata_cols = concat_df.columns[
        concat_df.columns.str.startswith("Metadata")
    ]
    concat_metadata_df = concat_df.loc[:, concat_metadata_cols]

    # concatenating metadata
    logging.info("Concatenating metadata")
    concat_df = concat_df.drop(concat_metadata_cols, axis="columns")
    concat_df = pd.concat([concat_metadata_df, concat_df], axis="columns")

    # dropping columns with na values
    logging.info("Dropping columns with 'na' values")
    na_cols = get_na_columns(concat_df, cutoff=0)
    concat_df = concat_df.drop(na_cols, axis="columns")

    # generating consensus profile
    logging.info("Generating consensus profile")
    x_consensus_df = consensus(
        concat_df,
        replicate_columns=consensus_config["replicate_columns"],
        operation=consensus_config["operation"],
        features=consensus_config["features"],
        compression_options=consensus_config["replicate_columns"],
        float_format=consensus_config["float_format"],
        modz_args=consensus_config["modz_args"],
    )

    # Saving consensus profile
    logging.info(f"Saving consensus profile: {consensus_file_out}")
    x_consensus_df.to_csv(consensus_file_out, sep="\t", index=False)


if __name__ in "__main__":

    # loading inputs
    inputs = [str(f_in) for f_in in snakemake.input]
    output = str(snakemake.output)
    config_path = str(snakemake.params["consensus_configs"])
    log_path = str(snakemake.log)

    # concatenated all Normalized aggregated profiles
    build_consensus(
        profile_list=inputs,
        consensus_file_out=output,
        config=config_path,
        log_file=log_path,
    )
