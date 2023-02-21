import yaml
from pathlib import Path

import pandas as pd
from pycytominer.operations import get_na_columns
from pycytominer import consensus


def build_dp_consensus(dp_profile: str, outname: str, config: str):
    """Builds a DeepProfiler consensus profile from aggregated DeepProfiler
    dataset

    Parameters
    ----------
    dp_profile : str
        Path to normalized aggregated dataset
    outname : str
        Name of generated consensus profile
    config : str
        Path pointing to configuration file
    """

    # loading config file
    consensus_path_obj = Path(config)
    if not consensus_path_obj.is_file():
        e_msg = "Unable to find consensus configuration file"
        raise FileNotFoundError(e_msg)

    consensus_config_path = consensus_path_obj.absolute()
    with open(consensus_config_path, "r") as yaml_contents:
        consensus_config = yaml.safe_load(yaml_contents)["consensus_config"]["params"]

    # drop na columns
    na_cols = get_na_columns(dp_profile, cutoff=0)
    concat_df = concat_df.drop(na_cols, axis="columns")

    # creating consensus
    dp_consensus_profile = consensus(
        concat_df,
        replicate_columns=consensus_config["replicate_columns"],
        operation=consensus_config["operation"],
        features=consensus_config["features"],
        compression_options=consensus_config["replicate_columns"],
        float_format=consensus_config["float_format"],
        modz_args=consensus_config["modz_args"],
    )

    dp_consensus_profile.to_csv(outname, sep="\t", index=False)


if __name__ == "__main__":

    # snakemake inputs
    norm_agg_dp_profile = str(snakemake.input)
    out_name = str(snakemake.output)
    config_path = str(snakemake.input)

    build_dp_consensus(
        dp_profile=norm_agg_dp_profile, outname=out_name, config=config_path
    )
