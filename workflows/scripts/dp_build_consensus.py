from pathlib import Path

import yaml
from pycytominer import consensus
from pycytominer.operations import get_na_columns


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


# building consensus profile from deep profiler features
if __name__ == "__main__":
    # snakemake inputs
    # more information how snakemake transfers workflow variables to scripts:
    # https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#python``
    norm_agg_dp_profile = str(snakemake.input)
    out_name = str(snakemake.output)
    config_path = str(snakemake.input)

    # building consensus profiles
    build_dp_consensus(
        dp_profile=norm_agg_dp_profile, outname=out_name, config=config_path
    )
