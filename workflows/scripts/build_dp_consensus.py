from pathlib import Path
from typing import Union

import pandas as pd
import snakemake
import yaml
from pycytominer import consensus
from pycytominer.operations import get_na_columns


def infer_dp_features(dp_profile) -> list[str]:
    """Returns a list of Deep profiler features found within the single cell
    dataframe

    Parameters
    ----------
    dp_profile : pd.DataFrame
        dataframe features captured from deep profiler

    Returns
    -------
    list[str]
        list of deep profiler features

    Raises
    ------
    ValueError
        Raised if no Deep profiler features are found within the given DataFrame
    """
    dp_features = []
    metadata_model = dp_profile["Metadata_Model"]
    for column in dp_profile.columns.tolist():
        if any([column.startswith(f"{meta_model}_") for meta_model in metadata_model]):
            dp_features.append(column)

    if len(dp_features) <= 0:
        raise ValueError(
            "No DP features found, Are you sure that these are DeepProfiler features?"
        )

    return dp_features


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
        e_msg = "Unable to find consensus configuration file"
        raise FileNotFoundError(e_msg)

    config_path = config_path_obj.absolute()
    with open(config_path, "r") as yaml_contents:
        loaded_configs = yaml.safe_load(yaml_contents)

    return loaded_configs


def build_dp_consensus(
    profiles: Union[list[str], str], outname: str, config: str
) -> None:
    """Build consensus profile from single or multiple aggregate profiles

    Parameters
    ----------
    profiles : Union[list[str], str]
        string or list of paths that points to normalized aggregate dp profiles
    outname : str
        Name of generated consensus output profile
    config : str
        path that points to configurations
    """

    # loading in configurations
    consensus_config = load_configs(config)
    consensus_params = consensus_config["consensus_config"]["params"]

    # type checking of profiles
    # -- if it is a list of normalized aggregate dp profiles, concatenate
    if isinstance(profiles, list):
        norm_agg_df = pd.concat(
            [pd.read_csv(profile_path) for profile_path in profiles], sort=True
        ).rename(
            {
                "Image_Metadata_Plate": "Metadata_Plate",
                "Image_Metadata_Well": "Metadata_Well",
            },
            axis="columns",
        )

        # realignment of the meta data column names
        concat_metadata_cols = norm_agg_df.columns[
            norm_agg_df.columns.str.startswith("Metadata")
        ]
        concat_metadata_df = norm_agg_df.loc[:, concat_metadata_cols]

        # concatenating metadata
        norm_agg_df = norm_agg_df.drop(concat_metadata_cols, axis="columns")
        norm_agg_df = pd.concat([concat_metadata_df, norm_agg_df], axis="columns")

    # -- if it is a single path, create DataFrame
    else:
        norm_agg_df = pd.read_csv(profiles)

    # check if fetures are infer, if so, extract deep profiler features
    if consensus_params["features"] == "infer":
        dp_features = infer_dp_features(norm_agg_df)
    else:
        dp_features = consensus_params["features"]

    # dropping columns with NA values
    na_cols = get_na_columns(norm_agg_df, cutoff=0, features=dp_features)
    norm_agg_df = norm_agg_df.drop(na_cols, axis="columns")

    # generating consensus profile
    dp_consensus_profile = consensus(
        norm_agg_df,
        replicate_columns=consensus_params["replicate_columns"],
        operation=consensus_params["operation"],
        features=dp_features,
        compression_options=consensus_params["replicate_columns"],
        float_format=consensus_params["float_format"],
        modz_args=consensus_params["modz_args"],
    )

    # saving consensus profile
    dp_consensus_profile.to_csv(outname, sep="\t", index=False)


if __name__ == "__main__":

    # snakemake inputs
    # more information how snakemake transfers workflow variables to scripts:
    # https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#python``
    norm_agg_data = str(snakemake.input)
    out_name = str(snakemake.output)
    consensus_config_path = str(snakemake.params["consensus_config"])

    # executing consensus function
    build_dp_consensus(
        profiles=norm_agg_data, outname=out_name, config=consensus_config_path
    )
