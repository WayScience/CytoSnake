import logging
from pathlib import Path
import yaml

import pandas as pd
from pycytominer import normalize


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
            "No DP features found, Are you sure that this dataframe is from DeepProfiler?"
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
        e_msg = "Unable to find Deep Profiler aggregation configuration file"
        raise FileNotFoundError(e_msg)

    config_path = config_path_obj.absolute()
    with open(config_path, "r") as yaml_contents:
        loaded_configs = yaml.safe_load(yaml_contents)

    return loaded_configs


def normalize_aggregate_dp_profiles(
    dp_aggregate_profile: str,
    outname: str,
    config: str,
):
    """Normalize aggregated DeepProfiler profiles. Generate
    Parameters
    ----------
    dp_aggregate_profile : str
        path to aggregate deep profiler features
    outname : str
        name of generated output
    config : str
        path to config file
    """

    # loading config file
    normalization_config = load_configs(config)
    normalization_params = normalization_config["normalize_configs"]["params"]

    # load in aggregated dataset
    dp_agg_df = pd.read_csv(dp_aggregate_profile)

    # normalizing aggregated dataset
    # -- this is to provide support for DeepProfiler aggregated datasets.
    # -- Currently pycytominer cannot infer DeepProfiler features
    if normalization_params["features"] == "infer":
        dp_features = infer_dp_features(dp_agg_df)
    else:
        dp_features = normalization_params["features"]

    # executing  normalization
    normalize(
        profiles=dp_agg_df,
        features=dp_features,
        image_features=normalization_params["features"],
        meta_features=normalization_params["meta_features"],
        samples=normalization_params["samples"],
        method=normalization_params["method"],
        compression_options=normalization_params["compression_options"],
        float_format=normalization_params["float_format"],
        mad_robustize_epsilon=normalization_params["mad_robustize_epsilon"],
        spherize_center=normalization_params["spherize_center"],
        spherize_method=normalization_params["spherize_method"],
        spherize_epsilon=normalization_params["spherize_epsilon"],
        output_file=outname,
    )


if __name__ == "__main__":

    # snakemake inputs
    agg_profile_path = str(snakemake.input)
    out_name = str(snakemake.output)
    config_path = str(snakemake.params["normalize_config"])

    # executing normalization function
    normalize_aggregate_dp_profiles(
        dp_aggregate_profile=agg_profile_path,
        outname=out_name,
        config=config_path,
    )
