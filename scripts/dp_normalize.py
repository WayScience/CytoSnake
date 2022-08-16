import logging
from pathlib import Path
import yaml

import pandas as pd
from pycytominer.cyto_utils.DeepProfiler_processing import AggregateDeepProfiler
from pycytominer import normalize
from cytopipe.utils.feature_utils import infer_dp_features


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
    # -- loading normalize config
    normalization_path_obj = Path(config)
    if not normalization_path_obj.is_file():
        e_msg = "Unable to find consensus configuration file"
        logging.error(e_msg)
        raise FileNotFoundError(e_msg)

    dp_normalization_config_path = normalization_path_obj.absolute()
    with open(dp_normalization_config_path, "r") as yaml_contents:
        normalization_config = yaml.safe_load(yaml_contents)["normalize_configs"][
            "params"
        ]
        logging.info("Normalization configuration loaded")

    # load in aggregated dataset
    dp_agg_df = pd.read_csv(dp_aggregate_profile)

    # normalizing aggregated dataset
    # -- this is to provide support for DeepProfiler aggregated datasets.
    # -- Currently pycytominer cannot infer DeepProfiler features
    if normalization_config["features"] == "infer":
        dp_features = infer_dp_features(dp_agg_df)
    else:
        dp_features = normalization_config["features"]

    # executing  normalization
    normalize(
        profiles=dp_agg_df,
        features=dp_features,
        image_features=normalization_config["features"],
        meta_features=normalization_config["meta_features"],
        samples=normalization_config["samples"],
        method=normalization_config["method"],
        compression_options=normalization_config["compression_options"],
        float_format=normalization_config["float_format"],
        mad_robustize_epsilon=normalization_config["mad_robustize_epsilon"],
        spherize_center=normalization_config["spherize_center"],
        spherize_method=normalization_config["spherize_method"],
        spherize_epsilon=normalization_config["spherize_epsilon"],
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
