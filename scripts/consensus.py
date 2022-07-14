import os
import yaml
from pathlib import Path

import pandas as pd
from pycytominer.operations import get_na_columns
from pycytominer import aggregate


def build_consensus(profile_list: list, consensus_file_out: str, config: str) -> None:
    """Concatenates all normalized aggregated features into one
    pandas DataFrame

    Parameters
    ----------
    profiles : list
        list of paths pointing to normalized aggregated features

    consensus_file_out : str
        path where to write consensus signature file

    Returns
    -------
    pd.DataFrame
        concatenated normalized aggregated features
    """
    # loading config file
    aggregate_path_obj = Path(config)
    aggregate_config_path = aggregate_path_obj.absolute()
    with open(aggregate_config_path, "r") as yaml_contents:
        consensus_config = yaml.safe_load(yaml_contents)["consensus_config"]["params"]

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

    concat_df = concat_df.drop(concat_metadata_cols, axis="columns")
    concat_df = pd.concat([concat_metadata_df, concat_df], axis="columns")

    # dropping columns with na values
    na_cols = get_na_columns(concat_df, cutoff=0)
    concat_df = concat_df.drop(na_cols, axis="columns")

    # aggregate
    consensus_aggregate_config = consensus_config["aggregate"]
    x_median_df = aggregate(
        concat_df,
        strata=consensus_aggregate_config["strata"],
        features=consensus_aggregate_config["features"],
        operation=consensus_aggregate_config["operation"],
        compute_object_count=consensus_aggregate_config["compute_object_count"],
        object_feature=consensus_aggregate_config["object_feature"],
        subset_data_df=consensus_aggregate_config["subset_data_df"],
        compression_options=consensus_aggregate_config["compression_options"],
        float_format=consensus_aggregate_config["float_format"],
    )

    # Output Profile Mapping for Downstream Analysis
    profile_id_mapping_df = x_median_df.loc[
        :, x_median_df.columns.str.startswith("Metadata")
    ]
    file = os.path.join("data", "profile_id_metadata_mapping.tsv")
    profile_id_mapping_df.to_csv(file, sep="\t", index=False)

    x_median_df.to_csv(consensus_file_out, sep="\t", index=False)


if __name__ in "__main__":

    # loading inputs
    inputs = [str(f_in) for f_in in snakemake.input]
    output = str(snakemake.output)
    config_path = str(snakemake.params["consensus_configs"])

    # concatenated all Normalized aggregated profiles
    build_consensus(profile_list=inputs, consensus_file_out=output, config=config_path)
