from calendar import c
from pathlib import Path
from re import A
import numpy as np
import pandas as pd

from pycytominer.consensus import modz
from pycytominer import get_na_columns, aggregate


def concatenate_data(profiles: list) -> pd.DataFrame:
    """Concatenates all normalized aggregated features into one
    pandas DataFrame

    Parameters
    ----------
    profiles : list
        list of paths pointing to normalized aggregated features

    Returns
    -------
    pd.DataFrame
        concatenated normalized aggregated features
    """
    concat_df = (
        pd.concat(profiles, sort=True)
        .rename(
            {
                "Image_Metadata_Plate": "Metadata_Plate",
                "Image_Metadata_Well": "Metadata_Well",
            },
            axis="column",
        )
        .drop(["Metadata_broad_sample"], axis="columns")
    )

    # realignment of the meta data column names
    concat_metadata_cols = concat_df.columns[
        concat_df.columns.str.startswith("Metadata")
    ]
    concat_metadata_df = concat_df.loc[:, concat_metadata_cols]
    concat_df = concat_df.drop(concat_metadata_cols, axis="columns")
    concat_df = pd.concat([concat_metadata_df, concat_df])

    # dropping columns with na values
    na_cols = get_na_columns(concat_df, cutoff=0)
    concat_df = concat_df.drop(na_cols, axis="columns")

    # droping costes features
    costes_cols = [x for x in concat_df.columns if "costes" in x.lower()]
    concat_df = concat_df.drop(costes_cols, axis="columns")

    return cocnat_df


if __name__ in "__main__":

    inputs = [str(f_in) for f_in in snakemake.input]
    output = None

    # concatenated all Normalized aggregated profiles
    concat_dataset = concatenate_data(inputs)
