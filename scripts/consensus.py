import pandas as pd

from pycytominer.operations import get_na_columns


def concatenate_data(profile_list: list) -> pd.DataFrame:
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
    concat_df = pd.concat([concat_metadata_df, concat_df])

    # dropping columns with na values
    na_cols = get_na_columns(concat_df, cutoff=0)
    concat_df = concat_df.drop(na_cols, axis="columns")

    return concat_df


if __name__ in "__main__":

    inputs = [str(f_in) for f_in in snakemake.input]
    output = str(snakemake.output)

    # concatenated all Normalized aggregated profiles
    concat_dataset = concatenate_data(inputs)
    concat_dataset.to_csv(output, compression="gzip")
