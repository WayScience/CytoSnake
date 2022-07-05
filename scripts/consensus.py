import os
import pandas as pd
from pycytominer.operations import get_na_columns
from pycytominer import aggregate


def build_consensus(profile_list: list, consensus_file_out) -> None:
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

    # dropping Costes Features
    costes_cols_to_drop = [x for x in concat_df.columns if "costes" in x.lower()]
    concat_df = concat_df.drop(costes_cols_to_drop, axis="columns")

    x_groupby_cols = ["Metadata_gene_name", "Metadata_pert_name", "Metadata_cell_line"]

    x_metacount_df = (
        concat_df.loc[:, x_groupby_cols]
        .assign(n_measurements=1)
        .groupby(x_groupby_cols)
        .count()
        .reset_index()
        .assign(data_type="cell_painting")
        .merge(
            concat_df.loc[:, x_groupby_cols + ["Metadata_Well", "Metadata_Plate"]],
            how="left",
            on=x_groupby_cols,
        )
    )

    # TODO: this path should be in the meta data
    y_df = pd.read_csv("./data/normalized_cell_health_labels.tsv", sep="\t").drop(
        ["plate_name", "well_col", "well_row"], axis="columns"
    )

    y_groupby_cols = ["guide", "cell_id"]
    y_metacount_df = (
        y_df.loc[:, y_groupby_cols]
        .assign(n_measurements=1)
        .groupby(y_groupby_cols)
        .count()
        .reset_index()
        .assign(data_type="cell_health")
    )

    x_groupby_cols = ["Metadata_gene_name", "Metadata_pert_name", "Metadata_cell_line"]
    x_metacount_df = (
        concat_df.loc[:, x_groupby_cols]
        .assign(n_measurements=1)
        .groupby(x_groupby_cols)
        .count()
        .reset_index()
        .assign(data_type="cell_painting")
        .merge(
            concat_df.loc[:, x_groupby_cols + ["Metadata_Well", "Metadata_Plate"]],
            how="left",
            on=x_groupby_cols,
        )
    )

    y_groupby_cols = ["guide", "cell_id"]
    y_metacount_df = (
        y_df.loc[:, y_groupby_cols]
        .assign(n_measurements=1)
        .groupby(y_groupby_cols)
        .count()
        .reset_index()
        .assign(data_type="cell_health")
    )

    all_measurements_df = (
        x_metacount_df.merge(
            y_metacount_df,
            left_on=["Metadata_pert_name", "Metadata_cell_line"],
            right_on=["guide", "cell_id"],
            suffixes=["_paint", "_health"],
            how="inner",
        )
        .sort_values(by=["Metadata_cell_line", "Metadata_pert_name"])
        .reset_index(drop=True)
        .drop(["Metadata_Well", "guide", "cell_id"], axis="columns")
    )

    # this should also be in the paramters as well
    file = os.path.join("results", "all_profile_metadata.tsv")
    all_measurements_df.to_csv(file, sep="\t", index=False)

    x_median_df = aggregate(
        concat_df,
        strata=["Metadata_cell_line", "Metadata_pert_name"],
        features="infer",
        operation="median",
    )

    x_median_df = (
        x_median_df.query(
            "Metadata_pert_name in @all_measurements_df.Metadata_pert_name.unique()"
        )
        .query("Metadata_cell_line in @all_measurements_df.Metadata_cell_line.unique()")
        .reset_index(drop=True)
        .reset_index()
        .rename({"index": "Metadata_profile_id"}, axis="columns")
    )
    x_median_df.Metadata_profile_id = [
        "profile_{}".format(x) for x in x_median_df.Metadata_profile_id
    ]

    # Output Profile Mapping for Downstream Analysis
    profile_id_mapping_df = x_median_df.loc[
        :, x_median_df.columns.str.startswith("Metadata")
    ]
    file = os.path.join("data", "profile_id_metadata_mapping.tsv")
    profile_id_mapping_df.to_csv(file, sep="\t", index=False)

    cell_health_meta_features = ["cell_id", "guide"]
    cell_health_features = y_df.drop(
        cell_health_meta_features, axis="columns"
    ).columns.tolist()
    y_meta_merge_cols = [
        "Metadata_profile_id",
        "Metadata_pert_name",
        "Metadata_cell_line",
    ]
    y_median_df = aggregate(
        y_df,
        strata=cell_health_meta_features,
        features=cell_health_features,
        operation="median",
    )

    y_median_df = y_median_df.reset_index(drop=True).merge(
        x_median_df.loc[:, y_meta_merge_cols],
        left_on=["guide", "cell_id"],
        right_on=["Metadata_pert_name", "Metadata_cell_line"],
        how="right",
    )

    # Get columns in correct order
    y_columns = (
        y_meta_merge_cols
        + y_median_df.loc[
            :, ~y_median_df.columns.str.startswith("Metadata_")
        ].columns.tolist()
    )

    y_median_df = y_median_df.loc[:, y_columns].drop(
        ["guide", "cell_id"], axis="columns"
    )

    # Confirm that matrices are aligned
    pd.testing.assert_series_equal(
        x_median_df.Metadata_profile_id,
        y_median_df.Metadata_profile_id,
        check_names=True,
    )

    # Are the guides aligned?
    pd.testing.assert_series_equal(
        x_median_df.Metadata_pert_name, y_median_df.Metadata_pert_name, check_names=True
    )

    # Are the cells aligned?
    pd.testing.assert_series_equal(
        x_median_df.Metadata_cell_line, y_median_df.Metadata_cell_line, check_names=True
    )

    x_median_df.to_csv(consensus_file_out, sep="\t", index=False)


if __name__ in "__main__":

    inputs = [str(f_in) for f_in in snakemake.input]
    output = str(snakemake.output)

    # concatenated all Normalized aggregated profiles
    concat_dataset = build_consensus(profile_list=inputs, consensus_file_out=output)
    # concat_dataset.to_csv(output, compression="gzip")
