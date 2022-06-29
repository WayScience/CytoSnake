import os
import yaml
import argparse
from pathlib import Path

import pandas as pd
from pycytominer.cyto_utils.cells import SingleCells


def aggregate(
    sql_file: str,
    metadata_dir: str,
    barcode_path: str,
    cell_count_out: str,
    aggregate_file_out: str,
    config: str,
):
    """aggregates single cell data into aggregate profiles

    Parameters:
    ----------
    sql_file: str
            SQL file that contains single cell data obtain from a single plate
    metadata_dir : str
        associated metadata file with the single cell data
    barcode_path : str
        file containing the barcode id of each plate data
    aggregate_file_out : str
        output file generated for aggregate profiles
    cell_count_out: str
        output file generating cell counts
    config: str
        Path to config file

    Returns:
    --------
        No object is returned
        Generates cell counts and aggregate profile output
    """

    # loading parameters
    aggregate_path_obj = Path(config)
    aggregate_config_path = aggregate_path_obj.absolute()
    with open(aggregate_config_path, "r") as yaml_contents:
        aggregate_configs = yaml.safe_load(yaml_contents)["single_cell_config"][
            "params"
        ]

    # Loading appropriate platemaps with given plate data
    plate = os.path.basename(sql_file).rsplit(".", 1)
    barcode_platemap_df = pd.read_csv(barcode_path)
    platemap = barcode_platemap_df.query(
        "Assay_Plate_Barcode == @plate"
    ).Plate_Map_Name.values[0]
    platemap_file = os.path.join(metadata_dir, "platemap", "{}.csv".format(platemap))
    platemap_df = pd.read_csv(platemap_file)

    sqlite_file = "sqlite:///{}".format(sql_file)
    # strata = ["Image_Metadata_Plate", "Image_Metadata_Well"]
    # image_cols = ["TableNumber", "ImageNumber"]
    single_cell_profile = SingleCells(
        sqlite_file,
        strata=aggregate_configs["strata"],
        image_cols=aggregate_configs["image_cols"],
        aggregation_operation=aggregate_configs["aggregation_operation"],
        output_file=aggregate_file_out,
        merge_cols=aggregate_configs["merge_cols"],
        add_image_features=aggregate_configs["add_image_features"],
        image_feature_categories=aggregate_configs["image_feature_categories"],
        features=aggregate_configs["features"],
        load_image_data=aggregate_configs["load_image_data"],
        subsample_frac=aggregate_configs["subsample_frac"],
        subsampling_random_state=aggregate_configs["subsampling_random_state"],
        fields_of_view=aggregate_configs["fields_of_view"],
        fields_of_view_feature="Image_Metadata_Well",
        object_feature=aggregate_configs["object_feature"],
    )

    # counting cells in each well and saving it as csv file
    print("Counting cells within each well")
    cell_count_df = single_cell_profile.count_cells()

    print("Saving cell counts in: {}".format(cell_count_out))
    cell_count_df = cell_count_df.merge(
        platemap_df, left_on="Image_Metadata_Well", right_on="well_position"
    ).drop(["WellRow", "WellCol", "well_position"], axis="columns")
    cell_count_df.to_csv(cell_count_out, sep="\t", index=False)

    print("aggregating cells")
    single_cell_profile.aggregate_profiles(
        output_file=aggregate_file_out, compression_options="gzip"
    )


if __name__ == "__main__":

    # CLI arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", type=str, required=True, help="cell plate feature data"
    )
    parser.add_argument(
        "-o", "--output", type=str, required=True, help="aggregated output files"
    )
    parser.add_argument(
        "-co", "--cellcount_out", type=str, required=True, help="cell out output file"
    )
    parser.add_argument(
        "-m",
        "--metadata_dir",
        type=str,
        required=True,
        help="Path to metadata directory",
    )
    parser.add_argument(
        "-b", "--barcode", type=str, required=True, help="Path to barcode labels"
    )
    parser.add_argument(
        "-c", "--config", type=str, required=True, help="Path to config file"
    )
    args = parser.parse_args()

    aggregate(
        sql_file=args.input,
        metadata_dir=args.metadata_dir,
        barcode_path=args.barcode,
        aggregate_file_out=args.output,
        cell_count_out=args.cellcount_out,
        config=args.config,
    )
