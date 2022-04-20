import os
import yaml
from pathlib import Path
import pandas as pd
from pycytominer.cyto_utils.cells import SingleCells


def aggregate(sql_file, metadata_dir, barcode_platemap):
    """aggregates single cell data into aggregate profiles

    Paramters:
    ----------
    sql_file: str
            SQl file that contains single cell data obtain from a single plate
    metadata_dir : str
        associated metadata file with the single cell data
    barcode_platemap : str
        file containing the barcode id of each platedata

    Returns:
    --------
        No object is returned

        Generates cell count, aggregate aggregate, augmented aggregate profiles, and normalized
        augmented aggregate profiles in the results/ directory.
    """

    # loading paramters
    aggregate_path_obj = Path(snakemake.params["aggregate_config"])
    aggregate_config_path = aggregate_path_obj.absolute()
    with open(aggregate_config_path, "r") as yaml_contents:
        aggregate_configs = yaml.safe_load(yaml_contents)["single_cell_config"][
            "params"
        ]

    # Loading appropriate platemaps with given plate data
    plate = os.path.basename(sql_file).rsplit(".", 1)
    barcode_platemap_df = pd.read_csv(barcode_platemap)
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
        output_file=aggregate_configs["output_file"],
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
    cell_count_outfile = str(snakemake.output["cell_counts"])
    cell_count_df = single_cell_profile.count_cells()

    print("Saving cell counts in: {}".format(cell_count_outfile))
    cell_count_df = cell_count_df.merge(
        platemap_df, left_on="Image_Metadata_Well", right_on="well_position"
    ).drop(["WellRow", "WellCol", "well_position"], axis="columns")
    cell_count_df.to_csv(cell_count_outfile, sep="\t", index=False)

    print("aggregating cells")
    aggregate_outfile = str(snakemake.output["aggregate_profile"])
    single_cell_profile.aggregate_profiles(
        output_file=aggregate_outfile, compression_options="gzip"
    )


if __name__ == "__main__":

    # transforming snakemake objects into python standard datatypes
    sqlfiles = [str(sqlfile) for sqlfile in snakemake.input["sql_files"]]
    meta_data_dir = str(snakemake.input["metadata"])
    barcode = str(snakemake.input["barcodes"])

    # running the aggregate algorithm
    for sqlfile in sqlfiles:
        aggregate(sqlfile, meta_data_dir, barcode)
