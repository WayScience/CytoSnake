import logging
import os
from pathlib import Path

import pandas as pd
from snakemake.script import Snakemake
import yaml
from pycytominer.cyto_utils.cells import SingleCells


def aggregate(
    sql_file: str,
    metadata_dir: str,
    barcode_path: str,
    cell_count_out: str,
    aggregate_file_out: str,
    config: str,
    log_file: str,
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
    log_file : str
        Path to log file

    Returns:
    --------
        No object is returned
        Generates cell counts and aggregate profile output
    """
    # opening logger
    log_path = Path(log_file).absolute()
    logging.basicConfig(
        filename=log_path,
        encoding="utf-8",
        level=logging.DEBUG,
        format="%(asctime)s.%(msecs)03d - %(levelname)s - %(thread)d - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    logging.info(f"Aggregating Single Cells in: {sql_file}")

    # loading parameters
    logging.info(f"Loading single cell aggregation configurations from: {config}")
    aggregate_path_obj = Path(config)
    aggregate_config_path = aggregate_path_obj.absolute()
    config_check = Path(aggregate_config_path).is_file()

    if not config_check:
        e_msg = "Unable to find aggregation configuration file"
        logging.error(e_msg)
        raise FileNotFoundError(e_msg)

    with open(aggregate_config_path, "r") as yaml_contents:
        aggregate_configs = yaml.safe_load(yaml_contents)["single_cell_config"][
            "params"
        ]
        logging.info("Aggregation configurations loaded.")

    # Loading appropriate platemaps with given plate data
    logging.info(f"Loading plate data from: {sql_file}")
    plate = os.path.basename(sql_file).rsplit(".", 1)
    plate_file_check = Path(sql_file).is_file()

    if not plate_file_check:
        e_msg = f"Unable to find plate file: {sql_file}"
        logging.error(e_msg)
        raise FileNotFoundError(e_msg)

    # generating cell count
    logging.info(f"Loading barcodes from: {barcode_path}")
    barcode_platemap_df = pd.read_csv(barcode_path)

    logging.info("Selecting associated plate map")
    try:
        platemap = barcode_platemap_df.query(
            "Assay_Plate_Barcode == @plate"
        ).Plate_Map_Name.values[0]
    except IndexError as e:
        logging.error(
            f"{e} raised. Unable to find associated platemap name with given plate barcode"
        )

    logging.info(f"Identified plate map: {platemap}")
    platemap_file = os.path.join(metadata_dir, "platemap", "{}.csv".format(platemap))
    platemap_df = pd.read_csv(platemap_file)

    logging.info(f"Loading Plate map data from: {sql_file}")
    sqlite_file = "sqlite:///{}".format(sql_file)
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
    logging.info("Counting cells within each well")
    cell_count_df = single_cell_profile.count_cells()

    logging.info(f"Saving cell counts in: {cell_count_out}")
    cell_count_df = cell_count_df.merge(
        platemap_df, left_on="Image_Metadata_Well", right_on="well_position"
    ).drop(["WellRow", "WellCol", "well_position"], axis="columns")
    cell_count_df.to_csv(cell_count_out, sep="\t", index=False)

    # aggregating singel cells into aggregate profile
    logging.info("Aggregating cells")
    single_cell_profile.aggregate_profiles(
        output_file=aggregate_file_out, compression_options="gzip"
    )
    logging.info(f"Aggregate profile saved in : {aggregate_file_out}")


if __name__ == "__main__":

    # snakemake inputs
    plate_data = str(snakemake.input["sql_files"])
    barcode_path = str(snakemake.input["barcodes"])
    metadata_dir_path = str(snakemake.input["metadata"])
    cell_count_output = str(snakemake.output["cell_counts"])
    aggregate_output = str(snakemake.output["aggregate_profile"])
    config_path = str(snakemake.params["aggregate_config"])
    log_path = str(snakemake.log)

    aggregate(
        sql_file=plate_data,
        metadata_dir=metadata_dir_path,
        barcode_path=barcode_path,
        aggregate_file_out=aggregate_output,
        cell_count_out=cell_count_output,
        config=config_path,
        log_file=log_path,
    )
