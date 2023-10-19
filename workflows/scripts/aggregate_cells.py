import logging
import pathlib

import memray
import pandas as pd
from pycytominer.cyto_utils.cells import SingleCells


def aggregate(
    sql_file: str,
    metadata_dir: str,
    barcode_path: str,
    cell_count_out: str,
    aggregate_file_out: str,
    single_cell_config: dict,
    aggregate_config: dict,
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
    single_cell_config: dict
        contains all single_cell parameters
    aggregate_config: dict
        contains all aggregate parameters
    log_file : str
        Path to log file

    Returns:
    --------
        No object is returned
        Generates cell counts and aggregate profile output
    """
    sqlite_file_path = pathlib.Path(sql_file).resolve(strict=True)

    # opening logger
    log_path = pathlib.Path(log_file).absolute()
    logging.basicConfig(
        filename=log_path,
        encoding="utf-8",
        level=logging.DEBUG,
        format="%(asctime)s.%(msecs)03d - %(levelname)s - %(thread)d - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    logging.info(f"Aggregating Single Cells in: {str(sqlite_file_path)}")

    # Loading appropriate platemaps with given plate data
    logging.info(f"Loading plate data from: {str(sqlite_file_path)}")

    # checking if the the sqli_file is found or it's a file type
    if not pathlib.Path(sql_file).is_file():
        e_msg = f"Unable to find plate file: {sql_file}"
        logging.error(e_msg)
        raise FileNotFoundError(e_msg)

    # loading in barcodes
    logging.info(f"Loading barcodes from: {barcode_path}")
    barcode_platemap_df = pd.read_csv(barcode_path)

    # plate map search
    logging.info("Selecting associated plate map")
    try:
        # plate variable is used in the query string below

        plate = pathlib.Path(plate_data).stem  # noqa
        platemap = barcode_platemap_df.query(
            "Assay_Plate_Barcode == @plate"
        ).Plate_Map_Name.values[0]
    except IndexError as e:
        logging.error(
            f"{e} raised."
            "Unable to find associated platemap name with given plate barcode"
        )

    # loading in identified platemap
    logging.info(f"Identified plate map: {platemap}")
    platemap_file = pathlib.Path(f"{metadata_dir}/platemap/{platemap}.csv")
    platemap_df = pd.read_csv(platemap_file)

    # Loading single-cell data into the SingleCell object
    logging.info(f"Loading Plate map data from: {sql_file}")
    sqlite_file = f"sqlite:///{plate_data}"
    single_cell_profile = SingleCells(
        sqlite_file,
        strata=single_cell_config["strata"],
        image_cols=single_cell_config["image_cols"],
        aggregation_operation=single_cell_config["aggregation_operation"],
        output_file=aggregate_file_out,
        merge_cols=single_cell_config["merge_cols"],
        add_image_features=single_cell_config["add_image_features"],
        image_feature_categories=single_cell_config["image_feature_categories"],
        features=single_cell_config["features"],
        load_image_data=single_cell_config["load_image_data"],
        subsample_frac=single_cell_config["subsample_frac"],
        subsampling_random_state=single_cell_config["subsampling_random_state"],
        fields_of_view=single_cell_config["fields_of_view"],
        fields_of_view_feature="Image_Metadata_Well",
        object_feature=single_cell_config["object_feature"],
    )

    # counting cells in each well and saving it as csv file
    logging.info("Counting cells within each well")
    cell_count_df = single_cell_profile.count_cells()

    logging.info(f"Saving cell counts in: {cell_count_out}")
    cell_count_df = cell_count_df.merge(
        platemap_df, left_on="Image_Metadata_Well", right_on="well_position"
    ).drop(["WellRow", "WellCol", "well_position"], axis="columns")
    cell_count_df.to_csv(cell_count_out, sep="\t", index=False)

    # Using the SingleCell object, aggregate single-cells to aggregate profiles
    logging.info("Aggregating cells")
    single_cell_profile.aggregate_profiles(
        output_file=aggregate_file_out,
        compute_subsample=aggregate_config["compute_subsample"],
        compression_options=aggregate_config["compression_options"],
        float_format=aggregate_config["float_format"],
        n_aggregation_memory_strata=aggregate_config["n_aggregation_memory_strata"],
    )
    logging.info(f"Aggregate profile saved in : {aggregate_file_out}")


# execute main code for aggregation
if __name__ == "__main__":
    # snakemake inputs
    # more information how snakemake transfers workflow variables to scripts:
    # https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#python
    plate_data = str(snakemake.input["sql_files"])
    barcode_path = str(snakemake.input["barcodes"])
    metadata_dir_path = str(snakemake.input["metadata"])
    cell_count_output = str(snakemake.output["cell_counts"])
    aggregate_output = str(snakemake.output["aggregate_profile"])
    single_cells_config = snakemake.params["single_cell_config"]["params"]
    aggregate_config = snakemake.params["aggregate_config"]["params"]
    log_path = str(snakemake.log)
    enable_profiling = snakemake.config["enable_profiling"]

    # exeucuting pycytominer aggregate function
    if enable_profiling:
        # anything below this context manager will be profiled
        root_name = pathlib.Path(plate_data).stem.split("_")[0]
        with memray.Tracker(f"{root_name}_aggregate_benchmark.bin"):
            aggregate(
                sql_file=plate_data,
                metadata_dir=metadata_dir_path,
                barcode_path=barcode_path,
                aggregate_file_out=aggregate_output,
                cell_count_out=cell_count_output,
                single_cell_config=single_cells_config,
                aggregate_config=aggregate_config,
                log_file=log_path,
            )
    else:
        aggregate(
            sql_file=plate_data,
            metadata_dir=metadata_dir_path,
            barcode_path=barcode_path,
            aggregate_file_out=aggregate_output,
            cell_count_out=cell_count_output,
            single_cell_config=single_cells_config,
            aggregate_config=aggregate_config,
            log_file=log_path,
        )
