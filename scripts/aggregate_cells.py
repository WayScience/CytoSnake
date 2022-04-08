import os
import faulthandler
import pandas as pd
from pycytominer.cyto_utils.cells import SingleCells


def aggregate(sql_file, metadata_dir, barcode_platemap):
    """_summary_

    Args:
        sql_file (str): _description_
        metadata_dir (str): _description_
        barcode_platemap (str): _description_
    """

    # Loading appropriate platemaps with given plate data
    plate = os.path.basename(sql_file).rsplit(".", 1)
    barcode_platemap_df = pd.read_csv(barcode_platemap)
    platemap = barcode_platemap_df.query(
        "Assay_Plate_Barcode == @plate"
    ).Plate_Map_Name.values[0]
    platemap_file = os.path.join(metadata_dir, "platemap", "{}.csv".format(platemap))
    platemap_df = pd.read_csv(platemap_file)

    # TODO: generate a configuration file that contains the default parameters
    # -- Configurations will be in the snakefile under the "param:" rule attribute
    sqlite_file = "sqlite:///{}".format(sql_file)
    strata = ["Image_Metadata_Plate", "Image_Metadata_Well"]
    image_cols = ["TableNumber", "ImageNumber"]
    single_cell_profile = SingleCells(
        sqlite_file,
        strata=strata,
        image_cols=image_cols,
        fields_of_view_feature="Image_Metadata_Well",
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

    # aggregating single cell profiles into well profiles
    # TODO: Figure out why SEGFAULTS are being raised here
    faulthandler.enable()
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
