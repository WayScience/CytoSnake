import os
import argparse 
from pathlib import Path

import pandas as pd
from pycytominer.cyto_utils.cells import SingleCells


def get_profiles(plate, sql_file, metadata_dir, barcode_platemap_df, strata, image_cols, compression):
    """
    Apply all profiling steps for a given plate.
    Output:
    Will write a series of processed files to disk
    """
    print("Processing {}.....".format(plate))

    if not os.path.isfile(sql_file):
        raise FileNotFoundError("Unable to find SQL file")

    sqlite_file = "sqlite:///{}.sqlite".format(plate)

    # Load specific platemap
    barcode_platemap_df = pd.read_csv(barcode_platemap_df)
    platemap = barcode_platemap_df.query(
        "Assay_Plate_Barcode == @plate"
    ).Plate_Map_Name.values[0]
    platemap_file = os.path.join(metadata_dir, "platemap", "{}.csv".format(platemap))
    platemap_df = pd.read_csv(platemap_file)

    # Prepare sql file for processing
    ap = SingleCells(
        sqlite_file, strata=strata, image_cols=image_cols, fields_of_view_feature=[]
    )

    # Count cells and output
    print("Counting cells")
    results_path = Path("results")
    results_path.mkdir(exist_ok=True)
    cell_count_file = os.path.join("results", "{}_cell_count.tsv".format(plate))
    cell_count_df = ap.count_cells()

    print("Saving cell counts in: {}".format(cell_count_file))
    cell_count_df = cell_count_df.merge(
        platemap_df, left_on="Image_Metadata_Well", right_on="well_position"
    ).drop(["WellRow", "WellCol", "well_position"], axis="columns")
    cell_count_df.to_csv(cell_count_file, sep="\t", index=False)

    # Begin processing profiles
    print("Compressing profiles...")
    data_path = Path("data")
    data_path.mkdir(exist_ok=True)
    output_dir = os.path.join("data", "profiles", plate)
    os.makedirs(output_dir, exist_ok=True)

    # Aggregate single cells into well profiles
    out_file = os.path.join(output_dir, "{}.csv.gz".format(plate))
    ap.aggregate_profiles(output_file=out_file, compression_options=compression)
    print("compressed file written on:")


if __name__ == "__main__":

    # CLI Arguments
    parser = argparse.ArgumentParser(description="Aggregation data")
    required = parser.add_argument_group("Required Arguments")
    optional = parser.add_argument_group("Optional Arguments")
    required.add_argument("-s", "--sql_file", type=str, 
                            help="SQL file")
    required.add_argument("-p", "--profile_platemap", type=str,
                            help="Single cell profile")
    required.add_argument("-m", "--metadata_dir", type=str, 
                            help="Cell metadata file")
    optional.add_argument("--strata", nargs="+", type=str, default=["Image_Metadata_Plate", "Image_Metadata_Well"], 
                          help="None")
    optional.add_argument("--image_cols", nargs="+", type=str, default=["TableNumber", "ImageNumber", "Metadata_Site"],
                          help="None")
    optional.add_argument("--compression", type=str, default="gzip", choices=["gizp"], 
                          help="Method of compression for aggregate profiles")
    args = parser.parse_args()


    #extracting platename name from sql file
    plate_name = os.path.splitext(os.path.basename(args.sql_file))[0]

    get_profiles(plate=plate_name,
                 sql_file=args.sql_file,
                 metadata_dir=args.metadata_dir,
                 barcode_platemap_df=args.profile_platemap,
                 strata=args.strata,
                 image_cols = args.image_cols,
                 compression=args.compression
                 )



