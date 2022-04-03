
import os
import pandas as pd
from pycytominer.cyto_utils.cells import SingleCells

# inputs
sql_file = "data/SQ00014613.sqlite"
metadata_dir = "data/metadata/"
barcode_platemap = "data/barcode_platemap.csv"

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
    sqlite_file, strata=strata, image_cols=image_cols, fields_of_view_feature="Image_Metadata_Well"
)

print("Counting cells within each well")
cell_count_outfile = "cellcount.csv"
cell_count_df = single_cell_profile.count_cells()
print(cell_count_df)

print("Saving cell counts in: {}".format(cell_count_outfile))
cell_count_df = cell_count_df.merge(
    platemap_df, left_on="Image_Metadata_Well", right_on="well_position"
).drop(["WellRow", "WellCol", "well_position"], axis="columns")
cell_count_df.to_csv(cell_count_outfile, sep="\t", index=False) 

aggregate_outfile = "hash.test.csv.gz"
single_cell_profile.aggregate_profiles(output_file=aggregate_outfile, compression_options="gzip")

# ap.to_csv("in_mem_data.csv.gz", compression="gzip")