"""
rule module: common.smk

common.smk is a workflow module that sets up the expected input and output paths
for the main analytical workflow.
"""
from snakemake.io import expand
from cytosnake.helpers import helper_funcs as hf
from cytosnake.utils.config_utils import load_data_path_configs, load_meta_path_configs

# Paths
DATA_DIR = str(load_data_path_configs())

# ------
# INPUTS
# ------
# generating a wild card list (just the file base names)
plate_name = hf.get_file_basenames(DATA_DIR, ext_target="sqlite")

# getting the rest of the input paths from helper functions

# level 2 data: Single cell dataused as inputs along with associated
# metadata and barcodes
#
# METADATA_DIR: contains information of what has been added to the cells
# BARCODES: Unique id that points to a specific level 2 dataset (plate data)
PLATE_DATA = hf.get_plate_data()
BARCODES = hf.get_barcodes()
METADATA_DIR = hf.get_metadata_dir()

# -------
# OUTPUTS
# -------
# Things to know:
# extended = list of the file names with a given wildcard
#
# To understand the level of data, please refere to PyCytominer documentation
# https://github.com/cytomining/pycytominer

# Level 2 data: converted into parquet format
CYTOTABLE_OUTPUT_DATA = hf.parquet_output()
CYTOTABLE_OUTPUT_DATA_EXTENDED = expand(CONVERTED_DATA, file_name=plate_name)

# level 2.5 data: annotated level 2 data based on given metadata (e.g treatments)
ANNOTATED_DATA = hf.annotated_output()
ANNOTATED_DATA_EXPAND = expand(ANNOTATED_DATA, file_name=plate_name)

# level 3 data: aggregated profile based on given aggregation level
# (e.g aggregating single-cell data to the well level)
AGGREGATE_DATA = hf.aggregate_output()
AGGREGATE_DATA_EXPAND = expand(AGGREGATE_DATA, file_name=plate_name)

# level 4a data: noramlzied profile
NORMALIZED_DATA = hf.normalized_output()
NORMALIZED_DATA_EXPAND = expand(NORMALIZED_DATA, file_name=plate_name)

# level 4b: selected features profile
SELECTED_FEATURE_DATA = hf.selected_features_output()
SELECTED_FEATURE_DATA_EXPAND = expand(SELECTED_FEATURE_DATA, file_name=plate_name)

# level 5: Consensus profile captures unique signatures that resulted from
# any external factor (e.g pertubations)
CONSENSUS_DATA = hf.consensus_output()

# other outputs
# Cell counts: cell counts per well in level 2 data
CELL_COUNTS = hf.cell_count_output()
CELL_COUNTS_EXPANDED = expand(CELL_COUNTS, file_name=plate_name)
