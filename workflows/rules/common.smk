"""
rule module: common.smk

common.smk is a workflow module that setups the expected input and output paths
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
# -- generating a wild card list (just the file base names)
plate_name = hf.get_file_basenames(DATA_DIR, ext_target="sqlite")

# -- getting the rest of the input paths from helper functions
PLATE_DATA = hf.get_plate_data()
BARCODES = hf.get_barcodes()
METADATA_DIR = hf.get_metadata_dir()

# -------
# OUTPUTS
# -------
# -- extended = list of the file names with a given wildcard
AGGREGATE_DATA = hf.aggregate_output()
CELL_COUNTS = hf.cell_count_output()

CELL_COUNTS_EXPANDED = expand(CELL_COUNTS, file_name=plate_name)
AGGREGATE_DATA_EXPAND = expand(AGGREGATE_DATA, file_name=plate_name)

ANNOTATED_DATA = hf.annotated_output()
ANNOTATED_DATA_EXPAND = expand(ANNOTATED_DATA, file_name=plate_name)

NORMALIZED_DATA = hf.normalized_output()
NORMALIZED_DATA_EXPAND = expand(NORMALIZED_DATA, file_name=plate_name)

SELECTED_FEATURE_DATA = hf.selected_features_output()
SELECTED_FEATURE_DATA_EXPAND = expand(SELECTED_FEATURE_DATA, file_name=plate_name)

CONSENSUS_DATA = hf.consensus_output()
