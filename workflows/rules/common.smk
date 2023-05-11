"""
rule module: common.smk

common.smk is a workflow module that sets up the expected input and output paths
for the main analytical workflow.
"""
from snakemake.io import expand
from cytosnake.helpers import datapaths
from cytosnake.utils.config_utils import load_data_path_configs

# Paths
DATA_DIR = str(load_data_path_configs())

# ------
# INPUTS
# ------
# generating a wild card list (just the file base names)
plate_name = datapaths.get_all_basenames(DATA_DIR, ext_target="sqlite")

# getting the rest of the input paths from helper functions

# level 2 data: Single cell dataused as inputs along with associated
# metadata and barcodes
#
# METADATA_DIR: contains information of what has been added to the cells
# BARCODES: Unique id that points to a specific level 2 dataset (plate data)
PLATE_DATA = datapaths.build_path(data_type="plate_data")
BARCODES = datapaths.get_barcodes()
METADATA_DIR = datapaths.get_metadata_dir()

# -------
# OUTPUTS
# -------
# Things to know:
# extended = list of the file names with a given wildcard
#
# To understand the level of data, please refere to PyCytominer documentation
# https://github.com/cytomining/pycytominer

# Level 2 data: converted into parquet format
CYTOTABLE_CONVERTED_PLATE_DATA = datapaths.build_path(
    data_type="plate_data", use_converted=True
)
CYTOTABLE_CONVERTED_PLATE_DATA_EXTENDED = expand(
    CYTOTABLE_CONVERTED_PLATE_DATA, basename=plate_name
)

# level 2.5 data: annotated level 2 data based on given metadata (e.g treatments)
ANNOTATED_DATA = datapaths.build_path(data_type="annotated")
ANNOTATED_DATA_EXPAND = expand(ANNOTATED_DATA, basename=plate_name)

# level 3 data: aggregated profile based on given aggregation level
# (e.g aggregating single-cell data to the well level)
AGGREGATE_DATA = datapaths.aggregate_output(data_type="aggregated")
AGGREGATE_DATA_EXPAND = expand(AGGREGATE_DATA, basename=plate_name)

# level 4a data: noramlzied profile
NORMALIZED_DATA = datapaths.normalized_output(data_type="normalized")
NORMALIZED_DATA_EXPAND = expand(NORMALIZED_DATA, basename=plate_name)

# level 4b: selected features profile
SELECTED_FEATURE_DATA = datapaths.build_path(data_type="feature_select")
SELECTED_FEATURE_DATA_EXPAND = expand(SELECTED_FEATURE_DATA, basename=plate_name)

# level 5: Consensus profile captures unique signatures that resulted from
# any external factor (e.g pertubations)
CONSENSUS_DATA = datapaths.consensus_output(data_type="consensus")
CONSENSUS_DATA_EXPAND = expand(CONSENSUS_DATA, basename=plate_name)

# other outputs
# Cell counts: cell counts per well in level 2 data
CELL_COUNTS = datapaths.build_path(data_type="consensus")
CELL_COUNTS_EXPANDED = expand(CELL_COUNTS, basename=plate_name)
