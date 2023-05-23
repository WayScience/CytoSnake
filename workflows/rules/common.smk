"""
rule module: common.smk

common.smk is a workflow module that sets up the expected input and output paths
for the main analytical workflow.
"""
from typing import Optional
from snakemake.io import expand
from cytosnake.helpers import datapaths
from cytosnake.utils.config_utils import load_data_path_configs, load_general_configs

# Paths
CYTOSNAKE_CONFIGS = load_general_configs()
DATA_DIR = str(load_data_path_configs())
DATA_CONFIGS = CYTOSNAKE_CONFIGS["data_configs"]

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
# EXTENDED is a list of file names that has been extended by using the wildcard.
# In this case, the wildcard placeholder is `basename` and is being extended by
# all basenames of the plate_data
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
AGGREGATE_DATA = datapaths.build_path(data_type="aggregated")
AGGREGATE_DATA_EXPAND = expand(AGGREGATE_DATA, basename=plate_name)

# level 4a data: noramlzied profile
NORMALIZED_DATA = datapaths.build_path(data_type="normalized")
NORMALIZED_DATA_EXPAND = expand(NORMALIZED_DATA, basename=plate_name)

# level 4b: selected features profile
SELECTED_FEATURE_DATA = datapaths.build_path(data_type="feature_select")
SELECTED_FEATURE_DATA_EXPAND = expand(SELECTED_FEATURE_DATA, basename=plate_name)

# level 5: Consensus profile captures unique signatures that resulted from
# any external factor (e.g pertubations)
CONSENSUS_DATA = datapaths.build_path(data_type="consensus")
CONSENSUS_DATA_EXPAND = expand(CONSENSUS_DATA, basename=plate_name)

# other outputsa
# Cell counts: cell counts per well in level 2 data
CELL_COUNTS = datapaths.build_path(data_type="consensus")
CELL_COUNTS_EXPANDED = expand(CELL_COUNTS, basename=plate_name)


def get_input(
    data_type: str,
    use_converted: Optional[bool] = False,
    tolist: Optional[bool] = False,
) -> str:
    """Returns absolute path of an input

    Parameters
    ----------
    data_type: str
        data

    use_conveted: Optional[bool]
        flag to return path conveted single-cell data

    tolist: Optional[bool]
        flag to return a list of paths of the desired data type

    Returns
    -------
    str
        returns path of select data type input path`
    """

    # checking
    if data_type not in DATA_CONFIGS["data_types"].keys():
        raise ValueError(f"`{data_type}` is not a supported datatype.")

    match data_type:
        case "cell_counts":
            return
        case "plate_data":
            if use_converted:
                return CYTOTABLE_CONVERTED_PLATE_DATA
            return PLATE_DATA
        case "aggregated":
            if tolist:
                return AGGREGATE_DATA_EXPAND
            return AGGREGATE_DATA
        case "annotated":
            if tolist:
                return ANNOTATED_DATA_EXPAND
            return ANNOTATED_DATA
        case "normalized":
            if tolist:
                return NORMALIZED_DATA_EXPAND
            return NORMALIZED_DATA
        case "feature_select":
            if tolist:
                return SELECTED_FEATURE_DATA_EXPAND
            return SELECTED_FEATURE_DATA
        case "consensus":
            if tolist:
                return CONSENSUS_DATA_EXPAND
            return CONSENSUS_DATA
