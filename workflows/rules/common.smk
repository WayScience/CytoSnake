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
BASENAME = datapaths.get_all_basenames(DATA_DIR, ext_target="sqlite")

# getting the rest of the input paths from helper functions

# level 2 data: Single cell dataused as inputs along with associated
# metadata and barcodes
#
# METADATA_DIR: contains information of what has been added to the cells
# BARCODES: Unique id that points to a specific level 2 dataset (plate data)
PLATE_DATA = datapaths.build_path(input_type="plate_data")
BARCODES = datapaths.get_barcodes()
METADATA_DIR = datapaths.get_metadata_dir()


# helper function
def get_data_path(
    input_type: str,
    use_converted: Optional[bool] = False,
    tolist: Optional[bool] = False,
) -> str:
    """Returns absolute path of an input.

    To understand the level of data, please refere to PyCytominer documentation
    https://github.com/cytomining/pycytominer

    The supported input types are:

    plate_data (level 2):
        Single-cell morphology datasets
    annotated (level 2.5):
        Datasets that have been augmented with metadata. This provides additional
        information to the profiles like perturbations types, cell types,
        well locations etc.
    aggregated (level 3):
        Aggregated single-cell morphology datasets.
    normalized (level 4a):
        Normalized profiles.
    feature_select (level 4b):
        selected features from given profile
    consensus (level 5)
        profile containing unique signatures pertaining to specific perturbations
        and/or external factors.

    Other paths generated:

    cell_counts:
        refers to the number of cells within a well

    Parameters
    ----------
    input_type: str
        data type to be used in order to generate a path

    use_conveted: Optional[bool]
        flag to return path converted single-cell data

    tolist: Optional[bool]
        flag to return a list of paths of the desired data type

    Returns
    -------
    str
        returns path of select data type input path`
    """
    # checking if provided `input_type` is valid
    if input_type not in DATA_CONFIGS["data_types"].keys():
        raise TypeError(f"`{input_type}` is not a supported datatype.")

    # checking if the user wants to use converted dataset
    # if plate data is selected, then check if the converted path is required
    if input_type == "plate_data" and use_converted:
        if tolist:
            data_path = datapaths.build_path(
                input_type=input_type, use_converted=use_converted
            )
            return expand(data_path, basename=BASENAME)

    # build path without converted data path
    data_path = datapaths.build_path(input_type=input_type, use_converted=use_converted)

    # check if the user want a list of paths or a single path
    if tolist:
        expanded_data = expand(data_path, basename=BASENAME)
        return expanded_data

    return data_path
