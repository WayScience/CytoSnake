import glob
from cytosnake.helpers import helper_funcs as hf

# obtaining plate_ids
PLATE_IDS = glob_wildcards("data/{id}.sqlite").id


# importing DAGs
include: "../rules/common.smk"
include: "../rules/preprocessing.smk"
include: "../rules/feature_select.smk"


# appending logs
rule all:
    input:
        AGGREGATE_DATA_EXPAND,
        CELL_COUNTS_EXPANDED,
        ANNOTATED_DATA_EXPAND,
        NORMALIZED_DATA_EXPAND,
        SELECTED_FEATURE_DATA_EXPAND,
