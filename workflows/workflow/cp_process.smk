import glob
from cytosnake.helpers import helper_funcs as hf


# importing rule modules
include: "../rules/common.smk"
include: "../rules/aggregate.smk"
include: "../rules/annotate.smk"
include: "../rules/normalize.smk"
include: "../rules/feature_select.smk"
include: "../rules/generate_consensus.smk"


# expected outputs from workflow
rule all:
    input:
        AGGREGATE_DATA_EXPAND,
        CELL_COUNTS_EXPANDED,
        ANNOTATED_DATA_EXPAND,
        NORMALIZED_DATA_EXPAND,
        SELECTED_FEATURE_DATA_EXPAND,
        CONSENSUS_DATA,
