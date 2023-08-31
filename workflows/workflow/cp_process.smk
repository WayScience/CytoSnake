"""
workflow: cp_process.smk

Description:
------------
Uses plate data and converts them to consensus profiles

The development of this workflow was heavily influenced this study:
https://github.com/broadinstitute/cell-health

Parameters:
-----------
input:
    PLATE_DATA: single-cell morphology plate datasets
    METADATA: metadata directory associated with plate data
    BARCODE (optional): file containing plate data to plate map pairings
output:
    AGGREGATE_DATA_EXPAND: Aggregated profiles
    CELL_COUNTS_EXPANDED: Cell counts per well
    ANNOTATED_DATA_EXPAND: Annotated profile that contains aggregated + metadata data
    NORMALIZED_DATA_EXPAND: Normalized aggregate datasets
    SELECTED_FEATURE_DATA_EXPAND: Selected features from normalized aggregate datasets
    CONSENSUS_DATA: Profile containing unique signatures

Returns
-------
    Consensus profile containing unique signatures associated with specific treatments.
"""


# import workflow configurations [genreal + workflow config]
configfile: "./configs/configuration.yaml"
configfile: "./configs/wf_configs/cp_process.yaml"


# importing rule modules
include: "../rules/common.smk"
include: "../rules/aggregate.smk"
include: "../rules/annotate.smk"
include: "../rules/normalize.smk"
include: "../rules/feature_select.smk"


# include: "../rules/generate_consensus.smk"


# set expected outputs from workflow
rule all:
    input:
        get_data_path(input_type="aggregated", tolist=True),
        get_data_path(input_type="cell_counts", tolist=True),
        get_data_path(input_type="annotated", tolist=True),
        get_data_path(input_type="normalized", tolist=True),
        get_data_path(input_type="feature_select", tolist=True),
        # get_data_path(input_type="consensus", tolist=True),
