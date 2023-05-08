"""
#

workflow: cp_process.smk


Description:
------------
Uses plate data and converts them to consensus profiles

The development of this workflow was heavily influced this study:
https://github.com/broadinstitute/cell-health

Parameters:
-----------
input:
    PLATE_DATA: single-cell morphology plate dataset
    METADATA: metadata directory asscociated with plate data
    BARCODE (optional): file containing plate data to plate map pairings
output:
    AGGREGATE_DATA_EXPAND: Aggregated profiles
    CELL_COUNTS_EXPANDED: Cell counts per well
    ANNOTATED_DATA_EXPAND: Annotated profile that contains aggregated + metadata data
    NORMALIZED_DATA_EXPAND: Normalized aggregate dataset
    SELECTED_FEATURE_DATA_EXPAND: Selected features from normalized aggregate dataset
    CONSENSUS_DATA: Profile containing unique signatures

Returns
-------
    Consensus profle containing unique signatures associated with specific treatments.
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
include: "../rules/generate_consensus.smk"


# set expected outputs from workflow
rule all:
    input:
        AGGREGATE_DATA_EXPAND,
        CELL_COUNTS_EXPANDED,
        ANNOTATED_DATA_EXPAND,
        NORMALIZED_DATA_EXPAND,
        SELECTED_FEATURE_DATA_EXPAND,
        CONSENSUS_DATA,
