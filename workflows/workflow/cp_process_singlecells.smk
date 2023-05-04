"""
workflow: cp_process_singlecells.smk

Description:
------------
Converts sqlite plate data into parquet and returns selected features in csv
format

Parameters:
----------
input:
    Plate data in sqlite format
output:
    Selected features in csv format

Returns
-------
    Selected morphological features
"""


# importing workflow configs [general + workflow config]
configfile: "./configs/configuration.yaml"
configfile: "./configs/wf_configs/cp_process_singlecells.yaml"


# importing modules
include: "../rules/common.smk"
include: "../rules/cytotable_convert.smk"
include: "../rules/normalize.smk"
include: "../rules/feature_select.smk"


rule all:
    input:
        CYTOTABLE_OUTPUT_DATA_EXTENDED,
        NORMALIZED_DATA_EXPAND,
        SELECTED_FEATURE_DATA_EXPAND,
