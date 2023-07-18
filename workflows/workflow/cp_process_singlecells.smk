"""
workflow: cp_process_singlecells.smk

Description:
------------
Converts sqlite plate data into parquet and returns selected features from
single-cell morphology profiles.

The development of this workflow was heavily influenced by:
https://github.com/WayScience/Benchmarking_NF1_data/tree/main/3_extracting_features

Parameters:
----------
input:
    PLATE_DATA: single-cell morphology plate in sqlite format
    METADATA: associated metadata
output:
    CYTOTABLE_OUTPUT_DATA_EXTENDED: Converted single-cell morphology datasets.
    ANNOTATED_DATA_EXPAND: Profile containing metadata.
    NORMALIZED_DATA_EXPAND: Normalized single-cell morphology datasets.
    SELECTED_FEATURE_DATA_EXPAND: Selected features from normalized single-cell
    morphology datasets.

Returns
-------
    Workflow generates selected features profile from single-cell morphology datasets
"""


# importing workflow configs [general + workflow config]
configfile: "./configs/configuration.yaml"
configfile: "./configs/wf_configs/cp_process_singlecells.yaml"


# importing modules
include: "../rules/common.smk"
include: "../rules/cytotable_convert.smk"
include: "../rules/annotate.smk"
include: "../rules/normalize.smk"
include: "../rules/feature_select.smk"


# define expected outputs
rule all:
    input:
        get_data_path(
            input_type="plate_data",
            tolist=True,
            use_converted=config["data_configs"]["use_converted_plate_data"],
        ),
        get_data_path(input_type="annotated", tolist=True),
        get_data_path(input_type="normalized", tolist=True),
        get_data_path(input_type="feature_select", tolist=True),
