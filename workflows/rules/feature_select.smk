"""
rule module: feature_select.smk

Utilizes pycytominer's feature select module:
https://github.com/cytomining/pycytominer/blob/master/pycytominer/feature_select.py

Performs feature selection based on this given profiles. PyCytominer contains
different operations to conduct its feature selection: variance_threshold,
correlation_threshold, drop_na_columns, drop_outliers, and noise_removal.

Parameters:
-----------
Input:
    Cell morphology profiles
Output:
    Selected features from profiles

Returns
-------
    CSV file containing selected features. Stored in the `results/` directory.
"""


rule feature_select:
    input:
        NORMALIZED_DATA_EXPAND,
    output:
        SELECTED_FEATURE_DATA_EXPAND,
    params:
        feature_select_config=config["feature_select_configs"],
    log:
        "logs/feature_select.log",
    conda:
        "../envs/cytominer_env.yaml"
    script:
        "../scripts/feature_select.py"
