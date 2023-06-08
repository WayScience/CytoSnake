rule feature_select:
    """
    Performs feature selection based on this given profiles. PyCytominer contains
    different operations to conduct its feature selection: variance_threshold,
    correlation_threshold, drop_na_columns, drop_outliers, and noise_removal.

    Utilizes pycytominer's feature select module:
    https://github.com/cytomining/pycytominer/blob/master/pycytominer/feature_select.py

    :input profile: single-cell morphology dataset.
    :input barcode: file containing unique barcodes that maps to a specific plate.
    :input metadata: metadata file associated with single-cell morphology dataset.

    :config: workflow config pointing to feature selection configs

    :output: selected features profile stored in `results/` directory
    """
    input:
        get_data_path(
            input_type=config["feature_select_configs"]["params"]["input_data"],
            tolist=True,
        ),
    output:
        get_data_path(input_type="feature_select", tolist=True),
    params:
        feature_select_config=config["feature_select_configs"],
    log:
        "logs/feature_select.log",
    conda:
        "../envs/cytominer_env.yaml"
    script:
        "../scripts/feature_select.py"
