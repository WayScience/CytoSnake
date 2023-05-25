rule normalize:
    """
    Normalizing single-cell or aggregate features. Current default normalization
    method is `standardize` other methods include:

    Utlizes pycytominer's normalization module:
    https://github.com/cytomining/pycytominer/blob/c90438fd7c11ad8b1689c21db16dab1a5280de6c/pycytominer/normalize.py

    :input profile: feature selected profiles.
    :input barcode: file containing unique barcodes that maps to a specific plate.
    :input metadata: metadata file associated with single-cell morphology dataset.

    :config: workflow config pointing to feature selection configs

    :output: normalized profile stored in the `results/` directory
    """
    input:
        get_data_path(
            input_type=config["normalize_configs"]["params"]["input_data"],
            use_converted=DATA_CONFIGS["use_converted_plate_data"],
        ),
    output:
        get_data_path(input_type="normalized"),
    conda:
        "../envs/cytominer_env.yaml"
    log:
        "logs/normalized_{basename}.log",
    params:
        normalize_config=config["normalize_configs"],
    script:
        "../scripts/normalize.py"
