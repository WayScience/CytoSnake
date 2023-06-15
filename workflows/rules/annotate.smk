rule annotate:
    """
    Generates an annotated profile with given metadata and is stored instored
    in the `results/` directory.

    Utilizes pycytominer's annotate module:
    https://github.com/cytomining/pycytominer/blob/master/pycytominer/annotate.py

    :input profiles: single-cell or aggregate profiles.
    :input barcode: file containing unique barcodes that maps to a specific plate.
    :input metadata: metadata file associated with single-cell morphology dataset.

    :config: workflow config pointing to annotate configs.

    :output annotated: annotated profile.
    """
    input:
        profile=get_data_path(
            input_type=config["annotate_configs"]["params"]["input_data"],
            use_converted=DATA_CONFIGS["use_converted_plate_data"],
        ),
        barcodes=BARCODES,
        metadata=METADATA_DIR,
    output:
        get_data_path(input_type="annotated"),
    conda:
        "../envs/cytominer_env.yaml"
    log:
        "logs/annotate_{basename}.log",
    params:
        annotate_config=config["config_paths"]["annotate"],
    script:
        "../scripts/annotate.py"
