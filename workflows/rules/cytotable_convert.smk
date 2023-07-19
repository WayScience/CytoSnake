rule convert:
    """
    Converts single-cell morphology dataset to parquet format.

    Utilizes CytoTable's convert workflow module:
    https://github.com/cytomining/CytoTable/blob/main/cytotable/convert.py

    :input profiles: single-cell profiles.
    :input barcode: file containing unique barcodes that maps to a specific plate.
    :input metadata: metadata file associated with single-cell morphology dataset.

    :config: workflow config pointing to convert configs.

    :output annotated: converted single-cell morphology dataset.
    """
    input:
        get_data_path(input_type=config["cytotable_convert"]["params"]["input_data"]),
    output:
        get_data_path(
            input_type="plate_data",
            use_converted=config["data_configs"]["use_converted_plate_data"],
        ),
    conda:
        "../envs/cytotable.yaml"
    params:
        data_configs=config["data_configs"]["data_types"]["plate_data"]["converted_ext"],
        cytotable_config=config["cytotable_convert"],
    script:
        "../scripts/convert.py"
