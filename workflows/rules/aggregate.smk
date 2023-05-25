rule aggregate:
    """
    Utilize's pycytominer's aggregate module:
    https://github.com/cytomining/pycytominer/blob/c90438fd7c11ad8b1689c21db16dab1a5280de6c/pycytominer/aggregate.py

    Aggregates single-cell profiles into aggregated profiles based on a given strata

    For example, users can configure `Metadata_Well` as their strata in order to
    aggregate single-cell data into the Well level.

    :input profile: single-cell morphology dataset.
    :input barcode: file containing unique barcodes that maps to a specific plate.
    :input metadata: metadata file associated with single-cell morphology dataset.

    :config: workflow config file pointing to aggregate configurations step

    :output aggregated-profile: aggregated datsaset
    :output cell-counts: csv file containg cell counts per well
    """
    input:
        sql_files=get_data_path(
            input_type=config["aggregate_configs"]["params"]["input_data"],
            use_converted=DATA_CONFIGS["use_converted_plate_data"],
        ),
        barcodes=BARCODES,
        metadata=METADATA_DIR,
    output:
        aggregate_profile=get_data_path(input_type="aggregated", tolist=True),
        cell_counts=CELL_COUNTS,
    log:
        "logs/aggregate_{file_name}.log",
    conda:
        "../envs/cytominer_env.yaml"
    params:
        aggregate_config=config["aggregate_configs"]["params"],
    script:
        "../scripts/aggregate_cells.py"
