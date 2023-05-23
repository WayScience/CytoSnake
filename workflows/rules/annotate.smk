"""
rule module: annotate.smk

Utilizes pycytominer's annotate module:
https://github.com/cytomining/pycytominer/blob/master/pycytominer/annotate.py

Annotates profiles with given metadata.

Parameters
----------
input:
    aggregate_profile: aggregated profile dataset

    barcodes: file containing unique barcodes that maps to a specific plate

    metadata: directory containing metadata associated with the aggregate
              profile

output:
    generates an annotated profile.


Returns:
--------
    Generates an annotated profile stored in the `results/` directory
"""
# prioritizing which input to use


rule annotate:
    input:
        profile=get_data_path(
            input_type=config["annotate_configs"]["params"]["input_data"],
            use_converted=DATA_CONFIGS["use_converted_plate_data"],
        ),
        barcodes=BARCODES,
        metadata=METADATA_DIR,
    output:
        ANNOTATED_DATA,
    conda:
        "../envs/cytominer_env.yaml"
    log:
        "logs/annotate_{basename}.log",
    params:
        annotate_config=config["config_paths"]["annotate"],
    script:
        "../scripts/annotate.py"
