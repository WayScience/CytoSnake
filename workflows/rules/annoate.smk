"""
rule module: annotate.smk

Utlizes pycytominer's annotate module:
https://github.com/cytomining/pycytominer/blob/c90438fd7c11ad8b1689c21db16dab1a5280de6c/pycytominer/annotate.py

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


rule annotate:
    input:
        aggregate_profile=AGGREGATE_DATA,
        barcodes=BARCODES,
        metadata=METADATA_DIR,
    output:
        ANNOTATED_DATA,
    conda:
        "../envs/cytominer_env.yaml"
    log:
        "logs/annotate_{file_name}.log",
    params:
        annotate_config=config["config_paths"]["annotate"],
    script:
        "../scripts/annotate.py"
