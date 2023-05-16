"""
rule module: normalize.smk

Utlizes pycytominer's normalization module:
https://github.com/cytomining/pycytominer/blob/c90438fd7c11ad8b1689c21db16dab1a5280de6c/pycytominer/normalize.py

Normalizing single-cell or aggregate features. Current default normalization
method is `standardize` other methods include:


parameters
----------
input
    single-cell or aggregated profiles

output
    normalized single-cell or aggregate dataset.

Output
------
    Generates an annotated profile stored in the `results/` directory
"""


rule normalize:
    input:
        CYTOTABLE_CONVERTED_PLATE_DATA
        if config["data_configs"]["use_converted_plate_data"]
        else PLATE_DATA,
    output:
        NORMALIZED_DATA,
    conda:
        "../envs/cytominer_env.yaml"
    log:
        "logs/normalized_{basename}.log",
    params:
        normalize_config=config["normalize_configs"],
    script:
        "../scripts/normalize.py"
