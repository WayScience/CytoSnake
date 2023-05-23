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
