"""
rule module: cytotable_convert.smk

Utilizes CytoTable's convert workflow module:
https://github.com/cytomining/CytoTable/blob/main/cytotable/convert.py

Parameters:
-----------


Returns:
--------
    parquet files stored within the data/ folder

"""


configfile: "configs/configuration.yaml"


rule convert:
    input:
        PLATE_DATA,
    output:
        CONVERTED_DATA,
    params:
        cytotable_config=config["config_paths"]["cytotable_config"],
    script:
        "../scripts/convert.py"
