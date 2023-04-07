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
        CYTOTABLE_OUTPUT_DATA,
    conda:
        "../envs/cytotable.yaml"
    params:
        data_configs=config["data_configs"]["plate_data_format"],
        cytotable_config=config["cytotable_convert"],
    script:
        "../scripts/convert.py"
