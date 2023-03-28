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


rule convert:
    input:
        PLATE_DATA,
        METADATA_DIR,
        BARCODES,
    output:
        OUTPUTS,
    script:
        "../scripts/"
