"""
rule module: aggregate.smk


Utilize's pycytominer's aggregate module:
https://github.com/cytomining/pycytominer/blob/c90438fd7c11ad8b1689c21db16dab1a5280de6c/pycytominer/aggregate.py

Aggregates single-cell profiles into aggregated profiles based on a given strata

For example, users can configure `Metadata_Well` as their strata in order to
aggregate single-cell data into the Well level.

Parameters:
-----------
input:
  sql_file: single-cell dataset

  barcodes: file containing unique barcodes that maps to a specific plate

  metadata: directory containing metadata associated with the aggregate
            profile
output:
  aggregated_profile: aggregated profiles
  cell_counts: CSV file that contains how many cells were counted per well

Returns
-------
  aggregated profiles and cell count data stored in the `results/` directory
# --------------------
"""


rule aggregate:
    input:
        sql_files=get_input(
            input_type=config["aggregate_configs"]["params"]["input_data"],
            use_converted=DATA_CONFIGS["use_converted_plate_data"],
        ),
        barcodes=BARCODES,
        metadata=METADATA_DIR,
    output:
        aggregate_profile=AGGREGATE_DATA,
        cell_counts=CELL_COUNTS,
    log:
        "logs/aggregate_{file_name}.log",
    conda:
        "../envs/cytominer_env.yaml"
    params:
        aggregate_config=config["aggregate_configs"]["params"],
    script:
        "../scripts/aggregate_cells.py"
