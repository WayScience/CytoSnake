"""
Documentation:
Workflow that involves preprocessing raw single-cell plate data and
transforming it into normalized aggregate profiles.

Parameters
----------
sql_files : List[str]
  List of SQL files containing plate data
barcodes : str
  path pointing to the barcode file storing platemap IDs
metadata : str
  path pointing to plate metadata

Generates
---------
cell_counts: .csv file
  csv file containing n_cells per well
augmented: csv.gz file
  Annotated aggregated profiles
normalized: csv.gz file
  Normalized annotated aggregate profiles
# --------------------
"""


# collecting all unique IDs from plate
configfile: "configs/configuration.yaml"


rule aggregate:
    input:
        sql_files=PLATE_DATA,
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
        aggregate_config=config["config_paths"]["single_cell"],
    script:
        "../scripts/aggregate_cells.py"


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


rule normalize:
    input:
        ANNOTATED_DATA,
    output:
        NORMALIZED_DATA,
    conda:
        "../envs/cytominer_env.yaml"
    log:
        "logs/normalized_{file_name}.log",
    params:
        normalize_config=config["config_paths"]["normalize"],
    script:
        "../scripts/normalize.py"
