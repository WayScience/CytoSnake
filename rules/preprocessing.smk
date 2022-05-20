# Documentation:
# Workflow that involves preprocessing raw single-cell plate data and
# transforming it into normalized aggregate profiles.
#
# Parameters
# ----------
# sql_files : List[str]
#   List of SQL files containing plate data
# barcodes : str
#   path pointing to the barcode file storing platemap IDs
# metadata : str
#   path pointing to plate metadata
#
# Generates
# ---------
# cell_counts: .csv file
#   csv file containing n_cells per well
# augmented: csv.gz file
#   Annotated aggregated profiles
# normalized: csv.gz file
#   Normalized annotated aggregate profiles
# --------------------


# collecting all unique IDs from plate
configfile: "configs/configuration.yaml"


rule aggregate:
    input:
        sql_files=expand("data/{plate_id}.sqlite", plate_id=PLATE_IDS),
        barcodes="data/barcode_platemap.csv",
        metadata="data/metadata",
    output:
        cell_counts=expand(
            "results/preprocessing/{plate_id}_cell_counts.tsv", plate_id=PLATE_IDS
        ),
        aggregate_profile=expand(
            "results/preprocessing/{plate_id}_aggregate.csv.gz", plate_id=PLATE_IDS
        ),
    conda:
        "../envs/cytominer_env.yaml"
    params:
        aggregate_config=config["config_paths"]["single_cell"],
    threads: config["analysis_configs"]["preprocessing"]["threads"]
    script:
        "../scripts/aggregate_cells.py"


rule annotate:
    input:
        barcodes="data/barcode_platemap.csv",
        aggregate_profile=expand(
            "results/preprocessing/{plate_id}_aggregate.csv.gz", plate_id=PLATE_IDS
        ),
        metadata="data/metadata",
    output:
        expand("results/preprocessing/{plate_id}_augmented.csv.gz", plate_id=PLATE_IDS),
    conda:
        "../envs/cytominer_env.yaml"
    params:
        annotate_config=config["config_paths"]["annotate"],
    script:
        "../scripts/annotate.py"


rule normalize:
    input:
        expand("results/preprocessing/{plate_id}_augmented.csv.gz", plate_id=PLATE_IDS),
    output:
        expand("results/preprocessing/{plate_id}_normalized.csv.gz", plate_id=PLATE_IDS),
    conda:
        "../envs/cytominer_env.yaml"
    params:
        normalize_config=config["config_paths"]["normalize"],
    script:
        "../scripts/normalize.py"
