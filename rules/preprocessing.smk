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
        sql_files="data/{plate_id}.sqlite",
        barcodes="data/barcode_platemap.csv",
        metadata="data/metadata",
    output:
        cell_counts="results/preprocessing/{plate_id}_cell_counts.tsv",
        aggregate_profile="results/preprocessing/{plate_id}_aggregate.csv.gz",
    conda:
        "../envs/cytominer_env.yaml"
    log:
        "logs/aggregate_{plate_id}.log",
    params:
        aggregate_config=config["config_paths"]["single_cell"],
    script:
        "../scripts/aggregate_cells.py"


rule annotate:
    input:
        barcodes="data/barcode_platemap.csv",
        aggregate_profile="results/preprocessing/{plate_id}_aggregate.csv.gz",
        metadata="data/metadata",
    output:
        "results/preprocessing/{plate_id}_augmented.csv.gz",
    conda:
        "../envs/cytominer_env.yaml"
    log:
        "logs/annotate_{plate_id}.log",
    params:
        annotate_config=config["config_paths"]["annotate"],
    script:
        "../scripts/annotate.py"


rule normalize:
    input:
        "results/preprocessing/{plate_id}_augmented.csv.gz",
    output:
        "results/preprocessing/{plate_id}_normalized.csv.gz",
    conda:
        "../envs/cytominer_env.yaml"
    log:
        "logs/normalized_{plate_id}.log",
    params:
        normalize_config=config["config_paths"]["normalize"],
    script:
        "../scripts/normalize.py"
