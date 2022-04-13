# Documentation:
# Workflow that involves preprocessing raw single-cell plate data and
# transforming it into normalized aggregate profiles.
#
# Parameters
# ---------
# sql_files : List[str]
#   List of SQL files containing plate data
# barcodes : str
#   path pointing to the
# metadata : str
#   path pointing to plate metadata
#
# Generates
# -------
# cell_counts: .csv file
#   csv file containing n_cells per well
# augmented: csv.gz file
#   Annotated aggregated profiles
# normalized: csv.gz file
#   Normalized annotated aggregate profiles
# --------------------


# collecting all unqieuids from plate
configfile: "configs/configuration.yaml"


rule aggregate:
    input:
        sql_files=expand("data/{plate_id}.sqlite", plate_id=PLATE_IDS),
        barcodes="data/barcode_platemap.csv",
        metadata="data/metadata",
    output:
        cell_counts=expand(
            "results/preprocessing/{plate_id}.cell_counts.tsv", plate_id=PLATE_IDS
        ),
        aggregate_profile=expand(
            "results/preprocessing/{plate_id}.aggregate.csv.gz", plate_id=PLATE_IDS
        ),
    conda:
        "../envs/cytominer_env.yaml"
    script:
        "../scripts/aggregate_cells.py"


rule annotate:
    input:
        barcodes="data/barcode_platemap.csv",
        aggregate_profiles=expand(
            "results/preprocessing/{plate_id}.aggregate.csv.gz", plate_id=PLATE_IDS
        ),
    output:
        expand("results/preprocessing/{plate_id}_augmented.csv.gz", plate_id=PLATE_IDS),
    conda:
        "../envs/cytominer_env.yaml"
    script:
        "../scripts/annotate.py"


rule normalize:
    input:
        expand("results/preprocessing/{plate_id}_augmented.csv.gz", plate_id=PLATE_IDS),
    output:
        expand("results/preprocessing/{plate_id}.normalized.csv.gz", plate_id=PLATE_IDS),
    params:
        norm_method=config["Normalization"]["parameters"]["method"],
    conda:
        "../envs/cytominer_env.yaml"
    script:
        "../scripts/normalize.py"
