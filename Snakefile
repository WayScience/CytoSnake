from pathlib import Path
import glob


# obtaining plate_ids
sql_paths = glob.glob("./data/*.sqlite")
PLATE_IDS = glob_wildcards("data/{id}.sqlite").id


# importing DAGs
include: "rules/preprocessing.smk"
include: "rules/feature_select.smk"


rule all:
    input:
        expand("results/preprocessing/{plate_id}_aggregate.csv.gz", plate_id=PLATE_IDS),  # expected outputs from the first DAG "Preprocessing"
        expand("results/preprocessing/{plate_id}_cell_counts.tsv", plate_id=PLATE_IDS),
        expand("results/preprocessing/{plate_id}_augmented.csv.gz", plate_id=PLATE_IDS),
        expand("results/preprocessing/{plate_id}_normalized.csv.gz", plate_id=PLATE_IDS),
        expand(
            "results/preprocessing/{plate_id}_feature_select.csv.gz",
            plate_id=PLATE_IDS,
        ),
        "results/preprocessing/consensus.tsv.gz",
