import os
import glob

sql_paths = glob.glob("./data/*.sqlite")
PLATE_IDS = [os.path.basename(fpath).rsplit(".", 1)[0] for fpath in sql_paths]

include: "rules/preprocessing.smk"

rule all:
    input:
        # expected outputs from the first DAG "Preprocessing"
        expand("results/preprocessing/{plate_id}.aggregate.csv.gz", plate_id=PLATE_IDS),
        expand("results/preprocessing/{plate_id}.cell_counts.tsv", plate_id=PLATE_IDS),
        expand("results/preprocessing/{plate_id}_augmented.csv.gz", plate_id=PLATE_IDS),
        expand("results/preprocessing/{plate_id}.normalized.csv.gz", plate_id=PLATE_IDS)

