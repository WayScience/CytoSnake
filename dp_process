"""
Analyzing morphological fetures extracted from DeepProfiler
"""
from pathlib import Path

# creating workflow varaible for index file
IDX_FILE_NAME = glob_wildcards("data/metadata/{indx_name}.csv").indx_name
if len(IDX_FILE_NAME) != 1:
    raise ValueError("Multiple index files were identified, please remove one")
IDX_FILE_NAME = IDX_FILE_NAME[0]

INDEX_FILE = f"data/metadata/{IDX_FILE_NAME}.csv"

# generating wildcard for directories containing dp features
DP_DIR = glob_wildcards("data/{dp_dir}_dp").dp_dir


# including processing modules
include: "rules/dp_process.smk"


# # expected output files
rule all:
    input:
        expand("results/processing/{dp_name}_dp_aggregated.csv.gz", dp_name=DP_DIR),
        expand(
            "results/processing/{dp_name}_dp_normalized_aggregated.csv.gz",
            dp_name=DP_DIR,
        ),
        expand("results/processing/{dp_name}_dp_consensus.csv.gz", dp_name=DP_DIR),
