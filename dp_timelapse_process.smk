"""
Analyzing time laps data captured in deep profiler features
"""
import os
# creating workflow varaible for index file
# INDEX_FILE="data/metadata/{index}.csv"
# IDX_NAME = glob_wildcards("data/metadata/{index}.csv").index

IDX_FILE = "data/metadata/index.csv"
# generating wildcard for directories containing dp features
DP_DIR = glob_wildcards("data/{dp_dir}_dp").dp_dir


# including processing modules
include: "rules/dp_process.smk"


# # expected output files
rule all:
    input:
        expand("results/processing/{dp_name}_dp_aggregated.csv.gz", dp_name=DP_DIR),
        expand("results/processing/{dp_name}_dp_normalized_aggregated.csv.gz", dp_name=DP_DIR),
        expand("results/processing/{dp_name}_dp_consensus.csv.gz", dp_name=DP_DIR)
