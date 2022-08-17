"""
This workflow focuses in processing dp profiler features.
"""


configfile: "configs/configuration.yaml"


rule dp_aggregate_profiles:
    input:
        idx_file = IDX_FILE,
        dp_features_dir="data/{dp_name}_dp",
    output:
        agg_out="results/processing/{dp_name}_dp_aggregated.csv.gz",
    conda:
        "../envs/dp_process.yaml"
    params:
        dp_data_config=config["config_paths"]["dp_data"],
        dp_agg_config=config["config_paths"]["dp_aggregation"],
    script:
        "../scripts/dp_aggregate.py"


rule normalize_aggregate:
    input:
        "results/processing/{dp_name}_dp_aggregated.csv.gz",
    output:
        "results/processing/{dp_name}_dp_normalized_aggregated.csv.gz",
    conda:
        "../envs/dp_process.yaml"
    params:
        normalize_config=config["config_paths"]["normalize"],
    script:
        "../scripts/dp_normalize.py"


rule build_consensus:
    input:
        "results/processing/{dp_name}_dp_normalized_aggregated.csv.gz",
    output:
        "results/processing/{dp_name}_dp_consensus.csv.gz"
    conda:
        "../envs/dp_process.yaml"
    params:
        consensus_config=config["config_paths"]["consensus_config"],
    script:
        "../scripts/dp_build_consensus.py"
