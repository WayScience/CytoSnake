"""
This workflow focuses in processing dp profiler features.
"""


configfile: "configs/configuration.yaml"


rule dp_aggregate_profiles:
    input:
        index_file=INDEX_FILE,
        dp_features_dir="data/{dp_name}_dp",
    output:
        "results/processing/{dp_name}_dp_aggregated.csv.gz",
    conda:
        "../envs/dp_process.yaml"
    params:
        dp_data_configs=config["config_paths"]["dp_data"],
        aggregator_configs=config["config_paths"]["dp_aggregator"],
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
        "results/processing/{dp_name}_dp_consensus.csv.gz",
    conda:
        "../envs/dp_process.yaml"
    params:
        consensus_config=config["config_paths"]["consensus_config"],
    script:
        "../scripts/build_dp_consensus.py"
