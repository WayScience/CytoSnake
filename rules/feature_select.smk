configfile: "configs/configuration.yaml"


rule feature_select:
    input:
        expand("results/preprocessing/{plate_id}_normalized.csv.gz", plate_id=PLATE_IDS),
    output:
        expand(
            "results/preprocessing/{plate_id}_feature_select.csv.gz",
            plate_id=PLATE_IDS,
        ),
    params:
        feature_select_config=config["config_paths"]["feature_select"],
    log:
        "logs/feature_select.log",
    conda:
        "../envs/cytominer_env.yaml"
    script:
        "../scripts/feature_select.py"


rule create_consensus:
    input:
        expand(
            "results/preprocessing/{plate_id}_feature_select.csv.gz",
            plate_id=PLATE_IDS,
        ),
    output:
        "results/preprocessing/consensus.tsv.gz",
    params:
        consensus_configs=config["config_paths"]["consensus_config"],
    log:
        "logs/create_consensus.log",
    conda:
        "../envs/cytominer_env.yaml"
    script:
        "../scripts/consensus.py"
