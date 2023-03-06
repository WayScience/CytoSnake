configfile: "configs/configuration.yaml"


rule feature_select:
    input:
        NORMALIZED_DATA_EXPAND,
    output:
        SELECTED_FEATURE_DATA_EXPAND,
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
        SELECTED_FEATURE_DATA_EXPAND,
    output:
        CONSENSUS_DATA,
    params:
        consensus_configs=config["config_paths"]["consensus_config"],
    log:
        "logs/create_consensus.log",
    conda:
        "../envs/cytominer_env.yaml"
    script:
        "../scripts/consensus.py"
