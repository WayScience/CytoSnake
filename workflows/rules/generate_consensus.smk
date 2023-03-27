"""
rule module: generate_consensus.smk

Utilize's pycytominer's consensus module:
https://github.com/cytomining/pycytominer/blob/master/pycytominer/consensus.py

Creates consensus profies that reflects unique signatures associated with
external factors.

Parameters:
----------
input:
    Selected features profile
output:
    Consensus profile

Return:
-------
    Consensus profile stored in the `results/` directory
"""


configfile: "configs/configuration.yaml"


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
