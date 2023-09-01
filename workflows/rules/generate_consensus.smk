rule create_consensus:
    """
    Creates consensus profiles that reflects unique signatures associated with
    external factors.

    Utilize's pycytominer's consensus module:
    https://github.com/cytomining/pycytominer/blob/master/pycytominer/consensus.py

    :input profile: selected features profiles
    :input barcode: file containing unique barcodes that maps to a specific plate.
    :input metadata: metadata file associated with single-cell morphology dataset.

    :config: workflow config pointing to feature selection configs

    :output: Consensus profile stored in the `results/` directory

    """
    input:
        get_data_path(
            input_type=config["consensus_config"]["params"]["input_data"], tolist=True
        ),
    output:
        get_data_path(input_type="consensus"),
    params:
        consensus_configs=config["config_paths"]["consensus_config"],
    log:
        "logs/create_consensus.log",
    conda:
        "../envs/cytominer_env.yaml"
    script:
        "../scripts/consensus.py"
