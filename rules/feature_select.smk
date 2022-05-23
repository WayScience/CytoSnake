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
        "results/preprocessing/median_consensus.csv.gz",
    params:
        aggregate_config=config["config_paths"]["aggregate"],
    conda:
        "../envs/cytominer_env.yaml"
    script:
        "../scripts/consensus.py"
