configfile: "configs/configuration.yaml"

rule feature_select:
    input:
        expand("results/preprocessing/{plate_id}_normalized.csv.gz", plate_id=PLATE_IDS),
    output:
        expand("results/preprocessing/{plate_id}_feature_select.csv.gz", plate_id=PLATE_IDS),
    params:
        feature_select_config=config["config_paths"]["feature_select"]
    conda:
        "../envs/cytominer_env.yaml"
    script:
        "../scripts/feature_select.py"


rule evaluate_features:
    input:
        expand("results/preprocessing/{plate_id}_feature_select.csv.gz", plate_id=PLATE_IDS),
    output:
        expand("results/preprocessing/{plate_id}_evaluated.csv", plate_id=PLATE_IDS),
    params:
        eval_config=config["config_paths"]["evaluate"]
    conda:
        "../envs/cytominer_env.yaml"
    script:
        "../scripts/evaluate_features.py"
