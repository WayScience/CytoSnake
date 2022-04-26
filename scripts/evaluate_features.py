import yaml
from cytominer_eval.evaluate import evaluate


def evaluate_features(profile, features, meta_features, replicate_groups, outname):

    # opening config file
    evaluate_ep = Path(snakemake.params["eval_config"])
    evaluate_config_path = evaluate_ep.absolute()
    with open(evaluate_config_path, "r") as yaml_contents:
        evaluate_config = yaml.safe_load(yaml_contents)["evaluate_configs"]["params"]

    evaluate_df = evaluate(
        profile,
        features,
        meta_features,
        replicate_groups,
        operation=evaluate_config["operation"],
        groupby_columns=evaluate_config["groupby_columns"],
        similarity_metric=evaluate_config["similarity_metric"],
        replicate_reproducibility_quantile=evaluate_config[
            "replicate_reproducibility_quantile"
        ],
        replicate_reproducibility_return_median_cor=evaluate_config[
            "replicate_reproducibility_return_median"
        ],
        precision_recall_k=evaluate_config["precision_recall_k"],
        grit_control_perts=evaluate_config["grit_control_perts"],
        grit_replicate_summary_method=evaluate_config["grit_replicate_summary_method"],
        mp_value_params=evaluate_config["mp_value_params"],
        enrichment_percentile=evaluate_config["enrichment_percentile"],
        hitk_percent_list=evaluate_config["hitk_percent_list"],
    )

    evaluate_df.to_csv(outname, index=False)


if __name__ == "__main__":

    profiles = [str(f_in) for f_in in snakemake.input]
    outnames = [(str(f_out) for f_out in snakemake.output)]
    io_files = zip(profiles, outnames)

    # TODO:
    # there are required inputs for the `evaluate()` function
    # Not implemented yet
    features = None
    meta_features = None
    replicate_groups = None

    for profile, outname in io_files:
        evaluate_features(profile, features, meta_features, replicate_groups, outname)
