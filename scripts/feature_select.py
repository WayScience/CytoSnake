from pathlib import Path
import yaml
from pycytominer.feature_select import feature_select


def feature_selection(normalized_profile, out_file):
    """Performs feature selection based on the given parameters explained
    in the configs/analysis_configs/feature_selection_configs.yaml file.

    Parameters
    ----------
    normalized_profile : str
        path that points to normalized profile

    Returns
    -------
    Generates output
    """

    # loading paramters
    feature_select_ep = Path(snakemake.params["feature_select_config"])
    feature_select_config_path = feature_select_ep.absolute()
    with open(feature_select_config_path, "r") as yaml_contents:
        feature_select_config = yaml.safe_load(yaml_contents)["feature_select_configs"][
            "params"
        ]

    feature_select(
        normalized_profile,
        features=feature_select_config["features"],
        image_features=feature_select_config["image_features"],
        samples=feature_select_config["samples"],
        operation=feature_select_config["operation"],
        na_cutoff=feature_select_config["na_cutoff"],
        corr_threshold=feature_select_config["corr_threshold"],
        corr_method=feature_select_config["corr_method"],
        freq_cut=feature_select_config["freq_cut"],
        unique_cut=feature_select_config["unique_cut"],
        compression_options=feature_select_config["compression_options"],
        float_format=feature_select_config["float_format"],
        blocklist_file=feature_select_config["blocklist_file"],
        outlier_cutoff=feature_select_config["outlier_cutoff"],
        noise_removal_perturb_groups=feature_select_config[
            "noise_removal_perturb_groups"
        ],
        noise_removal_stdev_cutoff=feature_select_config["noise_removal_stdev_cutoff"],
        output_file=out_file,
    )


if __name__ == "__main__":
    all_norm_profile = [str(f_in) for f_in in snakemake.input]
    out_files = [str(f_out) for f_out in snakemake.output]
    inputs = zip(all_norm_profile, out_files)

    # iteratively passing normalized data
    for norm_data, out_file in inputs:
        feature_selection(norm_data, out_file)
