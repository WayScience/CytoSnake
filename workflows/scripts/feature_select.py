import logging
from pathlib import Path

import yaml
from pycytominer.feature_select import feature_select


def feature_selection(
    normalized_profile: str, out_file: str, config: str, log_file: str
) -> None:
    """Performs feature selection based on the given parameters explained
    in the configs/analysis_configs/feature_selection_configs.yaml file.

    Parameters
    ----------
    normalized_profile : str
        Path that points to normalized profile
    out_file : str
        Name of generated outfile
    config: str
        Path pointing to config file
    log_file : str
        Path pointing to log file.

    Returns
    -------
        Generates a file containing the selected features.
    """

    # initiating logger
    log_path = Path(log_file).absolute()
    logging.basicConfig(
        filename=log_path,
        encoding="utf-8",
        level=logging.DEBUG,
        format="%(asctime)s.%(msecs)03d - %(levelname)s - %(thread)d - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logging.info("Starting Feature selection")

    # loading configs
    logging.info(f"Loading feature selection configuration from: {config}")

    feature_select_config = config["params"]
    logging.info("Feature Selection configuration loaded")

    # Feature selection
    logging.info("Conducting feature selection")
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
    logging.info(f"Selected features saved: {out_file}")


if __name__ == "__main__":
    all_norm_profile = [str(f_in) for f_in in snakemake.input]
    out_files = [str(f_out) for f_out in snakemake.output]
    config_path = snakemake.params["feature_select_config"]
    io_files = zip(all_norm_profile, out_files)
    log_path = str(snakemake.log)

    # iteratively passing normalized data
    for norm_data, feature_file_out in io_files:
        feature_selection(
            normalized_profile=norm_data,
            out_file=feature_file_out,
            config=config_path,
            log_file=log_path,
        )
