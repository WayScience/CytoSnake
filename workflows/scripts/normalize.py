import logging
import pathlib

from pycytominer.normalize import normalize


def normalization(
    anno_file: str, norm_outfile: str, config: str, log_file: str
) -> None:
    """Normalizes aggregate profiles

    Parameters
    ----------
    anno_file : str
        path leading to aggregate profiles file
    norm_outfile : str
        output name of the generated normalized file
    config : str
        Path to normalization config file
    log_file : str
        Path to log file

    Returns
    -------
    No python object is returned. Generates normalized aggregated profile in the
    results/ directory
    """

    # initiating logger
    log_path = pathlib.Path(log_file).absolute()
    logging.basicConfig(
        filename=log_path,
        encoding="utf-8",
        level=logging.DEBUG,
        format="%(asctime)s.%(msecs)03d - %(levelname)s - %(thread)d - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logging.info("Starting Annotation processes")

    # loading parameters
    logging.info(f"Loading Annotation configuration from: {config}")
    normalize_config = config["params"]

    # normalizing annotated aggregated profiles
    logging.info(f"Normalizing annotated aggregated profiles: {anno_file}")
    normalize(
        anno_file,
        features=normalize_config["features"],
        image_features=normalize_config["image_features"],
        meta_features=normalize_config["meta_features"],
        samples=normalize_config["samples"],
        method=normalize_config["method"],
        output_file=norm_outfile,
        compression_options=normalize_config["compression_options"],
        float_format=normalize_config["float_format"],
        mad_robustize_epsilon=normalize_config["mad_robustize_epsilon"],
        spherize_center=normalize_config["spherize_center"],
        spherize_method=normalize_config["spherize_method"],
        spherize_epsilon=normalize_config["spherize_epsilon"],
    )
    logging.info(f"Normalized aggregated profile saved: {norm_outfile}")


# executes normalization protocol
if __name__ == "__main__":
    # snakemake inputs
    # more information how snakemake transfers workflow variables to scripts:
    # https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#python
    annotated_data_path = str(snakemake.input)
    normalize_configs = snakemake.params["normalize_config"]
    normalized_data_output = str(snakemake.output)
    log_path = str(snakemake.log)

    # normalization step
    normalization(
        anno_file=annotated_data_path,
        norm_outfile=normalized_data_output,
        config=normalize_configs,
        log_file=log_path,
    )
