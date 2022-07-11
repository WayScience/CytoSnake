import argparse
from pathlib import Path

import yaml
from pycytominer.normalize import normalize


def normalization(anno_file: str, norm_outfile: str, config: str) -> None:
    """Normalizes aggregate profiles

    Parameters
    ----------
    anno_file : str
        path leading to aggregate profiles file
    norm_outfile : str
        output name of the generated normalized file
    config : str
        Path to normalization config file

    Returns
    -------
    No python object is returned. Generates normalized aggregated profile in the
    results/ directory
    """

    # loading parameters
    normalize_ep = Path(config)
    normalize_config_path = normalize_ep.absolute()
    with open(normalize_config_path, "r") as yaml_contents:
        normalize_config = yaml.safe_load(yaml_contents)["normalize_configs"]["params"]

    meta_features = [
        "Metadata_Plate",
        "Metadata_Well",
        "Metadata_WellRow",
        "Metadata_WellCol",
        "Metadata_gene_name",
        "Metadata_pert_name",
        "Metadata_broad_sample",
        "Metadata_cell_line",
    ]

    normalize(
        anno_file,
        features=normalize_config["features"],
        image_features=normalize_config["image_features"],
        # meta_features=normalize_config["meta_features"],
        meta_features=meta_features,
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


if __name__ == "__main__":

    # CLI Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", type=str, required=True, help="annotated aggregated profiles"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="normalized annotated aggregate profiles outputs",
    )
    parser.add_argument(
        "-c", "--config", type=str, required=True, help="path to config file"
    )
    args = parser.parse_args()

    # normalization step
    normalization(anno_file=args.input, norm_outfile=args.output, config=args.config)
