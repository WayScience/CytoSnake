from pathlib import Path
import yaml
from pycytominer.normalize import normalize


def normalization(anno_file, norm_outfile):
    """Normalizes aggregate profiles

    Parameters
    ----------
    anno_file : str
        path leading to aggregate profiles file
    norm_outfile : str
        output name of the generated normalized file
    norm_method : str
        Method of normalization

    Returns
    -------
    No python object is returned. Generates normalized aggregated profile in the
    results/ directory
    """

    # loading paramters
    normalize_ep = Path(snakemake.params["normalize_config"])
    normalize_config_path = normalize_ep.absolute()
    with open(normalize_config_path, "r") as yaml_contents:
        normalize_config = yaml.safe_load(yaml_contents)["normalize_configs"]["params"]

    meta_features = [
        "Metadata_Plate",
        "Metadata_Well",
        "Metadata_Plate_Map_Name",
        "Metadata_Plate",
        "Metadata_Well",
        "Metadata_Object_Count",
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

    # preprocessing converting snakemake objects into python strings
    annotated_files = [str(f_in) for f_in in snakemake.input]
    out_files = [str(f_out) for f_out in snakemake.output]
    io_files = zip(annotated_files, out_files)

    # iteratively normalizing annotated files
    for annotated_file, out_file in io_files:
        normalization(annotated_file, out_file)
