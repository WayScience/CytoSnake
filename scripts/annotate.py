from ast import arg
import os
from numpy import require
import yaml
import argparse
from pathlib import Path

import pandas as pd

# pycytominer imports
from pycytominer.annotate import annotate


def annotate_cells(
    aggregated_data: str,
    barcodes_path: str,
    metadata_dir: str,
    annotate_file_out: str,
    config: str,
) -> None:
    """Annoates

    Parameters
    ----------
    aggregated_data : str
        path pointing to aggregated dataset
    barcodes_path : str
        path pointing to platemaps
    metadata_dir : str
        path to metadata folder
    annotate_file_out: str
        name of generated annotated profile
    configs:
        path to configuration file
    """

    # loading in annotate configs
    annotate_path_obj = Path(config)
    annotate_config_path = annotate_path_obj.absolute()
    with open(annotate_config_path, "r") as yaml_contents:
        annotate_configs = yaml.safe_load(yaml_contents)["annotate_configs"]["params"]

    # loading in platmap
    plate = Path(aggregated_data).name.split("_")[0]
    barcode_platemap_df = pd.read_csv(barcodes_path)
    platemap = barcode_platemap_df.query(
        "Assay_Plate_Barcode == @plate"
    ).Plate_Map_Name.values[0]
    platemap_file = os.path.join(metadata_dir, "platemap", "{}.csv".format(platemap))
    platemap_df = pd.read_csv(platemap_file)

    # annotating the aggregated profiles
    annotate(
        profiles=aggregated_data,
        platemap=platemap_df,
        join_on=annotate_configs["join_on"],
        output_file=annotate_file_out,
        add_metadata_id_to_platemap=annotate_configs["add_metadata_id_to_platemap"],
        format_broad_cmap=annotate_configs["format_broad_cmap"],
        clean_cellprofiler=annotate_configs["clean_cellprofiler"],
        external_metadata=annotate_configs["external_metadata"],
        external_join_left=annotate_configs["external_join_left"],
        external_join_right=annotate_configs["external_join_right"],
        compression_options=annotate_configs["compression_options"],
        float_format=annotate_configs["float_format"],
        cmap_args=annotate_configs["cmap_args"],
    )


if __name__ == "__main__":

    # CLI
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", type=str, required=True, help="aggregated profiles"
    )
    parser.add_argument(
        "-o", "--output", type=str, required=True, help="Annotated profiles output"
    )
    parser.add_argument(
        "-b", "--barcode", type=str, required=True, help="Path to barcode file"
    )
    parser.add_argument(
        "-m",
        "--metadata_dir",
        type=str,
        required=True,
        help="Path of metadata directory",
    )
    parser.add_argument(
        "-c", "--config", type=str, required=True, help="Path to config file"
    )
    args = parser.parse_args()

    annotate_cells(
        aggregated_data=args.input,
        barcodes_path=args.barcode,
        metadata_dir=args.metadata_dir,
        annotate_file_out=args.output,
        config=args.config,
    )
