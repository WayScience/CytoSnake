import os
import logging
import yaml
from pathlib import Path

import snakemake
import pandas as pd
from pycytominer.annotate import annotate


def annotate_cells(
    aggregated_data: str,
    barcodes_path: str,
    metadata_dir: str,
    annotate_file_out: str,
    config: str,
    log_file: str,
) -> None:
    """Annotates aggregated profiles

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
    # initiating Logger
    log_path = Path(log_file).absolute()
    logging.basicConfig(
        filename=log_path,
        encoding="utf-8",
        level=logging.DEBUG,
        format="%(asctime)s.%(msecs)03d - %(levelname)s - %(thread)d - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logging.info("Starting Annotation processes")

    # loading in annotate configs
    logging.info(f"Loading Annotation configuration from: {config}")

    annotate_path_obj = Path(config)
    if not annotate_path_obj.is_file():
        e_msg = "Unable to find Annotation configuration file"
        logging.error(e_msg)
        raise FileNotFoundError(e_msg)

    annotate_config_path = annotate_path_obj.absolute()
    with open(annotate_config_path, "r") as yaml_contents:
        annotate_configs = yaml.safe_load(yaml_contents)["annotate_configs"]["params"]
        logging.info(f"Annotation configuration loaded")

    # loading in plate map
    logging.info(f"Loading plate data from: {aggregated_data}")

    logging.info(f"loadding barcodes from: {barcodes_path}")
    barcode_platemap_df = pd.read_csv(barcodes_path)

    logging.info("Searching plate map name")
    plate = Path(aggregated_data).name.split("_")[0]
    platemap = barcode_platemap_df.query(
        "Assay_Plate_Barcode == @plate"
    ).Plate_Map_Name.values[0]

    # checking if plate map is found
    if platemap == "" or platemap is None:
        e_msg = "Unable to find associated plate map"
        logging.error(e_msg)
        raise ValueError(e_msg)

    logging.info(f"Plate map found: {platemap}")

    # loading for plate file
    logging.info("loading plate map file")
    platemap_file = os.path.join(metadata_dir, "platemap", "{}.csv".format(platemap))
    if not Path(platemap_file).is_file():
        e_msg = "Unable to locate plate map file"
        logging.error(e_msg)
        raise FileNotFoundError(e_msg)

    # reading plate map
    platemap_df = pd.read_csv(platemap_file)
    logging.info(f"Loaded plate map: {platemap_file}")

    # annotating the aggregated profiles
    logging.info("Annotating Aggregated profiles")
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
    logging.info(f"Annotated files saved: {annotate_file_out}")


if __name__ == "__main__":

    # snakemake inputs
    aggregate_data_path = str(snakemake.input["aggregate_profile"])
    annotate_data_output = str(snakemake.output)
    barcode_path = str(snakemake.input["barcodes"])
    metadata_dir_path = str(snakemake.input["metadata"])
    config_path = str(snakemake.params["annotate_config"])
    log_path = str(snakemake.log)

    # annotating cells
    annotate_cells(
        aggregated_data=aggregate_data_path,
        barcodes_path=barcode_path,
        metadata_dir=metadata_dir_path,
        annotate_file_out=annotate_data_output,
        config=config_path,
        log_file=log_path,
    )
