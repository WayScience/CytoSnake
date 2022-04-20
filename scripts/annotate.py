import yaml
from pathlib import Path
import pandas as pd
from pycytominer.annotate import annotate


def annotate_cells(aggregated_data, barcodes_path):
    """Annoates

    Parameters
    ----------
    aggregated_data : str
        path pointing to aggregated dataset
    platemaps : str
        path pointing to platemaps

    """
    # loading in paramters
    annotate_ep = Path(snakemake.params["annotate_config"])
    annotate_config_path = annotate_ep.absolute()
    with open(annotate_config_path, "r") as yaml_contents:
        annotate_configs = yaml.safe_load(yaml_contents)["annotate_configs"]["params"]

    # input paths retrived by snakemake
    annotated_outfile = str(snakemake.output)
    platemap_df = pd.read_csv(barcodes_path)

    # annotating the aggregated profiles
    annotate(
        profiles=aggregated_data,
        platemap=platemap_df,
        join_on=annotate_configs["join_on"],
        output_file=annotated_outfile,
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

    # iterative loop
    barcodes = str(snakemake.input["barcodes"])
    aggregated_profiles = str(snakemake.input["aggregate_profile"]).split()
    print("Annotating profiles ...")
    for ap in aggregated_profiles:
        annotate_cells(ap, barcodes)
