import yaml
from pathlib import Path
import pandas as pd
from pycytominer.annotate import annotate


def annotate_cells(
    aggregated_data: str, barcodes_path: str, annotate_file_out: str, config: str
) -> None:
    """Annoates

    Parameters
    ----------
    aggregated_data : str
        path pointing to aggregated dataset
    barcodes_path : str
        path pointing to platemaps
    output: str
        name of generated annotated profile
    configs:
        path to configuration file
    """

    # loading in paramters
    annotate_path_obj = Path(config)
    annotate_config_path = annotate_path_obj.absolute()
    with open(annotate_config_path, "r") as yaml_contents:
        annotate_configs = yaml.safe_load(yaml_contents)["annotate_configs"]["params"]

    # input paths retrived by snakemake
    platemap_df = pd.read_csv(barcodes_path)

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

    # iterative loop
    barcodes = str(snakemake.input["barcodes"])
    config_path = Path(snakemake.params["annotate_config"])

    aggregated_profiles = [
        str(agg_prof) for agg_prof in snakemake.input["aggregate_profile"]
    ]
    outname = [str(f_out) for f_out in snakemake.output]
    inputs = zip(aggregated_profiles, outname)

    print("Annotating profiles ...")
    for agg_profile, out_name in inputs:
        annotate_cells(
            aggregated_data=agg_profile,
            barcodes_path=barcodes,
            annotate_file_out=out_name,
            config=config_path,
        )
