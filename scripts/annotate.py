import pandas as pd
from pycytominer.annotate import annotate


def annotate_cells(aggregated_data, platemaps):
    """Annoates

    Parameters
    ----------
    aggregated_data : str
        path pointing to aggregated dataset
    platemaps : str
        path pointing to platemaps

    """

    # input paths retrived by snakemake
    annotated_outfile = str(snakemake.output)
    platemap_df = pd.read_csv(barcodes)

    # annotating the aggregated profiles
    # NOTE: join_on should be in config file
    annotate(
        profiles=aggregated_profiles,
        platemap=platemap_df,
        join_on=["Metadata_Assay_Plate_Barcode", "Image_Metadata_Plate"],
        output_file=annotated_outfile,
        compression_options="gzip",
    )


if __name__ == "__main__":

    # iterative loop
    barcodes = str(snakemake.input["barcodes"])
    aggregated_profiles = str(snakemake.input["aggregate_profiles"])
    print("Annotaing profiles ...")
    for ap in aggregated_profiles:
        annotate_cells(ap, barcodes)
