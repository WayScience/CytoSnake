from pycytominer.normalize import normalize


def normalization(anno_file, norm_outfile, norm_method):
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
    meta_features = [
        "Metadata_Plate",
        "Metadata_Well",
        "Metadata_Plate_Map_Name",
        "Metadata_Plate",
        "Metadata_Well",
        "Metadata_Object_Count",
    ]

    normalize(
        profiles=anno_file,
        features="infer",
        meta_features=meta_features,
        samples="all",
        method=norm_method,
        output_file=norm_outfile,
        compression_options="gzip",
    )


if __name__ == "__main__":

    # preprocessing converting snakemake objects into python strings
    annotated_files = [str(f_in) for f_in in snakemake.input]
    out_files = [str(f_out) for f_out in snakemake.output]
    method = str(snakemake.params["norm_method"])

    io_files = zip(annotated_files, out_files)

    methods = ["standardize", "robustize", "mad_robustize", "spherize"]
    if method not in methods:
        raise ValueError(f"Unsupported method provided. Supported method: {methods}")

    # iteratively normalizing annotated files
    for annotated_file, out_file in io_files:
        normalization(annotated_file, out_file, norm_method=method)
