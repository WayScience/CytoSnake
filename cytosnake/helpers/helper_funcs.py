from cytosnake.utils.cyto_paths import get_data_path


def get_data() -> str:
    """Returns absolute path were the data is located

    Returns
    -------
    str
        returns a string that contains a snakemake wild card that collects all 
        plate data
        
    """

    # get the datapath
    data_path = get_data_path()
    return str(data_path) + "/{plate_name}.sqlite"


# def get_metadata_dir 