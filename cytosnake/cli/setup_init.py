"""setup_init.py

Contains functions for initializing different data structures
for processing
"""

from pathlib import Path
from typing import Union


def init_cp_data(data_fp: Union[list[str], str], metadata_fp: str, barcode_fp: str):
    """Sets up directory for CellProfiler datasets. Symlinks are created by taking
    the path of the inputs files and storing them into the "data" directory.

    Parameters
    ----------
    data_fp : Union[list[str], str]
        Plate data files
    metadata_fp : str
        metadata directory associated with plate data
    barcodes_fp : str
        barcode file associated with plate data
    """

    # checking that none are NoneTypes
    if not any([data_fp, metadata_fp, barcode_fp]):
        raise ValueError(
            "Inputs required do not much match cell_profiler datatype,"
            "please make sure to provide data plates, barcode and metadata"
        )

    # setting up paths
    # -- barcode
    barcode_path_obj = barcode_fp
    if barcode_fp is not None:
        barcode_path_obj = Path(barcode_fp)

    # -- metadata_path
    metadata_path_obj = Path(metadata_fp)

    # create data folder on working directory
    # raises error if data directory exists (prevents overwriting)
    data_dir_obj = Path("data")
    data_dir_obj.mkdir(exist_ok=False)

    # generating symlinks
    for data_file in data_fp:
        data_file_obj = Path(data_file)
        target_file = Path(f"../{data_file}")
        sym_link = Path(f"./{str(data_dir_obj)}/{data_file_obj.name}")
        sym_link.symlink_to(target_file)

    # generating symlinks of provided barcode file and metadata directory is
    if barcode_fp is not None:
        barcode_target = Path(f"../{barcode_path_obj}")
        barcode_symlink = Path(f"./{str(data_dir_obj)}/{barcode_path_obj.name}")
        barcode_symlink.symlink_to(barcode_target)

    metadata_target = Path(f"../{metadata_path_obj}")
    metadata_symlink = Path(f"./{str(data_dir_obj)}/{metadata_path_obj.name}")
    metadata_symlink.symlink_to(metadata_target)


def init_dp_data(data_fp: Union[list[str], str], metadata_fp: str):
    """Sets up directory for CellProfiler datasets. Symlinks are created by
    taking the path of the inputs files and storing them into the "data" directory.

    Parameters
    ----------
    data_fp : Union[list[str], str]
        directory containing extracted DeepProfiler features
    metadata_fp : str
        metadata directory associated with DeepProfiler data

    Returns
    -------
    None
        Generates a directory containing DeepProfiler Datasets
    """

    # setting up files
    metadata_path_obj = Path(metadata_fp)

    # create data folder on working directory
    # raises error if data directory exists (prevents overwriting)
    data_dir_obj = Path("data")
    data_dir_obj.mkdir(exist_ok=False)

    # generating symlink of input data files
    for data_file in data_fp:
        data_file_obj = Path(data_file)
        data_file_target = Path(f"../{data_file}")
        data_file_symlink = Path(f"./{str(data_dir_obj)}/{data_file_obj.name}_dp")
        data_file_symlink.symlink_to(data_file_target)

    # creating symlink of metadata dir to data directory
    metadata_target = Path(f"../{metadata_path_obj}")
    metadata_symlink = Path(f"./{str(data_dir_obj)}/{metadata_path_obj.name}")
    metadata_symlink.symlink_to(metadata_target)
