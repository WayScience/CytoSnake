"""setup_init.py

Contains functions for initializing different data structures
for processing
"""

import pathlib
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
        barcode_path_obj = pathlib.Path(barcode_fp).resolve(strict=True)

    # -- metadata_path
    metadata_path_obj = pathlib.Path(metadata_fp).resolve(strict=True)

    # create data folder on working directory
    # raises error if data directory exists (prevents overwriting)
    data_dir_obj = pathlib.Path("data").resolve()
    data_dir_obj.mkdir(exist_ok=False)

    # generating symlinks
    for data_file in data_fp:
        data_file_obj = pathlib.Path(data_file).resolve(strict=True)
        sym_link = data_dir_obj / data_file_obj.name
        sym_link.symlink_to(data_file_obj)

    # generating symlinks of provided barcode file and metadata directory is
    if barcode_fp is not None:
        barcode_symlink = pathlib.Path(data_dir_obj) / barcode_path_obj.name
        barcode_symlink.symlink_to(barcode_path_obj)

    metadata_symlink = pathlib.Path(data_dir_obj) / metadata_path_obj.name
    metadata_symlink.symlink_to(metadata_path_obj)


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
    metadata_path_obj = pathlib.Path(metadata_fp).resolve(strict=True)

    # create data folder on working directory
    # raises error if data directory exists (prevents overwriting)
    data_dir_obj = pathlib.Path("data")
    data_dir_obj.mkdir(exist_ok=False)

    # generating symlink of input data files
    for data_file in data_fp:
        data_file_obj = pathlib.Path(data_file).resolve(strict=True)
        data_file_symlink = data_dir_obj / f"{data_file_obj.name}_dp"
        data_file_symlink.symlink_to(data_file_obj)

    # creating symlink of metadata dir to data directory
    metadata_symlink = data_dir_obj / metadata_path_obj.name
    metadata_symlink.symlink_to(metadata_path_obj)
