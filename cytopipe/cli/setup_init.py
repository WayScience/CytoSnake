"""setup_init.py

Contains function of different variation of initializing data files
for processing
"""

import shutil
from typing import Union
from pathlib import Path


def init_cp_data(data_fp: Union[list[str], str], metadata_fp: str, barcode_fp: str):
    """Sets up directory for Cell Profiler datasets

    Parameters
    ----------
    data_fp : Union[list[str], str]
        Plate data files
    metadata_fp : str
        metadata directory associated with plate data
    barcodes_fp : str
        barcode file associated with plate data
    """

    # setting up paths
    barcode_path = str(Path(barcode_fp).absolute())
    metadata_path = str(Path(metadata_fp).absolute())

    # create data folder on working directory
    # raises error if data directory exists (prevents overwriting)
    data_dir_obj = Path("data")
    data_dir_obj.mkdir(exist_ok=False)
    data_dir_path = str(data_dir_obj.absolute())

    # moving all user provided plate data files into data folder
    for data_file in data_fp:
        f_path = str(Path(data_file).absolute())
        shutil.move(f_path, data_dir_path)

    # user provided barcode file and metadata directory is
    # moved to the data folder
    shutil.move(barcode_path, data_dir_path)
    shutil.move(metadata_path, data_dir_path)


def init_dp_data(data_fp: Union[list[str], str], metadata_fp: str):
    """Sets up directory for Deep Profiler datasets

    Parameters
    ----------
    data_fp : Union[list[str], str]
        directory containing extracted Deep Profiler features
    metadata_fp : str
        metadata directory associated with Deep Profiler data

    Returns
    -------
    None
        Generates a directory containing Deep Profiler Datasets
    """

    # setting up files
    metadata_path = str(Path(metadata_fp).absolute())

    # create data folder on working directory
    # raises error if data directory exists (prevents overwriting)
    data_dir_obj = Path("data")
    data_dir_obj.mkdir(exist_ok=False)
    data_dir_path = str(data_dir_obj.absolute())

    # moving all user provided plate data files into data folder
    for data_file in data_fp:

        # adding 'dp' suffix to the directory containing deep profiler datasets
        f_path_obj = Path(data_file).absolute()
        renamed_dirpath = f"{f_path_obj.parent}/{f_path_obj.name}_dp"

        # rename and move to data directory
        f_path_obj.rename(renamed_dirpath)
        shutil.move(renamed_dirpath, data_dir_path)

    # user metadata directory is
    # moved to the data folder
    shutil.move(metadata_path, data_dir_path)