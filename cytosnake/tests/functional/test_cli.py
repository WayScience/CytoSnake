"""
module: test_cli.py

Description:

The test_cli.py module is a crucial part of the CytoSnake project, responsible for
conducting functional tests on CytoSnake's Command Line Interface (CLI). These tests aim
to ensure that the CLI functions correctly, handling both positive and negative
scenarios effectively.

The purpose of this module is to verify that the CytoSnake CLI operates as expected and
provides accurate results when interacting with various user parameters and modes. It
validates that the CLI can handle different inputs, execute successfully in positive
cases, and properly report errors in negative cases.

The test scenarios are categorized into two main types:

Positive Cases:
These test cases validate the expected behavior of the CLI when
provided with valid user parameters. The primary goal is to ensure that the CLI executes
successfully without errors and produces the correct output.

Negative Cases:
These test cases simulate scenarios where invalid user parameters are provided.
The objective is to confirm that the CLI can identify and handle errors appropriately,
providing informative error messages to the user.
"""

import os
import pathlib
import shutil
import subprocess
from dataclasses import dataclass

import pytest

from cytosnake.common import errors
from cytosnake.tests.test_utils import get_raised_error


# --------------
# helper functions
# --------------
@dataclass
class DataFiles:
    """Structured datatype representation that contains all files in a selected dataset

    Attributes:
    -----------
    metdata: str | list[str]
        metadata directory name
    plate_data: list[str]
        list of plate data (parquet or sqlite files)
    barcode : str
        barcode file name

    Returns
    -------
    DataFiles
        DataStructure representation of test dataset files
    """

    # required parameters
    dataset_dir: str | pathlib.Path

    # extracted files
    metadata: str | list[str] = None
    plate_data: str | list[str] = None
    barcode: str = None

    # extracting file paths and setting into dataclass attributes
    def __post_init__(self):
        self._extract_content_files()

    def _extract_content_files(self):
        """extracts all files within given dataset folder and sets the DataFile dataclass
        attributes

        Raises
        ------
        TypeError
            raised if dataset_dir is not a str or pathlib.Path object.
            raised if plate data is not parquet or sqlite file
        """
        # get all top level files
        if not isinstance(self.dataset_dir, (str, pathlib.Path)):
            raise TypeError("dataset_dir must be a string or pathlib.Path object")
        if isinstance(self.dataset_dir, str):
            self.dataset_dir = pathlib.Path(self.dataset_dir).resolve(strict=True)

        # get all files
        all_files = list(self.dataset_dir.glob("*"))

        # get data files
        plate_data = [
            str(fpath.name)
            for fpath in all_files
            if fpath.suffix == ".parquet" or fpath.suffix == ".sqlite"
        ]
        self.plate_data = plate_data

        # get metadata_dir
        meta_data_path = [str(fpath.name) for fpath in all_files if fpath.is_dir()]
        self.metadata = (
            meta_data_path[0] if len(meta_data_path) == 1 else meta_data_path
        )

        # get barcode
        barcode_path = [
            str(fpath.name) for fpath in all_files if fpath.suffix == ".txt"
        ]
        self.barcode = barcode_path if len(barcode_path) == 1 else barcode_path


def get_test_data_folder(test_data_name: str) -> DataFiles:
    """Gets single or multiple datasets. Users provide the name of the datasets
    that will be used in their tests

    Parameters
    ----------
    test_data_name : str | list[str]
        name or names of datasets to be selected

    Returns
    -------
    DataFiles
        contains all files in a dataclass format

    Raises:
    -------
    FileNotFoundError
        Raised when the provided test_data_name is not a valid testing dataset.
    """

    # type checking
    if not isinstance(test_data_name, str):
        raise TypeError("`test_data_name` must be a string")

    # get testing dataset_path
    data_dir_path = pathlib.Path("./datasets")
    sel_test_data = (data_dir_path / test_data_name).resolve(strict=True)

    # convert to DataFiles content
    data_files = DataFiles(sel_test_data)

    return data_files


def prepare_dataset(
    test_data_name: str,
    test_dir_path: pathlib.Path,
) -> DataFiles:
    """Main function to prepare dataset into testing datafolder.

    Parameters
    ----------
    test_data_name : str
        name of the testing dataset you want to use

    test_dir_path : pathlib.Path
        path to testing directory

    Returns
    -------
    DataFiles
        Structured data object that contains all the files within the selected dataset
    """
    # get dataset and transfer to testing directory
    datafiles = get_test_data_folder(test_data_name=test_data_name)
    shutil.copytree(datafiles.dataset_dir, str(test_dir_path), dirs_exist_ok=True)

    # change directory to the testing directory
    os.chdir(str(test_dir_path))

    return datafiles


# ---------------
# PyTest Fixtures
# ---------------
@pytest.fixture
def testing_dir(tmp_path, request):
    """Creates a testing directory

    Note: Pytest will tear down tmp_path per test

    Parameters
    ----------
    tmp_path : pytest.fixture
        pytest default fixture value to be called when creating a temp dir.

    request : pytest.fixture
        Allows custom testing function to be added in the testing workflow

    returns:
    --------
    pathlib.Path
        Path pointing to temporary directory
    """
    # setting paths
    original_cwd = str(pathlib.Path(".").absolute())

    # creating a temporary path
    test_name = request.node.name
    tmp_dir = (tmp_path / test_name).absolute()
    tmp_dir.mkdir()

    # return the temporary directory path
    yield tmp_dir

    # teardown: remove the testing
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

    #  change the current working directory to the oroginal root
    os.chdir(str(original_cwd))


# ---------------
# Init Mode Tests
# ---------------
def test_multiplate_maps_no_barcode(testing_dir) -> None:
    """Negative case test: Expects `BarcodeMissingError` and Non-zero return code

    This tests checks if the CLI raises an error when a user provides multiple platemaps
    but no barcodes.

    Rational:
    ---------
        CytoSnake relies on a barcode file to accurately map each platemap to its
        corresponding plate. If users provide multiple platemaps without a barcode,
        it becomes impossible for CytoSnake to determine the correct association between
        platemaps and their respective plates.

    Parameters:
    -----------
    testing_dir: pytest.fixture
        Testing directory
    """

    # transfer data to testing folder
    prepare_dataset(test_data_name="standard_sqlite_multi", test_dir_path=testing_dir)

    # execute test
    cmd = "cytosnake init -d *.sqlite -m metadata".split()
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)

    # Grab raised exception
    raised_exception = get_raised_error(proc.stderr)

    # assert checking
    assert proc.returncode == 1
    assert raised_exception == errors.BarcodeMissingError.__name__


def test_one_plate_one_platemap(testing_dir) -> None:
    """Positive case test: Expects a 0 exit code

    This test checks if the CLI can handle one plate and one plate map without barcode
    as inputs.

    Rational:
    ---------
        As this input involves only one plate, CytoSnake does not require the process
        of identifying which platemap corresponds to that plate.

    Parameters:
    -----------
    test_dir: pytest.fixture
        Testing directory
    """

    # prepare testing files
    datafiles = prepare_dataset(
        test_data_name="standard_sqlite_single", test_dir_path=testing_dir
    )

    # Selecting one plate and meta data dir
    plate = datafiles.plate_data[0]
    metadata = datafiles.metadata

    # execute
    cmd = f"cytosnake init -d {plate} -m {metadata}".split()
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)

    # assert check
    assert proc.returncode == 0


def test_multiplates_with_multi_platemaps(testing_dir):
    """Positive case test: Expects a 0 exit code

    This test checks if the CLI can handle multiple plates and multiple plate maps
    with a barcode as inputs.

    Rational:
    ---------
        CytoSnake requires input from multiple platemaps and barcodes. To ensure
        accurate mapping of plate data, it necessitates the use of a barcode as an
        input. The barcode files are used by CytoSnake to correctly associate each
        plate map with its corresponding plate data.

    Parameters:
    -----------
    test_dir: pytest.fixture
        Testing directory
    """

    # prepare testing files
    datafiles = prepare_dataset(
        test_data_name="standard_sqlite_multi", test_dir_path=testing_dir
    )

    # Selecting one plate and meta data dir
    barcode = datafiles.barcode

    # execute
    cmd = f"cytosnake init -d *plate -m metadata -b {barcode}".split()
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)

    # assert checks
    assert proc.returncode == 0
