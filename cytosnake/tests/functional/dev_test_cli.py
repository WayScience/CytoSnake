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
These test cases simulate scenarios where invalid or inappropriate user parameters are
provided. The objective is to confirm that the CLI can identify and handle errors
appropriately, providing informative error messages to the user.
"""

import os
from dataclasses import dataclass
import pytest
import pathlib
import shutil
import subprocess


# --------------
# helper functions
# --------------
@dataclass
class DataFiles:
    """Structured datatype representation that contains all files in a selected dataset

    Attributes:
    -----------
    metdata: str
        path to metadata directory
    plate_data: list[str]
        list of plate data (parquet or sqlite files)


    Returns
    -------
    DataFiles
        DataStructure representation of test dataset files
    """

    # required parameters
    dataset_dir: str | pathlib.Path

    # extracted files
    metadata: str = None
    plate_data: list[str] = None
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
            raise TypeError(
                "dataset_dir must be a string or pathlib.Path object"
            )
        if isinstance(self.dataset_dir, str):
            self.dataset_dir = pathlib.Path(self.dataset_dir).resolve(
                strict=True
            )

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
        meta_data_path = [
            str(fpath.name) for fpath in all_files if fpath.is_dir()
        ]
        self.metadata = (
            meta_data_path[0] if len(meta_data_path) == 0 else meta_data_path
        )

        # get barcode
        barcode_path = [
            str(fpath.name) for fpath in all_files if fpath.suffix == ".txt"
        ]
        self.barcode = (
            barcode_path[0] if len(meta_data_path) == 0 else meta_data_path
        )


def prepare_dataset(
    source_dir: pathlib.Path, test_dir_path: pathlib.Path
) -> None:
    """Transports testing datasets into testing directory

    Parameters
    ----------
    test_dir_path : pathlib.Path
        _description_

    Returns
    -------
    None
        Testing files transported to testing directory
    """
    # un
    shutil.copytree(source_dir, test_dir_path, dirs_exist_ok=True)


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
    data_dir_path = pathlib.Path("datasets")
    sel_test_data = (data_dir_path / test_data_name).resolve(strict=True)

    # convert to DataFiles content
    data_files = DataFiles(sel_test_data)

    return data_files


def get_raised_error(traceback: str) -> str:
    """Parses traceback and attempts to obtain raised exception error.

    Traceback is parsed in this order:
    1. split by new lines
    2. grab the last line as it contains the raised exception and message
    3. split by ":" to separate exception name and exception message
    4. grab the first element since it contains that path to exception
    5. split by "." and grab last element, which is the exception name

    Parameters
    ----------
    traceback : str
        complete traceback generated by executing CLI

    Returns
    -------
    str
        return raised exception error
    """

    # returns exception name, refer to function documentation to understand
    # the order of parsing the traceback to obtain exception name.
    return traceback.splitlines()[-1].split(":")[0].split(".")[-1]


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
        pytest deafault fixture value to be called when creating a temp dir.

    request : pytest.fixture
        Allows custom testing function to be added in the testing workflow

    returns:
    --------
    pathlib.Path
        Path pointing to temporary directory
    """
    # setting paths
    original_cwd = pathlib.Path()

    # creating a temporary path and cd into it
    tmp_dir = tmp_path / "testing_dir"
    tmp_dir.mkdir()

    # go that directory
    os.chdir(str(tmp_dir))

    # create a custom tear down function
    def teardown():
        # go back to original working directory and remove testing directory
        os.chdir(original_cwd)
        shutil.rmtree(tmp_dir)

    # set finalizer with custom teardown. Will call teardown() at the ending of all
    # tests
    request.addfinalizer(teardown)

    # return the temporary directory path
    return tmp_dir


# ---------------
# Init Mode Tests
# ---------------
def test_barcode_logic_no_barcode_one_platemap(testing_dir) -> None:
    """Positive case test: Expects a succesfful run.

    Test emulates a user using CytoSnake and using multiple plate

    Parameters:
    -----------
    test_dir: pytest.fixture
        Testing directory
    """

    # select dataset to use, we are using standard sqlite
    testing_dataset = get_test_data_folder(test_data_name="standard_sqlite")
    print(os.getcwd())

    # transfer data to testing folder
    prepare_dataset(
        source_dir=testing_dataset.dataset_dir, test_dir_path=testing_dir
    )

    # execute test
    cmd = "cytosnake init -d *.sqlite -m metadata".split()
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)

    # assert checking
    assert proc.returncode == 0
