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

import pytest

from cytosnake.common import errors
from cytosnake.tests import test_utils


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
@pytest.mark.negative
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
    test_utils.prepare_dataset(
        test_data_name="standard_sqlite_multi", test_dir_path=testing_dir
    )

    # execute test
    proc = subprocess.run(
        "cytosnake init -d *.sqlite -m metadata".split(),
        capture_output=True,
        text=True,
        check=False,
    )

    # Grab raised exception
    raised_exception = test_utils.get_raised_error(proc.stderr)

    # assert checking
    assert proc.returncode == 1
    assert raised_exception == errors.BarcodeMissingError.__name__


@pytest.mark.postive
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
    datafiles = test_utils.prepare_dataset(
        test_data_name="standard_sqlite_single", test_dir_path=testing_dir
    )

    # Selecting one plate and meta data dir
    plate = datafiles.plate_data[0]
    metadata = datafiles.metadata

    # execute
    proc = subprocess.run(
        f"cytosnake init -d {plate} -m {metadata}".split(),
        capture_output=True,
        text=True,
        check=False,
    )

    # assert check
    assert proc.returncode == 0


@pytest.mark.postive
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
    datafiles = test_utils.prepare_dataset(
        test_data_name="standard_sqlite_multi", test_dir_path=testing_dir
    )

    # Selecting one plate and meta data dir
    barcode = datafiles.barcode

    # execute
    proc = subprocess.run(
        f"cytosnake init -d *plate -m metadata -b {barcode}".split(),
        capture_output=True,
        text=True,
        check=False,
    )

    # assert checks
    assert proc.returncode == 0
