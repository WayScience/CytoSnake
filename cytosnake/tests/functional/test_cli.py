"""
module: test_cli.py

This testing module composes of functional tests that contains checks for both positive
negative cases when using CytoSnake's CLI

A positive case indicates that given the user parameters we expect it to run
successfully.

A negative case indicates that given with the user parameters, our tests are able to
capture the errors.

Ultimately, test_cli.py will contains functional test to all modes that CytoSnake
contains.
"""

import os
import pathlib
import shutil
import subprocess
from typing import Optional

# import tempfile
# import pytest
# import subprocess
# import test_functions # This will contains helper functions for testing
from cytosnake.common import errors


# -----------------
# Helper functions
# -----------------
class CleanUpHandler:
    """Used to clean up directories in every single test run"""

    def __init__(self, tmp_path):
        self.tmp_path = tmp_path

    def __call__(self) -> None:
        shutil.rmtree(self.tmp_path)


def transfer_data(
    test_dir: pathlib.Path,
    n_plates: int,
    n_platemaps: int,
    metadata_dir_name: Optional[str] = "metadata",
    testing_data_dir="dummyfiles",
) -> None:
    """Wrapper function that transfer datasets found within the pytest module and
    transfers it to the assigned directory where pytest is conducting the functional
    tests.

    Parameters
    ----------
    test_dir : LocalPath
        `PyTest.LocalPath` object that contains the path were the test is being
        conducted

    Return
    ------
    None
        Transfers datafiles from the PyTest module to the testing directory
    """

    # get files to transfer
    dataset_dir = pathlib.Path(f"./datasets/{testing_data_dir}").resolve(strict=True)

    # grabbing all input paths
    sqlite_file_paths = list(dataset_dir.glob("*sqlite"))[:n_plates]
    platemaps_dir = dataset_dir / "metadata" / "platemap"
    plate_map_files = [
        str(_path.absolute()) for _path in platemaps_dir.glob("platemap*")
    ][:n_platemaps]
    barcode = dataset_dir / "barcode.txt"

    # create a metadata_dir in tmp_dir
    if not isinstance(metadata_dir_name, str):
        raise ValueError("metadata dir name must be a string")

    tmpdir_metadata_path = test_dir / metadata_dir_name / "platemap"
    tmpdir_metadata_path.mkdir(exist_ok=True, parents=True)

    # transferring all files to tmp dir
    for _path in plate_map_files:
        shutil.copy(_path, str(tmpdir_metadata_path))
    for _path in sqlite_file_paths:
        shutil.copy(_path, test_dir)
    shutil.copy(barcode, test_dir)


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


# --------------------------
# init mode functional tests
# --------------------------
# The tests below focuses on only executing the init mode.
def test_barcode_logic_no_barcode_one_platemap(tmp_path, request) -> None:
    """Positive case: This tests expects a successful run where the user provides
    multiple plate datasets, plate map, and no barcode. Since this is only one plate_map
    , this means that the generated dataset came from one experiment and multiple
    samples (plates) were used to generated the datasets.
    """
    # starting path
    test_module = str(pathlib.Path().absolute())

    # transfer dummy data to tmpdir
    transfer_data(test_dir=tmp_path, n_plates=2, n_platemaps=1)

    # change directory to tmpdir
    os.chdir(tmp_path)

    # execute CytoSnake
    cmd = "cytosnake init -d *.sqlite -m metadata".split()
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)

    # leave test directory
    os.chdir(test_module)

    # clean directory,
    cleanup_handler = CleanUpHandler(tmp_path)
    request.addfinalizer(cleanup_handler)

    # checking for success return code
    assert proc.returncode == 0


def test_barcode_logic_barcode_multi_platemaps(tmp_path, request) -> None:
    """Positive case: This tests expects a successful run where the user provides
    multiple plate datasets, multiple plate map, and barcode. Since this is only one
    plate_map , this means that the generated dataset came from one experiment and
    multiple samples (plates) were used to generated the datasets.
    """
    # PyTest module directory
    test_module = str(pathlib.Path().absolute())

    # transfer dummy data to tmpdir
    transfer_data(test_dir=tmp_path, n_plates=2, n_platemaps=2)

    # change directory to tmpdir
    os.chdir(tmp_path)

    # execute CytoSnake
    cmd = "cytosnake init -d *.sqlite -m metadata -b barcode.txt".split()
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)

    # leave testing dir
    os.chdir(test_module)

    # clean directory,
    cleanup_handler = CleanUpHandler(tmp_path)
    request.addfinalizer(cleanup_handler)

    # checking for success return code
    assert proc.returncode == 0


def test_barcode_logic_no_barcode_multi_platemaps(tmp_path, request) -> None:
    """Negative case: This test expects a failed run where the user provides multiple
    plate datasets, multiple plate maps (multi-experiments), and no barcode. Since
    there are plate maps, this indicates that the generated datasets came from multiple
    experiments.

    Checks:
    -------
        non-zero return code
        BarCodeRequiredError raised
    """
    # PyTest module directory
    test_module = str(pathlib.Path().absolute())

    # transfer dummy data to tmpdir
    transfer_data(test_dir=tmp_path, n_plates=2, n_platemaps=2)

    # change directory to tmpdir
    os.chdir(tmp_path)

    # execute CytoSnake
    cmd = "cytosnake init -d *.sqlite -m metadata"
    proc = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=False)
    raised_error = get_raised_error(proc.stderr)

    # leave testing dir
    os.chdir(test_module)

    # clean directory,
    cleanup_handler = CleanUpHandler(tmp_path)
    request.addfinalizer(cleanup_handler)

    # checking for success return code
    assert proc.returncode == 1
    assert raised_error == errors.BarcodeRequiredError.__name__
