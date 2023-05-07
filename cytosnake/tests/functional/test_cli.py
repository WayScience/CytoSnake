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

import pathlib
from dataclasses import dataclass
from typing import Optional

# import tempfile
# import pytest
# import subprocess
# import test_functions # This will contains helper functions for testing


# -----------------
# Helper functions
# -----------------
@dataclass
class DummyInputs:
    platedata: str
    metadata: str
    barcode: Optional[str] = None


def generate_input_paths(
    n_plates: int, n_platemaps: int, barcode: Optional[bool] = True
) -> DummyInputs:
    """Generates DummyInputs object that contains platedata, barcode and metadata paths.
    These paths point to dummy

    Parameters
    ----------
    n_plates : int
        number of dummy sqlite files
    n_platemaps : int
        number of dummy plate map files
    barcode : Optional[bool], optional
        True if a dummy barcode should be added into the DummyInputs, Else false
        , by default True
    """

    # getting root dir where all the dummy data is stored
    root_dummy_data_dir = pathlib.Path("./datasets/dummyfiles").resolve(strict=True)

    # now generated the amount of dummy data based on parameters
    # getting plate data
    plates_list = list(root_dummy_data_dir.glob("*.sqlite"))[: n_plates + 1]
    plates_list_str = " ".join(str(_path.absolute()) for _path in plates_list)

    # getting plate maps
    platemaps_dir = root_dummy_data_dir / "metadata" / "platemap"
    platemaps_list = list(platemaps_dir.glob("platemap*"))[: n_platemaps + 1]
    platemaps_list_str = " ".join(str(_path.absolute()) for _path in platemaps_list)

    # getting barcode
    if barcode is True:
        barcode_path = str(root_dummy_data_dir / "barcode.txt")
    else:
        barcode_path = None

    return DummyInputs(plates_list_str, platemaps_list_str, barcode_path)


# --------------------------
# init mode functional tests
# --------------------------
# The tests below focuses on only executing the init mode.


def test_barcode_logic_no_barcode_one_platemap() -> None:
    """Positive case: This tests expects a successful run where the user provides
    multiple plate datasets, plate map, and no barcode. Since this is only one plate_map
    , this means that the generated dataset came from one experiment and multiple
    samples (plates) were used to generated the datasets.
    """
    # generating input strings
    _inputs = generate_input_paths(n_plates=3, n_platemaps=1, barcode=False)


def test_barcode_logic_barcode_multi_platemaps() -> None:
    """Positive case: This tests expects a successful run where the user provides
    multiple plate datasets, plate map, and no barcode. Since this is only one plate_map
    , this means that the generated dataset came from one experiment and multiple
    samples (plates) were used to generated the datasets.
    """
    # generating input strings
    _inputs = generate_input_paths(n_plates=3, n_platemaps=2, barcode=True)

    # execute cli


def test_barcode_logic_no_barcode_multi_platemaps() -> None:
    """Negative case: This test expects a failed run where the user provides multiple
    plate datasets, multiple plate maps (multi-experiments), and no barcode. Since
    there are plate maps, this indicates that the generated datasets came from multiple
    experiments.
    """
    # generating input strings
    _inputs = generate_input_paths(n_plates=3, n_platemaps=2, barcode=False)

    # execute cli
