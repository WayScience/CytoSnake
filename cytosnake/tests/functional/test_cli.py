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

# import tempfile
# import pytest
# import subprocess

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
    pass


def test_barcode_logic_barcode_multi_platemaps() -> None:
    """Positive case: This tests expects a successful run where the user provides
    multiple plate datasets, plate map, and no barcode. Since this is only one plate_map
    , this means that the generated dataset came from one experiment and multiple
    samples (plates) were used to generated the datasets.
    """
    pass


def test_barcode_logic_no_barcode_multi_platemaps() -> None:
    """Negative case: This test expects a failed run where the user provides multiple
    plate datasets, multiple plate maps (multi-experiments), and no barcode. Since
    there are plate maps, this indicates that the generated datasets came from multiple
    experiments.
    """
    pass
