"""
module: input_guards.py

This module will handle the CytoSnake's CLI logic mostly interacting with user defined
parameters from CytoSnake's CLI.

There the logic establishes some rules of what inputs are required or what functionality
is or not allowed.
"""
import pathlib
from typing import TypeVar

from cytosnake.common.errors import BarcodeMissingError

# declaring user based type hinting
NameSpace = TypeVar("NameSpace")


def is_barcode_required(user_params: NameSpace) -> bool:
    """
    user_params: NameSpace
        Argparse.NameSpace object that contains all user provided parameters

    Returns
    -------
    bool
        With the given parameter inputs, True if barcodes are required else False
    """

    # getting both barcode and metadata from cli inputs
    barcode_param = user_params.barcode
    metadata_path = pathlib.Path(user_params.metadata).resolve(strict=True)

    # counting number of platemaps in metadata
    plate_maps_path = (metadata_path / "platemap").resolve(strict=True)
    n_platemaps = len(list(plate_maps_path.glob("*")))

    # if the metadata directory has more than 1 plate maps and no barcode file return
    # True.
    # This indicates that a barcode is required
    return n_platemaps > 1 and barcode_param is None


def check_init_parameter_inputs(user_params: NameSpace) -> bool:
    """Main wrapper to check `init` mode parameter logic.

    Parameters
    ----------
    user_params : NameSpace
        Argparse.NameSpace object that contains all user provided parameters.

    Returns
    -------
    bool
        True if all logic checks passed

    Raises
    ------
    BarcodeMissingError
        Raised if a multiple platemaps are found but no barcode file was provided
    """

    # checking if barcode is required
    if is_barcode_required(user_params=user_params):
        raise BarcodeMissingError("Barcode is required, multiple platemaps found")
