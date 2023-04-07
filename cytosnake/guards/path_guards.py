"""
Module: path_guards.py


Functions for checking paths.

This includes:
    - existing paths
    - valid path strings
"""

import pathlib
from typing import TypeGuard


def is_valid_path(val: object) -> TypeGuard[pathlib.Path]:
    """checks if provided value is a valid path"""

    # type checking
    if not isinstance(val, (str, pathlib.Path)):
        return False
    # convert to pathlib.Path
    if isinstance(val, str):
        val = pathlib.Path(val).resolve(strict=True)

    # check if the path exists
    return val.exists()
