"""
Module: path_guards.py


Functions for checking paths.

This includes:
    - existing paths
    - valid path strings
"""

from typing import TypeGuard
from pathlib import Path


def is_valid_path(val: object) -> TypeGuard[Path]:
    """checks if provided value is a valid path"""

    # check if the val is valid type
    # -- if string, convert to Path
    accepted_types = (str, Path)
    if not isinstance(val, accepted_types):
        return False
    if isinstance(val, str):
        val = Path(val)

    # check if the path exists
    return bool(val.exists())
