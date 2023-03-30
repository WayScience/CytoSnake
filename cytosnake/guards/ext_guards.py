"""
module: ext_guards.py

Checks if the correct extensions are provided
"""

import pathlib
from cytosnake.guards.path_guards import is_valid_path
from typing import TypeGuard


def has_parquet_ext(file_name: str | pathlib.Path) -> TypeGuard[str]:
    """Checks if the provided path is a parquet file.

    Parameters
    ----------
    file_name : str | pathlib.Path
        path to file

    Returns
    -------
    TypeGuard[str]
        return True if it is a parquet file, else False
    """
    return (
        file_name.suffix in [".parquet", ".parq", ".pq"]
        if is_valid_path(file_name)
        else False
    )


def has_sqlite_ext(file_name: str | pathlib.Path) -> TypeGuard[str]:
    """Checks if the provided path is a parquet file.

    Parameters
    ----------
    file_name : str | pathlib.Path
        path to file

    Returns
    -------
    TypeGuard[str]
        return True if it is a parquet file, else False
    """

    return (
        file_name.suffix in [".sqlite", ".sqlite3"]
        if is_valid_path(file_name)
        else False
    )
