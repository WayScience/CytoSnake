"""
Module: file_utils.py

Contains functions that involves file manipulations like:
    - searching 
    - transferring 
    - creating symbolic links 
    - creating and deleting files
"""
from pathlib import Path
from typing import Optional
from cytosnake.guards.path_guards import is_valid_path


def file_search(fpath: str | Path) -> dict:
    """Returns a list of all the files inside a directory

    Parameters
    ----------
    fpath : str | Path
        path to specific directory

    Returns
    -------
    dict
        Returns a list of dir names and file paths as key value pairs

    Raises
    ------
    """
    # type checking
    # -- convert str to Path object
    if isinstance(fpath, str):
        fpath = Path(fpath).absolute()

    # -- check if the path is valid
    if not is_valid_path(fpath):
        raise FileNotFoundError("Invalid path was provided")

    # search with given path
    list_of_files = fpath.glob("*")

    # collected_files = defaultdict(list)
    collected_files = {}
    for _file in list_of_files:
        name = _file.stem

        # if the path is a directory, get the files within directory (recursive)
        if _file.is_dir():
            nested_files = file_search(_file)
            collected_files[name] = nested_files
        elif _file.is_file():
            # collected_files[name].append(str(_file.absolute()))
            collected_files[name] = str(_file.absolute())
        else:
            raise RuntimeError(
                "Unexpected error captured: Non-path entry captured "
            )

    return collected_files


def find_project_dir(steps: Optional[int] = 10) -> Path:
    """Recursively searches for `.cytosnake` directory

    steps: Optional[int]
        number of recursive steps for searching. [default=10]

    Returns
    -------
    Path
        returns path that contains the `.cytosnake`

    Raises
    ------
    FileNotFoundError
        When recursive steps are exceeded and is unable to find `.cytosnake`
        folder
    """

    # type checking
    if not isinstance(steps, int):
        raise TypeError(f"'steps' must be an integer type. Not {type(steps)}")

    # using current working directory as starting pointg
    start_point = Path().absolute()
    for step in range(steps):
        # if the number of recursive steps is exceeded, raise error
        if step > 10:
            raise FileNotFoundError("unable to find project folder")

        # grab all files in directory
        all_files = start_point.glob("*")
        for _file in all_files:

            # check if the file is a directory and has the name `cytosnake`
            # -- if true, return the complete path
            if _file.is_dir() and _file.name == ".cytosnake":
                return _file.absolute()

        start_point = start_point.parent
