"""
workflow_exec.py

Module containing functions to execute workflows via cytopipe's
CLI interface
"""
from typing import Optional
from pathlib import Path
import snakemake


def __base_exec(
    n_cores: Optional[int] = 1,
    unlock: Optional[bool] = False,
    force: Optional[bool] = False,
) -> int:
    """Executes the cell profiler workflow

    Parameters
    ----------
    params_list : list
        list of par
    n_cores : int, optional
        max number of cores to use in the workflow, by default 1
    unlock : bool, optional
        if true, allows unlocking the directory where the workflow is being executed
        , by default False
    force : bool, optional
        if True, recreates all files produces from workflow


    Returns
    -------
    bool
        True for successful execution, otherwise False
    """

    # get workflow path
    snakefile = str(Path("Snakefile").absolute())

    # execute
    status = snakemake.snakemake(
        snakefile,
        cores=n_cores,
        unlock=unlock,
        forceall=force,
    )
    return status


def workflow_executor(
    workflow,
    n_cores: Optional[int] = 1,
    allow_unlock: Optional[bool] = False,
    force_run: Optional[bool] = False,
) -> int:
    """Wrapper for executing cytopipe workflows

    Parameters
    ----------
    n_cores : int, optional
        max number of cores to use in the workflow, by default 1
    allow_unlock : bool, optional
        Locks working directory when running cytopipe if
        interrupted. If set to True, cytopipe will automatically
        unlock the working directory. by default False
    force_run : bool, optional
        If set to True, when re-running cytopipe's workflow will
        start from the beginning. If False, the workflow will
        start where it last stopped. by default False
    """
    job = __base_exec(
        n_cores=n_cores, unlock=allow_unlock, force=force_run
    )
    if job is False:
        print(f"WARNING: {workflow} workflow failed")
        return 1
    return 0
