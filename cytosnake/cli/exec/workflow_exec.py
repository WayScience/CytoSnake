"""
workflow_exec.py

Module containing functions to execute workflows via cytopipe's
CLI interface
"""
from typing import Optional
import snakemake

from cytosnake.utils.config_utils import load_workflow_path


def __base_exec(
    workflow_file: str,
    n_cores: Optional[int] = 1,
    unlock: Optional[bool] = False,
    force: Optional[bool] = False,
) -> bool:
    """Executes the cell profiler workflow

    Parameters
    ----------
    workflow_file : str
        Workflow file path
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

    # execute
    status = snakemake.snakemake(
        workflow_file, cores=n_cores, unlock=unlock, forceall=force, use_conda=True
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

    # loading configurations that contains all workflow paths
    workflow_path = load_workflow_path(wf_name=workflow)
    
    # executing selected workflow
    job = __base_exec(
        workflow_file=workflow_path,
        n_cores=n_cores,
        unlock=allow_unlock,
        force=force_run,
    )

    # checking workflow job status.
    if job is False:
        print(f"ERROR: {workflow} workflow failed")
        return 1
    return 0
