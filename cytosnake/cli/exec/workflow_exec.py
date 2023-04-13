"""
workflow_exec.py

Module containing functions to execute workflows via cytopipe's
CLI interface
"""
from typing import Optional

import snakemake
from snakemake import shell

from cytosnake.utils.config_utils import load_general_configs


def __base_exec(
    workflow_file: str,
    n_cores: Optional[int] = 1,
    unlock: Optional[bool] = False,
    force: Optional[bool] = False,
    env_manager: Optional[str] = "conda",
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

    Raises
    ------
    TypeError
        Raised if env_manager is not a string type
    """
    # type checking
    if not isinstance(env_manager, str):
        raise TypeError("`env_manager` must be a string type")
    if env_manager not in ["conda", "mamba"]:
        raise ValueError(
            f"{env_manager} is invalid env_manager"
            "Only `conda` or `mamba` is supported."
        )

    # execute workflow with given configs
    return snakemake.snakemake(
        workflow_file,
        cores=n_cores,
        unlock=unlock,
        forceall=force,
        use_conda=True,
        conda_prefix=env_manager,
        conda_frontend=env_manager,
    )


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

    # load in general configs
    cytosnake_configs = load_general_configs()

    # executing selected workflow
    job = __base_exec(
        workflow_file=workflow,
        n_cores=n_cores,
        unlock=allow_unlock,
        force=force_run,
        env_manager=cytosnake_configs["env_manager"],
    )

    # checking workflow job status.
    if job is False:
        print(f"ERROR: {workflow} workflow failed")
        return 1
    return 0
