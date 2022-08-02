# ------------------------------------------------------------
# workflow_exec.py
#
# Module containing functions to execute workflows via cytopipe's
# CLI interface
# ------------------------------------------------------------
from pathlib import Path
import snakemake


def exec_preprocessing(n_cores=1) -> int:
    """Executes the cell profiler workflow

    Parameters
    ----------
    params_list : list
        list of par
    n_cores : int, optional
        max number of cores to use in the workflow, by default 1

    Returns
    -------
    bool
        True for successful execution, otherwise False
    """

    # get workflow path
    snakefile = str(Path("Snakefile").absolute())

    # execute
    status = snakemake.snakemake(snakefile, use_conda=True, cores=n_cores)
    return status
