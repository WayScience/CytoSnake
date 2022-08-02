# ------------------------------------------------------------
# cmd.py
#
# Generates CLI interfaces in order to interact with cytopipe.
# ------------------------------------------------------------
import sys
import shutil
from pathlib import Path

from .args import *
from .args_checker import CytoArgs
from ..exec.workflow_exec import exec_preprocessing


def run_cmd() -> None:
    """obtains all parameters and executes workflows

    Parameters
    ----------
    params : list
        list of user provided parameters

    Returns
    -------
    None
    """
    # collecting parameters
    cyto_args = CytoArgs()
    params = sys.argv[1:]
    if len(params) == 0:
        raise RuntimeError("please provide a mode")

    mode_type = params[0]

    # checking if mode type is supported
    # supported_modes = ["init", "run", "test"]
    if mode_type not in cyto_args.modes:
        raise ValueError(
            f"{mode_type} is not a mode supported modes: {cyto_args.modes}"
        )

    # checking if more than one mode is provided, raise error
    bool_mask = [param in cyto_args.modes for param in params]
    check = len([_bool for _bool in bool_mask if _bool == True])
    if check > 1:
        raise RuntimeError("two modes detected, please select one")

    # parsing arguments based on modes
    if mode_type == "init":

        # parsing inputs
        init_args = parse_init_args(params)

        # setting up file paths
        platemap_path = str(Path(init_args.platemap).absolute())
        barcode_path = str(Path(init_args.barcode).absolute())
        metadata_path = str(Path(init_args.metadata).absolute())

        # create data folder
        # raises error if data directory exists (prevents overwriting)
        data_dir_obj = Path("data")
        data_dir_obj.mkdir(exist_ok=False)
        data_dir_path = str(data_dir_obj.absolute())

        # moving all input files
        for data_file in init_args.data:
            f_path = str(Path(data_file).absolute())
            shutil.move(f_path, data_dir_path)

        shutil.move(platemap_path, data_dir_path)
        shutil.move(barcode_path, data_dir_path)
        shutil.move(metadata_path, data_dir_path)

    elif mode_type == "run":

        # selecting workflow process
        proc_sel = params[1]

        # checking if the workflow process exists
        if proc_sel not in cyto_args.workflows:
            raise ValueError(f"{proc_sel} is not a supported workflow")

        # executing process
        if proc_sel == "cp_process":
            cp_process_args = parse_cp_process_args(params)

            print("executing cell profiler preprocessing")
            status = exec_preprocessing(cores=cp_process_args.max_cores)

            # CLI exit based on status
            if status:
                sys.exit(0)
            else:
                print("Error: Unsuccessful workflow")
                sys.exit(1)

    elif mode_type == "help":
        # display documentation of all modes
        pass


if __name__ == "__main__":

    run_cmd()
