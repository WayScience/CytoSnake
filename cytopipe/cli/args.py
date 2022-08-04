# ------------------------------------------------------------
# args.py
#
# Responsible for centralizing and handling all user inputted
# parameters for cytopipe.
#
# Contains different modes with their respected parameters
# ------------------------------------------------------------
import argparse


def parse_init_args(args_list: list) -> argparse.Namespace:
    """Parses user inputs for setting up files

    Parameters
    ----------
    params : list
        User defined input parameters

    Returns
    -------
    argparse.Namespace
        parsed cli parameters

    Attributes:
    -----------
    mode : str
        name of the mode used
    data : list[str]
        list of plate data paths
    metadata : str
        Path to meta data directory
    barcode : str
        Path to barcode file
    platemaps : str
        Path to platemap file

    """

    # CLI inputs
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "mode",
        choices=["init"],
        help="Init mode sets up the all files for cytopipe processing",
    )
    parser.add_argument(
        "-d", "--data", nargs="+", required=True, help="list of plate data files"
    )
    parser.add_argument(
        "-m", "--metadata", type=str, required=True, help="path to metadata directory"
    )
    parser.add_argument(
        "-b", "--barcode", type=str, required=True, help="path to barcodes file"
    )
    args = parser.parse_args(args_list)

    return args


def parse_cp_process_args(args_list: list) -> argparse.Namespace:
    """Parses user inputs for cell profiler processing

    Parameters
    ----------
    args_list : list
        User defined input parameters

    Returns
    -------
    argparse.Namespace
        parsed cli parameters.

    Attributes
    ----------
    mode : str
        returns name of the mode
    workflow : str
        returns path of the workflow
    max_cores : int
        maximum number of cores that will be used

    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "mode", choices=["run"], type=str, help="Run mode executes specific workflows"
    )
    parser.add_argument("workflow", type=str, help="Name of desired workflow")
    parser.add_argument(
        "-c",
        "--max_cores",
        type=int,
        default=1,
        help="maximum number of cores to run the workflow",
    )
    args = parser.parse_args(args_list)

    return args
