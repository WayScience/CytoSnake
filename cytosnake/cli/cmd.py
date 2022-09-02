"""
Documentation

cmd.py Module

Generates CLI interface in order to interact with CytoSnake.
"""
import sys

from .args import CliControlPanel
from .exec.workflow_exec import workflow_executor
from .setup_init import init_cp_data, init_dp_data
from .cli_docs import init_doc, cli_docs, run_doc
from ..common.errors import WorkflowFailedException


def run_cmd() -> None:
    """obtains all parameters and executes workflows

    Parameters
    ----------
    params : list
        list of user provided parametersCytoSnake
    Returns
    -------
    None
    """
    # create args handler
    # -- Cli Control Panel
    args_handler = CliControlPanel(sys.argv)
    # print(args_handler.__repr__())

    # checking is user wanted to cli help
    if args_handler.cli_help is True:
        print(cli_docs)
        sys.exit(0)

    match args_handler.mode:

        case "init":

            # if mode help flag is shown, show cli help doc
            if args_handler.mode_help is True:
                print(init_doc)
                sys.exit(0)

            # setting up input files for cytosnake
            print("INFO: Formatting input files")
            init_args = args_handler.parse_init_args()

            # identifying which data type was added and how to set it up
            match init_args.datatype:
                case "cell_profiler":
                    init_cp_data(
                        data_fp=init_args.data,
                        metadata_fp=init_args.metadata,
                        barcode_fp=init_args.barcode,
                    )
                case "deep_profiler":
                    init_dp_data(data_fp=init_args.data, metadata_fp=init_args.metadata)
                case _:
                    raise RuntimeError(
                        "Unexpected error in identifying datatype. Did you specify `cell_profiler` or `deep_profiler`datatype?"
                    )

            print("INFO: Formatting complete!")

        case "run":
            # display run help documentation
            if args_handler.mode_help is True:
                print(run_doc)
                sys.exit(0)

            # parsing workflow parameters
            print(f"INFO: Executing {args_handler.workflow} workflow")
            wf_params = args_handler.parse_workflow_args()

            wf_executor = workflow_executor(
                workflow=wf_params.workflow,
                n_cores=wf_params.max_cores,
                allow_unlock=wf_params.lock,
                force_run=wf_params.force,
            )

            if wf_executor != 0:
                raise WorkflowFailedException(
                    f"Workflow encounter and error, please refer to the logs"
                )

        case "help":
            print(cli_docs)

        case _:
            raise RuntimeError("Unexpected error captured in mode selection")


if __name__ == "__main__":

    run_cmd()
