"""
Documentation

cmd.py Module

Generates CLI interface in order to interact with CytoSnake.
"""
import sys
import logging

# cytosnake imports
from cytosnake.cli.args import CliControlPanel
from cytosnake.cli.cli_docs import cli_docs, init_doc, run_doc
from cytosnake.cli.exec.workflow_exec import workflow_executor
from cytosnake.cli.setup_init import init_cp_data, init_dp_data
from cytosnake.common.errors import ProjectExistsError, WorkflowFailedException
from cytosnake.utils.cyto_paths import is_cytosnake_dir
from cytosnake.utils.cytosnake_setup import setup_cytosnake_env


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

    # set logging configurations
    logging.basicConfig(level="DEBUG")

    # create args handler
    # -- Cli Control Panel
    args_handler = CliControlPanel(sys.argv)
    # print(args_handler.__repr__())

    # checking is user wanted to cli help
    if args_handler.cli_help is True:
        print(cli_docs)
        sys.exit(0)

    # Main mode selection function. Each match
    match args_handler.mode:
        case "init":
            # if mode help flag is shown, show cli help doc
            if args_handler.mode_help is True:
                print(init_doc)
                sys.exit(0)

            # checking if current directory is a project folder
            # -- if True, raise error
            if is_cytosnake_dir():
                raise ProjectExistsError(
                    "This directory is already a cytosnake project directory"
                )

            # setting up input files for cytosnake
            logging.info(msg="Formatting input files")
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
                        "Unsupported datatype. Supported datatypes:"
                        "`cell_profiler` or `deep_profiler`datatype?"
                    )

            # now that the data is created, set up the current directory
            # into a project directory
            setup_cytosnake_env(args=init_args)
            logging.info("Initialization complete")

        # Executed if the user is using the `run` mode. This will execute the
        # workflow that are found within the `workflows` folder
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

            # exception error if the workflow failed
            if wf_executor != 0:
                raise WorkflowFailedException(
                    "Workflow encounter and error, please refer to the logs"
                )

        # if user execute `help` help mode. CLI documentation will appear.
        case "help":
            print(cli_docs)

        # Raise an error if an invalid mode is provided
        case _:
            raise RuntimeError("Unexpected error captured in mode selection")


if __name__ == "__main__":

    run_cmd()
