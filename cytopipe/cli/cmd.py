"""
Documentation

cmd.py Module

Generates CLI interface in order to interact with cytopipe.
"""
import sys
import shutil
from pathlib import Path

from .args import CliControlPanel
from .exec.workflow_exec import workflow_executor
from .cli_docs import init_doc, cli_docs, run_doc
from ..common.errors import WorkflowFailedException

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

            # setting up input files for cytopipe
            print("INFO: Formatting input files")
            init_args = args_handler.parse_init_args()

            # setting up file paths
            barcode_path = str(Path(init_args.barcode).absolute())
            metadata_path = str(Path(init_args.metadata).absolute())

            # create data folder on working directory
            # raises error if data directory exists (prevents overwriting)
            data_dir_obj = Path("data")
            data_dir_obj.mkdir(exist_ok=False)
            data_dir_path = str(data_dir_obj.absolute())

            # moving all user provided plate data files into data folder
            for data_file in init_args.data:
                f_path = str(Path(data_file).absolute())
                shutil.move(f_path, data_dir_path)

            # user provided barcode file and metadata directory is
            # moved to the data folder
            shutil.move(barcode_path, data_dir_path)
            shutil.move(metadata_path, data_dir_path)

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
                use_conda_env=wf_params.conda_env,
                allow_unlock=wf_params.lock,
                force_run=wf_params.force,
            )
            if not wf_executor:
                raise WorkflowFailedException(
                    f"Workflow encounter and error, please refer to the logs"
                )

        case "help":
            print(cli_docs)

        case _:
            raise RuntimeError("Unexpected error captured in mode selection")


if __name__ == "__main__":

    run_cmd()
