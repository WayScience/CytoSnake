"""
Documentation

args.py Module

Responsible for centralizing and handling all user inputted
parameters for cytosnake.

Contains different modes with their respected parameters
"""

import argparse
import shutil
from dataclasses import dataclass
from typing import Union

# cytosnake imports
from cytosnake.common.errors import (
    InvalidArgumentException,
    InvalidExecutableException,
    InvalidWorkflowException,
)
from cytosnake.utils.config_utils import load_workflow_path, load_workflow_paths_config


# CLI helper functions
def supported_workflows() -> tuple[str]:
    """Returns a tuple of supported workflows that `run` supports.

    Returns
    -------
    tuple[str]
        tuple containing the names of the workflow names
    """

    # load in workflow path config
    # -- only grab keys and add "help" into the supported workflows
    workflow_paths = list(load_workflow_paths_config().keys()) + ["help"]

    # returns as a tuple, we do not want to edit these names
    return tuple(workflow_paths)


# custom argparse actions
class WorkflowSearchPath(argparse.Action):
    """This class provides more functionality

    User provided workflow names are checked to see if they exist. If found,
    the parameter will return the absolute path to the workflow. If the name
    does not exist, raise an error indicating that the provided
    workflow name does not exist within cytosnake
    """

    def __call__(self, parser, args, values, option_string=None):

        # checking if user provided workflow exists
        supported_wf = supported_workflows()
        if values not in supported_wf:
            raise InvalidWorkflowException(
                f"Unable to find '{values}'. Please specify a supported workflow: {supported_wf}"
            )
        # grabbing and setting the new value with the extracted path
        values = str(load_workflow_path(values))

        # return new attributes of the `workflow` parameter
        setattr(args, self.dest, values)


# CLI Control Panel, controls all user inputs
@dataclass(repr=True)
class CliControlPanel:
    """

    Returns
    -------
    class CliControlPanel
        CliControlPanel Class contains all the parameters stored
        in the instances.
        The instance contains states of the provided parameters a user
        provides and provides easy logic flow for the CLI to work with.

    Raises
    ------
    InvalidExecutableException
        Raised if mismatching executable paths do not match
    InvalidArgumentException
        Raised if invalid arguments are provided
    InvalidModeException
        Raised if an invalid mode is provided
    """

    param_list: list[str]
    exec_path: Union[None, str] = None
    mode: Union[None, str] = None
    data_type: Union[None, str] = None
    workflow: Union[None, str] = None
    mode_check = False
    cli_help = False
    mode_help = False
    data_type: tuple[str] = ("cell_profiler", "deep_profiler")
    modes: tuple[str] = ("init", "run", "test", "help")
    __exec_name: str = "cytosnake"

    # ----------------------------------------
    # Metaclass class methods
    # ----------------------------------------

    # Used for setting up user parameters
    def __post_init__(self):
        self.__setup_args()

    # Class representation
    def __repr__(self) -> str:
        exec_mode_workflow = (
            f"exec_path={self.exec_path}, mode={self.mode}, data_type={self.data_type}"
        )
        help_checks = f"cli_help={self.cli_help}, mode_help={self.mode_help}"
        return f"CliArgsHandler({exec_mode_workflow}, {help_checks}"

    # making class printable
    def __str__(self) -> str:
        return self.__repr__()

    # ----------------------------------------
    # class in methods
    # ----------------------------------------
    def parse_init_args(self) -> argparse.Namespace:
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
        params = self.param_list[1:]
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "mode",
            choices=["init"],
            help="Init mode sets up the all files for cytosnake processing",
        )
        parser.add_argument(
            "-d",
            "--data",
            nargs="+",
            required=True,
            help="list of plate data files",
        )
        parser.add_argument(
            "-m",
            "--metadata",
            type=str,
            help="path to metadata directory",
            required=True,
        )
        parser.add_argument(
            "-b",
            "--barcode",
            type=str,
            required=False,
            default=None,
            help="path to barcodes file",
        )
        parser.add_argument(
            "--datatype",
            choices=list(self.data_type),
            default="cell_profiler",
            help="datatype used for processing",
        )
        args = parser.parse_args(params)

        return args

    def parse_workflow_args(self) -> argparse.Namespace:
        """Wrapper for parsing workflow parameters

        Returns
        -------
        argparse.Namespace
            returns parsed parameters based on workflow
        """
        args = self.__workflow_args_parser()
        return args

    # ----------------------------------------
    # Private formatting functions
    # ----------------------------------------

    def __setup_args(self):
        """Checks the input parameter list and updates parameter states


        Raises
        ------
        InvalidExecutableException
            Raised when issues in regards to the executable paths
        """

        # ------------------------------
        # check if path executing is the same as stored in env
        # ------------------------------
        check = shutil.which(self.__exec_name)
        if isinstance(check, str):
            if self.param_list[0] != check:
                raise InvalidExecutableException("Executable mismatches")
        if check is None:
            raise InvalidExecutableException("Unable to find executable")
        self.exec_path = self.param_list[0]

        # ------------------------------
        # checking the modes parameters
        # ------------------------------
        if not len(self.param_list) > 1:
            raise InvalidArgumentException(
                "No mode has been provided, please enter a mode"
            )

        mode = self.param_list[1]
        if mode in self.modes:
            self.mode_check = True

        if mode == "help":
            self.cli_help = True

            # checking any arguments after help
            self.__check_extra_help_args(help_flag_pos=1)

        # setting mode
        self.mode = mode

        # ------------------------------
        # checking the workflow parameters
        # ------------------------------
        # checking for mode help
        try:
            mode_help = self.param_list[2]
            if mode == "init" and mode_help == "help":
                self.__check_extra_help_args(help_flag_pos=2)
                self.mode_help = True

            elif mode == "run" and mode_help == "help":
                self.__check_extra_help_args(help_flag_pos=2)
                self.mode_help = True

        # indicates that the user only placed a mode
        except IndexError:
            self.mode_help = False
            self.workflow = False

    def __workflow_args_parser(self) -> argparse.Namespace:
        """Parses user inputs for CellProfiler processing.`

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

        # converting inputs as a list
        params = self.param_list[1:]
        parser = argparse.ArgumentParser()
        required = parser.add_argument_group("Required Arguments")
        cli_configs = parser.add_argument_group("Config Arguments")
        required.add_argument(
            "mode",
            choices=["run"],
            type=str,
            help="Run mode executes specific workflows",
        )
        required.add_argument(
            "workflow",
            action=WorkflowSearchPath,
            type=str,
            help="Name of desired workflow",
        )
        cli_configs.add_argument(
            "-c",
            "--max_cores",
            type=int,
            default=1,
            required=False,
            help="maximum number of cores to run the workflow",
        )
        cli_configs.add_argument(
            "--lock",
            action="store_true",
            default=False,
            required=False,
            help="Locking mechanism",
        )
        cli_configs.add_argument(
            "--force",
            action="store_true",
            default=False,
            required=False,
            help="Force run the entire workflow",
        )
        args = parser.parse_args(params)

        return args

    def __check_extra_help_args(self, help_flag_pos):
        """Checks arguemnts are help flag

        help_flag_pos : int
            index positions where the help flag is placed

        Returns
        -------
        None
            raises InvalidArgumentException if additional parameters are
            added after

        Raises
        ------
        IndexError
            Raised if the indicated help flag position is out of bounds
        InvalidArgumentException
            Raised if additional arguments are added after the help flag
        """
        extra_args = self.param_list[help_flag_pos + 1 :]
        if len(self.param_list) < help_flag_pos:
            raise IndexError(
                "Size of parameter inputs is smaller than the position of the help flag"
            )
        if len(extra_args) != 0:
            raise InvalidArgumentException(
                f"Unknown help parameters provided: {extra_args}"
            )
