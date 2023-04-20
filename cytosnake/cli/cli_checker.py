"""
Documentation
-------------

cli_checker.py Module

Argument checking module that checks and validates user input
arguments from CytoSnake's cli interface,
"""
import sys
from dataclasses import dataclass
from typing import Union

from cytosnake.cli.cli_docs import cli_docs, init_doc, run_doc
from cytosnake.common.errors import (
    InvalidWorkflowException,
    MultipleModesException,
    MultipleWorkflowsException,
    NoArgumentsException,
)


@dataclass
class CliProperties:
    """Struct object that contains information of all arguments in
    CytoSnake's cli interface

    Attributes
    -----------
    modes : tuple[str]
        Supported modes in CytoSnake's CLI
    workflows : tuple[str]
        Available CytoSnake workflows
    """

    modes: tuple[str] = ("init", "run", "test", "help")
    workflows: tuple[str] = ("cp_process", "dp_process")


def cli_check(args_list: list[Union[str, int, bool]]) -> bool:
    """Checks arguments

    Parameters
    ----------
    mode : str
        CytoSnake mode
    args_list : list[Union[str, int, bool]]
        list of user inputted arguments

    Returns
    -------
    None
        Passed the checking system, Raises errors if invalid
        arguments are passed.

    Raises
    ------
    NoArgumentsException
        Raised if the user provides no arguments
    InvalidModeException
        Raised if user inputs an invalid CytoSnake mode
    InvalidWorkflowException
        Raised if user inputs an invalid workflow name
    MultipleModesException
        Raised if user provides multiple modes

    """
    cli_props = CliProperties()

    # checking if arguments were passed
    if len(args_list) == 1:
        raise NoArgumentsException(
            "No Arguments were passed. Please provide a mode and required arguments"
        )

    # checking if supported mode was passed
    mode = args_list[1]
    if mode not in cli_props.modes:
        raise InvalidWorkflowException(f"{mode} is not a supported mode")

    # # checking if help message was called
    # _check_mode_help_arg(args_list)

    # checking if multiple modes were provided
    m_bool_mask = [
        param.lower() in cli_props.modes
        for param in args_list
        if param.lower() != "help"
    ]
    check = len([_bool for _bool in m_bool_mask if _bool])
    if check > 1:
        raise MultipleModesException("Multiple modes were declared, please select one")

    # checking arguments based on mode
    if mode == "run":
        workflow = args_list[2].lower()

        if workflow == "help":
            _check_mode_help_arg(args_list)

        # checking if workflow exists
        if workflow not in cli_props.workflows:
            raise InvalidWorkflowException(f"{workflow} is not a valid workflow")

        # checking for multiple workflows
        wf_bool_mask = [param in cli_props.workflows for param in args_list]
        check = len([_bool for _bool in wf_bool_mask if _bool])
        if check > 1:
            raise MultipleWorkflowsException(
                "Multiple workflows were declared, please selected one."
            )

    return True


def _check_mode_help_arg(args_list: list[Union[str, int, bool]]) -> None:
    """Checks if help documentation is called in the arguments

    Parameters
    ----------
    args_list : list[Union[str, int, bool]]
        list of user provided arguments
    """
    mode_opt = args_list[0].lower()

    # display cli documentation and exit
    if mode_opt == "help":
        print(cli_docs)
        sys.exit(0)

    # display mode help documentation and exit
    if len(args_list) == 2:
        mode_help = args_list[1]
        if mode_help == "help":
            match mode_opt.lower():
                case "init":
                    print(init_doc)
                    sys.exit(0)
                case "run":
                    print(run_doc)
                    sys.exit(0)
                case "test":
                    raise NotImplementedError("Documentation is not implemented yet")
                case _:
                    raise RuntimeError("Unexpected Error")
    else:
        return
