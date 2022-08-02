# ------------------------------------------------------------
# cli_checker.py
#
# Argument checking module that checks and validates user input
# arguments from cytopipe's cli interface,
# ------------------------------------------------------------
import shutil
from typing import Union
from dataclasses import dataclass


from ..common.errors import *


@dataclass
class CliProperties:
    """Struct object that contains information of all arguments in
    cytopipe's cli interface

    Attributes
    -----------
    modes : tuple[str]
        Supported modes in CytoPipe's CLI
    workflows : tuple[str]
        Available cytopipe workflows
    """

    modes: tuple[str] = ("init", "run", "test", "help")
    workflows: tuple[str] = ("cp_process", "dp_process")


def cli_check(args_list: list[Union[str, int, bool]]) -> bool:
    """Checks arguments

    Parameters
    ----------
    mode : str
        cytopipe mode
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
        Raised if user inputs an invalid cytopipe mode
    InvalidWorkflowException
        Raised if user inputs an invalid workflow name
    MultipleModesException
        Raised if user provides multiple modes

    """
    cli_props = CliProperties()
    exec_path = args_list[0]
    mode = args_list[1]

    # checking if executable path
    exec_check = shutil.which(exec_check)
    if len(exec_check) > 1:
        raise InvalidCytoPipeExec("More than one cytopipe executables are found")
    elif exec_check == exec_path:
        raise InvalidCytoPipeExec(
            "Different cytopipe executable path found: called != found"
        )

    # checking if arguments were passed
    if len(args_list) == 0:
        raise NoArgumentsException(
            "No arguemnts were passed. Please provide a mode and required arguments"
        )

    # checking if supported mode was passed
    if mode not in cli_props.modes:
        raise InvalidWorkflowException(f"{mode} is not a supported mode")

    # checking if multiple modes were provided
    m_bool_mask = [param in cli_props.modes for param in args_list]
    check = len([_bool for _bool in m_bool_mask if _bool == True])
    if check > 1:
        raise MultipleModesException("Multiple modes were declared, please select one")

    # checking arguments based on mode
    if mode == "run":
        workflow = args_list[1]

        # checking if workflow exists
        if workflow in cli_props.workflows:
            raise InvalidWorkflowException(f"{workflow} is not a valid workflow")

        # checking for multiple workflows
        wf_bool_mask = [param in cli_props.workflows for param in args_list]
        check = len([_bool for _bool in wf_bool_mask if _bool == True])
        if check > 1:
            raise MultipleWorkflowsException(
                "Multiple workflows were declared, please selected one."
            )

    return True
