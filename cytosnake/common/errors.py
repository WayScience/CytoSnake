# ------------------------------------------------------------
# errors.py
#
# Modules containing CytoPipe specific exceptions
#
# Base Exceptions:
# Focuses on basic error types like incorrect values or invalid
# operations.
#
# Cytopipe Specific Errors:
# Exceptions that specific to Cytopipe
# ------------------------------------------------------------

import sys
from typing import Union, Optional

# ------------------------------
# Base Exceptions
# ------------------------------
class BaseValueError(ValueError):
    """Base exception if incorrect values are passed"""


class BaseFileNotFound(FileNotFoundError):
    """Raised if a requested file is not found in cytopipe"""


class BaseFileExistsError(FileExistsError):
    """Raised if an existing file or directory is found in a directory"""


class BaseExecutorException(RuntimeError):
    """Base exception related to cytopipe execution"""


class BaseWorkflowException(RuntimeError):
    """Base exception related to cytopipe's workflow errors in runtime"""


# ------------------------------
# Cytopipe specific errors
# ------------------------------
class InvalidArgumentException(BaseValueError):
    """Raised when arguments requirements are not met"""


class InvalidCytoPipeExec(BaseExecutorException):
    """Raised if invalid cytopipe executable is being called"""


class InvalidModeException(BaseValueError):
    """Raised if in unsupported mode was passed"""


class InvalidWorkflowException(BaseValueError):
    """Raised if invalid workflows were specified"""


class MultipleModesException(BaseValueError):
    """Raised if multiple modes were provided"""


class NoArgumentsException(BaseValueError):
    """Raised if no arguments were passed in the CLI interface"""


class MultipleWorkflowsException(BaseValueError):
    """Raised if multiple workflows are declared"""


class InvalidExecutableException(BaseExecutorException):
    """Raised if any executables"""


class WorkflowNotFoundError(BaseFileNotFound):
    """Raised if workflow file is not found"""


class WorkflowFailedException(BaseWorkflowException):
    """Raised if a workflow fails during runtime"""


class ProjectExistsError(BaseFileExistsError):
    """Raised when `.cytosnake` file is found within a directory indicating
    that the current directory has already been set up for cytosnake analysis"""


# -----------------------
# Error handling functions
# -----------------------
def display_error(
    error: Union[BaseException, Exception],
    e_msg: Optional[Union[None, str]] = None,
    exit_code: Optional[int] = 1,
) -> None:
    """This function takes in an error type and a custom message. If the custom
    message is None, then the default string data found within the exception
    type will be used.
    The function will print the message of the error and ensures a non-zero
    exit code.

    Parameters
    ----------
    error : Union[BaseException, Exception]
        Takes in a python error object
    e_msg : Optional[Union[None, str]]
        Custom error message. Will overwrite default message from error objects.
        [default=None]
    exit_code : int
        Exit code when error is raised. Cannot be lower than 1. [Default=1]

    Return
    ------
        None

    Raises
    ------
    TypeError
        Raised if error object is not an error type\n
        Raised if e_msg is not a string type\n
        Raised if exit_code is not an integer type\n

    ValueError
        Raised if exit_code is a negative number or is equal to 0

    """

    # type checking (should learn more about tpye guards)
    if not isinstance(error, (BaseException, Exception)):
        raise TypeError("'error' must be and Exception or BaseException types")

    # type checking e_msg
    if not isinstance(e_msg, str) and e_msg is not None:
        raise TypeError("'e_msg' must be a string type")

    # error code checking
    if not isinstance(exit_code, int):
        raise TypeError("'exit_code', must be integer type")
    elif exit_code == 0:
        raise ValueError("exit error codes cannot be 0")
    elif exit_code < 0:
        raise ValueError("exit codes cannot be negative numbers")

    # formatting error message
    error_type = f"{error.__class__.__name__}:"

    # if no custom message provided, default to exception default mesasge
    if e_msg is None:
        e_msg = str(error)

    print(f"{error_type} {e_msg}")
    sys.exit(exit_code)
