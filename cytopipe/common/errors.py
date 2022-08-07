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


# ------------------------------
# Base Exceptions
# ------------------------------
class BaseValueError(ValueError):
    """Base  exception if incorrect values are passed"""


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


class WorkflowFailedException(BaseWorkflowException):
    """Raised if a workflow fails during runtime"""
