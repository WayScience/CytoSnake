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
    """Based Exception if incorrect values are passed"""


class BaseExecutorException(RuntimeError):
    """Based exception related to cytopipe execution"""


# ------------------------------
# Cytopipe specific errors
# ------------------------------
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
