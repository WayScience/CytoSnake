"""
module: test_utils.py

test_utils.py contain additional functions that enhance testing capabilities, providing
extra functionality for conducting comprehensive and robust tests.
"""

def get_raised_error(traceback: str) -> str:
    """Parses traceback and attempts to obtain raised exception error.

    Traceback is parsed in this order:
    1. split by new lines
    2. grab the last line as it contains the raised exception and message
    3. split by ":" to separate exception name and exception message
    4. grab the first element since it contains that path to exception
    5. split by "." and grab last element, which is the exception name

    Parameters
    ----------
    traceback : str
        complete traceback generated by executing CLI

    Returns
    -------
    str
        name of raised exception error
    """

    # returns exception name, refer to function documentation to understand
    # the order of parsing the traceback to obtain exception name.
    return traceback.splitlines()[-1].split(":")[0].split(".")[-1]