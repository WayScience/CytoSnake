#------------------------------------------------------------
# args_checker.py
#
# Argument checking module that checks and validates user input
# arguments from cytopipe's cli interface,
#------------------------------------------------------------
from dataclasses import dataclass

@dataclass
class CytoArgs:
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


def check_args(mode: str, args_list) -> bool:
    """Checks arguments 

    Parameters
    ----------
    mode : str
        _description_
    args_list : _type_
        _description_

    Returns
    -------
    bool
        _description_
    """
    pass