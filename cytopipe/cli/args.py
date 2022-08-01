#------------------------------------------------------------
# args.py
# 
# Responsible for centralizing and handling all user inputted 
# parameters for cytopipe.
#
# Contains different modes with their respected parameters
#------------------------------------------------------------
import argparse

def parse_init_args(args_list: list) -> argparse.Namespace:
    """Parses user inputs for 

    Parameters
    ----------
    params : list
        User defined input parameters

    Returns
    -------
    argparse.Namespace
        parsed cli parameters
    """
    
    # CLI inputs
    parser = argparse.ArgumentParser()
    parser.add_argument("init")
    parser.add_argument("-d", "--data", nargs="+", required=True, help="list of plate data files")
    parser.add_argument("-m", "--metadata", type=str, required=True, help="path to metadata directory")
    parser.add_argument("-b", "--barcode", type=str, required=True, help="path to barcodes file")
    parser.add_argument("-p", "--platemap", type=str, required=True, help="path to platemaps file")
    args = parser.parse_args(args_list)

    return args