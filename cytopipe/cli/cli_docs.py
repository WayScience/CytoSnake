from textwrap import dedent

run_doc = """
Command line::
Run Documentation

Cytopipe's run mode allows users to execute Cytopipes available workflows

Usage:
    cytopipe run workflow [-c]  
    cytopipe run help            

Required Arguments:
    workflow            Name of the workflow to execute

Optional Arguments:
    -c, --max_cores     Maximum number of cores to use for the workflow 
                        default=1

Help Arguments:
    help                Displays cytopipe's Run mode documentation j
"""

init_doc = """
Command line::
Init mode Documentation

Cytopipe's init mode allows user to setup up the require files for processing.

Usage:
    cytopipe init [-d] [-m] [-b] [-p]

Required Arguments:
    -d, --data          List of plate data files           
    -m, --metadata      Path to metadata directory
    -b, --barcode      Path to file containing barcode labeling
    -p, --platemap      Path to platemap file

Help Arguments:
    help                Displays cytopipe's init mode documentation                
"""

cli_docs = f"""
Cytopipe Documentation

Cytpipe's command line interface 

Usage:
    cytopipe [mode] [ mode options]
    cytopipe help

Required Arguments:
    mode                Cytopipe cli instruction on what to execute. There are
                        3 modes [init, run, help]. Init setups ups the required
                        files for processing. Run executes cytopipe's workflows.
                        Help displays the help message documentation. 
Help Argument:
    help                Displays cytopipe's CLI help and mode documentation

Mode Documentations:
{init_doc.replace("Command line::", "")}
{run_doc.replace("Command line::", "")}
"""
