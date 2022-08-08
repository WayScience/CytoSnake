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
    --conda_env         Enables workflow to use conda environment, default=True
    --lock              Enable auto unlocking mechanism for snakemake workflows.
                        Directory becomes locked when workflow is executed. if 
                        any interruptions has occurred, if True, the directory
                        will be automatically unlocked, else, it will remain
                        locked. Default is False.
    --force             Force re-run of the workflow. This means generated files
                        will be over-written with the outputs produced from the
                        forced re-run

Help Arguments:
    help                Displays cytopipe's Run mode documentation
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
    -b, --barcode       Path to file containing barcode labeling

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
