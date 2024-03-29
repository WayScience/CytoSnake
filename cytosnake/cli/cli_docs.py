run_doc = """
Command line::
Run Documentation

CytoSnake's run mode allows users to execute CytoSnake's available workflows

Usage:
    cytosnake run workflow [-c]
    cytosnake run help

Required Arguments:
    workflow            Name of the workflow to execute

Optional Arguments:
    -c, --max_cores     Maximum number of cores to use for the workflow
                        default=1
    --lock              Enable auto unlocking mechanism for snakemake workflows.
                        Directory becomes locked when workflow is executed. if
                        any interruptions has occurred, if True, the directory
                        will be automatically unlocked, else, it will remain
                        locked. Default is False.
    --force             Force re-run of the workflow. This means generated files
                        will be over-written with the outputs produced from the
                        forced re-run
Available Workflows:

cp_process              Workflow for analyzing cell morphology reads
                        obtained from CellProfiler
dp_process              Workflow for analyzing cell morphology reads
                        obtained from DeepProfiler


Help Arguments:
    help                Displays CytoSnake's Run mode documentation
"""

init_doc = """
Command line::
Init mode Documentation

CytoSnake's init mode allows user to setup up the require files for processing.

Usage:
    cytosnake init [-d] [-m] [-b] [-p]

Required Arguments:
    -d, --data          List of plate data files
    -m, --metadata      Path to metadata directory

Optional Arguments:
    -b, --barcode       Path to file containing barcode labeling. This is used
                        for cell morphology reads obtained by CellProfiler.
                        [Default=None]
    --datatype          Datatype flag helps CytoSnake in how to setup the input
                        files for processing.
                        [Choices = "cell_profiler", "deep_profiler"]
                        [Default="cell_profiler"]

Help Arguments:
    help                Displays CytoSnake's init mode documentation
"""

cli_docs = f"""
CytoSnake Documentation

CytoSnake's command line interface

Usage:
    cytosnake [mode] [ mode options]
    cyotsnake help

Required Arguments:
    mode                CytoSnake cli instruction on what to execute. There are
                        3 modes [init, run, help]. Init setups ups the required
                        files for processing. Run executes CytoSnake's workflows.
                        Help displays the help message documentation.
Help Argument:
    help                Displays CytoSnake's CLI help and mode documentation

Mode Documentations:
{init_doc.replace("Command line::", "")}
{run_doc.replace("Command line::", "")}
"""
