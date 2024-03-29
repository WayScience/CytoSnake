# CytoSnake Tutorial

```{toctree}
:maxdepth: 4
```

## About

This tutorial assumes that you have followed the installation steps and
you are ready to start taking off with CytoSnake!

Cytosnake is a command line interface (CLI) tool that contains a
multitude of workflows for analyzing morphological features obtained
from microscopy images of cells.

## Concepts

```{image} images/CytoSnake-cli_0.drawio.png
```

### Modes

Modes provide options on how the user can change the functionality of
CytoSnake. For example, if you would like to initialize your files for a
specific workflow, you can use the `init` mode:

```console
cytosnake init -d [<DATAFILES>] -m <METADATA> --data_type <DATATYPE>
```

- `DATAFILE` will refer to the raw data file(s) that you are going to
  analyze.
- `METADATA` refers to the associated metadata data directory that was
  generated along with the dataset
- `DATATYPE` flag tells cytosnake the origin of these morphology
  feature datasets (currently either CellProfiler or DeepProfiler).

The init mode initializes the provided input files into the appropriate
file structure that accommodates all the workflows available in
CytoSnake.

CytoSnake currently has three different types of modes, which are:

> 1. `init`: setup input files for workflows
> 2. `run`: execute a specific workflow
> 3. `help`: executes CytoSnake's CLI help documentation.

Example of using CytoSnake and its modes is added in the [Usage]
section.

### Configurations

CytoSnake has a configuration directory that allows users to change the
configurations for their specified workflows.

The configuration files are written in `.yaml` files, which contains
all the functions and its parameters used within the workflow. The
workflow's documentation provides information about the configuration
files involved within the workflow.

Each workflow possesses its own configurational file.

Below are the currently available workflows and the config files it accesses to conduct its processes.
accesses in order to conduct its processes.

#### **cp_process workflow docs**

| workflow name          | Path to config                                           | Documentation         |
| -------------- | ---------------------------------------------------------------- | --------------------- |
| cp_process     | ./CytoSnake/configs/analysis_configs/wf_configfs/cp_process.yaml |  [Workflow Docs](workflows.md#cp_process)|

**modules used:**
| Steps          | Path to config                                                   | Module Documentation  |
| -------------- | ---------------------------------------------------------------- | --------------------- |
| aggregate      | ./CytoSnake/configs/analysis_configs/aggregate_configs.yaml      | [aggregate docs](workflow-modules.md#aggregate)      |
| annotate       | ./CytoSnake/configs/analysis_configs/aggregate_configs.yaml      | [annotate docs](workflow-modules.md#annotate)      |
| normalize      | ./CytoSnake/configs/analysis_configs/normalize_configs.yaml      | [normalize docs](workflow-modules.md#normalize)      |
| feature_select | ./CytoSnake/configs/analysis_configs/feature_select_configs.yaml | [feature select docs](workflow-modules.md#feature-select) |
| consensus      | ./CytoSnake/configs/analysis_configs/consensus_configs.yaml      | [consensus docs](workflow-modules.md#generate-consensus)      |

#### **cp_process_singlecells workflow docs**

| workflow name          | Path to config                                           | Documentation         |
| -------------- | ---------------------------------------------------------------- | --------------------- |
| cp_process_singlecells | ./CytoSnake/configs/analysis_configs/wf_configfs/cp_process_singlecells.yaml |  [Workflow Docs](workflows.md#cp_process_singlecells)|

**modules used**:

| Steps          | Path to config                                                   | Module Documentation  |
| -------------- | ---------------------------------------------------------------- | --------------------- |
| convert      | ./CytoSnake/configs/analysis_configs/cytotable_convert.yaml        | [convert docs](workflow-modules.md#cytotable-convert)      |
| annotate       | ./CytoSnake/configs/analysis_configs/aggregate_configs.yaml      | [annotate docs](workflow-modules.md#annotate)       |
| normalize      | ./CytoSnake/configs/analysis_configs/normalize_configs.yaml      | [normalize docs](workflow-modules.md#normalize)      |
| feature_select | ./CytoSnake/configs/analysis_configs/feature_select_configs.yaml | [feature select docs](workflow-modules.md#feature-select) |

Users can easily find and change parameter values by accessing those
configurational files.

- **Steps** : instructions that the workflows
- **Path to config** : Location of the configurational files
- **Documentation** : Workflow documentation
- **Module Documentation** : Module documentation

## Documentation

To see CytoSnake's documentation, simply type:

```text
cytosnake help
```

This will display a large output into your terminal explaining all modes
and its parameters. If you are only interested in one, you can use the
`help` under any mode:

```text
# display help for run mode
cytosnake run help

# display help for init mode
cytosnake init run
```

### `init` mode documentation

Here are the list of parameters that CytoSnake's `init` mode currently
supports

| Parameters         | Documentation                                                                                                                                           |
| ------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Required Arguments |                                                                                                                                                         |
| `-d` `--data`      | List of plate data files                                                                                                                                |
| `-m` `--metadata`  | Path to metadata directory                                                                                                                              |
| Optional Arguments |                                                                                                                                                         |
| `-b` `--barcode`   | Path to file containing barcode labeling. This is used for cell morphology reads obtained by CellProfiler. \[Default=None\]                             |
| `--datatype`       | Datatype flag helps CytoSnake in how to setup the input files for processing. \[Choices = "cell_profiler" "deep_profiler"\] \[Default="cell_profiler"\] |

### `run` mode documentation

Here are the list of parameters that CytoSnake's `run` mode currently
supports

| Parameters         | Documentation                                                                                                                                                                                |
| ------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Required Arguments |                                                                                                                                                                                              |
| workflow           | Name of the workflow to execute                                                                                                                                                              |
| Optional Arguments |                                                                                                                                                                                              |
| `-c` `--max_core`  | Maximum number of cores to use for the workflow default=1                                                                                                                                    |
| `--lock`           | Directory becomes locked when workflow is executed. if any interruptions has occurred, if True, the directory will be automatically unlocked, else, it will remain locked. Default is False. |
| `--force`          | Force re-run of the workflow. This means generated files will be over-written with the outputs produced from the forced re-run                                                               |

## Usage

### Download data

In this usage tutorial, we will be using cell health datasets. ([Way, 2021](https://doi.org/10.1091/mbc.e20-12-0784))

You can download these datasets below (quite large files):

- [plate_data_1](https://nih.figshare.com/ndownloader/files/18506036): (10GB download)
- [plate_data_2](https://nih.figshare.com/ndownloader/files/18031619): (11GB download)
- [metadata_folder](https://minhaskamal.github.io/DownGit/#/home?url=https://github.com/broadinstitute/cell-health/tree/master/1.generate-profiles/data/metadata): Contains all associated perturbations per well
- [barcode](https://minhaskamal.github.io/DownGit/#/home?url=https://github.com/broadinstitute/cell-health/blob/master/1.generate-profiles/data/metadata/barcode_platemap.csv): Maps plate id with plate names

You can also use your dataset but some of the tasks that are being done
here are specific to the cell health dataset example.

### Setting up files

If you are using the downloaded datasets unzip the zip files in the
directory where the CytoSnake source is.

```text
unzip metadata.zip && unzip barcode_platemap.csv.zip
```

The first step is to prepare your files for analysis. This is simply
executed by typing:

```text
cytosnake init -d SQ00014613.sqlite SQ00014613.sqlite -d metadata -b barcode_platemap.csv
```

In instances where you may have a lot of data, CytoSnake supports
wildcard variables.

```text
cytosnake init -d *.sqlite -d metadata -b barcode_platemap..csv
```

Note that the data files under the `-d` parameter uses the wildcard `*` to
select all the `sqlite` files and places them into a single list of datafile
entries.

Wildcards are mainly suitable for selecting multiple files.

If there is an instance where you are going to use morphological
datasets obtained from DeepProfiler, then you must explicitly state the
datatype flag when using `init`:

Once entering the command, your out put should look like this:

```text
INFO: Formatting input files
INFO: Formatting complete!
```

If you are not receiving these messages, please refer to the [install.md] to see
if the installation process was done correctly.

So what just ocurred here?

`CytoSnake` now recognizes that your current directory is now known as the `ProjectDirectory`.

The `ProjectDirectory` is a way for `CytoSnake` to recognize that a project is occurring therefore preventing other projects from being created within the same directory.

The `ProjectDirectory` contains metadata information that also helps `CytoSnake` know what files have been initialized for downstream analysis.
You can find this information within the `.cytosnake` directory created after running the `init` mode.

### Running the workflow

In the `ProjectDirectory`,  a new folder called `./data` should appear in
your current directory. Inside the directory, it should contain symbolic
links of your data files that you have provided in the `init` mode. This
directory serves as a centralized location of data for the workflows to
have access too. Now that you have your data folder, you can simply
select which workflow to execute by using the `run` mode. Since the
cell-health dataset contains data extracted from CellProfiler, when we
will use the `cp_process` workflow.

```text
cytosnake run cp_process
```

These workflows contain their own environments, therefore there is no
need to download the dependencies that our workflows require. When the
the job is done, the last message you should see is:

```text
[Mon Sep 19 14:29:07 2022]
Finished job 0.
2 of 2 steps (100%) done
```

This indicates that all tasks within the workflow is complete.

### Accessing data

In your directory, `CytoSnake` produces a `results/` folder, which will
contain all the outputs generated from the workflow.
To list those
outputs, simply type:

```text
cd results/preprocessing/ && ls
```

This will take you to the directory where the generated outputs are and
lists all the files.

```text
consensus.tsv.gz                  SQ00014614_aggregate.csv.gz
SQ00014613_aggregate.csv.gz       SQ00014614_augmented.csv.gz
SQ00014613_augmented.csv.gz       SQ00014614_cell_counts.tsv
SQ00014613_cell_counts.tsv        SQ00014614_feature_select.csv.gz
SQ00014613_feature_select.csv.gz  SQ00014614_normalized.csv.gz
SQ00014613_normalized.csv.gz
```

These files contain different types of information that is denoted by
their suffix:

- `_cell_counts.tsv`: Number of cells in the dataset
- `_aggregate`: Refers to the aggregated dataset. Single cell dataset
  (your inputs) are aggregated into the “well” level.
- `_augmented`: A datasets contains metadata information in a per well
  level. For example, types of metadata can be: well position,
  treatments, controls, etc
- `_normalized` : normalized augmented dataset useful for feature
  selection.
- `_feature_select`: contains the selected morphological features that
  will be used to generate consensus profiles
- `_consensus`: is the consensus profile contains unique morphological
  signatures associated with a specific external treatment (drug,
  perturbations, controls (pos/neg), etc)
