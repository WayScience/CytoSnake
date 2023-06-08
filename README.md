# 🐍CytoSnake

"Reproducible pipelines for processing high-dimensional systems morphology data with snakemake 🐍 "

## Table of contents

- [🐍CytoSnake](#cytosnake)
  - [Table of contents](#table-of-contents)
  - [About](#about)
  - [Installation](#installation)
  - [Workflows](#workflows)
  - [💻 Development](#-development)
    - [Pre-commits](#pre-commits)
      - [how to use it](#how-to-use-it)
    - [Testing](#testing)
    - [Dry Runs](#dry-runs)

## About

CytoSnake is a command line interface (CLI) tool that contains reproducible workflows that analysis high dimensional single-cell morphology datasets.
CytoSnake's workflows are written in [`Snakemake`](https://github.com/snakemake/snakemake), which is a well established workflow manager that contains important features like data reproducibility, scalability and modularity.

`CytoSnake` make it easy for user to interact as it requires little inputs and parameters. below is an example on how to execute `CytoSnake` once installed

```text

```

## Installation

First, install `CytoSnake` into your local machine:

```text
git clone https://github.com/WayScience/CytoSnake.git

```

After cloning the repository, go into the `CytoSnake/` directory and create the `CytoSnake` environment

```text
conda env create -f cytosnake_env.yaml && conda activate cytosnake
```

Next is to install the `CytoSnake` module into your newly created environment

```text
pip install -e .
```

After this step, `CytoSnake` is installed. To check if `CytoSnake` is properly installed, simply type:

```text
cytosnake
```

you should see the CLI documentation pop out.

## Workflows

CytoSnake workflows are the main instructions on how your data is going be processed.
Each workflow comes with it's appropriate configuration file.

Here is an example below:

```yaml
annotate_configs:
  params:
    input_data: plate_data
    join_on:
      - Metadata_well_position
      - Image_Metadata_Well
    add_metadata_id_to_platemap: True
    format_broad_cmap: False
    clean_cellprofiler: True
    external_metadata: "none"
    external_join_left: "none"
    external_join_right: "none"
    compression_options:
      method: "gzip"
      mtime: 1
    float_format: null
    cmap_args: {}

aggregate_configs:
  params:
    input_data: annotated
    strata:
      - Metadata_Plate
      - Metadata_Well
    features: infer
    operation: median
    output_file: none
    compute_object_count: False
    object_feature: Metadata_ObjectNumber
    subset_data_df: none
    compression_options:
      method: gzip
      mtime: 1
    float_format: null

```

Here is a portion of the listed configs from the `cp_process` workflow.
Each block represents an analytical specific step that is conducted within the workflow.
In this example, `annotate_configs` and `aggregate_configs` are separate steps that occurs within the `cp_process` workflow.
Each block has the `params` parameter, which are the parameters associated with the analytical step.
Users can edit these parameters from the defaults if they want their workflow to analyze their data in a specific way.

Overall, each workflow will have a designated workflow config file.
It will contain all the steps conducted in the workflow and users have the options to change the preset that is specific to their dataset.

## 💻 Development

Below is the list of technologies used when developing CytoSnake.

### Pre-commits

[pre-commit](https://github.com/pre-commit/pre-commit) is a package that allows developers to run and check their code before commit changes.

The `pre-commit` configurations are found within `.pre-commit-config.yaml` file, and installs the technologies in the `.git/hooks` directory.
Within this directory, it contains a list of formatting technologies used in order to ensure that our code meets formatting standards like syntax and style.

#### How to use it

`pre-commits` and its dependencies already come included within the environment file.

Please make sure to active the conda environment when applying changes in order to activate our pre-commits and its dependencies.
When changes are applied to the source code, the `pre-commit` workflow will automatically executed on changed files.
However, if there is a need to run `pre-commit` on both edited and unedited files, you can directly execute the `pre-commit` directly with the `--all-files` parameter:

```sh
pre-commit run --all-files
```

This will let `pre-commit` know to execute the workflow on all files.

Below is the list of technologies used within the `pre-commit` workflow:

- [pycln](https://github.com/hadialqattan/pycln) : removes unused imports
- [isort](https://github.com/PyCQA/isort) : import sorter that sorts imports alphabetically, sections and by type
- [black](https://github.com/psf/black) : Python code formatter that supports [`PEP8`](https://peps.python.org/pep-0008/) styling.
- [black[Jupyter]](<https://github.com/drillan/jupyter-black>) : black formatter for notebooks
- [sourcery](https://github.com/sourcery-ai/sourcery) : AI-based code refactor. Provides suggestions to improve readability
- [ruff](https://github.com/charliermarsh/ruff) : Python Linter
- [snakefmt](https://github.com/snakemake/snakefmt) : Snakemake file formatter
- [pre-commit-hooks](https://github.com/pre-commit/pre-commit) : library of additional pre-commit technologies
  - remove trailing white spaces
  - remove space and tab mixtures
  - json formatter

### Testing

![CytoSnake-Testing-Framework-v1](https://github.com/axiomcura/CytoSnake/assets/31600622/010f1d8b-f5ce-4f99-b035-7a036dd12ae0)

Our testing frame are split into 3 testing modules to ensure that our Workflow is able to maintain reproducibility, speed and usability.

Our **functional** tests focuses on the user perspective to verify if the expected functionalities. `cytosnake/tests/functional`.
The functional testing modules contains different sets of user based parameter, which  weather our `CytoSnake` is able to conduct what is asked from the user. `

The next module is the **unit test** module, which individually tests every the core implementation of `CytoSnake`
The targeted audience is mainly for developers to verify that every single core component (modules, functions, objects) work as as expected.
Each of these tests are written in isolation meaning that the core implementation are being tested as if it was not part of the system.
This makes it extremely easy to tests edges cases for our core implementations and allows to quickly identify bugs.

Lastly, the **workflow testing** module tests all available workflows that `CytoSnake` contains.
The tests attempt to verify that each analytical steps within the workflow produces the expected input and output paths.
Pathing verification is also tested when testing for modularization, where specific workflows steps can be imported in other workflows.
Another important aspect of our workflow testing module is testing reproducibility.
These test contains expected outputs along with associated parameters.
The workflow tests takes in config files as inputs that is supposed to reproduce the known expected outputs.

### Dry Runs

`Snakemake`'s `dry-runs` is a simple test that ensures that all the required inputs and expected outputs are satisfied.

No execution is conducted, therefore the dry run **will not capture** any potential errors found within the executables.

Below is an example on how to conduct a `dry-run` on `snakemake`:

```text
snakemake -n
```

Using the `-n`, by default, will conduct a `dry-run` to all rules that exists within the workflow.

One can specify a specific rule to conduct a `dry-run` test by inputting the name of the rule after the `-n` flag.

This will only check all the required inputs and expected outputs that are specific to the specified rule.

The example below executes a `dry-run` test to only the `annotate` rule.

```text
snakemake -n annotate
```

output:

```text
Building DAG of jobs...
Job stats:
job          count    min threads    max threads
---------  -------  -------------  -------------
aggregate        1              1              1
annotate         1              1              1
total            2              1              1


[Mon Apr 18 10:42:55 2022]
rule aggregate:
    input: data/SQ00014614.sqlite, data/barcode_platemap.csv, data/metadata
    output: results/preprocessing/SQ00014614_cell_counts.tsv, results/preprocessing/SQ00014614_aggregate.csv.gz
    jobid: 1
    resources: tmpdir=/tmp

[Mon Apr 18 10:42:55 2022]
rule annotate:
    input: data/barcode_platemap.csv, results/preprocessing/SQ00014614_aggregate.csv.gz
    output: results/preprocessing/SQ00014614_augmented.csv.gz
    jobid: 0
    resources: tmpdir=/tmp

Job stats:
job          count    min threads    max threads
---------  -------  -------------  -------------
aggregate        1              1              1
annotate         1              1              1
total            2              1              1
```

The resulting output is a summary of the `dry-run` and it shows that two rules were tested despite only specifying one.
This is due to the `aggregate` rule creating an output that is a required input for the `annotate` rule; therefore, indicating that the `aggregate` rule is a requirement for the `annotate` rule.

The `job stats` gives a summary which rules were checked in the `dry-run` and the resources used.

More information of `snakemake`'s `dry-runs` can be found [here](https://snakemake.readthedocs.io/en/v5.1.4/executable.html#useful-command-line-arguments)
