# 🐍CytoSnake

"Reproducible pipelines for processing high-dimensional systems morphology data with snakemake 🐍 "

## Table of contents

- [🐍CytoSnake](#cytosnake)
  - [Table of contents](#table-of-contents)
  - [About](#about)
  - [Usage](#usage)
  - [💻 Development](#-development)
    - [Pre-commits](#pre-commits)
      - [how to use it](#how-to-use-it)
    - [Testing](#testing)
      - [Functional](#functional)
      - [Unit Testing](#unit-testing)
      - [workflow](#workflow)
    - [Dry Runs](#dry-runs)

## About
[About description here]


## Usage

[Usage description here]

## 💻 Development

Below is the list of technologies used when developing CytoSnake.

### Pre-commits

[pre-commit](https://github.com/pre-commit/pre-commit) is a package that allows developers to run and check their code before commit changes.

The `pre-commit` configurations are found within `.pre-commit-config.yaml` file and installs the technolgies in the `.git/hooks` directory where it contains a list of formatting technologies used in order to ensure that our code meets formatting standrds like: formatting, syntax and style.

#### how to use it

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

#### Functional

Our functional tests focuses on more on the user side `cytosnake/tests/functional`.

Users can provided different paramter inputs. Our tests

#### Unit Testing

[Unit testing docs ]

#### workflow

[Workflow testing Docs]

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
