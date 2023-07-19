# Contributing to CytoSnake

First of all, thank you for contributing to CytoSnake!

This document contains guidelines on how to most effectively contribute to the CytoSnake codebase.

If you are stuck, please feel free to ask any questions or ask for help.

## Table of contents

[Code of conduct](#code-of-conduct)

[Quick links](#quick-links)

[How can I contribute?](#how-can-i-contribute)

- [Bug reporting](#bug-reporting)
- [Suggesting enhancements](#suggesting-enhancements)
- [Your first code contribution](#your-first-code-contribution)
- [Pull requests](#pull-requests)
- [Documentation](#documentation)

[Style guides](#style-guides)

- [Git commit messages](#git-commit-messages)
- [Python style guide](#python-style-guide)
- [Documentation style guide](#documentation-style-guide)
- [Pre-commit](#pre-commit)

[Advanced development](#advanced-development)

- [Testing](#testing)
- [Dry runs](#dry-runs)

## Code of conduct

This project and everyone participating in it is governed by our [code of conduct](CODE_OF_CONDUCT.md).
By participating, you are expected to uphold this code.
Please report unacceptable behavior to cytodata.info@gmail.com.

## Quick links

- Documentation: https://cytosnake.readthedocs.io
- Issue tracker: https://github.com/WayScience/CytoSnake/issues
- Code coverage: https://app.codecov.io/gh/WayScience/CytoSnake

## How can I contribute?

### Bug reporting

We love hearing about use-cases when our software does not work.
This provides us an opportunity to improve.
However, in order for us to fix a bug, you need to tell us exactly what went wrong.

When you report a bug, please be prepared to tell us as much pertinent information as possible.
This information includes:

- The pycytominer version you’re using
- The format of input data
- Copy and paste two pieces of information: 1) your command and 2) the specific error message
- What you’ve tried to overcome the bug

Please provide this information as an issue in the repository: https://github.com/WayScience/CytoSnake/issues

Please also search the issues (and documentation) for an existing solution.
It’s possible we solved your bug already!
If you find an issue already describing your bug, please add a comment to the issue instead of opening a new one.

### Suggesting enhancements

We’re deeply committed to a simple, intuitive user experience, and to support core profiling pipeline data processing.
This commitment requires a good relationship, and open communication, with our users.

We encourage you to propose enhancements to improve the pycytominer package.

First, figure out if your proposal is already implemented, by reading the documentation!
Next, check the issues (https://github.com/WayScience/CytoSnake/issues) to see if someone else has already proposed the enhancement you have in mind.
If you do find the suggestion, please comment on the existing issue noting that you are also interested in this functionality.
If you do not find the suggestion, please open a new issue and clearly document the specific enhancement and why it would be helpful for your particular use case.

Please provide your enhancement suggestions as an issue in the repository:

### Your first code contribution

Contributing code for the first time can be a daunting task.
However, in our community, we strive to be as welcoming as possible to newcomers, while ensuring rigorous software development practices.

The first thing to figure out is exactly what you’re going to contribute! We describe all future work as individual [github issues](https://github.com/WayScience/CytoSnake/issues).

If you want to contribute code that we haven’t already outlined, please start a discussion in a new issue before actually writing any code.
A discussion will clarify the new code and reduce merge time.
Plus, it’s possible that your contribution belongs in a different code base, and we do not want to waste your time (or ours)!

### Pull requests

After you’ve decided to contribute code and have written it up, now it is time to file a pull request.
We specifically follow a [forked pull request model](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request-from-a-fork).
Please create a fork of the CytoSnake repository, clone the fork, and then create a new, feature-specific branch.
Once you make the necessary changes on this branch, you should file a pull request to incorporate your changes into the main CytoSnake repository.

The content and description of your pull request are directly related to the speed at which we are able to review, approve, and merge your contribution into pycytominer.
To ensure an efficient review process please perform the following steps:

1. Follow all instructions in the [pull request template](.github/PULL_REQUEST_TEMPLATE.md)
2. Triple check that your pull request is only adding _one_ specific feature. Small, bite-sized pull requests move so much faster than large pull requests.
3. After submitting your pull request, ensure that your contribution passes all status checks (e.g. passes all tests)

All pull requests must be reviewed and approved by at least one project maintainer in order to be merged.
We will do our best to review the code addition in a timely fashion.
Ensuring that you follow all steps above will increase our speed and ability to review.
We will check for accuracy, style, code coverage, and scope.

### Documentation

We use [sphinx](https://www.sphinx-doc.org/en/master/index.html) for autodocumentation of docstrings, using the [napoleon extenstion](https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html) to parse [NumPy style docstrings](https://numpydoc.readthedocs.io/en/latest/format.html), implemented with a [furo](https://pradyunsg.me/furo/) theme.
We host our documentation on [readthedocs.org](https://readthedocs.org/) at [https://cytosnake.readthedocs.io/](https://cytosnake.readthedocs.io/).

To build and test changes to the docs locally, run the following command:

```bash
sphinx-build -b html docs build
```

See [`docs/conf.py`](docs/conf.py) for full documentation configuration.

## Style guides

Please follow all style guides to the best of your abilities.

### Git commit messages

For all commit messages, please use a short phrase that describes the specific change.
For example, “Add feature to check normalization method string” is much preferred to “change code”.
When appropriate, reference issues (via `#` plus number) .

### Python style guide

For python code style, we use [black](https://github.com/psf/black).
Please use black before committing any code.
We will not accept code contributions that do not use black.
If you have set up your development environment using one of the dev container options specified above, the containers will install all required formatting tools, which will run automatically on any modified files before commits (using a tool called [pre-commit](https://pre-commit.com/)).

### Documentation style guide

We use the [numpy documentation style guide](https://numpydoc.readthedocs.io/en/latest/format.html).
We also use [prettier](https://prettier.io/) for automatic formatting of markdown, json and yaml files.

### Pre-commit

[pre-commit](https://github.com/pre-commit/pre-commit) is a package that allows developers to run and check their code before commit changes.

The `pre-commit` configurations are found within the `.pre-commit-config.yaml` file, and installs the technologies in the `.git/hooks` directory.
This directory contains a list of formatting technologies used in order to ensure that our code meets formatting syntax and style standards.

#### How to use it

`pre-commits` and its dependencies already come included within the environment file.

Please make sure to activate the conda environment when applying changes in order to activate our pre-commits and its dependencies.
When changes are applied to the source code, the `pre-commit` workflow will be automatically executed on changed files.
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

## Advanced development

### Testing

![CytoSnake-Testing-Framework-v1](https://github.com/axiomcura/CytoSnake/assets/31600622/010f1d8b-f5ce-4f99-b035-7a036dd12ae0)

Our testing framework is split into 3 testing modules to ensure that our Workflow can maintain reproducibility, speed, and usability.

Our **functional** tests focus on the user perspective to verify if the functions work as expected.
The **functional** testing modules contain different sets of user-based parameters, which determine if `CytoSnake` can conduct what is set by the user.

The next module is the **unit**  tests, which individually test every core implementation of `CytoSnake`.
The targeted audience for the **unit** tests is developers to verify that every core component (modules, functions, objects) works as expected.
Each of these tests are written in isolation, meaning that the core implementations are being tested as if it was not part of the system.
This makes it extremely easy to test edge cases for our core implementations and allows for identification of bugs quickly.

Lastly, the **workflow** tests will test all available workflows that `CytoSnake` contains.
The tests attempt to verify that each analytical step within the workflow produces the expected input and output paths.
Pathing verification is also tested when testing for modularization, where specific workflows steps can be imported in other workflows.
Another important aspect of our workflow testing module is testing reproducibility.
These tests contain expected outputs along with associated parameters.
The workflow tests take in config files as input that is supposed to reproduce the known expected outputs.

### Dry runs

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
