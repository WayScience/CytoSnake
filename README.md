# CytoPipe

"Reproducible pipelines for processing high-dimensional systems morphology data with snakemake 🐍 "

# Table of contents

- [CytoPipe](#cytopipe)
- [Table of contents](#table-of-contents)
  - [Formatting](#formatting)
    - [Black](#black)
    - [Snakefmt](#snakefmt)
- [Testing](#testing)
  - [Dry Runs](#dry-runs)

## Formatting

### Black
`Black` is a python code formatter that follows `PEP8` rules in order to maintain source code in a readable fromat.

To install `Black`:
```
pip install black
```
More in depth informaton can be found in the [Documentation](https://black.readthedocs.io/en/stable/)


### Snakefmt

`snakefmt` is a `snakemake` code formatter that is developed and maintained by the developers of `snakemake`.

In order to use `snakefmt`, `snakemake` must be installed, which can be found here: [snakemake-installation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

`snakefmt` installation processes can be found in the [snakefmt](https://github.com/snakemake/snakefmt#install) repository.

# Testing

## Dry Runs

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
