<img height="200" src="https://raw.githubusercontent.com/WayScience/CytoSnake/main/logo/with-text-for-light-bg.png?raw=true">

_CytoSnake: Orchestrating reproducible pipelines for processing high-dimensional systems morphology data with snakemake_

Erik Serrano, Gregory P. Way
University of Colorado Anschutz School of Medicine

## Table of contents

- [Table of contents](#table-of-contents)
- [About](#about)
- [Installation](#installation)
- [Workflows](#workflows)

## About

CytoSnake is a command line interface (CLI) tool that orchestrates reproducible workflows that process high-dimensional single-cell morphology features extracted from microscopy images.
CytoSnake's workflows are written in [`Snakemake`](https://github.com/snakemake/snakemake), which is a well established workflow manager that facilitates data reproducibility, scalability, and modularity.

`CytoSnake` makes it easy for user to process high-dimensional cell morphology data as it requires straightforward inputs and parameters.
Below is an example on how to execute `CytoSnake` once installed:

```bash
# setting up directory
cytosnake init -d <FILES or LIST OF FILES> -m <METADATA DIR> -b <BARCODE>

# executing workflow
cytosnake run <WORKFLOW>
```

**note**: `-b` is optional, it is used if there are multiple platemap files

## Installation

Install `CytoSnake` with `conda`.

```bash
conda install -c bioconda cytosnake
```

>**NOTE**: In case you don't currently have `conda` installed on your system, you can access the documentation [here](https://docs.conda.io/en/latest/miniconda.html).
We recommend using Miniconda, primarily due to its lightweight installation process.

To check if `CytoSnake` has been successfully installed, simply type `cytosnake help` to see the CLI documentation:

```bash
cytosnake help
```

## Workflows

CytoSnake workflows are the main instructions on how your data is going to be processed.
Each workflow comes with its appropriate configuration file.

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

Above is a portion of the listed configs from the `cp_process` workflow.
Each block represents an analytical specific step that is conducted within the workflow.
In this example, `annotate_configs` and `aggregate_configs` are separate steps that occur within the `cp_process` workflow.
Each block has the `params` parameter, which are the parameters associated with the analytical step.
Users can edit these parameters from the defaults if they want their workflow to analyze their data in a specific way.

Overall, each workflow will have a designated workflow config file.
It will contain all the steps conducted in the workflow, and users have the option to change the default parameters that are specific to their dataset.
