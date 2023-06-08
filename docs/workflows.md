# CytoSnake Workflows

```{toctree}
:maxdepth: 4
```

Below are the available workflow modules in `CytoSnake`

## cp_process

Generates consensus profiles.

Consensus profiles contain unique signatures that are associated with a specific state of the cell (e.g., control, DMSO, chemical/genetic perturbations, etc.).

The development of this workflow was heavily influenced by this study:
<https://github.com/broadinstitute/cell-health>

**included modules**:

- common
- aggregate
- annotate
- normalize
- feature_select
- generate_consensus

**parameters**

input:

- **profile**: single-cell morphology datasets
- **metadata**: metadata directory associated with plate data.
- **barcode** (optional): file containing plate data to plate map pairings. Default is `None`.

**outputs**:

- **profiles**: aggregated, annotated, normalized, selected features, and consensus profiles.
- **cell counts**: Cell counts per well

## cp_process_singlecells

Converts SQLite plate data into parquet and returns selected features from
single-cell morphology profiles.

The development of this workflow was heavily influenced by:
<https://github.com/WayScience/Benchmarking_NF1_data/tree/main/3_extracting_features>

**included modules**

- common
- cytotable_convert
- annotate
- normalize
- feature_select

**parameters**

inputs:

- **profiles**: single-cell morphology datasets in SQLite format.
- **metadata**: metadata file associated with single-cell morphology dataset.

outputs:

- **profiles**: annotated, normalized, and feature selected profiles.
