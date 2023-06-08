# Workflow Modules

```{toctree}
:maxdepth: 4
```

## Aggregate

Aggregates single-cell profiles into aggregated profiles based on the given strata.

For example, users can configure `Metadata_Well` as their strata in order to
aggregate single-cell data into the Well level.

Utilize's pycytominer's aggregate module:
<https://github.com/cytomining/pycytominer/blob/c90438fd7c11ad8b1689c21db16dab1a5280de6c/pycytominer/aggregate.py>

**Parameters** (params)

inputs:

- **profile**: profile single-cell morphology dataset.
- **barcode**: file containing unique barcodes that maps to a specific plate
- **metadata**: metadata file associated with single-cell morphology dataset.

outputs:

- **profile** : aggregated datsaset
- **cell-counts** : csv file containg cell counts per well

## Annotate

Generates an annotated profile with given metadata and is stored instored
in the `results/` directory.

Utilizes pycytominer's annotate module:
<https://github.com/cytomining/pycytominer/blob/master/pycytominer/annotate.py>

**parameters**

inputs:

- **profiles**: single-cell or aggregate profiles.
- **barcode**: file containing unique barcodes that maps to a specific plate.
- **metadata**: metadata file associated with single-cell morphology dataset.

output:

- profile: annotated profile.

## Common

common.smk is a workflow module that sets up the expected input and output paths
for the main analytical workflow.

## Cytotable Convert

Converts single-cell morphology dataset to parquet format.

Utilizes CytoTable's convert workflow module:
<https://github.com/cytomining/CytoTable/blob/main/cytotable/convert.py>

**paramters**

inputs:

- **profiles**: single-cell profiles.
- **barcode**: file containing unique barcodes that maps to a specific plate.
- **metadata**: metadata file associated with single-cell morphology dataset.

outputs:

- **annotated**: converted single-cell morphology dataset.

## Feature Select

Performs feature selection based on this given profiles.

PyCytominer contains different operations to conduct its feature selection: variance_threshold, correlation_threshold, drop_na_columns, drop_outliers, and noise_removal.

Utilizes pycytominer's feature select module:
<https://github.com/cytomining/pycytominer/blob/master/pycytominer/feature_select.py>

**paramters**

inputs:

- **profile**: single-cell
- **barcode**: file containing unique barcodes that maps to a specific plate.
- **metadata**: metadata file associated with single-cell morphology dataset.

outputs:

- **profiles**: selected features profiles.
"""

## Generate consensus

Creates consensus profiles that reflects unique signatures associated with external factors.

Utilize's pycytominer's consensus module:
<https://github.com/cytomining/pycytominer/blob/master/pycytominer/consensus.py>

**parameters**

inputs:

- profile: selected features profiles
- barcode: file containing unique barcodes that maps to a specific plate.
- metadata: metadata file associated with single-cell morphology dataset.

output:

- **profile**: Consensus profile

## Normalize

Normalizing single-cell or aggregate features. Current default normalization
method is `standardize` other methods include:

Utlizes pycytominer's normalization module:
<https://github.com/cytomining/pycytominer/blob/c90438fd7c11ad8b1689c21db16dab1a5280de6c/pycytominer/normalize.py>

- **profiles**: single-cell morphology or annotated profiles.
- barcode: file containing unique barcodes that maps to a specific plate.
- metadata: metadata file associated with single-cell morphology dataset.

output

- **profiles**: normalized profiles.
