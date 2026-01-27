# Hematopoietic stem cells population size estimate.

**\[stable\]** Hematopoietic stem cells population size estimate with
capture-recapture models.

## Usage

``` r
HSC_population_size_estimate(
  x,
  metadata,
  stable_timepoints = NULL,
  aggregation_key = c("SubjectID", "CellMarker", "Tissue", "TimePoint"),
  blood_lineages = blood_lineages_default(),
  timepoint_column = "TimePoint",
  seqCount_column = "seqCount_sum",
  fragmentEstimate_column = "fragmentEstimate_sum",
  seqCount_threshold = 3,
  fragmentEstimate_threshold = 3,
  nIS_threshold = 5,
  cell_type = "MYELOID",
  tissue_type = "PB",
  max_workers = 4
)
```

## Arguments

- x:

  An aggregated integration matrix. See details.

- metadata:

  An aggregated association file. See details.

- stable_timepoints:

  A numeric vector or NULL if there are no stable time points. NOTE: the
  vector is NOT intended as a sequence min-max, every stable time point
  has to be specified individually

- aggregation_key:

  A character vector indicating the key used for aggregating x and
  metadata. Note that x and metadata should always be aggregated with
  the same key.

- blood_lineages:

  A data frame containing information on the blood lineages. Users can
  supply their own, provided the columns `CellMarker` and `CellType` are
  present.

- timepoint_column:

  What is the name of the time point column to use? Note that this
  column must be present in the key.

- seqCount_column:

  What is the name of the column in x containing the values of sequence
  count quantification?

- fragmentEstimate_column:

  What is the name of the column in x containing the values of fragment
  estimate quantification? If fragment estimate is not present in the
  matrix, param should be set to `NULL`.

- seqCount_threshold:

  A single numeric value. After re-aggregating `x`, rows with a value
  greater or equal will be kept, the others will be discarded.

- fragmentEstimate_threshold:

  A single numeric value. Threshold value for fragment estimate, see
  details.

- nIS_threshold:

  A single numeric value. If a group (row) in the metadata data frame
  has a count of distinct integration sites strictly greater than this
  number it will be kept, otherwise discarded.

- cell_type:

  The cell types to include in the models. Note that the matching is
  case-insensitive.

- tissue_type:

  The tissue types to include in the models. Note that the matching is
  case-insensitive.

- max_workers:

  Maximum parallel workers allowed

## Value

A data frame with the results of the estimates

## Input formats

Both `x` and `metadata` should be supplied to the function in aggregated
format (ideally through the use of
[`aggregate_metadata`](https://calabrialab.github.io/ISAnalytics/dev/reference/aggregate_metadata.md)
and
[`aggregate_values_by_key`](https://calabrialab.github.io/ISAnalytics/dev/reference/aggregate_values_by_key.md)).
Note that the `aggregation_key`, aka the vector of column names used for
aggregation, must contain at least the columns associated with the tags
`subject`, `cell_marker`, `tissue` and a time point column (the user can
specify the name of the column in the argument `timepoint_column`).

## Specifying more than one group

Groups for the estimates are computed as a pair of cell type and tissue.
If the user wishes to compute estimates for more than one combination of
cell type and tissue, it is possible to specify them as character
vectors to the fields `cell_type` and `tissue_type` respectively, noting
that:

- Vectors must have the same length or one of the 2 has to be of length
  1

- It is a responsibility of the user to check whether the combination
  exists in the dataset provided.

Example:

    estimate <- HSC_population_size_estimate(
        x = aggreg,
        metadata = aggreg_meta,
        cell_type = c("MYELOID", "T", "B"),
        tissue_type = "PB"
    )

    # Evaluated groups will be:
    # - MYELOID PB
    # - T PB
    # - B PB

Note that estimates are computed individually for each group.

## On time points

If `stable_timepoints` is a vector with length \> 1, the function will
look for the first available stable time point and slice the data from
that time point onward. If `NULL` is supplied instead, it means there
are no stable time points available. Note that 0 time points are ALWAYS
discarded. Also, to be included in the analysis, a group must have at
least 2 distinct non-zero time points. NOTE: the vector passed has to
contain all individual time points, not just the minimum and maximum

## Setting a threshold for fragment estimate

If fragment estimate is present in the input matrix, the filtering logic
changes slightly: rows in the original matrix are kept if the sequence
count value is greater or equal than the `seqCount_threshold` AND the
fragment estimate value is greater or equal to the
`fragmentEstimate_threshold` IF PRESENT (non-zero value). This means
that for rows that miss fragment estimate, the filtering logic will be
applied only on sequence count. If the user wishes not to use the
combined filtering with fragment estimate, simply set
`fragmentEstimate_threshold = 0`.

## Required tags

The function will explicitly check for the presence of these tags:

- subject

- tissue

- cell_marker

## See also

Other Analysis functions:
[`CIS_grubbs()`](https://calabrialab.github.io/ISAnalytics/dev/reference/CIS_grubbs.md),
[`compute_abundance()`](https://calabrialab.github.io/ISAnalytics/dev/reference/compute_abundance.md),
[`cumulative_is()`](https://calabrialab.github.io/ISAnalytics/dev/reference/cumulative_is.md),
[`gene_frequency_fisher()`](https://calabrialab.github.io/ISAnalytics/dev/reference/gene_frequency_fisher.md),
[`is_sharing()`](https://calabrialab.github.io/ISAnalytics/dev/reference/is_sharing.md),
[`iss_source()`](https://calabrialab.github.io/ISAnalytics/dev/reference/iss_source.md),
[`sample_statistics()`](https://calabrialab.github.io/ISAnalytics/dev/reference/sample_statistics.md),
[`top_integrations()`](https://calabrialab.github.io/ISAnalytics/dev/reference/top_integrations.md),
[`top_targeted_genes()`](https://calabrialab.github.io/ISAnalytics/dev/reference/top_targeted_genes.md)

## Examples

``` r
data("integration_matrices", package = "ISAnalytics")
data("association_file", package = "ISAnalytics")
aggreg <- aggregate_values_by_key(
    x = integration_matrices,
    association_file = association_file,
    value_cols = c("seqCount", "fragmentEstimate")
)
aggreg_meta <- aggregate_metadata(association_file = association_file)
estimate <- HSC_population_size_estimate(
    x = aggreg,
    metadata = aggreg_meta,
    fragmentEstimate_column = NULL,
    stable_timepoints = c(90, 180, 360),
    cell_type = "Other"
)
#> Calculating number of IS for each group...
```
