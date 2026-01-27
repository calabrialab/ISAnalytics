# Scans input matrix to find and merge near integration sites.

**\[stable\]** This function scans the input integration matrix to
detect eventual integration sites that are too "near" to each other and
merges them into single integration sites adjusting their values if
needed.

## Usage

``` r
compute_near_integrations(
  x,
  threshold = 4,
  is_identity_tags = c("chromosome", "is_strand"),
  keep_criteria = c("max_value", "keep_first"),
  value_columns = c("seqCount", "fragmentEstimate"),
  max_value_column = "seqCount",
  sample_id_column = pcr_id_column(),
  additional_agg_lambda = list(.default = default_rec_agg_lambdas()),
  max_workers = 4,
  map_as_file = TRUE,
  file_path = default_report_path(),
  strand_specific = lifecycle::deprecated()
)
```

## Arguments

- x:

  An integration matrix

- threshold:

  A single integer that represents an absolute number of bases for which
  two integrations are considered distinct. If the threshold is set to 3
  it means, provided fields `chr` and `strand` are the same,
  integrations sites which have at least 3 bases in between them are
  considered distinct.

- is_identity_tags:

  Character vector of tags that identify the integration event as
  distinct (except for `"locus"`). See details.

- keep_criteria:

  While scanning, which integration should be kept? The 2 possible
  choices for this parameter are:

  - "max_value": keep the integration site which has the highest value
    (and collapse other values on that integration).

  - "keep_first": keeps the first integration

- value_columns:

  Character vector, contains the names of the numeric experimental
  columns

- max_value_column:

  The column that has to be considered for searching the maximum value

- sample_id_column:

  The name of the column containing the sample identifier

- additional_agg_lambda:

  A named list containing aggregating functions for additional columns.
  See details.

- max_workers:

  Maximum parallel workers allowed

- map_as_file:

  Produce recalibration map as a .tsv file?

- file_path:

  String representing the path were the file will be saved. Must be a
  folder. Relevant only if `map_as_file` is `TRUE`.

- strand_specific:

  **\[deprecated\]** Deprecated, use `is_identity_tags`

## Value

An integration matrix with same or less number of rows

## Details

### The concept of "near"

An integration event is uniquely identified by all fields specified in
the
[`mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/mandatory_IS_vars.md)
look-up table. It can happen to find IS that are formally distinct
(different combination of values in the fields), but that should not
considered distinct in practice, since they represent the same
integration event - this may be due to artefacts at the putative locus
of the IS in the merging of multiple sequencing libraries.

We say that an integration event IS1 is near to another integration
event IS2 if the absolute difference of their loci is strictly lower
than the set `threshold`.

### The IS identity

There is also another aspect to be considered. Since the algorithm is
based on a sliding window mechanism, on which groups of IS should we set
and slide the window?

By default, we have 3 fields in the
[`mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/mandatory_IS_vars.md):
chr, integration_locus, strand, and we assume that all the fields
contribute to the identity of the IS. This means that IS1 and IS2 can be
compared only if they have the same chromosome and the same strand.
However, if we would like to exclude the strand of the integration from
our considerations then IS1 and IS2 can be selected from all the events
that fall on the same chromosome. A practical example:

IS1 = `(chr = "1", strand = "+", integration_locus = 14568)`

IS2 = `(chr = "1", strand = "-", integration_locus = 14567)`

if `is_identity_tags = c("chromosome", "is_strand")` IS1 and IS2 are
considered distinct because they differ in strand, therefore no
correction will be applied to loci of either of the 2. If
`is_identity_tags = c("chromosome")` then IS1 and IS2 are considered
near, because the strand is irrelevant, hence one of the 2 IS will
change locus.

### Aggregating near IS

IS that fall in the same interval are evaluated according to the
criterion selected - if recalibration is necessary, rows with the same
sample ID are aggregated in a single row with a quantification value
that is the sum of all the merged rows.

If the input integration matrix contains annotation columns, that is
additional columns that are not

- part of the mandatory IS vars (see
  [`mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/mandatory_IS_vars.md))

- part of the annotation IS vars (see
  [`annotation_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/mandatory_IS_vars.md))

- the sample identifier column

- the quantification column

it is possible to specify how they should be aggregated. Defaults are
provided for each column type (character, integer, numeric...), but
custom functions can be specified as a named list, where names are
column names in `x` and values are functions to be applied. NOTE:
functions must be purrr-style lambdas and they must perform some kind of
aggregating operation, aka they must take a vector as input and return a
single value. The type of the output should match the type of the target
column. If you specify custom lambdas, provide defaults in the special
element `.defaults`. Example:

    list(
      numeric_col = ~ sum(.x),
      char_col = ~ paste0(.x, collapse = ", "),
      .defaults = default_rec_agg_lambdas()
    )

## Note

We do recommend to use this function in combination with
[comparison_matrix](https://calabrialab.github.io/ISAnalytics/dev/reference/comparison_matrix.md)
to automatically perform re-calibration on all quantification matrices.

## Required tags

The function will explicitly check for the presence of these tags:

- chromosome

- locus

- is_strand

- gene_symbol

## See also

Other Data cleaning and pre-processing:
[`aggregate_metadata()`](https://calabrialab.github.io/ISAnalytics/dev/reference/aggregate_metadata.md),
[`aggregate_values_by_key()`](https://calabrialab.github.io/ISAnalytics/dev/reference/aggregate_values_by_key.md),
[`default_meta_agg()`](https://calabrialab.github.io/ISAnalytics/dev/reference/default_meta_agg.md),
[`outlier_filter()`](https://calabrialab.github.io/ISAnalytics/dev/reference/outlier_filter.md),
[`outliers_by_pool_fragments()`](https://calabrialab.github.io/ISAnalytics/dev/reference/outliers_by_pool_fragments.md),
[`purity_filter()`](https://calabrialab.github.io/ISAnalytics/dev/reference/purity_filter.md),
[`realign_after_collisions()`](https://calabrialab.github.io/ISAnalytics/dev/reference/realign_after_collisions.md),
[`remove_collisions()`](https://calabrialab.github.io/ISAnalytics/dev/reference/remove_collisions.md),
[`threshold_filter()`](https://calabrialab.github.io/ISAnalytics/dev/reference/threshold_filter.md)

## Examples

``` r
data("integration_matrices", package = "ISAnalytics")
rec <- compute_near_integrations(
    x = integration_matrices, map_as_file = FALSE
)
head(rec)
#> # A tibble: 6 × 8
#>   chr   integration_locus strand GeneName GeneStrand CompleteAmplificationID    
#>   <chr>             <dbl> <chr>  <chr>    <chr>      <chr>                      
#> 1 1               8607357 +      RERE     -          PJ01_POOL01_LTR27LC94_PT00…
#> 2 1               8607357 +      RERE     -          PJ01_POOL01_LTR83LC46_PT00…
#> 3 1               8607357 +      RERE     -          PJ01_POOL01_LTR83LC66_PT00…
#> 4 1               8607357 +      RERE     -          PJ01_POOL01_LTR53LC32_PT00…
#> 5 1               8607357 +      RERE     -          PJ01_POOL02_LTR87LC74_PT00…
#> 6 1               8850362 +      RERE     -          PJ01_POOL03_LTR51LC86_PT00…
#> # ℹ 2 more variables: seqCount <dbl>, fragmentEstimate <dbl>
```
