# Import integration matrices from paths in the association file.

**\[stable\]** The function offers a convenient way of importing
multiple integration matrices in an automated or semi-automated way. For
more details see the "How to use import functions" vignette:
[`vignette("workflow_start", package = "ISAnalytics")`](https://calabrialab.github.io/ISAnalytics/dev/articles/workflow_start.md)

## Usage

``` r
import_parallel_Vispa2Matrices(
  association_file,
  quantification_type = c("seqCount", "fragmentEstimate"),
  matrix_type = c("annotated", "not_annotated"),
  workers = 2,
  multi_quant_matrix = TRUE,
  report_path = default_report_path(),
  patterns = NULL,
  matching_opt = matching_options(),
  mode = "AUTO",
  ...
)
```

## Arguments

- association_file:

  Data frame imported via
  [import_association_file](https://calabrialab.github.io/ISAnalytics/dev/reference/import_association_file.md)
  (with file system alignment)

- quantification_type:

  A vector of requested quantification_types. Possible choices are
  [quantification_types](https://calabrialab.github.io/ISAnalytics/dev/reference/quantification_types.md)

- matrix_type:

  A single string representing the type of matrices to be imported. Can
  only be one in "annotated" or "not_annotated".

- workers:

  A single integer representing the number of parallel workers to use
  for the import

- multi_quant_matrix:

  If set to `TRUE` will produce a multi-quantification matrix through
  [comparison_matrix](https://calabrialab.github.io/ISAnalytics/dev/reference/comparison_matrix.md)
  instead of a list.

- report_path:

  The path where the report file should be saved. Can be a folder or
  `NULL` if no report should be produced. Defaults to
  `{user_home}/ISAnalytics_reports`.

- patterns:

  A character vector of additional patterns to match on file names.
  Please note that patterns must be regular expressions. Can be `NULL`
  if no patterns need to be matched.

- matching_opt:

  A single value between
  [matching_options](https://calabrialab.github.io/ISAnalytics/dev/reference/matching_options.md)

- mode:

  Only `AUTO` is supported. As of `ISAnalytics 1.8.3`, the value
  `INTERACTIVE` is officially deprecated.

- ...:

  \<[`dynamic-dots`](https://rlang.r-lib.org/reference/dyn-dots.html)\>
  Additional named arguments to pass to `comparison_matrix` and
  `import_single_Vispa2_matrix`

## Value

Either a multi-quantification matrix or a list of integration matrices

## Required tags

The function will explicitly check for the presence of these tags:

- project_id

- vispa_concatenate

## See also

Other Import functions:
[`import_Vispa2_stats()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_Vispa2_stats.md),
[`import_association_file()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_association_file.md),
[`import_single_Vispa2Matrix()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_single_Vispa2Matrix.md)

## Examples

``` r
fs_path <- generate_default_folder_structure(type = "correct")
af <- import_association_file(fs_path$af,
    root = fs_path$root,
    report_path = NULL
)
matrices <- import_parallel_Vispa2Matrices(af,
    c("seqCount", "fragmentEstimate"),
    mode = "AUTO", report_path = NULL
)
head(matrices)
#> # A tibble: 6 × 8
#>   chr   integration_locus strand GeneName     GeneStrand CompleteAmplificationID
#>   <chr>             <int> <chr>  <chr>        <chr>      <chr>                  
#> 1 16             68164148 +      NFATC3       +          PJ01_POOL01_LTR75LC38_…
#> 2 4             129390130 +      LOC100507487 +          PJ01_POOL01_LTR75LC38_…
#> 3 5              84009671 -      EDIL3        -          PJ01_POOL01_LTR75LC38_…
#> 4 12             54635693 -      CBX5         -          PJ01_POOL01_LTR75LC38_…
#> 5 2             181930711 +      UBE2E3       +          PJ01_POOL01_LTR75LC38_…
#> 6 20             35920986 +      MANBAL       +          PJ01_POOL01_LTR75LC38_…
#> # ℹ 2 more variables: fragmentEstimate <dbl>, seqCount <int>
```
