# Check for genomic annotation problems in IS matrices.

**\[experimental\]** This helper function checks if each individual
integration site, identified by the
[`mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/mandatory_IS_vars.md),
has been annotated with two or more distinct gene symbols.

## Usage

``` r
annotation_issues(matrix)
```

## Arguments

- matrix:

  Either a single matrix or a list of matrices, ideally obtained via
  [`import_parallel_Vispa2Matrices()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_parallel_Vispa2Matrices.md)
  or
  [`import_single_Vispa2Matrix()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_single_Vispa2Matrix.md)

## Value

Either `NULL` if no issues were detected or 1 or more data frames with
genomic coordinates of the IS and the number of distinct genes
associated

## See also

Other Import functions helpers:
[`date_formats()`](https://calabrialab.github.io/ISAnalytics/dev/reference/date_formats.md),
[`default_af_transform()`](https://calabrialab.github.io/ISAnalytics/dev/reference/default_af_transform.md),
[`default_iss_file_prefixes()`](https://calabrialab.github.io/ISAnalytics/dev/reference/default_iss_file_prefixes.md),
[`matching_options()`](https://calabrialab.github.io/ISAnalytics/dev/reference/matching_options.md),
[`quantification_types()`](https://calabrialab.github.io/ISAnalytics/dev/reference/quantification_types.md)

## Examples

``` r
data("integration_matrices", package = "ISAnalytics")
annotation_issues(integration_matrices)
#> No annotation issues found
#> NULL
```
