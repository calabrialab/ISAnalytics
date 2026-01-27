# Resets dynamic vars to the default values.

Reverts all changes to dynamic vars to the default values. For more
details, refer to the dedicated vignette
[`vignette("workflow_start", package="ISAnalytics")`](https://calabrialab.github.io/ISAnalytics/dev/articles/workflow_start.md).

- `reset_mandatory_IS_vars()` re-sets the look-up table for mandatory IS
  vars.

&nbsp;

- `reset_annotation_IS_vars()` re-sets the look-up table for genomic
  annotation IS vars.

&nbsp;

- `reset_af_columns_def()` re-sets the look-up table for association
  file columns vars

&nbsp;

- `reset_iss_stats_specs()` re-sets the look-up table for VISPA2 pool
  statistics vars

&nbsp;

- `reset_matrix_file_suffixes()` re-sets the matrix file suffixes
  look-up table

&nbsp;

- `reset_dyn_vars_config()` re-sets all look-up tables

## Usage

``` r
reset_mandatory_IS_vars()

reset_annotation_IS_vars()

reset_af_columns_def()

reset_iss_stats_specs()

reset_matrix_file_suffixes()

reset_dyn_vars_config()
```

## Value

`NULL`

## See also

Other dynamic vars:
[`inspect_tags()`](https://calabrialab.github.io/ISAnalytics/dev/reference/inspect_tags.md),
[`mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/mandatory_IS_vars.md),
[`pcr_id_column()`](https://calabrialab.github.io/ISAnalytics/dev/reference/pcr_id_column.md),
[`set_mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/set_mandatory_IS_vars.md),
[`set_matrix_file_suffixes()`](https://calabrialab.github.io/ISAnalytics/dev/reference/set_matrix_file_suffixes.md)

## Examples

``` r
reset_mandatory_IS_vars()
#> Mandatory IS vars reset to default

reset_annotation_IS_vars()
#> Annotation IS vars reset to default

reset_af_columns_def()
#> Association file columns specs reset to default

reset_iss_stats_specs()
#> ISS stats specs reset to default

reset_matrix_file_suffixes()
#> Matrix suffixes specs reset to default

reset_dyn_vars_config()
#> Mandatory IS vars reset to default
#> Annotation IS vars reset to default
#> Association file columns specs reset to default
#> ISS stats specs reset to default
#> Matrix suffixes specs reset to default
```
