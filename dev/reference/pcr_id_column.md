# Easily retrieve the name of the pcr id column.

The function is a shortcut to retrieve the currently set pcr id column
name from the association file column tags look-up table. This column is
needed every time a joining operation with metadata needs to be
performed

## Usage

``` r
pcr_id_column()
```

## Value

The name of the column

## See also

Other dynamic vars:
[`inspect_tags()`](https://calabrialab.github.io/ISAnalytics/dev/reference/inspect_tags.md),
[`mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/mandatory_IS_vars.md),
[`reset_mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/reset_mandatory_IS_vars.md),
[`set_mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/set_mandatory_IS_vars.md),
[`set_matrix_file_suffixes()`](https://calabrialab.github.io/ISAnalytics/dev/reference/set_matrix_file_suffixes.md)

## Examples

``` r
pcr_id_column()
#> [1] "CompleteAmplificationID"
```
