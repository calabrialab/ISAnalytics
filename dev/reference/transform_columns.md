# Apply transformations to an arbitrary number of columns.

This function takes a named list of purr-style lambdas where names are
the names of the columns in the data frame that must be transformed.
NOTE: the columns are overridden, not appended.

## Usage

``` r
transform_columns(df, transf_list)
```

## Arguments

- df:

  The data frame on which transformations should be operated

- transf_list:

  A named list of purrr-style lambdas, where names are column names the
  function should be applied to.

## Value

A data frame with transformed columns

## Details

Lambdas provided in input must be transformations, aka functions that
take in input a vector and return a vector of the same length as the
input.

If the input transformation list contains column names that are not
present in the input data frame, they are simply ignored.

## See also

Other Utilities:
[`as_sparse_matrix()`](https://calabrialab.github.io/ISAnalytics/dev/reference/as_sparse_matrix.md),
[`comparison_matrix()`](https://calabrialab.github.io/ISAnalytics/dev/reference/comparison_matrix.md),
[`enable_progress_bars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/enable_progress_bars.md),
[`export_ISA_settings()`](https://calabrialab.github.io/ISAnalytics/dev/reference/export_ISA_settings.md),
[`generate_Vispa2_launch_AF()`](https://calabrialab.github.io/ISAnalytics/dev/reference/generate_Vispa2_launch_AF.md),
[`generate_blank_association_file()`](https://calabrialab.github.io/ISAnalytics/dev/reference/generate_blank_association_file.md),
[`generate_default_folder_structure()`](https://calabrialab.github.io/ISAnalytics/dev/reference/generate_default_folder_structure.md),
[`import_ISA_settings()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_ISA_settings.md),
[`separate_quant_matrices()`](https://calabrialab.github.io/ISAnalytics/dev/reference/separate_quant_matrices.md)

## Examples

``` r
df <- tibble::tribble(
    ~A, ~B, ~C, ~D,
    1, 2, "a", "aa",
    3, 4, "b", "bb",
    5, 6, "c", "cc"
)
lambdas <- list(A = ~ .x + 1, B = ~ .x + 2, C = ~ stringr::str_to_upper(.x))
transform_columns(df, lambdas)
#> # A tibble: 3 × 4
#>       A     B C     D    
#>   <dbl> <dbl> <chr> <chr>
#> 1     2     4 A     aa   
#> 2     4     6 B     bb   
#> 3     6     8 C     cc   
```
