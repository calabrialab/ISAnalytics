# Enable global progress bars for ISAnalytics functions.

This is a simple wrapper around functions from the package `progressr`.
To customize the appearance of the progress bar, please refer to
[progressr](https://progressr.futureverse.org/) documentation.

## Usage

``` r
enable_progress_bars()
```

## Value

`NULL`

## See also

Other Utilities:
[`as_sparse_matrix()`](https://calabrialab.github.io/ISAnalytics/dev/reference/as_sparse_matrix.md),
[`comparison_matrix()`](https://calabrialab.github.io/ISAnalytics/dev/reference/comparison_matrix.md),
[`export_ISA_settings()`](https://calabrialab.github.io/ISAnalytics/dev/reference/export_ISA_settings.md),
[`generate_Vispa2_launch_AF()`](https://calabrialab.github.io/ISAnalytics/dev/reference/generate_Vispa2_launch_AF.md),
[`generate_blank_association_file()`](https://calabrialab.github.io/ISAnalytics/dev/reference/generate_blank_association_file.md),
[`generate_default_folder_structure()`](https://calabrialab.github.io/ISAnalytics/dev/reference/generate_default_folder_structure.md),
[`import_ISA_settings()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_ISA_settings.md),
[`separate_quant_matrices()`](https://calabrialab.github.io/ISAnalytics/dev/reference/separate_quant_matrices.md),
[`transform_columns()`](https://calabrialab.github.io/ISAnalytics/dev/reference/transform_columns.md)

## Examples

``` r
enable_progress_bars()
#> Error in globalCallingHandlers(condition = global_progression_handler): should not be called with handlers on the stack
progressr::handlers(global = FALSE) # Deactivate
#> Error in globalCallingHandlers(NULL): should not be called with handlers on the stack
```
