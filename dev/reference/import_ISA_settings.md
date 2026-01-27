# Import a dynamic vars settings profile.

The function allows the import of an existing dynamic vars profile in
json format. This is a quick and convenient way to set up the workflow,
alternative to specifying lookup tables manually through the
corresponding setter functions. For more details, refer to the dedicated
vignette
[`vignette("workflow_start", package="ISAnalytics")`](https://calabrialab.github.io/ISAnalytics/dev/articles/workflow_start.md).

## Usage

``` r
import_ISA_settings(path)
```

## Arguments

- path:

  The path to the json file on disk

## Value

`NULL`

## See also

Other Utilities:
[`as_sparse_matrix()`](https://calabrialab.github.io/ISAnalytics/dev/reference/as_sparse_matrix.md),
[`comparison_matrix()`](https://calabrialab.github.io/ISAnalytics/dev/reference/comparison_matrix.md),
[`enable_progress_bars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/enable_progress_bars.md),
[`export_ISA_settings()`](https://calabrialab.github.io/ISAnalytics/dev/reference/export_ISA_settings.md),
[`generate_Vispa2_launch_AF()`](https://calabrialab.github.io/ISAnalytics/dev/reference/generate_Vispa2_launch_AF.md),
[`generate_blank_association_file()`](https://calabrialab.github.io/ISAnalytics/dev/reference/generate_blank_association_file.md),
[`generate_default_folder_structure()`](https://calabrialab.github.io/ISAnalytics/dev/reference/generate_default_folder_structure.md),
[`separate_quant_matrices()`](https://calabrialab.github.io/ISAnalytics/dev/reference/separate_quant_matrices.md),
[`transform_columns()`](https://calabrialab.github.io/ISAnalytics/dev/reference/transform_columns.md)

## Examples

``` r
tmp_folder <- tempdir()
export_ISA_settings(tmp_folder, "DEFAULT")
#> Settings profile correctly saved
#> ℹ Saved at: /tmp/RtmpNagWuo/DEFAULT_ISAsettings.json
import_ISA_settings(fs::path(tmp_folder, "DEFAULT_ISAsettings.json"))
#> Mandatory IS vars successfully changed
#> Annotation IS vars successfully changed
#> Association file columns specs successfully changed
#> ISS stats specs successfully changed
#> Matrix suffixes specs successfully changed
reset_dyn_vars_config()
#> Mandatory IS vars reset to default
#> Annotation IS vars reset to default
#> Association file columns specs reset to default
#> ISS stats specs reset to default
#> Matrix suffixes specs reset to default
```
