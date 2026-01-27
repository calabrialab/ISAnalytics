# Export a dynamic vars settings profile.

This function allows exporting the currently set dynamic vars in json
format so it can be quickly imported later. Dynamic variables need to be
properly set via the setter functions before calling the function. For
more details, refer to the dedicated vignette
[`vignette("workflow_start", package="ISAnalytics")`](https://calabrialab.github.io/ISAnalytics/dev/articles/workflow_start.md).

## Usage

``` r
export_ISA_settings(folder, setting_profile_name)
```

## Arguments

- folder:

  The path to the folder where the file should be saved. If the folder
  doesn't exist, it gets created automatically

- setting_profile_name:

  A name for the settings profile

## Value

`NULL`

## See also

Other Utilities:
[`as_sparse_matrix()`](https://calabrialab.github.io/ISAnalytics/dev/reference/as_sparse_matrix.md),
[`comparison_matrix()`](https://calabrialab.github.io/ISAnalytics/dev/reference/comparison_matrix.md),
[`enable_progress_bars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/enable_progress_bars.md),
[`generate_Vispa2_launch_AF()`](https://calabrialab.github.io/ISAnalytics/dev/reference/generate_Vispa2_launch_AF.md),
[`generate_blank_association_file()`](https://calabrialab.github.io/ISAnalytics/dev/reference/generate_blank_association_file.md),
[`generate_default_folder_structure()`](https://calabrialab.github.io/ISAnalytics/dev/reference/generate_default_folder_structure.md),
[`import_ISA_settings()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_ISA_settings.md),
[`separate_quant_matrices()`](https://calabrialab.github.io/ISAnalytics/dev/reference/separate_quant_matrices.md),
[`transform_columns()`](https://calabrialab.github.io/ISAnalytics/dev/reference/transform_columns.md)

## Examples

``` r
tmp_folder <- tempdir()
export_ISA_settings(tmp_folder, "DEFAULT")
#> Settings profile correctly saved
#> ℹ Saved at: /tmp/RtmpNagWuo/DEFAULT_ISAsettings.json
```
