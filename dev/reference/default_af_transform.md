# Default transformations to apply to association file columns.

A list of default transformations to apply to the association file
columns after importing it via
[`import_association_file()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_association_file.md)

## Usage

``` r
default_af_transform(convert_tp)
```

## Arguments

- convert_tp:

  The value of the argument `convert_tp` in the call to
  [`import_association_file()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_association_file.md)

## Value

A named list of lambdas

## See also

Other Import functions helpers:
[`annotation_issues()`](https://calabrialab.github.io/ISAnalytics/dev/reference/annotation_issues.md),
[`date_formats()`](https://calabrialab.github.io/ISAnalytics/dev/reference/date_formats.md),
[`default_iss_file_prefixes()`](https://calabrialab.github.io/ISAnalytics/dev/reference/default_iss_file_prefixes.md),
[`matching_options()`](https://calabrialab.github.io/ISAnalytics/dev/reference/matching_options.md),
[`quantification_types()`](https://calabrialab.github.io/ISAnalytics/dev/reference/quantification_types.md)

## Examples

``` r
default_af_transform(TRUE)
#> $TimepointMonths
#> ~dplyr::if_else(is.na(.x), NA_character_, stringr::str_pad(as.character(.x), 
#>     pad = "0", side = "left", width = max(nchar(as.character(.x[!is.na(.x)])), 
#>         na.rm = TRUE) + 1))
#> <environment: 0x55a17957d608>
#> 
#> $TimepointYears
#> ~dplyr::if_else(is.na(.x), NA_character_, stringr::str_pad(as.character(.x), 
#>     pad = "0", side = "left", width = max(nchar(as.character(.x[!is.na(.x)])), 
#>         na.rm = TRUE) + 1))
#> <environment: 0x55a17957d608>
#> 
```
