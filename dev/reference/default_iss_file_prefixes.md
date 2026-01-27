# Default regex prefixes for Vispa2 stats files.

Note that each element is a regular expression.

## Usage

``` r
default_iss_file_prefixes()
```

## Value

A character vector of regexes

## See also

Other Import functions helpers:
[`annotation_issues()`](https://calabrialab.github.io/ISAnalytics/dev/reference/annotation_issues.md),
[`date_formats()`](https://calabrialab.github.io/ISAnalytics/dev/reference/date_formats.md),
[`default_af_transform()`](https://calabrialab.github.io/ISAnalytics/dev/reference/default_af_transform.md),
[`matching_options()`](https://calabrialab.github.io/ISAnalytics/dev/reference/matching_options.md),
[`quantification_types()`](https://calabrialab.github.io/ISAnalytics/dev/reference/quantification_types.md)

## Examples

``` r
default_iss_file_prefixes()
#> [1] "stats\\.sequence." "stats\\.matrix."  
```
