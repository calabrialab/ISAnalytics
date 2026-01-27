# Possible choices for the `dates_format` parameter in `import_association_file`, `import_parallel_vispa2Matrices_interactive` and `import_parallel_vispa2Matrices_auto`.

All options correspond to `lubridate` functions, see more in the
dedicated package documentation.

## Usage

``` r
date_formats()
```

## Value

A character vector

## See also

[`import_association_file`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_association_file.md),
[`import_parallel_Vispa2Matrices_auto`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_parallel_Vispa2Matrices_auto.md)

Other Import functions helpers:
[`annotation_issues()`](https://calabrialab.github.io/ISAnalytics/dev/reference/annotation_issues.md),
[`default_af_transform()`](https://calabrialab.github.io/ISAnalytics/dev/reference/default_af_transform.md),
[`default_iss_file_prefixes()`](https://calabrialab.github.io/ISAnalytics/dev/reference/default_iss_file_prefixes.md),
[`matching_options()`](https://calabrialab.github.io/ISAnalytics/dev/reference/matching_options.md),
[`quantification_types()`](https://calabrialab.github.io/ISAnalytics/dev/reference/quantification_types.md)

## Examples

``` r
date_formats()
#>  [1] "ymd"     "ydm"     "mdy"     "myd"     "dmy"     "dym"     "yq"     
#>  [8] "ym"      "my"      "ymd_hms" "ymd_hm"  "ymd_h"   "dmy_hms" "dmy_hm" 
#> [15] "dmy_h"   "mdy_hms" "mdy_hm"  "mdy_h"   "ydm_hms" "ydm_hm"  "ydm_h"  
```
