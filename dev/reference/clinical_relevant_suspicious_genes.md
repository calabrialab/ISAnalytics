# Clinical relevant suspicious genes (for mouse and human).

Clinical relevant suspicious genes (for mouse and human).

## Usage

``` r
clinical_relevant_suspicious_genes()
```

## Value

A data frame

## See also

Other Plotting function helpers:
[`known_clinical_oncogenes()`](https://calabrialab.github.io/ISAnalytics/dev/reference/known_clinical_oncogenes.md)

## Examples

``` r
clinical_relevant_suspicious_genes()
#> # A tibble: 6 × 3
#>   GeneName ClinicalRelevance DOIReference                                
#>   <chr>    <lgl>             <chr>                                       
#> 1 DNMT3A   TRUE              https://doi.org/10.1182/blood-2018-01-829937
#> 2 TET2     TRUE              https://doi.org/10.1182/blood-2018-01-829937
#> 3 ASXL1    TRUE              https://doi.org/10.1182/blood-2018-01-829937
#> 4 JAK2     TRUE              https://doi.org/10.1182/blood-2018-01-829937
#> 5 CBL      TRUE              https://doi.org/10.1182/blood-2018-01-829937
#> 6 TP53     TRUE              https://doi.org/10.1182/blood-2018-01-829937
```
