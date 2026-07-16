# Defaults for column aggregations in `compute_near_integrations()`.

Defaults for column aggregations in
[`compute_near_integrations()`](https://calabrialab.github.io/ISAnalytics/dev/reference/compute_near_integrations.md).

## Usage

``` r
default_rec_agg_lambdas()
```

## Value

A named list of lambdas

## Examples

``` r
default_rec_agg_lambdas()
#> $character
#> ~paste0(.x, collapse = ";")
#> <environment: 0x563aae5f5b30>
#> 
#> $integer
#> ~sum(.x, na.rm = TRUE)
#> <environment: 0x563aae5f5b30>
#> 
#> $double
#> ~sum(.x, na.rm = TRUE)
#> <environment: 0x563aae5f5b30>
#> 
#> $logical
#> ~all(.x)
#> <environment: 0x563aae5f5b30>
#> 
```
