# A set of pre-defined functions for `sample_statistics`.

A set of pre-defined functions for `sample_statistics`.

## Usage

``` r
default_stats()
```

## Value

A named list of functions/purrr-style lambdas

## Examples

``` r
default_stats()
#> $sum
#> ~sum(.x, na.rm = TRUE)
#> <environment: 0x55cf393bca78>
#> 
#> $count
#> function (x)  .Primitive("length")
#> 
#> $shannon
#> ~vegan::diversity(.x, index = "shannon")
#> <environment: 0x55cf393bca78>
#> 
#> $simpson
#> ~vegan::diversity(.x, index = "simpson")
#> <environment: 0x55cf393bca78>
#> 
#> $invsimpson
#> ~vegan::diversity(.x, index = "invsimpson")
#> <environment: 0x55cf393bca78>
#> 
#> $describe
#> ~tibble::as_tibble(psych::describe(.x))
#> <environment: 0x55cf393bca78>
#> 
```
