# Top n targeted genes based on number of IS.

**\[experimental\]** Produces a summary of the number of integration
events per gene, orders the table in decreasing order and slices the
first n rows - either on all the data frame or by group.

## Usage

``` r
top_targeted_genes(
  x,
  n = 20,
  key = c("SubjectID", "CellMarker", "Tissue", "TimePoint"),
  consider_chr = TRUE,
  consider_gene_strand = TRUE,
  as_df = TRUE
)
```

## Arguments

- x:

  An integration matrix - must be annotated

- n:

  Number of rows to slice

- key:

  If slice has to be performed for each group, the character vector of
  column names that identify the groups. If `NULL` considers the whole
  input data frame.

- consider_chr:

  Logical, should the chromosome be taken into account? See details.

- consider_gene_strand:

  Logical, should the gene strand be taken into account? See details.

- as_df:

  If computation is performed by group, `TRUE` returns all groups merged
  in a single data frame with a column containing the group id. If
  `FALSE` returns a named list.

## Value

A data frame or a list of data frames

## Details

### Gene grouping

When producing a summary of IS by gene, there are different options that
can be chosen. The argument `consider_chr` accounts for the fact that
some genes (same gene symbol) may span more than one chromosome: if set
to `TRUE` counts of IS will be separated for those genes that span 2 or
more chromosomes - in other words they will be in 2 different rows of
the output table. On the contrary, if the argument is set to `FALSE`,
counts will be produced in a single row.

NOTE: the function counts **DISTINCT** integration events, which
logically corresponds to a union of sets. Be aware of the fact that
counts per group and counts with different arguments might be different:
if for example counts are performed by considering chromosome and there
is one gene symbol with 2 different counts, the sum of those 2 will
likely not be equal to the count obtained by performing the calculations
without considering the chromosome.

The same reasoning can be applied for the argument
`consider_gene_strand`, that takes into account the strand of the gene.

## Required tags

The function will explicitly check for the presence of these tags:

- chromosome

- locus

- gene_symbol

- gene_strand

Note that the tags "gene_strand" and "chromosome" are explicitly
required only if `consider_chr = TRUE` and/or
`consider_gene_strand = TRUE`.

## See also

Other Analysis functions:
[`CIS_grubbs()`](https://calabrialab.github.io/ISAnalytics/dev/reference/CIS_grubbs.md),
[`HSC_population_size_estimate()`](https://calabrialab.github.io/ISAnalytics/dev/reference/HSC_population_size_estimate.md),
[`compute_abundance()`](https://calabrialab.github.io/ISAnalytics/dev/reference/compute_abundance.md),
[`cumulative_is()`](https://calabrialab.github.io/ISAnalytics/dev/reference/cumulative_is.md),
[`gene_frequency_fisher()`](https://calabrialab.github.io/ISAnalytics/dev/reference/gene_frequency_fisher.md),
[`is_sharing()`](https://calabrialab.github.io/ISAnalytics/dev/reference/is_sharing.md),
[`iss_source()`](https://calabrialab.github.io/ISAnalytics/dev/reference/iss_source.md),
[`sample_statistics()`](https://calabrialab.github.io/ISAnalytics/dev/reference/sample_statistics.md),
[`top_integrations()`](https://calabrialab.github.io/ISAnalytics/dev/reference/top_integrations.md)

## Examples

``` r
data("integration_matrices", package = "ISAnalytics")
top_targ <- top_targeted_genes(
    integration_matrices,
    key = NULL
)
top_targ
#> # A tibble: 20 × 4
#>    GeneName GeneStrand chr    n_IS
#>    <chr>    <chr>      <chr> <int>
#>  1 KMT5B    -          11        5
#>  2 RERE     -          1         4
#>  3 ACOX1    -          17        3
#>  4 ADD1     +          4         3
#>  5 ANKFY1   -          17        3
#>  6 R3HDM2   -          12        3
#>  7 SPATS2   +          12        3
#>  8 ATF1     +          12        2
#>  9 ATF7     -          12        2
#> 10 AXIN1    -          16        2
#> 11 C6orf10  -          6         2
#> 12 CD163    -          12        2
#> 13 CDH2     -          18        2
#> 14 CLECL1   -          12        2
#> 15 CNOT6L   -          4         2
#> 16 DNMT1    -          19        2
#> 17 FCHSD2   -          11        2
#> 18 GRB2     -          17        2
#> 19 HBE1     -          11        2
#> 20 HNRNPM   +          19        2
```
