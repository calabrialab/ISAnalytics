# Sharing of integration sites between given groups.

**\[stable\]** Computes the amount of integration sites shared between
the groups identified in the input data.

## Usage

``` r
is_sharing(
  ...,
  group_key = c("SubjectID", "CellMarker", "Tissue", "TimePoint"),
  group_keys = NULL,
  n_comp = 2,
  is_count = TRUE,
  relative_is_sharing = TRUE,
  minimal = TRUE,
  include_self_comp = FALSE,
  keep_genomic_coord = FALSE,
  table_for_venn = FALSE
)
```

## Arguments

- ...:

  One or more integration matrices

- group_key:

  Character vector of column names which identify a single group. An
  associated group id will be derived by concatenating the values of
  these fields, separated by "\_"

- group_keys:

  A list of keys for asymmetric grouping. If not NULL the argument
  `group_key` is ignored

- n_comp:

  Number of comparisons to compute. This argument is relevant only if
  provided a single data frame and a single key.

- is_count:

  Logical, if `TRUE` returns also the count of IS for each group and the
  count for the union set

- relative_is_sharing:

  Logical, if `TRUE` also returns the relative sharing.

- minimal:

  Compute only combinations instead of all possible permutations? If
  `TRUE` saves time and excludes redundant comparisons.

- include_self_comp:

  Include comparisons with the same group?

- keep_genomic_coord:

  If `TRUE` keeps the genomic coordinates of the shared integration
  sites in a dedicated column (as a nested table)

- table_for_venn:

  Add column with truth tables for venn plots?

## Value

A data frame

## Details

An integration site is always identified by the combination of fields in
[`mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/mandatory_IS_vars.md),
thus these columns must be present in the input(s).

The function accepts multiple inputs for different scenarios, please
refer to the vignette
[`vignette("workflow_start", package = "ISAnalytics")`](https://calabrialab.github.io/ISAnalytics/dev/articles/workflow_start.md)
for a more in-depth explanation.

### Output

The function outputs a single data frame containing all requested
comparisons and optionally individual group counts, genomic coordinates
of the shared integration sites and truth tables for plotting venn
diagrams.

### Plotting sharing

The sharing data obtained can be easily plotted in a heatmap via the
function
[`sharing_heatmap`](https://calabrialab.github.io/ISAnalytics/dev/reference/sharing_heatmap.md)
or via the function
[`sharing_venn`](https://calabrialab.github.io/ISAnalytics/dev/reference/sharing_venn.md)

## Required tags

The function will explicitly check for the presence of these tags:

- All columns declared in
  [`mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/mandatory_IS_vars.md)

## See also

Other Analysis functions:
[`CIS_grubbs()`](https://calabrialab.github.io/ISAnalytics/dev/reference/CIS_grubbs.md),
[`HSC_population_size_estimate()`](https://calabrialab.github.io/ISAnalytics/dev/reference/HSC_population_size_estimate.md),
[`compute_abundance()`](https://calabrialab.github.io/ISAnalytics/dev/reference/compute_abundance.md),
[`cumulative_is()`](https://calabrialab.github.io/ISAnalytics/dev/reference/cumulative_is.md),
[`gene_frequency_fisher()`](https://calabrialab.github.io/ISAnalytics/dev/reference/gene_frequency_fisher.md),
[`iss_source()`](https://calabrialab.github.io/ISAnalytics/dev/reference/iss_source.md),
[`sample_statistics()`](https://calabrialab.github.io/ISAnalytics/dev/reference/sample_statistics.md),
[`top_integrations()`](https://calabrialab.github.io/ISAnalytics/dev/reference/top_integrations.md),
[`top_targeted_genes()`](https://calabrialab.github.io/ISAnalytics/dev/reference/top_targeted_genes.md)

## Examples

``` r
data("integration_matrices", package = "ISAnalytics")
data("association_file", package = "ISAnalytics")
aggreg <- aggregate_values_by_key(
    x = integration_matrices,
    association_file = association_file,
    value_cols = c("seqCount", "fragmentEstimate")
)
sharing <- is_sharing(aggreg)
sharing
#> # A tibble: 190 × 9
#>    g1            g2    shared count_g1 count_g2 count_union on_g1 on_g2 on_union
#>    <chr>         <chr>  <int>    <int>    <int>       <int> <dbl> <dbl>    <dbl>
#>  1 PT001_MNC_BM… PT00…     21       54      114         147 38.9   18.4    14.3 
#>  2 PT001_MNC_BM… PT00…     24       54       59          89 44.4   40.7    27.0 
#>  3 PT001_MNC_BM… PT00…     15       54       89         128 27.8   16.9    11.7 
#>  4 PT001_MNC_BM… PT00…     21       54       78         111 38.9   26.9    18.9 
#>  5 PT001_MNC_BM… PT00…      7       54       28          75 13.0   25       9.33
#>  6 PT001_MNC_BM… PT00…      8       54       59         105 14.8   13.6     7.62
#>  7 PT001_MNC_BM… PT00…     20       54       48          82 37.0   41.7    24.4 
#>  8 PT001_MNC_BM… PT00…      4       54       29          79  7.41  13.8     5.06
#>  9 PT001_MNC_BM… PT00…     15       54       43          82 27.8   34.9    18.3 
#> 10 PT001_MNC_BM… PT00…      0       54       98         152  0      0       0   
#> # ℹ 180 more rows
```
