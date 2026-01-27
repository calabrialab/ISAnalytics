# Compute CIS and Grubbs test over different time points and groups.

**\[experimental\]** Computes common insertion sites and Grubbs test for
each separate group and separating different time points among the same
group. The logic applied is the same as the function
[`CIS_grubbs()`](https://calabrialab.github.io/ISAnalytics/dev/reference/CIS_grubbs.md).

## Usage

``` r
CIS_grubbs_overtime(
  x,
  genomic_annotation_file = "hg19",
  grubbs_flanking_gene_bp = 1e+05,
  threshold_alpha = 0.05,
  group = "SubjectID",
  timepoint_col = "TimePoint",
  as_df = TRUE,
  return_missing_as_df = TRUE,
  max_workers = NULL
)
```

## Arguments

- x:

  An integration matrix, must include the
  [`mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/mandatory_IS_vars.md)
  columns and the
  [`annotation_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/mandatory_IS_vars.md)
  columns

- genomic_annotation_file:

  Database file for gene annotation, see details.

- grubbs_flanking_gene_bp:

  Number of base pairs flanking a gene

- threshold_alpha:

  Significance threshold

- group:

  A character vector of column names that identifies a group. Each group
  must contain one or more time points.

- timepoint_col:

  What is the name of the column containing time points?

- as_df:

  Choose the result format: if `TRUE` the results are returned as a
  single data frame containing a column for the group id and a column
  for the time point, if `FALSE` results are returned in the form of
  nested lists (one table for each time point and for each group), if
  `"group"` results are returned as a list separated for each group but
  containing a single table with all time points.

- return_missing_as_df:

  Returns those genes present in the input df but not in the refgenes as
  a data frame?

- max_workers:

  Maximum number of parallel workers. If `NULL` the maximum number of
  workers is calculated automatically.

## Value

A list with results and optionally missing genes info

## Details

### Genomic annotation file

A data frame containing genes annotation for the specific genome. From
version `1.5.4` the argument `genomic_annotation_file` accepts only data
frames or package provided defaults. The user is responsible for
importing the appropriate tabular files if customization is needed. The
annotations for the human genome (hg19 or hg38) and murine genome (mm9
or mm10) are already included in this package: to use one of them just
set the argument `genomic_annotation_file` to either `"hg19"`, `"hg38"`,
`"mm9"` or `"mm10"`. If for any reason the user is performing an
analysis on another genome, this file needs to be changed respecting the
USCS Genome Browser format, meaning the input file headers should
include:

name2, chrom, strand, min_txStart, max_txEnd, minmax_TxLen,
average_TxLen, name, min_cdsStart, max_cdsEnd, minmax_CdsLen,
average_CdsLen

## Examples

``` r
data("integration_matrices", package = "ISAnalytics")
data("association_file", package = "ISAnalytics")
aggreg <- aggregate_values_by_key(
    x = integration_matrices,
    association_file = association_file,
    value_cols = c("seqCount", "fragmentEstimate")
)
cis_overtime <- CIS_grubbs_overtime(aggreg)
#> Warning: Warning: missing genes in refgenes table
#> ℹ A total of 5 genes were found in the input data but not in the refgene table. This may be caused by a mismatch in the annotation phase of the matrix. Here is a summary: 
#> # A tibble: 5 × 3
#>   GeneName  GeneStrand chr  
#>   <chr>     <chr>      <chr>
#> 1 PLEKHG4B  -          14   
#> 2 CRELD2    -          15   
#> 3 UBE2D2    +          16   
#> 4 LINC01133 +          19   
#> 5 HTR4      +          6    
#> ℹ NOTE: missing genes will be removed from the final output! Review results carefully
#> ℹ A total of 25 IS will be removed because of missing genes ( 2.33 % of total IS in input)
cis_overtime
#> $missing_genes
#> # A tibble: 5 × 3
#>   GeneName  GeneStrand chr  
#>   <chr>     <chr>      <chr>
#> 1 PLEKHG4B  -          14   
#> 2 CRELD2    -          15   
#> 3 UBE2D2    +          16   
#> 4 LINC01133 +          19   
#> 5 HTR4      +          6    
#> 
#> $missing_is
#> $missing_is$absolute
#> [1] 25
#> 
#> $missing_is$perc
#> [1] 2.33
#> 
#> 
#> $cis
#> # A tibble: 932 × 39
#>    GeneName GeneStrand chr       n      mean    sd   median trimmed   mad    min
#>    <chr>    <chr>      <chr> <int>     <dbl> <dbl>    <dbl>   <dbl> <dbl>  <dbl>
#>  1 ABHD2    +          15        1  89642998    NA   8.96e7  8.96e7     0 8.96e7
#>  2 ACAP2    -          3         1 195107772    NA   1.95e8  1.95e8     0 1.95e8
#>  3 ACOX1    -          17        1  73953431    NA   7.40e7  7.40e7     0 7.40e7
#>  4 ADD1     +          4         1   2859237    NA   2.86e6  2.86e6     0 2.86e6
#>  5 ANKEF1   +          20        1  10028583    NA   1.00e7  1.00e7     0 1.00e7
#>  6 ANKRD52  -          12        1  56647382    NA   5.66e7  5.66e7     0 5.66e7
#>  7 ASAP1    -          8         1 131310762    NA   1.31e8  1.31e8     0 1.31e8
#>  8 ATF1     +          12        1  51189341    NA   5.12e7  5.12e7     0 5.12e7
#>  9 ATF7     -          12        1  53987984    NA   5.40e7  5.40e7     0 5.40e7
#> 10 C6orf106 -          6         1  34643226    NA   3.46e7  3.46e7     0 3.46e7
#> # ℹ 922 more rows
#> # ℹ 29 more variables: max <dbl>, range <dbl>, skew <dbl>, kurtosis <dbl>,
#> #   n_IS_perGene <int>, min_bp_integration_locus <dbl>,
#> #   max_bp_integration_locus <dbl>, IS_span_bp <dbl>,
#> #   avg_bp_integration_locus <dbl>, median_bp_integration_locus <dbl>,
#> #   distinct_orientations <int>, average_TxLen <dbl>,
#> #   raw_gene_integration_frequency <dbl>, …
#> 
```
