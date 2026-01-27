# Grubbs test for Common Insertion Sites (CIS).

**\[stable\]** Statistical approach for the validation of common
insertion sites significance based on the comparison of the integration
frequency at the CIS gene with respect to other genes contained in the
surrounding genomic regions. For more details please refer to this
paper:
<https://ashpublications.org/blood/article/117/20/5332/21206/Lentiviral-vector-common-integration-sites-in>

## Usage

``` r
CIS_grubbs(
  x,
  genomic_annotation_file = "hg19",
  grubbs_flanking_gene_bp = 1e+05,
  threshold_alpha = 0.05,
  by = NULL,
  return_missing_as_df = TRUE,
  results_as_list = TRUE
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

- by:

  Either `NULL` or a character vector of column names. If not NULL, the
  function will perform calculations for each group and return a list of
  data frames with the results. E.g. for `by = "SubjectID"`, CIS will be
  computed for each distinct SubjectID found in the table ("SubjectID"
  column must be included in the input data frame).

- return_missing_as_df:

  Returns those genes present in the input df but not in the refgenes as
  a data frame?

- results_as_list:

  If `TRUE` return the group computations as a named list, otherwise
  return a single df with an additional column containing the group id

## Value

A data frame

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

## Required tags

The function will explicitly check for the presence of these tags:

- chromosome

- locus

- is_strand

- gene_symbol

- gene_strand

## See also

Other Analysis functions:
[`HSC_population_size_estimate()`](https://calabrialab.github.io/ISAnalytics/dev/reference/HSC_population_size_estimate.md),
[`compute_abundance()`](https://calabrialab.github.io/ISAnalytics/dev/reference/compute_abundance.md),
[`cumulative_is()`](https://calabrialab.github.io/ISAnalytics/dev/reference/cumulative_is.md),
[`gene_frequency_fisher()`](https://calabrialab.github.io/ISAnalytics/dev/reference/gene_frequency_fisher.md),
[`is_sharing()`](https://calabrialab.github.io/ISAnalytics/dev/reference/is_sharing.md),
[`iss_source()`](https://calabrialab.github.io/ISAnalytics/dev/reference/iss_source.md),
[`sample_statistics()`](https://calabrialab.github.io/ISAnalytics/dev/reference/sample_statistics.md),
[`top_integrations()`](https://calabrialab.github.io/ISAnalytics/dev/reference/top_integrations.md),
[`top_targeted_genes()`](https://calabrialab.github.io/ISAnalytics/dev/reference/top_targeted_genes.md)

## Examples

``` r
data("integration_matrices", package = "ISAnalytics")
cis <- CIS_grubbs(integration_matrices)
#> Warning: Warning: missing genes in refgenes table
#> ℹ A total of 5 genes were found in the input data but not in the refgene table. This may be caused by a mismatch in the annotation phase of the matrix. Here is a summary: 
#>     GeneName GeneStrand    chr
#>       <char>     <char> <char>
#> 1:    CRELD2          -     15
#> 2:      HTR4          +      6
#> 3:    UBE2D2          +     16
#> 4: LINC01133          +     19
#> 5:  PLEKHG4B          -     14
#> ℹ NOTE: missing genes will be removed from the final output! Review results carefully
#> ℹ A total of 25 IS will be removed because of missing genes ( 1.48 % of total IS in input)
cis
#> $missing_genes
#>     GeneName GeneStrand    chr
#>       <char>     <char> <char>
#> 1:    CRELD2          -     15
#> 2:      HTR4          +      6
#> 3:    UBE2D2          +     16
#> 4: LINC01133          +     19
#> 5:  PLEKHG4B          -     14
#> 
#> $missing_is
#> $missing_is$absolute
#> [1] 25
#> 
#> $missing_is$perc
#> [1] 1.48
#> 
#> 
#> $cis
#> # A tibble: 501 × 37
#>    GeneName GeneStrand chr       n       mean    sd  median trimmed   mad    min
#>    <chr>    <chr>      <chr> <int>      <dbl> <dbl>   <dbl>   <dbl> <dbl>  <dbl>
#>  1 ABHD16A  -          6         3  31654876     0   3.17e7  3.17e7    0  3.17e7
#>  2 ABHD2    +          15        3  89642998     0   8.96e7  8.96e7    0  8.96e7
#>  3 ACAP2    -          3         2 195107772     0   1.95e8  1.95e8    0  1.95e8
#>  4 ACOX1    -          17       10  73960093. 8844.  7.40e7  7.40e7 4236. 7.40e7
#>  5 ACSM3    +          16        4  20777188     0   2.08e7  2.08e7    0  2.08e7
#>  6 ADD1     +          4         8   2865706. 6694.  2.87e6  2.87e6 8804. 2.86e6
#>  7 ADGRA3   -          4         2  22505019     0   2.25e7  2.25e7    0  2.25e7
#>  8 ADGRB3   +          6         4  69397949     0   6.94e7  6.94e7    0  6.94e7
#>  9 ADGRG3   +          16        3  57716345     0   5.77e7  5.77e7    0  5.77e7
#> 10 AHR      +          7         2  17227687     0   1.72e7  1.72e7    0  1.72e7
#> # ℹ 491 more rows
#> # ℹ 27 more variables: max <dbl>, range <dbl>, skew <dbl>, kurtosis <dbl>,
#> #   n_IS_perGene <int>, min_bp_integration_locus <dbl>,
#> #   max_bp_integration_locus <dbl>, IS_span_bp <dbl>,
#> #   avg_bp_integration_locus <dbl>, median_bp_integration_locus <dbl>,
#> #   distinct_orientations <int>, average_TxLen <dbl>,
#> #   raw_gene_integration_frequency <dbl>, …
#> 
```
