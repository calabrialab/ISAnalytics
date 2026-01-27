# Compute Fisher's exact test on gene frequencies.

**\[experimental\]** Provided 2 data frames with calculations for CIS,
via
[`CIS_grubbs()`](https://calabrialab.github.io/ISAnalytics/dev/reference/CIS_grubbs.md),
computes Fisher's exact test. Results can be plotted via
[`fisher_scatterplot()`](https://calabrialab.github.io/ISAnalytics/dev/reference/fisher_scatterplot.md).

## Usage

``` r
gene_frequency_fisher(
  cis_x,
  cis_y,
  min_is_per_gene = 3,
  gene_set_method = c("intersection", "union"),
  onco_db_file = "proto_oncogenes",
  tumor_suppressors_db_file = "tumor_suppressors",
  species = "human",
  known_onco = known_clinical_oncogenes(),
  suspicious_genes = clinical_relevant_suspicious_genes(),
  significance_threshold = 0.05,
  remove_unbalanced_0 = TRUE
)
```

## Arguments

- cis_x:

  A data frame obtained via
  [`CIS_grubbs()`](https://calabrialab.github.io/ISAnalytics/dev/reference/CIS_grubbs.md)

- cis_y:

  A data frame obtained via
  [`CIS_grubbs()`](https://calabrialab.github.io/ISAnalytics/dev/reference/CIS_grubbs.md)

- min_is_per_gene:

  Used for pre-filtering purposes. Genes with a number of distinct
  integration less than this number will be filtered out prior
  calculations. Single numeric or integer.

- gene_set_method:

  One between "intersection" and "union". When merging the 2 data
  frames, `intersection` will perform an inner join operation, while
  `union` will perform a full join operation.

- onco_db_file:

  Uniprot file for proto-oncogenes (see details). If different from
  default, should be supplied as a path to a file.

- tumor_suppressors_db_file:

  Uniprot file for tumor-suppressor genes. If different from default,
  should be supplied as a path to a file.

- species:

  One between `"human"`, `"mouse"` and `"all"`

- known_onco:

  Data frame with known oncogenes. See details.

- suspicious_genes:

  Data frame with clinical relevant suspicious genes. See details.

- significance_threshold:

  Significance threshold for the Fisher's test p-value

- remove_unbalanced_0:

  Remove from the final output those pairs in which there are no IS for
  one group or the other and the number of IS of the non-missing group
  are less than the mean number of IS for that group

## Value

A data frame

## Details

### Oncogene and tumor suppressor genes files

These files are included in the package for user convenience and are
simply UniProt files with gene annotations for human and mouse. For more
details on how this files were generated use the help
[`?tumor_suppressors`](https://calabrialab.github.io/ISAnalytics/dev/reference/proto_oncogenes.md),
[`?proto_oncogenes`](https://calabrialab.github.io/ISAnalytics/dev/reference/proto_oncogenes.md)

### Known oncogenes

The default values are included in this package and it can be accessed
by doing:

    known_clinical_oncogenes()

If the user wants to change this parameter the input data frame must
preserve the column structure. The same goes for the `suspicious_genes`
parameter (DOIReference column is optional):

    clinical_relevant_suspicious_genes()

## Required tags

The function will explicitly check for the presence of these tags:

- gene_symbol

## See also

Other Analysis functions:
[`CIS_grubbs()`](https://calabrialab.github.io/ISAnalytics/dev/reference/CIS_grubbs.md),
[`HSC_population_size_estimate()`](https://calabrialab.github.io/ISAnalytics/dev/reference/HSC_population_size_estimate.md),
[`compute_abundance()`](https://calabrialab.github.io/ISAnalytics/dev/reference/compute_abundance.md),
[`cumulative_is()`](https://calabrialab.github.io/ISAnalytics/dev/reference/cumulative_is.md),
[`is_sharing()`](https://calabrialab.github.io/ISAnalytics/dev/reference/is_sharing.md),
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
cis <- CIS_grubbs(aggreg, by = "SubjectID")
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
fisher <- gene_frequency_fisher(cis$cis$PT001, cis$cis$PT002,
    min_is_per_gene = 2
)
fisher
#> # A tibble: 1 × 28
#>   GeneName n_IS_perGene_1 average_TxLen_1 raw_gene_integration_frequency_1
#>   <chr>             <int>           <dbl>                            <dbl>
#> 1 KMT5B                 2          54306.                        0.0000368
#> # ℹ 24 more variables: IS_per_kbGeneLen_1 <dbl>, Sum_IS_per_kbGeneLen_1 <dbl>,
#> #   IS_per_kbGeneLen_perMDepth_TPM_1 <dbl>, n_IS_perGene_2 <int>,
#> #   average_TxLen_2 <dbl>, raw_gene_integration_frequency_2 <dbl>,
#> #   IS_per_kbGeneLen_2 <dbl>, Sum_IS_per_kbGeneLen_2 <dbl>,
#> #   IS_per_kbGeneLen_perMDepth_TPM_2 <dbl>, Onco1_TS2 <dbl>,
#> #   KnownClonalExpansion <lgl>, ClinicalRelevance <lgl>, DOIReference <chr>,
#> #   KnownGeneClass <chr>, CriticalForInsMut <lgl>, tot_n_IS_perGene_1 <int>, …
```
