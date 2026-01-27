# Retrieve description of a tag by name.

Given one or multiple tags, prints the associated description and
functions where the tag is explicitly used.

## Usage

``` r
inspect_tags(tags)
```

## Arguments

- tags:

  A character vector of tag names

## Value

`NULL`

## See also

Other dynamic vars:
[`mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/mandatory_IS_vars.md),
[`pcr_id_column()`](https://calabrialab.github.io/ISAnalytics/dev/reference/pcr_id_column.md),
[`reset_mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/reset_mandatory_IS_vars.md),
[`set_mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/set_mandatory_IS_vars.md),
[`set_matrix_file_suffixes()`](https://calabrialab.github.io/ISAnalytics/dev/reference/set_matrix_file_suffixes.md)

## Examples

``` r
inspect_tags(c("chromosome", "project_id", "x"))
#> * TAG: chromosome
#> ℹ Description: Number of the chromosome
#> ℹ Functions that use it: top_targeted_genes, CIS_grubbs, compute_near_integrations
#> * TAG: project_id
#> ℹ Description: Unique identifier of a project
#> ℹ Functions that use it: generate_default_folder_structure, import_Vispa2_stats, remove_collisions, generate_Vispa2_launch_AF, import_association_file, import_parallel_Vispa2Matrices
#> * TAG: x
#> ✖ Tag not found in available tags
```
