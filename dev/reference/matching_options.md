# Possible choices for the `matching_opt` parameter.

These are all the possible values for the `matching_opt` parameter in
`import_parallel_vispa2Matrices_auto`.

## Usage

``` r
matching_options()
```

## Value

A vector of characters for matching_opt

## Details

The values "ANY", "ALL" and "OPTIONAL", represent how the patterns
should be matched, more specifically

- ANY = look only for files that match AT LEAST one of the patterns
  specified

- ALL = look only for files that match ALL of the patterns specified

- OPTIONAL = look preferentially for files that match, in order, all
  patterns or any pattern and if no match is found return what is found
  (keep in mind that duplicates are discarded in automatic mode)

## See also

[`import_parallel_Vispa2Matrices_auto`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_parallel_Vispa2Matrices_auto.md)

Other Import functions helpers:
[`annotation_issues()`](https://calabrialab.github.io/ISAnalytics/dev/reference/annotation_issues.md),
[`date_formats()`](https://calabrialab.github.io/ISAnalytics/dev/reference/date_formats.md),
[`default_af_transform()`](https://calabrialab.github.io/ISAnalytics/dev/reference/default_af_transform.md),
[`default_iss_file_prefixes()`](https://calabrialab.github.io/ISAnalytics/dev/reference/default_iss_file_prefixes.md),
[`quantification_types()`](https://calabrialab.github.io/ISAnalytics/dev/reference/quantification_types.md)

## Examples

``` r
opts <- matching_options()
```
