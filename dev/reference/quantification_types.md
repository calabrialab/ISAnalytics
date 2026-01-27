# Possible choices for the `quantification_type` parameter.

These are all the possible values for the `quantification_type`
parameter in `import_parallel_vispa2Matrices_interactive` and
`import_parallel_vispa2Matrices_auto`.

## Usage

``` r
quantification_types()
```

## Value

A vector of characters for quantification types

## Details

The possible values are:

- fragmentEstimate

- seqCount

- barcodeCount

- cellCount

- ShsCount

## See also

[`import_parallel_Vispa2Matrices_interactive`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_parallel_Vispa2Matrices_interactive.md),
[`import_parallel_Vispa2Matrices_auto`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_parallel_Vispa2Matrices_auto.md)

Other Import functions helpers:
[`annotation_issues()`](https://calabrialab.github.io/ISAnalytics/dev/reference/annotation_issues.md),
[`date_formats()`](https://calabrialab.github.io/ISAnalytics/dev/reference/date_formats.md),
[`default_af_transform()`](https://calabrialab.github.io/ISAnalytics/dev/reference/default_af_transform.md),
[`default_iss_file_prefixes()`](https://calabrialab.github.io/ISAnalytics/dev/reference/default_iss_file_prefixes.md),
[`matching_options()`](https://calabrialab.github.io/ISAnalytics/dev/reference/matching_options.md)

## Examples

``` r
quant_types <- quantification_types()
```
