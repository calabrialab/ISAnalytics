\title{ISAnalytics News}

# ISAnalytics 0.99.15 (2020-10-22)

* Minor fix in import_association_file file function: added multiple strings
to be translated as NA

# ISAnalytics 0.99.14 (2020-10-21)

* Minor fixes in tests

# ISAnalytics 0.99.13 (2020-10-19)

## NEW FEATURES

* Added analysis functions `CIS_grubbs` and `cumulative_count_union`
* Added plotting functions `CIS_volcano_plot`

# ISAnalytics 0.99.12 (2020-10-04)

## NEW FEATURES

* Added analysis function `sample_statistics`

## SIGNIFICANT USER-VISIBLE CHANGES

* `aggregate_values_by_key` has a simplified interface and supports
multi-quantification matrices

## MINOR CHANGES

* Updated vignettes
* `import_parallel_Vispa2Matrices_interactive` and
`import_parallel_Vispa2Matrices_auto` now have an option to return 
a multi-quantification matrix directly after import instead of a list

# ISAnalytics 0.99.11 (2020-09-21)

## NEW FEATURES

* Added analysis functions `threshold_filter`, `top_integrations`
* Added support for multi-quantification matrices in `compute_abundance`

## MINOR FIXES

* Fixed bug in `comparison_matrix` that ignored custom column names
* Fixed issues in some documentation pages

# ISAnalytics 0.99.10 (2020-09-14)

ISanalytics is officially on bioconductor!

## NEW FEATURES

* Added analysis functions `comparison_matrix` and `separate_quant_matrices`
* Added utility function `as_sparse_matrix`
* Added package logo

## SIGNIFICANT USER-VISIBLE CHANGES

* Changed algorithm for `compute_near_integrations`
* Added support for multi-quantification matrices to `remove_collisions`
* Added usage of lifecycle badges in documentation: users can now see if 
a feature is experimental/maturing/stable etc

## MINOR FIXES

* Added fix for `import_single_Vispa2Matrix` to remove non significant 
0 values

# ISAnalytics 0.99.9 (2020-09-01)

## NEW FEATURES

* Added functionality: aggregate functions
* Added vignette on aggregate functions
* Added recalibration functions
* Added first analysis function (compute_abundance)

## SIGNIFICANT USER-VISIBLE CHANGES

* Dropped structure `ISADataFrame`: now the package only uses standard tibbles
* Modified package documentation

# ISAnalytics 0.99.8 (2020-08-12)

* Submitted to Bioconductor
