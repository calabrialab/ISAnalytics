# ISAnalytics News

## Changes in version 0.99.10 (2020-09-14)

ISanalytics is officially on bioconductor!

#### NEW FEATURES

* Added analysis functions `comparison_matrix` and `separate_quant_matrices`
* Added utility function `as_sparse_matrix`
* Added package logo

#### SIGNIFICANT USER-VISIBLE CHANGES

* Changed algorithm for `compute_near_integrations`
* Added support for multi-quantification matrices to `remove_collisions`
* Added usage of lifecycle badges in documentation: users can now see if 
a feature is experimental/maturing/stable etc

#### MINOR FIXES

* Added fix for `import_single_Vispa2Matrix` to remove non significant 
0 values

## Changes in version 0.99.9 (2020-09-01)

#### NEW FEATURES

* Added functionality: aggregate functions
* Added vignette on aggregate functions
* Added recalibration functions
* Added first analysis function (compute_abundance)

#### SIGNIFICANT USER-VISIBLE CHANGES

* Dropped structure `ISADataFrame`: now the package only uses standard tibbles
* Modified package documentation

## Changes in version 0.99.8 (2020-08-12)
* Submitted to Bioconductor
