\title{ISAnalytics News}

# ISAnalytics 1.0.7 (2021-02-06)

## UPGRADES

* It is now possible to save html reports to file from 
import_parallel_Vispa2Matrices_auto and 
import_parallel_Vispa2Matrices_interactive, remove_collisions and 
compute_near_integrations

## FIXES

* Fixed sample_statistics: now functions that have data frame output do not
produce nested tables. Flat tables are ready to be saved to file or can be
nested.
* Simplified association file check logic in remove_collisions: now 
function blocks only if the af doesn't contain the needed columns

# ISAnalytics 1.0.6 (2021-02-03)

## UPGRADES

* Upgraded import_association_file function: now file alignment is not
mandatory anymore and it is possible to save the html report to file
* Updated vignettes and documentation

# ISAnalytics 1.0.5 (2020-11-16)

## UPGRADES

* Greatly improved reports for collision removal functionality
* General improvements for all widget reports

# ISAnalytics 1.0.4 (2020-11-10)

## FIXES

* Further fixes for printing reports when widgets not available
* Added progress bar to collision processing in `remove_collisions`
* Updated vignettes

## NEW

* Added vignette "Using ISAnalytics without RStudio support"

# ISAnalytics 1.0.3 (2020-11-05)

## FIXES

* Fixed missing restarts for non-blocking widgets

# ISAnalytics 1.0.2 (2020-11-04)

## FIXES

* Functions that make use of widgets do not interrupt execution anymore if 
errors are thrown while producing or printing the widgets
* Optimized widget printing for importing functions
* If widgets can't be printed and verbose option is active, reports are now 
displayed on console instead (needed for usage in environments that do not 
have access to a browser)
* Other minor fixes (typos)

## IMPORTANT NOTES

* Vignette building might fail due to the fact that package "knitcitations" 
is temporarily unavailable through CRAN

# ISAnalytics 1.0.1 (2020-10-29)

* ISAnalytics is finally in release on bioconductor!
* Bug fixes: fixed a few bugs in importing and recalibration functions

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
