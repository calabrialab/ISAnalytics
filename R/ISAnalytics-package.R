#' ISAnalytics: Analyze gene therapy vector insertion sites data
#' identified from genomics next generation sequencing reads for
#' clonal tracking studies
#'
#' @description \lifecycle{stable}
#' In gene therapy, stem cells are modified using viral
#' vectors to deliver the therapeutic transgene and replace functional
#' properties since the genetic modification is stable and inherited in
#' all cell progeny. The retrieval and mapping of the sequences flanking
#' the virus-host DNA junctions allows the identification of insertion
#' sites (IS), essential for monitoring the evolution of genetically
#' modified cells in vivo. A comprehensive toolkit for the analysis of
#' IS is required to foster clonal trackign studies and supporting the
#' assessment of safety and long term efficacy in vivo. This package
#' is aimed at (1) supporting automation of IS workflow, (2) performing
#' base and advance analysis for IS tracking (clonal abundance, clonal
#' expansions and statistics for insertional mutagenesis, etc.),
#' (3) providing basic biology insights of transduced stem cells in vivo.
#'
#' @section Useful resources:
#' * \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5702242/}{VISPA2:
#' A Scalable Pipeline for High-Throughput Identification
#' and Annotation of Vector Integration Sites}
#'
#' @section ISAnalytics function families:
#' * Import functions:
#'   * \code{\link{import_single_Vispa2Matrix}}
#'   * \code{\link{import_association_file}}
#'   * \code{\link{import_Vispa2_stats}}
#'   * \code{\link{import_parallel_Vispa2Matrices}}
#' * Aggregation functions:
#'   * \code{\link{aggregate_metadata}}
#'   * \code{\link{aggregate_values_by_key}}
#' * Collision removal functions:
#'   * \code{\link{remove_collisions}}
#'   * \code{\link{realign_after_collisions}}
#' * Removal of outliers from raw reads
#'   * \code{\link{outlier_filter}}
#'   * \code{\link{outliers_by_pool_fragments}}
#' * Recalibration functions:
#'   * \code{\link{compute_near_integrations}}
#' * Analysis functions:
#'   * \code{\link{compute_abundance}}
#'   * \code{\link{comparison_matrix}}
#'   * \code{\link{separate_quant_matrices}}
#'   * \code{\link{threshold_filter}}
#'   * \code{\link{top_integrations}}
#'   * \code{\link{sample_statistics}}
#'   * \code{\link{CIS_grubbs}}
#'   * \code{\link{cumulative_count_union}}
#'   * \code{\link{is_sharing}}
#'   * \code{\link{cumulative_is}}
#' * HSC population size estimate:
#'   * \code{\link{HSC_population_size_estimate}}
#' * Plotting functions:
#'   * \code{\link{CIS_volcano_plot}}
#'   * \code{\link{HSC_population_plot}}
#'   * \code{\link{sharing_heatmap}}
#'   * \code{\link{integration_alluvial_plot}}
#'   * \code{\link{top_abund_tableGrob}}
#'   * \code{\link{circos_genomic_density}}
#'   * \code{\link{sharing_venn}}
#' * Utility functions:
#'   * \code{\link{generate_blank_association_file}}
#'   * \code{\link{generate_Vispa2_launch_AF}}
#'   * \code{\link{unzip_file_system}}
#'   * \code{\link{as_sparse_matrix}}
#'
#' @section Vignettes:
#' * \code{vignette("how_to_import_functions", package = "ISAnalytics")}
#' * \code{vignette("collision_removal",
#' package = "ISAnalytics")}
#' * \code{vignette("aggregate_function_usage",
#' package = "ISAnalytics")}
#' * \code{vignette("sharing_analyses",
#' package = "ISAnalytics")}
#'
#' @docType package
#' @name ISAnalytics
#' @keywords internal
NULL

## usethis namespace: start
#' @import lifecycle
## usethis namespace: end
NULL
