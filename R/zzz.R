
.onAttach <- function(libname, pkgname) {
    options(ISAnalytics.connection = stdin())
    options(ISAnalytics.verbose = TRUE)
    options(ISAnalytics.reports = TRUE)
    options(ISAnalytics.mandatory_is_vars = "default")
    options(ISAnalytics.genomic_annotation_vars = "default")
    options(ISAnalytics.af_specs = "default")
    options(ISAnalytics.iss_stats_specs = "default")
    options(ISAnalytics.matrix_file_suffix = "default")
    options(ISAnalytics.parallel_processing = TRUE)
}
