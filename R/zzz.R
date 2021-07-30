
.onAttach <- function(libname, pkgname) {
    options(ISAnalytics.connection = stdin())
    options(ISAnalytics.verbose = TRUE)
    options(ISAnalytics.reports = TRUE)
}
