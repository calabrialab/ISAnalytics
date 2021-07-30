#### ---- Internals for reports production ----####

# Name of the folder where templates are stored
.templates_folder <- function() {
    "rmd"
}

# All associated information to each type of report
.templates_info <- function() {
    list(
        collisions = list(
            template_name = "collision-report.Rmd",
            required_pkgs = c("flexdashboard", "plotly", "DT", "ggvenn"),
            def_filename = "collision_removal_report.html"
        ),
        vispa2_stats = list(
            template_name = "iss-import-report.Rmd",
            required_pkgs = c("flexdashboard", "DT"),
            def_filename = "vispa2_stats_import_report.html"
        ),
        asso_file = list(
            template_name = "af-report.Rmd",
            required_pkgs = c("flexdashboard", "DT"),
            def_filename = "association_file_import_report.html"
        ),
        matrix_imp = list(
            template_name = "matrix-import-report.Rmd",
            required_pkgs = c("flexdashboard", "DT"),
            def_filename = "matrices_import_report.html"
        ),
        outlier_flag = list(
            template_name = "outlier-report.Rmd",
            required_pkgs = c("flexdashboard", "DT"),
            def_filename = "raw_reads_outliers_report.html"
        )
    )
}

# Gets the default file name for the given report type
.get_default_rep_filename <- function(type) {
    return(.templates_info()[[type]]$def_filename)
}

# Retrieves the template file path
.get_template <- function(type) {
    filename <- .templates_info()[[type]]$template_name
    system.file(.templates_folder(), filename, package = "ISAnalytics")
}

# Retrieves all the required packages for the given report type
.get_sugg_packages <- function(type) {
    .templates_info()[[type]]$required_pkgs
}

# Renders the report with the appropriate parameters
.produce_report <- function(report_type, params, path) {
    if (!getOption("ISAnalytics.reports") == TRUE || is.null(path)) {
        return(NULL)
    }
    if (getOption("ISAnalytics.verbose") == TRUE) {
        rlang::inform("Producing report...")
    }
    if (!requireNamespace("rmarkdown")) {
        rlang::inform(.missing_pkg_error("rmarkdown"))
        return(NULL)
    }

    pkgs_present <- purrr::map_lgl(
        .get_sugg_packages(report_type),
        ~ requireNamespace(.x, quietly = TRUE)
    )
    if (any(pkgs_present == FALSE)) {
        missing_pkgs <- .get_sugg_packages(report_type)[!pkgs_present]
        rlang::inform(.missing_pkg_error(missing_pkgs[1]))
        return(NULL)
    }
    if (!is.character(path)) {
        not_str_path <- paste(
            "Provided report path is",
            "not a string, using default"
        )
        rlang::inform(not_str_path)
        path <- default_report_path()
    }
    template <- .get_template(report_type)
    path <- .clean_file_path(path, report_type)
    withRestarts(
        {
            rmarkdown::render(
                input = template,
                params = params,
                output_file = path,
                envir = new.env(),
                quiet = TRUE
            )
            rlang::inform(.report_saved_info(path))
        },
        report_fail = function(e) {
            rlang::inform(.report_fail_err(conditionMessage(e)))
        }
    )
}

# Gets a cleaned file path to the report file
.clean_file_path <- function(path, type) {
    if (fs::is_dir(path)) {
        if (!fs::dir_exists(path)) {
            fs::dir_create(path)
        }
        gen_filename <- .generate_report_filename(type)
        path <- fs::path(path, gen_filename)
    } else {
        if (fs::path_ext(path) == "") {
            path <- fs::path_ext_set(path, "html")
        }
    }
    return(path)
}

# Generates a default report filename if one is not provided in input
.generate_report_filename <- function(type) {
    def <- .get_default_rep_filename(type)
    date <- lubridate::today()
    return(paste0(date, "_", def))
}

.report_fail_err <- function(err) {
    c("Failure",
        x = "Report production failed, skipping",
        i = paste("Error: ", err)
    )
}

.report_saved_info <- function(file) {
    c("Report correctly saved",
        i = paste("Report saved to:", file)
    )
}

#' Default folder for saving ISAnalytics reports. Supplied as default
#' argument for several functions.
#'
#' @return A path
#' @importFrom fs path_home
#' @export
#'
#' @examples
#' default_report_path()
default_report_path <- function() {
    fs::path_home("ISAnalytics_reports")
}
