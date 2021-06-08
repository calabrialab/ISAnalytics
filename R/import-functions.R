#------------------------------------------------------------------------------#
# Importing functions
#------------------------------------------------------------------------------#

#' Import a single integration matrix from file
#'
#' @description \lifecycle{stable}
#' This function allows to read and import an integration matrix
#' produced as the output of Vispa2 pipeline and converts it to a tidy
#' format.
#'
#' @param path The path to the file on disk
#' @param to_exclude Either NULL or a character vector of column names that
#' should be ignored when importing
#' @param separator The column delimiter used
#'
#' @return A data.table object in tidy format
#' @family Import functions
#' @importFrom rlang abort inform
#' @importFrom fs path_ext
#' @importFrom readr read_delim cols
#' @importFrom tidyr separate
#' @importFrom magrittr `%>%`
#' @importFrom dplyr mutate
#' @importFrom stringr str_replace
#' @importFrom BiocParallel SnowParam MulticoreParam bplapply bpstop
#' @importFrom data.table melt.data.table rbindlist
#' @details The import series of functions is designed to work in combination
#' with the use of Vispa2 pipeline, please refer to this article for more
#' details: \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5702242/}{VISPA2:
#' A Scalable Pipeline for High-Throughput Identification and Annotation of
#' Vector Integration Sites}.
#' For more details on how to properly use these functions, refer to
#' \code{vignette("How to use import functions", package = "ISAnalytics")}
#' @export
#'
#' @examples
#' path_to_file <- system.file("extdata", "ex_annotated_ISMatrix.tsv.xz",
#'     package = "ISAnalytics"
#' )
#' isa_dataframe <- import_single_Vispa2Matrix(path_to_file)
import_single_Vispa2Matrix <- function(path,
    to_exclude = NULL,
    separator = "\t") {
    stopifnot(!missing(path) & is.character(path))
    if (!file.exists(path)) {
        rlang::abort(paste("File not found at", path))
    }
    if (!fs::is_file(path)) {
        rlang::abort(paste("Path exists but is not a file"))
    }
    mode <- "fread"
    ## Is the file compressed?
    is_compressed <- fs::path_ext(path) %in% .compressed_formats()
    if (is_compressed) {
        ## The compression type is supported by data.table::fread?
        compression_type <- fs::path_ext(path)
        if (!compression_type %in% .supported_fread_compression_formats()) {
            ### If not, switch to classic for reading
            mode <- "classic"
            if (getOption("ISAnalytics.verbose") == TRUE) {
                rlang::inform(.unsupported_comp_format_inf(),
                    class = "unsup_comp_format"
                )
            }
        }
    }
    ### Peak headers
    peek_headers <- readr::read_delim(path,
        delim = separator, n_max = 0,
        col_types = readr::cols()
    )
    ## - Detect type
    df_type <- .auto_detect_type(peek_headers)
    if (df_type == "MALFORMED") {
        rlang::abort(.malformed_ISmatrix_error(),
            class = "malformed_ism"
        )
    }
    is_annotated <- .is_annotated(peek_headers)
    ## - Start reading
    if (getOption("ISAnalytics.verbose") == TRUE) {
        rlang::inform(c("Reading file...", i = paste0("Mode: ", mode)))
    }
    df <- if (mode == "fread") {
        .read_with_fread(
            path = path, to_drop = to_exclude,
            df_type = df_type, annotated = is_annotated,
            sep = separator
        )
    } else {
        .read_with_readr(
            path = path, to_drop = to_exclude,
            df_type = df_type, annotated = is_annotated,
            sep = separator
        )
    }
    ## - Report summary
    if (getOption("ISAnalytics.verbose") == TRUE) {
        rlang::inform(.summary_ism_import_msg(
            df_type,
            .is_annotated(df),
            dim(df),
            mode
        ),
        class = "ism_import_summary"
        )
    }
    if (df_type == "OLD") {
        df <- df %>%
            tidyr::separate(
                col = .data$IS_genomicID,
                into = mandatory_IS_vars(),
                sep = "_", remove = TRUE,
                convert = TRUE
            ) %>%
            dplyr::mutate(chr = stringr::str_replace(
                .data$chr, "chr", ""
            ))
    }
    ## - Split in chunks
    if (getOption("ISAnalytics.verbose") == TRUE) {
        rlang::inform("Reshaping...")
    }
    chunks <- split(df,
        by = c("chr"),
        verbose = FALSE
    )
    ## - Melt in parallel
    p <- if (.Platform$OS.type == "windows") {
        BiocParallel::SnowParam(
            tasks = length(chunks),
            progressbar = getOption("ISAnalytics.verbose"),
            exportglobals = FALSE,
            stop.on.error = TRUE
        )
    } else {
        BiocParallel::MulticoreParam(
            tasks = length(chunks),
            progressbar = getOption("ISAnalytics.verbose"),
            exportglobals = FALSE,
            stop.on.error = TRUE
        )
    }
    mt <- function(data, annot) {
        id_vars <- if (annot) {
            c(
                mandatory_IS_vars(),
                annotation_IS_vars()
            )
        } else {
            mandatory_IS_vars()
        }
        data.table::melt.data.table(data,
            id.vars = id_vars,
            variable.name = "CompleteAmplificationID",
            value.name = "Value",
            na.rm = TRUE,
            verbose = FALSE
        )
    }
    tidy_chunks <- BiocParallel::bplapply(
        X = chunks,
        FUN = mt,
        annot = is_annotated,
        BPPARAM = p
    )
    BiocParallel::bpstop(p)
    tidy <- data.table::rbindlist(tidy_chunks)
    tidy <- tidy["Value" > 0]
    if (getOption("ISAnalytics.verbose") == TRUE) {
        rlang::inform("Done!")
    }
    return(tidy)
}


#' Import the association file from disk
#'
#' @description \lifecycle{maturing}
#' Imports the association file and immediately performs a check on
#' the file system starting from the root to assess the alignment between the
#' two.
#' @param path The path on disk to the association file.
#' @param root The path on disk of the root folder of Vispa2 output or NULL.
#' See details.
#' @param tp_padding Timepoint padding, indicates the number of digits of the
#' "Timepoint" column once imported. Fills the content with 0s up to the length
#' specified (ex: 1 becomes 0001 with a tp_padding of 4)
#' @param dates_format A single string indicating how dates should be parsed.
#' Must be a value in: \code{date_formats()}
#' @param separator The column separator used in the file
#' @param filter_for A named list where names represent column names that
#' must be filtered. For example: `list(ProjectID = c("PROJECT1", "PROJECT2))`
#' will filter the association file so that it contains only those rows
#' for which the value of the column "ProjectID" is one of the specified
#' values. If multiple columns are present in the list all filtering
#' conditions are applied as a logical AND.
#' @param import_iss Import Vispa2 stats and merge them with the
#' association file?
#' @param export_widget_path A path on disk to save produced widgets or NULL
#' if the user doesn't wish to save the html file
#' @param ... Additional arguments to pass to \code{\link{import_Vispa2_stats}}
#' @family Import functions
#' @return A tibble with the contents of the association file plus columns
#' containing the path in the file system for every project and pool if found.
#' @details The import series of functions is designed to work in combination
#' with the use of Vispa2 pipeline, please refer to this article for more
#' details: \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5702242/}{VISPA2:
#' A Scalable Pipeline for High-Throughput Identification
#' and Annotation of Vector Integration Sites}.\cr
#' The pipeline automatically produces an hierarchical structure in the file
#' system which follows this schema:
#' * /root_folder
#'   * Optional intermediate folders
#'     * ProjectID\cr
#'       |_bam\cr
#'       |_bcmuxall\cr
#'       |_bed\cr
#'       |_iss\cr
#'       |_quality\cr
#'       |_report\cr
#'       |_quantification\cr
#'        *|___concatenatePoolIDSeqRun\cr
#'
#' For each ProjectID there may be several nested PoolIDs. The alignment
#' function only looks for PoolIDs in the quantification folder, since it's
#' the location of the matrices to import.
#' For more details on how to properly use these functions, refer to the
#' vignette - vignette("how_to_import_functions").\cr
#' If 'NULL' the file system alignment step is skipped.
#' @export
#'
#' @importFrom purrr map_lgl set_names is_empty
#' @importFrom rlang inform abort dots_list exec
#' @importFrom htmltools tagList browsable
#' @importFrom tibble add_column
#' @import dplyr
#' @importFrom stringr str_pad
#' @importFrom magrittr `%>%`
#' @seealso \code{\link{date_formats}}
#' @examples
#' op <- options("ISAnalytics.widgets" = FALSE)
#' path <- system.file("extdata", "ex_association_file.tsv",
#'     package = "ISAnalytics"
#' )
#' root_pth <- system.file("extdata", "fs.zip", package = "ISAnalytics")
#' root <- unzip_file_system(root_pth, "fs")
#' association_file <- import_association_file(path, root, dates_format = "dmy")
#' options(op)
import_association_file <- function(path,
    root = NULL, tp_padding = 4, dates_format = "ymd",
    separator = "\t",
    filter_for = NULL,
    import_iss = FALSE,
    export_widget_path = NULL,
    convert_tp = TRUE,
    ...) {
    # Check parameters
    stopifnot(is.character(path) & length(path) == 1)
    stopifnot((is.character(root) & length(root) == 1) || (is.null(root)))
    stopifnot(file.exists(path))
    if (!is.null(root) && root != "") {
        stopifnot(file.exists(root))
    }
    stopifnot((is.numeric(tp_padding) |
        is.integer(tp_padding)) & length(tp_padding) == 1)
    stopifnot(length(dates_format) == 1 & dates_format %in% date_formats())
    stopifnot(is.character(separator) && length(separator) == 1)
    stopifnot(is.logical(import_iss) && length(import_iss) == 1)
    if (import_iss & is.null(root)) {
        rlang::abort(paste(
            "Can't import Vispa2 stats files without",
            "file system alignment. Provide the appropriate",
            "root."
        ))
    }
    # Check filter
    stopifnot(is.null(filter_for) ||
        (is.list(filter_for) && !is.null(names(filter_for))))
    # Read file and check the correctness
    af_checks <- .manage_association_file(
        path, root, tp_padding, dates_format,
        separator, filter_for
    )
    as_file <- af_checks$af
    parsing_problems <- af_checks$parsing_probs
    date_problems <- af_checks$date_probs
    checks <- af_checks$check
    col_probs <- NULL
    if (!.check_af_correctness(as_file)) {
        col_probs[["missing"]] <- association_file_columns()[
            !association_file_columns() %in% colnames(as_file)
        ]
    }
    non_standard <- colnames(as_file)[
        !colnames(as_file) %in% c(association_file_columns(), "Path")
    ]
    if (!purrr::is_empty(non_standard)) {
        col_probs[["non_standard"]] <- non_standard
    }
    missing_dates <- purrr::map_lgl(date_columns_coll(), function(date_col) {
        any(is.na(as_file[[date_col]]))
    }) %>% purrr::set_names(date_columns_coll())
    missing_dates <- names(missing_dates)[missing_dates == TRUE]
    something_to_report <- any(!is.null(c(
        parsing_problems,
        date_problems,
        checks,
        col_probs,
        missing_dates
    )))
    ## Fix timepoints
    if (convert_tp) {
        if (!"TimepointMonths" %in% colnames(as_file)) {
            as_file <- as_file %>%
                tibble::add_column(TimepointMonths = NA_real_)
        }
        if (!"TimepointYears" %in% colnames(as_file)) {
            as_file <- as_file %>%
                tibble::add_column(TimepointYears = NA_real_)
        }
        as_file <- as_file %>%
            dplyr::mutate(
                TimepointMonths = dplyr::if_else(
                        condition = as.numeric(.data$TimePoint) == 0,
                        true = 0,
                        false = dplyr::if_else(
                            condition = as.numeric(.data$TimePoint) > 0 &
                                as.numeric(.data$TimePoint) < 30,
                            true = ceiling(as.numeric(.data$TimePoint) / 30),
                            false = round(as.numeric(.data$TimePoint) / 30)
                        )
                    ),
                TimepointYears = dplyr::if_else(
                        condition = as.numeric(.data$TimePoint) == 0,
                        true = 0,
                        false = ceiling(as.numeric(.data$TimePoint) / 360)
                    )
                ) %>%
            dplyr::mutate(
                TimepointMonths = stringr::str_pad(
                    as.character(.data$TimepointMonths),
                    pad = "0", side = "left", width = 2
                ),
                TimepointYears = stringr::str_pad(
                    as.character(.data$TimepointYears),
                    pad = "0", side = "left", width = 2
                ),
            )
    }

    widget <- if (something_to_report) {
        summary_report <- .summary_af_import_msg(
            pars_prob = parsing_problems, dates_prob = date_problems,
            cols_prob = col_probs, crit_na = missing_dates,
            checks = ifelse(is.null(checks),
                yes = "skipped",
                no = ifelse(any(!checks$Found),
                    "problems detected",
                    "no problems detected"
                )
            )
        )
        if (getOption("ISAnalytics.widgets") == TRUE) {
            .produce_widget(".checker_widget",
                parsing_probs = parsing_problems,
                date_probs = date_problems,
                checker_df = checks,
                col_probs = col_probs,
                critical_nas = missing_dates
            )
        } else if (getOption("ISAnalytics.verbose") == TRUE) {
            rlang::inform(summary_report,
                class = "summary_report"
            )
            NULL
        } else {
            NULL
        }
    }
    if (import_iss) {
        dots <- rlang::dots_list(.named = TRUE)
        dots <- dots[!names(dots) %in% c(
            "association_file",
            "export_widget_path",
            "join_with_af"
        )]
        stats <- withCallingHandlers(
            {
                withRestarts(
                    {
                        rlang::exec(import_Vispa2_stats,
                            association_file = as_file,
                            export_widget_path = "INTERNAL",
                            join_with_af = TRUE,
                            !!!dots
                        )
                    },
                    fail_stats = function() {
                        rlang::inform("Issues in importing stats files, skipping")
                    }
                )
            },
            error = function(err) {
                rlang::inform(err$message)
                invokeRestart("fail_stats")
            }
        )
        if (!is.null(stats)) {
            as_file <- stats$stats
            if (getOption("ISAnalytics.widgets") == TRUE) {
                widget <- htmltools::browsable(
                    htmltools::tagList(
                        widget, stats$report_w
                    )
                )
            }
        }
    }
    if (!is.null(widget)) {
        .print_widget(widget, else_verbose = {
            rlang::inform(summary_report,
                class = "summary_report"
            )
        })
    } else if (getOption("ISAnalytics.verbose") == TRUE) {
        rlang::inform(summary_report,
            class = "summary_report"
        )
    }
    if (!is.null(export_widget_path)) {
        .export_widget_file(widget,
            path = export_widget_path,
            def_file_name = "af_import_report.html"
        )
    }
    as_file
}

#' Import Vispa2 stats given the aligned association file.
#'
#' \lifecycle{experimental}
#' Imports all the Vispa2 stats files for each pool provided the association
#' file has been aligned with the file system
#' (see \code{\link{import_association_file}}).
#'
#' @param association_file The file system aligned association file
#' (contains columns with absolute paths to the 'iss' folder)
#' @param file_prefixes A character vector with known file prefixes
#'  to match on file names. NOTE: the elements represent regular expressions.
#'  For defaults see \link{default_iss_file_prefixes}.
#' @param join_with_af Logical, if TRUE the imported stats files will be
#' merged with the association file, if false a single data frame holding
#' only the stats will be returned.
#' @param pool_col A single string. What is the name of the pool column
#' used in the Vispa2 run? This will be used as a key to perform a join
#' operation with the stats files `POOL` column.
#' @param export_widget_path Either NULL or the path on disk where the
#' widget report should be saved.
#'
#' @family Import functions
#' @importFrom rlang inform abort
#' @importFrom stats setNames
#' @import dplyr
#' @importFrom htmltools browsable tagList
#'
#' @return A data frame
#' @export
#'
#' @examples
#' op <- options("ISAnalytics.widgets" = FALSE)
#' path <- system.file("extdata", "ex_association_file.tsv",
#'     package = "ISAnalytics"
#' )
#' root_pth <- system.file("extdata", "fs.zip", package = "ISAnalytics")
#' root <- unzip_file_system(root_pth, "fs")
#' association_file <- import_association_file(path, root, dates_format = "dmy")
#' af_with_stats <- import_Vispa2_stats(association_file)
#' options(op)
import_Vispa2_stats <- function(association_file,
    file_prefixes = default_iss_file_prefixes(),
    join_with_af = TRUE,
    pool_col = "concatenatePoolIDSeqRun",
    export_widget_path = NULL) {
    ## Check param
    if (!is.data.frame(association_file)) {
        rlang::abort(.af_not_imported_err())
    }
    path_cols <- .path_cols_names()
    if (!path_cols$iss %in% colnames(association_file)) {
        rlang::abort(.af_not_aligned_err())
    }
    min_cols <- c("ProjectID", "concatenatePoolIDSeqRun", path_cols$iss)
    if (!all(min_cols %in% colnames(association_file))) {
        rlang::abort(
            .missing_needed_cols(
                min_cols[!min_cols %in% colnames(association_file)]
            )
        )
    }
    stopifnot(is.character(file_prefixes))
    stopifnot(is.logical(join_with_af))
    join_with_af <- join_with_af[1]
    if (join_with_af) {
        stopifnot(is.character(pool_col) & length(pool_col) == 1)
        stopifnot(pool_col %in% colnames(association_file))
    }
    ## export path has a special placeholder "INTERNAL" for calling
    ## the function inside import_association_file
    stopifnot(is.null(export_widget_path) || is.character(export_widget_path))
    ## Import
    stats <- .import_stats_iss(
        association_file = association_file,
        prefixes = file_prefixes
    )
    report <- stats$report
    stats <- stats$stats
    ## Produce widget report if requested
    widget_stats_import <- if (getOption("ISAnalytics.widgets") == TRUE) {
        .produce_widget(".iss_import_widget", report = report)
    } else {
        NULL
    }
    ## - IF NO STATS IMPORTED (STATS ARE NULL)
    if (is.null(stats)) {
        if (getOption("ISAnalytics.verbose") == TRUE) {
            rlang::inform(paste("No stats files imported"))
        }
        if (!is.null(widget_stats_import)) {
            if (is.null(export_widget_path)) {
                .print_widget(widget_stats_import)
            } else if (export_widget_path != "INTERNAL") {
                .print_widget(widget_stats_import)
            }
        }
        if (!is.null(export_widget_path) && export_widget_path != "INTERNAL") {
            .export_widget_file(
                widget_stats_import,
                export_widget_path,
                "vispa2_stats_import_report.html"
            )
        }
        if (!is.null(export_widget_path) && export_widget_path == "INTERNAL") {
            return(list(stats = NULL, report_w = widget_stats_import))
        } else {
            return(NULL)
        }
    }
    ## - IF STATS NOT NULL
    ## Merge if requested
    if (join_with_af) {
        association_file <- association_file %>%
            dplyr::left_join(stats, by = c(
                stats::setNames("POOL", pool_col),
                "TagSequence" = "TAG"
            ))
        ## Detect potential problems
        missing_stats <- association_file %>%
            dplyr::filter(is.na(.data$RUN_NAME)) %>%
            dplyr::select(
                .data$ProjectID,
                .data$concatenatePoolIDSeqRun,
                dplyr::all_of(pool_col),
                .data$CompleteAmplificationID
            ) %>%
            dplyr::distinct()
        widget_stats_import <- if (getOption("ISAnalytics.widgets") == TRUE) {
            miss_widget <- .produce_widget(".missing_iss_widget",
                missing_iss = missing_stats
            )
            htmltools::browsable(
                htmltools::tagList(widget_stats_import, miss_widget)
            )
        } else {
            NULL
        }
        if (!is.null(widget_stats_import)) {
            if (is.null(export_widget_path)) {
                .print_widget(widget_stats_import)
            } else if (export_widget_path != "INTERNAL") {
                .print_widget(widget_stats_import)
            }
        }
        if (!is.null(export_widget_path) && export_widget_path != "INTERNAL") {
            .export_widget_file(
                widget_stats_import,
                export_widget_path,
                "vispa2_stats_import_report.html"
            )
        }
        if (!is.null(export_widget_path) && export_widget_path == "INTERNAL") {
            return(list(
                stats = association_file,
                report_w = widget_stats_import
            ))
        } else {
            return(association_file)
        }
    }
    if (!is.null(widget_stats_import)) {
        if (is.null(export_widget_path)) {
            .print_widget(widget_stats_import)
        } else if (export_widget_path != "INTERNAL") {
            .print_widget(widget_stats_import)
        }
    }
    if (!is.null(export_widget_path) && export_widget_path != "INTERNAL") {
        .export_widget_file(
            widget_stats_import,
            export_widget_path,
            "vispa2_stats_import_report.html"
        )
    }
    if (!is.null(export_widget_path) && export_widget_path == "INTERNAL") {
        return(list(stats = stats, report_w = widget_stats_import))
    } else {
        return(stats)
    }
}

#' Import integration matrices based on the association file.
#'
#' @description \lifecycle{maturing}
#' These functions are designed to import the appropriate
#' integration matrix files given the association file and the root folder of
#' the file system where Vispa2 matrices are generated.
#' @details Import family functions are designed to work in combination with
#' Vispa2, for more details on this take a look here:
#' \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5702242/}{VISPA2:
#' A Scalable Pipeline for High-Throughput Identification
#' and Annotation of Vector Integration Sites}.\cr
#' For more details on how to properly use these functions, refer to the
#' vignette - vignette("how_to_import_functions")
#'
#' @section Interactive version:
#' The interactive version of import_parallel_Vispa2Matrices asks user for input
#' and allows a more detailed choice of projects to import, pools to import and,
#' if necessary, duplicate files. During the execution, a series of reports is
#' shown in html format.
#' @param association_file A single string containing the path to the
#' association file on disk, or a data frame resulting from a previous call to
#' `import_association_file`
#' @param quantification_type A vector of requested quantification_types. Must
#' be one in `quantification_types()`
#' @param matrix_type A single string representing the type of matrices to
#' be imported. Can only be one in `"annotated"` or `"not_annotated"`
#' @param workers A single integer representing the number
#' of parallel workers to use for the import
#' @param multi_quant_matrix If set to TRUE will produce a
#' multi-quantification matrix (data frame) through `comparison_matrix`
#' instead of a list.
#' @param export_report_path A path on disk to save produced import report
#'  or NULL if the user doesn't wish to save the html file
#' @param ... <[`dynamic-dots`][rlang::dyn-dots]> Additional named arguments
#' to pass to `ìmport_association_file` and `comparison_matrix`
#'
#' @seealso \code{\link{comparison_matrix}},
#' \code{\link{import_association_file}}
#'
#' @importFrom htmltools browsable tagList
#' @importFrom dplyr filter
#' @importFrom rlang dots_list inform abort call2 eval_tidy fn_fmls_names
#' @importFrom magrittr `%>%`
#' @family Import functions
#'
#' @return A named list of data frames containing data from
#' all imported integration
#' matrices, divided by quantification type or a multi-quantification matrix
#' @export
#'
#' @examples
#' \dontrun{
#' # Can't run because it's interactive and requires user input
#' matrices <- import_parallel_Vispa2Matrices_interactive(
#'     association_file,
#'     quantification_type,
#'     matrix_type = "annotated",
#'     workers = 2,
#'     multi_quant_matrix = FALSE,
#'     export_report_path = NULL,
#' )
#' }
import_parallel_Vispa2Matrices_interactive <- function(association_file,
    quantification_type,
    matrix_type = "annotated",
    workers = 2,
    multi_quant_matrix = TRUE,
    export_report_path = NULL,
    ...) {
    # Check parameters
    stopifnot((is.character(association_file) &
        length(association_file) == 1) ||
        is.data.frame(association_file))
    stopifnot(is.numeric(workers) & length(workers) == 1)
    workers <- floor(workers)
    stopifnot(!missing(quantification_type))
    stopifnot(all(quantification_type %in% quantification_types()))
    stopifnot(is.character(matrix_type) & matrix_type %in% c(
        "annotated",
        "not_annotated"
    ))
    stopifnot(is.logical(multi_quant_matrix) & length(multi_quant_matrix) == 1)
    ## Collect dot args
    if (is.character(association_file) || isTRUE(multi_quant_matrix)) {
        dots_args <- rlang::dots_list(..., .named = TRUE, .homonyms = "first")
        if (is.character(association_file)) {
            import_af_arg_names <- rlang::fn_fmls_names(import_association_file)
            import_af_arg_names <- import_af_arg_names[
                import_af_arg_names != "path"
            ]
            import_af_args <- dots_args[names(dots_args) %in%
                import_af_arg_names]
        }
        if (isTRUE(multi_quant_matrix)) {
            mult_arg_names <- rlang::fn_fmls_names(comparison_matrix)
            mult_arg_names <- mult_arg_names[mult_arg_names != "x"]
            mult_args <- dots_args[names(dots_args) %in%
                mult_arg_names]
        }
    }
    ## Import association file if provided a path
    if (is.character(association_file)) {
        association_file <- rlang::eval_tidy(
            rlang::call2("import_association_file",
                path = association_file,
                !!!import_af_args
            )
        )
    }
    ## Check there are the appropriate columns
    if (!"Path" %in% colnames(association_file)) {
        rlang::abort(.af_missing_path_error(), class = "missing_path_col")
    }
    association_file <- association_file %>%
        dplyr::filter(!is.na(.data$Path))
    if (nrow(association_file) == 0) {
        rlang::inform(c("The association file is empty, nothing to import",
            i = paste(
                "No projects left to import,",
                "absolute paths are all NA"
            )
        ))
        return(NULL)
    }
    ## User selects projects to keep
    association_file <- .interactive_select_projects_import(association_file)
    ## User selects pools to keep
    association_file <- .interactive_select_pools_import(association_file)
    ## Scan the appropriate file system paths and look for files
    files_found <- .lookup_matrices(
        association_file, quantification_type,
        matrix_type
    )
    ## Manage missing files and duplicates
    files_to_import <- .manage_anomalies_interactive(files_found)

    ## If files to import are 0 just terminate
    if (nrow(files_to_import) == 0) {
        rlang::abort("No files to import")
    }

    ## Import
    matrices <- .parallel_import_merge(files_to_import, workers)
    fimported <- matrices[[2]]
    if (getOption("ISAnalytics.widgets") == TRUE) {
        withCallingHandlers(
            {
                withRestarts(
                    {
                        import_widget <- .import_report_widget(
                            files_found,
                            files_to_import,
                            fimported
                        )
                        print(htmltools::browsable(htmltools::tagList(
                            import_widget
                        )))
                        if (!is.null(export_report_path)) {
                            .export_widget_file(
                                import_widget,
                                export_report_path,
                                "matrices_import_report.html"
                            )
                        }
                    },
                    print_err = function() {
                        rlang::inform(.widgets_error())
                        .summary_ism_par_import_msg(
                            fimported,
                            files_to_import,
                            files_found
                        )
                    }
                )
            },
            error = function(cnd) {
                rlang::inform(conditionMessage(cnd))
                invokeRestart("print_err")
            }
        )
    } else {
        .summary_ism_par_import_msg(
            fimported,
            files_to_import,
            files_found
        )
    }
    matrices <- matrices[[1]]
    if (multi_quant_matrix == TRUE) {
        matrices <- rlang::eval_tidy(rlang::call2(comparison_matrix,
            x = matrices,
            !!!mult_args
        ))
    }
    matrices
}


#' @inherit import_parallel_Vispa2Matrices_interactive
#'
#' @section Automatic version:
#' The automatic version of import_parallel_Vispa2Matrices doesn't interact with
#' the user directly, for this reason options in this modality are more limited
#' compared to the interactive version. In automatic version you can't:
#' * Choose single projects or pools: to have a selection import
#' the association file first and filter it according to your needs
#' before calling the function
#' (more details on this in the vignette)
#' * Choose duplicates: if, after filtering by the specified patterns,
#' duplicates are found they are automatically ignored
#'
#' @param patterns A character vector of additional patterns to match on file
#' names. Please note that patterns must be regular expressions. Can be NULL if
#' no patterns needs to be matched.
#' @param matching_opt A single value between `matching_options`
#' @seealso \code{\link{matching_options}},
#' \url{https://stringr.tidyverse.org/articles/regular-expressions.html}
#' @family Import functions
#' @importFrom htmltools browsable tagList
#' @importFrom dplyr filter
#' @importFrom rlang dots_list inform abort call2 eval_tidy fn_fmls_names
#' @importFrom magrittr `%>%`
#'
#' @export
#'
#' @examples
#' op <- options("ISAnalytics.widgets" = FALSE)
#' path <- system.file("extdata", "ex_association_file.tsv",
#'     package = "ISAnalytics"
#' )
#' root_pth <- system.file("extdata", "fs.zip", package = "ISAnalytics")
#' root <- unzip_file_system(root_pth, "fs")
#' matrices <- import_parallel_Vispa2Matrices_auto(
#'     association_file = path,
#'     quantification_type = c("fragmentEstimate", "seqCount"),
#'     patterns = NULL, matching_opt = "ANY",
#'     root = root,
#'     dates_format = "dmy",
#'     workers = 2
#' )
#' options(op)
import_parallel_Vispa2Matrices_auto <- function(association_file,
    quantification_type,
    matrix_type = "annotated",
    workers = 2,
    multi_quant_matrix = TRUE,
    export_report_path = NULL,
    patterns = NULL,
    matching_opt = matching_options(),
    ...) {
    # Check parameters
    stopifnot((is.character(association_file) &
        length(association_file) == 1) ||
        is.data.frame(association_file))
    stopifnot(is.numeric(workers) & length(workers) == 1)
    workers <- floor(workers)
    stopifnot(!missing(quantification_type))
    stopifnot(all(quantification_type %in% quantification_types()))
    stopifnot(is.character(matrix_type) & matrix_type %in% c(
        "annotated",
        "not_annotated"
    ))
    stopifnot(is.logical(multi_quant_matrix) & length(multi_quant_matrix) == 1)
    if (!is.null(patterns)) {
        stopifnot(is.character(patterns))
    }
    ### Evaluate matching_opt
    matching_option <- match.arg(matching_opt)
    stopifnot(is.character(matching_option))
    ## Collect dot args
    if (is.character(association_file) || isTRUE(multi_quant_matrix)) {
        dots_args <- rlang::dots_list(..., .named = TRUE, .homonyms = "first")
        if (is.character(association_file)) {
            import_af_arg_names <- rlang::fn_fmls_names(import_association_file)
            import_af_arg_names <- import_af_arg_names[
                import_af_arg_names != "path"
            ]
            import_af_args <- dots_args[names(dots_args) %in%
                import_af_arg_names]
        }
        if (isTRUE(multi_quant_matrix)) {
            mult_arg_names <- rlang::fn_fmls_names(comparison_matrix)
            mult_arg_names <- mult_arg_names[mult_arg_names != "x"]
            mult_args <- dots_args[names(dots_args) %in%
                mult_arg_names]
        }
    }
    ## Import association file if provided a path
    if (is.character(association_file)) {
        association_file <- rlang::eval_tidy(
            rlang::call2("import_association_file",
                path = association_file,
                !!!import_af_args
            )
        )
    }
    ## Check there are the appropriate columns
    if (!"Path" %in% colnames(association_file)) {
        rlang::abort(.af_missing_path_error(), class = "missing_path_col")
    }
    association_file <- association_file %>%
        dplyr::filter(!is.na(.data$Path))
    if (nrow(association_file) == 0) {
        rlang::inform(c("The association file is empty, nothing to import",
            i = paste(
                "No projects left to import,",
                "absolute paths are all NA"
            )
        ))
        return(NULL)
    }
    # Automatic workflow - limited options
    ## In automatic workflow all projects and pools contained in the association
    ## file are considered. If more precise selection is needed on this, user
    ## should use the interactive version or filter the association file
    ## appropriately before calling the function.

    ## Scan the appropriate file system paths and look for files
    files_found <- .lookup_matrices_auto(
        association_file, quantification_type,
        matrix_type, patterns, matching_option
    )
    ## Manage missing files and duplicates
    files_to_import <- .manage_anomalies_auto(files_found)
    ## If files to import are 0 just terminate
    if (nrow(files_to_import) == 0) {
        rlang::abort("No files to import")
    }

    ## Import
    matrices <- .parallel_import_merge(files_to_import, workers)
    fimported <- matrices[[2]]
    if (getOption("ISAnalytics.widgets") == TRUE) {
        withCallingHandlers(
            {
                withRestarts(
                    {
                        import_widget <- .import_report_widget(
                            files_found,
                            files_to_import,
                            fimported
                        )
                        print(htmltools::browsable(htmltools::tagList(
                            import_widget
                        )))
                        if (!is.null(export_report_path)) {
                            .export_widget_file(
                                import_widget,
                                export_report_path,
                                "matrices_import_report.html"
                            )
                        }
                    },
                    print_err = function() {
                        rlang::inform(.widgets_error())
                        .summary_ism_par_import_msg(
                            fimported,
                            files_to_import,
                            files_found
                        )
                    }
                )
            },
            error = function(cnd) {
                rlang::inform(conditionMessage(cnd))
                invokeRestart("print_err")
            }
        )
    } else {
        .summary_ism_par_import_msg(
            fimported,
            files_to_import,
            files_found
        )
    }
    matrices <- matrices[[1]]
    if (multi_quant_matrix == TRUE) {
        matrices <- rlang::eval_tidy(rlang::call2(comparison_matrix,
            x = matrices,
            !!!mult_args
        ))
    }
    matrices
}



#' Possible choices for the `quantification_type` parameter.
#'
#' These are all the possible values for the
#' `quantification_type` parameter in
#' `import_parallel_vispa2Matrices_interactive` and
#' `import_parallel_vispa2Matrices_auto`.
#'
#' @details The possible values are:
#' * fragmentEstimate
#' * seqCount
#' * barcodeCount
#' * cellCount
#' * ShsCount
#' @return A vector of characters for quantification types
#' @export
#' @seealso \code{\link{import_parallel_Vispa2Matrices_interactive}},
#' \code{\link{import_parallel_Vispa2Matrices_auto}}
#'
#' @family Import functions helpers
#'
#' @examples
#' quant_types <- quantification_types()
quantification_types <- function() {
    c(
        "fragmentEstimate", "seqCount",
        "barcodeCount", "cellCount",
        "ShsCount"
    )
}


#' Possible choices for the `matching_opt` parameter.
#'
#' These are all the possible values for the `matching_opt` parameter in
#' `import_parallel_vispa2Matrices_auto`.
#' @details The values "ANY", "ALL" and "OPTIONAL", represent how the patterns
#' should be matched, more specifically
#' * ANY = look only for files that match AT LEAST one of the
#' patterns specified
#' * ALL = look only for files that match ALL of the patterns specified
#' * OPTIONAL = look preferentially for files that match, in order, all
#' patterns or any pattern and if no match is found return what is found (keep
#' in mind that duplicates are discarded in automatic mode)
#' @return A vector of characters for matching_opt
#' @export
#' @family Import functions helpers
#' @seealso \code{\link{import_parallel_Vispa2Matrices_auto}}
#'
#' @examples
#' opts <- matching_options()
matching_options <- function() {
    c("ANY", "ALL", "OPTIONAL")
}


#' Possible choices for the `dates_format` parameter in
#' `import_association_file`,
#' `import_parallel_vispa2Matrices_interactive` and
#' `import_parallel_vispa2Matrices_auto`.
#'
#' All options correspond to `lubridate` functions:
#' * ymd: year, month, date
#' * ydm: year, day, month
#' * mdy: month, day, year
#' * myd: month, year, day
#' * dmy: day, month, year
#' * dym: day, year, month
#' * yq: year quantile
#'
#' **NOTE: use the same date format across the association file.**
#'
#' @return A character vector
#' @export
#' @seealso \code{\link{import_association_file}},
#' \code{\link{import_parallel_Vispa2Matrices_auto}}
#'
#' @examples
#' date_formats()
date_formats <- function() {
    c("ymd", "ydm", "mdy", "myd", "dmy", "dym", "yq")
}


#' Default regex prefixes for Vispa2 stats files.
#'
#' Note that each element is a regular expression.
#'
#' @return A character vector of regexes
#' @export
#'
#' @examples
#' default_iss_file_prefixes()
default_iss_file_prefixes <- function() {
    c("stats\\.sequence.", "stats\\.matrix.")
}
