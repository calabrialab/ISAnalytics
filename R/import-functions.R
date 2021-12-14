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
#' @param keep_excluded Keep the columns in `to_exclude` as additional
#' id columns?
#' @param separator The column delimiter used, defaults to `\t`
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
#' @details
#' For more details see the "How to use import functions" vignette:
#' \code{vignette("import_functions_howto", package = "ISAnalytics")}
#'
#' @export
#'
#' @examples
#' fs_path <- system.file("extdata", "fs.zip", package = "ISAnalytics")
#' fs <- unzip_file_system(fs_path, "fs")
#' matrix_path <- fs::path(
#'     fs,
#'     "PJ01",
#'     "quantification",
#'     "POOL01-1",
#'     "PJ01_POOL01-1_seqCount_matrix.no0.annotated.tsv.gz"
#' )
#' matrix <- import_single_Vispa2Matrix(matrix_path)
#' head(matrix)
import_single_Vispa2Matrix <- function(path,
    to_exclude = NULL,
    keep_excluded = FALSE,
    separator = "\t") {
    stopifnot(!missing(path) & is.character(path))
    stopifnot(is.null(to_exclude) || is.character(to_exclude))
    stopifnot(is.logical(keep_excluded))
    stopifnot(is.character(separator))
    if (!file.exists(path)) {
        not_found_msg <- paste("File not found at", path)
        rlang::abort(not_found_msg)
    }
    if (!fs::is_file(path)) {
        rlang::abort("Path exists but is not a file")
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
        if (!keep_excluded) {
            .read_with_fread(
                path = path, to_drop = to_exclude,
                df_type = df_type, annotated = is_annotated,
                sep = separator
            )
        } else {
            .read_with_fread(
                path = path, to_drop = NULL,
                df_type = df_type, annotated = is_annotated,
                sep = separator
            )
        }
    } else {
        if (!keep_excluded) {
            .read_with_readr(
                path = path, to_drop = to_exclude,
                df_type = df_type, annotated = is_annotated,
                sep = separator
            )
        } else {
            .read_with_readr(
                path = path, to_drop = NULL,
                df_type = df_type, annotated = is_annotated,
                sep = separator
            )
        }
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
            exportglobals = TRUE,
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
        if (!is.null(to_exclude)) {
            if (length(to_exclude) > 0 && keep_excluded) {
                id_vars <- c(id_vars, to_exclude)
            }
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
#' @description \lifecycle{stable}
#' Imports the association file and immediately performs a check on
#' the file system starting from the root to assess the alignment between the
#' two.
#' @param path The path on disk to the association file.
#' @param root The path on disk of the root folder of Vispa2 output or NULL.
#' See details.
#' @param tp_padding Timepoint padding, indicates the number of digits of the
#' "TimePoint" column once imported. Fills the content with 0s up to the length
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
#' @param convert_tp Should be time points be converted into months and
#' years?
#' @param report_path The path where the report file should be saved.
#' Can be a folder, a file or NULL if no report should be produced.
#' Defaults to `{user_home}/ISAnalytics_reports`.
#' @param ... Additional arguments to pass to
#' \code{\link{import_Vispa2_stats}}
#' @family Import functions
#' @return The data frame holding metadata
#' @details
#' If the `root` argument is set to `NULL` no file system alignment is
#' performed. This allows to import the basic file but it won't be
#' possible to perfom automated matrix and stats import.
#' For more details see the "How to use import functions" vignette:
#' \code{vignette("import_functions_howto", package = "ISAnalytics")}
#' @export
#'
#' @importFrom purrr map_lgl set_names is_empty
#' @importFrom rlang inform abort dots_list exec .data
#' @importFrom tibble add_column
#' @importFrom dplyr mutate if_else
#' @importFrom stringr str_pad
#' @importFrom magrittr `%>%`
#' @seealso \code{\link{date_formats}}
#' @examples
#' fs_path <- system.file("extdata", "fs.zip", package = "ISAnalytics")
#' fs <- unzip_file_system(fs_path, "fs")
#' af_path <- system.file("extdata", "asso.file.tsv.gz", package = "ISAnalytics")
#' af <- import_association_file(af_path, root = fs, report_path = NULL)
#' head(af)
import_association_file <- function(path,
    root = NULL,
    tp_padding = 4,
    dates_format = "ymd",
    separator = "\t",
    filter_for = NULL,
    import_iss = FALSE,
    convert_tp = TRUE,
    report_path = default_report_path(),
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
        rlang::abort(.no_stats_import_err())
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
    if (nrow(parsing_problems) == 0) {
        parsing_problems <- NULL
    }
    if (nrow(date_problems) == 0) {
        date_problems <- NULL
    }
    col_probs <- list(missing = NULL, non_standard = NULL)
    if (!.check_af_correctness(as_file)) {
        col_probs[["missing"]] <- association_file_columns()[
            !association_file_columns() %in% colnames(as_file)
        ]
    }
    non_standard <- colnames(as_file)[
        !colnames(as_file) %in% c(
            association_file_columns(), "Path",
            "Path_quant", "Path_iss"
        )
    ]
    if (!purrr::is_empty(non_standard)) {
        col_probs[["non_standard"]] <- non_standard
    }
    missing_dates <- purrr::map_lgl(date_columns_coll(), function(date_col) {
        any(is.na(as_file[[date_col]]))
    }) %>% purrr::set_names(date_columns_coll())
    missing_dates <- names(missing_dates)[missing_dates == TRUE]
    if (length(missing_dates) == 0) {
        missing_dates <- NULL
    }
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
    import_stats_rep <- NULL
    missing_stats_rep <- NULL
    if (import_iss) {
        dots <- rlang::dots_list(..., .named = TRUE)
        dots <- dots[!names(dots) %in% c(
            "association_file",
            "report_path",
            "join_with_af"
        )]
        stats <- NULL
        withCallingHandlers(
            {
                withRestarts(
                    {
                        stats <- rlang::exec(import_Vispa2_stats,
                            association_file = as_file,
                            report_path = "INTERNAL",
                            join_with_af = TRUE,
                            !!!dots
                        )
                    },
                    fail_stats = function() {
                        fail_stats_msg <- paste(
                            "Issues in importing stats",
                            "files, skipping"
                        )
                        rlang::inform(fail_stats_msg)
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
            if (!is.null(stats$report)) {
                import_stats_rep <- stats$report$import
                missing_stats_rep <- stats$report$miss
            }
        }
    }
    withCallingHandlers(
        {
            .produce_report("asso_file",
                params = list(
                    parsing_prob = parsing_problems,
                    dates_prob = date_problems,
                    col_prob = col_probs,
                    crit_nas = missing_dates,
                    fs_align = checks,
                    iss_stats = import_stats_rep,
                    iss_stats_miss = missing_stats_rep
                ),
                path = report_path
            )
        },
        error = function(cnd) {
            rest <- findRestart("report_fail")
            invokeRestart(rest, cnd)
        }
    )
    if (!getOption("ISAnalytics.reports") & getOption("ISAnalytics.verbose")) {
        summary_report <- .summary_af_import_msg(
            pars_prob = parsing_problems, dates_prob = date_problems,
            cols_prob = col_probs[[!is.null(col_probs)]],
            crit_na = missing_dates,
            checks = ifelse(is.null(checks),
                yes = "skipped",
                no = ifelse(any(!checks$Found),
                    "problems detected",
                    "no problems detected"
                )
            )
        )
        rlang::inform(summary_report, class = "summary_report")
    }
    as_file
}

#' Import Vispa2 stats given the aligned association file.
#'
#' \lifecycle{stable}
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
#' @param report_path The path where the report file should be saved.
#' Can be a folder, a file or NULL if no report should be produced.
#' Defaults to `{user_home}/ISAnalytics_reports`.
#'
#' @family Import functions
#' @importFrom rlang inform abort .data
#' @importFrom stats setNames
#' @importFrom dplyr left_join filter select all_of distinct
#'
#' @return A data frame
#' @export
#'
#' @examples
#' fs_path <- system.file("extdata", "fs.zip", package = "ISAnalytics")
#' fs <- unzip_file_system(fs_path, "fs")
#' af_path <- system.file("extdata", "asso.file.tsv.gz",
#'     package = "ISAnalytics"
#' )
#' af <- import_association_file(af_path,
#'     root = fs,
#'     import_iss = FALSE,
#'     report_path = NULL
#' )
#' stats_files <- import_Vispa2_stats(af,
#'     join_with_af = FALSE,
#'     report_path = NULL
#' )
#' head(stats_files)
import_Vispa2_stats <- function(association_file,
    file_prefixes = default_iss_file_prefixes(),
    join_with_af = TRUE,
    pool_col = "concatenatePoolIDSeqRun",
    report_path = default_report_path()) {
    ## Check param
    if (!is.data.frame(association_file)) {
        rlang::abort(.af_not_imported_err())
    }
    path_cols <- .path_cols_names()
    if (!path_cols$iss %in% colnames(association_file)) {
        rlang::abort(.af_not_aligned_err())
    }
    min_cols <- c(
        "ProjectID",
        "CompleteAmplificationID",
        pool_col,
        path_cols$iss,
        "TagSequence"
    )
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
    ## Import
    stats <- .import_stats_iss(
        association_file = association_file,
        prefixes = file_prefixes
    )
    report <- stats$report
    stats <- stats$stats
    ## - IF NO STATS IMPORTED (STATS ARE NULL)
    if (is.null(stats)) {
        if (getOption("ISAnalytics.verbose") == TRUE) {
            rlang::inform(.no_stat_files_imported())
        }
        if (report_path == "INTERNAL") {
            ## If function was called from import_association_file
            return(stats = association_file, report = report)
        }
        to_return <- if (join_with_af) {
            association_file
        } else {
            NULL
        }
        ## Produce report if it is requested
        withCallingHandlers(
            {
                .produce_report("vispa2_stats",
                    params = list(iss_stats = report),
                    path = report_path
                )
            },
            error = function(cnd) {
                rest <- findRestart("report_fail")
                invokeRestart(rest, cnd)
            }
        )
        ## Return nothing or the original af
        return(to_return)
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
        addit_columns <- c("SubjectID", "CellMarker", "Tissue", "TimePoint")
        addit_columns <- addit_columns[addit_columns %in%
            colnames(association_file)]
        missing_stats <- association_file %>%
            dplyr::filter(is.na(.data$RUN_NAME)) %>%
            dplyr::select(
                .data$ProjectID,
                dplyr::all_of(pool_col),
                .data$CompleteAmplificationID,
                .data$TagSequence,
                dplyr::all_of(addit_columns)
            ) %>%
            dplyr::distinct()
        if (report_path == "INTERNAL") {
            ## If function was called from import_association_file
            return(list(
                stats = association_file,
                report = list(
                    import = report,
                    miss = missing_stats
                )
            ))
        }
        ## If report was requested produce it (with missing df)
        withCallingHandlers(
            {
                .produce_report("vispa2_stats",
                    params = list(
                        iss_stats = report,
                        iss_stats_miss = missing_stats
                    ),
                    path = report_path
                )
            },
            error = function(cnd) {
                rest <- findRestart("report_fail")
                invokeRestart(rest, cnd)
            }
        )
        return(association_file)
    }
    ## If report was requested produce it
    withCallingHandlers(
        {
            .produce_report("vispa2_stats",
                params = list(iss_stats = report),
                path = report_path
            )
        },
        error = function(cnd) {
            rest <- findRestart("report_fail")
            invokeRestart(rest, cnd)
        }
    )
    return(stats)
}


#' Import integration matrices from paths in the association file.
#'
#' @description
#' \lifecycle{stable}
#' The function offers a convenient way of importing multiple integration
#' matrices in an automated or semi-automated way.
#' For more details see the "How to use import functions" vignette:
#' \code{vignette("import_functions_howto", package = "ISAnalytics")}
#'
#' @param association_file Data frame imported via
#' \link{import_association_file} (with file system alignment) or
#' a string containing the path to the association file on disk.
#' @param quantification_type A vector of requested quantification_types.
#' Possible choices are \link{quantification_types}
#' @param matrix_type A single string representing the type of matrices
#' to be imported. Can only be one in "annotated" or "not_annotated".
#' @param workers A single integer representing the number of parallel
#' workers to use for the import
#' @param multi_quant_matrix If set to TRUE will produce a
#' multi-quantification matrix through \link{comparison_matrix}
#' instead of a list.
#' @param report_path The path where the report file should be saved.
#' Can be a folder, a file or NULL if no report should be produced.
#' Defaults to `{user_home}/ISAnalytics_reports`.
#' @param patterns Relevant only if argument `mode` is set to `AUTO`.
#' A character vector of additional patterns to match on file
#' names. Please note that patterns must be regular expressions. Can be NULL if
#' no patterns need to be matched.
#' @param matching_opt Relevant only if argument `mode` is set to `AUTO`.
#' A single value between \link{matching_options}
#' @param mode A single value between `AUTO` and `INTERACTIVE`. If
#' INTERACTIVE, the function will ask for input from the user on console,
#' otherwise the process is fully automated (with limitations, see vignette).
#' @param ... <[`dynamic-dots`][rlang::dyn-dots]> Additional named arguments
#' to pass to `ìmport_association_file` and `comparison_matrix`
#'
#' @importFrom rlang fn_fmls_names dots_list arg_match inform abort
#' @importFrom rlang eval_tidy call2
#'
#' @family Import functions
#' @return Either a multi-quantification matrix or a list of integration
#' matrices
#' @export
#'
#' @examples
#' fs_path <- system.file("extdata", "fs.zip", package = "ISAnalytics")
#' fs <- unzip_file_system(fs_path, "fs")
#' af_path <- system.file("extdata", "asso.file.tsv.gz",
#'     package = "ISAnalytics"
#' )
#' af <- import_association_file(af_path,
#'     root = fs,
#'     import_iss = FALSE,
#'     report_path = NULL
#' )
#' matrices <- import_parallel_Vispa2Matrices(af,
#'     c("seqCount", "fragmentEstimate"),
#'     mode = "AUTO", report_path = NULL
#' )
#' head(matrices)
import_parallel_Vispa2Matrices <- function(association_file,
    quantification_type,
    matrix_type = "annotated",
    workers = 2,
    multi_quant_matrix = TRUE,
    report_path = default_report_path(),
    patterns = NULL,
    matching_opt = matching_options(),
    mode = c("AUTO", "INTERACTIVE"),
    ...) {
    .base_param_check(
        association_file, quantification_type, matrix_type,
        workers, multi_quant_matrix
    )
    rlang::arg_match(mode)
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
    association_file <- .pre_manage_af(association_file, import_af_args)
    if (nrow(association_file) == 0) {
        rlang::inform(.af_empty_msg())
        return(NULL)
    }
    ## Workflows
    ### --- Interactive
    if (mode == "INTERACTIVE") {
        ## User selects projects to keep
        association_file <- .interactive_select_projects_import(
            association_file
        )
        ## User selects pools to keep
        association_file <- .interactive_select_pools_import(association_file)
        ## Scan the appropriate file system paths and look for files
        files_found <- .lookup_matrices(
            association_file, quantification_type,
            matrix_type
        )
        ## Manage missing files and duplicates
        files_to_import <- .manage_anomalies_interactive(files_found)
    } else {
        ### --- Auto
        ## In automatic workflow all projects and pools contained
        ## in the association
        ## file are considered. If more precise selection is needed on this,
        ## user
        ## should use the interactive version or filter the association file
        ## appropriately before calling the function.
        ### Evaluate patterns
        stopifnot(is.logical(multi_quant_matrix) &
            length(multi_quant_matrix) == 1)
        if (!is.null(patterns)) {
            stopifnot(is.character(patterns))
        }
        ### Evaluate matching_opt
        matching_option <- match.arg(matching_opt)
        stopifnot(is.character(matching_option))
        ## Scan the appropriate file system paths and look for files
        files_found <- .lookup_matrices_auto(
            association_file, quantification_type,
            matrix_type, patterns, matching_option
        )
        ## Manage missing files and duplicates
        files_to_import <- .manage_anomalies_auto(files_found)
    }
    ## If files to import are 0 just terminate
    if (nrow(files_to_import) == 0) {
        rlang::abort("No files to import")
    }
    ## Import
    matrices <- .parallel_import_merge(files_to_import, workers)
    fimported <- matrices[[2]]
    if (nrow(fimported) == 0) {
        fimported <- NULL
    }
    matrices <- matrices[[1]]
    if (multi_quant_matrix == TRUE) {
        matrices <- rlang::eval_tidy(rlang::call2(comparison_matrix,
            x = matrices,
            !!!mult_args
        ))
    }
    annotation_problems <- if (getOption("ISAnalytics.reports") == TRUE &
        !is.null(report_path)) {
        tmp <- if (!multi_quant_matrix) {
            comparison_matrix(matrices)
        } else {
            matrices
        }
        annotation_issues(tmp)
    } else {
        NULL
    }
    withCallingHandlers(
        {
            .produce_report("matrix_imp",
                params = list(
                    files_found = files_found,
                    files_imp = fimported,
                    annot_prob = annotation_problems
                ),
                path = report_path
            )
        },
        error = function(cnd) {
            rest <- findRestart("report_fail")
            invokeRestart(rest, cnd)
        }
    )
    return(matrices)
}


#' Import integration matrices from association file.
#'
#' @description `r lifecycle::badge("deprecated")`
#' This function was deprecated to avoid redundancy.
#' Please refer to \code{\link{import_parallel_Vispa2Matrices}}.
#'
#' @importFrom lifecycle deprecate_warn
#' @importFrom rlang list2 `!!!`
#'
#' @export
#' @keywords internal
#' @return A data frame or a list
import_parallel_Vispa2Matrices_interactive <- function(association_file,
    quantification_type,
    matrix_type = "annotated",
    workers = 2,
    multi_quant_matrix = TRUE,
    export_report_path = NULL,
    ...) {
    lifecycle::deprecate_warn(
        when = "1.3.3",
        what = "import_parallel_Vispa2Matrices_interactive()",
        with = "import_parallel_Vispa2Matrices()",
        details = c(paste0(
            "Use import_parallel",
            "_Vispa2Matrices(mode = 'INTERACTIVE') ",
            "for interactive mode or ",
            "import_parallel",
            "_Vispa2Matrices(mode = 'AUTO') ",
            "for automatic mode"
        ))
    )
    dots <- rlang::list2(...)
    import_parallel_Vispa2Matrices(association_file,
        quantification_type,
        matrix_type = "annotated",
        workers = 2,
        multi_quant_matrix = TRUE,
        report_path = export_report_path,
        mode = "INTERACTIVE",
        !!!dots
    )
}

#' Import integration matrices from association file.
#'
#' @description `r lifecycle::badge("deprecated")`
#' This function was deprecated to avoid redundancy.
#' Please refer to \code{\link{import_parallel_Vispa2Matrices}}.
#' @importFrom lifecycle deprecate_warn
#' @importFrom rlang list2 `!!!`
#' @export
#' @keywords internal
#' @return A data frame or a list
import_parallel_Vispa2Matrices_auto <- function(association_file,
    quantification_type,
    matrix_type = "annotated",
    workers = 2,
    multi_quant_matrix = TRUE,
    patterns = NULL,
    matching_opt = matching_options(),
    export_report_path = NULL,
    ...) {
    lifecycle::deprecate_warn(
        when = "1.3.3",
        what = "import_parallel_Vispa2Matrices_auto()",
        with = "import_parallel_Vispa2Matrices()",
        details = c(paste0(
            "Use import_parallel",
            "_Vispa2Matrices(mode = 'INTERACTIVE') ",
            "for interactive mode or ",
            "import_parallel",
            "_Vispa2Matrices(mode = 'AUTO') ",
            "for automatic mode"
        ))
    )
    dots <- rlang::list2(...)
    import_parallel_Vispa2Matrices(association_file,
        quantification_type,
        matrix_type = "annotated",
        workers = 2,
        multi_quant_matrix = TRUE,
        mode = "AUTO",
        patterns = patterns,
        matching_opt = matching_opt,
        report_path = export_report_path,
        !!!dots
    )
}


#' Check for genomic annotation problems in IS matrices.
#'
#' \lifecycle{experimental}
#' This helper function checks if each individual integration site,
#' identified by the triplet (chr, integration locus, strand),
#' has been annotated with two or more distinct gene symbols.
#'
#' @param matrix Either a single matrix or a list of matrices, ideally obtained
#' via `import_parallel_Vispa2Matrices()` or `import_single_Vispa2Matrix()`
#'
#' @return Either `NULL` if no issues were detected or 1 or more data frames
#' with genomic coordinates of the IS and the number of distinct
#' genes associated
#' @export
#'
#' @family Import functions helpers
#'
#' @examples
#' data("integration_matrices", package = "ISAnalytics")
#' annotation_issues(integration_matrices)
annotation_issues <- function(matrix) {
    stopifnot(is.list(matrix))
    find_probs <- function(m) {
        needed <- c(mandatory_IS_vars(), annotation_IS_vars())
        if (!all(needed %in% colnames(m))) {
            rlang::inform(.missing_needed_cols(needed[!needed %in% colnames(m)]))
            return(NULL)
        }
        tmp <- m %>%
            dplyr::select(dplyr::all_of(c(
                mandatory_IS_vars(),
                annotation_IS_vars()
            ))) %>%
            dplyr::distinct() %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(mandatory_IS_vars()))) %>%
            dplyr::summarise(distinct_genes = dplyr::n())
        if (any(tmp$distinct_genes > 1)) {
            tmp %>%
                dplyr::filter(.data$distinct_genes > 1)
        } else {
            NULL
        }
    }
    if (is.data.frame(matrix)) {
        probs <- find_probs(matrix)
        if (is.null(probs) & getOption("ISAnalytics.verbose") == TRUE) {
            rlang::inform("No annotation issues found")
        }
        return(probs)
    }
    probs <- purrr::map(matrix, find_probs)
    if (all(is.null(probs)) & getOption("ISAnalytics.verbose") == TRUE) {
        rlang::inform("No annotation issues found")
    }
    return(probs)
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
