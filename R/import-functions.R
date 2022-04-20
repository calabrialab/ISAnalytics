#------------------------------------------------------------------------------#
# Importing functions
#------------------------------------------------------------------------------#
#' Import a single integration matrix from file
#'
#' @description
#' `r lifecycle::badge("stable")`
#' This function allows to read and import an integration matrix
#' (ideally produced by VISPA2) and converts it to a tidy
#' format.
#'
#' @param path The path to the file on disk
#' @param separator The column delimiter used, defaults to `\t`
#' @param additional_cols Either `NULL`, a named character vector or a named
#' list. See details.
#' @param sample_names_to Name of the output column holding the sample
#' identifier. Defaults to `pcr_id_column()`
#' @param values_to Name of the output column holding the quantification
#' values. Defaults to `Value`.
#' @param to_exclude `r lifecycle::badge("deprecated")`
#' Deprecated. Use `additonal_cols` instead
#' @param keep_excluded `r lifecycle::badge("deprecated")`
#' Deprecated. Use `additonal_cols` instead
#'
#' @details
#' ## Additional columns
#' Additional columns are annotation columns present in the integration matrix
#' to import that are not
#' * part of the mandatory IS vars (see `mandatory_IS_vars()`)
#' * part of the annotation IS vars (see `annotation_IS_vars()`)
#' * the sample identifier column
#' * the quantification column
#'
#' When specified they tell the function how to treat those columns in the
#' import phase, by providing a named character vector, where names correspond
#' to the additional column names and values are a choice of the following:
#'
#' * `"char"` for character (strings)
#' * `"int"` for integers
#' * `"logi"` for logical values (TRUE / FALSE)
#' * `"numeric"` for numeric values
#' * `"factor"` for factors
#' * `"date"` for generic date format - note that functions that
#' need to read and parse files will try to guess the format and parsing
#' may fail
#' * One of the accepted date/datetime formats by `lubridate`,
#' you can use `ISAnalytics::date_formats()` to view the accepted formats
#' * `"_"` to drop the column
#'
#' For more details see the "How to use import functions" vignette:
#' \code{vignette("import_functions_howto", package = "ISAnalytics")}
#'
#'
#' @template transform_list
#'
#' @section Required tags:
#' The function will explicitly check for the presence of these tags:
#'
#' * All columns declared in `mandatory_IS_vars()`
#'
#' @return A data.table object in tidy format
#'
#' @family Import functions
#' @export
#'
#' @examples
#' fs_path <- generate_default_folder_structure(type = "correct")
#' matrix_path <- fs::path(
#'     fs_path$root, "PJ01", "quantification",
#'     "POOL01-1", "PJ01_POOL01-1_seqCount_matrix.no0.annotated.tsv.gz"
#' )
#' matrix <- import_single_Vispa2Matrix(matrix_path)
#' head(matrix)
import_single_Vispa2Matrix <- function(path,
    separator = "\t",
    additional_cols = NULL,
    transformations = NULL,
    sample_names_to = pcr_id_column(),
    values_to = "Value",
    to_exclude = lifecycle::deprecated(),
    keep_excluded = lifecycle::deprecated()) {
    stopifnot(!missing(path) & is.character(path))
    stopifnot(is.character(separator))
    if (!file.exists(path)) {
        not_found_msg <- paste("File not found at", path)
        rlang::abort(not_found_msg)
    }
    if (!fs::is_file(path)) {
        rlang::abort("Path exists but is not a file")
    }
    stopifnot(is.null(transformations) ||
        (is.list(transformations) && !is.null(names(transformations))))
    stopifnot(is.character(sample_names_to))
    sample_names_to <- sample_names_to[1]
    stopifnot(is.character(values_to))
    values_to <- values_to[1]
    deprecation_details <- paste(
        "Arguments 'to_exclude' and 'keep_excluded'",
        "are deprecated in favor of a single argument",
        "which allows a more refined tuning. See",
        "`?import_single_Vispa2Matrix` for details"
    )
    if (lifecycle::is_present(to_exclude)) {
        lifecycle::deprecate_stop(
            when = "1.5.4",
            what = "import_single_Vispa2Matrix(to_exclude)",
            with = "import_single_Vispa2Matrix(additional_cols)",
            details = deprecation_details
        )
        return(NULL)
    }
    if (lifecycle::is_present(keep_excluded)) {
        lifecycle::deprecate_stop(
            when = "1.5.4",
            what = "import_single_Vispa2Matrix(keep_excluded)",
            with = "import_single_Vispa2Matrix(additional_cols)",
            details = deprecation_details
        )
        return(NULL)
    }

    tidy_df <- .import_single_matrix(
        path = path, separator = separator,
        additional_cols = additional_cols,
        transformations = transformations,
        call_mode = "EXTERNAL",
        id_col_name = sample_names_to,
        val_col_name = values_to
    )
    return(tidy_df)
}


#' Import the association file from disk
#'
#' @description
#' `r lifecycle::badge("stable")`
#' Imports the association file and optionally performs a check on
#' the file system starting from the root to assess the alignment between the
#' two.
#'
#' @param path The path on disk to the association file.
#' @param root The path on disk of the root folder of VISPA2 output or `NULL`.
#' See details.
#' @param dates_format A single string indicating how dates should be parsed.
#' Must be a value in: \code{date_formats()}
#' @param separator The column separator used in the file
#' @param filter_for A named list where names represent column names that
#' must be filtered. For example: `list(ProjectID = c("PROJECT1", "PROJECT2))`
#' will filter the association file so that it contains only those rows
#' for which the value of the column "ProjectID" is one of the specified
#' values. If multiple columns are present in the list all filtering
#' conditions are applied as a logical AND.
#' @param import_iss Import VISPA2 pool stats and merge them with the
#' association file? Logical value
#' @param convert_tp Should be time points be converted into months and
#' years? Logical value
#' @param ... Additional arguments to pass to
#' \code{\link{import_Vispa2_stats}}
#' @param tp_padding `r lifecycle::badge("deprecated")` Deprecated.
#' Use `transformations` instead.
#'
#' @template transform_list
#' @template report_path_param
#'
#' @family Import functions
#' @return The data frame containing metadata
#' @details
#' ## File system alignment
#' If the `root` argument is set to `NULL` no file system alignment is
#' performed. This allows to import the basic file but it won't be
#' possible to perform automated matrix and stats import.
#' For more details see the "How to use import functions" vignette:
#' \code{vignette("import_functions_howto", package = "ISAnalytics")}
#'
#' ## Time point conversion
#' The time point conversion is based on the following logic, given `TPD`
#' is the column containing the time point expressed in days and
#' `TPM` and `TPY` are respectively the time points expressed as month
#' and years
#' - If `TPD` is `NA` --> `NA` (for both months and years)
#' - `TPM` = 0, `TPY` = 0 if and only if `TPD` = 0
#' For conversion in months:
#' - `TPM` = ceiling(`TPD`/30) if `TPD` < 30 otherwise `TPM` = round(`TPD`/30)
#' For conversion in years:
#' - `TPY` = ceiling(`TPD`/360)
#'
#' @section Required tags:
#' The function will explicitly check for the presence of these tags:
#'
#' ```{r echo=FALSE, results="asis"}
#' all_tags <- available_tags()
#' af_needed <- unique(all_tags[purrr::map_lgl(eval(rlang::sym("needed_in")),
#'  ~ "import_association_file" %in% .x)][["tag"]])
#'  cat(paste0("* ", af_needed, collapse="\n"))
#' ```
#'
#' The function will use all the available specifications contained in
#' `association_file_columns(TRUE)` to read and parse the file.
#' If the specifications contain columns with a type `"date"`, the function
#' will parse the generic date with the format in the `dates_format` argument.
#'
#' @export
#'
#' @seealso \code{\link{date_formats}}
#' @examples
#' fs_path <- generate_default_folder_structure(type = "correct")
#' af <- import_association_file(fs_path$af,
#'     root = fs_path$root,
#'     report_path = NULL
#' )
#' head(af)
import_association_file <- function(path,
    root = NULL,
    dates_format = "ymd",
    separator = "\t",
    filter_for = NULL,
    import_iss = FALSE,
    convert_tp = TRUE,
    report_path = default_report_path(),
    transformations = default_af_transform(convert_tp),
    tp_padding = lifecycle::deprecated(),
    ...) {
    # Check parameters
    stopifnot(is.character(path))
    path <- path[1]
    stopifnot((is.character(root) & length(root) == 1) || (is.null(root)))
    stopifnot(file.exists(path))
    if (!is.null(root) && root != "") {
        stopifnot(file.exists(root))
    }
    stopifnot(length(dates_format) == 1 & dates_format %in% date_formats())
    stopifnot(is.character(separator))
    separator <- separator[1]
    stopifnot(is.logical(import_iss))
    import_iss <- import_iss[1]
    if (import_iss & is.null(root)) {
        rlang::abort(.no_stats_import_err())
    }
    stopifnot(is.null(transformations) ||
        (is.list(transformations) && !is.null(names(transformations))))
    if (lifecycle::is_present(tp_padding)) {
        lifecycle::deprecate_warn(
            when = "1.5.4",
            what = "import_association_file(tp_padding)",
            details = c(paste(
                "The argument is now deprecated in favor of custom",
                "column transformations"
            ),
            i = paste(
                "See the documentation of `transform_columns`",
                "or browse the package vignettes for more details"
            )
            )
        )
    }
    # Check filter
    stopifnot(is.null(filter_for) ||
        (is.list(filter_for) && !is.null(names(filter_for))))

    # Check presence of required tags
    required_tags <- list()
    req_tags_politic <- c()
    if (!is.null(root)) {
        ### Tags required for file system alignment
        required_tags <- append(
            required_tags,
            list(
                "project_id" = "char",
                "proj_folder" = "char",
                "vispa_concatenate" = "char"
            )
        )
        req_tags_politic <- c(req_tags_politic,
            "project_id" = "error",
            "proj_folder" = "error",
            "vispa_concatenate" = "error"
        )
    }
    if (convert_tp) {
        ### tags required for time point conversion
        required_tags <- append(
            required_tags,
            list("tp_days" = c("char", "int", "numeric"))
        )
        req_tags_politic <- c(req_tags_politic,
            "tp_days" = "first"
        )
    }
    tags_to_cols <- if (!purrr::is_empty(required_tags)) {
        .check_required_cols(
            required_tags = required_tags,
            vars_df = association_file_columns(TRUE),
            duplicate_politic = req_tags_politic
        )
    } else {
        NULL
    }
    # Read file and check the correctness
    af_checks <- .manage_association_file(
        af_path = path,
        root = root,
        format = dates_format,
        delimiter = separator,
        filter = filter_for,
        proj_fold_col = dplyr::if_else(!is.null(tags_to_cols),
            tags_to_cols %>%
                dplyr::filter(.data$tag == "proj_folder") %>%
                dplyr::pull(.data$names),
            NULL
        ),
        concat_pool_col = dplyr::if_else(!is.null(tags_to_cols),
            tags_to_cols %>%
                dplyr::filter(.data$tag == "vispa_concatenate") %>%
                dplyr::pull(.data$names),
            NULL
        ),
        project_id_col = dplyr::if_else(!is.null(tags_to_cols),
            tags_to_cols %>%
                dplyr::filter(.data$tag == "project_id") %>%
                dplyr::pull(.data$names),
            NULL
        )
    )
    as_file <- af_checks$af
    parsing_problems <- af_checks$parsing_probs
    date_problems <- af_checks$date_probs
    checks <- af_checks$check
    if (nrow(parsing_problems) == 0) {
        parsing_problems <- NULL
    }
    if (is.null(date_problems) || nrow(date_problems) == 0) {
        date_problems <- NULL
    }
    col_probs <- list(missing = NULL, non_standard = NULL)
    if (!.check_af_correctness(as_file)) {
        min_required_cols <- association_file_columns(TRUE) %>%
            dplyr::filter(.data$flag == "required") %>%
            dplyr::pull(.data$names)
        col_probs[["missing"]] <- min_required_cols[
            !min_required_cols %in% colnames(as_file)
        ]
    }
    non_standard <- colnames(as_file)[
        !colnames(as_file) %in% c(
            association_file_columns(), .path_cols_names()
        )
    ]
    if (!purrr::is_empty(non_standard)) {
        col_probs[["non_standard"]] <- non_standard
    }
    ## Fix timepoints
    if (convert_tp) {
        tp_col <- tags_to_cols %>%
            dplyr::filter(.data$tag == "tp_days") %>%
            dplyr::pull(.data$names)
        as_file <- as_file %>%
            dplyr::mutate(
                TimepointMonths = .timepoint_to_months(.data[[tp_col]]),
                TimepointYears = .timepoint_to_years(.data[[tp_col]])
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

    if (!is.null(transformations)) {
        as_file <- transform_columns(as_file, transformations)
    }
    crit_tags <- c(
        "project_id", "pool_id", "tag_seq", "subject", "tissue",
        "cell_marker", "pcr_replicate", "vispa_concatenate",
        "pcr_repl_id", "proj_folder"
    )
    crit_colnames <- association_file_columns(TRUE) %>%
        dplyr::filter(.data$tag %in% crit_tags) %>%
        dplyr::pull(.data$names)
    crit_colnames <- colnames(as_file)[colnames(as_file) %in% crit_colnames]
    crit_nas <- if (length(crit_colnames) > 0) {
        nas_crit <- purrr::map_lgl(crit_colnames, ~ {
            any(is.na(as_file[[.x]]))
        }) %>%
            purrr::set_names(crit_colnames)
        nas_crit <- names(purrr::keep(nas_crit, ~ .x == TRUE))
        if (length(nas_crit) == 0) {
            NULL
        } else {
            nas_crit
        }
    } else {
        NULL
    }
    withCallingHandlers(
        {
            .produce_report("asso_file",
                params = list(
                    parsing_prob = parsing_problems,
                    dates_prob = date_problems,
                    col_prob = col_probs,
                    crit_nas = crit_nas,
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
            crit_na = crit_nas,
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
#' @description
#' `r lifecycle::badge("stable")`
#' Imports all the Vispa2 stats files for each pool provided the association
#' file has been aligned with the file system
#' (see \code{\link{import_association_file}}).
#'
#' @param association_file The file system aligned association file
#' (contains columns with absolute paths to the 'iss' folder)
#' @param file_prefixes A character vector with known file prefixes
#'  to match on file names. NOTE: the elements represent regular expressions.
#'  For defaults see \link{default_iss_file_prefixes}.
#' @param join_with_af Logical, if `TRUE` the imported stats files will be
#' merged with the association file, if `FALSE` a single data frame holding
#' only the stats will be returned.
#' @param pool_col A single string. What is the name of the pool column
#' used in the Vispa2 run? This will be used as a key to perform a join
#' operation with the stats files `POOL` column.
#'
#' @template report_path_param
#' @section Required tags:
#' The function will explicitly check for the presence of these tags:
#'
#' ```{r echo=FALSE, results="asis"}
#' all_tags <- available_tags()
#' needed <- unique(all_tags[purrr::map_lgl(eval(rlang::sym("needed_in")),
#'  ~ "import_Vispa2_stats" %in% .x)][["tag"]])
#'  cat(paste0("* ", needed, collapse="\n"))
#' ```
#'
#'
#' @family Import functions
#' @importFrom rlang inform abort .data
#' @importFrom stats setNames
#'
#' @return A data frame
#' @export
#'
#' @examples
#' fs_path <- generate_default_folder_structure(type = "correct")
#' af <- import_association_file(fs_path$af,
#'     root = fs_path$root,
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
    required_tags <- list(
        "project_id" = "char",
        "tag_seq" = "char",
        "pcr_repl_id" = "char"
    )
    tag_politics <- list(
        "project_id" = "error",
        "tag_seq" = "error",
        "pcr_repl_id" = "error"
    )
    tags_to_cols <- .check_required_cols(
        required_tags = required_tags,
        vars_df = association_file_columns(TRUE),
        duplicate_politic = tag_politics
    )
    min_cols <- c(tags_to_cols$names, pool_col, path_cols$iss)
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
        prefixes = file_prefixes,
        pool_col = pool_col,
        path_iss_col = path_cols$iss,
        tags = tags_to_cols
    )
    report <- stats$report
    stats <- stats$stats
    ## - IF NO STATS IMPORTED (STATS ARE NULL)
    if (is.null(stats)) {
        if (getOption("ISAnalytics.verbose") == TRUE) {
            rlang::inform(.no_stat_files_imported())
        }
        if (!is.null(report_path) && report_path == "INTERNAL") {
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
        required_tags_for_join <- list(
            "tag_seq" = "char",
            "vispa_concatenate" = "char"
        )
        iss_tags_to_cols <- .check_required_cols(required_tags_for_join,
            iss_stats_specs(TRUE),
            duplicate_politic = "error"
        )
        if (any(!iss_tags_to_cols$names %in% colnames(stats))) {
            msg <- c("Error in joining VISPA2 stats with AF - skipping",
                .missing_needed_cols(iss_tags_to_cols$names[
                    !iss_tags_to_cols$names %in% colnames(stats)
                ]),
                i = paste(
                    "Needed columns for join were not found in",
                    "imported stats. Check your iss stats specs",
                    "with `iss_stats_specs(TRUE)` and check the files",
                    "are not malformed."
                ),
                i = paste("Returning imported data only")
            )
            rlang::inform(msg, class = "iss_join_missing")
            return(stats)
        }
        iss_pool_col <- iss_tags_to_cols %>%
            dplyr::filter(.data$tag == "vispa_concatenate") %>%
            dplyr::pull(.data$names)
        iss_tag_col <- iss_tags_to_cols %>%
            dplyr::filter(.data$tag == "tag_seq") %>%
            dplyr::pull(.data$names)
        association_file <- association_file %>%
            dplyr::left_join(stats, by = c(
                stats::setNames(iss_pool_col, pool_col),
                stats::setNames(iss_tag_col, tags_to_cols %>%
                    dplyr::filter(.data$tag == "tag_seq") %>%
                    dplyr::pull(.data$names))
            ))
        ## Detect potential problems
        addit_columns <- association_file_columns(TRUE) %>%
            dplyr::filter(.data$tag %in% c(
                "subject",
                "tissue",
                "cell_marker",
                "tp_days"
            ))
        addit_columns_names <- addit_columns %>%
            dplyr::pull(.data$names)
        addit_columns_names <- addit_columns_names[addit_columns_names %in%
            colnames(association_file)]
        iss_cols_in_af <- colnames(association_file)[
            colnames(association_file) %in% iss_stats_specs()
        ]

        missing_stats <- association_file %>%
            dplyr::filter(dplyr::if_all(
                dplyr::all_of(iss_cols_in_af),
                is.na
            )) %>%
            dplyr::select(dplyr::all_of(c(
                tags_to_cols$names,
                pool_col,
                addit_columns_names
            ))) %>%
            dplyr::distinct()
        all_af_tags <- tags_to_cols %>%
            dplyr::bind_rows(addit_columns)
        if (!is.null(report_path) && report_path == "INTERNAL") {
            ## If function was called from import_association_file
            return(list(
                stats = association_file,
                report = list(
                    import = report,
                    miss = missing_stats,
                    af_tag_map = all_af_tags
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
#' `r lifecycle::badge("stable")`
#' The function offers a convenient way of importing multiple integration
#' matrices in an automated or semi-automated way.
#' For more details see the "How to use import functions" vignette:
#' \code{vignette("import_functions_howto", package = "ISAnalytics")}
#'
#' @section Required tags:
#' The function will explicitly check for the presence of these tags:
#'
#' ```{r echo=FALSE, results="asis"}
#' all_tags <- available_tags()
#' needed <- unique(all_tags[purrr::map_lgl(eval(rlang::sym("needed_in")),
#'  ~ "import_parallel_Vispa2Matrices" %in% .x)][["tag"]])
#'  cat(paste0("* ", needed, collapse="\n"))
#' ```
#'
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
#' @param multi_quant_matrix If set to `TRUE` will produce a
#' multi-quantification matrix through \link{comparison_matrix}
#' instead of a list.
#' @param patterns Relevant only if argument `mode` is set to `AUTO`.
#' A character vector of additional patterns to match on file
#' names. Please note that patterns must be regular expressions. Can be `NULL`
#' if no patterns need to be matched.
#' @param matching_opt Relevant only if argument `mode` is set to `AUTO`.
#' A single value between \link{matching_options}
#' @param mode A single value between `AUTO` and `INTERACTIVE`. If
#' `INTERACTIVE`, the function will ask for input from the user on console,
#' otherwise the process is fully automated (with limitations, see vignette).
#' @param ... <[`dynamic-dots`][rlang::dyn-dots]> Additional named arguments
#' to pass to `ìmport_association_file`, `comparison_matrix` and
#' `import_single_Vispa2_matrix`
#'
#' @template report_path_param
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
#' fs_path <- generate_default_folder_structure(type = "correct")
#' af <- import_association_file(fs_path$af,
#'     root = fs_path$root,
#'     report_path = NULL
#' )
#' matrices <- import_parallel_Vispa2Matrices(af,
#'     c("seqCount", "fragmentEstimate"),
#'     mode = "AUTO", report_path = NULL
#' )
#' head(matrices)
import_parallel_Vispa2Matrices <- function(association_file,
    quantification_type = c("seqCount", "fragmentEstimate"),
    matrix_type = c("annotated", "not_annotated"),
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
    matrix_type <- rlang::arg_match(matrix_type)
    mode <- rlang::arg_match(mode)
    ## Collect dot args
    if (is.character(association_file) || isTRUE(multi_quant_matrix)) {
        dots_args <- rlang::dots_list(..., .named = TRUE, .homonyms = "first")
        if (is.character(association_file)) {
            import_af_arg_names <- rlang::fn_fmls_names(import_association_file)
            import_af_arg_names <- import_af_arg_names[
                !import_af_arg_names %in% c("path", "report_path")
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
    import_matrix_arg_names <- rlang::fn_fmls_names(
        import_single_Vispa2Matrix
    )
    import_matrix_arg_names <- import_matrix_arg_names[
        !import_matrix_arg_names %in% c(
            "path", "to_exclude",
            "keep_excluded"
        )
    ]
    import_matrix_args <- dots_args[names(dots_args) %in%
        import_matrix_arg_names]
    association_file <- .pre_manage_af(
        association_file,
        import_af_args,
        report_path
    )
    if (nrow(association_file) == 0) {
        rlang::inform(.af_empty_msg())
        return(NULL)
    }
    ## Workflows
    af_tags <- association_file_columns(TRUE)
    proj_col <- af_tags %>%
        dplyr::filter(.data$tag == "project_id") %>%
        dplyr::pull(.data$names)
    pool_col <- af_tags %>%
        dplyr::filter(.data$tag == "vispa_concatenate") %>%
        dplyr::pull(.data$names)
    ### --- Interactive
    if (mode == "INTERACTIVE") {
        matching_option <- NULL
        ## User selects projects to keep
        association_file <- .interactive_select_projects_import(
            association_file,
            proj_col = proj_col
        )
        ## User selects pools to keep
        association_file <- .interactive_select_pools_import(association_file,
            proj_col = proj_col,
            pool_col = pool_col
        )
        ## Scan the appropriate file system paths and look for files
        files_found <- .lookup_matrices(
            association_file, quantification_type,
            matrix_type, proj_col, pool_col
        )
        ## Manage missing files and duplicates
        files_to_import <- .manage_anomalies_interactive(
            files_found,
            proj_col,
            pool_col
        )
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
        matching_option <- rlang::arg_match(matching_opt)
        stopifnot(is.character(matching_option))
        ## Scan the appropriate file system paths and look for files
        files_found <- .lookup_matrices_auto(
            association_file, quantification_type,
            matrix_type, patterns, matching_option,
            proj_col, pool_col
        )
        ## Manage missing files and duplicates
        files_to_import <- .manage_anomalies_auto(
            files_found,
            proj_col, pool_col
        )
    }
    ## If files to import are 0 just terminate
    if (nrow(files_to_import) == 0) {
        rlang::abort("No files to import")
    }
    ## Import
    matrices <- .parallel_import_merge(
        files_to_import, workers,
        import_matrix_args
    )
    fimported <- matrices$summary
    if (nrow(fimported) == 0) {
        fimported <- NULL
    }
    matrices <- matrices$matrix
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
    launch_params <- list()
    if (!is.null(patterns)) {
        launch_params[["patterns"]] <- patterns
        launch_params[["matching_opt"]] <- matching_option
    }
    withCallingHandlers(
        {
            .produce_report("matrix_imp",
                params = list(
                    launch_params = launch_params,
                    set_vars = list(proj_col = proj_col, pool_col = pool_col),
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
#' @description `r lifecycle::badge("experimental")`
#' This helper function checks if each individual integration site,
#' identified by the `mandatory_IS_vars()`,
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
            rlang::inform(
                .missing_needed_cols(needed[!needed %in% colnames(m)])
            )
            return(NULL)
        }
        tmp <- m %>%
            dplyr::select(dplyr::all_of(c(
                mandatory_IS_vars(),
                annotation_IS_vars()
            ))) %>%
            dplyr::distinct() %>%
            dplyr::group_by(dplyr::across(
                dplyr::all_of(mandatory_IS_vars())
            )) %>%
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
#' All options correspond to `lubridate` functions, see more in the dedicated
#' package documentation.
#'
#' @family Import functions helpers
#' @return A character vector
#' @export
#' @seealso \code{\link{import_association_file}},
#' \code{\link{import_parallel_Vispa2Matrices_auto}}
#'
#' @examples
#' date_formats()
date_formats <- function() {
    c(
        "ymd", "ydm", "mdy", "myd", "dmy", "dym", "yq", "ym", "my",
        "ymd_hms", "ymd_hm", "ymd_h", "dmy_hms", "dmy_hm",
        "dmy_h", "mdy_hms", "mdy_hm", "mdy_h", "ydm_hms", "ydm_hm", "ydm_h"
    )
}


#' Default regex prefixes for Vispa2 stats files.
#'
#' Note that each element is a regular expression.
#'
#' @family Import functions helpers
#' @return A character vector of regexes
#' @export
#'
#' @examples
#' default_iss_file_prefixes()
default_iss_file_prefixes <- function() {
    c("stats\\.sequence.", "stats\\.matrix.")
}


#' Default transformations to apply to association file columns.
#'
#' @description
#' A list of default transformations to apply to the association file columns
#' after importing it via `import_association_file()`
#'
#' @param convert_tp The value of the argument `convert_tp` in the call
#' to `import_association_file()`
#'
#' @family Import functions helpers
#' @return A named list of lambdas
#' @export
#'
#' @examples
#' default_af_transform(TRUE)
default_af_transform <- function(convert_tp) {
    if (convert_tp) {
        return(list(
            TimepointMonths = ~ stringr::str_pad(
                as.character(.x),
                pad = "0", side = "left", width = 2
            ),
            TimepointYears = ~ stringr::str_pad(as.character(.x),
                pad = "0", side = "left", width = 2
            )
        ))
    }
    return(NULL)
}
