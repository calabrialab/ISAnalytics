#------------------------------------------------------------------------------#
# Raw reads filtering functions
#------------------------------------------------------------------------------#

#---- Filter outliers ---------------------------------------------------------#
#' Filter out outliers in metadata, identified by the chosen outlier test.
#'
#' \lifecycle{experimental}
#' Filter out outliers in metadata.
#'
#' @param metadata The metadata data frame
#' @param outlier_test A string representing a function name. The name
#' must be one of the available outlier tests, see
#' \link{available_outlier_tests}.
#' @param negate If TRUE will return only the metadata that was flagged to
#' be removed. If FALSE will return only the metadata that wasn't flagged
#' to be removed.
#' @param ... Additional named arguments passed to `outliers_test`
#'
#' @family Outliers filter
#' @importFrom rlang abort list2 exec as_function .data
#' @importFrom dplyr filter select
#' @return A data frame of metadata which has less or the same amount of rows
#'
#' @export
#' @examples
#' op <- options(ISAnalytics.widgets = FALSE)
#'
#' path_AF <- system.file("extdata", "ex_association_file.tsv",
#'     package = "ISAnalytics"
#' )
#' association_file <- import_association_file(path_AF, root = NULL,
#'     dates_format = "dmy"
#' )
#' filtered_af <- outlier_filter(association_file, key = "VCN")
#' options(op)
outlier_filter <- function(metadata,
                           outlier_test = "outliers_by_pool_fragments",
                           negate = FALSE,
                           ...) {
  stopifnot(is.data.frame(metadata))
  stopifnot(is.character(outlier_test))
  stopifnot(is.logical(negate))
  if (length(outlier_test) > 1) {
    outlier_test <- outlier_test[1]
  }
  if (!outlier_test %in% available_outlier_tests()) {
    rlang::abort(c("Unknown outlier test",
                   i = paste("To know what outliers tests are available",
                             "type `available_outlier_tests()`")
                   ))
  }
  dots <- rlang::list2(...)
  flagged_meta <- rlang::exec(rlang::as_function(outlier_test),
                metadata = metadata,
                !!!dots)
  filtered_meta <- if (negate) {
    flagged_meta %>%
      dplyr::filter(.data$to_remove == TRUE) %>%
      dplyr::select(-.data$to_remove)
  } else {
    flagged_meta %>%
      dplyr::filter(.data$to_remove == FALSE) %>%
      dplyr::select(-.data$to_remove)
  }
  filtered_meta
}


#---- Outlier tests -----------------------------------------------------------#
### Outlier tests perform some calculation on data and flag outliers according
### to parameters

#' Identify and flag outliers based on pool fragments.
#'
#' \lifecycle{experimental}
#' Identify and flag outliers
#'
#' @details
#' This particular test calculates for each column in the key
#'
#' * The zscore of the values
#' * The tstudent of the values
#' * The the distribution of the tstudent values
#'
#' Optionally the test can be performed for each pool and a normality test
#' can be run prior the actual calculations.
#' Samples are flagged if this condition is respected:
#'
#' * tdist < outlier_p_value_threshold & zscore < 0
#'
#' If the key contains more than one column an additional flag logic can be
#' specified for combining the results.
#' Example:
#' let's suppose the key contains the names of two columns, X and Y
#' `key = c("X", "Y")`
#' if we specify the the argument `flag_logic = "AND"` then the reads will
#' be flagged based on this global condition:
#' (tdist_X < outlier_p_value_threshold & zscore_X < 0) AND
#' (tdist_Y < outlier_p_value_threshold & zscore_Y < 0)
#'
#' The user can specify one or more logical operators that will be applied
#' in sequence.
#'
#' @param metadata The metadata data frame
#' @param key A character vector of numeric column names
#' @param outlier_p_value_threshold The p value threshold for a read to be
#' considered an outlier
#' @param normality_test Perform normality test? Normality is assessed for
#' each column in the key using Shapiro-Wilk test and if the values do not
#' follow a normal distribution, other calculations are skipped
#' @param normality_p_value_threshold Normality threshold
#' @param transform_log2 Perform a log2 trasformation on values prior the
#' actual calculations?
#' @param per_pool_test Perform the test for each pool?
#' @param pool_col A character vector of the names of the columns that
#' uniquely identify a pool
#' @param min_samples_per_pool The minimum number of samples that a pool
#' needs to contain in order to be processed - relevant only if
#' `per_pool_test = TRUE`
#' @param flag_logic A character vector of logic operators to obtain a
#' global flag formula - only relevant if the key is longer than one.
#' All operators must be chosen between:
#' `r paste0(flag_logics(), collapse = ", ")`
#' @param keep_calc_cols Keep the calculation columns in the output data frame?
#' @param save_widget_path Either null or a string containing the path
#' on disk where the report should be saved
#'
#' @return A data frame of metadata with the column `to_remove`
#'
#' @family Outlier tests
#' @export
#' @examples
#' op <- options(ISAnalytics.widgets = FALSE)
#'
#' path_AF <- system.file("extdata", "ex_association_file.tsv",
#'     package = "ISAnalytics"
#' )
#' association_file <- import_association_file(path_AF, root = NULL,
#'     dates_format = "dmy"
#' )
#' filtered_af <- outliers_by_pool_fragments(association_file, key = "VCN")
#' options(op)
outliers_by_pool_fragments <- function(metadata,
    key = "BARCODE_MUX",
    outlier_p_value_threshold = 0.05,
    normality_test = FALSE,
    normality_p_value_threshold = 0.05,
    transform_log2 = TRUE,
    per_pool_test = TRUE,
    pool_col = "PoolID",
    min_samples_per_pool = 5,
    flag_logic = "AND",
    keep_calc_cols = TRUE,
    save_widget_path = NULL) {
    ## Check
    stopifnot(is.data.frame(metadata))
    stopifnot(is.character(key))
    if (!all(key %in% colnames(metadata))) {
        rlang::abort(.missing_user_cols_error(key[!key %in%
            colnames(metadata)]),
        class = "missing_cols_key"
        )
    }
    stopifnot(is.numeric(outlier_p_value_threshold))
    stopifnot(is.logical(normality_test) && length(normality_test) == 1)
    if (normality_test) {
        stopifnot(is.numeric(normality_p_value_threshold) &&
            length(normality_p_value_threshold) == 1)
    }
    stopifnot(is.logical(transform_log2) && length(transform_log2) == 1)
    stopifnot(is.logical(per_pool_test) && length(per_pool_test) == 1)
    if (per_pool_test) {
        stopifnot(is.character(pool_col))
        stopifnot(is.numeric(min_samples_per_pool) &&
            length(min_samples_per_pool) == 1)
        if (!all(pool_col %in% colnames(metadata))) {
            rlang::abort(.missing_user_cols_error(
                pool_col[!pool_col %in% colnames(metadata)]
            ),
            class = "missing_cols_pool"
            )
        }
    }
    if (!"CompleteAmplificationID" %in% colnames(metadata)) {
        rlang::abort(.missing_complAmpID_error())
    }
    if (length(key) > 1) {
        ## Verify logic
        stopifnot(is.character(flag_logic))
        if (length(flag_logic) > length(key) - 1) {
            flag_logic <- flag_logic[seq_len(length(key) - 1)]
            if (getOption("ISAnalytics.verbose") == TRUE) {
                rlang::inform(c("'flag_logic' has more elements than expected",
                    i = paste(
                        "The vector will be trimmed to consider",
                        "the first", length(key) - 1, "elements",
                        "only."
                    )
                ), class = "flag_logic_long")
            }
        } else if (length(flag_logic) != length(key) - 1 &
            length(flag_logic) != 1) {
            flag_logic <- flag_logic[1]
            if (getOption("ISAnalytics.verbose") == TRUE) {
                rlang::inform(c(paste(
                    "'flag_logic' has an incorrect",
                    "amount of elements"
                ),
                paste(
                    "You should provide a minimum of 1 and a",
                    "maximum of",
                    length(key) - 1,
                    "logical operators"
                ),
                i = "Only the first parameter will be considered"
                ),
                class = "flag_logic_short"
                )
            }
        }
        flag_logic <- toupper(flag_logic)
        if (!all(flag_logic %in% flag_logics())) {
            rlang::abort(c(paste(
                "Unknown or unsupported logical operators:",
                paste0(flag_logic[!flag_logic %in% flag_logics()],
                    collapse = ", "
                )
            )))
        }
    }
    ## Cleaning
    if (getOption("ISAnalytics.verbose") == TRUE) {
        rlang::inform("Removing NAs from data...")
    }
    metadata_filtered <- metadata %>%
        dplyr::filter(dplyr::across(
            .cols = dplyr::all_of(key),
            .fns = ~ !is.na(.x)
        ))
    if (nrow(metadata_filtered) == 0) {
        rlang::inform(c("Metadata empty",
            paste(
                "Metadata has 0 rows after filtering,",
                "nothing to do"
            ),
            i = "Returning input data frame"
        ))
        return(metadata)
    }

    removed_nas <- metadata %>%
        dplyr::anti_join(metadata_filtered, by = "CompleteAmplificationID")

    ## Calculation
    calc_res <- .pool_frag_calc(
        meta = metadata_filtered,
        key = key,
        by_pool = per_pool_test,
        normality_test = normality_test,
        normality_threshold = normality_p_value_threshold,
        pool_col = pool_col,
        min_samples_per_pool = min_samples_per_pool,
        log2 = transform_log2
    )
    ## Extract dfs for report
    non_proc_samples <- if ("non_proc_samples" %in% names(calc_res)) {
        calc_res$non_proc_samples %>%
            dplyr::select(
                dplyr::all_of(pool_col),
                .data$CompleteAmplificationID
            ) %>%
            dplyr::distinct()
    } else {
        NULL
    }
    removed_zeros <- calc_res$removed_zeros
    calc_res <- calc_res$metadata

    ## Obtain final df
    calc_res <- calc_res %>% dplyr::bind_rows(
        removed_nas %>% dplyr::mutate(processed = FALSE)
    )

    ## Flag outliers
    for (k in key) {
        calc_res <- calc_res %>%
            dplyr::mutate("flag_{k}" := .flag_cond_outliers_pool_frag(
                proc = .data$processed,
                tdist = .data[[paste0("tdist_", k)]],
                zscore = .data[[paste0("zscore_", k)]],
                outlier_threshold = outlier_p_value_threshold
            ))
    }
    to_rem <- if (length(key) == 1) {
        dplyr::pull(calc_res, name = paste0("flag_", key))
    } else {
        purrr::pmap_lgl(calc_res[seq(
            from = length(calc_res) - length(key) + 1,
            to = length(calc_res)
        )], .apply_flag_logic, logic = flag_logic)
    }
    calc_res <- calc_res %>%
        dplyr::mutate(to_remove = to_rem)
    ## Produce report
    if (getOption("ISAnalytics.widgets") == TRUE) {
        rlang::inform("Producing report...")
        withCallingHandlers(
            {
                withRestarts(
                    {
                        flagged <- calc_res %>%
                            dplyr::select(
                                .data$processed,
                                .data$to_remove,
                                .data$ProjectID,
                                dplyr::all_of(c(pool_col)),
                                .data$CompleteAmplificationID,
                                dplyr::contains(key)
                            )
                        widget <- .outliers_report_widg(
                            by_pool = per_pool_test,
                            pool_col = pool_col,
                            norm_test = normality_test,
                            key = key,
                            flag_logic = flag_logic,
                            outlier_thresh = outlier_p_value_threshold,
                            log2_req = transform_log2,
                            removed_nas = removed_nas,
                            removed_zeros = removed_zeros,
                            flag_df = flagged,
                            non_proc_pools = non_proc_samples
                        )
                        print(widget)
                        # Export widget if requested
                        if (!is.null(save_widget_path)) {
                            if (getOption("ISAnalytics.verbose") == TRUE) {
                                rlang::inform("Saving report to file...")
                            }
                            .export_widget_file(
                                widget,
                                save_widget_path,
                                "raw_reads_report.html"
                            )
                        }
                    },
                    print_err = function() {
                        rlang::inform(.widgets_error())
                    }
                )
            },
            error = function(cnd) {
                rlang::inform(conditionMessage(cnd))
                invokeRestart("print_err")
            }
        )
    }
    if (!keep_calc_cols) {
      calc_res <- calc_res %>%
        dplyr::select(dplyr::all_of(colnames(metadata)), .data$to_remove)
    }
    calc_res
}

#' A character vector containing all the names of the currently supported
#' outliers tests that can be called in the function \link{outlier_filter}.
#'
#' @return A character vector
#'
#' @family Outlier tests
#' @export
#' @examples
#' available_outlier_tests()
available_outlier_tests <- function() {
  c("outliers_by_pool_fragments")
}
