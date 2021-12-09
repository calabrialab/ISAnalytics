#------------------------------------------------------------------------------#
# Collision removal functions
#------------------------------------------------------------------------------#

#' Identifies and removes collisions.
#'
#' \lifecycle{stable}
#' A collision is an integration (aka a unique combination of
#' `chr`, `integration_locus` and `strand`) which is observed in more than one
#' independent sample (a unique pair of `ProjectID` and `SubjectID`).
#' The function tries to decide to which subject an integration
#' should be assigned to and, if no
#' decision can be taken, the integration is completely removed from the data
#' frame.
#' For more details refer to the vignette "Collision removal functionality":
#' \code{vignette("collision_removal", package = "ISAnalytics")}
#'
#' @param x Either a multi-quantification matrix or a
#' named list of matrices (names must be quantification types)
#' @param association_file The association file imported via
#' `import_association_file()`
#' @param date_col The date column that should be considered.
#' Must be one value in `date_columns_coll()`
#' @param reads_ratio A single numeric value that represents the ratio that has
#' to be considered when deciding between `seqCount` value.
#' @param quant_cols A named character vector where names are
#' quantification types and
#' values are the names of the corresponding columns. The quantification
#' `seqCount` MUST be included in the vector.
#' @param report_path The path where the report file should be saved.
#' Can be a folder, a file or NULL if no report should be produced.
#' Defaults to `{user_home}/ISAnalytics_reports`.
#' @param max_workers Maximum number of parallel workers to distribute the
#' workload. If `NULL` (default) produces the maximum amount of workers allowed,
#' a numeric value is requested otherwise. WARNING: a higher number of workers
#' speeds up computation at the cost of memory consumption! Tune this parameter
#' accordingly.
#'
#' @family Collision removal
#' @importFrom magrittr `%>%`
#' @importFrom rlang inform abort exec .data
#' @importFrom purrr map2
#' @importFrom dplyr bind_rows select all_of group_by summarise across n
#' @importFrom dplyr distinct
#'
#' @seealso \code{\link{date_columns_coll}}
#'
#' @return Either a multi-quantification matrix or a list of data frames
#' @export
#'
#' @examples
#' data("integration_matrices", package = "ISAnalytics")
#' data("association_file", package = "ISAnalytics")
#' no_coll <- remove_collisions(
#'     x = integration_matrices,
#'     association_file = association_file,
#'     report_path = NULL
#' )
#' head(no_coll)
remove_collisions <- function(x,
    association_file,
    date_col = "SequencingDate",
    reads_ratio = 10,
    quant_cols = c(
        seqCount = "seqCount",
        fragmentEstimate = "fragmentEstimate"
    ),
    report_path = default_report_path(),
    max_workers = NULL) {

    # Check basic parameter correctness
    stopifnot(is.list(x) & !is.null(names(x)))
    stopifnot(is.character(quant_cols) && all(!is.na(names(quant_cols))))
    if (!all(names(quant_cols) %in% quantification_types())) {
        rlang::abort(
            .quantifications_names_err(
                quant_cols[!names(quant_cols) %in% quantification_types()]
            )
        )
    }
    stopifnot(is.null(max_workers[1]) || is.numeric(max_workers[1]))
    max_workers <- max_workers[1]
    mode <- NULL
    if (!is.data.frame(x)) {
        if (!all(names(x) %in% quantification_types())) {
            rlang::abort(.quantifications_names_err(
                names(x)[!names(x) %in% quantification_types()]
            ))
        }
        ## remove_collisions requires seqCount matrix, check if the list
        ## contains one
        if ((!"seqCount" %in% names(x)) ||
            nrow(x$seqCount) == 0) {
            rlang::abort(.seqCount_df_err())
        }
        all_correct <- purrr::map2(x, names(x), function(m, quant) {
            mand_cols <- .check_mandatory_vars(m)
            cAmp_col <- .check_complAmpID(m)
            if (mand_cols & cAmp_col) {
                return(NA_character_)
            } else {
                msgs <- c()
                if (!mand_cols) {
                    msgs <- .missing_mand_vars()
                }
                if (!cAmp_col) {
                    msgs <- c(msgs, .missing_cAmp_sub_msg())
                }
                msgs <- paste0(quant, " - ", paste0(msgs, collapse = ";\n"))
                return(msgs)
            }
        })
        if (!all(is.na(all_correct))) {
            message <- unlist(all_correct[!is.na(all_correct)])
            rlang::abort(c("Matrices miss required info, aborting", message))
        }
        ## Transform the list in a multi-quant matrix
        mode <- "LIST"
        quant_cols_lst <- as.list(quant_cols)
        args <- append(list(x = x), quant_cols_lst)
        x <- rlang::exec(comparison_matrix, !!!args)
    } else {
        if (.check_mandatory_vars(x) == FALSE) {
            rlang::abort(.missing_mand_vars())
        }
        if (.check_complAmpID(x) == FALSE) {
            rlang::abort(.missing_cAmp_sub_msg())
        }
        if (!all(quant_cols %in% colnames(x))) {
            rlang::abort(.missing_user_cols_error(
                quant_cols[!quant_cols %in% colnames(x)]
            ))
        }
        if (!"seqCount" %in% names(quant_cols)) {
            rlang::abort(.seqCount_col_err())
        }
        mode <- "MULTI"
    }

    stopifnot(is.data.frame(association_file))
    stopifnot(is.character(date_col))
    date_col <- date_col[1]
    stopifnot(date_col %in% date_columns_coll())
    stopifnot(is.integer(reads_ratio) || is.numeric(reads_ratio))
    reads_ratio <- reads_ratio[1]
    # Check association file correctness
    needed_af_cols <- c(
        "SubjectID", "ProjectID", "PoolID",
        "ReplicateNumber", date_col,
        "CompleteAmplificationID"
    )
    if (any(!needed_af_cols %in% colnames(association_file))) {
        missing_af_cols <- needed_af_cols[!needed_af_cols %in%
            colnames(association_file)]
        rlang::abort(.missing_af_needed_cols(missing_af_cols))
    }
    seq_count_col <- quant_cols["seqCount"]

    # Check date_col
    if (any(is.na(association_file[[date_col]]))) {
        rlang::abort(.na_in_date_col())
    }

    ## Check if association file contains all info relative to content the of
    ## the matrix
    verbose <- getOption("ISAnalytics.verbose")
    ### - Are all sample in the matrix present in the AF?
    missing_ind <- if (!all(x$CompleteAmplificationID %in%
        association_file$CompleteAmplificationID)) {
        which(!x$CompleteAmplificationID %in%
            association_file$CompleteAmplificationID)
    } else {
        NULL
    }
    pre_process <- if (!is.null(missing_ind)) {
        rlang::inform(.missing_af_samples_msg(
            length(unique(x$CompleteAmplificationID[missing_ind]))
        ))
        x[-missing_ind, ]
    } else {
        x
    }

    ## If after removing missing the matrix is empty, stop
    if (nrow(pre_process) == 0) {
        rlang::inform("Matrix is empty after removing missing samples")
        to_return <- if (mode == "LIST") {
            args <- append(list(x = x), as.list(quant_cols))
            rlang::exec(separate_quant_matrices, !!!args)
        } else {
            x
        }
        return(to_return)
    }

    ## Check if association file contains more info than matrix and
    ## keep only metadata that concerns projectIDs in the matrix
    add_samples <- .check_same_info(association_file, pre_process)
    association_file <- add_samples$reduced_af
    add_samples <- add_samples$miss
    if (nrow(add_samples) > 0) {
        if (verbose) {
            rlang::inform(.additional_ad_samples_msg())
        }
    } else {
        add_samples <- NULL
    }

    # Begin workflow
    ## - Join based on CompleteAmplificationID
    joined <- .join_matrix_af(pre_process, association_file, date_col)
    if (verbose) {
        rlang::inform("Identifying collisions...")
    }
    ## - Separate collisions from non-collisions
    splitted_df <- .identify_independent_samples(joined)

    if (nrow(splitted_df$collisions) == 0) {
        rlang::inform(.no_collisions_msg())
        to_return <- if (mode == "LIST") {
            args <- append(list(x = x), as.list(quant_cols))
            rlang::exec(separate_quant_matrices, !!!args)
        } else {
            x
        }
        return(to_return)
    }

    # Remove collisions
    if (verbose) {
        rlang::inform("Processing collisions...")
    }
    fixed_collisions <- .process_collisions(
        splitted_df$collisions,
        date_col, reads_ratio,
        seqCount_col = seq_count_col,
        max_workers = max_workers
    )
    reassigned <- fixed_collisions$reassigned
    removed <- fixed_collisions$removed
    fixed_collisions <- fixed_collisions$coll
    final_matr <- fixed_collisions %>%
        dplyr::bind_rows(splitted_df$non_collisions) %>%
        dplyr::select(dplyr::all_of(colnames(pre_process)))
    if (getOption("ISAnalytics.reports") == TRUE & !is.null(report_path)) {
        input_summary <- .summary_input(x, quant_cols)
        missing_smpl <- if (!is.null(missing_ind)) {
            x[missing_ind, ] %>%
                dplyr::group_by(.data$CompleteAmplificationID) %>%
                dplyr::summarise(
                    n_IS = dplyr::n(),
                    dplyr::across(
                        dplyr::all_of(quant_cols),
                        .fns = ~ sum(.x, na.rm = TRUE),
                        .names = "{.col}_tot"
                    )
                )
        } else {
            NULL
        }
        samples_info <- list(
            MATRIX = unique(x$CompleteAmplificationID),
            AF = unique(association_file$CompleteAmplificationID)
        )
        pre_summary <- .summary_input(pre_process, quant_cols)
        per_pool_stats_pre <- .per_pool_stats(joined, quant_cols)
        sharing_pre <- is_sharing(joined,
            group_key = c(
                "ProjectID",
                "SubjectID"
            ), n_comp = 2,
            is_count = FALSE,
            minimal = FALSE,
            include_self_comp = TRUE
        )
        coll_info <- list(
            coll_n = splitted_df$collisions %>%
                dplyr::distinct(
                    .data$chr,
                    .data$integration_locus,
                    .data$strand
                ) %>%
                nrow(),
            removed = removed,
            reassigned = reassigned
        )
        post_info <- .summary_input(final_matr, quant_cols)
        post_joined <- .join_matrix_af(final_matr,
            association_file,
            date_col = date_col
        )
        post_per_pool_stats <- .per_pool_stats(
            post_joined,
            quant_cols
        )
        sharing_post <- is_sharing(post_joined,
            group_key = c(
                "ProjectID",
                "SubjectID"
            ), n_comp = 2,
            is_count = FALSE,
            minimal = FALSE,
            include_self_comp = TRUE
        )
        summary_tbl <- .summary_table(
            before = joined, after = post_joined,
            seqCount_col = seq_count_col
        )
        withCallingHandlers(
            {
                .produce_report("collisions",
                    params = list(
                        input_info = input_summary,
                        missing_info = missing_smpl,
                        additional_info = add_samples,
                        samples_info = samples_info,
                        pre_info = pre_summary,
                        pre_stats = per_pool_stats_pre,
                        sharing_pre = sharing_pre,
                        coll_info = coll_info,
                        post_info = post_info,
                        post_stats = post_per_pool_stats,
                        summary_post = summary_tbl,
                        sharing_post = sharing_post
                    ), path = report_path
                )
            },
            error = function(cnd) {
                rest <- findRestart("report_fail")
                invokeRestart(rest, cnd)
            }
        )
    }

    ## If input was list, return the result in list form
    final <- if (mode == "LIST") {
        args <- append(list(x = final_matr), as.list(quant_cols))
        rlang::exec(separate_quant_matrices, !!!args)
    } else {
        final_matr
    }
    if (verbose) {
        rlang::inform("Finished!")
    }
    return(final)
}

#' Re-aligns matrices of other quantification types based on the processed
#' sequence count matrix.
#'
#' \lifecycle{stable}
#' This function should be used to keep data consistent among the same analysis:
#' if for some reason you removed the collisions by passing only the sequence
#' count matrix to `remove_collisions()`, you should call this
#' function afterwards, providing a list of other quantification matrices.
#' NOTE: if you provided a list of several quantification types to
#' `remove_collisions()` before, there is no need to call this function.
#'
#' @details For more details on how to use collision removal functionality:
#' \code{vignette("collision_removal", package = "ISAnalytics")}
#'
#' @param sc_matrix The sequence count matrix already processed for collisions
#' via `remove_collisions()`
#' @param other_matrices A named list of matrices to re-align. Names in the list
#' must be quantification types (\code{quantification_types()}) except
#' "seqCount".
#' @importFrom dplyr semi_join
#' @importFrom purrr map_lgl
#' @importFrom magrittr `%>%`
#' @family Collision removal
#' @seealso \code{\link{remove_collisions}}
#'
#' @return A named list with re-aligned matrices
#' @export
#'
#' @examples
#' data("integration_matrices", package = "ISAnalytics")
#' data("association_file", package = "ISAnalytics")
#' separated <- separate_quant_matrices(
#'     integration_matrices
#' )
#' no_coll <- remove_collisions(
#'     x = separated$seqCount,
#'     association_file = association_file,
#'     quant_cols = c(seqCount = "Value"),
#'     report_path = NULL
#' )
#' realigned <- realign_after_collisions(
#'     sc_matrix = no_coll,
#'     other_matrices = list(fragmentEstimate = separated$fragmentEstimate)
#' )
#' realigned
realign_after_collisions <- function(sc_matrix, other_matrices) {
    stopifnot(is.list(other_matrices) & !is.null(names(other_matrices)))
    stopifnot(all(names(other_matrices) %in% quantification_types()))
    all_ISm <- purrr::map_lgl(other_matrices, .check_mandatory_vars)
    if (!all(all_ISm)) {
        rlang::abort(.non_ISM_error())
    }
    all_campid <- purrr::map_lgl(other_matrices, .check_complAmpID)
    if (!all(all_campid)) {
        rlang::abort(.missing_complAmpID_error())
    }
    realigned <- purrr::map(other_matrices, function(x) {
        x %>% dplyr::semi_join(sc_matrix,
            by = c(
                mandatory_IS_vars(),
                "CompleteAmplificationID"
            )
        )
    })
    realigned
}

#' Possible choices for `date_col` parameter.
#'
#' @return A character vector of column names
#' @family Collision removal helpers
#' @seealso \code{\link{remove_collisions}}
#' @export
#'
#' @examples
#' dates <- date_columns_coll()
date_columns_coll <- function() {
    c("SequencingDate", "FusionPrimerPCRDate")
}
