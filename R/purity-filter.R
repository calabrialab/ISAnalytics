#' Filter integration sites based on purity.
#'
#' @description
#' `r lifecycle::badge("stable")`
#' Filter that targets possible contamination between cell lines based on
#' a numeric quantification (likely abundance or sequence count).
#'
#' @details
#' ## Setting input arguments
#'
#' The input matrix can be re-aggregated with the provided `group_key`
#' argument. This key contains the names of the columns to group on
#' (besides the columns holding genomic coordinates of the integration
#' sites) and must be contained in at least one of `x` or `lineages`
#' data frames. If the key is not found only in `x`, then a join operation
#' with the `lineages` data frame is performed on the common column(s)
#' `join_on`.
#'
#' ## Group selection
#' It is possible for the user to specify on which groups the logic of the
#' filter should be applied to. For example: if we have
#' `group_key = c("HematoLineage")` and we set
#' `selected_groups = c("CD34", "Myeloid","Lymphoid")`
#' it means that a single integration will be evaluated for the filter only
#' for groups that have the values of "CD34", "Myeloid" and "Lymphoid" in
#' the "HematoLineage" column.
#' If the same integration is present in other groups it is
#' kept as it is. `selected_groups` can be set to `NULL` if we want
#' the logic to apply to every group present in the data frame,
#' it can be set as a simple character vector as the example above if
#' the group key has length 1 (and there is no need to filter on time point).
#' If the group key is longer than 1 then the filter is applied only on the
#' first element of the key.
#'
#' If a more refined selection on groups is needed, a data frame can
#' be provided instead:
#'
#' ```
#' group_key = c("CellMarker", "Tissue")
#' selected_groups = tibble::tribble(
#' ~ CellMarker, ~ Tissue,
#' "CD34", "BM",
#' "CD14", "BM",
#' "CD14", "PB"
#' )
#' ```
#'
#' Columns in the data frame should be the same as group key (plus,
#' eventually, the time point column). In this example only those groups
#' identified by the rows in the provided data frame are processed.
#'
#' @family Data cleaning and pre-processing
#'
#' @param x An aggregated integration matrix, obtained via
#' `aggregate_values_by_key()`
#' @param lineages A data frame containing cell lineages information
#' @param aggregation_key The key used for aggregating `x`
#' @param group_key A character vector of column names for re-aggregation.
#' Column names must be either in `x` or in `lineages`. See details.
#' @param selected_groups Either NULL, a character vector or a
#' data frame for group selection. See details.
#' @param join_on Common columns to perform a join operation on
#' @param min_value A minimum value to filter the input matrix. Integrations
#' with a value strictly lower than `min_value` are excluded (dropped) from
#' the output.
#' @param impurity_threshold The ratio threshold for impurity in groups
#' @param by_timepoint Should filtering be applied on each time point? If
#' `FALSE`, all time points are merged together
#' @param timepoint_column Column in `x` containing the time point
#' @param value_column Column in `x` containing the numeric
#' quantification of interest
#'
#' @return A data frame
#' @export
#'
#' @examples
#' data("integration_matrices", package = "ISAnalytics")
#' data("association_file", package = "ISAnalytics")
#' aggreg <- aggregate_values_by_key(
#'     x = integration_matrices,
#'     association_file = association_file,
#'     value_cols = c("seqCount", "fragmentEstimate")
#' )
#' filtered_by_purity <- purity_filter(
#'     x = aggreg,
#'     value_column = "seqCount_sum"
#' )
#' head(filtered_by_purity)
purity_filter <- function(x,
    lineages = blood_lineages_default(),
    aggregation_key = c(
        "SubjectID", "CellMarker",
        "Tissue", "TimePoint"
    ),
    group_key = c("CellMarker", "Tissue"),
    selected_groups = NULL,
    join_on = "CellMarker",
    min_value = 3,
    impurity_threshold = 10,
    by_timepoint = TRUE,
    timepoint_column = "TimePoint",
    value_column = "seqCount_sum") {
    ## Checks
    #### - Base
    stopifnot(is.data.frame(x))
    stopifnot(is.character(aggregation_key))
    stopifnot(is.character(group_key))
    stopifnot(is.logical(by_timepoint))
    stopifnot(is.numeric(min_value) || is.integer(min_value))
    stopifnot(is.numeric(impurity_threshold) || is.integer(impurity_threshold))
    stopifnot(is.character(value_column))
    stopifnot(is.null(selected_groups) || is.character(selected_groups) ||
        is.data.frame(selected_groups))
    #### - Keys
    if (!all(aggregation_key %in% colnames(x))) {
        rlang::abort(.missing_user_cols_error(
            aggregation_key[!aggregation_key %in% colnames(x)]
        ))
    }
    if (!value_column[1] %in% colnames(x)) {
        rlang::abort(.missing_user_cols_error(value_column))
    }
    to_join <- if (all(group_key %in% colnames(x))) {
        ### If the groups rely only on attributes not related to lineages
        ### and are contained in the input matrix
        FALSE
    } else {
        ### If lineages info is needed
        stopifnot(is.data.frame(lineages))
        if (!all(group_key %in% unique(c(colnames(x), colnames(lineages))))) {
            missing_cols <- group_key[!group_key %in% unique(c(
                colnames(x),
                colnames(lineages)
            ))]
            rlang::abort(.missing_user_cols_error(missing_cols))
        }
        if (!(all(join_on %in% colnames(x)) &
            all(join_on %in% colnames(lineages)))
        ) {
            missing_common <- c("Missing common column(s) to join on",
                i = paste(
                    "The column(s) provided in argument",
                    "`join_on` is missing from one or",
                    "both data frames, aborting"
                )
            )
            rlang::abort(missing_common)
        }
        TRUE
    }
    if (by_timepoint) {
        stopifnot(is.character(timepoint_column))
        timepoint_column <- timepoint_column[1]
        if (!timepoint_column %in% colnames(x)) {
            rlang::abort(.missing_user_cols_error(timepoint_column))
        }
        group_key <- union(group_key, timepoint_column)
    }
    ## Pre-processing
    #### - Join if needed
    if (to_join) {
        x <- x %>%
            dplyr::left_join(lineages, by = join_on)
    }
    #### - Group and sum
    is_vars <- if (.is_annotated(x)) {
        c(mandatory_IS_vars(), annotation_IS_vars())
    } else {
        mandatory_IS_vars()
    }
    grouped <- x %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c(is_vars, group_key)))) %>%
        dplyr::summarise(
            Value = sum(.data[[value_column]]),
            .groups = "drop"
        )
    #### - value filter
    filtered_value <- threshold_filter(
        x = grouped,
        threshold = min_value,
        cols_to_compare = "Value",
        comparators = ">="
    )
    #### - Separating IS 1: group filtering
    pre_filt <- list()
    if (is.null(selected_groups) || purrr::is_empty(selected_groups)) {
        pre_filt[["process"]] <- filtered_value
        pre_filt[["keep"]] <- filtered_value[0, ]
    } else if (is.character(selected_groups)) {
        pre_filt[["process"]] <- filtered_value %>%
            dplyr::filter(.data[[group_key[1]]] %in% selected_groups)
        pre_filt[["keep"]] <- filtered_value %>%
            dplyr::filter(!.data[[group_key[1]]] %in% selected_groups)
    } else {
        ok_cols <- colnames(selected_groups)[colnames(selected_groups) %in%
            group_key]
        selected_groups <- selected_groups %>%
            dplyr::select(dplyr::all_of(ok_cols)) %>%
            dplyr::distinct()
        if (ncol(selected_groups) == 0 ||
            nrow(selected_groups) == 0) {
            pre_filt[["process"]] <- filtered_value
            pre_filt[["keep"]] <- filtered_value[0, ]
        } else {
            pre_filt[["process"]] <- dplyr::inner_join(filtered_value,
                selected_groups,
                by = ok_cols
            )
            pre_filt[["keep"]] <- dplyr::anti_join(filtered_value,
                selected_groups,
                by = ok_cols
            )
        }
    }
    if (nrow(pre_filt$process) == 0) {
        if (getOption("ISAnalytics.verbose", TRUE)) {
            rlang::inform("No iss to process, done")
        }
        return(filtered_value)
    }
    #### - Separating IS 2: iss that are shared between groups are going to be
    #### processed, unique iss are kept as they are
    vars_to_group <- if (by_timepoint) {
        c(is_vars, timepoint_column)
    } else {
        is_vars
    }
    by_is <- pre_filt$process %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(vars_to_group))) %>%
        dplyr::summarise(n = n(), .groups = "drop")
    to_process <- by_is %>%
        dplyr::filter(.data$n > 1) %>%
        dplyr::select(-dplyr::all_of("n")) %>%
        dplyr::inner_join(pre_filt$process, by = vars_to_group)
    if (nrow(to_process) == 0) {
        ## If there are no shared iss there is nothing to process,
        ## return just the filtered matrix
        return(filtered_value)
    }
    to_keep <- by_is %>%
        dplyr::filter(.data$n == 1) %>%
        dplyr::select(-dplyr::all_of("n")) %>%
        dplyr::inner_join(pre_filt$process, by = vars_to_group)
    #### - Process groups
    .filter_by_purity <- function(group) {
        max_val <- max(group$Value)
        processed <- group %>%
            dplyr::mutate(remove = (max_val / .data$Value) >
                impurity_threshold) %>%
            dplyr::filter(remove == FALSE) %>%
            dplyr::select(-dplyr::all_of("remove"))
        processed
    }
    processed_iss <- to_process %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(vars_to_group))) %>%
        dplyr::group_modify(~ .filter_by_purity(.x)) %>%
        dplyr::ungroup()
    #### - Re-compose matrix
    final <- processed_iss %>%
        dplyr::bind_rows(to_keep) %>%
        dplyr::bind_rows(pre_filt$keep)
    final
}
