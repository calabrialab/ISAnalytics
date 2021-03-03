#------------------------------------------------------------------------------#
# Re-calibration functions
#------------------------------------------------------------------------------#
#' Scans input matrix to find and merge near integration sites.
#'
#' \lifecycle{experimental}
#' This function scans the input integration matrix to detect eventual
#' integration sites that are too "near" to each other and merges them
#' into single integration sites adjusting their values if needed.
#'
#' @details The whole matrix is scanned with a sliding window mechanism:
#' for each row in the integration matrix an interval is calculated
#' based on the `threshold` value, then a "look ahead" operation is
#' performed to detect subsequent rows which integration locuses fall
#' in the interval. If CompleteAmplificationIDs of the near integrations
#' are different only the locus value (and optionally
#' GeneName and GeneStrand if the matrix is annotated) is modified,
#' otherwise rows with the same id are aggregated and values are summed.
#' If one of the map parameters is set to true the function will also
#' produce a re-calibration map: this data frame contains the reference
#' of pre-recalibration values for chr, strand and integration locus and
#' the value to which that integration was changed to after.
#'
#' @note We do recommend to use this function in combination with
#' \link{comparison_matrix} to automatically perform re-calibration on
#' all quantification matrices.
#'
#' @param x A single integration matrix, either with a single "Value"
#' column or multiple value columns corresponding to different
#' quantification types (obtained via \link{comparison_matrix})
#' @param threshold A single integer that represents an absolute
#' number of bases for which two integrations are considered distinct
#' @param keep_criteria While scanning, which integration should be kept?
#' The 2 possible choices for this parameter are:
#' * "max_value": keep the integration site which has the highest value
#' (and collapse other values on that integration).
#' * "keep_first": keeps the first integration
#' @param strand_specific Should strand be considered? If yes,
#' for example these two integration sites
#' c(chr = "1", strand = "+", integration_locus = 14568) and
#' c(chr = "1", strand = "-", integration_locus = 14568) are considered
#' different and not grouped together.
#' @param max_value_column The column that has to be considered for
#' searching the maximum value
#' @param map_as_widget Produce recalibration map as an HTML widget?
#' @param map_as_file Produce recalibration map as a .tsv file?
#' @param file_path String representing the path were the file will be
#' saved. By default the function produces a folder in the current
#' working directory and generates file names with time stamps.
#' @param export_widget_path A path on disk to save produced widgets or NULL
#' if the user doesn't wish to save the html file
#' @importFrom magrittr `%>%`
#' @importFrom tibble is_tibble tibble_row
#' @importFrom dplyr group_by group_split bind_rows
#' @importFrom purrr reduce map_lgl map_dfr
#' @importFrom rlang .data
#' @import BiocParallel
#'
#' @family Recalibration functions
#'
#' @return An integration matrix with same or less number of rows
#' @export
#'
#' @examples
#' path <- system.file("extdata", "ex_annotated_ISMatrix.tsv.xz",
#'     package = "ISAnalytics"
#' )
#' matrix <- import_single_Vispa2Matrix(path)
#' near <- compute_near_integrations(matrix,
#'     map_as_widget = FALSE,
#'     map_as_file = FALSE
#' )
compute_near_integrations <- function(x,
    threshold = 4,
    keep_criteria = "max_value",
    strand_specific = TRUE,
    max_value_column = "seqCount",
    map_as_widget = TRUE,
    map_as_file = TRUE,
    file_path = ".",
    export_widget_path = NULL) {
    # Check parameters
    stopifnot(tibble::is_tibble(x))
    ## Check tibble is an integration matrix
    if (.check_mandatory_vars(x) == FALSE) {
        stop(.non_ISM_error())
    }
    ## Check it contains CompleteAmplificationID column
    if (.check_complAmpID(x) == FALSE) {
        stop(.missing_complAmpID_error())
    }
    ## Check if it contains the "Value" column. If not find all numeric columns
    ## that are not default columns
    num_cols <- if (.check_value_col(x) == FALSE) {
        found <- .find_exp_cols(x)
        if (purrr::is_empty(found)) {
            stop(.missing_num_cols_error())
        }
        found
    } else {
        "Value"
    }
    stopifnot(is.numeric(threshold) | is.integer(threshold))
    stopifnot(length(threshold) == 1)
    criteria <- rlang::arg_match(keep_criteria, values = c(
        "max_value",
        "keep_first"
    ))
    stopifnot(is.logical(strand_specific) &
        length(strand_specific) == 1)
    stopifnot(is.character(max_value_column) & length(max_value_column) == 1)
    stopifnot(is.logical(map_as_widget) &
        length(map_as_widget) == 1)
    stopifnot(is.logical(map_as_file) &
        length(map_as_file) == 1)
    if (map_as_file == TRUE) {
        stopifnot(is.character(file_path) & length(file_path) == 1)
    }
    ## Notify user the name of columns which will be used
    if (getOption("ISAnalytics.verbose") == TRUE) {
        message(paste(c("Using the following numeric columns: ", num_cols),
            collapse = " "
        ))
    }
    ## Check that max_value_col parameter is set correctly if criteria is
    ## max_value
    if (criteria == "max_value") {
        if (!max_value_column %in% num_cols) {
            if (all(num_cols == "Value")) {
                warning(.using_val_col_warning(max_value_column),
                    immediate. = TRUE
                )
                max_value_column <- "Value"
            } else {
                stop(.max_val_stop_error(max_value_column))
            }
        }
    }
    ## Is x annotated?
    annotated <- .is_annotated(x)
    # Split data for parallel execution
    split <- if (strand_specific == TRUE) {
        x %>%
            dplyr::group_by(.data$chr, .data$strand) %>%
            dplyr::group_split()
    } else {
        x %>%
            dplyr::group_by(.data$chr) %>%
            dplyr::group_split()
    }
    ## Select only groups with 2 or more rows
    split_to_process <- split[purrr::map_lgl(split, function(x) {
        nrow(x) > 1
    })]
    # # Set up parallel workers
    if (.Platform$OS.type == "windows") {
        p <- BiocParallel::SnowParam(
            stop.on.error = FALSE,
            progressbar = TRUE,
            tasks = 20
        )
    } else {
        p <- BiocParallel::MulticoreParam(
            stop.on.error = FALSE, progressbar = TRUE, tasks = 20
        )
    }
    FUN <- function(x) {
        .sliding_window(x,
            threshold = threshold, keep_criteria = criteria,
            num_cols = num_cols, annotated = annotated,
            max_val_col = max_value_column,
            produce_map = any(map_as_file, map_as_widget)
        )
    }
    ## Obtain result: list of lists
    suppressWarnings({
        result <- BiocParallel::bptry(
            BiocParallel::bplapply(split_to_process,
                FUN = FUN,
                BPPARAM = p
            )
        )
    })
    BiocParallel::bpstop(p)
    ## Obtain single list
    result <- purrr::reduce(result, function(group1, group2) {
        recalibr_m <- dplyr::bind_rows(
            group1$recalibrated_matrix,
            group2$recalibrated_matrix
        )
        maps <- dplyr::bind_rows(
            group1$map,
            group2$map
        )
        list(recalibrated_matrix = recalibr_m, map = maps)
    })
    ## Add all rows that were not part of recalibration
    split_fine <- split[purrr::map_lgl(split, function(x) {
        !nrow(x) > 1
    })]
    result$recalibrated_matrix <- result$recalibrated_matrix %>%
        dplyr::bind_rows(split_fine)
    rows_to_add <- purrr::map_dfr(seq_along(nrow(split_fine)), function(ind) {
        row <- tibble::tibble_row(
            chr_before = split_fine$chr[ind],
            integration_locus_before =
                split_fine$integration_locus[ind],
            strand_before = split_fine$strand[ind],
            chr_after = split_fine$chr[ind],
            integration_locus_after =
                split_fine$integration_locus[ind],
            strand_after = split_fine$strand[ind]
        )
    })
    result$map <- result$map %>% dplyr::bind_rows(rows_to_add)

    if (map_as_widget == TRUE & getOption("ISAnalytics.widgets") == TRUE) {
        ## Produce widget
        withCallingHandlers(
            {
                withRestarts(
                    {
                        widget <- .recalibr_map_widget(result$map)
                        print(widget)
                        if (!is.null(export_widget_path)) {
                            .export_widget_file(
                                widget,
                                export_widget_path,
                                "recalibration_map.html"
                            )
                        }
                    },
                    print_err = function() {
                        message(.widgets_error())
                    }
                )
            },
            error = function(cnd) {
                message(conditionMessage(cnd))
                invokeRestart("print_err")
            }
        )
    }

    if (map_as_file == TRUE) {
        ### Manage file
        withCallingHandlers(
            {
                .write_recalibr_map(result$map, file_path)
            },
            error = function(e) {
                r <- findRestart("skip_write")
                if (is.null(r)) {
                    return()
                }
                invokeRestart(r)
            }
        )
    }
    return(result$recalibrated_matrix)
}
