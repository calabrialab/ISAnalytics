#------------------------------------------------------------------------------#
# Re-calibration functions
#------------------------------------------------------------------------------#
#' Scans input matrix to find and merge near integration sites.
#'
#' \lifecycle{stable}
#' This function scans the input integration matrix to detect eventual
#' integration sites that are too "near" to each other and merges them
#' into single integration sites adjusting their values if needed.
#'
#' @details The whole matrix is scanned with a sliding window mechanism:
#' for each row in the integration matrix an interval is calculated
#' based on the `threshold` value, then a "look ahead" operation is
#' performed to detect subsequent rows which integration locuses fall
#' in the interval. If `CompleteAmplificationID`s of the near integrations
#' are different only the locus value (and optionally
#' `GeneName` and `GeneStrand` if the matrix is annotated) is modified,
#' otherwise rows with the same id are aggregated and values are summed.
#' The function will also
#' produce a re-calibration map: this data frame contains the reference
#' of pre-recalibration values for `chr`, `strand` and `integration_locus` and
#' the value to which that integration was changed to.
#'
#' @note We do recommend to use this function in combination with
#' \link{comparison_matrix} to automatically perform re-calibration on
#' all quantification matrices.
#'
#' @param x An integration matrix
#' @param threshold A single integer that represents an absolute
#' number of bases for which two integrations are considered distinct.
#' If the threshold is set to 3 it means, provided fields `chr`
#' and `strand` are the same, integrations sites
#' which have at least 3 bases in between them are
#' considered distinct (e.g. (1, 14576, +) and (1, 14580, +) are
#' considered distinct)
#' @param keep_criteria While scanning, which integration should be kept?
#' The 2 possible choices for this parameter are:
#' * "max_value": keep the integration site which has the highest value
#' (and collapse other values on that integration).
#' * "keep_first": keeps the first integration
#' @param strand_specific Should strand be considered? If yes,
#' for example these two integration sites
#' `(chr = "1", strand = "+", integration_locus = 14568)` and
#' `(chr = "1", strand = "-", integration_locus = 14568)` are considered
#' different and not grouped together.
#' @param value_columns Character vector, contains the names of the numeric
#' experimental columns
#' @param max_value_column The column that has to be considered for
#' searching the maximum value
#' @param map_as_file Produce recalibration map as a .tsv file?
#' @param file_path String representing the path were the file will be
#' saved. Can be either a folder or a file. Relevant only if `map_as_file` is
#' `TRUE`.
#' @importFrom magrittr `%>%`
#' @importFrom data.table rbindlist setDT setnames
#' @importFrom dplyr group_by group_split
#' @importFrom purrr reduce map_lgl map_dfr
#' @importFrom rlang .data abort arg_match
#' @importFrom BiocParallel SnowParam MulticoreParam bptry bplapply bpstop
#'
#' @family Recalibration functions
#'
#' @return An integration matrix with same or less number of rows
#' @export
#'
#' @examples
#' data("integration_matrices", package = "ISAnalytics")
#' rec <- compute_near_integrations(
#'     x = integration_matrices, map_as_file = FALSE
#' )
#' head(rec)
compute_near_integrations <- function(x,
    threshold = 4,
    is_identity_tags = c("chromosome", "is_strand"),
    keep_criteria = c("max_value", "keep_first"),
    value_columns = c("seqCount", "fragmentEstimate"),
    max_value_column = "seqCount",
    sample_id_column = pcr_id_column(),
    additional_agg_lambda = list(.default = default_rec_agg_lambdas()),
    max_workers = 4,
    map_as_file = TRUE,
    file_path = default_report_path(),
    strand_specific = lifecycle::deprecated()) {
    ## --- Check parameters
    stopifnot(is.data.frame(x))
    stopifnot(is.numeric(threshold) || is.integer(threshold))
    threshold <- threshold[1]
    stopifnot(is.character(value_columns))
    stopifnot(is.character(max_value_column))
    num_cols <- unique(c(value_columns, max_value_column[1]))
    if (!all(num_cols %in% colnames(x))) {
      rlang::abort(.missing_needed_cols(
        num_cols[!num_cols %in% colnames(x)]
      ))
    }
    criteria <- rlang::arg_match(keep_criteria)
    stopifnot(is.character(sample_id_column))
    sample_id_column <- sample_id_column[1]
    stopifnot(is.logical(map_as_file))
    map_as_file <- map_as_file[1]
    if (map_as_file == TRUE) {
      stopifnot(is.character(file_path))
      file_path <- file_path[1]
    }
    if (lifecycle::is_present(strand_specific)) {
      lifecycle::deprecate_warn(
        when = "1.5.4",
        what = "compute_near_integrations(strand_specific)",
        with = "compute_near_integrations(is_identity_tags)",
        details = paste("See the documentation for details,",
                        "the argument will be likely dropped in the next",
                        "release cycle.")
      )
      stopifnot(is.logical(strand_specific))
      strand_specific <- strand_specific[1]
      if (!"is_strand" %in% is_identity_tags) {
        is_identity_tags <- c(is_identity_tags, "is_strand")
      }
    }
    ## --- Check tags
    stopifnot(is.null(is_identity_tags) || is.character(is_identity_tags))
    required_tags <- list("locus" = c("int", "numeric"))
    if (!is.null(is_identity_tags)) {
      is_identity_tags <- is_identity_tags[is_identity_tags != "locus"]
      id_tags_types <- purrr::map(is_identity_tags, ~ NULL) %>%
        purrr::set_names(is_identity_tags)
      required_tags <- append(required_tags, id_tags_types)
    }
    required_tag_cols <- .check_required_cols(
      required_tags = required_tags,
      vars_df = mandatory_IS_vars(TRUE),
      duplicate_politic = "error")
    if (!all(required_tag_cols$names %in% colnames(x))) {
      rlang::abort(.missing_needed_cols(required_tag_cols$names[
        !required_tag_cols$names %in% colnames(x)
      ]))
    }
    if (!sample_id_column %in% colnames(x)) {
      rlang::abort(.missing_needed_cols(sample_id_column))
    }
    ## Is x annotated?
    annotated <- .is_annotated(x)
    ## Any additional columns present?
    additional_cols <- colnames(x)[!colnames(x) %in% c(
      required_tag_cols$names, sample_id_column, value_columns,
      annotation_IS_vars()
    )]
    find_lambda <- function(.x) {
      if (.x %in% names(additional_agg_lambda)) {
        return(additional_agg_lambda[[.x]])
      }
      defaults <- additional_agg_lambda[[".defaults"]]
      if (is.null(defaults)) {
        return(NULL)
      }
      col_type <- typeof(x[[.x]])
      return(defaults[[col_type]])
    }
    add_cols_lambdas <- purrr::map(additional_cols, find_lambda) %>%
      purrr::set_names(additional_cols)
    add_cols_lambdas <- add_cols_lambdas[purrr::map_lgl(add_cols_lambdas,
                                                        ~ !is.null(.x))]
    # Process
    if (is.null(is_identity_tags) || purrr::is_empty(is_identity_tags)) {
      result <- .sliding_window(
        x = x,
        threshold = threshold,
        keep_criteria = criteria,
        annotated = annotated,
        num_cols = num_cols,
        max_val_col = max_value_column,
        sample_col = sample_id_column,
        req_tags = required_tag_cols,
        add_col_lambdas = add_cols_lambdas,
        produce_map = map_as_file
      )
      recalibr_m <- result$recalibrated_matrix
      maps <- result$map

    } else {
      stopifnot(is.numeric(max_workers))
      max_workers <- max_workers[1]
      is_identity_names <- required_tag_cols %>%
        dplyr::filter(.data$tag != "locus") %>%
        dplyr::pull(.data$names)
      # Split data for parallel execution
      split <- x %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(is_identity_names))) %>%
          dplyr::group_split()
      ## Select only groups with 2 or more rows
      split_to_process <- split[purrr::map_lgl(split, ~ nrow(.x) > 1)]
      # Set up parallel workers
      if (.Platform$OS.type == "windows") {
        p <- BiocParallel::SnowParam(
          workers = max_workers,
          progressbar = getOption("ISAnalytics.verbose"),
          tasks = length(split_to_process),
          exportglobals = TRUE
        )
      } else {
        p <- BiocParallel::MulticoreParam(
          workers = max_workers,
          progressbar = getOption("ISAnalytics.verbose"),
          tasks = length(split_to_process),
          exportglobals = FALSE
        )
      }
      ## Obtain result: list of lists
      result <- BiocParallel::bplapply(
        split_to_process,
        FUN = .sliding_window,
        threshold = threshold,
        keep_criteria = criteria,
        annotated = annotated,
        num_cols = num_cols,
        max_val_col = max_value_column,
        sample_col = sample_id_column,
        req_tags = required_tag_cols,
        add_col_lambdas = add_cols_lambdas,
        produce_map = map_as_file,
        BPPARAM = p
      )
      BiocParallel::bpstop(p)
      ## Obtain single list
      recalibr_m <- purrr::map(result, ~ .x$recalibrated_matrix)
      recalibr_m <- data.table::rbindlist(recalibr_m)
      maps <- purrr::map(result, ~ .x$map)
      maps <- data.table::rbindlist(maps)
      ## Add all rows that were not part of recalibration
      split_fine <- split[purrr::map_lgl(split, ~ !nrow(.x) > 1)]
      if (length(split_fine) > 0) {
        split_fine <- purrr::reduce(
          split_fine,
          function(g1, g2) {
            g1 <- data.table::setDT(g1)
            g2 <- data.table::setDT(g2)
            data.table::rbindlist(list(g1, g2))
          }
        )
        split_fine <- data.table::setDT(split_fine)
        map_fine <- unique(split_fine[, mget(required_tag_cols$names)])
        data.table::setnames(
          x = map_fine,
          old = mget(required_tag_cols$names),
          new = paste0(mget(required_tag_cols$names), "_before")
        )
        map_fine[, c(
          paste0(required_tag_cols$names, "_after")
        ) :=
          purrr::map(
            required_tag_cols$names,
            ~ eval(sym(paste0(.x, "_before")))
          )]
        recalibr_m <- data.table::rbindlist(list(recalibr_m, split_fine))
        maps <- data.table::rbindlist(list(maps, map_fine))
      }
    }
    if (map_as_file) {
        ### Manage file
        withCallingHandlers(
            {
                .write_recalibr_map(maps, file_path)
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
    return(recalibr_m)
}

default_rec_agg_lambdas <- function() {
  list(
    character = ~ paste0(.x, collapse = ";"),
    integer = ~ sum(.x, na.rm = TRUE),
    double = ~ sum(.x, na.rm = TRUE),
    logical = ~ all(.x)
  )
}
