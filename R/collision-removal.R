#------------------------------------------------------------------------------#
# Collision removal functions
#------------------------------------------------------------------------------#

#' Identifies and removes collisions based on the sequence count matrix.
#'
#' A collision is an integration (aka a unique combination of
#' chr, integration_locus and strand) which is observed in more than one
#' independent sample (a unique pair of ProjectID and SubjectID). The function
#' tries to decide to which subject an integration should be assigned and if no
#' decision can be taken, the integration is completely removed from the data
#' frame.
#'
#' @details If you don't want the function to show details and messages do:
#' \code{options(ISAnalitics.verbose = FALSE)}.
#' To restore to the original value:
#' \code{options(ISAnalitics.verbose = TRUE)}.
#' For more details on how to use collision removal functionality:
#' \code{vignette("Collision removal functionality", package = "ISAnalytics")}
#'
#' @param x A named list of matrices (names must be quantification types), or
#' a single ISADataFrame representing the sequence count matrix of interest.
#' @param association_file The association file imported via
#' `import_association_file`
#' @param date_col The date column that should be considered for the analysis.
#' Must be one value in `date_columns_coll()`
#' @param reads_ratio A single numeric value that represents the ratio that has
#' to be considered when deciding between seqCount value.
#' @family Collision removal
#' @importFrom dplyr bind_rows all_of select
#' @seealso \code{\link{date_columns_coll}}
#'
#' @return A list of ISADataframes with removed collisions
#' @export
#'
#' @examples
#' path <- system.file("extdata", "ex_association_file.tsv",
#'     package = "ISAnalytics"
#' )
#' root_pth <- system.file("extdata", "fs.zip", package = "ISAnalytics")
#' root <- unzip_file_system(root_pth, "fs")
#' association_file <- import_association_file(path, root)
#' matrices <- import_parallel_Vispa2Matrices_auto(
#'     association_file, NULL,
#'     c("fragmentEstimate", "seqCount"), "annotated", 2, NULL, "ANY"
#' )
#' matrices <- remove_collisions(matrices, association_file)
remove_collisions <- function(x,
    association_file,
    date_col = "SequencingDate",
    reads_ratio = 10) {
    # Check basic parameter correctness
    stopifnot(is.list(x) & !is.null(names(x)))
    if (is.ISADataFrame(x)) {
        x <- list(seqCount = x)
    }
    stopifnot(all(names(x) %in% quantification_types()))
    all_ISAdf <- purrr::map_lgl(x, is.ISADataFrame)
    if (!all(all_ISAdf)) {
        stop("x contains elements that are not ISADataFrame objects")
    }
    stopifnot(tibble::is_tibble(association_file))
    stopifnot(is.character(date_col) & length(date_col) == 1)
    stopifnot(date_col %in% date_columns_coll())
    stopifnot(is.integer(reads_ratio) | is.numeric(reads_ratio))
    stopifnot(length(reads_ratio) == 1)

    # Check association file correctness
    if (.check_af_correctness(association_file) == FALSE) {
        stop("Malformed association file: one or more columns are missing")
    }

    # Check date_col
    if (any(is.na(association_file[date_col]))) {
        stop(paste(
            "Selected date column contains NA values, please check",
            "and fix the association file"
        ))
    }

    # Check imported matrices vs association file
    ## remove_collisions requires seqCount matrix, check if the list contains
    ## one
    if ((!"seqCount" %in% names(x)) ||
        nrow(x$seqCount) == 0) {
        stop(paste(
            "Sequence count data frame is required for collision removal",
            "but none was detected in x"
        ))
    }
    seq_count_df <- x$seqCount
    ## Check if association file contains all info relative to content the of
    ## the matrix
    all_contained <- all(seq_count_df$CompleteAmplificationID %in%
        association_file$CompleteAmplificationID)
    if (!all_contained) {
        missing <- which(!seq_count_df$CompleteAmplificationID %in%
            association_file$CompleteAmplificationID)
        cat("Missing info for these observations: ", sep = "\n")
        print(seq_count_df[missing, ])
        stop("The association file is missing needed info on some experiments")
    }
    ## Check if association file contains more info than matrix
    verbose <- getOption("ISAnalytics.verbose")
    if (verbose == TRUE) {
        not_included <- .check_same_info(association_file, seq_count_df)
        if (nrow(not_included) > 0) {
            message(paste("Found additional data relative to some projects",
                          "that are not included in the imported matrices.",
                          "Here is a summary",
                collapse = "\n"
            ))
            print(not_included)
        }
    }

    # Begin workflow
    ## Join based on CompleteAmplificationID
    joined <- .join_matrix_af(seq_count_df, association_file, date_col)
    if (verbose == TRUE) {
        message("Identifying collisions...")
    }
    splitted_df <- .identify_independent_samples(joined)

    # Remove collisions
    if (verbose == TRUE) {
        message("Processing collisions...")
    }
    fixed_collisions <- .process_collisions(
        splitted_df$collisions,
        date_col, reads_ratio
    )
    reassigned <- fixed_collisions$reassigned
    removed <- fixed_collisions$removed
    fixed_collisions <- fixed_collisions$coll
    final_matr <- fixed_collisions %>%
        dplyr::bind_rows(splitted_df$non_collisions) %>%
        dplyr::select(dplyr::all_of(colnames(seq_count_df)))
    if (verbose == TRUE) {
        summary_tbl <- .summary_table(
            before = joined, after = final_matr,
            association_file = association_file
        )
        widget <- .summary_collisions_widget(removed, reassigned, summary_tbl,
            tot_rows = nrow(joined),
            collision_rows = nrow(
                splitted_df$collisions
            )
        )
        print(widget)
    }
    ## Align other matrices if present
    if (length(x) > 1) {
        if (verbose == TRUE) {
            message("Realigning matrices...")
        }
        other <- x[!names(x) %in% "seqCount"]
        other <- realign_after_collisions(final_matr, other)
        if (verbose == TRUE) {
            message("Finished!")
        }
        lst <- list(seqCount = final_matr)
        lst <- append(lst, other)
        return(lst)
    } else {
        if (verbose == TRUE) {
            message(paste(
                "Finished! You provided a single sequence count as matrix,",
                "to re-align other related matrices see",
                "?realign_after_collisions"
            ))
            return(list(seqCount = final_matr))
        }
    }
}

#' Re-aligns matrices of other quantification types based on the processed
#' sequence count matrix.
#'
#' This function should be used to keep data consistent among the same analysis:
#' if for some reason you removed the collisions by passing only the sequence
#' count matrix to the `remove_collisions` function, you should call this
#' function afterwards, providing a list of other quantification matrices.
#' NOTE: if you provided a list of several quantification types to
#' `remove_collisions` before, there is no need to call this function.
#'
#' @details For more details on how to use collision removal functionality:
#' \code{vignette("Collision removal functionality", package = "ISAnalytics")}
#'
#' @param sc_matrix The sequence count matrix already processed for collisions
#' via `remove_collisions`
#' @param other_matrices A named list of matrices to re-align. Names in the list
#' must be quantification types (\code{quantification_types()}) except
#' "seqCount".
#' @importFrom dplyr semi_join
#' @importFrom purrr map_lgl
#' @family Collision removal
#' @seealso \code{\link{remove_collisions}}
#'
#' @return A named list with re-aligned matrices
#' @export
#'
#' @examples
#' path <- system.file("extdata", "ex_association_file.tsv",
#'     package = "ISAnalytics"
#' )
#' root_pth <- system.file("extdata", "fs.zip", package = "ISAnalytics")
#' root <- unzip_file_system(root_pth, "fs")
#' association_file <- import_association_file(path, root)
#' matrices <- import_parallel_Vispa2Matrices_auto(
#'     association_file, NULL,
#'     c("fragmentEstimate", "seqCount"), "annotated", 2, NULL, "ANY"
#' )
#' sc_matrix <- remove_collisions(matrices$seqCount, association_file)
#' others <- matrices[!names(matrices) %in% "seqCount"]
#' aligned_matrices <- realign_after_collisions(sc_matrix$seqCount, others)
realign_after_collisions <- function(sc_matrix, other_matrices) {
    stopifnot(is.list(other_matrices) & !is.null(names(other_matrices)))
    stopifnot(all(names(other_matrices) %in% quantification_types()))
    all_ISAdf <- purrr::map_lgl(other_matrices, is.ISADataFrame)
    if (!all(all_ISAdf)) {
        stop(paste(
            "other_matrices list contains elements that are not",
            "ISADataFrame objects"
        ))
    }
    realigned <- purrr::map(other_matrices, function(x) {
        x %>% dplyr::semi_join(sc_matrix,
            by = c(
                "chr", "integration_locus",
                "strand", "CompleteAmplificationID"
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
