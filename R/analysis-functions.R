#------------------------------------------------------------------------------#
# Analysis functions
#------------------------------------------------------------------------------#
#' Computes the abundance of every integration in the sample.
#'
#' \lifecycle{maturing}
#' Abundance is obtained for every row by calculating the ratio
#' between the single value and the total value for the sample.
#'
#' @details Abundance will be computed upon the user selected columns
#' in the `columns` parameter. For each column a corresponding
#' relative abundance column (and optionally a percentage abundance
#' column) will be produced.
#'
#' @param x An integration matrix - aka a data frame that includes
#' the `mandatory_IS_vars()` as columns
#' @param columns A character vector of column names to process,
#' must be numeric or integer columns
#' @param percentage Add abundance as percentage?
#'
#' @family Analysis functions
#'
#' @importFrom magrittr `%>%`
#' @importFrom tibble is_tibble
#' @import dplyr
#' @importFrom rlang .data eval_tidy parse_expr
#' @importFrom stringr str_replace
#' @return An integration matrix
#' @export
#'
#' @examples
#' path <- system.file("extdata", "ex_annotated_ISMatrix.tsv.xz",
#'     package = "ISAnalytics"
#' )
#' matrix <- import_single_Vispa2Matrix(path)
#' abundance <- compute_abundance(matrix)
compute_abundance <- function(x, columns = "Value", percentage = TRUE) {
    ## Check parameters
    stopifnot(tibble::is_tibble(x))
    stopifnot(is.character(columns))
    if (.check_mandatory_vars(x) == FALSE) {
        stop(.non_ISM_error())
    }
    if (.check_complAmpID(x) == FALSE) {
        stop(.missing_complAmpID_error())
    }
    stopifnot(is.logical(percentage) & length(percentage) == 1)
    if (!all(columns %in% colnames(x))) {
        stop(.missing_user_cols_error())
    }
    purrr::walk(columns, function(col) {
        expr <- rlang::expr(`$`(x, !!col))
        if (!is.numeric(rlang::eval_tidy(expr)) &
            !is.numeric(rlang::eval_tidy(expr))) {
            stop(.non_num_user_cols_error())
        }
    })
    ## Computation
    totals <- x %>%
        dplyr::group_by(.data$CompleteAmplificationID) %>%
        dplyr::summarise(
            dplyr::across(dplyr::all_of(columns),
                sum,
                .names = "{.col}_sum"
            ),
            .groups = "drop"
        )
    abundance_df <- x %>%
        dplyr::left_join(totals, by = "CompleteAmplificationID") %>%
        dplyr::mutate(dplyr::across(dplyr::all_of(columns),
            list(ab = ~ .x / rlang::eval_tidy(
                rlang::parse_expr(
                    paste(
                        dplyr::cur_column(),
                        "sum",
                        sep = "_"
                    )
                )
            )),
            .names = "{.col}_RelAbundance"
        )) %>%
        dplyr::select(-c(dplyr::all_of(paste(columns, "sum", sep = "_")))) %>%
        dplyr::distinct()
    if (percentage == TRUE) {
        abundance_df <- abundance_df %>%
            dplyr::mutate(
                dplyr::across(dplyr::contains("RelAbundance"), ~ .x * 100,
                    .names = "{.col}_PercAbundance"
                )
            ) %>%
            dplyr::rename_with(
                ~ stringr::str_replace(.x, "_RelAbundance", ""),
                dplyr::contains("PercAbundance")
            )
    }
    abundance_df
}


#' obtain a single integration matrix from individual quantification
#' matrices.
#'
#' \lifecycle{maturing}
#' Takes a list of integration matrices referring to different qunatification
#' types and merges them in a single data frame that has multiple
#' value columns, each renamed according to their quantification type
#' of reference.
#'
#' @param x A named list of integration matrices, ideally obtained via
#' \link{import_parallel_Vispa2Matrices_interactive} or
#' \link{import_parallel_Vispa2Matrices_auto}. Names must be
#' quantification types.
#' @param fragmentEstimate The name of the output column for fragment
#' estimate values
#' @param seqCount The name of the output column for sequence
#' count values
#' @param barcodeCount The name of the output column for barcode count
#' values
#' @param cellCount The name of the output column for cell count values
#' @param ShsCount The name of the output column for Shs count values
#'
#' @importFrom purrr walk map2 reduce
#' @importFrom dplyr rename full_join intersect
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data `:=`
#'
#' @family Analysis functions
#'
#' @seealso \link{quantification_types}
#'
#' @return A tibble
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
#'     path, root,
#'     c("fragmentEstimate", "seqCount"), "annotated", 2, NULL, "ANY"
#' )
#' total_matrix <- comparison_matrix(matrices)
#' options(op)
comparison_matrix <- function(x,
    fragmentEstimate = "fragmentEstimate",
    seqCount = "seqCount",
    barcodeCount = "barcodeCount",
    cellCount = "cellCount",
    ShsCount = "ShsCount") {
    stopifnot(is.list(x) & !is.data.frame(x))
    stopifnot(all(names(x) %in% quantification_types()))
    purrr::walk(x, function(m) {
        mand <- .check_mandatory_vars(m)
        amp <- .check_complAmpID(m)
        val <- .check_value_col(m)
        if (any(c(mand, amp, val) == FALSE)) {
            stop(.non_ISM_error())
        }
    })
    stopifnot(is.character(fragmentEstimate) & length(fragmentEstimate) == 1)
    stopifnot(is.character(seqCount) & length(seqCount) == 1)
    stopifnot(is.character(barcodeCount) & length(barcodeCount) == 1)
    stopifnot(is.character(cellCount) & length(cellCount) == 1)
    stopifnot(is.character(ShsCount) & length(ShsCount) == 1)
    param_names <- c(
        fragmentEstimate = fragmentEstimate,
        seqCount = seqCount, barcodeCount = barcodeCount,
        cellCount = cellCount, ShsCount = ShsCount
    )
    x <- purrr::map2(x, names(x), function(matrix, quant_type) {
        quant_name <- param_names[names(param_names) %in% quant_type]
        matrix %>% dplyr::rename(!!quant_name := .data$Value)
    })
    result <- purrr::reduce(x, function(matrix1, matrix2) {
        commoncols <- dplyr::intersect(colnames(matrix1), colnames(matrix2))
        matrix1 %>%
            dplyr::full_join(matrix2, by = commoncols)
    })
    if (any(is.na(result)) & getOption("ISAnalytics.verbose") == TRUE) {
        message(.nas_introduced_msg())
    }
    result
}


#' Separate a multiple-quantification matrix into single quantification
#' matrices.
#'
#' \lifecycle{maturing}
#' The function separates a single multi-quantification integration
#' matrix, obtained via \link{comparison_matrix}, into single
#' quantification matrices as a named list of tibbles.
#'
#' @param x Single integration matrix with multiple quantification
#' value columns, likely obtained via \link{comparison_matrix}.
#' @param fragmentEstimate Name of the fragment estimate values column
#' in input
#' @param seqCount Name of the sequence count values column
#' in input
#' @param barcodeCount Name of the barcode count values column
#' in input
#' @param cellCount Name of the cell count values column
#' in input
#' @param ShsCount Name of the shs count values column
#' in input
#'
#' @importFrom purrr is_empty map set_names
#' @importFrom dplyr rename
#' @importFrom magrittr `%>%`
#'
#' @family Analysis functions
#'
#' @return A named list of tibbles, where names are quantification types
#' @seealso \link{quantification_types}
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
#'     path, root,
#'     c("fragmentEstimate", "seqCount"), "annotated", 2, NULL, "ANY"
#' )
#' total_matrix <- comparison_matrix(matrices)
#' separated_matrix <- separate_quant_matrices(total_matrix)
#' options(op)
separate_quant_matrices <- function(x, fragmentEstimate = "fragmentEstimate",
    seqCount = "seqCount",
    barcodeCount = "barcodeCount",
    cellCount = "cellCount",
    ShsCount = "ShsCount") {
    stopifnot(is.data.frame(x))
    if (.check_mandatory_vars(x) == FALSE) {
        stop(.non_ISM_error())
    }
    if (.check_complAmpID(x) == FALSE) {
        stop(.missing_complAmpID_error())
    }
    num_cols <- .find_exp_cols(x)
    if (purrr::is_empty(num_cols)) {
        stop(.missing_num_cols_error())
    }
    stopifnot(is.character(fragmentEstimate) & length(fragmentEstimate) == 1)
    stopifnot(is.character(seqCount) & length(seqCount) == 1)
    stopifnot(is.character(barcodeCount) & length(barcodeCount) == 1)
    stopifnot(is.character(cellCount) & length(cellCount) == 1)
    stopifnot(is.character(ShsCount) & length(ShsCount) == 1)
    param_col <- c(
        fragmentEstimate = fragmentEstimate,
        seqCount = seqCount, barcodeCount = barcodeCount,
        cellCount = cellCount,
        ShsCount = ShsCount
    )
    to_copy <- if (any(!num_cols %in% param_col)) {
        if (all(!num_cols %in% param_col)) {
            stop(.non_quant_cols_error())
        }
        num_cols[!num_cols %in% param_col]
    }
    num_cols <- param_col[param_col %in% num_cols]
    annot <- if (.is_annotated(x)) {
        annotation_IS_vars()
    } else {
        character(0)
    }
    if (!purrr::is_empty(to_copy) & getOption("ISAnalytics.verbose") == TRUE) {
        message(.non_quant_cols_msg(to_copy))
    }
    separated <- purrr::map(num_cols, function(quant) {
        x[c(
            mandatory_IS_vars(), annot, "CompleteAmplificationID",
            to_copy, quant
        )] %>% dplyr::rename(Value = quant)
    }) %>% purrr::set_names(names(num_cols))
    separated
}


#' Filter data frames with custom predicates
#'
#' @description
#' \lifecycle{experimental}
#' Filter a single data frame or a list of data frames with custom
#' predicates assembled from the function parameters.
#'
#' @details
#' ## A single data frame as input
#'
#' If the user chooses to operate on a single data frame, the other parameters
#' should only be vectors: numeric vector for `threshold` and character
#' vectors for both `cols_to_compare` and `comparators`.
#' A filtering condition is obtained by combining element by element
#' `cols_to_compare` + `comparators` + `threshold` (similarly to the
#' `paste` function). For example:
#'
#' \verb{
#' threshold = c(20, 35, 50)
#' cols_to_compare = c("a", "b", "c")
#' comparators = "<"
#' }
#'
#' given these vectors, the input data frame
#' will be filtered by checking which values in column "a" are less
#' than 20 **AND** which values in column "b" are less than 35 **AND**
#' which values in column "c" are less than 50.
#' Things the user should keep in mind are:
#' * The vectors of length 1 are going to be recycled if one or
#' more parameters are longer (in the example, the `comparators` value)
#' * If vectors are not of length 1 they must have the same length
#' * Columns to compare, of course, need to be included in the
#' input data frame and need to be numeric/integer
#' * The filtering will perform a logical "AND" on all the conditions,
#' only rows that satisfy ALL the conditions are preserved
#'
#' ## A list of data frames as input
#'
#' The input for the function may also be a list of data frames,
#' either named or unnamed.
#'
#' ### Unnamed list
#' If the input is a simple unnamed list, the other parameters should
#' be simple vectors (as for data frames). All the predicates will
#' simply be applied to every data frame in the list: this is useful
#' if it's desirable to filter for the same conditions different data frames
#' that have the same structure but different data.
#'
#' ### Named list
#' It is also possible to filter different data frames with different
#' sets of conditions. Besides having the possibility of defining the
#' other parameters as simple vector, which has the same results as
#' operating on an unnamed list, the user can define the parameters as
#' named lists containing vectors. For example:
#'
#' ```{r}
#'
#' example_df <- tibble::tibble(a = c(20, 30, 40),
#'                              b = c(40, 50, 60),
#'                              c = c("a", "b", "c"),
#'                              d = c(3L, 4L, 5L))
#' example_list <- list(first = example_df,
#'                      second = example_df,
#'                      third = example_df)
#' print(example_list)
#'
#' filtered <- threshold_filter(example_list,
#'                              threshold = list(first = c(20, 60),
#'                                               third = c(25)),
#'                              cols_to_compare = list(first = c("a", "b"),
#'                                                     third = c("a")),
#'                              comparators = list(first = c(">", "<"),
#'                                                 third = c(">=")))
#' print(filtered)
#'
#' ```
#' The above signature will roughly be translated as:
#' * Filter the element "first" in the list by checking that values in
#' column "a" are bigger than 20 AND values in column "b" are less than
#' 60
#' * Don't apply any filter to the element "second" (returns the
#' data frame as is)
#' * Filter the element "third" by checking that values in column "a"
#' are equal or bigger than 25.
#'
#' It is also possible to use some parameters as vectors and some as
#' lists: vectors will be recycled for every element filtered.
#'
#' ```r
#' filtered <- threshold_filter(example_list,
#'                              threshold = list(first = c(20, 60),
#'                                               third = c(25, 65)),
#'                              cols_to_compare = c("a", "b"),
#'                              comparators = list(first = c(">", "<"),
#'                                                 third = c(">=", "<=")))
#' ```
#' In this example, different threshold and comparators will be applied
#' to the same columns in all data frames.
#'
#' Things the user should keep in mind are:
#' * Names for the list parameters must be the same names in the
#' input list
#' * Only elements explicited in list parameters as names will
#' be filtered
#' * Lengths of both vectors and lists must be consistent
#'
#' @param x A data frame or a list of data frames
#' @param threshold A numeric/integer vector or a named list of
#' numeric/integer vectors
#' @param cols_to_compare A character vector or a named list of
#' character vectors
#' @param comparators A character vector or a named list of
#' character vectors. Must be one of the allowed values between
#' `c("<", ">", "==", "!=", ">=", "<=")`
#'
#' @family Analysis functions
#'
#' @return A data frame or a list of data frames
#' @export
#'
#' @examples
#' example_df <- tibble::tibble(
#'     a = c(20, 30, 40),
#'     b = c(40, 50, 60),
#'     c = c("a", "b", "c"),
#'     d = c(3L, 4L, 5L)
#' )
#' example_list <- list(
#'     first = example_df,
#'     second = example_df,
#'     third = example_df
#' )
#'
#' filtered <- threshold_filter(example_list,
#'     threshold = list(
#'         first = c(20, 60),
#'         third = c(25)
#'     ),
#'     cols_to_compare = list(
#'         first = c("a", "b"),
#'         third = c("a")
#'     ),
#'     comparators = list(
#'         first = c(">", "<"),
#'         third = c(">=")
#'     )
#' )
threshold_filter <- function(x,
    threshold,
    cols_to_compare = "Value",
    comparators = ">") {
    stopifnot(is.list(x))
    ### ---- If x is a data frame ---- ###
    if (is.data.frame(x)) {
        return(.tf_data_frame(x, threshold, cols_to_compare, comparators))
    }
    ### ---- If x is a list ---- ###
    return(.tf_list(x, threshold, cols_to_compare, comparators))
}


#' Sorts and keeps the top n integration sites in a data frame.
#'
#' \lifecycle{experimental}
#' The input data frame will be sorted by the highest values in
#' the columns specified and the top n rows will be returned as output.
#' The user can choose to keep additional columns in the output
#' by passing a vector of column names or passing 2 "shortcuts":
#' * `keep` = "everything" keeps all columns in the original data frame
#' * `keep` = "nothing" only keeps the mandatory columns
#' (`mandatory_IS_vars()`) plus the columns in the `columns` parameter.
#'
#' @param x An integration matrix (data frame containing
#' `mandatory_IS_vars()`)
#' @param n How many rows should the output have? Must be numeric
#' or integer and greater than 0
#' @param columns Columns to use for the sorting. If more than a column
#' is supplied primary ordering is done on the first column,
#' secondary ordering on all other columns
#' @param keep Names of the columns to keep besides `mandatory_IS_vars()`
#' and `columns`
#'
#' @family Analysis functions
#'
#' @importFrom dplyr arrange across all_of desc slice_head select
#' @importFrom magrittr `%>%`
#'
#' @return A data frame with `n` rows
#' @export
#'
#' @examples
#' smpl <- tibble::tibble(
#'     chr = c("1", "2", "3", "4", "5", "6"),
#'     integration_locus = c(14536, 14544, 14512, 14236, 14522, 14566),
#'     strand = c("+", "+", "-", "+", "-", "+"),
#'     CompleteAmplificationID = c("ID1", "ID2", "ID1", "ID1", "ID3", "ID2"),
#'     Value = c(3, 10, 40, 2, 15, 150),
#'     Value2 = c(456, 87, 87, 9, 64, 96),
#'     Value3 = c("a", "b", "c", "d", "e", "f")
#' )
#' top <- top_integrations(smpl,
#'     n = 3,
#'     columns = c("Value", "Value2"),
#'     keep = "nothing"
#' )
top_integrations <- function(x, n = 50, columns = "RelAbundance",
    keep = "everything") {
    stopifnot(is.data.frame(x))
    stopifnot(is.numeric(n) & length(n) == 1 & n > 0)
    stopifnot(is.character(keep))
    stopifnot(is.character(columns))
    if (!.check_mandatory_vars(x)) {
        stop(.non_ISM_error())
    }
    if (!all(columns %in% colnames(x))) {
        stop(.missing_user_cols_error())
    }
    if (!(all(keep == "everything") || all(keep == "nothing"))) {
        if (any(!keep %in% colnames(x))) {
            stop(.missing_user_cols_error())
        }
    }
    essential_cols <- c(mandatory_IS_vars(), columns)
    to_keep <- if (all(keep == "everything")) {
        colnames(x)[!colnames(x) %in% essential_cols]
    } else if (all(keep == "nothing")) {
        character(0)
    } else {
        keep[!keep %in% essential_cols]
    }
    result <- x %>%
        dplyr::arrange(dplyr::across(
            dplyr::all_of(columns),
            dplyr::desc
        )) %>%
        dplyr::slice_head(n = n) %>%
        dplyr::select(dplyr::all_of(c(essential_cols, to_keep)))
    return(result)
}


#' Computes user specified functions on numerical columns and updates
#' the metadata data frame accordingly.
#'
#' @description
#' \lifecycle{experimental}
#' The function operates on a data frame by grouping the content by
#' the sample key and computing every function specified on every
#' column in the `value_columns` parameter. After that the metadata
#' data frame is updated by including the computed results as columns
#' for the corresponding key.
#' For this reason it's required that both `x` and `metadata` have the
#' same sample key, and it's particularly important if the user is
#' working with previously aggregated data.
#' For example:
#'
#' ```r
#' ### Importing association file and matrices
#' path_AF <- system.file("extdata", "ex_association_file.tsv",
#' package = "ISAnalytics")
#' root_correct <- system.file("extdata", "fs.zip",
#' package = "ISAnalytics")
#' root_correct <- unzip_file_system(root_correct, "fs")
#'
#' association_file <- import_association_file(path_AF, root_correct)
#' matrices <- import_parallel_Vispa2Matrices_auto(
#' association_file = association_file , root = NULL,
#' quantification_type = c("seqCount","fragmentEstimate"),
#' matrix_type = "annotated", workers = 2, patterns = NULL,
#' matching_opt = "ANY")
#'
#' ### Aggregating data (both by same key)
#' aggreggated_x <- aggregate_values_by_key(matrices$seqCount,
#' association_file)
#' aggregated_meta <- aggregate_metadata(association_file)
#'
#' ### Sample statistics
#' sample_stats <- sample_statistics(x = aggregated_x,
#' metadata = aggregated_meta,
#' sample_key = c("SubjectID", "CellMarker","Tissue", "TimePoint"))
#'
#' ```
#' @param x A data frame
#' @param metadata The metadata data frame
#' @param sample_key Character vector representing the key for identifying
#' a sample
#' @param value_columns THe name of the columns to be computed,
#' must be numeric or integer
#' @param functions A named list of function or purrr-style lambdas
#'
#' @family Analysis functions
#' @importFrom rlang eval_tidy expr
#' @importFrom purrr is_function is_formula
#' @import dplyr
#' @importFrom magrittr `%>%`
#'
#' @return A list with modified x and metadata data frames
#' @export
#'
#' @examples
#' op <- options(ISAnalytics.widgets = FALSE)
#'
#' path_AF <- system.file("extdata", "ex_association_file.tsv",
#'     package = "ISAnalytics"
#' )
#' root_correct <- system.file("extdata", "fs.zip",
#'     package = "ISAnalytics"
#' )
#' root_correct <- unzip_file_system(root_correct, "fs")
#'
#' association_file <- import_association_file(path_AF, root_correct)
#' matrices <- import_parallel_Vispa2Matrices_auto(
#'     association_file = association_file, root = NULL,
#'     quantification_type = c("seqCount", "fragmentEstimate"),
#'     matrix_type = "annotated", workers = 2, patterns = NULL,
#'     matching_opt = "ANY"
#' )
#'
#' stats <- sample_statistics(matrices$seqCount, association_file)
#' options(op)
sample_statistics <- function(x, metadata,
    sample_key = "CompleteAmplificationID",
    value_columns = "Value",
    functions = default_stats()) {
    stopifnot(is.data.frame(x))
    stopifnot(is.data.frame(metadata))
    stopifnot(is.character(sample_key))
    stopifnot(is.character(value_columns))
    stopifnot(is.list(functions))
    if (!all(sample_key %in% colnames(x))) {
        stop(paste("Key columns not found in the data frame"))
    }
    if (!all(sample_key %in% colnames(metadata))) {
        stop(paste("Key columns not found in metadata"))
    }
    if (!all(value_columns %in% colnames(x))) {
        stop(paste("Value columns not found in the data frame"))
    }
    purrr::walk(value_columns, function(col) {
        expr <- rlang::expr(`$`(x, !!col))
        if (!is.numeric(rlang::eval_tidy(expr)) &&
            !is.integer(rlang::eval_tidy(expr))) {
            stop(paste("Some or all of value columns are not numeric"))
        }
    })
    purrr::walk(functions, function(f) {
        if (!(purrr::is_function(f) | purrr::is_formula(f))) {
            stop(paste(
                "The function parameter should contain a list",
                "of either functions or formulas.",
                "See ?sample_statistics for details"
            ))
        }
    })
    result <- x %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(sample_key))) %>%
        dplyr::summarise(dplyr::across(
            .cols = dplyr::all_of(value_columns),
            .fns = functions
        ),
        .groups = "drop"
        )

    updated_meta <- metadata %>% dplyr::left_join(result, by = sample_key)
    return(list(x = result, metadata = updated_meta))
}


#' A set of pre-defined functions for `sample_statistics`.
#'
#' @return A named list of functions/purrr-style lambdas
#' @export
#'
#' @family Analysis functions helpers
#'
#'
#' @examples
#' default_stats()
default_stats <- function() {
    list(
        diversity = vegan::diversity,
        sum = ~ sum(.x, na.rm = TRUE),
        count = length,
        describe = psych::describe
    )
}
