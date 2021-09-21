#------------------------------------------------------------------------------#
# Analysis functions
#------------------------------------------------------------------------------#

#' Computes the abundance for every integration event in the input data frame.
#'
#' \lifecycle{stable}
#' Abundance is obtained for every integration event by calculating the ratio
#' between the single value and the total value for the given group.
#'
#' @details Abundance will be computed upon the user selected columns
#' in the `columns` parameter. For each column a corresponding
#' relative abundance column (and optionally a percentage abundance
#' column) will be produced.
#'
#' @param x An integration matrix - aka a data frame that includes
#' the `mandatory_IS_vars()` as columns. The matrix can either be aggregated
#' (via `aggregate_values_by_key()`) or not.
#' @param columns A character vector of column names to process,
#' must be numeric or integer columns
#' @param percentage Add abundance as percentage?
#' @param key The key to group by when calculating totals
#' @param keep_totals A value between `TRUE`, `FALSE` or `df`. If `TRUE`,
#' the intermediate totals for each group will be kept in the output
#' data frame as a dedicated column with a trailing "_tot". If `FALSE`,
#' totals won't be included in the output data frame. If `df`, the totals
#' are returned to the user as a separate data frame, together with the
#' abundance data frame.
#'
#' @family Analysis functions
#'
#' @importFrom magrittr `%>%`
#' @importFrom dplyr group_by across all_of summarise left_join mutate
#' @importFrom dplyr cur_column distinct select contains rename_with
#' @importFrom rlang .data eval_tidy parse_expr abort
#' @importFrom purrr map_lgl
#' @importFrom stringr str_replace
#' @return Either a single data frame with computed abundance values or
#' a list of 2 data frames (abundance_df, quant_totals)
#' @export
#'
#' @examples
#' data("integration_matrices", package = "ISAnalytics")
#' abund <- compute_abundance(
#'     x = integration_matrices,
#'     columns = "fragmentEstimate",
#'     key = "CompleteAmplificationID"
#' )
#' head(abund)
compute_abundance <- function(x,
    columns = c("fragmentEstimate_sum"),
    percentage = TRUE,
    key = c("SubjectID", "CellMarker", "Tissue", "TimePoint"),
    keep_totals = FALSE) {
    ## Check parameters
    stopifnot(is.data.frame(x))
    stopifnot(is.character(columns))
    stopifnot(is.character(key))
    if (.check_mandatory_vars(x) == FALSE) {
        rlang::abort(.missing_mand_vars())
    }
    stopifnot(is.logical(percentage) & length(percentage) == 1)
    if (!all(columns %in% colnames(x)) | !all(key %in% colnames(x))) {
        missing_cols <- c(
            columns[!columns %in% colnames(x)],
            key[!key %in% colnames(x)]
        )
        rlang::abort(.missing_user_cols_error(missing_cols))
    }
    non_num_cols <- purrr::map_lgl(columns, function(col) {
        expr <- rlang::expr(`$`(x, !!col))
        if (is.numeric(rlang::eval_tidy(expr)) |
            is.integer(rlang::eval_tidy(expr))) {
            return(FALSE)
        } else {
            return(TRUE)
        }
    })
    if (any(non_num_cols)) {
        rlang::abort(.non_num_user_cols_error(columns[non_num_cols]))
    }
    stopifnot(is.logical(keep_totals) || keep_totals == "df")
    ## Computation
    ### Computes totals for each group defined by key
    totals <- x %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(key))) %>%
        dplyr::summarise(
            dplyr::across(dplyr::all_of(columns),
                sum,
                .names = "{.col}_tot"
            ),
            .groups = "drop"
        )
    ### Computes abundance as value (for each col) / total of the corresponding
    ### group (defined by key)
    abundance_df <- x %>%
        dplyr::left_join(totals, by = key) %>%
        dplyr::mutate(dplyr::across(dplyr::all_of(columns),
            list(ab = ~ .x / rlang::eval_tidy(
                rlang::parse_expr(
                    paste(
                        dplyr::cur_column(),
                        "tot",
                        sep = "_"
                    )
                )
            )),
            .names = "{.col}_RelAbundance"
        )) %>%
        dplyr::distinct()
    if (keep_totals == FALSE || keep_totals == "df") {
        abundance_df <- abundance_df %>%
            dplyr::select(-c(dplyr::all_of(paste(columns, "tot", sep = "_"))))
    }
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
    if (keep_totals == "df") {
        return(list(abundance_df = abundance_df, quant_totals = totals))
    } else {
        return(abundance_df)
    }
}


#' obtain a single integration matrix from individual quantification
#' matrices.
#'
#' \lifecycle{stable}
#' Takes a list of integration matrices referring to different quantification
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
#'     mode = "AUTO", report_path = NULL, multi_quant_matrix = FALSE
#' )
#' multi_quant <- comparison_matrix(matrices)
#' head(multi_quant)
comparison_matrix <- function(x,
    fragmentEstimate = "fragmentEstimate",
    seqCount = "seqCount",
    barcodeCount = "barcodeCount",
    cellCount = "cellCount",
    ShsCount = "ShsCount") {
    stopifnot(is.list(x) & !is.data.frame(x))
    stopifnot(all(names(x) %in% quantification_types()))
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
    na_introduced <- purrr::map_lgl(param_names, function(p) {
        any(is.na(result[[p]]))
    })
    if (any(na_introduced) & getOption("ISAnalytics.verbose") == TRUE) {
        rlang::inform(.nas_introduced_msg())
    }
    result
}


#' Separate a multiple-quantification matrix into single quantification
#' matrices.
#'
#' \lifecycle{stable}
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
#' @param key Key columns to perform the joining operation
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
#' data("integration_matrices", package = "ISAnalytics")
#' separated <- separate_quant_matrices(
#'     integration_matrices
#' )
#' separated
separate_quant_matrices <- function(x,
    fragmentEstimate = "fragmentEstimate",
    seqCount = "seqCount",
    barcodeCount = "barcodeCount",
    cellCount = "cellCount",
    ShsCount = "ShsCount",
    key = c(
        mandatory_IS_vars(),
        annotation_IS_vars(),
        "CompleteAmplificationID"
    )) {
    stopifnot(is.data.frame(x))
    if (!all(key %in% colnames(x))) {
        rlang::abort(.missing_user_cols_error(key[!key %in% colnames(x)]))
    }
    num_cols <- .find_exp_cols(x, key)
    if (purrr::is_empty(num_cols)) {
        rlang::abort(.missing_num_cols_error())
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
            rlang::abort(.non_quant_cols_error())
        }
        num_cols[!num_cols %in% param_col]
    }
    num_cols <- param_col[param_col %in% num_cols]
    if (!purrr::is_empty(to_copy) & getOption("ISAnalytics.verbose") == TRUE) {
        rlang::inform(.non_quant_cols_msg(to_copy))
    }
    separated <- purrr::map(num_cols, function(quant) {
        x %>%
            dplyr::select(dplyr::all_of(c(key, to_copy, quant))) %>%
            dplyr::rename(Value = quant)
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
#' threshold = list(first = c(20, 60),
#' third = c(25)),
#' cols_to_compare = list(first = c("a", "b"),
#' third = c("a")),
#' comparators = list(first = c(">", "<"),
#' third = c(">=")))
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
#' threshold = list(first = c(20, 60),
#' third = c(25, 65)),
#' cols_to_compare = c("a", "b"),
#' comparators = list(first = c(">", "<"),
#' third = c(">=", "<=")))
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
#' filtered
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


#' Sorts and keeps the top n integration sites based on the values
#' in a given column.
#'
#' \lifecycle{experimental}
#' The input data frame will be sorted by the highest values in
#' the columns specified and the top n rows will be returned as output.
#' The user can choose to keep additional columns in the output
#' by passing a vector of column names or passing 2 "shortcuts":
#' * `keep = "everything"` keeps all columns in the original data frame
#' * `keep = "nothing"` only keeps the mandatory columns
#' (`mandatory_IS_vars()`) plus the columns in the `columns` parameter.
#'
#' @param x An integration matrix (data frame containing
#' `mandatory_IS_vars()`)
#' @param n How many integrations should be sliced (in total or
#'  for each group)? Must be numeric
#' or integer and greater than 0
#' @param columns Columns to use for the sorting. If more than a column
#' is supplied primary ordering is done on the first column,
#' secondary ordering on all other columns
#' @param keep Names of the columns to keep besides `mandatory_IS_vars()`
#' and `columns`
#' @param key Either `NULL` or a character vector of column names to group
#' by. If not `NULL` the input will be grouped and the top fraction will
#' be extracted from each group.
#'
#' @family Analysis functions
#'
#' @importFrom magrittr `%>%`
#' @importFrom rlang abort
#'
#' @return Either a data frame with at most n rows or
#' a data frames with at most n*(number of groups) rows.
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
#' top_key <- top_integrations(smpl,
#'     n = 3,
#'     columns = "Value",
#'     keep = "Value2",
#'     key = "CompleteAmplificationID"
#' )
top_integrations <- function(x,
    n = 20,
    columns = "fragmentEstimate_sum_RelAbundance",
    keep = "everything",
    key = NULL) {
    stopifnot(is.data.frame(x))
    stopifnot(is.numeric(n) & length(n) == 1 & n > 0)
    stopifnot(is.character(keep))
    stopifnot(is.character(columns))
    stopifnot(is.null(key) || is.character(key))
    if (!.check_mandatory_vars(x)) {
        rlang::abort(.missing_mand_vars())
    }
    if (!all(columns %in% colnames(x))) {
        rlang::abort(.missing_user_cols_error(
            columns[!columns %in% colnames(x)]
        ))
    }
    if (!(all(keep == "everything") || all(keep == "nothing"))) {
        if (any(!keep %in% colnames(x))) {
            rlang::abort(.missing_user_cols_error(
                keep[!keep %in% colnames(x)]
            ))
        }
    }
    if (!is.null(key)) {
        if (!all(key %in% colnames(x))) {
            rlang::abort(.missing_user_cols_error(
                key[!key %in% colnames(x)]
            ))
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
    if (!is.null(key)) {
        result <- x %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(key))) %>%
            dplyr::arrange(dplyr::across(
                dplyr::all_of(columns),
                dplyr::desc
            ), .by_group = TRUE) %>%
            dplyr::slice_head(n = n) %>%
            dplyr::select(dplyr::all_of(c(key, essential_cols, to_keep))) %>%
            dplyr::ungroup()
        return(result)
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
#' data("integration_matrices", package = "ISAnalytics")
#' data("association_file", package = "ISAnalytics")
#' aggreg <- aggregate_values_by_key(
#'  x = integration_matrices,
#'  association_file = association_file,
#'  value_cols = c("seqCount", "fragmentEstimate")
#' )
#' aggreg_meta <- aggregate_metadata(association_file = association_file)
#'
#' sample_stats <- sample_statistics(x = aggreg,
#' metadata = aggreg_meta,
#' value_columns = c("seqCount", "fragmentEstimate"),
#' sample_key = c("SubjectID", "CellMarker","Tissue", "TimePoint"))
#'
#' ```
#' @param x A data frame
#' @param metadata The metadata data frame
#' @param sample_key Character vector representing the key for identifying
#' a sample
#' @param value_columns The name of the columns to be computed,
#' must be numeric or integer
#' @param functions A named list of function or purrr-style lambdas
#' @param add_integrations_count Add the count of distinct integration sites
#' for each group? Can be computed only if `x` contains the mandatory columns
#' `chr`, `integration_locus`, `strand`
#'
#' @family Analysis functions
#' @importFrom rlang eval_tidy expr abort .data sym inform
#' @importFrom purrr is_function is_formula map_lgl walk map set_names
#' @importFrom dplyr group_by across all_of summarise rename_with bind_cols
#' @importFrom dplyr n_distinct left_join
#' @importFrom magrittr `%>%`
#'
#' @return A list with modified x and metadata data frames
#' @export
#'
#' @examples
#' data("integration_matrices", package = "ISAnalytics")
#' data("association_file", package = "ISAnalytics")
#' stats <- sample_statistics(
#'     x = integration_matrices,
#'     metadata = association_file,
#'     value_columns = c("seqCount", "fragmentEstimate")
#' )
#' stats
sample_statistics <- function(x, metadata,
    sample_key = "CompleteAmplificationID",
    value_columns = "Value",
    functions = default_stats(),
    add_integrations_count = TRUE) {
    stopifnot(is.data.frame(x))
    stopifnot(is.data.frame(metadata))
    stopifnot(is.character(sample_key))
    stopifnot(is.character(value_columns))
    stopifnot(is.list(functions))
    stopifnot(is.logical(add_integrations_count))
    if (!all(c(sample_key, value_columns) %in% colnames(x))) {
        rlang::abort(.missing_user_cols_error(c(sample_key, value_columns)[
            !c(sample_key, value_columns) %in% colnames(x)
        ]))
    }
    if (!all(sample_key %in% colnames(metadata))) {
        rlang::abort(.missing_user_cols_meta_error(sample_key[
            !sample_key %in% colnames(x)
        ]))
    }
    vcols_are_numeric <- purrr::map_lgl(value_columns, function(col) {
        expr <- rlang::expr(`$`(x, !!col))
        is.numeric(rlang::eval_tidy(expr)) ||
            is.integer(rlang::eval_tidy(expr))
    })
    if (any(vcols_are_numeric == FALSE)) {
        rlang::abort(.non_num_user_cols_error(value_columns[!vcols_are_numeric]))
    }
    purrr::walk(functions, function(f) {
        if (!(purrr::is_function(f) | purrr::is_formula(f))) {
            rlang::abort(.non_function_elem_error())
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

    ## Flatten nested data frames
    df_cols <- purrr::map_lgl(result, ~ is.data.frame(.x))
    if (any(df_cols == TRUE)) {
        df_cols <- df_cols[df_cols]
        df_cols <- names(df_cols)
        dfss <- purrr::map(df_cols, function(col) {
            exp <- rlang::expr(`$`(result, !!col))
            df <- rlang::eval_tidy(exp)
            df <- df %>% dplyr::rename_with(.fn = ~ paste0(col, "_", .x))
            df
        }) %>% purrr::set_names(df_cols)
        for (dfc in df_cols) {
            result <- result %>%
                dplyr::select(-dplyr::all_of(dfc)) %>%
                dplyr::bind_cols(dfss[[dfc]])
        }
    }

    if (add_integrations_count) {
        if (all(mandatory_IS_vars() %in% colnames(x))) {
            mand_sym <- purrr::map(mandatory_IS_vars(), rlang::sym)
            nIS <- x %>%
                dplyr::group_by(dplyr::across(dplyr::all_of(sample_key))) %>%
                dplyr::summarise(
                    nIS = dplyr::n_distinct(!!!mand_sym),
                    .groups = "drop"
                )
            result <- result %>%
                dplyr::left_join(nIS, by = sample_key)
        } else {
            if (getOption("ISAnalytics.verbose")) {
                rlang::inform(.inform_skip_count_is())
            }
        }
    }

    updated_meta <- metadata %>% dplyr::left_join(result, by = sample_key)
    return(list(x = result, metadata = updated_meta))
}


#' Grubbs test for Common Insertion Sites (CIS).
#'
#' \lifecycle{stable}
#' Statistical approach for the validation of common insertion sites
#' significance based on the comparison of the integration frequency
#' at the CIS gene with respect to other genes contained in the
#' surrounding genomic regions. For more details please refer to
#' this paper:
#' <`r .lentiviral_CIS_paper()`>
#'
#' @details
#' ## Genomic annotation file
#' This file is a data base, or more simply a .tsv file to import, with
#' genes annotation for the specific genome. The annotations for the
#' human genome (hg19) and murine genome (mm9 and mm10) are already
#' included in this package: to use one of them just
#' set the argument `genomic_annotation_file` to either `"hg19"`,
#' `"mm9"` or `"mm10"`.
#' If for any reason the user is performing an analysis on another genome,
#' this file needs to be changed respecting the USCS Genome Browser
#' format, meaning the input file headers should include:
#'
#' `r refGene_table_cols()`
#'
#' @param x An integration matrix, must include the `mandatory_IS_vars()`
#' columns and the `annotation_IS_vars()` columns
#' @param genomic_annotation_file Database file for gene annotation,
#' see details.
#' @param grubbs_flanking_gene_bp Number of base pairs flanking a gene
#' @param threshold_alpha Significance threshold
#' @param by Either `NULL` or a character vector of column names. If not
#' NULL, the function will perform calculations for each group and return
#' a list of data frames with the results. E.g. for `by = "SubjectID"`,
#' CIS will be computed for each distinct SubjectID found in the table
#' (of course, "SubjectID" column must be included in the input data frame).
#'
#' @family Analysis functions
#'
#' @importFrom tibble as_tibble
#' @importFrom rlang .data abort current_env eval_tidy sym
#' @importFrom magrittr `%>%`
#' @importFrom utils read.csv
#' @importFrom stringr str_replace_all
#' @importFrom tidyr unite
#' @importFrom purrr set_names map
#'
#' @return A data frame
#' @export
#'
#' @examples
#' data("integration_matrices", package = "ISAnalytics")
#' cis <- CIS_grubbs(integration_matrices)
#' head(cis)
CIS_grubbs <- function(x,
    genomic_annotation_file = "hg19",
    grubbs_flanking_gene_bp = 100000,
    threshold_alpha = 0.05,
    by = NULL) {
    # Check x has the correct structure
    stopifnot(is.data.frame(x))
    if (!all(mandatory_IS_vars() %in% colnames(x))) {
        rlang::abort(.missing_mand_vars())
    }
    if (!.is_annotated(x)) {
        rlang::abort(.missing_annot())
    }
    # Check other parameters
    stopifnot(is.character(genomic_annotation_file))
    genomic_annotation_file <- genomic_annotation_file[1]
    if (genomic_annotation_file %in% c("hg19", "mm9", "mm10")) {
        gen_file <- paste0("refGenes_", genomic_annotation_file)
        utils::data(list = gen_file, envir = rlang::current_env())
        refgenes <- rlang::eval_tidy(rlang::sym(gen_file))
    } else {
        stopifnot(file.exists(genomic_annotation_file))
        # Determine file extension
        ext <- .check_file_extension(genomic_annotation_file)
        # Try to import annotation file
        if (ext == "tsv") {
            refgenes <- utils::read.csv(
                file = genomic_annotation_file,
                header = TRUE, fill = TRUE, sep = "\t",
                check.names = FALSE,
                na.strings = c("NONE", "NA", "NULL", "NaN", "")
            )
            # Check annotation file format
            if (!all(refGene_table_cols() %in% colnames(refgenes))) {
                rlang::abort(.non_standard_annotation_structure())
            }
            refgenes <- tibble::as_tibble(refgenes) %>%
                dplyr::mutate(chrom = stringr::str_replace_all(
                    .data$chrom,
                    "chr", ""
                ))
        } else if (ext == "csv") {
            refgenes <- utils::read.csv(
                file = genomic_annotation_file,
                header = TRUE, fill = TRUE,
                check.names = FALSE,
                na.strings = c("NONE", "NA", "NULL", "NaN", "")
            )
            # Check annotation file format
            if (!all(refGene_table_cols() %in% colnames(refgenes))) {
                rlang::abort(.non_standard_annotation_structure())
            }
            refgenes <- tibble::as_tibble(refgenes) %>%
                dplyr::mutate(chrom = stringr::str_replace_all(
                    .data$chrom,
                    "chr", ""
                ))
        } else {
            gen_file_err <- paste(
                "The genomic annotation file must be either in",
                ".tsv or .csv format (compressed or not)"
            )
            rlang::abort(gen_file_err)
        }
    }
    stopifnot(is.numeric(grubbs_flanking_gene_bp) |
        is.integer(grubbs_flanking_gene_bp))
    grubbs_flanking_gene_bp <- grubbs_flanking_gene_bp[1]
    stopifnot(is.numeric(threshold_alpha))
    threshold_alpha <- threshold_alpha[1]
    stopifnot(is.null(by) || is.character(by))
    if (!all(by %in% colnames(x))) {
        rlang::abort(.missing_user_cols_error(by[!by %in% colnames(x)]))
    }
    result <- if (is.null(by)) {
        .cis_grubb_calc(
            x = x,
            refgenes = refgenes,
            grubbs_flanking_gene_bp = grubbs_flanking_gene_bp,
            threshold_alpha = threshold_alpha
        )
    } else {
        grouped <- x %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(by)))
        group_keys <- grouped %>%
            dplyr::group_keys() %>%
            tidyr::unite(col = "id", dplyr::everything()) %>%
            dplyr::pull(.data$id)
        split <- x %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(by))) %>%
            dplyr::group_split() %>%
            purrr::set_names(group_keys)
        purrr::map(split, ~ .cis_grubb_calc(
            x = .x,
            refgenes = refgenes,
            grubbs_flanking_gene_bp = grubbs_flanking_gene_bp,
            threshold_alpha = threshold_alpha
        ))
    }
    return(result)
}

#' Integrations cumulative count in time by sample
#'
#' \lifecycle{experimental}
#' This function computes the cumulative number of integrations
#' observed in each sample at different time points by assuming that
#' if an integration is observed at time point "t" then it is also observed in
#' time point "t+1".
#'
#' @details
#' ## Input data frame
#' The user can provide as input for the `x` parameter both a simple
#' integration matrix AND setting the `aggregate` parameter to TRUE,
#' or provide an already aggregated matrix via
#' \link{aggregate_values_by_key}.
#' If the user supplies a matrix to be aggregated the `association_file`
#' parameter must not be NULL: aggregation will be done by an internal
#' call to the aggregation function.
#' If the user supplies an already aggregated matrix, the `key` parameter
#' is the key used for aggregation -
#' **NOTE: for this operation is mandatory
#' that the time point column is included in the key.**
#' ## Assumptions on time point format
#' By using the functions provided by this package, when imported,
#' an association file will be correctly formatted for future usage.
#' In the formatting process there is also a padding operation performed on
#' time points: this means the functions expects the time point column to
#' be of type character and to be correctly padded with 0s. If the
#' chosen column for time point is detected as numeric the function will
#' attempt the conversion to character and automatic padding.
#' If you choose to import the association file not using the
#' \link{import_association_file} function, be sure to check the format of
#' the chosen column to avoid undesired results.
#'
#' @param x A simple integration matrix or an aggregated matrix (see details)
#' @param association_file NULL or the association file for x if `aggregate`
#' is set to TRUE
#' @param timepoint_column What is the name of the time point column?
#' @param key The aggregation key - must always contain the `timepoint_column`
#' @param include_tp_zero Include timepoint 0?
#' @param zero How is 0 coded in the data frame?
#' @param aggregate Should x be aggregated?
#' @param ... Additional parameters to pass to `aggregate_values_by_key`
#'
#' @family Analysis functions
#'
#' @importFrom dplyr mutate filter across all_of select summarise group_by
#' @importFrom dplyr arrange group_split first full_join starts_with distinct
#' @importFrom dplyr semi_join n rename
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data abort inform `:=`
#' @importFrom stringr str_pad
#' @importFrom purrr reduce is_empty
#' @importFrom tidyr pivot_longer
#' @importFrom stats na.omit
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
#' cumulative_count <- cumulative_count_union(aggreg)
#' cumulative_count
cumulative_count_union <- function(x,
    association_file = NULL,
    timepoint_column = "TimePoint",
    key = c(
        "SubjectID",
        "CellMarker",
        "Tissue",
        "TimePoint"
    ),
    include_tp_zero = FALSE,
    zero = "0000",
    aggregate = FALSE,
    ...) {
    stopifnot(is.data.frame(x))
    stopifnot(is.data.frame(association_file) | is.null(association_file))
    stopifnot(is.character(timepoint_column) & length(timepoint_column) == 1)
    stopifnot(is.character(key))
    stopifnot(is.logical(include_tp_zero))
    stopifnot(is.character(zero) & length(zero) == 1)
    stopifnot(is.logical(aggregate))
    if (aggregate == TRUE & is.null(association_file)) {
        rlang::abort(.agg_with_null_meta_err())
    }
    if (!all(timepoint_column %in% key)) {
        rlang::abort(.key_without_tp_err())
    }
    if (aggregate == FALSE) {
        if (!all(key %in% colnames(x))) {
            rlang::abort(.key_not_found())
        }
    } else {
        x <- aggregate_values_by_key(
            x = x,
            association_file = association_file,
            key = key, ...
        )
        if (is.numeric(association_file[[timepoint_column]]) |
            is.integer(association_file[[timepoint_column]])) {
            max <- max(association_file[[timepoint_column]])
            digits <- floor(log10(x)) + 1
            association_file <- association_file %>%
                dplyr::mutate(
                    {{ timepoint_column }} := stringr::str_pad(
                        as.character(.data$TimePoint),
                        digits,
                        side = "left",
                        pad = "0"
                    )
                )
            zero <- paste0(rep_len("0", digits), collapse = "")
        }
    }
    if (include_tp_zero == FALSE) {
        x <- x %>%
            dplyr::filter(dplyr::across(
                dplyr::all_of(timepoint_column),
                ~ .x != zero
            ))
        if (nrow(x) == 0) {
            all_tp0_msg <- paste(
                "All time points zeros were excluded, the data",
                "frame is empty."
            )
            rlang::inform(all_tp0_msg)
            return(x)
        }
    }
    annot <- if (.is_annotated(x)) {
        annotation_IS_vars()
    } else {
        character(0)
    }
    x <- x %>% dplyr::select(dplyr::all_of(c(
        mandatory_IS_vars(),
        annot,
        key
    )))
    cols_for_join <- colnames(x)[!colnames(x) %in% timepoint_column]
    key_minus_tp <- key[!key %in% timepoint_column]
    distinct_tp_for_each <- x %>%
        dplyr::arrange(dplyr::across(dplyr::all_of(timepoint_column))) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(key_minus_tp))) %>%
        dplyr::summarise(
            Distinct_tp = unique(
                `$`(.data, !!timepoint_column)
            ),
            .groups = "drop"
        ) %>%
        dplyr::rename({{ timepoint_column }} := "Distinct_tp")

    splitted <- x %>%
        dplyr::arrange(dplyr::across(dplyr::all_of(timepoint_column))) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(timepoint_column))) %>%
        dplyr::group_split()

    mult_tp <- purrr::reduce(splitted, dplyr::full_join, by = cols_for_join)
    tp_indexes <- grep(timepoint_column, colnames(mult_tp))
    tp_indexes <- tp_indexes[-1]
    res <- if (!purrr::is_empty(tp_indexes)) {
        for (i in tp_indexes) {
            mod_col <- mult_tp[[i]]
            mod_col_prec <- mult_tp[[i - 1]]
            val <- dplyr::first(stats::na.omit(mod_col))
            mod_col[!is.na(mod_col_prec)] <- val
            mult_tp[i] <- mod_col
        }
        mult_tp %>%
            tidyr::pivot_longer(
                cols = dplyr::starts_with(timepoint_column),
                values_to = timepoint_column,
                values_drop_na = TRUE
            ) %>%
            dplyr::select(-c("name")) %>%
            dplyr::distinct()
    } else {
        mult_tp
    }
    res <- res %>% dplyr::semi_join(distinct_tp_for_each, by = key)
    res <- res %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(key))) %>%
        dplyr::summarise(count = dplyr::n(), .groups = "drop")
    return(res)
}

#' Expands integration matrix with the cumulative is union over time.
#'
#' @description \lifecycle{experimental}
#' Given an input integration matrix that can be grouped over time,
#' this function adds integrations in groups assuming that
#' if an integration is observed at time point "t" then it is also observed in
#' time point "t+1".
#'
#' @param x An integration matrix, ideally aggregated via
#' `aggregate_values_by_key()`
#' @param key The aggregation key used
#' @param timepoint_col The name of the time point column
#' @param include_tp_zero Should time point 0 be included?
#' @param keep_og_is Keep original set of integrations as a separate column?
#' @param expand If `FALSE`, for each group, the set of integration sites is
#' returned in a separate column as a nested table, otherwise the resulting
#' column is unnested.
#'
#' @family Analysis functions
#' @return A data frame
#' @export
#'
#' @importFrom rlang .data
#' @importFrom data.table .SD
#'
#' @examples
#' data("integration_matrices", package = "ISAnalytics")
#' data("association_file", package = "ISAnalytics")
#' aggreg <- aggregate_values_by_key(
#'     x = rlang::current_env()$integration_matrices,
#'     association_file = rlang::current_env()$association_file,
#'     value_cols = c("seqCount", "fragmentEstimate")
#' )
#' cumulated_is <- cumulative_is(aggreg)
#' cumulated_is
cumulative_is <- function(x,
    key = c(
        "SubjectID",
        "CellMarker",
        "Tissue",
        "TimePoint"
    ),
    timepoint_col = "TimePoint",
    include_tp_zero = FALSE,
    keep_og_is = TRUE,
    expand = FALSE) {
    stopifnot(is.data.frame(x))
    stopifnot(is.character(key))
    stopifnot(is.character(timepoint_col))
    timepoint_col <- timepoint_col[1]
    stopifnot(is.logical(include_tp_zero))
    include_tp_zero <- include_tp_zero[1]
    stopifnot(is.logical(keep_og_is))
    stopifnot(is.logical(expand))
    if (!timepoint_col %in% key) {
        rlang::abort(.key_without_tp_err())
    }
    if (!all(key %in% colnames(x))) {
        rlang::abort(.missing_user_cols_error(key[!key %in% colnames(x)]))
    }
    is_vars <- if (.is_annotated(x)) {
        c(mandatory_IS_vars(), annotation_IS_vars())
    } else {
        mandatory_IS_vars()
    }
    temp <- x %>%
        dplyr::select(dplyr::all_of(c(key, is_vars))) %>%
        dplyr::mutate(!!timepoint_col := as.numeric(.data[[timepoint_col]]))
    if (!include_tp_zero) {
        temp <- temp %>%
            dplyr::filter(.data[[timepoint_col]] != 0)
        if (nrow(temp) == 0) {
            rlang::inform(.only_zero_tp())
            return(NULL)
        }
    }
    temp <- temp %>%
        dplyr::group_by(dplyr::across({{ key }})) %>%
        dplyr::arrange(.data[[timepoint_col]], .by_group = TRUE) %>%
        dplyr::distinct(dplyr::across(dplyr::all_of(is_vars)),
            .keep_all = TRUE
        )
    temp <- data.table::setDT(temp)
    temp <- temp[, .(is = list(.SD)), by = key]
    no_tp_key <- key[key != timepoint_col]
    splitted <- split(temp, by = no_tp_key)
    cumulate <- purrr::map(splitted, function(x) {
        x[, cumulative_is := purrr::accumulate(
            is,
            ~ data.table::funion(.x, .y)
        )]
    })
    cumulate <- data.table::rbindlist(cumulate)
    if (!keep_og_is) {
        cumulate[, is := NULL]
    }
    if (expand) {
        cumulate <- tidyr::unnest(cumulate,
            cols = "cumulative_is"
        )
        cumulate <- data.table::setDT(cumulate)
    }
    cumulate
}

#' Sharing of integration sites between given groups.
#'
#' \lifecycle{experimental}
#' Computes the amount of integration sites shared between the groups identified
#' in the input data.
#'
#' @details
#' An integration site is always identified by the triple
#' `(chr, integration_locus, strand)`, thus these columns must be present
#' in the input(s).
#'
#' The function accepts multiple inputs for different scenarios, please refer
#' to the vignette
#' \code{vignette("sharing_analyses", package = "ISAnalytics")}
#' for a more in-depth explanation.
#'
#' ## Output
#' The function outputs a single data frame containing all requested
#' comparisons and optionally individual group counts, genomic coordinates
#' of the shared integration sites and truth tables for plotting venn diagrams.
#'
#' ## Plotting sharing
#' The sharing data obtained can be easily plotted in a heatmap via the
#' function \code{\link{sharing_heatmap}} or via the function
#' \code{\link{sharing_venn}}
#'
#' @param ... One or more integration matrices
#' @param group_key Character vector of column names which identify a
#' single group. An associated group id will be derived by concatenating
#' the values of these fields, separated by "_"
#' @param group_keys A list of keys for asymmetric grouping.
#' If not NULL the argument `group_key` is ignored
#' @param n_comp Number of comparisons to compute. This argument is relevant
#' only if provided a single data frame and a single key.
#' @param is_count Logical, if `TRUE` returns also the count of IS for
#' each group and the count for the union set
#' @param relative_is_sharing Logical, if `TRUE` also returns the relative
#' sharing.
#' @param minimal Compute only combinations instead of all possible
#' permutations? If `TRUE` saves time and excludes redundant comparisons.
#' @param include_self_comp Include comparisons with the same group?
#' @param keep_genomic_coord If `TRUE` keeps the genomic coordinates of the
#' shared integration sites in a dedicated column (as a nested table)
#' @param table_for_venn Add column with truth tables for venn plots?
#'
#' @family Analysis functions
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
#' sharing <- is_sharing(aggreg)
#' sharing
is_sharing <- function(...,
    group_key = c(
        "SubjectID",
        "CellMarker",
        "Tissue",
        "TimePoint"
    ),
    group_keys = NULL,
    n_comp = 2,
    is_count = TRUE,
    relative_is_sharing = TRUE,
    minimal = TRUE,
    include_self_comp = FALSE,
    keep_genomic_coord = FALSE,
    table_for_venn = FALSE) {
    ## Checks
    if (!requireNamespace("gtools", quietly = TRUE)) {
        rlang::abort(.missing_pkg_error("gtools"))
    }
    dots <- rlang::list2(...)
    if (is.null(dots) || purrr::is_empty(dots)) {
        rlang::abort(.no_data_supp())
    }
    all_dfs <- purrr::map_lgl(dots, ~ is.data.frame(.x))
    if (!all(all_dfs)) {
        rlang::abort(.non_df_input_err())
    }
    stopifnot(is.null(group_keys) || is.list(group_keys))
    stopifnot(is.null(group_key) || is.character(group_key))
    stopifnot(is.logical(minimal))
    stopifnot(is.logical(keep_genomic_coord))
    stopifnot(is.logical(table_for_venn))
    key_mode <- if (!is.null(group_keys)) {
        if (any(purrr::map_lgl(
            group_keys,
            ~ all(is.character(.x))
        ) == FALSE)) {
            rlang::abort(.keys_not_char_err())
        }
        if (length(unique(group_keys)) == 1) {
            if (getOption("ISAnalytics.verbose") == TRUE) {
                one_key_list <- c("Single key in list",
                    i = paste(
                        "Provided a single key in list,",
                        "automatically performing",
                        "group comparisons"
                    )
                )
                rlang::inform(one_key_list, class = "one_key_list")
            }
            group_key <- group_keys[[1]]
            "SINGLE_KEY"
        } else {
            if (is.null(names(group_keys))) {
                rlang::inform(.unnamed_keys_warn())
                def_keys <- paste0("g", seq_along(group_keys))
                names(group_keys) <- def_keys
            }
            "MULT_KEY"
        }
    } else {
        if (!is.character(group_key)) {
            rlang::abort(.keys_not_char_err())
        }
        "SINGLE_KEY"
    }
    if (key_mode == "SINGLE_KEY") {
        stopifnot(is.logical(include_self_comp))
    }
    df_mode <- if (length(dots) == 1) {
        "SINGLE_DF"
    } else {
        "MULT_DF"
    }
    stopifnot(is.logical(is_count))
    stopifnot(is.logical(relative_is_sharing))
    if (df_mode == "SINGLE_DF") {
        ## Single dataframe provided
        if (!all(mandatory_IS_vars() %in% colnames(dots[[1]]))) {
            rlang::abort(
                .missing_mand_vars()
            )
        }
        if (key_mode == "SINGLE_KEY") {
            ## Single df - Single key
            if (!all(group_key %in% colnames(dots[[1]]))) {
                rlang::abort(
                    .missing_user_cols_error(
                        group_key[!group_key %in% colnames(dots[[1]])]
                    )
                )
            }
            stopifnot(is.numeric(n_comp) || is.integer(n_comp))
            n_comp <- n_comp[1]
            if (n_comp < 2) {
                rlang::abort("`n_comp` must be at least 2")
            }
        } else {
            ## Single df - multiple keys
            all_cols <- unique(unlist(group_keys))
            if (!all(all_cols %in% colnames(dots[[1]]))) {
                rlang::abort(
                    .missing_user_cols_error(
                        all_cols[!all_cols %in% colnames(dots[[1]])]
                    )
                )
            }
        }
    } else {
        all_mand_vars <- purrr::map_lgl(
            dots,
            ~ all(mandatory_IS_vars() %in%
                colnames(.x))
        )
        if (!all(all_mand_vars)) {
            missing_mand_at <- c("Missing mandatory vars in data frames",
                i = paste(
                    "At positions: ",
                    paste0(which(!all_mand_vars),
                        collapse = ", "
                    )
                )
            )
            rlang::abort(missing_mand_at)
        }
        if (key_mode == "SINGLE_KEY") {
            ## Multiple df - single key
            key_found_df <- purrr::map_lgl(
                dots,
                ~ all(group_key %in% colnames(.x))
            )
            if (!all(key_found_df)) {
                err_msg_key_not_found <- paste(
                    "Key not found in data frames",
                    paste0(which(!key_found_df),
                        collapse = ", "
                    )
                )
                rlang::abort(err_msg_key_not_found)
            }
        } else {
            ## Multiple df - multiple keys
            if (length(dots) != length(group_keys)) {
                keys_length_err <- c("Wrong key length",
                    i = paste(
                        "When providing multiple",
                        "input data frames,",
                        "`group_keys` must have",
                        "the same length"
                    )
                )
                rlang::abort(keys_length_err)
            }
            keys_ok <- purrr::map2_lgl(
                dots, group_keys,
                ~ all(.y %in% colnames(.x))
            )
            if (!all(keys_ok)) {
                mult_key_err <- c("Some keys not found in corresponding df",
                    x = paste(
                        "Issues identified at positions:",
                        paste0(which(!keys_ok),
                            collapse = ", "
                        )
                    )
                )
                rlang::abort(mult_key_err)
            }
        }
    }
    sharing <- if (key_mode == "SINGLE_KEY" & df_mode == "SINGLE_DF") {
        .sharing_singledf_single_key(
            df = dots[[1]],
            key = group_key,
            minimal = minimal,
            n_comp = n_comp,
            is_count = is_count,
            rel_sharing = relative_is_sharing,
            include_self_comp = include_self_comp,
            keep_genomic_coord = keep_genomic_coord,
            venn = table_for_venn
        )
    } else if (key_mode == "SINGLE_KEY" & df_mode == "MULT_DF") {
        .sharing_multdf_single_key(
            dfs = dots, key = group_key,
            minimal = minimal, is_count = is_count,
            rel_sharing = relative_is_sharing,
            keep_genomic_coord = keep_genomic_coord,
            venn = table_for_venn
        )
    } else if (key_mode == "MULT_KEY" & df_mode == "SINGLE_DF") {
        .sharing_singledf_mult_key(
            df = dots[[1]],
            keys = group_keys,
            minimal = minimal,
            is_count = is_count,
            rel_sharing = relative_is_sharing,
            keep_genomic_coord = keep_genomic_coord,
            venn = table_for_venn
        )
    } else {
        .sharing_multdf_mult_key(
            dfs = dots, keys = group_keys,
            minimal = minimal,
            is_count = is_count,
            rel_sharing = relative_is_sharing,
            keep_genomic_coord = keep_genomic_coord,
            venn = table_for_venn
        )
    }
    if (getOption("ISAnalytics.verbose") == TRUE) {
        rlang::inform("Done!")
    }
    return(sharing)
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
        shannon = ~ vegan::diversity(.x, index = "shannon"),
        simpson = ~ vegan::diversity(.x, index = "simpson"),
        invsimpson = ~ vegan::diversity(.x, index = "invsimpson"),
        sum = ~ sum(.x, na.rm = TRUE),
        count = length,
        describe = ~ tibble::as_tibble(psych::describe(.x))
    )
}
