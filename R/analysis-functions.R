#------------------------------------------------------------------------------#
# Analysis functions
#------------------------------------------------------------------------------#

#' Computes the abundance for every integration event in the input data frame.
#'
#' @description
#' `r lifecycle::badge("stable")`
#' Abundance is obtained for every integration event by calculating the ratio
#' between the single value and the total value for the given group.
#'
#' @details Abundance will be computed upon the user selected columns
#' in the `columns` parameter. For each column a corresponding
#' relative abundance column (and optionally a percentage abundance
#' column) will be produced.
#'
#' @section Required tags:
#' The function will explicitly check for the presence of these tags:
#'
#' * All columns declared in `mandatory_IS_vars()`
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
    stopifnot(is.logical(percentage))
    percentage <- percentage[1]
    if (!all(c(columns, key) %in% colnames(x))) {
        missing_cols <- c(
            columns[!columns %in% colnames(x)],
            key[!key %in% colnames(x)]
        )
        rlang::abort(.missing_user_cols_error(missing_cols))
    }
    non_num_cols <- purrr::map_lgl(
        columns,
        ~ is.numeric(x[[.x]]) || is.integer(x[[.x]])
    )
    if (any(!non_num_cols)) {
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

#' Filter data frames with custom predicates
#'
#' @description
#' `r lifecycle::badge("stable")`
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
#' @family Data cleaning and pre-processing
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
#' @description
#' `r lifecycle::badge("stable")`
#' The input data frame will be sorted by the highest values in
#' the columns specified and the top n rows will be returned as output.
#' The user can choose to keep additional columns in the output
#' by passing a vector of column names or passing 2 "shortcuts":
#' * `keep = "everything"` keeps all columns in the original data frame
#' * `keep = "nothing"` only keeps the mandatory columns
#' (`mandatory_IS_vars()`) plus the columns in the `columns` parameter.
#'
#' @section Required tags:
#' The function will explicitly check for the presence of these tags:
#'
#' * All columns declared in `mandatory_IS_vars()`
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
# top_abundant_is
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


#' Top n targeted genes based on number of IS.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#' Produces a summary of the number of integration events per gene, orders
#' the table in decreasing order and slices the first n rows - either on
#' all the data frame or by group.
#'
#' @details
#' ## Gene grouping
#' When producing a summary of IS by gene, there are different options that
#' can be chosen.
#' The argument `consider_chr` accounts for the fact that some genes (same
#' gene symbol) may span more than one chromosome: if set to `TRUE`
#' counts of IS will be separated for those genes that span 2 or more
#' chromosomes - in other words they will be in 2 different rows of the
#' output table. On the contrary, if the argument is set to `FALSE`,
#' counts will be produced in a single row.
#'
#' NOTE: the function counts **DISTINCT** integration events, which logically
#' corresponds to a union of sets. Be aware of the fact that counts per group
#' and counts with different arguments might be different: if for example
#' counts are performed by considering chromosome and there is one gene symbol
#' with 2 different counts, the sum of those 2 will likely not be equal to
#' the count obtained by performing the calculations without
#' considering the chromosome.
#'
#' The same reasoning can be applied for the argument `consider_gene_strand`,
#' that takes into account the strand of the gene.
#'
#' @section Required tags:
#' The function will explicitly check for the presence of these tags:
#'
#' ```{r echo=FALSE, results="asis"}
#' all_tags <- available_tags()
#' needed <- unique(all_tags[purrr::map_lgl(eval(rlang::sym("needed_in")),
#'  ~ "top_targeted_genes" %in% .x)][["tag"]])
#'  cat(paste0("* ", needed, collapse="\n"))
#' ```
#'
#' Note that the tags "gene_strand" and "chromosome" are explicitly required
#' only if `consider_chr = TRUE` and/or `consider_gene_strand = TRUE`.
#'
#' @param x An integration matrix - must be annotated
#' @param n Number of rows to slice
#' @param key If slice has to be performed for each group, the character
#' vector of column names that identify the groups. If `NULL` considers the
#' whole input data frame.
#' @param consider_chr Logical, should the chromosome be taken into account?
#' See details.
#' @param consider_gene_strand Logical, should the gene strand be taken into
#' account? See details.
#' @param as_df If computation is performed by group, `TRUE` returns all
#' groups merged in a single data frame with a column containing the group id.
#' If `FALSE` returns a named list.
#'
#' @importFrom rlang sym
#'
#' @return A data frame or a list of data frames
#' @export
#' @family Analysis functions
#'
#' @examples
#' data("integration_matrices", package = "ISAnalytics")
#' top_targ <- top_targeted_genes(
#'     integration_matrices,
#'     key = NULL
#' )
#' top_targ
top_targeted_genes <- function(x,
    n = 20,
    key = c(
        "SubjectID", "CellMarker",
        "Tissue", "TimePoint"
    ),
    consider_chr = TRUE,
    consider_gene_strand = TRUE,
    as_df = TRUE) {
    stopifnot(is.data.frame(x))
    data.table::setDT(x)
    stopifnot(is.numeric(n) || is.integer(n))
    stopifnot(is.null(key) || is.character(key))
    stopifnot(is.logical(consider_chr))
    stopifnot(is.logical(consider_gene_strand))

    required_annot_tags <- c("gene_symbol")
    if (consider_gene_strand) {
        required_annot_tags <- c(required_annot_tags, "gene_strand")
    }
    annot_tag_cols <- .check_required_cols(
        required_annot_tags,
        annotation_IS_vars(TRUE),
        "error"
    )
    if (consider_chr) {
        chr_tag_col <- .check_required_cols(
            c("chromosome", "locus"),
            mandatory_IS_vars(TRUE),
            "error"
        )
        annot_tag_cols <- annot_tag_cols %>%
            dplyr::bind_rows(chr_tag_col)
    }
    data.table::setDT(annot_tag_cols)
    cols_to_check <- c(annot_tag_cols$names, key)
    if (!all(cols_to_check %in% colnames(x))) {
        rlang::abort(.missing_needed_cols(
            cols_to_check[!cols_to_check %in% colnames(x)]
        ))
    }

    df_with_is_counts <- if (is.null(key)) {
        .count_distinct_is_per_gene(
            x = x, include_chr = consider_chr,
            include_gene_strand = consider_gene_strand,
            gene_sym_col = annot_tag_cols[
                eval(sym("tag")) == "gene_symbol"
            ][["names"]],
            gene_strand_col = annot_tag_cols[
                eval(sym("tag")) == "gene_strand"
            ][["names"]],
            chr_col = annot_tag_cols[eval(sym("tag")) ==
                "chromosome"][["names"]],
            mand_vars_to_check = mandatory_IS_vars(TRUE)
        ) %>%
            dplyr::arrange(dplyr::desc(.data$n_IS)) %>%
            dplyr::slice_head(n = n)
    } else {
        tmp <- x[, .count_distinct_is_per_gene(
            x = .SD, include_chr = consider_chr,
            include_gene_strand = consider_gene_strand,
            gene_sym_col = annot_tag_cols[
                eval(sym("tag")) == "gene_symbol"
            ][["names"]],
            gene_strand_col = annot_tag_cols[
                eval(sym("tag")) == "gene_strand"
            ][["names"]],
            chr_col = annot_tag_cols[eval(sym("tag")) ==
                "chromosome"][["names"]],
            mand_vars_to_check = mandatory_IS_vars(TRUE)
        ), by = eval(key)]
        tmp[,
            .SD %>% dplyr::arrange(dplyr::desc(.data$n_IS)) %>%
                dplyr::slice_head(n = n),
            by = eval(key)
        ]
    }
    if (as_df) {
        return(df_with_is_counts)
    }
    return(split(df_with_is_counts, by = key))
}


#' Compute Fisher's exact test on gene frequencies.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#' Provided 2 data frames with calculations for CIS, via `CIS_grubbs()`,
#' computes Fisher's exact test.
#' Results can be plotted via `fisher_scatterplot()`.
#'
#' @param cis_x A data frame obtained via `CIS_grubbs()`
#' @param cis_y A data frame obtained via `CIS_grubbs()`
#' @param min_is_per_gene Used for pre-filtering purposes. Genes with a
#' number of distinct integration less than this number will be filtered out
#' prior calculations. Single numeric or integer.
#' @param gene_set_method One between "intersection" and "union". When merging
#' the 2 data frames, `intersection` will perform an inner join operation,
#' while `union` will perform a full join operation.
#' @param significance_threshold Significance threshold for the Fisher's
#' test p-value
#' @param remove_unbalanced_0 Remove from the final output those pairs in
#' which there are no IS for one group or the other and the number of
#' IS of the non-missing group are less than the mean number of IS for that
#' group
#'
#' @template genes_db
#'
#' @section Required tags:
#' The function will explicitly check for the presence of these tags:
#'
#' ```{r echo=FALSE, results="asis"}
#' all_tags <- available_tags()
#' needed <- unique(all_tags[purrr::map_lgl(eval(rlang::sym("needed_in")),
#'  ~ "gene_frequency_fisher" %in% .x)][["tag"]])
#'  cat(paste0("* ", needed, collapse="\n"))
#' ```
#'
#' @return A data frame
#' @export
#' @family Analysis functions
#'
#' @examples
#' data("integration_matrices", package = "ISAnalytics")
#' data("association_file", package = "ISAnalytics")
#' aggreg <- aggregate_values_by_key(
#'     x = integration_matrices,
#'     association_file = association_file,
#'     value_cols = c("seqCount", "fragmentEstimate")
#' )
#' cis <- CIS_grubbs(aggreg, by = "SubjectID")
#' fisher <- gene_frequency_fisher(cis$cis$PT001, cis$cis$PT002,
#'     min_is_per_gene = 2
#' )
#' fisher
gene_frequency_fisher <- function(cis_x,
    cis_y,
    min_is_per_gene = 3,
    gene_set_method = c("intersection", "union"),
    onco_db_file = "proto_oncogenes",
    tumor_suppressors_db_file = "tumor_suppressors",
    species = "human",
    known_onco = known_clinical_oncogenes(),
    suspicious_genes =
        clinical_relevant_suspicious_genes(),
    significance_threshold = 0.05,
    remove_unbalanced_0 = TRUE) {
    ## --- Input checks
    stopifnot(is.data.frame(cis_x) && is.data.frame(cis_y))
    stopifnot(is.integer(min_is_per_gene) || is.numeric(min_is_per_gene))
    gene_set_method <- rlang::arg_match(gene_set_method)
    stopifnot(is.character(onco_db_file))
    stopifnot(is.character(tumor_suppressors_db_file))
    stopifnot(is.character(species))
    stopifnot(is.data.frame(known_onco))
    stopifnot(is.data.frame(suspicious_genes))
    stopifnot(is.numeric(significance_threshold))
    stopifnot(is.logical(remove_unbalanced_0))
    ## -- Fetch gene symbol column
    gene_sym_col <- .check_required_cols(
        "gene_symbol", annotation_IS_vars(TRUE),
        duplicate_politic = "error"
    )[["names"]]
    req_cis_cols <- c(
        gene_sym_col, "n_IS_perGene", "average_TxLen",
        "raw_gene_integration_frequency"
    )
    quiet_expand <- purrr::quietly(.expand_cis_df)
    cols_for_join <- c(
        gene_sym_col,
        "Onco1_TS2", "ClinicalRelevance", "DOIReference",
        "KnownGeneClass", "KnownClonalExpansion",
        "CriticalForInsMut"
    )
    ## --- Calculations to perform on each df
    append_calc <- function(df, group_n) {
        if (!all(req_cis_cols %in% colnames(df))) {
            rlang::abort(
                .missing_needed_cols(req_cis_cols[!req_cis_cols %in%
                    colnames(df)])
            )
        }
        modified <- quiet_expand(
            df, gene_sym_col,
            onco_db_file, tumor_suppressors_db_file,
            species, known_onco, suspicious_genes
        )$result
        modified <- modified %>%
            dplyr::mutate(
                IS_per_kbGeneLen = .data$raw_gene_integration_frequency * 1000,
                Sum_IS_per_kbGeneLen = sum(.data$IS_per_kbGeneLen,
                    na.rm = TRUE
                ),
                IS_per_kbGeneLen_perMDepth_TPM = (.data$IS_per_kbGeneLen /
                    .data$Sum_IS_per_kbGeneLen) * 1e6
            ) %>%
            dplyr::filter(.data$n_IS_perGene >= min_is_per_gene) %>%
            dplyr::select(dplyr::all_of(c(
                req_cis_cols, cols_for_join,
                "IS_per_kbGeneLen",
                "Sum_IS_per_kbGeneLen",
                "IS_per_kbGeneLen_perMDepth_TPM"
            )))
        colnames(modified)[!colnames(modified) %in% cols_for_join] <- paste(
            colnames(modified)[!colnames(modified) %in% cols_for_join], group_n,
            sep = "_"
        )
        return(modified)
    }
    cis_mod <- purrr::map2(list(cis_x, cis_y), c(1, 2), append_calc)
    ## --- Merge the two in 1 df
    merged <- if (gene_set_method == "union") {
        purrr::reduce(cis_mod, ~ dplyr::full_join(.x, .y, by = cols_for_join))
    } else {
        purrr::reduce(cis_mod, ~ dplyr::inner_join(.x, .y, by = cols_for_join))
    }
    if (nrow(merged) == 0) {
        if (getOption("ISAnalytics.verbose", TRUE) == TRUE) {
            msg <- c("Data frame empty after filtering",
                i = paste(
                    "Data frame is empty after applying filter on IS,",
                    "is your filter too stringent?"
                ),
                x = "Nothing to return"
            )
            rlang::inform(msg, class = "empty_df_gene_freq")
        }
        return(NULL)
    }
    ## --- Actual computation of fisher test: test is applied on each row
    ## (each gene)
    merged <- merged %>%
        dplyr::mutate(
            tot_n_IS_perGene_1 = sum(cis_x$n_IS_perGene, na.rm = TRUE),
            tot_n_IS_perGene_2 = sum(cis_y$n_IS_perGene, na.rm = TRUE)
        )
    compute_fisher <- function(...) {
        row <- list(...)
        n_IS_perGene_1 <- row$n_IS_perGene_1
        n_IS_perGene_2 <- row$n_IS_perGene_2
        n_IS_perGene_1[which(is.na(n_IS_perGene_1))] <- 0
        n_IS_perGene_2[which(is.na(n_IS_perGene_2))] <- 0
        matrix <- matrix(
            data = c(
                n_IS_perGene_1,
                row$tot_n_IS_perGene_1 - n_IS_perGene_1,
                n_IS_perGene_2,
                row$tot_n_IS_perGene_2 - n_IS_perGene_2
            ),
            nrow = 2,
            dimnames = list(
                G1 = c("IS_of_gene", "TotalIS"),
                G2 = c("IS_of_gene", "TotalIS")
            )
        )
        ft <- stats::fisher.test(matrix)
        return(ft$p.value)
    }
    merged <- merged %>%
        dplyr::mutate(
            Fisher_p_value = purrr::pmap_dbl(., compute_fisher)
        ) %>%
        dplyr::mutate(
            Fisher_p_value_significant = dplyr::if_else(
                condition = .data$Fisher_p_value < significance_threshold,
                true = TRUE, false = FALSE
            )
        )
    ## --- Removing unbalanced 0s if requested - this scenario applies
    ## only if "union" is selected as method for join
    if (remove_unbalanced_0) {
        mean_is_per_gene_1 <- ceiling(mean(merged$n_IS_perGene_1, na.rm = TRUE))
        mean_is_per_gene_2 <- ceiling(mean(merged$n_IS_perGene_2, na.rm = TRUE))
        test_exclude <- function(...) {
            row <- list(...)
            if (is.na(row$n_IS_perGene_1) || is.na(row$n_IS_perGene_2)) {
                to_ex <- ifelse(
                    test = ((row$n_IS_perGene_1 < mean_is_per_gene_1) &
                        (is.na(row$n_IS_perGene_2))) |
                        ((is.na(row$n_IS_perGene_1)) &
                            (row$n_IS_perGene_2 < mean_is_per_gene_2)),
                    yes = TRUE,
                    no = FALSE
                )
                return(to_ex)
            }
            return(FALSE)
        }
        merged <- merged %>%
            dplyr::mutate(
                to_exclude_from_test = purrr::pmap(., test_exclude)
            ) %>%
            dplyr::filter(.data$to_exclude_from_test == FALSE) %>%
            dplyr::select(-dplyr::all_of("to_exclude_from_test"))
        if (nrow(merged) == 0) {
            if (getOption("ISAnalytics.verbose", TRUE) == TRUE) {
                msg <- c("Data frame empty after filtering",
                    i = paste(
                        "Data frame is after removing unbalanced IS,",
                        "nothing to return"
                    )
                )
                rlang::inform(msg, class = "empty_df_gene_freq_unbal")
            }
            return(NULL)
        }
    }
    ## --- Apply statistical corrections to p-value
    merged <- merged %>%
        dplyr::mutate(
            Fisher_p_value_fdr = stats::p.adjust(.data$Fisher_p_value,
                method = "fdr",
                n = length(.data$Fisher_p_value)
            ),
            Fisher_p_value_benjamini = stats::p.adjust(.data$Fisher_p_value,
                method = "BY",
                n = length(.data$Fisher_p_value)
            ),
            Fisher_p_value_bonferroni = stats::p.adjust(.data$Fisher_p_value,
                method = "bonferroni",
                n = length(.data$Fisher_p_value)
            ),
            minus_log10_pvalue = -log(.data$Fisher_p_value, base = 10)
        ) %>%
        dplyr::mutate(
            minus_log10_pvalue_fdr = -log(.data$Fisher_p_value_fdr, base = 10),
        )
    return(merged)
}


#' Computes user specified functions on numerical columns and updates
#' the metadata data frame accordingly.
#'
#' @description
#' `r lifecycle::badge("stable")`
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
#' `mandatory_IS_vars()`
#'
#' @section Required tags:
#' The function will explicitly check for the presence of these tags:
#'
#' * All columns declared in `mandatory_IS_vars()`
#'
#' These are checked only if `add_integrations_count = TRUE`.
#'
#' @family Analysis functions
#' @importFrom rlang .data sym
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
sample_statistics <- function(x,
    metadata,
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
        rlang::abort(.non_num_user_cols_error(
            value_columns[!vcols_are_numeric]
        ))
    }
    purrr::walk(functions, function(f) {
        if (!(purrr::is_function(f) | purrr::is_formula(f))) {
            rlang::abort(.non_function_elem_error())
        }
    })

    result <- x %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(sample_key))) %>%
        dplyr::summarise(
            dplyr::across(
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
            if (getOption("ISAnalytics.verbose", TRUE)) {
                rlang::inform(.inform_skip_count_is())
            }
        }
    }

    updated_meta <- metadata %>% dplyr::left_join(result, by = sample_key)
    return(list(x = result, metadata = updated_meta))
}


#' Grubbs test for Common Insertion Sites (CIS).
#'
#' @description
#' `r lifecycle::badge("stable")`
#' Statistical approach for the validation of common insertion sites
#' significance based on the comparison of the integration frequency
#' at the CIS gene with respect to other genes contained in the
#' surrounding genomic regions. For more details please refer to
#' this paper:
#' <`r .lentiviral_CIS_paper()`>
#'
#' @details
#' ## Genomic annotation file
#' A data frame containing
#' genes annotation for the specific genome.
#' From version `1.5.4` the argument `genomic_annotation_file` accepts only
#' data frames or package provided defaults.
#' The user is responsible for importing the appropriate tabular files if
#' customization is needed.
#' The annotations for the human genome (hg19) and
#' murine genome (mm9) are already
#' included in this package: to use one of them just
#' set the argument `genomic_annotation_file` to either `"hg19"` or
#' `"mm9"`.
#' If for any reason the user is performing an analysis on another genome,
#' this file needs to be changed respecting the USCS Genome Browser
#' format, meaning the input file headers should include:
#'
#' `r refGene_table_cols()`
#'
#' @section Required tags:
#' The function will explicitly check for the presence of these tags:
#'
#' ```{r echo=FALSE, results="asis"}
#' all_tags <- available_tags()
#' needed <- unique(all_tags[purrr::map_lgl(eval(rlang::sym("needed_in")),
#'  ~ "CIS_grubbs" %in% .x)][["tag"]])
#'  cat(paste0("* ", needed, collapse="\n"))
#' ```
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
#' ("SubjectID" column must be included in the input data frame).
#' @param return_missing_as_df Returns those genes present in the input df
#' but not in the refgenes as a data frame?
#' @param results_as_list If `TRUE`
#' return the group computations as a named list, otherwise return a single
#' df with an additional column containing the group id
#'
#' @family Analysis functions
#'
#' @importFrom rlang .data sym
#'
#' @return A data frame
#' @export
#'
#' @examples
#' data("integration_matrices", package = "ISAnalytics")
#' cis <- CIS_grubbs(integration_matrices)
#' cis
CIS_grubbs <- function(x,
    genomic_annotation_file = "hg19",
    grubbs_flanking_gene_bp = 100000,
    threshold_alpha = 0.05,
    by = NULL,
    return_missing_as_df = TRUE,
    results_as_list = TRUE) {
    res_checks <- .cis_param_check(
        x, genomic_annotation_file,
        grubbs_flanking_gene_bp,
        threshold_alpha,
        return_missing_as_df
    )
    stopifnot(is.null(by) || is.character(by))
    if (!all(by %in% colnames(x))) {
        rlang::abort(.missing_user_cols_error(by[!by %in% colnames(x)]))
    }
    result <- list()
    join_ref_res <- .cis_join_ref(x, res_checks)
    result <- append(result, join_ref_res[
        names(join_ref_res) %in% c("missing_genes", "missing_is")
    ])
    if (is.null(by)) {
        cis <- .cis_grubb_calc(
            x = join_ref_res$joint_ref,
            grubbs_flanking_gene_bp = res_checks$grubbs_flanking_gene_bp,
            threshold_alpha = res_checks$threshold_alpha,
            gene_symbol_col = res_checks$gene_symbol_col,
            gene_strand_col = res_checks$gene_strand_col,
            chr_col = res_checks$chrom_col, locus_col = res_checks$locus_col,
            strand_col = res_checks$strand_col
        )
    } else {
        grouped <- join_ref_res$joint_ref %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(by)))
        group_ks <- grouped %>%
            dplyr::group_keys() %>%
            tidyr::unite(col = "id", dplyr::everything()) %>%
            dplyr::pull(.data$id)
        split <- grouped %>%
            dplyr::group_split() %>%
            purrr::set_names(group_ks)
        cis <- purrr::map(split, ~ .cis_grubb_calc(
            x = .x,
            grubbs_flanking_gene_bp = res_checks$grubbs_flanking_gene_bp,
            threshold_alpha = res_checks$threshold_alpha,
            gene_symbol_col = res_checks$gene_symbol_col,
            gene_strand_col = res_checks$gene_strand_col,
            chr_col = res_checks$chrom_col, locus_col = res_checks$locus_col,
            strand_col = res_checks$strand_col
        ))
        if (!results_as_list) {
            cis <- purrr::map2(cis, names(cis), ~ {
                .x %>%
                    dplyr::mutate(group = .y)
            }) %>% purrr::reduce(dplyr::bind_rows)
        }
    }
    result$cis <- cis
    return(result)
}


#' Compute CIS and Grubbs test over different time points and groups.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#' Computes common insertion sites and Grubbs test for each separate group
#' and separating different time points among the same group. The logic
#' applied is the same as the function `CIS_grubbs()`.
#'
#' @inherit CIS_grubbs details
#'
#' @inheritParams CIS_grubbs
#' @param group A character vector of column names that identifies a group.
#' Each group must contain one or more time points.
#' @param timepoint_col What is the name of the column containing time points?
#' @param as_df Choose the result format: if `TRUE` the results are returned
#' as a single data frame containing a column for the group id and a column
#' for the time point, if `FALSE` results are returned in the form of nested
#' lists (one table for each time point and for each group), if `"group"`
#' results are returned as a list separated for each group but containing a
#' single table with all time points.
#' @param max_workers Maximum number of parallel workers. If `NULL` the
#' maximum number of workers is calculated automatically.
#'
#'
#' @return A list with results and optionally missing genes info
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
#' cis_overtime <- CIS_grubbs_overtime(aggreg)
#' cis_overtime
CIS_grubbs_overtime <- function(x,
    genomic_annotation_file = "hg19",
    grubbs_flanking_gene_bp = 100000,
    threshold_alpha = 0.05,
    group = "SubjectID",
    timepoint_col = "TimePoint",
    as_df = TRUE,
    return_missing_as_df = TRUE,
    max_workers = NULL) {
    result <- list()
    res_checks <- .cis_param_check(
        x, genomic_annotation_file,
        grubbs_flanking_gene_bp, threshold_alpha,
        return_missing_as_df
    )
    stopifnot(is.character(group))
    stopifnot(is.character(timepoint_col))
    stopifnot(is.null(max_workers) || is.numeric(max_workers))
    timepoint_col <- timepoint_col[1]
    if (!all(c(group, timepoint_col) %in% colnames(x))) {
        rlang::abort(.missing_user_cols_error(
            c(group, timepoint_col)[!c(group, timepoint_col) %in% colnames(x)]
        ))
    }
    stopifnot(is.logical(as_df) || is.character(as_df))
    result_as_df <- if (is.logical(as_df)) {
        if (as_df[1] == TRUE) {
            1
        } else {
            0
        }
    } else if (is.character(as_df) & as_df[1] == "group") {
        2
    } else {
        err_as_df <- c("The argument `as_df` must be one of the allowed values",
            x = paste(
                "Arg can be either logical (T/F) or",
                "equal to 'group'"
            ),
            i = paste("See `?CIS_grubbs_overtime`")
        )
        rlang::abort(err_as_df)
    }
    join_ref_res <- .cis_join_ref(x, res_checks)
    result <- append(result, join_ref_res[
        names(join_ref_res) %in% c("missing_genes", "missing_is")
    ])
    # --- Split according to group
    split <- join_ref_res$joint_ref %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(group)))
    keys_g <- split %>%
        dplyr::group_keys() %>%
        tidyr::unite(col = "id") %>%
        dplyr::pull(.data$id)
    split <- split %>%
        dplyr::group_split() %>%
        purrr::set_names(keys_g)

    # --- Calculate for each group and each tp
    tp_slice_cis <- function(df, timepoint_col,
    res_checks, result_as_df,
    progress) {
        tmp <- df %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(timepoint_col)))
        g_keys <- tmp %>%
            dplyr::group_keys() %>%
            dplyr::pull(!!timepoint_col)
        tmp <- tmp %>%
            dplyr::group_split() %>%
            purrr::set_names(g_keys)
        res <- if (result_as_df == 0) {
            purrr::map(tmp, ~ .cis_grubb_calc(
                x = .x,
                grubbs_flanking_gene_bp = res_checks$grubbs_flanking_gene_bp,
                threshold_alpha = res_checks$threshold_alpha,
                gene_symbol_col = res_checks$gene_symbol_col,
                gene_strand_col = res_checks$gene_strand_col,
                chr_col = res_checks$chrom_col,
                locus_col = res_checks$locus_col,
                strand_col = res_checks$strand_col
            ))
        } else {
            purrr::map2_df(tmp, names(tmp), ~ .cis_grubb_calc(
                x = .x,
                grubbs_flanking_gene_bp = res_checks$grubbs_flanking_gene_bp,
                threshold_alpha = res_checks$threshold_alpha,
                gene_symbol_col = res_checks$gene_symbol_col,
                gene_strand_col = res_checks$gene_strand_col,
                chr_col = res_checks$chrom_col,
                locus_col = res_checks$locus_col,
                strand_col = res_checks$strand_col
            ) %>% dplyr::mutate(!!timepoint_col := .y))
        }
        if (!is.null(progress)) {
            progress()
        }
        return(res)
    }
    cis_overtime <- .execute_map_job(
        data_list = split,
        fun_to_apply = tp_slice_cis,
        fun_args = list(
            timepoint_col = timepoint_col,
            result_as_df = result_as_df,
            res_checks = res_checks
        ),
        stop_on_error = TRUE,
        max_workers = max_workers
    )
    if (result_as_df == 1) {
        cis_overtime <- purrr::map2_df(
            cis_overtime$res, names(cis_overtime$res),
            ~ .x %>%
                dplyr::mutate(group = .y)
        )
    } else {
        cis_overtime <- cis_overtime$res
    }
    result$cis <- cis_overtime
    return(result)
}

#' Integrations cumulative count in time by sample
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' This function was deprecated in favour of a single function,
#' please use `cumulative_is` instead.
#'
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
#' @return A data frame
#' @export
#' @keywords internal
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
    lifecycle::deprecate_warn(
        when = "1.5.4",
        what = "cumulative_count_union()",
        with = "cumulative_is()",
        details = c(paste(
            "Use option `counts = TRUE`.",
            "Function will be likely dropped in the",
            "next release cycle"
        ))
    )
    cumulative_is(
        x = x, timepoint_col = timepoint_column, key = key,
        include_tp_zero = include_tp_zero, counts = TRUE
    )
}

#' Expands integration matrix with the cumulative IS union over time.
#'
#' @description
#' `r lifecycle::badge("experimental")`
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
#' @param counts Add cumulative counts? Logical
#'
#' @section Required tags:
#' The function will explicitly check for the presence of these tags:
#'
#' * All columns declared in `mandatory_IS_vars()`
#' * Checks if the matrix is annotated by assessing presence of
#' `annotation_IS_vars()`
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
#'     x = integration_matrices,
#'     association_file = association_file,
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
    counts = TRUE,
    keep_og_is = FALSE,
    expand = TRUE) {
    stopifnot(is.data.frame(x))
    stopifnot(is.character(key))
    stopifnot(is.character(timepoint_col))
    timepoint_col <- timepoint_col[1]
    stopifnot(is.logical(include_tp_zero))
    include_tp_zero <- include_tp_zero[1]
    stopifnot(is.logical(counts))
    counts <- counts[1]
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
            rlang::inform(.only_zero_tp(), class = "only_zero_tps")
            return(NULL)
        }
    }
    temp <- temp %>%
        dplyr::group_by(dplyr::across({{ key }})) %>%
        dplyr::arrange(.data[[timepoint_col]], .by_group = TRUE) %>%
        dplyr::distinct(dplyr::across(dplyr::all_of(is_vars)),
            .keep_all = TRUE
        )
    data.table::setDT(temp)
    temp <- temp[, list(is = list(.SD)), by = key]
    no_tp_key <- key[key != timepoint_col]
    split <- split(temp, by = no_tp_key)
    cumulate <- purrr::map(split, function(x) {
        x[, cumulative_is := purrr::accumulate(
            get("is"),
            ~ data.table::funion(.x, .y)
        )]
    })
    cumulate <- data.table::rbindlist(cumulate)
    if (!keep_og_is) {
        cumulate[, is := NULL]
    }
    counts_df <- if (counts) {
        cumulate[, list(is_n_cumulative = unlist(purrr::map(
            get("cumulative_is"), nrow
        ))), by = key]
    } else {
        NULL
    }
    if (expand) {
        cumulate <- tidyr::unnest(cumulate,
            cols = "cumulative_is"
        )
        data.table::setDT(cumulate)
    }
    to_return <- if (counts) {
        list(coordinates = cumulate, counts = counts_df)
    } else {
        cumulate
    }
    return(to_return)
}

#' Sharing of integration sites between given groups.
#'
#' @description
#' `r lifecycle::badge("stable")`
#' Computes the amount of integration sites shared between the groups identified
#' in the input data.
#'
#' @details
#' An integration site is always identified by the combination of fields in
#' `mandatory_IS_vars()`, thus these columns must be present
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
#' @section Required tags:
#' The function will explicitly check for the presence of these tags:
#'
#' * All columns declared in `mandatory_IS_vars()`
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
            if (getOption("ISAnalytics.verbose", TRUE) == TRUE) {
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
    if (getOption("ISAnalytics.verbose", TRUE) == TRUE) {
        rlang::inform("Done!")
    }
    return(sharing)
}


#' Find the source of IS by evaluating sharing.
#'
#' @description
#' `r lifecycle::badge("stable")`
#' The function computes the sharing between a reference group of interest
#' for each time point and a selection of groups of interest. In this way
#' it is possible to observe the percentage of shared integration sites between
#' reference and each group and identify in which time point a certain IS was
#' observed for the first time.
#'
#' @param reference A data frame containing one or more groups of reference.
#' Groups are identified by `ref_group_key`
#' @param selection A data frame containing one or more groups of interest
#' to compare.
#' Groups are identified by `selection_group_key`
#' @param ref_group_key Character vector of column names that identify a
#' unique group in the `reference` data frame
#' @param selection_group_key Character vector of column names that identify a
#' unique group in the `selection` data frame
#' @param timepoint_column Name of the column holding time point
#' info?
#' @param by_subject Should calculations be performed for each subject
#' separately?
#' @param subject_column Name of the column holding subjects information.
#' Relevant only if `by_subject = TRUE`
#'
#' @return A list of data frames or a data frame
#' @family Analysis functions
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
#' df1 <- aggreg %>%
#'     dplyr::filter(.data$Tissue == "BM")
#' df2 <- aggreg %>%
#'     dplyr::filter(.data$Tissue == "PB")
#' source <- iss_source(df1, df2)
#' source
#' ggplot2::ggplot(source$PT001, ggplot2::aes(
#'     x = as.factor(g2_TimePoint),
#'     y = sharing_perc, fill = g1
#' )) +
#'     ggplot2::geom_col() +
#'     ggplot2::labs(
#'         x = "Time point", y = "Shared IS % with MNC BM",
#'         title = "Source of is MNC BM vs MNC PB"
#'     )
iss_source <- function(reference,
    selection,
    ref_group_key = c(
        "SubjectID", "CellMarker",
        "Tissue", "TimePoint"
    ),
    selection_group_key = c(
        "SubjectID", "CellMarker",
        "Tissue", "TimePoint"
    ),
    timepoint_column = "TimePoint",
    by_subject = TRUE,
    subject_column = "SubjectID") {
    ## Checks
    stopifnot(is.data.frame(reference) & is.data.frame(selection))
    stopifnot(is.character(ref_group_key) & is.character(selection_group_key))
    stopifnot(is.character(timepoint_column))
    stopifnot(is.logical(by_subject))
    by_subject <- by_subject[1]
    if (!all(ref_group_key %in% colnames(reference))) {
        rlang::abort(.missing_user_cols_error(
            ref_group_key[!ref_group_key %in% colnames(reference)]
        ))
    }
    if (!all(selection_group_key %in% colnames(selection))) {
        rlang::abort(.missing_user_cols_error(
            selection_group_key[!selection_group_key %in% colnames(selection)]
        ))
    }
    timepoint_column <- timepoint_column[1]
    if (!timepoint_column %in% colnames(reference) |
        !timepoint_column %in% colnames(selection)) {
        rlang::abort(.missing_needed_cols(timepoint_column))
    }
    ## Workflow choice
    if (by_subject) {
        stopifnot(is.character(subject_column))
        subject_column <- subject_column[1]
        ref_split <- reference %>%
            dplyr::group_by(dplyr::across({{ subject_column }}))
        ref_subjs <- ref_split %>%
            dplyr::group_keys() %>%
            dplyr::pull(.data[[subject_column]])
        ref_split <- ref_split %>%
            dplyr::group_split() %>%
            purrr::set_names(ref_subjs)
        sel_split <- selection %>%
            dplyr::group_by(dplyr::across({{ subject_column }}))
        sel_subjs <- sel_split %>%
            dplyr::group_keys() %>%
            dplyr::pull(.data[[subject_column]])
        sel_split <- sel_split %>%
            dplyr::group_split() %>%
            purrr::set_names(sel_subjs)
        shared <- .sharing_for_source(ref_split,
            sel_split,
            ref_key = ref_group_key,
            sel_key = selection_group_key,
            tp_col = timepoint_column,
            subj_col = subject_column
        )
        shared <- purrr::map(
            shared,
            ~ .x %>%
                dplyr::select(
                    -dplyr::all_of(c("count_g1", "count_g2", "count_union"))
                ) %>%
                tidyr::unnest(dplyr::all_of("is_coord"), keep_empty = TRUE) %>%
                dplyr::mutate(sharing_perc = dplyr::if_else(
                    shared == 0, 0, .data$on_g2 / .data$shared
                )) %>%
                dplyr::select(
                    -dplyr::all_of(c("shared", "on_g1", "on_g2", "on_union"))
                )
        )
    } else {
        shared <- .sharing_for_source(reference,
            selection,
            ref_key = ref_group_key,
            sel_key = selection_group_key,
            tp_col = timepoint_column,
            subj_col = subject_column
        )
        shared <- shared %>%
            dplyr::select(
                -dplyr::all_of(c("count_g1", "count_g2", "count_union"))
            ) %>%
            tidyr::unnest(dplyr::all_of("is_coord"), keep_empty = TRUE) %>%
            dplyr::mutate(sharing_perc = dplyr::if_else(shared == 0,
                0,
                .data$on_g2 /
                    .data$shared
            )) %>%
            dplyr::select(
                -dplyr::all_of(c("shared", "on_g1", "on_g2", "on_union"))
            )
    }

    return(shared)
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
