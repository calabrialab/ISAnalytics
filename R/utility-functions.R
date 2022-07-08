#------------------------------------------------------------------------------#
# Utility functions
#------------------------------------------------------------------------------#

#' Retrieve description of a tag by name.
#'
#' @description
#' Given one or multiple tags, prints the associated description
#' and functions where the tag is explicitly used.
#'
#' @param tags A character vector of tag names
#'
#' @importFrom rlang sym
#'
#' @return `NULL`
#' @export
#' @family dynamic vars
#'
#' @examples
#' inspect_tags(c("chromosome", "project_id", "x"))
inspect_tags <- function(tags) {
    all_tags <- available_tags()
    desc_msg <- purrr::map(tags, ~ {
        tag_desc <- all_tags[eval(sym("tag")) == .x][["description"]]
        if (purrr::is_empty(tag_desc)) {
            return(
                c(paste("* TAG:", .x),
                    x = "Tag not found in available tags"
                )
            )
        }
        functions <- all_tags[eval(sym("tag")) == .x][["needed_in"]][[1]]
        functions <- paste0(functions, collapse = ", ")
        c(paste("* TAG:", .x),
            i = paste("Description:", tag_desc),
            i = paste("Functions that use it:", functions)
        )
    })
    purrr::walk(desc_msg, ~ rlang::inform(.x, class = "tag_inspect"))
}


#' Define custom dynamic vars.
#'
#' @description
#' This set of function allows users to specify custom look-up tables for
#' dynamic variables.
#' For more details, refer to the dedicated vignette
#' `vignette("workflow_start", package="ISAnalytics")`.
#'
#' * `set_mandatory_IS_vars()` sets the look-up table for mandatory IS vars.
#'
#'
#' @template vars_lookuptable_str
#'
#' @param specs Either a named vector or a data frame with specific format.
#' See details.
#'
#' @return `NULL`
#'
#' @export
#' @aliases vars_setters
#' @family dynamic vars
#'
#' @examples
#' tmp_mand_vars <- tibble::tribble(
#'     ~names, ~types, ~transform, ~flag, ~tag,
#'     "chrom", "char", ~ stringr::str_replace_all(.x, "chr", ""), "required",
#'     "chromosome",
#'     "position", "int", NULL, "required", "locus",
#'     "strand", "char", NULL, "required", "is_strand",
#'     "gap", "int", NULL, "required", NA_character_,
#'     "junction", "int", NULL, "required", NA_character_
#' )
#' set_mandatory_IS_vars(tmp_mand_vars)
#' print(mandatory_IS_vars(TRUE))
#' reset_mandatory_IS_vars()
#'
set_mandatory_IS_vars <- function(specs) {
    # Check structure
    err <- c("Input format error",
        i = paste(
            "The input must follow the specification",
            "in the documentation, see ?set_mandatory_IS_vars",
            "for details"
        )
    )
    new_vars <- .new_IS_vars_checks(specs, err, "mand_vars")
    options(ISAnalytics.mandatory_is_vars = new_vars)
    if (getOption("ISAnalytics.verbose", TRUE)) {
        rlang::inform("Mandatory IS vars successfully changed")
    }
}

#' Resets dynamic vars to the default values.
#'
#' @description
#' Reverts all changes to dynamic vars to the default values.
#' For more details, refer to the dedicated vignette
#' `vignette("workflow_start", package="ISAnalytics")`.
#'
#' * `reset_mandatory_IS_vars()` re-sets the look-up table for
#' mandatory IS vars.
#'
#' @return `NULL`
#' @export
#' @aliases vars_resetters
#' @family dynamic vars
#'
#' @examples
#' reset_mandatory_IS_vars()
#'
reset_mandatory_IS_vars <- function() {
    options(ISAnalytics.mandatory_is_vars = "default")
    if (getOption("ISAnalytics.verbose", TRUE)) {
        rlang::inform("Mandatory IS vars reset to default")
    }
}

#' Current dynamic vars specifications getters.
#'
#' @description Fetches the look-up tables for different categories of dynamic
#' vars.
#' For more details, refer to the dedicated vignette
#' `vignette("workflow_start", package="ISAnalytics")`.
#'
#' * `mandatory_IS_vars` returns the look-up table of variables that are
#' used to uniquely identify integration events
#'
#' @param include_types If set to `TRUE` returns both the names and the types
#' associated, otherwise returns only a character vector of names
#'
#' @return A character vector or a data frame
#' @export
#' @aliases vars_getters
#' @family dynamic vars
#'
#' @examples
#' # Names only
#' mandatory_IS_vars()
#'
#' # Names and types
#' mandatory_IS_vars(TRUE)
#'
mandatory_IS_vars <- function(include_types = FALSE) {
    opt <- getOption("ISAnalytics.mandatory_is_vars", "default")
    if (!include_types) {
        if (all(opt == "default")) {
            default_vars <- .default_mandatory_IS_vars()
            return(default_vars$names)
        }
        return(opt$names)
    }
    if (all(opt == "default")) {
        default_vars <- .default_mandatory_IS_vars()
        return(default_vars)
    }
    return(opt)
}



#' @rdname set_mandatory_IS_vars
#' @description
#' * `set_annotation_IS_vars()` sets the look-up table for genomic annotation
#' IS vars.
#' @export
#'
#' @examples
#' tmp_annot_vars <- tibble::tribble(
#'     ~names, ~types, ~transform, ~flag, ~tag,
#'     "gene", "char", NULL, "required",
#'     "gene_symbol",
#'     "gene_strand", "char", NULL, "required", "gene_strand"
#' )
#' print(annotation_IS_vars(TRUE))
#' reset_annotation_IS_vars()
#'
set_annotation_IS_vars <- function(specs) {
    # Check structure
    err <- c("Input format error",
        i = paste(
            "The input must follow the specification",
            "in the documentation, see ?set_annotation_IS_vars",
            "for details"
        )
    )
    new_vars <- .new_IS_vars_checks(specs, err, "annot_vars")
    options(ISAnalytics.genomic_annotation_vars = new_vars)
    if (getOption("ISAnalytics.verbose", TRUE)) {
        rlang::inform("Annotation IS vars successfully changed")
    }
}

#' @rdname reset_mandatory_IS_vars
#' @description
#' * `reset_annotation_IS_vars()` re-sets the look-up table for
#' genomic annotation IS vars.
#' @export
#' @examples
#' reset_annotation_IS_vars()
#'
reset_annotation_IS_vars <- function() {
    options(ISAnalytics.genomic_annotation_vars = "default")
    if (getOption("ISAnalytics.verbose", TRUE)) {
        rlang::inform("Annotation IS vars reset to default")
    }
}

#' @rdname mandatory_IS_vars
#' @description
#' * `annotation_IS_vars()` returns the look-up table of variables that
#' contain genomic annotations
#' @export
#'
#' @examples
#' # Names only
#' annotation_IS_vars()
#'
#' # Names and types
#' annotation_IS_vars(TRUE)
#'
annotation_IS_vars <- function(include_types = FALSE) {
    opt <- getOption("ISAnalytics.genomic_annotation_vars", "default")
    if (!include_types) {
        if (all(opt == "default")) {
            default_vars <- .default_annotation_IS_vars()
            return(default_vars$names)
        }
        return(opt$names)
    }
    if (all(opt == "default")) {
        default_vars <- .default_annotation_IS_vars()
        return(default_vars)
    }
    return(opt)
}


#' @rdname set_mandatory_IS_vars
#' @description
#' * `set_af_columns_def()` sets the look-up table for association file columns
#' vars
#' @export
#'
#' @examples
#' temp_af_cols <- tibble::tribble(
#'     ~names, ~types, ~transform, ~flag, ~tag,
#'     "project", "char", NULL, "required",
#'     "project_id",
#'     "pcr_id", "char", NULL, "required", "pcr_repl_id",
#'     "subject", "char", NULL, "required", "subject"
#' )
#' set_af_columns_def(temp_af_cols)
#' print(association_file_columns(TRUE))
#' reset_af_columns_def()
#'
set_af_columns_def <- function(specs) {
    # Check structure
    err <- c("Input format error",
        i = paste(
            "The input must follow the specification",
            "in the documentation, see ?set_af_columns_def",
            "for details"
        )
    )
    new_vars <- .new_IS_vars_checks(specs, err, "af_vars")
    options(ISAnalytics.af_specs = new_vars)
    if (getOption("ISAnalytics.verbose", TRUE)) {
        rlang::inform("Association file columns specs successfully changed")
    }
}

#' @rdname reset_mandatory_IS_vars
#' @export
#' @description
#' * `reset_af_columns_def()` re-sets the look-up table for
#' association file columns vars
#'
#' @examples
#' reset_af_columns_def()
#'
reset_af_columns_def <- function() {
    options(ISAnalytics.af_specs = "default")
    if (getOption("ISAnalytics.verbose", TRUE)) {
        rlang::inform("Association file columns specs reset to default")
    }
}

#' @rdname mandatory_IS_vars
#' @description
#' * `association_file_columns()` returns the look-up table of variables that
#' contains information on how  metadata is structured
#'
#' @export
#' @examples
#' # Names only
#' association_file_columns()
#'
#' # Names and types
#' association_file_columns(TRUE)
#'
association_file_columns <- function(include_types = FALSE) {
    opt <- getOption("ISAnalytics.af_specs", "default")
    if (!include_types) {
        if (all(opt == "default")) {
            default_vars <- .default_af_cols()
            return(default_vars$names)
        }
        return(opt$names)
    }
    if (all(opt == "default")) {
        default_vars <- .default_af_cols()
        return(default_vars)
    }
    return(opt)
}

#' @rdname set_mandatory_IS_vars
#' @description
#' * `set_iss_stats_specs()` sets the look-up table for VISPA2 pool statistics
#' vars
#' @export
#'
#' @examples
#' tmp_iss_vars <- tibble::tribble(
#'     ~names, ~types, ~transform, ~flag, ~tag,
#'     "pool", "char", NULL, "required",
#'     "vispa_concatenate",
#'     "tag", "char", NULL, "required", "tag_seq",
#'     "barcode", "int", NULL, "required", NA_character_
#' )
#' set_iss_stats_specs(tmp_iss_vars)
#' iss_stats_specs(TRUE)
#' reset_iss_stats_specs()
set_iss_stats_specs <- function(specs) {
    # Check structure
    err <- c("Input format error",
        i = paste(
            "The input must follow the specification",
            "in the documentation, see ?set_iss_stats_specs",
            "for details"
        )
    )
    new_vars <- .new_IS_vars_checks(specs, err, "iss_vars")
    options(ISAnalytics.iss_stats_specs = new_vars)
    if (getOption("ISAnalytics.verbose", TRUE)) {
        rlang::inform("ISS stats specs successfully changed")
    }
}

#' @rdname reset_mandatory_IS_vars
#' @description
#' * `reset_iss_stats_specs()` re-sets the look-up table for VISPA2 pool
#' statistics vars
#' @export
#'
#' @examples
#' reset_iss_stats_specs()
#'
reset_iss_stats_specs <- function() {
    options(ISAnalytics.iss_stats_specs = "default")
    if (getOption("ISAnalytics.verbose", TRUE)) {
        rlang::inform("ISS stats specs reset to default")
    }
}

#' @rdname mandatory_IS_vars
#' @description
#' * `iss_stats_specs()` returns the look-up table of variables that
#' contains information on the format of pool statistics files produced
#' automatically by VISPA2
#'
#' @export
#'
#' @examples
#' # Names only
#' iss_stats_specs()
#'
#' # Names and types
#' iss_stats_specs(TRUE)
#'
iss_stats_specs <- function(include_types = FALSE) {
    opt <- getOption("ISAnalytics.iss_stats_specs", "default")
    if (!include_types) {
        if (all(opt == "default")) {
            default_vars <- .default_iss_stats_specs()
            return(default_vars$names)
        }
        return(opt$names)
    }
    if (all(opt == "default")) {
        default_vars <- .default_iss_stats_specs()
        return(default_vars)
    }
    return(opt)
}


#' @rdname mandatory_IS_vars
#' @description
#' * `matrix_file_suffixes()` returns the look-up table of variables that
#' contains all default file names for each quantification type and it is
#' used by automated import functions
#'
#' @export
#'
#' @examples
#' # Names only
#' matrix_file_suffixes()
#'
matrix_file_suffixes <- function() {
    opt <- getOption("ISAnalytics.matrix_file_suffix", "default")

    if (all(opt == "default")) {
        defaults_params <- as.list(formals(set_matrix_file_suffixes))
        default_vars <- do.call(.generate_suffix_specs, args = defaults_params)
        return(default_vars)
    }

    return(opt)
}

#' Sets the look-up table for matrix file suffixes.
#'
#' @description
#' The function automatically produces and sets a look-up table of matrix file
#' suffixes based on user input.
#'
#' @param quantification_suffix A named list - names must be quantification
#' types in `quantification_types()`, and values must be single strings,
#' containing the associated suffix. Please note that ALL quantification
#' types must be specified or the function will produce an error.
#' @param annotation_suffix A named list - names must be `annotated` and
#' `not_annotated`, values must be single strings,
#' containing the associated suffix. Please note that both names must be
#' present in the list or the function will produce an error.
#' @param file_ext The file extension (e.g. `tsv`, `tsv.gz`)
#' @param glue_file_spec A string specifying the pattern used to form the
#' entire suffix, as per \code{\link[glue:glue]{glue::glue()}} requirements.
#' The string should contain the reference to `quantification_suffix`,
#' `annotation_suffix` and `file_ext`.
#'
#' @family dynamic vars
#' @return `NULL`
#' @export
#'
#' @examples
#' set_matrix_file_suffixes(
#'     quantification_suffix = list(
#'         seqCount = "sc",
#'         fragmentEstimate = "fe",
#'         barcodeCount = "barcodeCount",
#'         cellCount = "cellCount",
#'         ShsCount = "ShsCount"
#'     ),
#'     annotation_suffix = list(annotated = "annot", not_annotated = "")
#' )
#' matrix_file_suffixes()
#' reset_matrix_file_suffixes()
set_matrix_file_suffixes <- function(quantification_suffix = list(
        seqCount = "seqCount",
        fragmentEstimate = "fragmentEstimate",
        barcodeCount = "barcodeCount",
        cellCount = "cellCount",
        ShsCount = "ShsCount"
    ),
    annotation_suffix = list(
        annotated = ".no0.annotated",
        not_annotated = ""
    ),
    file_ext = "tsv.gz",
    glue_file_spec = "{quantification_suffix}_matrix{annotation_suffix}.{file_ext}") {
    stopifnot(is.list(quantification_suffix) &&
        !is.null(names(quantification_suffix)))
    stopifnot(is.list(annotation_suffix) &&
        !is.null(names(annotation_suffix)))
    stopifnot(is.character(file_ext))
    stopifnot(is.character(glue_file_spec))
    glue_file_spec <- glue_file_spec[1]
    if (!all(quantification_types() %in% names(quantification_suffix))) {
        warn_msg <- c("Warning: some quantification types are missing from specs",
            i = paste(
                "Missing: ",
                quantification_types()[!quantification_types() %in%
                    names(quantification_suffix)]
            ),
            i = paste(
                "Missing quantifications specs may cause problems",
                "when trying to import files.",
                "See documentation with `?set_matrix_file_suffixes`"
            )
        )
        rlang::inform(warn_msg, class = "missing_quant_specs")
    }
    quantification_suffix <- quantification_suffix[
        names(quantification_suffix) %in% quantification_types()
    ]
    if (!all(c("annotated", "not_annotated") %in% names(annotation_suffix))) {
        err_msg <- c(paste(
            "Error: a value for both 'annotated'",
            "and 'not_annotated' is requested. Quitting."
        ))
        rlang::abort(err_msg, class = "miss_annot_suff_specs")
    }

    final_specs <- .generate_suffix_specs(
        quantification_suffix,
        annotation_suffix,
        file_ext,
        glue_file_spec
    )

    options(ISAnalytics.matrix_file_suffix = final_specs)
    if (getOption("ISAnalytics.verbose", TRUE)) {
        rlang::inform("Matrix suffixes specs successfully changed")
    }
}

#' @rdname reset_mandatory_IS_vars
#' @description
#' * `reset_matrix_file_suffixes()` re-sets the matrix file suffixes look-up
#' table
#' @export
#'
#' @examples
#' reset_matrix_file_suffixes()
#'
reset_matrix_file_suffixes <- function() {
    options(ISAnalytics.matrix_file_suffix = "default")
    if (getOption("ISAnalytics.verbose", TRUE)) {
        rlang::inform("Matrix suffixes specs reset to default")
    }
}


#' @rdname reset_mandatory_IS_vars
#' @description
#' * `reset_dyn_vars_config()` re-sets all look-up tables
#' @export
#'
#' @examples
#' reset_dyn_vars_config()
reset_dyn_vars_config <- function() {
    reset_mandatory_IS_vars()
    reset_annotation_IS_vars()
    reset_af_columns_def()
    reset_iss_stats_specs()
    reset_matrix_file_suffixes()
}

#' Easily retrieve the name of the pcr id column.
#'
#' @description The function is a shortcut to retrieve the currently set
#' pcr id column name from the association file column tags look-up table.
#' This column is needed every time a joining operation with metadata
#' needs to be performed
#'
#' @return The name of the column
#' @family dynamic vars
#' @export
#'
#' @examples
#' pcr_id_column()
pcr_id_column <- function() {
    af_cols_specs <- association_file_columns(TRUE)
    selected_tags <- .check_required_cols(
        list("pcr_repl_id" = "char"),
        af_cols_specs, "error"
    )
    return(selected_tags$names)
}

#' Apply transformations to an arbitrary number of columns.
#'
#' @description
#' This function takes a named list of purr-style lambdas where names are the
#' names of the columns in the data frame that must be transformed.
#' NOTE: the columns are overridden, not appended.
#'
#' @details
#' Lambdas provided in input must be transformations, aka functions that take
#' in input a vector and return a vector of the same length as the input.
#'
#' If the input transformation list contains column names that are not present
#' in the input data frame, they are simply ignored.
#'
#' @param df The data frame on which transformations should be operated
#' @param transf_list A named list of purrr-style lambdas, where names
#' are column names the function should be applied to.
#'
#' @importFrom rlang sym
#'
#' @return A data frame with transformed columns
#' @export
#' @family Utilities
#'
#' @examples
#' df <- tibble::tribble(
#'     ~A, ~B, ~C, ~D,
#'     1, 2, "a", "aa",
#'     3, 4, "b", "bb",
#'     5, 6, "c", "cc"
#' )
#' lambdas <- list(A = ~ .x + 1, B = ~ .x + 2, C = ~ stringr::str_to_upper(.x))
#' transform_columns(df, lambdas)
transform_columns <- function(df, transf_list) {
    transf_list <- transf_list[names(transf_list) %in% colnames(df)]
    if (purrr::is_empty(transf_list)) {
        return(df)
    }
    transf_list_mod <- purrr::map2(names(transf_list), transf_list, ~ {
        rlang::expr(rlang::as_function(!!.y)(!!rlang::sym(.x)))
    }) %>%
        purrr::set_names(names(transf_list))
    withCallingHandlers(
        {
            withRestarts(
                {
                    return(dplyr::mutate(df, !!!transf_list_mod))
                },
                skip_transform = function() {
                    string_list <- paste(names(transf_list), transf_list,
                        sep = " = ", collapse = "; "
                    )
                    err_msg <- c(paste(
                        "Encountered a problem while applying transformations",
                        "to columns. Skipping this."
                    ),
                    i = paste(
                        "Check the correctness of your input,",
                        "see the documentation"
                    ),
                    i = paste("Your input: ", string_list)
                    )
                    if (getOption("ISAnalytics.verbose", TRUE)) {
                        rlang::inform(err_msg, class = "skip_col_transform")
                    }
                    return(df)
                }
            )
        },
        error = function(cnd) {
            res <- findRestart("skip_transform")
            invokeRestart(res)
        }
    )
}


#' Obtain a single integration matrix from individual quantification
#' matrices.
#'
#' @description
#' `r lifecycle::badge("stable")`
#' Takes a list of integration matrices referring to different quantification
#' types and merges them into a single data frame with multiple
#' value columns, each renamed according to their quantification type
#' of reference.
#'
#' @param x A named list of integration matrices, ideally obtained via
#' \link{import_parallel_Vispa2Matrices}. Names must be
#' quantification types in `quantification_types()`.
#' @param fragmentEstimate The name of the output column for fragment
#' estimate values
#' @param seqCount The name of the output column for sequence
#' count values
#' @param barcodeCount The name of the output column for barcode count
#' values
#' @param cellCount The name of the output column for cell count values
#' @param ShsCount The name of the output column for Shs count values
#' @param value_col_name Name of the column containing the corresponding
#' values in the single matrices
#'
#' @importFrom rlang .data `:=`
#'
#' @family Utilities
#'
#' @seealso \link{quantification_types}
#'
#' @return A single data frame
#' @export
#'
#' @examples
#' sc <- tibble::tribble(
#'     ~chr, ~integration_locus, ~strand, ~CompleteAmplificationID, ~Value,
#'     "1", 45324, "+", "ID1", 543,
#'     "2", 52423, "-", "ID1", 42,
#'     "6", 54623, "-", "ID2", 67,
#'     "X", 12314, "+", "ID3", 8
#' )
#' fe <- tibble::tribble(
#'     ~chr, ~integration_locus, ~strand, ~CompleteAmplificationID, ~Value,
#'     "1", 45324, "+", "ID1", 56.76,
#'     "2", 52423, "-", "ID1", 78.32,
#'     "6", 54623, "-", "ID2", 123.45,
#'     "X", 12314, "+", "ID3", 5.34
#' )
#' comparison_matrix(list(
#'     fragmentEstimate = fe,
#'     seqCount = sc
#' ))
comparison_matrix <- function(x,
    fragmentEstimate = "fragmentEstimate",
    seqCount = "seqCount",
    barcodeCount = "barcodeCount",
    cellCount = "cellCount",
    ShsCount = "ShsCount",
    value_col_name = "Value") {
    stopifnot(is.list(x) & !is.data.frame(x))
    stopifnot(all(names(x) %in% quantification_types()))
    stopifnot(is.character(fragmentEstimate) & length(fragmentEstimate) == 1)
    stopifnot(is.character(seqCount) & length(seqCount) == 1)
    stopifnot(is.character(barcodeCount) & length(barcodeCount) == 1)
    stopifnot(is.character(cellCount) & length(cellCount) == 1)
    stopifnot(is.character(ShsCount) & length(ShsCount) == 1)
    stopifnot(is.character(value_col_name) & length(value_col_name) == 1)
    param_names <- c(
        fragmentEstimate = fragmentEstimate,
        seqCount = seqCount, barcodeCount = barcodeCount,
        cellCount = cellCount, ShsCount = ShsCount
    )
    x <- purrr::map2(x, names(x), function(matrix, quant_type) {
        quant_name <- param_names[names(param_names) %in% quant_type]
        matrix %>% dplyr::rename(!!quant_name := .data[[value_col_name]])
    })
    result <- purrr::reduce(x, function(matrix1, matrix2) {
        commoncols <- dplyr::intersect(colnames(matrix1), colnames(matrix2))
        matrix1 %>%
            dplyr::full_join(matrix2, by = commoncols)
    })
    na_introduced <- purrr::map_lgl(param_names, function(p) {
        any(is.na(result[[p]]))
    })
    if (any(na_introduced) & getOption("ISAnalytics.verbose", TRUE) == TRUE) {
        rlang::inform(.nas_introduced_msg(), class = "comp_nas")
    }
    result
}


#' Separate a multiple-quantification matrix into single quantification
#' matrices.
#'
#' @description
#' `r lifecycle::badge("stable")`
#' The function separates a single multi-quantification integration
#' matrix, obtained via \link{comparison_matrix}, into single
#' quantification matrices as a named list of tibbles.
#'
#' @param x Single integration matrix with multiple quantification
#' value columns, obtained via \link{comparison_matrix}.
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
#' @family Utilities
#'
#' @return A named list of data frames, where names are quantification types
#' @seealso \link{quantification_types}
#' @export
#'
#' @examples
#' data("integration_matrices", package = "ISAnalytics")
#' separated <- separate_quant_matrices(
#'     integration_matrices
#' )
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
    if (!purrr::is_empty(to_copy) &
        getOption("ISAnalytics.verbose", TRUE) == TRUE) {
        rlang::inform(.non_quant_cols_msg(to_copy))
    }
    separated <- purrr::map(num_cols, function(quant) {
        x %>%
            dplyr::select(dplyr::all_of(c(key, to_copy, quant))) %>%
            dplyr::rename(Value = .data[[quant]])
    }) %>% purrr::set_names(names(num_cols))
    separated
}

#' Create a blank association file.
#'
#' @description
#' Produces a blank association file to start using
#' both VISPA2 and ISAnalytics
#'
#' @param path The path on disk where the file should be written - must be a
#' file
#'
#' @return `NULL`
#' @export
#' @family Utilities
#'
#' @examples
#' temp <- tempfile()
#' generate_blank_association_file(temp)
generate_blank_association_file <- function(path) {
    stopifnot(is.character(path))
    af <- data.frame(matrix(
        ncol = length(association_file_columns()),
        nrow = 0
    ))
    colnames(af) <- association_file_columns()
    path <- fs::as_fs_path(path)
    dir <- fs::path_dir(path)
    if (!fs::dir_exists(dir)) {
        fs::dir_create(dir, recurse = TRUE)
    }
    readr::write_tsv(af, file = path)
}

#' Creates a reduced association file for a VISPA2 run,
#' given project and pool
#'
#' @description
#' The function selects the appropriate columns and prepares a file for the
#' launch of VISPA2 pipeline for each project/pool pair specified.
#'
#' @details
#' Note: the function is vectorized, meaning you can specify more than
#' one project and more than one pool as vectors of characters, but you must
#' ensure that:
#'
#' * Both `project` and `pool` vectors have the same length
#' * You correclty type names in corresponding positions, for example
#' c("PJ01", "PJ01") - c("POOL01", "POOL02").
#' If you type a pool in the position of a corresponding
#' project that doesn't match no file will be produced since that pool doesn't
#' exist in the corresponding project.
#'
#' @section Required tags:
#' The function will explicitly check for the presence of these tags:
#'
#' ```{r echo=FALSE, results="asis"}
#' cat(paste0("* ", reduced_AF_columns()$tag, collapse = "\n"))
#' ```
#'
#' The names of the pools in the `pool` argument is checked against the
#' column corresponding to the `pool_id` tag.
#'
#' @param association_file The imported association file (via
#' `import_association_file()`)
#' @param project A vector of characters containing project names
#' @param pool A vector of characters containing pool names
#' @param path A single string representing the path to the folder where files
#' should be written. If the folder doesn't exist it will be created.
#'
#' @importFrom rlang sym
#'
#' @family Utilities
#'
#' @return `NULL`
#' @export
#'
#' @examples
#' temp <- tempdir()
#' data("association_file", package = "ISAnalytics")
#' generate_Vispa2_launch_AF(association_file, "PJ01", "POOL01", temp)
generate_Vispa2_launch_AF <- function(association_file,
    project,
    pool,
    path) {
    stopifnot(is.data.frame(association_file))
    stopifnot(is.character(project))
    stopifnot(is.character(pool))
    stopifnot(length(project) == length(pool))
    stopifnot(is.character(path) & length(path) == 1)
    path <- fs::as_fs_path(path)
    if (!fs::file_exists(path)) {
        fs::dir_create(path)
    }
    af_min_cols <- reduced_AF_columns()
    concat_col <- .check_required_cols("vispa_concatenate",
        vars_df = association_file_columns(TRUE),
        duplicate_politic = "error"
    )
    concat_col <- concat_col$names
    to_check <- c(af_min_cols$names, concat_col)
    if (!all(to_check %in% colnames(association_file))) {
        rlang::abort(.missing_af_needed_cols(to_check[!to_check %in%
            association_file_columns(TRUE)]))
    }
    proj_col <- af_min_cols[eval(sym("tag")) == "project_id", ][["names"]]
    pool_col <- af_min_cols[eval(sym("tag")) == "pool_id", ][["names"]]
    tag_id_col <- af_min_cols[eval(sym("tag")) == "tag_id", ][["names"]]
    tissue_col <- af_min_cols[eval(sym("tag")) == "tissue", ][["names"]]
    tp_col <- af_min_cols[eval(sym("tag")) == "tp_days", ][["names"]]
    fusion_col <- af_min_cols[eval(sym("tag")) == "fusion_id", ][["names"]]
    pcr_id_col <- af_min_cols[eval(sym("tag")) == "pcr_repl_id", ][["names"]]
    cm_col <- af_min_cols[eval(sym("tag")) == "cell_marker", ][["names"]]
    vec_id_col <- af_min_cols[eval(sym("tag")) == "vector_id", ][["names"]]
    process_proj_pool <- function(x, y) {
        selected_cols <- association_file %>%
            dplyr::select(dplyr::all_of(to_check)) %>%
            dplyr::filter(
                .data[[proj_col]] == x,
                .data[[pool_col]] == y
            ) %>%
            dplyr::mutate(TagID2 = .data[[tag_id_col]]) %>%
            dplyr::mutate(PoolName = dplyr::if_else(
                !is.na(.data[[concat_col]]),
                .data[[concat_col]],
                .data[[pool_col]]
            )) %>%
            dplyr::select(
                .data[[tag_id_col]], .data$TagID2, .data[[tissue_col]],
                .data[[tp_col]], .data[[fusion_col]], .data[[pcr_id_col]],
                .data[[cm_col]], .data[[proj_col]], .data[[vec_id_col]],
                .data[[pool_col]], .data$PoolName
            )
    }
    files <- purrr::map2(project, pool, process_proj_pool) %>%
        purrr::set_names(paste0(project, "-", pool, "_AF.tsv"))
    purrr::walk2(files, names(files), function(x, y) {
        complete_path <- fs::path(path, y)
        if (nrow(x) > 0) {
            readr::write_tsv(x, complete_path, col_names = FALSE)
        } else {
            msg <- c(paste("Nothing to write for ", y, ", skipping."))
            rlang::inform(msg, class = "launch_af_empty")
        }
    })
}


#' Converts tidy integration matrices in the original sparse matrix
#' form.
#'
#' @description
#' `r lifecycle::badge("stable")`
#' This function is particularly useful when a sparse matrix structure
#' is needed by a specific function (mainly from other packages).
#'
#' @param x A single tidy integration matrix or a list of integration
#' matrices. Supports also multi-quantification matrices
#' obtained via \link{comparison_matrix}
#' @param single_value_col Name of the column containing the values when
#' providing a single-quantification matrix
#' @param fragmentEstimate For multi-quantification matrix support:
#' the name of the fragment estimate values column
#' @param seqCount For multi-quantification matrix support:
#' the name of the sequence count values column
#' @param barcodeCount For multi-quantification matrix support:
#' the name of the barcode count values column
#' @param cellCount For multi-quantification matrix support:
#' the name of the cell count values column
#' @param ShsCount For multi-quantification matrix support:
#' the name of the Shs Count values column
#' @param key The name of the sample identifier fields (for aggregated
#' matrices can be a vector with more than 1 element)
#'
#' @family Utilities
#'
#' @return Depending on input, 2 possible outputs:
#' * A single sparse matrix (data frame) if input is a single quantification
#' matrix
#' * A list of sparse matrices divided by quantification if input
#' is a single multi-quantification matrix or a list of matrices
#' @export
#'
#' @examples
#' data("integration_matrices", package = "ISAnalytics")
#' sparse <- as_sparse_matrix(integration_matrices)
as_sparse_matrix <- function(x,
    single_value_col = "Value",
    fragmentEstimate = "fragmentEstimate",
    seqCount = "seqCount",
    barcodeCount = "barcodeCount",
    cellCount = "cellCount",
    ShsCount = "ShsCount",
    key = pcr_id_column()) {
    stopifnot(is.list(x))
    stopifnot(is.character(fragmentEstimate))
    stopifnot(is.character(seqCount))
    stopifnot(is.character(barcodeCount))
    stopifnot(is.character(cellCount))
    stopifnot(is.character(ShsCount))
    stopifnot(is.character(key))
    param_cols <- c(
        fragmentEstimate, seqCount, barcodeCount, cellCount,
        ShsCount, single_value_col
    )
    ## -- Internal for processing single df
    process_single_df <- function(df) {
        if (.check_mandatory_vars(df) == FALSE) {
            rlang::abort(.non_ISM_error())
        }
        if (!all(key %in% colnames(df))) {
            rlang::abort(.missing_needed_cols(key[!key %in% colnames(df)]))
        }
        id_cols <- c(mandatory_IS_vars())
        if (.is_annotated(df)) {
            id_cols <- c(id_cols, annotation_IS_vars())
        }
        id_cols <- c(id_cols, key)
        num_cols <- .find_exp_cols(df, id_cols)
        if (purrr::is_empty(num_cols)) {
            rlang::abort(.missing_num_cols_error())
        }
        pivot_dfs <- purrr::map(num_cols, ~ {
            if (.x %in% param_cols) {
                pivoted <- df %>%
                    dplyr::select(dplyr::all_of(
                        c(id_cols, .x)
                    )) %>%
                    tidyr::pivot_wider(
                        names_from = key,
                        values_from = .data[[.x]]
                    )
                return(pivoted)
            }
            return(NULL)
        }) %>%
            purrr::set_names(num_cols)
        pivot_dfs <- pivot_dfs[purrr::map_lgl(pivot_dfs, ~ !is.null(.x))]
        if (length(pivot_dfs) == 1) {
            pivot_dfs <- pivot_dfs[[1]]
        }
        return(pivot_dfs)
    }
    if (is.data.frame(x)) {
        return(process_single_df(x))
    } else {
        ## LIST
        sparse_m <- purrr::map(x, process_single_df)
        return(sparse_m)
    }
}

#### ---- Utilities for tests and examples ----####
#' A utility function to unzip and use example file systems included in the
#' package
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' From `ISAnalytics 1.5.4` this function is defunct, since the package
#' doesn't include example tabular files anymore. Use the function
#' `generate_default_folder_structure()` to generate a default folder
#' structure for running tests and play with the package import functions.
#' If you don't need to test import functions, you can simply load
#' package included data via `data("integration_matrices")` or
#' `data("association_file")`.
#'
#' @param zipfile The zipped file to decompress
#' @param name The name of the folder in the zipped archive ("fs" or "fserr")
#'
#' @export
#'
#' @return A path to reference
#' @keywords internal
unzip_file_system <- function(zipfile, name) {
    lifecycle::deprecate_stop(
        when = "1.5.4",
        what = "unzip_file_system()",
        with = "generate_default_folder_structure()",
        details = "Function will be removed in the next release cycle."
    )
}

#' Generate a default folder structure, following VISPA2 standards
#'
#' @description
#' The function produces a folder structure in the file system at the provided
#' path that respects VISPA2 standards, with package-included data.
#'
#' @section Required tags:
#' The function will explicitly check for the presence of these tags:
#'
#' ```{r echo=FALSE, results="asis"}
#' all_tags <- available_tags()
#' needed <- unique(all_tags[purrr::map_lgl(eval(sym("needed_in")),
#'  ~ "generate_default_folder_structure" %in% .x)][["tag"]])
#'  cat(paste0("* ", needed, collapse="\n"))
#' ```
#'
#' @param type One value between `"correct"`, `"incorrect"` and `"both"`.
#' Tells the function wheter to produce a correct structure or introduce some
#' errors (mainly for testing purposes).
#' @param dir Path to the folder in which the structure will be produced
#' @param af Either `"default"` for the association file provided as example
#' in the package or a custom association file as a data frame
#' @param matrices Either `"default"` for integration matrices
#' provided as example
#' in the package or a custom multi-quantification matrix
#'
#' @importFrom rlang sym
#'
#' @family Utilities
#' @return A named list containing the path to the association file and
#' the path to the top level folder(s) of the structure
#' @export
#'
#' @examples
#' fs_path <- generate_default_folder_structure(type = "correct")
#' fs_path
generate_default_folder_structure <- function(type = "correct",
    dir = tempdir(),
    af = "default",
    matrices = "default") {
    stopifnot(is.character(type))
    type <- type[1]
    stopifnot(type %in% c("correct", "incorrect", "both"))
    stopifnot(is.character(dir))
    dir <- dir[1]
    stopifnot(is.data.frame(af) || (is.character(af) & af == "default"))
    stopifnot(is.data.frame(matrices) || (is.character(matrices) &
        matrices == "default"))

    # Process association file
    association_file <- .process_af_for_gen(af)
    tag_list <- association_file$tag_list
    association_file <- association_file$af

    cols_selected <- colnames(
        association_file
    )[colnames(association_file) %in% association_file_columns()]
    association_file_reduced <- association_file %>%
        dplyr::select(dplyr::all_of(cols_selected))

    # Process matrices
    sep_matrices <- .process_m_for_gen(matrices, association_file, tag_list)

    # Process VISPA stats
    sep_stats <- .process_iss_for_gen(association_file, tag_list)

    if (is.null(sep_stats) & getOption("ISAnalytics.verbose", TRUE) == TRUE) {
        no_stats_msg <- c(paste("No VISPA2 stats column found in af,
                              skipping this step"))
        rlang::inform(no_stats_msg, class = "no_stats_warn")
    }
    # ---- Write files
    if (!dir_exists(dir)) {
        fs::dir_create(dir)
    }
    af_tmp_file <- fs::path(dir, "asso_file.tsv")
    readr::write_tsv(association_file_reduced, file = af_tmp_file, na = "")

    if (type == "both") {
        path_corr <- .generate_correct(dir, sep_matrices, sep_stats)
        path_incorr <- .generate_incorrect(dir, sep_matrices, sep_stats)
        return(list(
            af = af_tmp_file, root_corr = path_corr,
            root_inc = path_incorr
        ))
    }
    root <- if (type == "correct") {
        .generate_correct(dir, sep_matrices, sep_stats)
    } else {
        .generate_incorrect(dir, sep_matrices, sep_stats)
    }
    return(list(af = af_tmp_file, root = root))
}


#### ---- Configurations and others ----####
#' Export a dynamic vars settings profile.
#'
#' @description This function allows exporting the currently set dynamic
#' vars in json format so it can be quickly imported later. Dynamic
#' variables need to be properly set via the setter functions before calling
#' the function. For more details, refer to the dedicated vignette
#' `vignette("workflow_start", package="ISAnalytics")`.
#'
#' @param folder The path to the folder where the file should be saved. If
#' the folder doesn't exist, it gets created automatically
#' @param setting_profile_name A name for the settings profile
#'
#' @return `NULL`
#' @family Utilities
#' @export
#'
#' @examples
#' tmp_folder <- tempdir()
#' export_ISA_settings(tmp_folder, "DEFAULT")
export_ISA_settings <- function(folder, setting_profile_name) {
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
        rlang::abort(.missing_pkg_error("jsonlite"))
    }
    jsonify <- function(df) {
        tmp <- df %>%
            dplyr::mutate(transform = purrr::map(.data$transform, ~ {
                if (is.null(.x)) {
                    "NULL"
                } else {
                    rlang::f_text(.x)
                }
            }))
        tmp
    }
    all_specs_json <- list(
      mandatory_IS_vars = jsonify(mandatory_IS_vars(TRUE)),
      annotation_IS_vars = jsonify(annotation_IS_vars(TRUE)),
      association_file_columns = jsonify(association_file_columns(TRUE)),
      iss_stats_specs = jsonify(iss_stats_specs(TRUE)),
      matrix_file_suffixes = matrix_file_suffixes()
    )
    fs::dir_create(folder)
    file_name <- paste0(setting_profile_name, "_ISAsettings.json")
    all_specs_json %>%
        jsonlite::write_json(path = fs::path(folder, file_name))
    if (getOption("ISAnalytics.verbose", TRUE) == TRUE) {
        success_msg <- c("Settings profile correctly saved",
            i = paste("Saved at:", fs::path(folder, file_name))
        )
        rlang::inform(success_msg)
    }
}

#' Import a dynamic vars settings profile.
#'
#' @description The function allows the import of an existing dynamic vars
#' profile in json format. This is a quick and convenient way to set up
#' the workflow, alternative to specifying lookup tables manually through
#' the corresponding setter functions. For more details,
#' refer to the dedicated vignette
#' `vignette("workflow_start", package="ISAnalytics")`.
#'
#' @param path The path to the json file on disk
#'
#' @return `NULL`
#' @family Utilities
#' @export
#'
#' @examples
#' tmp_folder <- tempdir()
#' export_ISA_settings(tmp_folder, "DEFAULT")
#' import_ISA_settings(fs::path(tmp_folder, "DEFAULT_ISAsettings.json"))
#' reset_dyn_vars_config()
import_ISA_settings <- function(path) {
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
        rlang::abort(.missing_pkg_error("jsonlite"))
    }
    lookup_tbls <- jsonlite::fromJSON(path)
    unjsonify <- function(df, tbl_name) {
        tmp <- df %>%
            tibble::as_tibble() %>%
            mutate(transform = purrr::map(.data$transform, ~ {
                if (.x == "NULL") {
                    return(NULL)
                } else {
                    return(rlang::new_formula(NULL, rlang::parse_expr(.x)))
                }
            }))
        if (tbl_name == "association_file_columns") {
            set_af_columns_def(tmp)
        } else {
            function_name <- paste0("set_", tbl_name)
            rlang::exec(function_name, specs = tmp)
        }
    }
    matrix_suffixes_tbl <- lookup_tbls$matrix_file_suffixes
    lookup_tbls <- lookup_tbls[names(lookup_tbls) != "matrix_file_suffixes"]
    purrr::walk2(lookup_tbls, names(lookup_tbls), unjsonify)
    options(ISAnalytics.matrix_file_suffix = tibble::as_tibble(
        matrix_suffixes_tbl
    ))
    if (getOption("ISAnalytics.verbose", TRUE)) {
        rlang::inform("Matrix suffixes specs successfully changed")
    }
}

#' Launch the shiny application NGSdataExplorer.
#'
#' @return Nothing
#' @export
#'
#' @examples
#' \dontrun{
#' NGSdataExplorer()
#' }
NGSdataExplorer <- function() {
    required_pkgs <- c("shiny", "shinyWidgets", "datamods", "DT", "bslib")
    installed <- purrr::map_lgl(
        required_pkgs,
        ~ requireNamespace(.x, quietly = TRUE)
    )
    if (any(installed == FALSE)) {
        rlang::abort(.missing_pkg_error(required_pkgs[!installed]))
    }
    app <- shiny::shinyApp(ui = ui, server = server)
    shiny::runApp(app)
}
