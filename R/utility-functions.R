#------------------------------------------------------------------------------#
# Utility functions
#------------------------------------------------------------------------------#

#' Title
#'
#' @param tags
#'
#' @return
#' @export
#'
#' @examples
inspect_tags <- function(tags) {
  all_tags <- available_tags()
  desc_msg <- purrr::map(tags, ~ {
    tag_desc <- all_tags[eval(sym("tag")) == .x][["description"]]
    if (purrr::is_empty(tag_desc)) {
      return(
        c(paste("* TAG:", .x),
          x = "Tag not found in available tags")
      )
    }
    functions <- all_tags[eval(sym("tag")) == .x][["needed_in"]][[1]]
    functions <- paste0(functions, collapse = ", ")
    c(paste("* TAG:", .x),
      i = paste("Description:", tag_desc),
      i = paste("Functions that use it:", functions))
  })
  purrr::walk(desc_msg, ~ rlang::inform(.x, class = "tag_inspect"))
}


#' Define custom mandatory IS vars.
#'
#' This function lets the user specify the name and types of the mandatory
#' is vars: these variables uniquely identify an integration site within
#' a matrix.
#' Default values are `chr (char), integration_locus (int), strand (char)`.
#'
#' @details
#' The user can supply specifications in the form of a named vector or a
#' data frame.
#'
#' ## Named vector
#' When using a named vector, names should be the names of the columns,
#' values should be the type associated with each column in the form
#' of a string. The vector gets automatically converted into a data frame
#' with the right format (default values for the columns `transform` and
#' `flag` are `NULL` and `required` respectively).
#'
#' ## Data frame
#' When supplying a data frame, users must follow the format below:
#'
#' ```{r echo=FALSE}
#' print(tibble::tribble(
#' ~ names, ~ types, ~ transform, ~ flag,
#' "<col_name>", "<col_type>", ~ .x, "<required/optional>"
#' ))
#' ```
#'
#' Where the column `names` contains the names of the column,
#' the `types` column contains the column classes,
#' the column `transform` can hold a function or
#' a purrr-style lambda
#' that will be applied to the corresponding column
#' after import the column `flag` can contain for each column either the
#' value `required` or `optional`, and it is used in certain functions
#' to determine the minimal set of columns required.
#' Note that a valid function for the `transform` column
#'  takes as input a vector (the column)
#' and outputs a transformed vector that has the same length as the original
#' one. If no function should be applied, transform field should be set to
#' `NULL`. Example:
#'
#' ```{r}
#' tibble::tribble(
#' ~ names, ~ types, ~ transform, ~ flag,
#' "col1", "char", ~ stringr::str_to_upper(.x), "required",
#' "col2", "int", NULL, "required"
#' )
#' ```
#'
#' ## Types specification
#'
#' * Integer: `int`
#' * Character or string: `char`
#' * Numeric (double): `numeric`
#' * Logical (True/False): `logi`
#' * Factor: `factor`
#' * Date or date time: either generic `date` or any of the supported
#' date types in `date_formats()`.
#' See also the documentation of package \code{\link{lubridate}} for more info.
#' If `date` is supplied, the format of the date is either guessed or a
#' global specification supplied to the appropriate function is used for
#' parsing.
#'
#' @param specs Either a named vector or a data frame with specific format.
#' See details.
#'
#' @return Nothing
#' @seealso \code{\link{date_formats}}, \code{\link{reset_mandatory_IS_vars}},
#' \code{\link{mandatory_IS_vars}}
#'
#' @export
#'
#' @examples
#' set_mandatory_IS_vars(c(chr = "char"))
#' print(mandatory_IS_vars(TRUE))
#' reset_mandatory_IS_vars()
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
    if (getOption("ISAnalytics.verbose")) {
        rlang::inform("Mandatory IS vars successfully changed")
    }
}

#' Resets mandatory IS vars to the default value.
#'
#' @return Nothing
#' @export
#'
#' @seealso \code{\link{mandatory_IS_vars}},
#' \code{\link{set_mandatory_IS_vars}}
#'
#' @examples
#' set_mandatory_IS_vars(c(chr = "char"))
#' print(mandatory_IS_vars())
#' reset_mandatory_IS_vars()
#' print(mandatory_IS_vars())
reset_mandatory_IS_vars <- function() {
    options(ISAnalytics.mandatory_is_vars = "default")
    if (getOption("ISAnalytics.verbose")) {
        rlang::inform("Mandatory IS vars reset to default")
    }
}

#' Names of mandatory variables for an integration matrix.
#'
#' Contains the names of the columns that need to be present in order for a
#' tibble to be considered an integration matrix. To set
#' custom mandatory variables see the function
#' \code{\link{set_mandatory_IS_vars}}.
#'
#' @param include_types If set to `TRUE` returns both the names and the types
#' associated, otherwise returns only a character vector of names
#'
#' @return A character vector or a data frame
#' @export
#'
#' @examples
#' # Names only
#' mandatory_IS_vars()
#'
#' # Names and types
#' mandatory_IS_vars(TRUE)
mandatory_IS_vars <- function(include_types = FALSE) {
    opt <- getOption("ISAnalytics.mandatory_is_vars")
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

#' Define custom gene annotation vars.
#'
#' This function lets the user specify the name and types of the gene
#' annotation fields.
#' Default values are `GeneName (char), GeneStrand (char)`.
#'
#' @inherit set_mandatory_IS_vars details
#'
#' @param specs Either a named vector or a data frame with specific format.
#' See details.
#'
#' @return Nothing
#' @seealso \code{\link{date_formats}}, \code{\link{reset_annotation_IS_vars}},
#' \code{\link{annotation_IS_vars}}
#'
#' @export
#'
#' @examples
#' set_annotation_IS_vars(c(GeneName = "char"))
#' print(annotation_IS_vars(TRUE))
#' reset_annotation_IS_vars()
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
    if (getOption("ISAnalytics.verbose")) {
        rlang::inform("Annotation IS vars successfully changed")
    }
}

#' Resets annotation IS vars to the default value.
#'
#' @return Nothing
#' @export
#'
#' @seealso \code{\link{annotation_IS_vars}},
#' \code{\link{set_annotation_IS_vars}}
#'
#' @examples
#' set_annotation_IS_vars(c(GeneName = "char"))
#' print(annotation_IS_vars(TRUE))
#' reset_annotation_IS_vars()
#' print(annotation_IS_vars())
reset_annotation_IS_vars <- function() {
    options(ISAnalytics.genomic_annotation_vars = "default")
    if (getOption("ISAnalytics.verbose")) {
        rlang::inform("Annotation IS vars reset to default")
    }
}

#' Names of the genomic annotation variables for an integration matrix.
#'
#' Default values are `GeneName (char), GeneStrand (char)`. To set
#' custom genomic annotation variables see the function
#' \code{\link{set_annotation_IS_vars}}.
#'
#' @inheritParams mandatory_IS_vars
#'
#' @return A character vector or a data frame
#' @export
#'
#' @examples
#' # Names only
#' annotation_IS_vars()
#'
#' # Names and types
#' annotation_IS_vars(TRUE)
annotation_IS_vars <- function(include_types = FALSE) {
    opt <- getOption("ISAnalytics.genomic_annotation_vars")
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


#' Define custom association file columns.
#'
#' This function lets the user specify the name and types of association
#' file fields
#'
#' @inherit set_mandatory_IS_vars details
#'
#' @param specs Either a named vector or a data frame with specific format.
#' See details.
#'
#' @return Nothing
#' @seealso \code{\link{date_formats}}, \code{\link{reset_af_columns_def}},
#' \code{\link{association_file_columns}}
#'
#' @export
#'
#' @examples
#' temp_af_cols <- tibble::tribble(
#'     ~names, ~types, ~transform,
#'     "ProjectID", "char", NULL,
#'     "SubjectID", "char", NULL,
#'     "TimePoint", "int", NULL
#' )
#' set_af_columns_def(temp_af_cols)
#' print(association_file_columns(TRUE))
#' reset_af_columns_def()
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
    if (getOption("ISAnalytics.verbose")) {
        rlang::inform("Association file columns specs successfully changed")
    }
}

#' Resets the association file columns definitions to the default value.
#'
#' @return Nothing
#' @export
#'
#' @seealso \code{\link{association_file_columns}},
#' \code{\link{set_af_columns_def}}
#'
#' @examples
#' temp_af_cols <- tibble::tribble(
#'     ~names, ~types, ~transform,
#'     "ProjectID", "char", NULL,
#'     "SubjectID", "char", NULL,
#'     "TimePoint", "int", NULL
#' )
#' set_af_columns_def(temp_af_cols)
#' print(association_file_columns(TRUE))
#' reset_af_columns_def()
reset_af_columns_def <- function() {
    options(ISAnalytics.genomic_annotation_vars = "default")
    if (getOption("ISAnalytics.verbose")) {
        rlang::inform("Association file columns specs reset to default")
    }
}

#' Names of the columns in the association file.
#'
#' All the names of the columns present in the association file.
#'
#' @return A character vector or a data frame
#' @export
#'
#' @inheritParams mandatory_IS_vars
#'
#' @examples
#' # Names only
#' association_file_columns()
#'
#' # Names and types
#' association_file_columns(TRUE)
association_file_columns <- function(include_types = FALSE) {
    opt <- getOption("ISAnalytics.af_specs")
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

#' Title
#'
#' @param specs
#'
#' @return
#' @export
#'
#' @examples
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
    if (getOption("ISAnalytics.verbose")) {
        rlang::inform("ISS stats specs successfully changed")
    }
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
reset_iss_stats_specs <- function() {
    options(ISAnalytics.iss_stats_specs = "default")
    if (getOption("ISAnalytics.verbose")) {
        rlang::inform("ISS stats specs reset to default")
    }
}

#' Title
#'
#' @param include_types
#'
#' @return
#' @export
#'
#' @examples
iss_stats_specs <- function(include_types = FALSE) {
    opt <- getOption("ISAnalytics.iss_stats_specs")
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


#' Title
#'
#' @return
#' @export
#'
#' @examples
matrix_file_suffixes <- function() {
  opt <- getOption("ISAnalytics.matrix_file_suffix")

  if (all(opt == "default")) {
    defaults_params <- as.list(formals(set_matrix_file_suffixes))
    default_vars <- do.call(.generate_suffix_specs, args = defaults_params)
    return(default_vars)
  }

  return(opt)
}

#' Title
#'
#' @param quantification_suffix
#' @param annotation_suffix
#' @param file_ext
#' @param glue_file_spec
#'
#' @return
#' @export
#'
#' @examples
set_matrix_file_suffixes <- function(
  quantification_suffix = list(seqCount = "seqCount",
                               fragmentEstimate = "fragmentEstimate",
                               barcodeCount = "barcodeCount",
                               cellCount = "cellCount",
                               ShsCount = "ShsCount"),
  annotation_suffix = list(annotated = ".no0.annotated",
                           not_annotated = ""),
  file_ext = "tsv.gz",
  glue_file_spec = "{quantification_suffix}_matrix{annotation_suffix}.{file_ext}"
) {
  stopifnot(is.list(quantification_suffix) &&
              !is.null(names(quantification_suffix)))
  stopifnot(is.list(annotation_suffix) &&
              !is.null(names(annotation_suffix)))
  stopifnot(is.character(file_ext))
  stopifnot(is.character(glue_file_spec))
  glue_file_spec <- glue_file_spec[1]
  if (!all(quantification_types() %in% names(quantification_suffix))) {
    warn_msg <- c("Warning: some quantification types are missing from specs",
                  i = paste("Missing: ",
                            quantification_types()[!quantification_types() %in%
                                                     names(quantification_suffix)]),
                  i = paste("Missing quantifications specs may cause problems",
                            "when trying to import files.",
                            "See documentation with `?set_matrix_file_suffixes`")
                  )
    rlang::inform(warn_msg, class = "missing_quant_specs")
  }
  quantification_suffix <- quantification_suffix[
    names(quantification_suffix) %in% quantification_types()]
  if (!all(c("annotated", "not_annotated") %in% names(annotation_suffix))) {
    err_msg <- c(paste("Error: a value for both 'annotated'",
                       "and 'not_annotated' is requested. Quitting."))
    rlang::abort(err_msg, class = "miss_annot_suff_specs")
  }

  final_specs <- .generate_suffix_specs(quantification_suffix,
                                        annotation_suffix,
                                        file_ext,
                                        glue_file_spec)

  options(ISAnalytics.matrix_file_suffix = final_specs)
  if (getOption("ISAnalytics.verbose")) {
    rlang::inform("Matrix suffixes specs successfully changed")
  }
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
reset_matrix_file_suffixes <- function() {
  options(ISAnalytics.matrix_file_suffix = "default")
  if (getOption("ISAnalytics.verbose")) {
    rlang::inform("Matrix suffixes specs reset to default")
  }
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
pcr_id_column <- function() {
  af_cols_specs <- association_file_columns(TRUE)
  selected_tags <- .check_required_cols(list("pcr_repl_id" = "char"),
                                        af_cols_specs, "error")
  return(selected_tags$names)
}

# Function used to apply arbitrary transformations on columns.
# Expects a named list where names are names of columns and values are
# purrr style lambdas
#' Title
#'
#' @param df
#' @param transf_list
#'
#' @return
#' @export
#'
#' @examples
transform_columns <- function(df, transf_list) {
  transf_list <- transf_list[names(transf_list) %in% colnames(df)]
  if (purrr::is_empty(transf_list)) {
    return(df)
  }
  transf_list_mod <- purrr::map2(names(transf_list), transf_list, ~ {
    rlang::expr(rlang::as_function(!!.y)(!!rlang::sym(.x)))
  }) %>%
    purrr::set_names(names(transf_list))
  withCallingHandlers({
    withRestarts({
      return(dplyr::mutate(df, !!!transf_list_mod))
    },
    skip_transform = function() {
      string_list <- paste(names(transf_list), transf_list,
                           sep = " = ", collapse = "; " )
      err_msg <- c(paste("Encountered a problem while applying transformations",
                         "to columns. Skipping this."),
                   i = paste("Check the correctness of your input,",
                             "see the documentation"),
                   i = paste("Your input: ", string_list))
      if (getOption("ISAnalytics.verbose")) {
        rlang::inform(err_msg, class = "skip_col_transform")
      }
      return(df)
    })
  }, error = function(cnd) {
    res <- findRestart("skip_transform")
    invokeRestart(res)
  })
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
  if (any(na_introduced) & getOption("ISAnalytics.verbose") == TRUE) {
    rlang::inform(.nas_introduced_msg(), class = "comp_nas")
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
      dplyr::rename(Value = .data[[quant]])
  }) %>% purrr::set_names(names(num_cols))
  separated
}



#' Creates a blank association file.
#'
#' This function is useful if you want a blank association file to start using
#' both Vispa2 and this package or simply if you want a correct framework to
#' fix a malformed association file you have already.
#'
#' @param path The path on disk where the file should be written
#' @importFrom fs path_dir as_fs_path dir_create dir_exists
#' @importFrom readr write_tsv
#'
#' @return returns NULL
#' @export
#' @family Utility functions
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

#' Creates a reduced association file for Vispa2 run,
#' given project and pool
#'
#' The function selects the appropriate columns and prepares a file for the
#' launch of Vispa2 pipeline for each project/pool pair specified.
#'
#' @details Note: the function is vectorized, meaning you can specify more than
#' one project and more than one pool as vectors of characters, but you must
#' ensure that:
#' * Both `project` and `pool` vectors have the same length
#' * You correclty type names in corresponding positions, for example
#' c("CLOEXP", "PROJECT1100", "PROJECT1100") - c("POOL6", "ABX-LR-PL5-POOL14-1",
#' "ABX-LR-PL6-POOL15-1"). If you type a pool in the position of a corresponding
#' project that doesn't match no file will be produced since that pool doesn't
#' exist in the corresponding project.
#'
#' @param association_file The imported association file (via
#' `import_association_file`)
#' @param project A vector of characters containing project names
#' @param pool A vector of characters containing pool names.
#' **NOTE: the names should refer to the values contained in the
#' PoolID column of the association file and NOT the
#' concatenatePoolIDSeqRun column!**
#' @param path A single string representing the path to the folder where files
#' should be written. If the folder doesn't exist it will be created.
#' @importFrom fs as_fs_path file_exists dir_create path
#' @importFrom purrr map2 set_names walk2
#' @importFrom dplyr select filter mutate all_of
#' @importFrom readr write_tsv
#' @importFrom magrittr `%>%`
#' @family Utility functions
#'
#' @return returns NULL
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
                                       duplicate_politic = "error")
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
        dplyr::filter(.data[[proj_col]] == x,
                      .data[[pool_col]] == y) %>%
        dplyr::mutate(TagID2 = .data[[tag_id_col]]) %>%
        dplyr::mutate(PoolName = dplyr::if_else(
          !is.na(.data[[concat_col]]),
          .data[[concat_col]],
          .data[[pool_col]]
        )) %>%
        dplyr::select(.data[[tag_id_col]], .data$TagID2, .data[[tissue_col]],
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
#' \lifecycle{stable}
#' This function is particularly useful when a sparce matrix structure
#' is needed by a specific function (mainly from other packages).
#'
#' @param x A single tidy integration matrix or a list of integration
#' matrices. Supports also multi-quantification matrices
#' obtained via \link{comparison_matrix}
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
#'
#' @importFrom tidyr pivot_wider
#' @importFrom purrr is_empty set_names walk map
#' @importFrom dplyr select
#'
#' @family Utility functions
#'
#' @return Depending on input, 2 possible outputs:
#' * A single sparse matrix (tibble) if input is a single quantification
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
              names_from = .data[[key]],
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
#' This utility function is a simple shortcut to create a temporary directory,
#' unzip and reference the examples file systems included in the package for
#' testing purposes.
#'
#' @param zipfile The zipped file to decompress
#' @param name The name of the folder in the zipped archive ("fs" or "fserr")
#' @importFrom zip unzip
#' @family Utility functions
#' @export
#'
#' @return A path to reference
#' @examples
#' root_pth <- system.file("extdata", "fs.zip", package = "ISAnalytics")
#' root <- unzip_file_system(root_pth, "fs")
unzip_file_system <- function(zipfile, name) {
  lifecycle::deprecate_stop(
    when = "1.5.4",
    what = "unzip_file_system()",
    with = "generate_default_folder_structure()",
    details = "Function will be removed in the next release cycle."
  )
}

# TODO - add the incorrect version for testing (and both) +
# add possibility of supplying custom af and matrices (for tests)
#' Title
#'
#' @param type
#' @param dir
#' @param af
#' @param matrices
#'
#' @return
#' @export
#'
#' @examples
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

    if (is.null(sep_stats) & getOption("ISAnalytics.verbose") == TRUE) {
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
      return(list(af = af_tmp_file, root_corr = path_corr,
                  root_inc = path_incorr))
    }
    root <- if (type == "correct") {
      .generate_correct(dir, sep_matrices, sep_stats)
    } else {
      .generate_incorrect(dir, sep_matrices, sep_stats)
    }
    return(list(af = af_tmp_file, root = root))
}
