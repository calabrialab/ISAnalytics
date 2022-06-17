#------------------------------------------------------------------------------#
# Internal functions
#------------------------------------------------------------------------------#
## All functions in this file are NOT exported, to be used internally only.

#### ---- Internals for dynamic variables  ----####
# Internal for setting mandatory or annotation is vars
# tbl_type must be in {mand_vars, annot_vars, af_vars, iss_vars}
#' @importFrom rlang sym
.new_IS_vars_checks <- function(specs, err, tbl_type) {
    new_vars <- if (is.data.frame(specs)) {
        # if data frame supplied
        colnames_ok <- all(c("names", "types", "transform", "flag", "tag")
        %in% colnames(specs))
        col_types_ok <- if (colnames_ok) {
            all(purrr::map_lgl(
                specs[c("names", "types", "flag", "tag")],
                is.character
            ))
        } else {
            FALSE
        }
        types_ok <- if (colnames_ok & col_types_ok) {
            all(specs$types %in% .types_mapping()[["types"]])
        } else {
            FALSE
        }
        flags_ok <- if (colnames_ok & col_types_ok) {
            all(specs$flag %in% c("required", "optional"))
        } else {
            FALSE
        }
        transform_col_ok <- if (colnames_ok & col_types_ok) {
            is.list(specs$transform) &&
                all(purrr::map_lgl(
                    specs$transform,
                    ~ {
                        is.null(.x) || rlang::is_formula(.x) || is.function(.x)
                    }
                ))
        } else {
            FALSE
        }
        if (!(colnames_ok & col_types_ok & types_ok & transform_col_ok &
            flags_ok)) {
            rlang::abort(err)
        }
        specs
    } else {
        # if vector supplied
        if (is.null(names(specs)) || !all(specs %in%
            .types_mapping()[["types"]])) {
            rlang::abort(err)
        }
        purrr::map2_dfr(
            names(specs), specs,
            ~ tibble::tibble_row(
                names = .x,
                types = .y,
                transform = list(NULL),
                flag = "required",
                tag = NA_character_
            )
        )
    }
    # -- check critical tags
    if (!is.null(tbl_type)) {
        available_tags <- available_tags()
        to_check <- available_tags[eval(sym("dyn_vars_tbl")) == tbl_type &
            !purrr::map_lgl(
                eval(sym("needed_in")), ~ purrr::is_empty(.x)
            ), ]
        if (!all(to_check$tag %in% new_vars$tag)) {
            missing_tags <- to_check$tag[!to_check$tag %in% new_vars$tag]
            inspect_cmd <- paste0(
                "inspect_tags(c(", paste0("'", missing_tags,
                    "'",
                    collapse = ","
                ),
                "))"
            )
            warn <- c("Warning: important tags missing",
                i = paste(
                    "Some tags are required for proper execution",
                    "of some functions. If these tags are not",
                    "provided, execution of dependent functions",
                    "might fail.",
                    "Review your inputs carefully."
                ),
                i = paste("Missing tags:", paste0(missing_tags,
                    collapse = ", "
                )),
                i = paste0(
                    "To see where these are involved type `",
                    inspect_cmd, "`"
                )
            )
            rlang::warn(warn, class = "missing_crit_tags")
        }
    }
    new_vars
}

# Applies transformations on columns as specified in variables specs
# expects specs to be in data frame format
.apply_col_transform <- function(df, specs) {
    # Extract and associate names and transf
    non_null <- purrr::pmap(specs, function(names, types,
    transform, flags, tag) {
        if (is.null(transform)) {
            return(NULL)
        }
        return(transform)
    }) %>%
        purrr::set_names(specs$names)
    # Retain only non-null transform
    non_null <- non_null[purrr::map_lgl(non_null, ~ !is.null(.x))]
    if (length(non_null) > 0) {
        # if there are transf to apply
        apply_transform <- function(col, col_name) {
            if (!col_name %in% names(non_null)) {
                return(col)
            }
            transformation <- non_null[[col_name]]
            t <- if (rlang::is_formula(transformation)) {
                unlist(purrr::map(col, transformation))
            } else {
                # if it is a function
                do.call(what = transformation, args = list(col))
            }
        }
        df <- purrr::map2_dfc(df, colnames(df), apply_transform)
    }
    df
}

# Internal to quickly convert columns marked as dates with lubridate
.convert_dates <- function(date_vec, format) {
    if (format %in% date_formats()) {
        parsed <- lubridate::parse_date_time(date_vec, format)
        if (format %in% c(
            "ymd_hms", "ymd_hm", "ymd_h", "dmy_hms", "dmy_hm",
            "dmy_h", "mdy_hms", "mdy_hm", "mdy_h",
            "ydm_hms", "ydm_hm", "ydm_h"
        )) {
            return(parsed)
        } else {
            return(lubridate::as_date(parsed))
        }
    }
    if (format == "date") {
        # Guesses the format
        return(lubridate::as_date(as.character(date_vec)))
    }
}

# Generates a data frame holding file suffixes combinations for
# matrices
.generate_suffix_specs <- function(quantification_suffix,
    annotation_suffix,
    file_ext,
    glue_file_spec) {
    # Calculate combinations
    combinations <- purrr::cross(list(
        quantification_suffix = quantification_suffix,
        annotation_suffix = annotation_suffix,
        file_ext = file_ext
    ))
    glue_combo <- function(x) {
        quantif <- names(quantification_suffix[
            quantification_suffix == x$quantification_suffix
        ])
        matrix_type <- names(annotation_suffix[
            annotation_suffix == x$annotation_suffix
        ])
        quantification_suffix <- x$quantification_suffix
        annotation_suffix <- x$annotation_suffix
        file_ext <- x$file_ext
        glued <- glue::glue(glue_file_spec)
        tibble::tibble(
            quantification = quantif,
            matrix_type = matrix_type,
            file_suffix = glued
        )
    }
    final_specs <- purrr::map_df(combinations, glue_combo)
    final_specs
}

#### ---- Internals for utilities ----####
#---- USED IN : generate_default_folder_structure ----
.process_af_for_gen <- function(af) {
    association_file <- if (!is.data.frame(af) && all(af == "default")) {
        af_sym <- "association_file"
        utils::data(list = af_sym, envir = rlang::current_env())
        rlang::eval_tidy(rlang::sym(af_sym))
    } else {
        af
    }
    required_af_tags <- c(
        "pcr_repl_id" = "char",
        "vispa_concatenate" = "char",
        "tag_seq" = "char",
        "project_id" = "char"
    )
    tag_to_cols_af <- .check_required_cols(required_af_tags,
        vars_df = association_file_columns(TRUE),
        duplicate_politic = "error"
    )
    if (!all(tag_to_cols_af$names %in% colnames(association_file))) {
        rlang::abort(.missing_req_cols(
            tag_to_cols_af$names,
            tag_to_cols_af$names[!tag_to_cols_af$names %in%
                colnames(association_file)]
        ), class = "missing_req_col_err")
    }
    tag_to_cols_af_list <- purrr::map(tag_to_cols_af$tag, ~ tag_to_cols_af %>%
        dplyr::filter(.data$tag == .x) %>%
        dplyr::pull(.data$names)) %>%
        purrr::set_names(tag_to_cols_af$tag)

    return(list(af = association_file, tag_list = tag_to_cols_af_list))
}

.process_m_for_gen <- function(matrices, af, tag_list) {
    integration_matrices <- if (!is.data.frame(matrices) &&
        all(matrices == "default")) {
        matrices_sym <- "integration_matrices"
        utils::data(list = matrices_sym, envir = rlang::current_env())
        rlang::eval_tidy(rlang::sym(matrices_sym))
    } else {
        matrices
    }
    matrix_required_cols <- c(
        mandatory_IS_vars(),
        tag_list$pcr_repl_id
    )
    if (!all(matrix_required_cols %in% colnames(integration_matrices))) {
        rlang::abort(.missing_req_cols(
            matrix_required_cols,
            matrix_required_cols[!matrix_required_cols %in%
                colnames(integration_matrices)]
        ), class = "missing_req_col_err")
    }
    # First separate by project
    sep_matrices <- integration_matrices %>%
        dplyr::left_join(af %>%
            dplyr::select(
                dplyr::all_of(c(
                    tag_list$pcr_repl_id,
                    tag_list$vispa_concatenate,
                    tag_list$project_id
                ))
            ),
        by = tag_list$pcr_repl_id
        ) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(tag_list$project_id)))
    proj_names <- sep_matrices %>%
        dplyr::group_keys() %>%
        dplyr::pull(dplyr::all_of(tag_list$project_id))
    sep_matrices <- sep_matrices %>%
        dplyr::group_split(.keep = FALSE) %>%
        purrr::set_names(proj_names)
    # Then separate by pool
    sep_matrices <- purrr::map(sep_matrices, ~ {
        tmp <- .x %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(tag_list$vispa_concatenate)))
        pool_names <- tmp %>%
            dplyr::group_keys() %>%
            dplyr::pull(dplyr::all_of(tag_list$vispa_concatenate))
        tmp <- tmp %>%
            dplyr::group_split(.keep = FALSE) %>%
            purrr::set_names(pool_names)
        return(tmp)
    })
    return(sep_matrices)
}

# NOTE: the function expects vispa stats to be embedded in the association
# file, not provided as separate files
.process_iss_for_gen <- function(association_file, tag_list) {
    # Find VISPA stats minimal required columns
    minimal_iss_cols <- iss_stats_specs(TRUE) %>%
        dplyr::filter(
            .data$flag == "required",
            !.data$tag %in% c(
                "vispa_concatenate",
                "tag_seq"
            )
        ) %>%
        dplyr::pull(.data$names)
    # If the minimal cols are not found in af, nothing to do, return
    if (!all(minimal_iss_cols %in% colnames(association_file))) {
        return(NULL)
    }
    # Check required tags
    iss_specs <- .check_required_cols(
        required_tags = c(
            "vispa_concatenate" = "char",
            "tag_seq" = "char"
        ),
        vars_df = iss_stats_specs(TRUE),
        duplicate_politic = "error"
    )
    iss_tags <- purrr::map(iss_specs$tag, ~ iss_specs %>%
        dplyr::filter(.data$tag == .x) %>%
        dplyr::pull(.data$names)) %>%
        purrr::set_names(iss_specs$tag)
    iss_cols <- iss_stats_specs(TRUE) %>%
        dplyr::filter(!.data$tag %in% c("vispa_concatenate", "tag_seq")) %>%
        dplyr::pull(.data$names)
    iss_cols <- colnames(association_file)[colnames(association_file) %in%
        iss_cols]
    iss_data <- association_file %>%
        dplyr::select(
            dplyr::all_of(c(
                tag_list$project_id,
                tag_list$vispa_concatenate,
                tag_list$tag_seq,
                iss_cols
            ))
        ) %>%
        dplyr::rename(
            !!iss_tags$tag_seq := tag_list$tag_seq,
            !!iss_tags$vispa_concatenate :=
                tag_list$vispa_concatenate
        ) %>%
        dplyr::distinct() %>%
        dplyr::filter(!dplyr::if_all(dplyr::all_of(iss_cols), is.na))
    # Separate by project
    stats_split <- iss_data %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(tag_list$project_id)))
    proj_names <- stats_split %>%
        dplyr::group_keys() %>%
        dplyr::pull(dplyr::all_of(tag_list$project_id))
    stats_split <- stats_split %>%
        dplyr::group_split(.keep = FALSE) %>%
        purrr::set_names(proj_names)
    # Separate by pool
    stats_split <- purrr::map(stats_split, ~ {
        tmp <- .x %>%
            dplyr::group_by(dplyr::across(
                dplyr::all_of(iss_tags$vispa_concatenate)
            ))
        pool_names <- tmp %>%
            dplyr::group_keys() %>%
            dplyr::pull(dplyr::all_of(iss_tags$vispa_concatenate))
        tmp <- tmp %>%
            dplyr::group_split() %>%
            purrr::set_names(pool_names)
        return(tmp)
    })
    return(stats_split)
}

# Generates a correct standard file folder structure (either with package
# included data or provided)
.generate_correct <- function(dir, sep_matrices, sep_stats) {
    corr_fold <- fs::path(dir, "fs")
    fs::dir_create(corr_fold)
    ## -- for each project
    for (proj in names(sep_matrices)) {
        proj_fold <- fs::path(corr_fold, proj)
        quant_fold <- fs::path(proj_fold, "quantification")
        fs::dir_create(quant_fold)
        ## --- Write matrices
        purrr::walk2(sep_matrices[[proj]], names(sep_matrices[[proj]]), ~ {
            sparse <- as_sparse_matrix(.x)
            pool_fold <- fs::path(quant_fold, .y)
            fs::dir_create(pool_fold)
            fe_suffix <- dplyr::if_else(.is_annotated(.x),
                "fragmentEstimate_matrix.no0.annotated.tsv.gz",
                "fragmentEstimate_matrix.tsv.gz"
            )
            sc_suffix <- dplyr::if_else(.is_annotated(.x),
                "seqCount_matrix.no0.annotated.tsv.gz",
                "seqCount_matrix.tsv.gz"
            )
            prefix <- paste(proj, .y, sep = "_")
            readr::write_tsv(sparse$fragmentEstimate,
                file = fs::path(pool_fold, paste(prefix,
                    fe_suffix,
                    sep = "_"
                )),
                na = ""
            )
            readr::write_tsv(sparse$seqCount,
                file = fs::path(pool_fold, paste(prefix,
                    sc_suffix,
                    sep = "_"
                )),
                na = ""
            )
        })
    }
    # --- Write iss
    if (any(purrr::map_lgl(sep_stats, ~ !is.null(.x)))) {
        for (proj in names(sep_stats)) {
            proj_fold <- fs::path(corr_fold, proj)
            iss_fold <- fs::path(proj_fold, "iss")
            fs::dir_create(iss_fold)
            purrr::walk2(sep_stats[[proj]], names(sep_stats[[proj]]), ~ {
                pool_fold <- fs::path(iss_fold, .y)
                fs::dir_create(pool_fold)
                if (nrow(.x) > 0) {
                    filename <- paste0(
                        "stats.sequence_", proj,
                        ".", .y, ".tsv"
                    )
                    readr::write_tsv(.x,
                        file = fs::path(pool_fold, filename),
                        na = ""
                    )
                } else {
                    if (getOption("ISAnalytics.verbose", TRUE) == TRUE) {
                        warn_empty_iss <- c(paste0(
                            "Warning: empty stats for project '",
                            proj, "', pool '", .y, "'"
                        ),
                        i = paste(
                            "No files will be produced for",
                            "this pool, if this behaviour is",
                            "not desired check the association",
                            "file provided in input"
                        )
                        )
                        rlang::inform(warn_empty_iss, class = "warn_empty_iss")
                    }
                }
            })
        }
    }
    return(corr_fold)
}

# Generates an incorrect standard file folder structure (either with package
# included data or provided) - intended for testing purposes only
.generate_incorrect <- function(dir, sep_matrices, sep_stats) {
    err_fold <- fs::path(dir, "fserr")
    fs::dir_create(err_fold)
    # -- Generate error for projects
    if (length(sep_matrices) > 1) {
        names(sep_matrices)[1] <- paste0(names(sep_matrices)[1], "-err")
    }
    ## -- for each project
    for (proj in names(sep_matrices)) {
        proj_fold <- fs::path(err_fold, proj)
        quant_fold <- fs::path(proj_fold, "quantification")
        fs::dir_create(quant_fold)
        ### -- Generate pool error
        if (length(sep_matrices[[proj]]) > 1) {
            names(sep_matrices[[proj]])[1] <- paste0(
                names(sep_matrices[[proj]])[1], "-err"
            )
        }
        ### -- Choose a pool for matrix error
        pool_to_mod <- purrr::detect(
            names(sep_matrices[[proj]]),
            ~ stringr::str_detect(.x,
                "-err",
                negate = TRUE
            )
        )
        ## --- Write matrices
        purrr::walk2(sep_matrices[[proj]], names(sep_matrices[[proj]]), ~ {
            sparse <- as_sparse_matrix(.x)
            pool_fold <- fs::path(quant_fold, .y)
            fs::dir_create(pool_fold)
            fe_suffix <- dplyr::if_else(.is_annotated(.x),
                "fragmentEstimate_matrix.no0.annotated.tsv.gz",
                "fragmentEstimate_matrix.tsv.gz"
            )
            sc_suffix <- dplyr::if_else(.is_annotated(.x),
                "seqCount_matrix.no0.annotated.tsv.gz",
                "seqCount_matrix.tsv.gz"
            )
            prefix <- paste(proj, .y, sep = "_")
            if (!.y == pool_to_mod) {
                readr::write_tsv(sparse$fragmentEstimate,
                    file = fs::path(pool_fold, paste(prefix,
                        fe_suffix,
                        sep = "_"
                    )),
                    na = ""
                )
            }
            readr::write_tsv(sparse$seqCount,
                file = fs::path(pool_fold, paste(prefix,
                    sc_suffix,
                    sep = "_"
                )),
                na = ""
            )
        })
    }
    # --- Write iss
    if (any(purrr::map_lgl(sep_stats, ~ !is.null(.x)))) {
        for (proj in names(sep_stats)) {
            proj_fold <- fs::path(err_fold, proj)
            iss_fold <- fs::path(proj_fold, "iss")
            fs::dir_create(iss_fold)
            ### -- Generate pool error
            if (length(sep_stats[[proj]]) > 1) {
                names(sep_stats[[proj]])[1] <- paste0(
                    names(sep_stats[[proj]])[1], "-err"
                )
            }
            purrr::walk2(sep_stats[[proj]], names(sep_stats[[proj]]), ~ {
                pool_fold <- fs::path(iss_fold, .y)
                fs::dir_create(pool_fold)
                if (!.y == pool_to_mod) {
                    if (nrow(.x) > 0) {
                        filename <- paste0(
                            "stats.sequence_", proj,
                            ".", .y, ".tsv"
                        )
                        readr::write_tsv(.x,
                            file = fs::path(pool_fold, filename),
                            na = ""
                        )
                    } else {
                        if (getOption("ISAnalytics.verbose", TRUE) == TRUE) {
                            warn_empty_iss <- c(paste0(
                                "Warning: empty stats for project '",
                                proj, "', pool '", .y, "'"
                            ),
                            i = paste(
                                "No files will be produced for",
                                "this pool, if this behaviour is",
                                "not desired check the association",
                                "file provided in input"
                            )
                            )
                            rlang::inform(warn_empty_iss,
                                class = "warn_empty_iss"
                            )
                        }
                    }
                }
            })
        }
    }
    return(err_fold)
}


#### ---- Internals for multiple functions/general purpose ----####

# Returns the file format for each of the file paths passed as a parameter.
#' @importFrom fs path_ext
#' @importFrom tools file_path_sans_ext
#' @importFrom purrr is_empty
.check_file_extension <- function(file_path) {
    ## Get last portion of file name
    last <- fs::path_ext(file_path)
    compressed <- which(last %in% .compressed_formats())
    if (!purrr::is_empty(compressed)) {
        ## for compressed files try extracting the true extension
        file_path[compressed] <- tools::file_path_sans_ext(
            file_path[compressed]
        )
        last <- fs::path_ext(file_path)
    }
    return(last)
}

# Internal to check if required tags are in the specified variables.
# * required_tags is a vector of char, if supplied with names
#   in the form of c("tag" = "col_type") also types are checked
# * duplicate_politic is used to determine what to do in case a tag
#   is found more than once:
#     - "error": raises an error and returns
#     - "keep": keeps all columns, no error or msgs
#     - "first": keeps the first row it finds
#     - named vector with names as tags and values either "keep" for
#       keeping all duplicates, "error" to raise an error if duplicates
#       are found for that specific column, "first" to keep only first
.check_required_cols <- function(required_tags, vars_df, duplicate_politic) {
    error_msg <- function(missing_cols) {
        e <- c("Missing required tags or wrong types",
            x = paste(
                "Missing tags:",
                paste0(missing_cols, collapse = ", ")
            )
        )
    }
    # Check just colnames
    tags_names <- if (!is.null(names(required_tags))) {
        names(required_tags)
    } else {
        required_tags
    }
    cols_in_tags <- vars_df %>%
        dplyr::filter(.data$tag %in% tags_names)
    # No tags found
    if (nrow(cols_in_tags) == 0) {
        rlang::abort(error_msg(tags_names), class = "missing_tags_err")
    }
    missing_tags <- tags_names[!tags_names %in% cols_in_tags$tag]
    if (length(missing_tags) > 0) {
        rlang::abort(error_msg(missing_tags), class = "missing_tags_err")
    }
    # ---- If all tags are found
    ### Evaluate types if necessary
    if (!is.null(names(required_tags))) {
        cols_in_tags <- .eval_tag_types(cols_in_tags, required_tags)
    }

    ### Evaluate duplicates
    cols_in_tags <- .eval_tag_duplicates(cols_in_tags,
        duplicate_politic = duplicate_politic
    )
    cols_in_tags
}

# * duplicate_politic is used to determine what to do in case a tag
#   is found more than once:
#     - "error": raises an error and returns
#     - "keep": keeps all columns, no error or msgs
#     - "first": keeps the first row it finds
#     - named vector with names as tags and values either "keep" for
#       keeping all duplicates, "error" to raise an error if duplicates
#       are found for that specific column, "first" to keep only first
.eval_tag_duplicates <- function(df, duplicate_politic) {
    if (any(duplicate_politic != "keep")) {
        tags_split <- df %>%
            dplyr::group_by(.data$tag) %>%
            dplyr::group_split()

        first_politic <- function(sub_df) {
            first_row <- sub_df[1, ]
            warn <- c(paste0(
                "--> Duplicates found for tag '",
                sub_df$tag[1], "', only first match kept"
            ),
            i = paste(
                "Found:",
                paste0(sub_df$names, collapse = ", ")
            ),
            i = paste("Kept:", first_row$names)
            )
            return(list(type = "MOD", result = first_row, warning = warn))
        }

        error_politic <- function(sub_df) {
            err <- c(
                paste0(
                    "Duplicates not allowed for tag '", sub_df$tag[1],
                    "', found ", nrow(sub_df)
                ),
                utils::capture.output(print((sub_df)))
            )
            return(list(type = "ERROR", result = err))
        }
        # If named vector
        single_tag_eval <- function(sub_df) {
            if (nrow(sub_df) == 1) {
                return(list(type = "DF", result = sub_df))
            }
            choice <- duplicate_politic[sub_df$tag[1]]
            if (choice == "error") {
                return(error_politic(sub_df))
            }
            if (choice == "first") {
                return(first_politic(sub_df))
            }
            if (choice == "keep") {
                return(list(type = "DF", result = sub_df))
            }
        }

        dupl_eval <- if (!is.null(names(duplicate_politic))) {
            purrr::map(tags_split, single_tag_eval)
        } else if (all(duplicate_politic == "error")) {
            purrr::map(tags_split, ~ {
                if (nrow(.x) == 1) {
                    return(list(type = "DF", result = .x))
                }
                return(error_politic(.x))
            })
        } else {
            purrr::map(tags_split, ~ {
                if (nrow(.x) == 1) {
                    return(list(type = "DF", result = .x))
                }
                return(first_politic(.x))
            })
        }

        errors <- purrr::map(dupl_eval, ~ {
            if (.x$type == "ERROR") {
                return(.x$result)
            }
            NULL
        })
        if (any(purrr::map_lgl(errors, ~ !is.null(.x)))) {
            single_err <- c("Duplicates found for some tags",
                purrr::reduce(errors, c),
                i = paste(
                    "Check the function documentation",
                    "to know more on how to set variables",
                    "correctly"
                )
            )
            rlang::abort(single_err, class = "tag_dupl_err")
        }
        warnings <- purrr::map(dupl_eval, ~ {
            if (.x$type == "MOD") {
                return(.x$warning)
            } else {
                NULL
            }
        })
        if (any(purrr::map_lgl(warnings, ~ !is.null(.x)))) {
            single_warn <- c("Warning: duplicates found for some tags",
                purrr::reduce(warnings, c),
                i = paste(
                    "Check the function documentation",
                    "to know more on how to set variables",
                    "correctly"
                )
            )
            rlang::inform(single_warn, class = "tag_dupl_warn")
        }
        out_df <- purrr::map(dupl_eval, ~ .x$result) %>%
            purrr::reduce(dplyr::bind_rows)
        return(out_df)
    }

    return(df)
}

# Types is a named vec in the form of c("tag" = "col_type")
.eval_tag_types <- function(df, types) {
    tag_split <- df %>%
        dplyr::group_by(.data$tag) %>%
        dplyr::group_split()
    single_type_eval <- function(sub_df) {
        tag_type <- types[[sub_df$tag[1]]]
        if (is.null(tag_type) || purrr::is_empty(tag_type)) {
            return(list(type = "DF", result = sub_df))
        }
        out <- sub_df %>%
            dplyr::filter(.data$types %in% tag_type)
        if (nrow(out) == 0) {
            err <- c(paste("Wrong col class for tag '", sub_df$tag[1], "'"),
                x = paste0(
                    "Expected: ",
                    paste0(tag_type, collapse = ", "),
                    " - Found: ",
                    paste0(sub_df$types, collapse = ", ")
                )
            )
            return(list(type = "ERROR", result = err))
        }
        return(list(type = "DF", result = out))
    }
    result <- purrr::map(tag_split, single_type_eval)
    errors <- purrr::map(result, ~ {
        if (.x$type == "ERROR") {
            return(.x$result)
        }
        return(NULL)
    })
    if (any(purrr::map_lgl(errors, ~ !is.null(.x)))) {
        compact_msg <- c(
            "Wrong column classes for some tags",
            purrr::reduce(errors, c)
        )
        rlang::abort(compact_msg, class = "tag_type_err")
    }
    result <- purrr::map(result, ~ .x$result) %>%
        purrr::reduce(dplyr::bind_rows)
    return(result)
}

#### ---- Internals for checks on integration matrices----####

# Internal helper function for checking mandatory vars presence in x.
#
# Checks if the elements of `mandatory_IS_vars` are present as column names
# in the data frame.
# @param x A data.frame object (or any extending class)
# @keywords internal
#
# @return FALSE if all or some elements are not found in the data frame, TRUE
# otherwise
.check_mandatory_vars <- function(x) {
    stopifnot(is.data.frame(x))
    res <- if (all(mandatory_IS_vars() %in% colnames(x))) {
        TRUE
    } else {
        FALSE
    }
    return(res)
}


# Internal helper function for checking that the pcr replicate
# column presence in x.
#
# @param x A data.frame object (or any extending class)
# @keywords internal
#
# @return FALSE if not found, TRUE otherwise
.check_sample_col <- function(x) {
    stopifnot(is.data.frame(x))
    if (pcr_id_column() %in% colnames(x)) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

# Finds experimental columns in an integration matrix.
#
# The function checks if there are numeric columns which are not
# standard integration matrix columns, if there are returns their names.
#
# @param x A data.frame
#' @importFrom purrr map_lgl
#' @importFrom rlang expr eval_tidy
#' @keywords internal
#
# @return A character vector of column names
.find_exp_cols <- function(x, to_exclude) {
    stopifnot(is.data.frame(x))
    remaining <- colnames(x)[!colnames(x) %in% to_exclude]
    remaining_numeric <- purrr::map_lgl(remaining, function(y) {
        exp <- rlang::expr(`$`(x, !!y))
        exp <- rlang::eval_tidy(exp)
        is.numeric(rlang::eval_tidy(exp)) | is.integer(rlang::eval_tidy(exp))
    })
    remaining[remaining_numeric]
}

# Checks if the integration matrix is annotated or not.
#
# @param x A data.frame
# @keywords internal
#
# @return A logical value
.is_annotated <- function(x) {
    stopifnot(is.data.frame(x))
    if (all(annotation_IS_vars() %in% colnames(x))) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

#### ---- Internals for matrix import ----####

#---- USED IN : import_single_Vispa2Matrix ----

# Reads an integration matrix using data.table::fread
#' @importFrom data.table fread
.read_with_fread <- function(path, additional_cols, annotated, sep) {
    col_types <- .mandatory_IS_types("fread")
    if (annotated) {
        annot_types <- .annotation_IS_types("fread")
        col_types <- purrr::map2(
            col_types, names(col_types),
            ~ c(.x, annot_types[[.y]])
        )
        add_types <- annot_types[!names(col_types) %in% names(annot_types)]
        if (!all(purrr::map_lgl(add_types, is.null))) {
            col_types <- append(col_types, add_types)
        }
    }
    to_drop <- NULL
    if (!is.null(additional_cols) && !purrr::is_empty(additional_cols)) {
        to_drop <- names(additional_cols[additional_cols == "_"])
        others <- additional_cols[!names(additional_cols) %in% to_drop]
        others <- purrr::map(others, ~ {
            .types_mapping() %>%
                dplyr::filter(.data$types == .x) %>%
                dplyr::pull(.data$fread)
        })
        col_types <- purrr::map2(
            col_types, names(col_types),
            ~ c(.x, names(others[others == .y]))
        )
    }

    tmp <- data.table::fread(
        file = path,
        sep = sep,
        na.strings = c("NONE", "NA", "NULL", "NaN", ""),
        verbose = FALSE,
        drop = to_drop,
        colClasses = col_types,
        showProgress = getOption("ISAnalytics.verbose", TRUE),
        data.table = TRUE
    )

    vars <- mandatory_IS_vars(TRUE)
    if (annotated) {
        vars <- vars %>%
            dplyr::bind_rows(annotation_IS_vars(TRUE))
    }
    # Check if there are dates to convert
    dates <- vars %>%
        dplyr::filter(.data$types %in% c(date_formats(), "date"))
    if (nrow(dates) > 0) {
        as_list <- setNames(dates$types, dates$names)
        for (col in names(as_list)) {
            tmp <- tmp %>%
                dplyr::mutate(!!col := .convert_dates(
                    .data[[col]],
                    as_list[[col]]
                ))
        }
    }
    ## Apply transformations if present
    tmp <- .apply_col_transform(tmp, vars)
    tmp <- data.table::setDT(tmp)
    return(tmp)
}

# Reads an integration matrix using readr::read_delim
#' @importFrom readr read_delim cols
#' @importFrom data.table setDT
.read_with_readr <- function(path, additional_cols, annotated, sep) {
    col_types <- .mandatory_IS_types("classic")
    if (annotated) {
        col_types <- append(
            col_types,
            .annotation_IS_types("classic")
        )
    }
    if (!is.null(additional_cols) && !purrr::is_empty(additional_cols)) {
        to_drop <- additional_cols[additional_cols == "_"]
        others <- additional_cols[!names(additional_cols) %in% names(to_drop)]
        others <- purrr::map(others, ~ {
            .types_mapping() %>%
                dplyr::filter(.data$types == .x) %>%
                dplyr::pull(.data$mapping)
        })
        for (x in names(to_drop)) {
            col_types[[x]] <- "_"
        }
        for (x in names(others)) {
            col_types[[x]] <- others[[x]]
        }
    }
    col_types[[".default"]] <- "n"
    df <- readr::read_delim(
        file = path,
        delim = sep,
        col_types = do.call(readr::cols, col_types),
        na = c("NONE", "NA", "NULL", "NaN", ""),
        trim_ws = TRUE,
        progress = getOption("ISAnalytics.verbose", TRUE)
    )
    vars <- mandatory_IS_vars(TRUE)
    if (annotated) {
        vars <- vars %>%
            dplyr::bind_rows(annotation_IS_vars(TRUE))
    }
    # Check if there are dates to convert
    dates <- vars %>%
        dplyr::filter(.data$types %in% c(date_formats(), "date"))
    if (nrow(dates) > 0) {
        as_list <- setNames(dates$types, dates$names)
        for (col in names(as_list)) {
            df <- df %>%
                dplyr::mutate(!!col := .convert_dates(
                    .data[[col]],
                    as_list[[col]]
                ))
        }
    }
    ## Apply transformations if present
    df <- .apply_col_transform(df, vars)
    df <- data.table::setDT(df)
    return(df)
}

# - call_mode: either INTERNAL or EXTERNAL. First is used when calling
#              the function
#              from parallel importing function, second one when calling
#              import_single_Vispa2_matrix
.import_single_matrix <- function(path,
    separator = "\t",
    additional_cols = NULL,
    transformations = NULL,
    call_mode = "INTERNAL",
    id_col_name = pcr_id_column(),
    val_col_name = "Value") {
    ## - Check additional cols specs: must be named vector or list or NULL
    stopifnot(is.null(additional_cols) ||
        (is.character(additional_cols) &
            !is.null(names(additional_cols))) ||
        (is.list(additional_cols) & !is.null(names(additional_cols))))
    correct_types <- additional_cols %in% c(.types_mapping()$types, "_")
    if (any(correct_types == FALSE)) {
        add_types_err <- c("Unknown column type in specified additional columns",
            x = paste("Types must be in the allowed types"),
            i = paste(
                "Unknown formats: ",
                paste0(unique(
                    additional_cols[correct_types == FALSE]
                ), collapse = ", ")
            ),
            i = paste(
                "See documentation for details",
                "?import_single_Vispa2Matrix"
            )
        )
        rlang::abort(add_types_err, class = "add_types_err")
    }
    mode <- "fread"
    ## Is the file compressed?
    is_compressed <- fs::path_ext(path) %in% .compressed_formats()
    if (is_compressed) {
        ## The compression type is supported by data.table::fread?
        compression_type <- fs::path_ext(path)
        if (!compression_type %in% .supported_fread_compression_formats()) {
            ### If not, switch to classic for reading
            mode <- "classic"
            if (call_mode == "EXTERNAL" &&
                getOption("ISAnalytics.verbose", TRUE) == TRUE) {
                rlang::inform(.unsupported_comp_format_inf(),
                    class = "unsup_comp_format"
                )
            }
        }
    }
    ### Peak headers
    peek_headers <- readr::read_delim(path,
        delim = separator, n_max = 0,
        col_types = readr::cols()
    )
    ## - Mandatory vars should always be present (as declared in specs)
    if (!all(mandatory_IS_vars() %in% colnames(peek_headers))) {
        missing_mand <- mandatory_IS_vars()[!mandatory_IS_vars() %in%
            colnames(peek_headers)]
        rlang::abort(.missing_req_cols(mandatory_IS_vars(), missing_mand),
            class = "im_single_miss_mand_vars"
        )
    }
    is_annotated <- .is_annotated(peek_headers)
    additional_cols <- additional_cols[names(additional_cols) %in%
        colnames(peek_headers)]
    ## - Start reading
    if (call_mode == "EXTERNAL" &&
        getOption("ISAnalytics.verbose", TRUE) == TRUE) {
        rlang::inform(c("Reading file...", i = paste0("Mode: ", mode)))
    }
    df <- if (mode == "fread") {
        .read_with_fread(
            path = path, additional_cols = additional_cols,
            annotated = is_annotated, sep = separator
        )
    } else {
        .read_with_readr(
            path = path, additional_cols = additional_cols,
            annotated = is_annotated, sep = separator
        )
    }
    df_dim <- dim(df)
    ## - Melt
    id_vars <- c(
        mandatory_IS_vars(),
        names(additional_cols[additional_cols != "_"])
    )
    if (is_annotated) {
        id_vars <- c(id_vars, annotation_IS_vars())
    }
    mt <- function(data, id_vars, sample_col, value_col) {
        data.table::melt.data.table(data,
            id.vars = id_vars,
            variable.name = sample_col,
            value.name = value_col,
            na.rm = TRUE,
            verbose = FALSE
        )
    }
    if (call_mode == "EXTERNAL" &&
        getOption("ISAnalytics.verbose", TRUE) == TRUE) {
        rlang::inform("Reshaping...")
    }
    max_workers <- trunc(BiocParallel::snowWorkers() / 3)
    if (call_mode == "INTERNAL" || nrow(df) <= 500000 || max_workers <= 1) {
        tidy <- mt(df, id_vars, id_col_name, val_col_name)
    } else if (nrow(df) > 500000 & max_workers > 1) {
        ## - Melt in parallel
        p <- if (.Platform$OS.type == "windows") {
            BiocParallel::SnowParam(
                workers = max_workers,
                tasks = trunc(nrow(df) / max_workers),
                progressbar = getOption("ISAnalytics.verbose", TRUE),
                exportglobals = TRUE,
                stop.on.error = TRUE
            )
        } else {
            BiocParallel::MulticoreParam(
                workers = max_workers,
                tasks = trunc(nrow(df) / max_workers),
                progressbar = getOption("ISAnalytics.verbose", TRUE),
                exportglobals = FALSE,
                stop.on.error = TRUE
            )
        }
        ## - Split in chunks
        chunk_id_vec <- rep(seq_len(max_workers - 1),
            each = trunc(nrow(df) / max_workers)
        )
        chunk_id_vec <- c(
            chunk_id_vec,
            rep_len(chunk_id_vec[length(chunk_id_vec)] + 1,
                length.out = nrow(df) - length(chunk_id_vec)
            )
        )
        df[, c("chunk_id") := chunk_id_vec]
        chunks <- split(df,
            by = c("chunk_id"),
            verbose = FALSE,
            keep.by = FALSE
        )
        tidy_chunks <- BiocParallel::bplapply(
            X = chunks,
            FUN = mt,
            id_vars = id_vars,
            sample_col = id_col_name,
            value_col = val_col_name,
            BPPARAM = p
        )
        BiocParallel::bpstop(p)
        tidy <- data.table::rbindlist(tidy_chunks)
    }
    tidy <- tidy[val_col_name > 0]
    ## Transform cols
    if (call_mode == "EXTERNAL" &&
        getOption("ISAnalytics.verbose", TRUE) == TRUE &&
        !is.null(transformations)) {
        rlang::inform("Applying column transformations...")
    }
    if (!is.null(transformations)) {
        tidy <- transform_columns(tidy, transf_list = transformations)
    }
    ## - Report summary
    if (call_mode == "EXTERNAL" &&
        getOption("ISAnalytics.verbose", TRUE) == TRUE) {
        rlang::inform(.summary_ism_import_msg(
            is_annotated,
            df_dim,
            mode,
            length(unique(tidy[[id_col_name]]))
        ),
        class = "ism_import_summary"
        )
    }
    return(tidy)
}

#---- USED IN : import_association_file ----

# Checks if the association file contains at least the the columns flagged
# as required in `association_columns(TRUE)`
#
# @param df The imported association file
.check_af_correctness <- function(df) {
    af_cols_specs <- association_file_columns(TRUE)
    ### Retrieves the column names flagged as required
    min_required <- af_cols_specs %>%
        dplyr::filter(.data$flag == "required") %>%
        dplyr::pull(.data$names)
    if (length(min_required) == 0) {
        return(TRUE) # No required columns
    }
    if (all(min_required %in% colnames(df))) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

# Imports association file from disk. Converts dates and pads timepoints,
# reporting parsing problems.
.read_af <- function(path, date_format, delimiter) {
    mode <- "readr"
    ## - Check file extension
    file_ext <- .check_file_extension(path)
    if (file_ext %in% c("xls", "xlsx")) {
        mode <- "readxl"
        if (!requireNamespace("readxl", quietly = TRUE)) {
            rlang::abort(.missing_pkg_error("readxl"))
        }
        if (getOption("ISAnalytics.verbose", TRUE) == TRUE) {
            rlang::inform(.xls_file_warning(),
                class = "xls_file"
            )
        }
    }
    ## - Peek headers to know column types
    col_types <- if (mode != "readxl") {
        headers_peek <- readr::read_delim(path,
            delim = delimiter,
            n_max = 0,
            col_types = readr::cols()
        )
        t_readr <- .af_col_types("readr")
        t_readr[names(t_readr) %in% colnames(headers_peek)]
    } else {
        headers_peek <- readxl::read_excel(path, n_max = 0)
        t_readr <- .af_col_types("readr")
        ordered <- purrr::map_chr(
            colnames(headers_peek),
            function(x) {
                el <- getElement(t_readr, x)
                el <- if (is.null(el)) {
                    "guess"
                } else if (el == "c") {
                    "text"
                } else if (el %in% c("i", "d")) {
                    "numeric"
                } else {
                    "guess"
                }
            }
        )
        ordered
    }
    suppressWarnings({
        as_file <- if (mode == "readr") {
            df <- readr::read_delim(path,
                delim = delimiter,
                na = c("NONE", "NA", "NULL", "NaN", ""),
                col_types = do.call(readr::cols, col_types),
                trim_ws = TRUE,
                progress = getOption("ISAnalytics.verbose", TRUE)
            )
            problems <- readr::problems(df)
            df
        } else {
            df <- readxl::read_excel(path,
                col_types = col_types,
                na = c("NONE", "NA", "NULL", "NaN", ""),
                trim_ws = TRUE,
                progress = getOption("ISAnalytics.verbose", TRUE)
            )
            problems <- NULL
            df
        }

        vars <- association_file_columns(TRUE)
        # Check if there are dates to convert
        dates <- vars %>%
            dplyr::filter(.data$types %in% c(date_formats(), "date")) %>%
            dplyr::mutate(types = dplyr::if_else(.data$types == "date",
                true = date_format,
                false = .data$types
            ))
        dates <- dates %>%
            dplyr::filter(.data$names %in% colnames(as_file))
        date_failures <- NULL
        if (nrow(dates) > 0) {
            before <- as_file %>%
                dplyr::select(dplyr::all_of(dates$names))
            as_list <- setNames(dates$types, dates$names)
            for (col in names(as_list)) {
                as_file <- as_file %>%
                    dplyr::mutate(!!col := .convert_dates(
                        .data[[col]],
                        as_list[[col]]
                    ))
            }
            date_failures <- purrr::map_dfr(dates$names, function(col) {
                before_col <- purrr::pluck(before, col)
                row_failed <- which(
                    purrr::map2_lgl(
                        before_col, as_file[[col]],
                        function(d1, d2) {
                            if (!is.na(d1) & is.na(d2)) {
                                TRUE
                            } else {
                                FALSE
                            }
                        }
                    )
                )
                if (!is.null(row_failed) && !purrr::is_empty(row_failed)) {
                    tibble::tibble(
                        row = row_failed, col = col,
                        original = before_col[row_failed],
                        parsed = NA, format = date_format
                    )
                } else {
                    return(NULL)
                }
            })
        }
        ## Apply transformations if present
        as_file <- .apply_col_transform(as_file, vars)
        as_file <- data.table::setDT(as_file)
    })
    return(list(af = as_file, probs = problems, date_fail = date_failures))
}

# Internal function to check alignment between association file and file
# system starting from the root. The alignment is checked at folder level.
#
# - df: The imported association file (data.frame or tibble)
# - root_folder: Path to the root folder
#' @importFrom rlang .data `:=`
# Returns a data frame with this format:
# ProjectID - ConcatenatePoolIDSeqRun - PathToFolderProjectID - Found - Path -
# Path_quant - Path_iss (NOTE: headers are dynamic!)
.check_file_system_alignment <- function(df,
    root_folder,
    proj_fold_col,
    concat_pool_col,
    project_id_col) {
    if (!proj_fold_col %in% colnames(df)) {
        rlang::abort(.af_missing_pathfolder_error(proj_fold_col))
    }
    temp_df <- df %>%
        dplyr::select(dplyr::all_of(c(
            project_id_col,
            concat_pool_col,
            proj_fold_col
        ))) %>%
        dplyr::distinct()
    path_cols <- .path_cols_names()
    proj_folders_exist <- temp_df %>%
        dplyr::select(.data[[proj_fold_col]]) %>%
        dplyr::distinct() %>%
        dplyr::mutate(Found = !is.na(
            .data[[proj_fold_col]]
        ) &
            unname(fs::dir_exists(
                fs::path(
                    fs::path(root_folder),
                    .data[[proj_fold_col]]
                )
            )))
    partial_check <- temp_df %>%
        dplyr::left_join(proj_folders_exist, by = proj_fold_col)
    temp_df <- partial_check %>%
        dplyr::filter(.data$Found == TRUE)
    partial_check <- partial_check %>%
        dplyr::filter(.data$Found == FALSE) %>%
        dplyr::mutate(
            !!path_cols$project := NA_character_,
            !!path_cols$quant := NA_character_,
            !!path_cols$iss := NA_character_
        )
    FUN <- function(...) {
        cur <- tibble::tibble(...)
        project_folder <- fs::path(
            fs::path(root_folder),
            cur[[proj_fold_col]]
        )
        quant_folder <- if (!is.na(cur[[concat_pool_col]])) {
            paste0(fs::path(
                "quantification",
                fs::path(cur[[concat_pool_col]])
            ), "$")
        } else {
            if (getOption("ISAnalytics.verbose", TRUE) == TRUE) {
                msg <- c(paste(
                    "Warning: found NA",
                    paste0(concat_pool_col, " field")
                ))
                rlang::inform(msg,
                    i = "Check association file for possible issues",
                    class = "na_concat"
                )
            }
            NA_character_
        }
        iss_folder <- if (!is.na(cur[[concat_pool_col]])) {
            paste0(fs::path(
                "iss",
                fs::path(cur[[concat_pool_col]])
            ), "$")
        } else {
            NA_character_
        }
        quant_found <- unique(fs::dir_ls(
            path = project_folder, recurse = TRUE,
            type = "directory", fail = FALSE,
            regexp = quant_folder
        ))
        if (length(quant_found) == 0 || all(is.na(quant_found))) {
            quant_found <- NA_character_
        }
        iss_found <- unique(fs::dir_ls(
            path = project_folder, recurse = TRUE,
            type = "directory", fail = FALSE,
            regexp = iss_folder
        ))
        if (length(iss_found) == 0 || all(is.na(iss_found))) {
            iss_found <- NA_character_
        }
        return(
            cur %>%
                dplyr::mutate(
                    !!path_cols$project := project_folder,
                    !!path_cols$quant := quant_found,
                    !!path_cols$iss := iss_found
                )
        )
    }
    results_df <- purrr::pmap_dfr(
        temp_df,
        FUN
    )
    checker_df <- results_df %>%
        dplyr::bind_rows(partial_check)
    checker_df
}

# Updates the association file after the alignment check --> adds the columns
# containing paths to project folder, quant folders and iss folders
.update_af_after_alignment <- function(as_file, checks,
    proj_fold_col,
    concat_pool_col,
    project_id_col) {
    as_file <- as_file %>%
        dplyr::left_join(checks %>%
            dplyr::select(
                -.data[[proj_fold_col]],
                -.data$Found
            ),
        by = c(project_id_col, concat_pool_col)
        )
    as_file
}

# Helper function to be used internally to treat association file.
.manage_association_file <- function(af_path,
    root,
    format,
    delimiter,
    filter,
    proj_fold_col,
    concat_pool_col,
    project_id_col) {
    # Import the association file
    association_file <- .read_af(
        path = af_path,
        date_format = format,
        delimiter = delimiter
    )
    parsing_probs <- association_file$probs
    date_probs <- association_file$date_fail
    association_file <- association_file$af
    if (!is.null(filter)) {
        # Pre-filtering of association file
        if (!all(names(filter) %in% colnames(association_file))) {
            if (getOption("ISAnalytics.verbose", TRUE) == TRUE) {
                missing_filt <- c(
                    "Some or all the names in the filter not found",
                    i = paste(
                        "Columns not found: ",
                        paste0(
                            names(filter)[names(filter) %in%
                                colnames(association_file)],
                            collapse = ", "
                        )
                    ),
                    "Ignoring the missing columns"
                )
                rlang::inform(
                    missing_filt,
                    class = "filter_warn"
                )
            }
            filter <- filter[names(filter) %in% colnames(association_file)]
        }
        if (!purrr::is_empty(filter)) {
            predicate <- purrr::map2_chr(
                filter, names(filter),
                function(x, y) {
                    if (is.character(x)) {
                        paste0(
                            y, " %in% c(", paste0("'",
                                rlang::enexpr(x),
                                "'",
                                collapse = ", "
                            ),
                            ")"
                        )
                    } else {
                        paste0(
                            y, " %in% c(", paste0(rlang::enexpr(x),
                                collapse = ", "
                            ),
                            ")"
                        )
                    }
                }
            )
            predicate <- paste0(c(
                "dplyr::filter(association_file, ",
                paste0(predicate, collapse = ", "),
                ")"
            ), collapse = "")
            association_file <- rlang::eval_tidy(
                rlang::parse_expr(predicate)
            )
        }
    }
    checks <- NULL
    if (!is.null(root) & nrow(association_file) > 0) {
        checks <- .check_file_system_alignment(
            df = association_file,
            root_folder = root,
            proj_fold_col = proj_fold_col,
            concat_pool_col = concat_pool_col,
            project_id_col = project_id_col
        )
        association_file <- .update_af_after_alignment(
            as_file = association_file,
            checks = checks,
            proj_fold_col = proj_fold_col,
            concat_pool_col = concat_pool_col,
            project_id_col = project_id_col
        )
    }
    res <- list(
        af = association_file, check = checks,
        parsing_probs = parsing_probs, date_probs = date_probs
    )
    return(res)
}

# Converts a vector of timepoints expressed as days to a vector
# of timepoints expressed as months with this logic:
# - Month is 0 if and only if timepoint is 0
# - If timepoint is strictly less than 30 returns the smallest integer
#   NOT LESS THAN tp / 30
#   ---> ex: 29 becomes ceiling(29/30) = ceiling(0.9666667) = 1
# - If timepoint is greater than or equal to 30 returns the rounding of
#   tp / 30 to the closest integer
#   ---> ex: 35 becomes round(35/30) = round(1.166667) = 1
# - If for any reason the timepoint is negative NA is returned instead
.timepoint_to_months <- function(tp) {
    dplyr::if_else(
        condition = as.numeric(tp) > 0,
        true = dplyr::if_else(
            condition = tp < 30,
            true = ceiling(as.numeric(tp) / 30),
            false = round(as.numeric(tp) / 30)
        ),
        false = dplyr::if_else(
            condition = as.numeric(tp) != 0,
            true = NA_real_,
            false = 0
        )
    )
}

# Converts a vector of timepoints expressed as days to a vector
# of timepoints expressed as years with this logic:
# - Year is 0 if and only if tp is 0
# - otherwise take the ceiling of tp / 360
.timepoint_to_years <- function(tp) {
    dplyr::if_else(
        condition = as.numeric(tp) == 0,
        true = 0,
        false = dplyr::if_else(
            condition = as.numeric(tp) < 0,
            true = NA_real_,
            false = ceiling(as.numeric(tp) / 360)
        )
    )
}

#---- USED IN : import_Vispa2_stats ----

# Finds automatically the path on disk to each stats file.
#' @importFrom purrr pmap_dfr detect_index
#' @importFrom tibble tibble
#' @importFrom fs dir_ls
#' @importFrom rlang .data
# @keywords internal
# @return A tibble with columns: ProjectID, concatenatePoolIDSeqRun,
# Path_iss (or designated dynamic name), stats_files, info
.stats_report <- function(association_file,
    prefixes,
    proj_col,
    pool_col,
    path_iss_col) {
    temp <- association_file %>%
        dplyr::select(
            dplyr::all_of(c(proj_col, pool_col, path_iss_col))
        ) %>%
        dplyr::distinct() %>%
        dplyr::mutate(
            !!proj_col := as.factor(.data[[proj_col]]),
            !!pool_col := as.factor(.data[[pool_col]])
        )
    # If paths are all NA return
    if (all(is.na(temp[[path_iss_col]]))) {
        return(temp %>% dplyr::mutate(
            stats_files = NA_character_,
            info = list("NO FOLDER FOUND")
        ))
    }
    match_pattern <- function(pattern, temp_row) {
        if (is.na(temp_row[[path_iss_col]])) {
            return(tibble::tibble(pattern = pattern, file = NA_character_))
        }
        # For each prefix pattern search in the iss folder
        # Note: there can be
        # - a single file matching the prefix (ideal)
        # - multiple files matching
        # - no file matching
        files <- fs::dir_ls(temp_row[[path_iss_col]],
            type = "file", fail = FALSE,
            regexp = pattern
        )
        if (length(files) == 0) {
            files <- NA_character_
        }
        tibble::tibble(pattern = pattern, file = files)
    }
    stats_paths <- purrr::pmap_dfr(temp, function(...) {
        temp_row <- tibble::tibble(...)
        matches <- purrr::map_dfr(prefixes, ~ match_pattern(.x, temp_row))
        if (all(is.na(matches$file))) {
            # No stats files found for all prefixes
            return(temp_row %>%
                dplyr::mutate(
                    stats_files = NA_character_,
                    info = list("NOT FOUND")
                ))
        }
        if (length(which(!is.na(matches$file))) == 1) {
            return(temp_row %>%
                dplyr::mutate(
                    stats_files = (matches %>%
                        dplyr::filter(!is.na(.data$file)) %>%
                        dplyr::pull(.data$file)),
                    info = list("NULL")
                ))
        }
        # Get the count of files found for each pattern and order them
        # as preference - preferred pattern is the first in the prefixes
        # parameter
        pattern_counts <- matches %>%
            dplyr::mutate(pattern = factor(.data$pattern,
                levels = prefixes
            )) %>%
            dplyr::filter(!is.na(.data$file)) %>%
            dplyr::group_by(.data$pattern) %>%
            dplyr::tally() %>%
            dplyr::arrange(.data$pattern)
        one_index <- purrr::detect_index(pattern_counts$n, ~ .x == 1)
        if (one_index == 0) {
            ## Means there are duplicates for all patterns, gets reported
            ## but no files will be imported
            return(temp_row %>%
                dplyr::mutate(
                    stats_files = NA_character_,
                    info = list(matches %>%
                        dplyr::filter(!is.na(.data$file)) %>%
                        dplyr::mutate(info = "DUPLICATE"))
                ))
        } else {
            ## Pick the file with the pattern that has a count of 1 (no
            ## duplicates)
            chosen_pattern <- pattern_counts$pattern[one_index]
            return(temp_row %>%
                dplyr::mutate(
                    stats_files = matches %>%
                        dplyr::filter(.data$pattern == chosen_pattern) %>%
                        dplyr::pull(.data$file),
                    info = list("NULL")
                ))
        }
    })
    stats_paths
}


# Checks if the stats file contains the minimal set of columns.
.check_stats <- function(x, required_cols) {
    if (all(required_cols %in% colnames(x))) {
        TRUE
    } else {
        FALSE
    }
}

# Imports all found Vispa2 stats files.
#
# - association_file: imported and aligned association file
# - prefixes: vector of regex to match file names
# - pool_col: the column containing the pool as specified in input
# - path_iss_col: name of the column that contains the path
# - tags: injected from parent function - collection of required tags and
#   corresponding column names
#
# @return A list with the imported stats and a report of imported files. If no
# files were imported returns NULL instead
.import_stats_iss <- function(association_file,
    prefixes,
    pool_col,
    path_iss_col,
    tags) {
    proj_col <- tags %>%
        dplyr::filter(.data$tag == "project_id") %>%
        dplyr::pull(.data$names)
    # Obtain paths
    stats_paths <- .stats_report(association_file,
        prefixes,
        proj_col = proj_col,
        pool_col = pool_col,
        path_iss_col = path_iss_col
    )
    stats_paths <- stats_paths %>%
        tidyr::unnest(.data$info)
    if (all(is.na(stats_paths$stats_files))) {
        stats_paths <- stats_paths %>%
            dplyr::mutate(
                Imported = FALSE,
                reason = as.factor(.data$info)
            ) %>%
            dplyr::mutate(-.data$info)
        return(list(stats = NULL, report = stats_paths))
    }
    vispa_stats_req <- iss_stats_specs(TRUE) %>%
        dplyr::filter(.data$flag == "required") %>%
        dplyr::pull(.data$names)
    # Setup parallel workers and import
    # Register backend according to platform
    if (.Platform$OS.type == "windows") {
        p <- BiocParallel::SnowParam(
            stop.on.error = FALSE,
            tasks = length(stats_paths$stats_files),
            progressbar = getOption("ISAnalytics.verbose", TRUE)
        )
    } else {
        p <- BiocParallel::MulticoreParam(
            stop.on.error = FALSE,
            tasks = length(stats_paths$stats_files),
            progressbar = getOption("ISAnalytics.verbose", TRUE)
        )
    }
    FUN <- function(x, req_cols) {
        if (is.na(x)) {
            return(NULL)
        }
        stats <- data.table::fread(
            file = x, sep = "\t",
            na.strings = c("", "NA", "na", "NONE"),
            data.table = TRUE
        )
        ok <- .check_stats(stats, req_cols)
        if (ok == TRUE) {
            # Apply column transformations if present
            stats <- .apply_col_transform(stats, iss_stats_specs(TRUE))
            return(stats)
        } else {
            return("MALFORMED")
        }
    }
    stats_dfs <- BiocParallel::bptry(
        BiocParallel::bplapply(stats_paths$stats_files,
            FUN,
            req_cols = vispa_stats_req,
            BPPARAM = p
        )
    )
    BiocParallel::bpstop(p)
    correct <- purrr::map_chr(stats_dfs, function(x) {
        if (all(is.data.frame(x))) {
            "TRUE"
        } else if (is.null(x)) {
            "FALSE"
        } else {
            x
        }
    })
    stats_paths <- stats_paths %>%
        dplyr::mutate(Imported = correct)
    stats_paths <- purrr::pmap_dfr(stats_paths, function(...) {
        row <- tibble::tibble(...)
        condition <- if (row$Imported == "MALFORMED") {
            "MALFORMED"
        } else if (row$Imported == "TRUE") {
            NA_character_
        } else {
            row$info
        }
        row %>%
            dplyr::mutate(reason = as.factor(condition)) %>%
            dplyr::mutate(Imported = dplyr::if_else(
                condition = (.data$Imported == "MALFORMED"),
                true = FALSE,
                false = as.logical(.data$Imported)
            ))
    })
    stats_paths <- stats_paths %>%
        dplyr::select(-.data$info)
    stats_dfs <- stats_dfs[stats_paths$Imported]
    # Bind rows in single tibble for all files
    if (purrr::is_empty(stats_dfs)) {
        return(list(stats = NULL, report = stats_paths))
    }
    stats_dfs <- purrr::reduce(stats_dfs, function(x, y) {
        x %>%
            dplyr::bind_rows(y) %>%
            dplyr::distinct()
    })
    list(stats = stats_dfs, report = stats_paths)
}

#---- USED IN : import_parallel_Vispa2Matrices ----
# Base arg check (valid for both versions)
.base_param_check <- function(association_file,
    quantification_type,
    matrix_type,
    workers,
    multi_quant_matrix) {
    # Check parameters
    stopifnot((is.character(association_file) &
        length(association_file) == 1) ||
        is.data.frame(association_file))
    stopifnot(is.numeric(workers) & length(workers) == 1)
    stopifnot(!missing(quantification_type))
    stopifnot(all(quantification_type %in% quantification_types()))
    stopifnot(is.logical(multi_quant_matrix) & length(multi_quant_matrix) == 1)
}

# Import AF if necessary, if already imported check it is aligned with
# file system
.pre_manage_af <- function(association_file, import_af_args, report_path) {
    if (!is.null(report_path) && !fs::is_dir(report_path)) {
        report_path <- fs::path_dir(report_path)
    }
    ## Import association file if provided a path
    if (is.character(association_file)) {
        association_file <- rlang::eval_tidy(
            rlang::call2("import_association_file",
                path = association_file,
                report_path = report_path,
                !!!import_af_args
            )
        )
    }
    ## Check there are the appropriate columns
    if (!.path_cols_names()$quant %in% colnames(association_file)) {
        rlang::abort(.af_missing_path_error(.path_cols_names()$quant),
            class = "missing_path_col"
        )
    }
    association_file <- association_file %>%
        dplyr::filter(!is.na(.data[[.path_cols_names()$quant]]))
    return(association_file)
}


# Allows the user to choose interactively the projects
# to consider for import.
#
# @param association_file The tibble representing the imported association file
# @keywords internal
#' @importFrom dplyr distinct select filter
#' @importFrom rlang .data
#' @importFrom stringr str_split
#' @importFrom purrr map_dbl
#
# @return A modified version of the association file where only selected
# projects are present
.interactive_select_projects_import <- function(association_file,
    proj_col) {
    repeat {
        cat("Which projects would you like to import?\n")
        cat("[1] ALL", "\n", "[2] ONLY SOME\n", "[0] QUIT\n", sep = "")
        # Keep asking user until the choice is valid
        repeat {
            cat("Your choice: ")
            connection <- getOption("ISAnalytics.connection", stdin())
            if (length(connection) > 1) {
                con <- connection[1]
            } else {
                con <- connection
            }
            n_projects_to_import <- readLines(con = con, n = 1)
            n_projects_to_import <- as.numeric(n_projects_to_import)
            if (!is.na(n_projects_to_import) &
                is.numeric(n_projects_to_import) &
                n_projects_to_import %in% c(0, 1, 2)) {
                if (n_projects_to_import == 0) {
                    rlang::abort("Quitting")
                } else {
                    break
                }
            } else {
                rlang::inform(paste(
                    "Invalid choice,",
                    "please choose between the above."
                ))
            }
        }
        # If only some projects selected
        if (n_projects_to_import == 2) {
            cat(
                "Here is the list of available projects, type the indexes",
                "separated by a comma:\n",
                sep = ""
            )
            project_list <- association_file %>%
                dplyr::distinct(.data[[proj_col]]) %>%
                dplyr::pull(.data[[proj_col]])
            cat(paste0("[", seq_along(project_list), "] ", project_list),
                "[0] QUIT\n",
                sep = "\n"
            )
            # Keep asking user until a valid choice
            repeat {
                cat("Your choice: ")
                if (length(connection) > 1) {
                    con <- connection[2]
                } else {
                    con <- connection
                }
                projects_to_import <- readLines(con = con, n = 1)
                projects_to_import <- purrr::map_dbl(unlist(
                    stringr::str_split(projects_to_import, ",")
                ), as.numeric)
                if (!any(is.na(projects_to_import)) &
                    all(is.numeric(projects_to_import)) &
                    all(projects_to_import %in% c(
                        seq_along(project_list),
                        0
                    ))) {
                    if (any(projects_to_import == 0)) {
                        stop("Quitting", call. = FALSE)
                    } else {
                        break
                    }
                } else {
                    message(paste(
                        "Invalid choice, please select valid indexes separated",
                        "by a comma (example: 1, 2, 3) "
                    ))
                }
            }
            # Ask confirm
            cat("\nYour choices: ",
                paste(project_list[projects_to_import], collapse = ","),
                "\n",
                "Confirm your choices? [y/n]",
                sep = ""
            )
            repeat {
                if (length(connection) > 1) {
                    con <- connection[3]
                } else {
                    con <- connection
                }
                confirm <- readLines(con = con, n = 1)
                if (!confirm %in% c("y", "n")) {
                    message(paste(
                        "Invalid choice, please select one between y",
                        "or n"
                    ))
                } else {
                    break
                }
            }
            if (confirm == "y") {
                result <- association_file %>%
                    dplyr::filter(.data[[proj_col]] %in%
                        project_list[projects_to_import])
                return(result)
            } else {
                next
            }
        } else {
            # Ask confirm
            cat("\nYour choices: ", "import all projects", "\n",
                "Confirm your choices? [y/n]",
                sep = ""
            )
            repeat {
                if (length(connection) > 1) {
                    con <- connection[3]
                } else {
                    con <- connection
                }
                confirm <- readLines(con = con, n = 1)
                if (!confirm %in% c("y", "n")) {
                    message(paste(
                        "Invalid choice, please select one between y",
                        "or n"
                    ))
                } else {
                    break
                }
            }
            if (confirm == "y") {
                return(association_file)
            } else {
                next
            }
        }
    }
}

# Simple internal helper function to handle user input for selection
# of number of pools.
# @keywords internal
#
# @return Numeric representing user selection (1 for all pools, 2 for
# only some pools, 0 to exit)
.pool_number_IN <- function() {
    cat("Which pools for each project would you like to import?\n")
    cat("[1] ALL", "[2] ONLY SOME", "[0] QUIT", sep = "\n")
    repeat {
        cat("Your choice: ")
        connection <- getOption("ISAnalytics.connection", stdin())
        if (length(connection) > 1) {
            connection <- connection[4]
        }
        n_pools_to_import <- readLines(con = connection, n = 1)
        n_pools_to_import <- as.numeric(n_pools_to_import)
        if (!is.na(n_pools_to_import) & is.numeric(n_pools_to_import) &
            n_pools_to_import %in% c(0, 1, 2)) {
            if (n_pools_to_import == 0) {
                stop("Quitting", call. = FALSE)
            } else {
                break
            }
        } else {
            message("Invalid choice, please choose between the above.")
        }
    }
    n_pools_to_import
}

# Simple helper interal function to handle user input for actual
# pool choices.
#
# @param indexes A vector of integer indexes available
# @keywords internal
#' @importFrom stringr str_split
#' @importFrom purrr map_dbl
#
# @return The user selection as a numeric vector
.pool_choices_IN <- function(indexes) {
    repeat {
        cat("\nYour choice: ")
        connection <- getOption("ISAnalytics.connection", stdin())
        if (length(connection) > 1) {
            connection <- connection[5]
        }
        to_imp <- readLines(
            con = connection,
            n = 1
        )
        to_imp <- purrr::map_dbl(unlist(
            stringr::str_split(to_imp, ",")
        ), as.numeric)
        if (!any(is.na(to_imp)) & all(is.numeric(to_imp))) {
            if (all(to_imp %in%
                c(indexes, 0))) {
                if (any(to_imp == 0)) {
                    stop("Quitting", call. = FALSE)
                } else {
                    break
                }
            } else {
                message(paste(
                    "Invalid choice, please select valid indexes",
                    "separated by a comma (example: 1, 2, 3) "
                ))
                next
            }
        } else {
            message(paste(
                "Invalid choice, please select valid indexes separated",
                "by a comma (example: 1, 2, 3) "
            ))
        }
    }
    to_imp
}

# Allows the user to choose interactively
# the pools to consider for import.
#' @importFrom rlang .data
#'
# @return A modified version of the association file where only selected
# pools for each project are present
.interactive_select_pools_import <- function(association_file,
    proj_col,
    pool_col) {
    repeat {
        n_pools_to_import <- .pool_number_IN()
        if (n_pools_to_import == 2) {
            cat("Here is the list of available pools for each project,",
                "type the indexes separated by a comma or 0 to quit: \n\n",
                sep = ""
            )
            available <- association_file %>%
                dplyr::select(
                    .data[[proj_col]],
                    .data[[pool_col]]
                ) %>%
                dplyr::distinct() %>%
                tidyr::nest(data = c(.data[[pool_col]]))
            pools_to_import <- purrr::pmap(available, function(...) {
                l <- list(...)
                current <- l[proj_col]
                current_data <- tibble::as_tibble(
                    purrr::flatten(l["data"])
                )
                indexes <- seq_along(current_data[[pool_col]])
                cat("Project: ", current[[proj_col]], "\n\n")
                cat("[0] QUIT\n")
                purrr::walk2(
                    current_data[[pool_col]],
                    indexes,
                    function(x, y) {
                        cat(paste0("[", y, "] ", x),
                            sep = "\n"
                        )
                    }
                )
                to_imp <- .pool_choices_IN(indexes)
                tibble::tibble(
                    !!proj_col := current[[proj_col]],
                    !!pool_col := current_data[[pool_col]][to_imp]
                )
            })
            pools_to_import <- purrr::reduce(pools_to_import, function(x, y) {
                dplyr::bind_rows(x, y)
            })
            cat("\nYour choices: ", sep = "\n")
            print(pools_to_import)
            cat("\nConfirm your choices? [y/n]", sep = "")
            repeat {
                connection <- getOption("ISAnalytics.connection", stdin())
                if (length(connection) > 1) {
                    connection <- connection[6]
                }
                confirm <- readLines(con = connection, n = 1)
                if (!confirm %in% c("y", "n")) {
                    message(paste(
                        "Invalid choice, please select one between y",
                        "or n"
                    ))
                } else {
                    break
                }
            }
            if (confirm == "y") {
                association_file <- association_file %>%
                    dplyr::inner_join(pools_to_import,
                        by = c(proj_col, pool_col)
                    )
                return(association_file)
            } else {
                next
            }
        } else {
            cat("\nYour choices: ", sep = "\n")
            print(association_file %>%
                dplyr::select(
                    .data[[proj_col]],
                    .data[[pool_col]]
                ) %>%
                dplyr::distinct())
            cat("\nConfirm your choices? [y/n]", sep = "")
            repeat {
                connection <- getOption("ISAnalytics.connection", stdin())
                if (length(connection) > 1) {
                    connection <- connection[6]
                }
                confirm <- readLines(con = connection, n = 1)
                if (!confirm %in% c("y", "n")) {
                    message(paste(
                        "Invalid choice, please select one between y",
                        "or n"
                    ))
                } else {
                    break
                }
            }
            if (confirm == "y") {
                return(association_file)
            } else {
                next
            }
        }
    }
}

# Updates a files_found tibble with Files_count and Anomalies column.
#
# @param lups A files_found tibble obtained in a lookup function. Must contain
# the Files column (nested table quantification type and files)
# @keywords internal
#
# @return Updated files_found with Anomalies and Files_count columns
.trace_anomalies <- function(lups) {
    files_count <- purrr::pmap(lups, function(...) {
        temp <- tibble::tibble(...)
        an <- temp$Files
        an <- an %>%
            dplyr::group_by(.data$Quantification_type) %>%
            dplyr::summarise(
                Found = sum(!is.na(.data$Files_found)),
                .groups = "drop_last"
            )
        an
    })

    lups <- lups %>%
        dplyr::mutate(Files_count = files_count) %>%
        dplyr::mutate(
            Anomalies = purrr::map_lgl(
                .data$Files_count,
                ~ any(.x["Found"] != 1)
            ),
            .before = .data$Files
        )
    lups
}

# Looks up matrices to import given the association file and the
# root of the file system.
#
# @param association_file Tibble representing the association file
# @param quantification_type The type of quantification matrices to look for
# (one in `quantification_types()`)
# @param matrix_type The matrix_type to lookup (one between "annotated" or
# "not_annotated")
# @keywords internal
#' @importFrom tibble tibble
#' @importFrom fs dir_ls as_fs_path
#' @importFrom purrr pmap_dfr cross_df
#' @importFrom stringr str_replace_all
#' @importFrom tidyr nest
#
# @return A tibble containing all found files, including duplicates and missing
.lookup_matrices <- function(association_file,
    quantification_type,
    matrix_type,
    proj_col,
    pool_col) {
    path_col_names <- .path_cols_names()
    temp <- association_file %>%
        dplyr::select(
            .data[[proj_col]],
            .data[[pool_col]],
            .data[[path_col_names$quant]]
        ) %>%
        dplyr::distinct()

    ## Obtain a df with all possible combination of suffixes for each
    ## quantification
    suffixes <- matrix_file_suffixes() %>%
        dplyr::filter(
            .data$matrix_type == matrix_type,
            .data$quantification %in% quantification_type
        ) %>%
        dplyr::mutate(suffix_regex = paste0(stringr::str_replace_all(
            .data$file_suffix, "\\.", "\\\\."
        ), "$"))

    ## For each row in temp (aka for each ProjectID and
    ## concatenatePoolIDSeqRun) scan the quantification folder
    lups <- purrr::pmap_dfr(temp, function(...) {
        temp_row <- tibble::tibble(...)
        found <- purrr::pmap_dfr(suffixes, function(...) {
            ## For each quantification scan the folder for the
            ## corresponding suffixes
            cross_row <- tibble::tibble(...)
            matches <- fs::dir_ls(temp_row[[path_col_names$quant]],
                type = "file", fail = FALSE,
                regexp = cross_row$suffix_regex
            )
            if (length(matches) == 0) {
                matches <- NA_character_
            }
            tibble::tibble(
                Quantification_type = cross_row$quantification,
                Files_found = matches
            )
        })
        found <- found %>%
            dplyr::group_by(.data$Quantification_type) %>%
            dplyr::distinct() %>%
            dplyr::group_modify(~ {
                if (nrow(.x) > 1) {
                    if (any(is.na(.x$Files_found)) &
                        !all(is.na(.x$Files_found))) {
                        .x %>%
                            dplyr::filter(!is.na(.data$Files_found))
                    } else {
                        .x
                    }
                } else {
                    .x
                }
            })
        tibble::tibble(
            !!proj_col := temp_row[[proj_col]],
            !!pool_col := temp_row[[pool_col]],
            found
        )
    })
    lups <- lups %>%
        tidyr::nest(Files = c(.data$Quantification_type, .data$Files_found))
    lups <- .trace_anomalies(lups)
    lups
}


# Simple function to manage user input for duplicate file choices for each
# quantification type in a single project/pool pair (use internally in
# `.manage_anomalies_interactive`).
#
# @param q_types Vector of characters containing the unique quantification
# types that are detected as duplicates
# @param dupl The tibble containing quantification types and path to the files
# found for a single project/pool pair
# @keywords internal
#' @importFrom dplyr filter slice bind_rows
#' @importFrom purrr map reduce
#' @importFrom rlang .data
#
# @return An updated tibble containing files chosen for each type
.choose_duplicates_files_interactive <- function(q_types, dupl) {
    non_dup <- dupl %>% dplyr::filter(!.data$Quantification_type %in% q_types)
    type_choice <- purrr::map(q_types, function(x) {
        cat("---Quantification type: ", x, "---\n\n")
        temp <- dupl %>% dplyr::filter(.data$Quantification_type == x)
        cat("[0] QUIT\n",
            paste0("[", seq_along(temp$Files_found), "] ", temp$Files_found,
                collapse = "\n\n"
            ),
            sep = "\n"
        )
        repeat {
            cat("\n\nType the index number: ")
            connection <- getOption("ISAnalytics.connection", stdin())
            if (length(connection) > 1) {
                connection <- connection[7]
            }
            ch <- readLines(con = connection, n = 1)
            ch <- as.numeric(ch)
            if (!is.na(ch) &
                is.numeric(ch) &
                ch >= 0 & ch <= length(temp$Files_found)) {
                if (ch == 0) {
                    stop("Quitting", call. = FALSE)
                } else {
                    break
                }
            }
        }
        temp %>%
            dplyr::slice(ch) %>%
            dplyr::bind_rows(non_dup)
    })
    type_choice <- purrr::reduce(type_choice, dplyr::bind_rows)
}

# Manages anomalies for files found.
#
# The function manages anomalies found for files after scanning appropriate
# folders by:
# * Removing files not found (files for which Files_count$Found == 0 and Path
# is NA) and printing a message to notify user
# * Removing duplicates by asking the user which files to keep for each
# quantification type, pool and project
#
# @param files_found The tibble obtained via calling `.lookup_matrices`
# @keywords internal
#' @importFrom dplyr filter select rename bind_rows arrange
#' @importFrom tidyr unnest
#' @importFrom tibble as_tibble tibble
#' @importFrom purrr pmap flatten reduce
#
# @return A tibble containing for each project, pool and quantification type
# the files chosen (ideally 1 for each quantification type if found, no more
# than 1 per type)
.manage_anomalies_interactive <- function(files_found, proj_col, pool_col) {
    # Isolate anomalies in files found
    anomalies <- files_found %>% dplyr::filter(.data$Anomalies == TRUE)
    files_found <- files_found %>% dplyr::filter(.data$Anomalies == FALSE)
    # If there are no anomalies to fix, return chosen files
    if (nrow(anomalies) == 0) {
        files_to_import <- files_found %>%
            dplyr::select(
                .data[[proj_col]],
                .data[[pool_col]], .data$Files
            ) %>%
            tidyr::unnest(.data$Files) %>%
            dplyr::rename(Files_chosen = "Files_found")
        return(files_to_import)
    }
    repeat {
        # For each ProjectID and PoolID
        files_to_import <-
            purrr::pmap(anomalies, function(...) {
                l <- list(...)
                current <- tibble::as_tibble(l[c(
                    proj_col,
                    pool_col
                )])
                current_files <- tibble::as_tibble(purrr::flatten(l["Files"]))
                current_count <- tibble::as_tibble(
                    purrr::flatten(l["Files_count"])
                )
                # Find missing files first
                missing <- current_count %>% dplyr::filter(.data$Found == 0)
                if (nrow(missing) > 0) {
                    rlang::inform("Some files are missing and will be ignored")
                    current_files <- current_files %>%
                        dplyr::filter(!.data$Quantification_type %in%
                            missing$Quantification_type)
                }
                # Manage duplicates
                duplicate <- current_count %>% dplyr::filter(.data$Found > 1)
                if (nrow(duplicate) > 0) {
                    rlang::inform("Duplicates found for some files")
                    cat("Plese select one file for each group, type the index ",
                        "as requested\n\n",
                        sep = ""
                    )
                    cat("#### Project: ", current[[proj_col]], "####\n\n")
                    cat(
                        "#### Pool: ", current[[pool_col]],
                        "####\n\n"
                    )
                    current_files <- .choose_duplicates_files_interactive(
                        duplicate$Quantification_type, current_files
                    )
                }
                to_import <- tibble::tibble(
                    !!proj_col := current[[proj_col]],
                    !!pool_col := current[[pool_col]],
                    current_files
                )
                to_import <- to_import %>%
                    dplyr::rename(Files_chosen = "Files_found")
            })
        # Obtain single tibble for all anomalies
        files_to_import <- purrr::reduce(files_to_import, dplyr::bind_rows)
        cat("\nYour choices: ", sep = "\n")
        print(files_to_import)
        cat("\nFiles chosen details:\n")
        cat(paste0(
            seq_along(files_to_import$Files_chosen), "-",
            files_to_import$Files_chosen
        ), sep = "\n\n")
        cat("\nConfirm your choices? [y/n]", sep = "")
        repeat {
            connection <- getOption("ISAnalytics.connection", stdin())
            if (length(connection) > 1) {
                connection <- connection[8]
            }
            confirm <- readLines(con = connection, n = 1)
            if (!confirm %in% c("y", "n")) {
                message(paste(
                    "Invalid choice, please select one between",
                    "y or n"
                ))
            } else {
                break
            }
        }
        if (confirm == "y") {
            # Add non-anomalies
            if (nrow(files_found) > 0) {
                files_found <- files_found %>%
                    dplyr::select(
                        .data[[proj_col]],
                        .data[[pool_col]], .data$Files
                    ) %>%
                    tidyr::unnest(.data$Files) %>%
                    dplyr::rename(Files_chosen = "Files_found")
                files_to_import <- files_to_import %>%
                    dplyr::bind_rows(files_found) %>%
                    dplyr::arrange(.data[[proj_col]])
            }
            return(files_to_import)
        } else {
            next
        }
    }
}


# Internal function to match user defined patterns on a vector
# of file names.
#
# For each pattern specified by the user, the function tries to find a match
# on all the file names and combines the results as a tibble in which column
# names are the patterns and the values are TRUE if the pattern matched on the
# element at that index or FALSE otherwise.
#
# @param filenames A character vector of file names
# @param patterns A character vector of patterns to be matched
#
# @return A tibble
.pattern_matching <- function(filenames, patterns) {
    p_matches <- purrr::map(patterns, function(x) {
        mtc <- stringr::str_detect(filenames, x)
        tibble::as_tibble_col(mtc, column_name = x)
    })
    p_matches <- purrr::reduce(p_matches, function(x, y) {
        x %>% dplyr::bind_cols(y)
    })
    p_matches
}

# Helper function for checking if any of the elements of the list is true.
#
# @param ... A list of logical values
# @keywords internal
#
# @return TRUE if any of the parameters is true
.any_match <- function(...) {
    l <- unlist(list(...))
    any(l)
}

# Helper function for checking if all of the elements of the list is true.
#
# @param ... A list of logical values
# @keywords internal
#
# @return TRUE if all of the parameters is true
.all_match <- function(...) {
    l <- unlist(list(...))
    all(l)
}

# Updates files_found tibble according to pattern and matching options.
#
# @param files_nested The tibble containing Quantification_type and Files_found
# columns relative to a single project/pool pair
# @param p_matches The tibble representing the pattern matchings resulting from
# `pattern_matching`
# @param matching_opt The matching option
#' @importFrom purrr map reduce
#' @importFrom rlang .data
#' @importFrom fs as_fs_path
#
# @return An updated files_found tibble according to the matching option
.update_as_option <- function(files_nested, p_matches, matching_opt) {
    patterns <- colnames(p_matches)
    # Bind columns of Files nested tbl with pattern matches results
    p_matches <- files_nested %>% dplyr::bind_cols(p_matches)
    # Find the different quantification types in the files_found
    types <- p_matches %>%
        dplyr::select(.data$Quantification_type) %>%
        dplyr::distinct() %>%
        dplyr::pull(.data$Quantification_type)
    # For each quantification type
    select_matches <- function(x) {
        # Obtain a summary of matches (matches ALL patterns or ANY pattern)
        temp <- p_matches %>%
            dplyr::filter(.data$Quantification_type == x) %>%
            dplyr::rowwise() %>%
            dplyr::mutate(
                ALL = .all_match(
                    dplyr::c_across(dplyr::all_of(patterns))
                ),
                ANY = .any_match(dplyr::c_across(dplyr::all_of(patterns)))
            ) %>%
            dplyr::ungroup()
        # If the matching option is "ANY" - preserve files that match any of the
        # given patterns
        if (matching_opt == "ANY") {
            if (!all(is.na(temp$ANY)) & any(temp$ANY)) {
                # If there is at least 1 match in the group, files that match
                # are kept
                files_keep <- temp %>%
                    dplyr::filter(.data$ANY == TRUE) %>%
                    dplyr::select(.data$Quantification_type, .data$Files_found)
            } else {
                # If none of the files match, none is preserved, file is
                # converted to NA
                files_keep <- temp %>%
                    dplyr::slice(1) %>%
                    dplyr::mutate(Files_found = fs::as_fs_path(
                        NA_character_
                    )) %>%
                    dplyr::select(.data$Quantification_type, .data$Files_found)
            }
            files_keep
            # If the matching option is "ALL" - preserve only files that match
            # all the given patterns
        } else if (matching_opt == "ALL") {
            if (!all(is.na(temp$ALL)) & any(temp$ALL)) {
                # If there is at least 1 match in the group, files that match
                # are kept
                files_keep <- temp %>%
                    dplyr::filter(.data$ALL == TRUE) %>%
                    dplyr::select(.data$Quantification_type, .data$Files_found)
            } else {
                # If none of the files match none is preserved, file is
                # converted to NA
                files_keep <- temp %>%
                    dplyr::slice(1) %>%
                    dplyr::mutate(Files_found = fs::as_fs_path(
                        NA_character_
                    )) %>%
                    dplyr::select(.data$Quantification_type, .data$Files_found)
            }
            files_keep
            # If the matching option is "OPTIONAL" - preserve preferentially
            # files that match all or any of the given patterns, if none is
            # matching preserves file anyway
        } else {
            # Check first if some files match all the patterns
            if (!all(is.na(temp$ALL)) & any(temp$ALL)) {
                # Preserve only the ones that match
                files_keep <- temp %>%
                    dplyr::filter(.data$ALL == TRUE) %>%
                    dplyr::select(.data$Quantification_type, .data$Files_found)
            } else if (!all(is.na(temp$ANY)) & any(temp$ANY)) {
                # If none match all the patterns check files that match any of
                # the patterns and preserve only those
                files_keep <- temp %>%
                    dplyr::filter(.data$ANY == TRUE) %>%
                    dplyr::select(.data$Quantification_type, .data$Files_found)
            } else {
                # If there are no files that match any of the patterns simply
                # preserve the files found
                files_keep <- temp %>%
                    dplyr::select(.data$Quantification_type, .data$Files_found)
            }
            files_keep
        }
    }
    to_keep <- purrr::map(types, select_matches)
    to_keep <- purrr::reduce(to_keep, dplyr::bind_rows)
    to_keep
}

# Looks up matrices to import given the association file and the
# root of the file system for AUTO mode.
#
# @return A tibble containing all found files, including duplicates and missing
.lookup_matrices_auto <- function(association_file,
    quantification_type,
    matrix_type,
    patterns,
    matching_opt,
    proj_col,
    pool_col) {
    files_found <- .lookup_matrices(
        association_file,
        quantification_type,
        matrix_type,
        proj_col,
        pool_col
    )
    if (nrow(files_found) > 0 & !is.null(patterns)) {
        pre_filtering <- function(...) {
            l <- list(...)
            files_nested <- tibble::as_tibble(purrr::flatten(l["Files"]))
            split <- stringr::str_split(files_nested$Files_found, "\\/")
            filenames <- unlist(purrr::map(split, function(x) {
                utils::tail(x, n = 1)
            }))
            p_matches <- .pattern_matching(filenames, patterns)
            to_keep <- .update_as_option(files_nested, p_matches, matching_opt)
            to_keep
        }
        matching_patterns <- purrr::pmap(files_found, pre_filtering)
        files_found <- files_found %>%
            dplyr::mutate(Files = matching_patterns) %>%
            dplyr::select(
                .data[[proj_col]], .data[[pool_col]],
                .data$Files
            )
        files_found <- .trace_anomalies(files_found)
    }
    files_found
}

# Manages anomalies for files found.
#
# The function manages anomalies found for files after scanning appropriate
# folders by:
# * Removing files not found (files for which Files_count$Found == 0 and Path
# is NA) and printing a message to notify user
# * Removing duplicates
#
# @param files_found The tibble obtained via calling `.lookup_matrices_auto`
# @keywords internal
#' @importFrom dplyr filter select rename bind_rows arrange
#' @importFrom rlang .data
#' @importFrom tidyr unnest
#' @importFrom tibble as_tibble tibble
#' @importFrom purrr pmap flatten
#
# @return A tibble containing for each project, pool and quantification type
# the files chosen
.manage_anomalies_auto <- function(files_found, proj_col,
    pool_col) {
    # Isolate anomalies in files found
    anomalies <- files_found %>% dplyr::filter(.data$Anomalies == TRUE)
    files_found <- files_found %>% dplyr::filter(.data$Anomalies == FALSE)
    # If there are no anomalies to fix, return chosen files
    if (nrow(anomalies) == 0) {
        files_to_import <- files_found %>%
            dplyr::select(
                .data[[proj_col]],
                .data[[pool_col]], .data$Files
            ) %>%
            tidyr::unnest(.data$Files) %>%
            dplyr::rename(Files_chosen = "Files_found")
        return(files_to_import)
    }
    # For each ProjectID and PoolID
    process_anomalie <- function(...) {
        l <- list(...)
        current <- tibble::as_tibble(l[c(
            proj_col,
            pool_col
        )])
        current_files <- tibble::as_tibble(purrr::flatten(l["Files"]))
        current_count <- tibble::as_tibble(purrr::flatten(l["Files_count"]))
        reported_missing_or_dupl <- tibble::tibble(
            !!proj_col := character(),
            !!pool_col := character(),
            Quantification = character(),
            Missing = logical(),
            Duplicated = logical()
        )
        # Find missing files first
        missing <- current_count %>% dplyr::filter(.data$Found == 0)
        if (nrow(missing) > 0) {
            current_files <- current_files %>%
                dplyr::filter(!.data$Quantification_type %in%
                    missing$Quantification_type)
            reported_missing_or_dupl <- dplyr::bind_rows(
                reported_missing_or_dupl,
                missing %>%
                    dplyr::select(.data$Quantification_type) %>%
                    dplyr::mutate(
                        !!proj_col := current[[proj_col]],
                        !!pool_col := current[[pool_col]],
                        Missing = TRUE,
                        Duplicated = FALSE
                    ) %>%
                    dplyr::rename(Quantification = "Quantification_type")
            )
        }
        # Manage duplicates
        duplicate <- current_count %>% dplyr::filter(.data$Found > 1)
        if (nrow(duplicate) > 0) {
            current_files <- current_files %>%
                dplyr::filter(!.data$Quantification_type %in%
                    duplicate$Quantification_type)
            reported_missing_or_dupl <- dplyr::bind_rows(
                reported_missing_or_dupl,
                duplicate %>%
                    dplyr::select(.data$Quantification_type) %>%
                    dplyr::mutate(
                        !!proj_col := current[[proj_col]],
                        !!pool_col := current[[pool_col]],
                        Missing = FALSE,
                        Duplicated = TRUE
                    ) %>%
                    dplyr::rename(Quantification = "Quantification_type")
            )
        }
        to_import <- tibble::tibble(
            ProjectID = current[[proj_col]],
            concatenatePoolIDSeqRun =
                current[[pool_col]],
            current_files
        )
        to_import <- to_import %>% dplyr::rename(Files_chosen = "Files_found")
        return(list(to_import = to_import, reported = reported_missing_or_dupl))
    }
    files_to_import <- purrr::pmap(anomalies, process_anomalie)
    reported_anomalies <- purrr::map_df(files_to_import, ~ .x$reported)
    files_to_import <- purrr::map_df(files_to_import, ~ .x$to_import)
    # Report missing or duplicated
    if (getOption("ISAnalytics.verbose", TRUE) == TRUE &&
        nrow(reported_anomalies) > 0) {
        if (any(reported_anomalies$Missing)) {
            missing_msg <- c("Some files are missing and will be ignored",
                i = paste(
                    "Here is a summary of missing files:\n",
                    paste0(utils::capture.output(print(
                        reported_anomalies %>%
                            dplyr::filter(.data$Missing == TRUE) %>%
                            dplyr::select(
                                -.data$Duplicated,
                                -.data$Missing
                            )
                    )),
                    collapse = "\n"
                    )
                )
            )
            rlang::inform(missing_msg, class = "auto_mode_miss")
        }
        if (any(reported_anomalies$Duplicated)) {
            dupl_msg <- c(paste("Duplicates found for some files"),
                i = paste(
                    "In automatic mode duplicates are not preserved.",
                    "Use interactive mode for more accurate",
                    "file selection or provide appropriate file patterns",
                    "(see `?import_parallel_Vispa2Matrices`)"
                ),
                i = paste(
                    "Here is a summary of duplicated files:\n",
                    paste0(utils::capture.output(
                        print(
                            reported_anomalies %>%
                                dplyr::filter(.data$Duplicated == TRUE) %>%
                                dplyr::select(
                                    -.data$Duplicated,
                                    -.data$Missing
                                )
                        )
                    ), collapse = "\n")
                )
            )
            rlang::inform(dupl_msg, class = "auto_mode_dupl")
        }
    }
    # Add non-anomalies
    if (nrow(files_found) > 0) {
        files_found <- files_found %>%
            dplyr::select(
                .data[[proj_col]],
                .data[[pool_col]], .data$Files
            ) %>%
            tidyr::unnest(.data$Files) %>%
            dplyr::rename(Files_chosen = "Files_found")
        files_to_import <- files_to_import %>%
            dplyr::bind_rows(files_found) %>%
            dplyr::arrange(.data[[proj_col]], .data[[pool_col]])
    }
    files_to_import
}


# Internal function for parallel import of a single quantification
# type files.
#
# @param q_type The quantification type (single string)
# @param files Files_found table where absolute paths of chosen files
# are stored
# @param workers Number of parallel workers
# @keywords internal
#' @importFrom dplyr filter mutate bind_rows distinct
#' @importFrom BiocParallel SnowParam MulticoreParam bptry bplapply bpstop bpok
#' @importFrom purrr is_empty reduce
#
# @return A single tibble with all data from matrices of same quantification
# type in tidy format
.import_type <- function(q_type, files, cluster, import_matrix_args) {
    files <- files %>%
        dplyr::filter(.data$Quantification_type == q_type)
    sample_col_name <- if ("id_col_name" %in% names(import_matrix_args)) {
        import_matrix_args[["id_col_name"]]
    } else {
        pcr_id_column()
    }
    FUN <- function(x, arg_list) {
        do.call(.import_single_matrix, args = append(
            list(path = x),
            arg_list
        ))
    }
    # Import every file
    suppressMessages(suppressWarnings({
        matrices <- BiocParallel::bptry(
            BiocParallel::bplapply(files$Files_chosen,
                FUN = FUN,
                arg_list = import_matrix_args,
                BPPARAM = cluster
            )
        )
    }))
    correct <- BiocParallel::bpok(matrices)
    samples_count <- purrr::map2_int(matrices, correct, ~ {
        if (.y) {
            .x %>%
                dplyr::distinct(.data[[sample_col_name]]) %>%
                dplyr::pull(.data[[sample_col_name]]) %>%
                length()
        } else {
            NA_integer_
        }
    })
    distinct_is <- purrr::map2_int(matrices, correct, ~ {
        if (.y) {
            .x %>%
                dplyr::distinct(dplyr::across(
                    dplyr::all_of(mandatory_IS_vars())
                )) %>%
                nrow()
        } else {
            NA_integer_
        }
    })
    imported_files <- files %>%
        dplyr::mutate(
            Imported = correct,
            Number_of_samples = samples_count,
            Distinct_is = distinct_is
        )
    matrices <- matrices[correct]
    # Bind rows in single tibble for all files
    if (purrr::is_empty(matrices)) {
        return(list(matrix = NULL, imported_files = imported_files))
    }
    matrices <- purrr::reduce(matrices, function(x, y) {
        x %>%
            dplyr::bind_rows(y) %>%
            dplyr::distinct()
    })
    list(matrix = matrices, imported_files = imported_files)
}

# Internal function for importing all files for each quantification type.
#
# @param files_to_import The tibble containing the files to import
# @param workers Number of parallel workers
# @keywords internal
#' @importFrom dplyr select distinct bind_rows
#' @importFrom purrr map set_names reduce flatten
#' @importFrom tibble as_tibble
#
# @return A named list of tibbles
.parallel_import_merge <- function(files_to_import, workers,
    import_matrix_args) {
    # Find the actual quantification types included
    q_types <- files_to_import %>%
        dplyr::distinct(.data$Quantification_type) %>%
        dplyr::pull(.data$Quantification_type)
    # Register backend according to platform
    if (.Platform$OS.type == "windows") {
        p <- BiocParallel::SnowParam(
            workers = workers,
            stop.on.error = FALSE,
            progressbar = getOption("ISAnalytics.verbose", TRUE),
            tasks = nrow(files_to_import)
        )
    } else {
        p <- BiocParallel::MulticoreParam(
            workers = workers,
            stop.on.error = FALSE,
            progressbar = getOption("ISAnalytics.verbose", TRUE),
            tasks = nrow(files_to_import)
        )
    }
    # Import and merge for every quantification type
    imported_matrices <- purrr::map(q_types,
        .f = ~ .import_type(
            .x,
            files_to_import,
            p,
            import_matrix_args
        )
    ) %>%
        purrr::set_names(q_types)
    BiocParallel::bpstop(p)
    summary_files <- purrr::map(imported_matrices, ~ .x$imported_files) %>%
        purrr::reduce(dplyr::bind_rows)
    imported_matrices <- purrr::map(imported_matrices, ~ .x$matrix)

    list(matrix = imported_matrices, summary = summary_files)
}

#### ---- Internals for collision removal ----####

#---- USED IN : remove_collisions ----

# Used as a wrapper for basic correctness check on collisions input
# - Returns the mode (LIST or MULTI)
.collisions_check_input_matrices <- function(x, quant_cols) {
    stopifnot(is.list(x) & !is.null(names(x)))
    stopifnot(is.character(quant_cols) && all(!is.na(names(quant_cols))))
    if (!all(names(quant_cols) %in% quantification_types())) {
        rlang::abort(
            .quantifications_names_err(
                quant_cols[!names(quant_cols) %in% quantification_types()]
            )
        )
    }
    mode <- NULL
    if (!is.data.frame(x)) {
        #---- If input x is provided as list of data frames
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
            cAmp_col <- .check_sample_col(m)
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
            names(message) <- NULL
            rlang::abort(c("Matrices miss required info, aborting", message),
                class = "coll_matrix_issues"
            )
        }
        ## Transform the list in a multi-quant matrix
        mode <- "LIST"
        quant_cols_lst <- as.list(quant_cols)
        args <- append(list(x = x), quant_cols_lst)
        x <- rlang::exec(comparison_matrix, !!!args)
        rlang::env_poke(env = rlang::caller_env(), nm = "x", value = x)
    } else {
        #---- If input x is provided as data frame
        if (.check_mandatory_vars(x) == FALSE) {
            rlang::abort(.missing_mand_vars())
        }
        if (.check_sample_col(x) == FALSE) {
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
    return(mode)
}

# Used as a wrapper for basic correctness check on collisions input -
# for association file. Checks if all columns provided in input are present
# and in correct format and checks if required tags are present
# - Returns the selected and found tags if no errors are raised
.collisions_check_input_af <- function(association_file,
    date_col,
    independent_sample_id) {
    stopifnot(is.data.frame(association_file))
    stopifnot(is.character(independent_sample_id) &&
        !purrr::is_empty(independent_sample_id))
    af_specs <- association_file_columns(TRUE)
    ## Check if input columns are present actual af
    user_input_cols <- unique(c(independent_sample_id, date_col))
    if (!all(user_input_cols %in%
        colnames(association_file))) {
        missing_cols <- user_input_cols[!user_input_cols %in%
            colnames(association_file)]
        rlang::abort(.missing_af_needed_cols(missing_cols))
    }
    ## Check tags
    required_tags <- list(
        project_id = "char",
        pool_id = "char",
        pcr_replicate = c("int", "numeric"),
        pcr_repl_id = "char"
    )
    req_tag_cols <- .check_required_cols(required_tags, af_specs, "error")
    if (!all(req_tag_cols$names %in% colnames(association_file))) {
        rlang::abort(.missing_af_needed_cols(req_tag_cols$names[
            !req_tag_cols$names %in% colnames(association_file)
        ]))
    }
    ## Check date_col
    if (!lubridate::is.Date(association_file[[date_col]])) {
        not_date_err <- c(paste0("'", date_col, "'", " is not a date vector"))
        rlang::abort(not_date_err, class = "not_date_coll_err")
    }
    if (any(is.na(association_file[[date_col]]))) {
        rlang::abort(.na_in_date_col())
    }
    return(req_tag_cols)
}


# Checks if association file contains more information than the matrix.
#
# Used to notify the user that wants to know if for the projects contained in
# the examined matrix there are additional CompleteAmplificationIDs contained
# in the association file that weren't included in the integration matrix (for
# example failed runs).
#' @importFrom rlang .data sym
.check_same_info <- function(association_file, df, req_tag_cols,
    indep_sample_id) {
    data.table::setDT(req_tag_cols)
    joined <- association_file %>%
        dplyr::left_join(df, by = pcr_id_column())
    proj_col_name <- req_tag_cols[eval(rlang::expr(
        !!sym("tag") == "project_id"
    ))][["names"]]
    projects <- joined %>%
        dplyr::filter(!dplyr::if_all(
            .cols = dplyr::all_of(mandatory_IS_vars()),
            .fns = is.na
        )) %>%
        dplyr::distinct(.data[[proj_col_name]]) %>%
        dplyr::pull(.data[[proj_col_name]])
    wanted_tags <- c(
        "subject", "tissue", "cell_marker", "tp_days"
    )
    wanted <- association_file_columns(TRUE)
    data.table::setDT(wanted)
    wanted <- wanted[eval(sym("tag")) %in% wanted_tags &
        eval(sym("names")) %in% colnames(association_file), ]
    wanted <- data.table::rbindlist(list(req_tag_cols, wanted))
    missing_in_df <- joined %>%
        dplyr::filter(
            dplyr::if_all(
                .cols = dplyr::all_of(mandatory_IS_vars()),
                .fns = is.na
            ) & .data[[proj_col_name]] %in% projects
        ) %>%
        dplyr::distinct(dplyr::across(
            dplyr::all_of(c(wanted$names, indep_sample_id))
        ))
    reduced_af <- association_file %>%
        dplyr::filter(
            .data[[proj_col_name]] %in% projects
        )
    data.table::setDT(reduced_af)
    data.table::setDT(missing_in_df)
    list(miss = missing_in_df, reduced_af = reduced_af)
}


# Identifies independent samples and separates the joined_df in
# collisions and non-collisions
#' @importFrom data.table .SD .N
#' @importFrom rlang sym
# @return A named list containing the splitted joined_df for collisions and
# non-collisions
.identify_independent_samples <- function(joined, indep_sample_id) {
    indep_syms <- purrr::map(indep_sample_id, rlang::sym)
    temp <- joined %>%
        dplyr::select(dplyr::all_of(
            c(mandatory_IS_vars(), indep_sample_id)
        )) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(mandatory_IS_vars()))) %>%
        dplyr::summarise(
            N = dplyr::n_distinct(!!!indep_syms),
            .groups = "drop"
        ) %>%
        dplyr::filter(.data$N > 1)
    non_collisions <- joined %>%
        dplyr::anti_join(temp, by = mandatory_IS_vars())
    collisions <- joined %>%
        dplyr::semi_join(temp, by = mandatory_IS_vars())
    data.table::setDT(non_collisions)
    data.table::setDT(collisions)
    list(collisions = collisions, non_collisions = non_collisions)
}

# Internal for date discrimination in collision removal.
#
# It's the first of 4 steps in the algorithm for collision removal: it tries to
# find a single sample who has an associated date which is earlier than any
# other. If comparison is not possible the analysis fails and returns the
# original input data.
#
# @param x - All metadata associated with a single collision event in
#            a data.table format (does not contain mandatory IS vars)
# @param date_col - The name of the date column chosen for collision removal
# @param ind_sample_key - character vector containing the columns that identify
#                       independent samples
# @return A named list with:
# * $data: a data.table, containing the data (unmodified or modified)
# * $check: a logical value indicating whether the analysis was successful or
# not (and therefore there is the need to perform the next step)
.discriminate_by_date <- function(x, date_col, ind_sample_key) {
    # Test for all dates equality
    all_equal_dates <- length(unique(x[[date_col]])) == 1
    if (all_equal_dates == TRUE) {
        return(list(data = x, check = FALSE))
    }
    # If not all are equal sort them asc
    x <- dplyr::arrange(x, .data[[date_col]])
    # Test if first date is unique
    dates <- x[[date_col]]
    if (length(dates[dates %in% dates[1]]) > 1) {
        return(list(data = x, check = FALSE))
    }
    winning_sample <- x[1, mget(ind_sample_key)]
    # Filter the winning rows
    x <- x %>%
        dplyr::semi_join(winning_sample, by = ind_sample_key)
    data.table::setDT(x)
    return(list(data = x, check = TRUE))
}

# Internal for replicate discrimination in collision removal.
#
# It's the second of 4 steps in the algorithm for collision removal:
# grouping by independent sample it counts the number of
# rows (replicates) found for each group, orders them from biggest to smallest
# and, if a single group has more rows than any other group the integration is
# assigned to that sample, otherwise the analysis fails and returns the
# original input to be submitted to the next step.
#
# @param x - All metadata associated with a single collision event in
#            a data.table format (does not contain mandatory IS vars)
# @param repl_col - The name of the column containing the replicate number
# @param ind_sample_key - character vector containing the columns that identify
#                       independent samples
#' @importFrom data.table .N
# @return A named list with:
# * $data: a data.table, containing the data (unmodified or modified)
# * $check: a logical value indicating whether the analysis was successful or
# not (and therefore there is the need to perform the next step)
.discriminate_by_replicate <- function(x, repl_col, ind_sample_key) {
    temp <- x %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(ind_sample_key))) %>%
        dplyr::summarise(N = dplyr::n(), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(.data$N))
    data.table::setDT(temp)
    if (length(temp$N) != 1 & !temp$N[1] > temp$N[2]) {
        return(list(data = x, check = FALSE))
    }
    temp <- temp[1, mget(ind_sample_key)]
    x <- x %>%
        dplyr::semi_join(temp, by = ind_sample_key)
    data.table::setDT(x)
    return(list(data = x, check = TRUE))
}

# Internal for sequence count discrimination in collision removal.
#
# It's the third of 4 steps in the algorithm for collision removal:
# grouping by independent sample, it sums the value of
# the sequence count for each group, sorts the groups from highest to lowest
# cumulative value and then checks the ratio between the first element and the
# second: if the ratio is > `reads_ratio` the integration is assigned to
# the first group, otherwise the analysis fails and returns the original input.
#
# @param x - All metadata associated with a single collision event in
#            a data.table format (does not contain mandatory IS vars)
# @param reads_ratio The value of the ratio between sequence count values to
# check
# @param seqCount_col The name of the sequence count column (support
# for multi quantification matrix)
# @param ind_sample_key - character vector containing the columns that identify
#                       independent samples
# @return A named list with:
# * $data: a tibble, containing the data (unmodified or modified)
# * $check: a logical value indicating whether the analysis was successful or
# not (and therefore there is the need to perform the next step)
.discriminate_by_seqCount <- function(x,
    reads_ratio,
    seqCount_col,
    ind_sample_key) {
    temp <- x %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(ind_sample_key))) %>%
        dplyr::summarise(sum = sum(.data[[seqCount_col]]), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(.data$sum))
    data.table::setDT(temp)
    ratio <- temp$sum[1] / temp$sum[2]
    if (!ratio > reads_ratio) {
        return(list(data = x, check = FALSE))
    }
    temp <- temp[1, mget(ind_sample_key)]
    x <- x %>%
        dplyr::semi_join(temp, by = ind_sample_key)
    data.table::setDT(x)
    return(list(data = x, check = TRUE))
}

# Internal function that performs four-step-check of collisions for
# a single integration.
#
# @param x - Represents the sub-table (data.table) containing coordinates
#            and metadata about a SINGLE collision event (same values for
#            mandatory IS vars)
# @param date_col The date column to consider
# @param repl_col - The name of the column containing the replicate number
# @param reads_ratio The value of the ratio between sequence count values to
# check
# @param seqCount_col The name of the sequence count column (support
# for multi quantification matrix)
# @param ind_sample_key - character vector containing the columns that identify
#                         independent samples
# @return A list with:
# * $data: an updated tibble with processed collisions or NULL if no
# criteria was sufficient
# * $reassigned: 1 if the integration was successfully reassigned, 0 otherwise
# * $removed: 1 if the integration is removed entirely because no criteria was
# met, 0 otherwise
.four_step_check <- function(x,
    date_col,
    repl_col,
    reads_ratio,
    seqCount_col,
    ind_sample_key) {
    current_data <- x[, !mandatory_IS_vars()]
    # Try to discriminate by date
    result <- .discriminate_by_date(current_data, date_col, ind_sample_key)
    current_data <- result$data
    if (result$check == TRUE) {
        coordinates <- x[seq_len(nrow(current_data)), mget(mandatory_IS_vars())]
        res <- cbind(coordinates, current_data)
        return(list(data = res, reassigned = 1, removed = 0))
    }
    # If first check fails try to discriminate by replicate
    result <- .discriminate_by_replicate(current_data, repl_col, ind_sample_key)
    current_data <- result$data
    if (result$check == TRUE) {
        coordinates <- x[seq_len(nrow(current_data)), mget(mandatory_IS_vars())]
        res <- cbind(coordinates, current_data)
        return(list(data = res, reassigned = 1, removed = 0))
    }
    # If second check fails try to discriminate by seqCount
    result <- .discriminate_by_seqCount(
        current_data, reads_ratio, seqCount_col,
        ind_sample_key
    )
    current_data <- result$data
    if (result$check == TRUE) {
        coordinates <- x[seq_len(nrow(current_data)), mget(mandatory_IS_vars())]
        res <- cbind(coordinates, current_data)
        return(list(data = res, reassigned = 1, removed = 0))
    }
    # If all check fails remove the integration from all subjects
    return(list(data = NULL, reassigned = 0, removed = 1))
}


# Internal function to process collisions on multiple integrations,
# parallelized.
# @return A list containing the updated collisions, a numeric value
# representing the number of integrations removed and a numeric value
# representing the number of integrations reassigned
.process_collisions <- function(collisions,
    date_col,
    repl_col,
    reads_ratio,
    seqCount_col,
    ind_sample_key,
    max_workers) {
    # Split by IS
    split_data <- split(collisions, by = mandatory_IS_vars())
    # Manage workers
    if (.Platform$OS.type == "windows" & is.null(max_workers)) {
        max_workers <- BiocParallel::snowWorkers()
    } else if (.Platform$OS.type != "windows" & is.null(max_workers)) {
        max_workers <- BiocParallel::multicoreWorkers()
    }
    # Register backend according to platform
    if (.Platform$OS.type == "windows") {
        p <- BiocParallel::SnowParam(
            stop.on.error = TRUE,
            progressbar = getOption("ISAnalytics.verbose", TRUE),
            tasks = length(split_data),
            workers = max_workers
        )
    } else {
        p <- BiocParallel::MulticoreParam(
            stop.on.error = TRUE,
            progressbar = getOption("ISAnalytics.verbose", TRUE),
            tasks = length(split_data),
            workers = max_workers
        )
    }
    # For each chunk process collisions in parallel
    result <- BiocParallel::bplapply(split_data,
        FUN = .four_step_check,
        date_col = date_col,
        repl_col = repl_col,
        reads_ratio = reads_ratio,
        seqCount_col = seqCount_col,
        ind_sample_key = ind_sample_key,
        BPPARAM = p
    )
    BiocParallel::bpstop(p)
    # For each element of result extract and bind the 3 components
    processed_collisions_df <- purrr::map(result, ~ .x$data) %>%
        purrr::reduce(~ data.table::rbindlist(list(.x, .y)))
    removed_total <- purrr::map(result, ~ .x$removed) %>%
        purrr::reduce(sum)
    reassigned_total <- purrr::map(result, ~ .x$reassigned) %>%
        purrr::reduce(sum)
    # Return the summary
    list(
        coll = processed_collisions_df, removed = removed_total,
        reassigned = reassigned_total
    )
}

# Internal function to obtain a summary table after collision processing.
#
# @param before The joined table before collision removing
# (obtained via `.join_matrix_af`)
# @param after The final matrix obtained after collision processing
# @param association_file The association file
#' @importFrom rlang .data
# @keywords internal
#
# @return A tibble with a summary containing for each SubjectID the number of
# integrations found before and after, the sum of the value of the sequence
# count for each subject before and after and the corresponding deltas.
.summary_table <- function(before, after, seqCount_col) {
    after_matr <- after %>%
        dplyr::select(dplyr::all_of(colnames(after)), .data$SubjectID)

    n_int_after <- after_matr %>%
        dplyr::group_by(.data$SubjectID) %>%
        dplyr::tally(name = "count_integrations_after")
    sumReads_after <- after_matr %>%
        dplyr::group_by(.data$SubjectID) %>%
        dplyr::summarise(
            sumSeqReads_after = sum(.data[[seqCount_col]]),
            .groups = "drop_last"
        )
    temp_aft <- n_int_after %>% dplyr::left_join(sumReads_after,
        by = c("SubjectID")
    )

    before_matr <- before %>%
        dplyr::select(dplyr::all_of(colnames(after)), .data$SubjectID)
    n_int_before <- before_matr %>%
        dplyr::group_by(.data$SubjectID) %>%
        dplyr::tally(name = "count_integrations_before")
    sumReads_before <- before_matr %>%
        dplyr::group_by(.data$SubjectID) %>%
        dplyr::summarise(
            sumSeqReads_before = sum(.data[[seqCount_col]]),
            .groups = "drop_last"
        )
    temp_bef <- n_int_before %>% dplyr::left_join(sumReads_before,
        by = c("SubjectID")
    )
    summary <- temp_bef %>%
        dplyr::left_join(temp_aft, by = c("SubjectID")) %>%
        dplyr::mutate(
            delta_integrations =
                .data$count_integrations_before -
                    .data$count_integrations_after,
            delta_seqReads = .data$sumSeqReads_before - .data$sumSeqReads_after
        )
    summary
}

# Internal for obtaining summary info about the input sequence count matrix
# as is (no pre-processing).
# Info contain
# - Number of distinct integration sites
# - Total for each quantification present (if multi-quant) or seqCount
.summary_input <- function(input_df, quant_cols) {
    names(quant_cols) <- NULL
    # Distinct integration sites
    n_IS <- input_df %>%
        dplyr::distinct(dplyr::across(mandatory_IS_vars())) %>%
        nrow()
    quant_totals <- input_df %>%
        dplyr::select(dplyr::all_of(quant_cols)) %>%
        dplyr::summarise(dplyr::across(
            .cols = dplyr::all_of(quant_cols),
            .fns = list(
                sum = ~ sum(.x, na.rm = TRUE)
            ),
            .names = "{.col}"
        ),
        .groups = "drop"
        ) %>%
        as.list()
    list(total_iss = n_IS, quant_totals = quant_totals)
}

.per_pool_stats <- function(joined, quant_cols, pool_col) {
    ## Joined is the matrix already joined with metadata
    df_by_pool <- joined %>%
        dplyr::group_by(.data[[pool_col]]) %>%
        dplyr::summarise(dplyr::across(
            .cols = dplyr::all_of(quant_cols),
            .fns = list(
                sum = ~ sum(.x, na.rm = TRUE),
                count = length,
                describe = ~ psych::describe(.x)
            )
        ))
    df_cols <- purrr::map_lgl(df_by_pool, ~ is.data.frame(.x))
    if (any(df_cols == TRUE)) {
        df_cols <- df_cols[df_cols]
        df_cols <- names(df_cols)
        dfss <- purrr::map(df_cols, function(col) {
            exp <- rlang::expr(`$`(df_by_pool, !!col))
            df <- rlang::eval_tidy(exp)
            df <- df %>% dplyr::rename_with(.fn = ~ paste0(col, "_", .x))
            df
        }) %>% purrr::set_names(df_cols)
        for (dfc in df_cols) {
            df_by_pool <- df_by_pool %>%
                dplyr::select(-dplyr::all_of(dfc)) %>%
                dplyr::bind_cols(dfss[[dfc]])
        }
    }
    df_by_pool
}

.collisions_obtain_report_summaries <- function(x,
    association_file,
    quant_cols, missing_ind,
    pcr_col, pre_process,
    collisions,
    removed, reassigned,
    joined, post_joined, pool_col,
    final_matr,
    replicate_n_col,
    independent_sample_id,
    seq_count_col) {
    input_summary <- .summary_input(x, quant_cols)
    missing_smpl <- if (!is.null(missing_ind)) {
        x[missing_ind, ] %>%
            dplyr::group_by(.data[[pcr_col]]) %>%
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
        MATRIX = unique(x[[pcr_col]]),
        AF = unique(association_file[[pcr_col]])
    )
    pre_summary <- .summary_input(pre_process, quant_cols)
    per_pool_stats_pre <- .per_pool_stats(joined, quant_cols, pool_col)
    coll_info <- list(
        coll_n = collisions %>%
            dplyr::distinct(
                dplyr::across(dplyr::all_of(mandatory_IS_vars()))
            ) %>%
            nrow(),
        removed = removed,
        reassigned = reassigned
    )
    post_info <- .summary_input(final_matr, quant_cols)
    post_per_pool_stats <- .per_pool_stats(
        post_joined,
        quant_cols,
        pool_col
    )
    summary_tbl <- .summary_table(
        before = joined, after = post_joined,
        seqCount_col = seq_count_col
    )
    return(
        list(
            input_info = input_summary,
            missing_info = missing_smpl,
            samples_info = samples_info,
            pre_info = pre_summary,
            pre_stats = per_pool_stats_pre,
            coll_info = coll_info,
            post_info = post_info,
            post_stats = post_per_pool_stats,
            summary_post = summary_tbl
        )
    )
}

.collisions_obtain_sharing_heatmaps <- function(joined,
    independent_sample_id,
    post_joined,
    report_path) {
    sharing_heatmaps_pre <- sharing_heatmaps_post <- NULL
    withCallingHandlers(
        {
            withRestarts(
                {
                    sharing_pre <- is_sharing(joined,
                        group_key = independent_sample_id,
                        n_comp = 2,
                        is_count = FALSE,
                        minimal = FALSE,
                        include_self_comp = TRUE
                    )
                    sharing_post <- is_sharing(post_joined,
                        group_key = independent_sample_id,
                        n_comp = 2,
                        is_count = FALSE,
                        minimal = FALSE,
                        include_self_comp = TRUE
                    )
                    too_many_samples <- nrow(dplyr::distinct(
                        joined,
                        dplyr::across(
                            dplyr::all_of(independent_sample_id)
                        )
                    )) > 25
                    sharing_heatmaps_pre <- if (too_many_samples) {
                        sharing_heatmap(sharing_pre, interactive = FALSE)
                    } else {
                        sharing_heatmap(sharing_pre, interactive = TRUE)
                    }
                    sharing_heatmaps_post <- if (too_many_samples) {
                        sharing_heatmap(sharing_post, interactive = FALSE)
                    } else {
                        sharing_heatmap(sharing_post, interactive = TRUE)
                    }
                    if (too_many_samples) {
                        skip_sharing_msg <- c(
                            "Too many independent samples",
                            i = paste(
                                "To avoid issues with report size,",
                                "sharing heatmaps will be saved",
                                "as separate files in the same",
                                "folder"
                            )
                        )
                        rlang::inform(skip_sharing_msg)
                        fs::dir_create(report_path)
                        prefix_filename <- paste(lubridate::today(),
                            "ISAnalytics_collision_removal",
                            sep = "_"
                        )
                        for (map_type in names(sharing_heatmaps_pre)) {
                            filename <- paste0(paste(prefix_filename,
                                "pre-processing-sharing",
                                map_type,
                                sep = "_"
                            ), ".pdf")
                            ggplot2::ggsave(
                                plot = sharing_heatmaps_pre[[map_type]],
                                filename = filename,
                                path = report_path,
                                width = 15, height = 15
                            )
                        }
                        for (map_type in names(sharing_heatmaps_post)) {
                            filename <- paste0(
                                paste(prefix_filename,
                                    "post-processing-sharing",
                                    map_type,
                                    sep = "_"
                                ),
                                ".pdf"
                            )
                            ggplot2::ggsave(
                                plot = sharing_heatmaps_post[[map_type]],
                                filename = filename,
                                path = report_path,
                                width = 15, height = 15
                            )
                        }
                        sharing_heatmaps_pre <- NULL
                        sharing_heatmaps_post <- NULL
                    }
                },
                sharing_err = function(e) {
                    rlang::inform(c(paste(
                        "Unable to compute sharing:",
                        conditionMessage(e)
                    ),
                    i = "Skipping"
                    ))
                }
            )
        },
        error = function(cnd) {
            rest <- findRestart("sharing_err")
            invokeRestart(rest, cnd)
        }
    )
    return(
        list(
            sharing_pre = sharing_heatmaps_pre,
            sharing_post = sharing_heatmaps_post
        )
    )
}

#### ---- Internals for aggregate functions ----####

#---- USED IN : aggregate_metadata ----

# Aggregates the association file based on the function table.
#' @importFrom rlang .data is_function is_formula
#' @importFrom tibble tibble
#' @importFrom purrr pmap_df
# @keywords internal
#
# @return A tibble
.aggregate_meta <- function(association_file, grouping_keys, function_tbl) {
    ## Discard columns that are not present in the af
    function_tbl <- function_tbl %>%
        dplyr::filter(.data$Column %in% colnames(association_file))
    ## If no columns are left return
    if (nrow(function_tbl) == 0) {
        return(NULL)
    }
    function_tbl <- function_tbl %>%
        tidyr::nest(cols = .data$Column)
    apply_function <- function(Function, Args, Output_colname, cols) {
        Function <- list(Function)
        Args <- list(Args)
        row <- tibble::tibble(Function, Args,
            Output_colname,
            cols = list(cols)
        )
        if (rlang::is_formula(row$Function[[1]])) {
            res <- association_file %>%
                dplyr::group_by(dplyr::across(dplyr::all_of(grouping_keys))) %>%
                dplyr::summarise(
                    dplyr::across(
                        dplyr::all_of(row$cols[[1]]$Column),
                        .fns = row$Function[[1]],
                        .names = row$Output_colname
                    ),
                    .groups = "drop"
                )
            return(res)
        }
        if (rlang::is_function(row$Function[[1]])) {
            args <- if (is.na(row$Args[[1]]) | is.null(is.na(row$Args[[1]]))) {
                list()
            } else {
                row$Args[[1]]
            }
            res <- association_file %>%
                dplyr::group_by(dplyr::across(dplyr::all_of(grouping_keys))) %>%
                dplyr::summarise(
                    dplyr::across(
                        .cols = dplyr::all_of(row$cols[[1]]$Column),
                        .fns = row$Function[[1]],
                        .names = row$Output_colname,
                        !!!args
                    ),
                    .groups = "drop"
                )
            return(res)
        }
        return(NULL)
    }
    agg <- purrr::pmap(function_tbl, apply_function)
    agg <- purrr::reduce(agg, ~ dplyr::full_join(.x, .y, by = grouping_keys))
    return(agg)
}

#---- USED IN : aggregate_values_by_key ----

# Internal function that performs aggregation on values with a lambda
# operation.
# - x is a single matrix
.aggregate_lambda <- function(x, af, key, value_cols, lambda, group,
    join_af_by) {
    cols_to_check <- c(group, key, value_cols, join_af_by)
    joint <- x %>%
        dplyr::left_join(af, by = dplyr::all_of(join_af_by))
    if (any(!cols_to_check %in% colnames(joint))) {
        missing <- cols_to_check[!cols_to_check %in% colnames(joint)]
        rlang::abort(.missing_user_cols_error(missing))
    }
    cols <- c(colnames(x), key)
    aggregated_matrices <- joint %>%
        dplyr::select(dplyr::all_of(cols)) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c(group, key)))) %>%
        dplyr::summarise(dplyr::across(
            .cols = dplyr::all_of(value_cols),
            .fns = lambda
        ), .groups = "drop")
    aggregated_matrices
}


#### ---- Internals for re-calibration functions ----####

#---- USED IN : compute_near_integrations ----

# Finds a unique maximum value in a vector of numbers.
#
# Unlike the standard `max` function, this function returns a single maximum
# value if and only if a maximum exists, more in details that means:
# * It ignores NA values by default
# * If the identified maximum value is not unique in the vector, it returns
# an empty vector instead
#
# @param x A numeric or integer vector
#' @importFrom purrr is_empty
#
# @return A single numeric value or an empty numeric vector
# @keywords internal
.find_unique_max <- function(x) {
    if (any(is.na(x))) {
        x <- x[!is.na(x)]
    }
    uniques <- unique(x)
    if (purrr::is_empty(uniques)) {
        return(uniques)
    }
    if (length(uniques) == 1) {
        if (length(x[x == uniques]) > 1) {
            return(numeric(0))
        } else {
            return(uniques)
        }
    }
    max_val <- max(uniques)
    if (length(x[x == max_val]) > 1) {
        return(numeric(0))
    }
    max_val
}

# Internal function that implements the sliding window algorithm.
#
# **NOTE: this function is meant to be called on a SINGLE GROUP,
# meaning a subset of an integration matrix in which all rows
# share the same chr (and optionally same strand).**
#' @importFrom dplyr arrange
#' @importFrom purrr is_empty map_dfr
#' @importFrom rlang expr eval_tidy .data sym
#' @importFrom data.table setDT setnames `%chin%`
# @keywords internal
#
# @return A named list with recalibrated matrix and recalibration map.
.sliding_window <- function(x,
    threshold,
    keep_criteria,
    annotated,
    num_cols,
    max_val_col,
    sample_col,
    req_tags,
    add_col_lambdas,
    produce_map) {
    locus_col <- req_tags %>%
        dplyr::filter(.data$tag == "locus") %>%
        dplyr::pull(.data$names)
    ## Order by integration_locus
    x <- x %>% dplyr::arrange(.data[[locus_col]])
    x <- data.table::setDT(x)
    map_recalibr <- if (produce_map == TRUE) {
        req_vars <- req_tags$names
        rec_map <- unique(x[, mget(req_vars)])
        data.table::setnames(
            x = rec_map,
            old = req_vars,
            new = paste0(req_vars, "_before")
        )
        na_type <- purrr::map(req_vars, ~ {
            col_type <- typeof(x[, get(.x)])
            fun_name <- paste0("as.", col_type)
            return(do.call(what = fun_name, args = list(NA)))
        })
        rec_map[, c(paste(req_vars, "after", sep = "_")) := na_type]
    } else {
        NULL
    }
    index <- 1
    repeat {
        # If index is the last row in the data frame return
        if (index == nrow(x)) {
            if (produce_map == TRUE) {
                is_predicate <- purrr::map(req_vars, ~ {
                    var_before <- sym(paste0(.x, "_before"))
                    value <- x[index, get(.x)]
                    rlang::expr(!!var_before == !!value)
                }) %>%
                    purrr::reduce(~ rlang::expr(!!.x & !!.y))
                map_recalibr[
                    eval(is_predicate),
                    c(
                        paste0(req_vars, "_after")
                    ) :=
                        purrr::map(
                            req_vars,
                            ~ eval(sym(paste0(.x, "_before")))
                        )
                ]
            }
            return(list(recalibrated_matrix = x, map = map_recalibr))
        }
        ## Compute interval for row
        interval <- x[index, ][[locus_col]] + 1 + threshold
        ## Look ahead for every integration that falls in the interval
        near <- numeric()
        k <- index
        repeat {
            if (k == nrow(x)) {
                break
            }
            k <- k + 1
            if (x[k, ][[locus_col]] < interval) {
                # Saves the indexes of the rows that are in the interval
                near <- append(near, k)
            } else {
                break
            }
        }
        window <- c(index, near)
        row_to_keep <- index
        if (!purrr::is_empty(near)) {
            ## Change loci according to criteria
            ######## CRITERIA PROCESSING
            if (keep_criteria == "max_value") {
                expr <- rlang::expr(`$`(x[window, ], !!max_val_col))
                to_check <- rlang::eval_tidy(expr)
                max <- .find_unique_max(to_check)
                if (!purrr::is_empty(max)) {
                    row_to_keep <- window[which(to_check == max)]
                    near <- window[!window == row_to_keep]
                }
            }
            # Fill map if needed
            if (produce_map == TRUE) {
                is_predicate <- purrr::map(req_vars, ~ {
                    var_before <- sym(paste0(.x, "_before"))
                    value <- x[window, get(.x)]
                    if (!is.character(value)) {
                        return(rlang::expr(!!var_before %in% !!value))
                    }
                    rlang::expr(!!var_before %chin% !!value)
                }) %>%
                    purrr::reduce(~ rlang::expr(!!.x & !!.y))
                map_recalibr[
                    eval(is_predicate),
                    c(
                        paste0(req_vars, "_after")
                    ) :=
                        purrr::map(
                            req_vars,
                            ~ x[row_to_keep][[.x]]
                        )
                ]
            }
            # Change loci and other ids of near integrations
            fields_to_change <- req_tags$names
            if (annotated) {
                fields_to_change <- c(fields_to_change, annotation_IS_vars())
                x[near, c(
                    req_tags$names,
                    annotation_IS_vars()
                ) := x[
                    row_to_keep,
                    mget(fields_to_change)
                ]]
            } else {
                x[near, c(
                    req_tags$names
                ) := x[
                    row_to_keep,
                    mget(fields_to_change)
                ]]
            }
            ## Aggregate same IDs
            starting_rows <- nrow(x)
            repeat {
                t <- x[[sample_col]][window]
                d <- unique(t[duplicated(t)])
                if (purrr::is_empty(d)) {
                    break
                }
                dupl_indexes <- which(t == d[1])
                values_sum <- colSums(x[window[dupl_indexes], mget(num_cols)],
                    na.rm = TRUE
                )
                x[window[dupl_indexes[1]], c(num_cols) := as.list(values_sum)]
                ### Aggregate additional cols
                if (!is.null(add_col_lambdas) &&
                    !purrr::is_empty(add_col_lambdas)) {
                    agg_cols <- purrr::map(names(add_col_lambdas), ~ {
                        if (is.null(add_col_lambdas[[.x]])) {
                            return(x[window[dupl_indexes[1]], ][[.x]])
                        }
                        values <- x[window[dupl_indexes], ][[.x]]
                        return(rlang::as_function(
                            add_col_lambdas[[.x]]
                        )(values))
                    })
                    x[window[dupl_indexes[1]], c(names(add_col_lambdas)) :=
                        agg_cols]
                }

                x <- x[-window[dupl_indexes[-1]], ]
                to_drop <- seq(
                    from = length(window),
                    length.out = length(dupl_indexes[-1]),
                    by = -1
                )
                window <- window[-to_drop]
            }
            ## Increment index
            if (nrow(x) == starting_rows) {
                index <- k
            } else {
                index <- utils::tail(window, n = 1) + 1
            }
        } else {
            if (produce_map == TRUE) {
                is_predicate <- purrr::map(req_vars, ~ {
                    var_before <- sym(paste0(.x, "_before"))
                    value <- x[window, get(.x)]
                    if (!is.character(value)) {
                        return(rlang::expr(!!var_before %in% !!value))
                    }
                    rlang::expr(!!var_before %chin% !!value)
                }) %>%
                    purrr::reduce(~ rlang::expr(!!.x & !!.y))
                map_recalibr[
                    eval(is_predicate),
                    c(
                        paste0(req_vars, "_after")
                    ) :=
                        purrr::map(
                            req_vars,
                            ~ x[row_to_keep][[.x]]
                        )
                ]
            }
            index <- index + 1
        }
    }
    list(recalibrated_matrix = x, map = map_recalibr)
}


# Generates a file name for recalibration maps.
#
# Unique file names include current date and timestamp.
# @keywords internal
#' @importFrom stringr str_replace_all
#
# @return A string
.generate_rec_map_filename <- function() {
    def <- "recalibration_map.tsv.gz"
    date <- lubridate::today()
    return(paste0(date, "_", def))
}

# Tries to write the recalibration map to a tsv file.
#
# @param map The recalibration map
# @param file_path The file path as a string
#' @importFrom fs dir_create path_wd path
#' @importFrom readr write_tsv
# @keywords internal
#
# @return Nothing
.write_recalibr_map <- function(map, file_path) {
    if (!fs::file_exists(file_path)) {
        ext <- fs::path_ext(file_path)
        ## if ext is empty assume it's a folder
        if (ext == "") {
            fs::dir_create(file_path)
            gen_filename <- .generate_rec_map_filename()
            tmp_filename <- fs::path(file_path, gen_filename)
        } else {
            tmp_filename <- fs::path_ext_remove(file_path)
            if (ext %in% .compressed_formats()) {
                ext <- paste(fs::path_ext(tmp_filename), ext, sep = ".")
                tmp_filename <- fs::path_ext_remove(tmp_filename)
            }
            if (!ext %in% c(
                "tsv", paste("tsv", .compressed_formats(), sep = "."),
                "csv", paste("csv", .compressed_formats(), sep = "."),
                "txt", paste("txt", .compressed_formats(), sep = ".")
            )) {
                if (getOption("ISAnalytics.verbose", TRUE) == TRUE) {
                    warn <- c("Recalibration file format unsupported",
                        i = "Writing file in 'tsv.gz' format"
                    )
                    rlang::inform(warn, class = "rec_unsupp_ext")
                }
                tmp_filename <- paste0(tmp_filename, ".tsv.gz")
            } else {
                tmp_filename <- paste0(tmp_filename, ".", ext)
            }
        }
    } else if (fs::is_dir(file_path)) {
        gen_filename <- .generate_rec_map_filename()
        tmp_filename <- fs::path(file_path, gen_filename)
    } else {
        tmp_filename <- file_path
    }
    withRestarts(
        {
            data.table::fwrite(map, file = tmp_filename, sep = "\t", na = "")
            saved_msg <- paste(
                "Recalibration map saved to: ",
                tmp_filename
            )
            if (getOption("ISAnalytics.verbose", TRUE) == TRUE) {
                rlang::inform(saved_msg)
            }
        },
        skip_write = function() {
            skip_msg <- paste(
                "Could not write recalibration map file.",
                "Skipping."
            )
            rlang::inform(skip_msg)
        }
    )
}

#### ---- Internals for analysis functions ----####


#---- USED IN : threshold_filter ----

# @keywords internal
.check_threshold_param <- function(threshold) {
    if (is.list(threshold)) {
        # List must have names
        if (is.null(names(threshold))) {
            stop(.names_list_param_err())
        }
        all_nums <- purrr::map_lgl(threshold, function(t) {
            is.numeric(t) || is.integer(t)
        })
        if (!all(all_nums)) {
            stop(.threshold_err(), call. = FALSE)
        }
    } else {
        if (!(is.numeric(threshold) || is.integer(threshold))) {
            stop(.threshold_err(), call. = FALSE)
        }
    }
}

# @keywords internal
.check_comparators_param <- function(comparators) {
    if (is.list(comparators)) {
        # List must have names
        if (is.null(names(comparators))) {
            stop(.names_list_param_err())
        }
        chk <- purrr::map_lgl(comparators, function(comp) {
            is.character(comp) &
                all(comp %in% c("<", ">", "==", "!=", ">=", "<="))
        })
        if (!all(chk)) {
            stop(.comparators_err(), call. = FALSE)
        }
    } else {
        if (!(is.character(comparators) &
            all(comparators %in% c("<", ">", "==", "!=", ">=", "<=")))) {
            stop(.comparators_err(), call. = FALSE)
        }
    }
}

# @keywords internal
.check_cols_to_compare_param <- function(cols_to_compare) {
    if (is.list(cols_to_compare)) {
        # List must have names
        if (is.null(names(cols_to_compare))) {
            stop(.names_list_param_err())
        }
        chk <- purrr::map_lgl(cols_to_compare, function(col) {
            is.character(col)
        })
        if (!all(chk)) {
            stop(paste("`cols_to_compare` elements must be character"),
                call. = FALSE
            )
        }
    } else {
        if (!is.character(cols_to_compare)) {
            stop(paste("`cols_to_compare` must be character"), call. = FALSE)
        }
    }
}

# @keywords internal
.list_names_to_mod <- function(x, list_params) {
    ### list_params is a list with 1-3 elements with names as parameters
    ### single elements are also lists
    names_equal_list <- purrr::map_lgl(list_params, function(p) {
        all(names(p) %in% names(x))
    })
    if (!all(names_equal_list)) {
        stop(.param_names_not_in_list_err())
    }
    names_same <- purrr::reduce(list_params, function(p1, p2) {
        names_same <- all(names(p1) == names(p2))
        if (names_same == FALSE) {
            stop(.param_names_not_equal_err())
        }
    })
    return(names(list_params[[1]]))
}

# @keywords internal
.check_col_correct <- function(cols_to_compare, x, names_to_mod) {
    if (is.list(cols_to_compare)) {
        purrr::walk2(
            cols_to_compare,
            names(cols_to_compare),
            function(set, name) {
                df <- x[[name]]
                if (!all(set %in% colnames(df))) {
                    rlang::abort(.missing_user_cols_list_error(set, name))
                }
                purrr::walk(set, function(col) {
                    expr <- rlang::expr(`$`(df, !!col))
                    if (!is.numeric(rlang::eval_tidy(expr)) &&
                        !is.integer(rlang::eval_tidy(expr))) {
                        rlang::abort(.non_num_user_cols_error(
                            rlang::eval_tidy(expr)
                        ))
                    }
                })
            }
        )
    } else {
        if (is.data.frame(x)) {
            if (!all(cols_to_compare %in% colnames(x))) {
                rlang::abort(.missing_user_cols_error(
                    cols_to_compare[!cols_to_compare %in% colnames(x)]
                ))
            }
            ### Check they're numeric
            all_numeric <- purrr::map_lgl(cols_to_compare, function(col) {
                expr <- rlang::expr(`$`(x, !!col))
                is.numeric(rlang::eval_tidy(expr)) ||
                    is.integer(rlang::eval_tidy(expr))
            })
            if (!all(all_numeric)) {
                rlang::abort(
                    .non_num_user_cols_error(
                        cols_to_compare[!all_numeric]
                    )
                )
            }
        } else {
            purrr::walk(x[names_to_mod], function(ddf) {
                if (!all(cols_to_compare %in% colnames(ddf))) {
                    rlang::abort(
                        .missing_user_cols_error(
                            cols_to_compare[!cols_to_compare %in% colnames(ddf)]
                        )
                    )
                }
                purrr::walk(cols_to_compare, function(col) {
                    expr <- rlang::expr(`$`(ddf, !!col))
                    if (!is.numeric(rlang::eval_tidy(expr)) &&
                        !is.integer(rlang::eval_tidy(expr))) {
                        rlang::abort(
                            .non_num_user_cols_error(rlang::eval_tidy(expr))
                        )
                    }
                })
            })
        }
    }
}

# @keywords internal
.check_length_list_params <- function(list_params) {
    purrr::pwalk(list(
        list_params$threshold, list_params$cols_to_compare,
        list_params$comparators
    ), function(x, y, z) {
        if (length(x) != length(y) || length(x) != length(z) ||
            length(y) != length(z)) {
            stop(.diff_leng_args())
        }
    })
}

# @keywords internal
.check_length_vector_params <- function(vec_params, list_params) {
    ### Check lengths between vector params first
    vec_lengths <- purrr::map_dbl(vec_params, length)
    max <- max(vec_lengths)
    vec_lengths_ok <- purrr::map_lgl(vec_lengths, function(l) {
        l == 1 || l == max
    })
    if (any(vec_lengths_ok == FALSE)) {
        stop(.diff_leng_args(), call. = FALSE)
    }
    if (!is.null(list_params) || !purrr::is_empty(list_params)) {
        list_lengths_ok <- purrr::map_lgl(list_params, function(p) {
            l <- purrr::map_dbl(p, length)
            all(l == max)
        })
        if (any(list_lengths_ok == FALSE)) {
            stop(.diff_leng_args(), call. = FALSE)
        }
    }
}

# Threshold filter for data frames.
#
# @param x A data frame
# @param threshold Vector of threshold values
# @param cols_to_compare Vector of columns to compare
# @param comparators Vector of comparators
#
#' @importFrom purrr pmap_chr reduce
#' @importFrom rlang eval_tidy parse_expr .data
#' @importFrom magrittr `%>%`
#
# @keywords internal
# @return A data frame
.tf_data_frame <- function(x, threshold, cols_to_compare, comparators) {
    ### Check all columns in input are present
    if (is.list(threshold) || is.list(comparators) ||
        is.list(cols_to_compare)) {
        stop(.list_params_err())
    }
    .check_threshold_param(threshold)
    .check_comparators_param(comparators)
    .check_cols_to_compare_param(cols_to_compare)
    .check_col_correct(cols_to_compare, x, NULL)
    .check_length_vector_params(
        list(threshold, cols_to_compare, comparators),
        NULL
    )
    ### Obtain single expression for all columns
    conditions <- purrr::pmap_chr(
        list(
            cols_to_compare,
            comparators,
            threshold
        ),
        function(q, c, t) {
            temp <- paste(q, c, t)
            paste0(".data$", temp)
        }
    )
    single_predicate_string <- purrr::reduce(conditions, function(c1, c2) {
        paste0(c1, ", ", c2)
    })
    single_predicate_string <- paste0(
        "x %>% dplyr::filter(",
        single_predicate_string,
        ")"
    )
    filtered <- rlang::eval_tidy(rlang::parse_expr(single_predicate_string))
    return(filtered)
}

# Threshold filter for lists.
#
# @param x The list of data frames
# @param threshold Vector or named list for threshold values
# @param cols_to_compare Vector or named list for columns to
# compare
# @param comparators Vector or named list for comparators
#
#' @importFrom purrr map_lgl is_empty map pmap_chr reduce map2
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data parse_expr eval_tidy
# @keywords internal
#
# @return A list of data frames
.tf_list <- function(x, threshold, cols_to_compare, comparators) {
    ### Check list contains data frames
    contain_df <- purrr::map_lgl(x, function(el) {
        is.data.frame(el)
    })
    if (!all(contain_df)) {
        stop(paste("Elements of the list must be data frames"))
    }
    ### Check the nature of the parameters (list or vector)
    list_has_names <- !is.null(names(x))
    params <- list(
        threshold = threshold, comparators = comparators,
        cols_to_compare = cols_to_compare
    )
    param_are_lists <- purrr::map_lgl(params, is.list)
    ### With unnamed lists only vector parameters allowed
    if (any(param_are_lists) & list_has_names == FALSE) {
        stop(.unnamed_list_err())
    }
    ### Check single params
    .check_threshold_param(threshold)
    .check_cols_to_compare_param(cols_to_compare)
    .check_comparators_param(comparators)
    list_params <- params[param_are_lists]
    names_to_mod <- if (!purrr::is_empty(list_params)) {
        .list_names_to_mod(x, list_params)
    } else {
        names(x)
    }
    if (all(param_are_lists)) {
        .check_length_list_params(params)
    }
    if (any(param_are_lists == FALSE)) {
        .check_length_vector_params(params[!param_are_lists], list_params)
    }
    .check_col_correct(cols_to_compare, x, names_to_mod)
    ### Filter each
    if (!list_has_names) {
        filtered <- purrr::map(x, function(df) {
            conditions <- purrr::pmap_chr(
                list(
                    cols_to_compare,
                    comparators,
                    threshold
                ),
                function(q, c, t) {
                    temp <- paste(q, c, t)
                    paste0(".data$", temp)
                }
            )
            single_predicate_string <- purrr::reduce(
                conditions,
                function(c1, c2) {
                    paste0(c1, ", ", c2)
                }
            )
            single_predicate_string <- paste0(
                "df %>% dplyr::filter(",
                single_predicate_string,
                ")"
            )
            rlang::eval_tidy(
                rlang::parse_expr(single_predicate_string)
            )
        })
        return(filtered)
    }
    filtered <- purrr::map2(x, names(x), function(df, nm) {
        if (!nm %in% names_to_mod) {
            df
        } else {
            col_set <- if (is.list(cols_to_compare)) {
                cols_to_compare[[nm]]
            } else {
                cols_to_compare
            }
            comp <- if (is.list(comparators)) {
                comparators[[nm]]
            } else {
                comparators
            }
            thr <- if (is.list(threshold)) {
                threshold[[nm]]
            } else {
                threshold
            }
            conditions <- purrr::pmap_chr(
                list(
                    col_set,
                    comp,
                    thr
                ),
                function(q, c, t) {
                    temp <- paste(q, c, t)
                    paste0(".data$", temp)
                }
            )
            single_predicate_string <- purrr::reduce(
                conditions,
                function(c1, c2) {
                    paste0(c1, ", ", c2)
                }
            )
            single_predicate_string <- paste0(
                "df %>% dplyr::filter(",
                single_predicate_string,
                ")"
            )
            rlang::eval_tidy(
                rlang::parse_expr(single_predicate_string)
            )
        }
    })
    return(filtered)
}

#---- USED IN : top_targeted_genes ----
# Internal to count the number of distinct integration sites for each gene
# - x: input df, must include mand vars and annotation vars
# - include_chr: should the distinction be made by taking into account that
#   some genes span over different chromosomes? If TRUE keeps same
#   (gene_sym, gene_strand) but different chr separated (different rows)
# - include_gene_strand: same but for gene strand
#' @importFrom data.table %chin%
#' @importFrom rlang sym
.count_distinct_is_per_gene <- function(x, include_chr,
    include_gene_strand,
    gene_sym_col,
    gene_strand_col,
    chr_col,
    mand_vars_to_check) {
    data.table::setDT(mand_vars_to_check)
    data.table::setDT(x)
    present_mand_vars <- mand_vars_to_check[eval(sym("names")) %chin% colnames(x)]
    group_key <- c(gene_sym_col)
    if (include_gene_strand) {
        group_key <- c(group_key, gene_strand_col)
    }
    if (include_chr) {
        group_key <- c(group_key, chr_col)
        present_mand_vars <- present_mand_vars[!eval(sym("names")) %chin% chr_col]
    }
    is_vars <- present_mand_vars$names
    count_by_gene <- x[, list(n_IS = dplyr::n_distinct(
        .SD[, mget(is_vars)]
    )), by = eval(group_key)]
    count_by_gene
}

#---- USED IN : CIS_grubbs, CIS_grubbs_overtime ----
# Param check for CIS functions
.cis_param_check <- function(x, genomic_annotation_file,
    grubbs_flanking_gene_bp,
    threshold_alpha,
    return_missing_as_df) {
    ## Check x has the correct structure
    stopifnot(is.data.frame(x))
    check_res <- list()
    ## Check dyn vars for required tags
    check_res$req_mand_vars <- .check_required_cols(c("chromosome", "locus"),
        vars_df = mandatory_IS_vars(TRUE),
        duplicate_politic = "error"
    )
    check_res$chrom_col <- check_res$req_mand_vars %>%
        dplyr::filter(.data$tag == "chromosome") %>%
        dplyr::pull(.data$names)
    check_res$locus_col <- check_res$req_mand_vars %>%
        dplyr::filter(.data$tag == "locus") %>%
        dplyr::pull(.data$names)
    check_res$strand_col <- if ("is_strand" %in% mandatory_IS_vars(TRUE)$tag) {
        col <- mandatory_IS_vars(TRUE) %>%
            dplyr::filter(.data$tag == "is_strand") %>%
            dplyr::pull(.data$names)
        if (col %in% colnames(x)) {
            col
        } else {
            NULL
        }
    } else {
        NULL
    }
    check_res$req_annot_col <- .check_required_cols(
        list(gene_symbol = "char", gene_strand = "char"),
        vars_df = annotation_IS_vars(TRUE),
        "error"
    )
    check_res$gene_symbol_col <- check_res$req_annot_col %>%
        dplyr::filter(.data$tag == "gene_symbol") %>%
        dplyr::pull(.data$names)
    check_res$gene_strand_col <- check_res$req_annot_col %>%
        dplyr::filter(.data$tag == "gene_strand") %>%
        dplyr::pull(.data$names)
    check_res$cols_required <- c(
        check_res$chrom_col, check_res$locus_col,
        check_res$gene_symbol_col,
        check_res$gene_strand_col
    )
    if (!all(check_res$cols_required %in% colnames(x))) {
        rlang::abort(.missing_user_cols_error(
            check_res$cols_required[!check_res$cols_required %in% colnames(x)]
        ))
    }
    # Check other parameters
    stopifnot(is.data.frame(genomic_annotation_file) ||
        is.character(genomic_annotation_file))
    genomic_annotation_file <- genomic_annotation_file[1]
    if (is.character(genomic_annotation_file) &&
        !genomic_annotation_file %in% c("hg19", "mm9")) {
        err_msg <- c("Genomic annotation file unknown",
            x = paste(
                "Since ISAnalytics 1.5.4, if provided as",
                "character vector, `genomic_annotation_file`",
                "parameter must be one of 'hg19' or 'mm9'"
            ),
            i = paste(
                "For using other genome reference files",
                "import them in the R environment and pass",
                "them to the function"
            )
        )
        rlang::abort(err_msg, class = "genomic_file_char")
    }
    if (is.character(genomic_annotation_file)) {
        gen_file <- paste0("refGenes_", genomic_annotation_file)
        utils::data(list = gen_file, envir = rlang::current_env())
        check_res$refgenes <- rlang::eval_tidy(rlang::sym(gen_file))
    } else {
        # Check annotation file format
        refgenes <- genomic_annotation_file
        if (!all(refGene_table_cols() %in% colnames(refgenes))) {
            rlang::abort(.non_standard_annotation_structure())
        }
        check_res$refgenes <- tibble::as_tibble(refgenes) %>%
            dplyr::mutate(chrom = stringr::str_replace_all(
                .data$chrom,
                "chr", ""
            ))
    }
    stopifnot(is.numeric(grubbs_flanking_gene_bp) ||
        is.integer(grubbs_flanking_gene_bp))
    check_res$grubbs_flanking_gene_bp <- grubbs_flanking_gene_bp[1]
    stopifnot(is.numeric(threshold_alpha))
    check_res$threshold_alpha <- threshold_alpha[1]
    stopifnot(is.logical(return_missing_as_df))
    check_res$return_missing_as_df <- return_missing_as_df[1]
    return(check_res)
}

# Join with refgenes
.cis_join_ref <- function(x, res_checks) {
    result <- list()
    # Join input with refgenes - inner join only
    result$joint_ref <- x %>%
        dplyr::inner_join(res_checks$refgenes %>%
            dplyr::select(
                dplyr::all_of(
                    c("name2", "chrom", "strand", "average_TxLen")
                )
            ), by = setNames(
            object = c("chrom", "strand", "name2"),
            nm = c(
                res_checks$chrom_col,
                res_checks$gene_strand_col,
                res_checks$gene_symbol_col
            )
        ))
    # Retrieve eventual missing genes
    missing_genes <- x %>%
        dplyr::anti_join(result$joint_ref,
            by = c(
                res_checks$gene_symbol_col,
                res_checks$gene_strand_col, res_checks$chrom_col
            )
        )
    missing_is_tot <- nrow(missing_genes)
    missing_genes <- missing_genes %>%
        dplyr::distinct(dplyr::across(dplyr::all_of(
            c(
                res_checks$gene_symbol_col, res_checks$gene_strand_col,
                res_checks$chrom_col
            )
        )))
    if (nrow(missing_genes) > 0 &
        getOption("ISAnalytics.verbose", TRUE) == TRUE) {
        warn_miss <- c("Warning: missing genes in refgenes table",
            i = paste(paste(
                "A total of", nrow(missing_genes),
                "genes",
                "were found in the input data but not",
                "in the refgene table. This may be caused by",
                "a mismatch in the annotation phase of",
                "the matrix. Here is a summary: "
            ),
            paste0(utils::capture.output({
                print(missing_genes, n = Inf)
            }), collapse = "\n"),
            sep = "\n"
            ),
            i = paste(
                "NOTE: missing genes will be removed from",
                "the final output! Review results carefully"
            ),
            i = paste(
                "A total of", missing_is_tot,
                "IS will be removed because of missing genes (",
                round((missing_is_tot / nrow(x)) * 100, 2), "% of",
                "total IS in input)"
            )
        )
        rlang::warn(warn_miss, class = "warn_miss_genes")
    }
    if (res_checks$return_missing_as_df) {
        result$missing_genes <- missing_genes
        result$missing_is <- list(
            absolute = missing_is_tot,
            perc = round(
                (missing_is_tot / nrow(x)) * 100, 2
            )
        )
    }
    return(result)
}

# Internal computation of CIS_grubbs test
.cis_grubb_calc <- function(x,
    grubbs_flanking_gene_bp,
    threshold_alpha,
    gene_symbol_col,
    gene_strand_col,
    chr_col,
    locus_col,
    strand_col) {
    ## -- Grouping by gene
    df_by_gene <- if (requireNamespace("psych", quietly = TRUE)) {
        x %>%
            dplyr::group_by(
                dplyr::across(
                    dplyr::all_of(c(gene_symbol_col, gene_strand_col, chr_col))
                )
            ) %>%
            dplyr::summarise(
                n = dplyr::n(),
                mean = mean(.data[[locus_col]]),
                sd = stats::sd(.data[[locus_col]], na.rm = TRUE),
                median = stats::median(.data[[locus_col]]),
                trimmed = mean(.data[[locus_col]], trim = .1),
                mad = stats::mad(.data[[locus_col]]),
                min = min(.data[[locus_col]]),
                max = max(.data[[locus_col]]),
                range = max(.data[[locus_col]]) - min(.data[[locus_col]]),
                skew = psych::skew(.data[[locus_col]]),
                kurtosis = psych::kurtosi(.data[[locus_col]]),
                n_IS_perGene = dplyr::n_distinct(
                    .data[[locus_col]]
                ),
                min_bp_integration_locus =
                    min(.data[[locus_col]]),
                max_bp_integration_locus =
                    max(.data[[locus_col]]),
                IS_span_bp = (max(.data[[locus_col]]) -
                    min(.data[[locus_col]])),
                avg_bp_integration_locus =
                    mean(.data[[locus_col]]),
                median_bp_integration_locus =
                    stats::median(.data[[locus_col]]),
                distinct_orientations = dplyr::if_else(
                    condition = is.null(strand_col),
                    true = NA_integer_,
                    false = dplyr::n_distinct(.data[[strand_col]])
                ),
                average_TxLen = .data$average_TxLen[1],
                .groups = "drop"
            )
    } else {
        x %>%
            dplyr::group_by(
                dplyr::across(
                    dplyr::all_of(c(gene_symbol_col, gene_strand_col, chr_col))
                )
            ) %>%
            dplyr::summarise(
                n = dplyr::n(),
                mean = mean(.data[[locus_col]]),
                sd = stats::sd(.data[[locus_col]], na.rm = TRUE),
                median = stats::median(.data[[locus_col]]),
                trimmed = mean(.data[[locus_col]], trim = .1),
                mad = stats::mad(.data[[locus_col]]),
                min = min(.data[[locus_col]]),
                max = max(.data[[locus_col]]),
                range = max(.data[[locus_col]]) - min(.data[[locus_col]]),
                n_IS_perGene = dplyr::n_distinct(
                    .data[[locus_col]]
                ),
                min_bp_integration_locus =
                    min(.data[[locus_col]]),
                max_bp_integration_locus =
                    max(.data[[locus_col]]),
                IS_span_bp = (max(.data[[locus_col]]) -
                    min(.data[[locus_col]])),
                avg_bp_integration_locus =
                    mean(.data[[locus_col]]),
                median_bp_integration_locus =
                    stats::median(.data[[locus_col]]),
                distinct_orientations = dplyr::if_else(
                    condition = is.null(strand_col),
                    true = NA_integer_,
                    false = dplyr::n_distinct(.data[[strand_col]])
                ),
                average_TxLen = .data$average_TxLen[1],
                .groups = "drop"
            )
    }

    n_elements <- nrow(df_by_gene)
    ### Grubbs test
    ### --- Gene Frequency
    df_bygene_withannotation <- df_by_gene %>%
        dplyr::mutate(
            raw_gene_integration_frequency =
                .data$n_IS_perGene / .data$average_TxLen,
            integration_frequency_withtolerance = (.data$n_IS_perGene /
                (.data$average_TxLen + grubbs_flanking_gene_bp)) * 1000,
            minus_log2_integration_freq_withtolerance =
                -log(x = .data$integration_frequency_withtolerance, base = 2)
        )
    ### --- z score
    z_mlif <- function(x) {
        sqrt((n_elements * (n_elements - 2) * x^2) /
            (((n_elements - 1)^2) - (n_elements * x^2)))
    }
    df_bygene_withannotation <- df_bygene_withannotation %>%
        dplyr::mutate(
            zscore_minus_log2_int_freq_tolerance =
                scale(.data$minus_log2_integration_freq_withtolerance)[, 1],
            neg_zscore_minus_log2_int_freq_tolerance =
                -.data$zscore_minus_log2_int_freq_tolerance,
            t_z_mlif = z_mlif(
                .data$neg_zscore_minus_log2_int_freq_tolerance
            )
        )
    ### --- tdist
    t_dist_2t <- function(x, deg) {
        return((1 - stats::pt(x, deg)) * 2)
    }
    df_bygene_withannotation <- df_bygene_withannotation %>%
        dplyr::mutate(
            tdist2t = t_dist_2t(.data$t_z_mlif, n_elements - 2),
            tdist_pt = stats::pt(
                q = .data$t_z_mlif,
                df = n_elements - 2
            ),
            tdist_bonferroni_default = ifelse(
                .data$tdist2t * n_elements > 1, 1,
                .data$tdist2t * n_elements
            ),
            tdist_bonferroni = stats::p.adjust(
                .data$tdist2t,
                method = "bonferroni",
                n = length(.data$tdist2t)
            ),
            tdist_fdr = stats::p.adjust(
                .data$tdist2t,
                method = "fdr",
                n = length(.data$tdist2t)
            ),
            tdist_benjamini = stats::p.adjust(
                .data$tdist2t,
                method = "BY",
                n = length(.data$tdist2t)
            )
        )
    df_bygene_withannotation <- df_bygene_withannotation %>%
        dplyr::mutate(
            tdist_positive_and_corrected =
                ifelse(
                    (.data$tdist_bonferroni_default < threshold_alpha &
                        .data$neg_zscore_minus_log2_int_freq_tolerance > 0),
                    .data$tdist_bonferroni_default,
                    NA
                ),
            tdist_positive = ifelse(
                (.data$tdist2t < threshold_alpha &
                    .data$neg_zscore_minus_log2_int_freq_tolerance > 0),
                .data$tdist2t,
                NA
            )
        )
    EM_correction_N <- length(
        df_bygene_withannotation$tdist_positive[
            !is.na(df_bygene_withannotation$tdist_positive)
        ]
    )
    df_bygene_withannotation <- df_bygene_withannotation %>%
        dplyr::mutate(
            tdist_positive_and_correctedEM =
                ifelse(
                    (.data$tdist2t * EM_correction_N <
                        threshold_alpha &
                        .data$neg_zscore_minus_log2_int_freq_tolerance > 0),
                    .data$tdist2t * EM_correction_N,
                    NA
                )
        )
    return(df_bygene_withannotation)
}

#---- USED IN : CIS_volcano_plot ----
#' @importFrom rlang arg_match abort inform current_env
#' @importFrom utils read.delim
.load_onco_ts_genes <- function(onco_db_file,
    tumor_suppressors_db_file,
    species) {
    onco_db <- if (onco_db_file == "proto_oncogenes") {
        utils::data("proto_oncogenes", envir = rlang::current_env())
        rlang::current_env()$proto_oncogenes
    } else {
        if (!file.exists(onco_db_file)) {
            onco_db_not_found <- c(paste("`onco_db_file` was not found"),
                x = "Did you provide the correct path?"
            )
            rlang::abort(onco_db_not_found)
        }
        utils::read.delim(
            file = onco_db_file, header = TRUE,
            fill = TRUE, sep = "\t",
            check.names = FALSE,
            na.strings = c("NONE", "NA", "NULL", "NaN", "")
        )
    }
    tumsup_db <- if (tumor_suppressors_db_file == "tumor_suppressors") {
        utils::data("tumor_suppressors", envir = rlang::current_env())
        rlang::current_env()$tumor_suppressors
    } else {
        if (!file.exists(tumor_suppressors_db_file)) {
            tum_db_not_found <- c("`tumor_suppressors_db_file` was not found",
                x = "Did you provide the correct path?"
            )
            rlang::abort(tum_db_not_found)
        }
        utils::read.delim(
            file = tumor_suppressors_db_file,
            header = TRUE,
            fill = TRUE, sep = "\t",
            check.names = FALSE,
            na.strings = c("NONE", "NA", "NULL", "NaN", "")
        )
    }
    specie <- rlang::arg_match(species, values = c("human", "mouse", "all"))
    specie_name <- switch(specie,
        "human" = "Homo sapiens (Human)",
        "mouse" = "Mus musculus (Mouse)",
        "all" = c(
            "Homo sapiens (Human)",
            "Mus musculus (Mouse)"
        )
    )
    if (getOption("ISAnalytics.verbose", TRUE) == TRUE) {
        load_onco_msg <- c(
            "Loading annotated genes -  species selected: ",
            paste(c(specie_name), collapse = ", ")
        )
        rlang::inform(load_onco_msg)
    }
    # Filter and merge
    onco_df <- .filter_db(onco_db, specie_name, "OncoGene")
    tumsup_df <- .filter_db(tumsup_db, specie_name, "TumorSuppressor")
    oncots_df_to_use <- .merge_onco_tumsup(onco_df, tumsup_df)
    if (getOption("ISAnalytics.verbose", TRUE) == TRUE) {
        rlang::inform(c(
            "Loading annotated genes -  done"
        ))
    }
    return(oncots_df_to_use)
}

#' @importFrom rlang .data
#' @importFrom magrittr `%>%`
.filter_db <- function(onco_db, species, output_col_name) {
    filtered_db <- onco_db %>%
        dplyr::filter(.data$Organism %in% species) %>%
        dplyr::filter(.data$Status == "reviewed" &
            !is.na(.data$`Gene names`)) %>%
        dplyr::select(.data$`Gene names`) %>%
        dplyr::distinct()

    filtered_db <- filtered_db %>%
        dplyr::mutate(`Gene names` = stringr::str_split(
            .data$`Gene names`, " "
        )) %>%
        tidyr::unnest("Gene names") %>%
        dplyr::mutate({{ output_col_name }} := .data$`Gene names`) %>%
        dplyr::rename("GeneName" = "Gene names")

    return(filtered_db)
}

#' @importFrom rlang .data
#' @importFrom magrittr `%>%`
.merge_onco_tumsup <- function(onco_df, tum_df) {
    onco_tumsup <- onco_df %>%
        dplyr::full_join(tum_df, by = "GeneName") %>%
        dplyr::mutate(Onco1_TS2 = ifelse(
            !is.na(.data$OncoGene) & !is.na(.data$TumorSuppressor),
            yes = 3,
            no = ifelse(!is.na(.data$OncoGene),
                yes = 1,
                no = ifelse(!is.na(.data$TumorSuppressor),
                    yes = 2,
                    no = NA
                )
            )
        )) %>%
        dplyr::select(-c("OncoGene", "TumorSuppressor")) %>%
        dplyr::distinct()
    return(onco_tumsup)
}

.expand_cis_df <- function(cis_grubbs_df,
    gene_sym_col,
    onco_db_file,
    tumor_suppressors_db_file,
    species,
    known_onco,
    suspicious_genes) {
    ## Load onco and ts
    oncots_to_use <- .load_onco_ts_genes(
        onco_db_file,
        tumor_suppressors_db_file,
        species
    )
    ## Join all dfs by gene
    cis_grubbs_df <- cis_grubbs_df %>%
        dplyr::left_join(oncots_to_use, by = gene_sym_col) %>%
        dplyr::left_join(known_onco, by = gene_sym_col) %>%
        dplyr::left_join(suspicious_genes, by = gene_sym_col)
    ## Add info
    cis_grubbs_df <- cis_grubbs_df %>%
        dplyr::mutate(
            KnownGeneClass = ifelse(
                is.na(.data$Onco1_TS2),
                yes = "Other",
                no = ifelse(.data$Onco1_TS2 == 1,
                    yes = "OncoGene",
                    no = "TumSuppressor"
                )
            ),
            CriticalForInsMut = ifelse(!is.na(.data$KnownClonalExpansion),
                yes = TRUE, no = FALSE
            )
        )
    return(cis_grubbs_df)
}


#---- USED IN : outliers_by_pool_fragments ----
.outlier_pool_frag_base_checks <- function(metadata,
    key,
    outlier_p_value_threshold,
    normality_test,
    normality_p_value_threshold,
    transform_log2,
    min_samples_per_pool,
    per_pool_test,
    pool_col, pcr_id_col) {
    stopifnot(is.data.frame(metadata))
    stopifnot(is.character(key))
    if (!all(key %in% colnames(metadata))) {
        rlang::abort(.missing_user_cols_error(key[!key %in%
            colnames(metadata)]),
        class = "missing_cols_key"
        )
    }
    stopifnot(is.numeric(outlier_p_value_threshold))
    stopifnot(is.logical(normality_test))
    if (normality_test) {
        stopifnot(is.numeric(normality_p_value_threshold))
    }
    stopifnot(is.logical(transform_log2))
    stopifnot(is.logical(per_pool_test))
    if (per_pool_test) {
        stopifnot(is.character(pool_col))
        stopifnot(is.numeric(min_samples_per_pool) &&
            length(min_samples_per_pool) == 1)
        if (!all(pool_col %in% colnames(metadata))) {
            rlang::abort(.missing_user_cols_error(
                pool_col[!pool_col %in% colnames(metadata)]
            ),
            class = "missing_cols_pool"
            )
        }
    }
    if (!pcr_id_col %in% colnames(metadata)) {
        rlang::abort(.missing_user_cols_error(pcr_id_col))
    }
}

.outlier_test_verify_logiop <- function(key, flag_logic, symbol_name) {
    if (length(key) > 1) {
        ## Verify logic
        stopifnot(is.character(flag_logic))
        if (length(flag_logic) > length(key) - 1) {
            flag_logic_new <- flag_logic[seq_len(length(key) - 1)]
            rlang::env_poke(
                env = rlang::caller_env(), nm = symbol_name,
                value = flag_logic_new
            )
            if (getOption("ISAnalytics.verbose", TRUE) == TRUE) {
                flag_msg <- c(paste0(
                    "'", symbol_name,
                    "' has more elements than expected"
                ),
                i = paste(
                    "The vector will be trimmed to consider",
                    "the first", length(key) - 1, "elements",
                    "only."
                )
                )
                rlang::inform(flag_msg, class = "flag_logic_long")
            }
        } else if (length(flag_logic) != length(key) - 1 &
            length(flag_logic) != 1) {
            flag_logic_new <- flag_logic[1]
            rlang::env_poke(
                env = rlang::caller_env(), nm = symbol_name,
                value = flag_logic_new
            )
            if (getOption("ISAnalytics.verbose", TRUE) == TRUE) {
                flag_logic_err <- c(paste(
                    "'", symbol_name, "' has an incorrect",
                    "amount of elements"
                ),
                paste(
                    "You should provide 1 or",
                    length(key) - 1,
                    "logical operators"
                ),
                i = "Only the first parameter will be considered"
                )
                rlang::inform(flag_logic_err,
                    class = "flag_logic_short"
                )
            }
        } else {
            flag_logic_new <- flag_logic
        }
        flag_logic_new <- toupper(flag_logic_new)
        rlang::env_poke(
            env = rlang::caller_env(), nm = symbol_name,
            value = flag_logic_new
        )
        if (!all(flag_logic_new %in% flag_logics())) {
            unknown_logi_op_err <- c(paste(
                "Unknown or unsupported logical operators:",
                paste0(flag_logic_new[!flag_logic_new %in% flag_logics()],
                    collapse = ", "
                )
            ))
            rlang::abort(unknown_logi_op_err, class = "unsupp_logi_op")
        }
    }
}

# Actual computation of statistical test on pre-filtered metadata (no NA values)
#' @importFrom data.table .N
#' @importFrom rlang sym
.pool_frag_calc <- function(meta,
    key,
    by_pool,
    normality_test,
    normality_threshold,
    pool_col,
    min_samples_per_pool,
    log2,
    pcr_id_col) {
    log2_removed_report <- NULL
    if (log2) {
        ## discard values < = 0 and report removed
        if (getOption("ISAnalytics.verbose", TRUE) == TRUE) {
            rlang::inform("Log2 transformation, removing values <= 0")
        }
        old_meta <- meta
        predicate_g0 <- purrr::map(key, ~ {
            rlang::expr(!!sym(.x) > 0)
        }) %>% purrr::reduce(~ rlang::expr(!!.x & !!.y))
        meta <- meta[eval(predicate_g0), ]
        log2_removed_report <- old_meta[!meta, on = pcr_id_col]
        log2_removed_report <- unique(
            log2_removed_report[, mget(c(pool_col, pcr_id_col, key))]
        )
    }
    if (by_pool) {
        # Group by pool
        pool_sample_count <- meta[, .N, by = eval(pool_col)]
        to_process <- pool_sample_count[
            eval(sym("N")) >= min_samples_per_pool,
            get(pool_col)
        ]
        pools_to_process <- meta[eval(sym(pool_col)) %in% to_process, ]
        split <- split(pools_to_process, by = pool_col)
        test_res <- purrr::map(
            split,
            ~ .process_pool_frag(
                .x, key, normality_test,
                normality_threshold, log2
            )
        )
        test_res <- purrr::reduce(
            test_res,
            ~ data.table::rbindlist(list(.x, .y))
        )
        test_res <- test_res[, c("processed") := TRUE]
        # Groups not processed because not enough samples
        not_to_process <- pool_sample_count[
            eval(sym("N")) < min_samples_per_pool, get(pool_col)
        ]
        non_proc <- meta[eval(sym(pool_col)) %in% not_to_process, ]
        non_proc <- non_proc[, c("processed") := FALSE]
        final <- data.table::rbindlist(list(test_res, non_proc), fill = TRUE)
        # Rows not processed because log2 requested and value <= 0
        if (exists("old_meta")) {
            log2_removed <- old_meta[!meta, on = pcr_id_col]
            log2_removed <- log2_removed[, c("processed") := FALSE]
            final <- data.table::rbindlist(list(final, log2_removed), fill = TRUE)
        }
        return(list(
            metadata = final,
            removed_zeros = log2_removed_report,
            non_proc_samples = non_proc
        ))
    } else {
        test_res <- .process_pool_frag(
            chunk = data.table::copy(meta),
            key = key,
            normality_test = normality_test,
            normality_threshold = normality_threshold,
            log2 = log2
        )
        final <- test_res[, c("processed") := TRUE]
        # Rows not processed because log2 requested and value <= 0
        if (exists("old_meta")) {
            negate_predicate_g0 <- rlang::expr(!(!!predicate_g0))
            log2_removed <- old_meta[eval(negate_predicate_g0), ]
            log2_removed <- log2_removed[, c("processed") := FALSE]
            final <- data.table::rbindlist(list(final, log2_removed), fill = TRUE)
        }
        return(list(
            metadata = final,
            removed_zeros = log2_removed_report
        ))
    }
}

# Only computes for each key the actual calculations (either on pool chunk
# or entire df)
#' @importFrom purrr map reduce flatten
#' @importFrom dplyr bind_cols
.process_pool_frag <- function(chunk, key, normality_test,
    normality_threshold,
    log2) {
    res <- purrr::map(key, ~ .tests_pool_frag(
        x = chunk[[.x]],
        suffix = .x,
        normality_test = normality_test,
        normality_threshold = normality_threshold,
        log2 = log2
    ))
    res <- purrr::reduce(
        list(og = chunk, data.table::as.data.table(purrr::flatten(res))),
        cbind
    )
    res
}

# Computation on single vectors (columns)
#' @importFrom stats dt shapiro.test
#' @importFrom tibble tibble
.tests_pool_frag <- function(x,
    suffix, normality_test,
    normality_threshold,
    log2) {
    if (log2) {
        x <- log2(x)
    }
    norm <- NA
    res <- list(NA_real_, NA_real_, NA_real_)
    z_statistic <- function() {
        zscore <- scale(x)
        zscore <- zscore[, 1]
        count <- length(x)
        tstudent <- sqrt(
            ((count * (count - 2) * (zscore)^2) /
                ((count - 1)^2 - count * (zscore)^2))
        )
        tdist <- stats::dt(tstudent, df = count - 2)
        return(list(zscore = zscore, tstudent = tstudent, tdist = tdist))
    }
    if (normality_test) {
        withCallingHandlers(
            {
                withRestarts(
                    {
                        shapiro_test <- stats::shapiro.test(x)
                        norm <- shapiro_test$p.value >= normality_threshold
                        if (norm) {
                            res <- z_statistic()
                        }
                    },
                    test_err = function() {
                        rlang::inform(
                            c("Unable to perform normality test, skipping")
                        )
                    }
                )
            },
            error = function(cnd) {
                rlang::inform(cnd$message)
                invokeRestart("test_err")
            }
        )
    } else {
        res <- z_statistic()
    }
    results <- if (log2) {
        res_colnames <- c("log2", "normality", "zscore", "tstudent", "tdist")
        res_colnames <- paste(res_colnames, suffix, sep = "_")
        res <- setNames(res, res_colnames[seq(3, 5)])
        res <- append(setNames(list(x, norm), res_colnames[c(1, 2)]), res)
        data.table::as.data.table(res)
    } else {
        res_colnames <- c("normality", "zscore", "tstudent", "tdist")
        res_colnames <- paste(res_colnames, suffix, sep = "_")
        res <- setNames(res, res_colnames[seq(2, 4)])
        res <- append(setNames(list(norm), res_colnames[1]), res)
        data.table::as.data.table(res)
    }
    return(results)
}

# Flag logic to use on a single key component
# returns a logical vector
#' @importFrom purrr pmap_lgl
.flag_cond_outliers_pool_frag <- function(proc,
    tdist,
    zscore,
    outlier_threshold) {
    # basic flag condition is tdist < out.thresh AND zscore < 0
    purrr::pmap_lgl(list(proc, tdist, zscore), function(p, t, z) {
        if (p == FALSE) {
            return(FALSE)
        }
        if (is.na(t) | is.na(z)) {
            return(FALSE)
        }
        return(t < outlier_threshold & z < 0)
    })
}

# Obtain and apply global logic on multiple columns based on the single
# key logic
#' @importFrom tibble tibble
.apply_flag_logic <- function(..., logic) {
    all_flags <- rlang::list2(...)
    partial <- NULL
    index <- 1
    for (op in logic) {
        switch(op,
            "AND" = {
                if (is.null(partial)) {
                    partial <- all_flags[[1]] & all_flags[[2]]
                    index <- 3
                } else {
                    partial <- partial & all_flags[[index]]
                    index <- index + 1
                }
            },
            "OR" = {
                if (is.null(partial)) {
                    partial <- all_flags[[1]] | all_flags[[2]]
                    index <- 3
                } else {
                    partial <- partial | all_flags[[index]]
                    index <- index + 1
                }
            },
            "XOR" = {
                if (is.null(partial)) {
                    partial <- xor(all_flags[[1]], all_flags[[2]])
                    index <- 3
                } else {
                    partial <- xor(partial, all_flags[[index]])
                    index <- index + 1
                }
            },
            "NAND" = {
                if (is.null(partial)) {
                    partial <- !(all_flags[[1]] & all_flags[[2]])
                    index <- 3
                } else {
                    partial <- !(partial & all_flags[[index]])
                    index <- index + 1
                }
            },
            "NOR" = {
                if (is.null(partial)) {
                    partial <- !(all_flags[[1]] | all_flags[[2]])
                    index <- 3
                } else {
                    partial <- !(partial | all_flags[[index]])
                    index <- index + 1
                }
            },
            "XNOR" = {
                if (is.null(partial)) {
                    partial <- !xor(all_flags[[1]], all_flags[[2]])
                    index <- 3
                } else {
                    partial <- !xor(partial, all_flags[[index]])
                    index <- index + 1
                }
            }
        )
    }
    return(partial)
}

.validate_outlier_output_format <- function(out, meta_rows, pcr_id_col) {
    format_err <- c("Wrong outlier test result format",
        x = paste(
            "Test results should follow the format",
            "described in the documentation"
        ),
        i = "See `?outlier_filter`"
    )
    if (any(!c("to_remove", pcr_id_col) %in% colnames(out))) {
        rlang::abort(format_err, class = "outlier_format_err")
    }
    if (!all(meta_rows %in% out[[pcr_id_col]])) {
        rlang::abort(format_err, class = "outlier_format_err")
    }
}

#---- USED IN : HSC_population_size_estimate ----
# Calculates population estimates (all)
#' @importFrom purrr map_lgl detect_index
#' @importFrom tidyr pivot_wider
.estimate_pop <- function(df,
    seqCount_column,
    fragmentEstimate_column,
    timepoint_column,
    annotation_cols,
    stable_timepoints,
    subject,
    tissue_col) {
    quant_cols <- if (is.null(fragmentEstimate_column)) {
        seqCount_column
    } else {
        c(seqCount_column, fragmentEstimate_column)
    }
    # --- STABLE TIMEPOINTS?
    found_stable <- purrr::map_lgl(
        stable_timepoints,
        ~ .x %in%
            as.numeric(df[[timepoint_column]])
    )
    first_stable_index <- if (length(found_stable) > 0) {
        purrr::detect_index(found_stable, ~ .x == TRUE)
    } else {
        0
    }
    stable_tps <- if (first_stable_index > 0) {
        TRUE
    } else {
        FALSE
    }
    # --- OBTAIN MATRIX (ALL TPs)
    matrix_desc <- df %>%
        dplyr::mutate(bin = 1) %>%
        dplyr::select(-dplyr::all_of(quant_cols)) %>%
        tidyr::pivot_wider(
            names_from = dplyr::all_of(c(
                "CellType",
                tissue_col,
                timepoint_column
            )),
            values_from = .data$bin,
            names_sort = TRUE,
            values_fill = 0
        ) %>%
        dplyr::select(-dplyr::all_of(c(
            mandatory_IS_vars(),
            annotation_cols
        ))) %>%
        as.matrix()
    # --- OBTAIN MATRIX (STABLE TPs)
    patient_slice_stable <- if (first_stable_index > 0) {
        first_stable <- stable_timepoints[first_stable_index]
        df %>%
            dplyr::filter(as.numeric(.data[[timepoint_column]]) >=
                first_stable) %>%
            dplyr::mutate(bin = 1) %>%
            dplyr::select(-dplyr::all_of(quant_cols)) %>%
            tidyr::pivot_wider(
                names_from = dplyr::all_of(c(
                    "CellType",
                    tissue_col,
                    timepoint_column
                )),
                values_from = .data$bin,
                names_sort = TRUE,
                values_fill = 0
            ) %>%
            dplyr::select(-dplyr::all_of(c(
                mandatory_IS_vars(),
                annotation_cols
            ))) %>%
            as.matrix()
    } else {
        NULL
    }

    timecaptures <- length(colnames(matrix_desc))
    cols_estimate_mcm <- c("Model", "abundance", "stderr")
    # --- ESTIMATES
    ## models0
    ### Estimate abundance without bias correction nor het
    #### For all tps
    models0_all <- .closed_m0_est(
        matrix_desc = matrix_desc,
        timecaptures = timecaptures,
        cols_estimate_mcm = cols_estimate_mcm,
        subject = subject, stable = FALSE
    )
    #### For the slice (first stable - last tp)
    models0_stable <- if (!is.null(patient_slice_stable) &&
        ncol(patient_slice_stable) > 1) {
        .closed_m0_est(
            matrix_desc = patient_slice_stable,
            timecaptures = length(colnames(patient_slice_stable)),
            cols_estimate_mcm = cols_estimate_mcm,
            subject = subject, stable = TRUE
        )
    } else {
        NULL
    }
    ## BC models
    ### Estimate abundance without bias correction nor het
    models_bc_all <- .closed_mthchaobc_est(
        matrix_desc = matrix_desc,
        timecaptures = timecaptures,
        cols_estimate_mcm = cols_estimate_mcm,
        subject = subject, stable = FALSE
    )
    models_bc_stable <- if (!is.null(patient_slice_stable) &&
        ncol(patient_slice_stable) > 1) {
        .closed_mthchaobc_est(
            matrix_desc = patient_slice_stable,
            timecaptures = length(colnames(patient_slice_stable)),
            cols_estimate_mcm = cols_estimate_mcm,
            subject = subject, stable = TRUE
        )
    } else {
        NULL
    }
    ## Consecutive timepoints
    ### Model m0 bc
    estimate_consecutive_m0 <- if (ncol(matrix_desc) > 1) {
        .consecutive_m0bc_est(
            matrix_desc = matrix_desc,
            cols_estimate_mcm = cols_estimate_mcm,
            subject = subject
        )
    } else {
        NULL
    }
    ### Model Mth
    estimate_consecutive_mth <- if (stable_tps & ncol(matrix_desc) > 2) {
        # - Note: pass the whole matrix, not only stable slice
        .consecutive_mth_est(
            matrix_desc = matrix_desc,
            cols_estimate_mcm = cols_estimate_mcm,
            subject = subject
        )
    } else {
        NULL
    }
    results <- models0_all %>%
        dplyr::bind_rows(models0_stable) %>%
        dplyr::bind_rows(models_bc_all) %>%
        dplyr::bind_rows(models_bc_stable) %>%
        dplyr::bind_rows(estimate_consecutive_m0) %>%
        dplyr::bind_rows(estimate_consecutive_mth)
    results
}

#' @importFrom Rcapture closedp.0
#' @importFrom purrr flatten_chr
#' @importFrom stringr str_split
#' @importFrom tibble as_tibble add_column
#' @importFrom dplyr select all_of
.closed_m0_est <- function(matrix_desc, timecaptures, cols_estimate_mcm,
    subject, stable) {
    ### Estimate abundance without bias correction nor het
    models0 <- Rcapture::closedp.0(matrix_desc,
        dfreq = FALSE,
        dtype = "hist",
        t = timecaptures,
        t0 = timecaptures,
        neg = FALSE
    )
    splitted_col_nms1 <- purrr::flatten_chr(stringr::str_split(
        colnames(matrix_desc)[1],
        pattern = "_"
    ))
    splitted_col_nmsn <- purrr::flatten_chr(stringr::str_split(
        colnames(matrix_desc)[ncol(matrix_desc)],
        pattern = "_"
    ))
    tp_string <- if (stable) {
        "Stable"
    } else {
        "All"
    }
    estimate_results <- tibble::as_tibble(models0$results,
        rownames = "Model"
    ) %>%
        dplyr::select(dplyr::all_of(cols_estimate_mcm)) %>%
        tibble::add_column(
            SubjectID = subject,
            Timepoints = tp_string,
            CellType = splitted_col_nms1[1],
            Tissue = splitted_col_nms1[2],
            TimePoint_from = as.numeric(splitted_col_nms1[3]),
            TimePoint_to = as.numeric(splitted_col_nmsn[3]),
            ModelType = "ClosedPopulation",
            ModelSetUp = "models0"
        )
    estimate_results
}

#' @importFrom Rcapture closedp.bc
#' @importFrom purrr flatten_chr
#' @importFrom stringr str_split
#' @importFrom tibble as_tibble add_column
#' @importFrom dplyr select all_of
.closed_mthchaobc_est <- function(matrix_desc, timecaptures, cols_estimate_mcm,
    subject, stable) {
    mthchaobc <- Rcapture::closedp.bc(matrix_desc,
        dfreq = FALSE,
        dtype = "hist",
        t = timecaptures,
        t0 = timecaptures,
        m = c("M0", "Mt", "Mth")
    )
    splitted_col_nms1 <- purrr::flatten_chr(stringr::str_split(
        colnames(matrix_desc)[1],
        pattern = "_"
    ))
    splitted_col_nmsn <- purrr::flatten_chr(stringr::str_split(
        colnames(matrix_desc)[ncol(matrix_desc)],
        pattern = "_"
    ))
    tp_string <- if (stable) {
        "Stable"
    } else {
        "All"
    }
    estimate_results_chaobc <- tibble::as_tibble(mthchaobc$results,
        rownames = "Model"
    ) %>%
        dplyr::select(dplyr::all_of(cols_estimate_mcm)) %>%
        tibble::add_column(
            SubjectID = subject,
            Timepoints = tp_string,
            CellType = splitted_col_nms1[1],
            Tissue = splitted_col_nms1[2],
            TimePoint_from = as.numeric(splitted_col_nms1[3]),
            TimePoint_to = as.numeric(splitted_col_nmsn[3]),
            ModelType = "ClosedPopulation",
            ModelSetUp = "mthchaobc"
        )
    estimate_results_chaobc
}

#' @importFrom Rcapture closedp.bc
#' @importFrom purrr flatten_chr map2_dfr
#' @importFrom stringr str_split
#' @importFrom tibble as_tibble add_column
#' @importFrom dplyr select all_of
.consecutive_m0bc_est <- function(matrix_desc, cols_estimate_mcm, subject) {
    ## Consecutive timepoints
    indexes <- seq(from = 1, to = ncol(matrix_desc) - 1, by = 1)
    purrr::map2_dfr(
        indexes,
        indexes + 1,
        function(t1, t2) {
            sub_matrix <- matrix_desc[, seq(from = t1, to = t2, by = 1)]
            splitted_col_nms1 <- purrr::flatten_chr(stringr::str_split(
                colnames(sub_matrix)[1],
                pattern = "_"
            ))
            splitted_col_nmsn <- purrr::flatten_chr(stringr::str_split(
                colnames(sub_matrix)[ncol(sub_matrix)],
                pattern = "_"
            ))
            patient_trend_M0 <- Rcapture::closedp.bc(sub_matrix,
                dfreq = FALSE,
                dtype = "hist",
                t = 2,
                t0 = 2,
                m = c("M0")
            )
            results <- tibble::as_tibble(patient_trend_M0$results,
                rownames = "Model"
            ) %>%
                dplyr::select(dplyr::all_of(cols_estimate_mcm)) %>%
                tibble::add_column(
                    SubjectID = subject,
                    Timepoints = "Consecutive",
                    CellType = splitted_col_nms1[1],
                    Tissue = splitted_col_nms1[2],
                    TimePoint_from = as.numeric(splitted_col_nms1[3]),
                    TimePoint_to = as.numeric(splitted_col_nmsn[3]),
                    ModelType = "ClosedPopulation",
                    ModelSetUp = "models0BC"
                )
            results
        }
    )
}

#' @importFrom Rcapture closedp.bc
#' @importFrom purrr flatten_chr map2_dfr
#' @importFrom stringr str_split
#' @importFrom tibble as_tibble add_column
#' @importFrom dplyr select all_of
.consecutive_mth_est <- function(matrix_desc, cols_estimate_mcm, subject) {
    indexes_s <- seq(from = 1, to = ncol(matrix_desc) - 2, by = 1)
    purrr::map2_dfr(
        indexes_s,
        indexes_s + 2,
        function(t1, t2) {
            sub_matrix <- matrix_desc[, seq(from = t1, to = t2, by = 1)]
            splitted_col_nms1 <- purrr::flatten_chr(stringr::str_split(
                colnames(sub_matrix)[1],
                pattern = "_"
            ))
            splitted_col_nmsn <- purrr::flatten_chr(stringr::str_split(
                colnames(sub_matrix)[ncol(sub_matrix)],
                pattern = "_"
            ))
            patient_trend_Mth <- Rcapture::closedp.bc(sub_matrix,
                dfreq = FALSE,
                dtype = "hist",
                t = 3, t0 = 3,
                m = c("Mth"),
                h = "Chao"
            )
            result <- tibble::as_tibble(patient_trend_Mth$results,
                rownames = "Model"
            ) %>%
                dplyr::select(dplyr::all_of(cols_estimate_mcm)) %>%
                tibble::add_column(
                    SubjectID = subject,
                    Timepoints = "Consecutive",
                    CellType = splitted_col_nms1[1],
                    Tissue = splitted_col_nms1[2],
                    TimePoint_from = as.numeric(splitted_col_nms1[3]),
                    TimePoint_to = as.numeric(splitted_col_nmsn[3]),
                    ModelType = "ClosedPopulation",
                    ModelSetUp = "modelMTHBC"
                )
            result
        }
    )
}

## Internal to call on each slice (a slice typically corresponding to 1 patient)
.re_agg_and_estimate <- function(df,
    metadata,
    fragmentEstimate_column,
    seqCount_column,
    tissue_col,
    timepoint_column,
    aggregation_key,
    seqCount_threshold,
    fragmentEstimate_threshold,
    cell_type,
    tissue_type,
    annotation_cols,
    subj_col,
    stable_timepoints) {
    # --- RE-AGGREGATION - by cell type, tissue and time point
    val_cols <- if (!is.null(fragmentEstimate_column)) {
        c(seqCount_column, fragmentEstimate_column)
    } else {
        seqCount_column
    }
    re_agg <- aggregate_values_by_key(
        x = df,
        association_file = metadata,
        value_cols = val_cols,
        key = c("CellType", tissue_col, timepoint_column),
        join_af_by = aggregation_key
    )
    re_agg <- re_agg %>%
        dplyr::filter(!is.na(.data$CellType))
    seqCount_column <- colnames(re_agg)[stringr::str_detect(
        colnames(re_agg),
        seqCount_column
    )]
    fragmentEstimate_column <- if (!is.null(fragmentEstimate_column)) {
        colnames(re_agg)[stringr::str_detect(
            colnames(re_agg),
            fragmentEstimate_column
        )]
    } else {
        NULL
    }
    ### Keep only sequence count value greater or equal to
    ### seqCount_threshold (or put it in AND with fe filter)
    ### (for each is)
    per_int_sums_low <- if (is.null(fragmentEstimate_column)) {
        re_agg %>%
            dplyr::group_by(dplyr::across(
                dplyr::all_of(mandatory_IS_vars())
            )) %>%
            dplyr::summarise(
                sum_sc = sum(.data[[seqCount_column]]),
                .groups = "drop"
            ) %>%
            dplyr::filter(.data$sum_sc >= seqCount_threshold)
    } else {
        re_agg %>%
            dplyr::group_by(dplyr::across(
                dplyr::all_of(mandatory_IS_vars())
            )) %>%
            dplyr::summarise(
                sum_sc = sum(.data[[seqCount_column]]),
                sum_fe = sum(.data[[fragmentEstimate_column]]),
                .groups = "drop"
            ) %>%
            dplyr::filter(
                .data$sum_sc >= seqCount_threshold,
                .data$sum_fe >= fragmentEstimate_threshold |
                    .data$sum_fe == 0
            )
    }
    re_agg <- re_agg %>%
        dplyr::semi_join(per_int_sums_low, by = mandatory_IS_vars())
    ### Keep values whose cell type is in cell_type, tissue
    ### in tissue_type and
    ### have timepoint greater than zero
    patient_slice <- re_agg %>%
        dplyr::filter(
            stringr::str_to_upper(.data$CellType) %in% cell_type,
            stringr::str_to_upper(.data[[tissue_col]]) %in% tissue_type,
            as.numeric(.data[[timepoint_column]]) > 0
        )
    ### If selected patient does not have at least 3 distinct timepoints
    ### simply return (do not consider for population estimate)
    if (length(unique(patient_slice[[timepoint_column]])) <= 2) {
        return(NULL)
    }
    estimate <- .estimate_pop(
        df = patient_slice,
        seqCount_column = seqCount_column,
        fragmentEstimate_column = fragmentEstimate_column,
        timepoint_column = timepoint_column,
        annotation_cols = annotation_cols,
        subject = df[[subj_col]][1],
        stable_timepoints = stable_timepoints,
        tissue_col = tissue_col
    )
    return(estimate)
}

#---- USED IN : is_sharing ----
## Internal to find absolute shared number of is between an arbitrary
## number of groups
## Dots are group names in sharing df (actually a data.table)
#' @importFrom rlang sym
.find_in_common <- function(..., lookup_tbl, keep_genomic_coord) {
    groups <- as.list(...)
    in_common <- purrr::pmap(groups, function(...) {
        grps <- list(...)
        filt <- lookup_tbl[eval(sym("group_id")) %in% grps, ]
        common <- purrr::reduce(filt$is, function(l, r) {
            l[r, on = mandatory_IS_vars(), nomatch = 0]
        })
        common
    })
    if (!keep_genomic_coord) {
        return(purrr::map(in_common, ~ nrow(.x)))
    }
    return(list(purrr::map(in_common, ~ nrow(.x)), in_common))
}

## Internal, to use on each row of the combinations df.
## Expands the row with all its permutations keeping same absolute shared is
## and counts if present
## - ... : row passed as a list
## - g_names: names of the groups (g1, g2...)
## - counts: are counts present? TRUE/FALSE
.sh_row_permut <- function(..., g_names, counts) {
    og_row <- list(...)
    coord <- if ("is_coord" %in% names(og_row)) {
        og_row$is_coord <- list(og_row$is_coord)
        TRUE
    } else {
        FALSE
    }
    ids <- unlist(og_row[g_names])
    og_row <- data.table::setDT(og_row)
    # If row of all equal elements no need for permutations
    if (length(unique(ids)) == 1) {
        return(og_row)
    }
    # If elements are different
    shared_is <- og_row$shared
    if (counts) {
        count_union <- og_row$count_union
    }
    perm <- gtools::permutations(
        n = length(g_names),
        r = length(g_names),
        v = g_names,
        set = TRUE,
        repeats.allowed = FALSE
    )
    colnames(perm) <- g_names
    perm <- data.table::setDT(as.data.frame(perm))
    perm[, c("shared") := shared_is]
    if (coord) {
        perm[, c("is_coord") := list(og_row$is_coord)]
    }
    if (counts) {
        for (g in g_names) {
            count_col <- paste0("count_", g)
            perm[, c(count_col) := paste0("count_", get(g))]
        }
        perm[, count_union := count_union]
    }
    sub_with_val <- function(val) {
        unlist(purrr::map(val, ~ og_row[[.x]]))
    }
    for (g in g_names) {
        perm[, c(g) := list(sub_with_val(get(g)))]
        if (counts) {
            count_col <- paste0("count_", g)
            perm[, c(count_col) := .(sub_with_val(get(count_col)))]
        }
    }
    return(perm)
}

# Counts the number of integrations in the union of all groups of a row
.count_group_union <- function(..., col_groups, lookup_tbl) {
    dots <- list(...)
    if ("is_coord" %in% names(dots)) {
        dots$is_coord <- list(dots$is_coord)
    }
    row <- data.table::setDT(dots)
    groups_in_row <- row[1, mget(col_groups)]
    sub_lookup <- lookup_tbl[eval(rlang::sym("group_id")) %in% groups_in_row]
    count_union <- nrow(purrr::reduce(sub_lookup$is, data.table::funion))
    row[, count_union := count_union]
    return(row)
}

# Obtains lookup table for groups
.sh_obtain_lookup <- function(key, df) {
    temp <- df %>%
        dplyr::select(dplyr::all_of(c(key, mandatory_IS_vars()))) %>%
        dplyr::group_by(dplyr::across({{ key }})) %>%
        dplyr::distinct(dplyr::across(dplyr::all_of(mandatory_IS_vars())),
            .keep_all = TRUE
        ) %>%
        tidyr::unite(col = "group_id", dplyr::all_of(key))
    temp <- data.table::setDT(temp)
    temp <- temp[, list(is = list(.SD)), by = "group_id"]
    return(temp)
}

# Obtains truth table for venn diagrams. To apply to each row with pmap.
# - groups : g1, g2, g3...
# - lookup : either list (mult/mult) or df
.sh_truth_tbl_venn <- function(..., lookup, groups) {
    row <- list(...)
    if (!is.data.frame(lookup)) {
        # If lookup is a list and not single df
        retrieve_is <- function(group_name) {
            label <- row[[group_name]] # a string
            retrieved <- (lookup[[group_name]])[
                eval(rlang::sym("group_id")) == label
            ]
            retrieved[, c("group_id") := paste0(
                get("group_id"), "(", group_name, ")"
            )]
            retrieved <- retrieved %>%
                tidyr::unnest(.data$is) %>%
                tidyr::unite(
                    col = "int_id",
                    dplyr::all_of(mandatory_IS_vars())
                )
            data.table::setDT(retrieved)
        }
        retrieved_iss <- purrr::map(groups, retrieve_is)
        retrieved_iss <- data.table::rbindlist(retrieved_iss)
        retrieved_iss[, c("observed") := TRUE]
        truth_tbl <- data.table::dcast(retrieved_iss, int_id ~ group_id,
            value.var = "observed",
            fill = FALSE
        )
        return(truth_tbl)
    } else {
        labels <- unlist(row[groups])
        retrieved <- lookup[eval(rlang::sym("group_id")) %in% labels]
        retrieved <- retrieved %>%
            tidyr::unnest(.data$is) %>%
            tidyr::unite(
                col = "int_id",
                dplyr::all_of(mandatory_IS_vars())
            )
        data.table::setDT(retrieved)
        retrieved[, c("observed") := TRUE]
        truth_tbl <- data.table::dcast(retrieved, int_id ~ group_id,
            value.var = "observed",
            fill = FALSE
        )
        return(truth_tbl)
    }
}

## Computes sharing table for single input df and single key
## -key: name of the columns to do a group by
## -minimal: true or false, if true computes ONLY combinations and not
## all permutations
## -n_comp: number of comparisons (2-way sharing, 3-way sharing...).
## Should be less or equal to the distinct number of groups
## -is_count: keep the counts in the table?
## -rel_sharing: compute relative sharing? (on g_i & on_union)
## -include_self_comp: include rows with the same group for each comparison?
## useful for heatmaps. Sharing for this rows is always 100%
.sharing_singledf_single_key <- function(df, key, minimal, n_comp,
    is_count, rel_sharing,
    include_self_comp,
    keep_genomic_coord,
    venn) {
    temp <- .sh_obtain_lookup(key, df)
    # Check n and k
    if (nrow(temp) < n_comp) {
        n_comp <- nrow(temp)
        if (getOption("ISAnalytics.verbose", TRUE) == TRUE) {
            warn_msg <- c("Number of requested comparisons too big",
                i = paste(
                    "The number of requested comparisons",
                    "is greater than the number of groups.",
                    "Reducing comparisons to the biggest value",
                    "allowed"
                )
            )
            rlang::inform(warn_msg)
        }
    }
    if (getOption("ISAnalytics.verbose", TRUE) == TRUE) {
        rlang::inform("Calculating combinations...")
    }
    group_comb <- gtools::combinations(
        n = length(temp$group_id),
        r = n_comp,
        v = temp$group_id,
        set = TRUE,
        repeats.allowed = FALSE
    )
    cols <- paste0("g", seq_len(ncol(group_comb)))
    colnames(group_comb) <- cols
    group_comb <- data.table::setDT(as.data.frame(group_comb))
    is_counts <- temp %>%
        dplyr::mutate(count = purrr::map_int(.data$is, ~ nrow(.x))) %>%
        dplyr::select(-.data$is)
    sharing_df <- if (!keep_genomic_coord) {
        group_comb[, c("shared") := .find_in_common(.SD,
            lookup_tbl = temp,
            keep_genomic_coord = keep_genomic_coord
        ),
        .SDcols = cols
        ]
    } else {
        group_comb[, c("shared", "is_coord") := .find_in_common(.SD,
            lookup_tbl = temp,
            keep_genomic_coord = keep_genomic_coord
        ),
        .SDcols = cols
        ]
    }

    if (include_self_comp) {
        ## Calculate groups with equal components
        if (getOption("ISAnalytics.verbose", TRUE) == TRUE) {
            rlang::inform("Calculating self groups (requested)...")
        }
        self_rows <- purrr::map2_df(
            is_counts$group_id, is_counts$count,
            function(x, y) {
                row_ls <- as.list(setNames(rep_len(
                    x,
                    length.out = n_comp
                ), nm = cols))
                row <- data.table::setDT(row_ls)
                row[, c("shared") := y]
                if (keep_genomic_coord) {
                    row[, c("is_coord") := list(temp[
                        eval(rlang::sym("group_id")) == x
                    ]$is)]
                }
                return(row)
            }
        )
        sharing_df <- data.table::rbindlist(list(
            self_rows,
            sharing_df
        ))
    }
    group_cols <- colnames(sharing_df)[!colnames(sharing_df) %in% c(
        "shared",
        "is_coord"
    )]
    if (is_count || rel_sharing) {
        ## Add counts -groups
        for (col in group_cols) {
            count_col_name <- paste0("count_", col)
            sharing_df <- sharing_df[
                dplyr::rename(
                    is_counts,
                    !!col := .data$group_id,
                    !!count_col_name := .data$count
                ),
                on = col,
                nomatch = 0
            ]
        }
        ## Add counts -union
        sharing_df <- purrr::pmap_df(sharing_df, .count_group_union,
            col_groups = group_cols,
            lookup_tbl = temp
        )
    }
    if (!minimal) {
        if (getOption("ISAnalytics.verbose", TRUE) == TRUE) {
            rlang::inform("Calculating permutations (requested)...")
        }
        sharing_df <- purrr::pmap_df(sharing_df, .sh_row_permut,
            g_names = group_cols,
            counts = is_count || rel_sharing
        )
    }
    if (rel_sharing) {
        # for groups
        for (col in group_cols) {
            rel_col_name <- paste0("on_", col)
            count_col_name <- paste0("count_", col)
            sharing_df <- sharing_df %>%
                dplyr::mutate(!!rel_col_name := (.data$shared /
                    .data[[count_col_name]]) * 100)
        }
        # for union
        sharing_df <- sharing_df %>%
            dplyr::mutate(on_union = (.data$shared /
                .data$count_union) * 100)
    }
    if (!is_count) {
        sharing_df <- sharing_df %>%
            dplyr::select(!dplyr::contains("count_"))
    }
    if (venn) {
        sharing_df <- sharing_df %>%
            dplyr::mutate(truth_tbl_venn = purrr::pmap(.,
                .sh_truth_tbl_venn,
                lookup = temp,
                groups = group_cols
            ))
    }
    sharing_df
}

## Computes sharing table for single input df and multiple keys
## -keys: name of the columns to do a group by (it is a named list)
## -minimal: true or false, if true computes ONLY combinations and not
## all permutations
## -is_count: keep the counts in the table?
## -rel_sharing: compute relative sharing? (on g_i & on_union)
.sharing_singledf_mult_key <- function(df, keys,
    minimal, is_count,
    rel_sharing, keep_genomic_coord,
    venn) {
    g_names <- names(keys)
    ## Obtain lookup table for each key
    lookup <- purrr::map(keys, ~ .sh_obtain_lookup(.x, df))
    group_labels <- purrr::map(lookup, ~ .x$group_id)
    unique_keys <- names(keys[!duplicated(keys)])
    lookup <- data.table::rbindlist(lookup[unique_keys])
    ## Obtain combinations
    combin <- data.table::setDT(purrr::cross_df(group_labels))
    sharing_df <- if (!keep_genomic_coord) {
        combin[, c("shared") := .find_in_common(.SD,
            lookup_tbl = lookup,
            keep_genomic_coord = keep_genomic_coord
        ),
        .SDcols = g_names
        ]
    } else {
        combin[, c("shared", "is_coord") := .find_in_common(.SD,
            lookup_tbl = lookup,
            keep_genomic_coord = keep_genomic_coord
        ),
        .SDcols = g_names
        ]
    }

    if (is_count || rel_sharing) {
        is_counts <- lookup %>%
            dplyr::mutate(count = purrr::map_int(.data$is, ~ nrow(.x))) %>%
            dplyr::select(-.data$is)
        ## Add counts -groups
        for (col in g_names) {
            count_col_name <- paste0("count_", col)
            sharing_df <- sharing_df[
                dplyr::rename(
                    is_counts,
                    !!col := .data$group_id,
                    !!count_col_name := .data$count
                ),
                on = col,
                nomatch = 0
            ]
        }
        ## Add counts -union
        sharing_df <- purrr::pmap_df(sharing_df, .count_group_union,
            col_groups = g_names,
            lookup_tbl = lookup
        )
    }
    if (!minimal) {
        if (getOption("ISAnalytics.verbose", TRUE) == TRUE) {
            rlang::inform("Calculating permutations (requested)...")
        }
        sharing_df <- purrr::pmap_df(sharing_df, .sh_row_permut,
            g_names = g_names,
            counts = is_count || rel_sharing
        )
    }
    if (rel_sharing) {
        # for groups
        for (col in g_names) {
            rel_col_name <- paste0("on_", col)
            count_col_name <- paste0("count_", col)
            sharing_df <- sharing_df %>%
                dplyr::mutate(!!rel_col_name := (.data$shared /
                    .data[[count_col_name]]) * 100)
        }
        # for union
        sharing_df <- sharing_df %>%
            dplyr::mutate(on_union = (.data$shared /
                .data$count_union) * 100)
    }
    if (!is_count) {
        sharing_df <- sharing_df %>%
            dplyr::select(!dplyr::contains("count_"))
    }
    if (venn) {
        sharing_df <- sharing_df %>%
            dplyr::mutate(truth_tbl_venn = purrr::pmap(.,
                .sh_truth_tbl_venn,
                lookup = lookup,
                groups = g_names
            ))
    }
    sharing_df
}

## Computes sharing table for mult input dfs and single key
.sharing_multdf_single_key <- function(dfs, key, minimal,
    is_count, rel_sharing, keep_genomic_coord, venn) {
    if (is.null(names(dfs))) {
        g_names <- paste0("g", seq_len(length(dfs)))
        names(dfs) <- g_names
    }
    lookups <- purrr::map(dfs, ~ .sh_obtain_lookup(key, .x))
    group_labels <- purrr::map(lookups, ~ .x$group_id)
    lookup <- data.table::rbindlist(lookups)
    ## Obtain combinations
    combin <- data.table::setDT(purrr::cross_df(group_labels))
    sharing_df <- if (!keep_genomic_coord) {
        combin[, c("shared") := .find_in_common(.SD,
            lookup_tbl = lookup,
            keep_genomic_coord = keep_genomic_coord
        ),
        .SDcols = names(dfs)
        ]
    } else {
        combin[, c("shared", "is_coord") := .find_in_common(.SD,
            lookup_tbl = lookup,
            keep_genomic_coord = keep_genomic_coord
        ),
        .SDcols = names(dfs)
        ]
    }

    if (is_count || rel_sharing) {
        is_counts <- purrr::map2(lookups, names(lookups), function(x, y) {
            count_col_name <- paste0("count_", y)
            x %>%
                dplyr::mutate(!!count_col_name := purrr::map_int(
                    .data$is, ~ nrow(.x)
                )) %>%
                dplyr::select(-.data$is) %>%
                dplyr::rename(!!y := .data$group_id)
        })
        ## Add counts -groups
        for (col in names(dfs)) {
            cnt_tbl <- is_counts[[col]]
            sharing_df <- sharing_df[
                cnt_tbl,
                on = col,
                nomatch = 0
            ]
        }
        ## Add counts -union
        sharing_df <- purrr::pmap_df(sharing_df, .count_group_union,
            col_groups = names(dfs),
            lookup_tbl = lookup
        )
    }
    if (!minimal) {
        if (getOption("ISAnalytics.verbose", TRUE) == TRUE) {
            rlang::inform("Calculating permutations (requested)...")
        }
        sharing_df <- purrr::pmap_df(sharing_df, .sh_row_permut,
            g_names = names(dfs),
            counts = is_count || rel_sharing
        )
    }
    if (rel_sharing) {
        # for groups
        for (col in names(dfs)) {
            rel_col_name <- paste0("on_", col)
            count_col_name <- paste0("count_", col)
            sharing_df <- sharing_df %>%
                dplyr::mutate(!!rel_col_name := (.data$shared /
                    .data[[count_col_name]]) * 100)
        }
        # for union
        sharing_df <- sharing_df %>%
            dplyr::mutate(on_union = (.data$shared /
                .data$count_union) * 100)
    }
    if (!is_count) {
        sharing_df <- sharing_df %>%
            dplyr::select(!dplyr::contains("count_"))
    }
    if (venn) {
        sharing_df <- sharing_df %>%
            dplyr::mutate(truth_tbl_venn = purrr::pmap(.,
                .sh_truth_tbl_venn,
                lookup = lookups,
                groups = names(dfs)
            ))
    }
    sharing_df
}

## Computes sharing table for mult input dfs and mult key
.sharing_multdf_mult_key <- function(dfs, keys, minimal,
    is_count, rel_sharing, keep_genomic_coord, venn) {
    lookups <- purrr::map2(dfs, keys, ~ .sh_obtain_lookup(.y, .x)) %>%
        purrr::set_names(names(keys))
    group_labels <- purrr::map(lookups, ~ .x$group_id)
    lookup <- data.table::rbindlist(lookups)
    ## Obtain combinations
    combin <- data.table::setDT(purrr::cross_df(group_labels))
    sharing_df <- if (!keep_genomic_coord) {
        combin[, c("shared") := .find_in_common(.SD,
            lookup_tbl = lookup,
            keep_genomic_coord = keep_genomic_coord
        ),
        .SDcols = names(keys)
        ]
    } else {
        combin[, c("shared", "is_coord") := .find_in_common(.SD,
            lookup_tbl = lookup,
            keep_genomic_coord = keep_genomic_coord
        ),
        .SDcols = names(keys)
        ]
    }
    if (is_count || rel_sharing) {
        is_counts <- purrr::map2(lookups, names(lookups), function(x, y) {
            count_col_name <- paste0("count_", y)
            x %>%
                dplyr::mutate(!!count_col_name := purrr::map_int(
                    .data$is, ~ nrow(.x)
                )) %>%
                dplyr::select(-.data$is) %>%
                dplyr::rename(!!y := .data$group_id)
        })
        ## Add counts -groups
        for (col in names(keys)) {
            cnt_tbl <- is_counts[[col]]
            sharing_df <- sharing_df[
                cnt_tbl,
                on = col,
                nomatch = 0
            ]
        }
        ## Add counts -union
        sharing_df <- purrr::pmap_df(sharing_df, .count_group_union,
            col_groups = names(keys),
            lookup_tbl = lookup
        )
    }
    if (!minimal) {
        if (getOption("ISAnalytics.verbose", TRUE) == TRUE) {
            rlang::inform("Calculating permutations (requested)...")
        }
        sharing_df <- purrr::pmap_df(sharing_df, .sh_row_permut,
            g_names = names(keys),
            counts = is_count || rel_sharing
        )
    }
    if (rel_sharing) {
        # for groups
        for (col in names(keys)) {
            rel_col_name <- paste0("on_", col)
            count_col_name <- paste0("count_", col)
            sharing_df <- sharing_df %>%
                dplyr::mutate(!!rel_col_name := (.data$shared /
                    .data[[count_col_name]]) * 100)
        }
        # for union
        sharing_df <- sharing_df %>%
            dplyr::mutate(on_union = (.data$shared /
                .data$count_union) * 100)
    }
    if (!is_count) {
        sharing_df <- sharing_df %>%
            dplyr::select(!dplyr::contains("count_"))
    }
    if (venn) {
        sharing_df <- sharing_df %>%
            dplyr::mutate(truth_tbl_venn = purrr::pmap(.,
                .sh_truth_tbl_venn,
                lookup = lookups,
                groups = names(keys)
            ))
    }
    sharing_df
}

#---- USED IN : iss_source ----
.assign_iss_by_tp <- function(df, timepoint_column) {
    is_vars <- if (.is_annotated(df)) {
        c(mandatory_IS_vars(), annotation_IS_vars())
    } else {
        mandatory_IS_vars()
    }
    return(df %>%
        dplyr::mutate(!!timepoint_column := as.numeric(
            .data[[timepoint_column]]
        )) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(is_vars))) %>%
        dplyr::group_modify(~ {
            if (nrow(.x) == 1) {
                return(.x)
            }
            min_tp <- min(.x[[timepoint_column]])
            return(.x %>%
                dplyr::filter(.data[[timepoint_column]] == min_tp))
        }) %>%
        dplyr::ungroup())
}

.sharing_for_source <- function(ref,
    sel,
    ref_key,
    sel_key,
    tp_col,
    subj_col) {
    if (!is.data.frame(ref)) {
        ### Workflow by subject
        common_names <- if (!all(names(ref) == names(sel))) {
            intersect(names(ref), names(sel))
        } else {
            names(ref)
        }
        if (length(common_names) == 0) {
            no_common_err <- c("No common subjects",
                x = paste(
                    "No common subjects between",
                    "reference and selection"
                ),
                i = paste(
                    "In reference:",
                    paste0(names(ref), collapse = ", "),
                    "\n",
                    "In selection: ",
                    paste0(names(sel), collapse = ", ")
                )
            )
            rlang::abort(no_common_err)
        }
        if (getOption("ISAnalytics.verbose", TRUE) == TRUE &&
            (length(common_names) < length(names(ref)) ||
                length(common_names) < length(names(sel)))) {
            all_names <- union(names(ref), names(sel))
            common_warn_msg <- c("Mismatch in subjects found",
                i = paste(
                    "Some subjects were excluded from",
                    "computations because they were",
                    "absent from reference or from",
                    "selection"
                ),
                paste(
                    "Excluded: ",
                    paste0(all_names[!all_names %in%
                        common_names],
                    collapse = ", "
                    )
                )
            )
        }
        if (.Platform$OS.type == "windows") {
            p <- BiocParallel::SnowParam(
                tasks = length(common_names),
                progressbar = getOption("ISAnalytics.verbose", TRUE),
                exportglobals = TRUE
            )
        } else {
            p <- BiocParallel::MulticoreParam(
                tasks = length(common_names),
                progressbar = getOption("ISAnalytics.verbose", TRUE),
                exportglobals = FALSE
            )
        }
        FUN <- function(subj,
    ref, sel, ref_key, sel_key,
    tp_col) {
            ref_df <- ref[[subj]]
            sel_df <- sel[[subj]]
            ref_key_min <- ref_key[ref_key != tp_col]
            ref_split <- ref_df %>%
                dplyr::group_by(dplyr::across(dplyr::all_of(ref_key_min))) %>%
                dplyr::group_split()
            quiet_sharing <- purrr::quietly(is_sharing)
            sharing <- purrr::map_df(ref_split, ~ {
                assigned_ref <- .assign_iss_by_tp(.x,
                    timepoint_column = tp_col
                )
                sh <- (quiet_sharing(assigned_ref,
                    sel_df,
                    group_keys = list(
                        g1 = ref_key,
                        g2 = sel_key
                    ),
                    keep_genomic_coord = TRUE
                ))$result
                sh <- sh %>%
                    tidyr::separate(
                        col = "g1",
                        into = paste0("g1_", ref_key),
                        remove = FALSE,
                        convert = TRUE
                    ) %>%
                    tidyr::separate(
                        col = "g2",
                        into = paste0("g2_", sel_key),
                        remove = FALSE,
                        convert = TRUE
                    )
            })
        }

        shared <- BiocParallel::bplapply(common_names, FUN,
            ref = ref,
            sel = sel,
            ref_key = ref_key,
            sel_key = sel_key,
            tp_col = tp_col,
            BPPARAM = p
        ) %>%
            purrr::set_names(common_names)
        return(shared)
    }
    ref_key_min <- ref_key[ref_key != tp_col]
    ref_split <- ref %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(ref_key_min))) %>%
        dplyr::group_split()
    quiet_sharing <- purrr::quietly(is_sharing)
    sharing <- purrr::map_df(ref_split, ~ {
        assigned_ref <- .assign_iss_by_tp(.x,
            timepoint_column = tp_col
        )
        (quiet_sharing(assigned_ref,
            sel,
            group_keys = list(
                g1 = ref_key,
                g2 = sel_key
            ),
            keep_genomic_coord = TRUE
        ))$result
    })
    sharing %>%
        tidyr::separate(
            col = "g1",
            into = paste0("g1_", ref_key),
            remove = FALSE,
            convert = TRUE
        ) %>%
        tidyr::separate(
            col = "g2",
            into = paste0("g2_", sel_key),
            remove = FALSE,
            convert = TRUE
        )
}
