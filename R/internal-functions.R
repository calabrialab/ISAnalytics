#------------------------------------------------------------------------------#
# Internal functions
#------------------------------------------------------------------------------#
## All functions in this file are NOT exported, to be used internally only.

#### ---- Internals for ISADataFrame construction ----####

#' Internal helper function for checking metadata fields.
#'
#' Checks if the specified metadata are present in the data frame.
#' @param x An ISADataFrame object
#'
#' @return FALSE if some metadata are not found in the data frame, TRUE
#' otherwise
.check_metadata <- function(x) {
    stopifnot(is.ISADataFrame(x))
    if (!all(vapply(
        X = metadata(x), FUN = is.element,
        set = colnames(x), FUN.VALUE = logical(1)
    ))) {
        return(FALSE)
    } else {
        return(TRUE)
    }
}

#' Internal helper function for checking presence of experimental data columns.
#'
#' Checks if there is at least one experimental data column
#' @param x An ISADataFrame object
#'
#' @return FALSE if no experimental columns were detected, TRUE if at least
#' one was detected
.check_atLeastOneExp <- function(x) {
    stopifnot(is.ISADataFrame(x))
    mandAndMeta <- c(mandatoryVars(x), metadata(x))
    if (length(colnames(x)) <= length(mandAndMeta)) {
        return(FALSE)
    } else {
        return(TRUE)
    }
}

#' Internal helper function for checking presence of experimental data columns
#' and verifying they're numeric.
#'
#' @param x An ISADataFrame object
#'
#' @return TRUE if at least one experimental data column of type numeric (or
#' integer) was detected, 'Warning' if one or more non-numeric
#' columns were detected which are not metadata, FALSE in all other cases
.check_nonNumdata <- function(x) {
    # checks if there are non numeric experimental data columns
    atleastone <- .check_atLeastOneExp(x)
    if (atleastone == FALSE) {
        return(FALSE)
    } else {
        nd <- .find_nonNumData(x)
        if (length(nd) > 0) {
            return("Warning")
        } else {
            return(TRUE)
        }
    }
}

#' Internal helper function to find the indexes of detected non numeric columns.
#'
#' @param x An ISADataFrame object
#'
#' @return A numeric vector
.find_nonNumData <- function(x) {
    mandAndMeta <- c(mandatoryVars(x), metadata(x))
    remCols <- colnames(x)[which(!(colnames(x) %in% mandAndMeta))]
    which(vapply(x[remCols], function(.x) {
        !(is.numeric(.x) || is.integer(.x))
    }, FUN.VALUE = logical(1)))
}

#' Internal utility function to use for fixing the metadata attribute after
#' structural changes to an ISADataFrame object.
#'
#' @param x An ISADataFrame object
#' @importFrom rlang expr eval_tidy
#' @importFrom purrr map_lgl
#'
#' @return The same ISADataFrame as input with fixed metadata
.fix_metadata <- function(x) {
    col_names <- colnames(x)[!colnames(x) %in% mandatoryVars(x)]
    non_num <- purrr::map_lgl(col_names, function(y) {
        cur_col <- rlang::expr(`$`(x, !!y))
        cur_col <- rlang::eval_tidy(cur_col)
        if (is.numeric(cur_col) | is.integer(cur_col)) {
            FALSE
        } else {
            TRUE
        }
    })
    new_meta <- col_names[which(non_num)]
    attr(x, "metadata") <- new_meta
    x
}

#### ---- Internals for matrix import ----####

#---- USED IN : import_single_Vispa2Matrix ----

#' Internal function to convert a messy matrix to a tidy data frame
#'
#' @description Uses the suite of functions provided by the
#' tidyverse to prooduce a more dense and ordered structure.
#' This function is not exported and should be called in other importing
#' functions.
#'
#' @param df ISADataFrame to convert to tidy
#'
#' @return a tidy ISADataFrame
#' @importFrom rlang .data
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr arrange all_of
#' @importFrom forcats fct_inseq as_factor
.messy_to_tidy <- function(df) {
    stopifnot(is.ISADataFrame(df))
    exp_cols <- which(!(colnames(df) %in% c(mandatoryVars(df), metadata(df))))
    isadf_tidy <- df %>%
        tidyr::pivot_longer(
            cols = dplyr::all_of(exp_cols),
            names_to = "CompleteAmplificationID",
            values_to = "Value",
            values_drop_na = TRUE
        ) %>%
        dplyr::arrange(forcats::fct_inseq(forcats::as_factor(.data$chr)))
    new_meta <- c(metadata(df), "CompleteAmplificationID")
    attr(isadf_tidy, "metadata") <- new_meta
    isadf_tidy
}

#' Internal function to auto-detect the type of IS based on the headers.
#'
#' @param df the data frame to inspect
#'
#' @return one value among:
#' * "OLD" : for old-style matrices that had only one column holding all genomic
#' coordinates
#' * "NEW_ANNOTATED" :  for the classic Vispa2 annotated matrices
#' * "NEW_NOTANN" : for Vispa2 not annotated matrices
#' * "MALFORMED" : in any other case
.auto_detect_type <- function(df) {
    if ("IS_genomicID" %in% colnames(df) &
        all(!(c("chr", "integration_locus", "strand") %in% colnames(df)))) {
        return("OLD")
    }
    if (all(c("chr", "integration_locus", "strand") %in% colnames(df))) {
        if (all(c("GeneName", "GeneStrand") %in% colnames(df))) {
            return("NEW_ANNOTATED")
        } else {
            return("NEW_NOTANN")
        }
    }
    return("MALFORMED")
}

#---- USED IN : import_association_file ----

#' Checks if the association file has the right format (standard headers).
#'
#' @param df The imported association file
#'
#' @return TRUE if the check passes, FALSE otherwise
.check_af_correctness <- function(df) {
    if (all(association_file_columns %in% colnames(df))) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}


#' Reads association file and checks if it's correct or not.
#'
#' @param path The path to the association file on disk
#' @importFrom tibble as_tibble
#'
#' @return A tibble containing the association file.
.read_and_correctness_af <- function(path) {
    stopifnot(is.character(path))
    stopifnot(file.exists(path))

    as_file <- read.csv(path,
        header = TRUE, check.names = FALSE,
        stringsAsFactors = FALSE, sep = "\t"
    )
    as_file <- tibble::as_tibble(as_file)
    # Checks if association file is correct
    correct <- .check_af_correctness(as_file)
    if (!correct) {
        stop(paste(
            "Malformed association file, could not import.",
            "Check the file or generate a new blank one with",
            "generate_blank_association_file()"
        ), call. = FALSE)
    }
    as_file
}

#' Internal function to check alignment between association file and file
#' system starting from the root. The alignment is checked at folder level.
#'
#' @param df The imported association file (data.frame or tibble)
#' @param root_folder Path to the root folder
#' @importFrom dplyr select distinct mutate bind_rows
#' @importFrom fs dir_ls
#' @importFrom purrr pmap is_empty reduce map_dbl
#' @importFrom stringr str_replace_all str_extract_all
#' @importFrom tibble tibble
#'
#' @return A data frame containing, for each ProjectID and
#' concatenatePoolIDSeqRun the
#' corresponding path on disk if found, NA otherwise.
.check_file_system_alignment <- function(df, root_folder) {
    temp_df <- df %>%
        dplyr::select(
            .data$ProjectID,
            .data$concatenatePoolIDSeqRun,
            .data$PathToFolderProjectID
        ) %>%
        dplyr::distinct()
    tree_struct <- fs::dir_ls(
        path = root_folder, recurse = TRUE,
        type = "directory"
    )
    root_regexp <- stringr::str_replace_all(root_folder, "\\/", "\\\\/")
    results_df <- purrr::pmap(
        temp_df,
        function(ProjectID, concatenatePoolIDSeqRun,
    PathToFolderProjectID) {
            pattern <- paste0(
                "^", root_regexp,
                stringr::str_replace_all(
                    PathToFolderProjectID, "\\/", "\\\\/"
                ),
                "\\/quantification\\/",
                concatenatePoolIDSeqRun, "$"
            )
            pattern <- stringr::str_replace_all(pattern,
                                                pattern = "\\\\\\/\\\\\\/",
                                                replacement = "\\\\/")
            found <- stringr::str_extract_all(tree_struct, pattern)
            found <- unlist(found)
            value <- if (purrr::is_empty(found)) {
                NA_character_
            } else {
                found
            }
            tibble::tibble(
                ProjectID = ProjectID,
                concatenatePoolIDSeqRun =
                    concatenatePoolIDSeqRun,
                Path = value
            )
        }
    )
    checker_df <- purrr::reduce(results_df, dplyr::bind_rows) %>%
        dplyr::mutate(
            Found = ifelse(!is.na(.data$Path), TRUE, FALSE),
            .before = .data$Path
        )
    checker_df
}

#' Updates the association file after the alignment check.
#'
#' The function checks if there are missing folders and updates the association
#' file by adding a column `Path` where the absolute path on disk for the
#' project and pool is found, if no path is found NA is inserted instead.
#'
#' @param as_file The tibble representing the read association_file
#' @param checks The tibble representing the results
#' of `.check_file_system_alignment`
#' @param root The root folder
#'
#' @return An updated association file with absolute paths
.update_af_after_alignment <- function(as_file, checks, root) {
    # If some some folders are missing a warning is thrown
    not_found <- checks %>% dplyr::filter(.data$Found == FALSE)
    if (nrow(not_found) > 0) {
        warning(paste0("One or more projects were not found in the file",
            "system starting from ", root, ", please check your",
            "association file for errors and/or your file system.",
            "Until you re-import the association file these ",
            "missing files will be ignored.",
            collapse = ""
        ),
        immediate. = TRUE, call. = FALSE
        )
    }
    # Finally import modified association file
    checks <- checks %>%
        dplyr::select(
            .data$ProjectID, .data$concatenatePoolIDSeqRun,
            .data$Path
        )
    as_file <- as_file %>%
        dplyr::left_join(checks, by = c("ProjectID", "concatenatePoolIDSeqRun"))
    as_file
}

#---- USED IN : import_parallel_Vispa2Matrices_interactive ----

#' Helper function to be used when importing matrices in parallel.
#'
#' @param association_file Either the path to the association file or the tibble
#' representing the imported association file (done via
#' `import_association_file`)
#' @param root If `association_file` is the path to the file, root is a single
#' string holding the path to the root folder, otherwise root is `NULL`
#'
#' @return A list of two elements: the first element is an updated version of
#' the association file with NAs removed, the second element is a widget showing
#' results of alignment checks.
.manage_association_file <- function(association_file, root) {
    # Manage association file
    if (is.character(association_file)) {
        # If it's a path to file import the association file
        association_file <- .read_and_correctness_af(association_file)
        checks <- .check_file_system_alignment(association_file, root)
        widget_checks <- .checker_widget(checks)
        association_file <- .update_af_after_alignment(
            association_file,
            checks, root
        )
        association_file <- association_file %>% dplyr::filter(!is.na(
            .data$Path
        ))
        return(list(association_file, widget_checks))
    } else {
        # If it's a tibble (file already imported) check the correctness
        correct <- ifelse(.check_af_correctness(association_file),
            ifelse("Path" %in% colnames(association_file), TRUE,
                FALSE
            ), FALSE
        )
        if (isFALSE(correct)) {
            stop(paste("Malformed association file"))
        }
        association_file <- association_file %>% dplyr::filter(!is.na(
            .data$Path
        ))
        return(list(association_file, NULL))
    }
}

#' Allows the user to choose interactively the projects to consider for import.
#'
#' @param association_file The tibble representing the imported association file
#' @importFrom dplyr distinct select filter
#' @importFrom rlang .data
#' @importFrom stringr str_split
#'
#' @return A modified version of the association file where only selected
#' projects are present
.interactive_select_projects_import <- function(association_file) {
    repeat {
        cat("Which projects would you like to import?\n")
        cat("[1] ALL", "\n", "[2] ONLY SOME\n", "[0] QUIT\n", sep = "")
        # Keep asking user until the choice is valid
        repeat {
            cat("Your choice: ")
            connection <- getOption("ISAnalytics.connection")
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
                    stop("Quitting", call. = FALSE)
                } else {
                    break
                }
            } else {
                message("Invalid choice, please choose between the above.")
            }
        }
        # If only some projects selected
        if (n_projects_to_import == 2) {
            cat(
                "Here is the list of available projects, type the indexes",
                "separated by a comma:\n",
                sep = ""
            )
            project_list <- dplyr::distinct(
                dplyr::select(association_file, .data$ProjectID)
            )$ProjectID
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
                projects_to_import <- sapply(unlist(
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
                    dplyr::filter(.data$ProjectID %in%
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

#' Simple internal helper function to handle user input for selection of number
#' of pools.
#'
#' @return Numeric representing user selection (1 for all pools, 2 for
#' only some pools, 0 to exit)
.pool_number_IN <- function() {
    cat("Which pools for each project would you like to import?\n")
    cat("[1] ALL", "[2] ONLY SOME", "[0] QUIT", sep = "\n")
    repeat {
        cat("Your choice: ")
        connection <- getOption("ISAnalytics.connection")
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

#' Simple helper interal function to handle user input for actual pool choices.
#'
#' @param indexes A vector of integer indexes available
#' @importFrom stringr str_split
#'
#' @return The user selection as a numeric vector
.pool_choices_IN <- function(indexes) {
    repeat {
        cat("\nYour choice: ")
        connection <- getOption("ISAnalytics.connection")
        if (length(connection) > 1) {
            connection <- connection[5]
        }
        to_imp <- readLines(
            con = connection,
            n = 1
        )
        to_imp <- sapply(unlist(
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

#' Allows the user to choose interactively the pools to consider for import.
#'
#' @param association_file The tibble representing the imported association file
#' @importFrom dplyr select distinct group_by bind_rows inner_join
#' @importFrom tibble tibble
#' @importFrom tidyr nest
#' @importFrom purrr map pmap reduce
#' @importFrom rlang .data
#'
#' @return A modified version of the association file where only selected
#' pools for each project are present
.interactive_select_pools_import <- function(association_file) {
    repeat {
        n_pools_to_import <- .pool_number_IN()
        if (n_pools_to_import == 2) {
            cat("Here is the list of available pools for each project,",
                "type the indexes separated by a comma or 0 to quit: \n\n",
                sep = ""
            )
            available <- association_file %>%
                dplyr::select(
                    .data$ProjectID,
                    .data$concatenatePoolIDSeqRun
                ) %>%
                dplyr::distinct() %>%
                tidyr::nest(data = c(.data$concatenatePoolIDSeqRun))
            pools_to_import <- purrr::pmap(available, function(...) {
                l <- list(...)
                current <- l["ProjectID"]
                current_data <- tibble::as_tibble(
                    purrr::flatten(l["data"])
                )
                indexes <- seq_along(current_data$concatenatePoolIDSeqRun)
                cat("ProjectID: ", current$ProjectID, "\n\n")
                cat("[0] QUIT\n")
                purrr::walk2(
                    current_data$concatenatePoolIDSeqRun,
                    indexes,
                    function(x, y) {
                        cat(paste0("[", y, "] ", x),
                            sep = "\n"
                        )
                    }
                )
                to_imp <- .pool_choices_IN(indexes)
                tibble::tibble(
                    ProjectID = current$ProjectID,
                    concatenatePoolIDSeqRun =
                        current_data$concatenatePoolIDSeqRun[to_imp]
                )
            })
            pools_to_import <- purrr::reduce(pools_to_import, function(x, y) {
                dplyr::bind_rows(x, y)
            })
            cat("\nYour choices: ", sep = "\n")
            print(pools_to_import)
            cat("\nConfirm your choices? [y/n]", sep = "")
            repeat {
                connection <- getOption("ISAnalytics.connection")
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
                        by = c("ProjectID", "concatenatePoolIDSeqRun")
                    )
                return(association_file)
            } else {
                next
            }
        } else {
            cat("\nYour choices: ", sep = "\n")
            print(association_file %>%
                dplyr::select(
                    .data$ProjectID,
                    .data$concatenatePoolIDSeqRun
                ) %>%
                dplyr::distinct())
            cat("\nConfirm your choices? [y/n]", sep = "")
            repeat {
                connection <- getOption("ISAnalytics.connection")
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

#' Updates a files_found tibble with Files_count and Anomalies column.
#'
#' @param lups A files_found tibble obtained in a lookup function. Must contain
#' the Files column (nested table quantification type and files)
#'
#' @return Updated files_found with Anomalies and Files_count columns
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

#' Looks up matrices to import given the association file and the root of the
#' file system.
#'
#' @param association_file Tibble representing the association file
#' @param quantification_type The type of quantification matrices to look for
#' (one in `quantification_types()`)
#' @param matrix_type The matrix_type to lookup (one between "annotated" or
#' "not_annotated")
#' @importFrom tibble tibble
#' @importFrom fs dir_ls as_fs_path
#' @importFrom purrr map reduce map_dbl
#' @importFrom stringr str_detect
#' @importFrom dplyr select distinct bind_rows mutate
#' @importFrom tidyr nest
#'
#' @return A tibble containing all found files, including duplicates and missing
.lookup_matrices <- function(association_file,
    quantification_type,
    matrix_type) {
    temp <- association_file %>%
        dplyr::select(
            .data$ProjectID,
            .data$concatenatePoolIDSeqRun,
            .data$Path
        ) %>%
        dplyr::distinct()
    # Regex for matrix type matching
    matrix_type_regexp <- paste0("\\.no0\\.annotated")
    # Map, for each row in temp
    lups <- purrr::pmap(temp, function(...) {
        # Obtain a tibble with all the columns in temp but a single row
        current <- tibble::tibble(...)
        # Scan file system starting from Path
        files_found <- fs::dir_ls(path = current$Path)
        # Map for each quantification type requested
        matching <- purrr::map(quantification_type, function(x) {
            # Are there files with the requested quantification type in the
            # scanned folder?
            file_regexp <- paste0("\\_", x, "\\_", "matrix")
            detected <- unlist(stringr::str_detect(files_found,
                pattern = file_regexp
            ))
            found <- if (length(files_found[detected]) > 0) {
                # If yes subset and in the subset find only the ones that match
                # the matrix type regex
                files_found <- files_found[detected]
                detected <- unlist(stringr::str_detect(files_found,
                    pattern = matrix_type_regexp
                ))
                type_match <- if (matrix_type == "annotated") {
                    # If the wanted matrix type is annotated
                    if (length(detected) > 0) {
                        # If some paths were found to contain the type regex
                        # then return the subsetted vector of file paths
                        files_found[detected]
                    } else {
                        # If no paths were found to contain the type regex
                        # return NA
                        NA_character_
                    }
                } else {
                    # If the wanted matrix is NOT annotated
                    if (length(files_found[detected]) > 0) {
                        # If some paths were found to contain the type regex
                        if (length(files_found[!detected]) > 0) {
                            # If the files_found MINUS the detected ones is
                            # not an empty set then return the set
                            files_found[!detected]
                        } else {
                            # If the files_found MINUS the detected ones is an
                            # empty set then return NA
                            NA_character_
                        }
                    } else {
                        # If there were no paths that matched the type regex
                        # simply return the original files_found
                        files_found
                    }
                }
            } else {
                NA_character_
            }
            tibble::tibble(
                Quantification_type = x,
                Files_found = fs::as_fs_path(found)
            )
        })
        matching <- purrr::reduce(matching, dplyr::bind_rows) %>%
            dplyr::mutate(Files_found = fs::as_fs_path(.data$Files_found))

        tibble::tibble(
            ProjectID = current$ProjectID,
            concatenatePoolIDSeqRun = current$concatenatePoolIDSeqRun,
            matching
        )
    })

    lups <- purrr::reduce(lups, dplyr::bind_rows) %>%
        tidyr::nest(Files = c(.data$Quantification_type, .data$Files_found))
    lups <- .trace_anomalies(lups)
    lups
}


#' Simple function to manage user input for duplicate file choices for each
#' quantification type in a single project/pool pair (use internally in
#' `.manage_anomalies_interactive`).
#'
#' @param q_types Vector of characters containing the unique quantification
#' types that are detected as duplicates
#' @param dupl The tibble containing quantification types and path to the files
#' found for a single project/pool pair
#' @importFrom dplyr filter slice bind_rows
#' @importFrom purrr map reduce
#' @importFrom rlang .data
#'
#' @return An updated tibble containing files chosen for each type
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
            connection <- getOption("ISAnalytics.connection")
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

#' Manages anomalies for files found.
#'
#' The function manages anomalies found for files after scanning appropriate
#' folders by:
#' * Removing files not found (files for which Files_count$Found == 0 and Path
#' is NA) and printing a message to notify user
#' * Removing duplicates by asking the user which files to keep for each
#' quantification type, pool and project
#'
#' @param files_found The tibble obtained via calling `.lookup_matrices`
#' @importFrom dplyr filter select rename bind_rows arrange
#' @importFrom tidyr unnest
#' @importFrom tibble as_tibble tibble
#' @importFrom purrr pmap flatten reduce
#'
#' @return A tibble containing for each project, pool and quantification type
#' the files chosen (ideally 1 for each quantification type if found, no more
#' than 1 per type)
.manage_anomalies_interactive <- function(files_found) {
    # Isolate anomalies in files found
    anomalies <- files_found %>% dplyr::filter(.data$Anomalies == TRUE)
    files_found <- files_found %>% dplyr::filter(.data$Anomalies == FALSE)
    # If there are no anomalies to fix, return chosen files
    if (nrow(anomalies) == 0) {
        files_to_import <- files_found %>%
            dplyr::select(
                .data$ProjectID,
                .data$concatenatePoolIDSeqRun, .data$Files
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
                    "ProjectID",
                    "concatenatePoolIDSeqRun"
                )])
                current_files <- tibble::as_tibble(purrr::flatten(l["Files"]))
                current_count <- tibble::as_tibble(
                    purrr::flatten(l["Files_count"])
                )
                # Find missing files first
                missing <- current_count %>% dplyr::filter(.data$Found == 0)
                if (nrow(missing) > 0) {
                    message("Some files are missing and will be ignored")
                    current_files <- current_files %>%
                        dplyr::filter(!.data$Quantification_type %in%
                            missing$Quantification_type)
                }
                # Manage duplicates
                duplicate <- current_count %>% dplyr::filter(.data$Found > 1)
                if (nrow(duplicate) > 0) {
                    message("Duplicates found for some files")
                    cat("Plese select one file for each group, type the index ",
                        "as requested\n\n",
                        sep = ""
                    )
                    cat("#### ProjectID: ", current$ProjectID, "####\n\n")
                    cat(
                        "#### Pool: ", current$concatenatePoolIDSeqRun,
                        "####\n\n"
                    )
                    current_files <- .choose_duplicates_files_interactive(
                        duplicate$Quantification_type, current_files
                    )
                }
                to_import <- tibble::tibble(
                    ProjectID = current$ProjectID,
                    concatenatePoolIDSeqRun =
                        current$concatenatePoolIDSeqRun,
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
            connection <- getOption("ISAnalytics.connection")
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
                        .data$ProjectID,
                        .data$concatenatePoolIDSeqRun, .data$Files
                    ) %>%
                    tidyr::unnest(.data$Files) %>%
                    dplyr::rename(Files_chosen = "Files_found")
                files_to_import <- files_to_import %>%
                    dplyr::bind_rows(files_found) %>%
                    dplyr::arrange(.data$ProjectID)
            }
            return(files_to_import)
        } else {
            next
        }
    }
}

#' Internal function for parallel import of a single quantification type files.
#'
#' @param q_type The quantification type (single string)
#' @param files Files_found table were absolute paths of chosen files are stored
#' @param workers Number of parallel workers
#' @importFrom dplyr filter mutate bind_rows distinct
#' @importFrom BiocParallel SnowParam MulticoreParam bptry bplapply bpstop bpok
#' @importFrom purrr is_empty reduce
#'
#' @return A single tibble with all data from matrices of same quantification
#' type in tidy format
.import_type <- function(q_type, files, workers) {
    files <- files %>% dplyr::filter(.data$Quantification_type == q_type)
    # Register backend according to platform
    if (.Platform$OS.type == "windows") {
        p <- BiocParallel::SnowParam(workers = workers, stop.on.error = FALSE)
    } else {
        p <- BiocParallel::MulticoreParam(
            workers = workers,
            stop.on.error = FALSE
        )
    }
    # Import every file
    FUN <- function(x) {
        matrix <- import_single_Vispa2Matrix(x)
    }
    suppressMessages(suppressWarnings({
        matrices <- BiocParallel::bptry(
            BiocParallel::bplapply(files$Files_chosen, FUN, BPPARAM = p)
        )
    }))
    BiocParallel::bpstop(p)
    correct <- BiocParallel::bpok(matrices)
    imported_files <- files %>% dplyr::mutate(Imported = correct)
    matrices <- matrices[correct]
    # Bind rows in single tibble for all files
    if (purrr::is_empty(matrices)) {
        return(NULL)
    }
    matrices <- purrr::reduce(matrices, function(x, y) {
        x %>%
            dplyr::bind_rows(y) %>%
            dplyr::distinct()
    })
    list(matrices, imported_files)
}

#' Internal function for importing all files for each quantification type.
#'
#' @param files_to_import The tibble containing the files to import
#' @param workers Number of parallel workers
#' @importFrom dplyr select distinct bind_rows
#' @importFrom purrr map set_names reduce flatten
#' @importFrom tibble as_tibble
#'
#' @return A named list of ISADataFrames
.parallel_import_merge <- function(files_to_import, workers) {
    # Find the actual quantification types included
    q_types <- files_to_import %>%
        dplyr::select(.data$Quantification_type) %>%
        dplyr::distinct()
    q_types <- q_types$Quantification_type
    # Import and merge for every quantification type
    imported_matrices <- purrr::map(q_types,
        .f = ~ .import_type(
            .x,
            files_to_import,
            workers
        )
    ) %>%
        purrr::set_names(q_types)
    summary_files <- purrr::map(imported_matrices, function(x) {
        tibble::as_tibble(purrr::flatten(x[2]))
    }) %>% purrr::reduce(dplyr::bind_rows)
    imported_matrices <- purrr::map(imported_matrices, function(x) {
        suppressMessages({
            ISADataFrame(purrr::flatten(x[1]))
        })
    })

    list(imported_matrices, summary_files)
}

#---- USED IN : import_parallel_Vispa2Matrices_auto ----

#' Internal function to match user defined patterns on a vector of file names.
#'
#' For each pattern specified by the user, the function tries to find a match
#' on all the file names and combines the results as a tibble in which column
#' names are the patterns and the values are TRUE if the pattern matched on the
#' element at that index or FALSE otherwise.
#'
#' @param filenames A character vector of file names
#' @param patterns A character vector of patterns to be matched
#' @importFrom tibble as_tibble_col
#' @importFrom stringr str_detect
#' @importFrom purrr map reduce
#' @importFrom dplyr bind_cols
#'
#' @return A tibble
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

#' Helper function for checking if any of the elements of the list is true.
#'
#' @param ... A list of logical values
#'
#' @return TRUE if any of the parameters is true
.any_match <- function(...) {
    l <- unlist(list(...))
    any(l)
}

#' Helper function for checking if all of the elements of the list is true.
#'
#' @param ... A list of logical values
#'
#' @return TRUE if all of the parameters is true
.all_match <- function(...) {
    l <- unlist(list(...))
    all(l)
}

#' Updates files_found tibble according to pattern and matching options.
#'
#' @param files_nested The tibble containing Quantification_type and Files_found
#' columns relative to a single project/pool pair
#' @param p_matches The tibble representing the pattern matchings resulting from
#' `pattern_matching`
#' @param matching_opt The matching option
#' @import dplyr
#' @importFrom purrr map reduce
#' @importFrom rlang .data
#' @importFrom fs as_fs_path
#'
#' @return An updated files_found tibble according to the matching option
.update_as_option <- function(files_nested, p_matches, matching_opt) {
    patterns <- colnames(p_matches)
    # Bind columns of Files nested tbl with pattern matches results
    p_matches <- files_nested %>% dplyr::bind_cols(p_matches)
    # Find the different quantification types in the files_found
    types <- p_matches %>%
        dplyr::select(.data$Quantification_type) %>%
        dplyr::distinct()
    # For each quantification type
    to_keep <- purrr::map(types$Quantification_type, function(x) {
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
    })
    to_keep <- purrr::reduce(to_keep, dplyr::bind_rows)
    to_keep
}

#' Looks up matrices to import given the association file and the root of the
#' file system.
#'
#' @inheritParams .lookup_matrices
#' @param patterns A character vector of patterns to be matched
#' @param matching_opt A single character representing the matching option (one
#' of "ANY", "ALL" or "OPTIONAL")
#' @importFrom purrr pmap flatten map
#' @importFrom tibble as_tibble
#' @importFrom stringr str_split
#' @importFrom dplyr mutate select
#' @importFrom utils tail
#'
#' @return A tibble containing all found files, including duplicates and missing
.lookup_matrices_auto <- function(association_file,
    quantification_type,
    matrix_type,
    patterns,
    matching_opt) {
    files_found <- .lookup_matrices(
        association_file,
        quantification_type,
        matrix_type
    )
    if (nrow(files_found) > 0 & !is.null(patterns)) {
        matching_patterns <- purrr::pmap(files_found, function(...) {
            l <- list(...)
            files_nested <- tibble::as_tibble(purrr::flatten(l["Files"]))
            splitted <- stringr::str_split(files_nested$Files_found, "\\/")
            filenames <- unlist(purrr::map(splitted, function(x) {
                utils::tail(x, n = 1)
            }))
            p_matches <- .pattern_matching(filenames, patterns)
            to_keep <- .update_as_option(files_nested, p_matches, matching_opt)
        })
        files_found <- files_found %>%
            dplyr::mutate(Files = matching_patterns) %>%
            dplyr::select(
                .data$ProjectID, .data$concatenatePoolIDSeqRun,
                .data$Files
            )
        files_found <- .trace_anomalies(files_found)
    }
    files_found
}

#' Manages anomalies for files found.
#'
#' The function manages anomalies found for files after scanning appropriate
#' folders by:
#' * Removing files not found (files for which Files_count$Found == 0 and Path
#' is NA) and printing a message to notify user
#' * Removing duplicates
#'
#' @param files_found The tibble obtained via calling `.lookup_matrices_auto`
#' @importFrom dplyr filter select rename bind_rows arrange
#' @importFrom rlang .data
#' @importFrom tidyr unnest
#' @importFrom tibble as_tibble tibble
#' @importFrom purrr pmap flatten
#'
#' @return A tibble containing for each project, pool and quantification type
#' the files chosen
.manage_anomalies_auto <- function(files_found) {
    # Isolate anomalies in files found
    anomalies <- files_found %>% dplyr::filter(.data$Anomalies == TRUE)
    files_found <- files_found %>% dplyr::filter(.data$Anomalies == FALSE)
    # If there are no anomalies to fix, return chosen files
    if (nrow(anomalies) == 0) {
        files_to_import <- files_found %>%
            dplyr::select(
                .data$ProjectID,
                .data$concatenatePoolIDSeqRun, .data$Files
            ) %>%
            tidyr::unnest(.data$Files) %>%
            dplyr::rename(Files_chosen = "Files_found")
        return(files_to_import)
    }
    # For each ProjectID and PoolID
    files_to_import <- purrr::pmap(anomalies, function(...) {
        l <- list(...)
        current <- tibble::as_tibble(l[c(
            "ProjectID",
            "concatenatePoolIDSeqRun"
        )])
        current_files <- tibble::as_tibble(purrr::flatten(l["Files"]))
        current_count <- tibble::as_tibble(purrr::flatten(l["Files_count"]))
        # Find missing files first
        missing <- current_count %>% dplyr::filter(.data$Found == 0)
        if (nrow(missing) > 0) {
            message("Some files are missing and will be ignored")
            current_files <- current_files %>%
                dplyr::filter(!.data$Quantification_type %in%
                    missing$Quantification_type)
        }
        # Manage duplicates
        duplicate <- current_count %>% dplyr::filter(.data$Found > 1)
        if (nrow(duplicate) > 0) {
            message(paste(
                "Duplicates found for some files: in automatic mode",
                "duplicates are not preserved - use interactive mode for more",
                "accurate file selection"
            ))
            current_files <- current_files %>%
                dplyr::filter(!.data$Quantification_type %in%
                    duplicate$Quantification_type)
        }
        to_import <- tibble::tibble(
            ProjectID = current$ProjectID,
            concatenatePoolIDSeqRun =
                current$concatenatePoolIDSeqRun,
            current_files
        )
        to_import <- to_import %>% dplyr::rename(Files_chosen = "Files_found")
    })
    # Obtain single tibble for all anomalies
    files_to_import <- purrr::reduce(files_to_import, dplyr::bind_rows)
    # Add non-anomalies
    if (nrow(files_found) > 0) {
        files_found <- files_found %>%
            dplyr::select(
                .data$ProjectID,
                .data$concatenatePoolIDSeqRun, .data$Files
            ) %>%
            tidyr::unnest(.data$Files) %>%
            dplyr::rename(Files_chosen = "Files_found")
        files_to_import <- files_to_import %>%
            dplyr::bind_rows(files_found) %>%
            dplyr::arrange(.data$ProjectID)
    }
    files_to_import
}

#### ---- Internals for HTML widgets construction ----####

#' Builds the html widget for the checker table.
#'
#' @param checker_df Tibble obtained via `.check_file_system_alignment`
#' @importFrom reactable reactable reactableTheme colDef
#' @importFrom htmltools div span h2 css browsable
#'
#' @return An html widget
.checker_widget <- function(checker_df) {
    styled_df <- reactable::reactable(
        checker_df,
        striped = TRUE,
        sortable = TRUE,
        showSortable = TRUE,
        bordered = FALSE,
        outlined = TRUE,
        searchable = TRUE,
        pagination = TRUE,
        paginationType = "simple",
        showPageSizeOptions = TRUE,
        pageSizeOptions = c(4, 8, 12),
        defaultPageSize = 5,
        showPagination = TRUE,
        resizable = TRUE,
        defaultSorted = list(Found = "asc"),
        theme = reactable::reactableTheme(
            style = list(
                fontFamily = "Calibri"
            ),
            cellStyle = list(
                display = "flex",
                flexDirection = "column",
                justifyContent = "center"
            )
        ),
        defaultColDef = reactable::colDef(
            headerStyle = list(fontSize = "18px", paddingLeft = "15px"),
            align = "left"
        ),
        columns = list(
            ProjectID = reactable::colDef(
                align = "right",
                style = list(
                    color = "#9e9e9e",
                    fontWeight = "800",
                    borderRight = "2px solid #E6E6E6"
                ),
                minWidth = 60
            ),
            concatenatePoolIDSeqRun = reactable::colDef(
                minWidth = 100,
                style = list(paddingLeft = "15px")
            ),
            Found = reactable::colDef(
                maxWidth = 100,
                align = "center",
                style = function(value) {
                    color <- if (value == TRUE) {
                        "#6afc21"
                    } else {
                        "#d61e1e"
                    }
                    list(
                        color = color, paddingLeft = "15px",
                        fontWeight = "bold"
                    )
                },
                cell = function(value) {
                    if (value == TRUE) "\u2713" else "\u2718"
                }
            ),
            Path = reactable::colDef(
                minWidth = 200,
                style = list(paddingLeft = "15px")
            )
        )
    )
    widget <- htmltools::div(
        htmltools::div(
            htmltools::h2("ALIGNMENT RESULTS"),
            htmltools::span(paste(
                "Results of alignment between file system and",
                "association file. If some folders are not found",
                "they will be ignored until the problem is fixed",
                "and the association file re-imported."
            )),
            style = htmltools::css(font.family = "Calibri"),
            styled_df
        )
    )
    htmltools::browsable(widget)
}

#' Builds the html widget for the files_found table.
#'
#' @param files_found Tibble obtained via `.lookup_matrices` or
#' `.lookup_matrices_auto`
#' @importFrom reactable reactable reactableTheme colDef
#' @importFrom htmltools div span h2 css h3 browsable
#' @importFrom dplyr select
#' @importFrom rlang .data
#'
#' @return An html widget
.files_found_widget <- function(files_found) {
    main_cols <- files_found %>% dplyr::select(
        .data$ProjectID,
        .data$concatenatePoolIDSeqRun,
        .data$Anomalies
    )
    styled_df <- reactable::reactable(
        main_cols,
        striped = TRUE,
        sortable = TRUE,
        showSortable = TRUE,
        bordered = FALSE,
        outlined = TRUE,
        searchable = TRUE,
        pagination = TRUE,
        paginationType = "simple",
        showPageSizeOptions = TRUE,
        pageSizeOptions = c(4, 8, 12),
        defaultPageSize = 5,
        showPagination = TRUE,
        resizable = TRUE,
        theme = reactable::reactableTheme(
            style = list(
                fontFamily = "Calibri"
            ),
            cellStyle = list(
                display = "flex",
                flexDirection = "column",
                justifyContent = "center"
            )
        ),
        defaultColDef = reactable::colDef(
            headerStyle = list(fontSize = "18px", paddingLeft = "15px"),
            align = "left"
        ),
        columns = list(
            ProjectID = reactable::colDef(
                align = "right",
                style = list(
                    color = "#9e9e9e",
                    fontWeight = "800",
                    borderRight = "2px solid #E6E6E6"
                ),
                minWidth = 60
            ),
            concatenatePoolIDSeqRun = reactable::colDef(
                minWidth = 100,
                style = list(paddingLeft = "15px")
            ),
            Anomalies = reactable::colDef(
                minWidth = 150,
                align = "center",
                style = function(value) {
                    color <- if (value == TRUE) {
                        "#f2cd29"
                    } else {
                        "#6afc21"
                    }
                    list(
                        color = color, paddingLeft = "15px",
                        fontWeight = "bold",
                        fontSize = "20px"
                    )
                },
                cell = function(value) {
                    if (value == TRUE) "\u26A0" else "\u2713"
                }
            )
        ),
        details = function(index) {
            files <- files_found$Files[[index]]
            count <- files_found$Files_count[[index]]
            styled_files <- reactable::reactable(
                files,
                bordered = FALSE,
                outlined = TRUE,
                resizable = TRUE,
                striped = TRUE,
                pagination = TRUE,
                defaultPageSize = 4,
                showPagination = TRUE,
                paginationType = "simple",
                defaultColDef = reactable::colDef(
                    headerStyle = list(fontSize = "18px", paddingLeft = "15px"),
                    align = "left",
                    style = list(paddingLeft = "15px"),
                    header = function(value) gsub("_", " ", value, fixed = TRUE)
                ),
                columns = list(
                    Quantification_type = reactable::colDef(
                        minWidth = 200,
                        maxWidth = 200,
                        filterable = TRUE
                    )
                )
            )
            styled_count <- reactable::reactable(
                count,
                bordered = FALSE,
                outlined = TRUE,
                resizable = TRUE,
                striped = TRUE,
                defaultColDef = reactable::colDef(
                    headerStyle = list(fontSize = "18px", paddingLeft = "15px"),
                    align = "left",
                    style = list(paddingLeft = "15px")
                ),
                columns = list(
                    Quantification_type = reactable::colDef(
                        header = function(value) {
                            gsub("_", " ", value,
                                fixed = TRUE
                            )
                        }
                    ),
                    Found = reactable::colDef(
                        cell = function(value) {
                            if (value == 1) {
                                value
                            } else {
                                if (value > 1) {
                                    paste(value, "\u2691")
                                } else {
                                    paste(value, "\u2716")
                                }
                            }
                        },
                        style = function(value) {
                            if (value == 1) {
                                color <- "black"
                                weight <- "normal"
                            } else {
                                if (value > 1) {
                                    color <- "#f2cd29"
                                    weight <- "bold"
                                } else {
                                    color <- "#d61e1e"
                                    weight <- "bold"
                                }
                            }
                            list(
                                color = color, fontWeight = weight,
                                paddingLeft = "15px"
                            )
                        }
                    )
                )
            )
            htmltools::div(
                style = paste(
                    "padding-left: 40px;",
                    "padding-right: 40px;",
                    "padding-bottom: 20px"
                ),
                htmltools::h3(paste(
                    "Summary of files count for each",
                    "quantification type"
                )),
                styled_count,
                htmltools::h3(paste(
                    "Summary of files found for each",
                    "quantification type"
                )),
                styled_files
            )
        }
    )
    widget <- htmltools::div(
        htmltools::h2("INTEGRATION MATRICES FOUND REPORT"),
        htmltools::span(paste(
            "Report of all files found for each quantification",
            "type. Click on the arrow on the left side of each",
            "row to see details."
        )),
        style = htmltools::css(font.family = "Calibri"),
        styled_df
    )
    htmltools::browsable(widget)
}

#' Builds the html widget for the files_to_import table.
#'
#' @param files_to_import Tibble obtained via `.manage_anomalies_interactive` or
#' `.manage_anomalies_auto`
#' @importFrom reactable reactable reactableTheme colDef
#' @importFrom htmltools div span h2 css browsable
#'
#' @return An html widget
.files_to_import_widget <- function(files_to_import) {
    styled_df <- reactable::reactable(
        files_to_import,
        striped = TRUE,
        sortable = TRUE,
        showSortable = TRUE,
        bordered = FALSE,
        outlined = TRUE,
        searchable = TRUE,
        pagination = TRUE,
        paginationType = "simple",
        showPageSizeOptions = TRUE,
        pageSizeOptions = c(4, 8, 12),
        defaultPageSize = 5,
        showPagination = TRUE,
        resizable = TRUE,
        theme = reactable::reactableTheme(
            style = list(
                fontFamily = "Calibri"
            ),
            cellStyle = list(
                display = "flex",
                flexDirection = "column",
                justifyContent = "center"
            )
        ),
        defaultColDef = reactable::colDef(
            headerStyle = list(
                fontSize = "18px", paddingLeft = "15px",
                display = "flex",
                flexDirection = "column",
                justifyContent = "center"
            ),
            style = list(paddingLeft = "15px"),
            align = "left",
            header = function(value) gsub("_", " ", value, fixed = TRUE)
        ),
        columns = list(
            Files_chosen = reactable::colDef(
                minWidth = 250
            ),
            ProjectID = reactable::colDef(
                filterable = TRUE
            ),
            concatenatePoolIDSeqRun = reactable::colDef(
                filterable = TRUE
            ),
            Quantification_type = reactable::colDef(
                filterable = TRUE,
                align = "center"
            )
        )
    )
    widget <- htmltools::div(
        style = htmltools::css(font.family = "Calibri"),
        htmltools::h2("SUMMARY OF FILES CHOSEN FOR IMPORT"),
        htmltools::span("Here is a summary of all files chosen for import"),
        styled_df
    )
    htmltools::browsable(widget)
}

#' Builds the html widget for the files_imported table.
#'
#' @param files_imported Tibble obtained via `.parallel_import_merge`
#' @importFrom reactable reactable reactableTheme colDef
#' @importFrom htmltools div span h2 css browsable
#'
#' @return An html widget
.files_imported_widget <- function(files_imported) {
    styled_df <- reactable::reactable(
        files_imported,
        striped = TRUE,
        sortable = TRUE,
        showSortable = TRUE,
        bordered = FALSE,
        outlined = TRUE,
        searchable = TRUE,
        pagination = TRUE,
        paginationType = "simple",
        showPageSizeOptions = TRUE,
        pageSizeOptions = c(4, 8, 12),
        defaultPageSize = 5,
        showPagination = TRUE,
        resizable = TRUE,
        theme = reactable::reactableTheme(
            style = list(
                fontFamily = "Calibri"
            ),
            cellStyle = list(
                display = "flex",
                flexDirection = "column",
                justifyContent = "center"
            )
        ),
        defaultColDef = reactable::colDef(
            headerStyle = list(
                fontSize = "18px", paddingLeft = "15px",
                display = "flex",
                flexDirection = "column",
                justifyContent = "center"
            ),
            style = list(paddingLeft = "15px"),
            align = "left",
            header = function(value) gsub("_", " ", value, fixed = TRUE)
        ),
        columns = list(
            Files_chosen = reactable::colDef(
                minWidth = 250
            ),
            ProjectID = reactable::colDef(
                filterable = TRUE
            ),
            concatenatePoolIDSeqRun = reactable::colDef(
                filterable = TRUE
            ),
            Quantification_type = reactable::colDef(
                filterable = TRUE,
                align = "center"
            ),
            Imported = reactable::colDef(
                style = function(value) {
                    color <- if (value == TRUE) {
                        "#6afc21"
                    } else {
                        "#d61e1e"
                    }
                    list(
                        paddingLeft = "15px",
                        textTransform = "uppercase",
                        color = color,
                        fontWeight = "bold"
                    )
                },
                align = "center"
            )
        )
    )
    widget <- htmltools::div(
        style = htmltools::css(font.family = "Calibri"),
        htmltools::h2("REPORT: FILES IMPORTED"),
        htmltools::span("Here is a summary of all files actually imported for
        each quantification type. If you see 'false' in the column Imported,
        some errors might have occurred and the function was unable to import
                        that matrix."),
        styled_df
    )
    htmltools::browsable(widget)
}

#' Builds the html widget for the summary table.
#'
#' @param removed Number of removed collisions
#' @param reassigned Number of re-assigned collisions
#' @param summary Summary table
#' @param tot_rows Total number rows of sequence count matrix before processing
#' @param collision_rows Total number of rows of collisions
#' @importFrom reactable reactable reactableTheme colDef
#' @importFrom htmltools div h2 h4 css browsable
#'
#' @return A widget
.summary_collisions_widget <- function(removed,
    reassigned,
    summary,
    tot_rows,
    collision_rows) {
    theme <- reactable::reactableTheme(
        style = list(
            fontFamily = "Calibri"
        ),
        cellStyle = list(
            display = "flex",
            flexDirection = "column",
            justifyContent = "center"
        )
    )

    styled_df <- reactable::reactable(
        summary,
        striped = TRUE,
        sortable = TRUE,
        showSortable = TRUE,
        bordered = FALSE,
        outlined = TRUE,
        searchable = TRUE,
        pagination = TRUE,
        paginationType = "simple",
        showPageSizeOptions = TRUE,
        pageSizeOptions = c(4, 8, 12),
        defaultPageSize = 5,
        showPagination = TRUE,
        resizable = TRUE,
        theme = theme,
        defaultColDef = reactable::colDef(
            headerStyle = list(
                fontSize = "18px", paddingLeft = "15px",
                display = "flex",
                flexDirection = "column",
                justifyContent = "center"
            ),
            style = list(paddingLeft = "15px"),
            align = "left",
            header = function(value) gsub("_", " ", value, fixed = TRUE),
            filterable = TRUE
        )
    )

    collisions_found <- tibble::tibble(
        abs_number = collision_rows,
        percentage_on_total =
            collision_rows / tot_rows
    )
    collisions_removed <- tibble::tibble(
        abs_number = removed,
        percentage_on_collisions =
            removed / collision_rows,
        percentage_on_total =
            removed / tot_rows
    )
    collisions_reassigned <- tibble::tibble(
        abs_number = reassigned,
        percentage_on_collisions =
            reassigned / collision_rows,
        percentage_on_total =
            reassigned / tot_rows
    )

    collisions_found_styled <- reactable::reactable(
        collisions_found,
        fullWidth = FALSE,
        bordered = FALSE,
        outlined = TRUE,
        theme = theme,
        defaultColDef = reactable::colDef(
            header = function(value) gsub("_", " ", value, fixed = TRUE)
        ),
        columns = list(
            percentage_on_total = reactable::colDef(
                format = reactable::colFormat(
                    percent = TRUE,
                    digits = 2
                )
            )
        )
    )

    collisions_removed_styled <- reactable::reactable(
        collisions_removed,
        fullWidth = FALSE,
        bordered = FALSE,
        outlined = TRUE,
        theme = theme,
        defaultColDef = reactable::colDef(
            header = function(value) gsub("_", " ", value, fixed = TRUE)
        ),
        columns = list(
            percentage_on_total = reactable::colDef(
                format = reactable::colFormat(
                    percent = TRUE,
                    digits = 2
                )
            ),
            percentage_on_collisions = reactable::colDef(
                format = reactable::colFormat(
                    percent = TRUE,
                    digits = 2
                )
            )
        )
    )

    collisions_reassigned_styled <- reactable::reactable(
        collisions_reassigned,
        fullWidth = FALSE,
        bordered = FALSE,
        outlined = TRUE,
        theme = theme,
        defaultColDef = reactable::colDef(
            header = function(value) gsub("_", " ", value, fixed = TRUE)
        ),
        columns = list(
            percentage_on_total = reactable::colDef(
                format = reactable::colFormat(
                    percent = TRUE,
                    digits = 2
                )
            ),
            percentage_on_collisions = reactable::colDef(
                format = reactable::colFormat(
                    percent = TRUE,
                    digits = 2
                )
            )
        )
    )

    widget <- htmltools::div(
        style = htmltools::css(font.family = "Calibri"),
        htmltools::h2("COLLISION REMOVAL SUMMARY"),
        htmltools::h4("TOTAL READS:"),
        htmltools::div(tot_rows),
        htmltools::h4("COLLISIONS FOUND:"),
        collisions_found_styled,
        htmltools::h4("REMOVED:"),
        collisions_removed_styled,
        htmltools::h4("REASSIGNED:"),
        collisions_reassigned_styled,
        htmltools::h4("SUMMARY:"),
        styled_df
    )
    htmltools::browsable(widget)
}

#### ---- Internals for collision removal ----####

#---- USED IN : remove_collisions ----

#' Checks if association file contains more information than the matrix.
#'
#' Used to notify the user that wants to know if for the projects contained in
#' the examined matrix there are additional CompleteAmplificationIDs contained
#' in the association file that weren't included in the integration matrix (for
#' example failed runs).
#'
#' @param association_file The imported association file
#' @param df The sequence count matrix to examine
#' @import dplyr
#' @importFrom purrr map reduce
#' @importFrom rlang .data
#'
#' @return A tibble containing ProjectID, PoolID and CompleteAmplificationID
#' only for additional info and only for projects already present in df
.check_same_info <- function(association_file, df) {
    joined <- association_file %>%
        dplyr::left_join(df, by = "CompleteAmplificationID")
    projects <- (joined %>% dplyr::select(.data$ProjectID) %>%
        dplyr::distinct())$ProjectID
    joined1 <- purrr::map(projects, function(x) {
        temp <- joined %>% dplyr::filter(.data$ProjectID == x)
        if (!all(is.na(temp$chr))) {
            temp %>%
                dplyr::filter(is.na(.data$chr)) %>%
                dplyr::select(
                    .data$ProjectID,
                    .data$PoolID,
                    .data$CompleteAmplificationID
                )
        }
    })
    purrr::reduce(joined1, dplyr::bind_rows)
}

#' Produces a joined tibble between the sequence count matrix and the
#' association file
#'
#' @param seq_count_df The sequence count tibble
#' @param association_file The association file tibble
#' @param date_col The date column chosen
#' @importFrom lubridate dmy
#' @import dplyr
#' @importFrom rlang .data
#'
#' @return A tibble
.join_matrix_af <- function(seq_count_df, association_file, date_col) {
    joined <- seq_count_df %>%
        dplyr::left_join(association_file, by = "CompleteAmplificationID") %>%
        dplyr::select(
            dplyr::all_of(colnames(seq_count_df)), all_of(date_col),
            .data$ProjectID, .data$PoolID, .data$SubjectID,
            .data$ReplicateNumber
        ) %>%
        dplyr::mutate(!!date_col := lubridate::dmy(`$`(.data, !!date_col)))
    joined
}

#' Identifies independent samples and separates the joined_df in collisions and
#' non-collisions
#'
#' @param joined_df The joined tibble obtained via `.join_matrix_af`
#' @import dplyr
#' @importFrom rlang .data
#'
#' @return A named list containing the splitted joined_df for collisions and
#' non-collisions
.identify_independent_samples <- function(joined_df) {
    temp <- joined_df %>%
        dplyr::select(
            dplyr::all_of(mandatoryVars(joined_df)),
            .data$ProjectID, .data$SubjectID
        ) %>%
        dplyr::group_by(.data$chr, .data$integration_locus, .data$strand) %>%
        dplyr::distinct(.data$ProjectID, .data$SubjectID, .keep_all = TRUE) %>%
        dplyr::summarise(n = dplyr::n(), .groups = "drop_last") %>%
        dplyr::ungroup() %>%
        dplyr::filter(.data$n > 1) %>%
        dplyr::select(-c(.data$n))

    non_collisions <- joined_df %>%
        dplyr::anti_join(temp, by = c("chr", "integration_locus", "strand"))
    collisions <- joined_df %>%
        dplyr::right_join(temp, by = c("chr", "integration_locus", "strand"))
    list(collisions = collisions, non_collisions = non_collisions)
}

#' Returns a polished version of collisions table for analysis.
#'
#' Polishing consists in nesting all information relative to a single
#' integration, in other words, obtaining a linked data frame for each
#' integration.
#'
#' @param collisions The collisions table obtained via
#' `.identify_independent_samples`
#' @importFrom dplyr group_by
#' @importFrom rlang .data
#' @importFrom tidyr nest
#'
#' @return A nested ISADataFrame
.obtain_nested <- function(collisions) {
    collisions %>%
        dplyr::group_by(
            .data$chr,
            .data$integration_locus, .data$strand
        ) %>%
        tidyr::nest()
}


#' Internal for date discrimination in collision removal.
#'
#' It's the first of 4 steps in the algorithm for collision removal: it tries to
#' find a single sample who has an associated date which is earlier than any
#' other. If comparison is not possible the analysis fails and returns the
#' orginal input data.
#'
#' @details NOTE: this function is meant to be used inside a mapping function
#' such as `purrr::pmap`. The function only works on data regarding a SINGLE
#' integration (triplet chr, integration_locus, strand).
#'
#' @param nest The nested table associated with a single integration
#' @param date_col The name of the date column chosen for collision removal,
#' one in \code{date_columns_coll()}
#' @importFrom rlang expr eval_tidy .data
#' @importFrom dplyr filter arrange
#'
#' @return A named list with:
#' * $data: a tibble, containing the data (unmodified or modified)
#' * $check: a logical value indicating whether the analysis was successful or
#' not (and therefore there is the need to perform the next step)
.discriminate_by_date <- function(nest, date_col) {
    dates <- rlang::expr(`$`(nest, !!date_col))
    dates <- rlang::eval_tidy(dates)
    # Test for all dates equality
    all_equal_dates <- length(unique(dates)) == 1
    if (all_equal_dates == TRUE) {
        return(list(data = nest, check = FALSE))
    }
    # If not all are equal sort them asc
    nest <- nest %>% dplyr::arrange(`$`(.data, !!date_col))
    # Test if first date is unique
    dates <- rlang::expr(`$`(nest, !!date_col))
    dates <- rlang::eval_tidy(dates)
    if (length(dates[dates %in% dates[1]]) > 1) {
        return(list(data = nest, check = FALSE))
    }
    winning_subj <- nest$SubjectID[1]
    # Filter the winning rows
    nest <- nest %>% dplyr::filter(.data$SubjectID == winning_subj)
    return(list(data = nest, check = TRUE))
}

#' Internal for replicate discrimination in collision removal.
#'
#' It's the second of 4 steps in the algorithm for collision removal:
#' grouping by independent sample (ProjectID, SubjectID) it counts the number of
#' rows (replicates) found for each group, orders them from biggest to smallest
#' and, if a single group has more rows than any other group the integration is
#' assigned to that sample, otherwise the analysis fails and returns the
#' original input to be submitted to the next step.
#'
#' @details NOTE: this function is meant to be used inside a mapping function
#' such as `purrr::pmap`. The function only works on data regarding a SINGLE
#' integration (triplet chr, integration_locus, strand).
#'
#' @param nest The nested table associated with a single integration
#' @import dplyr
#' @importFrom rlang .data
#'
#' @return A named list with:
#' * $data: a tibble, containing the data (unmodified or modified)
#' * $check: a logical value indicating whether the analysis was successful or
#' not (and therefore there is the need to perform the next step)
.discriminate_by_replicate <- function(nest) {
    temp <- nest %>%
        dplyr::group_by(.data$ProjectID, .data$SubjectID) %>%
        dplyr::summarise(n_rows = dplyr::n(), .groups = "keep") %>%
        dplyr::arrange(dplyr::desc(.data$n_rows))
    if (length(temp$n_rows) != 1 & !temp$n_rows[1] > temp$n_rows[2]) {
        return(list(data = nest, check = FALSE))
    }
    nest <- nest %>%
        dplyr::filter(
            .data$ProjectID == temp$ProjectID[1],
            .data$SubjectID == temp$SubjectID[1]
        )
    return(list(data = nest, check = TRUE))
}

#' Internal for sequence count discrimination in collision removal.
#'
#' It's the third of 4 steps in the algorithm for collision removal:
#' grouping by independent sample (ProjectID, SubjectID), it sums the value of
#' the sequence count for each group, sorts the groups from highest to lowest
#' cumulative value and then checks the ratio between the first element and the
#' second: if the ratio is > `reads_ratio` the integration is assigned to
#' the first group, otherwise the analysis fails and returns the original input.
#'
#' @details NOTE: this function is meant to be used inside a mapping function
#' such as `purrr::pmap`. The function only works on data regarding a SINGLE
#' integration (triplet chr, integration_locus, strand).
#'
#' @param nest The nested table associated with a single integration
#' @param reads_ratio The value of the ratio between sequence count values to
#' check
#' @import dplyr
#' @importFrom rlang .data
#'
#' @return A named list with:
#' * $data: a tibble, containing the data (unmodified or modified)
#' * $check: a logical value indicating whether the analysis was successful or
#' not (and therefore there is the need to perform the next step)
.discriminate_by_seqCount <- function(nest, reads_ratio) {
    temp <- nest %>%
        dplyr::group_by(.data$ProjectID, .data$SubjectID) %>%
        dplyr::summarise(sum_seqCount = sum(.data$Value), .groups = "keep") %>%
        dplyr::arrange(dplyr::desc(.data$sum_seqCount))

    ratio <- temp$sum_seqCount[1] / temp$sum_seqCount[2]
    if (!ratio > reads_ratio) {
        return(list(data = nest, check = FALSE))
    }
    nest <- nest %>%
        dplyr::filter(
            .data$ProjectID == temp$ProjectID[1],
            .data$SubjectID == temp$SubjectID[1]
        )
    return(list(data = nest, check = TRUE))
}

#' Internal function that performs four-step-check of collisions for a single
#' integration.
#'
#' @details NOTE: this function is meant to be used inside a mapping function
#' such as `purrr::pmap`. The function only works on data regarding a SINGLE
#' integration (triplet chr, integration_locus, strand).
#'
#' @param ... Represents a single row of a tibble obtained via `obtain_nested`.
#' One row contains 4 variables: chr, integration_locus, strand and data, where
#' data is a nested table that contains all rows that share that integration
#' (collisions).
#' @param date_col The date column to consider
#' @param reads_ratio The value of the ratio between sequence count values to
#' check
#'
#' @importFrom tibble as_tibble tibble
#' @importFrom purrr flatten
#' @importFrom tidyr unnest
#' @importFrom rlang .data env_bind
#' @return A list with:
#' * $data: an updated ISADataFrame with processed collisions or NULL if no
#' criteria was sufficient
#' * $reassigned: 1 if the integration was successfully reassigned, 0 otherwise
#' * $removed: 1 if the integration is removed entirely because no criteria was
#' met, 0 otherwise
.four_step_check <- function(..., date_col, reads_ratio) {
    l <- list(...)
    current <- tibble::as_tibble(l[c("chr", "integration_locus", "strand")])
    current_data <- tibble::as_tibble(purrr::flatten(l["data"]))
    # Try to discriminate by date
    result <- .discriminate_by_date(current_data, date_col)
    if (result$check == TRUE) {
        current_data <- result$data
        res <- tibble::tibble(
            chr = current$chr,
            integration_locus = current$integration_locus,
            strand = current$strand,
            data = list(current_data)
        )
        res <- res %>% tidyr::unnest(.data$data)
        suppressMessages({
            res <- ISADataFrame(res)
        })
        return(list(data = res, reassigned = 1, removed = 0))
    }
    current_data <- result$data
    # If first check fails try to discriminate by replicate
    result <- .discriminate_by_replicate(current_data)
    if (result$check == TRUE) {
        current_data <- result$data
        res <- tibble::tibble(
            chr = current$chr,
            integration_locus = current$integration_locus,
            strand = current$strand,
            data = list(current_data)
        )
        res <- res %>% tidyr::unnest(.data$data)
        suppressMessages({
            res <- ISADataFrame(res)
        })
        return(list(data = res, reassigned = 1, removed = 0))
    }
    current_data <- result$data
    # If second check fails try to discriminate by seqCount
    result <- .discriminate_by_seqCount(current_data, reads_ratio)
    if (result$check == TRUE) {
        current_data <- result$data
        res <- tibble::tibble(
            chr = current$chr,
            integration_locus = current$integration_locus,
            strand = current$strand,
            data = list(current_data)
        )
        res <- res %>% tidyr::unnest(.data$data)
        suppressMessages({
            res <- ISADataFrame(res)
        })
        return(list(data = res, reassigned = 1, removed = 0))
    }
    # If all check fails remove the integration from all subjects
    return(list(data = NULL, reassigned = 0, removed = 1))
}

#' Internal function to process collisions on multiple integrations.
#'
#' @param x The nested collision data frame (or chunks for parallel execution)
#' @param date_col The date column to consider
#' @param reads_ratio The value of the ratio between sequence count values to
#' check
#'
#' @importFrom purrr pmap reduce
#' @importFrom dplyr bind_rows
#' @return A list containing the updated collisions, a numeric value
#' representing the number of integrations removed and a numeric value
#' representing the number of integrations reassigned
.coll_mapping <- function(x, date_col, reads_ratio) {
    result <- purrr::pmap(x,
        .f = .four_step_check,
        date_col = date_col,
        reads_ratio = reads_ratio
    )
    proc_collisions <- purrr::map(result, function(x) {
        x$data
    })
    proc_collisions <- purrr::reduce(proc_collisions, dplyr::bind_rows)
    removed_tot <- purrr::map(result, function(x) {
        x$removed
    })
    removed_tot <- purrr::reduce(removed_tot, sum)
    reassigned_tot <- purrr::map(result, function(x) {
        x$reassigned
    })
    reassigned_tot <- purrr::reduce(reassigned_tot, sum)
    list(
        coll = proc_collisions,
        removed = removed_tot,
        reassigned = reassigned_tot
    )
}

#' Internal function to process collisions on multiple integrations,
#' parallelized.
#'
#' @param collisions The collision data frame
#' @param date_col The date column to consider
#' @param reads_ratio The value of the ratio between sequence count values to
#' check
#'
#' @import BiocParallel
#' @importFrom purrr map reduce
#' @importFrom tibble as_tibble_col add_column
#' @importFrom dplyr group_by group_split
#' @return A list containing the updated collisions, a numeric value
#' representing the number of integrations removed and a numeric value
#' representing the number of integrations reassigned
.process_collisions <- function(collisions, date_col, reads_ratio) {
    # Obtain nested version of collisions
    nested <- .obtain_nested(collisions)
    # Register backend according to platform
    if (.Platform$OS.type == "windows") {
        p <- BiocParallel::SnowParam(stop.on.error = FALSE)
    } else {
        p <- BiocParallel::MulticoreParam(stop.on.error = FALSE)
    }
    # Split the data frame in chunks according to number of workers
    workers <- BiocParallel::bpworkers(p)
    exceeding <- nrow(nested) %% workers
    ampl <- trunc(nrow(nested) / workers)
    chunks <- unlist(lapply(seq_len(workers), FUN = rep_len, length.out = ampl))
    if (exceeding > 0) {
        chunks <- c(chunks, rep_len(
            x = tail(workers, n = 1),
            length.out = exceeding
        ))
    }
    chunks <- tibble::as_tibble_col(chunks, column_name = "chunk")
    nested <- nested %>% tibble::add_column(chunks)
    split_data <- nested %>%
        dplyr::group_by(.data$chunk) %>%
        dplyr::group_split(.keep = FALSE)
    # For each chunk process collisions in parallel
    suppressMessages(suppressWarnings({
        ### result is a list with n elements, each element is a list of 3,
        ### containing processed collisions, number of removed collisions and
        ### number of reassigned collisions
        result <- BiocParallel::bptry(
            BiocParallel::bplapply(split_data,
                FUN = .coll_mapping,
                date_col = date_col,
                reads_ratio = reads_ratio,
                BPPARAM = p
            )
        )
    }))
    BiocParallel::bpstop(p)
    # For each element of result extract and bind the 3 components
    processed_collisions_df <- purrr::map(result, function(x) {
        x$coll
    })
    processed_collisions_df <- purrr::reduce(
        processed_collisions_df,
        dplyr::bind_rows
    )
    removed_total <- purrr::map(result, function(x) {
        x$removed
    })
    removed_total <- purrr::reduce(removed_total, sum)
    reassigned_total <- purrr::map(result, function(x) {
        x$reassigned
    })
    reassigned_total <- purrr::reduce(reassigned_total, sum)
    # Return the summary
    list(
        coll = processed_collisions_df, removed = removed_total,
        reassigned = reassigned_total
    )
}

#' Internal function to obtain a summary table after collision processing.
#'
#' @param before The joined table before collision removing
#' (obtained via `.join_matrix_af`)
#' @param after The final matrix obtained after collision processing
#' @param association_file The association file
#' @import dplyr
#' @importFrom rlang .data
#'
#' @return A tibble with a summary containing for each SubjectID the number of
#' integrations found before and after, the sum of the value of the sequence
#' count for each subject before and after and the corresponding deltas.
.summary_table <- function(before, after, association_file) {
    after_matr <- after %>%
        dplyr::left_join(association_file,
            by = c("CompleteAmplificationID")
        ) %>%
        dplyr::select(dplyr::all_of(colnames(after)), .data$SubjectID)

    n_int_after <- after_matr %>%
        dplyr::group_by(.data$SubjectID) %>%
        dplyr::tally(name = "count_integrations_after")
    sumReads_after <- after_matr %>%
        dplyr::group_by(.data$SubjectID) %>%
        dplyr::summarise(
            sumSeqReads_after = sum(.data$Value),
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
            sumSeqReads_before = sum(.data$Value),
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
