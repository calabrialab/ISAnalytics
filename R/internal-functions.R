#------------------------------------------------------------------------------#
# Internal functions
#------------------------------------------------------------------------------#
## All functions in this file are NOT exported, to be used internally only.

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

# Internal helper function for checking `Value` column presence in x.
#
# Checks if the column `Value` is present in the data frame and also
# checks if the column is numeric or integer.
# @param x A data.frame object (or any extending class)
# @keywords internal
#
# @return FALSE if not found or contains non-numeric data, TRUE otherwise
.check_value_col <- function(x) {
    stopifnot(is.data.frame(x))
    present <- if ("Value" %in% colnames(x)) {
        TRUE
    } else {
        FALSE
    }
    if (present == TRUE) {
        return(is.numeric(x$Value) | is.integer(x$Value))
    } else {
        return(FALSE)
    }
}

# Internal helper function for checking `CompleteAmplifcationID`
# column presence in x.
#
# Checks if the column `CompleteAmplifcationID` is present in the data frame.
#
# @param x A data.frame object (or any extending class)
# @keywords internal
#
# @return FALSE if not found, TRUE otherwise
.check_complAmpID <- function(x) {
    stopifnot(is.data.frame(x))
    if ("CompleteAmplificationID" %in% colnames(x)) {
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
.find_exp_cols <- function(x) {
    stopifnot(is.data.frame(x))
    default_cols <- c(
        mandatory_IS_vars(), annotation_IS_vars(),
        "CompleteAmplificationID"
    )
    remaining <- colnames(x)[!colnames(x) %in% default_cols]
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

# Internal function to convert a messy matrix to a tidy data frame
#
# @description Uses the suite of functions provided by the
# tidyverse to produce a more dense and ordered structure.
# This function is not exported and should be called in other importing
# functions.
#
# @param df Messy tibble to convert to tidy
# @keywords internal
#
# @return a tidy tibble
#' @importFrom rlang .data
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr arrange all_of filter
#' @importFrom forcats fct_inseq as_factor
.messy_to_tidy <- function(df) {
    exp_cols <- which(!(colnames(df) %in% c(
        mandatory_IS_vars(),
        annotation_IS_vars()
    )))
    isadf_tidy <- df %>%
        tidyr::pivot_longer(
            cols = dplyr::all_of(exp_cols),
            names_to = "CompleteAmplificationID",
            values_to = "Value",
            values_drop_na = TRUE
        ) %>%
        dplyr::filter(.data$Value > 0) %>%
        dplyr::arrange(forcats::fct_inseq(forcats::as_factor(.data$chr)))
    isadf_tidy
}

# Internal function to auto-detect the type of IS based on the headers.
#
# @param df the data frame to inspect
# @keywords internal
#
# @return one value among:
# * "OLD" : for old-style matrices that had only one column holding
# all genomic coordinates
# * "NEW" :  for the classic Vispa2 annotated/not annotated matrices
# * "MALFORMED" : in any other case
.auto_detect_type <- function(df) {
    if ("IS_genomicID" %in% colnames(df) &
        all(!mandatory_IS_vars() %in% colnames(df))) {
        return("OLD")
    }
    if (all(mandatory_IS_vars() %in% colnames(df))) {
        return("NEW")
    }
    return("MALFORMED")
}

#---- USED IN : import_association_file ----

# Checks if the association file has the right format (standard headers).
#
# @param df The imported association file
# @keywords internal
#
# @return TRUE if the check passes, FALSE otherwise
.check_af_correctness <- function(df) {
    if (all(association_file_columns() %in% colnames(df))) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}


# Reads association file and checks if it's correct or not.
#
# @param path The path to the association file on disk
# @param padding The padding for TimePoint field
# @param date_format The date format of date columns
# @keywords internal
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate across contains
#' @importFrom rlang .data
#' @importFrom stringr str_pad
#' @import lubridate
#
# @return A tibble containing the association file.
.read_and_correctness_af <- function(path, padding, date_format) {
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
    as_file <- as_file %>%
        dplyr::mutate(TimePoint = stringr::str_pad(
            as.character(.data$TimePoint),
            padding,
            side = "left",
            pad = "0"
        )) %>%
        dplyr::mutate(dplyr::across(dplyr::contains("Date"), ~ do.call(
            getExportedValue("lubridate", date_format),
            list(.x),
            quote = TRUE
        )))
    as_file
}

# Internal function to check alignment between association file and file
# system starting from the root. The alignment is checked at folder level.
#
# @param df The imported association file (data.frame or tibble)
# @param root_folder Path to the root folder
# @keywords internal
#' @importFrom dplyr select distinct mutate bind_rows
#' @importFrom fs dir_ls
#' @importFrom purrr pmap is_empty reduce map_dbl
#' @importFrom stringr str_replace_all str_extract_all
#' @importFrom tibble tibble
#
# @return A data frame containing, for each ProjectID and
# concatenatePoolIDSeqRun the
# corresponding path on disk if found, NA otherwise.
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
                replacement = "\\\\/"
            )
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

# Updates the association file after the alignment check.
#
# The function checks if there are missing folders and updates the association
# file by adding a column `Path` where the absolute path on disk for the
# project and pool is found, if no path is found NA is inserted instead.
#
# @param as_file The tibble representing the read association_file
# @param checks The tibble representing the results
# of `.check_file_system_alignment`
# @param root The root folder
# @keywords internal
#
# @return An updated association file with absolute paths
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

# Helper function to be used when importing matrices in parallel.
#
# @param association_file Either the path to the association file or the tibble
# representing the imported association file (done via
# `import_association_file`)
# @param root If `association_file` is the path to the file, root is a single
# string holding the path to the root folder, otherwise root is `NULL`
# @param padding The padding for TimePoint field
# @param format The date format of date columns
# @keywords internal
#
# @return A list of two elements: the first element is an updated version of
# the association file with NAs removed, the second element is a widget showing
# results of alignment checks.
.manage_association_file <- function(association_file, root, padding, format) {
    # Manage association file
    if (is.character(association_file)) {
        # If it's a path to file import the association file
        association_file <- .read_and_correctness_af(
            association_file,
            padding, format
        )
        checks <- .check_file_system_alignment(association_file, root)
        association_file <- .update_af_after_alignment(
            association_file,
            checks, root
        )
        association_file <- association_file %>% dplyr::filter(!is.na(
            .data$Path
        ))
        res <- list(association_file, checks)
        return(res)
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
        connection <- getOption("ISAnalytics.connection")
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
#
# @param association_file The tibble representing the imported
# association file
#' @importFrom dplyr select distinct group_by bind_rows inner_join
#' @importFrom tibble tibble
#' @importFrom tidyr nest
#' @importFrom purrr map pmap reduce
#' @importFrom rlang .data
# @keywords internal
#
# @return A modified version of the association file where only selected
# pools for each project are present
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
#' @importFrom purrr map reduce map_dbl
#' @importFrom stringr str_detect
#' @importFrom dplyr select distinct bind_rows mutate
#' @importFrom tidyr nest
#
# @return A tibble containing all found files, including duplicates and missing
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

# Internal function for parallel import of a single quantification
# type files.
#
# @param q_type The quantification type (single string)
# @param files Files_found table were absolute paths of chosen files
# are stored
# @param workers Number of parallel workers
# @keywords internal
#' @importFrom dplyr filter mutate bind_rows distinct
#' @importFrom BiocParallel SnowParam MulticoreParam bptry bplapply bpstop bpok
#' @importFrom purrr is_empty reduce
#
# @return A single tibble with all data from matrices of same quantification
# type in tidy format
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
        tibble::as_tibble(purrr::flatten(x[1]))
    })

    list(imported_matrices, summary_files)
}

#---- USED IN : import_parallel_Vispa2Matrices_auto ----

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
# @keywords internal
#' @importFrom tibble as_tibble_col
#' @importFrom stringr str_detect
#' @importFrom purrr map reduce
#' @importFrom dplyr bind_cols
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
# @keywords internal
#' @import dplyr
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

# Looks up matrices to import given the association file and the
# root of the file system.
#
# @inheritParams .lookup_matrices
# @param patterns A character vector of patterns to be matched
# @param matching_opt A single character representing the matching option (one
# of "ANY", "ALL" or "OPTIONAL")
# @keywords internal
#' @importFrom purrr pmap flatten map
#' @importFrom tibble as_tibble
#' @importFrom stringr str_split
#' @importFrom dplyr mutate select
#' @importFrom utils tail
#
# @return A tibble containing all found files, including duplicates and missing
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

#### ---- Internals for collision removal ----####

#---- USED IN : remove_collisions ----

# Checks if association file contains more information than the matrix.
#
# Used to notify the user that wants to know if for the projects contained in
# the examined matrix there are additional CompleteAmplificationIDs contained
# in the association file that weren't included in the integration matrix (for
# example failed runs).
#
# @param association_file The imported association file
# @param df The sequence count matrix to examine
# @keywords internal
#' @import dplyr
#' @importFrom purrr map reduce
#' @importFrom rlang .data
#
# @return A tibble containing ProjectID, PoolID and CompleteAmplificationID
# only for additional info and only for projects already present in df
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

# Produces a joined tibble between the sequence count matrix and the
# association file
#
# @param seq_count_df The sequence count tibble
# @param association_file The association file tibble
# @param date_col The date column chosen
# @keywords internal
#' @import dplyr
#' @importFrom rlang .data
#
# @return A tibble
.join_matrix_af <- function(seq_count_df, association_file, date_col) {
    joined <- seq_count_df %>%
        dplyr::left_join(association_file, by = "CompleteAmplificationID") %>%
        dplyr::select(
            dplyr::all_of(colnames(seq_count_df)), all_of(date_col),
            .data$ProjectID, .data$PoolID, .data$SubjectID,
            .data$ReplicateNumber
        )
    joined
}

# Identifies independent samples and separates the joined_df in
# collisions and non-collisions
#
# @param joined_df The joined tibble obtained via `.join_matrix_af`
#' @import dplyr
#' @importFrom rlang .data
# @keywords internal
#
# @return A named list containing the splitted joined_df for collisions and
# non-collisions
.identify_independent_samples <- function(joined_df) {
    temp <- joined_df %>%
        dplyr::select(
            dplyr::all_of(mandatory_IS_vars()),
            .data$ProjectID, .data$SubjectID
        ) %>%
        dplyr::group_by(dplyr::across(mandatory_IS_vars())) %>%
        dplyr::distinct(.data$ProjectID, .data$SubjectID, .keep_all = TRUE) %>%
        dplyr::summarise(n = dplyr::n(), .groups = "drop_last") %>%
        dplyr::ungroup() %>%
        dplyr::filter(.data$n > 1) %>%
        dplyr::select(-c(.data$n))

    non_collisions <- joined_df %>%
        dplyr::anti_join(temp, by = mandatory_IS_vars())
    collisions <- joined_df %>%
        dplyr::right_join(temp, by = mandatory_IS_vars())
    list(collisions = collisions, non_collisions = non_collisions)
}

# Returns a polished version of collisions table for analysis.
#
# Polishing consists in nesting all information relative to a single
# integration, in other words, obtaining a linked data frame for each
# integration.
#
# @param collisions The collisions table obtained via
# `.identify_independent_samples`
# @keywords internal
#' @importFrom dplyr group_by across
#' @importFrom rlang .data
#' @importFrom tidyr nest
#
# @return A nested tibble
.obtain_nested <- function(collisions) {
    collisions %>%
        dplyr::group_by(dplyr::across(mandatory_IS_vars())) %>%
        tidyr::nest()
}


# Internal for date discrimination in collision removal.
#
# It's the first of 4 steps in the algorithm for collision removal: it tries to
# find a single sample who has an associated date which is earlier than any
# other. If comparison is not possible the analysis fails and returns the
# orginal input data.
#
# @details NOTE: this function is meant to be used inside a mapping function
# such as `purrr::pmap`. The function only works on data regarding a SINGLE
# integration (triplet chr, integration_locus, strand).
#
# @param nest The nested table associated with a single integration
# @param date_col The name of the date column chosen for collision removal,
# one in \code{date_columns_coll()}
#' @importFrom rlang expr eval_tidy .data
#' @importFrom dplyr filter arrange
# @keywords internal
#
# @return A named list with:
# * $data: a tibble, containing the data (unmodified or modified)
# * $check: a logical value indicating whether the analysis was successful or
# not (and therefore there is the need to perform the next step)
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

# Internal for replicate discrimination in collision removal.
#
# It's the second of 4 steps in the algorithm for collision removal:
# grouping by independent sample (ProjectID, SubjectID) it counts the number of
# rows (replicates) found for each group, orders them from biggest to smallest
# and, if a single group has more rows than any other group the integration is
# assigned to that sample, otherwise the analysis fails and returns the
# original input to be submitted to the next step.
#
# @details NOTE: this function is meant to be used inside a mapping function
# such as `purrr::pmap`. The function only works on data regarding a SINGLE
# integration (triplet chr, integration_locus, strand).
#
# @param nest The nested table associated with a single integration
#' @import dplyr
#' @importFrom rlang .data
# @keywords internal
#
# @return A named list with:
# * $data: a tibble, containing the data (unmodified or modified)
# * $check: a logical value indicating whether the analysis was successful or
# not (and therefore there is the need to perform the next step)
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

# Internal for sequence count discrimination in collision removal.
#
# It's the third of 4 steps in the algorithm for collision removal:
# grouping by independent sample (ProjectID, SubjectID), it sums the value of
# the sequence count for each group, sorts the groups from highest to lowest
# cumulative value and then checks the ratio between the first element and the
# second: if the ratio is > `reads_ratio` the integration is assigned to
# the first group, otherwise the analysis fails and returns the original input.
#
# @details NOTE: this function is meant to be used inside a mapping function
# such as `purrr::pmap`. The function only works on data regarding a SINGLE
# integration (triplet chr, integration_locus, strand).
#
# @param nest The nested table associated with a single integration
# @param reads_ratio The value of the ratio between sequence count values to
# check
# @param seqCount_col The name of the sequence count column (support
# for multi quantification matrix)
#' @import dplyr
#' @importFrom rlang .data eval_tidy expr
# @keywords internal
#
# @return A named list with:
# * $data: a tibble, containing the data (unmodified or modified)
# * $check: a logical value indicating whether the analysis was successful or
# not (and therefore there is the need to perform the next step)
.discriminate_by_seqCount <- function(nest, reads_ratio, seqCount_col) {
    temp <- nest %>%
        dplyr::group_by(.data$ProjectID, .data$SubjectID) %>%
        dplyr::summarise(sum_seqCount = sum(
            rlang::eval_tidy(rlang::expr(`$`(.data, !!seqCount_col)))
        ), .groups = "keep") %>%
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

# Internal function that performs four-step-check of collisions for
# a single integration.
#
# @details NOTE: this function is meant to be used inside a mapping function
# such as `purrr::pmap`. The function only works on data regarding a SINGLE
# integration (triplet chr, integration_locus, strand).
#
# @param ... Represents a single row of a tibble obtained via
# `obtain_nested`. One row contains 4 variables: chr, integration_locus,
# strand and data, where data is a nested table that contains all rows
# that share that integration (collisions).
# @param date_col The date column to consider
# @param reads_ratio The value of the ratio between sequence count values to
# check
# @param seqCount_col The name of the sequence count column (support
# for multi quantification matrix)
# @keywords internal
#
#' @importFrom tibble as_tibble tibble
#' @importFrom purrr flatten
#' @importFrom tidyr unnest
#' @importFrom rlang .data env_bind
# @return A list with:
# * $data: an updated tibble with processed collisions or NULL if no
# criteria was sufficient
# * $reassigned: 1 if the integration was successfully reassigned, 0 otherwise
# * $removed: 1 if the integration is removed entirely because no criteria was
# met, 0 otherwise
.four_step_check <- function(..., date_col, reads_ratio, seqCount_col) {
    l <- list(...)
    current <- tibble::as_tibble(l[mandatory_IS_vars()])
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
        return(list(data = res, reassigned = 1, removed = 0))
    }
    current_data <- result$data
    # If second check fails try to discriminate by seqCount
    result <- .discriminate_by_seqCount(current_data, reads_ratio, seqCount_col)
    if (result$check == TRUE) {
        current_data <- result$data
        res <- tibble::tibble(
            chr = current$chr,
            integration_locus = current$integration_locus,
            strand = current$strand,
            data = list(current_data)
        )
        res <- res %>% tidyr::unnest(.data$data)
        return(list(data = res, reassigned = 1, removed = 0))
    }
    # If all check fails remove the integration from all subjects
    return(list(data = NULL, reassigned = 0, removed = 1))
}

# Internal function to process collisions on multiple integrations.
#
# @param x The nested collision data frame (or chunks for parallel execution)
# @param date_col The date column to consider
# @param reads_ratio The value of the ratio between sequence count values to
# check
# @param seqCount_col The name of the sequence count column (support
# for multi quantification matrix)
# @keywords internal
#' @importFrom purrr pmap reduce
#' @importFrom dplyr bind_rows
# @return A list containing the updated collisions, a numeric value
# representing the number of integrations removed and a numeric value
# representing the number of integrations reassigned
.coll_mapping <- function(x, date_col, reads_ratio, seqCount_col) {
    result <- purrr::pmap(x,
        .f = .four_step_check,
        date_col = date_col,
        reads_ratio = reads_ratio,
        seqCount_col = seqCount_col
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

# Internal function to process collisions on multiple integrations,
# parallelized.
#
# @param collisions The collision data frame
# @param date_col The date column to consider
# @param reads_ratio The value of the ratio between sequence count values to
# check
# @param seqCount_col The name of the sequence count column (support
# for multi quantification matrix)
# @keywords internal
#
#' @import BiocParallel
#' @importFrom purrr map reduce
#' @importFrom tibble as_tibble_col add_column
#' @importFrom dplyr group_by group_split
# @return A list containing the updated collisions, a numeric value
# representing the number of integrations removed and a numeric value
# representing the number of integrations reassigned
.process_collisions <- function(collisions, date_col, reads_ratio,
    seqCount_col) {
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
                seqCount_col = seqCount_col,
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

# Internal function to obtain a summary table after collision processing.
#
# @param before The joined table before collision removing
# (obtained via `.join_matrix_af`)
# @param after The final matrix obtained after collision processing
# @param association_file The association file
#' @import dplyr
#' @importFrom rlang .data
# @keywords internal
#
# @return A tibble with a summary containing for each SubjectID the number of
# integrations found before and after, the sum of the value of the sequence
# count for each subject before and after and the corresponding deltas.
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

#### ---- Internals for aggregate functions ----####

#---- USED IN : aggregate_metadata ----

# Minimal association_file variable set.
#
# Contains the names of the columns of the association file that are a minimum
# requirement to perform aggregation.
# @keywords internal
#
# @return A character vector
.min_var_set <- function() {
    c(
        "FusionPrimerPCRDate", "LinearPCRDate", "VCN", "DNAngUsed", "Kapa",
        "ulForPool", "Path"
    )
}

# Minimal stats column set.
#
# Contains the name of the columns that are a minimum requirement for
# aggregation.
# @keywords internal
#
# @return A character vector
.stats_columns_min <- function() {
    c(
        "POOL", "TAG", "BARCODE_MUX", "TRIMMING_FINAL_LTRLC",
        "LV_MAPPED", "BWA_MAPPED_OVERALL", "ISS_MAPPED_PP"
    )
}

# Checks if the stats file contains the minimal set of columns.
#
# @param x The stats df
# @keywords internal
#
# @return TRUE or FALSE
.check_stats <- function(x) {
    if (all(.stats_columns_min() %in% colnames(x))) {
        TRUE
    } else {
        FALSE
    }
}

# Finds automatically the path on disk to each stats file.
#
# @param association_file The association file
#' @import dplyr
#' @importFrom purrr pmap_dfr pmap_df
#' @importFrom tibble as_tibble tibble add_column
#' @importFrom stringr str_extract str_extract_all str_detect
#' @importFrom fs dir_exists dir_ls
#' @importFrom tidyr unnest unite
#' @importFrom rlang .data
# @keywords internal
# @return A tibble with ProjectID and the absolute path on disk to each file
# if found
.stats_report <- function(association_file) {
    # Obtain unique projectID and path
    temp <- association_file %>%
        dplyr::select(.data$ProjectID, .data$Path) %>%
        dplyr::distinct()
    # If paths are all NA return
    if (all(is.na(temp$Path))) {
        return(NULL)
    }
    pattern <- "^stats\\."
    file_reg <- "(?<=iss\\/).*"
    # Obtain paths to iss folder
    stats_paths <- purrr::pmap_dfr(temp, function(...) {
        current <- tibble::tibble(...)
        if (is.na(current$Path)) {
            l <- list(ProjectID = current$ProjectID, stats_path = NA_character_)
            tibble::as_tibble(l)
        } else {
            pj_path <-
                stringr::str_extract(current$Path,
                    pattern = paste0(
                        ".*",
                        current$ProjectID,
                        "(?=\\/quantification)"
                    )
                )
            pj_path <- paste0(pj_path, "/iss")
            l <- list(ProjectID = current$ProjectID, stats_path = pj_path)
            tibble::as_tibble(l)
        }
    })
    stats_paths <- stats_paths %>% dplyr::distinct()
    # Find all stats files in iss folders
    stats_paths <- purrr::pmap_df(stats_paths, function(...) {
        cur <- tibble::tibble(...)
        # Check if folder exists
        # Set to NA the iss folders not found
        if (!fs::dir_exists(cur$stats_path)) {
            cur$stats_path <- NA_character_
            cur %>% tibble::add_column(stats_files = list(NA_character_))
        } else {
            files_in_iss <- unlist(fs::dir_ls(cur$stats_path))
            files_in_iss <- unlist(stringr::str_extract_all(
                files_in_iss, file_reg
            ))
            which_match <- unlist(stringr::str_detect(files_in_iss, pattern))
            files_in_iss <- files_in_iss[which_match]
            cur %>% tibble::add_column(stats_files = list(files_in_iss))
        }
    })
    stats_paths <- stats_paths %>%
        tidyr::unnest(.data$stats_files) %>%
        tidyr::unite(
            col = "files", .data$stats_path, .data$stats_files,
            sep = "/", na.rm = TRUE
        )
}

# Imports all found Vispa2 stats files.
#
# @param association_file The association file
#' @import BiocParallel
#' @importFrom tibble as_tibble
#' @importFrom purrr map2_lgl reduce is_empty
#' @importFrom dplyr mutate bind_rows distinct
# @keywords internal
#
# @return A list with the imported stats and a report of imported files. If no
# files were imported returns NULL instead
.import_stats_iss <- function(association_file) {
    # Obtain paths
    stats_paths <- .stats_report(association_file)
    if (is.null(stats_paths)) {
        return(NULL)
    }
    # Setup parallel workers and import
    # Register backend according to platform
    if (.Platform$OS.type == "windows") {
        p <- BiocParallel::SnowParam(stop.on.error = FALSE)
    } else {
        p <- BiocParallel::MulticoreParam(
            stop.on.error = FALSE
        )
    }
    FUN <- function(x) {
        stats <- read.csv(
            file = x, sep = "\t", stringsAsFactors = FALSE,
            check.names = FALSE
        )
        stats <- tibble::as_tibble(stats)
        ok <- .check_stats(stats)
        if (ok == TRUE) {
            return(stats)
        } else {
            return(NULL)
        }
    }
    suppressMessages(suppressWarnings({
        stats_dfs <- BiocParallel::bptry(
            BiocParallel::bplapply(stats_paths$files, FUN, BPPARAM = p)
        )
    }))
    BiocParallel::bpstop(p)
    correct <- BiocParallel::bpok(stats_dfs)
    correct <- purrr::map2_lgl(stats_dfs, correct, function(x, y) {
        if (is.null(x)) {
            FALSE
        } else {
            y
        }
    })
    stats_paths <- stats_paths %>% dplyr::mutate(Imported = correct)
    stats_dfs <- stats_dfs[correct]
    # Bind rows in single tibble for all files
    if (purrr::is_empty(stats_dfs)) {
        return(NULL)
    }
    stats_dfs <- purrr::reduce(stats_dfs, function(x, y) {
        x %>%
            dplyr::bind_rows(y) %>%
            dplyr::distinct()
    })
    list(stats_dfs, stats_paths)
}

# Performs eventual join with stats and aggregation.
#
# @param association_file The imported association file
# @param stats The imported stats df
# @param grouping_keys The names of the columns to group by
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom stringr str_replace_all
#' @importFrom tidyr unite
# @keywords internal
#
# @return A tibble
.join_and_aggregate <- function(association_file, stats, grouping_keys) {
    # If stats not null join with af
    if (!is.null(stats)) {
        stats <- stats %>% dplyr::mutate(TAG = stringr::str_replace_all(
            .data$TAG,
            pattern = "\\.", replacement = ""
        ))
        association_file <- association_file %>%
            dplyr::left_join(stats,
                by = c(
                    "concatenatePoolIDSeqRun" = "POOL",
                    "TagSequence" = "TAG"
                )
            )
        aggregated <- association_file %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(grouping_keys))) %>%
            dplyr::summarise(
                FusionPrimerPCRDate = suppressWarnings({
                    min(.data$FusionPrimerPCRDate,
                        na.rm = TRUE
                    )
                }),
                LinearPCRDate = suppressWarnings({
                    min(.data$LinearPCRDate, na.rm = TRUE)
                }),
                VCN = suppressWarnings({
                    mean(.data$VCN, na.rm = TRUE)
                }),
                Avg_DNAngUsed = suppressWarnings({
                    mean(.data$DNAngUsed, na.rm = TRUE)
                }),
                Kapa = suppressWarnings({
                    mean(.data$Kapa, na.rm = TRUE)
                }),
                DNAngUsed = suppressWarnings({
                    sum(.data$DNAngUsed, na.rm = TRUE)
                }),
                ulForPool = suppressWarnings({
                    sum(.data$ulForPool, na.rm = TRUE)
                }),
                BARCODE_MUX = suppressWarnings({
                    sum(.data$BARCODE_MUX, na.rm = TRUE)
                }),
                TRIMMING_FINAL_LTRLC = suppressWarnings({
                    sum(.data$TRIMMING_FINAL_LTRLC, na.rm = TRUE)
                }),
                LV_MAPPED = suppressWarnings({
                    sum(.data$LV_MAPPED, na.rm = TRUE)
                }),
                BWA_MAPPED_OVERALL = suppressWarnings({
                    sum(.data$BWA_MAPPED_OVERALL, na.rm = TRUE)
                }),
                ISS_MAPPED_PP = suppressWarnings({
                    sum(.data$ISS_MAPPED_PP, na.rm = TRUE)
                })
            ) %>%
            dplyr::ungroup() %>%
            tidyr::unite(
                col = "AggregateMeta", dplyr::all_of(grouping_keys),
                sep = "_", remove = FALSE
            )
        return(aggregated)
    }
    # Aggregate
    aggregated <- association_file %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(grouping_keys))) %>%
        dplyr::summarise(
            FusionPrimerPCRDate = suppressWarnings({
                min(.data$FusionPrimerPCRDate,
                    na.rm = TRUE
                )
            }),
            LinearPCRDate = suppressWarnings({
                min(.data$LinearPCRDate, na.rm = TRUE)
            }),
            VCN = suppressWarnings({
                mean(.data$VCN, na.rm = TRUE)
            }),
            Avg_DNAngUsed = suppressWarnings({
                mean(.data$DNAngUsed, na.rm = TRUE)
            }),
            Kapa = suppressWarnings({
                mean(.data$Kapa, na.rm = TRUE)
            }),
            DNAngUsed = suppressWarnings({
                sum(.data$DNAngUsed, na.rm = TRUE)
            }),
            ulForPool = suppressWarnings({
                sum(.data$ulForPool, na.rm = TRUE)
            })
        ) %>%
        dplyr::ungroup() %>%
        tidyr::unite(
            col = "AggregateMeta", dplyr::all_of(grouping_keys),
            sep = "_", remove = FALSE
        )
}

#---- USED IN : aggregate_values_by_key ----

# Internal function that performs aggregation on values with a lambda
# operation.
#
# @param x The list of matrices to aggregate. If a single matrix has to be
# supplied it must be enclosed in a list. For example `x = list(matrix)`.
# @param af The association file
# @param key A string or character vector to use as key
# @param lambda The aggregating operation to apply to values. Must take as
# input a numeric/integer vector and return a single value
# @param group The additional variables to add to grouping
# @param args Additional arguments passed on to lambda (named list)
# @param namespace The namespace from which the lambda is exported
# @param envir The environment in which symbols must be evaluated
#' @import dplyr
#' @importFrom rlang .data
# @keywords internal
#
# @return A list of tibbles with aggregated values
.aggregate_lambda <- function(x, af, key, value_cols, lambda, group) {
    # Vectorize
    aggregated_matrices <- purrr::map(x, function(y) {
        cols <- c(colnames(y), key)
        agg <- y %>%
            dplyr::left_join(af, by = "CompleteAmplificationID") %>%
            dplyr::select(dplyr::all_of(cols)) %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(c(group, key)))) %>%
            dplyr::summarise(dplyr::across(
                .cols = dplyr::all_of(value_cols),
                .fns = lambda
            ), .groups = "drop")
        agg
    })
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
#
# @param x An integration matrix subset (see description)
# @param threshold The numeric value representing an absolute number
# of bases for which two integrations are considered distinct
# @param keep_criteria The string with selecting criteria
# @param annotated Is `x` annotated? Logical value
# @param num_cols A character vector with the names of the numeric columns
# @param max_val_col The column to consider if criteria is max_value
# @param produce_map Produce recalibration map?
#' @importFrom tibble tibble tibble_row add_row
#' @importFrom dplyr arrange bind_rows distinct
#' @importFrom purrr is_empty map_dfr
#' @importFrom rlang expr eval_tidy
# @keywords internal
#
# @return A named list with recalibrated matrix and recalibration map.
.sliding_window <- function(x,
    threshold,
    keep_criteria,
    annotated,
    num_cols,
    max_val_col,
    produce_map) {
    ## Order by integration_locus
    x <- x %>% dplyr::arrange(.data$integration_locus)
    map_recalibr <- if (produce_map == TRUE) {
        tibble::tibble(
            chr_before = character(0),
            integration_locus_before = integer(0),
            strand_before = character(0),
            chr_after = character(0),
            integration_locus_after = integer(0),
            strand_after = character(0)
        )
    } else {
        NULL
    }
    index <- 1
    repeat {
        # If index is the last row in the data frame return
        if (index == nrow(x)) {
            if (produce_map == TRUE) {
                map_recalibr <- tibble::add_row(map_recalibr,
                    chr_before = x$chr[index],
                    integration_locus_before =
                        x$integration_locus[index],
                    strand_before = x$strand[index],
                    chr_after = x$chr[index],
                    integration_locus_after =
                        x$integration_locus[index],
                    strand_after = x$strand[index]
                )
            }
            if (!is.null(map_recalibr)) {
                map_recalibr <- dplyr::distinct(
                    map_recalibr
                )
            }
            return(list(recalibrated_matrix = x, map = map_recalibr))
        }
        ## Compute interval for row
        interval <- x[index, ]$integration_locus + 1 + threshold
        ## Look ahead for every integration that falls in the interval
        near <- numeric()
        k <- index
        repeat {
            if (k == nrow(x)) {
                break
            }
            k <- k + 1
            if (x[k, ]$integration_locus < interval) {
                # Saves the indexes of the rows that are in the interval
                near <- append(near, k)
            } else {
                break
            }
        }
        window <- c(index, near)
        if (!purrr::is_empty(near)) {
            ## Change loci according to criteria
            ######## CRITERIA PROCESSING
            row_to_keep <- index
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
                recalib_rows <- purrr::map_dfr(window, function(cur_row) {
                    tibble::tibble_row(
                        chr_before = x$chr[cur_row],
                        integration_locus_before =
                            x$integration_locus[cur_row],
                        strand_before = x$strand[cur_row],
                        chr_after = x$chr[row_to_keep],
                        integration_locus_after =
                            x$integration_locus[row_to_keep],
                        strand_after = x$strand[row_to_keep]
                    )
                })
                map_recalibr <- dplyr::bind_rows(map_recalibr, recalib_rows)
            }
            # Change loci and strand of near integrations
            x[near, ]$integration_locus <- x[row_to_keep, ]$integration_locus
            x[near, ]$strand <- x[row_to_keep, ]$strand

            if (annotated == TRUE) {
                x[near, ]$GeneName <- x[row_to_keep, ]$GeneName
                x[near, ]$GeneStrand <- x[row_to_keep, ]$GeneStrand
            }
            ## Aggregate same IDs
            starting_rows <- nrow(x)
            repeat {
                t <- x$CompleteAmplificationID[window]
                d <- unique(t[duplicated(t)])
                if (purrr::is_empty(d)) {
                    break
                }
                dupl_indexes <- which(t == d[1])
                values_sum <- colSums(x[window[dupl_indexes], num_cols],
                    na.rm = TRUE
                )
                x[window[dupl_indexes[1]], num_cols] <- as.list(values_sum)
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
                index <- tail(window, n = 1) + 1
            }
        } else {
            if (produce_map == TRUE) {
                map_recalibr <- tibble::add_row(map_recalibr,
                    chr_before = x$chr[index],
                    integration_locus_before =
                        x$integration_locus[index],
                    strand_before = x$strand[index],
                    chr_after = x$chr[index],
                    integration_locus_after =
                        x$integration_locus[index],
                    strand_after = x$strand[index]
                )
            }
            index <- index + 1
        }
    }
    if (!is.null(map_recalibr)) {
        map_recalibr <- dplyr::distinct(
            map_recalibr
        )
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
.generate_filename <- function() {
    time <- as.character(Sys.time())
    time <- stringr::str_replace_all(
        string = time,
        pattern = "-",
        replacement = "_"
    )
    time <- stringr::str_replace_all(
        string = time,
        pattern = " ",
        replacement = "_"
    )
    time <- stringr::str_replace_all(
        string = time,
        pattern = ":",
        replacement = ""
    )
    filename <- paste0("recalibr_map_", time, ".tsv")
    filename
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
    withRestarts(
        {
            if (file_path == ".") {
                filename <- .generate_filename()
                fs::dir_create(fs::path_wd("recalibration_maps"))
                file_path <- fs::path_wd("recalibration_maps", filename)
                readr::write_tsv(map, file_path)
                if (getOption("ISAnalytics.verbose") == TRUE) {
                    message(paste(
                        "Recalibration map saved to: ",
                        file_path
                    ))
                }
            } else {
                file_path <- fs::path(file_path)
                readr::write_tsv(map, file_path)
                if (getOption("ISAnalytics.verbose") == TRUE) {
                    message(paste(
                        "Recalibration map saved to: ",
                        file_path
                    ))
                }
            }
        },
        skip_write = function() {
            message(paste(
                "Could not write recalibration map file.",
                "Skipping."
            ))
        }
    )
}

#### ---- Internals for analysis functions ----####
#---- USED IN : CIS_grubbs ----
### link to article in documentation
.lentiviral_CIS_paper <- function() {
    paste0(
        "https://ashpublications.org/blood/article/117/20/5332/21206/",
        "Lentiviral-vector-common-integration-sites-in"
    )
}

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
                    stop(.missing_user_cols_list_error(set, name),
                        call. = FALSE
                    )
                }
                purrr::walk(set, function(col) {
                    expr <- rlang::expr(`$`(df, !!col))
                    if (!is.numeric(rlang::eval_tidy(expr)) &&
                        !is.integer(rlang::eval_tidy(expr))) {
                        stop(.non_num_user_cols_error(),
                            call. = FALSE
                        )
                    }
                })
            }
        )
    } else {
        if (is.data.frame(x)) {
            if (!all(cols_to_compare %in% colnames(x))) {
                stop(.missing_user_cols_error(), call. = FALSE)
            }
            ### Check they're numeric
            all_numeric <- purrr::map_lgl(cols_to_compare, function(col) {
                expr <- rlang::expr(`$`(x, !!col))
                is.numeric(rlang::eval_tidy(expr)) ||
                    is.integer(rlang::eval_tidy(expr))
            })
            if (!all(all_numeric)) {
                stop(.non_num_user_cols_error(), call. = FALSE)
            }
        } else {
            purrr::walk(x[names_to_mod], function(ddf) {
                if (!all(cols_to_compare %in% colnames(ddf))) {
                    stop(.missing_user_cols_error(), call. = FALSE)
                }
                purrr::walk(cols_to_compare, function(col) {
                    expr <- rlang::expr(`$`(ddf, !!col))
                    if (!is.numeric(rlang::eval_tidy(expr)) &&
                        !is.integer(rlang::eval_tidy(expr))) {
                        stop(.non_num_user_cols_error(),
                            call. = FALSE
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
#' @import dplyr
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

#---- USED IN : CIS_volcano_plot ----
#' @importFrom rlang arg_match
#' @importFrom utils read.delim
.load_onco_ts_genes <- function(onco_db_file,
    tumor_suppressors_db_file,
    species) {
    if (!file.exists(onco_db_file)) {
        stop(paste(
            "`onco_db_file` was not found, check you provided the",
            "correct path for the file"
        ))
    }
    if (!file.exists(tumor_suppressors_db_file)) {
        stop(paste(
            "`tumor_suppressors_db_file` was not found,",
            "check you provided the",
            "correct path for the file"
        ))
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
    # Acquire DB
    onco_db <- utils::read.delim(
        file = onco_db_file, header = TRUE,
        fill = TRUE, sep = "\t", check.names = FALSE
    )
    tumsup_db <- utils::read.delim(
        file = tumor_suppressors_db_file,
        header = TRUE,
        fill = TRUE, sep = "\t",
        check.names = FALSE
    )
    if (getOption("ISAnalytics.verbose") == TRUE) {
        message(paste(c(
            "Loading annotated genes -  species selected: ",
            paste(c(specie_name), collapse = ", ")
        )))
    }
    # Filter and merge
    onco_df <- .filter_db(onco_db, specie_name, "OncoGene")
    tumsup_df <- .filter_db(tumsup_db, specie_name, "TumorSuppressor")
    oncots_df_to_use <- .merge_onco_tumsup(onco_df, tumsup_df)
    if (getOption("ISAnalytics.verbose") == TRUE) {
        message(paste(c(
            "Loading annotated genes -  done"
        )))
    }
    return(oncots_df_to_use)
}

#' @import dplyr
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

#' @import dplyr
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
