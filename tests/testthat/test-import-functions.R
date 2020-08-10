context("Importing IS from files")

library(ISAnalytics)

#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
# Annotated new matrix
example_path <- system.file("extdata", "ex_annotated_ISMatrix.tsv.xz",
    package = "ISAnalytics"
)

# Old style not annotated matrix
example_path2 <- system.file("extdata", "ex_old_style_ISMatrix.tsv.xz",
    package = "ISAnalytics"
)

# New not annotated matrix
example_path3 <- system.file("extdata", "ex_notann_ISMatrix.tsv.xz",
    package = "ISAnalytics"
)

# Malformed matrix
example_path4 <- system.file("extdata", "ex_malformed_ISMatrix.tsv.xz",
    package = "ISAnalytics"
)

# Path to example association file
path_af <- system.file("extdata", "ex_association_file.tsv",
    package = "ISAnalytics"
)

# Path to correct file system example
path_root_correct <- system.file("extdata", "fs.zip",
    package = "ISAnalytics"
)
root_correct <- unzip_file_system(path_root_correct, "fs")

# Path to incorrect file system example
path_root_err <- system.file("extdata", "fserr.zip",
    package = "ISAnalytics"
)
root_err <- unzip_file_system(path_root_err, "fserr")

#------------------------------------------------------------------------------#
# Tests .messy_to_tidy
#------------------------------------------------------------------------------#
test_that(".messy_to_tidy fails if input is not an ISADataFrame", {
    expect_error(.messy_to_tidy(data.frame(a = seq_len(10), b = seq_len(10))))
    expect_error(.messy_to_tidy(tibble(list(a = seq_len(10), b = seq_len(10)))))
})

test_that(".messy_to_tidy produces correct tidy ISADataFrame - annotated new", {
    raw_isadf <- read.csv(example_path, sep = "\t", check.names = FALSE)
    raw_isadf <- new_ISADataFrame(raw_isadf, meta = c("GeneName", "GeneStrand"))
    tidy_isadf <- .messy_to_tidy(raw_isadf)
    # Resulting data frame must preserve ISADataFrame class
    expect_true(is.ISADataFrame(tidy_isadf))
    # Resulting ISADataFrame must have the columns "ExperimentID" and "Value"
    expect_true(all(c(
        "CompleteAmplificationID",
        "Value"
    ) %in% colnames(tidy_isadf)))
    # The metadata field updates and contains the new "CompleteAmplificationID"
    # column
    expect_equal(metadata(tidy_isadf), c(
        "GeneName", "GeneStrand",
        "CompleteAmplificationID"
    ))
    # The data frame must not contain null values
    expect_true(all(!is.na(tidy_isadf$Value)))
})

test_that(".messy_to_tidy produces correct tidy ISADataFrame - not annotated", {
    raw_isadf <- read.csv(example_path3, sep = "\t", check.names = FALSE)
    raw_isadf <- new_ISADataFrame(raw_isadf)
    tidy_isadf <- .messy_to_tidy(raw_isadf)
    # Resulting data frame must preserve ISADataFrame class
    expect_true(is.ISADataFrame(tidy_isadf))
    # Resulting ISADataFrame must have the columns "ExperimentID" and "Value"
    expect_true(all(c(
        "CompleteAmplificationID",
        "Value"
    ) %in% colnames(tidy_isadf)))
    # The metadata field updates and contains the new "CompleteAmplificationID"
    # column
    expect_equal(metadata(tidy_isadf), c("CompleteAmplificationID"))
    # The data frame must not contain null values
    expect_true(all(!is.na(tidy_isadf$Value)))
})

#------------------------------------------------------------------------------#
# Tests .auto_detect_type
#------------------------------------------------------------------------------#
test_that(".auto_detect_type detects new style annotated matrices", {
    raw_isadf <- read.csv(example_path, sep = "\t", check.names = FALSE)
    res <- .auto_detect_type(raw_isadf)
    expect_equal(res, "NEW_ANNOTATED")
})

test_that(".auto_detect_type detects old style matrices", {
    raw_isadf <- read.csv(example_path2, sep = "\t", check.names = FALSE)
    res <- .auto_detect_type(raw_isadf)
    expect_equal(res, "OLD")
})

test_that(".auto_detect_type detects new not annotated matrices", {
    raw_isadf <- read.csv(example_path, sep = "\t", check.names = FALSE)
    raw_isadf <- dplyr::select(raw_isadf, -c("GeneName", "GeneStrand"))
    res <- .auto_detect_type(raw_isadf)
    expect_equal(res, "NEW_NOTANN")
})

test_that(".auto_detect_type detects malformed matrices", {
    raw_isadf <- read.csv(example_path, sep = "\t", check.names = FALSE)
    raw_isadf1 <- dplyr::select(
        raw_isadf,
        -c("GeneName", "GeneStrand", "integration_locus")
    )
    res <- .auto_detect_type(raw_isadf1)
    expect_equal(res, "MALFORMED")
    raw_isadf2 <- read.csv(example_path2, sep = "\t", check.names = FALSE)
    raw_isadf3 <- tibble::add_column(raw_isadf2, chr = raw_isadf$chr)
    res <- .auto_detect_type(raw_isadf3)
    expect_equal(res, "MALFORMED")
})

#------------------------------------------------------------------------------#
# Tests import_single_Vispa2Matrix
#------------------------------------------------------------------------------#
test_that("import_single_Vispa2Matrix succeeds with standard annotated", {
    isadf <- import_single_Vispa2Matrix(example_path)
    expect_s3_class(isadf, "ISADataFrame")
    expect_equal(metadata(isadf), c(
        "GeneName", "GeneStrand",
        "CompleteAmplificationID"
    ))
    expect_true(all(c("chr", "integration_locus", "strand") %in% colnames(isadf)))
    expect_true(all(!is.na(isadf$Value)))
})

test_that("import_single_Vispa2Matrix succeeds with old", {
    isadf <- import_single_Vispa2Matrix(example_path2)
    expect_s3_class(isadf, "ISADataFrame")
    expect_equal(metadata(isadf), c("CompleteAmplificationID"))
    expect_true(all(c("chr", "integration_locus", "strand") %in% colnames(isadf)))
    expect_true(all(!is.na(isadf$Value)))
})

test_that("import_single_Vispa2Matrix succeeds with new not annotated", {
    isadf <- import_single_Vispa2Matrix(example_path3)
    expect_s3_class(isadf, "ISADataFrame")
    expect_equal(metadata(isadf), c("CompleteAmplificationID"))
    expect_true(all(c("chr", "integration_locus", "strand") %in% colnames(isadf)))
    expect_true(all(!is.na(isadf$Value)))
})

test_that("import_single_Vispa2Matrix fails with malformed matrix", {
    expect_error({
        isadf <- import_single_Vispa2Matrix(example_path4)
    })
})

test_that("import_single_Vispa2Matrix fails if file does not exist", {
    expect_error(import_single_Vispa2Matrix(system.file("extdata",
        "matrix.tsv",
        package = "ISAnalytics"
    )))
})


#------------------------------------------------------------------------------#
# Tests .check_af_correctness
#------------------------------------------------------------------------------#
test_that(".check_af_correctness returns TRUE for correct file", {
    af <- tibble::as_tibble(read.csv(path_af,
        check.names = FALSE, sep = "\t",
        stringsAsFactors = FALSE
    ))
    check <- .check_af_correctness(af)
    expect_true(check)
})

test_that(".check_af_correctness returns FALSE for malformed file", {
    af <- tibble::as_tibble(read.csv(path_af,
        check.names = FALSE, sep = "\t",
        stringsAsFactors = FALSE
    ))
    af <- af %>% dplyr::select(-c(.data$ProjectID)) # Missing a column
    check <- .check_af_correctness(af)
    expect_false(check)
})

#------------------------------------------------------------------------------#
# Tests .read_and_correctness_af
#------------------------------------------------------------------------------#
test_that(".read_and_correctness_af succeeds when association file is correct", {
    af <- .read_and_correctness_af(path_af)
    expect_s3_class(af, "tbl_df")
})

test_that(".read_and_correctness_af fails when association file is incorrect", {
    af <- tibble::as_tibble(read.csv(path_af,
        check.names = FALSE, sep = "\t",
        stringsAsFactors = FALSE
    ))
    af <- af %>% dplyr::select(-c(.data$ProjectID)) # Missing a column
    temp_path <- tempfile(fileext = ".tsv")
    readr::write_tsv(af, temp_path)
    expect_error({
        af <- .read_and_correctness_af(temp_path)
    })
})

#### OTHER VARS ####
af_df <- .read_and_correctness_af(path_af)
missing_project <- list(
    ProjectID = "PROJECT1101",
    concatenatePoolIDSeqRun = "ABY-LR-PL4-POOL54-2"
)


#------------------------------------------------------------------------------#
# Tests .check_file_system_alignment
#------------------------------------------------------------------------------#
test_that(".check_file_system_alignment finds all projects for correct fs", {
    checks <- .check_file_system_alignment(af_df, root_correct)
    expect_true(all(checks$Found == TRUE))
    expect_true(all(!is.na(checks$Path)))
})

test_that(".check_file_system_alignment finds missing projects for incorrect
          fs", {
    checks <- .check_file_system_alignment(af_df, root_err)
    expect_false(all(checks$Found == TRUE))
    expect_false(all(!is.na(checks$Path)))
    expect_equal(
        (checks %>%
            dplyr::filter(
                .data$ProjectID == missing_project$ProjectID,
                .data$concatenatePoolIDSeqRun ==
                    missing_project$concatenatePoolIDSeqRun
            ))$Found, FALSE
    )
    expect_equal(
        (checks %>%
            dplyr::filter(
                .data$ProjectID == missing_project$ProjectID,
                .data$concatenatePoolIDSeqRun ==
                    missing_project$concatenatePoolIDSeqRun
            ))$Path,
        NA_character_
    )
})


#------------------------------------------------------------------------------#
# Tests .update_af_after_alignment
#------------------------------------------------------------------------------#
test_that(".update_af_after_alignment updates correct fs - no NAs", {
    checks <- .check_file_system_alignment(af_df, root_correct)
    updated_af <- .update_af_after_alignment(af_df,
        checks = checks,
        root = root_correct
    )
    expect_true(all(!is.na(updated_af$Path)))
    expect_equal(colnames(updated_af), c(association_file_columns, "Path"))
})

test_that(".update_af_after_alignment updates correct fserr - with NAs", {
    checks <- .check_file_system_alignment(af_df, root_err)
    expect_warning({
        updated_af <- .update_af_after_alignment(af_df,
            checks = checks,
            root = root_err
        )
    })
    expect_false(all(!is.na(updated_af$Path)))
    expect_equal(colnames(updated_af), c(association_file_columns, "Path"))
    missing <- updated_af %>%
        dplyr::filter(
            .data$ProjectID == missing_project$ProjectID,
            .data$concatenatePoolIDSeqRun ==
                missing_project$concatenatePoolIDSeqRun
        )
    expect_true(all(is.na(missing$Path)))
})

#### OTHER VARS ####
warning_update_after_alignment <-
    paste0("One or more projects were not found in the file",
           "system starting from ", root_err, ", please check your",
           "association file for errors and/or your file system.",
           "Until you re-import the association file these ",
           "missing files will be ignored.",
           collapse = ""
    )

#------------------------------------------------------------------------------#
# Tests import_association_file
#------------------------------------------------------------------------------#
## Testing input
test_that("import_association_file fails if path is not a single string", {
    # Not a string
    expect_error(
        import_association_file(1, root_correct)
    )
    # Vector of strings
    expect_error(
        import_association_file(c("a", "b", "c"), root_correct)
    )
})
test_that("import_association_file fails if root is not a single string", {
    # Not a string
    expect_error(
        import_association_file(path_af, 1)
    )
    # Vector of strings
    expect_error(
        import_association_file(path_af, c("a", "b", "c"))
    )
})
test_that("import_association_file fails if path or root don't exist", {
    # Path doesn't exist
    expect_error(
        import_association_file("a", root_correct)
    )
    # Root doesn't exist
    expect_error(
        import_association_file(path_af, "a")
    )
})
## Testing results
test_that("import_association_file imports af with no warnings for fs", {
    op <- options("viewer" = NULL)
    on.exit(options(op), add = TRUE, after = FALSE)
    expect_warning(
        {
            af <- import_association_file(path_af, root_correct)
        },
        regexp = NA
    )
    expect_equal(colnames(af), c(association_file_columns, "Path"))
    expect_true(all(!is.na(af$Path)))
})

test_that("import_association_file imports af with warnings for fserr", {
    op <- options("viewer" = NULL)
    on.exit(options(op), add = TRUE, after = FALSE)
    expect_warning(
        {
            af <- import_association_file(path_af, root_err)
        },
        regexp = warning_update_after_alignment
    )
    expect_equal(colnames(af), c(association_file_columns, "Path"))
    expect_false(all(!is.na(af$Path)))
    missing <- af %>% dplyr::filter(
        .data$ProjectID == missing_project$ProjectID,
        .data$concatenatePoolIDSeqRun ==
            missing_project$concatenatePoolIDSeqRun
    )
    expect_true(all(is.na(missing$Path)))
})

#------------------------------------------------------------------------------#
# Tests .manage_association_file
#------------------------------------------------------------------------------#

#### OTHER VARS ####
all_proj_pools <- list(
    ProjectID = c("CLOEXP", "PROJECT1100", "PROJECT1101"),
    PoolID = c(
        "POOL6-1", "ABX-LR-PL5-POOL14-1",
        "ABX-LR-PL6-POOL15-1", "ABY-LR-PL4-POOL54-2"
    )
)

test_that(".manage_association_file succeeds with fs - path version", {
    expect_warning(
        {
            af <- .manage_association_file(path_af, root_correct)
        },
        regexp = NA
    )
    expect_true(all(!is.na(af[[1]]$Path)))
    expect_true(!is.null(af[[2]]))
    expect_true(all(all_proj_pools$ProjectID %in% af[[1]]$ProjectID))
    expect_true(all(all_proj_pools$PoolID %in% af[[1]]$concatenatePoolIDSeqRun))
})

test_that(".manage_association_file succeeds with fs - tibble version", {
    op <- options("viewer" = NULL)
    on.exit(options(op), add = TRUE, after = FALSE)
    expect_warning(
        {
            af <- import_association_file(path_af, root_correct)
        },
        regexp = NA
    )
    expect_warning(
        {
            af <- .manage_association_file(af, NULL)
        },
        regexp = NA
    )
    expect_true(all(!is.na(af[[1]]$Path)))
    expect_true(is.null(af[[2]]))
    expect_true(all(all_proj_pools$ProjectID %in% af[[1]]$ProjectID))
    expect_true(all(all_proj_pools$PoolID %in% af[[1]]$concatenatePoolIDSeqRun))
})

test_that(".manage_association_file succeeds with fserr - path version", {
    expect_warning(
        {
            af <- .manage_association_file(path_af, root_err)
        },
        regexp = warning_update_after_alignment
    )
    expect_true(all(!is.na(af[[1]]$Path)))
    expect_true(!is.null(af[[2]]))
    expect_true(all(!missing_project$ProjectID %in% af[[1]]$ProjectID))
    expect_true(all(!missing_project$concatenatePoolIDSeqRun %in%
        af[[1]]$concatenatePoolIDSeqRun))
    not_missing_proj <- all_proj_pools$ProjectID[!all_proj_pools$ProjectID %in%
        missing_project$ProjectID]
    not_missing_pool <- all_proj_pools$PoolID[!all_proj_pools$PoolID %in%
        missing_project$concatenatePoolIDSeqRun]
    expect_true(all(not_missing_proj %in% af[[1]]$ProjectID))
    expect_true(all(not_missing_pool %in% af[[1]]$concatenatePoolIDSeqRun))
})

test_that(".manage_association_file succeeds with fserr - tibble version", {
    op <- options("viewer" = NULL)
    on.exit(options(op), add = TRUE, after = FALSE)
    expect_warning(
        {
            af <- import_association_file(path_af, root_err)
        },
        regexp = warning_update_after_alignment
    )
    expect_warning(
        {
            af <- .manage_association_file(af, NULL)
        },
        regexp = NA
    )
    expect_true(all(!is.na(af[[1]]$Path)))
    expect_true(is.null(af[[2]]))
    expect_true(all(!missing_project$ProjectID %in% af[[1]]$ProjectID))
    expect_true(all(!missing_project$concatenatePoolIDSeqRun %in%
        af[[1]]$concatenatePoolIDSeqRun))
    not_missing_proj <- all_proj_pools$ProjectID[!all_proj_pools$ProjectID %in%
        missing_project$ProjectID]
    not_missing_pool <- all_proj_pools$PoolID[!all_proj_pools$PoolID %in%
        missing_project$concatenatePoolIDSeqRun]
    expect_true(all(not_missing_proj %in% af[[1]]$ProjectID))
    expect_true(all(not_missing_pool %in% af[[1]]$concatenatePoolIDSeqRun))
})

test_that(".manage_association_file fails if association file is malformed", {
    op <- options("viewer" = NULL)
    on.exit(options(op), add = TRUE, after = FALSE)
    expect_warning(
        {
            af <- import_association_file(path_af, root_correct)
        },
        regexp = NA
    )
    af <- af %>% dplyr::select(-c("Path"))
    expect_error({
        af <- .manage_association_file(af, NULL)
    })
    expect_warning(
        {
            af <- import_association_file(path_af, root_correct)
        },
        regexp = NA
    )
    af <- af %>% dplyr::select(-c("ProjectID"))
    expect_error({
        af <- .manage_association_file(af, NULL)
    })
})

#------------------------------------------------------------------------------#
# Tests .interactive_select_projects_import
#------------------------------------------------------------------------------#

#### OTHER VARS ####
as_file_correct <- .manage_association_file(path_af, root_correct)[[1]]
suppressWarnings({
    as_file_err <- .manage_association_file(path_af, root_err)[[1]]
})

test_that(".interactive_select_projects_import stops if input 1 is 0", {
    input1 <- tempfile()
    input_value1 <- "0"
    op <- options(ISAnalytics.connection = input1)
    on.exit(options(op), add = TRUE, after = FALSE)
    write(input_value1, input1)
    expect_error(
        {
            invisible(capture_output({
                .interactive_select_projects_import(as_file_correct)
            }))
        },
        regexp = "Quitting"
    )
})

test_that(".interactive_select_projects_import selects all projects
          when inputs are '1, y'", {
    input1 <- tempfile()
    input2 <- tempfile()
    input3 <- tempfile()
    input_value1 <- "1"
    input_value3 <- "y"
    op <- options(ISAnalytics.connection = c(input1, input2, input3))
    on.exit(options(op), add = TRUE, after = FALSE)
    write(input_value1, input1)
    write(input_value3, input3)
    invisible(capture_output({
        updated_af <- .interactive_select_projects_import(as_file_correct)
    }))
    expect_true(all(all_proj_pools$ProjectID %in% updated_af$ProjectID))
})

test_that(".interactive_select_projects_import selects only some projects
          when input 1 is '2'", {
    # Single choice
    input1 <- tempfile()
    input2 <- tempfile()
    input3 <- tempfile()
    input_value1 <- "2"
    input_value2 <- "1"
    input_value3 <- "y"
    op <- options(ISAnalytics.connection = c(input1, input2, input3))
    on.exit(options(op), add = TRUE, after = FALSE)
    write(input_value1, input1)
    write(input_value2, input2)
    write(input_value3, input3)
    project_list <- dplyr::distinct(
        dplyr::select(as_file_correct, .data$ProjectID)
    )$ProjectID
    invisible(capture_output({
        updated_af <- .interactive_select_projects_import(as_file_correct)
    }))
    selected_proj <- project_list[as.numeric(input_value2)]
    expect_true(all(updated_af$ProjectID %in% selected_proj))

    input_value1 <- "2"
    input_value2 <- "2"
    input_value3 <- "y"
    write(input_value1, input1)
    write(input_value2, input2)
    write(input_value3, input3)
    invisible(capture_output({
        updated_af <- .interactive_select_projects_import(as_file_correct)
    }))
    selected_proj <- project_list[as.numeric(input_value2)]
    expect_true(all(updated_af$ProjectID %in% selected_proj))

    input_value1 <- "2"
    input_value2 <- "3"
    input_value3 <- "y"
    write(input_value1, input1)
    write(input_value2, input2)
    write(input_value3, input3)
    invisible(capture_output({
        updated_af <- .interactive_select_projects_import(as_file_correct)
    }))
    selected_proj <- project_list[as.numeric(input_value2)]
    expect_true(all(updated_af$ProjectID %in% selected_proj))

    # Multiple choice
    input_value1 <- "2"
    input_value2 <- "1,2"
    input_value3 <- "y"
    write(input_value1, input1)
    write(input_value2, input2)
    write(input_value3, input3)
    invisible(capture_output({
        updated_af <- .interactive_select_projects_import(as_file_correct)
    }))
    indexes <- unlist(stringr::str_split(input_value2, ","))
    indexes <- as.numeric(indexes)
    selected_proj <- project_list[indexes]
    expect_true(all(updated_af$ProjectID %in% selected_proj))

    input_value1 <- "2"
    input_value2 <- "2,3"
    input_value3 <- "y"
    write(input_value1, input1)
    write(input_value2, input2)
    write(input_value3, input3)
    invisible(capture_output({
        updated_af <- .interactive_select_projects_import(as_file_correct)
    }))
    indexes <- unlist(stringr::str_split(input_value2, ","))
    indexes <- as.numeric(indexes)
    selected_proj <- project_list[indexes]
    expect_true(all(updated_af$ProjectID %in% selected_proj))
})

test_that(".interactive_select_projects_import stops if input 2 is 0", {
    input1 <- tempfile()
    input2 <- tempfile()
    input_value1 <- "2"
    input_value2 <- "0"
    op <- options(ISAnalytics.connection = c(input1, input2))
    on.exit(options(op), add = TRUE, after = FALSE)
    write(input_value1, input1)
    write(input_value2, input2)
    expect_error(
        {
            invisible(capture_output({
                .interactive_select_projects_import(as_file_correct)
            }))
        },
        regexp = "Quitting"
    )
})

#------------------------------------------------------------------------------#
# Tests .pool_number_IN
#------------------------------------------------------------------------------#
test_that(".pool_number_IN stops when input is 0", {
    input4 <- tempfile()
    op <- options(ISAnalytics.connection = c(NA, NA, NA, input4))
    on.exit(options(op), add = TRUE, after = FALSE)
    input_value4 <- "0"
    write(input_value4, input4)
    expect_error(
        {
            invisible(capture_output({
                .pool_number_IN()
            }))
        },
        regexp = "Quitting"
    )
})

test_that(".pool_number_IN returns user choice correctly when not 0", {
    input4 <- tempfile()
    op <- options(ISAnalytics.connection = c(NA, NA, NA, input4))
    on.exit(options(op), add = TRUE, after = FALSE)
    input_value4 <- "1"
    write(input_value4, input4)
    invisible(capture_output({
        choice <- .pool_number_IN()
    }))
    expect_equal(choice, 1)

    input_value4 <- "2"
    write(input_value4, input4)
    invisible(capture_output({
        choice <- .pool_number_IN()
    }))
    expect_equal(choice, 2)
})

#------------------------------------------------------------------------------#
# Tests .pool_choices_IN
#------------------------------------------------------------------------------#
test_that("pool_choices_IN stops when input is 0", {
    input5 <- tempfile()
    op <- options(ISAnalytics.connection = c(NA, NA, NA, NA, input5))
    on.exit(options(op), add = TRUE, after = FALSE)
    input_value5 <- "0"
    write(input_value5, input5)
    expect_error(
        {
            invisible(capture_output({
                .pool_choices_IN(c(1, 2, 3))
            }))
        },
        regexp = "Quitting"
    )
})

test_that(".pool_number_IN returns user choice correctly when not 0 - single", {
    input5 <- tempfile()
    op <- options(ISAnalytics.connection = c(NA, NA, NA, NA, input5))
    on.exit(options(op), add = TRUE, after = FALSE)
    input_value5 <- "1"
    write(input_value5, input5)
    invisible(capture_output({
        choice <- .pool_choices_IN(c(1, 2, 3))
    }))
    expect_true(choice == 1)

    input_value5 <- "2"
    write(input_value5, input5)
    invisible(capture_output({
        choice <- .pool_choices_IN(c(1, 2, 3))
    }))
    expect_true(choice == 2)

    input_value5 <- "3"
    write(input_value5, input5)
    invisible(capture_output({
        choice <- .pool_choices_IN(c(1, 2, 3))
    }))
    expect_true(choice == 3)
})

test_that(".pool_number_IN returns user choice correctly when not 0 - multi", {
    input5 <- tempfile()
    op <- options(ISAnalytics.connection = c(NA, NA, NA, NA, input5))
    on.exit(options(op), add = TRUE, after = FALSE)
    input_value5 <- "1,2"
    write(input_value5, input5)
    invisible(capture_output({
        choice <- .pool_choices_IN(c(1, 2, 3))
    }))
    expect_equal(choice, c("1" = 1, "2" = 2))

    input_value5 <- "2,3"
    write(input_value5, input5)
    invisible(capture_output({
        choice <- .pool_choices_IN(c(1, 2, 3))
    }))
    expect_equal(choice, c("2" = 2, "3" = 3))
})

#------------------------------------------------------------------------------#
# Tests .interactive_select_pools_import
#------------------------------------------------------------------------------#
test_that(".interactive_select_pools_import selects all pools if input is 1", {
    input4 <- tempfile()
    input6 <- tempfile()
    op <- options(ISAnalytics.connection = c(NA, NA, NA, input4, NA, input6))
    on.exit(options(op), add = TRUE, after = FALSE)
    input_value4 <- "1"
    input_value6 <- "y"
    write(input_value4, input4)
    write(input_value6, input6)
    invisible(capture_output({
        updated_af <- .interactive_select_pools_import(as_file_correct)
    }))
    expect_true(all(all_proj_pools$PoolID %in%
        updated_af$concatenatePoolIDSeqRun))
})

test_that(".interactive_select_pools_import selects only the first pools
          if input is '2, 1, y'", {
    input4 <- tempfile()
    input5 <- tempfile()
    input6 <- tempfile()
    op <- options(ISAnalytics.connection = c(NA, NA, NA, input4, input5, input6))
    on.exit(options(op), add = TRUE, after = FALSE)
    input_value4 <- "2"
    input_value5 <- "1"
    input_value6 <- "y"
    write(input_value4, input4)
    write(input_value5, input5)
    write(input_value6, input6)
    invisible(capture_output({
        updated_af <- .interactive_select_pools_import(as_file_correct)
    }))
    available <- as_file_correct %>%
        dplyr::select(.data$ProjectID, .data$concatenatePoolIDSeqRun) %>%
        dplyr::distinct()
    pools_selected <- purrr::map(
        unique(available$ProjectID),
        function(x) {
            pools <- (available %>%
                dplyr::filter(.data$ProjectID == x))$concatenatePoolIDSeqRun[1]
            pools
        }
    )
    pools_selected <- unlist(pools_selected)
    expect_true(all(updated_af$concatenatePoolIDSeqRun %in% pools_selected))
})

#------------------------------------------------------------------------------#
# Tests .trace_anomalies
#------------------------------------------------------------------------------#
smpl_ff <- function() {
    tb <- tibble::tibble(
        ProjectID = c(
            "PROJ1", "PROJ1", "PROJ2", "PROJ2",
            "PROJ3", "PROJ3"
        ),
        concatenatePoolIDSeqRun = c(
            "POOL1", "POOL1", "POOL2", "POOL2",
            "POOL3", "POOL3"
        ),
        Quantification_type = c(
            "fragmentEstimate", "seqCount",
            "fragmentEstimate", "seqCount",
            "fragmentEstimate", "seqCount"
        ),
        Files_found = c(
            "pth1", "pth2", "pth3", "pth4", "pth5",
            "pth6"
        )
    )
    tb <- tb %>% tidyr::nest(Files = c("Quantification_type", "Files_found"))
    tb
}

test_that(".trace_anomalies updates tibble correctly", {
    example_df <- smpl_ff()
    updated_df <- .trace_anomalies(example_df)
    expect_true(all(c("Anomalies", "Files_count") %in% colnames(updated_df)))
    expect_true(all(updated_df$Anomalies == FALSE))
    expect_true(all(updated_df$Files_count[[1]]$Found == 1))
    expect_true(all(updated_df$Files_count[[2]]$Found == 1))
    expect_true(all(updated_df$Files_count[[3]]$Found == 1))
})

#------------------------------------------------------------------------------#
# Tests .lookup_matrices
#------------------------------------------------------------------------------#

#### OTHER VARS ####
quant_types <- c("fragmentEstimate", "seqCount")

test_that(".lookup_matrices finds all files for fs - no duplicates", {
    files_found <- .lookup_matrices(as_file_correct, quant_types, "annotated")
    expect_true(all(c(
        "ProjectID", "concatenatePoolIDSeqRun", "Anomalies",
        "Files", "Files_count"
    ) %in% colnames(files_found)))
    expect_true(all(files_found$Anomalies == FALSE))
    expect_true(all(files_found$Files_count[[1]]$Found == 1))
    expect_true(all(files_found$Files_count[[2]]$Found == 1))
    expect_true(all(files_found$Files_count[[3]]$Found == 1))
    expect_true(all(files_found$Files_count[[4]]$Found == 1))

    files_found <- .lookup_matrices(as_file_correct, quant_types, "not_annotated")
    expect_true(all(c(
        "ProjectID", "concatenatePoolIDSeqRun", "Anomalies",
        "Files", "Files_count"
    ) %in% colnames(files_found)))
    expect_true(all(files_found$Anomalies == FALSE))
    expect_true(all(files_found$Files_count[[1]]$Found == 1))
    expect_true(all(files_found$Files_count[[2]]$Found == 1))
    expect_true(all(files_found$Files_count[[3]]$Found == 1))
    expect_true(all(files_found$Files_count[[4]]$Found == 1))
})

test_that(".lookup_matrices finds all files for fserr - duplicates", {
    # Annotated matrices
    files_found <- .lookup_matrices(as_file_err, quant_types, "annotated")
    expect_true(all(c(
        "ProjectID", "concatenatePoolIDSeqRun", "Anomalies",
        "Files", "Files_count"
    ) %in% colnames(files_found)))
    expect_false(all(files_found$Anomalies == FALSE))
    expect_true((files_found %>%
        dplyr::filter(
            .data$ProjectID == "CLOEXP",
            .data$concatenatePoolIDSeqRun == "POOL6-1"
        )
    )$Anomalies == TRUE)
    cloexp_anomalies <- (files_found %>%
        dplyr::filter(
            .data$ProjectID == "CLOEXP",
            .data$concatenatePoolIDSeqRun == "POOL6-1"
        ))$Files_count[[1]]
    expect_equal(dplyr::filter(
        cloexp_anomalies, .data$Quantification_type == "fragmentEstimate"
    )$Found, 2)
    expect_equal(dplyr::filter(
        cloexp_anomalies, .data$Quantification_type == "seqCount"
    )$Found, 2)
    proj1100_1_anomalies <- (files_found %>%
        dplyr::filter(
            .data$ProjectID == "PROJECT1100",
            .data$concatenatePoolIDSeqRun == "ABX-LR-PL5-POOL14-1"
        ))$Files_count[[1]]
    expect_equal(dplyr::filter(
        proj1100_1_anomalies,
        .data$Quantification_type == "fragmentEstimate"
    )$Found, 3)
    expect_equal(dplyr::filter(
        proj1100_1_anomalies, .data$Quantification_type == "seqCount"
    )$Found, 1)
    # Non annotated matrices
    files_found <- .lookup_matrices(as_file_err, quant_types, "not_annotated")
    expect_true(all(c(
        "ProjectID", "concatenatePoolIDSeqRun", "Anomalies",
        "Files", "Files_count"
    ) %in% colnames(files_found)))
    expect_false(all(files_found$Anomalies == FALSE))
    expect_true((files_found %>%
        dplyr::filter(
            .data$ProjectID == "CLOEXP",
            .data$concatenatePoolIDSeqRun == "POOL6-1"
        )
    )$Anomalies == FALSE)
    expect_true((files_found %>%
        dplyr::filter(
            .data$ProjectID == "PROJECT1100",
            .data$concatenatePoolIDSeqRun ==
                "ABX-LR-PL5-POOL14-1"
        )
    )$Anomalies == TRUE)
    expect_true((files_found %>%
        dplyr::filter(
            .data$ProjectID == "PROJECT1100",
            .data$concatenatePoolIDSeqRun ==
                "ABX-LR-PL6-POOL15-1"
        )
    )$Anomalies == TRUE)
    proj1100_1_anomalies <- (files_found %>%
        dplyr::filter(
            .data$ProjectID == "PROJECT1100",
            .data$concatenatePoolIDSeqRun == "ABX-LR-PL5-POOL14-1"
        ))$Files_count[[1]]
    expect_equal(dplyr::filter(
        proj1100_1_anomalies,
        .data$Quantification_type == "fragmentEstimate"
    )$Found, 0)
    expect_equal(dplyr::filter(
        proj1100_1_anomalies, .data$Quantification_type == "seqCount"
    )$Found, 0)
    proj1100_2_anomalies <- (files_found %>%
        dplyr::filter(
            .data$ProjectID == "PROJECT1100",
            .data$concatenatePoolIDSeqRun == "ABX-LR-PL6-POOL15-1"
        ))$Files_count[[1]]
    expect_equal(dplyr::filter(
        proj1100_2_anomalies,
        .data$Quantification_type == "fragmentEstimate"
    )$Found, 0)
    expect_equal(dplyr::filter(
        proj1100_2_anomalies, .data$Quantification_type == "seqCount"
    )$Found, 0)
})

#------------------------------------------------------------------------------#
# Tests .choose_duplicates_files_interactive
#------------------------------------------------------------------------------#
smpl_dupl <- function() {
    tb <- tibble::tibble(
        Quantification_type = c(
            "fragmentEstimate", "seqCount",
            "fragmentEstimate"
        ),
        Files_found = c("pth1", "pth2", "pth3")
    )
    tb
}

test_that(".choose_duplicates_files_interactive modifies tbl according to user
          choices", {
    input7 <- tempfile()
    op <- options(ISAnalytics.connection = c(
        NA, NA, NA, NA, NA, NA,
        input7
    ))
    on.exit(options(op), add = TRUE, after = FALSE)
    input_value7 <- "1"
    write(input_value7, input7)
    smpl_dupl_t <- smpl_dupl()
    dupl_quant <- "fragmentEstimate"
    invisible(capture_output({
        updated_tbl <- .choose_duplicates_files_interactive(
            dupl_quant,
            smpl_dupl_t
        )
    }))
    expect_true(all(c("fragmentEstimate", "seqCount") %in%
        updated_tbl$Quantification_type))
    filtered <- updated_tbl %>%
        dplyr::filter(.data$Quantification_type == dupl_quant)
    expect_equal(nrow(filtered), 1)
    expect_equal(filtered$Files_found, c("pth1"))
})

test_that(".choose_duplicates_files_interactive stops if input is 0", {
    input7 <- tempfile()
    op <- options(ISAnalytics.connection = c(
        NA, NA, NA, NA, NA, NA,
        input7
    ))
    on.exit(options(op), add = TRUE, after = FALSE)
    input_value7 <- "0"
    write(input_value7, input7)
    smpl_dupl_t <- smpl_dupl()
    dupl_quant <- "fragmentEstimate"
    expect_error(
        {
            invisible(capture_output({
                updated_tbl <- .choose_duplicates_files_interactive(
                    dupl_quant,
                    smpl_dupl_t
                )
            }))
        },
        regexp = "Quitting"
    )
})

#------------------------------------------------------------------------------#
# Tests .manage_anomalies_interactive
#------------------------------------------------------------------------------#

#### OTHER VARS ####
files_found_fs_ann <- .lookup_matrices(
    as_file_correct,
    quant_types, "annotated"
)
files_found_fs_notann <- .lookup_matrices(
    as_file_correct,
    quant_types, "not_annotated"
)
files_found_fs_err_ann <- .lookup_matrices(
    as_file_err,
    quant_types, "annotated"
)
files_found_fs_err_notann <- .lookup_matrices(
    as_file_err,
    quant_types, "not_annotated"
)

test_that(".manage_anomalies_interactive succeeds with no anomalies", {
    # Annotated matrices
    expect_message(
        {
            invisible(capture_output({
                to_import <- .manage_anomalies_interactive(files_found_fs_ann)
            }))
        },
        regexp = NA
    )
    expect_true(all(c("Files_chosen") %in% colnames(to_import)))
    expect_true(all(!is.na(to_import$Files_chosen)))
    expect_true(all(quant_types %in% to_import$Quantification_type))
    present <- purrr::map2_lgl(
        to_import$ProjectID,
        to_import$concatenatePoolIDSeqRun,
        function(x, y) {
            temp <- to_import %>%
                dplyr::filter(
                    .data$ProjectID == x,
                    .data$concatenatePoolIDSeqRun ==
                        y
                )
            all(quant_types %in% temp$Quantification_type)
        }
    )
    expect_true(all(present))

    # Not annotated matrices
    expect_message(
        {
            invisible(capture_output({
                to_import <- .manage_anomalies_interactive(files_found_fs_notann)
            }))
        },
        regexp = NA
    )
    expect_true(all(c("Files_chosen") %in% colnames(to_import)))
    expect_true(all(!is.na(to_import$Files_chosen)))
    expect_true(all(quant_types %in% to_import$Quantification_type))
    present <- purrr::map2_lgl(
        to_import$ProjectID,
        to_import$concatenatePoolIDSeqRun,
        function(x, y) {
            temp <- to_import %>%
                dplyr::filter(
                    .data$ProjectID == x,
                    .data$concatenatePoolIDSeqRun ==
                        y
                )
            all(quant_types %in% temp$Quantification_type)
        }
    )
    expect_true(all(present))
})

test_that(".manage_anomalies_interactive succeeds with anomalies", {
    input7 <- tempfile()
    input8 <- tempfile()
    op <- options(ISAnalytics.connection = c(
        NA, NA, NA, NA, NA, NA, input7,
        input8
    ))
    on.exit(options(op), add = TRUE, after = FALSE)
    input_value7 <- "1"
    input_value8 <- "y"
    write(input_value7, input7)
    write(input_value8, input8)
    # Annotated matrices
    expect_message(
        {
            invisible(capture_output({
                to_import <- .manage_anomalies_interactive(files_found_fs_err_ann)
            }))
        },
        regexp = "Duplicates found for some files"
    )
    expect_true(all(c("Files_chosen") %in% colnames(to_import)))
    expect_true(all(!is.na(to_import$Files_chosen)))
    expect_true(all(quant_types %in% to_import$Quantification_type))
    present <- purrr::map2_lgl(
        to_import$ProjectID,
        to_import$concatenatePoolIDSeqRun,
        function(x, y) {
            temp <- to_import %>%
                dplyr::filter(
                    .data$ProjectID == x,
                    .data$concatenatePoolIDSeqRun ==
                        y
                )
            all(quant_types %in% temp$Quantification_type)
        }
    )
    expect_true(all(present))

    # Not annotated matrices
    expect_message(
        {
            invisible(capture_output({
                to_import <- .manage_anomalies_interactive(files_found_fs_err_notann)
            }))
        },
        regexp = "Some files are missing and will be ignored"
    )
    expect_true(all(c("Files_chosen") %in% colnames(to_import)))
    expect_true(all(!is.na(to_import$Files_chosen)))
    expect_true(all(
        quant_types %in% (to_import %>%
            dplyr::filter(
                .data$ProjectID == "CLOEXP",
                .data$concatenatePoolIDSeqRun ==
                    "POOL6-1"
            ))$Quantification_type
    ))
    expect_true(all( # Missing not annotated matrices from PROJECT1100
        !quant_types %in% (to_import %>%
            dplyr::filter(
                .data$ProjectID == "PROJECT1100",
                .data$concatenatePoolIDSeqRun ==
                    "ABX-LR-PL5-POOL14-1"
            )
        )$Quantification_type
    ))
    expect_true(all(
        !quant_types %in% (to_import %>%
            dplyr::filter(
                .data$ProjectID == "PROJECT1100",
                .data$concatenatePoolIDSeqRun ==
                    "ABX-LR-PL6-POOL15-1"
            )
        )$Quantification_type
    ))
})

#------------------------------------------------------------------------------#
# Tests .import_type
#------------------------------------------------------------------------------#
#### OTHER VARS ####
invisible(capture_output({
    smpl_files_to_import <- .manage_anomalies_interactive(files_found_fs_ann)
}))

test_that(".import_type correctly imports all files", {
    imports <- .import_type("fragmentEstimate", files = smpl_files_to_import, 2)
    matrices <- imports[[1]]
    report <- imports[[2]]
    expect_true(all(report$Imported == TRUE))
    expect_true(nrow(matrices) == 4 * 1476)

    imports <- .import_type("seqCount", files = smpl_files_to_import, 2)
    matrices <- imports[[1]]
    report <- imports[[2]]
    expect_true(all(report$Imported == TRUE))
    expect_true(nrow(matrices) == 4 * 1476)
})

#------------------------------------------------------------------------------#
# Tests .parallel_import_merge
#------------------------------------------------------------------------------#
test_that(".parallel_import_merge correctly imports files for each q type", {
    imported <- .parallel_import_merge(smpl_files_to_import, 2)
    matrices <- imported[[1]]
    report <- imported[[2]]
    expect_true(all(report$Imported == TRUE))
    expect_true(nrow(matrices$fragmentEstimate) == 4 * 1476)
    expect_true(nrow(matrices$seqCount) == 4 * 1476)
})

#------------------------------------------------------------------------------#
# Tests import_parallel_Vispa2Matrices_interactive
#------------------------------------------------------------------------------#
## Testing input
test_that("import_parallel_Vispa2Matrices_interactive stops if af is missing", {
    expect_error({
        import_parallel_Vispa2Matrices_interactive(
            root = "x",
            quantification_type =
                quant_types,
            matrix_type = "annotated",
            workers = 2
        )
    })
})

test_that("import_parallel_Vispa2Matrices_interactive stops if af is not
          character or tibble", {
    expect_error({
        import_parallel_Vispa2Matrices_interactive(
            association_file = 1,
            root = "x",
            quantification_type =
                quant_types,
            matrix_type = "annotated",
            workers = 2
        )
    })
})

test_that("import_parallel_Vispa2Matrices_interactive stops if af is a character
          vector longer than 1", {
    expect_error({
        import_parallel_Vispa2Matrices_interactive(
            association_file = c(
                path_af,
                path_af
            ),
            root = "x",
            quantification_type =
                quant_types,
            matrix_type = "annotated",
            workers = 2
        )
    })
})

test_that("import_parallel_Vispa2Matrices_interactive stops if af is a character
          but root is incorrect", {
    expect_error({
        import_parallel_Vispa2Matrices_interactive(
            association_file = path_af,
            quantification_type =
                quant_types,
            matrix_type = "annotated",
            workers = 2
        )
    })
    expect_error({
        import_parallel_Vispa2Matrices_interactive(
            association_file = path_af,
            root = 1,
            quantification_type =
                quant_types,
            matrix_type = "annotated",
            workers = 2
        )
    })
    expect_error({
        import_parallel_Vispa2Matrices_interactive(
            association_file = path_af,
            root = c(root_correct, root_err),
            quantification_type =
                quant_types,
            matrix_type = "annotated",
            workers = 2
        )
    })
})

test_that("import_parallel_Vispa2Matrices_interactive stops if workers is not
         numeric or have multiple values", {
    expect_error({
        import_parallel_Vispa2Matrices_interactive(
            association_file = path_af,
            root = root_correct,
            quantification_type =
                quant_types,
            matrix_type = "annotated",
            workers = "x"
        )
    })
    expect_error({
        import_parallel_Vispa2Matrices_interactive(
            association_file = path_af,
            root = root_correct,
            quantification_type =
                quant_types,
            matrix_type = "annotated",
            workers = c(1, 2, 3)
        )
    })
})


test_that("import_parallel_Vispa2Matrices_interactive stops if
         quantification_types is incorrect or missing", {
    expect_error({
        import_parallel_Vispa2Matrices_interactive(
            association_file = path_af,
            root = root_correct,
            matrix_type = "annotated",
            workers = 2
        )
    })
    expect_error({
        import_parallel_Vispa2Matrices_interactive(
            association_file = path_af,
            root = root_correct,
            quantification_type =
                c("a"),
            matrix_type = "annotated",
            workers = 2
        )
    })
})

test_that("import_parallel_Vispa2Matrices_interactive stops if
          matrix_type is incorrect", {
    expect_error({
        import_parallel_Vispa2Matrices_interactive(
            association_file = path_af,
            root = root_correct,
            quantification_type =
                quant_types,
            matrix_type = 1,
            workers = 2
        )
    })
    expect_error({
        import_parallel_Vispa2Matrices_interactive(
            association_file = path_af,
            root = root_correct,
            quantification_type =
                quant_types,
            matrix_type = "ann",
            workers = 2
        )
    })
})
## Testing results
test_that("import_parallel_Vispa2Matrices_interactive succeeds for fs - path", {
    # Annotated
    input1 <- tempfile()
    input3 <- tempfile()
    input4 <- tempfile()
    input6 <- tempfile()
    op1 <- options(ISAnalytics.connection = c(
        input1, NA, input3, input4, NA,
        input6, NA, NA
    ))
    op2 <- options("viewer" = NULL)
    on.exit(
        {
            options(op1)
            options(op2)
        },
        add = TRUE,
        after = FALSE
    )
    input_value1 <- "1"
    input_value3 <- "y"
    input_value4 <- "1"
    input_value6 <- "y"
    write(input_value1, input1)
    write(input_value3, input3)
    write(input_value4, input4)
    write(input_value6, input6)
    invisible(capture_output({
        matrices <-
            import_parallel_Vispa2Matrices_interactive(path_af,
                root_correct,
                quantification_type =
                    quant_types,
                "annotated", 2
            )
    }))
    expect_equal(nrow(matrices$fragmentEstimate), 4 * 1476)
    expect_equal(nrow(matrices$seqCount), 4 * 1476)
    # Not annotated
    input_value1 <- "1"
    input_value3 <- "y"
    input_value4 <- "1"
    input_value6 <- "y"
    write(input_value1, input1)
    write(input_value3, input3)
    write(input_value4, input4)
    write(input_value6, input6)
    invisible(capture_output({
        matrices <-
            import_parallel_Vispa2Matrices_interactive(path_af,
                root_correct,
                quantification_type =
                    quant_types,
                "not_annotated", 2
            )
    }))
    expect_equal(nrow(matrices$fragmentEstimate), 4 * 1476)
    expect_equal(nrow(matrices$seqCount), 4 * 1476)
})

test_that("import_parallel_Vispa2Matrices_interactive succeeds
          for fs - tibble", {
    # Annotated
    input1 <- tempfile()
    input3 <- tempfile()
    input4 <- tempfile()
    input6 <- tempfile()
    op1 <- options(ISAnalytics.connection = c(
        input1, NA, input3, input4, NA,
        input6, NA, NA
    ))
    op2 <- options("viewer" = NULL)
    on.exit(
        {
            options(op1)
            options(op2)
        },
        add = TRUE,
        after = FALSE
    )
    input_value1 <- "1"
    input_value3 <- "y"
    input_value4 <- "1"
    input_value6 <- "y"
    write(input_value1, input1)
    write(input_value3, input3)
    write(input_value4, input4)
    write(input_value6, input6)
    invisible(capture_output({
        matrices <-
            import_parallel_Vispa2Matrices_interactive(as_file_correct,
                NULL,
                quantification_type =
                    quant_types,
                "annotated", 2
            )
    }))
    expect_equal(nrow(matrices$fragmentEstimate), 4 * 1476)
    expect_equal(nrow(matrices$seqCount), 4 * 1476)
    # Not annotated
    input_value1 <- "1"
    input_value3 <- "y"
    input_value4 <- "1"
    input_value6 <- "y"
    write(input_value1, input1)
    write(input_value3, input3)
    write(input_value4, input4)
    write(input_value6, input6)
    invisible(capture_output({
        matrices <-
            import_parallel_Vispa2Matrices_interactive(as_file_correct,
                NULL,
                quantification_type =
                    quant_types,
                "not_annotated", 2
            )
    }))
    expect_equal(nrow(matrices$fragmentEstimate), 4 * 1476)
    expect_equal(nrow(matrices$seqCount), 4 * 1476)
})

test_that("import_parallel_Vispa2Matrices_interactive succeeds for
          fserr - path", {
    # Annotated
    input1 <- tempfile()
    input3 <- tempfile()
    input4 <- tempfile()
    input6 <- tempfile()
    input7 <- tempfile()
    input8 <- tempfile()
    op1 <- options(ISAnalytics.connection = c(
        input1, NA, input3, input4, NA,
        input6, input7, input8
    ))
    op2 <- options("viewer" = NULL)
    on.exit(
        {
            options(op1)
            options(op2)
        },
        add = TRUE,
        after = FALSE
    )
    input_value1 <- "1"
    input_value3 <- "y"
    input_value4 <- "1"
    input_value6 <- "y"
    input_value7 <- "1"
    input_value8 <- "y"
    write(input_value1, input1)
    write(input_value3, input3)
    write(input_value4, input4)
    write(input_value6, input6)
    write(input_value7, input7)
    write(input_value8, input8)
    expect_warning({
        invisible(capture_output({
            matrices <-
                import_parallel_Vispa2Matrices_interactive(path_af,
                    root_err,
                    quantification_type =
                        quant_types,
                    "annotated", 2
                )
        }))
    })
    expect_equal(nrow(matrices$fragmentEstimate), 3 * 1476)
    expect_equal(nrow(matrices$seqCount), 3 * 1476)
    # Not annotated
    input_value1 <- "1"
    input_value3 <- "y"
    input_value4 <- "1"
    input_value6 <- "y"
    input_value7 <- "1"
    input_value8 <- "y"
    write(input_value1, input1)
    write(input_value3, input3)
    write(input_value4, input4)
    write(input_value6, input6)
    write(input_value7, input7)
    write(input_value8, input8)
    expect_warning({
        invisible(capture_output({
            matrices <-
                import_parallel_Vispa2Matrices_interactive(path_af,
                    root_err,
                    quantification_type =
                        quant_types,
                    "not_annotated", 2
                )
        }))
    })
    expect_equal(nrow(matrices$fragmentEstimate), 1 * 1476)
    expect_equal(nrow(matrices$seqCount), 1 * 1476)
})

#------------------------------------------------------------------------------#
# Tests .pattern_matching
#------------------------------------------------------------------------------#
test_that(".pattern_matching matches single pattern", {
    p_matches <- .pattern_matching(c("file1", "file2", "file3test"), "file")
    expect_true(all(colnames(p_matches) %in% c("file")))
    expect_true(all(p_matches$file == TRUE))
})

test_that(".pattern_matching matches multiple patterns", {
    p_matches <- .pattern_matching(c("file1", "file2", "file3test"), c(
        "file",
        "1",
        "test"
    ))
    expect_true(all(colnames(p_matches) %in% c("file", "1", "test")))
    expect_true(all(p_matches$file == TRUE))
    expect_equal(p_matches$`1`, c(TRUE, FALSE, FALSE))
    expect_equal(p_matches$test, c(FALSE, FALSE, TRUE))
})

#------------------------------------------------------------------------------#
# Tests .any_match & .all_match
#------------------------------------------------------------------------------#
test_that(".any_match provides correct output", {
    test_values <- list(TRUE, FALSE, FALSE)
    any <- .any_match(test_values)
    all <- .all_match(test_values)
    expect_true(any)
    expect_false(all)
})

#------------------------------------------------------------------------------#
# Tests .update_as_option
#------------------------------------------------------------------------------#
test_that(".update_as_option works for fs err - single pattern", {
    cloexp_pool6 <- (files_found_fs_err_ann %>%
        dplyr::filter(
            .data$ProjectID == "CLOEXP",
            .data$concatenatePoolIDSeqRun == "POOL6-1"
        ))$Files[[1]]
    cloexp_matches <- .pattern_matching(cloexp_pool6$Files_found, "NoMate")
    updated <- .update_as_option(cloexp_pool6, cloexp_matches, "ANY")
    expect_true(nrow(updated) == 2)
    expect_true(all(quant_types %in% updated$Quantification_type))
    expect_true(all(!is.na(updated$Files_found)))
    updated <- .update_as_option(cloexp_pool6, cloexp_matches, "ALL")
    expect_true(nrow(updated) == 2)
    expect_true(all(quant_types %in% updated$Quantification_type))
    expect_true(all(!is.na(updated$Files_found)))
    updated <- .update_as_option(cloexp_pool6, cloexp_matches, "OPTIONAL")
    expect_true(nrow(updated) == 2)
    expect_true(all(quant_types %in% updated$Quantification_type))
    expect_true(all(!is.na(updated$Files_found)))
})

test_that(".update_as_option works for fs err - multiple patterns", {
    cloexp_pool6 <- (files_found_fs_err_ann %>%
        dplyr::filter(
            .data$ProjectID == "CLOEXP",
            .data$concatenatePoolIDSeqRun == "POOL6-1"
        )
    )$Files[[1]]
    cloexp_matches <- .pattern_matching(cloexp_pool6$Files_found, c(
        "NoMate",
        "test"
    ))
    updated <- .update_as_option(cloexp_pool6, cloexp_matches, "ANY")
    expect_true(nrow(updated) == 2)
    expect_true(all(quant_types %in% updated$Quantification_type))
    expect_true(all(!is.na(updated$Files_found)))
    updated <- .update_as_option(cloexp_pool6, cloexp_matches, "ALL")
    expect_true(nrow(updated) == 2)
    expect_true(all(is.na(updated$Files_found)))
    expect_true(all(quant_types %in% updated$Quantification_type))
    updated <- .update_as_option(cloexp_pool6, cloexp_matches, "OPTIONAL")
    expect_true(nrow(updated) == 2)
    expect_true(all(quant_types %in% updated$Quantification_type))
    expect_true(all(!is.na(updated$Files_found)))

    proj1100_pool1 <- (files_found_fs_err_ann %>%
        dplyr::filter(
            .data$ProjectID == "PROJECT1100",
            .data$concatenatePoolIDSeqRun ==
                "ABX-LR-PL5-POOL14-1"
        )
    )$Files[[1]]
    proj1100_pool1_matches <- .pattern_matching(
        proj1100_pool1$Files_found,
        c("NoMate", "test")
    )
    updated <- .update_as_option(proj1100_pool1, proj1100_pool1_matches, "ANY")
    expect_true(nrow(updated) == 3)
    expect_true(all(quant_types %in% updated$Quantification_type))
    expect_true(all(!is.na((updated %>%
        dplyr::filter(
            .data$Quantification_type == "fragmentEstimate"
        )
    )$Files_found)))
    expect_true(all(is.na((updated %>%
        dplyr::filter(
            .data$Quantification_type == "seqCount"
        )
    )$Files_found)))
    updated <- .update_as_option(proj1100_pool1, proj1100_pool1_matches, "ALL")
    expect_true(nrow(updated) == 2)
    expect_true(all(quant_types %in% updated$Quantification_type))
    expect_true(all(is.na(updated$Files_found)))
    proj1100_pool1_matches <- .pattern_matching(
        proj1100_pool1$Files_found,
        c("NoMate", "test")
    )
    updated <- .update_as_option(
        proj1100_pool1, proj1100_pool1_matches,
        "OPTIONAL"
    )
    expect_true(nrow(updated) == 3)
    expect_true(all(quant_types %in% updated$Quantification_type))
    expect_true(all(!is.na(updated$Files_found)))

    proj1100_pool2 <- (files_found_fs_err_ann %>%
        dplyr::filter(
            .data$ProjectID == "PROJECT1100",
            .data$concatenatePoolIDSeqRun ==
                "ABX-LR-PL6-POOL15-1"
        )
    )$Files[[1]]
    proj1100_pool2_matches <- .pattern_matching(
        proj1100_pool2$Files_found,
        c("NoMate", "test")
    )
    updated <- .update_as_option(proj1100_pool2, proj1100_pool2_matches, "ANY")
    expect_true(nrow(updated) == 2)
    expect_true(all(quant_types %in% updated$Quantification_type))
    expect_true(all(is.na(updated$Files_found)))
    updated <- .update_as_option(proj1100_pool2, proj1100_pool2_matches, "ALL")
    expect_true(nrow(updated) == 2)
    expect_true(all(quant_types %in% updated$Quantification_type))
    expect_true(all(is.na(updated$Files_found)))
    updated <- .update_as_option(
        proj1100_pool2, proj1100_pool2_matches,
        "OPTIONAL"
    )
    expect_true(nrow(updated) == 2)
    expect_true(all(quant_types %in% updated$Quantification_type))
    expect_true(all(!is.na(updated$Files_found)))
})

#------------------------------------------------------------------------------#
# Tests .lookup_matrices_auto
#------------------------------------------------------------------------------#
test_that(".lookup_matrices_auto returns correctly with patterns null", {
    to_import <- .lookup_matrices_auto(as_file_correct, quant_types, "annotated",
        patterns = NULL, matching_opt = "ANY"
    )
    expect_true(all(to_import$Anomalies == FALSE))
    to_import <- .lookup_matrices_auto(as_file_correct, quant_types, "annotated",
        patterns = NULL, matching_opt = "ALL"
    )
    expect_true(all(to_import$Anomalies == FALSE))
    to_import <- .lookup_matrices_auto(as_file_correct, quant_types, "annotated",
        patterns = NULL, matching_opt = "OPTIONAL"
    )
    expect_true(all(to_import$Anomalies == FALSE))
})

test_that(".lookup_matrices_auto returns correctly for fs", {
    to_import <- .lookup_matrices_auto(as_file_correct, quant_types, "annotated",
        patterns = "NoMate", matching_opt = "ANY"
    )
    expect_true(all(to_import$Anomalies == TRUE))
    to_import <- .lookup_matrices_auto(as_file_correct, quant_types, "annotated",
        patterns = "NoMate", matching_opt = "ALL"
    )
    expect_true(all(to_import$Anomalies == TRUE))
    to_import <- .lookup_matrices_auto(as_file_correct, quant_types, "annotated",
        patterns = "NoMate",
        matching_opt = "OPTIONAL"
    )
    expect_true(all(to_import$Anomalies == FALSE))
})

test_that(".lookup_matrices_auto returns correctly for fserr", {
    to_import <- .lookup_matrices_auto(as_file_err, quant_types, "annotated",
        patterns = "NoMate", matching_opt = "ANY"
    )
    expect_true(all((to_import %>%
        dplyr::filter(.data$ProjectID == "PROJECT1100")
    )$Anomalies) == TRUE)
    expect_true(all((to_import %>%
        dplyr::filter(.data$ProjectID == "CLOEXP")
    )$Anomalies) == FALSE)
    to_import <- .lookup_matrices_auto(as_file_err, quant_types, "annotated",
        patterns = "NoMate", matching_opt = "ALL"
    )
    expect_true(all((to_import %>%
        dplyr::filter(.data$ProjectID == "PROJECT1100")
    )$Anomalies) == TRUE)
    expect_true(all((to_import %>%
        dplyr::filter(.data$ProjectID == "CLOEXP")
    )$Anomalies) == FALSE)
    to_import <- .lookup_matrices_auto(as_file_err, quant_types, "annotated",
        patterns = "NoMate",
        matching_opt = "OPTIONAL"
    )
    expect_true(all((to_import %>%
        dplyr::filter(
            .data$ProjectID == "PROJECT1100",
            .data$concatenatePoolIDSeqRun ==
                "ABX-LR-PL5-POOL14-1"
        )
    )$Anomalies) == TRUE)
})

#------------------------------------------------------------------------------#
# Tests .manage_anomalies_auto
#------------------------------------------------------------------------------#
test_that(".manage_anomalies_auto correctly manages for fs", {
    files_found <- .lookup_matrices_auto(as_file_correct,
        quant_types, "annotated",
        patterns = NULL, matching_opt = "ANY"
    )
    expect_message(
        {
            files_to_import <- .manage_anomalies_auto(files_found)
        },
        regexp = NA
    )
    expect_true(nrow(files_to_import) == 4 * length(quant_types))
    files_found <- .lookup_matrices_auto(as_file_correct,
        quant_types, "annotated",
        patterns = "NoMate",
        matching_opt = "ANY"
    )
    expect_message(
        {
            files_to_import <- .manage_anomalies_auto(files_found)
        },
        regexp = "Some files are missing and will be ignored"
    )
    expect_true(nrow(files_to_import) == 0)
    files_found <- .lookup_matrices_auto(as_file_correct,
        quant_types, "annotated",
        patterns = "NoMate",
        matching_opt = "ALL"
    )
    expect_message(
        {
            files_to_import <- .manage_anomalies_auto(files_found)
        },
        regexp = "Some files are missing and will be ignored"
    )
    expect_true(nrow(files_to_import) == 0)
    files_found <- .lookup_matrices_auto(as_file_correct,
        quant_types, "annotated",
        patterns = "NoMate",
        matching_opt = "OPTIONAL"
    )
    expect_message(
        {
            files_to_import <- .manage_anomalies_auto(files_found)
        },
        regexp = NA
    )
    expect_true(nrow(files_to_import) == 4 * length(quant_types))
})

test_that(".manage_anomalies_auto correctly manages for fserr", {
    dupl_message <- paste(
        "Duplicates found for some files: in automatic mode",
        "duplicates are not preserved - use interactive mode",
        "for more accurate",
        "file selection"
    )
    missing_message <- "Some files are missing and will be ignored"
    files_found <- .lookup_matrices_auto(as_file_err,
        quant_types, "annotated",
        patterns = NULL, matching_opt = "ANY"
    )
    expect_message(
        {
            files_to_import <- .manage_anomalies_auto(files_found)
        },
        regexp = dupl_message
    )
    expect_true(nrow(files_to_import) == 3)
    files_found <- .lookup_matrices_auto(as_file_err,
        quant_types, "annotated",
        patterns = "NoMate",
        matching_opt = "ANY"
    )
    expect_message(
        {
            files_to_import <- .manage_anomalies_auto(files_found)
        },
        regexp = c(missing_message)
    )
    expect_true(nrow(files_to_import) == 2)
    files_found <- .lookup_matrices_auto(as_file_err,
        quant_types, "annotated",
        patterns = "NoMate",
        matching_opt = "ALL"
    )
    expect_message(
        {
            files_to_import <- .manage_anomalies_auto(files_found)
        },
        regexp = c(missing_message)
    )
    expect_true(nrow(files_to_import) == 2)
    files_found <- .lookup_matrices_auto(as_file_err,
        quant_types, "annotated",
        patterns = "NoMate",
        matching_opt = "OPTIONAL"
    )
    expect_message(
        {
            files_to_import <- .manage_anomalies_auto(files_found)
        },
        regexp = c(dupl_message)
    )
    expect_true(nrow(files_to_import) == 5)
})

#------------------------------------------------------------------------------#
# Tests import_parallel_Vispa2Matrices_auto
#------------------------------------------------------------------------------#
## Testing input
test_that("import_parallel_Vispa2Matrices_auto stops if af is missing", {
    expect_error({
        import_parallel_Vispa2Matrices_auto(
            root = "x",
            quantification_type =
                quant_types,
            matrix_type = "annotated",
            workers = 2,
            patterns = NULL,
            matching_opt = "ANY"
        )
    })
})

test_that("import_parallel_Vispa2Matrices_auto stops if af is not
          character or tibble", {
    expect_error({
        import_parallel_Vispa2Matrices_auto(
            association_file = 1,
            root = "x",
            quantification_type =
                quant_types,
            matrix_type = "annotated",
            workers = 2,
            patterns = NULL,
            matching_opt = "ANY"
        )
    })
})

test_that("import_parallel_Vispa2Matrices_auto stops if af is a character
          vector longer than 1", {
    expect_error({
        import_parallel_Vispa2Matrices_auto(
            association_file = c(
                path_af,
                path_af
            ),
            root = "x",
            quantification_type =
                quant_types,
            matrix_type = "annotated",
            workers = 2,
            patterns = NULL,
            matching_opt = "ANY"
        )
    })
})

test_that("import_parallel_Vispa2Matrices_auto stops if af is a character
          but root is incorrect", {
    expect_error({
        import_parallel_Vispa2Matrices_auto(
            association_file = path_af,
            quantification_type =
                quant_types,
            matrix_type = "annotated",
            workers = 2,
            patterns = NULL,
            matching_opt = "ANY"
        )
    })
    expect_error({
        import_parallel_Vispa2Matrices_auto(
            association_file = path_af,
            root = 1,
            quantification_type =
                quant_types,
            matrix_type = "annotated",
            workers = 2,
            patterns = NULL,
            matching_opt = "ANY"
        )
    })
    expect_error({
        import_parallel_Vispa2Matrices_auto(
            association_file = path_af,
            root = c(root_correct, root_err),
            quantification_type =
                quant_types,
            matrix_type = "annotated",
            workers = 2,
            patterns = NULL,
            matching_opt = "ANY"
        )
    })
})

test_that("import_parallel_Vispa2Matrices_auto stops if workers is not
         numeric or have multiple values", {
    expect_error({
        import_parallel_Vispa2Matrices_auto(
            association_file = path_af,
            root = root_correct,
            quantification_type =
                quant_types,
            matrix_type = "annotated",
            workers = "x",
            patterns = NULL,
            matching_opt = "ANY"
        )
    })
    expect_error({
        import_parallel_Vispa2Matrices_auto(
            association_file = path_af,
            root = root_correct,
            quantification_type =
                quant_types,
            matrix_type = "annotated",
            workers = c(1, 2, 3),
            patterns = NULL,
            matching_opt = "ANY"
        )
    })
})


test_that("import_parallel_Vispa2Matrices_auto stops if
         quantification_types is incorrect or missing", {
    expect_error({
        import_parallel_Vispa2Matrices_auto(
            association_file = path_af,
            root = root_correct,
            matrix_type = "annotated",
            workers = 2,
            patterns = NULL,
            matching_opt = "ANY"
        )
    })
    expect_error({
        import_parallel_Vispa2Matrices_auto(
            association_file = path_af,
            root = root_correct,
            quantification_type =
                c("a"),
            matrix_type = "annotated",
            workers = 2,
            patterns = NULL,
            matching_opt = "ANY"
        )
    })
})

test_that("import_parallel_Vispa2Matrices_auto stops if
          matrix_type is incorrect", {
    expect_error({
        import_parallel_Vispa2Matrices_auto(
            association_file = path_af,
            root = root_correct,
            quantification_type =
                quant_types,
            matrix_type = 1,
            workers = 2,
            patterns = NULL,
            matching_opt = "ANY"
        )
    })
    expect_error({
        import_parallel_Vispa2Matrices_auto(
            association_file = path_af,
            root = root_correct,
            quantification_type =
                quant_types,
            matrix_type = "ann",
            workers = 2,
            patterns = NULL,
            matching_opt = "ANY"
        )
    })
})

test_that("import_parallel_Vispa2Matrices_auto stops if
          patterns is incorrect", {
    expect_error({
        import_parallel_Vispa2Matrices_auto(
            association_file = path_af,
            root = root_correct,
            quantification_type =
                quant_types,
            matrix_type = 1,
            workers = 2,
            patterns = 1,
            matching_opt = "ANY"
        )
    })
})

test_that("import_parallel_Vispa2Matrices_auto stops if
           matching_option is incorrect", {
    expect_error({
        import_parallel_Vispa2Matrices_auto(
            association_file = path_af,
            root = root_correct,
            quantification_type =
                quant_types,
            matrix_type = 1,
            workers = 2,
            patterns = NULL,
            matching_opt = "SOME"
        )
    })
    expect_error({
        import_parallel_Vispa2Matrices_auto(
            association_file = path_af,
            root = root_correct,
            quantification_type =
                quant_types,
            matrix_type = 1,
            workers = 2,
            patterns = NULL,
            matching_opt = c("ANY", "ALL")
        )
    })
    expect_error({
        import_parallel_Vispa2Matrices_auto(
            association_file = path_af,
            root = root_correct,
            quantification_type =
                quant_types,
            matrix_type = 1,
            workers = 2,
            patterns = NULL
        )
    })
})
## Testing results
test_that("import_parallel_Vispa2Matrices_auto succeeds for fs - path", {
    op <- options("viewer" = NULL)
    on.exit(options(op), add = TRUE, after = FALSE)
    matrices <- import_parallel_Vispa2Matrices_auto(path_af,
        root_correct,
        quant_types,
        "annotated",
        2,
        patterns = NULL,
        matching_opt = "ANY"
    )
    expect_true(nrow(matrices$fragmentEstimate) == 4 * 1476)
    expect_true(nrow(matrices$seqCount) == 4 * 1476)
})

test_that("import_parallel_Vispa2Matrices_auto succeeds for fs - tibble", {
    op <- options("viewer" = NULL)
    on.exit(options(op), add = TRUE, after = FALSE)
    af <- as_file_correct %>% dplyr::filter(.data$ProjectID == "CLOEXP")
    matrices <- import_parallel_Vispa2Matrices_auto(af,
        NULL,
        quant_types,
        "annotated",
        2,
        patterns = NULL,
        matching_opt = "ANY"
    )
    expect_true(nrow(matrices$fragmentEstimate) == 1476)
    expect_true(nrow(matrices$seqCount) == 1476)
})

test_that("import_parallel_Vispa2Matrices_auto succeeds for fserr - tibble", {
    op <- options("viewer" = NULL)
    on.exit(options(op), add = TRUE, after = FALSE)
    af <- as_file_err %>% dplyr::filter(.data$ProjectID == "CLOEXP")
    matrices <- import_parallel_Vispa2Matrices_auto(af,
        NULL,
        quant_types,
        "annotated",
        2,
        patterns = "NoMate",
        matching_opt = "ANY"
    )
    expect_true(nrow(matrices$fragmentEstimate) == 1476)
    expect_true(nrow(matrices$seqCount) == 1476)
    matrices <- import_parallel_Vispa2Matrices_auto(af,
        NULL,
        quant_types,
        "annotated",
        2,
        patterns = "NoMate",
        matching_opt = "ALL"
    )
    expect_true(nrow(matrices$fragmentEstimate) == 1476)
    expect_true(nrow(matrices$seqCount) == 1476)
    matrices <- import_parallel_Vispa2Matrices_auto(af,
        NULL,
        quant_types,
        "annotated",
        2,
        patterns = "NoMate",
        matching_opt = "OPTIONAL"
    )
    expect_true(nrow(matrices$fragmentEstimate) == 1476)
    expect_true(nrow(matrices$seqCount) == 1476)
    af <- as_file_err %>%
        dplyr::filter(
            .data$ProjectID == "PROJECT1100",
            .data$concatenatePoolIDSeqRun == "ABX-LR-PL5-POOL14-1"
        )
    expect_error({
        matrices <- import_parallel_Vispa2Matrices_auto(af,
            NULL,
            quant_types,
            "annotated",
            2,
            patterns = "NoMate",
            matching_opt = "ANY"
        )
    })
    expect_error({
        matrices <- import_parallel_Vispa2Matrices_auto(af,
            NULL,
            quant_types,
            "annotated",
            2,
            patterns = "NoMate",
            matching_opt = "ALL"
        )
    })
    matrices <- import_parallel_Vispa2Matrices_auto(af,
        NULL,
        quant_types,
        "annotated",
        2,
        patterns = "NoMate",
        matching_opt = "OPTIONAL"
    )
    expect_true(nrow(matrices$seqCount) == 1476)
    af <- as_file_err %>%
        dplyr::filter(
            .data$ProjectID == "PROJECT1100",
            .data$concatenatePoolIDSeqRun == "ABX-LR-PL6-POOL15-1"
        )
    expect_error({
        matrices <- import_parallel_Vispa2Matrices_auto(af,
            NULL,
            quant_types,
            "annotated",
            2,
            patterns = "NoMate",
            matching_opt = "ANY"
        )
    })
    expect_error({
        matrices <- import_parallel_Vispa2Matrices_auto(af,
            NULL,
            quant_types,
            "annotated",
            2,
            patterns = "NoMate",
            matching_opt = "ALL"
        )
    })
    matrices <- import_parallel_Vispa2Matrices_auto(af,
        NULL,
        quant_types,
        "annotated",
        2,
        patterns = "NoMate",
        matching_opt = "OPTIONAL"
    )
    expect_true(nrow(matrices$fragmentEstimate) == 1476)
    expect_true(nrow(matrices$seqCount) == 1476)
})
