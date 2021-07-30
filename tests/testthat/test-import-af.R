library(ISAnalytics)
func_name <- c(
    ".read_af", ".check_file_system_alignment",
    ".update_af_after_alignment",
    "import_association_file"
)

#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
withr::local_options(list(ISAnalytics.reports = FALSE))
test_af_path <- system.file("testdata", "ex_association_file.tsv.gz",
    package = "ISAnalytics"
)

# Path to correct file system example
path_root_correct <- system.file("testdata", "fs.zip",
    package = "ISAnalytics"
)
root_correct <- unzip_file_system(path_root_correct, "fs")

# Path to incorrect file system example
path_root_err <- system.file("testdata", "fserr.zip",
    package = "ISAnalytics"
)
root_err <- unzip_file_system(path_root_err, "fserr")

path_cols <- .path_cols_names()
#------------------------------------------------------------------------------#
# Tests .read_af
#------------------------------------------------------------------------------#
test_that(paste(func_name[1], "works for tsv file"), {
    read_res <- .read_af(test_af_path,
        padding = 4, date_format = "dmy",
        delimiter = "\t"
    )
    expect_equal(dim(read_res$af), c(49, 66))
    expect_true(nrow(read_res$probs) == 0)
    expect_true(purrr::is_empty(read_res$date_fail))
})

#------------------------------------------------------------------------------#
# Tests .check_file_system_alignment
#------------------------------------------------------------------------------#
test_that(paste(func_name[2], "finds all projects for correct fs"), {
    read_res <- .read_af(test_af_path,
        padding = 4, date_format = "dmy",
        delimiter = "\t"
    )
    checks <- .check_file_system_alignment(read_res$af, root_correct)
    expect_true(all(checks$Found == TRUE))
    expect_true(all(!is.na(checks[[path_cols$project]])))
    expect_true(all(!is.na(checks[[path_cols$quant]])))
    expect_true(all(!is.na(checks[[path_cols$iss]])))
})

test_that(paste(func_name[2], "finds missing projects for incorrect fs"), {
    missing_project <- list(
        ProjectID = "PROJECT1101",
        concatenatePoolIDSeqRun = "ABY-LR-PL4-POOL54-2"
    )
    missing_iss <- list(
        ProjectID = "PROJECT1100",
        concatenatePoolIDSeqRun = "ABX-LR-PL6-POOL15-1"
    )
    read_res <- .read_af(test_af_path,
        padding = 4, date_format = "dmy",
        delimiter = "\t"
    )
    checks <- .check_file_system_alignment(read_res$af, root_err)
    expect_false(all(checks$Found == TRUE))
    expect_false(all(!is.na(checks[[path_cols$project]])))
    missing_row <- checks %>%
        dplyr::filter(
            .data$ProjectID == missing_project$ProjectID,
            .data$concatenatePoolIDSeqRun ==
                missing_project$concatenatePoolIDSeqRun
        )
    expect_equal(
        missing_row$Found, FALSE
    )
    expect_true(is.na(missing_row[[path_cols$project]]))
    expect_true(
        is.na(missing_row[[path_cols$quant]])
    )
    expect_true(
        is.na(missing_row[[path_cols$iss]])
    )
    missing_row_iss <- checks %>%
        dplyr::filter(
            .data$ProjectID == missing_iss$ProjectID,
            .data$concatenatePoolIDSeqRun ==
                missing_iss$concatenatePoolIDSeqRun
        )
    expect_true(
        is.na(missing_row_iss[[path_cols$iss]])
    )
})


#------------------------------------------------------------------------------#
# Tests .update_af_after_alignment
#------------------------------------------------------------------------------#
test_that(paste(func_name[3], "updates correct fs - no NAs"), {
    read_res <- .read_af(test_af_path,
        padding = 4, date_format = "dmy",
        delimiter = "\t"
    )
    checks <- .check_file_system_alignment(read_res$af, root_correct)
    updated_af <- .update_af_after_alignment(read_res$af,
        checks = checks,
        root = root_correct
    )
    expect_true(all(!is.na(updated_af$Path)))
})


test_that(paste(func_name[3], "updates correct fserr - with NAs"), {
    missing_project <- list(
        ProjectID = "PROJECT1101",
        concatenatePoolIDSeqRun = "ABY-LR-PL4-POOL54-2"
    )
    read_res <- .read_af(test_af_path,
        padding = 4, date_format = "dmy",
        delimiter = "\t"
    )
    checks <- .check_file_system_alignment(read_res$af, root_err)
    updated_af <- .update_af_after_alignment(read_res$af,
        checks = checks,
        root = root_err
    )

    expect_false(all(!is.na(updated_af$Path)))
    missing <- updated_af %>%
        dplyr::filter(
            .data$ProjectID == missing_project$ProjectID,
            .data$concatenatePoolIDSeqRun ==
                missing_project$concatenatePoolIDSeqRun
        )
    expect_true(all(is.na(missing$Path)))
})

#------------------------------------------------------------------------------#
# Tests import_association_file
#------------------------------------------------------------------------------#
## Testing input
test_that(paste(func_name[4], "fails if path is not a single string"), {
    # Not a string
    expect_error(
        import_association_file(1, root_correct)
    )
    # Vector of strings
    expect_error(
        import_association_file(c("a", "b", "c"), root_correct)
    )
})

test_that(paste(func_name[4], "fails if root is not a single string"), {
    # Not a string
    expect_error(
        import_association_file(test_af_path, 1)
    )
    # Vector of strings
    expect_error(
        import_association_file(test_af_path, c("a", "b", "c"))
    )
})

test_that(paste(func_name[4], "fails if path or root don't exist"), {
    # Path doesn't exist
    expect_error(
        import_association_file("a", root_correct)
    )
    # Root doesn't exist
    expect_error(
        import_association_file(test_af_path, "a")
    )
})

test_that(paste(func_name[4], "fails if padding is incorrect"), {
    # Padding is not numeric or integer
    expect_error(
        import_association_file(test_af_path, root_correct, tp_padding = "3")
    )
    # Padding is vector
    expect_error(
        import_association_file(test_af_path, root_correct, tp_padding = c(1, 2))
    )
})

test_that(paste(func_name[4], "fails if date format is incorrect"), {
    # dates_format is not a character
    expect_error(
        import_association_file(test_af_path, root_correct,
            tp_padding = 4, dates_format = 1
        )
    )
    # dates_format has not length 1
    expect_error(
        import_association_file(test_af_path, root_correct,
            tp_padding = 4, dates_format = c("dmy", "myd")
        )
    )
    # dates_format is not one of the allowed
    expect_error(
        import_association_file(test_af_path, root_correct,
            tp_padding = 4, dates_format = "ggmmaaa"
        )
    )
})

test_that(paste(func_name[4], "fails if filter is not a named list"), {
    expect_error(
        import_association_file(test_af_path, filter_for = c("a", "b"))
    )
    expect_error(
        import_association_file(test_af_path, filter_for = list("a", "b"))
    )
})

test_that(paste(func_name[4], "fails if separator is wrong"), {
    expect_error(
        import_association_file(test_af_path, separator = NULL)
    )
    expect_error(
        import_association_file(test_af_path, separator = c("\t", ","))
    )
    expect_error(
        import_association_file(test_af_path, separator = 1)
    )
})

## Testing results
test_that(paste(func_name[4], "imports with defaults"), {
    expect_message(
        {
            af1 <- import_association_file(test_af_path)
        },
        class = "summary_report"
    )
    expect_message(
        {
            af2 <- import_association_file(test_af_path, root = root_correct)
        },
        class = "summary_report"
    )
    expect_message(
        {
            af3 <- import_association_file(test_af_path, root = root_err)
        },
        class = "summary_report"
    )
    ### Check all have same dim
    expect_true(all(dim(af1) == c(49, 68)))
    expect_true(all(c(dim(af2), dim(af3)) == c(49, 71)))
})

test_that(paste(func_name[4], "imports with filtering"), {
    expect_message(
        {
            af1 <- import_association_file(test_af_path,
                filter_for = list(ProjectID = "CLOEXP")
            )
        },
        class = "summary_report"
    )
    expect_message(
        {
            af2 <- import_association_file(test_af_path,
                root = root_correct,
                filter_for = list(ProjectID = "CLOEXP")
            )
        },
        class = "summary_report"
    )
    expect_message(
        {
            af3 <- import_association_file(test_af_path,
                root = root_err,
                filter_for = list(ProjectID = "CLOEXP")
            )
        },
        class = "summary_report"
    )
    ### Check all have same dim
    expect_true(unique(af1$ProjectID) == "CLOEXP")
    expect_true(unique(af2$ProjectID) == "CLOEXP")
    expect_true(unique(af3$ProjectID) == "CLOEXP")
})
