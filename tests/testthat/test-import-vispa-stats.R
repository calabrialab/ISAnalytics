#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
withr::local_options(
    list(
        ISAnalytics.reports = FALSE,
        ISAnalytics.verbose = FALSE
    )
)

#------------------------------------------------------------------------------#
# Test .stats_report
#------------------------------------------------------------------------------#
test_that(".stats_report produces correct output", {
    output <- .stats_report(local_af_corr,
        prefixes = default_iss_file_prefixes(),
        proj_col = "ProjectID",
        pool_col = "concatenatePoolIDSeqRun",
        path_iss_col = .path_cols_names()$iss
    )
    expect_true(unique(output$ProjectID) == "PJ01")
    expect_true(nrow(output) == 4)
    expect_true(all(is.na(output[[.path_cols_names()$iss]]) == FALSE))
    expect_true(all(is.na(output$stats_files) == FALSE))

    output <- .stats_report(local_af_inc,
        prefixes = default_iss_file_prefixes(),
        proj_col = "ProjectID",
        pool_col = "concatenatePoolIDSeqRun",
        path_iss_col = .path_cols_names()$iss
    )
    expect_true(unique(output$ProjectID) == "PJ01")
    expect_true(nrow(output) == 4)
    expect_true(is.na(output %>%
        dplyr::filter(.data$concatenatePoolIDSeqRun == "POOL01-1") %>%
        dplyr::pull(.data[[.path_cols_names()$iss]])))
    expect_true(all(is.na(
        output %>%
            dplyr::filter(.data$concatenatePoolIDSeqRun %in% c(
                "POOL01-1",
                "POOL02-1"
            )) %>%
            dplyr::pull(.data$stats_files)
    )))
})

#------------------------------------------------------------------------------#
# Test .import_stats_iss
#------------------------------------------------------------------------------#
test_that(".import_stats_iss works as expected", {
    tags_to_cols <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "CompleteAmplificationID", "char", NULL, "required", "pcr_repl_id",
        "ProjectID", "char", NULL, "required", "project_id",
        "TagSequence", "char", NULL, "required", "tag_seq"
    )
    imported_iss <- .import_stats_iss(
        local_af_corr,
        default_iss_file_prefixes(),
        "concatenatePoolIDSeqRun",
        "Path_iss",
        tags_to_cols
    )
    expect_true(nrow(imported_iss$stats) == 53 & ncol(imported_iss$stats) == 14)
    expect_true(all(imported_iss$report$Imported == TRUE))
})

#------------------------------------------------------------------------------#
# Test import_Vispa2_stats
#------------------------------------------------------------------------------#
test_that("import_Vispa2_stats works with no join - default", {
    imported_stats <- import_Vispa2_stats(local_af_corr,
        join_with_af = FALSE
    )
    expect_true(nrow(imported_stats) == 53 & ncol(imported_stats) == 14)

    imported_stats <- import_Vispa2_stats(local_af_inc,
        join_with_af = FALSE
    )
    expect_true(nrow(imported_stats) == 27 & ncol(imported_stats) == 14)
})

test_that("import_Vispa2_stats works with join - default", {
    imported_stats <- import_Vispa2_stats(local_af_corr,
        join_with_af = TRUE
    )
    expect_true(nrow(imported_stats) == 53 & ncol(imported_stats) == 86)
    expect_true(all(!is.na(imported_stats$RUN_NAME)))
    imported_stats <- import_Vispa2_stats(local_af_inc,
        join_with_af = TRUE
    )
})
