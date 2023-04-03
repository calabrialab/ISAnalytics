#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
path_cols <- .path_cols_names()
withr::local_options(list(ISAnalytics.reports = FALSE))

#------------------------------------------------------------------------------#
# .read_af
#------------------------------------------------------------------------------#
test_that(".read_af works for tsv file", {
    read_res <- .read_af(
        path = fs_path$af, date_format = "ymd",
        delimiter = "\t"
    )
    expect_equal(dim(read_res$af), c(53, 71))
    expect_true(nrow(read_res$probs) == 0)
    expect_true(purrr::is_empty(read_res$date_fail))
})

test_that(".read_af works for xslx file", {
    skip_if_not_installed("readxl")
    skip_if_not_installed("openxlsx")
    cols_selected <- colnames(
        association_file
    )[colnames(association_file) %in% association_file_columns()]
    association_file_mod <- association_file |>
        dplyr::select(dplyr::all_of(cols_selected))
    withr::with_tempfile(
        new = "af", fileext = ".xlsx", pattern = "asso_file",
        code = {
            wb <- openxlsx::createWorkbook()
            openxlsx::addWorksheet(wb, "AF")
            options(openxlsx.dateFormat = "yyyy-mm-dd")
            openxlsx::writeData(wb, "AF", association_file_mod)
            openxlsx::saveWorkbook(wb, af, overwrite = TRUE)
            expect_message(
                {
                    read_res <- .read_af(
                        path = af,
                        date_format = "ymd",
                        delimiter = "\t"
                    )
                },
                class = "xls_file"
            )
            expect_equal(dim(read_res$af), c(53, 71))
        }
    )
})

#------------------------------------------------------------------------------#
# .check_file_system_alignment
#------------------------------------------------------------------------------#
test_that(".check_file_system_alignment works for correct fs", {
    checks <- .check_file_system_alignment(
        df = association_file,
        root_folder = fs_path$root_corr,
        proj_fold_col = "PathToFolderProjectID",
        concat_pool_col = "concatenatePoolIDSeqRun",
        project_id_col = "ProjectID"
    )
    expect_true(all(checks$Found))
    expect_true(all(!is.na(checks$Path)))
    expect_true(all(!is.na(checks$Path_quant)))
    expect_true(all(!is.na(checks$Path_iss)))
})

test_that(".check_file_system_alignment works for incorrect fs", {
    checks <- .check_file_system_alignment(
        df = association_file,
        root_folder = fs_path$root_inc,
        proj_fold_col = "PathToFolderProjectID",
        concat_pool_col = "concatenatePoolIDSeqRun",
        project_id_col = "ProjectID"
    )
    expect_true(all(checks$Found))
    expect_true(all(!is.na(checks$Path)))
    expect_true(is.na(checks$Path_quant[1]))
    expect_true(is.na(checks$Path_iss[1]))
})

test_that(".check_file_system_alignment no error for NA concat", {
    withr::with_dir(tempdir(), {
        fake_pj <- fs::path("Project")
        fs::dir_create(fake_pj)
        fs::dir_create(fs::path(fake_pj, "fold_1_NA"))
        fs::dir_create(fs::path(fake_pj, "fold_2_NA"))
        fs::dir_create(fs::path(fake_pj, "fold_3_NA"))
        test_df <- tibble::tribble(
            ~ProjectID, ~PoolID, ~concatenatePoolIDSeqRun,
            ~PathToFolderProjectID,
            "Project", "POOL1", NA_character_, "/Project",
            "Project", "POOL2", NA_character_, "/Project",
            "Project", "POOL3", NA_character_, "/Project",
        )
        expect_message(
            {
                checks <- .check_file_system_alignment(
                    df = test_df,
                    root_folder = ".",
                    proj_fold_col = "PathToFolderProjectID",
                    concat_pool_col = "concatenatePoolIDSeqRun",
                    project_id_col = "ProjectID"
                )
            },
            class = "na_concat"
        )
        expect_true(is.data.frame(checks))
    })
})

#------------------------------------------------------------------------------#
# .update_af_after_alignment
#------------------------------------------------------------------------------#
test_that(".update_af_after_alignment updates correct fs - no NAs", {
    checks <- .check_file_system_alignment(
        df = association_file,
        root_folder = fs_path$root_corr,
        proj_fold_col = "PathToFolderProjectID",
        concat_pool_col = "concatenatePoolIDSeqRun",
        project_id_col = "ProjectID"
    )
    updated_af <- .update_af_after_alignment(association_file,
        checks = checks,
        proj_fold_col = "PathToFolderProjectID",
        concat_pool_col = "concatenatePoolIDSeqRun",
        project_id_col = "ProjectID"
    )
    expect_true(all(!is.na(updated_af[[path_cols$project]])))
    expect_true(all(!is.na(updated_af[[path_cols$quant]])))
    expect_true(all(!is.na(updated_af[[path_cols$iss]])))
})


test_that(".update_af_after_alignment updates correct fserr - with NAs", {
    checks <- .check_file_system_alignment(
        df = association_file,
        root_folder = fs_path$root_inc,
        proj_fold_col = "PathToFolderProjectID",
        concat_pool_col = "concatenatePoolIDSeqRun",
        project_id_col = "ProjectID"
    )
    updated_af <- .update_af_after_alignment(association_file,
        checks = checks,
        proj_fold_col = "PathToFolderProjectID",
        concat_pool_col = "concatenatePoolIDSeqRun",
        project_id_col = "ProjectID"
    )
    expect_true(all(
        is.na(updated_af |>
            dplyr::filter(.data$concatenatePoolIDSeqRun == "POOL01-1") |>
            dplyr::pull(.data[[path_cols$quant]]))
    ))
    expect_true(all(
        is.na(updated_af |>
            dplyr::filter(.data$concatenatePoolIDSeqRun == "POOL01-1") |>
            dplyr::pull(.data[[path_cols$iss]]))
    ))
})

#------------------------------------------------------------------------------#
# .manage_association_file
#------------------------------------------------------------------------------#
test_that(".manage_association_file works as expected", {
    managed_af <- .manage_association_file(
        af_path = fs_path$af,
        root = fs_path$root_corr,
        delimiter = "\t",
        format = "ymd",
        filter = NULL,
        proj_fold_col = "PathToFolderProjectID",
        concat_pool_col = "concatenatePoolIDSeqRun",
        project_id_col = "ProjectID"
    )
    expect_true(nrow(managed_af$date_probs) == 0)
    expect_true(nrow(managed_af$parsing_probs) == 0)
    expect_true(is.data.frame(managed_af$check) & nrow(managed_af$check) == 4)
    expect_true(all(path_cols %in% colnames(managed_af$af)))
    expect_message(
        {
            managed_af <- .manage_association_file(
                af_path = fs_path$af,
                root = fs_path$root_corr,
                delimiter = "\t",
                format = "ymd",
                filter = list(X = "A"),
                proj_fold_col = "PathToFolderProjectID",
                concat_pool_col = "concatenatePoolIDSeqRun",
                project_id_col = "ProjectID"
            )
        },
        class = "filter_warn"
    )
})

#------------------------------------------------------------------------------#
# import_association_file
#------------------------------------------------------------------------------#
## Testing input

test_that("import_association_file fails if date format is incorrect", {
    # dates_format is not a character
    expect_error(
        import_association_file(fs_path$af, fs_path$root_corr,
            dates_format = 1
        )
    )
    # dates_format has not length 1
    expect_error(
        import_association_file(fs_path$af, fs_path$root_corr,
            dates_format = c("dmy", "myd")
        )
    )
    # dates_format is not one of the allowed
    expect_error(
        import_association_file(fs_path$af, fs_path$root_corr,
            dates_format = "ggmmaaa"
        )
    )
})

test_that("import_association_file fails if filter is not a named list", {
    expect_error(
        import_association_file(fs_path$af, filter_for = c("a", "b"))
    )
    expect_error(
        import_association_file(fs_path$af, filter_for = list("a", "b"))
    )
})

test_that("import_association_file fails if separator is wrong", {
    expect_error(
        import_association_file(fs_path$af, separator = NULL)
    )
    expect_error(
        import_association_file(fs_path$af, separator = 1)
    )
})

## Testing results
test_that("import_association_file default - no fs align", {
    expect_message(
        {
            af <- import_association_file(fs_path$af,
                convert_tp = FALSE
            )
        },
        class = "summary_report"
    )
    expect_true(nrow(af) == 53 & ncol(af) == 71)

    # With date conversion but no transformation
    expect_message(
        {
            af <- import_association_file(fs_path$af,
                convert_tp = TRUE,
                transformations = NULL
            )
        },
        class = "summary_report"
    )
    expect_true(nrow(af) == 53 & ncol(af) == 71)
    expect_true(all(c("TimepointMonths", "TimepointYears") %in% colnames(af)))
    expect_true(is.numeric(af$TimepointMonths) & is.numeric(af$TimepointYears))

    # With date conversion and transformation
    expect_message(
        {
            af <- import_association_file(fs_path$af,
                convert_tp = TRUE
            )
        },
        class = "summary_report"
    )
    expect_true(nrow(af) == 53 & ncol(af) == 71)
    expect_true(all(c("TimepointMonths", "TimepointYears") %in% colnames(af)))
    expect_true(is.character(af$TimepointMonths) &
        is.character(af$TimepointYears))
})

test_that("import_association_file default - with fs align", {
    expect_message(
        {
            af <- import_association_file(fs_path$af,
                root = fs_path$root_corr,
                convert_tp = FALSE
            )
        },
        class = "summary_report"
    )
    expect_true(nrow(af) == 53 & ncol(af) == 74)
    expect_true(all(.path_cols_names() %in% colnames(af)))

    expect_message(
        {
            af <- import_association_file(fs_path$af,
                root = fs_path$root_inc,
                convert_tp = FALSE
            )
        },
        class = "summary_report"
    )
    expect_true(nrow(af) == 53 & ncol(af) == 74)
    expect_true(all(.path_cols_names() %in% colnames(af)))
})

### --- With dynamic vars
test_that("import_association_file custom - no fs align - no missing tags", {
    temp_specs <- .default_af_cols() |>
        dplyr::mutate(names = dplyr::if_else(
            condition = .data$names == "ProjectID",
            true = "Proj",
            false = .data$names
        ))
    new_af_path <- fs::path(fs::path_dir(fs_path$af), "asso_file_mod.tsv")
    readr::write_tsv(
        association_file |>
            dplyr::select(dplyr::all_of(.default_af_cols()$names)) |>
            dplyr::rename(Proj = "ProjectID"),
        file = new_af_path, na = ""
    )
    withr::local_options(.new = list(ISAnalytics.af_specs = temp_specs))
    expect_message(
        {
            af <- import_association_file(new_af_path,
                convert_tp = FALSE
            )
        },
        class = "summary_report"
    )
    expect_true(nrow(af) == 53 & ncol(af) == 71)
    expect_true("Proj" %in% colnames(af))
})

test_that("import_association_file custom - no fs align - missing tags", {
    temp_specs <- .default_af_cols() |> # missing project tag
        dplyr::filter(.data$tag != "project_id" | is.na(.data$tag))
    withr::local_options(.new = list(ISAnalytics.af_specs = temp_specs))
    expect_message(
        { # Everything should work, no alignment requested
            af <- import_association_file(fs_path$af,
                convert_tp = FALSE
            )
        },
        class = "summary_report"
    )
})

test_that("import_association_file custom - no fs align - missing cols", {
    new_af_path <- fs::path(fs::path_dir(fs_path$af), "asso_file_mod.tsv")
    readr::write_tsv(
        association_file |>
            dplyr::select(
                dplyr::all_of(.default_af_cols()$names),
                -dplyr::all_of("ProjectID")
            ),
        file = new_af_path, na = ""
    )
    expect_message(
        { ## Should work without errors since no fs alignment
            af <- import_association_file(new_af_path,
                convert_tp = FALSE
            )
        },
        class = "summary_report"
    )
})

test_that("import_association_file custom - fs align - no missing tags", {
    temp_specs <- .default_af_cols() |>
        dplyr::mutate(names = dplyr::if_else(
            condition = .data$names == "ProjectID",
            true = "Proj",
            false = .data$names
        ))
    new_af_path <- fs::path(fs::path_dir(fs_path$af), "asso_file_mod.tsv")
    readr::write_tsv(
        association_file |>
            dplyr::select(dplyr::all_of(.default_af_cols()$names)) |>
            dplyr::rename(Proj = "ProjectID"),
        file = new_af_path, na = ""
    )
    withr::local_options(.new = list(ISAnalytics.af_specs = temp_specs))
    expect_message(
        {
            af <- import_association_file(new_af_path,
                root = fs_path$root_corr,
                convert_tp = FALSE
            )
        },
        class = "summary_report"
    )
    expect_true(nrow(af) == 53 & ncol(af) == 74)
    expect_true("Proj" %in% colnames(af))
})

test_that("import_association_file custom - fs align - missing tags", {
    temp_specs <- .default_af_cols() |>
        dplyr::filter(.data$tag != "project_id") # missing project tag
    withr::local_options(.new = list(ISAnalytics.af_specs = temp_specs))
    expect_error(
        {
            af <- import_association_file(fs_path$af,
                root = fs_path$root_corr,
                convert_tp = FALSE
            )
        },
        class = "missing_tags_err"
    )
})

test_that("import_association_file custom - fs align - missing cols", {
    new_af_path <- fs::path(fs::path_dir(fs_path$af), "asso_file_mod.tsv")
    readr::write_tsv(
        association_file |>
            dplyr::select(
                dplyr::all_of(.default_af_cols()$names),
                -dplyr::all_of("ProjectID")
            ),
        file = new_af_path, na = ""
    )
    expect_error({ # Should error, missing ProjectID col
        af <- import_association_file(new_af_path,
            root = fs_path$root_corr,
            convert_tp = FALSE
        )
    })
})

test_that("import_association_file signals deprecation of tp_padding", {
    expect_message(
        {
            expect_deprecated({
                af <- import_association_file(fs_path$af,
                    convert_tp = FALSE,
                    tp_padding = 4
                )
            })
        },
        class = "summary_report"
    )
    expect_true(nrow(af) == 53 & ncol(af) == 71)
})

test_that("import_association_file works with import_iss = TRUE", {
    withr::local_options(list(ISAnalytics.verbose = FALSE))
    af <- import_association_file(fs_path$af,
        root = fs_path$root_corr,
        convert_tp = TRUE,
        import_iss = TRUE
    )
    expect_true(all(c(
        "RUN_NAME", "PHIX_MAPPING", "PLASMID_MAPPED_BYPOOL",
        "BARCODE_MUX", "LTR_IDENTIFIED", "TRIMMING_FINAL_LTRLC",
        "LV_MAPPED", "BWA_MAPPED_OVERALL", "ISS_MAPPED_OVERALL",
        "RAW_READS", "QUALITY_PASSED", "ISS_MAPPED_PP"
    ) %in%
        colnames(af)))
})
