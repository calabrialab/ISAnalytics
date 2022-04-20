library(ISAnalytics)
func_name <- c(
    ".interactive_select_projects_import",
    ".pool_number_IN",
    ".pool_choices_IN",
    ".interactive_select_pools_import",
    ".trace_anomalies",
    ".lookup_matrices",
    ".choose_duplicates_files_interactive",
    ".manage_anomalies_interactive",
    ".import_type",
    ".parallel_import_merge",
    "import_parallel_Vispa2Matrices_interactive",
    ".pattern_matching",
    ".update_as_option",
    ".lookup_matrices_auto",
    ".manage_anomalies_auto",
    ".base_param_check",
    ".pre_manage_af",
    "import_parallel_Vispa2Matrices"
)

#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
withr::local_options(
    ISAnalytics.reports = FALSE,
    ISAnalytics.verbose = FALSE
)

all_proj_pools <- list(
    ProjectID = unique(association_file$ProjectID),
    PoolID = unique(association_file$PoolID)
)

copy_path <- fs::path(tempdir(), "fs_dupl")
## Generate file system structure containing duplicates
fs::dir_copy(fs_path$root_corr, copy_path, overwrite = TRUE)
to_dupl <- fs::path(
    copy_path, "PJ01", "quantification", "POOL01-1",
    paste0(
        "PJ01_POOL01-1_fragmentEstimate_",
        "matrix.no0.annotated.tsv.gz"
    )
)
fs::file_copy(
    to_dupl,
    fs::path(
        copy_path, "PJ01", "quantification", "POOL01-1",
        paste0(
            "PJ01_POOL01-1_FOO_fragmentEstimate_",
            "matrix.no0.annotated.tsv.gz"
        )
    ),
    overwrite = TRUE
)

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

quant_types <- c("fragmentEstimate", "seqCount")

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

#------------------------------------------------------------------------------#
# Tests .interactive_select_projects_import
#------------------------------------------------------------------------------#
test_that(paste(func_name[1], "stops if input 1 is 0"), {
    input1 <- withr::local_tempfile()
    input_value1 <- "0"
    withr::local_options(ISAnalytics.connection = input1)
    write(input_value1, input1)
    expect_error(
        {
            invisible(capture_output(capture_output({
                .interactive_select_projects_import(local_af_corr, "ProjectID")
            })))
        },
        regexp = "Quitting"
    )
})

test_that(paste(func_name[1], "selects all projects when inputs are '1, y'"), {
    input1 <- withr::local_tempfile()
    input2 <- withr::local_tempfile()
    input3 <- withr::local_tempfile()
    input_value1 <- "1"
    input_value3 <- "y"
    op <- withr::local_options(ISAnalytics.connection = c(
        input1,
        input2,
        input3
    ))
    write(input_value1, input1)
    write(input_value3, input3)
    invisible(capture_output(capture_output({
        updated_af <- .interactive_select_projects_import(
            local_af_corr,
            "ProjectID"
        )
    })))
    expect_true(all(all_proj_pools$ProjectID %in% updated_af$ProjectID))
})

test_that(paste(
    func_name[1], "selects only some projects",
    "when input 1 is '2'"
), {
    # Single choice
    input1 <- withr::local_tempfile()
    input2 <- withr::local_tempfile()
    input3 <- withr::local_tempfile()
    input_value1 <- "2"
    input_value2 <- "1"
    input_value3 <- "y"
    op <- withr::local_options(ISAnalytics.connection = c(
        input1,
        input2,
        input3
    ))
    write(input_value1, input1)
    write(input_value2, input2)
    write(input_value3, input3)
    project_list <- dplyr::distinct(
        dplyr::select(local_af_corr, .data$ProjectID)
    )$ProjectID
    invisible(capture_output(capture_output({
        updated_af <- .interactive_select_projects_import(
            local_af_corr,
            "ProjectID"
        )
    })))
    selected_proj <- project_list[as.numeric(input_value2)]
    expect_true(all(updated_af$ProjectID %in% selected_proj))
})

test_that(paste(func_name[1], "stops if input 2 is 0"), {
    input1 <- withr::local_tempfile()
    input2 <- withr::local_tempfile()
    input_value1 <- "2"
    input_value2 <- "0"
    op <- withr::local_options(ISAnalytics.connection = c(input1, input2))
    write(input_value1, input1)
    write(input_value2, input2)
    expect_error(
        {
            invisible(capture_output(capture_output({
                .interactive_select_projects_import(
                    local_af_corr,
                    "ProjectID"
                )
            })))
        },
        regexp = "Quitting"
    )
})

#------------------------------------------------------------------------------#
# Tests .pool_number_IN
#------------------------------------------------------------------------------#
test_that(paste(func_name[2], "stops when input is 0"), {
    input4 <- withr::local_tempfile()
    op <- withr::local_options(ISAnalytics.connection = c(NA, NA, NA, input4))
    input_value4 <- "0"
    write(input_value4, input4)
    expect_error(
        {
            invisible(capture_output(capture_output({
                .pool_number_IN()
            })))
        },
        regexp = "Quitting"
    )
})

test_that(paste(func_name[2], "returns user choice correctly when not 0"), {
    input4 <- withr::local_tempfile()
    op <- withr::local_options(ISAnalytics.connection = c(NA, NA, NA, input4))
    input_value4 <- "1"
    write(input_value4, input4)
    invisible(capture_output(capture_output({
        choice <- .pool_number_IN()
    })))
    expect_equal(choice, 1)

    input_value4 <- "2"
    write(input_value4, input4)
    invisible(capture_output(capture_output({
        choice <- .pool_number_IN()
    })))
    expect_equal(choice, 2)
})

#------------------------------------------------------------------------------#
# Tests .pool_choices_IN
#------------------------------------------------------------------------------#
test_that(paste(func_name[3], "stops when input is 0"), {
    input5 <- withr::local_tempfile()
    op <- withr::local_options(ISAnalytics.connection = c(
        NA, NA, NA, NA,
        input5
    ))
    input_value5 <- "0"
    write(input_value5, input5)
    expect_error(
        {
            invisible(capture_output(capture_output({
                .pool_choices_IN(c(1, 2, 3))
            })))
        },
        regexp = "Quitting"
    )
})

test_that(paste(
    func_name[3], "returns user choice correctly when not",
    "0 - single"
), {
    input5 <- withr::local_tempfile()
    op <- withr::local_options(ISAnalytics.connection = c(
        NA, NA, NA, NA,
        input5
    ))
    input_value5 <- "1"
    write(input_value5, input5)
    invisible(capture_output(capture_output({
        choice <- .pool_choices_IN(c(1, 2, 3))
    })))
    expect_true(choice == 1)

    input_value5 <- "2"
    write(input_value5, input5)
    invisible(capture_output(capture_output({
        choice <- .pool_choices_IN(c(1, 2, 3))
    })))
    expect_true(choice == 2)

    input_value5 <- "3"
    write(input_value5, input5)
    invisible(capture_output(capture_output({
        choice <- .pool_choices_IN(c(1, 2, 3))
    })))
    expect_true(choice == 3)
})

test_that(paste(
    func_name[3], "returns user choice correctly",
    "when not 0 - multi"
), {
    input5 <- withr::local_tempfile()
    op <- withr::local_options(ISAnalytics.connection = c(
        NA, NA, NA, NA,
        input5
    ))
    input_value5 <- "1,2"
    write(input_value5, input5)
    invisible(capture_output(capture_output({
        choice <- .pool_choices_IN(c(1, 2, 3))
    })))
    expect_equal(choice, c(1, 2))

    input_value5 <- "2,3"
    write(input_value5, input5)
    invisible(capture_output(capture_output({
        choice <- .pool_choices_IN(c(1, 2, 3))
    })))
    expect_equal(choice, c(2, 3))
})

#------------------------------------------------------------------------------#
# Tests .interactive_select_pools_import
#------------------------------------------------------------------------------#
test_that(paste(
    func_name[4], "selects all pools if input is 1"
), {
    input4 <- withr::local_tempfile()
    input6 <- withr::local_tempfile()
    op <- withr::local_options(ISAnalytics.connection = c(
        NA, NA, NA,
        input4, NA, input6
    ))
    input_value4 <- "1"
    input_value6 <- "y"
    write(input_value4, input4)
    write(input_value6, input6)
    invisible(capture_output(capture_output({
        updated_af <- .interactive_select_pools_import(
            local_af_corr,
            "ProjectID",
            "concatenatePoolIDSeqRun"
        )
    })))
    expect_true(all(all_proj_pools$concatenatePoolIDSeqRun %in%
        updated_af$concatenatePoolIDSeqRun))
})

test_that(paste(func_name[4], "selects only the first pools
          if input is '2, 1, y'"), {
    input4 <- withr::local_tempfile()
    input5 <- withr::local_tempfile()
    input6 <- withr::local_tempfile()
    op <- withr::local_options(ISAnalytics.connection = c(
        NA, NA, NA, input4,
        input5, input6
    ))
    input_value4 <- "2"
    input_value5 <- "1"
    input_value6 <- "y"
    write(input_value4, input4)
    write(input_value5, input5)
    write(input_value6, input6)
    invisible(capture_output(capture_output({
        updated_af <- .interactive_select_pools_import(
            local_af_corr,
            "ProjectID",
            "concatenatePoolIDSeqRun"
        )
    })))
    available <- local_af_corr %>%
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
test_that(paste(func_name[5], "updates tibble correctly"), {
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
test_that(paste(func_name[6], "finds all files for fs - no duplicates"), {
    files_found <- .lookup_matrices(
        local_af_corr,
        quant_types,
        "annotated",
        "ProjectID",
        "concatenatePoolIDSeqRun"
    )
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

test_that(paste(func_name[6], "finds all files for fserr - duplicates"), {
    aligned_af_mod <- local_af_corr %>%
        dplyr::mutate(!!.path_cols_names()$quant :=
            stringr::str_replace(
                .data[[.path_cols_names()$quant]],
                "/fs/", "/fs_dupl/"
            ))
    files_found <- .lookup_matrices(
        aligned_af_mod,
        quant_types,
        "annotated",
        "ProjectID",
        "concatenatePoolIDSeqRun"
    )
    expect_true(files_found %>%
        dplyr::filter(.data$concatenatePoolIDSeqRun == "POOL01-1") %>%
        dplyr::pull(.data$Anomalies))
    expect_true((files_found %>%
        dplyr::filter(.data$concatenatePoolIDSeqRun == "POOL01-1") %>%
        dplyr::pull(.data$Files_count))[[1]] %>%
        dplyr::filter(.data$Quantification_type == "fragmentEstimate") %>%
        dplyr::pull(.data$Found) == 2)
})

#------------------------------------------------------------------------------#
# Tests .choose_duplicates_files_interactive
#------------------------------------------------------------------------------#
test_that(paste(func_name[7], "modifies tbl according to user choices"), {
    input7 <- withr::local_tempfile()
    op <- withr::local_options(ISAnalytics.connection = c(
        NA, NA, NA, NA, NA, NA,
        input7
    ))
    input_value7 <- "1"
    write(input_value7, input7)
    smpl_dupl_t <- smpl_dupl()
    dupl_quant <- "fragmentEstimate"
    invisible(capture_output(capture_output({
        updated_tbl <- .choose_duplicates_files_interactive(
            dupl_quant,
            smpl_dupl_t
        )
    })))
    expect_true(all(c("fragmentEstimate", "seqCount") %in%
        updated_tbl$Quantification_type))
    filtered <- updated_tbl %>%
        dplyr::filter(.data$Quantification_type == dupl_quant)
    expect_equal(nrow(filtered), 1)
    expect_equal(filtered$Files_found, c("pth1"))
})

test_that(paste(func_name[7], "stops if input is 0"), {
    input7 <- withr::local_tempfile()
    op <- withr::local_options(ISAnalytics.connection = c(
        NA, NA, NA, NA, NA, NA,
        input7
    ))
    input_value7 <- "0"
    write(input_value7, input7)
    smpl_dupl_t <- smpl_dupl()
    dupl_quant <- "fragmentEstimate"
    expect_error(
        {
            invisible(capture_output(capture_output({
                updated_tbl <- .choose_duplicates_files_interactive(
                    dupl_quant,
                    smpl_dupl_t
                )
            })))
        },
        regexp = "Quitting"
    )
})

#------------------------------------------------------------------------------#
# Tests .manage_anomalies_interactive
#------------------------------------------------------------------------------#
test_that(paste(func_name[8], "succeeds with no anomalies"), {
    files_found_fs_ann <- .lookup_matrices(
        local_af_corr,
        quant_types,
        "annotated",
        "ProjectID",
        "concatenatePoolIDSeqRun"
    )
    files_found_fs_err_ann <- .lookup_matrices(
        local_af_inc %>%
            dplyr::filter(!is.na(.data[[.path_cols_names()$quant]])),
        quant_types,
        "annotated",
        "ProjectID",
        "concatenatePoolIDSeqRun"
    )

    expect_message(
        {
            invisible(capture_output(capture_output({
                to_import <- .manage_anomalies_interactive(
                    files_found_fs_ann,
                    "ProjectID",
                    "concatenatePoolIDSeqRun"
                )
            })))
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

test_that(paste(func_name[8], "succeeds with anomalies"), {
    aligned_af_mod <- local_af_corr %>%
        dplyr::mutate(!!.path_cols_names()$quant :=
            stringr::str_replace(
                .data[[.path_cols_names()$quant]],
                "/fs/", "/fs_dupl/"
            ))
    files_found <- .lookup_matrices(
        aligned_af_mod,
        quant_types,
        "annotated",
        "ProjectID",
        "concatenatePoolIDSeqRun"
    )
    input7 <- withr::local_tempfile()
    input8 <- withr::local_tempfile()
    op <- withr::local_options(ISAnalytics.connection = c(
        NA, NA, NA, NA, NA, NA, input7,
        input8
    ))
    input_value7 <- "1"
    input_value8 <- "y"
    write(input_value7, input7)
    write(input_value8, input8)
    # Annotated matrices
    expect_message({
        invisible(capture_output(capture_output({
            to_import <- .manage_anomalies_interactive(
                files_found,
                "ProjectID",
                "concatenatePoolIDSeqRun"
            )
        })))
    })
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
    expect_true(to_import %>%
        dplyr::filter(
            .data$concatenatePoolIDSeqRun == "POOL01-1",
            .data$Quantification_type == "fragmentEstimate"
        ) %>%
        dplyr::pull(.data$Files_chosen) == ((files_found %>%
        dplyr::filter(.data$concatenatePoolIDSeqRun ==
            "POOL01-1"))$Files[[1]] %>%
        dplyr::filter(.data$Quantification_type == "fragmentEstimate") %>%
        dplyr::pull(.data$Files_found))[1])
})

#------------------------------------------------------------------------------#
# Tests .pattern_matching
#------------------------------------------------------------------------------#
test_that(paste(func_name[12], "matches single pattern"), {
    p_matches <- .pattern_matching(c("file1", "file2", "file3test"), "file")
    expect_true(all(colnames(p_matches) %in% c("file")))
    expect_true(all(p_matches$file == TRUE))
})

test_that(paste(func_name[12], "matches multiple patterns"), {
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
test_that(paste(func_name[13], "works as expected"), {
    # ---- PREP
    aligned_af_mod <- local_af_corr %>%
        dplyr::mutate(!!.path_cols_names()$quant :=
            stringr::str_replace(
                .data[[.path_cols_names()$quant]],
                "/fs/", "/fs_dupl/"
            ))
    files_found_fs_dupl <- .lookup_matrices(
        aligned_af_mod,
        quant_types, "annotated",
        "ProjectID",
        "concatenatePoolIDSeqRun"
    )
    pool_with_dupl <- (files_found_fs_dupl %>%
        dplyr::filter(.data$concatenatePoolIDSeqRun == "POOL01-1") %>%
        dplyr::pull(.data$Files))[[1]]

    # ---- TEST SINGLE PATTERN
    pattern <- c("FOO")
    pattern_match <- .pattern_matching(pool_with_dupl$Files_found, pattern)
    expect_true(pattern_match$FOO[1] == TRUE &
        all(pattern_match$FOO[c(2, 3)] == FALSE))
    updated <- .update_as_option(pool_with_dupl, pattern_match, "ANY")
    expect_true(nrow(updated) == 2)
    expect_true(all(quant_types %in% updated$Quantification_type))
    expect_true(is.na(updated %>%
        dplyr::filter(.data$Quantification_type == "seqCount") %>%
        dplyr::pull(.data$Files_found)))
    updated <- .update_as_option(pool_with_dupl, pattern_match, "ALL")
    expect_true(nrow(updated) == 2)
    expect_true(all(quant_types %in% updated$Quantification_type))
    expect_true(is.na(updated %>%
        dplyr::filter(.data$Quantification_type == "seqCount") %>%
        dplyr::pull(.data$Files_found)))
    updated <- .update_as_option(pool_with_dupl, pattern_match, "OPTIONAL")
    expect_true(nrow(updated) == 2)
    expect_true(all(quant_types %in% updated$Quantification_type))
    expect_true(all(!is.na(updated$Files_found)))

    # ---- TEST MULTIPLE PATTERN
    pattern <- c("FOO", "BAR", "X")
    pattern_match <- .pattern_matching(pool_with_dupl$Files_found, pattern)
    updated <- .update_as_option(pool_with_dupl, pattern_match, "ANY")
    expect_true(nrow(updated) == 2)
    expect_true(all(quant_types %in% updated$Quantification_type))
    expect_true(is.na(updated %>%
        dplyr::filter(.data$Quantification_type == "seqCount") %>%
        dplyr::pull(.data$Files_found)))
    updated <- .update_as_option(pool_with_dupl, pattern_match, "ALL")
    expect_true(nrow(updated) == 2)
    expect_true(all(quant_types %in% updated$Quantification_type))
    expect_true(all(is.na(updated$Files_found)))
    updated <- .update_as_option(pool_with_dupl, pattern_match, "OPTIONAL")
    expect_true(nrow(updated) == 2)
    expect_true(all(quant_types %in% updated$Quantification_type))
    expect_true(all(!is.na(updated$Files_found)))

    # ---- TEST NEGATED PATTERNS
    # * Keep all the files that DO NOT contain "FOO"
    pattern <- c("^(?!.*FOO).*")
    pattern_match <- .pattern_matching(pool_with_dupl$Files_found, pattern)
    updated <- .update_as_option(pool_with_dupl, pattern_match, "ANY")
    expect_true(nrow(updated) == 2)
    expect_true(all(quant_types %in% updated$Quantification_type))
    expect_true(all(!is.na(updated$Files_found)))
    updated <- .update_as_option(pool_with_dupl, pattern_match, "ALL")
    expect_true(nrow(updated) == 2)
    expect_true(all(quant_types %in% updated$Quantification_type))
    expect_true(all(!is.na(updated$Files_found)))
    updated <- .update_as_option(pool_with_dupl, pattern_match, "OPTIONAL")
    expect_true(nrow(updated) == 2)
    expect_true(all(quant_types %in% updated$Quantification_type))
    expect_true(all(!is.na(updated$Files_found)))
})

#------------------------------------------------------------------------------#
# Tests .lookup_matrices_auto
#------------------------------------------------------------------------------#
test_that(paste(func_name[14], "returns correctly with patterns null"), {
    to_import <- .lookup_matrices_auto(local_af_corr,
        quant_types,
        "annotated",
        patterns = NULL, matching_opt = "ANY",
        "ProjectID",
        "concatenatePoolIDSeqRun"
    )
    expect_true(all(to_import$Anomalies == FALSE))
    to_import <- .lookup_matrices_auto(local_af_corr, quant_types,
        "annotated",
        patterns = NULL, matching_opt = "ALL",
        "ProjectID",
        "concatenatePoolIDSeqRun"
    )
    expect_true(all(to_import$Anomalies == FALSE))
    to_import <- .lookup_matrices_auto(local_af_corr, quant_types,
        "annotated",
        patterns = NULL, matching_opt = "OPTIONAL",
        "ProjectID",
        "concatenatePoolIDSeqRun"
    )
    expect_true(all(to_import$Anomalies == FALSE))
})

test_that(paste(func_name[14], "returns correctly for fs"), {
    to_import <- .lookup_matrices_auto(local_af_corr, quant_types,
        "annotated",
        patterns = "FOO", matching_opt = "ANY",
        "ProjectID",
        "concatenatePoolIDSeqRun"
    )
    expect_true(all(to_import$Anomalies == TRUE))
    to_import <- .lookup_matrices_auto(local_af_corr, quant_types,
        "annotated",
        patterns = "FOO", matching_opt = "ALL",
        "ProjectID",
        "concatenatePoolIDSeqRun"
    )
    expect_true(all(to_import$Anomalies == TRUE))
    to_import <- .lookup_matrices_auto(local_af_corr, quant_types,
        "annotated",
        patterns = "FOO",
        matching_opt = "OPTIONAL",
        "ProjectID",
        "concatenatePoolIDSeqRun"
    )
    expect_true(all(to_import$Anomalies == FALSE))
})

test_that(paste(func_name[14], "returns correctly for fserr"), {
    aligned_af_mod <- local_af_corr %>%
        dplyr::mutate(!!.path_cols_names()$quant :=
            stringr::str_replace(
                .data[[.path_cols_names()$quant]],
                "/fs/", "/fs_dupl/"
            ))
    to_import <- .lookup_matrices_auto(aligned_af_mod,
        quant_types, "annotated",
        patterns = "FOO", matching_opt = "ANY",
        "ProjectID",
        "concatenatePoolIDSeqRun"
    )
    expect_true(all(to_import$Anomalies == TRUE))
    na_pools <- c(
        "POOL02-1",
        "POOL03-1",
        "POOL04-1"
    )
    na_files <- purrr::map_lgl(na_pools, ~ {
        all(
            is.na((to_import %>%
                dplyr::filter(
                    .data$concatenatePoolIDSeqRun == .x
                )
            )$Files[[1]]$Files_found)
        )
    })
    expect_true(all(na_files))
    non_na_files <- !is.na((to_import %>%
        dplyr::filter(
            .data$concatenatePoolIDSeqRun == "POOL01-1"
        ))$Files[[1]] %>%
        dplyr::filter(.data$Quantification_type == "fragmentEstimate") %>%
        dplyr::pull(.data$Files_found))
    expect_true(non_na_files)

    to_import <- .lookup_matrices_auto(aligned_af_mod,
        quant_types, "annotated",
        patterns = "FOO", matching_opt = "ALL",
        "ProjectID",
        "concatenatePoolIDSeqRun"
    )
    expect_true(all(to_import$Anomalies == TRUE))
    na_files <- purrr::map_lgl(na_pools, ~ {
        all(
            is.na((to_import %>%
                dplyr::filter(
                    .data$concatenatePoolIDSeqRun == .x
                )
            )$Files[[1]]$Files_found)
        )
    })
    expect_true(all(na_files))
    non_na_files <- !is.na((to_import %>%
        dplyr::filter(
            .data$concatenatePoolIDSeqRun == "POOL01-1"
        ))$Files[[1]] %>%
        dplyr::filter(.data$Quantification_type == "fragmentEstimate") %>%
        dplyr::pull(.data$Files_found))
    expect_true(non_na_files)

    to_import <- .lookup_matrices_auto(aligned_af_mod,
        quant_types, "annotated",
        patterns = "FOO", matching_opt = "OPTIONAL",
        "ProjectID",
        "concatenatePoolIDSeqRun"
    )
    expect_true(all(to_import$Anomalies == FALSE))
})

#------------------------------------------------------------------------------#
# Tests .manage_anomalies_auto
#------------------------------------------------------------------------------#
test_that(paste(func_name[15], "correctly manages for fs"), {
    files_found <- .lookup_matrices_auto(local_af_corr,
        quant_types, "annotated",
        patterns = NULL, matching_opt = "ANY",
        "ProjectID",
        "concatenatePoolIDSeqRun"
    )
    files_to_import <- .manage_anomalies_auto(
        files_found,
        "ProjectID",
        "concatenatePoolIDSeqRun"
    )
    expect_true(nrow(files_to_import) == 4 * length(quant_types))
})

test_that(paste(func_name[15], "correctly manages for fserr"), {
    withr::local_options(list(ISAnalytics.verbose = TRUE))
    aligned_af_mod <- local_af_corr %>%
        dplyr::mutate(!!.path_cols_names()$quant :=
            stringr::str_replace(
                .data[[.path_cols_names()$quant]],
                "/fs/", "/fs_dupl/"
            ))
    files_found <- .lookup_matrices_auto(aligned_af_mod,
        quant_types, "annotated",
        patterns = NULL, matching_opt = "ANY",
        "ProjectID",
        "concatenatePoolIDSeqRun"
    )
    expect_message(
        {
            files_to_import <- .manage_anomalies_auto(
                files_found,
                "ProjectID",
                "concatenatePoolIDSeqRun"
            )
        },
        class = "auto_mode_dupl"
    )
    expect_true(nrow(files_to_import) == 7)
    files_found <- .lookup_matrices_auto(aligned_af_mod,
        quant_types, "annotated",
        patterns = "FOO",
        matching_opt = "ANY",
        "ProjectID",
        "concatenatePoolIDSeqRun"
    )
    expect_message(
        {
            files_to_import <- .manage_anomalies_auto(
                files_found,
                "ProjectID",
                "concatenatePoolIDSeqRun"
            )
        },
        class = "auto_mode_miss"
    )
    expect_true(nrow(files_to_import) == 1)
})

#------------------------------------------------------------------------------#
# Tests .import_type
#------------------------------------------------------------------------------#
test_that(paste(func_name[9], "reports errors silently"), {
    # Expecting errors when mand vars do not mach
    sample_df <- tibble::tribble(
        ~IS_GenomicID, ~GeneName, ~GeneStrand, ~id1, ~id2,
        ~id3, ~id4,
        "1_140546_+", "GENE1", "-", 4, NA, NA, 1,
        "14_43567_-", "GENE2", "+", 231, NA, 2, NA,
        "5_214676_-", "GENE3", "-", NA, NA, NA, NA,
        "7_66778_-", "GENE4", "-", NA, 355, NA, NA,
        "1_75687_+", "GENE5", "+", NA, NA, NA, 65,
        "5_64576_+", "GENE6", "-", 1, 667, NA, NA,
        "20_57587_-", "GENE7", "-", NA, 13, 1, NA,
        "X_457658_+", "GENE8", "+", NA, NA, NA, 768
    )
    temp_path <- withr::local_tempfile(fileext = "tsv.gz")
    readr::write_tsv(sample_df, file = temp_path, na = "")
    sample_file_to_import <- tibble::tribble(
        ~ProjectID, ~concatenatePoolIDSeqRun, ~Quantification_type,
        ~Files_chosen,
        "P1", "POOL1", "seqCount", temp_path
    )
    p <- if (.Platform$OS.type == "windows") {
        BiocParallel::SnowParam(
            workers = 2,
            stop.on.error = FALSE
        )
    } else {
        BiocParallel::MulticoreParam(
            workers = 2,
            stop.on.error = FALSE
        )
    }
    import_result <- .import_type(
        q_type = "seqCount",
        files = sample_file_to_import,
        cluster = p, import_matrix_args = list()
    )
    expect_true(is.null(import_result$matrix))
    expect_true(import_result$imported_files$Imported[1] == FALSE)
    new_mand <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "IS_GenomicID", "char", NULL, "required", NA_character_
    )
    withr::with_options(list(ISAnalytics.mandatory_is_vars = new_mand), {
        import_result <- .import_type(
            q_type = "seqCount",
            files = sample_file_to_import,
            cluster = p, import_matrix_args = list()
        )
        expect_true(!is.null(import_result$matrix))
        expect_true(nrow(import_result$matrix) == 11 &
            ncol(import_result$matrix) == 5)
        expect_true(import_result$imported_files$Imported[1] == TRUE)
        expect_true(import_result$imported_files$Number_of_samples[1] == 4)
        expect_true(import_result$imported_files$Distinct_is[1] == 7)
    })
    BiocParallel::bpstop(p)
})

#------------------------------------------------------------------------------#
# Tests .parallel_import_merge
#------------------------------------------------------------------------------#
test_that(paste(func_name[10], "works as expected"), {
    sample_df <- tibble::tribble(
        ~IS_GenomicID, ~GeneName, ~GeneStrand, ~id1, ~id2,
        ~id3, ~id4,
        "1_140546_+", "GENE1", "-", 4, NA, NA, 1,
        "14_43567_-", "GENE2", "+", 231, NA, 2, NA,
        "5_214676_-", "GENE3", "-", NA, NA, NA, NA,
        "7_66778_-", "GENE4", "-", NA, 355, NA, NA,
        "1_75687_+", "GENE5", "+", NA, NA, NA, 65,
        "5_64576_+", "GENE6", "-", 1, 667, NA, NA,
        "20_57587_-", "GENE7", "-", NA, 13, 1, NA,
        "X_457658_+", "GENE8", "+", NA, NA, NA, 768
    )
    temp_path <- withr::local_tempfile(fileext = "tsv.gz")
    readr::write_tsv(sample_df, file = temp_path, na = "")
    sample_df_corr <- sample_df %>%
        tidyr::separate(.data$IS_GenomicID,
            into = c("chr", "integration_locus", "strand"),
            sep = "_", remove = TRUE, convert = TRUE
        )
    temp_path_corr <- withr::local_tempfile(fileext = "tsv.gz")
    readr::write_tsv(sample_df_corr, file = temp_path_corr, na = "")
    sample_file_to_import <- tibble::tribble(
        ~ProjectID, ~concatenatePoolIDSeqRun, ~Quantification_type,
        ~Files_chosen,
        "P1", "POOL1", "seqCount", temp_path,
        "P1", "POOL1", "fragmentEstimate", temp_path_corr
    )
    results <- .parallel_import_merge(sample_file_to_import,
        workers = 2,
        import_matrix_args = list()
    )
    expect_true(is.null(results$matrix$seqCount) &
        !is.null(results$matrix$fragmentEstimate))
    expect_equal(results$summary$Imported, c(FALSE, TRUE))
})

#------------------------------------------------------------------------------#
# Tests .base_param_check
#------------------------------------------------------------------------------#
test_that(paste(func_name[16], "executes correctly"), {
    expect_error(
        {
            .base_param_check(
                local_af_corr, c("seqCount", "fragmentEstimate"),
                "annotated", 2, TRUE
            )
        },
        regexp = NA
    )
    expect_error(
        {
            .base_param_check(
                "abc", c("seqCount", "fragmentEstimate"),
                "annotated", 2, TRUE
            )
        },
        regexp = NA
    )
    expect_error({
        .base_param_check(
            local_af_corr, c("a", "b"),
            "annotated", 2, TRUE
        )
    })
    expect_error(
        {
            .base_param_check(
                local_af_corr, c("seqCount", "fragmentEstimate"),
                "a", 2, TRUE
            )
        },
        regexp = NA
    )
    expect_error({
        .base_param_check(
            local_af_corr, c("seqCount", "fragmentEstimate"),
            "annotated", "blabla", TRUE
        )
    })
    expect_error({
        .base_param_check(
            local_af_corr, c("seqCount", "fragmentEstimate"),
            "annotated", 2, "TRUE"
        )
    })
})

#------------------------------------------------------------------------------#
# Tests .pre_manage_af
#------------------------------------------------------------------------------#
test_that(paste(func_name[17], "executes correctly"), {
    # Already imported af
    checked_af <- .pre_manage_af(local_af_corr, list(), report_path = NULL)
    expect_equal(local_af_corr, checked_af)

    # To import
    checked_af <- .pre_manage_af(fs_path$af, list(root = fs_path$root_corr),
        report_path = NULL
    )
    expect_equal(local_af_corr, checked_af)
})

test_that(paste(func_name[17], "produces report in right location"), {
    withr::local_options(list(ISAnalytics.reports = TRUE))
    # When provided only a folder
    temp_fold <- tempdir()
    complete_filename <- .clean_file_path(temp_fold, "asso_file")
    expect_message({
        checked_af <- .pre_manage_af(fs_path$af, list(root = fs_path$root_corr),
            report_path = temp_fold
        )
    })
    expect_equal(local_af_corr, checked_af)
    expect_true(fs::file_exists(complete_filename))

    # When provided folder and file for matrix import
    matrix_report_filname <- fs::path(temp_fold, "matrix_report.html")
    expect_message({
        checked_af <- .pre_manage_af(fs_path$af, list(root = fs_path$root_corr),
            report_path = matrix_report_filname
        )
    })
    expect_equal(local_af_corr, checked_af)
    expect_true(fs::file_exists(complete_filename))
})

#------------------------------------------------------------------------------#
# Tests import_parallel_Vispa2Matrices
#------------------------------------------------------------------------------#
test_that(paste(func_name[18], "executes correctly"), {
    matrices <- import_parallel_Vispa2Matrices(
        association_file = local_af_corr,
        quantification_type = c("seqCount", "fragmentEstimate"),
        mode = "AUTO", report_path = NULL
    )
    expect_true(nrow(matrices) == 1689 & ncol(matrices) == 8)
    expect_true(all(c("seqCount", "fragmentEstimate") %in% colnames(matrices)))
})

test_that("deprecation of old functions gets signaled - AUTO", {
    expect_deprecated({
        matrices <- import_parallel_Vispa2Matrices_auto(
            local_af_corr,
            quantification_type = c("seqCount", "fragmentEstimate")
        )
    })
    matrices_new_fun <- import_parallel_Vispa2Matrices(
        association_file = local_af_corr,
        quantification_type = c("seqCount", "fragmentEstimate"),
        mode = "AUTO", report_path = NULL
    )
    expect_equal(matrices %>%
        dplyr::arrange(.data$seqCount, .data$fragmentEstimate),
    matrices_new_fun %>%
        dplyr::arrange(.data$seqCount, .data$fragmentEstimate),
    ignore_attr = TRUE
    )
})

test_that("deprecation of old functions gets signaled - INTERACTIVE", {
    input1 <- withr::local_tempfile()
    input2 <- withr::local_tempfile()
    input3 <- withr::local_tempfile()
    input4 <- withr::local_tempfile()
    input5 <- withr::local_tempfile()
    input6 <- withr::local_tempfile()
    input_value1 <- "1"
    input_value3 <- "y"
    input_value4 <- "1"
    input_value6 <- "y"
    op <- withr::local_options(ISAnalytics.connection = c(
        input1,
        input2,
        input3,
        input4,
        input5,
        input6
    ))
    write(input_value1, input1)
    write(input_value3, input3)
    write(input_value4, input4)
    write(input_value6, input6)
    expect_deprecated({
        invisible(capture.output({
            matrices <- import_parallel_Vispa2Matrices_interactive(
                local_af_corr,
                quantification_type = c("seqCount", "fragmentEstimate")
            )
        }))
    })
    invisible(capture.output({
        matrices_new_fun <- import_parallel_Vispa2Matrices(
            association_file = local_af_corr,
            quantification_type = c("seqCount", "fragmentEstimate"),
            mode = "INTERACTIVE", report_path = NULL
        )
    }))
    expect_equal(matrices %>%
        dplyr::arrange(.data$seqCount, .data$fragmentEstimate),
    matrices_new_fun %>%
        dplyr::arrange(.data$seqCount, .data$fragmentEstimate),
    ignore_attr = TRUE
    )
})
