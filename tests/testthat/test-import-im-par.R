library(ISAnalytics)
func_name <- c(
    ".interactive_select_projects_import", ".pool_number_IN",
    ".pool_choices_IN", ".interactive_select_pools_import", ".trace_anomalies",
    ".lookup_matrices", ".choose_duplicates_files_interactive",
    ".manage_anomalies_interactive", ".import_type", ".parallel_import_merge",
    "import_parallel_Vispa2Matrices_interactive",
    ".pattern_matching",
    ".update_as_option",
    ".lookup_matrices_auto", ".manage_anomalies_auto",
    "import_parallel_Vispa2Matrices_auto"
)

#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
withr::local_options(
    ISAnalytics.widgets = FALSE,
    ISAnalytics.verbose = FALSE
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

withr::with_options(list(
    ISAnalytics.widgets = FALSE,
    ISAnalytics.verbose = FALSE
), {
    as_file_correct <- import_association_file(path_af, root_correct)
})

withr::with_options(list(
    ISAnalytics.widgets = FALSE,
    ISAnalytics.verbose = FALSE
), {
    as_file_errors <- import_association_file(path_af, root_err)
})

all_proj_pools <- list(
    ProjectID = c("CLOEXP", "PROJECT1100", "PROJECT1101"),
    PoolID = c(
        "POOL6-1", "ABX-LR-PL5-POOL14-1",
        "ABX-LR-PL6-POOL15-1", "ABY-LR-PL4-POOL54-2"
    )
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
                .interactive_select_projects_import(as_file_correct)
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
        updated_af <- .interactive_select_projects_import(as_file_correct)
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
        dplyr::select(as_file_correct, .data$ProjectID)
    )$ProjectID
    invisible(capture_output(capture_output({
        updated_af <- .interactive_select_projects_import(as_file_correct)
    })))
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
    invisible(capture_output(capture_output({
        updated_af <- .interactive_select_projects_import(as_file_correct)
    })))
    selected_proj <- project_list[as.numeric(input_value2)]
    expect_true(all(updated_af$ProjectID %in% selected_proj))

    # Multiple choice
    input_value1 <- "2"
    input_value2 <- "1,2"
    input_value3 <- "y"
    write(input_value1, input1)
    write(input_value2, input2)
    write(input_value3, input3)
    invisible(capture_output(capture_output({
        updated_af <- .interactive_select_projects_import(as_file_correct)
    })))
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
    invisible(capture_output(capture_output({
        updated_af <- .interactive_select_projects_import(as_file_correct)
    })))
    indexes <- unlist(stringr::str_split(input_value2, ","))
    indexes <- as.numeric(indexes)
    selected_proj <- project_list[indexes]
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
                .interactive_select_projects_import(as_file_correct)
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
        updated_af <- .interactive_select_pools_import(as_file_correct)
    })))
    expect_true(all(all_proj_pools$PoolID %in%
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
        updated_af <- .interactive_select_pools_import(as_file_correct)
    })))
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

    files_found <- .lookup_matrices(
        as_file_correct, quant_types,
        "not_annotated"
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
    # Annotated matrices
    as_file_err <- as_file_errors %>% dplyr::filter(!is.na(.data$Path))
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
    # Annotated matrices
    as_file_err <- as_file_errors %>% dplyr::filter(!is.na(.data$Path))
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
    expect_message(
        {
            invisible(capture_output(capture_output({
                to_import <- .manage_anomalies_interactive(files_found_fs_ann)
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

    # Not annotated matrices
    expect_message(
        {
            invisible(capture_output(capture_output({
                to_import <- .manage_anomalies_interactive(
                    files_found_fs_notann
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
    as_file_err <- as_file_errors %>% dplyr::filter(!is.na(.data$Path))
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
    expect_message(
        {
            invisible(capture_output(capture_output({
                to_import <- .manage_anomalies_interactive(
                    files_found_fs_err_ann
                )
            })))
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
            invisible(capture_output(capture_output({
                to_import <- .manage_anomalies_interactive(
                    files_found_fs_err_notann
                )
            })))
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
test_that(paste(func_name[9], "correctly imports all files"), {
    as_file_err <- as_file_errors %>% dplyr::filter(!is.na(.data$Path))
    files_found_fs_ann <- .lookup_matrices(
        as_file_correct,
        quant_types, "annotated"
    )
    invisible(capture_output(capture_output({
        smpl_files_to_import <- .manage_anomalies_interactive(
            files_found_fs_ann
        )
    })))
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
test_that(paste(func_name[10], "correctly imports files for each q type"), {
    as_file_err <- as_file_errors %>% dplyr::filter(!is.na(.data$Path))
    files_found_fs_ann <- .lookup_matrices(
        as_file_correct,
        quant_types, "annotated"
    )
    invisible(capture_output(capture_output({
        smpl_files_to_import <- .manage_anomalies_interactive(
            files_found_fs_ann
        )
    })))
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
test_that(paste(func_name[11], "stops if af is missing"), {
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

test_that(paste(func_name[11], "stops if af is not character or tibble"), {
    expect_error({
        import_parallel_Vispa2Matrices_interactive(
            association_file = 1,
            quantification_type =
                quant_types
        )
    })
})

test_that(paste(
    func_name[11], "stops if af is a character",
    "vector longer than 1"
), {
    expect_error({
        import_parallel_Vispa2Matrices_interactive(
            association_file = c(
                path_af,
                path_af
            ),
            quantification_type =
                quant_types
        )
    })
})

test_that(paste(
    func_name[11], "stops if af is a character",
    "but root is incorrect"
), {
    expect_error({
        import_parallel_Vispa2Matrices_interactive(
            association_file = path_af,
            quantification_type =
                quant_types
        )
    })
    expect_error({
        import_parallel_Vispa2Matrices_interactive(
            association_file = path_af,
            root = 1,
            quantification_type =
                quant_types
        )
    })
    expect_error({
        import_parallel_Vispa2Matrices_interactive(
            association_file = path_af,
            root = c(root_correct, root_err),
            quantification_type =
                quant_types
        )
    })
})

test_that(paste(
    func_name[11], "stops if workers is not",
    "numeric or have multiple values"
), {
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


test_that(paste(
    func_name[11], "stops if quantification_types",
    "is incorrect or missing"
), {
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

test_that(paste(func_name[11], "stops if matrix_type is incorrect"), {
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
test_that(paste(func_name[11], "succeeds for fs - path"), {
    # Annotated
    input1 <- withr::local_tempfile()
    input3 <- withr::local_tempfile()
    input4 <- withr::local_tempfile()
    input6 <- withr::local_tempfile()
    op1 <- withr::local_options(ISAnalytics.connection = c(
        input1, NA, input3, input4, NA,
        input6, NA, NA
    ), ISAnalytics.widgets = FALSE, ISAnalytics.verbose = FALSE)
    input_value1 <- "1"
    input_value3 <- "y"
    input_value4 <- "1"
    input_value6 <- "y"
    write(input_value1, input1)
    write(input_value3, input3)
    write(input_value4, input4)
    write(input_value6, input6)
    invisible(capture_output(capture_output({
        matrices <-
            import_parallel_Vispa2Matrices_interactive(
                association_file = path_af,
                root = root_correct,
                quantification_type =
                    quant_types,
                matrix_type = "annotated",
                workers = 2,
                dates_format = "dmy",
                multi_quant_matrix = FALSE
            )
    })))
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
    invisible(capture_output(capture_output({
        matrices <-
            import_parallel_Vispa2Matrices_interactive(path_af,
                root = root_correct,
                quantification_type =
                    quant_types,
                matrix_type = "not_annotated",
                workers = 2,
                dates_format = "dmy",
                multi_quant_matrix = FALSE
            )
    })))
    expect_equal(nrow(matrices$fragmentEstimate), 4 * 1476)
    expect_equal(nrow(matrices$seqCount), 4 * 1476)
})

test_that(paste(func_name[11], "succeeds for fs - tibble"), {
    # Annotated
    input1 <- withr::local_tempfile()
    input3 <- withr::local_tempfile()
    input4 <- withr::local_tempfile()
    input6 <- withr::local_tempfile()
    op1 <- withr::local_options(ISAnalytics.connection = c(
        input1, NA, input3, input4, NA,
        input6, NA, NA
    ), ISAnalytics.widgets = FALSE, ISAnalytics.verbose = FALSE)
    input_value1 <- "1"
    input_value3 <- "y"
    input_value4 <- "1"
    input_value6 <- "y"
    write(input_value1, input1)
    write(input_value3, input3)
    write(input_value4, input4)
    write(input_value6, input6)
    invisible(capture_output(capture_output({
        matrices <-
            import_parallel_Vispa2Matrices_interactive(as_file_correct,
                quantification_type =
                    quant_types,
                matrix_type = "annotated",
                workers = 2,
                multi_quant_matrix = FALSE
            )
    })))
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
    invisible(capture_output(capture_output({
        matrices <-
            import_parallel_Vispa2Matrices_interactive(as_file_correct,
                NULL,
                quantification_type =
                    quant_types,
                matrix_type = "not_annotated",
                workers = 2,
                multi_quant_matrix = FALSE
            )
    })))
    expect_equal(nrow(matrices$fragmentEstimate), 4 * 1476)
    expect_equal(nrow(matrices$seqCount), 4 * 1476)
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
test_that(paste(func_name[13], "works for fs err - single pattern"), {
    as_file_err <- as_file_errors %>% dplyr::filter(!is.na(.data$Path))
    files_found_fs_ann <- .lookup_matrices(
        as_file_correct,
        quant_types, "annotated"
    )
    files_found_fs_err_ann <- .lookup_matrices(
        as_file_err,
        quant_types, "annotated"
    )
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

test_that(paste(func_name[13], "works for fs err - multiple patterns"), {
    as_file_err <- as_file_errors %>% dplyr::filter(!is.na(.data$Path))
    files_found_fs_ann <- .lookup_matrices(
        as_file_correct,
        quant_types, "annotated"
    )
    files_found_fs_err_ann <- .lookup_matrices(
        as_file_err,
        quant_types, "annotated"
    )
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
test_that(paste(func_name[14], "returns correctly with patterns null"), {
    to_import <- .lookup_matrices_auto(as_file_correct, quant_types,
        "annotated",
        patterns = NULL, matching_opt = "ANY"
    )
    expect_true(all(to_import$Anomalies == FALSE))
    to_import <- .lookup_matrices_auto(as_file_correct, quant_types,
        "annotated",
        patterns = NULL, matching_opt = "ALL"
    )
    expect_true(all(to_import$Anomalies == FALSE))
    to_import <- .lookup_matrices_auto(as_file_correct, quant_types,
        "annotated",
        patterns = NULL, matching_opt = "OPTIONAL"
    )
    expect_true(all(to_import$Anomalies == FALSE))
})

test_that(paste(func_name[14], "returns correctly for fs"), {
    to_import <- .lookup_matrices_auto(as_file_correct, quant_types,
        "annotated",
        patterns = "NoMate", matching_opt = "ANY"
    )
    expect_true(all(to_import$Anomalies == TRUE))
    to_import <- .lookup_matrices_auto(as_file_correct, quant_types,
        "annotated",
        patterns = "NoMate", matching_opt = "ALL"
    )
    expect_true(all(to_import$Anomalies == TRUE))
    to_import <- .lookup_matrices_auto(as_file_correct, quant_types,
        "annotated",
        patterns = "NoMate",
        matching_opt = "OPTIONAL"
    )
    expect_true(all(to_import$Anomalies == FALSE))
})

test_that(paste(func_name[14], "returns correctly for fserr"), {
    as_file_err <- as_file_errors %>% dplyr::filter(!is.na(.data$Path))
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
test_that(paste(func_name[15], "correctly manages for fs"), {
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

test_that(paste(func_name[15], "correctly manages for fserr"), {
    as_file_err <- as_file_errors %>% dplyr::filter(!is.na(.data$Path))
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
test_that(paste(func_name[16], "stops if af is missing"), {
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

test_that(paste(
    func_name[16], "stops if af is not",
    "character or tibble"
), {
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

test_that(paste(
    func_name[16], "stops if af is a character",
    "vector longer than 1"
), {
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

test_that(paste(func_name[16], "stops if patterns is incorrect"), {
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

test_that(paste(func_name[16], "stops if matching_option is incorrect"), {
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
test_that(paste(func_name[16], "succeeds for fs - path"), {
    matrices <- import_parallel_Vispa2Matrices_auto(path_af,
        root = root_correct,
        quantification_type = quant_types,
        matrix_type = "annotated",
        workers = 2,
        patterns = NULL,
        matching_opt = "ANY",
        dates_format = "dmy",
        multi_quant_matrix = FALSE
    )
    expect_true(nrow(matrices$fragmentEstimate) == 4 * 1476)
    expect_true(nrow(matrices$seqCount) == 4 * 1476)
})

test_that(paste(func_name[16], "succeeds for fs - tibble"), {
    af <- as_file_correct %>% dplyr::filter(.data$ProjectID == "CLOEXP")
    matrices <- import_parallel_Vispa2Matrices_auto(af,
        root = NULL,
        quantification_type = quant_types,
        matrix_type = "annotated",
        workers = 2,
        patterns = NULL,
        matching_opt = "ANY",
        multi_quant_matrix = FALSE
    )
    expect_true(nrow(matrices$fragmentEstimate) == 1476)
    expect_true(nrow(matrices$seqCount) == 1476)
})

test_that(paste(func_name[16], "succeeds for fserr - tibble"), {
    as_file_err <- as_file_errors %>% dplyr::filter(!is.na(.data$Path))
    af <- as_file_err %>% dplyr::filter(.data$ProjectID == "CLOEXP")
    matrices <- import_parallel_Vispa2Matrices_auto(af,
        quantification_type = quant_types,
        matrix_type = "annotated",
        workers = 2,
        patterns = "NoMate",
        matching_opt = "ANY",
        multi_quant_matrix = FALSE
    )
    expect_true(nrow(matrices$fragmentEstimate) == 1476)
    expect_true(nrow(matrices$seqCount) == 1476)
    matrices <- import_parallel_Vispa2Matrices_auto(af,
        quantification_type = quant_types,
        matrix_type = "annotated",
        workers = 2,
        patterns = "NoMate",
        matching_opt = "ALL",
        multi_quant_matrix = FALSE
    )
    expect_true(nrow(matrices$fragmentEstimate) == 1476)
    expect_true(nrow(matrices$seqCount) == 1476)
    matrices <- import_parallel_Vispa2Matrices_auto(af,
        quantification_type = quant_types,
        matrix_type = "annotated",
        workers = 2,
        patterns = "NoMate",
        matching_opt = "OPTIONAL",
        multi_quant_matrix = FALSE
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
            quantification_type = quant_types,
            matrix_type = "annotated",
            workers = 2,
            patterns = "NoMate",
            matching_opt = "ANY",
            multi_quant_matrix = FALSE
        )
    })
    expect_error({
        matrices <- import_parallel_Vispa2Matrices_auto(af,
            Nquantification_type = quant_types,
            matrix_type = "annotated",
            workers = 2,
            patterns = "NoMate",
            matching_opt = "ALL",
            multi_quant_matrix = FALSE
        )
    })
    matrices <- import_parallel_Vispa2Matrices_auto(af,
        quantification_type = quant_types,
        matrix_type = "annotated",
        workers = 2,
        patterns = "NoMate",
        matching_opt = "OPTIONAL",
        multi_quant_matrix = FALSE
    )
    expect_true(nrow(matrices$seqCount) == 1476)
    af <- as_file_err %>%
        dplyr::filter(
            .data$ProjectID == "PROJECT1100",
            .data$concatenatePoolIDSeqRun == "ABX-LR-PL6-POOL15-1"
        )
    expect_error({
        matrices <- import_parallel_Vispa2Matrices_auto(af,
            quantification_type = quant_types,
            matrix_type = "annotated",
            workers = 2,
            patterns = "NoMate",
            matching_opt = "ANY",
            multi_quant_matrix = FALSE
        )
    })
    expect_error({
        matrices <- import_parallel_Vispa2Matrices_auto(af,
            quantification_type = quant_types,
            matrix_type = "annotated",
            workers = 2,
            patterns = "NoMate",
            matching_opt = "ALL",
            multi_quant_matrix = FALSE
        )
    })
    matrices <- import_parallel_Vispa2Matrices_auto(af,
        quantification_type = quant_types,
        matrix_type = "annotated",
        workers = 2,
        patterns = "NoMate",
        matching_opt = "OPTIONAL",
        multi_quant_matrix = FALSE
    )
    expect_true(nrow(matrices$fragmentEstimate) == 1476)
    expect_true(nrow(matrices$seqCount) == 1476)
})
