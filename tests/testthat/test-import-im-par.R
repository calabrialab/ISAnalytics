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
to_dupl_case1 <- fs::path(
    copy_path, "PJ01", "quantification", "POOL01-1",
    paste0(
        "PJ01_POOL01-1_fragmentEstimate_",
        "matrix.no0.annotated.tsv.gz"
    )
)
to_dupl_case2 <- fs::path(
    copy_path, "PJ01", "quantification", "POOL02-1",
    paste0(
        "PJ01_POOL02-1_fragmentEstimate_",
        "matrix.no0.annotated.tsv.gz"
    )
)
fs::file_copy(
    to_dupl_case1,
    fs::path(
        copy_path, "PJ01", "quantification", "POOL01-1",
        paste0(
            "PJ01_POOL01-1_FOO_fragmentEstimate_",
            "matrix.no0.annotated.tsv.gz"
        )
    ),
    overwrite = TRUE
)
fs::file_copy(
    to_dupl_case2,
    fs::path(
        copy_path, "PJ01", "quantification", "POOL02-1",
        paste0(
            "PJ01_POOL02-1_fragmentEstimate_",
            "matrix.tsv.gz"
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
    tb <- tb |> tidyr::nest(Files = c("Quantification_type", "Files_found"))
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
# Tests .trace_anomalies
#------------------------------------------------------------------------------#
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
test_that(".lookup_matrices finds all files for fs - no duplicates", {
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

test_that(".lookup_matrices finds all files for fserr - duplicates", {
    aligned_af_mod <- local_af_corr |>
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
    expect_true(files_found |>
        dplyr::filter(.data$concatenatePoolIDSeqRun == "POOL01-1") |>
        dplyr::pull(.data$Anomalies))
    expect_true((files_found |>
        dplyr::filter(.data$concatenatePoolIDSeqRun == "POOL01-1") |>
        dplyr::pull(.data$Files_count))[[1]] |>
        dplyr::filter(.data$Quantification_type == "fragmentEstimate") |>
        dplyr::pull(.data$Found) == 2)
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
test_that(".update_as_option works as expected", {
    # ---- PREP
    aligned_af_mod <- local_af_corr |>
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
    pool_with_dupl <- (files_found_fs_dupl |>
        dplyr::filter(.data$concatenatePoolIDSeqRun == "POOL01-1") |>
        dplyr::pull(.data$Files))[[1]]

    # ---- TEST SINGLE PATTERN
    pattern <- c("FOO")
    pattern_match <- .pattern_matching(pool_with_dupl$Files_found, pattern)
    expect_true(pattern_match$FOO[1] == TRUE &
        all(pattern_match$FOO[c(2, 3)] == FALSE))
    updated <- .update_as_option(pool_with_dupl, pattern_match, "ANY")
    expect_true(nrow(updated) == 2)
    expect_true(all(quant_types %in% updated$Quantification_type))
    expect_true(is.na(updated |>
        dplyr::filter(.data$Quantification_type == "seqCount") |>
        dplyr::pull(.data$Files_found)))
    updated <- .update_as_option(pool_with_dupl, pattern_match, "ALL")
    expect_true(nrow(updated) == 2)
    expect_true(all(quant_types %in% updated$Quantification_type))
    expect_true(is.na(updated |>
        dplyr::filter(.data$Quantification_type == "seqCount") |>
        dplyr::pull(.data$Files_found)))
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
test_that(".lookup_matrices_auto returns correctly with patterns null", {
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

test_that(".lookup_matrices_auto returns correctly for fs", {
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

test_that(".lookup_matrices_auto returns correctly for fserr", {
    aligned_af_mod <- local_af_corr |>
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
            is.na((to_import |>
                dplyr::filter(
                    .data$concatenatePoolIDSeqRun == .x
                )
            )$Files[[1]]$Files_found)
        )
    })
    expect_true(all(na_files))
    non_na_files <- !is.na((to_import |>
        dplyr::filter(
            .data$concatenatePoolIDSeqRun == "POOL01-1"
        ))$Files[[1]] |>
        dplyr::filter(.data$Quantification_type == "fragmentEstimate") |>
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
            is.na((to_import |>
                dplyr::filter(
                    .data$concatenatePoolIDSeqRun == .x
                )
            )$Files[[1]]$Files_found)
        )
    })
    expect_true(all(na_files))
    non_na_files <- !is.na((to_import |>
        dplyr::filter(
            .data$concatenatePoolIDSeqRun == "POOL01-1"
        ))$Files[[1]] |>
        dplyr::filter(.data$Quantification_type == "fragmentEstimate") |>
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
test_that(".manage_anomalies_auto correctly manages for fs", {
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

test_that(".manage_anomalies_auto correctly manages for fserr", {
    withr::local_options(list(ISAnalytics.verbose = TRUE))
    aligned_af_mod <- local_af_corr |>
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
test_that(".import_type works as expected", {
    sample_df_1 <- tibble::tribble(
        ~chr, ~integration_locus, ~strand, ~GeneName, ~GeneStrand, ~id1, ~id2,
        ~id3, ~id4,
        "1", 140546, "+", "GENE1", "-", 4, NA, NA, 1,
        "14", 43567, "-", "GENE2", "+", 231, NA, 2, NA,
        "5", 214676, "-", "GENE3", "-", NA, NA, NA, NA,
        "7", 66778, "-", "GENE4", "-", NA, 355, NA, NA,
        "1", 75687, "+", "GENE5", "+", NA, NA, NA, 65,
        "5", 64576, "+", "GENE6", "-", 1, 667, NA, NA,
        "20", 57587, "-", "GENE7", "-", NA, 13, 1, NA,
        "X", 457658, "+", "GENE8", "+", NA, NA, NA, 768
    )
    sample_df_2 <- tibble::tribble(
        ~chr, ~integration_locus, ~strand, ~GeneName, ~GeneStrand, ~id5, ~id6,
        ~id7, ~id8,
        "1", 324167, "+", "GENE9", "-", 4, NA, NA, 1,
        "14", 86223, "-", "GENE10", "+", 231, NA, 2, NA,
        "5", 96423, "-", "GENE11", "-", NA, NA, NA, NA,
        "7", 34211, "-", "GENE12", "-", NA, 355, NA, NA,
        "1", 75687, "+", "GENE5", "+", NA, NA, NA, 65,
        "5", 562123, "+", "GENE13", "-", 1, 667, NA, NA,
        "20", 452134, "-", "GENE14", "-", NA, 13, 1, NA,
        "X", 457658, "+", "GENE8", "+", NA, NA, NA, 768
    )
    temp_path_1 <- withr::local_tempfile(fileext = ".tsv.gz")
    readr::write_tsv(sample_df_1, file = temp_path_1, na = "")
    temp_path_2 <- withr::local_tempfile(fileext = ".tsv.gz")
    readr::write_tsv(sample_df_2, file = temp_path_2, na = "")
    sample_file_to_import <- tibble::tribble(
        ~ProjectID, ~concatenatePoolIDSeqRun, ~Quantification_type,
        ~Files_chosen,
        "P1", "POOL1", "seqCount", temp_path_1,
        "P1", "POOL2", "seqCount", temp_path_2
    )
    output <- .import_type(
        type_df = sample_file_to_import,
        import_matrix_args = list(
            call_mode = "INTERNAL",
            id_col_name = "Sample"
        ),
        progr = NULL, max_workers = 2
    )
    expect_true(is.data.frame(output$matrix) & nrow(output$matrix) == 22)
    expect_true(is.data.frame(output$imported_files))
    expect_true(all(output$imported_files$Imported))
    expect_true(all(output$imported_files$Number_of_samples == 4))
    expect_true(all(output$imported_files$Distinct_is == 7))
})

test_that(".import_type reports errors silently", {
    sample_df_1 <- tibble::tribble(
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
    sample_df_2 <- tibble::tribble(
        ~chr, ~integration_locus, ~strand, ~GeneName, ~GeneStrand, ~id5, ~id6,
        ~id7, ~id8,
        "1", 324167, "+", "GENE9", "-", 4, NA, NA, 1,
        "14", 86223, "-", "GENE10", "+", 231, NA, 2, NA,
        "5", 96423, "-", "GENE11", "-", NA, NA, NA, NA,
        "7", 34211, "-", "GENE12", "-", NA, 355, NA, NA,
        "1", 75687, "+", "GENE5", "+", NA, NA, NA, 65,
        "5", 562123, "+", "GENE13", "-", 1, 667, NA, NA,
        "20", 452134, "-", "GENE14", "-", NA, 13, 1, NA,
        "X", 457658, "+", "GENE8", "+", NA, NA, NA, 768
    )
    temp_path_1 <- withr::local_tempfile(fileext = ".tsv.gz")
    readr::write_tsv(sample_df_1, file = temp_path_1, na = "")
    temp_path_2 <- withr::local_tempfile(fileext = ".tsv.gz")
    readr::write_tsv(sample_df_2, file = temp_path_2, na = "")
    sample_file_to_import <- tibble::tribble(
        ~ProjectID, ~concatenatePoolIDSeqRun, ~Quantification_type,
        ~Files_chosen,
        "P1", "POOL1", "seqCount", temp_path_1,
        "P1", "POOL2", "seqCount", temp_path_2
    )
    output <- .import_type(
        type_df = sample_file_to_import,
        import_matrix_args = list(
            call_mode = "INTERNAL",
            id_col_name = "Sample"
        ),
        progr = NULL, max_workers = 2
    )
    expect_true(is.data.frame(output$matrix) & nrow(output$matrix) == 11)
    expect_true(is.data.frame(output$imported_files))
    expect_false(output$imported_files$Imported[1])
    expect_true(output$imported_files$Imported[2])
})

#------------------------------------------------------------------------------#
# Tests .parallel_import_merge
#------------------------------------------------------------------------------#
test_that(".parallel_import_merge works as expected", {
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
    temp_path <- withr::local_tempfile(fileext = ".tsv.gz")
    readr::write_tsv(sample_df, file = temp_path, na = "")
    sample_df_corr <- sample_df |>
        tidyr::separate(.data$IS_GenomicID,
            into = c("chr", "integration_locus", "strand"),
            sep = "_", remove = TRUE, convert = TRUE
        )
    temp_path_corr <- withr::local_tempfile(fileext = ".tsv.gz")
    readr::write_tsv(sample_df_corr, file = temp_path_corr, na = "")
    sample_file_to_import <- tibble::tribble(
        ~ProjectID, ~concatenatePoolIDSeqRun, ~Quantification_type,
        ~Files_chosen,
        "P1", "POOL1", "seqCount", temp_path,
        "P1", "POOL1", "fragmentEstimate", temp_path_corr
    )
    results <- .parallel_import_merge(sample_file_to_import,
        workers = 2,
        import_matrix_args = list(
            call_mode = "INTERNAL",
            id_col_name = "Sample"
        )
    )
    expect_true(is.null(results$matrix$seqCount) &
        !is.null(results$matrix$fragmentEstimate))
    expect_equal(results$summary |>
        dplyr::filter(
            .data$Quantification_type == "fragmentEstimate"
        ) |>
        dplyr::pull(.data$Imported), TRUE)
    expect_equal(results$summary |>
        dplyr::filter(
            .data$Quantification_type == "seqCount"
        ) |>
        dplyr::pull(.data$Imported), FALSE)
})

#------------------------------------------------------------------------------#
# Tests .base_param_check
#------------------------------------------------------------------------------#
test_that(".base_param_check executes correctly", {
    expect_error(
        {
            .base_param_check(
                local_af_corr, c("seqCount", "fragmentEstimate"),
                "annotated", 2, TRUE
            )
        },
        regexp = NA
    )
    expect_error({
        .base_param_check(
            "abc", c("seqCount", "fragmentEstimate"),
            "annotated", 2, TRUE
        )
    })
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
test_that(".pre_manage_af executes correctly", {
    # Already imported af
    checked_af <- .pre_manage_af(local_af_corr)
    expect_equal(local_af_corr, checked_af)
})

#------------------------------------------------------------------------------#
# Tests import_parallel_Vispa2Matrices
#------------------------------------------------------------------------------#
test_that("import_parallel_Vispa2Matrices executes correctly", {
    matrices <- import_parallel_Vispa2Matrices(
        association_file = local_af_corr,
        quantification_type = c("seqCount", "fragmentEstimate"),
        mode = "AUTO", report_path = NULL
    )
    expect_true(nrow(matrices) == 1689 & ncol(matrices) == 8)
    expect_true(all(c("seqCount", "fragmentEstimate") %in% colnames(matrices)))
})

test_that("deprecation of old functions gets signaled", {
    expect_defunct({
        matrices <- import_parallel_Vispa2Matrices_auto(
            local_af_corr,
            quantification_type = c("seqCount", "fragmentEstimate")
        )
    })
    expect_defunct({
        matrices <- import_parallel_Vispa2Matrices_interactive(
            local_af_corr,
            quantification_type = c("seqCount", "fragmentEstimate")
        )
    })
})

test_that("import_parallel_Vispa2Matrices signals deprecation of interactive", {
    expect_error({
        matrices <- import_parallel_Vispa2Matrices(
            association_file = local_af_corr,
            quantification_type = c("seqCount", "fragmentEstimate"),
            mode = "INTERACTIVE", report_path = NULL
        )
    })
})
