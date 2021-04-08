library(ISAnalytics)

#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
op <- withr::local_options(
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
root_err <- system.file("extdata", "fserr.zip",
    package = "ISAnalytics"
)

root_err <- unzip_file_system(root_err, "fserr")

# Association file
association_file <- import_association_file(path_af, root_correct,
    dates_format = "dmy"
)
suppressWarnings({
    association_file_err <- import_association_file(path_af, root_err,
        dates_format = "dmy"
    )
})

# Matrices
matrices <- import_parallel_Vispa2Matrices_auto(
    association_file = path_af, root = root_correct,
    quantification_type = c("seqCount", "fragmentEstimate"),
    matrix_type = "annotated", workers = 2, patterns = NULL,
    matching_opt = "ANY",
    dates_format = "dmy", multi_quant_matrix = FALSE
)

#------------------------------------------------------------------------------#
# Tests .stats_report
#------------------------------------------------------------------------------#
test_that(".stats_report returns NULL if all paths are NA", {
    association_file_mod <- association_file
    association_file_mod$Path <- rep_len(
        NA_character_,
        length(association_file_mod$Path)
    )
    stats_rep <- .stats_report(association_file_mod)
    expect_null(stats_rep)
    # Or filtering af
    af_filtered <- association_file_err %>%
        dplyr::filter(.data$ProjectID == "PROJECT1101")
    stats_rep <- .stats_report(af_filtered)
    expect_null(stats_rep)
})

test_that(".stats_report returns correctly for both af", {
    stats_rep <- .stats_report(association_file)
    expect_true(all(!is.na(stats_rep$stats_files)))
    stats_rep <- .stats_report(association_file_err)
    expect_true(any(is.na(stats_rep$stats_path)))
})

#------------------------------------------------------------------------------#
# Tests .import_stats_iss
#------------------------------------------------------------------------------#
test_that(".import_stats_iss returns null if nothing to import", {
    af_filtered <- association_file_err %>%
        dplyr::filter(.data$ProjectID == "PROJECT1101")
    stats <- .import_stats_iss(af_filtered)
    expect_null(stats)
})

test_that(".import_stats_iss detects missing/malformed files", {
    stats <- .import_stats_iss(association_file_err)
    report <- stats[[2]]
    expect_true(all((report %>%
        dplyr::filter(.data$ProjectID == "PROJECT1100")
    )$Imported == FALSE))
    expect_true(all((report %>%
        dplyr::filter(.data$ProjectID == "CLOEXP")
    )$Imported == TRUE))
})

test_that(".import_stats_iss returns NULL if no files were imported", {
    af_filtered <- association_file_err %>%
        dplyr::filter(.data$ProjectID == "PROJECT1100")
    stats <- .import_stats_iss(af_filtered)
    expect_null(stats)
})

#------------------------------------------------------------------------------#
# Tests .join_and_aggregate
#------------------------------------------------------------------------------#
test_that(".join_and_aggregate correctly join and aggregate for stats", {
    af_filtered <- association_file %>%
        dplyr::filter(.data$ProjectID == "CLOEXP")
    stats <- .import_stats_iss(af_filtered)
    agg <- .join_and_aggregate(
        af_filtered, stats[[1]],
        c(
            "SubjectID", "CellMarker",
            "Tissue", "TimePoint"
        )
    )
    expect_equal(agg$VCN_avg, c(0.7, 2.8, 4.3, 9.89), ignore_attr = TRUE)
    expect_true(all(agg$DNAngUsed_avg == 100))
    expect_equal(agg$Kapa_avg, c(26.54667, 59.81000, 79.80333, 85.87000),
        tolerance = .5, ignore_attr = TRUE
    )
    expect_true(all(agg$DNAngUsed_sum == 300))
    expect_equal(agg$ulForPool_sum, c(25.7, 12.29, 7.98, 7.46),
        ignore_attr = TRUE
    )
    expect_equal(agg$BARCODE_MUX_sum, c(4872448, 8573674, 5841068, 6086489),
        ignore_attr = TRUE
    )
    expect_equal(
        agg$TRIMMING_FINAL_LTRLC_sum,
        c(4863932, 8549992, 5818598, 6062299),
        ignore_attr = TRUE
    )
    expect_equal(
        agg$LV_MAPPED_sum,
        c(2114395, 3485464, 2923561, 2484034),
        ignore_attr = TRUE
    )
    expect_equal(
        agg$BWA_MAPPED_OVERALL_sum,
        c(2619004, 4751887, 2631197, 3297007),
        ignore_attr = TRUE
    )
    expect_equal(
        agg$ISS_MAPPED_PP_sum,
        c(2386173, 3746655, 1938140, 2330286),
        ignore_attr = TRUE
    )
})

test_that(".join_and_aggregate correctly join and aggregate for af only", {
    af_filtered <- association_file %>%
        dplyr::filter(.data$ProjectID == "CLOEXP")
    agg <- .join_and_aggregate(
        af_filtered, NULL,
        c(
            "SubjectID", "CellMarker",
            "Tissue", "TimePoint"
        )
    )
    expect_equal(agg$VCN_avg, c(0.7, 2.8, 4.3, 9.89), ignore_attr = TRUE)
    expect_true(all(agg$DNAngUsed_avg == 100))
    expect_equal(agg$Kapa_avg, c(26.54667, 59.81000, 79.80333, 85.87000),
        tolerance = .5, ignore_attr = TRUE
    )
    expect_true(all(agg$DNAngUsed_sum == 300))
    expect_equal(agg$ulForPool_sum, c(25.7, 12.29, 7.98, 7.46),
        ignore_attr = TRUE
    )
    expect_true(all(!.stats_columns_min() %in% colnames(agg)))
})

#------------------------------------------------------------------------------#
# Tests aggregate_metadata
#------------------------------------------------------------------------------#
# Test input
test_that("aggregate_metadata stops if af is missing mandatory columns", {
    association_file <- association_file %>%
        dplyr::select(-c(.data$FusionPrimerPCRDate))
    expect_error({
        agg_meta <- aggregate_metadata(association_file = association_file)
    })
})
test_that("aggregate_metadata stops if grouping keys is null", {
    expect_error({
        agg_meta <- aggregate_metadata(
            association_file = association_file,
            grouping_keys = NULL
        )
    })
})
test_that("aggregate_metadata stops if grouping keys is not a char vector", {
    expect_error({
        agg_meta <- aggregate_metadata(
            association_file = association_file,
            grouping_keys = 1
        )
    })
    expect_error({
        agg_meta <- aggregate_metadata(
            association_file = association_file,
            grouping_keys = c(1, 2)
        )
    })
})
test_that("aggregate_metadata stops if grouping keys are missing from af", {
    expect_error({
        agg_meta <- aggregate_metadata(
            association_file = association_file,
            grouping_keys = "a"
        )
    })
    expect_error({
        agg_meta <- aggregate_metadata(
            association_file = association_file,
            grouping_keys = c("a", "ProjectID")
        )
    })
})

test_that("aggregate_metadata stops if import_stats is not logical", {
    expect_error({
        agg_meta <- aggregate_metadata(
            association_file = association_file,
            import_stats = 1
        )
    })
    expect_error({
        agg_meta <- aggregate_metadata(
            association_file = association_file,
            import_stats = c(TRUE, TRUE, FALSE)
        )
    })
})

#------------------------------------------------------------------------------#
# Tests aggregate_values_by_key
#------------------------------------------------------------------------------#

# Test input
test_that("aggregate_values_by_key stops if x incorrect", {
    expect_error({
        agg <- aggregate_values_by_key(
            x = 1,
            association_file = association_file
        )
    })
    expect_error(
        {
            agg <- aggregate_values_by_key(
                x = tibble::tibble(a = c(1, 2, 3)),
                association_file = association_file
            )
        },
        regexp = .non_ISM_error()
    )
    expect_error({
        agg <- aggregate_values_by_key(matrices$seqCount, association_file,
            value_cols = "x"
        )
    })
    expect_error({
        agg <- aggregate_values_by_key(matrices$seqCount, association_file,
            value_cols = c("chr")
        )
    })
    expect_error({
        agg <- aggregate_values_by_key(
            matrices$seqCount %>%
                dplyr::select(-c("CompleteAmplificationID")),
            association_file
        )
    })
    expect_error({
        agg <- aggregate_values_by_key(matrices, association_file,
            value_cols = "x"
        )
    })
    expect_error({
        agg <- aggregate_values_by_key(matrices, association_file,
            value_cols = c("chr")
        )
    })
    mod_list <- purrr::map(matrices, function(m) {
        m %>% dplyr::select(-c("CompleteAmplificationID"))
    })
    expect_error(
        {
            agg <- aggregate_values_by_key(
                mod_list,
                association_file
            )
        },
        regexp = .missing_complAmpID_error()
    )
})

test_that("aggregate_values_by_key stops if key is incorrect", {
    expect_error({
        agg <- aggregate_values_by_key(
            x = matrices$seqCount,
            association_file = association_file,
            key = 1
        )
    })
    expect_error(
        {
            agg <- aggregate_values_by_key(
                x = matrices$seqCount,
                association_file = association_file,
                key = "x"
            )
        },
        regexp = "Key fields are missing from association file"
    )
    expect_error(
        {
            agg <- aggregate_values_by_key(
                x = matrices$seqCount,
                association_file = association_file,
                key = c("SubjectID", "x")
            )
        },
        regexp = "Key fields are missing from association file"
    )
})

test_that("aggregate_values_by_key stops if lambda is incorrect", {
    expect_error({
        agg <- aggregate_values_by_key(
            x = isadf,
            association_file = association_file,
            key = "SubjectID",
            lambda = 1
        )
    })
    expect_error({
        agg <- aggregate_values_by_key(
            x = isadf,
            association_file = association_file,
            key = "SubjectID",
            lambda = sum
        )
    })
    expect_error({
        agg <- aggregate_values_by_key(
            x = isadf,
            association_file = association_file,
            key = "SubjectID",
            lambda = c("sum", "mean")
        )
    })
})

test_that("aggregate_values_by_key stops if group is incorrect", {
    expect_error({
        agg <- aggregate_values_by_key(
            x = matrices$seqCount,
            association_file = association_file,
            group = c(1, 2)
        )
    })
    expect_error(
        {
            agg <- aggregate_values_by_key(
                x = matrices$seqCount,
                association_file = association_file,
                group = c("xfsf", "fwre")
            )
        },
        regexp = "Grouping variables not found"
    )
})

# Test Values
test_that("aggregate_values_by_key succeeds for single IS matrix", {
    expect_error(
        {
            agg <- aggregate_values_by_key(matrices$seqCount, association_file)
        },
        regexp = NA
    )
    expect_true(all(c(
        "chr", "integration_locus", "strand", "GeneName",
        "GeneStrand", "SubjectID",
        "Value_sum"
    ) %in% colnames(agg)))
    comp <- comparison_matrix(matrices)
    agg <- aggregate_values_by_key(comp, association_file,
        value_cols = c(
            "fragmentEstimate",
            "seqCount"
        )
    )
    expect_true(all(c(
        "chr", "integration_locus", "strand", "GeneName",
        "GeneStrand", "SubjectID",
        "fragmentEstimate_sum", "seqCount_sum"
    ) %in% colnames(agg)))
})

test_that("aggregate_values_by_key succeeds for list of IS matrices", {
    expect_error(
        {
            agg <- aggregate_values_by_key(matrices, association_file)
        },
        regexp = NA
    )
    expect_true(all(c(
        "chr", "integration_locus", "strand", "GeneName",
        "GeneStrand", "SubjectID",
        "Value_sum"
    ) %in% colnames(agg[[1]])))
    expect_true(all(c(
        "chr", "integration_locus", "strand", "GeneName",
        "GeneStrand", "SubjectID",
        "Value_sum"
    ) %in% colnames(agg[[2]])))
})
