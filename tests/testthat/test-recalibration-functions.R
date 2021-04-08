library(ISAnalytics)

#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
op <- withr::local_options(
    ISAnalytics.widgets = FALSE,
    ISAnalytics.verbose = FALSE
)
# Samples
sample_group1 <- tibble::tibble(
    chr = c(rep_len("1", 6)),
    integration_locus = c(
        14572, 14572,
        14575, 14571,
        14577, 14581
    ),
    strand = c(rep_len("+", 6)),
    CompleteAmplificationID = paste0("ID", 1:6),
    Value = c(70, 150, 120, 1400, 36, 15)
)

sample_group2 <- tibble::tibble(
    chr = rep_len("3", 7),
    integration_locus = c(
        16380,
        16396,
        16402,
        16395,
        16378,
        16399,
        16387
    ),
    strand = rep_len("+", 7),
    CompleteAmplificationID = c(
        "ID1", "ID4",
        "ID5", "ID4",
        "ID1", "ID3",
        "ID1"
    ),
    Value = c(
        1846, 64, 543, 89, 123, 886,
        48
    )
)
sample_group_mult1 <- tibble::tibble(
    chr = c(rep_len("1", 6)),
    integration_locus = c(
        14572, 14572,
        14575, 14571,
        14577, 14581
    ),
    strand = c(rep_len("+", 6)),
    CompleteAmplificationID = paste0("ID", 1:6),
    seqCount = c(70, 150, 120, 1400, 36, 15),
    fragmentEstimate = c(
        83.4, 125.9, 1.656,
        64.545, 6.4, 564.69
    )
)

sample_group_mult2 <- tibble::tibble(
    chr = rep_len("3", 7),
    integration_locus = c(
        16380,
        16396,
        16402,
        16395,
        16378,
        16399,
        16387
    ),
    strand = rep_len("+", 7),
    CompleteAmplificationID = c(
        "ID1", "ID4",
        "ID5", "ID4",
        "ID1", "ID3",
        "ID1"
    ),
    seqCount = c(
        1846, 64, 543, 89, 123, 886,
        48
    ),
    fragmentEstimate = c(
        334.54, 5456.45, 12.55,
        64.65, 654.5, 453.6, 1.36
    )
)

#------------------------------------------------------------------------------#
# Tests .find_unique_max
#------------------------------------------------------------------------------#
test_that(".find_unique_max returns empty if values incomparable", {
    # Max is not unique, with NAs
    t1 <- c(NA, NA, 170, 170)
    # All NAs
    t2 <- c(NA_real_, NA_real_, NA_real_)
    # The max is not unique
    t3 <- c(50, 60, 70, 70)
    max <- .find_unique_max(t1)
    expect_equal(max, numeric(0))
    max <- .find_unique_max(t2)
    expect_equal(max, numeric(0))
    max <- .find_unique_max(t3)
    expect_equal(max, numeric(0))
})

test_that(".find_unique_max returns value if values comparable", {
    # Max is unique, with NAs
    t1 <- c(NA_real_, NA_real_, 170)
    # The max is unique
    t2 <- c(50, 50, 60, 70)
    max <- .find_unique_max(t1)
    expect_equal(max, 170)
    max <- .find_unique_max(t2)
    expect_equal(max, 70)
})

#------------------------------------------------------------------------------#
# Tests .sliding_window
#------------------------------------------------------------------------------#
### OTHER VARS ###
expected_for_smpl1_kf <- tibble::tibble(
    chr = c(rep_len("1", 6)),
    integration_locus = c(
        14571, 14571,
        14571, 14571,
        14577, 14577
    ),
    strand = c(rep_len("+", 6)),
    CompleteAmplificationID = c(
        "ID4",
        "ID1",
        "ID2",
        "ID3",
        "ID5",
        "ID6"
    ),
    Value = c(1400, 70, 150, 120, 36, 15)
)
recalibr_map_smpl1_kf <- tibble::tibble(
    chr_before = rep_len("1", 5),
    integration_locus_before = c(14571, 14572, 14575, 14577, 14581),
    strand_before = rep_len("+", 5),
    chr_after = rep_len("1", 5),
    integration_locus_after = c(14571, 14571, 14571, 14577, 14577),
    strand_after = rep_len("+", 5)
)

expected_for_smpl2_kf <- tibble::tibble(
    chr = c(rep_len("3", 5)),
    integration_locus = c(
        16378, 16387,
        16395, 16395,
        16402
    ),
    strand = c(rep_len("+", 5)),
    CompleteAmplificationID = c(
        "ID1",
        "ID1",
        "ID4",
        "ID3",
        "ID5"
    ),
    Value = c(1969, 48, 153, 886, 543)
)
recalibr_map_smpl2_kf <- tibble::tibble(
    chr_before = rep_len("3", 7),
    integration_locus_before = c(16378, 16380, 16387, 16395, 16396, 16399, 16402),
    strand_before = rep_len("+", 7),
    chr_after = rep_len("3", 7),
    integration_locus_after = c(16378, 16378, 16387, 16395, 16395, 16395, 16402),
    strand_after = rep_len("+", 7)
)

expected_for_smpl2_mv <- tibble::tibble(
    chr = c(rep_len("3", 5)),
    integration_locus = c(
        16380, 16387,
        16399, 16399,
        16402
    ),
    strand = c(rep_len("+", 5)),
    CompleteAmplificationID = c(
        "ID1",
        "ID1",
        "ID4",
        "ID3",
        "ID5"
    ),
    Value = c(1969, 48, 153, 886, 543)
)
recalibr_map_smpl2_mv <- tibble::tibble(
    chr_before = rep_len("3", 7),
    integration_locus_before = c(16378, 16380, 16387, 16395, 16396, 16399, 16402),
    strand_before = rep_len("+", 7),
    chr_after = rep_len("3", 7),
    integration_locus_after = c(16380, 16380, 16387, 16399, 16399, 16399, 16402),
    strand_after = rep_len("+", 7)
)

expected_for_smplmult1_kf <- tibble::tibble(
    chr = c(rep_len("1", 6)),
    integration_locus = c(
        14571, 14571,
        14571, 14571,
        14577, 14577
    ),
    strand = c(rep_len("+", 6)),
    CompleteAmplificationID = c(
        "ID4",
        "ID1",
        "ID2",
        "ID3",
        "ID5",
        "ID6"
    ),
    seqCount = c(
        1400, 70, 150, 120, 36,
        15
    ),
    fragmentEstimate = c(
        64.545, 83.4,
        125.9, 1.656,
        6.4, 564.69
    )
)

expected_for_smplmult2_kf <- tibble::tibble(
    chr = c(rep_len("3", 5)),
    integration_locus = c(
        16378, 16387,
        16395, 16395,
        16402
    ),
    strand = c(rep_len("+", 5)),
    CompleteAmplificationID = c(
        "ID1",
        "ID1",
        "ID4",
        "ID3",
        "ID5"
    ),
    seqCount = c(1969, 48, 153, 886, 543),
    fragmentEstimate = c(
        989.04, 1.36,
        5521.1, 453.60,
        12.55
    )
)
expected_for_smplmult2_mv <- tibble::tibble(
    chr = c(rep_len("3", 5)),
    integration_locus = c(
        16380, 16387,
        16399, 16399,
        16402
    ),
    strand = c(rep_len("+", 5)),
    CompleteAmplificationID = c(
        "ID1",
        "ID1",
        "ID4",
        "ID3",
        "ID5"
    ),
    seqCount = c(
        1969, 48, 153, 886,
        543
    ),
    fragmentEstimate = c(
        989.04, 1.36,
        5521.1, 453.60,
        12.55
    )
)

test_that(".sliding_window produces correct output for sample1", {
    result <- .sliding_window(
        x = sample_group1, threshold = 4,
        keep_criteria = "keep_first", annotated = FALSE,
        num_cols = "Value", max_val_col = "Value",
        produce_map = TRUE
    )
    expect_equal(result$recalibrated_matrix, expected_for_smpl1_kf)
    expect_equal(result$map, recalibr_map_smpl1_kf)
    result <- .sliding_window(
        x = sample_group1, threshold = 4,
        keep_criteria = "max_value", annotated = FALSE,
        num_cols = "Value", max_val_col = "Value",
        produce_map = TRUE
    )
    expect_equal(result$recalibrated_matrix, expected_for_smpl1_kf)
    expect_equal(result$map, recalibr_map_smpl1_kf)
})

test_that(".sliding_window produces correct output for sample2", {
    result <- .sliding_window(
        x = sample_group2, threshold = 4,
        keep_criteria = "keep_first", annotated = FALSE,
        num_cols = "Value", max_val_col = "Value",
        produce_map = TRUE
    )
    expect_equal(result$recalibrated_matrix, expected_for_smpl2_kf)
    expect_equal(result$map, recalibr_map_smpl2_kf)
    result <- .sliding_window(
        x = sample_group2, threshold = 4,
        keep_criteria = "max_value", annotated = FALSE,
        num_cols = "Value", max_val_col = "Value",
        produce_map = TRUE
    )
    expect_equal(result$recalibrated_matrix, expected_for_smpl2_mv)
    expect_equal(result$map, recalibr_map_smpl2_mv)
})

test_that(".sliding_window produces correct output for sample1 - mult column", {
    result <- .sliding_window(
        x = sample_group_mult1, threshold = 4,
        keep_criteria = "keep_first", annotated = FALSE,
        num_cols = c("seqCount", "fragmentEstimate"),
        max_val_col = "seqCount",
        produce_map = TRUE
    )
    expect_equal(result$recalibrated_matrix, expected_for_smplmult1_kf)
    expect_equal(result$map, recalibr_map_smpl1_kf)
    result <- .sliding_window(
        x = sample_group_mult1, threshold = 4,
        keep_criteria = "max_value", annotated = FALSE,
        num_cols = c("seqCount", "fragmentEstimate"),
        max_val_col = "seqCount",
        produce_map = TRUE
    )
    expect_equal(result$recalibrated_matrix, expected_for_smplmult1_kf)
    expect_equal(result$map, recalibr_map_smpl1_kf)
})

test_that(".sliding_window produces correct output for sample2 - mult column", {
    result <- .sliding_window(
        x = sample_group_mult2, threshold = 4,
        keep_criteria = "keep_first", annotated = FALSE,
        num_cols = c("seqCount", "fragmentEstimate"),
        max_val_col = "seqCount",
        produce_map = TRUE
    )
    expect_equal(result$recalibrated_matrix, expected_for_smplmult2_kf)
    expect_equal(result$map, recalibr_map_smpl2_kf)
    result <- .sliding_window(
        x = sample_group_mult2, threshold = 4,
        keep_criteria = "max_value", annotated = FALSE,
        num_cols = c("seqCount", "fragmentEstimate"),
        max_val_col = "seqCount",
        produce_map = TRUE
    )
    expect_equal(result$recalibrated_matrix, expected_for_smplmult2_mv)
    expect_equal(result$map, recalibr_map_smpl2_mv)
})

#------------------------------------------------------------------------------#
# Tests compute_near_integrations
#------------------------------------------------------------------------------#
## Test input
test_that("compute_near_integrations stops if x is not single tibble", {
    expect_error({
        res <- compute_near_integrations(x = 1)
    })
    expect_error({
        res <- compute_near_integrations(x = c(sample_group1, sample_group2))
    })
})
test_that("compute_near_integrations stops if x is not IM", {
    smpl <- sample_group1 %>% dplyr::select(-c("chr"))
    expect_error(
        {
            res <- compute_near_integrations(x = smpl)
        },
        regexp = .non_ISM_error()
    )
})
test_that("compute_near_integrations stops if missing IDs columns", {
    smpl <- sample_group1 %>% dplyr::select(-c("CompleteAmplificationID"))
    expect_error(
        {
            res <- compute_near_integrations(x = smpl)
        },
        regexp = .missing_complAmpID_error()
    )
})
test_that("compute_near_integrations stops if missing any num cols", {
    smpl <- sample_group1 %>% dplyr::select(-c("Value"))
    expect_error(
        {
            res <- compute_near_integrations(x = smpl)
        },
        regexp = .missing_num_cols_error()
    )
})
test_that("compute_near_integrations stops if treshold is incorrect", {
    expect_error({
        res <- compute_near_integrations(x = sample_group1, threshold = "5")
    })
    expect_error({
        res <- compute_near_integrations(x = sample_group1, threshold = c(1, 2))
    })
})
test_that("compute_near_integrations stops if keep_criteria is incorrect", {
    expect_error({
        res <- compute_near_integrations(
            x = sample_group1, threshold = 4,
            keep_criteria = c(1, 2)
        )
    })
    expect_error({
        res <- compute_near_integrations(
            x = sample_group1, threshold = 4,
            keep_criteria = c("a", "b")
        )
    })
    expect_error({
        res <- compute_near_integrations(
            x = sample_group1, threshold = 4,
            keep_criteria = c("max_value", "first")
        )
    })
})
test_that("compute_near_integrations stops if strand_specific is incorrect", {
    expect_error({
        res <- compute_near_integrations(
            x = sample_group1,
            strand_specific = "TRUE"
        )
    })
    expect_error({
        res <- compute_near_integrations(
            x = sample_group1,
            strand_specific = c(TRUE, FALSE)
        )
    })
})
test_that("compute_near_integrations stops if max_value_col is incorrect", {
    expect_error({
        res <- compute_near_integrations(
            x = sample_group1,
            max_value_column = 1
        )
    })
    expect_error({
        res <- compute_near_integrations(
            x = sample_group1,
            max_value_column = c("a", "b")
        )
    })
})
test_that("compute_near_integrations stops if map parameters are incorrect", {
    expect_error({
        res <- compute_near_integrations(
            x = sample_group1,
            map_as_widget = "TRUE"
        )
    })
    expect_error({
        res <- compute_near_integrations(
            x = sample_group1,
            map_as_widget = c(TRUE, FALSE)
        )
    })
    expect_error({
        res <- compute_near_integrations(
            x = sample_group1,
            map_as_file = "TRUE"
        )
    })
    expect_error({
        res <- compute_near_integrations(
            x = sample_group1,
            map_as_file = c(TRUE, FALSE)
        )
    })
})
test_that("compute_near_integrations stops if max_col not found", {
    expect_error(
        {
            res <- compute_near_integrations(
                x = sample_group_mult1,
                keep_criteria = "max_value",
                max_value_column = "a"
            )
        },
        regexp = .max_val_stop_error("a")
    )
})
## TEST VALUES
test_that("compute_near_integrations produces warning for max col", {
    expect_warning(invisible(capture.output({
        res <- compute_near_integrations(
            x = sample_group1,
            keep_criteria = "max_value",
            max_value_column = "seqCount",
            map_as_file = FALSE,
            file_path = NULL
        )
    })), regexp = .using_val_col_warning("seqCount"))
    expect_equal(res, expected_for_smpl1_kf)
})
test_that("compute_near_integrations produces correct output for total", {
    total_simple <- sample_group1 %>% dplyr::bind_rows(sample_group2)
    total_mult <- sample_group_mult1 %>% dplyr::bind_rows(sample_group_mult2)
    expect_warning(
        {
            invisible(capture.output({
                res <- compute_near_integrations(
                    x = total_simple,
                    keep_criteria = "keep_first",
                    max_value_column = "seqCount",
                    map_as_file = FALSE
                )
            }))
        },
        regexp = NA
    )
    expected_simple <- expected_for_smpl1_kf %>%
        dplyr::bind_rows(expected_for_smpl2_kf)
    map_simple_exp <- recalibr_map_smpl1_kf %>%
        dplyr::bind_rows(recalibr_map_smpl2_kf)
    expect_equal(res, expected_simple)
    expect_warning(
        {
            invisible(capture.output({
                res <- compute_near_integrations(
                    x = total_mult,
                    keep_criteria = "keep_first",
                    max_value_column = "seqCount",
                    map_as_file = FALSE
                )
            }))
        },
        regexp = NA
    )
    expected_mult <- expected_for_smplmult1_kf %>%
        dplyr::bind_rows(expected_for_smplmult2_kf)
    expect_equal(res, expected_mult)
})
