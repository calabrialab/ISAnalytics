library(ISAnalytics)

#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
op <- withr::local_options(
    ISAnalytics.reports = FALSE,
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
    expect_equal(result$recalibrated_matrix, expected_for_smpl1_kf,
        ignore_attr = TRUE
    )
    expect_equal(result$map, recalibr_map_smpl1_kf,
        ignore_attr = TRUE
    )
    result <- .sliding_window(
        x = sample_group1, threshold = 4,
        keep_criteria = "max_value", annotated = FALSE,
        num_cols = "Value", max_val_col = "Value",
        produce_map = TRUE
    )
    expect_equal(result$recalibrated_matrix, expected_for_smpl1_kf,
        ignore_attr = TRUE
    )
    expect_equal(result$map, recalibr_map_smpl1_kf,
        ignore_attr = TRUE
    )
})

test_that(".sliding_window produces correct output for sample2", {
    result <- .sliding_window(
        x = sample_group2, threshold = 4,
        keep_criteria = "keep_first", annotated = FALSE,
        num_cols = "Value", max_val_col = "Value",
        produce_map = TRUE
    )
    expect_equal(result$recalibrated_matrix, expected_for_smpl2_kf,
        ignore_attr = TRUE
    )
    expect_equal(result$map, recalibr_map_smpl2_kf,
        ignore_attr = TRUE
    )
    result <- .sliding_window(
        x = sample_group2, threshold = 4,
        keep_criteria = "max_value", annotated = FALSE,
        num_cols = "Value", max_val_col = "Value",
        produce_map = TRUE
    )
    expect_equal(result$recalibrated_matrix, expected_for_smpl2_mv,
        ignore_attr = TRUE
    )
    expect_equal(result$map, recalibr_map_smpl2_mv,
        ignore_attr = TRUE
    )
})

test_that(".sliding_window produces correct output for sample1 - mult column", {
    result <- .sliding_window(
        x = sample_group_mult1, threshold = 4,
        keep_criteria = "keep_first", annotated = FALSE,
        num_cols = c("seqCount", "fragmentEstimate"),
        max_val_col = "seqCount",
        produce_map = TRUE
    )
    expect_equal(result$recalibrated_matrix, expected_for_smplmult1_kf,
        ignore_attr = TRUE
    )
    expect_equal(result$map, recalibr_map_smpl1_kf,
        ignore_attr = TRUE
    )
    result <- .sliding_window(
        x = sample_group_mult1, threshold = 4,
        keep_criteria = "max_value", annotated = FALSE,
        num_cols = c("seqCount", "fragmentEstimate"),
        max_val_col = "seqCount",
        produce_map = TRUE
    )
    expect_equal(result$recalibrated_matrix, expected_for_smplmult1_kf,
        ignore_attr = TRUE
    )
    expect_equal(result$map, recalibr_map_smpl1_kf,
        ignore_attr = TRUE
    )
})

test_that(".sliding_window produces correct output for sample2 - mult column", {
    result <- .sliding_window(
        x = sample_group_mult2, threshold = 4,
        keep_criteria = "keep_first", annotated = FALSE,
        num_cols = c("seqCount", "fragmentEstimate"),
        max_val_col = "seqCount",
        produce_map = TRUE
    )
    expect_equal(result$recalibrated_matrix, expected_for_smplmult2_kf,
        ignore_attr = TRUE
    )
    expect_equal(result$map, recalibr_map_smpl2_kf,
        ignore_attr = TRUE
    )
    result <- .sliding_window(
        x = sample_group_mult2, threshold = 4,
        keep_criteria = "max_value", annotated = FALSE,
        num_cols = c("seqCount", "fragmentEstimate"),
        max_val_col = "seqCount",
        produce_map = TRUE
    )
    expect_equal(result$recalibrated_matrix, expected_for_smplmult2_mv,
        ignore_attr = TRUE
    )
    expect_equal(result$map, recalibr_map_smpl2_mv,
        ignore_attr = TRUE
    )
})

#------------------------------------------------------------------------------#
# Tests compute_near_integrations
#------------------------------------------------------------------------------#
## TEST VALUES
test_that("compute_near_integrations produces correct output for total", {
    total_simple <- sample_group1 %>% dplyr::bind_rows(sample_group2)
    total_mult <- sample_group_mult1 %>% dplyr::bind_rows(sample_group_mult2)
    res <- compute_near_integrations(
        x = total_simple,
        keep_criteria = "keep_first",
        max_value_column = "Value",
        value_columns = c("Value"),
        map_as_file = FALSE
    )
    expected_simple <- expected_for_smpl1_kf %>%
        dplyr::bind_rows(expected_for_smpl2_kf)
    map_simple_exp <- recalibr_map_smpl1_kf %>%
        dplyr::bind_rows(recalibr_map_smpl2_kf)
    expect_equal(res, expected_simple,
        ignore_attr = TRUE
    )
    res <- compute_near_integrations(
        x = total_mult,
        keep_criteria = "keep_first",
        max_value_column = "seqCount",
        value_columns = c("seqCount", "fragmentEstimate"),
        map_as_file = FALSE
    )
    expected_mult <- expected_for_smplmult1_kf %>%
        dplyr::bind_rows(expected_for_smplmult2_kf)
    expect_equal(res, expected_mult,
        ignore_attr = TRUE
    )
})
