#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
withr::local_options(
    ISAnalytics.verbose = FALSE,
    ISAnalytics.reports = FALSE
)

minimal_agg_example <- tibble::tibble(
    chr = c("1", "2", "2", "5", "5", "3", "3", "10", "11"),
    integration_locus = c(
        14505,
        1005,
        15513,
        4561,
        10167,
        5247,
        10951,
        23403,
        25611
    ),
    strand = c(
        "+", "-", "+", "+", "-", "-",
        "+", "-", "-"
    ),
    SubjectID = c(
        rep_len("S1", 5),
        rep_len("S2", 4)
    ),
    CellMarker = c(
        rep_len("CM1", 5),
        rep_len("CM2", 4)
    ),
    Tissue = c(
        rep_len("T1", 5),
        rep_len("T2", 4)
    ),
    TimePoint = c(
        rep_len("0001", 3),
        "0010", "0015",
        rep_len("0001", 2),
        rep_len("0020", 2)
    ),
    Value_sum = c(
        50, 80, 650, 46, 79, 633,
        875, 99, 123
    )
)

#------------------------------------------------------------------------------#
# Test sample_statistics
#------------------------------------------------------------------------------#
test_that("sample_statistics works as expected with default functions", {
    val_cols <- c("seqCount", "fragmentEstimate")
    stats <- sample_statistics(
        x = integration_matrices,
        metadata = association_file,
        value_columns = val_cols
    )
    fun_names <- names(default_stats())
    fun_names <- fun_names[fun_names != "describe"]
    desc_fun_names <- paste("describe", c(
        "vars", "n", "mean", "sd",
        "median", "trimmed", "mad",
        "min", "max", "range", "skew",
        "kurtosis", "se"
    ), sep = "_")
    fun_names <- c(fun_names, desc_fun_names)
    expected_name_cols <- unlist(purrr::map(
        val_cols,
        ~ paste(.x, fun_names, sep = "_")
    ))
    expect_true(all(c(expected_name_cols, "nIS") %in% colnames(stats$x)))
    expect_true(all(c(expected_name_cols, "nIS") %in% colnames(stats$metadata)))
})

#------------------------------------------------------------------------------#
# Test top_integrations
#------------------------------------------------------------------------------#
test_that("top_integrations works as expected", {
    example <- tibble::tribble(
        ~chr, ~integration_locus, ~strand, ~SubjectID, ~Tissue, ~abundance,
        "1", 34254, "+", "S1", "T1", 0.27,
        "1", 56423, "+", "S1", "T1", 0.02,
        "2", 85834, "-", "S1", "T1", 0.12,
        "5", 12334, "-", "S1", "T1", 0.05,
        "6", 64721, "+", "S1", "T1", 0.16,
        "1", 96473, "-", "S1", "T1", 0.29,
        "7", 31434, "+", "S1", "T1", 0.09,
        "1", 34254, "+", "S2", "T1", 0.56,
        "X", 31246, "+", "S2", "T1", 0.01,
        "8", 12354, "-", "S2", "T1", 0.08,
        "9", 13468, "+", "S2", "T1", 0.15,
        "5", 86534, "-", "S2", "T1", 0.19,
        "1", 74567, "+", "S2", "T1", 0.01
    )
    top_3 <- top_integrations(example,
        n = 3,
        columns = "abundance"
    )
    expected_top_3 <- tibble::tribble(
        ~chr, ~integration_locus, ~strand, ~abundance, ~SubjectID, ~Tissue,
        "1", 34254, "+", 0.56, "S2", "T1",
        "1", 96473, "-", 0.29, "S1", "T1",
        "1", 34254, "+", 0.27, "S1", "T1"
    )
    expect_equal(top_3, expected_top_3)
    top_3 <- top_integrations(example,
        n = 3,
        columns = "abundance",
        key = c("SubjectID", "Tissue")
    )
    expected_top_3 <- tibble::tribble(
        ~SubjectID, ~Tissue, ~chr, ~integration_locus, ~strand, ~abundance,
        "S1", "T1", "1", 96473, "-", 0.29,
        "S1", "T1", "1", 34254, "+", 0.27,
        "S1", "T1", "6", 64721, "+", 0.16,
        "S2", "T1", "1", 34254, "+", 0.56,
        "S2", "T1", "5", 86534, "-", 0.19,
        "S2", "T1", "9", 13468, "+", 0.15,
    )
    expect_equal(top_3, expected_top_3)
})

#------------------------------------------------------------------------------#
# Test cumulative_is
#------------------------------------------------------------------------------#
test_that("cumulative_is produces correct output", {
    c_is <- cumulative_is(minimal_agg_example,
        expand = FALSE,
        keep_og_is = TRUE
    )
    expect_true(nrow(c_is$coordinates) == 5)
    expect_true(all(c("is", "cumulative_is") %in% colnames(c_is$coordinates)))
    expect_equal(c_is$counts$is_n_cumulative, c(3, 4, 5, 2, 4))
    c_is <- cumulative_is(minimal_agg_example,
        keep_og_is = FALSE,
        expand = FALSE
    )
    expect_true(nrow(c_is$coordinates) == 5)
    expect_true(all(c("cumulative_is") %in% colnames(c_is$coordinates)))
    expect_equal(c_is$counts$is_n_cumulative, c(3, 4, 5, 2, 4))
    c_is <- cumulative_is(minimal_agg_example,
        keep_og_is = FALSE,
        expand = TRUE
    )
    expect_true(all(mandatory_IS_vars() %in% colnames(c_is$coordinates)))
})

test_that("cumulative_is works as expected with tp 0", {
    all_zeros <- minimal_agg_example |>
        dplyr::mutate(TimePoint = "0000")
    expect_message(
        {
            c_is <- cumulative_is(all_zeros)
        },
        class = "only_zero_tps"
    )
    expect_null(c_is)

    some_zeros <- minimal_agg_example
    some_zeros[c(1, 2, 3), ]$TimePoint <- "0000"
    c_is <- cumulative_is(some_zeros, expand = FALSE)
    expect_true(c_is$counts[1, ]$TimePoint == 10 &
        c_is$counts[1, ]$is_n_cumulative == 1)
})

#------------------------------------------------------------------------------#
# Test .assign_iss_by_tp
#------------------------------------------------------------------------------#
test_that(".assign_iss_by_tp assigns iss correctly", {
    test_df <- tibble::tribble(
        ~chr, ~integration_locus, ~strand, ~CellMarker, ~Tissue, ~TimePoint,
        "1", 12345, "+", "CD34", "BM", "03",
        "1", 12345, "+", "CD34", "BM", "01",
        "1", 12345, "+", "CD34", "BM", "06",
        "2", 54321, "-", "CD34", "BM", "01",
        "3", 62135, "+", "CD34", "BM", "02",
        "3", 62135, "+", "CD34", "BM", "08"
    )

    expected_assign <- tibble::tribble(
        ~chr, ~integration_locus, ~strand, ~CellMarker, ~Tissue, ~TimePoint,
        "1", 12345, "+", "CD34", "BM", 1,
        "2", 54321, "-", "CD34", "BM", 1,
        "3", 62135, "+", "CD34", "BM", 2
    )
    re_assigned <- .assign_iss_by_tp(
        df = test_df,
        timepoint_column = "TimePoint"
    )
    expect_equal(re_assigned, expected_assign)
})

#------------------------------------------------------------------------------#
# Test iss_source
#------------------------------------------------------------------------------#
test_df2 <- tibble::tibble(
    chr = c(
        "1", "1", "1", "1", "2", "2", "3", "3",
        "1", "1", "1", "1", "1", "2", "3", "3"
    ),
    integration_locus = c(
        12345, 43524, 12345, 12345, 54321,
        76835, 62135, 62135, 12345, 12345,
        56832, 12345, 12345, 54321, 62135,
        62135
    ),
    strand = c(
        "+", "+", "+", "+", "-", "-", "+", "+", "+",
        "+", "-", "+", "+", "-", "+", "+"
    ),
    SubjectID = c(
        "PT01", "PT02", "PT01", "PT01", "PT01",
        "PT02", "PT01", "PT01", "PT01", "PT01",
        "PT02", "PT01", "PT02", "PT01", "PT01",
        "PT01"
    ),
    CellMarker = c(
        "CD34", "CD34", "CD34", "CD34", "CD34",
        "CD34", "CD34", "CD34", "Whole", "Whole",
        "Whole", "Whole", "Whole", "Whole", "Whole",
        "Whole"
    ),
    Tissue = c(
        "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM",
        "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM"
    ),
    TimePoint = c(
        "03", "03", "01", "06", "01", "01", "02",
        "08", "03", "01", "01", "06", "06", "01",
        "02", "08"
    )
)

test_df3 <- tibble::tibble(
    chr = c(
        "1", "1", "1", "1", "2", "2", "3", "3", "1",
        "1", "1", "1", "1", "2", "3", "3"
    ),
    integration_locus = c(
        56252, 43524, 12345, 86435, 54321,
        76835, 62135, 432245, 12345, 534587,
        56832, 12345, 12345, 578635, 62135,
        62135
    ),
    strand = c(
        "+", "+", "+", "+", "-", "-", "+", "+", "+",
        "+", "-", "+", "+", "-", "+", "+"
    ),
    SubjectID = c(
        "PT02", "PT02", "PT01", "PT01", "PT01",
        "PT02", "PT01", "PT01", "PT01", "PT01",
        "PT02", "PT01", "PT02", "PT01", "PT02",
        "PT01"
    ),
    CellMarker = c(
        "CD14", "CD13", "CD15", "CD14", "CD14",
        "CD13", "CD14", "CD13", "CD14", "CD15",
        "CD15", "CD14", "CD15", "CD13", "CD14",
        "CD13"
    ),
    Tissue = c(
        "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB",
        "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB"
    ),
    TimePoint = c(
        "03", "03", "01", "06", "01", "01", "02",
        "08", "03", "01", "01", "06", "06", "01",
        "02", "08"
    )
)
expected_res1 <- list(
    PT01 = tibble::tibble(
        g1 = c(
            "PT01_CD34_BM_1", "PT01_CD34_BM_2",
            "PT01_CD34_BM_1", "PT01_CD34_BM_2",
            "PT01_CD34_BM_1", "PT01_CD34_BM_2",
            "PT01_CD34_BM_1", "PT01_CD34_BM_2",
            "PT01_CD34_BM_1", "PT01_CD34_BM_2",
            "PT01_CD34_BM_1", "PT01_CD34_BM_2",
            "PT01_CD34_BM_1", "PT01_CD34_BM_2",
            "PT01_Whole_BM_1", "PT01_Whole_BM_2",
            "PT01_Whole_BM_1", "PT01_Whole_BM_2",
            "PT01_Whole_BM_1", "PT01_Whole_BM_2",
            "PT01_Whole_BM_1", "PT01_Whole_BM_2",
            "PT01_Whole_BM_1", "PT01_Whole_BM_2",
            "PT01_Whole_BM_1", "PT01_Whole_BM_2",
            "PT01_Whole_BM_1", "PT01_Whole_BM_2"
        ) |> paste0("-1"),
        g1_SubjectID = c(
            "PT01", "PT01", "PT01", "PT01",
            "PT01", "PT01", "PT01", "PT01",
            "PT01", "PT01", "PT01", "PT01",
            "PT01", "PT01", "PT01", "PT01",
            "PT01", "PT01", "PT01", "PT01",
            "PT01", "PT01", "PT01", "PT01",
            "PT01", "PT01", "PT01", "PT01"
        ),
        g1_CellMarker = c(
            "CD34", "CD34", "CD34", "CD34",
            "CD34", "CD34", "CD34", "CD34",
            "CD34", "CD34", "CD34", "CD34",
            "CD34", "CD34", "Whole", "Whole",
            "Whole", "Whole", "Whole", "Whole",
            "Whole", "Whole", "Whole", "Whole",
            "Whole", "Whole", "Whole", "Whole"
        ),
        g1_Tissue = c(
            "BM", "BM", "BM", "BM", "BM", "BM",
            "BM", "BM", "BM", "BM", "BM", "BM",
            "BM", "BM", "BM", "BM", "BM", "BM",
            "BM", "BM", "BM", "BM", "BM", "BM",
            "BM", "BM", "BM", "BM"
        ),
        g1_TimePoint = c(
            1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2,
            1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2,
            1, 2, 1, 2
        ),
        g2 = c(
            "PT01_CD15_PB_01", "PT01_CD15_PB_01",
            "PT01_CD14_PB_06", "PT01_CD14_PB_06",
            "PT01_CD14_PB_01", "PT01_CD14_PB_01",
            "PT01_CD14_PB_02", "PT01_CD14_PB_02",
            "PT01_CD13_PB_08", "PT01_CD13_PB_08",
            "PT01_CD14_PB_03", "PT01_CD14_PB_03",
            "PT01_CD13_PB_01", "PT01_CD13_PB_01",
            "PT01_CD15_PB_01", "PT01_CD15_PB_01",
            "PT01_CD14_PB_06", "PT01_CD14_PB_06",
            "PT01_CD14_PB_01", "PT01_CD14_PB_01",
            "PT01_CD14_PB_02", "PT01_CD14_PB_02",
            "PT01_CD13_PB_08", "PT01_CD13_PB_08",
            "PT01_CD14_PB_03", "PT01_CD14_PB_03",
            "PT01_CD13_PB_01", "PT01_CD13_PB_01"
        ) |> paste0("-2"),
        g2_SubjectID = c(
            "PT01", "PT01", "PT01", "PT01",
            "PT01", "PT01", "PT01", "PT01",
            "PT01", "PT01", "PT01", "PT01",
            "PT01", "PT01", "PT01", "PT01",
            "PT01", "PT01", "PT01", "PT01",
            "PT01", "PT01", "PT01", "PT01",
            "PT01", "PT01", "PT01", "PT01"
        ),
        g2_CellMarker = c(
            "CD15", "CD15", "CD14", "CD14",
            "CD14", "CD14", "CD14", "CD14",
            "CD13", "CD13", "CD14", "CD14",
            "CD13", "CD13", "CD15", "CD15",
            "CD14", "CD14", "CD14", "CD14",
            "CD14", "CD14", "CD13", "CD13",
            "CD14", "CD14", "CD13", "CD13"
        ),
        g2_Tissue = c(
            "PB", "PB", "PB", "PB", "PB", "PB",
            "PB", "PB", "PB", "PB", "PB", "PB",
            "PB", "PB", "PB", "PB", "PB", "PB",
            "PB", "PB", "PB", "PB", "PB", "PB",
            "PB", "PB", "PB", "PB"
        ),
        g2_TimePoint = c(
            1, 1, 6, 6, 1, 1, 2, 2, 8, 8, 3, 3,
            1, 1, 1, 1, 6, 6, 1, 1, 2, 2, 8, 8,
            3, 3, 1, 1
        ),
        chr = c(
            "1", NA, "1", NA, "2", NA, NA, "3",
            NA, "3", "1", NA, NA, NA, "1", NA,
            "1", NA, "2", NA, NA, "3", NA, "3",
            "1", NA, NA, NA
        ),
        integration_locus = c(
            12345, NA, 12345, NA, 54321,
            NA, NA, 62135, NA, 62135, 12345,
            NA, NA, NA, 12345, NA, 12345,
            NA, 54321, NA, NA, 62135, NA,
            62135, 12345, NA, NA, NA
        ),
        strand = c(
            "+", NA, "+", NA, "-", NA, NA,
            "+", NA, "+", "+", NA, NA, NA,
            "+", NA, "+", NA, "-", NA, NA,
            "+", NA, "+", "+", NA, NA, NA
        ),
        sharing_perc = c(
            50, 0, 50, 0, 100, 0, 0, 100, 0, 50,
            100, 0, 0, 0, 50, 0, 50, 0, 100, 0,
            0, 100, 0, 50, 100, 0, 0, 0
        )
    ) |>
        dplyr::arrange(.data$g1, .data$g2),
    PT02 = tibble::tibble(
        g1 = c(
            "PT02_CD34_BM_3", "PT02_CD34_BM_1",
            "PT02_CD34_BM_3", "PT02_CD34_BM_1",
            "PT02_CD34_BM_3", "PT02_CD34_BM_1",
            "PT02_CD34_BM_3", "PT02_CD34_BM_1",
            "PT02_CD34_BM_3", "PT02_CD34_BM_1",
            "PT02_CD34_BM_3", "PT02_CD34_BM_1",
            "PT02_Whole_BM_6", "PT02_Whole_BM_1",
            "PT02_Whole_BM_6", "PT02_Whole_BM_1",
            "PT02_Whole_BM_6", "PT02_Whole_BM_1",
            "PT02_Whole_BM_6", "PT02_Whole_BM_1",
            "PT02_Whole_BM_6", "PT02_Whole_BM_1",
            "PT02_Whole_BM_6", "PT02_Whole_BM_1"
        ) |> paste0("-1"),
        g1_SubjectID = c(
            "PT02", "PT02", "PT02", "PT02",
            "PT02", "PT02", "PT02", "PT02",
            "PT02", "PT02", "PT02", "PT02",
            "PT02", "PT02", "PT02", "PT02",
            "PT02", "PT02", "PT02", "PT02",
            "PT02", "PT02", "PT02", "PT02"
        ),
        g1_CellMarker = c(
            "CD34", "CD34", "CD34", "CD34",
            "CD34", "CD34", "CD34", "CD34",
            "CD34", "CD34", "CD34", "CD34",
            "Whole", "Whole", "Whole", "Whole",
            "Whole", "Whole", "Whole", "Whole",
            "Whole", "Whole", "Whole", "Whole"
        ),
        g1_Tissue = c(
            "BM", "BM", "BM", "BM", "BM", "BM",
            "BM", "BM", "BM", "BM", "BM", "BM",
            "BM", "BM", "BM", "BM", "BM", "BM",
            "BM", "BM", "BM", "BM", "BM", "BM"
        ),
        g1_TimePoint = c(
            3, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1,
            6, 1, 6, 1, 6, 1, 6, 1, 6, 1, 6, 1
        ),
        g2 = c(
            "PT02_CD14_PB_03", "PT02_CD14_PB_03",
            "PT02_CD13_PB_03", "PT02_CD13_PB_03",
            "PT02_CD13_PB_01", "PT02_CD13_PB_01",
            "PT02_CD15_PB_01", "PT02_CD15_PB_01",
            "PT02_CD15_PB_06", "PT02_CD15_PB_06",
            "PT02_CD14_PB_02", "PT02_CD14_PB_02",
            "PT02_CD14_PB_03", "PT02_CD14_PB_03",
            "PT02_CD13_PB_03", "PT02_CD13_PB_03",
            "PT02_CD13_PB_01", "PT02_CD13_PB_01",
            "PT02_CD15_PB_01", "PT02_CD15_PB_01",
            "PT02_CD15_PB_06", "PT02_CD15_PB_06",
            "PT02_CD14_PB_02", "PT02_CD14_PB_02"
        ) |> paste0("-2"),
        g2_SubjectID = c(
            "PT02", "PT02", "PT02", "PT02",
            "PT02", "PT02", "PT02", "PT02",
            "PT02", "PT02", "PT02", "PT02",
            "PT02", "PT02", "PT02", "PT02",
            "PT02", "PT02", "PT02", "PT02",
            "PT02", "PT02", "PT02", "PT02"
        ),
        g2_CellMarker = c(
            "CD14", "CD14", "CD13", "CD13",
            "CD13", "CD13", "CD15", "CD15",
            "CD15", "CD15", "CD14", "CD14",
            "CD14", "CD14", "CD13", "CD13",
            "CD13", "CD13", "CD15", "CD15",
            "CD15", "CD15", "CD14", "CD14"
        ),
        g2_Tissue = c(
            "PB", "PB", "PB", "PB", "PB", "PB",
            "PB", "PB", "PB", "PB", "PB", "PB",
            "PB", "PB", "PB", "PB", "PB", "PB",
            "PB", "PB", "PB", "PB", "PB", "PB"
        ),
        g2_TimePoint = c(
            3, 3, 3, 3, 1, 1, 1, 1, 6, 6, 2, 2,
            3, 3, 3, 3, 1, 1, 1, 1, 6, 6, 2, 2
        ),
        chr = c(
            NA, NA, "1", NA, NA, "2", NA, NA,
            NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, 1, 1, NA, NA,
            NA
        ),
        integration_locus = c(
            NA, NA, 43524, NA, NA, 76835,
            NA, NA, NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, 56832,
            12345, NA, NA, NA
        ),
        strand = c(
            NA, NA, "+", NA, NA, "-", NA,
            NA, NA, NA, NA, NA, NA,
            NA, NA, NA, NA, NA, NA,
            "-", "+", NA, NA, NA
        ),
        sharing_perc = c(
            0, 0, 100, 0, 0, 100, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 100,
            100, 0, 0, 0
        )
    ) |>
        dplyr::arrange(.data$g1, .data$g2)
)

expected_res2 <- tibble::tibble(
    g1 = c(
        "PT01_CD34_BM_1", "PT01_CD34_BM_2", "PT01_CD34_BM_1", "PT01_CD34_BM_2",
        "PT01_CD34_BM_1", "PT01_CD34_BM_2", "PT01_CD34_BM_1", "PT01_CD34_BM_2",
        "PT01_CD34_BM_1", "PT01_CD34_BM_2", "PT01_CD34_BM_1", "PT01_CD34_BM_2",
        "PT01_CD34_BM_1", "PT01_CD34_BM_2", "PT01_CD34_BM_1", "PT01_CD34_BM_2",
        "PT01_CD34_BM_1", "PT01_CD34_BM_2", "PT01_CD34_BM_1", "PT01_CD34_BM_2",
        "PT01_CD34_BM_1", "PT01_CD34_BM_2", "PT01_CD34_BM_1", "PT01_CD34_BM_2",
        "PT01_CD34_BM_1", "PT01_CD34_BM_2", "PT01_Whole_BM_1",
        "PT01_Whole_BM_2", "PT01_Whole_BM_1", "PT01_Whole_BM_2",
        "PT01_Whole_BM_1", "PT01_Whole_BM_2", "PT01_Whole_BM_1",
        "PT01_Whole_BM_2", "PT01_Whole_BM_1", "PT01_Whole_BM_2",
        "PT01_Whole_BM_1", "PT01_Whole_BM_2", "PT01_Whole_BM_1",
        "PT01_Whole_BM_2", "PT01_Whole_BM_1", "PT01_Whole_BM_2",
        "PT01_Whole_BM_1", "PT01_Whole_BM_2", "PT01_Whole_BM_1",
        "PT01_Whole_BM_2", "PT01_Whole_BM_1", "PT01_Whole_BM_2",
        "PT01_Whole_BM_1", "PT01_Whole_BM_2", "PT01_Whole_BM_1",
        "PT01_Whole_BM_2", "PT02_CD34_BM_3", "PT02_CD34_BM_1",
        "PT02_CD34_BM_3", "PT02_CD34_BM_1", "PT02_CD34_BM_3",
        "PT02_CD34_BM_1", "PT02_CD34_BM_3", "PT02_CD34_BM_1",
        "PT02_CD34_BM_3", "PT02_CD34_BM_1", "PT02_CD34_BM_3",
        "PT02_CD34_BM_1", "PT02_CD34_BM_3", "PT02_CD34_BM_1",
        "PT02_CD34_BM_3", "PT02_CD34_BM_1", "PT02_CD34_BM_3",
        "PT02_CD34_BM_1", "PT02_CD34_BM_3", "PT02_CD34_BM_1",
        "PT02_CD34_BM_3", "PT02_CD34_BM_1", "PT02_CD34_BM_3",
        "PT02_CD34_BM_1", "PT02_CD34_BM_3", "PT02_CD34_BM_1",
        "PT02_Whole_BM_6", "PT02_Whole_BM_1", "PT02_Whole_BM_6",
        "PT02_Whole_BM_1", "PT02_Whole_BM_6", "PT02_Whole_BM_1",
        "PT02_Whole_BM_6", "PT02_Whole_BM_1", "PT02_Whole_BM_6",
        "PT02_Whole_BM_1", "PT02_Whole_BM_6", "PT02_Whole_BM_1",
        "PT02_Whole_BM_6", "PT02_Whole_BM_1", "PT02_Whole_BM_6",
        "PT02_Whole_BM_1", "PT02_Whole_BM_6", "PT02_Whole_BM_1",
        "PT02_Whole_BM_6", "PT02_Whole_BM_1", "PT02_Whole_BM_6",
        "PT02_Whole_BM_1", "PT02_Whole_BM_6", "PT02_Whole_BM_1",
        "PT02_Whole_BM_6", "PT02_Whole_BM_1"
    ) |> paste0("-1"),
    g1_SubjectID = c(
        "PT01", "PT01", "PT01", "PT01", "PT01", "PT01", "PT01",
        "PT01", "PT01", "PT01", "PT01", "PT01", "PT01", "PT01",
        "PT01", "PT01", "PT01", "PT01", "PT01", "PT01", "PT01",
        "PT01", "PT01", "PT01", "PT01", "PT01", "PT01", "PT01",
        "PT01", "PT01", "PT01", "PT01", "PT01", "PT01", "PT01",
        "PT01", "PT01", "PT01", "PT01", "PT01", "PT01", "PT01",
        "PT01", "PT01", "PT01", "PT01", "PT01", "PT01", "PT01",
        "PT01", "PT01", "PT01", "PT02", "PT02", "PT02", "PT02",
        "PT02", "PT02", "PT02", "PT02", "PT02", "PT02", "PT02",
        "PT02", "PT02", "PT02", "PT02", "PT02", "PT02", "PT02",
        "PT02", "PT02", "PT02", "PT02", "PT02", "PT02", "PT02",
        "PT02", "PT02", "PT02", "PT02", "PT02", "PT02", "PT02",
        "PT02", "PT02", "PT02", "PT02", "PT02", "PT02", "PT02",
        "PT02", "PT02", "PT02", "PT02", "PT02", "PT02", "PT02",
        "PT02", "PT02", "PT02", "PT02", "PT02", "PT02"
    ),
    g1_CellMarker = c(
        "CD34", "CD34", "CD34", "CD34", "CD34", "CD34", "CD34",
        "CD34", "CD34", "CD34", "CD34", "CD34", "CD34", "CD34",
        "CD34", "CD34", "CD34", "CD34", "CD34", "CD34", "CD34",
        "CD34", "CD34", "CD34", "CD34", "CD34", "Whole", "Whole",
        "Whole", "Whole", "Whole", "Whole", "Whole", "Whole",
        "Whole", "Whole", "Whole", "Whole", "Whole", "Whole",
        "Whole", "Whole", "Whole", "Whole", "Whole", "Whole",
        "Whole", "Whole", "Whole", "Whole", "Whole", "Whole",
        "CD34", "CD34", "CD34", "CD34", "CD34", "CD34", "CD34",
        "CD34", "CD34", "CD34", "CD34", "CD34", "CD34", "CD34",
        "CD34", "CD34", "CD34", "CD34", "CD34", "CD34", "CD34",
        "CD34", "CD34", "CD34", "CD34", "CD34", "Whole", "Whole",
        "Whole", "Whole", "Whole", "Whole", "Whole", "Whole",
        "Whole", "Whole", "Whole", "Whole", "Whole", "Whole",
        "Whole", "Whole", "Whole", "Whole", "Whole", "Whole",
        "Whole", "Whole", "Whole", "Whole", "Whole", "Whole"
    ),
    g1_Tissue = c(
        "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM",
        "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM",
        "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM",
        "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM",
        "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM",
        "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM",
        "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM",
        "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM",
        "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM",
        "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM",
        "BM", "BM", "BM", "BM"
    ),
    g1_TimePoint = c(
        1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1,
        2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2,
        1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 3, 1, 3, 1, 3,
        1, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1,
        3, 1, 6, 1, 6, 1, 6, 1, 6, 1, 6, 1, 6, 1, 6, 1, 6, 1, 6,
        1, 6, 1, 6, 1, 6, 1, 6, 1
    ),
    g2 = c(
        "PT02_CD14_PB_03", "PT02_CD14_PB_03", "PT02_CD13_PB_03",
        "PT02_CD13_PB_03", "PT01_CD15_PB_01", "PT01_CD15_PB_01",
        "PT01_CD14_PB_06", "PT01_CD14_PB_06", "PT01_CD14_PB_01",
        "PT01_CD14_PB_01", "PT02_CD13_PB_01", "PT02_CD13_PB_01",
        "PT01_CD14_PB_02", "PT01_CD14_PB_02", "PT01_CD13_PB_08",
        "PT01_CD13_PB_08", "PT01_CD14_PB_03", "PT01_CD14_PB_03",
        "PT02_CD15_PB_01", "PT02_CD15_PB_01", "PT02_CD15_PB_06",
        "PT02_CD15_PB_06", "PT01_CD13_PB_01", "PT01_CD13_PB_01",
        "PT02_CD14_PB_02", "PT02_CD14_PB_02", "PT02_CD14_PB_03",
        "PT02_CD14_PB_03", "PT02_CD13_PB_03", "PT02_CD13_PB_03",
        "PT01_CD15_PB_01", "PT01_CD15_PB_01", "PT01_CD14_PB_06",
        "PT01_CD14_PB_06", "PT01_CD14_PB_01", "PT01_CD14_PB_01",
        "PT02_CD13_PB_01", "PT02_CD13_PB_01", "PT01_CD14_PB_02",
        "PT01_CD14_PB_02", "PT01_CD13_PB_08", "PT01_CD13_PB_08",
        "PT01_CD14_PB_03", "PT01_CD14_PB_03", "PT02_CD15_PB_01",
        "PT02_CD15_PB_01", "PT02_CD15_PB_06", "PT02_CD15_PB_06",
        "PT01_CD13_PB_01", "PT01_CD13_PB_01", "PT02_CD14_PB_02",
        "PT02_CD14_PB_02", "PT02_CD14_PB_03", "PT02_CD14_PB_03",
        "PT02_CD13_PB_03", "PT02_CD13_PB_03", "PT01_CD15_PB_01",
        "PT01_CD15_PB_01", "PT01_CD14_PB_06", "PT01_CD14_PB_06",
        "PT01_CD14_PB_01", "PT01_CD14_PB_01", "PT02_CD13_PB_01",
        "PT02_CD13_PB_01", "PT01_CD14_PB_02", "PT01_CD14_PB_02",
        "PT01_CD13_PB_08", "PT01_CD13_PB_08", "PT01_CD14_PB_03",
        "PT01_CD14_PB_03", "PT02_CD15_PB_01", "PT02_CD15_PB_01",
        "PT02_CD15_PB_06", "PT02_CD15_PB_06", "PT01_CD13_PB_01",
        "PT01_CD13_PB_01", "PT02_CD14_PB_02", "PT02_CD14_PB_02",
        "PT02_CD14_PB_03", "PT02_CD14_PB_03", "PT02_CD13_PB_03",
        "PT02_CD13_PB_03", "PT01_CD15_PB_01", "PT01_CD15_PB_01",
        "PT01_CD14_PB_06", "PT01_CD14_PB_06", "PT01_CD14_PB_01",
        "PT01_CD14_PB_01", "PT02_CD13_PB_01", "PT02_CD13_PB_01",
        "PT01_CD14_PB_02", "PT01_CD14_PB_02", "PT01_CD13_PB_08",
        "PT01_CD13_PB_08", "PT01_CD14_PB_03", "PT01_CD14_PB_03",
        "PT02_CD15_PB_01", "PT02_CD15_PB_01", "PT02_CD15_PB_06",
        "PT02_CD15_PB_06", "PT01_CD13_PB_01", "PT01_CD13_PB_01",
        "PT02_CD14_PB_02", "PT02_CD14_PB_02"
    ) |> paste0("-2"),
    g2_SubjectID = c(
        "PT02", "PT02", "PT02", "PT02", "PT01", "PT01", "PT01",
        "PT01", "PT01", "PT01", "PT02", "PT02", "PT01", "PT01",
        "PT01", "PT01", "PT01", "PT01", "PT02", "PT02", "PT02",
        "PT02", "PT01", "PT01", "PT02", "PT02", "PT02", "PT02",
        "PT02", "PT02", "PT01", "PT01", "PT01", "PT01", "PT01",
        "PT01", "PT02", "PT02", "PT01", "PT01", "PT01", "PT01",
        "PT01", "PT01", "PT02", "PT02", "PT02", "PT02", "PT01",
        "PT01", "PT02", "PT02", "PT02", "PT02", "PT02", "PT02",
        "PT01", "PT01", "PT01", "PT01", "PT01", "PT01", "PT02",
        "PT02", "PT01", "PT01", "PT01", "PT01", "PT01", "PT01",
        "PT02", "PT02", "PT02", "PT02", "PT01", "PT01", "PT02",
        "PT02", "PT02", "PT02", "PT02", "PT02", "PT01", "PT01",
        "PT01", "PT01", "PT01", "PT01", "PT02", "PT02", "PT01",
        "PT01", "PT01", "PT01", "PT01", "PT01", "PT02", "PT02",
        "PT02", "PT02", "PT01", "PT01", "PT02", "PT02"
    ),
    g2_CellMarker = c(
        "CD14", "CD14", "CD13", "CD13", "CD15", "CD15", "CD14",
        "CD14", "CD14", "CD14", "CD13", "CD13", "CD14", "CD14",
        "CD13", "CD13", "CD14", "CD14", "CD15", "CD15", "CD15",
        "CD15", "CD13", "CD13", "CD14", "CD14", "CD14", "CD14",
        "CD13", "CD13", "CD15", "CD15", "CD14", "CD14", "CD14",
        "CD14", "CD13", "CD13", "CD14", "CD14", "CD13", "CD13",
        "CD14", "CD14", "CD15", "CD15", "CD15", "CD15", "CD13",
        "CD13", "CD14", "CD14", "CD14", "CD14", "CD13", "CD13",
        "CD15", "CD15", "CD14", "CD14", "CD14", "CD14", "CD13",
        "CD13", "CD14", "CD14", "CD13", "CD13", "CD14", "CD14",
        "CD15", "CD15", "CD15", "CD15", "CD13", "CD13", "CD14",
        "CD14", "CD14", "CD14", "CD13", "CD13", "CD15", "CD15",
        "CD14", "CD14", "CD14", "CD14", "CD13", "CD13", "CD14",
        "CD14", "CD13", "CD13", "CD14", "CD14", "CD15", "CD15",
        "CD15", "CD15", "CD13", "CD13", "CD14", "CD14"
    ),
    g2_Tissue = c(
        "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB",
        "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB",
        "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB",
        "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB",
        "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB",
        "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB",
        "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB",
        "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB",
        "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB",
        "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB",
        "PB", "PB", "PB", "PB"
    ),
    g2_TimePoint = c(
        3, 3, 3, 3, 1, 1, 6, 6, 1, 1, 1, 1, 2, 2, 8, 8, 3, 3, 1,
        1, 6, 6, 1, 1, 2, 2, 3, 3, 3, 3, 1, 1, 6, 6, 1, 1, 1, 1,
        2, 2, 8, 8, 3, 3, 1, 1, 6, 6, 1, 1, 2, 2, 3, 3, 3, 3, 1,
        1, 6, 6, 1, 1, 1, 1, 2, 2, 8, 8, 3, 3, 1, 1, 6, 6, 1, 1,
        2, 2, 3, 3, 3, 3, 1, 1, 6, 6, 1, 1, 1, 1, 2, 2, 8, 8, 3,
        3, 1, 1, 6, 6, 1, 1, 2, 2
    ),
    chr = c(
        NA, NA, NA, NA, "1", NA, "1", NA, "2", NA, NA, NA, NA, "3", NA,
        "3", "1", NA, NA, NA, "1", NA, NA, NA, NA, "3", NA, NA, NA, NA,
        "1", NA, "1", NA, "2", NA, NA, NA, NA, "3", NA, "3", "1", NA, NA,
        NA, "1", NA, NA, NA, NA, "3", NA, NA, "1", NA, NA, NA, NA, NA, NA,
        NA, NA, "2", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
        NA, NA, NA, NA, NA, "1", NA, "1", NA, NA, NA, NA, NA, NA, NA, NA,
        NA, "1", NA, NA, "1", "1", NA, NA, NA, NA, NA
    ),
    integration_locus = c(
        NA, NA, NA, NA, 12345, NA, 12345, NA, 54321, NA, NA,
        NA, NA, 62135, NA, 62135, 12345, NA, NA, NA, 12345,
        NA, NA, NA, NA, 62135, NA, NA, NA, NA, 12345, NA,
        12345, NA, 54321, NA, NA, NA, NA, 62135, NA, 62135,
        12345, NA, NA, NA, 12345, NA, NA, NA, NA, 62135, NA,
        NA, 43524, NA, NA, NA, NA, NA, NA, NA, NA, 76835, NA,
        NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
        NA, NA, NA, NA, 12345, NA, 12345, NA, NA, NA, NA, NA,
        NA, NA, NA, NA, 12345, NA, NA, 56832, 12345, NA, NA,
        NA, NA, NA
    ),
    strand = c(
        NA, NA, NA, NA, "+", NA, "+", NA, "-", NA, NA, NA, NA, "+", NA,
        "+", "+", NA, NA, NA, "+", NA, NA, NA, NA, "+", NA, NA, NA, NA,
        "+", NA, "+", NA, "-", NA, NA, NA, NA, "+", NA, "+", "+", NA,
        NA, NA, "+", NA, NA, NA, NA, "+", NA, NA, "+", NA, NA, NA, NA,
        NA, NA, NA, NA, "-", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
        NA, NA, NA, NA, NA, NA, NA, "+", NA, "+", NA, NA, NA, NA, NA, NA,
        NA, NA, NA, "+", NA, NA, "-", "+", NA, NA, NA, NA, NA
    ),
    sharing_perc = c(
        0, 0, 0, 0, 50, 0, 50, 0, 100, 0, 0, 0, 0, 100, 0, 50,
        100, 0, 0, 0, 100, 0, 0, 0, 0, 100, 0, 0, 0, 0, 50, 0, 50,
        0, 100, 0, 0, 0, 0, 100, 0, 50, 100, 0, 0, 0, 100, 0, 0,
        0, 0, 100, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0, 50,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0, 100, 100, 0, 0, 0,
        0, 0
    )
) |> dplyr::arrange(.data$g1, .data$g2)

test_that("iss_source produces expected output - per patient", {
    withr::local_options(list(
        ISAnalytics.parallel_processing = FALSE
    ))
    res <- iss_source(test_df2, test_df3) |>
        purrr::map(~ dplyr::arrange(.x, .data$g1, .data$g2))
    expect_equal(res, expected_res1)
})

test_that("iss_source produces expected output - NO per patient", {
    withr::local_options(list(
        ISAnalytics.parallel_processing = FALSE
    ))
    res <- iss_source(test_df2, test_df3, by_subject = FALSE) |>
        dplyr::arrange(.data$g1, .data$g2)
    expect_equal(res, expected_res2)
})


#------------------------------------------------------------------------------#
# Test gene_frequency_fisher
#------------------------------------------------------------------------------#
test_that("gene_frequency_fisher produces expected output", {
    withr::local_options(list(ISAnalytics.verbose = TRUE))
    test_cis <- readRDS(fs::path(testdata_path, "test_cis_for_fisher.Rds"))
    ## -- Testing intersection
    fisher_df <- gene_frequency_fisher(test_cis[[1]], test_cis[[2]],
        min_is_per_gene = 1,
        gene_set_method = "intersection"
    )
    common_genes <- intersect(test_cis$PT001$GeneName, test_cis$PT002$GeneName)
    expect_true(all(fisher_df$GeneName %in% common_genes) &
        all(common_genes %in% fisher_df$GeneName))
    ## -- Testing union
    fisher_df <- gene_frequency_fisher(test_cis[[1]], test_cis[[2]],
        min_is_per_gene = 1,
        gene_set_method = "union",
        remove_unbalanced_0 = FALSE
    )
    union_genes <- union(test_cis$PT001$GeneName, test_cis$PT002$GeneName)
    expect_true(all(fisher_df$GeneName %in% union_genes) &
        all(union_genes %in% fisher_df$GeneName))
    ## -- Testing default filtering
    expect_message(
        {
            fisher_df <- gene_frequency_fisher(test_cis[[1]], test_cis[[2]],
                min_is_per_gene = 3,
                gene_set_method = "intersection"
            )
        },
        class = "empty_df_gene_freq"
    )
    expect_null(fisher_df)
    fisher_df <- gene_frequency_fisher(test_cis[[1]], test_cis[[2]],
        min_is_per_gene = 3,
        gene_set_method = "union"
    )
    expect_true(nrow(fisher_df) == 2)
})

#------------------------------------------------------------------------------#
# Test top_targeted_genes
#------------------------------------------------------------------------------#
test_that("top_targeted_genes produces expected output", {
    annot_ex <- minimal_agg_example |>
        dplyr::mutate(
            GeneName = c(
                "GENE1", "GENE1", "GENE1", "GENE2", "GENE2", "GENE3", "GENE4",
                "GENE4", "GENE4"
            ),
            GeneStrand = c("+", "+", "+", "-", "-", "-", "+", "+", "+"),
            .after = "strand"
        )

    # --- top genes overall
    top_5 <- top_targeted_genes(annot_ex, n = 5, key = NULL)
    expected <- tibble::tibble(
        GeneName = c("GENE1", "GENE2", "GENE1", "GENE3", "GENE4"),
        GeneStrand = c("+", "-", "+", "-", "+"),
        chr = c("2", "5", "1", "3", "10"),
        n_IS = c(2, 2, 1, 1, 1)
    )
    expect_equal(top_5, expected)
    top_5 <- top_targeted_genes(annot_ex,
        n = 5, key = NULL,
        consider_chr = FALSE
    )
    expected <- tibble::tibble(
        GeneName = c("GENE1", "GENE4", "GENE2", "GENE3"),
        GeneStrand = c("+", "+", "-", "-"),
        n_IS = c(3, 3, 2, 1)
    )
    expect_equal(top_5, expected)
    top_5 <- top_targeted_genes(annot_ex,
        n = 5, key = NULL,
        consider_chr = FALSE,
        consider_gene_strand = FALSE
    )
    expected <- tibble::tibble(
        GeneName = c("GENE1", "GENE4", "GENE2", "GENE3"),
        n_IS = c(3, 3, 2, 1)
    )
    expect_equal(top_5, expected)

    # --- top genes per group
    top_5 <- top_targeted_genes(annot_ex,
        n = 5,
        key = c("SubjectID", "CellMarker")
    )
    expected <- tibble::tibble(
        SubjectID = c("S1", "S1", "S1", "S2", "S2", "S2", "S2"),
        CellMarker = c("CM1", "CM1", "CM1", "CM2", "CM2", "CM2", "CM2"),
        GeneName = c(
            "GENE1", "GENE2", "GENE1", "GENE3", "GENE4", "GENE4",
            "GENE4"
        ),
        GeneStrand = c("+", "-", "+", "-", "+", "+", "+"),
        chr = c("2", "5", "1", "3", "3", "10", "11"),
        n_IS = c(2, 2, 1, 1, 1, 1, 1)
    )
    top_5 <- top_targeted_genes(annot_ex,
        n = 5,
        key = c("SubjectID", "CellMarker"),
        consider_chr = FALSE
    )
    expected <- tibble::tibble(
        SubjectID = c("S1", "S1", "S2", "S2"),
        CellMarker = c("CM1", "CM1", "CM2", "CM2"),
        GeneName = c("GENE1", "GENE2", "GENE4", "GENE3"),
        GeneStrand = c("+", "-", "+", "-"),
        n_IS = c(3, 2, 3, 1)
    )
})
