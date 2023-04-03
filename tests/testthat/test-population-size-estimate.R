withr::local_options(list(ISAnalytics.verbose = FALSE))
#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
test_meta <- readRDS(fs::path(testdata_path, "test_population_meta.Rds"))

test_data <- readRDS(fs::path(testdata_path, "test_population_iss.Rds"))

test_expected <- readRDS(fs::path(testdata_path,
                                  "test_population_expected.Rds"))

test_filter_expected <- readRDS(
  fs::path(testdata_path, "test_population_filter_df_expected.Rds")
)

meta_with_lineage <- test_meta |>
    dplyr::filter(.data$NumIS > 5) |>
    dplyr::left_join(blood_lineages_default(), by = "CellMarker")

#------------------------------------------------------------------------------#
# .match_celltypes_tissues
#------------------------------------------------------------------------------#
test_that(".match_celltypes_tissues works as expected", {
    # Error for incompatible lengths
    cell_types <- c("MYELOID", "B", "T")
    tissues <- c("PB", "PB")
    expect_error(
        {
            .match_celltypes_tissues(cell_types, tissues)
        },
        class = "pop_size_err_length"
    )
    # Works if either is length 1
    cell_types <- c("MYELOID", "B", "T")
    tissues <- c("PB")
    matched <- .match_celltypes_tissues(cell_types, tissues)
    expect_true(length(matched) == 3)
    expect_equal(matched[[1]], c("MYELOID", "PB"))
    expect_equal(matched[[2]], c("B", "PB"))
    expect_equal(matched[[3]], c("T", "PB"))

    cell_types <- c("MYELOID")
    tissues <- c("PB", "BM")
    matched <- .match_celltypes_tissues(cell_types, tissues)
    expect_true(length(matched) == 2)
    expect_equal(matched[[1]], c("MYELOID", "PB"))
    expect_equal(matched[[2]], c("MYELOID", "BM"))

    cell_types <- c("MYELOID")
    tissues <- c("PB")
    matched <- .match_celltypes_tissues(cell_types, tissues)
    expect_true(length(matched) == 1)
    expect_equal(matched[[1]], c("MYELOID", "PB"))
})

test_that(".match_celltypes_tissues drops duplicated groups", {
    cell_types <- c("MYELOID", "MYELOID", "MYELOID")
    tissues <- c("PB", "PB", "PB")
    matched <- .match_celltypes_tissues(cell_types, tissues)
    expect_true(length(matched) == 1)
    expect_equal(matched[[1]], c("MYELOID", "PB"))
    cell_types <- c("MYELOID", "MYELOID", "B", "T", "B", "T")
    tissues <- c("PB", "PB", "PB", "PB", "PB", "PB")
    matched <- .match_celltypes_tissues(cell_types, tissues)
    expect_true(length(matched) == 3)
    expect_equal(matched[[1]], c("MYELOID", "PB"))
    expect_equal(matched[[2]], c("B", "PB"))
    expect_equal(matched[[3]], c("T", "PB"))
})


#------------------------------------------------------------------------------#
# .re_agg_and_filter
#------------------------------------------------------------------------------#
test_that(".re_agg_and_filter works as expected", {
    re_agg <- .re_agg_and_filter(
        test_data,
        metadata = meta_with_lineage,
        fragmentEstimate_column = NULL,
        seqCount_column = "seqCount_sum",
        tissue_col = "Tissue",
        timepoint_column = "TimePoint",
        aggregation_key = c(
            "SubjectID",
            "CellMarker",
            "Tissue",
            "TimePoint"
        ),
        seqCount_threshold = 3,
        fragmentEstimate_threshold = 3,
        groups_to_proc = list(
            c("MYELOID", "PB"),
            c("B", "PB"),
            c("T", "PB")
        ),
        annotation_cols = c(
            "GeneName",
            "GeneStrand"
        ),
        subj_col = "SubjectID"
    )
    expect_equal(
        re_agg |> dplyr::arrange(.data$chr, .data$integration_locus),
        test_filter_expected[[1]] |>
            dplyr::arrange(.data$chr, .data$integration_locus)
    )
    re_agg <- .re_agg_and_filter(
        test_data,
        metadata = meta_with_lineage,
        fragmentEstimate_column = "fragmentEstimate_sum",
        seqCount_column = "seqCount_sum",
        tissue_col = "Tissue",
        timepoint_column = "TimePoint",
        aggregation_key = c(
            "SubjectID",
            "CellMarker",
            "Tissue",
            "TimePoint"
        ),
        seqCount_threshold = 3,
        fragmentEstimate_threshold = 3,
        groups_to_proc = list(
            c("MYELOID", "PB"),
            c("B", "PB"),
            c("T", "PB")
        ),
        annotation_cols = c(
            "GeneName",
            "GeneStrand"
        ),
        subj_col = "SubjectID"
    )
    expect_equal(
        re_agg |> dplyr::arrange(.data$chr, .data$integration_locus),
        test_filter_expected[[2]] |>
            dplyr::arrange(.data$chr, .data$integration_locus)
    )
})

#------------------------------------------------------------------------------#
# .convert_to_captures_matrix
#------------------------------------------------------------------------------#
test_that(".convert_to_captures_matrix gives correct column order", {
    ex_df <- tibble::tribble(
        ~chr, ~integration_locus, ~strand, ~GeneName, ~GeneStrand, ~CellType,
        ~Tissue, ~TimePoint, ~seqCount_sum,
        "1", 1034971, "-", "C1orf159", "-", "Myeloid", "BM", 60, 3,
        "1", 1034971, "-", "C1orf159", "-", "Myeloid", "BM", 140, 3,
        "1", 1034971, "-", "C1orf159", "-", "Myeloid", "BM", 360, 3,
    )
    captures_matrix <- .convert_to_captures_matrix(
        ex_df,
        quant_cols = "seqCount_sum", tissue_col = "Tissue",
        timepoint_column = "TimePoint", annotation_cols = annotation_IS_vars()
    )
    expect_equal(
        colnames(captures_matrix),
        c("Myeloid_BM_60", "Myeloid_BM_140", "Myeloid_BM_360")
    )
})

#------------------------------------------------------------------------------#
# .closed_m0_est
#------------------------------------------------------------------------------#
test_that(".closed_m0_est works as expected", {
    skip_if_not_installed("Rcapture")
    data_selection <- test_filter_expected[[2]] |>
        dplyr::filter(.data$CellType == "Myeloid", .data$Tissue == "PB")
    data_captures <- .convert_to_captures_matrix(
        data_selection,
        quant_cols = c("seqCount_sum", "fragmentEstimate_sum"),
        tissue_col = "Tissue",
        timepoint_column = "TimePoint",
        annotation_cols = annotation_IS_vars()
    )
    estimates_df <- .closed_m0_est(data_captures,
        timecaptures = ncol(data_captures),
        cols_estimate_mcm = c(
            "Model", "abundance",
            "stderr"
        ),
        subject = "PATIENT01", stable = FALSE
    )
    expected_df <- tibble::tibble(
      Model = c('M0', 'Mh Chao (LB)', 'Mh Poisson2', 'Mh Darroch',
                'Mh Gamma3.5'),
      abundance = c(3728.42791320413, 3812.05951417004,
                    3791.87510620205, 3852.04875115893,
                    3919.18254958108),
      stderr = c(16.9534089361061, 26.382222136728, 24.6483695631177,
                 35.820450572057, 51.6639000626023),
      SubjectID = c('PATIENT01', 'PATIENT01', 'PATIENT01', 'PATIENT01',
                    'PATIENT01'),
      Timepoints = c('All', 'All', 'All', 'All', 'All'),
      CellType = c('Myeloid', 'Myeloid', 'Myeloid', 'Myeloid',
                   'Myeloid'),
      Tissue = c('PB', 'PB', 'PB', 'PB', 'PB'),
      TimePoint_from = c(6, 6, 6, 6, 6),
      TimePoint_to = c(24, 24, 24, 24, 24),
      Timepoints_included = c(4, 4, 4, 4, 4),
      ModelType = c('ClosedPopulation', 'ClosedPopulation',
                    'ClosedPopulation', 'ClosedPopulation',
                    'ClosedPopulation'),
      ModelSetUp = c('models0', 'models0', 'models0', 'models0',
                     'models0'))
    expect_equal(estimates_df, expected_df)
})

test_that(".closed_m0_est works with 2 tps", {
    skip_if_not_installed("Rcapture")
    data_selection <- test_filter_expected[[2]] |>
        dplyr::filter(.data$CellType == "Myeloid", .data$Tissue == "PB")
    data_captures <- .convert_to_captures_matrix(
        data_selection,
        quant_cols = c("seqCount_sum", "fragmentEstimate_sum"),
        tissue_col = "Tissue",
        timepoint_column = "TimePoint",
        annotation_cols = annotation_IS_vars()
    )
    data_captures <- data_captures[, c(1, 2)]
    estimates_df <- .closed_m0_est(data_captures,
        timecaptures = ncol(data_captures),
        cols_estimate_mcm = c(
            "Model", "abundance",
            "stderr"
        ),
        subject = "PATIENT01", stable = FALSE
    )
    expected_df <- tibble::tibble(Model = c('M0'),
                                  abundance = c(4237.52083333333),
                                  stderr = c(138.781572395047),
                                  SubjectID = c('PATIENT01'),
                                  Timepoints = c('All'),
                                  CellType = c('Myeloid'),
                                  Tissue = c('PB'),
                                  TimePoint_from = c(6),
                                  TimePoint_to = c(12),
                                  Timepoints_included = c(2),
                                  ModelType = c('ClosedPopulation'),
                                  ModelSetUp = c('models0'))
    expect_equal(estimates_df, expected_df)
})

#------------------------------------------------------------------------------#
# .closed_mthchaobc_est
#------------------------------------------------------------------------------#
test_that(".closed_mthchaobc_est works as expected", {
    skip_if_not_installed("Rcapture")
    data_selection <- test_filter_expected[[2]] |>
        dplyr::filter(.data$CellType == "Myeloid", .data$Tissue == "PB")
    data_captures <- .convert_to_captures_matrix(
        data_selection,
        quant_cols = c("seqCount_sum", "fragmentEstimate_sum"),
        tissue_col = "Tissue",
        timepoint_column = "TimePoint",
        annotation_cols = annotation_IS_vars()
    )
    estimates_df <- .closed_mthchaobc_est(data_captures,
        timecaptures = ncol(data_captures),
        cols_estimate_mcm = c(
            "Model", "abundance",
            "stderr"
        ),
        subject = "PATIENT01", stable = FALSE
    )
    expected_df <- tibble::tibble(
      Model = c('M0', 'Mt', 'Mth Chao (LB)'),
      abundance = c(3728.39194399291, 3611.286068073, 3703.79712337948),
      stderr = c(16.9512894097357, 10.6507907639934, 19.058670892136),
      SubjectID = c('PATIENT01', 'PATIENT01', 'PATIENT01'),
      Timepoints = c('All', 'All', 'All'),
      CellType = c('Myeloid', 'Myeloid', 'Myeloid'),
      Tissue = c('PB', 'PB', 'PB'),
      TimePoint_from = c(6, 6, 6),
      TimePoint_to = c(24, 24, 24),
      Timepoints_included = c(4, 4, 4),
      ModelType = c('ClosedPopulation', 'ClosedPopulation',
                    'ClosedPopulation'),
      ModelSetUp = c('mthchaobc', 'mthchaobc', 'mthchaobc'))
    expect_equal(estimates_df, expected_df)
})

test_that(".closed_mthchaobc_est works with 2 tps", {
    skip_if_not_installed("Rcapture")
    data_selection <- test_filter_expected[[2]] |>
        dplyr::filter(.data$CellType == "Myeloid", .data$Tissue == "PB")
    data_captures <- .convert_to_captures_matrix(
        data_selection,
        quant_cols = c("seqCount_sum", "fragmentEstimate_sum"),
        tissue_col = "Tissue",
        timepoint_column = "TimePoint",
        annotation_cols = annotation_IS_vars()
    )
    data_captures <- data_captures[, c(1, 2)]
    estimates_df <- .closed_mthchaobc_est(
        data_captures,
        timecaptures = ncol(data_captures),
        cols_estimate_mcm = c(
            "Model", "abundance",
            "stderr"
        ),
        subject = "PATIENT01", stable = FALSE
    )
    expected_df <- tibble::tibble(
      Model = c('M0', 'Mt'),
      abundance = c(4231.92277713625, 3019.07159353349),
      stderr = c(138.340184265819, 72.0050170227728),
      SubjectID = c('PATIENT01', 'PATIENT01'),
      Timepoints = c('All', 'All'),
      CellType = c('Myeloid', 'Myeloid'),
      Tissue = c('PB', 'PB'),
      TimePoint_from = c(6, 6),
      TimePoint_to = c(12, 12),
      Timepoints_included = c(2, 2),
      ModelType = c('ClosedPopulation', 'ClosedPopulation'),
      ModelSetUp = c('mthchaobc', 'mthchaobc'))
    expect_equal(estimates_df, expected_df)
})

#------------------------------------------------------------------------------#
# .consecutive_m0bc_est
#------------------------------------------------------------------------------#
test_that(".consecutive_m0bc_est works as expected", {
    skip_if_not_installed("Rcapture")
    data_selection <- test_filter_expected[[2]] |>
        dplyr::filter(.data$CellType == "Myeloid", .data$Tissue == "PB")
    data_captures <- .convert_to_captures_matrix(
        data_selection,
        quant_cols = c("seqCount_sum", "fragmentEstimate_sum"),
        tissue_col = "Tissue",
        timepoint_column = "TimePoint",
        annotation_cols = annotation_IS_vars()
    )
    estimates_df <- .consecutive_m0bc_est(
        matrix_desc = data_captures,
        cols_estimate_mcm = c(
            "Model", "abundance",
            "stderr"
        ),
        subject = "PATIENT01"
    )
    expected_df <- tibble::tibble(
      Model = c('M0', 'M0', 'M0'),
      abundance = c(4231.92277713625, 3213.31381578947, 3639.2570458528),
      stderr = c(138.340184265819, 34.7602331515537, 27.6367496532401),
      SubjectID = c('PATIENT01', 'PATIENT01', 'PATIENT01'),
      Timepoints = c('Consecutive', 'Consecutive', 'Consecutive'),
      CellType = c('Myeloid', 'Myeloid', 'Myeloid'),
      Tissue = c('PB', 'PB', 'PB'),
      TimePoint_from = c(6, 12, 18),
      TimePoint_to = c(12, 18, 24),
      Timepoints_included = c(2, 2, 2),
      ModelType = c('ClosedPopulation', 'ClosedPopulation', 'ClosedPopulation'),
      ModelSetUp = c('models0BC', 'models0BC', 'models0BC'))
    expect_equal(estimates_df, expected_df)
})

test_that(".consecutive_m0bc_est works for 2 tps", {
    skip_if_not_installed("Rcapture")
    data_selection <- test_filter_expected[[2]] |>
        dplyr::filter(.data$CellType == "Myeloid", .data$Tissue == "PB")
    data_captures <- .convert_to_captures_matrix(
        data_selection,
        quant_cols = c("seqCount_sum", "fragmentEstimate_sum"),
        tissue_col = "Tissue",
        timepoint_column = "TimePoint",
        annotation_cols = annotation_IS_vars()
    )
    data_captures <- data_captures[, c(1, 2)]
    estimates_df <- .consecutive_m0bc_est(
        matrix_desc = data_captures,
        cols_estimate_mcm = c(
            "Model", "abundance",
            "stderr"
        ),
        subject = "PATIENT01"
    )
    expected_df <- tibble::tibble(Model = c('M0'),
                                  abundance = c(4231.92277713625),
                                  stderr = c(138.340184265819),
                                  SubjectID = c('PATIENT01'),
                                  Timepoints = c('Consecutive'),
                                  CellType = c('Myeloid'),
                                  Tissue = c('PB'),
                                  TimePoint_from = c(6),
                                  TimePoint_to = c(12),
                                  Timepoints_included = c(2),
                                  ModelType = c('ClosedPopulation'),
                                  ModelSetUp = c('models0BC'))
    expect_equal(estimates_df, expected_df)
})

#------------------------------------------------------------------------------#
# .consecutive_mth_est
#------------------------------------------------------------------------------#
test_that(".consecutive_mth_est works as expected", {
    skip_if_not_installed("Rcapture")
    data_selection <- test_filter_expected[[2]] |>
        dplyr::filter(.data$CellType == "Myeloid", .data$Tissue == "PB")
    data_captures <- .convert_to_captures_matrix(
        data_selection,
        quant_cols = c("seqCount_sum", "fragmentEstimate_sum"),
        tissue_col = "Tissue",
        timepoint_column = "TimePoint",
        annotation_cols = annotation_IS_vars()
    )
    estimates_df <- .consecutive_mth_est(
        matrix_desc = data_captures,
        cols_estimate_mcm = c(
            "Model", "abundance",
            "stderr"
        ),
        subject = "PATIENT01"
    )
    expected_df <- tibble::tibble(
      Model = c('Mth Chao (LB)', 'Mth Chao (LB)'),
      abundance = c(3389.66173863073, 3642.4428335005),
      stderr = c(38.9612415429626, 18.5479954820153),
      SubjectID = c('PATIENT01', 'PATIENT01'),
      Timepoints = c('Consecutive', 'Consecutive'),
      CellType = c('Myeloid', 'Myeloid'),
      Tissue = c('PB', 'PB'),
      TimePoint_from = c(6, 12),
      TimePoint_to = c(18, 24),
      Timepoints_included = c(3, 3),
      ModelType = c('ClosedPopulation', 'ClosedPopulation'),
      ModelSetUp = c('modelMTHBC', 'modelMTHBC'))
    expect_equal(estimates_df, expected_df)
})

#------------------------------------------------------------------------------#
# .estimate_pop
#------------------------------------------------------------------------------#
test_that(".estimate_pop works as expected", {
    skip_if_not_installed("Rcapture")
    data_selection <- test_filter_expected[[2]] |>
        dplyr::filter(.data$CellType == "Myeloid", .data$Tissue == "PB")
    estimates <- .estimate_pop(
        df = data_selection,
        seqCount_column = "seqCount_sum",
        fragmentEstimate_column = "fragmentEstimate_sum",
        timepoint_column = "TimePoint",
        annotation_cols = annotation_IS_vars(),
        stable_timepoints = NULL,
        subject = "PATIENT01",
        tissue_col = "Tissue"
    )
    expected_df <- tibble::tibble(
      Model = c('M0', 'Mh Chao (LB)', 'Mh Poisson2', 'Mh Darroch',
                'Mh Gamma3.5', 'M0', 'Mt', 'Mth Chao (LB)', 'M0',
                'M0', 'M0'),
      abundance = c(3728.42791320413, 3812.05951417004,
                    3791.87510620205, 3852.04875115893,
                    3919.18254958108, 3728.39194399291,
                    3611.286068073, 3703.79712337948,
                    4231.92277713625, 3213.31381578947,
                    3639.2570458528),
      stderr = c(16.9534089361061, 26.382222136728, 24.6483695631177,
                 35.820450572057, 51.6639000626023, 16.9512894097357,
                 10.6507907639934, 19.058670892136, 138.340184265819,
                 34.7602331515537, 27.6367496532401),
      SubjectID = c('PATIENT01', 'PATIENT01', 'PATIENT01', 'PATIENT01',
                    'PATIENT01', 'PATIENT01', 'PATIENT01', 'PATIENT01',
                    'PATIENT01', 'PATIENT01', 'PATIENT01'),
      Timepoints = c('All', 'All', 'All', 'All', 'All', 'All', 'All',
                     'All', 'Consecutive', 'Consecutive', 'Consecutive'),
      CellType = c('Myeloid', 'Myeloid', 'Myeloid', 'Myeloid', 'Myeloid',
                   'Myeloid', 'Myeloid', 'Myeloid', 'Myeloid', 'Myeloid',
                   'Myeloid'),
      Tissue = c('PB', 'PB', 'PB', 'PB', 'PB', 'PB', 'PB', 'PB', 'PB',
                 'PB', 'PB'),
      TimePoint_from = c(6, 6, 6, 6, 6, 6, 6, 6, 6, 12, 18),
      TimePoint_to = c(24, 24, 24, 24, 24, 24, 24, 24, 12, 18, 24),
      Timepoints_included = c(4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2),
      ModelType = c('ClosedPopulation', 'ClosedPopulation',
                    'ClosedPopulation', 'ClosedPopulation',
                    'ClosedPopulation', 'ClosedPopulation',
                    'ClosedPopulation', 'ClosedPopulation',
                    'ClosedPopulation', 'ClosedPopulation',
                    'ClosedPopulation'),
      ModelSetUp = c('models0', 'models0', 'models0', 'models0',
                     'models0', 'mthchaobc', 'mthchaobc', 'mthchaobc',
                     'models0BC', 'models0BC', 'models0BC'))
    expect_equal(estimates$est, expected_df)
    expect_true(length(estimates$logger) == 3)

    estimates <- .estimate_pop(
        df = data_selection,
        seqCount_column = "seqCount_sum",
        fragmentEstimate_column = "fragmentEstimate_sum",
        timepoint_column = "TimePoint",
        annotation_cols = annotation_IS_vars(),
        stable_timepoints = c(12, 18, 24),
        subject = "PATIENT01",
        tissue_col = "Tissue"
    )
    expected_df <- tibble::tibble(
      Model = c('M0', 'Mh Chao (LB)', 'Mh Poisson2', 'Mh Darroch',
                'Mh Gamma3.5', 'M0', 'Mh Chao (LB)', 'Mh Poisson2',
                'Mh Darroch', 'Mh Gamma3.5', 'M0', 'Mt',
                'Mth Chao (LB)', 'M0', 'Mt', 'Mth Chao (LB)',
                'M0', 'M0', 'M0', 'Mth Chao (LB)', 'Mth Chao (LB)'),
      abundance = c(3728.42791320413, 3812.05951417004, 3791.87510620205,
                    3852.04875115893, 3919.18254958108, 3624.15261606024,
                    3738.0582174137, 3859.94234449885, 4037.04032298801,
                    4298.38860033931, 3728.39194399291, 3611.286068073,
                    3703.79712337948, 3624.09305004649, 3561.16510345987,
                    3642.4428335005, 4231.92277713625, 3213.31381578947,
                    3639.2570458528, 3389.66173863073, 3642.4428335005),
      stderr = c(16.9534089361061, 26.382222136728, 24.6483695631177,
                 35.820450572057, 51.6639000626023, 14.6949984853492,
                 24.6645237658986, 42.7923810992611, 77.030081870443,
                 137.914165992544, 16.9512894097357, 10.6507907639934,
                 19.058670892136, 14.6913550528532, 11.1577713011559,
                 18.5479954820153, 138.340184265819, 34.7602331515537,
                 27.6367496532401, 38.9612415429626, 18.5479954820153),
      SubjectID = c('PATIENT01', 'PATIENT01', 'PATIENT01', 'PATIENT01',
                    'PATIENT01', 'PATIENT01', 'PATIENT01', 'PATIENT01',
                    'PATIENT01', 'PATIENT01', 'PATIENT01', 'PATIENT01',
                    'PATIENT01', 'PATIENT01', 'PATIENT01', 'PATIENT01',
                    'PATIENT01', 'PATIENT01', 'PATIENT01', 'PATIENT01',
                    'PATIENT01'),
      Timepoints = c('All', 'All', 'All', 'All', 'All', 'Stable',
                     'Stable', 'Stable', 'Stable', 'Stable', 'All',
                     'All', 'All', 'Stable', 'Stable', 'Stable',
                     'Consecutive', 'Consecutive', 'Consecutive',
                     'Consecutive', 'Consecutive'),
      CellType = c('Myeloid', 'Myeloid', 'Myeloid', 'Myeloid', 'Myeloid',
                   'Myeloid', 'Myeloid', 'Myeloid', 'Myeloid', 'Myeloid',
                   'Myeloid', 'Myeloid', 'Myeloid', 'Myeloid', 'Myeloid',
                   'Myeloid', 'Myeloid', 'Myeloid', 'Myeloid', 'Myeloid',
                   'Myeloid'),
      Tissue = c('PB', 'PB', 'PB', 'PB', 'PB', 'PB', 'PB', 'PB', 'PB',
                 'PB', 'PB', 'PB', 'PB', 'PB', 'PB', 'PB', 'PB', 'PB',
                 'PB', 'PB', 'PB'),
      TimePoint_from = c(6, 6, 6, 6, 6, 12, 12, 12, 12, 12, 6, 6, 6, 12,
                         12, 12, 6, 12, 18, 6, 12),
      TimePoint_to = c(24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24,
                       24, 24, 24, 12, 18, 24, 18, 24),
      Timepoints_included = c(4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 4, 4, 4, 3, 3,
                              3, 2, 2, 2, 3, 3),
      ModelType = c('ClosedPopulation', 'ClosedPopulation',
                    'ClosedPopulation', 'ClosedPopulation',
                    'ClosedPopulation', 'ClosedPopulation',
                    'ClosedPopulation', 'ClosedPopulation',
                    'ClosedPopulation', 'ClosedPopulation',
                    'ClosedPopulation', 'ClosedPopulation',
                    'ClosedPopulation', 'ClosedPopulation',
                    'ClosedPopulation', 'ClosedPopulation',
                    'ClosedPopulation', 'ClosedPopulation',
                    'ClosedPopulation', 'ClosedPopulation',
                    'ClosedPopulation'),
      ModelSetUp = c('models0', 'models0', 'models0', 'models0',
                     'models0', 'models0', 'models0', 'models0',
                     'models0', 'models0', 'mthchaobc', 'mthchaobc',
                     'mthchaobc', 'mthchaobc', 'mthchaobc',
                     'mthchaobc', 'models0BC', 'models0BC',
                     'models0BC', 'modelMTHBC', 'modelMTHBC'))
    expect_equal(estimates$est, expected_df)
    expect_true(length(estimates$logger) == 0)
})

#------------------------------------------------------------------------------#
# .re_agg_and_estimate
#------------------------------------------------------------------------------#
test_that(".re_agg_and_estimate works as expected - base case", {
    # Base case: all groups of interest possess at least 2 time points
    est_slice <- .re_agg_and_estimate(
        test_data,
        metadata = meta_with_lineage,
        fragmentEstimate_column = "fragmentEstimate_sum",
        seqCount_column = "seqCount_sum",
        tissue_col = "Tissue",
        timepoint_column = "TimePoint",
        aggregation_key = c(
            "SubjectID",
            "CellMarker",
            "Tissue",
            "TimePoint"
        ),
        seqCount_threshold = 3,
        fragmentEstimate_threshold = 3,
        groups_to_proc = list(
            c("MYELOID", "PB"),
            c("B", "PB"),
            c("T", "PB")
        ),
        annotation_cols = c(
            "GeneName",
            "GeneStrand"
        ),
        subj_col = "SubjectID",
        stable_timepoints = c(12, 18, 24)
    )
    expect_equal(est_slice$est, test_expected)
    expect_true(length(est_slice$log) == 0)
})

test_that(".re_agg_and_estimate works as expected - less tps", {
    # One group of interest has only 1 time point (excluded)
    filtered_test_data <- test_data |>
        dplyr::filter(!(.data$CellMarker == "CD19" & .data$TimePoint %in% c(
            "12", "18", "24"
        ) & .data$Tissue == "PB"))
    est_slice <- .re_agg_and_estimate(
        filtered_test_data,
        metadata = meta_with_lineage,
        fragmentEstimate_column = "fragmentEstimate_sum",
        seqCount_column = "seqCount_sum",
        tissue_col = "Tissue",
        timepoint_column = "TimePoint",
        aggregation_key = c(
            "SubjectID",
            "CellMarker",
            "Tissue",
            "TimePoint"
        ),
        seqCount_threshold = 3,
        fragmentEstimate_threshold = 3,
        groups_to_proc = list(
            c("MYELOID", "PB"),
            c("B", "PB"),
            c("T", "PB")
        ),
        annotation_cols = c(
            "GeneName",
            "GeneStrand"
        ),
        subj_col = "SubjectID",
        stable_timepoints = c(12, 18, 24)
    )
    expect_equal(
        est_slice$log[[1]],
        "PATIENT01, B, PB - Has less than 2 time points, won't be processed"
    )
    expect_true(
        all(unique(est_slice$est$CellType) == c("Myeloid", "T")) &
            unique(est_slice$est$Tissue) == c("PB")
    )
})

#------------------------------------------------------------------------------#
# HSC_population_size_estimate
#------------------------------------------------------------------------------#
test_that("HSC_population_size_estimate works as expected", {
    data_two_patients <- test_data |>
        dplyr::bind_rows(
            test_data |>
                dplyr::mutate(SubjectID = "PATIENT02")
        )
    meta_with_two_patients <- test_meta |>
        dplyr::bind_rows(
            test_meta |>
                dplyr::mutate(SubjectID = "PATIENT02")
        )
    expected_final <- test_expected |>
        dplyr::mutate(
            PopSize = round(.data$abundance - .data$stderr)
        )
    expected_final <- expected_final |>
        dplyr::bind_rows(
            expected_final |>
                dplyr::mutate(SubjectID = "PATIENT02")
        )
    pop_est <- HSC_population_size_estimate(
        data_two_patients,
        metadata = meta_with_two_patients,
        stable_timepoints = c(12, 18, 24),
        aggregation_key = c("SubjectID", "CellMarker", "Tissue", "TimePoint"),
        cell_type = c("MYELOID", "B", "T"),
        tissue_type = c("PB")
    )
    expect_true(length(pop_est$log[[1]]) == 0)
    expect_true(length(pop_est$log[[2]]) == 0)
    expect_equal(pop_est$est, expected_final)
})

test_that("HSC_population_size_estimate produces output missing NumIS", {
    data_two_patients <- test_data |>
        dplyr::bind_rows(
            test_data |>
                dplyr::mutate(SubjectID = "PATIENT02")
        )
    meta_with_two_patients <- test_meta |>
        dplyr::bind_rows(
            test_meta |>
                dplyr::mutate(SubjectID = "PATIENT02")
        ) |>
        dplyr::select(-dplyr::all_of("NumIS"))
    expected_final <- test_expected |>
        dplyr::mutate(
            PopSize = round(.data$abundance - .data$stderr)
        )
    expected_final <- expected_final |>
        dplyr::bind_rows(
            expected_final |>
                dplyr::mutate(SubjectID = "PATIENT02")
        )
    pop_est <- HSC_population_size_estimate(
        data_two_patients,
        metadata = meta_with_two_patients,
        stable_timepoints = c(12, 18, 24),
        aggregation_key = c("SubjectID", "CellMarker", "Tissue", "TimePoint"),
        cell_type = c("MYELOID", "B", "T"),
        tissue_type = c("PB")
    )
    expect_true(length(pop_est$log[[1]]) == 0)
    expect_true(length(pop_est$log[[2]]) == 0)
    expect_equal(pop_est$est, expected_final)
})
