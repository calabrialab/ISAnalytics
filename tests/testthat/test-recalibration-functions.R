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
  withr::local_options(list(ISAnalytics.mandatory_is_vars = "default"))
    result <- .sliding_window(
        x = sample_group1, threshold = 4,
        keep_criteria = "keep_first", annotated = FALSE,
        num_cols = "Value", max_val_col = "Value",
        sample_col = "CompleteAmplificationID",
        req_tags = mandatory_IS_vars(TRUE),
        add_col_lambdas = NULL,
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
        sample_col = "CompleteAmplificationID",
        req_tags = mandatory_IS_vars(TRUE),
        add_col_lambdas = NULL,
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
  withr::local_options(list(ISAnalytics.mandatory_is_vars = "default"))
    result <- .sliding_window(
        x = sample_group2, threshold = 4,
        keep_criteria = "keep_first", annotated = FALSE,
        num_cols = "Value", max_val_col = "Value",
        sample_col = "CompleteAmplificationID",
        req_tags = mandatory_IS_vars(TRUE),
        add_col_lambdas = NULL,
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
        sample_col = "CompleteAmplificationID",
        req_tags = mandatory_IS_vars(TRUE),
        add_col_lambdas = NULL,
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
  withr::local_options(list(ISAnalytics.mandatory_is_vars = "default"))
    result <- .sliding_window(
        x = sample_group_mult1, threshold = 4,
        keep_criteria = "keep_first", annotated = FALSE,
        num_cols = c("seqCount", "fragmentEstimate"),
        max_val_col = "seqCount",
        sample_col = "CompleteAmplificationID",
        req_tags = mandatory_IS_vars(TRUE),
        add_col_lambdas = NULL,
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
        sample_col = "CompleteAmplificationID",
        req_tags = mandatory_IS_vars(TRUE),
        add_col_lambdas = NULL,
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
  withr::local_options(list(ISAnalytics.mandatory_is_vars = "default"))
    result <- .sliding_window(
        x = sample_group_mult2, threshold = 4,
        keep_criteria = "keep_first", annotated = FALSE,
        num_cols = c("seqCount", "fragmentEstimate"),
        max_val_col = "seqCount",
        sample_col = "CompleteAmplificationID",
        req_tags = mandatory_IS_vars(TRUE),
        add_col_lambdas = NULL,
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
        sample_col = "CompleteAmplificationID",
        req_tags = mandatory_IS_vars(TRUE),
        add_col_lambdas = NULL,
        produce_map = TRUE
    )
    expect_equal(result$recalibrated_matrix, expected_for_smplmult2_mv,
        ignore_attr = TRUE
    )
    expect_equal(result$map, recalibr_map_smpl2_mv,
        ignore_attr = TRUE
    )
})

test_that(".sliding_window works on annotated", {
  withr::local_options(list(ISAnalytics.mandatory_is_vars = "default"))
  annot_group1 <- sample_group1 %>%
    dplyr::mutate(GeneName = paste0("GENE", seq_len(nrow(sample_group1))),
                  GeneStrand = "+", .after = "strand")
  expected_matrix <- expected_for_smpl1_kf %>%
    dplyr::mutate(GeneName = c("GENE4", "GENE4", "GENE4", "GENE4",
                               "GENE5", "GENE5"),
                  GeneStrand = "+", .after = "strand")
  result <- .sliding_window(
    x = annot_group1, threshold = 4,
    keep_criteria = "max_value", annotated = TRUE,
    num_cols = "Value",
    max_val_col = "Value",
    sample_col = "CompleteAmplificationID",
    req_tags = mandatory_IS_vars(TRUE),
    add_col_lambdas = NULL,
    produce_map = TRUE
  )
  expect_equal(result$recalibrated_matrix, expected_matrix,
               ignore_attr = TRUE
  )
  expect_equal(result$map, recalibr_map_smpl1_kf,
               ignore_attr = TRUE
  )
  annot_group2 <- sample_group2 %>%
    dplyr::mutate(GeneName = paste0("GENE", seq_len(nrow(sample_group2))),
                  GeneStrand = "-", .after = "strand")
  expected_matrix <- expected_for_smpl2_mv %>%
    dplyr::mutate(GeneName = c("GENE1", "GENE7", "GENE6", "GENE6", "GENE3"),
                  GeneStrand = "-", .after = "strand")
  result <- .sliding_window(
    x = annot_group2, threshold = 4,
    keep_criteria = "max_value", annotated = TRUE,
    num_cols = "Value",
    max_val_col = "Value",
    sample_col = "CompleteAmplificationID",
    req_tags = mandatory_IS_vars(TRUE),
    add_col_lambdas = NULL,
    produce_map = TRUE
  )
  expect_equal(result$recalibrated_matrix, expected_matrix,
               ignore_attr = TRUE
  )
  expect_equal(result$map, recalibr_map_smpl2_mv,
               ignore_attr = TRUE
  )
})

test_that(".sliding_window aggregates add columns correctly", {
  withr::local_options(list(ISAnalytics.mandatory_is_vars = "default"))
  s2_add_cols <- sample_group2 %>%
    dplyr::mutate(ann1 = seq_len(nrow(sample_group2)),
                  ann2 = c("a", "b", "c", "d", "e", "f", "g"),
                  ann3 = c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE),
                  .after = "strand")
  result <- .sliding_window(
    x = s2_add_cols, threshold = 4,
    keep_criteria = "max_value", annotated = FALSE,
    num_cols = "Value",
    max_val_col = "Value",
    sample_col = "CompleteAmplificationID",
    req_tags = mandatory_IS_vars(TRUE),
    add_col_lambdas = list(ann1 = ~sum(.x, na.rm = TRUE),
                           ann2 = ~paste0(.x, collapse = ";"),
                           ann3 = ~any(.x)),
    produce_map = TRUE
  )
  expected_v1 <- expected_for_smpl2_mv %>%
    dplyr::mutate(ann1 = c(6, 7, 6, 6, 3),
                  ann2 = c("e;a", "g", "d;b", "f", "c"),
                  ann3 = c(TRUE, FALSE, TRUE, TRUE, FALSE),
                  .after = "strand")
  expect_equal(result$recalibrated_matrix, expected_v1,
               ignore_attr = TRUE
  )
  expect_equal(result$map, recalibr_map_smpl2_mv,
               ignore_attr = TRUE
  )
  ## NULL lambdas correspond to keep any value
  result <- .sliding_window(
    x = s2_add_cols, threshold = 4,
    keep_criteria = "max_value", annotated = FALSE,
    num_cols = "Value",
    max_val_col = "Value",
    sample_col = "CompleteAmplificationID",
    req_tags = mandatory_IS_vars(TRUE),
    add_col_lambdas = list(ann1 = NULL,
                           ann2 = NULL,
                           ann3 = NULL),
    produce_map = TRUE
  )
  expected_v2 <- expected_for_smpl2_mv %>%
    dplyr::mutate(ann1 = c(5, 7, 4, 6, 3),
                  ann2 = c("e", "g", "d", "f", "c"),
                  ann3 = c(FALSE, FALSE, TRUE, TRUE, FALSE),
                  .after = "strand")
  expect_equal(result$recalibrated_matrix, expected_v2,
               ignore_attr = TRUE
  )
  expect_equal(result$map, recalibr_map_smpl2_mv,
               ignore_attr = TRUE
  )
})

test_that(".sliding_window works with custom vars", {
  withr::local_options(list(ISAnalytics.mandatory_is_vars = "default"))
  customized <- sample_group2 %>%
    dplyr::rename(chrom = "chr", locus = "integration_locus")
  temp_vars <- mandatory_IS_vars(TRUE)
  temp_vars[1, ]$names <- "chrom"
  temp_vars[2, ]$names <- "locus"
  withr::with_options(list(ISAnalytics.mandatory_is_vars = temp_vars),
                      {
                        result <- .sliding_window(
                          x = customized, threshold = 4,
                          keep_criteria = "max_value", annotated = FALSE,
                          num_cols = "Value", max_val_col = "Value",
                          sample_col = "CompleteAmplificationID",
                          req_tags = mandatory_IS_vars(TRUE),
                          add_col_lambdas = NULL,
                          produce_map = TRUE
                        )
                      })
  expected_matrix <- expected_for_smpl2_mv %>%
    dplyr::rename(chrom = "chr", locus = "integration_locus")
  expect_equal(result$recalibrated_matrix, expected_matrix,
               ignore_attr = TRUE
  )
  expect_equal(result$map, recalibr_map_smpl2_mv %>%
                 dplyr::rename(chrom_before = "chr_before",
                               locus_before = "integration_locus_before",
                               chrom_after = "chr_after",
                               locus_after = "integration_locus_after"),
               ignore_attr = TRUE
  )
})

#------------------------------------------------------------------------------#
# Tests .write_recalibr_map
#------------------------------------------------------------------------------#
test_that(".write_recalibr_map works if path provided is dir", {
  tmp_dir <- fs::path(tempdir(), "ISAtest")
  ## Works if folder doesn't exist (creates it and writes the file in it)
  withr::with_file(tmp_dir, {
    .write_recalibr_map(recalibr_map_smpl1_kf, tmp_dir)
    expect_true(fs::dir_exists(tmp_dir))
    expect_true(length(fs::dir_ls(tmp_dir)) == 1)
  })
  ## Works if folder already exists
  withr::with_file(tmp_dir, {
    fs::dir_create(tmp_dir)
    .write_recalibr_map(recalibr_map_smpl1_kf, tmp_dir)
    expect_true(fs::dir_exists(tmp_dir))
    expect_true(length(fs::dir_ls(tmp_dir)) == 1)
  })
})

test_that(".write_recalibr_map works if path provided is file", {
  ## Works for accepted formats
  for (ext in c("tsv", "csv", "txt",
                      paste("tsv", .compressed_formats(), sep = "."),
                      paste("csv", .compressed_formats(), sep = "."),
                      paste("txt", .compressed_formats(), sep = "."))) {
    file_name <- paste0("recalibration_map.", ext)
    tmp_file <- fs::path(tempdir(), file_name)
    withr::with_file(tmp_file, {
      .write_recalibr_map(recalibr_map_smpl1_kf, tmp_file)
      expect_true(fs::file_exists(tmp_file))
    })
  }
  ## Changes extension for unsupported extension
  file_name <- "recalibration_map.xslx"
  expected_filename <- "recalibration_map.tsv.gz"
  tmp_file <- fs::path(tempdir(), file_name)
  withr::with_options(list(ISAnalytics.verbose = TRUE), {
    expect_message({
      expect_message({
        withr::with_file(tmp_file, {
          .write_recalibr_map(recalibr_map_smpl1_kf, tmp_file)
          expect_false(fs::file_exists(tmp_file))
          expect_true(fs::file_exists(fs::path(tempdir(), expected_filename)))
        })
      }, class = "rec_unsupp_ext")
    })
  })
  file_name <- "recalibration_map.gz"
  tmp_file <- fs::path(tempdir(), file_name)
  withr::with_options(list(ISAnalytics.verbose = TRUE), {
    expect_message({
    expect_message({
      withr::with_file(tmp_file, {
        .write_recalibr_map(recalibr_map_smpl1_kf, tmp_file)
        expect_false(fs::file_exists(tmp_file))
        expect_true(fs::file_exists(fs::path(tempdir(), expected_filename)))
      })
    }, class = "rec_unsupp_ext")
    })
  })
})

#------------------------------------------------------------------------------#
# Tests compute_near_integrations
#------------------------------------------------------------------------------#
test_that("compute_near_integrations produces correct output for total", {
  withr::local_options(list(ISAnalytics.mandatory_is_vars = "default"))
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

    res <- compute_near_integrations(
      x = total_mult,
      keep_criteria = "keep_first",
      max_value_column = "seqCount",
      value_columns = c("seqCount", "fragmentEstimate"),
      is_identity_tags = NULL,
      map_as_file = FALSE
    )
    expected <- data.table::data.table(
      chr = c('1', '1', '1', '1', '1', '1', '3', '3', '3', '3', '3'),
      integration_locus = c(14571, 14571, 14571, 14571, 14577, 14577,
                            16378, 16387, 16395, 16395, 16402),
      strand = c('+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+'),
      CompleteAmplificationID = c('ID4', 'ID1', 'ID2', 'ID3', 'ID5', 'ID6',
                                  'ID1', 'ID1', 'ID4', 'ID3', 'ID5'),
      seqCount = c(1400, 70, 150, 120, 36, 15, 1969, 48, 153, 886, 543),
      fragmentEstimate = c(64.545, 83.4, 125.9, 1.656, 6.4, 564.69, 989.04,
                           1.36, 5521.1, 453.6, 12.55))
    expect_equal(res, expected)
})

test_that("compute_near_integrations warns deprecation", {
  expect_deprecated({
    res <- compute_near_integrations(
      x = sample_group1,
      keep_criteria = "keep_first",
      max_value_column = "Value",
      value_columns = "Value",
      strand_specific = TRUE,
      map_as_file = FALSE
    )
  })
})

test_that("compute_near_integrations works for package examples", {
  tmp_dir <- withr::local_tempdir()
  test_with_fine <- sample_group_mult1 %>%
    dplyr::bind_rows(sample_group_mult2) %>%
    tibble::add_case(chr = "5", integration_locus = 45213, strand = "-",
                     CompleteAmplificationID = "ID1",
                     seqCount = 45, fragmentEstimate = 56.45)

  recalibr <- compute_near_integrations(test_with_fine,
                                        map_as_file = TRUE,
                                        file_path = tmp_dir)
  expect_true(nrow(recalibr) == 12 & ncol(recalibr) == 6)
  expect_true(fs::file_exists(fs::path(tmp_dir, .generate_rec_map_filename())))
})
