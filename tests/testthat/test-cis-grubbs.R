
withr::local_options(list(ISAnalytics.verbose = FALSE))
agg <- aggregate_values_by_key(integration_matrices, association_file,
                               value_cols = c("seqCount", "fragmentEstimate"))

#------------------------------------------------------------------------------#
# Test CIS_grubbs
#------------------------------------------------------------------------------#
test_that("CIS_grubbs produces correct df", {
  output_cols <- c(
    "GeneName", "GeneStrand", "chr", "n_IS_perGene",
    "min_bp_integration_locus", "max_bp_integration_locus",
    "IS_span_bp", "avg_bp_integration_locus",
    "median_bp_integration_locus", "distinct_orientations",
    "vars", "n", "mean", "sd",
    "median", "trimmed", "mad",
    "min", "max", "range",
    "skew", "kurtosis", "se",
    "average_TxLen",
    "raw_gene_integration_frequency",
    "integration_frequency_withtolerance",
    "minus_log2_integration_freq_withtolerance",
    "zscore_minus_log2_int_freq_tolerance",
    "neg_zscore_minus_log2_int_freq_tolerance", "t_z_mlif",
    "tdist2t", "tdist_pt", "tdist_bonferroni_default",
    "tdist_bonferroni", "tdist_fdr", "tdist_benjamini",
    "tdist_positive_and_corrected", "tdist_positive",
    "tdist_positive_and_correctedEM"
  )
  result <- CIS_grubbs(agg)
  expect_true(ncol(result$cis) == 39 & nrow(result$cis) == 501)
  expect_true(all(output_cols %in% colnames(result$cis)))
  expect_equal(nrow(result$missing_genes), 5)
  ## -- Group
  result <- CIS_grubbs(agg, by = "SubjectID")
  expect_true(is.list(result$cis) & length(result$cis) == 2)
  expect_true(nrow(result$cis$PT001) == 250 & ncol(result$cis$PT001) == 39)
  expect_true(nrow(result$cis$PT002) == 292 & ncol(result$cis$PT002) == 39)
  expect_equal(nrow(result$missing_genes), 5)

  result <- CIS_grubbs(agg, by = "SubjectID", results_as_list = FALSE)
  expect_true(ncol(result$cis) == 40 & nrow(result$cis) == 542)
  expect_true("group" %in% colnames(result$cis))
  expect_equal(nrow(result$missing_genes), 5)
})

test_that("CIS_grubbs errors for char non-standard refgene", {
  expect_error({
    result <- CIS_grubbs(agg, genomic_annotation_file = "my_file.tsv")
  }, class = "genomic_file_char")
})

test_that("CIS_grubbs warns missing genes if verbose", {
  withr::local_options(list(ISAnalytics.verbose = TRUE))
  expect_message({
    result <- CIS_grubbs(agg)
  }, class = "warn_miss_genes")
})
