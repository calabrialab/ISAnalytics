context("Aggregate functions")

library(ISAnalytics)

#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
op <- options(ISAnalytics.widgets = FALSE)
on.exit(options(op))

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
                          package = "ISAnalytics")

root_err <- unzip_file_system(root_err, "fserr")

# Association file
association_file <- import_association_file(path_af, root_correct)
suppressWarnings({
  association_file_err <- import_association_file(path_af, root_err)
})

# Matrices
matrices <- matrices <- import_parallel_Vispa2Matrices_auto(
  association_file = path_af, root = root_correct,
  quantification_type = c("seqCount", "fragmentEstimate"),
  matrix_type = "annotated", workers = 2, patterns = NULL,
  matching_opt = "ANY"
)

#------------------------------------------------------------------------------#
# Tests .stats_report
#------------------------------------------------------------------------------#
test_that(".stats_report returns NULL if all paths are NA", {
  association_file_mod <- association_file
  association_file_mod$Path <- rep_len(NA_character_,
                                       length(association_file_mod$Path))
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
  expect_true(all(!stats_rep$files == ""))
  stats_rep <- .stats_report(association_file_err)
  expect_false(all(!stats_rep$files == ""))
  expect_true((stats_rep %>%
                 dplyr::filter(.data$ProjectID == "PROJECT1101"))$files == "")
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
  agg <- .join_and_aggregate(af_filtered, stats[[1]],
                             c("SubjectID", "CellMarker",
                               "Tissue", "TimePoint"))
  expect_equivalent(agg$VCN, c(0.7, 2.8, 4.3, 9.89))
  expect_true(all(agg$Avg_DNAngUsed == 100))
  expect_equivalent(agg$Kapa, c(26.54667, 59.81000, 79.80333, 85.87000),
                    tolerance = .5)
  expect_true(all(agg$DNAngUsed == 300))
  expect_equivalent(agg$ulForPool, c(25.7, 12.29, 7.98, 7.46))
  expect_equivalent(agg$BARCODE_MUX, c(4872448, 8573674, 5841068, 6086489))
  expect_equivalent(agg$TRIMMING_FINAL_LTRLC,
                    c(4863932, 8549992, 5818598, 6062299))
  expect_equivalent(agg$LV_MAPPED,
                    c(2114395, 3485464, 2923561, 2484034))
  expect_equivalent(agg$BWA_MAPPED_OVERALL,
                    c(2619004, 4751887, 2631197, 3297007))
  expect_equivalent(agg$ISS_MAPPED_PP,
                    c(2386173, 3746655, 1938140, 2330286))
})

test_that(".join_and_aggregate correctly join and aggregate for af only", {
  af_filtered <- association_file %>%
    dplyr::filter(.data$ProjectID == "CLOEXP")
  agg <- .join_and_aggregate(af_filtered, NULL,
                             c("SubjectID", "CellMarker",
                               "Tissue", "TimePoint"))
  expect_equivalent(agg$VCN, c(0.7, 2.8, 4.3, 9.89))
  expect_true(all(agg$Avg_DNAngUsed == 100))
  expect_equivalent(agg$Kapa, c(26.54667, 59.81000, 79.80333, 85.87000),
                    tolerance = .5)
  expect_true(all(agg$DNAngUsed == 300))
  expect_equivalent(agg$ulForPool, c(25.7, 12.29, 7.98, 7.46))
  expect_true(all(!.stats_columns_min() %in% colnames(agg)))
})

#------------------------------------------------------------------------------#
# Tests aggregate_metadata
#------------------------------------------------------------------------------#
# Test input
test_that("aggregate_metadata stops if association file is not a tibble", {
  expect_error({
    agg_meta <- aggregate_metadata(association_file = as.data.frame(
      association_file))
  })
  expect_error({
    agg_meta <- aggregate_metadata(association_file = 1)
  })
})
test_that("aggregate_metadata stops if af is missing mandatory columns", {
  association_file <- association_file %>%
    dplyr::select(-c(.data$FusionPrimerPCRDate))
  expect_error({
    agg_meta <- aggregate_metadata(association_file = association_file)
  })
})
test_that("aggregate_metadata stops if grouping keys is null", {
  expect_error({
    agg_meta <- aggregate_metadata(association_file = association_file,
                                   grouping_keys = NULL)
  })
})
test_that("aggregate_metadata stops if grouping keys is not a char vector", {
  expect_error({
    agg_meta <- aggregate_metadata(association_file = association_file,
                                   grouping_keys = 1)
  })
  expect_error({
    agg_meta <- aggregate_metadata(association_file = association_file,
                                   grouping_keys = c(1,2))
  })
})
test_that("aggregate_metadata stops if grouping keys are missing from af", {
  expect_error({
    agg_meta <- aggregate_metadata(association_file = association_file,
                                   grouping_keys = "a")
  })
  expect_error({
    agg_meta <- aggregate_metadata(association_file = association_file,
                                   grouping_keys = c("a", "ProjectID"))
  })
})

test_that("aggregate_metadata stops if import_stats is not logical", {
  expect_error({
    agg_meta <- aggregate_metadata(association_file = association_file,
                                   import_stats = 1)
  })
  expect_error({
    agg_meta <- aggregate_metadata(association_file = association_file,
                                   import_stats = c(TRUE, TRUE, FALSE))
  })
})

# Test Values
test_that("aggregate_metadata succeeds if all params are correct", {
  op <- options(ISAnalytics.verbose = FALSE)
  on.exit({
    options(op)
  })
  af_filtered <- association_file %>%
    dplyr::filter(.data$ProjectID == "CLOEXP")
  agg <- aggregate_metadata(af_filtered)
  expect_equivalent(agg$VCN, c(0.7, 2.8, 4.3, 9.89))
  expect_true(all(agg$Avg_DNAngUsed == 100))
  expect_equivalent(agg$Kapa, c(26.54667, 59.81000, 79.80333, 85.87000),
                    tolerance = .5)
  expect_true(all(agg$DNAngUsed == 300))
  expect_equivalent(agg$ulForPool, c(25.7, 12.29, 7.98, 7.46))
  expect_equivalent(agg$BARCODE_MUX, c(4872448, 8573674, 5841068, 6086489))
  expect_equivalent(agg$TRIMMING_FINAL_LTRLC,
                    c(4863932, 8549992, 5818598, 6062299))
  expect_equivalent(agg$LV_MAPPED,
                    c(2114395, 3485464, 2923561, 2484034))
  expect_equivalent(agg$BWA_MAPPED_OVERALL,
                    c(2619004, 4751887, 2631197, 3297007))
  expect_equivalent(agg$ISS_MAPPED_PP,
                    c(2386173, 3746655, 1938140, 2330286))
})

#------------------------------------------------------------------------------#
# Tests .wrapper
#------------------------------------------------------------------------------#
test_that(".wrapper returns simple vector if atomic", {
  res <- .wrapper(c(1,2,3), lambda = "sum", args = NULL, namespace = "base",
                  env = .GlobalEnv)
  expect_false(is.list(res))
  expect_true(is.numeric(res))
  expect_true(length(res) == 1)
})

test_that(".wrapper returns list if not atomic", {
  res <- .wrapper(c(1,2,3), lambda = "summary", args = NULL, namespace = "base",
                  env = .GlobalEnv)
  expect_true(is.list(res))
  expect_s3_class(res[[1]], "summaryDefault")
  expect_true(length(res) == 1)
  res <- .wrapper(c(1,2,3), lambda = "as.data.frame",
                  args = NULL, namespace = "base",
                  env = .GlobalEnv)
  expect_true(is.list(res))
  expect_true(tibble::is_tibble(res[[1]]))
  expect_true(length(res) == 1)
})

#------------------------------------------------------------------------------#
# Tests .aggregate_lambda
#------------------------------------------------------------------------------#
### OTHER VARS ###
ex_matrix <- function() {
  tibble::tibble(chr = c(rep_len(1, 8)),
                 integration_locus = c(rep_len(56382099, 4),
                                       rep_len(57531752, 4)),
                 strand = c(rep_len("-", 4),
                            rep_len("+", 4)),
                 GeneName = c(rep_len("PLPP3", 4),
                              rep_len("DAB1", 4)),
                 GeneStrand = c(rep_len("-", 8)),
                 CompleteAmplificationID = c(paste("id", 1:8, sep = "_")),
                 Value = c(rep_len(4, 4), rep_len(1, 4)))
}

ex_af <- function() {
  tibble::tibble(CompleteAmplificationID = c(paste("id", 1:8, sep = "_")),
                 SubjectID = c(rep_len("subj1", 4), rep_len("subj2", 4)))
}

example_m <- ex_matrix()
example_af <- ex_af()

test_that(".aggregate_lambda produces correct output", {
  agg <- .aggregate_lambda(x = list(example_m), af = example_af,
                           key = "SubjectID", lambda = "sum",
                           group = c("chr", "integration_locus","strand",
                                     "GeneName", "GeneStrand"),
                           args = NULL, namespace = "base", envir = .GlobalEnv)
  expect_true((agg[[1]][1,]$Aggregated_value == 16))
  expect_true((agg[[1]][2,]$Aggregated_value == 4))
})

test_that(".aggregate_lambda produces correct output - custom", {
  foo <- function(x) {
    sum(x)
  }
  env <- rlang::current_env()
  agg <- .aggregate_lambda(x = list(example_m), af = example_af,
                           key = "SubjectID", lambda = "foo",
                           group = c("chr", "integration_locus","strand",
                                     "GeneName", "GeneStrand"),
                           args = NULL, namespace = NULL, envir = env)
  expect_true((agg[[1]][1,]$Aggregated_value == 16))
  expect_true((agg[[1]][2,]$Aggregated_value == 4))
})

#------------------------------------------------------------------------------#
# Tests aggregate_values_by_key
#------------------------------------------------------------------------------#
### OTHER VARS ###
path_sm <- system.file("extdata", "ex_annotated_ISMatrix.tsv.xz",
                       package = "ISAnalytics")
isadf <- import_single_Vispa2Matrix(path_sm)

# Test input
test_that("aggregate_values_by_key stops if x is not a list or an IS matrix", {
  expect_error({
    agg <- aggregate_values_by_key(x = 1, association_file = association_file)
  })
  expect_error({
    agg <- aggregate_values_by_key(x = tibble::tibble(a = c(1,2,3)),
                                   association_file = association_file)
  }, regexp = .non_ISM_error())
  isadf_mod1 <- isadf %>% dplyr::select(-c(.data$Value))
  isadf_mod2 <- isadf %>% dplyr::mutate(Value = as.character(.data$Value))
  expect_error({
    agg <- aggregate_values_by_key(x = isadf_mod1,
                                   association_file = association_file)
  }, regexp = .missing_value_col_error())
  expect_error({
    agg <- aggregate_values_by_key(x = isadf_mod2,
                                   association_file = association_file)
  }, regexp = .missing_value_col_error())
  expect_error({
    agg <- aggregate_values_by_key(x = list(isadf_mod1, isadf),
                                   association_file = association_file)
  }, regexp = .missing_value_col_error())
})

test_that("aggregate_values_by_key stops if af is not a tibble", {
  expect_error({
    agg <- aggregate_values_by_key(x = isadf, association_file = 1)
  })
})
test_that("aggregate_values_by_key stops if key is incorrect", {
  expect_error({
    agg <- aggregate_values_by_key(x = isadf,
                                   association_file = association_file,
                                   key = 1)
  })
  expect_error({
    agg <- aggregate_values_by_key(x = isadf,
                                   association_file = association_file,
                                   key = "x")
  })
  expect_error({
    agg <- aggregate_values_by_key(x = isadf,
                                   association_file = association_file,
                                   key = c("SubjectID", "x"))
  })
})
test_that("aggregate_values_by_key stops if lambda is incorrect", {
  expect_error({
    agg <- aggregate_values_by_key(x = isadf,
                                   association_file = association_file,
                                   key = "SubjectID",
                                   lambda = 1)
  })
  expect_error({
    agg <- aggregate_values_by_key(x = isadf,
                                   association_file = association_file,
                                   key = "SubjectID",
                                   lambda = sum)
  })
  expect_error({
    agg <- aggregate_values_by_key(x = isadf,
                                   association_file = association_file,
                                   key = "SubjectID",
                                   lambda = c("sum", "mean"))
  })
})
test_that("aggregate_values_by_key stops if args is incorrect", {
  expect_error({
    agg <- aggregate_values_by_key(x = isadf,
                                   association_file = association_file,
                                   key = "SubjectID",
                                   lambda = "sum",
                                   args = c(1,2))
  })
  expect_error({
    agg <- aggregate_values_by_key(x = isadf,
                                   association_file = association_file,
                                   key = "SubjectID",
                                   lambda = "sum",
                                   args = list(1,2))
  })
  expect_error({
    agg <- aggregate_values_by_key(x = isadf,
                                   association_file = association_file,
                                   key = "SubjectID",
                                   lambda = "sum",
                                   args = NULL)
  }, regexp = NA)
})
test_that("aggregate_values_by_key stops if group is incorrect", {
  expect_error({
    agg <- aggregate_values_by_key(x = isadf,
                                   association_file = association_file,
                                   key = "SubjectID",
                                   lambda = "sum",
                                   args = NULL,
                                   group = c(1,2))
  })
  expect_error({
    agg <- aggregate_values_by_key(x = isadf,
                                   association_file = association_file,
                                   key = "SubjectID",
                                   lambda = "sum",
                                   args = NULL,
                                   group = c("xfsf", "fwre"))
  })
  expect_error({
    agg <- aggregate_values_by_key(x = isadf,
                                   association_file = association_file,
                                   key = "SubjectID",
                                   lambda = "sum",
                                   args = NULL,
                                   group = NULL)
  }, regexp = NA)
})
test_that("aggregate_values_by_key stops if namespace is incorrect", {
  expect_error({
    agg <- aggregate_values_by_key(x = isadf,
                                   association_file = association_file,
                                   key = "SubjectID",
                                   lambda = "sum",
                                   args = list(na.rm = TRUE),
                                   namespace = 1)
  })
  expect_error({
    agg <- aggregate_values_by_key(x = isadf,
                                   association_file = association_file,
                                   key = "SubjectID",
                                   lambda = "sum",
                                   args = list(na.rm = TRUE),
                                   namespace = c("base", "ISAnalytics"))
  })
  expect_error({
    foo <- function(x) {
      sum(x)
    }
    env <- rlang::current_env()
    agg <- aggregate_values_by_key(x = isadf,
                                   association_file = association_file,
                                   key = "SubjectID",
                                   lambda = "foo",
                                   args = NULL,
                                   namespace = NULL,
                                   env = env)
  }, regexp = NA)
})
test_that("aggregate_values_by_key stops if env is incorrect", {
  expect_error({
    agg <- aggregate_values_by_key(x = isadf,
                                   association_file = association_file,
                                   key = "SubjectID",
                                   lambda = "sum",
                                   args = list(na.rm = TRUE),
                                   namespace = "base",
                                   env = NULL)
  })
  expect_error({
    agg <- aggregate_values_by_key(x = isadf,
                                   association_file = association_file,
                                   key = "SubjectID",
                                   lambda = "sum",
                                   args = list(na.rm = TRUE),
                                   namespace = "base",
                                   env = "global")
  })
})
test_that("aggregate_values_by_key stops if problems with function chosen", {
  expect_error({
    agg <- aggregate_values_by_key(x = isadf,
                                   association_file = association_file,
                                   key = "SubjectID",
                                   lambda = "sum",
                                   args = list(na.rm = TRUE),
                                   namespace = "ertfweg"
                                  )
  })
  expect_error({
    agg <- aggregate_values_by_key(x = isadf,
                                   association_file = association_file,
                                   key = "SubjectID",
                                   lambda = "seefwegfe",
                                   args = list(na.rm = TRUE),
                                   namespace = "base"
    )
  })
  env <- rlang::current_env()
  n <- 4
  expect_error({
    agg <- aggregate_values_by_key(x = isadf,
                                   association_file = association_file,
                                   key = "SubjectID",
                                   lambda = "foo",
                                   args = list(na.rm = TRUE),
                                   namespace = NULL,
                                   env = env)
  })
  expect_error({
    agg <- aggregate_values_by_key(x = isadf,
                                   association_file = association_file,
                                   key = "SubjectID",
                                   lambda = "n",
                                   args = list(na.rm = TRUE),
                                   namespace = NULL,
                                   env = env)
  })
})

# Test Values
test_that("aggregate_values_by_key succeeds for single IS matrix", {
  expect_error({
    agg <- aggregate_values_by_key(matrices$seqCount, association_file)
  }, regexp = NA)
  expect_true(all(c("chr", "integration_locus", "strand", "GeneName",
                    "GeneStrand", "SubjectID",
                    "Aggregated_value") %in% colnames(agg)))
})

test_that("aggregate_values_by_key succeeds for list of IS matrices", {
  expect_error({
    agg <- aggregate_values_by_key(matrices, association_file)
  }, regexp = NA)
  expect_true(all(c("chr", "integration_locus", "strand", "GeneName",
                    "GeneStrand", "SubjectID",
                    "Aggregated_value") %in% colnames(agg[[1]])))
  expect_true(all(c("chr", "integration_locus", "strand", "GeneName",
                    "GeneStrand", "SubjectID",
                    "Aggregated_value") %in% colnames(agg[[2]])))
})

