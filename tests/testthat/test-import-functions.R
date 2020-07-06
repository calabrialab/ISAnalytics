context("Importing IS from files")

library(ISAnalytics)
library(tibble)

#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
# Annotated new matrix
example_path <- system.file("extdata", "ex_annotated_ISMatrix.tsv.xz",
                            package = "ISAnalytics")

# Old style not annotated matrix
example_path2 <- system.file("extdata", "ex_old_style_ISMatrix.tsv.xz",
                            package = "ISAnalytics")

# New not annotated matrix
example_path3 <- system.file("extdata", "ex_notann_ISMatrix.tsv.xz",
                             package = "ISAnalytics")

# Malformed matrix
example_path4 <- system.file("extdata", "ex_malformed_ISMatrix.tsv.xz",
                             package = "ISAnalytics")
#------------------------------------------------------------------------------#
# Tests .messy_to_tidy
#------------------------------------------------------------------------------#
test_that(".messy_to_tidy fails if input is not an ISADataFrame", {
  expect_error(.messy_to_tidy(data.frame(a = seq_len(10), b = seq_len(10))))
  expect_error(.messy_to_tidy(tibble(list(a = seq_len(10), b = seq_len(10)))))
})

test_that(".messy_to_tidy produces correct tidy ISADataFrame", {
  raw_isadf <- read.csv(example_path, sep = "\t", check.names = FALSE)
  raw_isadf <- new_ISADataFrame(raw_isadf, meta = c("GeneName", "GeneStrand"))
  tidy_isadf <- .messy_to_tidy(raw_isadf)
  # Resulting data frame must preserve ISADataFrame class
  expect_true(is.ISADataFrame(tidy_isadf))
  # Resulting ISADataFrame must have the columns "ExperimentID" and "Value"
  expect_true(all(c("ExperimentID", "Value") %in% colnames(tidy_isadf)))
  # The metadata field updates and contains the new ExperimentID column
  expect_equal(metadata(tidy_isadf), c("GeneName", "GeneStrand",
                                        "ExperimentID"))
  # The data frame must not contain null values
  expect_true(all(!is.na(tidy_isadf$Value)))
})

#------------------------------------------------------------------------------#
# Tests .auto_detect_type
#------------------------------------------------------------------------------#
test_that(".auto_detect_type detects new style annotated matrices", {
  raw_isadf <- read.csv(example_path, sep = "\t", check.names = FALSE)
  res <- .auto_detect_type(raw_isadf)
  expect_equal(res, "NEW_ANNOTATED")
})

test_that(".auto_detect_type detects old style matrices", {
  raw_isadf <- read.csv(example_path2, sep = "\t", check.names = FALSE)
  res <- .auto_detect_type(raw_isadf)
  expect_equal(res, "OLD")
})

test_that(".auto_detect_type detects new not annotated matrices", {
  raw_isadf <- read.csv(example_path, sep = "\t", check.names = FALSE)
  raw_isadf <- dplyr::select(raw_isadf, -c("GeneName", "GeneStrand"))
  res <- .auto_detect_type(raw_isadf)
  expect_equal(res, "NEW_NOTANN")
})

test_that(".auto_detect_type detects malformed matrices", {
  raw_isadf <- read.csv(example_path, sep = "\t", check.names = FALSE)
  raw_isadf1 <- dplyr::select(raw_isadf,
                              -c("GeneName", "GeneStrand", "integration_locus"))
  res <- .auto_detect_type(raw_isadf1)
  expect_equal(res, "MALFORMED")
  raw_isadf2 <- read.csv(example_path2, sep = "\t", check.names = FALSE)
  raw_isadf3 <- add_column(raw_isadf2, chr = raw_isadf$chr)
  res <- .auto_detect_type(raw_isadf3)
  expect_equal(res, "MALFORMED")
})

#------------------------------------------------------------------------------#
# Tests import_single_Vispa2Matrix
#------------------------------------------------------------------------------#
test_that("import_single_Vispa2Matrix succeeds with standard annotated", {
  isadf <- import_single_Vispa2Matrix(example_path)
  expect_s3_class(isadf, "ISADataFrame")
  expect_equal(metadata(isadf), c("GeneName", "GeneStrand", "ExperimentID"))
  expect_true(all(c("chr", "integration_locus", "strand") %in% colnames(isadf)))
  expect_true(all(!is.na(isadf$Value)))
})

test_that("import_single_Vispa2Matrix succeeds with old", {
  isadf <- import_single_Vispa2Matrix(example_path2)
  expect_s3_class(isadf, "ISADataFrame")
  expect_equal(metadata(isadf), c("ExperimentID"))
  expect_true(all(c("chr", "integration_locus", "strand") %in% colnames(isadf)))
  expect_true(all(!is.na(isadf$Value)))
})

test_that("import_single_Vispa2Matrix succeeds with new not annotated", {
  isadf <- import_single_Vispa2Matrix(example_path3)
  expect_s3_class(isadf, "ISADataFrame")
  expect_equal(metadata(isadf), c("ExperimentID"))
  expect_true(all(c("chr", "integration_locus", "strand") %in% colnames(isadf)))
  expect_true(all(!is.na(isadf$Value)))
})

test_that("import_single_Vispa2Matrix fails with malformed matrix", {
  expect_error(isadf <- import_single_Vispa2Matrix(example_path4))
})

test_that("import_single_Vispa2Matrix fails if file does not exist", {
  expect_error(import_single_Vispa2Matrix(system.file("extdata",
                                                      "matrix.tsv",
                                                      package = "ISAnalytics")))
})
