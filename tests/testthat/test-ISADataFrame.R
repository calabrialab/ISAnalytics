context("Building ISADataFrames")
library(ISAnalytics)
library(tibble)

listToTibble_casual <- list(a = 1:10, b = 10:1)
listToTibble_correct <- list(chr = c(as.character(1:10)), integration_locus = runif(10, min=100, max =10000),
                             strand = sample(c("+", "-"), 10, replace = TRUE),
                             exp_1 = runif(10, min = 0, max = 10000),
                             exp_2 = runif(10, min = 0, max = 10000),
                             exp_3 = runif(10, min = 0, max = 10000))
test_df_casual <- as.data.frame(listToTibble_casual)
test_tibble_casual <- as_tibble(listToTibble_casual)
test_df_correct <- as.data.frame(listToTibble_correct)
test_tibble_correct <- as_tibble(listToTibble_correct)

#--------------------------------------------------------------------------------------------------------------------------#
# Tests new_ISADataFrame
#--------------------------------------------------------------------------------------------------------------------------#

test_that("new_ISADataFrame works when input is a generic named list", {
  df <- new_ISADataFrame(listToTibble_casual)
  df2 <- new_ISADataFrame(listToTibble_correct)
  expect_s3_class(df, "ISADataFrame")
  expect_equal(attr(df, "mandatoryVars"), c("chr", "integration_locus", "strand"))
  expect_equal(attr(df, "metadata"), character())
  expect_s3_class(df2, "ISADataFrame")
  expect_equal(attr(df2, "mandatoryVars"), c("chr", "integration_locus", "strand"))
  expect_equal(attr(df2, "metadata"), character())
})

test_that("new_ISADataFrame fails when input is not a named list or a data frame", {
  expect_error(df <- new_ISADataFrame(1:10)) #with numeric vector
  expect_error(df <- new_ISADataFrame(list(1:10, 10:1))) #with unnamed list
})

test_that("new_ISADataFrame works when the input is a dataframe", {
  isadf <- new_ISADataFrame(test_df_casual) # data.frame with arbitrary values
  expect_s3_class(isadf, "ISADataFrame")
  expect_equal(attr(isadf, "mandatoryVars"), c("chr", "integration_locus", "strand"))
  isadf <- new_ISADataFrame(test_tibble_casual) # tibble with arbitrary values
  expect_s3_class(isadf, "ISADataFrame")
  expect_equal(attr(isadf, "mandatoryVars"), c("chr", "integration_locus", "strand"))
  isadf <- new_ISADataFrame(test_df_correct) # data.frame with correct values
  expect_s3_class(isadf, "ISADataFrame")
  expect_equal(attr(isadf, "mandatoryVars"), c("chr", "integration_locus", "strand"))
  isadf <- new_ISADataFrame(test_tibble_correct) # tibble with correct values
  expect_s3_class(isadf, "ISADataFrame")
  expect_equal(attr(isadf, "mandatoryVars"), c("chr", "integration_locus", "strand"))
})

#--------------------------------------------------------------------------------------------------------------------------#
# Tests validate_ISADataFrame
#--------------------------------------------------------------------------------------------------------------------------#
test_that("validate_ISADataFrame passes with correctly built ISAdf", {
  isadf <- new_ISADataFrame(listToTibble_correct) #With list
  expect_true(validate_ISADataFrame(isadf))
  isadf <- new_ISADataFrame(test_tibble_correct) #With tibble
  expect_true(validate_ISADataFrame(isadf))
  isadf <- new_ISADataFrame(test_df_correct) #With data.frame
  expect_true(validate_ISADataFrame(isadf))
})

test_that("validate_ISADataFrame fails if input is not an ISAdf", {
  expect_error(validate_ISADataFrame(test_tibble_correct))
  expect_error(validate_ISADataFrame(test_df_correct))
})

test_that("validate_ISADataFrame fails if input does not contain mandatory variables", {
  isadf <- new_ISADataFrame(listToTibble_casual) # with casual list
  expect_error(validate_ISADataFrame(isadf),
               regexp = "Validation of ISADataFrame failed: the input data doesn't contain the mandatory variables")
  isadf <- new_ISADataFrame(test_tibble_casual) # with casual tibble
  expect_error(validate_ISADataFrame(isadf),
               regexp = "Validation of ISADataFrame failed: the input data doesn't contain the mandatory variables")
  isadf <- new_ISADataFrame(test_df_casual) # with casual data.frame
  expect_error(validate_ISADataFrame(isadf),
               regexp = "Validation of ISADataFrame failed: the input data doesn't contain the mandatory variables")
})

test_that("validate_ISADataFrame throws a warning if declared metadata variables aren't present in the dataset", {
  isadf <- new_ISADataFrame(listToTibble_correct) # with correct list
  attr(isadf, "metadata") <- c("meta1, meta2")
  expect_warning(validate_ISADataFrame(isadf),
                 regexp = "Validation of ISADataFrame - warning: the input data doesn't contain the specified metadata columns")
  isadf <- new_ISADataFrame(test_tibble_correct) # with correct tibble
  attr(isadf, "metadata") <- c("meta1, meta2")
  expect_warning(validate_ISADataFrame(isadf),
                 regexp = "Validation of ISADataFrame - warning: the input data doesn't contain the specified metadata columns")
  isadf <- new_ISADataFrame(test_df_correct) # with correct data.frame
  attr(isadf, "metadata") <- c("meta1, meta2")
  expect_warning(validate_ISADataFrame(isadf),
                 regexp = "Validation of ISADataFrame - warning: the input data doesn't contain the specified metadata columns")
})

test_that("Validate_ISADataFrame fails if no experimental columns are detected", {
  #---Creating structures with NO exp data---#
  noExpList <- list(chr = c(as.character(1:10)), integration_locus = runif(10, min=100, max =10000),
                    strand = sample(c("+", "-"), 10, replace = TRUE))
  noExpListWithMeta <- list(chr = c(as.character(1:10)), integration_locus = runif(10, min=100, max =10000),
                            strand = sample(c("+", "-"), 10, replace = TRUE), meta1 = rep_len("met1", 10))
  noExpTibble <- as_tibble(noExpList)
  noExpTibbleWithMeta <- as_tibble(noExpListWithMeta)
  noExpDf <- as.data.frame(noExpList)
  noExpDfWithMeta <- as.data.frame(noExpListWithMeta)

  isadf <- new_ISADataFrame(noExpList) # with list
  isadf2 <- new_ISADataFrame(noExpListWithMeta, meta = "meta1") # with list and meta
  expect_error(validate_ISADataFrame(isadf), regexp = "Validation of ISADataFrame failed: no experimental variables found")
  expect_error(validate_ISADataFrame(isadf2), regexp = "Validation of ISADataFrame failed: no experimental variables found")

  isadf <- new_ISADataFrame(noExpTibble) # with tibble
  isadf2 <- new_ISADataFrame(noExpTibbleWithMeta, meta = "meta1") # with tibble and meta
  expect_error(validate_ISADataFrame(isadf), regexp = "Validation of ISADataFrame failed: no experimental variables found")
  expect_error(validate_ISADataFrame(isadf2), regexp = "Validation of ISADataFrame failed: no experimental variables found")

  isadf <- new_ISADataFrame(noExpDf) # with df
  isadf2 <- new_ISADataFrame(noExpDfWithMeta, meta = "meta1") # with df and meta
  expect_error(validate_ISADataFrame(isadf), regexp = "Validation of ISADataFrame failed: no experimental variables found")
  expect_error(validate_ISADataFrame(isadf2), regexp = "Validation of ISADataFrame failed: no experimental variables found")
})

test_that("Validate_ISADataFrame throws a warning if exp data is detected but it's not numeric", {
  #---Creating structures with non numeric exp data---#
  nonNumExpList <- list(chr = c(as.character(1:10)), integration_locus = runif(10, min=100, max =10000),
                        strand = sample(c("+", "-"), 10, replace = TRUE), exp_1 = rep_len("exp1", 10))
  nonNumExpListWithMeta <- list(chr = c(as.character(1:10)), integration_locus = runif(10, min=100, max =10000),
                      strand = sample(c("+", "-"), 10, replace = TRUE), exp_1 = rep_len("exp1", 10), meta_1 = rep_len("met1", 10))
  nonNumExpTibble <- as_tibble(nonNumExpList)
  nonNumExpTibbleWithMeta <- as_tibble(nonNumExpListWithMeta)
  nonNumExpDf <- as.data.frame(nonNumExpList)
  nonNumExpDfWithMeta <- as.data.frame(nonNumExpListWithMeta)

  isadf <- new_ISADataFrame(nonNumExpList) # with list
  isadf2 <- new_ISADataFrame(nonNumExpListWithMeta) # with list and meta
  expect_warning(validate_ISADataFrame(isadf),
                 regexp = "Validation of ISADataFrame - warning: found experimental columns with non numeric type")
  expect_warning(validate_ISADataFrame(isadf2),
                 regexp = "Validation of ISADataFrame - warning: found experimental columns with non numeric type")

  isadf <- new_ISADataFrame(nonNumExpTibble) # with tibble
  isadf2 <- new_ISADataFrame(nonNumExpTibbleWithMeta) # with tibble and meta
  expect_warning(validate_ISADataFrame(isadf),
                 regexp = "Validation of ISADataFrame - warning: found experimental columns with non numeric type")
  expect_warning(validate_ISADataFrame(isadf2),
                 regexp = "Validation of ISADataFrame - warning: found experimental columns with non numeric type")

  isadf <- new_ISADataFrame(nonNumExpDf) # with df
  isadf2 <- new_ISADataFrame(nonNumExpDfWithMeta) # with df and meta
  expect_warning(validate_ISADataFrame(isadf),
                 regexp = "Validation of ISADataFrame - warning: found experimental columns with non numeric type")
  expect_warning(validate_ISADataFrame(isadf2),
                 regexp = "Validation of ISADataFrame - warning: found experimental columns with non numeric type")
})
