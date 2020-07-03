context("Building ISADataFrames")
library(ISAnalytics)
library(tibble)

listToTibble_casual <- list(a = 1:10, b = 10:1)
listToTibble_correct <- list(chr = c(as.character(1:10)),
                             integration_locus = runif(10,
                                                       min = 100,
                                                       max = 10000),
                             strand = sample(c("+","-"),
                                             10, replace = TRUE),
                             exp_1 = runif(10, min = 0, max = 10000),
                             exp_2 = runif(10, min = 0, max = 10000),
                             exp_3 = runif(10, min = 0, max = 10000))
test_df_casual <- as.data.frame(listToTibble_casual)
test_tibble_casual <- tibble::as_tibble(listToTibble_casual)
test_df_correct <- as.data.frame(listToTibble_correct)
test_tibble_correct <- tibble::as_tibble(listToTibble_correct)

#------------------------------------------------------------------------------#
# Tests new_ISADataFrame
#------------------------------------------------------------------------------#

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
    expect_error(df <- new_ISADataFrame(1:10)) # with numeric vector
    expect_error(df <- new_ISADataFrame(list(1:10, 10:1))) # with unnamed list
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

test_that("ISADataFrame is extensible", {
    subClass <- function(x, y) {
        new_ISADataFrame(x, y = y, class = "subClass")
    }
    testSub <- subClass(test_tibble_correct, c("one, two"))
    expect_s3_class(testSub, c("subClass", "ISADataFrame"))
    expect_s3_class(testSub[1], c("subClass", "ISADataFrame"))
})

#------------------------------------------------------------------------------#
# Tests validate_ISADataFrame
#------------------------------------------------------------------------------#
#---Creating structures with NO exp data---#
noExpList <- list(chr = c(as.character(1:10)), integration_locus = runif(10, min = 100, max = 10000), strand = sample(c("+", "-"),
    10,
    replace = TRUE
))
noExpListWithMeta <- list(chr = c(as.character(1:10)), integration_locus = runif(10, min = 100, max = 10000), strand = sample(c(
    "+",
    "-"
), 10, replace = TRUE), meta1 = rep_len("met1", 10))
noExpTibble <- tibble::as_tibble(noExpList)
noExpTibbleWithMeta <- tibble::as_tibble(noExpListWithMeta)
noExpDf <- as.data.frame(noExpList)
noExpDfWithMeta <- as.data.frame(noExpListWithMeta)

test_that("validate_ISADataFrame passes with correctly built ISAdf", {
    isadf <- new_ISADataFrame(listToTibble_correct) # With list
    expect_true(validate_ISADataFrame(isadf))
    isadf <- new_ISADataFrame(test_tibble_correct) # With tibble
    expect_true(validate_ISADataFrame(isadf))
    isadf <- new_ISADataFrame(test_df_correct) # With data.frame
    expect_true(validate_ISADataFrame(isadf))
})

test_that("validate_ISADataFrame fails if input is not an ISAdf", {
    expect_error(validate_ISADataFrame(test_tibble_correct))
    expect_error(validate_ISADataFrame(test_df_correct))
})

test_that("validate_ISADataFrame fails if input does not contain mandatory variables", {
    isadf <- new_ISADataFrame(listToTibble_casual) # with casual list
    expect_error(validate_ISADataFrame(isadf), regexp = "Validation of ISADataFrame failed: the input data doesn't contain the mandatory variables")
    isadf <- new_ISADataFrame(test_tibble_casual) # with casual tibble
    expect_error(validate_ISADataFrame(isadf), regexp = "Validation of ISADataFrame failed: the input data doesn't contain the mandatory variables")
    isadf <- new_ISADataFrame(test_df_casual) # with casual data.frame
    expect_error(validate_ISADataFrame(isadf), regexp = "Validation of ISADataFrame failed: the input data doesn't contain the mandatory variables")
})

test_that("validate_ISADataFrame throws a warning if declared metadata variables aren't present in the dataset", {
    isadf <- new_ISADataFrame(listToTibble_correct) # with correct list
    attr(isadf, "metadata") <- c("meta1, meta2")
    expect_warning(validate_ISADataFrame(isadf), regexp = "Validation of ISADataFrame - warning: the input data doesn't contain the specified metadata columns")
    isadf <- new_ISADataFrame(test_tibble_correct) # with correct tibble
    attr(isadf, "metadata") <- c("meta1, meta2")
    expect_warning(validate_ISADataFrame(isadf), regexp = "Validation of ISADataFrame - warning: the input data doesn't contain the specified metadata columns")
    isadf <- new_ISADataFrame(test_df_correct) # with correct data.frame
    attr(isadf, "metadata") <- c("meta1, meta2")
    expect_warning(validate_ISADataFrame(isadf), regexp = "Validation of ISADataFrame - warning: the input data doesn't contain the specified metadata columns")
})

test_that("Validate_ISADataFrame fails if no experimental columns are detected", {
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
    nonNumExpList <- list(chr = c(as.character(1:10)), integration_locus = runif(10, min = 100, max = 10000), strand = sample(c(
        "+",
        "-"
    ), 10, replace = TRUE), exp_1 = rep_len("exp1", 10))
    nonNumExpListWithMeta <- list(chr = c(as.character(1:10)), integration_locus = runif(10, min = 100, max = 10000), strand = sample(c(
        "+",
        "-"
    ), 10, replace = TRUE), exp_1 = rep_len("exp1", 10), meta_1 = rep_len("met1", 10))
    nonNumExpTibble <- tibble::as_tibble(nonNumExpList)
    nonNumExpTibbleWithMeta <- tibble::as_tibble(nonNumExpListWithMeta)
    nonNumExpDf <- as.data.frame(nonNumExpList)
    nonNumExpDfWithMeta <- as.data.frame(nonNumExpListWithMeta)

    isadf <- new_ISADataFrame(nonNumExpList) # with list
    isadf2 <- new_ISADataFrame(nonNumExpListWithMeta) # with list and meta
    expect_warning(validate_ISADataFrame(isadf), regexp = "Validation of ISADataFrame - warning: found experimental columns with non numeric type")
    expect_warning(validate_ISADataFrame(isadf2), regexp = "Validation of ISADataFrame - warning: found experimental columns with non numeric type")

    isadf <- new_ISADataFrame(nonNumExpTibble) # with tibble
    isadf2 <- new_ISADataFrame(nonNumExpTibbleWithMeta) # with tibble and meta
    expect_warning(validate_ISADataFrame(isadf), regexp = "Validation of ISADataFrame - warning: found experimental columns with non numeric type")
    expect_warning(validate_ISADataFrame(isadf2), regexp = "Validation of ISADataFrame - warning: found experimental columns with non numeric type")

    isadf <- new_ISADataFrame(nonNumExpDf) # with df
    isadf2 <- new_ISADataFrame(nonNumExpDfWithMeta) # with df and meta
    expect_warning(validate_ISADataFrame(isadf), regexp = "Validation of ISADataFrame - warning: found experimental columns with non numeric type")
    expect_warning(validate_ISADataFrame(isadf2), regexp = "Validation of ISADataFrame - warning: found experimental columns with non numeric type")
})

#--------------------------------------------------------------------------------------------------------------------------#
# Tests ISADataFrame()
#--------------------------------------------------------------------------------------------------------------------------#
listToTibble_correctMeta <- list(
    chr = c(as.character(1:10)), integration_locus = runif(10, min = 100, max = 10000), strand = sample(c(
        "+",
        "-"
    ), 10, replace = TRUE), meta1 = rep_len("m1", 10), exp_1 = runif(10, min = 0, max = 10000), exp_2 = runif(10, min = 0, max = 10000),
    exp_3 = runif(10, min = 0, max = 10000)
)
test_df_correctMeta <- as.data.frame(listToTibble_correctMeta)
test_tibble_correctMeta <- tibble::as_tibble(listToTibble_correctMeta)

listToTibble_diffLeng <- list(
    chr = c(as.character(1:10)), integration_locus = runif(10, min = 100, max = 10000), strand = sample(c(
        "+",
        "-"
    ), 10, replace = TRUE), meta1 = rep_len("m1", 10), exp_1 = runif(5, min = 0, max = 10000), exp_2 = runif(10, min = 0, max = 10000),
    exp_3 = runif(8, min = 0, max = 10000)
)

listToTibble_correctNonNum <- list(chr = c(as.character(1:10)), integration_locus = runif(10, min = 100, max = 10000), strand = sample(c(
    "+",
    "-"
), 10, replace = TRUE), nonNumericdata = rep_len("random", 10), exp_1 = runif(10, min = 0, max = 10000), exp_2 = runif(10,
    min = 0,
    max = 10000
), exp_3 = runif(10, min = 0, max = 10000))
test_df_correctNonNum <- as.data.frame(listToTibble_correctNonNum)
test_tibble_correctNonNum <- tibble::as_tibble(listToTibble_correctNonNum)

listToTibble_correctAll <- list(
    chr = c(as.character(1:10)), integration_locus = runif(10, min = 100, max = 10000), strand = sample(c(
        "+",
        "-"
    ), 10, replace = TRUE), meta1 = rep_len("m1", 10), nonNumericdata = rep_len("random", 10), exp_1 = runif(5, min = 0, max = 10000),
    exp_2 = runif(10, min = 0, max = 10000), exp_3 = runif(8, min = 0, max = 10000)
)

test_that("ISADataFrame with try.correct is able to fix list element lengths", {
    expect_message(isadf <- ISADataFrame(listToTibble_diffLeng, metadata = "meta1"), regexp = "Warning - introduced NAs to fix issues in provided list")
    isadf <- ISADataFrame(listToTibble_diffLeng, metadata = "meta1")
    expect_s3_class(isadf, "ISADataFrame")
    expect_true(any(is.na(isadf$exp_1)))
    expect_true(any(is.na(isadf$exp_3)))
    expect_false(any(is.na(isadf$exp_2)))
})

test_that("ISADataFrame without try.correct is not able to fix list element lengths and fails", {
    expect_error(isadf <- ISADataFrame(listToTibble_diffLeng, metadata = "meta1", try.correct = FALSE))
})

test_that("ISADataFrame() fails if there are no mandatory vars", {
    expect_error(isadf <- ISADataFrame(listToTibble_casual))
    expect_error(isadf <- ISADataFrame(test_tibble_casual))
    expect_error(isadf <- ISADataFrame(test_df_casual))
})

test_that("ISADataFrame() fails if there is no experimental data", {
    expect_error(isadf <- ISADataFrame(noExpList))
    expect_error(isadf <- ISADataFrame(noExpListWithMeta, metadata = "meta1"))
    expect_error(isadf <- ISADataFrame(noExpDf))
    expect_error(isadf <- ISADataFrame(noExpDfWithMeta, metadata = "meta1"))
    expect_error(isadf <- ISADataFrame(noExpTibble))
    expect_error(isadf <- ISADataFrame(noExpTibbleWithMeta, metadata = "meta1"))
})

test_that("ISADataFrame() with try.correct corrects missing metadata", {
    isadf <- ISADataFrame(listToTibble_correctMeta, metadata = c("meta1", "meta2"))
    expect_s3_class(isadf, "ISADataFrame")
    expect_equal(attr(isadf, "metadata"), c("meta1"))
    expect_message(isadf <- ISADataFrame(listToTibble_correctMeta, metadata = c("meta1", "meta2")), regexp = paste("Auto-corrected: the input data does not",
                                                                                                                   "contain the specified metadata columns"))
    isadf <- ISADataFrame(test_df_correctMeta, metadata = c("meta1", "meta2"))
    expect_s3_class(isadf, "ISADataFrame")
    expect_equal(attr(isadf, "metadata"), c("meta1"))
    expect_message(isadf <- ISADataFrame(test_df_correctMeta, metadata = c("meta1", "meta2")), regexp = paste("Auto-corrected: the input data does not",
                                                                                                              "contain the specified metadata columns"))
    isadf <- ISADataFrame(test_tibble_correctMeta, metadata = c("meta1", "meta2"))
    expect_s3_class(isadf, "ISADataFrame")
    expect_equal(attr(isadf, "metadata"), c("meta1"))
    expect_message(isadf <- ISADataFrame(test_tibble_correctMeta, metadata = c("meta1", "meta2")), regexp = paste("Auto-corrected: the input data does not",
                                                                                                                  "contain the specified metadata columns"))
})

test_that("ISADataFrame() without try.correct fails and doesn't correct missing metadata", {
    expect_error(isadf <- ISADataFrame(listToTibble_correctMeta, metadata = c("meta1", "meta2"), try.correct = FALSE))
    expect_error(isadf <- ISADataFrame(test_df_correctMeta, metadata = c("meta1", "meta2"), try.correct = FALSE))
    expect_error(isadf <- ISADataFrame(test_tibble_correctMeta, metadata = c("meta1", "meta2"), try.correct = FALSE))
})

test_that("ISADataFrame() with try.correct corrects non numeric experimental data", {
    isadf <- ISADataFrame(listToTibble_correctNonNum)
    expect_s3_class(isadf, "ISADataFrame")
    expect_equal(attr(isadf, "metadata"), c("nonNumericdata"))
    expect_message(isadf <- ISADataFrame(listToTibble_correctNonNum), regexp = paste("Auto-corrected: found experimental",
                                                                                     "columns with non numeric type"))
    isadf <- ISADataFrame(test_df_correctNonNum)
    expect_s3_class(isadf, "ISADataFrame")
    expect_equal(attr(isadf, "metadata"), c("nonNumericdata"))
    expect_message(isadf <- ISADataFrame(test_df_correctNonNum), regexp = paste("Auto-corrected: found experimental",
                                                                                "columns with non numeric type"))
    isadf <- ISADataFrame(test_tibble_correctNonNum)
    expect_s3_class(isadf, "ISADataFrame")
    expect_equal(attr(isadf, "metadata"), c("nonNumericdata"))
    expect_message(isadf <- ISADataFrame(test_tibble_correctNonNum), regexp = paste("Auto-corrected: found experimental",
                                                                                    "columns with non numeric type"))
})

test_that("ISADataFrame() without try.correct fails and doesn't correct non numeric experimental data", {
    expect_error(isadf <- ISADataFrame(listToTibble_correctNonNum, try.correct = FALSE))
    expect_error(isadf <- ISADataFrame(test_df_correctNonNum, try.correct = FALSE))
    expect_error(isadf <- ISADataFrame(test_tibble_correctNonNum, try.correct = FALSE))
})

test_that("ISADataFrame() with try.correct is successful for all warnings", {
    isadf <- ISADataFrame(listToTibble_correctAll, metadata = c("meta1", "meta2", "meta3"))
    expect_s3_class(isadf, "ISADataFrame")
    expect_true(any(is.na(isadf$exp_1)))
    expect_true(any(is.na(isadf$exp_3)))
    expect_false(any(is.na(isadf$exp_2)))
    expect_equal(attr(isadf, "metadata"), c("meta1", "nonNumericdata"))
    expect_message(isadf <- ISADataFrame(listToTibble_correctAll, metadata = c("meta1", "meta2", "meta3")), regexp = paste("Auto-corrected: found experimental",
                                                                                                                           "columns with non numeric type"))
    expect_message(isadf <- ISADataFrame(listToTibble_correctAll, metadata = c("meta1", "meta2", "meta3")), regexp = paste("Auto-corrected: the input data does not",
                                                                                                                           "contain the specified metadata columns"))
    expect_message(isadf <- ISADataFrame(listToTibble_correctAll, metadata = c("meta1", "meta2", "meta3")), regexp = paste("Warning - introduced NAs to fix issues in",
                                                                                                                           "provided list"))
})

listToTibble_noExpOneNonNum <- list(chr = c(as.character(1:10)), integration_locus = runif(10, min = 100, max = 10000), strand = sample(c(
    "+",
    "-"
), 10, replace = TRUE), nonNumericdata = rep_len("random", 10))

test_that("ISADataFrame() with try.correct is able to identify missing exp data after correction", {
    expect_error(isadf <- ISADataFrame(listToTibble_noExpOneNonNum))
})
