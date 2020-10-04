context("Analysis functions")

library(ISAnalytics)

#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
op <- options(ISAnalytics.verbose = FALSE, ISAnalytics.widgets = FALSE)
on.exit(options(op))

path_AF <- system.file("extdata", "ex_association_file.tsv",
    package = "ISAnalytics"
)
root_correct <- system.file("extdata", "fs.zip",
    package = "ISAnalytics"
)
root_correct <- unzip_file_system(root_correct, "fs")
root_error <- system.file("extdata", "fserr.zip",
    package = "ISAnalytics"
)

root_error <- unzip_file_system(root_error, "fserr")

matrices_correct <- import_parallel_Vispa2Matrices_auto(
    association_file = path_AF, root = root_correct,
    quantification_type = c("seqCount", "fragmentEstimate"),
    matrix_type = "annotated", workers = 2, patterns = NULL,
    matching_opt = "ANY"
)

matrices_incorr <- suppressWarnings({
    import_parallel_Vispa2Matrices_auto(
        association_file = path_AF, root = root_error,
        quantification_type = c("fragmentEstimate", "seqCount"),
        matrix_type = "annotated", workers = 2, patterns = NULL,
        matching_opt = "ANY"
    )
})

#------------------------------------------------------------------------------#
# Test compute_abundance
#------------------------------------------------------------------------------#
smpl <- tibble::tibble(
    chr = c("1", "2", "3", "4", "5", "6"),
    integration_locus = c(
        14536,
        14544,
        14512,
        14236,
        14522,
        14566
    ),
    strand = c("+", "+", "-", "+", "-", "+"),
    CompleteAmplificationID = c(
        "ID1", "ID2", "ID1",
        "ID1", "ID3", "ID2"
    ),
    Value = c(3, 10, 40, 2, 15, 150),
    Value2 = c(456, 87, 87, 9, 64, 96),
    Value3 = c("a", "b", "c", "d", "e", "f")
)

test_that("compute_abundance stops if x is not tibble", {
    expect_error({
        ab <- compute_abundance(list(
            a = matrices_correct$seqCount,
            b = matrices_correct$fragmentEstimate
        ))
    })
})

test_that("compute_abundance stops if x is not integration matrix", {
    expect_error(
        {
            ab <- compute_abundance(tibble::tibble(a = c(1, 2, 3)))
        },
        regexp = .non_ISM_error()
    )
    expect_error(
        {
            ab <- compute_abundance(tibble::tibble(
                chr = c("1", "2", "3", "4", "5", "6"),
                integration_locus = c(
                    14536,
                    14544,
                    14512,
                    14236,
                    14522,
                    14566
                ),
                strand = c(
                    "+", "+", "-", "+",
                    "-", "+"
                )
            ))
        },
        regexp = .missing_complAmpID_error()
    )
    expect_error(
        {
            ab <- compute_abundance(smpl, columns = c("a"))
        },
        regexp = .missing_user_cols_error()
    )
    expect_error(
        {
            ab <- compute_abundance(smpl, columns = "Value3")
        },
        regexp = .non_num_user_cols_error()
    )
})

test_that("compute_abundance produces expected output", {
    ab <- compute_abundance(smpl)
    expected <- tibble::tibble(
        chr = c("1", "2", "3", "4", "5", "6"),
        integration_locus = c(
            14536,
            14544,
            14512,
            14236,
            14522,
            14566
        ),
        strand = c("+", "+", "-", "+", "-", "+"),
        CompleteAmplificationID = c(
            "ID1", "ID2", "ID1",
            "ID1", "ID3", "ID2"
        ),
        Value = c(3, 10, 40, 2, 15, 150),
        Value2 = c(456, 87, 87, 9, 64, 96),
        Value3 = c("a", "b", "c", "d", "e", "f"),
        Value_RelAbundance = c(
            0.06667, 0.0625, 0.88889,
            0.04445, 1.0, 0.9375
        ),
        Value_PercAbundance = Value_RelAbundance * 100
    )
    expect_equal(ab, expected, tolerance = 0.05)

    ab <- compute_abundance(smpl, columns = c("Value", "Value2"))
    expected <- tibble::tibble(
        chr = c("1", "2", "3", "4", "5", "6"),
        integration_locus = c(
            14536,
            14544,
            14512,
            14236,
            14522,
            14566
        ),
        strand = c("+", "+", "-", "+", "-", "+"),
        CompleteAmplificationID = c(
            "ID1", "ID2", "ID1",
            "ID1", "ID3", "ID2"
        ),
        Value = c(3, 10, 40, 2, 15, 150),
        Value2 = c(456, 87, 87, 9, 64, 96),
        Value3 = c("a", "b", "c", "d", "e", "f"),
        Value_RelAbundance = c(
            0.06667, 0.0625, 0.88889,
            0.04445, 1.0, 0.9375
        ),
        Value2_RelAbundance = c(
            0.82608, 0.475409, 0.157608,
            0.016304, 1.0, 0.52459
        ),
        Value_PercAbundance = Value_RelAbundance * 100,
        Value2_PercAbundance = Value2_RelAbundance * 100,
    )
    expect_equal(ab, expected, tolerance = 0.05)
})

#------------------------------------------------------------------------------#
# Test comparison_matrix
#------------------------------------------------------------------------------#
test_that("comparison_matrix throws error if input is incorrect", {
    # Input must be a named list, not single data frame
    expect_error({
        comp <- comparison_matrix(matrices_correct$seqCount)
    })
    # List names must be quantification types
    og_names <- names(matrices_correct)
    names(matrices_correct) <- c("a", "b")
    expect_error({
        comp <- comparison_matrix(matrices_correct)
    })
    # Every element of the list must be an integration matrix
    names(matrices_correct) <- og_names
    expect_error({
        comp <- comparison_matrix(list(seqCount = "a"))
    })
    expect_error(
        {
            comp <- comparison_matrix(
                list(
                    seqCount = matrices_correct$seqCount,
                    fragmentEstimate = matrices_correct$fragmentEstimate %>%
                        dplyr::select(-.data$chr)
                )
            )
        },
        regexp = .non_ISM_error()
    )
    expect_error(
        {
            comp <- comparison_matrix(
                list(
                    seqCount = matrices_correct$seqCount,
                    fragmentEstimate = matrices_correct$fragmentEstimate %>%
                        dplyr::select(-.data$CompleteAmplificationID)
                )
            )
        },
        regexp = .non_ISM_error()
    )
    expect_error(
        {
            comp <- comparison_matrix(
                list(
                    seqCount = matrices_correct$seqCount,
                    fragmentEstimate = matrices_correct$fragmentEstimate %>%
                        dplyr::select(-.data$Value)
                )
            )
        },
        regexp = .non_ISM_error()
    )
})

test_that("comparison_matrix notifies NA introduction", {
    op <- options(ISAnalytics.verbose = TRUE)
    on.exit(options(op))
    expect_message(
        {
            comp <- comparison_matrix(matrices_incorr)
        },
        regexp = .nas_introduced_msg()
    )
})

test_that("comparison_matrix produces correct output", {
    expect_message(
        {
            comp <- comparison_matrix(matrices_correct)
        },
        regexp = NA
    )
    expect_true(all(c("fragmentEstimate", "seqCount") %in% colnames(comp)))
    expect_true(nrow(comp) == nrow(matrices_correct$seqCount))
})

test_that("comparison_matrix supports custom names", {
    expect_message(
        {
            comp <- comparison_matrix(matrices_correct,
                seqCount = "sc",
                fragmentEstimate = "fe"
            )
        },
        regexp = NA
    )
    expect_true(all(c("fe", "sc") %in% colnames(comp)))
    expect_true(nrow(comp) == nrow(matrices_correct$seqCount))
    expect_true(is.integer(comp$sc))
    expect_true(is.numeric(comp$fe))
})

#------------------------------------------------------------------------------#
# Test separate_quant_matrices
#------------------------------------------------------------------------------#
smpl <- tibble::tibble(
    chr = c("1", "2", "3"),
    integration_locus = c(1263, 1264, 1265),
    strand = c("+", "+", "+"),
    CompleteAmplificationID = c("ID1", "ID2", "ID3"),
    fragmentEstimate = c(432.43, 532.43, 23.43),
    seqCount = c(34, 435, 65),
    random_col = c(6483, 486, 873)
)

test_that("separate_quant_matrices stops if param incorrect", {
    # Must be single data frame
    expect_error({
        sep <- separate_quant_matrices(list(a = smpl, b = smpl))
    })
    # Must be an integration matrix
    expect_error(
        {
            sep <- separate_quant_matrices(tibble::tibble(a = c(1, 2, 3)))
        },
        regexp = .non_ISM_error()
    )
    # Must contain CompleteAmplificationID column
    expect_error(
        {
            sep <- separate_quant_matrices(tibble::tibble(
                chr = c("1", "2", "3"),
                integration_locus = c(1263, 1264, 1265),
                strand = c("+", "+", "+")
            ))
        },
        regexp = .missing_complAmpID_error()
    )
    # Must contain numeric columns
    expect_error(
        {
            sep <- separate_quant_matrices(tibble::tibble(
                chr = c("1", "2", "3"),
                integration_locus = c(1263, 1264, 1265),
                strand = c("+", "+", "+"),
                CompleteAmplificationID = c("ID1", "ID2", "ID3")
            ))
        },
        regexp = .missing_num_cols_error()
    )
    # Must contain at least one quantification
    expect_error(
        {
            sep <- separate_quant_matrices(tibble::tibble(
                chr = c("1", "2", "3"),
                integration_locus = c(1263, 1264, 1265),
                strand = c("+", "+", "+"),
                CompleteAmplificationID = c("ID1", "ID2", "ID3"),
                random_col = c(6483, 486, 873)
            ))
        },
        regexp = .non_quant_cols_error()
    )
})

test_that("separate_quant_matrices warns if additional columns", {
    op <- options(ISAnalytics.verbose = TRUE)
    on.exit(options(op))
    expect_message(
        {
            sep <- separate_quant_matrices(smpl)
        },
        regexp = .non_quant_cols_msg("random_col")
    )
    expected_output <- list(
        fragmentEstimate =
            tibble::tibble(
                chr = c("1", "2", "3"),
                integration_locus = c(
                    1263,
                    1264,
                    1265
                ),
                strand = c("+", "+", "+"),
                CompleteAmplificationID = c(
                    "ID1",
                    "ID2",
                    "ID3"
                ),
                random_col = c(6483, 486, 873),
                Value = c(
                    432.43, 532.43,
                    23.43
                )
            ),
        seqCount =
            tibble::tibble(
                chr = c("1", "2", "3"),
                integration_locus = c(
                    1263,
                    1264,
                    1265
                ),
                strand = c("+", "+", "+"),
                CompleteAmplificationID = c(
                    "ID1",
                    "ID2",
                    "ID3"
                ),
                random_col = c(6483, 486, 873),
                Value = c(34, 435, 65)
            )
    )
    expect_equal(sep, expected_output)
})

test_that("separate_quant_matrices supports custom names", {
    colnames(smpl)[c(5, 6)] <- c("fe", "sc")
    sep <- separate_quant_matrices(smpl,
        fragmentEstimate = "fe",
        seqCount = "sc"
    )
    expected_output <- list(
        fragmentEstimate =
            tibble::tibble(
                chr = c("1", "2", "3"),
                integration_locus = c(
                    1263,
                    1264,
                    1265
                ),
                strand = c("+", "+", "+"),
                CompleteAmplificationID = c(
                    "ID1",
                    "ID2",
                    "ID3"
                ),
                random_col = c(6483, 486, 873),
                Value = c(
                    432.43, 532.43,
                    23.43
                )
            ),
        seqCount =
            tibble::tibble(
                chr = c("1", "2", "3"),
                integration_locus = c(
                    1263,
                    1264,
                    1265
                ),
                strand = c("+", "+", "+"),
                CompleteAmplificationID = c(
                    "ID1",
                    "ID2",
                    "ID3"
                ),
                random_col = c(6483, 486, 873),
                Value = c(34, 435, 65)
            )
    )
    expect_equal(sep, expected_output)
})

#------------------------------------------------------------------------------#
# Test threshold_filter
#------------------------------------------------------------------------------#
example_df <- tibble::tibble(
    a = c(20, 30, 40), b = c(40, 50, 60),
    c = c("a", "b", "c"), d = c(3L, 4L, 5L)
)

example_list <- list(
    first = example_df,
    second = example_df,
    third = example_df
)

example_list2 <- example_list
names(example_list2) <- NULL

test_that("threshold_filter gives errors with df - params wrong", {
    expect_error(
        {
            threshold_filter(example_df, threshold = list(1, 2))
        },
        regexp = .list_params_err(),
        fixed = TRUE
    )
    expect_error(
        {
            threshold_filter(example_df, threshold = "1")
        },
        regexp = .threshold_err(),
        fixed = TRUE
    )
    expect_error(
        {
            threshold_filter(example_df, threshold = 1, cols_to_compare = c(1, 2))
        },
        regexp = paste("`cols_to_compare` must be character"),
        fixed = TRUE
    )
    expect_error(
        {
            threshold_filter(example_df, threshold = 1, comparators = 1)
        },
        regexp = .comparators_err(),
        fixed = TRUE
    )
    expect_error(
        {
            threshold_filter(example_df, threshold = 1, comparators = "%")
        },
        regexp = .comparators_err(),
        fixed = TRUE
    )
    expect_error(
        {
            threshold_filter(example_df,
                threshold = 1,
                cols_to_compare = c("a", "f")
            )
        },
        regexp = .missing_user_cols_error(),
        fixed = TRUE
    )
    expect_error(
        {
            threshold_filter(example_df,
                threshold = 1,
                cols_to_compare = c("a", "c")
            )
        },
        regexp = .non_num_user_cols_error(),
        fixed = TRUE
    )
    expect_error(
        {
            threshold_filter(example_df,
                threshold = 1,
                cols_to_compare = c("a", "b", "d"),
                comparators = c(">", "<")
            )
        },
        regexp = .diff_leng_args(),
        fixed = TRUE
    )
})

test_that("threshold_filter gives result with df", {
    filtered <- threshold_filter(example_df,
        threshold = c(25, 55),
        cols_to_compare = c("a", "b"),
        comparators = ">"
    )
    expected <- tibble::tibble(a = c(40), b = c(60), c = c("c"), d = c(5))
    expect_equal(filtered, expected)
})

test_that("threshold_filter gives errors with list - params wrong", {
    ### List doesn't contain data frames
    expect_error(
        {
            threshold_filter(list(one = 1, two = 2), threshold = 1)
        },
        regexp = paste("Elements of the list must be data frames")
    )
    expect_error(
        {
            threshold_filter(list(one = example_df, two = 2), threshold = 1)
        },
        regexp = paste("Elements of the list must be data frames")
    )

    ### All param vectors but wrong types
    expect_error(
        {
            threshold_filter(example_list, threshold = c("1", "2"))
        },
        regexp = .threshold_err()
    )
    expect_error(
        {
            threshold_filter(example_list,
                threshold = c(1, 2),
                cols_to_compare = c(1, 2)
            )
        },
        regexp = paste("`cols_to_compare` must be character")
    )
    expect_error(
        {
            threshold_filter(example_list,
                threshold = c(1, 2),
                cols_to_compare = c("a", "b"), comparators = c(2, 4)
            )
        },
        regexp = .comparators_err()
    )
    expect_error(
        {
            threshold_filter(example_list,
                threshold = c(1, 2),
                cols_to_compare = c("a", "b"),
                comparators = c("%in%")
            )
        },
        regexp = .comparators_err()
    )

    ### Params as lists with unnamed list
    expect_error(
        {
            threshold_filter(example_list2,
                threshold = list(a = c(1, 2), b = c("a", "b"))
            )
        },
        regexp = .unnamed_list_err(),
        fixed = TRUE
    )
    expect_error(
        {
            threshold_filter(example_list2,
                threshold = c(1, 2),
                cols_to_compare = list(a = c("a", "b"))
            )
        },
        regexp = .unnamed_list_err(),
        fixed = TRUE
    )
    expect_error(
        {
            threshold_filter(example_list2,
                threshold = c(1, 2),
                cols_to_compare = c("a", "b"),
                comparators = list(a = c("<", ">"))
            )
        },
        regexp = .unnamed_list_err(),
        fixed = TRUE
    )

    ### List params without names
    expect_error(
        {
            threshold_filter(example_list, threshold = list(c(1, 2), c(2, 3)))
        },
        regexp = .names_list_param_err()
    )
    expect_error(
        {
            threshold_filter(example_list,
                threshold = list(
                    a = c(1, 2),
                    b = c(2, 3)
                ),
                cols_to_compare = list(c("a", "b"), c("b", "c"))
            )
        },
        regexp = .names_list_param_err()
    )
    expect_error(
        {
            threshold_filter(example_list,
                threshold = list(
                    a = c(1, 2),
                    b = c(2, 3)
                ),
                cols_to_compare = list(
                    a = c("a", "b"),
                    b = c("b", "c")
                ),
                comparators = list(c("<", ">"), c(">", "<"))
            )
        },
        regexp = .names_list_param_err()
    )

    ### List params with internal wrong types
    expect_error(
        {
            threshold_filter(example_list, threshold = list(
                a = c("1", "2"),
                b = c(2, 3)
            ))
        },
        regexp = .threshold_err()
    )
    expect_error(
        {
            threshold_filter(example_list,
                threshold = list(
                    a = c(1, 2),
                    b = c(2, 3)
                ),
                cols_to_compare = list(
                    a = c(1, 2),
                    b = c("b", "c")
                ),
                comparators = list(c("<", ">"), c(">", "<"))
            )
        },
        regexp = paste("`cols_to_compare` elements must be character")
    )
    expect_error(
        {
            threshold_filter(example_list,
                threshold = list(
                    a = c(1, 2),
                    b = c(2, 3)
                ),
                cols_to_compare = list(
                    a = c("a", "b"),
                    b = c("b", "c")
                ),
                comparators = list(
                    a = c("%", ">"),
                    b = c(">", "<")
                )
            )
        },
        regexp = .comparators_err()
    )

    ### Names in param not same names as list
    expect_error(
        {
            threshold_filter(example_list,
                threshold = list(
                    a = c(1, 2),
                    b = c(2, 3)
                ),
                cols_to_compare = list(
                    a = c("a", "b"),
                    b = c("b", "c")
                ),
                comparators = list(
                    a = c("==", ">"),
                    b = c(">", "<")
                )
            )
        },
        regexp = .param_names_not_in_list_err(),
        fixed = TRUE
    )

    ### Names not all equal
    expect_error(
        {
            threshold_filter(example_list,
                threshold = list(
                    first = c(1, 2),
                    second = c(2, 3)
                ),
                cols_to_compare = list(
                    second = c("a", "b"),
                    third = c("b", "c")
                ),
                comparators = list(first = c("==", ">"))
            )
        },
        regexp = .param_names_not_equal_err(),
        fixed = TRUE
    )
    expect_error(
        {
            threshold_filter(example_list,
                threshold = list(
                    first = c(1, 2),
                    second = c(2, 3)
                ),
                cols_to_compare = list(
                    first = c("a", "b"),
                    second = c("b", "c")
                ),
                comparators = list(first = c("==", ">"))
            )
        },
        regexp = .param_names_not_equal_err(),
        fixed = TRUE
    )

    ### Length params not equal
    expect_error(
        {
            threshold_filter(example_list,
                threshold = list(
                    first = c(1, 2, 3),
                    second = c(2, 3)
                ),
                cols_to_compare = list(
                    first = c("a", "b"),
                    second = c("b", "c")
                ),
                comparators = list(
                    first = c("==", ">"),
                    second = c("==", ">")
                )
            )
        },
        regexp = .diff_leng_args(),
        fixed = TRUE
    )
    expect_error(
        {
            threshold_filter(example_list,
                threshold = list(
                    first = c(1, 2),
                    second = c(2, 3)
                ),
                cols_to_compare = list(
                    first = c("a", "b"),
                    second = c("b")
                ),
                comparators = list(
                    first = c("==", ">"),
                    second = c("==")
                )
            )
        },
        regexp = .diff_leng_args(),
        fixed = TRUE
    )
    expect_error(
        {
            threshold_filter(example_list,
                threshold = c(1, 2),
                cols_to_compare = c("V1", "V2", "V3"),
                comparators = c("<", ">")
            )
        },
        regexp = .diff_leng_args(),
        fixed = TRUE
    )
    expect_error(
        {
            threshold_filter(example_list,
                threshold = c(1),
                cols_to_compare = c("V1", "V2", "V3"),
                comparators = c("<", ">")
            )
        },
        regexp = .diff_leng_args(),
        fixed = TRUE
    )

    ### Cols not in df
    expect_error(
        {
            threshold_filter(example_list,
                threshold = 1,
                cols_to_compare = list(
                    first = c("a", "b"),
                    second = c("aa", "bb")
                ),
                comparators = c("<", ">")
            )
        },
        regexp = .missing_user_cols_list_error(c("aa", "bb"), "second"),
        fixed = TRUE
    )
    expect_error(
        {
            threshold_filter(example_list,
                threshold = 1,
                cols_to_compare = c("aa", "bb"),
                comparators = c("<", ">")
            )
        },
        regexp = .missing_user_cols_error(),
        fixed = TRUE
    )
    expect_error(
        {
            threshold_filter(example_list,
                threshold = list(first = c(1, 2), second = c(2, 3)),
                cols_to_compare = c("aa", "bb"),
                comparators = c("<", ">")
            )
        },
        regexp = .missing_user_cols_error(),
        fixed = TRUE
    )

    ### Cols not numeric
    expect_error(
        {
            threshold_filter(example_list,
                threshold = 1,
                cols_to_compare = list(
                    first = c("a", "b"),
                    second = c("c", "d")
                ),
                comparators = c("<", ">")
            )
        },
        regexp = .non_num_user_cols_error(),
        fixed = TRUE
    )
    expect_error(
        {
            threshold_filter(example_list,
                threshold = 1,
                cols_to_compare = c("c", "b"),
                comparators = c("<", ">")
            )
        },
        regexp = .non_num_user_cols_error(),
        fixed = TRUE
    )
    expect_error(
        {
            threshold_filter(example_list,
                threshold = list(first = c(1, 2), second = c(2, 3)),
                cols_to_compare = c("c", "b"),
                comparators = c("<", ">")
            )
        },
        regexp = .non_num_user_cols_error(),
        fixed = TRUE
    )
})

test_that("threshold_filter gives result with list", {
    ### Unnamed list
    filtered <- threshold_filter(example_list2,
        threshold = c(25, 55),
        cols_to_compare = c("a", "b"),
        comparators = ">"
    )
    expected <- tibble::tibble(a = c(40), b = c(60), c = c("c"), d = c(5))
    expect_equal(filtered, list(
        expected,
        expected,
        expected
    ))
    ### All vector params, named list
    filtered <- threshold_filter(example_list,
        threshold = c(25, 55),
        cols_to_compare = c("a", "b"),
        comparators = ">"
    )
    expected <- tibble::tibble(a = c(40), b = c(60), c = c("c"), d = c(5))
    expect_equal(filtered, list(
        first = expected,
        second = expected,
        third = expected
    ))
    ### All list params, named list
    filtered <- threshold_filter(example_list,
        threshold = list(
            first = c(20, 40),
            second = c(25)
        ),
        cols_to_compare = list(
            first = c("a", "b"),
            second = c("a")
        ),
        comparators = list(
            first = c("==", "=="),
            second = c(">")
        )
    )
    expected1 <- tibble::tibble(a = c(20), b = c(40), c = c("a"), d = c(3))
    expected2 <- tibble::tibble(
        a = c(30, 40), b = c(50, 60),
        c = c("b", "c"), d = c(4, 5)
    )
    expect_equal(filtered, list(
        first = expected1,
        second = expected2,
        third = example_df
    ))
    ### Mixed params, named list
    filtered <- threshold_filter(example_list,
        threshold = list(
            first = c(20, 40),
            second = c(25, 60)
        ),
        cols_to_compare = c("a", "b"),
        comparators = list(
            first = c("==", "=="),
            second = c(">", "<")
        )
    )
    expected1 <- tibble::tibble(a = c(20), b = c(40), c = c("a"), d = c(3))
    expected2 <- tibble::tibble(
        a = c(30), b = c(50),
        c = c("b"), d = c(4)
    )
    expect_equal(filtered, list(
        first = expected1,
        second = expected2,
        third = example_df
    ))
})

#------------------------------------------------------------------------------#
# Test top_integrations
#------------------------------------------------------------------------------#
test_that("top_integrations stops if params incorrect", {
    ## Missing value cols
    expect_error(
        {
            top_integrations(smpl)
        },
        regexp = .missing_user_cols_error()
    )
    ## Missing cols to keep
    expect_error(
        {
            top_integrations(smpl, columns = "Value", keep = "a")
        },
        regexp = .missing_user_cols_error()
    )
    ## No mand vars
    expect_error(
        {
            top_integrations(smpl[2:7], columns = "Value")
        },
        regexp = .non_ISM_error()
    )
})

smpl <- tibble::tibble(
    chr = c("1", "2", "3", "4", "5", "6"),
    integration_locus = c(
        14536,
        14544,
        14512,
        14236,
        14522,
        14566
    ),
    strand = c("+", "+", "-", "+", "-", "+"),
    CompleteAmplificationID = c(
        "ID1", "ID2", "ID1",
        "ID1", "ID3", "ID2"
    ),
    Value = c(3, 10, 40, 2, 15, 150),
    Value2 = c(456, 87, 87, 9, 64, 96),
    Value3 = c("a", "b", "c", "d", "e", "f")
)

test_that("top_integrations works correctly", {
    top <- top_integrations(smpl,
        n = 3, columns = c("Value", "Value2"),
        keep = "nothing"
    )
    expected <- tibble::tibble(
        chr = c("6", "3", "5"),
        integration_locus = c(14566, 14512, 14522),
        strand = c("+", "-", "-"),
        Value = c(150, 40, 15),
        Value2 = c(96, 87, 64)
    )
    expect_equal(top, expected)
    top <- top_integrations(smpl,
        n = 3, columns = c("Value", "Value2"),
        keep = c("Value3", "CompleteAmplificationID")
    )
    expected <- tibble::tibble(
        chr = c("6", "3", "5"),
        integration_locus = c(14566, 14512, 14522),
        strand = c("+", "-", "-"),
        Value = c(150, 40, 15),
        Value2 = c(96, 87, 64),
        Value3 = c("f", "c", "e"),
        CompleteAmplificationID = c("ID2", "ID1", "ID3")
    )
    expect_equal(top, expected)
})

#------------------------------------------------------------------------------#
# Test sample_statistics
#------------------------------------------------------------------------------#
association_file <- import_association_file(path_AF, root_correct)
## Test input
test_that("sample_statistics stops if param incorrect", {
    # x must be df
    expect_error({
        stats <- sample_statistics(matrices_correct,
            metadata = association_file
        )
    })
    # metadata must be df
    expect_error({
        stats <- sample_statistics(matrices_correct$seqCount,
            metadata = 1
        )
    })
    # Sample key
    expect_error({
        stats <- sample_statistics(matrices_correct$seqCount,
            metadata = association_file,
            sample_key = 1
        )
    })
    expect_error(
        {
            stats <- sample_statistics(matrices_correct$seqCount,
                metadata = association_file,
                sample_key = c("sdgfs", "fwre")
            )
        },
        regexp = "Key columns not found in the data frame"
    )
    # Value cols
    expect_error(
        {
            stats <- sample_statistics(matrices_correct$seqCount,
                metadata = association_file,
                value_columns = c("x", "y")
            )
        },
        regexp = "Value columns not found in the data frame"
    )
    expect_error(
        {
            stats <- sample_statistics(matrices_correct$seqCount,
                metadata = association_file,
                value_columns = "chr"
            )
        },
        regexp = "Some or all of value columns are not numeric"
    )
    # Functions
    expect_error({
        stats <- sample_statistics(matrices_correct$seqCount,
            metadata = association_file,
            functions = c("sum")
        )
    })
    expect_error(
        {
            stats <- sample_statistics(matrices_correct$seqCount,
                metadata = association_file,
                functions = list(sum = "sum")
            )
        },
        regexp = paste(
            "The function parameter should contain a list",
            "of either functions or formulas.",
            "See ?sample_statistics for details"
        ),
        fixed = TRUE
    )
})

## Test values
test_that("sample_statistics returns correctly", {
    res <- sample_statistics(matrices_correct$seqCount,
        metadata = association_file,
        functions = list(sum = sum)
    )
    expect_true(is.list(res))
    expect_true(all(c("CompleteAmplificationID", "Value_sum") %in%
        colnames(res$x)))
    expect_true("Value_sum" %in% colnames(res$metadata))
})
