library(ISAnalytics)

#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
op <- withr::local_options(
    ISAnalytics.verbose = FALSE,
    ISAnalytics.widgets = FALSE
)

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
    matching_opt = "ANY", dates_format = "dmy", multi_quant_matrix = FALSE
)

matrices_incorr <- suppressWarnings({
    import_parallel_Vispa2Matrices_auto(
        association_file = path_AF, root = root_error,
        quantification_type = c("fragmentEstimate", "seqCount"),
        matrix_type = "annotated", workers = 2, patterns = NULL,
        matching_opt = "ANY", dates_format = "dmy", multi_quant_matrix = FALSE
    )
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
})

test_that("comparison_matrix notifies NA introduction", {
    withr::local_options(list(ISAnalytics.verbose = TRUE))
    expect_message({
        comp <- comparison_matrix(matrices_incorr)
    })
})

test_that("comparison_matrix produces correct output", {
    comp <- comparison_matrix(matrices_correct)
    expect_true(all(c("fragmentEstimate", "seqCount") %in% colnames(comp)))
    expect_true(nrow(comp) == nrow(matrices_correct$seqCount))
})

test_that("comparison_matrix supports custom names", {
    comp <- comparison_matrix(matrices_correct,
        seqCount = "sc",
        fragmentEstimate = "fe"
    )
    expect_true(all(c("fe", "sc") %in% colnames(comp)))
    expect_true(nrow(comp) == nrow(matrices_correct$seqCount))
    expect_true(is.numeric(comp$sc))
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
    expect_error({
        sep <- separate_quant_matrices(tibble::tibble(a = c(1, 2, 3)))
    })
    # Must contain at least one quantification
    expect_error(
        {
            sep <- separate_quant_matrices(tibble::tibble(
                chr = c("1", "2", "3"),
                integration_locus = c(1263, 1264, 1265),
                strand = c("+", "+", "+"),
                CompleteAmplificationID = c("ID1", "ID2", "ID3"),
                random_col = c(6483, 486, 873)
            ), key = c(mandatory_IS_vars(), "CompleteAmplificationID"))
        },
        regexp = .non_quant_cols_error()
    )
})

test_that("separate_quant_matrices warns if additional columns", {
    withr::local_options(list(ISAnalytics.verbose = TRUE))
    expect_message({
        sep <- separate_quant_matrices(smpl,
            key = c(
                mandatory_IS_vars(),
                "CompleteAmplificationID"
            )
        )
    })
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
        seqCount = "sc",
        key = c(mandatory_IS_vars(), "CompleteAmplificationID")
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
    expect_error({
        threshold_filter(example_df,
            threshold = 1,
            cols_to_compare = c("a", "f")
        )
    })
    expect_error({
        threshold_filter(example_df,
            threshold = 1,
            cols_to_compare = c("a", "c")
        )
    })
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
    expect_error({
        threshold_filter(example_list,
            threshold = 1,
            cols_to_compare = list(
                first = c("a", "b"),
                second = c("aa", "bb")
            ),
            comparators = c("<", ">")
        )
    })
    expect_error({
        threshold_filter(example_list,
            threshold = 1,
            cols_to_compare = c("aa", "bb"),
            comparators = c("<", ">")
        )
    })
    expect_error({
        threshold_filter(example_list,
            threshold = list(first = c(1, 2), second = c(2, 3)),
            cols_to_compare = c("aa", "bb"),
            comparators = c("<", ">")
        )
    })

    ### Cols not numeric
    expect_error({
        threshold_filter(example_list,
            threshold = 1,
            cols_to_compare = list(
                first = c("a", "b"),
                second = c("c", "d")
            ),
            comparators = c("<", ">")
        )
    })
    expect_error({
        threshold_filter(example_list,
            threshold = 1,
            cols_to_compare = c("c", "b"),
            comparators = c("<", ">")
        )
    })
    expect_error({
        threshold_filter(example_list,
            threshold = list(first = c(1, 2), second = c(2, 3)),
            cols_to_compare = c("c", "b"),
            comparators = c("<", ">")
        )
    })
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
# Test sample_statistics
#------------------------------------------------------------------------------#
association_file <- import_association_file(path_AF, root_correct,
    dates_format = "dmy"
)
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

#------------------------------------------------------------------------------#
# Test CIS_grubbs
#------------------------------------------------------------------------------#
test_that("CIS_grubbs fails if x is not standard matrix", {
    wrong <- tibble::tibble(a = c(1, 2, 3), b = c(4, 5, 6))
    expect_error(
        {
            CIS_grubbs(wrong)
        },
        regexp = .non_ISM_error()
    )
    expect_error(
        {
            CIS_grubbs(smpl)
        },
        regexp = .missing_annot()
    )
})

test_that("CIS_grubbs fails if file doesn't exist", {
    expect_error({
        CIS_grubbs(matrices_correct$seqCount,
            genomic_annotation_file = "myfile.tsv"
        )
    })
})

test_that("CIS_grubbs produces correct df", {
    output_cols <- c(
        "GeneName", "GeneStrand", "chr", "n_IS_perGene",
        "min_bp_integration_locus", "max_bp_integration_locus",
        "IS_span_bp", "avg_bp_integration_locus",
        "median_bp_integration_locus", "distinct_orientations",
        "describe_vars", "describe_n", "describe_mean", "describe_sd",
        "describe_median", "describe_trimmed", "describe_mad",
        "describe_min", "describe_max", "describe_range",
        "describe_skew", "describe_kurtosis", "describe_se",
        "average_TxLen", "geneIS_frequency_byHitIS",
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
    result <- CIS_grubbs(matrices_correct$seqCount)
    expect_true(ncol(result) == 40)
    expect_true(all(output_cols %in% colnames(result)))
    result <- CIS_grubbs(matrices_correct$seqCount,
        add_standard_padjust = FALSE
    )
    expect_true(ncol(result) == 37)
    expect_true(all(output_cols[!output_cols %in% c(
        "tdist_bonferroni",
        "tdist_fdr",
        "tdist_benjamini"
    )] %in%
        colnames(result)))
})

#------------------------------------------------------------------------------#
# Test cumulative_count_union
#------------------------------------------------------------------------------#
minimal_agg_example <- tibble::tibble(
    chr = c(1, 2, 2, 5, 5, 3, 3, 10, 11),
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

test_that("cumulative_count_union stops if tp not in key", {
    expect_error(
        {
            cumulative_count_union(minimal_agg_example, key = c(
                "SubjectID",
                "CellMarker",
                "Tissue"
            ))
        },
        regexp = .key_without_tp_err()
    )
})

test_that("cumulative_count_union produces correct output", {
    cum_count <- cumulative_count_union(minimal_agg_example)
    expected <- tibble::tibble(
        SubjectID = c(
            rep_len("S1", 3),
            rep_len("S2", 2)
        ),
        CellMarker = c(
            rep_len("CM1", 3),
            rep_len("CM2", 2)
        ),
        Tissue = c(
            rep_len("T1", 3),
            rep_len("T2", 2)
        ),
        TimePoint = c(
            "0001", "0010", "0015",
            "0001", "0020"
        ),
        count = c(3, 4, 5, 2, 4)
    )
    expect_equal(cum_count, expected)
})

#------------------------------------------------------------------------------#
# Test is_sharing
#------------------------------------------------------------------------------#
test_sharing_input <- tibble::tribble(
    ~chr, ~integration_locus, ~strand, ~SubjectID,
    "1", 10032, "+", "S1",
    "10", 20162, "-", "S1",
    "5", 45612, "-", "S1",
    "2", 21206, "+", "S1",
    "1", 10102, "+", "S1",
    "3", 38167, "-", "S1",
    "1", 10032, "+", "S2",
    "5", 45612, "-", "S2",
    "6", 12542, "-", "S2",
    "1", 42571, "-", "S2",
    "4", 21054, "+", "S2",
    "4", 12072, "-", "S2",
    "11", 25722, "-", "S2",
    "10", 11725, "+", "S2",
    "10", 12247, "+", "S2",
    "11", 25722, "-", "S2",
    "1", 10032, "+", "S3",
    "2", 21206, "+", "S3",
    "6", 12542, "-", "S3",
    "1", 42571, "-", "S3",
    "2", 51232, "+", "S3"
)

test_expected_shared_counts <- tibble::tribble(
    ~group_id, ~num_IS,
    "S1", 6,
    "S2", 9,
    "S3", 5
)

test_expected_sharing_df <- tibble::tribble(
    ~group1, ~group2, ~shared, ~on_g1, ~on_g2, ~on_union,
    "S1", "S1", 6, 100, 100, 100,
    "S1", "S2", 2, (2 / 6) * 100, (2 / 9) * 100, (2 / 13) * 100,
    "S2", "S1", 2, (2 / 9) * 100, (2 / 6) * 100, (2 / 13) * 100,
    "S1", "S3", 2, (2 / 6) * 100, (2 / 5) * 100, (2 / 9) * 100,
    "S3", "S1", 2, (2 / 5) * 100, (2 / 6) * 100, (2 / 9) * 100,
    "S2", "S2", 9, 100, 100, 100,
    "S2", "S3", 3, (3 / 9) * 100, (3 / 5) * 100, (3 / 11) * 100,
    "S3", "S2", 3, (3 / 5) * 100, (3 / 9) * 100, (3 / 11) * 100,
    "S3", "S3", 5, 100, 100, 100
)

test_that("is_sharing expects mandatory cols", {
    expect_error({
        sharing <- is_sharing(test_sharing_input %>%
            dplyr::select(-.data$chr))
    })
})

test_that("is_sharing produces expected output", {
    sharing_res <- is_sharing(test_sharing_input,
        group_key = "SubjectID"
    )
    expect_equal(sharing_res$is_count, test_expected_shared_counts)
    expect_equal(sharing_res$sharing, test_expected_sharing_df)
})
