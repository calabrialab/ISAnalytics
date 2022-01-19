library(ISAnalytics)

#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
op <- withr::local_options(
    ISAnalytics.verbose = FALSE,
    ISAnalytics.reports = FALSE
)

path_AF <- system.file("testdata", "ex_association_file.tsv.gz",
    package = "ISAnalytics"
)
root_correct <- system.file("testdata", "fs.zip",
    package = "ISAnalytics"
)
root_correct <- unzip_file_system(root_correct, "fs")
root_error <- system.file("testdata", "fserr.zip",
    package = "ISAnalytics"
)

root_error <- unzip_file_system(root_error, "fserr")

matrices_correct <- import_parallel_Vispa2Matrices(
    association_file = path_AF, root = root_correct,
    quantification_type = c("seqCount", "fragmentEstimate"),
    matrix_type = "annotated", workers = 2, patterns = NULL,
    matching_opt = "ANY", dates_format = "dmy", multi_quant_matrix = FALSE,
    mode = "AUTO"
)

matrices_incorr <- suppressWarnings({
    import_parallel_Vispa2Matrices(
        association_file = path_AF, root = root_error,
        quantification_type = c("fragmentEstimate", "seqCount"),
        matrix_type = "annotated", workers = 2, patterns = NULL,
        matching_opt = "ANY", dates_format = "dmy", multi_quant_matrix = FALSE,
        mode = "AUTO"
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
    expect_error({
        stats <- sample_statistics(matrices_correct$seqCount,
            metadata = association_file,
            sample_key = c("sdgfs", "fwre")
        )
    })
    # Value cols
    expect_error({
        stats <- sample_statistics(matrices_correct$seqCount,
            metadata = association_file,
            value_columns = c("x", "y")
        )
    })
    expect_error({
        stats <- sample_statistics(matrices_correct$seqCount,
            metadata = association_file,
            value_columns = "chr"
        )
    })
    # Functions
    expect_error({
        stats <- sample_statistics(matrices_correct$seqCount,
            metadata = association_file,
            functions = c("sum")
        )
    })
    expect_error({
        stats <- sample_statistics(matrices_correct$seqCount,
            metadata = association_file,
            functions = list(sum = "sum")
        )
    })
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
    expect_error({
        CIS_grubbs(wrong)
    })
    expect_error({
        CIS_grubbs(smpl)
    })
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
        "vars", "n", "mean", "sd",
        "median", "trimmed", "mad",
        "min", "max", "range",
        "skew", "kurtosis", "se",
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
# Test cumulative_is
#------------------------------------------------------------------------------#

test_that("cumulative_is produces correct output", {
    c_is <- cumulative_is(minimal_agg_example)
    expect_true(nrow(c_is) == 5)
    expect_true(all(c("is", "cumulative_is") %in% colnames(c_is)))
    counts <- purrr::map_int(c_is$cumulative_is, ~ nrow(.x))
    expect_equal(counts, c(3, 4, 5, 2, 4))
    c_is <- cumulative_is(minimal_agg_example, keep_og_is = FALSE)
    expect_true(nrow(c_is) == 5)
    expect_true(all(c("cumulative_is") %in% colnames(c_is)))
    counts <- purrr::map_int(c_is$cumulative_is, ~ nrow(.x))
    expect_equal(counts, c(3, 4, 5, 2, 4))
    c_is <- cumulative_is(minimal_agg_example,
        keep_og_is = FALSE,
        expand = TRUE
    )
    expect_true(all(mandatory_IS_vars() %in% colnames(c_is)))
    expect_true(nrow(c_is) == sum(counts))
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
    PT01 = data.table::data.table(
        g1 = c(
            "PT01_CD34_BM_1", "PT01_CD34_BM_1",
            "PT01_CD34_BM_1", "PT01_CD34_BM_1",
            "PT01_CD34_BM_1", "PT01_CD34_BM_1",
            "PT01_CD34_BM_1", "PT01_CD34_BM_2",
            "PT01_CD34_BM_2", "PT01_CD34_BM_2",
            "PT01_CD34_BM_2", "PT01_CD34_BM_2",
            "PT01_CD34_BM_2", "PT01_CD34_BM_2",
            "PT01_Whole_BM_1", "PT01_Whole_BM_1",
            "PT01_Whole_BM_1", "PT01_Whole_BM_1",
            "PT01_Whole_BM_1", "PT01_Whole_BM_1",
            "PT01_Whole_BM_1", "PT01_Whole_BM_2",
            "PT01_Whole_BM_2", "PT01_Whole_BM_2",
            "PT01_Whole_BM_2", "PT01_Whole_BM_2",
            "PT01_Whole_BM_2", "PT01_Whole_BM_2"
        ),
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
            1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
            2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2,
            2, 2, 2, 2
        ),
        g2 = c(
            "PT01_CD15_PB_01", "PT01_CD14_PB_06",
            "PT01_CD14_PB_01", "PT01_CD14_PB_02",
            "PT01_CD13_PB_08", "PT01_CD14_PB_03",
            "PT01_CD13_PB_01", "PT01_CD15_PB_01",
            "PT01_CD14_PB_06", "PT01_CD14_PB_01",
            "PT01_CD14_PB_02", "PT01_CD13_PB_08",
            "PT01_CD14_PB_03", "PT01_CD13_PB_01",
            "PT01_CD15_PB_01", "PT01_CD14_PB_06",
            "PT01_CD14_PB_01", "PT01_CD14_PB_02",
            "PT01_CD13_PB_08", "PT01_CD14_PB_03",
            "PT01_CD13_PB_01", "PT01_CD15_PB_01",
            "PT01_CD14_PB_06", "PT01_CD14_PB_01",
            "PT01_CD14_PB_02", "PT01_CD13_PB_08",
            "PT01_CD14_PB_03", "PT01_CD13_PB_01"
        ),
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
            "CD15", "CD14", "CD14", "CD14",
            "CD13", "CD14", "CD13", "CD15",
            "CD14", "CD14", "CD14", "CD13",
            "CD14", "CD13", "CD15", "CD14",
            "CD14", "CD14", "CD13", "CD14",
            "CD13", "CD15", "CD14", "CD14",
            "CD14", "CD13", "CD14", "CD13"
        ),
        g2_Tissue = c(
            "PB", "PB", "PB", "PB", "PB", "PB",
            "PB", "PB", "PB", "PB", "PB", "PB",
            "PB", "PB", "PB", "PB", "PB", "PB",
            "PB", "PB", "PB", "PB", "PB", "PB",
            "PB", "PB", "PB", "PB"
        ),
        g2_TimePoint = c(
            1, 6, 1, 2, 8, 3, 1, 1, 6, 1, 2,
            8, 3, 1, 1, 6, 1, 2, 8, 3, 1, 1,
            6, 1, 2, 8, 3, 1
        ),
        shared = c(
            1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0,
            0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1,
            0, 0
        ),
        count_g1 = c(
            2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1,
            1, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1,
            1, 1
        ),
        count_g2 = c(
            2, 2, 1, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1,
            1, 2, 2, 1, 1, 2, 1, 1, 2, 2, 1, 1, 2,
            1, 1
        ),
        count_union = c(
            3, 3, 2, 3, 4, 2, 3, 3, 3, 2, 1, 2,
            2, 2, 3, 3, 2, 3, 4, 2, 3, 3, 3, 2,
            1, 2, 2, 2
        ),
        on_g1 = c(
            50, 50, 50, 0, 0, 50, 0, 0, 0, 0, 100, 100,
            0, 0, 50, 50, 50, 0, 0, 50, 0, 0, 0, 0,
            100, 100, 0, 0
        ),
        on_g2 = c(
            50, 50, 100, 0, 0, 100, 0, 0, 0, 0, 100,
            50, 0, 0, 50, 50, 100, 0, 0, 100, 0, 0, 0,
            0, 100, 50, 0, 0
        ),
        on_union = c(
            33.3333333333333, 33.3333333333333, 50,
            0, 0, 50, 0, 0, 0, 0, 100, 50, 0, 0,
            33.3333333333333, 33.3333333333333,
            50, 0, 0, 50, 0, 0, 0, 0, 100, 50,
            0, 0
        )
    ) %>%
        dplyr::arrange(.data$g1),
    PT02 = data.table::data.table(
        g1 = c(
            "PT02_CD34_BM_1", "PT02_CD34_BM_1", "PT02_CD34_BM_1",
            "PT02_CD34_BM_1", "PT02_CD34_BM_1", "PT02_CD34_BM_1",
            "PT02_CD34_BM_3", "PT02_CD34_BM_3", "PT02_CD34_BM_3",
            "PT02_CD34_BM_3", "PT02_CD34_BM_3", "PT02_CD34_BM_3",
            "PT02_Whole_BM_1", "PT02_Whole_BM_1", "PT02_Whole_BM_1",
            "PT02_Whole_BM_1", "PT02_Whole_BM_1", "PT02_Whole_BM_1",
            "PT02_Whole_BM_6", "PT02_Whole_BM_6", "PT02_Whole_BM_6",
            "PT02_Whole_BM_6", "PT02_Whole_BM_6", "PT02_Whole_BM_6"
        ),
        g1_SubjectID = c(
            "PT02", "PT02", "PT02", "PT02", "PT02", "PT02",
            "PT02", "PT02", "PT02", "PT02", "PT02", "PT02",
            "PT02", "PT02", "PT02", "PT02", "PT02", "PT02",
            "PT02", "PT02", "PT02", "PT02", "PT02", "PT02"
        ),
        g1_CellMarker = c(
            "CD34", "CD34", "CD34", "CD34", "CD34", "CD34",
            "CD34", "CD34", "CD34", "CD34", "CD34", "CD34",
            "Whole", "Whole", "Whole", "Whole", "Whole",
            "Whole", "Whole", "Whole", "Whole", "Whole",
            "Whole", "Whole"
        ),
        g1_Tissue = c(
            "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM",
            "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM",
            "BM", "BM", "BM", "BM", "BM", "BM"
        ),
        g1_TimePoint = c(
            1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1,
            6, 6, 6, 6, 6, 6
        ),
        g2 = c(
            "PT02_CD14_PB_03", "PT02_CD13_PB_03", "PT02_CD13_PB_01",
            "PT02_CD15_PB_01", "PT02_CD15_PB_06", "PT02_CD14_PB_02",
            "PT02_CD14_PB_03", "PT02_CD13_PB_03", "PT02_CD13_PB_01",
            "PT02_CD15_PB_01", "PT02_CD15_PB_06", "PT02_CD14_PB_02",
            "PT02_CD14_PB_03", "PT02_CD13_PB_03", "PT02_CD13_PB_01",
            "PT02_CD15_PB_01", "PT02_CD15_PB_06", "PT02_CD14_PB_02",
            "PT02_CD14_PB_03", "PT02_CD13_PB_03", "PT02_CD13_PB_01",
            "PT02_CD15_PB_01", "PT02_CD15_PB_06", "PT02_CD14_PB_02"
        ),
        g2_SubjectID = c(
            "PT02", "PT02", "PT02", "PT02", "PT02", "PT02",
            "PT02", "PT02", "PT02", "PT02", "PT02", "PT02",
            "PT02", "PT02", "PT02", "PT02", "PT02", "PT02",
            "PT02", "PT02", "PT02", "PT02", "PT02", "PT02"
        ),
        g2_CellMarker = c(
            "CD14", "CD13", "CD13", "CD15", "CD15", "CD14",
            "CD14", "CD13", "CD13", "CD15", "CD15", "CD14",
            "CD14", "CD13", "CD13", "CD15", "CD15", "CD14",
            "CD14", "CD13", "CD13", "CD15", "CD15", "CD14"
        ),
        g2_Tissue = c(
            "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB",
            "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB",
            "PB", "PB", "PB", "PB", "PB", "PB"
        ),
        g2_TimePoint = c(
            3, 3, 1, 1, 6, 2, 3, 3, 1, 1, 6, 2, 3, 3, 1, 1, 6, 2,
            3, 3, 1, 1, 6, 2
        ),
        shared = c(
            0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0
        ),
        count_g1 = c(
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1
        ),
        count_g2 = c(
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1
        ),
        count_union = c(
            2, 2, 1, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2,
            2, 2, 2, 1, 2
        ),
        on_g1 = c(
            0, 0, 100, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0,
            0, 0, 0, 100, 0
        ),
        on_g2 = c(
            0, 0, 100, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0,
            0, 0, 0, 100, 0
        ),
        on_union = c(
            0, 0, 100, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0,
            0, 0, 0, 0, 100, 0
        )
    ) %>%
        dplyr::arrange(.data$g1)
)

expected_res2 <- data.table::data.table(
    g1 = c("PT01_CD34_BM_1", "PT01_CD34_BM_2", "PT01_CD34_BM_1",
           "PT01_CD34_BM_2", "PT01_CD34_BM_1", "PT01_CD34_BM_2",
           "PT01_CD34_BM_1", "PT01_CD34_BM_2", "PT01_CD34_BM_1",
           "PT01_CD34_BM_2", "PT01_CD34_BM_1", "PT01_CD34_BM_2",
           "PT01_CD34_BM_1", "PT01_CD34_BM_2", "PT01_CD34_BM_1",
           "PT01_CD34_BM_2", "PT01_CD34_BM_1", "PT01_CD34_BM_2",
           "PT01_CD34_BM_1", "PT01_CD34_BM_2", "PT01_CD34_BM_1",
           "PT01_CD34_BM_2", "PT01_CD34_BM_1", "PT01_CD34_BM_2",
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
           "PT02_Whole_BM_6", "PT02_Whole_BM_1"),
    g1_SubjectID = c("PT01", "PT01", "PT01", "PT01", "PT01", "PT01",
                     "PT01", "PT01", "PT01", "PT01", "PT01", "PT01",
                     "PT01", "PT01", "PT01", "PT01", "PT01", "PT01",
                     "PT01", "PT01", "PT01", "PT01", "PT01", "PT01",
                     "PT01", "PT01", "PT01", "PT01", "PT01", "PT01",
                     "PT01", "PT01", "PT01", "PT01", "PT01", "PT01",
                     "PT01", "PT01", "PT01", "PT01", "PT01", "PT01",
                     "PT01", "PT01", "PT01", "PT01", "PT01", "PT01",
                     "PT01", "PT01", "PT01", "PT01", "PT02", "PT02",
                     "PT02", "PT02", "PT02", "PT02", "PT02", "PT02",
                     "PT02", "PT02", "PT02", "PT02", "PT02", "PT02",
                     "PT02", "PT02", "PT02", "PT02", "PT02", "PT02",
                     "PT02", "PT02", "PT02", "PT02", "PT02", "PT02",
                     "PT02", "PT02", "PT02", "PT02", "PT02", "PT02",
                     "PT02", "PT02", "PT02", "PT02", "PT02", "PT02",
                     "PT02", "PT02", "PT02", "PT02", "PT02", "PT02",
                     "PT02", "PT02", "PT02", "PT02", "PT02", "PT02",
                     "PT02", "PT02"),
    g1_CellMarker = c("CD34", "CD34", "CD34", "CD34", "CD34", "CD34",
                      "CD34", "CD34", "CD34", "CD34", "CD34", "CD34",
                      "CD34", "CD34", "CD34", "CD34", "CD34", "CD34",
                      "CD34", "CD34", "CD34", "CD34", "CD34", "CD34",
                      "CD34", "CD34", "Whole", "Whole", "Whole", "Whole",
                      "Whole", "Whole", "Whole", "Whole", "Whole", "Whole",
                      "Whole", "Whole", "Whole", "Whole", "Whole", "Whole",
                      "Whole", "Whole", "Whole", "Whole", "Whole", "Whole",
                      "Whole", "Whole", "Whole", "Whole", "CD34", "CD34",
                      "CD34", "CD34", "CD34", "CD34", "CD34", "CD34", "CD34",
                      "CD34", "CD34", "CD34", "CD34", "CD34", "CD34", "CD34",
                      "CD34", "CD34", "CD34", "CD34", "CD34", "CD34", "CD34",
                      "CD34", "CD34", "CD34", "Whole", "Whole", "Whole",
                      "Whole", "Whole", "Whole", "Whole", "Whole", "Whole",
                      "Whole", "Whole", "Whole", "Whole", "Whole", "Whole",
                      "Whole", "Whole", "Whole", "Whole", "Whole", "Whole",
                      "Whole", "Whole", "Whole", "Whole", "Whole"),
    g1_Tissue = c("BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM",
                  "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM",
                  "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM",
                  "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM",
                  "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM",
                  "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM",
                  "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM",
                  "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM",
                  "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM",
                  "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM",
                  "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM", "BM",
                  "BM", "BM", "BM", "BM", "BM"),
    g1_TimePoint = c(1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1,
                     2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2,
                     1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1,
                     2, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1,
                     3, 1, 3, 1, 3, 1, 3, 1, 3, 1, 6, 1, 6, 1, 6, 1, 6,
                     1, 6, 1, 6, 1, 6, 1, 6, 1, 6, 1, 6, 1, 6, 1, 6, 1, 6, 1),
    g2 = c("PT02_CD14_PB_03", "PT02_CD14_PB_03", "PT02_CD13_PB_03",
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
           "PT02_CD14_PB_02", "PT02_CD14_PB_02"),
    g2_SubjectID = c("PT02", "PT02", "PT02", "PT02", "PT01", "PT01",
                     "PT01", "PT01", "PT01", "PT01", "PT02", "PT02",
                     "PT01", "PT01", "PT01", "PT01", "PT01", "PT01",
                     "PT02", "PT02", "PT02", "PT02", "PT01", "PT01",
                     "PT02", "PT02", "PT02", "PT02", "PT02", "PT02",
                     "PT01", "PT01", "PT01", "PT01", "PT01", "PT01",
                     "PT02", "PT02", "PT01", "PT01", "PT01", "PT01",
                     "PT01", "PT01", "PT02", "PT02", "PT02", "PT02",
                     "PT01", "PT01", "PT02", "PT02", "PT02", "PT02",
                     "PT02", "PT02", "PT01", "PT01", "PT01", "PT01",
                     "PT01", "PT01", "PT02", "PT02", "PT01", "PT01",
                     "PT01", "PT01", "PT01", "PT01", "PT02", "PT02",
                     "PT02", "PT02", "PT01", "PT01", "PT02", "PT02",
                     "PT02", "PT02", "PT02", "PT02", "PT01", "PT01",
                     "PT01", "PT01", "PT01", "PT01", "PT02", "PT02",
                     "PT01", "PT01", "PT01", "PT01", "PT01", "PT01",
                     "PT02", "PT02", "PT02", "PT02", "PT01", "PT01",
                     "PT02", "PT02"),
    g2_CellMarker = c("CD14", "CD14", "CD13", "CD13", "CD15", "CD15",
                      "CD14", "CD14", "CD14", "CD14", "CD13", "CD13",
                      "CD14", "CD14", "CD13", "CD13", "CD14", "CD14",
                      "CD15", "CD15", "CD15", "CD15", "CD13", "CD13",
                      "CD14", "CD14", "CD14", "CD14", "CD13", "CD13",
                      "CD15", "CD15", "CD14", "CD14", "CD14", "CD14",
                      "CD13", "CD13", "CD14", "CD14", "CD13", "CD13",
                      "CD14", "CD14", "CD15", "CD15", "CD15", "CD15",
                      "CD13", "CD13", "CD14", "CD14", "CD14", "CD14",
                      "CD13", "CD13", "CD15", "CD15", "CD14", "CD14",
                      "CD14", "CD14", "CD13", "CD13", "CD14", "CD14",
                      "CD13", "CD13", "CD14", "CD14", "CD15", "CD15",
                      "CD15", "CD15", "CD13", "CD13", "CD14", "CD14",
                      "CD14", "CD14", "CD13", "CD13", "CD15", "CD15",
                      "CD14", "CD14", "CD14", "CD14", "CD13", "CD13",
                      "CD14", "CD14", "CD13", "CD13", "CD14", "CD14",
                      "CD15", "CD15", "CD15", "CD15", "CD13", "CD13",
                      "CD14", "CD14"),
    g2_Tissue = c("PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB",
                  "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB",
                  "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB",
                  "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB",
                  "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB",
                  "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB",
                  "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB",
                  "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB",
                  "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB",
                  "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB",
                  "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB", "PB",
                  "PB", "PB", "PB", "PB", "PB"),
    g2_TimePoint = c(3, 3, 3, 3, 1, 1, 6, 6, 1, 1, 1, 1, 2, 2, 8, 8, 3, 3,
                     1, 1, 6, 6, 1, 1, 2, 2, 3, 3, 3, 3, 1, 1, 6, 6, 1, 1,
                     1, 1, 2, 2, 8, 8, 3, 3, 1, 1, 6, 6, 1, 1, 2, 2, 3, 3,
                     3, 3, 1, 1, 6, 6, 1, 1, 1, 1, 2, 2, 8, 8, 3, 3, 1, 1,
                     6, 6, 1, 1, 2, 2, 3, 3, 3, 3, 1, 1, 6, 6, 1, 1, 1, 1,
                     2, 2, 8, 8, 3, 3, 1, 1, 6, 6, 1, 1, 2, 2),
    shared = c(0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0,
               1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1,
               0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0,
               0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0,
               0, 0, 0, 0),
    count_g1 = c(2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2,
                 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1,
                 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 1, 1, 1, 1),
    count_g2 = c(1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1,
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1,
                 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2,
                 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 2, 1,
                 1, 1, 1, 1, 1, 1, 1, 1, 1),
    count_union = c(3, 2, 3, 2, 3, 3, 3, 3, 2, 2, 3, 2, 3, 1, 4, 2, 2, 2,
                    3, 2, 2, 2, 3, 2, 3, 1, 3, 2, 3, 2, 3, 3, 3, 3, 2, 2,
                    3, 2, 3, 1, 4, 2, 2, 2, 3, 2, 2, 2, 3, 2, 3, 1, 2, 2,
                    1, 2, 3, 3, 3, 3, 2, 2, 2, 1, 2, 2, 3, 3, 2, 2, 2, 2,
                    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 3, 2, 2, 2, 2,
                    2, 2, 3, 3, 1, 2, 2, 1, 1, 2, 2, 2, 2, 2),
    on_g1 = c(0, 0, 0, 0, 50, 0, 50, 0, 50, 0, 0, 0, 0, 100, 0, 100, 50, 0,
              0, 0, 50, 0, 0, 0, 0, 100, 0, 0, 0, 0, 50, 0, 50, 0, 50, 0, 0,
              0, 0, 100, 0, 100, 50, 0, 0, 0, 50, 0, 0, 0, 0, 100, 0, 0, 100,
              0, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 100, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              100, 0, 0, 100, 100, 0, 0, 0, 0, 0),
    on_g2 = c(0, 0, 0, 0, 50, 0, 50, 0, 100, 0, 0, 0, 0, 100, 0, 50, 100,
              0, 0, 0, 100, 0, 0, 0, 0, 100, 0, 0, 0, 0, 50, 0, 50, 0, 100,
              0, 0, 0, 0, 100, 0, 50, 100, 0, 0, 0, 100, 0, 0, 0, 0, 100,
              0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0, 50, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 100, 0, 0, 100, 100, 0, 0, 0, 0, 0),
    on_union = c(0, 0, 0, 0, 33.3333333333333, 0, 33.3333333333333, 0, 50,
                 0, 0, 0, 0, 100, 0, 50, 50, 0, 0, 0, 50, 0, 0, 0, 0, 100,
                 0, 0, 0, 0, 33.3333333333333, 0, 33.3333333333333, 0, 50,
                 0, 0, 0, 0, 100, 0, 50, 50, 0, 0, 0, 50, 0, 0, 0, 0, 100,
                 0, 0, 100, 0, 0, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 0, 50, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 100, 0, 0, 100, 100, 0, 0, 0, 0, 0)
)

expected_intra_marker <- list(
  PT01 = data.table::data.table(
    g1 = c('PT01_CD34_BM_1', 'PT01_CD34_BM_1', 'PT01_CD34_BM_1',
           'PT01_CD34_BM_1', 'PT01_CD34_BM_1', 'PT01_CD34_BM_2',
           'PT01_CD34_BM_2', 'PT01_CD34_BM_2', 'PT01_CD34_BM_2',
           'PT01_CD34_BM_2', 'PT01_CD34_BM_3', 'PT01_CD34_BM_3',
           'PT01_CD34_BM_3', 'PT01_CD34_BM_3', 'PT01_CD34_BM_3'),
    g1_SubjectID = c('PT01', 'PT01', 'PT01', 'PT01', 'PT01', 'PT01',
                     'PT01', 'PT01', 'PT01', 'PT01', 'PT01', 'PT01',
                     'PT01', 'PT01', 'PT01'),
    g1_CellMarker = c('CD34', 'CD34', 'CD34', 'CD34', 'CD34', 'CD34',
                      'CD34', 'CD34', 'CD34', 'CD34', 'CD34', 'CD34',
                      'CD34', 'CD34', 'CD34'),
    g1_Tissue = c('BM', 'BM', 'BM', 'BM', 'BM', 'BM', 'BM', 'BM', 'BM',
                  'BM', 'BM', 'BM', 'BM', 'BM', 'BM'),
    g1_TimePoint = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3),
    g2 = c('PT01_CD34_BM_03', 'PT01_CD34_BM_01', 'PT01_CD34_BM_06',
           'PT01_CD34_BM_02', 'PT01_CD34_BM_08', 'PT01_CD34_BM_03',
           'PT01_CD34_BM_01', 'PT01_CD34_BM_06', 'PT01_CD34_BM_02',
           'PT01_CD34_BM_08', 'PT01_CD34_BM_03', 'PT01_CD34_BM_01',
           'PT01_CD34_BM_06', 'PT01_CD34_BM_02', 'PT01_CD34_BM_08'),
    g2_SubjectID = c('PT01', 'PT01', 'PT01', 'PT01', 'PT01', 'PT01',
                     'PT01', 'PT01', 'PT01', 'PT01', 'PT01', 'PT01',
                     'PT01', 'PT01', 'PT01'),
    g2_CellMarker = c('CD34', 'CD34', 'CD34', 'CD34', 'CD34', 'CD34',
                      'CD34', 'CD34', 'CD34', 'CD34', 'CD34', 'CD34',
                      'CD34', 'CD34', 'CD34'),
    g2_Tissue = c('BM', 'BM', 'BM', 'BM', 'BM', 'BM', 'BM', 'BM', 'BM',
                  'BM', 'BM', 'BM', 'BM', 'BM', 'BM'),
    g2_TimePoint = c(3, 1, 6, 2, 8, 3, 1, 6, 2, 8, 3, 1, 6, 2, 8),
    shared = c(1, 3, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0),
    count_g1 = c(3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    count_g2 = c(2, 3, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 1, 1),
    count_union = c(4, 3, 3, 4, 4, 3, 4, 2, 1, 1, 2, 4, 2, 2, 2),
    on_g1 = c(33.3333333333333, 100, 33.3333333333333, 0, 0, 0, 0, 0, 100,
              100, 100, 0, 0, 0, 0),
    on_g2 = c(50, 100, 100, 0, 0, 0, 0, 0, 100, 100, 50, 0, 0, 0, 0),
    on_union = c(25, 100, 33.3333333333333, 0, 0, 0, 0, 0, 100, 100, 50, 0,
                 0, 0, 0)
))

test_that("iss_source produces expected output - per patient", {
    res <- iss_source(test_df2, test_df3) %>%
        purrr::map(~ {
            .x %>%
                dplyr::select(-.data$is_coord) %>%
                dplyr::arrange(.data$g1)
        })
    expect_equal(res, expected_res1)
})

test_that("iss_source behaves right with intra-marker", {
  sub <- test_df2 %>%
    dplyr::filter(.data$CellMarker == "CD34") %>%
    dplyr::mutate(SubjectID = "PT01")
  res <- iss_source(sub, sub) %>%
    purrr::map(~ {
      .x %>%
        dplyr::select(-.data$is_coord) %>%
        dplyr::arrange(.data$g1)
    })
  expect_equal(res, expected_intra_marker)
})

test_that("iss_source produces expected output - NO per patient", {
    res <- iss_source(test_df2, test_df3, by_subject = FALSE) %>%
      select(-is_coord)
    expect_equal(res, expected_res2)
})
