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
    Value = c(3, 10, 40, 2, 15, 150)
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
                ),
                CompleteAmplificationID = c(
                    "ID1",
                    "ID2",
                    "ID1",
                    "ID1",
                    "ID3",
                    "ID2"
                )
            ))
        },
        regexp = .missing_value_col_error()
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
        AbsAbundance = c(
            0.06667, 0.0625, 0.88889,
            0.04445, 1.0, 0.9375
        ),
        PercAbundance = AbsAbundance * 100
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
