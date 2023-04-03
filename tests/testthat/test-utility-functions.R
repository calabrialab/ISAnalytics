withr::local_options(list(ISAnalytics.verbose = FALSE))

#------------------------------------------------------------------------------#
# Tests .process_af_for_gen
#------------------------------------------------------------------------------#
test_that(".process_af_for_gen works for default af", {
    res <- .process_af_for_gen(af = "default")
    af_out <- res$af
    tag_out <- res$tag_list
    expected_tags <- list(
        "pcr_repl_id" = "CompleteAmplificationID",
        "project_id" = "ProjectID",
        "tag_seq" = "TagSequence",
        "vispa_concatenate" = "concatenatePoolIDSeqRun"
    )
    af_sym <- "association_file"
    utils::data(list = af_sym, envir = rlang::current_env())
    association_file <- rlang::eval_tidy(rlang::sym(af_sym))
    expect_equal(tag_out, expected_tags)
    expect_equal(af_out, association_file)
})

test_that(".process_af_for_gen works for custom af", {
    af_custom <- tibble::tribble(
        ~ProjectID, ~TagSequence, ~CompleteAmplificationID,
        ~concatenatePoolIDSeqRun, ~X,
        "PJ01", "LTR75LC38", "ID1", "POOL1", "a",
        "PJ01", "LTR53LC32", "ID2", "POOL1", "b",
        "PJ01", "LTR27LC94", "ID3", "POOL2", "c",
        "PJ01", "LTR69LC52", "ID4", "POOL3", "d"
    )
    res <- .process_af_for_gen(af = af_custom)
    af_out <- res$af
    tag_out <- res$tag_list
    expected_tags <- list(
        "pcr_repl_id" = "CompleteAmplificationID",
        "project_id" = "ProjectID",
        "tag_seq" = "TagSequence",
        "vispa_concatenate" = "concatenatePoolIDSeqRun"
    )
    expect_equal(tag_out, expected_tags)
    expect_equal(af_out, af_custom)

    af_custom_err <- af_custom |>
        dplyr::select(-dplyr::all_of(c("ProjectID")))
    expect_error(
        {
            res <- .process_af_for_gen(af = af_custom_err)
        },
        class = "missing_req_col_err"
    )
})

#------------------------------------------------------------------------------#
# Tests .process_m_for_gen
#------------------------------------------------------------------------------#
test_that(".process_m_for_gen works with default matrices", {
    af <- .process_af_for_gen(af = "default")
    tags <- af$tag_list
    af <- af$af
    res <- .process_m_for_gen(matrices = "default", af, tags)
    expect_equal(names(res), "PJ01")
    expect_equal(names(res$PJ01), c(
        "POOL01-1", "POOL02-1",
        "POOL03-1", "POOL04-1"
    ))
})

#------------------------------------------------------------------------------#
# Tests .process_iss_for_gen
#------------------------------------------------------------------------------#
test_that(".process_iss_for_gen works with default stats", {
    af <- .process_af_for_gen(af = "default")
    tags <- af$tag_list
    af <- af$af
    res <- .process_iss_for_gen(af, tags)
    expect_equal(names(res), "PJ01")
    expect_equal(names(res$PJ01), c(
        "POOL01-1", "POOL02-1",
        "POOL03-1", "POOL04-1"
    ))
})

#------------------------------------------------------------------------------#
# Tests .generate_correct
#------------------------------------------------------------------------------#
test_that(".generate_correct works correctly with defaults", {
    af <- .process_af_for_gen(af = "default")
    tags <- af$tag_list
    af <- af$af
    mat <- .process_m_for_gen(matrices = "default", af, tags)
    iss <- .process_iss_for_gen(af, tags)
    temp_dir <- tempdir()
    expect_error(
        {
            path_out <- .generate_correct(temp_dir,
                sep_matrices = mat,
                sep_stats = iss
            )
        },
        regexp = NA
    )
    proj_fold <- fs::path(path_out, "PJ01")
    quant_fold <- fs::path(proj_fold, "quantification")
    iss_fold <- fs::path(proj_fold, "iss")
    pools <- paste0(paste0("POOL0", 1:4), "-1")
    pools_quant <- fs::path(quant_fold, pools)
    pools_iss <- fs::path(iss_fold, pools)
    expect_true(all(pools_quant %in% fs::dir_ls(path_out, recurse = TRUE)))
    expect_true(all(pools_iss %in% fs::dir_ls(path_out, recurse = TRUE)))
    expect_true(all(purrr::map_lgl(
        pools_quant,
        ~ length(fs::dir_ls(.x)) == 2
    )))
    expect_true(all(purrr::map_lgl(
        pools_iss,
        ~ length(fs::dir_ls(.x)) == 1
    )))
})

#------------------------------------------------------------------------------#
# Tests .generate_incorrect
#------------------------------------------------------------------------------#
test_that(".generate_incorrect works correctly with defaults", {
    af <- .process_af_for_gen(af = "default")
    tags <- af$tag_list
    af <- af$af
    mat <- .process_m_for_gen(matrices = "default", af, tags)
    iss <- .process_iss_for_gen(af, tags)
    temp_dir <- tempdir()
    expect_error(
        {
            path_out <- .generate_incorrect(temp_dir,
                sep_matrices = mat,
                sep_stats = iss
            )
        },
        regexp = NA
    )
    proj_fold <- fs::path(path_out, "PJ01")
    quant_fold <- fs::path(proj_fold, "quantification")
    iss_fold <- fs::path(proj_fold, "iss")
    pools <- c("POOL01-1-err", paste0(paste0("POOL0", 2:4), "-1"))
    pools_quant <- fs::path(quant_fold, pools)
    pools_iss <- fs::path(iss_fold, pools)
    expect_true(all(pools_quant %in% fs::dir_ls(path_out, recurse = TRUE)))
    expect_true(all(pools_iss %in% fs::dir_ls(path_out, recurse = TRUE)))
    expect_true(all(purrr::map_lgl(pools_quant, ~ {
        split <- fs::path_split(.x)
        pool_fold <- split[[1]][length(split[[1]])]
        if (pool_fold == "POOL02-1" & length(fs::dir_ls(.x)) == 1) {
            return(TRUE)
        }
        return(length(fs::dir_ls(.x)) == 2)
    })))
    expect_true(all(purrr::map_lgl(pools_iss, ~ {
        split <- fs::path_split(.x)
        pool_fold <- split[[1]][length(split[[1]])]
        if (pool_fold == "POOL02-1" & length(fs::dir_ls(.x)) == 0) {
            return(TRUE)
        }
        return(length(fs::dir_ls(.x)) == 1)
    })))
})

#------------------------------------------------------------------------------#
# Tests transform_columns
#------------------------------------------------------------------------------#
test_that("transform_columns applies transform as expected", {
    test_df <- tibble::tribble(
        ~A, ~B, ~C,
        1, 2, 3,
        4, 5, 6,
        7, 8, 9
    )
    to_apply <- list(A = ~ stringr::str_pad(
        as.character(.x),
        pad = "0", side = "left", width = 2
    ), B = ~ .x + 2)
    result <- transform_columns(test_df, to_apply)
    expect_equal(result$A, c("01", "04", "07"))
    expect_equal(result$B, c(4, 7, 10))
})

#------------------------------------------------------------------------------#
# Tests pcr_id_column
#------------------------------------------------------------------------------#
test_that("pcr_id_column errors if no pcr_repl_id in vars", {
    withr::local_options(list(ISAnalytics.af_specs = "default"))
    af_specs <- association_file_columns(TRUE)
    af_specs_mod <- af_specs |>
        dplyr::filter(!.data$tag %in% c("pcr_repl_id"))
    withr::with_options(list(ISAnalytics.af_specs = af_specs_mod), {
        expect_error({
            col <- pcr_id_column()
        })
    })
})

#------------------------------------------------------------------------------#
# Test comparison_matrix
#------------------------------------------------------------------------------#
m1 <- tibble::tribble(
    ~chr, ~integration_locus, ~strand, ~CompleteAmplificationID, ~Value,
    "1", 42522, "+", "ID1", 453,
    "3", 43521, "-", "ID1", 564,
    "2", 12425, "-", "ID2", 67,
    "1", 25756, "+", "ID3", 123,
    "7", 12345, "-", "ID3", 897,
)
m2 <- tibble::tribble(
    ~chr, ~integration_locus, ~strand, ~CompleteAmplificationID, ~Value,
    "1", 42522, "+", "ID1", 123.65,
    "3", 43521, "-", "ID1", 786.12,
    "2", 12425, "-", "ID2", 12.35,
    "1", 25756, "+", "ID3", 78.546,
    "7", 12345, "-", "ID3", 93.456
)
example <- list(seqCount = m1, fragmentEstimate = m2)

test_that("comparison_matrix throws error if input is incorrect", {
    # Input must be a named list, not single data frame
    expect_error({
        comp <- comparison_matrix(m1)
    })
    # List names must be quantification types
    expect_error({
        comp <- comparison_matrix(setNames(example, c("a", "b")))
    })
})

test_that("comparison_matrix notifies NA introduction", {
    withr::local_options(list(ISAnalytics.verbose = TRUE))
    mod_m2 <- m2[c(1, 2, 3), ]
    expect_message(
        {
            comp <- comparison_matrix(list(
                seqCount = m1,
                fragmentEstimate = mod_m2
            ))
        },
        class = "comp_nas"
    )
})

test_that("comparison_matrix produces correct output", {
    comp <- comparison_matrix(example)
    expected <- tibble::tribble(
        ~chr, ~integration_locus, ~strand, ~CompleteAmplificationID,
        ~seqCount, ~fragmentEstimate,
        "1", 42522, "+", "ID1", 453, 123.65,
        "3", 43521, "-", "ID1", 564, 786.12,
        "2", 12425, "-", "ID2", 67, 12.35,
        "1", 25756, "+", "ID3", 123, 78.546,
        "7", 12345, "-", "ID3", 897, 93.456
    )
    expect_equal(comp |> dplyr::arrange(.data$integration_locus),
        expected |> dplyr::arrange(.data$integration_locus),
        ignore_attr = TRUE
    )
})

test_that("comparison_matrix supports custom names", {
    comp <- comparison_matrix(example,
        seqCount = "sc",
        fragmentEstimate = "fe"
    )
    expect_true(all(c("fe", "sc") %in% colnames(comp)))
    expect_true(nrow(comp) == nrow(example$seqCount))
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
        key = c(
            mandatory_IS_vars(),
            "CompleteAmplificationID"
        )
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
# Tests generate_blank_association_file
#------------------------------------------------------------------------------#
test_that("generate_blank_association_file stops if path is not char", {
    expect_error({
        generate_blank_association_file(1)
    })
})

test_that("generate_blank_association_file works correctly", {
    temp <- withr::local_tempfile()
    generate_blank_association_file(temp)
    af <- read.csv(temp, sep = "\t", check.names = FALSE, header = TRUE)
    expect_true(all(colnames(af) %in% association_file_columns()))
})

#------------------------------------------------------------------------------#
# Tests generate_Vispa2_launch_AF
#------------------------------------------------------------------------------#
## Testing input
test_that("generate_Vispa2_launch_AF stops if association_file is not df", {
    expect_error({
        generate_Vispa2_launch_AF("af.tsv", "x", "y", "z")
    })
})

test_that("generate_Vispa2_launch_AF stops if project is not char", {
    expect_error({
        generate_Vispa2_launch_AF(association_file, 1, "y", "z")
    })
})

test_that("generate_Vispa2_launch_AF stops if pool is not char", {
    expect_error({
        generate_Vispa2_launch_AF(association_file, "PJ01", 1, "z")
    })
})

test_that("generate_Vispa2_launch_AF stops if lengths of projects and pool is
          not the same", {
    expect_error({
        generate_Vispa2_launch_AF(
            association_file, c("PJ01"),
            c("POOL01", "POOL02"), "z"
        )
    })
})

test_that("generate_Vispa2_launch_AF stops if af is malformed", {
    af <- association_file |>
        dplyr::select(-dplyr::all_of(c("ProjectID")))
    expect_error({
        generate_Vispa2_launch_AF(af, c("PJ01"), c("POOL01"), "test")
    })
})

## Testing output
test_that("generate_Vispa2_launch_AF works for single pair", {
    temp <- withr::local_tempdir()
    project <- c("PJ01")
    pool <- c("POOL01")
    name <- paste0(project, "-", pool, "_AF.tsv")
    complete_path <- fs::path(temp, name)
    generate_Vispa2_launch_AF(association_file, project, pool, temp)
    df <- read.csv(complete_path, sep = "\t", header = FALSE)
    expect_true(ncol(df) == 11)
    expect_equal(df[, 1], df[, 2])
    expect_message(
        {
            generate_Vispa2_launch_AF(
                association_file,
                c("PJ01"), c("POOL05"), temp
            )
        },
        class = "launch_af_empty"
    )
})

test_that("generate_Vispa2_launch_AF works for multiple pair", {
    temp <- withr::local_tempdir()
    project <- c("PJ01", "PJ01")
    pool <- c("POOL01", "POOL02")
    name <- paste0(project, "-", pool, "_AF.tsv")
    complete_path <- fs::path(temp, name)
    generate_Vispa2_launch_AF(association_file, project, pool, temp)
    df <- read.csv(complete_path[1], sep = "\t", header = FALSE)
    expect_true(ncol(df) == 11)
    expect_equal(df[, 1], df[, 2])
    df <- read.csv(complete_path[2], sep = "\t", header = FALSE)
    expect_true(ncol(df) == 11)
    expect_equal(df[, 1], df[, 2])
})

#------------------------------------------------------------------------------#
# Tests as_sparse_matrix
#------------------------------------------------------------------------------#
smpl <- tibble::tibble(
    chr = c(1, 2, 3), integration_locus = c(1354, 5634, 4765),
    strand = c("+", "+", "+"),
    GeneName = c("GENE1", "GENE2", "GENE3"),
    GeneStrand = c("+", "+", "+"),
    CompleteAmplificationID = c("ID1", "ID2", "ID3"),
    Value = c(46, 546, 587)
)

smpl_multi <- tibble::tibble(
    chr = c(1, 2, 3),
    integration_locus = c(1354, 5634, 4765),
    strand = c("+", "+", "+"),
    GeneName = c("GENE1", "GENE2", "GENE3"),
    GeneStrand = c("+", "+", "+"),
    CompleteAmplificationID = c("ID1", "ID2", "ID3"),
    seqCount = c(46, 546, 587),
    fragmentEstimate = c(4234.5, 533.45, 5431.43),
    barcodeCount = c(46, 6456, 456)
)

test_that("as_sparse_matrix works with single matrix", {
    sparse <- as_sparse_matrix(smpl)
    expected <- tibble::tibble(
        chr = c(1, 2, 3),
        integration_locus = c(1354, 5634, 4765),
        strand = c("+", "+", "+"),
        GeneName = c("GENE1", "GENE2", "GENE3"),
        GeneStrand = c("+", "+", "+"),
        ID1 = c(46, NA, NA),
        ID2 = c(NA, 546, NA),
        ID3 = c(NA, NA, 587)
    )
    expect_equal(sparse, expected, ignore_attr = TRUE)
})

test_that("as_sparse_matrix works with multi-quant matrix", {
    sparse <- as_sparse_matrix(smpl_multi)
    expect_true(is.list(sparse) & !is.data.frame(sparse))
    expect_true(all(c("seqCount", "fragmentEstimate", "barcodeCount") %in%
        names(sparse)))
    expect_true(all(c(
        "chr", "integration_locus", "strand",
        "GeneName", "GeneStrand", "ID1", "ID2", "ID3"
    ) %in% colnames(sparse$seqCount)))
    expect_equal(nrow(sparse$seqCount), 3)
    expect_true(all(c(
        "chr", "integration_locus", "strand",
        "GeneName", "GeneStrand", "ID1", "ID2", "ID3"
    ) %in% colnames(sparse$fragmentEstimate)))
    expect_equal(nrow(sparse$fragmentEstimate), 3)
    expect_true(all(c(
        "chr", "integration_locus", "strand",
        "GeneName", "GeneStrand", "ID1", "ID2", "ID3"
    ) %in% colnames(sparse$barcodeCount)))
    expect_equal(nrow(sparse$barcodeCount), 3)
})

test_that("as_sparse_matrix works with list of matrices", {
    sparse <- as_sparse_matrix(list(smpl, smpl))
    expect_true(!tibble::is_tibble(sparse) & is.list(sparse))
    expect_true(all(c(
        "chr", "integration_locus", "strand",
        "GeneName", "GeneStrand", "ID1", "ID2", "ID3"
    ) %in%
        colnames(sparse[[1]])))
    expect_equal(nrow(sparse[[1]]), 3)
    expect_true(all(c(
        "chr", "integration_locus", "strand",
        "GeneName", "GeneStrand", "ID1", "ID2", "ID3"
    ) %in%
        colnames(sparse[[2]])))
    expect_equal(nrow(sparse[[2]]), 3)
})

test_that("as_sparse_matrix works with aggreg matrix", {
    smpl_agg <- tibble::tibble(
        chr = c(1, 2, 3),
        integration_locus = c(1354, 5634, 4765),
        strand = c("+", "+", "+"),
        GeneName = c("GENE1", "GENE2", "GENE3"),
        GeneStrand = c("+", "+", "+"),
        SubjectID = c("S1", "S2", "S2"),
        CellMarker = c("C1", "C1", "C2"),
        Value_sum = c(46, 546, 587)
    )
    sparse <- as_sparse_matrix(smpl_agg,
        single_value_col = "Value_sum",
        key = c("SubjectID", "CellMarker")
    )
    expected <- tibble::tibble(
        chr = c(1, 2, 3),
        integration_locus = c(1354, 5634, 4765),
        strand = c("+", "+", "+"),
        GeneName = c("GENE1", "GENE2", "GENE3"),
        GeneStrand = c("+", "+", "+"),
        S1_C1 = c(46, NA, NA),
        S2_C1 = c(NA, 546, NA),
        S2_C2 = c(NA, NA, 587)
    )
    expect_equal(sparse, expected)
})

#------------------------------------------------------------------------------#
# Tests annotation_issues
#------------------------------------------------------------------------------#
test_df_issues <- tibble::tribble(
    ~chr, ~integration_locus, ~strand, ~GeneName, ~GeneStrand,
    ~CompleteAmplificationID, ~Value,
    "1", 123456, "+", "ABCDE", "+", "ID1", 56,
    "1", 123456, "+", "ABCDE", "-", "ID2", 675,
    "1", 123456, "+", "FGHIL", "-", "ID3", 67,
    "2", 5674653, "-", "FGHIL", "-", "ID2", 873,
    "1", 4578768, "-", "RSPQX", "-", "ID3", 983,
)

test_df_no_issues <- tibble::tribble(
    ~chr, ~integration_locus, ~strand, ~GeneName, ~GeneStrand,
    ~CompleteAmplificationID, ~Value,
    "1", 123456, "+", "ABCDE", "+", "ID1", 56,
    "1", 123456, "+", "ABCDE", "+", "ID2", 675,
    "1", 123456, "+", "ABCDE", "+", "ID3", 67,
    "2", 5674653, "-", "FGHIL", "-", "ID2", 873,
    "1", 4578768, "-", "RSPQX", "-", "ID3", 983,
)

test_that("annotation_issues returns df if issues", {
    res <- annotation_issues(test_df_issues)
    expect_true(!is.null(res))
    expect_true(nrow(res) == 1 & res$chr[1] == "1" &
        res$integration_locus[1] == 123456 & res$strand == "+" &
        res$distinct_genes == 3)
})

test_that("annotation_issues returns null if no issues", {
    withr::local_options(list(ISAnalytics.verbose = TRUE))
    expect_message({
        res <- annotation_issues(test_df_no_issues)
    })
    expect_null(res)
})

test_that("annotation_issues works with lists", {
    res <- annotation_issues(list(a = test_df_issues, b = test_df_no_issues))
    expect_true(!is.null(res))
    expect_true(is.null(res[[2]]))
    expect_true(nrow(res[[1]]) == 1)
})

#------------------------------------------------------------------------------#
# Tests inspect_tags
#------------------------------------------------------------------------------#
test_that("inspect_tags prints tag info", {
    expect_message(
        {
            inspect_tags("chromosome")
        },
        class = "tag_inspect"
    )
    expect_message(
        {
            expect_message(
                {
                    inspect_tags(c("chromosome", "locus"))
                },
                class = "tag_inspect"
            )
        },
        class = "tag_inspect"
    )
})

#------------------------------------------------------------------------------#
# Tests unzip_file_system
#------------------------------------------------------------------------------#
test_that("unzip_file_system signals deprecation", {
    expect_defunct({
        unzip_file_system("fs", "fs")
    })
})

#------------------------------------------------------------------------------#
# .split_df_in_chunks
#------------------------------------------------------------------------------#
test_that(".split_df_in_chunks correctly splits", {
    df_1 <- tibble::tribble(
        ~A, ~B,
        1, "a",
        2, "b",
        3, "c",
        4, "d",
        5, "e",
        6, "f"
    )
    split_df <- .split_df_in_chunks(df_1, 3)
    expect_equal(length(split_df), 3)
    expect_true(all(purrr::map_lgl(split_df, ~ nrow(.x) == 2)))

    df_2 <- tibble::tribble(
        ~A, ~B,
        1, "a",
        2, "b",
        3, "c",
        4, "d",
        5, "e",
        6, "f",
        7, "g",
        8, "h"
    )
    split_df <- .split_df_in_chunks(df_2, 3)
    expect_equal(length(split_df), 3)
    expect_true(all(purrr::map_lgl(split_df[c(1, 2)], ~ nrow(.x) == 3)))
    expect_true(nrow(split_df[[3]]) == 2)

    df_3 <- tibble::tribble(
        ~A, ~B,
        1, "a",
        2, "b",
        3, "c",
        4, "d",
        5, "e",
        6, "f",
        7, "g",
        8, "h",
        9, "i",
        10, "j",
        11, "k"
    )
    split_df <- .split_df_in_chunks(df_3, 7)
    expect_equal(length(split_df), 7)
    expect_true(all(purrr::map_lgl(split_df[c(1, 2, 3, 4)], ~ nrow(.x) == 2)))
    expect_true(all(purrr::map_lgl(split_df[c(5, 6, 7)], ~ nrow(.x) == 1)))

    rejoint <- purrr::reduce(split_df, dplyr::bind_rows) |>
        dplyr::distinct()

    expect_true(nrow(rejoint) == nrow(df_3))
})

#------------------------------------------------------------------------------#
# .execute_map_job
#------------------------------------------------------------------------------#
test_that(".execute_map_job does not stop on error", {
    data_list <- list(1, 3, "2", -1)
    fun_to_app <- function(x, progress) {
        suppressWarnings({
            res <- sqrt(x)
        })
        if (!is.null(progress)) {
            progress()
        }
        return(res)
    }
    computed <- .execute_map_job(data_list,
        fun_to_apply = fun_to_app,
        stop_on_error = FALSE, fun_args = list()
    )
    expect_true(!is.null(computed$err[[3]]))
    expect_true(all(is.null(c(
        computed$err[[1]], computed$err[[2]],
        computed$err[[4]]
    ))))
    expect_true(is.nan(computed$res[[4]]))
    expect_true(is.null(computed$res[[3]]))
})

test_that(".execute_map_job stops on error", {
    data_list <- list(1, 3, "2", -1)
    fun_to_app <- function(x, progress) {
        suppressWarnings({
            res <- sqrt(x)
        })
        if (!is.null(progress)) {
            progress()
        }
        return(res)
    }
    expect_error({
        computed <- .execute_map_job(data_list,
            fun_to_apply = fun_to_app,
            stop_on_error = TRUE, fun_args = list()
        )
    })
})
