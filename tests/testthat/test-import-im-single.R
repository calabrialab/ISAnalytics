library(ISAnalytics)
func_name <- "import_single_Vispa2Matrix"

#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
sample_df <- tibble::tribble(
    ~chr, ~integration_locus, ~strand, ~GeneName, ~GeneStrand, ~id1, ~id2,
    ~id3, ~id4,
    "1", 140546, "+", "GENE1", "-", 4, NA, NA, 1,
    "14", 43567, "-", "GENE2", "+", 231, NA, 2, NA,
    "5", 214676, "-", "GENE3", "-", NA, NA, NA, NA,
    "7", 66778, "-", "GENE4", "-", NA, 355, NA, NA,
    "1", 75687, "+", "GENE5", "+", NA, NA, NA, 65,
    "5", 64576, "+", "GENE6", "-", 1, 667, NA, NA,
    "20", 57587, "-", "GENE7", "-", NA, 13, 1, NA,
    "X", 457658, "+", "GENE8", "+", NA, NA, NA, 768
)

sample_df_old <- tidyr::unite(sample_df,
    col = "IS_genomicID",
    .data$chr, .data$integration_locus,
    .data$strand, sep = "_", remove = FALSE
)

sample_df_add_columns <- sample_df %>%
    dplyr::mutate(add1 = "a", add2 = 56)

#------------------------------------------------------------------------------#
# Tests
#------------------------------------------------------------------------------#
## --- Test failures
test_that(paste(func_name, "fails if path doesn't exist"), {
    expect_error(
        {
            df <- import_single_Vispa2Matrix("my_file")
        },
        regexp = paste("File not found at", "my_file")
    )
})

test_that(paste(func_name, "fails if path is a dir"), {
    tmpdir <- withr::local_tempdir()
    expect_error(
        {
            df <- import_single_Vispa2Matrix(tmpdir)
        },
        regexp = paste("Path exists but is not a file")
    )
})

test_that(paste(func_name, "stops if malformed"), {
    tmp <- withr::local_tempfile(fileext = ".tsv")
    readr::write_tsv(sample_df_old, tmp)
    expect_error(
        {
            df <- import_single_Vispa2Matrix(tmp)
        },
        class = "malformed_ism"
    )
})

## --- Test different matrices types
test_that(paste(func_name, "reads type NEW standard"), {
    tf <- withr::local_tempfile(fileext = ".tsv")
    readr::write_tsv(sample_df, tf)
    expected_summary_msg <- .summary_ism_import_msg(
        "NEW",
        TRUE,
        c(8, 9),
        "fread"
    )
    expected_summary_msg <- c(
        expected_summary_msg[1],
        paste(
            "*",
            expected_summary_msg[seq(
                from = 2,
                to = length(expected_summary_msg)
            )]
        )
    )
    expected_summary_msg <- paste0(expected_summary_msg, collapse = "\n")
    capture_output(
        suppressMessages(
            expect_message(
                expect_message(
                    withr::with_file(
                        file = tf,
                        code = {
                            df <- import_single_Vispa2Matrix(tf)
                        }
                    )
                ),
                class = "ism_import_summary",
                regexp = expected_summary_msg,
                fixed = TRUE
            )
        )
    )
    expect_equal(dim(df), c(11, 7))
    expect_type(df$chr, "character")
    expect_type(df$strand, "character")
    expect_type(df$GeneName, "character")
    expect_type(df$GeneStrand, "character")
    expect_type(df$integration_locus, "integer")
    expect_true(typeof(df$Value) %in% c("double", "integer"))
})

test_that(paste(func_name, "reads type NEW different params"), {
    tf <- withr::local_tempfile(fileext = ".tsv")
    readr::write_tsv(sample_df, tf)
    expected_summary_msg <- .summary_ism_import_msg(
        "NEW",
        TRUE,
        c(8, 7),
        "fread"
    )
    expected_summary_msg <- c(
        expected_summary_msg[1],
        paste(
            "*",
            expected_summary_msg[seq(
                from = 2,
                to = length(expected_summary_msg)
            )]
        )
    )
    expected_summary_msg <- paste0(expected_summary_msg, collapse = "\n")
    capture_output(
        suppressMessages(
            expect_message(
                expect_message(
                    withr::with_file(
                        file = tf,
                        code = {
                            df <- import_single_Vispa2Matrix(tf,
                                to_exclude = c("id1", "id2")
                            )
                        }
                    )
                ),
                class = "ism_import_summary",
                regexp = expected_summary_msg,
                fixed = TRUE
            )
        )
    )
    expect_equal(dim(df), c(5, 7))
    expect_type(df$chr, "character")
    expect_type(df$strand, "character")
    expect_type(df$GeneName, "character")
    expect_type(df$GeneStrand, "character")
    expect_type(df$integration_locus, "integer")
    expect_true(typeof(df$Value) %in% c("double", "integer"))
    expect_true(all(df$CompleteAmplificationID %in% c("id3", "id4")))

    readr::write_tsv(sample_df %>% dplyr::select(-GeneName, -GeneStrand), tf)
    expected_summary_msg <- .summary_ism_import_msg(
        "NEW",
        FALSE,
        c(8, 7),
        "fread"
    )
    expected_summary_msg <- c(
        expected_summary_msg[1],
        paste(
            "*",
            expected_summary_msg[seq(
                from = 2,
                to = length(expected_summary_msg)
            )]
        )
    )
    expected_summary_msg <- paste0(expected_summary_msg, collapse = "\n")
    capture_output(
        suppressMessages(
            expect_message(
                expect_message(
                    withr::with_file(
                        file = tf,
                        code = {
                            df <- import_single_Vispa2Matrix(tf)
                        }
                    )
                ),
                class = "ism_import_summary",
                regexp = expected_summary_msg,
                fixed = TRUE
            )
        )
    )
    expect_equal(dim(df), c(11, 5))
    expect_type(df$chr, "character")
    expect_type(df$strand, "character")
    expect_type(df$integration_locus, "integer")
    expect_true(typeof(df$Value) %in% c("double", "integer"))
})

test_that(paste(func_name, "reads type OLD"), {
    tf <- withr::local_tempfile(fileext = ".tsv")
    readr::write_tsv(sample_df_old %>%
        dplyr::select(-chr, -integration_locus, -strand), tf)
    expected_summary_msg <- .summary_ism_import_msg(
        "OLD",
        TRUE,
        c(8, 7),
        "fread"
    )
    expected_summary_msg <- c(
        expected_summary_msg[1],
        paste(
            "*",
            expected_summary_msg[seq(
                from = 2,
                to = length(expected_summary_msg)
            )]
        )
    )
    expected_summary_msg <- paste0(expected_summary_msg, collapse = "\n")
    capture_output(
        suppressMessages(
            expect_message(
                expect_message(
                    withr::with_file(
                        file = tf,
                        code = {
                            df <- import_single_Vispa2Matrix(tf)
                        }
                    )
                ),
                class = "ism_import_summary",
                regexp = expected_summary_msg,
                fixed = TRUE
            )
        )
    )
    expect_equal(dim(df), c(11, 7))
    expect_type(df$chr, "character")
    expect_type(df$strand, "character")
    expect_type(df$GeneName, "character")
    expect_type(df$GeneStrand, "character")
    expect_type(df$integration_locus, "integer")
    expect_true(typeof(df$Value) %in% c("double", "integer"))
})

## --- Test different delimiters
test_that(paste(func_name, "reads comma delimited"), {
    tf <- withr::local_tempfile(fileext = ".csv")
    readr::write_csv(sample_df, tf)
    expected_summary_msg <- .summary_ism_import_msg(
        "NEW",
        TRUE,
        c(8, 9),
        "fread"
    )
    expected_summary_msg <- c(
        expected_summary_msg[1],
        paste(
            "*",
            expected_summary_msg[seq(
                from = 2,
                to = length(expected_summary_msg)
            )]
        )
    )
    expected_summary_msg <- paste0(expected_summary_msg, collapse = "\n")
    capture_output(
        suppressMessages(
            expect_message(
                expect_message(
                    withr::with_file(
                        file = tf,
                        code = {
                            df <- import_single_Vispa2Matrix(tf,
                                separator = ","
                            )
                        }
                    )
                ),
                class = "ism_import_summary",
                regexp = expected_summary_msg,
                fixed = TRUE
            )
        )
    )
    expect_equal(dim(df), c(11, 7))
    expect_type(df$chr, "character")
    expect_type(df$strand, "character")
    expect_type(df$GeneName, "character")
    expect_type(df$GeneStrand, "character")
    expect_type(df$integration_locus, "integer")
    expect_true(typeof(df$Value) %in% c("double", "integer"))
})

test_that(paste(func_name, "reads semicolon delimited"), {
    tf <- withr::local_tempfile(fileext = ".csv")
    readr::write_csv2(sample_df, tf)
    expected_summary_msg <- .summary_ism_import_msg(
        "NEW",
        TRUE,
        c(8, 9),
        "fread"
    )
    expected_summary_msg <- c(
        expected_summary_msg[1],
        paste(
            "*",
            expected_summary_msg[seq(
                from = 2,
                to = length(expected_summary_msg)
            )]
        )
    )
    expected_summary_msg <- paste0(expected_summary_msg, collapse = "\n")
    capture_output(
        suppressMessages(
            expect_message(
                expect_message(
                    withr::with_file(
                        file = tf,
                        code = {
                            df <- import_single_Vispa2Matrix(tf,
                                separator = ";"
                            )
                        }
                    )
                ),
                class = "ism_import_summary",
                regexp = expected_summary_msg,
                fixed = TRUE
            )
        )
    )
    expect_equal(dim(df), c(11, 7))
    expect_type(df$chr, "character")
    expect_type(df$strand, "character")
    expect_type(df$GeneName, "character")
    expect_type(df$GeneStrand, "character")
    expect_type(df$integration_locus, "integer")
    expect_true(typeof(df$Value) %in% c("double", "integer"))
})

## --- Test reading compressed files
test_that(paste(func_name, "reads correctly compressed .xz"), {
    tf <- withr::local_tempfile(fileext = ".tsv.xz")
    readr::write_tsv(sample_df, tf)
    expected_summary_msg <- .summary_ism_import_msg(
        "NEW",
        TRUE,
        c(8, 9),
        "classic"
    )
    expected_summary_msg <- c(
        expected_summary_msg[1],
        paste(
            "*",
            expected_summary_msg[seq(
                from = 2,
                to = length(expected_summary_msg)
            )]
        )
    )
    expected_summary_msg <- paste0(expected_summary_msg, collapse = "\n")
    capture_output(
        suppressMessages(
            expect_message(
                expect_message(
                    expect_message(withr::with_file(
                        file = tf,
                        code = {
                            df <- import_single_Vispa2Matrix(tf)
                        }
                    ), class = "unsup_comp_format")
                ),
                class = "ism_import_summary",
                regexp = expected_summary_msg,
                fixed = TRUE
            )
        )
    )
    expect_equal(dim(df), c(11, 7))
    expect_type(df$chr, "character")
    expect_type(df$strand, "character")
    expect_type(df$GeneName, "character")
    expect_type(df$GeneStrand, "character")
    expect_type(df$integration_locus, "integer")
    expect_true(typeof(df$Value) %in% c("double", "integer"))
})

test_that(paste(func_name, "reads correctly compressed .gz"), {
    tf <- withr::local_tempfile(fileext = ".tsv.gz")
    readr::write_tsv(sample_df, tf)
    expected_summary_msg <- .summary_ism_import_msg(
        "NEW",
        TRUE,
        c(8, 9),
        "fread"
    )
    expected_summary_msg <- c(
        expected_summary_msg[1],
        paste(
            "*",
            expected_summary_msg[seq(
                from = 2,
                to = length(expected_summary_msg)
            )]
        )
    )
    expected_summary_msg <- paste0(expected_summary_msg, collapse = "\n")
    capture_output(
        suppressMessages(
            expect_message(
                expect_message(
                    withr::with_file(
                        file = tf,
                        code = {
                            df <- import_single_Vispa2Matrix(tf)
                        }
                    )
                ),
                class = "ism_import_summary",
                regexp = expected_summary_msg,
                fixed = TRUE
            )
        )
    )
    expect_equal(dim(df), c(11, 7))
    expect_type(df$chr, "character")
    expect_type(df$strand, "character")
    expect_type(df$GeneName, "character")
    expect_type(df$GeneStrand, "character")
    expect_type(df$integration_locus, "integer")
    expect_true(typeof(df$Value) %in% c("double", "integer"))
})

test_that(paste(func_name, "reads correctly compressed .bz2"), {
    testthat::skip_if_not_installed(pkg = "R.utils")
    tf <- withr::local_tempfile(fileext = ".tsv.bz2")
    readr::write_tsv(sample_df, tf)
    expected_summary_msg <- .summary_ism_import_msg(
        "NEW",
        TRUE,
        c(8, 9),
        "fread"
    )
    expected_summary_msg <- c(
        expected_summary_msg[1],
        paste(
            "*",
            expected_summary_msg[seq(
                from = 2,
                to = length(expected_summary_msg)
            )]
        )
    )
    expected_summary_msg <- paste0(expected_summary_msg, collapse = "\n")
    capture_output(
        suppressMessages(
            expect_message(
                expect_message(
                    withr::with_file(
                        file = tf,
                        code = {
                            df <- import_single_Vispa2Matrix(tf)
                        }
                    )
                ),
                class = "ism_import_summary",
                regexp = expected_summary_msg,
                fixed = TRUE
            )
        )
    )
    expect_equal(dim(df), c(11, 7))
    expect_type(df$chr, "character")
    expect_type(df$strand, "character")
    expect_type(df$GeneName, "character")
    expect_type(df$GeneStrand, "character")
    expect_type(df$integration_locus, "integer")
    expect_true(typeof(df$Value) %in% c("double", "integer"))
})

test_that(paste(func_name, "reads correctly compressed .zip"), {
    tf <- withr::local_tempfile(fileext = ".tsv")
    readr::write_tsv(sample_df, tf)
    tz <- withr::local_tempfile(fileext = ".zip")
    utils::zip(tz, files = tf)
    expected_summary_msg <- .summary_ism_import_msg(
        "NEW",
        TRUE,
        c(8, 9),
        "classic"
    )
    expected_summary_msg <- c(
        expected_summary_msg[1],
        paste(
            "*",
            expected_summary_msg[seq(
                from = 2,
                to = length(expected_summary_msg)
            )]
        )
    )
    expected_summary_msg <- paste0(expected_summary_msg, collapse = "\n")
    capture_output(
        suppressMessages(
            expect_message(
                expect_message(
                    withr::with_file(
                        file = tf,
                        code = {
                            df <- import_single_Vispa2Matrix(tz)
                        }
                    )
                ),
                class = "ism_import_summary",
                regexp = expected_summary_msg,
                fixed = TRUE
            )
        )
    )
    expect_equal(dim(df), c(11, 7))
    expect_type(df$chr, "character")
    expect_type(df$strand, "character")
    expect_type(df$GeneName, "character")
    expect_type(df$GeneStrand, "character")
    expect_type(df$integration_locus, "integer")
    expect_true(typeof(df$Value) %in% c("double", "integer"))
})

test_that(paste(func_name, "keeps additional cols"), {
    tf <- withr::local_tempfile(fileext = ".tsv.gz")
    readr::write_tsv(sample_df_add_columns, tf)
    capture_output(
        suppressMessages(
            expect_message(
                expect_message(
                    withr::with_file(
                        file = tf,
                        code = {
                            df <- import_single_Vispa2Matrix(tf,
                                to_exclude = c("add1", "add2"),
                                keep_excluded = TRUE
                            )
                        }
                    )
                )
            )
        )
    )
    expect_true(all(c(
        "chr", "integration_locus", "strand",
        "GeneName", "GeneStrand", "add1", "add2",
        "CompleteAmplificationID", "Value"
    ) %in% colnames(df)))
})
