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

sample_df_path <- withr::local_tempfile(fileext = "tsv.gz")
readr::write_tsv(x = sample_df, na = "", file = sample_df_path)

sample_df_old <- tidyr::unite(sample_df,
    col = "IS_genomicID",
    dplyr::all_of(c("chr", "integration_locus", "strand")),
    sep = "_", remove = FALSE
)

sample_df_old_path <- withr::local_tempfile(fileext = "tsv.gz")
readr::write_tsv(x = sample_df_old, na = "", file = sample_df_old_path)

sample_df_old_2 <- tidyr::unite(sample_df,
    col = "IS_genomicID",
    dplyr::all_of(c("chr", "integration_locus", "strand")),
    sep = "_", remove = TRUE
)
sample_df_old_2_path <- withr::local_tempfile(fileext = "tsv.gz")
readr::write_tsv(x = sample_df_old_2, na = "", file = sample_df_old_2_path)

sample_df_add_columns <- sample_df |>
    dplyr::mutate(add1 = "a", add2 = 56, .after = "GeneStrand")
sample_df_add_path <- withr::local_tempfile(fileext = "tsv.gz")
readr::write_tsv(x = sample_df_add_columns, na = "", file = sample_df_add_path)

#------------------------------------------------------------------------------#
# .read_with_fread
#------------------------------------------------------------------------------#
test_that(".read_with_fread works with all matrices - standard mand", {
    skip_if_not_installed("data.table")
    df <- .read_with_fread(
        path = sample_df_path, additional_cols = NULL,
        annotated = TRUE, sep = "\t"
    )
    expect_true(nrow(df) == 8 & ncol(df) == 9)
    expect_true(is.character(df$chr) & is.integer(df$integration_locus) &
        is.character(df$strand) & is.character(df$GeneName) &
        is.character(df$GeneStrand) & all(is.integer(c(
        df$id1,
        df$id2,
        df$id3,
        df$id4
    ))))
    df <- .read_with_fread(
        path = sample_df_old_path,
        additional_cols = list(IS_genomicID = "char"),
        annotated = TRUE, sep = "\t"
    )
    expect_true(nrow(df) == 8 & ncol(df) == 10)
    expect_true(is.character(df$chr) & is.integer(df$integration_locus) &
        is.character(df$strand) & is.character(df$GeneName) &
        is.character(df$GeneStrand) & all(is.integer(c(
        df$id1,
        df$id2,
        df$id3,
        df$id4
    ))) &
        is.character(df$IS_genomicID))
    df <- .read_with_fread(
        path = sample_df_add_path,
        additional_cols = list(add1 = "char", add2 = "int"),
        annotated = TRUE, sep = "\t"
    )
    expect_true(nrow(df) == 8 & ncol(df) == 11)
    expect_true(is.character(df$chr) & is.integer(df$integration_locus) &
        is.character(df$strand) & is.character(df$GeneName) &
        is.character(df$GeneStrand) & all(is.integer(c(
        df$id1,
        df$id2,
        df$id3,
        df$id4
    ))) &
        is.character(df$add1) & is.integer(df$add2))
    df <- .read_with_fread(
        path = sample_df_add_path,
        additional_cols = list(add1 = "char", add2 = "_"),
        annotated = TRUE, sep = "\t"
    )
    expect_true(nrow(df) == 8 & ncol(df) == 10)
    expect_true(!c("add2") %in% colnames(df))
})

#------------------------------------------------------------------------------#
# .read_with_readr
#------------------------------------------------------------------------------#
test_that(".read_with_readr works with all matrices - standard mand", {
    withr::local_options(list(ISAnalytics.verbose = FALSE))
    df <- .read_with_readr(
        path = sample_df_path,
        additional_cols = NULL,
        annotated = TRUE,
        sep = "\t"
    )
    expect_true(nrow(df) == 8 & ncol(df) == 9)
    expect_true(is.character(df$chr) & is.integer(df$integration_locus) &
        is.character(df$strand) & is.character(df$GeneName) &
        is.character(df$GeneStrand) & all(is.numeric(c(
        df$id1,
        df$id2,
        df$id3,
        df$id4
    ))))
    df <- .read_with_readr(
        path = sample_df_old_path,
        additional_cols = list("IS_genomicID" = "char"),
        annotated = TRUE, sep = "\t"
    )
    expect_true(nrow(df) == 8 & ncol(df) == 10)
    expect_true(is.character(df$chr) & is.integer(df$integration_locus) &
        is.character(df$strand) & is.character(df$GeneName) &
        is.character(df$GeneStrand) & all(is.numeric(c(
        df$id1,
        df$id2,
        df$id3,
        df$id4
    ))) &
        is.character(df$IS_genomicID))
    df <- .read_with_readr(
        path = sample_df_add_path,
        additional_cols = c("add1" = "char", "add2" = "int"),
        annotated = TRUE, sep = "\t"
    )
    expect_true(nrow(df) == 8 & ncol(df) == 11)
    expect_true(is.character(df$chr) & is.integer(df$integration_locus) &
        is.character(df$strand) & is.character(df$GeneName) &
        is.character(df$GeneStrand) & all(is.numeric(c(
        df$id1,
        df$id2,
        df$id3,
        df$id4
    ))) &
        is.character(df$add1) & is.integer(df$add2))
    df <- .read_with_readr(
        path = sample_df_add_path,
        additional_cols = c("add2" = "_", "add1" = "char"),
        annotated = TRUE, sep = "\t"
    )
    expect_true(nrow(df) == 8 & ncol(df) == 10)
    expect_true(!c("add2") %in% colnames(df))
})

#------------------------------------------------------------------------------#
# .melt_datatable
#------------------------------------------------------------------------------#
test_that(".melt_datatable produces expected output", {
    skip_if_not_installed("data.table")
    as_dt <- data.table::as.data.table(sample_df)
    melted <- .melt_datatable(as_dt,
        id_vars = c(
            mandatory_IS_vars(),
            annotation_IS_vars()
        ),
        sample_col = "Sample",
        value_col = "Value",
        progress = NULL
    )
    expect_true(nrow(melted) == 11 & ncol(melted) == 7)
    expect_true(all(is.character(melted$Sample)))
    expect_true(all(is.numeric(melted$Value)) || all(is.integer(melted$Value)))
    expect_false(any(is.na(melted$Value)))
})

#------------------------------------------------------------------------------#
# .melt_tidyr
#------------------------------------------------------------------------------#
test_that(".melt_tidyr produces expected output", {
    melted <- .melt_tidyr(sample_df,
        id_vars = c(
            mandatory_IS_vars(),
            annotation_IS_vars()
        ),
        sample_col = "Sample",
        value_col = "Value",
        progress = NULL
    )
    expect_true(nrow(melted) == 11 & ncol(melted) == 7)
    expect_true(all(is.character(melted$Sample)))
    expect_true(all(is.numeric(melted$Value)) || all(is.integer(melted$Value)))
    expect_false(any(is.na(melted$Value)))
})

#------------------------------------------------------------------------------#
# .melt_in_parallel
#------------------------------------------------------------------------------#
test_that(".melt_in_parallel produces expected output", {
    melted <- .melt_in_parallel(sample_df,
        n_chunks = 3, melt_fun = .melt_tidyr,
        id_vars = c(
            mandatory_IS_vars(),
            annotation_IS_vars()
        ),
        id_col_name = "Sample", val_col_name = "Value"
    )
    expect_true(nrow(melted) == 11 & ncol(melted) == 7)
    expect_true(all(is.character(melted$Sample)))
    expect_true(all(is.numeric(melted$Value)) || all(is.integer(melted$Value)))
    expect_false(any(is.na(melted$Value)))
})

#------------------------------------------------------------------------------#
# .import_single_matrix
#------------------------------------------------------------------------------#
test_that(".import_single_matrix for EXTERNAL - error msgs", {
    ## Errors on additional columns
    expect_error({
        df <- .import_single_matrix(
            path = sample_df_path, separator = "\t",
            additional_cols = c(1, 2, 3),
            transformations = NULL,
            call_mode = "EXTERNAL",
            id_col_name = "CompleteAmplificationID",
            val_col_name = "Value"
        )
    })
    expect_error({
        df <- .import_single_matrix(
            path = sample_df_path, separator = "\t",
            additional_cols = c("a", "b", "c"),
            transformations = NULL,
            call_mode = "EXTERNAL",
            id_col_name = "CompleteAmplificationID",
            val_col_name = "Value"
        )
    })
    expect_error({
        df <- .import_single_matrix(
            path = sample_df_path, separator = "\t",
            additional_cols = list("a", "b", "c"),
            transformations = NULL,
            call_mode = "EXTERNAL",
            id_col_name = "CompleteAmplificationID",
            val_col_name = "Value"
        )
    })
    expect_error(
        {
            df <- .import_single_matrix(
                path = sample_df_path, separator = "\t",
                additional_cols = list(
                    a = "c", b = "char",
                    c = "i"
                ),
                transformations = NULL,
                call_mode = "EXTERNAL",
                id_col_name = "CompleteAmplificationID",
                val_col_name = "Value"
            )
        },
        class = "add_types_err"
    )
    ## Missing mandatory vars
    expect_error(
        {
            df <- .import_single_matrix(
                path = sample_df_old_2_path, separator = "\t",
                additional_cols = NULL,
                transformations = NULL,
                call_mode = "EXTERNAL",
                id_col_name = "CompleteAmplificationID",
                val_col_name = "Value"
            )
        },
        class = "im_single_miss_mand_vars"
    )
})

test_that(".import_single_matrix for EXTERNAL - works", {
    ## Standard df
    expect_message(
        {
            expect_message({
                expect_message({
                    df <- .import_single_matrix(
                        path = sample_df_path, separator = "\t",
                        additional_cols = NULL,
                        transformations = NULL,
                        call_mode = "EXTERNAL",
                        id_col_name = "CompleteAmplificationID",
                        val_col_name = "Value"
                    )
                })
            })
        },
        class = "ism_import_summary"
    )
    ## With additional cols
    expect_message(
        {
            expect_message({
                expect_message({
                    df <- .import_single_matrix(
                        path = sample_df_add_path,
                        separator = "\t",
                        additional_cols = list(
                            add1 = "char",
                            add2 = "int"
                        ),
                        transformations = NULL,
                        call_mode = "EXTERNAL",
                        id_col_name = "CompleteAmplificationID",
                        val_col_name = "Value"
                    )
                })
            })
        },
        class = "ism_import_summary"
    )
    ## With additional cols and transformations
    expect_message(
        {
            expect_message({
                expect_message({
                    expect_message({
                        df <- .import_single_matrix(
                            path = sample_df_add_path,
                            separator = "\t",
                            additional_cols = list(
                                add1 = "char",
                                add2 = "int"
                            ),
                            transformations = list(
                                add1 = ~ stringr::str_to_upper(.x)
                            ),
                            call_mode = "EXTERNAL",
                            id_col_name = "CompleteAmplificationID",
                            val_col_name = "Value"
                        )
                    })
                })
            })
        },
        class = "ism_import_summary"
    )
    expect_true(all(df$add1 == "A"))
})

#------------------------------------------------------------------------------#
# import_single_Vispa2Matrix
#------------------------------------------------------------------------------#
test_that("import_single_Vispa2Matrix warns deprecation", {
    expect_defunct({
        df <- import_single_Vispa2Matrix(sample_df_add_path,
            to_exclude = "chr"
        )
    })
    expect_defunct({
        df <- import_single_Vispa2Matrix(sample_df_add_path,
            keep_excluded = TRUE
        )
    })
})

### ---- Testing with custom vars -------------------------------------------###
test_that("import_single_Vispa2Matrix reads correctly with custom vars", {
    sample_df_mod <- sample_df |>
        dplyr::mutate(gap = 100, junction = 100, .after = "GeneStrand")
    tmp_mand_vars <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "chr", "char", NULL, "required", "chromosome",
        "integration_locus", "int", NULL, "required", "locus",
        "strand", "char", NULL, "required", "is_strand",
        "gap", "int", NULL, "required", NA_character_,
        "junction", "int", NULL, "required", NA_character_
    )
    withr::with_options(
        list(ISAnalytics.mandatory_is_vars = tmp_mand_vars),
        {
            tf <- withr::local_tempfile(fileext = ".tsv.gz")
            readr::write_tsv(sample_df_mod, tf)
            expected_summary_msg <- .summary_ism_import_msg(
                TRUE,
                c(8, 11),
                "fread",
                4
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
            expected_summary_msg <- paste0(expected_summary_msg,
                collapse = "\n"
            )
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
            expect_equal(dim(df), c(11, 9))
            expect_type(df$chr, "character")
            expect_type(df$strand, "character")
            expect_type(df$GeneName, "character")
            expect_type(df$GeneStrand, "character")
            expect_type(df$integration_locus, "integer")
            expect_type(df$gap, "integer")
            expect_type(df$junction, "integer")
            expect_true(typeof(df$Value) %in% c("double", "integer"))
        }
    )
})

test_that(paste(
    "import_single_Vispa2Matrix",
    "reads correctly with custom vars and dates"
), {
    sample_df_mod <- sample_df |>
        dplyr::mutate(custom_date = c(
            "2016-11-15", "2017-06-23",
            "2017-06-23", "2017-06-23",
            "2018-03-15", "2018-03-15",
            "2018-04-14", "2018-04-14"
        ))
    tmp_mand_vars <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "chr", "char", NULL, "required", "chromosome",
        "integration_locus", "int", NULL, "required", "locus",
        "strand", "char", NULL, "required", "is_strand",
        "custom_date", "ymd", NULL, "required", NA_character_
    )
    withr::with_options(
        list(ISAnalytics.mandatory_is_vars = tmp_mand_vars),
        {
            tf <- withr::local_tempfile(fileext = ".tsv.gz")
            readr::write_tsv(sample_df_mod, tf)
            expected_summary_msg <- .summary_ism_import_msg(
                TRUE,
                c(8, 10),
                "fread",
                4
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
            expected_summary_msg <- paste0(expected_summary_msg,
                collapse = "\n"
            )
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
            expect_equal(dim(df), c(11, 8))
            expect_type(df$chr, "character")
            expect_type(df$strand, "character")
            expect_type(df$GeneName, "character")
            expect_type(df$GeneStrand, "character")
            expect_type(df$integration_locus, "integer")
            expect_true(typeof(df$Value) %in% c("double", "integer"))
            expect_true(all(lubridate::is.Date(df$custom_date)))
        }
    )
})
