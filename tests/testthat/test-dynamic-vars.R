#------------------------------------------------------------------------------#
# Tests for general use of dynamic variables
# * Internals
# * Getters/setters
#------------------------------------------------------------------------------#

### ---- Internals ----------------------------------------------------------###
test_that(".new_IS_vars_checks errors if input not in format", {
    specs <- 5
    expect_error({
        vars <- .new_IS_vars_checks(specs, "err", NULL)
    })
    specs <- tibble::tribble(
        ~a, ~b,
        "var1", "int"
    )
    expect_error({
        vars <- .new_IS_vars_checks(specs, "err", NULL)
    })
    specs <- tibble::tribble(
        ~names, ~types,
        "var1", integer()
    )
    expect_error({
        vars <- .new_IS_vars_checks(specs, "err", NULL)
    })
    specs <- tibble::tribble(
        ~names, ~types, ~transform, ~flag,
        "var1", "unknown", NULL, "no"
    )
    expect_error({
        vars <- .new_IS_vars_checks(specs, "err", NULL)
    })
    specs <- c(var1 = "unknown")
    expect_error({
        vars <- .new_IS_vars_checks(specs, "err", NULL)
    })
})

test_that(".new_IS_vars_checks accepts named vec", {
    specs <- c(var1 = "char", var2 = "int", var3 = "numeric")
    vars <- .new_IS_vars_checks(specs, "err", NULL)
    expected <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "var1", "char", NULL, "required", NA_character_,
        "var2", "int", NULL, "required", NA_character_,
        "var3", "numeric", NULL, "required", NA_character_
    )
    expect_equal(vars, expected)
})

test_that(".new_IS_vars_checks accepts named list", {
    specs <- list(var1 = "char", var2 = "int", var3 = "numeric")
    vars <- .new_IS_vars_checks(specs, "err", NULL)
    expected <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "var1", "char", NULL, "required", NA_character_,
        "var2", "int", NULL, "required", NA_character_,
        "var3", "numeric", NULL, "required", NA_character_
    )
    expect_equal(vars, expected)
})

test_that(".new_IS_vars_checks accepts data frame", {
    specs <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "var1", "char", NULL, "required", NA_character_,
        "var2", "int", NULL, "required", NA_character_,
        "var3", "numeric", NULL, "required", NA_character_
    )
    vars <- .new_IS_vars_checks(specs, "err", NULL)
    expect_equal(vars, specs)
})

test_that(".apply_col_transform works with purrr lambdas", {
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
    specs <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "chr", "char", ~ paste0("chr", .x), "required", NA_character_
    )
    expected <- sample_df |>
        dplyr::mutate(chr = paste0("chr", .data$chr))
    res <- .apply_col_transform(sample_df, specs)
    expect_equal(res, expected)
})

test_that(".apply_col_transform works with functions single arg", {
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
    foo <- function(x) {
        paste0("chr", x)
    }
    specs <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "chr", "char", foo, "required", NA_character_
    )
    expected <- sample_df |>
        dplyr::mutate(chr = foo(.data$chr))
    res <- .apply_col_transform(sample_df, specs)
    expect_equal(res, expected)
})

test_that(".eval_tag_duplicates ok with all err politic", {
    tmp_specs <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "var1", "char", NULL, "required", "project_id",
        "var2", "char", NULL, "required", "project_id",
        "var3", "int", NULL, "required", "pcr_replicate",
        "var4", "char", NULL, "required", "tp_days",
        "var5", "int", NULL, "required", "tp_days"
    )
    expect_error(
        {
            res <- .eval_tag_duplicates(tmp_specs, duplicate_politic = "error")
        },
        class = "tag_dupl_err"
    )
})

test_that(".eval_tag_duplicates ok with all first politic", {
    tmp_specs <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "var1", "char", NULL, "required", "project_id",
        "var2", "char", NULL, "required", "project_id",
        "var3", "int", NULL, "required", "pcr_replicate",
        "var4", "char", NULL, "required", "tp_days",
        "var5", "int", NULL, "required", "tp_days"
    )
    expected_out <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "var3", "int", NULL, "required", "pcr_replicate",
        "var1", "char", NULL, "required", "project_id",
        "var4", "char", NULL, "required", "tp_days"
    )
    expect_message(
        {
            res <- .eval_tag_duplicates(tmp_specs, duplicate_politic = "first")
        },
        class = "tag_dupl_warn"
    )
    expect_equal(res, expected_out)
})

test_that(".eval_tag_duplicates ok with all keep politic", {
    tmp_specs <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "var1", "char", NULL, "required", "project_id",
        "var2", "char", NULL, "required", "project_id",
        "var3", "int", NULL, "required", "pcr_replicate",
        "var4", "char", NULL, "required", "tp_days",
        "var5", "int", NULL, "required", "tp_days"
    )
    res <- .eval_tag_duplicates(tmp_specs, duplicate_politic = "keep")

    expect_equal(res, tmp_specs)
})

test_that(".eval_tag_duplicates ok with custom tag politic", {
    tmp_specs <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "var1", "char", NULL, "required", "project_id",
        "var2", "char", NULL, "required", "project_id",
        "var3", "int", NULL, "required", "pcr_replicate",
        "var4", "char", NULL, "required", "tp_days",
        "var5", "int", NULL, "required", "tp_days"
    )
    politics <- c(
        "project_id" = "error", "pcr_replicate" = "error",
        "tp_days" = "keep"
    )
    expect_error(
        {
            res <- .eval_tag_duplicates(tmp_specs, duplicate_politic = politics)
        },
        class = "tag_dupl_err"
    )

    politics <- c(
        "project_id" = "first", "pcr_replicate" = "error",
        "tp_days" = "keep"
    )

    expected_out <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "var3", "int", NULL, "required", "pcr_replicate",
        "var1", "char", NULL, "required", "project_id",
        "var4", "char", NULL, "required", "tp_days",
        "var5", "int", NULL, "required", "tp_days"
    )

    expect_message(
        {
            res <- .eval_tag_duplicates(tmp_specs, duplicate_politic = politics)
        },
        class = "tag_dupl_warn"
    )
    expect_equal(res, expected_out)
})

test_that(".eval_tag_types works as expected", {
    tmp_specs <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "var1", "char", NULL, "required", "project_id",
        "var2", "char", NULL, "required", "project_id",
        "var3", "int", NULL, "required", "pcr_replicate",
        "var4", "char", NULL, "required", "tp_days",
        "var5", "int", NULL, "required", "tp_days"
    )
    types_specs <- c(
        "project_id" = "char", "pcr_replicate" = "int",
        "tp_days" = "char"
    )
    expected_out <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "var3", "int", NULL, "required", "pcr_replicate",
        "var1", "char", NULL, "required", "project_id",
        "var2", "char", NULL, "required", "project_id",
        "var4", "char", NULL, "required", "tp_days",
    )
    result <- .eval_tag_types(df = tmp_specs, types = types_specs)
    expect_equal(result, expected_out)

    types_specs <- c(
        "project_id" = "int", "pcr_replicate" = "int",
        "tp_days" = "char"
    )
    expect_error(
        {
            result <- .eval_tag_types(df = tmp_specs, types = types_specs)
        },
        class = "tag_type_err"
    )
})

test_that(".check_required_cols errors for missing tags", {
    tmp_specs <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "var1", "char", NULL, "required", "project_id",
        "var2", "char", NULL, "required", "project_id",
        "var3", "int", NULL, "required", "pcr_replicate",
        "var4", "char", NULL, "required", "tp_days",
        "var5", "int", NULL, "required", "tp_days"
    )

    req_tags <- c("project_id", "pcr_replicate", "tp_days", "subject")
    expect_error(
        {
            result <- .check_required_cols(
                required_tags = req_tags,
                vars_df = tmp_specs,
                duplicate_politic = "error"
            )
        },
        class = "missing_tags_err"
    )

    req_tags <- c("subject")
    expect_error(
        {
            result <- .check_required_cols(
                required_tags = req_tags,
                vars_df = tmp_specs,
                duplicate_politic = "error"
            )
        },
        class = "missing_tags_err"
    )

    req_tags <- c(
        "project_id" = "char",
        "pcr_replicate" = "int",
        "tp_days" = "char",
        "subject" = "char"
    )
    expect_error(
        {
            result <- .check_required_cols(
                required_tags = req_tags,
                vars_df = tmp_specs,
                duplicate_politic = "error"
            )
        },
        class = "missing_tags_err"
    )

    req_tags <- c("subject" = "char")
    expect_error(
        {
            result <- .check_required_cols(
                required_tags = req_tags,
                vars_df = tmp_specs,
                duplicate_politic = "error"
            )
        },
        class = "missing_tags_err"
    )
})

test_that(".check_required_cols errors for duplicates", {
    tmp_specs <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "var1", "char", NULL, "required", "project_id",
        "var2", "char", NULL, "required", "project_id",
        "var3", "int", NULL, "required", "pcr_replicate",
        "var4", "char", NULL, "required", "tp_days",
        "var5", "int", NULL, "required", "tp_days"
    )

    req_tags <- c("project_id", "pcr_replicate", "tp_days")
    expect_error(
        {
            result <- .check_required_cols(
                required_tags = req_tags,
                vars_df = tmp_specs,
                duplicate_politic = "error"
            )
        },
        class = "tag_dupl_err"
    )

    custom_politic <- c(
        "project_id" = "error",
        "pcr_replicate" = "first",
        "tp_days" = "keep"
    )
    expect_error(
        {
            result <- .check_required_cols(
                required_tags = req_tags,
                vars_df = tmp_specs,
                duplicate_politic = custom_politic
            )
        },
        class = "tag_dupl_err"
    )

    req_tags <- c(
        "project_id" = "char", "pcr_replicate" = "int",
        "tp_days" = "char"
    )
    expect_error(
        {
            result <- .check_required_cols(
                required_tags = req_tags,
                vars_df = tmp_specs,
                duplicate_politic = custom_politic
            )
        },
        class = "tag_dupl_err"
    )
})

test_that(".check_required_cols errors for types", {
    tmp_specs <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "var1", "char", NULL, "required", "project_id",
        "var2", "char", NULL, "required", "project_id",
        "var3", "int", NULL, "required", "pcr_replicate",
        "var4", "char", NULL, "required", "tp_days",
        "var5", "int", NULL, "required", "tp_days"
    )

    req_tags <- c(
        "project_id" = "char", "pcr_replicate" = "char",
        "tp_days" = "char"
    )
    expect_error(
        {
            result <- .check_required_cols(
                required_tags = req_tags,
                vars_df = tmp_specs,
                duplicate_politic = "keep"
            )
        },
        class = "tag_type_err"
    )
})

test_that(".check_required_cols works as expected", {
    tmp_specs <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "var1", "char", NULL, "required", "project_id",
        "var2", "char", NULL, "required", "project_id",
        "var3", "int", NULL, "required", "pcr_replicate",
        "var4", "char", NULL, "required", "tp_days",
        "var5", "int", NULL, "required", "tp_days"
    )
    req_tags <- list(
        "project_id" = "char", "pcr_replicate" = "int",
        "tp_days" = c("char", "int")
    )
    result <- .check_required_cols(
        required_tags = req_tags,
        vars_df = tmp_specs,
        duplicate_politic = "keep"
    )
    expected_out <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "var3", "int", NULL, "required", "pcr_replicate",
        "var1", "char", NULL, "required", "project_id",
        "var2", "char", NULL, "required", "project_id",
        "var4", "char", NULL, "required", "tp_days",
        "var5", "int", NULL, "required", "tp_days"
    )
    expect_equal(result, expected_out)
    req_tags <- list(
        "project_id" = "char", "pcr_replicate" = "int",
        "tp_days" = c("char", "logi")
    )
    expected_out <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "var3", "int", NULL, "required", "pcr_replicate",
        "var1", "char", NULL, "required", "project_id",
        "var2", "char", NULL, "required", "project_id",
        "var4", "char", NULL, "required", "tp_days",
    )
    result <- .check_required_cols(
        required_tags = req_tags,
        vars_df = tmp_specs,
        duplicate_politic = "keep"
    )
    expect_equal(result, expected_out)
})

### ---- Getters/setters ----------------------------------------------------###
test_that("mandatory_IS_vars works as expected", {
    withr::with_options(list(ISAnalytics.mandatory_is_vars = "default"), {
        vars <- mandatory_IS_vars()
        expect_equal(vars, .default_mandatory_IS_vars()$names)
        vars <- mandatory_IS_vars(TRUE)
        expect_equal(vars, .default_mandatory_IS_vars())
    })
    temp_specs <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "var1", "char", NULL, "required", NA_character_,
        "var2", "int", NULL, "required", NA_character_,
        "var3", "numeric", NULL, "required", NA_character_
    )
    withr::with_options(list(ISAnalytics.mandatory_is_vars = temp_specs), {
        vars <- mandatory_IS_vars()
        expect_equal(vars, temp_specs$names)
        vars <- mandatory_IS_vars(TRUE)
        expect_equal(vars, temp_specs)
    })
})

test_that("annotation_IS_vars works as expected", {
    withr::with_options(list(ISAnalytics.genomic_annotation_vars = "default"), {
        vars <- annotation_IS_vars()
        expect_equal(vars, .default_annotation_IS_vars()$names)
        vars <- annotation_IS_vars(TRUE)
        expect_equal(vars, .default_annotation_IS_vars())
    })
    temp_specs <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "var1", "char", NULL, "required", NA_character_,
        "var2", "int", NULL, "required", NA_character_,
        "var3", "numeric", NULL, "required", NA_character_
    )
    withr::with_options(list(ISAnalytics.genomic_annotation_vars = temp_specs), {
        vars <- annotation_IS_vars()
        expect_equal(vars, temp_specs$names)
        vars <- annotation_IS_vars(TRUE)
        expect_equal(vars, temp_specs)
    })
})

test_that("association_file_columns works as expected", {
    withr::with_options(list(ISAnalytics.af_specs = "default"), {
        vars <- association_file_columns()
        expect_equal(vars, .default_af_cols()$names)
        vars <- association_file_columns(TRUE)
        expect_equal(vars, .default_af_cols(), ignore_attr = TRUE)
    })
    temp_specs <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "var1", "char", NULL, "required", NA_character_,
        "var2", "int", NULL, "required", NA_character_,
        "var3", "numeric", NULL, "required", NA_character_
    )
    withr::with_options(list(ISAnalytics.af_specs = temp_specs), {
        vars <- association_file_columns()
        expect_equal(vars, temp_specs$names)
        vars <- association_file_columns(TRUE)
        expect_equal(vars, temp_specs)
    })
})

test_that("set_mandatory_IS_vars signals missing tags", {
    temp_specs <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "var1", "char", NULL, "required", NA_character_,
        "var2", "int", NULL, "required", NA_character_,
        "var3", "numeric", NULL, "required", NA_character_
    )
    before <- getOption("ISAnalytics.mandatory_is_vars")
    expect_message(
        expect_warning(
            {
                set_mandatory_IS_vars(temp_specs)
            },
            class = "missing_crit_tags"
        )
    )
    after <- getOption("ISAnalytics.mandatory_is_vars")
    expect_equal(after, temp_specs)
    options(ISAnalytics.mandatory_is_vars = before)
})

test_that("set_mandatory_IS_vars works as expected", {
    temp_specs <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "var1", "char", NULL, "required", "chromosome",
        "var2", "int", NULL, "required", "locus",
        "var3", "numeric", NULL, "required", "is_strand"
    )
    before <- getOption("ISAnalytics.mandatory_is_vars")
    expect_message(
        set_mandatory_IS_vars(temp_specs)
    )
    after <- getOption("ISAnalytics.mandatory_is_vars")
    expect_equal(after, temp_specs)
    options(ISAnalytics.mandatory_is_vars = before)
})

test_that("set_annotation_IS_vars signals missing tags", {
    temp_specs <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "var1", "char", NULL, "required", NA_character_,
        "var2", "int", NULL, "required", NA_character_,
        "var3", "numeric", NULL, "required", NA_character_
    )
    before <- getOption("ISAnalytics.genomic_annotation_vars")
    expect_message(
        expect_message(
            {
                set_annotation_IS_vars(temp_specs)
            },
            class = "missing_crit_tags"
        )
    )
    after <- getOption("ISAnalytics.genomic_annotation_vars")
    expect_equal(after, temp_specs)
    options(ISAnalytics.genomic_annotation_vars = before)
})

test_that("set_annotation_IS_vars works as expected", {
    temp_specs <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "var1", "char", NULL, "required", "gene_symbol",
        "var2", "int", NULL, "required", "gene_strand",
        "var3", "numeric", NULL, "required", NA_character_
    )
    before <- getOption("ISAnalytics.genomic_annotation_vars")
    expect_message(
        set_annotation_IS_vars(temp_specs)
    )
    after <- getOption("ISAnalytics.genomic_annotation_vars")
    expect_equal(after, temp_specs)
    options(ISAnalytics.genomic_annotation_vars = before)
})

test_that("set_af_columns_def works as expected", {
    temp_specs <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "var1", "char", NULL, "required", NA_character_,
        "var2", "int", NULL, "required", NA_character_,
        "var3", "numeric", NULL, "required", NA_character_
    )
    before <- getOption("ISAnalytics.af_specs")
    expect_message(
        expect_message(
            {
                set_af_columns_def(temp_specs)
            },
            class = "missing_crit_tags"
        )
    )
    after <- getOption("ISAnalytics.af_specs")
    expect_equal(after, temp_specs)
    options(ISAnalytics.af_specs = before)
})

test_that("iss_stats_specs works as expected", {
    withr::with_options(list(ISAnalytics.iss_stats_specs = "default"), {
        vars <- iss_stats_specs()
        expect_equal(vars, .default_iss_stats_specs()$names)
        vars <- iss_stats_specs(TRUE)
        expect_equal(vars, .default_iss_stats_specs())
    })
    temp_specs <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "var1", "char", NULL, "required", "vispa_concatenate",
        "var2", "int", NULL, "required", "tag_seq",
        "var3", "numeric", NULL, "required", NA_character_
    )
    withr::with_options(list(ISAnalytics.iss_stats_specs = temp_specs), {
        vars <- iss_stats_specs()
        expect_equal(vars, temp_specs$names)
        vars <- iss_stats_specs(TRUE)
        expect_equal(vars, temp_specs)
    })
})

test_that("set_iss_stats_specs works as expected", {
    temp_specs <- tibble::tribble(
        ~names, ~types, ~transform, ~flag, ~tag,
        "var1", "char", NULL, "required", "vispa_concatenate",
        "var2", "int", NULL, "required", "tag_seq",
        "var3", "numeric", NULL, "required", NA_character_
    )
    before <- getOption("ISAnalytics.iss_stats_specs")
    expect_message(
        set_iss_stats_specs(temp_specs)
    )
    after <- getOption("ISAnalytics.iss_stats_specs")
    expect_equal(after, temp_specs)
    options(ISAnalytics.iss_stats_specs = before)
})

test_that("matrix_file_suffixes returns expected", {
    expected <- tibble::tibble(
        quantification = c(
            "seqCount", "fragmentEstimate", "barcodeCount",
            "cellCount", "ShsCount", "seqCount",
            "fragmentEstimate", "barcodeCount",
            "cellCount", "ShsCount"
        ),
        matrix_type = c(
            "annotated", "annotated", "annotated",
            "annotated", "annotated", "not_annotated",
            "not_annotated", "not_annotated", "not_annotated",
            "not_annotated"
        ),
        file_suffix = c(
            "seqCount_matrix.no0.annotated.tsv.gz",
            "fragmentEstimate_matrix.no0.annotated.tsv.gz",
            "barcodeCount_matrix.no0.annotated.tsv.gz",
            "cellCount_matrix.no0.annotated.tsv.gz",
            "ShsCount_matrix.no0.annotated.tsv.gz",
            "seqCount_matrix.tsv.gz", "fragmentEstimate_matrix.tsv.gz",
            "barcodeCount_matrix.tsv.gz", "cellCount_matrix.tsv.gz",
            "ShsCount_matrix.tsv.gz"
        )
    )
    matrix_suffixes <- matrix_file_suffixes() |>
        dplyr::arrange(.data$quantification)
    expect_equal(matrix_suffixes, expected |>
        dplyr::arrange(.data$quantification))
})

test_that("set_matrix_file_suffixes works as expected", {
    before <- getOption("ISAnalytics.matrix_file_suffix")
    expect_error(
        {
            set_matrix_file_suffixes(
                annotation_suffix = list(annotated = "ann")
            )
        },
        class = "miss_annot_suff_specs"
    )
    expect_message({
        expect_message(
            {
                set_matrix_file_suffixes(quantification_suffix = list(
                    fragmentEstimate = "fe", seqCount = "sc"
                ))
            },
            class = "missing_quant_specs"
        )
    })
    expected_out <- tibble::tibble(
        quantification = c(
            "fragmentEstimate", "seqCount", "fragmentEstimate",
            "seqCount"
        ),
        matrix_type = c(
            "annotated", "annotated",
            "not_annotated", "not_annotated"
        ),
        file_suffix = c(
            "fe_matrix.no0.annotated.tsv.gz",
            "sc_matrix.no0.annotated.tsv.gz",
            "fe_matrix.tsv.gz", "sc_matrix.tsv.gz"
        )
    ) |>
        dplyr::arrange(.data$quantification)
    expect_equal(getOption("ISAnalytics.matrix_file_suffix"), expected_out)
    options(ISAnalytics.matrix_file_suffix = before)
})
