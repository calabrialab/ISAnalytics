library(ISAnalytics)
func_name <- c(
    ".aggregate_meta", "aggregate_metadata",
    "aggregate_values_by_key"
)

#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
# Path to example association file
path_af <- system.file("extdata", "ex_association_file.tsv",
    package = "ISAnalytics"
)

# Path to correct file system example
path_root_correct <- system.file("extdata", "fs.zip",
    package = "ISAnalytics"
)
root_correct <- unzip_file_system(path_root_correct, "fs")

# Association file
association_file <- withr::with_options(
    list(ISAnalytics.widgets = FALSE, ISAnalytics.verbose = FALSE),
    {
        import_association_file(path_af, root_correct,
            dates_format = "dmy", filter_for = list(ProjectID = "CLOEXP"),
            import_iss = TRUE
        )
    }
)

# Matrices
matrices <- withr::with_options(
    list(ISAnalytics.widgets = FALSE, ISAnalytics.verbose = FALSE),
    {
        import_parallel_Vispa2Matrices_auto(
            association_file = association_file,
            quantification_type = c("seqCount", "fragmentEstimate"),
            matrix_type = "annotated", workers = 2, patterns = NULL,
            matching_opt = "ANY",
            dates_format = "dmy", multi_quant_matrix = FALSE
        )
    }
)

# Standard grouping keys
key <- c(
    "SubjectID",
    "CellMarker",
    "Tissue",
    "TimePoint"
)

#------------------------------------------------------------------------------#
# Tests .aggregate_meta
#------------------------------------------------------------------------------#
test_that(paste(func_name[1], "aggregates correct default"), {
    aggreg_meta <- .aggregate_meta(association_file,
        grouping_keys = key,
        function_tbl = default_meta_agg()
    )
    expect_true(aggreg_meta %>% nrow() == 4)
})

test_that(paste(func_name[1], "ignores column if not present"), {
    func_table <- default_meta_agg() %>%
        tibble::add_row(
            Column = "A", Function = list(~ sum(.x)), Args = NA,
            Output_colname = "{.col}_sum"
        )
    aggreg_meta <- .aggregate_meta(association_file,
        grouping_keys = key,
        function_tbl = func_table
    )
    expect_true(!"A_sum" %in% colnames(aggreg_meta))
})

test_that(paste(func_name[1], "works for mixed function formula"), {
    func_table <- tibble::tribble(
        ~Column, ~Function, ~Args, ~Output_colname,
        "FusionPrimerPCRDate", ~ min(.x, na.rm = TRUE), NA, "{.col}_min",
        "FusionPrimerPCRDate", min, list(na.rm = TRUE), "{.col}_min_fun"
    )
    aggreg_meta <- .aggregate_meta(association_file,
        grouping_keys = key,
        function_tbl = func_table
    )
    expect_true(all(aggreg_meta$FusionPrimerPCRDate_min ==
        aggreg_meta$FusionPrimerPCRDate_min_fun))
})

#------------------------------------------------------------------------------#
# Tests aggregate_metadata
#------------------------------------------------------------------------------#
# Test input
test_that(paste(func_name[2], "stops if grouping keys is null"), {
    expect_error({
        agg_meta <- aggregate_metadata(
            association_file = association_file,
            grouping_keys = NULL
        )
    })
})
test_that(paste(func_name[2], "stops if grouping keys is not a char vector"), {
    expect_error({
        agg_meta <- aggregate_metadata(
            association_file = association_file,
            grouping_keys = 1
        )
    })
    expect_error({
        agg_meta <- aggregate_metadata(
            association_file = association_file,
            grouping_keys = c(1, 2)
        )
    })
})
test_that(paste(func_name[2], "stops if grouping keys are missing from af"), {
    expect_error({
        agg_meta <- aggregate_metadata(
            association_file = association_file,
            grouping_keys = "a"
        )
    })
    expect_error({
        agg_meta <- aggregate_metadata(
            association_file = association_file,
            grouping_keys = c("a", "ProjectID")
        )
    })
})

test_that(paste(func_name[2], "warns deprecation of import_stats"), {
    expect_deprecated({
        agg_meta <- aggregate_metadata(
            association_file = association_file,
            import_stats = TRUE
        )
    })
})

test_that(paste(func_name[2], "warns if no columns found"), {
    func_table <- tibble::tribble(
        ~Column, ~Function, ~Args, ~Output_colname,
        "A", ~ min(.x, na.rm = TRUE), NA, "{.col}_min",
        "B", min, list(na.rm = TRUE), "{.col}_min_fun"
    )
    expect_message({
        aggreg_meta <- aggregate_metadata(
            association_file = association_file,
            grouping_keys = key,
            aggregating_functions = func_table
        )
    })
    expect_null(aggreg_meta)
})

#------------------------------------------------------------------------------#
# Tests aggregate_values_by_key
#------------------------------------------------------------------------------#
# Test input
test_that(paste(func_name[3], "stops if x incorrect"), {
    expect_error({
        agg <- aggregate_values_by_key(
            x = 1,
            association_file = association_file
        )
    })
    expect_error(
        {
            agg <- aggregate_values_by_key(
                x = tibble::tibble(a = c(1, 2, 3)),
                association_file = association_file
            )
        },
        regexp = .non_ISM_error()
    )
    expect_error({
        agg <- aggregate_values_by_key(matrices$seqCount, association_file,
            value_cols = "x"
        )
    })
    expect_error({
        agg <- aggregate_values_by_key(matrices$seqCount, association_file,
            value_cols = c("chr")
        )
    })
    expect_error({
        agg <- aggregate_values_by_key(
            matrices$seqCount %>%
                dplyr::select(-c("CompleteAmplificationID")),
            association_file
        )
    })
    expect_error({
        agg <- aggregate_values_by_key(matrices, association_file,
            value_cols = "x"
        )
    })
    expect_error({
        agg <- aggregate_values_by_key(matrices, association_file,
            value_cols = c("chr")
        )
    })
    mod_list <- purrr::map(matrices, function(m) {
        m %>% dplyr::select(-c("CompleteAmplificationID"))
    })
    expect_error({
        agg <- aggregate_values_by_key(
            mod_list,
            association_file
        )
    })
})

test_that(paste(func_name[3], "stops if key is incorrect"), {
    expect_error({
        agg <- aggregate_values_by_key(
            x = matrices$seqCount,
            association_file = association_file,
            key = 1
        )
    })
    expect_error(
        {
            agg <- aggregate_values_by_key(
                x = matrices$seqCount,
                association_file = association_file,
                key = "x"
            )
        },
        regexp = "Key fields are missing from association file"
    )
    expect_error(
        {
            agg <- aggregate_values_by_key(
                x = matrices$seqCount,
                association_file = association_file,
                key = c("SubjectID", "x")
            )
        },
        regexp = "Key fields are missing from association file"
    )
})

test_that(paste(func_name[3], "stops if lambda is incorrect"), {
    expect_error({
        agg <- aggregate_values_by_key(
            x = isadf,
            association_file = association_file,
            key = "SubjectID",
            lambda = 1
        )
    })
    expect_error({
        agg <- aggregate_values_by_key(
            x = isadf,
            association_file = association_file,
            key = "SubjectID",
            lambda = sum
        )
    })
    expect_error({
        agg <- aggregate_values_by_key(
            x = isadf,
            association_file = association_file,
            key = "SubjectID",
            lambda = c("sum", "mean")
        )
    })
})

test_that(paste(func_name[3], "stops if group is incorrect"), {
    expect_error({
        agg <- aggregate_values_by_key(
            x = matrices$seqCount,
            association_file = association_file,
            group = c(1, 2)
        )
    })
    expect_error(
        {
            agg <- aggregate_values_by_key(
                x = matrices$seqCount,
                association_file = association_file,
                group = c("xfsf", "fwre")
            )
        },
        regexp = "Grouping variables not found"
    )
})

# Test Values
test_that(paste(func_name[3], "succeeds for single IS matrix"), {
    expect_error(
        {
            agg <- aggregate_values_by_key(matrices$seqCount, association_file)
        },
        regexp = NA
    )
    expect_true(all(c(
        "chr", "integration_locus", "strand", "GeneName",
        "GeneStrand", "SubjectID",
        "Value_sum"
    ) %in% colnames(agg)))
    comp <- comparison_matrix(matrices)
    agg <- aggregate_values_by_key(comp, association_file,
        value_cols = c(
            "fragmentEstimate",
            "seqCount"
        )
    )
    expect_true(all(c(
        "chr", "integration_locus", "strand", "GeneName",
        "GeneStrand", "SubjectID",
        "fragmentEstimate_sum", "seqCount_sum"
    ) %in% colnames(agg)))
})

test_that(paste(func_name[3], "succeeds for list of IS matrices"), {
    expect_error(
        {
            agg <- aggregate_values_by_key(matrices, association_file)
        },
        regexp = NA
    )
    expect_true(all(c(
        "chr", "integration_locus", "strand", "GeneName",
        "GeneStrand", "SubjectID",
        "Value_sum"
    ) %in% colnames(agg[[1]])))
    expect_true(all(c(
        "chr", "integration_locus", "strand", "GeneName",
        "GeneStrand", "SubjectID",
        "Value_sum"
    ) %in% colnames(agg[[2]])))
})
