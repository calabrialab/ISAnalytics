#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
min_example <- tibble::tribble(
    ~ProjectID, ~SubjectID, ~Tissue, ~FusionPrimerPCRDate, ~Value,
    "PJ01", "S1", "BM", lubridate::ymd("2020-02-23"), 34646,
    "PJ01", "S1", "BM", lubridate::ymd("2020-02-23"), 56762,
    "PJ01", "S1", "BM", lubridate::ymd("2020-02-23"), 5678,
    "PJ01", "S1", "PB", lubridate::ymd("2020-02-23"), 7934,
    "PJ01", "S1", "PB", lubridate::ymd("2020-02-23"), 786,
    "PJ01", "S1", "PB", lubridate::ymd("2020-03-09"), 5322,
    "PJ01", "S2", "BM", lubridate::ymd("2020-04-05"), 8933,
    "PJ01", "S2", "BM", lubridate::ymd("2020-01-05"), 894,
    "PJ01", "S2", "PB", lubridate::ymd("2020-01-05"), 1356,
    "PJ01", "S2", "PB", lubridate::ymd("2020-01-05"), 784,
)
agg_key_meta <- c("SubjectID", "Tissue")

#------------------------------------------------------------------------------#
# Tests .aggregate_meta
#------------------------------------------------------------------------------#
test_that(".aggregate_meta aggregates correct default", {
    to_apply <- tibble::tribble(
        ~Column, ~Function, ~Args, ~Output_colname,
        "FusionPrimerPCRDate", ~ suppressWarnings(min(.x, na.rm = TRUE)),
        NA, "{col}_min",
        "Value", ~ sum(.x, na.rm = TRUE), NA, "{col}_sum"
    )
    aggreg_meta <- .aggregate_meta(min_example,
        grouping_keys = agg_key_meta,
        function_tbl = to_apply
    )
    expected <- tibble::tribble(
        ~SubjectID, ~Tissue, ~FusionPrimerPCRDate_min, ~Value_sum,
        "S1", "BM", lubridate::ymd("2020-02-23"), 97086,
        "S1", "PB", lubridate::ymd("2020-02-23"), 14042,
        "S2", "BM", lubridate::ymd("2020-01-05"), 9827,
        "S2", "PB", lubridate::ymd("2020-01-05"), 2140
    )
    expect_equal(aggreg_meta |> dplyr::arrange(dplyr::desc(.data$Value_sum)),
        expected,
        ignore_attr = TRUE
    )
})

test_that(".aggregate_meta ignores column if not present", {
    to_apply <- tibble::tribble(
        ~Column, ~Function, ~Args, ~Output_colname,
        "FusionPrimerPCRDate", ~ suppressWarnings(min(.x, na.rm = TRUE)),
        NA, "{col}_min",
        "Value", ~ sum(.x, na.rm = TRUE), NA, "{col}_sum",
        "A", ~ sum(.x), NA, "{.col}_sum"
    )
    aggreg_meta <- .aggregate_meta(min_example,
        grouping_keys = agg_key_meta,
        function_tbl = to_apply
    )
    expect_true(!"A_sum" %in% colnames(aggreg_meta))
})

test_that(".aggregate_meta works for mixed function formula", {
    func_table <- tibble::tribble(
        ~Column, ~Function, ~Args, ~Output_colname,
        "FusionPrimerPCRDate", ~ min(.x, na.rm = TRUE), NA, "{.col}_min",
        "FusionPrimerPCRDate", min, list(na.rm = TRUE), "{.col}_min_fun"
    )
    aggreg_meta <- .aggregate_meta(min_example,
        grouping_keys = agg_key_meta,
        function_tbl = func_table
    )
    expect_true(all(aggreg_meta$FusionPrimerPCRDate_min ==
        aggreg_meta$FusionPrimerPCRDate_min_fun))
})

#------------------------------------------------------------------------------#
# Tests aggregate_metadata
#------------------------------------------------------------------------------#
test_that(".aggregate_meta warns deprecation of import_stats", {
    expect_deprecated({
        agg_meta <- aggregate_metadata(
            association_file = min_example,
            grouping_keys = agg_key_meta,
            import_stats = TRUE
        )
    })
})

test_that(".aggregate_meta warns if no columns found", {
    func_table <- tibble::tribble(
        ~Column, ~Function, ~Args, ~Output_colname,
        "A", ~ min(.x, na.rm = TRUE), NA, "{.col}_min",
        "B", min, list(na.rm = TRUE), "{.col}_min_fun"
    )
    expect_message({
        aggreg_meta <- aggregate_metadata(
            association_file = min_example,
            grouping_keys = agg_key_meta,
            aggregating_functions = func_table
        )
    })
    expect_null(aggreg_meta)
})

#------------------------------------------------------------------------------#
# .aggregate_lambda
#------------------------------------------------------------------------------#
test_that(".aggregate_lambda aggregates correctly", {
    df <- tibble::tribble(
        ~chr, ~integration_locus, ~strand, ~CompleteAmplificationID, ~Value,
        "1", 34255, "+", "ID1", 567,
        "2", 35675, "+", "ID1", 678,
        "1", 34255, "+", "ID2", 34,
        "1", 87953, "-", "ID3", 234,
        "5", 65789, "+", "ID3", 213,
        "1", 87953, "-", "ID4", 1233,
        "1", 32311, "-", "ID4", 52,
    )
    af <- tibble::tribble(
        ~CompleteAmplificationID, ~SubjectID,
        "ID1", "S1",
        "ID2", "S1",
        "ID3", "S2",
        "ID4", "S2"
    )
    agg <- aggregate_values_by_key(
        x = df, association_file = af, value_cols = "Value",
        key = "SubjectID", group = c(mandatory_IS_vars())
    )
    expected <- tibble::tribble(
        ~chr, ~integration_locus, ~strand, ~SubjectID, ~Value_sum,
        "1", 32311, "-", "S2", 52,
        "1", 34255, "+", "S1", 601,
        "1", 87953, "-", "S2", 1467,
        "2", 35675, "+", "S1", 678,
        "5", 65789, "+", "S2", 213
    )
    expect_equal(
        agg |> dplyr::arrange(.data$Value_sum),
        expected |> dplyr::arrange(.data$Value_sum)
    )
})

#------------------------------------------------------------------------------#
# Tests aggregate_values_by_key
#------------------------------------------------------------------------------#
test_that("aggregate_values_by_key succeeds for single IS matrix", {
    agg <- aggregate_values_by_key(
        integration_matrices,
        association_file,
        value_cols = c("seqCount", "fragmentEstimate")
    )
    expect_true(all(c(
        "chr", "integration_locus", "strand", "GeneName",
        "GeneStrand", "SubjectID", "CellMarker", "Tissue", "TimePoint",
        "fragmentEstimate_sum", "seqCount_sum"
    ) %in% colnames(agg)))
    expect_equal(nrow(agg), 1074)
    expect_equal(ncol(agg), 11)
})

test_that("aggregate_values_by_key succeeds for list of IS matrices", {
    matrix_list <- separate_quant_matrices(integration_matrices)
    agg <- aggregate_values_by_key(matrix_list, association_file)
    expect_true(all(c(
        "chr", "integration_locus", "strand", "GeneName",
        "GeneStrand", "SubjectID", "CellMarker", "Tissue", "TimePoint",
        "Value_sum"
    ) %in% colnames(agg[[1]])))
    expect_true(all(c(
        "chr", "integration_locus", "strand", "GeneName",
        "GeneStrand", "SubjectID", "CellMarker", "Tissue", "TimePoint",
        "Value_sum"
    ) %in% colnames(agg[[2]])))
})
