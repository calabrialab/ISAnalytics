library(ISAnalytics)

#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
op <- withr::local_options(
    ISAnalytics.reports = FALSE,
    ISAnalytics.verbose = FALSE
)

minimal_test_coll <- tibble::tribble(
    ~chr, ~integration_locus, ~strand, ~CompleteAmplificationID, ~seqCount,
    ~fragmentEstimate,
    "1", 62432, "+", "SAMPLE1", 654, 543.12,
    "6", 54632, "-", "SAMPLE1", 64, 56.65,
    "4", 23435, "+", "SAMPLE1", 865, 422.56,
    "4", 23435, "+", "SAMPLE2", 65, 15.98,
    "3", 57678, "-", "SAMPLE2", 14, 11.43,
    "4", 23435, "+", "SAMPLE3", 874, 123.4,
)

minimal_test_coll_meta <- tibble::tribble(
    ~ProjectID, ~PoolID, ~SubjectID, ~SequencingDate, ~ReplicateNumber,
    ~CompleteAmplificationID,
    "PJ1", "POOL1", "SJ001", lubridate::as_date(x = "2020-03-11"), 1, "SAMPLE1",
    "PJ1", "POOL1", "SJ001", lubridate::as_date(x = "2020-03-11"), 1, "SAMPLE2",
    "PJ1", "POOL1", "SJ002", lubridate::as_date(x = "2020-03-11"), 1, "SAMPLE3"
)

minimal_test_coll_meta_probs <- minimal_test_coll_meta %>%
    tibble::add_case(
        ProjectID = "PJ1",
        PoolID = "POOL1",
        SubjectID = "SJ002",
        SequencingDate = lubridate::as_date(x = "2020-03-11"),
        ReplicateNumber = 1,
        CompleteAmplificationID = "SAMPLE4"
    )

#------------------------------------------------------------------------------#
# Tests .check_same_info
#------------------------------------------------------------------------------#
test_that(".check_same_info returns empty if no problems", {
    check <- .check_same_info(minimal_test_coll_meta, minimal_test_coll)
    expect_true(nrow(check$miss) == 0)
})

test_that(".check_same_info returns non empty if more info", {
    check <- .check_same_info(minimal_test_coll_meta_probs, minimal_test_coll)
    expect_true(nrow(check$miss) == 1)
    expect_true(check$miss$ProjectID == "PJ1" &
        check$miss$PoolID == "POOL1" &
        check$miss$CompleteAmplificationID == "SAMPLE4" &
        check$miss$SubjectID == "SJ002")
})

#------------------------------------------------------------------------------#
# Tests .identify_independent_samples
#------------------------------------------------------------------------------#
test_that(".identify_independent_samples splits joined_df", {
    joined <- .join_matrix_af(minimal_test_coll,
        association_file = minimal_test_coll_meta,
        date_col = "SequencingDate"
    )
    splitted <- .identify_independent_samples(joined)
    expect_true(nrow(splitted$collisions) +
        nrow(splitted$non_collisions) == nrow(joined))
    expect_true(nrow(dplyr::inner_join(splitted[[1]], splitted[[2]],
        by = c(
            "chr",
            "integration_locus",
            "strand"
        )
    )) == 0)
    expect_true(unique(splitted$collisions$chr) == "4" &
        unique(splitted$collisions$integration_locus) == 23435 &
        unique(splitted$collisions$strand) == "+")
})

#------------------------------------------------------------------------------#
# Tests .discriminate_by_date
#------------------------------------------------------------------------------#
### OTHER VARS ###
example_df_nested <- function(ProjectID,
    Value,
    SequencingDate,
    PoolID,
    SubjectID,
    Replicate) {
    t <- tibble::tibble(
        Value = Value,
        SequencingDate = SequencingDate,
        ProjectID = ProjectID,
        PoolID = PoolID,
        SubjectID = SubjectID,
        ReplicateNumber = Replicate
    )
    t <- t %>%
        dplyr::mutate(
            CompleteAmplificationID = paste(.data$ProjectID,
                .data$PoolID,
                .data$SubjectID,
                .data$ReplicateNumber,
                sep = "_"
            ),
            .before = .data$Value
        )
}

example_df_nested_multi <- function(ProjectID,
    seqCount,
    fragmentEstimate,
    SequencingDate,
    PoolID,
    SubjectID,
    Replicate) {
    t <- tibble::tibble(
        seqCount = seqCount,
        fragmentEstimate = fragmentEstimate,
        SequencingDate = SequencingDate,
        ProjectID = ProjectID,
        PoolID = PoolID,
        SubjectID = SubjectID,
        ReplicateNumber = Replicate
    )
    t <- t %>%
        dplyr::mutate(
            CompleteAmplificationID = paste(.data$ProjectID,
                .data$PoolID,
                .data$SubjectID,
                .data$ReplicateNumber,
                sep = "_"
            ),
            .before = .data$seqCount
        )
}

test_that(".discriminate_by_date returns as expected for all equal dates", {
    dates <- lubridate::dmy("05/08/2020")
    df <- example_df_nested(
        ProjectID = "CLOEXP", Value = c(125, 25, 2),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c("subj1", "subj2", "subj3"),
        Replicate = c(1, 1, 2)
    )
    result <- .discriminate_by_date(df, "SequencingDate")
    expect_true(result$check == FALSE)
    expect_equal(result$data, df)
    df <- example_df_nested_multi(
        ProjectID = "CLOEXP", seqCount = c(125, 25, 2),
        fragmentEstimate = c(1354.54, 484.5, 386),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c("subj1", "subj2", "subj3"),
        Replicate = c(1, 1, 2)
    )
    result <- .discriminate_by_date(df, "SequencingDate")
    expect_true(result$check == FALSE)
    expect_equal(result$data, df)
})

test_that(".discriminate_by_date returns as expected for paired equal dates", {
    dates <- c(
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("01/08/2020"),
        lubridate::dmy("01/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("03/08/2020")
    )
    df <- example_df_nested(
        ProjectID = "CLOEXP",
        Value = c(125, 25, 2, 50, 800, 6),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c(
            "subj1", "subj2", "subj3", "subj4",
            "subj2", "subj5"
        ),
        Replicate = c(1, 1, 2, 1, 2, 3)
    )
    result <- .discriminate_by_date(df, "SequencingDate")
    expect_true(result$check == FALSE)
    df <- df %>% dplyr::arrange(.data$SequencingDate)
    expect_equal(result$data, df)
    df <- example_df_nested_multi(
        ProjectID = "CLOEXP",
        seqCount = c(125, 25, 2, 50, 800, 6),
        fragmentEstimate = c(456, 56, 786.87, 644.56, 4.857, 86.563),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c(
            "subj1", "subj2", "subj3", "subj4",
            "subj2", "subj5"
        ),
        Replicate = c(1, 1, 2, 1, 2, 3)
    )
    result <- .discriminate_by_date(df, "SequencingDate")
    expect_true(result$check == FALSE)
    df <- df %>% dplyr::arrange(.data$SequencingDate)
    expect_equal(result$data, df)
})

test_that(".discriminate_by_date returns as expected for single min date", {
    dates <- c(
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("01/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("03/08/2020")
    )
    df <- example_df_nested(
        ProjectID = "CLOEXP",
        Value = c(125, 25, 2, 50, 800, 6),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c(
            "subj1", "subj2", "subj3", "subj4",
            "subj2", "subj5"
        ),
        Replicate = c(1, 1, 2, 1, 2, 3)
    )
    result <- .discriminate_by_date(df, "SequencingDate")
    expect_true(result$check == TRUE)
    df <- df %>%
        dplyr::filter(.data$SequencingDate == lubridate::dmy("01/08/2020"))
    expect_equal(result$data, df)

    df <- example_df_nested_multi(
        ProjectID = "CLOEXP",
        seqCount = c(125, 25, 2, 50, 800, 6),
        fragmentEstimate = c(456, 56, 786.87, 644.56, 4.857, 86.563),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c(
            "subj1", "subj2", "subj3", "subj4",
            "subj2", "subj5"
        ),
        Replicate = c(1, 1, 2, 1, 2, 3)
    )
    result <- .discriminate_by_date(df, "SequencingDate")
    expect_true(result$check == TRUE)
    df <- df %>%
        dplyr::filter(.data$SequencingDate == lubridate::dmy("01/08/2020"))
    expect_equal(result$data, df)
})

#------------------------------------------------------------------------------#
# Tests .discriminate_by_replicate
#------------------------------------------------------------------------------#
test_that(".discriminate_by_replicate returns as expected for single max", {
    dates <- c(
        lubridate::dmy("01/08/2020"),
        lubridate::dmy("01/08/2020"),
        lubridate::dmy("03/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020")
    )
    df <- example_df_nested(
        ProjectID = "CLOEXP",
        Value = c(125, 25, 2, 50, 800, 6),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c(
            "subj1", "subj6", "subj3", "subj2",
            "subj2", "subj5"
        ),
        Replicate = c(1, 1, 2, 1, 2, 3)
    )
    result <- .discriminate_by_replicate(df)
    expect_true(result$check == TRUE)
    df <- df %>% dplyr::filter(.data$SubjectID == "subj2")
    expect_equal(result$data, df)

    df <- example_df_nested_multi(
        ProjectID = "CLOEXP",
        seqCount = c(125, 25, 2, 50, 800, 6),
        fragmentEstimate = c(456, 56, 786.87, 644.56, 4.857, 86.563),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c(
            "subj1", "subj6", "subj3", "subj2",
            "subj2", "subj5"
        ),
        Replicate = c(1, 1, 2, 1, 2, 3)
    )
    result <- .discriminate_by_replicate(df)
    expect_true(result$check == TRUE)
    df <- df %>% dplyr::filter(.data$SubjectID == "subj2")
    expect_equal(result$data, df)
})

test_that(".discriminate_by_replicate returns as expected for not max", {
    dates <- c(
        lubridate::dmy("01/08/2020"),
        lubridate::dmy("01/08/2020"),
        lubridate::dmy("03/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020")
    )
    df <- example_df_nested(
        ProjectID = "CLOEXP",
        Value = c(125, 25, 2, 50, 800, 6),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c(
            "subj1", "subj1", "subj3", "subj2",
            "subj2", "subj5"
        ),
        Replicate = c(1, 2, 2, 1, 2, 3)
    )
    result <- .discriminate_by_replicate(df)
    expect_true(result$check == FALSE)
    expect_equal(result$data, df)

    df <- example_df_nested_multi(
        ProjectID = "CLOEXP",
        seqCount = c(125, 25, 2, 50, 800, 6),
        fragmentEstimate = c(456, 56, 786.87, 644.56, 4.857, 86.563),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c(
            "subj1", "subj1", "subj3", "subj2",
            "subj2", "subj5"
        ),
        Replicate = c(1, 2, 2, 1, 2, 3)
    )
    result <- .discriminate_by_replicate(df)
    expect_true(result$check == FALSE)
    expect_equal(result$data, df)
})

#------------------------------------------------------------------------------#
# Tests .discriminate_by_seqCount
#------------------------------------------------------------------------------#
test_that(".discriminate_by_seqCount returns as expected for ratio < 10", {
    dates <- c(
        lubridate::dmy("01/08/2020"),
        lubridate::dmy("01/08/2020"),
        lubridate::dmy("03/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020")
    )
    df <- example_df_nested(
        ProjectID = "CLOEXP",
        Value = c(125, 25, 2, 50, 800, 6),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c(
            "subj1", "subj1", "subj3", "subj2",
            "subj2", "subj5"
        ),
        Replicate = c(1, 2, 2, 1, 2, 3)
    )
    result <- .discriminate_by_seqCount(df, 10, "Value")
    expect_true(result$check == FALSE)
    expect_equal(result$data, df)

    df <- example_df_nested_multi(
        ProjectID = "CLOEXP",
        seqCount = c(125, 25, 2, 50, 800, 6),
        fragmentEstimate = c(456, 56, 786.87, 644.56, 4.857, 86.563),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c(
            "subj1", "subj1", "subj3", "subj2",
            "subj2", "subj5"
        ),
        Replicate = c(1, 2, 2, 1, 2, 3)
    )
    result <- .discriminate_by_seqCount(df, 10, "seqCount")
    expect_true(result$check == FALSE)
    expect_equal(result$data, df)
})

test_that(".discriminate_by_seqCount returns as expected for ratio > 10", {
    dates <- c(
        lubridate::dmy("01/08/2020"),
        lubridate::dmy("01/08/2020"),
        lubridate::dmy("03/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020")
    )
    df <- example_df_nested(
        ProjectID = "CLOEXP",
        Value = c(20, 25, 2, 50, 800, 6),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c(
            "subj1", "subj1", "subj3", "subj2",
            "subj2", "subj5"
        ),
        Replicate = c(1, 2, 2, 1, 2, 3)
    )
    result <- .discriminate_by_seqCount(df, 10, "Value")
    expect_true(result$check == TRUE)
    df <- df %>% dplyr::filter(.data$SubjectID == "subj2")
    expect_equal(result$data, df)

    df <- example_df_nested_multi(
        ProjectID = "CLOEXP",
        seqCount = c(20, 25, 2, 50, 800, 6),
        fragmentEstimate = c(456, 56, 786.87, 644.56, 4.857, 86.563),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c(
            "subj1", "subj1", "subj3", "subj2",
            "subj2", "subj5"
        ),
        Replicate = c(1, 2, 2, 1, 2, 3)
    )
    result <- .discriminate_by_seqCount(df, 10, "seqCount")
    expect_true(result$check == TRUE)
    df <- df %>% dplyr::filter(.data$SubjectID == "subj2")
    expect_equal(result$data, df)
})

#------------------------------------------------------------------------------#
# Tests .four_step_check
#------------------------------------------------------------------------------#
test_that(".four_step_check returns as expected for first step", {
    dates <- c(
        lubridate::dmy("01/08/2020"),
        lubridate::dmy("03/08/2020"),
        lubridate::dmy("03/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020")
    )
    df <- example_df_nested(
        ProjectID = "CLOEXP",
        Value = c(20, 25, 2, 50, 800, 6),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c(
            "subj1", "subj6", "subj3", "subj2",
            "subj2", "subj5"
        ),
        Replicate = c(1, 2, 2, 1, 2, 3)
    )
    ex <- tibble::tibble(
        chr = 1, integration_locus = 10493, strand = "+",
        data = list(df)
    )
    res <- purrr::pmap(ex,
        .f = .four_step_check, date_col = "SequencingDate",
        reads_ratio = 10, seqCount_col = "Value"
    )
    rem_col <- purrr::map(res, function(x) {
        x$data
    })
    rem_col <- purrr::reduce(rem_col, dplyr::bind_rows)
    reassigned <- purrr::map(res, function(x) {
        x$reassigned
    })
    reassigned <- purrr::reduce(reassigned, sum)
    removed <- purrr::map(res, function(x) {
        x$removed
    })
    removed <- purrr::reduce(removed, sum)
    expect_equal(rem_col$SequencingDate, c(lubridate::dmy("01/08/2020")))
    expect_equal(rem_col$SubjectID, c("subj1"))
    expect_true(nrow(rem_col) == 1)
    expect_true(reassigned == 1)
    expect_true(removed == 0)

    df <- example_df_nested_multi(
        ProjectID = "CLOEXP",
        seqCount = c(20, 25, 2, 50, 800, 6),
        fragmentEstimate = c(456, 56, 786.87, 644.56, 4.857, 86.563),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c(
            "subj1", "subj6", "subj3", "subj2",
            "subj2", "subj5"
        ),
        Replicate = c(1, 2, 2, 1, 2, 3)
    )
    ex <- tibble::tibble(
        chr = 1, integration_locus = 10493, strand = "+",
        data = list(df)
    )
    res <- purrr::pmap(ex,
        .f = .four_step_check, date_col = "SequencingDate",
        reads_ratio = 10, seqCount_col = "seqCount"
    )
    rem_col <- purrr::map(res, function(x) {
        x$data
    })
    rem_col <- purrr::reduce(rem_col, dplyr::bind_rows)
    reassigned <- purrr::map(res, function(x) {
        x$reassigned
    })
    reassigned <- purrr::reduce(reassigned, sum)
    removed <- purrr::map(res, function(x) {
        x$removed
    })
    removed <- purrr::reduce(removed, sum)
    expect_equal(rem_col$SequencingDate, c(lubridate::dmy("01/08/2020")))
    expect_equal(rem_col$SubjectID, c("subj1"))
    expect_true(nrow(rem_col) == 1)
    expect_true(reassigned == 1)
    expect_true(removed == 0)
})

test_that(".four_step_check returns as expected for second step", {
    dates <- c(
        lubridate::dmy("01/08/2020"),
        lubridate::dmy("01/08/2020"),
        lubridate::dmy("03/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020")
    )
    df <- example_df_nested(
        ProjectID = "CLOEXP",
        Value = c(20, 25, 2, 50, 800, 6),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c(
            "subj1", "subj1", "subj3", "subj2",
            "subj4", "subj5"
        ),
        Replicate = c(1, 2, 2, 1, 2, 3)
    )
    ex <- tibble::tibble(
        chr = 1, integration_locus = 10493, strand = "+",
        data = list(df)
    )
    res <- purrr::pmap(ex,
        .f = .four_step_check, date_col = "SequencingDate",
        reads_ratio = 10, seqCount_col = "Value"
    )
    rem_col <- purrr::map(res, function(x) {
        x$data
    })
    rem_col <- purrr::reduce(rem_col, dplyr::bind_rows)
    reassigned <- purrr::map(res, function(x) {
        x$reassigned
    })
    reassigned <- purrr::reduce(reassigned, sum)
    removed <- purrr::map(res, function(x) {
        x$removed
    })
    removed <- purrr::reduce(removed, sum)
    expect_equal(rem_col$ReplicateNumber, c(1, 2))
    expect_equal(rem_col$SubjectID, c("subj1", "subj1"))
    expect_true(nrow(rem_col) == 2)
    expect_true(reassigned == 1)
    expect_true(removed == 0)

    df <- example_df_nested_multi(
        ProjectID = "CLOEXP",
        seqCount = c(20, 25, 2, 50, 800, 6),
        fragmentEstimate = c(456, 56, 786.87, 644.56, 4.857, 86.563),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c(
            "subj1", "subj1", "subj3", "subj2",
            "subj4", "subj5"
        ),
        Replicate = c(1, 2, 2, 1, 2, 3)
    )
    ex <- tibble::tibble(
        chr = 1, integration_locus = 10493, strand = "+",
        data = list(df)
    )
    res <- purrr::pmap(ex,
        .f = .four_step_check, date_col = "SequencingDate",
        reads_ratio = 10, seqCount_col = "seqCount"
    )
    rem_col <- purrr::map(res, function(x) {
        x$data
    })
    rem_col <- purrr::reduce(rem_col, dplyr::bind_rows)
    reassigned <- purrr::map(res, function(x) {
        x$reassigned
    })
    reassigned <- purrr::reduce(reassigned, sum)
    removed <- purrr::map(res, function(x) {
        x$removed
    })
    removed <- purrr::reduce(removed, sum)
    expect_equal(rem_col$ReplicateNumber, c(1, 2))
    expect_equal(rem_col$SubjectID, c("subj1", "subj1"))
    expect_true(nrow(rem_col) == 2)
    expect_true(reassigned == 1)
    expect_true(removed == 0)
})

test_that(".four_step_check returns as expected for third step", {
    dates <- c(
        lubridate::dmy("01/08/2020"),
        lubridate::dmy("01/08/2020"),
        lubridate::dmy("03/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020")
    )
    df <- example_df_nested(
        ProjectID = "CLOEXP",
        Value = c(20, 25, 2, 50, 800, 6),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c(
            "subj1", "subj1", "subj3", "subj2",
            "subj2", "subj5"
        ),
        Replicate = c(1, 2, 2, 1, 2, 3)
    )
    ex <- tibble::tibble(
        chr = 1, integration_locus = 10493, strand = "+",
        data = list(df)
    )
    res <- purrr::pmap(ex,
        .f = .four_step_check, date_col = "SequencingDate",
        reads_ratio = 10, seqCount_col = "Value"
    )
    rem_col <- purrr::map(res, function(x) {
        x$data
    })
    rem_col <- purrr::reduce(rem_col, dplyr::bind_rows)
    reassigned <- purrr::map(res, function(x) {
        x$reassigned
    })
    reassigned <- purrr::reduce(reassigned, sum)
    removed <- purrr::map(res, function(x) {
        x$removed
    })
    removed <- purrr::reduce(removed, sum)
    expect_equal(rem_col$ReplicateNumber, c(1, 2))
    expect_equal(rem_col$Value, c(50, 800))
    expect_equal(rem_col$SubjectID, c("subj2", "subj2"))
    expect_true(nrow(rem_col) == 2)
    expect_true(reassigned == 1)
    expect_true(removed == 0)

    df <- example_df_nested_multi(
        ProjectID = "CLOEXP",
        seqCount = c(20, 25, 2, 50, 800, 6),
        fragmentEstimate = c(456, 56, 786.87, 644.56, 4.857, 86.563),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c(
            "subj1", "subj1", "subj3", "subj2",
            "subj2", "subj5"
        ),
        Replicate = c(1, 2, 2, 1, 2, 3)
    )
    ex <- tibble::tibble(
        chr = 1, integration_locus = 10493, strand = "+",
        data = list(df)
    )
    res <- purrr::pmap(ex,
        .f = .four_step_check, date_col = "SequencingDate",
        reads_ratio = 10, seqCount_col = "seqCount"
    )
    rem_col <- purrr::map(res, function(x) {
        x$data
    })
    rem_col <- purrr::reduce(rem_col, dplyr::bind_rows)
    reassigned <- purrr::map(res, function(x) {
        x$reassigned
    })
    reassigned <- purrr::reduce(reassigned, sum)
    removed <- purrr::map(res, function(x) {
        x$removed
    })
    removed <- purrr::reduce(removed, sum)
    expect_equal(rem_col$ReplicateNumber, c(1, 2))
    expect_equal(rem_col$seqCount, c(50, 800))
    expect_equal(rem_col$SubjectID, c("subj2", "subj2"))
    expect_true(nrow(rem_col) == 2)
    expect_true(reassigned == 1)
    expect_true(removed == 0)
})

test_that(".four_step_check returns as expected for fourth step", {
    dates <- c(
        lubridate::dmy("01/08/2020"),
        lubridate::dmy("01/08/2020"),
        lubridate::dmy("03/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020")
    )
    df <- example_df_nested(
        ProjectID = "CLOEXP",
        Value = c(20, 25, 2, 50, 80, 6),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c(
            "subj1", "subj1", "subj3", "subj2",
            "subj2", "subj5"
        ),
        Replicate = c(1, 2, 2, 1, 2, 3)
    )
    ex <- tibble::tibble(
        chr = 1, integration_locus = 10493, strand = "+",
        data = list(df)
    )
    res <- purrr::pmap(ex,
        .f = .four_step_check, date_col = "SequencingDate",
        reads_ratio = 10, seqCount_col = "Value"
    )
    rem_col <- purrr::map(res, function(x) {
        x$data
    })
    rem_col <- purrr::reduce(rem_col, dplyr::bind_rows)
    reassigned <- purrr::map(res, function(x) {
        x$reassigned
    })
    reassigned <- purrr::reduce(reassigned, sum)
    removed <- purrr::map(res, function(x) {
        x$removed
    })
    removed <- purrr::reduce(removed, sum)
    expect_null(rem_col)
    expect_true(reassigned == 0)
    expect_true(removed == 1)

    df <- example_df_nested_multi(
        ProjectID = "CLOEXP",
        seqCount = c(20, 25, 2, 50, 80, 6),
        fragmentEstimate = c(456, 56, 786.87, 644.56, 4.857, 86.563),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c(
            "subj1", "subj1", "subj3", "subj2",
            "subj2", "subj5"
        ),
        Replicate = c(1, 2, 2, 1, 2, 3)
    )
    ex <- tibble::tibble(
        chr = 1, integration_locus = 10493, strand = "+",
        data = list(df)
    )
    res <- purrr::pmap(ex,
        .f = .four_step_check, date_col = "SequencingDate",
        reads_ratio = 10, seqCount_col = "seqCount"
    )
    rem_col <- purrr::map(res, function(x) {
        x$data
    })
    rem_col <- purrr::reduce(rem_col, dplyr::bind_rows)
    reassigned <- purrr::map(res, function(x) {
        x$reassigned
    })
    reassigned <- purrr::reduce(reassigned, sum)
    removed <- purrr::map(res, function(x) {
        x$removed
    })
    removed <- purrr::reduce(removed, sum)
    expect_null(rem_col)
    expect_true(reassigned == 0)
    expect_true(removed == 1)
})

#------------------------------------------------------------------------------#
# Tests .coll_mapping
#------------------------------------------------------------------------------#
### OTHER VARS ###
example_collisions <- function() {
    cols <- list(
        chr = c(1, 2, 3, 4),
        integration_locus = c(103948, 14390, 12453, 12353),
        strand = c("+", "+", "-", "-")
    )
    # Can discriminate by date
    dates <- c(
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("01/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("03/08/2020")
    )
    smpl1 <- example_df_nested(
        ProjectID = "CLOEXP",
        Value = c(125, 25, 2, 50, 800, 6),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c(
            "subj1", "subj2", "subj3", "subj4",
            "subj2", "subj5"
        ),
        Replicate = c(1, 1, 2, 1, 2, 3)
    )
    # Can discriminate by replicate
    dates <- c(
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020")
    )
    smpl2 <- example_df_nested(
        ProjectID = "CLOEXP",
        Value = c(125, 25, 2, 50, 800, 6),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c(
            "subj1", "subj2", "subj3", "subj4",
            "subj2", "subj5"
        ),
        Replicate = c(1, 1, 2, 1, 2, 3)
    )
    # Can discriminate by seqCount
    smpl3 <- example_df_nested(
        ProjectID = "CLOEXP",
        Value = c(25, 25, 2, 50, 800, 6),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c(
            "subj1", "subj2", "subj3", "subj4",
            "subj2", "subj1"
        ),
        Replicate = c(1, 1, 2, 1, 2, 3)
    )
    # Can't discriminate
    smpl4 <- example_df_nested(
        ProjectID = "CLOEXP",
        Value = c(125, 25, 2, 50, 800, 6),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c(
            "subj1", "subj2", "subj3", "subj4",
            "subj2", "subj1"
        ),
        Replicate = c(1, 1, 2, 1, 2, 3)
    )
    t <- tibble::as_tibble(cols)
    t <- t %>%
        tibble::add_column(tibble::as_tibble_col(list(
            smpl1, smpl2,
            smpl3, smpl4
        ),
        column_name = "data"
        ))
    t
}

example_collisions_multi <- function() {
    cols <- list(
        chr = c(1, 2, 3, 4),
        integration_locus = c(103948, 14390, 12453, 12353),
        strand = c("+", "+", "-", "-")
    )
    # Can discriminate by date
    dates <- c(
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("01/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("03/08/2020")
    )
    smpl1 <- example_df_nested_multi(
        ProjectID = "CLOEXP",
        seqCount = c(125, 25, 2, 50, 800, 6),
        fragmentEstimate = c(456, 56, 786.87, 644.56, 4.857, 86.563),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c(
            "subj1", "subj2", "subj3", "subj4",
            "subj2", "subj5"
        ),
        Replicate = c(1, 1, 2, 1, 2, 3)
    )
    # Can discriminate by replicate
    dates <- c(
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020"),
        lubridate::dmy("05/08/2020")
    )
    smpl2 <- example_df_nested_multi(
        ProjectID = "CLOEXP",
        seqCount = c(125, 25, 2, 50, 800, 6),
        fragmentEstimate = c(456, 56, 786.87, 644.56, 4.857, 86.563),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c(
            "subj1", "subj2", "subj3", "subj4",
            "subj2", "subj5"
        ),
        Replicate = c(1, 1, 2, 1, 2, 3)
    )
    # Can discriminate by seqCount
    smpl3 <- example_df_nested_multi(
        ProjectID = "CLOEXP",
        seqCount = c(25, 25, 2, 50, 800, 6),
        fragmentEstimate = c(456, 56, 786.87, 644.56, 4.857, 86.563),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c(
            "subj1", "subj2", "subj3", "subj4",
            "subj2", "subj1"
        ),
        Replicate = c(1, 1, 2, 1, 2, 3)
    )
    # Can't discriminate
    smpl4 <- example_df_nested_multi(
        ProjectID = "CLOEXP",
        seqCount = c(125, 25, 2, 50, 800, 6),
        fragmentEstimate = c(456, 56, 786.87, 644.56, 4.857, 86.563),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c(
            "subj1", "subj2", "subj3", "subj4",
            "subj2", "subj1"
        ),
        Replicate = c(1, 1, 2, 1, 2, 3)
    )
    t <- tibble::as_tibble(cols)
    t <- t %>%
        tibble::add_column(tibble::as_tibble_col(list(
            smpl1, smpl2,
            smpl3, smpl4
        ),
        column_name = "data"
        ))
    t
}

ex_collisions <- example_collisions()
ex_collisions_multi <- example_collisions_multi()

test_that(".coll_mapping produces expected result for example collisions", {
    result <- .coll_mapping(ex_collisions, "SequencingDate", 10, "Value")
    coll <- result$coll
    removed <- result$removed
    reassigned <- result$reassigned
    expect_true(removed == 1)
    expect_true(reassigned == 3)
    # Removed integration should not be in the result
    rem_integration <- coll %>%
        dplyr::filter(
            .data$chr == 4, .data$integration_locus == 12353,
            .data$strand == "-"
        )
    expect_true(nrow(rem_integration) == 0)
    # First integration expected result
    first <- coll %>%
        dplyr::filter(
            .data$chr == 1, .data$integration_locus == 103948,
            .data$strand == "+"
        )
    expect_true(all(first$SubjectID == "subj3"))
    expect_true(nrow(first) == 1)
    # Second integration expected result
    second <- coll %>%
        dplyr::filter(
            .data$chr == 2, .data$integration_locus == 14390,
            .data$strand == "+"
        )
    expect_true(all(second$SubjectID == "subj2"))
    expect_true(nrow(second) == 2)
    # Third integration expected result
    third <- coll %>%
        dplyr::filter(
            .data$chr == 3, .data$integration_locus == 12453,
            .data$strand == "-"
        )
    expect_true(all(third$SubjectID == "subj2"))
    expect_true(nrow(third) == 2)

    result <- .coll_mapping(
        ex_collisions_multi, "SequencingDate", 10,
        "seqCount"
    )
    coll <- result$coll
    removed <- result$removed
    reassigned <- result$reassigned
    expect_true(removed == 1)
    expect_true(reassigned == 3)
    # Removed integration should not be in the result
    rem_integration <- coll %>%
        dplyr::filter(
            .data$chr == 4, .data$integration_locus == 12353,
            .data$strand == "-"
        )
    expect_true(nrow(rem_integration) == 0)
    # First integration expected result
    first <- coll %>%
        dplyr::filter(
            .data$chr == 1, .data$integration_locus == 103948,
            .data$strand == "+"
        )
    expect_true(all(first$SubjectID == "subj3"))
    expect_true(nrow(first) == 1)
    # Second integration expected result
    second <- coll %>%
        dplyr::filter(
            .data$chr == 2, .data$integration_locus == 14390,
            .data$strand == "+"
        )
    expect_true(all(second$SubjectID == "subj2"))
    expect_true(nrow(second) == 2)
    # Third integration expected result
    third <- coll %>%
        dplyr::filter(
            .data$chr == 3, .data$integration_locus == 12453,
            .data$strand == "-"
        )
    expect_true(all(third$SubjectID == "subj2"))
    expect_true(nrow(third) == 2)
})

#------------------------------------------------------------------------------#
# Tests .process_collisions
#------------------------------------------------------------------------------#
test_that(".process_collisions returns updated collisions", {
    result <- .process_collisions(
        ex_collisions %>%
            tidyr::unnest(.data$data),
        "SequencingDate", 10, "Value",
        NULL
    )
    coll <- result$coll
    removed <- result$removed
    reassigned <- result$reassigned
    expect_true(removed == 1)
    expect_true(reassigned == 3)
    # Removed integration should not be in the result
    rem_integration <- coll %>%
        dplyr::filter(
            .data$chr == 4, .data$integration_locus == 12353,
            .data$strand == "-"
        )
    expect_true(nrow(rem_integration) == 0)
    # First integration expected result
    first <- coll %>%
        dplyr::filter(
            .data$chr == 1, .data$integration_locus == 103948,
            .data$strand == "+"
        )
    expect_true(all(first$SubjectID == "subj3"))
    expect_true(nrow(first) == 1)
    # Second integration expected result
    second <- coll %>%
        dplyr::filter(
            .data$chr == 2, .data$integration_locus == 14390,
            .data$strand == "+"
        )
    expect_true(all(second$SubjectID == "subj2"))
    expect_true(nrow(second) == 2)
    # Third integration expected result
    third <- coll %>%
        dplyr::filter(
            .data$chr == 3, .data$integration_locus == 12453,
            .data$strand == "-"
        )
    expect_true(all(third$SubjectID == "subj2"))
    expect_true(nrow(third) == 2)

    result <- .process_collisions(
        ex_collisions_multi %>%
            tidyr::unnest(.data$data),
        "SequencingDate", 10, "seqCount",
        NULL
    )
    coll <- result$coll
    removed <- result$removed
    reassigned <- result$reassigned
    expect_true(removed == 1)
    expect_true(reassigned == 3)
    # Removed integration should not be in the result
    rem_integration <- coll %>%
        dplyr::filter(
            .data$chr == 4, .data$integration_locus == 12353,
            .data$strand == "-"
        )
    expect_true(nrow(rem_integration) == 0)
    # First integration expected result
    first <- coll %>%
        dplyr::filter(
            .data$chr == 1, .data$integration_locus == 103948,
            .data$strand == "+"
        )
    expect_true(all(first$SubjectID == "subj3"))
    expect_true(nrow(first) == 1)
    # Second integration expected result
    second <- coll %>%
        dplyr::filter(
            .data$chr == 2, .data$integration_locus == 14390,
            .data$strand == "+"
        )
    expect_true(all(second$SubjectID == "subj2"))
    expect_true(nrow(second) == 2)
    # Third integration expected result
    third <- coll %>%
        dplyr::filter(
            .data$chr == 3, .data$integration_locus == 12453,
            .data$strand == "-"
        )
    expect_true(all(third$SubjectID == "subj2"))
    expect_true(nrow(third) == 2)
})

#------------------------------------------------------------------------------#
# Tests remove_collisions
#------------------------------------------------------------------------------#
test_that("remove_collisions succeeds", {
    invisible(capture_output({
        coll_rem <- remove_collisions(minimal_test_coll,
            minimal_test_coll_meta_probs,
            date_col = "SequencingDate",
            reads_ratio = 10
        )
        expected <- coll_rem %>%
            dplyr::filter(
                .data$chr == "4",
                .data$integration_locus == 23435,
                .data$strand == "+"
            ) %>%
            dplyr::distinct(.data$CompleteAmplificationID)
        expect_true(
            all(expected$CompleteAmplificationID %in% c("SAMPLE1", "SAMPLE2"))
        )
    }))
})

#------------------------------------------------------------------------------#
# Tests realign_after_collisions
#------------------------------------------------------------------------------#
test_that("realign_after_collisions correctly re-aligns", {
    separated <- separate_quant_matrices(minimal_test_coll,
        key = c(
            mandatory_IS_vars(),
            "CompleteAmplificationID"
        )
    )
    coll_single <- remove_collisions(
        x = separated$seqCount,
        association_file = minimal_test_coll_meta,
        quant_cols = c(seqCount = "Value")
    )
    realigned <- realign_after_collisions(
        coll_single,
        list(fragmentEstimate = separated$fragmentEstimate)
    )
    expect_true(nrow(coll_single) == nrow(realigned$fragmentEstimate))
    expect_true(all(realigned$fragmentEstimate$chr %in% coll_single$chr))
    expect_true(all(realigned$fragmentEstimate$integration_locus
        %in% coll_single$integration_locus))
    expect_true(all(realigned$fragmentEstimate$CompleteAmplificationID
        %in% coll_single$CompleteAmplificationID))
})
