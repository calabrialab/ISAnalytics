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

sample_tag_cols <- tibble::tribble(
    ~names, ~types, ~transform, ~flag, ~tag,
    "CompleteAmplificationID", "char", NULL, "required", "pcr_repl_id",
    "ReplicateNumber", "int", NULL, "required", "pcr_replicate",
    "PoolID", "char", NULL, "required", "pool_id",
    "ProjectID", "char", NULL, "required", "project_id"
)

#------------------------------------------------------------------------------#
# Tests .check_same_info
#------------------------------------------------------------------------------#
test_that(".check_same_info returns empty if no problems", {
    check <- .check_same_info(minimal_test_coll_meta,
        minimal_test_coll,
        req_tag_cols = sample_tag_cols,
        indep_sample_id = c("ProjectID", "SubjectID")
    )
    expect_true(nrow(check$miss) == 0)
})

test_that(".check_same_info returns non empty if more info", {
    check <- .check_same_info(minimal_test_coll_meta_probs,
        minimal_test_coll,
        req_tag_cols = sample_tag_cols,
        indep_sample_id = c("ProjectID", "SubjectID")
    )
    expect_true(nrow(check$miss) == 1)
    expect_true(check$miss$ProjectID == "PJ1" &
        check$miss$PoolID == "POOL1" &
        check$miss$CompleteAmplificationID == "SAMPLE4" &
        check$miss$SubjectID == "SJ002")
})

test_that(".check_same_info reports only projects of interest", {
    mod_af <- minimal_test_coll_meta_probs %>%
        tibble::add_case(
            ProjectID = "PJ2", PoolID = "POOL2",
            SubjectID = "SJ003",
            SequencingDate = lubridate::as_date(x = "2021-04-12"),
            ReplicateNumber = 1,
            CompleteAmplificationID = "SAMPLE5"
        ) %>%
        tibble::add_case(
            ProjectID = "PJ2", PoolID = "POOL2",
            SubjectID = "SJ004",
            SequencingDate = lubridate::as_date(x = "2021-04-12"),
            ReplicateNumber = 1,
            CompleteAmplificationID = "SAMPLE6"
        ) %>%
        tibble::add_case(
            ProjectID = "PJ3", PoolID = "POOL3",
            SubjectID = "SJ005",
            SequencingDate = lubridate::as_date(x = "2021-04-12"),
            ReplicateNumber = 1,
            CompleteAmplificationID = "SAMPLE7"
        )
    mod_matrix <- minimal_test_coll %>%
        tibble::add_case(
            chr = "7", integration_locus = 56753, strand = "+",
            CompleteAmplificationID = "SAMPLE5", seqCount = 432,
            fragmentEstimate = 453.5
        )
    check <- .check_same_info(mod_af,
        mod_matrix,
        req_tag_cols = sample_tag_cols,
        indep_sample_id = c("ProjectID", "SubjectID")
    )
    expect_false("SAMPLE7" %in% check$miss$CompleteAmplificationID)
    expect_true(all(c("SAMPLE4", "SAMPLE6") %in%
        check$miss$CompleteAmplificationID))
})

#------------------------------------------------------------------------------#
# Tests .identify_independent_samples
#------------------------------------------------------------------------------#
test_that(".identify_independent_samples splits joined_df", {
    joined <- minimal_test_coll %>%
      dplyr::left_join(minimal_test_coll_meta,
                       by = "CompleteAmplificationID") %>%
      dplyr::select(dplyr::all_of(c(
        colnames(minimal_test_coll), "SequencingDate",
        "ReplicateNumber", "ProjectID", "SubjectID"
      )))

    split <- .identify_independent_samples(joined,
        indep_sample_id = c(
            "ProjectID",
            "SubjectID"
        )
    )
    expect_true(nrow(split$collisions) +
        nrow(split$non_collisions) == nrow(joined))
    expect_true(nrow(dplyr::inner_join(split[[1]], split[[2]],
        by = c(
            "chr",
            "integration_locus",
            "strand"
        )
    )) == 0)
    expect_true(unique(split$collisions$chr) == "4" &
        unique(split$collisions$integration_locus) == 23435 &
        unique(split$collisions$strand) == "+")
})

test_that(".identify_independent_samples works with custom indep sample", {
    custom_matrix <- tibble::tribble(
        ~chr, ~integration_locus, ~strand, ~CompleteAmplificationID, ~seqCount,
        ~fragmentEstimate, ~Field1, ~Field2, ~Field3,
        "5", 32424, "+", "SAMPLE1", 564, 565.43, "A01", "B01", "C01",
        "5", 32424, "+", "SAMPLE2", 56, 13.4, "A01", "B01", "C02",
        "6", 43544, "-", "SAMPLE3", 564, 324.5, "A01", "B01", "C01",
        "6", 43544, "-", "SAMPLE4", 53, 46.9, "A01", "B01", "C01",
        "1", 54354, "+", "SAMPLE5", 67, 54.6, "A01", "B02", "C01",
        "2", 67544, "+", "SAMPLE6", 732, 563.8, "A01", "B02", "C01"
    )
    independent_sample_key <- c("Field1", "Field2", "Field3")
    split <- .identify_independent_samples(data.table::setDT(custom_matrix),
        indep_sample_id = independent_sample_key
    )
    expect_true(nrow(split$collisions) +
        nrow(split$non_collisions) == nrow(custom_matrix))
    expect_true(unique(split$collisions$chr) == "5" &
        unique(split$collisions$integration_locus) == 32424 &
        unique(split$collisions$strand) == "+")
})

#------------------------------------------------------------------------------#
# Tests .discriminate_by_date
#------------------------------------------------------------------------------#
example_df_multi <- function(ProjectID,
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
    t <- data.table::setDT(t)
    return(t)
}

test_that(".discriminate_by_date returns as expected for all equal dates", {
    dates <- lubridate::dmy("05/08/2020")
    df <- example_df_multi(
        ProjectID = "CLOEXP", seqCount = c(125, 25, 2),
        fragmentEstimate = c(1354.54, 484.5, 386),
        SequencingDate = dates, PoolID = "POOL6",
        SubjectID = c("subj1", "subj2", "subj3"),
        Replicate = c(1, 1, 2)
    )
    result <- .discriminate_by_date(
        df, "SequencingDate",
        c("ProjectID", "SubjectID")
    )
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
    df <- example_df_multi(
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
    result <- .discriminate_by_date(
        df, "SequencingDate",
        c("ProjectID", "SubjectID")
    )
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
    df <- example_df_multi(
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
    result <- .discriminate_by_date(
        df, "SequencingDate",
        c("ProjectID", "SubjectID")
    )
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
    df <- example_df_multi(
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
    result <- .discriminate_by_replicate(
        df, "ReplicateNumber",
        c("ProjectID", "SubjectID")
    )
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
    df <- example_df_multi(
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
    result <- .discriminate_by_replicate(
        df, "ReplicateNumber",
        c("ProjectID", "SubjectID")
    )
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

    df <- example_df_multi(
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
    result <- .discriminate_by_seqCount(
        df, 10, "seqCount",
        c("ProjectID", "SubjectID")
    )
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
    df <- example_df_multi(
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
    result <- .discriminate_by_seqCount(
        df, 10, "seqCount",
        c("ProjectID", "SubjectID")
    )
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
    df <- example_df_multi(
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
    ex <- df %>%
        dplyr::mutate(
            chr = "1",
            integration_locus = 10493,
            strand = "+",
            .before = "CompleteAmplificationID"
        )
    res <- .four_step_check(ex,
        date_col = "SequencingDate",
        repl_col = "ReplicateNumber",
        reads_ratio = 10,
        seqCount_col = "seqCount",
        ind_sample_key = c("ProjectID", "SubjectID")
    )
    rem_col <- res$data
    reassigned <- res$reassigned
    removed <- res$removed
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
    df <- example_df_multi(
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
    ex <- df %>%
        dplyr::mutate(
            chr = "1",
            integration_locus = 10493,
            strand = "+",
            .before = "CompleteAmplificationID"
        )

    res <- .four_step_check(ex,
        date_col = "SequencingDate",
        repl_col = "ReplicateNumber",
        reads_ratio = 10,
        seqCount_col = "seqCount",
        ind_sample_key = c("ProjectID", "SubjectID")
    )
    rem_col <- res$data
    reassigned <- res$reassigned
    removed <- res$removed
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
    df <- example_df_multi(
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
    ex <- df %>%
        dplyr::mutate(
            chr = "1",
            integration_locus = 10493,
            strand = "+",
            .before = "CompleteAmplificationID"
        )

    res <- .four_step_check(ex,
        date_col = "SequencingDate",
        repl_col = "ReplicateNumber",
        reads_ratio = 10,
        seqCount_col = "seqCount",
        ind_sample_key = c("ProjectID", "SubjectID")
    )
    rem_col <- res$data
    reassigned <- res$reassigned
    removed <- res$removed
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
    df <- example_df_multi(
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
    ex <- df %>%
        dplyr::mutate(
            chr = "1",
            integration_locus = 10493,
            strand = "+",
            .before = "CompleteAmplificationID"
        )

    res <- .four_step_check(ex,
        date_col = "SequencingDate",
        repl_col = "ReplicateNumber",
        reads_ratio = 10,
        seqCount_col = "seqCount",
        ind_sample_key = c("ProjectID", "SubjectID")
    )
    rem_col <- res$data
    reassigned <- res$reassigned
    removed <- res$removed
    expect_null(rem_col)
    expect_true(reassigned == 0)
    expect_true(removed == 1)
})

#------------------------------------------------------------------------------#
# Tests .process_collisions
#------------------------------------------------------------------------------#
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
    smpl1 <- example_df_multi(
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
    smpl2 <- example_df_multi(
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
    smpl3 <- example_df_multi(
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
    smpl4 <- example_df_multi(
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
        )) %>%
        tidyr::unnest(.data$data)
    return(data.table::setDT(t))
}

ex_collisions_multi <- example_collisions_multi()

test_that(".process_collisions returns updated collisions", {
    result <- .process_collisions(
        collisions = ex_collisions_multi,
        date_col = "SequencingDate",
        reads_ratio = 10,
        seqCount_col = "seqCount",
        repl_col = "ReplicateNumber",
        ind_sample_key = c("ProjectID", "SubjectID"),
        max_workers = 4
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
# Tests .collisions_check_input_matrices
#------------------------------------------------------------------------------#
test_that(paste(
    ".collisions_check_input_matrices throws errors",
    "for input list - matrix issues"
), {
    separated_m <- separate_quant_matrices(minimal_test_coll,
        key = c(
            mandatory_IS_vars(),
            "CompleteAmplificationID"
        )
    )
    issues1 <- separated_m
    issues1$fragmentEstimate <- issues1$fragmentEstimate %>%
        dplyr::rename(chrom = "chr")
    issues1$seqCount <- issues1$seqCount %>%
        dplyr::rename(sample_id = "CompleteAmplificationID")
    expect_error(
        {
            mode <- .collisions_check_input_matrices(issues1, c(
                seqCount = "seqCount",
                fragmentEstimate = "fragmentEstimate"
            ))
        },
        class = "coll_matrix_issues"
    )
})

test_that(paste(".collisions_check_input_matrices works"), {
    separated_m <- separate_quant_matrices(minimal_test_coll,
        key = c(
            mandatory_IS_vars(),
            "CompleteAmplificationID"
        )
    )
    mode <- .collisions_check_input_matrices(separated_m, c(
        seqCount = "seqCount",
        fragmentEstimate = "fragmentEstimate"
    ))
    expect_equal(mode, "LIST")
    mode <- .collisions_check_input_matrices(minimal_test_coll, c(
        seqCount = "seqCount",
        fragmentEstimate = "fragmentEstimate"
    ))
    expect_equal(mode, "MULTI")
})

#------------------------------------------------------------------------------#
# Tests .collisions_check_input_af
#------------------------------------------------------------------------------#
test_that(".collisions_check_input_af works as expected", {
    ## NULL or empty independent sample key throws error
    expect_error({
        checks <- .collisions_check_input_af(minimal_test_coll_meta,
            date_col = "SequencingDate",
            independent_sample_id = NULL
        )
    })
    expect_error({
        checks <- .collisions_check_input_af(minimal_test_coll_meta,
            date_col = "SequencingDate",
            independent_sample_id = c()
        )
    })
    ## Throws error if user input cols are not in af
    expect_error({
        checks <- .collisions_check_input_af(minimal_test_coll_meta,
            date_col = "SequencingDate",
            independent_sample_id = c("A", "B")
        )
    })
    ## Throws error if required tags are ok but actual names are not present
    ## in the data frame
    expect_error({
        checks <- .collisions_check_input_af(minimal_test_coll_meta %>%
            dplyr::rename(Project = "ProjectID"),
        date_col = "SequencingDate",
        independent_sample_id = c("SubjectID")
        )
    })
    ## Throws error if date column is not a date
    expect_error(
        {
            checks <- .collisions_check_input_af(minimal_test_coll_meta %>%
                dplyr::mutate(SequencingDate = as.character(
                  .data$SequencingDate)),
            date_col = "SequencingDate",
            independent_sample_id = c("SubjectID")
            )
        },
        class = "not_date_coll_err"
    )
    ## Throws error if date col contains NA
    expect_error({
        mod_af <- minimal_test_coll_meta
        mod_af[1, ]$SequencingDate <- NA
        checks <- .collisions_check_input_af(mod_af,
            date_col = "SequencingDate",
            independent_sample_id = c("ProjectID", "SubjectID")
        )
    })
})

#------------------------------------------------------------------------------#
# Tests .summary_input
#------------------------------------------------------------------------------#
test_that(".summary_input returns info correctly", {
    summary <- .summary_input(
        ex_collisions_multi,
        c("seqCount", "fragmentEstimate")
    )
    expect_equal(summary$total_iss, 4)
    expect_equal(summary$quant_totals$seqCount, 3932)
    expect_equal(summary$quant_totals$fragmentEstimate, 8139.4)
})

#------------------------------------------------------------------------------#
# Tests .per_pool_stats
#------------------------------------------------------------------------------#
test_that(".per_pool_stats works as expected", {
    pool_stats <- .per_pool_stats(
        ex_collisions_multi,
        c("seqCount", "fragmentEstimate"),
        "PoolID"
    )
    expect_true(nrow(pool_stats) == 1)
    expect_true(pool_stats$PoolID == "POOL6")
    expect_true(ncol(pool_stats) == 31)
})

#------------------------------------------------------------------------------#
# Tests remove_collisions
#------------------------------------------------------------------------------#
test_that("remove_collisions succeeds", {
    withr::local_options(list(ISAnalytics.verbose = FALSE))
    coll_rem <- remove_collisions(
        x = minimal_test_coll,
        association_file = minimal_test_coll_meta_probs,
        report_path = NULL, max_workers = 4
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
})

test_that("remove_collisions produces report", {
    withr::local_options(list(ISAnalytics.reports = TRUE))
    tmp_dir <- withr::local_tempdir()
    coll_rem <- remove_collisions(
        integration_matrices, association_file,
        report_path = tmp_dir
    )
    path_to_file <- fs::path(tmp_dir, .generate_report_filename("collisions"))
    expect_true(fs::file_exists(path_to_file))
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
