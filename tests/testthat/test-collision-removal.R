context("Collision Removal")

library(ISAnalytics)

#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
op <- options(ISAnalytics.widgets = FALSE)
on.exit(options(op))

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
af_missing_info <- function(af) {
    af %>%
        dplyr::filter(
            .data$CompleteAmplificationID != paste0(
                "CLOEXP_POOL6_LTR10LC10_VA2020-",
                "mix10_VA2020-mix10_lenti_737.",
                "pCCLsin.PPT.SFFV.eGFP.Wpre_",
                "PLATE1_NONE_1_NONE_SLiM_0000_",
                "NONE_NONE_OSR_NONE_NONE_NONE"
            )
        )
}

af_no_prob <- function(af) {
    af %>%
        dplyr::filter(
            .data$CompleteAmplificationID != paste0(
                "PROJECT1101_ABY-LR",
                "-PL4-POOL54_LTR51LC51_VA2020-",
                "mix14_VA2020-mix14_lenti_737.",
                "pCCLsin.PPT.SFFV.eGFP.Wpre_PL4_",
                "NONE_2_NONE_SLiM_0_NONE_NONE_",
                "OSR_NONE_NONE_NONE"
            )
        )
}

association_file_more <- import_association_file(path_af, root_correct)
association_file_np <- af_no_prob(association_file_more)
association_file_miss <- af_missing_info(association_file_more)

# Matrices
import_matr_silent <- function(type) {
    op <- options(ISAnalytics.verbose = FALSE)
    on.exit(options(op))
    matrices <- import_parallel_Vispa2Matrices_auto(
        association_file = path_af, root = root_correct,
        quantification_type = type,
        matrix_type = "annotated", workers = 2, patterns = NULL,
        matching_opt = "ANY"
    )
    ex <- rlang::expr(`$`(matrices, !!type))
    matr <- rlang::eval_tidy(ex)
    matr
}

seq_count_m <- import_matr_silent("seqCount")
fe_m <- import_matr_silent("fragmentEstimate")

#------------------------------------------------------------------------------#
# Tests .check_same_info
#------------------------------------------------------------------------------#
test_that(".check_same_info returns empty if no problems", {
    check <- .check_same_info(association_file_np, seq_count_m)
    expect_true(nrow(check) == 0)
})

test_that(".check_same_info returns non empty if more info", {
    check <- .check_same_info(association_file_more, seq_count_m)
    expect_true(nrow(check) == 1)
    expect_true(check$ProjectID == "PROJECT1101")
    expect_true(check$PoolID == "ABY-LR-PL4-POOL54")
})

#------------------------------------------------------------------------------#
# Tests .join_matrix_af
#------------------------------------------------------------------------------#
test_that(".join_matrix_af produces the right table", {
    joined_df <- .join_matrix_af(
        seq_count_m, association_file_np,
        "SequencingDate"
    )
    expect_true(!is.null(joined_df))
    expect_true(nrow(joined_df) == nrow(seq_count_m))
    expect_true(ncol(joined_df) == 12)
})

#------------------------------------------------------------------------------#
# Tests .identify_independent_samples
#------------------------------------------------------------------------------#
### OTHER VARS ###
joined_df <- .join_matrix_af(seq_count_m, association_file_np, "SequencingDate")

test_that(".identify_independent_samples splits joined_df", {
    splitted <- .identify_independent_samples(joined_df)
    expect_true(nrow(splitted$collisions) +
        nrow(splitted$non_collisions) == nrow(joined_df))
    expect_true(nrow(dplyr::inner_join(splitted[[1]], splitted[[2]],
        by = c(
            "chr",
            "integration_locus",
            "strand"
        )
    )) == 0)
})

#------------------------------------------------------------------------------#
# Tests .obtain_nested
#------------------------------------------------------------------------------#
test_that(".obtain_nested correctly returns nested table", {
    splitted <- .identify_independent_samples(joined_df)
    nested <- .obtain_nested(splitted$collisions)
    expect_true(ncol(nested) == 4)
    expect_equal(colnames(nested), c(
        "chr", "integration_locus", "strand",
        "data"
    ))
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
    result <- .discriminate_by_seqCount(df, 10)
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
    result <- .discriminate_by_seqCount(df, 10)
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
        reads_ratio = 10
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
        reads_ratio = 10
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
        reads_ratio = 10
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
        reads_ratio = 10
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
        tibble::add_column(tibble::as_tibble_col(list(smpl1, smpl2,
                                                      smpl3, smpl4),
            column_name = "data"
        ))
    t
}

ex_collisions <- example_collisions()

test_that(".coll_mapping produces expected result for example collisions", {
    result <- .coll_mapping(ex_collisions, "SequencingDate", 10)
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
        "SequencingDate", 10
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
## TEST INPUT
test_that("remove_collisions stops if x is not a named list", {
    expect_error({
        rc <- remove_collisions(1, association_file_np)
    })
    expect_error({
        rc <- remove_collisions(list(1, 2, 3), association_file_np)
    })
})

test_that("remove_collisions stops if names of x are not in quant types", {
    # All
    expect_error({
        rc <- remove_collisions(list(a = "ab", b = "cd"), association_file_np)
    })
    # Only some
    expect_error({
        rc <- remove_collisions(
            list(seqCount = "ab", b = "cd"),
            association_file_np
        )
    })
})

test_that("remove_collisions stops if elements are not IS matrices", {
    non_ism <- tibble::tibble(a = c(1,2,3))
    miss_value <- seq_count_m %>% dplyr::select(-c(.data$Value))
    miss_camp <- seq_count_m %>%
        dplyr::select(-c(.data$CompleteAmplificationID))
    expect_error(
        {
            rc <- remove_collisions(x = non_ism, association_file_np)
        }, regexp = .non_ISM_error()
    )
    expect_error(
        {
            rc <- remove_collisions(x = miss_value, association_file_np)
        }, regexp = .missing_value_col_error()
    )
    expect_error(
        {
            rc <- remove_collisions(x = miss_camp, association_file_np)
        }, regexp = .missing_complAmpID_error()
    )
    expect_error(
        {
            rc <- remove_collisions(x = list(seqCount = non_ism,
                                             fragmentEstimate = seq_count_m),
                                    association_file_np)
        }, regexp = .non_ISM_error()
    )
    expect_error(
        {
            rc <- remove_collisions(x = list(seqCount = miss_value,
                                             fragmentEstimate = seq_count_m),
                                    association_file_np)
        }, regexp = .missing_value_col_error()
    )
    expect_error(
        {
            rc <- remove_collisions(x = list(seqCount = miss_camp,
                                             fragmentEstimate = seq_count_m),
                                    association_file_np)
        }, regexp = .missing_complAmpID_error()
    )
})

test_that("remove_collisions stops if association file is not a tibble", {
    expect_error({
        rc <- remove_collisions(
            list(seqCount = seq_count_m),
            path_af
        )
    })
})

test_that("remove_collisions stops if date_col is incorrect", {
    expect_error({
        rc <- remove_collisions(list(seqCount = seq_count_m),
            association_file_np,
            date_col = 2
        )
    })
    expect_error({
        rc <- remove_collisions(list(seqCount = seq_count_m),
            association_file_np,
            date_col = c(
                "SequencingDate", "FusionPrimerPCRDate"
            )
        )
    })
    expect_error({
        rc <- remove_collisions(list(seqCount = seq_count_m),
            association_file_np,
            date_col = "Date"
        )
    })
})

test_that("remove_collisions stops if reads_ratio is not a number", {
    expect_error({
        rc <- remove_collisions(list(seqCount = seq_count_m),
            association_file_np,
            date_col = "SequencingDate",
            reads_ratio = "20"
        )
    })
})

test_that("remove_collisions stops if reads_ratio is not a single number", {
    expect_error({
        rc <- remove_collisions(list(seqCount = seq_count_m),
            association_file_np,
            date_col = "SequencingDate",
            reads_ratio = c(10, 20)
        )
    })
})

test_that("remove_collisions stops if af is malformed", {
    af <- association_file_np %>% dplyr::select(-c(.data$ProjectID))
    expect_error(
        {
            rc <- remove_collisions(
                list(seqCount = seq_count_m),
                af
            )
        },
        regexp = "Malformed association file: one or more columns are missing"
    )
})

test_that("remove_collisions stops if date_col contains NA", {
    expect_error(
        {
            rc <- remove_collisions(list(seqCount = seq_count_m),
                association_file_np,
                date_col = "FusionPrimerPCRDate"
            )
        },
        regexp = paste(
            "Selected date column contains NA values, please check",
            "and fix the association file"
        )
    )
})

test_that("remove_collisions stops if there's no seqCount matrix", {
    expect_error(
        {
            rc <- remove_collisions(
                list(fragmentEstimate = seq_count_m),
                association_file_np
            )
        },
        regexp = paste(
            "Sequence count data frame is required for collision",
            "removal but none was detected in x"
        )
    )
    expect_error(
        {
            sq <- seq_count_m %>% dplyr::filter(.data$chr == "zzz")
            rc <- remove_collisions(
                list(seqCount = sq),
                association_file_np
            )
        },
        regexp = paste(
            "Sequence count data frame is required for collision",
            "removal but none was detected in x"
        )
    )
})

test_that("remove_collisions stops if missing info from af", {
    expect_error(
        {
            invisible(capture.output({
                rc <- remove_collisions(
                    list(seqCount = seq_count_m),
                    association_file_miss
                )
            }))
        },
        regexp = "The association file is missing needed info on some experiments"
    )
})

test_that("remove_collisions notifies additional data and succeeds", {
    op <- options(ISAnalytics.verbose = TRUE)
    on.exit(options(op))
    expect_message(
        {
            invisible(capture_output({
                coll_rem <- remove_collisions(seq_count_m, association_file_more,
                                              date_col = "SequencingDate",
                                              reads_ratio = 10
                )
            }))
        },
        regexp = paste("Found additional data relative to some projects",
                       "that are not included in the imported matrices.",
                       "Here is a summary",
                       collapse = "\n"
        ), fixed = TRUE
    )
})

#------------------------------------------------------------------------------#
# Tests realign_after_collisions
#------------------------------------------------------------------------------#

### OTHER VARS ###
silent_coll <- function() {
    op <- options(ISAnalytics.verbose = FALSE)
    on.exit(options(op))
    invisible(capture_output({
        coll_rem <- remove_collisions(seq_count_m, association_file_more,
            date_col = "SequencingDate",
            reads_ratio = 10
        )
    }))
    coll_rem
}

coll_rem <- silent_coll()

test_that("realign_after_collisions fails if others is not a named list", {
    expect_error({
        realigned <- realign_after_collisions(coll_rem$seqCount, 2)
    })
    expect_error({
        realigned <- realign_after_collisions(coll_rem$seqCount, NULL)
    })
    expect_error({
        realigned <- realign_after_collisions(
            coll_rem$seqCount,
            coll_rem$fragmentEstimate
        )
    })
})

test_that("realign_after_collisions fails if names not quant types", {
    expect_error({
        realigned <- realign_after_collisions(coll_rem$seqCount, fe_m)
    })
    expect_error({
        realigned <- realign_after_collisions(coll_rem$seqCount, list(a = fe_m))
    })
})

test_that("realign_after_collisions fails if list contains non-IS matrices", {
    fe1 <- fe_m %>% dplyr::select(-c(.data$CompleteAmplificationID))
    expect_error({
        realigned <- realign_after_collisions(
            coll_rem$seqCount,
            list(
                fragmentEstimate =
                    tibble::tibble(a = c(1,2,3))
            )
        )
    }, regexp = .non_ISM_error())
    expect_error({
        realigned <- realign_after_collisions(
            coll_rem$seqCount,
            list(
                fragmentEstimate = fe_m,
                barcodeCount = fe1
            )
        )
    }, regexp = .missing_complAmpID_error())
})

test_that("realign_after_collisions correctly re-aligns", {
    realigned <- realign_after_collisions(
        coll_rem$seqCount,
        list(fragmentEstimate = fe_m)
    )
    expect_true(nrow(coll_rem$seqCount) == nrow(realigned$fragmentEstimate))
    expect_true(all(realigned$fragmentEstimate$chr %in% coll_rem$seqCount$chr))
    expect_true(all(realigned$fragmentEstimate$integration_locus
        %in% coll_rem$seqCount$integration_locus))
    expect_true(all(realigned$fragmentEstimate$CompleteAmplificationID
        %in% coll_rem$seqCount$CompleteAmplificationID))
})
