withr::local_options(list(ISAnalytics.verbose = FALSE))
#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
test_sharing_input <- tibble::tribble(
    ~chr, ~integration_locus, ~strand, ~SubjectID, ~CellMarker, ~Tissue,
    "1", 10032, "+", "S1", "CD34", "BM",
    "10", 20162, "-", "S1", "CD34", "BM",
    "5", 45612, "-", "S1", "CD34", "BM",
    "2", 21206, "+", "S1", "CD13", "PB",
    "1", 10102, "+", "S1", "CD13", "PB",
    "3", 38167, "-", "S1", "CD14", "PB",
    "1", 10032, "+", "S2", "CD14", "PB",
    "5", 45612, "-", "S2", "CD34", "BM",
    "6", 12542, "-", "S2", "CD13", "PB",
    "1", 42571, "-", "S2", "CD14", "PB",
    "4", 21054, "+", "S2", "CD14", "PB",
    "4", 12072, "-", "S2", "CD13", "PB",
    "11", 25722, "-", "S2", "CD34", "BM",
    "10", 11725, "+", "S2", "CD14", "PB",
    "10", 12247, "+", "S2", "CD34", "BM",
    "11", 25722, "-", "S2", "CD13", "PB",
    "1", 10032, "+", "S3", "CD14", "PB",
    "2", 21206, "+", "S3", "CD14", "PB",
    "6", 12542, "-", "S3", "CD14", "PB",
    "1", 42571, "-", "S3", "CD34", "BM",
    "2", 51232, "+", "S3", "CD13", "PB"
)

#------------------------------------------------------------------------------#
# Test .sh_obtain_lookup
#------------------------------------------------------------------------------#
test_that(".sh_obtain_lookup correct output", {
    key <- c("SubjectID", "Tissue", "CellMarker")
    lu <- .sh_obtain_lookup(key, test_sharing_input)
    expect_true(length(lu$group_id) == 9)
    counts <- purrr::map_int(lu$is, ~ nrow(.x))
    expect_equal(counts, c(3, 2, 1, 4, 3, 3, 3, 1, 1))

    # Testing duplicated is
    key <- c("SubjectID")
    lu <- .sh_obtain_lookup(key, test_sharing_input)
    expect_true(length(lu$group_id) == 3)
    counts <- purrr::map_int(lu$is, ~ nrow(.x))
    expect_equal(
        counts, c(6, 9, 5)
    ) # S2 has duplicated is that is removed --> 9
})

#------------------------------------------------------------------------------#
# Test .find_in_common
#------------------------------------------------------------------------------#
test_that(".find_in_common works as expected", {
    key <- c("SubjectID", "CellMarker")
    lu <- .sh_obtain_lookup(key, test_sharing_input)
    labels <- tibble::tibble(g1 = "S1_CD34", g2 = "S2_CD14")
    common_is <- purrr::pmap(labels, .find_in_common,
        lookup_tbl = lu,
        keep_genomic_coord = FALSE
    )
    expect_equal(common_is[[1]], tibble::tibble(shared = 1))
    common_is <- purrr::pmap(labels, .find_in_common,
        lookup_tbl = lu,
        keep_genomic_coord = TRUE
    )
    is_1 <- (lu |>
        dplyr::filter(.data$group_id == "S1_CD34") |>
        dplyr::pull("is"))[[1]]
    is_2 <- (lu |>
        dplyr::filter(.data$group_id == "S2_CD14") |>
        dplyr::pull("is"))[[1]]
    expect_equal(
        common_is[[1]],
        tibble::tibble(
            shared = 1,
            is_coord = list(
                is_1 |>
                    dplyr::inner_join(is_2, by = mandatory_IS_vars())
            )
        )
    )
})

#------------------------------------------------------------------------------#
# Test .count_group_union
#------------------------------------------------------------------------------#
test_that(".count_group_union works as expected", {
    key <- c("SubjectID", "CellMarker")
    lu <- .sh_obtain_lookup(key, test_sharing_input)
    labels <- c("S1_CD34", "S2_CD14")
    c_union <- .count_group_union(
        g1 = labels[1], g2 = labels[2],
        col_groups = c("g1", "g2"), lookup_tbl = lu
    )
    expect_equal(c_union$count_union, 6)
})

#------------------------------------------------------------------------------#
# Test .sh_truth_tbl_venn
#------------------------------------------------------------------------------#
test_that(".sh_truth_tbl_venn works as expected", {
    key <- c("SubjectID", "CellMarker")
    lu <- .sh_obtain_lookup(key, test_sharing_input)
    truth_tbl <- .sh_truth_tbl_venn(
        g1 = "S2_CD13", g2 = "S2_CD34", lookup = lu,
        groups = c("g1", "g2")
    )
    expect_true(
        all(truth_tbl |>
            dplyr::filter(.data$S2_CD34 == TRUE & .data$S2_CD13 == TRUE) |>
            dplyr::pull(.data$int_id) == "11_25722_-")
    )
})

#------------------------------------------------------------------------------#
# Test .sh_row_permut
#------------------------------------------------------------------------------#
test_that(".sh_row_permut works as expected", {
    test_row <- list(
        g1 = "S2_CD13",
        g2 = "S2_CD34",
        shared = 1,
        is_coord = list(tibble::tibble(
            chr = "11",
            integration_locus = 25722,
            strand = "-"
        )),
        count_g1 = 3,
        count_g2 = 3,
        count_union = 5,
        truth_tbl_venn = list(tibble::tibble(
            int_id = c(
                "5_45612_-", "11_25722_-", "10_12247_+", "6_12542_-",
                "4_12072_-"
            ),
            S2_CD34 = c(TRUE, TRUE, TRUE, FALSE, FALSE),
            S2_CD13 = c(FALSE, TRUE, FALSE, TRUE, TRUE)
        ))
    )
    perm <- .sh_row_permut(
        !!!test_row,
        g_names = c("g1", "g2"),
        counts = TRUE
    )
    expect_equal(perm$truth_tbl_venn[[1]], test_row$truth_tbl_venn[[1]])
    expect_equal(perm$truth_tbl_venn[[2]], test_row$truth_tbl_venn[[1]])
    expect_equal(perm$is_coord[[1]], test_row$is_coord[[1]])
    expect_equal(perm$is_coord[[2]], test_row$is_coord[[1]])
    expect_equal(nrow(perm), 2)
})

test_that(".sh_row_permut does nothing for same ids", {
    test_row <- list(
        g1 = "S2_CD34",
        g2 = "S2_CD34",
        shared = 3,
        count_g1 = 3,
        count_g2 = 3,
        count_union = 3
    )
    perm <- .sh_row_permut(
        !!!test_row,
        g_names = c("g1", "g2"),
        counts = TRUE
    )
    expect_equal(nrow(perm), 1)
})

test_that(".sh_row_permut works as expected - 3 groups", {
    test_row <- tibble::tibble(
        g1 = "S1_CD34",
        g2 = "S2_CD14",
        g3 = "S2_CD13",
        shared = 0,
        count_g1 = 3,
        count_g2 = 4,
        count_g3 = 3,
        count_union = 9
    )
    perm <- .sh_row_permut(
        g1 = test_row$g1,
        g2 = test_row$g2,
        g3 = test_row$g3,
        shared = test_row$shared,
        count_g1 = test_row$count_g1,
        count_g2 = test_row$count_g2,
        count_g3 = test_row$count_g3,
        count_union = test_row$count_union,
        g_names = c("g1", "g2", "g3"),
        counts = TRUE
    )
    expect_true(nrow(perm) == 6)
    expect_true(all(perm$shared == 0))
    expect_true(all(perm$count_union == 9))
    expect_true(all(perm |>
        dplyr::filter(.data[["g1"]] == "S1_CD34") |>
        dplyr::pull("count_g1") == 3))
    expect_true(all(perm |>
        dplyr::filter(.data[["g1"]] == "S2_CD14") |>
        dplyr::pull("count_g1") == 4))
    expect_true(all(perm |>
        dplyr::filter(.data[["g1"]] == "S2_CD13") |>
        dplyr::pull("count_g1") == 3))
    expect_true(all(perm |>
        dplyr::filter(.data[["g2"]] == "S1_CD34") |>
        dplyr::pull("count_g2") == 3))
    expect_true(all(perm |>
        dplyr::filter(.data[["g2"]] == "S2_CD14") |>
        dplyr::pull("count_g2") == 4))
    expect_true(all(perm |>
        dplyr::filter(.data[["g2"]] == "S2_CD13") |>
        dplyr::pull("count_g2") == 3))
    expect_true(all(perm |>
        dplyr::filter(.data[["g3"]] == "S1_CD34") |>
        dplyr::pull("count_g3") == 3))
    expect_true(all(perm |>
        dplyr::filter(.data[["g3"]] == "S2_CD14") |>
        dplyr::pull("count_g3") == 4))
    expect_true(all(perm |>
        dplyr::filter(.data[["g3"]] == "S2_CD13") |>
        dplyr::pull("count_g3") == 3))
    # No counts
    test_row <- tibble::tibble(
        g1 = "S1_CD34",
        g2 = "S2_CD14",
        g3 = "S2_CD13",
        shared = 0
    )
    perm <- .sh_row_permut(
        g1 = test_row$g1,
        g2 = test_row$g2,
        g3 = test_row$g3,
        shared = test_row$shared,
        g_names = c("g1", "g2", "g3"),
        counts = FALSE
    )
    expect_true(nrow(perm) == 6)
    expect_true(all(perm$shared == 0))
    # Self row
    test_row <- tibble::tibble(
        g1 = "S1_CD34",
        g2 = "S1_CD34",
        g3 = "S1_CD34",
        shared = 3,
        count_g1 = 3,
        count_g2 = 3,
        count_g3 = 3,
        count_union = 3
    )
    perm <- .sh_row_permut(
        g1 = test_row$g1,
        g2 = test_row$g2,
        g3 = test_row$g3,
        shared = test_row$shared,
        count_g1 = test_row$count_g1,
        count_g2 = test_row$count_g2,
        count_g3 = test_row$count_g3,
        count_union = test_row$count_union,
        g_names = c("g1", "g2", "g3"),
        counts = TRUE
    )
    expect_true(nrow(perm) == 1)
})

#------------------------------------------------------------------------------#
# Test .single_row_sharing
#------------------------------------------------------------------------------#
test_that(".single_row_sharing works as expected", {
    test_row <- tibble::tibble(
        g1 = c("S1_CD13"),
        g2 = c("S1_CD14"),
        g3 = c("S1_CD34")
    )
    key <- c("SubjectID", "CellMarker")
    lu <- .sh_obtain_lookup(key, test_sharing_input)
    is_counts <- tibble::tibble(
        group_id = c(
            "S1_CD34", "S1_CD13", "S1_CD14", "S2_CD14", "S2_CD34",
            "S2_CD13", "S3_CD14", "S3_CD34", "S3_CD13"
        ),
        count = c(3, 2, 1, 4, 3, 3, 3, 1, 1)
    )
    result <- .single_row_sharing(
        test_row,
        lookup = lu,
        is_counts = is_counts,
        add_is_count = TRUE,
        keep_genomic_coord = TRUE,
        rel_sharing = TRUE, venn = TRUE,
        minimal = FALSE, progress = NULL
    )
    expected <- tibble::tibble(
        g1 = c("S1_CD13", "S1_CD13", "S1_CD14", "S1_CD14", "S1_CD34", "S1_CD34"),
        g2 = c("S1_CD14", "S1_CD34", "S1_CD13", "S1_CD34", "S1_CD13", "S1_CD14"),
        g3 = c("S1_CD34", "S1_CD14", "S1_CD34", "S1_CD13", "S1_CD14", "S1_CD13"),
        shared = c(0, 0, 0, 0, 0, 0),
        count_g1 = c(2, 2, 1, 1, 3, 3),
        count_g2 = c(1, 3, 2, 3, 2, 1),
        count_g3 = c(3, 1, 3, 2, 1, 2),
        count_union = c(6, 6, 6, 6, 6, 6),
        on_g1 = c(0, 0, 0, 0, 0, 0),
        on_g2 = c(0, 0, 0, 0, 0, 0),
        on_g3 = c(0, 0, 0, 0, 0, 0),
        on_union = c(0, 0, 0, 0, 0, 0)
    )
    expect_equal(
        result |>
            dplyr::select(-dplyr::all_of(c("is_coord", "truth_tbl_venn"))) |>
            dplyr::arrange(.data$g1),
        expected |>
            dplyr::arrange(.data$g1)
    )
})

#------------------------------------------------------------------------------#
# Test .sharing_singledf_single_key
#------------------------------------------------------------------------------#
test_that(".sharing_singledf_single_key ok", {
    withr::local_options(list(
        ISAnalytics.parallel_processing = FALSE
    ))
    key <- c("SubjectID", "CellMarker")
    sh <- .sharing_singledf_single_key(
        df = test_sharing_input,
        key = key,
        minimal = TRUE,
        n_comp = 2,
        is_count = FALSE,
        rel_sharing = FALSE,
        include_self_comp = FALSE,
        keep_genomic_coord = FALSE,
        venn = FALSE
    )
    expect_true(nrow(sh) == 36 & ncol(sh) == 3)
    sh <- .sharing_singledf_single_key(
        df = test_sharing_input,
        key = key,
        minimal = TRUE,
        n_comp = 2,
        is_count = FALSE,
        rel_sharing = FALSE,
        include_self_comp = TRUE,
        keep_genomic_coord = FALSE,
        venn = FALSE
    )
    expect_true(nrow(sh) == 36 + 9 & ncol(sh) == 3)
    sh <- .sharing_singledf_single_key(
        df = test_sharing_input,
        key = key,
        minimal = FALSE,
        n_comp = 2,
        is_count = FALSE,
        rel_sharing = FALSE,
        include_self_comp = TRUE,
        keep_genomic_coord = FALSE,
        venn = FALSE
    )
    expect_true(nrow(sh) == 36 * 2 + 9 & ncol(sh) == 3)
    sh <- .sharing_singledf_single_key(
        df = test_sharing_input,
        key = key,
        minimal = FALSE,
        n_comp = 2,
        is_count = TRUE,
        rel_sharing = FALSE,
        include_self_comp = TRUE,
        keep_genomic_coord = FALSE,
        venn = FALSE
    )
    expect_true(nrow(sh) == 36 * 2 + 9 & ncol(sh) == 6)
    sh <- .sharing_singledf_single_key(
        df = test_sharing_input,
        key = key,
        minimal = FALSE,
        n_comp = 2,
        is_count = TRUE,
        rel_sharing = TRUE,
        include_self_comp = TRUE,
        keep_genomic_coord = FALSE,
        venn = FALSE
    )
    expect_true(nrow(sh) == 36 * 2 + 9 & ncol(sh) == 9)
    sh <- .sharing_singledf_single_key(
        df = test_sharing_input,
        key = key,
        minimal = FALSE,
        n_comp = 2,
        is_count = TRUE,
        rel_sharing = TRUE,
        include_self_comp = TRUE,
        keep_genomic_coord = TRUE,
        venn = FALSE
    )
    expect_true(nrow(sh) == 36 * 2 + 9 & "is_coord" %in% colnames(sh))
    sh <- .sharing_singledf_single_key(
        df = test_sharing_input,
        key = key,
        minimal = FALSE,
        n_comp = 2,
        is_count = TRUE,
        rel_sharing = TRUE,
        include_self_comp = TRUE,
        keep_genomic_coord = TRUE,
        venn = TRUE
    )
    expect_true(nrow(sh) == 36 * 2 + 9 & "truth_tbl_venn" %in% colnames(sh))
})

test_that(".sharing_singledf_single_key works with 1 group", {
    withr::local_options(list(
        ISAnalytics.parallel_processing = FALSE
    ))
    sub <- test_sharing_input |>
        dplyr::filter(.data$SubjectID == "S1")
    key <- "SubjectID"
    sh <- .sharing_singledf_single_key(
        df = sub,
        key = key,
        minimal = FALSE,
        n_comp = 2,
        is_count = TRUE,
        rel_sharing = TRUE,
        include_self_comp = TRUE,
        keep_genomic_coord = TRUE,
        venn = TRUE
    )
    expect_true(nrow(sh) == 1)
    expect_equal(sh$g1, "S1")
    expect_equal(sh$g2, "S1")
    expect_equal(sh$shared, 6)
    expect_equal(sh$count_g1, 6)
    expect_equal(sh$count_g2, 6)
    expect_equal(sh$count_union, 6)
    expect_equal(sh$on_g1, 100)
    expect_equal(sh$on_g2, 100)
    expect_equal(sh$on_union, 100)

    sh_1 <- .sharing_singledf_single_key(
        df = sub,
        key = key,
        minimal = FALSE,
        n_comp = 2,
        is_count = TRUE,
        rel_sharing = TRUE,
        include_self_comp = FALSE,
        keep_genomic_coord = TRUE,
        venn = TRUE
    )
    expect_null(sh_1)
})

#------------------------------------------------------------------------------#
# Test .sharing_singledf_mult_key
#------------------------------------------------------------------------------#
test_that(".sharing_singledf_mult_key ok", {
    withr::local_options(list(
        ISAnalytics.parallel_processing = FALSE
    ))
    keys <- list(g1 = c("SubjectID", "CellMarker"), g2 = c("SubjectID"))
    sh <- .sharing_singledf_mult_key(
        df = test_sharing_input,
        keys = keys, minimal = TRUE,
        is_count = TRUE,
        rel_sharing = TRUE,
        keep_genomic_coord = FALSE,
        venn = FALSE
    )
    expect_true(nrow(sh) == 27 & ncol(sh) == 9)
    expect_true(all(
        stringr::str_detect(sh$g1, "^S(1|2|3)_CD[1-9]{2}$") &
            stringr::str_detect(sh$g2, "^S[1-3]{1}$")
    ))
    sh <- .sharing_singledf_mult_key(
        df = test_sharing_input,
        keys = keys, minimal = FALSE,
        is_count = TRUE,
        rel_sharing = TRUE,
        keep_genomic_coord = FALSE,
        venn = FALSE
    )
    expect_true(nrow(sh) == 54 & ncol(sh) == 9)
    sh <- .sharing_singledf_mult_key(
        df = test_sharing_input,
        keys = keys, minimal = FALSE,
        is_count = TRUE,
        rel_sharing = TRUE,
        keep_genomic_coord = TRUE,
        venn = TRUE
    )
    expect_true("is_coord" %in% colnames(sh))
    expect_true("truth_tbl_venn" %in% colnames(sh))
})

#------------------------------------------------------------------------------#
# Test .sharing_multdf_single_key
#------------------------------------------------------------------------------#
test_that(".sharing_multdf_single_key ok", {
    withr::local_options(list(
        ISAnalytics.parallel_processing = FALSE
    ))
    df1 <- test_sharing_input |> dplyr::filter(CellMarker == "CD34")
    df2 <- test_sharing_input |> dplyr::filter(CellMarker == "CD13")
    key <- c("SubjectID", "CellMarker")
    sh <- .sharing_multdf_single_key(
        dfs = list(df1, df2),
        key = key, minimal = TRUE,
        is_count = TRUE,
        rel_sharing = TRUE,
        keep_genomic_coord = FALSE,
        venn = FALSE
    )
    expect_true(nrow(sh) == 9)
    sh <- .sharing_multdf_single_key(
        dfs = list(df1, df2),
        key = key, minimal = FALSE,
        is_count = TRUE,
        rel_sharing = TRUE,
        keep_genomic_coord = FALSE,
        venn = FALSE
    )
    expect_true(nrow(sh) == 18)
})

test_that(".sharing_multdf_single_key handles limit cases", {
    withr::local_options(list(
        ISAnalytics.parallel_processing = FALSE
    ))
    df1 <- test_sharing_input |>
        dplyr::filter(CellMarker == "CD34", SubjectID == "S1") |>
        dplyr::filter(!.data$chr == "1")
    df2 <- test_sharing_input |>
        dplyr::filter(CellMarker == "CD34", SubjectID == "S1")
    key <- c("SubjectID", "CellMarker")
    sh <- .sharing_multdf_single_key(
        dfs = list(df1, df2),
        key = key, minimal = FALSE,
        is_count = TRUE,
        rel_sharing = TRUE,
        keep_genomic_coord = TRUE,
        venn = TRUE
    )
    expect_true(all(sh$shared == 2))
    expect_true(nrow(sh) == 2)
    expect_true(sh[1, ]$count_g1 == sh[2, ]$count_g2 &
        sh[2, ]$count_g1 == sh[1, ]$count_g2)
    expect_equal(sh[1, ]$truth_tbl_venn[[1]], sh[2, ]$truth_tbl_venn[[1]])
    expect_equal(sh[1, ]$is_coord[[1]], sh[2, ]$is_coord[[1]])
})

#------------------------------------------------------------------------------#
# Test .sharing_multdf_mult_key
#------------------------------------------------------------------------------#
test_that(".sharing_multdf_mult_key ok", {
    withr::local_options(list(
        ISAnalytics.parallel_processing = FALSE
    ))
    df1 <- test_sharing_input |> dplyr::filter(CellMarker == "CD34")
    df2 <- test_sharing_input |> dplyr::filter(CellMarker == "CD13")
    keys <- list(g1 = c("SubjectID"), g2 = c("SubjectID", "CellMarker"))
    sh <- .sharing_multdf_mult_key(
        dfs = list(df1, df2),
        keys = keys,
        minimal = TRUE,
        is_count = TRUE,
        rel_sharing = TRUE,
        keep_genomic_coord = FALSE,
        venn = FALSE
    )
    expect_true(all(
        stringr::str_detect(sh$g2, "^S(1|2|3)_CD13-2$") &
            stringr::str_detect(sh$g1, "^S[1-3]{1}-1$")
    ))
    expect_true(nrow(sh) == 9)
    expect_true(
        sh |>
            dplyr::filter(.data$g1 == "S2-1", .data$g2 == "S2_CD13-2") |>
            dplyr::pull("shared") == 1
    )

    sh <- .sharing_multdf_mult_key(
        dfs = list(df1, df2),
        keys = keys,
        minimal = FALSE,
        is_count = TRUE,
        rel_sharing = TRUE,
        keep_genomic_coord = TRUE,
        venn = TRUE
    )
    expect_true(nrow(sh) == 18)
})

#------------------------------------------------------------------------------#
# Test is_sharing
#------------------------------------------------------------------------------#
test_that("is_sharing detects single key in list", {
    keys <- list(c("SubjectID", "CellMarker"))
    withr::local_options(list(
        ISAnalytics.parallel_processing = FALSE
    ))
    withr::with_options(list(ISAnalytics.verbose = TRUE), {
        expect_message(
            {
                sharing <- is_sharing(test_sharing_input,
                    group_keys = keys
                )
            },
            class = "one_key_list"
        )
    })
    keys <- list(
        c("SubjectID", "CellMarker"), c("SubjectID", "CellMarker"),
        c("SubjectID", "CellMarker")
    )
    withr::with_options(list(ISAnalytics.verbose = TRUE), {
        expect_message(
            {
                sharing <- is_sharing(test_sharing_input,
                    group_keys = keys
                )
            },
            class = "one_key_list"
        )
    })
})
