library(ISAnalytics)
library(data.table, include.only = "%like%")
withr::local_options(list(ISAnalytics.verbose = FALSE))
#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
withr::local_options(ISAnalytics.verbose = FALSE)
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
    expect_equal(counts, c(6, 9, 5)) # S2 has duplicated is that is removed --> 9
})

#------------------------------------------------------------------------------#
# Test .find_in_common
#------------------------------------------------------------------------------#
test_that(".find_in_common works as expected", {
    key <- c("SubjectID", "CellMarker")
    lu <- .sh_obtain_lookup(key, test_sharing_input)
    labels <- data.table::data.table(g1 = "S1_CD34", g2 = "S2_CD14")
    common_is <- .find_in_common(labels,
        lookup_tbl = lu,
        keep_genomic_coord = FALSE
    )
    expect_equal(common_is, list(1))
    common_is <- .find_in_common(labels,
        lookup_tbl = lu,
        keep_genomic_coord = TRUE
    )
    is_1 <- lu[group_id == "S1_CD34"]$is[[1]]
    is_2 <- lu[group_id == "S2_CD14"]$is[[1]]
    expect_equal(common_is, list(
        list(1),
        list(is_2[is_1,
            on = mandatory_IS_vars(),
            nomatch = 0
        ])
    ))
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
# Test .sh_row_permut
#------------------------------------------------------------------------------#
test_that(".sh_row_permut works as expected", {
    test_row <- data.table::data.table(
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
    expect_true(all(perm[g1 == "S1_CD34", count_g1] == 3))
    expect_true(all(perm[g1 == "S2_CD14", count_g1] == 4))
    expect_true(all(perm[g1 == "S2_CD13", count_g1] == 3))
    expect_true(all(perm[g2 == "S1_CD34", count_g2] == 3))
    expect_true(all(perm[g2 == "S2_CD14", count_g2] == 4))
    expect_true(all(perm[g2 == "S2_CD13", count_g2] == 3))
    expect_true(all(perm[g3 == "S1_CD34", count_g3] == 3))
    expect_true(all(perm[g3 == "S2_CD14", count_g3] == 4))
    expect_true(all(perm[g3 == "S2_CD13", count_g3] == 3))
    # No counts
    test_row <- data.table::data.table(
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
    test_row <- data.table::data.table(
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
# Test .sharing_singledf_single_key
#------------------------------------------------------------------------------#
test_that(".sharing_singledf_single_key ok", {
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

#------------------------------------------------------------------------------#
# Test .sharing_singledf_mult_key
#------------------------------------------------------------------------------#
test_that(".sharing_singledf_mult_key ok", {
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
    expect_true(nrow(sh[g1 %like% "^S(1|2|3){1}_CD(34|14|13)$" &
        g2 %like% "^S(1|2|3){1}$"]) == 27)
    sh <- .sharing_singledf_mult_key(
        df = test_sharing_input,
        keys = keys, minimal = FALSE,
        is_count = TRUE,
        rel_sharing = TRUE,
        keep_genomic_coord = FALSE,
        venn = FALSE
    )
    expect_true(nrow(sh) == 54 & ncol(sh) == 9)
})

#------------------------------------------------------------------------------#
# Test .sharing_multdf_single_key
#------------------------------------------------------------------------------#
test_that(".sharing_multdf_single_key ok", {
    df1 <- test_sharing_input %>% dplyr::filter(CellMarker == "CD34")
    df2 <- test_sharing_input %>% dplyr::filter(CellMarker == "CD13")
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

#------------------------------------------------------------------------------#
# Test .sharing_multdf_mult_key
#------------------------------------------------------------------------------#
test_that(".sharing_multdf_mult_key ok", {
    df1 <- test_sharing_input %>% dplyr::filter(CellMarker == "CD34")
    df2 <- test_sharing_input %>% dplyr::filter(CellMarker == "CD13")
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
    expect_true(nrow(sh[g1 %like% "^S(1|2|3){1}$" &
        g2 %like% "^S(1|2|3){1}_CD(13)$"]) == 9)
    sh <- .sharing_multdf_mult_key(
        dfs = list(df1, df2),
        keys = keys,
        minimal = FALSE,
        is_count = TRUE,
        rel_sharing = TRUE,
        keep_genomic_coord = FALSE,
        venn = FALSE
    )
    expect_true(nrow(sh) == 18)
    expect_true(nrow(sh[g1 %like% "^S(1|2|3){1}$" &
        g2 %like% "^S(1|2|3){1}_CD(13)$"]) == 9)
    expect_true(nrow(sh[g2 %like% "^S(1|2|3){1}$" &
        g1 %like% "^S(1|2|3){1}_CD(13)$"]) == 9)
})

#------------------------------------------------------------------------------#
# Test is_sharing
#------------------------------------------------------------------------------#
test_that("is_sharing detects single key in list", {
    keys <- list(c("SubjectID", "CellMarker"))
    withr::with_options(list(ISAnalytics.verbose = TRUE), {
      expect_message(
      expect_message(
        expect_message(
            {
                sharing <- is_sharing(test_sharing_input, group_keys = keys)
            },
            class = "one_key_list"
        )))
    })
    keys <- list(
        c("SubjectID", "CellMarker"), c("SubjectID", "CellMarker"),
        c("SubjectID", "CellMarker")
    )
    withr::with_options(list(ISAnalytics.verbose = TRUE), {
      expect_message(
      expect_message(
        expect_message(
            {
                sharing <- is_sharing(test_sharing_input, group_keys = keys)
            },
            class = "one_key_list"
        )))
    })
})
