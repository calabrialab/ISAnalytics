
func_name <- "compute_abundance"
#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
test_df_simple <- tibble::tibble(
    chr = c("18", "8", "15", "15", "9", "15", "17", "7"),
    integration_locus = c(2247, 2524, 3917, 7227, 4499, 2200, 4305, 3792),
    strand = c("+", "+", "-", "+", "-", "+", "+", "-"),
    CompleteAmplificationID = c(
        "ID5", "ID6", "ID7", "ID5", "ID9",
        "ID10", "ID1", "ID6"
    ),
    Value = c(8, 51, 379, 450, 118, 294, 269, 402)
)

expected_totals_simple <- tibble::tibble(
    CompleteAmplificationID = c("ID5", "ID6", "ID7", "ID9", "ID10", "ID1"),
    Value_tot = c(458, 453, 379, 118, 294, 269)
)

expected_ab_simple <- test_df_simple %>%
    tibble::add_column(
        Value_tot = c(458, 453, 379, 458, 118, 294, 269, 453),
        Value_RelAbundance = c(
            0.017467249,
            0.112582781,
            1,
            0.982532751,
            1,
            1,
            1,
            0.887417219
        ),
        Value_PercAbundance = c(
            1.746724891,
            11.25827815,
            100,
            98.25327511,
            100,
            100,
            100,
            88.74172185
        )
    )

test_df_agg <- tibble::tribble(
    ~chr, ~integration_locus, ~strand, ~Key1, ~Key2, ~Val1, ~Val2, ~Val3,
    "14", 4999, "+", "KEY10", "key14", 176, 14, 312,
    "15", 4363, "+", "KEY15", "key15", 170, 42, 180,
    "7", 7079, "-", "KEY7", "key7", 475, 465, 43,
    "21", 6962, "+", "KEY21", "key21", 452, 83, 337,
    "21", 5329, "-", "KEY21", "key23", 147, 143, 225,
    "1", 2804, "+", "KEY10", "key1", 9, 203, 137,
    "10", 7471, "+", "KEY10", "key10", 116, 348, 442,
    "2", 2955, "-", "KEY2", "key21", 202, 117, 408
)

expected_ab_agg <- tibble::tribble(
    ~chr, ~integration_locus, ~strand, ~Key1, ~Key2, ~Val1, ~Val2, ~Val3,
    ~Val1_tot, ~Val2_tot, ~Val3_tot, ~Val1_RelAbundance, ~Val2_RelAbundance,
    ~Val3_RelAbundance, ~Val1_PercAbundance, ~Val2_PercAbundance,
    ~Val3_PercAbundance,
    "14", 4999, "+", "KEY10", "key14", 176, 14, 312, 301, 565, 891,
    0.584717608, 0.024778761, 0.35016835, 58.4717608, 2.477876106, 35.01683502,
    "15", 4363, "+", "KEY15", "key15", 170, 42, 180, 170, 42, 180, 1, 1, 1,
    100, 100, 100,
    "7", 7079, "-", "KEY7", "key7", 475, 465, 43, 475, 465, 43, 1, 1, 1,
    100, 100, 100,
    "21", 6962, "+", "KEY21", "key21", 452, 83, 337, 599, 226, 562, 0.754590985,
    0.367256637, 0.599644128, 75.4590985, 36.72566372, 59.96441281,
    "21", 5329, "-", "KEY21", "key23", 147, 143, 225, 599, 226, 562, 0.245409015,
    0.632743363, 0.400355872, 24.5409015, 63.27433628, 40.03558719,
    "1", 2804, "+", "KEY10", "key1", 9, 203, 137, 301, 565, 891, 0.029900332,
    0.359292035, 0.15375982, 2.990033223, 35.92920354, 15.37598204,
    "10", 7471, "+", "KEY10", "key10", 116, 348, 442, 301, 565, 891, 0.38538206,
    0.615929204, 0.496071829, 38.53820598, 61.59292035, 49.60718294,
    "2", 2955, "-", "KEY2", "key21", 202, 117, 408, 202, 117, 408, 1, 1, 1,
    100, 100, 100
)

expected_totals_agg <- tibble::tribble(
    ~Key1, ~Val1_tot, ~Val2_tot, ~Val3_tot,
    "KEY10", 301, 565, 891,
    "KEY21", 599, 226, 562,
    "KEY15", 170, 42, 180,
    "KEY7", 475, 465, 43,
    "KEY2", 202, 117, 408
)
#------------------------------------------------------------------------------#
# Test compute_abundance
#------------------------------------------------------------------------------#
### Test on inputs
### - x
test_that(paste(func_name, "fails if x is not a data frame"), {
    expect_error(compute_abundance(x = list(a = c(1, 2, 3))))
})

### - Columns
test_that(paste(func_name, "fails if columns is not chr"), {
    expect_error(compute_abundance(
        x = test_df_simple, columns = c(1, 2),
        key = "CompleteAmplificationID"
    ))
    expect_error(compute_abundance(
        x = test_df_simple, columns = NULL,
        key = "CompleteAmplificationID"
    ))
})
test_that(paste(func_name, "fails if columns is empty"), {
    expect_error(compute_abundance(
        x = test_df_simple,
        columns = "",
        key = "CompleteAmplificationID"
    ))
})
test_that(paste(func_name, "fails if columns not found"), {
    err <- expect_error(compute_abundance(
        x = test_df_simple,
        columns = c("a", "b"),
        key = "CompleteAmplificationID"
    ))
    expect_equal(
        err$message,
        stringr::str_replace(
            rlang::format_error_bullets(.missing_user_cols_error(c("a", "b"))),
            "\\* ", ""
        )
    )
})
test_that(paste(func_name, "fails if any columns not numeric"), {
    expect_error(compute_abundance(
        x = test_df_simple, columns = c("chr"),
        key = "CompleteAmplificationID"
    ))
})

### - Percentage
test_that(paste(func_name, "fails if percentage is not logical"), {
    expect_error(compute_abundance(
        x = test_df_simple,
        columns = "Value",
        percentage = "true",
        key = "CompleteAmplificationID"
    ))
})

### - Key
test_that(paste(func_name, "fails if key is not chr"), {
    expect_error(compute_abundance(
        x = test_df_simple,
        key = 1,
        columns = "Value",
    ))
    expect_error(compute_abundance(
        x = test_df_simple,
        key = NULL,
        columns = "Value"
    ))
})
test_that(paste(func_name, "fails if key is empty"), {
    expect_error(compute_abundance(
        x = test_df_simple,
        key = "",
        columns = "Value"
    ))
})
test_that(paste(func_name, "fails if key not found"), {
    err <- expect_error(compute_abundance(
        x = test_df_simple,
        key = "a",
        columns = "Value"
    ))
    expect_equal(
        err$message,
        stringr::str_replace(
            rlang::format_error_bullets(.missing_user_cols_error(c("a"))),
            "\\* ", ""
        )
    )
})

### - Keep_totals
test_that(paste(func_name, "fails if keep_totals is not one of values"), {
    expect_error(compute_abundance(
        x = test_df_simple,
        keep_totals = "true",
        columns = "Value",
        key = "CompleteAmplificationID"
    ))
    expect_error(compute_abundance(
        x = test_df_simple,
        keep_totals = NULL,
        columns = "Value",
        key = "CompleteAmplificationID"
    ))
})

# Output correctness
### For simple non-agg matrix
test_that(paste(func_name, "has correct output for simple is matrix"), {
    abund <- compute_abundance(test_df_simple,
        columns = "Value",
        key = "CompleteAmplificationID"
    )
    expect_equal(abund, expected_ab_simple %>% dplyr::select(-Value_tot))
    abund <- compute_abundance(test_df_simple,
        keep_totals = TRUE,
        columns = "Value",
        key = "CompleteAmplificationID"
    )
    expect_equal(abund, expected_ab_simple)
    abund <- compute_abundance(test_df_simple,
        keep_totals = "df",
        columns = "Value",
        key = "CompleteAmplificationID"
    )
    expect_equal(abund$abundance_df, expected_ab_simple
    %>% dplyr::select(-Value_tot))
    expect_equal(
        abund$quant_totals %>%
            dplyr::arrange(CompleteAmplificationID),
        expected_totals_simple %>%
            dplyr::arrange(CompleteAmplificationID)
    )
    abund <- compute_abundance(test_df_simple,
        percentage = FALSE,
        columns = "Value",
        key = "CompleteAmplificationID"
    )
    expect_equal(abund, expected_ab_simple
    %>% dplyr::select(-Value_tot, -Value_PercAbundance))
})

## For agg matrix
test_that(paste(func_name, "has correct output for agg is matrix"), {
    abund <- compute_abundance(test_df_agg,
        columns = c("Val1", "Val2", "Val3"),
        key = "Key1", keep_totals = TRUE
    )
    expect_equal(abund, expected_ab_agg)
    abund <- compute_abundance(test_df_agg,
        columns = c("Val1", "Val2", "Val3"),
        key = "Key1",
        keep_totals = FALSE
    )
    expect_equal(abund, expected_ab_agg %>% dplyr::select(
        -Val1_tot,
        -Val2_tot,
        -Val3_tot
    ))
    abund <- compute_abundance(test_df_agg,
        columns = c("Val1", "Val2", "Val3"),
        key = "Key1",
        keep_totals = "df"
    )
    expect_equal(abund$abundance_df, expected_ab_agg %>%
        dplyr::select(
            -Val1_tot,
            -Val2_tot,
            -Val3_tot
        ))
    expect_equal(
        abund$quant_totals %>%
            dplyr::arrange(Key1),
        expected_totals_agg %>%
            dplyr::arrange(Key1)
    )
    abund <- compute_abundance(test_df_agg,
        columns = c("Val1", "Val2", "Val3"),
        key = "Key1",
        percentage = FALSE
    )
    expect_equal(abund, expected_ab_agg %>% dplyr::select(
        -Val1_tot,
        -Val2_tot,
        -Val3_tot,
        -Val1_PercAbundance,
        -Val2_PercAbundance,
        -Val3_PercAbundance,
    ))
})
