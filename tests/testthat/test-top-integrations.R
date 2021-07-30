library(ISAnalytics)
func_name <- "top_integrations"
#------------------------------------------------------------------------------#
# Global vars
#------------------------------------------------------------------------------#
test_df <- tibble::tibble(
    chr = c("1", "2", "3", "4", "5", "6"),
    integration_locus = c(14536, 14544, 14512, 14236, 14522, 14566),
    strand = c("+", "+", "-", "+", "-", "+"),
    CompleteAmplificationID = c("ID1", "ID2", "ID1", "ID1", "ID3", "ID2"),
    Value = c(3, 10, 40, 2, 15, 150),
    Value2 = c(456, 87, 87, 9, 64, 96),
    Value3 = c("a", "b", "c", "d", "e", "f")
)

expected_no_key <- tibble::tibble(
    chr = c("6", "3", "5"),
    integration_locus = c(14566, 14512, 14522),
    strand = c("+", "-", "-"),
    Value = c(150, 40, 15),
    CompleteAmplificationID = c("ID2", "ID1", "ID3"),
    Value2 = c(96, 87, 64),
    Value3 = c("f", "c", "e")
)

expected_no_key2 <- tibble::tibble(
    chr = c("6", "3", "5"),
    integration_locus = c(14566, 14512, 14522),
    strand = c("+", "-", "-"),
    Value = c(150, 40, 15),
    Value2 = c(96, 87, 64),
    CompleteAmplificationID = c("ID2", "ID1", "ID3"),
    Value3 = c("f", "c", "e")
)

expected_key <- tibble::tibble(
    CompleteAmplificationID = c("ID1", "ID2", "ID3"),
    chr = c("3", "6", "5"),
    integration_locus = c(14512, 14566, 14522),
    strand = c("-", "+", "-"),
    Value = c(40, 150, 15),
    Value2 = c(87, 96, 64),
    Value3 = c("c", "f", "e")
)

#------------------------------------------------------------------------------#
# Test top_integrations
#------------------------------------------------------------------------------#
### Test on inputs
test_that(paste0(func_name, "stops if params incorrect"), {
    ## Missing value cols
    err <- expect_error({
        top_integrations(test_df)
    })
    expect_equal(
        err$message,
        stringr::str_replace(
            rlang::format_error_bullets(
                .missing_user_cols_error("fragmentEstimate_sum_RelAbundance")
            ),
            "\\* ", ""
        )
    )
    ## Missing cols to keep
    err <- expect_error({
        top_integrations(test_df, columns = "Value", keep = "a")
    })
    expect_equal(
        err$message,
        stringr::str_replace(
            rlang::format_error_bullets(
                .missing_user_cols_error("a")
            ),
            "\\* ", ""
        )
    )
    ## No mand vars
    err <- expect_error({
        top_integrations(test_df[2:7], columns = "Value")
    })
    ## Missing non null key
    err <- expect_error({
        top_integrations(test_df, columns = "Value", key = "a")
    })
    expect_equal(
        err$message,
        stringr::str_replace(
            rlang::format_error_bullets(
                .missing_user_cols_error("a")
            ),
            "\\* ", ""
        )
    )
})

### Test on outputs
test_that(paste0(func_name, "works on df no key"), {
    top <- top_integrations(test_df, n = 3, columns = "Value")
    expect_equal(top, expected_no_key)
})

test_that(paste0(func_name, "works on df no key secondary ordering"), {
    top <- top_integrations(test_df, n = 3, columns = c("Value", "Value2"))
    expect_equal(top, expected_no_key2)
})

test_that(paste0(func_name, "works on df with key"), {
    top <- top_integrations(test_df,
        n = 1, columns = "Value",
        key = "CompleteAmplificationID"
    )
    expect_equal(top, expected_key)
})

test_that(paste0(func_name, "drops unwanted columns"), {
    top <- top_integrations(test_df, n = 3, columns = "Value", keep = "nothing")
    expect_equal(top, expected_no_key %>%
        dplyr::select(chr, integration_locus, strand, Value))
    top <- top_integrations(test_df, n = 3, columns = "Value", keep = "Value3")
    expect_equal(top, expected_no_key %>%
        dplyr::select(chr, integration_locus, strand, Value, Value3))
    top <- top_integrations(test_df,
        n = 1, columns = "Value",
        key = "CompleteAmplificationID", keep = "Value2"
    )
    expect_equal(top, expected_key %>%
        dplyr::select(
            CompleteAmplificationID, chr,
            integration_locus, strand,
            Value, Value2
        ))
})
