RecalibrationMap <- R6::R6Class(
    classname = "RecalibrationMap",
    private = list(
        .map = NULL,
        .req_tags = NULL,
        .generate_rec_map_filename = function() {
            def <- "recalibration_map.tsv.gz"
            date <- lubridate::today()
            return(paste0(date, "_", def))
        }
    ),
    public = list(
        initialize = function(req_tags) {
            private$.req_tags <- req_tags
            with_types <- req_tags |>
                dplyr::left_join(.types_mapping(), by = "types") |>
                purrr::pmap(function(...) {
                    row <- list(...)
                    type_mapping <- row$fread
                    names(type_mapping) <- row$names
                    return(type_mapping)
                }) |>
                purrr::list_c()
            col_declaration <- purrr::map2(with_types, names(with_types), ~ {
                content <- rlang::exec(.fn = .x)
                return(content)
            })
            col_declaration_before <- col_declaration
            names(col_declaration_before) <- paste0(
                names(col_declaration_before),
                "_before"
            )
            col_declaration_after <- col_declaration
            names(col_declaration_after) <- paste0(
                names(col_declaration_after),
                "_after"
            )
            rec_map <- tibble::tibble(
                !!!col_declaration_before,
                !!!col_declaration_after
            )
            private$.map <- rec_map
        },
        get_map = function() {
            return(private$.map)
        },
        update = function(before, after) {
            if (nrow(before) > nrow(after)) {
                # In this case after is composed of only 1 row, recycle it to match
                # the length
                after <- after |>
                    dplyr::slice(
                        rep(1, each = nrow(before))
                    )
            }
            # Rename cols
            before <- before |>
                dplyr::rename_with(.fn = ~ paste0(.x, "_before"))
            after <- after |>
                dplyr::rename_with(.fn = ~ paste0(.x, "_after"))
            # Obtain a single df
            update_df <- before |>
                dplyr::bind_cols(after)
            # Update map
            private$.map <- private$.map |>
                dplyr::bind_rows(update_df) |>
                dplyr::distinct()
        },
        write_recalibr_map = function(file_path) {
            if (!fs::file_exists(file_path)) {
                ext <- fs::path_ext(file_path)
                ## if ext is empty assume it's a folder
                if (ext == "") {
                    fs::dir_create(file_path)
                    gen_filename <- private$.generate_rec_map_filename()
                    tmp_filename <- fs::path(file_path, gen_filename)
                } else {
                    tmp_filename <- fs::path_ext_remove(file_path)
                    if (ext %in% .compressed_formats()) {
                        ext <- paste(fs::path_ext(tmp_filename), ext, sep = ".")
                        tmp_filename <- fs::path_ext_remove(tmp_filename)
                    }
                    if (!ext %in% c(
                        "tsv", paste("tsv", .compressed_formats(), sep = "."),
                        "csv", paste("csv", .compressed_formats(), sep = "."),
                        "txt", paste("txt", .compressed_formats(), sep = ".")
                    )) {
                        if (getOption("ISAnalytics.verbose", TRUE) == TRUE) {
                            warn <- c("Recalibration file format unsupported",
                                i = "Writing file in 'tsv.gz' format"
                            )
                            rlang::inform(warn, class = "rec_unsupp_ext")
                        }
                        tmp_filename <- paste0(tmp_filename, ".tsv.gz")
                    } else {
                        tmp_filename <- paste0(tmp_filename, ".", ext)
                    }
                }
            } else if (fs::is_dir(file_path)) {
                gen_filename <- private$.generate_rec_map_filename()
                tmp_filename <- fs::path(file_path, gen_filename)
            } else {
                tmp_filename <- file_path
            }
            withRestarts(
                {
                    readr::write_delim(
                        private$.map,
                        file = tmp_filename, delim = "\t", na = ""
                    )
                    saved_msg <- paste(
                        "Recalibration map saved to: ",
                        tmp_filename
                    )
                    if (getOption("ISAnalytics.verbose", TRUE) == TRUE) {
                        rlang::inform(saved_msg)
                    }
                },
                skip_write = function() {
                    skip_msg <- paste(
                        "Could not write recalibration map file.",
                        "Skipping."
                    )
                    rlang::inform(skip_msg)
                }
            )
        }
    )
)
