# ISAnalytics maintenance notes

Working notes for potential fixes, documentation updates, and behaviour that
may deserve follow-up while reviewing the package.

## CIS_volcano_plot: `annotation_threshold_ontots` may not affect labels

- Area: `R/plotting-functions.R`, `CIS_volcano_plot()`.
- Status: needs verification / likely documentation-code mismatch.
- Observation: the argument `annotation_threshold_ontots` is documented as a
  value above which genes are annotated with labels. In the current
  implementation it is transformed into `annotation_threshold_ontots_log` and
  mentioned in the plot subtitle, but label selection appears to be driven by
  `tdist_fdr < significance_threshold`.
- User-visible impact: changing `annotation_threshold_ontots` may not change
  which genes are labelled, even though the documentation suggests it should.
- Possible follow-up:
  - decide whether labels should use `annotation_threshold_ontots`;
  - update the filtering logic if the argument is intended to control labels;
  - otherwise simplify or correct the documentation/subtitle.

## CIS_volcano_plot: unused `ggplot()` arguments warning

- Area: `R/plotting-functions.R`, `CIS_volcano_plot()`.
- Status: likely small bug / warning cleanup.
- Observation: `ggplot2::ggplot()` is called with `na.rm = TRUE` and `se = TRUE`.
  These arguments are not consumed by `ggplot()` and produce a warning from
  `fortify()`.
- Example warning:

```text
Warning message:
In fortify(data, ...) :
  Arguments in `...` must be used.
Problematic arguments:
* na.rm = TRUE
* se = TRUE
```

- User-visible impact: the plot is produced, but users see a confusing warning.
- Possible follow-up:
  - remove `na.rm = TRUE, se = TRUE` from the `ggplot()` call;
  - if NA handling is needed, pass `na.rm = TRUE` to the relevant geom layer;
  - `se` does not appear relevant for the current layers unless a smoothing
    layer is added.

## HSC_population_plot: unused `ggplot()` arguments warning

- Area: `R/plotting-functions.R`, `HSC_population_plot()`.
- Status: likely small bug / warning cleanup.
- Observation: the function also passes `na.rm = TRUE` and `se = TRUE` to
  `ggplot2::ggplot()`, which does not consume those arguments.
- User-visible impact: the same warning seen in `CIS_volcano_plot()` may appear
  when drawing HSC population estimates.
- Possible follow-up:
  - remove `na.rm = TRUE, se = TRUE` from the `ggplot()` call;
  - place `na.rm = TRUE` on `geom_point()` / `geom_line()` only if needed.

## HSC_population_plot: labels are partly hardcoded

- Area: `R/plotting-functions.R`, `HSC_population_plot()`.
- Status: likely documentation/UX improvement.
- Observation: axis labels and subtitle refer to a specific interpretation,
  including "Time Point (months after GT)", "Chao model with bias correction",
  and "IS from Myeloid PB cells as surrogate of HSC." The function arguments
  allow users to plot different `timepoints` and `models`, and the estimates
  data can contain different `CellType` / `Tissue` combinations.
- User-visible impact: plots can be misleading when the user selected a model
  other than the default or estimated a group other than Myeloid PB.
- Possible follow-up:
  - derive subtitle components from the filtered `estimates` data;
  - make axis/model labels configurable;
  - avoid saying "months" unless the selected `timepoint_column` is known to
    be expressed in months.

## HSC_population_size_estimate: empty output when `cell_type` matches marker, not `CellType`

- Area: `R/population-size-estimate.R`, `HSC_population_size_estimate()`.
- Status: documentation/logging improvement.
- Observation: the `cell_type` argument is matched against the `CellType`
  column after joining `blood_lineages_default()`, not against `CellMarker`.
  For example, the default lineage table maps `CellMarker == "MNC"` to
  `CellType == "Other"`. Calling the function with `cell_type = "MNC"` can
  therefore return `est = NULL` even though the input data has many MNC rows.
- User-visible impact: users may naturally pass their observed marker value
  and receive an empty result with little diagnostic information.
- Possible follow-up:
  - clarify in the parameter documentation that `cell_type` means lineage
    `CellType`, not raw `CellMarker`;
  - log the available `CellMarker`/`CellType`/`Tissue` combinations when no
    estimates are produced;
  - consider an explicit warning when requested `cell_type` values are not
    present after lineage mapping.

## circos_genomic_density: chromosome prefix is always added

- Area: `R/plotting-functions.R`, `circos_genomic_density()`.
- Status: potential robustness improvement.
- Observation: the internal data preparation step builds circos coordinates
  with `chr = paste0("chr", .data$chr)`. This works for package-standard
  matrices where chromosomes are stored as `"1"`, `"2"`, `"X"`, etc., but
  custom inputs that already contain UCSC-style names such as `"chr1"` would
  become `"chrchr1"`.
- User-visible impact: circos initialization or density plotting can fail or
  silently miss data for custom matrices with pre-prefixed chromosome names.
- Possible follow-up:
  - normalize chromosome names with a helper that adds `chr` only when absent;
  - document the expected chromosome naming convention for this plot.

## generate_Vispa2_launch_AF: missing-column error may report the wrong set

- Area: `R/utility-functions.R`, `generate_Vispa2_launch_AF()`.
- Status: needs verification.
- Observation: after building `to_check`, the missing-column branch calls
  `.missing_af_needed_cols(to_check[!to_check %in% association_file_columns(TRUE)])`.
  `association_file_columns(TRUE)` returns the full dynamic-vars table, not just
  a character vector of column names. The intended comparison may have been
  against `colnames(association_file)` or `association_file_columns()`.
- User-visible impact: when required columns are missing, the error message may
  list an unexpected or overly broad set of columns.
- Possible follow-up:
  - add a test with a deliberately incomplete association file;
  - ensure the missing-column error reports exactly `to_check[!to_check %in%
    colnames(association_file)]`.

## top_cis_overtime_heatmap: re-drawing can overlay existing grid output

- Area: `R/plotting-functions.R`, `top_cis_overtime_heatmap()`.
- Status: mostly expected `grid`/`pheatmap` behaviour, but documentation could
  be clearer.
- Observation: returned heatmaps are `pheatmap` objects containing a `gtable`.
  Re-drawing one with `grid::grid.draw(hmaps$PT001$gtable)` draws on the current
  grid page. If a previous plot is still present and `grid::grid.newpage()` is
  not called first, output can be overlaid.
- User-visible impact: heatmaps may appear plotted on top of a previous plot
  when manually re-drawn.
- Possible follow-up:
  - document the recommended re-draw pattern more explicitly:

```r
grid::grid.newpage()
grid::grid.draw(hmaps$PT001$gtable)
```

  - consider examples with `silent = TRUE` when users want to store heatmaps
    without drawing them immediately.
