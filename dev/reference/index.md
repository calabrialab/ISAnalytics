# Package index

## Import functions

Importing integration matrices and metadata from files

- [`import_association_file()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_association_file.md)
  **\[stable\]** : Import the association file from disk
- [`import_Vispa2_stats()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_Vispa2_stats.md)
  **\[stable\]** : Import Vispa2 stats given the aligned association
  file.
- [`import_single_Vispa2Matrix()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_single_Vispa2Matrix.md)
  **\[stable\]** : Import a single integration matrix from file
- [`import_parallel_Vispa2Matrices()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_parallel_Vispa2Matrices.md)
  **\[stable\]** : Import integration matrices from paths in the
  association file.

## Data cleaning and pre-processing

- [`annotation_issues()`](https://calabrialab.github.io/ISAnalytics/dev/reference/annotation_issues.md)
  **\[experimental\]** : Check for genomic annotation problems in IS
  matrices.
- [`compute_near_integrations()`](https://calabrialab.github.io/ISAnalytics/dev/reference/compute_near_integrations.md)
  **\[stable\]** : Scans input matrix to find and merge near integration
  sites.
- [`remove_collisions()`](https://calabrialab.github.io/ISAnalytics/dev/reference/remove_collisions.md)
  **\[stable\]** : Identifies and removes collisions
- [`realign_after_collisions()`](https://calabrialab.github.io/ISAnalytics/dev/reference/realign_after_collisions.md)
  **\[stable\]** : Re-aligns matrices of other quantification types
  based on the processed sequence count matrix.
- [`outlier_filter()`](https://calabrialab.github.io/ISAnalytics/dev/reference/outlier_filter.md)
  **\[experimental\]** : Filter out outliers in metadata, identified by
  the chosen outlier test.
- [`outliers_by_pool_fragments()`](https://calabrialab.github.io/ISAnalytics/dev/reference/outliers_by_pool_fragments.md)
  **\[stable\]** : Identify and flag outliers based on pool fragments.
- [`aggregate_metadata()`](https://calabrialab.github.io/ISAnalytics/dev/reference/aggregate_metadata.md)
  **\[stable\]** : Performs aggregation on metadata contained in the
  association file.
- [`aggregate_values_by_key()`](https://calabrialab.github.io/ISAnalytics/dev/reference/aggregate_values_by_key.md)
  **\[stable\]** : Aggregates matrices values based on specified key.

## Analysis functions

- [`compute_abundance()`](https://calabrialab.github.io/ISAnalytics/dev/reference/compute_abundance.md)
  **\[stable\]** : Computes the abundance for every integration event in
  the input data frame.
- [`top_integrations()`](https://calabrialab.github.io/ISAnalytics/dev/reference/top_integrations.md)
  **\[stable\]** : Sorts and keeps the top n integration sites based on
  the values in a given column.
- [`top_targeted_genes()`](https://calabrialab.github.io/ISAnalytics/dev/reference/top_targeted_genes.md)
  **\[experimental\]** : Top n targeted genes based on number of IS.
- [`sample_statistics()`](https://calabrialab.github.io/ISAnalytics/dev/reference/sample_statistics.md)
  **\[stable\]** : Computes user specified functions on numerical
  columns and updates the metadata data frame accordingly.
- [`CIS_grubbs()`](https://calabrialab.github.io/ISAnalytics/dev/reference/CIS_grubbs.md)
  **\[stable\]** : Grubbs test for Common Insertion Sites (CIS).
- [`CIS_grubbs_overtime()`](https://calabrialab.github.io/ISAnalytics/dev/reference/CIS_grubbs_overtime.md)
  **\[experimental\]** : Compute CIS and Grubbs test over different time
  points and groups.
- [`gene_frequency_fisher()`](https://calabrialab.github.io/ISAnalytics/dev/reference/gene_frequency_fisher.md)
  **\[experimental\]** : Compute Fisher's exact test on gene
  frequencies.
- [`cumulative_is()`](https://calabrialab.github.io/ISAnalytics/dev/reference/cumulative_is.md)
  **\[experimental\]** : Expands integration matrix with the cumulative
  IS union over time.
- [`is_sharing()`](https://calabrialab.github.io/ISAnalytics/dev/reference/is_sharing.md)
  **\[stable\]** : Sharing of integration sites between given groups.
- [`purity_filter()`](https://calabrialab.github.io/ISAnalytics/dev/reference/purity_filter.md)
  **\[stable\]** : Filter integration sites based on purity.
- [`iss_source()`](https://calabrialab.github.io/ISAnalytics/dev/reference/iss_source.md)
  **\[stable\]** : Find the source of IS by evaluating sharing.
- [`HSC_population_size_estimate()`](https://calabrialab.github.io/ISAnalytics/dev/reference/HSC_population_size_estimate.md)
  **\[stable\]** : Hematopoietic stem cells population size estimate.

## Plotting functions

- [`CIS_volcano_plot()`](https://calabrialab.github.io/ISAnalytics/dev/reference/CIS_volcano_plot.md)
  **\[stable\]** : Trace volcano plot for computed CIS data.
- [`top_cis_overtime_heatmap()`](https://calabrialab.github.io/ISAnalytics/dev/reference/top_cis_overtime_heatmap.md)
  **\[experimental\]** : Heatmaps for the top N common insertion sites
  over time.
- [`HSC_population_plot()`](https://calabrialab.github.io/ISAnalytics/dev/reference/HSC_population_plot.md)
  : Plot of the estimated HSC population size for each patient.
- [`integration_alluvial_plot()`](https://calabrialab.github.io/ISAnalytics/dev/reference/integration_alluvial_plot.md)
  **\[stable\]** : Alluvial plots for IS distribution in time.
- [`top_abund_tableGrob()`](https://calabrialab.github.io/ISAnalytics/dev/reference/top_abund_tableGrob.md)
  : Summary top abundant tableGrobs for plots.
- [`sharing_heatmap()`](https://calabrialab.github.io/ISAnalytics/dev/reference/sharing_heatmap.md)
  **\[stable\]** : Plot IS sharing heatmaps.
- [`sharing_venn()`](https://calabrialab.github.io/ISAnalytics/dev/reference/sharing_venn.md)
  **\[stable\]** : Produce tables to plot sharing venn or euler
  diagrams.
- [`circos_genomic_density()`](https://calabrialab.github.io/ISAnalytics/dev/reference/circos_genomic_density.md)
  **\[stable\]** : Trace a circos plot of genomic densities.
- [`fisher_scatterplot()`](https://calabrialab.github.io/ISAnalytics/dev/reference/fisher_scatterplot.md)
  **\[stable\]** : Plot results of gene frequency Fisher's exact test.

## Utility functions

- [`as_sparse_matrix()`](https://calabrialab.github.io/ISAnalytics/dev/reference/as_sparse_matrix.md)
  **\[stable\]** : Converts tidy integration matrices in the original
  sparse matrix form.
- [`generate_blank_association_file()`](https://calabrialab.github.io/ISAnalytics/dev/reference/generate_blank_association_file.md)
  : Create a blank association file.
- [`generate_Vispa2_launch_AF()`](https://calabrialab.github.io/ISAnalytics/dev/reference/generate_Vispa2_launch_AF.md)
  : Creates a reduced association file for a VISPA2 run, given project
  and pool
- [`inspect_tags()`](https://calabrialab.github.io/ISAnalytics/dev/reference/inspect_tags.md)
  : Retrieve description of a tag by name.
- [`available_tags()`](https://calabrialab.github.io/ISAnalytics/dev/reference/available_tags.md)
  : All available tags for dynamic vars look-up tables.
- [`set_mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/set_mandatory_IS_vars.md)
  [`set_annotation_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/set_mandatory_IS_vars.md)
  [`set_af_columns_def()`](https://calabrialab.github.io/ISAnalytics/dev/reference/set_mandatory_IS_vars.md)
  [`set_iss_stats_specs()`](https://calabrialab.github.io/ISAnalytics/dev/reference/set_mandatory_IS_vars.md)
  : Define custom dynamic vars.
- [`reset_mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/reset_mandatory_IS_vars.md)
  [`reset_annotation_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/reset_mandatory_IS_vars.md)
  [`reset_af_columns_def()`](https://calabrialab.github.io/ISAnalytics/dev/reference/reset_mandatory_IS_vars.md)
  [`reset_iss_stats_specs()`](https://calabrialab.github.io/ISAnalytics/dev/reference/reset_mandatory_IS_vars.md)
  [`reset_matrix_file_suffixes()`](https://calabrialab.github.io/ISAnalytics/dev/reference/reset_mandatory_IS_vars.md)
  [`reset_dyn_vars_config()`](https://calabrialab.github.io/ISAnalytics/dev/reference/reset_mandatory_IS_vars.md)
  : Resets dynamic vars to the default values.
- [`mandatory_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/mandatory_IS_vars.md)
  [`annotation_IS_vars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/mandatory_IS_vars.md)
  [`association_file_columns()`](https://calabrialab.github.io/ISAnalytics/dev/reference/mandatory_IS_vars.md)
  [`iss_stats_specs()`](https://calabrialab.github.io/ISAnalytics/dev/reference/mandatory_IS_vars.md)
  [`matrix_file_suffixes()`](https://calabrialab.github.io/ISAnalytics/dev/reference/mandatory_IS_vars.md)
  : Current dynamic vars specifications getters.
- [`set_matrix_file_suffixes()`](https://calabrialab.github.io/ISAnalytics/dev/reference/set_matrix_file_suffixes.md)
  : Sets the look-up table for matrix file suffixes.
- [`export_ISA_settings()`](https://calabrialab.github.io/ISAnalytics/dev/reference/export_ISA_settings.md)
  : Export a dynamic vars settings profile.
- [`import_ISA_settings()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_ISA_settings.md)
  : Import a dynamic vars settings profile.
- [`pcr_id_column()`](https://calabrialab.github.io/ISAnalytics/dev/reference/pcr_id_column.md)
  : Easily retrieve the name of the pcr id column.
- [`transform_columns()`](https://calabrialab.github.io/ISAnalytics/dev/reference/transform_columns.md)
  : Apply transformations to an arbitrary number of columns.
- [`comparison_matrix()`](https://calabrialab.github.io/ISAnalytics/dev/reference/comparison_matrix.md)
  **\[stable\]** : Obtain a single integration matrix from individual
  quantification matrices.
- [`separate_quant_matrices()`](https://calabrialab.github.io/ISAnalytics/dev/reference/separate_quant_matrices.md)
  **\[stable\]** : Separate a multiple-quantification matrix into single
  quantification matrices.
- [`generate_default_folder_structure()`](https://calabrialab.github.io/ISAnalytics/dev/reference/generate_default_folder_structure.md)
  : Generate a default folder structure, following VISPA2 standards
- [`enable_progress_bars()`](https://calabrialab.github.io/ISAnalytics/dev/reference/enable_progress_bars.md)
  : Enable global progress bars for ISAnalytics functions.
- [`NGSdataExplorer()`](https://calabrialab.github.io/ISAnalytics/dev/reference/NGSdataExplorer.md)
  : Launch the shiny application NGSdataExplorer.

## Exported variables and helpers

- [`date_formats()`](https://calabrialab.github.io/ISAnalytics/dev/reference/date_formats.md)
  :

  Possible choices for the `dates_format` parameter in
  `import_association_file`,
  `import_parallel_vispa2Matrices_interactive` and
  `import_parallel_vispa2Matrices_auto`.

- [`matching_options()`](https://calabrialab.github.io/ISAnalytics/dev/reference/matching_options.md)
  :

  Possible choices for the `matching_opt` parameter.

- [`quantification_types()`](https://calabrialab.github.io/ISAnalytics/dev/reference/quantification_types.md)
  :

  Possible choices for the `quantification_type` parameter.

- [`reduced_AF_columns()`](https://calabrialab.github.io/ISAnalytics/dev/reference/reduced_AF_columns.md)
  : Names of the columns of the association file to consider for Vispa2
  launch.

- [`default_stats()`](https://calabrialab.github.io/ISAnalytics/dev/reference/default_stats.md)
  :

  A set of pre-defined functions for `sample_statistics`.

- [`clinical_relevant_suspicious_genes()`](https://calabrialab.github.io/ISAnalytics/dev/reference/clinical_relevant_suspicious_genes.md)
  : Clinical relevant suspicious genes (for mouse and human).

- [`known_clinical_oncogenes()`](https://calabrialab.github.io/ISAnalytics/dev/reference/known_clinical_oncogenes.md)
  : Known clinical oncogenes (for mouse and human).

- [`available_outlier_tests()`](https://calabrialab.github.io/ISAnalytics/dev/reference/available_outlier_tests.md)
  : A character vector containing all the names of the currently
  supported outliers tests that can be called in the function
  outlier_filter.

- [`blood_lineages_default()`](https://calabrialab.github.io/ISAnalytics/dev/reference/blood_lineages_default.md)
  : Default blood lineages info

- [`default_iss_file_prefixes()`](https://calabrialab.github.io/ISAnalytics/dev/reference/default_iss_file_prefixes.md)
  : Default regex prefixes for Vispa2 stats files.

- [`default_meta_agg()`](https://calabrialab.github.io/ISAnalytics/dev/reference/default_meta_agg.md)
  : Default metadata aggregation function table

- [`default_report_path()`](https://calabrialab.github.io/ISAnalytics/dev/reference/default_report_path.md)
  : Default folder for saving ISAnalytics reports. Supplied as default
  argument for several functions.

- [`refGene_table_cols()`](https://calabrialab.github.io/ISAnalytics/dev/reference/refGene_table_cols.md)
  : Required columns for refGene file.

- [`default_af_transform()`](https://calabrialab.github.io/ISAnalytics/dev/reference/default_af_transform.md)
  : Default transformations to apply to association file columns.

- [`default_rec_agg_lambdas()`](https://calabrialab.github.io/ISAnalytics/dev/reference/default_rec_agg_lambdas.md)
  :

  Defaults for column aggregations in
  [`compute_near_integrations()`](https://calabrialab.github.io/ISAnalytics/dev/reference/compute_near_integrations.md).

## Deprecated

- [`import_parallel_Vispa2Matrices_interactive()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_parallel_Vispa2Matrices_interactive.md)
  **\[defunct\]** : Import integration matrices from association file.
- [`import_parallel_Vispa2Matrices_auto()`](https://calabrialab.github.io/ISAnalytics/dev/reference/import_parallel_Vispa2Matrices_auto.md)
  **\[defunct\]** : Import integration matrices from association file.
- [`cumulative_count_union()`](https://calabrialab.github.io/ISAnalytics/dev/reference/cumulative_count_union.md)
  **\[defunct\]** : Integrations cumulative count in time by sample
- [`unzip_file_system()`](https://calabrialab.github.io/ISAnalytics/dev/reference/unzip_file_system.md)
  **\[deprecated\]** : A utility function to unzip and use example file
  systems included in the package
- [`threshold_filter()`](https://calabrialab.github.io/ISAnalytics/dev/reference/threshold_filter.md)
  **\[deprecated\]** : Filter data frames with custom predicates

## Data

- [`association_file`](https://calabrialab.github.io/ISAnalytics/dev/reference/association_file.md)
  : Example of association file.
- [`integration_matrices`](https://calabrialab.github.io/ISAnalytics/dev/reference/integration_matrices.md)
  : Example of imported multi-quantification integration matrices.
- [`proto_oncogenes`](https://calabrialab.github.io/ISAnalytics/dev/reference/proto_oncogenes.md)
  [`tumor_suppressors`](https://calabrialab.github.io/ISAnalytics/dev/reference/proto_oncogenes.md)
  : Data frames for proto-oncogenes (human and mouse) and
  tumor-suppressor genes from UniProt.
- [`refGenes_hg19`](https://calabrialab.github.io/ISAnalytics/dev/reference/refGenes_hg19.md)
  [`refGenes_mm9`](https://calabrialab.github.io/ISAnalytics/dev/reference/refGenes_hg19.md)
  : Gene annotation files for hg19, mm9.
- [`refGenes_hg38`](https://calabrialab.github.io/ISAnalytics/dev/reference/refGenes_hg38.md)
  [`refGenes_mm10`](https://calabrialab.github.io/ISAnalytics/dev/reference/refGenes_hg38.md)
  : Reference gene annotation for hg38 or mm10.
