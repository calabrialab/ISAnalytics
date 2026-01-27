# A utility function to unzip and use example file systems included in the package

**\[deprecated\]** From `ISAnalytics 1.5.4` this function is defunct,
since the package doesn't include example tabular files anymore. Use the
function
[`generate_default_folder_structure()`](https://calabrialab.github.io/ISAnalytics/dev/reference/generate_default_folder_structure.md)
to generate a default folder structure for running tests and play with
the package import functions. If you don't need to test import
functions, you can simply load package included data via
`data("integration_matrices")` or `data("association_file")`.

## Usage

``` r
unzip_file_system(zipfile, name)
```

## Arguments

- zipfile:

  The zipped file to decompress

- name:

  The name of the folder in the zipped archive ("fs" or "fserr")

## Value

A path to reference
