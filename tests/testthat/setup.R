library(ISAnalytics)

# Load AF
af_sym <- "association_file"
utils::data(list = af_sym, envir = rlang::current_env())
association_file <- rlang::eval_tidy(rlang::sym(af_sym))

# Load matrices
matrices_sym <- "integration_matrices"
utils::data(list = matrices_sym, envir = rlang::current_env())
integration_matrices <- rlang::eval_tidy(rlang::sym(matrices_sym))

# Generate file folder structure
fs_path <- generate_default_folder_structure(type = "both")

# Align AF
local_af_corr <- import_association_file(fs_path$af,
    root = fs_path$root_corr,
    report_path = NULL
)
local_af_inc <- import_association_file(fs_path$af,
    root = fs_path$root_inc,
    report_path = NULL
)

# Unzip testdata
testdata_path <- tempdir()
utils::unzip(zipfile = system.file("testdata", "testdata.zip",
                                   package = "ISAnalytics"),
             exdir = testdata_path)
