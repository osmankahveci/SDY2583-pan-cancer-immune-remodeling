# Install packages used by the recovered and reconstructed SDY2583 workflows.
# Exact versions are captured after execution with scripts/99_session_info.R.
cran_packages <- c(
  "broom", "dplyr", "forcats", "ggplot2", "patchwork", "purrr",
  "readr", "rlang", "scales", "stringr", "tibble", "tidyr"
)
missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_cran) > 0L) install.packages(missing_cran)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc_packages <- c("flowCore")
for (package in bioc_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    BiocManager::install(package, ask = FALSE, update = FALSE)
  }
}
message("All-panel dependency check complete.")
