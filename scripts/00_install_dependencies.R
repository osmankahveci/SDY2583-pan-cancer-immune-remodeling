# Install packages used by the recovered SDY2583 analysis scripts.
# Package versions are intentionally not asserted because the original session
# metadata were not preserved in the analysis archive.

cran_packages <- c(
  "broom",
  "dplyr",
  "forcats",
  "ggplot2",
  "patchwork",
  "purrr",
  "readr",
  "rlang",
  "scales",
  "stringr",
  "tibble",
  "tidyr"
)

missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("flowCore", quietly = TRUE)) {
  BiocManager::install("flowCore", ask = FALSE, update = FALSE)
}

message("Dependency check complete.")
