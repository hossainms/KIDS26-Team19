## Installs every R package used in this project.
## Run from the repo root:  Rscript R_Scripts/utils/install.pkgs.R

cran_pkgs <- c(
  # Shiny client
  "shiny",
  "DT",
  "DBI",
  "duckdb",
  # R_Scripts data pulls / reporting
  "dplyr",
  "data.table",
  "xml2",
  "readxl",
  "writexl"
)

bioc_pkgs <- c("GEOquery")

missing_cran <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(missing_cran)) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

missing_bioc <- setdiff(bioc_pkgs, rownames(installed.packages()))
if (length(missing_bioc)) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

message("All packages installed.")
