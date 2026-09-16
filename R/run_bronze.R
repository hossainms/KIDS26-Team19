# Download every study in the GSE inventory into bronze/.
# Run from the repo root: Rscript R/run_bronze.R
source("R/bronze.R")
gse.ids <- unique(readxl::read_excel("GSE_Metadata_Inventory.xlsx", sheet = "Metadata")$SeriesAccession)
build_bronze(gse.ids)
