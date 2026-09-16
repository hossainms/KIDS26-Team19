# Opt-in network test; run from the repository root:
# Rscript R_Scripts/tests/test_geo_matrix_live.R
# Downloads and CSV reports are kept under downloads/ (ignored by Git).
source("R_Scripts/05_Download_Metadata_Inventory.R")
out <- "downloads/geo_matrix_smoke_test"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Bounded keyword search: verify metadata-to-file discovery without downloading
# arbitrary large matrices from keyword results.
query <- '(AML OR "acute myeloid leukemia") AND (DMSO OR vehicle) AND "Homo sapiens"[Organism]'
keyword_result <- discover_geo_by_keywords(query, keywords = c("AML", "acute myeloid leukemia", "DMSO", "vehicle"),
                                           max_studies_per_query = 2, retries = 1)
for (name in names(keyword_result)) {
  write.csv(keyword_result[[name]], file.path(out, paste0("keyword_", name, ".csv")), row.names = FALSE)
}
print(keyword_result$query_log)
stopifnot(keyword_result$query_log$Status %in% c("complete", "truncated"),
          nrow(keyword_result$studies) > 0,
          all(keyword_result$studies$MetadataStatus == "OK"),
          all(keyword_result$matrices$QueryStrings == query),
          any(keyword_result$matrices$Status == "discovered"))

inventory <- discover_geo_matrices(
  c("GSE1000", "GSE101791", "GSE999999999", "not-an-accession"),
  retries = 1, timeout_seconds = 60
)
write.csv(inventory, file.path(out, "discovery.csv"), row.names = FALSE)
print(inventory[, c("SeriesAccession", "Platform", "MatrixFile", "Status")])
stopifnot(sum(inventory$SeriesAccession == "GSE1000" & inventory$Status == "discovered") == 1L,
          sum(inventory$SeriesAccession == "GSE101791" & inventory$Status == "discovered") == 2L,
          all(c("GPL11154", "GPL20301") %in% inventory$Platform),
          inventory$Status[inventory$SeriesAccession == "GSE999999999"] == "discovery_failed",
          inventory$Status[inventory$SeriesAccession == "NOT-AN-ACCESSION"] == "invalid_accession")

result <- download_geo_matrices(inventory, file.path(out, "matrices"), retries = 1)
write.csv(result, file.path(out, "download_results.csv"), row.names = FALSE)
print(result[, c("MatrixFile", "Status", "ErrorMessage")])
good <- inventory$Status == "discovered"
stopifnot(all(result$Status[good] %in% c("downloaded", "already_exists")),
          all(file.exists(result$LocalFile[good])),
          identical(result$Status[!good], inventory$Status[!good]))

# Reruns should reuse valid files without modifying them.
before <- file.info(result$LocalFile[good])$mtime
rerun <- download_geo_matrices(result, file.path(out, "matrices"), retries = 0)
stopifnot(all(rerun$Status[good] == "already_exists"),
          identical(before, file.info(result$LocalFile[good])$mtime))
write.csv(rerun, file.path(out, "rerun_results.csv"), row.names = FALSE)

# Exercise a real missing-file response, followed by a valid cached row.
missing <- .geo_row("GSE1000", "discovered", file = "GSE1000-GPL999999999_series_matrix.txt.gz",
                    url = paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1000/matrix/",
                                 "GSE1000-GPL999999999_series_matrix.txt.gz"))
failure <- download_geo_matrices(rbind(missing, result[1, ]), file.path(out, "matrices"),
                                 retries = 0, timeout_seconds = 60)
write.csv(failure, file.path(out, "failure_results.csv"), row.names = FALSE)
stopifnot(identical(failure$Status, c("download_failed", "already_exists")),
          !file.exists(file.path(out, "matrices", missing$MatrixFile)))
cat("All live GEO checks passed. Reports and matrices:", normalizePath(out), "\n")
