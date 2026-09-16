# Run from the repository root:
# Rscript R_Scripts/tests/test_geo_matrix_download.R
# Offline fixtures exercise failures without contacting NCBI.
env <- new.env()
sys.source("R_Scripts/05_Download_Metadata_Inventory.R", envir = env)
scratch <- tempfile("geo-matrix-tests-")
dir.create(scratch)
on.exit_cleanup <- function() unlink(scratch, recursive = TRUE)

tryCatch({
  urls <- character()
  env$.geo_listing <- function(url) {
    urls <<- c(urls, url)
    if (grepl("/GSE982/", url)) return(c("../", "GSE982_series_matrix.txt.gz"))
    if (grepl("/GSE1000/", url)) return(c(
      "GSE1000-GPL96_series_matrix.txt.gz", "GSE1000-GPL97_series_matrix.txt.gz",
      "GSE1000-GPL97_series_matrix.txt.gz", "unrelated.txt", "../bad_series_matrix.txt.gz"
    ))
    if (grepl("/GSE1001/", url)) return("../")
    stop("HTTP 404 fixture")
  }
  inventory <- env$discover_geo_matrices(
    c(" gse982 ", "GSE982", "GSE1000", "GSE1001", "GSE1002", "bad", NA),
    retries = 0, delay_seconds = 0
  )
  stopifnot(nrow(inventory) == 7L,
            sum(inventory$Status == "discovered") == 3L,
            sum(inventory$Status == "invalid_accession") == 2L,
            "no_matrix_files" %in% inventory$Status,
            "discovery_failed" %in% inventory$Status,
            identical(inventory$Platform[2:3], c("GPL96", "GPL97")),
            is.na(inventory$Platform[1]),
            any(grepl("/GSEnnn/GSE982/", urls)),
            any(grepl("/GSE1nnn/GSE1000/", urls)),
            nrow(env$discover_geo_matrices(character())) == 0L)

  attempts <- 0L
  env$.geo_fetch <- function(url, destination) {
    if (grepl("GPL97", url)) {
      writeLines("partial response", destination)
      stop("Network interruption fixture")
    }
    if (grepl("GPL96", url)) {
      attempts <<- attempts + 1L
      if (attempts == 1L) stop("Temporary failure fixture")
    }
    con <- gzfile(destination, "wt")
    on.exit(close(con))
    writeLines(c('!Series_title\t"Fixture"', '!Series_platform_id\t"GPL96"', '!series_matrix_table_begin',
                 '"ID_REF"\t"GSM1"', '"probe1"\t1', '!series_matrix_table_end'), con)
  }
  result <- env$download_geo_matrices(inventory, scratch, retries = 1, delay_seconds = 0)
  stopifnot(identical(result$Status[1:3], c("downloaded", "downloaded", "download_failed")),
            result$Platform[1] == "GPL96", result$PlatformSource[1] == "matrix_series_header",
            attempts == 2L, all(file.exists(result$LocalFile[1:2])),
            is.na(result$LocalFile[3]), identical(result$Status[4:7], inventory$Status[4:7]),
            length(list.files(scratch, all.files = TRUE, no.. = TRUE)) == 2L)
  again <- env$download_geo_matrices(result, scratch, retries = 0, delay_seconds = 0)
  stopifnot(all(again$Status[1:2] == "already_exists"))

  # Corrupt cached files are replaced, not silently skipped.
  writeLines("broken cache", result$LocalFile[1])
  repaired <- env$download_geo_matrices(result[1, ], scratch, retries = 0, delay_seconds = 0)
  stopifnot(repaired$Status == "downloaded")
  env$.geo_validate_matrix(repaired$LocalFile)

  # A failed refresh must preserve the existing valid file.
  before <- tools::md5sum(repaired$LocalFile)
  env$.geo_fetch <- function(url, destination) writeLines("HTML error page", destination)
  failed <- env$download_geo_matrices(repaired, scratch, overwrite = TRUE,
                                     retries = 0, delay_seconds = 0)
  stopifnot(failed$Status == "download_failed", identical(before, tools::md5sum(repaired$LocalFile)),
            length(list.files(scratch, all.files = TRUE, no.. = TRUE)) == 2L)
  cat("All GEO matrix discovery/download tests passed.\n")
}, finally = on.exit_cleanup())
