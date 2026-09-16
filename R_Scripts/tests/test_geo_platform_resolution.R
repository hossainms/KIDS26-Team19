# Offline tests: no GEO requests or matrix downloads.
env <- new.env()
sys.source("R_Scripts/05_Download_Metadata_Inventory.R", env)
scratch <- tempfile("geo-platform-tests-")
dir.create(scratch)
tryCatch({
  fixture <- function(name, header) {
    path <- file.path(scratch, name)
    con <- gzfile(path, "wt")
    writeLines(c(header, "!series_matrix_table_begin", "!series_matrix_table_end"), con)
    close(con)
    path
  }
  files <- c(
    fixture("series.gz", c(rep('!Series_title\t"Padding"', 300), '!Series_platform_id\t"GPL96"')),
    fixture("sample.gz", c('!Series_platform_id = GPL96', '!Series_platform_id = GPL97',
                           '!Sample_platform_id\t"GPL97"\t"GPL97"')),
    fixture("fallback.gz", '!Sample_geo_accession\t"GSM1"'),
    fixture("fallback2.gz", '!Sample_geo_accession\t"GSM2"'),
    fixture("ambiguous.gz", '!Series_title\t"Missing sample identities"'),
    fixture("failure.gz", '!Series_title\t"Failure fixture"'),
    fixture("single.gz", '!Series_title\t"Single platform"')
  )
  queries <- character()
  env$.geo_query_platform_metadata <- function(gse) {
    queries <<- c(queries, gse)
    if (gse == "GSE999") stop("Network failure")
    if (gse == "GSE3") return(list(series = "GPL570", samples = list()))
    list(series = c("GPL96", "GPL97"), samples = list(GSM1 = "GPL96", GSM2 = "GPL97"))
  }
  input <- data.frame(SeriesAccession = c("GSE1", "GSE1", "GSE2", "GSE2", "GSE2", "GSE999", "GSE3"),
                      Platform = NA_character_, LocalFile = files, Status = "downloaded")
  output <- env$fill_geo_platforms(input, retries = 0, delay_seconds = 0)
  stopifnot(identical(output$Platform, c("GPL96", "GPL97", "GPL96", "GPL97", NA, NA, "GPL570")),
            identical(output$PlatformSource, c("matrix_series_header", "matrix_sample_header",
              "GEOquery_samples", "GEOquery_samples", "unresolved", "unresolved", "GEOquery_series")),
            sum(queries == "GSE2") == 1L, !"GSE1" %in% queries,
            all(output$Status == "downloaded"), grepl("Network failure", output$PlatformError[6]),
            grepl("unambiguous", output$PlatformError[5]))
  # Existing IDs are kept; unsuccessful download rows trigger no metadata requests.
  input$Platform <- "GPL123"
  count <- length(queries)
  kept <- env$fill_geo_platforms(input)
  input$Platform <- NA_character_
  input$Status <- "download_failed"
  skipped <- env$fill_geo_platforms(input)
  stopifnot(all(kept$Platform == "GPL123"), length(queries) == count,
            all(is.na(skipped$Platform)), nrow(env$fill_geo_platforms(input[FALSE, ])) == 0L)
  cat("All platform resolution tests passed.\n")
}, finally = unlink(scratch, recursive = TRUE))
