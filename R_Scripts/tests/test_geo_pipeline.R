# Offline end-to-end reporting and interrupted-rerun checks.
env <- new.env()
sys.source("R_Scripts/05_Download_Metadata_Inventory.R", env)
scratch <- tempfile("geo-pipeline-")
dir.create(scratch)
tryCatch({
  studies <- c("GSE1", "GSE2")
  make_search <- function(ids) {
    matrices <- do.call(rbind, lapply(ids, function(gse) {
      env$.geo_row(gse, "discovered", file = paste0(gse, "_series_matrix.txt.gz"),
                   url = paste0("https://fixture/", gse))
    }))
    if (is.null(matrices)) matrices <- env$.geo_empty()
    list(studies = data.frame(SeriesAccession = ids, Title = rep("Old title", length(ids)),
                              MetadataStatus = rep("OK", length(ids))), matrices = matrices,
         query_log = data.frame(Status = if (length(ids)) "complete" else "no_hits"))
  }
  env$discover_geo_by_keywords <- function(...) make_search(studies)
  env$.geo_fetch <- function(url, destination) {
    con <- gzfile(destination, "wt")
    on.exit(close(con))
    writeLines(c('!Series_title\t"Actual title"', '!Sample_platform_id\t"GPL96"\t"GPL96"',
                 '!Sample_organism_ch1\t"Homo sapiens"\t"Homo sapiens"',
                 '!Sample_geo_accession\t"GSM1"\t"GSM2"', '!series_matrix_table_begin',
                 '"ID_REF"\t"GSM1"\t"GSM2"', '"probe1"\t1\t2', '!series_matrix_table_end'), con)
  }
  original_download <- env$download_geo_matrices
  env$download_geo_matrices <- function(inventory, dest_dir, ...) {
    original_download(inventory, dest_dir, retries = 0, delay_seconds = 0)
  }
  env$run_geo_pipeline(out_dir = scratch, query_strings = "fixture")
  read <- function(name) read.csv(file.path(scratch, name), stringsAsFactors = FALSE)
  success <- read("downloaded_matrices.csv")
  stopifnot(nrow(success) == 2L, all(success$Platform == "GPL96"), all(success$SampleCount == 2L),
            all(success$SeriesTitle == "Actual title"), all(success$Organism == "Homo sapiens"),
            nrow(read("failures.csv")) == 0L, nrow(read("platform_failures.csv")) == 0L)

  # A rerun adds GSE3, then is interrupted on GSE4. Prior rows must survive.
  studies <- c("GSE3", "GSE4")
  env$download_geo_matrices <- function(inventory, dest_dir, ...) {
    if (inventory$SeriesAccession == "GSE4") stop("Simulated interruption")
    original_download(inventory, dest_dir, retries = 0, delay_seconds = 0)
  }
  interrupted <- tryCatch({env$run_geo_pipeline(out_dir = scratch, query_strings = "fixture"); FALSE},
                           error = function(e) grepl("Simulated interruption", conditionMessage(e)))
  success <- read("downloaded_matrices.csv")
  stopifnot(interrupted, setequal(success$SeriesAccession, c("GSE1", "GSE2", "GSE3")))

  # Reprocessing cached files does not duplicate report rows.
  studies <- "GSE1"
  env$run_geo_pipeline(out_dir = scratch, query_strings = "fixture")
  success <- read("downloaded_matrices.csv")
  stopifnot(nrow(success) == 3L, !anyDuplicated(success$LocalFile),
            success$Status[success$SeriesAccession == "GSE1"] == "already_exists")
  studies <- character()
  env$run_geo_pipeline(out_dir = scratch, query_strings = "fixture")
  stopifnot(nrow(read("downloaded_matrices.csv")) == 3L)
  empty_dir <- file.path(scratch, "empty")
  env$run_geo_pipeline(out_dir = empty_dir, query_strings = "fixture")
  stopifnot(nrow(read.csv(file.path(empty_dir, "downloaded_matrices.csv"))) == 0L)
  cat("All standalone pipeline tests passed.\n")
}, finally = unlink(scratch, recursive = TRUE))
