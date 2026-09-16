# Bronze layer: download GEO series matrix files exactly as published.
# Files land in bronze/geo/<GSE>/ and each one gets a row in bronze/manifest.csv.

GEO_SERIES_URL <- "https://ftp.ncbi.nlm.nih.gov/geo/series"

# GEO groups studies by all but the last three digits: GSE100446 -> GSE100nnn.
gse_folder <- function(gse.id) {
  digits <- sub("^GSE", "", gse.id)
  prefix <- if (nchar(digits) > 3) substr(digits, 1, nchar(digits) - 3) else ""
  paste0("GSE", prefix, "nnn")
}

matrix_dir_url <- function(gse.id) {
  paste(GEO_SERIES_URL, gse_folder(gse.id), gse.id, "matrix", "", sep = "/")
}

# Multi-platform studies publish one file per platform, so list the folder.
parse_matrix_listing <- function(html) {
  hrefs <- unlist(regmatches(html, gregexpr('href="[^"]+"', html)))
  files <- gsub('^href="|"$', "", hrefs)
  unique(files[grepl("_series_matrix\\.txt\\.gz$", files)])
}

platform_from_file <- function(file.name) {
  m <- regmatches(file.name, regexpr("GPL[0-9]+", file.name))
  if (length(m) == 0) NA_character_ else m
}

list_matrix_files <- function(gse.id) {
  html <- tryCatch(readLines(matrix_dir_url(gse.id), warn = FALSE), error = function(e) NULL)
  if (is.null(html)) return(character())
  parse_matrix_listing(html)
}

# Download every matrix file for one study. Existing files are kept.
download_study <- function(gse.id, bronze.dir = "bronze") {
  files <- list_matrix_files(gse.id)
  if (length(files) == 0) {
    return(data.frame(SeriesAccession = gse.id, Platform = NA, MatrixFile = NA, URL = matrix_dir_url(gse.id),
                      Bytes = NA, MD5 = NA, DownloadedAt = NA, Status = "no matrix files listed"))
  }
  out.dir <- file.path(bronze.dir, "geo", gse.id)
  dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
  rows <- lapply(files, function(f) {
    url <- paste0(matrix_dir_url(gse.id), f)
    dest <- file.path(out.dir, f)
    status <- "OK"
    if (!file.exists(dest)) {
      tmp <- paste0(dest, ".part")
      ok <- tryCatch(download.file(url, tmp, mode = "wb", quiet = TRUE) == 0, error = function(e) FALSE)
      if (ok) file.rename(tmp, dest) else { unlink(tmp); status <- "download failed" }
    }
    have <- file.exists(dest)
    data.frame(SeriesAccession = gse.id, Platform = platform_from_file(f), MatrixFile = f, URL = url,
               Bytes = if (have) file.size(dest) else NA,
               MD5 = if (have) unname(tools::md5sum(dest)) else NA,
               DownloadedAt = if (have) format(file.mtime(dest), "%Y-%m-%dT%H:%M:%S") else NA,
               Status = status)
  })
  do.call(rbind, rows)
}

# Download a list of studies and write the manifest.
build_bronze <- function(gse.ids, bronze.dir = "bronze", delay.seconds = 0.4) {
  rows <- vector("list", length(gse.ids))
  for (i in seq_along(gse.ids)) {
    message(i, "/", length(gse.ids), " ", gse.ids[i])
    rows[[i]] <- download_study(gse.ids[i], bronze.dir)
    Sys.sleep(delay.seconds)
  }
  manifest <- do.call(rbind, rows)
  write.csv(manifest, file.path(bronze.dir, "manifest.csv"), row.names = FALSE, na = "")
  manifest
}
