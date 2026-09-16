# Offline tests; run from repository root.
env <- new.env()
sys.source("R_Scripts/05_Download_Metadata_Inventory.R", env)
stopifnot(env$geo_keyword_query(c("AML", "DMSO")) ==
            '("AML"[All Fields] AND "DMSO"[All Fields])')
calls <- list()
env$.geo_eutils <- function(endpoint, params, ...) {
  calls[[length(calls) + 1L]] <<- params
  if (endpoint == "esearch") {
    stopifnot(grepl("gse[ETYP]", params$term, fixed = TRUE))
    if (grepl("broken", params$term)) stop("Service unavailable")
    count <- if (grepl("empty", params$term)) 0 else 2
    ids <- if (count == 0) character() else c("100", "101")
    ids <- head(ids[seq_along(ids) > params$retstart], params$retmax)
    return(xml2::read_xml(paste0("<eSearchResult><Count>", count,
      "</Count><QueryTranslation>fixture</QueryTranslation><IdList>",
      paste(paste0("<Id>", ids, "</Id>")[seq_along(ids)], collapse = ""),
      "</IdList></eSearchResult>")))
  }
  xml2::read_xml(paste0('<eSummaryResult><DocSum><Id>100</Id>',
    '<Item Name="Accession">GSE1000</Item><Item Name="title">AML response</Item>',
    '<Item Name="summary">Cells exposed to dmso.</Item></DocSum>',
    '<DocSum><Id>101</Id><Item Name="Accession">GSE1001</Item>',
    '<Item Name="title">Other experiment</Item><Item Name="summary">Vehicle</Item>',
    '</DocSum></eSummaryResult>'))
}
env$.geo_listing <- function(url) {
  if (grepl("GSE1000/", url)) "GSE1000_series_matrix.txt.gz" else character()
}
result <- env$discover_geo_by_keywords(c("AML DMSO", "AML vehicle", "empty", "broken"),
                                       keywords = c("aml", "DMSO", "vehicle"), page_size = 1,
                                       retries = 0)
stopifnot(nrow(result$studies) == 2L, nrow(result$query_study_links) == 4L,
          identical(result$query_log$Status, c("complete", "complete", "no_hits", "search_failed")),
          result$studies$TitleKeywordHits[1] == "aml",
          result$studies$SummaryKeywordHits[1] == "DMSO",
          result$studies$QueryStrings[1] == "AML DMSO; AML vehicle",
          result$matrices$QueryStrings[1] == result$studies$QueryStrings[1],
          "no_matrix_files" %in% result$matrices$Status,
          any(vapply(calls, function(x) !is.null(x$retstart) && x$retstart == 1, logical(1))))
capped <- env$discover_geo_by_keywords("AML", max_studies_per_query = 1, retries = 0)
stopifnot(capped$query_log$Status == "truncated", nrow(capped$studies) == 1L)
empty <- env$discover_geo_by_keywords("empty", retries = 0)
stopifnot(nrow(empty$studies) == 0L, nrow(empty$matrices) == 0L)
# Summary failures remain visible even when no GSE could be resolved.
original <- env$.geo_eutils
env$.geo_eutils <- function(endpoint, params, ...) {
  if (endpoint == "esummary") stop("Summary failed")
  original(endpoint, params, ...)
}
failed <- env$discover_geo_by_keywords("AML", retries = 0)
stopifnot(all(failed$studies$MetadataStatus == "metadata_failed"),
          all(failed$studies$MetadataError == "Summary failed"), nrow(failed$matrices) == 0L)
cat("All keyword discovery tests passed.\n")
