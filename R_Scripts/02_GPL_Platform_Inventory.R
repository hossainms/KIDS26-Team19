#Setup
rm(list = ls())
options(stringsAsFactors = FALSE)

library(readxl)
library(dplyr)
library(data.table)
library(writexl)

#Settings

candidate.file <- "/home/nmakonne/AML_GEO_Project/GEO_DataSets/Inputs/AML_GEO_candidates_datasets_2026-07-24.xlsx"
geo.dir <- "/home/nmakonne/AML_GEO_Project/GEO_DataSets"

#Read workbook
geo.dt <- read_excel(candidate.file, sheet = "human_expression_methylation")

#Keep only expression studies
expression.dt <- geo.dt %>%
  filter(
    grepl(
      "Expression profiling",
      Type,
      ignore.case = TRUE
    )
  )
#Unique GSE list:
gse.ids <- unique(expression.dt$SeriesAccession)

####Helper Functions
##Build FTP URL
build_gse_url <- function(gse.id) {
  paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/",
    substring(
      gse.id,
      1,
      nchar(gse.id) - 3
    ),
    "nnn/",
    gse.id,
    "/matrix/",
    gse.id,
    "_series_matrix.txt.gz"
  )
}
##Download if needed:

download_gse_if_needed <- function(gse.id, geo.dir) {
  file.name <- file.path(geo.dir, paste0(gse.id,"_series_matrix.txt.gz"))
  if (!file.exists(file.name)) {
    message("Downloading", gse.id, "...")
    download.file(url = build_gse_url(gse.id),
                  destfile = file.name,
                  mode = "wb",
                  quiet = FALSE)
  } else {
    message(gse.id, " already exists.")
  }
  return(file.name)
}

##Read only the header:

read_series_header <- function(file.name) {
  con <- gzfile(file.name, open = "rt")
  header <- readLines(con, n = 250, warn = FALSE)  #Reading first 250 lines to include series info
  close(con)
  return(header)
}

##Extracting the metadata:
extract_header_metadata <- function(
    header, gse.id) {
  ###platform ids
  platform.lines <-
    grep(
      "^!Series_platform_id",
      header,
      value = T
    )
  platforms <- sub(
    "!Series_platform_id = ",
    "",
    platform.lines,
    fixed = T
  )
  ###study title
  title.line <-
    grep(
      "^!Series_title",
      header,
      value = T
    )
  title <- if (length(title.line) == 0) {
    NA_character_
  } else {
    sub(
      "!Series_title = ",
      "",
      title.line[1],
      fixed = T
    )
  }
  ##Sample count
  sample.count <-
    sum(
      grepl(
        "^!Sample_geo_accession",
        header
      )
    )
  #Return one row per platform
  data.table(
    SeriesAccession = gse.id,
    Platform = platforms,
    SeriesTitle = title,
    SampleCount = sample.count
  )
}

### Main wrapper function

get_gse_header_metadata <- function(gse.id, geo.dir) {
  file.name <- download_gse_if_needed(gse.id, geo.dir)
  header <- read_series_header(file.name)
  result <- extract_header_metadata(header, gse.id)
  return(result)
}

###Test

test.result <- get_gse_header_metadata(gse.id = gse.id[1], geo.dir = geo.dir)
print(test.result)






















