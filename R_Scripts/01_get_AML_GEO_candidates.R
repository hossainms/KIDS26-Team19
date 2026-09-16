### Setup
setwd("C:/Users/nmakonne/Box/Nobel/Rscratchwork/R Scripts")
rm(list = ls())
options(stringsAsFactors = FALSE)
library(writexl)
library(dplyr)
library(xml2) #(Dora suggested I read xml instead of html)
library(data.table) #Same thing
library(GEOquery) #Getting the platforms from GEOquery instead of Eutils

#User Settings
out.dir <- "C:/Users/nmakonne/Box/Nobel/Rscratchwork/outputs"

if(!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = TRUE)
}

#Diagnosis abbreviation to search for. Only "AML" and "ALL" are supported for now.
diagnosis <- "AML"

#Map of supported diagnosis abbreviations to their full names:
diagnosis.map <- c(
  AML = "acute myeloid leukemia",
  ALL = "acute lymphocytic leukemia"
)

out.path <- file.path(
  out.dir,
  paste0(diagnosis, "_GEO_candidates_datasets_", Sys.Date(), ".xlsx")
)

#Build the search terms for a given diagnosis abbreviation:
build_query_strings <- function(abbreviation, diagnosis.map) {
  if (!abbreviation %in% names(diagnosis.map)) {
    stop(
      "Unsupported diagnosis '", abbreviation, "'. Supported values: ",
      paste(names(diagnosis.map), collapse = ", ")
    )
  }
  full.name <- diagnosis.map[[abbreviation]]
  c(
    paste(abbreviation, "DMSO"),
    paste0("\"", full.name, "\" DMSO"),
    paste(abbreviation, "vehicle"),
    paste0("\"", full.name, "\" vehicle")
  )
}

#The Search Terms

query.strings <- build_query_strings(diagnosis, diagnosis.map)

#Optional NCBI API Key 
ncbi.api.key <- Sys.getenv("NCBI_API_KEY")

###Helper functions

#Add API key to a URL if available
add_api_key <- function(url, api.key = ""){
  if (nzchar(api.key)) {
    url <- paste0(url, "&api_key=", api.key)
  }
  return(url)
}

#Safely read xml from a url
safe_read_xml <- function(url) {
  result <- tryCatch(
    read_xml(url),
    error = function(e) {
      warning("Could not read XML from: ", url)
      return(NULL)
    }
  )
  return(result)
}

###Function #1, Searching GEO/GDS and returning NCBI UIDs

get_gds_ids <- function(query.string, retmax = 9999, api.key = "") {
  message("Searching GEO/GDS for:", query.string) # Searching message to be fancy
  #Building the NCBI Esearch URL:
  url <- paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?",
    "db=gds",
    "&term=", URLencode(query.string, reserved = TRUE),
    "&retmode=xml",
    "&retmax=", retmax
  )
  #Adding the NCBI API key if one is available:
  url <- add_api_key(url, api.key)
  #Reading the XML returned by NCBI:
  xml_doc <- safe_read_xml(url)
  #Returning an empty data table if the XML couldn't be read:
  if (is.null(xml_doc)) {
    return(
      data.table(
        query_string = query.string,
        ID = character()
      )
    )
  }
  
  #looking for any nodes named Id in the doc:
  id_nodes <- xml_find_all(xml_doc, ".//Id")
  #Extracting text:
  ids <- xml_text(id_nodes, trim = TRUE)
  # removing duplicate ids:
  ids <- unique(ids)
  
  #Returning a table with the query and all the IDs it found:
  return(
    data.table(
      query_string = query.string,
      ID = ids
    )
  )
}

###Another Helper function to extract one named item from an Esummary record:
get_item <- function(record, item.name) {
  #Building an XPath query to find an Item with a specific name:
  xpath <- sprintf('.//Item[@Name="%s"]', item.name)
  node <- xml_find_first(record, xpath) #returns first matching node
  if (inherits(node, "xml_missing")) {   #If the item is missing, returning a xml_missing object
    return(NA_character_)
  }
  return(xml_text(node, trim = TRUE)) #If not missing then return the text inside the node:
}

### Another Helper function: trying multiple possible Item names:

get_first_item <- function(record, item.names) {
  #Looping through the possible item names
  for (item.name in item.names) {
    value <- get_item(record, item.name)
    if (!is.na(value) && nzchar(value)) {
      return(value)
    }
  }
  #If nothing works then return missing text:
  return(NA_character_)
}

### Another Helper funciton: extracting PubMed IDs

get_pubmed_ids <- function(record) {
  #Finding item nodes who has a Pubmed or PMID in the name
  pubmed_nodes <- xml_find_all(
    record,
    './/Item[
    contains(translate(@Name, "ABCDEFGHIJKLMNOPQRSTUVWXYZ", "abcdefghijklmnopqrstuvwxyz"), "pubmed")
    or
    contains(translate(@Name, "ABCDEFGHIJKLMNOPQRSTUVWXYZ", "abcdefghijklmnopqrstuvwxyz"), "pmid")
    ]'
  )
  #If no Pubmed-related nodes exist then return NA
  if (length(pubmed_nodes) == 0) {
    return(NA_character_)
  }
  
  #Extracting text from the matching nodes:
  pubmed_text <- xml_text(pubmed_nodes, trim = TRUE)
  #Pulling out digit sequences from the text:
  pubmed_ids <- unlist(
    regmatches(
      pubmed_text,
      gregexpr('[0-9]+', pubmed_text)
    )
  )
  #Removing the duplicates:
  pubmed_ids <- unique(pubmed_ids)
  
  #If no numeric Ids are found then return NA:
  if(length(pubmed_ids) == 0) {
    return(NA_character_)
  }
  #Collapse multiple Pubmed Ids into a sing semicolon-separated string:
  return(paste(pubmed_ids, collapse = "; "))
}

# Helper function: Extract one Docsum record into a Data.table row:

extract_docsum_record <- function(record) {
  id_node <- xml_find_first(record, "./Id")
  
  id_value <- if (inherits(id_node, "xml_missing")) {
    NA_character_
  } else {
    xml_text(id_node, trim = TRUE)
  }
  
  # Build one row of metadata.
  result <- data.table(
    ID = id_value, #NCBI UID for GEO record
    SeriesAccession = get_first_item(               #Using get_first_item to try multiple possible names
      record,
      c("Accession", "GSE", "SeriesAccession")
    ),
    # Study title
    Title = get_first_item(
      record,
      c("title", "Title")
    ),
    #Study Summary
    Summary = get_first_item(
      record,
      c("summary", "Summary")
    ),
    #Organism
    Organism = get_first_item(
      record,
      c("taxon", "Organism")
    ),
    #Geodata type (ie. Expression values/methylation profiles)
    Type = get_first_item(
      record,
      c("gdsType", "Type")
    ),
    #GPL ids
    Platform = get_first_item(
      record,
      c("GPL", "Platform")
    ),
    #Number of samples in the record
    SampleCount = suppressWarnings(
      as.integer(
        get_first_item(
          record,
          c("n_samples", "SampleCount")
        )
      )
    ),
    #Supplementary info if provided:
    DownloadData = get_first_item(
      record,
      c("suppFile", "SupplementaryFiles", "DownloadData")
    ),
    #Pubmed IDs linked to the record:
    PubMedIDs = get_pubmed_ids(record)
  )
  #Return the one-row data table:
  return(result)
  
}

#Function #2, Extracting GEO metadata for a batch of NCBI UIDs:

extract_geo_records_batch <- function(uid.vec, api.key = "") {
  #Convert Ids to character values in case they came in as numbers:
  uid.vec <- as.character(uid.vec)
  #Remove missing/blank Ids
  uid.vec <- uid.vec[!is.na(uid.vec) & nzchar(uid.vec)]
  #IF no valid ids are available then return an empty table:
  if(length(uid.vec) == 0) {
    return(
      data.table(
        ID = character(),
        SeriesAccession = character(),
        Title = character(),
        Summary = character(),
        Organism = character(),
        Type = character(),
        Platform = character(),
        SampleCount = integer(),
        DownloadData = character(),
        PubMedIDs = character()
      )
    )
  }
  #Combine IDs into a comma separated list for a ESummary request (Dora suggested this)
  uid.string <- paste(uid.vec, collapse = ",")
  url <- paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?",
    "db=gds",
    "&id=", uid.string,
    "&retmode=xml"
  )
  #Add API key if present:
  url <- add_api_key(url, api.key)
  #Read the XML result:
  xml_doc <- safe_read_xml(url)
  
  #If the XML request fails, return one missing row per requested ID:
  if(is.null(xml_doc)) {
    return(
      data.table(
        ID = uid.vec,
        SeriesAccession = NA_character_,
        Title = NA_character_,
        Summary = NA_character_,
        Organism = NA_character_,
        Type = NA_character_,
        Platform = NA_character_,
        SampleCount = NA_integer_,
        DownloadData = NA_character_,
        PubMedIDs = NA_character_
      )
    )
  }
  # Storing each GEO summary record as a Docsum node:
  records <- xml_find_all(xml_doc, ".//DocSum")
  
  #If NCBI returned no records then return the missing rows:
  if (length(records) == 0) {
    return(
      data.table(
        ID = uid.vec,
        SeriesAccession = NA_character_,
        Title = NA_character_,
        Summary = NA_character_,
        Organism = NA_character_,
        Type = NA_character_,
        Platform = NA_character_,
        SampleCount = NA_integer_,
        DownloadData = NA_character_,
        PubMedIDs = NA_character_
      )
    )
  }
  #Convert each XML DocSum record into a single datatable row:
  record_list <- lapply(records, extract_docsum_record)
  #Stacking all the data.tables into one table:
  result <- rbindlist(record_list, fill = TRUE)
  
  return(result)
}

###Function #3: Extract GEO metadata for many IDs in batches:

extract_geo_records <- function(uid.vec, batch.size = 200, delay.seconds = 0.4, api.key = "") {
  #Converting IDs to unique character values:
  uid.vec <- unique(as.character(uid.vec))
  #Removing the missing or blank Ids:
  uid.vec <- uid.vec[!is.na(uid.vec) & nzchar(uid.vec)]
  
  #Split Ids into batches:
  batches <- split(
    uid.vec,
    ceiling(seq_along(uid.vec) / batch.size)
  )
  #Create an empty list to store results from each batch:
  batch_results <- vector("list", length(batches))
  #loop over batches:
  for (i in seq_along(batches)) {
    #Progress message for clarity:
    message("Getting GEO metadata batch ", i, " of ", length(batches))
    
    #Extract metadata for this batch:
    batch_results[[i]] <- extract_geo_records_batch(
      uid.vec = batches[[i]],
      api.key = api.key
    )
    #Courtesy pause to not ban St Jude again:
    Sys.sleep(delay.seconds)
  }
  #Combine all batch results into a single table:
  result <- rbindlist(batch_results, fill = TRUE)
  #Return the finished and combined metadata table:
  return(result)
}

### Another Helper Function: finding keyword hits inside text:

find_hits <- function(text, terms) {
  #Combine possible multiple text fields into a single searchable string:
  text <- paste(text, collapse = " ")
  
  #For each term, checking if it appears in the text:
  hits <- terms[
    vapply(
      terms,
      function(term) {
        grepl(term, text, ignore.case = TRUE, perl = TRUE)
      },
      logical(1)
    )
  ]
  
  #If there are no hits, then return NA:
  if(length(hits) == 0) {
    return(NA_character_)
  }
  #Otherwise return the unique hits semicolon separated:
  return(paste(unique(hits), collapse = "; "))
}

### Another Helper Function to classify one GEO record:

classify_one_record <- function(title, summary, type) {
  #Put the searchable record text together:
  text <- paste(title, summary, type, sep = " ")
  
  #Terms suggesting a clinical trial or clinical cohort:
  clinical.terms <- c(            # I used AI to come up with these terms
    "clinical trial",
    "phase I",
    "phase II",
    "phase III",
    "NCT[0-9]+",
    "AAML[0-9]+",
    "COG",
    "HOVON",
    "trial",
    "protocol",
    "overall survival",
    "event-free survival",
    "induction",
    "remission",
    "NCT id"
  )
  
  # Terms suggesting a control condition.
  control.terms <- c(
    "DMSO",
    "vehicle",
    "untreated",
    "mock",
    "placebo",
    "control"
  )
  
  # Terms suggesting cell-line samples.
  cell.line.terms <- c(
    "cell line",
    "MOLM",
    "MOLM13",
    "MV4",
    "MV4-11",
    "MV4;11",
    "THP-1",
    "HL-60",
    "HL60",
    "K562",
    "OCI-AML",
    "KG-1",
    "KG1",
    "Kasumi",
    "NB4"
  )
  
  # Terms suggesting patient derived material.
  patient.terms <- c(
    "patient",
    "primary AML",
    "primary sample",
    "bone marrow",
    "peripheral blood",
    "blast"
  )
  #Finding hits for each category:
  clinical.hits <- find_hits(text, clinical.terms)
  control.hits <- find_hits(text, control.terms)
  cell.hits <- find_hits(text, cell.line.terms)
  patient.hits <- find_hits(text, patient.terms)
  #Determine sample type:
  sample.type.guess <- if (!is.na(cell.hits) && !is.na(patient.hits)) {
    "Mixed/Unclear"
  } else if (!is.na(cell.hits)) {
    "Cell line"
  } else if (!is.na(patient.hits)) {
    "Patient/sample"
  } else {
    "Unknown"
  }
  
  #Determine the broad study group:
  group.guess <- if (!is.na(clinical.hits)) {
    "Clinical trial/cohort"
  } else if (!is.na(control.hits) || !is.na(cell.hits) || !is.na(patient.hits)) {
    "Cell line/patient sample experiment"
  } else {
    "Needs manual review"
  }
  
  #REturn one row explaining the classification:
  return(
    data.table(
      GroupGuess = group.guess,
      SampleTypeGuess = sample.type.guess,
      ClinicalHits = clinical.hits,
      ControlHits = control.hits,
      CellLineHits = cell.hits,
      PatientHits = patient.hits
    )
  )
}

###Function #4, Classifying a GEO metadata table

classify_geo_table <- function(geo.dt) {
  #If the table has no rows, return it unchanged"
  if (nrow(geo.dt) == 0) {
    return (geo.dt)
  }
  #Apply the classify_one_record() function to each row:
  class_list <- lapply(
    seq_len(nrow(geo.dt)),
    function(i) {
      classify_one_record(
        title = geo.dt$Title[i],
        summary = geo.dt$Summary[i],
        type = geo.dt$Type[i]
      )
    }
  )
  #Stack the classification rows into one table:
  
  class_dt <- rbindlist(class_list, fill = TRUE)
  #Bind the metadata and classification colulmsn together:
  result <- cbind(geo.dt, class_dt)
  return(result)
}

#####!!! Main Function: Search the GEO database and build candidate workbook data:

get_GEO_candidates <- function(query.strings, retmax = 9999, batch.size = 200, delay.seconds = 0.4, api.key = "") {
  ##Step 1. run ESearch for each query string
  
  search_results <- lapply(
    query.strings,
    function(q){
      get_gds_ids(
        query.string = q,
        retmax = retmax,
        api.key = api.key
      )
    }
  )
  #Combine all search result tables:
  query_dt <- rbindlist(search_results, fill = TRUE)
  #Remove the duplicate query/ID rows:
  query_dt <- unique(query_dt)
  #Message for clarity:
  message("Unique GEO/GDS UIDs found: ", length(unique(query_dt$ID)))
  
  ##STEP 2: Summarize which query strings found each ID:
  
  #For each ID, combine all the queries:
  query_by_id <- query_dt[
    ,
    .(
      QueryString = paste(unique(query_string), collapse = "; ")
    ), 
    by = ID
  ]
  
  ##STEP 3: Use Esummary XML to get metadata:
  
  geo_dt <- extract_geo_records(
    uid.vec = unique(query_dt$ID),
    batch.size = batch.size,
    delay.seconds = delay.seconds,
    api.key = api.key
  )
  
  ##STEP 4: Attach the query-string info
  
  #Merging query strings onto the GEO metadata table we made:
  geo_dt <- merge(
    geo_dt,
    query_by_id,
    by = "ID",
    all.x = TRUE
  )
  
  ##STEP 5: classify studies using the metadata text
  geo_dt <- classify_geo_table(geo_dt)
  ##STEP 6: create useful filtered sheets:
  
  # Expression/Methylation records:
  human_expression_methylation <- geo_dt[
    grepl("Homo sapiens|human", Organism, ignore.case = TRUE) &
      grepl("Expression profiling|Methylation profiling", Type, ignore.case = TRUE)
  ]
  # Controlled Experiments/Patient samples
  controlled_exp <- human_expression_methylation[
    GroupGuess == "Cell line/patient sample experiment"
  ]
  #Clinical trial/cohort datasets:
  clinical_cohort <- human_expression_methylation[
    GroupGuess == "Clinical trial/cohort"
  ]
  # Records that would need manual reviewing:
  needs_review <- human_expression_methylation[
    GroupGuess == "Needs manual review"
  ]
  ### Step 7: return everything as a list:
  
  return(
    list(
      all_records = geo_dt,
      human_expression_methylation = human_expression_methylation,
      controlled_exp = controlled_exp,
      clinical_cohort = clinical_cohort,
      needs_review = needs_review,
      query_log = query_dt
    )
  )
}

############### Running and testing the workflow::

geo_results <- get_GEO_candidates(
  query.strings = query.strings,
  retmax = 9999,
  batch.size = 200,
  delay.seconds = 0.4,
  api.key = ncbi.api.key
)

#Write the results to one Excel file:
write_xlsx(geo_results, out.path)





"#TODO maybe find pubmed ids as well

# ----------------------------------------------------------------------------------------------------------



#############################################
# Functions to get information on GEO study data sets
# currently written as a script
# convert to a function that accepts query string like "AML+DMSO" as input
# has a lot of webscraping of HTML/XML, string operations, etc
# we'd like also capture pubmed ID, title, summary
# explore using XML R package to make webscraping easier

# E-utils link for the GEO GDS ids
gds.ids.link="https://eutils.ncbi.nlm.nih.gov/eutils/esearch.fcgi?db=gds&term=AML+DMSO&rettype=uilist&retmax=9999"

# get the GEO data set GDS IDs
gds.ids=scan(gds.ids.link,what=character(),sep="\n")
gds.ids=grep("<Id>",gds.ids,value=T,fixed=T)
gds.ids=gsub("\t<Id>","",gds.ids,fixed=T)
gds.ids=gsub("</Id>","",gds.ids,fixed=T)



library(xml2)
library(data.table)

extract_geo_record = function(uid) {
  
  url = paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
    "esummary.fcgi?db=gds&id=", uid, "&retmode=xml"
  )
  
  xml_doc = read_xml(url)
  records = xml_find_all(xml_doc, ".//DocSum")
  
  # Return an empty/error row if no records were found
  if (length(records) == 0) {
    return(
      data.table(
        ID = as.character(uid),
        SeriesAccession = NA_character_,
        Title = NA_character_,
        Summary = NA_character_,
        Organism = NA_character_,
        Type = NA_character_,
        Platform = NA_character_,
        SampleCount = NA_integer_,
        DownloadData = NA_character_
      )
    )
  }
  
  record_list = lapply(records, function(record) {
    
    get_item = function(name) {
      node = xml_find_first(
        record,
        sprintf('.//Item[@Name="%s"]', name)
      )
      
      if (inherits(node, "xml_missing")) {
        NA_character_
      } else {
        xml_text(node, trim = TRUE)
      }
    }
    
    data.table(
      ID = xml_text(
        xml_find_first(record, "./Id"),
        trim = TRUE
      ),
      SeriesAccession = get_item("Accession"),
      Title = get_item("title"),
      Summary = get_item("summary"),
      Organism = get_item("taxon"),
      Type = get_item("gdsType"),
      Platform = get_item("GPL"),
      SampleCount = as.integer(get_item("n_samples")),
      DownloadData = get_item("suppFile")
    )
  })
  
  rbindlist(record_list, fill = TRUE)
}


geo_dt = rbindlist(
  lapply(gds.ids[1:10], function(uid) {
    
    result = extract_geo_record(uid)
    
    # Avoid sending requests too quickly
    Sys.sleep(0.4)
    
    result
  }),
  fill = TRUE
)

