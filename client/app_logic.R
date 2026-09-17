run_command_capture <- function(command, args) {
  output <- suppressWarnings(system2(command, args, stdout = TRUE, stderr = TRUE))
  command_status <- attr(output, "status")

  list(
    status = if (is.null(command_status)) 0L else as.integer(command_status),
    output = paste(output, collapse = "\n")
  )
}

diagnosis_out_dir <- function(diagnosis, repo_root) {
  file.path(repo_root, "downloads", paste0("geo_", tolower(diagnosis)))
}

diagnosis_matrices_dir <- function(diagnosis, repo_root) {
  file.path(diagnosis_out_dir(diagnosis, repo_root), "matrices")
}

has_cached_matrices <- function(diagnosis, repo_root) {
  matrices_dir <- diagnosis_matrices_dir(diagnosis, repo_root)
  if (!dir.exists(matrices_dir)) {
    return(FALSE)
  }
  length(list.files(matrices_dir, pattern = "_series_matrix\\.txt\\.gz$", full.names = FALSE)) > 0L
}

run_python_build_db <- function(repo_root, python_bin, matrices_dir, db_path) {
  if (file.exists(db_path)) {
    unlink(db_path)
  }
  run_command_capture(
    python_bin,
    c(
      file.path(repo_root, "parsing", "buildDb.py"),
      matrices_dir,
      "--db-path", db_path
    )
  )
}

load_status_message <- function(diagnosis, repo_root) {
  if (has_cached_matrices(diagnosis, repo_root)) {
    return(paste0("Loading '", diagnosis, "' from cached matrices (buildDb.py)..."))
  }
  paste0(
    "No matrices under downloads/geo_", tolower(diagnosis), "/matrices. ",
    "Run parsing/buildDb.py upstream after downloading GEO data."
  )
}

#' Build data/geo.duckdb from cached matrices for a diagnosis (no GEO download in the app).
build_database_from_cache <- function(diagnosis, repo_root, python_bin, db_path) {
  matrices_dir <- diagnosis_matrices_dir(diagnosis, repo_root)
  cached <- has_cached_matrices(diagnosis, repo_root)

  if (!cached) {
    stop(
      "No cached matrices for '", diagnosis, "' in ", matrices_dir, ". ",
      "Download GEO matrices and run parsing/buildDb.py outside this app."
    )
  }

  build_result <- run_python_build_db(repo_root, python_bin, matrices_dir, db_path)
  if (!identical(build_result$status, 0L)) {
    stop(
      "Building the DuckDB for '", diagnosis, "' failed (exit code ", build_result$status,
      ").\n", build_result$output
    )
  }

  list(
    db_path = db_path,
    output = c(
      "Built from cached matrices:",
      matrices_dir,
      "",
      "Python build output:",
      build_result$output
    )
  )
}

SAMPLES_TABLE_COLUMNS <- c(
  "series_accession",
  "series_platform_id",
  "sample_geo_accession",
  "sample_organism_ch1",
  "sample_data_row_count",
  "sample_characteristics_ch1",
  "sample_molecule_ch1"
)

SAMPLES_TABLE_QUERY <- paste(
  "SELECT",
  paste(SAMPLES_TABLE_COLUMNS, collapse = ", "),
  "FROM samples",
  "ORDER BY series_accession, sample_geo_accession"
)

samples_table_column_labels <- c(
  "Series",
  "Platform",
  "Sample",
  "Organism",
  "Row count",
  "Sample characteristics",
  "Sample molecule"
)

fetch_samples_table <- function(connection) {
  DBI::dbGetQuery(connection, SAMPLES_TABLE_QUERY)
}

#' Open an existing DuckDB file read-only if it has a samples table.
connect_existing_database <- function(db_path) {
  if (!file.exists(db_path)) {
    return(NULL)
  }
  con <- tryCatch(
    DBI::dbConnect(duckdb::duckdb(), dbdir = db_path, read_only = TRUE),
    error = function(e) NULL
  )
  if (is.null(con)) {
    return(NULL)
  }
  ok <- tryCatch({
    DBI::dbGetQuery(con, "SELECT 1 FROM samples LIMIT 1")
    TRUE
  }, error = function(e) FALSE)
  if (!ok) {
    DBI::dbDisconnect(con, shutdown = TRUE)
    return(NULL)
  }
  con
}
