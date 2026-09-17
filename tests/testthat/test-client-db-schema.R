test_that("buildDb samples table uses per-sample columns when inventory CSV exists", {
  repo_root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
  python_bin <- if (file.exists(file.path(repo_root, "parsing", ".venv", "bin", "python3"))) {
    file.path(repo_root, "parsing", ".venv", "bin", "python3")
  } else {
    Sys.which("python3")
  }
  if (!nzchar(python_bin)) {
    testthat::skip("python3 not available")
  }

  toy_matrix <- file.path(repo_root, "Toy-Datasets", "GSE128103_series_matrix.txt.gz")
  if (!file.exists(toy_matrix)) {
    testthat::skip("Toy-Datasets matrix fixture missing")
  }

  work <- tempfile("geo_build_test_")
  dir.create(work)
  matrices_dir <- file.path(work, "matrices")
  dir.create(matrices_dir)
  file.copy(toy_matrix, file.path(matrices_dir, basename(toy_matrix)))

  writeLines(
    c(
      '"SeriesAccession","Platform","Status"',
      '"GSE99999","GPL00000","downloaded"'
    ),
    file.path(work, "downloaded_matrices.csv")
  )

  db_file <- file.path(work, "test.duckdb")
  build_script <- file.path(repo_root, "parsing", "buildDb.py")
  result <- system2(
    python_bin,
    c(build_script, matrices_dir, "--db-path", db_file),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(result, "status")
  if (!is.null(status) && status != 0L) {
    testthat::fail(paste(c(result, collapse = "\n")))
  }
  if (!file.exists(db_file)) {
    testthat::fail("buildDb did not create database file")
  }

  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = db_file, read_only = TRUE)
  sample_cols <- DBI::dbGetQuery(con, "DESCRIBE samples")$column_name
  expect_true("series_accession" %in% sample_cols)
  expect_true("sample_geo_accession" %in% sample_cols)
  expect_true("sample_characteristics_ch1" %in% sample_cols)
  expect_true("sample_molecule_ch1" %in% sample_cols)
  expect_false("SeriesAccession" %in% sample_cols)

  study_tables <- DBI::dbGetQuery(
    con,
    "SELECT table_name FROM information_schema.tables WHERE table_schema = 'main'"
  )$table_name
  expect_true("studies" %in% study_tables)

  DBI::dbDisconnect(con, shutdown = TRUE)
  unlink(work, recursive = TRUE)
})
