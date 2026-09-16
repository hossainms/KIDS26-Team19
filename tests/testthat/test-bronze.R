source(file.path("..", "..", "R", "bronze.R"))

test_that("series folder drops the last three digits", {
  expect_equal(gse_folder("GSE100446"), "GSE100nnn")
  expect_equal(gse_folder("GSE982"), "GSEnnn")
  expect_equal(gse_folder("GSE14"), "GSEnnn")
})

test_that("matrix listing finds single and multi-platform files", {
  html <- c(
    '<a href="GSE995-GPL80_series_matrix.txt.gz">GSE995-GPL80_series_matrix.txt.gz</a>',
    '<a href="GSE995-GPL96_series_matrix.txt.gz">GSE995-GPL96_series_matrix.txt.gz</a>',
    '<a href="/geo/series/GSEnnn/">Parent</a>'
  )
  expect_equal(
    parse_matrix_listing(html),
    c("GSE995-GPL80_series_matrix.txt.gz", "GSE995-GPL96_series_matrix.txt.gz")
  )
})

test_that("platform comes from the file name or is blank", {
  expect_equal(platform_from_file("GSE995-GPL96_series_matrix.txt.gz"), "GPL96")
  expect_equal(platform_from_file("GSE982_series_matrix.txt.gz"), NA_character_)
})
