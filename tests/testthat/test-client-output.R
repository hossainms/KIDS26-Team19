source(file.path("client", "app_logic.R"))

test_that("run_command_capture returns combined command output", {
  result <- run_command_capture(
    "Rscript",
    c("-e", shQuote("cat('R output\\n'); message('R diagnostic')"))
  )

  expect_identical(result$status, 0L)
  expect_match(result$output, "R output")
  expect_match(result$output, "R diagnostic")
})