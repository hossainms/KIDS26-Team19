REQUIRED <- c("date", "time", "time_zone", "person", "coding_agent", "model", "machine", "branch", "commits")

test_that("every log has a Central time name and full frontmatter", {
  logs <- list.files(file.path("..", "..", "logs"), pattern = "\\.md$", full.names = TRUE)
  expect_gt(length(logs), 0)
  for (log in logs) {
    expect_match(basename(log), "^\\d{4}-\\d{2}-\\d{2}-\\d{4}-[a-z0-9-]+\\.md$", info = log)
    lines <- readLines(log, warn = FALSE)
    ends <- which(lines == "---")
    expect_true(length(ends) >= 2 && ends[1] == 1, info = log)
    keys <- sub(":.*", "", lines[(ends[1] + 1):(ends[2] - 1)])
    expect_true(all(REQUIRED %in% keys), info = paste(log, "missing", paste(setdiff(REQUIRED, keys), collapse = ", ")))
    expect_true(any(grepl("^time_zone: America/Chicago", lines)), info = log)
  }
})
