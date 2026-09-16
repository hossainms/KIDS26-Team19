## Shiny app: fetch GEO matrices for the selected diagnosis, rebuild data/geo.duckdb
## from scratch (parsing/buildDb.py), then browse the resulting sample metadata.
## R setup:      install.packages(c("shiny", "DT", "DBI", "duckdb"))
## Python setup: pip install -r parsing/requirements.txt

library(shiny)
library(DT)
library(DBI)
library(duckdb)

#Locate the repo root regardless of whether the app is launched from "." or "client/":
find_repo_root <- function() {
  for (candidate in c(".", "..")) {
    if (file.exists(file.path(candidate, "parsing", "buildDb.py"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }
  stop("Cannot locate repository root (expected parsing/buildDb.py).")
}

#Prefer the parsing/.venv Python if it exists, otherwise fall back to python3 on PATH:
find_python <- function(repo_root) {
  venv_python <- file.path(repo_root, "parsing", ".venv", "bin", "python3")
  if (file.exists(venv_python)) venv_python else "python3"
}

repo_root <- find_repo_root()
python_bin <- find_python(repo_root)
db_path <- file.path(repo_root, "data", "geo.duckdb")
source(file.path(repo_root, "client", "app_logic.R"))

#Diagnosis abbreviations supported by the dropdown, mapped to full names:
diagnosis_map <- c(
  AML = "acute myeloid leukemia",
  ALL = "acute lymphocytic leukemia"
)

#Re-download matrices for a diagnosis and rebuild data/geo.duckdb from scratch:
rebuild_database <- function(diagnosis) {
  out_dir <- file.path(repo_root, "downloads", paste0("geo_", tolower(diagnosis)))
  fetch_result <- run_command_capture(
    "Rscript",
    c(shQuote(file.path(repo_root, "R_Scripts", "05_Download_Metadata_Inventory.R")),
      "--diagnosis", diagnosis, "--out", shQuote(out_dir))
  )
  if (!identical(fetch_result$status, 0L)) {
    stop(
      "Fetching GEO matrices for '", diagnosis, "' failed (exit code ", fetch_result$status,
      ").\n", fetch_result$output
    )
  }

  if (file.exists(db_path)) unlink(db_path) #destroy any previous diagnosis's database

  build_result <- run_command_capture(
    python_bin,
    c(shQuote(file.path(repo_root, "parsing", "buildDb.py")),
      shQuote(file.path(out_dir, "matrices")), "--db-path", shQuote(db_path))
  )
  if (!identical(build_result$status, 0L)) {
    stop(
      "Building the DuckDB for '", diagnosis, "' failed (exit code ", build_result$status,
      ").\n", build_result$output
    )
  }

  list(
    db_path = db_path,
    output = c("R download output:", fetch_result$output, "", "Python build output:", build_result$output)
  )
}

ui <- fluidPage(
  titlePanel("GEO Sample Explorer"),
  sidebarLayout(
    sidebarPanel(
      selectInput("diagnosis", "Diagnosis", choices = names(diagnosis_map), selected = "AML"),
      helpText("Changing the diagnosis re-downloads GEO matrices and rebuilds the database.")
    ),
    mainPanel(
      verbatimTextOutput("message"),
      DTOutput("results")
    )
  )
)

server <- function(input, output, session) {
  con <- reactiveVal(NULL)
  status <- reactiveVal("")

  rebuild <- function(diagnosis) {
    old_con <- con()
    if (!is.null(old_con)) dbDisconnect(old_con, shutdown = TRUE)
    con(NULL)
    status(paste0("Loading '", diagnosis, "'..."))
    outcome <- tryCatch({
      rebuild_result <- rebuild_database(diagnosis)
      list(ok = TRUE, output = rebuild_result$output)
    }, error = function(e) list(ok = FALSE, error = conditionMessage(e)))
    if (outcome$ok) {
      con(dbConnect(duckdb::duckdb(), dbdir = db_path, read_only = TRUE))
      status(paste(outcome$output, collapse = "\n"))
    } else {
      status(paste0("Failed to load '", diagnosis, "': ", outcome$error))
    }
  }

  observeEvent(input$diagnosis, rebuild(input$diagnosis), ignoreNULL = FALSE)

  session$onSessionEnded(function() {
    active <- isolate(con())
    if (!is.null(active)) dbDisconnect(active, shutdown = TRUE)
  })

  results <- reactive({
    active <- con()
    req(active)
    dbGetQuery(active, "SELECT * FROM samples ORDER BY series_accession, sample_geo_accession")
  })

  output$message <- renderText({ status() })

  output$results <- renderDT({
    req(con())
    datatable(results(), rownames = FALSE, filter = "top",
              options = list(pageLength = 25, scrollX = TRUE))
  })
}

shinyApp(ui, server)

shinyApp(ui, server)
