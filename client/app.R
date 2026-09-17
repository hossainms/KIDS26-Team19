## Shiny app: browse sample metadata in data/geo.duckdb built upstream via parsing/buildDb.py.
## R setup:      install.packages(c("shiny", "DT", "DBI", "duckdb"))
## Python setup: pip install -r parsing/requirements.txt (for buildDb.py outside the app)

library(shiny)
library(DT)
library(DBI)
library(duckdb)

find_repo_root <- function() {
  for (candidate in c(".", "..")) {
    if (file.exists(file.path(candidate, "parsing", "buildDb.py"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }
  stop("Cannot locate repository root (expected parsing/buildDb.py).")
}

find_python <- function(repo_root) {
  venv_python <- file.path(repo_root, "parsing", ".venv", "bin", "python3")
  if (file.exists(venv_python)) venv_python else "python3"
}

repo_root <- find_repo_root()
python_bin <- find_python(repo_root)
db_path <- file.path(repo_root, "data", "geo.duckdb")
source(file.path(repo_root, "client", "app_logic.R"))

diagnosis_map <- c(
  AML = "acute myeloid leukemia",
  ALL = "acute lymphocytic leukemia"
)

ui <- fluidPage(
  titlePanel("GEO Sample Explorer"),
  sidebarLayout(
    sidebarPanel(
      selectInput("diagnosis", "Diagnosis", choices = names(diagnosis_map), selected = "AML"),
      helpText(
        "Matrices and DuckDB are prepared outside this app ",
        "(R_Scripts/05_Download_Metadata_Inventory.R and parsing/buildDb.py). ",
        "Startup opens data/geo.duckdb if present. Changing diagnosis rebuilds from that diagnosis's cached matrices."
      )
    ),
    mainPanel(
      h4(textOutput("state", inline = TRUE)),
      verbatimTextOutput("message"),
      DTOutput("results")
    )
  )
)

server <- function(input, output, session) {
  con <- reactiveVal(NULL)
  status <- reactiveVal("Ready.")
  query_error <- reactiveVal("")
  loading <- reactiveVal(FALSE)
  db_env <- new.env(parent = emptyenv())
  db_env$connection <- NULL

  disconnect_db <- function() {
    active <- db_env$connection
    if (!is.null(active)) {
      dbDisconnect(active, shutdown = TRUE)
    }
    db_env$connection <- NULL
    con(NULL)
  }

  attach_connection <- function(connection, message_text) {
    db_env$connection <- connection
    con(connection)
    query_error("")
    status(message_text)
    loading(FALSE)
  }

  open_existing_database <- function() {
    existing <- connect_existing_database(db_path)
    if (is.null(existing)) {
      return(FALSE)
    }
    attach_connection(
      existing,
      paste0("Opened existing database: ", db_path)
    )
    TRUE
  }

  rebuild_diagnosis_database <- function(diagnosis) {
    loading(TRUE)
    disconnect_db()
    query_error("")
    status(load_status_message(diagnosis, repo_root))

    outcome <- tryCatch(
      {
        result <- build_database_from_cache(diagnosis, repo_root, python_bin, db_path)
        list(ok = TRUE, output = result$output)
      },
      error = function(e) list(ok = FALSE, error = conditionMessage(e))
    )

    if (outcome$ok) {
      attach_connection(
        connect_existing_database(db_path),
        paste(outcome$output, collapse = "\n")
      )
    } else {
      loading(FALSE)
      status(paste0("Failed to load '", diagnosis, "': ", outcome$error))
    }
  }

  diagnosis_initialized <- reactiveVal(FALSE)

  observeEvent(input$diagnosis, {
    if (!diagnosis_initialized()) {
      diagnosis_initialized(TRUE)
      if (open_existing_database()) {
        return()
      }
      rebuild_diagnosis_database(input$diagnosis)
      return()
    }
    rebuild_diagnosis_database(input$diagnosis)
  }, ignoreInit = FALSE)

  session$onSessionEnded(function() {
    active <- db_env$connection
    if (!is.null(active)) {
      dbDisconnect(active, shutdown = TRUE)
    }
    db_env$connection <- NULL
  })

  results <- reactive({
    active <- con()
    req(active)
    tryCatch(
      fetch_samples_table(active),
      error = function(e) {
        query_error(conditionMessage(e))
        NULL
      }
    )
  })

  output$state <- renderText({
    if (isTRUE(isolate(loading()))) {
      return(isolate(status()))
    }
    if (nzchar(isolate(query_error()))) {
      return(paste0("Database query error: ", isolate(query_error())))
    }
    if (is.null(isolate(con()))) {
      return(isolate(status()))
    }
    "Database loaded. Browse samples below."
  })

  output$message <- renderText({ status() })

  output$results <- renderDT({
    tbl <- results()
    query_err <- query_error()
    if (!is.null(query_err) && nzchar(query_err)) {
      validate(need(FALSE, paste0("Could not read samples table:\n", query_err)))
    }
    validate(need(!is.null(tbl), "No sample rows to display yet."))
    characteristics_idx <- which(names(tbl) == "sample_characteristics_ch1") - 1L
    column_defs <- list()
    if (length(characteristics_idx) == 1L && characteristics_idx >= 0L) {
      column_defs <- list(
        list(width = "320px", targets = characteristics_idx)
      )
    }
    datatable(
      tbl,
      rownames = FALSE,
      filter = "top",
      colnames = samples_table_column_labels,
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        columnDefs = column_defs
      )
    )
  })
}

shinyApp(ui, server)
