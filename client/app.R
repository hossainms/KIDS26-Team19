## Baseline Shiny app: query the DuckDB sample metadata by diagnosis.
## Build the database first:  python parsing/buildDb.py
## One-time setup: install.packages(c("shiny", "DT", "DBI", "duckdb"))

library(shiny)
library(DT)
library(DBI)
library(duckdb)

db_path <- file.path("..", "data", "geo.duckdb")
if (!file.exists(db_path)) {
  db_path <- file.path("data", "geo.duckdb")
}
db_path <- normalizePath(db_path, mustWork = TRUE)

con <- dbConnect(duckdb::duckdb(), dbdir = db_path, read_only = TRUE)
onStop(function() dbDisconnect(con, shutdown = TRUE))

ui <- fluidPage(
  titlePanel("GEO Sample Explorer"),
  sidebarLayout(
    sidebarPanel(
      textInput("diagnosis", "Diagnosis", value = "AML"),
      helpText("Only 'AML' is supported for now."),
      actionButton("search", "Search", class = "btn-primary")
    ),
    mainPanel(
      textOutput("message"),
      DTOutput("results")
    )
  )
)

server <- function(input, output, session) {
  ## TODO: two upstream fetch requests (supplied by another group) will populate
  ## the database before this query runs.
  results <- eventReactive(input$search, ignoreNULL = FALSE, {
    if (!identical(toupper(trimws(input$diagnosis)), "AML")) {
      return(NULL)
    }
    dbGetQuery(
      con,
      "SELECT * FROM samples ORDER BY series_accession, sample_geo_accession"
    )
  })

  output$message <- renderText({
    if (is.null(results())) "No data for that diagnosis. Try 'AML'." else ""
  })

  output$results <- renderDT({
    req(results())
    datatable(results(), rownames = FALSE, filter = "top",
              options = list(pageLength = 25, scrollX = TRUE))
  })
}

shinyApp(ui, server)
