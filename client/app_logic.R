run_command_capture <- function(command, args) {
  output <- suppressWarnings(system2(command, args, stdout = TRUE, stderr = TRUE))
  command_status <- attr(output, "status")

  list(
    status = if (is.null(command_status)) 0L else as.integer(command_status),
    output = paste(output, collapse = "\n")
  )
}