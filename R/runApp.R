#' Run Shiny App
#'
#' @export
runShinyApp <- function() {
  appDir <- system.file("shinyapp", package = "PepFlow")
  if (appDir == "") {
    stop("The App couldn´t be located. Please install the package another time.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}
