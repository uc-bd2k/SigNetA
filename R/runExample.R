#' @export
#' 
#' @importFrom shiny runApp
runExample <- function() {
  appDir <- system.file("examples", "SigNetA-app", package = "SigNetA")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}