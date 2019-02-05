#' Launch Chordomics
#' 
#' Function launches Chordomics server and UI
#' @export
#' 

launchApp <- function(){
  shiny::shinyApp(server = shinyAppServer, ui = shinyAppUI)
}
