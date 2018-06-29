#' check_package
#'
#' This will build documentation for the package without errors due to new data file
#'
#' @export


check_package <- function(){
  devtools::load_all()
  roxygen2::roxygenise()
}


