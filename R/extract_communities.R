#' Identifies all R / L interactions
#'
#' This function will map all RL interactions
#'
#' @param input the input graph analysis object
#' @param select the communities to select
#' @export
#' @details
#' This will use the calc_agg_bulk results to ID networks
#' @examples
#' ex_sc_example <- id_rl(input = ex_sc_example)

extract_communities <- function(input, select){
  members <- input$Communities$membership
  keep <- c()
  for (i in 1:length(select)) {
    int <- select[i]
    ind <- which(members == int)
    keep <- c(keep, ind)
  }
  subgraph <- induced_subgraph(input$input, keep)
  return(subgraph)
}

