#' Identifies all R / L interactions
#'
#' This function will map all RL interactions
#'
#' @param input the input graph analysis object
#' @param select the communities to select
#' @param expand if true will expand to grab all communities that share an edge with the input community
#' @export
#' @details
#' This will use the calc_agg_bulk results to ID networks
#' @examples
#' ex_sc_example <- id_rl(input = ex_sc_example)

extract_communities <- function(input, select, expand = TRUE){

  members <- input$Clusters_Results$membership
  keep <- c()
  for (i in 1:length(select)) {
    int <- select[i]
    ind <- which(members == int)
    keep <- c(keep, ind)
  }

  if(expand == TRUE){
    keep2 <- unlist(adjacent_vertices(input$igraph_Network, keep, mode = "all"))
    select2 <- unique(members[keep2])
    keep <- c()
    for (i in 1:length(select2)) {
      int <- select2[i]
      ind <- which(members == int)
      keep <- c(keep, ind)
    }
  }

  subgraph <- induced_subgraph(input$igraph_Network, keep)
  ind2 <- match(names(V(subgraph)), names(V(input$igraph_Network)))
  l  <- input$layout[ind2,]
  results <- vector(mode = "list", length = 2)
  results[[1]] <- subgraph
  results[[2]] <- l
  names(results) <- c("igraph_Network", "layout")
  return(results)

}

