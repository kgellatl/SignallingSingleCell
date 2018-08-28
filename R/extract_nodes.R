#' Identifies all R / L interactions
#'
#' This function will map all RL interactions
#'
#' @param input the input graph analysis object
#' @param nodes the nodes to extract
#' @param expand if true will expand to grab all nodes within distance n of a node of interest
#' @export
#' @details
#' This will use the calc_agg_bulk results to ID networks
#' @examples
#' ex_sc_example <- id_rl(input = ex_sc_example)

extract_nodes <- function(input, nodes, expand = 0){
  path <- c()
  for (i in 1:length(nodes)) {
    if(i < length(nodes)){
      int1 <- nodes[i]
      int2 <- nodes[i+1]
      walk <- shortest_paths(graph = input, from = int1, to = int2, mode = "all")
      walk <- unlist(walk$vpath)
      path <- c(path, walk)
    } else {
      int1 <- nodes[i]
      int2 <- nodes[1]
      walk <- shortest_paths(graph = input, from = int1, to = int2, mode = "all")
      walk <- unlist(walk$vpath)
      path <- c(path, walk)
      path <- unique(path)
      path <- names(V(input))[path]
    }
  }

  if(expand > 0){
    output_graph <- make_ego_graph(input, order = expand, nodes = nodes, mode = "all")
    keep_v <- c()
    for (j in 1:length(output_graph)) {
      gtmp <- output_graph[[j]]
      ed <- names(V(gtmp))
      keep_v <- c(keep_v, ed)
      keep_v <- unique(keep_v)
    }
    path <- unique(c(path, keep_v))
  }

  keep <- match(path, names(V(input)))
  output <- delete.vertices(input, V(input)[-keep])

  return(output)

}

