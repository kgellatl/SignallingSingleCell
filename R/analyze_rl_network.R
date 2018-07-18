#' Identifies all R / L interactions
#'
#' This function will map all RL interactions
#'
#' @param input the input network from plot_rl_network
#' @export
#' @details
#' This will use the calc_agg_bulk results to ID networks
#' @examples
#' ex_sc_example <- id_rl(input = ex_sc_example)

analyze_rl_network <- function(input, h = 5, w = 20){

  tmp_net <- input
  l <- layout_nicely(tmp_net)

  cfg <- cluster_edge_betweenness(as.undirected(tmp_net))

  pdf("Fullnetwork_communities.pdf", h = 8, w = 8, useDingbats = FALSE)
  plot(cfg, tmp_net, layout = l, edge.curved=curve_multiple(net_graph), vertex.frame.color = NA, cex.col= "black", rescale = TRUE)
  dev.off()

  pdf("Dendrogram.pdf", height = 10, width = 50)
  dendPlot(cfg, mode="hclust")
  dev.off()

  deg <- degree(tmp_net, mode="all")
  deg <- sort(deg, decreasing = T)

  deg_rec <- degree(tmp_net, mode="in")
  deg_rec <- sort(deg_rec, decreasing = T)

  deg_lig <- degree(tmp_net, mode="out")
  deg_lig <- sort(deg_lig, decreasing = T)

  vert_rank <- (betweenness(as.undirected(tmp_net))) # VERTICES
  vert_rank <- sort(vert_rank, decreasing = T)

  edg_rank <- (edge_betweenness(as.undirected(tmp_net))) # EDGES
  tmp <- igraph::as_edgelist(tmp_net)
  tmp2 <- apply( tmp[ , 1:2 ] , 1 , paste , collapse = "-" )
  names(edg_rank) <- tmp2
  edg_rank <- sort(edg_rank, decreasing = T)

  hs <- hub_score(tmp_net, weights=NA)$vector #Outgoing (ligands)
  hs <- sort(hs, decreasing = T)

  as <- authority_score(tmp_net, weights=NA)$vector # Incoming (receptors)
  as <- sort(as, decreasing = T)

  results <- vector(mode = "list", length = 7)
  results[[1]] <- deg
  results[[2]] <- deg_rec
  results[[3]] <- vert_rank
  results[[4]] <- as
  results[[5]] <- deg_lig
  results[[6]] <- edg_rank
  results[[7]] <- hs
  results[[8]] <- cfg

  names(results) <- c("Degree", "Node_degree", "Node_betweeness", "Node_authority",
                      "Edge_degree", "Edge_betweeness", "Edge_hub", "Communities")
 return(results)
}






