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

analyze_rl_network <- function(input){

  tmp_net <- input
  l <- layout_nicely(tmp_net)

  pdf("Fullnetwork.pdf", h = 8, w = 8, useDingbats = FALSE)
  plot(tmp_net, layout = l, edge.curved=curve_multiple(tmp_net), vertex.frame.color = NA, cex.col= "black", rescale = TRUE)
  dev.off()

  pdf("Fullnetwork_noname.pdf", h = 8, w = 8, useDingbats = FALSE)
  net2 <- tmp_net
  V(net2)$name <- ""
  plot(net2, layout = l, edge.curved=curve_multiple(tmp_net), vertex.frame.color = NA, cex.col= "black", rescale = TRUE)
  dev.off()

  tmp_net_cp <- tmp_net
  cfg <- cluster_edge_betweenness(as.undirected(tmp_net_cp))

  pdf("Fullnetwork_communities.pdf", h = 8, w = 8, useDingbats = FALSE)
  plot(cfg, tmp_net_cp, layout = l, edge.curved=curve_multiple(tmp_net), vertex.frame.color = NA, cex.col= "black", rescale = TRUE)
  dev.off()

  pdf("Dendrogram.pdf", height = 10, width = 50)
  dendPlot(cfg, mode="hclust")
  dev.off()

  #####
  nodes <- igraph::as_data_frame(tmp_net, what = "vertices")
  links <- igraph::as_data_frame(tmp_net, what = "edges")
  links$arrows <- "to"

  colnames(nodes)[1] <- "id"
  nodes <- nodes[,c("id", "color")]

  nodes$label <- V(tmp_net)$name
  links$value <- E(tmp_net)$width

  deg <- degree(tmp_net, mode="all")
  deg <- rank(deg)
  deg <- (20/max(deg)*deg)
  deg[which(deg < 0.5)] <- 0.5
  nodes$value <- deg

  nodes$community <- cfg$membership

  links$width <- 3

  vit_net <- visNetwork::visNetwork(nodes, links, width="100%", height="1000px")

  vit_net <- visNetwork::visOptions(vit_net, highlightNearest = TRUE, selectedBy = "community")

  visNetwork::visSave(vit_net, file="Interactive_Network_analyzed.html")

  #####

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
  results[[9]] <- vit_net

  names(results) <- c("Degree", "Node_degree", "Node_betweeness", "Node_authority",
                      "Edge_degree", "Edge_betweeness", "Edge_hub", "Communities", "Interactive")
  return(results)
}






