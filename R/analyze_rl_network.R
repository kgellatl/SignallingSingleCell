#' Identifies all R / L interactions
#'
#' This function will map all RL interactions
#'
#' @param input the input network from plot_rl_network
#' @param h the height of images
#' @param w the height of images
#' @param prefix a character to be appended to start of file names
#' @export
#' @details
#' This will use the calc_agg_bulk results to ID networks
#' @examples
#' ex_sc_example <- id_rl(input = ex_sc_example)

analyze_rl_network <- function(input, h = 8, w = 8, prefix = ""){

  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }

  tmp_net <- input
  l <- layout_nicely(tmp_net)

  tmp_net_cp <- tmp_net
  cfg <- cluster_edge_betweenness(as.undirected(tmp_net_cp))

  cs2 <- crossing(cfg, tmp_net_cp)


  cols_clust <- gg_color_hue(length(unique(cfg$membership)))
  cols_clust2 <- adjustcolor(cols_clust, alpha.f = 1)
  clusts <- as.vector(cfg$membership)

  for (i in 1:length(cols_clust)) {
    cl <- cols_clust[i]
    clusts[which(clusts == i)] <- cl
  }

  V(tmp_net_cp)$color <- clusts

  #####

  deg <- degree(tmp_net, mode="all")
  deg <- sort(deg, decreasing = T)

  deg_rec <- degree(tmp_net, mode="in")
  deg_rec <- sort(deg_rec, decreasing = T)

  deg_lig <- degree(tmp_net, mode="out")
  deg_lig <- sort(deg_lig, decreasing = T)

  vert_rank <- (betweenness(as.undirected(tmp_net))) # VERTICES
  vert_rank_srt <- sort(vert_rank, decreasing = T)

  edg_rank <- (edge_betweenness(as.undirected(tmp_net))) # EDGES
  tmp <- igraph::as_edgelist(tmp_net)
  tmp2 <- apply( tmp[ , 1:2 ] , 1 , paste , collapse = "-" )
  names(edg_rank) <- tmp2
  edg_rank <- sort(edg_rank, decreasing = T)

  hs <- hub_score(tmp_net, weights=NA)$vector #Outgoing (ligands)
  hs <- sort(hs, decreasing = T)

  as <- authority_score(tmp_net, weights=NA)$vector # Incoming (receptors)
  as <- sort(as, decreasing = T)

  cluster_edge_btw <- vector(mode = "list", length = length(unique(cfg$membership)))
  sizes <- rep(0, length(V(tmp_net)))

  for (i in 1:length(unique(cfg$membership))) {
    int <- which(cfg$membership == i)
    vals <- vert_rank[int]
    cluster_edge_btw[[i]] <- sort(vals, decreasing = T)
    vals <- rank(vals)
    size <- (3/max(vals)*vals)
    sizes[int] <- size
  }

  cluster_centers <- c()
  for(i in 1:length(cluster_edge_btw)){
    name_center <- cluster_edge_btw[[i]][1]
    cluster_centers <- c(cluster_centers,name_center)
  }

  V(tmp_net_cp)$size <- sizes

  coGrph <- delete_edges(tmp_net, E(tmp_net)[crossing(cfg, tmp_net)])
  comm_ind <- igraph::decompose.graph(coGrph)

  #####

  pdf(paste0(prefix, "Analyzed_Network.pdf"), h = h, w = w, useDingbats = FALSE)
  plot(tmp_net, layout = l, edge.curved=curve_multiple(tmp_net), vertex.frame.color = NA, cex.col= "black", rescale = TRUE)
  dev.off()

  pdf(paste0(prefix, "Analyzed_Network_noname.pdf"), h = h, w = w, useDingbats = FALSE)
  net2 <- tmp_net
  V(net2)$name <- ""
  plot(net2, layout = l, edge.curved=curve_multiple(tmp_net), cex.col= "black", vertex.frame.color = 'white', rescale = TRUE)
  dev.off()

  pdf(paste0(prefix, "Analyzed_Network_communities.pdf"), h = h, w = w, useDingbats = FALSE)
  plot(tmp_net_cp, mark.groups = cfg,  layout = l, edge.curved=curve_multiple(tmp_net), vertex.frame.color = 'white', cex.col= "black", rescale = TRUE)
  dev.off()


  pdf(paste0(prefix, "Analyzed_Network_communities_nonames.pdf"), h = h, w = w, useDingbats = FALSE)
  plot(tmp_net_cp, mark.groups = cfg, vertex.label = "",  layout = l, edge.curved=curve_multiple(tmp_net), vertex.frame.color = 'white', vertex.frame.color = 'white', cex.col= "black", rescale = TRUE)
  dev.off()

  pdf(paste0(prefix, "Analyzed_Network_communities_crossing.pdf"), h = h, w = w, useDingbats = FALSE)
  tmp_net_cp2 <- tmp_net_cp
  E(tmp_net_cp2)$color <- "gray"
  E(tmp_net_cp2)$color[crossing(cfg,tmp_net_cp2)] <- "red"
  plot(tmp_net_cp2, mark.groups = cfg, vertex.label = "",  layout = l, edge.curved=curve_multiple(tmp_net), vertex.frame.color = 'white', cex.col= "black", rescale = TRUE)
  dev.off()

  pdf(paste0(prefix, "Analyzed_Network_communities_crossing_named.pdf"), h = h, w = w, useDingbats = FALSE)
  keep_name <- unique(ends(tmp_net_cp2, names(cs2)[which(cs2 == TRUE)], names =  T))
  V(tmp_net_cp2)$name[-match(keep_name, V(tmp_net_cp2)$name)] <- ""
  plot(tmp_net_cp2, layout = l, edge.curved=curve_multiple(tmp_net), vertex.frame.color = 'white', cex.col= "black", rescale = TRUE)
  dev.off()

  pdf(paste0(prefix, "Analyzed_Network_communities_centers_named.pdf"), h = h, w = w, useDingbats = FALSE)
  tmp_net3 <-  tmp_net_cp
  keep_center <- match(names(cluster_centers), names(V(tmp_net3)))
  V(tmp_net3)$name[-keep_center] <- ""
  plot(tmp_net3,  layout = l, edge.curved=curve_multiple(tmp_net), vertex.frame.color = 'white', cex.col= "black", rescale = TRUE)
  dev.off()


  pdf(paste0(prefix, "Analyzed_Network_communities_numbered.pdf"), h = h, w = w, useDingbats = FALSE)
  V(tmp_net_cp)$name <- cfg$membership
  plot(tmp_net_cp, layout = l, edge.curved=curve_multiple(tmp_net), vertex.frame.color = 'white', cex.col= "black", rescale = TRUE)
  dev.off()

  # pdf(paste0(prefix, "Dendrogram.pdf"), height = 10, width = 50)
  # dendPlot(cfg, mode="hclust")
  # dev.off()

  #####

  nodes <- igraph::as_data_frame(tmp_net, what = "vertices")
  links <- igraph::as_data_frame(tmp_net, what = "edges")
  links$arrows <- "to"

  colnames(nodes)[1] <- "id"
  nodes <- nodes[,c("id", "color")]

  nodes$label <- V(tmp_net)$name
  links$value <- E(tmp_net)$width

  size2 <- (20/max(sizes)*sizes)
  nodes$value <- size2

  nodes$community <- cfg$membership

  links$width <- 3

  nodes <- nodes[order(nodes$id),]

  vit_net <- visNetwork::visNetwork(nodes, links, width="100%", height="1000px")

  vit_net <- visNetwork::visOptions(vit_net, highlightNearest = TRUE, selectedBy = "community", nodesIdSelection = TRUE)

  visNetwork::visSave(vit_net, file=paste0(prefix, "Interactive_Network_analyzed.html"))

  #####

  results <- vector(mode = "list", length = 10)
  results[[1]] <- deg
  results[[2]] <- deg_rec
  results[[3]] <- vert_rank_srt
  results[[4]] <- as
  results[[5]] <- deg_lig
  results[[6]] <- edg_rank
  results[[7]] <- hs
  results[[8]] <- cfg
  results[[9]] <- cs2
  results[[10]] <- cluster_edge_btw
  results[[11]] <- as.data.frame(matrix(c(names(cluster_centers), cluster_centers), ncol = 2))
  results[[12]] <- comm_ind
  results[[13]] <- vit_net
  results[[14]] <- input

  names(results) <- c("Degree", "Node_degree", "Node_betweeness", "Node_authority",
                      "Edge_degree", "Edge_betweeness", "Edge_hub", "Communities", "Crossing",
                      "Communities_betweeness", "Community_centers", "Communities_individual", "Interactive", "input")
  return(results)
}






