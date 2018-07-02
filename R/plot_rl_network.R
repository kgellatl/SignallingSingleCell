#' Identifies all R / L interactions
#'
#' This function will map all RL interactions
#'
#' @param input the input full_network from calc_rc_network
#' @param group_by the pData columns calc_rl_network was calculated on to split the networks
#' @param layout nicely, kk, circle
#' @param write_interactive whether or not to write an interactive visNetwork html object
#' @param interactive_groups the dropdown menu for selection nodes, either "nodes", "group_by", or "community"
#' into independent networks
#' @export
#' @details
#' This will use the calc_agg_bulk results to ID networks
#' @examples
#' ex_sc_example <- id_rl(input = ex_sc_example)

plot_rl_network <- function(input, group_by = FALSE, layout = "nicely", write_interactive = TRUE, interactive_groups = "nodes", nodesize = 3, textsize = 0.5){
  ##### Colors to match ggplot #####
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  n = length(unique(as.character(input$Summary[,"Lig_produce"])))
  dynamic_colors = gg_color_hue(n)
  ##### Arrowhead Sizing #####
  eqarrowPlot <- function(graph, layout, edge.lty=rep(1, ecount(graph)),
                          edge.arrow.size=rep(1, ecount(graph)),
                          vertex.shape="circle",
                          edge.curved=rep(0.0, length(E(graph))), ...) {
    plot(graph, edge.lty=0, edge.arrow.size=0, layout=layout,
         vertex.shape="none", edge.color=edge.col)
    for (e in seq_len(ecount(graph))) {
      graph2 <- delete.edges(graph, E(graph)[(1:ecount(graph))[-e]])
      plot(graph2, edge.lty=edge.lty[e], edge.arrow.size=edge.arrow.size[e],
           edge.curved=edge.curved[e], edge.color=edge.col[e], layout=layout, vertex.shape="none",
           vertex.label=NA, add=TRUE, ...)
    }
    plot(graph, edge.lty=0, edge.arrow.size=0, layout=layout,
         vertex.shape=vertex.shape, edge.color=edge.col, add=TRUE, ...)
    invisible(NULL)
  }

  ##### Group by  #####
  if(group_by != FALSE){
    tmpdat <- input$full_network[,c(2,4,1,3,5)]
    for (i in 1:nrow(tmpdat)) {
      int <- tmpdat[i,]
      for (j in 1:2) {
        if(j == 1){
          tmpdat[i,6] <- paste0(unlist(int[,c(j,c(3,5))]), collapse = "_")
        } else {
          tmpdat[i,7] <- paste0(unlist(int[,c(j,c(4,5))]), collapse = "_")
        }
      }
    }
    net_graph <- graph_from_data_frame(tmpdat[,c("V6", "V7")], directed = TRUE)
  }
  ##### No Group  by  #####
  if(group_by == FALSE){
    tmpdat <- input$full_network[,c(2,4,1,3)]
    for (i in 1:nrow(tmpdat)) {
      int <- tmpdat[i,]
      for (j in 1:2) {
        if(j == 1){
          tmpdat[i,5] <- paste0(unlist(int[,c(j,c(3))]), collapse = "_")
        } else {
          tmpdat[i,6] <- paste0(unlist(int[,c(j,c(4))]), collapse = "_")
        }
      }
    }
    net_graph <- graph_from_data_frame(tmpdat[,c("V5", "V6")], directed = TRUE)
  }
  ##### Color Edges  #####
  cols <- sort(unique(c(tmpdat[,3], tmpdat[,4])))
  dynamic_colors = gg_color_hue(length(cols))
  cols <- matrix(c(cols, dynamic_colors), ncol = 2)
  colors_edge <- c()
  for (i in 1:nrow(tmpdat)) {
    ind <- match(tmpdat[i,3], cols[,1])
    col <- cols[ind,2]
    colors_edge <- c(colors_edge, col)
  }
  ##### Color Vertices and get groups (nodes by default) #####
  colors_vert <- c()
  vertcol <- names(V(net_graph))
  names <- c()
  groups <- c()
  for (i in 1:length(vertcol)) {
    int <- unlist(strsplit(vertcol[i], split = "_"))[2]
    ind <- match(int, cols[,1])
    col <- cols[ind,2]
    colors_vert <- c(colors_vert, col)
    names <- c(names, unlist(strsplit(vertcol[i], split = "_"))[1])
    groups <- c(groups, int)
  }
  V(net_graph)$group <- groups
  name_backup <- V(net_graph)$name
  V(net_graph)$name <- names
  if(group_by != FALSE){
    V(net_graph)$group_by <- NA
    for (i in 1:length(V(net_graph))) {
      int <- names(V(net_graph))[i]
      sk <- strsplit(int, "-")[[1]][length(strsplit(int, "-")[[1]])]
      V(net_graph)$group_by[i] <- sk
    }
  }
  ##### Graphing parameters #####
  V(net_graph)$size <- nodesize
  V(net_graph)$label.cex <- textsize
  V(net_graph)$label.color <- "black"
  V(net_graph)$vertex.frame.color <- NA
  V(net_graph)$color <- colors_vert
  edge_widths <- input$full_network$log10_Connection
  edge_widths <- scale(edge_widths)
  edge_widths <- (edge_widths + abs(min(edge_widths)))
  E(net_graph)$width <- (2/max(edge_widths)*edge_widths)
  E(net_graph)$width <- E(net_graph)$width+1
  E(net_graph)$arrow.size <- 0.1
  E(net_graph)$color <- colors_edge

  if(layout == "nicely"){
    l <- layout_nicely(net_graph)
  }
  if(layout == "kk"){
    l <- layout_with_kk(net_graph)
  }
  if(layout == "circle"){
    l <- layout_in_circle(net_graph)
  }

  l <- norm_coords(l, ymin=0, ymax=1, xmin=0, xmax=1)
  edge.curved=
  pdf("Fullnetwork.pdf", h = 8, w = 8, useDingbats = FALSE)
  plot(net_graph, layout = l, edge.curved=curve_multiple(net_graph), vertex.frame.color = NA, cex.col= "black", rescale = TRUE)
  cell_legend <- sort(unique(tmpdat[,3]))
  legend(x=-1.5, y=0, cell_legend, pch=21,
         col="#777777", pt.bg=rev(dynamic_colors), pt.cex=2, cex=.8, bty="n", ncol=1)
  dev.off()

  pdf("Fullnetwork_communities.pdf", h = 8, w = 8, useDingbats = FALSE)
  cfg <- cluster_edge_betweenness(as.undirected(net_graph))
  plot(cfg, net_graph, layout = l, edge.curved=curve_multiple(net_graph), vertex.frame.color = NA, cex.col= "black", rescale = TRUE)
  dev.off()

  #####

}


#####

if(write_interactive == TRUE){
  V(net_graph)$name <- name_backup
  nodes <- igraph::as_data_frame(net_graph, what = "vertices")
  links <- igraph::as_data_frame(net_graph, what = "edges")
  links$arrows <- "to"

  colnames(nodes)[1] <- "id"
  nodes <- nodes[,c("id", "color")]

  nodes$label <- V(net_graph)$name
  links$width <- E(net_graph)$width

  cfg <- cluster_edge_betweenness(as.undirected(net_graph))
  nodes$community <- cfg$membership

  nodes$nodes <- V(net_graph)$group

  if(group_by != FALSE){
    nodes$condition <- V(net_graph)$skin
  }

  vit_net <- visNetwork::visNetwork(nodes, links, width="100%", height="1000px")

  vit_net <- visOptions(vit_net, highlightNearest = TRUE, selectedBy = "nodes")

  if(interactive_groups == "condition"){
    vit_net <- visOptions(vit_net, highlightNearest = TRUE, selectedBy = "condition")
  }

  if(interactive_groups == "community"){
    vit_net <- visOptions(vit_net, highlightNearest = TRUE, selectedBy = "community")
  }

  visNetwork::visSave(vit_net, file="Interactive_Network.html")
}




