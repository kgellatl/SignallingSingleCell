#' Identifies all R / L interactions
#'
#' This function will map all RL interactions
#'
#' @param input the input networked (usually filtered)
#' @param input_full the input full network
#' @param group_by the pData columns calc_rl_network was calculated on to split the networks
#' @param comparitive If true will look at connections across group_by and color edges by their change
#' @param from the comparitive reference
#' @param to the condition to compare to reference
#' @param value the column of the full network to calculate foldChange on
#' @param write_interactive whether or not to write an interactive visNetwork html object
#' @param interactive_groups the dropdown menu for selection nodes, either "nodes", "group_by", or "cluster"
#' @param nodesize The size of nodes
#' @param size_by_connections if true will override node size and size by the degree of the node
#' @param textsize The size of text
#' @param h pdf height
#' @param w pdf width
#' @param prefix a character to be appeneded to the start of file names
#' into independent networks
#' @export
#' @details
#' This will use the calc_agg_bulk results to ID networks
#' @examples
#' ex_sc_example <- id_rl(input = ex_sc_example)

plot_rl_network <- function(input, value = "log10_Connection_product", group_by = FALSE, comparitive = FALSE, from = FALSE, to = FALSE, input_full = NA,
                            write_interactive = TRUE, interactive_groups = "nodes", nodesize = 3, size_by_connections = TRUE,
                            textsize = 0.5, h = 8, w = 8, prefix = ""){
  ###############################################################################################
  ##### Colors to match ggplot #####
  ###############################################################################################

  plot_rl_results <- list()
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  n = length(unique(as.character(input$Summary[,"Lig_produce"])))
  dynamic_colors = gg_color_hue(n)

  ###############################################################################################
  ##### Group by  #####
  ###############################################################################################

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
      tmpdat$expression <- input$full_network[,value]
    }
    net_graph <- graph_from_data_frame(tmpdat[,c("V6", "V7")], directed = TRUE)
  }

  ###############################################################################################
  ##### No Group  by  #####
  ###############################################################################################

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
      tmpdat$expression <- input$full_network[,value]
    }
    net_graph <- graph_from_data_frame(tmpdat[,c("V5", "V6")], directed = TRUE)
  }

  ###############################################################################################
  ##### Color Edges  #####
  ###############################################################################################

  cols <- sort(unique(c(tmpdat[,3], tmpdat[,4])))
  dynamic_colors = gg_color_hue(length(cols))
  cols <- matrix(c(cols, dynamic_colors), ncol = 2)
  colors_edge <- c()
  for (i in 1:nrow(tmpdat)) {
    ind <- match(tmpdat[i,3], cols[,1])
    col <- cols[ind,2]
    colors_edge <- c(colors_edge, col)
  }

  ###############################################################################################
  ##### Color Edges  FOR COMPARITIVE MUCH MORE COMPLICATED #####
  ###############################################################################################

  if(comparitive == TRUE){
    if(group_by == FALSE || from == FALSE || to == FALSE || value == FALSE){
      stop("Comparisons are only valid across groups, provide a group_by, from, and to")
    }
    search <- apply(input$full_network[ , 1:4 ] , 1 , paste , collapse = "_" )
    search <- unique(search)
    new_dat <- as.data.frame(matrix(unlist(strsplit(search, split = "_")), ncol = 4, byrow = T))
    new_dat$FC <- 0
    new_dat$expression <- 0
    if(!is.null(input$full_network$keep)){
      new_dat$significant <- FALSE
    }
    search_full <- apply(input_full$full_network[ , 1:4 ] , 1 , paste , collapse = "_" )
    remove <- c()
    for (i in 1:length(search)) {
      int <- search[i]
      ind <- grep(paste0("^", int, "$"), search_full)
      int2 <- input_full$full_network[ind, ]
      vals <- int2[,group_by]
      pos1 <- match(from, vals)
      pos2 <- match(to, vals)
      from_val <- int2[pos1,value]
      to_val <- int2[pos2,value]
      FC <- from_val - to_val
      exp <- max(int2[,value])
      new_dat$expression[i] <- exp
      if(all(is.na(c(from_val, to_val)))){
        remove <- c(remove, i)
      }
      if(length(which(is.na(c(from_val, to_val)))) == 1){
        pos <- which(is.na(c(from_val, to_val)))
        if(pos == 1){
          FC <- "ON"
        } else {
          FC <- "OFF"
        }
      }
      new_dat$FC[i] <- FC
      if(!is.null(input$full_network$keep)){
        tmpvals <- input$full_network$keep[ind]
        state <- unique(tmpvals)
        new_dat$significant[i] <- state
      }
    }
    if(!is.null(remove)){
      new_dat <- new_dat[-remove,]
    }
    new_dat_keep <- new_dat
    new_dat <- new_dat[,c(2,4,1,3,5,6)]
    tmpdat <- new_dat
    for (i in 1:nrow(tmpdat)) {
      int <- tmpdat[i,]
      for (j in 1:2) {
        if(j == 1){
          tmpdat[i,7] <- paste0(unlist(int[,c(j,c(3))]), collapse = "_")
        } else {
          tmpdat[i,8] <- paste0(unlist(int[,c(j,c(4))]), collapse = "_")
        }
      }
    }
    net_graph <- graph_from_data_frame(tmpdat[,c("V7", "V8")], directed = TRUE)
    new_dat_BACKUP <- new_dat
    if(length(which(new_dat_BACKUP$FC == "ON") >= 1) || length(which(new_dat_BACKUP$FC == "OFF") >= 1)){
      new_dat_BACKUP <- new_dat_BACKUP[-which(new_dat_BACKUP$FC == "ON"),]
      new_dat_BACKUP <- new_dat_BACKUP[-which(new_dat_BACKUP$FC == "OFF"),]
    }
    max_val <- max(as.numeric(new_dat_BACKUP$FC))
    min_val <- min(as.numeric(new_dat_BACKUP$FC))
    new_dat$FC[which(new_dat$FC == "ON")] <- min_val
    new_dat$FC[which(new_dat$FC == "OFF")] <- max_val
    new_dat$FC <- as.numeric(new_dat$FC)
    new_dat$color <- NA

    ##### Need to take this and split into positive and negative sections
    col2s <- (viridis::cividis(15))
    # plot(seq(1:15), col = col2s, pch = 20)

    colfunc_low <- colorRampPalette(col2s[1:5])
    lowind <- which(new_dat$FC > 0)
    cols2_down <- colfunc_low(length(lowind))
    new_dat$color[lowind][order(-new_dat$FC[lowind])] <- cols2_down

    colfunc_high <- colorRampPalette(col2s[11:15])
    highind <- which(new_dat$FC < 0)
    cols2_up <- colfunc_high(length(highind))
    new_dat$color[highind][order(-new_dat$FC[highind])] <- cols2_up

    alternative_colors <- new_dat$color
    #####

    if(!is.null(input$full_network$keep)){
      grayout <- which(new_dat_keep$significant == FALSE)
      alternative_colors[grayout] <- col2s[8]
    }

    colors_edge <- c()
    for (i in 1:nrow(tmpdat)) {
      ind <- match(tmpdat[i,3], cols[,1])
      col <- cols[ind,2]
      colors_edge <- c(colors_edge, col)
    }

  }

  ###############################################################################################
  ##### Color Vertices and get groups (nodes by default) #####
  ###############################################################################################

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
  # V(net_graph)$name <- names
  if(group_by != FALSE){
    V(net_graph)$group_by <- NA
    for (i in 1:length(V(net_graph))) {
      int <- names(V(net_graph))[i]
      sk <- strsplit(int, "-")[[1]][length(strsplit(int, "-")[[1]])]
      V(net_graph)$group_by[i] <- sk
    }
  }

  ###############################################################################################
  ##### Graphing parameters #####
  ###############################################################################################

  V(net_graph)$size <- nodesize
  V(net_graph)$label.cex <- textsize
  V(net_graph)$label.color <- "black"
  V(net_graph)$vertex.frame.color <- "white"
  V(net_graph)$color_celltype <- colors_vert

  E(net_graph)$arrow.size <- 0.1
  E(net_graph)$color_celltype <- colors_edge

  l <- layout_components(net_graph, layout = layout_with_kk)

  if(comparitive!=FALSE){
    E(net_graph)$color_compare <- alternative_colors
  }

  if(comparitive!=FALSE){
    newsize <- new_dat$expression
    newsize <- rank(-newsize)
    newsize <- (3/max(newsize)*newsize)
    E(net_graph)$width <- newsize
  } else {
    newsize <- tmpdat$expression
    newsize <- rank(-newsize)
    newsize <- (3/max(newsize)*newsize)
    E(net_graph)$width <- newsize
  }

  if(size_by_connections == TRUE){
    deg <- degree(net_graph, mode="all")
    deg <- rank(deg)
    deg <- (3/max(deg)*deg)
    deg[which(deg < 0.5)] <- 0.5
    V(net_graph)$size <- deg
  }

  ###############################################################################################
  ##### Write out results #####
  ###############################################################################################

  l <- norm_coords(l, ymin=0, ymax=1, xmin=0, xmax=1)
  pdf(paste0(prefix, "Fullnetwork.pdf"), h = h, w = w, useDingbats = FALSE)
  plot(net_graph, layout = l, rescale = TRUE,
       vertex.color = V(net_graph)$color_celltype,
       edge.color = E(net_graph)$color_celltype)
  cell_legend <- sort(unique(tmpdat[,3]))
  legend(x=-1.5, y=0, cell_legend, pch=21,
         col="#777777", pt.bg=rev(dynamic_colors), pt.cex=2, cex=.8, bty="n", ncol=1)
  dev.off()

  pdf(paste0(prefix, "Fullnetwork_noname.pdf"), h = h, w = w, useDingbats = FALSE)
  V(net_graph)$name_blank <- ""
  plot(net_graph, layout = l, rescale = TRUE,
       vertex.color = V(net_graph)$color_celltype,
       edge.color = E(net_graph)$color_celltype,
       vertex.label = V(net_graph)$name_blank)
  dev.off()

  if(comparitive == TRUE){
    pdf(paste0(prefix, "Fullnetwork_noname_compare.pdf"), h = h, w = w, useDingbats = FALSE)
    plot(net_graph, layout = l, rescale = TRUE,
         vertex.color = V(net_graph)$color_celltype,
         edge.color = E(net_graph)$color_compare,
         vertex.label = V(net_graph)$name_blank)
    dev.off()
  }

  subgraphs <- igraph::clusters(net_graph)
  cols_clust <- gg_color_hue(length(unique(subgraphs$membership)))
  clusts <- as.vector(subgraphs$membership)
  for (i in 1:length(cols_clust)) {
    cl <- cols_clust[i]
    clusts[which(clusts == i)] <- cl
  }

  pdf(paste0(prefix, "Fullnetwork_clusters.pdf"), h = h, w = w, useDingbats = FALSE)
  V(net_graph)$name_membership <- subgraphs$membership
  V(net_graph)$color_membership <- clusts
  plot(net_graph, layout = l, rescale = TRUE,
       vertex.color = V(net_graph)$color_membership,
       edge.color = E(net_graph)$color_celltype,
       vertex.label = V(net_graph)$name_membership)
  dev.off()

  plot_rl_results[[1]] <- net_graph
  plot_rl_results[[2]] <- l
  plot_rl_results[[3]] <- subgraphs
  plot_rl_results[[4]] <- decompose.graph(net_graph)
  names(plot_rl_results) <- c("igraph_Network", "layout", "clusters", "clusters_subgraphs")

  if(comparitive!= FALSE){
    plot_rl_results[[5]] <- new_dat_BACKUP
    names(plot_rl_results) <- c("igraph_Network", "layout", "clusters", "clusters_subgraphs", "comparitive_table")
  }

  ###############################################################################################
  ##### Interactive #####
  ###############################################################################################

  if(write_interactive == TRUE){
    V(net_graph)$name <- name_backup
    nodes <- igraph::as_data_frame(net_graph, what = "vertices")
    links <- igraph::as_data_frame(net_graph, what = "edges")
    links$arrows <- "to"

    colnames(nodes)[1] <- "id"
    nodes <- nodes[,c("id", "color_celltype")]

    nodes$label <- V(net_graph)$name
    links$value <- E(net_graph)$width

    if(size_by_connections == TRUE){
      deg <- degree(net_graph, mode="all")
      deg <- rank(deg)
      deg <- (20/max(deg)*deg)
      deg[which(deg < 0.5)] <- 0.5
      nodes$value <- deg
    }

    nodes$cluster <- as.vector(plot_rl_results$clusters$membership)

    nodes$nodes <- V(net_graph)$group
    links$width <- 3

    nodes$color <-  V(net_graph)$color_celltype
    links$color <- E(net_graph)$color_celltype

        if(group_by != FALSE){
      nodes$condition <- V(net_graph)$skin
    }

    nodes <- nodes[order(nodes$id),]

    vit_net <- visNetwork::visNetwork(nodes, links, width="100%", height="1000px")

    vit_net <- visNetwork::visOptions(vit_net, highlightNearest = TRUE, selectedBy = "nodes", nodesIdSelection = TRUE)
    if(interactive_groups == "condition"){
      vit_net <- visNetwork::visOptions(vit_net, highlightNearest = TRUE, selectedBy = "condition")
    }
    if(interactive_groups == "cluster"){
      vit_net <- visNetwork::visOptions(vit_net, highlightNearest = TRUE, selectedBy = "cluster")
    }
    visNetwork::visSave(vit_net, file=paste0(prefix, "Fullnetwork_interactive.html"))

    if(comparitive!= FALSE){
      plot_rl_results[[6]] <- vit_net
      names(plot_rl_results) <- c("igraph_Network", "layout", "clusters", "clusters_subgraphs", "comparitive_table", "interactive")
    } else {
      plot_rl_results[[5]] <- vit_net
      names(plot_rl_results) <- c("igraph_Network", "layout", "clusters", "clusters_subgraphs", "interactive")
    }
  }
  return(plot_rl_results)
}
