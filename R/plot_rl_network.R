#' Identifies all R / L interactions
#'
#' This function will map all RL interactions
#'
#' @param input the input full_network from calc_rc_network
#' @param group_by the pData columns calc_rl_network was calculated on to split the networks
#' @param mode the network plot type
#' @param edge_weight the network edge_weight types. See column names of calc_rl_network$Summary
#' into independent networks
#' @export
#' @details
#' This will use the calc_agg_bulk results to ID networks
#' @examples
#' ex_sc_example <- id_rl(input = ex_sc_example)

plot_rl_network <- function(input, group_by = FALSE, mode = "Summary", edge_weight = "num_connections", rank = TRUE){
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
  ##### Summary plot type is edges weighted by the number of connections #####
  if(mode == "Summary"){
    net_dat <- input$Summary
    for (i in 1:ncol(net_dat)) {
      net_dat[,i] <- as.character(net_dat[,i])
    }
    if(group_by!= FALSE){
      ##### GROUPED DATA #####
      net_dat_final <- net_dat
      for (i in 1:nrow(net_dat_final)) {
        tmp <- net_dat_final[i,]
        tmp[,"Lig_produce"] <- paste0(c(tmp[,group_by], tmp[,"Lig_produce"]), collapse = "_")
        tmp[,"Rec_receive"] <- paste0(c(tmp[,group_by], tmp[,"Rec_receive"]), collapse = "_")
        net_dat_final[i,] <- tmp
      }
      for (i in 1:nrow(net_dat_final)) {
        tmp <- net_dat_final[i,]
        tmp[,"Lig_produce"] <- paste0(c("Lig", tmp[,"Lig_produce"]), collapse = "_")
        tmp[,"Rec_receive"] <- paste0(c("Rec", tmp[,"Rec_receive"]), collapse = "_")
        net_dat_final[i,] <- tmp
      }
      order <- c()
      groups  <- unique(net_dat[,group_by])
      for (i in 1:length(groups)) {
        int <- groups[i]
        tab <- grep(int, net_dat[,group_by])
        order <- c(order , tab)
      }
      net_dat_final <- net_dat_final[order,]
      net_graph <- graph_from_data_frame(net_dat_final)
      ##### Layout #####
      l <- layout_in_circle(net_graph)
      rownames(l) <- names(V(net_graph))
      y_pos <- seq(1:length(unique(net_dat[,"Lig_produce"])))
      l[,2] <- rep(y_pos)
      xpos <- seq(1:(length(groups)*2))
      for (i in 1:length(xpos)) {
        if(i == 1){
          # do nothing
        } else {
          if(i %% 2 == 0){
            xpos[i] <- xpos[i-1]+5
          } else {
            xpos[i] <- xpos[i-1]+1
          }
        }
      }
      chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
      grouped_pos <- chunk2(xpos, length(groups))
      for (i in 1:length(groups)) {
        int <- groups[i]
        ind <- grep(int, rownames(l))
        posits <- chunk2(ind, 2)
        for (j in 1:length(posits)) {
          l[posits[[j]],1] <- grouped_pos[[i]][j]
        }
      }
      ##### Color Nodes  #####
      V(net_graph)$color <- rep(dynamic_colors, length(unique(net_dat_final[,group_by]))*2)
      ##### Color Edges  #####
      edge.col <- rep(rep(dynamic_colors, each = length(unique(net_dat[,"Lig_produce"]))), length(groups))
      ##### Size Edge and Arrows #####
      if(rank == TRUE){
        E(net_graph)$width <- rank(as.numeric(net_dat_final[,edge_weight]))
        E(net_graph)$width <- (10/max(E(net_graph)$width)*E(net_graph)$width)
      }
      if(rank == FALSE){
        E(net_graph)$width <- as.numeric(net_dat_final[,edge_weight])
        E(net_graph)$width <- (10/max(E(net_graph)$width)*E(net_graph)$width)
      }
      ##### Size Vertices #####
      verts <- names(V(net_graph))
      vert_weights <- c()
      for (i in 1:length(verts)) {
        ind <- grep(verts[i], net_dat_final[,"Lig_produce"])
        sum <- net_dat_final$num_connections[ind]
        sum <- sum(as.numeric(sum))
        vert_weights <- c(vert_weights, sum)
      }
      vert_weights <- vert_weights[which(vert_weights > 0)]
      vert_weights <- rank(as.numeric(vert_weights))
      V(net_graph)$size <- (10/max(vert_weights)*vert_weights)
      ##### Remove Names #####
      V(net_graph)$name <- ""
      ##### Plot #####
      eqarrowPlot(net_graph, layout = l, edge.arrow.size=(2/max(E(net_graph)$width)*E(net_graph)$width),
                  edge.width=E(net_graph)$size)
      ##### Add Legend ####
      cell_legend <- unique(net_dat[,"Lig_produce"])
      cell_legend <-rev(cell_legend)
      legend(x=-1.5, y=0, cell_legend, pch=21,
             col="#777777", pt.bg=rev(dynamic_colors), pt.cex=2, cex=.8, bty="n", ncol=1)
    } else {
      ##### NON GROUPED DATA #####
      net_dat_final <- net_dat
      for (i in 1:nrow(net_dat_final)) {
        tmp <- net_dat_final[i,]
        tmp[,"Lig_produce"] <- paste0(c("Lig", tmp[,"Lig_produce"]), collapse = "_")
        tmp[,"Rec_receive"] <- paste0(c("Rec", tmp[,"Rec_receive"]), collapse = "_")
        net_dat_final[i,] <- tmp
      }
      net_graph <- graph_from_data_frame(net_dat_final)
      ##### Layout #####
      l <- layout_in_circle(net_graph)
      rownames(l) <- names(V(net_graph))
      ligs <- grep("Lig", rownames(l))
      recs <- grep("Rec", rownames(l))
      l[ligs,1] <- 1
      l[recs,1] <- 2
      y_pos <- seq(1:length(unique(net_dat[,"Lig_produce"])))
      l[,2] <- rep(y_pos, 2)
      ##### Color Nodes  #####
      V(net_graph)$color <- rep(dynamic_colors, 2)
      ##### Color Edges  #####
      ecol <- c()
      for (i in 1:length(dynamic_colors)) {
        col <- rep(dynamic_colors[i], n )
        ecol <- c(ecol, col)
      }
      edge.col <- ecol
      ##### Size Edge and Arrows #####
      if(rank == TRUE){
        E(net_graph)$width <- rank(as.numeric(net_dat_final[,edge_weight]))
        E(net_graph)$width <- (10/max(E(net_graph)$width)*E(net_graph)$width)
      }
      if(rank == FALSE){
        E(net_graph)$width <- as.numeric(net_dat_final[,edge_weight])
        E(net_graph)$width <- (10/max(E(net_graph)$width)*E(net_graph)$width)
      }
      ##### Size Vertices #####
      verts <- names(V(net_graph))
      vert_weights <- c()
      for (i in 1:length(verts)) {
        ind <- grep(verts[i], net_dat_final[,"Lig_produce"])
        sum <- net_dat_final$num_connections[ind]
        sum <- sum(as.numeric(sum))
        vert_weights <- c(vert_weights, sum)
      }
      vert_weights <- vert_weights[which(vert_weights > 0)]
      vert_weights <- rank(as.numeric(vert_weights))
      V(net_graph)$size <- (10/max(vert_weights)*vert_weights)
      ##### Remove Names #####
      V(net_graph)$name <- ""
      ##### Plot #####
      eqarrowPlot(net_graph, layout = l, edge.arrow.size=(2/max(E(net_graph)$width)*E(net_graph)$width),
                  edge.width=E(net_graph)$size)
      ##### Add Legend ####
      cell_legend <- unique(net_dat[,"Lig_produce"])
      cell_legend <-rev(cell_legend)
      legend(x=-1.5, y=0, cell_legend, pch=21,
             col="#777777", pt.bg=rev(dynamic_colors), pt.cex=2, cex=.8, bty="n", ncol=1)
    }
  }

  #####

  if(mode == "RLPair"){
  }

  #####

  #####

  if(mode == "Hairball"){

    ##### Group by  #####

    # tmpdat <- val[,c(2,4,1,3,5)]
    #
    # for (i in 1:nrow(tmpdat)) {
    #   int <- tmpdat[i,]
    #   for (j in 1:2) {
    #     if(j == 1){
    #       tmpdat[i,6] <- paste0(unlist(int[,c(j,c(3,5))]), collapse = "_")
    #     } else {
    #       tmpdat[i,7] <- paste0(unlist(int[,c(j,c(4,5))]), collapse = "_")
    #     }
    #   }
    # }
    # net_graph <- graph_from_data_frame(tmpdat[,c("V6", "V7")], directed = TRUE)
    #
    # ##### No Group  by  #####
    #
    # tmpdat <- val[,c(2,4,1,3)]
    #
    # for (i in 1:nrow(tmpdat)) {
    #   int <- tmpdat[i,]
    #   for (j in 1:2) {
    #     if(j == 1){
    #       tmpdat[i,5] <- paste0(unlist(int[,c(j,c(3))]), collapse = "_")
    #     } else {
    #       tmpdat[i,6] <- paste0(unlist(int[,c(j,c(4))]), collapse = "_")
    #     }
    #   }
    # }
    # net_graph <- graph_from_data_frame(tmpdat[,c("V5", "V6")], directed = TRUE)
    #
    # head(tmpdat)
    # length(E(net_graph))
    # length(V(net_graph))
    #
    # cols <- sort(unique(Vitiligo_Network$full_network[,1]))
    # dynamic_colors = gg_color_hue(length(cols))
    # cols <- matrix(c(cols, dynamic_colors), ncol = 2)
    # colors_edge <- c()
    # for (i in 1:nrow(tmpdat)) {
    #   ind <- match(tmpdat[i,3], cols[,1])
    #   col <- cols[ind,2]
    #   colors_edge <- c(colors_edge, col)
    # }
    #
    # colors_vert <- c()
    # vertcol <- names(V(net_graph))
    # names <- c()
    # groups <- c()
    # for (i in 1:length(vertcol)) {
    #   int <- unlist(strsplit(vertcol[i], split = "_"))[2]
    #   ind <- match(int, cols[,1])
    #   col <- cols[ind,2]
    #   colors_vert <- c(colors_vert, col)
    #   names <- c(names, unlist(strsplit(vertcol[i], split = "_"))[1])
    #   groups <- c(groups, int)
    # }
    # V(net_graph)$group <- groups
    #
    #
    # V(net_graph)$size <- 1
    # V(net_graph)$label.cex <- 0.5
    # V(net_graph)$label.color <- "black"
    # V(net_graph)$vertex.frame.color <- NA
    # V(net_graph)$color <- colors_vert
    # # V(net_graph)$name <- names
    #
    # V(net_graph)$skin <- NA
    # # for (i in 1:length(V(net_graph))) {
    # #   int <- names(V(net_graph))[i]
    # #   sk <- strsplit(int, "-")[[1]][length(strsplit(int, "-")[[1]])]
    # #   V(net_graph)$skin[i] <- sk
    # # }
    #
    # E(net_graph)$width <- 0.1
    # E(net_graph)$arrow.size <- 1
    # E(net_graph)$color <- colors_edge
    # E(net_graph)$color <- colors_edge
    #
    # E(net_graph)
    # V(net_graph)
    #
    # l <- layout_nicely(net_graph)
    # # l <- layout_with_kk(net_graph)
    #
    # # l <- layout_on_sphere(net_graph)
    # # l <- layout_in_circle(net_graph)
    #
    # l <- norm_coords(l, ymin=0, ymax=1, xmin=0, xmax=1)
    #
    #
    #
    # pdf("Fullnetwork_4.pdf", h = 20, w = 20, useDingbats = FALSE)
    # plot(net_graph, layout = l, edge.curved=0.2, vertex.frame.color = NA, cex.col= "black", rescale = TRUE)
    # cell_legend <- sort(unique(tmpdat[,3]))
    # legend(x=-1.5, y=0, cell_legend, pch=21,
    #        col="#777777", pt.bg=rev(dynamic_colors), pt.cex=2, cex=.8, bty="n", ncol=1)
    # dev.off()
    #
    # pdf("Fullnetwork_communities_edge.pdf", h = 6, w = 6, useDingbats = FALSE)
    # ceb <- cluster_edge_betweenness(as.undirected(net_graph))
    # plot(ceb, net_graph, layout = l, edge.curved=0.2, vertex.frame.color = NA, cex.col= "black", rescale = TRUE)
    # dev.off()
    #
    # pdf("Fullnetwork_communities_propogation.pdf", h = 20, w = 20, useDingbats = FALSE)
    # clp <- cluster_label_prop(as.undirected(net_graph))
    # plot(clp, net_graph, layout = l, edge.curved=0.2, vertex.frame.color = NA, cex.col= "black", rescale = TRUE)
    # dev.off()
    #
    # pdf("Fullnetwork_communities_greedy.pdf", h = 6, w = 6, useDingbats = FALSE)
    # cfg <- cluster_fast_greedy(as.undirected(net_graph))
    # plot(cfg, net_graph, layout = l, edge.curved=0.2, vertex.frame.color = NA, cex.col= "black", rescale = TRUE)
    # dev.off()


  }

  #####

}


#####

#
# gjs <- threejs::graphjs(g = net_graph, vertex.label = V(net_graph)$name, vertex.shape =  V(net_graph)$name, layout = layout_with_kk(net_graph, dim = 3),
#                         width = 1000, height = 1000, edge.width = 3, vertex.size = 0.1, edge.alpha = 1)
# htmlwidgets::saveWidget(gjs, file="Interactive_network.html")
#
#
# #####
#
# nodes <- igraph::as_data_frame(net_graph, what = "vertices")
# links <- igraph::as_data_frame(net_graph, what = "edges")
# links$arrows <- "to"
#
# colnames(nodes)[1] <- "id"
# nodes <- nodes[,c("id", "color")]
#
# nodes$label <- V(net_graph)$name
#
# cfg <- cluster_fast_greedy(as.undirected(net_graph))
# nodes$community <- cfg$membership
#
# vit_net <- visNetwork::visNetwork(nodes, links, width="500px", height="500px")
#
# vit_net <- visOptions(vit_net, highlightNearest = TRUE, selectedBy = "community")
#
# # nodes$group <- V(net_graph)$group
# # vit_net <- visOptions(vit_net, highlightNearest = TRUE, selectedBy = "group")
#
# # nodes$skin <- V(net_graph)$skin
# # vit_net <- visOptions(vit_net, highlightNearest = TRUE, selectedBy = "skin")
#
# # vit_net
#
# visNetwork::visSave(vit_net, file="Net_log6_z1.5_allskin.html")
#
# #####
