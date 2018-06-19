#' Identifies all R / L interactions
#'
#' This function will map all RL interactions
#'
#' @param input the input full_network from calc_rc_network
#' @param group_by the pData columns calc_rl_network was calculated on to split the networks
#' @param mode the network plot type
#' @param normalization the network normalization types. See column names of calc_rl_network$Summary
#' into independent networks
#' @export
#' @details
#' This will use the calc_agg_bulk results to ID networks
#' @examples
#' ex_sc_example <- id_rl(input = ex_sc_example)

plot_rl_network <- function(input, group_by = FALSE, mode = "Summary", normalization = "num_connections"){
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





      if(normalization == "num_connections"){
        E(net_graph)$width <- rank(as.numeric(net_dat_final$num_connections))
        E(net_graph)$width <- (10/max(E(net_graph)$width)*E(net_graph)$width)
      }
      if(normalization == "fraction_connections"){
        E(net_graph)$width <- rank(as.numeric(net_dat_final$fraction_connections))
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





      if(normalization == "num_connections"){
        E(net_graph)$width <- rank(as.numeric(net_dat_final$num_connections))
        E(net_graph)$width <- (10/max(E(net_graph)$width)*E(net_graph)$width)
      }
      if(normalization == "fraction_connections"){
        E(net_graph)$width <- rank(as.numeric(net_dat_final$fraction_connections))
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
  }
}
