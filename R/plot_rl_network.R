#' Identifies all R / L interactions
#'
#' This function will map all RL interactions
#'
#' @param input the input full_network from calc_rc_network
#' @param group_by the pData columns calc_rl_network was calculated on to split the networks
#' @param mode the network plot type
#' into independent networks
#' @export
#' @details
#' This will use the calc_agg_bulk results to ID networks
#' @examples
#' ex_sc_example <- id_rl(input = ex_sc_example)

plot_rl_network <- function(input, group_by = FALSE, mode = "Summary"){

  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }

  n = length(unique(as.character(input$Summary$V1)))
  dynamic_colors = gg_color_hue(n)

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

  ######

  if(mode == "Summary"){
    net_dat <- input$Summary
    for (i in 1:ncol(net_dat)) {
      net_dat[,i] <- as.character(net_dat[,i])
    }
    network <- c()
    if(group_by!= FALSE){
      groups <- unique(net_dat$V3)
      for (i in 1:length(groups)) {
        int <- groups[i]
        tab <- net_dat[grep(int, net_dat$V3),1:3]
        for(j in 1:nrow(tab)){
          int2 <- tab[j,]
          tab[j,1] <- paste0(c(int2$V1, int2$V3), collapse = "_")
          tab[j,2] <- paste0(c(int2$V2, int2$V3), collapse = "_")
          net_dat[grep(int, net_dat$V3),1:3] <- tab
        }
      }
      for (i in 1:nrow(net_dat)) {
        add <- as.vector(net_dat[i,1:2])
        network <- append(network, add)
      }
      net_graph <- igraph::graph(c(network))

    } else {
      net_dat_final <- net_dat
      for (i in 1:nrow(net_dat_final)) {
        tmp <- net_dat_final[i,]
        tmp$V1 <- paste0(c("Lig", tmp$V1), collapse = "_")
        tmp$V2 <- paste0(c("Rec", tmp$V2), collapse = "_")
        net_dat_final[i,] <- tmp
      }
      net_graph <- igraph::graph_from_data_frame(net_dat_final)

      ##### Layout #####
      l <- layout_in_circle(net_graph)
      # rownames(l) <- names(V(net_graph))
      # ligs <- grep("Lig", rownames(l))
      # recs <- grep("Rec", rownames(l))
      # l[ligs,1] <- 1
      # l[recs,1] <- 2
      # y_pos <- seq(1:length(unique(net_dat$V1)))
      # l[,2] <- rep(y_pos, 2)
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
      E(net_graph)$width <- rank(as.numeric(net_dat_final$freq))
      E(net_graph)$width <- (10/max(E(net_graph)$width)*E(net_graph)$width)
      ##### Remove Names #####
      V(net_graph)$name <- ""
      #####
      eqarrowPlot(net_graph, layout = l, edge.arrow.size=(2/max(E(net_graph)$width)*E(net_graph)$width),
                  edge.width=E(net_graph)$size)
      ##### Add Legend ####
      cell_legend <- unique(net_dat$V1)
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
