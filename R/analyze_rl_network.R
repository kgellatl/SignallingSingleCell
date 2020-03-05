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

analyze_rl_network <- function(input, h = 8, w = 8, prefix = "", mult = 1, layout = FALSE, subset = FALSE, subset_on = "connected", merge_singles = T, cluster_type  = "louvain"){

  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }

  ###############################################################################################
  ##### Subsetting and determining layout #####
  ###############################################################################################

  if(subset != FALSE){
    if(subset_on == "connected"){
      tmp_net <- input$clusters_subgraphs[[subset]]
      ind2 <- match(names(V(tmp_net)), names(V(input$igraph_Network)))
      l  <- input$layout[ind2,]
    }
    if(subset_on == "clusters"){
      tmp_net <- input$Clusters_individual[[subset]]
      ind2 <- match(names(V(tmp_net)), names(V(input$igraph_Network)))
      l  <- input$layout[ind2,]
    }
  } else {
    tmp_net <- input$igraph_Network
    l <- input$layout
  }

  if(layout){
    l <- layout_with_kk(tmp_net)
  }

    ###############################################################################################
    ##### Subsetting and determining layout #####
    ###############################################################################################

  if(cluster_type == "between"){
    cfg <- cluster_edge_betweenness(as.undirected(tmp_net), weights = NULL)

  }

  if(cluster_type == "louvain"){
    cfg <- cluster_louvain(as.undirected(tmp_net), weights = NULL)

  }
    cs2 <- crossing(cfg, tmp_net)
    # coGrph <- delete_edges(tmp_net, E(tmp_net)[crossing(cfg, tmp_net)])
    # comm_ind <- igraph::decompose.graph(coGrph)

    found_clust <- sort(unique(cfg$membership))
    comm_ind <- vector(mode = "list", length = length(found_clust))
    for (i in 1:length(found_clust)) {
      nodes_int_clust <- names(V(tmp_net))[which(cfg$membership == found_clust[i])]
      comm_ind[[i]] <- induced_subgraph(tmp_net, nodes_int_clust)
    }

    if(merge_singles){
      members <-  as.vector(cfg$membership)
      cluster_sizes <- table(cfg$membership)
      check_singles <- which(cluster_sizes == 1)
      if(length(check_singles) > 0){
        ind_singles <- as.vector(which(cluster_sizes == 1))
        ind_multiples <- as.vector(which(cluster_sizes > 1))
        multiples_length <- seq(1:length(ind_multiples))
        for (i in 1:length(ind_multiples)) {
          int_clust <- ind_multiples[i]
          ind_clust <- which(members == int_clust)
          members[ind_clust] <- multiples_length[i]
        }
        members[which(members > max(multiples_length))] <- max(multiples_length)+1
        all_new_cs <- sort(unique(members))
        comm_ind <- vector(mode = "list", length = max(multiples_length)+1)
        for (i in 1:length(all_new_cs)) {
          new_cluster <- V(tmp_net)[which(members == all_new_cs[i])]
          new_cluster <- induced_subgraph(tmp_net, new_cluster)
          comm_ind[[i]] <- new_cluster
        }
        cfg$membership <- members
      }
    }

    cols_clust <- gg_color_hue(length(unique(cfg$membership)))
    cols_clust2 <- adjustcolor(cols_clust, alpha.f = 1)
    clusts <- as.vector(cfg$membership)

    for (i in 1:length(cols_clust)) {
      cl <- cols_clust[i]
      clusts[which(clusts == i)] <- cl
    }

    V(tmp_net)$color_clusters <- clusts

    ###############################################################################################
    ##### Calc stats #####
    ###############################################################################################

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
    if(nrow(tmp) > 1){
      tmp2 <- apply(tmp[ , 1:2 ] , 1 , paste , collapse = "-" )
      names(edg_rank) <- tmp2
    }

    edg_rank <- sort(edg_rank, decreasing = T)

    hs <- hub_score(tmp_net, weights=NA)$vector #Outgoing (ligands)
    hs_srt <- sort(hs, decreasing = T)

    as <- authority_score(tmp_net, weights=NA)$vector # Incoming (receptors)
    as_srt <- sort(as, decreasing = T)

    cluster_edge_btw <- vector(mode = "list", length = length(unique(cfg$membership)))
    cluster_edge_hub <- vector(mode = "list", length = length(unique(cfg$membership)))
    cluster_node_authority <- vector(mode = "list", length = length(unique(cfg$membership)))

    sizes <- rep(0, length(V(tmp_net)))

    for (i in 1:length(unique(cfg$membership))) {
      int <- which(cfg$membership == i)
      vals <- vert_rank[int]
      cluster_edge_hub[[i]] <- sort(hs[int], decreasing = T)
      cluster_node_authority[[i]] <- sort(as[int], decreasing = T)
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

    V(tmp_net)$size <- sizes

    ###############################################################################################
    ##### Format outputs #####
    ###############################################################################################

    E(tmp_net)$width <- E(tmp_net)$width*mult
    V(tmp_net)$size <- V(tmp_net)$size*mult

    pdf(paste0(prefix, "Analyzed_Network.pdf"), h = h, w = w, useDingbats = FALSE)
    plot(tmp_net, layout = l, rescale = TRUE,
         vertex.color = V(tmp_net)$color_celltype,
         edge.color = E(tmp_net)$color_celltype)
    dev.off()

    pdf(paste0(prefix, "Analyzed_Network_noname.pdf"), h = h, w = w, useDingbats = FALSE)
    plot(tmp_net, layout = l, rescale = TRUE,
         vertex.color = V(tmp_net)$color_celltype,
         edge.color = E(tmp_net)$color_celltype,
         vertex.label = V(tmp_net)$name_blank)
    dev.off()

    pdf(paste0(prefix, "Analyzed_Network_communities.pdf"), h = h, w = w, useDingbats = FALSE)
    plot(tmp_net, layout = l, rescale = TRUE,
         mark.groups = cfg,
         vertex.color = V(tmp_net)$color_celltype,
         edge.color = E(tmp_net)$color_celltype)
    dev.off()

    pdf(paste0(prefix, "Analyzed_Network_communities_noname.pdf"), h = h, w = w, useDingbats = FALSE)
    plot(tmp_net, layout = l, rescale = TRUE,
         mark.groups = cfg,
         vertex.color = V(tmp_net)$color_clusters,
         edge.color = E(tmp_net)$color_celltype,
         vertex.label = V(tmp_net)$name_blank)
    dev.off()

    if(!is.null(E(tmp_net)$color_compare)){
      pdf(paste0(prefix, "Analyzed_Network_compared.pdf"), h = h, w = w, useDingbats = FALSE)
      plot(tmp_net, layout = l, rescale = TRUE,
           vertex.color = V(tmp_net)$color_celltype,
           edge.color = E(tmp_net)$color_compare)
      dev.off()
    }

    if(!is.null(E(tmp_net)$color_compare)){
      pdf(paste0(prefix, "Analyzed_Network_compared_noname.pdf"), h = h, w = w, useDingbats = FALSE)
      plot(tmp_net, layout = l, rescale = TRUE,
           vertex.color = V(tmp_net)$color_celltype,
           edge.color = E(tmp_net)$color_compare,
           vertex.label = V(tmp_net)$name_blank)
      dev.off()
    }

    if(!is.null(E(tmp_net)$color_compare)){
      pdf(paste0(prefix, "Analyzed_Network_communities_compare.pdf"), h = h, w = w, useDingbats = FALSE)
      plot(tmp_net, layout = l, rescale = TRUE,
           mark.groups = cfg,
           vertex.color = V(tmp_net)$color_clusters,
           edge.color = E(tmp_net)$color_compare)
      dev.off()
    }

    if(!is.null(E(tmp_net)$color_compare)){
      pdf(paste0(prefix, "Analyzed_Network_communities_compare_noname.pdf"), h = h, w = w, useDingbats = FALSE)
      plot(tmp_net, layout = l, rescale = TRUE,
           mark.groups = cfg,
           vertex.color = V(tmp_net)$color_clusters,
           edge.color = E(tmp_net)$color_compare,
           vertex.label = V(tmp_net)$name_blank)
      dev.off()
    }

    pdf(paste0(prefix, "Analyzed_Network_communities_crossing_noname.pdf"), h = h, w = w, useDingbats = FALSE)
    E(tmp_net)$color_crossing <- "gray"
    E(tmp_net)$color_crossing[crossing(cfg,tmp_net)] <- "red"
    plot(tmp_net, layout = l, rescale = TRUE,
         mark.groups = cfg,
         vertex.color = V(tmp_net)$color_clusters,
         edge.color = E(tmp_net)$color_crossing,
         vertex.label = V(tmp_net)$name_blank)
    dev.off()


    pdf(paste0(prefix, "Analyzed_Network_communities_crossing.pdf"), h = h, w = w, useDingbats = FALSE)
    keep_name <- unique(ends(tmp_net, names(cs2)[which(cs2 == TRUE)], names =  T))
    V(tmp_net)$name_crossing <- V(tmp_net)$name
    V(tmp_net)$name_crossing[-match(keep_name, V(tmp_net)$name_crossing)] <- ""
    plot(tmp_net, layout = l, rescale = TRUE,
         mark.groups = cfg,
         vertex.color = V(tmp_net)$color_clusters,
         edge.color = E(tmp_net)$color_crossing,
         vertex.label = V(tmp_net)$name_crossing)
    dev.off()

    pdf(paste0(prefix, "Analyzed_Network_communities_centers_named.pdf"), h = h, w = w, useDingbats = FALSE)
    keep_center <- match(names(cluster_centers), names(V(tmp_net)))
    V(tmp_net)$name_center <-V(tmp_net)$name
    V(tmp_net)$name_center[-keep_center] <- ""
    plot(tmp_net, layout = l, rescale = TRUE,
         mark.groups = cfg,
         vertex.color = V(tmp_net)$color_clusters,
         edge.color = E(tmp_net)$color_crossing,
         vertex.label = V(tmp_net)$name_center)
    dev.off()


    pdf(paste0(prefix, "Analyzed_Network_communities_numbered.pdf"), h = h, w = w, useDingbats = FALSE)
    V(tmp_net)$name_clusters <- cfg$membership
    plot(tmp_net, layout = l, rescale = TRUE,
         mark.groups = cfg,
         vertex.color = V(tmp_net)$color_clusters,
         edge.color = E(tmp_net)$color_crossing,
         vertex.label = V(tmp_net)$name_clusters)
    dev.off()


    ###############################################################################################
    ##### Interactive #####
    ###############################################################################################

    nodes <- igraph::as_data_frame(tmp_net, what = "vertices")
    links <- igraph::as_data_frame(tmp_net, what = "edges")
    if(dim(links)[1] > 0){
    links$arrows <- "to"

    colnames(nodes)[1] <- "id"
    nodes <- nodes[,c("id", "color_celltype")]

    nodes$color <-  V(tmp_net)$color_celltype
    links$color <- E(tmp_net)$color_celltype

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
    } else {
      vit_net <- NULL
    }

    ###############################################################################################
    ##### Outputs #####
    ###############################################################################################

    results <- vector(mode = "list", length = 10)
    results[[1]] <- deg
    results[[2]] <- deg_rec
    results[[3]] <- vert_rank_srt
    results[[4]] <- as_srt
    results[[5]] <- deg_lig
    results[[6]] <- edg_rank
    results[[7]] <- hs_srt
    results[[8]] <- cs2
    results[[9]] <- cfg
    results[[10]] <- comm_ind
    results[[11]] <- cluster_edge_btw
    results[[12]] <- cluster_edge_hub
    results[[13]] <- cluster_node_authority
    results[[14]] <- vit_net
    results[[15]] <- tmp_net
    results[[16]] <- l

    names(results) <- c("Degree", "Node_degree", "Node_betweeness", "Node_authority",
                        "Edge_degree", "Edge_betweeness", "Edge_hub", "Crossing", "Clusters_Results", "Clusters_individual",
                        "Clusters_betweeness", "Clusters_edge_hub","Clusters_node_authority", "Interactive", "igraph_Network", "layout")
    return(results)
}






