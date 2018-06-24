#' Identifies all R / L interactions
#'
#' This function will map all RL interactions
#'
#' @param input the input ex_sc
#' @param nodes the pData variable used in calc_agg_bulk that will be used to place nodes (such as cluster or cell type)
#' @param group_by the pData columns calc_agg_bulk was calculated on to split the networks
#' into independent networks (such as an experimental condition)
#' @param weight_by_proportion if true will multiple ligand mean expression within a node by the proportion of that group
#' @export
#' @details
#' This will use the calc_agg_bulk results to ID networks
#' @examples
#' ex_sc_example <- id_rl(input = ex_sc_example)

calc_rl_network <- function(input, nodes, group_by = FALSE, weight_by_proportion = FALSE, print_progress = TRUE){
  ##### Get all Receptor Ligand Pairs expressed in the data into long format
  if(print_progress == TRUE){
    print("Getting all RL Pairs")
  }
  all_pairs <- which(!is.na(fData(input)[,"networks_ligs_to_receptor"]))
  all_pairs <- fData(input)[all_pairs,]
  all_pairs <- all_pairs[,grep("networks_", colnames(all_pairs))]
  ligs <- strsplit(all_pairs$networks_ligs_to_receptor, "_")
  all_pairs_long <- matrix()
  all_pairs_long <- as.data.frame(all_pairs_long)
  for (i in 1:length(ligs)) {
    ligand <- ligs[[i]]
    receptor <- rep(rownames(all_pairs)[i], length(ligand))
    dat <- cbind(receptor, ligand)
    if(i == 1){
      all_pairs_long <- cbind(all_pairs_long, dat)
      all_pairs_long <- all_pairs_long[,2:3]
    } else {
      all_pairs_long <- rbind(all_pairs_long, dat )
    }
  }
  ##### Extract the expression data for the expressed ligand receptor interactions
  if(print_progress == TRUE){
    print("Extracting RL Pairs Expression Information")
  }
  all_expr <- fData(input)[,grep("_bulk", colnames(fData(input)))]
  receptor_bulks <- as.character(unique(all_pairs_long$receptor))
  ligand_bulks <- as.character(unique(all_pairs_long$ligand))
  all_expr <- all_expr[c(receptor_bulks, ligand_bulks),]
  ##### Construct the possible nodes in the network
  if(print_progress == TRUE){
    print("Constructing Nodes")
  }
  node <- unique(pData(input)[,nodes])
  nod <- vector(mode = "list", 2)
  nod[[1]] <- node
  nod[[2]] <- node
  bulks <- expand.grid(nod, stringsAsFactors = FALSE)
  interactions <- matrix(c(bulks$Var2, bulks$Var1), ncol= 2)
  ##### Break those nodes by a variable across break by
  if(group_by != FALSE){
    breaks <- sort(unique(pData(input)[,group_by]))
    col1 <- rep(interactions[,1], length(breaks))
    col2 <- rep(interactions[,2], length(breaks))
    col3 <- c()
    for (i in 1:length(breaks)) {
      bras <- breaks[i]
      bras <- rep(bras, nrow(interactions))
      col3 <- c(col3, bras)
    }
    interactions <- matrix(c(col1, col2, col3), ncol = 3)
  }
  ##### Now construct the full network
  if(print_progress == TRUE){
    print("Constructing Interactome")
  }
  full_network <- data.frame()
  for (i in 1:nrow(interactions)) {
    inter <- interactions[i,]
    dat <- rep(inter, nrow(all_pairs_long))
    if(group_by != FALSE){
      dat <- matrix(dat, ncol = 3, byrow = T)
      dat <- as.data.frame(dat)
      dat[,4] <- all_pairs_long[,1]
      dat[,5] <- all_pairs_long[,2]
      full_network <- rbind(full_network, dat)
    } else {
      dat <- matrix(dat, ncol = 2, byrow = T)
      dat <- as.data.frame(dat)
      dat[,3] <- all_pairs_long[,1]
      dat[,4] <- all_pairs_long[,2]
      full_network <- rbind(full_network, dat)
    }
  }
  if(group_by != FALSE){
    full_network <- full_network[,c(1,5,2,4,3)]
    colnames(full_network) <- c(nodes, "Ligand", nodes, "Receptor", group_by)
  } else {
    full_network <- full_network[,c(1,4,2,3)]
    colnames(full_network) <- c(nodes, "Ligand", nodes, "Receptor")
  }
  full_network[,1] <- as.character(full_network[,1])
  full_network[,2] <- as.character(full_network[,2])
  full_network[,3] <- as.character(full_network[,3])
  full_network[,4] <- as.character(full_network[,4])
  if(group_by != FALSE){
    full_network[,5] <- as.character(full_network[,5])
  }
  ##### Write in the expression values
  if(print_progress == TRUE){
    print("Calculating Interactions")
    alerts <- c()
    for (i in 1:20) {
      printi <- floor(nrow(full_network)/20)*i
      alerts <- c(alerts, printi)
    }
  }
  full_network$Ligand_expression <- 0
  full_network$Receptor_expression <- 0
  remove <- c()
  for (i in 1:nrow(full_network)) {
    if(print_progress == TRUE){
      if(i %in% alerts){
        ind <- match(i, alerts)
        print(paste0(ind*5, "% Complete"))
      }
    }
    int <- full_network[i,]
    int_exp <- all_expr[c(int$Ligand, int$Receptor),]
    ctr <- int[,3]
    ctl <- int[,1]
    lig <- int$Ligand
    rec <- int$Receptor
    if(group_by != FALSE){
      ctr <- paste0(int[,group_by], "_", ctr)
      ctl <- paste0(int[,group_by], "_", ctl)
    } else {
      ctr <- paste0(ctr)
      ctl <- paste0(ctl)
    }
    rec_ex <- int_exp[rec,grep(ctr, colnames(int_exp))]
    lig_ex <- int_exp[lig,grep(ctl, colnames(int_exp))]
    if(weight_by_proportion == TRUE){
      to_parse <- colnames(int_exp)[grep(ctl, colnames(int_exp))]
      to_parse <- unlist(strsplit(to_parse, "_"))
      to_parse <- to_parse[match("percent", to_parse)+1]
      to_parse <- as.numeric(to_parse)/100
      lig_ex <- lig_ex*to_parse
    }
    full_network[i, "Ligand_expression"] <- lig_ex
    full_network[i, "Receptor_expression"] <- rec_ex
    if(full_network[i, "Ligand_expression"] == 0 || full_network[i, "Receptor_expression"] == 0 ){
      remove <- c(remove, i)
    }
  }
  full_network <- full_network[-remove,]
  full_network$Ligand_expression <- full_network$Ligand_expression*1E6
  full_network$Receptor_expression <- full_network$Receptor_expression*1E6
  full_network$Connection <- full_network$Ligand_expression*full_network$Receptor_expression
  ##### Summarize Interactions #####
  if(print_progress == TRUE){
    print("Summarizing Interactions")
  }
  vec1 <- full_network[,1]
  vec2 <- full_network[,3]
  if(group_by != FALSE){
    vec3 <- full_network[,5]
    dat <- matrix(c(vec1, vec2, vec3), ncol = 3)
  } else {
    dat <- matrix(c(vec1, vec2), ncol = 2)
  }
  dat <- as.data.frame(dat)
  ##### Number of connections #####
  summary <- plyr::count(dat)
  summary[,1:2] <- data.frame(lapply(summary[,1:2], as.character), stringsAsFactors=FALSE)
  ##### Number of connections divided by num genes #####
  if(group_by != FALSE){
    tmp_dat <- matrix(c(as.character(summary$V3), as.character(summary$V1)), ncol = 2)
    tmp_dat <- as.data.frame(tmp_dat)
    tmp_dat <- data.frame(lapply(tmp_dat, as.character), stringsAsFactors=FALSE)
    grp <- c()
    for (i in 1:nrow(tmp_dat)) {
      topaste <- tmp_dat[i,1:2]
      topaste <- as.character(as.vector(topaste))
      int <- paste0(topaste, collapse = "_")
      grp <- c(grp, int)
    }
    summary$freq_frac <- NA
    for (i in 1:length(grp)) {
      int <- grp[i]
      int <- colnames(fData(input))[grep(int, colnames(fData(input)))]
      int <- unlist(strsplit(int, "_"))
      pos <- match("genes", int)
      frac <- summary$freq[i]/as.numeric(int[pos+1])
      summary$freq_frac[i] <- frac
    }
  } else {
    tmp_dat <- data.frame(lapply(summary, as.character), stringsAsFactors=FALSE)
    grp <- unique(tmp_dat$V1)
    summary$freq_frac <- NA
    for (i in 1:length(grp)) {
      int <- grp[i]
      int <- colnames(fData(input))[grep(int, colnames(fData(input)))]
      int <- unlist(strsplit(int, "_"))
      pos <- match("genes", int)
      frac <- summary$freq[grep(grp[i], summary$V1)]/as.numeric(int[pos+1])
      summary$freq_frac[grep(grp[i], summary$V1)] <- frac
    }
  }
  ##### Number of connections divided by total outgoing connections #####
  summary$prop_freq <- NA
  if(group_by != FALSE){
    nodes <- unique(summary$V1)
    conds <- unique(summary$V3)
    mat <- expand.grid(nodes, conds)
    for (i in 1:nrow(mat)) {
      int <- mat[i,]
      ind1 <- grep(int$Var1, summary$V1)
      ind2 <- grep(int$Var2, summary$V3)
      pos <- intersect(ind1, ind2)
      tot <- sum(summary[pos,"freq"])
      prop <- summary[pos,"freq"]/tot
      summary[pos,"prop_freq"] <- prop
    }
  } else {
    tmp_dat <- data.frame(lapply(summary, as.character), stringsAsFactors=FALSE)
    grp <- unique(tmp_dat$V1)
    summary$prop_freq <- NA
    for (i in 1:length(grp)) {
      total <- sum(summary[grep(grp[i], summary$V1),"freq"])
      summary[grep(grp[i], summary$V1),"prop_freq"] <- summary[grep(grp[i], summary$V1),"freq"]/total
    }
  }
  ##### Rank Connections by L * R coarse (within cell type A to anyone!) #####
  full_network$Connection_rank_coarse <- NA
  full_network$Connection_Z_coarse <- NA
  if(print_progress == TRUE){
    print("Calculating Coarse Ranks")
  }
  if(group_by != FALSE){
    full_network$Connection_rank_coarse_grouped <- NA
    full_network$Connection_Z_coarse_grouped <- NA
    tmpdat <- expand.grid(unique(summary$V1), unique(summary$V3))
    for (i in 1:nrow(tmpdat)) {
      int <- tmpdat[i,]
      ind1 <- grep(int$Var1, full_network[,1])
      ind2 <- grep(int$Var2, full_network[,5])
      ind <- intersect(ind1,ind2)
      vals <- full_network[ind,"Connection"]
      zval <- scale(vals)
      vals <- rank(-vals)
      full_network[ind,"Connection_rank_coarse_grouped"] <- vals
      full_network[ind,"Connection_Z_coarse_grouped"] <- zval

    }
  }
  tmpdat <- unique(summary$V1)
  for (i in 1:length(tmpdat)) {
    int <- tmpdat[i]
    ind1 <- grep(int, full_network[,1])
    vals <- full_network[ind1,"Connection"]
    zval <- scale(vals)
    vals <- rank(-vals)
    full_network[ind1,"Connection_rank_coarse"] <- vals
    full_network[ind1,"Connection_Z_coarse"] <- zval
  }
  ##### Rank Connections by L * R fine (within cell type A to cell type B!) #####
  full_network$Connection_rank_fine <- NA
  full_network$Connection_Z_fine <- NA
  if(print_progress == TRUE){
    print("Calculating Fine Ranks")
  }
  tmp_dat <- data.frame(lapply(summary, as.character), stringsAsFactors=FALSE)
  if(group_by != FALSE){
    full_network$Connection_rank_fine_grouped <- NA
    full_network$Connection_Z_fine_grouped <- NA
    for (i in 1:nrow(tmp_dat)) {
      int <- tmp_dat[i,]
      ind1 <- grep(int$V1, full_network[,1])
      ind2 <- grep(int$V2, full_network[,3])
      ind3 <- grep(int$V3, full_network[,5])
      ind <- intersect(intersect(ind1,ind2),ind3)
      vals <- full_network[ind,"Connection"]
      vals <- full_network[ind,"Connection"]
      zval <- scale(vals)
      vals <- rank(-vals)
      full_network[ind,"Connection_rank_fine_grouped"] <- vals
      full_network[ind,"Connection_Z_fine_grouped"] <- zval
    }
  }
  for (i in 1:nrow(tmp_dat)) {
    int <- tmp_dat[i,]
    ind1 <- grep(int$V1, full_network[,1])
    ind2 <- grep(int$V2, full_network[,3])
    ind <- intersect(ind1,ind2)
    vals <- full_network[ind,"Connection"]
    zval <- scale(vals)
    vals <- rank(-vals)
    full_network[ind,"Connection_rank_fine"] <- vals
    full_network[ind,"Connection_Z_fine"] <- vals
  }
  ##### Rank Connections by Zscore #####
  full_network$Zscores_genes <- NA
  if(print_progress == TRUE){
    print("Calculating Z Scores Grouped")
    alerts <- c()
    for (i in 1:20) {
      printi <- floor(nrow(all_pairs_long)/20)*i
      alerts <- c(alerts, printi)
    }
  }
  tmp_dat <- data.frame(lapply(all_pairs_long, as.character), stringsAsFactors=FALSE)
  if(group_by != FALSE){
    full_network$Zscores_genes_grouped <- NA
    for (i in 1:nrow(tmp_dat)) {
      if(print_progress == TRUE){
        if(i %in% alerts){
          ind <- match(i, alerts)
          print(paste0(ind*5, "% Complete"))
        }
      }
      int <- tmp_dat[i,]
      ind1 <- grep(int$receptor, full_network$Receptor)
      ind2 <- grep(int$ligand, full_network$Ligand)
      ind <- intersect(ind1, ind2)
      tmp <- full_network[ind,]
      vals <- scale(tmp$Connection)
      full_network[ind,"Zscores_genes_grouped"] <- vals
      for (j in 1:length(breaks)) {
        int <- breaks[j]
        ind3 <- grep(int, tmp[,group_by])
        tmp2 <- tmp[ind3,]
        vals <- scale(tmp2$Connection)
        ind4 <- match(rownames(tmp2), rownames(full_network))
        full_network[ind4,"Zscores_genes_grouped"] <- vals
      }
    }
  }
  if(print_progress == TRUE){
    print("Calculating Z Scores")
  }
  for (i in 1:nrow(tmp_dat)) {
    if(print_progress == TRUE){
      if(i %in% alerts){
        ind <- match(i, alerts)
        print(paste0(ind*5, "% Complete"))
      }
    }
    int <- tmp_dat[i,]
    ind1 <- grep(int$receptor, full_network$Receptor)
    ind2 <- grep(int$ligand, full_network$Ligand)
    ind <- intersect(ind1, ind2)
    tmp <- full_network[ind,]
    vals <- scale(tmp$Connection)
    full_network[ind,"Zscores_genes"] <- vals
  }
  #####
  if(group_by == FALSE){
    colnames(summary) <- c("Lig_produce", "Rec_receive", "num_connections", "fraction_connections", "proportion_connections")
  } else {
    colnames(summary) <- c("Lig_produce", "Rec_receive", group_by, "num_connections", "fraction_connections", "proportion_connections")
  }
  #####
  results <- vector(mode = "list", 2)
  results[[1]] <- summary
  results[[2]] <- full_network
  names(results) <- c("Summary", "full_network")
  return(results)
  #####
}
