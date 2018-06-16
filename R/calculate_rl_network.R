#' Identifies all R / L interactions
#'
#' This function will map all RL interactions
#'
#' @param input the input ex_sc
#' @param nodes the pData variable used in calc_agg_bulk that will be used to place nodes (such as cluster or cell type)
#' @param group_by the pData columns calc_agg_bulk was calculated on to split the netoworks into independent networks
#'(such as an experimental condition)
#' @param weight_by_proportion if true will multiple ligand mean expression within a node by the proportion of that group
#' @export
#' @details
#' This will use the calc_agg_bulk results to ID networks
#' @examples
#' ex_sc_example <- id_rl(input = ex_sc_example)

calculate_rl_network <- function(input, nodes, group_by = FALSE, print_progress = TRUE, weight_by_proportion = TRUE){
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
  if(break_by != FALSE){
    breaks <- sort(unique(pData(input)[,break_by]))
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
    if(break_by != FALSE){
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
  if(break_by != FALSE){
    full_network <- full_network[,c(1,5,2,4,3)]
    colnames(full_network) <- c(nodes, "Ligand", nodes, "Receptor", break_by)
  } else {
    full_network <- full_network[,c(1,3,2,4)]
    colnames(full_network) <- c(nodes, "Ligand", nodes, "Receptor")
  }
  full_network[,1] <- as.character(full_network[,1])
  full_network[,2] <- as.character(full_network[,2])
  full_network[,3] <- as.character(full_network[,3])
  full_network[,4] <- as.character(full_network[,4])
  if(break_by != FALSE){
    full_network[,5] <- as.character(full_network[,5])
  }
  ##### Write in the expression values
  if(print_progress == TRUE){
    print("Calculating Interactions")
  }
  full_network$Ligand_expression <- 0
  full_network$Receptor_expression <- 0
  remove <- c()
  if(print_progress == TRUE){
    alerts <- c()
    for (i in 1:20) {
      printi <- floor(nrow(full_network)/20)*i
      alerts <- c(alerts, printi)
    }
  }
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
    if(break_by != FALSE){
      ctr <- paste0(int[,break_by], "_", ctr)
      ctl <- paste0(int[,break_by], "_", ctl)
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
  plot(density(log10(full_network$Connection)))
  ##### Summarize Interactions #####

  #####
  return(full_network)
  #####
}
