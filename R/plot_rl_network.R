#' Identifies all R / L interactions
#'
#' This function will map all RL interactions
#'
#' @param input the input ex_sc
#' @param nodes the pData variable used in calc_agg_bulk that will be used to place nodes (such as cluster or cell type)
#' @param break_by the pData columns calc_agg_bulk was calculated on to split the netoworks into independent networks
#'(such as an experimental condition)
#' @export
#' @details
#' This will use the calc_agg_bulk results to ID networks
#' @examples
#' ex_sc_example <- id_rl(input = ex_sc_example)

plot_rl_network <- function(input, nodes, break_by = FALSE){
  ##### Get all Receptor Ligand Pairs expressed in the data into long format
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
  all_expr <- fData(input)[,grep("_bulk", colnames(fData(input)))]
  receptor_bulks <- as.character(unique(all_pairs_long$receptor))
  ligand_bulks <- as.character(unique(all_pairs_long$ligand))
  all_expr <- all_expr[c(receptor_bulks, ligand_bulks),]
  ##### Construct the possible nodes in the network
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
  full_network$Ligand_expression <- 0
  full_network$Receptor_expression <- 0
  str(full_network)
  ##### Write in the expression values
  head(all_expr)
  head(full_network)
  for (i in 1:nrow(full_network)) {
    int <- full_network[i,]
    int_exp <-all_expr[c(int$Ligand, int$Receptor),]

  }
}
