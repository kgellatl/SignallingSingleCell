#' Identifies all R / L interactions
#'
#' This function will map all RL interactions
#'
#' @param input the output of calc_rl_network
#' @param filter_type "network" or "DE"
#' @param filter_by the column to summarize by
#' @param DEfolder the output location of findDEgenes
#' @param absol whether or not to take the absolute value of the column
#' @param cutoff value to filter by
#' @param direction >, <, =
#' @param keep if true will retain all edges and annotate them as passing the filter or not
#' @export
#' @details
#' This will use the calc_agg_bulk results to ID networks
#' @examples
#' ex_sc_example <- id_rl(input = ex_sc_example)
#'
filter_rl_network <- function(input, filter_by, filter_type = "network", DEfolder = NULL, cutoff, direction, absol = FALSE,
                              group_by = FALSE, keep = FALSE){
  rownames(input$full_network) <- seq(1:nrow(input$full_network))
  if(filter_type != "DE" && filter_type != "network"){
    stop("Filter types can be 'DE' or 'network'")
  }
  if(filter_type == "network"){
    bools <- c()
    val <- input$full_network[,filter_by]
    if(absol == TRUE){
      val <- abs(val)
    }
    for (j in 1:length(val)) {
      qr <- val[j]
      bool <- f(qr, cutoff)
      bools <- c(bools, bool)
    }
    ind_keep <- which(bools == TRUE)
    val <- input$full_network[ind_keep,]
    if(group_by != FALSE){
      summary <- plyr::count(val[,c(1,3,5)])
      data1 <- apply( summary[ , 1:3 ] , 1 , paste , collapse = "-" )
      data2 <- apply( input$Summary[ , 1:3 ] , 1 , paste , collapse = "-" )
    } else {
      summary <- plyr::count(val[,c(1,3)])
      data1 <- apply( summary[ , 1:2 ] , 1 , paste , collapse = "-" )
      data2 <- apply( input$Summary[ , 1:2 ] , 1 , paste , collapse = "-" )
    }
    input$Summary[,filter_by] <- 0
    colnames(input$Summary)[ncol(input$Summary)] <- paste0(c(filter_by, cutoff, direction), collapse = "_")
    for (i in 1:length(data1)) {
      int <- data1[i]
      ind <- match(int, data2)
      input$Summary[ind,ncol(input$Summary)] <- summary$freq[i]
    }
    if (keep == TRUE) {
      if(is.null(input$full_network$keep)){
        input$full_network$keep <- FALSE
        input$full_network$keep[ind_keep] <- TRUE
      } else {
        keep2 <- which(input$full_network$keep == TRUE)
        ind <- intersect(ind_keep, keep2)
        input$full_network$keep <- FALSE
        input$full_network$keep[ind_keep] <- TRUE
      }
    } else {
      input$full_network <- val
    }
  }
  if(filter_type == "DE"){
    if(is.null(DEfolder)){
      stop("Requires a DE folder from findDEgenes output")
    }
    l1 <- length(filter_by)
    l2 <-length(absol)
    l3 <- length(cutoff)
    l4 <- length(direction)
    if(length(unique(c(l1,l2,l3,l4))) > 1){
      stop("filter_by, absol, cutoff, and direction must be vectors of equal length")
    }

    #####
    # First list all the files in the DE folder and find ones that correspond to the nodes in network
    #####
    filelist <- list.files(DEfolder)
    filelist <- paste0(DEfolder, "/", filelist)
    types <- unique(input$full_network[,1])
    rows_keep <- c()
    for (i in 1:length(filelist)) {
      int <- filelist[i]
      interested <- c()
      for (j in 1:length(types)) {
        match <- grep(types[j], int)
        if(length(match > 0)){
          if(match == 1){
            interested <- types[j]
          }
        }
      }

      #####
      # Find ones that pass the filter criteria
      #####

      if(length(filter_by) > 1){
        full_gene_set <- c()
        for (l in 1:length(filter_by)) {
          int_filt <- filter_by[l]
          f <- function(a, b, op=direction[l]) {
            call <- call(op, a, b)
            result <- eval(call)
            result
          } # ENDS F DEFINITION
          tmp <- read.table(int, sep = "\t", row.names = 1, header = T)
          vec <- tmp[,int_filt]

          if(absol[l] == TRUE){
            vec <- abs(vec)
          }
          bool <- f(vec, cutoff[l])
          int_gene <- rownames(tmp)[bool]
          full_gene_set <- c(full_gene_set, int_gene)
        } # FINISHES FILTER BASED ON PARAMTER

        int_gene <- which(table(full_gene_set) == length(filter_by))
        int_gene <- names(int_gene)
        #####
        # Find now the network rows corresponding to the input celltype
        #####
        indices <-  grep(interested, input$full_network[,1]) # Find from position
        indices2 <-  grep(interested, input$full_network[,3]) # Find to position
        indices <- unique(c(indices, indices2)) ### Is this correct???
        kint2 <- input$full_network[indices,]
        rownames(kint2) <- rownames(input$full_network)[indices]
        rec <- kint2$Receptor
        ligs <- kint2$Ligand
        mfin <- c()
        for (k in 1:length(int_gene)) {
          gene2 <- int_gene[k]
          rec_match <- grep(paste0("^", gene2, "$"), rec)
          lig_match <- grep(paste0("^", gene2, "$"), ligs)
          fin <- c(rec_match, lig_match)
          mfin <- c(mfin, fin)
          mfin <- unique(mfin)
        } # ENDS GENE LOOP
        #####
        # Keep the network rows based on this match
        #####
        genes_keep <- kint2[mfin,]
        rows_keep <- c(rows_keep, rownames(genes_keep))
        rows_keep <- unique(rows_keep)
        rows_keep <- as.numeric(rows_keep)

      } else { # STARTS SINGLE PARAMETER
        ########################################
        ##### SINGLE FILTER BY PARAMETER #####
        ########################################
        f <- function(a, b, op=direction) {
          call <- call(op, a, b)
          result <- eval(call)
          result
        }
        tmp <- read.table(int, sep = "\t", row.names = 1, header = T)
        vec <- tmp[,filter_by]
        if(absol == TRUE){
          vec <- abs(vec)
        }
        bool <- f(vec, cutoff)
        int_gene <- rownames(tmp)[bool]
        #####
        # Find now the network rows corresponding to the input celltype
        #####
        indices <-  grep(interested, input$full_network[,1]) # Find from position
        indices2 <-  grep(interested, input$full_network[,3]) # Find to position
        indices <- unique(c(indices, indices2)) ### Is this correct???
        kint2 <- input$full_network[indices,]
        rownames(kint2) <- rownames(input$full_network)[indices]
        rec <- kint2$Receptor
        ligs <- kint2$Ligand
        mfin <- c()
        for (k in 1:length(int_gene)) {
          gene2 <- int_gene[k]
          rec_match <- grep(paste0("^", gene2, "$"), rec)
          lig_match <- grep(paste0("^", gene2, "$"), ligs)
          fin <- c(rec_match, lig_match)
          mfin <- c(mfin, fin)
          mfin <- unique(mfin)
        } # ENDS GENE LOOP
        #####
        # Keep the network rows based on this match
        #####
        genes_keep <- kint2[mfin,]
        rows_keep <- c(rows_keep, rownames(genes_keep))
        rows_keep <- unique(rows_keep)
        rows_keep <- as.numeric(rows_keep)

      } # ENDS SINGLE PARAMETER LOOP

    } # ENDS FILE LIST LOOP

    if (keep == TRUE) {
      if(is.null(input$full_network$keep)){
        input$full_network$keep <- FALSE
        input$full_network$keep[(rows_keep)] <- TRUE
      } else {
        keep2 <- which(input$full_network$keep == TRUE)
        ind <- intersect((rows_keep), keep2)
        input$full_network$keep <- FALSE
        input$full_network$keep[(rows_keep)] <- TRUE
      }
    } else {
      input$full_network <- input$full_network[sort(rows_keep),]
    }
  }
  return(input)
}

