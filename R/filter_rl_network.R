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
#' @export
#' @details
#' This will use the calc_agg_bulk results to ID networks
#' @examples
#' ex_sc_example <- id_rl(input = ex_sc_example)
filter_rl_network <- function(input, filter_by, filter_type = "network", DEfolder = NULL, cutoff, direction, absol = FALSE, group_by = FALSE){
  f <- function(a, b, op=direction) {
    call <- call(op, a, b)
    result <- eval(call)
    result
  }
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
    ind <- which(bools == TRUE)
    val <- input$full_network[ind,]
    input$full_network <- val
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
  }
  if(filter_type == "DE"){
    if(is.null(DEfolder)){
      stop("Requires a DE folder from findDEgenes output")
    }
    filelist <- list.files(DEfolder)
    filelist <- paste0(DEfolder, "/", filelist)
    keep_row <- c()
    genes <- c()
    for (i in 1:length(filelist)) {
      int <- filelist[i]
      tmp <- read.table(int, sep = "\t", row.names = 1, header = T)
      vec <- tmp[,filter_by]
      if(absol == TRUE){
        vec <- abs(vec)
      }
      bool <- f(vec, cutoff)
      int_gene <- rownames(tmp)[bool]
      genes <- unique(c(genes, int_gene))
    }
    pairs <- paste0(input$full_network$Ligand, "_", input$full_network$Receptor)
    for (i in 1:length(genes)) {
      int <- genes[i]
      keep <- grep(int, pairs)
      keep_row <- unique(keep_row)
      if(length(keep)>0){
        keep_row <- c(keep_row, keep)
      }
    }
    input$full_network <- input$full_network[sort(keep_row),]
  }
  return(input)
}

