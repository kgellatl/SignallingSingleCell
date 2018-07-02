#' Identifies all R / L interactions
#'
#' This function will map all RL interactions
#'
#' @param input the output of calc_rl_network
#' @param filter_by the column of full_network to summarize by
#' @param cutoff value to filter by
#' @param direction >, <, =
#' @export
#' @details
#' This will use the calc_agg_bulk results to ID networks
#' @examples
#' ex_sc_example <- id_rl(input = ex_sc_example)
filter_rl_network <- function(input, filter_by, cutoff, direction, group_by = FALSE){
  f <- function(a, b, op=direction) {
    call <- call(op, a, b)
    result <- eval(call)
    result
  }
  bools <- c()
  val <- input$full_network[,filter_by]
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
  return(input)
}

