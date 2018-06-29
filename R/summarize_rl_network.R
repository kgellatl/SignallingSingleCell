#' Identifies all R / L interactions
#'
#' This function will map all RL interactions
#'
#' @param input the output of calc_rl_network
#' @param summarize_by the column of full_network to summarize by
#' @param cutoff value to filter by
#' @param direction >, <, =
#' @export
#' @details
#' This will use the calc_agg_bulk results to ID networks
#' @examples
#' ex_sc_example <- id_rl(input = ex_sc_example)

summarize_rl_network <- function(input, summarize_by, cutoff, direction){
  f <- function(a, b, op=direction) {
    call <- call(op, a, b)
    result <- eval(call)
    result
  }
  bools <- c()
  val <- input$full_network[,summarize_by]
  for (j in 1:length(val)) {
    qr <- val[j]
    bool <- f(qr, cutoff)
    bools <- c(bools, bool)
  }
  ind <- which(bools == TRUE)
  val <- input$full_network[ind,]
  input$full_network <- val


  gene_list <- unique(c(val$Ligand, val$Receptor))
  summary <- plyr::count(val[,c(1,3,5)])
  input$Summary[,summarize_by] <- NA
  input$Summary[,summarize_by] <- summary$freq
  head(input$Summary)
}

