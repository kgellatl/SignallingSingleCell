#' Search Gene
#'
#' This function will search for a gene in your data by some input
#'
#' @param input the input ex_sc
#' @param search The keyword or partial word to match
#' @export
#' @details
#' This will calculate total UMIs on a per cell basis.
#' @examples
#' ex_sc_example <- calc_libsize(input = ex_sc_example)

search_gene <- function(input, search){
  if(class(input)=="ExpressionSet"){
    gene_set <- rownames(exprs(input))
    matches <- grep(search, gene_set, value = T)
    return(matches)
  }
  if(length(input$igraph_Network) > 0){
    matches <- names(V(input$igraph_Network))[grep(search, names(V(input$igraph_Network)))]
    return(matches)
  }
  if(length(input$input) > 0){
    matches <- names(V(input$input))[grep(search, names(V(input$input)))]
    return(matches)
  }
}
