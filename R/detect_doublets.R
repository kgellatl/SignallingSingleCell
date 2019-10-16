#' Calculate UPM expression values across pData values
#'
#' This will calculate UMIs per million UPM expression values across pData columns. Useful for heatmaps and networking analysis.
#'
#' @param input the input ex_sc.
#' @param id_by the pData variable to operate on for doublet detection
#' @param num_markers the number of markers to calculate similarity to other cells
#' @param remove the percentage of data to remove

#' @export
#' @details
#' Utilize information stored in pData to control the plot display.
#' @examples
#' detect_doublets()

detect_doublets <- function(input, id_by, num_markers, remove){
  questionable <- c()
  search <- paste0("marker_score_", id_by)
  ind <- grep(search,(colnames(fData(input))))
  if(length(ind) == 0 ){
    stop("Run id_markers on the id_by provided here first!")
  }
  marks <- return_markers(input, return_by = id_by, num_markers = num_markers)
  marks_all <- unique(unlist(marks))
  tmp <- input[marks_all,]
  groups <- sort(unique(pData(input)[,id_by]))
  for (i in 1:length(groups)) {
    int_group <- groups[i]
    out_groups <- groups[-i]
    ind <- grep(int_group,pData(input)[,id_by])
    group_cells <- tmp[,ind]
    int_genes <- marks[[grep(paste0("^", int_group, "$"), matrix(c(unlist(strsplit(names(marks), "_"))), ncol = 2, byrow = T)[,1])]]
    internal_sums <- apply(exprs(group_cells)[int_genes,],2,sum)
    # if(min(internal_sums) == 0){
    #   stop("Some Cells have zero expression, increase num_markers")
    # }
    for (j in 1:length(out_groups)) {
      out_group <- out_groups[j]
      out_genes <- marks[[grep(paste0("^", out_group, "$"), matrix(c(unlist(strsplit(names(marks), "_"))), ncol = 2, byrow = T)[,1])]]
      external_sums <- apply(exprs(group_cells)[out_genes,],2,sum)
      ratio <- external_sums/internal_sums
      questionable <- c(questionable, ratio)
    }
  }
  # vals <- matrix(ncol = length(marks)-1, nrow = ncol(input))
  # rownames(vals) <- colnames(input)
  #
  # alerts <- floor(ncol(input) / 20)
  # alerts <-  as.numeric(alerts)
  # when_alert <- seq(from = 1, to = ncol(input), by = alerts)
  # when_alert <- when_alert[-1]
  #
  # for (i in 1:nrow(vals)) {
  #   if(i %in% when_alert){
  #     ind <- match(i, when_alert)
  #     print(paste0(5*ind, " Percent Done"))
  #   }
  #   int <- rownames(vals)[i]
  #   ind <- grep(int, names(questionable))
  #   thevals <- questionable[ind]
  #   vals[i,] <- thevals
  # }
  pData(input)$Doublets <- "Good"
  bad <- questionable
  very_bad1 <- which(is.na(bad))
  very_bad2 <- which(!is.finite(bad))
  removed <- c(very_bad1, very_bad2)
  very_bad <- c(unique(names(c(very_bad1, very_bad2))))
  pData(input)[very_bad,]$Doublets <- "Questionable"
  bad <- bad[-removed]
  outlier <- sort(-bad)
  num_remove <- round(remove*ncol(input))
  outlier <- outlier[(1:num_remove)]
  total <- unique(c(names(outlier),very_bad))
  pData(input)[total,]$Doublets <- "Doublet"
  return(input)
}
