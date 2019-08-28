#' Calculate Library Size
#'
#' This function will count the number of reads / UMIs per cell
#'
#' @param input the input DIRECTORY
#' @export
#' @details
#' This will calculate the total number of UMIs on a per cell basis.
#' @examples
#' ex_sc_example <- calc_libsize(input = ex_sc_example)

plot_volcano <- function(de_path, de_file, title = NA, fdr_cut = 0.001, logfc_cut = 0.75){
  volcano_dat <- read.delim(paste0(de_path, de_file))

  ind_fc_up <- which(volcano_dat$logFC > logfc_cut)
  ind_fc_down <- which(volcano_dat$logFC < -logfc_cut)

  ind_fdr <- which(volcano_dat$FDR < fdr_cut)

  ind_fc_up <- intersect(ind_fdr, ind_fc_up)
  ind_fc_down <- intersect(ind_fdr, ind_fc_down)


  volcano_dat$color <- "black"
  volcano_dat[ind_fc_up,"color"] <- "red"
  volcano_dat[ind_fc_down,"color"] <- "blue"

  volcano_dat$FDR[which(volcano_dat$FDR == 0)] <- min(volcano_dat$FDR[which(volcano_dat$FDR != 0)])
  volcano_dat$y <- -log10(volcano_dat$FDR)
  volcano_dat$x <- volcano_dat$logFC

  volcano_dat$label <- rownames(volcano_dat)
  keep_labels <- intersect(ind_fdr, c(ind_fc_down, ind_fc_up))
  volcano_dat[-keep_labels,"label"] <- ""

  g <- ggplot(volcano_dat)
  g <- g + geom_point(aes(x = volcano_dat$x, y = volcano_dat$y), colour = volcano_dat$color)
  g <- g + ggrepel::geom_text_repel(data = volcano_dat[keep_labels,], mapping = ggplot2::aes(x = x , y = y, label = label) )
  g <- g + theme_classic()
  g <- g + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  g <- g + xlab("log2(FC)")
  g <- g + ylab("-log10(FDR)")
  if(!is.na(title)){
    g <- g + labs(title= title)
  } else {
    g <- g + labs(title= de_file)
  }
  plot(g)
}
