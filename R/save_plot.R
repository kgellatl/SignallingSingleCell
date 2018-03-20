#' Save plot
#'
#' This will save the last ggplot call
#'
#' @param filename The file name to be written
#' @param format File extension can be "pdf", "svg", "tiff", "png", "jpg"
#' @param dpi resolution for raste output types
#' @param dingbats FALSE or TRUE (for pdf only)
#' @param h The height
#' @param w The Width
#' @param units "in", "cm"
#' @param dpi resolution for raster output types
#' @export
#' @details
#' This will save the save ggplot object that was written to the device.
#' @examples
#' save_plot(filename = "Cluster plot", format = "pdf")

save_ggplot <- function(filename = "Rplot", format = "pdf", dingbats = FALSE, dpi = 300, h = 6, w = 6, units = "in"){
  ggplot2::ggsave(filename = paste0(filename, ".", format), height = h, width = w)
}


