#' This will create a violin plot
#'
#' This will plot a given gene via violin plot
#'
#' @param input The input data
#' @param title The title
#' @param gene if provided will color_by the gene
#' @param color_by a pData variable
#' @param facet_by a pData variable
#' @param plot_mean If true will create a secondary axis to plot the mean expression value
#' @param ncol How many columns if faceting
#' @param size The size of the points
#' @param colors What colors to utilize for categorial data. Be sure it is of the proper length!
#' @export
#' @details
#' Utilize information stored in pData to control the plot display.
#' @examples
#' plot_violin(input, title = "Actb across clusters", gene = "Actb", color_by = "Timepoint", facet_by = "Cluster", size = 1, ncol = 3)

plot_violin <- function(input, title = "", color_by, gene, facet_by = "NA", ncol = "NA", size = 1, colors = "NA", plot_mean = TRUE, theme = "classic"){
  if(facet_by == "NA"){
    geneColored1 <- pData(input)[,c("x", "y", color_by)]
    geneColored1 <- cbind(geneColored1, log2(exprs(input)[gene,]+2)-1)
    colnames(geneColored1)<-c("x","y", color_by, gene)
  } else {
    geneColored1 <- pData(input)[,c("x", "y", color_by, facet_by)]
    geneColored1 <- cbind(geneColored1, log2(exprs(input)[gene,]+2)-1)
    colnames(geneColored1)<-c("x","y", color_by, facet_by, gene)
  }
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  cols <- gg_color_hue(length(unique(pData(input)[,color_by])))
  label.n = function(x) {
    return(c(y=-0.1, label=length(x)))
  }
  fracSC = function(x){
    return(c(y = -.2, label = round(mean(x), 2)))
  }
  colnames(geneColored1) <- gsub("-", "", colnames(geneColored1))
  gene <- gsub("-", "", gene)
  geneColored1$frac <- 0
  frac_cells <- "frac"
  geneColored1[which(geneColored1[,gene] >0),"frac"] <- 1
  g <- ggplot(geneColored1)
  if(theme == "bw") {
    g <- g + theme_bw();
  } else {
    g <- g + theme_classic()
  }
  if(title == ""){
    title <- gene
    g <- g + labs(title= title, y = gene)

  } else {
    g <- g + labs(title= title, y = gene)
  }
  g <- g + theme(plot.title = element_text(size = 20), axis.title = element_text(size = 10), legend.title = element_text(size = 15), legend.text=element_text(size=10))
  g <- g + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  g <- g + geom_jitter(aes_string(x=color_by, y=gene, col = color_by), width = 0.2, size = size, alpha = 0.25)
  g <- g + geom_violin(aes_string(x=color_by, y=gene, col = color_by), trim = T, fill = NA)
  g <- g + stat_summary(aes_string(x = color_by, y = gene), fun.data = label.n, fun.y = "mean", colour = "black", geom = "text", size=3)
  g <- g + stat_summary(aes_string(x = color_by, y = frac_cells), fun.data = fracSC, fun.y = "mean", colour = "black", geom = "text", size=3)

  #####
  if(plot_mean == TRUE){
    breaks <- sort(unique(geneColored1[,color_by]))
    means <- c()
    for (i in 1:length(breaks)) {
      means <- c(means, mean(geneColored1[grep(breaks[i], geneColored1[,color_by]),gene]))
    }
    summary_dat <- as.data.frame(matrix(c(breaks, means), ncol = 2))
    summary_dat[,2] <- as.numeric(as.character(summary_dat[,2]))
    colnames(summary_dat) <- c(color_by, gene)
    scale <- max(geneColored1[,gene] / max(summary_dat[,2]))
    summary_dat[,2] <- summary_dat[,2]*(scale*.5)
    if(facet_by == "NA"){
      g <- g + geom_point(data = summary_dat, aes_string(x = color_by, y = gene))
      g <- g + scale_y_continuous(sec.axis = sec_axis(~./(scale*.5), name = "Mean Expression"))
    }
  }

  #####

  if(facet_by != "NA"){
    if(plot_mean == TRUE){
      breaks <- sort(unique(geneColored1[,color_by]))
      facets <- sort(unique(geneColored1[,facet_by]))
      all_dat <- expand.grid(breaks, facets)
      all_dat <- as.data.frame(all_dat)
      all_dat[,3] <- 0
      rrow <- c()
      for (i in 1:nrow(all_dat)) {
        ind1 <- which(geneColored1[,color_by] == all_dat[i,1])
        ind2 <- which(geneColored1[,facet_by] == all_dat[i,2])
        ind <- intersect(ind1, ind2)
        val <- mean(geneColored1[ind,gene])
        if(is.nan(val) == TRUE){
          rrow <- c(rrow, i)
        } else {
          all_dat[i,3] <- val
        }
      }
      if(length(rrow > 0)){
        all_dat <- all_dat[-rrow,]
      }
      all_dat[,3] <- as.numeric(as.character(all_dat[,3]))
      colnames(all_dat) <- c(color_by, facet_by, gene)
      scale <- max(geneColored1[,gene] / max(all_dat[,3]))
      all_dat[,3] <- all_dat[,3]*(scale*.5)
      g <- g + geom_point(data = all_dat, aes_string(x = color_by, y = gene))
      g <- g + scale_y_continuous(sec.axis = sec_axis(~./(scale*.5), name = "Mean Expression"))
    }
    if(ncol != "NA"){
      g <- g +  facet_wrap(facets = reformulate(facet_by), ncol = ncol, scales = "free_x")
    } else {
      g <- g +  facet_grid(facets = reformulate(facet_by), scales = "free_x", space = "free_x")
    }
  }
  #####
  g <- g + theme(axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank())
  return(g)
}

