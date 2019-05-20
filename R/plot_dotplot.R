#' This will create a bar plot with dots
#'
#' @param input Bioconductor’s ExpressionSet Class with bulk value stored in fData.
#' @param title The title
#' @param gene to plot the expression level of
#' @param color_by a pData variable
#' @param colors What colors to utilize for categorial data. Be sure it is of the proper length!
#' @param facet_by a pData variable
#' @param point_by a pData variable.
#' @param ncol How many columns if faceting
#' @param text_sizes a vector of title_size, axis_title, axis_text, legend_title, legend_text, facet_text, faults too c(20,10,5,10,5,5)
#' @param theme the plot theme
#' @param number_labels to show the cell numbers and cell percentage of each bar.
#' @param stackratio the overlap of dots.
#' @param dotsize the size of dots.
#' @param bar the weighted mean.
#' @param binwidth average the values when dots are within the range of (max-min)*binwidth.
#' @export
#' @details
#' Utilize information stored in pData to control the plot display. Each point_by as a dot with a bar showing the weighted mean of all point_by dots.
#' @examples
#' plot_dotplot(ex_sc, gene = "ADCY7", color_by = "Skin", facet_by = "subCellType", point_by = "Patient")


plot_dotplot <- function(input, gene, color_by, facet_by = "NA", point_by, title = "", colors = NA, ncol = "NA", number_labels = T, text_sizes = c(20, 10, 5, 10, 5, 5), theme = "classic", alpha = 0.5, stackratio = 0.4, dotsize = 3, bar = T, binwidth = 0.005)
{
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  cols <- gg_color_hue(length(unique(pData(input)[, color_by])))
  geneColored1 <- fData(input)
  ind <- grep("bulk", colnames(fData(input)))
  if (length(ind) == 0) {
    stop("Must calc_agg_bulk() before using this function")
  }
  geneColored1 <- geneColored1[, ind]
  color_bys <- sort(unique(pData(input)[, color_by]))
  if (facet_by == "NA") {
    ind <- grep(color_bys[1], colnames(geneColored1))
    if (length(ind) == 0) {
      stop("Provided color_by argument was not used to calc_agg_bulk()")
    }
  }
  point_bys <- sort(unique(pData(input)[, point_by]))

  geneColored1 <- as.data.frame(t(geneColored1[gene, ]))
  geneColored1[, "tmp_val"] <- sub("_num_.*", "", rownames(geneColored1))

  ncol_tmp <- length(strsplit(geneColored1[1, "tmp_val"], split = "_")[[1]])
  tmpmat <- matrix(unlist(strsplit(geneColored1[, "tmp_val"], split = "_")), byrow = T, ncol = ncol_tmp)
  if (facet_by != "NA") {facet_bys <- sort(unique(pData(input)[, facet_by]))}else{facet_bys <- "NA"}

  ind <- apply(tmpmat, 2, function(x) all(x %in% color_bys) | all(x %in% facet_bys) | all(x %in% point_bys))
  if (!all(ind)) {
    stop("Bulk values are calculated with more variables. You may want to provide a facet_by.")
  }
  tmpmat <- tmpmat[, ind]
  names_tmp <- c()
  for (i in 1:ncol(tmpmat)) {
    ind <- c(all(tmpmat[,i] %in% color_bys), all(tmpmat[,i] %in% facet_bys), all(tmpmat[,i] %in% point_bys))
    names_tmp <- c(names_tmp, c(color_by, facet_by, point_by)[ind])
  }
  colnames(tmpmat) <- names_tmp

  geneColored1 <- cbind(geneColored1, tmpmat)

  for (i in 1:nrow(geneColored1)) {
    int_string <- rownames(geneColored1)[i]
    int_string <- strsplit(int_string, split = "_")[[1]]
    geneColored1$num[i] <- int_string[match("cells", int_string) +
                                        1]
  }
  for (i in 1:nrow(geneColored1)) {
    int_string <- rownames(geneColored1)[i]
    int_string <- strsplit(int_string, split = "_")[[1]]
    geneColored1$frac[i] <- int_string[match("percent",
                                             int_string) + 1]
  }

  geneColored1$weighted_mean <- c()
  geneColored1$num_sum <- c()
  geneColored1$frac_sum <- c()
  if (bar) {
    for (i in 1:length(color_bys)) {
      if (facet_by == "NA") {
        ind <- which(geneColored1[,color_by] == color_bys[i])
        if (length(ind) == 0) {next}
        mean_tmp <- sum(as.numeric(geneColored1[ind,gene])*as.numeric(geneColored1[ind,"num"]))/sum(as.numeric(geneColored1[ind,"num"]))
        geneColored1[ind[1],"weighted_mean"] <- mean_tmp
        geneColored1[ind[1],"num_sum"] <- sum(as.numeric(geneColored1[ind,"num"]))
        geneColored1[ind[1],"frac_sum"] <- sum(as.numeric(geneColored1[ind,"frac"]))
      }else{
        for (j in 1:length(facet_bys)){
          ind <- which(geneColored1[,color_by] == color_bys[i] & geneColored1[,facet_by] == facet_bys[j])
          if (length(ind) == 0) {next}
          mean_tmp <- sum(as.numeric(geneColored1[ind,gene])*as.numeric(geneColored1[ind,"num"]))/sum(as.numeric(geneColored1[ind,"num"]))
          geneColored1[ind[1],"weighted_mean"] <- mean_tmp
          geneColored1[ind[1],"num_sum"] <- sum(as.numeric(geneColored1[ind,"num"]))
          geneColored1[ind[1],"frac_sum"] <- sum(as.numeric(geneColored1[ind,"frac"]))
        }
      }
    }
  }

  genename <- gene

  colnames(geneColored1) <- gsub("-", "", colnames(geneColored1))
  gene <- gsub("-", "", gene)
  if (facet_by != "NA") {
    if (length(unique(geneColored1[, facet_by])) == 1) {
      stop("facet_by is provided, but was not used to calculate_agg_bulk")
    }
    if (length(unique(geneColored1[, color_by])) == 1) {
      stop("color_by is provided, but was not used to calculate_agg_bulk")
    }
    if (round(sum(as.numeric(geneColored1$frac[grep(facet_bys[1], rownames(geneColored1))]))) != 100) {
      warning("The proportions reported are internal to the group_by argument used to calc_agg_bulk")
    }
  }

  g <- ggplot(geneColored1)
  if (number_labels == T) {
    ytextposnum <- -max(geneColored1[, gene])/30
    ytextposfrac <- -max(geneColored1[, gene])/10
    num = "num"
    frac = "frac"
    g <- g + geom_text(aes_string(x = color_by, y = ytextposnum, label = "num_sum"), size = 2)
    g <- g + geom_text(aes_string(x = color_by, y = ytextposfrac, label = "frac_sum"), size = 2)
  }
  if (all(!is.na(colors))) {
    g <- g + scale_color_manual(values = c(colors))
    g <- g + scale_fill_manual(values = c(colors))
  }
  if (theme == "bw") {
    g <- g + theme_bw()
  }else {
    g <- g + theme_classic()
  }

  if (title == "") {
    title <- genename
    g <- g + labs(title = title, y = genename)
  }else {
    g <- g + labs(title = title, y = genename)
  }

  g <- g + theme(plot.title = element_text(size = text_sizes[1]),
                 axis.title = element_text(size = text_sizes[2]), axis.text = element_text(size = text_sizes[3]),
                 legend.title = element_text(size = text_sizes[4]), legend.text = element_text(size = text_sizes[5]))
  g <- g + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))

  g <- g + geom_col(aes_string(x = color_by, y = "weighted_mean", fill = color_by), col = "black", alpha = alpha)
  g <- g + geom_dotplot(aes_string(x = color_by, y = gene, fill = color_by), binaxis='y', stackdir = 'center', stackratio = stackratio, dotsize = dotsize, alpha = 0.8, binwidth = (max(geneColored1[,gene])-min(geneColored1[,gene]))*binwidth)

  if (facet_by != "NA") {
    if (ncol != "NA") {
      g <- g + facet_wrap(facets = reformulate(facet_by),
                          ncol = ncol, scales = "free_x")
      g <- g + theme(strip.text.x = element_text(size = text_sizes[6]))
    }
    else {
      g <- g + facet_grid(facets = reformulate(facet_by),
                          scales = "free_x", space = "free_x")
      g <- g + theme(strip.text.x = element_text(size = text_sizes[6]))
    }
  }
  g <- g + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  return(g)
}
