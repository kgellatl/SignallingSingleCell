#' This will create a bar plot
#'
#' This will plot a given gene via violin plot
#'
#' @param input The input data
#' @param title The title
#' @param gene if provided will color_by the gene
#' @param color_by a pData variable
#' @param colors What colors to utilize for categorial data. Be sure it is of the proper length!
#' @param facet_by a pData variable
#' @param ncol How many columns if faceting
#' @param text_sizes a vector of title_size, axis_title, axis_text, legend_title, legend_text, facet_text, faults too c(20,10,5,10,5,5)
#' @param theme the plot theme
#' @export
#' @details
#' Utilize information stored in pData to control the plot display.
#' @examples
#' plot_violin(input, title = "Actb across clusters", gene = "Actb", color_by = "Timepoint", facet_by = "Cluster", size = 1, ncol = 3)

plot_barplot <- function(input, title = "", gene, color_by, facet_by = "NA", colors = NA, ncol = "NA", number_labels = T,
                         text_sizes = c(20,10,5,10,5,5), theme = "classic", alpha = 0.5){

  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }

  cols <- gg_color_hue(length(unique(pData(input)[,color_by])))
  geneColored1 <- fData(input)

  ### Check that aggregate bulk was calculated
  ind <- grep("bulk", colnames(fData(input)))
  if(length(ind) == 0){
    stop("Must calc_agg_bulk() before using this function")
  }
  geneColored1 <- geneColored1[,ind]
  color_bys <- sort(unique(pData(input)[,color_by]))


  if(facet_by == "NA"){
    ind <- grep(color_bys[1], colnames(geneColored1))
    if(length(ind) == 0){
      stop("Provided color_by argument was not used to calc_agg_bulk()")
    }
  }


  if(facet_by == "NA"){
    geneColored1 <- as.data.frame(t(geneColored1[gene,]))
    geneColored1[,color_by] <- 0
    for (i in 1:length(color_bys)) {
      int <- color_bys[i]
      ind <- grep(int, rownames(geneColored1))
      geneColored1[ind,color_by] <- color_bys[i]
    }
  } else {
    geneColored1 <- as.data.frame(t(geneColored1[gene,]))
    geneColored1[,color_by] <- 0
    geneColored1[,facet_by] <- 0
    facet_bys <- sort(unique(pData(input)[,facet_by]))
    for (i in 1:length(color_bys)) {
      int <- color_bys[i]
      ind <- grep(int, rownames(geneColored1))
      geneColored1[ind,color_by] <- int

    }
    for (i in 1:length(facet_bys)) {
      int <- facet_bys[i]
      ind <- grep(int, rownames(geneColored1))
      geneColored1[ind,facet_by] <- int
    }
  }

  ### Grab number of cells in group
  for (i in 1:nrow(geneColored1)) {
    int_string <- rownames(geneColored1)[i]
    int_string <- strsplit(int_string, split = "_")[[1]]
    geneColored1$num[i] <- int_string[match("cells", int_string)+1]
  }

  for (i in 1:nrow(geneColored1)) {
    int_string <- rownames(geneColored1)[i]
    int_string <- strsplit(int_string, split = "_")[[1]]
    geneColored1$frac[i] <- int_string[match("percent", int_string)+1]
  }



  colnames(geneColored1) <- gsub("-", "", colnames(geneColored1))
  gene <- gsub("-", "", gene)

  if(facet_by == "NA"){
    if(length(unique(geneColored1[,color_by])) < nrow(geneColored1)){
      stop("calc_agg_bulk was provided multiple aggregate_by values, you are only providing color_by. Please provide a facet_by argument")
    }
  } else {
    if(length(unique(geneColored1[,facet_by])) == 1){
      stop("facet_by is provided, but was not used to calculate_agg_bulk")
    }
    if(length(unique(geneColored1[,color_by])) == 1){
      stop("facet_by is provided, but was not used to calculate_agg_bulk")
    }
  }




  ###
  g <- ggplot(geneColored1)
  if(number_labels == T){
    ytextposnum <- -max(geneColored1[,gene])/30
    ytextposfrac <- -max(geneColored1[,gene])/10
    num = "num"
    frac = "frac"
    g <- g + geom_text(aes_string(x = color_by, y = ytextposnum, label = num))
    g <- g + geom_text(aes_string(x = color_by, y = ytextposfrac, label = frac))

  }
  if(all(!is.na(colors))){
    g <- g + scale_color_manual(values = c(colors))
  }
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
  g <- g + theme(plot.title = element_text(size = text_sizes[1]), axis.title = element_text(size = text_sizes[2]), axis.text = element_text(size = text_sizes[3]), legend.title = element_text(size = text_sizes[4]), legend.text=element_text(size=text_sizes[5]))
  g <- g + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  g <- g + geom_col(aes_string(x=color_by, y=gene, fill = color_by), col = "black", alpha = alpha)
  if(facet_by != "NA"){
    if(ncol != "NA"){
      g <- g +  facet_wrap(facets = reformulate(facet_by), ncol = ncol, scales = "free_x")
      g <- g + theme(strip.text.x = element_text(size = text_sizes[6]))
    } else {
      g <- g +  facet_grid(facets = reformulate(facet_by), scales = "free_x", space = "free_x")
      g <- g + theme(strip.text.x = element_text(size = text_sizes[6]))
    }
  }
  #####
  g <- g + theme(axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank())
  return(g)
}
