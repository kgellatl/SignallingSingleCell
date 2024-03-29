% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_heatmap.R
\name{plot_heatmap}
\alias{plot_heatmap}
\title{Plots a heatmap}
\usage{
plot_heatmap(
  input,
  genes,
  type,
  facet_by = FALSE,
  scale_group = F,
  title = "Heatmap",
  scale_by = "row",
  cluster_by = "row",
  cluster_type = "hierarchical",
  k = NULL,
  show_k = F,
  ceiling = FALSE,
  color_pal = viridis::magma(256),
  color_facets = FALSE,
  group_names = TRUE,
  gene_names = TRUE,
  text_angle = 90,
  pdf_format = "tile",
  interactive = FALSE,
  text_sizes = c(20, 10, 5, 10, 5, 5),
  gene_labels = NULL,
  gene_labels_size = 2,
  gene_labels_nudge = -0.5,
  gene_labels_col = 1,
  gene_labels_force = 1,
  return_results = F
)
}
\arguments{
\item{input}{the input ex_sc}

\item{genes}{a vector of genes to go into the heatmap}

\item{type}{can either be "bulk" or "single_cell. If type is "bulk" it will utilize the bulk info stored in fData}

\item{facet_by}{will create breaks in the heatmap by some pData Variable}

\item{title}{A heatmap title}

\item{scale_by}{scale across "row" (genes), "col" (groups), or FALSE}

\item{cluster_by}{either "row", col, or both}

\item{cluster_type}{"kmeans" or "hierarchical"}

\item{k}{if cluster type is kmeans must provide k}

\item{ceiling}{A value above which to truncate}

\item{color_pal}{The color pallete to be used}

\item{color_facets}{if true will use colors instead of text labels for the facets}

\item{group_names}{whether groups should be labelled}

\item{gene_names}{whether genes should be labelled}

\item{text_angle}{The desired angle for text on the group labels}

\item{pdf_format}{can be "tile" or "raster." tile is generally higher quality while raster is more efficient}

\item{text_sizes}{a vector of title_size, axis_title, axis_text, legend_title, legend_text, facet_text, faults too c(20,10,5,10,5,5)}

\item{plotly}{if true will be interactive
# note that this option cannot be saved with save_ggplot(), and also is time consuming for single cell heatmaps}
}
\description{
Makes either a single cell or bulk heatmap,
}
\details{
Utilize information stored in pData to control the plot display.
}
\examples{
plot_tsne_metadata(ex_sc_example, color_by = "UMI_sum", title = "UMI_sum across clusters", facet_by = "Cluster", ncol = 3)
}
