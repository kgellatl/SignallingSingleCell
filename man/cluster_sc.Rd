% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_sc.R
\name{cluster_sc}
\alias{cluster_sc}
\title{Cluster Single Cell}
\usage{
cluster_sc(
  input,
  dimension,
  method,
  num_clust = NA,
  s = 2,
  xcol = "x",
  ycol = "y",
  set.seed = F
)
}
\arguments{
\item{input}{the input ex_sc}

\item{dimension}{either "Comp" or "2d"}

\item{method}{can either be "spectral" or "density" which is on 2d}

\item{num_clust}{the number of clusters}

\item{s}{the number of standard deviations from the curve to select cluster centers}

\item{xcol}{first column to use with dimentions for the 2d method}

\item{ycol}{second column to use with dimentions for the 2d method}
}
\description{
This will perform clustering on your single cell data.
}
\details{
This will perform clustering on either the high dimensional PCA / ICA components if dimension = Comp,
or the 2d tsne result if method = density. Typically spectral clustering works much better on higher dimensional data,
which density based clustering works better on 2d data.
}
\examples{
ex_sc_example <- cluster_sc(input = ex_sc_example, dimension = "Comp", method = "spectral", num_clust = 6)
}
