#' Cluster Single Cell
#'
#' This will perform clustering on your single cell data.
#'
#' @param input the input ex_sc
#' @param dimension either "Comp" or "2d"
#' @param method can either be "spectral" or "density" which is on 2d
#' @param num_clust the number of clusters
#' @param s the number of standard deviations from the curve to select cluster centers
#' @param xcol first column to use with dimentions for the 2d method
#' @param ycol second column to use with dimentions for the 2d method
#' @export
#' @details
#' This will perform clustering on either the high dimensional PCA / ICA components if dimension = Comp,
#' or the 2d tsne result if method = density. Typically spectral clustering works much better on higher dimensional data,
#' which density based clustering works better on 2d data.
#' @examples
#' ex_sc_example <- cluster_sc(input = ex_sc_example, dimension = "Comp", method = "spectral", num_clust = 6)

cluster_sc <- function(input,
                       dimension,
                       method,
                       num_clust = NA,
                       s=2,
                       xcol="x",
                       ycol="y") {
  if(dimension == "Comp"){
    tocluster = pData(input)[,grep("Comp", colnames(pData(input)))]
  }
  if(dimension == "2d"){
    tocluster = pData(input)[,c(xcol, ycol)]
  }
  if(method == "spectral"){
    spec <- kknn::specClust(tocluster, centers = num_clust, method = 'random-walk')
    cluster <- spec$cluster
    cluster <- paste0("Cluster", cluster)
    pData(input)$Cluster <- cluster
  }
  if(method == "density"){
    dist_tSNE = dist(tocluster)
    clust_tSNE = densityClust::densityClust(dist_tSNE, gaussian = T) # calculate density
    comb = as.data.frame(clust_tSNE$rho*clust_tSNE$delta) # combine rho and delta values
    colnames(comb) = "gamma"
    comb = comb[order(comb$gamma, decreasing = T), ,drop=F]
    comb$index = seq(nrow(comb))
    if (is.na(num_clust)) {
      # chose the max k from gamma distribution
      fit = mgcv::gam(formula = gamma ~ s(index, bs="cs"), data = log10(comb[floor(0.01*nrow(comb)):nrow(comb),]+1))
      test = log10(comb[,"index", drop=F])
      p = predict(fit, test, type = "link", se.fit = T)
      comb$pred = (10^predict(fit, test))-1
      comb$residual = comb$gamma-comb$pred
      comb$predsd = comb$pred+(s*sd(comb$residual))
      print(ggplot(comb, aes(index, gamma)) +
              geom_point(size = 0.5) +
              theme_bw() +
              scale_x_log10() +
              scale_y_log10() +
              geom_line(data=comb, aes(index, pred), colour="red") +
              geom_line(data=comb, aes(index, predsd), colour="blue"))
      cellcut = rownames(comb[comb$gamma>comb$predsd,])
    } else {
      cellcut = rownames(comb)[1:num_clust]
    }
    cellidx = which(names(clust_tSNE$rho)%in%cellcut)
    clust_tSNE = densityClust::findClusters(clust_tSNE, peaks = cellidx)
    pData(input)$Cluster = paste0("Cluster", clust_tSNE$clusters)
  }
  return(input)
}
