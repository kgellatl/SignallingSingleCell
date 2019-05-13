### Setup

setMethod("show",
          signature = "SignalSet",
          definition = function(object){
            dims <- dim(object@assays@counts)

            cat("An object of class ", class(object), "\n", sep = "")

            cat("", dims[1], "features by ", dims[2], "cells \n")

            if(is.null(dim(object@DE))){
              cat("", "DE has not been performed \n")
            } else {
              cat("", "DE Results present \n")
            }

            if(nrow(object@AggBulk) == 0){
              cat("", "calc_agg_bulk has not been performed \n")
            } else {
              cat("", "AggBulk Results present \n")
            }

            if(nrow(object@network_dataframe) == 0){
              cat("", "calc_rl_connections has not been performed \n")
            } else {
              cat("", "network_dataframe Results present \n")
            }

            if(length(object@network_igraph) == 0){
              cat("", "build_rl_network has not been performed \n")
            } else {
              cat("", "network_igraph Results present \n")
            }

            invisible(NULL)
          })

setGeneric("pData", function(object, ...) standardGeneric("pData"))
setGeneric("fData", function(object, ...) standardGeneric("fData"))
setMethod("pData", "SignalSet", function(object) object@assays@pD)
setMethod("fData", "SignalSet", function(object) object@assays@fD)

setGeneric("pData<-", function(object, value) standardGeneric("pData<-"))
setGeneric("fData<-", function(object, value) standardGeneric("fData<-"))
setMethod("pData<-", "SignalSet", function(object, value){
  object@assays@pD <- value
  if(validObject(object))
    return(object)
})
setMethod("fData<-", "SignalSet", function(object, value){
  object@assays@fD <- value
  if(validObject(object))
    return(object)
})


setMethod("[", "SignalSet",
          function(x,i,j,drop="missing") {
            .marray <- x@assays@counts[i, j]
            .mnormarray <- x@assays@norm_counts[i, j]
            .pmeta <- x@assays@pD[j, ]
            .fmeta <- x@assays@fD[i, ]
            SignalSet(
              assays = ExSc(
                counts = .marray,
                norm_counts = .mnormarray,
                fD = .fmeta,
                pD = .pmeta),
              AggBulk = x@AggBulk,
              DE = x@DE,
              network_dataframe = x@network_dataframe,
              network_igraph = x@network_igraph,
              commands = x@commands)
          })



