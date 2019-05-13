### Setup

setMethod("show",
          signature = "ExSc",
          definition = function(object){

            cat("An object of class ", class(object), "\n", sep = "")

            cat("", nrow(object@counts), "features by ",
                ncol(object@counts), " cells")

            invisible(NULL)
          })

setMethod("show",
          signature = "SignalSet",
          definition = function(object){
            dims <- dim(object@assays@counts)

            cat("An object of class ", class(object), "\n", sep = "")

            cat("", dims[1], "features by ", dims[2], "cells \n")

            if(is.null(dim(object@DimReduc))){
              cat("", "dim_reduce has not been performed \n")
            } else {
              cat("", "Dimension Reduction Results present \n")
            }

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
setMethod("pData", "SignalSet", function(object) object@assays@pD)

setGeneric("fData", function(object, ...) standardGeneric("fData"))
setMethod("fData", "SignalSet", function(object) object@assays@fD)

