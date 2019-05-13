### Setup

setMethod("show",
          signature = "Ex_Sc",
          definition = function(object){
            cat("An object of class ", class(object), "\n", sep = "")
            cat("", nrow(object@counts), "features by ",
                ncol(object@counts), " cells")
            invisible(NULL)
          })

setMethod("show",
          signature = "Signal_Set",
          definition = function(object){
            dims <- dim(object@assays@counts)
            cat("An object of class ", class(object), "\n", sep = "")
            cat("", dims[1], "features by ", dims[2], "cells")
            invisible(NULL)
          })


exsc
signal_set

dim(signal_set@assays@counts)[[1]]
