### Setup

setMethod("show",
          signature = "SignalSet",
          definition = function(object){
            dims <- dim(object@counts)

            cat("An object of class ", class(object), "\n", sep = "")

            cat("", dims[1], "features by ", dims[2], "cells \n")

            cat("", "This object contains ", length(object@group_ids), "groups\n")

            cat("", "The active group is ", object@active_group@group, "\n")


            # if(length(object@active_group@DE) == 0){
            #   cat("  DE has not been performed on the active group \n")
            # } else {
            #   cat("  DE COMPLETE on active group\n")
            #
            # }
            #
            # if(length(object@active_group@AggBulk) == 0){
            #   cat("  DE has not been performed on the active group \n")
            # } else {
            #   cat("  DE COMPLETE on active group\n")
            #
            # }
            #
            # if(length(object@active_group@DE) == 0){
            #   cat("  DE has not been performed on the active group \n")
            # } else {
            #   cat("  DE COMPLETE on active group\n")
            #
            # }

            invisible(NULL)
          })



setGeneric("pData", function(object, ...) standardGeneric("pData"))
setMethod("pData", "SignalSet", function(object) object@active_group@pD)
setGeneric("pData<-", function(object, value) standardGeneric("pData<-"))
setMethod("pData<-", "SignalSet", function(object, value){
  object@active_group@pD <- value
  if(validObject(object))
    return(object)
})


setGeneric("fData", function(object, ...) standardGeneric("fData"))
setMethod("fData", "SignalSet", function(object) object@active_group@fD)
setGeneric("fData<-", function(object, value) standardGeneric("fData<-"))
setMethod("fData<-", "SignalSet", function(object, value){
  object@active_group@fD <- value
  if(validObject(object))
    return(object)
})


# setMethod("[", "SignalSet",
#           function(x,i,j,drop="missing") {
#             .marray <- x@assays@counts[i, j]
#             .mnormarray <- x@assays@norm_counts[i, j]
#             .pmeta <- x@assays@pD[j, ]
#             .fmeta <- x@assays@fD[i, ]
#             SignalSet(
#               assays = ExSc(
#                 counts = .marray,
#                 norm_counts = .mnormarray,
#                 fD = .fmeta,
#                 pD = .pmeta),
#               AggBulk = x@AggBulk,
#               DE = x@DE,
#               network_dataframe = x@network_dataframe,
#               network_igraph = x@network_igraph,
#               commands = x@commands)
#           })



