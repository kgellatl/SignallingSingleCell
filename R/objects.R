### Setup

ExSc <- setClass(Class = "ExSc",

                  slots = c(
                    counts = "dgCMatrix",
                    norm_counts = "dgCMatrix",
                    pD = 'data.frame',
                    fD = 'data.frame'),

                  prototype = list(counts = Matrix::Matrix(c(0,1), nrow = 100, ncol = 50, sparse = T),
                                   norm_counts = Matrix::Matrix(c(0,1), nrow = 100, ncol = 50, sparse = T),
                                   pD = as.data.frame(matrix(NA, nrow = 50, ncol =2)),
                                   fD = as.data.frame(matrix(NA, nrow = 100, ncol =2)))
)

SignalSetCommand <- setClass(
  Class = 'SignalSetCommand',
  slots = c(
    name = 'character',
    time.stamp = 'POSIXct',
    call.string = 'character',
    params = 'ANY'
  )
)

SignalSet  <- setClass(Class = "SignalSet",

                        slots = c(
                          assays = "ExSc",
                          AggBulk = "matrix",
                          DE = "list",
                          network_dataframe = "data.frame",
                          network_igraph = "list",
                          commands = "list"),

                       validity = function(object){
                         errors <- character()
                         genes <- nrow(object@assays@counts)
                         cells <- ncol(object@assays@counts)
                         if(genes != nrow(object@assays@fD)){
                           msg <- c("fD must contain as many rows as the input matrix has rows")
                           errors <- c(errors,msg)
                         }
                         cells <- ncol(object@assays@counts)
                         if(cells != nrow(object@assays@pD)){
                           msg <- ("pD must contain as many rows as the input matrix has cols")
                           errors <- c(errors,msg)
                         }
                         if (length(errors) == 0 ) TRUE else errors
                       }
)




