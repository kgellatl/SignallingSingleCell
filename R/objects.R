### Setup

ExSc <- setClass(Class = "ExSc",

                  slots = c(
                    counts = "dgCMatrix",
                    norm_counts = "dgCMatrix",
                    pD = 'data.frame',
                    fD = 'data.frame'),

                  prototype = list(counts = Matrix::Matrix(c(0,1), nrow = 100, ncol = 10, sparse = T),
                                   norm_counts = Matrix::Matrix(c(0,1), nrow = 100, ncol = 10, sparse = T),
                                   pD = as.data.frame(matrix(NA, nrow = 10, ncol =2)),
                                   fD = as.data.frame(matrix(NA, nrow = 100, ncol =2))),

                  validity = function(object){
                    errors <- character()
                    genes <- nrow(object@counts)
                    cells <- ncol(object@counts)
                    if(genes != nrow(object@fD)){
                      msg <- c("fD must contain as many rows as the input matrix has rows")
                      errors <- c(errors,msg)
                    }
                    cells <- ncol(object@counts)
                    if(cells != nrow(object@pD)){
                      msg <- ("pD must contain as many rows as the input matrix has cols")
                      errors <- c(errors,msg)
                    }
                    if (length(errors) == 0 ) TRUE else errors
                  }
)

DimReduc <- setClass(Class = "DimReduc",

                  slots = c(
                    cells = "matrix",
                    x_y = "matrix")
)

SignalSet  <- setClass(Class = "SignalSet",

                        slots = c(
                          assays = "ExSc",
                          DimReduc = "DimReduc",
                          AggBulk = "matrix",
                          DE = "list",
                          network_dataframe = "data.frame",
                          network_igraph = "list")
)




