### GNN feedforward generic ####################################################

ffGNN <- function(x, data, verbose = 0) UseMethod("ffGNN")


### GNN feedfoward method ######################################################

##' @title Feedforward Method for Objects of Class "gnn_GNN"
##' @param x object of S3 class "gnn_GNN" to be sampled from (input layer is
##'        d-dimensional)
##' @param data (n, d)-matrix of data to be fed forward through 'x'
##' @param verbose verbosity level of the underlying predict.keras.engine.training.Model()
##' @return the output (matrix) of the GNN 'x'
##' @author Marius Hofert
ffGNN.gnn_GNN <- function(x, data, verbose = 0)
{
    stopifnot(inherits(x, "gnn_GNN"))
    if(!is.matrix(data))
        data <- rbind(data)
    predict(x[["model"]], x = data, verbose = verbose)
}
