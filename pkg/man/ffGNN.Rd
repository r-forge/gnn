\name{ffGNN}
\alias{ffGNN}
\alias{ffGNN.gnn_GNN}
\title{Feedforward for Generative Neural Networks}
\description{
  Feedforward method for objects of \code{\link{S3}} class \code{"gnn_GNN"}.
}
\usage{
\method{ffGNN}{gnn_GNN}(x, data)
}
\arguments{
  \item{x}{object of \code{\link{S3}} class \code{"gnn_GNN"}.}
  \item{data}{\code{\link{matrix}} to be fed forward through \code{x}.}
}
\value{
  The output \code{\link{matrix}} of \code{x} when fed with \code{data}.
}
\author{Marius Hofert}
\examples{
if(TensorFlow_available()) { # rather restrictive (due to R-Forge, winbuilder)
library(gnn) # for being standalone

## Define dummy model
d <- 2 # bivariate case
GMMN <- FNN(c(d, 300, d)) # Feedforward NN with MMD loss (a GMMN; random weights)

## Feedforward
n <- 3
set.seed(271)
X <- ffGNN(GMMN, data = matrix(runif(n * d), ncol = d))
stopifnot(dim(X) == c(n, d))

}
}
\keyword{methods}
