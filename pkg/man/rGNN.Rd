\name{rGNN}
\alias{rGNN}
\alias{rGNN.gnn_GNN}
\title{Sampling from a Generative Neural Network}
\description{
  Sampling method for objects of \code{\link{S3}} class \code{"gnn_GNN"}.
}
\usage{
\method{rGNN}{gnn_GNN}(x, size, prior = NULL, pobs = FALSE, ...)
}
\arguments{
  \item{x}{object of \code{\link{S3}} class \code{"gnn_GNN"}.}
  \item{size}{sample size, a positive \code{\link{integer}}. Ignored
    if \code{prior} is a \code{\link{matrix}}.}
  \item{prior}{one of
    \describe{
      \item{\code{NULL}}{the default, generates independent N(0,1)
	realizations as prior sample.}
      \item{\code{\link{matrix}}}{passes the
	given sample through the GNN \code{x}. Such a matrix is
        internally (if \code{prior = NULL}) and typically obtained
	via \code{\link{rPrior}()}.}
  }}
  \item{pobs}{\code{\link{logical}} indicating whether the
    pseudo-observations of the generated samples should be returned.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link{rPrior}()} if \code{prior = NULL}.}
}
\value{
  \code{(size, dim(x)[1])}-\code{\link{matrix}} of samples.
}
\author{Marius Hofert}
\examples{
if(TensorFlow_available()) { # rather restrictive (due to R-Forge, winbuilder)
library(gnn) # for being standalone

## Define dummy model
d <- 2 # bivariate case
GMMN <- FNN(c(d, 300, d)) # Feedforward NN with MMD loss (a GMMN; random weights)

## Sampling
n <- 3
(X1 <- rGNN(GMMN, size = n)) # default (independent N(0,1) samples as prior)
(X2 <- rGNN(GMMN, size = n, # passing additional arguments to rPrior()
            qmargins = qexp, method = "sobol", seed = 271))
(X3 <- rGNN(GMMN, prior = matrix(rexp(n * d), ncol = d))) # providing 'prior'
stopifnot(dim(X1) == c(n, d), dim(X2) == c(n, d), dim(X3) == c(n, d))

}
}
\keyword{methods}
