\name{find_box}
\alias{find_box}
\title{Box Numbers (Multivariate) Points Fall Into}
\description{
  Finding the numbers of boxes that given (multivariate) points fall
  into (the default is similar to \code{\link{findInterval}()} but
  other methods are provided, too).
}
\usage{
find_box(x, endpoints = NULL,
         method = c("per.dim", "diagonal", "lexicographic"),
         rightmost.closed = TRUE, left.open = TRUE, ...)
}
\arguments{
  \item{x}{\eqn{(n, d)}-matrix containing \eqn{n}
    \eqn{d}-dimensional data points, \eqn{n \ge 1}{n >= 1},
    \eqn{d \ge 1}{d >= 1}.}
  \item{endpoints}{\eqn{d}-\code{\link{list}} containing
    numeric vectors of endpoints of the intervals in each dimension
    (each of the \eqn{d} elements is an argument \code{vec}
    as required by \code{\link{findInterval}()}).}
  \item{method}{\code{\link{character}} string indicating the
    method to be used. Available are:
    \describe{
      \item{\code{"per.dim"}}{the default. Each row \code{x[i,]}
	of \code{x} produces \eqn{d} numbers, where the \eqn{j}th
	indicates in which interval \code{x[i,j]} falls. This is
	essentially \code{\link{findInterval}()} applied to the \eqn{d}
	coordinate samples (the columns of \code{x}).
      }
      \item{\code{"diagonal"}}{Each row \code{x[i,]} produces one
	number, indicating in which box \code{x[i,]} falls along
	the \eqn{d}-dimensional diagonal (note that this method
        requires all elements of \code{endpoints} to have the
	same length, so that the diagonal is well-defined).}
      \item{\code{"lexicographic"}}{Each row \code{x[i,]} produces one
	number, indicating in which box \code{x[i,]} falls
	if all \emph{nonempty} boxes are numbered consecutively
	in lexicographic order.}
    }
  }
  \item{rightmost.closed}{see \code{\link{findInterval}()} (note
    the different default here).}
  \item{left.open}{see \code{\link{findInterval}()} (note
    the different default here).}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link{findInterval}()}.}
}
\value{
  \describe{
    \item{\code{"per.dim"}}{\eqn{(n,d)}-matrix of box numbers per
      dimension.}
    \item{\code{"diagonal"}, \code{"lexicographic"}}{\eqn{n}-vector
      with box numbers.}
  }
  Note that, as \code{\link{findInterval}()}, \eqn{0} means `in no box'.
}
\details{
  The box numbers can be used, for example, to color points; see the
  examples below.
}
\author{
  Marius Hofert
}
\examples{
### 1 Basic example calls ######################################################

## Example data
n <- 1000
d <- 2
set.seed(271)
U <- matrix(runif(n * d), ncol = d) # (n, d)-matrix of data (here: in [0,1]^d)

## Define endpoints and evaluate for different methods
epts <- lapply(1:d, function(j) seq(0, 1, by = 1/5)) # 5 boxes per dimension
find_box(U, endpoints = epts)[1:10,] # default "per.dim" (first 10 points only)
boxes.diag <- find_box(U, endpoints = epts, method = "diagonal")
boxes.lexi <- find_box(U, endpoints = epts, method = "lexicographic")

## Special cases
## First row of U (n = 1)
U[1,] # ~= (0.25, 0.14)
stopifnot(find_box(U[1, 1:2], endpoints = epts) == c(2, 1))
stopifnot(find_box(U[1, 1:2], endpoints = epts, method = "diagonal") == 0)
stopifnot(find_box(U[1, 1:2], endpoints = epts, method = "lexicographic") == 1)
## Note concerning the last line: It's 1 because all other boxes are empty
## Single number U[1,1] (d = 1)
U[1,1] # ~= 0.25
stopifnot(find_box(U[1,1], endpoints = epts[[1]]) == 2)
stopifnot(find_box(U[1,1], endpoints = epts[[1]], method = "lexicographic") == 1)


### 2 Coloring points along the diagonal #######################################

## Colors
library(RColorBrewer)
basecols <- c("#000000", brewer.pal(8, name = "Dark2")[c(7, 3, 5, 4, 6)])
mypal <- colorRampPalette(basecols, space = "Lab") # palette for interpolation
ncols <- diff(range(boxes.diag)) + 1 # number of colors we need
palette(mypal(ncols)) # set palette accordingly

## Boxes of equal size
plot(U, pch = 20, xlab = expression(U[1]), ylab = expression(U[2]),
     col = if(min(boxes.diag) == 0) boxes.diag + 1 else boxes.diag)
abline(v = epts[[1]], h = epts[[2]], col = "gray50") # guides

## Boxes of different sizes (have to be the same number of
## boxes per dimension, otherwise there is no obvious 'diagonal')
epts. <- lapply(1:d, function(j) c(0, 0.1, 0.3, 0.6, 1)) # 4 boxes per dimension
boxes.diag <- find_box(U, endpoints = epts., method = "diagonal")
plot(U, pch = 20, xlab = expression(U[1]), ylab = expression(U[2]),
     col = if(min(boxes.diag) == 0) boxes.diag + 1 else boxes.diag)
abline(v = epts.[[1]], h = epts.[[2]], col = "gray50") # guides


### 3 Coloring points in lexicographic ordering ################################

## Boxes of equal size
ncols <- diff(range(boxes.lexi)) + 1 # number of colors we need
palette(mypal(ncols)) # set palette accordingly
plot(U, pch = 20, xlab = expression(U[1]), ylab = expression(U[2]),
     col = if(min(boxes.lexi) == 0) boxes.lexi + 1 else boxes.lexi)
abline(v = epts[[1]], h = epts[[2]], col = "gray50") # guides

## Boxes of different sizes and numbers
epts. <- list(seq(0.2, 1, by = 1/5), seq(1/3, 1, by = 1/3))
boxes.lexi <- find_box(U, endpoints = epts., method = "lexicographic")
ncols <- diff(range(boxes.lexi)) + 1
palette(mypal(ncols))
plot(U, pch = 20, xlab = expression(U[1]), ylab = expression(U[2]),
     col = if(min(boxes.lexi) == 0) boxes.lexi + 1 else boxes.lexi)
abline(v = epts.[[1]], h = epts.[[2]], col = "gray50")
}
\keyword{programming}
