\name{save_load_rda}
\alias{save_gnn}
\alias{load_gnn}
\title{Save and Load .rda Files with Conversion to Raw and Keras Models}
\description{
  Save and load \code{.rda} files with conversion to objects of class
  \code{raw} (for \code{save_gnn()}) or \code{"keras.engine.training.Model"}
  (for \code{load_gnn()}).
}
\usage{
save_gnn(..., file, name = NULL)
load_gnn(file)
}
\arguments{
  \item{\dots}{objects to be saved in \code{file} (under the provided
  names if \code{name} was provided). Those objects which are of class
  \code{"gnn_GNN"} are converted with \code{\link{as.raw}()} before they
  are saved.}
  \item{file}{file name; see the underlying \code{\link{save}()}
  and \code{\link{load}()}.}
  \item{name}{\code{\link{character}} (vector) of name(s) under
  which the objects in \code{\dots} are to be saved in \code{file}. If
  \code{\link{NULL}}, the names of the objects provided by \dots are
  taken as \code{name}.}
}
\value{
  \describe{
    \item{save_gnn()}{nothing (generates an \code{.rda} by side-effect).}
    \item{load_gnn()}{the loaded object(s). Those of class \code{"gnn_GNN"}
    are converted with \code{\link{as.keras}()} before they are
    returned; this also applies to a component of a loaded object of class
    \code{\link{list}}.}
  }
}
\author{Marius Hofert}
\seealso{
  See the underlying functions \code{\link{load}()}
  and \code{\link{save}()} (among others).
}
\examples{
if(TensorFlow_available()) { # rather restrictive (due to R-Forge, winbuilder)
library(gnn) # for being standalone

file <- tempfile("foo", fileext = ".rda")
GMMN1 <- FNN()
save_gnn(GMMN1, file = file) # converts GMMN via as.raw()
GMMN2 <- load_gnn(file) # converts loaded object via as.keras()
stopifnot(is.GNN(GMMN2), inherits(GMMN2[["model"]], "keras.engine.training.Model"))
rm(GMMN1, GMMN2) # clean-up
stopifnot(file.remove(file))

}
}
\keyword{manip}
