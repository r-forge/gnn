\name{FNN}
\alias{FNN}
\title{Generative Moment Matching Network}
\description{
  Constructor for a generative feedforward neural network (FNN) model,
  an object of \code{\link{S3}} class \code{"gnn_FNN"}.
}
\usage{
FNN(dim = c(2, 2), activation = c(rep("relu", length(dim) - 2), "sigmoid"),
    batch.norm = FALSE, dropout.rate = 0, loss.fun = "MMD", n.GPU = 0, ...)
}
\arguments{
  \item{dim}{\code{\link{integer}} vector of length at least two, giving
    the dimensions of the input layer, the hidden layer(s) (if any) and
    the output layer (in this order).}
  \item{activation}{\code{\link{character}} vector of length
    \code{length(dim) - 1} specifying the activation functions
    for all hidden layers and the output layer (in this order);
    note that the input layer does not have an activation function.}
  \item{loss.fun}{loss function specified as \code{\link{character}}
    or \code{\link{function}}.}
  \item{batch.norm}{\code{\link{logical}} indicating whether batch
    normalization layers are to be added after each hidden layer.}
  \item{dropout.rate}{\code{\link{numeric}} value in [0,1] specifying
    the fraction of input to be dropped; see the rate parameter of
    \code{\link{layer_dropout}()}. Note that only if positive, dropout
    layers are added after each hidden layer.}
  \item{n.GPU}{non-negative \code{\link{integer}} specifying the number of GPUs
    available if the GPU version of TensorFlow is installed.
    If positive, a (special) multiple GPU model for data
    parallelism is instantiated. Note that for multi-layer perceptrons
    on a few GPUs, this model does not yet yield any scale-up computational
    factor (in fact, currently very slightly negative scale-ups are likely due
    to overhead costs).}
  \item{\dots}{additional arguments passed to \code{\link{loss}()}.}
}
\value{
  \code{FNN()} returns an object of \code{\link{S3}} class \code{"gnn_FNN"}
  with components
  \describe{
    \item{\code{model}}{FNN model (a \pkg{keras} object inheriting from
      the R6 classes \code{"keras.engine.training.Model"},
      \code{"keras.engine.network.Network"},
      \code{"keras.engine.base_layer.Layer"}
      and \code{"python.builtin.object"}, or a \code{\link{raw}}
      object).}
    \item{\code{type}}{\code{\link{character}} string indicating
      the type of model.}
    \item{\code{dim}}{see above.}
    \item{\code{activation}}{see above.}
    \item{\code{batch.norm}}{see above.}
    \item{\code{dropout.rate}}{see above.}
    \item{\code{n.param}}{number of trainable, non-trainable and total
      number of parameters.}
    \item{\code{loss.type}}{type of loss function (\code{\link{character}}).}
    \item{\code{n.train}}{number of training samples (\code{\link{NA_integer_}}
      unless trained).}
    \item{\code{batch.size}}{batch size (\code{\link{NA_integer_}} unless trained).}
    \item{\code{n.epoch}}{number of epochs (\code{\link{NA_integer_}}
      unless trained).}
    \item{\code{loss}}{\code{\link{numeric}(n.epoch)} containing the
      loss function values per epoch.}
    \item{\code{time}}{object of S3 class \code{"proc_time"}
      containing the training time (if trained).}
    \item{\code{prior}}{\code{\link{matrix}} containing a (sub-)sample
      of the prior (if trained).}
  }
}
\details{
  The \code{\link{S3}} class \code{"gnn_FNN"} is a subclass of the
  \code{\link{S3}} class \code{"gnn_GNN"} which in turn is a subclass of
  \code{"gnn_Model"}.
}
\author{Marius Hofert and Avinash Prasad}
\references{
  Li, Y., Swersky, K. and Zemel, R. (2015).
  Generative moment matching networks.
  \emph{Proceedings of Machine Learning Research}, \bold{37}
  (International Conference on Maching Learning), 1718--1727.
  See http://proceedings.mlr.press/v37/li15.pdf (2019-08-24)

  Dziugaite, G. K., Roy, D. M. and Ghahramani, Z. (2015).
  Training generative neural networks via maximum mean discrepancy
  optimization. \emph{AUAI Press}, 258--267.
  See http://www.auai.org/uai2015/proceedings/papers/230.pdf (2019-08-24)

  Hofert, M., Prasad, A. and Zhu, M. (2020).
  Quasi-random sampling for multivariate distributions via generative
  neural networks. \emph{Journal of Computational and Graphical
  Statistics}, \doi{10.1080/10618600.2020.1868302}.

  Hofert, M., Prasad, A. and Zhu, M. (2020).
  Multivariate time-series modeling with generative neural networks.
  See \url{https://arxiv.org/abs/2002.10645}.

  Hofert, M. Prasad, A. and Zhu, M. (2020).
  Applications of multivariate quasi-random sampling with neural
  networks. See \url{https://arxiv.org/abs/2012.08036}.
}
% Note:
% - The following runs much longer than the allowed time limit but it
%   doesn't matter as it's not run by R-Forge or CRAN checks.
% - Also note that the local check with --as-cran produces a note
%   saying that __pycache__, tmphfypfw3y.py, tmppvcp8ymi.py,
%   tmpz6_axkch.py are found. It's unclear why, these files/folders
%   don't exists locally and there is also no such warning on R-Forge
%   or winbuilder.
\examples{
if(TensorFlow_available()) { # rather restrictive (due to R-Forge, winbuilder)
library(gnn) # for being standalone

## Training data
d <- 2 # bivariate case
P <- matrix(0.9, nrow = d, ncol = d); diag(P) <- 1 # correlation matrix
ntrn <- 60000 # training data sample size
set.seed(271)
library(nvmix)
X <- abs(rNorm(ntrn, scale = P)) # componentwise absolute values of N(0,P) sample

## Plot a subsample
m <- 2000 # subsample size for plots
opar <- par(pty = "s")
plot(X[1:m,], xlab = expression(X[1]), ylab = expression(X[2])) # plot |X|
U <- apply(X, 2, rank) / (ntrn + 1) # pseudo-observations of |X|
plot(U[1:m,], xlab = expression(U[1]), ylab = expression(U[2])) # visual check

## Model 1: A basic feedforward neural network (FNN) with MSE loss function
fnn <- FNN(c(d, 300, d), loss.fun = "MSE") # define the FNN
fnn <- fitGNN(fnn, data = U, n.epoch = 40) # train with batch optimization
plot(fnn, kind = "loss") # plot the loss after each epoch

## Model 2: A GMMN (FNN with MMD loss function)
gmmn <- FNN(c(d, 300, d)) # define the GMMN (initialized with random weights)
## For training we need to use a mini-batch optimization (batch size < nrow(U)).
## For a fair comparison (same number of gradient steps) to NN, we use 500
## samples (25% = 4 gradient steps/epoch) for 10 epochs for GMMN.
library(keras) # for callback_early_stopping()
## We monitor the loss function and stop earlier if the loss function
## over the last patience-many epochs has changed by less than min_delta
## in absolute value. Then we keep the weights that led to the smallest
## loss seen throughout training.
gmmn <- fitGNN(gmmn, data = U, batch.size = 500, n.epoch = 10,
               callbacks = callback_early_stopping(monitor = "loss",
                                                   min_delta = 1e-3, patience = 3,
                                                   restore_best_weights = TRUE))
plot(gmmn, kind = "loss") # plot the loss after each epoch
## Note:
## - Obviously, in a real-world application, batch.size and n.epoch
##   should be (much) larger (e.g., batch.size = 5000, n.epoch = 300).
## - Training is not reproducible (due to keras).

## Model 3: A FNN with CvM loss function
fnnCvM <- FNN(c(d, 300, d), loss.fun = "CvM")
fnnCvM <- fitGNN(fnnCvM, data = U, batch.size = 500, n.epoch = 10,
                 callbacks = callback_early_stopping(monitor = "loss",
                                                     min_delta = 1e-3, patience = 3,
                                                     restore_best_weights = TRUE))
plot(fnnCvM, kind = "loss") # plot the loss after each epoch

## Sample from the different models
set.seed(271)
V.fnn <- rGNN(fnn, size = m)
set.seed(271)
V.gmmn <- rGNN(gmmn, size = m)
set.seed(271)
V.fnnCvM <- rGNN(fnnCvM, size = m)

## Joint plot of training subsample with GMMN PRNs. Clearly, the MSE
## cannot be used to learn the distribution correctly.
layout(matrix(1:4, ncol = 2, byrow = TRUE))
plot(U[1:m,], xlab = expression(U[1]), ylab = expression(U[2]), cex = 0.2)
mtext("Training subsample", side = 4, line = 0.4, adj = 0)
plot(V.fnn,    xlab = expression(V[1]), ylab = expression(V[2]), cex = 0.2)
mtext("Trained NN with MSE loss", side = 4, line = 0.4, adj = 0)
plot(V.gmmn,  xlab = expression(V[1]), ylab = expression(V[2]), cex = 0.2)
mtext("Trained NN with MMD loss", side = 4, line = 0.4, adj = 0)
plot(V.fnnCvM,    xlab = expression(V[1]), ylab = expression(V[2]), cex = 0.2)
mtext("Trained NN with CvM loss", side = 4, line = 0.4, adj = 0)

## Joint plot of training subsample with GMMN QRNs
library(qrng) # for sobol()
V.fnn.    <- rGNN(fnn,    size = m, method = "sobol", randomize = "Owen")
V.gmmn.   <- rGNN(gmmn,   size = m, method = "sobol", randomize = "Owen")
V.fnnCvM. <- rGNN(fnnCvM, size = m, method = "sobol", randomize = "Owen")
plot(U[1:m,], xlab = expression(U[1]), ylab = expression(U[2]), cex = 0.2)
mtext("Training subsample", side = 4, line = 0.4, adj = 0)
plot(V.fnn., xlab = expression(V[1]), ylab = expression(V[2]), cex = 0.2)
mtext("Trained NN with MSE loss", side = 4, line = 0.4, adj = 0)
plot(V.gmmn., xlab = expression(V[1]), ylab = expression(V[2]), cex = 0.2)
mtext("Trained NN with MMD loss", side = 4, line = 0.4, adj = 0)
plot(V.fnnCvM., xlab = expression(V[1]), ylab = expression(V[2]), cex = 0.2)
mtext("Trained NN with CvM loss", side = 4, line = 0.4, adj = 0)
layout(1)
par(opar)
}
}
\keyword{models}
