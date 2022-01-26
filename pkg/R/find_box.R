### Box numbers of multivariate points #########################################

##' @title Box Numbers (Multivariate) Points Fall Into
##' @param x (n, d)-matrix of n d-dimensional points, n >= 1, d >= 1
##' @param endpoints d-list of endpoints e of the boxes, with endpoints[[j]]
##'        being the argument 'vec' of findInterval(). Defaults to the
##'        empirical quantiles of k = 5 intervals per dimension.
##' @param method character string indicating the method of how boxes
##'        are determined:
##'        - "per.dim": d box numbers per point x[i,], indicating, in each
##'          dimensions, in which box x[i,] falls.
##'        - "diagonal": one box number per point x[i,], indicating in
##'          which box (along the d-dimensional diagonal) x[i,] falls.
##'          This method requires all elements of 'endpoints' to have the
##'          same length (so the same number k of intervals per dimensions).
##'        - "lexicographic": one box number per point x[i,], indicating
##'          in which box (numbered consecutively over all *nonempty* boxes
##'          in lexicographic order) x[i,] falls.
##' @param rightmost.closed see ?findInterval (but with different default here)
##' @param left.open see ?findInterval (but with different default here)
##' @param ... additional arguments passed to the underlying findInterval()
##' @return - If method = "per.dim": (n, d)-matrix of box numbers
##'           per dimension (0 if 'in no box'); so the (i,j)th entry y[i,j] = l
##'           <=> e[[j]][l] < x[i,j] <= e[[j]][l+1] (with '<' being '<=' if l = 1)
##'           and y[i,j] = 0 if x[i,j] is in no box of e[[j]].
##'         - If method = "lexicographic" or "diagonal": n-vector
##'           with box numbers (again 0 if 'in no box').
##' @author Marius Hofert
##' @note - method != "per.dim" can be used to index a vector of colors
##'       - language: box = hyperrectangle, cube = hyperrectangle
##'         with equal sides
find_box <- function(x, endpoints = NULL,
                     method = c("per.dim", "diagonal", "lexicographic"),
                     rightmost.closed = TRUE, left.open = TRUE, ...)
{
    ## Checks
    if(!is.matrix(x))
        x <- rbind(x)
    d <- ncol(x)
    if(is.null(endpoints)) # use k = 5 intervals per dimension
        endpoints <- lapply(1:d, function(j)
            quantile(x[,j], probs = seq(0, 1, 0.2), names = FALSE))
    if(!is.list(endpoints)) endpoints <- list(endpoints) # for d = 1
    stopifnot(length(endpoints) == d,
              sapply(endpoints, function(e.j) all(diff(e.j) > 0)))

    ## For each dimension, determine in which box coordinate samples fall
    ## Note: - left.open = TRUE & rightmost.closed = TRUE => boxes are of the form
    ##         [e_{j,1}, e_{j,2}], (e_{j,2}, e_{j,3}], ..., (e_{j,k}, e_{j,k+1}]
    ##       - findInterval(1, vec = c(2, 3)) correctly returns 0
    B <- sapply(1:d, function(j) # (k, d)-matrix containing the box number in each dimension
        findInterval(x[,j], vec = endpoints[[j]],
                     rightmost.closed = rightmost.closed, left.open = left.open, ...))

    ## Switch
    method <- match.arg(method)
    switch(method,
           "per.dim" = {
               B # return B
           },
           "diagonal" = {
               if(d == 1) stop("method = \"diagonal\" requires d > 1")
               ## Require the number of boxes in all dimensions to be equal (conservative)
               elens <- sapply(endpoints, length)
               if(!all(diff(elens) == 0))
                   stop("If method = \"diagonal\" the number of endpoints have to be the same for all dimensions.")

               ## For those points on the diagonal, return the (same, first)
               ## box number (per dimension).
               if(!is.matrix(B)) B <- rbind(B) # for case n = 1
               apply(B, 1, function(b) if(all(diff(b) == 0)) b[1] else 0)
           },
           "lexicographic" = {
               ## Order boxes in lexicographic order (and keep the ordering for unsorting afterwards)
               if(!is.matrix(B)) B <- rbind(B) # for case n = 1
               ord <- do.call(order, split(B, col(B))) # (n, d)-matrix so that B[ord,] is lexicographically ordered; see https://stat.ethz.ch/pipermail/r-help/2003-July/036560.html
               B. <- B[ord,, drop = FALSE] # sorted B

               ## Iterate over all points' boxes and assign consecutive box numbers over
               ## all *nonempty*  boxes
               box <- 0 # running box number
               n <- nrow(x)
               res <- rep(NA, n)
               res[1] <- if(any(B.[1,] == 0)) 0 else 1
               if(n > 1)
                   for(i in 2:n)
                       res[i] <- if(any(B.[i,] == 0)) {
                                     0 # set to 0 as the point is in no box
                                 } else if(all(B.[i,] == B.[i-1,])) {
                                     res[i-1] # set to previous number as the two points are in the same box
                                 } else {
                                     ## Note: Increment of 1 provides consecutive box numbers for all
                                     ##       *nonempty* boxes. If at least one box is empty, the proper
                                     ##       increments among box numbers from 1 to \prod_{j=1}^d k_j
                                     ##       (with k_j = number of boxes in dimension j) are
                                     ##       not easy to determine and so the maximal number we assign
                                     ##       may be smaller than \prod_{j=1}^d k_j.
                                     box <- box + 1
                                     box
                                 }
                       stopifnot(!is.na(res)) # sanity check

               ## Undo sorting and return
               res[order(ord)] # sanity check via: stopifnot(cbind(B., res)[order(ord), 1:d] == B)
           },
           stop("Wrong 'method'"))
}
