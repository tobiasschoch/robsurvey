\name{robsvyreg}
\alias{robsvyreg}
\alias{svyreg_control}
\title{Internal Function for the Regression GM-Estimator}
\usage{
robsvyreg(x, y, w, k, psi, type, xwgt, var = NULL, verbose = TRUE, ...)
svyreg_control(tol = 1e-5, maxit = 100, k_Inf = 1e6, init = NULL,
               mad_center = TRUE, ...)
}
\arguments{
    \item{x}{\code{[numeric matrix]} design matrix (\code{NA} values not
        allowed).}
    \item{y}{\code{[numeric vector]} dependent variable (\code{NA} values not
        allowed).}
    \item{w}{\code{[numeric vector]} weights (no \code{NA}'s allowed).}
    \item{k}{\code{[double]} robustness tuning constant
        (\eqn{0 < k \leq \infty}{0 < k <= Inf}).}
    \item{psi}{\code{[integer]} psi-functions: \code{0}: Huber, \code{1}:
        asymmetric Huber, and \code{2}: Tukey biweight.}
    \item{type}{\code{[integer]} type of estimator; \code{0}: M-estimator;
        \code{1}: Mallows and \code{2}: Schweppe type GM-estimator.}
    \item{xwgt}{\code{[numeric vector]} weights for design space used in
        GM-estimators (default: \code{NULL}, \code{NA} values not allowed).}
    \item{var}{\code{[numeric vector]} heteroscedastic variance
        (default: \code{NULL}).}
    \item{verbose}{\code{[logical]} indicating whether additional information
        is printed to the console (default: \code{TRUE}).}
    \item{\dots}{additional arguments passed to the method
        (see \code{svyreg_control}).}
    \item{tol}{\code{[double]} numerical tolerance criterion to stop the
        iterations (default: \code{1e-05}).}
    \item{maxit}{\code{[integer]} maximum number of iterations to use
        (default: \code{100}).}
    \item{k_Inf}{\code{[integer]} numerical value that represents \code{Inf}
        (default: \code{1e+06}).}
    \item{init}{either \code{NULL} or \code{[numeric vector]}, if
        \code{init = NULL} the regression estimator is initialized by
        weighted least squares; otherwise, \code{init} can be specified as
        the estimate (i.e., \emph{p}-vector) to initialize the iteratively
        re-weighted least squares method (default: \code{NULL}).}
    \item{mad_center}{\code{[logical]} if \code{TRUE}, the weighted MAD is
        centered about the (weighted) median, otherwise the weighted MAD is
        centered about zero (default: \code{TRUE}).}
}
\description{
    \strong{Internal} function to call the robust survey regression
    \emph{GM}-estimator; this function is \strong{only} intended for internal
    use.  The function does \strong{not} check or validate the arguments.
    In particular, missing values in the data may make the function crash.
}
\details{
    Not documented
}
\value{
\code{[list]}
}
