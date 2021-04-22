\name{robsvyreg}
\alias{robsvyreg}
\alias{svyreg_control}
\title{Internal Function for the Regression GM-Estimator}
\usage{
robsvyreg(x, y, w, k, psi, type, xwgt = NULL, var = NULL, ...)
svyreg_control(tol = 1e-05, maxit = 100, k_Inf = 1e+05, ...)
}
\arguments{
\item{x}{\code{[array]} design matrix (no \code{NA}'s allowed).}
\item{y}{\code{[numeric vector]} dependent variable (no \code{NA}'s allowed).}
\item{w}{\code{[numeric vector]} weights (no \code{NA}'s allowed).}
\item{k}{\code{[double]} robustness tuning constant
    (\eqn{0 < k \leq \infty}{0 < k <= Inf}).}
\item{psi}{\code{[integer]} psi-functions: \code{0}: Huber, \code{1}:
    asymmetric Huber, \code{2}: Tukey biweight.}
\item{type}{\code{[integer]} type of estimator; \code{0}: M-estimator;
    \code{1}: Mallows and \code{2}: Schweppe type GM-estimator.}
\item{xwgt}{\code{[numeric vector]} weights for design space used in
    GM-estimators (default: \code{NULL}, (no \code{NA}'s allowed).}
\item{var}{\code{[numeric vector]} heteroscedastic variance
    (default: \code{NULL}).}
\item{...}{additional arguments passed to the method
    (see \code{svyreg_control}).}
\item{tol}{\code{[double]} numerical tolerance criterion to stop the
    iterations (default: \code{1e-05}).}
\item{maxit}{\code{[integer]} maximum number of iterations to use
    (default: \code{100}).}
\item{k_Inf}{\code{[integer]} numerical value that represents \code{Inf}
    (default: \code{1e+05}).}
}
\value{
\code{[list]}
}
\description{
\strong{Internal} function for robust survey regression GM-estimator;
this function is \strong{only} intended for internal use. The function
does \strong{not} check or validate the arguments. In particular,
missing values in the data may make the function crash.
}
\details{
Not documented
}
