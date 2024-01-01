\name{svyreg_huber-deprecated}
\alias{svyreg_huber}
\title{Deprecated Huber Robust Survey Regression M-Estimator}
\usage{
svyreg_huber(formula, design, k, var = NULL, na.rm = FALSE, asym = FALSE,
             verbose = TRUE, ...)
}
\arguments{
    \item{formula}{a \code{[formula]} object (i.e., symbolic description
        of the model)}
    \item{design}{an object of class \code{survey.design}; see
        \code{\link[survey]{svydesign}}.}
    \item{k}{\code{[double]} robustness tuning constant
        (\eqn{0 < k \leq \infty}{0 < k <= Inf}).}
    \item{var}{a one-sided \code{[formula]} object or variable name
        (\code{[character]}) that defines the heteroscedastic variance
        or \code{[NULL]} indicating homoscedastic variance (default:
        \code{NULL}).}
    \item{na.rm}{\code{[logical]} indicating whether \code{NA} values should
        be removed before the computation proceeds (default: \code{FALSE}).}
    \item{asym}{\code{[logical]} toggle for asymmetric Huber psi-function
        (default: \code{FALSE}).}
    \item{verbose}{\code{[logical]} indicating whether additional information
            is printed to the console (default: \code{TRUE}).}
    \item{\dots}{additional arguments passed to the method (e.g., \code{maxit}:
        maxit number of iterations, etc.).}
}
\description{
    The function \code{svyreg_huber} is \bold{deprecated}; use instead
    \code{\link{svyreg_huberM}}.
}
\details{
    See \code{\link{svyreg_huberM}}.
}
\value{
    Object of class \code{svyreg.rob}
}
\seealso{
    \link{robsurvey-deprecated}
}
