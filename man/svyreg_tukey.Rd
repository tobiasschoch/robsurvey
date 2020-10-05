\name{svyreg_tukey}
\alias{svyreg_tukey}
\alias{svyreg_tukeyGM}
\title{Huber robust survey regression M- and GM-estimator}
\usage{
svyreg_tukey(formula, design, k, var = NULL, na.rm = FALSE, ...)

svyreg_tukeyGM(formula, design, k, type, xwgt, var = NULL, na.rm = FALSE, ...)
}
\arguments{
\item{formula}{a \code{[formula]} object (i.e., symbolic description of the model)}

\item{design}{an object of class \code{survey.design} or \code{survey.design2}.}

\item{k}{\code{[double]} robustness tuning constant (\eqn{0 < k \leq \infty}{0 < k <= Inf}).}

\item{var}{\code{[numeric vector]} heteroscedastic variance (default: \code{NULL}, i.e., homoscedastic variance).}

\item{na.rm}{\code{[logical]} indicating whether \code{NA} values should be removed before the computation proceeds (default: \code{FALSE}).}

\item{...}{additional arguments passed to the method (e.g., \code{maxit}: maxit number of iterations, etc.).}

\item{type}{\code{[character]} \code{"Mallows"} or \code{"Schweppe"}.}

\item{xwgt}{\code{[numerical vector]} of weights in the design space (default: \code{NULL}); \code{xwgt} is only relevant if \code{type = "Mallows"} or \code{type = "Schweppe"}.}
}
\value{
object of class \code{svyreg.rob}
}
\description{
Huber robust survey regression M- and GM-estimator (Mallows and Schweppe type)
}
\details{
Details
}
\section{Failure of convergence}{
 
By default, the method assumes a maximum number of \code{maxit = 100} iterations and a numerical tolerance criterion to stop the iterations of \code{tol = 1e-05}. You can run the code with specifications other than the default values by specifying the arguments \code{maxit} and/or \code{tol} in the function call; see also \code{\link{svyreg_control}}.
}

