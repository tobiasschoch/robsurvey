\name{huber2}
\alias{huber2}
\title{Weighted Huber Proposal 2 Estimator}
\usage{
huber2(x, w, k = 1.5, na.rm = FALSE, maxit = 50, tol = 1e-04, info = FALSE,
    k_Inf = 1e5, df_cor = TRUE)
}
\arguments{
\item{x}{\code{[numeric vector]} observations.}
\item{w}{\code{[numeric vector]} weights (same length as vector \code{x}).}
\item{k}{\code{[double]} robustness tuning constant
    (\eqn{0 < k \leq \infty}{0 < k <= Inf}).}
\item{na.rm}{\code{[logical]} indicating whether \code{NA} values should
    be removed before the computation proceeds (default: \code{FALSE}).}
\item{maxit}{\code{[integer]} maximum number of iterations to use
    (default: \code{50}).}
\item{tol}{\code{[double]} numerical tolerance criterion to stop the
    iterations (default: \code{1e-04}).}
\item{info}{\code{[logical]} observations (default: \code{FALSE}).}
\item{k_Inf}{\code{[integer]} numerical value that represents \code{Inf}
    (default: \code{1e+05}).}
\item{df_cor}{\code{[logical]} toggle for the adjustment of the degrees of
    freedom for the estimate of scale (default: \code{TRUE}).}
}
\value{
The return value depends on \code{info}: \describe{
    \item{\code{info = FALSE}:}{estimate of mean or total \code{[double]}}
    \item{\code{info = TRUE}:}{a \code{[list]} with items:
        \itemize{
            \item \code{characteristic} \code{[character]},
            \item \code{estimator} \code{[character]},
            \item \code{estimate} \code{[double]},
            \item \code{variance} (default: \code{NA}),
            \item \code{robust} \code{[list]},
            \item \code{residuals} \code{[numeric vector]},
            \item \code{model} \code{[list]},
            \item \code{design} (default: \code{NA}),
            \item \code{[call]}
        }
    }
}
}
\description{
Weighted Huber Proposal 2 estimator
}
\details{
The function \code{huber2} computes the weighted Huber (1964) Proposal 2
estimates of location and scale.

The method is initialized by the weighted median (location) and the
weighted interquartile range (scale).
}
\examples{
data(workplace)

# Weighted Proposal 2 M-estimator of the mean
huber2(workplace$employment, workplace$weight, k = 8)
}
\references{
Huber, P. J. (1964). Robust Estimation of a Location Parameter,
\emph{Ann. Math. Statist.} 35, pp. 73--101.
}
