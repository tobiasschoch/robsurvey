\name{huber2}
\alias{huber2}
\title{Weighted Huber Proposal 2 Estimator}
\usage{
huber2(x, w, k = 1.5, na.rm = FALSE, maxit = 50, tol = 1e-04, info = FALSE,
       k_Inf = 1e6, df_cor = TRUE)
}
\arguments{
    \item{x}{\code{[numeric vector]} data.}
    \item{w}{\code{[numeric vector]} weights (same length as \code{x}).}
    \item{k}{\code{[double]} robustness tuning constant
        (\eqn{0 < k \leq \infty}{0 < k <= Inf}).}
    \item{na.rm}{\code{[logical]} indicating whether \code{NA} values should
        be removed before the computation proceeds (default: \code{FALSE}).}
    \item{maxit}{\code{[integer]} maximum number of iterations to use
        (default: \code{50}).}
    \item{tol}{\code{[double]} numerical tolerance criterion to stop the
        iterations (default: \code{1e-04}).}
    \item{info}{\code{[logical]} indicating whether additional information
        should be returned (default: \code{FALSE}).}
    \item{k_Inf}{\code{[integer]} numerical value that represents \code{Inf}
        (default: \code{1e+06}).}
    \item{df_cor}{\code{[logical]} if \code{TRUE}, the degrees of freedom of
        the estimate of scale is adjusted (default: \code{TRUE}).}
}
\description{
    Weighted Huber Proposal 2 estimator of location and scatter.
}
\details{
    Function \code{huber2} computes the weighted Huber (1964) Proposal 2
    estimates of location and scale.

    The method is initialized by the weighted median (location) and the
    weighted interquartile range (scale).
}
\section{Comparison}{
    The \code{huber2} estimator is initialized by the weighted median and
    the weighted (scaled) interquartile range. For unweighted data, this
    estimator \emph{differs} from \code{\link[MASS]{hubers}} in \pkg{MASS},
    which is initialized by \code{\link{mad}}.

    The difference between the estimators is usually negligible (for
    sufficiently small values of \code{tol}). See examples.
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
\references{
Huber, P. J. (1964). Robust Estimation of a Location Parameter.
    \emph{Annals of Mathematical Statistics} \bold{35}, 73--101.
    \doi{10.1214/aoms/1177703732}
}
\examples{
head(workplace)

# Weighted "Proposal 2" estimator of the mean
huber2(workplace$employment, workplace$weight, k = 8)

# More information on the estimate, i.e., info = TRUE
m <- huber2(workplace$employment, workplace$weight, k = 8, info = TRUE)

# Estimate of scale
m$scale

# Comparison with MASS::hubers (without weights). We make a copy of MASS::hubers
library(MASS)
hubers_mod <- hubers

# Then we replace mad by the (scaled) IQR as initial scale estimator
body(hubers_mod)[[7]][[3]][[2]] <- substitute(s0 <- IQR(y, type = 2) * 0.7413)

# Define the numerical tolerance
TOLERANCE <- 1e-8

# Comparison
m1 <- huber2(workplace$payroll, rep(1, 142), tol = TOLERANCE)
m2 <- hubers_mod(workplace$payroll, tol = TOLERANCE)$mu
m1 / m2 - 1

# The absolute relative difference is < 4.0-09 (smaller than TOLERANCE)
}
