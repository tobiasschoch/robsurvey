\name{weighted_mean_huber}
\alias{weighted_mean_huber}
\alias{weighted_total_huber}
\title{Weighted Huber Mean and Total (bare-bone functions)}
\usage{
weighted_mean_huber(x, w, k = 1.5, type = "rwm", asym = FALSE, info = FALSE,
    na.rm = FALSE, ...)
weighted_total_huber(x, w, k = 1.5, type = "rwm", asym = FALSE, info = FALSE,
    na.rm = FALSE, ...)
}
\arguments{
\item{x}{\code{[numeric vector]} observations.}
\item{w}{\code{[numeric vector]} weights (same length as vector \code{x}).}
\item{k}{\code{[double]} robustness tuning constant
    (\eqn{0 < k \leq \infty}{0 < k <= Inf}; default: \code{k = 1.5}).}
\item{type}{\code{[character]} type of method: \code{"rwm"} or \code{"rht"};
    see below (default: \code{"rwm"}).}
\item{asym}{\code{[logical]} toggle for asymmetric Huber psi-function
    (default: \code{FALSE}).}
\item{info}{\code{[logical]} indicating whether additional information
    should be returned (default: \code{FALSE}).}
\item{na.rm}{\code{[logical]} indicating whether \code{NA} values should
    be removed before the computation proceeds (default: \code{FALSE}).}
\item{...}{additional arguments passed to the method (e.g., \code{maxit}:
    maxit number of iterations, etc.).}
}
\value{
The return value depends on \code{info}:
\describe{
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
Weighted Huber mean and total (bare-bone functions with limited
functionality; see \code{\link{svymean_huber}} and
\code{\link{svytotal_huber}} for more capable methods)
}
\details{
\describe{
    \item{Characteristic.}{Population mean or total. Let \eqn{\mu} denote
        the estimated population mean; then, the estimated total is given
        by \eqn{N \mu} with \eqn{N =\sum w_i}, where summation is over all
        observations in the sample.}
    \item{Type.}{Two methods/types are available for estimating the location
        \eqn{\mu} (and the scale \eqn{\sigma}; see models, below):
        \describe{
            \item{\code{type = "rht"}:}{a}
            \item{\code{type = "rmw"}:}{a}
        }
    }
    \item{Variance estimation.}{See the related but more capable functions:
        \itemize{
            \item \code{\link{svymean_huber}},
            \item \code{\link{svytotal_huber}}.
        }
    }
    \item{Psi-function.}{By default, the \code{Huber} psi-function is used
        in the specification of the M-estimator. An asymmetric version of
        the Huber psi-function can be used by calling the functions with
        the additional argument \code{psi = "asymHuber"}.}
}
}
\section{Failure of convergence}{
By default, the method assumes a maximum number of \code{maxit = 100}
iterations and a numerical tolerance criterion to stop the iterations of
\code{tol = 1e-05}. You can run the code with specifications other than
the default values by specifying the arguments \code{maxit} and/or
\code{tol} in the function call; see also \code{\link{svyreg_control}}.
}

\examples{
data(workplace)

# Robust Horvitz-Thompson M-estimator of the population total
weighted_total_huber(workplace$employment, workplace$weight, k = 9,
    type = "rht")

# Robust weighted M-estimator of the population mean
weighted_mean_huber(workplace$employment, workplace$weight, k = 12,
    type = "rwm")
}
\references{
Hulliger, B. (1995). Outlier Robust Horvitz-Thompson Estimators,
\emph{Surv. Methodol.} 21, pp. 79-87.
}
\seealso{
\code{\link{weighted_mean_tukey}} and \code{\link{weighted_total_tukey}}
}
