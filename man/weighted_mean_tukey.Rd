\name{weighted_mean_tukey}
\alias{weighted_mean_tukey}
\alias{weighted_total_tukey}
\title{Weighted Tukey Biweight Mean and Total (bare-bone functions)}
\usage{
weighted_mean_tukey(x, w, k, type = "rwm", info = FALSE, na.rm = FALSE,
    verbose = TRUE, ...)
weighted_total_tukey(x, w, k, type = "rwm", info = FALSE, na.rm = FALSE,
    verbose = TRUE, ...)
}
\arguments{
\item{x}{\code{[numeric vector]} data.}
\item{w}{\code{[numeric vector]} weights (same length as vector \code{x}).}
\item{k}{\code{[double]} robustness tuning constant
    (\eqn{0 < k \leq \infty}{0 < k <= Inf}).}
\item{type}{\code{[character]} type of method: \code{"rwm"} or \code{"rht"};
    see below (default: \code{"rwm"}).}
\item{info}{\code{[logical]} indicating whether additional information
    should be returned (default: \code{FALSE}).}
\item{na.rm}{\code{[logical]} indicating whether \code{NA} values should
    be removed before the computation proceeds (default: \code{FALSE}).}
\item{verbose}{\code{[logical]} indicating whether additional information
		is printed to the console (default: \code{TRUE}).}
\item{\dots}{additional arguments passed to the method (e.g., \code{maxit}:
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
Weighted Tukey biweight M-estimator of the mean and total (bare-bone
functions with limited functionality; see \code{\link{svymean_tukey}} and
\code{\link{svytotal_tukey}} for more capable methods)
}
\details{
\describe{
    \item{Characteristic.}{Population mean or total. Let \eqn{\mu} denote
        the estimated population mean; then, the estimated total is
        given by \eqn{\hat{N} \mu}{Nhat \mu} with
        \eqn{\hat{N} =\sum w_i}{Nhat = sum(w[i])}, where
        summation is over all observations in the sample.}
    \item{Type.}{Two methods/types are available for estimating the
        location \eqn{\mu} (and the scale \eqn{\sigma}; see models, below):
        \describe{
            \item{\code{type = "rht"}:}{a}
            \item{\code{type = "rmw"}:}{a}
        }
    }
    \item{Variance estimation.}{See survey methods:
        \itemize{
            \item \code{\link{svymean_tukey}},
            \item \code{\link{svytotal_tukey}}.
        }
    }
    \item{Psi-function.}{Tukey biweight psi-function with tuning parameter
        \code{k}}
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
weighted_total_tukey(workplace$employment, workplace$weight, k = 9,
    type = "rht")

# Robust weighted M-estimator of the population mean
weighted_mean_tukey(workplace$employment, workplace$weight, k = 12,
    type = "rwm")
}
\references{
Hulliger, B. (1995). Outlier Robust Horvitz-Thompson Estimators.
\emph{Survey Methodology} 21, pp. 79-87.
}
\seealso{
\code{\link{weighted_mean_huber}} and \code{\link{weighted_total_huber}}
}
