\name{weighted-m-estimator}
\alias{weighted_mean_huber}
\alias{weighted_total_huber}
\alias{weighted_mean_tukey}
\alias{weighted_total_tukey}
\title{Weighted Huber and Tukey Mean and Total (bare-bone functions)}
\usage{
weighted_mean_huber(x, w, k, type = "rwm", asym = FALSE, info = FALSE,
                    na.rm = FALSE, verbose = TRUE, ...)
weighted_total_huber(x, w, k, type = "rwm", asym = FALSE, info = FALSE,
                     na.rm = FALSE, verbose = TRUE, ...)
weighted_mean_tukey(x, w, k, type = "rwm", info = FALSE, na.rm = FALSE,
                    verbose = TRUE, ...)
weighted_total_tukey(x, w, k, type = "rwm", info = FALSE, na.rm = FALSE,
                     verbose = TRUE, ...)
}
\arguments{
    \item{x}{\code{[numeric vector]} data.}
    \item{w}{\code{[numeric vector]} weights (same length as \code{x}).}
    \item{k}{\code{[double]} robustness tuning constant
        (\eqn{0 < k \leq \infty}{0 < k <= Inf}).}
    \item{type}{\code{[character]} type of method: \code{"rwm"} or
        \code{"rht"}; see below (default: \code{"rwm"}).}
    \item{asym}{\code{[logical]} toggle for asymmetric Huber psi-function
        (default: \code{FALSE}).}
    \item{info}{\code{[logical]} indicating whether additional information
        should be returned (default: \code{FALSE}).}
    \item{na.rm}{\code{[logical]} indicating whether \code{NA} values should
        be removed before the computation proceeds (default: \code{FALSE}).}
    \item{verbose}{\code{[logical]} indicating whether additional
        information is printed to the console (default: \code{TRUE}).}
    \item{\dots}{additional arguments passed to the method (e.g.,
        \code{maxit}: maxit number of iterations, etc.).}
}
\description{
    Weighted Huber and Tukey \emph{M}-estimator of the mean and total
    (bare-bone function with limited functionality; see
    \code{\link{svymean_huber}},  \code{\link{svymean_tukey}},
    \code{\link{svytotal_huber}}, and \code{\link{svytotal_tukey}} for more
    capable methods)
}
\details{
    \describe{
        \item{Characteristic.}{Population mean or total. Let \eqn{\mu}
            denote the estimated population mean; then, the estimated
            total is given by \eqn{\hat{N} \mu}{Nhat \mu} with
            \eqn{\hat{N} =\sum w_i}{Nhat = sum(w[i])}, where
            summation is over all observations in the sample.}
        \item{Type.}{Two methods/types are available for estimating the
            location \eqn{\mu}:
            \describe{
                \item{\code{type = "rwm" (default)}:}{robust weighted
                    \emph{M}-estimator of the population mean and total,
                    respectively. This estimator is recommended for sampling
                    designs whose inclusion probabilities are \emph{not}
                    proportional to some measure of size. [Legacy note: In an
                    earlier version, the method \code{type = "rwm"} was called
                    \code{"rhj"}; the type \code{"rhj"} is now silently
                    converted to \code{"rwm"}]
                }
                \item{\code{type = "rht"}:}{robust Horvitz-Thompson
                    \emph{M}-estimator of the population mean and total,
                    respectively. This estimator is recommended for
                    proportional-to-size sampling designs.}
            }
        }
        \item{Variance estimation.}{See the related but more capable
            functions:
            \itemize{
                \item \code{\link{svymean_huber}} and
                    \code{\link{svymean_tukey}},
                \item \code{\link{svytotal_huber}} and
                    \code{\link{svytotal_tukey}}.
            }
        }
        \item{Psi-function.}{By default, the \code{Huber} or \code{Tukey}
            psi-function are used in the specification of the M-estimators. For
            the Huber estimator, an asymmetric version of the Huber
            psi-function can be used by setting the argument
            \code{asym = TRUE} in the function call.}
    }
}
\section{Failure of convergence}{
    By default, the method assumes a maximum number of \code{maxit = 100}
    iterations and a numerical tolerance criterion to stop the iterations of
    \code{tol = 1e-05}. If the algorithm fails to converge, you may
    consider changing the default values; see \code{\link{svyreg_control}}.
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
\references{
Hulliger, B. (1995). Outlier Robust Horvitz-Thompson Estimators.
    \emph{Survey Methodology} \bold{21}, 79--87.
}
\seealso{
\link[=robsurvey-package]{Overview} (of all implemented functions)
}
\examples{
head(workplace)

# Robust Horvitz-Thompson M-estimator of the population total
weighted_total_huber(workplace$employment, workplace$weight, k = 9,
    type = "rht")

# Robust weighted M-estimator of the population mean
weighted_mean_huber(workplace$employment, workplace$weight, k = 12,
    type = "rwm")
}
