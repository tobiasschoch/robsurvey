\name{svymean_winsorized}
\alias{svymean_winsorized}
\alias{svymean_k_winsorized}
\alias{svytotal_winsorized}
\alias{svytotal_k_winsorized}
\title{Weighted Winsorized Mean and Total}
\usage{
svymean_winsorized(x, design, LB = 0.05, UB = 1 - LB, na.rm = FALSE,
                   trim_var = FALSE, ...)
svymean_k_winsorized(x, design, k, na.rm = FALSE, trim_var = FALSE, ...)
svytotal_winsorized(x, design, LB = 0.05, UB = 1 - LB, na.rm = FALSE,
                    trim_var = FALSE, ...)
svytotal_k_winsorized(x, design, k, na.rm = FALSE, trim_var = FALSE, ...)
}
\arguments{
    \item{x}{a one-sided \code{[formula]}, e.g., \code{~myVariable}.}
    \item{design}{an object of class \code{survey.design}; see
        \code{\link[survey]{svydesign}}.}
    \item{LB}{\code{[double]} lower bound of winsorization such that
        \eqn{0 \leq} \code{LB} \eqn{<} \code{UB} \eqn{\leq 1}.}
    \item{UB}{\code{[double]} upper bound of winsorization such that
        \eqn{0 \leq} \code{LB} \eqn{<} \code{UB} \eqn{\leq 1}.}
    \item{na.rm}{\code{[logical]} indicating whether \code{NA} values should
        be removed before the computation proceeds (default: \code{FALSE}).}
    \item{trim_var}{\code{[logical]} indicating whether the variance
        should be approximated by the variance estimator of the trimmed mean/
        total (default: \code{FALSE}).}
    \item{k}{\code{[integer]} number of observations to be winsorized at the
        top of the distribution.}
    \item{\dots}{additional arguments (currently not used).}
}
\description{
    Weighted winsorized mean and total
}
\details{
    Package \pkg{survey} must be attached to the search path in order to use
    the functions (see \code{\link{library}} or \code{\link{require}}).

    \describe{
        \item{Characteristic.}{Population mean or total. Let \eqn{\mu}
            denote the estimated winsorized population mean; then, the
            estimated winsorized total is given by
            \eqn{\hat{N} \mu}{Nhat \mu} with
            \eqn{\hat{N} =\sum w_i}{Nhat = sum(w[i])}, where summation
            is over all observations in the sample.}
        \item{Modes of winsorization.}{The amount of winsorization can be
            specified in relative or absolute terms:
        \itemize{
            \item \emph{Relative:} By specifying \code{LB} and \code{UB},
                the method winsorizes the \code{LB}\eqn{~\cdot 100\%}
                of the smallest observations and the
                (1 - \code{UB})\eqn{~\cdot 100\%} of the largest
                observations from the data.
            \item \emph{Absolute:} By specifying argument \code{k} in the
                functions with the "infix" \code{_k_} in their name (e.g.,
                \code{svymean_k_winsorized}), the
                largest \eqn{k} observations are winsorized, \eqn{0<k<n},
                where \eqn{n} denotes the sample size. E.g., \code{k = 2}
                implies that the largest and the second largest observation
                are winsorized.
            }
        }
        \item{Variance estimation.}{Large-sample approximation based on the
            influence function; see Huber and Ronchetti (2009, Chap. 3.3)
            and Shao (1994). Two estimators are available:
            \describe{
                \item{\code{simple_var = FALSE}}{Variance estimator of
                    the winsorized mean/ total. The estimator depends on
                    the estimated probability density function evaluated at
                    the winsorization thresholds, which can be -- depending
                    on the context -- numerically unstable. As a remedy,
                    a simplified variance estimator is available by
                    setting \code{simple_var = TRUE}.}
                \item{\code{simple_var = TRUE}}{Variance is approximated
                    using the variance estimator of the trimmed mean/ total.}
            }
        }
        \item{Utility functions.}{\code{\link[=svystat_rob]{summary}},
            \code{\link[=svystat_rob]{coef}}, \code{\link[=svystat_rob]{SE}},
            \code{\link[=svystat_rob]{vcov}},
            \code{\link[=svystat_rob]{residuals}},
            \code{\link[=svystat_rob]{fitted}} and
            \code{\link[=svystat_rob]{robweights}}.}
        \item{Bare-bone functions.}{See:
            \itemize{
                \item \code{\link{weighted_mean_winsorized}},
                \item \code{\link{weighted_mean_k_winsorized}},
                \item \code{\link{weighted_total_winsorized}},
                \item \code{\link{weighted_total_k_winsorized}}.}
        }
    }
}
\value{
    Object of class \code{\link{svystat_rob}}
}
\references{
Huber, P. J. and Ronchetti, E. (2009). \emph{Robust Statistics},
    New York: John Wiley and Sons, 2nd edition.
    \doi{10.1002/9780470434697}

Shao, J. (1994). L-Statistics in Complex Survey Problems.
    \emph{The Annals of Statistics} \bold{22}, 976--967.
    \doi{10.1214/aos/1176325505}
}
\seealso{
\link[=robsurvey-package]{Overview} (of all implemented functions)

\code{\link{weighted_mean_winsorized}},
\code{\link{weighted_mean_k_winsorized}},
\code{\link{weighted_total_winsorized}} and
\code{\link{weighted_total_k_winsorized}}
}
\examples{
head(workplace)

library(survey)
# Survey design for stratified simple random sampling without replacement
dn <- if (packageVersion("survey") >= "4.2") {
        # survey design with pre-calibrated weights
        svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight,
                  data = workplace, calibrate.formula = ~-1 + strat)
    } else {
        # legacy mode
        svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight,
                  data = workplace)
    }

# Estimated winsorized population mean (5\% symmetric winsorization)
svymean_winsorized(~employment, dn, LB = 0.05)

# Estimated one-sided k winsorized population total (2 observations are
# winsorized at the top of the distribution)
svytotal_k_winsorized(~employment, dn, k = 2)
}
