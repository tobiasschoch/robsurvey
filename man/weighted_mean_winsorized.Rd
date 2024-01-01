\name{weighted_mean_winsorized}
\alias{weighted_mean_winsorized}
\alias{weighted_mean_k_winsorized}
\alias{weighted_total_winsorized}
\alias{weighted_total_k_winsorized}
\title{Weighted Winsorized Mean and Total (bare-bone functions)}
\usage{
weighted_mean_winsorized(x, w, LB = 0.05, UB = 1 - LB, info = FALSE,
                         na.rm = FALSE)
weighted_mean_k_winsorized(x, w, k, info = FALSE, na.rm = FALSE)
weighted_total_winsorized(x, w, LB = 0.05, UB = 1 - LB, info = FALSE,
                          na.rm = FALSE)
weighted_total_k_winsorized(x, w, k, info = FALSE, na.rm = FALSE)
}
\arguments{
    \item{x}{\code{[numeric vector]} data.}
    \item{w}{\code{[numeric vector]} weights (same length as \code{x}).}
    \item{LB}{\code{[double]} lower bound of winsorization such that
        \eqn{0 \leq} \code{LB} \eqn{<} \code{UB} \eqn{\leq 1}.}
    \item{UB}{\code{[double]} upper bound of winsorization such that
        \eqn{0 \leq} \code{LB} \eqn{<} \code{UB} \eqn{\leq 1}.}
    \item{info}{\code{[logical]} indicating whether additional information
        should be returned (default: \code{FALSE}).}
    \item{na.rm}{\code{[logical]} indicating whether \code{NA} values should
        be removed before the computation proceeds (default: \code{FALSE}).}
    \item{k}{\code{[integer]} number of observations to be winsorized at the
        top of the distribution.}
}
\description{
    Weighted winsorized mean and total (bare-bone functions with limited
    functionality; see \code{\link{svymean_winsorized}} and
    \code{\link{svytotal_winsorized}} for more capable methods)
}
\details{
    \describe{
        \item{Characteristic.}{Population mean or total. Let \eqn{\mu}
            denote the estimated winsorized population mean; then, the
            estimated population total is given by \eqn{\hat{N} \mu}{Nhat \mu}
            with \eqn{\hat{N} =\sum w_i}{Nhat = sum(w[i])}, where summation
            is over all observations in the sample.}
        \item{Modes of winsorization.}{The amount of winsorization can be
            specified in relative or absolute terms:
            \itemize{
                \item \emph{Relative:} By specifying \code{LB} and \code{UB},
                    the methods winsorizes the \code{LB}\eqn{~\cdot 100\%}
                    of the smallest observations and the
                    (1 - \code{UB})\eqn{~\cdot 100\%} of the largest
                    observations from the data.
                \item \emph{Absolute:} By specifying argument \code{k} in the
                    functions with the "infix" \code{_k_} in their name, the
                    largest \eqn{k} observations are winsorized, \eqn{0<k<n},
                    where \eqn{n} denotes the sample size. E.g., \code{k = 2}
                    implies that the largest and the second largest
                    observation are winsorized.
            }
        }
        \item{Variance estimation.}{See survey methods:
            \itemize{
                \item \code{\link{svymean_winsorized}},
                \item \code{\link{svytotal_winsorized}},
                \item \code{\link{svymean_k_winsorized}},
                \item \code{\link{svytotal_k_winsorized}}.
            }
        }
    }
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
\seealso{
\link[=robsurvey-package]{Overview} (of all implemented functions)

\code{\link{svymean_winsorized}}, \code{\link{svymean_k_winsorized}},
\code{\link{svytotal_winsorized}} and \code{\link{svytotal_k_winsorized}}
}
\examples{
head(workplace)

# Estimated winsorized population mean (5\% symmetric winsorization)
weighted_mean_winsorized(workplace$employment, workplace$weight, LB = 0.05)

# Estimated one-sided k winsorized population total (2 observations are
# winsorized at the top of the distribution)
weighted_total_k_winsorized(workplace$employment, workplace$weight, k = 2)
}
