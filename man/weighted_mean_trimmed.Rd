\name{weighted_mean_trimmed}
\alias{weighted_mean_trimmed}
\alias{weighted_total_trimmed}
\title{Weighted Trimmed Mean and Total (bare-bone functions)}
\usage{
weighted_mean_trimmed(x, w, LB = 0.05, UB = 1 - LB, info = FALSE,
                      na.rm = FALSE)
weighted_total_trimmed(x, w, LB = 0.05, UB = 1 - LB, info = FALSE,
                       na.rm = FALSE)
}
\arguments{
    \item{x}{\code{[numeric vector]} data.}
    \item{w}{\code{[numeric vector]} weights (same length as \code{x}).}
    \item{LB}{\code{[double]} lower bound of trimming such that
        \eqn{0 \leq} \code{LB} \eqn{<} \code{UB} \eqn{\leq 1}.}
    \item{UB}{\code{[double]} upper bound of trimming such that
        \eqn{0 \leq} \code{LB} \eqn{<} \code{UB} \eqn{\leq 1}.}
    \item{info}{\code{[logical]} indicating whether additional information
        should be returned (default: \code{FALSE}).}
    \item{na.rm}{\code{[logical]} indicating whether \code{NA} values should
        be removed before the computation proceeds (default: \code{FALSE}).}
}
\description{
    Weighted trimmed mean and total (bare-bone functions with limited
    functionality; see \code{\link{svymean_trimmed}} and
    \code{\link{svytotal_trimmed}} for more capable methods)
}
\details{
    \describe{
        \item{Characteristic.}{Population mean or total. Let \eqn{\mu} denote
            the estimated trimmed population mean; then, the estimated trimmed
            population total is given by \eqn{\hat{N} \mu}{Nhat \mu} with
            \eqn{\hat{N} =\sum w_i}{Nhat = sum(w[i])}, where
            summation is over all observations in the sample.}
        \item{Trimming.}{The methods trims the \code{LB}\eqn{~\cdot 100\%}
            of the smallest observations and the
            (1 - \code{UB})\eqn{~\cdot 100\%} of the largest observations
            from the data.}
        \item{Variance estimation.}{See survey methods:
            \itemize{
                \item \code{\link{svymean_trimmed}},
                \item \code{\link{svytotal_trimmed}}.
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

\code{\link{svymean_trimmed}} and \code{\link{svytotal_trimmed}}
}
\examples{
head(workplace)

# Estimated trimmed population total (5\% symmetric trimming)
weighted_total_trimmed(workplace$employment, workplace$weight, LB = 0.05,
    UB = 0.95)

# Estimated trimmed population mean (5\% trimming at the top of the distr.)
weighted_mean_trimmed(workplace$employment, workplace$weight, UB = 0.95)
}
