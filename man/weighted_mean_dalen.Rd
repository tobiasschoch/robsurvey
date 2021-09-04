\name{weighted_mean_dalen}
\alias{weighted_mean_dalen}
\alias{weighted_total_dalen}
\title{Dalen Estimator of the Mean and Total}
\usage{
weighted_mean_dalen(x, w, censoring, verbose = TRUE, info = FALSE,
    na.rm = FALSE)
weighted_total_dalen( x, w, censoring, verbose = TRUE, info = FALSE,
    na.rm = FALSE)
}
\arguments{
\item{x}{\code{[numeric vector]} data.}
\item{w}{\code{[numeric vector]} weights (same length as vector \code{x}).}
\item{censoring}{\code{[double]} threshold above which all observations are
    censored.}
\item{verbose}{\code{[logical]} indicating whether additional information
    should be printed to the console (default: \code{FALSE}).}
\item{info}{\code{[logical]} indicating whether additional information
    should be returned (default: \code{FALSE}).}
\item{na.rm}{\code{[logical]} indicating whether \code{NA} values should be
    removed before the computation proceeds (default: \code{FALSE}).}
}
\description{
Dalén's estimator of the mean and total (bare-bone functions with limited
functionality)
}
\details{
Let \eqn{sum_{i \in s} w_i x_i}{sum(w[i]x[i])} denote the estimator of
the \eqn{x}-total (summation is over all elements \eqn{i} in sample \eqn{s}).
The estimator of Dalén (1987) sums over \eqn{\min(c, w_ix_i)}{min(c, w[i]x[i]};
hence, it censors the products \eqn{w_ix_i}{w[i]x[i]} to the censoring
constant \eqn{c} (\code{censoring}). The estimator of the population
\eqn{x}-mean is is defined similarly.
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
\examples{
data(workplace)

# Dalen's estimator of the total (with censoring threshold: 100000)
weighted_total_dalen(workplace$employment, workplace$weight, 100000)
}
\references{
Dalén, J. (1987). Practical Estimators of a Population Total Which Reduce
the Impact of Large Observations, Research Report, Statistics Sweden.
}
