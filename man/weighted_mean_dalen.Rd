\name{weighted_mean_dalen}
\alias{weighted_mean_dalen}
\alias{weighted_total_dalen}
\title{Dalen Estimators of the Mean and Total}
\usage{
weighted_mean_dalen(x, w, censoring, type = "Z2", info = FALSE,
    na.rm = FALSE, verbose = TRUE)
weighted_total_dalen( x, w, censoring, type = "Z2", info = FALSE,
    na.rm = FALSE, verbose = TRUE)
}
\arguments{
    \item{x}{\code{[numeric vector]} data.}
    \item{w}{\code{[numeric vector]} weights (same length as \code{x}).}
    \item{censoring}{\code{[double]} cutoff threshold above which the
        observations are censored.}
    \item{type}{\code{[character]} type of estimator; either \code{"Z2"} or
        \code{"Z3"} (default: \code{"Z2"}).}
    \item{info}{\code{[logical]} indicating whether additional information
        should be returned (default: \code{FALSE}).}
    \item{na.rm}{\code{[logical]} indicating whether \code{NA} values should be
        removed before the computation proceeds (default: \code{FALSE}).}
    \item{verbose}{\code{[logical]} indicating whether additional information
        should be printed to the console (default: \code{FALSE}).}
}
\description{
    Dalén's estimators of the population mean and the population total
    (bare-bone functions with limited functionality).
}
\details{
    Let \eqn{\sum_{i \in s} w_ix_i}{sum(w[i]x[i])} denote the expansion
    estimator of the \eqn{x}-total (summation is over all elements \eqn{i}
    in sample \eqn{s}). The estimators Z2 and Z3 of Dalén (1987) are
    defined as follows.
    \describe{
        \item{Estimator Z2}{
            The estimator Z2 of the population total sums over
            \eqn{\min(c, w_ix_i)}{min(c, w[i]x[i]}; hence, it
            censors the products \eqn{w_ix_i}{w[i]x[i]} to the
            censoring constant \eqn{c} (\code{censoring}). The estimator of
            the population \eqn{x}-mean is is defined as the total divided
            by the population size.
        }
        \item{Estimator Z3}{
            The estimator Z3 of the population total is defined as the sum
            over the elements \eqn{z_i}{z[i]}, which is equal to
            \eqn{z_i = w_ix_i}{z[i] = w[i]x[i]}
            if \eqn{w_iy_i \leq c}{w[i]x[i] <= c} and
            \eqn{z_i = c + (y_i - c/w_i)}{z[i] = c + (y[i] - c/w[i])}
            otherwise.
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
\examples{
data(workplace)

# Dalen's estimator of the total (with censoring threshold: 100000)
weighted_total_dalen(workplace$employment, workplace$weight, 100000)
}
\references{
Dalén, J. (1987). Practical Estimators of a Population Total Which Reduce
    the Impact of Large Observations. R & D Report U/STM 1987:32, Statistics
    Sweden, Stockholm.
}
