\name{svymean_trimmed}
\alias{svymean_trimmed}
\alias{svytotal_trimmed}
\title{Weighted Trimmed Mean and Total}
\usage{
svymean_trimmed(x, design, LB = 0.05, UB = 1 - LB, na.rm = FALSE)
svytotal_trimmed(x, design, LB = 0.05, UB = 1 - LB, na.rm = FALSE)
}
\arguments{
\item{x}{a one-sided \code{[formula]}, e.g., \code{~myVariable}.}
\item{design}{an object of class \code{survey.design} or \code{survey.design2}.}
\item{LB}{\code{[double]} lower bound of trimming such that
    \eqn{0 \leq} \code{LB} \eqn{<} \code{UB} \eqn{\leq 1}.}
\item{UB}{\code{[double]} upper bound of trimming such that
    \eqn{0 \leq} \code{LB} \eqn{<} \code{UB} \eqn{\leq 1}.}
\item{na.rm}{\code{[logical]} indicating whether \code{NA} values should
    be removed before the computation proceeds (default: \code{FALSE}).}
}
\value{
object of class \code{\link{svystat_rob}}
}
\description{
Weighted trimmed population mean and total.
}
\details{
\describe{
    \item{Characteristic.}{Population mean or total. Let \eqn{\mu} denote
        the estimated trimmed population mean; then, the estimated trimmed
        total is given by \eqn{\hat{N} \mu}{Nhat \mu} with
        \eqn{\hat{N} =\sum w_i}{Nhat = sum(w[i])}, where
        summation is over all observations in the sample.}
    \item{Trimming.}{The methods trims the \code{LB}\eqn{~\cdot 100\%}
        percentage of the smallest observations and the
        (1 - \code{UB})\eqn{~\cdot 100\%} percentage of the largest
        observations from the data.}
    \item{Variance estimation.}{Taylor linearization.}
    \item{Utility functions.}{\code{\link[=svystat_rob]{summary}},
        \code{\link[=svystat_rob]{coef}}, \code{\link[=svystat_rob]{SE}},
        \code{\link[=svystat_rob]{vcov}}, \code{\link[=svystat_rob]{residuals}},
        \code{\link[=svystat_rob]{fitted}}, and
        \code{\link[=svystat_rob]{robweights}}.}
    \item{Bare-bone functions.}{See \code{\link{weighted_mean_trimmed}} and
        \code{\link{weighted_total_trimmed}}.}
}
}
\references{
Huber, P.J. (1981). \emph{Robust Statistics}, New York: John Wiley and Sons.
Shao, J. (1994). L-Statistics in Complex Survey Problems.
\emph{Annals of Statistics} 22, pp. 976-967.
}
\examples{
data(workplace)

library(survey)
# Survey design for simple random sampling without replacement
dn <- svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight,
    data = workplace)

# Estimated trimmed population total (5\% symmetric trimming)
svytotal_trimmed(~employment, dn, LB = 0.05, UB = 0.95)

# Estimated trimmed population mean (5\% trimming at the top of the distr.)
svymean_trimmed(~employment, dn, UB = 0.95)
}
\seealso{
\code{\link{weighted_mean_trimmed}} and \code{\link{weighted_total_trimmed}}
}
