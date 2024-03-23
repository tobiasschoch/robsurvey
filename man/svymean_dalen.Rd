\name{svymean_dalen}
\alias{svymean_dalen}
\alias{svytotal_dalen}
\title{Dalen's Estimators of the Population Mean and Total}
\usage{
svymean_dalen(x, design, censoring, type = "Z2", na.rm = FALSE,
              verbose = TRUE, ...)
svytotal_dalen(x, design, censoring, type = "Z2", na.rm = FALSE,
               verbose = TRUE, ...)
}
\arguments{
    \item{x}{a one-sided \code{[formula]}, e.g., \code{~myVariable}.}
    \item{design}{an object of class \code{survey.design}; see
        \code{\link[survey]{svydesign}}.}
    \item{censoring}{\code{[double]} cutoff threshold above which the
        observations are censored.}
    \item{type}{\code{[character]} type of estimator; either \code{"Z2"} or
        \code{"Z3"} (default: \code{"Z2"}).}
    \item{na.rm}{\code{[logical]} indicating whether \code{NA} values should
        be removed before the computation proceeds (default: \code{FALSE}).}
    \item{verbose}{\code{[logical]} indicating whether additional information
            is printed to the console (default: \code{TRUE}).}
    \item{\dots}{additional arguments (currently not used).}
}
\description{
    Dalen's estimators Z2 and Z3 of the population mean and total; see
    \code{\link{weighted_mean_dalen}} for further details.
}
\details{
    Package \pkg{survey} must be attached to the search path in order to use
    the functions (see \code{\link{library}} or \code{\link{require}}).

    \describe{
        \item{Methods/ types}{\code{type = "Z2"} or \code{type = "Z3"}; see
            \code{\link{weighted_mean_dalen}} for more details.}
        \item{Utility functions}{\code{\link[=svystat_rob]{summary}},
            \code{\link[=svystat_rob]{coef}}, \code{\link[=svystat_rob]{SE}},
            \code{\link[=svystat_rob]{vcov}},
            \code{\link[=svystat_rob]{residuals}},
            \code{\link[=svystat_rob]{fitted}},
            \code{\link[=svystat_rob]{robweights}}.}
        \item{Bare-bone functions}{See \code{\link{weighted_mean_dalen}} and
            \code{\link{weighted_total_dalen}}.}
    }
}
\value{
    Object of class \code{\link{svystat_rob}}
}
\references{
Dalén, J. (1987). Practical Estimators of a Population Total Which Reduce
    the Impact of Large Observations. R & D Report U/STM 1987:32, Statistics
    Sweden, Stockholm.
}
\seealso{
\link[=robsurvey-package]{Overview} (of all implemented functions)

\code{\link{svymean_trimmed}}, \code{\link{svytotal_trimmed}},
\code{\link{svymean_winsorized}}, \code{\link{svytotal_winsorized}},
\code{\link{svymean_huber}} and \code{\link{svytotal_huber}}
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

# Dalen's estimator Z3 of the population total
svytotal_dalen(~employment, dn, censoring = 20000, type = "Z3")

# Dalen's estimator Z3 of the population mean
m <- svymean_dalen(~employment, dn, censoring = 20000, type = "Z3")

# Summarize
summary(m)

# Extract estimate
coef(m)

# Extract estimated standard error
SE(m)
}
