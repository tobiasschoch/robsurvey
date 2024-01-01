\name{svysummary}
\alias{svysummary}
\title{Weighted Five-Number Summary of a Variable}
\usage{
svysummary(object, design, na.rm = FALSE, ...)
}
\arguments{
    \item{object}{one-sided \code{[formula]} for which a summary is desired,
        e.g., \code{~payroll}.}
    \item{design}{an object of class \code{survey.design}; see
        \code{\link{svydesign}}.}
    \item{na.rm}{\code{[logical]} indicating whether \code{NA} values should
        be removed before the computation proceeds (default: \code{FALSE}).}
    \item{\dots}{additional arguments.}
}
\description{
    Weighted five-number summary used for \code{survey.design} and
    \code{survey.design2} objects (similar to \code{base::summary}
    for \code{[numeric vectors]}).
}
\value{
    A weighted five-number summary (numeric variable) or a frequency table
    (factor variable).
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

svysummary(~payroll, dn)
}
