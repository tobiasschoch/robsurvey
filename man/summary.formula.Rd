\name{summary.formula}
\alias{summary.formula}
\title{Weighted Five-Number Summary of a Variable}
\usage{
\method{summary}{formula}(object, design, na.rm = FALSE, ...)
}
\arguments{
\item{object}{one-sided \code{[formula]} for which a summary is desired, e.g.,
    \code{~payroll}.}
\item{design}{an object of class \code{survey.design} or \code{survey.design2}.}
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
data(workplace)

library(survey)
# Survey design for simple random sampling without replacement
dn <- svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight,
    data = workplace)

summary(~payroll, dn)
}
