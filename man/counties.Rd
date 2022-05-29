\name{counties}
\alias{counties}
\docType{data}
\title{Data on a Simple Random Sample of 100 Counties in the U.S.}
\description{
    Data from a simple random sample (without replacement) of 100 of the 3141
    counties in the United Stated (U.S. Bureau of the Census, 1994).
}
\usage{data(counties)}
\format{
    A \code{data.frame} with 100 observations on the following variables:
    \describe{
        \item{\code{state}}{state, \code{[character]}.}
        \item{\code{county}}{county, \code{[character]}.}
        \item{\code{landarea}}{land area, 1990 (square miles),
            \code{[double]}.}
        \item{\code{totpop}}{population total, 1992, \code{[double]}.}
        \item{\code{unemp}}{number of unemployed persons, 1991,
            \code{[double]}.}
        \item{\code{farmpop}}{farm population, 1990, \code{[double]}.}
        \item{\code{numfarm}}{number of farms, 1987, \code{[double]}.}
        \item{\code{farmacre}}{acreage in farms, 1987, \code{[double]}.}
        \item{\code{weights}}{sampling weight, \code{[double]}.}
        \item{\code{fpc}}{finite population corretion, \code{[double]}.}
    }
}
\details{
    The data (and 10 additional variables) are published in Lohr (1999,
    Appendix C).
}
\source{
Lohr, S. L. (1999). \emph{Sampling: Design and Analysis},
    Pacific Grove (CA): Duxbury Press.
}
\examples{
data(counties)

\dontrun{
# survey design for counties data (pkg survey is required)
library(survey)
dn <- svydesign(ids = ~1, fpc = ~fpc, weights = ~weights, data = counties)
}
}
\keyword{datasets}
