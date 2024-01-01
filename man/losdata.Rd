\name{losdata}
\alias{losdata}
\docType{data}
\title{Length-of-Stay (LOS) Hospital Data}
\description{
    A simple random sample of 70 patients in inpatient hospital treatment.
}
\usage{data(losdata)}
\format{
    A \code{data.frame} with data on the following variables:
    \describe{
        \item{\code{los}}{length of stay (days) \code{[integer]}.}
        \item{\code{weight}}{sampling weight \code{[double]}.}
        \item{\code{fpc}}{finite population correction \code{[double]}.}
    }
}
\details{
    The \code{losdata} are a simple random sample without replacement (SRSWOR)
    of size \eqn{n = 70} patients from the (fictive) population of
    \eqn{N = 2479} patients in inpatient hospital treatment. We have
    constructed the \code{losdata} as a showcase; though, the LOS
    measurements are real data that we have taken from the 201 observations
    in Ruffieux et al.  (2000). The original LOS data of Ruffieux et al.
    (2000) are available in the R package \pkg{robustbase}; see
    \code{robustbase::data(los)}.  Our \code{losdata} are a SRSWOR of
    size \eqn{n = 70} from the 201 original observations.

Ruffieux et al. (2000) and data.frame \code{los} in the R
    package \pkg{robustbase}.
}
\source{
Ruffieux, C., Paccaud, F. and Marazzi, A. (2000). Comparing rules for
    truncating hospital length of stay. \emph{Casemix Quarterly} \bold{2}.
}
\examples{
head(losdata)

library(survey)
# Survey design for simple random sampling without replacement
dn <- if (packageVersion("survey") >= "4.2") {
        # survey design with pre-calibrated weights
        svydesign(ids = ~1, fpc = ~fpc, weights = ~weight, data = losdata,
                  calibrate.formula = ~1)
    } else {
        # legacy mode
        svydesign(ids = ~1, fpc = ~fpc, weights = ~weight, data = losdata)
    }
}
\keyword{datasets}
