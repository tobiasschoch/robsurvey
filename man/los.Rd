\name{losdata}
\alias{losdata}
\docType{data}
\title{Length-of-stay (LOS) hospital data}
\description{
A simple random sample of 70 patients in inpatient hospital treatment.
}
\usage{data(workplace)}
\format{
   A \code{data.frame} with data on 70 patients 
   \describe{
      \item{\code{los}}{length of stay (days) \code{[integer]}.}
      \item{\code{weight}}{sampling weight \code{[double]}.}
      \item{\code{fpc}}{finite population correction \code{[double]}.}
   }
}
\details{
The \code{losdata} are a simple random sample without replacement (SRSWOR) of \eqn{n = 70} patients from the (fictive) population of \eqn{N = 2479} patients in inpatient hospital treatment. We have constructed the \code{losdata} as a showcase; though, the LOS measurements are real data that we have taken from the 201 observations in Ruffieux et al. (2000). The original LOS data of Ruffieux et al. (2000) are available in the R package \pkg{robustbase}; see \code{robustbase::data(los)}. Our \code{losdata} are a SRSWOR of size \eqn{n = 70} from the 201 original observations.
}

\source{
Ruffieux et al. (2000) and data.frame \code{los} in the R package \pkg{robustbase}.
}

\references{
Ruffieux, C., F. Paccaud, and A. Marazzi (2000). Comparing rules for truncating hospital length of stay, \emph{Casemix Quarterly} 2.
}

\examples{
data(losdata)
}
\keyword{datasets}
