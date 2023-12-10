\name{flour}
\alias{flour}
\docType{data}
\title{Measurement of Copper Content in Wholemeal Flour}
\description{
    Measurement of copper content in wholemeal flour (measured in parts per
    million).
}
\usage{data(flour)}
\format{
    A \code{data.frame} with 24 observations (sorted in ascending order)
    on the following variables:
    \describe{
        \item{\code{copper}}{copper content \code{[double]}.}
        \item{\code{weight}}{weight \code{[double]}.}
    }
}
\details{
    The data are published in Maronna et al. (2019, p. 2).
}
\source{
Maronna, R. A., Martin, R. D., Yohai, V. J. and Salibián-Barrera, M.
    (2019). \emph{Robust Statistics: Theory and Methods (with R)},
    Hoboken (NJ): John Wiley and Sons, 2nd edition.
    \doi{10.1002/9781119214656}
}
\examples{
head(flour)
}
\keyword{datasets}
