\name{weighted_quantile}
\alias{weighted_quantile}
\title{Weighted Quantile}
\usage{
weighted_quantile(x, w, probs, na.rm = FALSE)
}
\arguments{
    \item{x}{\code{[numeric vector]} data.}
    \item{w}{\code{[numeric vector]} weights (same length as \code{x}).}
    \item{probs}{\code{[numeric vector]} vector of probabilities with
        values in \code{[0,1]}.}
    \item{na.rm}{\code{[logical]} indicating whether \code{NA} values
        should be removed before the computation proceeds
        (default: \code{FALSE}).}
}
\description{
    Weighted population quantile.
}
\details{
    \describe{
        \item{Overview.}{\code{weighted_quantile} computes the weighted
            sample quantiles; argument \code{probs} allows vector inputs.}
        \item{Implementation.}{The function is based on a weighted version
            of the quickselect/Find algorithm with the Bentley and McIlroy
            (1993) 3-way partitioning scheme. For very small arrays,
            we use insertion sort.}
       \item{Compatibility.}{For equal weighting, i.e., when all elements in
            \code{w} are equal, \code{weighted_quantile} is identical with
            \code{type = 2} of \code{stats::quantile}; see also
            Hyndman and Fan (1996).}
    }
}
\value{
    Weighted estimate of the population quantiles
}
\references{
Bentley, J. L. and McIlroy, D. M. (1993). Engineering a Sort Function,
    \emph{Software - Practice and Experience} \bold{23}, 1249--1265.
    \doi{10.1002/spe.4380231105}

Hyndman, R.J. and Fan, Y. (1996). Sample Quantiles in Statistical Packages,
    \emph{The American Statistician} \bold{50}, 361--365.
    \doi{10.1080/00031305.1996.10473566}
}
\seealso{
\link[=robsurvey-package]{Overview} (of all implemented functions)

\code{\link{weighted_median}}
}
\examples{
head(workplace)

# Weighted 25\% quantile (1st quartile)
weighted_quantile(workplace$employment, workplace$weight, 0.25)
}
