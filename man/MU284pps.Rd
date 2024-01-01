\name{MU284pps}
\alias{MU284pps}
\docType{data}
\title{PPS Sample From the MU284 Population}
\description{
    Probability-proportional-to-size sample (PPS) without replacement of
    municipalities from the MU284 population in Särndal et al. (1992).
    The sample inclusion probabilities are proportional to the
    population size in 1975 (variable P75).
}
\usage{data(MU284pps)}
\format{
    A \code{data.frame} with 60 observations on the following variables:
    \describe{
        \item{\code{LABEL}}{identifier variable, \code{[integer]}.}
        \item{\code{P85}}{1985 population size (in thousands),
            \code{[double]}.}
        \item{\code{P75}}{1975 population size (in thousands),
            \code{[double]}.}
        \item{\code{RMT85}}{Revenues from the 1985 municipal taxation
            (in millions of kronor), \code{[double]}.}
        \item{\code{CS82}}{number of Conservative seats in municipal council,
            \code{[double]}.}
        \item{\code{SS82}}{number of Social-Democrat seats in municipal
            council (1982), \code{[double]}.}
        \item{\code{S82}}{total number of seats in municipal council (1982),
            \code{[double]}.}
        \item{\code{ME84}}{number of municipal employees in 1984,
            \code{[double]}.}
        \item{\code{REV84}}{real estate values according to 1984 assessment
            (in millions of kronor), \code{[double]}.}
        \item{\code{REG}}{geographic region indicator, \code{[integer]}.}
        \item{\code{CL}}{cluster indicator (a cluster consists of a set of
            neighbouring municipalities), \code{[integer]}.}
        \item{\code{weights}}{sampling weights, \code{[double]}.}
        \item{\code{pi}}{sample inclusion probability, \code{[double]}.}
    }
}
\details{
    The MU284 population of Särndal et al. (1992, Appendix B) is a
    dataset with observations on the 284 municipalities in Sweden in the
    late 1970s and early 1980s. The \code{MU284} \emph{population} data
    are available in the \pkg{sampling} package of Tillé and Matei (2021).

    The data frame \code{MU284pps} is a probability-proportional-to-size
    sample (PPS) without replacement from the MU284 population.
    The sample inclusion probabilities are proportional to the
    population size in 1975 (variable P75). The sample has been
    selected by Brewer’s method; see Tillé (2006, Chap. 7).
    The sampling weight (inclusion probabilities) are calibrated to
    the population size and the population total of P75.
}
\source{
Särndal, C.-E., Swensson, B. and Wretman, J. (1992).
    \emph{Model Assisted Survey Sampling}, New York: Springer-Verlag.

Tillé, Y. and Matei, A. (2021). \emph{sampling: Survey Sampling}.
    R package version 2.9. \url{https://CRAN.R-project.org/package=sampling}

Tillé, Y. (2006). \emph{Sampling Algorithms}. New York: Springer-Verlag.
}
\seealso{
    \code{\link{MU284strat}}
}
\examples{
head(MU284pps)

library(survey)
# Survey design with inclusion probabilities proportional to size
dn <- if (packageVersion("survey") >= "4.2") {
        # survey design with pre-calibrated weights
        svydesign(ids = ~LABEL, fpc = ~pi, data = MU284pps, pps = "brewer",
                  calibrate.formula = ~1)
    } else {
        # legacy mode
        svydesign(ids = ~LABEL, fpc = ~pi, data = MU284pps, pps = "brewer")
    }
}
\keyword{datasets}
