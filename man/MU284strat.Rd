\name{MU284strat}
\alias{MU284strat}
\docType{data}
\title{Stratified Sample from the MU284 Population}
\description{
    Stratified simple random sample (without replacement) of municipalities
    from the MU284 population in Särndal et al. (1992). Stratification is by
    geographic region and a take-all stratum (by 1975 population size), which
    includes the big cities Stockholm, Göteborg, and Malmö.
}
\usage{data(MU284strat)}
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
        \item{\code{CL}}{cluster indicator (a cluster consists of a set of
            neighbouring municipalities), \code{[integer]}.}
        \item{\code{REG}}{geographic region indicator, \code{[integer]}.}
        \item{\code{Stratum}}{stratum indicator, \code{[integer]}.}
        \item{\code{weights}}{sampling weights, \code{[double]}.}
        \item{\code{fpc}}{finite population correction, \code{[double]}.}
    }
}
\details{
    The MU284 population of Särndal et al. (1992, Appendix B) is a
    dataset with observations on the 284 municipalities in Sweden in the
    late 1970s and early 1980s. The \code{MU284} \emph{population} data
    are available in the \pkg{sampling} package of Tillé and Matei (2021).

    The population is divided into two parts based on 1975 population
    size (\code{P75}):
        \itemize{
            \item the MU281 population, which consists of the 281 smallest
                municipalities;
            \item the MU3 population of the three biggest municipalities/
                cities in Sweden (Stockholm, Göteborg, and Malmö).
        }
    The three biggest cities take exceedingly large values (representative
    outliers) on almost all of the variables. To account for this, a stratified
    sample has been drawn from the MU284 population using a take-all stratum.
    The sample data, \code{MU284strat}, (of size \eqn{n=60}) consists of
        \itemize{
            \item a stratified simple random sample (without replacement)
                from the MU281 population, where stratification is by
                geographic region (\code{REG}) with proportional sample
                size allocation;
            \item a take-all stratum that includes the three biggest
                cities/ municipalities (population M3).
        }
}
\source{
Särndal, C.-E., Swensson, B. and Wretman, J. (1992).
    \emph{Model Assisted Survey Sampling}, New York: Springer-Verlag.

Tillé, Y. and Matei, A. (2021). \emph{sampling: Survey Sampling}.
    R package version 2.9. \url{https://CRAN.R-project.org/package=sampling}
}
\seealso{
    \code{\link{MU284pps}}
}
\examples{
head(MU284strat)

library(survey)
# Survey design for stratified simple random sampling without replacement
dn <- if (packageVersion("survey") >= "4.2") {
        # survey design with pre-calibrated weights
        svydesign(ids = ~LABEL, strata = ~Stratum, fpc = ~fpc,
                  weights = ~weights, data = MU284strat,
                  calibrate.formula = ~-1 + Stratum)
    } else {
        # legacy mode
        svydesign(ids = ~LABEL, strata = ~Stratum, fpc = ~fpc,
                  weights = ~weights, data = MU284strat)
    }
}
\keyword{datasets}
