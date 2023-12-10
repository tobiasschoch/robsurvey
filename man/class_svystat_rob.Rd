\name{class_svystat_rob}
\alias{class_svystat_rob}
\alias{svystat_rob}
\alias{summary.svystat_rob}
\alias{coef.svystat_rob}
\alias{SE.svystat_rob}
\alias{vcov.svystat_rob}
\alias{residuals.svystat_rob}
\alias{fitted.svystat_rob}
\alias{robweights}
\alias{robweights.svystat_rob}
\alias{print.svystat_rob}
\alias{scale.svystat_rob}
\alias{mse}
\title{Utility Functions for Objects of Class svystat_rob}
\usage{
mse(object)
\method{summary}{svystat_rob}(object, digits = max(3L,
    getOption("digits") - 3L), ...)
\method{coef}{svystat_rob}(object, ...)
\method{SE}{svystat_rob}(object, ...)
\method{vcov}{svystat_rob}(object, ...)
\method{scale}{svystat_rob}(x, ...)
\method{residuals}{svystat_rob}(object, ...)
\method{fitted}{svystat_rob}(object, ...)
robweights(object)
\method{robweights}{svystat_rob}(object)
\method{print}{svystat_rob}(x, digits = max(3L, getOption("digits") - 3L), ...)
}
\arguments{
    \item{object}{object of class \code{svystat_rob}.}
    \item{digits}{\code{[integer]} minimal number of significant digits.}
    \item{\dots}{additional arguments passed to the method.}
    \item{x}{object of class \code{svystat_rob}.}
}
\description{
    Methods and utility functions for objects of class \code{svystat_rob}.
}
\details{
    Package \pkg{survey} must be attached to the search path in order to use
    the functions (see \code{\link{library}} or \code{\link{require}}).

    Utility functions:
    \itemize{
        \item \code{mse} computes the estimated risk (mean square
            error) in presence of representative outliers; see also
            \code{\link{mer}}
        \item \code{summary} gives a summary of the estimation properties
        \item \code{robweights} extracts the robustness weights
        \item \code{coef} extracts the estimate of location
        \item \code{SE} extracts the (estimated) standard error
        \item \code{vcov} extracts the (estimated) covariance matrix
        \item \code{residuals} extracts the residuals
        \item \code{fitted} extracts the fitted values
    }
}
