\name{class_svyreg_rob}
\alias{class_svyreg_rob}
\alias{svyreg_rob}
\alias{print.svyreg_rob}
\alias{summary.svyreg_rob}
\alias{coef.svyreg_rob}
\alias{vcov.svyreg_rob}
\alias{residuals.svyreg_rob}
\alias{fitted.svyreg_rob}
\alias{robweights.svyreg_rob}
\alias{plot.svyreg_rob}
\title{Utility Functions for Objects of Class svyreg_rob}
\usage{
\method{print}{svyreg_rob}(x, digits = max(3L, getOption("digits") - 3L), ...)

\method{summary}{svyreg_rob}(object, digits = max(3L,
    getOption("digits") - 3L), ...)

\method{coef}{svyreg_rob}(object, ...)

\method{vcov}{svyreg_rob}(object, ...)

\method{residual}{svyreg_rob}(object, ...)

\method{fitted}{svyreg_rob}(object, ...)

\method{robweights}{svyreg_rob}(object)

\method{plot}{svyreg_rob}(x, which = 1:5, ...)
}
\arguments{
\item{x}{object of class \code{svyreg_rob}.}

\item{digits}{\code{[integer]} minimal number of significant digits.}

\item{...}{additional arguments passed to the method.}

\item{object, }{x object of class \code{svyreg_rob}.}

\item{which}{indicating which plots to be drawn; if a subset of the plots
    is required, you can specify a subset of the numbers \code{1:5}.}
}
\description{
Methods and utility functions for objects of classs code{svyreg_rob}
}
\details{
Utility functions
\itemize{
   \item \code{summary} gives a summary of the estimation properties
   \item \code{robweights} extracts the robustness weights
   \item \code{coef} extracts the estimates
   \item \code{SE} extracts the (estimated) standard error
   \item \code{vcov} extracts the (estimated) covariance matrix
   \item \code{residuals} extracts the residuals
   \item \code{fitted} extracts the fitted values
}
}
