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

\method{summary}{svyreg_rob}(object, var = c("design", "model", "compound"),
    digits = max(3L, getOption("digits") - 3L), ...)

\method{coef}{svyreg_rob}(object, ...)

\method{vcov}{svyreg_rob}(object, var = c("design", "model", "compound"), ...)

\method{residuals}{svyreg_rob}(object, ...)

\method{fitted}{svyreg_rob}(object, ...)

\method{robweights}{svyreg_rob}(object)

\method{plot}{svyreg_rob}(x, which = 1:5, ...)
}
\arguments{
\item{x}{object of class \code{svyreg_rob}.}
\item{digits}{\code{[integer]} minimal number of significant digits.}
\item{\dots}{additional arguments passed to the method.}
\item{object}{object of class \code{svyreg_rob}.}
\item{var}{type of variance estimator: \code{"design"}, \code{"model"} or
    \code{"compound"}, default: \code{"design"}.}
\item{which}{indicating which plots to be drawn; if a subset of the plots
    is required, you can specify a subset of the numbers \code{1:5}.}
}
\description{
Methods and utility functions for objects of class \code{svyreg_rob}.

Variance estimation (\code{summary} and \code{vcov}):
\itemize{
    \item \code{"design"}: design-based variance estimator using
        linearization; see Binder (1983)
    \item \code{"model"}: model-based weighted variance estimator (the
        sampling design is ignored)
    \item \code{"compound"}: design-model-based variance estimator; see
        Godambe and Thompson (2009)
}
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
\references{
Binder, D.A. (1983).  On the Variances of Asymptotically Normal Estimators
from Complex Surveys, \emph{International Statistical Review} 51, pp. 279-292.

Godambe, V.P. and M.E. Thompson (2009). Estimating Functions and Survey
Sampling, in: D. Pfeffermann and C.R. Rao (eds.),
\emph{Handbook of Statistics}, vol. 29B, Sample Surveys: Inference and
Analysis, Chapter 26, pp. 83-101, Amsterdam: Elsevier.
}
