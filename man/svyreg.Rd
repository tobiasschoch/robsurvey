\name{svyreg}
\alias{svyreg}
\title{Survey Regression Estimator}
\usage{
svyreg(formula, design, var = NULL, na.rm = FALSE)
}
\arguments{
\item{formula}{a \code{[formula]} object (i.e., symbolic description of the
    model)}
\item{design}{an object of class \code{survey.design} or \code{survey.design2}.}
\item{var}{\code{[numeric vector]} heteroscedastic variance (default:
    \code{NULL}, i.e., homoscedastic variance).}
\item{na.rm}{\code{[logical]} indicating whether \code{NA} values should be
    removed before the computation proceeds (default: \code{FALSE}).}
}
\value{
object of class \code{svyreg}
}
\description{
Weighted regression estimator
%FIXME
Details
}
