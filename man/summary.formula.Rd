\name{summary.formula}
\alias{summary.formula}
\title{Weighted five-number summary of a variable}
\usage{
\method{summary}{formula}(object, design, na.rm = FALSE, ...)
}
\arguments{
\item{object}{\code{[character]} name of the variable for which a summary is desired.}

\item{design}{an object of class \code{survey.design} or \code{survey.design2}.}

\item{na.rm}{\code{[logical]} indicating whether \code{NA} values should be removed before the computation proceeds (default: \code{FALSE}).}

\item{...}{additional arguments.}
}
\value{
A numerical summary
}
\description{
Weighted five-number summary of a variable (similar to \code{base::summary} for \code{[numeric vectors]})
}
\details{
A weighted five-number summary (numeric variable) or a frequency table (factor variable).
}
