\name{svymean_reg}
\alias{svymean_reg}
\alias{svytotal_reg}
\title{Regression Predictor of the Mean (Generalized Regression Estimator)}
\usage{
svymean_reg(object, auxiliary, k = Inf,
    psi = c("Huber", "Huberasym", "Tukey"), check.names = TRUE)
svytotal_reg()
}
\arguments{
\item{object}{fitted regression model (object of class \code{svyreg}).}
\item{auxiliary}{\code{[numeric vector]} population means of the
    auxiliary variables.}
\item{k}{tuning constant.}
\item{psi}{psi-function.}
\item{check.names}{\code{[logical]} should the names be checked.}
}
\description{
%FIXME
description
}
\details{
%FIXME
details
}
\value{
object of class \code{svystat_rob}
}
\examples{
}
