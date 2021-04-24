\name{svymean_reg}
\alias{svymean_reg}
\title{Regression Predictor of the Mean (Generalized Regression Estimator)}
\usage{
svymean_reg(object, auxiliary, k = NULL, check.names = TRUE)
}
\arguments{
\item{object}{fitted regression model (object of class \code{svyreg}).}
\item{auxiliary}{\code{[numeric vector]} population means of the
    auxiliary variables.}
\item{k}{\code{[double]} robustness tuning constant
    (\eqn{0 < k \leq \infty}{0 < k <= Inf}; default: \code{k = 1.5}).}
\item{check.names}{\code{[logical]} should the names be checked.}
}
\value{
object of class \code{svystat.rob}
}
\description{
Regression predictor of the mean
}
\details{
Details
}
