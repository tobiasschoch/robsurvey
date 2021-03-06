\name{svymean_reg_huber}
\alias{svymean_reg_huber}
\alias{svytotal_reg_huber}
\title{Robust Regression Predictor of the Mean}
\usage{
svymean_reg_huber(object, mean_auxiliary, k)
svytotal_reg_huber()
}
\arguments{
\item{object}{fitted regression model (object of class \code{svyreg}).}
\item{mean_auxiliary}{\code{[numeric vector]} population means of the
    auxiliary variables.}
\item{k}{\code{[double]} robustness tuning constant
    (\eqn{0 < k \leq \infty}{0 < k <= Inf}).}
}
\value{
object of class \code{svystat.rob}
}
\description{
regression estimator of the mean
}
\details{
Details
}
\examples{
}
