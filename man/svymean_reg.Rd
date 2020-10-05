\name{svymean_reg}
\alias{svymean_reg}
\alias{svymean_reg.svyreg_rob}
\title{Regression predictor of the mean (generalized regression estimator)}
\usage{
svymean_reg(object, ...)

\method{svymean_reg}{svyreg_rob}(object, mean_auxiliary, ...)
}
\arguments{
\item{object}{fitted regression model (object of class \code{svyreg}).}

\item{...}{additional arguments passed to the method.}

\item{mean_auxiliary}{\code{[numeric vector]} population means of the auxiliary variables.}
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
