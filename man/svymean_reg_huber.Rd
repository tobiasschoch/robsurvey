\name{svymean_reg_huber}
\alias{svymean_reg_huber}
\title{Robust regression predictor of the mean}
\usage{
svymean_reg_huber(object, mean_auxiliary, k)
}
\arguments{
\item{object}{fitted regression model (object of class \code{svyreg}).}

\item{mean_auxiliary}{\code{[numeric vector]} population means of the auxiliary variables.}

\item{k}{\code{[double]} robustness tuning constant (\eqn{0 < k \leq \infty}{0 < k <= Inf}).}
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
\section{Failure of convergence}{
 
By default, the method assumes a maximum number of \code{maxit = 100} iterations and a numerical tolerance criterion to stop the iterations of \code{tol = 1e-05}. You can run the code with specifications other than the default values by specifying the arguments \code{maxit} and/or \code{tol} in the function call; see also \code{\link{svyreg_control}}.
}

