\name{wBACON}
\alias{wBACON}
\alias{print.robmv}
\alias{distance}
\title{Weighted BACON algorithm: Multivariate robust estimation of location and scatter}
\usage{
wBACON(x, w = NULL, alpha = 0.95, intercept = FALSE, na.rm = FALSE, maxiter = 50, verbose = TRUE) 

distance(object)

\method{print}{robmv}(x, digits = max(3L, getOption("digits") - 3L), ...)
}
\arguments{
\item{x}{\code{[matrix]} or \code{[data.frame]}.}

\item{w}{\code{[numeric]} sampling weight (default \code{w = NULL}).}

\item{alpha}{\code{[double]} tuning constant, \eqn{0 < \alpha < 1}{0 < alpha < 1}; (default: \code{alpha = 0.95}).}

\item{intercept}{\code{[logical]} indicating whether \code{x} contains a regression intercept (default: \code{FALSE}). If \code{TRUE}, \code{wBACON} drops the constant variable.}

\item{na.rm}{\code{[logical]} indicating whether \code{NA} values should be removed before the computation proceeds (default: \code{FALSE}).}

\item{maxiter}{\code{[integer]} maximal number of iterations (default: \code{maxiter = 50}).}

\item{verbose}{\code{[logical]} indicating whether additional information is printed to the console (default: \code{TRUE}).}

\item{digits}{\code{[integer]} minimal number of significant digits.}

\item{...}{additional arguments passed to the method.}

\item{object}{object of class \code{robmv}.}
}
\description{
\code{wBACON} is an iterative method for the computation of multivariate location and scatter (under the assumption of a Gaussian distribution). The algorithm is initialized from a set of uncontaminated data. Then the subset is iteratively refined; i.e., additional observations are included into the subset if their Mahalanobis distance is below some threshold (likewise, observations are removed from the subset if their distance larger than the threshold). This process iterates until the set of good data remain stable. Observations not among the good data are outliers; see Billor et al. (2000). The weighted Bacon algorithm is due to Béguin and Hulliger (2008). 
}
\details{
\itemize{
\item The threshold for the (squared) Mahalanobis distances is defined as the standardized chi-square \eqn{alpha}{alpha} quantile. All observations whose squared Mahalanobis distances is larger than the threshold are regarded as outliers. 
\item The \code{wBACON} cannot deal with missing values. In contrast, function \code{\link[modi]{BEM}} in package \pkg{modi} implements the BACON-EEM algorithm of Béguin and Hulliger (2008), which is tailored to work with outlying and missing values. 
\item If the sampling weights \code{w} are not explicitly specified (i.e., \code{w = NULL}), they are taken to be 1.0.
}

}
\examples{
data(workplace) 

# estimation multivariate location and scatter
wb <- wBACON(workplace[, c("payroll", "employment")], workplace$weights)

# summary of the robust distances
summary(distance(wb))
}
\references{
Billor N., Hadi A.S., Vellemann P.F. (2000). BACON: Blocked Adaptative Computationally efficient Outlier Nominators. \emph{Computational Statistics and Data Analysis} 34, pp. 279-298.

Béguin C., Hulliger B. (2008). The BACON-EEM Algorithm for Multivariate Outlier Detection in Incomplete Survey Data. \emph{Survey Methodology} 34, pp. 91-103.
} 
