\name{wgt_functions}
\alias{huberWgt}
\alias{tukeyWgt}
\alias{simpsonWgt}

\title{Weight Functions (for the M- and GM-Estimators)}
\usage{
huberWgt(x, k = 1.345)
tukeyWgt(x, k = 4.685)
simpsonWgt(x, a, b)
}
\arguments{
\item{x}{\code{[numeric vector]} variable.}
\item{k}{\code{[double]} robustness tuning constant
    (\eqn{0 < k \leq \infty}{0 < k <= Inf}).}
\item{a}{\code{[double]} robustness tuning constant
    (\eqn{0 <= a \leq \infty}{0 <= a <= Inf}; default: \code{a = 1});
    see details below.}
\item{b}{\code{[double]} robustness tuning constant
    (\eqn{0 < b \leq \infty}{0 < b <= Inf}; see details below.}
}
\value{
Numerical vector of weights
}
\description{
Weight functions associated with the Huber and the Tukey biweight
psi-functions; and the weight function of Simpson et al. (1992)
for GM-estimators.
}
\details{
The functions \code{huberWgt} and \code{tukeyWgt} return the weights
associated with the respective psi-function.

The function \code{simpsonWgt} is used (in regression GM-estimators)
to downweight leverage observations (i.e., outliers in the model's design
space). Let \eqn{d_i}{d[i]} denote the (robust) squared Mahalanobis
distance of the i-th observation. The Simpson et al. type of weight is
defined as \eqn{\min \{1, (b/d_i)^{a/2}\}}{min[1, (b/d[i])^(a/2)]}, where
\code{a} and \code{b} are tuning constants.
\itemize{
    \item By default, \code{a = 1}; this choice implies that the weights
        are computed on the basis of the robust Mahalanobis distances.
        Alternative: \code{a = Inf} implies a weight of zero for all
        observations whose (robust) squared Mahalanobis is larger than
        \code{b}.
    \item The tuning constants \code{b} is a threshold on the distances.
}
}
\references{
Simpson, D.G., Ruppert, D., and Carroll, R.J. (1992). On One-Step GM
    Estimates and Stability of Inferences in Linear Regression,
    \emph{Journal of the American Statistical Association} 87, pp. 439-450.
}
\seealso{
\code{\link{svyreg_huber}}, \code{\link{svyreg_huberGM}},
\code{\link{svyreg_tukey}}, and \code{\link{svyreg_tukeyGM}}
}
