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
    \item{x}{\code{[numeric vector]} data.}
    \item{k}{\code{[double]} robustness tuning constant
        (\eqn{0 < k \leq \infty}{0 < k <= Inf}).}
    \item{a}{\code{[double]} robustness tuning constant
        (\eqn{0 \leq a \leq \infty}{0 <= a <= Inf}); see details below.}
    \item{b}{\code{[double]} robustness tuning constant
        (\eqn{0 < b \leq \infty}{0 < b <= Inf}; see details below.}
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
    distance of the i-th observation. The Simpson et al. (1992) type of
    weight is defined as
    \eqn{\min \{1, (b/d_i)^{a/2}\}}{min[1, (b/d[i])^(a/2)]}, where
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
\value{
    Numerical vector of weights
}
\references{
Simpson, D. G., Ruppert, D. and Carroll, R.J. (1992). On One-Step GM
    Estimates and Stability of Inferences in Linear Regression.
    \emph{Journal of the American Statistical Association} \bold{87}, 439--450.
    \doi{10.2307/2290275}
}
\seealso{
\link[=robsurvey-package]{Overview} (of all implemented functions)

\code{\link{svyreg_huberM}}, \code{\link{svyreg_huberGM}},
\code{\link{svyreg_tukeyM}} and \code{\link{svyreg_tukeyGM}}
}
\examples{
head(flour)

# standardized distance from median (copper content in wholemeal flour)
x <- flour$copper
z <- abs(x - median(x)) / mad(x)

# plot of weight functions vs. distance
plot(z, huberWgt(z, k = 3), ylim = c(0, 1), xlab = "distance",
     ylab = "weight")
points(z, tukeyWgt(z, k = 6), pch = 2, col = 2)
points(z, simpsonWgt(z, a = Inf, b = 3), pch = 3, col = 4)
legend("topright", c("huberWgt(k = 3)", "tukeyWgt(k = 6)",
       "simpsonWgt(a = Inf, b = 3)"), pch = 1:3, col = c(1, 2, 4))
}
