\name{svymean_reg}
\alias{svymean_reg}
\alias{svytotal_reg}
\title{Generalized Regression Estimator (GREG) of the Mean and Total}
\usage{
svymean_reg(object, auxiliary, type, k = NULL, check.names = TRUE,
    na.rm = FALSE, keep_object = TRUE)
svytotal_reg(object, auxiliary, type, k = NULL, check.names = TRUE,
    na.rm = FALSE, keep_object = TRUE)
}
\arguments{
    \item{object}{of class \code{[svyreg_rob]}, e.g., result of the Huber
        regression \eqn{M}-estimator \code{svyreg_huberM}.}
    \item{auxiliary}{\code{[numeric]} mean or total of the auxiliary
        variables.}
    \item{type}{\code{[character]} type of estimator.}
    \item{k}{\code{[numeric]} robustness tuning constant of the
        \code{psi}-function used in the bias-correction term in the GREG
        (default: \code{Inf}, which implies a non-robust GREG).}
    \item{check.names}{\code{[logical]} if \code{TRUE}, the names of
        \code{auxiliary} is checked against the independent variables of the
        fitted model \code{object} (default: \code{TRUE}.}
    \item{na.rm}{\code{[logical]} indicating whether \code{NA} values should
        be removed before the computation proceeds (default: \code{FALSE}).}
    \item{keep_object}{\code{[logical]} if \code{TRUE}, \code{object} is
        returned.}
}
\description{
(Robust) Generalized Regression Estimator (GREG) of the Mean and Total
}
\details{
    \describe{
        \item{How to}{
            Text
        }
        \item{Types}{
            \itemize{
                \item \code{"projective"} (Särndal et al., 1992)
                \item \code{"ADU"} (Särndal et al., 1992; Gwet and Rivest, 1992)
                \item \code{"robust"}
                \item \code{"BR"} (Beaumont and Rivest, 2009)
                \item \code{"lee"} (Lee, 1991; 1992)
                \item \code{"duchesne"} (Duchesne, 1999)
            }
        }
    }
}
\value{
Object of class \code{\link{svystat_rob}}
}
\references{
Beaumont, J.-F. and L.-P. Rivest (2009). Dealing with outliers in survey data,
    in \emph{Sample Surveys: Theory, Methods and Inference},
    ed. by D. Pfeffermann and C. R. Rao,
    Amsterdam: Elsevier, vol. 29A of \emph{Handbook of Statistics},
    chap. 11, 247--280.

Chambers, R. (1986). Outlier Robust Finite Population Estimation,
    \emph{Journal of the American Statistical Association} \bold{81},
    1063--1069.

Duchesne, P. (1999). Robust calibration estimators,
    \emph{Survey Methodology} \bold{25}, 43--56.

Gwet, J.-P. and L.-P. Rivest (1992).
    Outlier Resistant Alternatives to the Ratio Estimator,
    \emph{Journal of the American Statistical Association} \bold{87},
    1174--1182.

Lee, H. (1991). Model-Based Estimators That Are Robust to Outliers,
    in \emph{Proceedings of the 1991 Annual Research Conference},
    Bureau of the Census, 178--202. Washington, DC, Department of Commerce.

Lee, H. (1995). Outliers in business surveys,
    in \emph{Business survey methods},
    ed. by B. G. Cox, D. A. Binder, B. N. Chinnappa, A. Christianson,
    M. J. Colledge, and P. S. Kott,
    New York: John Wiley and Sons, chap. 26, 503--526.

Särndal, C.-E., B. Swensson and J. Wretman (1992).
    \emph{Model Assisted Survey Sampling}, New York: Springer.

}
\examples{
data(workplace)

library(survey)
# Survey design for simple random sampling without replacement
dn <- svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight,
    data = workplace)

# Robust regression M-estimator with Huber psi-function
reg <- svyreg_huberM(payroll ~ employment, dn, k = 3)

# Robust GREG estimator of the mean; the population means of the auxiliary
# variables are from a register
m <- svymean_reg(reg, auxiliary = c(1, 11.03), type = "ADU")
m

# ADU
m <- svytotal_reg(reg, auxiliary = c(86514, 953555), type = "ADU")
m

m <- svyreg(payroll ~ employment, dn)
svytotal_reg(m, auxiliary = c(86514, 953555), type = "ADU")


# Summarize
summary(m)

# Extract estimate
coef(m)

# Extract estimated standard error
SE(m)
}
\seealso{
\code{\link{svyreg}}, \code{\link{svyreg_huberM}}, \code{\link{svyreg_huberGM}},
\code{\link{svyreg_tukeyM}}, and \code{\link{svyreg_tukeyGM}}
}
