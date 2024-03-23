\name{svymean_reg}
\alias{svymean_reg}
\alias{svytotal_reg}
\title{Robust Generalized Regression Predictor (GREG) of the Mean and Total}
\usage{
svytotal_reg(object, totals, N = NULL, type, k = NULL, check.names = TRUE,
             keep_object = TRUE, ...)
svymean_reg(object, totals, N = NULL, type, k = NULL, check.names = TRUE,
            keep_object = TRUE, N_unknown = FALSE, ...)
}
\arguments{
    \item{object}{an object of class \code{[svyreg_rob]}, e.g., result of
        the Huber regression \eqn{M}-estimator \code{svyreg_huberM}.}
    \item{totals}{\code{[numeric]} vector of population totals of the
        auxiliary variables.}
    \item{N}{\code{[numeric]} population size (see also \code{N_unknown}.}
    \item{type}{\code{[character]} type of predictor; see Details Section.}
    \item{k}{\code{[numeric]} robustness tuning constant of the
        \code{psi}-function used in the bias-correction term of the GREG.
        The definition of \code{k} depends on the \code{type} of predictor
        and is discussed in the Details Section.}
    \item{check.names}{\code{[logical]} if \code{TRUE}, the names of
        \code{auxiliary} are checked against the names of the independent
        variables of the fitted model \code{object} (default: \code{TRUE}).}
    \item{keep_object}{\code{[logical]} if \code{TRUE}, \code{object} is
        returned as an additional slot of the return value (default:
        \code{TRUE}).}
    \item{N_unknown}{\code{[logical]} if \code{TRUE}, it is assumed that
        the population size is unknown; thus, it is estimated
        (default: \code{FALSE}).}
    \item{\dots}{additional arguments (currently not used).}
}
\description{
    Generalized regression estimator (GREG) predictor of the mean and total,
    and robust GREG \eqn{M}-estimator predictor
}
\details{
    Package \pkg{survey} must be attached to the search path in order to use
    the functions (see \code{\link{library}} or \code{\link{require}}).

    The (robust) GREG predictor of the population total or mean is
    computed in two steps.
    \itemize{
        \item Step 1: Fit the regression model associated with the GREG
            predictor by one of the functions \code{\link{svyreg}},
            \code{\link{svyreg_huberM}}, \code{\link{svyreg_huberGM}},
            \code{\link{svyreg_tukeyM}} or \code{\link{svyreg_tukeyGM}}.
            The fitted model is called \code{object}.
        \item Step 2: Based on the fitted model obtained in the first step,
            we predict the population total and mean, respectively, by
            the predictors \code{svytotal_reg} and \code{svymean_reg},
            where \code{object} is the fitted regression model.
    }
    The following GREG predictors are available:
    \describe{
        \item{GREG (not robust, \code{k = NULL})}{
            The following \emph{non-robust} GREG predictors are available:
            \itemize{
                \item \code{type = "projective"} ignores the bias correction
                    term of the GREG predictor; see Särndal and Wright (1984).
                \item \code{type = "ADU"} is the "standard" GREG,
                    which is an asymptotically design unbiased (ADU)
                    predictor; see Särndal et al.(1992, Chapter 6).
            }
            If the fitted regression model (\code{object}) does include
            a regression intercept, the predictor types \code{"projective"}
            and \code{"ADU"} are identical because the bias correction
            of the GREG is zero by design.
        }
        \item{Robust GREG}{
            The following \emph{robust} GREG predictors are available:
            \itemize{
                \item \code{type = "huber"} and \code{type = "tukey"} are,
                    respectively, the robust GREG predictors with Huber
                    and Tukey bisquare (biweight) psi-function. The tuning
                    constant must satisfy \code{0 < k <= Inf}.
                    We can use the Huber-type GREG predictor although the
                    model has been fitted by the regression estimator
                    with Tukey psi-function (and vice versa).
                \item \code{type = "BR"} is the bias-corrected robust GREG
                    predictor of Beaumont and Rivest (2009), which is
                    inspired by the bias-corrected robust predictor of
                    Chambers (1986). The tuning constant must satisfy
                    \code{0 < k <= Inf}.
                \item \code{type = "lee"} is the bias-corrected predictor
                    of Lee (1991; 1992). Tthe tuning constant \code{k} must
                    satisfy \code{0 <= k <= 1}.
                \item \code{type = "duchesne"} is the bias-corrected,
                    calibration-type estimator/ predictor of Duchesne (1999).
                    The tuning constant \code{k} must be specified as a
                    vector \code{k = c(a, b)}, where \code{a} and \code{b}
                    are the tuning constants of Duchesne's modified Huber
                    psi-function (default values: \code{a = 9} and
                    \code{b = 0.25}).
            }
        }
        \item{Auxiliary data}{
            Two types of auxiliary variables are distinguished: (1)
            population size \eqn{N} and (2) population totals of the
            auxiliary variables used in the regression model (i.e.,
            non-constant explanatory variables).

            The option \code{N_unknown = TRUE} can be used in the predictor
            of the population mean if \eqn{N} is unknown.

            The names of the entries of \code{totals} are checked against
            the names of the regression fit (\code{object}), unless we specify
            \code{check.names = FALSE}.
        }
        \item{Utility functions}{
            The return value is an object of class \code{\link{svystat_rob}}.
            Thus, the utility functions \code{\link[=svystat_rob]{summary}},
                \code{\link[=svystat_rob]{coef}},
                \code{\link[=svystat_rob]{SE}},
                \code{\link[=svystat_rob]{vcov}},
                \code{\link[=svystat_rob]{residuals}},
                \code{\link[=svystat_rob]{fitted}}, and
                \code{\link[=svystat_rob]{robweights}} are available.
        }
    }
}
\value{
    Object of class \code{\link{svystat_rob}}
}
\references{
Beaumont, J.-F. and Rivest, L.-P. (2009). Dealing with outliers in survey
    data. In: \emph{Sample Surveys: Theory, Methods and Inference}
    ed. by Pfeffermann, D. and Rao, C. R. Volume 29A of
    \emph{Handbook of Statistics}, Amsterdam: Elsevier, Chap. 11, 247--280.
    \doi{10.1016/S0169-7161(08)00011-4}

Chambers, R. (1986). Outlier Robust Finite Population Estimation.
    \emph{Journal of the American Statistical Association} \bold{81},
    1063--1069. \doi{10.1080/01621459.1986.10478374}

Duchesne, P. (1999). Robust calibration estimators, \emph{Survey Methodology}
    \bold{25}, 43--56.

Gwet, J.-P. and Rivest, L.-P. (1992). Outlier Resistant Alternatives to
    the Ratio Estimator. \emph{Journal of the American Statistical Association}
    \bold{87}, 1174--1182. \doi{10.1080/01621459.1992.10476275}

Lee, H. (1991). Model-Based Estimators That Are Robust to Outliers,
    in \emph{Proceedings of the 1991 Annual Research Conference},
    Bureau of the Census, 178--202. Washington, DC, Department of Commerce.

Lee, H. (1995). Outliers in business surveys. In:
    \emph{Business survey methods} ed. by Cox, B. G., Binder, D. A.,
    Chinnappa, B. N., Christianson, A., Colledge, M. J. and Kott, P. S.
    New York: John Wiley and Sons, Chap. 26, 503--526.
    \doi{10.1002/9781118150504.ch26}

Särndal, C.-E., Swensson, B. and Wretman, J. (1992).
    \emph{Model Assisted Survey Sampling}, New York: Springer.

Särndal, C.-E. and Wright, R. L. (1984). Cosmetic Form of Estimators
    in Survey Sampling. \emph{Scandinavian Journal of Statistics} \bold{11},
    146--156.
}
\seealso{
\link[=robsurvey-package]{Overview} (of all implemented functions)

\code{\link{svymean_ratio}} and \code{\link{svytotal_ratio}} for (robust)
ratio predictors

\code{\link{svymean_huber}}, \code{\link{svytotal_huber}},
\code{\link{svymean_tukey}} and \code{\link{svytotal_tukey}} for
\eqn{M}-estimators

\code{\link{svyreg}}, \code{\link{svyreg_huberM}}, \code{\link{svyreg_huberGM}},
\code{\link{svyreg_tukeyM}} and \code{\link{svyreg_tukeyGM}} for robust
regression \eqn{M}- and \eqn{GM}-estimators
}
\examples{
head(workplace)

library(survey)
# Survey design for stratified simple random sampling without replacement
dn <- if (packageVersion("survey") >= "4.2") {
        # survey design with pre-calibrated weights
        svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight,
                  data = workplace, calibrate.formula = ~-1 + strat)
    } else {
        # legacy mode
        svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight,
                  data = workplace)
    }

# Robust regression M-estimator with Huber psi-function
reg <- svyreg_huberM(payroll ~ employment, dn, k = 3)

# Summary of the regression M-estimate
summary(reg)

# Diagnostic plots of the regression M-estimate (e.g., standardized
# residuals against fitted values)
plot(reg, which = 1L)

# Plot of the robustness weights of the regression M-estimate against
# its residuals
plot(residuals(reg), robweights(reg))

# ADU (asymptotically design unbiased) estimator
m <- svytotal_reg(reg, totals = 1001233, 90840, type = "ADU")
m

# Robust GREG estimator of the mean; the population means of the auxiliary
# variables are from a register
m <- svymean_reg(reg, totals = 1001233, 90840, type = "huber", k = 20)
m

# Summary of the robust GREG estimate
summary(m)

# Extract estimate
coef(m)

# Extract estimated standard error
SE(m)

# Approximation of the estimated mean square error
mse(m)
}
