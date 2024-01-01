\name{svymean_ratio}
\alias{svymean_ratio}
\alias{svytotal_ratio}
\title{Robust Ratio Predictor of the Mean and Total}
%FIXME: N and total
\usage{
svytotal_ratio(object, total, variance = "wu", keep_object = TRUE)
svymean_ratio(object, total, N = NULL, variance = "wu",
              keep_object = TRUE, N_unknown = FALSE)
}
\arguments{
    \item{object}{an object of class \code{[ratio]}, e.g., result of
        the Huber ratio \eqn{M}-estimator \code{svyratio_huber}.}
    \item{total}{\code{[numeric]} vector of population totals of the
        auxiliary variables.}
    \item{N}{\code{[numeric]} population size (see also \code{N_unknown}.}
    \item{variance}{\code{[character]} type of variance estimator (default:
        "wu"); see Details Section.}
    \item{keep_object}{\code{[logical]} if \code{TRUE}, \code{object} is
        returned as an additional slot of the return value (default:
        \code{TRUE}).}
    \item{N_unknown}{\code{[logical]} if \code{TRUE}, it is assumed that
        the population size is unknown; thus, it is estimated
        (default: \code{FALSE}).}
}
\description{
    Robust ratio predictor (\emph{M}-estimator) of the population mean and
    total with Huber and Tukey biweight (bisquare) psi-function.
}
\details{
    Package \pkg{survey} must be attached to the search path in order to use
    the functions (see \code{\link{library}} or \code{\link{require}}).

    The (robust) ratio predictor of the population total or mean is
    computed in two steps.
    \itemize{
        \item Step 1: Fit the ratio model associated with the predictor
            by one of the functions \code{\link{svyratio_huber}}
            or \code{\link{svyratio_tukey}}. The fitted model is called
            \code{object}.
        \item Step 2: Based on the fitted model obtained in the first step,
            we predict the population total and mean, respectively, by
            the predictors \code{svytotal_ratio} and \code{svymean_ratio},
            where \code{object} is the fitted ratio model.
    }
    \describe{
        \item{Auxiliary data}{
            Two types of auxiliary variables are distinguished: (1)
            population size \eqn{N} and (2) the population total of the
            auxiliary variable (denominator) used in the ratio model.

            The option \code{N_unknown = TRUE} can be used in the predictor
            of the population mean if \eqn{N} is unknown.
        }
        \item{Variance estimation}{
            Three variance estimators are implemented (argument
            \code{variance}): \code{"base"}, \code{"wu"}, and \code{"hajek"}.
            These estimators correspond to the estimators \code{v0}, \code{v1},
            and \code{v2} in Wu (1982).
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
Wu, C.-F. (1982). Estimation of Variance of the Ratio Estimator.
    \emph{Biometrika} \bold{69}, 183--189.
}
\seealso{
\link[=robsurvey-package]{Overview} (of all implemented functions)

\code{\link{svymean_reg}} and \code{\link{svytotal_reg}} for (robust) GREG
regression predictors

\code{\link{svyreg_huberM}}, \code{\link{svyreg_huberGM}},
\code{\link{svyreg_tukeyM}} and \code{\link{svyreg_tukeyGM}} for robust
regression \eqn{M}- and \eqn{GM}-estimators

\code{\link{svymean_huber}}, \code{\link{svytotal_huber}},
\code{\link{svymean_tukey}} and \code{\link{svytotal_tukey}} for
\eqn{M}-estimators
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

# Robust ratio M-estimator with Huber psi-function
rat <- svyratio_huber(~payroll, ~ employment, dn, k = 5)

# Summary of the ratio estimate
summary(rat)

# Diagnostic plots of the ration/regression M-estimate (e.g.,
# standardized residuals against fitted values)
plot(rat, which = 1L)

# Plot of the robustness weights of the ratio/regression M-estimate
# against its residuals
plot(residuals(rat), robweights(rat))

# Robust ratio predictor of the population mean
m <- svymean_ratio(rat, total = 1001233, N = 90840)
m

# Summary of the ratio estimate of the population mean
summary(m)

# Extract estimate
coef(m)

# Extract estimate of scale
scale(m)

# Extract estimated standard error
SE(m)
}
