\name{class_svyreg_rob}
\alias{class_svyreg_rob}
\alias{svyreg_rob}
\alias{print.svyreg_rob}
\alias{summary.svyreg_rob}
\alias{coef.svyreg_rob}
\alias{vcov.svyreg_rob}
\alias{residuals.svyreg_rob}
\alias{fitted.svyreg_rob}
\alias{robweights.svyreg_rob}
\alias{plot.svyreg_rob}
\title{Utility Functions for Objects of Class svyreg_rob}
\usage{
\method{print}{svyreg_rob}(x, digits = max(3L, getOption("digits") - 3L), ...)

\method{summary}{svyreg_rob}(object, mode = c("design", "model", "compound"),
    digits = max(3L, getOption("digits") - 3L), ...)

\method{coef}{svyreg_rob}(object, ...)

\method{vcov}{svyreg_rob}(object, mode = c("design", "model", "compound"), ...)

\method{residuals}{svyreg_rob}(object, ...)

\method{fitted}{svyreg_rob}(object, ...)

\method{robweights}{svyreg_rob}(object)

\method{plot}{svyreg_rob}(x, which = 1L:4L,
    hex = FALSE, caption = c("Standardized residuals vs. Fitted Values",
    "Normal Q-Q", "Response vs. Fitted values",
    "Sqrt of abs(Residuals) vs. Fitted Values"),
	panel = if (add.smooth) function(x, y, ...) panel.smooth(x, y,
    iter = iter.smooth, ...) else points, sub.caption = NULL, main = "",
	ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...,
	id.n = 3, labels.id = names(residuals(x)), cex.id = 0.75, qqline = TRUE,
	add.smooth = getOption("add.smooth"), iter.smooth = 3,
	label.pos = c(4, 2), cex.caption = 1, cex.oma.main = 1.25)
}
\arguments{
    \item{x}{object of class \code{svyreg_rob}.}
    \item{digits}{\code{[integer]} minimal number of significant digits.}
    \item{\dots}{additional arguments passed to the method.}
    \item{object}{object of class \code{svyreg_rob}.}
    \item{mode}{mode of variance estimator: \code{"design"}, \code{"model"} or
        \code{"compound"} (default: \code{"design"}).}
    \item{which}{indicating which plots to be drawn; if a subset of the plots
        is required, you can specify a subset of the numbers \code{1:4}.}
    \item{hex}{\code{[logical]} if \code{TRUE}, a hexagonally binned plot
        is shown in place of a scatterplot.}
    \item{caption}{captions to appear above the plots;
		\code{[character]} vector of valid graphics annotations. It
		can be set to \code{""} or \code{NA} to suppress all captions.}
	\item{panel}{panel function. The useful alternative to
		\code{\link{points}}, \code{\link{panel.smooth}} can be chosen
		by \code{add.smooth = TRUE}.}
	\item{sub.caption}{common title \code{[character]}---above the figures
		if there are more than one; used as \code{sub} (s.\code{\link{title}})
		otherwise. If \code{NULL}, as by default, a possible abbreviated
		version of \code{deparse(x$call)} is used.}
	\item{main}{title to each plot \code{[character]}---in addition
		to \code{caption}.}
	\item{ask}{\code{[logical]}; if \code{TRUE}, the user is \emph{ask}ed
        before each plot, see \code{\link{par}(ask=.)}.}
	\item{id.n}{number of points to be labelled in each plot, starting
		with the most extreme, \code{[integer]}.}
	\item{labels.id}{vector of labels \code{[character]}, from which the
		labels for extreme points will be chosen. \code{NULL} uses
		observation numbers.}
	\item{cex.id}{magnification of point labels, \code{[numeric]}.}
	\item{qqline}{\code{[logical]} indicating if a \code{\link{qqline}()}
		should be added to the normal Q-Q plot.}
	\item{add.smooth}{\code{[logical]} indicating if a smoother should be
		added to most plots; see also \code{panel} above.}
	\item{iter.smooth}{the number of robustness iterations \code{[integer]},
		the argument \code{iter} in \code{\link{panel.smooth}()}.}
	\item{label.pos}{positioning of labels \code{[numeric]}, for the left
		half and right half of the graph respectively.}
	\item{cex.caption}{controls the size of \code{caption}, \code{[numeric]}.}
	\item{cex.oma.main}{controls the size of the \code{sub.caption} only if
		that is \emph{above} the figures when there is more than one,
		\code{[numeric]}.}
}
\description{
Methods and utility functions for objects of class \code{svyreg_rob}.

}
\details{
    \describe{
        \item{Variance}{
            For variance estimation (\code{summary} and \code{vcov}),
            three modes are available:
            \itemize{
                \item \code{"design"}: design-based variance estimator using
                    linearization; see Binder (1983)
                \item \code{"model"}: model-based weighted variance estimator
                    (the sampling design is ignored)
                \item \code{"compound"}: design-model-based variance
                    estimator; see Godambe and Thompson (2009)
            }
        }
        \item{Utility functions}{
            \itemize{
                \item \code{summary} gives a summary of the estimation
                    properties
                \item \code{plot} shows diagnostic plots for the estimated
                    regression model
                \item \code{robweights} extracts the robustness weights
                    (if available)
                \item \code{coef} extracts the estimated regression coefficients
                \item \code{vcov} extracts the (estimated) covariance matrix
                \item \code{residuals} extracts the residuals
                \item \code{fitted} extracts the fitted values
            }
        }
    }
}
\references{
Binder, D. A. (1983).  On the Variances of Asymptotically Normal Estimators
from Complex Surveys. \emph{International Statistical Review}, \bold{51},
279--292.

Godambe, V.P. and Thompson, M.E. (2009). Estimating Functions and Survey
Sampling, in: D. Pfeffermann and C.R. Rao (eds.),
\emph{Handbook of Statistics}, \bold{vol. 29B}, Sample Surveys: Inference and
Analysis, Chapter 26, 83--101, Amsterdam: Elsevier.
}
\seealso{
Weighted least squares: \code{\link{svyreg}}; robust weighted regression
\code{\link{svyreg_huber}}, \code{\link{svyreg_huberGM}},
\code{\link{svyreg_tukey}} and \code{\link{svyreg_tukeyGM}}
}
\examples{
data(workplace)

library(survey)
# Survey design for simple random sampling without replacement
dn <- svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight,
    data = workplace)

# Compute the regression estimate
m <- svyreg_huber(payroll ~ employment, dn, k = 8)

# utility functions
summary(m)
coef(m)
SE(m)
vcov(m)
residuals(m)
fitted(m)
}
