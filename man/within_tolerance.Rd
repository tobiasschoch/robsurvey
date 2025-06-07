\name{within_tolerance}
\alias{within_tolerance}
\title{Tolerance Interval}
\usage{
within_tolerance(x, w, method = c("quartile", "modified", "boxplot"),
                 constants, lambda = 0.05, info = FALSE,
                 boxplot_coef = 1.5)
}
\arguments{
    \item{x}{\code{[numeric vector]} data vector.}
    \item{w}{\code{[numeric vector]} design weights (same length as \code{x}).}
    \item{method}{[character] one of the methods: \code{"quartile"},
        \code{"modified"} (quartile method), or \code{"boxplot"}.}
    \item{constants}{\code{[numeric vector]} a vector of size two with
        nonnegative tuning constants; it is only used by the methods
        \code{"quartile"} and \code{"modified"}.}
    \item{lambda}{\code{[numeric]} a tuning constant that takes values in the
        closed unit interval; it is only used by method \code{"modified"},
        default: \code{lambda = 0.05}.}
    \item{info}{\code{[logical]} if \code{TRUE}, the tolerance interval is
        printed out.}
    \item{boxplot_coef}{\code{[numeric]} determines how far the whiskers of the
       boxplot extend out from the box; the default is 1.5.}
}
\description{
    The function flags observations that fall within the tolerance interval.
    Observations that fall outside the interval are regarded as (potential)
    outliers.
}
\details{
    Three methods are available.

    \describe{
        \item{Quartile method (\code{"quartile"})}{
            For the quartile method, the tolerance interval is given by
            \deqn{[m - c_l \cdot L_l, \; m + c_u \cdot L_u]}{
                [m - c_l * L_l,  m + c_u * L_u]}
            with
            \deqn{L_l = m - q_1 \quad \text{and} \quad L_u = q_3 - m,}{
                L_l = m - q_1,  and  L_u = q_3 - m,}
            where \eqn{m}{m} denotes the (weighted) median; \eqn{q_1}{q_1} and
            \eqn{q_3}{q_3} are, respectively, the first and third (weighted)
            quartiles. The tuning constants \eqn{c_l}{c_l} and \eqn{c_u}{c_u}
            are combined into the vector \eqn{(c_l, c_u)}{(c_l, c_u)}, which is
            available as argument \code{constants}; both constants must be
            nonnegative numbers.

            The quartiles are calculated using design weights.
        }
        \item{Modified quartile method (\code{"modified"})}{
            For the modified quartile method (Lee, 1995), the tolerance
            interval is given by replacing \eqn{L_l}{L_l} and \eqn{L_u}{L_u}
            with, respectively,
            \deqn{L_l = \max\big(m - q_1, \vert \lambda \cdot m\vert\big),}{
            L_l = max[m - q_1, |lambda * m|],}
            and
            \deqn{L_u = \max\big(q_3 - m, \vert \lambda \cdot m \vert\big)}{
            L_u = max[q_3 - m, |lambda * m|].}
            The tuning constant \eqn{\lambda}{lambda} can only take values in
            the closed unit interval and is available as argument \code{lambda}.

            The quartiles are calculated using design weights.
        }
        \item{Boxplot (box-and-whisker plot) method (\code{"boxplot"})}{
            The tolerance interval for the boxplot method extends from the
            lower whisker to the upper whisker. By default, the length of the
            whiskers is set to 1.5 times the interquartile range; see argument
            \code{boxplot_coef}. For more details, see
            \link[graphics]{boxplot}.

            The quartiles, and therefore the interquartile range, are
            calculated using design weights.
        }
    }
}
\value{
A vector of logicals, where \code{TRUE} indicates that an observation is within
the tolerance limits and \code{FALSE} indicates a (potential) outlier.

If \code{info = TRUE}, the function prints the tolerance interval. The
endpoints of the interval can be numbers or the symbols \sQuote{min.} and
\sQuote{max.}, which denote the minimum and maximum values in the data,
respectively.
}
\seealso{
\link[=robsurvey-package]{Overview} (of all implemented functions)
}
\references{
Lee, H. (1995). Outliers in Business Surveys, in: Cox, B. G. et al. (eds.),
    \emph{Business Survey Methods}, p. 503--526. New York: John Wiley and Sons.
}
\examples{
head(workplace)
attach(workplace)

# Show the tolerance limits
within_tolerance(payroll, weight, method = "boxplot", info = TRUE)

# Observations that fall outside the tolerance limits are (potential) outliers
outlier <- !within_tolerance(payroll, weight, method = "boxplot")
outlier[1:10]
}
