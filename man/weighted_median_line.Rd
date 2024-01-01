\name{weighted_median_line}
\alias{weighted_median_line}
\title{Robust Simple Linear Regression Based on Medians}
\usage{
weighted_median_line(x, y = NULL, w, type = "slopes", na.rm = FALSE)
}
\arguments{
    \item{x}{\code{[numeric vector]} explanatory variable.}
    \item{y}{\code{[numeric vector]} response variable (default: \code{NULL}).}
    \item{w}{\code{[numeric vector]} weights (same length as \code{x}).}
    \item{type}{\code{[character]} \code{"slopes"} or \code{"products"}
        (default: \code{"slopes"}).}
    \item{na.rm}{\code{[logical]} indicating whether \code{NA} values should
        be removed before the computation proceeds (default: \code{FALSE}).}
}
\description{
    Robust simple linear regression based on medians: two methods are available:
    \code{"slopes"} and \code{"product"}.
}
\details{
    \describe{
        \item{Overview.}{Robust simple linear regression based on medians}
        \item{Type.}{Two methods/ types are available. Let \eqn{m(x,w)}
            denote the weighted median of variable \code{x} with weights
            \code{w}:
            \describe{
                \item{\code{type = "slopes"}:}{The slope is computed as
                    \deqn{b1 = m\left( \frac{y - m(y,w)}{x - m(x,w)}, w\right).
                        }{m[(y - m[y, w]) / (x - m[x, w]), w].}
                }
                \item{\code{type = "products"}:}{The slope is computed as
                    \deqn{b1 = \frac{m\big([y - m(y,w)][x - m(x,w)], w\big)}
                        {m\big([x - m(x,w)]^2, w\big)}.}{%
                        m([y - m(y, w)][x - m(x, w)], w) / m([x - m(x, w)]^2, w).}
                }
            }
        }
    }
}
\value{
    A vector with two components: intercept and slope
}
\seealso{
\link[=robsurvey-package]{Overview} (of all implemented functions)

\code{\link[stats]{line}}, \code{\link{weighted_line}} and
\code{\link{weighted_median_ratio}}
}
\examples{
x <- c(1, 2, 4, 5)
y <- c(3, 2, 7, 4)
weighted_line(y ~ x, w = rep(1, length(x)))
weighted_median_line(y ~ x, w = rep(1, length(x)))
m <- weighted_median_line(y ~ x, w = rep(1, length(x)), type = "prod")
m
coef(m)
fitted(m)
residuals(m)

# cars data
head(cars)
with(cars, weighted_median_line(dist ~ speed, w = rep(1, length(dist))))
with(cars, weighted_median_line(dist ~ speed, w = rep(1, length(dist)),
type = "prod"))

# weighted
w <- c(rep(1,20), rep(2,20), rep(5, 10))
with(cars, weighted_median_line(dist ~ speed, w = w))
with(cars, weighted_median_line(dist ~ speed, w = w, type = "prod"))

# outlier in y
cars$dist[49] <- 360
with(cars, weighted_median_line(dist ~ speed, w = w))
with(cars, weighted_median_line(dist ~ speed, w = w, type = "prod"))

# outlier in x
data(cars)
cars$speed[49] <- 72
with(cars, weighted_median_line(dist ~ speed, w = w))
with(cars, weighted_median_line(dist ~ speed, w = w, type = "prod"))
}
