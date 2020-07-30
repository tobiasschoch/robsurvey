#' Weighted robust line fitting
#'
#' \code{weighted_line} fits a robust line and allows weights.
#'
#' \code{weighted_line} uses different quantiles for splitting the sample than \code{stats::line()}.
#'
#' @param x \code{[numeric vector]} explanatory variable.
#' @param y \code{[numeric vector]} response variable (default: \code{NULL}).
#' @param w \code{[numeric vector]} weights (same length as vector \code{x}).
#' @param iter \code{[integer]} number of iterations to use (default: \code{1}). 
#' @param na.rm \code{[logical]} indicating whether \code{NA} values should be removed before the computation proceeds (default: \code{FALSE}). 
#' @return intercept and slope of the fitted line
#' @examples
#' data(cars)
#' weighted_line(cars$speed, cars$dist, w = rep(1, length(cars$speed)))
#' weighted_line(cars$speed, cars$dist, w = rep(1:10, each = 5))
#' @seealso \code{\link[stats]{line}}
#' @export weighted_line
weighted_line <- function(x, y = NULL, w, na.rm = FALSE, iter = 1)
{
   if(missing(w)) stop("Argument 'w' (weights) is missing, with no default.\n")

   # quantiles as implemented in line() but with weights x and w sorted 
   # according to x
   line.quantile <- function(x, w, prob)
   {
      n <- length(x)
      cumw <- n * cumsum(w) / sum(w)
      low <- min(which(cumw >= (n - 1) * prob))
      high <- max(which(cumw <= (n - 1) * prob))
      return((x[low + 1] + x[high + 1]) /2)
   }

   if (inherits(x, "formula")) {
   mf <- stats::model.frame(x)
   y <- mf[, 1]
   x <- mf[, -1]
   }
   if (NCOL(x) > 1) stop("Argument 'x' contains more than 1 explanatory 
      variables.")
   dat <- data.frame(x, y, w)
   ok <- stats::complete.cases(dat$x, dat$y, dat$w)
   n <- sum(ok)
   if (n < length(x)) {
      if (na.rm) {
	 x <- dat$x[ok]; y <- dat$y[ok]; w <- dat$w[ok]
      } else {
	 stop("There are missing values in 'x', 'y' or 'w'.\n")
      }
   }

   # Sort data according to x
   ord <- order(x)
   x <- x[ord]
   y <- y[ord]

   # standardise weights to n
   w <- n * w[ord] / sum(w)

   # Groups
   lim1 <- line.quantile(x, w, 1 / 3)
   lim2 <- line.quantile(x, w, 2 / 3)
   groups <- (x <= lim1) * 1 + (x > lim1 & x < lim2) * 2 + (x >= lim2) * 3

   # Medians for x
   wmedx <- c(weighted_median(x[groups == 1], w[groups == 1], na.rm), 
	      weighted_median(x[groups == 3], w[groups == 3], na.rm))

   # polishing (affects only the slope)
   slope <- 0; r <- y; j <- 0
   while (j <= iter - 1) {
      wmedr <- c(weighted_median(r[groups == 1], w[groups == 1], na.rm),
	         weighted_median(r[groups == 3], w[groups == 3], na.rm))
      slope <- slope + (wmedr[2] - wmedr[1]) / (wmedx[2] - wmedx[1])
      r <- y - slope * x
      j <- j + 1
   }

   # intercept and predicted values
   intercept <- weighted_median(r, w, na.rm)
   yhat <- intercept + slope * x

   structure(list(call = match.call(), coefficients = c(intercept, slope),
      residuals = y - yhat, fitted.values = yhat), class = "tukeyline")
}

#' Robust simple linear regression based on medians
#'
#' Robust simple linear regression based on medians: two methods are available. 
#'
#' \describe{
#'    \item{Overview.}{Robust simple linear regression based on medians}
#'    \item{Type.}{Two methods/ types are available (let \eqn{m(x,w)} denote the weighted median of variable \code{x} with weights \code{w}):
#'	 \describe{
#'	    \item{\code{type = "slopes"}:}{The slope is computed as
#'	       \deqn{b1 = m\left( \frac{y - m(y,w)}{x - m(x,w)}, w\right).}{m[(y - m[y, w]) / (x - m[x, w]), w].}
#'	    }
#'	    \item{\code{type = "products"}:}{The slope is computed as
#'	       \deqn{b1 = \frac{m\big([y - m(y,w)][x - m(x,w)], w\big)}{m\big([x - m(x,w)]^2, w\big)}.}{m([y - m(y, w)][x - m(x, w)], w) / m([x - m(x, w)]^2, w).} 
#'	    }
#'	 }
#'    }
#' }
#'
#' @param x \code{[numeric vector]} explanatory variable.
#' @param y \code{[numeric vector]} response variable (default: \code{NULL}).
#' @param w \code{[numeric vector]} weights (same length as vector \code{x}).
#' @param type \code{[character]} either \code{slopes} or \code{products} (default: \code{slopes}). 
#' @param na.rm \code{[logical]} indicating whether \code{NA} values should be removed before the computation proceeds (default: \code{FALSE}). 
#' @return a vector with two components: intercept and slope
#'
#' @examples
#' x <- c(1, 2, 4, 5)
#' y <- c(3, 2, 7, 4)
#' weighted_line(y~x, w=rep(1, length(x)))
#' weighted_median_line(y~x, w = rep(1, length(x)))
#' weighted_median_line(y~x, w = rep(1, length(x)), type = "prod")
#'
#' data(cars)
#' with(cars, weighted_median_line(dist ~ speed, w = rep(1, length(dist))))
#' with(cars, weighted_median_line(dist ~ speed, w = rep(1, length(dist)), 
#' type = "prod"))
#'
#' # weighted
#' w <- c(rep(1,20), rep(2,20), rep(5, 10))
#' with(cars, weighted_median_line(dist ~ speed, w = w))
#' with(cars, weighted_median_line(dist ~ speed, w = w, type = "prod"))
#'
#' # outlier in y
#' cars$dist[49] <- 360
#' with(cars, weighted_median_line(dist ~ speed, w = w))
#' with(cars, weighted_median_line(dist ~ speed, w = w, type = "prod"))
#'
#' # outlier in x
#' data(cars)
#' cars$speed[49] <- 72
#' with(cars, weighted_median_line(dist ~ speed, w = w))
#' with(cars, weighted_median_line(dist ~ speed, w = w, type = "prod"))
#' @seealso \code{\link[stats]{line}}, \code{\link{weighted_line}}, \code{\link{weighted_median_ratio}}
#' @export weighted_median_line
weighted_median_line <- function(x, y = NULL, w, type = "slopes", 
   na.rm = FALSE)
{
   if(missing(w)) stop("Argument 'w' (weights) is missing, with no default.\n")

   if (inherits(x, "formula")) {
      dat <- grDevices::xy.coords(x)
      x <- dat$x
      y <- dat$y
   }
   if (NCOL(x) > 1) stop("Argument 'x' contains more than 1 explanatory 
      variables.\n")

   # Robustification type
   stype <- pmatch(type, c("slopes"), nomatch = 2)

   # Ensure case-wise completeness
   ok <- stats::complete.cases(data.frame(x, y, w))
   n <- sum(ok)
   if (n < length(x)) {
      if (na.rm) {
	 x <- x[ok]; y <- y[ok]; w <- w[ok]
      } else {
	 stop("There are missing values in 'x', 'y' or 'w'.\n")
      }
   }

   # univariate medians
   wmedx <- weighted_median(x, w, na.rm)
   wmedy <- weighted_median(y, w, na.rm)

   # slope (remove NA created due to division by 0)
   slope <- switch(stype,
      weighted_median((y - wmedy) / (x - wmedx), w, na.rm = TRUE),
      weighted_median((y - wmedy) * (x - wmedx), w, na.rm = TRUE) / 
	 weighted_median((x - wmedx)^2, w, na.rm = TRUE)
   )

   # residuals
   r <- y - slope * x

   # intercept
   intercept <- weighted_median(r, w, na.rm)
   yhat <- intercept + slope * x

   structure(list(call = sys.call(), coefficients = c(intercept, slope),
	      residuals = y - yhat, fitted.values = yhat), class = "medline")
}

#' Weighted robust ratio based on median
#'
#' A weighted median of the ratios y/x determines the slope of a regression through the origin.
#'
#' @param x \code{[numeric vector]} explanatory variable.
#' @param y \code{[numeric vector]} response variable (default: \code{NULL}).
#' @param w \code{[numeric vector]} weights (same length as vector \code{x}).
#' @param na.rm \code{[logical]} indicating whether \code{NA} values should be removed before the computation proceeds (default: \code{FALSE}). 
#' @return a vector with two components: intercept and slope
#' @examples
#' x <- c(1,2,4,5)
#' y <- c(1,0,5,2)
#' weighted_median_ratio(y~x, w = rep(1, length(y)))
#' @seealso \code{\link[stats]{line}}, \code{\link{weighted_line}}, \code{\link{weighted_median_line}}
#' @export weighted_median_ratio
weighted_median_ratio <- function(x, y = NULL, w, na.rm = FALSE)
{
   if(missing(w)) stop("Argument 'w' (weights) is missing, with no default.\n")

   if (inherits(x, "formula")) {
      dat <- grDevices::xy.coords(x)
      y <- dat$y
      x <- dat$x
   }
   if (NCOL(x) > 1) stop("Argument 'x' contains more than 1 explanatory 
   variables.\n")

   # Ensure case-wise completeness
   ok <- stats::complete.cases(data.frame(x, y, w))
   n <- sum(ok)
   if (n < length(x)) {
      if (na.rm) {
	 x <- x[ok]; y <- y[ok]; w <- w[ok]
      } else {
	 stop("There are missing values in 'x', 'y' or 'w'.\n")
      }
   }

   # ratio
   ratio <- weighted_median(y / x, w, na.rm)

   # fitted.varlues and residuals
   yhat <- ratio * x
   r <- y - yhat

   structure(list(call = match.call(), coefficients = ratio,
	   residuals = y - yhat, fitted.values = yhat), class = "medline")
} 
