#' Weighted winsorized mean and total 
#'
#' Weighted winsorized mean and total 
#'
#' \describe{
#'    \item{Characteristic.}{Population mean or total. Let \eqn{\mu} denote the estimated winsorized population mean; then, the estimated winsorized total is given by \eqn{N \mu} with \eqn{N =\sum w_i}, where summation is over all observations in the sample.}
#'    \item{Modes of winsorization.}{The amount of winsorization can be specified in relative or absolute terms:
#'	 \itemize{
#'	    \item \emph{relative:} By specifying \code{LB} and \code{UB}, the method winsorizes the \code{LB}\eqn{~\cdot 100\%} percentage of the smallest observations and the (1 - \code{UB})\eqn{~\cdot 100\%} percentage of the largest observations from the data. 
#'	    \item \emph{absolute:} By specifying argument \code{k} in the functions with the "infix" \code{_k_} in their name, the largest \eqn{k} observations are winsorized, \eqn{0<k<n}, where \eqn{n} denotes the sample size.
#'	 }
#'    } 
#'    \item{Variance estimation.}{Taylor linearization; two estimators are available: 
#'	 \describe{
#'	    \item{\code{simple_var = FALSE}:}{Variance estimator of the winsorized mean/ total. The estimator depends on the estimated probability density function evaluated at the winsorization thresholds, which can be -- depending on the context -- numerically unstable. As a remedy, a simplified variance estimator is available by setting \code{simple_var = TRUE}.}
#'	    \item{\code{simple_var = TRUE}:}{Variance is approximated using the variance estimator of the trimmed mean/ total.}
#'	 }
#'    }
#'    \item{Utility functions.}{\code{\link[=svystat_rob]{summary}}, \code{\link[=svystat_rob]{coef}}, \code{\link[=svystat_rob]{SE}}, \code{\link[=svystat_rob]{vcov}}, \code{\link[=svystat_rob]{residuals}}, \code{\link[=svystat_rob]{fitted}}, and \code{\link[=svystat_rob]{robweights}}.}
#'    \item{Bare-bone functions.}{See: 
#'	 \itemize{
#'	    \item \code{\link{weighted_mean_winsorized}}, 
#'	    \item \code{\link{weighted_mean_k_winsorized}}, 
#'	    \item \code{\link{weighted_total_winsorized}}, 
#'	    \item \code{\link{weighted_total_k_winsorized}}.
#'	 }
#'    }
#' } 
#'
#' @param x a one-sided \code{[formula]}, e.g., \code{~myVariable}. 
#' @param design an object of class \code{survey.design} or \code{survey.design2}.
#' @param LB \code{[double]} lower bound of winsorization, such that \eqn{0 \leq} \code{LB} \eqn{<} \code{UB} \eqn{\leq 1}.
#' @param UB \code{[double]} upper bound of winsorization, such that \eqn{0 \leq} \code{LB} \eqn{<} \code{UB} \eqn{\leq 1}.
#' @param na.rm \code{[logical]} indicating whether \code{NA} values should be removed before the computation proceeds (default: \code{FALSE}). 
#' @param simple_var \code{[logical]} indicating whether a simplified variance estimator should be used (default: \code{FALSE}).
#' @param k \code{[integer]} number of observations to be winsorized at the top of the distribution. 
#' @return object of class \code{\link{svystat_rob}} 
#' @examples
#' data(workplace) 
#'
#' library(survey)
#' # Survey design for simple random sampling without replacement
#' dn <- svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight, 
#'                 data = workplace)
#'  
#' # Estimated winsorized population mean (5% symmetric winsorization)
#' svymean_winsorized(~employment, dn, LB = 0.05)
#'
#' # Estimated one-sided k winsorized population total (2 observations are 
#' # winsorized at the top of the distribution)
#' svytotal_k_winsorized(~employment, dn, k = 2)
#' @seealso \code{\link{weighted_mean_winsorized}}, \code{\link{weighted_mean_k_winsorized}}, \code{\link{weighted_total_winsorized}}, and \code{\link{weighted_total_k_winsorized}} 
#' @export
svymean_winsorized <- function(x, design, LB = 0.05, UB = 1 - LB, na.rm = FALSE,
   simple_var = FALSE)
{
   dat <- .checkformula(x, design)
   res <- weighted_mean_winsorized(dat$y, dat$w, LB, UB, info = TRUE, na.rm)
   # influence function 
   if (simple_var){
      infl <- .infl_trimmed(dat$y, dat$w, LB, UB, res$estimate)
   }else{
      infl <- .infl_winsorized(dat$y, dat$w, LB, UB, res$estimate)
   }
   # variance 
   infl <- infl * dat$w / sum(dat$w) 
   res$variance <- survey::svyrecvar(infl, design$cluster, design$strata, 
        design$fpc, postStrata = design$postStrata)
   names(res$estimate) <- dat$yname
   res$call <- match.call()
   res$design <- design
   class(res) <- "svystat_rob"
   res
}

#' @rdname svymean_winsorized
#' @export
svymean_k_winsorized <- function(x, design, k, na.rm = FALSE, 
   simple_var = FALSE)
{
   dat <- .checkformula(x, design)
   res <- weighted_mean_k_winsorized(dat$y, dat$w, k, info = TRUE, na.rm)
   w <- res$model$w; x <- res$model$y
   # influence function 
   if (simple_var){
      infl <- .infl_trimmed(dat$y, dat$w, 0, res$robust$UB, res$estimate)
   }else{
      infl <- .infl_winsorized(dat$y, dat$w, 0, res$robust$UB, res$estimate)
   }
   # variance 
   infl <- infl * w / sum(w) 
   res$variance <- survey::svyrecvar(infl, design$cluster, design$strata, 
        design$fpc, postStrata = design$postStrata)
   names(res$estimate) <- dat$yname
   res$call <- match.call()
   res$design <- design
   class(res) <- "svystat_rob"
   res
}

#' @rdname svymean_winsorized
#' @export
svytotal_winsorized <- function(x, design, LB = 0.05, UB = 1 - LB, 
   na.rm = FALSE, simple_var = FALSE)
{
   res <- svymean_winsorized(x, design, LB, UB, na.rm, simple_var)
   sum_w <- sum(res$model$w) 
   res$estimate <- res$estimate * sum_w 
   res$variance <- res$variance * sum_w^2 
   res$characteristic <- "total"
   res$call <- match.call()
   res
}

#' @rdname svymean_winsorized
#' @export
svytotal_k_winsorized <- function(x, design, k, na.rm = FALSE, 
   simple_var = FALSE)
{
   res <- svymean_k_winsorized(x, design, k, na.rm, simple_var)
   sum_w <- sum(res$model$w) 
   res$estimate <- res$estimate * sum_w 
   res$variance <- res$variance * sum_w^2 
   res$characteristic <- "total"
   res$call <- match.call()
   res
} 
