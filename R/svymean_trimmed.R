#' Weighted trimmed mean and total 
#'
#' Weighted trimmed mean and total 
#'
#' \describe{
#'    \item{Characteristic.}{Population mean or total. Let \eqn{\mu} denote the estimated trimmed population mean; then, the estimated trimmed total is given by \eqn{N \mu} with \eqn{N =\sum w_i}, where summation is over all observations in the sample.}
#'    \item{Trimming.}{The methods trims the \code{LB}\eqn{~\cdot 100\%} percentage of the smallest observations and the (1 - \code{UB})\eqn{~\cdot 100\%} percentage of the largest observations from the data.} 
#'    \item{Variance estimation.}{Taylor linearization.}
#'    \item{Utility functions.}{\code{\link[=svystat_rob]{summary}}, \code{\link[=svystat_rob]{coef}}, \code{\link[=svystat_rob]{SE}}, \code{\link[=svystat_rob]{vcov}}, \code{\link[=svystat_rob]{residuals}}, \code{\link[=svystat_rob]{fitted}}, and \code{\link[=svystat_rob]{robweights}}.}
#'    \item{Bare-bone functions.}{See \code{\link{weighted_mean_trimmed}} and \code{\link{weighted_total_trimmed}}.}
#' } 
#'
#' @param x a one-sided \code{[formula]}, e.g., \code{~myVariable}. 
#' @param design an object of class \code{survey.design} or \code{survey.design2}.
#' @param LB \code{[double]} lower bound of trimming, such that \eqn{0 \leq} \code{LB} \eqn{<} \code{UB} \eqn{\leq 1}.
#' @param UB \code{[double]} upper bound of trimming, such that \eqn{0 \leq} \code{LB} \eqn{<} \code{UB} \eqn{\leq 1}.
#' @param na.rm \code{[logical]} indicating whether \code{NA} values should be removed before the computation proceeds (default: \code{FALSE}). 
#' @return object of class \code{\link{svystat_rob}} 
#'
#' @examples
#' data(workplace) 
#'
#' library(survey)
#' # Survey design for simple random sampling without replacement
#' dn <- svydesign(ids = ~ID, strata = ~strat, fpc = ~fpc, weights = ~weight, 
#'                 data = workplace)
#'  
#' # Estimated trimmed population total (5% symmetric trimming) 
#' svytotal_trimmed(~employment, dn, LB = 0.05, UB = 0.95) 
#' 
#' # Estimated trimmed population mean (5% trimming at the top of the distr.) 
#' svymean_trimmed(~employment, dn, UB = 0.95) 
#' @seealso \code{\link{weighted_mean_trimmed}} and \code{\link{weighted_total_trimmed}}
#' @export
svymean_trimmed <- function(x, design, LB = 0.05, UB = 1 - LB, na.rm = FALSE)
{
   dat <- .checkformula(x, design)
   res <- weighted_mean_trimmed(dat$y, dat$w, LB, UB, info = TRUE, na.rm)
   # influence function 
   infl <- .infl_trimmed(dat$y, dat$w, LB, UB, res$estimate)
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

#' @rdname svymean_trimmed
#' @export
svytotal_trimmed <- function(x, design, LB = 0.05, UB = 1 - LB, na.rm = FALSE)
{
   res <- svymean_trimmed(x, design, LB, UB, na.rm)
   sum_w <- sum(res$model$w) 
   res$estimate <- res$estimate * sum_w 
   res$variance <- res$variance * sum_w^2 
   res$characteristic <- "total"
   res$call <- match.call()
   res
} 
