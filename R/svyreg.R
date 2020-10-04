#' Survey regression estimator 
#' 
#' Weighted regression estimator
#'  
#' Details
#'
#' @param formula a \code{[formula]} object (i.e., symbolic description of the model)
#' @param design an object of class \code{survey.design} or \code{survey.design2}.
#' @param var \code{[numeric vector]} heteroscedastic variance (default: \code{NULL}, i.e., homoscedastic variance).
#' @param na.rm \code{[logical]} indicating whether \code{NA} values should be removed before the computation proceeds (default: \code{FALSE}). 
#' @return object of class \code{svyreg} 
#' @export 
svyreg <- function(formula, design, var = NULL, na.rm = FALSE)
{
   res <- svyreg_huber(formula, design, var, k = Inf, na.rm) 
   res$estimator <- "Survey regression estimator"
   res$call <- match.call()
   res$robust <- NULL 
   res$optim <- NULL
   class(res) <- "svyreg_rob"
   res
} 
