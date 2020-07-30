#' Robust regression predictor of the mean 
#' 
#' regression estimator of the mean
#'
#' Details
#'  
#' @section Failure of convergence: 
#' By default, the method assumes a maximum number of \code{maxit = 100} iterations and a numerical tolerance criterion to stop the iterations of \code{tol = 1e-05}. You can run the code with specifications other than the default values by specifying the arguments \code{maxit} and/or \code{tol} in the function call; see also \code{\link{svyreg_control}}. 
#'
#' @param object fitted regression model (object of class \code{svyreg}). 
#' @param mean_auxiliary \code{[numeric vector]} population means of the auxiliary variables. 
#' @return object of class \code{svystat.rob} 
#' @param k \code{[double]} robustness tuning constant (\eqn{0 < k \leq \infty}{0 < k <= Inf}). 
#' @export 
svymean_reg_huber <- function(object, mean_auxiliary, k)
{
# FIXME: Do we need the S3 formalism?
   stopifnot(k > 0)
   call <- match.call()
   call[[1]] <- substitute(svymean_reg_huber)
   # check dimensions
   if (length(mean_auxiliary) != object$model$p){
      stop("Dimension of argument 'mean_auxiliary' is not correct\n")
   }
   # robust greg estimate
   w <- object$model$w; sum_w <- sum(w)
   resid_winsorized <- pmin.int(k, pmax.int(-k, object$residuals))
   est <- sum(mean_auxiliary * object$estimate) + sum(w * resid_winsorized) / 
      sum_w
   names(est) <- object$model$yname
   # compute variance
   design <- object$design
   v <- survey::svyrecvar(resid_winsorized * w / sum_w, design$cluster, 
      design$strata, design$fpc, postStrata = design$postStrata)
   res <- list(characteristic = "mean", estimator = "GREG M-estiamtor", 
      estimate = est, variance = v, design = design, call = call)
   class(res) <- "svystat.rob"
   res
}


#FIXME: total
 
