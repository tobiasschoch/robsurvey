#' Regression predictor of the mean (generalized regression estimator) 
#' 
#' Regression predictor of the mean 
#'
#' Details
#'  
#' @param object fitted regression model (object of class \code{svyreg}). 
#' @param mean_auxiliary \code{[numeric vector]} population means of the auxiliary variables. 
#' @param ... additional arguments passed to the method. 
#' @return object of class \code{svystat.rob} 
#' @export 
svymean_reg <- function(object, ...)
{
   UseMethod("svymean_reg", object)
}

#' @rdname svymean_reg 
#' @export 
svymean_reg.svyreg_rob <- function(object, mean_auxiliary, ...)
{
   call <- match.call()
   call[[1]] <- substitute(svymean.reg)
   # check dimensions
   if (length(mean_auxiliary) != object$model$p){
      stop("Dimension of argument 'mean_auxiliary' is not correct\n")
   }
   # greg estimate
   w <- object$model$w; sum_w <- sum(w)
   est <- sum(mean_auxiliary * object$estimate) 
   if (object$model$intercept == 0){
      est <- est + sum(w * object$residuals) / sum_w
   }
   names(est) <- object$model$yname

   # variance estimate
   design <- object$design 
   v <- survey::svyrecvar(object$residuals * w / sum_w , design$cluster, 
      design$strata, design$fpc, postStrata = design$postStrata)
   res <- list(characteristic = "mean", estimator = "GREG", estimate = est, 
      variance = v, design = design, call = call)
   class(res) <- "svystat.rob"
   res
}

# FIXME: total 
