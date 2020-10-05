weighted_mean_huber <- function(x, w, k = 1.5, type = "rwm", asym = FALSE, 
   info = FALSE, na.rm = FALSE, ...)
{
   dat <- .check(x, w, na.rm); if (is.null(dat)) return(NA)
   psi <- ifelse(asym, 1, 0)
   if (type == "rwm"){
      res <- robsvyreg(rep(1, dat$n), dat$x, dat$w, k, psi, 0, NULL, NULL, ...)
   }else if (type == "rht"){
      res <- robsvyreg(mean(dat$w) / dat$w, dat$x, dat$w, k, psi, 0, NULL, x, 
	 ...)
   }else{
      stop(paste0("Method '", type, "' does not exist\n"), call. = FALSE)
   }
   if (length(res) == 1) return(NA)
   if (info){
      res$model[c("n", "p")] <- NULL 
      if (type == "rwm") res$model$x <- NULL
      res$characteristic <- "mean"
      res$estimator = paste0("Huber M-estimator (type = ", type, ifelse(asym, 
	 "; asym. psi", ""), ")")
      res$robust[c("Epsi2", "Epsiprime")] <- NULL
      res$call <- match.call()
      return(res)
   }else{
      return(res$estimate)
   }
}

weighted_total_huber <- function(x, w, k = 1.5, type = "rwm", asym = FALSE, 
   info = FALSE, na.rm = FALSE, ...)
{
   res <- weighted_mean_huber(x, w, k, type, asym, info, na.rm, ...)
   if (length(res) == 1){
      res <- res * sum(w)
   }else{
      res$characteristic <- "total"
      res$estimate <- res$estimate * sum(w)
      res$call <- match.call()
   }
   return(res)
} 
