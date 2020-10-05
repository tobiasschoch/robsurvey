svymean_huber <- function(x, design, k = 1.5, type = "rwm", asym = FALSE, 
   na.rm = FALSE, ...)
{
   dat <- .checkformula(x, design)
   res <- weighted_mean_huber(dat$y, dat$w, k, type, asym, info = TRUE, na.rm, 
      ...)
   # modify residuals for type 'rht' (only for variance estimation)
   if (type == "rht"){
      r <- sqrt(res$model$var) * res$model$y - res$estimate 
   } else {
      r <- res$residuals
   }
   # compute variance

# FIXME: take y and from res
   infl <- res$robust$robweights * r * dat$w / sum(dat$w) 
   res$variance <- survey::svyrecvar(infl, design$cluster, design$strata, 
        design$fpc, postStrata = design$postStrata)
   names(res$estimate) <- dat$yname
   res$call <- match.call()
   res$design <- design
   class(res) <- "svystat_rob"
   res
}

svytotal_huber <- function(x, design, k = 1.5, type = "rwm", asym = FALSE, 
   na.rm = FALSE, ...)
{
   res <- svymean_huber(x, design, k, type, asym, na.rm, ...)  
   sum_w <- sum(res$model$w) 
   res$estimate <- res$estimate * sum_w 
   res$variance <- res$variance * sum_w^2 
   res$characteristic <- "total"
   res$call <- match.call()
   res
}

# fixed_downweighting <- function(object, at = 0.95){
#    if (class(object) != "svystat_rob") stop("method not supported for class: ", 
#       class(object), "\n") 
#    if (object$call[[1]] == "svymean_huber" && 
#       object$call[[1]] == "svytotal_huber"){
#       stop("method: '", object$call[[1]], "' is not supported\n") 
#    }
#    stopifnot(at > 0.5, at < 1)
#    ctrl <- svyreg_control()
#
#    cl <- object$call 
#    foo <- function(k, at){
#       cl$k <- k
#       tmp <- eval(cl)
#       mean(tmp$robust$robweights) - at
#    }
#    res <- uniroot(foo, interval = c(0.0001, ctrl$k_Inf), at = at)
#    res$root
# }
# 
