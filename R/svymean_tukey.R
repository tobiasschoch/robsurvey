svymean_tukey <- function(x, design, k = 1.5, type = "rwm", na.rm = FALSE, ...)
{
   dat <- .checkformula(x, design)
   res <- weighted_mean_tukey(dat$y, dat$w, k, type, info = TRUE, na.rm, ...)
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

svytotal_tukey <- function(x, design, k = 1.5, type = "rwm", na.rm = FALSE, ...)
{
   res <- svymean_tukey(x, design, k, type, na.rm, ...)  
   sum_w <- sum(res$model$w) 
   res$estimate <- res$estimate * sum_w 
   res$variance <- res$variance * sum_w^2 
   res$characteristic <- "total"
   res$call <- match.call()
   res
} 
