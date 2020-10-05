svyreg_tukey <- function(formula, design, k, var = NULL, na.rm = FALSE, 
   ...)
{
   dat <- .checkreg(formula, design, var, na.rm)
   res <- robsvyreg(dat$x, dat$y, dat$w, k, 2, 0, NULL, dat$var, ...)
   res$design <- design
   res$call <- match.call()
   res$model$intercept <- dat$intercept
   res$model$yname <- dat$yname 
   class(res) <- "svyreg_rob"
   res
}

svyreg_tukeyGM <- function(formula, design, k, type, xwgt, var = NULL, 
   na.rm = FALSE, ...)
{
   type_int <- switch(toupper(type), "MALLOWS" = 1, "SCHWEPPE" = 2)
   if (is.null(type_int))
      stop("Type '", type,"' is not defined\n")

   dat <- .checkreg(formula, design, var, na.rm)

   if (NCOL(xwgt) > 1) {
      xwgt <- as.numeric(xwgt[, 1])
      warning("Only first column of argument 'xwgt' is used\n")
   }
   if (length(xwgt) != length(dat$y))
      stop("Argument 'xwgt' is not of length n\n")

   res <- robsvyreg(dat$x, dat$y, dat$w, k, 2, type_int, xwgt, dat$var, ...)
   res$design <- design
   res$call <- match.call()
   res$model$intercept <- dat$intercept
   res$model$yname <- dat$yname 
   res$model$xwgt <- xwgt 
   class(res) <- "svyreg_rob"
   res
} 
