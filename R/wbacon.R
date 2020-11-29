wBACON <- function(x, w = NULL, alpha = 0.95, intercept = FALSE, na.rm = FALSE, 
   maxiter = 50, verbose = FALSE) 
{ 
   if (!is.matrix(x))
      x <- as.matrix(x)

   if (intercept) {
      x <- x[, - which(apply(x, 2, var) == 0)]
   }

   n <- nrow(x); p <- ncol(x)

   if (is.null(w)) w <- rep(1, n)

   stopifnot(n > p, p > 1,  0 < alpha, alpha < 1, n == length(w))

   # NA treatment
   cc <- stats::complete.cases(x, w) 
   if (sum(cc) != n) {
      if (na.rm) { 
	 x <- x[cc, ]
	 w <- w[cc] 
      } else 
	 stop("Data must not contain missing values; see argument 'na.rm'\n", 
	    call. = FALSE)
   } 
   n <- nrow(x)

   # check if any element is not finite 
   chk <- sum(is.finite(c(x, w))) != (1 + p) * n 
   if (chk) 
      stop("Some observations are not finite\n", call. = FALSE)
   
   # compute weighted BACON algorithm
   tmp <- .C("wbacon", x = as.double(x), w = as.double(w), 
      center = as.double(numeric(p)), scatter = as.double(numeric(p * p)),
      dist = as.double(numeric(n)), n = as.integer(n), p = as.integer(p),
      alpha = as.double(alpha), subset = as.integer(rep(0, n)), 
      cutoff = as.double(numeric(1)), maxiter = as.integer(abs(maxiter)),
      verbose = as.integer(verbose), PACKAGE = "robsurvey")

   tmp$cov <- matrix(tmp$scatter, ncol = p)
   tmp$verbose <- NULL
   tmp$converged <- ifelse(tmp$maxiter < maxiter, TRUE, FALSE)
   tmp$call <- match.call()
   names(tmp$center) <- colnames(x)
   colnames(tmp$cov) <- colnames(x)
   rownames(tmp$cov) <- colnames(x)
   class(tmp) <- "robmv"
   tmp
}

print.robmv <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
   cat("\nWeighted BACON: Robust location, covariance, and distances\n")
   if (x$converged) {
      cat(paste0("Converged in ", x$maxiter, " iterations (alpha = ", x$alpha, 
	 ")\n\n"))
      cat("Robust estimate of location:\n")
      print(x$center, digits)
      cat("\nRobust estimate of covariance:\n")
      print(x$cov, digits = digits)
      cat(paste0("\nDistances (cutoff: ", format(x$cutoff, digits = digits), 
	 "):\n"))
      print(summary(x$dist), digits = digits) 
   } else
      cat(paste0("Algorithm did not converge in ", x$maxiter," iterations!\n"))

   cat("\n")
}

distance <- function(x){
   if (!inherits(x, "robmv"))
      cat("not defined for this type of argument\n")
   else
      x$dist      
} 
