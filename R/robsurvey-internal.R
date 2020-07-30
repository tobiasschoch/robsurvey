# some sanity checks
.check <- function(x, w, na.rm)
{
   if (is.factor(x) || is.factor(w) || is.data.frame(x)){
      stop("Arguments data and weights must be numeric vectors\n")
   }
   n <- length(x); nw <- length(w)
   if (nw != n) stop("Data vector and weights are not of the same dimension\n", 
      call. = FALSE)
   if (n == 0){
      return(NA)
   }
   # check for missing values
   cc <- stats::complete.cases(x, w)
   if (sum(cc) != n) {
      if (na.rm) {
	 x <- x[cc]; w <- w[cc]
      } else{
	 return(NULL)
      }
   }
   n <- length(x) 
   # check if data vector and weights are finite
   if (sum(is.finite(c(x, w))) != 2 * n) {
      warning("Some observations are not finite\n", call. = FALSE, 
	 immediate. = TRUE)
      return(NULL)      
   }
   list(x = x, w = w, n = n)
}

.checkformula <- function(x, design)
{
   if (class(x) == "formula"){
      mf <- stats::model.frame(x, design$variables, na.action = stats::na.pass)
      n <- nrow(mf)
      if (ncol(mf) > 1){
	 stop("Argument 'y' must be a formula of one single variable", 
	    call. = FALSE)
      }
      xname <- names(mf)
      xdat <- mf[[1]]
   }else{
      if (is.character(x)){
	 xname <- x
	 xdat <- design$variables[, x]
	 n <- length(xdat)
	 if (any(is.na(xdat))) stop(paste0("Variable '", xname, 
	    "' must not contain NA's\n"), call. = FALSE)
      }else{
	 stop("svymean_winsorized is not defined for object of class: ", 
	    class(x), "\n", call. = FALSE)
      }
   }
   list(xname = xname, x = xdat, w = as.numeric(1 / design$prob))   
}

# influence function winsorized mean
.infl_winsorized <- function(x, w, LB, UB, wm, ngrid = 401)
{
   qs <- weighted_quantile(x, w, probs = c(LB, UB))
   # density estimates at 'LB' and '1 - UB'
   bwd <- KernSmooth::dpik(x, scalest = "minim", level = 2, kernel = "normal", 
      canonical = FALSE, gridsize = ngrid, range.x = range(x), truncate = TRUE)
   at <- seq(min(x), max(x), length = ngrid)
   nx <- rowsum(c(rep(0, ngrid), w), c(1:ngrid, findInterval(x, at)))
   dens <- KernSmooth::locpoly(rep(1, ngrid), nx * ngrid / (diff(range(x)) * 
      sum(w)), binned = TRUE, bandwidth = bwd, range.x = range(x))
   f_LB <- dens$y[floor(ngrid * LB) + 1]
   f_UB <- dens$y[floor(ngrid * UB)]
   # influence function
   infl <- pmin.int(qs[2] + (1 - UB) / f_UB, pmax.int(qs[1] - LB / f_LB, x))
   w <- (UB - LB) * wm + LB * qs[1] + (1 - UB) * qs[2] 
   # return
   infl - w + LB^2 / f_LB + (1 - UB)^2 / f_UB
}

# influence function trimmed mean
.infl_trimmed <- function(x, w, LB, UB, tm)
{
   qs <- weighted_quantile(x, w, probs = c(LB, UB))
   # influence function
   infl <- pmin.int(qs[2], pmax.int(qs[1], x))
   w <- (UB - LB) * tm + LB * qs[1] + (1 - UB) * qs[2] 
   (infl - w) / (UB - LB)
}

# onAttach function
.onAttach <- function(libname, pkgname)
{
   quietly <- TRUE

   # this code is taken from Henrik Bengtsson
   # https://github.com/HenrikBengtsson/R.methodsS3/ => pkgStartupMessage.R 
   tryCatch({
      # The default, if not found
      quietly <- formals(base::library)$quietly

      # Identify the environment/frame of interest by making sure
      # it at least contains all the arguments of source().
      argsToFind <- names(formals(base::library))

      # Scan the call frames/environments backwards...
      srcfileList <- list()
      for (ff in sys.nframe():0) {
	 env <- sys.frame(ff)

	 # Does the environment look like a library() environment?
	 exist <- sapply(argsToFind, FUN = exists, envir = env, 
	    inherits = FALSE)
	 if (!all(exist)) {
	 # Nope, then skip to the next one
	 next
	 }

	 # Was argument 'quietly' specified?
	 missing <- eval(expression(missing(quietly)), envir = env)
	 if (!missing) {
	    quietly <- get("quietly", envir = env, inherits = FALSE)
	    break
	 }

	 # ...otherwise keep searching due to nested library() calls.
      } 
   }, error = function() {})

   if (!quietly){
      packageStartupMessage("\n                  88
		  88
		  88
     e8d88 .d8b.  8888b.   ___ _   _ _ ____   _____ _   _
     8P'  d8' '8b 88 '8b  / __| | | | '__\\ \\ / / _ \\ | | |
     88   Y8. .8P 88  dP  \\__ \\ |_| | |   \\ V /  __/ |_| |
     88    'Y8P'  88e8P'  |___/\\__,_|_|    \\_/ \\___|\\__, |
						     __/ |
					version 0.2 |___/\n")
   }
} 
