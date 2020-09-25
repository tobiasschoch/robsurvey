#' Internal function for regression M-estimator 
#' 
#' \strong{Internal} function for Huber type robust survey regression M-estimator; this function is \strong{not} intended for external use  
#' 
#' Not documented
#' 
#' @param x \code{[array]} design matrix. 
#' @param y \code{[numeric vector]} dependent variable. 
#' @param w \code{[numeric vector]} weights.
#' @param psi \code{[integer]} psi-functions: \code{0}: Huber (default), \code{1}: asymmetric Huber, \code{2}: Tukey biweight.
#' @param k \code{[double]} robustness tuning constant (\eqn{0 < k \leq \infty}{0 < k <= Inf}).
#' @param type \code{[integer]} type of estimator; \code{0}: M-estimator; \code{1}: Mallows and \code{2}: Schweppe type GM-estimator
#' @param intercept \code{[integer]} \code{1} if has intercept, otherwise \code{0}.  
#' @param xwgt \code{[numeric vector]} weights for design space used in GM-estimators (default: \code{NULL}).
#' @param var \code{[numeric vector]} heteroscedastic variance (default: \code{NULL}).
#' @param na.rm \code{[logical]} indicating whether \code{NA} values should be removed before the computation proceeds (default: \code{FALSE}). 
#' @param maxit \code{[integer]} maximum number of iterations to use (default: \code{100}). 
#' @param tol \code{[double]} numerical tolerance criterion to stop the iterations (default: \code{1e-05}). 
#' @param psi \code{[character]} \code{"Huber"} or \code{"asymHuber"} (default: \code{psi = "Huber"}). 
#' @param k_Inf \code{[integer]} numerical value that represents \code{Inf} (default: \code{1e+05}).
#' @param ... additional arguments passed to the method (see \code{svyreg_control}).  
#' @return \code{[list]} 
#' @useDynLib robsurvey rwlslm 
#' @export 
robsvyreg <- function(x, y, w, k, psi, type, intercept, xwgt = NULL, var = NULL,
   na.rm, ...)
{
   if (!is.matrix(x)) x <- as.matrix(x)

#FIXME: NA's

   ctrl <- svyreg_control(...)
   if (k <= 0) stop("Argument 'k' must be > 0\n", call. = FALSE)
   if (k == Inf) k <- ctrl$k_Inf 

   if (is.null(xwgt)) { 
      if(type %in% 1:2) 
	 stop("Argument 'xwgt' is missing\n", call. = FALSE)
      else
	 xwgt <- rep(1, length(y))
   }

   n <- length(y); p <- ifelse(is.null(ncol(x)), 1, ncol(x))
   # account for heteroscedasticity
   if (!is.null(var)){
      x <- x / sqrt(var); y <- y / sqrt(var)
   }
   #
   tmp <- .C("rwlslm", x = as.double(x), y = as.double(y), w = as.double(w), 
      resid = as.double(numeric(n)), robwgt = as.double(numeric(n)), 
      xwgt = as.double(xwgt), n = as.integer(n), p = as.integer(p), 
      k = as.double(k), beta = as.double(numeric(p)), 
      scale = as.double(numeric(1)), maxit = as.integer(ctrl$maxit), 
      tol = as.double(ctrl$tol), psi = as.integer(psi), 
      type = as.integer(type), Epsi2 = as.double(numeric(1)), 
      Epsiprime = as.double(numeric(1)))#FIXME: , PACKAGE = "robsurvey")

   psi_fun <- switch(ctrl$psi + 1, "Huber", "asymHuber", "Tukey")
   names(tmp$beta) <- colnames(x) 

   list(
      characteristic = "regression", 
      #FIXME:
      estimator = paste0("Survey regression M-estimator (", psi_fun,", k = ", 
	 k, ")"), 
      estimate = tmp$beta, 
      variance = NA,
      robust = list(psifunction = psi_fun, 
	 k = k, robweights = tmp$robwgt, outliers = 1 * (tmp$robwgt < 1), 
	 Epsi2 = tmp$Epsi2, Epsiprime = tmp$Epsiprime, scale = tmp$scale^2), 
      optim = list(converged = (tmp$maxit != 0), niter = ifelse(tmp$maxit == 0, 
	 ctrl$maxit, tmp$maxit)), 
      residuals = tmp$resid, 
      model = list(y = y, x = x, w = w, var = var, n = n, p = p, 
	 intercept = intercept, yname = NA), 
      design = NA, 
      call = NA)
}

#' @rdname robsvyreg 
#' @export 
svyreg_control <- function(tol = 1e-5, maxit = 100, psi = "Huber", k_Inf = 1e5,
   ...)
{
   #FIXME:
   if(!(psi %in% c("Huber", "asymHuber"))) stop("Function 'psi' must be 
      either 'Huber' or 'asymHuber'\n")
   psi0 <- switch(psi, "Huber" = 0L, "asymHuber" = 1L)
   list(tol= unname(tol), maxit = unname(maxit), psi = unname(psi0), 
      k_Inf = unname(k_Inf))
} 
