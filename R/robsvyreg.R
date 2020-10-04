#' Internal function for regression GM-estimator 
#' 
#' \strong{Internal} function for robust survey regression GM-estimator; this function is \strong{only} intended for internal use. The function does \strong{not} check or validate the arguments. In particular, missing values in the data may make the function crash. 
#' 
#' Not documented
#' 
#' @param x \code{[array]} design matrix (no \code{NA}'s allowed). 
#' @param y \code{[numeric vector]} dependent variable (no \code{NA}'s allowed). 
#' @param w \code{[numeric vector]} weights (no \code{NA}'s allowed).
#' @param psi \code{[integer]} psi-functions: \code{0}: Huber, \code{1}: asymmetric Huber, \code{2}: Tukey biweight.
#' @param k \code{[double]} robustness tuning constant (\eqn{0 < k \leq \infty}{0 < k <= Inf}).
#' @param type \code{[integer]} type of estimator; \code{0}: M-estimator; \code{1}: Mallows and \code{2}: Schweppe type GM-estimator.
#' @param xwgt \code{[numeric vector]} weights for design space used in GM-estimators (default: \code{NULL}, (no \code{NA}'s allowed).
#' @param var \code{[numeric vector]} heteroscedastic variance (default: \code{NULL}).
#' @param maxit \code{[integer]} maximum number of iterations to use (default: \code{100}). 
#' @param tol \code{[double]} numerical tolerance criterion to stop the iterations (default: \code{1e-05}). 
#' @param k_Inf \code{[integer]} numerical value that represents \code{Inf} (default: \code{1e+05}).
#' @param ... additional arguments passed to the method (see \code{svyreg_control}).  
#' @return \code{[list]} 
#' @useDynLib robsurvey rwlslm 
#' @export 
robsvyreg <- function(x, y, w, k, psi, type, xwgt = NULL, var = NULL,
   ...)
{
   ctrl <- svyreg_control(...)
   if (k <= 0) stop("Argument 'k' must be > 0\n", call. = FALSE)
   if (k == Inf) k <- ctrl$k_Inf 
   n <- length(y); p <- NCOL(x) 

   if (is.null(xwgt)) 
      xwgt <- rep(1, n)

   # account for heteroscedasticity
   if (!is.null(var)){
      x <- x / sqrt(var); y <- y / sqrt(var)
   }

   tmp <- .C("rwlslm", x = as.double(x), y = as.double(y), w = as.double(w), 
      resid = as.double(numeric(n)), robwgt = as.double(numeric(n)), 
      xwgt = as.double(xwgt), n = as.integer(n), p = as.integer(p), 
      k = as.double(k), beta = as.double(numeric(p)), 
      scale = as.double(numeric(1)), maxit = as.integer(ctrl$maxit), 
      tol = as.double(ctrl$tol), psi = as.integer(psi), 
      type = as.integer(type), Epsi2 = as.double(numeric(1)), 
      Epsiprime = as.double(numeric(1)), PACKAGE = "robsurvey")

   psi_fun <- switch(psi + 1, "Huber", "asymHuber", "Tukey")
   names(tmp$beta) <- colnames(x) 
   list(
      characteristic = "regression", 
      estimator = paste0("Survey regression ", switch(type + 1, "", 
	 "Mallows G", "Schweppe G") ,"M-estimator (", psi_fun," psi, k = ", 
	 k, ")"),
      estimate = tmp$beta, 
      variance = NA,
      robust = list(psifunction = psi_fun, 
	 k = k, robweights = tmp$robwgt, outliers = 1 * (tmp$robwgt < 1), 
	 Epsi2 = tmp$Epsi2, Epsiprime = tmp$Epsiprime, scale = tmp$scale^2), 
      optim = list(converged = (tmp$maxit != 0), niter = ifelse(tmp$maxit == 0, 
	 ctrl$maxit, tmp$maxit), tol = ctrl$tol), 
      residuals = tmp$resid, 
      model = list(x = x, y = y, w = w, var = var, n = n, p = p), 
      design = NA, 
      call = NA)
}

#' @rdname robsvyreg 
#' @export 
svyreg_control <- function(tol = 1e-5, maxit = 100, k_Inf = 1e5, ...)
{
   list(tol = unname(tol), maxit = unname(maxit), k_Inf = unname(k_Inf))
} 
