robsvyreg <- function(x, y, w, k, psi, type, xwgt = NULL, var = NULL, ...)
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
      scale = as.double(numeric(1)), tol = as.double(ctrl$tol), 
      maxit = as.integer(ctrl$maxit), psi = as.integer(psi), 
      type = as.integer(type), PACKAGE = "robsurvey")

   psi_fun <- switch(psi + 1, "Huber", "asymHuber", "Tukey")
   names(tmp$beta) <- colnames(x) 
   list(
      characteristic = "regression", 
      estimator = list(string = paste0("Survey regression ", switch(type + 1, 
	 "", "Mallows G", "Schweppe G") ,"M-estimator (", psi_fun," psi, k = ", 
	 k, ")"), type = type, psi = psi, psi_fun = psi_fun, k = k),
      estimate = tmp$beta, 
      scale = tmp$scale, 
      robust = list(robweights = tmp$robwgt, outliers = 1 * (tmp$robwgt < 1)), 
      optim = list(converged = (tmp$maxit != 0), niter = ifelse(tmp$maxit == 0, 
	 ctrl$maxit, tmp$maxit), tol = ctrl$tol), 
      residuals = tmp$resid, 
      model = list(x = x, y = y, w = w, var = var, xwgt = xwgt, n = n, p = p), 
      design = NA, 
      call = NA)
}

svyreg_control <- function(tol = 1e-5, maxit = 100, k_Inf = 1e5, ...)
{
   list(tol = unname(tol), maxit = unname(maxit), k_Inf = unname(k_Inf))
} 
