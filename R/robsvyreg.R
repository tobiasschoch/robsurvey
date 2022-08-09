# workhorse function for robust regression
robsvyreg <- function(x, y, w, k, psi, type, xwgt, var = NULL,
    verbose = TRUE, ...)
{
    stopifnot(is.numeric(k))
    n <- length(y); p <- NCOL(x)
    model <- list(x = x, y = y, w = w, var = var, xwgt = xwgt, n = n, p = p)

    ctrl <- svyreg_control(...)
    if (k <= 0)
        stop("Argument 'k' must be > 0\n", call. = FALSE)
    if (k == Inf)
        k <- ctrl$k_Inf

    # account for heteroscedasticity
    if (!is.null(var)) {
        x <- x / sqrt(var)
        y <- y / sqrt(var)
    }

    # initialization
    if (is.null(ctrl$init)) {
        init_flag <- 1 # initialization by weighted least squares
        beta <- numeric(p)
    } else {
        init_flag <- 0 # regression is initialized by ctrl$init
        beta <- ctrl$init
        if (length(as.vector(beta)) != p || !is.numeric(beta))
            stop("Argument 'init' must be a numerical p-vector\n",
                call. = FALSE)
        if (sqrt(sum(beta^2)) < sqrt(.Machine$double.eps))
            stop("Euclidean norm of 'init' is zero (or nearly so)\n",
                call. = FALSE)
    }

    tmp <- .C("rwlslm", x = as.double(x), y = as.double(y), w = as.double(w),
        resid = as.double(numeric(n)), robwgt = as.double(numeric(n)),
        xwgt = as.double(xwgt), n = as.integer(n), p = as.integer(p),
        k = as.double(k), beta = as.double(beta),
        scale = as.double(numeric(1)), tol = as.double(ctrl$tol),
        maxit = as.integer(ctrl$maxit), psi = as.integer(psi),
        type = as.integer(type), init = as.integer(init_flag),
        mad_center = as.integer(ctrl$mad_center), verbose = as.integer(verbose),
        used_iqr = as.integer(0), PACKAGE = "robsurvey")
    # Note: The residuals are (y - x*beta) / sqrt(var)

    converged <- (tmp$maxit != 0)

    psi_fun <- switch(psi + 1, "Huber", "asymHuber", "Tukey")
    names(tmp$beta) <- colnames(x)
    list(characteristic = "regression",
        estimator = list(string = paste0("Survey regression ",
            switch(type + 1, "", "Mallows G", "Schweppe G"),
            "M-estimator (", psi_fun," psi, k = ",
	        k, ")"), type = type, psi = psi, psi_fun = psi_fun, k = k),
        estimate = tmp$beta, scale = tmp$scale,
        robust = list(robweights = tmp$robwgt, outliers = 1 * (tmp$robwgt < 1)),
        optim = list(converged = converged,
            niter = ifelse(tmp$maxit == 0, ctrl$maxit, tmp$maxit),
            tol = ctrl$tol, used_iqr = tmp$used_iqr),
        residuals = tmp$resid, model = model, design = NA, call = match.call())
}
# control function for robust regression
svyreg_control <- function(tol = 1e-5, maxit = 100, k_Inf = 1e6, init = NULL,
    mad_center = TRUE, ...)
{
    if (!is.numeric(tol))
        stop("Argument 'tol' must be of type 'numeric'\n", call. = FALSE)
    if (maxit <= 0 || !is.numeric(maxit))
        stop("Argument 'maxit' must be a positive integer\n", call. = FALSE)
    if (!is.logical(mad_center))
        stop("Argument 'mad_center' must be of type 'logical'\n",
            call. = FALSE)
    if (!is.numeric(k_Inf))
        stop("Argument 'k_Inf' must be of type 'numeric'\n",
            call. = FALSE)
    list(tol = unname(tol), maxit = unname(as.integer(maxit)),
        k_Inf = unname(k_Inf), init = unname(init),
        mad_center = unname(mad_center))
}
