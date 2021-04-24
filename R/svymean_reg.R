svymean_reg <- function(object, auxiliary, k = NULL, check.names = TRUE)
{
    if (!inherits(object, "svyreg_rob"))
        stop("svymean_reg is not defined for this object\n", call. = FALSE)

    beta <- object$estimate
    p <- object$model$p

    if (is.data.frame(auxiliary))
        auxiliary <- as.matrix(auxiliary)

    if (is.matrix(auxiliary)) {
        if (NCOL(auxiliary) < NROW(auxiliary))
	        auxiliary <- t(auxiliary)

        if (NROW(auxiliary) > 1)
            stop("Dimension of argument 'auxiliary' is wrong\n", call. = FALSE)

        cnames <- colnames(auxiliary)
        auxiliary <- c(auxiliary)
        names(auxiliary) <- cnames
    }

    if (length(auxiliary) != p)
        stop("Dim. of 'auxiliary' is wrong; it should be of dim. ", p, " not ",
            length(auxiliary), "\n", call. = FALSE)

    cnames <- names(auxiliary)
    intercept <- object$model$intercept
    if (check.names && !is.null(cnames)) {
        at <- ifelse(intercept, 2, 1)
        m <- match(cnames[at:p], names(beta)[at:p])

        if (is.na(any(m))) {
	        warning("Names of 'auxiliary' do not match (check.names = TRUE)\n",
	            immediate. = TRUE, call. = FALSE)
        } else {
	        auxiliary <- if (at == 1)
	            auxiliary[m]
	        else
	            auxiliary[c(1, m + 1)]
        }
   }

    # estimate
    est <- sum(auxiliary * beta)
    names(est) <- object$model$yname

   # GREG correction
    w <- object$model$w; sum_w <- sum(w)
    if (intercept == 0) {
        est <- if (is.null(k))
            est +  sum(w * object$residuals) / sum_w
        else
	        est	#FIXME:
    }

    # variance estimate  # NOTE:: check
    design <- object$design
    v <- survey::svyrecvar(object$residuals * w / sum_w , design$cluster,
        design$strata, design$fpc, postStrata = design$postStrata)

    object$characteristic <- "mean"
    object$estimate <- est
    object$variance <- v
    object$call <- match.call()
    class(object) <- "svystat_rob"
    object
}


# auxiliary <- c(100, 200)
#
# auxiliary <- matrix(c(100, 200), nrow = 2)
# rownames(auxiliary) <- c("(Intercept)", "employment")
#
#
#
# auxiliary <- matrix(c(100, 200), ncol = 2)
# colnames(auxiliary) <- c("", "employment")
#
# svymean_reg(object, auxiliary)
#
#setwd("C:/My/code/robsurvey/R")
#source("svymean_reg.R")
