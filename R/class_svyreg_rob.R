# Slots of objects of class 'svyreg_rob'
#  + characteristic: "regression"
#  + estimator [list]
#     + string: [char]
#     + type: [int]
#     + psi: [int]
#     + psi_fun: [char]
#     + k: [numeric]
#  + estimate: [numeric] vector of estimated regression coefficients
#  + scale: [numeric] scale estimate
#  + robust [list]
#     + robweights: [numeric] robustness weights
#     + outliers: [numeric] indicator variable
#  + optim [list]
#     + converged: [logical]
#     + niter: [int] number of IRWLS iterations
#     + tol: [numeric] numerical tolerance criterion (IRLWS)
#     + used_iqr: [int] 1 = scale estimated by IQR not MAD
#  + residuals: [numeric]
#  + model [list]
#     + x: [matrix] design matrix
#     + y: [numeric] response variable
#     + w: [numeric] sampling weights
#     + var: [numeric] heteroscedasticity variances
#     + xwgt: [numeric] weights in the model's design space (GM-estimator)
#     + n [int] number of observations
#     + p [int] number of independent variables
#  + design: [survey.design object without 'variables']
#  + call: [call object]
#  + terms: [terms object]

# print method for robust regression object
print.svyreg_rob <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
    if (any(is.na(x$estimate))) {
        cat(paste0("\n", x$estimator$string, "\n"))
        cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
            "\n", sep = "")
        cat("\nCoefficients:\n")
        print.default(format(x$estimate, digits = digits), print.gap = 2L,
                      quote = FALSE)
        cat("\nScale estimate:", NA, "\n")
        return(invisible(x))
    }
    converged <- x$optim$converged
    if (is.null(converged) || converged) {
        cat(paste0("\n", x$estimator$string, "\n"))
        cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
            "\n", sep = "")
        if (is.finite(x$estimator$k))
            cat("\nIRWLS converged in", x$optim$niter, "iterations\n\n")
        else
            cat("\n")

        if (length(x$estimate)) {
            cat("Coefficients:\n")
            print.default(format(x$estimate, digits = digits), print.gap = 2L,
                          quote = FALSE)
        }
        cat("\nScale estimate:", format(signif(x$scale, digits)),
            if (!is.null(x$optim$used_iqr))
                if (x$optim$used_iqr) "(weighted IQR)" else "(weighted MAD)"
            else
                "", "\n")
    } else {
        cat("\nIRWLS not converged in", x$optim$niter, "iterations (with tol =",
            x$optim$tol, ")\n")
    }
    # return
    invisible(x)
}
# summary method for robust regression object
summary.svyreg_rob <- function(object, mode = c("design", "model", "compound"),
                               digits = max(3L, getOption("digits") - 3L), ...)
{
    if (any(is.na(object$estimate))) {
        cat("\nCall:\n")
        print(object$call)
        cat("\nResiduals:\n")
        print(NA)
        cat("\nCoefficients:\n")
        NAs <- object$estimate
        printCoefmat(cbind(Estimate = NAs, 'Std. Error' = NAs, 't value' = NAs,
                           'Pr(>|t|)' = NAs), digits = digits)
        cat("\nResidual standard error:", NA, "on", NA, "degrees of freedom\n")
        return(invisible(object))
    }
    converged <- object$optim$converged
    if (is.null(converged) || converged) {
        n <- object$model$n; p <- object$model$p
        r <- object$residuals; N <- sum(object$model$w)

        cat("\nCall:\n")
        print(object$call)
        cat("\nResiduals:\n")
        print(summary(r))
        cat("\nCoefficients:\n")

        # covariance matrix
        tmp <- switch(match.arg(mode),
            "model" = .cov_reg_model(object),
            "design" = .cov_reg_design(object),
            "compound" = .cov_reg_compound(object))
        if (tmp$ok == 0)
            stop("\nCovariance estimation failed\n", call. = FALSE)

        res <- list(stddev = tmp$scale, covmat = tmp$cov, n = n, p = p, N = N)
        ns <- colnames(object$model$x)
        dimnames(res$covmat) <- list(ns, ns)

        # print out
        est <- object$estimate
        se <- sqrt(diag(res$covmat))
        tval <- est / se
        printCoefmat(cbind(Estimate = est, 'Std. Error' = se, 't value' = tval,
            'Pr(>|t|)' = 2 * pt(abs(tval), N - p, lower.tail = FALSE)),
            digits = digits)

        cat("\nResidual standard error:", format(signif(res$stddev, digits)),
            "on", N - p, "degrees of freedom\n")

        if (is.finite(object$estimator$k)) {
            cat("\nRobustness weights:\n")
            print(summary(object$robust$robweights))
        }
        invisible(res)
    } else {
        warning(" covariance is not available because",
            "\n regression algorithm did not converge\n", call. = FALSE)
    }
}
# extract variance from robust regression object
vcov.svyreg_rob <- function(object, mode = c("design", "model", "compound"),
                            ...)
{
    converged <- object$optim$converged
    if (is.null(converged) || converged) {
        mat <- switch(match.arg(mode),
            "model" = .cov_reg_model(object)$cov,
            "design" = .cov_reg_design(object)$cov,
            "compound" = .cov_reg_compound(object)$cov)
    } else {
        warning(" covariance is not available because",
                "\n regression algorithm did not converge\n", call. = FALSE)
        mat <- matrix(NA, object$model$p, object$model$p)
    }
    ns <- colnames(object$model$x)
    if (length(ns) == 1)
        ns <- names(object$estimate)
    dimnames(mat) <- list(ns, ns)
    mat
}

# extract standard error from robust regression object
SE.svyreg_rob <- function(object, mode = c("design", "model", "compound"),
                          ...)
{
    sqrt(diag(vcov.svyreg_rob(object, mode, ...)))
}

# model-based covariance matrix of M- and GM-regression estimators
.cov_reg_model <- function(object)
{
    # robustness tuning constant
    k <- object$estimator$k
    if (is.infinite(k))
        k <- svyreg_control()$k

    # account for heteroscedasticity
    x <- object$model$x
    if (!is.null(object$model$var))
        x <- x / sqrt(object$model$var)

    # robustness weights
    ui <- robweights(object)

    tmp <- .C(C_cov_reg_model, resid = as.double(object$residuals),
        x = as.double(object$model$x), xwgt = as.double(object$model$xwgt),
        robwgt = as.double(ui), w = as.double(object$model$w), k = as.double(k),
        scale = as.double(object$scale), scale = as.double(numeric(1)),
        n = as.integer(object$model$n), p = as.integer(object$model$p),
        psi = as.integer(object$estimator$psi),
        type = as.integer(object$estimator$type), ok = as.integer(0))
    p <- object$model$p
    if (tmp$ok == 0) {
        warning("Covariance estimation failed", call. = FALSE)
        cov_mat <- matrix(NA, p, p)
    } else {
        cov_mat <- matrix(tmp$x[1:(p * p)], ncol = p)
    }
    # return
    list(ok = tmp$ok, cov = cov_mat, scale = tmp$scale)
}

# design-based covariance matrix of M- and GM-regression estimators
.cov_reg_design <- function(object)
{
    # robustness tuning constant
    k <- object$estimator$k
    if (is.infinite(k))
        k <- svyreg_control()$k

    x <- object$model$x
    n <- NROW(x); p <- NCOL(x)
    # account for heteroscedasticity
    if (!is.null(object$model$var))
        x <- x / sqrt(object$model$var)

    # weights in the model's design space
    xwgt <- object$model$xwgt
    if (is.null(xwgt))
        xwgt <- rep(1, n)

    # scale the weights (prevent overflow)
    w <- object$model$w / sum(object$model$w)

    # robustness weights
    ui <- robweights(object)

    # Q matrix
    r <- object$residuals
    Q <- svyrecvar(w * xwgt * x * ui * r, object$design$cluster,
                   object$design$strata, object$design$fpc)

    # covariance matrix (takes Q matrix as one of its argument)
    tmp <- .C(C_cov_reg_design, x = as.double(x), w = as.double(w),
        xwgt = as.double(xwgt), resid = as.double(r),
        scale = as.double(object$scale), k = as.double(k),
        psi = as.integer(object$estimator$psi),
        type = as.integer(object$estimator$type), n = as.integer(n),
        p = as.integer(p), ok = as.integer(0), mat = as.double(Q))

    if (tmp$ok == 0) {
        warning("Covariance estimation failed", call. = FALSE)
        cov_mat <- matrix(NA, p, p)
    } else {
        cov_mat <- matrix(tmp$mat, ncol = p, nrow = p)
    }
    # return
    list(ok = 1, cov = cov_mat, scale = object$scale)
}

# compound design-model-based covariance matrix of M- and GM-regression est.
.cov_reg_compound <- function(object)
{
    gamma_m <- .cov_reg_model(object)
    gamma_d <- .cov_reg_design(object)
    # sampling fraction
    f <- object$model$n / sum(object$model$w)
    if (gamma_d$ok * gamma_m$ok == 0)
        cov_mat <- matrix(NA, object$model$p, object$model$p)
    else
        cov_mat <- gamma_d$cov + f * gamma_m$cov
    # return
    list(ok = 1, cov = cov_mat, scale = object$scale)
}

# extract coefficients from robust regression object
coef.svyreg_rob <- function(object, ...)
{
    object$estimate
}

# extract residuals from robust regression object
residuals.svyreg_rob <- function(object, ...)
{
    object$residuals
}
# extract fitted values from robust regression object
fitted.svyreg_rob <- function(object, ...)
{
    object$model$y - object$residuals
}
# extract robustness weights from robust regression object
robweights.svyreg_rob <- function(object)
{
    tmp <- object$robust$robweights
    if (is.null(tmp))
        rep(1, object$model$n)
    else
        tmp
}
# plot method for robust regression object
plot.svyreg_rob <- function(x, which = 1L:4L, hex = FALSE,
    caption = c("Standardized residuals vs. Fitted Values",
        "Normal Q-Q", "Response vs. Fitted values",
        "Sqrt of abs(Residuals) vs. Fitted Values"),
	panel = if (add.smooth) function(x, y, ...) panel.smooth(x, y,
        iter = iter.smooth, ...) else points,
    sub.caption = NULL, main = "", ask = prod(par("mfcol")) < length(which) &&
        dev.interactive(),
    ..., id.n = 3, labels.id = names(residuals(x)), cex.id = 0.75,
    qqline = TRUE, add.smooth = getOption("add.smooth"), iter.smooth = 3,
	label.pos = c(4, 2), cex.caption = 1, cex.oma.main = 1.25)
{
    if (!is.numeric(which) || any(which < 1) || any(which > 4))
        stop("'which' must be in 1:4")

    show <- rep(FALSE, 4)
    show[which] <- TRUE

    r <- residuals(x)
    n <- length(r)
    yh <- fitted(x)
    w <- x$model$w
    y <- x$model$y

    # Standardized residuals
    rs <- r / x$scale

    if (is.null(id.n)) {
        id.n <- 0
    } else {
        id.n <- as.integer(id.n)
        if (id.n < 0L || id.n > n)
            stop(gettextf("'id.n' must be in {1,..,%d}", n), domain = NA)
    }
    if (id.n > 0L) {
        if (is.null(labels.id))
            labels.id <- paste(1L:n)
        iid <- 1L:id.n
        show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
        if (any(show[2L:3L]))
            show.rs <- sort.list(abs(rs), decreasing = TRUE)[iid]
        text.id <- function(x, y, ind, adj.x = TRUE)
		{
            labpos <- if (adj.x)
                label.pos[1 + as.numeric(x > mean(range(x)))]
            else
				3

            text(x, y, labels.id[ind], cex = cex.id, xpd = TRUE,
                pos = labpos, offset = 0.25)
        }
    }
    getCaption <- function(k)
	{
		if (length(caption) < k)
        	NA_character_
	    else
			as.graphicsAnnot(caption[[k]])
	}
    if (is.null(sub.caption)) {
        cal <- x$call
        if (!is.na(m.f <- match("formula", names(cal)))) {
            cal <- cal[c(1, m.f)]
            names(cal)[2L] <- ""
        }
        cc <- deparse(cal, 80)
        nc <- nchar(cc[1L], "c")
        abbr <- length(cc) > 1 || nc > 75
        sub.caption <- if (abbr)
            paste(substr(cc[1L], 1L, min(75L, nc)), "...")
        else
			cc[1L]
    }
    one.fig <- prod(par("mfcol")) == 1
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
	# 1 Standardized residuals vs. Fitted values
    if (show[1L]) {
        ylim <- range(rs, na.rm = TRUE)
        if (id.n > 0)
            ylim <- extendrange(r = ylim, f = 0.08)
        dev.hold()
        if (hex) {
            requireNamespace("hexbin")
            hb <- hexbin::hexbin(yh, rs, ybnds = ylim)
            hvp <- hexbin::plot(hb, xlab = "Fitted values",
                                ylab = "Standardized residuals", main = main)
            hexbin::hexVP.abline(hvp$plot, h = 0, lty = 3, col = "gray")
            hexbin::hexVP.loess(hb, hvp = hvp$plot, span = 2 / 3, ...)
        } else {
            plot(yh, rs, xlab = "Fitted values",
                 ylab = "Standardized residuals", main = main, ylim = ylim,
                 type = "n", ...)
            panel(yh, rs, ...)
            if (one.fig)
                title(sub = sub.caption, ...)
            mtext(getCaption(1), 3, 0.25, cex = cex.caption)
            if (id.n > 0) {
                y.id <- rs[show.r]
                y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ") / 3
                text.id(yh[show.r], y.id, show.r)
            }
            abline(h = 0, lty = 3, col = "gray")
        }
        dev.flush()
    }
	# 2 Normal Q-Q
    if (show[2L]) {
        ylim <- range(rs, na.rm = TRUE)
        ylim[2L] <- ylim[2L] + diff(ylim) * 0.075
        dev.hold()
        qq <- qqnorm(rs, main = main, ylab = "Standardized residuals",
                     ylim = ylim, ...)
        if (qqline)
            qqline(rs, lty = 3, col = "gray50")
        if (one.fig)
            title(sub = sub.caption, ...)
        mtext(getCaption(2), 3, 0.25, cex = cex.caption)
        if (id.n > 0)
            text.id(qq$x[show.rs], qq$y[show.rs], show.rs)
        dev.flush()
    }
    # 3 Response vs. Fitted values
    if (show[3L]) {
        ylim <- range(y, na.rm = TRUE)
        if (id.n > 0)
            ylim <- extendrange(r = ylim, f = 0.08)
        dev.hold()
        if (hex) {
            requireNamespace("hexbin")
            hb <- hexbin::hexbin(yh, y, ybnds = ylim)
            hvp <- hexbin::plot(hb, xlab = "Fitted values", ylab = "Response",
                                main = main)
            hexbin::hexVP.abline(hvp$plot, h = 0, lty = 3, col = "gray")
            hexbin::hexVP.loess(hb, hvp = hvp$plot, span = 2 / 3, ...)
        } else {
            plot(yh, y, xlab = "Fitted values", ylab = "Response",
                 main = main, ylim = ylim, type = "n", ...)
            panel(yh, y, ...)
            if (one.fig)
                title(sub = sub.caption, ...)
            mtext(getCaption(3), 3, 0.25, cex = cex.caption)
            if (id.n > 0) {
                y.id <- y[show.r]
                y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ") / 3
                text.id(yh[show.r], y.id, show.r)
            }
            abline(h = 0, lty = 3, col = "gray")
            abline(0, 1, lty = 2, col = "gray")
        }

        dev.flush()
    }
    # 4 Sqrt of abs(Residuals) vs. Fitted values
    if (show[4L]) {
        sqrtabsr <- sqrt(abs(r))
        ylim <- c(0, max(sqrtabsr, na.rm = TRUE))
        yhn0 <- yh
        dev.hold()
        if (hex) {
            requireNamespace("hexbin")
            hb <- hexbin::hexbin(yhn0, sqrtabsr, ybnds = ylim)
            hvp <- plot(hb, xlab = "Fitted values",
                        ylab = "Sqrt of abs(Residuals)", main = main)
            hexbin::hexVP.loess(hb, hvp = hvp$plot, span = 2 / 3, ...)
        } else {
            plot(yhn0, sqrtabsr, xlab = "Fitted values",
                 ylab = "Sqrt of abs(Residuals)", main = main, ylim = ylim,
                 type = "n", ...)
            panel(yhn0, sqrtabsr, ...)
            if (one.fig)
                title(sub = sub.caption, ...)
            mtext(getCaption(4), 3, 0.25, cex = cex.caption)
            if (id.n > 0)
                text.id(yhn0[show.r], sqrtabsr[show.r], show.r)
        }
        dev.flush()
    }
    if (!one.fig && par("oma")[3L] >= 1)
        mtext(sub.caption, outer = TRUE, cex = cex.oma.main)
    invisible()
}
