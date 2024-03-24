# Style

Tobias Schoch

## General

In terms of style (`R` and `C`code), we mostly follow the Linux Kernel Coding Style; see https://www.kernel.org/doc/html/v4.10/process/coding-style.html. But tabs are 4 characters in our case, not 8.



## Indentation

Function signature and body

```R
svymean_huber <- function(x, design, k, type = "rwm", asym = FALSE,
                          na.rm = FALSE, verbose = TRUE, ...)
{
    
}
```

General function calls

```R
res <- weighted_mean_huber(dat$y, dat$w, k, type, asym, TRUE, FALSE,
                           verbose, ...)
```

Call of `C` and `Fortran` symbols

```R
tmp <- .C(C_rwlslm, x = as.double(x), y = as.double(y), w = as.double(w),
    resid = as.double(numeric(n)), robwgt = as.double(numeric(n)),
    xwgt = as.double(xwgt), n = as.integer(n), p = as.integer(p),
    k = as.double(k), beta = as.double(beta),
    scale = as.double(numeric(1)), tol = as.double(ctrl$tol),
    maxit = as.integer(ctrl$maxit), psi = as.integer(psi),
    type = as.integer(type), init = as.integer(init_flag),
    mad_center = as.integer(ctrl$mad_center), verbose = as.integer(verbose),
    used_iqr = as.integer(0))
```

List

```R
res <- list(characteristic = "total",
	estimator = list(string = paste0("Dalen ", type,
		" estimator (censored at ", censoring, ")"), censoring = censoring),
    estimate = estimate, variance = NA,
    robust = list(xw = xw),
    residuals = NA,
    model = list(y = dat$x, w = dat$w),
    design = NA, call = match.call())
```

Switch

```R
# simple
tmp <- switch(variance, "base" = 0, "wu" = 0.5, "hajek" = 1)

# more elaborate
tmp <- switch(match.arg(mode),
    "model" = .cov_reg_model(object),
    "design" = .cov_reg_design(object),
    "compound" = .cov_reg_compound(object))
```

With common sense

```R
# like this
return(.new_svystat_rob("mean", dat$yname,
    paste0("Huber M-estimator (type = ", type,
    if (asym) "; asym. psi" else "", ")"), match.call(),
    design, "mest", type = type, psi = if (asym) 1 else 0,
    psi_fun = "Huber", k = k))

# not
return(.new_svystat_rob("mean", dat$yname, paste0("Huber M-estimator
                                                  (type = ", type, if
                                                   (asym) "; asym. psi"
                                                 else "", ")"),
                                             match.call(), design,
                                             "mest", type = type, psi =
                                                 if (asym) 1 else 0,
                                             psi_fun = "Huber", k = k))
```

