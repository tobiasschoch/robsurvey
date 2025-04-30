# Handling of NA and domains

Tobias Schoch – 31-01-2025 

## 1 Intro



## 2 Basic robust estimators

### 2.1 Bare-bone functions (`R/robsurvey-interal.R`)

`.check_data_weights(x, w, na.rm = FALSE, check_NA = TRUE)`, file: `robsurvey-interal.R`

* if `na.rm = TRUE`, missing values are removed using `complete.cases()`
* observations with zero weights kept untreated



### 2.2 Survey methods (`R/robsurvey-interal.R`)

`.check_formula(f, design, na.rm = FALSE, check_NA = TRUE)`, file: `robsurvey-interal.R`

* if `na.rm = TRUE`, missing values are removed using `complete.cases()`
  * the design is subset to complete cases
* domains are identified with $w_i = 0$
  * flag: `domain = TRUE`  
  * `in_domain` is a Boolean vector of length $n$, respectively, `in_domain = NULL` if `domain = FALSE`



**The prototype of a call** (for the example `svymean_trimmed()`; unimportant code has been left away)

```R
    dat <- .check_formula(x, design, na.rm)
```

* if `na.rm = TRUE` in the call of `svymean_`, the argument `na.rm = TRUE` is passed on to `.check_formula()` and subsequently missing values are removed using `complete.cases()`
* also, the design object in `dat$design` is a subset using `complete.cases()`

* domains are identified with $w_i = 0$

  * if any $w_i=0$, the flag `domain = TRUE` is set

  * The observations of the domains are identified using `in_domain`, which is a Boolean vector of length $n$ (without domains: `in_domain = NULL`)



Population- vs. domain-level estimate using the `dat` object

```R
    res <- if (dat$domain)
        weighted_mean_trimmed(dat$y[dat$in_domain], dat$w[dat$in_domain], LB,
                              UB, info = TRUE, na.rm = FALSE)
    else
        weighted_mean_trimmed(dat$y, dat$w, LB, UB, info = TRUE, na.rm = FALSE)
```

* `weighted_mean_trimmed()` is called with `na.rm = FALSE` because the NA treatment took place earlier



Variance estimation

```R
    # influence function
    infl <- .infl_trimmed(res$model$y, res$model$w, LB, UB, res$estimate) *
                res$model$w / sum(res$model$w)
    if (dat$domain) {
        tmp <- numeric(dat$n)
        tmp[dat$in_domain] <- infl
        infl <- tmp
    }
    # variance
    design <- dat$design
    res$variance <- svyrecvar(infl, design$cluster, design$strata, design$fpc,
                              postStrata = design$postStrata)
```

* The influence function is computed using the `res` object, which results from `weighted_mean_trimmed()`, which is the result of `.check_data_weights()`
* `design` is taken from the `dat` object



==ISSUE== With NA present in the data, `infl` is vector of length $\vert r\vert$ (response set), whereas all vectors of `dat$design` are of length $\vert s\vert$ (sample).



```R
	# return
	res$estimator$domain <- dat$domain
```

* the domain flag is added to the returned list



## 3 Regression

