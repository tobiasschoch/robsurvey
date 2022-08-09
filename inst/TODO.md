# TODO - Road Map

---

## MILESTONE ver. 0.4

### DOCUMENTATION

* `packge.Rd`link to \pkg{survey} does not work
* Improve documentation of `weighted_line()`  and  `weighted_median_lines()`
* Documentation, examples, and vignette of `svymean_reg()` and `svytotal_reg()`

### NEW FEATURE

* `svyreg()`: predict method
* robust calibration

## BUG FIX

* GREG tests + heteroscedasticity (of calibrate)
* Mallows type GM-estimator



```R
library("survey")
library("robsurvey")
# qi's in Wright's QR-estimator
.qi_weights <- function(object)
{
    # robustness and design weights
    qi <- robweights(object) * object$model$w
    # account for heteroscedasticity
    if (!is.null(object$model$var))
        qi <- qi / object$model$var
    # multiplicative factor for Mallows type GM-estimator
    if (object$estimator$type == 1)
        qi <- qi * object$model$xwgt
    qi
}

dataset <- "MU284strat"
dn <- svydesign(ids = ~LABEL, strata = ~Stratum, fpc = ~fpc,
         weights = ~weights, data = MU284strat)
f <- RMT85 ~ P85 + S82 + CS82
set.seed(1)
xwgt <- runif(60, 0.8, 1)
m <- svyreg_huberGM(f, dn, k = 3, xwgt = xwgt, type = "Mallows")
coef(m)
lm.wfit(m$model$x, m$model$y, .qi_weights(m))$coefficients

```







---

## UNSCHEDULED
