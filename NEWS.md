**If you are viewing this file on CRAN, please check [latest news on GitHub](https://github.com/tobiasschoch/robsurvey/blob/master/NEWS.md) where the formatting is also better.**

# robsurvey VERSION 0.5 (2022-10-07)

## CHANGES

* Files in `/doc` folder are now in `*.pdf` format which takes less space compared with the `*.html` format. Thus, the warning `checking installed package size ... NOTE
   installed size is 5.4Mb
   sub-directories of 1Mb or more:
    doc 4.9Mb` disappeared

* The print method for objects of class `svystat_rob` now correctly prints: [Estimator] `of the population` [mean/total].

* The default value of argument `type` in the functions `weighted_mean_huber()`, `weighted_mean_tukey()`, `svymean_huber()` and `svymean_tukey()` is now `"rwm"`. Type `"rhj"` is still available (and will be supported in the future) but is silently converted to `"rwm"`.

# robsurvey VERSION 0.4 (2022-09-08)

## NEW FEATURES

* Robust estimators of the ratio of two variables (`svyratio_huber()` and `svyratio_tukey()`); these functions are robust alternatives to `survey::svyratio()`.
* Robust ratio estimators of the population mean and total, see `svymean_ratio()` and `svytotal_ratio()`.
* Example data `MU284pps`: A pps sample without replacement of size 50 from the MU284 population in Särndal et al. (1992).

## CHANGES

* Functions `svymean_reg()` and `svytotal_reg()` are not flagged as "experimental" anymore. Several changes took place (in fact, the functions have undergone a complete code refactoring):
  * Argument `auxiliary` has been replaced by the two arguments `N` (population size) and `totals` (i.e., population totals of non-constant explanatory variables). **Important:** `svymean_reg()` is now called with `totals` not the population means.
  * The arguments `na.rm` and `verbose` have been dropped (not needed).
  * The variance estimators in `svymean_reg()` and `svytotal_reg()` are now implemented as *g*-weighted residual variance estimators.


* Added documentation for variable `strat` in the `workplace` data and updated description of variable `payroll`.

* Added 45-degree line in the diagnostic `plot` method for "3 Response vs. Fitted values" (`which = 3`) of class `svyreg_rob`.
* Method `SE()` for class `svyreg_rob` is now exported to the namespace.

## BUG FIX

* Slot `estimator$string` in the return value of function `mer()` indicates the name of the underlying estimator correctly.
* Fixed annotation of observations in diagnostic plot "Sqrt of abs(Residuals) vs. Fitted values" (`which = 4` in `plot`) for class `svyreg_rob`.

# robsurvey VERSION 0.3 (2022-06-04)

## NEW FEATURES

* Diagnostic plots for fitted regression model, i.e., objects of class `svyreg_rob`
* Robust regression: If the estimated regression scale (by default weighted MAD) is zero (or nearly so), the weighted IQR is tried instead. If the weighted IQR ist also zero, the function returns with an error.
* Function `mer()` for minimum estimated risk estimation of location gained two new arguments:
  * `method`: the method used in the search for a minimum, e.g., `"Brent"`, `"BFGS"`, see `stats::optim()` for more details
  * `init` determines the left side of the search interval and the initial value in the minimization approach

* Function `mse()` computes/ extracts the estimated mean square error/ estimated risk in presence of representative outliers; see also `mer()`

* Robust generalized regression estimation (GREG) of the mean and total; see `svymean_reg()` and `svytotal_reg()`. The current implementation of the functions is **EXPERIMENTAL** and a warning is issued when calling the functions (unless `verbose = FALSE`). Experimental features may:
  * have undergone less extensive testing than is normal for standard features
  * interact with unstable (external) dependencies
  * be subject to change
  * not be directly supported by the developers in the event
    issues arise
  


## CHANGES

* The default functions for regression M-estimators are now called `svyreg_huberM()` and `svyreg_tukeyM()`; the old functions `svyreg_huber()` and `svyreg_tukey()` are deprecated but are kept for compatibility reasons.
* Documentation files in folder `inst/doc` and test cases in `tests` have been updated
* Several internal changes (e.g., default value of `k_Inf` is now `1e06` not `1e05`; see function `svyreg_control()`).

## BUG FIX

For designs with unequal probability sampling, the variance estimates of the robust estimators of mean and total are now identical with the estimates of `survey::svymean()` and `survey::svytotal()` if the tuning constant is `k = Inf` or `LB = 0` and `UB = 1`

## MISC

Added `DOI` to all references (where available).

# robsurvey VERSION 0.2 (2022-01-17)

## NEW FEATURES

* Weighted regression GM-estimators; see e.g., `svyreg_huberGM()`
* M- and GM-estimator of regression with Tukey biweight psi-function; see `svyreg_tukey()` and `svyreg_tukeyGM()`
* Weighted `k` winsorized mean and total; see `weighted_mean_k_winsorized()` and `svymean_k_winsorized()`
* Dalen's weight reduction estimator of the mean and total; see `weighted_mean_dalen()` and `svymean_dalen()`
* Weighted M-estimator of the mean and total with Tukey biweight psi-function
* Weighted Huber Proposal 2 estimator; see `huber2()`
* Data sets `counties`, `flour`, `losdata`, and `MU284strat`

## BUG FIX

The original C implementation of `wquantile` was buggy (with implications for the R function `weighted_quantile()` and also the iterative re-weighted least squares algorithm). The new C implementation of `wquantile` is sound.

## CHANGE

Argument `type = "rwm"` of `weighted_mean_huber()` is not used anymore (deprecated); instead, the type is now called `"rhj"`.

## CHANGE of LICENSE and MAINTAINER

* The authors agreed on November 25, 2020, to license/ re-license the package version 0.2 under the GPL-2 resp.  GPL-3 license. (Version 0.1 has been licensed under the MIT license).
* Tobias Schoch is now the maintainer of the package.
