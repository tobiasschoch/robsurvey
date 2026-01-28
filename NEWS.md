**If you are viewing this file on CRAN, please check [latest news on GitHub](https://github.com/tobiasschoch/robsurvey/blob/master/NEWS.md) where the formatting is also better.**

# robsurvey VERSION 0.7-1 (2026-01-28)

## NEW FEATURES

Added function `within_tolerance()` which returns the value `TRUE` if an observation is within tolerance limits. Observations that fall outside are declared (potential) outliers.

Added a `confint()` method for objects of class `svystat_rob`.

## CHANGES

The (internal) handling of subdomains has been changed; this should go unnoticed by the user.

# robsurvey VERSION 0.7 (2024-08-22)

## NEW FEATURES

* The survey methods with prefix `svymean_` or `svytotal_` can now be used with the function `survey::svyby()`.
* Method `mse()` is now also available for objects of class `svystat` (which is defined in pkg `survey`).
* Functions to draw pps samples (without replacement) by Brewer's method.

# robsurvey VERSION 0.6 (2024-01-14)

This release of the package fixes a number of issues that were brought to our attention by an anonymous reviewer. It also fixes some other problems.

## BUG FIX

Function `huber2()` erroneously returned `NA` if the initial scale estimate was less than `DBL_EPSILON` (i.e., zero in floating-point arithmetics) although not zero in real arithmetics. After removing the unnecessary  `scale < DBL_EPSILON` check, the function behaves as expected.

## CHANGES

* `R-core` has been added to the list of intellectual property owners with the roles `c("cph", ctb")` for `zeroin2.c` (see `DESCRIPTION` file).
* Automated tests have been part of the [GitHub](https://github.com/tobiasschoch/robsurvey) repo since the first commit. However, the tests have not been shipped with the source package. Now, the tests are included in the source package.
* Function `robsvyreg()` is not exported to the namespace anymore because it is regarded as an internal function (not to be called by users).
* Since **version 4.2** of the **survey** package (released in May 2023), the `survey` package allows the definition of pre-calibrated weights (see argument `calibrate.formula` of the function `survey::svydesign()`; see also vignette [Pre-calibrated weights](https://CRAN.R-project.org/package=survey/vignettes/precalibrated.pdf) of the `survey` package). From now on, we will use this functionality by default in the examples, vignettes and documentation. Our code automatically reverts/falls back to calling `svydesign()` without pre-calibrated weights (**legacy mode**) on R installations with an **earlier** version of the `survey` package. As a consequence, some of the variance and standard error estimates in legacy mode may differ from those with pre-calibrated weights.

## MISC

* Many of the help files have been updated and expanded.
* Since version 0.3, the functions `svyreg_huber()` and `svyreg_tukey()` are deprecated but have been kept for compatibility reasons. The deprecated functions are now equipped with a call to `.Deprecated()` and are documented separately in `help("robsurvey-deprecated")`.
* Some of the R and C code has been cleaned up; the `NAMESPACE` file has been consolidated; symbols of shared objects are now registered using `src/init.c` (this file has previously been called `robsurvey_init.c`) and loaded by `useDynLib()` as R objects (not character strings), whose names are pre-fixed by `"C_"` for reasons of transparency.
* The stratum variable `strat` in the dataset `workplace` is now a `factor`; the same applies to the variables `REG`, `CL` and `Stratum` in the datasets `MU284strat` and `MU284pps`.

# robsurvey VERSION 0.5-2 (2022-12-04)

## BUG FIX

Fixed a bug in the C function `wquant0`. For the special case of samples of size 2, the weighted quantile (other than the median) was wrong if the data were sorted in descending order. (Thanks to Ryota Suzuki, who detected the bug, [Issue #1](https://github.com/tobiasschoch/robsurvey/issues/1)).

# robsurvey VERSION 0.5-1 (2022-11-17)

## CHANGE

The `summary()` method for objects of class `formula` has been replaced by `svysummary()` because it did not handle non-standard cases correctly (Thanks to the editorial office of the Journal of Statistical Software for pointing this out).

## MISC

* Fixed defunct links in vignettes and added `requireNamespace()` as a guard for suggested packages in the vignettes.
* In the help files of functions that depend on the `survey` package or that extend its functionality, we added the following note: "Package `survey` must be loaded in order to use this function." to the Details section.

# robsurvey VERSION 0.5 (2022-10-07)

## CHANGES

* Files in `/doc` folder are now in `*.pdf` format which takes less space compared with the `*.html` format. Thus, the warning `checking installed package size... NOTE installed size is 5.4Mb sub-directories of 1Mb or more: doc 4.9Mb` disappeared
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
  * not be directly supported by the developers in the event issues arise


## CHANGES

* The default functions for regression M-estimators are now called `svyreg_huberM()` and `svyreg_tukeyM()`; the old functions `svyreg_huber()` and `svyreg_tukey()` are deprecated but are kept for compatibility reasons.
* Documentation files in folder `inst/doc` and test cases in `tests` have been updated
* Several internal changes (e.g., default value of `k_Inf` is now `1e06` not `1e05`; see function `svyreg_control()`).

## BUG FIX

For designs with unequal probability sampling, the variance estimates of the robust estimators of mean and total are now identical with the estimates of `survey::svymean()` and `survey::svytotal()` if the tuning constant is `k = Inf` or `LB = 0` and `UB = 1`.

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
