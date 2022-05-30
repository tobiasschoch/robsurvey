# robsurvey<img src="inst/varia/logo.svg" align="right" width=120 height=139 alt="" />

<!-- badges: start -->
[![CRAN](https://www.r-pkg.org/badges/version/robsurvey)](https://cran.r-project.org/package=robsurvey)
<!-- badges: end -->

## Summary

Functions to compute robust (outlier-resistant) estimates of
finite population characteristics. The package supports the computations of robust means, totals, ratios, etc. Available methods are regression
*M*- and *GM*-estimators, trimming, and winsorization, etc. The package complements the [survey](https://cran.r-project.org/package=survey) package.

## 1 What the package offers

### 1.1 Basic Robust Estimators

A key *design pattern* of the package is that the majority of the estimating methods is available in two "flavors":

- bare-bone methods
- survey methods

Bare-bone methods are stripped-down versions of the survey methods in terms of functionality and informativeness. These functions may serve users and other package developers as building blocks. In particular, bare-bone functions *cannot compute* variances. The survey methods are much more capable and depend, for variance estimation, on the [survey](https://CRAN.R-project.org/package=survey) package.

**Trimming**

| Bare-bone methods          | Survey methods       |
| -------------------------- | -------------------- |
| `weighted_mean_trimmed()`  | `svymean_trimmed()`  |
| `weighted_total_trimmed()` | `svytotal_trimmed()` |

**Winsorization**

| Bare-bone methods               | Survey methods            |
| ------------------------------- | ------------------------- |
| `weighted_mean_winsorized()`    | `svymean_winsorized()`    |
| `weighted_mean_k_winsorized()`  | `svymean_k_winsorized()`  |
| `weighted_total_winsorized()`   | `svytotal_winsorized()`   |
| `weighted_total_k_winsorized()` | `svytotal_k_winsorized()` |

**Dalen's estimators (weight reduction methods)**

| Bare-bone methods        | Survey methods     |
| ------------------------ | ------------------ |
| `weighted_mean_dalen()`  | `svymean_dalen()`  |
| `weighted_total_dalen()` | `svytotal_dalen()` |

**M-estimators**

| Bare-bone methods        | Survey methods     |
| ------------------------ | ------------------ |
| `weighted_mean_huber()`  | `svymean_huber()`  |
| `weighted_mean_tukey()`  | `svymean_tukey()`  |
| `weighted_total_huber()` | `svytotal_huber()` |
| `weighted_total_tukey()` | `svytotal_tukey()` |

The *M*-estimators have a `type` argument taking the values `"rht"` or `"rhj"` to specify, respectively, the robust Horvitz–Thompson (RHT) or the robust Hajek (RHJ) estimator. 

In addition, `huber2()` implements a weighted Huber "Proposal 2" estimator (only bare-bone function). Function `mer()` (minimum estimated risk estimator) implements an adaptive *M*-estimator.

**Utility functions**

- `weighted_quantile()` and `weighted_median()`
- `weighted_mad()` and `weighted_IQR()`
- `weighted_mean()` and `weighted_total()`

### 1.2 Robust weighted regression

**Weighted least squares**

* `svyreg()`

**Weighted regression *M*-estimator**

* `svyreg_huberM()`
* `svyreg_tukeyM()`

**Weighted regression *GM*-estimator** (Mallows and Schweppe type)

* `svyreg_huberGM()`
* `svyreg_tukeyGM()`

### 1.3 Robust generalized regression estimator (GREG)

> The functions are experimental; they may change!

* `svymean_reg()`
* `svytotal_reg()`

### 1.4 Weighted resistant line

- `weighted_line()`
- `weighted_median_line()`
- `weighted_median_ratio()`

## 2 Installation

The package can be installed from CRAN using
```
install.packages("robsurvey")
```

## 3 Building

Make sure that the R package `devtools` is installed. Then, the `robsurvey` package can be pulled from this GitHub repository and installed by
```
devtools::install_github("tobiasschoch/robsurvey")
```

## 4 Community guidelines

#### Submitting an issue

If you have any suggestions for feature additions or any problems with the software that you would like addressed with the development community, please submit an issue on the Issues tab of the project GitHub repository. You may want to search the existing issues before submitting, to avoid asking a question or requesting a feature that has already been discussed.

#### How to contribute

If you are interested in modifying the code, you may fork the project for your own use, as detailed in the FreeBSD License we have adopted for the project. In order to contribute, please contact the developer by Tobias Schoch at gmail dot com (the names are separated by a dot) after making the desired changes.

#### Asking for help

If you have questions about how to use the software, or would like to seek out collaborations related to this project, you may contact Tobias Schoch (see contact details above).
