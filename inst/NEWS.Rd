\name{NEWS}
\title{News for \R Package \pkg{robsurvey}}
\encoding{UTF-8}
\section{CHANGES in robsurvey VERSION 0.3 (2022-04-02)}{
    \subsection{NEW FEATURES}{
        \itemize{
            \item robust generalized regression estimation (GREG) of the mean
                and total; see \code{svymean_reg} and \code{svytotal_reg}
            \item Diagnostic plots for regression, i.e., objects of class
                \code{svyreg_rob}
            \item Implementation of robust regression: If the estimated
                regression scale (by default weighted MAD) is zero (or nearly
                so), the weighted IQR is tried instead
        }
    }
    \subsection{CHANGES}{
        \itemize{
            \item The default functions for regression M-estimators are now
                called \code{svyreg_huberM} and \code{svyreg_tukeyM}; the
                old functions \code{svyreg_huber} and \code{svyreg_tukey}
                are deprecated but are kept for compatibility reasons.
            \item Documentation files in folder \code{inst/doc} have been
                updated.
            \item Several internal changes (e.g., default value of
                \code{k_Inf} is now \code{1e06} not \code{1e05}; see
                function \code{svyreg_control()}).
        }
    }
    \subsection{BUG FIXES}{
        \itemize{
            \item For designs with unequal probability sampling, the variance
                estimates of the robust estimators of mean and total
                are now identical with the estimates of
                \code{survey::svymean()} and \code{survey::svytotal()} if
                the tuning constant is \code{k = Inf} or \code{LB = 0}
                and \code{UB = 1}
        }
    }
}
\section{CHANGES in robsurvey VERSION 0.2 (2020-11-25)}{
    \subsection{NEW FEATURES}{
        \itemize{
            \item Weighted regression GM-estimators (e.g.,
                \code{svyreg_huberGM})
            \item M- and GM-estimator of regression with Tukey biweight
                psi-function (\code{svyreg_tukey} and \code{svyreg_tukeyGM})
            \item Weighted \code{k} winsorized mean and total (e.g.,
            \code{weighted_mean_k_winsorized} and \code{svymean_k_winsorized})
            \item Dalen's weight reduction estimator of the mean and total
                (e.g., \code{weighted_mean_dalen} and \code{svymean_dalen})
            \item Weighted M-estimator of the mean and total with Tukey
                biweight psi-function
            \item Weighted Huber Proposal 2 estimator (\code{huber2})
            \item Data sets \code{counties}, \code{flour}, \code{losdata}, and
                \code{MU284strat}
        }
    }
    \subsection{BUG FIXES}{
        \itemize{
            \item The original C implementation of \code{wquantile} was buggy
                (with implications for \code{weighted_quantile} and also the
                iterative re-weighted least squares algorithm). The new C
                implementation of of \code{wquantile} is sound.
        }
    }
    \subsection{CHANGES}{
        \itemize{
            \item Argument \code{type = "rwm"} of \code{weighted_mean_huber}
                is not used anymore (deprecated); instead, the type is now
                called \code{"rhj"}.
        }
    }
    \subsection{CHANGE of LICENSE and MAINTAINER}{
        \itemize{
            \item The authors agreed on November 25, 2020, to license/
                re-license the package version 0.2 under the GPL-2 resp.
                GPL-3 license.  (Version 0.1 has been licensed under the
                MIT license).
            \item Tobias Schoch is now the maintainer of the package.
        }
    }
}
