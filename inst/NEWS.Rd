\name{NEWS}
\title{News for \R Package \pkg{robsurvey}}
\encoding{UTF-8}

\section{CHANGES in robsurvey VERSION 0.2 (2020-11-25)}{
\subsection{NEW FEATURES}{
\itemize{
    \item Sample weighted regression GM-estimators (e.g.,
    \code{svyreg_huberGM})
    \item M- and GM-estimator of regression with Tukey biweight
    psi-function (\code{svyreg_tukey} and \code{svyreg_tukeyGM})
    \item (Robust) generalized regression estimation (GREG and robust GREG)
    \item Weighted \code{k} winsorized mean and total (e.g.,
    \code{weighted_mean_k_winsorized})
    \item Dalen's weight reduction estimator of the mean and total (e.g.,
    \code{weighted_mean_dalen})
    \item Weighted M-estimator of the mean and total with Tukey biweight
    psi-function
    \item Weighted Huber Proposal 2 estimator (\code{huber2})
    \item Data sets \code{flour} and \code{losdata}
}
}
\subsection{BUG FIXES}{
\itemize{
    \item The original C implementation of \code{wquantile} was buggy (with
    implications for \code{weighted_quantile} and also the iterative
    re-weighted least squares algorithm). The new C implementation of
    of \code{wquantile} is sound.
}
}
\subsection{CHANGES}{
\itemize{
    \item Argument \code{type = "rwm"} of \code{weighted_mean_huber} is not
    used anymore (deprecated); instead, the type is now called \code{"rhj"}.
}
}
\subsection{CHANGE of LICENSE and MAINTAINER}{
\itemize{
    \item The authors agreed on November 25, 2020, to license/ re-license
    the package version 0.2 under the GPL-2 resp. GPL-3 license.
    (Version 0.1 has been licensed under the MIT license).
    \item Tobias Schoch is now the maintainer of the package.
}
}
}
