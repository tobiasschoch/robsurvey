\docType{package}
\name{robsurvey-package}
\alias{robsurvey-package}
\alias{robsurvey}
\title{Package Overview}
\description{
    A key \emph{design pattern} of the package is that the majority
    of the estimating methods is available in two "flavors":
    \itemize{
        \item bare-bone methods
        \item survey methods
    }
    Bare-bone methods are stripped-down versions of the survey methods in
    terms of functionality and informativeness. These functions may serve
    users and package developers as building blocks. In particular,
    bare-bone functions \emph{cannot compute} variances.

    The survey methods are much more capable and depend, for variance
    estimation, on the \pkg{survey} package.
}
\section{Basic Robust Estimators}{
    \subsection{Trimming}{
        \itemize{
            \item Bare-bone methods: \code{\link{weighted_mean_trimmed}} and
                \code{\link{weighted_total_trimmed}}
            \item Survey methods: \code{\link{svymean_trimmed}} and
                \code{\link{svytotal_trimmed}}
        }
    }
    \subsection{Winsorization}{
        \itemize{
            \item Bare-bone methods:
            \itemize{
                \item \code{\link{weighted_mean_winsorized}} and
                    \code{\link{weighted_total_winsorized}}
                \item \code{\link{weighted_mean_k_winsorized}} and
                    \code{\link{weighted_total_k_winsorized}}
            }
            \item Survey methods:
            \itemize{
                \item \code{\link{svymean_winsorized}} and
                    \code{\link{svytotal_winsorized}}
                \item \code{\link{svymean_k_winsorized}} and
                    \code{\link{svytotal_k_winsorized}}
            }
        }
    }
    \subsection{Dalen's estimators (weight reduction methods)}{
        \itemize{
            \item Bare-bone methods: \code{\link{weighted_mean_dalen}} and
                \code{\link{weighted_total_dalen}}
            \item Survey methods: \code{\link{svymean_dalen}} and
                \code{\link{svytotal_dalen}}
        }
    }
    \subsection{M-estimators}{
        \itemize{
            \item Bare-bone methods:
            \itemize{
                \item \code{\link{weighted_mean_huber}} and
                    \code{\link{weighted_total_huber}}
                \item \code{\link{weighted_mean_tukey}} and
                    \code{\link{weighted_total_tukey}}
                \item \code{\link{huber2}} (weighted Huber Proposal 2
                    estimator)
            }
            \item Survey methods:
            \itemize{
                \item \code{\link{svymean_huber}} and
                    \code{\link{svytotal_huber}}
                \item \code{\link{svymean_tukey}} and
                    \code{\link{svytotal_tukey}}
                \item \code{\link{mer}} (minimum estimated risk estimator)
            }
        }
    }
}
\section{Survey Regression (weighted least squares)}{
    \code{\link{svyreg}}
}
\section{Robust Regression and Ratio Estimation (weighted)}{
    \itemize{
        \item Regression M-estimators: \code{\link{svyreg_huberM}} and
            \code{\link{svyreg_tukeyM}}
        \item Regression GM-estimators (Mallows and Schweppe):
            \code{\link{svyreg_huberGM}} and \code{\link{svyreg_tukeyGM}}
        \item Ratio M-estimators:
            \code{\link{svyratio_huber}} and \code{\link{svyratio_tukey}}

    }
}
\section{Robust Generalized Regression (GREG) and Ratio Prediction of the
        Population Mean and Total}{
    \itemize{
        \item Regression predictors: \code{\link{svymean_reg}} and
            \code{\link{svytotal_reg}}
        \item Ratio predictors: \code{\link{svymean_ratio}} and
            \code{\link{svytotal_ratio}}
    }
}
\section{Utility functions}{
    \itemize{
        \item \code{\link{weighted_quantile}} and \code{\link{weighted_median}}
        \item \code{\link{weighted_mad}} and \code{\link{weighted_IQR}}
        \item \code{\link{weighted_mean}} and \code{\link{weighted_total}}
        \item \code{\link{weighted_line}}, \code{\link{weighted_median_line}},
            and \code{\link{weighted_median_ratio}}
    }
}
