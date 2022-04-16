# Robust Generalized Regression Estimator

Tobias Schoch – April 15, 2022

[TOC]

## 1 Introduction

The *population* regression model is given by
$$
\begin{equation}\label{eq:reg_model}
	\xi: \quad Y_i = \bold x_i^T \boldsymbol{\theta} + \sigma \sqrt{v_i}E_i, \qquad \boldsymbol{\theta} \in \R^p, \quad \sigma > 0, \quad i \in U,
\end{equation}
$$
where the population $U$ is of size $N$; the parameters $\boldsymbol{\theta}$ and $\sigma$ are unknown; the $\bold x_i$'s are known values (possibly containing outliers), $\bold x_i \in \R^p$, $1 \leq p < N$; the $v_i$'s are known positive (heteroscedasticity) constants; the errors $E_i$ are independent and identically distributed (i.i.d.) random variables with zero expectation and unit variance; it is assumed that $\sum_{i \in U} \bold x_i \bold x_i^T / v_i$ is a non-singular $(p \times p)$ matrix.

It is assumed that a sample $s$ is drawn from $U$ with sampling design $p(s)$ such that the independence structure of the model in ($\ref{eq:reg_model}$) is maintained. 

The (sample) weighted least squares estimator of $\boldsymbol \theta$ is denoted by $\widehat{\boldsymbol \theta}_{n}^{LS}$. The sample regression $GM$-estimator of $\boldsymbol{\theta}$ is defined as the root to the estimating equation $\widehat{\boldsymbol{\Psi}}_n(\boldsymbol{\theta}, \sigma) = \bold 0$ (for all $\sigma > 0$), where
$$
\begin{equation*}
    \widehat{\boldsymbol{\Psi}}_n(\boldsymbol{\theta}, \sigma) = \sum_{i \in s} w_i \boldsymbol{\Psi}_i(\boldsymbol{\theta}, \sigma) \qquad \text{with} \qquad \boldsymbol{\Psi}_i(\boldsymbol{\theta}, \sigma) = \eta\left(\frac{y_i - \bold x_i^T \boldsymbol{\theta}}{\sigma \sqrt{v_i}}, \; \bold x_i\right) \frac{\bold x_i}{\sigma \sqrt{v_i}},
\end{equation*}
$$
where the function $\eta: \R \times \R^p \rightarrow \R$ parametrizes the estimators; see file `regression_covariance.html` for more details. Let $(\widehat{\boldsymbol \theta}_n,\widehat{\sigma}_n)$ denote the tuple of sample-based weighted $GM$-estimates of the (super-) population parameters $(\boldsymbol \theta, \sigma)$ or the finite population parameters $(\widehat{\boldsymbol \theta}_N, \widehat{\sigma}_N)$.

## 2 Generalized regression estimator

Consider the availability of auxiliary data $\bold x_i$, $i \in U$. The
goal is to estimate the $y$-total (or $y$-mean), where the study variable
$Y_i$ is supposed to be related to the $\bold x_i$'s by the heteroscedastic
(super-) population model in ($\ref{eq:reg_model}$). Under the
model, the generalized regression estimator (GREG) of the $y$-total is
[(Särndal et al.](#literature), 1992, Chap. 6.4)
$$
\begin{equation}\label{eq:greg}
    \widehat{t}_y^{\;reg} = \bold t_{\bold x}^T \widehat{\boldsymbol \theta}_n^{LS} +
        \sum_{i \in s} \frac{y_i - \bold x_i^T \widehat{\boldsymbol \theta}_n^{LS}}{\pi_i},
\end{equation}
$$

where the population $x$-total, $\bold t_{\bold x} = \sum_{i \in U} \bold x_i$, is
a known quantity.

**Special case**. A special case of $\widehat{t}_y^{\;reg}$ obtains
if the following condition holds.

> *CS condition.* The vector $(v_1, \ldots, v_N)^T$ of variance constants is in the *column space* of the
> design matrix $(\bold x_1^T, \ldots, \bold x_N^T)^T$. That is, there
> exists $\boldsymbol \lambda \in \R^p$ such that $v_i = \boldsymbol \lambda^T \boldsymbol x_i$
> for all $i \in U$.

For model specifications which satisfy the CS condition, the second
summand on the r.h.s. of ($\ref{eq:greg}$) vanishes ([Särndal et al.](#literature), 1992, Result 6.5.1). Thus, the GREG estimator takes the
simple projective form
$$
\begin{equation*}
	\widehat{t}_y^{\;reg} = \bold t_{\bold x}^T \widehat{\boldsymbol \theta}_{LS}.
\end{equation*}
$$

## 2 Robust generalized regression estimator

The GREG estimator is an asymptotically design unbiased
(ADU) estimator of $t_y$ whether the model in ($\ref{eq:reg_model}$)
is correctly specified or not
([Särndal et al.](#literature), 1992); hence it is
commonly said to be robust. However, the situation is more
intricate because we need to distinguish representative from nonrepresentative outliers ([Chambers](#literature), 1986). In the presence
of nonrepresentative outliers, $\widehat{\boldsymbol \theta}_{LS}$  tends
to be biased and very inefficient as an estimator of $\boldsymbol \theta$.
The same applies to $\widehat{t}_y^{\;reg}$ as an estimator of $t_y$.
In the presence of representative outliers, the GREG
estimator owns the ADUness property and is said bias-robust. But it
is not variance-robust because representative outliers tend to
inflate the design variance ([Hulliger](#literature), 1995; [Beaumont and Alavi](#literature), 2004).

Robust
$$
\begin{equation}\label{eq:rob_greg}
    \widehat{t}_y^{\;rob} = \bold t_{\bold x}^T \widehat{\boldsymbol \theta}_{n} + \sum_{i \in s} d_i (y_i - \bold x_i^T \widehat{\boldsymbol \theta}_{n}),
\end{equation}
$$
where the $d_i$'s are real-valued constants. Six sets of constants are available.

* `type = "projective"`. This $d_i \equiv 0$

* `type = "ADU"`.  $d_i \equiv 1 / \pi_i$ 

* `type = "robust"` $d_i \equiv w_i^*(k)$ 

* `type = "BC"` $d_i \equiv w_i^*(k_p)$ 

* `type = "Lee"`$d_i \equiv \vartheta / \pi_i$ 

* `type = "Duchesne"`

$$
\begin{equation*}
d_i \equiv \frac{h(\boldsymbol x_i)\psi_{mod}(r_i)}{\pi_i v_i r_i h_1(\bold x_i)}
\end{equation*}
$$











## Literature

BEAUMONT, J.-F. AND A. ALAVI (2004). Robust Generalized Regression Estimation, *Survey Methodology* **30**, 195–208.

BEAUMONT, J.-F. AND L.-P. RIVEST (2009). Dealing with outliers in survey data, in *Sample Surveys: Theory, Methods and Inference*, ed. by D. Pfeffermann and C. R. Rao, Amsterdam: Elsevier, vol. 29A of *Handbook of Statistics*, chap. 11, 247–280.

CHAMBERS, R. (1986). Outlier Robust Finite Population Estimation, *Journal of the American Statistical Association* **81**, 1063–1069.

DUCHESNE, P. (1999). Robust calibration estimators, *Survey Methodology* **25**, 43–56.

LEE, H. (1995a). Outliers in business surveys, in *Business survey methods*, ed. by B. G. Cox, D. A. Binder, B. N. Chinnappa, A. Christianson, M. J. Colledge, and P. S. Kott, New York: John Wiley and Sons, chap. 26, 503–526.

SÄRNDAL, C.-E., B. SWENSSON, AND J. WRETMAN (1992). Model Assisted Survey Sampling, New York: Springer.

