# Model-Based Covariance Estimation for Regression *M*- and *GM*-Estimators

Tobias Schoch – April 2, 2022

[TOC]

## 1 Introduction

The *population* regression model is given by
$$
\begin{equation*}
	\xi: \quad Y_i = \bold x_i^T \boldsymbol{\theta} + \sigma \sqrt{v_i}E_i, \qquad \boldsymbol{\theta} \in \R^p, \quad \sigma > 0, \quad i \in U,
\end{equation*}
$$
where the population $U$ is of size $N$; the parameters $\boldsymbol{\theta}$ and $\sigma$ are unknown; the $\bold x_i$'s are known values (possibly containing outliers), $\bold x_i \in \R^p$, $1 \leq p < N$; the $v_i$'s are known positive (heteroscedasticity) constants; the errors $E_i$ are independent and identically distributed (i.i.d.) random variables with zero expectation and unit variance; it is assumed that $\sum_{i \in U} \bold x_i \bold x_i^T / v_i$ is a non-singular $(p \times p)$ matrix.

It is assumed that a sample $s$ is drawn from $U$ with sampling design $p(s)$ such that the independence structure of model $\xi$ is maintained. The sample regression $GM$-estimator of $\boldsymbol{\theta}$ is defined as the root to the estimating equation $\widehat{\boldsymbol{\Psi}}_n(\boldsymbol{\theta}, \sigma) = \bold 0$ (for all $\sigma > 0$), where

$$
\begin{equation*}
    \widehat{\boldsymbol{\Psi}}_n(\boldsymbol{\theta}, \sigma) = \sum_{i \in s} w_i \boldsymbol{\Psi}_i(\boldsymbol{\theta}, \sigma) \qquad \text{with} \qquad \boldsymbol{\Psi}_i(\boldsymbol{\theta}, \sigma) = \eta\left(\frac{y_i - \bold x_i^T \boldsymbol{\theta}}{\sigma \sqrt{v_i}}, \; \bold x_i\right) \frac{\bold x_i}{\sigma \sqrt{v_i}},
\end{equation*}
$$
where the function $\eta: \R \times \R^p \rightarrow \R$ parametrizes the following estimators

| M-estimator                 | Mallows GM-estimator                        | Schweppe GM-Estimator                                        |
| --------------------------- | ------------------------------------------- | ------------------------------------------------------------ |
| $\eta(r,\bold x) = \psi(r)$ | $\eta(r,\bold x) =\psi(r) \cdot h(\bold x)$ | $\eta(r,\bold x) = \psi\left(\displaystyle{\frac{r}{h(\bold x)}}\right) \cdot h(\bold x)$ |

where $\psi:\R \rightarrow \R$ is a continuous, bounded, and odd (possibly redescending) function, and $h: \R^p \rightarrow \R_+$ is a weight function.

## 2 Covariance estimation

The model-based covariance matrix of $\boldsymbol{\theta}$ is ([Hampel et al.](#literature), 1986. Chap. 6.3)
$$
\begin{equation}\label{eq:cov}
    \mathrm{cov}_{\xi}(\boldsymbol{\theta}, \sigma) = \bold M^{-1}(\boldsymbol{\theta}, \sigma) \cdot \bold Q(\boldsymbol{\theta}, \sigma) \cdot \bold M^{-T}(\boldsymbol{\theta}, \sigma) \qquad \text{for known} \; \sigma > 0,
\end{equation}
$$
where
$$
\begin{equation*}
    \bold M(\boldsymbol{\theta}, \sigma) = \sum_{i=1}^N \mathrm{E}_{\xi} \big\{ \boldsymbol{\Psi}_i'(\boldsymbol{\theta}, \sigma) \big\}, \quad \text{where} \quad
    \boldsymbol{\Psi}_i'(\boldsymbol{\theta}, \sigma) = -\frac{\partial}{\partial \boldsymbol{\theta}^*} \boldsymbol{\Psi}_i(Y_i, \bold x_i; \boldsymbol{\theta}^*, \sigma) \bigg\vert_{\boldsymbol{\theta}^* = \boldsymbol{\theta}},
\end{equation*}
$$

$$
\begin{equation*}
	\bold Q(\boldsymbol{\theta}, \sigma) = \frac{1}{N} \sum_{i=1}^N \mathrm{E}_{\xi} \big\{ \boldsymbol{\Psi}_i(Y_i, \bold x_i; \boldsymbol{\theta}, \sigma) \boldsymbol{\Psi}_i(Y_i, \bold x_i;\boldsymbol{\theta}, \sigma)^T \big\},
\end{equation*}
$$

and $\mathrm{E}_{\xi}$ denotes expectation with respect to model $\xi$. For the sample regression $GM$-estimator $\widehat{\boldsymbol{\theta}}_n$, the matrices $\bold M$ and $\bold Q$ must be estimated. In place of $\bold M$ and $\bold Q$ in ($\ref{eq:cov}$), we have

| M-estimator                                                  | Mallows GM-estimator                                         | Schweppe GM-Estimator                                        |
| ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| $\widehat{\bold M}_M = - \overline{\psi'} \cdot \bold X^T \bold W \bold X$ | $\widehat{\bold M}_{Mal} = - \overline{\psi'} \cdot \bold X^T \bold W \bold H \bold X$ | $\widehat{\bold M}_{Sch} = - \bold X^T \bold W \bold S_1 \bold X$ |
| $\widehat{\bold Q}_M = \overline{\psi^2} \cdot \bold X^T \bold W \bold X$ | $\widehat{\bold Q}_{Mal} = \overline{\psi^2} \cdot \bold X^T \bold W \bold H^2 \bold X$ | $\widehat{\bold Q}_{Sch} = \bold X^T \bold W \bold S_2 \bold X$ |

where
$$
\begin{equation*}
    \bold W = \mathrm{diag}_{i=1,\ldots,n}\{w_i\} \qquad \text{and} \qquad
    \bold H = \mathrm{diag}_{i=1,\ldots,n}\{h(\bold x_i)\},
\end{equation*}
$$

$$
\begin{equation*}
	\overline{\psi'} = \frac{1}{\widehat{N}}\sum_{i \in s} w_i \psi' \left( \frac{r_i}{\widehat{\sigma} \sqrt{v_i}} \right) \qquad \text{and} \qquad
\overline{\psi^2} = \frac{1}{\widehat{N}}\sum_{i \in s} w_i \psi^2 \left( \frac{r_i}{\widehat{\sigma} \sqrt{v_i}} \right),
\end{equation*}
$$

$$
\begin{equation*}
    \bold S_1 = \mathrm{diag}_{i=1,\ldots,n} \big\{ s_1^i \big\}, \quad \text{with} \quad s_1^i = \frac{1}{\widehat{N}}\sum_{j \in s} w_j \psi'\left(\frac{r_j}{h(\bold x_i)\widehat{\sigma} \sqrt{v_j}}\right),
\end{equation*}
$$

and
$$
\begin{equation*}
    \bold S_2 = \mathrm{diag}_{i=1,\ldots,n} \big\{ s_2^i \big\}, \quad \text{with} \quad s_2^i = \frac{1}{\widehat{N}}\sum_{j \in s} w_j \psi^2\left(\frac{r_j}{h(\bold x_i)\widehat{\sigma} \sqrt{v_j}}\right).
\end{equation*}
$$
**Remarks**

* The $i$-th diagonal element of $\bold S_1$ and $\bold S_2$ depends on $h(\bold x_i)$, but the summation is over $j \in s$; see also [Marazzi](#literature) (1987, Chap. 6).

* When $\bold W$ is equal to the identity matrix $\bold I$, the asymptotic covariance of $\widehat{\boldsymbol{\theta}}_M$ is equal to the expression in  [Huber](#literature) (1981, Eq. 6.5), which is implemented in the R packages `MASS` [Venables and Ripley](#literature) (2002) and `robeth` [Marazzi](#literature) (2020).

* For the Mallows and Schweppe type $GM$-estimators and given that $\bold W = \bold I$, the asymptotic covariance coincides with the one implemented in package/ library `robeth` for the option *averaged*; see [Marazzi](#literature) (1993, Chap. 4) and [Marazzi](#literature) (1987, Chap. 2.6) on the earlier `ROBETH-85` implementation.

## 3 Implementation

The main function—which is only a wrapper function— is `cov_reg_model` . The following display shows  pseudo code of the main function.

````c
cov_reg_model()
{
   get_psi_function() 			// get psi function (function pointer)
	get_psi_prime_function() 	// get psi-prime function (function pointer)
	switch(type) {
		case 0: cov_m_est()					// M-estimator
     	case 1: cov_mallows_gm_est() 		// Mallows GM-estimator
	 	case 2: cov_schweppe_gm_est() 	// Schweppe GM-estimator
	}
   robsurvey_error()			// terminate with error in case of failure
}
````

The functions `cov_m_est()`, `cov_mallows_gm_est()`, and `cov_schweppe_gm_est()` implement the covariance estimators; see below. All functions are based on the subroutines in `BLAS` ([Anderson et al.](#literature), 1999) and `LAPACK` ([Blackford et al.](#literature), 2002).

To fix notation, denote the Hadamard product of the matrices $\bold A$ and $\bold B$ by $\bold A\circ \bold B$ and suppose that $\sqrt{\cdot}$ is applied element by element.

**M-estimator** (`cov_m_est`). The covariance matrix is (up to a scalar) equal to
$$
\begin{equation}\label{eq:cov_m}
   (\bold X^T \bold W \bold X)^{-1}
\end{equation}
$$
and is computed as follows:

* Compute the factorization $\sqrt{\bold w} \circ \bold X := \bold Q \bold R$ (LAPACK: `dgeqrf`).
* Invert the upper triangular matrix $\bold R$ by backward substitution to get $\bold R^{-1}$ (LAPACK: `dtrtri`).
* Compute $\bold R^{-1} \bold R^{-T}$, which is equal to ($\ref{eq:cov_m}$); taking advantage of the triangular shape of $\bold R^{-1}$ and $\bold R^{-T}$ (LAPACK: `dtrmm`).

**Mallows GM-estimator** (`cov_mallows_gm_est`). The covariance matrix is (up to a scalar) equal to
$$
\begin{equation}\label{eq:cov_mallows}
   \big(\bold X^T \bold W \bold H \bold X\big)^{-1} \bold X^T \bold W \bold H^2 \bold X \big(\bold X^T \bold W \bold H \bold X\big)^{-1}
\end{equation}
$$
and is computed as follows:

* Compute the QR factorization: $\sqrt{\bold w \cdot \bold h} \circ \bold X := \bold Q \bold R$ (LAPACK: `dgeqrf`).
* Invert the upper triangular matrix $\bold R$ by backward substitution to get $\bold R^{-1}$ (LAPACK: `dtrtri`).
* Define a new matrix: $\bold A \leftarrow \sqrt{\bold h} \circ \bold Q$ (extraction of $\bold Q$ matrix with LAPACK: `dorgqr`).
* Update the matrix: $\bold A \leftarrow \bold A \bold R^{-T}$ (taking advantage of the triangular shape of $\bold R^{-1}$; LAPACK: `dtrmm`).
* Compute $\bold A \bold A^T$, which corresponds to the expression in ($\ref{eq:cov_mallows}$); (LAPACK: `dgemm`).

**Schweppe GM-estimator** (`cov_schweppe_gm_est`). The covariance matrix is (up to a scalar) equal to
$$
\begin{equation}\label{eq:cov_schweppe}
   \big(\bold X^T \bold W \bold S_1 \bold X\big)^{-1} \bold X^T \bold W \bold S_2 \bold X \big(\bold X^T \bold W \bold S_1 \bold X\big)^{-1}.
\end{equation}
$$
 Put $\bold s_1 = \mathrm{diag}(\bold S_1)$, $\bold s_2 = \mathrm{diag}(\bold S_2)$, and let $\cdot / \cdot $ denote elemental division (i.e., the inverse of the Hadamard product). The covariance matrix in ($\ref{eq:cov_schweppe}$) is computed as follows

* Compute the factorization $\sqrt{\bold w \circ \bold s_1} \circ \bold X := \bold Q \bold R$ (LAPACK: `dgeqrf`).
* Invert the upper triangular matrix $\bold R$ by backward substitution to get $\bold R^{-1}$ (LAPACK: `dtrtri`).
* Define a new matrix: $\bold A \leftarrow \sqrt{\bold s_2 / \bold s_1 } \circ \bold Q$ (extraction of $\bold Q$ matrix with LAPACK: `dorgqr`).
* Update the matrix: $\bold A \leftarrow \bold A \bold R^{-T}$ (taking advantage of the triangular shape of $\bold R^{-1}$; LAPACK: `dtrmm`).
* Compute $\bold A \bold A^T$, which corresponds to the expression in ($\ref{eq:cov_schweppe}$); (LAPACK: ` dgemm`).

**Note**. [Marazzi](#literature) (1987) uses the Cholesky factorization (see subroutines `RTASKV` and `RTASKW`) which is computationally a bit cheaper than our QR factorization.

## Literature

ANDERSON, E., BAI, Z., BISCHOF, C., BLACKFORD, L. S. , DEMMEL, J., DONGARRA, J., CROZ, J. D. , GREENHAUM, A., HAMMARLING, S., MCKENNEY, A., AND D. SORENSEN (1999). LAPACK Users’ Guide. 3rd edition.
Society for Industrial and Applied Mathematics (SIAM), Philadelphia.

BLACKFORD, L. S., PETITET, A., POZO, R., REMINGTON, K., WHALEY, R. C., DEMMEL, J., DONGARRA, J., DUFF,
I., HAMMARLING, S., HENRY, G., HEROUX, M., KAUFMAN, L. AND A. Lumsdaine (2002). An Updated
Set of Basic Linear Algebra Subprograms (BLAS). ACM Transactions on Mathematical
Software, 28, 135–151.

HAMPEL, F. R., E. M. RONCHETTI, P. J. ROUSSEEUW, AND W. A. STAHEL (1986). Robust Statistics: The Approach Based on Influence Functions, New York: John Wiley & Sons.

HUBER, P. J. (1981). Robust Statistics, New York: John Wiley & Sons.

MARAZZI, A. (1987). Subroutines for robust and bounded influence regression in ROBETH, Cahiers
de Recherches et de Documentation, 3 ROBETH 2, Division de Statistique et Informatique, Institut Universitaire de Médecine Sociale et Préventive, Lausanne, ROBETH-85 Document No. 2, August
1985, revised April 1987.

MARAZZI, A. (1993). Algorithms, Routines, and S Functions for Robust Statistics: The FORTRAN Library
ROBETH with an interface to S-PLUS, New York: Chapman & Hall.

MARAZZI, A. (2020). robeth: R Functions for Robust Statistics, R package version 2.7-6.

VENABLES, W. N. AND B. D. RIPLEY (2002). Modern Applied Statistics with S, New York: Springer, 4th ed.
