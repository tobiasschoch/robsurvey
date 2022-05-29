# Robust Generalized Regression Estimator

Tobias Schoch

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

The (sample) weighted least squares estimator of $\boldsymbol \theta$ is denoted by $\widehat{\boldsymbol \theta}_{n}^{LS}$. The sample regression *M*- and *GM*-estimator of $\boldsymbol{\theta}$ are defined as the root to the following estimating equations ([Hampel et al.](#literature), 1986, Chap. 6.3)
$$
\begin{align*}
	\sum_{i \in s} \frac{w_i}{\sqrt{v_i}} \psi_k (r_i)\bold x_i &= \bold 0 & \text{($M$-estimator)},\\
	\sum_{i \in s} \frac{w_i}{\sqrt{v_i}} h(\bold x_i) \psi_k (r_i) \bold x_i &= \bold 0 & \text{(Mallows $GM$-estimator)},\\
	\sum_{i \in s} \frac{w_i}{\sqrt{v_i}} h(\bold x_i) \psi_k\left(\frac{r_i}{h(\bold x_i)}\right) \bold x_i &= \bold 0& \text{(Schweppe $GM$-estimator)},
\end{align*}
$$
where $\psi_k$ is a $\psi$-function indexed by the robustness tuning constant $k$, $r_i = (y_i - \bold x_i^T\boldsymbol \theta) / (\sigma\sqrt{v_i})$ is the standardized residual, $w_i$ denotes the sampling weight, and $h:\R^p \rightarrow \R_+ $ is a weight function.

Let $(\widehat{\boldsymbol \theta}_n,\widehat{\sigma}_n)$ denote the tuple of estimates The estimators can be written as a weighted least squares problem
$$
\begin{equation*}
\sum_{i \in s} \frac{w_i}{v_i} u_i (y_i - \bold x_i^T \widehat{\boldsymbol \theta}_n) \bold x_i = \bold 0,
\end{equation*}
$$
where the $u_i$'s are defined in the following table.

| M-estimator               | Mallows GM-estimator                  | Schweppe GM-estimator                                        |
| ------------------------- | ------------------------------------- | ------------------------------------------------------------ |
| $u_i = \psi_k(r_i) / r_i$ | $u_i = h(\bold x_i)\psi_k(r_i) / r_i$ | $u_i = \psi_k(r_i^*) / r_i^*$, where $r_i^* = r_i / h(\bold x_i)$ |


## 2 Robust generalized regression estimator

Following [Duchsene](#literature) (1999) the robust GREG estimator of the $y$-total can be written as a QR-estimator ([Wright](#literature), 1983)
$$
\begin{equation*}
    \widehat{t}_y^{\,rob} = \bold t_{\bold x}^T \widehat{\boldsymbol \theta}_{n} + \sum_{i \in s} b_i (y_i - \bold x_i^T \widehat{\boldsymbol \theta}_{n}),
\end{equation*}
$$
where $\widehat{\boldsymbol \theta}_n$ is a weighted regression *GM*-estimate, $\bold t_x$ is the population $x$-total, and the $b_i$'s are real-valued constants. Six sets of constants $\{b_i\}$ are available that define different types of estimators.

**`type = "projective"`** ignores the bias correction term. Therefore, $\widehat{t}_y^{\,rob}$ is equal to the simple projective estimator. It obtains for $b_i \equiv 0$. The projective estimator might be very biased; see [Gwet and Rivest](#literature) (1992).

**`type = "ADU"`** yields an estimator that is asymptotically design unbiased (ADU) but not robust. This estimator reflects the finding of [Gwet and Rivest](#literature) (1992) that an estimator which is asymptotically design unbiased cannot be resistant to outliers (and vice versa). This estimator obtains for $b_i \equiv w_i / v_i$.

**`type = "robust"`** gives a robust GREG estimator without bias correction. It obtains for  $b_i \equiv u_i w_i / v_i$. The choice of $k$ in $\psi_k$ can be different from the one used in the computation of $\widehat{\boldsymbol \theta}_n$.

**`type = "lee"`** has been proposed by [Lee](#literature) (1991, 1995). He observed that the estimator has an asymptotic bias of $(k - 1)\sum_{i \in U} r_i(\boldsymbol \theta_N)$, where $r_i(\boldsymbol \theta_N)$ is the population residual that obtains by fitting a robust regression model for the entire population. By choosing $0 \leq k \leq 1$, the statistician can control the bias of the estimator. The estimator obtains for $b_i \equiv k w_i/v_i$.

**`type = "BR"`** of [Beaumont and Rivest](#literature) (2009) is a bias-corrected estimator. It obtains for $b_i \equiv u_i w_i/v_i$, where $u_i$ is defined with $\psi_k$ replaced by 
$$
\begin{equation*}
	\psi_k^*(x) = \frac{x}{w_i} + \frac{w_i - 1}{w_i}\psi_k^{hub}(x),
\end{equation*}
$$
where $\psi_k^{hub}$ is the Huber $\psi$-function with tuning constant $k$. The resulting estimator $\widehat{t}_y^{\;rob}$ is census consistent, i.e., it is equal to $t_y$ if $s = U$.

**`type = "duchesne"`** has been proposed by [Duchesne](#literature) (1999). It is a bias-corrected estimator and obtains for $b_i \equiv u_i w_i / v_i$, where $u_i$ is defined with $\psi_k$ replaced by
$$
\begin{equation*}
	\psi_{a,b}^*(x) = \begin{cases} x & \text{if} \vert x \vert \leq a, \\ a \cdot \mathrm{sign}(x) & \text{if} \vert x \vert > a \; \text{and} \; \vert x \vert < a/b, \\ b \cdot x & \text{if} \vert x \vert > a/b, 
        \end{cases}
\end{equation*}
$$

where $\psi_{a,b}$ is a modified Huber $\psi$-function with tuning constants $a$ and $b$.  [Duchesne](#literature) (1999) suggested the default parametrization $a=9$ and $b=0.25$.

```R


svymean_reg(object, auxiliary, theta, k, check.names, na.rm)



```



## 3 Variance

By a first-order Taylor series expansion, the GREG estimator of the $y$-total, $\widehat{t}_y^{\,greg}$, can be approximated by ([Särndal et al.](#literature), 1992, Result 6.6.1)
$$
\begin{align}
	\widehat{t}_y^{\,greg} & \approx t_y^0 = \sum_{i \in s}w_i y_i + (\bold t_x - \widehat{\bold t}_x)^T \boldsymbol \theta_N \nonumber \\
	& = \bold t_x^T\boldsymbol \theta_N +\sum_{i \in s}w_iE_i,\label{eq:greg_taylor}
\end{align}
$$
where $\boldsymbol \theta_N$ is the census fit parameter and $E_i$ is the census residual, $E_i = y_i - \bold x_i^T \boldsymbol \theta_N$. The first term on the r.h.s. of ($\ref{eq:greg_taylor}$) is a constant. Since $\widehat{t}_y^{\,greg}$ is an approximately unbiased estimator for $t_y$, its variance can be approximated by 
$$
\begin{equation*}
\mathrm{var}\left(\widehat{t}_y^{\,greg}\right) \approx \mathrm{var}\left(t_y^0\right) = \mathrm{var}\left(\sum_{i \in s}w_i E_i\right)
\end{equation*}
$$
which can be estimated—upon replacing $E_i$ by the sample residual $e_i$— by
$$
\begin{equation*}
	\widehat{\mathrm{var}}\left(\widehat{t}_y^{\,greg}\right) = \sum_{i \in s}\sum_{j \in s} \frac{\Delta_{ij}}{\pi_i}\frac{e_ie_j}{\pi_i\pi_j}.
\end{equation*}
$$

---

The GREG estimator of the $y$-total, $\widehat{t}_y^{\,greg}$, can be written using the $g$-weights as ([Särndal et al.](#literature), 1992, p. 232–233) $\widehat{t}_y^{\,greg} = \sum_{i \in s} g_i w_i y_i$, where
$$
\begin{equation*}
g_i = 1 + (\bold t_x - \widehat{\bold t}_x)^T \left(\sum_{i \in s} \frac{w_i}{v_i} \bold x_i \bold x_i^T\right)^{-1}\frac{\bold x_i}{v_i}.
\end{equation*}
$$
The census residuals are $E_i = y_i - \bold x_i^T \boldsymbol \theta_N$ for all $i \in U$. 
$$
\begin{align}
	\widehat{t}_y^{\,greg} &= \sum_{i \in s} g_i w_i y_i = \sum_{i \in s} g_i w_i(\bold x_i^T \boldsymbol \theta_N + E_i) \nonumber \\
	&= \sum_{i \in U} \bold x_i^T \boldsymbol \theta_N + \sum_{i \in s} g_iw_i E_i, \label{eq:greg_approx}
\end{align}
$$
where we have used the calibration property of the $g$-weights in ($\ref{eq:greg_approx}$), $\sum_{i \in s} g_i w_i \bold x_i = \sum_{i \in U} \bold x_i$. The first term on the r.h.s of ($\ref{eq:greg_approx}$) is a population quantity and does therefore not contribute to the variance of $\widehat{t}_y^{\;greg}$. Thus, we have

$$
\begin{equation}\label{eq:var_greg2}
	\mathrm{var}\left(\widehat{t}_y^{\,greg}\right) = \mathrm{var}\left(\sum_{i \in s} g_i w_i E_i\right),
\end{equation}
$$
which resembles the variance of a weighted total except that it also depends on $g_i$. Disregarding the fact that the $g$-weights are sample dependent and substituting the sample residual, $e_i$, for $E_i$ in ($\ref{eq:var_greg2}$), the variance of the GREG estimator of the $y$-total can be computed by the $g$-weighted variance ([Särndal et al.](#literature), 1992, Result 6.6.1)
$$
\begin{equation*}
	\widehat{\mathrm{var}}\left(\widehat{t}_y^{\;greg}\right) = \sum_{i \in s}\sum_{j \in s} \frac{\Delta_{ij}}{\pi_i}\frac{g_ie_ig_je_j}{\pi_i\pi_j}.
\end{equation*}
$$

## Variance

$$
\begin{equation*}
	\widehat{t}_y^{\,rob} = \bold t_x^T \boldsymbol \theta_{N} + \sum_{i \in s} b_i (y_i - \bold x_i^T \boldsymbol \theta_{N}) + (\bold t_x - \widehat{\bold t}_{bx})^T \big(\widehat{\boldsymbol \theta}_n - \boldsymbol \theta_N\big)
\end{equation*}
$$

where $\widehat{\bold t}_{bx} = \sum_{i \in s} b_i\bold x_i$.



## Variance QR

Let $(q_i, r_i)$ denote constants. The QR-predictor of the $y$-total is ([Wright](#literature), 1983)
$$
\begin{equation*}
    \widehat{t}_y^{\,qr} = \widehat{t}_{ry} + \big(\bold t_x - \widehat{\bold t}_{rx} \big)^T \widehat{\boldsymbol \theta}_{n}^{\,qr} = \widehat{t}_{ry} + \big(\bold t_x - \widehat{\bold t}_{rx} \big)^T \widehat{\bold T}_{qxx}^{-1} \widehat{\bold t}_{qxy},
\end{equation*}
$$
where $\widehat{t}_{ry} = \sum_{i \in s} r_i y_i$, $\widehat{\bold t}_{rx} = \sum_{i \in s}r_i \bold x_i$, $\widehat{\bold T}_{qxx} = \sum_{i \in s}q_i\bold x_i \bold x_i^T$, and $\widehat{\bold t}_{qxy} = \sum_{i \in s}q_i \bold x_i y_i$. The GREG estimator of the $y$-total obtains for the choice $(q_i = w_i/v_i, r_i = w_i)$. A QR-predictor is ADU if $r_i = w_i$. It can be written as $\widehat{t}_y^{\,qr} = \sum_{i \in s}w_i g_i^{qr}y_i$, where ([Duchesne](#literature), 1999)
$$
\begin{equation*}
	w_ig_i^{\,qr} = r_i + \big(\bold t_x - \widehat{\bold t}_{rx}\big)^T\left(\sum_{i \in s} q_i \bold x_i \bold x_i^T\right)^{-1} q_i \bold x_i.
\end{equation*}
$$

Let $\boldsymbol \theta_N$ be the census GM-estimate of $\boldsymbol \theta$ under model $\xi$, and denote its residual by $E_i = y_i - \bold x_i^T \boldsymbol \theta_N$. The GREG estimator of the $y$-total can be written as
$$
\begin{equation}\label{eq:qr_gweight}
	\widehat{t}_y^{\,qr} = \sum_{i \in s} w_i g_i^{\,qr} y_i = \sum_{i \in s} w_i g_i^{\,qr} (\bold x_i^T \boldsymbol \theta_N + E_i) = \sum_{i \in U} \bold x_i^T \boldsymbol \theta_N + \sum_{i \in s} w_i g_i^{\,qr} E_i,
\end{equation}
$$
where we have used the calibration property $\sum_{i \in s} w_i g_i^{\,qr} \bold x_i = \sum_{i \in U} \bold x_i$. Using ($\ref{eq:qr_gweight}$), we have
$$
\begin{equation*}
	\mathrm{var}\big(\widehat{t}_y^{\,qr}\big) = \mathrm{var}\left(\sum_{i \in s}w_ig_i^{\,qr}E_i\right)
\end{equation*}
$$
and can be estimated by
$$
\begin{equation*}
	\widehat{\mathrm{var}}\big(\widehat{t}_y^{\,qr}\big) = \mathrm{var}\left(\sum_{i \in s}w_ig_i^{\,qr}e_i\right),
\end{equation*}
$$
where $e_i$ is the sample residual.



---



The population fit is




Put
$$
\begin{equation*}
	\widehat{t}_y^{\,rob} = f\big(\widehat{t}_{by}, \widehat{\bold t}_{bx}, \widehat{\bold T}_{cxx}, \widehat{\bold t}_{cxy}\big)
\end{equation*}
$$
with partial derivatives evaluated at
$$
\begin{equation*}
	\tau = \big(t_y, \bold t_x, \bold T_{xx}, \bold t_{xy}\big)
\end{equation*}
$$

$$
\begin{equation*}
	\frac{\partial}{\partial \,\widehat{t}_{by}} f \big\vert_{\tau} = 1, \qquad \frac{\partial}{\partial \, \widehat{\bold t}_{bx}} f \big\vert_{\tau} = - \widehat{\boldsymbol \theta}_n \big\vert_{\tau} = -\widehat{\bold T}_{cxx}^{-1} \widehat{\bold t}_{cxy}\big\vert_{\tau}= -\boldsymbol \theta_N 
\end{equation*}
$$

$$
\begin{equation*}
	\frac{\partial}{\partial \,\widehat{\bold T}_{bxx}} f \big \vert_{\tau} = -\big(\bold t_x - \widehat{\bold t}_{bx} \big)^T \widehat{\bold T}_{cxx}^{-1}\widehat{\boldsymbol \theta}_n \big\vert_{\tau} = \bold 0, \qquad \frac{\partial}{\partial \,\widehat{\bold t}_{bxy}} f \big\vert_{\tau} = \big(\bold t_x - \widehat{\bold t}_{bx} \big)^T \widehat{\bold T}_{cxx}^{-1} \Big\vert_{\tau} = \bold 0
\end{equation*}
$$

So the first-order Taylor series expansion about $\tau$ is
$$
\begin{equation*}
	t_y^0 = f(\tau) + (\widehat{t}_{by} - t_y) - (\widehat{\bold t}_{bx} - \bold t_x)^T \boldsymbol \theta_N = \widehat{t}_{by} + \big(\bold t_x - \widehat{\bold t}_{bx}\big)^T \boldsymbol \theta_N
\end{equation*}
$$

## 9 New

### 9.1 Only representative outliers

* reference: $t_y$
* census parameter is estimated by weighted LS, $\boldsymbol \theta_{LS}$ (weighting is required because of heteroscedasticity)
* robust estimation with sample data incurs a bias => mse not variance

### 9.2 Both

* reference: $\sum_{i \in U} u_i y_i$, where $u_i$ is a population-level robustness weight
* census parameter is estimated by robust methods, $\boldsymbol \theta_{rob}$

* symmetric contamination distribution: robust estimation does not incur bias (if $E_i$  also have a symmetric distribution)
* otherwise (more realistic): bias

## GREG: first principles

Consider model
$$
\begin{equation*}
	\xi: \quad Y_i = \bold x_i^T \boldsymbol{\theta} + \sigma \sqrt{v_i}E_i, \qquad \boldsymbol{\theta} \in \R^p, \quad \sigma > 0, \quad i \in U,
\end{equation*}
$$
and note that $\widehat{t}_y=\sum_{i \in s} w_i y_i$ is a design-unbiased estimator of $t_y$, but it is not model-unbiased under model $\xi$. The model bias of $\widehat{t}_y$ is
$$
\begin{equation*}
	\mathrm{E}_{\xi}(\widehat{t}_y - t_y) = \sum_{i \in s} w_i \bold x_i^T\boldsymbol \theta - \sum_{i \in U} \bold x_i^T \boldsymbol \theta,
\end{equation*}
$$
given that $\mathrm{E}_{\xi}(E_i)=0$. If $\boldsymbol \theta$ were known, we could define the "corrected" estimator
$$
\begin{equation*}
	\widehat{t}_y^* = \sum_{i \in s}w_i y_i - \sum_{i \in s} w_i \bold x_i^T\boldsymbol \theta + \sum_{i \in U} \bold x_i^T \boldsymbol \theta,
\end{equation*}
$$
which is design- and model-unbiased as an estimator of $t_y$.

* estimator $\widehat{\boldsymbol \theta}$ must be model-unbiased
* estimator $\widehat{\boldsymbol \theta}$ must be design-consistent for $\boldsymbol \theta_N$



## Literature

==must be modified==

BEAUMONT, J.-F. AND A. ALAVI (2004). Robust Generalized Regression Estimation, *Survey Methodology* **30**, 195–208.

BEAUMONT, J.-F. AND L.-P. RIVEST (2009). Dealing with outliers in survey data, in *Sample Surveys: Theory, Methods and Inference*, ed. by D. Pfeffermann and C. R. Rao, Amsterdam: Elsevier, vol. 29A of *Handbook of Statistics*, chap. 11, 247–280.

CHAMBERS, R. (1986). Outlier Robust Finite Population Estimation, *Journal of the American Statistical Association* **81**, 1063–1069.

DUCHESNE, P. (1999). Robust calibration estimators, *Survey Methodology* **25**, 43–56.

HAMPEL, F. R., E. M. RONCHETTI, P. J. ROUSSEEUW, AND W. A. STAHEL (1986). *Robust Statistics: The Approach Based on Influence Functions*, New York: John Wiley & Sons.

LEE, H. (1995a). Outliers in business surveys, in *Business survey methods*, ed. by B. G. Cox, D. A. Binder, B. N. Chinnappa, A. Christianson, M. J. Colledge, and P. S. Kott, New York: John Wiley and Sons, chap. 26, 503–526.

SÄRNDAL, C.-E., B. SWENSSON, AND J. WRETMAN (1992). Model Assisted Survey Sampling, New York: Springer.

WRIGHT, R. L. (1983). Finite population sampling with multivariate auxiliary information, *Journal of the American Statistical Association* **78**, 879–884.



