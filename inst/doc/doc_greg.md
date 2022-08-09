

# Robust Generalized Regression Predictor

Tobias Schoch

[toc]

## 1 Introduction

The *population* regression model is given by
$$
\begin{equation}\label{eq:reg_model}
	\xi: \quad Y_i = \bold x_i^T \boldsymbol{\theta} + \sigma \sqrt{v_i}E_i, \qquad \boldsymbol{\theta} \in \R^p, \quad \sigma > 0, \quad i \in U,
\end{equation}
$$
where

* the population $U$ is of size $N$,

* the parameters $\boldsymbol{\theta}$ and $\sigma$ are unknown,
* the $\bold x_i$'s are known values (possibly containing outliers), $\bold x_i \in \R^p$, $1 \leq p < N$; we denote the design matrix by $\bold X = (\bold x_1, \ldots, \bold x_n)^T$,
* the $v_i$'s are known positive (heteroscedasticity) constants,
* the errors $E_i$ are independent and identically distributed (i.i.d.) random variables (r.v.) with zero expectation and unit variance,
* it is assumed that $\sum_{i \in U} \bold x_i \bold x_i^T / v_i$ is a non-singular $(p \times p)$ matrix.

*Remarks.* The i.i.d. assumption on the errors $E_i$ is rather strict. This assumption can be replaced by the assumption that the $E_i$ are identically distributed
r.v. such that
$\mathbb{E}_{\xi} (E_i \mid \bold x_i, \ldots, \bold x_N) = 0$ and
$\mathbb{E}_{\xi} (E_i E_j \mid \bold x_i, \ldots, \bold x_N) = 1$ if $i = j$
and zero otherwise for all $i,j \in U$, where $\mathbb{E}_{\xi}$ denotes expectation w.r.t. model $\xi$ in ($\ref{eq:reg_model}$). Another generalization obtains by requiring that $\mathbb{E}_{\xi}(E_i\bold x_i)
    = \mathbb{E}_{\xi}(\bold x_i E_i) = \bold 0$ in place of the conditional
expectation. If the distribution of the errors $E_i$ is asymmetric with non-zero mean,
the regression intercept and the errors are confounded. The slope parameters,
however, are identifiable with asymmetric distributions
[(Carroll and Welsh, 1988)](#literature). In the context of GREG prediction, however, we deal with prediction under
the model. Thus, identifiability is not an issue.

It is assumed that a sample $s$ is drawn from $U$ with sampling design $p(s)$ such that the independence (orthogonality) structure of the model errors in ($\ref{eq:reg_model}$) is maintained. The sample regression *M*- and *GM*-estimator of $\boldsymbol{\theta}$ are defined as the root to the following estimating equations (cf. [Hampel et al., 1986](#literature), Chapter 6.3)
$$
\begin{align*}
	\sum_{i \in s} \frac{w_i}{\sqrt{v_i}} \psi_k (r_i)\bold x_i &= \bold 0 & \text{($M$-estimator)},\\
	\sum_{i \in s} \frac{w_i}{\sqrt{v_i}} h(\bold x_i) \psi_k (r_i) \bold x_i &= \bold 0 & \text{(Mallows $GM$-estimator)},\\
	\sum_{i \in s} \frac{w_i}{\sqrt{v_i}} h(\bold x_i) \psi_k\left(\frac{r_i}{h(\bold x_i)}\right) \bold x_i &= \bold 0& \text{(Schweppe $GM$-estimator)},
\end{align*}
$$
where

* $w_i$ is the sampling weight,
* $\psi_k$ is a *generic* $\psi$-function indexed by the robustness tuning constant $k$,

* $r_i$ is the standardized residual, defined as

$$
\begin{equation}\label{eq:residuals}
	r_i = \frac{y_i - \bold x_i^T\boldsymbol \theta}{\sigma\sqrt{v_i}},
\end{equation}
$$

* $h:\R^p \rightarrow \R_+ $ is a weight function,
* $\sigma$ is the regression scale which is estimated by the (normalized) weighted median of the absolute deviations from the weighted median of the residuals.

The Huber and Tukey bisquare (biweight) $\psi$-functions are denoted by, respectively, $\psi_k^{hub}$ and $\psi_k^{tuk}$. The sample-based estimators of $\boldsymbol \theta$ can be written as a weighted least squares problem
$$
\begin{equation*}
\sum_{i \in s} \frac{w_i}{v_i} u_i(r_i,k) (y_i - \bold x_i^T \widehat{\boldsymbol \theta}_n) \bold x_i = \bold 0,
\end{equation*}
$$
where 
$$
\begin{equation}\label{eq:ui}
    u_i(r_i,k) = \begin{cases}
        \displaystyle{\frac{\psi_k(r_i)}{r_i}} & \qquad \text{$M$-estimator},\\
        h(\bold x_i)\displaystyle{\frac{\psi_k(r_i)}{r_i}} &
            \qquad \text{Mallows $GM$-estimator},\\
        \displaystyle{\frac{\psi_k(r_i^*)}{r_i^*}}, \quad \text{where} \quad
            r_i^* = \frac{r_i}{h(\bold x_i)} & \qquad \text{Schweppe $GM$-estimator},
    \end{cases}
\end{equation}
$$

and $k$ denotes the robustness tuning constant.



## 2 Representation of the robust GREG as a QR-predictor

The robust GREG predictor of the population $y$-total can be written in terms of the $g$-weights [(see e.g., Särndal et al., 1992, Chapter 6)](#literature) as
$$
\begin{equation}\label{eq:gweighted_total}
	\widehat{t}_{y}^{\,rob} = \sum_{i \in s}g_i y_i,
\end{equation}
$$
where the $g$-weights are defined as [(Duchesne, 1999)](#literature)
$$
\begin{equation}\label{eq:gweights}
	g_i = b_i + \big(\bold t_{\bold x} - \widehat{\bold t}_{b\bold x}\big)^T \left(\sum_{i \in s}q_i \bold x_i \bold x_i^T\right)^{-1}q_i\bold x_i,
\end{equation}
$$
where $\widehat{\bold t}_{b\bold x} = \sum_{i \in s}b_i \bold x_i$ and $\bold t_{\bold x} = \sum_{i \in U}\bold x_i$. The sampling weights, $w_i$, are "embedded" into the $g$-weights in ($\ref{eq:gweighted_total}$).

In contrast to the non-robust "standard" GREG predictor, the $g$-weights in ($\ref{eq:gweights}$) depend on the study variable, $y_i$, through the choice of the constants $(q_i, b_i) = \{(q_i,b_i): i \in s\}$. This will be easily recognized once we define the set of constants. The predictors of the population $y$-total that are defined in terms of the constants $(q_i,r_i)$ form the class of *QR-predictor* due to [Wright (1983)](#literature). 

<div style="padding-left: 1.5rem;
    padding-top: 0.25rem;
    padding-bottom: 0.25rem;
    margin-top: 0.25rem;
    margin-bottom: 0.25rem;
    border: 0px solid #eee;
    border-left-width: 0.75rem;
    border-radius: .25rem;
    border-left-color: #ce5b00;">
   <b style="color: #ce5b00;">Important.</b> We denote the constants by <i>(q<sub>i</sub>, b<sub>i</sub>)</i> instead of <i>(q<sub>i</sub>, r<sub>i</sub>)</i> because <i>r<sub>i</sub></i> is our notation for the residuals.
</div>

In passing we note that $\widehat{t}_y^{\,rob}$ can be expressed in a "standard" GREG representation. Let
$$
\begin{equation*}
	\widehat{\boldsymbol \theta} = \left(\sum_{i \in s}q_i \bold x_i \bold x_i^T\right)^{-1} \sum_{i \in s} q_i\bold x_iy_i,
\end{equation*}
$$
then $\widehat{t}_y^{\,rob}$ in ($\ref{eq:gweighted_total}$) can be written as
$$
\begin{equation*}
	\widehat{t}_{y}^{\,rob} = \sum_{i \in s} b_i y_i + \big(\bold t_{\bold x} - \widehat{\bold t}_{b\bold x}\big)^T \widehat{\boldsymbol \theta}
	= \bold t_{\bold x}^T\widehat{\boldsymbol \theta} + \sum_{i \in s}b_i (y_i - \bold x_i^T\widehat{\boldsymbol \theta}).
\end{equation*}
$$

In the next two sections, we define the constants $(q_i, b_i)$ of the QR-predictor.

### 2.1 Constants $q_i$ of the QR-predictor 

The set of constants $\{q_i\}$ is defined as 
$$
\begin{equation}\label{eq:qi}
    q_i =\frac{w_i \cdot u_i(r_i,k)}{v_i}, \qquad i = 1,\ldots,n,
\end{equation}
$$
where $v_i$ is given in ($\ref{eq:reg_model}$) and $u_i(r_i,k)$ is defined in ($\ref{eq:ui}$). The tuning constant $k$ in $u_i(r_i,k)$ is the one that is used to estimate $\boldsymbol \theta$.

### 2.2 Constants $b_i$ of the QR-predictor

The constants $\{b_i\}$ are predictor-specific. They depend on the argument `type`. Moreover, the $b_i$'s depend on the robustness tuning constant `k` —which is an argument of `svymean_reg()` and `svytotal_reg()`—to control the robustness of the prediction. To distinguish it from the tuning constant $k$, which is used in fitting model $\xi$ in ($\ref{eq:reg_model}$), it will be denoted by $\kappa$. Seven sets $\{b_i\}$ are available.

* `type = "projective"`: $b_i \equiv 0$ [(Särndal and Wright, 1984)](#literature),

* `type = "ADU"`: $b_i \equiv w_i$ [(Särndal et al., 1992, Chapter 6)](#literature),

  `type = "huber"`: $b_i \equiv w_i \cdot u_i(r_i,\kappa)$, where $u_i$ is defined in ($\ref{eq:ui}$) with $\psi_k \equiv \psi_k^{hub}$ [(Lee, 1995; Hulliger, 1995; Beaumont and Alavi, 2004)](#literature),

* `type = "tukey"`: $b_i \equiv w_i \cdot u_i(r_i, \kappa)$, where $u_i$ is defined in ($\ref{eq:ui}$) with $\psi_k \equiv \psi_k^{tuk}$[(Lee, 1995; Hulliger, 1995; Beaumont and Alavi, 2004)](#literature)  , 

* `type = "lee"`: $b_i \equiv \kappa \cdot w_i$, where $0 \leq \kappa \leq 1$ [(Lee, 1991, 1995)](#literature), 

* `type = "BR"`: $b_i \equiv w_i \cdot u_i(r_i,\kappa)$, where $u_i$ is defined in ($\ref{eq:ui}$) with $\psi_k$ replaced by [(Beaumont and Rivest, 2009)](#literature)
  $$
  \begin{equation*}
  	\psi_k^{mod}(x) = \frac{x}{w_i} + \frac{w_i - 1}{w_i}\psi_k^{hub}(x),
  \end{equation*}
  $$

* `type = "duchesne"`: $b_i \equiv w_i \cdot u_i(r_i; a,b)$, where $u_i$ is defined in ($\ref{eq:ui}$) with $\psi_k$ replaced by [(Duchesne, 1999)](#literature)
  $$
  \begin{equation*}
  	\psi_{a,b}^{hub}(x) = \begin{cases} x & \text{if} \quad \vert x \vert \leq a, \\ a \cdot \mathrm{sign}(x) & \text{if} \quad \vert x \vert > a \quad \text{and} \quad \vert x \vert < a/b, \\ b \cdot x & \text{if} \quad \vert x \vert > a/b, 
          \end{cases}
  \end{equation*}
  $$
  where $\psi_{a,b}^{hub}$ is a modified Huber $\psi$-function with tuning constants $a$ and $b$ (in place of $\kappa$). [Duchesne](#literature) (1999) suggested the default parametrization $a=9$ and $b=0.25$.

### 2.3 Implementation

Let $\bold q = (q_1, \ldots, q_n)^T$ and $\bold b = (b_1, \ldots, b_n)^T$, where $q_i$ and $b_i$ are defined in, respectively, ($\ref{eq:qi}$) and Section 2.2. Put $\bold Z = \sqrt{\bold q} \circ \bold X$, where $\circ$ denotes Hadamard multiplication and the square root is applied element by element. The vector-valued $g$-weights, $\bold g = (g_1, \ldots, g_n)^T$, in ($\ref{eq:gweights}$) can be written as
$$
\begin{equation*}
	\bold g^T = \bold b^T + \big(\bold t_{\bold x} - \widehat{\bold t}_{b\bold x}\big)^T \underbrace{\big(\bold Z^T \bold Z\big)^{-1} \bold Z^T}_{=\bold H, \; \text{say}} \circ (\sqrt{\bold q})^T.
\end{equation*}
$$
Define the QR factorization $\bold Z = \bold Q \bold R$, where $\bold Q$ is an orthogonal matrix and $\bold R$ is an upper triangular matrix (both of conformable size). Note that the matrix QR-factorization and Wright's QR-estimators have nothing in common besides the name; in particular, $\bold q$ and $\bold Q$ are unrelated. With this we have
$$
\begin{equation*}
	\bold H = \big(\bold Z^T\bold Z\big)^{-1} \bold Z^T =
   \bold R^{-1}\bold Q^T
\end{equation*}
$$
and multiplying both sides by $\bold R$, we get $\bold R\bold H = \bold Q^T$ which can be solved easily for $\bold H$ since $\bold R$ is an upper triangular matrix (see `base::backsolve()`). Thus, the $g$-weights can be computed as
$$
\begin{equation*}
	\bold g = \bold b + \bold H^T \big(\bold t_{\bold x} - \widehat{\bold t}_{b\bold x}\big) \circ \sqrt{\bold q},
\end{equation*}
$$

where the $(p \times n)$ matrix $\bold H$ need not be explicitly transposed when using `base::crossprod()`. The terms $\bold b$ and $\widehat{\bold t}_{b \bold x}$ are easy to compute. Thus,
$$
\begin{equation*}
	\widehat{t}_y^{\;rob} = \bold g^T\bold y, \qquad \text{where} \quad \bold y = (y_1, \ldots, y_n)^T.
\end{equation*}
$$

## 3 Variance estimation

<div style="padding-left: 1.5rem;
    padding-top: 0.25rem;
    padding-bottom: 0.25rem;
    margin-top: 0.25rem;
    margin-bottom: 0.25rem;
    border: 0px solid #eee;
    border-left-width: 0.75rem;
    border-radius: .25rem;
    border-left-color: #ce5b00;">
   <b style="color: #ce5b00;">Important.</b> Inference of the regression estimator is only implemented under the assumption of representative outliers (in the sense of Chambers, 1986). We do not cover inference in presence of nonrepresentative
outliers.
</div>

Our discussion for variance estimation follows the line of reasoning in [Särndal et al., (1992](#literature), 233-234) on the variance of the non-robust GREG estimator. To this end, denote by $E_i = y_i - \bold x_i^T \boldsymbol \theta_N$, $i \in U$, the census residuals, where $\boldsymbol \theta_N$ is the census parameter. With this, any $g$-weighted predictor can be written as
$$
\begin{equation}\label{eq:census_fit}
	\widehat{t}_y^{\,rob} = \sum_{i \in s} g_i y_i = \sum_{i \in s} g_i (\bold x_i^T \boldsymbol \theta_N + E_i)
	= \sum_{i \in U} \bold x_i^T \boldsymbol \theta_N + \sum_{i \in s} g_i E_i, 
\end{equation}
$$
where we have used the fact that the $g$-weights in ($\ref{eq:gweights}$) satisfy the calibration property
$$
\begin{equation*}
	\sum_{i \in s} g_i \bold x_i = \sum_{i \in U} \bold x_i.
\end{equation*}
$$
The first term on the r.h.s. of the last equality in ($\ref{eq:census_fit}$) is a population quantity and does therefore not contribute to the variance of $\widehat{t}_y^{\,rob}$. Thus, we can calculate the variance of the robust GREG predictor by
$$
\begin{equation}\label{eq:var_censusfit}
	\mathrm{var}\left(\widehat{t}_y^{\,rob}\right) = \mathrm{var}\left(\sum_{i \in s} g_i E_i\right)
\end{equation}
$$
under the assumptions that (1) the $E_i$ are known quantities and (2) the $g_i$ do not depend on the $y_i$. 

Disregarding the fact that the $g$-weights are sample dependent and substituting the sample residual $r_i$ for $E_i$ in ($\ref{eq:var_censusfit}$), [Särndal et al. (1992](#literature), 233-234 and Result 6.6.1) propose to estimate the variance of the GREG predictor by the $g$-weighted variance of the total $\sum_{i \in s}g_ir_i$. Following the same train of thought and disregarding in addition that the $g_i$ depend on $y_i$, the variance of $\widehat{t}_y^{\,rob}$ can be approximated by 
$$
\begin{equation*}
	\widehat{\mathrm{var}}\big(\widehat{t}_y^{\,rob}\big) \approx \widehat{\mathrm{var}}\left(\sum_{i \in s}g_i r_i\right),
\end{equation*}
$$
where $\widehat{\mathrm{var}}(\cdot)$ denotes a variance estimator of a total for the sampling design $p(s)$.

## Literature

BEAUMONT, J.-F. AND A. ALAVI (2004). Robust Generalized Regression Estimation, *Survey Methodology* **30**, 195–208.

BEAUMONT, J.-F. AND L.-P. RIVEST (2009). Dealing with outliers in survey data, in *Sample Surveys: Theory, Methods and Inference*, ed. by D. Pfeffermann and C. R. Rao, Amsterdam: Elsevier, vol. 29A of *Handbook of Statistics*, chap. 11, 247–280. [DOI:10.1016/ S0169-7161(08)00011-4](https://doi.org/10.1016/S0169-7161(08)00011-4)

CARROLL, R. J. AND A. H. WELSH (1988). A Note on Asymmetry and Robustness in Linear Regression.
*The American Statistician* **42**, 285–287. [DOI:10.1080/00031305.1988.10475591](https://doi.org/10.1080/00031305.1988.10475591)

CHAMBERS, R. (1986). Outlier Robust Finite Population Estimation. *Journal of the American Statistical Association* **81**, 1063–1069, [DOI:10.1080/01621459.1986.10478374](https://doi.org/10.1080/01621459.1986.10478374).

DUCHESNE, P. (1999). Robust calibration estimators, *Survey Methodology* **25**, 43–56.

HAMPEL, F. R., E. M. RONCHETTI, P. J. ROUSSEEUW, AND W. A. STAHEL (1986). *Robust Statistics: The Approach Based on Influence Functions*, New York: John Wiley & Sons. [DOI:10.1002/9781118186435](https://doi.org/10.1002/9781118186435)

HULLIGER, B. (1995). Outlier Robust Horvitz–Thompson Estimators. *Survey Methodology* **21**, 79–87.

LEE, H. (1995). Outliers in business surveys, in *Business survey methods*, ed. by B. G. Cox, D. A. Binder, B. N. Chinnappa, A. Christianson, M. J. Colledge, and P. S. Kott, New York: John Wiley and Sons, chap. 26, 503–526. [DOI:10.1002/9781118150504.ch26](https://doi.org/10.1002/9781118150504.ch26)

SÄRNDAL, C.-E., B. SWENSSON, AND J. WRETMAN (1992). *Model Assisted Survey Sampling*, New York: Springer. 

SÄRNDAL, C.-E. AND R. L. WRIGHT (1984). Cosmetic Form of Estimators in Survey Sampling. *Scandinavian Journal of Statistics* **11**, 146–156.

WRIGHT, R. L. (1983). Finite population sampling with multivariate auxiliary information, *Journal of the American Statistical Association* **78**, 879–884. [DOI:10.1080/01621459.1983.10477035](https://doi.org/10.1080/01621459.1983.10477035).

