# Chapter 20 — Statistical Analysis

[← Previous](19-probability.md) · [Table of Contents](../00-front/table-of-contents.md) · [Next →](21-random-processes.md)

> **Prerequisites:** Ch. 5 (Integral Calculus), Ch. 6 (Multivariable Calculus), Ch. 8 (Matrices and Linear Systems), Ch. 19 (Probability)
> **Learning outcomes:** After this chapter you can …
> - summarise a data set with location, spread and shape statistics, and choose robust statistics when outliers are present;
> - compute covariance, correlation coefficients and covariance matrices, and interpret them correctly;
> - recognise the standard distributions (Bernoulli, binomial, Poisson, uniform, exponential, Gaussian, Rayleigh, chi-square, $t$, $F$), their moments and the relations between them;
> - derive sampling distributions, apply the central limit theorem, and build estimators by the method of moments and maximum likelihood;
> - construct confidence intervals and carry out $z$, $t$, chi-square and $F$ tests with correct interpretation of $p$-values and error probabilities;
> - fit and diagnose a linear regression, and evaluate and report measurement uncertainty in the GUM framework.

Probability theory (Ch. 19) starts from a model — a distribution with known parameters — and deduces what data should look like. **Statistics** runs the argument backwards: it starts from data and infers something about the model that produced them. An engineer who measures the noise floor of a receiver, calibrates a temperature sensor, qualifies a batch of capacitors or reports "$4.702\ \mathrm{k\Omega} \pm 0.012\ \mathrm{k\Omega}$" is doing statistics, whether or not the word is used.

The chapter follows the natural order of a data analysis. We first describe data (Sections 20.1–20.3), then recall the distributions used as models (20.4). Sections 20.5–20.8 form the core of *inference*: how statistics computed from a sample are themselves random (sampling distributions), how to estimate parameters, how to attach intervals, and how to test claims. Section 20.9 fits relations between variables by regression, and Section 20.10 assembles everything into the internationally agreed method for evaluating measurement uncertainty (the *GUM*).

Throughout, $X_1,\dots,X_n$ denotes a random sample (random variables, uppercase) and $x_1,\dots,x_n$ the observed numbers (lowercase), following the notation of Ch. 19.

## 20.1 Descriptive statistics

### Why describe before inferring

Before any model is fitted, a data set should be *looked at*. A single miswired channel, a saturated ADC, a unit error between millivolts and volts, or a bimodal population (two suppliers mixed in one reel) will invalidate every later calculation, and each is obvious from a summary or a plot but invisible in a single average. Descriptive statistics compress $n$ numbers into a handful of meaningful ones and are the first defence against garbage-in, garbage-out.

> [!definition] Data, population and sample
> A **population** is the complete set of items or outcomes about which we want to draw conclusions (all resistors produced by a process, all possible readings of a voltmeter on a fixed source). A **sample** is the finite subset $x_1,\dots,x_n$ actually observed; $n$ is the **sample size**. A **statistic** is any function of the sample, $t = g(x_1,\dots,x_n)$, that does not depend on unknown parameters.

Data come in several types, and the type dictates the admissible summaries:

| Type | Example | Meaningful operations |
|---|---|---|
| Nominal (categorical) | failure mode: open / short / drift | counting, mode |
| Ordinal | inspection grade A/B/C | ordering, median, ranks |
| Interval | temperature in °C | differences, mean, standard deviation |
| Ratio | resistance in Ω, time to failure | all of the above plus ratios, coefficient of variation |

Averaging category codes ("open = 1, short = 2") is meaningless; quoting a coefficient of variation for a temperature in °C is meaningless because the zero is arbitrary.

### Order statistics and the five-number summary

Sorting the sample gives the **order statistics** $x_{(1)}\le x_{(2)}\le\cdots\le x_{(n)}$. Many robust summaries are built from them.

> [!definition] Quantiles, quartiles and the five-number summary
> For $0<p<1$, a **$p$-quantile** $q_p$ is a value with (approximately) a fraction $p$ of the data below it and $1-p$ above. The **median** is $q_{0.5}$; the **lower and upper quartiles** are $Q_1=q_{0.25}$ and $Q_3=q_{0.75}$. The **interquartile range** is $\mathrm{IQR}=Q_3-Q_1$. The **five-number summary** is $\bigl(x_{(1)},\,Q_1,\,\tilde x,\,Q_3,\,x_{(n)}\bigr)$, where $\tilde x$ is the median.

Several conventions exist for quartiles of a finite sample (software packages offer up to nine); for hand work we use the *median of halves*: $Q_1$ is the median of the lower half of the sorted data and $Q_3$ the median of the upper half (excluding the overall median when $n$ is odd). Differences between conventions shrink as $n$ grows and do not matter for the decisions quartiles are used for.

The **box plot** draws a box from $Q_1$ to $Q_3$ with a line at the median and whiskers to the most extreme points within the **fences**

$$
Q_1 - 1.5\,\mathrm{IQR}, \qquad Q_3 + 1.5\,\mathrm{IQR}.
$$

Points beyond the fences are flagged as potential outliers. For Gaussian data, about $0.7\%$ of points fall outside the fences, so a flagged point is a reason to *investigate*, not automatically to delete.

### Frequency tables and histograms

For larger $n$ the data are grouped into $K$ class intervals (bins) of width $h$, with **frequencies** $f_k$ (counts), **relative frequencies** $f_k/n$ and **frequency density** $f_k/(nh)$. A histogram plotted with frequency density has unit total area and is directly comparable with a probability density function $f_X(x)$ from Ch. 19 — this is the bridge between data and model.

The choice of bin width matters: too wide hides structure, too narrow shows noise. Two common rules are

$$
K \approx 1+\log_2 n \quad(\text{Sturges}), \qquad h \approx 2\,\mathrm{IQR}\,n^{-1/3}\quad(\text{Freedman–Diaconis}).
$$

Sturges' rule assumes near-Gaussian data and undersmooths large samples; Freedman–Diaconis uses the IQR and is robust to outliers. The **empirical cumulative distribution function**

$$
\hat F_n(x) = \frac{1}{n}\,\#\{i : x_i \le x\}
$$

needs no binning at all: it is a staircase rising by $1/n$ at each data point, and by the Glivenko–Cantelli theorem it converges uniformly to the true CDF $F_X(x)$ as $n\to\infty$.

### Shape: skewness and kurtosis

Location and spread (Section 20.2) are not the whole story. The *shape* of a distribution is summarised by standardised third and fourth moments. With $m_k = \frac1n\sum_i (x_i-\bar x)^k$ the $k$th sample central moment,

$$
g_1 = \frac{m_3}{m_2^{3/2}} \quad(\text{sample skewness}), \qquad g_2 = \frac{m_4}{m_2^{2}} - 3 \quad(\text{sample excess kurtosis}).
$$

Skewness is zero for symmetric data, positive when the right tail is long (time-to-failure, power of noise, lognormal component values) and negative when the left tail is long. Excess kurtosis is zero for a Gaussian, positive for heavy-tailed data with more extreme values (impulsive noise, interference), and negative for flat, light-tailed data (a uniform distribution has $g_2=-1.2$). Because they involve high powers of deviations, $g_1$ and $g_2$ are very sensitive to outliers and need $n$ of at least several tens to be informative.

> [!example] Worked Example 20.1 — Summarising a batch of resistors
> Twelve resistors of nominal value $1\ \mathrm{k\Omega}$ are measured (in Ω):
> $$
> 998,\ 1003,\ 1001,\ 995,\ 1007,\ 1002,\ 999,\ 1004,\ 1000,\ 1012,\ 997,\ 1001 .
> $$
> Find the five-number summary, the IQR, the fences, and comment on the shape.
>
> **Step 1 — sort.** $995,\ 997,\ 998,\ 999,\ 1000,\ 1001,\ 1001,\ 1002,\ 1003,\ 1004,\ 1007,\ 1012$.
>
> **Step 2 — median.** With $n=12$ (even) the median is the average of the 6th and 7th values: $\tilde x = (1001+1001)/2 = 1001\ \Omega$.
>
> **Step 3 — quartiles (median of halves).** Lower half $995,\dots,1001$: $Q_1 = (998+999)/2 = 998.5\ \Omega$. Upper half $1001,\dots,1012$: $Q_3=(1003+1004)/2 = 1003.5\ \Omega$. Hence $\mathrm{IQR}=5\ \Omega$.
>
> **Step 4 — five-number summary.** $(995,\ 998.5,\ 1001,\ 1003.5,\ 1012)\ \Omega$; range $=17\ \Omega$.
>
> **Step 5 — fences.** $Q_1-1.5\,\mathrm{IQR} = 991\ \Omega$, $Q_3+1.5\,\mathrm{IQR}=1011\ \Omega$. The value $1012\ \Omega$ lies beyond the upper fence and is flagged.
>
> **Step 6 — shape.** Working with deviations $d_i = x_i-1000$ keeps the arithmetic small: $\sum d_i = 19$, so $\bar x = 1001.583\ \Omega$. The central moments are $m_2=19.41\ \Omega^2$, $m_3=71.42\ \Omega^3$, $m_4=1267\ \Omega^4$, giving $g_1 = 71.42/19.41^{1.5}=0.84$ and $g_2 = 1267/19.41^2-3=0.36$. The positive skewness is driven almost entirely by the flagged point. Before reporting statistics, one should check whether the $1012\ \Omega$ reading is a genuine part (perhaps from a different reel) or a contact-resistance error in the fixture.

> [!example] Worked Example 20.2 — Statistics from a frequency table
> The rise times of 40 logic buffers are grouped into 1 ns classes:
>
> | Class (ns) | 4–5 | 5–6 | 6–7 | 7–8 | 8–9 |
> |---|---|---|---|---|---|
> | Frequency $f_k$ | 3 | 9 | 15 | 10 | 3 |
>
> Estimate the mean and standard deviation.
>
> **Step 1 — class marks.** Represent each class by its midpoint $c_k = 4.5,\,5.5,\,6.5,\,7.5,\,8.5$ ns.
>
> **Step 2 — mean.** $\bar x \approx \frac{1}{n}\sum_k f_k c_k = \frac{3(4.5)+9(5.5)+15(6.5)+10(7.5)+3(8.5)}{40} = \frac{261}{40}=6.525\ \mathrm{ns}$.
>
> **Step 3 — variance.** $s^2 \approx \frac{1}{n-1}\sum_k f_k (c_k-\bar x)^2 = \frac{42.975}{39} = 1.102\ \mathrm{ns^2}$, so $s\approx 1.05\ \mathrm{ns}$.
>
> **Step 4 — grouping error.** Replacing every value by its class mark adds an error that is roughly uniform over a class of width $h$; Sheppard's correction subtracts $h^2/12 = 0.083\ \mathrm{ns^2}$ from the variance, giving $s\approx1.01\ \mathrm{ns}$. The correction is small when $h \ll s$ and is often ignored, but it shows why grouped statistics are approximate.

> [!tip] Engineering Note — Look at the raw data first
> Production test logs routinely contain sentinel values ($-999$, $9.9\times10^{37}$ "overload" from a DMM, zero for "no reading"). A five-number summary exposes them instantly because they appear as $x_{(1)}$ or $x_{(n)}$. Make computing the minimum and maximum the reflex first step of every analysis.

> [!warning] Common Pitfall — Deleting outliers to improve the statistics
> Removing points only because they are inconvenient biases every later result (spread is underestimated, confidence intervals become too narrow). Delete a point only with a documented physical cause (fixture fault, logging error). Otherwise keep it and, if needed, use robust statistics (Section 20.2) whose value it cannot dominate.

## 20.2 Mean, median, variance, and standard deviation

### Measures of location

Three numbers compete to describe "the typical value".

> [!definition] Sample mean, median and mode
> - **Sample mean:** $\displaystyle \bar x = \frac1n\sum_{i=1}^n x_i$.
> - **Sample median:** $\tilde x = x_{((n+1)/2)}$ for $n$ odd, $\tilde x = \tfrac12\bigl(x_{(n/2)}+x_{(n/2+1)}\bigr)$ for $n$ even.
> - **Mode:** the most frequent value (or, for continuous grouped data, the class with the highest density).

Each is the solution of a minimisation problem, which explains its behaviour:

$$
\bar x = \arg\min_{c}\sum_{i=1}^n (x_i-c)^2, \qquad \tilde x = \arg\min_{c}\sum_{i=1}^n \lvert x_i-c\rvert .
$$

*Proof sketch.* Differentiating $\sum (x_i-c)^2$ gives $-2\sum(x_i-c)=0$, i.e. $c=\bar x$; the second derivative $2n>0$ confirms a minimum. For $\sum\lvert x_i-c\rvert$ the derivative (where it exists) is $\#\{x_i<c\}-\#\{x_i>c\}$, which changes sign exactly when half the points lie on each side — at the median. $\square$

Squared loss punishes large deviations heavily, so a single wild point drags the mean; absolute loss treats all deviations linearly, so the median moves at most to the next data point. This is formalised by the **breakdown point**, the fraction of data that can be made arbitrarily bad before the statistic becomes arbitrarily bad: $0$ (strictly $1/n$) for the mean, $50\%$ for the median.

The mean has compensating virtues: it is linear ($\overline{ax+b}=a\bar x+b$), easy to update, and for Gaussian data it is the most efficient estimator of the centre (Section 20.6) — the median needs about $\pi/2\approx1.57$ times as many samples to achieve the same variance.

A useful compromise is the **trimmed mean**, which discards a fixed fraction (say $10\%$) from each end before averaging. It is common in sensor fusion and in reading-averaging routines in instrument firmware.

### Weighted mean

When observations have different reliabilities, they should not count equally. If $x_i$ has standard uncertainty $u_i$, the **inverse-variance weighted mean** is

$$
\bar x_w = \frac{\sum_i w_i x_i}{\sum_i w_i}, \qquad w_i = \frac{1}{u_i^2}, \qquad u(\bar x_w) = \Bigl(\sum_i w_i\Bigr)^{-1/2}.
$$

It is the minimum-variance unbiased linear combination of independent measurements of the same quantity: minimising $\operatorname{Var}\bigl(\sum a_i x_i\bigr)=\sum a_i^2u_i^2$ subject to $\sum a_i=1$ with a Lagrange multiplier (Ch. 6) gives $a_i\propto 1/u_i^2$.

### Measures of spread

> [!definition] Sample variance and standard deviation
> $$
> s^2 = \frac{1}{n-1}\sum_{i=1}^n (x_i-\bar x)^2, \qquad s = \sqrt{s^2}.
> $$
> The **range** is $x_{(n)}-x_{(1)}$; the **mean absolute deviation** is $\frac1n\sum\lvert x_i-\bar x\rvert$; the **median absolute deviation** is $\mathrm{MAD}=\operatorname{median}_i\lvert x_i-\tilde x\rvert$. The **coefficient of variation** (ratio data only) is $\mathrm{CV}=s/\bar x$.

The standard deviation $s$ carries the units of the data and is the usual measure of spread. For Gaussian data, $1.4826\,\mathrm{MAD}$ and $\mathrm{IQR}/1.349$ are consistent estimates of $\sigma$ that are robust to outliers.

### Why divide by $n-1$?

The divisor $n-1$ (Bessel's correction) makes $s^2$ an *unbiased* estimator of the population variance $\sigma^2$. The derivation is short and worth knowing.

> [!theorem] Unbiasedness of the sample variance
> If $X_1,\dots,X_n$ are independent with common mean $\mu$ and variance $\sigma^2$, then $E[S^2]=\sigma^2$, where $S^2=\frac{1}{n-1}\sum (X_i-\bar X)^2$.

*Proof.* Write $X_i-\bar X = (X_i-\mu)-(\bar X-\mu)$ and expand:

$$
\sum_{i=1}^n (X_i-\bar X)^2 = \sum_{i=1}^n (X_i-\mu)^2 - n(\bar X-\mu)^2 ,
$$

because the cross term is $-2(\bar X-\mu)\sum(X_i-\mu) = -2n(\bar X-\mu)^2$. Taking expectations, $E[(X_i-\mu)^2]=\sigma^2$ and $E[(\bar X-\mu)^2]=\operatorname{Var}(\bar X)=\sigma^2/n$ (Ch. 19, variance of a sum of independent variables). Hence

$$
E\Bigl[\sum (X_i-\bar X)^2\Bigr] = n\sigma^2 - n\cdot\frac{\sigma^2}{n} = (n-1)\sigma^2 . \qquad\square
$$

The intuition: deviations are measured from $\bar x$, which is itself fitted to the data and therefore sits closer to the points than the true $\mu$ does. One "degree of freedom" is used up in estimating $\mu$ — the $n$ deviations satisfy one linear constraint $\sum(x_i-\bar x)=0$, so only $n-1$ are free. The same counting gives the divisor $n-p$ when $p$ parameters are fitted in regression (Section 20.9).

Note that $s$ itself is *not* unbiased for $\sigma$: by Jensen's inequality $E[S]<\sigma$. For Gaussian data $E[S]=c_4\sigma$ with $c_4 = \sqrt{2/(n-1)}\,\Gamma(n/2)/\Gamma((n-1)/2)$, e.g. $c_4=0.798$ for $n=2$ and $0.973$ for $n=10$. The bias is usually ignored except in quality-control charts.

### Computational formulas and numerical stability

Expanding the square gives the "computational" formula

$$
s^2 = \frac{1}{n-1}\Bigl(\sum x_i^2 - n\bar x^2\Bigr) = \frac{1}{n-1}\Bigl(\sum x_i^2 - \frac{(\sum x_i)^2}{n}\Bigr),
$$

which needs only running sums. It is exact in real arithmetic but disastrous in floating point when $\lvert\bar x\rvert\gg s$: two large, nearly equal numbers are subtracted, and the relative rounding error of each (about $10^{-16}$ in double precision) is multiplied by $\bar x^2/s^2$ (catastrophic cancellation, Ch. 22).

The stable alternatives are (i) the **two-pass** algorithm (compute $\bar x$ first, then $\sum(x_i-\bar x)^2$), (ii) **shifting** the data by any value close to the mean before using the one-pass formula, or (iii) **Welford's online algorithm**, which updates mean and sum of squared deviations sample by sample:

$$
\begin{aligned}
\delta &= x_k - \bar x_{k-1}, \qquad \bar x_k = \bar x_{k-1} + \frac{\delta}{k},\\
M_k &= M_{k-1} + \delta\,(x_k-\bar x_k), \qquad s_k^2 = \frac{M_k}{k-1}.
\end{aligned}
$$

Welford's recursion is what a well-written DMM "statistics" mode or a streaming data logger should use.

### Variance, RMS and power

For an electrical signal sampled as $v_1,\dots,v_n$, the mean is the DC component and the (population-divisor) variance is the AC power in a $1\ \Omega$ reference. Since

$$
\frac1n\sum v_i^2 = \bar v^{\,2} + \frac1n\sum (v_i-\bar v)^2 ,
$$

the squared RMS value equals the squared DC value plus the squared standard deviation:

$$
V_\mathrm{rms}^2 = V_\mathrm{DC}^2 + \sigma_v^2 .
$$

A "true-RMS AC" meter that blocks DC reports $\sigma_v$; a DC-coupled true-RMS meter reports $\sqrt{V_\mathrm{DC}^2+\sigma_v^2}$. This identity is the statistical form of Parseval's theorem (Ch. 15) restricted to the zero-frequency term.

> [!example] Worked Example 20.3 — Mean and standard deviation of ADC readings
> Eight readings of a nominal $2.5\ \mathrm{V}$ reference are
> $2.47,\ 2.53,\ 2.51,\ 2.49,\ 2.55,\ 2.45,\ 2.50,\ 2.52\ \mathrm{V}$. Find $\bar x$, $s$, and the RMS value.
>
> **Step 1 — mean.** $\sum x_i = 20.02$, so $\bar x = 2.5025\ \mathrm{V}$.
>
> **Step 2 — deviations.** $x_i-\bar x = -0.0325,\ 0.0275,\ 0.0075,\ -0.0125,\ 0.0475,\ -0.0525,\ -0.0025,\ 0.0175$ (they sum to zero — a free arithmetic check).
>
> **Step 3 — sum of squares.** $\sum(x_i-\bar x)^2 = 0.00735\ \mathrm{V^2}$.
>
> **Step 4 — variance and standard deviation.** $s^2 = 0.00735/7 = 1.05\times10^{-3}\ \mathrm{V^2}$, $s = 32.4\ \mathrm{mV}$. (With divisor $n$ one would get $30.3\ \mathrm{mV}$; the difference is $\sqrt{8/7}$, i.e. about $7\%$ for $n=8$.)
>
> **Step 5 — RMS.** $\frac1n\sum x_i^2 = 50.1074/8 = 6.263425\ \mathrm{V^2}$, and indeed $\bar x^2 + \frac{1}{n}\sum(x_i-\bar x)^2 = 6.26250625+0.00091875 = 6.263425$. So $V_\mathrm{rms}=2.50268\ \mathrm{V}$, dominated by the DC value as expected.

> [!example] Worked Example 20.4 — When the computational formula fails
> A frequency counter logs four readings $10^9+4,\ 10^9+7,\ 10^9+13,\ 10^9+16$ (in Hz). The true sample variance is that of $4,7,13,16$: $\bar d=10$, $\sum(d_i-\bar d)^2 = 36+9+9+36=90$, $s^2 = 30\ \mathrm{Hz^2}$.
>
> Evaluating $\bigl(\sum x_i^2-(\sum x_i)^2/n\bigr)/(n-1)$ in IEEE double precision returns $-170.7\ \mathrm{Hz^2}$ — a *negative* variance. Each $x_i^2\approx10^{18}$ is stored with an absolute rounding error of about $10^{18}\times 1.1\times10^{-16}\approx 100$, comparable to the 90 we are trying to recover. Subtracting $10^9$ first (or using the two-pass or Welford algorithm) returns exactly $30\ \mathrm{Hz^2}$.

> [!example] Worked Example 20.5 — Combining three calibrations
> Three laboratories report the value of a standard resistor as $10.12\pm0.02$, $10.05\pm0.05$ and $10.09\pm0.03\ \Omega$ (standard uncertainties). Find the weighted mean and its uncertainty.
>
> **Step 1 — weights.** $w_i=1/u_i^2 = 2500,\ 400,\ 1111.1\ \Omega^{-2}$; $\sum w_i=4011.1$.
>
> **Step 2 — weighted mean.** $\bar x_w = (2500\times10.12+400\times10.05+1111.1\times10.09)/4011.1 = 10.1047\ \Omega$.
>
> **Step 3 — uncertainty.** $u(\bar x_w) = 1/\sqrt{4011.1} = 0.0158\ \Omega$ — smaller than the best single uncertainty, as it should be.
>
> **Step 4 — consistency check.** $\chi^2 = \sum w_i(x_i-\bar x_w)^2 = 0.585+1.197+0.240 = 2.02$ with $n-1=2$ degrees of freedom. A value near the degrees of freedom indicates the three results agree within their stated uncertainties (Section 20.8); a value much larger would suggest an unrecognised error in at least one laboratory.

> [!warning] Common Pitfall — Quoting $s$ when you mean the uncertainty of the mean
> The standard deviation $s$ describes the scatter of *individual* readings. The uncertainty of their *average* is the standard error $s/\sqrt n$ (Section 20.5). In Example 20.3, the readings scatter by $32\ \mathrm{mV}$ but the mean is known to about $32/\sqrt8\approx11\ \mathrm{mV}$. Confusing the two either overstates or understates the result by a factor $\sqrt n$.

> [!tip] Engineering Note — Robust statistics in firmware
> Embedded systems reading noisy sensors often use a *median-of-$k$* filter (typically $k=3$ or $5$) to reject single-sample glitches from ESD or relay switching, followed by an averaging filter to reduce Gaussian noise. The median removes impulses that would corrupt a mean; the subsequent mean recovers the efficiency the median loses.

## 20.3 Covariance and correlation

### Motivation

Engineering data are rarely one-dimensional. A sensor's output depends on temperature; two receiver channels share a local oscillator and therefore share its phase noise; the gain and offset of an ADC drift together with supply voltage. Quantifying how two variables vary *together* is needed to decide whether one predicts the other, to propagate uncertainties correctly when inputs are not independent (Section 20.10), and to design systems that cancel common-mode effects.

### Population and sample covariance

Recall from Ch. 19 that for random variables $X$ and $Y$ with means $\mu_X,\mu_Y$,

$$
\operatorname{Cov}(X,Y) = E\bigl[(X-\mu_X)(Y-\mu_Y)\bigr] = E[XY]-\mu_X\mu_Y, \qquad \rho_{XY} = \frac{\operatorname{Cov}(X,Y)}{\sigma_X\sigma_Y}.
$$

> [!definition] Sample covariance and Pearson correlation coefficient
> For paired data $(x_i,y_i)$, $i=1,\dots,n$, define the sums of squares and products
> $$
> S_{xx}=\sum (x_i-\bar x)^2,\qquad S_{yy}=\sum (y_i-\bar y)^2,\qquad S_{xy}=\sum (x_i-\bar x)(y_i-\bar y).
> $$
> The **sample covariance** is $s_{xy}=S_{xy}/(n-1)$ and the **sample (Pearson) correlation coefficient** is
> $$
> r = \frac{s_{xy}}{s_x s_y} = \frac{S_{xy}}{\sqrt{S_{xx}S_{yy}}} .
> $$

The sign of $s_{xy}$ says whether large $x$ tends to go with large $y$ (positive) or small $y$ (negative). Its magnitude depends on units (Ω·°C, say), which is why the dimensionless $r$ is preferred for describing the *strength* of association.

> [!theorem] Bounds on the correlation coefficient
> $-1\le r\le 1$, with $\lvert r\rvert=1$ if and only if all points lie exactly on a straight line $y=a+bx$ with $b\neq0$ (sign of $r$ = sign of $b$). The same holds for $\rho_{XY}$, with "exactly" replaced by "with probability one".

*Proof.* Let $\mathbf{u}=(x_i-\bar x)_i$ and $\mathbf{v}=(y_i-\bar y)_i$ be the centred data vectors in $\mathbb{R}^n$. Then $r=\langle\mathbf{u},\mathbf{v}\rangle/(\lVert\mathbf{u}\rVert\lVert\mathbf{v}\rVert)$ is the cosine of the angle between them. The Cauchy–Schwarz inequality (Ch. 9) gives $\lvert\langle\mathbf{u},\mathbf{v}\rangle\rvert\le\lVert\mathbf{u}\rVert\lVert\mathbf{v}\rVert$, with equality iff $\mathbf{v}=b\mathbf{u}$, i.e. $y_i-\bar y=b(x_i-\bar x)$ for all $i$. $\square$

This geometric picture — **correlation is the cosine of the angle between centred data vectors** — is worth remembering: $r=0$ means the centred vectors are orthogonal, and $r^2$ is the squared cosine, the fraction of $\lVert\mathbf v\rVert^2$ captured by projecting $\mathbf v$ onto $\mathbf u$. This is exactly the $r^2$ of simple linear regression (Section 20.9).

Further properties follow directly from the bilinearity of covariance:

- $\operatorname{Cov}(aX+b,\,cY+d)=ac\operatorname{Cov}(X,Y)$, so $\rho$ is invariant to changes of units and origin (up to sign if $ac<0$).
- $\operatorname{Var}(X\pm Y)=\operatorname{Var}X+\operatorname{Var}Y\pm2\operatorname{Cov}(X,Y)$.
- Independent $\Rightarrow$ uncorrelated; the converse fails in general but holds for jointly Gaussian variables (Ch. 19).

### What correlation does not measure

Pearson's $r$ measures **linear** association only. If $X$ is symmetric about zero and $Y=X^2$, then $Y$ is a deterministic function of $X$, yet $\operatorname{Cov}(X,Y)=E[X^3]-E[X]E[X^2]=0$, so $\rho=0$. Conversely, a strong $r$ can be produced by a single extreme point, by mixing two groups with different means, or by a common cause (both $x$ and $y$ rise with ambient temperature). The four data sets of *Anscombe's quartet* share identical means, variances, $r=0.816$ and regression lines while looking completely different when plotted. The rule is simple: never interpret $r$ without a scatter plot.

> [!warning] Common Pitfall — Correlation is not causation, and not agreement
> A correlation between solder-joint failures and production day of week may reflect a third variable (which operator was on shift). Separately, two instruments measuring the same voltage can have $r=0.999$ while one reads consistently $2\%$ high: correlation measures association, not agreement. For method comparison, analyse the *differences* $y_i-x_i$ (mean bias and its spread) instead.

### Rank correlation

When the relation is monotonic but not linear, or the data contain outliers or are ordinal, replace values by their ranks.

> [!definition] Spearman's rank correlation
> Let $R_i$ and $Q_i$ be the ranks of $x_i$ among the $x$'s and of $y_i$ among the $y$'s (average ranks for ties). **Spearman's** $r_s$ is the Pearson correlation of the ranks. With no ties it simplifies to
> $$
> r_s = 1 - \frac{6\sum_i d_i^2}{n(n^2-1)}, \qquad d_i = R_i - Q_i .
> $$

The shortcut follows from $\sum R_i=\sum Q_i=n(n+1)/2$ and $\sum R_i^2=\sum Q_i^2=n(n+1)(2n+1)/6$, which fix the means and variances of the ranks; expanding $\sum d_i^2=\sum R_i^2+\sum Q_i^2-2\sum R_iQ_i$ and substituting gives the formula. $r_s=\pm1$ for any strictly monotonic relation, linear or not.

### The covariance matrix

For a random vector $\mathbf{X}=(X_1,\dots,X_p)^\mathsf{T}$ with mean $\boldsymbol{\mu}$, the **covariance matrix** is

$$
\boldsymbol{\Sigma} = E\bigl[(\mathbf X-\boldsymbol\mu)(\mathbf X-\boldsymbol\mu)^\mathsf{T}\bigr], \qquad \Sigma_{ij}=\operatorname{Cov}(X_i,X_j).
$$

Its sample counterpart, from $n$ observation vectors $\mathbf{x}_k$, is $\mathbf{S}=\frac{1}{n-1}\sum_k(\mathbf x_k-\bar{\mathbf x})(\mathbf x_k-\bar{\mathbf x})^\mathsf T$.

> [!theorem] Properties of the covariance matrix
> 1. $\boldsymbol\Sigma$ is symmetric and positive semidefinite: $\mathbf a^\mathsf T\boldsymbol\Sigma\mathbf a=\operatorname{Var}(\mathbf a^\mathsf T\mathbf X)\ge0$ for every $\mathbf a$.
> 2. For a linear map $\mathbf Y=\mathbf A\mathbf X+\mathbf b$: $\boldsymbol\Sigma_Y=\mathbf A\boldsymbol\Sigma_X\mathbf A^\mathsf T$.
> 3. Its eigenvectors (Ch. 10) are orthogonal **principal axes**; the eigenvalues are the variances along those axes, and $\operatorname{tr}\boldsymbol\Sigma=\sum_i\sigma_i^2$ is the total variance.

*Proof of 1 and 2.* $\mathbf a^\mathsf T\mathbf X$ is a scalar random variable with mean $\mathbf a^\mathsf T\boldsymbol\mu$, so its variance is $E[(\mathbf a^\mathsf T(\mathbf X-\boldsymbol\mu))^2] = \mathbf a^\mathsf TE[(\mathbf X-\boldsymbol\mu)(\mathbf X-\boldsymbol\mu)^\mathsf T]\mathbf a$, which is non-negative because it is a variance. For 2, $\mathbf Y-E\mathbf Y=\mathbf A(\mathbf X-\boldsymbol\mu)$, and substituting into the definition gives $\mathbf A\boldsymbol\Sigma_X\mathbf A^\mathsf T$. Property 3 is the spectral theorem for real symmetric matrices. $\square$

Property 2 is the matrix form of uncertainty propagation used in Section 20.10. Property 3 underlies principal component analysis: rotating to the eigenvector basis produces *uncorrelated* components, and discarding those with small eigenvalues compresses the data with minimal loss of variance. The **correlation matrix** is $\mathbf{R}=\mathbf D^{-1}\boldsymbol\Sigma\mathbf D^{-1}$ with $\mathbf D=\operatorname{diag}(\sigma_1,\dots,\sigma_p)$; it has ones on the diagonal.

> [!example] Worked Example 20.6 — Correlation of an RTD calibration
> A platinum RTD is calibrated at six temperatures:
>
> | $T$ (°C) | 0 | 20 | 40 | 60 | 80 | 100 |
> |---|---|---|---|---|---|---|
> | $R$ (Ω) | 100.0 | 107.9 | 115.5 | 123.2 | 130.9 | 138.5 |
>
> Find $s_{TR}$ and $r$.
>
> **Step 1 — means.** $\bar T=50\ \mathrm{^\circ C}$, $\bar R=716/6=119.333\ \Omega$.
>
> **Step 2 — sums of squares.** $S_{TT}=2(50^2+30^2+10^2)=7000\ \mathrm{(^\circ C)^2}$; $S_{RR}=1035.29\ \Omega^2$.
>
> **Step 3 — sum of products.** Because $\sum(T_i-\bar T)=0$, $S_{TR}=\sum (T_i-\bar T)R_i = -50(100.0)-30(107.9)-10(115.5)+10(123.2)+30(130.9)+50(138.5) = 2692\ \Omega\,\mathrm{^\circ C}$.
>
> **Step 4 — covariance and correlation.** $s_{TR}=2692/5=538.4\ \Omega\,\mathrm{^\circ C}$; $r = 2692/\sqrt{7000\times1035.29}=0.99999$.
>
> **Interpretation.** $r$ is so close to 1 that it no longer discriminates between "linear" and "slightly curved". The physically known curvature of platinum (Callendar–Van Dusen equation) is invisible in $r$; it shows up only in the residuals of a fit (Section 20.9). A high $r$ is necessary but not sufficient for a good linear model.

> [!example] Worked Example 20.7 — Rank correlation of capacitor ageing
> Six electrolytic capacitors have equivalent series resistance (ESR, in mΩ) and measured remaining life (kh):
>
> | ESR | 3.1 | 2.4 | 4.8 | 1.9 | 3.6 | 5.2 |
> |---|---|---|---|---|---|---|
> | Life | 12 | 15 | 7 | 18 | 9 | 8 |
>
> **Step 1 — rank each variable.** ESR ranks: $3,2,5,1,4,6$. Life ranks: $4,5,1,6,3,2$.
>
> **Step 2 — rank differences.** $d_i = -1,-3,4,-5,1,4$; $\sum d_i^2 = 1+9+16+25+1+16=68$.
>
> **Step 3 — Spearman.** $r_s = 1-6(68)/(6\cdot35) = 1-408/210 = -0.943$.
>
> **Step 4 — compare.** Pearson's $r$ on the raw values is $-0.939$. The two agree here because the relation is roughly linear; had one capacitor failed early with a huge ESR, $r$ would be pulled strongly while $r_s$ would change by at most one rank.

> [!example] Worked Example 20.8 — Sum and difference of correlated noise
> Two amplifier channels have output noise with $\sigma_1 = 2\ \mathrm{mV}$, $\sigma_2=3\ \mathrm{mV}$ and correlation $\rho=0.5$ (shared supply ripple). Find the noise of the sum and of the difference, and the principal axes.
>
> **Step 1 — covariance matrix.** $\operatorname{Cov}=\rho\sigma_1\sigma_2 = 3\ \mathrm{mV^2}$, so
> $$
> \boldsymbol\Sigma=\begin{pmatrix}4&3\\3&9\end{pmatrix}\ \mathrm{mV^2}.
> $$
>
> **Step 2 — sum.** $\mathbf a=(1,1)^\mathsf T$: $\operatorname{Var}=\mathbf a^\mathsf T\boldsymbol\Sigma\mathbf a = 4+9+2(3)=19\ \mathrm{mV^2}$, $\sigma_{1+2}=4.36\ \mathrm{mV}$.
>
> **Step 3 — difference.** $\mathbf a=(1,-1)^\mathsf T$: $\operatorname{Var}=4+9-6=7\ \mathrm{mV^2}$, $\sigma_{1-2}=2.65\ \mathrm{mV}$. A differential measurement partially cancels the common component — the statistical basis of common-mode rejection.
>
> **Step 4 — principal axes.** The eigenvalues solve $\lambda^2-13\lambda+27=0$: $\lambda=(13\pm\sqrt{61})/2 = 10.41,\ 2.59\ \mathrm{mV^2}$. The unit-norm combination with the least variance is the eigenvector of $\lambda=2.59$, approximately $(0.91,-0.42)^\mathsf T$, with $\sigma=1.61\ \mathrm{mV}$. (For comparison, the normalised plain difference $(1,-1)^\mathsf T/\sqrt2$ has variance $3.5\ \mathrm{mV^2}$.) Weighting the quieter channel more heavily cancels more of the shared ripple.

> [!tip] Engineering Note — Covariance in multi-channel systems
> In antenna arrays, the sample covariance matrix of the element signals is the starting point of adaptive beamforming and direction-of-arrival estimation: its large eigenvalues correspond to signals and its small ones to noise. In instrumentation, the off-diagonal terms decide whether averaging two channels helps (uncorrelated noise, $\sigma/\sqrt2$) or barely helps (strongly correlated noise).

## 20.4 Common distributions

### Why a catalogue

Statistical inference needs models, and a surprisingly small set of distributions covers most engineering situations. Each arises from a specific mechanism — counting successes, counting rare events, waiting for the next event, adding many small effects, taking the magnitude of a 2-D Gaussian vector, summing squared Gaussian errors — and recognising the mechanism tells you which model to use. The last three entries ($\chi^2$, $t$, $F$) are less physical models of data than *sampling distributions*: they describe statistics computed from Gaussian samples and are the workhorses of Sections 20.7–20.8.

General definitions (pmf, pdf, CDF, expectation, moment generating function $M_X(s)=E[e^{sX}]$) are in Ch. 19. For each distribution below we give the mechanism, the pmf or pdf, the mean and variance, and where it appears in EE.

### Bernoulli distribution

> [!definition] Bernoulli$(p)$
> $X\in\{0,1\}$ with $P(X=1)=p$, $P(X=0)=1-p$. Then $E[X]=p$, $\operatorname{Var}(X)=p(1-p)$, $M_X(s)=1-p+pe^{s}$.

The variance follows from $X^2=X$: $E[X^2]=p$, so $\operatorname{Var}X=p-p^2$. It is maximal, $1/4$, at $p=1/2$ — the worst case used for sample-size planning of proportions (Section 20.7). Models: one transmitted bit is in error or not; one device passes or fails test; one photon is detected or not.

### Binomial distribution

> [!definition] Binomial$(n,p)$
> The number of successes in $n$ independent Bernoulli$(p)$ trials:
> $$
> P(X=k)=\binom{n}{k}p^k(1-p)^{n-k},\qquad k=0,1,\dots,n,
> $$
> with $E[X]=np$, $\operatorname{Var}(X)=np(1-p)$, $M_X(s)=(1-p+pe^s)^n$.

*Derivation.* Any particular sequence with $k$ successes has probability $p^k(1-p)^{n-k}$, and there are $\binom nk$ such sequences. Writing $X=\sum_{i=1}^n B_i$ with independent Bernoulli $B_i$, linearity gives $E[X]=np$ and independence gives $\operatorname{Var}X=np(1-p)$; the MGF is the $n$th power of the Bernoulli MGF because MGFs of independent sums multiply. The pmf is unimodal with mode at $\lfloor (n+1)p\rfloor$.

### Poisson distribution

> [!definition] Poisson$(\lambda)$
> $$
> P(X=k)=\frac{\lambda^k e^{-\lambda}}{k!},\qquad k=0,1,2,\dots,
> $$
> with $E[X]=\operatorname{Var}(X)=\lambda$ and $M_X(s)=\exp\bigl(\lambda(e^s-1)\bigr)$.

> [!theorem] Poisson limit of the binomial
> If $n\to\infty$ and $p\to0$ with $np=\lambda$ fixed, then $\binom nk p^k(1-p)^{n-k}\to \lambda^ke^{-\lambda}/k!$.

*Proof sketch.* $\binom nk p^k = \frac{n(n-1)\cdots(n-k+1)}{n^k}\frac{\lambda^k}{k!}\to\frac{\lambda^k}{k!}$, and $(1-\lambda/n)^{n-k}\to e^{-\lambda}$. Equivalently, the binomial MGF $(1+p(e^s-1))^n=(1+\lambda(e^s-1)/n)^n\to e^{\lambda(e^s-1)}$. $\square$

So the Poisson law describes the count of *rare, independent* events in a fixed interval: electrons crossing a junction (shot noise), photons on a detector, radioactive decays, packet arrivals on a lightly loaded link, soft errors in a memory over a day. The equality of mean and variance is a useful diagnostic: count data whose variance greatly exceeds the mean are *overdispersed* (clustered events), and a plain Poisson model will underestimate their variability.

```desmos-graph
left=-0.5; right=10.5; bottom=-0.02; top=0.3
---
(0,0.0282)|blue
(1,0.1211)|blue
(2,0.2335)|blue
(3,0.2668)|blue
(4,0.2001)|blue
(5,0.1029)|blue
(6,0.0368)|blue
(7,0.0090)|blue
(8,0.0014)|blue
(9,0.0001)|blue
(10,0.0000)|blue
(0,0.0498)|red
(1,0.1494)|red
(2,0.2240)|red
(3,0.2240)|red
(4,0.1680)|red
(5,0.1008)|red
(6,0.0504)|red
(7,0.0216)|red
(8,0.0081)|red
(9,0.0027)|red
(10,0.0008)|red
```

*Binomial$(10,0.3)$ pmf (blue) and Poisson$(3)$ pmf (red): equal means, but the Poisson has variance 3 against the binomial's 2.1, so it is flatter with a longer right tail.*

### Uniform distribution

> [!definition] Uniform$(a,b)$
> $$
> f_X(x)=\frac{1}{b-a}\quad (a\le x\le b),\qquad E[X]=\frac{a+b}{2},\qquad \operatorname{Var}(X)=\frac{(b-a)^2}{12}.
> $$

*Derivation of the variance.* Centre at zero with half-width $c=(b-a)/2$: $\operatorname{Var}X=\frac{1}{2c}\int_{-c}^{c}x^2\,dx=\frac{c^2}{3}=\frac{(b-a)^2}{12}$.

The uniform law models complete ignorance within known bounds. Two cases matter enormously in EE. **Quantisation error** of an ideal ADC with step $\Delta$ is (under mild conditions) uniform on $[-\Delta/2,\Delta/2]$, with variance $\Delta^2/12$ — the source of the familiar $6.02N+1.76\ \mathrm{dB}$ SNR formula. **Type B uncertainty** from a manufacturer's "$\pm a$" specification is modelled as uniform on $[-a,a]$, giving standard uncertainty $a/\sqrt3$ (Section 20.10). The random phase of an unsynchronised oscillator is uniform on $[0,2\pi)$.

### Exponential distribution

> [!definition] Exponential$(\lambda)$
> $$
> f_X(x)=\lambda e^{-\lambda x}\ (x\ge0),\qquad F_X(x)=1-e^{-\lambda x},\qquad E[X]=\frac1\lambda,\qquad \operatorname{Var}(X)=\frac{1}{\lambda^2}.
> $$

It is the waiting time to the first event of a Poisson process of rate $\lambda$: $P(X>x)=P(\text{no events in }[0,x])=e^{-\lambda x}$. It is the only continuous distribution with the **memoryless property**

$$
P(X>s+t\mid X>s)=\frac{e^{-\lambda(s+t)}}{e^{-\lambda s}}=e^{-\lambda t}=P(X>t).
$$

In reliability engineering, $\lambda$ is a constant **failure rate** (failures per hour, or FIT = failures per $10^9$ device-hours) and $1/\lambda$ is the **MTBF** (mean time between failures). Memorylessness means a working part is "as good as new" — appropriate for the flat bottom of the bathtub curve (random failures), wrong for wear-out (electrolytic capacitors, fans, relays), which is modelled by the Weibull distribution instead.

```desmos-graph
left=0; right=5; bottom=0; top=2.1
---
y=0.5e^{-0.5x}|blue|0<x<5
y=e^{-x}|red|0<x<5
y=2e^{-2x}|green|0<x<5
```

*Exponential pdfs for $\lambda=0.5$ (blue), $1$ (red) and $2$ (green). The pdf starts at $\lambda$ and decays with time constant $1/\lambda$ — the same shape as an RC discharge.*

### Gaussian (normal) distribution

> [!definition] Gaussian $\mathcal N(\mu,\sigma^2)$
> $$
> f_X(x)=\frac{1}{\sigma\sqrt{2\pi}}\exp\Bigl(-\frac{(x-\mu)^2}{2\sigma^2}\Bigr),\qquad E[X]=\mu,\quad\operatorname{Var}(X)=\sigma^2,\quad M_X(s)=e^{\mu s+\sigma^2s^2/2}.
> $$
> The **standard** Gaussian $Z\sim\mathcal N(0,1)$ has CDF $\Phi(z)$; engineers usually use the tail function $Q(z)=1-\Phi(z)=P(Z>z)$.

Any Gaussian is standardised by $Z=(X-\mu)/\sigma$, so $P(X>x)=Q\bigl((x-\mu)/\sigma\bigr)$. Key values: $P(\lvert Z\rvert\le1)=0.6827$, $P(\lvert Z\rvert\le2)=0.9545$, $P(\lvert Z\rvert\le3)=0.9973$ (the "68–95–99.7 rule"); $Q(1.645)=0.05$, $Q(1.960)=0.025$, $Q(2.576)=0.005$, $Q(3)=1.35\times10^{-3}$. For large $z$, $Q(z)\approx \frac{1}{z\sqrt{2\pi}}e^{-z^2/2}$.

Three facts make the Gaussian central:

1. **Closure under linear operations.** Any linear combination of jointly Gaussian variables is Gaussian; the output of an LTI system driven by Gaussian noise is Gaussian (Ch. 21).
2. **The central limit theorem** (Section 20.5): sums of many small independent effects are approximately Gaussian, which is why thermal noise, aggregate manufacturing variation and many measurement errors are Gaussian.
3. **Maximum entropy.** Among all distributions with given mean and variance, the Gaussian has the largest entropy — it is the "least committal" model consistent with knowing only two moments.

```desmos-graph
left=-6; right=6; bottom=0; top=0.85
---
y=\frac{1}{\sqrt{2\pi}}e^{-\frac{x^{2}}{2}}|blue
y=\frac{1}{0.5\sqrt{2\pi}}e^{-\frac{x^{2}}{2\cdot0.25}}|red
y=\frac{1}{2\sqrt{2\pi}}e^{-\frac{(x-1)^{2}}{8}}|green
```

*Gaussian pdfs: $\mathcal N(0,1)$ (blue), $\mathcal N(0,0.25)$ (red) and $\mathcal N(1,4)$ (green). Halving $\sigma$ doubles the peak height; the area stays 1.*

### Rayleigh distribution

> [!definition] Rayleigh$(\sigma)$
> If $X,Y$ are independent $\mathcal N(0,\sigma^2)$, then $R=\sqrt{X^2+Y^2}$ has
> $$
> f_R(r)=\frac{r}{\sigma^2}e^{-r^2/(2\sigma^2)}\ (r\ge0),\qquad P(R>r)=e^{-r^2/(2\sigma^2)},
> $$
> $$
> E[R]=\sigma\sqrt{\pi/2}\approx1.253\sigma,\qquad E[R^2]=2\sigma^2,\qquad\operatorname{Var}(R)=\Bigl(2-\frac\pi2\Bigr)\sigma^2\approx0.429\sigma^2 .
> $$

*Derivation.* In polar coordinates $(r,\phi)$ the joint pdf $\frac{1}{2\pi\sigma^2}e^{-(x^2+y^2)/(2\sigma^2)}$ becomes, with the Jacobian $r$ (Ch. 6), $\frac{r}{2\pi\sigma^2}e^{-r^2/(2\sigma^2)}$. Integrating out $\phi\in[0,2\pi)$ gives $f_R$; the phase is uniform and independent of $R$. The moments follow from $\int_0^\infty r^2e^{-r^2/(2\sigma^2)}dr=\sigma^3\sqrt{\pi/2}$.

The Rayleigh law is the **envelope of narrowband Gaussian noise** (in-phase and quadrature components are independent Gaussians) and the amplitude of a **multipath-faded** radio signal with no line-of-sight path. $R^2$ is exponential with mean $2\sigma^2$, so the *power* of a Rayleigh-faded signal is exponentially distributed. Adding a deterministic component gives the Rician distribution.

```desmos-graph
left=0; right=8; bottom=0; top=1.3
---
y=xe^{-\frac{x^{2}}{2}}|blue|0<x<8
y=\frac{x}{4}e^{-\frac{x^{2}}{8}}|red|0<x<8
y=\frac{x}{0.25}e^{-\frac{x^{2}}{0.5}}|green|0<x<8
```

*Rayleigh pdfs for $\sigma=1$ (blue), $\sigma=2$ (red) and $\sigma=0.5$ (green). The mode is at $r=\sigma$; the long right tail is the occasional noise peak that triggers a false alarm.*

### Chi-square distribution

> [!definition] Chi-square $\chi^2_k$
> If $Z_1,\dots,Z_k$ are independent $\mathcal N(0,1)$, then $Q=\sum_{i=1}^k Z_i^2$ has the chi-square distribution with $k$ **degrees of freedom**:
> $$
> f_Q(q)=\frac{q^{k/2-1}e^{-q/2}}{2^{k/2}\Gamma(k/2)}\ (q>0),\qquad E[Q]=k,\qquad\operatorname{Var}(Q)=2k,\qquad M_Q(s)=(1-2s)^{-k/2}.
> $$

*Derivation of the moments.* For one term, $E[Z^2]=1$ and $E[Z^4]=3$, so $\operatorname{Var}(Z^2)=2$; adding $k$ independent terms gives mean $k$ and variance $2k$. The MGF of $Z^2$ is $\int\frac{1}{\sqrt{2\pi}}e^{sz^2-z^2/2}dz=(1-2s)^{-1/2}$, and the $k$th power gives $M_Q$, which is the MGF of a gamma distribution with shape $k/2$ and scale 2 — hence the pdf.

Special cases tie the catalogue together: $\chi^2_2$ is exponential with mean 2, and $\sqrt{\chi^2_2}\cdot\sigma$ is Rayleigh. $\chi^2_k$ is skewed right for small $k$ and approaches $\mathcal N(k,2k)$ as $k\to\infty$ (CLT). In EE it describes **noise power estimates** (a sum of squared Gaussian samples), **energy detectors** in spectrum sensing, and, most importantly for this chapter, the sample variance of Gaussian data (Section 20.5) and goodness-of-fit statistics (Section 20.8).

```desmos-graph
left=0; right=16; bottom=0; top=0.55
---
y=\frac{1}{\sqrt{2\pi x}}e^{-\frac{x}{2}}|blue|0<x<16
y=\frac{1}{2}e^{-\frac{x}{2}}|red|0<x<16
y=\frac{x}{4}e^{-\frac{x}{2}}|green|0<x<16
y=\frac{x^{2}}{16}e^{-\frac{x}{2}}|purple|0<x<16
```

*Chi-square pdfs for $k=1$ (blue), $2$ (red), $4$ (green) and $6$ (purple). For $k\ge3$ the mode is at $k-2$; the distributions spread and become more symmetric as $k$ grows.*

### Student's $t$ distribution

> [!definition] Student's $t_\nu$
> If $Z\sim\mathcal N(0,1)$ and $V\sim\chi^2_\nu$ are independent, then $T=Z/\sqrt{V/\nu}$ has the $t$ distribution with $\nu$ degrees of freedom:
> $$
> f_T(t)=\frac{\Gamma\bigl(\frac{\nu+1}{2}\bigr)}{\sqrt{\nu\pi}\,\Gamma\bigl(\frac\nu2\bigr)}\Bigl(1+\frac{t^2}{\nu}\Bigr)^{-(\nu+1)/2},
> $$
> with $E[T]=0$ ($\nu>1$) and $\operatorname{Var}(T)=\nu/(\nu-2)$ ($\nu>2$).

The $t$ distribution is what a standardised Gaussian becomes when the unknown $\sigma$ is replaced by an estimate from $\nu$ degrees of freedom. The extra randomness of the estimate fattens the tails: with $\nu=3$, $P(\lvert T\rvert>3)=0.058$, versus $0.0027$ for the Gaussian — a factor of 21. As $\nu\to\infty$, $V/\nu\to1$ and $t_\nu\to\mathcal N(0,1)$. The two-sided $95\%$ critical value $t_{0.025,\nu}$ is $12.71$ for $\nu=1$, $3.18$ for $\nu=3$, $2.57$ for $\nu=5$, $2.23$ for $\nu=10$ and $2.04$ for $\nu=30$, against $1.96$ for the Gaussian. The case $\nu=1$ is the Cauchy distribution, which has no mean at all: averaging Cauchy samples does not reduce their spread.

```desmos-graph
left=-5; right=5; bottom=0; top=0.45
---
y=\frac{1}{\sqrt{2\pi}}e^{-\frac{x^{2}}{2}}|blue
y=\frac{1}{\pi(1+x^{2})}|red
y=0.36755\left(1+\frac{x^{2}}{3}\right)^{-2}|green
```

*Student $t$ pdfs for $\nu=1$ (red, Cauchy) and $\nu=3$ (green) against $\mathcal N(0,1)$ (blue). Lower peaks, heavier tails: small samples must use wider intervals.*

### The $F$ distribution

> [!definition] Fisher–Snedecor $F_{d_1,d_2}$
> If $U\sim\chi^2_{d_1}$ and $V\sim\chi^2_{d_2}$ are independent, $F=\dfrac{U/d_1}{V/d_2}$ has pdf
> $$
> f(x)=\frac{(d_1/d_2)^{d_1/2}}{B(d_1/2,d_2/2)}\,x^{d_1/2-1}\Bigl(1+\frac{d_1}{d_2}x\Bigr)^{-(d_1+d_2)/2}\quad(x>0),
> $$
> with $E[F]=d_2/(d_2-2)$ for $d_2>2$; $B$ is the beta function.

$F$ is the distribution of a **ratio of two independent variance estimates** from Gaussian data, normalised by their degrees of freedom; it is used to compare the noise of two instruments (Section 20.8) and to test the significance of regression models (Section 20.9). Useful identities: $1/F_{d_1,d_2}\sim F_{d_2,d_1}$, so lower critical values follow from upper ones, and $t_\nu^2\sim F_{1,\nu}$.

```desmos-graph
left=0; right=5; bottom=0; top=1
---
y=10.3684x^{1.5}\left(1+0.5x\right)^{-7.5}|blue|0<x<5
y=312.8125x^{4}\left(1+0.5x\right)^{-15}|red|0<x<5
```

*$F$ pdfs for $(d_1,d_2)=(5,10)$ (blue) and $(10,20)$ (red). Both concentrate near 1 — the value expected when two variance estimates agree — and tighten as the degrees of freedom grow.*

### Relations between the distributions

| From | Operation | Gives |
|---|---|---|
| $n$ independent Bernoulli$(p)$ | sum | Binomial$(n,p)$ |
| Binomial$(n,p)$ | $n\to\infty$, $np=\lambda$ | Poisson$(\lambda)$ |
| Binomial, Poisson | large $np(1-p)$, large $\lambda$ | Gaussian (CLT) |
| Poisson process rate $\lambda$ | time to first event | Exponential$(\lambda)$ |
| Two independent $\mathcal N(0,\sigma^2)$ | magnitude $\sqrt{X^2+Y^2}$ | Rayleigh$(\sigma)$ |
| $k$ independent $\mathcal N(0,1)$ | sum of squares | $\chi^2_k$ |
| $\chi^2_2/2$ | — | Exponential$(1)$ |
| $\mathcal N(0,1)$ and $\chi^2_\nu$ | $Z/\sqrt{V/\nu}$ | $t_\nu$ |
| $\chi^2_{d_1}$ and $\chi^2_{d_2}$ | ratio of normalised | $F_{d_1,d_2}$ |
| Uniform$(0,1)$ | $-\ln(U)/\lambda$ | Exponential$(\lambda)$ |

The last line is the **inverse-transform method** for generating random numbers (Ch. 24): if $F$ is a continuous CDF and $U$ is uniform, $F^{-1}(U)$ has CDF $F$.

> [!example] Worked Example 20.9 — Packet errors: binomial versus Poisson
> A 100-bit packet is sent over a channel with independent bit errors of probability $p=0.01$. The error-correcting code fixes up to 2 bit errors. Find the probability that a packet is lost.
>
> **Step 1 — model.** The number of bit errors is $X\sim$ Binomial$(100,0.01)$; a packet is lost if $X\ge3$.
>
> **Step 2 — exact.** $P(X=0)=0.99^{100}=0.36603$; $P(X=1)=100(0.01)(0.99)^{99}=0.36973$; $P(X=2)=\binom{100}{2}(0.01)^2(0.99)^{98}=0.18486$. So $P(X\ge3)=1-0.92063=0.0794$.
>
> **Step 3 — Poisson approximation.** $\lambda=np=1$: $P(X\ge3)\approx 1-e^{-1}(1+1+\tfrac12)=1-2.5e^{-1}=0.0803$.
>
> **Step 4 — comment.** The approximation is within $1.2\%$ (relative). Without the code, a packet would be lost whenever $X\ge1$, with probability $0.634$; the code reduces the loss rate by a factor of 8.

> [!example] Worked Example 20.10 — Photon counting
> A photodetector registers on average $\lambda=4$ photons per $1\ \mathrm{ms}$ gate. Find (a) the probability of an empty gate, (b) the probability of 8 or more counts, (c) the signal-to-noise ratio of a count.
>
> **(a)** $P(X=0)=e^{-4}=0.0183$.
>
> **(b)** $P(X\ge8)=1-\sum_{k=0}^{7}e^{-4}4^k/k! = 0.0511$. Doubling the mean count happens in 5% of gates purely by chance — the irreducible **shot noise** of light.
>
> **(c)** Mean $4$, standard deviation $\sqrt4=2$, so $\mathrm{SNR}=\lambda/\sqrt\lambda=\sqrt\lambda=2$. To improve SNR by a factor of 10 one must collect 100 times more photons (longer gate or larger aperture): shot-noise-limited SNR grows only as the square root of signal.

> [!example] Worked Example 20.11 — Reliability with exponential lifetimes
> A power module has constant failure rate $\lambda = 2\times10^{-5}\ \mathrm{h^{-1}}$ (MTBF $=50\,000\ \mathrm{h}$).
>
> **(a) One-year survival.** $P(T>8760\ \mathrm h)=e^{-8760/50000}=e^{-0.1752}=0.839$.
>
> **(b) Survival to the MTBF.** $P(T>50\,000)=e^{-1}=0.368$: only $37\%$ of units reach the MTBF — the MTBF is a mean, not a typical life.
>
> **(c) Memorylessness.** Given a module has already run $3$ years, the probability it survives one more year is still $0.839$. Burn-in therefore does nothing for exponential lifetimes; it helps only when early failures (decreasing hazard rate) are present.
>
> **(d) Series system.** Ten independent modules, all required: the system lifetime is the minimum of ten exponentials, which is exponential with rate $10\lambda$, so MTBF$_\text{sys}=5000\ \mathrm h$.

> [!example] Worked Example 20.12 — Gaussian yield and Rayleigh false alarms
> **(a) Yield.** Resistor values are $\mathcal N(1000\ \Omega,\ (4\ \Omega)^2)$. The fraction inside a $\pm1\%$ tolerance band ($990$–$1010\ \Omega$) is $P(\lvert Z\rvert\le 10/4)=1-2Q(2.5)=1-2(0.00621)=0.9876$, i.e. $1.24\%$ rejects.
>
> **(b) False alarms.** A radar detector compares the envelope of receiver noise with a threshold. The I and Q noise components are each $\mathcal N(0,\sigma^2)$ with $\sigma=1\ \mathrm{mV}$, so the envelope is Rayleigh$(1\ \mathrm{mV})$. For a threshold of $3\ \mathrm{mV}$, $P_\mathrm{FA}=e^{-9/2}=0.0111$. To achieve $P_\mathrm{FA}=10^{-6}$ the threshold must be $\sigma\sqrt{2\ln10^6}=5.26\ \mathrm{mV}$. Note how slowly the threshold grows: a millionfold reduction in false alarms costs less than a doubling of threshold.

> [!warning] Common Pitfall — Assuming Gaussian tails far from the mean
> Real component distributions are often truncated (parts outside tolerance were screened out) or have heavier tails than Gaussian. Predicting a $10^{-9}$ failure rate from a Gaussian fitted to 50 samples extrapolates far beyond the data; at $6\sigma$ the answer depends entirely on the assumed tail shape. Tail probabilities should be supported by data or by physics, not by the convenience of the Gaussian formula.

> [!tip] Engineering Note — Matching mechanism to model
> Ask what generates the data. Counting events in a fixed window: binomial (fixed $n$) or Poisson (no natural $n$). Time between events: exponential (no ageing) or Weibull (ageing). Sum of many small effects: Gaussian. Magnitude of a 2-D Gaussian: Rayleigh. Product of many small factors (component values after several process steps): lognormal, since the logarithm is a sum.

## 20.5 Sampling

### Statistics are random variables

A sample mean computed today from 16 readings will differ from one computed tomorrow from another 16. Before the data are taken, every statistic $T=g(X_1,\dots,X_n)$ is a random variable, and its distribution — the **sampling distribution** — tells us how far $T$ is likely to be from the quantity it is meant to estimate. All of inference rests on this idea: a confidence interval is a statement about the sampling distribution of an estimator, and a hypothesis test compares an observed statistic with its sampling distribution under a hypothesis.

> [!definition] Random sample
> $X_1,\dots,X_n$ is a **random sample** (i.i.d. sample) from a distribution $F$ if the $X_i$ are independent and each has CDF $F$. The **standard error** of a statistic is the standard deviation of its sampling distribution.

Independence is a physical assumption, not a mathematical gift. Readings taken faster than the correlation time of the noise (Ch. 21), units taken consecutively from one reel, or repeated measurements without re-seating a probe are *not* independent, and the formulas below then understate the true variability.

### Sampling distribution of the mean

> [!theorem] Mean and variance of the sample mean
> If $X_1,\dots,X_n$ is a random sample from a distribution with mean $\mu$ and variance $\sigma^2$, then
> $$
> E[\bar X]=\mu,\qquad \operatorname{Var}(\bar X)=\frac{\sigma^2}{n},\qquad \operatorname{SE}(\bar X)=\frac{\sigma}{\sqrt n}.
> $$
> If the sample is from $\mathcal N(\mu,\sigma^2)$, then $\bar X\sim\mathcal N(\mu,\sigma^2/n)$ exactly.

*Proof.* Linearity of expectation gives $E[\bar X]=\frac1n\sum E[X_i]=\mu$. Independence lets variances add: $\operatorname{Var}(\sum X_i)=n\sigma^2$, and dividing by $n^2$ gives $\sigma^2/n$. A linear combination of independent Gaussians is Gaussian. $\square$

The **$1/\sqrt n$ law** is the most important scaling in measurement: averaging four readings halves the noise, averaging a hundred reduces it tenfold. It is the statistical reason a lock-in amplifier, a signal averager on an oscilloscope, or a long integration time on a DMM improves resolution — and the reason improvement becomes expensive: each extra decimal digit costs a factor of 100 in measurement time.

If the readings are correlated with correlation coefficient $\rho_k$ between samples $k$ apart, the variance of the mean becomes

$$
\operatorname{Var}(\bar X)=\frac{\sigma^2}{n}\Bigl[1+2\sum_{k=1}^{n-1}\Bigl(1-\frac kn\Bigr)\rho_k\Bigr],
$$

so positive correlation inflates the standard error. The ratio $n_\text{eff}=n/[\,\cdot\,]$ is the **effective number of independent samples**. For first-order (exponentially correlated) noise with $\rho_k=a^k$ and large $n$, $n_\text{eff}\approx n(1-a)/(1+a)$.

**Finite populations.** When sampling *without replacement* $n$ items from a finite lot of $N$, the variance of the mean is multiplied by the finite-population correction $(N-n)/(N-1)$. It matters only when the sample is a sizeable fraction of the lot (more than about 5%).

### The central limit theorem

The Gaussian result above needs Gaussian data. The CLT says that for *any* distribution with finite variance, the mean of a large sample is approximately Gaussian anyway.

> [!theorem] Central limit theorem (Lindeberg–Lévy)
> Let $X_1,X_2,\dots$ be i.i.d. with mean $\mu$ and finite variance $\sigma^2>0$. Then the standardised mean
> $$
> Z_n=\frac{\bar X_n-\mu}{\sigma/\sqrt n}=\frac{\sum_{i=1}^n X_i-n\mu}{\sigma\sqrt n}
> $$
> converges in distribution to $\mathcal N(0,1)$: $P(Z_n\le z)\to\Phi(z)$ for every $z$.

*Proof sketch (via characteristic functions).* Let $Y_i=(X_i-\mu)/\sigma$, so $E[Y_i]=0$, $E[Y_i^2]=1$, and let $\varphi(\omega)=E[e^{j\omega Y}]$ be their common characteristic function (the Fourier transform of the pdf, Ch. 16). A Taylor expansion about $\omega=0$ gives

$$
\varphi(\omega)=1+j\omega E[Y]-\frac{\omega^2}{2}E[Y^2]+o(\omega^2)=1-\frac{\omega^2}{2}+o(\omega^2).
$$

Since $Z_n=\frac{1}{\sqrt n}\sum Y_i$ and characteristic functions of independent sums multiply,

$$
\varphi_{Z_n}(\omega)=\Bigl[\varphi\Bigl(\frac{\omega}{\sqrt n}\Bigr)\Bigr]^n=\Bigl[1-\frac{\omega^2}{2n}+o\Bigl(\frac1n\Bigr)\Bigr]^n\ \longrightarrow\ e^{-\omega^2/2},
$$

which is the characteristic function of $\mathcal N(0,1)$. Lévy's continuity theorem converts convergence of characteristic functions into convergence of CDFs. $\square$

In signal terms: the pdf of a sum is the *convolution* of the individual pdfs (Ch. 19), and repeated convolution of any reasonable pulse with itself tends to a Gaussian shape. The figure shows this for uniform pdfs, whose repeated convolutions (the Irwin–Hall densities) are piecewise polynomials.

```desmos-graph
left=-2; right=2; bottom=0; top=1.1
---
y=1|blue|-0.5<x<0.5
y=1-\operatorname{abs}\left(x\right)|red|-1<x<1
y=\frac{(x+1.5)^{2}}{2}|green|-1.5<x<-0.5
y=0.75-x^{2}|green|-0.5<x<0.5
y=\frac{(1.5-x)^{2}}{2}|green|0.5<x<1.5
y=\frac{1}{0.5\sqrt{2\pi}}e^{-\frac{x^{2}}{0.5}}|black|dashed
```

*Pdf of the centred sum of $n$ independent Uniform$(0,1)$ variables: $n=1$ (blue), $n=2$ (red, triangle), $n=3$ (green, piecewise quadratic), against the Gaussian with the $n=3$ variance $3/12=0.25$ (dashed). Three uniforms are already close to Gaussian.*

How large must $n$ be? The **Berry–Esseen theorem** bounds the error: $\sup_z\lvert P(Z_n\le z)-\Phi(z)\rvert\le C\rho/(\sigma^3\sqrt n)$ with $\rho=E\lvert X-\mu\rvert^3$ and $C<0.48$. The convergence is fast for symmetric, light-tailed data (uniform: $n\approx3$–$5$ suffices for the centre of the distribution), slow for skewed data (exponential, lognormal: $n$ of 30–50), and fails for distributions without finite variance (Cauchy). The approximation is always worst in the far tails — the CLT is a statement about the *bulk* of the distribution.

The CLT also explains the **normal approximation to the binomial and Poisson**: a Binomial$(n,p)$ count is a sum of $n$ Bernoullis, so it is approximately $\mathcal N(np,np(1-p))$ when $np(1-p)\gtrsim 10$; a Poisson$(\lambda)$ count is approximately $\mathcal N(\lambda,\lambda)$ for $\lambda\gtrsim10$. A **continuity correction** (evaluating at $k\pm\frac12$) improves the approximation for integer counts.

### Sampling distribution of the variance

For Gaussian samples the sample variance has an exact distribution, and it is independent of the sample mean.

> [!theorem] Sampling theorem for Gaussian samples
> If $X_1,\dots,X_n$ is a random sample from $\mathcal N(\mu,\sigma^2)$, then
> 1. $\bar X$ and $S^2$ are independent;
> 2. $\dfrac{(n-1)S^2}{\sigma^2}\sim\chi^2_{n-1}$;
> 3. $T=\dfrac{\bar X-\mu}{S/\sqrt n}\sim t_{n-1}$.

*Proof sketch.* Choose an orthogonal matrix $\mathbf Q$ (Ch. 9) whose first row is $(1,\dots,1)/\sqrt n$, and set $\mathbf Y=\mathbf Q(\mathbf X-\mu\mathbf 1)/\sigma$. Because $\mathbf Q$ is orthogonal and the $(X_i-\mu)/\sigma$ are i.i.d. $\mathcal N(0,1)$, the $Y_i$ are again i.i.d. $\mathcal N(0,1)$ (the covariance $\mathbf Q\mathbf I\mathbf Q^\mathsf T=\mathbf I$). Now $Y_1=\sqrt n(\bar X-\mu)/\sigma$, and since $\lVert\mathbf Y\rVert^2=\lVert(\mathbf X-\mu\mathbf1)/\sigma\rVert^2$,

$$
\sum_{i=2}^n Y_i^2=\frac{1}{\sigma^2}\Bigl[\sum_i(X_i-\mu)^2-n(\bar X-\mu)^2\Bigr]=\frac{(n-1)S^2}{\sigma^2}.
$$

So $\bar X$ depends only on $Y_1$ and $S^2$ only on $Y_2,\dots,Y_n$, giving independence (1), and $(n-1)S^2/\sigma^2$ is a sum of $n-1$ squared standard Gaussians, giving (2). Then $T=Y_1/\sqrt{\chi^2_{n-1}/(n-1)}$, which is $t_{n-1}$ by definition (3). $\square$

The rotation argument shows geometrically where the "lost" degree of freedom goes: one dimension of the $n$-dimensional data space (the direction of $\mathbf 1$) is used by the mean.

Consequences: $\operatorname{Var}(S^2)=2\sigma^4/(n-1)$, so the relative standard uncertainty of $S^2$ is $\sqrt{2/(n-1)}$ and that of $S$ is about $1/\sqrt{2(n-1)}$. With $n=10$ readings, $s$ is uncertain by about $24\%$; with $n=50$, by $10\%$. Standard deviations estimated from a handful of readings are themselves very uncertain — which is precisely why the $t$ distribution, rather than the Gaussian, is used for small samples.

> [!warning] Common Pitfall — Point 2 needs Gaussian data
> The $\chi^2$ result for $S^2$ is *not* protected by the CLT. For heavy-tailed data the variance of $S^2$ is $\sigma^4\bigl(\frac{2}{n-1}+\frac{\kappa}{n}\bigr)$ with $\kappa$ the excess kurtosis, which can be much larger than the Gaussian value. Intervals and tests for variances (Sections 20.7–20.8) are therefore far less robust than those for means.

### Two-sample statistics

For independent Gaussian samples of sizes $n_1,n_2$ with variances $\sigma_1^2,\sigma_2^2$:

$$
\bar X_1-\bar X_2\sim\mathcal N\Bigl(\mu_1-\mu_2,\ \frac{\sigma_1^2}{n_1}+\frac{\sigma_2^2}{n_2}\Bigr),\qquad \frac{S_1^2/\sigma_1^2}{S_2^2/\sigma_2^2}\sim F_{n_1-1,\,n_2-1}.
$$

The second result follows directly from the definition of $F$ as a ratio of independent normalised chi-squares.

### Sample proportions

If $X\sim$ Binomial$(n,p)$ counts defectives in a sample, the sample proportion $\hat p=X/n$ has $E[\hat p]=p$ and $\operatorname{SE}(\hat p)=\sqrt{p(1-p)/n}$, and is approximately Gaussian for $np(1-p)\gtrsim10$.

> [!example] Worked Example 20.13 — How much averaging is enough?
> An ADC channel has Gaussian reading noise with $\sigma=5\ \mathrm{mV}$.
>
> **(a)** For a 16-reading average, $\operatorname{SE}=5/\sqrt{16}=1.25\ \mathrm{mV}$, and $P(\lvert\bar X-\mu\rvert>2\ \mathrm{mV})=2Q(2/1.25)=2Q(1.6)=0.110$.
>
> **(b)** To make the standard error at most $0.5\ \mathrm{mV}$ we need $\sigma/\sqrt n\le0.5$, i.e. $n\ge(5/0.5)^2=100$ readings.
>
> **(c)** If the readings are correlated with $\rho_k=0.6^k$ (the noise is lowpass-filtered before sampling), the effective sample size is only $n_\text{eff}\approx n(0.4/1.6)=n/4$. One then needs about 400 readings — or, better, should sample more slowly so that successive readings become independent.

> [!example] Worked Example 20.14 — CLT for skewed lifetimes
> Component lifetimes are exponential with mean $1000\ \mathrm h$ (so $\sigma=1000\ \mathrm h$). A life test runs 40 parts. Find the probability that the sample mean lifetime is below $800\ \mathrm h$.
>
> **Step 1 — CLT.** $\operatorname{SE}=1000/\sqrt{40}=158.1\ \mathrm h$. $P(\bar X<800)\approx\Phi\bigl((800-1000)/158.1\bigr)=\Phi(-1.265)=0.103$.
>
> **Step 2 — exact.** A sum of 40 independent exponentials is gamma distributed (shape 40, scale $1000$), so $\bar X$ is gamma with shape 40 and scale 25; its CDF at 800 gives $0.0956$.
>
> **Step 3 — comment.** The CLT answer is about $8\%$ high in relative terms because the exponential is strongly skewed; the true distribution of $\bar X$ is still slightly right-skewed at $n=40$. For a quick engineering estimate the CLT is adequate; for an acceptance test with contractual consequences, use the exact distribution.

> [!example] Worked Example 20.15 — Distribution of a noise-variance estimate
> A noise generator is specified as $\sigma^2=4\ \mathrm{mV^2}$. From $n=10$ Gaussian samples, what is the probability that the estimate $s^2$ exceeds $6\ \mathrm{mV^2}$, even though the generator meets specification?
>
> **Step 1 — pivot.** $(n-1)S^2/\sigma^2=9S^2/4\sim\chi^2_9$.
>
> **Step 2 — translate.** $S^2>6$ is equivalent to $\chi^2_9>9(6)/4=13.5$.
>
> **Step 3 — evaluate.** From chi-square tables or software, $P(\chi^2_9>13.5)=0.141$.
>
> **Comment.** In one test out of seven, a perfectly good generator would appear to be $50\%$ noisier than specified. The central $95\%$ range of $\chi^2_9$ is $[2.70,\ 19.02]$, so $s^2$ ranges over $[1.20,\ 8.45]\ \mathrm{mV^2}$ — a factor of seven. Ten samples are far too few to verify a variance.

> [!example] Worked Example 20.16 — Sampling a production lot
> A process has defect rate $p=2\%$. In a sample of $n=500$ parts, what is the probability that the observed defect rate is $3\%$ or more (15 or more defectives)?
>
> **Step 1 — mean and SE.** $E[X]=np=10$, $\operatorname{Var}X=np(1-p)=9.8$, $\operatorname{SD}=3.13$.
>
> **Step 2 — normal approximation with continuity correction.** $P(X\ge15)\approx Q\bigl((14.5-10)/3.13\bigr)=Q(1.44)=0.075$.
>
> **Step 3 — exact.** The binomial tail $P(X\ge15)=0.081$. The approximation is serviceable; the remaining discrepancy reflects the right skew of a binomial with small $p$. (Without the continuity correction one would get $Q(1.60)=0.055$, noticeably worse.)

> [!tip] Engineering Note — Dithering and the CLT
> Adding a small random dither before quantisation and then averaging many readings lets an ADC resolve DC levels far finer than one LSB: the dither makes the quantisation errors of successive readings independent and zero-mean, so their average converges at the $1/\sqrt n$ rate. Without dither, a constant input always produces the same code and averaging gains nothing.

## 20.6 Estimation

### The estimation task

A model has unknown parameters $\theta$ — the mean and standard deviation of a noise source, the failure rate of a component, the Rayleigh parameter of a fading channel, the slope of a sensor. An **estimator** $\hat\theta=g(X_1,\dots,X_n)$ is a rule that turns data into a value for $\theta$; the number it produces for a particular data set is an **estimate**. Because $\hat\theta$ is a statistic, it has a sampling distribution, and "good" estimators are those whose sampling distribution is concentrated near the true $\theta$.

### Properties of estimators

> [!definition] Bias, variance, MSE, consistency, efficiency
> - **Bias:** $b(\hat\theta)=E[\hat\theta]-\theta$. $\hat\theta$ is **unbiased** if $b=0$ for all $\theta$.
> - **Mean squared error:** $\operatorname{MSE}(\hat\theta)=E\bigl[(\hat\theta-\theta)^2\bigr]$.
> - **Consistency:** $\hat\theta_n\to\theta$ in probability as $n\to\infty$.
> - **Efficiency:** among unbiased estimators, one with smaller variance is more efficient; the relative efficiency of $\hat\theta_1$ to $\hat\theta_2$ is $\operatorname{Var}\hat\theta_2/\operatorname{Var}\hat\theta_1$.

> [!theorem] Bias–variance decomposition
> $$
> \operatorname{MSE}(\hat\theta)=\operatorname{Var}(\hat\theta)+b(\hat\theta)^2 .
> $$

*Proof.* Add and subtract $E\hat\theta$: $E[(\hat\theta-E\hat\theta+E\hat\theta-\theta)^2]=E[(\hat\theta-E\hat\theta)^2]+2(E\hat\theta-\theta)E[\hat\theta-E\hat\theta]+(E\hat\theta-\theta)^2$, and the middle expectation is zero. $\square$

Unbiasedness is not sacred. A slightly biased estimator with much smaller variance can have lower MSE — for Gaussian data, the variance estimator with divisor $n+1$ has the smallest MSE of all estimators of the form $c\sum(X_i-\bar X)^2$, even though it is biased. A sufficient condition for consistency is that both bias and variance tend to zero; $\bar X$ is consistent for $\mu$ (by the law of large numbers, Ch. 19) and $S^2$ is consistent for $\sigma^2$.

### The Cramér–Rao lower bound

How good can an unbiased estimator possibly be? For a model with pdf $f(x;\theta)$ that is smooth in $\theta$, define the **score** $\partial_\theta\ln f(X;\theta)$ and the **Fisher information** per observation

$$
I(\theta)=E\Bigl[\Bigl(\frac{\partial}{\partial\theta}\ln f(X;\theta)\Bigr)^2\Bigr]=-E\Bigl[\frac{\partial^2}{\partial\theta^2}\ln f(X;\theta)\Bigr].
$$

> [!theorem] Cramér–Rao lower bound (CRLB)
> Under regularity conditions (support independent of $\theta$, differentiation under the integral allowed), every unbiased estimator from an i.i.d. sample of size $n$ satisfies
> $$
> \operatorname{Var}(\hat\theta)\ge\frac{1}{n\,I(\theta)} .
> $$
> An unbiased estimator attaining the bound is called **efficient**.

*Proof sketch.* Let $L(\theta)=\prod f(x_i;\theta)$ be the joint pdf and $U=\partial_\theta\ln L$ the total score. Differentiating $\int L\,d\mathbf x=1$ gives $E[U]=0$; differentiating $\int\hat\theta L\,d\mathbf x=\theta$ gives $E[\hat\theta U]=1$, so $\operatorname{Cov}(\hat\theta,U)=1$. The Cauchy–Schwarz inequality for covariances gives $1\le\operatorname{Var}(\hat\theta)\operatorname{Var}(U)$, and $\operatorname{Var}(U)=nI(\theta)$ by independence. $\square$

For a DC level $\theta$ in white Gaussian noise of variance $\sigma^2$, $\ln f=-\frac{(x-\theta)^2}{2\sigma^2}+\text{const}$, so $I=1/\sigma^2$ and the bound is $\sigma^2/n$ — exactly the variance of $\bar X$. The sample mean is efficient: no unbiased processing of the samples can do better. The CRLB is the standard benchmark in estimation theory for radar, communications and sensor arrays.

### Method of moments

The oldest general recipe: equate sample moments to population moments and solve.

> [!definition] Method of moments (MoM)
> If the model has $p$ parameters and its first $p$ moments are known functions $\mu_k'(\theta)=E[X^k]$, the MoM estimator solves
> $$
> \frac1n\sum_{i=1}^n x_i^k=\mu_k'(\hat\theta),\qquad k=1,\dots,p .
> $$
> (Central moments may be used instead of raw moments.)

MoM is simple, needs only moment formulas, and gives consistent estimators (sample moments converge to population moments, and continuous functions preserve convergence). It is often a good *starting value* for an iterative maximum-likelihood search. Its weakness is that it can be inefficient and can produce impossible values (outside the parameter space or the data range).

### Maximum likelihood

> [!definition] Likelihood and maximum likelihood estimator (MLE)
> For observed data $x_1,\dots,x_n$, the **likelihood** is the joint pdf (or pmf) viewed as a function of the parameter:
> $$
> L(\theta)=\prod_{i=1}^n f(x_i;\theta),\qquad \ell(\theta)=\ln L(\theta)=\sum_{i=1}^n\ln f(x_i;\theta).
> $$
> The **MLE** $\hat\theta$ maximises $L$ (equivalently $\ell$). When the maximum is interior and $\ell$ is smooth, it solves the **likelihood equation** $\partial\ell/\partial\theta=0$.

The MLE chooses the parameter value under which the observed data were most probable. It is the dominant estimation method in engineering because of its large-sample properties:

> [!theorem] Asymptotic properties of the MLE
> Under regularity conditions, as $n\to\infty$,
> 1. $\hat\theta$ is **consistent**;
> 2. $\hat\theta$ is **asymptotically Gaussian and efficient**: $\hat\theta\approx\mathcal N\bigl(\theta,\,1/(nI(\theta))\bigr)$;
> 3. **invariance:** the MLE of $g(\theta)$ is $g(\hat\theta)$.

Point 2 provides standard errors almost for free: $\operatorname{SE}(\hat\theta)\approx 1/\sqrt{-\ell''(\hat\theta)}$, the inverse square root of the curvature of the log-likelihood at its peak (the **observed information**). A sharply peaked log-likelihood means the data pin $\theta$ down well. With several parameters, the inverse of the matrix of second derivatives (the negative Hessian of $\ell$) approximates the covariance matrix of $\hat{\boldsymbol\theta}$. Point 3 means, for example, that the MLE of a standard deviation is the square root of the MLE of the variance, and the MLE of an MTBF is the reciprocal of the MLE of the failure rate.

> [!theorem] MLE in Gaussian noise is least squares
> If $x_i=s_i(\boldsymbol\theta)+w_i$ with $w_i$ i.i.d. $\mathcal N(0,\sigma^2)$, then
> $$
> \ell(\boldsymbol\theta)=-\frac n2\ln(2\pi\sigma^2)-\frac{1}{2\sigma^2}\sum_{i=1}^n\bigl(x_i-s_i(\boldsymbol\theta)\bigr)^2,
> $$
> so maximising the likelihood over $\boldsymbol\theta$ is the same as minimising the sum of squared residuals.

This single observation justifies least-squares curve fitting (Section 20.9, Ch. 23) whenever the noise is Gaussian, and suggests weighted least squares (weights $1/\sigma_i^2$) when the noise variances differ.

### Standard MLEs

| Model | MLE | Notes |
|---|---|---|
| Bernoulli$(p)$ | $\hat p=\bar x$ (proportion of successes) | unbiased |
| Poisson$(\lambda)$ | $\hat\lambda=\bar x$ | unbiased, efficient |
| Exponential$(\lambda)$ | $\hat\lambda=1/\bar x$ | biased: $E\hat\lambda=n\lambda/(n-1)$ |
| $\mathcal N(\mu,\sigma^2)$ | $\hat\mu=\bar x$, $\hat\sigma^2=\frac1n\sum(x_i-\bar x)^2$ | $\hat\sigma^2$ biased by factor $(n-1)/n$ |
| Uniform$(0,\theta)$ | $\hat\theta=x_{(n)}$ | biased low; support depends on $\theta$ |
| Rayleigh$(\sigma)$ | $\hat\sigma^2=\frac{1}{2n}\sum r_i^2$ | unbiased for $\sigma^2$ |

*Derivation for the Gaussian.* $\ell(\mu,\sigma^2)=-\frac n2\ln(2\pi\sigma^2)-\frac{1}{2\sigma^2}\sum(x_i-\mu)^2$. Setting $\partial\ell/\partial\mu=\frac{1}{\sigma^2}\sum(x_i-\mu)=0$ gives $\hat\mu=\bar x$; setting $\partial\ell/\partial\sigma^2=-\frac{n}{2\sigma^2}+\frac{1}{2\sigma^4}\sum(x_i-\hat\mu)^2=0$ gives $\hat\sigma^2=\frac1n\sum(x_i-\bar x)^2$. The MLE of the variance is biased because it does not account for the degree of freedom used by $\hat\mu$ — the MLE optimises fit, not unbiasedness. The bias vanishes as $n\to\infty$, consistent with the asymptotic theory.

```desmos-graph
left=0; right=0.012; bottom=-6; top=0.5
---
y=8\ln\left(\frac{x}{0.00396}\right)-2020\left(x-0.00396\right)|blue|0.0005<x<0.012
y=-0.5|red|dashed
x=0.00396|grey|dashed
```

*Log-likelihood of the exponential failure rate in Example 20.18, shifted so that its maximum is 0. The MLE sits at the peak; the values of $\lambda$ where $\ell$ has fallen by $\tfrac12$ (red line) give an approximate $68\%$ likelihood interval, which is visibly asymmetric.*

> [!example] Worked Example 20.17 — Uniform: method of moments versus maximum likelihood
> Timing jitter is modelled as Uniform$(0,\theta)$ (in ns). Five observations: $0.8,\ 2.9,\ 1.7,\ 3.6,\ 2.2$. Estimate $\theta$ by both methods and compare their quality.
>
> **Step 1 — MoM.** $E[X]=\theta/2$, so $\hat\theta_\text{MoM}=2\bar x=2(2.24)=4.48\ \mathrm{ns}$. It is unbiased, with $\operatorname{Var}=4\operatorname{Var}(X)/n=4\theta^2/(12n)=\theta^2/(3n)$.
>
> **Step 2 — MLE.** $L(\theta)=\theta^{-n}$ if $\theta\ge x_{(n)}$, and $0$ otherwise (some observation would be impossible). $L$ decreases in $\theta$, so it is maximised at the smallest admissible value: $\hat\theta_\text{ML}=x_{(n)}=3.6\ \mathrm{ns}$. The likelihood equation $\ell'=0$ has no solution — the maximum is at a boundary, and the regularity conditions of the CRLB fail.
>
> **Step 3 — bias correction.** The maximum of $n$ uniforms has $E[X_{(n)}]=n\theta/(n+1)$, so $\hat\theta_\text{U}=\frac{n+1}{n}x_{(n)}=4.32\ \mathrm{ns}$ is unbiased, with variance $\theta^2/(n(n+2))$.
>
> **Step 4 — compare.** For $n=5$: $\operatorname{Var}\hat\theta_\text{MoM}=\theta^2/15$, $\operatorname{Var}\hat\theta_\text{U}=\theta^2/35$. The order-statistic estimator is more than twice as efficient, and its advantage grows linearly with $n$ (variance $\propto1/n^2$ instead of $1/n$). Note also that MoM could produce an estimate *smaller* than the largest observation — a logically impossible value — for unlucky samples.

> [!example] Worked Example 20.18 — Failure rate from a life test
> Eight units run to failure at $120,\ 340,\ 75,\ 510,\ 260,\ 195,\ 430,\ 90\ \mathrm h$. Assume exponential lifetimes and find the MLE of $\lambda$, its standard error and the MLE of the MTBF.
>
> **Step 1 — log-likelihood.** $\ell(\lambda)=\sum\bigl(\ln\lambda-\lambda t_i\bigr)=n\ln\lambda-\lambda\sum t_i = 8\ln\lambda-2020\lambda$.
>
> **Step 2 — maximise.** $\ell'(\lambda)=8/\lambda-2020=0$ gives $\hat\lambda=8/2020=3.96\times10^{-3}\ \mathrm{h^{-1}}$. $\ell''=-8/\lambda^2<0$ confirms a maximum.
>
> **Step 3 — standard error.** Observed information $-\ell''(\hat\lambda)=n/\hat\lambda^2$, so $\operatorname{SE}(\hat\lambda)\approx\hat\lambda/\sqrt n=1.40\times10^{-3}\ \mathrm{h^{-1}}$ — a $35\%$ relative uncertainty from eight failures.
>
> **Step 4 — invariance.** $\widehat{\mathrm{MTBF}}=1/\hat\lambda=\bar t=252.5\ \mathrm h$.
>
> **Comment.** The relative uncertainty $1/\sqrt n$ depends only on the number of *failures*. Reliability demonstrations are expensive precisely because small failure rates produce few failures; Section 20.7 gives an exact interval using the chi-square distribution.

> [!example] Worked Example 20.19 — Rayleigh parameter of a noise envelope
> Envelope samples of narrowband receiver noise (mV): $1.2,\ 0.7,\ 2.1,\ 1.5,\ 0.9,\ 1.8$. Estimate $\sigma$ by ML and by the method of moments.
>
> **Step 1 — log-likelihood.** $\ell(\sigma^2)=\sum\ln r_i-n\ln\sigma^2-\frac{1}{2\sigma^2}\sum r_i^2$.
>
> **Step 2 — maximise.** $\partial\ell/\partial\sigma^2=-\frac{n}{\sigma^2}+\frac{1}{2\sigma^4}\sum r_i^2=0\Rightarrow\hat\sigma^2=\frac{1}{2n}\sum r_i^2=\frac{12.64}{12}=1.053\ \mathrm{mV^2}$, so $\hat\sigma_\text{ML}=1.026\ \mathrm{mV}$ (invariance).
>
> **Step 3 — method of moments.** $E[R]=\sigma\sqrt{\pi/2}$, so $\hat\sigma_\text{MoM}=\bar r/1.2533=1.367/1.2533=1.090\ \mathrm{mV}$.
>
> **Step 4 — interpretation.** The MLE uses $\sum r_i^2$, which is the estimated *noise power* ($E[R^2]=2\sigma^2$): it is the natural statistic, and $2n\hat\sigma^2/\sigma^2\sim\chi^2_{2n}$ exactly, which gives exact intervals. The two estimates differ by $6\%$ — well within the sampling scatter for $n=6$.

> [!example] Worked Example 20.20 — Cramér–Rao bound for a DC level
> A DC voltage $A$ is measured in white Gaussian noise of $\sigma=2\ \mathrm{mV}$ with $n=25$ samples.
>
> **Step 1 — Fisher information.** $\ln f(x;A)=-\frac{(x-A)^2}{2\sigma^2}+\text{const}$, $\partial^2_A\ln f=-1/\sigma^2$, so $I(A)=1/\sigma^2$.
>
> **Step 2 — bound.** $\operatorname{Var}(\hat A)\ge\sigma^2/n=4/25=0.16\ \mathrm{mV^2}$, i.e. $\operatorname{SD}\ge0.4\ \mathrm{mV}$.
>
> **Step 3 — efficiency of the mean.** $\operatorname{Var}(\bar X)=\sigma^2/n$ equals the bound, so the sample mean is efficient. The median, with asymptotic variance $\pi\sigma^2/(2n)=0.251\ \mathrm{mV^2}$, reaches only $64\%$ efficiency.
>
> **Step 4 — non-Gaussian noise.** For Laplacian (double-exponential) noise with the same variance, the Fisher information doubles to $2/\sigma^2$ and the MLE becomes the *median*, which is then efficient while the mean is only $50\%$ efficient. The best estimator depends on the noise distribution — a key lesson of estimation theory.

> [!warning] Common Pitfall — Maximum likelihood with few data
> The attractive properties of the MLE are *asymptotic*. With small $n$ it can be noticeably biased ($\hat\sigma^2$ with divisor $n$; $\hat\lambda=1/\bar t$), and the Gaussian approximation to its distribution can be poor (the log-likelihood in Example 20.18 is visibly asymmetric). Quote small-sample results with exact intervals where available, or check by simulation (Ch. 24).

> [!tip] Engineering Note — Estimators inside equipment
> A spectrum analyser's "noise marker" estimates power by averaging squared magnitudes (the Rayleigh MLE in disguise); its "video averaging" in log mode instead averages decibels, which is biased low by $2.51\ \mathrm{dB}$ for Gaussian noise — a correction built into good instruments. A GPS receiver computes position by iterative least squares, i.e. maximum likelihood under Gaussian range errors, and reports a dilution-of-precision that is the CRLB-like covariance of the fix.

## 20.7 Confidence intervals

### From a point to an interval

A point estimate without an indication of its reliability is of limited use: "the noise floor is $-112.3\ \mathrm{dBm}$" means one thing if the estimate is good to $0.1\ \mathrm{dB}$ and quite another if it is good to $3\ \mathrm{dB}$. A **confidence interval** reports a range of plausible parameter values together with a stated probability that the *procedure* captures the true value.

> [!definition] Confidence interval
> An interval $[L,U]$ whose endpoints are statistics $L=g_1(X_1,\dots,X_n)$, $U=g_2(X_1,\dots,X_n)$ is a $100(1-\alpha)\%$ **confidence interval** for $\theta$ if
> $$
> P(L\le\theta\le U)=1-\alpha\quad\text{for every value of }\theta .
> $$
> $1-\alpha$ is the **confidence level** (coverage probability). One-sided bounds $(-\infty,U]$ or $[L,\infty)$ are defined similarly.

The probability statement concerns the random endpoints, not $\theta$, which is fixed. Once the data are in and the interval is, say, $[2.475,\ 2.530]\ \mathrm V$, the true value either is or is not in it; the $95\%$ describes the long-run success rate of the method. If many laboratories each compute a $95\%$ interval from their own data, about $95\%$ of the intervals will contain the true value.

### The pivotal method

Most intervals in practice are built from a **pivot**: a function $Q(X_1,\dots,X_n;\theta)$ whose distribution is known and does not depend on $\theta$. If $q_{\alpha/2}$ and $q_{1-\alpha/2}$ are quantiles of that distribution, then $P(q_{\alpha/2}\le Q\le q_{1-\alpha/2})=1-\alpha$, and solving the inequalities for $\theta$ yields the interval. The sampling theorems of Section 20.5 supply the pivots:

| Parameter | Conditions | Pivot | Distribution |
|---|---|---|---|
| $\mu$ | Gaussian (or large $n$), $\sigma$ known | $(\bar X-\mu)/(\sigma/\sqrt n)$ | $\mathcal N(0,1)$ |
| $\mu$ | Gaussian, $\sigma$ unknown | $(\bar X-\mu)/(S/\sqrt n)$ | $t_{n-1}$ |
| $\sigma^2$ | Gaussian | $(n-1)S^2/\sigma^2$ | $\chi^2_{n-1}$ |
| $\sigma_1^2/\sigma_2^2$ | two Gaussian samples | $(S_1^2/\sigma_1^2)/(S_2^2/\sigma_2^2)$ | $F_{n_1-1,n_2-1}$ |
| $\lambda$ (exponential) | complete sample | $2\lambda\sum T_i$ | $\chi^2_{2n}$ |

Below, $z_{\alpha}$, $t_{\alpha,\nu}$ and $\chi^2_{\alpha,\nu}$ denote **upper** $\alpha$-points: $P(Z>z_\alpha)=\alpha$, and so on.

### Interval for a mean

> [!theorem] Confidence intervals for a mean
> **$\sigma$ known** (or $n$ large, with $s$ in place of $\sigma$):
> $$
> \bar x\pm z_{\alpha/2}\frac{\sigma}{\sqrt n}.
> $$
> **$\sigma$ unknown, Gaussian data:**
> $$
> \bar x\pm t_{\alpha/2,\,n-1}\frac{s}{\sqrt n}.
> $$

*Derivation of the $t$ interval.* $P\bigl(-t_{\alpha/2,n-1}\le\frac{\bar X-\mu}{S/\sqrt n}\le t_{\alpha/2,n-1}\bigr)=1-\alpha$; multiplying through by $S/\sqrt n$ and rearranging for $\mu$ gives $P\bigl(\bar X-t\,S/\sqrt n\le\mu\le\bar X+t\,S/\sqrt n\bigr)=1-\alpha$. $\square$

The half-width $E=z_{\alpha/2}\sigma/\sqrt n$ is the **margin of error**. Inverting it gives the sample size needed for a target margin:

$$
n=\Bigl(\frac{z_{\alpha/2}\,\sigma}{E}\Bigr)^2 .
$$

Because of the $1/\sqrt n$ law, halving the interval requires four times the data. The $t$ interval is robust to moderate non-normality (thanks to the CLT) once $n$ exceeds about 15–20, unless the data are strongly skewed or contain outliers.

### Interval for a variance

> [!theorem] Confidence interval for $\sigma^2$ (Gaussian data)
> $$
> \frac{(n-1)s^2}{\chi^2_{\alpha/2,\,n-1}}\ \le\ \sigma^2\ \le\ \frac{(n-1)s^2}{\chi^2_{1-\alpha/2,\,n-1}} .
> $$
> Taking square roots gives the interval for $\sigma$.

The interval is not symmetric about $s^2$, because the chi-square distribution is skewed; the upper limit lies further from $s^2$ than the lower. Unlike the $t$ interval, it is *not* robust to non-Gaussian data (Section 20.5).

### Interval for a proportion

For $x$ successes in $n$ trials, $\hat p=x/n$. The textbook **Wald interval** $\hat p\pm z_{\alpha/2}\sqrt{\hat p(1-\hat p)/n}$ is simple but performs poorly when $np$ is small: its actual coverage can fall well below the nominal level, and it can extend below zero. The **Wilson score interval** inverts the test statistic $(\hat p-p)/\sqrt{p(1-p)/n}$ exactly, solving a quadratic in $p$:

$$
\frac{\hat p+\dfrac{z^2}{2n}\ \pm\ z\sqrt{\dfrac{\hat p(1-\hat p)}{n}+\dfrac{z^2}{4n^2}}}{1+\dfrac{z^2}{n}},\qquad z=z_{\alpha/2}.
$$

It is recommended for routine use. The exact (Clopper–Pearson) interval, computed from binomial tail probabilities, is conservative (coverage at least $1-\alpha$) and is preferred for contractual or safety-related claims.

### Interval for a difference of means

For two independent samples, the **Welch interval** is

$$
(\bar x_1-\bar x_2)\pm t_{\alpha/2,\nu}\sqrt{\frac{s_1^2}{n_1}+\frac{s_2^2}{n_2}},\qquad
\nu=\frac{\bigl(s_1^2/n_1+s_2^2/n_2\bigr)^2}{\dfrac{(s_1^2/n_1)^2}{n_1-1}+\dfrac{(s_2^2/n_2)^2}{n_2-1}} ,
$$

where $\nu$ is the Welch–Satterthwaite effective degrees of freedom — the same formula reappears in the GUM (Section 20.10). If the two variances are known to be equal, the **pooled** variance $s_p^2=\frac{(n_1-1)s_1^2+(n_2-1)s_2^2}{n_1+n_2-2}$ with $n_1+n_2-2$ degrees of freedom may be used instead. For paired data (the same units measured twice), compute the differences $d_i$ and use the one-sample $t$ interval on them — pairing removes unit-to-unit variation and usually gives a much narrower interval.

### Confidence, prediction and tolerance intervals

Three different questions are often confused:

| Interval | Question | Gaussian form | As $n\to\infty$ |
|---|---|---|---|
| Confidence for $\mu$ | where is the population mean? | $\bar x\pm t\,s/\sqrt n$ | shrinks to a point |
| Prediction | where will the *next single* observation fall? | $\bar x\pm t\,s\sqrt{1+1/n}$ | tends to $\mu\pm z\sigma$ |
| Tolerance | where do $P\%$ of all units lie (with confidence $\gamma$)? | $\bar x\pm k\,s$ (tabulated $k$) | tends to $\mu\pm z\sigma$ |

The prediction interval follows because $X_\text{new}-\bar X$ has variance $\sigma^2(1+1/n)$. A datasheet "typical" range or a production limit is a tolerance-interval question, not a confidence-interval question.

> [!example] Worked Example 20.21 — Interval for a reference voltage
> Use the eight readings of Example 20.3: $\bar x=2.5025\ \mathrm V$, $s=32.40\ \mathrm{mV}$, $n=8$.
>
> **(a) $\sigma$ unknown.** $\nu=7$, $t_{0.025,7}=2.365$. Half-width $=2.365\times32.40/\sqrt8=27.1\ \mathrm{mV}$. The $95\%$ CI is $[2.4754,\ 2.5296]\ \mathrm V$.
>
> **(b) If $\sigma=32.4\ \mathrm{mV}$ were known** from long experience with this ADC, the half-width would be $1.960\times32.4/\sqrt8=22.5\ \mathrm{mV}$. The $t$ interval is $21\%$ wider: the price of estimating $\sigma$ from only eight readings.
>
> **(c) Sample size.** To obtain a $95\%$ half-width of $5\ \mathrm{mV}$ (with $\sigma\approx32.4\ \mathrm{mV}$): $n=(1.96\times32.4/5)^2=161.3$, so $n=162$ readings.
>
> **(d) Prediction.** A single future reading will lie in $2.5025\pm2.365\times32.40\sqrt{1+1/8}=2.5025\pm0.0813\ \mathrm V$ with $95\%$ probability — three times wider than the CI for the mean.

> [!example] Worked Example 20.22 — Interval for the noise standard deviation
> From the same data, find a $95\%$ CI for $\sigma$.
>
> **Step 1 — quantiles.** For $\chi^2_7$: $\chi^2_{0.975,7}=1.690$, $\chi^2_{0.025,7}=16.013$.
>
> **Step 2 — interval for $\sigma^2$.** $(n-1)s^2=7\times1.05\times10^{-3}=7.35\times10^{-3}\ \mathrm{V^2}$. The interval is $[7.35\times10^{-3}/16.013,\ 7.35\times10^{-3}/1.690]=[4.59\times10^{-4},\ 4.35\times10^{-3}]\ \mathrm{V^2}$.
>
> **Step 3 — interval for $\sigma$.** $[21.4,\ 66.0]\ \mathrm{mV}$, around the point estimate $32.4\ \mathrm{mV}$. The upper limit is double the estimate: a standard deviation from eight readings is known only to within a factor of about two.

> [!example] Worked Example 20.23 — Defect rate: Wald, Wilson and exact
> A sample of $n=200$ boards contains $x=6$ with solder defects. Find a $95\%$ CI for the defect rate $p$.
>
> **Step 1 — estimate.** $\hat p=0.030$, $\sqrt{\hat p(1-\hat p)/n}=0.01206$.
>
> **Step 2 — Wald.** $0.030\pm1.96(0.01206)=[0.0064,\ 0.0536]$.
>
> **Step 3 — Wilson.** With $z^2/n=0.01921$: centre $=(0.030+0.00960)/1.01921=0.0389$; half-width $=1.96\sqrt{1.455\times10^{-4}+2.40\times10^{-5}}/1.01921=0.0250$. Interval $[0.0138,\ 0.0639]$.
>
> **Step 4 — exact (Clopper–Pearson).** $[0.0111,\ 0.0642]$.
>
> **Comment.** With only six defects, the Wald interval is shifted too far to the left and is too short on the right; Wilson is close to the exact interval. Note the asymmetry: the data are compatible with a defect rate more than twice the point estimate.

> [!example] Worked Example 20.24 — Exact interval for a failure rate, and zero-failure tests
> **(a)** In Example 20.18, $n=8$ failures in total time $\sum t_i=2020\ \mathrm h$. The pivot $2\lambda\sum T_i\sim\chi^2_{16}$, with $\chi^2_{0.975,16}=6.908$ and $\chi^2_{0.025,16}=28.845$. Hence
> $$
> \frac{6.908}{2\times2020}\le\lambda\le\frac{28.845}{2\times2020}\quad\Longrightarrow\quad 1.71\times10^{-3}\le\lambda\le7.14\times10^{-3}\ \mathrm{h^{-1}},
> $$
> and the MTBF lies in $[140,\ 585]\ \mathrm h$ with $95\%$ confidence. The interval is strongly asymmetric about $\hat\lambda=3.96\times10^{-3}$ — compare the Gaussian approximation $\hat\lambda\pm1.96\,\mathrm{SE}=[1.2,\ 6.7]\times10^{-3}$.
>
> **(b) Zero failures.** A supplier claims $\lambda\le10^{-4}\ \mathrm{h^{-1}}$. How many device-hours $T$ of testing with no failures demonstrate this at $90\%$ confidence? For a time-terminated test with zero failures, the one-sided upper bound is $\lambda_U=\chi^2_{0.10,2}/(2T)=4.605/(2T)=2.303/T$. Setting $\lambda_U=10^{-4}$ gives $T=23\,026$ device-hours — e.g. 100 units for 230 h. (The same "rule of three" in another guise: with zero failures, the $95\%$ upper bound is $3/T$.)

> [!warning] Common Pitfall — "There is a 95% probability that $\mu$ is in this interval"
> In the frequentist framework used here, $\mu$ is not random, so this sentence is not strictly meaningful; the $95\%$ is a property of the procedure. A related error is to judge two results "different" only if their $95\%$ intervals do not overlap: two intervals can overlap substantially while the difference between the means is still significant. Compare means with an interval or test for the *difference* instead.

> [!tip] Engineering Note — Coverage factors in measurement reports
> Calibration certificates quote an expanded uncertainty $U=ku_c$ with $k=2$ "corresponding to a coverage probability of approximately $95\%$". This is the GUM version of a confidence interval: $k=2$ replaces $z_{0.025}=1.96$, and when the effective degrees of freedom are small, $k$ is taken from the $t$ distribution instead (Section 20.10).

## 20.8 Hypothesis testing

### Decisions under uncertainty

Many engineering questions have yes/no answers: Does this batch meet the $1\%$ tolerance? Did the firmware update change the mean current draw? Is the new supplier's capacitor noisier? Is this random number generator uniform? Is there a target in this radar cell? A **hypothesis test** is a rule for answering such questions from data, with controlled probabilities of the two possible mistakes.

> [!definition] Hypotheses, test statistic, rejection region
> - The **null hypothesis** $H_0$ is the default claim (no change, no effect, meets specification); the **alternative** $H_1$ is what we look for evidence of.
> - A **test statistic** $T$ is a statistic whose distribution under $H_0$ is known.
> - The **rejection (critical) region** $\mathcal R$ is the set of values of $T$ for which $H_0$ is rejected.
> - A hypothesis such as $\mu=\mu_0$ that fixes the distribution completely is **simple**; one such as $\mu>\mu_0$ is **composite**.

| | $H_0$ true | $H_0$ false |
|---|---|---|
| Reject $H_0$ | **Type I error** (false alarm), probability $\alpha$ | correct (detection), probability $1-\beta$ = **power** |
| Do not reject $H_0$ | correct | **Type II error** (miss), probability $\beta$ |

The **significance level** $\alpha=P(\text{reject }H_0\mid H_0)$ is chosen in advance, conventionally $0.05$ or $0.01$. The **power** $1-\beta$ depends on how far the truth is from $H_0$, on $n$, and on $\alpha$. For a fixed sample size, $\alpha$ and $\beta$ trade against each other: moving the threshold to reduce false alarms increases misses. Only more data (or less noise) reduces both.

The radar vocabulary is exact, not metaphorical: $\alpha$ is the probability of false alarm $P_\text{FA}$ and $1-\beta$ the probability of detection $P_\text{D}$; plotting $P_\text{D}$ against $P_\text{FA}$ as the threshold varies gives the **receiver operating characteristic** (ROC).

```desmos-graph
left=-4; right=6.5; bottom=0; top=0.45
---
y=\frac{1}{\sqrt{2\pi}}e^{-\frac{x^{2}}{2}}|blue
y=\frac{1}{\sqrt{2\pi}}e^{-\frac{(x-2.5)^{2}}{2}}|red
x=1.96|black|dashed
0\le y\le\frac{1}{\sqrt{2\pi}}e^{-\frac{x^{2}}{2}}|blue|1.96<x<6.5
0\le y\le\frac{1}{\sqrt{2\pi}}e^{-\frac{(x-2.5)^{2}}{2}}|red|-4<x<1.96
```

*Distribution of the standardised statistic under $H_0$ (blue, mean 0) and under an alternative (red, mean 2.5), with threshold $1.96$ (dashed). The blue area to the right is $\alpha=0.025$; the red area to the left is $\beta=0.295$. Moving the threshold right shrinks $\alpha$ and grows $\beta$.*

### The $p$-value

> [!definition] $p$-value
> The **$p$-value** of an observed statistic $t_\text{obs}$ is the probability, computed under $H_0$, of obtaining a value of $T$ at least as extreme (in the direction of $H_1$) as $t_\text{obs}$:
> $$
> p=P_{H_0}(T\ge t_\text{obs})\ \ \text{(upper one-sided)},\qquad p=P_{H_0}(\lvert T\rvert\ge\lvert t_\text{obs}\rvert)\ \ \text{(two-sided, symmetric }T).
> $$
> Rejecting $H_0$ when $p\le\alpha$ is equivalent to the fixed-$\alpha$ test.

The $p$-value measures how surprising the data would be *if $H_0$ were true*. It is not the probability that $H_0$ is true, and it is not the probability that the result is due to chance. A small $p$ says "either $H_0$ is false or something unlikely happened"; a large $p$ says the data are compatible with $H_0$ — which is not the same as evidence *for* $H_0$, particularly when $n$ is small and the test has little power.

> [!theorem] Duality of tests and confidence intervals
> A two-sided level-$\alpha$ test of $H_0:\theta=\theta_0$ rejects exactly when $\theta_0$ lies outside the $100(1-\alpha)\%$ confidence interval built from the same pivot.

*Proof.* Both statements are the same inequality on the pivot: $\lvert(\hat\theta-\theta_0)/\operatorname{SE}\rvert>c$ (reject) is equivalent to $\theta_0\notin[\hat\theta-c\operatorname{SE},\ \hat\theta+c\operatorname{SE}]$. $\square$

The interval is usually more informative than the test result: it says not only *whether* the parameter differs from $\theta_0$ but *by how much* it plausibly differs.

### The general procedure

1. State $H_0$ and $H_1$ (one- or two-sided) *before* looking at the data.
2. Choose $\alpha$, and ideally the sample size needed for adequate power at the smallest effect that matters.
3. Choose the test statistic and check its assumptions (independence, normality, equal variances).
4. Compute the statistic and the $p$-value (or compare with the critical value).
5. State the conclusion in the engineering context, with the size of the effect (a confidence interval), not just "significant".

### Tests for means

> [!theorem] One-sample $z$ and $t$ tests
> For $H_0:\mu=\mu_0$:
> $$
> Z=\frac{\bar X-\mu_0}{\sigma/\sqrt n}\sim\mathcal N(0,1)\ \ (\sigma\text{ known}),\qquad T=\frac{\bar X-\mu_0}{S/\sqrt n}\sim t_{n-1}\ \ (\sigma\text{ unknown, Gaussian data}).
> $$
> Reject at level $\alpha$ when $\lvert Z\rvert>z_{\alpha/2}$ (two-sided), $Z>z_\alpha$ ($H_1:\mu>\mu_0$) or $Z<-z_\alpha$ ($H_1:\mu<\mu_0$); likewise for $T$ with $t_{\alpha/2,n-1}$ or $t_{\alpha,n-1}$.

**Power and sample size.** If the true mean is $\mu_0+\delta$, the $z$ statistic has mean $\delta\sqrt n/\sigma$, so for the one-sided test

$$
1-\beta=Q\Bigl(z_\alpha-\frac{\delta\sqrt n}{\sigma}\Bigr),\qquad n=\Bigl(\frac{(z_\alpha+z_\beta)\,\sigma}{\delta}\Bigr)^2 ,
$$

and for a two-sided test $z_\alpha$ is replaced by $z_{\alpha/2}$ (neglecting the tiny opposite tail). The ratio $\delta/\sigma$ is the **effect size**; detecting an effect of a tenth of a standard deviation needs about a hundred times more data than detecting one of a full standard deviation.

> [!theorem] Two-sample and paired tests
> - **Welch $t$ test** ($H_0:\mu_1=\mu_2$, independent samples, unequal variances): $T=\dfrac{\bar X_1-\bar X_2}{\sqrt{S_1^2/n_1+S_2^2/n_2}}$, approximately $t_\nu$ with the Welch–Satterthwaite $\nu$ of Section 20.7.
> - **Pooled $t$ test** (equal variances assumed): denominator $S_p\sqrt{1/n_1+1/n_2}$, exactly $t_{n_1+n_2-2}$.
> - **Paired $t$ test** (each unit measured under both conditions): one-sample $t$ test on the differences $D_i=X_{1i}-X_{2i}$, with $n-1$ degrees of freedom.

### Chi-square tests

**Test for a variance.** For Gaussian data and $H_0:\sigma^2=\sigma_0^2$, the statistic $\chi^2=(n-1)S^2/\sigma_0^2\sim\chi^2_{n-1}$ under $H_0$. Large values indicate $\sigma^2>\sigma_0^2$.

**Goodness of fit.** Data are classified into $K$ categories with observed counts $O_k$; under $H_0$ the category probabilities are $p_k$, so the expected counts are $E_k=np_k$. **Pearson's statistic**

$$
\chi^2=\sum_{k=1}^K\frac{(O_k-E_k)^2}{E_k}
$$

is approximately $\chi^2_{K-1-m}$ under $H_0$, where $m$ is the number of parameters estimated from the data to compute the $p_k$. *Why chi-square:* each $O_k$ is approximately Gaussian with variance close to $E_k$ (Poisson-like counting), so each term is roughly a squared standard Gaussian; the constraint $\sum O_k=n$ removes one degree of freedom and each fitted parameter removes another. The approximation requires expected counts that are not too small — the usual rule is $E_k\ge5$; merge sparse categories otherwise.

**Independence in a contingency table.** For an $r\times c$ table of counts $O_{ij}$ with row totals $R_i$, column totals $C_j$ and grand total $n$, the expected counts under independence are $E_{ij}=R_iC_j/n$, and $\sum_{i,j}(O_{ij}-E_{ij})^2/E_{ij}$ is approximately $\chi^2_{(r-1)(c-1)}$.

### $F$ test for equality of variances

For two independent Gaussian samples and $H_0:\sigma_1^2=\sigma_2^2$, $F=S_1^2/S_2^2\sim F_{n_1-1,n_2-1}$. For a two-sided test, put the larger variance in the numerator and double the upper-tail probability. Like the chi-square variance test, it is sensitive to non-normality; Levene's test is a robust alternative.

### The Neyman–Pearson lemma

Which test is *best*? For two simple hypotheses the answer is known exactly.

> [!theorem] Neyman–Pearson lemma
> For testing $H_0:\theta=\theta_0$ against $H_1:\theta=\theta_1$, the test that rejects $H_0$ when the **likelihood ratio**
> $$
> \Lambda(\mathbf x)=\frac{L(\theta_1;\mathbf x)}{L(\theta_0;\mathbf x)}>\eta
> $$
> (with $\eta$ chosen so that the size is $\alpha$) has the greatest power among all tests of size at most $\alpha$.

For a Gaussian mean with known $\sigma$, $\ln\Lambda$ is an increasing linear function of $\bar x$ when $\mu_1>\mu_0$, so the optimal test is "reject when $\bar x$ is large" — the $z$ test. The likelihood-ratio test is the foundation of optimal detection in communications and radar: the matched filter followed by a threshold is the Neyman–Pearson detector for a known signal in white Gaussian noise.

> [!example] Worked Example 20.25 — $z$ test and power for a regulator
> A $5.000\ \mathrm V$ regulator line has known output noise between units of $\sigma=20\ \mathrm{mV}$. After a component change, $n=25$ units give $\bar x=5.009\ \mathrm V$. Has the mean shifted? Use $\alpha=0.05$.
>
> **Step 1 — hypotheses.** $H_0:\mu=5.000$, $H_1:\mu\ne5.000$ (two-sided: a shift either way matters).
>
> **Step 2 — statistic.** $\operatorname{SE}=20/\sqrt{25}=4\ \mathrm{mV}$; $z=9/4=2.25$.
>
> **Step 3 — decision.** $\lvert z\rvert>1.96$; $p=2Q(2.25)=0.0244$. Reject $H_0$: the shift is statistically significant. The $95\%$ CI for the shift is $9\pm7.8\ \mathrm{mV}$, i.e. $[1.2,\ 16.8]\ \mathrm{mV}$.
>
> **Step 4 — power.** What is the probability that this test detects a true shift of $\delta=10\ \mathrm{mV}$? The statistic then has mean $10/4=2.5$: power $=Q(1.96-2.5)+\Phi(-1.96-2.5)=0.705$. A $10\ \mathrm{mV}$ shift would be missed $30\%$ of the time.
>
> **Step 5 — sample size for 90% power.** $n=\bigl((1.960+1.282)\times20/10\bigr)^2=42.0$, so 43 units.

> [!example] Worked Example 20.26 — One-sample $t$ test and the danger of low power
> For the eight reference readings of Example 20.3, test $H_0:\mu=2.500\ \mathrm V$.
>
> $t=(2.5025-2.5000)/(0.03240/\sqrt8)=0.0025/0.01146=0.218$ with 7 degrees of freedom; $p=0.83$. There is no evidence of an offset. But the $95\%$ CI is $[2.475,\ 2.530]\ \mathrm V$: offsets as large as $\pm27\ \mathrm{mV}$ are entirely compatible with the data. "Not significant" here means "not enough data to tell", not "no offset".

> [!example] Worked Example 20.27 — Paired versus unpaired comparison of drift
> Six precision resistors are measured (kΩ) before and after a 1000 h burn-in:
>
> | Unit | 1 | 2 | 3 | 4 | 5 | 6 |
> |---|---|---|---|---|---|---|
> | Before | 10.02 | 9.98 | 10.05 | 10.01 | 9.97 | 10.03 |
> | After | 10.06 | 10.01 | 10.08 | 10.06 | 9.99 | 10.07 |
>
> **Step 1 — paired analysis.** Differences $d_i=0.04,\ 0.03,\ 0.03,\ 0.05,\ 0.02,\ 0.04$; $\bar d=0.0350$, $s_d=0.01049$. $t=0.0350/(0.01049/\sqrt6)=8.17$ with 5 degrees of freedom, $p=0.00045$. The drift is unmistakable.
>
> **Step 2 — (wrong) unpaired analysis.** Treating the columns as independent samples gives $t=1.82$, $p=0.10$ — "not significant".
>
> **Step 3 — explanation.** Unit-to-unit variation (about $30\ \Omega$) is much larger than the drift ($35\ \Omega$ on average, but with spread only $10\ \Omega$). Pairing subtracts each unit's own value and removes that variation. Always use the paired test when the same units are measured twice; the unpaired test is not only less powerful but answers a different question.

> [!example] Worked Example 20.28 — Comparing two suppliers: Welch $t$ and $F$ tests
> Capacitance (nF) of samples from two suppliers:
>
> A: $47.2,\ 46.8,\ 47.5,\ 47.1,\ 46.9,\ 47.4,\ 47.0,\ 47.3$ ($n_1=8$, $\bar x_1=47.150$, $s_1=0.245$)
>
> B: $46.1,\ 47.9,\ 46.5,\ 48.2,\ 45.8,\ 47.6,\ 46.9$ ($n_2=7$, $\bar x_2=47.000$, $s_2=0.924$)
>
> **Step 1 — means (Welch).** $\operatorname{SE}=\sqrt{0.245^2/8+0.924^2/7}=0.360$; $t=0.150/0.360=0.417$. Welch–Satterthwaite: $\nu=6.74$. $p=0.69$. No evidence that the means differ.
>
> **Step 2 — variances ($F$).** $F=s_2^2/s_1^2=0.8533/0.0600=14.2$ with $(6,7)$ degrees of freedom. The upper $2.5\%$ point of $F_{6,7}$ is $5.12$, so $H_0:\sigma_1=\sigma_2$ is rejected; the two-sided $p=0.0026$.
>
> **Step 3 — conclusion.** The suppliers agree on average, but supplier B's parts scatter almost four times as widely — exactly the difference that matters for a filter with tight tolerance. Note that the pooled $t$ test would have been inappropriate here, since its equal-variance assumption is clearly violated; Welch's test does not need it.

> [!example] Worked Example 20.29 — Chi-square tests for a bit generator and for shift dependence
> **(a) Goodness of fit.** A hardware random-bit generator produces 400 two-bit symbols with counts $00{:}\ 112$, $01{:}\ 91$, $10{:}\ 85$, $11{:}\ 112$. Under $H_0$ (uniform) each expected count is 100:
> $$
> \chi^2=\frac{12^2+9^2+15^2+12^2}{100}=\frac{594}{100}=5.94,\qquad \nu=4-1=3 .
> $$
> The $5\%$ critical value of $\chi^2_3$ is $7.81$; $p=0.115$. No evidence of bias at this sample size (but a real generator would be tested with millions of symbols and a full test suite).
>
> **(b) Independence.** Defective boards by shift: day $18$ of $400$, night $34$ of $400$. Under independence the expected defectives are $26$ per shift and the expected good boards $374$. Then
> $$
> \chi^2=\frac{8^2}{26}+\frac{8^2}{374}+\frac{8^2}{26}+\frac{8^2}{374}=5.27,\qquad \nu=(2-1)(2-1)=1,
> $$
> exceeding $3.84$; $p=0.022$. The defect rate depends on shift (whether the cause is the shift itself or something correlated with it is a separate, engineering question).

> [!example] Worked Example 20.30 — A borderline variance test
> A sensor's noise is specified as $\sigma\le2.0\ \mathrm{mV}$. Fifteen readings give $s=2.6\ \mathrm{mV}$. Test $H_0:\sigma=2.0$ against $H_1:\sigma>2.0$ at $\alpha=0.05$.
>
> $\chi^2=14\times2.6^2/2.0^2=23.66$; the upper $5\%$ point of $\chi^2_{14}$ is $23.68$; $p=0.0503$. Formally, $H_0$ is not rejected. Practically, a result this close to the threshold says the evidence is marginal; the sensible engineering response is to collect more readings, not to treat $p=0.0503$ and $p=0.0497$ as opposite conclusions.

> [!warning] Common Pitfall — Multiple testing and data dredging
> Testing 20 independent parameters of a product at $\alpha=0.05$ when nothing has changed produces at least one "significant" result with probability $1-0.95^{20}=0.64$. When many tests are made, control the family-wise error, e.g. by the **Bonferroni** rule (use $\alpha/m$ for each of $m$ tests), or the false-discovery rate. Choosing which hypothesis to test *after* seeing which looks interesting invalidates the stated $\alpha$ entirely.

> [!warning] Common Pitfall — Statistical versus practical significance
> With $10^6$ samples, a $0.01\ \mathrm{mV}$ offset in a $5\ \mathrm V$ supply can be "highly significant" ($p<10^{-6}$) and completely irrelevant. With 5 samples, a $50\ \mathrm{mV}$ offset may be "not significant" and very relevant. Always report the estimated effect and its confidence interval alongside the $p$-value.

> [!tip] Engineering Note — Production test limits as hypothesis tests
> A pass/fail limit on a test station is a hypothesis test performed on every unit. Measurement noise causes good units near the limit to fail (false rejects, Type I from the producer's view) and bad units to pass (escapes, Type II). **Guard-banding** moves the test limit inside the specification by a multiple of the measurement uncertainty to cut escapes at the cost of more false rejects — the $\alpha$–$\beta$ trade-off in hardware.

## 20.9 Regression

### Fitting relations between variables

Calibration curves, sensor characterisation, extracting a time constant from a step response, finding a device's temperature coefficient, modelling battery capacity against cycle count: all of these fit a functional relation $y\approx g(x;\boldsymbol\beta)$ to noisy data. **Regression analysis** estimates the parameters $\boldsymbol\beta$, quantifies their uncertainty, and — just as importantly — checks whether the assumed form of $g$ is adequate. This section treats the statistical side; the optimisation side (nonlinear least squares, regularisation) is in Ch. 23.

### The simple linear model

> [!definition] Simple linear regression model
> $$
> Y_i=\beta_0+\beta_1x_i+\varepsilon_i,\qquad i=1,\dots,n,
> $$
> where the $x_i$ are known (not random, or measured without appreciable error) and the errors $\varepsilon_i$ are independent with $E[\varepsilon_i]=0$ and constant variance $\operatorname{Var}(\varepsilon_i)=\sigma^2$. For exact intervals and tests we further assume $\varepsilon_i\sim\mathcal N(0,\sigma^2)$.

### Least-squares estimates

The **least-squares** estimates minimise the sum of squared vertical deviations

$$
S(\beta_0,\beta_1)=\sum_{i=1}^n\bigl(y_i-\beta_0-\beta_1x_i\bigr)^2 .
$$

Setting both partial derivatives to zero (Ch. 6) gives the **normal equations**

$$
\frac{\partial S}{\partial\beta_0}=-2\sum(y_i-\beta_0-\beta_1x_i)=0,\qquad\frac{\partial S}{\partial\beta_1}=-2\sum x_i(y_i-\beta_0-\beta_1x_i)=0 .
$$

> [!theorem] Least-squares line
> $$
> \hat\beta_1=\frac{S_{xy}}{S_{xx}}=r\,\frac{s_y}{s_x},\qquad \hat\beta_0=\bar y-\hat\beta_1\bar x .
> $$
> The fitted line passes through $(\bar x,\bar y)$, and the residuals $e_i=y_i-\hat y_i$ satisfy $\sum e_i=0$ and $\sum x_ie_i=0$.

*Derivation.* The first normal equation gives $\bar y=\beta_0+\beta_1\bar x$, i.e. the line passes through the centroid. Substituting $\beta_0=\bar y-\beta_1\bar x$ into the second and using $\sum x_i(y_i-\bar y)=S_{xy}$, $\sum x_i(x_i-\bar x)=S_{xx}$ gives $S_{xy}-\beta_1S_{xx}=0$. The residual identities are the normal equations themselves. $\square$

Under Gaussian errors, least squares is maximum likelihood (Section 20.6). Without the Gaussian assumption, the **Gauss–Markov theorem** guarantees that least squares gives the *best linear unbiased estimators* (minimum variance among all unbiased estimators that are linear in the $y_i$), provided the errors are uncorrelated with equal variance.

### Sampling properties and inference

Since $\hat\beta_1=\sum(x_i-\bar x)Y_i/S_{xx}$ is a linear combination of the independent $Y_i$,

$$
E[\hat\beta_1]=\beta_1,\qquad\operatorname{Var}(\hat\beta_1)=\frac{\sigma^2}{S_{xx}},\qquad\operatorname{Var}(\hat\beta_0)=\sigma^2\Bigl(\frac1n+\frac{\bar x^2}{S_{xx}}\Bigr),\qquad\operatorname{Cov}(\hat\beta_0,\hat\beta_1)=-\frac{\sigma^2\bar x}{S_{xx}} .
$$

The error variance is estimated from the residuals with $n-2$ degrees of freedom (two parameters fitted):

$$
s^2=\frac{\mathrm{SSE}}{n-2}=\frac{1}{n-2}\sum_{i=1}^n e_i^2 .
$$

> [!theorem] Inference for the regression line (Gaussian errors)
> 1. $\dfrac{\hat\beta_1-\beta_1}{s/\sqrt{S_{xx}}}\sim t_{n-2}$; a $100(1-\alpha)\%$ CI for the slope is $\hat\beta_1\pm t_{\alpha/2,n-2}\,s/\sqrt{S_{xx}}$.
> 2. CI for the **mean response** at $x_0$: $\hat y_0\pm t_{\alpha/2,n-2}\,s\sqrt{\dfrac1n+\dfrac{(x_0-\bar x)^2}{S_{xx}}}$.
> 3. **Prediction interval** for a new observation at $x_0$: $\hat y_0\pm t_{\alpha/2,n-2}\,s\sqrt{1+\dfrac1n+\dfrac{(x_0-\bar x)^2}{S_{xx}}}$.

The term $(x_0-\bar x)^2/S_{xx}$ makes both bands narrowest at $\bar x$ and flare outwards: the line is pinned down best at the centre of the data, and uncertainty in the slope is amplified with distance. This is the quantitative argument against extrapolation. The formula $\operatorname{Var}(\hat\beta_1)=\sigma^2/S_{xx}$ also guides experiment design: to measure a slope precisely, spread the $x$ values as widely as the model remains valid — placing half the points at each end maximises $S_{xx}$, although it leaves no way to detect curvature.

### Analysis of variance and $r^2$

Each observation's deviation from the mean splits into an explained part and a residual: $y_i-\bar y=(\hat y_i-\bar y)+(y_i-\hat y_i)$. Squaring and summing, the cross term vanishes (by the residual identities), giving the **ANOVA decomposition**

$$
\underbrace{\sum(y_i-\bar y)^2}_{\mathrm{SST}}=\underbrace{\sum(\hat y_i-\bar y)^2}_{\mathrm{SSR}}+\underbrace{\sum(y_i-\hat y_i)^2}_{\mathrm{SSE}} .
$$

> [!definition] Coefficient of determination
> $$
> r^2=\frac{\mathrm{SSR}}{\mathrm{SST}}=1-\frac{\mathrm{SSE}}{\mathrm{SST}} ,
> $$
> the fraction of the variation in $y$ explained by the fitted line. In simple linear regression it equals the square of the Pearson correlation $r$ (Section 20.3). With $p$ parameters, the **adjusted** $\bar r^2=1-\dfrac{\mathrm{SSE}/(n-p)}{\mathrm{SST}/(n-1)}$ penalises extra parameters.

The overall significance of the regression is tested by $F=\dfrac{\mathrm{SSR}/1}{\mathrm{SSE}/(n-2)}\sim F_{1,n-2}$ under $H_0:\beta_1=0$; in simple regression $F=t^2$ for the slope.

### Residual analysis

The model assumptions — correct functional form, independent errors, constant variance, Gaussian errors — are checked by examining the residuals $e_i$, ideally the **standardised residuals** $e_i/s$. A good fit leaves residuals that look like structureless noise. Characteristic patterns diagnose specific failures:

| Residual pattern | Plotted against | Diagnosis | Remedy |
|---|---|---|---|
| U-shape or arch | $x$ or $\hat y$ | missing curvature | add $x^2$ term, transform, physical model |
| Funnel (spread grows with $x$) | $x$ or $\hat y$ | non-constant variance | weighted least squares, log transform |
| Runs, waves | time/order | correlated errors, drift | model the drift, randomise run order |
| One very large $\lvert e_i/s\rvert$ ($>3$) | — | outlier | check the point; robust fit |
| Curved normal-probability plot | normal quantiles | non-Gaussian errors | robust methods; intervals approximate |

A point can also be **influential** without having a large residual: a point at extreme $x$ (high **leverage** $h_{ii}=\frac1n+\frac{(x_i-\bar x)^2}{S_{xx}}$) pulls the line towards itself. Cook's distance combines residual and leverage to measure how much the fit changes when a point is deleted.

### Matrix form and multiple regression

With $p$ parameters, write the model as $\mathbf y=\mathbf X\boldsymbol\beta+\boldsymbol\varepsilon$, where the $n\times p$ **design matrix** $\mathbf X$ has a column of ones and one column per regressor (for a polynomial fit: $1,\ x_i,\ x_i^2,\dots$). Minimising $\lVert\mathbf y-\mathbf X\boldsymbol\beta\rVert^2$ gives the normal equations and their solution (Ch. 8, Ch. 23):

$$
\mathbf X^\mathsf T\mathbf X\,\hat{\boldsymbol\beta}=\mathbf X^\mathsf T\mathbf y,\qquad\hat{\boldsymbol\beta}=(\mathbf X^\mathsf T\mathbf X)^{-1}\mathbf X^\mathsf T\mathbf y,\qquad\operatorname{Cov}(\hat{\boldsymbol\beta})=\sigma^2(\mathbf X^\mathsf T\mathbf X)^{-1},\qquad s^2=\frac{\mathrm{SSE}}{n-p}.
$$

Geometrically, $\hat{\mathbf y}=\mathbf X\hat{\boldsymbol\beta}$ is the orthogonal projection of $\mathbf y$ onto the column space of $\mathbf X$, and the residual vector is orthogonal to every column — the matrix form of $\sum e_i=0$, $\sum x_ie_i=0$. The diagonal of the "hat" matrix $\mathbf H=\mathbf X(\mathbf X^\mathsf T\mathbf X)^{-1}\mathbf X^\mathsf T$ gives the leverages. Numerically, solve by QR factorisation rather than by forming $\mathbf X^\mathsf T\mathbf X$, whose condition number is the square of that of $\mathbf X$ (Ch. 22).

If the errors have different variances $\sigma_i^2$, **weighted least squares** minimises $\sum w_i(y_i-\mathbf x_i^\mathsf T\boldsymbol\beta)^2$ with $w_i=1/\sigma_i^2$, giving $\hat{\boldsymbol\beta}=(\mathbf X^\mathsf T\mathbf W\mathbf X)^{-1}\mathbf X^\mathsf T\mathbf W\mathbf y$ — the regression version of the weighted mean of Section 20.2.

### Linearising transformations

Many physical laws become linear after a transformation, and ordinary linear regression can then be used:

| Model | Transformation | Linear form |
|---|---|---|
| $y=Ae^{-t/\tau}$ (RC discharge) | $\ln y$ | $\ln y=\ln A-t/\tau$ |
| $y=Ax^b$ (power law) | $\ln y,\ \ln x$ | $\ln y=\ln A+b\ln x$ |
| $I=I_s\bigl(e^{V/(nV_T)}-1\bigr)$, $I\gg I_s$ (diode) | $\ln I$ | $\ln I=\ln I_s+V/(nV_T)$ |
| $R=R_0e^{B(1/T-1/T_0)}$ (NTC thermistor) | $\ln R,\ 1/T$ | $\ln R=\text{const}+B/T$ |

A transformation also transforms the errors. If $y$ has constant *absolute* noise $\sigma$, then $\ln y$ has noise $\approx\sigma/y$, which grows as $y$ decays; the unweighted log fit then over-weights the small, noisy late values. Use weights $w_i\propto y_i^2$ in the log domain, or fit the nonlinear model directly (Ch. 23).

```desmos-graph
left=-0.5; right=6.5; bottom=-1; top=20
---
(0,1.05)|black
(1,3.12)|black
(2,5.63)|black
(3,8.31)|black
(4,11.45)|black
(5,14.70)|black
(6,18.42)|black
y=0.2661+2.8961x|red|0<x<6
y=1.0321+1.9768x+0.15321x^{2}|blue|dashed|0<x<6
```

*Data of Example 20.33 with the straight-line fit (red, $r^2=0.992$) and the quadratic fit (blue, dashed). Both look acceptable on this scale — the residual plot below shows the difference.*

```desmos-graph
left=-0.5; right=6.5; bottom=-1; top=1
---
(0,0.7839)|red
(1,-0.0421)|red
(2,-0.4282)|red
(3,-0.6443)|red
(4,-0.4004)|red
(5,-0.0464)|red
(6,0.7775)|red
(0,0.0179)|blue
(1,-0.0421)|blue
(2,0.0314)|blue
(3,-0.0314)|blue
(4,0.0593)|blue
(5,-0.0464)|blue
(6,0.0114)|blue
y=0|black|dashed
```

*Residuals of the straight-line fit (red) form a clear U: the model is missing curvature. Residuals of the quadratic fit (blue) are small and patternless.*

> [!example] Worked Example 20.31 — Calibrating an RTD
> Fit $R=\beta_0+\beta_1T$ to the data of Example 20.6 and give the sensitivity with a $95\%$ CI.
>
> **Step 1 — slope and intercept.** From Example 20.6, $S_{TT}=7000$, $S_{TR}=2692$, $\bar T=50$, $\bar R=119.333$. Hence $\hat\beta_1=2692/7000=0.384571\ \Omega/^\circ\mathrm C$ and $\hat\beta_0=119.333-0.384571(50)=100.105\ \Omega$.
>
> **Step 2 — residuals.** $e_i=R_i-\hat R_i=-0.105,\ +0.104,\ +0.012,\ +0.021,\ +0.030,\ -0.062\ \Omega$ (they sum to zero). $\mathrm{SSE}=0.02705\ \Omega^2$.
>
> **Step 3 — error variance.** $s^2=0.02705/4=6.76\times10^{-3}\ \Omega^2$, $s=0.0822\ \Omega$.
>
> **Step 4 — standard errors.** $\operatorname{SE}(\hat\beta_1)=0.0822/\sqrt{7000}=9.83\times10^{-4}\ \Omega/^\circ\mathrm C$; $\operatorname{SE}(\hat\beta_0)=0.0822\sqrt{1/6+2500/7000}=0.0595\ \Omega$.
>
> **Step 5 — CI for the slope.** $t_{0.025,4}=2.776$: $\beta_1\in0.38457\pm0.00273=[0.3818,\ 0.3873]\ \Omega/^\circ\mathrm C$. Relative to $R_0\approx100\ \Omega$, this is a temperature coefficient of about $0.385\%/^\circ\mathrm C$, consistent with the IEC value $0.00385\ /^\circ\mathrm C$ for industrial platinum.
>
> **Step 6 — $r^2$ and the residual pattern.** $r^2=1-0.02705/1035.29=0.99997$. The residuals are negative at both ends and positive in the middle — a faint arch consistent with the known negative quadratic term of platinum. A quadratic fit gives a curvature coefficient $-4.9\times10^{-5}\ \Omega/(^\circ\mathrm C)^2$, close to the Callendar–Van Dusen value $-5.8\times10^{-5}$ for a $100\ \Omega$ sensor. With six points and residuals near the noise level, the evidence for curvature is weak in isolation, but physics supports the quadratic model.

> [!example] Worked Example 20.32 — Using the calibration: inverse prediction and bands
> Using the fit of Example 20.31:
>
> **(a) Mean response at $T=50\ ^\circ\mathrm C$.** $\hat R=119.333\ \Omega$; half-width $2.776\times0.0822\sqrt{1/6}=0.093\ \Omega$.
>
> **(b) Prediction for a single new reading at $50\ ^\circ\mathrm C$.** Half-width $2.776\times0.0822\sqrt{1+1/6}=0.247\ \Omega$.
>
> **(c) Inverse prediction.** The sensor reads $R=120.00\ \Omega$. Solving the fitted line for $T$: $\hat T=(120.00-100.105)/0.384571=51.73\ ^\circ\mathrm C$. Its uncertainty is roughly the prediction half-width converted to temperature: $0.247/0.3846\approx0.64\ ^\circ\mathrm C$ ($95\%$). This "calibrate $y$ on $x$, then invert" procedure is the standard calibration method when the reference temperatures are known much better than the resistance readings.
>
> **(d) Extrapolation.** At $T=300\ ^\circ\mathrm C$ the linear fit predicts $215.48\ \Omega$; the true Pt100 value is about $212.05\ \Omega$. The curvature invisible over $0$–$100\ ^\circ\mathrm C$ causes an error of $3.4\ \Omega$ ($9\ ^\circ\mathrm C$) — far outside any band computed from the fit, because those bands assume the model is correct.

> [!example] Worked Example 20.33 — Detecting curvature from the residuals
> A sensor gives $y=1.05,\ 3.12,\ 5.63,\ 8.31,\ 11.45,\ 14.70,\ 18.42$ at $x=0,1,\dots,6$.
>
> **Step 1 — straight line.** $\bar x=3$, $S_{xx}=28$, $S_{xy}=81.09$: $\hat\beta_1=2.8961$, $\hat\beta_0=0.2661$. $\mathrm{SSE}=1.982$, $\mathrm{SST}=236.82$, $r^2=0.9916$.
>
> **Step 2 — residuals.** $0.78,\ -0.04,\ -0.43,\ -0.64,\ -0.40,\ -0.05,\ 0.78$: positive at the ends, negative in the middle. Despite $r^2>0.99$ the line is systematically wrong.
>
> **Step 3 — quadratic.** With design matrix columns $(1,\ x,\ x^2)$, solving $\mathbf X^\mathsf T\mathbf X\hat{\boldsymbol\beta}=\mathbf X^\mathsf T\mathbf y$ gives $\hat y=1.032+1.977x+0.1532x^2$. Residuals are now all below $0.06$ in magnitude with no pattern, and $\mathrm{SSE}=0.00987$ — two hundred times smaller.
>
> **Step 4 — is the improvement real?** The partial $F$ statistic for adding one parameter is $F=\frac{(1.982-0.00987)/1}{0.00987/4}=799$ on $(1,4)$ degrees of freedom, overwhelmingly significant. Adjusted $\bar r^2$ rises from $0.9900$ to $0.99994$.

> [!example] Worked Example 20.34 — Time constant from an RC discharge
> A capacitor discharges through a resistor; the voltage is sampled at $t=0,1,\dots,5\ \mathrm{ms}$: $5.02,\ 3.01,\ 1.86,\ 1.12,\ 0.67,\ 0.41\ \mathrm V$. Estimate $\tau$.
>
> **Step 1 — linearise.** $y_i=\ln v_i=1.6134,\ 1.1019,\ 0.6206,\ 0.1133,\ -0.4005,\ -0.8916$.
>
> **Step 2 — fit.** $\bar t=2.5$, $S_{tt}=17.5$, $S_{ty}=-8.7698$: slope $\hat\beta_1=-0.50113\ \mathrm{ms^{-1}}$, intercept $\hat\beta_0=1.6124$.
>
> **Step 3 — parameters.** $\hat\tau=-1/\hat\beta_1=1.9955\ \mathrm{ms}$, $\hat V_0=e^{1.6124}=5.015\ \mathrm V$.
>
> **Step 4 — uncertainty.** $s=0.00849$ (in $\ln$ units), $\operatorname{SE}(\hat\beta_1)=0.00849/\sqrt{17.5}=0.00203$. By first-order propagation (Section 20.10), $\operatorname{SE}(\hat\tau)=\operatorname{SE}(\hat\beta_1)/\hat\beta_1^2=0.0081\ \mathrm{ms}$. So $\tau=1.996\pm0.022\ \mathrm{ms}$ ($95\%$, $t_{0.025,4}=2.776$).
>
> **Step 5 — check.** The residuals in the log domain ($+0.001,\ -0.009,\ +0.011,\ +0.004,\ -0.008,\ +0.002$) show no trend, so a single exponential is adequate. If the noise were constant in volts, the later points (small $v$) would show larger log-residuals and a weighted fit would be preferable.

> [!warning] Common Pitfall — $r^2$ is not a test of model adequacy
> A high $r^2$ says the line explains most of the variation *in this data set*; it says nothing about whether the functional form is right (Examples 20.31 and 20.33), whether it can be extrapolated, or whether the errors are independent. $r^2$ also depends on the range of $x$: the same physical relation sampled over a wider range gives a higher $r^2$. Judge fits by residual plots and by the size of $s$ relative to the required accuracy.

> [!warning] Common Pitfall — Errors in $x$ and regression dilution
> Least squares assumes $x$ is known exactly. If $x$ is measured with noise of variance $\sigma_u^2$, the fitted slope is biased towards zero by the factor $S_{xx}/(S_{xx}+(n-1)\sigma_u^2)$ (**regression dilution**). When both variables are noisy and comparably so, use errors-in-variables methods (orthogonal/Deming regression) or make sure the reference variable is far more accurate than the response.

> [!tip] Engineering Note — Linearity specifications
> An ADC's or sensor's "integral non-linearity" is the maximum residual from a fitted straight line, expressed in LSB or % of full scale. Datasheets specify whether the line is the *end-point* line (through the first and last codes) or the *best-fit* (least-squares) line; the best-fit figure is always smaller, sometimes by a factor of two, so the two are not comparable.

## 20.10 Experimental uncertainty

### Why a common framework

A measurement result is incomplete without a statement of its uncertainty: without one, results cannot be compared with each other, with specifications, or with theory. Before 1993, laboratories used incompatible conventions ("random" and "systematic" errors combined linearly by some, in quadrature by others, at $1\sigma$ by some and $3\sigma$ by others). The *Guide to the Expression of Uncertainty in Measurement* (**GUM**, JCGM 100), issued jointly by BIPM, IEC, ISO and other bodies, established one method that is now used by every accredited calibration laboratory. It is the statistics of this chapter applied systematically to measurement.

### Vocabulary

> [!definition] Measurand, error and uncertainty
> - The **measurand** $Y$ is the quantity intended to be measured (the resistance of *this* resistor at $23\ ^\circ\mathrm C$).
> - The **measurement model** $Y=f(X_1,\dots,X_N)$ expresses the measurand in terms of **input quantities** $X_i$ (readings, calibration corrections, environmental quantities).
> - The **error** is the difference between the measured value and the true value; it is unknowable.
> - The **standard uncertainty** $u(x_i)$ of an input estimate $x_i$ is the standard deviation of the probability distribution that describes our knowledge of $X_i$.
> - The **combined standard uncertainty** $u_c(y)$ is the standard uncertainty of the output estimate $y=f(x_1,\dots,x_N)$.
> - The **expanded uncertainty** $U=k\,u_c(y)$ defines an interval $y\pm U$ expected to contain a large fraction (the **coverage probability**, typically $95\%$) of the values that could reasonably be attributed to $Y$; $k$ is the **coverage factor**.

Two further terms are often confused. **Accuracy** (closeness to the true value) is a qualitative concept; **precision** refers to the scatter of repeated readings (repeatability under the same conditions, reproducibility under changed conditions). Uncertainty is the quantitative statement that covers both. Known systematic effects are *corrected* — the GUM assumes the result has been corrected for all recognised significant systematic effects — and the uncertainty of each correction becomes an input to the budget.

### Type A and Type B evaluation

The GUM classifies uncertainty components not by their physical nature ("random" vs "systematic") but by the **method used to evaluate them**.

> [!definition] Type A evaluation
> Evaluation by statistical analysis of a series of observations. For $n$ independent repeated readings $q_1,\dots,q_n$ of an input quantity, the estimate is $\bar q$ and
> $$
> u(\bar q)=\frac{s(q)}{\sqrt n},\qquad \nu=n-1\ \text{degrees of freedom}.
> $$

This is the standard error of the mean of Section 20.5. If a well-characterised *pooled* standard deviation $s_p$ from earlier extensive experiments is available, it can replace $s(q)$ (with its larger degrees of freedom), which is valuable when only a few readings are taken in the measurement itself.

> [!definition] Type B evaluation
> Evaluation by any other means: manufacturer's specifications, calibration certificates, data from handbooks, resolution of a display, experience with similar instruments. The available information is converted to a standard deviation by assuming a probability distribution.

| Information available | Assumed distribution | Standard uncertainty |
|---|---|---|
| Certificate: "$U$ at $k=2$" (or "95%") | Gaussian | $u=U/k$ (or $U/1.96$) |
| Limits $\pm a$, no other knowledge (specification, tolerance) | rectangular (uniform) | $u=a/\sqrt3\approx0.577a$ |
| Limits $\pm a$, values near centre more likely | triangular | $u=a/\sqrt6\approx0.408a$ |
| Quantity varies sinusoidally between $\pm a$ (e.g. RF mismatch, temperature cycling) | U-shaped (arcsine) | $u=a/\sqrt2\approx0.707a$ |
| Digital resolution $\delta$ (one count) | rectangular, half-width $\delta/2$ | $u=\delta/\sqrt{12}\approx0.289\delta$ |

The factors come from the variances of the distributions: uniform on $[-a,a]$ has variance $a^2/3$ (Section 20.4); the symmetric triangular distribution has variance $a^2/6$; the arcsine distribution of $a\sin\theta$ with $\theta$ uniform has variance $E[a^2\sin^2\theta]=a^2/2$.

Type A and Type B uncertainties are both standard deviations and are combined in the same way. A Type B component with well-founded limits may be assigned effectively infinite degrees of freedom. When the limits themselves are uncertain, the GUM gives the approximation $\nu\approx\frac12\bigl(\Delta u/u\bigr)^{-2}$, where $\Delta u/u$ is the relative uncertainty of the uncertainty (e.g. limits believed reliable to $25\%$ give $\nu=8$).

### The law of propagation of uncertainty

Linearise the model about the input estimates (Ch. 4, Ch. 6):

$$
Y\approx f(x_1,\dots,x_N)+\sum_{i=1}^N c_i\,(X_i-x_i),\qquad c_i=\left.\frac{\partial f}{\partial X_i}\right\rvert_{\mathbf x} .
$$

The $c_i$ are **sensitivity coefficients**; they convert the uncertainty of each input into units of the output. Taking the variance of the linear approximation, with $\mathbf c=(c_1,\dots,c_N)^\mathsf T$ and the input covariance matrix $\boldsymbol\Sigma_x$, gives Section 20.3's rule $\operatorname{Var}(\mathbf c^\mathsf T\mathbf X)=\mathbf c^\mathsf T\boldsymbol\Sigma_x\mathbf c$:

> [!theorem] Law of propagation of uncertainty (GUM)
> $$
> u_c^2(y)=\sum_{i=1}^N c_i^2u^2(x_i)+2\sum_{i=1}^{N-1}\sum_{j=i+1}^{N}c_ic_j\,u(x_i,x_j),
> $$
> where $u(x_i,x_j)=r(x_i,x_j)\,u(x_i)\,u(x_j)$ is the estimated covariance. For **uncorrelated** inputs,
> $$
> u_c(y)=\sqrt{\sum_{i=1}^N\bigl(c_i\,u(x_i)\bigr)^2} ,
> $$
> the root-sum-of-squares (RSS) of the individual contributions $u_i(y)=\lvert c_i\rvert u(x_i)$.

Two special cases are used constantly:

- **Sums and differences** $Y=\sum a_iX_i$: $c_i=a_i$, so $u_c^2=\sum a_i^2u^2(x_i)$ (plus covariance terms).
- **Products and quotients** $Y=A\,X_1^{p_1}X_2^{p_2}\cdots X_N^{p_N}$: $c_i=p_iy/x_i$, so the **relative** uncertainties combine:
$$
\Bigl(\frac{u_c(y)}{y}\Bigr)^2=\sum_{i=1}^N\Bigl(p_i\frac{u(x_i)}{x_i}\Bigr)^2 .
$$

The RSS combination has an important practical consequence: the largest contributions dominate. A component one third the size of the largest adds only about $5\%$ to $u_c$. Effort spent reducing small components is wasted; the **uncertainty budget** shows where improvement is worthwhile.

**Correlated inputs.** If two inputs are fully correlated ($r=+1$, e.g. two voltages measured with the same uncalibrated meter on the same range), their contributions add *linearly*: $u_c=\lvert c_1u_1+c_2u_2\rvert$. For a sum of two such readings this increases the uncertainty; for a difference or a ratio it can cancel it — the reason ratio and substitution methods (bridges, comparators, ratiometric ADC references) are so accurate.

**Validity.** The law of propagation is first order. It is accurate when the model is nearly linear over the range $x_i\pm u(x_i)$, i.e. when the second-order terms $\frac12\sum\frac{\partial^2f}{\partial X_i^2}u^2(x_i)$ are negligible compared with $u_c$. It fails, for example, for $Y=X^2$ when $x\approx0$ (the linear term vanishes) or for strongly nonlinear functions such as decibel conversions of very noisy quantities. GUM Supplement 1 then prescribes the **Monte Carlo method**: draw $M$ (typically $10^6$) samples of each input from its assigned distribution, evaluate the model for each draw, and take the mean, standard deviation and the $2.5\%$/$97.5\%$ percentiles of the resulting output sample (Ch. 24).

### Effective degrees of freedom and the coverage factor

If $u_c$ is dominated by Type A components with few degrees of freedom, $u_c$ is itself uncertain and $k=2$ gives less than $95\%$ coverage. The output is then treated as $t$-distributed with an effective number of degrees of freedom given by the **Welch–Satterthwaite formula**:

$$
\nu_\text{eff}=\frac{u_c^4(y)}{\displaystyle\sum_{i=1}^N\frac{u_i^4(y)}{\nu_i}} ,\qquad u_i(y)=\lvert c_i\rvert u(x_i) .
$$

It is obtained by matching the variance of the estimated $u_c^2$ to that of a scaled $\chi^2_\nu$ variable: each $u_i^2$ estimated with $\nu_i$ degrees of freedom has $\operatorname{Var}(u_i^2)\approx2u_i^4/\nu_i$ (Section 20.5), these add, and a $\chi^2_\nu$-based estimate of $u_c^2$ would have variance $2u_c^4/\nu$. The coverage factor is then $k=t_{0.025,\nu_\text{eff}}$ (truncating $\nu_\text{eff}$ down to an integer): for $\nu_\text{eff}\ge30$, $k\approx2$; for $\nu_\text{eff}=5$, $k=2.57$.

### Reporting

The GUM recommends stating the result as, for example,

$$
R=4.9748\ \mathrm{k\Omega}\pm0.0045\ \mathrm{k\Omega},
$$

"where the number following the symbol $\pm$ is the numerical value of an expanded uncertainty $U=ku_c$, with $U$ determined from a combined standard uncertainty $u_c=2.3\ \Omega$ and a coverage factor $k=2$, corresponding to a coverage probability of approximately $95\%$." Uncertainties are rounded to at most two significant digits, and the result is rounded to the same decimal position. A compact equivalent with standard uncertainty is $R=4.9748(23)\ \mathrm{k\Omega}$, the digits in parentheses being $u_c$ in units of the last digit.

> [!example] Worked Example 20.35 — Type B evaluations from a datasheet and a certificate
> **(a) DMM specification.** On the $10\ \mathrm V$ range a DMM is specified to $\pm(0.02\%\ \text{of reading}+5\ \text{counts})$, with resolution $0.1\ \mathrm{mV}$. For a reading of $4.9852\ \mathrm V$: $a=0.0002(4.9852)+0.0005=1.497\ \mathrm{mV}$. With no further knowledge, rectangular: $u=1.497/\sqrt3=0.864\ \mathrm{mV}$.
>
> **(b) Resolution.** One count $\delta=0.1\ \mathrm{mV}$: $u=0.1/\sqrt{12}=0.029\ \mathrm{mV}$ — negligible beside (a). (The DMM specification usually already includes resolution; it is listed separately here only to show its size.)
>
> **(c) Certificate.** A standard resistor's certificate states $U=0.8\ \mu\Omega/\Omega$ at $k=2$. Then $u=0.4\ \mu\Omega/\Omega$, with $\nu$ taken from the certificate or as large.
>
> **(d) RF mismatch.** A power meter's mismatch uncertainty has limits $\pm0.10\ \mathrm{dB}$, arising from an unknown relative phase of reflections, so the U-shaped distribution applies: $u=0.10/\sqrt2=0.071\ \mathrm{dB}$ — larger than the rectangular value $0.058\ \mathrm{dB}$, because the error spends most of its time near the limits.

> [!example] Worked Example 20.36 — Propagation for a product: power in a resistor
> $P=V^2/R$ with $V=12.00\ \mathrm V$, $u(V)=0.02\ \mathrm V$ and $R=47.0\ \Omega$, $u(R)=0.2\ \Omega$ (uncorrelated).
>
> **Step 1 — estimate.** $P=144/47=3.0638\ \mathrm W$.
>
> **Step 2 — relative uncertainties.** $u(V)/V=0.167\%$, $u(R)/R=0.426\%$.
>
> **Step 3 — combine with exponents $p_V=2$, $p_R=-1$.**
> $$
> \frac{u_c(P)}{P}=\sqrt{(2\times0.167\%)^2+(0.426\%)^2}=\sqrt{0.111+0.181}\,\%=0.541\% .
> $$
> So $u_c(P)=0.0166\ \mathrm W$ and $U=0.033\ \mathrm W$ ($k=2$): $P=3.064\pm0.033\ \mathrm W$.
>
> **Step 4 — sensitivity coefficients (the same thing, explicitly).** $c_V=2V/R=0.5106\ \mathrm{W/V}$, $c_R=-V^2/R^2=-0.06519\ \mathrm{W/\Omega}$; contributions $0.0102\ \mathrm W$ and $0.0130\ \mathrm W$; RSS $=0.0166\ \mathrm W$.
>
> **Step 5 — Monte Carlo check.** Drawing $10^6$ Gaussian samples of $V$ and $R$ and evaluating $V^2/R$ gives mean $3.0639\ \mathrm W$, standard deviation $0.0165\ \mathrm W$ and a $95\%$ interval $[3.0316,\ 3.0964]\ \mathrm W$, confirming the linearised result: the model is nearly linear over $\pm u$.

> [!example] Worked Example 20.37 — Correlated inputs: series and ratio
> Two nominally equal $1\ \mathrm{k\Omega}$ resistors are measured with the same ohmmeter, whose calibration uncertainty $0.5\ \Omega$ is common to both readings (correlation $r=1$); random reading noise is negligible.
>
> **(a) Series combination** $R_s=R_1+R_2$: $c_1=c_2=1$, so $u_c^2=0.25+0.25+2(1)(1)(0.5)(0.5)=1.0\ \Omega^2$, $u_c=1.0\ \Omega$. Ignoring the correlation would give $0.71\ \Omega$ — an underestimate.
>
> **(b) Ratio** $Q=R_1/R_2$: write the calibration error as a common scale factor. With $c_1=1/R_2$ and $c_2=-R_1/R_2^2$ (both $\approx\pm10^{-3}\ \Omega^{-1}$), $u_c^2\approx(10^{-3}\times0.5)^2+(10^{-3}\times0.5)^2-2(10^{-3})(10^{-3})(0.25)=0$. The common error cancels completely in the ratio. This is why a resistance bridge or a ratiometric measurement (sensor and ADC sharing one reference) can be far more accurate than the absolute accuracy of its parts.

> [!example] Worked Example 20.38 — A complete uncertainty budget: $R=V/I$
> A resistor is measured by passing a current and reading voltage and current with two DMMs. Ten simultaneous readings give $\bar V=4.9852\ \mathrm V$ with $s(V)=1.2\ \mathrm{mV}$ and $\bar I=1.00210\ \mathrm{mA}$ with $s(I)=0.30\ \mu\mathrm A$. The voltmeter is specified as in Example 20.35(a); the ammeter as $\pm(0.05\%\ \text{of reading}+0.2\ \mu\mathrm A)$. Assume the two meters' errors are uncorrelated, as is the reading noise.
>
> **Step 1 — model and estimate.** $R=V/I=4.9852/1.00210=4.97475\ \mathrm{k\Omega}$.
>
> **Step 2 — sensitivity coefficients.** $c_V=1/I=0.99790\ \mathrm{k\Omega/V}$; $c_I=-V/I^2=-4.9643\ \mathrm{k\Omega/mA}$.
>
> **Step 3 — input uncertainties.**
> - $V$, Type A: $1.2/\sqrt{10}=0.379\ \mathrm{mV}$, $\nu=9$.
> - $V$, Type B (specification): $0.864\ \mathrm{mV}$, $\nu=\infty$.
> - $I$, Type A: $0.30/\sqrt{10}=0.0949\ \mu\mathrm A$, $\nu=9$.
> - $I$, Type B: $a=0.0005(1.00210\ \mathrm{mA})+0.2\ \mu\mathrm A=0.701\ \mu\mathrm A$; $u=0.701/\sqrt3=0.405\ \mu\mathrm A$, $\nu=\infty$.
>
> **Step 4 — budget.**
>
> | Source | $u(x_i)$ | $c_i$ | $\lvert u_i(R)\rvert$ (Ω) | $\nu_i$ |
> |---|---|---|---|---|
> | $V$ repeatability (A) | $0.379\ \mathrm{mV}$ | $0.9979\ \mathrm{k\Omega/V}$ | 0.38 | 9 |
> | $V$ meter spec (B, rect.) | $0.864\ \mathrm{mV}$ | $0.9979\ \mathrm{k\Omega/V}$ | 0.86 | $\infty$ |
> | $I$ repeatability (A) | $0.0949\ \mu\mathrm A$ | $-4.964\ \mathrm{k\Omega/mA}$ | 0.47 | 9 |
> | $I$ meter spec (B, rect.) | $0.405\ \mu\mathrm A$ | $-4.964\ \mathrm{k\Omega/mA}$ | 2.01 | $\infty$ |
> | **Combined** $u_c(R)$ | | | **2.27** | $\nu_\text{eff}\approx3400$ |
>
> **Step 5 — combine.** $u_c=\sqrt{0.38^2+0.86^2+0.47^2+2.01^2}=2.27\ \Omega$ (relative $0.046\%$).
>
> **Step 6 — effective degrees of freedom.** $\nu_\text{eff}=2.27^4/\bigl(0.38^4/9+0.47^4/9\bigr)\approx3.4\times10^3$; Type B components dominate, so $k=2$ is appropriate.
>
> **Step 7 — report.** $R=4.9748\ \mathrm{k\Omega}\pm0.0045\ \mathrm{k\Omega}$ ($k=2$, $\approx95\%$).
>
> **Step 8 — read the budget.** The ammeter specification contributes $2.01^2/2.27^2=78\%$ of the variance. Taking more readings (reducing the Type A terms) would barely help; a better ammeter, or measuring the current as the voltage across a calibrated standard resistor (turning the measurement into a ratio, Example 20.37), would.

> [!example] Worked Example 20.39 — When few readings dominate: Welch–Satterthwaite
> A frequency is measured four times; the Type A standard uncertainty of the mean is $u_A=0.8\ \mathrm{Hz}$ ($\nu=3$). The counter's timebase contributes a Type B $u_B=0.5\ \mathrm{Hz}$ ($\nu=\infty$).
>
> **Step 1 — combine.** $u_c=\sqrt{0.8^2+0.5^2}=0.943\ \mathrm{Hz}$.
>
> **Step 2 — degrees of freedom.** $\nu_\text{eff}=0.943^4/(0.8^4/3)=0.792/0.1365=5.8$, truncated to 5.
>
> **Step 3 — coverage factor.** $k=t_{0.025,5}=2.57$, so $U=2.57\times0.943=2.4\ \mathrm{Hz}$. Using $k=2$ would give $1.9\ \mathrm{Hz}$ and an actual coverage of only about $90\%$.
>
> **Step 4 — improvement.** Taking 16 readings instead of 4 would halve $u_A$ to $0.4\ \mathrm{Hz}$ (if $s$ is unchanged) and raise $\nu$ to 15; then $u_c=0.64\ \mathrm{Hz}$, $\nu_\text{eff}\approx98$, $k\approx1.98$ and $U\approx1.3\ \mathrm{Hz}$ — nearly half the previous value.

> [!warning] Common Pitfall — Adding uncertainties linearly, or double counting
> Adding all contributions linearly ("worst case") overstates the uncertainty for independent inputs; RSS is correct unless the inputs are correlated. Conversely, listing the same physical effect twice — e.g. a DMM's resolution in addition to a specification that already includes it, or the scatter of readings that already reflects the temperature fluctuation you also list as Type B — inflates the budget. Each physical effect should appear once.

> [!warning] Common Pitfall — Confusing tolerance with uncertainty
> A $1\%$ resistor's tolerance is the limit of its deviation from nominal, a Type B input with $u=1\%/\sqrt3=0.58\%$ if nothing else is known. Once the resistor has been *measured*, its value is known to the measurement uncertainty, which is usually far smaller; the tolerance is then irrelevant. Using tolerance limits as if they were $1\sigma$ values, or as if they were $95\%$ limits, are both common errors.

> [!tip] Engineering Note — Uncertainty and conformance decisions
> Deciding whether a measured value conforms to a specification must account for the measurement uncertainty: a reading $0.3\%$ inside a $\pm1\%$ limit with $U=0.5\%$ does not prove conformance. Standards such as ISO 14253-1 and ILAC-G8 formalise this with **guard bands** of width related to $U$ — the measurement-uncertainty version of the production-test trade-off in Section 20.8. A common rule of thumb is a **test uncertainty ratio** of at least $4{:}1$ between the specification tolerance and the measurement's expanded uncertainty.

## Summary

| Topic | Key formula or fact |
|---|---|
| Sample mean, variance | $\bar x=\frac1n\sum x_i$, $\ s^2=\frac{1}{n-1}\sum(x_i-\bar x)^2$ (unbiased) |
| Robust spread | $\mathrm{IQR}=Q_3-Q_1$; fences $Q_1-1.5\,\mathrm{IQR}$, $Q_3+1.5\,\mathrm{IQR}$; $\sigma\approx1.4826\,\mathrm{MAD}$ |
| Weighted mean | $\bar x_w=\sum w_ix_i/\sum w_i$, $w_i=1/u_i^2$, $u(\bar x_w)=(\sum w_i)^{-1/2}$ |
| RMS and variance | $V_\mathrm{rms}^2=V_\mathrm{DC}^2+\sigma_v^2$ |
| Correlation | $r=S_{xy}/\sqrt{S_{xx}S_{yy}}$, $-1\le r\le1$; Spearman $r_s=1-6\sum d_i^2/(n(n^2-1))$ |
| Covariance propagation | $\boldsymbol\Sigma_Y=\mathbf A\boldsymbol\Sigma_X\mathbf A^\mathsf T$; $\operatorname{Var}(X\pm Y)=\sigma_X^2+\sigma_Y^2\pm2\operatorname{Cov}$ |
| Binomial / Poisson | $np,\ np(1-p)$ / $\lambda,\ \lambda$ |
| Uniform / exponential | $\frac{a+b}{2},\ \frac{(b-a)^2}{12}$ / $\frac1\lambda,\ \frac1{\lambda^2}$; memoryless |
| Rayleigh | $E[R]=\sigma\sqrt{\pi/2}$, $E[R^2]=2\sigma^2$, $P(R>r)=e^{-r^2/(2\sigma^2)}$ |
| $\chi^2_k$, $t_\nu$, $F_{d_1,d_2}$ | sum of $k$ squared $\mathcal N(0,1)$; $Z/\sqrt{V/\nu}$; ratio of normalised $\chi^2$ |
| Sample mean | $E\bar X=\mu$, $\operatorname{SE}=\sigma/\sqrt n$; CLT: $\frac{\bar X-\mu}{\sigma/\sqrt n}\to\mathcal N(0,1)$ |
| Gaussian sampling | $(n-1)S^2/\sigma^2\sim\chi^2_{n-1}$; $\frac{\bar X-\mu}{S/\sqrt n}\sim t_{n-1}$; $\bar X\perp S^2$ |
| Estimator quality | $\operatorname{MSE}=\operatorname{Var}+b^2$; CRLB $\operatorname{Var}\hat\theta\ge1/(nI(\theta))$ |
| MLE | maximise $\ell(\theta)=\sum\ln f(x_i;\theta)$; $\operatorname{SE}\approx1/\sqrt{-\ell''(\hat\theta)}$; Gaussian noise $\Rightarrow$ least squares |
| CI for mean | $\bar x\pm z_{\alpha/2}\sigma/\sqrt n$; $\ \bar x\pm t_{\alpha/2,n-1}s/\sqrt n$; $n=(z_{\alpha/2}\sigma/E)^2$ |
| CI for variance | $\bigl[(n-1)s^2/\chi^2_{\alpha/2,n-1},\ (n-1)s^2/\chi^2_{1-\alpha/2,n-1}\bigr]$ |
| Tests | $z$, $t$ (one-sample, Welch, paired); $\chi^2=\sum(O-E)^2/E$ with $K-1-m$ d.o.f.; $F=s_1^2/s_2^2$ |
| Errors | $\alpha=P(\text{reject}\mid H_0)$, $\beta=P(\text{accept}\mid H_1)$; power $=1-\beta$; $n=((z_\alpha+z_\beta)\sigma/\delta)^2$ |
| Regression | $\hat\beta_1=S_{xy}/S_{xx}$, $\hat\beta_0=\bar y-\hat\beta_1\bar x$, $s^2=\mathrm{SSE}/(n-2)$, $\operatorname{SE}(\hat\beta_1)=s/\sqrt{S_{xx}}$ |
| $r^2$, matrix LS | $r^2=1-\mathrm{SSE}/\mathrm{SST}$; $\hat{\boldsymbol\beta}=(\mathbf X^\mathsf T\mathbf X)^{-1}\mathbf X^\mathsf T\mathbf y$, $\operatorname{Cov}=\sigma^2(\mathbf X^\mathsf T\mathbf X)^{-1}$ |
| Type A / Type B | $u=s/\sqrt n$ ($\nu=n-1$); rectangular $a/\sqrt3$, triangular $a/\sqrt6$, U-shaped $a/\sqrt2$, resolution $\delta/\sqrt{12}$ |
| Propagation (GUM) | $u_c^2=\sum c_i^2u^2(x_i)+2\sum_{i<j}c_ic_ju(x_i,x_j)$, $c_i=\partial f/\partial x_i$ |
| Products | $(u_c/y)^2=\sum(p_iu(x_i)/x_i)^2$ for $y=A\prod x_i^{p_i}$ |
| Expanded uncertainty | $U=ku_c$; $\nu_\text{eff}=u_c^4/\sum(u_i^4/\nu_i)$; $k=t_{0.025,\nu_\text{eff}}$ ($\approx2$ for large $\nu$) |

## Further reading

- D. C. Montgomery and G. C. Runger, *Applied Statistics and Probability for Engineers*, 7th ed., Wiley, 2018 — the standard engineering statistics text; thorough on intervals, tests and regression.
- A. Papoulis and S. U. Pillai, *Probability, Random Variables and Stochastic Processes*, 4th ed., McGraw-Hill, 2002 — distributions, estimation and the bridge to random processes (Ch. 21).
- S. M. Kay, *Fundamentals of Statistical Signal Processing, Vol. I: Estimation Theory*, Prentice Hall, 1993 — Cramér–Rao bounds and maximum likelihood with signal-processing examples.
- JCGM 100:2008, *Evaluation of Measurement Data — Guide to the Expression of Uncertainty in Measurement* (GUM), and JCGM 101:2008 (Supplement 1, Monte Carlo method), freely available from the BIPM — the authoritative source for Section 20.10.
- I. G. Hughes and T. P. A. Hase, *Measurements and their Uncertainties*, Oxford University Press, 2010 — a concise practical guide to error analysis and curve fitting.
