# Chapter 21 — Random Processes

[← Previous](20-statistical-analysis.md) · [Table of Contents](../00-front/table-of-contents.md) · [Next →](../part-7-numerical-computational/22-numerical-methods.md)

> **Prerequisites:** Ch. 16 (Fourier Transforms), Ch. 17 (Z-Transforms), Ch. 19 (Probability), Ch. 20 (Statistical Analysis)
> **Learning outcomes:** After this chapter you can …
> - describe a random process as an ensemble of sample functions, distinguish ensemble averages from time averages, and compute mean, autocorrelation and autocovariance functions;
> - test a process for strict- and wide-sense stationarity and use the properties of autocorrelation and cross-correlation functions;
> - move between autocorrelation and power spectral density with the Wiener–Khinchin theorem, and compute the output mean, autocorrelation and PSD $S_{YY}=\lvert H\rvert^2S_{XX}$ of an LTI system with random input, in continuous and discrete time (including AR(1) sequences);
> - work with Gaussian processes, white and band-limited white noise, and the noise-equivalent bandwidth of a filter;
> - decide when time averages may replace ensemble averages (ergodicity) and how long a record must be;
> - carry out standard noise calculations: thermal noise $4kTRB$, $kT/C$ noise, SNR after filtering and averaging, and the matched filter.

A random variable (Ch. 19) assigns a number to each outcome of an experiment. Many engineering signals — the hiss at a receiver output, a speech waveform, the bit stream on a serial link, the drift of a gyroscope, the error of an ADC over time — are not single numbers but *functions of time* whose exact shape cannot be predicted in advance. Yet they are far from lawless: the noise at the output of a given amplifier always has about the same RMS value, its spectrum always has about the same shape, and a filter always reduces it by about the same amount. A **random process** is the mathematical object that captures this combination of unpredictable detail and predictable statistics. It assigns a whole waveform to each outcome of an experiment.

This chapter develops the **second-order theory** of random processes — means, correlation functions and power spectra. Second-order theory is remarkably economical: it is exactly what is preserved and transformed by linear systems, it completely characterizes Gaussian processes, and it answers the questions an engineer asks most often (What is the noise power? How is it distributed over frequency? How much does this filter remove? How long must I average?). We use the Fourier transform conventions of [Ch. 16](../part-5-transforms/16-fourier-transforms.md) ($X(\omega)=\int x(t)e^{-j\omega t}dt$) and the discrete-time transforms of [Ch. 17](../part-5-transforms/17-z-transforms.md) throughout, and the probability language of [Ch. 19](19-probability.md).

## 21.1 Random sequences and processes

### Motivation

Suppose you connect an oscilloscope across a resistor at room temperature, set a large vertical gain, and press "single" a hundred times. Every trace is different — no two records of thermal noise ever coincide — but every trace has the same "look": the same typical amplitude, the same graininess. If you sample all hundred traces at the instant $t_1=3\ \mu\text{s}$ after the trigger, the hundred numbers you collect form a histogram that looks Gaussian with zero mean. If you instead take *one* trace and histogram its values along time, you get, for this kind of noise, the same Gaussian. These two procedures — averaging *across* records at a fixed time, and averaging *along* one record — are the two fundamental operations of this chapter, and the question of when they agree (ergodicity, §21.7) is one of its central themes.

### Definition and the two viewpoints

> [!definition] Random process
> Let $(\Omega,\mathcal{F},P)$ be a probability space and $T$ an index set. A **random process** (stochastic process) is a family of random variables $\{X(t),\ t\in T\}$ all defined on the same $\Omega$. For each outcome $\zeta\in\Omega$ the function $t\mapsto X(t,\zeta)$ is a **sample function** (realization, sample path); the collection of all sample functions is the **ensemble**.
> - $T\subseteq\mathbb{R}$ (an interval): **continuous-time** process $X(t)$.
> - $T\subseteq\mathbb{Z}$: **discrete-time** process or **random sequence** $X[n]$.
>
> Independently, the *values* may be continuous (a noise voltage) or discrete (a bit, a count).

A random process is thus a function of two variables, $X(t,\zeta)$, and there are two ways to slice it:

- **Fix $t$:** $X(t_1,\cdot)$ is an ordinary random variable with CDF $F_X(x;t_1)=P\{X(t_1)\le x\}$ and pdf $f_X(x;t_1)=\partial F_X/\partial x$. Everything in Ch. 19 applies to it.
- **Fix $\zeta$:** $X(\cdot,\zeta_1)$ is a deterministic waveform — the single trace you see on the oscilloscope. Everything in Chs. 15–17 (Fourier series, transforms) applies to it, at least formally.

We follow the usual convention of suppressing $\zeta$ and writing $X(t)$; lowercase $x(t)$ denotes one particular sample function.

The four combinations of continuous/discrete time and continuous/discrete values all occur in practice:

| | continuous values | discrete values |
|---|---|---|
| **continuous time** | thermal noise voltage $V_n(t)$ | random telegraph signal $\pm1$; photon count $N(t)$ |
| **discrete time** | sampled noise $V_n[n]=V_n(nT_s)$ | data bits $b[n]\in\{0,1\}$; quantized samples |

### Finite-dimensional distributions

A single random variable is described by its CDF. A random process involves infinitely many random variables, and their joint behaviour matters: the probability that a noise voltage exceeds a threshold at two instants 1 ns apart is very different from the probability for two instants 1 s apart.

> [!definition] Finite-dimensional distributions
> The **$k$-th order distribution** of $X(t)$ is the joint CDF
> $$
> F_X(x_1,\dots,x_k;\,t_1,\dots,t_k) = P\{X(t_1)\le x_1,\dots,X(t_k)\le x_k\},
> $$
> with joint pdf $f_X(x_1,\dots,x_k;t_1,\dots,t_k)$. The process is **specified in distribution** when these are known for every $k$ and every choice of times $t_1,\dots,t_k$.

Kolmogorov's extension theorem guarantees that any *consistent* family of finite-dimensional distributions (marginalizing the $(k+1)$-th order distribution over one variable gives the $k$-th order one, and permuting the times permutes the variables) defines a random process. In practice we almost never know all of these distributions. Two situations are tractable:

1. **The process is built explicitly** from a few random variables, e.g. $X(t)=A\cos(\omega_0t+\Theta)$. Then all distributions follow from those of $A$ and $\Theta$.
2. **The process is Gaussian** (§21.5). Then the first two moments determine everything.

Otherwise we settle for the first two moments, which is what the rest of this section introduces.

### Ensemble averages: mean, autocorrelation, autocovariance

> [!definition] Mean, autocorrelation, autocovariance
> $$
> \begin{aligned}
> \mu_X(t) &= E[X(t)]=\int_{-\infty}^{\infty}x\,f_X(x;t)\,dx \\
> R_{XX}(t_1,t_2) &= E[X(t_1)X(t_2)]=\iint x_1x_2\,f_X(x_1,x_2;t_1,t_2)\,dx_1\,dx_2 \\
> C_{XX}(t_1,t_2) &= E\big[(X(t_1)-\mu_X(t_1))(X(t_2)-\mu_X(t_2))\big] = R_{XX}(t_1,t_2)-\mu_X(t_1)\mu_X(t_2)
> \end{aligned}
> $$
> The **variance function** is $\sigma_X^2(t)=C_{XX}(t,t)$ and the **mean-square value** (average power at time $t$, in a $1\ \Omega$ normalization) is $E[X^2(t)]=R_{XX}(t,t)$. The **correlation coefficient** is $\rho_{XX}(t_1,t_2)=C_{XX}(t_1,t_2)/(\sigma_X(t_1)\sigma_X(t_2))\in[-1,1]$.
> For a complex process, $R_{XX}(t_1,t_2)=E[X(t_1)X^*(t_2)]$ and $C_{XX}(t_1,t_2)=R_{XX}(t_1,t_2)-\mu_X(t_1)\mu_X^*(t_2)$.

These are **ensemble averages**: at each time they average over all outcomes $\zeta$, i.e. vertically across the ensemble of traces. The mean function $\mu_X(t)$ is a deterministic function of time (for example, a DC offset that drifts as an amplifier warms up). The autocorrelation is a deterministic function of *two* times: it measures how strongly the value at $t_1$ is (linearly) related to the value at $t_2$.

Two facts follow immediately from the definition and are used constantly:

- **Symmetry:** $R_{XX}(t_1,t_2)=R_{XX}(t_2,t_1)$ for a real process (and $R_{XX}(t_1,t_2)=R_{XX}^*(t_2,t_1)$ for a complex one).
- **Non-negative definiteness:** for any times $t_1,\dots,t_n$ and real constants $a_1,\dots,a_n$,
$$
\sum_{i=1}^n\sum_{k=1}^n a_ia_k\,R_{XX}(t_i,t_k)=E\Big[\Big(\sum_i a_iX(t_i)\Big)^2\Big]\ge0 .
$$
So the matrix $[R_{XX}(t_i,t_k)]$ is a correlation matrix (Ch. 19) for every choice of times. Not every symmetric function of two variables can be an autocorrelation.

For two processes on the same probability space we also define the **cross-correlation** $R_{XY}(t_1,t_2)=E[X(t_1)Y(t_2)]$ and **cross-covariance** $C_{XY}(t_1,t_2)=R_{XY}(t_1,t_2)-\mu_X(t_1)\mu_Y(t_2)$; they are studied in §21.3.

### Time averages

For a single sample function $x(t)$ we can also average *along time*:

> [!definition] Time averages of a sample function
> $$
> \langle x\rangle=\lim_{T\to\infty}\frac{1}{2T}\int_{-T}^{T}x(t)\,dt,\qquad
> \mathcal{R}_x(\tau)=\langle x(t+\tau)x(t)\rangle=\lim_{T\to\infty}\frac{1}{2T}\int_{-T}^{T}x(t+\tau)x(t)\,dt .
> $$
> For a sequence, $\langle x\rangle=\lim_{N\to\infty}\frac{1}{2N+1}\sum_{n=-N}^{N}x[n]$, and similarly for $\mathcal{R}_x[k]$.

The two kinds of average are conceptually very different:

| | Ensemble average $E[\cdot]$ | Time average $\langle\cdot\rangle$ |
|---|---|---|
| averages over | outcomes $\zeta$ at fixed $t$ | time $t$ for fixed $\zeta$ |
| result | deterministic function of $t$ | a number (or a function of $\tau$) that may depend on $\zeta$, i.e. a random variable |
| what you need | many independent records | one long record |
| what a DC voltmeter reads | — | $\langle x\rangle$ (over its averaging time) |
| what a true-RMS meter reads | — | $\sqrt{\langle x^2\rangle}$ |

Laboratory instruments compute time averages; theory is most easily done with ensemble averages. A process for which the two coincide (with probability 1) is called **ergodic** — §21.7 makes this precise and shows that it is a genuine assumption, not a tautology.

### A catalogue of standard processes

| Process | Construction | Typical EE use |
|---|---|---|
| Random-phase sinusoid | $X(t)=A\cos(\omega_0 t+\Theta)$, $\Theta\sim U[0,2\pi)$ | oscillator with unknown phase; carrier at a receiver |
| i.i.d. sequence | $X[n]$ independent, identically distributed | data symbols; ideal discrete white noise |
| Bernoulli / binary sequence | $X[n]\in\{-1,+1\}$ i.i.d., equally likely | random data bits, PRBS test patterns |
| Random walk | $S[n]=\sum_{k=1}^{n}X[k]$, $X[k]$ i.i.d. $\pm1$ | phase-noise accumulation, drift, integrated error |
| Moving average MA($q$) | $X[n]=\sum_{k=0}^{q}b_kW[n-k]$, $W$ i.i.d. | FIR-filtered noise |
| Autoregressive AR(1) | $X[n]=aX[n-1]+W[n]$ | IIR-filtered noise, sensor bias, channel fading (§21.4) |
| Random binary waveform | $X(t)=\sum_k a_k\,p(t-kT-D)$, $a_k=\pm A$, $D\sim U[0,T)$ | NRZ data on a cable |
| Random telegraph signal | $X(t)=\pm1$, sign flips at Poisson times of rate $\lambda$ | burst (RTS) noise in MOSFETs, on/off traffic |
| Poisson counting process | $N(t)$ = number of events in $(0,t]$, rate $\lambda$ | photon and electron arrivals, shot noise |
| Wiener process | limit of a scaled random walk (§21.5) | Brownian motion, oscillator phase noise |
| Gaussian process | all finite-dimensional distributions Gaussian | thermal noise, aggregate interference |

A useful structural property shared by the random walk, the Poisson process and the Wiener process is:

> [!definition] Independent increments
> $X(t)$ has **independent increments** if for any $t_1<t_2<\dots<t_k$ the increments $X(t_2)-X(t_1),\ X(t_3)-X(t_2),\dots,X(t_k)-X(t_{k-1})$ are mutually independent. It has **stationary increments** if the distribution of $X(t+s)-X(t)$ depends only on $s$.

For a process with independent increments starting from $X(0)=0$, the autocovariance has the universal form $C_{XX}(t_1,t_2)=\sigma_X^2(\min(t_1,t_2))$: for $t_1\le t_2$, write $X(t_2)=X(t_1)+[X(t_2)-X(t_1)]$; the bracket is independent of $X(t_1)$, so $\operatorname{Cov}(X(t_1),X(t_2))=\operatorname{Var}X(t_1)$.

### Worked examples

> [!example] Worked Example 21.1 — Mean and autocorrelation of a random-phase sinusoid
> Let $X(t)=A\cos(\omega_0 t+\Theta)$ with constants $A,\omega_0$ and $\Theta\sim U[0,2\pi)$.
>
> **Step 1 — mean.** The pdf of $\Theta$ is $1/(2\pi)$ on $[0,2\pi)$, so
> $$
> \mu_X(t)=\frac{1}{2\pi}\int_0^{2\pi}A\cos(\omega_0 t+\theta)\,d\theta = 0\quad\text{for every }t,
> $$
> because a cosine integrates to zero over a full period of its argument.
>
> **Step 2 — autocorrelation.** Using $\cos a\cos b=\tfrac12[\cos(a-b)+\cos(a+b)]$,
> $$
> R_{XX}(t_1,t_2) = \frac{A^2}{2}E\big[\cos\omega_0(t_1-t_2)\big] + \frac{A^2}{2}E\big[\cos(\omega_0(t_1+t_2)+2\Theta)\big] = \frac{A^2}{2}\cos\omega_0(t_1-t_2).
> $$
> The first expectation is of a constant; the second is zero because $2\Theta$ is uniform over two full periods.
>
> **Step 3 — interpretation.** $R_{XX}$ depends only on $t_1-t_2$, not on absolute time — the hallmark of stationarity (§21.2). The power $R_{XX}(t,t)=A^2/2$ is the familiar mean-square value of a sinusoid of amplitude $A$.
>
> **Step 4 — first-order pdf.** Since $\omega_0t+\Theta$ (mod $2\pi$) is uniform, $X(t)$ has the arcsine pdf $f_X(x)=1/(\pi\sqrt{A^2-x^2})$ for $\lvert x\rvert<A$, the same for every $t$. A numerical check: $\int_{-A}^{A}x^2/(\pi\sqrt{A^2-x^2})\,dx=A^2/2$. ✓
>
> **Step 5 — time averages.** Take one sample function, $x(t)=A\cos(\omega_0t+\theta_1)$. Its time average is $0$ and its time autocorrelation is $\langle x(t+\tau)x(t)\rangle=\frac{A^2}{2}\cos\omega_0\tau$ — the same as the ensemble results, whatever $\theta_1$ is. This process is ergodic (§21.7).

> [!example] Worked Example 21.2 — Random walk: variance that grows
> Let $S[0]=0$ and $S[n]=\sum_{k=1}^{n}X[k]$ for $n\ge1$, with $X[k]$ i.i.d., $P(X=\pm1)=\tfrac12$, so $E[X]=0$ and $E[X^2]=1$.
>
> **Mean.** $\mu_S[n]=\sum_{k=1}^n E[X[k]]=0$.
>
> **Autocorrelation.** $R_{SS}[n_1,n_2]=E\Big[\sum_{k=1}^{n_1}\sum_{l=1}^{n_2}X[k]X[l]\Big]=\sum_k\sum_l E[X[k]X[l]]$. By independence $E[X[k]X[l]]=\delta[k-l]$, so only the terms with $k=l\le\min(n_1,n_2)$ survive:
> $$
> R_{SS}[n_1,n_2]=\min(n_1,n_2).
> $$
> This agrees with the independent-increments formula $C(t_1,t_2)=\sigma^2(\min(t_1,t_2))$.
>
> **Consequences.** $\sigma_S^2[n]=n$: the ensemble spreads like $\sqrt n$. After $n=10^4$ steps the typical excursion is about $100$, and by the central limit theorem $P\{\lvert S[10^4]\rvert>200\}\approx 2Q(2)\approx0.046$. The correlation coefficient between $S[n]$ and $S[4n]$ is $n/\sqrt{n\cdot4n}=1/2$ for every $n$: the walk "remembers" its past indefinitely. Integrated noise — e.g. oscillator phase driven by white frequency noise, or the position estimate of an inertial navigator integrating noisy accelerometer readings — behaves the same way.

> [!example] Worked Example 21.3 — The Poisson counting process
> $N(t)$, $t\ge0$, counts events in $(0,t]$; it has independent, stationary increments and $N(t+s)-N(t)\sim\text{Poisson}(\lambda s)$. Hence $\mu_N(t)=\lambda t$ and $\sigma_N^2(t)=\lambda t$.
>
> **Autocorrelation.** For $t_1\le t_2$, by the independent-increments argument $C_{NN}(t_1,t_2)=\operatorname{Var}N(t_1)=\lambda t_1$. So in general
> $$
> C_{NN}(t_1,t_2)=\lambda\min(t_1,t_2),\qquad R_{NN}(t_1,t_2)=\lambda\min(t_1,t_2)+\lambda^2t_1t_2 .
> $$
>
> **Numbers.** A photodiode receives photons at $\lambda=10^{9}\ \text{s}^{-1}$. In a $1\ \mu$s window the expected count is $\lambda t=1000$ with standard deviation $\sqrt{1000}=31.6$, a relative fluctuation of $3.2\%$. This $1/\sqrt{\lambda t}$ fluctuation is **shot noise**; the corresponding photocurrent noise has one-sided PSD $2qI$ (§21.8). Note the mean grows with $t$, so $N(t)$ is not stationary — but its *increments over fixed-length windows* are.

> [!example] Worked Example 21.4 — Sampling a process gives a random sequence
> A continuous-time process $X(t)$ with $\mu_X=0$ and $R_{XX}(t_1,t_2)=\sigma^2e^{-\alpha\lvert t_1-t_2\rvert}$ is sampled every $T_s$ seconds: $X[n]=X(nT_s)$. Then
> $$
> R_{XX}[n_1,n_2]=E[X(n_1T_s)X(n_2T_s)]=\sigma^2e^{-\alpha T_s\lvert n_1-n_2\rvert}=\sigma^2a^{\lvert n_1-n_2\rvert},\qquad a=e^{-\alpha T_s}.
> $$
> With $\alpha=1000\ \text{s}^{-1}$ (correlation time 1 ms) and $T_s=100\ \mu$s, $a=e^{-0.1}=0.905$: adjacent samples are 90% correlated, and samples are essentially uncorrelated ($\rho<0.05$) only when $\lvert n_1-n_2\rvert\ge30$. Sampling at $T_s=5$ ms instead gives $a=e^{-5}=0.0067$: practically uncorrelated samples. This geometric autocorrelation is exactly that of the AR(1) sequence we study in §21.4 — sampling a first-order Gauss–Markov process produces an AR(1) sequence.

> [!warning] Common Pitfall — confusing the two averages
> $E[X(t)]$ is a function of $t$ computed across the ensemble; $\langle x\rangle$ is a number computed along one record. Writing "the mean of the noise" without saying which is ambiguous. For a process like $X(t)=A$ (a random constant, see §21.7), $E[X(t)]=E[A]$ but $\langle x\rangle=A$: they differ for every outcome except those with $A=E[A]$.

> [!tip] Engineering Note
> Measurement noise, quantization error, channel noise, phase noise and interference are all modelled as random processes. What we need for design is rarely the full set of finite-dimensional distributions; it is the mean (DC offset), the power (RMS value) and *how that power is spread over time lag and frequency*. The autocorrelation (§21.3) and the power spectral density (§21.4) carry exactly this information, and linear systems transform them by simple rules.

## 21.2 Stationarity

### Motivation

Most of the convenient machinery of signal analysis — transfer functions, spectra, filters designed once and used forever — rests on **time invariance**. For random signals the corresponding idea is that the *statistics* of the process do not change with time: the noise from a resistor in thermal equilibrium looks statistically the same today as tomorrow, even though no two sample functions agree. Stationarity is what makes the autocorrelation a function of one variable (the lag) instead of two, and that in turn is what allows the Fourier transform to enter (§21.4). Without it, the notion of a "noise spectrum" is not even well defined.

### Strict-sense stationarity

> [!definition] Strict-sense stationarity (SSS)
> $X(t)$ is **strict-sense stationary** if all its finite-dimensional distributions are invariant under a time shift: for every $k$, every set of times $t_1,\dots,t_k$ and every shift $c$,
> $$
> F_X(x_1,\dots,x_k;\,t_1+c,\dots,t_k+c)=F_X(x_1,\dots,x_k;\,t_1,\dots,t_k).
> $$
> If this holds only for $k\le N$, the process is **stationary of order $N$**.

Consequences for an SSS process:

- **First order** ($k=1$): $f_X(x;t)=f_X(x)$ does not depend on $t$. In particular $\mu_X(t)=\mu_X$ and $\sigma_X^2(t)=\sigma_X^2$ are constants.
- **Second order** ($k=2$): $f_X(x_1,x_2;t_1,t_2)$ depends on the times only through $\tau=t_1-t_2$ (choose $c=-t_2$). Hence $R_{XX}(t_1,t_2)$ is a function of $t_1-t_2$ alone.

> [!theorem] i.i.d. sequences are strictly stationary
> If $X[n]$ are independent and identically distributed with CDF $F$, then $X[n]$ is SSS.

*Proof.* For any distinct $n_1,\dots,n_k$ and shift $c$, independence gives $F_X(x_1,\dots,x_k;n_1+c,\dots,n_k+c)=\prod_iF(x_i)$, which does not involve $c$. $\square$

A memoryless function of an SSS process, $Y(t)=g(X(t))$, is again SSS (the joint distribution of $Y(t_1+c),\dots,Y(t_k+c)$ is determined by that of $X(t_1+c),\dots,X(t_k+c)$). So is the output of any time-invariant (possibly nonlinear) system with memory driven by an SSS input, once start-up transients have died away.

SSS is a strong requirement and is hard to verify from data: it constrains distributions of every order. For linear-system analysis a much weaker condition, involving only the first two moments, suffices.

### Wide-sense stationarity

> [!definition] Wide-sense stationarity (WSS)
> $X(t)$ is **wide-sense stationary** (weakly stationary) if
> 1. $E[X^2(t)]<\infty$ for all $t$,
> 2. $\mu_X(t)=\mu_X$ is constant, and
> 3. $R_{XX}(t_1,t_2)$ depends only on $\tau=t_1-t_2$.
>
> We then write, with a slight abuse of notation,
> $$
> R_{XX}(\tau)=E[X(t+\tau)X(t)],\qquad C_{XX}(\tau)=R_{XX}(\tau)-\mu_X^2,\qquad \sigma_X^2=C_{XX}(0).
> $$
> For a random sequence: $E[X[n]]=\mu_X$ and $R_{XX}[k]=E[X[n+k]X[n]]$ independent of $n$.

> [!definition] Joint wide-sense stationarity
> $X(t)$ and $Y(t)$ are **jointly WSS** if each is WSS and $R_{XY}(t_1,t_2)=E[X(t_1)Y(t_2)]$ depends only on $t_1-t_2$. We write $R_{XY}(\tau)=E[X(t+\tau)Y(t)]$.

Relations between the notions:

- **SSS with finite second moment $\Rightarrow$ WSS.** Shown above.
- **WSS $\not\Rightarrow$ SSS** in general (Example 21.8 below).
- **For Gaussian processes, WSS $\Leftrightarrow$ SSS** (§21.5), because the mean and autocovariance determine all finite-dimensional distributions.
- Each of $X$ and $Y$ being WSS does **not** make them jointly WSS. For example, with $X$ WSS and $Y(t)=X(-t)$ (also WSS), $R_{XY}(t_1,t_2)=R_{XX}(t_1+t_2)$, which depends on the sum of the times.

> [!tip] Engineering Note — why "WSS" is the working assumption
> Every result about filters and spectra in this chapter needs only WSS. Since estimating higher-order distributions takes enormous amounts of data, practising engineers almost always mean "WSS" when they say "stationary". When a noise source is additionally Gaussian — the usual case for thermal noise and many aggregate interference sources — nothing is lost.

### Testing for stationarity: worked examples

The procedure is always the same: compute $\mu_X(t)$ and $R_{XX}(t+\tau,t)$ and check whether $t$ drops out.

> [!example] Worked Example 21.5 — Random phase, random amplitude, random pair
> **(a)** $X(t)=A\cos(\omega_0 t+\Theta)$, $\Theta\sim U[0,2\pi)$: from Example 21.1, $\mu_X=0$ and $R_{XX}(t+\tau,t)=\tfrac{A^2}{2}\cos\omega_0\tau$. **WSS.** (In fact SSS: shifting time by $c$ is the same as shifting $\Theta$ by $\omega_0c$, which leaves its uniform distribution modulo $2\pi$ unchanged.)
>
> **(b)** $Y(t)=B\cos\omega_0 t$ with $B$ a random variable, $E[B]=\mu_B$, $E[B^2]=m_2$:
> $$
> \mu_Y(t)=\mu_B\cos\omega_0 t,\qquad R_{YY}(t+\tau,t)=m_2\cos\omega_0(t+\tau)\cos\omega_0t=\frac{m_2}{2}\big[\cos\omega_0\tau+\cos\omega_0(2t+\tau)\big].
> $$
> If $\mu_B\neq0$ the mean depends on $t$; even with $\mu_B=0$, $R_{YY}$ depends on $t$ through the second term. **Not WSS.** Physically: every sample function passes through zero at the same instants $t=(k+\tfrac12)\pi/\omega_0$, so the ensemble power $R_{YY}(t,t)=m_2\cos^2\omega_0t$ is zero there.
>
> **(c)** $Z(t)=B\cos\omega_0t + C\sin\omega_0t$ with $B,C$ zero mean, uncorrelated, both of variance $\sigma^2$:
> $$
> \begin{aligned}
> \mu_Z(t)&=0,\\
> R_{ZZ}(t_1,t_2)&=E[B^2]\cos\omega_0t_1\cos\omega_0t_2+E[C^2]\sin\omega_0t_1\sin\omega_0t_2+E[BC]\sin\omega_0(t_1+t_2)\\
> &=\sigma^2\cos\omega_0(t_1-t_2).
> \end{aligned}
> $$
> **WSS.** The random phase is supplied by the random pair $(B,C)$. If $E[B^2]\ne E[C^2]$ or $E[BC]\ne0$, the $t$-dependence returns — both conditions are necessary. $Z$ is WSS for *any* distribution of $(B,C)$ with these moments, but it is SSS only if $(B,C)$ is circularly symmetric (e.g. jointly Gaussian with equal variances).

> [!example] Worked Example 21.6 — A phase that is not uniform enough
> Let $X(t)=A\cos(\omega_0t+\Theta)$ with $\Theta\sim U[0,\pi]$ (only half a cycle). Then
> $$
> \mu_X(t)=\frac{A}{\pi}\int_0^{\pi}\cos(\omega_0t+\theta)\,d\theta=\frac{A}{\pi}\big[\sin(\omega_0t+\pi)-\sin\omega_0t\big]=-\frac{2A}{\pi}\sin\omega_0t .
> $$
> The mean oscillates with amplitude $2A/\pi\approx0.637A$, so the process is **not WSS**. Curiously, the autocorrelation is still $\tfrac{A^2}{2}\cos\omega_0(t_1-t_2)$: the term $E[\cos(\omega_0(t_1+t_2)+2\Theta)]$ vanishes because $2\Theta$ *is* uniform on $[0,2\pi)$. Conditions 2 and 3 of the WSS definition are independent, and both must be checked. (The autocovariance, $R_{XX}-\mu_X(t_1)\mu_X(t_2)$, does depend on both times.)

> [!example] Worked Example 21.7 — Random binary (NRZ) waveform
> Let $X(t)=\sum_{k=-\infty}^{\infty}a_k\,p(t-kT-D)$, where $p(t)=1$ for $0\le t<T$ and $0$ otherwise, the $a_k$ are i.i.d. with $P(a_k=\pm A)=\tfrac12$, and the timing offset $D\sim U[0,T)$ is independent of the $a_k$. This models a data stream whose clock phase is unknown to the observer.
>
> **Step 1 — mean.** $E[X(t)]=\sum_kE[a_k]\,E[p(t-kT-D)]=0$.
>
> **Step 2 — the product.** $X(t)$ and $X(t+\tau)$ are each $\pm A$. If $t$ and $t+\tau$ fall in the **same** bit interval, $X(t+\tau)X(t)=a_k^2=A^2$. If they fall in **different** intervals, the product is $a_ka_l$ with $k\ne l$, which has zero mean. So
> $$
> R_{XX}(t+\tau,t)=A^2\,P\{t\text{ and }t+\tau\text{ lie in the same bit}\}.
> $$
>
> **Step 3 — the probability.** For $\lvert\tau\rvert\ge T$ it is $0$. For $0\le\tau<T$: the bit boundaries are at $kT+D$, and since $D$ is uniform, the position of $t$ within its bit, $U=(t-D)\bmod T$, is uniform on $[0,T)$. Both instants lie in the same bit iff $U+\tau<T$, which has probability $(T-\tau)/T$. By symmetry the same holds for $\tau<0$ with $\lvert\tau\rvert$.
>
> **Result.**
> $$
> R_{XX}(\tau)=\begin{cases}A^2\left(1-\dfrac{\lvert\tau\rvert}{T}\right),&\lvert\tau\rvert<T,\\[2mm]0,&\lvert\tau\rvert\ge T.\end{cases}
> $$
> The process is **WSS** with a triangular autocorrelation. Without the random offset $D$ (a fixed, known clock) the probability in Step 3 depends on where $t$ sits within its bit, and $R_{XX}(t+\tau,t)$ is periodic in $t$ with period $T$ — the process is then cyclostationary (below). The triangle reappears in §21.4 as the $\operatorname{sinc}^2$ spectrum of NRZ data.

```desmos-graph
left=-2.5; right=2.5; bottom=-0.2; top=1.3
---
y=\max(1-\left|x\right|,0)|blue
y=e^{-\left|x\right|}|red|dashed
```

*Triangular autocorrelation of the random NRZ waveform ($A=1$, $T=1$, blue) and the exponential autocorrelation $e^{-2\lambda\lvert\tau\rvert}$ of a random telegraph signal with the same average switching rate $\lambda=1/(2T)$ (red, dashed). The NRZ signal is completely decorrelated after one bit; the telegraph signal, whose switching instants are not tied to a clock, decorrelates only gradually.*

> [!example] Worked Example 21.8 — WSS but not strictly stationary
> Build a sequence of independent variables with $X[n]\sim\mathcal{N}(0,1)$ for even $n$ and $X[n]=\pm1$ with probability $\tfrac12$ each for odd $n$. Then $E[X[n]]=0$ and $E[X^2[n]]=1$ for all $n$, and $E[X[n+k]X[n]]=0$ for $k\ne0$ by independence, so
> $$
> \mu_X=0,\qquad R_{XX}[k]=\delta[k].
> $$
> The sequence is **WSS** (indeed white, §21.6). But the first-order distribution alternates between a Gaussian and two point masses, so $P\{X[n]>1.5\}$ equals $Q(1.5)\approx0.067$ for even $n$ and $0$ for odd $n$: **not SSS**. A second-order description cannot see this difference — a reminder that "white" and "i.i.d." are not the same thing.

### Cyclostationarity

> [!definition] Wide-sense cyclostationarity
> $X(t)$ is **wide-sense cyclostationary** with period $T_0$ if
> $$
> \mu_X(t+T_0)=\mu_X(t),\qquad R_{XX}(t_1+T_0,t_2+T_0)=R_{XX}(t_1,t_2)\quad\text{for all }t,t_1,t_2 .
> $$

Digitally modulated signals $X(t)=\sum_k a_k\,p(t-kT_b)$ with a fixed clock, sampled-and-held noise, and noise in switched-capacitor and mixer circuits (whose operating point is periodically pumped by a clock or local oscillator) are all cyclostationary. The standard way of bringing them into the WSS framework is to randomize the time origin.

> [!theorem] Stationarizing a cyclostationary process
> If $X(t)$ is wide-sense cyclostationary with period $T_0$ and $D\sim U[0,T_0)$ is independent of $X$, then $Y(t)=X(t-D)$ is WSS with
> $$
> \mu_Y=\frac{1}{T_0}\int_0^{T_0}\mu_X(t)\,dt,\qquad R_{YY}(\tau)=\bar R_{XX}(\tau)=\frac{1}{T_0}\int_0^{T_0}R_{XX}(t+\tau,t)\,dt .
> $$

*Proof.* Condition on $D=d$: $E[Y(t+\tau)Y(t)\mid D=d]=R_{XX}(t+\tau-d,\,t-d)$. Averaging over $d$ uniform on one period gives $\frac{1}{T_0}\int_{t-T_0}^{t}R_{XX}(s+\tau,s)\,ds$ with $s=t-d$. Because $s\mapsto R_{XX}(s+\tau,s)$ is $T_0$-periodic, the integral over any interval of length $T_0$ is the same, so the result does not depend on $t$. The mean is handled identically. $\square$

The **period-averaged autocorrelation** $\bar R_{XX}(\tau)$ is what a spectrum analyzer effectively measures for a cyclostationary signal, because the analyzer's sweep is not synchronized to the data clock. Example 21.7 is exactly this construction applied to NRZ data.

### Non-stationary processes in practice

Many processes are not stationary but are close enough for practical purposes:

- **Asymptotically stationary:** a stable filter switched on at $t=0$ with a WSS input produces an output whose statistics converge to stationary values as the transient decays (the AR(1) sequence started from $X[0]=0$ in §21.4 is an example).
- **Locally stationary:** speech is commonly treated as stationary over 20–30 ms frames; a radio channel's fading statistics are stationary over distances of tens of wavelengths.
- **Stationary increments:** the random walk, Poisson process and Wiener process are not stationary, but their increments over fixed lags are. Oscillator phase noise and clock-jitter accumulation are analyzed this way (Allan variance).

> [!warning] Common Pitfall — stationarity belongs to the ensemble
> Stationarity is a property of the *ensemble*, not of one record. A single sample function of a stationary process can look "non-stationary" over a short window (a burst of large values happens by chance), and a non-stationary process such as slow thermal drift can look stationary over a short window. Before invoking stationarity, state the time scale over which the assumption is meant to hold, and check that the record is long compared with the correlation time but short compared with the time scale of drift.

## 21.3 Autocorrelation and cross-correlation

### Motivation

The autocorrelation function answers a simple question: *if I know the value of the signal now, how much do I know about its value $\tau$ seconds from now?* For thermal noise in a wideband amplifier the answer is "almost nothing after a few nanoseconds"; for the slowly wandering bias of a MEMS gyroscope, "a great deal, even after minutes". The **cross-correlation** asks the same question about two different signals, and it is the workhorse of synchronization, ranging, system identification and interference analysis: whenever an engineer needs to know whether, and with what delay, one signal is contained in another, a correlator is involved.

Throughout this section $X(t)$ and $Y(t)$ are real and (jointly) WSS unless stated otherwise.

### Properties of the autocorrelation function

> [!theorem] Properties of the autocorrelation function
> 1. **Power:** $R_{XX}(0)=E[X^2(t)]\ge0$ is the average power (mean-square value).
> 2. **Even symmetry:** $R_{XX}(-\tau)=R_{XX}(\tau)$.
> 3. **Peak at the origin:** $\lvert R_{XX}(\tau)\rvert\le R_{XX}(0)$ for all $\tau$.
> 4. **Non-negative definiteness:** for all $n$, times $t_1,\dots,t_n$ and real $a_1,\dots,a_n$,
> $$
> \sum_{i=1}^n\sum_{k=1}^n a_ia_k\,R_{XX}(t_i-t_k)\ge0 .
> $$
> Equivalently (Bochner's theorem, §21.4), its Fourier transform is non-negative.
> 5. **Periodicity:** if $R_{XX}(T_0)=R_{XX}(0)$ for some $T_0\neq0$, then $R_{XX}$ is periodic with period $T_0$, and $X(t+T_0)=X(t)$ with probability 1 for every $t$.
> 6. **DC component:** if $X$ has no periodic components and $\mu_X\ne0$, then $R_{XX}(\tau)\to\mu_X^2$ as $\lvert\tau\rvert\to\infty$; i.e. $C_{XX}(\tau)\to0$.
> 7. **Continuity:** if $R_{XX}$ is continuous at $\tau=0$, it is continuous everywhere.

*Proofs.*

1. Definition with $\tau=0$.
2. $R_{XX}(-\tau)=E[X(t-\tau)X(t)]$; substitute $s=t-\tau$ and use stationarity: $E[X(s)X(s+\tau)]=R_{XX}(\tau)$.
3. Cauchy–Schwarz for random variables (Ch. 19): $\big(E[X(t+\tau)X(t)]\big)^2\le E[X^2(t+\tau)]\,E[X^2(t)]=R_{XX}(0)^2$. Alternatively, $0\le E[(X(t+\tau)\pm X(t))^2]=2R_{XX}(0)\pm2R_{XX}(\tau)$.
4. $0\le E\big[(\sum_i a_iX(t_i))^2\big]$ expands to the double sum.
5. $E[(X(t+T_0)-X(t))^2]=2R_{XX}(0)-2R_{XX}(T_0)=0$, so $X(t+T_0)=X(t)$ with probability 1. Then $R_{XX}(\tau+T_0)=E[X(t+\tau+T_0)X(t)]=E[X(t+\tau)X(t)]=R_{XX}(\tau)$.
6. Write $X=\mu_X+\tilde X$ with $\tilde X$ zero mean; $R_{XX}=\mu_X^2+C_{XX}$, and "no periodic components" is the (informal) condition that $C_{XX}(\tau)\to0$ — values far apart become uncorrelated.
7. $\lvert R_{XX}(\tau+h)-R_{XX}(\tau)\rvert=\lvert E[(X(t+\tau+h)-X(t+\tau))X(t)]\rvert\le\sqrt{E[(X(t+h)-X(t))^2]\,R_{XX}(0)}=\sqrt{2(R_{XX}(0)-R_{XX}(h))R_{XX}(0)}$, which tends to $0$ as $h\to0$ if $R_{XX}$ is continuous at $0$. $\square$

Property 7 has a physical reading: $E[(X(t+h)-X(t))^2]=2[R_{XX}(0)-R_{XX}(h)]$. So $R_{XX}$ continuous at $0$ is exactly the statement that the process is **mean-square continuous** — its sample values do not jump on average. A large curvature of $R_{XX}$ at the origin means a rapidly fluctuating process.

**Reading an autocorrelation.** A typical autocorrelation of a zero-mean noise process is a peaked, even function:

- its **height** $R_{XX}(0)$ is the power;
- its **width** is the correlation time — how long the process "remembers";
- its **shape at the origin** tells how smooth the sample functions are: a cusp ($e^{-\alpha\lvert\tau\rvert}$) means continuous but non-differentiable paths; a rounded top ($e^{-\alpha\tau^2}$, sinc) means smooth paths;
- a **constant floor** $\mu_X^2$ at large lags is DC; a **persistent oscillation** is a sinusoidal component.

### Correlation coefficient and correlation time

The normalized autocovariance

$$
\rho_{XX}(\tau)=\frac{C_{XX}(\tau)}{C_{XX}(0)}\in[-1,1]
$$

is the **correlation coefficient function**. Several definitions of a single **correlation time** are in use; the two most common are

$$
\tau_c=\int_0^\infty \rho_{XX}(\tau)\,d\tau\qquad\text{or}\qquad \tau_e:\ \rho_{XX}(\tau_e)=e^{-1}.
$$

For $\rho=e^{-\alpha\lvert\tau\rvert}$ both give $1/\alpha$. Samples spaced much farther apart than $\tau_c$ are nearly uncorrelated; this is what makes long averages effective (§21.7).

### Worked examples on autocorrelation

> [!example] Worked Example 21.9 — Random telegraph signal
> $X(t)$ takes values $\pm1$ and switches sign at the event times of a Poisson process of rate $\lambda$; $P(X(0)=1)=\tfrac12$, independent of the switching times.
>
> **Step 1 — mean.** By symmetry $P(X(t)=1)=\tfrac12$ for all $t$, so $\mu_X=0$.
>
> **Step 2 — the product.** $X(t+\tau)X(t)=+1$ if an even number of switches occurs in the interval between $t$ and $t+\tau$ and $-1$ if odd. The number of switches $N$ in an interval of length $\lvert\tau\rvert$ is Poisson with mean $\lambda\lvert\tau\rvert$.
>
> **Step 3 — sum the series.**
> $$
> R_{XX}(\tau)=P(N\text{ even})-P(N\text{ odd})=\sum_{k=0}^{\infty}(-1)^k\frac{(\lambda\lvert\tau\rvert)^k}{k!}e^{-\lambda\lvert\tau\rvert}=e^{-\lambda\lvert\tau\rvert}e^{-\lambda\lvert\tau\rvert}=e^{-2\lambda\lvert\tau\rvert}.
> $$
>
> **Interpretation.** Power $R_{XX}(0)=1$ (obvious, since $X^2=1$), correlation time $\tau_c=1/(2\lambda)$. With $\lambda=10^4$ switches per second (a typical MOSFET random-telegraph-noise trap), $\tau_c=50\ \mu$s; values 200 μs apart have correlation $e^{-4}=0.018$. The cusp at $\tau=0$ reflects the jumps of the sample paths.

```desmos-graph
left=-4; right=4; bottom=-1.2; top=1.3
---
y=e^{-\left|x\right|}|blue
y=\cos(2x)|red|dashed
y=\frac{\sin(\pi x)}{\pi x}|green
y=e^{-x^2}|purple|dotted
```

*Four valid autocorrelation shapes: exponential (telegraph or RC-filtered noise, blue — note the cusp), cosine (random-phase sinusoid, red — never decays because the process is periodic), sinc (band-limited white noise, green — negative lobes are allowed) and Gaussian (smooth, differentiable noise, purple). All are even and peak at $\tau=0$.*

> [!example] Worked Example 21.10 — Is this a valid autocorrelation?
> Consider the rectangular function $g(\tau)=1$ for $\lvert\tau\rvert\le1$ and $0$ otherwise. It is even, non-negative and peaks at $0$ — properties 1–3 hold. Is it an autocorrelation?
>
> **Test property 4** with times $t_1=0,\ t_2=1,\ t_3=2$. The matrix $[g(t_i-t_k)]$ is
> $$
> \mathbf{G}=\begin{bmatrix}1&1&0\\1&1&1\\0&1&1\end{bmatrix},\qquad \mathbf{a}=\begin{bmatrix}1\\-1\\1\end{bmatrix}:\quad \mathbf{a}^\mathsf{T}\mathbf{G}\mathbf{a}=3+2(-1-1+0)=-1<0 .
> $$
> A "variance" of $-1$ is impossible, so **$g$ is not an autocorrelation**. (The eigenvalues of $\mathbf{G}$ are $1$ and $1\pm\sqrt2$; one is negative.) In the frequency domain the same conclusion is immediate: the Fourier transform of $g$ is $2\sin\omega/\omega$, which goes negative, violating $S_{XX}\ge0$ (§21.4).
>
> By contrast the triangle $1-\lvert\tau\rvert$ (Example 21.7) has transform $\sin^2(\omega/2)/(\omega/2)^2\ge0$ and is valid. Being even and peaked is necessary but nowhere near sufficient.

> [!example] Worked Example 21.11 — Moving-average (MA(1)) sequence
> Let $W[n]$ be i.i.d. with zero mean and variance $\sigma_W^2$, and define the two-tap FIR output $X[n]=W[n]+bW[n-1]$.
>
> **Step 1 — mean:** $E[X[n]]=0$.
>
> **Step 2 — autocorrelation.** $R_{XX}[k]=E\big[(W[n+k]+bW[n+k-1])(W[n]+bW[n-1])\big]$. Using $E[W[i]W[l]]=\sigma_W^2\delta[i-l]$, the four cross terms survive only when their indices match:
> $$
> R_{XX}[0]=\sigma_W^2(1+b^2),\qquad R_{XX}[\pm1]=b\,\sigma_W^2,\qquad R_{XX}[k]=0\ \ (\lvert k\rvert\ge2).
> $$
>
> **Step 3 — numbers.** With $b=0.5$, $\sigma_W^2=1$: $R_{XX}[0]=1.25$, $R_{XX}[\pm1]=0.5$, $\rho_{XX}[1]=0.4$. Notice $\lvert\rho_{XX}[1]\rvert=\lvert b\rvert/(1+b^2)\le\tfrac12$ for every $b$: an MA(1) sequence can never have lag-1 correlation exceeding $0.5$ in magnitude. If a measured sequence has $\hat\rho[1]=0.8$ and $\hat\rho[k]\approx0$ for $k\ge2$, it cannot be MA(1) — such constraints are exactly the non-negative-definiteness of property 4.
>
> **Generalization.** For an FIR filter $X[n]=\sum_{m=0}^{q}b_mW[n-m]$ driven by white noise, $R_{XX}[k]=\sigma_W^2\sum_m b_mb_{m+k}$ (the deterministic autocorrelation of the coefficient sequence), which vanishes for $\lvert k\rvert>q$.

### Derivatives of a WSS process

Many physical quantities are derivatives of others (current in a capacitor, velocity from position), so it is useful to know the correlation of $X'(t)$. If $R_{XX}(\tau)$ is twice differentiable at $\tau=0$, then $X(t)$ is **mean-square differentiable** and, interchanging derivative and expectation,

$$
R_{X'X}(\tau)=E[X'(t+\tau)X(t)]=\frac{d}{d\tau}R_{XX}(\tau),\qquad
R_{X'X'}(\tau)=E[X'(t+\tau)X'(t)]=-\frac{d^2}{d\tau^2}R_{XX}(\tau).
$$

*Derivation.* $R_{X'X}(\tau)=\frac{\partial}{\partial t_1}E[X(t_1)X(t_2)]\big|_{t_1=t+\tau,t_2=t}=R_{XX}'(\tau)$. Then $R_{X'X'}(\tau)=\frac{\partial}{\partial t_2}R_{X'X}(t_1-t_2)=-R_{XX}''(\tau)$.

Two consequences: since $R_{XX}$ is even, $R'_{XX}(0)=0$ and so $E[X'(t)X(t)]=0$ — **a WSS process and its derivative are uncorrelated at the same instant**. And the derivative's power is $E[X'^2]=-R_{XX}''(0)$, which is finite only if $R_{XX}$ is smooth at the origin.

> [!example] Worked Example 21.12 — Smooth versus rough noise
> **(a)** $R_{XX}(\tau)=\sigma^2e^{-\alpha\tau^2}$. Then $R'_{XX}=-2\alpha\tau\sigma^2e^{-\alpha\tau^2}$ and
> $$
> R_{X'X'}(\tau)=-R''_{XX}(\tau)=2\alpha\sigma^2(1-2\alpha\tau^2)e^{-\alpha\tau^2},\qquad E[X'^2]=2\alpha\sigma^2 .
> $$
> With $\sigma=1$ mV and $\alpha=10^{8}\ \text{s}^{-2}$, the RMS slope is $\sqrt{2\times10^8}\times10^{-3}=14.1$ V/s. This value is needed, for instance, to convert voltage noise into timing jitter at a comparator: $\sigma_t\approx\sigma_v/\lvert\text{slope of the signal}\rvert$, and the expected rate of zero crossings (Rice's formula, §21.5).
>
> **(b)** $R_{XX}(\tau)=\sigma^2e^{-\alpha\lvert\tau\rvert}$ has a cusp at $0$: $R''_{XX}(0)$ does not exist (formally, $-R''_{XX}$ contains $2\alpha\sigma^2\delta(\tau)$). The process is mean-square continuous but **not differentiable** — the derivative would have infinite power. RC-filtered white noise and the telegraph signal are of this kind; their paths are rough at every scale down to where the idealized model breaks down.

### Cross-correlation

> [!definition] Cross-correlation and cross-covariance
> For jointly WSS $X$ and $Y$:
> $$
> R_{XY}(\tau)=E[X(t+\tau)Y(t)],\qquad C_{XY}(\tau)=R_{XY}(\tau)-\mu_X\mu_Y .
> $$
> $X$ and $Y$ are **uncorrelated** if $C_{XY}(\tau)=0$ for all $\tau$, and **orthogonal** if $R_{XY}(\tau)=0$ for all $\tau$. Independent processes are uncorrelated; zero-mean uncorrelated processes are orthogonal.

The ordering convention matters: here the *first* subscript carries the time shift $+\tau$. Some texts use $E[X(t)Y(t+\tau)]$, which reverses the sign of $\tau$. Check the convention before using a published formula.

> [!theorem] Properties of the cross-correlation
> 1. $R_{XY}(\tau)=R_{YX}(-\tau)$ — **not** even in general; its peak may be anywhere.
> 2. $\lvert R_{XY}(\tau)\rvert\le\sqrt{R_{XX}(0)R_{YY}(0)}$.
> 3. $\lvert R_{XY}(\tau)\rvert\le\tfrac12\big[R_{XX}(0)+R_{YY}(0)\big]$.
> 4. If $Z=X+Y$: $R_{ZZ}(\tau)=R_{XX}(\tau)+R_{YY}(\tau)+R_{XY}(\tau)+R_{YX}(\tau)$. For uncorrelated $X,Y$ with at least one of them zero mean, $R_{ZZ}=R_{XX}+R_{YY}$: **powers of uncorrelated zero-mean processes add**.

*Proofs.* (1) Substitute $s=t+\tau$ in $E[X(t+\tau)Y(t)]=E[Y(s-\tau)X(s)]=R_{YX}(-\tau)$. (2) Cauchy–Schwarz. (3) The geometric mean is at most the arithmetic mean. (4) Expand $E[(X(t+\tau)+Y(t+\tau))(X(t)+Y(t))]$. $\square$

Property 4 is the reason noise powers, not noise voltages, are added when independent sources are combined: two independent 1 mV-rms noise sources in series give $\sqrt2$ mV rms, not 2 mV.

> [!example] Worked Example 21.13 — Signal plus independent noise
> A WSS signal $S(t)$ with $R_{SS}(\tau)=4e^{-100\lvert\tau\rvert}\ \mathrm{V^2}$ is corrupted by independent zero-mean noise $N(t)$ with $R_{NN}(\tau)=0.25e^{-10^5\lvert\tau\rvert}\ \mathrm{V^2}$. The observed process is $X=S+N$.
>
> **Autocorrelation.** By property 4 with $R_{SN}=R_{NS}=E[S]E[N]=0$:
> $$
> R_{XX}(\tau)=4e^{-100\lvert\tau\rvert}+0.25e^{-10^5\lvert\tau\rvert}.
> $$
> Total power $R_{XX}(0)=4.25\ \mathrm{V^2}$ (2.06 V rms); SNR $=4/0.25=16$ (12.0 dB).
>
> **Separating them by lag.** At $\tau=50\ \mu$s the noise term is $0.25e^{-5}=0.0017$ while the signal term is $4e^{-0.005}=3.98$. A correlator evaluated at a lag longer than the noise correlation time but shorter than the signal's sees almost pure signal: the "spike" at the origin of a measured autocorrelation is the signature of wideband noise, and its height is the noise power. This is how noise power is often estimated in practice without switching the signal off.
>
> **Cross-correlation with the clean signal.** $R_{XS}(\tau)=E[(S(t+\tau)+N(t+\tau))S(t)]=R_{SS}(\tau)$: independent noise drops out of a cross-correlation entirely (on average).

> [!example] Worked Example 21.14 — Time-delay estimation
> A transmitted signal $X(t)$ (WSS, zero mean, $R_{XX}(\tau)=e^{-1000\lvert\tau\rvert}$) arrives as $Y(t)=0.5\,X(t-D)+N(t)$ with $D=2\ \text{ms}$ and noise $N(t)$ zero mean and independent of $X$.
>
> **Step 1.**
> $$
> R_{YX}(\tau)=E[Y(t+\tau)X(t)]=0.5\,E[X(t+\tau-D)X(t)]+E[N(t+\tau)]E[X(t)]=0.5\,R_{XX}(\tau-D).
> $$
>
> **Step 2.** So $R_{YX}(\tau)=0.5\,e^{-1000\lvert\tau-0.002\rvert}$, which peaks with value $0.5$ at $\tau=D=2$ ms. By property 1, $R_{XY}(\tau)=R_{YX}(-\tau)$ peaks at $\tau=-2$ ms.
>
> **Step 3 — resolution.** The peak's half-width (where $R_{YX}$ falls to $0.5\,e^{-1}$) is $\pm1$ ms. Two echoes closer than about one correlation time of $X$ merge into one peak. Hence radar and GPS use signals with very **short correlation times** (wide bandwidth): a GPS C/A code with chip duration $T_c\approx0.98\ \mu$s has a triangular autocorrelation of half-width $T_c$, corresponding to about 300 m of range, and receivers interpolate the peak to a small fraction of that.
>
> Independent noise does not shift the peak; it only adds random fluctuations to the *estimated* correlation, whose size decreases with averaging time (§21.7).

```desmos-graph
left=-1; right=6; bottom=-0.1; top=0.6
---
y=0.5e^{-\left|x-2\right|}|blue
y=0.5e^{-\left|x-2\right|}+0.04\sin(23x)e^{-0.1x}|grey|dotted
x=2|red|dashed
```

*Cross-correlation $R_{YX}(\tau)$ of Example 21.14 with $\tau$ in milliseconds (blue): the peak locates the delay $D=2$ ms. The grey dotted curve suggests what a finite-length estimate looks like — fluctuations from the noise ride on top, but the peak location is preserved.*

> [!tip] Engineering Note — correlators everywhere
> In discrete time, $R_{XX}[k]=E[X[n+k]X[n]]$ and $R_{XY}[k]=E[X[n+k]Y[n]]$. Correlators in GNSS receivers (code phase search), preamble synchronizers in Wi-Fi and Ethernet PHYs, lock-in amplifiers (correlation with a reference sinusoid), time-domain reflectometers and acoustic source localizers all compute sample cross-correlations against a known reference and then look for the peak. The mathematics of this section — especially "independent noise drops out of $R_{YX}$" — is the reason they work.

> [!warning] Common Pitfall — correlation is not causation, and zero correlation is not independence
> $R_{XY}\equiv0$ says only that there is no *linear* relation at any lag. $Y(t)=X^2(t)-E[X^2]$ for a zero-mean Gaussian $X$ is uncorrelated with $X$ (odd moments of a Gaussian vanish) but completely determined by it. For jointly Gaussian processes, and only then, uncorrelated implies independent (§21.5).

## 21.4 Power spectral density

### Motivation

Filters, amplifiers and channels are characterized by their frequency response, so the natural question about noise is *how its power is distributed over frequency*. For a deterministic energy signal $x(t)$, Parseval's theorem (Ch. 16) says $\lvert X(\omega)\rvert^2/(2\pi)$ is the energy per unit angular frequency. A sample function of a WSS process, however, has the same statistics forever, so its energy is infinite and its Fourier transform does not exist in the ordinary sense. The way out is to work with **power** per unit bandwidth, averaged over the ensemble. The remarkable fact — the Wiener–Khinchin theorem — is that this spectral power density is simply the Fourier transform of the autocorrelation function. Correlation in the time domain and power distribution in the frequency domain are two views of the same second-order information.

### Definition via the truncated transform

> [!definition] Power spectral density
> For a WSS process $X(t)$, let $X_T(t)=X(t)$ for $\lvert t\rvert\le T$ and $0$ otherwise, with Fourier transform $X_T(\omega)=\int_{-T}^{T}X(t)e^{-j\omega t}dt$. The **power spectral density** (PSD) is
> $$
> S_{XX}(\omega)=\lim_{T\to\infty}\frac{E\big[\lvert X_T(\omega)\rvert^2\big]}{2T}.
> $$

The quantity $\lvert X_T(\omega)\rvert^2/(2T)$ for one sample function is the **periodogram**: energy density of the truncated record divided by its duration, i.e. power density. The PSD is its ensemble average in the limit of a long record.

### The Wiener–Khinchin theorem

> [!theorem] Wiener–Khinchin theorem
> If $X(t)$ is WSS and $\int_{-\infty}^{\infty}\lvert\tau R_{XX}(\tau)\rvert\,d\tau<\infty$ (so that the limit below may be taken; the result holds in the distributional sense under much weaker conditions), then the PSD and the autocorrelation are a Fourier-transform pair:
> $$
> S_{XX}(\omega)=\int_{-\infty}^{\infty}R_{XX}(\tau)\,e^{-j\omega\tau}\,d\tau,\qquad
> R_{XX}(\tau)=\frac{1}{2\pi}\int_{-\infty}^{\infty}S_{XX}(\omega)\,e^{j\omega\tau}\,d\omega .
> $$

*Proof.* Write the squared magnitude as a double integral and take the expectation inside:

$$
E\big[\lvert X_T(\omega)\rvert^2\big]=E\Big[\int_{-T}^{T}\!\!\int_{-T}^{T}X(t)X(s)e^{-j\omega(t-s)}\,dt\,ds\Big]=\int_{-T}^{T}\!\!\int_{-T}^{T}R_{XX}(t-s)e^{-j\omega(t-s)}\,dt\,ds .
$$

The integrand depends only on $\tau=t-s$. In the $(t,s)$ square of side $2T$, the line $t-s=\tau$ has length proportional to $2T-\lvert\tau\rvert$ for $\lvert\tau\rvert\le2T$; changing variables to $(\tau,s)$ gives

$$
E\big[\lvert X_T(\omega)\rvert^2\big]=\int_{-2T}^{2T}(2T-\lvert\tau\rvert)\,R_{XX}(\tau)e^{-j\omega\tau}\,d\tau .
$$

Dividing by $2T$:

$$
\frac{E\lvert X_T(\omega)\rvert^2}{2T}=\int_{-2T}^{2T}\Big(1-\frac{\lvert\tau\rvert}{2T}\Big)R_{XX}(\tau)e^{-j\omega\tau}\,d\tau .
$$

As $T\to\infty$ the triangular weight tends to $1$ for each fixed $\tau$, and the error term $\frac{1}{2T}\int\lvert\tau R_{XX}(\tau)\rvert d\tau\to0$ under the stated condition. The limit is $\int R_{XX}(\tau)e^{-j\omega\tau}d\tau$. The inverse relation is the Fourier inversion theorem. $\square$

The triangular weight in the proof is worth remembering: a finite record of length $2T$ "sees" the autocorrelation through a triangular window, which is why spectral estimates from finite data are always smeared (§21.7).

A deep converse is **Bochner's theorem**: a continuous function $R(\tau)$ is non-negative definite (property 4 of §21.3) if and only if it is the inverse Fourier transform of a non-negative (possibly generalized) spectral measure. So a candidate autocorrelation is valid exactly when its transform is non-negative — the easiest test in practice (Example 21.10).

### Properties of the PSD

> [!theorem] Properties of the PSD
> 1. $S_{XX}(\omega)$ is **real** and **non-negative**: $S_{XX}(\omega)\ge0$.
> 2. For real $X$, $S_{XX}(-\omega)=S_{XX}(\omega)$ (even).
> 3. **Total power:**
> $$
> E[X^2(t)]=R_{XX}(0)=\frac{1}{2\pi}\int_{-\infty}^{\infty}S_{XX}(\omega)\,d\omega=\int_{-\infty}^{\infty}S_{XX}(2\pi f)\,df .
> $$
> 4. **Power in a band.** The power in the band $\omega_1\le\lvert\omega\rvert\le\omega_2$ (both signs of frequency counted) is
> $$
> P_{[\omega_1,\omega_2]}=\frac{1}{\pi}\int_{\omega_1}^{\omega_2}S_{XX}(\omega)\,d\omega .
> $$
> 5. **DC and periodic components.** A non-zero mean contributes $2\pi\mu_X^2\,\delta(\omega)$; a sinusoidal component of amplitude $A$ at $\omega_0$ contributes $\frac{\pi A^2}{2}[\delta(\omega-\omega_0)+\delta(\omega+\omega_0)]$.
> 6. **Zero-frequency value:** $S_{XX}(0)=\int R_{XX}(\tau)\,d\tau$ — the area under the autocorrelation. For a zero-mean process this governs the variance of long time averages (§21.7).

*Proofs.* (1) Non-negativity is immediate from the definition as a limit of $E\lvert X_T\rvert^2/(2T)\ge0$. Reality follows from $R_{XX}$ being even: $S_{XX}(\omega)=\int R_{XX}(\tau)\cos\omega\tau\,d\tau$. (2) Cosine is even. (3) Set $\tau=0$ in the inverse transform. (4) Integrate the PSD over $\omega_1\le\lvert\omega\rvert\le\omega_2$ and use evenness. That this really is "the power passed by an ideal bandpass filter" follows from the LTI result $S_{YY}=\lvert H\rvert^2S_{XX}$ below. (5) Transform of a constant and of a cosine. (6) Set $\omega=0$. $\square$

> [!warning] Common Pitfall — $\omega$ versus $f$, one-sided versus two-sided
> With the angular-frequency convention used here, power is $\frac{1}{2\pi}\int S_{XX}(\omega)\,d\omega$. Expressed per hertz, $S_X(f)=S_{XX}(2\pi f)$ and power is $\int_{-\infty}^{\infty}S_X(f)\,df$, with units $\mathrm{V^2/Hz}$. Instrument data sheets and spectrum analyzers usually quote the **one-sided** PSD $G_X(f)=2S_X(f)$ for $f\ge0$, so that power is $\int_0^\infty G_X(f)\,df$. Noise densities in $\mathrm{nV/\sqrt{Hz}}$ are $\sqrt{G_X(f)}$. Mixing these conventions produces factor-of-2 or factor-of-$2\pi$ errors; write the convention next to every PSD you quote.

### Standard transform pairs

With $\alpha>0$, $\operatorname{sinc}$-type expressions written out explicitly:

| $R_{XX}(\tau)$ | $S_{XX}(\omega)$ | origin |
|---|---|---|
| $\sigma^2e^{-\alpha\lvert\tau\rvert}$ | $\dfrac{2\alpha\sigma^2}{\alpha^2+\omega^2}$ | RC-filtered white noise, telegraph signal |
| $\dfrac{A^2}{2}\cos\omega_0\tau$ | $\dfrac{\pi A^2}{2}\big[\delta(\omega-\omega_0)+\delta(\omega+\omega_0)\big]$ | random-phase sinusoid |
| $\mu^2$ (constant) | $2\pi\mu^2\delta(\omega)$ | DC |
| $\dfrac{N_0}{2}\delta(\tau)$ | $\dfrac{N_0}{2}$ | white noise |
| $N_0B\,\dfrac{\sin 2\pi B\tau}{2\pi B\tau}$ | $\dfrac{N_0}{2}$ for $\lvert\omega\rvert<2\pi B$, else $0$ | band-limited white noise |
| $\sigma^2\big(1-\lvert\tau\rvert/T\big)$ for $\lvert\tau\rvert\le T$, else $0$ | $\sigma^2T\,\dfrac{\sin^2(\omega T/2)}{(\omega T/2)^2}$ | random NRZ data |
| $\sigma^2e^{-\alpha\tau^2}$ | $\sigma^2\sqrt{\pi/\alpha}\;e^{-\omega^2/(4\alpha)}$ | smooth (Gaussian-shaped) noise |
| $\sigma^2e^{-\alpha\lvert\tau\rvert}\cos\omega_0\tau$ | $\alpha\sigma^2\Big[\dfrac{1}{\alpha^2+(\omega-\omega_0)^2}+\dfrac{1}{\alpha^2+(\omega+\omega_0)^2}\Big]$ | narrowband (bandpass) noise |

The last pair follows from the first by the modulation theorem: multiplying $R_{XX}$ by $\cos\omega_0\tau$ shifts half the spectrum to $+\omega_0$ and half to $-\omega_0$.

### Worked examples

> [!example] Worked Example 21.15 — Exponential autocorrelation and the Lorentzian spectrum
> Let $R_{XX}(\tau)=\sigma^2e^{-\alpha\lvert\tau\rvert}$.
>
> **Step 1 — transform.** Split the integral at $\tau=0$:
> $$
> S_{XX}(\omega)=\sigma^2\left[\int_0^\infty e^{-(\alpha+j\omega)\tau}d\tau+\int_0^\infty e^{-(\alpha-j\omega)\tau}d\tau\right]=\sigma^2\left[\frac{1}{\alpha+j\omega}+\frac{1}{\alpha-j\omega}\right]=\frac{2\alpha\sigma^2}{\alpha^2+\omega^2}.
> $$
>
> **Step 2 — check the power.** $\frac{1}{2\pi}\int_{-\infty}^{\infty}\frac{2\alpha\sigma^2}{\alpha^2+\omega^2}d\omega=\frac{2\alpha\sigma^2}{2\pi}\cdot\frac{\pi}{\alpha}=\sigma^2=R_{XX}(0)$. ✓
>
> **Step 3 — power in a band.** By property 4, the power in $\lvert\omega\rvert\le\omega_1$ is
> $$
> \frac{1}{\pi}\int_0^{\omega_1}\frac{2\alpha\sigma^2}{\alpha^2+\omega^2}d\omega=\frac{2\sigma^2}{\pi}\arctan\frac{\omega_1}{\alpha}.
> $$
> With $\alpha=1000\ \text{s}^{-1}$ and a 1 kHz bandwidth ($\omega_1=2\pi\times1000$), the fraction is $\frac{2}{\pi}\arctan(6.283)=0.900$: 90% of the power lies below 1 kHz. Half the power lies below $\omega=\alpha$ (159 Hz), which is also the half-power (−3 dB) point of the spectrum.
>
> **Step 4 — interpretation.** Shorter correlation time (larger $\alpha$) means a wider, flatter spectrum with the same total power: time scale and bandwidth trade off exactly as for deterministic signals (Ch. 16). The value at the origin, $S_{XX}(0)=2\sigma^2/\alpha=\int R_{XX}d\tau$, confirms property 6.

```desmos-graph
left=-8; right=8; bottom=-0.1; top=2.2
---
y=\frac{2}{1+x^2}|blue
y=\frac{6}{9+x^2}|red
y=1|grey|dotted
```

*Lorentzian PSDs $S_{XX}(\omega)=2\alpha/(\alpha^2+\omega^2)$ for $\alpha=1$ (blue) and $\alpha=3$ (red), both with $\sigma^2=1$. The areas under the two curves are equal (same power); faster decorrelation spreads that power over a wider band. The grey line marks half the peak of the blue curve, reached at $\omega=\alpha=1$.*

> [!example] Worked Example 21.16 — Spectrum of random NRZ data
> From Example 21.7, $R_{XX}(\tau)=A^2(1-\lvert\tau\rvert/T)$ for $\lvert\tau\rvert<T$.
>
> **Step 1 — transform.** The triangle of half-width $T$ and height $1$ is the self-convolution of a rectangle of width $T$ and height $1/\sqrt T$; its transform is therefore the square of the rectangle's transform:
> $$
> S_{XX}(\omega)=A^2T\,\frac{\sin^2(\omega T/2)}{(\omega T/2)^2}.
> $$
> In hertz, $S_X(f)=A^2T\big[\sin(\pi fT)/(\pi fT)\big]^2$, with nulls at $f=k/T$, $k\ne0$.
>
> **Step 2 — power in the main lobe.** Total power is $A^2$. The fraction in $\lvert f\rvert<1/T$ is
> $$
> \int_{-1/T}^{1/T}T\Big[\frac{\sin\pi fT}{\pi fT}\Big]^2df=\int_{-1}^{1}\Big[\frac{\sin\pi u}{\pi u}\Big]^2du\approx0.903 .
> $$
> So a 1 Gb/s NRZ stream ($T=1$ ns) has 90.3% of its power below 1 GHz, while the remaining 9.7% in the sidelobes is what causes crosstalk and EMI — the motivation for pulse shaping.
>
> **Step 3 — check at DC.** $S_{XX}(0)=A^2T=\int R_{XX}\,d\tau$ (area of the triangle $=A^2T$). ✓ There is no impulse at DC because the data are equally likely $\pm A$ (zero mean). With unipolar data $\{0,A\}$ the mean is $A/2$ and an impulse $2\pi(A/2)^2\delta(\omega)$ appears.

> [!example] Worked Example 21.17 — Separating a sinusoid from noise in the spectrum
> $X(t)=2\cos(2\pi\cdot 50\,t+\Theta)+N(t)$, $\Theta\sim U[0,2\pi)$, with $N$ independent zero-mean noise of $R_{NN}(\tau)=0.5e^{-200\lvert\tau\rvert}$. Then
> $$
> R_{XX}(\tau)=2\cos(100\pi\tau)+0.5e^{-200\lvert\tau\rvert},\qquad
> S_{XX}(\omega)=2\pi\big[\delta(\omega-100\pi)+\delta(\omega+100\pi)\big]+\frac{200}{200^2+\omega^2}.
> $$
> The sinusoid's power $2$ sits in two impulses of weight $2\pi$ each (check: $\frac{1}{2\pi}\cdot2\cdot2\pi=2$); the noise's power $0.5$ is spread as a Lorentzian. In a measured spectrum with resolution bandwidth $\Delta f$, a line appears as a peak of height $\approx P_{\text{line}}/\Delta f$ while the noise floor stays at its density; halving $\Delta f$ raises the line 3 dB above the floor. In the autocorrelation the same separation appears in the lag domain: the cosine persists forever, the noise dies out within about $5/200=25$ ms.

### Cross-power spectral density

> [!definition] Cross-power spectral density
> For jointly WSS $X,Y$:
> $$
> S_{XY}(\omega)=\int_{-\infty}^{\infty}R_{XY}(\tau)e^{-j\omega\tau}d\tau,\qquad S_{YX}(\omega)=S_{XY}^*(\omega)\ \ (\text{real }X,Y).
> $$
> $S_{XY}$ is complex in general; its phase carries delay information. The **magnitude-squared coherence** is
> $$
> \gamma_{XY}^2(\omega)=\frac{\lvert S_{XY}(\omega)\rvert^2}{S_{XX}(\omega)S_{YY}(\omega)}\in[0,1],
> $$
> which measures how linearly related $X$ and $Y$ are at each frequency (Cauchy–Schwarz in the frequency domain).

For the delayed-echo model of Example 21.14, $R_{YX}(\tau)=0.5R_{XX}(\tau-D)$ gives $S_{YX}(\omega)=0.5\,e^{-j\omega D}S_{XX}(\omega)$: a linear phase of slope $-D$. Measuring that slope estimates the delay; this is the frequency-domain twin of peak-picking the cross-correlation.

### LTI systems with random inputs

This is the single most-used result of the chapter. Let a WSS process $X(t)$ drive a stable LTI system with impulse response $h(t)$ and frequency response $H(\omega)=\int h(t)e^{-j\omega t}dt$, producing

$$
Y(t)=\int_{-\infty}^{\infty}h(\alpha)X(t-\alpha)\,d\alpha .
$$

> [!theorem] Input–output relations for LTI systems
> $Y$ is WSS, jointly WSS with $X$, and
> $$
> \begin{aligned}
> \mu_Y &= \mu_X\,H(0), \\
> R_{YX}(\tau) &= h(\tau)*R_{XX}(\tau), & S_{YX}(\omega)&=H(\omega)\,S_{XX}(\omega),\\
> R_{XY}(\tau) &= h(-\tau)*R_{XX}(\tau), & S_{XY}(\omega)&=H^*(\omega)\,S_{XX}(\omega),\\
> R_{YY}(\tau) &= h(\tau)*h(-\tau)*R_{XX}(\tau), & S_{YY}(\omega)&=\lvert H(\omega)\rvert^2S_{XX}(\omega).
> \end{aligned}
> $$
> The output power is
> $$
> E[Y^2]=R_{YY}(0)=\frac{1}{2\pi}\int_{-\infty}^{\infty}\lvert H(\omega)\rvert^2S_{XX}(\omega)\,d\omega .
> $$
> If $X$ is Gaussian, so is $Y$ (§21.5).

*Derivation.*

- **Mean.** $E[Y(t)]=\int h(\alpha)E[X(t-\alpha)]d\alpha=\mu_X\int h(\alpha)d\alpha=\mu_XH(0)$ — the DC gain applied to the DC value.
- **Cross-correlation.** $R_{YX}(\tau)=E[Y(t+\tau)X(t)]=\int h(\alpha)E[X(t+\tau-\alpha)X(t)]\,d\alpha=\int h(\alpha)R_{XX}(\tau-\alpha)\,d\alpha=(h*R_{XX})(\tau)$. Independent of $t$, so $X,Y$ are jointly WSS.
- **Autocorrelation.** $R_{YY}(\tau)=E[Y(t+\tau)Y(t)]=\int h(\beta)E[Y(t+\tau)X(t-\beta)]\,d\beta=\int h(\beta)R_{YX}(\tau+\beta)\,d\beta$. With $\beta\to-\beta$ this is the convolution of $R_{YX}$ with $h(-\tau)$.
- **Spectra.** Transform the convolutions, using $\mathcal{F}\{h(-\tau)\}=H(-\omega)=H^*(\omega)$ for real $h$. Then $S_{YY}=H\,H^*\,S_{XX}=\lvert H\rvert^2S_{XX}$. $\square$

Notice that only $\lvert H\rvert$ affects the output PSD: the **phase response of a filter has no effect on the output noise power or spectrum**, although it does affect the cross-spectrum $S_{YX}$ and hence the shape of the noise waveforms.

More generally, if $Y_1=h_1*X$ and $Y_2=h_2*X$ are two filtered versions of the same input,

$$
S_{Y_1Y_2}(\omega)=H_1(\omega)H_2^*(\omega)S_{XX}(\omega).
$$

If the passbands of $H_1$ and $H_2$ do not overlap, $S_{Y_1Y_2}\equiv0$: **outputs of non-overlapping filters are uncorrelated** — the basis of filter-bank spectrum analyzers and of treating different frequency bins of noise as independent.

Two immediate uses of the theorem:

- **System identification with white noise.** If $S_{XX}=N_0/2$ then $R_{YX}(\tau)=\frac{N_0}{2}h(\tau)$: cross-correlating the output with a white input measures the impulse response directly. Pseudo-random binary sequences (PRBS), whose autocorrelation is very nearly an impulse, are used for exactly this in acoustics and in-circuit test.
- **Shaping filters.** Any rational PSD $S(\omega)=\lvert H(\omega)\rvert^2N_0/2$ can be generated by filtering white noise through $H$ — the basis of noise simulation and of Kalman-filter noise models. Finding a stable, minimum-phase $H$ from a given $S$ is called **spectral factorization**.

> [!example] Worked Example 21.18 — White noise through an RC lowpass
> $H(\omega)=1/(1+j\omega RC)$ and $S_{XX}=N_0/2$ (white noise, §21.6).
>
> **Step 1 — output PSD.** $S_{YY}=\lvert H\rvert^2S_{XX}=\dfrac{N_0/2}{1+\omega^2R^2C^2}$. Writing $\alpha=1/(RC)$, $S_{YY}=\dfrac{N_0}{2}\cdot\dfrac{\alpha^2}{\alpha^2+\omega^2}=\dfrac{N_0\alpha}{4}\cdot\dfrac{2\alpha}{\alpha^2+\omega^2}$.
>
> **Step 2 — output autocorrelation.** Matching with Example 21.15 ($\sigma^2=N_0\alpha/4$):
> $$
> R_{YY}(\tau)=\frac{N_0}{4RC}\,e^{-\lvert\tau\rvert/RC},\qquad E[Y^2]=\frac{N_0}{4RC}.
> $$
>
> **Step 3 — cross-check in the time domain.** $h(t)=\alpha e^{-\alpha t}u(t)$, so $R_{YY}(0)=\frac{N_0}{2}\int h^2dt=\frac{N_0}{2}\cdot\frac{\alpha^2}{2\alpha}=\frac{N_0\alpha}{4}$. ✓ (For white input, $E[Y^2]=\frac{N_0}{2}\int h^2(t)\,dt$ — the output power is the input density times the filter's impulse-response energy.)
>
> **Step 4 — numbers.** $R=1$ kΩ, $C=1\ \mu$F ($RC=1$ ms), $N_0/2=10^{-6}\ \mathrm{V^2/Hz}$: $E[Y^2]=2\times10^{-6}/(4\times10^{-3})=5\times10^{-4}\ \mathrm{V^2}$, i.e. $22.4$ mV rms. The output samples decorrelate on the time scale $RC$: white noise has become exponentially correlated (Lorentzian) noise.

```desmos-graph
left=-6; right=6; bottom=-0.1; top=1.3
---
y=1|grey|dashed
y=\frac{1}{1+x^2}|blue
(1,0.5)|label:half-power at ωRC=1
```

*White input PSD (normalized to 1, dashed) and the output PSD of a first-order RC lowpass, $\lvert H(\omega)\rvert^2=1/(1+\omega^2R^2C^2)$, plotted against $\omega RC$. The filter turns white noise into Lorentzian, exponentially correlated noise.*

> [!example] Worked Example 21.19 — Coloured noise through a differentiator-like highpass
> The noise of Example 21.15 ($\sigma^2=1\ \mathrm{V^2}$, $\alpha=1000\ \text{s}^{-1}$) is applied to a first-order highpass $H(\omega)=\dfrac{j\omega\tau_0}{1+j\omega\tau_0}$ with $\tau_0=1$ ms (corner at $\omega=1000$ rad/s, equal to $\alpha$). Find the output power.
>
> **Step 1.** $\lvert H\rvert^2=\dfrac{\omega^2\tau_0^2}{1+\omega^2\tau_0^2}=\dfrac{\omega^2}{\alpha^2+\omega^2}$ (since $1/\tau_0=\alpha$).
>
> **Step 2.** $S_{YY}=\dfrac{\omega^2}{\alpha^2+\omega^2}\cdot\dfrac{2\alpha\sigma^2}{\alpha^2+\omega^2}=\dfrac{2\alpha\sigma^2\omega^2}{(\alpha^2+\omega^2)^2}$.
>
> **Step 3.** Using $\int_{-\infty}^{\infty}\frac{\omega^2}{(\alpha^2+\omega^2)^2}d\omega=\frac{\pi}{2\alpha}$,
> $$
> E[Y^2]=\frac{1}{2\pi}\cdot2\alpha\sigma^2\cdot\frac{\pi}{2\alpha}=\frac{\sigma^2}{2}=0.5\ \mathrm{V^2}.
> $$
> A highpass whose corner coincides with the noise's own corner removes exactly half the noise power. The complementary lowpass $1/(1+j\omega\tau_0)$ passes $\frac{1}{2\pi}\int\frac{\alpha^2}{\alpha^2+\omega^2}\frac{2\alpha\sigma^2}{\alpha^2+\omega^2}d\omega=\sigma^2/2$ as well: the two outputs' powers add to the input power here because $\lvert H_{LP}\rvert^2+\lvert H_{HP}\rvert^2=1$ — a power-complementary pair.

### Discrete-time random sequences

All of the above has an exact discrete-time counterpart, which is what is actually computed in DSP systems.

> [!theorem] Discrete-time Wiener–Khinchin and LTI relations
> For a WSS sequence with $\sum_k\lvert R_{XX}[k]\rvert<\infty$,
> $$
> S_{XX}(e^{j\Omega})=\sum_{k=-\infty}^{\infty}R_{XX}[k]\,e^{-j\Omega k},\qquad R_{XX}[k]=\frac{1}{2\pi}\int_{-\pi}^{\pi}S_{XX}(e^{j\Omega})e^{j\Omega k}\,d\Omega,
> $$
> a $2\pi$-periodic, real, non-negative spectrum in the normalized frequency $\Omega$ (rad/sample). The power is $E[X^2[n]]=\frac{1}{2\pi}\int_{-\pi}^{\pi}S_{XX}(e^{j\Omega})d\Omega$. If $Y[n]=\sum_mh[m]X[n-m]$ with $H(e^{j\Omega})=\sum_mh[m]e^{-j\Omega m}$, then
> $$
> \mu_Y=\mu_XH(e^{j0}),\qquad S_{YX}(e^{j\Omega})=H(e^{j\Omega})S_{XX}(e^{j\Omega}),\qquad S_{YY}(e^{j\Omega})=\lvert H(e^{j\Omega})\rvert^2S_{XX}(e^{j\Omega}).
> $$
> In $z$-transform form, $S_{YY}(z)=H(z)H(z^{-1})S_{XX}(z)$.

The derivation is word-for-word the continuous one with sums for integrals. When a continuous process is sampled at rate $f_s=1/T_s$, $\Omega=\omega T_s$, and the sampled PSD is the aliased sum $S_{X_d}(e^{j\Omega})=\frac{1}{T_s}\sum_m S_{XX}\big(\frac{\Omega-2\pi m}{T_s}\big)$; for a process band-limited below $f_s/2$ there is no aliasing.

#### The AR(1) sequence

> [!definition] First-order autoregressive (AR(1)) sequence
> $$
> X[n]=aX[n-1]+W[n],\qquad\lvert a\rvert<1,
> $$
> where $W[n]$ is zero-mean white noise with variance $\sigma_W^2$ ("innovations"), with $W[n]$ uncorrelated with past values $X[n-1],X[n-2],\dots$.

This is the discrete-time analogue of RC-filtered white noise and the simplest model of a correlated sequence. We derive its statistics in three independent ways.

**(i) As a filtered white sequence.** The recursion is an LTI filter with $H(z)=\dfrac{1}{1-az^{-1}}$, $h[n]=a^nu[n]$. In steady state,

$$
S_{XX}(e^{j\Omega})=\frac{\sigma_W^2}{\lvert1-ae^{-j\Omega}\rvert^2}=\frac{\sigma_W^2}{1-2a\cos\Omega+a^2},
$$

and the variance is $\sigma_W^2\sum_{n\ge0}h^2[n]=\sigma_W^2\sum a^{2n}=\dfrac{\sigma_W^2}{1-a^2}$.

**(ii) Yule–Walker recursion.** Multiply the recursion by $X[n-k]$ ($k\ge1$) and take expectations. Since $W[n]$ is uncorrelated with the past, $E[W[n]X[n-k]]=0$, giving

$$
R_{XX}[k]=aR_{XX}[k-1]\quad(k\ge1)\qquad\Longrightarrow\qquad R_{XX}[k]=R_{XX}[0]\,a^{\lvert k\rvert}.
$$

For $k=0$: $R_{XX}[0]=E[(aX[n-1]+W[n])^2]=a^2R_{XX}[0]+\sigma_W^2$, so $R_{XX}[0]=\sigma_W^2/(1-a^2)$. Hence

$$
R_{XX}[k]=\frac{\sigma_W^2}{1-a^2}\,a^{\lvert k\rvert}.
$$

**(iii) Wiener–Khinchin check.** Summing the two geometric series,

$$
\sum_k\frac{\sigma_W^2}{1-a^2}a^{\lvert k\rvert}e^{-j\Omega k}=\frac{\sigma_W^2}{1-a^2}\Big[\frac{1}{1-ae^{-j\Omega}}+\frac{1}{1-ae^{j\Omega}}-1\Big]=\frac{\sigma_W^2}{1-2a\cos\Omega+a^2}. \ ✓
$$

The spectrum is lowpass for $0<a<1$ (peak $\sigma_W^2/(1-a)^2$ at $\Omega=0$) and highpass for $-1<a<0$ (peak at $\Omega=\pi$, with alternating-sign correlation $R[k]\propto a^{\lvert k\rvert}$). As $a\to1$ the correlation time $-1/\ln a\approx1/(1-a)$ samples grows and the spectrum concentrates near DC; at $a=1$ the recursion becomes the random walk of Example 21.2, which has no stationary state.

```desmos-graph
left=-3.1416; right=3.1416; bottom=0; top=12
---
y=\frac{1}{1-2(0.5)\cos(x)+0.25}|blue
y=\frac{1}{1-2(0.7)\cos(x)+0.49}|red
y=\frac{1}{1+2(0.7)\cos(x)+0.49}|green|dashed
```

*AR(1) power spectra $S_{XX}(e^{j\Omega})=1/(1-2a\cos\Omega+a^2)$ over one period $-\pi\le\Omega\le\pi$ ($\sigma_W^2=1$): $a=0.5$ (blue) and $a=0.7$ (red) are lowpass, sharper for larger $a$; $a=-0.7$ (green, dashed) is the mirror-image highpass spectrum.*

> [!example] Worked Example 21.20 — AR(1) statistics and start-up transient
> Take $a=0.9$, $\sigma_W^2=1$.
>
> **Step 1 — steady state.** $\sigma_X^2=1/(1-0.81)=5.263$; $\rho[k]=0.9^{\lvert k\rvert}$, so $\rho[10]=0.349$ and $\rho[22]=0.098$. The spectrum ranges from $S(e^{j0})=1/(0.1)^2=100$ down to $S(e^{j\pi})=1/(1.9)^2=0.277$ — a 25.6 dB lowpass. Its half-power point solves $1-1.8\cos\Omega+0.81=2(0.01)$, giving $\Omega_{3\text{dB}}=0.105$ rad/sample $\approx1-a$.
>
> **Step 2 — power check by integration.** $\frac{1}{2\pi}\int_{-\pi}^{\pi}\frac{d\Omega}{1.81-1.8\cos\Omega}=\frac{1}{\sqrt{1.81^2-1.8^2}}=\frac{1}{\sqrt{0.0361}}=5.263$. ✓ (Using $\frac{1}{2\pi}\int_{-\pi}^{\pi}\frac{d\Omega}{A-B\cos\Omega}=1/\sqrt{A^2-B^2}$ for $A>\lvert B\rvert$.)
>
> **Step 3 — transient.** If the recursion starts at $X[0]=0$, then $X[n]=\sum_{m=0}^{n-1}a^mW[n-m]$ and
> $$
> \operatorname{Var}X[n]=\sigma_W^2\frac{1-a^{2n}}{1-a^2}.
> $$
> After $n=10$ steps the variance is $4.62$ (88% of the steady state); after $n=50$, $5.263$ to four digits. The process is **asymptotically WSS**. In simulation, discard the first several correlation times of output ("burn-in"), or start from $X[0]\sim\mathcal{N}(0,\sigma_W^2/(1-a^2))$, which makes the sequence stationary from the start.

> [!example] Worked Example 21.21 — Fitting an AR(1) model to measured correlations
> The bias of a temperature sensor, sampled once per second, has estimated autocorrelation $\hat R[0]=2.0\ \mathrm{K^2}$ and $\hat R[1]=1.6\ \mathrm{K^2}$. Fit an AR(1) model and predict the bias one and ten seconds ahead.
>
> **Step 1 — parameters (Yule–Walker).** $\hat a=\hat R[1]/\hat R[0]=0.8$ and $\hat\sigma_W^2=\hat R[0](1-\hat a^2)=2.0\times0.36=0.72\ \mathrm{K^2}$.
>
> **Step 2 — one-step prediction.** Given the present value $X[n]=x$, the best linear predictor of $X[n+1]$ is $\hat X[n+1]=ax=0.8x$, with error $W[n+1]$ of variance $0.72$ — compared with $2.0$ without using $x$. Knowledge of the present removes 64% ($=a^2$) of the uncertainty.
>
> **Step 3 — ten steps ahead.** Iterating, $\hat X[n+10]=a^{10}x=0.107x$, error variance $\sigma_X^2(1-a^{20})=2.0(1-0.0115)=1.977\ \mathrm{K^2}$: after ten correlation times the present value is nearly useless for prediction.
>
> **Step 4 — continuous-time equivalent.** Matching $a=e^{-\alpha T_s}$ with $T_s=1$ s gives $\alpha=-\ln0.8=0.223\ \text{s}^{-1}$, i.e. a Gauss–Markov process (§21.5) with correlation time 4.5 s. This is the standard bias model in a Kalman filter.

> [!tip] Engineering Note — why AR models are everywhere
> AR($p$) models $X[n]=\sum_{i=1}^pa_iX[n-i]+W[n]$ generalize the above and are fitted by solving the Yule–Walker equations, a Toeplitz linear system in $R_{XX}[0],\dots,R_{XX}[p]$ (solved efficiently by the Levinson–Durbin recursion). They underlie linear-predictive speech coding (LPC in mobile-phone codecs), high-resolution "parametric" spectrum estimation, channel-fading simulators and the process models of Kalman filters.

## 21.5 Gaussian processes

### Motivation

Second-order theory tells us means, correlations and spectra, but an engineer often needs *probabilities*: How often will the noise exceed the decision threshold? What is the chance that a sensor error exceeds its specification? In general such questions need the full distribution. For one enormously important class — Gaussian processes — the mean and autocorrelation already determine every distribution, so everything computed in §§21.3–21.4 converts directly into probabilities. And Gaussian processes are not a mathematical convenience only: by the central limit theorem (Ch. 19), any signal formed as the sum of a very large number of small, independent contributions — the thermal agitation of $10^{20}$ electrons, the superposition of many weak interferers, the accumulated effect of many small manufacturing variations — is very nearly Gaussian.

### Definition

> [!definition] Gaussian process
> $X(t)$ is a **Gaussian process** if for every $k$ and every set of times $t_1,\dots,t_k$ the random vector $\mathbf{X}=[X(t_1),\dots,X(t_k)]^\mathsf{T}$ is jointly Gaussian (Ch. 19):
> $$
> f_{\mathbf{X}}(\mathbf{x})=\frac{1}{(2\pi)^{k/2}(\det\mathbf{C})^{1/2}}\exp\!\Big(-\tfrac12(\mathbf{x}-\boldsymbol{\mu})^\mathsf{T}\mathbf{C}^{-1}(\mathbf{x}-\boldsymbol{\mu})\Big),
> $$
> with $\mu_i=\mu_X(t_i)$ and $C_{ik}=C_{XX}(t_i,t_k)$. (If $\mathbf{C}$ is singular the vector is a degenerate Gaussian, defined through its characteristic function $\Phi_{\mathbf{X}}(\boldsymbol{\omega})=\exp\big(j\boldsymbol{\omega}^\mathsf{T}\boldsymbol{\mu}-\tfrac12\boldsymbol{\omega}^\mathsf{T}\mathbf{C}\boldsymbol{\omega}\big)$.)
> An equivalent definition: **every finite linear combination** $\sum_ia_iX(t_i)$ is a (possibly degenerate) Gaussian random variable.

Two processes $X,Y$ are **jointly Gaussian** if every vector formed from samples of both is jointly Gaussian.

> [!warning] Common Pitfall — Gaussian marginals are not enough
> A process whose first-order pdf $f_X(x;t)$ is Gaussian for every $t$ need **not** be a Gaussian process. Example: let $Z\sim\mathcal{N}(0,1)$ and $S=\pm1$ independent and equally likely; the pair $(Z,SZ)$ has Gaussian marginals, but $Z+SZ$ equals $0$ with probability $\tfrac12$, which no Gaussian variable does. The joint Gaussian property must hold for all finite collections of samples.

### Key properties

> [!theorem] Properties of Gaussian processes
> 1. **Second-order specification:** a Gaussian process is completely specified, in distribution, by $\mu_X(t)$ and $C_{XX}(t_1,t_2)$. Conversely, for any function $\mu(t)$ and any symmetric non-negative-definite $C(t_1,t_2)$ there exists a Gaussian process with these moments.
> 2. **WSS $\Leftrightarrow$ SSS** for Gaussian processes.
> 3. **Uncorrelated $\Rightarrow$ independent:** jointly Gaussian samples (or processes) that are uncorrelated are independent.
> 4. **Closure under linear operations:** if $X$ is Gaussian, then the output of any linear system — samples, sums, derivatives, integrals, LTI filtering — is Gaussian and jointly Gaussian with $X$.
> 5. **Linear conditional means:** for jointly Gaussian $X$ and $Y$, $E[Y\mid X]$ is an affine function of $X$, so the best (minimum mean-square error) estimator of one Gaussian quantity from others is linear, and the error variance does not depend on the observed values.

*Proof sketches.*

1. The joint pdf above involves only $\boldsymbol\mu$ and $\mathbf{C}$. Existence follows from Kolmogorov's extension theorem, because the Gaussian distributions with these parameters are consistent (marginals of Gaussians are Gaussian with the corresponding sub-vector and sub-matrix).
2. If $X$ is WSS, the mean vector and covariance matrix of $[X(t_1+c),\dots,X(t_k+c)]$ are $\mu_X\mathbf{1}$ and $[C_{XX}(t_i-t_k)]$, which do not involve $c$. By property 1 the joint distribution is therefore shift-invariant for every $k$: SSS.
3. If $\mathbf{C}$ is diagonal, the quadratic form in the exponent separates into a sum and the joint pdf factors into the product of the marginals.
4. A sample, sum or finite difference of Gaussian samples is a linear combination, hence Gaussian by the equivalent definition. An integral $\int h(\alpha)X(t-\alpha)d\alpha$ is a mean-square limit of Riemann sums, each Gaussian; mean-square limits of Gaussian variables are Gaussian (their characteristic functions $e^{j\omega\mu_n-\omega^2\sigma_n^2/2}$ converge to one of the same form). Derivatives are limits of difference quotients, likewise.
5. For jointly Gaussian $(X,Y)$ with correlation coefficient $\rho$, the conditional distribution of $Y$ given $X=x$ is $\mathcal{N}\big(\mu_Y+\rho\frac{\sigma_Y}{\sigma_X}(x-\mu_X),\ \sigma_Y^2(1-\rho^2)\big)$ (Ch. 19). $\square$

Property 4 is what makes Gaussian noise analysis tractable: to find the distribution of a filtered Gaussian noise sample we need only its mean and variance, both of which follow from $\mu_Y=\mu_XH(0)$ and $S_{YY}=\lvert H\rvert^2S_{XX}$ (§21.4). Non-Gaussian inputs, by contrast, generally produce outputs whose distributions are hard to compute — although narrowband filtering of almost any wideband noise tends to make the output *more* Gaussian (a central-limit effect: each output sample is a weighted sum of many nearly independent input samples).

### Worked examples

> [!example] Worked Example 21.22 — Probability of a large increment
> $X(t)$ is a zero-mean Gaussian WSS process with $R_{XX}(\tau)=4e^{-\lvert\tau\rvert}$ (units $\mathrm{V^2}$, $\tau$ in seconds). Find $P\{X(1)-X(0)>2\}$.
>
> **Step 1 — distribution of the increment.** $D=X(1)-X(0)$ is a linear combination of jointly Gaussian samples, hence Gaussian, with mean $0$ and
> $$
> \operatorname{Var}(D)=R_{XX}(0)+R_{XX}(0)-2R_{XX}(1)=8-8e^{-1}=5.057\ \mathrm{V^2},\qquad\sigma_D=2.249\ \mathrm{V}.
> $$
>
> **Step 2 — probability.** $P\{D>2\}=Q(2/2.249)=Q(0.889)\approx0.187$, where $Q(x)=1-\Phi(x)$.
>
> **Step 3 — comparison.** $P\{X(0)>2\}=Q(2/2)=Q(1)\approx0.159$. The correlation between $X(0)$ and $X(1)$ is $\rho=e^{-1}=0.368$; if the samples were independent, $\operatorname{Var}(D)$ would be $8$, and if perfectly correlated it would be $0$. The general formula $\operatorname{Var}[X(t+\tau)-X(t)]=2[R_{XX}(0)-R_{XX}(\tau)]$ — the **structure function** — is used to specify jitter and drift over a given time interval.

> [!example] Worked Example 21.23 — Predicting a Gaussian process
> Same process as Example 21.22. The value $X(0)=3$ V is observed. Find the distribution of $X(0.5)$ and the probability that it still exceeds 3 V.
>
> **Step 1 — correlation.** $\rho=R_{XX}(0.5)/R_{XX}(0)=e^{-0.5}=0.607$.
>
> **Step 2 — conditional distribution** (property 5 with $\mu=0$, $\sigma_X=\sigma_Y=2$):
> $$
> X(0.5)\mid X(0)=3\ \sim\ \mathcal{N}\big(\rho\cdot3,\ 4(1-\rho^2)\big)=\mathcal{N}(1.820,\ 2.528).
> $$
> The conditional standard deviation is $1.590$ V.
>
> **Step 3 — probability.** $P\{X(0.5)>3\mid X(0)=3\}=Q\big((3-1.820)/1.590\big)=Q(0.742)\approx0.229$, compared with the unconditional $P\{X(0.5)>3\}=Q(1.5)=0.067$. Knowing that the noise is high now makes it more than three times as likely to be high half a second later.
>
> **Remark.** Because the conditional variance $4(1-\rho^2)$ does not depend on the observed value, the prediction error of this optimal linear predictor is known in advance — the fact that makes Kalman-filter covariance computations possible before any data arrive. For the exponential autocorrelation, observing $X$ at earlier times $t<0$ in addition to $X(0)$ would not improve the prediction: this is the Markov property of the Gauss–Markov process (below).

> [!example] Worked Example 21.24 — Zero-crossing rate (Rice's formula)
> For a zero-mean stationary Gaussian process with twice-differentiable $R_{XX}$, the expected number of crossings of the level $u$ per unit time (both directions) is **Rice's formula**
> $$
> \nu_u=\frac{1}{\pi}\sqrt{\frac{-R_{XX}''(0)}{R_{XX}(0)}}\;\exp\!\Big(-\frac{u^2}{2R_{XX}(0)}\Big).
> $$
> *Sketch.* A crossing occurs in $(t,t+dt)$ if $X(t)$ is within $\lvert X'(t)\rvert dt$ of $u$. Since $X(t)$ and $X'(t)$ are uncorrelated (§21.3), hence independent (Gaussian), the expected count is $dt\,f_X(u)\,E\lvert X'\rvert=dt\,f_X(u)\sqrt{2/\pi}\,\sigma_{X'}$. With $\sigma_{X'}^2=-R''_{XX}(0)$ and $f_X(u)=e^{-u^2/2\sigma^2}/(\sqrt{2\pi}\sigma)$ the formula follows.
>
> **Check on a sinusoid.** For $R_{XX}=\cos\omega_0\tau$, $\nu_0=\omega_0/\pi=2f_0$: two zero crossings per period. ✓
>
> **Band-limited noise.** For flat noise on $\lvert f\rvert<B$, $R_{XX}(\tau)=\sigma^2\frac{\sin2\pi B\tau}{2\pi B\tau}$ and $-R''_{XX}(0)/R_{XX}(0)=(2\pi B)^2/3$, so $\nu_0=2B/\sqrt3$. With $B=1$ kHz: $\nu_0=1155$ crossings per second. Crossings of $u=3\sigma$ (both directions) occur at $1155\,e^{-4.5}=12.8$ per second, i.e. about 6.4 upward excursions per second. This is how false-trigger rates of comparators and false-alarm rates of threshold detectors are estimated.

### Gauss–Markov and Wiener processes

> [!definition] Gauss–Markov (Ornstein–Uhlenbeck) process
> The stationary Gaussian process with
> $$
> R_{XX}(\tau)=\sigma^2e^{-\alpha\lvert\tau\rvert},\qquad S_{XX}(\omega)=\frac{2\alpha\sigma^2}{\alpha^2+\omega^2}
> $$
> is the only stationary Gaussian process (up to trivial cases) that is also **Markov**: given $X(t)$, the future is independent of the past. It is the output of a first-order lowpass filter driven by white Gaussian noise, $\dot X=-\alpha X+W$, and its samples at spacing $T_s$ form an AR(1) sequence with $a=e^{-\alpha T_s}$ (Example 21.4).

It is the standard model for slowly varying sensor biases and drifts (gyroscope bias, clock frequency error, the bias of Example 21.21) and for the velocity of a Brownian particle.

> [!definition] Wiener process (Brownian motion)
> $W(t)$, $t\ge0$, is the Gaussian process with $W(0)=0$, independent stationary increments, $W(t+s)-W(t)\sim\mathcal{N}(0,\sigma^2s)$. Equivalently, it is the zero-mean Gaussian process with
> $$
> R_{WW}(t_1,t_2)=\sigma^2\min(t_1,t_2).
> $$

The Wiener process is the continuous-time limit of the random walk: scale steps of size $\pm\sigma\sqrt{\Delta t}$ every $\Delta t$ seconds; by the CLT, the walk at time $t$ is approximately $\mathcal{N}(0,\sigma^2t)$, and its autocorrelation $\sigma^2\Delta t\min(n_1,n_2)\to\sigma^2\min(t_1,t_2)$. It is not stationary ($\operatorname{Var}W(t)=\sigma^2t$), its paths are continuous but nowhere differentiable, and — formally — its derivative is white Gaussian noise of PSD $\sigma^2$ (§21.6): $\frac{\partial^2}{\partial t_1\partial t_2}\sigma^2\min(t_1,t_2)=\sigma^2\delta(t_1-t_2)$. Free-running oscillator phase noise is modelled as a Wiener process, which is why phase error grows like $\sqrt t$ and why the phase-noise spectrum of an oscillator falls as $1/f^2$ away from the carrier.

### Bandpass Gaussian noise

At the output of a narrowband IF or RF filter centred at $\omega_0$, noise looks like a sinusoid at $\omega_0$ with slowly wandering amplitude and phase. This is made precise by the quadrature representation.

> [!theorem] In-phase/quadrature representation of narrowband noise
> Let $N(t)$ be zero-mean, WSS, Gaussian with PSD concentrated in $\lvert\omega\mp\omega_0\rvert<W$, $W<\omega_0$. Then
> $$
> N(t)=N_c(t)\cos\omega_0t-N_s(t)\sin\omega_0t=R(t)\cos\big(\omega_0t+\Phi(t)\big),
> $$
> where the in-phase and quadrature components $N_c,N_s$ are jointly Gaussian, zero-mean, lowpass processes, each with the same power as $N$:
> $$
> E[N_c^2]=E[N_s^2]=E[N^2]=\sigma^2,\qquad S_{N_cN_c}(\omega)=S_{N_sN_s}(\omega)=\big[S_{NN}(\omega-\omega_0)+S_{NN}(\omega+\omega_0)\big]_{\lvert\omega\rvert<W}.
> $$
> If $S_{NN}$ is symmetric about $\pm\omega_0$, $N_c$ and $N_s$ are independent. At any instant the **envelope** $R=\sqrt{N_c^2+N_s^2}$ is **Rayleigh** distributed and the phase $\Phi$ is uniform on $[0,2\pi)$, independent of $R$:
> $$
> f_R(r)=\frac{r}{\sigma^2}e^{-r^2/(2\sigma^2)}\ (r\ge0),\qquad P\{R>r\}=e^{-r^2/(2\sigma^2)} .
> $$

The Rayleigh result follows from transforming the pdf of two independent $\mathcal{N}(0,\sigma^2)$ variables to polar coordinates (Ch. 19): $f_{R\Phi}(r,\phi)=\frac{r}{2\pi\sigma^2}e^{-r^2/2\sigma^2}$, which factors.

> [!example] Worked Example 21.25 — Envelope-detector false-alarm threshold
> A radar receiver's IF noise is Gaussian narrowband noise with power $\sigma^2$. An envelope detector declares a target when $R(t)$ exceeds a threshold $r_0$. Choose $r_0$ so that the probability of a false alarm on any single sample is $P_{fa}=10^{-6}$.
>
> **Step 1.** $P\{R>r_0\}=e^{-r_0^2/(2\sigma^2)}=10^{-6}$, so $r_0=\sigma\sqrt{2\ln10^6}=\sigma\sqrt{27.63}=5.26\,\sigma$.
>
> **Step 2 — compare with a baseband threshold.** For a real Gaussian sample, $P\{N>5.26\sigma\}=Q(5.26)\approx7.2\times10^{-8}$, about 14 times smaller. The envelope combines two independent Gaussian components, so large values are more probable; using the Gaussian tail here would make the false-alarm rate 14 times worse than designed.
>
> **Step 3 — in decibels.** The threshold is $20\log_{10}5.26=14.4$ dB above the RMS noise of each component, which sets the minimum detectable signal of the radar.

> [!tip] Engineering Note — where Gaussianity fails
> The Gaussian model is excellent for thermal noise and usually for the sum of many independent interferers, but it fails for impulsive noise (lightning, ignition, switching transients), for a small number of dominant interferers, for quantization error of a slowly varying signal (closer to uniform), and for clipped or rectified signals. Heavy-tailed noise makes rare large errors far more likely than $Q(\cdot)$ predicts; robust detectors and median filters are the response. Always check a histogram, ideally on a logarithmic probability axis where the Gaussian tail is a parabola.

## 21.6 White noise

### Motivation

Every receiver, amplifier and sensor adds noise, and in most cases that noise is spread almost uniformly over a frequency range far wider than the circuit that processes it. The thermal noise of a resistor is flat to within 1% up to about 100 GHz at room temperature; the circuits that see it have bandwidths of hertz to gigahertz. In such cases the detailed spectrum of the noise beyond the circuit's passband is irrelevant — only its *density* inside the passband matters. **White noise** is the idealization that makes this explicit: a flat PSD at all frequencies, characterized by one number. It is the "impulse" of random-signal analysis — unphysical on its own, indispensable as an input.

### Definition and interpretation

> [!definition] White noise
> A zero-mean WSS process $W(t)$ is **white** if its PSD is constant over all frequencies:
> $$
> S_{WW}(\omega)=\frac{N_0}{2}\quad\Longleftrightarrow\quad R_{WW}(\tau)=\frac{N_0}{2}\,\delta(\tau).
> $$
> If in addition $W$ is Gaussian, it is **white Gaussian noise** (WGN). The factor $\tfrac12$ is a convention that makes $N_0$ the **one-sided** PSD, in $\mathrm{V^2/Hz}$ (or $\mathrm{W/Hz}$ for available power).

Consequences of the definition:

- **Infinite power.** $R_{WW}(0)=\infty$; equivalently $\frac{1}{2\pi}\int\frac{N_0}{2}d\omega$ diverges. White noise is not a process with finite-variance samples; it is a generalized process that only makes sense after filtering or integration, just as $\delta(t)$ only makes sense under an integral.
- **Uncorrelated at every separation.** $R_{WW}(\tau)=0$ for every $\tau\neq0$, however small: the process has no memory at all.
- **Filtered white noise is well defined.** By §21.4, a filter $h$ with finite energy turns it into a finite-power process with
$$
S_{YY}(\omega)=\frac{N_0}{2}\lvert H(\omega)\rvert^2,\qquad R_{YY}(\tau)=\frac{N_0}{2}\int h(t+\tau)h(t)\,dt,\qquad E[Y^2]=\frac{N_0}{2}\int_{-\infty}^{\infty}h^2(t)\,dt .
$$
The output autocorrelation is the *deterministic autocorrelation of the impulse response*: the filter "imprints" its own memory on the noise.
- **Integrated white noise is a Wiener process.** $\int_0^tW(s)\,ds$ has variance $\frac{N_0}{2}t$ and independent increments (§21.5).

White noise is useful precisely *because* physical noise is flat over a band far wider than any system that processes it, so only $N_0$ matters and the system's own bandwidth sets the output power.

### Band-limited white noise

> [!definition] Band-limited (lowpass) white noise
> $$
> S_{XX}(\omega)=\begin{cases}N_0/2,&\lvert\omega\rvert<2\pi B\\0,&\text{otherwise}\end{cases}
> \quad\Longrightarrow\quad
> R_{XX}(\tau)=N_0B\,\frac{\sin(2\pi B\tau)}{2\pi B\tau}.
> $$

*Derivation.* $R_{XX}(\tau)=\frac{1}{2\pi}\int_{-2\pi B}^{2\pi B}\frac{N_0}{2}e^{j\omega\tau}d\omega=\frac{N_0}{2}\cdot\frac{1}{2\pi}\cdot\frac{2\sin(2\pi B\tau)}{\tau}=N_0B\,\frac{\sin(2\pi B\tau)}{2\pi B\tau}$.

The power is $R_{XX}(0)=N_0B$: **noise power = one-sided density × bandwidth**. This is the single most-used formula in link budgets.

```desmos-graph
left=-2.5; right=2.5; bottom=-0.35; top=1.2
---
y=\frac{\sin(2\pi x)}{2\pi x}|blue
(0.5,0)|label:1/(2B)
(1,0)|label:2/(2B)
(-0.5,0)
(-1,0)
```

*Normalized autocorrelation of band-limited white noise with $B=1$: zeros at $\tau=k/(2B)$. Samples taken at the Nyquist rate $2B$ are uncorrelated — and, if the noise is Gaussian, independent.*

Similarly, **bandpass white noise** with density $N_0/2$ on $\lvert f\mp f_0\rvert<B/2$ has $R_{XX}(\tau)=N_0B\frac{\sin\pi B\tau}{\pi B\tau}\cos2\pi f_0\tau$ and power $N_0B$ — a sinc envelope on a carrier, consistent with the narrowband representation of §21.5.

### Discrete-time white noise

> [!theorem] Discrete-time white noise and sampling
> A sequence $W[n]$ with $E[W[n]]=0$ and $R_{WW}[k]=\sigma_W^2\,\delta[k]$ has the flat spectrum $S_{WW}(e^{j\Omega})=\sigma_W^2$ for all $\Omega$, and — unlike its continuous-time counterpart — finite power $\sigma_W^2$. Sampling band-limited white noise of bandwidth $B$ at exactly $f_s=2B$ produces discrete white noise with $\sigma_W^2=N_0B$.

*Proof of the sampling statement.* $R_{XX}[k]=R_{XX}(k/(2B))=N_0B\frac{\sin\pi k}{\pi k}=N_0B\,\delta[k]$. $\square$

Sampling at a higher rate $f_s>2B$ gives correlated samples (sinc values at non-integer multiples of $\pi$); sampling unfiltered wideband noise at $f_s$ aliases all of its power into $[-f_s/2,f_s/2]$. This is why an anti-aliasing filter is needed even when the *signal* is band-limited: without it, wideband noise is folded into the band of interest and the SNR degrades.

> [!example] Worked Example 21.26 — Noise power in a receiver bandwidth
> White noise with two-sided PSD $N_0/2=10^{-8}\ \mathrm{V^2/Hz}$ passes through an ideal lowpass filter of bandwidth $B=5$ kHz.
>
> **Step 1.** $P=N_0B=2\times10^{-8}\times5\times10^{3}=10^{-4}\ \mathrm{V^2}$, so $V_{\text{rms}}=10\ \mathrm{mV}$.
>
> **Step 2 — scaling.** Doubling the bandwidth doubles the noise power (+3 dB) and raises the RMS noise by $\sqrt2$ to $14.1$ mV. Halving it gives $7.07$ mV. In general $V_{\text{rms}}\propto\sqrt B$, which is why noise densities are quoted in $\mathrm{V/\sqrt{Hz}}$: here $\sqrt{N_0}=\sqrt{2\times10^{-8}}=141\ \mu\mathrm{V}/\sqrt{\mathrm{Hz}}$, and $141\ \mu\text{V}\times\sqrt{5000}=10$ mV. ✓
>
> **Step 3 — samples.** If this output is sampled at $f_s=10$ kHz $=2B$, the samples are uncorrelated with variance $10^{-4}\ \mathrm{V^2}$; if the noise is Gaussian, $P\{\lvert v\rvert>30\ \text{mV}\}=2Q(3)=0.0027$ per sample, i.e. about 27 samples per second exceed $\pm30$ mV.

### Noise-equivalent bandwidth

Real filters do not have brick-wall responses. To keep the convenient formula "power = density × bandwidth", we define the bandwidth that makes it exact.

> [!definition] Noise-equivalent bandwidth
> For a lowpass filter with DC gain $H(0)\neq0$, the **noise-equivalent bandwidth** is
> $$
> B_N=\frac{1}{\lvert H(0)\rvert^2}\int_0^\infty\lvert H(2\pi f)\rvert^2\,df ,
> $$
> the width (in Hz) of an ideal rectangular filter of the same peak gain that passes the same white-noise power. For a bandpass filter, replace $\lvert H(0)\rvert$ by the peak gain $\lvert H(2\pi f_0)\rvert$ and integrate over $f>0$.

With this definition the output power of a filter driven by white noise of one-sided density $N_0$ is exactly

$$
E[Y^2]=N_0\,B_N\,\lvert H(0)\rvert^2 .
$$

By Parseval's theorem there is an equivalent **time-domain formula**:

$$
B_N=\frac{\int_{-\infty}^{\infty}h^2(t)\,dt}{2\left(\int_{-\infty}^{\infty}h(t)\,dt\right)^2},
$$

since $\int_0^\infty\lvert H(2\pi f)\rvert^2df=\frac12\int_{-\infty}^\infty\lvert H\rvert^2df=\frac12\int h^2dt$ and $H(0)=\int h\,dt$.

```desmos-graph
left=0; right=4; bottom=0; top=1.2
---
y=\frac{1}{1+x^2}|blue
y=1|red|dashed|0<x<1.5708
x=1.5708|red|dashed|0<y<1
```

*Noise-equivalent bandwidth of a first-order lowpass, with frequency normalized to $f_{3\text{dB}}$: the area under $\lvert H\rvert^2=1/(1+x^2)$ (blue) over $x>0$ is $\pi/2$, equal to the area of the red rectangle of height $1$ and width $\pi/2\approx1.571$.*

> [!example] Worked Example 21.27 — Noise bandwidth of RC filters
> **(a) Single RC lowpass.** $\lvert H(2\pi f)\rvert^2=1/(1+(2\pi fRC)^2)$, $H(0)=1$:
> $$
> B_N=\int_0^\infty\frac{df}{1+(2\pi fRC)^2}=\frac{1}{2\pi RC}\cdot\frac{\pi}{2}=\frac{1}{4RC}=\frac{\pi}{2}f_{3\text{dB}},
> $$
> where $f_{3\text{dB}}=1/(2\pi RC)$. The time-domain formula agrees: $h=\frac{1}{RC}e^{-t/RC}u(t)$, $\int h^2=\frac{1}{2RC}$, $\int h=1$, so $B_N=\frac{1}{4RC}$. ✓ With $RC=1$ ms, $f_{3\text{dB}}=159$ Hz and $B_N=250$ Hz: using $f_{3\text{dB}}$ in place of $B_N$ would underestimate the noise power by 36% (1.96 dB).
>
> **(b) Two identical RC sections, buffered.** $\lvert H\rvert^2=1/(1+(2\pi fRC)^2)^2$. With $u=2\pi fRC$ and $\int_0^\infty\frac{du}{(1+u^2)^2}=\frac{\pi}{4}$,
> $$
> B_N=\frac{1}{2\pi RC}\cdot\frac{\pi}{4}=\frac{1}{8RC}.
> $$
> The 3 dB frequency of the cascade solves $(1+u^2)^2=2$, $u=\sqrt{\sqrt2-1}=0.644$, so $f_{3\text{dB}}=0.1024/RC$ and $B_N/f_{3\text{dB}}=1.22$. Steeper skirts bring $B_N$ closer to $f_{3\text{dB}}$.

The ratio $B_N/f_{3\text{dB}}$ for an $n$-th order Butterworth lowpass is $\dfrac{\pi/(2n)}{\sin(\pi/(2n))}$:

| order $n$ | 1 | 2 | 3 | 4 | $\infty$ |
|---|---|---|---|---|---|
| $B_N/f_{3\text{dB}}$ | 1.571 | 1.111 | 1.047 | 1.026 | 1 |

> [!example] Worked Example 21.28 — Averaging as filtering
> A **moving-average** (integrate-and-dump) filter outputs $Y(t)=\frac1T\int_{t-T}^{t}X(s)\,ds$, i.e. $h(t)=\frac1T$ on $[0,T]$. Then $\int h=1$, $\int h^2=1/T$ and
> $$
> B_N=\frac{1/T}{2\cdot1^2}=\frac{1}{2T}.
> $$
> Averaging white noise of density $N_0$ for $T$ seconds gives output variance $N_0/(2T)$. A digital voltmeter integrating over one power-line cycle ($T=20$ ms at 50 Hz) has $B_N=25$ Hz; increasing the integration time to 10 cycles reduces the noise bandwidth to 2.5 Hz and the RMS noise by $\sqrt{10}$. The frequency response, $\lvert H(2\pi f)\rvert=\lvert\sin\pi fT/(\pi fT)\rvert$, also has nulls at multiples of $1/T$ — which is why integrating over exact line cycles rejects hum.

> [!example] Worked Example 21.29 — Noise bandwidth of a tuned circuit
> A second-order bandpass (e.g. a parallel RLC tank driven by a current, with voltage output) has
> $$
> \lvert H(2\pi f)\rvert^2=\frac{1}{1+Q^2\big(f/f_0-f_0/f\big)^2}
> $$
> with peak gain $1$ at $f_0$ and $B_{3\text{dB}}=f_0/Q$. The integral (by contour methods, Ch. 18, or numerically) is
> $$
> B_N=\int_0^\infty\lvert H\rvert^2df=\frac{\pi f_0}{2Q}=\frac{\pi}{2}B_{3\text{dB}} .
> $$
> For a 10.7 MHz IF filter with $Q=50$: $B_{3\text{dB}}=214$ kHz, $B_N=336$ kHz. A numerical check with $f_0=1$, $Q=50$ gives $0.03138$ against $\pi/100=0.03142$ (the small difference is the truncated tail). The same $\pi/2$ factor as the first-order lowpass appears because the bandpass response is the lowpass response mapped onto the frequency axis.

### Coloured noise

Noise whose PSD is not flat is called **coloured**. The common types are named by their spectral slope:

| Name | One-sided PSD | Autocorrelation behaviour | Where it occurs |
|---|---|---|---|
| white | $G(f)=N_0$ | $\delta(\tau)$ | thermal, shot noise |
| pink / flicker | $\propto1/f$ | very long memory, no finite variance over $0<f<\infty$ | MOSFETs, carbon resistors, oscillator close-in phase noise |
| brown / red | $\propto1/f^2$ | integral of white noise (Wiener) | random walk, free-running oscillator phase |
| Lorentzian | $\propto1/(1+(f/f_c)^2)$ | $e^{-2\pi f_c\lvert\tau\rvert}$ | RC-filtered white, RTS noise |
| blue / violet | $\propto f$, $\propto f^2$ | derivative-like | differentiated noise, noise-shaped quantization |

> [!example] Worked Example 21.30 — Integrating a white-plus-flicker noise density
> An op-amp has input voltage noise with one-sided density $G(f)=e_w^2\big(1+f_c/f\big)$, white level $e_w=10\ \mathrm{nV/\sqrt{Hz}}$ and **1/f corner** $f_c=100$ Hz. Find the RMS noise from $f_1=0.1$ Hz to $f_2=10$ kHz.
>
> **Step 1 — integrate.**
> $$
> \overline{v^2}=\int_{f_1}^{f_2}e_w^2\Big(1+\frac{f_c}{f}\Big)df=e_w^2\Big[(f_2-f_1)+f_c\ln\frac{f_2}{f_1}\Big].
> $$
>
> **Step 2 — numbers.** $e_w^2=10^{-16}\ \mathrm{V^2/Hz}$; $(f_2-f_1)=9999.9$ Hz; $f_c\ln(10^5)=100\times11.51=1151$ Hz. Total $\overline{v^2}=10^{-16}\times11151=1.115\times10^{-12}\ \mathrm{V^2}$, so $v_{\text{rms}}=1.06\ \mu$V.
>
> **Step 3 — interpretation.** The 1/f part contributes the equivalent of 1151 Hz of white bandwidth, i.e. 10% of the total power. Its contribution grows only logarithmically with the observation time ($f_1\approx1/T_{\text{obs}}$): each additional decade at low frequency adds $f_c\ln10=230$ Hz. The variance does not converge as $f_1\to0$ — flicker noise has no well-defined stationary variance, a symptom of its very long memory.

> [!warning] Common Pitfall — "white" is not "Gaussian"
> "White" refers to the *spectrum* (uncorrelated samples); "Gaussian" refers to the *amplitude distribution*. They are independent properties. A random $\pm1$ i.i.d. sequence is white but not Gaussian; the Gauss–Markov process is Gaussian but not white; thermal noise is both. Statements such as "white noise has a Gaussian distribution" or "Gaussian noise has a flat spectrum" are false in general.

> [!tip] Engineering Note — quantization noise as discrete white noise
> An ADC with step $\Delta$, digitizing a busy signal, adds an error that is well modelled as discrete white noise, uniformly distributed on $[-\Delta/2,\Delta/2]$, with variance $\Delta^2/12$ spread evenly over $\lvert f\rvert<f_s/2$. Its one-sided density is therefore $\Delta^2/(6f_s)$: **oversampling** by a factor $M$ (then filtering to the signal band) reduces the in-band quantization noise by $M$ — 3 dB, or half a bit, per doubling of $f_s$. Noise-shaping ($\Sigma\Delta$) converters go further by making the quantization noise "blue", pushing it out of band.

## 21.7 Ergodicity

### Motivation

Everything so far has been defined by ensemble averages — expectations over infinitely many hypothetical repetitions of the experiment. In the laboratory we usually have **one** long record: one amplifier, one data capture, one night of sensor logging. A spectrum analyzer, an RMS voltmeter and an oscilloscope's "average" mode all compute *time* averages over that single record. The property that justifies replacing $E[\cdot]$ by $\langle\cdot\rangle$ is **ergodicity**. It is not automatic, and the theory below tells us both *when* it holds and *how long* the record must be for the time average to be accurate — the practical question behind every measurement-time specification.

### Time averages as random variables

For a finite record, define the time averages

$$
\langle X\rangle_T=\frac{1}{2T}\int_{-T}^{T}X(t)\,dt,\qquad
\hat R_T(\tau)=\frac{1}{2T}\int_{-T}^{T}X(t+\tau)X(t)\,dt .
$$

Because they are computed from a random sample function, these are **random variables**. For a WSS process, interchanging expectation and integration,

$$
E\big[\langle X\rangle_T\big]=\frac{1}{2T}\int_{-T}^{T}\mu_X\,dt=\mu_X ,
$$

so $\langle X\rangle_T$ is an **unbiased estimator** of $\mu_X$ (Ch. 20). The question is whether its variance shrinks to zero as $T$ grows.

> [!definition] Ergodicity in the mean and in autocorrelation
> A WSS process is **mean-ergodic** if $\langle X\rangle_T\to\mu_X$ in mean square as $T\to\infty$:
> $$
> \lim_{T\to\infty}E\big[(\langle X\rangle_T-\mu_X)^2\big]=0 .
> $$
> It is **autocorrelation-ergodic** if $\hat R_T(\tau)\to R_{XX}(\tau)$ in mean square for every $\tau$, and **ergodic** (in the strict sense) if *every* time average of every measurable function of the process converges to the corresponding ensemble average. For a sequence, replace $\frac{1}{2T}\int_{-T}^{T}$ by $\frac{1}{N}\sum_{n=0}^{N-1}$.

### The mean-ergodic theorem

> [!theorem] Mean-ergodic theorem (Slutsky)
> For a WSS process,
> $$
> \operatorname{Var}\big(\langle X\rangle_T\big)=\frac{1}{2T}\int_{-2T}^{2T}C_{XX}(\tau)\Big(1-\frac{\lvert\tau\rvert}{2T}\Big)d\tau ,
> $$
> and $X$ is mean-ergodic if and only if this tends to zero. Equivalently (Slutsky's condition), $\frac{1}{2T}\int_{-2T}^{2T}C_{XX}(\tau)\,d\tau\to0$. A **sufficient** condition is $\int_{-\infty}^{\infty}\lvert C_{XX}(\tau)\rvert\,d\tau<\infty$, in which case, for large $T$,
> $$
> \operatorname{Var}\big(\langle X\rangle_T\big)\approx\frac{1}{2T}\int_{-\infty}^{\infty}C_{XX}(\tau)\,d\tau=\frac{S_{CC}(0)}{2T},
> $$
> where $S_{CC}(0)$ is the PSD of the zero-mean part $X-\mu_X$ at zero frequency.
> For a WSS sequence with $N$ samples: $\operatorname{Var}(\bar X_N)=\frac1N\sum_{\lvert k\rvert<N}\big(1-\frac{\lvert k\rvert}{N}\big)C_{XX}[k]\approx\frac1N\sum_kC_{XX}[k]$.

*Derivation.*

$$
\operatorname{Var}\big(\langle X\rangle_T\big)=E\Big[\Big(\frac{1}{2T}\int_{-T}^{T}(X(t)-\mu_X)dt\Big)^2\Big]=\frac{1}{4T^2}\int_{-T}^{T}\!\!\int_{-T}^{T}C_{XX}(t-s)\,dt\,ds .
$$

This is the same double integral over a square that appeared in the proof of Wiener–Khinchin (§21.4), with $\omega=0$. Collapsing it along lines $t-s=\tau$ gives $\frac{1}{4T^2}\int_{-2T}^{2T}(2T-\lvert\tau\rvert)C_{XX}(\tau)\,d\tau$, which is the stated formula. If $\int\lvert C_{XX}\rvert<\infty$, the triangular weight tends to $1$ and the integral is bounded, so the variance behaves like $(1/2T)\int C_{XX}$ and tends to zero. $\square$

Two important interpretations:

- **Effective number of independent samples.** If $X$ had independent samples every $\Delta$ seconds, a record of length $2T$ would contain $2T/\Delta$ of them and the variance of the mean would be $\sigma_X^2\Delta/(2T)$. Comparing with the formula, a correlated process behaves as if it had independent samples spaced $\Delta_{\text{eff}}=\int C_{XX}\,d\tau/\sigma_X^2=2\tau_c$ apart (twice the correlation time). The **effective sample size** is $N_{\text{eff}}=2T/(2\tau_c)$.
- **Frequency-domain view.** Averaging over $2T$ is a lowpass filter of noise bandwidth $B_N=1/(4T)$ (Example 21.28 with $T\to2T$). The residual fluctuation is the noise power in that tiny band around DC, $2\cdot S_{CC}(0)\cdot B_N=S_{CC}(0)/(2T)$ — exactly the formula. Only the noise density at DC matters, which is why 1/f noise, whose density diverges at DC, defeats long averaging.

> [!example] Worked Example 21.31 — How long must we average?
> A sensor output is $X(t)=\mu+N(t)$, where $N$ has $C_{NN}(\tau)=\sigma^2e^{-\alpha\lvert\tau\rvert}$ with $\sigma=1$ mV and $\alpha=1000\ \mathrm{s^{-1}}$. We estimate $\mu$ by averaging a record of total length $2T=1$ s.
>
> **Step 1 — asymptotic formula.**
> $$
> \operatorname{Var}(\hat\mu)\approx\frac{1}{2T}\int_{-\infty}^{\infty}\sigma^2e^{-\alpha\lvert\tau\rvert}d\tau=\frac{1}{2T}\cdot\frac{2\sigma^2}{\alpha}=\frac{2\times10^{-6}}{1\times1000}=2\times10^{-9}\ \mathrm{V^2},
> $$
> so the standard error is $\sqrt{2\times10^{-9}}\approx44.7\ \mu\mathrm{V}$.
>
> **Step 2 — exact formula.** Evaluating the triangular-weighted integral exactly with $L=2T$ gives $\frac{2\sigma^2}{\alpha L}\Big[1-\frac{1-e^{-\alpha L}}{\alpha L}\Big]=2\times10^{-9}(1-0.001)=1.998\times10^{-9}\ \mathrm{V^2}$: the asymptotic formula is excellent once $L\gg1/\alpha$.
>
> **Step 3 — effective samples.** $\tau_c=1/\alpha=1$ ms, so $N_{\text{eff}}=L/(2\tau_c)=500$, and $\sigma/\sqrt{500}=44.7\ \mu$V. ✓ Sampling the record at 1 MHz (a million samples) does not help: the samples are highly correlated, and the error is set by the 500 effectively independent ones.
>
> **Step 4 — design.** To reach a standard error of $10\ \mu$V, we need $N_{\text{eff}}=(1\ \text{mV}/10\ \mu\text{V})^2=10^4$, i.e. $L=2\tau_cN_{\text{eff}}=20$ s of data.

> [!example] Worked Example 21.32 — Sample mean of an AR(1) sequence
> Take the AR(1) sequence of Example 21.20 ($a=0.9$, $\sigma_W^2=1$, $\sigma_X^2=5.263$) and $N=1000$ samples.
>
> **Step 1 — sum of covariances.** $\sum_kC_{XX}[k]=\sigma_X^2\sum_ka^{\lvert k\rvert}=\sigma_X^2\frac{1+a}{1-a}=5.263\times19=100$. (This equals $S_{XX}(e^{j0})=\sigma_W^2/(1-a)^2=100$, as it must.)
>
> **Step 2 — variance of the mean.** $\operatorname{Var}(\bar X_N)\approx100/1000=0.1$; the exact finite-$N$ sum gives $0.0991$. Standard error $0.316$.
>
> **Step 3 — comparison with i.i.d.** Treating the samples as independent would give $\sigma_X^2/N=0.00526$ — **19 times too small**, an error bar $\sqrt{19}=4.4$ times too optimistic. The effective sample size is $N_{\text{eff}}=N\frac{1-a}{1+a}=52.6$.
>
> This is a common error when applying the confidence intervals of Ch. 20 to time-series data: those formulas assume independent observations.

### Non-ergodic processes

> [!example] Worked Example 21.33 — Stationary processes that are not ergodic
> **(a) Random constant.** Let $X(t)=A$ for all $t$, where $A\sim\mathcal{N}(0,1)$. Then $\mu_X=0$ and $R_{XX}(\tau)=1$: the process is SSS. But each sample function is a constant, so $\langle x\rangle_T=A$, which is not $0$ in general. Here $C_{XX}(\tau)=1$ does not decay, and the variance formula gives $\operatorname{Var}\langle X\rangle_T=1\not\to0$. One record never reveals the ensemble mean. In the PSD this shows up as an impulse $2\pi\delta(\omega)$ in the spectrum of the zero-mean part — a random DC component.
>
> **(b) Random-phase sinusoid.** $A\cos(\omega_0t+\Theta)$ with fixed $A$ *is* ergodic in mean and autocorrelation: every sample function has time average $0$ and time autocorrelation $\tfrac{A^2}{2}\cos\omega_0\tau$ (Example 21.1, Step 5).
>
> **(c) Random amplitude and phase.** $A\cos(\omega_0t+\Theta)$ with *random* $A$ (independent of $\Theta$) is still mean-ergodic, but **not autocorrelation-ergodic**: the time autocorrelation of one record is $\tfrac{A^2}{2}\cos\omega_0\tau$, while the ensemble one is $\tfrac{E[A^2]}{2}\cos\omega_0\tau$. A power meter on one record reports that record's $A^2/2$.
>
> **(d) Random choice of source.** A switch, set once at random, connects either a 1 mV-rms or a 10 mV-rms noise source to the output. The resulting process is stationary, with ensemble power $\tfrac12(1+100)=50.5\ \mathrm{mV^2}$, but every individual record has time-average power $1$ or $100\ \mathrm{mV^2}$ — never $50.5$. Non-ergodicity typically means "the ensemble is a mixture of different sub-populations, and each record sees only one".

The general principle: a stationary process fails to be ergodic when some random quantity is fixed for the whole record (an amplitude, a DC offset, a choice of state). In engineering terms, ergodicity is a statement that **one device, observed for a long time, is representative of the population of devices**. This is plausible for the thermal noise of a resistor; it is false for, say, the offset voltage of op-amps, which varies from unit to unit but is constant in time for each unit. Unit-to-unit variation needs ensemble statistics (many devices), not long records.

### Autocorrelation ergodicity

Autocorrelation-ergodicity is mean-ergodicity of the product process $Z_\tau(t)=X(t+\tau)X(t)$, whose autocovariance involves **fourth** moments of $X$. For a zero-mean Gaussian process, the Isserlis (Wick) moment theorem gives $E[X_1X_2X_3X_4]=R_{12}R_{34}+R_{13}R_{24}+R_{14}R_{23}$, and one finds

$$
C_{Z_\tau Z_\tau}(s)=R_{XX}^2(s)+R_{XX}(s+\tau)R_{XX}(s-\tau).
$$

> [!theorem] Autocorrelation ergodicity of Gaussian processes
> A zero-mean stationary Gaussian process with $\int_{-\infty}^{\infty}R_{XX}^2(\tau)\,d\tau<\infty$ (equivalently, finite-energy correlation, which holds whenever its spectrum has no impulses and is square-integrable) is autocorrelation-ergodic. More generally, a stationary Gaussian process is ergodic (in every sense) if and only if its spectrum has no discrete (impulsive) part.

So for Gaussian noise without spectral lines — the typical receiver-noise case — time averages of any kind converge to ensemble averages, which is why the ergodic assumption is so rarely questioned in noise measurements.

### Estimating correlations and spectra from data

From $N$ samples $x[0],\dots,x[N-1]$ of a zero-mean sequence, the standard estimators are

$$
\hat R_{XX}[k]=\frac{1}{N}\sum_{n=0}^{N-1-k}x[n+k]x[n]\quad(0\le k<N),\qquad \hat R_{XX}[-k]=\hat R_{XX}[k],
$$

$$
\hat S_{XX}(e^{j\Omega})=\frac{1}{N}\Big\lvert\sum_{n=0}^{N-1}x[n]e^{-j\Omega n}\Big\rvert^2=\sum_{\lvert k\rvert<N}\hat R_{XX}[k]e^{-j\Omega k}\quad\text{(periodogram)}.
$$

Properties (all resting on an ergodicity assumption):

- **Bias of $\hat R$.** $E[\hat R_{XX}[k]]=\big(1-\frac{\lvert k\rvert}{N}\big)R_{XX}[k]$: biased (a triangular taper), but asymptotically unbiased, and guaranteed non-negative definite. Dividing by $N-\lvert k\rvert$ instead removes the bias but can produce invalid (non-positive-definite) estimates at large lags.
- **Bias of the periodogram.** $E[\hat S]=\frac{1}{2\pi}\big(S_{XX}*F_N\big)(e^{j\Omega})$, the true spectrum convolved with the Fejér kernel $F_N(\Omega)=\frac1N\frac{\sin^2(N\Omega/2)}{\sin^2(\Omega/2)}$: finite records smear spectral detail finer than about $2\pi/N$ rad/sample and leak power from strong peaks into their sidelobes. Windowing (Hann, Blackman) lowers the sidelobes at the expense of resolution.
- **Variance of the periodogram.** For Gaussian noise, at each frequency $\hat S\approx S\cdot\chi^2_2/2$ — an exponential random variable whose standard deviation equals its mean, **however large $N$ is**. The raw periodogram is not a consistent estimator: more data gives more frequency points, not smoother ones.
- **Welch's method.** Split the record into $K$ segments of length $L$, window each, compute $K$ periodograms and average them. The variance drops by about $K$ (more with 50% overlap), while the resolution coarsens to $\sim f_s/L$. This is what every FFT spectrum analyzer's "averaging" does.

> [!example] Worked Example 21.34 — Planning a Welch spectral estimate
> A noise signal is captured at $f_s=100$ kHz for 1 s ($N=10^5$ samples). We want a PSD estimate with about 100 Hz resolution and better than 10% relative standard deviation.
>
> **Step 1 — segment length.** Resolution $\approx f_s/L$; choose $L=1024$, giving $\Delta f=97.7$ Hz.
>
> **Step 2 — number of segments.** Non-overlapping: $K=\lfloor10^5/1024\rfloor=97$. Relative standard deviation of the averaged estimate $\approx1/\sqrt K=1/\sqrt{97}=10.2\%$ — just at the target. With 50% overlap and a Hann window, $K\approx194$ segments whose effective independence is about 0.9, bringing the relative error to roughly 7.5%.
>
> **Step 3 — trade-off.** Halving $L$ to 512 doubles $K$ (relative error $\div\sqrt2$) but doubles $\Delta f$ to 195 Hz. Resolution × (relative error)$^2$ × record length is roughly constant — a time–bandwidth uncertainty relation for spectral estimation: to get both fine resolution and low variance, you need a long record.

> [!warning] Common Pitfall — ergodicity cannot be checked from one record
> Ergodicity is a property relating one record to the ensemble; a single record, however long, cannot reveal whether other records would behave differently (Example 21.33(d)). It must be justified physically ("this noise is generated by the same mechanism in every realization, with no hidden parameters fixed per record") or checked with several independent records — for instance, several devices or several power-up cycles.

## 21.8 Applications to communications and measurement

### Motivation

The tools of §§21.1–21.7 combine into a short, repeatable procedure that settles most noise questions in circuits, instruments and communication links: model each noise source by its PSD, carry it through the system with $\lvert H\rvert^2$, add the PSDs of independent sources, integrate to get power, and — if the noise is Gaussian — convert powers into probabilities with $Q(\cdot)$. This section applies that procedure to the physical noise sources an EE meets most often (thermal and shot noise), to the figures of merit used to specify them (noise bandwidth, SNR, noise figure), and to the central design question of detection: which filter makes a known signal most visible in white noise? The aim is to show the mathematics at work; the device physics and communication-system design belong to other courses.

### Thermal (Johnson–Nyquist) noise

The random thermal motion of charge carriers in a resistor $R$ at absolute temperature $T$ produces a fluctuating open-circuit voltage. Because it is the sum of the effects of an enormous number of independent carrier collisions, each lasting about $10^{-14}$ s, it is to an excellent approximation zero-mean, **Gaussian** (central limit theorem) and **white** up to frequencies of order $kT/h\approx6$ THz at 300 K.

> [!theorem] Johnson–Nyquist formula
> The open-circuit noise voltage of a resistor $R$ in thermal equilibrium at temperature $T$ has
> $$
> G_V(f)=4kTR\ \ \mathrm{V^2/Hz}\ (\text{one-sided}),\qquad S_V(f)=2kTR\ \ (\text{two-sided, per Hz}),
> $$
> where $k=1.380649\times10^{-23}$ J/K is Boltzmann's constant. In a bandwidth $B$ (or noise bandwidth $B_N$):
> $$
> \overline{v_n^2}=4kTRB .
> $$
> Equivalently (Norton form), a noise current source of one-sided PSD $G_I=4kT/R$ in parallel with a noiseless $R$.

*Sketch of Nyquist's argument.* Connect two equal resistors $R$ by a lossless transmission line of characteristic impedance $R$ and length $L$; everything is at temperature $T$. Each resistor delivers noise power into the line, matched, so it is absorbed by the other. Now short both ends: the line becomes a resonator with modes at $f_m=mc/(2L)$, i.e. $2L\,B/c$ modes in a band $B$. By the equipartition theorem each mode holds average energy $kT$, so the energy in band $B$ is $kT\cdot2LB/c$; the line delivers it in a transit time $L/c$ in each direction, giving power $kTB$ from each resistor into the other. A resistor with open-circuit noise voltage $v_n$ delivers $\overline{v_n^2}/(4R)$ into a matched load, so $\overline{v_n^2}/(4R)=kTB$, i.e. $\overline{v_n^2}=4kTRB$. $\square$

Two corollaries are used constantly:

- **Available noise power** $P_n=kTB$, independent of $R$: $kT_0=4.00\times10^{-21}$ W/Hz at the reference temperature $T_0=290$ K, i.e. **−174 dBm/Hz**.
- **Combining resistors.** Independent noise sources add in power. Two resistors in series: $4kT(R_1+R_2)B$ — the same as one resistor $R_1+R_2$. In parallel: $4kT(R_1\parallel R_2)B$. More generally, the noise of any passive network of resistors, capacitors and inductors at a uniform temperature is $4kT\operatorname{Re}\{Z(f)\}$ per hertz, where $Z$ is the network's impedance (the fluctuation–dissipation theorem).

> [!example] Worked Example 21.35 — Resistor noise in audio and RF bands
> **(a) Audio.** $R=10$ kΩ, $T=300$ K, $B=20$ kHz:
> $$
> v_{n,\text{rms}}=\sqrt{4kTRB}=\sqrt{4\times1.381\times10^{-23}\times300\times10^4\times2\times10^4}=\sqrt{3.31\times10^{-12}}\approx1.82\ \mu\mathrm{V}.
> $$
> The often-quoted figure $\sqrt{4kTR}\approx12.9\ \mathrm{nV/\sqrt{Hz}}$ for 10 kΩ is the same result per root-hertz ($12.9\ \text{nV}\times\sqrt{20000}=1.82\ \mu$V ✓).
>
> **(b) RF.** $R=50\ \Omega$, $T=290$ K, $B=1$ MHz: $\overline{v_n^2}=8.0\times10^{-13}\ \mathrm{V^2}$, $v_{n,\text{rms}}=0.895\ \mu$V open-circuit. The available power is $kTB=4.0\times10^{-15}$ W $=-114$ dBm — the noise floor of a 1 MHz-wide receiver before any amplifier noise is added.
>
> **(c) Parallel combination.** Two 1 kΩ resistors in parallel at 300 K, $B=10$ kHz: equivalent $500\ \Omega$, $v_{n,\text{rms}}=\sqrt{4kT\cdot500\cdot10^4}=0.288\ \mu$V. Check via Norton sources: each contributes current density $4kT/1000$; the two independent currents add in power to $8kT/1000$, and across $500\ \Omega$ give voltage density $8kT/1000\times500^2=4kT\times500$. ✓

> [!example] Worked Example 21.36 — $kT/C$ noise
> Model a resistor $R$ (with its noise source) connected across a capacitor $C$. The capacitor voltage is the RC-lowpass output of white noise with one-sided density $4kTR$, i.e. $N_0/2=2kTR$ two-sided.
>
> **Step 1.** By Example 21.18, $\overline{v_C^2}=\dfrac{N_0}{4RC}=\dfrac{4kTR}{4RC}=\dfrac{kT}{C}$.
>
> **Step 2 — the resistance cancels.** A larger $R$ produces more noise density but a proportionally narrower noise bandwidth $B_N=1/(4RC)$. Equivalently, by equipartition, the capacitor's mean stored energy $\tfrac12C\overline{v_C^2}$ equals $\tfrac12kT$.
>
> **Step 3 — numbers.** For $C=1$ pF at 300 K, $v_{\text{rms}}=\sqrt{4.14\times10^{-21}/10^{-12}}\approx64\ \mu$V. For a 12-bit ADC with 1 V full scale (LSB $=244\ \mu$V), the quantization noise is $244/\sqrt{12}=70\ \mu$V rms, so a 1 pF sampling capacitor already contributes noise comparable to the quantization error. This sets a floor on the noise of sample-and-hold circuits, switched-capacitor filters and CMOS image-sensor pixels, and is why low-noise ADC front ends use large sampling capacitors.

### Shot noise

Current in a diode or photodetector consists of discrete charges $q$ crossing a barrier at random (Poisson) times. Modelling the current as $I(t)=q\sum_k\delta(t-t_k)$ with Poisson arrivals of rate $\lambda=\bar I/q$ (Example 21.3), one finds that its fluctuation about the mean $\bar I$ is white with one-sided density

$$
G_I(f)=2q\bar I\qquad\Longrightarrow\qquad\overline{i_n^2}=2q\bar IB .
$$

*Sketch.* For a Poisson impulse train of rate $\lambda$, the autocovariance is $\lambda q^2\delta(\tau)$ (increments over disjoint intervals are independent), so the two-sided PSD of the fluctuation is $\lambda q^2=q\bar I$, and the one-sided density is $2q\bar I$.

For example, a photodiode current of $\bar I=1\ \mu$A in a 1 MHz bandwidth has $\sqrt{2\times1.602\times10^{-19}\times10^{-6}\times10^{6}}=0.566$ nA rms of shot noise — a relative fluctuation of $5.7\times10^{-4}$, the $1/\sqrt{\lambda t}$ law of Example 21.3 in another form.

### Signal-to-noise ratio, noise figure and averaging

> [!definition] Signal-to-noise ratio
> For a signal component with power $P_s$ and an additive noise component with power $P_n$ at the same point, $\text{SNR}=P_s/P_n$, usually quoted as $10\log_{10}(\text{SNR})$ dB. If signal $s(t)$ (deterministic or random, with PSD $S_{ss}$) and independent noise with PSD $S_{NN}$ pass through $H$,
> $$
> \text{SNR}_{\text{out}}=\frac{\frac{1}{2\pi}\int\lvert H\rvert^2S_{ss}\,d\omega}{\frac{1}{2\pi}\int\lvert H\rvert^2S_{NN}\,d\omega}.
> $$

> [!definition] Noise factor and noise figure
> For a two-port with available gain $G$, driven by a source at $T_0=290$ K, the **noise factor** is
> $$
> F=\frac{\text{SNR}_{\text{in}}}{\text{SNR}_{\text{out}}}=1+\frac{T_e}{T_0}\ \ge1,
> $$
> where $T_e$ is the **equivalent input noise temperature**: the amplifier's internal noise referred to its input as if it were extra thermal noise $kT_eB$. The **noise figure** is $NF=10\log_{10}F$ dB.

For a cascade of stages with gains $G_1,G_2,\dots$ and noise factors $F_1,F_2,\dots$, refer each stage's added noise to the input of the chain by dividing by the gain preceding it. Because the noise contributions are independent, their powers add, and one obtains **Friis's formula**:

$$
F=F_1+\frac{F_2-1}{G_1}+\frac{F_3-1}{G_1G_2}+\cdots
$$

> [!example] Worked Example 21.37 — Receiver noise floor and cascade
> A receiver front end has an LNA with $NF_1=1.5$ dB ($F_1=1.413$) and gain $15$ dB ($G_1=31.6$), followed by a mixer with $NF_2=10$ dB ($F_2=10$).
>
> **Step 1 — cascade.** $F=1.413+\dfrac{10-1}{31.6}=1.413+0.285=1.697$, i.e. $NF=2.30$ dB. The noisy mixer costs only 0.8 dB because the LNA's gain divides its contribution.
>
> **Step 2 — sensitivity.** The noise floor referred to the input in a noise bandwidth $B_N=200$ kHz is $-174+10\log_{10}(2\times10^5)+NF=-174+53.0+2.3=-118.7$ dBm. If the demodulator needs $\text{SNR}=10$ dB, the minimum detectable signal is $-108.7$ dBm.
>
> **Step 3 — reversing the order.** Put the mixer first; suppose it is passive with a conversion loss of 7 dB ($G=0.2$). Then $F=10+\dfrac{1.413-1}{0.2}=12.07$, i.e. $NF=10.8$ dB — 8.5 dB worse, and the sensitivity drops to about $-100$ dBm. The first stage dominates; this is the mathematical reason the lowest-noise, highest-gain stage goes first.

> [!theorem] SNR gain from averaging
> If a deterministic waveform $s(t)$ is observed $M$ times in independent zero-mean noise of variance $\sigma^2$ (per sample), and the $M$ records are averaged point by point (synchronous or **coherent averaging**), the signal is unchanged while the noise variance becomes $\sigma^2/M$. The SNR improves by a factor $M$ (10 dB per decade of $M$). If the noise in successive records is correlated with correlation coefficients $\rho_{ik}$, the improvement is $M^2/\sum_{i,k}\rho_{ik}$ instead.

*Proof.* The average of the noise samples has variance $\frac{1}{M^2}\sum_{i,k}\operatorname{Cov}(N_i,N_k)=\frac{\sigma^2}{M^2}\sum_{i,k}\rho_{ik}$, which is $\sigma^2/M$ for uncorrelated records. $\square$

This is how oscilloscopes, evoked-potential recorders and radar pulse integrators dig repetitive signals out of noise; Example 21.31 is the same principle for a constant signal.

### The matched filter

A pulse $s(t)$ of known shape, nonzero only on $[0,T]$, is received in additive white noise of two-sided PSD $N_0/2$. It is passed through a filter $h$ and the output is sampled at time $T_0$. Which $h$ makes the pulse most visible — i.e. maximizes the output SNR at $T_0$?

At the sampling instant the signal output is $y_s(T_0)=\int h(\alpha)s(T_0-\alpha)\,d\alpha$, and the noise output has variance $\frac{N_0}{2}\int h^2(\alpha)\,d\alpha$ (§21.6). So

$$
\text{SNR}(T_0)=\frac{\Big(\int h(\alpha)s(T_0-\alpha)\,d\alpha\Big)^2}{\frac{N_0}{2}\int h^2(\alpha)\,d\alpha}.
$$

> [!theorem] Matched filter
> The output SNR at $T_0$ is maximized by
> $$
> h_{\text{opt}}(t)=c\,s(T_0-t)\qquad(\text{any constant }c\neq0),\qquad H_{\text{opt}}(\omega)=c\,S^*(\omega)e^{-j\omega T_0},
> $$
> and the maximum is
> $$
> \text{SNR}_{\max}=\frac{2E}{N_0},\qquad E=\int s^2(t)\,dt ,
> $$
> which depends only on the pulse **energy**, not its shape. Choosing $T_0=T$ makes $h_{\text{opt}}$ causal.

*Proof.* By the Cauchy–Schwarz inequality for integrals,

$$
\Big(\int h(\alpha)s(T_0-\alpha)\,d\alpha\Big)^2\le\int h^2(\alpha)\,d\alpha\cdot\int s^2(T_0-\alpha)\,d\alpha=E\int h^2(\alpha)\,d\alpha,
$$

with equality iff $h(\alpha)\propto s(T_0-\alpha)$. Substituting, $\text{SNR}\le E/(N_0/2)=2E/N_0$. $\square$

The filter's impulse response is the pulse **reversed in time**; equivalently, its output at $T_0$ is the **correlation** $\int r(t)s(t)\,dt$ of the received waveform with the known pulse — a correlator receiver. Its frequency response emphasizes the frequencies where the signal is strong ($\lvert H\rvert\propto\lvert S\rvert$) and aligns all their phases to peak simultaneously at $T_0$. For coloured noise of PSD $S_{NN}$, the same argument applied after a whitening filter gives $H_{\text{opt}}\propto S^*(\omega)e^{-j\omega T_0}/S_{NN}(\omega)$.

> [!example] Worked Example 21.38 — Matched filter for a rectangular pulse, and what a simple RC filter loses
> **Pulse.** $s(t)=A$ on $[0,T]$, $A=1$ mV, $T=1$ ms; noise $N_0/2=10^{-10}\ \mathrm{V^2/Hz}$.
>
> **Step 1 — matched filter.** $h(t)=s(T-t)$ is again a rectangle on $[0,T]$: the **integrate-and-dump** filter. $E=A^2T=10^{-9}\ \mathrm{V^2s}$ and
> $$
> \text{SNR}_{\max}=\frac{2E}{N_0}=\frac{2\times10^{-9}}{2\times10^{-10}}=10\ \ (10\ \text{dB}).
> $$
>
> **Step 2 — detection.** For antipodal signalling ($\pm s(t)$) with a threshold at zero and Gaussian noise, the error probability is $P_e=Q\big(\sqrt{\text{SNR}_{\max}}\big)=Q\big(\sqrt{2E/N_0}\big)=Q(3.16)=7.8\times10^{-4}$. Raising $E/N_0$ to $9.6$ dB gives the textbook $P_e\approx10^{-5}$.
>
> **Step 3 — RC lowpass instead.** Use $H=1/(1+j\omega RC)$ and sample at $T$. The signal output is $A(1-e^{-T/RC})$; the noise variance is $N_0/(4RC)$ (Example 21.18). With $x=T/RC$,
> $$
> \text{SNR}_{RC}=\frac{A^2(1-e^{-x})^2}{N_0/(4RC)}=\frac{2E}{N_0}\cdot\frac{2(1-e^{-x})^2}{x}.
> $$
> The factor $g(x)=2(1-e^{-x})^2/x$ is maximized numerically at $x=1.26$ (i.e. $RC=0.80\,T$), where $g=0.815$: the best RC filter is only $0.89$ dB worse than the matched filter. Too narrow a filter ($x\to0$) cuts signal; too wide ($x\to\infty$) admits noise.

```desmos-graph
left=0; right=6; bottom=0; top=1.1
---
y=\frac{2(1-e^{-x})^{2}}{x}|blue
y=1|red|dashed
(1.256,0.815)|label:best RC, −0.89 dB
```

*SNR of an RC lowpass relative to the matched filter for a rectangular pulse, against $x=T/RC$ (blue). The matched filter (red, dashed) is the upper bound set by Cauchy–Schwarz; the RC filter comes within 0.89 dB of it at $T/RC\approx1.26$.*

> [!tip] Engineering Note — estimation filters
> The matched filter maximizes SNR for *detecting* a known pulse. When instead a random signal $S(t)$ must be *estimated* from $X=S+N$, the filter minimizing $E[(S-\hat S)^2]$ is the (non-causal) **Wiener filter** $H(\omega)=S_{SX}(\omega)/S_{XX}(\omega)=S_{SS}/(S_{SS}+S_{NN})$ for independent $S,N$: it passes frequencies where the signal dominates and attenuates those where the noise does. Its recursive, time-domain counterpart for state-space models is the Kalman filter. Both are built entirely from the correlation functions and spectra of this chapter.

> [!tip] Engineering Note — the noise-calculation workflow
> Nearly every noise calculation in circuits and communications follows the same steps: (1) model each noise source by its PSD (thermal $4kTR$, shot $2qI$, amplifier $e_n^2$, $i_n^2$, flicker $\propto1/f$); (2) transfer each to the output with its own $\lvert H_i(\omega)\rvert^2$; (3) add the PSDs of independent sources; (4) integrate (or multiply by a noise bandwidth) to get power; (5) refer to the input by dividing by the signal gain if an input-referred figure is wanted; (6) if Gaussian, convert to probabilities with $Q(\cdot)$; (7) if the result comes from a measurement, check the averaging time with the ergodic formulas of §21.7.

> [!warning] Common Pitfall — adding RMS values
> Independent noise sources add **in power** (mean-square), never in RMS: $1\ \mu$V and $1\ \mu$V rms of independent noise give $1.41\ \mu$V, not $2\ \mu$V. Only fully correlated sources (the same noise reaching the output by two paths) add in amplitude — and then with their relative phase, possibly cancelling. Also remember that the PSD conventions (one-sided $G$, two-sided $S$, per hertz or per rad/s) must be the same for all terms before adding.

## Summary

| Concept | Key formula |
|---|---|
| Mean, autocorrelation | $\mu_X(t)=E[X(t)]$, $R_{XX}(t_1,t_2)=E[X(t_1)X(t_2)]$ |
| Autocovariance | $C_{XX}(t_1,t_2)=R_{XX}(t_1,t_2)-\mu_X(t_1)\mu_X(t_2)$ |
| Time averages | $\langle x\rangle=\lim\frac{1}{2T}\int_{-T}^{T}x\,dt$, $\langle x(t+\tau)x(t)\rangle$ |
| Independent increments | $C_{XX}(t_1,t_2)=\sigma_X^2(\min(t_1,t_2))$; random walk $\min(n_1,n_2)$; Poisson $\lambda\min(t_1,t_2)$ |
| WSS | $\mu_X$ constant, $R_{XX}(\tau)=E[X(t+\tau)X(t)]$ |
| Cyclostationary → WSS | $\bar R_{XX}(\tau)=\frac{1}{T_0}\int_0^{T_0}R_{XX}(t+\tau,t)\,dt$ |
| Autocorrelation properties | $R_{XX}(0)=E[X^2]$, even, $\lvert R_{XX}(\tau)\rvert\le R_{XX}(0)$, non-negative definite |
| Derivative | $R_{X'X'}(\tau)=-R_{XX}''(\tau)$; $E[X'X]=0$ |
| Cross-correlation | $R_{XY}(\tau)=E[X(t+\tau)Y(t)]=R_{YX}(-\tau)$ |
| Wiener–Khinchin | $S_{XX}(\omega)=\int R_{XX}(\tau)e^{-j\omega\tau}d\tau$ |
| Power | $E[X^2]=\frac{1}{2\pi}\int S_{XX}(\omega)\,d\omega=\int S_X(f)\,df=\int_0^\infty G_X(f)\,df$ |
| Exponential pair | $\sigma^2e^{-\alpha\lvert\tau\rvert}\leftrightarrow 2\alpha\sigma^2/(\alpha^2+\omega^2)$ |
| NRZ data | $A^2(1-\lvert\tau\rvert/T)\leftrightarrow A^2T\,\sin^2(\omega T/2)/(\omega T/2)^2$ |
| LTI output | $\mu_Y=\mu_XH(0)$, $S_{YX}=HS_{XX}$, $S_{YY}=\lvert H\rvert^2S_{XX}$ |
| Discrete WK, LTI | $S_{XX}(e^{j\Omega})=\sum_kR_{XX}[k]e^{-j\Omega k}$; $S_{YY}=\lvert H(e^{j\Omega})\rvert^2S_{XX}$ |
| AR(1) | $R_{XX}[k]=\frac{\sigma_W^2}{1-a^2}a^{\lvert k\rvert}$, $S_{XX}=\frac{\sigma_W^2}{1-2a\cos\Omega+a^2}$ |
| Gaussian conditional | $X(t+\tau)\mid X(t)=x\sim\mathcal{N}(\rho x,\ \sigma^2(1-\rho^2))$ (zero mean) |
| Rice's formula | $\nu_u=\frac{1}{\pi}\sqrt{-R''(0)/R(0)}\,e^{-u^2/2R(0)}$ |
| Narrowband envelope | Rayleigh: $P\{R>r\}=e^{-r^2/2\sigma^2}$ |
| White noise | $S=N_0/2$, $R=\frac{N_0}{2}\delta(\tau)$; filtered: $E[Y^2]=\frac{N_0}{2}\int h^2dt$ |
| Band-limited white noise | $R=N_0B\,\dfrac{\sin 2\pi B\tau}{2\pi B\tau}$, power $N_0B$ |
| Noise-equivalent bandwidth | $B_N=\frac{1}{\lvert H(0)\rvert^2}\int_0^\infty\lvert H(2\pi f)\rvert^2df$; RC: $1/(4RC)$; averager: $1/(2T)$ |
| Variance of time average | $\operatorname{Var}\langle X\rangle_T\approx\frac{1}{2T}\int C_{XX}(\tau)d\tau$; sequence: $\frac1N\sum_kC_{XX}[k]$ |
| Thermal noise | $\overline{v_n^2}=4kTRB$; available $kTB$ (−174 dBm/Hz); on a capacitor $kT/C$ |
| Shot noise | $\overline{i_n^2}=2q\bar IB$ |
| Friis | $F=F_1+(F_2-1)/G_1+\cdots$ |
| Matched filter | $h(t)=s(T_0-t)$, $\text{SNR}_{\max}=2E/N_0$ |

## Further reading

- A. Papoulis and S. U. Pillai, *Probability, Random Variables, and Stochastic Processes*, 4th ed., McGraw-Hill — the standard reference for stationarity, ergodicity and spectral theory.
- A. Leon-Garcia, *Probability, Statistics, and Random Processes for Electrical Engineering*, 3rd ed., Pearson — gentler, with many EE examples, including AR models and estimation.
- P. Z. Peebles, *Probability, Random Variables, and Random Signal Principles*, 4th ed., McGraw-Hill — concise treatment of noise in linear systems, noise bandwidth and noise figure.
- R. G. Gallager, *Stochastic Processes: Theory for Applications*, Cambridge University Press — rigorous treatment of Gaussian and Poisson processes.
