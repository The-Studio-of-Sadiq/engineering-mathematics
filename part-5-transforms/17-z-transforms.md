# Chapter 17 — Z-Transforms

[← Previous](16-fourier-transforms.md) · [Table of Contents](../00-front/table-of-contents.md) · [Next →](18-complex-analysis.md)

> **Prerequisites:** Ch. 2 (complex numbers), Ch. 4–5 (sequences, series and integral calculus), Ch. 14 (Laplace transforms), Ch. 16 (Fourier transforms and sampling)
> **Learning outcomes:** After this chapter you can
> - describe discrete-time signals as sequences, manipulate them with shifts, folds and convolution sums, and classify them (energy/power, periodic, causal);
> - compute bilateral and unilateral Z-transforms of standard sequences and apply the transform properties instead of summing series;
> - determine and sketch the region of convergence (ROC) for right-sided, left-sided, two-sided and finite sequences, and read causality and stability from it;
> - invert Z-transforms by partial fractions, power-series (long) division and the residue theorem;
> - solve linear constant-coefficient difference equations with initial conditions using the unilateral transform;
> - represent discrete-time LTI systems by their system function $H(z)$, pole–zero plots and block diagrams, and analyse their stability and frequency response $H(e^{j\omega})$, relating them to continuous-time systems through $z=e^{sT}$.

The Laplace transform of [Ch. 14](14-laplace-transforms.md) turns linear differential equations into algebra; the Fourier transform of [Ch. 16](16-fourier-transforms.md) describes signals as superpositions of sinusoids. Digital electronics, however, does not deal in functions of a continuous variable $t$. An analogue-to-digital converter delivers a list of numbers $x[0], x[1], x[2],\dots$, one per clock tick; a digital filter, a microcontroller control loop or a sigma-delta modulator computes each new output from a handful of past inputs and outputs. The governing equations are *difference* equations rather than differential equations, and the natural transform is the **Z-transform**, which converts a sequence into a function of a complex variable $z$ by using the sequence values as the coefficients of a power series in $z^{-1}$.

Everything in this chapter has a continuous-time twin: the Z-transform plays the role of the Laplace transform, the unit circle $\lvert z\rvert=1$ plays the role of the imaginary axis $s=j\omega$, the inside of the unit circle plays the role of the left half-plane, and the region of convergence is an annulus instead of a vertical strip. The dictionary between the two worlds is the single relation $z=e^{sT}$, derived in [§17.7](#177-discrete-time-system-analysis). Because $X(z)$ is a Laurent series, the deeper theory — why the inverse transform is a contour integral, why the ROC is an annulus, why residues compute the inverse — belongs to complex analysis and is completed in [Ch. 18](18-complex-analysis.md). Here we state what we need and use it.

## 17.1 Sequences and discrete-time functions

### Why sequences?

A continuous-time signal $x(t)$ sampled every $T$ seconds gives the numbers $x(nT)$, $n\in\mathbb{Z}$. Once in a processor the sampling interval is no longer visible: the processor sees only the index $n$. It is therefore convenient to define discrete-time signals abstractly, as functions on the integers, and to carry $T$ along only when we need to translate back to seconds and hertz.

> [!definition] Discrete-time signal (sequence)
> A **discrete-time signal** is a function $x:\mathbb{Z}\to\mathbb{C}$, written $x[n]$ with square brackets. The integer $n$ is the **sample index**. When $x[n]=x_a(nT)$ comes from sampling an analogue signal $x_a(t)$, $T$ is the **sampling interval** and $f_s=1/T$ the **sampling rate**.

Square brackets are a deliberate visual cue: $x[n]$ is defined only at integers, and $x[1.5]$ is meaningless (not zero). A finite sequence is often listed with an arrow, or here with a subscript marker, under the $n=0$ entry:

$$
x[n]=\{\,2,\ \underset{\uparrow}{1},\ -1,\ 3\,\}\quad\Longleftrightarrow\quad x[-1]=2,\ x[0]=1,\ x[1]=-1,\ x[2]=3,
$$

with $x[n]=0$ for all other $n$.

### Elementary sequences

> [!definition] Unit impulse, unit step and unit ramp
> $$
> \delta[n]=\begin{cases}1,&n=0\\0,&n\ne0\end{cases}\qquad
> u[n]=\begin{cases}1,&n\ge0\\0,&n<0\end{cases}\qquad
> r[n]=n\,u[n].
> $$

Unlike the Dirac delta of continuous time, the **Kronecker** (unit) impulse $\delta[n]$ is an ordinary function: it has the value 1 at one point. There is no limiting process, no distribution theory and no problem with $\delta[n]^2=\delta[n]$. The three sequences are linked by differences and running sums, the discrete analogues of derivatives and integrals:

$$
\delta[n]=u[n]-u[n-1],\qquad u[n]=\sum_{k=-\infty}^{n}\delta[k]=\sum_{k=0}^{\infty}\delta[n-k],\qquad r[n+1]-r[n]=u[n].
$$

**The sifting property** is the discrete statement that any sequence is a weighted sum of shifted impulses:

$$
x[n]=\sum_{k=-\infty}^{\infty}x[k]\,\delta[n-k].
$$

This identity, trivial as it looks, is the foundation of the convolution sum below.

> [!definition] Exponential and sinusoidal sequences
> The **complex exponential sequence** is $x[n]=A\,\alpha^n$ with $A,\alpha\in\mathbb{C}$. Writing $\alpha=r e^{j\omega_0}$,
> $$
> A\alpha^n=A\,r^n e^{j\omega_0 n}=A\,r^n\left(\cos\omega_0 n+j\sin\omega_0 n\right).
> $$
> For $r=1$ this is a **discrete-time phasor** of frequency $\omega_0$ (radians per sample); its real part $\cos(\omega_0 n+\theta)$ is a **sinusoidal sequence**.

The modulus $r=\lvert\alpha\rvert$ controls growth: $r<1$ gives a decaying envelope, $r=1$ a constant envelope, $r>1$ a growing one. The angle $\omega_0$ controls oscillation. For real $\alpha$, $0<\alpha<1$ gives monotone decay, and $-1<\alpha<0$ gives decay with alternating sign — something a real continuous exponential $e^{at}$ can never do, because $e^{aT}>0$. A negative real pole is the discrete-time signature of oscillation at exactly half the sampling rate, $\omega_0=\pi$.

```desmos-graph
left=-1; right=16; bottom=-1.2; top=1.2
---
(\left[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15\right],0.8^{\left[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15\right]})|blue
(\left[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15\right],(-0.8)^{\left[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15\right]})|red
y=0.8^{x}|blue|dashed|0<x<15
y=-0.8^{x}|red|dashed|0<x<15
```

*The sequences $0.8^n u[n]$ (blue) and $(-0.8)^n u[n]$ (red) sampled at integers; both decay inside the envelope $\pm0.8^n$ (dashed), but the negative base alternates in sign — oscillation at $\omega_0=\pi$.*

### Discrete-time frequency is periodic

Two features of discrete-time sinusoids have no continuous-time counterpart, and both follow from the fact that $n$ is an integer.

> [!theorem] Frequency is only defined modulo $2\pi$
> For every integer $k$,
> $$
> e^{j(\omega_0+2\pi k)n}=e^{j\omega_0 n}\,e^{j2\pi kn}=e^{j\omega_0 n}.
> $$
> Hence the sequences with frequencies $\omega_0$ and $\omega_0+2\pi k$ are *identical*. All distinct discrete-time frequencies lie in an interval of length $2\pi$, conventionally $-\pi<\omega\le\pi$ or $0\le\omega<2\pi$.

In consequence, "high frequency" in discrete time means $\omega$ near $\pm\pi$ (the sequence alternates sign as fast as possible), and "low frequency" means $\omega$ near $0$ or near any multiple of $2\pi$. The sampled version of $\cos(\Omega t)$ with sampling interval $T$ is $\cos(\Omega T n)$, so the discrete frequency is

$$
\omega=\Omega T=2\pi\frac{f}{f_s}\quad\text{rad/sample},
$$

and the periodicity of $\omega$ is precisely the aliasing of [§16.8](16-fourier-transforms.md#168-sampling-related-mathematical-foundations): analogue frequencies $f$ and $f+kf_s$ produce the same samples.

> [!theorem] Periodicity of a discrete sinusoid
> The sequence $e^{j\omega_0 n}$ is periodic with period $N$ (a positive integer) if and only if $\omega_0 N=2\pi m$ for some integer $m$, i.e. if and only if $\omega_0/2\pi=m/N$ is **rational**. The fundamental period is the smallest such $N$, obtained with $m/N$ in lowest terms.

*Proof.* $e^{j\omega_0(n+N)}=e^{j\omega_0 n}$ for all $n$ iff $e^{j\omega_0N}=1$ iff $\omega_0N\in2\pi\mathbb{Z}$. $\square$

So $\cos(0.3\pi n)$ has $\omega_0/2\pi=3/20$ and fundamental period $N=20$ (it completes $m=3$ cycles of the underlying cosine in 20 samples), whereas $\cos(n)$ — frequency 1 rad/sample — is *not periodic at all*, because $1/2\pi$ is irrational, even though its samples lie on a perfectly periodic continuous curve.

### Operations on sequences

| Operation | Definition | Remarks |
|---|---|---|
| Shift (delay by $n_0>0$) | $y[n]=x[n-n_0]$ | moves the sequence $n_0$ samples right |
| Advance | $y[n]=x[n+n_0]$ | moves left; non-causal if used in real time |
| Fold (time reversal) | $y[n]=x[-n]$ | mirror about $n=0$ |
| Decimation (downsampling by $M$) | $y[n]=x[Mn]$ | keeps every $M$th sample; information is lost |
| Expansion (upsampling by $L$) | $y[n]=x[n/L]$ if $L\mid n$, else $0$ | inserts $L-1$ zeros between samples |
| First difference | $\nabla x[n]=x[n]-x[n-1]$ | analogue of $dx/dt$ (backward) |
| Accumulation | $y[n]=\sum_{k=-\infty}^{n}x[k]$ | analogue of $\int_{-\infty}^t x\,d\tau$ |

When shift and fold are combined, the order matters: $x[-n+2]=x[-(n-2)]$ is $x[-n]$ shifted right by 2 — equivalently, $x[n]$ advanced by 2 and then folded. The safe method is to substitute: to sketch $y[n]=x[2-n]$, tabulate $n$, compute $2-n$, and read off $x$.

> [!definition] Classes of sequences
> - **Causal** (right-sided from 0): $x[n]=0$ for $n<0$. **Anticausal**: $x[n]=0$ for $n>0$.
> - **Right-sided**: $x[n]=0$ for $n<N_1$ for some $N_1$. **Left-sided**: $x[n]=0$ for $n>N_2$. **Two-sided**: neither. **Finite-length**: both right- and left-sided.
> - **Even**: $x[-n]=x[n]$; **odd**: $x[-n]=-x[n]$. Every sequence splits uniquely as $x_e[n]=\tfrac12(x[n]+x[-n])$ plus $x_o[n]=\tfrac12(x[n]-x[-n])$.
> - **Energy**: $E_x=\sum_n\lvert x[n]\rvert^2$. **Average power**: $P_x=\lim_{N\to\infty}\frac{1}{2N+1}\sum_{n=-N}^{N}\lvert x[n]\rvert^2$. An **energy signal** has $0<E_x<\infty$ (so $P_x=0$); a **power signal** has $0<P_x<\infty$ (so $E_x=\infty$).
> - **Absolutely summable**: $\sum_n\lvert x[n]\rvert<\infty$ (written $x\in\ell^1$). **Bounded**: $\lvert x[n]\rvert\le B$ for all $n$.

The right/left/two-sided classification looks bookkeeping-level here, but in [§17.3](#173-region-of-convergence) it will determine the *shape* of the region of convergence.

> [!example] Worked Example 17.1 — Energy and power of standard sequences
> Classify (a) $x_1[n]=(0.5)^n u[n]$, (b) $x_2[n]=u[n]$, (c) $x_3[n]=\cos(0.3\pi n)$, (d) $x_4[n]=n\,u[n]$.
>
> **(a)** $E=\sum_{n=0}^\infty 0.25^n=\dfrac{1}{1-0.25}=\dfrac43$. Finite, so an energy signal; $P=0$.
>
> **(b)** $E=\sum_{n\ge0}1=\infty$. For power, $\sum_{n=-N}^{N}u[n]^2=N+1$, so $P=\lim (N+1)/(2N+1)=\tfrac12$. A power signal.
>
> **(c)** $\cos^2\theta=\tfrac12(1+\cos2\theta)$, and over a full period of 20 samples $\sum\cos(0.6\pi n)=0$, so the average of $\cos^2(0.3\pi n)$ over one period is exactly $\tfrac12$. Hence $P=\tfrac12$: a power signal (as for any sinusoid of unit amplitude).
>
> **(d)** $\sum_{n=0}^{N}n^2\approx N^3/3$, so $P=\lim N^3/(3(2N+1))=\infty$. Neither an energy nor a power signal.
>
> *Check:* numerically, $\sum_{n=0}^{50}0.25^n=1.33333$ and $\frac{1}{40}\sum_{n=0}^{39}\cos^2(0.3\pi n)=0.5000$.

### Discrete-time systems and the convolution sum

A **discrete-time system** maps an input sequence $x[n]$ to an output sequence $y[n]=\mathcal{T}\{x[n]\}$. The properties of [Ch. 14](14-laplace-transforms.md) carry over verbatim:

- **Linear**: $\mathcal{T}\{a x_1+b x_2\}=a\mathcal{T}\{x_1\}+b\mathcal{T}\{x_2\}$.
- **Time-invariant** (shift-invariant): if $x[n]\mapsto y[n]$, then $x[n-n_0]\mapsto y[n-n_0]$ for every integer $n_0$.
- **Causal**: $y[n]$ depends only on $x[k]$ with $k\le n$.
- **BIBO stable**: every bounded input produces a bounded output.

For a linear time-invariant (LTI) system, define the **impulse response** $h[n]=\mathcal{T}\{\delta[n]\}$. Applying $\mathcal{T}$ to the sifting identity, linearity moves $\mathcal{T}$ inside the sum and time invariance turns $\mathcal{T}\{\delta[n-k]\}$ into $h[n-k]$:

> [!theorem] Convolution sum
> The output of an LTI system with impulse response $h[n]$ is
> $$
> y[n]=\sum_{k=-\infty}^{\infty}x[k]\,h[n-k]=\sum_{k=-\infty}^{\infty}h[k]\,x[n-k]\;\equiv\;(x*h)[n].
> $$
> Discrete convolution is commutative, associative and distributive over addition, and $\delta[n]$ is its identity: $x*\delta=x$, $x*\delta[n-n_0]=x[n-n_0]$.

*Derivation.* $y[n]=\mathcal{T}\big\{\sum_k x[k]\delta[n-k]\big\}=\sum_k x[k]\,\mathcal{T}\{\delta[n-k]\}=\sum_k x[k]h[n-k]$. The substitution $m=n-k$ gives the second form. $\square$

Unlike the convolution integral, the convolution sum can be evaluated exactly by hand for short sequences. Two practical methods:

1. **Flip-and-slide.** Fold $h[k]$ to $h[-k]$, shift it to $h[n-k]$, multiply pointwise with $x[k]$, and add.
2. **Polynomial multiplication.** For finite sequences starting at $n=0$, $y$ has the coefficients of the product of the polynomials $\sum x[n]w^n$ and $\sum h[n]w^n$. This is not a coincidence: it is the convolution theorem of the Z-transform in disguise (§17.2).

> [!example] Worked Example 17.2 — Convolution of finite sequences
> Compute $y=x*h$ for $x[n]=\{\underset{\uparrow}{1},2,3\}$ and $h[n]=\{\underset{\uparrow}{1},-1,2\}$.
>
> **Polynomial method.** $(1+2w+3w^2)(1-w+2w^2)=1+(2-1)w+(3-2+2)w^2+(-3+4)w^3+6w^4$, so
> $$
> y[n]=\{\underset{\uparrow}{1},\,1,\,3,\,1,\,6\}.
> $$
> **Check by flip-and-slide at $n=2$:** $y[2]=x[0]h[2]+x[1]h[1]+x[2]h[0]=2-2+3=3$. ✓
>
> **Length rule.** A length-$L_1$ sequence convolved with a length-$L_2$ sequence has length $L_1+L_2-1=5$; if the first starts at $N_1$ and the second at $N_2$, the result starts at $N_1+N_2$.

> [!example] Worked Example 17.3 — Step response of a first-order system by convolution
> An LTI system has $h[n]=a^n u[n]$ with $a\ne1$. Find its step response $s[n]=(u*h)[n]$.
>
> For $n<0$, every term of $\sum_k h[k]u[n-k]$ has either $k<0$ (so $h[k]=0$) or $n-k<0$ (so $u[n-k]=0$); thus $s[n]=0$. For $n\ge0$ the product is non-zero only for $0\le k\le n$:
> $$
> s[n]=\sum_{k=0}^{n}a^k=\frac{1-a^{n+1}}{1-a},\qquad n\ge0.
> $$
> If $\lvert a\rvert<1$, $s[n]\to1/(1-a)$, the **DC gain**. With $a=0.5$: $s[0..4]=1,\ 1.5,\ 1.75,\ 1.875,\ 1.9375\to2$.
>
> *Check:* direct recursion $s[n]=0.5\,s[n-1]+1$ from $s[-1]=0$ gives the same values.

> [!tip] Engineering Note — Sampled circuits are sequences
> When an ADC samples the capacitor voltage of an RC circuit with time constant $\tau$ every $T$ seconds, the free response $v(t)=V_0e^{-t/\tau}$ becomes $v[n]=V_0\left(e^{-T/\tau}\right)^n$: a geometric sequence with base $a=e^{-T/\tau}$. The ratio between consecutive samples, not the time constant, is what the processor sees. This is our first glimpse of $z=e^{sT}$: the continuous pole $s=-1/\tau$ becomes the discrete pole $a=e^{-T/\tau}$.

> [!warning] Common Pitfall — "Sampling a periodic signal gives a periodic sequence"
> Only if $f/f_s$ is rational. A 1 kHz tone sampled at $f_s=44.1$ kHz gives $\omega_0/2\pi=10/441$, periodic with $N=441$ samples; a tone at $1000\sqrt2$ Hz sampled at the same rate gives a sequence that never repeats exactly. Algorithms that assume an integer number of samples per period (e.g. some RMS meters) then show leakage errors.

## 17.2 Definition and properties

### Motivation: eigenfunctions of LTI systems

Feed the geometric sequence $x[n]=z^n$, with $z$ any non-zero complex number, into an LTI system with impulse response $h[n]$. The convolution sum gives

$$
y[n]=\sum_{k}h[k]\,z^{\,n-k}=z^n\sum_{k}h[k]\,z^{-k}=H(z)\,z^n,
$$

provided the sum converges. The output is the input multiplied by a complex constant $H(z)$: geometric sequences are **eigenfunctions** of every LTI system, and $H(z)$ is the corresponding eigenvalue. Exactly as $e^{st}$ led to the Laplace transform in continuous time, $z^n$ leads to the Z-transform. Decomposing an input into geometric sequences and multiplying each by $H(z)$ is the whole strategy of this chapter.

### Definition

> [!definition] Bilateral Z-transform
> The (bilateral, or two-sided) **Z-transform** of a sequence $x[n]$ is the function of the complex variable $z$
> $$
> X(z)=\mathcal{Z}\{x[n]\}=\sum_{n=-\infty}^{\infty}x[n]\,z^{-n},
> $$
> defined on the set of $z$ for which the series converges absolutely. That set is the **region of convergence** (ROC). We write $x[n]\overset{\mathcal{Z}}{\longleftrightarrow}X(z)$, ROC $=R_x$.

> [!definition] Unilateral Z-transform
> The **unilateral** (one-sided) Z-transform is
> $$
> \mathcal{X}(z)=\mathcal{Z}_u\{x[n]\}=\sum_{n=0}^{\infty}x[n]\,z^{-n}.
> $$
> It ignores $x[n]$ for $n<0$. For causal sequences the two transforms coincide. The unilateral form is the tool for difference equations with initial conditions (§17.5), just as the unilateral Laplace transform is for differential equations.

Writing $z=re^{j\omega}$,

$$
X(re^{j\omega})=\sum_n\left(x[n]\,r^{-n}\right)e^{-j\omega n},
$$

so the Z-transform on the circle of radius $r$ is the discrete-time Fourier transform (DTFT) of the exponentially weighted sequence $x[n]r^{-n}$. The weighting $r^{-n}$ plays the role of the factor $e^{-\sigma t}$ in the Laplace transform: choosing $r$ large tames sequences that grow as $n\to+\infty$; choosing $r$ small tames those that grow as $n\to-\infty$. On the unit circle, $r=1$, the Z-transform *is* the DTFT:

$$
X(e^{j\omega})=\sum_n x[n]\,e^{-j\omega n}\quad(\text{if the unit circle lies in the ROC}).
$$

The notation $X(e^{j\omega})$ — rather than $X(\omega)$ — is standard and deliberate: it reminds us that the DTFT is $X(z)$ evaluated at $z=e^{j\omega}$ and that it is $2\pi$-periodic in $\omega$.

### The fundamental pair: the geometric sequence

> [!example] Worked Example 17.4 — Transform of $a^n u[n]$
> Let $x[n]=a^n u[n]$. Then
> $$
> X(z)=\sum_{n=0}^{\infty}a^nz^{-n}=\sum_{n=0}^{\infty}\left(az^{-1}\right)^n=\frac{1}{1-az^{-1}}=\frac{z}{z-a},
> $$
> by the geometric series, which converges absolutely if and only if $\lvert az^{-1}\rvert<1$, i.e. $\lvert z\rvert>\lvert a\rvert$.
>
> **ROC:** the exterior of the circle of radius $\lvert a\rvert$. $X(z)$ has a **zero** at $z=0$ and a **pole** at $z=a$; the pole lies on the inner boundary of the ROC.
>
> *Check* ($a=0.5$, $z=2$): $\sum_{n=0}^{60}(0.25)^n=1.33333=2/(2-0.5)$. ✓

> [!example] Worked Example 17.5 — Transform of $-a^n u[-n-1]$
> Let $x[n]=-a^n u[-n-1]$, which is non-zero only for $n\le-1$. Substituting $m=-n$,
> $$
> X(z)=-\sum_{n=-\infty}^{-1}a^nz^{-n}=-\sum_{m=1}^{\infty}\left(a^{-1}z\right)^m=-\frac{a^{-1}z}{1-a^{-1}z}=\frac{z}{z-a},
> $$
> converging if and only if $\lvert a^{-1}z\rvert<1$, i.e. $\lvert z\rvert<\lvert a\rvert$.
>
> **The same algebraic expression** $z/(z-a)$ as in Example 17.4, but with ROC the *interior* of the circle $\lvert z\rvert=\lvert a\rvert$.

These two examples carry the most important lesson of the chapter: **a Z-transform is the pair (formula, ROC)**. The formula alone does not determine the sequence. In practice the ROC is often implied by context — causal signals in real-time processing have exterior ROCs — but whenever two-sided or non-causal sequences are possible it must be stated.

### Table of standard pairs

Each pair below follows from the geometric series and the properties derived next.

| $x[n]$ | $X(z)$ | ROC |
|---|---|---|
| $\delta[n]$ | $1$ | all $z$ |
| $\delta[n-n_0]$, $n_0>0$ | $z^{-n_0}$ | $z\ne0$ |
| $\delta[n+n_0]$, $n_0>0$ | $z^{n_0}$ | all finite $z$ ($z\ne\infty$) |
| $u[n]$ | $\dfrac{1}{1-z^{-1}}=\dfrac{z}{z-1}$ | $\lvert z\rvert>1$ |
| $-u[-n-1]$ | $\dfrac{z}{z-1}$ | $\lvert z\rvert<1$ |
| $a^nu[n]$ | $\dfrac{z}{z-a}$ | $\lvert z\rvert>\lvert a\rvert$ |
| $-a^nu[-n-1]$ | $\dfrac{z}{z-a}$ | $\lvert z\rvert<\lvert a\rvert$ |
| $n\,u[n]$ | $\dfrac{z}{(z-1)^2}$ | $\lvert z\rvert>1$ |
| $n\,a^nu[n]$ | $\dfrac{az}{(z-a)^2}$ | $\lvert z\rvert>\lvert a\rvert$ |
| $\binom{n+m-1}{m-1}a^nu[n]$ | $\dfrac{z^m}{(z-a)^m}$ | $\lvert z\rvert>\lvert a\rvert$ |
| $\cos(\omega_0n)\,u[n]$ | $\dfrac{z(z-\cos\omega_0)}{z^2-2z\cos\omega_0+1}$ | $\lvert z\rvert>1$ |
| $\sin(\omega_0n)\,u[n]$ | $\dfrac{z\sin\omega_0}{z^2-2z\cos\omega_0+1}$ | $\lvert z\rvert>1$ |
| $r^n\cos(\omega_0n)\,u[n]$ | $\dfrac{z(z-r\cos\omega_0)}{z^2-2rz\cos\omega_0+r^2}$ | $\lvert z\rvert>r$ |
| $r^n\sin(\omega_0n)\,u[n]$ | $\dfrac{zr\sin\omega_0}{z^2-2rz\cos\omega_0+r^2}$ | $\lvert z\rvert>r$ |
| $a^n$ for $0\le n\le N-1$, else 0 | $\dfrac{1-a^Nz^{-N}}{1-az^{-1}}$ | $z\ne0$ |

> [!tip] Engineering Note — Write transforms in $z^{-1}$ for implementation, in $z$ for algebra
> DSP texts and filter-design software (for example the coefficient vectors `b`, `a` passed to MATLAB's `filter` or SciPy's `lfilter`) write rational transforms as ratios of polynomials in $z^{-1}$, because $z^{-1}$ is a one-sample delay that can be built from a register. Pole–zero and partial-fraction work is cleaner in positive powers of $z$. Converting is just multiplying numerator and denominator by the same power of $z$ — but keep track of the poles or zeros at $z=0$ that this creates or removes.

### Linearity

> [!theorem] Linearity
> If $x_1\leftrightarrow X_1$ (ROC $R_1$) and $x_2\leftrightarrow X_2$ (ROC $R_2$), then
> $$
> a\,x_1[n]+b\,x_2[n]\;\longleftrightarrow\;a\,X_1(z)+b\,X_2(z),\qquad\text{ROC}\supseteq R_1\cap R_2.
> $$

The ROC can be larger than the intersection when poles cancel. For example $a^nu[n]-a^nu[n-1]=\delta[n]$: each term has ROC $\lvert z\rvert>\lvert a\rvert$, but the difference $\frac{z}{z-a}-\frac{a}{z-a}=1$ converges everywhere.

Linearity immediately gives the sinusoid pairs. By Euler, $\cos(\omega_0 n)u[n]=\tfrac12\left(e^{j\omega_0n}+e^{-j\omega_0n}\right)u[n]$, so with $a=e^{\pm j\omega_0}$:

$$
\frac12\left(\frac{z}{z-e^{j\omega_0}}+\frac{z}{z-e^{-j\omega_0}}\right)=\frac{z}{2}\cdot\frac{2z-2\cos\omega_0}{z^2-2z\cos\omega_0+1}=\frac{z(z-\cos\omega_0)}{z^2-2z\cos\omega_0+1},
$$

with ROC $\lvert z\rvert>1$. The two poles $e^{\pm j\omega_0}$ sit *on* the unit circle at angle $\pm\omega_0$: a pole on the unit circle at angle $\omega_0$ is an undamped oscillation at frequency $\omega_0$.

### Time shifting

> [!theorem] Time shift (bilateral)
> $$
> x[n-n_0]\;\longleftrightarrow\;z^{-n_0}X(z),
> $$
> with the same ROC as $X$, except possibly for the addition or removal of $z=0$ or $z=\infty$.

*Proof.* $\sum_n x[n-n_0]z^{-n}=\sum_m x[m]z^{-(m+n_0)}=z^{-n_0}X(z)$ with $m=n-n_0$. $\square$

So a delay of one sample is multiplication by $z^{-1}$ — the **unit-delay operator** in a block diagram. For the *unilateral* transform the shift rule must account for samples that cross $n=0$; it is derived in §17.5 where it is needed.

### Multiplication by an exponential (scaling in $z$)

> [!theorem] Exponential multiplication
> $$
> a^n x[n]\;\longleftrightarrow\;X\!\left(\frac{z}{a}\right),\qquad\text{ROC}=\lvert a\rvert R_x.
> $$
> Every pole and zero $p$ of $X$ moves to $ap$: the pole–zero pattern is scaled by $\lvert a\rvert$ and rotated by $\angle a$.

*Proof.* $\sum_n a^nx[n]z^{-n}=\sum_n x[n](z/a)^{-n}$. $\square$

With $a=e^{j\omega_0}$ this is **modulation**: $e^{j\omega_0n}x[n]\leftrightarrow X(e^{-j\omega_0}z)$, a rotation of the pole–zero plot by $\omega_0$, i.e. a frequency shift of the DTFT. With $a=r$ real it gives the damped sinusoid pairs in the table from the undamped ones: replace $z$ by $z/r$ in the $\cos(\omega_0 n)u[n]$ pair and simplify.

### Differentiation in $z$

> [!theorem] Multiplication by $n$
> $$
> n\,x[n]\;\longleftrightarrow\;-z\frac{dX(z)}{dz},\qquad\text{same ROC (up to } 0,\infty).
> $$

*Proof.* Inside its ROC a power series may be differentiated term by term ([Ch. 18](18-complex-analysis.md)): $\dfrac{dX}{dz}=\sum_n(-n)x[n]z^{-n-1}$; multiply by $-z$. $\square$

Applied to $a^nu[n]$:

$$
-z\frac{d}{dz}\frac{z}{z-a}=-z\cdot\frac{(z-a)-z}{(z-a)^2}=\frac{az}{(z-a)^2},
$$

the $na^nu[n]$ entry in the table. Applying the rule repeatedly produces $n^2a^nu[n]$, and so on, and hence the double, triple, … pole pairs needed in partial-fraction inversion.

### Time reversal and conjugation

> [!theorem] Time reversal and conjugation
> $$
> x[-n]\;\longleftrightarrow\;X(1/z),\quad\text{ROC}=1/R_x;\qquad x^*[n]\;\longleftrightarrow\;X^*(z^*),\quad\text{ROC}=R_x.
> $$

Reversal maps a pole at $p$ to a pole at $1/p$ and turns an exterior ROC $\lvert z\rvert>r$ into an interior one $\lvert z\rvert<1/r$. The conjugation property implies that for **real** sequences $X(z^*)=X^*(z)$: poles and zeros are real or occur in complex-conjugate pairs, and on the unit circle $X(e^{-j\omega})=X^*(e^{j\omega})$ — the magnitude response is even in $\omega$ and the phase is odd.

### Convolution

> [!theorem] Convolution theorem
> $$
> (x*h)[n]\;\longleftrightarrow\;X(z)\,H(z),\qquad\text{ROC}\supseteq R_x\cap R_h.
> $$

*Proof.* Take $z\in R_x\cap R_h$, so both series converge absolutely and the order of summation may be exchanged:
$$
\sum_n\sum_k x[k]h[n-k]z^{-n}=\sum_k x[k]z^{-k}\sum_n h[n-k]z^{-(n-k)}=X(z)H(z).\ \square
$$

The polynomial-multiplication trick of Example 17.2 is this theorem with $w=z^{-1}$. More important, the theorem turns the output of an LTI system into a product, $Y(z)=H(z)X(z)$, which is the starting point of §17.6.

### Accumulation and first difference

Because $\sum_{k=-\infty}^nx[k]=(x*u)[n]$ and $\mathcal{Z}\{u\}=z/(z-1)$:

$$
\sum_{k=-\infty}^{n}x[k]\;\longleftrightarrow\;\frac{z}{z-1}X(z)=\frac{1}{1-z^{-1}}X(z),\qquad x[n]-x[n-1]\;\longleftrightarrow\;(1-z^{-1})X(z).
$$

The accumulator has a pole at $z=1$ — the discrete analogue of the integrator's pole at $s=0$ — and the first difference has a zero there, like the differentiator's zero at $s=0$.

### Parseval's relation

For two sequences whose ROCs both contain the unit circle, the convolution theorem applied to $x[n]$ and $y^*[-n]$, evaluated at $n=0$, gives

$$
\sum_{n=-\infty}^{\infty}x[n]\,y^*[n]=\frac{1}{2\pi}\int_{-\pi}^{\pi}X(e^{j\omega})\,Y^*(e^{j\omega})\,d\omega,
\qquad
E_x=\sum_n\lvert x[n]\rvert^2=\frac{1}{2\pi}\int_{-\pi}^{\pi}\lvert X(e^{j\omega})\rvert^2d\omega.
$$

It is the discrete-time counterpart of the energy theorem of [§16.7](16-fourier-transforms.md#167-energy-spectral-density); the only difference is that the frequency integral runs over one period, $2\pi$, rather than over the whole real line. For $x[n]=0.5^nu[n]$, $\lvert X(e^{j\omega})\rvert^2=1/(1.25-\cos\omega)$ and the integral evaluates to $\frac{1}{2\pi}\cdot\frac{2\pi}{\sqrt{1.25^2-1}}=\frac{1}{0.75}=\frac43$, matching Example 17.1(a).

### Initial and final value theorems

> [!theorem] Initial value theorem
> If $x[n]$ is causal, then
> $$
> x[0]=\lim_{z\to\infty}X(z).
> $$
> More generally $x[1]=\lim_{z\to\infty}z\left(X(z)-x[0]\right)$, and so on.

*Proof.* $X(z)=x[0]+x[1]z^{-1}+x[2]z^{-2}+\cdots$; every term but the first vanishes as $z\to\infty$. $\square$

> [!theorem] Final value theorem
> If $x[n]$ is causal and $(z-1)X(z)$ has all its poles strictly inside the unit circle, then
> $$
> \lim_{n\to\infty}x[n]=\lim_{z\to1}(z-1)X(z)=\lim_{z\to1}\left(1-z^{-1}\right)X(z).
> $$

*Proof sketch.* $\mathcal{Z}\{x[n]-x[n-1]\}=(1-z^{-1})X(z)$. Summing the differences, $x[N]=\sum_{n=0}^{N}(x[n]-x[n-1])$; letting $N\to\infty$, $\lim x[N]=\lim_{z\to1}\sum_{n\ge0}(x[n]-x[n-1])z^{-n}$, the interchange of limits being justified by the pole condition, which guarantees the differences decay geometrically. $\square$

The hypothesis is essential: $X(z)=z/(z+1)$ (the sequence $(-1)^n u[n]$) gives $\lim_{z\to1}(z-1)X(z)=0$, yet $(-1)^n$ has no limit. The pole of $(z-1)X(z)$ at $z=-1$ lies on the unit circle, so the theorem does not apply.

### Summary table of properties

| Property | Sequence | Transform | ROC |
|---|---|---|---|
| Linearity | $ax_1+bx_2$ | $aX_1+bX_2$ | $\supseteq R_1\cap R_2$ |
| Time shift | $x[n-n_0]$ | $z^{-n_0}X(z)$ | $R_x$ (except $0,\infty$) |
| Exponential scaling | $a^nx[n]$ | $X(z/a)$ | $\lvert a\rvert R_x$ |
| Modulation | $e^{j\omega_0n}x[n]$ | $X(e^{-j\omega_0}z)$ | $R_x$ |
| Multiplication by $n$ | $nx[n]$ | $-z\,X'(z)$ | $R_x$ |
| Time reversal | $x[-n]$ | $X(1/z)$ | $1/R_x$ |
| Conjugation | $x^*[n]$ | $X^*(z^*)$ | $R_x$ |
| Convolution | $x*h$ | $XH$ | $\supseteq R_x\cap R_h$ |
| Accumulation | $\sum_{k\le n}x[k]$ | $X/(1-z^{-1})$ | $\supseteq R_x\cap\{\lvert z\rvert>1\}$ |
| First difference | $x[n]-x[n-1]$ | $(1-z^{-1})X$ | $\supseteq R_x\cap\{\lvert z\rvert>0\}$ |
| Initial value (causal) | $x[0]$ | $\lim_{z\to\infty}X(z)$ | — |
| Final value (causal) | $x[\infty]$ | $\lim_{z\to1}(z-1)X(z)$ | poles of $(z-1)X$ inside $\lvert z\rvert=1$ |

### Worked examples

> [!example] Worked Example 17.6 — A finite sequence
> Find $X(z)$ for $x[n]=\{2,\underset{\uparrow}{1},-1,3\}$.
>
> Directly from the definition:
> $$
> X(z)=2z+1-z^{-1}+3z^{-2}=\frac{2z^3+z^2-z+3}{z^2}.
> $$
> **ROC.** The positive power $2z$ diverges as $z\to\infty$, and the negative powers diverge at $z=0$, so the ROC is $0<\lvert z\rvert<\infty$. A finite sequence always has an ROC of the entire plane, except possibly $0$ and/or $\infty$.
>
> *Check:* at $z=1$, $X(1)=2+1-1+3=5=\sum_n x[n]$, as it must be: the Z-transform at $z=1$ is the sum of the sequence.

> [!example] Worked Example 17.7 — A damped sinusoid via scaling
> Find the Z-transform of $x[n]=(0.8)^n\cos(\pi n/3)\,u[n]$.
>
> **Step 1.** Start from the cosine pair with $\cos(\pi/3)=\tfrac12$:
> $$
> \cos\left(\tfrac{\pi}{3}n\right)u[n]\longleftrightarrow\frac{z(z-0.5)}{z^2-z+1},\quad\lvert z\rvert>1.
> $$
> **Step 2.** Multiply by $0.8^n$, i.e. replace $z$ by $z/0.8$:
> $$
> X(z)=\frac{(z/0.8)(z/0.8-0.5)}{(z/0.8)^2-z/0.8+1}=\frac{z(z-0.4)}{z^2-0.8z+0.64},\quad\lvert z\rvert>0.8.
> $$
> **Poles** at $z=0.4\pm j0.4\sqrt3=0.8e^{\pm j\pi/3}$: radius $0.8$ (the damping), angle $\pm\pi/3$ (the oscillation frequency). **Zeros** at $0$ and $0.4$.
>
> *Check* at $z=2$: the formula gives $2(1.6)/(4-1.6+0.64)=3.2/3.04=1.052632$; summing $\sum_{n=0}^{80}0.8^n\cos(\pi n/3)\,2^{-n}$ numerically gives $1.052632$. ✓

> [!example] Worked Example 17.8 — Using several properties together
> Find $X(z)$ for $x[n]=n\,(0.5)^{n-2}u[n-2]$.
>
> **Step 1.** Let $g[n]=(0.5)^{n-2}u[n-2]$, the sequence $0.5^nu[n]$ delayed by 2. Then $G(z)=z^{-2}\cdot\dfrac{z}{z-0.5}=\dfrac{1}{z(z-0.5)}$, $\lvert z\rvert>0.5$.
>
> **Step 2.** $x[n]=n\,g[n]$, so
> $$
> X(z)=-z\,\frac{d}{dz}\left[\frac{1}{z^2-0.5z}\right]=-z\cdot\frac{-(2z-0.5)}{(z^2-0.5z)^2}=\frac{2z-0.5}{z(z-0.5)^2},\quad\lvert z\rvert>0.5.
> $$
> *Check by leading coefficients.* $x[0]=x[1]=0$, $x[2]=2$, $x[3]=3(0.5)=1.5$. Expanding $X(z)$ for large $z$: $X=2z^{-2}\left(1-\tfrac{0.25}{z}\right)\left(1+\tfrac{1}{z}+\cdots\right)=2z^{-2}+1.5z^{-3}+\cdots$ ✓.

> [!example] Worked Example 17.9 — Initial and final values without inversion
> A causal sequence has $Y(z)=\dfrac{z^2}{(z-1)(z-0.5)}$. Find $y[0]$ and $\lim_{n\to\infty}y[n]$.
>
> **Initial value.** $\lim_{z\to\infty}\dfrac{z^2}{z^2-1.5z+0.5}=1$, so $y[0]=1$.
>
> **Final value.** $(z-1)Y(z)=\dfrac{z^2}{z-0.5}$ has its only pole at $0.5$, inside the unit circle, so the theorem applies: $\lim_{z\to1}\dfrac{z^2}{z-0.5}=\dfrac{1}{0.5}=2$.
>
> This $Y$ is $\dfrac{z}{z-1}\cdot\dfrac{z}{z-0.5}$, the step response of $h[n]=0.5^nu[n]$ (Example 17.3), whose values $1, 1.5, 1.75,\dots\to2$ confirm both results.

> [!warning] Common Pitfall — Dropping the ROC
> Writing "$\mathcal{Z}\{x\}=z/(z-a)$" without an ROC is incomplete: Examples 17.4 and 17.5 show two different sequences with this transform. Tables in most texts list only causal sequences and silently assume the exterior ROC; when you invert, you must know whether that assumption is right for your problem.

## 17.3 Region of convergence

### Why the ROC deserves its own section

Section 17.2 showed that $z/(z-a)$ is the transform of two different sequences. The ROC is what tells them apart, and it carries physical information: whether the sequence is right- or left-sided, whether the system it describes is causal, and whether it is stable. In continuous time the Laplace ROC is a vertical strip $\sigma_1<\operatorname{Re}s<\sigma_2$; in discrete time it is an **annulus** $r_1<\lvert z\rvert<r_2$, because the convergence of $\sum x[n]z^{-n}$ depends only on $\lvert z\rvert$.

### The ROC is an annulus

Absolute convergence of $\sum_n x[n]z^{-n}$ means

$$
\sum_{n=-\infty}^{\infty}\lvert x[n]\rvert\,r^{-n}<\infty,\qquad r=\lvert z\rvert.
$$

The angle of $z$ does not appear, so if one point on a circle centred at the origin is in the ROC, the whole circle is. Split the sum into the **causal part** $n\ge0$ and the **anticausal part** $n<0$:

$$
\underbrace{\sum_{n=0}^{\infty}\lvert x[n]\rvert\,r^{-n}}_{\text{converges for } r>r_1}\;+\;\underbrace{\sum_{m=1}^{\infty}\lvert x[-m]\rvert\,r^{m}}_{\text{converges for } r<r_2}.
$$

The first is a power series in $1/r$; by the root test (Ch. 4) it converges for $1/r<1/r_1$, i.e. $r>r_1$, where

$$
r_1=\limsup_{n\to\infty}\lvert x[n]\rvert^{1/n}.
$$

The second is a power series in $r$ and converges for $r<r_2$, where $1/r_2=\limsup_{m\to\infty}\lvert x[-m]\rvert^{1/m}$. Both must converge, so:

> [!theorem] Shape of the region of convergence
> The ROC of a Z-transform is an annulus centred at the origin,
> $$
> R_x=\{\,z:\ r_1<\lvert z\rvert<r_2\,\},\qquad 0\le r_1,\ r_2\le\infty,
> $$
> possibly including one or both boundary circles, and possibly empty (when $r_1\ge r_2$). The rate of growth of $x[n]$ as $n\to+\infty$ sets the inner radius $r_1$; the rate of growth as $n\to-\infty$ sets the outer radius $r_2$.

On the open annulus $X(z)$ is a convergent Laurent series and is therefore analytic ([Ch. 18](18-complex-analysis.md)). An analytic function has no poles, which gives the next property.

### Properties of the ROC

> [!theorem] Properties of the ROC
> 1. **No poles inside.** The ROC contains no poles of $X(z)$. For rational $X(z)$ the ROC is bounded by poles (or extends to $0$ or $\infty$).
> 2. **Finite-length sequences:** the ROC is the whole $z$-plane, except possibly $z=0$ (if $x[n]\ne0$ for some $n>0$) and $z=\infty$ (if $x[n]\ne0$ for some $n<0$).
> 3. **Right-sided sequences:** the ROC is the **exterior** of a circle, $\lvert z\rvert>r_1$ — for rational $X$, outside the outermost pole. It includes $z=\infty$ if and only if the sequence is causal.
> 4. **Left-sided sequences:** the ROC is the **interior** of a circle, $\lvert z\rvert<r_2$ — for rational $X$, inside the innermost (non-zero) pole. It includes $z=0$ if and only if the sequence is anticausal ($x[n]=0$ for $n>0$).
> 5. **Two-sided sequences:** the ROC is a **ring** $r_1<\lvert z\rvert<r_2$ between two consecutive pole radii, or empty.
> 6. **Connectedness:** the ROC is a single connected region.

*Proof of 3.* If $x[n]=0$ for $n<N_1$, the anticausal part has finitely many terms, so $r_2=\infty$ (if $N_1\ge0$, even $z=\infty$ is allowed, because every term $x[n]z^{-n}$ with $n\ge0$ stays finite as $z\to\infty$). Only $r_1$ is constrained: the ROC is $\lvert z\rvert>r_1$. Since there can be no poles in it, every pole has $\lvert p\rvert\le r_1$, and the boundary must pass through the outermost pole. Property 4 is the mirror image. $\square$

### ROC sketches

A sketch of the ROC consists of the pole–zero plot (poles $\times$, zeros $\circ$), the unit circle (for reference) and the shaded convergence region. Consider the rational function with poles at $z=0.5$ and $z=2$,

$$
X(z)=\frac{z}{(z-0.5)(z-2)}.
$$

Its ROC must be an annulus bounded by pole radii, so there are exactly three possibilities.

```desmos-graph
left=-3; right=3; bottom=-3; top=3
---
x^2+y^2>4|blue
x^2+y^2=1|black|dashed
(0.5,0)|red|label:pole 0.5
(2,0)|red|label:pole 2
(0,0)|green|label:zero
```

*ROC $\lvert z\rvert>2$ (shaded): the right-sided (here causal) sequence $x[n]=\tfrac23\left(2^n-0.5^n\right)u[n]$. The ROC lies outside the outermost pole and excludes the unit circle (dashed), so the sequence grows without bound.*

```desmos-graph
left=-3; right=3; bottom=-3; top=3
---
0.25<x^2+y^2<4|blue
x^2+y^2=1|black|dashed
(0.5,0)|red|label:pole 0.5
(2,0)|red|label:pole 2
(0,0)|green|label:zero
```

*ROC $0.5<\lvert z\rvert<2$ (shaded ring): the two-sided sequence $x[n]=-\tfrac23\,0.5^n u[n]-\tfrac23\,2^n u[-n-1]$. The ring contains the unit circle, so this sequence is absolutely summable and has a DTFT.*

```desmos-graph
left=-3; right=3; bottom=-3; top=3
---
x^2+y^2<0.25|blue
x^2+y^2=1|black|dashed
(0.5,0)|red|label:pole 0.5
(2,0)|red|label:pole 2
(0,0)|green|label:zero
```

*ROC $\lvert z\rvert<0.5$ (shaded disc): the left-sided sequence $x[n]=\tfrac23\left(0.5^n-2^n\right)u[-n-1]$, which blows up as $n\to-\infty$.*

The three sequences are derived in Example 17.12 below; the point here is that one rational function has three inverse transforms, selected by the ROC.

### Causality and stability read from the ROC

When $X(z)=H(z)$ is the transform of an impulse response, two system properties translate into ROC statements.

> [!theorem] Causality and stability from the ROC
> Let $h[n]$ have transform $H(z)$ with ROC $R_h$.
> - The system is **causal** ($h[n]=0$ for $n<0$) if and only if $R_h$ is the exterior of a circle *including* $z=\infty$. For rational $H(z)$: the ROC is outside the outermost pole, *and* the numerator degree does not exceed the denominator degree (in positive powers of $z$).
> - The system is **BIBO stable** ($\sum\lvert h[n]\rvert<\infty$) if and only if $R_h$ **includes the unit circle** $\lvert z\rvert=1$.
> - Consequently, a causal rational system is stable if and only if **all its poles lie strictly inside the unit circle**.

*Proof of the stability statement.* On $\lvert z\rvert=1$, $\sum_n\lvert h[n]z^{-n}\rvert=\sum_n\lvert h[n]\rvert$. So the unit circle is in the ROC precisely when $h$ is absolutely summable, which (§17.7) is precisely BIBO stability. $\square$

The degree condition in the causality statement is what excludes $z=\infty$ from being a pole: $H(z)=z^2/(z-0.5)$ has ROC $\lvert z\rvert>0.5$ as a right-sided sequence, but it equals $z+0.5+\frac{0.25}{z-0.5}$, whose polynomial part $z$ is an advance $\delta[n+1]$. The sequence starts at $n=-1$: right-sided, but not causal.

> [!tip] Engineering Note — The unit circle is the frequency axis
> Every real-time digital filter is causal, so its ROC is the exterior of its outermost pole. Stability then requires that exterior to contain the unit circle, i.e. every pole inside $\lvert z\rvert<1$. This single picture — poles inside the circle — is the discrete-time version of "poles in the left half-plane". The distance of the outermost pole from the circle is the stability margin: a pole at radius $0.999$ in a fixed-point implementation may cross the circle after coefficient quantisation.

### Worked examples

> [!example] Worked Example 17.10 — ROC of a sum of right-sided sequences
> Find $X(z)$ and its ROC for $x[n]=(0.5)^nu[n]+\left(-\tfrac13\right)^nu[n]$.
>
> **Step 1.** Each term is right-sided:
> $$
> (0.5)^nu[n]\leftrightarrow\frac{z}{z-0.5},\ \lvert z\rvert>0.5;\qquad\left(-\tfrac13\right)^nu[n]\leftrightarrow\frac{z}{z+\frac13},\ \lvert z\rvert>\tfrac13.
> $$
> **Step 2.** Add, over a common denominator:
> $$
> X(z)=\frac{z(z+\frac13)+z(z-0.5)}{(z-0.5)(z+\frac13)}=\frac{2z\left(z-\frac{1}{12}\right)}{(z-0.5)(z+\frac13)}.
> $$
> **Step 3.** ROC is the intersection, $\lvert z\rvert>0.5$: outside the outermost pole, as it must be for a right-sided sequence. It contains the unit circle and $z=\infty$: the sequence is causal and absolutely summable.
>
> *Check:* $X(\infty)=2=x[0]=1+1$ ✓ (initial value theorem).

> [!example] Worked Example 17.11 — The two-sided exponential $a^{\lvert n\rvert}$
> Find the Z-transform of $x[n]=a^{\lvert n\rvert}$ for real $a$, and state when it exists.
>
> **Step 1. Split.** $a^{\lvert n\rvert}=a^nu[n]+a^{-n}u[-n-1]$.
>
> **Step 2. Causal part.** $a^nu[n]\leftrightarrow\dfrac{z}{z-a}$, ROC $\lvert z\rvert>\lvert a\rvert$.
>
> **Step 3. Anticausal part.** $a^{-n}u[-n-1]=(1/a)^nu[-n-1]$. By Example 17.5 with base $1/a$, $-(1/a)^nu[-n-1]\leftrightarrow\dfrac{z}{z-1/a}$ for $\lvert z\rvert<1/\lvert a\rvert$; hence $(1/a)^nu[-n-1]\leftrightarrow-\dfrac{z}{z-1/a}$.
>
> **Step 4. Combine.**
> $$
> X(z)=\frac{z}{z-a}-\frac{z}{z-1/a}=\frac{(a-1/a)\,z}{(z-a)(z-1/a)},\qquad\lvert a\rvert<\lvert z\rvert<\frac{1}{\lvert a\rvert}.
> $$
> **Existence.** The ring is non-empty only if $\lvert a\rvert<1/\lvert a\rvert$, i.e. $\lvert a\rvert<1$. For $\lvert a\rvert\ge1$ the sequence grows (or fails to decay) in both directions and has **no** Z-transform. When it exists the ring always contains the unit circle, as it must for an absolutely summable sequence.
>
> *Check* ($a=0.5$, $z=1$): $\sum_n0.5^{\lvert n\rvert}=1+2\sum_{n\ge1}0.5^n=3$; the formula gives $(-1.5)(1)/((0.5)(-1))=3$ ✓.

> [!example] Worked Example 17.12 — One formula, three sequences
> For $X(z)=\dfrac{z}{(z-0.5)(z-2)}$, find the sequence for each possible ROC and say which is causal and which is stable.
>
> **Step 1. Partial fractions of $X(z)/z$** (method of §17.4):
> $$
> \frac{X(z)}{z}=\frac{1}{(z-0.5)(z-2)}=\frac{-2/3}{z-0.5}+\frac{2/3}{z-2}\quad\Longrightarrow\quad X(z)=-\frac23\,\frac{z}{z-0.5}+\frac23\,\frac{z}{z-2}.
> $$
> **Step 2. Assign each pole's term according to the ROC.** A pole *inside* the inner boundary of the ROC (the ROC lies outside it) gives a right-sided term $p^nu[n]$; a pole *outside* the outer boundary gives a left-sided term $-p^nu[-n-1]$.
>
> | ROC | pole 0.5 | pole 2 | $x[n]$ | causal | stable |
> |---|---|---|---|---|---|
> | $\lvert z\rvert>2$ | right | right | $\tfrac23\left(2^n-0.5^n\right)u[n]$ | yes | no |
> | $0.5<\lvert z\rvert<2$ | right | left | $-\tfrac23\,0.5^nu[n]-\tfrac23\,2^nu[-n-1]$ | no | yes |
> | $\lvert z\rvert<0.5$ | left | left | $\tfrac23\left(0.5^n-2^n\right)u[-n-1]$ | no | no |
>
> **Step 3. Interpret.** No ROC makes this system both causal and stable, because one pole is outside the unit circle. The stable choice is two-sided and equals $-\tfrac23\,0.5^{\lvert n\rvert}$ — compare Example 17.11 with $a=0.5$, where $X$ was $-1.5$ times this one.
>
> *Check* (middle ROC): evaluating the inverse contour integral (§17.4) numerically on $\lvert z\rvert=1$ gives $x[0]=-0.6667$, $x[1]=x[-1]=-0.3333$, $x[3]=x[-3]=-0.0833$, matching $-\tfrac23(0.5)^{\lvert n\rvert}$.

> [!example] Worked Example 17.13 — An empty ROC
> Does $x[n]=2^nu[n]-(0.5)^nu[-n-1]$ have a Z-transform?
>
> The causal part needs $\lvert z\rvert>2$; the anticausal part needs $\lvert z\rvert<0.5$. The intersection is empty, so **no Z-transform exists** — even though the formal sum of the two rational expressions, $\frac{z}{z-2}+\frac{z}{z-0.5}$, can be written down. Likewise $x[n]=1$ for all $n$ (or $\cos\omega_0 n$ for all $n$) has no bilateral Z-transform: the causal part needs $\lvert z\rvert>1$, the anticausal part $\lvert z\rvert<1$. Such sequences are handled by the DTFT with impulses (as in [§16.4](16-fourier-transforms.md#164-frequency-domain-representation)) or, for $n\ge0$ only, by the unilateral transform.

> [!warning] Common Pitfall — Choosing the ROC from the formula alone
> Students often pick "$\lvert z\rvert>$ largest pole" automatically. That is right only if the sequence is known to be right-sided (e.g. a causal system). If instead you are told the system is *stable*, the ROC must contain the unit circle — which, with poles on both sides of it, produces a two-sided impulse response. Read the problem's physical constraints first, then choose the ROC.

## 17.4 Inverse Z-transform

### The inversion integral

Recovering $x[n]$ from $X(z)$ is the step that turns an algebraic answer into a sequence you can plot or implement. The general formula is a contour integral, and every practical method is a shortcut for evaluating it.

> [!theorem] Inverse Z-transform
> If $X(z)$ has ROC $R_x$ and $C$ is any counter-clockwise circle $\lvert z\rvert=r$ lying in $R_x$, then
> $$
> x[n]=\frac{1}{2\pi j}\oint_C X(z)\,z^{\,n-1}\,dz,\qquad n\in\mathbb{Z}.
> $$

*Derivation.* Substitute the series for $X$ and integrate term by term (legitimate on a circle inside the ROC, where convergence is uniform):
$$
\frac{1}{2\pi j}\oint_C\sum_k x[k]z^{-k}z^{n-1}\,dz=\sum_k x[k]\,\frac{1}{2\pi j}\oint_C z^{\,n-k-1}\,dz.
$$
With $z=re^{j\theta}$, $dz=jz\,d\theta$, the inner integral is $\frac{r^{n-k}}{2\pi}\int_0^{2\pi}e^{j(n-k)\theta}d\theta=\delta[n-k]$. Only the $k=n$ term survives, leaving $x[n]$. $\square$

On the unit circle ($r=1$, $z=e^{j\omega}$, $dz=jz\,d\omega$) the formula becomes the inverse DTFT,

$$
x[n]=\frac{1}{2\pi}\int_{-\pi}^{\pi}X(e^{j\omega})\,e^{j\omega n}\,d\omega,
$$

the discrete-time counterpart of the inverse Fourier integral of [§16.2](16-fourier-transforms.md#162-inverse-transform). The Z-transform and its inverse are simply the Laurent-series expansion of $X$ in the annulus $R_x$ and the formula for its coefficients — which is why the choice of annulus (the ROC) changes the answer.

There are four practical methods:

1. **Inspection** using the table of pairs and properties;
2. **Partial-fraction expansion**, for rational $X(z)$ — the workhorse;
3. **Power-series expansion**, including **long division** — gives numerical values sample by sample and handles non-rational $X$;
4. **Residues** — evaluates the inversion integral directly; elegant for closed-form answers, especially with repeated poles.

### Method 1: partial fractions

For a rational $X(z)=N(z)/D(z)$, the aim is to write $X$ as a sum of terms $\frac{z}{z-p}$ (and $\frac{z}{(z-p)^m}$ for repeated poles), each of which inverts by the table. Because the table entries carry a factor $z$ in the numerator, the cleanest procedure expands $X(z)/z$, not $X(z)$.

> [!theorem] Partial-fraction procedure (expansion of $X(z)/z$)
> 1. Write $X(z)$ in positive powers of $z$ and factor the denominator: poles $p_1,\dots,p_K$.
> 2. Form $F(z)=X(z)/z$. If $F$ is improper (numerator degree $\ge$ denominator degree), first divide out the polynomial part.
> 3. Expand $F(z)=\sum_i\dfrac{A_i}{z-p_i}$ (+ repeated-pole terms), with **cover-up** residues for simple poles
> $$
> A_i=\left.(z-p_i)F(z)\right|_{z=p_i},
> $$
> and for a pole of order $m$ at $p$, the coefficients of $\frac{1}{(z-p)^{k}}$,
> $$
> C_{k}=\frac{1}{(m-k)!}\left.\frac{d^{\,m-k}}{dz^{\,m-k}}\Big[(z-p)^mF(z)\Big]\right|_{z=p},\quad k=1,\dots,m.
> $$
> 4. Multiply back by $z$: $X(z)=\sum_iA_i\dfrac{z}{z-p_i}+\cdots$.
> 5. Invert each term using the ROC: for $\lvert z\rvert>\lvert p\rvert$, $\frac{z}{z-p}\to p^nu[n]$; for $\lvert z\rvert<\lvert p\rvert$, $\frac{z}{z-p}\to-p^nu[-n-1]$. Repeated poles use $\dfrac{z}{(z-p)^2}\to n\,p^{\,n-1}u[n]$ and $\dfrac{z}{(z-p)^3}\to\frac{n(n-1)}{2}p^{\,n-2}u[n]$ (right-sided case).

Dividing by $z$ first has two advantages: a pole at $z=0$, if created, simply produces a $\delta[n]$ term; and no separate polynomial division is needed when the numerator and denominator of $X$ have equal degree — the common case for causal systems.

An equivalent route, preferred in DSP texts, expands $X$ in powers of $w=z^{-1}$: $X=\sum_iA_i/(1-p_iz^{-1})$ with $A_i=\left.(1-p_iz^{-1})X\right|_{z=p_i}$. Both give the same coefficients when $X$ is proper in $z^{-1}$.

> [!example] Worked Example 17.14 — Distinct real poles
> Find the causal sequence with $X(z)=\dfrac{z}{(z-0.5)(z-0.25)}$.
>
> **Step 1.** $F(z)=\dfrac{X}{z}=\dfrac{1}{(z-0.5)(z-0.25)}$.
>
> **Step 2. Cover-up.** $A_1=\dfrac{1}{0.5-0.25}=4$, $A_2=\dfrac{1}{0.25-0.5}=-4$.
>
> **Step 3.** $X(z)=4\dfrac{z}{z-0.5}-4\dfrac{z}{z-0.25}$; causal means $\lvert z\rvert>0.5$, so both terms are right-sided:
> $$
> x[n]=4\left[(0.5)^n-(0.25)^n\right]u[n].
> $$
> *Check:* $x[0]=0$ (and $X(\infty)=0$ ✓); $x[1]=4(0.25)=1$, $x[2]=4(0.25-0.0625)=0.75$. Running the recursion $x[n]=0.75x[n-1]-0.125x[n-2]+\delta[n-1]$ gives $0,\ 1,\ 0.75,\ 0.4375,\ 0.2344,\dots$, identical to the formula.

> [!example] Worked Example 17.15 — A repeated pole
> Find the causal inverse of $X(z)=\dfrac{z^2}{(z-1)(z-0.5)^2}$.
>
> **Step 1.** $F(z)=\dfrac{z}{(z-1)(z-0.5)^2}=\dfrac{A}{z-1}+\dfrac{B}{z-0.5}+\dfrac{C}{(z-0.5)^2}$.
>
> **Step 2. Simple pole.** $A=\left.\dfrac{z}{(z-0.5)^2}\right|_{z=1}=\dfrac{1}{0.25}=4$.
>
> **Step 3. Double pole.** With $G(z)=(z-0.5)^2F(z)=\dfrac{z}{z-1}$:
> $$
> C=G(0.5)=\frac{0.5}{-0.5}=-1,\qquad B=G'(0.5)=\left.\frac{-1}{(z-1)^2}\right|_{0.5}=-4.
> $$
> **Step 4.** $X(z)=4\dfrac{z}{z-1}-4\dfrac{z}{z-0.5}-\dfrac{z}{(z-0.5)^2}$, and with $\dfrac{z}{(z-0.5)^2}\leftrightarrow n(0.5)^{n-1}u[n]=2n(0.5)^nu[n]$:
> $$
> x[n]=\left[4-4(0.5)^n-2n(0.5)^n\right]u[n].
> $$
> *Check:* $x[0..5]=0,\ 1,\ 2,\ 2.75,\ 3.25,\ 3.5625$ from both the formula and a direct recursion; $x[n]\to4$, which the final value theorem confirms: $\lim_{z\to1}(z-1)X=1/0.25=4$.

> [!example] Worked Example 17.16 — Complex-conjugate poles
> Find the causal inverse of $X(z)=\dfrac{z}{z^2-z+0.5}$.
>
> **Step 1. Poles.** $z=\tfrac12\pm\tfrac12j=\tfrac{1}{\sqrt2}e^{\pm j\pi/4}$: radius $r=1/\sqrt2\approx0.707$, angle $\omega_0=\pi/4$.
>
> **Step 2. Match the damped-sine pair** $\dfrac{zr\sin\omega_0}{z^2-2rz\cos\omega_0+r^2}$. Here $2r\cos\omega_0=2\cdot\tfrac{1}{\sqrt2}\cdot\tfrac{1}{\sqrt2}=1$ ✓ and $r^2=0.5$ ✓, while $r\sin\omega_0=0.5$. So $X(z)=2\times\dfrac{0.5\,z}{z^2-z+0.5}$ and
> $$
> x[n]=2\left(\tfrac{1}{\sqrt2}\right)^n\sin\left(\tfrac{\pi n}{4}\right)u[n].
> $$
> **Alternative (complex partial fractions).** $F=X/z=\dfrac{1}{(z-p)(z-p^*)}$, $A=\dfrac{1}{p-p^*}=\dfrac{1}{j}=-j$, and the conjugate pole has $A^*=j$. Then $x[n]=-jp^n+jp^{*n}=2\operatorname{Im}\{p^n\}=2r^n\sin(\omega_0n)$, the same answer. In general a conjugate pair with residues $A, A^*$ contributes $2\lvert A\rvert r^n\cos(\omega_0n+\angle A)\,u[n]$.
>
> *Check:* $x[0..7]=0,\ 1,\ 1,\ 0.5,\ 0,\ -0.25,\ -0.25,\ -0.125$, agreeing with the recursion $x[n]=x[n-1]-0.5x[n-2]+\delta[n-1]$.

> [!example] Worked Example 17.17 — Equal-degree numerator and a pole at the origin
> Find the causal inverse of $X(z)=\dfrac{z^2+1}{z^2-0.25}$.
>
> **Step 1.** $F(z)=\dfrac{z^2+1}{z(z-0.5)(z+0.5)}$: dividing by $z$ creates a pole at $0$, and $F$ is now proper.
>
> **Step 2. Cover-up.** $A_0=\dfrac{1}{-0.25}=-4$, $A_{0.5}=\dfrac{1.25}{0.5\cdot1}=2.5$, $A_{-0.5}=\dfrac{1.25}{(-0.5)(-1)}=2.5$.
>
> **Step 3.** $X(z)=-4+2.5\dfrac{z}{z-0.5}+2.5\dfrac{z}{z+0.5}$, so
> $$
> x[n]=-4\,\delta[n]+2.5\left[(0.5)^n+(-0.5)^n\right]u[n].
> $$
> The two exponentials cancel at odd $n$ and add at even $n$: $x=1,\ 0,\ 1.25,\ 0,\ 0.3125,\ 0,\dots$.
>
> *Check:* in powers of $z^{-1}$, $X=\dfrac{1+z^{-2}}{1-0.25z^{-2}}=(1+z^{-2})(1+0.25z^{-2}+0.0625z^{-4}+\cdots)=1+1.25z^{-2}+0.3125z^{-4}+\cdots$ ✓.

### Method 2: power series and long division

Since $X(z)=\sum x[n]z^{-n}$, any method that expands $X$ as a power series in the ROC reads off $x[n]$ as the coefficients. For rational $X$ the expansion is obtained by **polynomial long division**:

- For a **right-sided** (causal) ROC $\lvert z\rvert>r_1$, expand in *descending* powers of $z$, i.e. ascending powers of $z^{-1}$: arrange numerator and denominator in decreasing powers of $z$ and divide.
- For a **left-sided** ROC $\lvert z\rvert<r_2$, expand in *ascending* powers of $z$: arrange both polynomials in increasing powers of $z$ and divide.

Long division does not normally produce a closed form, but it is foolproof for the first few samples, it is how a computer "inverts" (it is the recursion of the difference equation), and it is the method of choice for checking a partial-fraction answer.

> [!example] Worked Example 17.18 — Causal long division
> Find the first four samples of the causal sequence with $X(z)=\dfrac{z}{z^2-0.75z+0.125}$ (Example 17.14).
>
> Divide $z$ by $z^2-0.75z+0.125$ in descending powers:
>
> | Step | Quotient term | Remainder after subtraction |
> |---|---|---|
> | $z\div z^2$ | $z^{-1}$ | $z-(z-0.75+0.125z^{-1})=0.75-0.125z^{-1}$ |
> | $0.75\div z^2$ | $0.75z^{-2}$ | $-0.125z^{-1}-(-0.5625z^{-1}+0.09375z^{-2})=0.4375z^{-1}-0.09375z^{-2}$ |
> | $0.4375z^{-1}\div z^2$ | $0.4375z^{-3}$ | $\dots$ |
>
> So $X=0+z^{-1}+0.75z^{-2}+0.4375z^{-3}+\cdots$, giving $x[0..3]=0,\ 1,\ 0.75,\ 0.4375$ — exactly the partial-fraction result.

> [!example] Worked Example 17.19 — Anticausal long division
> Expand $X(z)=\dfrac{z}{z-a}$ with ROC $\lvert z\rvert<\lvert a\rvert$.
>
> Arrange in ascending powers of $z$: $\dfrac{z}{-a+z}$. Dividing, the first quotient term is $z/(-a)=-a^{-1}z$; the remainder is $z-(-a^{-1}z)(-a+z)=a^{-1}z^2$; the next quotient term is $-a^{-2}z^2$, and so on:
> $$
> X(z)=-a^{-1}z-a^{-2}z^2-a^{-3}z^3-\cdots=\sum_{n=-\infty}^{-1}\left(-a^{n}\right)z^{-n}.
> $$
> Hence $x[n]=-a^nu[-n-1]$, in agreement with Example 17.5.

Non-rational transforms are inverted by known series. The logarithm is the classic case.

> [!example] Worked Example 17.20 — A logarithmic transform
> Find the right-sided $x[n]$ with $X(z)=\ln\left(1+az^{-1}\right)$, $\lvert z\rvert>\lvert a\rvert$.
>
> Use $\ln(1+w)=\sum_{n=1}^{\infty}\dfrac{(-1)^{n+1}}{n}w^n$ for $\lvert w\rvert<1$ with $w=az^{-1}$ (valid exactly when $\lvert z\rvert>\lvert a\rvert$):
> $$
> X(z)=\sum_{n=1}^{\infty}\frac{(-1)^{n+1}a^n}{n}z^{-n}\quad\Longrightarrow\quad x[n]=\frac{(-1)^{n+1}a^n}{n}\,u[n-1].
> $$
> *Cross-check via the differentiation property:* $-z\,X'(z)=\dfrac{az^{-1}}{1+az^{-1}}=\dfrac{a}{z+a}$, whose inverse is $a(-a)^{n-1}u[n-1]$; dividing by $n$ reproduces $x[n]$. Sequences of this type appear in cepstral analysis of speech and echo signals.

### Method 3: residues

When $X(z)$ is rational, the integrand $X(z)z^{n-1}$ of the inversion integral is analytic everywhere except at finitely many poles, and the residue theorem ([Ch. 18](18-complex-analysis.md)) evaluates the integral as a sum.

> [!theorem] Residue formula for the inverse Z-transform
> $$
> x[n]=\sum_{\text{poles } p_k\text{ inside } C}\operatorname{Res}_{z=p_k}\left[X(z)\,z^{n-1}\right],
> $$
> where $C$ is a circle in the ROC. For a pole of order $m$ at $p$,
> $$
> \operatorname{Res}_{z=p}F(z)=\frac{1}{(m-1)!}\left.\frac{d^{\,m-1}}{dz^{\,m-1}}\Big[(z-p)^mF(z)\Big]\right|_{z=p};\qquad m=1:\ \ \operatorname{Res}=\left.(z-p)F(z)\right|_{z=p}.
> $$

Two practical remarks:

- For a **causal** $X$ and $n\ge0$, the poles inside $C$ are just the poles of $X$; but for $n<0$ the factor $z^{n-1}$ adds a pole of order $1-n$ at the origin. Rather than computing those awkward residues, use the fact that the sum of all residues, including the one at infinity, is zero: $x[n]=-\sum(\text{residues outside }C)$. For a causal $X$ there are no poles outside $C$ and, for $n<0$, no residue at infinity, so $x[n]=0$ — confirming causality.
- Even for $n\ge0$, a factor $z^{n-1}$ with $n=0$ creates a pole at $z=0$ if $X(0)\ne0$. That is the $n=0$ case requiring care; the $X(z)/z$ partial-fraction method handles it automatically.

> [!example] Worked Example 17.21 — Residues for simple poles
> Invert $X(z)=\dfrac{z}{(z-0.5)(z-0.25)}$, $\lvert z\rvert>0.5$, by residues.
>
> For $n\ge0$, $X(z)z^{n-1}=\dfrac{z^n}{(z-0.5)(z-0.25)}$ has poles only at $0.5$ and $0.25$ (no pole at the origin, because the factor $z$ in $X$ cancels the $z^{-1}$):
> $$
> \operatorname{Res}_{0.5}=\frac{0.5^n}{0.25}=4(0.5)^n,\qquad\operatorname{Res}_{0.25}=\frac{0.25^n}{-0.25}=-4(0.25)^n.
> $$
> Sum: $x[n]=4(0.5^n-0.25^n)$ for $n\ge0$. For $n<0$ there are no poles outside $C$ (taken just larger than $0.5$), so $x[n]=0$. Same as Example 17.14.

> [!example] Worked Example 17.22 — Residues for a repeated pole, and the $n=0$ subtlety
> (a) Invert $X(z)=\dfrac{z}{(z-a)^2}$, $\lvert z\rvert>\lvert a\rvert$.
>
> For $n\ge0$, $X z^{n-1}=\dfrac{z^n}{(z-a)^2}$ has a double pole at $a$:
> $$
> x[n]=\left.\frac{d}{dz}z^n\right|_{z=a}=n\,a^{n-1},\qquad n\ge0,
> $$
> confirming the table entry $z/(z-a)^2\leftrightarrow na^{n-1}u[n]$ with no differentiation of transforms needed.
>
> (b) Invert $X(z)=\dfrac{1}{z-0.5}$, $\lvert z\rvert>0.5$. Here $Xz^{n-1}=\dfrac{z^{n-1}}{z-0.5}$.
> - $n\ge1$: only the pole at $0.5$; $x[n]=0.5^{n-1}$.
> - $n=0$: $\dfrac{1}{z(z-0.5)}$ has poles at $0$ and $0.5$, with residues $\dfrac{1}{-0.5}=-2$ and $\dfrac{1}{0.5}=2$. Sum: $x[0]=0$.
>
> So $x[n]=0.5^{n-1}u[n-1]$ — the delayed geometric sequence, as the shift property predicts ($z^{-1}\cdot\frac{z}{z-0.5}$). Forgetting the pole at the origin would give the wrong value $x[0]=2$.

### Comparing the methods

| Method | Gives | Best for | Weakness |
|---|---|---|---|
| Inspection / table | closed form | standard shapes, combined with properties | needs recognition |
| Partial fractions | closed form | any rational $X$; two-sided ROCs | repeated and complex poles need care |
| Long division | numbers, sample by sample | checking, first few samples, computer implementation | no closed form; errors accumulate by hand |
| Power series | closed form or numbers | non-rational $X$ (logs, exponentials) | requires a known series |
| Residues | closed form | repeated poles; theory | pole at origin for small $n$ |

> [!warning] Common Pitfall — Expanding $X(z)$ instead of $X(z)/z$
> Expanding $X(z)=\frac{z}{(z-0.5)(z-0.25)}$ directly gives $\frac{2}{z-0.5}-\frac{1}{z-0.25}$, whose terms are $\frac{1}{z-p}=z^{-1}\frac{z}{z-p}\leftrightarrow p^{n-1}u[n-1]$ — correct, but only if you remember the delay. Inverting $\frac{2}{z-0.5}$ as "$2(0.5)^n$" is a classic one-sample error. Dividing by $z$ first avoids it.

> [!tip] Engineering Note — Numerical inversion
> In software, a rational $X(z)$ in powers of $z^{-1}$ is inverted by running the difference equation it defines with an impulse input (MATLAB `impz`, SciPy `signal.dimpulse` or `lfilter(b, a, impulse)`), and partial fractions are computed by `residuez`. Root-finding for high-order denominators is ill-conditioned; for orders above about ten, factor the filter into second-order sections rather than computing all poles of the full polynomial.

## 17.5 Difference equations

### Where difference equations come from

A **linear constant-coefficient difference equation** (LCCDE) relates the present output of a discrete-time system to past outputs and present and past inputs:

$$
\sum_{k=0}^{N}a_k\,y[n-k]=\sum_{k=0}^{M}b_k\,x[n-k],\qquad a_0\ne0.
$$

The **order** is $N$, the number of past outputs the system remembers. Difference equations arise in three ways that an electrical engineer meets constantly:

- **By design**: a digital filter *is* its difference equation; the coefficients $a_k, b_k$ are what the designer chooses and the processor stores.
- **By discretisation**: replacing derivatives in a differential equation by finite differences (Ch. 22) — for instance a forward-Euler model of an RC circuit — produces an LCCDE.
- **By counting**: recurrences in sampled-data control, compound-interest-like growth, ladder networks of identical sections and combinatorial counts (Fibonacci) are naturally indexed by an integer.

Normalising so that $a_0=1$, the equation can be solved for $y[n]$:

$$
y[n]=-\sum_{k=1}^{N}a_k\,y[n-k]+\sum_{k=0}^{M}b_k\,x[n-k].
$$

This is a **recursion**: given $N$ initial values $y[-1],\dots,y[-N]$ and the input, it computes $y[0], y[1],\dots$ one after another. Unlike a differential equation, a difference equation can therefore always be "solved" numerically with nothing more than multiplication and addition. What the Z-transform supplies is the **closed form** — the formula valid for every $n$ at once — and with it insight into the modes, stability and steady state.

### The unilateral shift theorem

The bilateral shift rule $x[n-1]\leftrightarrow z^{-1}X(z)$ assumes knowledge of the whole sequence. In an initial-value problem we know $y[n]$ for $n\ge0$ only through the equation, plus a few stored values $y[-1], y[-2],\dots$. The unilateral transform keeps exactly this information.

> [!theorem] Unilateral time-shift (delay and advance)
> With $\mathcal{Y}(z)=\sum_{n\ge0}y[n]z^{-n}$:
> $$
> \begin{aligned}
> y[n-1]&\;\longleftrightarrow\;z^{-1}\mathcal{Y}(z)+y[-1],\\
> y[n-2]&\;\longleftrightarrow\;z^{-2}\mathcal{Y}(z)+z^{-1}y[-1]+y[-2],\\
> y[n-k]&\;\longleftrightarrow\;z^{-k}\mathcal{Y}(z)+\sum_{m=1}^{k}y[-m]\,z^{-(k-m)},\\
> y[n+1]&\;\longleftrightarrow\;z\,\mathcal{Y}(z)-z\,y[0],\\
> y[n+2]&\;\longleftrightarrow\;z^2\mathcal{Y}(z)-z^2y[0]-z\,y[1].
> \end{aligned}
> $$

*Derivation for one delay.* $\sum_{n\ge0}y[n-1]z^{-n}=y[-1]+\sum_{n\ge1}y[n-1]z^{-n}=y[-1]+z^{-1}\sum_{m\ge0}y[m]z^{-m}$. *For one advance:* $\sum_{n\ge0}y[n+1]z^{-n}=z\sum_{m\ge1}y[m]z^{-m}=z(\mathcal{Y}-y[0])$. The others follow by repeating. $\square$

The parallel with the Laplace derivative rule $\mathcal{L}\{\dot y\}=sY-y(0^-)$ is exact: delays bring in *past* values ($y[-1]$, $y[-2]$), advances bring in *initial* values ($y[0]$, $y[1]$). Use the delay form when the equation is written with $y[n-k]$ and the data are $y[-1],y[-2],\dots$; use the advance form when it is written with $y[n+k]$ and the data are $y[0],y[1],\dots$.

For a **causal input** applied at $n=0$, $x[-1]=x[-2]=\dots=0$, so the input terms shift without extra terms.

### The solution procedure

> [!theorem] Solving an LCCDE with the unilateral Z-transform
> 1. Take the unilateral transform of every term, using the shift theorem to bring in the initial conditions.
> 2. Solve the resulting **algebraic** equation for $\mathcal{Y}(z)$.
> 3. Invert by partial fractions of $\mathcal{Y}(z)/z$ (ROC: exterior of the outermost pole, since $y$ is considered for $n\ge0$).
> 4. Check $y[0]$ (and $y[1]$) against a direct run of the recursion.

Step 2 always produces the structure

$$
\mathcal{Y}(z)=\underbrace{\frac{B(z)}{A(z)}\,\mathcal{X}(z)}_{\text{zero-state response}}+\underbrace{\frac{C(z)}{A(z)}}_{\text{zero-input response}},
\qquad A(z)=\sum_{k=0}^{N}a_kz^{-k},\quad B(z)=\sum_{k=0}^{M}b_kz^{-k},
$$

where $C(z)$ is a polynomial built from the initial conditions. The **zero-state response** (ZSR) is what the input produces from rest; the **zero-input response** (ZIR) is what the stored initial conditions produce with no input. Both share the denominator $A(z)$, whose roots — the **characteristic roots** or **natural modes** — are the same as those of the classical characteristic equation $\lambda^N+a_1\lambda^{N-1}+\dots+a_N=0$. The ratio $B/A$ is the system function $H(z)$ of §17.6.

A second, equally useful split is into **natural** (transient) and **forced** (steady-state) parts: the natural response collects the terms from the poles of $H$, the forced response those from the poles of $\mathcal{X}$.

### Worked examples

> [!example] Worked Example 17.23 — First order with a step input
> Solve $y[n]-0.5\,y[n-1]=u[n]$ with $y[-1]=1$.
>
> **Step 1. Transform.** $\mathcal{Y}-0.5\left(z^{-1}\mathcal{Y}+y[-1]\right)=\dfrac{z}{z-1}$.
>
> **Step 2. Solve.** $\mathcal{Y}\left(1-0.5z^{-1}\right)=0.5+\dfrac{z}{z-1}$, so
> $$
> \mathcal{Y}(z)=\underbrace{\frac{0.5z}{z-0.5}}_{\text{ZIR}}+\underbrace{\frac{z^2}{(z-1)(z-0.5)}}_{\text{ZSR}}.
> $$
> **Step 3. Invert.** The ZIR is $0.5(0.5)^n=(0.5)^{n+1}$. For the ZSR, $\dfrac{z}{(z-1)(z-0.5)}=\dfrac{2}{z-1}-\dfrac{1}{z-0.5}$, so ZSR $=2-(0.5)^n$. Adding,
> $$
> y[n]=2-(0.5)^n+(0.5)^{n+1}=2-(0.5)^{n+1},\qquad n\ge0.
> $$
> **Step 4. Check.** Recursion: $y[0]=0.5(1)+1=1.5$, $y[1]=0.5(1.5)+1=1.75$, $y[2]=1.875$; formula: $2-0.5=1.5$, $2-0.25=1.75$, $2-0.125=1.875$ ✓. The forced response is the constant $2=H(1)$, the DC gain; the natural response $-(0.5)^{n+1}$ decays.

> [!example] Worked Example 17.24 — Second order: ZIR, ZSR, natural and forced parts
> Solve $y[n]-0.75\,y[n-1]+0.125\,y[n-2]=u[n]$ with $y[-1]=1$, $y[-2]=2$.
>
> **Step 1. Transform.**
> $$
> \mathcal{Y}-0.75\left(z^{-1}\mathcal{Y}+1\right)+0.125\left(z^{-2}\mathcal{Y}+z^{-1}\cdot1+2\right)=\frac{z}{z-1}.
> $$
> **Step 2. Collect.** The initial-condition terms give $0.75-0.125z^{-1}-0.25=0.5-0.125z^{-1}$, so
> $$
> \mathcal{Y}\,(1-0.75z^{-1}+0.125z^{-2})=\frac{z}{z-1}+0.5-0.125z^{-1}.
> $$
> With $1-0.75z^{-1}+0.125z^{-2}=(1-0.5z^{-1})(1-0.25z^{-1})$:
> $$
> \mathcal{Y}_{\text{zi}}=\frac{0.5(1-0.25z^{-1})}{(1-0.5z^{-1})(1-0.25z^{-1})}=\frac{0.5z}{z-0.5},\qquad
> \mathcal{Y}_{\text{zs}}=\frac{z^3}{(z-1)(z-0.5)(z-0.25)}.
> $$
> The initial conditions happen to cancel the $0.25$ mode in the ZIR.
>
> **Step 3. Invert.** ZIR: $0.5(0.5)^n$. ZSR: expanding $\mathcal{Y}_{\text{zs}}/z$ by cover-up,
> $$
> A_1=\frac{1}{0.5\cdot0.75}=\frac83,\quad A_{0.5}=\frac{0.25}{(-0.5)(0.25)}=-2,\quad A_{0.25}=\frac{0.0625}{(-0.75)(-0.25)}=\frac13,
> $$
> so ZSR $=\tfrac83-2(0.5)^n+\tfrac13(0.25)^n$. Total:
> $$
> y[n]=\underbrace{\tfrac83}_{\text{forced}}\;\underbrace{-\,1.5(0.5)^n+\tfrac13(0.25)^n}_{\text{natural}},\qquad n\ge0.
> $$
> **Step 4. Check.** Formula: $y[0]=2.6667-1.5+0.3333=1.5$; $y[1]=2.6667-0.75+0.0833=2.0$; $y[2]=2.6667-0.375+0.0208=2.3125$. Recursion: $y[0]=0.75(1)-0.125(2)+1=1.5$, $y[1]=0.75(1.5)-0.125(1)+1=2.0$, $y[2]=0.75(2)-0.125(1.5)+1=2.3125$ ✓. The steady state $8/3=1/(1-0.75+0.125)=H(1)$.

> [!example] Worked Example 17.25 — The Fibonacci recurrence (advance form)
> Solve $y[n+2]=y[n+1]+y[n]$ with $y[0]=0$, $y[1]=1$.
>
> **Step 1. Advance rule.** $\left(z^2\mathcal{Y}-z^2y[0]-zy[1]\right)=\left(z\mathcal{Y}-zy[0]\right)+\mathcal{Y}$, i.e. $z^2\mathcal{Y}-z=z\mathcal{Y}+\mathcal{Y}$.
>
> **Step 2.** $\mathcal{Y}(z)=\dfrac{z}{z^2-z-1}=\dfrac{z}{(z-\varphi)(z-\psi)}$, with $\varphi=\dfrac{1+\sqrt5}{2}\approx1.618$ and $\psi=\dfrac{1-\sqrt5}{2}\approx-0.618$.
>
> **Step 3.** $\dfrac{\mathcal{Y}}{z}=\dfrac{1}{\varphi-\psi}\left(\dfrac{1}{z-\varphi}-\dfrac{1}{z-\psi}\right)$ with $\varphi-\psi=\sqrt5$, so
> $$
> y[n]=\frac{\varphi^n-\psi^n}{\sqrt5}\qquad\text{(Binet's formula).}
> $$
> *Check:* $y[10]=(122.992-0.00813)/2.23607=55.000$ ✓. One characteristic root lies outside the unit circle, so the sequence grows geometrically, by the factor $\varphi$ per step: as a "system" this recursion is unstable.

> [!example] Worked Example 17.26 — Input at a natural frequency (discrete resonance)
> Solve $y[n]-0.5\,y[n-1]=(0.5)^nu[n]$ from rest.
>
> $\mathcal{Y}(1-0.5z^{-1})=\dfrac{z}{z-0.5}$, so $\mathcal{Y}=\dfrac{z^2}{(z-0.5)^2}$, a **double** pole: the input mode coincides with the natural mode. From the table ($m=2$): $\dfrac{z^2}{(z-a)^2}\leftrightarrow(n+1)a^nu[n]$, so
> $$
> y[n]=(n+1)(0.5)^n,\qquad n\ge0.
> $$
> The factor $(n+1)$ is the discrete counterpart of the $te^{-at}$ term in a continuous system driven at its own pole. *Check:* $y[0]=1$, $y[1]=0.5(1)+0.5=1$, $y[2]=0.5(1)+0.25=0.75$; formula $1,\ 1,\ 0.75$ ✓.

> [!example] Worked Example 17.27 — Forward-Euler model of an RC low-pass circuit
> An RC low-pass filter obeys $\tau\,\dot v=x-v$, $\tau=RC$. Approximating $\dot v(nT)\approx\left(v[n+1]-v[n]\right)/T$ gives
> $$
> v[n+1]=\left(1-\frac{T}{\tau}\right)v[n]+\frac{T}{\tau}\,x[n].
> $$
> Find the step response for $\tau=1$ ms, $T=0.25$ ms, $v[0]=0$, and compare with the exact samples $1-e^{-nT/\tau}$.
>
> **Step 1.** With $\alpha=T/\tau=0.25$ and $x[n]=u[n]$, the advance rule gives $z\mathcal{V}-zv[0]=0.75\mathcal{V}+0.25\dfrac{z}{z-1}$, so
> $$
> \mathcal{V}(z)=\frac{0.25z}{(z-1)(z-0.75)}.
> $$
> **Step 2.** $\mathcal{V}/z=\dfrac{1}{z-1}-\dfrac{1}{z-0.75}$, so $v[n]=1-(0.75)^n$.
>
> **Step 3. Compare.** The exact samples are $1-(e^{-0.25})^n=1-(0.7788)^n$. At $n=4$ ($t=\tau$): Euler gives $0.6836$, exact $0.6321$. The discrete pole $0.75$ approximates the exact $e^{-T/\tau}=0.7788$ — a first-order Taylor approximation of $z=e^{sT}$ (§17.7).
>
> **Step 4. Numerical stability.** The pole $1-T/\tau$ lies inside the unit circle only if $0<T/\tau<2$. With $T>2\tau$ the simulation diverges although the circuit is perfectly stable — the origin of the step-size limits studied in [Ch. 22](../part-7-numerical-computational/22-numerical-methods.md).

### Classical solution, for comparison

The method of undetermined coefficients works for difference equations exactly as for differential equations: the homogeneous solution is $\sum c_i\lambda_i^n$ over the characteristic roots $\lambda_i$ (with $n^k\lambda^n$ for repeated roots), a particular solution is guessed in the form of the input, and the constants are fitted to the initial conditions. For Example 17.24: characteristic roots $0.5, 0.25$; particular solution $y_p=K$ with $K(1-0.75+0.125)=1$, $K=8/3$; then $y=\tfrac83+c_1(0.5)^n+c_2(0.25)^n$ fitted to $y[0]=1.5$, $y[1]=2$ gives $c_1=-1.5$, $c_2=\tfrac13$. The Z-transform does the same work but builds the initial conditions in automatically and handles resonant inputs (Example 17.26) without a separate rule.

> [!warning] Common Pitfall — Using $y[0]$ where the delay rule needs $y[-1]$
> The delay rule needs the values *before* the input is applied ($y[-1], y[-2]$); the advance rule needs $y[0], y[1]$. If a problem gives $y[0]$ for an equation written in delays, first back-solve the recursion for $y[-1]$ (with the known $x[0]$) or rewrite the equation in advance form. Mixing them silently shifts the answer by one sample.

> [!tip] Engineering Note — Initial conditions are register contents
> In an implemented IIR filter, $y[-1]$ and $y[-2]$ are literally the numbers sitting in the delay registers when processing starts. Clearing them (zero state) at power-up gives the zero-state response; a filter that is switched between coefficient sets on the fly carries its old register contents as initial conditions, and the resulting zero-input transient is the "click" heard when audio equaliser settings change abruptly.

## 17.6 System-function representation

### The system function

For an LTI system with impulse response $h[n]$, the convolution theorem gives $Y(z)=H(z)X(z)$. The transform of the impulse response therefore characterises the system completely, in the same way the transfer function $H(s)$ does in continuous time.

> [!definition] System function (transfer function)
> The **system function** of a discrete-time LTI system is
> $$
> H(z)=\mathcal{Z}\{h[n]\}=\frac{Y(z)}{X(z)}\quad(\text{zero initial state}),
> $$
> together with its ROC. Evaluated on the unit circle, $H(e^{j\omega})$ is the **frequency response**.

For a system described by the LCCDE $\sum_{k=0}^Na_ky[n-k]=\sum_{k=0}^Mb_kx[n-k]$, transforming with zero initial conditions (bilateral shift rule) gives

$$
H(z)=\frac{\displaystyle\sum_{k=0}^{M}b_kz^{-k}}{\displaystyle\sum_{k=0}^{N}a_kz^{-k}}=\frac{b_0+b_1z^{-1}+\cdots+b_Mz^{-M}}{a_0+a_1z^{-1}+\cdots+a_Nz^{-N}}.
$$

The coefficients of the difference equation are the coefficients of the numerator and denominator polynomials, read off directly. Going the other way, any rational $H(z)$ written in powers of $z^{-1}$ is a difference equation: cross-multiply $Y\cdot A(z)=X\cdot B(z)$ and interpret each $z^{-k}$ as a $k$-sample delay.

An LCCDE alone does not fix the ROC. If the system is known to be causal (the usual case), the ROC is the exterior of the outermost pole; if it is known to be stable, the ROC contains the unit circle. These two pieces of information, plus $H(z)$, determine $h[n]$ uniquely.

### Pole–zero (factored) form

Factoring numerator and denominator,

$$
H(z)=b_0\,z^{N-M}\,\frac{\prod_{k=1}^{M}\left(z-q_k\right)}{\prod_{k=1}^{N}\left(z-p_k\right)}\qquad(a_0=1),
$$

with **zeros** $q_k$, **poles** $p_k$, and $\lvert N-M\rvert$ extra zeros or poles at $z=0$. The **pole–zero plot** — poles marked $\times$, zeros $\circ$, with the unit circle drawn for reference — together with the gain $b_0$ is a complete description. It shows at a glance:

- **stability** of a causal system: all poles inside the unit circle;
- **modes** of the natural response: each pole $p=re^{j\theta}$ contributes $r^n$-weighted oscillation at $\theta$ rad/sample;
- **frequency selectivity**: poles near the unit circle at angle $\theta$ produce gain peaks near $\omega=\theta$; zeros on or near the circle produce nulls (§17.7).

> [!example] Worked Example 17.28 — From difference equation to pole–zero plot
> A causal system obeys $y[n]=0.9\,y[n-1]-0.81\,y[n-2]+x[n]-x[n-2]$. Find $H(z)$, its poles and zeros, and its ROC.
>
> **Step 1.** Rearranged, $y[n]-0.9y[n-1]+0.81y[n-2]=x[n]-x[n-2]$, so
> $$
> H(z)=\frac{1-z^{-2}}{1-0.9z^{-1}+0.81z^{-2}}=\frac{z^2-1}{z^2-0.9z+0.81}.
> $$
> **Step 2. Zeros:** $z^2=1$, i.e. $z=\pm1$ — on the unit circle at $\omega=0$ and $\omega=\pi$.
>
> **Step 3. Poles:** $z=\dfrac{0.9\pm\sqrt{0.81-3.24}}{2}=0.45\pm j0.7794=0.9e^{\pm j\pi/3}$ (check: $2r\cos\theta=2(0.9)(0.5)=0.9$ ✓, $r^2=0.81$ ✓).
>
> **Step 4. ROC.** Causal $\Rightarrow\lvert z\rvert>0.9$, which contains the unit circle: the system is stable. The zeros block DC and $\omega=\pi$; the poles near the circle at $\pm60^\circ$ boost $\omega\approx\pi/3$. This is a **second-order band-pass resonator**; its magnitude response is drawn in §17.7.

```desmos-graph
left=-1.6; right=1.6; bottom=-1.6; top=1.6
---
x^2+y^2=1|black|dashed
x^2+y^2>0.81|blue
(0.45,0.7794)|red|label:pole r=0.9 at 60°
(0.45,-0.7794)|red|label:pole r=0.9 at −60°
(1,0)|green|label:zero +1
(-1,0)|green|label:zero −1
```

*Pole–zero plot of $H(z)=(z^2-1)/(z^2-0.9z+0.81)$ with the causal ROC $\lvert z\rvert>0.9$ shaded. Poles (red) at radius 0.9 and angle $\pm\pi/3$; zeros (green) on the unit circle at $\pm1$.*

### FIR and IIR systems

> [!definition] FIR and IIR
> A system is **finite impulse response (FIR)** if $h[n]$ has finitely many non-zero samples. For a causal FIR system $H(z)=\sum_{k=0}^{M}b_kz^{-k}$: all poles are at $z=0$, so an FIR system is always stable. Its difference equation has no feedback: $y[n]=\sum_k b_kx[n-k]$, and $h[n]=b_n$.
>
> A system is **infinite impulse response (IIR)** if $h[n]$ has infinitely many non-zero samples — for rational $H$, when there is at least one pole not cancelled by a zero away from the origin. IIR systems need feedback ($a_k\ne0$ for some $k\ge1$) and may be unstable.

> [!example] Worked Example 17.29 — The moving-average filter: FIR with a recursive form
> The 4-point moving average is $y[n]=\tfrac14\left(x[n]+x[n-1]+x[n-2]+x[n-3]\right)$.
>
> **Step 1.** $H(z)=\tfrac14\left(1+z^{-1}+z^{-2}+z^{-3}\right)$. By the geometric sum,
> $$
> H(z)=\frac14\,\frac{1-z^{-4}}{1-z^{-1}}=\frac14\,\frac{z^4-1}{z^3(z-1)}.
> $$
> **Step 2. Zeros:** the fourth roots of unity $1, j, -1, -j$; **poles:** three at $z=0$ and one at $z=1$. The pole at $1$ cancels the zero at $1$, leaving zeros at $\omega=\pm\pi/2$ and $\pi$ on the unit circle — frequencies that are completely rejected.
>
> **Step 3. Recursive implementation.** The ratio form corresponds to the difference equation $y[n]=y[n-1]+\tfrac14\left(x[n]-x[n-4]\right)$: one addition and one subtraction per sample regardless of length (the **running-sum** or CIC structure), instead of $M$ additions. Its impulse response is still finite, because the pole is cancelled.
>
> *Check:* $\lvert H(e^{j0})\rvert=1$ (unit DC gain), $\lvert H(e^{j\pi/2})\rvert=0$, $\lvert H(e^{j\pi/4})\rvert=0.653$.

### From impulse response to system function, and back

> [!example] Worked Example 17.30 — Identifying a system from its impulse response
> A causal system has $h[n]=\left[2(0.5)^n-(0.25)^n\right]u[n]$. Find $H(z)$ and the difference equation.
>
> **Step 1.**
> $$
> H(z)=\frac{2z}{z-0.5}-\frac{z}{z-0.25}=\frac{z\left[2(z-0.25)-(z-0.5)\right]}{(z-0.5)(z-0.25)}=\frac{z^2}{z^2-0.75z+0.125}.
> $$
> **Step 2.** In powers of $z^{-1}$: $H=\dfrac{1}{1-0.75z^{-1}+0.125z^{-2}}$, so
> $$
> y[n]=0.75\,y[n-1]-0.125\,y[n-2]+x[n].
> $$
> *Check:* $h[0]=2-1=1$ and the recursion gives $h[0]=1$; $h[1]=1-0.25=0.75$ and the recursion gives $0.75(1)=0.75$ ✓. This is the system of Example 17.24, whose step response had steady state $H(1)=1/0.375=8/3$.

### Block diagrams and realisations

Three primitives build every LCCDE: the **adder**, the **multiplier** (gain) and the **unit delay** $z^{-1}$. Because the same $H(z)$ can be wired in many ways, the choice of **realisation** (structure) matters for memory, arithmetic and — above all in fixed-point hardware — sensitivity to coefficient rounding.

> [!definition] Direct forms
> Write $H(z)=B(z)\cdot\dfrac{1}{A(z)}$.
> - **Direct form I** implements $B(z)$ (a tapped delay line on $x$, $M$ delays) followed by $1/A(z)$ (a tapped delay line on $y$ with feedback, $N$ delays): $M+N$ delays.
> - **Direct form II** reverses the order, $\dfrac{1}{A(z)}$ then $B(z)$. Both sections then delay the same intermediate signal $w[n]$, so the delay lines merge: $w[n]=x[n]-\sum_{k=1}^Na_kw[n-k]$, $y[n]=\sum_{k=0}^Mb_kw[n-k]$ — only $\max(M,N)$ delays, the minimum possible (**canonical**).
> - **Transposed direct form II** reverses every arrow and swaps inputs and outputs (transposition preserves $H$); it is often preferred in floating point because the adders see smaller dynamic range.

Writing the DF-II equations for Example 17.28 ($b=[1,0,-1]$, $a=[1,-0.9,0.81]$):

$$
w[n]=x[n]+0.9\,w[n-1]-0.81\,w[n-2],\qquad y[n]=w[n]-w[n-2],
$$

two delays (storing $w[n-1]$, $w[n-2]$), three multiplications and three additions per sample.

> [!definition] Cascade and parallel forms
> - **Cascade (series) form:** factor $H(z)=G\prod_kH_k(z)$ into first- and second-order sections ("biquads"), each with real coefficients — complex-conjugate pole and zero pairs are kept together. This is the standard implementation of IIR filters of order above two.
> - **Parallel form:** a partial-fraction expansion in $z^{-1}$, $H(z)=C+\sum_k\dfrac{\beta_{0k}+\beta_{1k}z^{-1}}{1+\alpha_{1k}z^{-1}+\alpha_{2k}z^{-2}}$, sums the outputs of independent sections.

The reason for cascades is numerical: the roots of a high-degree polynomial are extremely sensitive to its coefficients (Ch. 22). A 10th-order narrow-band filter stored as one denominator polynomial may become unstable after rounding its coefficients to 16 bits, while the same filter as five biquads, each with its own well-conditioned quadratic, stays stable.

### Interconnections

The algebra of system functions mirrors that of Laplace transfer functions:

| Connection | System function |
|---|---|
| Cascade $H_1$ then $H_2$ | $H_1(z)H_2(z)$ (order irrelevant for LTI) |
| Parallel | $H_1(z)+H_2(z)$ |
| Negative feedback, forward $G$, feedback $F$ | $\dfrac{G(z)}{1+G(z)F(z)}$ |
| Inverse system | $H^{-1}(z)=1/H(z)$ |

> [!example] Worked Example 17.31 — A digital feedback loop
> A digital integrating controller with a one-sample computation delay has forward path $G(z)=\dfrac{Kz^{-1}}{1-z^{-1}}=\dfrac{K}{z-1}$ in a unity negative-feedback loop. Find the closed-loop system function, the range of $K$ for stability, and the step response for $K=0.5$.
>
> **Step 1.**
> $$
> T(z)=\frac{G}{1+G}=\frac{K/(z-1)}{1+K/(z-1)}=\frac{K}{z-1+K}=\frac{Kz^{-1}}{1-(1-K)z^{-1}}.
> $$
> **Step 2. Stability.** The single closed-loop pole is $z=1-K$; stability requires $\lvert1-K\rvert<1$, i.e. $0<K<2$. The loop is unstable for too much gain, as a sampled control loop always is — the computation delay limits the achievable bandwidth.
>
> **Step 3. Step response** ($K=0.5$, pole $0.5$): $Y=\dfrac{0.5}{z-0.5}\cdot\dfrac{z}{z-1}$, $Y/z=\dfrac{1}{z-1}-\dfrac{1}{z-0.5}$, so $y[n]=\left[1-(0.5)^n\right]u[n]$, with zero steady-state error ($T(1)=1$) thanks to the integrator.
>
> *Check:* simulated outputs $0,\ 0.5,\ 0.75,\ 0.875,\ 0.9375$ ✓. With $K=1$ the pole is at $0$ and the output reaches 1 in one step and stays there (a **deadbeat** response, impossible in continuous time); with $K=1.5$ the pole is at $-0.5$ and the response overshoots with alternating sign: $0,\ 1.5,\ 0.75,\ 1.125,\ 0.9375,\dots$.

### Inverse systems and minimum phase

The inverse of $H(z)$ has the zeros of $H$ as its poles. For the inverse to be **both causal and stable**, those poles — the zeros of $H$ — must lie inside the unit circle.

> [!definition] Minimum-phase system
> A causal, stable rational system is **minimum phase** if all its zeros are also strictly inside the unit circle. Exactly the minimum-phase systems have causal, stable inverses.

For example $H(z)=\dfrac{1-0.5z^{-1}}{1-0.8z^{-1}}$ (zero $0.5$, pole $0.8$) is minimum phase; its inverse $\dfrac{1-0.8z^{-1}}{1-0.5z^{-1}}$ is causal and stable. By contrast $H(z)=1-2z^{-1}$ (zero at $2$) has the inverse $1/(1-2z^{-1})$, which is either causal and unstable (ROC $\lvert z\rvert>2$) or stable and anticausal (ROC $\lvert z\rvert<2$). Channel equalisers in communication receivers face exactly this choice, and solve it by allowing a delay.

> [!warning] Common Pitfall — Pole–zero cancellation hides instability
> Cancelling a factor $(z-p)$ between numerator and denominator is legitimate algebra, but if $\lvert p\rvert\ge1$ and the cancellation occurs between two separately built subsystems (e.g. an unstable plant and a controller zero placed on its pole), the internal signal still contains the growing mode $p^n$. Any mismatch — and there always is one in hardware — lets it appear at the output. Cancel unstable poles on paper only.

> [!tip] Engineering Note — Reading a filter specification sheet
> Coefficients exported by filter-design tools come either as `(b, a)` polynomials, as `(z, p, k)` zeros, poles and gain, or as a matrix of second-order sections `sos`, one row $[b_0,b_1,b_2,1,a_1,a_2]$ per biquad. All three describe the same $H(z)$; the `sos` form is the one to implement.

## 17.7 Discrete-time system analysis

This section collects the analysis questions that the system function answers: Is the system stable? How does it respond to sinusoids? What is its frequency response, and how can it be read from the pole–zero plot? And how does all of this relate to the continuous-time system that a digital system usually imitates or controls?

### BIBO stability

> [!theorem] BIBO stability of LTI systems
> A discrete-time LTI system is bounded-input bounded-output stable if and only if its impulse response is absolutely summable:
> $$
> S=\sum_{n=-\infty}^{\infty}\lvert h[n]\rvert<\infty.
> $$
> Equivalently, the ROC of $H(z)$ contains the unit circle. For a **causal** system with rational $H(z)$: stable $\iff$ every pole satisfies $\lvert p_k\rvert<1$.

*Proof.* **Sufficiency.** If $\lvert x[n]\rvert\le B$ then $\lvert y[n]\rvert=\left\lvert\sum_kh[k]x[n-k]\right\rvert\le B\sum_k\lvert h[k]\rvert=BS$.
**Necessity.** Suppose $S=\infty$. For real $h$, choose the bounded input $x[k]=\operatorname{sgn}(h[-k])$ (with $\operatorname{sgn}0=0$). Then $y[0]=\sum_kh[k]x[-k]=\sum_k\lvert h[k]\rvert=\infty$, so a bounded input produces an unbounded output.
**Pole form.** A causal rational $h[n]$ is a finite sum of terms $c\,n^mp^nu[n]$. Each is absolutely summable iff $\lvert p\rvert<1$ (the ratio test: $\lvert(n+1)^mp^{n+1}\rvert/\lvert n^mp^n\rvert\to\lvert p\rvert$). $\square$

Three cases are worth naming for causal systems:

| Pole location | Natural response | Classification |
|---|---|---|
| all $\lvert p\rvert<1$ | decays geometrically | **stable** |
| simple poles on $\lvert z\rvert=1$, none outside | bounded, non-decaying (constant, $(-1)^n$, sinusoid) | **marginally stable** — *not* BIBO stable |
| any $\lvert p\rvert>1$, or repeated poles on $\lvert z\rvert=1$ | grows | **unstable** |

The accumulator $H(z)=1/(1-z^{-1})$ is the canonical marginal case: its impulse response $u[n]$ is bounded, but the bounded step input produces the unbounded ramp $(n+1)u[n]$ — the double pole of $Y=z^2/(z-1)^2$ on the unit circle.

### Stability tests without factoring

For high-order denominators, root finding is avoidable. Write the denominator in positive powers as $D(z)=z^N+a_1z^{N-1}+\dots+a_N$.

> [!theorem] Necessary conditions (any order)
> If all roots of the real monic polynomial $D(z)$ lie strictly inside the unit circle, then
> $$
> D(1)>0,\qquad(-1)^ND(-1)>0,\qquad\lvert a_N\rvert<1.
> $$

*Reason.* $D(1)=\prod(1-p_k)$; for real roots in $(-1,1)$ every factor is positive, and conjugate pairs contribute $\lvert1-p\rvert^2>0$. Likewise $(-1)^ND(-1)=\prod(1+p_k)>0$. And $\lvert a_N\rvert=\prod\lvert p_k\rvert<1$. $\square$

> [!theorem] Second-order stability triangle
> The polynomial $z^2+a_1z+a_2$ has both roots strictly inside the unit circle if and only if
> $$
> \lvert a_2\rvert<1\qquad\text{and}\qquad\lvert a_1\rvert<1+a_2.
> $$
> For $N=2$ the necessary conditions above are also sufficient. The roots are complex for $a_2>a_1^2/4$, in which case their radius is $\sqrt{a_2}$.

```desmos-graph
left=-2.5; right=2.5; bottom=-1.5; top=1.5
---
y<1|blue
y>x-1|blue
y>-x-1|blue
y=\frac{x^2}{4}|red|dashed
(-0.9,0.81)|label:resonator 17.28
(-0.75,0.125)|label:Ex. 17.24
```

*The stability triangle in the $(a_1,a_2)$ plane for $z^2+a_1z+a_2$: the region where all three shadings overlap is stable. Above the dashed parabola $a_2=a_1^2/4$ the poles are complex with radius $\sqrt{a_2}$; below it they are real. Both systems marked are stable.*

For third and higher order, the **Jury** (or **Schur–Cohn**) table extends the test by a recursive construction; for a cubic $z^3+a_1z^2+a_2z+a_3$ it adds to the three necessary conditions the single extra inequality

$$
\left\lvert a_3^2-1\right\rvert>\left\lvert a_3a_1-a_2\right\rvert.
$$

The **bilinear map** $z=(1+w)/(1-w)$ sends the unit disc onto the left half $w$-plane, so applying the Routh–Hurwitz test of [Ch. 14](14-laplace-transforms.md) to the transformed polynomial is an alternative route.

> [!example] Worked Example 17.32 — Testing cubics
> (a) $D(z)=z^3-1.2z^2+0.5z-0.1$. Conditions: $D(1)=1-1.2+0.5-0.1=0.2>0$ ✓; $-D(-1)=-(-1-1.2-0.5-0.1)=2.8>0$ ✓; $\lvert a_3\rvert=0.1<1$ ✓; $\lvert0.01-1\rvert=0.99>\lvert0.12-0.5\rvert=0.38$ ✓. **Stable.** *Check:* the roots are $0.682$ and $0.259\pm j0.282$ (radius $0.383$).
>
> (b) $D(z)=z^3-1.5z^2+1.1z-0.6$. $D(1)=1-1.5+1.1-0.6=0$: the first condition fails with equality, so there is a root *on* the unit circle at $z=1$. Dividing it out, $D=(z-1)(z^2-0.5z+0.6)$, whose quadratic factor has roots of radius $\sqrt{0.6}=0.775$. **Marginally stable** (a simple pole at $z=1$ makes the system an accumulator in disguise).

### Frequency response

> [!theorem] Sinusoidal steady state
> Let a stable LTI system with real impulse response have frequency response $H(e^{j\omega})=\lvert H(e^{j\omega})\rvert e^{j\angle H(e^{j\omega})}$. Then:
> 1. The eternal input $x[n]=e^{j\omega_0n}$ produces exactly $y[n]=H(e^{j\omega_0})\,e^{j\omega_0n}$.
> 2. The eternal input $x[n]=A\cos(\omega_0n+\phi)$ produces
> $$
> y[n]=A\left\lvert H(e^{j\omega_0})\right\rvert\cos\!\left(\omega_0n+\phi+\angle H(e^{j\omega_0})\right).
> $$
> 3. A causal input $A\cos(\omega_0n+\phi)u[n]$ produces this response plus a transient made of the system's natural modes, which decays because the poles are inside the unit circle.

*Derivation.* Statement 1 is the eigenfunction property of §17.2 with $z=e^{j\omega_0}$ — the sum converges because the unit circle is in the ROC. For 2, write the cosine as $\tfrac A2\left(e^{j\phi}e^{j\omega_0n}+e^{-j\phi}e^{-j\omega_0n}\right)$ and use $H(e^{-j\omega_0})=H^*(e^{j\omega_0})$ for real $h$. For 3, the partial-fraction expansion of $Y(z)=H(z)X(z)$ separates the terms with the input's poles $e^{\pm j\omega_0}$ (the steady state) from those with the poles of $H$ (the transient). $\square$

Properties of $H(e^{j\omega})$ that follow from §17.2:

- **Periodicity:** $H(e^{j(\omega+2\pi)})=H(e^{j\omega})$; only $-\pi<\omega\le\pi$ needs plotting.
- **Symmetry** (real $h$): $\lvert H\rvert$ is even and $\angle H$ is odd, so $0\le\omega\le\pi$ suffices. In physical units $\omega=\pi$ corresponds to the Nyquist frequency $f_s/2$.
- **DC and Nyquist gains:** $H(e^{j0})=H(1)=\sum h[n]$ and $H(e^{j\pi})=H(-1)=\sum(-1)^nh[n]$.
- **Group delay:** $\tau_g(\omega)=-\dfrac{d\angle H(e^{j\omega})}{d\omega}$, in samples. A **linear-phase** FIR filter with symmetric coefficients $b_k=b_{M-k}$ has $\angle H=-\omega M/2$ (plus jumps of $\pi$ at zeros), i.e. a constant delay of $M/2$ samples for all frequencies.

> [!example] Worked Example 17.33 — First-order low-pass: the frequency response
> The "exponential smoother" $y[n]=a\,y[n-1]+(1-a)\,x[n]$, $0<a<1$, has $H(z)=\dfrac{1-a}{1-az^{-1}}$.
>
> **Step 1. Magnitude.** On $z=e^{j\omega}$: $\lvert1-ae^{-j\omega}\rvert^2=(1-a\cos\omega)^2+a^2\sin^2\omega=1-2a\cos\omega+a^2$, so
> $$
> \lvert H(e^{j\omega})\rvert=\frac{1-a}{\sqrt{1-2a\cos\omega+a^2}},\qquad\angle H(e^{j\omega})=-\arctan\frac{a\sin\omega}{1-a\cos\omega}.
> $$
> **Step 2. Values for $a=0.8$.** DC: $\lvert H\rvert=0.2/0.2=1$. Nyquist ($\omega=\pi$): $0.2/1.8=0.111$, i.e. $-19.1$ dB.
>
> **Step 3. 3-dB frequency.** Setting $\lvert H\rvert^2=\tfrac12$: $1-2a\cos\omega_c+a^2=2(1-a)^2$, so
> $$
> \cos\omega_c=1-\frac{(1-a)^2}{2a}=1-\frac{0.04}{1.6}=0.975,\qquad\omega_c=0.2241\ \text{rad/sample}.
> $$
> At $f_s=10$ kHz this is $f_c=\omega_cf_s/2\pi=357$ Hz.

```desmos-graph
left=-3.1416; right=3.1416; bottom=0; top=1.15
---
y=\frac{0.5}{\sqrt{1.25-\cos(x)}}|green
y=\frac{0.2}{\sqrt{1.64-1.6\cos(x)}}|blue
y=\frac{0.05}{\sqrt{1.9025-1.9\cos(x)}}|red
y=0.7071|grey|dashed
```

*Magnitude responses $\lvert H(e^{j\omega})\rvert=(1-a)/\sqrt{1-2a\cos\omega+a^2}$ of the first-order smoother for $a=0.5$ (green), $0.8$ (blue) and $0.95$ (red), plotted over one period $-\pi<\omega\le\pi$. As the pole approaches $z=1$ the pass-band narrows; the dashed line is the $-3$ dB level.*

### Geometric evaluation from the pole–zero plot

Write $H$ in factored form and set $z=e^{j\omega}$:

$$
\lvert H(e^{j\omega})\rvert=\lvert b_0\rvert\,\frac{\prod_k\lvert e^{j\omega}-q_k\rvert}{\prod_k\lvert e^{j\omega}-p_k\rvert},\qquad
\angle H(e^{j\omega})=\angle b_0+(N-M)\omega+\sum_k\angle\left(e^{j\omega}-q_k\right)-\sum_k\angle\left(e^{j\omega}-p_k\right).
$$

Each factor $e^{j\omega}-c$ is the vector from the pole or zero $c$ to the point $e^{j\omega}$ on the unit circle. As $\omega$ sweeps from $0$ to $\pi$, the point travels along the upper half of the circle, and:

- the **magnitude** is (gain) × (product of zero-vector lengths) / (product of pole-vector lengths);
- passing close to a **pole** makes a denominator length small: a **peak**, sharper the closer the pole is to the circle;
- passing through a **zero** on the circle makes the magnitude exactly **zero**;
- poles and zeros at the origin have vectors of length 1 for all $\omega$: they affect only the phase (a pure delay).

For a pole at radius $r$ close to 1, the 3-dB bandwidth of the peak it creates is approximately

$$
\Delta\omega_{3\,\text{dB}}\approx2(1-r)\quad\text{rad/sample},
$$

because near the pole the unit circle is almost straight, and the pole-vector length grows by $\sqrt2$ when the point has moved a distance $1-r$ either side of the closest approach. This rule makes pole placement a design tool.

> [!example] Worked Example 17.34 — The resonator of Example 17.28
> Analyse $H(z)=\dfrac{z^2-1}{z^2-0.9z+0.81}$, poles $0.9e^{\pm j\pi/3}$, zeros $\pm1$.
>
> **Step 1. Geometric reading.** At $\omega=0$ and $\omega=\pi$ the point $e^{j\omega}$ sits on a zero: $\lvert H\rvert=0$. Near $\omega=\pi/3$ it passes within $1-r=0.1$ of a pole: a peak, of approximate width $2(0.1)=0.2$ rad.
>
> **Step 2. Peak value.** At $\omega=\pi/3$: the zero vectors have lengths $\lvert e^{j\pi/3}-1\rvert=1$ and $\lvert e^{j\pi/3}+1\rvert=\sqrt3$; the near pole vector has length $0.1$; the far pole vector is $\lvert e^{j\pi/3}-0.9e^{-j\pi/3}\rvert=\lvert0.05+j1.6454\rvert=1.6462$. So
> $$
> \lvert H(e^{j\pi/3})\rvert=\frac{1\cdot\sqrt3}{0.1\times1.6462}=10.52.
> $$
> **Step 3. Steady-state response.** For $x[n]=\cos(\pi n/3)$, a numerical evaluation gives $H(e^{j\pi/3})=10.52\,e^{j0.030}$, so in steady state $y[n]=10.52\cos(\pi n/3+0.030)$. The exact maximum of $\lvert H\rvert$ is $10.53$, at $\omega=0.334\pi$, very slightly above $\pi/3$ (the far pole and the zeros tilt the peak).

```desmos-graph
left=0; right=3.1416; bottom=0; top=12
---
y=\frac{2\left|\sin(x)\right|}{\sqrt{\left(\cos(2x)-0.9\cos(x)+0.81\right)^2+\left(\sin(2x)-0.9\sin(x)\right)^2}}|blue
x=1.0472|grey|dashed
(1.0472,10.52)|label:peak ≈ 10.5 at π/3
```

*Magnitude response of the resonator $H(z)=(1-z^{-2})/(1-0.9z^{-1}+0.81z^{-2})$ for $0\le\omega\le\pi$: zeros at $z=\pm1$ force nulls at DC and at Nyquist; the poles at $0.9e^{\pm j\pi/3}$ produce a peak at $\omega\approx\pi/3$ with bandwidth $\approx0.2$ rad.*

> [!example] Worked Example 17.35 — Designing a 50 Hz notch filter by pole–zero placement
> Remove 50 Hz mains hum from an ECG sampled at $f_s=1$ kHz, with a notch about 15 Hz wide.
>
> **Step 1. Notch frequency.** $\omega_0=2\pi\cdot50/1000=0.1\pi$ rad/sample.
>
> **Step 2. Zeros on the circle** at $e^{\pm j\omega_0}$ give exact nulls: numerator $1-2\cos\omega_0z^{-1}+z^{-2}$.
>
> **Step 3. Poles just inside, at the same angle,** $re^{\pm j\omega_0}$, restore unit gain away from $\omega_0$ (each pole nearly cancels the effect of its zero except near the notch). Bandwidth $2(1-r)\approx2\pi\cdot15/1000=0.094$ rad gives $r\approx0.95$:
> $$
> H(z)=\frac{1-1.9021z^{-1}+z^{-2}}{1-1.8070z^{-1}+0.9025z^{-2}}.
> $$
> **Step 4. Check.** Numerically, $\lvert H\rvert=0$ at $50$ Hz, $1.025$ at DC and $1.052$ at Nyquist (scale by $1/1.025$ if exact unit DC gain is wanted), and the $-3$ dB notch width is $0.0936$ rad $=14.9$ Hz. The poles have radius $0.95$: stable.

> [!example] Worked Example 17.36 — Moving average and linear phase
> For the 4-point moving average of Example 17.29,
> $$
> H(e^{j\omega})=\frac14\,\frac{1-e^{-j4\omega}}{1-e^{-j\omega}}=\frac14\,\frac{e^{-j2\omega}\left(2j\sin2\omega\right)}{e^{-j\omega/2}\left(2j\sin\frac{\omega}{2}\right)}=e^{-j3\omega/2}\,\frac{\sin2\omega}{4\sin(\omega/2)}.
> $$
> The magnitude is the **Dirichlet kernel** $\left\lvert\dfrac{\sin2\omega}{4\sin(\omega/2)}\right\rvert$, with zeros at $\omega=\pi/2,\pi$; the phase is $-1.5\omega$ (plus $\pi$ jumps where the real factor changes sign) — **linear phase**, a constant group delay of $M/2=1.5$ samples, because the coefficients are symmetric.

```desmos-graph
left=-3.1416; right=3.1416; bottom=0; top=1.1
---
y=\left|\frac{\sin(2x)}{4\sin(x/2)}\right||blue
y=\left|\frac{\sin(4x)}{8\sin(x/2)}\right||red|dashed
```

*Magnitude responses of the 4-point (blue) and 8-point (red, dashed) moving averages: crude low-pass filters with nulls at multiples of $2\pi/M$. Doubling the length halves the main-lobe width.*

### All-pass systems

The geometric picture shows how to build a filter with *flat* magnitude: place each zero at the **conjugate reciprocal** $1/p^*$ of a pole $p$. Then for $\lvert z\rvert=1$,

$$
H_{\text{ap}}(z)=\frac{z^{-1}-p^*}{1-pz^{-1}},\qquad\lvert H_{\text{ap}}(e^{j\omega})\rvert=\frac{\lvert e^{-j\omega}-p^*\rvert}{\lvert1-pe^{-j\omega}\rvert}=\frac{\lvert e^{-j\omega}\rvert\,\lvert1-p^*e^{j\omega}\rvert}{\lvert1-pe^{-j\omega}\rvert}=1,
$$

since the last numerator is the conjugate of the denominator. All-pass sections change only the phase; they are used as phase equalisers, as building blocks of fractional delays, and in the theorem that **every rational stable system factors as (minimum phase) × (all-pass)**: reflecting a zero from outside the circle to its conjugate reciprocal inside leaves $\lvert H\rvert$ unchanged and reduces the phase lag, hence "minimum phase".

### Relation to the Laplace transform: $z=e^{sT}$

Suppose $x[n]=x_a(nT)$ are samples of a continuous signal. Represent the sampled signal, as in [§16.8](16-fourier-transforms.md#168-sampling-related-mathematical-foundations), by the impulse train

$$
x_s(t)=\sum_{n}x_a(nT)\,\delta(t-nT).
$$

Its Laplace transform, using $\mathcal{L}\{\delta(t-nT)\}=e^{-snT}$, is

$$
X_s(s)=\sum_nx[n]\,e^{-snT}=\left.X(z)\right|_{z=e^{sT}}.
$$

> [!theorem] The $s$-to-$z$ map
> The Z-transform of the samples equals the Laplace transform of the impulse-sampled signal under the substitution
> $$
> z=e^{sT},\qquad s=\sigma+j\Omega\ \Longrightarrow\ \lvert z\rvert=e^{\sigma T},\quad\angle z=\Omega T.
> $$
> In particular:
> - the imaginary axis $\sigma=0$ maps onto the unit circle, with digital frequency $\omega=\Omega T$;
> - the open left half-plane $\sigma<0$ maps into the open unit disc $\lvert z\rvert<1$; the right half-plane maps outside it;
> - a continuous pole at $s=-\alpha\pm j\Omega_d$ appears in the sampled signal as a pole at $z=e^{-\alpha T}e^{\pm j\Omega_dT}$;
> - the map is **not one-to-one**: $s$ and $s+jk\,2\pi/T$ give the same $z$. Each horizontal strip of height $\Omega_s=2\pi/T$ covers the whole $z$-plane once.

The last statement is aliasing seen from the $s$-plane: frequencies $\Omega$ and $\Omega+k\Omega_s$ are indistinguishable after sampling, and so are poles that differ by $jk\Omega_s$. The **primary strip** $-\pi/T<\Omega\le\pi/T$ is the set of continuous frequencies that sampling represents faithfully, and its edges $\Omega=\pm\pi/T$ map to the negative real axis $z<0$ — the Nyquist frequency.

```desmos-graph
left=-1.6; right=1.6; bottom=-1.6; top=1.6
---
x^2+y^2=1|black
x^2+y^2=0.3679|blue|dashed
x^2+y^2=1.8221|red|dashed
y=x|green|0<x<1.6
y=-x|green|0<x<1.6
y=0|purple|-1.6<x<0
```

*Images under $z=e^{sT}$: the $j\Omega$ axis maps to the unit circle (black); the vertical lines $\sigma T=-0.5$ and $\sigma T=+0.3$ map to circles of radius $e^{-0.5}=0.607$ (blue) and $e^{0.3}=1.35$ (red); the horizontal lines $\Omega T=\pm\pi/4$ map to rays at $\pm45^\circ$ (green); $\Omega T=\pm\pi$ maps to the negative real axis (purple).*

Useful landmarks of the map: constant damping $\sigma$ ↔ circles; constant frequency $\Omega$ ↔ rays; the point $s=0$ (DC) ↔ $z=1$; $s=-\infty$ ↔ $z=0$. Fast continuous poles (large $\alpha T$) crowd towards the origin; slow ones (small $\alpha T$) crowd towards $z=1$, which is why high sampling rates push all poles of a digital controller into a small, numerically delicate neighbourhood of $z=1$.

> [!example] Worked Example 17.37 — Mapping a damped oscillation
> A continuous-time system has poles $s=-200\pm j1000$ s⁻¹ and is sampled at $f_s=4$ kHz ($T=0.25$ ms). Where are the poles of the sampled system, and what are its digital damping and frequency?
>
> **Step 1.** $\sigma T=-200\times0.25\times10^{-3}=-0.05$, $\Omega T=1000\times0.25\times10^{-3}=0.25$ rad.
>
> **Step 2.** $z=e^{-0.05}e^{\pm j0.25}=0.9512\,e^{\pm j0.25}=0.9217\pm j0.2353$.
>
> **Step 3.** The pole radius $0.951<1$ confirms stability (left half-plane → inside the circle). The digital oscillation frequency $0.25$ rad/sample corresponds to $0.25\,f_s/2\pi=159$ Hz $=1000/2\pi$ Hz, as it must. If the same system were sampled at only $f_s=150$ Hz, the pole angle $\Omega T=6.667$ rad would wrap to $6.667-2\pi=0.384$ rad — an alias at $9.2$ Hz.

### Discretising continuous-time systems

The substitution $z=e^{sT}$ cannot be applied directly to a rational $H(s)$ — it would produce a non-rational $H(z)$. Practical discretisation methods are rational approximations that preserve different features.

**Impulse invariance.** Sample the analogue impulse response: $h[n]=T\,h_a(nT)$. For $H_a(s)=\sum_k\dfrac{A_k}{s-s_k}$, $h_a(t)=\sum_kA_ke^{s_kt}u(t)$, so

$$
H(z)=\sum_k\frac{T A_k}{1-e^{s_kT}z^{-1}}.
$$

Each analogue pole $s_k$ maps exactly to $e^{s_kT}$, but the frequency response is the aliased sum $H(e^{j\omega})=\sum_kH_a\!\left(j\frac{\omega+2\pi k}{T}\right)$, so the method suits only band-limited (low-pass, band-pass) prototypes.

**Bilinear transform.** Replace $s$ by the first-order rational approximation of $s=\frac1T\ln z$ (from the trapezoidal rule of integration):

$$
s=\frac{2}{T}\,\frac{1-z^{-1}}{1+z^{-1}},\qquad z=\frac{1+sT/2}{1-sT/2}.
$$

It maps the left half-plane *exactly* onto the unit disc (so stable analogue filters always give stable digital ones) and the whole $j\Omega$ axis onto the unit circle *once* — no aliasing. The price is **frequency warping**: with $z=e^{j\omega}$,

$$
s=\frac2T\,\frac{1-e^{-j\omega}}{1+e^{-j\omega}}=j\,\frac2T\tan\frac{\omega}{2}\quad\Longrightarrow\quad\Omega=\frac{2}{T}\tan\frac{\omega}{2},\qquad\omega=2\arctan\frac{\Omega T}{2}.
$$

Critical frequencies are therefore **pre-warped**: to place a digital cutoff at $\omega_c$, design the analogue prototype with cutoff $\Omega_c=\frac2T\tan\frac{\omega_c}{2}$.

```desmos-graph
left=0; right=6; bottom=0; top=3.5
---
y=2\arctan\left(\frac{x}{2}\right)|blue
y=x|grey|dashed|0<x<3.1416
y=3.1416|red|dotted
```

*Bilinear frequency warping, $\omega=2\arctan(\Omega T/2)$ (blue) against the ideal $\omega=\Omega T$ (dashed): the two agree at low frequency, but the whole infinite analogue axis is compressed into $0\le\omega<\pi$ (red dotted line).*

> [!example] Worked Example 17.38 — Impulse invariance versus bilinear for an RC low-pass
> Discretise $H_a(s)=\dfrac{1}{1+s\tau}$ with $\tau=1$ ms and $T=0.25$ ms.
>
> **Impulse invariance.** $h_a(t)=\frac1\tau e^{-t/\tau}u(t)$, so $h[n]=\frac{T}{\tau}\left(e^{-T/\tau}\right)^nu[n]$ and
> $$
> H_{\text{ii}}(z)=\frac{0.25}{1-0.7788z^{-1}}.
> $$
> The pole $e^{-T/\tau}=0.7788$ is exact (compare the Euler pole $0.75$ of Example 17.27), but the DC gain is $0.25/(1-0.7788)=1.130$, not 1: the analogue response is not band-limited, so aliased spectral copies add to the gain. In practice one rescales by $1/1.130$.
>
> **Bilinear.** With $\frac{2}{T}=8000$ s⁻¹ and $\tau=10^{-3}$ s: $1+s\tau=1+8\,\dfrac{1-z^{-1}}{1+z^{-1}}=\dfrac{9-7z^{-1}}{1+z^{-1}}$, so
> $$
> H_{\text{bl}}(z)=\frac{1+z^{-1}}{9-7z^{-1}}=\frac{0.1111\,(1+z^{-1})}{1-0.7778z^{-1}}.
> $$
> DC gain $H(1)=2/2=1$ exactly; $H(-1)=0$ — the Nyquist frequency is a transmission zero, the image of $s=\infty$.

> [!example] Worked Example 17.39 — Bilinear design with pre-warping
> Design a first-order digital low-pass with 3-dB cutoff $f_c=1$ kHz at $f_s=8$ kHz.
>
> **Step 1. Digital cutoff.** $\omega_c=2\pi f_c/f_s=\pi/4$.
>
> **Step 2. Pre-warp.** $\Omega_c=\dfrac2T\tan\dfrac{\omega_c}{2}=16000\tan(\pi/8)=6627$ rad/s (the naive value would be $2\pi\cdot1000=6283$ rad/s). Then $\Omega_cT=2\tan(\pi/8)=0.8284$.
>
> **Step 3. Prototype** $H_a(s)=\dfrac{\Omega_c}{s+\Omega_c}$. Substituting the bilinear map and multiplying through by $T(1+z^{-1})$:
> $$
> H(z)=\frac{\Omega_cT\,(1+z^{-1})}{(2+\Omega_cT)+(\Omega_cT-2)z^{-1}}=\frac{0.2929\,(1+z^{-1})}{1-0.4142z^{-1}}.
> $$
> **Step 4. Check.** $\lvert H(e^{j\pi/4})\rvert=0.7071$ exactly, i.e. $-3$ dB at 1 kHz ✓; $H(1)=0.5858/0.5858=1$ ✓. Without pre-warping the cutoff would land at $\omega=2\arctan(6283/16000)=0.748$ rad, i.e. $952$ Hz.

> [!warning] Common Pitfall — Confusing rad/sample with rad/s
> $\omega$ (digital, rad/sample, periodic with period $2\pi$) and $\Omega$ (analogue, rad/s) differ by the factor $T$. A digital filter "with cutoff $0.2$" is meaningless until $f_s$ is known, and design software variously normalises frequency to $\pi$ rad/sample, to $f_s/2$ ("1.0 = Nyquist") or to $f_s$. Always convert with $\omega=2\pi f/f_s$ and check that $\omega_c<\pi$.

> [!tip] Engineering Note — What the unit-circle picture buys a designer
> Every specification of a digital filter translates into a pole–zero picture: stopband frequencies become zeros on the circle, passband peaks become poles near the circle at the corresponding angle, the sharpness of a feature becomes the distance $1-r$, and stability is the single constraint $r<1$. The magnitude plot is then the shadow of these poles and zeros on the circle. Sophisticated design methods (Butterworth, Chebyshev, elliptic prototypes transformed bilinearly; Parks–McClellan FIR design) are simply systematic ways of choosing those positions.

## Summary

| Concept | Formula |
|---|---|
| Discrete frequency | $\omega=\Omega T=2\pi f/f_s$ rad/sample; $e^{j(\omega+2\pi)n}=e^{j\omega n}$ |
| Convolution sum | $y[n]=\sum_kx[k]h[n-k]$ |
| Bilateral Z-transform | $X(z)=\sum_{n=-\infty}^{\infty}x[n]z^{-n}$, with ROC |
| Unilateral Z-transform | $\mathcal{X}(z)=\sum_{n=0}^{\infty}x[n]z^{-n}$ |
| Geometric pairs | $a^nu[n]\leftrightarrow\frac{z}{z-a}$, $\lvert z\rvert>\lvert a\rvert$; $-a^nu[-n-1]\leftrightarrow\frac{z}{z-a}$, $\lvert z\rvert<\lvert a\rvert$ |
| Repeated pole | $na^nu[n]\leftrightarrow\frac{az}{(z-a)^2}$ |
| Damped cosine | $r^n\cos(\omega_0n)u[n]\leftrightarrow\frac{z(z-r\cos\omega_0)}{z^2-2rz\cos\omega_0+r^2}$ |
| Shift / scaling / multiply by $n$ | $z^{-n_0}X(z)$; $a^nx[n]\leftrightarrow X(z/a)$; $nx[n]\leftrightarrow-zX'(z)$ |
| Reversal / convolution | $x[-n]\leftrightarrow X(1/z)$; $x*h\leftrightarrow XH$ |
| Initial / final value | $x[0]=\lim_{z\to\infty}X$; $x[\infty]=\lim_{z\to1}(z-1)X$ |
| Parseval | $\sum\lvert x[n]\rvert^2=\frac1{2\pi}\int_{-\pi}^{\pi}\lvert X(e^{j\omega})\rvert^2d\omega$ |
| ROC shapes | right-sided: $\lvert z\rvert>r_1$; left-sided: $\lvert z\rvert<r_2$; two-sided: $r_1<\lvert z\rvert<r_2$; finite: all $z$ (except $0,\infty$) |
| Inversion integral | $x[n]=\frac1{2\pi j}\oint_CX(z)z^{n-1}dz=\sum\operatorname{Res}\left[Xz^{n-1}\right]$ |
| Partial fractions | expand $X(z)/z$; $A_i=(z-p_i)X(z)/z\,\rvert_{p_i}$ |
| Unilateral delay | $y[n-1]\leftrightarrow z^{-1}\mathcal{Y}+y[-1]$; $y[n-2]\leftrightarrow z^{-2}\mathcal{Y}+z^{-1}y[-1]+y[-2]$ |
| Unilateral advance | $y[n+1]\leftrightarrow z\mathcal{Y}-zy[0]$ |
| System function | $H(z)=\dfrac{\sum b_kz^{-k}}{\sum a_kz^{-k}}=\dfrac{Y}{X}$ |
| Feedback | $\dfrac{G}{1+GF}$ |
| Stability (causal) | all poles $\lvert p_k\rvert<1$ $\iff\sum\lvert h[n]\rvert<\infty$ |
| Second-order test | $z^2+a_1z+a_2$: $\lvert a_2\rvert<1$, $\lvert a_1\rvert<1+a_2$ |
| Frequency response | $H(e^{j\omega})$; $\cos\omega_0n\mapsto\lvert H\rvert\cos(\omega_0n+\angle H)$ |
| Pole bandwidth | $\Delta\omega_{3\,\text{dB}}\approx2(1-r)$ |
| $s$–$z$ map | $z=e^{sT}$: LHP $\to\lvert z\rvert<1$, $j\Omega$ axis $\to$ unit circle |
| Impulse invariance | $\frac{A}{s-s_k}\to\frac{TA}{1-e^{s_kT}z^{-1}}$ |
| Bilinear transform | $s=\frac2T\frac{1-z^{-1}}{1+z^{-1}}$, $\Omega=\frac2T\tan\frac\omega2$ |

## Further reading

- A. V. Oppenheim & R. W. Schafer, *Discrete-Time Signal Processing*, 3rd ed., Pearson — Chs. 2–5: sequences, the Z-transform, ROC, system functions, structures and the bilinear transform; the standard reference.
- A. V. Oppenheim, A. S. Willsky & S. H. Nawab, *Signals and Systems*, 2nd ed., Prentice Hall — Ch. 10 develops the Z-transform in parallel with the Laplace transform of Ch. 9.
- J. G. Proakis & D. G. Manolakis, *Digital Signal Processing: Principles, Algorithms and Applications*, 4th ed., Pearson — many worked inversions and the unilateral transform for difference equations.
- K. Ogata, *Discrete-Time Control Systems*, 2nd ed., Prentice Hall — the $s$–$z$ mapping, Jury stability test and sampled-data feedback loops.
