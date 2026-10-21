# Chapter 16 — Fourier Transforms

[← Previous](15-fourier-series.md) · [Table of Contents](../00-front/table-of-contents.md) · [Next →](17-z-transforms.md)

> **Prerequisites:** Ch. 2 (complex numbers), Ch. 5 (integral calculus), Ch. 14 (Laplace transforms), Ch. 15 (Fourier series and Fourier analysis)
> **Learning outcomes:** After this chapter you can
> - compute forward and inverse Fourier transforms of standard signals, including generalized transforms of impulses, steps, constants and sinusoids;
> - derive and apply the transform properties (linearity, shift, scaling, duality, differentiation, integration, modulation) so that most transforms are found without integrating;
> - use the convolution and multiplication theorems to analyse linear time-invariant systems and products of signals in the frequency domain;
> - interpret magnitude and phase spectra, frequency responses, group delay and bandwidth, including the time–bandwidth (uncertainty) trade-off;
> - analyse periodic and finite-length signals spectrally, including windowing, leakage and resolution;
> - compute energy and energy spectral density with Parseval's theorem, and state the mathematics of ideal sampling: impulse-train spectra, aliasing, the Nyquist condition and sinc reconstruction.

The Fourier series of [Ch. 15](15-fourier-series.md) represents a *periodic* signal as a sum of harmonics at the discrete frequencies $k\omega_0$. Most signals an engineer meets are not periodic: a single radar pulse, a transient on a power line, a spoken word, a step applied to an amplifier, a bit of data on a cable. To describe such signals in terms of frequency we let the period grow without bound. The harmonics crowd together, their spacing $\omega_0$ shrinks to $d\omega$, the sum becomes an integral, and we obtain the **Fourier transform**: a representation of an *aperiodic* signal as a continuum of complex exponentials $e^{j\omega t}$.

The Fourier transform sits between the two neighbouring chapters. It is the Laplace transform of [Ch. 14](14-laplace-transforms.md) restricted to the imaginary axis $s=j\omega$ (whenever that axis lies in the region of convergence), and when a signal is sampled it becomes, after a change of variable, the discrete-time Fourier transform and the Z-transform of [Ch. 17](17-z-transforms.md) on the unit circle. Compared with the Laplace transform it gives up the ability to handle growing signals directly, but gains a perfectly symmetric inverse, a clean physical meaning (the frequency content of the signal), and an energy theorem (Parseval) that makes "how much of the signal lies in this band?" a well-posed question.

Throughout this chapter $x(t)$ is a real or complex function of continuous time $t$ (in seconds), $\omega$ is angular frequency in rad/s, and $f=\omega/2\pi$ is frequency in Hz. Wherever a voltage example helps, $x(t)$ is a voltage in volts.

## 16.1 Continuous Fourier transform

### Motivation: from series to integral

A sinusoid $\cos(\omega t)$ passing through a linear circuit comes out as a sinusoid of the same frequency, scaled and phase-shifted; that is the whole content of phasor analysis. If we can write an arbitrary input as a superposition of sinusoids — or, more conveniently, of complex exponentials $e^{j\omega t}$ — then the response of any linear time-invariant circuit to that input is the same superposition of scaled exponentials. For periodic inputs, the Fourier series supplies the superposition. For aperiodic inputs we need a *continuous* superposition, and the weight given to each frequency is the Fourier transform.

**Heuristic derivation.** Let $x(t)$ be zero outside $\lvert t\rvert<T_0/2$, and let $\tilde x(t)$ be its periodic extension with period $T_0$ and fundamental $\omega_0=2\pi/T_0$. The exponential Fourier coefficients of $\tilde x$ ([Ch. 15](15-fourier-series.md)) are

$$
c_k=\frac{1}{T_0}\int_{-T_0/2}^{T_0/2}x(t)\,e^{-jk\omega_0t}\,dt .
$$

Define the function

$$
X(\omega)=\int_{-\infty}^{\infty}x(t)\,e^{-j\omega t}\,dt ,
$$

which, because $x$ vanishes outside one period, gives $c_k=X(k\omega_0)/T_0$. The series becomes

$$
\tilde x(t)=\sum_{k=-\infty}^{\infty}\frac{1}{T_0}X(k\omega_0)\,e^{jk\omega_0t}
=\frac{1}{2\pi}\sum_{k=-\infty}^{\infty}X(k\omega_0)\,e^{jk\omega_0t}\,\omega_0 .
$$

Now let $T_0\to\infty$. Then $\tilde x\to x$ (the neighbouring copies move off to infinity), the spacing $\omega_0$ becomes the differential $d\omega$, the discrete frequencies $k\omega_0$ become a continuous variable $\omega$, and the sum is a Riemann sum for

$$
x(t)=\frac{1}{2\pi}\int_{-\infty}^{\infty}X(\omega)\,e^{j\omega t}\,d\omega .
$$

Two facts come out of this argument and are worth keeping:

1. **The Fourier series coefficients of a periodic signal are samples of the transform of one period**, divided by the period: $c_k=X(k\omega_0)/T_0$.
2. As $T_0$ grows each individual coefficient $c_k$ shrinks like $1/T_0$, but the *density* $T_0c_k$ stays finite. The transform is that density — it is an amplitude *per unit frequency*, not an amplitude.

### Definition

> [!definition] Fourier transform (angular-frequency convention)
> For a signal $x(t)$, the Fourier transform is
>
> $$
> X(\omega) = \mathcal{F}\{x(t)\} = \int_{-\infty}^{\infty} x(t)\,e^{-j\omega t}\,dt ,
> $$
>
> where $\omega$ is angular frequency in rad/s. We write $x(t) \leftrightarrow X(\omega)$ and call $x$ and $X$ a **transform pair**.

$X(\omega)$ is in general complex even when $x(t)$ is real, and we write it in polar form $X(\omega) = \lvert X(\omega)\rvert e^{j\angle X(\omega)}$. If $x(t)$ is a voltage in V, then $X(\omega)$ has units V·s = V/Hz: it is a **spectral density**. Two immediate special values are useful checks on every computation:

$$
X(0)=\int_{-\infty}^{\infty}x(t)\,dt\quad(\text{the area of } x),\qquad
x(0)=\frac{1}{2\pi}\int_{-\infty}^{\infty}X(\omega)\,d\omega\quad(\text{the area of } X\text{, divided by }2\pi).
$$

The second follows from the inverse formula (Section 16.2) at $t=0$.

**Relation to the Laplace transform.** The bilateral Laplace transform is $X_L(s)=\int x(t)e^{-st}dt$ with $s=\sigma+j\omega$. Setting $\sigma=0$ gives the Fourier transform: $X(\omega)=X_L(j\omega)$, *provided the region of convergence (ROC) of $X_L$ contains the imaginary axis*. For a causal signal this means all poles of $X_L(s)$ lie strictly in the left half-plane — the same condition as BIBO stability of an impulse response. If the ROC does not contain the $j\omega$-axis, the Fourier transform either does not exist or exists only as a generalized function (Section 16.2).

### Existence

The integral defining $X(\omega)$ is improper on both ends, so we need conditions under which it makes sense.

> [!theorem] Sufficient conditions for existence (Dirichlet conditions)
> Suppose
> 1. $x(t)$ is **absolutely integrable**: $\int_{-\infty}^{\infty}\lvert x(t)\rvert\,dt<\infty$;
> 2. in every finite interval $x(t)$ has finitely many maxima and minima and finitely many finite discontinuities.
>
> Then $X(\omega)$ exists for every $\omega$, is continuous and bounded, $\lvert X(\omega)\rvert\le\int\lvert x\rvert\,dt$, and $X(\omega)\to0$ as $\lvert\omega\rvert\to\infty$ (the **Riemann–Lebesgue lemma**).

*Proof sketch.* Since $\lvert e^{-j\omega t}\rvert=1$, $\lvert X(\omega)\rvert\le\int\lvert x(t)\rvert dt<\infty$, which is the bound. Continuity follows because $\lvert X(\omega+h)-X(\omega)\rvert\le\int\lvert x(t)\rvert\,\lvert e^{-jht}-1\rvert dt$, and the integrand tends to zero for each $t$ while being dominated by $2\lvert x(t)\rvert$. For Riemann–Lebesgue, approximate $x$ by a step function; for a single rectangle $\int_a^b e^{-j\omega t}dt=(e^{-j\omega a}-e^{-j\omega b})/(j\omega)\to0$, and the approximation error can be made arbitrarily small uniformly in $\omega$. $\square$

There is a second, equally important class.

> [!theorem] Finite-energy signals (Plancherel)
> If $x(t)$ has finite energy, $\int_{-\infty}^{\infty}\lvert x(t)\rvert^2dt<\infty$, then the truncated transforms $X_T(\omega)=\int_{-T}^{T}x(t)e^{-j\omega t}dt$ converge **in mean square** to a function $X(\omega)$ as $T\to\infty$, the inverse formula holds in the same mean-square sense, and Parseval's theorem (Section 16.7) holds.

Absolute integrability and finite energy are different conditions. The sinc pulse $\sin(t)/t$ has finite energy but is not absolutely integrable (its envelope decays only as $1/\lvert t\rvert$); its transform exists in the mean-square sense and turns out to be a rectangle — discontinuous, which is consistent with the Dirichlet theorem *not* applying. Conversely $1/\sqrt{\lvert t\rvert}$ on $\lvert t\rvert<1$ is absolutely integrable but has infinite energy.

Signals such as $1$, $u(t)$, $\cos\omega_0 t$ and $\operatorname{sgn}(t)$ are neither absolutely integrable nor of finite energy. They are nevertheless indispensable — a DC supply, a switch closing, a carrier wave — and they have **generalized transforms** containing Dirac impulses $\delta(\omega)$, developed in Section 16.2.

### Standard functions: rect, sinc, triangle, signum

We use the following standard shapes throughout the chapter.

> [!definition] Rectangle, sinc, triangle and signum
>
> $$
> \operatorname{rect}(t) = \begin{cases}1, & \lvert t\rvert<\tfrac12\\ \tfrac12, & \lvert t\rvert=\tfrac12\\ 0, & \lvert t\rvert>\tfrac12\end{cases},
> \qquad
> \operatorname{sinc}(x) = \frac{\sin(\pi x)}{\pi x},\ \ \operatorname{sinc}(0)=1,
> $$
>
> $$
> \Lambda(t)=\max\big(0,\,1-\lvert t\rvert\big),\qquad
> \operatorname{sgn}(t)=\begin{cases}1,&t>0\\0,&t=0\\-1,&t<0\end{cases}.
> $$

So $\operatorname{rect}(t/T)$ is a unit-height pulse of width $T$ centred on $t=0$, and $\Lambda(t/T)$ is a unit-height triangle of base $2T$. This is the **normalized** sinc: its zeros are at the nonzero integers and $\int_{-\infty}^{\infty}\operatorname{sinc}(x)\,dx=1$.

> [!warning] Common Pitfall — two sinc conventions
> Mathematics texts and some circuits books define $\operatorname{sinc}(x)=\sin(x)/x$ (unnormalized, zeros at $x=k\pi$). This book uses $\sin(\pi x)/(\pi x)$, the convention of numerical libraries such as NumPy and MATLAB. Whenever a formula containing "sinc" is copied from another source, check which one it uses; the safest habit is to write the explicit ratio $2\sin(\omega T/2)/\omega$ when there is any doubt.

### Worked examples

> [!example] Worked Example 16.1 — the one-sided exponential
> Find the transform of $x(t)=e^{-at}u(t)$, $a>0$, and its magnitude and phase.
>
> **Step 1 — integrate.** Only $t>0$ contributes:
>
> $$
> X(\omega)=\int_0^{\infty} e^{-at}e^{-j\omega t}\,dt=\int_0^\infty e^{-(a+j\omega)t}dt=\left[\frac{-e^{-(a+j\omega)t}}{a+j\omega}\right]_0^{\infty}=\frac{1}{a+j\omega}.
> $$
>
> The upper limit vanishes because $\lvert e^{-(a+j\omega)t}\rvert=e^{-at}\to0$; this is where $a>0$ is used.
>
> **Step 2 — polar form.** Multiply numerator and denominator by the conjugate, or read off directly:
>
> $$
> \lvert X(\omega)\rvert=\frac{1}{\sqrt{a^2+\omega^2}},\qquad \angle X(\omega)=-\arctan\frac{\omega}{a}.
> $$
>
> **Step 3 — interpret.** $X(0)=1/a$, which equals $\int_0^\infty e^{-at}dt$ ✓. At $\omega=a$ the magnitude is $1/(a\sqrt2)$, i.e. $1/\sqrt2$ of its DC value (−3 dB), and the phase is $-45^\circ$. For $\omega\gg a$, $\lvert X\rvert\approx1/\lvert\omega\rvert$ and the phase approaches $-90^\circ$. The result equals the Laplace transform $1/(s+a)$ at $s=j\omega$, which is legitimate because the ROC $\operatorname{Re}\{s\}>-a$ contains the $j\omega$-axis.
>
> **Numerical check** ($a=2$, $\omega=2$): $X=1/(2+2j)=0.25-0.25j$, $\lvert X\rvert=0.3536=1/(2\sqrt2)$ ✓, $\angle X=-45^\circ$ ✓.

> [!example] Worked Example 16.2 — the rectangular pulse
> Find the transform of $x(t)=A\operatorname{rect}(t/T)$, a pulse of height $A$ and width $T$.
>
> $$
> X(\omega)=A\int_{-T/2}^{T/2} e^{-j\omega t}\,dt=A\,\frac{e^{j\omega T/2}-e^{-j\omega T/2}}{j\omega}=\frac{2A\sin(\omega T/2)}{\omega}=AT\operatorname{sinc}\!\left(\frac{\omega T}{2\pi}\right).
> $$
>
> The last step writes $2\sin(\omega T/2)/\omega=T\sin(\pi u)/(\pi u)$ with $u=\omega T/2\pi$.
>
> - $X(0)=AT$, the pulse area ✓.
> - Zeros occur where $\omega T/2=k\pi$, i.e. $\omega=2\pi k/T$, $k=\pm1,\pm2,\dots$ (in Hz: $f=k/T$).
> - The main lobe runs from $-2\pi/T$ to $2\pi/T$: width $4\pi/T$ rad/s. **Halving the pulse width doubles the spectral width.**
> - $X$ is real because $x$ is real and even; its sign alternates between lobes, so the phase alternates between $0$ and $\pm\pi$.
> - The envelope decays as $2A/\lvert\omega\rvert$ — slowly, because $x$ has jump discontinuities.
>
> For a 1 μs pulse ($T=10^{-6}$ s) the first null is at $f=1$ MHz: a rule of thumb used everywhere in pulse and digital-signal work is that a pulse of duration $T$ needs a bandwidth of order $1/T$ Hz.

```desmos-graph
left=-12; right=12; bottom=-1.5; top=4.5
---
y=2\frac{\sin(x)}{x}|blue
y=4\frac{\sin(2x)}{2x}|red
y=0|black|dotted
```

*Spectra $T\operatorname{sinc}(\omega T/2\pi)=2\sin(\omega T/2)/\omega$ of unit-height rectangular pulses of width $T=2$ (blue) and $T=4$ (red): the longer pulse has twice the peak and half the main-lobe width.*

> [!example] Worked Example 16.3 — the two-sided exponential
> Find the transform of $x(t)=e^{-a\lvert t\rvert}$, $a>0$.
>
> Split the integral at $t=0$:
>
> $$
> X(\omega)=\int_{-\infty}^{0}e^{(a-j\omega)t}dt+\int_0^\infty e^{-(a+j\omega)t}dt=\frac{1}{a-j\omega}+\frac{1}{a+j\omega}=\frac{2a}{a^2+\omega^2}.
> $$
>
> The result is real and even, as it must be for a real even signal, and positive for all $\omega$. Its value at the origin is $X(0)=2/a=\int e^{-a\lvert t\rvert}dt$ ✓. The decay is $1/\omega^2$, one power faster than for the one-sided exponential: $e^{-a\lvert t\rvert}$ is continuous (only its derivative jumps), whereas $e^{-at}u(t)$ itself jumps at $t=0$. This shape is called a **Lorentzian**; it appears again as the power spectrum of a first-order Gauss–Markov process in [Ch. 21](../part-6-probability-statistics/21-random-processes.md).

> [!example] Worked Example 16.4 — the Gaussian pulse
> Find the transform of $x(t)=e^{-t^2/(2\sigma^2)}$, $\sigma>0$.
>
> **Step 1 — complete the square** in the exponent:
>
> $$
> -\frac{t^2}{2\sigma^2}-j\omega t=-\frac{1}{2\sigma^2}\left(t+j\sigma^2\omega\right)^2-\frac{\sigma^2\omega^2}{2}.
> $$
>
> (Expand the square: $-(t^2+2j\sigma^2\omega t-\sigma^4\omega^2)/(2\sigma^2)=-t^2/(2\sigma^2)-j\omega t+\sigma^2\omega^2/2$; subtracting $\sigma^2\omega^2/2$ restores the original.)
>
> **Step 2 — shift the contour.** Then
>
> $$
> X(\omega)=e^{-\sigma^2\omega^2/2}\int_{-\infty}^{\infty}e^{-(t+j\sigma^2\omega)^2/(2\sigma^2)}\,dt .
> $$
>
> The remaining integral is a Gaussian integral along a line parallel to the real axis. Because $e^{-z^2/(2\sigma^2)}$ is analytic and decays rapidly as $\lvert\operatorname{Re}z\rvert\to\infty$ in any horizontal strip, Cauchy's theorem ([Ch. 18](18-complex-analysis.md)) lets us move the line back to the real axis without changing the value, giving $\int e^{-t^2/(2\sigma^2)}dt=\sigma\sqrt{2\pi}$.
>
> $$
> e^{-t^2/(2\sigma^2)}\ \leftrightarrow\ \sigma\sqrt{2\pi}\,e^{-\sigma^2\omega^2/2}.
> $$
>
> **Interpretation.** The transform of a Gaussian is a Gaussian. Its width in $\omega$ is $1/\sigma$: a pulse of time-width $\sigma$ has frequency-width $1/\sigma$. With $\sigma=1$ the pair is self-reciprocal up to the factor $\sqrt{2\pi}$. Because the Gaussian is infinitely smooth, its spectrum decays faster than any power of $\omega$.
>
> **Numerical check** ($\sigma=0.5$, $\omega=2$): formula gives $0.5\sqrt{2\pi}e^{-0.5}=0.7602$; numerical integration of $\int e^{-2t^2}\cos(2t)\,dt$ (the sine part vanishes by oddness) also gives $0.7602$ ✓.

> [!tip] Engineering Note — what the transform tells a circuit designer
> The magnitude $\lvert X(\omega)\rvert$ shows which frequency bands the signal occupies, and therefore what bandwidth an amplifier, cable or filter must pass to reproduce it. The rectangle's $1/\omega$ tail is why sharp-edged digital signals radiate interference far above their clock rate; rounding the edges (Gaussian or raised-cosine shaping) converts that tail into a much faster decay. The phase $\angle X(\omega)$ records *when* each component occurs: a delay shows up purely in the phase (Section 16.3).

## 16.2 Inverse transform

### Statement

The derivation in Section 16.1 already suggested how to recover $x(t)$ from $X(\omega)$. The precise statement is the **Fourier integral theorem**.

> [!theorem] Inverse Fourier transform (Fourier integral theorem)
> If $x(t)$ satisfies the Dirichlet conditions of Section 16.1, then
>
> $$
> x(t)=\mathcal{F}^{-1}\{X(\omega)\}=\frac{1}{2\pi}\int_{-\infty}^{\infty}X(\omega)\,e^{j\omega t}\,d\omega
> $$
>
> at every point of continuity of $x$, where the integral is understood as the limit $\lim_{W\to\infty}\frac{1}{2\pi}\int_{-W}^{W}$. At a jump discontinuity the integral converges to the midpoint $\tfrac12[x(t^-)+x(t^+)]$.

The inverse has exactly the same form as the forward transform except for the sign of the exponent and the factor $1/(2\pi)$. The factor appears because we integrate over angular frequency: $d\omega=2\pi\,df$. In words, the inverse formula says that $x(t)$ is a superposition of complex exponentials $e^{j\omega t}$, the one at frequency $\omega$ having complex weight $X(\omega)\,d\omega/2\pi$.

*Proof sketch.* Insert the definition of $X$ into the truncated inverse and exchange the order of integration (justified by absolute integrability of $x$ and the finite range in $\omega$):

$$
x_W(t)=\frac{1}{2\pi}\int_{-W}^{W}\left[\int_{-\infty}^{\infty}x(\tau)e^{-j\omega\tau}d\tau\right]e^{j\omega t}d\omega
=\int_{-\infty}^{\infty}x(\tau)\,D_W(t-\tau)\,d\tau,\qquad D_W(s)=\frac{\sin Ws}{\pi s}.
$$

The kernel $D_W$ is the continuous-time analogue of the Dirichlet kernel of Fourier series. For every $W$ it has total area $\int D_W(s)\,ds=1$, and as $W\to\infty$ it oscillates ever faster away from $s=0$ while its central peak (height $W/\pi$, width about $2\pi/W$) grows. By the Riemann–Lebesgue lemma the contribution away from $s=0$ tends to zero, so $D_W$ acts in the limit like $\delta(s)$ and $x_W(t)\to x(t)$ at points of continuity. At a jump, half of the kernel's area lies on each side of $s=0$ (it is even), so the limit is the average of the one-sided limits. $\square$

The truncation $x_W$ is exactly what an ideal low-pass filter of bandwidth $W$ does to $x$, and near a jump it overshoots by about 9 % of the jump no matter how large $W$ is — the **Gibbs phenomenon**, inherited unchanged from Fourier series.

> [!tip] The hertz convention
> Many communications and optics texts (for example Bracewell) use frequency $f$ in Hz instead of $\omega$:
>
> $$
> X_f(f)=\int_{-\infty}^{\infty}x(t)\,e^{-j2\pi ft}\,dt,\qquad x(t)=\int_{-\infty}^{\infty}X_f(f)\,e^{j2\pi ft}\,df .
> $$
>
> The two are related by $X_f(f)=X(2\pi f)$. The Hz form is symmetric (no $2\pi$ outside either integral), which removes $2\pi$ factors from duality, the multiplication theorem and Parseval's theorem; the $\omega$ form matches the Laplace variable $s=\sigma+j\omega$ and circuit phasors. **This book uses $\omega$ throughout.** To convert a result from a Hz-convention source, substitute $f=\omega/2\pi$ in the transform and, in any integral over frequency, $df=d\omega/2\pi$. (A third, "unitary" convention puts $1/\sqrt{2\pi}$ in front of both integrals; it is common in physics and is not used here.)

### Methods of inversion

1. **Table lookup plus properties.** This is the usual method. A short table of pairs (below) combined with the properties of Section 16.3 covers almost every signal met in practice.
2. **Partial fractions** for rational functions of $j\omega$, exactly as for Laplace transforms. Treat $j\omega$ as the variable, expand, and use

   $$
   \frac{1}{(a+j\omega)^m}\ \leftrightarrow\ \frac{t^{m-1}}{(m-1)!}\,e^{-at}u(t),\qquad \operatorname{Re}\{a\}>0 .
   $$

   If a term has $\operatorname{Re}\{a\}<0$, the corresponding absolutely integrable time function is anticausal: $1/(a+j\omega)\leftrightarrow-e^{-at}u(-t)$ for $\operatorname{Re}\{a\}<0$. A Fourier transform, unlike a unilateral Laplace transform, needs no separately stated ROC because the $j\omega$-axis *is* the ROC; the sign of $\operatorname{Re}\{a\}$ decides causal versus anticausal.
3. **Direct integration**, often by contour integration and residues ([Ch. 18](18-complex-analysis.md)), which is how rational spectra with mixed left- and right-half-plane poles are inverted systematically.

### Worked examples

> [!example] Worked Example 16.5 — ideal low-pass spectrum
> Let $X(\omega)=\operatorname{rect}(\omega/2W)$, i.e. $X=1$ for $\lvert\omega\rvert<W$ and $0$ otherwise. Find $x(t)$.
>
> $$
> x(t)=\frac{1}{2\pi}\int_{-W}^{W}e^{j\omega t}\,d\omega=\frac{1}{2\pi}\cdot\frac{e^{jWt}-e^{-jWt}}{jt}=\frac{\sin(Wt)}{\pi t}=\frac{W}{\pi}\operatorname{sinc}\!\left(\frac{Wt}{\pi}\right).
> $$
>
> Checks: $x(0)=W/\pi$, which equals $\frac{1}{2\pi}\int X\,d\omega=2W/2\pi$ ✓. The zeros are at $t=k\pi/W$, $k\ne0$.
>
> If $X$ is the frequency response of a filter, $x$ is its impulse response. The ideal ("brick-wall") low-pass filter therefore has a sinc impulse response that extends over all $t<0$: it responds before the impulse arrives. It is **non-causal** and hence not physically realisable — a fact read directly from the inverse transform. Realisable filters approximate the brick wall with a delayed, truncated sinc or with rational responses.

> [!example] Worked Example 16.6 — partial fractions with distinct poles
> Invert $X(\omega)=\dfrac{j\omega+3}{(j\omega+1)(j\omega+2)}$.
>
> Write $X=\dfrac{A}{1+j\omega}+\dfrac{B}{2+j\omega}$. By cover-up,
>
> $$
> A=\left.\frac{j\omega+3}{j\omega+2}\right|_{j\omega=-1}=2,\qquad B=\left.\frac{j\omega+3}{j\omega+1}\right|_{j\omega=-2}=-1 .
> $$
>
> Both poles are in the left half-plane ($a=1,2>0$), so both terms are causal:
>
> $$
> x(t)=\left(2e^{-t}-e^{-2t}\right)u(t).
> $$
>
> Check at $\omega=0$: $X(0)=3/2$ and $\int_0^\infty(2e^{-t}-e^{-2t})dt=2-\tfrac12=\tfrac32$ ✓. Check at $t=0^+$: $x(0^+)=1$, and the initial-value rule $x(0^+)=\lim_{\omega\to\infty}j\omega X(\omega)=1$ ✓.

> [!example] Worked Example 16.7 — complex poles, and a two-sided signal
> *(a)* Invert $X(\omega)=\dfrac{1}{(j\omega)^2+2j\omega+5}$.
>
> Complete the square: $(j\omega)^2+2j\omega+5=(j\omega+1)^2+4$. From the pair $e^{-at}\sin(bt)u(t)\leftrightarrow b/[(a+j\omega)^2+b^2]$ (the Laplace pair at $s=j\omega$, valid since $a=1>0$),
>
> $$
> x(t)=\tfrac12e^{-t}\sin(2t)\,u(t).
> $$
>
> Check: $X(0)=1/5$, and $\int_0^\infty\tfrac12e^{-t}\sin2t\,dt=\tfrac12\cdot\frac{2}{1+4}=\tfrac15$ ✓.
>
> *(b)* Invert $Y(\omega)=\dfrac{6}{(2+j\omega)(3-j\omega)}$.
>
> Partial fractions in $j\omega$: $Y=\dfrac{A}{2+j\omega}+\dfrac{B}{3-j\omega}$ with $A=6/(3+2)=6/5$ (set $j\omega=-2$) and $B=6/(2+3)=6/5$ (set $j\omega=3$). The first term is causal. The second, $1/(3-j\omega)$, is the transform of $e^{3t}u(-t)$ (replace $t\to-t$ in $e^{-3t}u(t)\leftrightarrow1/(3+j\omega)$, which replaces $\omega\to-\omega$). Hence
>
> $$
> y(t)=\tfrac65\left[e^{-2t}u(t)+e^{3t}u(-t)\right],
> $$
>
> a two-sided signal that decays in both directions. Check: $Y(0)=1$ and $\int y\,dt=\tfrac65(\tfrac12+\tfrac13)=1$ ✓. A unilateral Laplace table cannot produce this answer; the Fourier transform handles it naturally because its "ROC" is fixed to be the $j\omega$-axis.

### Generalized transforms

Many important signals — constants, steps, sinusoids — are not absolutely integrable. They acquire transforms when we allow the Dirac impulse in the frequency domain. The cleanest route is to compute transforms of impulses directly and then use the inverse formula.

**Impulse.** By the sifting property, $\mathcal{F}\{\delta(t-t_0)\}=\int\delta(t-t_0)e^{-j\omega t}dt=e^{-j\omega t_0}$. In particular $\delta(t)\leftrightarrow1$: an impulse contains all frequencies with equal weight, which is why an impulse response characterises an LTI system completely.

**Complex exponential and constant.** Apply the inverse formula to $X(\omega)=2\pi\delta(\omega-\omega_0)$:

$$
\frac{1}{2\pi}\int_{-\infty}^{\infty} 2\pi\delta(\omega-\omega_0)\,e^{j\omega t}\,d\omega=e^{j\omega_0t}.
$$

Hence

$$
e^{j\omega_0t}\leftrightarrow2\pi\delta(\omega-\omega_0),\qquad 1\leftrightarrow2\pi\delta(\omega),
$$

and by Euler's formula

$$
\cos\omega_0t\leftrightarrow\pi[\delta(\omega-\omega_0)+\delta(\omega+\omega_0)],\qquad
\sin\omega_0t\leftrightarrow-j\pi[\delta(\omega-\omega_0)-\delta(\omega+\omega_0)].
$$

A pure tone has all its spectral content concentrated at $\pm\omega_0$, with impulse areas $\pi$ (not $\tfrac12$: the areas are $2\pi$ times the Fourier-series coefficients $\tfrac12$).

**A limiting argument.** The same pair can be obtained without the inverse formula by viewing $1$ as the limit of $e^{-a\lvert t\rvert}$ as $a\to0^+$. Its transform $2a/(a^2+\omega^2)$ has area $\int2a/(a^2+\omega^2)\,d\omega=2\pi$ for every $a$, and as $a\to0$ it becomes an ever taller, narrower peak at $\omega=0$. Its limit is therefore $2\pi\delta(\omega)$ ✓.

```desmos-graph
left=-6; right=6; bottom=-0.5; top=8.5
---
y=\frac{2}{1+x^2}|blue
y=\frac{1}{0.25+x^2}|red
y=\frac{0.5}{0.0625+x^2}|green
```

*The transforms $2a/(a^2+\omega^2)$ of $e^{-a\lvert t\rvert}$ for $a=1$ (blue), $0.5$ (red) and $0.25$ (green): every curve has area $2\pi$, and as $a\to0$ they converge to the impulse $2\pi\delta(\omega)$, the transform of the constant 1.*

**Signum.** Write $\operatorname{sgn}(t)=\lim_{a\to0^+}\left[e^{-at}u(t)-e^{at}u(-t)\right]$. The bracket has transform

$$
\frac{1}{a+j\omega}-\frac{1}{a-j\omega}=\frac{-2j\omega}{a^2+\omega^2}\ \xrightarrow{a\to0}\ \frac{-2j\omega}{\omega^2}=\frac{2}{j\omega}.
$$

There is no impulse at $\omega=0$ because the bracket is odd, so its transform is odd and integrates to zero across the origin. Thus $\operatorname{sgn}(t)\leftrightarrow 2/(j\omega)$.

**Unit step.** Since $u(t)=\tfrac12[1+\operatorname{sgn}(t)]$,

$$
u(t)\leftrightarrow\pi\delta(\omega)+\frac{1}{j\omega}.
$$

The impulse $\pi\delta(\omega)$ carries the average (DC) value $\tfrac12$ of the step; the $1/(j\omega)$ term carries the jump.

> [!warning] Common Pitfall — substituting $s=j\omega$ outside the ROC
> The Fourier transform of $u(t)$ is **not** obtained by putting $s=j\omega$ in $1/s$: the ROC of $\mathcal{L}\{u(t)\}$ is $\operatorname{Re}\{s\}>0$, which excludes the $j\omega$-axis. The extra $\pi\delta(\omega)$ accounts for the DC value $\tfrac12$ of $u(t)$. Similarly $\mathcal{F}\{\cos(\omega_0t)u(t)\}=\dfrac{j\omega}{\omega_0^2-\omega^2}+\dfrac{\pi}{2}[\delta(\omega-\omega_0)+\delta(\omega+\omega_0)]$, not just $s/(s^2+\omega_0^2)$ at $s=j\omega$. Substitute $s=j\omega$ only when the ROC strictly contains the imaginary axis. When simple poles lie *on* the axis, each pole at $s=j\omega_p$ with residue $r$ contributes an extra $\pi r\,\delta(\omega-\omega_p)$.

> [!example] Worked Example 16.8 — a DC-plus-tone signal and a switched constant
> *(a)* $x(t)=3+2\cos(5t)-\sin(5t)$. Term by term:
>
> $$
> X(\omega)=6\pi\delta(\omega)+2\pi[\delta(\omega-5)+\delta(\omega+5)]+j\pi[\delta(\omega-5)-\delta(\omega+5)]
> =6\pi\delta(\omega)+(2+j)\pi\,\delta(\omega-5)+(2-j)\pi\,\delta(\omega+5).
> $$
>
> The impulse weights at $\pm5$ are complex conjugates, as they must be for a real signal. Dividing by $2\pi$ gives the Fourier-series coefficients $c_{\pm1}=1\pm j/2$, whose magnitude $\sqrt5/2$ is half the amplitude $\sqrt{2^2+1^2}=\sqrt5$ of the combined sinusoid $\sqrt5\cos(5t+26.6^\circ)$ ✓.
>
> *(b)* $y(t)=V_0\,u(t-t_0)$, a supply switched on at $t_0$. By the time-shift property (Section 16.3),
>
> $$
> Y(\omega)=V_0e^{-j\omega t_0}\left[\pi\delta(\omega)+\frac{1}{j\omega}\right]=\pi V_0\delta(\omega)+\frac{V_0e^{-j\omega t_0}}{j\omega},
> $$
>
> where we used $e^{-j\omega t_0}\delta(\omega)=\delta(\omega)$ (the exponential equals 1 at $\omega=0$).

### Table of standard pairs

| $x(t)$ | $X(\omega)$ | Condition |
|---|---|---|
| $\delta(t)$ | $1$ | |
| $\delta(t-t_0)$ | $e^{-j\omega t_0}$ | |
| $1$ | $2\pi\delta(\omega)$ | |
| $u(t)$ | $\pi\delta(\omega)+1/(j\omega)$ | |
| $\operatorname{sgn}(t)$ | $2/(j\omega)$ | |
| $e^{-at}u(t)$ | $1/(a+j\omega)$ | $\operatorname{Re}\{a\}>0$ |
| $te^{-at}u(t)$ | $1/(a+j\omega)^2$ | $\operatorname{Re}\{a\}>0$ |
| $\dfrac{t^{m-1}}{(m-1)!}e^{-at}u(t)$ | $1/(a+j\omega)^m$ | $\operatorname{Re}\{a\}>0$ |
| $e^{-at}\sin(bt)\,u(t)$ | $b/[(a+j\omega)^2+b^2]$ | $a>0$ |
| $e^{-at}\cos(bt)\,u(t)$ | $(a+j\omega)/[(a+j\omega)^2+b^2]$ | $a>0$ |
| $e^{-a\lvert t\rvert}$ | $2a/(a^2+\omega^2)$ | $a>0$ |
| $1/(a^2+t^2)$ | $(\pi/a)\,e^{-a\lvert\omega\rvert}$ | $a>0$ |
| $e^{-t^2/(2\sigma^2)}$ | $\sigma\sqrt{2\pi}\,e^{-\sigma^2\omega^2/2}$ | $\sigma>0$ |
| $\operatorname{rect}(t/T)$ | $T\operatorname{sinc}(\omega T/2\pi)=2\sin(\omega T/2)/\omega$ | |
| $\Lambda(t/T)$ | $T\operatorname{sinc}^2(\omega T/2\pi)$ | |
| $\sin(Wt)/(\pi t)$ | $\operatorname{rect}(\omega/2W)$ | |
| $e^{j\omega_0t}$ | $2\pi\delta(\omega-\omega_0)$ | |
| $\cos\omega_0t$ | $\pi[\delta(\omega-\omega_0)+\delta(\omega+\omega_0)]$ | |
| $\sin\omega_0t$ | $-j\pi[\delta(\omega-\omega_0)-\delta(\omega+\omega_0)]$ | |
| $\sum_n\delta(t-nT)$ | $\omega_s\sum_k\delta(\omega-k\omega_s)$ | $\omega_s=2\pi/T$ |

## 16.3 Transform properties

### Why properties matter

Evaluating $\int x(t)e^{-j\omega t}dt$ from scratch for every signal is slow and error-prone. In practice one memorises a dozen pairs and a dozen properties; together they generate hundreds of transforms by pure algebra. The properties are also where the physics lives: "a delay is a linear phase", "differentiation boosts high frequencies", "compressing a pulse spreads its spectrum", "multiplying by a carrier shifts the spectrum" are each one line of this section.

Throughout, $x(t)\leftrightarrow X(\omega)$ and $y(t)\leftrightarrow Y(\omega)$.

| Property | Time domain | Frequency domain |
|---|---|---|
| Linearity | $a x(t)+b y(t)$ | $aX(\omega)+bY(\omega)$ |
| Time shift | $x(t-t_0)$ | $e^{-j\omega t_0}X(\omega)$ |
| Frequency shift | $e^{j\omega_0t}x(t)$ | $X(\omega-\omega_0)$ |
| Modulation | $x(t)\cos\omega_0t$ | $\tfrac12[X(\omega-\omega_0)+X(\omega+\omega_0)]$ |
| Scaling | $x(at)$, $a\ne0$ real | $\dfrac{1}{\lvert a\rvert}X(\omega/a)$ |
| Reversal | $x(-t)$ | $X(-\omega)$ |
| Conjugation | $x^*(t)$ | $X^*(-\omega)$ |
| Duality | $X(t)$ | $2\pi x(-\omega)$ |
| Differentiation in time | $d^nx/dt^n$ | $(j\omega)^nX(\omega)$ |
| Integration | $\int_{-\infty}^t x(\tau)d\tau$ | $X(\omega)/(j\omega)+\pi X(0)\delta(\omega)$ |
| Multiplication by $t$ | $t^n\,x(t)$ | $j^n\,d^nX/d\omega^n$ |
| Convolution | $x*y$ | $X(\omega)Y(\omega)$ |
| Multiplication | $x(t)y(t)$ | $\dfrac{1}{2\pi}(X*Y)(\omega)$ |
| Parseval | $\int\lvert x\rvert^2dt$ | $\dfrac{1}{2\pi}\int\lvert X\rvert^2d\omega$ |

The convolution, multiplication and Parseval rows are the subject of Sections 16.5 and 16.7. The others are derived now.

### Derivations

**Linearity** is inherited from linearity of the integral.

**Time shift.** Substitute $\tau=t-t_0$:

$$
\int x(t-t_0)e^{-j\omega t}dt=\int x(\tau)e^{-j\omega(\tau+t_0)}d\tau=e^{-j\omega t_0}X(\omega).
$$

A delay multiplies the transform by a unit-magnitude factor: $\lvert X\rvert$ is unchanged and the phase acquires a term $-\omega t_0$, **linear in frequency**. Each frequency component is delayed by the same time $t_0$, which requires a phase lag proportional to its frequency.

**Frequency shift.** $\int e^{j\omega_0t}x(t)e^{-j\omega t}dt=\int x(t)e^{-j(\omega-\omega_0)t}dt=X(\omega-\omega_0)$. Multiplying by a complex exponential slides the whole spectrum by $\omega_0$. The **modulation** property follows by writing $\cos\omega_0t=\tfrac12(e^{j\omega_0t}+e^{-j\omega_0t})$; similarly $x(t)\sin\omega_0t\leftrightarrow\frac{1}{2j}[X(\omega-\omega_0)-X(\omega+\omega_0)]$.

**Scaling.** For $a>0$ substitute $\tau=at$, $dt=d\tau/a$:

$$
\int x(at)e^{-j\omega t}dt=\frac1a\int x(\tau)e^{-j(\omega/a)\tau}d\tau=\frac1aX\!\left(\frac{\omega}{a}\right).
$$

For $a<0$ the substitution reverses the limits, producing a minus sign that turns $1/a$ into $1/\lvert a\rvert$. With $a=-1$ this gives the **reversal** property. Compressing in time ($\lvert a\rvert>1$) expands the spectrum horizontally and lowers it by $1/\lvert a\rvert$; the factor keeps the area $\frac{1}{2\pi}\int X\,d\omega=x(0)$ unchanged, as it must be since $x(at)$ and $x(t)$ have the same value at $t=0$.

**Conjugation.** $\int x^*(t)e^{-j\omega t}dt=\left[\int x(t)e^{j\omega t}dt\right]^*=X^*(-\omega)$. For real $x$ this gives $X(-\omega)=X^*(\omega)$, the conjugate symmetry studied in Section 16.4.

**Duality.** Write the inverse formula with $t$ replaced by $-t$ and multiply by $2\pi$:

$$
2\pi x(-t)=\int_{-\infty}^{\infty}X(\omega)e^{-j\omega t}d\omega .
$$

Now rename the integration variable $\omega\to\tau$ and the free variable $t\to\omega$: $2\pi x(-\omega)=\int X(\tau)e^{-j\omega\tau}d\tau=\mathcal{F}\{X(t)\}$. So **if $x(t)\leftrightarrow X(\omega)$ then $X(t)\leftrightarrow2\pi x(-\omega)$**. Every pair in the table therefore yields a second pair for free; the rectangle–sinc and sinc–rectangle pairs are duals.

**Differentiation.** Differentiate the inverse formula under the integral sign:

$$
\frac{dx}{dt}=\frac{d}{dt}\frac{1}{2\pi}\int X(\omega)e^{j\omega t}d\omega=\frac{1}{2\pi}\int j\omega X(\omega)e^{j\omega t}d\omega ,
$$

so $dx/dt\leftrightarrow j\omega X(\omega)$; repeating gives $(j\omega)^n$. Differentiation multiplies the spectrum by $\lvert\omega\rvert$: it **emphasises high frequencies** and destroys DC. (Equivalently, integrate $\int\dot xe^{-j\omega t}dt$ by parts, the boundary terms vanishing because $x\to0$ at $\pm\infty$.)

**Multiplication by $t$.** Differentiate the forward formula with respect to $\omega$: $dX/d\omega=\int(-jt)x(t)e^{-j\omega t}dt$, so $t\,x(t)\leftrightarrow j\,dX/d\omega$. This is the dual of the differentiation property.

**Integration.** Let $y(t)=\int_{-\infty}^tx(\tau)d\tau$. Then $y=x*u$, and the convolution theorem (Section 16.5) with $u\leftrightarrow\pi\delta(\omega)+1/(j\omega)$ gives

$$
Y(\omega)=X(\omega)\left[\pi\delta(\omega)+\frac{1}{j\omega}\right]=\frac{X(\omega)}{j\omega}+\pi X(0)\delta(\omega).
$$

The impulse term is present whenever $x$ has nonzero area, because then $y(t)\to X(0)$ as $t\to\infty$ and $y$ has a DC component $X(0)/2$ (it rises from $0$ to $X(0)$).

**Moment theorems.** Setting $\omega=0$ in the forward formula and in the $t^n$ property gives

$$
\int_{-\infty}^{\infty}t^n x(t)\,dt=j^n\left.\frac{d^nX}{d\omega^n}\right|_{\omega=0}.
$$

The area, centroid and spread of a pulse can therefore be read from the value, slope and curvature of its transform at DC.

### Worked examples

> [!example] Worked Example 16.9 — shift and duality
> Find the transforms of $x(t)=e^{-2(t-1)}u(t-1)$ and of $g(t)=\dfrac{1}{1+t^2}$.
>
> *(a)* $x$ is $e^{-2t}u(t)$ delayed by $t_0=1$. From $e^{-2t}u(t)\leftrightarrow 1/(2+j\omega)$ and the time-shift property:
>
> $$
> X(\omega)=\frac{e^{-j\omega}}{2+j\omega},\qquad \lvert X\rvert=\frac{1}{\sqrt{4+\omega^2}},\qquad \angle X=-\omega-\arctan\frac{\omega}{2}.
> $$
>
> The magnitude is that of the undelayed pulse; only the phase gained the linear term $-\omega$.
>
> *(b)* From $e^{-\lvert t\rvert}\leftrightarrow 2/(1+\omega^2)$, duality gives $\dfrac{2}{1+t^2}\leftrightarrow 2\pi e^{-\lvert-\omega\rvert}=2\pi e^{-\lvert\omega\rvert}$. Dividing by 2:
>
> $$
> \frac{1}{1+t^2}\leftrightarrow \pi e^{-\lvert\omega\rvert}.
> $$
>
> Check: $G(0)=\pi=\int_{-\infty}^{\infty}\frac{dt}{1+t^2}=[\arctan t]_{-\infty}^{\infty}$ ✓. Scaling then gives the general table entry $1/(a^2+t^2)\leftrightarrow(\pi/a)e^{-a\lvert\omega\rvert}$. Evaluating this integral directly requires residues ([Ch. 18](18-complex-analysis.md)); duality did it in one line.

> [!example] Worked Example 16.10 — combining scaling and shifting
> *(a)* Find the transform of a unit pulse that starts at $t=0$ and ends at $t=T$.
>
> The pulse is $\operatorname{rect}\big((t-T/2)/T\big)$, the centred pulse delayed by $T/2$:
>
> $$
> X(\omega)=T\operatorname{sinc}\!\left(\frac{\omega T}{2\pi}\right)e^{-j\omega T/2}=\frac{1-e^{-j\omega T}}{j\omega}.
> $$
>
> (The second form is the direct integral $\int_0^Te^{-j\omega t}dt$; factoring $e^{-j\omega T/2}$ out of $1-e^{-j\omega T}$ leaves $2j\sin(\omega T/2)$, confirming the first form.)
>
> *(b)* Given $x(t)\leftrightarrow X(\omega)$, find the transform of $y(t)=x(2t-3)$.
>
> Write $y(t)=x\big(2(t-\tfrac32)\big)$. First scale: $x(2t)\leftrightarrow\tfrac12X(\omega/2)$. Then delay *the scaled signal* by $\tfrac32$:
>
> $$
> Y(\omega)=\tfrac12X\!\left(\frac{\omega}{2}\right)e^{-j3\omega/2}.
> $$
>
> Alternatively shift first: $v(t)=x(t-3)\leftrightarrow e^{-j3\omega}X(\omega)$, and $y(t)=v(2t)\leftrightarrow\tfrac12e^{-j3\omega/2}X(\omega/2)$ — the same answer, because scaling acts on *every* $\omega$ in the expression, including the one in the exponential.
>
> **Check** with $x=\operatorname{rect}(t)$: $y=\operatorname{rect}(2t-3)$ is a pulse of width $\tfrac12$ centred at $t=\tfrac32$, so directly $Y=\tfrac12\operatorname{sinc}(\omega/4\pi)e^{-j3\omega/2}$; the formula gives $\tfrac12\operatorname{sinc}(\omega/4\pi)e^{-j3\omega/2}$ ✓.

> [!warning] Common Pitfall — the delay in $x(at-b)$
> The delay of $x(at-b)$ is $b/a$, not $b$. Writing $e^{-j\omega b}$ instead of $e^{-j\omega b/a}$ is the single most common error with these two properties. Always rewrite the argument as $a(t-b/a)$ first.

> [!example] Worked Example 16.11 — the triangle by differentiation
> Find the transform of $x(t)=\Lambda(t/T)$ without integrating.
>
> **Step 1 — differentiate until only impulses remain.** The triangle has slope $1/T$ on $(-T,0)$ and $-1/T$ on $(0,T)$, so $\dot x=\frac1T[\operatorname{rect}\text{ on }(-T,0)]-\frac1T[\operatorname{rect}\text{ on }(0,T)]$. Differentiating again, each jump of the slope becomes an impulse:
>
> $$
> \ddot x(t)=\frac1T\left[\delta(t+T)-2\delta(t)+\delta(t-T)\right].
> $$
>
> **Step 2 — transform.** By the shift property,
>
> $$
> (j\omega)^2X(\omega)=\frac1T\left[e^{j\omega T}-2+e^{-j\omega T}\right]=\frac{2\cos\omega T-2}{T}=-\frac{4\sin^2(\omega T/2)}{T}.
> $$
>
> **Step 3 — divide.** Since $(j\omega)^2=-\omega^2$,
>
> $$
> X(\omega)=\frac{4\sin^2(\omega T/2)}{T\omega^2}=T\operatorname{sinc}^2\!\left(\frac{\omega T}{2\pi}\right).
> $$
>
> Division by $\omega^2$ is safe here because $x$ has finite area, so $X$ has no impulse at the origin; the value $X(0)=T$ (by L'Hôpital or the sinc limit) equals the triangle's area $\tfrac12\cdot2T\cdot1$ ✓. The method works for any **piecewise-polynomial** signal: differentiate until the result is a sum of impulses, transform, divide by $(j\omega)^n$, and fix the DC value by the area.

> [!example] Worked Example 16.12 — multiplication by $t$
> *(a)* From $e^{-at}u(t)\leftrightarrow(a+j\omega)^{-1}$,
>
> $$
> te^{-at}u(t)\leftrightarrow j\frac{d}{d\omega}(a+j\omega)^{-1}=j\cdot(-1)(a+j\omega)^{-2}\cdot j=\frac{1}{(a+j\omega)^2},
> $$
>
> confirming the table entry. Repeating gives the $t^{m-1}/(m-1)!$ pairs used for repeated poles.
>
> *(b)* Find the transform of the "Gaussian monocycle" $x(t)=t\,e^{-t^2/2}$ (an ultra-wideband pulse shape).
>
> From $e^{-t^2/2}\leftrightarrow\sqrt{2\pi}e^{-\omega^2/2}$:
>
> $$
> X(\omega)=j\frac{d}{d\omega}\left(\sqrt{2\pi}e^{-\omega^2/2}\right)=-j\sqrt{2\pi}\,\omega\,e^{-\omega^2/2}.
> $$
>
> The result is purely imaginary and odd, as it must be for a real odd signal, and it vanishes at DC (the monocycle has zero area). Its magnitude peaks where $\frac{d}{d\omega}(\omega e^{-\omega^2/2})=0$, i.e. at $\omega=1$ rad/s, where $\lvert X\rvert=\sqrt{2\pi}e^{-1/2}=1.520$. Equivalently, $x=-\frac{d}{dt}e^{-t^2/2}$, and the differentiation property gives $-j\omega\sqrt{2\pi}e^{-\omega^2/2}$ — the same answer by a second route ✓.

> [!example] Worked Example 16.13 — an RF pulse
> A radar transmitter emits $x(t)=A\operatorname{rect}(t/T)\cos(\omega_ct)$ with $\omega_c\gg2\pi/T$. Find and sketch its spectrum.
>
> By the modulation property applied to $A\operatorname{rect}(t/T)\leftrightarrow AT\operatorname{sinc}(\omega T/2\pi)$:
>
> $$
> X(\omega)=\frac{AT}{2}\left[\operatorname{sinc}\!\left(\frac{(\omega-\omega_c)T}{2\pi}\right)+\operatorname{sinc}\!\left(\frac{(\omega+\omega_c)T}{2\pi}\right)\right].
> $$
>
> The spectrum consists of two copies of the pulse spectrum, each of half the height, centred on $\pm\omega_c$. The main lobe around $\omega_c$ has null-to-null width $4\pi/T$ rad/s, or $2/T$ Hz. For $A=1$ V, $T=1$ μs, $f_c=10$ MHz: peaks of $0.5$ μV·s at $\pm10$ MHz, first nulls at $9$ and $11$ MHz. The two copies interfere slightly with each other through their sidelobes; the condition $\omega_c\gg2\pi/T$ makes this negligible.

```desmos-graph
left=-16; right=16; bottom=-0.5; top=1.2
---
y=\frac{\sin(x-8)}{x-8}+\frac{\sin(x+8)}{x+8}|blue
y=\frac{\sin(x)}{x}|grey|dashed
```

*Modulation shifts the spectrum: the baseband pulse spectrum (grey dashed, $T=2$, peak normalised to 1) is copied to $\pm\omega_c=\pm8$ rad/s with half the height (blue).*

> [!tip] Engineering Note — differentiation turns ODEs into algebra
> For a circuit described by $a_2\ddot y+a_1\dot y+a_0y=b_1\dot x+b_0x$, the differentiation property gives $[a_2(j\omega)^2+a_1j\omega+a_0]Y=(b_1j\omega+b_0)X$, so
>
> $$
> Y(\omega)=H(\omega)X(\omega),\qquad H(\omega)=\frac{b_0+jb_1\omega}{a_0-a_2\omega^2+ja_1\omega}.
> $$
>
> This is the phasor method of AC circuit analysis, extended from single sinusoids to arbitrary spectra: an inductor's impedance $j\omega L$ and a capacitor's $1/(j\omega C)$ are the differentiation and integration properties applied to $v=L\,di/dt$ and $v=\frac1C\int i\,dt$. The Fourier method gives the zero-state response of a stable system; initial conditions and unstable systems are the territory of the Laplace transform ([Ch. 14](14-laplace-transforms.md)).

## 16.4 Frequency-domain representation

A signal and its transform are two complete descriptions of the same object. This section is about reading the frequency-domain description: what the magnitude and phase spectra mean, what symmetry a real signal forces on them, how an LTI system acts on them, and how "width" in one domain constrains width in the other.

### Magnitude and phase spectra

The pair $\lvert X(\omega)\rvert$ (the **magnitude** or **amplitude spectrum**) and $\angle X(\omega)$ (the **phase spectrum**) is the frequency-domain representation of $x(t)$. Equivalently one may use the real and imaginary parts $X_R(\omega)=\operatorname{Re}\{X\}$, $X_I(\omega)=\operatorname{Im}\{X\}$.

- Magnitudes spanning many decades are plotted in decibels, $20\log_{10}\lvert X(\omega)\rvert$ (relative to a stated reference), against $\log\omega$ — the Bode format. A factor of 2 in magnitude is about 6 dB; a factor of 10 is 20 dB.
- Phase is defined only modulo $2\pi$. Numerical routines return the **principal value** in $(-\pi,\pi]$, producing artificial jumps of $2\pi$; **unwrapping** removes them by adding multiples of $2\pi$ to make the phase continuous wherever $\lvert X\rvert\ne0$. Where $X$ passes through zero and changes sign (the nulls of a sinc), the phase genuinely jumps by $\pi$.
- Phase matters as much as magnitude. Two signals with the same $\lvert X\rvert$ but different phases can look entirely different in time: an impulse and a "chirp" can share a flat magnitude spectrum.

### Symmetry of real signals

> [!theorem] Conjugate symmetry
> If $x(t)$ is real then
>
> $$
> X(-\omega)=X^*(\omega),
> $$
>
> so $\lvert X(\omega)\rvert$ and $\operatorname{Re}\{X(\omega)\}$ are **even**, while $\angle X(\omega)$ and $\operatorname{Im}\{X(\omega)\}$ are **odd**. Moreover:
> - $x$ real and even $\Rightarrow$ $X$ real and even;
> - $x$ real and odd $\Rightarrow$ $X$ purely imaginary and odd;
> - for any real $x$ with even part $x_e(t)=\tfrac12[x(t)+x(-t)]$ and odd part $x_o(t)=\tfrac12[x(t)-x(-t)]$: $\ x_e\leftrightarrow\operatorname{Re}\{X\}$ and $x_o\leftrightarrow j\operatorname{Im}\{X\}$.

*Proof.* For real $x$, $X^*(\omega)=\int x(t)e^{+j\omega t}dt=X(-\omega)$. Writing $e^{-j\omega t}=\cos\omega t-j\sin\omega t$,

$$
X(\omega)=\int x(t)\cos\omega t\,dt-j\int x(t)\sin\omega t\,dt .
$$

The cosine integral sees only the even part of $x$ (the odd part integrates to zero against an even function), and the sine integral sees only the odd part. This gives the even/odd statements; the cosine integral is even in $\omega$, the sine integral odd. $\square$

Consequently, for real signals the negative-frequency half of the spectrum carries no new information. Engineers often plot only $\omega\ge0$ (a **one-sided** spectrum), remembering to double energies when integrating over $\omega>0$ only (Section 16.7). Negative frequencies are not a physical mystery: $\cos\omega_0t$ is the sum of two counter-rotating phasors $e^{\pm j\omega_0t}/2$, and the transform simply accounts for both.

> [!example] Worked Example 16.14 — even and odd parts of the one-sided exponential
> For $x(t)=e^{-at}u(t)$, $X=\dfrac{1}{a+j\omega}=\dfrac{a}{a^2+\omega^2}-j\dfrac{\omega}{a^2+\omega^2}$.
>
> The even part of $x$ is $x_e=\tfrac12e^{-a\lvert t\rvert}$, whose transform is $\tfrac12\cdot\dfrac{2a}{a^2+\omega^2}=\dfrac{a}{a^2+\omega^2}=\operatorname{Re}\{X\}$ ✓.
> The odd part is $x_o=\tfrac12\operatorname{sgn}(t)e^{-a\lvert t\rvert}$, which must therefore have transform $j\operatorname{Im}\{X\}=-j\omega/(a^2+\omega^2)$. As $a\to0$ this tends to $-j/\omega=1/(j\omega)$, which is indeed the transform of $\tfrac12\operatorname{sgn}(t)$ ✓.
>
> This example shows a deeper fact: for a **causal** signal, $x_o(t)=\operatorname{sgn}(t)\,x_e(t)$, so the odd part is determined by the even part and hence $\operatorname{Im}\{X\}$ is determined by $\operatorname{Re}\{X\}$. (The precise relation is a Hilbert transform; for minimum-phase systems it also links the magnitude and phase of a frequency response.)

### Frequency response of LTI systems

Complex exponentials are **eigenfunctions** of linear time-invariant (LTI) systems. If the impulse response is $h(t)$, the input $e^{j\omega t}$ produces

$$
y(t)=\int_{-\infty}^{\infty}h(\tau)e^{j\omega(t-\tau)}d\tau = \left[\int h(\tau)e^{-j\omega\tau}d\tau\right]e^{j\omega t}=H(\omega)\,e^{j\omega t},\qquad H(\omega)=\mathcal{F}\{h(t)\}.
$$

The output is the input multiplied by a complex constant; the constant (the eigenvalue) is the **frequency response** $H(\omega)$. For a real system ($h$ real), the input $\cos(\omega_0t+\theta)$ — the real part of $e^{j\theta}e^{j\omega_0t}$ — produces

$$
y(t)=\lvert H(\omega_0)\rvert\cos\big(\omega_0t+\theta+\angle H(\omega_0)\big).
$$

The **magnitude response** $\lvert H\rvert$ is the gain and the **phase response** $\angle H$ the phase shift at each frequency. For a general input, superposition over all frequencies gives $Y(\omega)=H(\omega)X(\omega)$, proved formally as the convolution theorem in Section 16.5.

> [!example] Worked Example 16.15 — RC low-pass filter
> For a series $R$ with shunt $C$, output across $C$: $RC\dot y+y=x$, so by the differentiation property
>
> $$
> H(\omega)=\frac{1}{1+j\omega RC},\qquad h(t)=\frac{1}{RC}e^{-t/RC}u(t).
> $$
>
> With $R=1\ \text{k}\Omega$, $C=1\ \mu\text{F}$: $RC=1$ ms, so the cutoff is $\omega_c=1/RC=1000$ rad/s ($f_c\approx159$ Hz).
>
> For the input $\cos(2000t)$: $\lvert H\rvert=1/\sqrt{1+2^2}=0.447$ ($-7.0$ dB), $\angle H=-\arctan2=-63.4^\circ$, so in steady state $y(t)=0.447\cos(2000t-63.4^\circ)$.
>
> For the input $\cos(100t)$: $\lvert H\rvert=1/\sqrt{1.01}=0.995$, $\angle H=-5.7^\circ$: frequencies well below cutoff pass almost unaltered. Far above cutoff $\lvert H\rvert\approx\omega_c/\omega$, a fall of 20 dB per decade.

```desmos-graph
left=-6; right=6; bottom=-1.8; top=1.4
---
y=\frac{1}{\sqrt{1+x^2}}|blue
y=-\arctan(x)|red
y=\frac{1}{\sqrt{2}}|grey|dashed
(1,0.7071)|label:-3 dB
```

*Magnitude $\lvert H\rvert=1/\sqrt{1+\omega^2}$ (blue, even) and phase $\angle H=-\arctan\omega$ in radians (red, odd) for $H=1/(1+j\omega)$, with $\omega$ normalised to $\omega_c$.*

> [!example] Worked Example 16.16 — a resonant second-order response
> A series RLC circuit with output across $C$ satisfies $LC\ddot y+RC\dot y+y=x$, so
>
> $$
> H(\omega)=\frac{1}{1-\omega^2LC+j\omega RC}=\frac{1}{1-(\omega/\omega_0)^2+j(\omega/\omega_0)/Q},\qquad \omega_0=\frac{1}{\sqrt{LC}},\ \ Q=\frac{\omega_0L}{R}.
> $$
>
> Take $L=10$ mH, $C=1$ μF, $R=50\ \Omega$: $\omega_0=1/\sqrt{10^{-8}}=10^4$ rad/s ($f_0=1.59$ kHz) and $Q=10^4\cdot0.01/50=2$.
>
> | $\omega$ | $\lvert H\rvert$ | $\angle H$ |
> |---|---|---|
> | $\omega_0/2=5000$ | $1.265$ | $-18.4^\circ$ |
> | $\omega_0=10^4$ | $2.000=Q$ | $-90^\circ$ |
> | $2\omega_0=2\times10^4$ | $0.316$ | $-161.6^\circ$ |
>
> (At $\omega_0/2$: denominator $=0.75+0.25j$, magnitude $0.7906$, so $\lvert H\rvert=1.265$ and $\angle H=-\arctan(1/3)$. At $2\omega_0$: denominator $=-3+j$, magnitude $\sqrt{10}$.) At resonance the real part of the denominator vanishes, the gain equals $Q$, and the phase passes through $-90^\circ$. The magnitude response peaks near $\omega_0$ and falls as $1/\omega^2$ (40 dB/decade) above it.

```desmos-graph
left=0; right=4; bottom=0; top=2.4
---
y=\frac{1}{\sqrt{(1-x^2)^2+(x/2)^2}}|blue
y=\frac{1}{\sqrt{(1-x^2)^2+(x/0.707)^2}}|red
y=\frac{1}{\sqrt{(1-x^2)^2+(x/0.5)^2}}|green
```

*Magnitude of the second-order low-pass response against $\omega/\omega_0$ for $Q=2$ (blue, resonant peak), $Q=0.707$ (red, maximally flat) and $Q=0.5$ (green, critically damped).*

### Distortionless transmission, phase delay and group delay

> [!definition] Distortionless system
> A system is **distortionless** if every input emerges scaled and delayed but otherwise unchanged: $y(t)=Kx(t-t_d)$. By the time-shift property this is equivalent to
>
> $$
> H(\omega)=Ke^{-j\omega t_d}:\qquad \lvert H(\omega)\rvert=K\ \ (\text{flat}),\qquad \angle H(\omega)=-\omega t_d\ \ (\text{linear phase}).
> $$
>
> It suffices that these hold over the band occupied by the signal.

Departures from flat magnitude are **amplitude distortion**; departures from linear phase are **phase (delay) distortion**. Two delays characterise a non-ideal phase response $\phi(\omega)=\angle H(\omega)$:

$$
\tau_p(\omega)=-\frac{\phi(\omega)}{\omega}\quad(\text{phase delay}),\qquad
\tau_g(\omega)=-\frac{d\phi}{d\omega}\quad(\text{group delay}).
$$

The phase delay is the delay of a single sinusoid at $\omega$. The group delay is the delay of the **envelope** of a narrow-band signal centred at $\omega$. To see why, let $x(t)=m(t)\cos\omega_ct$ with $m$ slowly varying, so that $X$ is concentrated near $\pm\omega_c$. Expand $\phi(\omega)\approx\phi(\omega_c)+\phi'(\omega_c)(\omega-\omega_c)$ near $\omega_c$ and take $\lvert H\rvert\approx\lvert H(\omega_c)\rvert$ there. Applying the shift and modulation properties to each term gives

$$
y(t)\approx\lvert H(\omega_c)\rvert\,m\big(t-\tau_g(\omega_c)\big)\cos\big(\omega_c(t-\tau_p(\omega_c))\big).
$$

The envelope travels at the group delay and the carrier at the phase delay. Both are equal (to $t_d$) for linear phase.

> [!example] Worked Example 16.17 — group delay of an RC filter and of an all-pass section
> *(a)* For $H=1/(1+j\omega RC)$, $\phi=-\arctan(\omega RC)$, so
>
> $$
> \tau_g(\omega)=\frac{RC}{1+\omega^2R^2C^2}.
> $$
>
> At DC the group delay is $RC$ (1 ms for the values of Worked Example 16.15); at the cutoff it has halved to $RC/2$. A pulse whose spectrum lies well inside the passband is delayed by about $RC$ with little distortion; a pulse that extends past cutoff suffers both amplitude and phase distortion.
>
> *(b)* The first-order **all-pass** section $H(\omega)=\dfrac{a-j\omega}{a+j\omega}$, $a>0$, has $\lvert H\rvert=1$ for all $\omega$ — numerator and denominator are conjugates — so there is no amplitude distortion at all. Its phase is $\phi=-2\arctan(\omega/a)$ and
>
> $$
> \tau_g(\omega)=\frac{2a}{a^2+\omega^2},
> $$
>
> which is $2/a$ at DC and falls towards zero at high frequency. With $a=10^3$ s$^{-1}$: $\tau_g(0)=2$ ms, $\tau_g(10^3)=1$ ms, $\tau_g(10^4)=19.8$ μs. A broadband pulse through this "flat" filter is dispersed, its low-frequency components arriving later than its high-frequency ones. All-pass sections are used to *equalise* group delay, i.e. to cancel the delay variation of another filter.

```desmos-graph
left=0; right=5; bottom=-3.4; top=2.4
---
y=-2\arctan(x)|red
y=\frac{2}{1+x^2}|blue
y=1|grey|dashed
```

*First-order all-pass with $a=1$: flat magnitude 1 (grey dashed), phase $-2\arctan\omega$ in radians (red) and group delay $2/(1+\omega^2)$ in seconds (blue). Flat magnitude does not imply distortionless transmission.*

### Bandwidth and the time–bandwidth trade-off

**Bandwidth measures.** For a low-pass spectrum, common definitions are:
- the **3-dB bandwidth**, where $\lvert X\rvert$ first falls to $\lvert X(0)\rvert/\sqrt2$;
- the **null-to-null** or first-null bandwidth (for a rectangular pulse of width $T$: $2\pi/T$ rad/s one-sided);
- the **essential bandwidth**, containing a stated fraction (say 95 % or 99 %) of the energy (Section 16.7);
- the **RMS bandwidth** $\Delta\omega$ defined below.

For band-pass spectra the same measures are applied to the band around the centre frequency. No single definition is "correct"; a specification must say which one it uses.

**The trade-off.** The scaling property already implies that duration and bandwidth vary inversely: replacing $x(t)$ by $x(at)$ multiplies every width in time by $1/a$ and every width in frequency by $a$, leaving their product unchanged. The product depends only on the *shape*. How small can it be? With the root-mean-square widths

$$
(\Delta t)^2=\frac{\int t^2\lvert x(t)\rvert^2dt}{\int\lvert x(t)\rvert^2dt},\qquad (\Delta\omega)^2=\frac{\int\omega^2\lvert X(\omega)\rvert^2d\omega}{\int\lvert X(\omega)\rvert^2d\omega},
$$

for a signal centred at $t=0$ and $\omega=0$:

> [!theorem] Uncertainty principle (time–bandwidth inequality)
> If $\sqrt{\lvert t\rvert}\,x(t)\to0$ as $\lvert t\rvert\to\infty$ and the widths are finite, then
>
> $$
> \Delta t\,\Delta\omega\ \ge\ \frac12 ,
> $$
>
> with equality if and only if $x(t)$ is a Gaussian, $x(t)=Ce^{-t^2/(2\sigma^2)}$.

*Proof sketch* (real $x$). Let $E=\int x^2dt$. By Parseval (Section 16.7) and $\dot x\leftrightarrow j\omega X$, $\int\omega^2\lvert X\rvert^2d\omega/2\pi=\int\dot x^2dt$, so $(\Delta\omega)^2=\int\dot x^2dt/E$. The Cauchy–Schwarz inequality gives

$$
\left(\int t\,x\,\dot x\,dt\right)^2\le\int t^2x^2dt\int\dot x^2dt=E^2(\Delta t)^2(\Delta\omega)^2 .
$$

Integrating by parts, $\int t\,x\dot x\,dt=\int t\,\tfrac12\tfrac{d}{dt}(x^2)\,dt=\big[\tfrac12tx^2\big]-\tfrac12\int x^2dt=-E/2$, the boundary term vanishing by hypothesis. Hence $E^2/4\le E^2(\Delta t\,\Delta\omega)^2$. Equality in Cauchy–Schwarz requires $\dot x=-ctx$, whose solutions are Gaussians. $\square$

> [!example] Worked Example 16.18 — time–bandwidth products
> *(a) Gaussian.* For $x=e^{-t^2/(2\sigma^2)}$, $x^2=e^{-t^2/\sigma^2}$ is proportional to a normal density of variance $\sigma^2/2$, so $\Delta t=\sigma/\sqrt2$. Similarly $\lvert X\rvert^2\propto e^{-\sigma^2\omega^2}$ has variance $1/(2\sigma^2)$, so $\Delta\omega=1/(\sigma\sqrt2)$. The product is $\tfrac12$ exactly ✓.
>
> *(b) Two-sided exponential.* For $x=e^{-\lvert t\rvert}$: $\int t^2e^{-2\lvert t\rvert}dt=\tfrac12$ and $\int e^{-2\lvert t\rvert}dt=1$, so $(\Delta t)^2=\tfrac12$. In frequency, $(\Delta\omega)^2=\int\dot x^2dt/E=1$ (because $\dot x^2=e^{-2\lvert t\rvert}$). So $\Delta t\,\Delta\omega=1/\sqrt2\approx0.707>\tfrac12$ ✓.
>
> *(c) Rectangular pulse.* $\lvert X\rvert^2\propto\sin^2(\omega T/2)/\omega^2$, so $\omega^2\lvert X\rvert^2$ does not decay and $\Delta\omega=\infty$: by the RMS measure a signal with jumps has infinite bandwidth, even though 90 % of its energy lies in the main lobe (Section 16.7). The RMS measure is sensitive to tails; the essential bandwidth is not.

```desmos-graph
left=-6; right=6; bottom=-0.2; top=2.8
---
y=e^{-x^2/2}|blue
y=\sqrt{2\pi}e^{-x^2/2}|blue|dashed
y=e^{-x^2/(2\cdot 4)}|red
y=2\sqrt{2\pi}e^{-4x^2/2}|red|dashed
```

*Gaussian pulses (solid) and their transforms (dashed): $\sigma=1$ (blue) gives a transform of the same width; $\sigma=2$ (red) is twice as long in time and half as wide in frequency.*

### Smoothness and spectral decay

The examples so far display a pattern: the rectangle (jumps) has a spectrum decaying as $1/\lvert\omega\rvert$; the triangle and two-sided exponential (continuous, with jumps in the first derivative) decay as $1/\omega^2$; the Gaussian (infinitely smooth) decays faster than any power.

> [!theorem] Smoothness and decay
> If $x,\dot x,\dots,x^{(n-1)}$ are continuous and absolutely integrable and $x^{(n)}$ is absolutely integrable (possibly with finitely many jumps), then $\lvert X(\omega)\rvert\le C/\lvert\omega\rvert^n$ for large $\lvert\omega\rvert$. If $x^{(n)}$ has jumps, $X(\omega)$ decays *exactly* as $1/\lvert\omega\rvert^{n+1}$.

*Reason.* $x^{(n)}\leftrightarrow(j\omega)^nX$, and the transform of an absolutely integrable function is bounded; if $x^{(n)}$ has a jump, its derivative contains an impulse whose transform does not decay, so $(j\omega)^{n+1}X$ stays bounded but nonzero at infinity.

> [!tip] Engineering Note — smooth edges, narrow spectra
> The spectral occupancy of a digital waveform is set by its edges, not its bit rate alone. Trapezoidal edges (continuous signal, jumps in slope) give $1/\omega^2$ (40 dB/decade) decay beyond a corner frequency of roughly $1/(\pi t_r)$ for rise time $t_r$; raised-cosine edges give $1/\omega^3$. Electromagnetic-compatibility estimates and pulse-shaping filters in data links are applications of this theorem.

## 16.5 Convolution and multiplication

### Motivation

The output of an LTI system with impulse response $h(t)$ to an input $x(t)$ is the **convolution**

$$
y(t)=(x*h)(t)=\int_{-\infty}^{\infty}x(\tau)\,h(t-\tau)\,d\tau .
$$

Each input value $x(\tau)\,d\tau$ launches a scaled, delayed copy of the impulse response, and the output is their sum. Convolution integrals are tedious to evaluate directly — the integrand changes form as $t$ moves the flipped $h$ across $x$ — and the single most useful fact about the Fourier transform is that it turns convolution into multiplication. The dual statement, that multiplication in time becomes convolution in frequency, describes modulation, windowing and sampling.

### The two theorems

> [!theorem] Convolution theorem
> If $x$ and $h$ have Fourier transforms $X$ and $H$ (for instance both absolutely integrable), and $y=x*h$, then
>
> $$
> Y(\omega)=X(\omega)H(\omega).
> $$

*Proof.* By definition,

$$
Y(\omega)=\int_{-\infty}^{\infty}\left[\int_{-\infty}^{\infty}x(\tau)h(t-\tau)\,d\tau\right]e^{-j\omega t}\,dt .
$$

Exchange the order of integration (absolutely convergent when $x,h$ are absolutely integrable) and in the inner integral put $t=\tau+\sigma$, $dt=d\sigma$:

$$
Y(\omega)=\int x(\tau)\left[\int h(\sigma)e^{-j\omega(\tau+\sigma)}d\sigma\right]d\tau=\int x(\tau)e^{-j\omega\tau}d\tau\cdot\int h(\sigma)e^{-j\omega\sigma}d\sigma=X(\omega)H(\omega).\ \square
$$

Algebraic properties of convolution follow at once from those of multiplication: $x*h=h*x$ (commutative), $(x*h)*g=x*(h*g)$ (associative — cascaded systems multiply their frequency responses), $x*(h+g)=x*h+x*g$ (parallel systems add), and $x*\delta=x$ (because $\mathcal{F}\{\delta\}=1$). Also $x(t-t_1)*h(t-t_2)=y(t-t_1-t_2)$, since delays multiply as $e^{-j\omega t_1}e^{-j\omega t_2}$.

> [!theorem] Multiplication (modulation) theorem
>
> $$
> x(t)\,y(t)\ \leftrightarrow\ \frac{1}{2\pi}\int_{-\infty}^{\infty}X(\lambda)\,Y(\omega-\lambda)\,d\lambda=\frac{1}{2\pi}(X*Y)(\omega).
> $$

*Proof.* Write $x$ via the inverse formula and compute the transform of the product:

$$
\int\left[\frac{1}{2\pi}\int X(\lambda)e^{j\lambda t}d\lambda\right]y(t)e^{-j\omega t}dt=\frac{1}{2\pi}\int X(\lambda)\left[\int y(t)e^{-j(\omega-\lambda)t}dt\right]d\lambda=\frac{1}{2\pi}\int X(\lambda)Y(\omega-\lambda)d\lambda.\ \square
$$

The $1/(2\pi)$ is the price of the angular-frequency convention; in the Hz convention it disappears. The two theorems are duals: **convolution in one domain is multiplication in the other.** The modulation property of Section 16.3 is the special case $y=\cos\omega_0t$, for which $Y$ is a pair of impulses and convolution with an impulse is a shift: $X(\omega)*\delta(\omega-\omega_0)=X(\omega-\omega_0)$.

> [!warning] Common Pitfall — the $2\pi$ in the multiplication theorem
> Forgetting the $1/(2\pi)$ is the commonest error with this theorem. A quick test: take $y(t)=1$. Then $x\cdot1$ must transform to $X$, and indeed $\frac{1}{2\pi}X*2\pi\delta=X$ ✓. Without the factor you would get $2\pi X$.

### Worked examples

> [!example] Worked Example 16.19 — convolution by transforms
> Find $y=x*h$ for $x(t)=e^{-t}u(t)$ and $h(t)=e^{-2t}u(t)$.
>
> $$
> Y(\omega)=\frac{1}{(1+j\omega)(2+j\omega)}=\frac{1}{1+j\omega}-\frac{1}{2+j\omega}
> \ \Rightarrow\ y(t)=\left(e^{-t}-e^{-2t}\right)u(t).
> $$
>
> Direct check: for $t>0$, $\int_0^t e^{-\tau}e^{-2(t-\tau)}d\tau=e^{-2t}\int_0^te^{\tau}d\tau=e^{-2t}(e^{t}-1)=e^{-t}-e^{-2t}$ ✓. Also $Y(0)=\tfrac12=X(0)H(0)=1\cdot\tfrac12$, consistent with the general rule that **areas multiply under convolution**.

> [!example] Worked Example 16.20 — triangle from rectangles
> Show that $\operatorname{rect}(t/T)*\operatorname{rect}(t/T)=T\,\Lambda(t/T)$ and deduce the transform of the triangle.
>
> **Time domain.** The flipped rectangle $\operatorname{rect}((t-\tau)/T)$ overlaps $\operatorname{rect}(\tau/T)$ on an interval of length $T-\lvert t\rvert$ when $\lvert t\rvert<T$ and not at all otherwise. So the convolution is $\max(0,T-\lvert t\rvert)=T\Lambda(t/T)$, a triangle of height $T$ and base $2T$.
>
> **Frequency domain.** By the convolution theorem, $T\Lambda(t/T)\leftrightarrow\big[T\operatorname{sinc}(\omega T/2\pi)\big]^2$, hence
>
> $$
> \Lambda(t/T)\leftrightarrow T\operatorname{sinc}^2(\omega T/2\pi),
> $$
>
> confirming Worked Example 16.11 by a third route. The sinc² decays as $1/\omega^2$ rather than $1/\omega$ because the triangle is continuous while the rectangle has jumps — **convolution smooths, and smoother signals have faster-decaying spectra.** Each further convolution with the rectangle multiplies the spectrum by another sinc and adds one more order of smoothness; the $n$-fold convolution tends to a Gaussian shape (the central limit theorem of [Ch. 19](../part-6-probability-statistics/19-probability.md) in transform form).

> [!example] Worked Example 16.21 — an RC filter driven by a pulse
> A unit pulse of width $T$ starting at $t=0$, $x(t)=u(t)-u(t-T)$, drives the RC low-pass filter $H(\omega)=1/(1+j\omega\tau)$, $\tau=RC$. Find the output.
>
> **Step 1 — transforms.** From Worked Example 16.10, $X(\omega)=(1-e^{-j\omega T})/(j\omega)$. So
>
> $$
> Y(\omega)=\left(1-e^{-j\omega T}\right)G_0(\omega),\qquad G_0(\omega)=\frac{1}{j\omega(1+j\omega\tau)} .
> $$
>
> **Step 2 — recognise the step response.** The step response $g(t)=h*u$ has transform $H(\omega)[\pi\delta(\omega)+1/(j\omega)]=\pi\delta(\omega)+G_0(\omega)$. In $Y$ the factor $1-e^{-j\omega T}$ vanishes at $\omega=0$, so adding the impulse terms changes nothing: $Y=(1-e^{-j\omega T})[\pi\delta+G_0]=G(\omega)-e^{-j\omega T}G(\omega)$, i.e.
>
> $$
> y(t)=g(t)-g(t-T).
> $$
>
> **Step 3 — find $g$.** Partial fractions, $G_0=\dfrac{1}{j\omega}-\dfrac{\tau}{1+j\omega\tau}$, so $G=\pi\delta+\dfrac{1}{j\omega}-\dfrac{1}{1/\tau+j\omega}$, giving $g(t)=u(t)-e^{-t/\tau}u(t)=(1-e^{-t/\tau})u(t)$.
>
> **Result.**
>
> $$
> y(t)=\left(1-e^{-t/\tau}\right)u(t)-\left(1-e^{-(t-T)/\tau}\right)u(t-T).
> $$
>
> With $\tau=1$ ms and $T=2$ ms the output reaches $1-e^{-2}=0.865$ at the end of the pulse and then decays as $0.865\,e^{-(t-T)/\tau}$. Check: the output area is $X(0)H(0)=T\cdot1=2$ ms, and integrating the result gives $\int_0^T(1-e^{-t/\tau})dt+0.865\tau=(T-\tau+\tau e^{-2})+0.865\tau=T$ ✓.

```desmos-graph
left=-0.5; right=6; bottom=-0.2; top=1.2
---
y=\left\{0<x<2:1,0\right\}|grey|dashed
y=\left\{0<x<2:1-e^{-x}\right\}|blue
y=\left\{x>2:(1-e^{-2})e^{-(x-2)}\right\}|blue
```

*Rectangular input of width $T=2\tau$ (grey dashed) and RC filter output (blue), time in units of $\tau$: the output rises to $0.865$ and decays exponentially — the low-pass filter has rounded the edges, removing the slow-decaying high-frequency content of the rectangle.*

> [!example] Worked Example 16.22 — the multiplication theorem
> Find the transform of $x(t)=\left(\dfrac{\sin Wt}{\pi t}\right)^2$.
>
> From Worked Example 16.5, $\sin(Wt)/(\pi t)\leftrightarrow\operatorname{rect}(\omega/2W)$, a rectangle of width $2W$ in $\omega$. By the multiplication theorem and Worked Example 16.20 (with $T=2W$, variable $\omega$),
>
> $$
> X(\omega)=\frac{1}{2\pi}\operatorname{rect}\!\left(\frac{\omega}{2W}\right)*\operatorname{rect}\!\left(\frac{\omega}{2W}\right)=\frac{2W}{2\pi}\,\Lambda\!\left(\frac{\omega}{2W}\right)=\frac{W}{\pi}\Lambda\!\left(\frac{\omega}{2W}\right).
> $$
>
> The squared sinc occupies **twice** the bandwidth of the sinc ($\lvert\omega\rvert<2W$): squaring a signal doubles its bandwidth, which is why nonlinear devices generate new frequencies.
>
> **Check** by $X(0)=\int x\,dt$: $\int\sin^2(Wt)/(\pi^2t^2)\,dt=\pi W/\pi^2=W/\pi$ (using $\int\sin^2(u)/u^2du=\pi$) ✓.

> [!example] Worked Example 16.23 — modulation and coherent demodulation
> A message $m(t)$ band-limited to $\lvert\omega\rvert<\omega_M$ is transmitted as $s(t)=m(t)\cos\omega_ct$ with $\omega_c>\omega_M$ (double-sideband suppressed-carrier modulation). At the receiver it is multiplied again by $\cos\omega_ct$. Show that a low-pass filter recovers $m$.
>
> **Transmitter.** $S(\omega)=\tfrac12[M(\omega-\omega_c)+M(\omega+\omega_c)]$: two half-height copies of $M$ centred on $\pm\omega_c$, occupying $\omega_c-\omega_M<\lvert\omega\rvert<\omega_c+\omega_M$.
>
> **Receiver.** $r(t)=s(t)\cos\omega_ct=m(t)\cos^2\omega_ct=\tfrac12m(t)+\tfrac12m(t)\cos2\omega_ct$, so
>
> $$
> R(\omega)=\tfrac12M(\omega)+\tfrac14\left[M(\omega-2\omega_c)+M(\omega+2\omega_c)\right].
> $$
>
> The baseband term occupies $\lvert\omega\rvert<\omega_M$ and the other two occupy $2\omega_c-\omega_M<\lvert\omega\rvert<2\omega_c+\omega_M$. Since $\omega_c>\omega_M$, the gap $\omega_M<\lvert\omega\rvert<2\omega_c-\omega_M$ separates them, and an ideal low-pass filter with gain 2 and any cutoff in that gap outputs exactly $m(t)$.
>
> If the receiver's oscillator has a phase error, $\cos(\omega_ct+\theta)$, the baseband term becomes $\tfrac12m(t)\cos\theta$: the output fades as $\theta\to90^\circ$. This is why coherent receivers need carrier synchronisation — a conclusion reached entirely by the multiplication theorem.

> [!tip] Engineering Note — where the two theorems appear
> - *Filtering* is multiplication of spectra: cascaded amplifier stages multiply gains, and in dB they add.
> - *Mixing* (multiplying by a local oscillator) shifts spectra; superheterodyne receivers, frequency synthesisers and lock-in amplifiers are all the modulation property.
> - *Windowing* a measured record (Section 16.6) and *sampling* (Section 16.8) are both multiplications in time and hence convolutions in frequency.
> - Fast convolution in digital filters computes $x*h$ as $\mathcal{F}^{-1}\{XH\}$ with the FFT ([Ch. 24](../part-7-numerical-computational/24-computational-mathematics.md)), because a product is far cheaper than a convolution sum.

## 16.6 Spectral analysis

Spectral analysis asks: *which frequencies does a signal contain, and how strongly?* The Fourier transform answers it exactly for idealised signals defined for all time. A measured signal, however, is only ever observed for a finite time, and that limitation is itself a multiplication in time, so it changes the spectrum we see in a precisely calculable way. This section unifies periodic and aperiodic signals under one transform and then quantifies the effect of finite observation: resolution, leakage and windows.

### Transforms of periodic signals

A periodic signal with Fourier series $x(t)=\sum_kc_ke^{jk\omega_0t}$ ([Ch. 15](15-fourier-series.md)) has, term by term from $e^{jk\omega_0t}\leftrightarrow2\pi\delta(\omega-k\omega_0)$, the generalized transform

$$
X(\omega)=2\pi\sum_{k=-\infty}^{\infty}c_k\,\delta(\omega-k\omega_0),
$$

a **line spectrum** of impulses at the harmonics, with areas $2\pi c_k$.

There is a second, illuminating derivation. Let $x_1(t)$ be one period of $x$ (zero outside $\lvert t\rvert<T_0/2$), with transform $X_1(\omega)$. Then $x$ is $x_1$ convolved with an impulse train, $x=x_1*\sum_n\delta(t-nT_0)$, and the impulse-train pair (proved in Section 16.8) gives

$$
X(\omega)=X_1(\omega)\cdot\omega_0\sum_k\delta(\omega-k\omega_0)=2\pi\sum_k\frac{X_1(k\omega_0)}{T_0}\,\delta(\omega-k\omega_0).
$$

Comparing, $c_k=X_1(k\omega_0)/T_0$ — the result found heuristically in Section 16.1, now as a theorem. **Periodicity in time is sampling in frequency.** The single-period spectrum $X_1(\omega)$ is the *envelope* of the line spectrum.

Aperiodic signals have continuous spectra; periodic ones have discrete spectra; a signal with both parts (a pulse plus a sinusoid, or a periodic signal plus a transient) has a mixed spectrum, and the transform handles all three with the same formulas.

> [!example] Worked Example 16.24 — a rectangular pulse train
> A clock signal consists of unit pulses of width $\tau=1$ ms repeating every $T_0=4$ ms (duty cycle $25\,\%$). Find its line spectrum.
>
> One period is $x_1=\operatorname{rect}(t/\tau)$, with $X_1(\omega)=\tau\operatorname{sinc}(\omega\tau/2\pi)$. With $\omega_0=2\pi/T_0$,
>
> $$
> c_k=\frac{X_1(k\omega_0)}{T_0}=\frac{\tau}{T_0}\operatorname{sinc}\!\left(\frac{k\tau}{T_0}\right)=0.25\operatorname{sinc}(0.25k).
> $$
>
> | $k$ | $0$ | $\pm1$ | $\pm2$ | $\pm3$ | $\pm4$ | $\pm5$ |
> |---|---|---|---|---|---|---|
> | $c_k$ | $0.2500$ | $0.2251$ | $0.1592$ | $0.0750$ | $0$ | $-0.0450$ |
>
> The DC value $c_0=0.25$ is the average of the waveform ✓. Every fourth harmonic vanishes because it falls on a null of the envelope ($k\tau/T_0$ an integer). The transform is $X(\omega)=2\pi\sum_kc_k\delta(\omega-2\pi k/T_0)$: lines every 250 Hz under a sinc envelope whose first null is at $1/\tau=1$ kHz. Making the pulses narrower (smaller duty cycle) spreads the envelope and makes the lines denser relative to it; making the period longer packs the lines closer together under the same envelope shape — and in the limit $T_0\to\infty$ they merge into the continuous spectrum of a single pulse.

### Finite observation: windowing and leakage

Observing $x(t)$ only for $\lvert t\rvert<T/2$ is multiplication by a **window** $w(t)=\operatorname{rect}(t/T)$. By the multiplication theorem the spectrum we actually compute is

$$
X_w(\omega)=\mathcal{F}\{x(t)w(t)\}=\frac{1}{2\pi}(X*W)(\omega),\qquad W(\omega)=T\operatorname{sinc}\!\left(\frac{\omega T}{2\pi}\right).
$$

The true spectrum is blurred by the window's transform. A single spectral line, $2\pi c\,\delta(\omega-\omega_1)$, becomes $c\,W(\omega-\omega_1)$: an impulse is smeared into a sinc.

> [!example] Worked Example 16.25 — a truncated cosine
> Find the spectrum of $x(t)=\cos(\omega_0t)\operatorname{rect}(t/T)$.
>
> With $\cos\omega_0t\leftrightarrow\pi[\delta(\omega-\omega_0)+\delta(\omega+\omega_0)]$ and convolution with an impulse being a shift,
>
> $$
> X_w(\omega)=\frac{T}{2}\left[\operatorname{sinc}\!\left(\frac{(\omega-\omega_0)T}{2\pi}\right)+\operatorname{sinc}\!\left(\frac{(\omega+\omega_0)T}{2\pi}\right)\right]
> =\frac{\sin\big((\omega-\omega_0)T/2\big)}{\omega-\omega_0}+\frac{\sin\big((\omega+\omega_0)T/2\big)}{\omega+\omega_0}.
> $$
>
> Instead of two impulses we see two sinc lobes of peak $T/2$ (growing with observation time, as an impulse is approached) and main-lobe width $4\pi/T$, plus **sidelobes** whose largest is $0.217$ of the peak ($-13.3$ dB), located about $1.43$ bins ($1.43\times2\pi/T$ rad/s) from the centre. This spreading of energy away from the true frequency is called **spectral leakage**.

```desmos-graph
left=-12; right=12; bottom=-1; top=2.5
---
y=\frac{\sin(2(x-5))}{x-5}+\frac{\sin(2(x+5))}{x+5}|blue
y=0|black|dotted
x=5|red|dashed
x=-5|red|dashed
```

*Spectrum of $\cos(5t)$ observed for $T=4$ s: the ideal lines at $\omega=\pm5$ (red dashed) become sinc lobes of width $4\pi/T\approx3.1$ rad/s with sidelobe leakage.*

Two consequences govern every practical spectrum analysis.

- **Resolution.** Two sinusoids closer than about one main-lobe half-width, $2\pi/T$ rad/s (i.e. $1/T$ Hz), merge into a single lobe. Finer resolution requires a longer record: **frequency resolution is the reciprocal of observation time.** This is the time–bandwidth trade-off of Section 16.4 applied to measurement.
- **Leakage.** Sidelobes of a strong component can mask a weak component nearby. The rectangular window's sidelobes decay only as $1/\lvert\omega\rvert$ because the window has jumps.

### Tapered windows

Leakage is reduced by replacing the abrupt rectangle with a window that goes smoothly to zero at its ends. By the smoothness theorem of Section 16.4, a smoother window has faster-decaying sidelobes. The standard example is the **Hann (raised-cosine) window**

$$
w_H(t)=\tfrac12\left[1+\cos\frac{2\pi t}{T}\right]\operatorname{rect}\!\left(\frac tT\right).
$$

Its transform follows from the modulation property: the rectangle contributes $\tfrac12W(\omega)$ and the cosine term contributes two shifted half-height copies $\tfrac14W(\omega\mp2\pi/T)$. With $u=\omega T/2\pi$,

$$
W_H(\omega)=\frac T2\operatorname{sinc}(u)+\frac T4\left[\operatorname{sinc}(u-1)+\operatorname{sinc}(u+1)\right]=\frac T2\,\frac{\operatorname{sinc}(u)}{1-u^2}.
$$

The three shifted sincs are arranged so that their sidelobes largely cancel. Compared with the rectangle:

- the first null moves from $u=1$ to $u=2$ (the main lobe is **twice as wide**);
- the highest sidelobe falls from $-13.3$ dB to $-31.5$ dB;
- the sidelobes decay as $1/\lvert u\rvert^3$ instead of $1/\lvert u\rvert$ (the Hann window is continuous with a continuous first derivative);
- the peak is $T/2$ instead of $T$ (the window's area is halved).

| Window | First null ($u$) | 3-dB half-width ($u$) | Highest sidelobe | Sidelobe decay |
|---|---|---|---|---|
| Rectangular | 1 | 0.44 | $-13.3$ dB | $1/\lvert u\rvert$ |
| Hann | 2 | 0.72 | $-31.5$ dB | $1/\lvert u\rvert^3$ |
| Hamming ($0.54+0.46\cos$) | 2 | 0.65 | $-42.7$ dB | $1/\lvert u\rvert$ |
| Blackman | 3 | 0.82 | $-58.1$ dB | $1/\lvert u\rvert^3$ |

The trade-off is universal: **lower sidelobes cost a wider main lobe**, i.e. better dynamic range costs resolution, because a window cannot be both short in time and narrow in frequency.

```desmos-graph
left=-5; right=5; bottom=-0.3; top=1.1
---
y=\frac{\sin(\pi x)}{\pi x}|blue
y=\frac{\sin(\pi x)}{\pi x(1-x^2)}|red
y=0|black|dotted
```

*Normalised window transforms against $u=\omega T/2\pi$: rectangular $\operatorname{sinc}(u)$ (blue) and Hann $\operatorname{sinc}(u)/(1-u^2)$ (red). The Hann main lobe is twice as wide, but its sidelobes are almost invisible on this scale.*

> [!example] Worked Example 16.26 — choosing a record length and window
> A signal contains a strong tone at $1000$ Hz and a weak tone, $40$ dB below it, at $1200$ Hz. The analysis record is $T=0.1$ s long. Can the weak tone be seen (a) with a rectangular window, (b) with a Hann window?
>
> The frequency offset is $\Delta f=200$ Hz, i.e. $u=\Delta f\,T=20$ bins from the strong tone.
>
> *(a) Rectangular.* The sidelobe envelope is $\lvert\operatorname{sinc}(u)\rvert\le1/(\pi u)$; at $u=20$ this is $1/(20\pi)=0.0159$, i.e. $-36.0$ dB. The leakage from the strong tone at $1200$ Hz is up to $-36$ dB, *above* the weak tone at $-40$ dB. The weak tone is masked.
>
> *(b) Hann.* The envelope is $1/(\pi u(u^2-1))$; at $u=20$ this is $4.0\times10^{-5}$, i.e. $-88.0$ dB. The weak tone stands about 48 dB above the leakage floor and is clearly visible. Its main lobe (null-to-null $4/T=40$ Hz) is far narrower than the 200 Hz separation, so the wider Hann lobe costs nothing here.
>
> *(c) Resolution check.* If instead the second tone were at $1015$ Hz with equal amplitude, $u=1.5$. With the rectangle the two main lobes (half-width 1 bin) are just separable; with Hann (half-width 2 bins) they merge. Either lengthen the record to $T\ge0.2$ s, or accept the rectangle's leakage. Resolution and dynamic range compete, and the record length is the only free lunch.

> [!warning] Common Pitfall — reading peak heights from windowed spectra
> The height of a windowed line depends on the window: it is $c\,W(0)$, i.e. proportional to the window's area ($T$ for rectangular, $T/2$ for Hann). Amplitude estimates must be divided by the window's area (its "coherent gain"), and the frequency of a tone lying between computed frequency points is underestimated in height ("scalloping"). These are consequences of $X_w=\frac1{2\pi}X*W$, not instrument faults.

### Time-varying spectra

The Fourier transform integrates over all time, so it cannot say *when* a frequency occurred: a chirp and a time-reversed chirp have the same magnitude spectrum. Sliding a short window along the signal and transforming each segment gives the **short-time Fourier transform**

$$
X(\tau,\omega)=\int_{-\infty}^{\infty}x(t)\,w(t-\tau)\,e^{-j\omega t}\,dt ,
$$

whose magnitude squared (the spectrogram) is the standard picture of speech, radar returns and vibration. Its time and frequency resolutions are bounded by the window's duration and bandwidth, and the uncertainty principle forbids making both arbitrarily fine.

> [!tip] Engineering Note — from continuous theory to the FFT
> A spectrum analyser's resolution bandwidth, and an FFT's bin spacing $1/T$ Hz, are both manifestations of $W(\omega)$. The discrete computation (DFT/FFT) evaluates samples of $X_w$, and its windowing behaviour is exactly the continuous theory above; the algorithms themselves are treated in [Ch. 24](../part-7-numerical-computational/24-computational-mathematics.md).

## 16.7 Energy spectral density

### Energy signals and power signals

How much energy does a pulse deliver, and how is that energy distributed over frequency? These questions decide how much bandwidth a channel must have to pass "most" of a pulse, how much noise a filter admits, and how much interference a transmitter puts into a neighbouring band.

> [!definition] Energy and power of a signal
> The **energy** of $x(t)$ is
>
> $$
> E=\int_{-\infty}^{\infty}\lvert x(t)\rvert^2\,dt .
> $$
>
> If $x$ is a voltage across (or a current through) a $1\ \Omega$ resistor, $E$ is the energy dissipated in joules; for a real load $R$ divide (or multiply) by $R$. The **average power** is $P=\lim_{T\to\infty}\frac1T\int_{-T/2}^{T/2}\lvert x(t)\rvert^2dt$. A signal with $0<E<\infty$ is an **energy signal** (and has $P=0$); a signal with $0<P<\infty$ is a **power signal** (and has $E=\infty$).

Pulses and transients are energy signals; sinusoids, periodic waveforms, and stationary noise are power signals. This section treats energy signals; the power-signal theory is outlined at the end and developed for random signals in [Ch. 21](../part-6-probability-statistics/21-random-processes.md).

### Parseval's theorem

> [!theorem] Parseval's (Rayleigh's) energy theorem
> If $x$ has finite energy, then
>
> $$
> E=\int_{-\infty}^{\infty}\lvert x(t)\rvert^2dt=\frac{1}{2\pi}\int_{-\infty}^{\infty}\lvert X(\omega)\rvert^2d\omega .
> $$
>
> More generally, for two finite-energy signals,
>
> $$
> \int_{-\infty}^{\infty}x(t)\,y^*(t)\,dt=\frac{1}{2\pi}\int_{-\infty}^{\infty}X(\omega)\,Y^*(\omega)\,d\omega .
> $$

*Proof.* Substitute the inverse formula for $x$ and exchange the order of integration:

$$
\int x(t)y^*(t)\,dt=\int\left[\frac1{2\pi}\int X(\omega)e^{j\omega t}d\omega\right]y^*(t)\,dt=\frac{1}{2\pi}\int X(\omega)\left[\int y(t)e^{-j\omega t}dt\right]^*d\omega=\frac{1}{2\pi}\int X(\omega)Y^*(\omega)\,d\omega .
$$

Setting $y=x$ gives the energy theorem. (Alternatively, apply the multiplication theorem to $x(t)x^*(t)$, using $x^*\leftrightarrow X^*(-\omega)$, and evaluate the result at $\omega=0$.) $\square$

The general form says the Fourier transform preserves inner products (up to the factor $2\pi$): it is a *rotation* of function space, in the language of [Ch. 15](15-fourier-series.md)'s orthogonal expansions. In particular, two signals are orthogonal in time if and only if their transforms are orthogonal in frequency; signals occupying disjoint frequency bands are automatically orthogonal, which is the mathematical basis of frequency-division multiplexing.

### Energy spectral density

> [!definition] Energy spectral density
> The **energy spectral density (ESD)** of $x$ is
>
> $$
> \Psi_x(\omega)=\lvert X(\omega)\rvert^2 .
> $$
>
> By Parseval, $\Psi_x(\omega)\,d\omega/2\pi$ is the energy contained in the frequency band $d\omega$ (equivalently, $\Psi_x\,df$ in a band $df$ Hz).

Units: for a voltage, $\lvert X\rvert^2$ is in V²·s² = (V²·s)/Hz, i.e. joules (into 1 Ω) per hertz. The ESD discards phase, so many different signals share the same ESD; a pulse and its delayed or time-reversed version, for example.

For **real** $x$ the ESD is even, so the energy in the band $\omega_1<\lvert\omega\rvert<\omega_2$ — counting both the positive and negative frequencies that make up the real signal components — is

$$
E_{[\omega_1,\omega_2]}=\frac{1}{2\pi}\cdot2\int_{\omega_1}^{\omega_2}\lvert X(\omega)\rvert^2d\omega=\frac{1}{\pi}\int_{\omega_1}^{\omega_2}\lvert X(\omega)\rvert^2d\omega .
$$

**ESD through an LTI system.** Since $Y=HX$,

$$
\Psi_y(\omega)=\lvert H(\omega)\rvert^2\,\Psi_x(\omega),\qquad E_y=\frac{1}{2\pi}\int\lvert H(\omega)\rvert^2\lvert X(\omega)\rvert^2d\omega .
$$

The output ESD depends only on the **magnitude** response; phase distortion reshapes a pulse but does not change its energy or the energy in any band. The quantity $\lvert H\rvert^2$ is the **power gain**, and an ideal band-pass filter with $\lvert H\rvert^2=1$ on a band and $0$ elsewhere measures the energy in that band — which is how a spectrum analyser physically measures the ESD.

### Worked examples

> [!example] Worked Example 16.27 — energy and essential bandwidth of an exponential pulse
> For $x(t)=e^{-at}u(t)$: directly, $E=\int_0^\infty e^{-2at}dt=\dfrac{1}{2a}$. By Parseval,
>
> $$
> E=\frac{1}{2\pi}\int_{-\infty}^{\infty}\frac{d\omega}{a^2+\omega^2}=\frac{1}{2\pi}\cdot\frac{\pi}{a}=\frac{1}{2a}\ ✓
> $$
>
> The energy in $\lvert\omega\rvert<W$ is $\frac{1}{\pi}\int_0^W\frac{d\omega}{a^2+\omega^2}=\frac{1}{\pi a}\arctan\frac{W}{a}$, a fraction
>
> $$
> \frac{E_W}{E}=\frac{2}{\pi}\arctan\frac{W}{a}.
> $$
>
> At the 3-dB frequency $W=a$ this is exactly $50\,\%$. For $95\,\%$: $\arctan(W/a)=0.95\pi/2$, so $W=a\tan(85.5^\circ)\approx12.7a$. For $99\,\%$, $W=a\tan(89.1^\circ)\approx63.7a$. The slow $1/\omega^2$ decay of the ESD (caused by the jump at $t=0$) makes the essential bandwidth very much larger than the 3-dB bandwidth.
>
> With $a=10^4$ s$^{-1}$ (a 100 μs time constant), the 95 % essential bandwidth is $1.27\times10^5$ rad/s, about 20 kHz.

```desmos-graph
left=0; right=15; bottom=0; top=1.1
---
y=\frac{1}{1+x^2}|blue
y=\frac{2}{\pi}\arctan(x)|red
y=0.95|grey|dashed
(12.7,0.95)|label:95% energy
```

*For $x=e^{-at}u(t)$ with $\omega$ in units of $a$: normalised ESD $a^2\lvert X\rvert^2=1/(1+\omega^2)$ (blue) and cumulative energy fraction $\tfrac{2}{\pi}\arctan\omega$ (red).*

> [!example] Worked Example 16.28 — the rectangular pulse and a sinc integral
> *(a)* Use Parseval to evaluate $\displaystyle\int_{-\infty}^{\infty}\frac{\sin^2u}{u^2}\,du$.
>
> For $x=\operatorname{rect}(t)$, $E=1$ and $X=2\sin(\omega/2)/\omega$. Parseval gives $1=\frac1{2\pi}\int\frac{4\sin^2(\omega/2)}{\omega^2}d\omega$. Put $u=\omega/2$, $d\omega=2du$: $1=\frac{1}{2\pi}\int\frac{4\sin^2u}{4u^2}2\,du$, so
>
> $$
> \int_{-\infty}^{\infty}\frac{\sin^2u}{u^2}\,du=\pi .
> $$
>
> *(b)* What fraction of the energy of a rectangular pulse of width $T$ lies within its main lobe $\lvert\omega\rvert<2\pi/T$?
>
> $$
> \frac{E_{\text{main}}}{E}=\frac{1}{T}\cdot\frac{1}{\pi}\int_0^{2\pi/T}\frac{4\sin^2(\omega T/2)}{\omega^2}d\omega=\frac{2}{\pi}\int_0^{\pi}\frac{\sin^2v}{v^2}dv\approx0.903,
> $$
>
> after substituting $v=\omega T/2$. The integral has no elementary closed form; numerical quadrature gives $\int_0^\pi\sin^2v/v^2\,dv=1.4182$, hence the fraction $0.903$. So **about 90 % of the energy of a rectangular pulse lies in its main lobe**, even though its RMS bandwidth is infinite (Worked Example 16.18). A channel of bandwidth $1/T$ Hz passes 90 % of the energy of a pulse of duration $T$; it is the source of the "bandwidth ≈ 1/duration" rule.

> [!example] Worked Example 16.29 — energy at the output of a filter
> The pulse $x(t)=e^{-at}u(t)$ drives an RC low-pass filter $H=b/(b+j\omega)$, $b=1/RC$. Find the output energy by Parseval and check it in the time domain.
>
> **Frequency domain.** $\lvert H\rvert^2\lvert X\rvert^2=\dfrac{b^2}{(b^2+\omega^2)(a^2+\omega^2)}$. With the partial-fraction identity
>
> $$
> \frac{1}{(a^2+\omega^2)(b^2+\omega^2)}=\frac{1}{b^2-a^2}\left[\frac{1}{a^2+\omega^2}-\frac{1}{b^2+\omega^2}\right]
> $$
>
> and $\int d\omega/(c^2+\omega^2)=\pi/c$, the integral is $\dfrac{\pi}{b^2-a^2}\left(\dfrac1a-\dfrac1b\right)=\dfrac{\pi}{ab(a+b)}$. Hence
>
> $$
> E_y=\frac{1}{2\pi}\cdot b^2\cdot\frac{\pi}{ab(a+b)}=\frac{b}{2a(a+b)} .
> $$
>
> **Time domain.** For $a\ne b$, $y=x*h=\dfrac{b}{b-a}\left(e^{-at}-e^{-bt}\right)u(t)$, and
>
> $$
> E_y=\frac{b^2}{(b-a)^2}\left[\frac{1}{2a}-\frac{2}{a+b}+\frac{1}{2b}\right]=\frac{b^2}{(b-a)^2}\cdot\frac{(b-a)^2}{2ab(a+b)}=\frac{b}{2a(a+b)}\ ✓
> $$
>
> **Numbers.** $a=500$ s$^{-1}$, $b=1000$ s$^{-1}$: $E_x=1/(2a)=1.000\times10^{-3}$, $E_y=1000/(2\cdot500\cdot1500)=6.67\times10^{-4}$ (V²·s). The filter passes $E_y/E_x=b/(a+b)=2/3$ of the input energy. As $b\to\infty$ (no filtering) the ratio tends to 1, and as $b\to0$ to 0, as expected.

> [!example] Worked Example 16.30 — essential bandwidth of a Gaussian pulse
> For $x=e^{-t^2/(2\sigma^2)}$, $\lvert X\rvert^2=2\pi\sigma^2e^{-\sigma^2\omega^2}$. The fraction of energy in $\lvert\omega\rvert<W$ is
>
> $$
> \frac{E_W}{E}=\frac{\int_{-W}^{W}e^{-\sigma^2\omega^2}d\omega}{\int_{-\infty}^{\infty}e^{-\sigma^2\omega^2}d\omega}=\operatorname{erf}(\sigma W).
> $$
>
> For $95\,\%$, $\sigma W=\operatorname{erf}^{-1}(0.95)=1.386$, so $W=1.386/\sigma$. Compare the exponential pulse (Worked Example 16.27), whose 95 % bandwidth is 12.7 times its 3-dB bandwidth: the Gaussian's 3-dB bandwidth (where $e^{-\sigma^2\omega^2/2}=1/\sqrt2$) is $\sqrt{\ln2}/\sigma=0.833/\sigma$, only a factor $1.66$ below its 95 % bandwidth. Smooth pulses have compact spectra.

### Autocorrelation and the energy Wiener–Khinchin relation

> [!definition] Deterministic autocorrelation
> The autocorrelation of a finite-energy signal is
>
> $$
> R_x(\tau)=\int_{-\infty}^{\infty}x(t+\tau)\,x^*(t)\,dt .
> $$
>
> It measures the similarity of $x$ with a copy of itself shifted by $\tau$. For real $x$ it is real and even, and $\lvert R_x(\tau)\rvert\le R_x(0)=E$ by the Cauchy–Schwarz inequality.

Substituting $s=-t$ shows $R_x(\tau)=x(\tau)*x^*(-\tau)$. Since $x^*(-t)\leftrightarrow X^*(\omega)$, the convolution theorem gives

> [!theorem] Autocorrelation–ESD pair
>
> $$
> R_x(\tau)\ \leftrightarrow\ \lvert X(\omega)\rvert^2=\Psi_x(\omega),\qquad R_x(0)=E=\frac{1}{2\pi}\int\Psi_x(\omega)\,d\omega .
> $$

This is the energy-signal form of the **Wiener–Khinchin theorem**: the ESD is the Fourier transform of the autocorrelation. Because $\Psi_x$ carries no phase information, neither does $R_x$; any two signals with the same magnitude spectrum have the same autocorrelation. Similarly the **cross-correlation** $R_{xy}(\tau)=\int x(t+\tau)y^*(t)\,dt$ transforms to the **cross-energy spectral density** $X(\omega)Y^*(\omega)$.

> [!example] Worked Example 16.31 — autocorrelation of an exponential pulse
> For $x=e^{-at}u(t)$ and $\tau\ge0$: $R_x(\tau)=\int_0^\infty e^{-a(t+\tau)}e^{-at}dt=e^{-a\tau}/(2a)$. By evenness $R_x(\tau)=e^{-a\lvert\tau\rvert}/(2a)$. Its transform is $\dfrac{1}{2a}\cdot\dfrac{2a}{a^2+\omega^2}=\dfrac{1}{a^2+\omega^2}=\lvert X(\omega)\rvert^2$ ✓, and $R_x(0)=1/(2a)=E$ ✓.

### Power signals

For a power signal the ESD is infinite, and we work with power per unit bandwidth instead. Truncate the signal to $x_T(t)=x(t)\operatorname{rect}(t/T)$, compute its ESD $\lvert X_T(\omega)\rvert^2$, divide by the duration and let $T\to\infty$:

$$
S_x(\omega)=\lim_{T\to\infty}\frac{\lvert X_T(\omega)\rvert^2}{T},\qquad P=\frac{1}{2\pi}\int_{-\infty}^{\infty}S_x(\omega)\,d\omega .
$$

This is the **power spectral density**. For a periodic signal it reduces to Parseval's theorem for Fourier series, $P=\sum_k\lvert c_k\rvert^2$, with $S_x(\omega)=2\pi\sum_k\lvert c_k\rvert^2\delta(\omega-k\omega_0)$. For random signals, where the limit must be combined with an expectation, the PSD $S_{XX}(\omega)$ and its relation to the autocorrelation $R_{XX}(\tau)$ are the subject of [Ch. 21](../part-6-probability-statistics/21-random-processes.md).

> [!warning] Common Pitfall — one-sided versus two-sided spectra
> Integrating $\lvert X\rvert^2$ over $\omega>0$ only and forgetting the factor 2 halves every energy figure. Conversely, many instruments display a *one-sided* density that already includes the factor 2. Before comparing numbers, establish whether a density is per rad/s or per Hz (they differ by $2\pi$) and one- or two-sided (they differ by 2).

> [!tip] Engineering Note — matched filtering
> The generalized Parseval theorem underlies the **matched filter**: among all filters, the one with $H(\omega)=X^*(\omega)e^{-j\omega t_0}$ maximises the output peak relative to white noise, and its output at $t_0$ is $\frac{1}{2\pi}\int\lvert X\rvert^2d\omega=E$, the pulse energy. Its impulse response $x(t_0-t)$ makes the filter output the autocorrelation $R_x(t-t_0)$. Radar range resolution and digital-receiver detection both rest on these two facts.

## 16.8 Sampling-related mathematical foundations

Every digital system — a data converter, a software radio, a digital oscilloscope, a control loop running on a microcontroller — replaces a continuous signal $x(t)$ by its samples $x(nT)$. Whether anything is lost, and if so what, is a question about Fourier transforms. This section gives the mathematics of ideal (instantaneous) sampling: the spectrum of an impulse train, the spectral replication caused by sampling, the Nyquist condition under which the replicas do not overlap, aliasing when they do, and the interpolation formula that reconstructs $x(t)$. Converter hardware is outside the scope of this book; the discrete-time transforms that take over once the samples exist are the subject of [Ch. 17](17-z-transforms.md).

Notation: $T$ is the sampling interval, $f_s=1/T$ the sampling rate in Hz and $\omega_s=2\pi/T$ in rad/s.

### The impulse train and its transform

> [!theorem] Impulse train (Dirac comb)
> For $p(t)=\displaystyle\sum_{n=-\infty}^{\infty}\delta(t-nT)$,
>
> $$
> P(\omega)=\omega_s\sum_{k=-\infty}^{\infty}\delta(\omega-k\omega_s),\qquad \omega_s=\frac{2\pi}{T}.
> $$
>
> An impulse train in time, of spacing $T$ and unit areas, transforms to an impulse train in frequency, of spacing $\omega_s=2\pi/T$ and areas $\omega_s$.

*Proof.* $p(t)$ is periodic with period $T$. Its Fourier series coefficients are

$$
c_k=\frac1T\int_{-T/2}^{T/2}\delta(t)e^{-jk\omega_st}dt=\frac1T\quad\text{for every }k,
$$

so $p(t)=\frac1T\sum_ke^{jk\omega_st}$. Transforming term by term with $e^{jk\omega_st}\leftrightarrow2\pi\delta(\omega-k\omega_s)$ gives $P(\omega)=\frac{2\pi}{T}\sum_k\delta(\omega-k\omega_s)$. $\square$

The comb is (up to scale) its own transform, and its spacings are reciprocal: **dense sampling in time means widely spaced lines in frequency.** The series $p(t)=\frac1T\sum_ke^{jk\omega_st}$ does not converge in the ordinary sense; it is an identity between generalized functions, meaning that both sides give the same result when integrated against any smooth, rapidly decaying test function.

> [!theorem] Poisson summation formula
> For a sufficiently smooth and rapidly decaying $x(t)$ with transform $X(\omega)$,
>
> $$
> \sum_{n=-\infty}^{\infty}x(nT)\,e^{-j\omega nT}=\frac1T\sum_{k=-\infty}^{\infty}X(\omega-k\omega_s),
> $$
>
> and in particular (at $\omega=0$) $\displaystyle T\sum_nx(nT)=\sum_kX(k\omega_s)$.

This is the impulse-train theorem integrated against $x$, and it is the analytic heart of everything that follows: the left side is computed from samples, the right side from the continuous spectrum.

### Ideal sampling and spectral replication

Ideal (impulse) sampling of $x(t)$ at interval $T$ is modelled as multiplication by the comb:

$$
x_s(t)=x(t)\,p(t)=\sum_{n=-\infty}^{\infty}x(nT)\,\delta(t-nT).
$$

The product keeps only the values $x(nT)$, each as the area of an impulse. By the multiplication theorem,

$$
X_s(\omega)=\frac{1}{2\pi}X(\omega)*\omega_s\sum_k\delta(\omega-k\omega_s)=\frac{1}{T}\sum_{k=-\infty}^{\infty}X(\omega-k\omega_s).
$$

> [!theorem] Spectrum of a sampled signal
>
> $$
> X_s(\omega)=\frac1T\sum_{k=-\infty}^{\infty}X(\omega-k\omega_s)=\sum_{n=-\infty}^{\infty}x(nT)\,e^{-j\omega nT}.
> $$
>
> **Sampling replicates the spectrum at every multiple of $\omega_s$, scaled by $1/T$.** The sampled spectrum is periodic in $\omega$ with period $\omega_s$.

The second form follows by transforming $x_s$ term by term ($\delta(t-nT)\leftrightarrow e^{-j\omega nT}$); the equality of the two forms is the Poisson summation formula. The second form is the bridge to discrete time: with the normalised frequency $\Omega=\omega T$ (radians per sample), $X_s=\sum_nx[n]e^{-j\Omega n}$ is the **discrete-time Fourier transform** of the sequence $x[n]=x(nT)$, and with $z=e^{sT}$ (so $z=e^{j\Omega}$ on the $j\omega$-axis) it is the Z-transform of [Ch. 17](17-z-transforms.md) evaluated on the unit circle.

Periodicity in frequency is the dual of Section 16.6's statement that periodicity in time gives a line spectrum: **sampling in one domain is periodic repetition in the other.**

> [!example] Worked Example 16.32 — the sampled spectrum of a two-sided exponential
> Sample $x(t)=e^{-a\lvert t\rvert}$ (transform $2a/(a^2+\omega^2)$, not band-limited) at interval $T$. Find $X_s(\omega)$ in closed form and compare $T X_s$ with $X$.
>
> **From the samples.** With $r=e^{-aT}<1$, $x(nT)=r^{\lvert n\rvert}$, and two geometric series give
>
> $$
> X_s(\omega)=\sum_{n=0}^{\infty}r^ne^{-j\omega nT}+\sum_{n=1}^{\infty}r^ne^{j\omega nT}=\frac{1}{1-re^{-j\omega T}}+\frac{re^{j\omega T}}{1-re^{j\omega T}}=\frac{1-r^2}{1-2r\cos\omega T+r^2}.
> $$
>
> **From the replication formula** the same function must equal $\frac1T\sum_k\dfrac{2a}{a^2+(\omega-k\omega_s)^2}$. With $a=1$, $T=0.5$ s ($r=0.6065$) and $\omega=1.3$ rad/s, the closed form gives $1.57173$; summing the replicas numerically for $\lvert k\rvert\le20000$ gives $1.57173$ ✓.
>
> **Aliasing error.** Ideally $TX_s=X$ in the band $\lvert\omega\rvert<\omega_s/2$. Here:
>
> | $\omega$ | $X(\omega)$ | $TX_s(\omega)$ |
> |---|---|---|
> | $0$ | $2.000$ | $T\frac{1+r}{1-r}=2.041$ |
> | $\omega_s/2=2\pi$ rad/s | $0.0494$ | $T\frac{1-r}{1+r}=0.1225$ |
>
> Near DC the error is 2 %, but at the folding frequency the sampled spectrum is about 2.5 times the true one: the tails of the replicas centred at $\pm\omega_s$ have folded onto the band edge. A signal whose spectrum decays slowly suffers aliasing at any finite sampling rate.

### Nyquist condition and reconstruction

If $X(\omega)=0$ for $\lvert\omega\rvert>\omega_M$, the replicas $X(\omega-k\omega_s)$ occupy the intervals $\lvert\omega-k\omega_s\rvert<\omega_M$. Adjacent replicas do not overlap provided $\omega_s-\omega_M>\omega_M$, i.e. $\omega_s>2\omega_M$. In that case the baseband replica ($k=0$) can be isolated by a low-pass filter, and $x(t)$ recovered exactly.

> [!theorem] Sampling theorem (Shannon–Nyquist)
> If $x(t)$ is band-limited, $X(\omega)=0$ for $\lvert\omega\rvert>\omega_M$, and
>
> $$
> \omega_s>2\omega_M\qquad(\text{equivalently } f_s>2f_M),
> $$
>
> then $x(t)$ is uniquely determined by its samples $x(nT)$ and
>
> $$
> x(t)=\sum_{n=-\infty}^{\infty}x(nT)\operatorname{sinc}\!\left(\frac{t-nT}{T}\right)=\sum_{n=-\infty}^{\infty}x(nT)\,\frac{\sin\big(\pi(t-nT)/T\big)}{\pi(t-nT)/T}.
> $$
>
> $2\omega_M$ (or $2f_M$) is the **Nyquist rate**; $\omega_s/2$ (or $f_s/2$) is the **Nyquist (folding) frequency**.

*Proof.* Because the replicas do not overlap, $X(\omega)=T\,X_s(\omega)$ for $\lvert\omega\rvert<\omega_s/2$ and $X(\omega)=0$ for $\omega_M<\lvert\omega\rvert$. So multiplying $X_s$ by the ideal low-pass filter

$$
H_r(\omega)=T\operatorname{rect}\!\left(\frac{\omega}{\omega_s}\right)\qquad(\text{gain }T,\ \text{cutoff }\omega_s/2)
$$

returns exactly $X$: $X(\omega)=H_r(\omega)X_s(\omega)$. By Worked Example 16.5 with $W=\omega_s/2=\pi/T$, the impulse response of $H_r$ is

$$
h_r(t)=T\cdot\frac{\sin(\pi t/T)}{\pi t}=\operatorname{sinc}\!\left(\frac tT\right).
$$

In the time domain, $x=h_r*x_s=h_r*\sum_nx(nT)\delta(t-nT)=\sum_nx(nT)h_r(t-nT)$, which is the interpolation formula. $\square$

Two features of the interpolation formula deserve attention.

- **Interpolation property.** $\operatorname{sinc}((t-nT)/T)$ equals 1 at $t=nT$ and 0 at every other sample instant $t=mT$, $m\ne n$. So the formula passes exactly through every sample, and between samples each sinc contributes in proportion to its sample value.
- **Orthogonality.** The shifted sincs $\operatorname{sinc}((t-nT)/T)$ are mutually orthogonal, each with energy $T$. By the generalized Parseval theorem, with transforms $T\operatorname{rect}(\omega/\omega_s)e^{-j\omega nT}$, $\int\operatorname{sinc}\big(\tfrac{t-nT}{T}\big)\operatorname{sinc}\big(\tfrac{t-mT}{T}\big)dt=\frac{1}{2\pi}\int_{-\pi/T}^{\pi/T}T^2e^{-j\omega(n-m)T}d\omega=T\,\delta_{nm}$, where $\delta_{nm}$ is the Kronecker delta. The sampling theorem is therefore an orthogonal expansion, like a Fourier series, in which the coefficients are simply the sample values: $x(nT)=\frac1T\int x(t)\operatorname{sinc}((t-nT)/T)\,dt$. It follows that $E=T\sum_n\lvert x(nT)\rvert^2$ for a band-limited signal sampled above the Nyquist rate.

```desmos-graph
left=-3; right=5; bottom=-0.4; top=1.4
---
y=0.5\frac{\sin(\pi(x+1))}{\pi(x+1)}|grey|dashed
y=\frac{\sin(\pi x)}{\pi x}|grey|dashed
y=0.8\frac{\sin(\pi(x-1))}{\pi(x-1)}|grey|dashed
y=0.3\frac{\sin(\pi(x-2))}{\pi(x-2)}|grey|dashed
y=0.5\frac{\sin(\pi(x+1))}{\pi(x+1)}+\frac{\sin(\pi x)}{\pi x}+0.8\frac{\sin(\pi(x-1))}{\pi(x-1)}+0.3\frac{\sin(\pi(x-2))}{\pi(x-2)}|blue
(-1,0.5)|red
(0,1)|red
(1,0.8)|red
(2,0.3)|red
```

*Sinc interpolation with $T=1$: four samples (red), the weighted sincs $x(nT)\operatorname{sinc}(t-n)$ (grey dashed), and their sum (blue), which passes through every sample because each sinc vanishes at the other sample instants.*

```desmos-graph
left=-7; right=7; bottom=-0.3; top=1.3
---
y=\max(0,1-\left|x\right|)|blue
y=\max(0,1-\left|x-3\right|)|blue|dashed
y=\max(0,1-\left|x+3\right|)|blue|dashed
y=\max(0,1-\left|x-6\right|)|blue|dashed
y=\max(0,1-\left|x+6\right|)|blue|dashed
x=1.5|red|dotted
x=-1.5|red|dotted
```

*Replicated spectrum $TX_s(\omega)$ for a triangular $X$ with $\omega_M=1$ and $\omega_s=3>2\omega_M$: the replicas are separated by guard bands, so a low-pass filter with cutoff $\omega_s/2=1.5$ (red) recovers $X$ exactly.*

> [!example] Worked Example 16.33 — Nyquist rates of derived signals
> Let $x(t)=\dfrac{\sin(200\pi t)}{\pi t}$, band-limited to $\omega_M=200\pi$ rad/s ($f_M=100$ Hz). Find the Nyquist rate of each signal.
>
> *(a)* $x(t)$ itself: $f_s>2f_M=200$ Hz.
>
> *(b)* $x^2(t)$: by Worked Example 16.22 squaring doubles the bandwidth to $f_M=200$ Hz, so $f_s>400$ Hz.
>
> *(c)* $x(t)\cos(1000\pi t)$: modulation moves the spectrum to $400<\lvert f\rvert<600$ Hz. Treated as a low-pass signal with highest frequency 600 Hz, $f_s>1200$ Hz. (Band-pass sampling theory allows lower rates for such signals by exploiting the empty band below 400 Hz; the condition is that no replica lands on the occupied band.)
>
> *(d)* $x(t)*y(t)$ with $y$ band-limited to 50 Hz: the product $X(\omega)Y(\omega)$ vanishes outside the *narrower* band, so $f_s>100$ Hz.
>
> *(e)* $x(t)+x(t-t_0)$: delay does not change the band, so $f_s>200$ Hz.
>
> Multiplication widens bands (bandwidths add, by convolution of spectra); filtering and delay do not.

> [!example] Worked Example 16.34 — reconstruction from samples, and truncation
> The tone $x(t)=\cos(6\pi t)$ ($f_0=3$ Hz) is sampled at $f_s=10$ Hz ($T=0.1$ s), above its Nyquist rate of 6 Hz. Reconstruct $x(0.05)$, midway between two samples.
>
> The samples are $x(nT)=\cos(0.6\pi n)$ and
>
> $$
> x(0.05)=\sum_n\cos(0.6\pi n)\operatorname{sinc}(0.5-n).
> $$
>
> The exact value is $\cos(0.3\pi)=0.587785$. Summing over $\lvert n\rvert\le N$:
>
> | $N$ | 5 | 50 | 500 | 5000 |
> |---|---|---|---|---|
> | partial sum | 0.584465 | 0.587725 | 0.587785 | 0.587785 |
>
> The sum converges to the correct value, but slowly — the sinc tails decay only as $1/\lvert n\rvert$, so even 11 samples leave an error of $3\times10^{-3}$. Practical interpolators therefore use a faster-decaying kernel, which in turn requires sampling somewhat above the Nyquist rate so that a gentler filter edge fits into the guard band.

### Aliasing

If $\omega_s<2\omega_M$, adjacent replicas overlap. In the baseband $\lvert\omega\rvert<\omega_s/2$, $TX_s(\omega)=X(\omega)+X(\omega-\omega_s)+X(\omega+\omega_s)+\cdots$, and the extra terms are the high-frequency parts of $X$ **folded** about $\pm\omega_s/2$. The folded content is indistinguishable from genuine low-frequency content: no processing of the samples can separate them. This is **aliasing**.

For a single sinusoid the effect is explicit. A tone of frequency $f_0$ sampled at $f_s$ produces the same samples as one of frequency

$$
f_a=\lvert f_0-kf_s\rvert,\qquad k=\text{the integer nearest } f_0/f_s ,
$$

which lies in $[0,f_s/2]$. Indeed $\cos(2\pi f_0nT+\theta)=\cos(2\pi(f_0-kf_s)nT+\theta)$ because $2\pi kf_snT=2\pi kn$. When $f_0-kf_s<0$, writing the result in terms of the positive frequency $f_a$ reverses the sign of the phase: a sine aliases to a *negated* sine.

> [!example] Worked Example 16.35 — aliasing of a tone
> A 7 kHz tone is sampled at $f_s=10$ kHz. Here $f_0/f_s=0.7$, so $k=1$ and $f_a=\lvert7-10\rvert=3$ kHz.
>
> The samples are $\cos(2\pi\cdot0.7n)$ and $\cos(2\pi\cdot0.3n)$, which are identical since $0.7n=n-0.3n$ and cosine is even: both sequences are $1,\,-0.309,\,-0.809,\,0.809,\,0.309,\,-1,\dots$ For a sine, $\sin(2\pi\cdot0.7n)=-\sin(2\pi\cdot0.3n)$: the alias has inverted phase.
>
> To avoid this, the signal must be band-limited to below $f_s/2=5$ kHz *before* sampling, or $f_s$ must exceed $14$ kHz. A tone at 13 kHz would alias to $\lvert13-10\rvert=3$ kHz as well, and one at 17 kHz to $\lvert17-20\rvert=3$ kHz: every frequency $kf_s\pm3$ kHz lands on the same alias.

```desmos-graph
left=-0.3; right=6.3; bottom=-1.3; top=1.3
---
y=\cos(1.4\pi x)|blue
y=\cos(0.6\pi x)|red|dashed
(0,1)|black
(1,-0.309)|black
(2,-0.809)|black
(3,0.809)|black
(4,0.309)|black
(5,-1)|black
(6,0.309)|black
```

*Aliasing in time, with $t$ measured in sample periods: the 7 kHz tone (blue) and the 3 kHz tone (red dashed) pass through exactly the same samples (black) at $f_s=10$ kHz.*

```desmos-graph
left=-4; right=4; bottom=-0.3; top=1.6
---
y=\max(0,1-\left|x\right|)+\max(0,1-\left|x-1.5\right|)+\max(0,1-\left|x+1.5\right|)+\max(0,1-\left|x-3\right|)+\max(0,1-\left|x+3\right|)|blue
y=\max(0,1-\left|x\right|)|red|dashed
x=0.75|grey|dotted
x=-0.75|grey|dotted
```

*Undersampling ($\omega_M=1$, $\omega_s=1.5<2\omega_M$): the sum of overlapping replicas (blue) no longer matches the original triangle (red dashed) inside $\lvert\omega\rvert<\omega_s/2$ — aliasing.*

> [!warning] Common Pitfall — sampling exactly at twice the highest frequency
> Sampling at $f_s=2f_M$ is not enough for a sinusoid at $f_M$: $\sin(2\pi f_MnT)=\sin(\pi n)=0$ for all $n$, so the tone is invisible, and $\cos(\pi n+\theta)=\cos\theta\,(-1)^n$ cannot distinguish amplitude from phase. The inequality in the sampling theorem is strict. (For a signal whose spectrum contains no impulse at $\pm\omega_M$, equality is harmless, but nothing is gained by relying on it.)

### Non-ideal reconstruction: the zero-order hold

The interpolation formula needs the ideal, non-causal sinc kernel. The simplest causal alternative holds each sample constant for one interval. Mathematically this is the convolution of $x_s$ with a unit pulse on $[0,T)$:

$$
x_{\text{ZOH}}(t)=\sum_nx(nT)\operatorname{rect}\!\left(\frac{t-nT-T/2}{T}\right)\ \leftrightarrow\ X_s(\omega)\cdot T\operatorname{sinc}\!\left(\frac{\omega T}{2\pi}\right)e^{-j\omega T/2}.
$$

Compared with the ideal reconstruction filter (flat gain $T$ in the baseband), the hold introduces (i) a magnitude **droop** $\operatorname{sinc}(\omega T/2\pi)$, which at the folding frequency is $\operatorname{sinc}(\tfrac12)=2/\pi=0.637$ ($-3.9$ dB); (ii) a pure delay of $T/2$ (linear phase); and (iii) incomplete suppression of the replicas, attenuated only by the sinc sidelobes. Each is read directly from the convolution theorem.

> [!tip] Engineering Note — anti-aliasing and oversampling
> Real systems place an analog anti-aliasing low-pass filter before the sampler and oversample (for example $f_s=2.5f_M$, as in audio at 44.1 kHz for a 20 kHz band) to leave a guard band for a realisable filter roll-off. Choosing the filter order and the sampling rate is a design application of the replication formula $X_s=\frac1T\sum_kX(\omega-k\omega_s)$: the attenuation the filter provides at $f_s-f_M$ is the attenuation of the worst alias that can land in band.

## Summary

| Concept | Formula |
|---|---|
| Forward transform | $X(\omega)=\int_{-\infty}^{\infty}x(t)e^{-j\omega t}dt$ |
| Inverse transform | $x(t)=\frac{1}{2\pi}\int_{-\infty}^{\infty}X(\omega)e^{j\omega t}d\omega$ (midpoint at jumps) |
| Hz convention | $X_f(f)=X(2\pi f)$, $x(t)=\int X_f(f)e^{j2\pi ft}df$ |
| Relation to Laplace | $X(\omega)=X_L(j\omega)$ if the ROC contains the $j\omega$-axis |
| DC values | $X(0)=\int x\,dt$, $x(0)=\frac{1}{2\pi}\int X\,d\omega$ |
| Exponential | $e^{-at}u(t)\leftrightarrow1/(a+j\omega)$; $e^{-a\lvert t\rvert}\leftrightarrow2a/(a^2+\omega^2)$ |
| Rectangle / sinc | $\operatorname{rect}(t/T)\leftrightarrow T\operatorname{sinc}(\omega T/2\pi)$; $\sin(Wt)/(\pi t)\leftrightarrow\operatorname{rect}(\omega/2W)$ |
| Triangle | $\Lambda(t/T)\leftrightarrow T\operatorname{sinc}^2(\omega T/2\pi)$ |
| Gaussian | $e^{-t^2/2\sigma^2}\leftrightarrow\sigma\sqrt{2\pi}e^{-\sigma^2\omega^2/2}$ |
| Generalized pairs | $\delta\leftrightarrow1$, $1\leftrightarrow2\pi\delta(\omega)$, $u\leftrightarrow\pi\delta(\omega)+1/(j\omega)$, $\operatorname{sgn}\leftrightarrow2/(j\omega)$ |
| Sinusoids | $\cos\omega_0t\leftrightarrow\pi[\delta(\omega-\omega_0)+\delta(\omega+\omega_0)]$ |
| Time shift / frequency shift | $x(t-t_0)\leftrightarrow e^{-j\omega t_0}X$; $e^{j\omega_0t}x\leftrightarrow X(\omega-\omega_0)$ |
| Modulation | $x\cos\omega_0t\leftrightarrow\frac12[X(\omega-\omega_0)+X(\omega+\omega_0)]$ |
| Scaling | $x(at)\leftrightarrow X(\omega/a)/\lvert a\rvert$ |
| Duality | $X(t)\leftrightarrow2\pi x(-\omega)$ |
| Differentiation / multiplication by $t$ | $\dot x\leftrightarrow j\omega X$; $tx\leftrightarrow j\,dX/d\omega$ |
| Integration | $\int_{-\infty}^tx\,d\tau\leftrightarrow X/(j\omega)+\pi X(0)\delta(\omega)$ |
| Convolution / multiplication | $x*h\leftrightarrow XH$; $xy\leftrightarrow\frac1{2\pi}X*Y$ |
| Real-signal symmetry | $X(-\omega)=X^*(\omega)$; $\lvert X\rvert$ even, $\angle X$ odd |
| Frequency response | $e^{j\omega t}\to H(\omega)e^{j\omega t}$, $H=\mathcal{F}\{h\}$ |
| Distortionless system / group delay | $H=Ke^{-j\omega t_d}$; $\tau_g=-d\angle H/d\omega$ |
| Uncertainty principle | $\Delta t\,\Delta\omega\ge\frac12$, equality for the Gaussian |
| Periodic signal | $X=2\pi\sum_kc_k\delta(\omega-k\omega_0)$, $c_k=X_1(k\omega_0)/T_0$ |
| Windowing | $X_w=\frac{1}{2\pi}X*W$; resolution $\approx1/T$ Hz |
| Hann window | $W_H=\frac T2\operatorname{sinc}(u)/(1-u^2)$, $u=\omega T/2\pi$ |
| Parseval | $\int\lvert x\rvert^2dt=\frac1{2\pi}\int\lvert X\rvert^2d\omega$; $\int xy^*dt=\frac1{2\pi}\int XY^*d\omega$ |
| ESD | $\Psi_x=\lvert X\rvert^2\leftrightarrow R_x(\tau)$; $\Psi_y=\lvert H\rvert^2\Psi_x$ |
| Band energy (real $x$) | $E_{[\omega_1,\omega_2]}=\frac1\pi\int_{\omega_1}^{\omega_2}\lvert X\rvert^2d\omega$ |
| Impulse train | $\sum_n\delta(t-nT)\leftrightarrow\omega_s\sum_k\delta(\omega-k\omega_s)$ |
| Sampled spectrum | $X_s(\omega)=\frac1T\sum_kX(\omega-k\omega_s)=\sum_nx(nT)e^{-j\omega nT}$ |
| Nyquist condition | $\omega_s>2\omega_M$ ($f_s>2f_M$) |
| Interpolation | $x(t)=\sum_nx(nT)\operatorname{sinc}((t-nT)/T)$ |
| Alias frequency | $f_a=\lvert f_0-kf_s\rvert$, $k$ nearest integer to $f_0/f_s$ |

## Further reading

- A. V. Oppenheim, A. S. Willsky & S. H. Nawab, *Signals and Systems*, 2nd ed., Prentice Hall — Chs. 4 (continuous-time Fourier transform), 6 (time and frequency characterisation) and 7 (sampling); uses the angular-frequency convention of this chapter.
- R. N. Bracewell, *The Fourier Transform and Its Applications*, 3rd ed., McGraw-Hill — the classic treatment in the Hz convention, with extensive pictorial dictionaries of pairs and a thorough account of the impulse train and sampling.
- B. P. Lathi & R. Green, *Signal Processing and Linear Systems*, 2nd ed., Oxford University Press — energy spectral density, essential bandwidth and filtering, with many circuit examples.
- A. Papoulis, *The Fourier Integral and Its Applications*, McGraw-Hill — rigorous yet engineering-oriented treatment of existence, generalized functions, windows and the uncertainty principle.
