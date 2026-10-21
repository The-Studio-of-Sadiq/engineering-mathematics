# Chapter 15 — Fourier Series and Fourier Analysis

[← Previous](14-laplace-transforms.md) · [Table of Contents](../00-front/table-of-contents.md) · [Next →](16-fourier-transforms.md)

> **Prerequisites:** Ch. 2 (complex numbers), Ch. 3 (trigonometry and periodic functions), Ch. 5 (integral calculus), Ch. 9 (vector spaces and inner products), Ch. 12 (linear differential equations)
> **Learning outcomes:** After this chapter you can
> - compute trigonometric and complex (exponential) Fourier series of periodic waveforms of any period, and convert between the coefficient forms $(a_n,b_n)$, $(A_n,\phi_n)$ and $c_n$;
> - use even, odd, half-wave and quarter-wave symmetry to predict which harmonics are present before integrating, and build half-range expansions;
> - state the Dirichlet conditions and the pointwise, uniform and mean-square convergence theorems, and explain and quantify the Gibbs phenomenon;
> - obtain coefficients quickly with the shift, differentiation, integration and jump methods, and predict how fast they decay from the smoothness of the waveform;
> - apply Parseval's theorem to compute RMS values, harmonic powers, power factor and total harmonic distortion, and to sum numerical series;
> - find the steady-state response of a linear time-invariant system to any periodic input, harmonic by harmonic.

A periodic waveform that is not a sinusoid — a clock square wave, the triangular ripple current of a buck converter, the humps of a rectified sine, the stepped output of an inverter — is still, in a precise sense, *made of* sinusoids. [Chapter 3](../part-1-foundations/03-trigonometry-and-periodic-functions.md) introduced this idea informally: it established the orthogonality of harmonics (§3.7), the power of a finite sum of harmonics, the notion of total harmonic distortion, and the shapes and symmetries of the standard waveforms (§3.8). This chapter turns those observations into a complete theory. We ask, and answer, four questions:

1. **Representation.** Given a periodic function, which sinusoids does it contain, and with what amplitudes and phases? (The Fourier coefficients, §15.2–15.3, §15.7.)
2. **Structure.** How do properties of the waveform — symmetry, smoothness, discontinuities — show up in its coefficients? (§15.4, §15.5, §15.7.)
3. **Validity.** In what sense does the infinite sum of sinusoids actually equal the function, and what happens at jumps? (§15.6.)
4. **Use.** How do power, RMS values and distortion follow from the coefficients, and how does a linear system respond to a periodic input? (§15.8, §15.9.)

The [next chapter](16-fourier-transforms.md) lets the period tend to infinity and turns the discrete set of harmonics into the continuous Fourier transform. Everything here — orthogonality, Parseval, convolution, the trade-off between smoothness and spectral spread — carries over almost word for word.

**Notation.** A periodic function $x(t)$ has fundamental period $T_0$, fundamental frequency $f_0=1/T_0$ and fundamental angular frequency $\omega_0=2\pi/T_0$ (Ch. 3, §3.4). The harmonic index is $n$, so the $n$th harmonic has angular frequency $n\omega_0$. (Chapter 16 uses $k$ for the same index when the discrete-time index $n$ is also in play; the formulas are identical.) Many examples use the normalised period $T_0=2\pi$, $\omega_0=1$, in which case the time variable is written $t$ and harmonics are $\cos nt$, $\sin nt$. The notation $\int_{T_0}$ means an integral over *any* interval of length $T_0$.

## 15.1 Periodic function decomposition

### Motivation: why sinusoids, and why decompose at all?

There are infinitely many ways to write a function as a sum of simpler functions — powers of $t$ (Taylor series, Ch. 4), polynomials, piecewise constants, wavelets. Engineers single out sums of sinusoids for one overriding reason:

> **A complex exponential $e^{j\omega t}$ is an eigenfunction of every linear time-invariant (LTI) system.**

If an LTI system with impulse response $h(t)$ is driven by $x(t)=e^{j\omega t}$, its output is

$$
y(t)=\int_{-\infty}^{\infty}h(\tau)\,e^{j\omega(t-\tau)}\,d\tau=\Big[\int_{-\infty}^{\infty}h(\tau)e^{-j\omega\tau}\,d\tau\Big]e^{j\omega t}=H(j\omega)\,e^{j\omega t},
$$

provided the integral converges. The input passes through *unchanged in shape*, multiplied by a single complex number $H(j\omega)$ — the frequency response, which for a circuit is the phasor transfer function of Ch. 2 and the Laplace transfer function $H(s)$ of [Ch. 14](14-laplace-transforms.md) evaluated at $s=j\omega$. No other family of signals has this property for *all* LTI systems at once.

Now suppose a periodic input can be written as a sum of exponentials at the harmonic frequencies,

$$
x(t)=\sum_{n=-\infty}^{\infty}c_n\,e^{jn\omega_0t}.
$$

By linearity (and, for an infinite sum, a convergence argument we postpone to §15.9) the steady-state output is

$$
y(t)=\sum_{n=-\infty}^{\infty}c_n\,H(jn\omega_0)\,e^{jn\omega_0t}.
$$

A problem that in the time domain requires solving a differential equation with a discontinuous forcing function becomes, in the harmonic domain, a list of multiplications. That is the whole engineering case for Fourier series: **decompose, multiply harmonic by harmonic, recombine**. Section 15.9 carries this programme out in detail; the intervening sections build the mathematics that makes it trustworthy.

A second, independent motivation is *measurement*. Power-quality standards, audio distortion specifications, EMC emission limits and spectrum-analyser displays are all stated harmonic by harmonic. To read or write such a specification you must be able to move fluently between a waveform and its list of harmonics.

### Periodic functions and their building blocks

> [!definition] Periodic function
> A function $x:\mathbb{R}\to\mathbb{C}$ is **periodic with period $T>0$** if $x(t+T)=x(t)$ for all $t$. The smallest such $T$, if it exists, is the **fundamental period** $T_0$. A function defined only on an interval of length $T_0$ is extended to all of $\mathbb{R}$ by **periodic extension**, $x(t+kT_0)=x(t)$, $k\in\mathbb{Z}$.

Every function in the set

$$
\mathcal{B}=\{\,1,\ \cos\omega_0t,\ \sin\omega_0t,\ \cos2\omega_0t,\ \sin2\omega_0t,\ \dots\,\}
$$

has period $T_0$ (the $n$th harmonic has period $T_0/n$, hence also period $T_0$). So does every finite linear combination of them, a **trigonometric polynomial of degree $N$**:

$$
S_N(t)=\frac{a_0}{2}+\sum_{n=1}^{N}\big(a_n\cos n\omega_0t+b_n\sin n\omega_0t\big).
$$

The decomposition problem is: *given a $T_0$-periodic $x(t)$, choose the coefficients so that $S_N$ is as close as possible to $x$, and decide whether $S_N\to x$ as $N\to\infty$.*

### The inner-product viewpoint

Chapter 9 treated vectors abstractly as elements of a space with an inner product. Functions of period $T_0$ form such a space. For real or complex functions that are square-integrable over a period, define

$$
\langle x,y\rangle=\frac{1}{T_0}\int_{T_0}x(t)\,y^*(t)\,dt,\qquad \lVert x\rVert=\sqrt{\langle x,x\rangle}=\sqrt{\frac1{T_0}\int_{T_0}\lvert x(t)\rvert^2dt}.
$$

The normalisation $1/T_0$ is the engineer's choice: it makes $\lVert x\rVert$ the **RMS value** of $x$ and $\lVert x\rVert^2$ its **mean-square value**, the average power it would dissipate in a 1 Ω resistor. (Mathematics texts often omit $1/T_0$, which changes constants but nothing else.) The space of all such functions is denoted $L^2(T_0)$.

> [!theorem] Orthogonality of the harmonic basis
> With $\omega_0=2\pi/T_0$ and integers $m,n$:
> $$
> \langle e^{jm\omega_0t},e^{jn\omega_0t}\rangle=\frac1{T_0}\int_{T_0}e^{j(m-n)\omega_0t}\,dt=\begin{cases}1,&m=n,\\0,&m\ne n.\end{cases}
> $$
> Equivalently, for the real functions with $m,n\ge0$:
> $$
> \begin{aligned}
> \frac{2}{T_0}\int_{T_0}\cos m\omega_0t\,\cos n\omega_0t\,dt&=\delta_{mn}\quad(m,n\ge1), & \frac{2}{T_0}\int_{T_0}\sin m\omega_0t\,\sin n\omega_0t\,dt&=\delta_{mn},\\
> \int_{T_0}\sin m\omega_0t\,\cos n\omega_0t\,dt&=0, & \frac1{T_0}\int_{T_0}1\cdot 1\,dt&=1,
> \end{aligned}
> $$
> where $\delta_{mn}$ is the Kronecker delta.

**Proof (complex form).** For $m=n$ the integrand is 1. For $m\ne n$, let $k=m-n\ne0$; then $\int_{t_1}^{t_1+T_0}e^{jk\omega_0t}dt=\big[e^{jk\omega_0t}/(jk\omega_0)\big]_{t_1}^{t_1+T_0}=e^{jk\omega_0t_1}(e^{j2\pi k}-1)/(jk\omega_0)=0$. The real relations follow by writing $\cos$ and $\sin$ as combinations of $e^{\pm j\theta}$, or directly by product-to-sum as in Ch. 3, §3.7. $\square$

So the set $\{e^{jn\omega_0t}\}_{n\in\mathbb{Z}}$ is **orthonormal**, and the real set $\{1,\sqrt2\cos n\omega_0t,\sqrt2\sin n\omega_0t\}$ is orthonormal too. Each harmonic is a unit vector at "right angles" to every other.

### Coefficients as projections

In $\mathbb{R}^3$ the component of a vector $\mathbf{v}$ along a unit vector $\hat{\mathbf{e}}$ is the projection $\langle\mathbf{v},\hat{\mathbf{e}}\rangle$, and the closest point to $\mathbf{v}$ in a plane spanned by orthonormal $\hat{\mathbf{e}}_1,\hat{\mathbf{e}}_2$ is $\langle\mathbf{v},\hat{\mathbf{e}}_1\rangle\hat{\mathbf{e}}_1+\langle\mathbf{v},\hat{\mathbf{e}}_2\rangle\hat{\mathbf{e}}_2$. Exactly the same holds in the function space.

> [!theorem] Best approximation (least squares)
> Let $x\in L^2(T_0)$ and let $\{\varphi_n\}_{n=-N}^{N}$ with $\varphi_n(t)=e^{jn\omega_0t}$. Among all trigonometric polynomials $p(t)=\sum_{n=-N}^{N}d_n\varphi_n(t)$, the mean-square error $\lVert x-p\rVert^2$ is minimised uniquely by
> $$
> d_n=c_n=\langle x,\varphi_n\rangle=\frac1{T_0}\int_{T_0}x(t)\,e^{-jn\omega_0t}\,dt,
> $$
> and the minimum error is
> $$
> \lVert x-S_N\rVert^2=\lVert x\rVert^2-\sum_{n=-N}^{N}\lvert c_n\rvert^2 .
> $$

**Proof.** Expand, using orthonormality:

$$
\begin{aligned}
\lVert x-p\rVert^2&=\lVert x\rVert^2-\sum_n\big(d_n^*\langle x,\varphi_n\rangle+d_n\langle x,\varphi_n\rangle^*\big)+\sum_n\lvert d_n\rvert^2\\
&=\lVert x\rVert^2-\sum_n\lvert c_n\rvert^2+\sum_n\lvert d_n-c_n\rvert^2 .
\end{aligned}
$$

(The second line follows by completing the square: $\lvert d_n-c_n\rvert^2=\lvert d_n\rvert^2-d_n^*c_n-d_nc_n^*+\lvert c_n\rvert^2$.) The last sum is non-negative and vanishes only when every $d_n=c_n$. $\square$

Three consequences are worth stating separately, because they recur throughout the chapter.

1. **Coefficients do not change when $N$ grows.** Because the basis is orthogonal, the best coefficient for $\cos3\omega_0t$ is the same whether we approximate with 3 harmonics or 300. Contrast this with fitting a polynomial by least squares in the monomials $1,t,t^2,\dots$, where adding a term changes all the earlier coefficients. This *finality of coefficients* is why a spectrum analyser can report "the 3rd harmonic" without reference to how many others it measures.
2. **Bessel's inequality.** Since the error is non-negative, $\sum_{n=-N}^{N}\lvert c_n\rvert^2\le\lVert x\rVert^2$ for every $N$, hence $\sum_{n=-\infty}^{\infty}\lvert c_n\rvert^2\le\lVert x\rVert^2<\infty$. In particular $c_n\to0$ as $\lvert n\rvert\to\infty$.
3. **Completeness and Parseval.** If the error tends to zero as $N\to\infty$, Bessel's inequality becomes an equality — Parseval's theorem (§15.8). That it does tend to zero for *every* square-integrable $x$ is a deep theorem (the completeness of the trigonometric system) whose proof sketch appears in §15.6.

> [!definition] Fourier series
> The **Fourier series** of a $T_0$-periodic function $x$ is the formal series whose coefficients are the projections of $x$ on the harmonics:
> $$
> x(t)\sim\sum_{n=-\infty}^{\infty}c_n\,e^{jn\omega_0t}=\frac{a_0}{2}+\sum_{n=1}^{\infty}\big(a_n\cos n\omega_0t+b_n\sin n\omega_0t\big).
> $$
> The symbol $\sim$ means "has the Fourier series"; replacing it by $=$ requires the convergence theorems of §15.6.

> [!example] Worked Example 15.1 — Least-squares approximation of a sawtooth
> Let $x(t)=t$ on $(-\pi,\pi)$, extended with period $2\pi$ (a sawtooth). Find the best mean-square approximation by (a) a single term $b_1\sin t$ and (b) $\sum_{n=1}^{3}b_n\sin nt$, and the residual mean-square error in each case. (Cosine terms and the constant are zero because $x$ is odd; see §15.4.)
>
> **Solution.** Mean-square value: $\lVert x\rVert^2=\frac1{2\pi}\int_{-\pi}^{\pi}t^2dt=\frac{\pi^2}{3}=3.2899$.
>
> The projection on $\sin nt$ in real form is $b_n=\frac{2}{2\pi}\int_{-\pi}^{\pi}t\sin nt\,dt$. Integrating by parts,
> $$
> \int_{-\pi}^{\pi}t\sin nt\,dt=\Big[-\frac{t\cos nt}{n}\Big]_{-\pi}^{\pi}+\frac1n\int_{-\pi}^{\pi}\cos nt\,dt=-\frac{2\pi\cos n\pi}{n}=\frac{2\pi(-1)^{n+1}}{n},
> $$
> so $b_n=2(-1)^{n+1}/n$: $b_1=2$, $b_2=-1$, $b_3=2/3$.
>
> Each term $b_n\sin nt$ has mean-square value $b_n^2/2$, and by orthogonality the error is $\lVert x\rVert^2$ minus the sum of these.
>
> (a) Error $=3.2899-\tfrac{4}{2}=1.2899$: the single sine captures $2/3.2899=60.8\%$ of the mean-square value.
>
> (b) Error $=3.2899-2-0.5-0.2222=0.5676$: three terms capture $82.7\%$.
>
> **Check.** Numerically integrating $(t-\sum_{n\le3}b_n\sin nt)^2$ over a period and dividing by $2\pi$ gives $0.56765$. ✓ Note that $b_1=2$ is the best coefficient in both (a) and (b) — finality of coefficients. The slow capture (still 17% error with three terms) reflects the jump of the sawtooth at $t=\pm\pi$; §15.7 explains why jumps force $\lvert b_n\rvert\propto1/n$.

> [!example] Worked Example 15.2 — Why a non-orthogonal basis is awkward
> Approximate $x(t)=\lvert t\rvert$ on $(-\pi,\pi)$ in mean square by (a) a constant $d_0$, and (b) $d_0+d_1\cos t$. Compare with what happens if we instead used the non-orthogonal pair $\{1,\cos^2(t/2)\}$.
>
> **Solution.** (a) The best constant is the projection on $1$: $d_0=\frac{1}{2\pi}\int_{-\pi}^{\pi}\lvert t\rvert dt=\frac{\pi}{2}=1.5708$, the mean value.
>
> (b) Because $1$ and $\cos t$ are orthogonal, $d_0$ is unchanged and $d_1=\frac{2}{2\pi}\int_{-\pi}^{\pi}\lvert t\rvert\cos t\,dt=\frac{2}{\pi}\int_0^{\pi}t\cos t\,dt=\frac{2}{\pi}[t\sin t+\cos t]_0^\pi=\frac{2}{\pi}(-2)=-\frac{4}{\pi}=-1.2732$.
>
> Now note $\cos^2(t/2)=\tfrac12+\tfrac12\cos t$ spans the same two-dimensional space as $\{1,\cos t\}$, so the best approximation is the same function $1.5708-1.2732\cos t=e_0+e_1\cos^2(t/2)$ with $e_1=-2.5465$ and $e_0=1.5708+1.2732=2.8440$. The best constant *by itself* was $1.5708$, but in the combination it must change to $2.8440$: with a non-orthogonal basis, adding a function alters the coefficients already found, and each coefficient is no longer a single integral but the solution of a linear system (the normal equations of Ch. 8). Orthogonality makes each coefficient an independent, final, one-line integral.

> [!tip] Engineering note — a spectrum analyser computes projections
> A swept or FFT spectrum analyser estimates $\lvert c_n\rvert$ by correlating the input with $e^{-jn\omega_0t}$ over a finite window: exactly the projection integral above. A lock-in amplifier does the same for one chosen harmonic, multiplying by a reference $\cos n\omega_0t$ and $\sin n\omega_0t$ and low-pass filtering (averaging) the products. The orthogonality relations are why the other harmonics average out of the measurement.

> [!warning] Pitfall — the decomposition needs exact periodicity
> Orthogonality holds only over an *exact* period, or a whole number of periods. Integrating over a window that is not a whole number of periods (a common situation when analysing measured data) leaks energy from each harmonic into its neighbours. The Fourier series is a statement about an exactly periodic idealisation; finite-record effects belong to the Fourier transform and the DFT (Ch. 16, Ch. 24).

## 15.2 Trigonometric Fourier series

### Motivation

The trigonometric form uses only real functions — cosines and sines — and so is the natural form when the waveform is real, when we want to see its symmetry directly, and when we want amplitudes and phases that can be read on an oscilloscope or written into a phasor calculation. It is also the form that arises from separation of variables in the heat, wave and Laplace equations ([Ch. 13](../part-4-differential-equations/13-partial-differential-equations.md)), where sines and cosines appear as eigenfunctions of $d^2/dx^2$ with boundary conditions.

### The Euler (Euler–Fourier) formulas

Suppose a $T_0$-periodic function can be written as

$$
x(t)=\frac{a_0}{2}+\sum_{n=1}^{\infty}\big(a_n\cos n\omega_0t+b_n\sin n\omega_0t\big),
$$

and that the series may be integrated term by term after multiplication by a bounded function (true, for example, when it converges uniformly; §15.6 shows the formulas below are correct in far greater generality). Multiply both sides by $\cos m\omega_0t$, $m\ge1$, and integrate over a period. By orthogonality every term vanishes except $a_m\int_{T_0}\cos^2m\omega_0t\,dt=a_mT_0/2$. Similarly with $\sin m\omega_0t$. Integrating without a multiplier kills every harmonic and leaves $\frac{a_0}{2}T_0$.

> [!theorem] Euler formulas for the trigonometric Fourier coefficients
> For a $T_0$-periodic function integrable over a period, with $\omega_0=2\pi/T_0$,
> $$
> a_n=\frac{2}{T_0}\int_{T_0}x(t)\cos n\omega_0t\,dt\quad(n\ge0),\qquad b_n=\frac{2}{T_0}\int_{T_0}x(t)\sin n\omega_0t\,dt\quad(n\ge1).
> $$
> The constant term is $\dfrac{a_0}{2}=\dfrac1{T_0}\displaystyle\int_{T_0}x(t)\,dt$, the **mean (DC) value** of $x$.

The convention of writing the constant as $a_0/2$ is purely so that the formula for $a_n$ also gives $a_0$; some books write $a_0$ for the mean value instead and adjust the formula. Always check which convention a table uses.

For the normalised period $T_0=2\pi$ ($\omega_0=1$) the formulas read

$$
a_n=\frac1\pi\int_{-\pi}^{\pi}x(t)\cos nt\,dt,\qquad b_n=\frac1\pi\int_{-\pi}^{\pi}x(t)\sin nt\,dt .
$$

For a function given on $(-L,L)$, so $T_0=2L$ and $\omega_0=\pi/L$ — the usual form for spatial problems in Ch. 13:

$$
x(s)\sim\frac{a_0}2+\sum_{n=1}^\infty\Big(a_n\cos\frac{n\pi s}{L}+b_n\sin\frac{n\pi s}{L}\Big),\qquad a_n=\frac1L\int_{-L}^{L}x(s)\cos\frac{n\pi s}{L}\,ds,\quad b_n=\frac1L\int_{-L}^{L}x(s)\sin\frac{n\pi s}{L}\,ds.
$$

**Choice of interval.** The integrand $x(t)\cos n\omega_0t$ is itself $T_0$-periodic, and the integral of a periodic function over one period does not depend on where the period starts:

$$
\int_{\alpha}^{\alpha+T_0}g(t)\,dt=\int_{0}^{T_0}g(t)\,dt\qquad\text{for every }\alpha,
$$

because the piece cut from one end is exactly the piece added at the other. Choose the interval that makes the integration easiest — usually the one in which $x$ has the fewest pieces or a symmetric form.

**Dimensions.** $a_n$ and $b_n$ carry the units of $x$ (the factor $2/T_0$ cancels the seconds from $dt$), and $n\omega_0t$ is dimensionless. A coefficient in volt-seconds signals a forgotten $1/T_0$.

### Amplitude–phase (harmonic) form

Grouping each cosine–sine pair into one sinusoid (Ch. 3, §3.2):

$$
a_n\cos n\omega_0t+b_n\sin n\omega_0t=A_n\cos(n\omega_0t+\phi_n),\qquad A_n=\sqrt{a_n^2+b_n^2},\quad \phi_n=\operatorname{atan2}(-b_n,\,a_n).
$$

So

$$
x(t)\sim A_0+\sum_{n=1}^{\infty}A_n\cos(n\omega_0t+\phi_n),\qquad A_0=\frac{a_0}{2}.
$$

The pair $(A_n,\phi_n)$ is the **amplitude** and **phase** of the $n$th harmonic; plotted against $n\omega_0$ they form the one-sided amplitude and phase spectra (§15.5). Note the sign in $\phi_n$: with the cosine reference, a pure $+\sin$ term has phase $-90^\circ$. A sine-reference form $A_n\sin(n\omega_0t+\psi_n)$ with $\psi_n=\phi_n+90^\circ$ is common in power-systems literature; state the reference whenever phases are quoted.

### A first full example: the square wave

> [!example] Worked Example 15.3 — Square wave
> Find the Fourier series of the square wave $x(t)=A$ for $0<t<T_0/2$, $x(t)=-A$ for $T_0/2<t<T_0$, period $T_0$.
>
> **Solution.** Mean value: the two halves cancel, so $a_0=0$.
>
> Cosine coefficients:
> $$
> a_n=\frac2{T_0}\Big[\int_0^{T_0/2}A\cos n\omega_0t\,dt-\int_{T_0/2}^{T_0}A\cos n\omega_0t\,dt\Big]=\frac{2A}{T_0n\omega_0}\Big[\sin n\omega_0t\Big|_0^{T_0/2}-\sin n\omega_0t\Big|_{T_0/2}^{T_0}\Big]=0,
> $$
> since $\sin(n\pi)=\sin(2n\pi)=0$.
>
> Sine coefficients, with $n\omega_0T_0/2=n\pi$ and $T_0\omega_0=2\pi$:
> $$
> b_n=\frac{2A}{T_0n\omega_0}\Big[-\cos n\omega_0t\Big|_0^{T_0/2}+\cos n\omega_0t\Big|_{T_0/2}^{T_0}\Big]=\frac{A}{n\pi}\big[(1-\cos n\pi)+(1-\cos n\pi)\big]=\frac{2A}{n\pi}\big(1-(-1)^n\big).
> $$
> So $b_n=4A/(n\pi)$ for odd $n$ and $0$ for even $n$:
> $$
> x(t)=\frac{4A}{\pi}\Big(\sin\omega_0t+\frac13\sin3\omega_0t+\frac15\sin5\omega_0t+\cdots\Big)=\frac{4A}{\pi}\sum_{k=1}^{\infty}\frac{\sin(2k-1)\omega_0t}{2k-1}.
> $$
>
> **Checks.** (i) Only odd harmonics — predicted by half-wave symmetry (§15.4). (ii) Only sines — predicted by odd symmetry. (iii) At $t=T_0/4$: $\frac{4A}{\pi}(1-\frac13+\frac15-\cdots)=\frac{4A}{\pi}\cdot\frac\pi4=A$ by the Leibniz series. ✓ (iv) The fundamental amplitude $4A/\pi=1.273A$ *exceeds* the square wave's peak — a reminder that harmonics can add destructively.

```desmos-graph
left=-1; right=13; bottom=-1.6; top=1.6
---
y=\operatorname{sign}(\sin(x))|black|dashed
y=\frac{4}{\pi}\sin(x)|red
y=\sum_{k=1}^{2}\frac{4}{\pi(2k-1)}\sin((2k-1)x)|green
y=\sum_{k=1}^{10}\frac{4}{\pi(2k-1)}\sin((2k-1)x)|blue
```
*Partial sums of the square-wave series with $\omega_0=1$, $A=1$: fundamental only (red), harmonics 1 and 3 (green), harmonics up to 19 (blue). The sum approaches the square wave (dashed) everywhere except near the jumps, where an overshoot persists (§15.6).*

### More standard waveforms

> [!example] Worked Example 15.4 — Ramp on a non-symmetric interval
> A sawtooth generator produces $x(t)=t$ volts for $0<t<2$ s, repeating with period $T_0=2$ s. Find its trigonometric Fourier series.
>
> **Solution.** $\omega_0=2\pi/2=\pi$ rad/s. The interval $(0,2)$ is the natural one here.
>
> Mean value: $\frac{a_0}2=\frac12\int_0^2t\,dt=1$ V.
>
> For $n\ge1$ integrate by parts, using $\sin2n\pi=0$, $\cos2n\pi=1$:
> $$
> a_n=\frac22\int_0^2t\cos n\pi t\,dt=\Big[\frac{t\sin n\pi t}{n\pi}\Big]_0^2-\frac{1}{n\pi}\int_0^2\sin n\pi t\,dt=0+\frac{1}{n^2\pi^2}\big[\cos n\pi t\big]_0^2=0,
> $$
> $$
> b_n=\int_0^2t\sin n\pi t\,dt=\Big[-\frac{t\cos n\pi t}{n\pi}\Big]_0^2+\frac{1}{n\pi}\int_0^2\cos n\pi t\,dt=-\frac{2}{n\pi}.
> $$
> Hence
> $$
> x(t)=1-\frac{2}{\pi}\sum_{n=1}^\infty\frac{\sin n\pi t}{n}\ \text{V}.
> $$
>
> **Interpretation.** After removing the mean, the periodic function $x(t)-1$ is odd: $x(-t)-1=-(x(t)-1)$ (check: for $0<t<2$, $x(-t)=x(2-t)=2-t$, so $x(-t)-1=1-t$). That is why only sines survive. Every harmonic, odd and even, is present, with amplitude $2/(n\pi)$: 0.637 V, 0.318 V, 0.212 V, … — falling only as $1/n$ because of the jump at $t=0$.
>
> **Check.** Numerical quadrature gives $a_1=0.0000$, $b_1=-0.6366$, $b_2=-0.3183$, $b_3=-0.2122$. ✓

> [!example] Worked Example 15.5 — Half-wave rectified sine
> A half-wave rectifier fed by $V_m\sin\omega_0t$ produces $x(t)=V_m\sin\omega_0t$ for $0\le\omega_0t\le\pi$ and $x(t)=0$ for $\pi<\omega_0t<2\pi$. Find the Fourier series.
>
> **Solution.** Work in $\theta=\omega_0t$, period $2\pi$, and set $V_m=1$ (multiply by $V_m$ at the end).
>
> Mean: $\frac{a_0}2=\frac1{2\pi}\int_0^\pi\sin\theta\,d\theta=\frac1\pi$.
>
> Cosine coefficients. Use $\sin\theta\cos n\theta=\frac12[\sin(1+n)\theta+\sin(1-n)\theta]$. For $n\ne1$,
> $$
> a_n=\frac1\pi\int_0^\pi\sin\theta\cos n\theta\,d\theta=\frac1{2\pi}\Big[\frac{1-\cos(1+n)\pi}{1+n}+\frac{1-\cos(1-n)\pi}{1-n}\Big].
> $$
> For odd $n\ne1$, $1\pm n$ is even, both cosines equal 1 and $a_n=0$. For even $n$, $\cos(1\pm n)\pi=-1$ and
> $$
> a_n=\frac1{2\pi}\Big[\frac2{1+n}+\frac2{1-n}\Big]=\frac1\pi\cdot\frac{2}{1-n^2}=-\frac{2}{\pi(n^2-1)}.
> $$
> For $n=1$: $a_1=\frac1\pi\int_0^\pi\sin\theta\cos\theta\,d\theta=0$.
>
> Sine coefficients. $b_1=\frac1\pi\int_0^\pi\sin^2\theta\,d\theta=\frac12$. For $n\ge2$, $\sin\theta\sin n\theta=\frac12[\cos(n-1)\theta-\cos(n+1)\theta]$ integrates over $(0,\pi)$ to $\frac12\big[\frac{\sin(n-1)\pi}{n-1}-\frac{\sin(n+1)\pi}{n+1}\big]=0$. So $b_n=0$ for $n\ge2$.
>
> Result (writing $n=2k$):
> $$
> x(t)=V_m\Big[\frac1\pi+\frac12\sin\omega_0t-\frac2\pi\sum_{k=1}^\infty\frac{\cos2k\omega_0t}{4k^2-1}\Big]=V_m\Big[0.3183+0.5\sin\omega_0t-0.2122\cos2\omega_0t-0.0424\cos4\omega_0t-\cdots\Big].
> $$
>
> **Interpretation.** The waveform contains DC (the useful output of a rectifier), *half* of the original fundamental, and even harmonics only. The structure is explained in §15.4: $x(t)=\frac{V_m}{2}\big(\sin\omega_0t+\lvert\sin\omega_0t\rvert\big)$, the sum of half the input sinusoid and half of a full-wave rectified sine, which contains only DC and even harmonics.
>
> **Check.** Numerical quadrature: $a_0=0.63662$, $b_1=0.5$, $a_2=-0.212207$, $a_4=-0.042441$, $a_3=b_3=0$. ✓

> [!example] Worked Example 15.6 — Rectangular pulse train (trigonometric form)
> A pulse train has height $A$, width $\tau$ and period $T_0$, centred on $t=0$: $x(t)=A$ for $\lvert t\rvert<\tau/2$ and $0$ elsewhere in $(-T_0/2,T_0/2)$. Let $D=\tau/T_0$ be the duty cycle. Find $a_n,b_n$.
>
> **Solution.** The function is even about $t=0$, so every $b_n=\frac{2}{T_0}\int x(t)\sin n\omega_0t\,dt=0$ (odd integrand on a symmetric interval). Mean value $a_0/2=A\tau/T_0=AD$. For $n\ge1$,
> $$
> a_n=\frac2{T_0}\int_{-\tau/2}^{\tau/2}A\cos n\omega_0t\,dt=\frac{4A}{T_0n\omega_0}\sin\frac{n\omega_0\tau}{2}=\frac{2A}{n\pi}\sin(n\pi D)=2AD\operatorname{sinc}(nD),
> $$
> with the normalised $\operatorname{sinc}(u)=\sin(\pi u)/(\pi u)$ used in Ch. 16. So
> $$
> x(t)=AD+2AD\sum_{n=1}^{\infty}\operatorname{sinc}(nD)\cos n\omega_0t .
> $$
>
> **Special cases.** $D=\tfrac12$ gives $a_n=\frac{2A}{n\pi}\sin\frac{n\pi}2$: $\frac{2A}{\pi},0,-\frac{2A}{3\pi},0,\frac{2A}{5\pi},\dots$ — a square wave of levels $0$ and $A$, i.e. $\frac A2$ plus half of Worked Example 15.3 shifted by a quarter period (the sign alternation is that shift; §15.7). $D=\tfrac14$ gives harmonics 4, 8, 12, … missing, because $\sin(n\pi/4)=0$ at those $n$. In general, whenever $1/D$ is an integer, the harmonics at multiples of $1/D$ vanish.
>
> The envelope of $a_n$ is $2AD\operatorname{sinc}(fD/f_0)$ with first null at $f=1/\tau$: **narrower pulses spread their power over more harmonics.** This is the discrete version of the pulse-width/bandwidth reciprocity of Ch. 16.

> [!tip] Engineering note — table conventions
> Tables of Fourier series differ in (i) whether the constant is $a_0$ or $a_0/2$, (ii) the time origin of the waveform (a shift changes the sine/cosine mix but not the amplitudes $A_n$), (iii) whether amplitudes are peak or RMS values, and (iv) whether frequency is in rad/s or Hz. When a datasheet says "3rd harmonic 4%" it almost always means the ratio $A_3/A_1$, which is the same for peak and RMS values.

> [!warning] Pitfall — special values and the forgotten $n=1$
> Two errors account for most wrong coefficients. (1) Writing $\cos n\pi=1$ instead of $(-1)^n$, or $\sin(n\pi/2)=1$ for all $n$. Keep a table: $\sin n\pi=0$, $\cos n\pi=(-1)^n$, $\sin\frac{n\pi}{2}=0,1,0,-1,\dots$ for $n=0,1,2,3$, $\cos\frac{n\pi}2=1,0,-1,0,\dots$. (2) Dividing by $n-1$ or $n^2-1$ and forgetting to treat $n=1$ separately, as in Worked Example 15.5. Whenever a formula has a factor that vanishes for a particular $n$, compute that coefficient directly.

## 15.3 Complex Fourier series

### Motivation

The trigonometric form needs two families of coefficients and two families of formulas, and its algebra — shifting, differentiating, multiplying by a frequency response — keeps generating products of sines and cosines that must be untangled with identities. The complex exponential form packs each harmonic into a *single* complex number, turns differentiation into multiplication by $jn\omega_0$, turns a time shift into multiplication by a phase factor, and makes the LTI response $c_nH(jn\omega_0)$ a one-line statement. It is also the form that generalises directly to the Fourier transform (Ch. 16), the Z-transform (Ch. 17) and the DFT. In signals and systems work it is the default.

### From real to complex form

Euler's formula gives $\cos\theta=\tfrac12(e^{j\theta}+e^{-j\theta})$ and $\sin\theta=\tfrac1{2j}(e^{j\theta}-e^{-j\theta})$. Substituting into the trigonometric series,

$$
a_n\cos n\omega_0t+b_n\sin n\omega_0t=\frac{a_n-jb_n}{2}e^{jn\omega_0t}+\frac{a_n+jb_n}{2}e^{-jn\omega_0t}.
$$

Defining

$$
c_0=\frac{a_0}{2},\qquad c_n=\frac{a_n-jb_n}{2},\qquad c_{-n}=\frac{a_n+jb_n}{2}\qquad(n\ge1),
$$

the series becomes a single two-sided sum.

> [!definition] Complex (exponential) Fourier series
> The complex Fourier series of a $T_0$-periodic function $x$ is
> $$
> x(t)\sim\sum_{n=-\infty}^{\infty}c_n\,e^{jn\omega_0t},\qquad c_n=\frac1{T_0}\int_{T_0}x(t)\,e^{-jn\omega_0t}\,dt,\quad n\in\mathbb{Z}.
> $$
> The partial sums are taken symmetrically, $S_N(t)=\sum_{n=-N}^{N}c_ne^{jn\omega_0t}$, so that $S_N$ is the same function as the trigonometric partial sum of degree $N$.

The coefficient formula can be derived directly — it is the projection $\langle x,e^{jn\omega_0t}\rangle$ of §15.1 — or checked by combining the Euler formulas: $\frac{a_n-jb_n}{2}=\frac1{T_0}\int x(t)(\cos n\omega_0t-j\sin n\omega_0t)\,dt$. Note the single prefactor $1/T_0$, not $2/T_0$: each real harmonic of amplitude $A_n$ is split equally between two lines of magnitude $A_n/2$ at $\pm n\omega_0$.

> [!theorem] Conversions between coefficient forms
> For $n\ge1$:
> $$
> \begin{aligned}
> a_n&=c_n+c_{-n}, & b_n&=j(c_n-c_{-n}),\\
> A_n&=2\lvert c_n\rvert\ (\text{if }x\text{ real}), & \phi_n&=\angle c_n\ (\text{if }x\text{ real}),
> \end{aligned}
> $$
> and $c_0=a_0/2=A_0$. For a real $x$, $a_n=2\operatorname{Re}\{c_n\}$ and $b_n=-2\operatorname{Im}\{c_n\}$.

The phase relation follows because $c_n=\frac{a_n-jb_n}{2}=\frac{A_n}{2}e^{j\phi_n}$ with $\phi_n=\operatorname{atan2}(-b_n,a_n)$, exactly the phase of the amplitude–phase form in §15.2. Then $c_ne^{jn\omega_0t}+c_{-n}e^{-jn\omega_0t}=2\operatorname{Re}\{c_ne^{jn\omega_0t}\}=A_n\cos(n\omega_0t+\phi_n)$: **$c_n$ is half the phasor of the $n$th harmonic.** This is the bridge to Ch. 2 phasor analysis.

### Real signals: conjugate symmetry

> [!theorem] Conjugate symmetry
> $x(t)$ is real if and only if $c_{-n}=c_n^*$ for all $n$. Then $\lvert c_{-n}\rvert=\lvert c_n\rvert$ (even magnitude spectrum) and $\angle c_{-n}=-\angle c_n$ (odd phase spectrum).

**Proof.** $c_{-n}=\frac1{T_0}\int x(t)e^{jn\omega_0t}dt$ and $c_n^*=\frac1{T_0}\int x^*(t)e^{jn\omega_0t}dt$; these agree for all $n$ iff $x=x^*$ (by uniqueness of coefficients, §15.6). $\square$

The negative-frequency coefficients of a real signal therefore carry no new information; they are the "clockwise-rotating" partners of Ch. 3, §3.5, needed so that the imaginary parts cancel. For a *complex* signal — the analytic signal or the baseband I/Q signal of a radio receiver, $x=I+jQ$ — the positive- and negative-frequency coefficients are independent, and a one-sided description loses information.

### Line spectra, two-sided

Plotting $\lvert c_n\rvert$ and $\angle c_n$ against $\omega=n\omega_0$ (or $f=nf_0$) for $n\in\mathbb{Z}$ gives the **two-sided line spectrum**: a set of vertical lines at the harmonic frequencies. Because the frequencies are discrete, a periodic signal has a *line* (discrete) spectrum, while an aperiodic signal has a continuous spectrum (Ch. 16). Section 15.5 compares one- and two-sided spectra and their conventions.

### Computing $c_n$ directly

The complex formula usually needs fewer integrations than the real pair, because $\int e^{at}dt=e^{at}/a$ handles exponential, sinusoidal and constant pieces at once.

> [!example] Worked Example 15.7 — Pulse train in complex form, with numbers
> A clock signal is a rectangular pulse train of height $A=5$ V, width $\tau=1\ \mu$s and period $T_0=4\ \mu$s, with a pulse centred on $t=0$. Find $c_n$, and list $c_0,\dots,c_4$ and the frequencies at which they occur.
>
> **Solution.** $\omega_0=2\pi/T_0$, $D=\tau/T_0=1/4$. For $n\ne0$,
> $$
> c_n=\frac1{T_0}\int_{-\tau/2}^{\tau/2}Ae^{-jn\omega_0t}\,dt=\frac{A}{T_0}\cdot\frac{e^{-jn\omega_0\tau/2}-e^{jn\omega_0\tau/2}}{-jn\omega_0}=\frac{A}{T_0}\cdot\frac{2\sin(n\omega_0\tau/2)}{n\omega_0}=AD\operatorname{sinc}(nD),
> $$
> using $n\omega_0\tau/2=n\pi D$. For $n=0$ the mean is $AD$, consistent with $\operatorname{sinc}0=1$. The coefficients are real (because $x$ is real and even; §15.4).
>
> | $n$ | $f=n/T_0$ | $c_n=1.25\operatorname{sinc}(n/4)$ (V) |
> |---|---|---|
> | 0 | 0 | 1.2500 |
> | ±1 | 250 kHz | 1.1254 |
> | ±2 | 500 kHz | 0.7958 |
> | ±3 | 750 kHz | 0.3751 |
> | ±4 | 1 MHz | 0 |
>
> The first null is at $f=1/\tau=1$ MHz. Compare with Worked Example 15.6: $a_n=2c_n$, e.g. the 250 kHz harmonic has peak amplitude $2(1.1254)=2.251$ V.
>
> **Check.** $c_1=1.25\cdot\sin(\pi/4)/(\pi/4)=1.25\times0.9003=1.1254$. ✓

> [!example] Worked Example 15.8 — Periodic exponential decay
> The waveform $x(t)=e^{-at}$ for $0\le t<T_0$, repeated with period $T_0$, models a periodically recharged and discharged RC node. Find $c_n$, and evaluate for $a=1\ \text{s}^{-1}$, $T_0=2\pi$ s.
>
> **Solution.**
> $$
> c_n=\frac1{T_0}\int_0^{T_0}e^{-at}e^{-jn\omega_0t}\,dt=\frac1{T_0}\Big[\frac{e^{-(a+jn\omega_0)t}}{-(a+jn\omega_0)}\Big]_0^{T_0}=\frac{1-e^{-aT_0}e^{-jn2\pi}}{T_0(a+jn\omega_0)}=\frac{1-e^{-aT_0}}{T_0\,(a+jn\omega_0)} .
> $$
> With $a=1$, $T_0=2\pi$, $\omega_0=1$: $\frac{1-e^{-2\pi}}{2\pi}=0.158858$, so $c_n=0.158858/(1+jn)$:
> $c_0=0.1589$, $c_1=0.0794-0.0794j=0.1123\angle{-45^\circ}$, $c_2=0.0318-0.0635j=0.0710\angle{-63.4^\circ}$.
>
> **Interpretation.** $\lvert c_n\rvert\approx0.159/n$ for large $n$: a $1/n$ decay, because the waveform jumps from $e^{-aT_0}$ back to 1 at each period boundary. The factor $1/(a+jn\omega_0)$ is the Laplace transform $1/(s+a)$ at $s=jn\omega_0$ — a connection made general next.
>
> **Check.** Numerical quadrature of the defining integral reproduces $c_1=0.079429-0.079429j$ and $c_2=0.031772-0.063543j$. ✓

### Coefficients from the Laplace transform of one period

Let $x_1(t)$ be one period of $x$, i.e. $x_1(t)=x(t)$ on $[0,T_0)$ and $0$ elsewhere, and let $X_1(s)=\mathcal{L}\{x_1\}=\int_0^{T_0}x(t)e^{-st}dt$, an entire function of $s$ because the integral is over a finite interval. Comparing with the definition of $c_n$:

> [!theorem] Coefficients via the one-period transform
> $$
> c_n=\frac{1}{T_0}\,X_1(jn\omega_0),\qquad X_1(s)=\int_0^{T_0}x(t)\,e^{-st}\,dt .
> $$

Chapter 14 showed that the Laplace transform of the periodic extension (starting at $t=0$) is $X(s)=X_1(s)/(1-e^{-sT_0})$. Its poles are at the zeros of $1-e^{-sT_0}$, i.e. at $s=jn\omega_0$, and the residue there is $X_1(jn\omega_0)/T_0=c_n$. So the Fourier coefficients are exactly the residues of the Laplace transform at its imaginary-axis poles — the steady-state part of the inverse Laplace transform. Chapter 16 obtains the same relation as $c_n=X(n\omega_0)/T_0$ from the Fourier transform of one period.

> [!example] Worked Example 15.9 — Half-wave rectifier via the one-period transform
> Recompute the coefficients of Worked Example 15.5 ($V_m=1$, $\omega_0=1$) using $X_1(s)$.
>
> **Solution.** One period is $\sin t$ on $[0,\pi]$ and 0 on $(\pi,2\pi)$. Writing $\sin t\,[u(t)-u(t-\pi)]=\sin t\,u(t)+\sin(t-\pi)u(t-\pi)$ and using the shift theorem of Ch. 14,
> $$
> X_1(s)=\frac{1+e^{-\pi s}}{s^2+1}.
> $$
> For $n\ne\pm1$: $c_n=\frac{1}{2\pi}\cdot\frac{1+(-1)^n}{1-n^2}$, which is $0$ for odd $n$ and $\frac{1}{\pi(1-n^2)}$ for even $n$. Then $a_n=2\operatorname{Re}c_n=-\frac{2}{\pi(n^2-1)}$ — as before. For $n=0$, $c_0=1/\pi$. ✓
>
> For $n=1$, $X_1(s)$ has a removable $0/0$ at $s=j$: write $X_1(s)=\frac{1+e^{-\pi s}}{(s-j)(s+j)}$ and apply l'Hôpital to the numerator over $(s-j)$: $\lim_{s\to j}\frac{-\pi e^{-\pi s}}{1}\cdot\frac{1}{s+j}=\frac{-\pi e^{-j\pi}}{2j}=\frac{\pi}{2j}$. So $c_1=\frac{1}{2\pi}\cdot\frac{\pi}{2j}=-\frac{j}{4}$ and $b_1=-2\operatorname{Im}c_1=\frac12$. ✓
>
> **Lesson.** The one-period transform is often the quickest route when $x_1$ is built from standard Laplace pairs and shifts; the special harmonics where the denominator vanishes need a limit.

> [!example] Worked Example 15.10 — Reading coefficients off without integration
> Find the complex Fourier coefficients of $x(t)=\cos^3\omega_0t+2\sin2\omega_0t$, and of the sawtooth $x(t)=t$ on $(-\pi,\pi)$ from its known sine series.
>
> **Solution.** (a) By the power-reduction formula of Ch. 3, $\cos^3\theta=\frac34\cos\theta+\frac14\cos3\theta$. So
> $$
> x(t)=\tfrac38\big(e^{j\omega_0t}+e^{-j\omega_0t}\big)+\tfrac18\big(e^{j3\omega_0t}+e^{-j3\omega_0t}\big)+\tfrac1j\big(e^{j2\omega_0t}-e^{-j2\omega_0t}\big).
> $$
> Hence $c_{\pm1}=\tfrac38$, $c_{\pm3}=\tfrac18$, $c_2=-j$, $c_{-2}=j$, all others zero. Uniqueness of coefficients (§15.6) guarantees that *any* way of writing $x$ as a finite sum of harmonics gives its Fourier series.
>
> (b) From Worked Example 15.1, $a_n=0$, $b_n=2(-1)^{n+1}/n$. Then $c_n=\frac{a_n-jb_n}{2}=-j\frac{(-1)^{n+1}}{n}=j\frac{(-1)^n}{n}$ for $n\ne0$, $c_0=0$. Check conjugate symmetry: $c_{-n}=j\frac{(-1)^{n}}{-n}=-j\frac{(-1)^n}{n}=c_n^*$. ✓

> [!warning] Pitfall — factor-of-two errors between the forms
> The most frequent mistake in this whole subject is mixing the $1/T_0$ of $c_n$ with the $2/T_0$ of $a_n,b_n$, or reading $\lvert c_n\rvert$ as the amplitude of the $n$th harmonic. For a real signal the amplitude of the harmonic at $n\omega_0$ is $A_n=2\lvert c_n\rvert$; only the DC term has $A_0=c_0$ with no factor 2. A quick sanity test: $\cos\omega_0t$ has $c_{\pm1}=\tfrac12$ and $a_1=1$.

> [!tip] Engineering note — phasors are positive-frequency coefficients
> In AC circuit analysis a voltage $V_m\cos(\omega t+\phi)$ is represented by the phasor $\mathbf{V}=V_me^{j\phi}$ (peak) or $\frac{V_m}{\sqrt2}e^{j\phi}$ (RMS). In Fourier language $\mathbf{V}_{\text{peak}}=2c_1$. A multi-harmonic waveform is then a *list* of phasors $\mathbf{V}_n=2c_n$, one per harmonic, each analysed with the impedances $Z(jn\omega_0)$ at its own frequency — the method of §15.9.

## 15.4 Even and odd functions

### Motivation

Most of the integrals in §15.2 turned out to be zero, and in every case the zero could have been predicted from a symmetry of the waveform. Recognising symmetry before integrating halves (or quarters) the work, provides an independent check on the answer, and — more importantly for engineers — links the *shape* of a waveform to its *harmonic content*: a symmetrical inverter output has no even harmonics, a balanced push–pull amplifier cancels even-order distortion, and an asymmetrical clipping characteristic generates a second harmonic and a DC shift. This section collects the symmetry rules and the half-range expansions that follow from them.

### Even and odd functions

> [!definition] Even and odd functions
> A function is **even** if $x(-t)=x(t)$ for all $t$ (graph symmetric about the vertical axis) and **odd** if $x(-t)=-x(t)$ (graph symmetric under a 180° rotation about the origin). An odd function defined at $0$ has $x(0)=0$.

The rules for products follow the rules for signs: even × even and odd × odd are even; even × odd is odd. For integrals over a symmetric interval,

$$
\int_{-a}^{a}x_{\text{even}}(t)\,dt=2\int_0^ax_{\text{even}}(t)\,dt,\qquad \int_{-a}^{a}x_{\text{odd}}(t)\,dt=0 .
$$

Since $\cos n\omega_0t$ is even and $\sin n\omega_0t$ is odd, applying these rules to the Euler formulas over $(-T_0/2,T_0/2)$ gives:

> [!theorem] Fourier series of even and odd functions
> If $x$ is even, $b_n=0$ for all $n$ and $x$ has a **Fourier cosine series**
> $$
> x(t)\sim\frac{a_0}2+\sum_{n=1}^\infty a_n\cos n\omega_0t,\qquad a_n=\frac{4}{T_0}\int_0^{T_0/2}x(t)\cos n\omega_0t\,dt .
> $$
> If $x$ is odd, $a_n=0$ for all $n$ (including $a_0$) and $x$ has a **Fourier sine series**
> $$
> x(t)\sim\sum_{n=1}^\infty b_n\sin n\omega_0t,\qquad b_n=\frac{4}{T_0}\int_0^{T_0/2}x(t)\sin n\omega_0t\,dt .
> $$
> In complex form: $x$ real and even $\iff$ $c_n$ real and even in $n$; $x$ real and odd $\iff$ $c_n$ purely imaginary and odd in $n$.

**Even/odd decomposition.** Every function splits uniquely into even and odd parts,

$$
x(t)=x_e(t)+x_o(t),\qquad x_e(t)=\frac{x(t)+x(-t)}2,\qquad x_o(t)=\frac{x(t)-x(-t)}2 .
$$

The cosine terms of the Fourier series of $x$ are exactly the Fourier series of $x_e$, and the sine terms are the series of $x_o$. For real $x$: $x_e\leftrightarrow\operatorname{Re}\{c_n\}$ and $x_o\leftrightarrow j\operatorname{Im}\{c_n\}$.

Symmetry is a property of the waveform *and the chosen time origin*. The square wave of Worked Example 15.3 is odd about $t=0$; shifting the origin by $T_0/4$ makes it even, and its series becomes a cosine series with the same amplitudes (§15.7). Many waveforms have **hidden symmetry**: they become even or odd after subtracting their mean value (Worked Example 15.4) or after a shift of origin. Look for both before integrating.

### Half-wave symmetry

> [!definition] Half-wave (odd-harmonic) symmetry
> A $T_0$-periodic function has **half-wave symmetry** if
> $$
> x\Big(t+\frac{T_0}{2}\Big)=-x(t)\qquad\text{for all }t .
> $$
> Graphically, the second half-cycle is the first half-cycle shifted by $T_0/2$ and turned upside down. The mean value is necessarily zero.

> [!theorem] Half-wave symmetry ⇔ odd harmonics only
> $x$ has half-wave symmetry if and only if $c_n=0$ (equivalently $a_n=b_n=0$) for every even $n$, including $n=0$. In that case, for odd $n$,
> $$
> a_n=\frac{4}{T_0}\int_0^{T_0/2}x(t)\cos n\omega_0t\,dt,\qquad b_n=\frac{4}{T_0}\int_0^{T_0/2}x(t)\sin n\omega_0t\,dt .
> $$

**Proof.** Split the coefficient integral into two half-periods and substitute $t=t'+T_0/2$ in the second:

$$
c_n=\frac1{T_0}\int_0^{T_0/2}x(t)e^{-jn\omega_0t}dt+\frac1{T_0}\int_0^{T_0/2}x\Big(t'+\frac{T_0}2\Big)e^{-jn\omega_0t'}e^{-jn\pi}\,dt'.
$$

With half-wave symmetry the second integrand is $-x(t')e^{-jn\omega_0t'}(-1)^n$, so $c_n=\big(1-(-1)^n\big)\frac1{T_0}\int_0^{T_0/2}x(t)e^{-jn\omega_0t}dt$: zero for even $n$ and twice the half-period integral for odd $n$. Conversely, if only odd harmonics are present then each term satisfies $e^{jn\omega_0(t+T_0/2)}=(-1)^ne^{jn\omega_0t}=-e^{jn\omega_0t}$, so the sum does too. $\square$

The physical reading is powerful: **any system that treats positive and negative half-cycles identically (apart from sign) produces only odd harmonics.** A symmetric clipper, a push–pull output stage, a transformer core with a symmetric B–H loop and a full-bridge inverter switched symmetrically all have half-wave-symmetric outputs when driven by a sinusoid. Even harmonics are the signature of *asymmetry*: a half-wave rectifier, a single-ended amplifier stage with a square-law term, or DC magnetisation of a transformer core.

The opposite condition $x(t+T_0/2)=+x(t)$ simply means that $T_0/2$ is a period; then only even harmonics of $\omega_0$ appear — they are all the harmonics of the true fundamental $2\omega_0$. The full-wave rectified sine $\lvert\sin\omega_0t\rvert$ is the standard example.

```desmos-graph
left=-1; right=13; bottom=-1.8; top=1.8
---
y=\sin(x)+0.4\sin(3x+1)|blue
y=\sin(x)+0.4\sin(2x+1)|red
y=0|black|dashed
```
*Blue: fundamental plus a third harmonic — each negative half-cycle is the inverted copy of the preceding positive one (half-wave symmetry). Red: fundamental plus a second harmonic — the positive and negative half-cycles have different shapes, and no shift-and-invert maps one onto the other.*

### Quarter-wave symmetry

> [!definition] Quarter-wave symmetry
> A function has **quarter-wave symmetry** if it has half-wave symmetry *and* is either even or odd (about a suitable origin).

Then only odd harmonics are present and, in addition, only cosines (even case) or only sines (odd case). Combining the two reductions, one quarter-period suffices:

> [!theorem] Quarter-wave formulas
> If $x$ is odd with half-wave symmetry, then $a_n=0$ for all $n$, $b_n=0$ for even $n$, and for odd $n$
> $$
> b_n=\frac{8}{T_0}\int_0^{T_0/4}x(t)\sin n\omega_0t\,dt .
> $$
> If $x$ is even with half-wave symmetry, then $b_n=0$ for all $n$, $a_n=0$ for even $n$, and for odd $n$
> $$
> a_n=\frac{8}{T_0}\int_0^{T_0/4}x(t)\cos n\omega_0t\,dt .
> $$

**Proof sketch (odd case).** Odd symmetry gives $b_n=\frac4{T_0}\int_0^{T_0/2}x\sin n\omega_0t\,dt$. Odd plus half-wave implies $x(T_0/2-t)=-x(t-T_0/2)=x(t)$: on $(0,T_0/2)$ the graph is symmetric about $t=T_0/4$. For odd $n$, $\sin n\omega_0(T_0/2-t)=\sin(n\pi-n\omega_0t)=\sin n\omega_0t$, so the integrand is symmetric about $T_0/4$ too and the half-period integral is twice the quarter-period integral. $\square$

| Symmetry | Condition | Harmonics present | Coefficient integral |
|---|---|---|---|
| none | — | all, $a_n$ and $b_n$ | full period, $2/T_0$ |
| even | $x(-t)=x(t)$ | DC and cosines | half period, $4/T_0$ |
| odd | $x(-t)=-x(t)$ | sines only, no DC | half period, $4/T_0$ |
| half-wave | $x(t+T_0/2)=-x(t)$ | odd $n$ only, no DC | half period, $4/T_0$ |
| even quarter-wave | even and half-wave | cosines, odd $n$ | quarter period, $8/T_0$ |
| odd quarter-wave | odd and half-wave | sines, odd $n$ | quarter period, $8/T_0$ |

> [!example] Worked Example 15.11 — Triangle wave by the quarter-wave formula
> The triangle wave of peak $A$ and period $T_0$ has $x(t)=A\big(1-4\lvert t\rvert/T_0\big)$ on $(-T_0/2,T_0/2)$. Find its Fourier series using symmetry.
>
> **Solution.** Symmetry: $x$ is even; and $x(t+T_0/2)=-x(t)$ (check at $t=0$: $x(T_0/2)=A(1-2)=-A=-x(0)$; the graph is a straight line falling from $A$ to $-A$ on $(0,T_0/2)$ and rising back on $(T_0/2,T_0)$). So it is even quarter-wave: only $a_n$ with odd $n$.
>
> Let $\theta=\omega_0t$; on $(0,\pi/2)$, $x=A(1-2\theta/\pi)$. For odd $n$,
> $$
> a_n=\frac{8}{T_0}\int_0^{T_0/4}x\cos n\omega_0t\,dt=\frac4\pi\int_0^{\pi/2}A\Big(1-\frac{2\theta}{\pi}\Big)\cos n\theta\,d\theta .
> $$
> By parts, $\int_0^{\pi/2}(1-\tfrac{2\theta}{\pi})\cos n\theta\,d\theta=\big[(1-\tfrac{2\theta}\pi)\tfrac{\sin n\theta}{n}\big]_0^{\pi/2}+\frac{2}{\pi n}\int_0^{\pi/2}\sin n\theta\,d\theta=0+\frac{2}{\pi n^2}\Big(1-\cos\frac{n\pi}2\Big)=\frac{2}{\pi n^2}$ for odd $n$. Hence
> $$
> a_n=\frac{8A}{\pi^2n^2}\ (n\text{ odd}),\qquad x(t)=\frac{8A}{\pi^2}\Big(\cos\omega_0t+\frac19\cos3\omega_0t+\frac1{25}\cos5\omega_0t+\cdots\Big).
> $$
>
> **Checks.** At $t=0$: $\frac{8A}{\pi^2}\sum_{\text{odd}}\frac1{n^2}=\frac{8A}{\pi^2}\cdot\frac{\pi^2}{8}=A$. ✓ The coefficients fall as $1/n^2$ — the triangle is continuous with jumps only in its slope (§15.7). The fundamental, $0.811A$, is already very close to the waveform: the third harmonic is only $1/9$ of it.

> [!example] Worked Example 15.12 — Quasi-square (120° conduction) waveform
> A three-phase six-step inverter produces a line-to-line voltage that is $+V$ for $30^\circ<\theta<150^\circ$, $-V$ for $210^\circ<\theta<330^\circ$ and zero otherwise ($\theta=\omega_0t$). Find its harmonic content.
>
> **Solution.** With the origin as given, $x$ is odd ($x(-\theta)=-x(\theta)$, since $-V$ occupies $(-150^\circ,-30^\circ)\equiv(210^\circ,330^\circ)$) and half-wave symmetric. So only sine terms with odd $n$; use the quarter-wave formula with $T_0=2\pi$:
> $$
> b_n=\frac{8}{2\pi}\int_0^{\pi/2}x\sin n\theta\,d\theta=\frac{4V}{\pi}\int_{\pi/6}^{\pi/2}\sin n\theta\,d\theta=\frac{4V}{n\pi}\Big[\cos\frac{n\pi}6-\cos\frac{n\pi}{2}\Big]=\frac{4V}{n\pi}\cos\frac{n\pi}{6}\quad(n\text{ odd}).
> $$
> - $n=3,9,15,\dots$ (triplen): $\cos(n\pi/6)=\cos(k\pi/2)$ with odd $k$, which is $0$. **No triplen harmonics** — the line-to-line voltage is the difference of two phase waveforms displaced by $120^\circ$, and triplens are in phase in all three phases (Ch. 3, §3.7), so they cancel.
> - $n=6k\pm1$: $\lvert\cos(n\pi/6)\rvert=\sqrt3/2$, so $\lvert b_n\rvert=\frac{2\sqrt3V}{n\pi}=\frac{1.1027V}{n}$.
>
> Numbers: $b_1=1.1027V$, $b_5=-0.2205V$, $b_7=-0.1575V$, $b_{11}=0.1002V$, $b_{13}=0.0848V$ — each exactly $1/n$ of the fundamental. **Check:** direct quadrature over a full period reproduces all five values to $10^{-12}$. ✓ The distortion of this waveform is computed with Parseval's theorem in §15.8.

### Half-range expansions

In boundary-value problems ([Ch. 13](../part-4-differential-equations/13-partial-differential-equations.md)) a function is given only on a half-interval $0<x<L$ — the initial temperature along a rod, the initial displacement of a string, the voltage on one face of a box — and we are free to extend it to $(-L,0)$ however suits the boundary conditions. Two extensions are standard.

> [!definition] Half-range expansions
> Let $f$ be given on $(0,L)$.
>
> The **even extension** ($f(-x)=f(x)$, period $2L$) yields the **half-range cosine series**
> $$
> f(x)\sim\frac{a_0}{2}+\sum_{n=1}^{\infty}a_n\cos\frac{n\pi x}{L},\qquad a_n=\frac2L\int_0^Lf(x)\cos\frac{n\pi x}{L}\,dx .
> $$
> The **odd extension** ($f(-x)=-f(x)$, period $2L$) yields the **half-range sine series**
> $$
> f(x)\sim\sum_{n=1}^{\infty}b_n\sin\frac{n\pi x}{L},\qquad b_n=\frac2L\int_0^Lf(x)\sin\frac{n\pi x}{L}\,dx .
> $$

Both series represent $f$ on $(0,L)$; they differ outside it and at the endpoints. The sine series vanishes at $x=0$ and $x=L$, suiting Dirichlet conditions $u=0$ at the ends; the cosine series has zero slope there (each term does), suiting Neumann (insulated) conditions. Note the harmonics are $\pi/L$ apart, half the spacing $2\pi/L$ of an ordinary series of period $L$ — the extension doubles the period.

> [!example] Worked Example 15.13 — Sine and cosine expansions of a ramp
> Expand $f(x)=x$ on $(0,L)$ in (a) a half-range sine series and (b) a half-range cosine series, and compare them.
>
> **Solution.** (a) Odd extension is the sawtooth $x$ on $(-L,L)$. By parts,
> $$
> b_n=\frac2L\int_0^Lx\sin\frac{n\pi x}{L}\,dx=\frac2L\Big[-\frac{Lx}{n\pi}\cos\frac{n\pi x}{L}\Big]_0^L+\frac2L\cdot\frac{L}{n\pi}\int_0^L\cos\frac{n\pi x}{L}dx=\frac{2L(-1)^{n+1}}{n\pi}.
> $$
> So $x=\frac{2L}{\pi}\big(\sin\frac{\pi x}L-\frac12\sin\frac{2\pi x}L+\frac13\sin\frac{3\pi x}{L}-\cdots\big)$ on $(0,L)$.
>
> (b) Even extension is the triangle $\lvert x\rvert$ on $(-L,L)$. $a_0=\frac2L\int_0^Lx\,dx=L$, and for $n\ge1$
> $$
> a_n=\frac2L\int_0^Lx\cos\frac{n\pi x}{L}\,dx=\frac2L\Big[\frac{L^2}{n^2\pi^2}\cos\frac{n\pi x}{L}\Big]_0^L=\frac{2L}{n^2\pi^2}\big((-1)^n-1\big),
> $$
> i.e. $-\frac{4L}{n^2\pi^2}$ for odd $n$, $0$ for even $n$. So $x=\frac L2-\frac{4L}{\pi^2}\big(\cos\frac{\pi x}{L}+\frac19\cos\frac{3\pi x}L+\cdots\big)$ on $(0,L)$.
>
> **Comparison.** The odd extension jumps from $L$ to $-L$ at $x=L$, so the sine series converges slowly ($1/n$), shows a Gibbs overshoot near $x=L$ (§15.6) and converges to $0$ (the midpoint) at $x=L$. The even extension is continuous, so the cosine series converges uniformly with $1/n^2$ coefficients and gives the right value $L$ at $x=L$. With $L=1$ and three non-zero harmonics each (sine: $n=1,2,3$; cosine: $n=1,3,5$): at $x=0.9$ the sine series gives $0.556$ and the cosine series $0.912$ (exact $0.9$); at $x=1$ the sine series gives $0$, the cosine series $0.967$ (exact $1$). Choose the extension that is continuous whenever the boundary conditions leave you the choice.

> [!example] Worked Example 15.14 — Even and odd parts from the coefficients
> For the periodic exponential of Worked Example 15.8 ($e^{-t}$ on $[0,2\pi)$), find the Fourier series of its even part without further integration.
>
> **Solution.** $c_n=\frac{0.158858}{1+jn}=0.158858\frac{1-jn}{1+n^2}$. The even part has coefficients $\operatorname{Re}\{c_n\}=\frac{0.158858}{1+n^2}$, so
> $$
> x_e(t)=0.158858\Big(1+2\sum_{n=1}^\infty\frac{\cos nt}{1+n^2}\Big),
> $$
> and the odd part is $x_o(t)=0.317716\sum_{n\ge1}\frac{n\sin nt}{1+n^2}$. Directly, $x_e(t)=\frac12\big(e^{-t}+e^{-(2\pi-t)}\big)$ on $(0,2\pi)$, which is $\frac12(1+e^{-2\pi})=0.500934$ at $t=0^+$. **Check:** summing the cosine series at $t=0$ to 200 000 terms gives $0.158858\,(1+2\times1.076674)=0.50093$. ✓ (The continuous even part's series has $1/n^2$ coefficients; the odd part carries the jump and has $1/n$ coefficients.)

> [!warning] Pitfall — symmetry depends on the origin, and on the whole period
> A waveform that "looks odd" over the half-period you sketched may not be odd as a periodic function; check $x(-t)=-x(t)$ using the periodic extension (for instance $x(-t)=x(T_0-t)$). And half-wave symmetry is a statement about shifting by exactly $T_0/2$: a waveform whose negative half-cycle is a *mirror image in time* (rather than a shifted copy) of the positive one need not have half-wave symmetry.

> [!tip] Engineering note — even-order distortion and push–pull cancellation
> If an amplifier's transfer characteristic is $y=a_1x+a_2x^2+a_3x^3$, a push–pull pair of matched devices produces $y_1-y_2$ with $y_2$ the same characteristic driven by $-x$: the even-order terms cancel, leaving $2a_1x+2a_3x^3$, an odd function of $x$. An odd function of a half-wave-symmetric input is half-wave symmetric, so the output has no even harmonics. Measured second-harmonic distortion in such an amplifier is therefore a direct measure of device mismatch.

## 15.5 Harmonics

### Motivation

Once a waveform has been decomposed, engineers stop thinking about the waveform and start thinking about its **harmonics**: "the 5th harmonic current exceeds the limit", "the third harmonic of the oscillator is 40 dB down", "the triplen harmonics overload the neutral". This section fixes the vocabulary, defines line spectra precisely in their one-sided and two-sided forms, and records how simple operations on a waveform — shifting it, compressing it in time, scaling it — move and reshape its spectrum.

### Vocabulary

> [!definition] Harmonics
> For a periodic signal with fundamental angular frequency $\omega_0$:
> - the **fundamental** (first harmonic) is the component at $\omega_0$;
> - the **$n$th harmonic** is the component $A_n\cos(n\omega_0t+\phi_n)$ at $n\omega_0$, $n\ge2$;
> - the **DC component** is $A_0=c_0$, sometimes called the zeroth harmonic;
> - **odd** and **even** harmonics have odd and even $n$; harmonics with $n$ a multiple of 3 are called **triplen** harmonics.

Three related terms are often confused:

- An **overtone** (acoustics) is any component above the fundamental; the first overtone is the *second* harmonic.
- An **interharmonic** is a component at a frequency that is *not* an integer multiple of the supply frequency, for instance a 175 Hz component in a 50 Hz system produced by a cycloconverter or an arc furnace. A signal containing interharmonics is periodic only with a longer period (if the frequency ratio is rational) or not at all (Ch. 3, §3.4), and its Fourier series with respect to 50 Hz does not exist.
- A **subharmonic** has a frequency $\omega_0/m$ below the fundamental. Its presence means the true period is $mT_0$; relative to that longer period, it is simply a low-order harmonic.

### One-sided and two-sided line spectra

For a real signal there are two equivalent ways to draw the spectrum.

> [!definition] Line spectra
> The **one-sided amplitude and phase spectra** are the sets of lines $A_n$ and $\phi_n$ at frequencies $n\omega_0\ge0$, from the form $x(t)=A_0+\sum_{n\ge1}A_n\cos(n\omega_0t+\phi_n)$.
>
> The **two-sided magnitude and phase spectra** are the lines $\lvert c_n\rvert$ and $\angle c_n$ at $n\omega_0$ for all $n\in\mathbb{Z}$, from $x(t)=\sum_nc_ne^{jn\omega_0t}$.
>
> The **power spectrum** consists of lines $\lvert c_n\rvert^2$ (two-sided) or $A_0^2$, $A_n^2/2$ (one-sided), whose sum is the mean-square value (Parseval, §15.8).

The conversion is: one-sided line height $=2\times$ two-sided height for $n\ge1$; equal for DC. For phase, the one-sided spectrum is the right half of the two-sided one, and the left half is its negative. Lines with $A_n=0$ have undefined phase and are omitted from phase plots.

Phase conventions matter when lines have negative real coefficients. A coefficient $a_n=-0.2$ can be plotted as amplitude $-0.2$ with phase $0$, or (properly, since amplitudes are non-negative) as $A_n=0.2$ with phase $\pm180^\circ$. Both appear in textbooks; the non-negative-amplitude convention is the one used by instruments.

Amplitudes are frequently given in decibels relative to the fundamental, $20\log_{10}(A_n/A_1)$ dBc, or relative to a reference level such as 1 µV (dBµV, used in EMC). Because harmonics often span several orders of magnitude, a log scale is the only practical display.

```desmos-graph
left=-9; right=9; bottom=-0.4; top=1.4
---
y=1.25\frac{\sin(\pi x/4)}{\pi x/4}|grey|dashed
(0,1.25)|blue
(1,1.1254)|blue
(-1,1.1254)|blue
(2,0.7958)|blue
(-2,0.7958)|blue
(3,0.3751)|blue
(-3,0.3751)|blue
(4,0)|blue
(-4,0)|blue
(5,-0.2251)|blue
(-5,-0.2251)|blue
(6,-0.2653)|blue
(-6,-0.2653)|blue
(7,-0.1608)|blue
(-7,-0.1608)|blue
(8,0)|blue
(-8,0)|blue
```
*Two-sided line spectrum $c_n=1.25\operatorname{sinc}(n/4)$ of the 5 V, 25%-duty pulse train of Worked Example 15.7 (blue points, horizontal axis = harmonic number $n$), lying on the continuous sinc envelope (dashed). Every 4th harmonic is missing. Plotted as magnitudes, lines 5–7 would be drawn upwards with phase $180^\circ$.*

### How operations on the waveform change the spectrum

Let $x(t)\leftrightarrow c_n$ denote that $x$ has coefficients $c_n$ with fundamental $\omega_0$. Each of the following follows in one line from the coefficient integral.

**Amplitude scaling and adding a constant.** $Kx(t)+B\leftrightarrow Kc_n+B\,\delta_{n0}$: scaling scales every line; adding a DC level changes only the $n=0$ line.

**Time shift.** $x(t-t_0)\leftrightarrow c_ne^{-jn\omega_0t_0}$. Substitute $t'=t-t_0$ in $\frac1{T_0}\int x(t-t_0)e^{-jn\omega_0t}dt$. The magnitudes are unchanged; each phase shifts by $-n\omega_0t_0$, **linearly in frequency**. A delay is invisible on an amplitude spectrum and is a straight line of slope $-t_0$ on a phase-versus-frequency plot.

**Time scaling.** $x(at)$, $a>0$, has period $T_0/a$ and fundamental $a\omega_0$, with the *same* coefficients $c_n$: compressing a waveform in time moves every line outward in frequency in proportion, without changing heights. A 1 kHz square wave and a 1 MHz square wave have identical coefficient lists.

**Time reversal.** $x(-t)\leftrightarrow c_{-n}$; for real $x$ this is $c_n^*$, which negates every phase.

**Harmonic content versus shape.** The amplitude spectrum alone does not determine the waveform: phases matter. Ch. 3 showed a fundamental plus a third harmonic in phase and in anti-phase, with identical $A_n$, giving a flattened and a peaked waveform. Peak value, crest factor and the risk of clipping depend on phases; power and RMS do not (§15.8).

> [!example] Worked Example 15.15 — One- and two-sided spectra of a rectified sine
> Draw up the one-sided amplitude and phase spectrum and the two-sided spectrum of the half-wave rectified 230 V (RMS), 50 Hz mains, $V_m=325.3$ V, up to the 6th harmonic, with amplitudes also in dBc.
>
> **Solution.** From Worked Example 15.5, $x=V_m\big[\frac1\pi+\frac12\sin\omega_0t-\frac2\pi\sum_k\frac{\cos2k\omega_0t}{4k^2-1}\big]$. Convert to cosine-reference amplitude–phase form: $\frac12\sin\omega_0t=\frac12\cos(\omega_0t-90^\circ)$, and $-\frac{2}{\pi(4k^2-1)}\cos2k\omega_0t=\frac{2}{\pi(4k^2-1)}\cos(2k\omega_0t+180^\circ)$.
>
> | $n$ | $f$ (Hz) | $A_n/V_m$ | $A_n$ (V) | $\phi_n$ | dBc | $c_n$ (V) |
> |---|---|---|---|---|---|---|
> | 0 | 0 | 0.3183 | 103.5 | — | — | 103.5 |
> | 1 | 50 | 0.5000 | 162.6 | $-90^\circ$ | 0 | $81.3\angle{-90^\circ}$ |
> | 2 | 100 | 0.2122 | 69.0 | $180^\circ$ | $-7.4$ | $34.5\angle180^\circ$ |
> | 3, 5 | 150, 250 | 0 | 0 | — | — | 0 |
> | 4 | 200 | 0.0424 | 13.8 | $180^\circ$ | $-21.4$ | $6.9\angle180^\circ$ |
> | 6 | 300 | 0.0182 | 5.9 | $180^\circ$ | $-28.8$ | $3.0\angle180^\circ$ |
>
> The two-sided spectrum has $c_{-n}=c_n^*$: lines of the same heights at $-50,-100,\dots$ Hz, with phases $+90^\circ$ at $-50$ Hz and $180^\circ$ (equivalently $-180^\circ$) at the even negative harmonics.
>
> **Check.** $A_2/A_1=0.2122/0.5=0.4244$, $20\log_{10}0.4244=-7.44$ dB. ✓ $A_6/V_m=2/(35\pi)=0.01819$. ✓

> [!example] Worked Example 15.16 — Spectrum of a delayed pulse train
> The pulse train of Worked Example 15.7 is delayed so that each pulse starts, rather than is centred, at $t=0$: $x(t)=5$ V for $0<t<1\ \mu$s, period $4\ \mu$s. Find the new coefficients.
>
> **Solution.** The new waveform is the old one delayed by $t_0=\tau/2=0.5\ \mu$s. By the shift rule,
> $$
> c_n^{\text{new}}=1.25\operatorname{sinc}(n/4)\,e^{-jn\omega_0t_0}=1.25\operatorname{sinc}(n/4)\,e^{-jn\pi/4},
> $$
> since $\omega_0t_0=2\pi\cdot0.5/4=\pi/4$. Magnitudes are unchanged; the phase of harmonic $n$ is $-45^\circ n$ (plus $180^\circ$ where sinc is negative). For example $c_1=1.1254\angle{-45^\circ}$, $c_2=0.7958\angle{-90^\circ}$, $c_3=0.3751\angle{-135^\circ}$.
>
> **Check by direct integration.** $c_1=\frac{5}{4\,\mu\text{s}}\int_0^{1\,\mu\text{s}}e^{-j\omega_0t}dt=\frac{5}{4}\cdot\frac{1-e^{-j\pi/2}}{j\pi/2}=\frac54\cdot\frac{1+j}{j\pi/2}=\frac{5}{2\pi}(1-j)=0.7958-0.7958j=1.1254\angle{-45^\circ}$. ✓

### Harmonics in three-phase systems

In a balanced three-phase system the phase-$b$ and phase-$c$ waveforms are the phase-$a$ waveform delayed by $T_0/3$ and $2T_0/3$. By the shift rule the $n$th harmonic of phase $b$ is that of phase $a$ shifted by $-n\cdot120^\circ$. Hence:

- $n=1,4,7,\dots$ ($n=3k+1$): shifts $-120^\circ,-240^\circ$ — **positive sequence** $a$-$b$-$c$, like the fundamental;
- $n=2,5,8,\dots$ ($n=3k+2$): shifts $-240^\circ\equiv+120^\circ$ — **negative sequence** $a$-$c$-$b$, producing reverse torque in induction motors;
- $n=3,6,9,\dots$ (triplen): shifts that are multiples of $360^\circ$ — **zero sequence**, in phase in all three lines.

Zero-sequence currents add in the neutral instead of cancelling, and cancel in any line-to-line quantity (Worked Example 15.12).

> [!example] Worked Example 15.17 — Neutral current from triplen harmonics
> Each phase of a balanced four-wire system supplies a nonlinear load drawing $I_1=10$ A (RMS) at 50 Hz and a third harmonic $I_3=3$ A (RMS), plus negligible higher harmonics. Find the RMS neutral current, and compare with the phase current.
>
> **Solution.** The neutral current is $i_a+i_b+i_c$. Fundamental components are $120^\circ$ apart and sum to zero; third-harmonic components are in phase and add: $i_N=3\times\sqrt2\,(3)\cos(3\omega_0t+\phi_3)$. So $I_N=9$ A RMS at 150 Hz.
>
> The phase RMS current is $\sqrt{10^2+3^2}=10.44$ A (orthogonality, §15.8). The neutral carries 86% of the phase current even though the load is perfectly balanced — the reason modern installations with many switch-mode loads size the neutral conductor at or above the phase conductors.

> [!warning] Pitfall — "the 3rd harmonic" of what?
> Harmonic numbers refer to a stated fundamental. A 12-pulse rectifier's characteristic harmonics are $12k\pm1$ of the *supply* frequency; the ripple of a full-wave rectifier is at the 2nd harmonic of the supply but the *fundamental* of the ripple waveform. When writing or reading a specification, say which fundamental the harmonic numbers refer to.

> [!tip] Engineering note — harmonic limits
> Power-quality standards (e.g. IEC 61000-3-2 for equipment currents, IEEE 519 for system voltages and currents) specify limits harmonic by harmonic, typically up to the 40th or 50th, and in total via THD (§15.8). Measurement standards (IEC 61000-4-7) prescribe exactly the Fourier analysis of this chapter applied to a 10-cycle (50 Hz) or 12-cycle (60 Hz) window, which gives 5 Hz resolution and therefore also resolves interharmonics.

## 15.6 Convergence

### Motivation

So far we have computed coefficients and written "$\sim$". An engineer is entitled to ask whether the sum really *is* the function. The question is not pedantic. The square-wave series converges to $0$ at the jumps, where the square wave is $\pm A$ — so at some points the answer is "no". Its partial sums overshoot by 9% near every jump no matter how many terms are taken — so even "approximately equal" needs qualification. And the identities we want to use — Parseval's theorem, term-by-term differentiation, the LTI response formula — each require a particular kind of convergence. This section states the three standard modes of convergence (pointwise, uniform, mean-square), the conditions under which each holds, and the Gibbs phenomenon that marks the boundary between them.

### Dirichlet conditions and pointwise convergence

> [!definition] Dirichlet conditions
> A $T_0$-periodic function $x$ satisfies the **Dirichlet conditions** if, over one period,
> 1. $x$ is absolutely integrable: $\int_{T_0}\lvert x(t)\rvert\,dt<\infty$;
> 2. $x$ has a finite number of maxima and minima;
> 3. $x$ has a finite number of discontinuities, each of them finite (a jump).
>
> Conditions 2 and 3 together say that $x$ is of bounded variation on a period. A simpler sufficient set used in many texts is that $x$ be **piecewise smooth**: the period splits into finitely many intervals on each of which $x$ and $x'$ are continuous and have finite one-sided limits at the ends.

> [!theorem] Dirichlet's theorem (pointwise convergence)
> If $x$ satisfies the Dirichlet conditions, its Fourier series converges at every $t$, and
> $$
> \lim_{N\to\infty}S_N(t)=\frac{x(t^+)+x(t^-)}{2},
> $$
> where $x(t^\pm)$ are the one-sided limits. At points of continuity the series converges to $x(t)$; at a jump it converges to the **midpoint** of the jump.

All waveforms met in practice — square, triangle, sawtooth, rectified sines, pulse trains, PWM signals, exponential segments — satisfy the Dirichlet conditions. Functions that violate them are mathematical constructions such as $\sin(1/t)$ near $t=0$ (infinitely many extrema) or $\tan t$ (not integrable).

**Proof sketch.** The argument has three steps and is worth knowing because the same kernel reappears in the DFT and in antenna-array theory.

*Step 1 — the Dirichlet kernel.* Use $T_0=2\pi$. Substituting the coefficient integral into the partial sum and exchanging the finite sum with the integral,

$$
S_N(t)=\sum_{n=-N}^{N}\Big[\frac1{2\pi}\int_{-\pi}^{\pi}x(s)e^{-jns}ds\Big]e^{jnt}=\frac1{2\pi}\int_{-\pi}^{\pi}x(s)\,D_N(t-s)\,ds,
$$

where, summing the geometric series as in Ch. 3, §3.6,

$$
D_N(u)=\sum_{n=-N}^{N}e^{jnu}=\frac{\sin\big((N+\tfrac12)u\big)}{\sin(u/2)} .
$$

$D_N$ is even, $2\pi$-periodic, has $\frac{1}{2\pi}\int_{-\pi}^{\pi}D_N(u)\,du=1$ (only the $n=0$ term survives), and has a tall central peak of height $2N+1$ and width about $2\pi/(N+\frac12)$ flanked by oscillating side lobes. So $S_N(t)$ is a weighted local average of $x$ near $t$, with weights that become more sharply concentrated as $N$ grows — but that also oscillate in sign.

*Step 2 — split at the point.* Changing variable to $u=t-s$ and using evenness of $D_N$,

$$
S_N(t)-\frac{x(t^+)+x(t^-)}{2}=\frac{1}{2\pi}\int_{0}^{\pi}\Big[x(t+u)-x(t^+)+x(t-u)-x(t^-)\Big]D_N(u)\,du,
$$

because $\frac{1}{2\pi}\int_0^\pi D_N=\frac12$.

*Step 3 — Riemann–Lebesgue.* Write $D_N(u)=\sin\big((N+\tfrac12)u\big)/\sin(u/2)$. The function $g(u)=\big[x(t+u)-x(t^+)+x(t-u)-x(t^-)\big]/\sin(u/2)$ is integrable on $(0,\pi)$: near $u=0$ the bracket is $O(u)$ when $x$ is piecewise smooth (one-sided derivatives exist), cancelling the $u/2$ in the denominator. The **Riemann–Lebesgue lemma** states that $\int_a^bg(u)\sin(\lambda u)\,du\to0$ as $\lambda\to\infty$ for every integrable $g$ — rapid oscillation averages any fixed integrable function to zero (proof: true for step functions by direct integration; approximate $g$ by step functions). With $\lambda=N+\frac12$ the right-hand side tends to zero. $\square$

The Riemann–Lebesgue lemma applied to the coefficient integral itself gives the companion fact $c_n\to0$ as $\lvert n\rvert\to\infty$ for every absolutely integrable $x$.

> [!theorem] Uniqueness
> If two piecewise-continuous $T_0$-periodic functions have the same Fourier coefficients, they are equal at every point where both are continuous. In particular, a function whose coefficients are all zero vanishes except at isolated points.

Uniqueness is what justifies reading off coefficients from any convenient expansion (Worked Example 15.10) and solving for unknown waveforms by matching harmonics.

### Uniform convergence

Pointwise convergence can be very uneven: near a jump, the number of terms needed for a given accuracy grows without bound as $t$ approaches the jump. **Uniform convergence** means $\max_t\lvert S_N(t)-x(t)\rvert\to0$: a single $N$ achieves a given accuracy everywhere.

> [!theorem] Uniform convergence
> (a) If $\sum_n\lvert c_n\rvert<\infty$, the Fourier series converges absolutely and uniformly to a continuous function (Weierstrass M-test), which equals $x$ wherever $x$ is continuous.
>
> (b) If $x$ is continuous everywhere (including across the period boundary, $x(t_1+T_0)=x(t_1)$) and piecewise smooth, then $\sum\lvert c_n\rvert<\infty$ and the series converges uniformly to $x$.
>
> (c) If $x$ has a jump, the convergence cannot be uniform near that jump: a uniform limit of continuous functions is continuous.

**Proof of (b).** Let $x'$ have coefficients $d_n$. Integrating by parts, with the boundary terms cancelling by continuity and periodicity, $d_n=jn\omega_0c_n$ (§15.7). By Cauchy–Schwarz and Bessel's inequality applied to $x'$,

$$
\sum_{n\ne0}\lvert c_n\rvert=\sum_{n\ne0}\frac{\lvert d_n\rvert}{\lvert n\rvert\omega_0}\le\frac{1}{\omega_0}\Big(\sum_{n\ne0}\frac{1}{n^2}\Big)^{1/2}\Big(\sum_{n}\lvert d_n\rvert^2\Big)^{1/2}\le\frac1{\omega_0}\sqrt{\frac{\pi^2}{3}}\,\lVert x'\rVert<\infty. \quad\square
$$

So the triangle wave, the full-wave rectified sine and the trapezoid converge uniformly; the square wave, sawtooth, pulse trains and the periodic exponential of Worked Example 15.8 do not.

### Mean-square convergence

The weakest but most useful mode is convergence in the RMS sense, $\lVert x-S_N\rVert\to0$. It ignores what happens at individual points and measures only the power of the error — exactly the quantity that matters for energy, noise and distortion.

> [!theorem] Mean-square convergence (completeness)
> If $\frac1{T_0}\int_{T_0}\lvert x\rvert^2dt<\infty$, then
> $$
> \lim_{N\to\infty}\frac1{T_0}\int_{T_0}\big\lvert x(t)-S_N(t)\big\rvert^2dt=0 .
> $$
> Conversely (Riesz–Fischer), every sequence with $\sum\lvert c_n\rvert^2<\infty$ is the coefficient sequence of some square-integrable $x$.

**Proof sketch.** By the best-approximation theorem of §15.1, $\lVert x-S_N\rVert\le\lVert x-p\rVert$ for every trigonometric polynomial $p$ of degree $N$. So it suffices that trigonometric polynomials can approximate $x$ arbitrarily well in RMS. First approximate $x$ in RMS by a continuous periodic function $g$ (possible for any square-integrable $x$). Then approximate $g$ *uniformly* — hence also in RMS — by the **Fejér (Cesàro) means**
$$
\sigma_N(t)=\frac{S_0(t)+S_1(t)+\cdots+S_{N}(t)}{N+1}=\frac1{2\pi}\int_{-\pi}^{\pi}g(t-u)F_N(u)\,du,\qquad F_N(u)=\frac{1}{N+1}\Big[\frac{\sin\big((N+1)u/2\big)}{\sin(u/2)}\Big]^2 .
$$
The **Fejér kernel** $F_N$ is *non-negative*, has unit mean, and concentrates at $u=0$; a positive averaging kernel of this kind converges uniformly for continuous $g$ (Fejér's theorem). The triangle inequality completes the argument. $\square$

Mean-square convergence is what makes Parseval's theorem (§15.8) true for *every* finite-power periodic signal, square wave included, and it is the natural setting for the Fourier-series solutions of PDEs in Ch. 13.

| Mode | Measure of error | Holds when | Square wave? |
|---|---|---|---|
| pointwise | $\lvert x(t)-S_N(t)\rvert$ at each $t$ | Dirichlet conditions (limit = midpoint at jumps) | yes, with midpoint at jumps |
| uniform | $\max_t\lvert x-S_N\rvert$ | $x$ continuous and piecewise smooth, or $\sum\lvert c_n\rvert<\infty$ | no |
| mean-square | $\lVert x-S_N\rVert$ | $x$ square-integrable | yes |

> [!warning] Pitfall — continuity alone is not enough
> There exist continuous periodic functions whose Fourier series diverge at a point (du Bois-Reymond, 1873). Continuity guarantees mean-square convergence and Fejér-mean convergence, but pointwise convergence of the ordinary partial sums needs something more — piecewise smoothness, bounded variation (Dirichlet), or a Hölder condition. For engineering waveforms this is never an issue, but it explains why theorems are stated with conditions that look fussier than "continuous".

> [!example] Worked Example 15.18 — Using pointwise convergence to sum series
> Use the convergence theorem to evaluate (a) $1-\frac13+\frac15-\cdots$, (b) $\sum_{k\ge1}\frac{1}{(2k-1)^2}$, (c) the value of the sawtooth series $\sum_{n\ge1}\frac{2(-1)^{n+1}}{n}\sin nt$ at $t=\pi$.
>
> **Solution.** (a) The square wave ($A=1$, $\omega_0=1$) is continuous at $t=\pi/2$ with value $1$: $\frac4\pi\big(\sin\frac\pi2+\frac13\sin\frac{3\pi}2+\frac15\sin\frac{5\pi}{2}+\cdots\big)=\frac4\pi\big(1-\frac13+\frac15-\cdots\big)=1$, so the sum is $\pi/4$.
>
> (b) The triangle wave of Worked Example 15.11 is continuous at $t=0$ with value $A$: $\frac{8A}{\pi^2}\sum_{\text{odd }n}\frac1{n^2}=A$, so $\sum\frac{1}{(2k-1)^2}=\frac{\pi^2}{8}=1.2337$.
>
> (c) Every term contains $\sin n\pi=0$, so the series sums to $0$. The sawtooth jumps from $\pi$ (left limit) to $-\pi$ (right limit) at $t=\pi$; the midpoint is $0$. ✓ This is not an approximation failure: the series *converges*, to the midpoint.
>
> **Check.** $4\sum_{k=0}^{10^6-1}\frac{(-1)^k}{2k+1}=3.1415917$ (error $10^{-6}$, the slow $1/N$ convergence of an alternating series). The first three terms of (b) give $1+\frac19+\frac1{25}=1.1511$; the exact $\pi^2/8=1.2337$. ✓

### The Gibbs phenomenon

Near a jump, the partial sums of a Fourier series overshoot. As $N$ grows, the overshoot moves closer to the jump and becomes narrower, but its height does **not** decrease. This is the **Gibbs phenomenon** (observed by Michelson's harmonic analyser in 1898, explained by Gibbs in 1899, and in fact found by Wilbraham in 1848).

**Derivation for the square wave.** Take $A=1$, $\omega_0=1$ and let $S_M(t)=\frac4\pi\sum_{k=1}^{M}\frac{\sin(2k-1)t}{2k-1}$ contain the first $M$ odd harmonics (highest harmonic $N=2M-1$). Differentiate and sum the cosines with the formula of Ch. 3, §3.6:

$$
S_M'(t)=\frac4\pi\sum_{k=1}^{M}\cos(2k-1)t=\frac{2}{\pi}\cdot\frac{\sin2Mt}{\sin t}.
$$

The first zero of $S_M'$ to the right of the jump at $t=0$ is at $t^*=\pi/(2M)$, the location of the first and largest maximum. Since $S_M(0)=0$,

$$
S_M(t^*)=\frac2\pi\int_0^{\pi/(2M)}\frac{\sin2Mt}{\sin t}\,dt .
$$

Substitute $u=2Mt$: $S_M(t^*)=\frac2\pi\int_0^{\pi}\frac{\sin u}{2M\sin(u/2M)}\,du$. As $M\to\infty$, $2M\sin(u/2M)\to u$, so

$$
\lim_{M\to\infty}S_M(t^*)=\frac2\pi\int_0^\pi\frac{\sin u}{u}\,du=\frac{2}{\pi}\operatorname{Si}(\pi)=1.178980,
$$

where $\operatorname{Si}(z)=\int_0^z\frac{\sin u}{u}du$ is the sine integral, $\operatorname{Si}(\pi)=1.851937$.

> [!theorem] Gibbs phenomenon
> Let $x$ be piecewise smooth with a jump of size $J=x(t_0^+)-x(t_0^-)$ at $t_0$. Then near $t_0$ the partial sums $S_N$ overshoot $x(t_0^+)$ (and undershoot $x(t_0^-)$) by an amount tending to
> $$
> J\Big(\frac{1}{\pi}\operatorname{Si}(\pi)-\frac12\Big)=0.08949\,J,
> $$
> i.e. about **8.95% of the jump**, on each side. The overshoot occurs at a distance of order $T_0/(2N)$ from the jump, which tends to zero; its height does not.

The general case reduces to the square wave: subtract from $x$ a multiple of a sawtooth or square wave with the same jump at $t_0$; the remainder is continuous near $t_0$ and its series converges uniformly there, so it contributes no overshoot.

Gibbs is fully consistent with the convergence theorems. At each *fixed* $t\ne t_0$ the overshoot eventually moves past $t$, so pointwise convergence holds; the overshoot region has width $O(1/N)$, so its contribution to the mean-square error tends to zero; but the *maximum* error never falls below about $0.09J$, so convergence is not uniform.

```desmos-graph
left=-0.1; right=1.2; bottom=-0.2; top=1.3
---
y=\operatorname{sign}(x)|black|dashed
y=\sum_{k=1}^{10}\frac{4}{\pi(2k-1)}\sin((2k-1)x)|red
y=\sum_{k=1}^{40}\frac{4}{\pi(2k-1)}\sin((2k-1)x)|blue
y=1.17898|grey|dotted
```
*Zoom on the jump of the unit square wave at $t=0$: partial sums with harmonics up to 19 (red) and up to 79 (blue). Increasing $N$ squeezes the overshoot towards the jump but both peaks reach the same height, about $1.179$ (dotted).*

> [!example] Worked Example 15.19 — Size and position of the Gibbs overshoot
> A 1 kHz, $\pm5$ V square wave is passed through an ideal (brick-wall) low-pass filter that keeps harmonics up to 19 kHz and removes the rest. Find the peak output voltage and the time after each edge at which it occurs. Repeat for a 79 kHz cutoff.
>
> **Solution.** The filter output is the partial sum with harmonics $1,3,\dots,19$, i.e. $M=10$ odd harmonics. The peak is at $\theta^*=\pi/(2M)=\pi/20$ after the edge, i.e. $t^*=\frac{\theta^*}{\omega_0}=\frac{\pi/20}{2\pi\times1000}=25\ \mu$s. Its value for finite $M$: numerically $S_{10}(\pi/20)=1.17981$, so the peak is $5\times1.1798=5.90$ V — an overshoot of $0.90$ V, which is $9.0\%$ of the $10$ V jump.
>
> With harmonics to 79 kHz ($M=40$): $t^*=\frac{1}{4\cdot40\cdot1000}=6.25\ \mu$s, and $S_{40}(t^*)=1.17903$, peak $5.895$ V. The overshoot has moved four times closer to the edge but is essentially the same size; the limit is $5\times1.17898=5.8949$ V.
>
> **Check.** A brute-force search of $S_M$ on a fine grid gives maxima $1.17981$ at $t=0.15708$ ($M=10$) and $1.17903$ at $t=0.03927$ ($M=40$), matching $\pi/20$ and $\pi/80$. ✓

### Taming Gibbs: smoothing the partial sums

When a truncated series must be *used* — to synthesise a waveform, to design an FIR filter by truncating an ideal response, to display a reconstructed signal — the overshoot is usually unwanted. It comes from the abrupt truncation of the coefficients at $n=N$; tapering them removes it at the cost of a slower edge.

- **Fejér (Cesàro) means**: multiply $c_n$ by $1-\frac{\lvert n\rvert}{N+1}$. Since the Fejér kernel is non-negative, $\sigma_N$ never leaves the range $[\min x,\max x]$: no overshoot at all, but a noticeably rounded edge.
- **Lanczos $\sigma$-factors**: multiply $c_n$ by $\operatorname{sinc}\big(n/(N+1)\big)$. Equivalent to averaging $S_N$ over a window one period of the highest harmonic wide; it reduces the overshoot to about 1% of the jump while keeping a reasonably sharp edge.

> [!example] Worked Example 15.20 — Smoothed partial sums of a square wave
> For the unit square wave with harmonics up to $N=19$, compare the maximum of the raw partial sum with the Lanczos-smoothed and Fejér-smoothed sums.
>
> **Solution.** The weights on harmonic $n$ are $1$, $\operatorname{sinc}(n/20)$ and $1-n/20$ respectively; e.g. for $n=19$ they are $1$, $0.0524$ and $0.05$. Evaluating on a fine grid over a half-period, and measuring the 10%–90% rise time from the edge:
>
> | Partial sum | Maximum | Overshoot | 10–90% rise (rad) |
> |---|---|---|---|
> | raw $S_{19}$ | 1.1798 | 18.0% of $A$ (9.0% of jump) | 0.074 |
> | Lanczos $\sigma$ | 1.0235 | 2.3% of $A$ | 0.143 |
> | Fejér | 0.9682 (at $t=\pi/2$) | none | 0.263 |
>
> The Fejér sum stays below 1 everywhere, as the positivity of its kernel guarantees, but at the price of a rise time three and a half times that of the raw sum; the Lanczos sum is the usual compromise.

```desmos-graph
left=-0.1; right=1.6; bottom=-0.2; top=1.3
---
y=\operatorname{sign}(x)|black|dashed
y=\sum_{k=1}^{10}\frac{4}{\pi(2k-1)}\sin((2k-1)x)|red
y=\sum_{k=1}^{10}\frac{4}{\pi(2k-1)}\frac{\sin(\pi(2k-1)/20)}{\pi(2k-1)/20}\sin((2k-1)x)|blue
y=\sum_{k=1}^{10}\frac{4}{\pi(2k-1)}(1-\frac{2k-1}{20})\sin((2k-1)x)|green
```
*Harmonics up to 19 of the unit square wave: raw partial sum (red) with the Gibbs overshoot, Lanczos-smoothed (blue) with almost no overshoot, Fejér-smoothed (green) with none but a slower rise.*

> [!tip] Engineering note — Gibbs in hardware
> Any system with a sharp cutoff in frequency rings at sharp edges in time. The step response of an ideal low-pass filter overshoots by 8.95% for the same mathematical reason (Ch. 16); steep anti-aliasing filters and brick-wall digital filters show this ringing on square-wave tests; and FIR filters designed by truncating an ideal impulse response have ripples of about 9% (−21 dB) near the band edge whatever their length. Tapered windows (Hann, Hamming, Kaiser) are the filter-design version of the $\sigma$-factors above: they trade transition sharpness for suppressed ripple.

> [!warning] Pitfall — more terms do not remove the overshoot
> It is tempting to blame a 9% overshoot in a reconstructed waveform on "too few harmonics". Adding harmonics only narrows it. If the overshoot matters, taper the coefficients (or, in hardware, use a filter with a gentler roll-off) — do not just raise $N$.

## 15.7 Fourier coefficients

### Motivation

The Euler integrals of §15.2 always work, but integrating by parts over several pieces is slow and error-prone. Experienced analysts rarely do it. Instead they (i) build new series from known ones using the *properties* of the coefficient map — shifting, scaling, differentiating, integrating, multiplying; (ii) use the **jump method**, which reads the coefficients of any piecewise-polynomial waveform straight off the sizes and positions of its discontinuities; and (iii) compute coefficients numerically with the FFT when no formula is available, knowing exactly what error that introduces. This section develops all three, together with the result that ties them together: **the smoothness of a waveform determines how fast its coefficients decay**.

### Properties of the coefficients

Write $x\leftrightarrow c_n$ and $y\leftrightarrow d_n$ for two signals with the same period $T_0$.

> [!theorem] Fourier-coefficient properties
> | Operation | Signal | Coefficients |
> |---|---|---|
> | linearity | $\alpha x+\beta y$ | $\alpha c_n+\beta d_n$ |
> | time shift | $x(t-t_0)$ | $e^{-jn\omega_0t_0}c_n$ |
> | frequency shift (modulation by a harmonic) | $e^{jm\omega_0t}x(t)$ | $c_{n-m}$ |
> | time reversal | $x(-t)$ | $c_{-n}$ |
> | conjugation | $x^*(t)$ | $c_{-n}^*$ |
> | time scaling ($a>0$; period $T_0/a$) | $x(at)$ | $c_n$ (at frequencies $na\omega_0$) |
> | differentiation | $x'(t)$ | $jn\omega_0\,c_n$ |
> | integration (requires $c_0=0$) | $\int^tx(\tau)\,d\tau$ | $\dfrac{c_n}{jn\omega_0}$ ($n\ne0$), plus a constant |
> | periodic convolution | $\dfrac1{T_0}\displaystyle\int_{T_0}x(\tau)y(t-\tau)\,d\tau$ | $c_nd_n$ |
> | multiplication | $x(t)\,y(t)$ | $\displaystyle\sum_{m=-\infty}^{\infty}c_md_{n-m}$ |

**Proofs.** The shift, modulation, reversal and conjugation rules are substitutions in the coefficient integral. *Differentiation:* if $x$ is continuous and periodic with piecewise-continuous $x'$, integrate by parts:

$$
\frac1{T_0}\int_{T_0}x'(t)e^{-jn\omega_0t}dt=\frac1{T_0}\Big[x(t)e^{-jn\omega_0t}\Big]_{t_1}^{t_1+T_0}+\frac{jn\omega_0}{T_0}\int_{T_0}x(t)e^{-jn\omega_0t}dt=jn\omega_0c_n,
$$

the boundary term vanishing by periodicity. *Integration* is differentiation read backwards; the condition $c_0=0$ is needed because the integral of a signal with non-zero mean grows linearly and is not periodic. *Periodic convolution:* substitute the series of $y$ and use orthogonality,

$$
\frac1{T_0}\int_{T_0}x(\tau)\sum_md_me^{jm\omega_0(t-\tau)}d\tau=\sum_md_me^{jm\omega_0t}\,\frac1{T_0}\int_{T_0}x(\tau)e^{-jm\omega_0\tau}d\tau=\sum_mc_md_me^{jm\omega_0t}.
$$

*Multiplication:* the product of $\sum_mc_me^{jm\omega_0t}$ and $\sum_kd_ke^{jk\omega_0t}$ collects, at frequency $n\omega_0$, all pairs with $m+k=n$. $\square$

The last two rules are the duality at the heart of Fourier analysis: **convolution in time ↔ multiplication of coefficients; multiplication in time ↔ (discrete) convolution of coefficients.** Chapter 16 states the same pair for the transform. The multiplication rule explains why a nonlinearity spreads the spectrum: squaring a signal convolves its line spectrum with itself, producing lines at all sums and differences of the original frequencies (Ch. 3, §3.7).

> [!theorem] Term-by-term differentiation and integration
> (a) **Integration is always allowed.** If $x$ is piecewise continuous, the Fourier series of $x$ — whether or not it converges — may be integrated term by term between any limits, and the resulting series converges (uniformly) to the integral of $x$.
>
> (b) **Differentiation needs care.** If $x$ is continuous, periodic and piecewise smooth (in particular, no jumps, *including at the period boundary*), then the term-by-term derivative of its series is the Fourier series of $x'$, converging to $\frac12[x'(t^+)+x'(t^-)]$ wherever $x'$ satisfies the Dirichlet conditions.

Integration improves convergence (coefficients gain a factor $1/n$); differentiation worsens it (coefficients gain a factor $n$). Differentiating the square-wave series term by term gives $\frac{4A\omega_0}{\pi}\sum_{\text{odd}}\cos n\omega_0t$, which diverges at every $t$ — the derivative of a function with jumps contains impulses, which no ordinary function represents. The jump method below makes this precise.

> [!example] Worked Example 15.21 — Triangle wave from the square wave by integration
> Obtain the series of the triangle wave of Worked Example 15.11 by integrating the square wave, and confirm the result.
>
> **Solution.** Let $s(t)$ be the square wave of amplitude $A$ (Worked Example 15.3), $s(t)=\frac{4A}{\pi}\sum_{\text{odd }n}\frac{\sin n\omega_0t}{n}$. Its integral from $0$ is a triangle: $\int_0^ts=At$ for $0<t<T_0/2$, rising from $0$ to $AT_0/2$, then falling back — a triangle of peak-to-peak $AT_0/2$ with mean $AT_0/4$. Integrating term by term:
> $$
> \int_0^ts(\tau)\,d\tau=\frac{4A}{\pi\omega_0}\sum_{\text{odd }n}\frac{1-\cos n\omega_0t}{n^2}=\frac{4A}{\pi\omega_0}\cdot\frac{\pi^2}{8}-\frac{4A}{\pi\omega_0}\sum_{\text{odd }n}\frac{\cos n\omega_0t}{n^2},
> $$
> using $\sum_{\text{odd}}1/n^2=\pi^2/8$ (Worked Example 15.18). The constant is $\frac{A\pi}{2\omega_0}=\frac{AT_0}{4}$, the mean. ✓ The triangle of Worked Example 15.11 has peak $B$ at $t=0$ and is $B-\frac{4B}{AT_0}\int_0^ts$ (scaled and inverted so that it falls from $B$ to $-B$). Its cosine coefficients are therefore $\frac{4B}{AT_0}\cdot\frac{4A}{\pi\omega_0n^2}=\frac{16B}{\pi T_0\omega_0n^2}=\frac{8B}{\pi^2n^2}$ for odd $n$. ✓

### Smoothness and the rate of decay

> [!theorem] Decay of Fourier coefficients
> Let $x$ be $T_0$-periodic and let $k\ge0$ be an integer.
> - If $x$ has a jump discontinuity, then $c_n$ decays like $1/n$ and no faster (as $\lvert n\rvert\to\infty$, typically $\lvert c_n\rvert\approx C/\lvert n\rvert$ along a subsequence).
> - If $x,x',\dots,x^{(k-1)}$ are continuous (as periodic functions) and $x^{(k)}$ is piecewise smooth with jumps, then $\lvert c_n\rvert=O\big(1/\lvert n\rvert^{k+1}\big)$, with the leading term set by the jumps of $x^{(k)}$.
> - If $x$ has continuous derivatives of every order, $c_n$ decays faster than any power of $1/n$; if $x$ is analytic (Ch. 18), $c_n$ decays exponentially.

**Proof.** Apply the differentiation rule $k$ times: $c_n=\frac{\text{coefficient of }x^{(k)}}{(jn\omega_0)^{k}}$, and the coefficients of the piecewise-smooth $x^{(k)}$ with jumps are $O(1/n)$ by one more integration by parts (the boundary terms at the jumps give exactly $1/n$; the remaining integral is $o(1/n)$ by Riemann–Lebesgue). $\square$

| Waveform | Smoothness | Decay of $\lvert c_n\rvert$ | Coefficients |
|---|---|---|---|
| square, sawtooth, pulse train | jumps in $x$ | $1/n$ | square: $b_n=4A/(n\pi)$ |
| triangle, trapezoid, rectified sine | $x$ continuous, jumps in $x'$ | $1/n^2$ | triangle: $a_n=8A/(\pi^2n^2)$ |
| parabolic wave $t^2$ on $(-\pi,\pi)$ | $x$ continuous, jumps in $x'$ at $\pm\pi$ | $1/n^2$ | $a_n=4(-1)^n/n^2$ |
| cubic wave $t(\pi^2-t^2)$ on $(-\pi,\pi)$ | $x,x'$ continuous, jumps in $x''$ | $1/n^3$ | $b_n=12(-1)^{n+1}/n^3$ |
| $e^{\cos t}$ | analytic | exponential | $c_n=I_n(1)$ (modified Bessel) |

(The parabolic wave is continuous at $\pm\pi$ but its slope jumps from $2\pi$ to $-2\pi$ there. The cubic wave vanishes at $\pm\pi$ and has slope $-2\pi^2$ at both ends, so it and its derivative are continuous; its second derivative $-6t$ jumps from $-6\pi$ to $6\pi$.)

This theorem has direct engineering consequences. **Edges cost bandwidth.** A clock with ideal vertical edges has harmonics falling at only 20 dB/decade, so its spectrum extends far beyond the clock frequency and causes interference; rounding the edges (finite rise time) turns the decay into 40 dB/decade above a corner frequency of about $1/(\pi t_r)$. Conversely, the number of harmonics needed to reproduce a waveform to a given accuracy can be estimated before computing anything.

### The jump method

For piecewise-polynomial waveforms, the coefficients can be read off from the jumps. Suppose that over one period $x$ is piecewise smooth with jumps $J_i=x(t_i^+)-x(t_i^-)$ at points $t_i$. Then, in the generalised sense of Ch. 14 and Ch. 16, $x'$ is its ordinary (piecewise) derivative $x'_{\text{reg}}$ plus impulses:

$$
x'(t)=x'_{\text{reg}}(t)+\sum_iJ_i\,\delta_{T_0}(t-t_i),\qquad \delta_{T_0}(t)=\sum_{k}\delta(t-kT_0),
$$

and the periodic impulse train has coefficients $\frac1{T_0}\int_{T_0}\delta(t-t_i)e^{-jn\omega_0t}dt=\frac{1}{T_0}e^{-jn\omega_0t_i}$. Applying $c_n(x')=jn\omega_0c_n(x)$:

> [!theorem] Jump formula
> For $n\ne0$,
> $$
> c_n=\frac{1}{jn\omega_0}\Big[\frac1{T_0}\sum_iJ_i\,e^{-jn\omega_0t_i}+c_n(x'_{\text{reg}})\Big]=\frac{1}{jn\omega_0T_0}\sum_iJ_ie^{-jn\omega_0t_i}+\frac{1}{jn\omega_0}\,c_n(x'_{\text{reg}}).
> $$
> Apply the same formula to $x'_{\text{reg}}$, using the jumps $J_i'$ of the derivative, and so on. For a piecewise-polynomial waveform the process terminates, because some derivative is piecewise constant and its regular derivative is zero.

The same result follows without impulses by integrating by parts piece by piece; the impulses simply keep the bookkeeping clean.

> [!example] Worked Example 15.22 — Square wave and sawtooth by the jump method
> (a) The square wave of Worked Example 15.3 has jumps $+2A$ at $t=0$ and $-2A$ at $t=T_0/2$, and is constant in between ($x'_{\text{reg}}=0$). So
> $$
> c_n=\frac{1}{jn\omega_0T_0}\big(2A-2Ae^{-jn\pi}\big)=\frac{2A\,(1-(-1)^n)}{j2\pi n}=\frac{2A}{j\pi n}\ (n\text{ odd}),\ 0\ (n\text{ even}),
> $$
> and $b_n=-2\operatorname{Im}c_n=\frac{4A}{\pi n}$. ✓
>
> (b) The sawtooth $x=t$ on $(-\pi,\pi)$ has one jump, $-2\pi$ at $t=\pi$, and $x'_{\text{reg}}=1$, whose coefficients vanish for $n\ne0$. So $c_n=\frac{1}{jn\cdot2\pi}(-2\pi)e^{-jn\pi}=\frac{-(-1)^n}{jn}=\frac{j(-1)^n}{n}$ — Worked Example 15.10(b) in one line. ✓

> [!example] Worked Example 15.23 — Trapezoidal wave
> A trapezoidal voltage (a square wave with finite rise time) rises linearly from $0$ to $A$ over $0<\theta<\alpha$, stays at $A$ until $\theta=\pi-\alpha$, falls linearly to $0$ at $\theta=\pi$, and repeats with opposite sign on $(\pi,2\pi)$ ($\theta=\omega_0t$, $0<\alpha\le\pi/2$). Find its coefficients by the jump method and evaluate them for $\alpha=\pi/6$ ($30^\circ$ ramps).
>
> **Solution.** $x$ is continuous, so it has no jumps. Its derivative (with respect to $\theta$, so $\omega_0=1$, $T_0=2\pi$) is piecewise constant: $+A/\alpha$ on $(-\alpha,\alpha)$, $-A/\alpha$ on $(\pi-\alpha,\pi+\alpha)$, and $0$ elsewhere. The jumps of $x'$ are
>
> | position | $-\alpha$ | $\alpha$ | $\pi-\alpha$ | $\pi+\alpha$ |
> |---|---|---|---|---|
> | jump $J'_i$ | $+A/\alpha$ | $-A/\alpha$ | $-A/\alpha$ | $+A/\alpha$ |
>
> and $x''_{\text{reg}}=0$. Applying the jump formula twice ($x$ has no jumps, so only the second level contributes):
> $$
> c_n=\frac{1}{jn}\cdot\frac{1}{jn\,2\pi}\sum_iJ'_ie^{-jn\theta_i}=-\frac{A}{2\pi\alpha n^2}\Big[e^{jn\alpha}-e^{-jn\alpha}-e^{-jn(\pi-\alpha)}+e^{-jn(\pi+\alpha)}\Big].
> $$
> The bracket is $2j\sin n\alpha-(-1)^n\big(e^{jn\alpha}-e^{-jn\alpha}\big)=2j\sin n\alpha\,\big(1-(-1)^n\big)$. So $c_n=0$ for even $n$ (half-wave symmetry ✓) and, for odd $n$, $c_n=-\frac{2jA\sin n\alpha}{\pi\alpha n^2}$, giving
> $$
> b_n=-2\operatorname{Im}c_n=\frac{4A}{\pi\alpha}\cdot\frac{\sin n\alpha}{n^2}\quad(n\text{ odd}),\qquad a_n=0 .
> $$
> With $\alpha=\pi/6$: $b_1=1.2159A$, $b_3=0.2702A$, $b_5=0.0486A$, $b_7=-0.0248A$.
>
> **Checks.** (i) Numerical quadrature reproduces all four values to $10^{-12}$. ✓ (ii) As $\alpha\to0$, $\frac{\sin n\alpha}{\alpha}\to n$ and $b_n\to\frac{4A}{\pi n}$ — the square wave. ✓ (iii) The ramps make the waveform continuous, so $b_n=O(1/n^2)$: compared with the square wave ($b_7=0.1819A$) the 7th harmonic is 7 times smaller. Writing $b_n=\frac{4A}{\pi n}\operatorname{sinc}\!\big(\frac{n\alpha}{\pi}\big)$ shows the trapezoid's spectrum as the square-wave spectrum multiplied by a sinc envelope whose first null is at $n=\pi/\alpha$: the frequency-domain picture of "finite rise time costs less bandwidth".

```desmos-graph
left=-0.5; right=6.8; bottom=-1.4; top=1.4
---
y=\min(1,\max(-1,\frac{6}{\pi}\arcsin(\sin(x))))|black|dashed
y=\sum_{k=1}^{4}\frac{24}{\pi^2}\frac{\sin((2k-1)\pi/6)}{(2k-1)^2}\sin((2k-1)x)|blue
y=\sum_{k=1}^{4}\frac{4}{\pi(2k-1)}\sin((2k-1)x)|red
```
*The trapezoidal wave with $30^\circ$ ramps (dashed) and its partial sum through the 7th harmonic (blue), against the square-wave partial sum through the 7th harmonic (red). With the same four harmonics the trapezoid is reproduced closely because its coefficients fall as $1/n^2$.*

### Computing coefficients numerically: the DFT

When $x$ is known only by samples — a measured current, a simulation output — the coefficient integral is approximated by a sum. Take $N$ equally spaced samples $x_k=x(kT_0/N)$, $k=0,\dots,N-1$, and apply the rectangle (equivalently, for periodic integrands, trapezoidal) rule:

$$
\tilde c_n=\frac1N\sum_{k=0}^{N-1}x_k\,e^{-j2\pi nk/N}.
$$

This is the **discrete Fourier transform** divided by $N$, computed by the FFT in $O(N\log N)$ operations (Ch. 24). Its error is described exactly by the following identity.

> [!theorem] Aliasing of Fourier coefficients
> If the Fourier series of $x$ converges at the sample points (with the midpoint value used at any jump), then
> $$
> \tilde c_n=\sum_{m=-\infty}^{\infty}c_{n+mN}.
> $$
> In particular $\tilde c_n=c_n$ exactly for $\lvert n\rvert<N/2$ if $x$ is band-limited with $c_n=0$ for $\lvert n\rvert\ge N/2$.

**Proof.** Substitute $x_k=\sum_mc_me^{j2\pi mk/N}$ and use $\frac1N\sum_{k=0}^{N-1}e^{j2\pi(m-n)k/N}=1$ if $m\equiv n\pmod N$ and $0$ otherwise (a finite geometric series). $\square$

Each computed coefficient is contaminated by every true coefficient whose index differs from it by a multiple of $N$ — the coefficient-domain version of aliasing (Ch. 16, §16.8). The error therefore depends on how fast $c_n$ decays: for a smooth periodic waveform it is tiny even for modest $N$ (the trapezoidal rule is *spectrally accurate* for smooth periodic integrands), while for a waveform with jumps it is in general $O(1/N)$ (improving to $O(1/N^2)$ for a fixed harmonic when every jump falls on a sample instant and the midpoint value is used there).

> [!example] Worked Example 15.24 — FFT estimate of sawtooth coefficients
> Sample the sawtooth $x=t$ on $(-\pi,\pi)$ at $N=16$ points $t_k=-\pi+2\pi k/N$ (using the midpoint value $0$ at the jump, $k=0$), compute $\tilde c_1$ and $\tilde c_4$, and compare with the exact $c_n=j(-1)^n/n$. Repeat with $N=64$.
>
> **Solution.** Exact: $c_1=-j$, $c_4=0.25j$. The FFT (with the phase factor $e^{jn\pi}$ that accounts for the samples starting at $-\pi$ instead of $0$) gives
>
> | $N$ | $\tilde c_1$ | error | $\tilde c_4$ | error |
> |---|---|---|---|---|
> | 16 | $-0.98712j$ | 1.3% | $0.19635j$ | 21% |
> | 64 | $-0.99920j$ | 0.08% | $0.24678j$ | 1.3% |
>
> **Explanation by the aliasing identity.** For $N=16$, $\tilde c_4=\sum_mc_{4+16m}=j\sum_m\frac{1}{4+16m}=\frac{j}{16}\sum_m\frac{1}{m+\frac14}=\frac{j}{16}\,\pi\cot\frac\pi4=\frac{j\pi}{16}=0.19635j$, using the partial-fraction expansion $\pi\cot\pi z=\sum_m\frac{1}{z+m}$ (summed symmetrically; Ch. 18). ✓ For $N=64$ the same formula gives $\frac{\pi}{64}\cot\frac{\pi}{16}=0.24678$. ✓ The harmonic at a quarter of the sampling rate is badly aliased because the $1/n$ tail of the sawtooth is heavy, and the error of a fixed harmonic falls only like $1/N^2$ (here $21\%\to1.3\%$ for a fourfold increase in $N$). A practical rule: to measure harmonics up to $n_{\max}$ to within a few per cent when the waveform has jumps, sample at many times $2n_{\max}$ per period, or low-pass filter before sampling.

> [!warning] Pitfall — coefficients of a waveform with jumps need an explicit midpoint
> When sampling a waveform that has a jump exactly at a sample instant, use the midpoint value. Using either one-sided value adds a spurious constant $\pm J/(2N)$ to *every* $\tilde c_n$ — an error that does not decrease with $n$ and that looks like a broadband noise floor.

> [!example] Worked Example 15.25 — Full-wave rectified sine by the multiplication rule
> Obtain the Fourier series of $\lvert\sin t\rvert$ from those of $\sin t$ and the square wave $s(t)=\operatorname{sgn}(\sin t)$, using $\lvert\sin t\rvert=\sin t\cdot s(t)$.
>
> **Solution.** $\sin t$ has $d_{1}=\frac{1}{2j}$, $d_{-1}=-\frac1{2j}$. The square wave ($A=1$) has $c_m=\frac{2}{j\pi m}$ for odd $m$, zero for even $m$. By the multiplication rule the product has
> $$
> e_n=\sum_mc_md_{n-m}=c_{n-1}d_1+c_{n+1}d_{-1}=\frac{1}{2j}\big(c_{n-1}-c_{n+1}\big).
> $$
> For odd $n$, $n\pm1$ is even, so $e_n=0$. For even $n$, $e_n=\frac{1}{2j}\cdot\frac{2}{j\pi}\Big(\frac1{n-1}-\frac1{n+1}\Big)=-\frac{1}{\pi}\cdot\frac{2}{n^2-1}$. In particular $e_0=\frac{2}{\pi}$. Then $a_n=2e_n$:
> $$
> \lvert\sin t\rvert=\frac{2}{\pi}-\frac{4}{\pi}\sum_{k=1}^{\infty}\frac{\cos2kt}{4k^2-1}.
> $$
> **Check.** This is twice the even-harmonic part of Worked Example 15.5, as $\lvert\sin t\rvert=2\times$(half-wave rectified) $-\sin t$ requires. ✓ At $t=0$: $\sum_k\frac{1}{4k^2-1}=\frac12$ (telescoping), so the series gives $\frac2\pi-\frac4\pi\cdot\frac12=0=\lvert\sin0\rvert$. ✓

> [!tip] Engineering note — estimating EMI from edge rates
> The jump formula applied to a trapezoidal clock of period $T_0$, duty 50% and rise time $t_r$ gives the familiar EMC envelope: flat to $f_0$, $-20$ dB/decade from $f_0$ to about $1/(\pi\tau)$ (with $\tau$ the pulse width), and $-40$ dB/decade above $1/(\pi t_r)$. Slowing edges from 1 ns to 5 ns moves the second corner from 318 MHz to 64 MHz and reduces the harmonic amplitudes above 318 MHz by a factor 5 (14 dB), with a smaller reduction between the two corners — without changing the logic function.

## 15.8 Power and energy relationships

### Motivation

The most frequently used result about Fourier series is not a formula for coefficients but a statement about *power*: the average power of a periodic signal equals the sum of the powers of its harmonics. It underlies every RMS meter calibrated for distorted waveforms, every harmonic-distortion specification, the definition of power factor in non-sinusoidal circuits, and the estimate of how much bandwidth a signal really needs. Chapter 3 (§3.7) proved it for a finite sum of harmonics; here we prove it for every finite-power periodic signal, extend it to products of two different signals (voltage and current), and use it to define and compute distortion measures.

### Power signals

A periodic signal does not decay, so its total energy $\int_{-\infty}^\infty\lvert x\rvert^2dt$ is infinite. The meaningful quantities are the **energy per period** $E_{T_0}=\int_{T_0}\lvert x\rvert^2dt$ and the **average (normalised) power**

$$
P_x=\frac1{T_0}\int_{T_0}\lvert x(t)\rvert^2dt=X_{\text{rms}}^2 ,
$$

the power the signal would dissipate in a $1\ \Omega$ resistor if it were a voltage or a current. Periodic signals are **power signals**; the pulses and transients of Ch. 16, with finite total energy, are **energy signals**. For a voltage $v$ across a resistance $R$ the actual power is $V_{\text{rms}}^2/R$.

### Parseval's theorem

> [!theorem] Parseval's theorem for Fourier series
> If $x$ is $T_0$-periodic with $\int_{T_0}\lvert x\rvert^2dt<\infty$, then
> $$
> \frac1{T_0}\int_{T_0}\lvert x(t)\rvert^2dt=\sum_{n=-\infty}^{\infty}\lvert c_n\rvert^2 .
> $$
> For real $x$ in trigonometric and amplitude–phase form,
> $$
> \frac1{T_0}\int_{T_0}x^2(t)\,dt=\Big(\frac{a_0}{2}\Big)^2+\frac12\sum_{n=1}^\infty\big(a_n^2+b_n^2\big)=A_0^2+\sum_{n=1}^\infty\frac{A_n^2}{2}=A_0^2+\sum_{n=1}^{\infty}X_{n,\text{rms}}^2 ,
> $$
> where $X_{n,\text{rms}}=A_n/\sqrt2$ is the RMS value of the $n$th harmonic.

**Proof.** By the best-approximation theorem (§15.1), $\lVert x-S_N\rVert^2=\lVert x\rVert^2-\sum_{\lvert n\rvert\le N}\lvert c_n\rvert^2$. By the mean-square convergence theorem (§15.6), the left side tends to $0$ as $N\to\infty$. The trigonometric form follows from $\lvert c_n\rvert^2+\lvert c_{-n}\rvert^2=2\lvert c_n\rvert^2=\frac{A_n^2}{2}=\frac{a_n^2+b_n^2}{2}$ for real $x$. $\square$

In words: **the mean-square value is the sum of the mean-square values of the harmonics, whatever their phases.** The harmonics are orthogonal "vectors" and Parseval's theorem is Pythagoras' theorem in the function space. It is a statement about power only: the peak value, the crest factor, and how hard an amplifier clips all depend on the phases.

The squared magnitudes $\lvert c_n\rvert^2$ form the **power spectrum**; the partial sums of Parseval give the **cumulative power**, and the error of a truncated series is

$$
\lVert x-S_N\rVert^2=\sum_{\lvert n\rvert>N}\lvert c_n\rvert^2,
$$

the power in the discarded harmonics. This is how one decides how many harmonics, or how much bandwidth, a waveform "needs".

> [!theorem] Generalised Parseval (power of a product)
> If $x\leftrightarrow c_n$ and $y\leftrightarrow d_n$ are $T_0$-periodic and square-integrable, then
> $$
> \frac1{T_0}\int_{T_0}x(t)\,y^*(t)\,dt=\sum_{n=-\infty}^\infty c_n\,d_n^* .
> $$
> For a real voltage $v=V_0+\sum_{n\ge1}V_n\cos(n\omega_0t+\alpha_n)$ and current $i=I_0+\sum_{n\ge1}I_n\cos(n\omega_0t+\beta_n)$ (peak amplitudes), the **average power** is
> $$
> P=\frac1{T_0}\int_{T_0}v\,i\,dt=V_0I_0+\sum_{n=1}^\infty\frac{V_nI_n}{2}\cos(\alpha_n-\beta_n)=V_0I_0+\sum_{n=1}^{\infty}V_{n,\text{rms}}I_{n,\text{rms}}\cos\theta_n .
> $$

**Proof.** Apply Parseval to $x+y$ and $x+jy$ and use polarisation, or argue directly: the cross terms $c_md_n^*$ with $m\ne n$ average to zero by orthogonality. For the power formula take $c_n=\frac{V_n}2e^{j\alpha_n}$, $d_n=\frac{I_n}{2}e^{j\beta_n}$ and add the $\pm n$ terms: $2\operatorname{Re}\{c_nd_n^*\}=\frac{V_nI_n}{2}\cos(\alpha_n-\beta_n)$. $\square$

**Only voltage and current harmonics of the same frequency transfer average power.** A sinusoidal supply voltage delivers power only through the fundamental of the current; the current harmonics circulate without net energy transfer, yet they contribute to the RMS current and hence to losses in every conductor they pass through.

### Distortion measures

> [!definition] Total harmonic distortion and related measures
> For a periodic signal with zero mean, harmonic RMS values $X_1,X_2,\dots$ and total RMS value $X=\sqrt{\sum_{n\ge1}X_n^2}$:
> - **THD (fundamental-referenced)**, the standard power-engineering definition (IEEE 519, IEC 61000):
> $$
> \mathrm{THD}_F=\frac{\sqrt{\sum_{n\ge2}X_n^2}}{X_1}=\sqrt{\frac{X^2}{X_1^2}-1};
> $$
> - **THD (RMS-referenced)**, used by some audio analysers ("THD-R", or THD+N when noise is included):
> $$
> \mathrm{THD}_R=\frac{\sqrt{\sum_{n\ge2}X_n^2}}{X}=\frac{\mathrm{THD}_F}{\sqrt{1+\mathrm{THD}_F^2}};
> $$
> - **distortion factor** (fundamental fraction) $\;X_1/X=1/\sqrt{1+\mathrm{THD}_F^2}$;
> - **crest factor** $\;X_{\text{peak}}/X$ (equal to $\sqrt2$ for a sinusoid);
> - **form factor** $\;X/\overline{\lvert x\rvert}$ (RMS over rectified mean; $\pi/(2\sqrt2)=1.1107$ for a sinusoid).

$\mathrm{THD}_F$ can exceed 100%; $\mathrm{THD}_R$ cannot. For small distortion they agree. Parseval's theorem lets THD be computed from the total RMS value and the fundamental alone, without summing the infinite series of harmonics — the method used in all the examples below. A DC component, if present, is excluded from both numerator and denominator.

**Power factor with a sinusoidal voltage.** If $v$ is a pure sinusoid at $\omega_0$ and $i$ is distorted, then $P=V I_1\cos\theta_1$, with $\theta_1$ the phase angle between $v$ and the fundamental of $i$, while the apparent power is $S=VI$ (RMS values). Hence

$$
\mathrm{PF}=\frac PS=\underbrace{\frac{I_1}{I}}_{\text{distortion factor}}\times\underbrace{\cos\theta_1}_{\text{displacement factor}}=\frac{\cos\theta_1}{\sqrt{1+\mathrm{THD}_i^2}} .
$$

A load can have unity displacement factor and still a poor power factor.

> [!example] Worked Example 15.26 — Summing series with Parseval
> Use Parseval's theorem to evaluate $\sum_{n\ge1}\frac{1}{n^2}$, $\sum_{n\ge1}\frac1{n^4}$ and $\sum_{k\ge1}\frac{1}{(4k^2-1)^2}$.
>
> **Solution.** (a) Sawtooth $x=t$ on $(-\pi,\pi)$: $b_n=\frac{2(-1)^{n+1}}{n}$, mean square $\frac{\pi^2}{3}$. Parseval: $\frac12\sum\frac{4}{n^2}=\frac{\pi^2}{3}$, so $\sum\frac1{n^2}=\frac{\pi^2}{6}=1.644934$.
>
> (b) Triangle of peak 1: $a_n=\frac{8}{\pi^2n^2}$ (odd $n$), mean square $\frac{1}{3}$ (a linear ramp from $-1$ to $1$). Parseval: $\frac12\sum_{\text{odd}}\frac{64}{\pi^4n^4}=\frac13$, so $\sum_{\text{odd}}\frac1{n^4}=\frac{\pi^4}{96}$. The even terms contribute $\sum\frac{1}{(2m)^4}=\frac{1}{16}\sum\frac1{m^4}$, so $\sum\frac{1}{n^4}\big(1-\frac1{16}\big)=\frac{\pi^4}{96}$ and $\sum\frac1{n^4}=\frac{\pi^4}{90}=1.082323$.
>
> (c) Full-wave rectified sine $\lvert\sin t\rvert$: mean square $\frac12$, DC $\frac2\pi$, $a_{2k}=-\frac{4}{\pi(4k^2-1)}$. Parseval: $\frac{4}{\pi^2}+\frac12\sum_k\frac{16}{\pi^2(4k^2-1)^2}=\frac12$, so $\sum_k\frac{1}{(4k^2-1)^2}=\frac{\pi^2}{16}\big(\frac12-\frac4{\pi^2}\big)=\frac{\pi^2-8}{16}=0.116850$.
>
> **Check.** Direct summation to $10^6$ terms: $1.644933$, $1.082323$, $0.116850$. ✓

> [!example] Worked Example 15.27 — THD of standard waveforms
> Compute $\mathrm{THD}_F$ for the square wave, triangle, sawtooth, six-step (120° quasi-square) waveform of Worked Example 15.12 and the $30^\circ$ trapezoid of Worked Example 15.23.
>
> **Solution.** In each case compute the total mean square directly from the waveform and the fundamental's mean square $A_1^2/2$, then $\mathrm{THD}_F=\sqrt{X^2/X_1^2-1}$.
>
> - *Square wave* ($\pm A$): $X^2=A^2$; $X_1^2=\frac12\big(\frac{4A}{\pi}\big)^2=\frac{8A^2}{\pi^2}$. $\mathrm{THD}=\sqrt{\frac{\pi^2}{8}-1}=48.34\%$.
> - *Triangle* (peak $A$): $X^2=\frac{A^2}3$; $X_1^2=\frac12\cdot\frac{64A^2}{\pi^4}$. $\mathrm{THD}=\sqrt{\frac{\pi^4}{96}-1}=12.12\%$.
> - *Sawtooth* ($t$ on $(-\pi,\pi)$): $X^2=\frac{\pi^2}{3}$; $X_1^2=\frac{2^2}{2}=2$. $\mathrm{THD}=\sqrt{\frac{\pi^2}{6}-1}=80.31\%$.
> - *Six-step*: the waveform is $\pm V$ for $\frac23$ of the period, so $X^2=\frac23V^2$; $b_1=\frac{2\sqrt3}{\pi}V$, $X_1^2=\frac{6}{\pi^2}V^2$. $\mathrm{THD}=\sqrt{\frac{2\pi^2}{18}-1}=\sqrt{\frac{\pi^2}{9}-1}=31.08\%$.
> - *Trapezoid* ($\alpha=\pi/6$, peak $A$): each ramp of duration $\alpha$ contributes mean square $\frac13A^2$, so over a half-period $X^2=\frac{(\pi-2\alpha)A^2+2\alpha A^2/3}{\pi}=\big(1-\frac{4\alpha}{3\pi}\big)A^2=0.7778A^2$; $X_1^2=\frac12(1.2159A)^2=0.7392A^2$. $\mathrm{THD}=22.86\%$.
>
> | Waveform | $X/A$ | $X_1/A$ | $\mathrm{THD}_F$ | $\mathrm{THD}_R$ | Crest factor |
> |---|---|---|---|---|---|
> | sine | 0.7071 | 0.7071 | 0 | 0 | 1.414 |
> | triangle | 0.5774 | 0.5732 | 12.1% | 12.0% | 1.732 |
> | trapezoid (30°) | 0.8819 | 0.8597 | 22.9% | 22.3% | 1.134 |
> | six-step (120°) | 0.8165 | 0.7797 | 31.1% | 29.7% | 1.225 |
> | square | 1 | 0.9003 | 48.3% | 43.5% | 1 |
> | sawtooth ($A=\pi$) | 0.5774 | 0.4502 | 80.3% | 62.6% | 1.732 |
>
> **Check.** Summing $b_n^2/2$ numerically for the six-step wave over $n\le 20\,001$ and comparing with $\frac23V^2$ reproduces $31.08\%$. ✓ Note how THD tracks smoothness (§15.7): jumps give large THD, kinks small THD.

> [!example] Worked Example 15.28 — Power in a circuit with harmonics
> A supply voltage $v(t)=325\cos\omega_0t+30\cos3\omega_0t$ V drives a nonlinear load drawing $i(t)=10\cos(\omega_0t-30^\circ)+4\cos(3\omega_0t-60^\circ)+2\cos(5\omega_0t+20^\circ)$ A. Find the RMS voltage and current, the average power, the apparent power, the power factor, and the voltage and current THD.
>
> **Solution.** RMS values by Parseval:
> $V=\sqrt{\frac{325^2+30^2}{2}}=230.79$ V, $\quad I=\sqrt{\frac{10^2+4^2+2^2}{2}}=\sqrt{60}=7.746$ A.
>
> Average power — only matching harmonics contribute; the 5th harmonic current has no 5th harmonic voltage to work against:
> $$
> P=\frac{325\times10}{2}\cos30^\circ+\frac{30\times4}{2}\cos60^\circ+0=1407.3+30.0=1437.3\ \text{W}.
> $$
> Apparent power $S=VI=1787.7$ VA; power factor $P/S=0.804$.
>
> $\mathrm{THD}_v=30/325=9.23\%$; $\mathrm{THD}_i=\sqrt{4^2+2^2}/10=44.7\%$.
>
> **Interpretation.** The fundamental displacement factor alone is $\cos30^\circ=0.866$; harmonics pull the true power factor down to $0.804$. The third-harmonic product contributes a small amount of *real* power (30 W) because the supply itself is distorted.
>
> **Check.** Averaging $v(t)\,i(t)$, $v^2$ and $i^2$ numerically over $10^6$ points per period gives $1437.29$ W, $230.787$ V and $7.7460$ A. ✓

> [!example] Worked Example 15.29 — Power factor of a rectifier with a large DC inductor
> A single-phase bridge rectifier with a large DC-side inductor draws a square-wave current of amplitude $I_d$ in phase with the sinusoidal supply voltage $V\sqrt2\sin\omega_0t$. Find the power factor.
>
> **Solution.** The current is a square wave $\pm I_d$, so $I=I_d$ and $I_1=\frac{4I_d}{\pi\sqrt2}=0.9003I_d$, in phase with the voltage ($\cos\theta_1=1$). Hence
> $$
> \mathrm{PF}=\frac{I_1}{I}\cos\theta_1=\frac{2\sqrt2}{\pi}=0.9003,\qquad \mathrm{THD}_i=48.3\%.
> $$
> The load is "resistive" as far as the fundamental is concerned, yet about 10% of the conductor capacity is wasted on harmonic current. Active power-factor-correction front ends shape the current into a sinusoid to push PF above 0.99.

> [!example] Worked Example 15.30 — How much bandwidth does a pulse train need?
> For the 5 V, 25%-duty pulse train of Worked Example 15.7, what fraction of the total power lies in the harmonics inside the first spectral null ($\lvert f\rvert<1/\tau$)? How much power lies in the DC term? For the $\pm A$ square wave, what fraction lies in the fundamental, in harmonics up to the 5th, and up to the 19th?
>
> **Solution.** Total power of the pulse train: the signal is $5$ V for a quarter of the time, so $P=25\times\frac14=6.25\ \text{V}^2$. Inside the null, $\lvert n\rvert\le3$:
> $$
> \sum_{n=-3}^{3}\lvert c_n\rvert^2=1.25^2+2\big(1.1254^2+0.7958^2+0.3751^2\big)=1.5625+2(1.2665+0.6333+0.1407)=5.6435\ \text{V}^2 ,
> $$
> which is $90.3\%$ of the total. The DC term alone carries $1.5625/6.25=25\%$.
>
> For the square wave, the fraction in the odd harmonics up to $N$ is $\frac{8}{\pi^2}\sum_{\text{odd }n\le N}\frac{1}{n^2}$: fundamental $81.06\%$; up to the 5th $93.31\%$; up to the 19th $97.98\%$.
>
> **Interpretation.** Power converges far faster than the waveform does pointwise: the 19-harmonic square-wave approximation of Worked Example 15.19 has a 9% overshoot at every edge, yet it carries 98% of the power. "Bandwidth containing 90% of the power" is a standard practical definition of signal bandwidth.

> [!warning] Pitfall — adding RMS values linearly, or adding powers of the same frequency
> RMS values of *different* harmonics add in quadrature, never linearly: $10$ A at 50 Hz plus $3$ A at 150 Hz is $10.44$ A, not $13$ A. But two components at the *same* frequency (two sources' 5th-harmonic currents meeting at a bus, say) must first be combined as phasors — their cross term does not average to zero, and the result can be anywhere between the difference and the sum of their amplitudes.

> [!tip] Engineering note — "true RMS" meters
> An average-responding meter measures $\overline{\lvert x\rvert}$ and multiplies by the sinusoidal form factor $1.1107$. On a square wave ($\overline{\lvert x\rvert}=A$, true RMS $A$) it reads $1.11A$, 11% high; on the narrow, peaky current of a capacitor-input rectifier (crest factor 2.5–3) it typically reads 20–40% low. A true-RMS meter computes $\sqrt{\overline{x^2}}$ directly, which by Parseval is the quadrature sum of all harmonics within its bandwidth.

## 15.9 Applications to periodic signals

### Motivation

Section 15.1 promised that Fourier series would turn the response of a linear system to a periodic input into a list of multiplications. This section keeps that promise, states the conditions under which it is valid, and works through the standard situations: filtering a switching waveform, resonance with one harmonic of a non-sinusoidal drive, and ripple at the output of a rectifier. The applications are motivation; the method is general and applies equally to mechanical, thermal and acoustic systems described by linear differential equations with constant coefficients.

### Steady-state response of an LTI system

> [!theorem] Periodic steady state of an LTI system
> Let a stable LTI system with frequency response $H(j\omega)$ be driven by a $T_0$-periodic input $x(t)=\sum_nc_ne^{jn\omega_0t}$ with finite power. Then the output, after transients have decayed, is the $T_0$-periodic signal
> $$
> y_{\text{ss}}(t)=\sum_{n=-\infty}^{\infty}c_n\,H(jn\omega_0)\,e^{jn\omega_0t},
> $$
> i.e. the output coefficients are $d_n=c_nH(jn\omega_0)$. For real signals, in amplitude–phase form,
> $$
> y_{\text{ss}}(t)=A_0H(0)+\sum_{n=1}^{\infty}A_n\,\lvert H(jn\omega_0)\rvert\cos\big(n\omega_0t+\phi_n+\angle H(jn\omega_0)\big).
> $$
> The output power is $P_y=\sum_n\lvert c_n\rvert^2\lvert H(jn\omega_0)\rvert^2$.

**Justification.** Stability means $\int\lvert h\rvert<\infty$, so $H(j\omega)$ is bounded and continuous. The periodic steady state is $y_{\text{ss}}(t)=\int_{-\infty}^{\infty}h(\tau)x(t-\tau)\,d\tau$ (the input applied for all past time). Writing the infinite convolution as a periodic convolution with the periodised impulse response $h_{T_0}(t)=\sum_kh(t-kT_0)$, whose Fourier coefficients are $\frac1{T_0}H(jn\omega_0)$, the periodic-convolution property of §15.7 gives $d_n=T_0\cdot c_n\cdot\frac1{T_0}H(jn\omega_0)=c_nH(jn\omega_0)$. Because $H$ is bounded, $\sum\lvert d_n\rvert^2\le\max\lvert H\rvert^2\sum\lvert c_n\rvert^2<\infty$, so the output series converges in mean square; if $H(j\omega)$ falls off at high frequency, it converges faster than the input series. $\square$

The **complete** response to an input switched on at $t=0$ is $y_{\text{ss}}$ plus a transient determined by the system's poles and initial conditions; for a stable system the transient decays and $y\to y_{\text{ss}}$. The Laplace transform (Ch. 14) gives both together; the Fourier series gives the steady state directly, which is usually what is wanted. The procedure is:

1. Find the input's Fourier coefficients (§15.2–15.7), using symmetry and tables.
2. Find $H(j\omega)$ — from impedances, from the differential equation by $d/dt\to j\omega$, or from $H(s)$ at $s=j\omega$.
3. Multiply harmonic by harmonic: $d_n=c_nH(jn\omega_0)$.
4. Recombine, or — more often — stop at step 3 and read off amplitudes, power, THD or ripple directly.

For a linear constant-coefficient ODE $\sum_k\alpha_ky^{(k)}=\sum_m\beta_mx^{(m)}$ the frequency response is $H(j\omega)=\sum_m\beta_m(j\omega)^m/\sum_k\alpha_k(j\omega)^k$, because $d/dt$ acts on $e^{jn\omega_0t}$ as multiplication by $jn\omega_0$ (the differentiation rule of §15.7).

> [!example] Worked Example 15.31 — RC low-pass filter driven by a square wave
> A $\pm1$ V, 1 kHz square wave drives an RC low-pass filter with $RC=159.15\ \mu$s (corner frequency $f_c=1/(2\pi RC)=1$ kHz). Find the first four non-zero output harmonics, the output RMS value and THD, and the peak output voltage. Check the peak against a time-domain solution.
>
> **Solution.** *Input:* $x(t)=\frac{4}{\pi}\sum_{\text{odd}}\frac{\sin n\omega_0t}{n}$ (Worked Example 15.3), $\omega_0=2\pi\times1000$ rad/s.
>
> *Frequency response:* $H(j\omega)=\frac{1}{1+j\omega RC}$; at the harmonics $\omega_0RC=1$, so $H(jn\omega_0)=\frac{1}{1+jn}$, $\lvert H\rvert=\frac1{\sqrt{1+n^2}}$, $\angle H=-\arctan n$.
>
> | $n$ | input amplitude $\frac{4}{n\pi}$ (V) | $\lvert H\rvert$ | output amplitude (V) | phase shift |
> |---|---|---|---|---|
> | 1 | 1.2732 | 0.7071 | 0.9003 | $-45.0^\circ$ |
> | 3 | 0.4244 | 0.3162 | 0.1342 | $-71.6^\circ$ |
> | 5 | 0.2546 | 0.1961 | 0.0499 | $-78.7^\circ$ |
> | 7 | 0.1819 | 0.1414 | 0.0257 | $-81.9^\circ$ |
>
> So $y_{\text{ss}}(t)=0.9003\sin(\omega_0t-45^\circ)+0.1342\sin(3\omega_0t-71.6^\circ)+0.0499\sin(5\omega_0t-78.7^\circ)+\cdots$ V. For large $n$ the output amplitude is $\approx\frac{4}{\pi n^2}$: the filter converts the $1/n$ decay of the square wave into the $1/n^2$ decay of a waveform with kinks but no jumps — as it must, since a capacitor voltage cannot jump.
>
> *RMS and THD by Parseval:* $Y^2=\sum_{\text{odd}}\frac12\big(\frac{4}{\pi n}\big)^2\frac{1}{1+n^2}=0.41612$ V² (summing to $n=2\times10^5$), so $Y=0.6451$ V. The fundamental contributes $\frac12(0.9003)^2=0.40528$ V², hence $\mathrm{THD}_y=\sqrt{0.41612/0.40528-1}=16.4\%$, down from $48.3\%$ at the input.
>
> *Peak value:* the output peaks at the end of each positive half-cycle, $t=T_0/2$. Summing the series there to $10^5$ harmonics gives $0.91715$ V.
>
> **Time-domain check.** In the periodic steady state the capacitor charges towards $+1$ for half a period from $-V_p$ and must end at $+V_p$: $V_p=1-(1+V_p)e^{-T_0/(2RC)}$, so $V_p=\frac{1-e^{-T_0/2RC}}{1+e^{-T_0/2RC}}=\tanh\frac{T_0}{4RC}=\tanh\frac{\pi}{2}=0.917152$ V. ✓ The exact steady-state waveform $1-(1+V_p)e^{-t/RC}$ on the half-period also has RMS $0.64508$ V. ✓

```desmos-graph
left=-0.5; right=12.8; bottom=-1.4; top=1.4
---
y=\operatorname{sign}(\sin(x))|grey|dashed
y=\sum_{k=1}^{60}\frac{4}{\pi(2k-1)}\frac{\sin((2k-1)x)-(2k-1)\cos((2k-1)x)}{1+(2k-1)^2}|blue
y=1-(1+\tanh(\pi/2))e^{-x}|red|dotted|0<x<3.1416
```
*Steady-state RC filter output for a square-wave input with $\omega_0RC=1$ (horizontal axis $\omega_0t$): the Fourier sum through the 119th harmonic (blue) and the exact exponential charging curve on the first half-period (red, dotted) coincide.*

> [!example] Worked Example 15.32 — Resonance with a harmonic
> A lightly damped second-order system $\ddot y+0.2\dot y+25y=x(t)$ (natural frequency $\omega_n=5$ rad/s, $\zeta=0.02$) is driven by the unit square wave of period $2\pi$ s ($\omega_0=1$ rad/s). Find the steady-state response.
>
> **Solution.** $H(j\omega)=\frac{1}{25-\omega^2+0.2j\omega}$. With $b_n=\frac{4}{\pi n}$ for odd $n$:
>
> | $n$ | $b_n$ | $\lvert H(jn)\rvert$ | output amplitude | phase of $H$ |
> |---|---|---|---|---|
> | 1 | 1.2732 | 0.04167 | 0.05305 | $-0.5^\circ$ |
> | 3 | 0.4244 | 0.06246 | 0.02651 | $-2.1^\circ$ |
> | 5 | 0.2546 | 1.0000 | **0.25465** | $-90.0^\circ$ |
> | 7 | 0.1819 | 0.04160 | 0.00757 | $-176.7^\circ$ |
> | 9 | 0.1415 | 0.01786 | 0.00252 | $-178.2^\circ$ |
>
> At $n=5$ the denominator is $25-25+0.2j\cdot5=j$, so $H(j5)=-j$. The steady state is
> $$
> y_{\text{ss}}(t)\approx0.0531\sin t+0.0265\sin3t-0.2546\cos5t+0.0076\sin(7t-176.7^\circ)+\cdots .
> $$
> The 5th harmonic, only 20% of the input fundamental, is amplified by the resonance until it dominates the output: Parseval gives the output mean square $0.034214$, of which $94.8\%$ is in the 5th harmonic. The system responds at *five times* the drive frequency.
>
> **The undamped limit.** With no damping, $H(j5)=\frac1{25-25}$ is infinite: there is no periodic steady state, and the true solution contains a term growing like $t\cos5t$ (Ch. 12, resonance). In general, $\ddot y+\omega_n^2y=x(t)$ has a $T_0$-periodic solution if and only if the input has **no component at the resonant frequency**: either $\omega_n/\omega_0$ is not an integer, or $c_{\pm n}=0$ for $n=\omega_n/\omega_0$. Had the system been tuned to $\omega_n=4$, an even harmonic, the square wave (which has none) would drive it without resonance.

```desmos-graph
left=-0.5; right=12.8; bottom=-0.4; top=0.4
---
y=\sum_{k=1}^{30}\frac{4}{\pi(2k-1)}\frac{(25-(2k-1)^2)\sin((2k-1)x)-0.2(2k-1)\cos((2k-1)x)}{(25-(2k-1)^2)^2+0.04(2k-1)^2}|blue
y=0.25\operatorname{sign}(\sin(x))|grey|dashed
```
*Steady-state response of $\ddot y+0.2\dot y+25y=x(t)$ to a unit square wave of period $2\pi$ (blue; the input is shown scaled by $0.25$, dashed). The output oscillates five times per input period: the system has picked out the resonant 5th harmonic.*

> [!example] Worked Example 15.33 — Ripple at the output of a rectifier filter
> A full-wave rectifier fed from 230 V RMS, 50 Hz ($V_m=325$ V) supplies an LC filter: series inductor $L=0.5$ H, shunt capacitor $C=1000\ \mu$F, with a light load that can be neglected in computing the ripple. Assuming the rectifier output is the ideal $V_m\lvert\sin\omega t\rvert$, find the DC output and the ripple amplitudes.
>
> **Solution.** *Input series* (Worked Example 15.25), fundamental of the rectified wave at 100 Hz:
> $$
> v_r(t)=\frac{2V_m}{\pi}-\frac{4V_m}{\pi}\sum_{k\ge1}\frac{\cos(2k\omega t)}{4k^2-1}=206.90-137.93\cos(2\omega t)-27.59\cos(4\omega t)-11.82\cos(6\omega t)-\cdots\ \text{V}.
> $$
> *Filter* (unloaded LC divider): $H(j\Omega)=\frac{1/(j\Omega C)}{j\Omega L+1/(j\Omega C)}=\frac{1}{1-\Omega^2LC}$, with $H(0)=1$.
>
> | component | $\Omega$ (rad/s) | input (V) | $\lvert H\rvert$ | output (V) |
> |---|---|---|---|---|
> | DC | 0 | 206.90 | 1 | 206.90 |
> | 100 Hz | 628.3 | 137.93 | 0.005092 | 0.702 |
> | 200 Hz | 1256.6 | 27.59 | 0.001268 | 0.035 |
> | 300 Hz | 1885.0 | 11.82 | 0.000563 | 0.0067 |
>
> At 100 Hz, $\Omega^2LC=628.3^2\times0.5\times10^{-3}=197.4$. The output is $206.9$ V DC with $0.70$ V of 100 Hz ripple; the RMS ripple by Parseval is $0.497$ V, a ripple factor of $0.24\%$.
>
> **Interpretation.** The second-order filter attenuates by $\approx1/(\Omega^2LC)$, i.e. 40 dB/decade, while the rectified-sine harmonics themselves fall as $1/k^2$; the combination makes the 200 Hz and higher ripple negligible. Only the lowest ripple harmonic matters for the design, which is why rectifier filter formulas are written for the 100 Hz (or 120 Hz) component alone.

> [!example] Worked Example 15.34 — Choosing a filter from a harmonic specification
> A 1 kHz, $\pm1$ V square wave is to be converted to a sine wave by a filter that leaves the fundamental essentially unchanged. What attenuation is needed at 3 kHz and above for the output THD to be below 1%? Is a second-order Butterworth low-pass with $f_c=1$ kHz sufficient?
>
> **Solution.** The output harmonics are $\frac{4}{\pi n}\lvert H(jn\omega_0)\rvert$. If the filter's attenuation grows at least as fast as $n^2$ beyond the 3rd harmonic, the 3rd dominates the THD and the requirement is approximately $\frac{1}{3}\frac{\lvert H(j3\omega_0)\rvert}{\lvert H(j\omega_0)\rvert}\le0.01$, i.e. a relative attenuation of at least $33.3$, or $30.5$ dB, at 3 kHz.
>
> A second-order Butterworth filter has $\lvert H(j\omega)\rvert=1/\sqrt{1+(\omega/\omega_c)^4}$. With $f_c=1$ kHz: $\lvert H(j\omega_0)\rvert=0.7071$ and $\lvert H(j3\omega_0)\rvert=1/\sqrt{82}=0.1104$, a relative attenuation of only $6.4$ (16 dB). Output harmonic ratios: 3rd $\frac13\cdot\frac{0.1104}{0.7071}=5.2\%$, 5th $\frac15\cdot\frac{0.0400}{0.7071}=1.1\%$, 7th $0.4\%$. THD $\approx\sqrt{5.2^2+1.1^2+0.4^2+\cdots}\%\approx5.3\%$. Not sufficient.
>
> Moving the corner down to $f_c=400$ Hz: $\lvert H(j\omega_0)\rvert=1/\sqrt{1+2.5^4}=0.1580$ and $\lvert H(j3\omega_0)\rvert=1/\sqrt{1+7.5^4}=0.01778$; 3rd-harmonic ratio $\frac13\cdot\frac{0.01778}{0.1580}=3.75\%$ — better, but the fundamental is now attenuated by 16 dB, and the filter still falls short. In fact, for this Butterworth filter the relative attenuation $\lvert H(j\omega_0)\rvert/\lvert H(j3\omega_0)\rvert=\sqrt{(1+81r^4)/(1+r^4)}$, $r=f_0/f_c$, can never exceed $9$ (19 dB, the 40 dB/decade asymptote over a factor of 3 in frequency), whatever the corner. Since the 3rd harmonic needs $30.5$ dB more attenuation than the fundamental, a much steeper roll-off is required: a filter of order 4 or more (e.g. a 4th-order Butterworth with $f_c\approx1$ kHz gives a relative attenuation of $\sqrt{(1+3^8)/2}=57.3$, 35 dB) or a band-pass tuned to 1 kHz.
>
> **Lesson.** Harmonic specifications translate directly into frequency-response requirements at a few discrete frequencies — a far easier design problem than specifying a time-domain waveform.

### Other uses of the same idea

- **Partial differential equations.** Separation of variables for the heat, wave and Laplace equations ([Ch. 13](../part-4-differential-equations/13-partial-differential-equations.md)) expands the initial or boundary data in a half-range sine or cosine series (§15.4); each spatial harmonic then evolves independently, decaying like $e^{-\kappa(n\pi/L)^2t}$ in the heat equation. This is the same "multiply harmonic by harmonic" principle with the role of $H(jn\omega_0)$ played by the time evolution of each mode.
- **Periodic signals and the Fourier transform.** A periodic signal has the generalised Fourier transform $X(\omega)=2\pi\sum_nc_n\delta(\omega-n\omega_0)$ — a train of impulses whose areas are the line-spectrum coefficients ([Ch. 16](16-fourier-transforms.md), §16.6). The ideal impulse-train sampler of §16.8 is analysed exactly this way: its Fourier series has $c_n=1/T_s$ for all $n$.
- **Nonlinear devices with periodic inputs.** A memoryless nonlinearity maps a periodic input to a periodic output with the same fundamental, whose coefficients follow from the multiplication rule (§15.7) or directly from the Euler integrals of the distorted waveform. Describing-function analysis of nonlinear feedback loops keeps only the fundamental of the output.
- **Switching power converters.** Their waveforms are piecewise linear (§15.7, jump method), and harmonic analysis of the switch-node voltage predicts ripple and EMI; the "average model" of a converter is exactly the $n=0$ term.

> [!warning] Pitfall — the steady-state formula hides the transient and assumes stability
> The harmonic-by-harmonic result describes only the periodic steady state of a *stable* system. It says nothing about the start-up transient (which can be much larger, e.g. inrush currents), and it fails outright for marginally stable systems driven at a resonance, as in the undamped case of Worked Example 15.32. Use the Laplace transform (Ch. 14) when the transient matters.

> [!tip] Engineering note — harmonic power flow studies
> Power-system harmonic studies apply this section's theorem at network scale: each nonlinear load is modelled as a set of harmonic current sources $I_n$, the network as an impedance matrix $\mathbf{Z}(jn\omega_0)$ computed separately for each harmonic, and the harmonic bus voltages as $\mathbf{V}_n=\mathbf{Z}(jn\omega_0)\mathbf{I}_n$. Parallel resonances between line inductance and power-factor-correction capacitors, at which $\lvert Z\rvert$ peaks near the 5th or 7th harmonic, are the network-scale analogue of Worked Example 15.32.

## Summary

| Topic | Key formula or fact |
|---|---|
| Fundamental | $T_0$, $f_0=1/T_0$, $\omega_0=2\pi/T_0$; $n$th harmonic at $n\omega_0$ |
| Inner product, RMS | $\langle x,y\rangle=\frac1{T_0}\int_{T_0}xy^*dt$; $\lVert x\rVert=X_{\text{rms}}$ |
| Orthogonality | $\frac1{T_0}\int_{T_0}e^{jm\omega_0t}e^{-jn\omega_0t}dt=\delta_{mn}$ |
| Trigonometric series | $x\sim\frac{a_0}2+\sum_{n\ge1}(a_n\cos n\omega_0t+b_n\sin n\omega_0t)$ |
| Euler formulas | $a_n=\frac2{T_0}\int_{T_0}x\cos n\omega_0t\,dt$, $b_n=\frac2{T_0}\int_{T_0}x\sin n\omega_0t\,dt$ |
| Amplitude–phase form | $A_n=\sqrt{a_n^2+b_n^2}$, $\phi_n=\operatorname{atan2}(-b_n,a_n)$, $x\sim A_0+\sum A_n\cos(n\omega_0t+\phi_n)$ |
| Complex series | $x\sim\sum_nc_ne^{jn\omega_0t}$, $c_n=\frac1{T_0}\int_{T_0}xe^{-jn\omega_0t}dt$ |
| Conversions | $c_n=\frac{a_n-jb_n}2$, $c_0=\frac{a_0}{2}$; real $x$: $c_{-n}=c_n^*$, $A_n=2\lvert c_n\rvert$, $\phi_n=\angle c_n$ |
| One-period transform | $c_n=\frac1{T_0}X_1(jn\omega_0)$, $X_1(s)=\int_0^{T_0}xe^{-st}dt$ |
| Best approximation | $S_N$ minimises $\lVert x-p\rVert$; error $=\lVert x\rVert^2-\sum_{\lvert n\rvert\le N}\lvert c_n\rvert^2$ |
| Even / odd | even: cosines, $a_n=\frac4{T_0}\int_0^{T_0/2}x\cos$; odd: sines, $b_n=\frac4{T_0}\int_0^{T_0/2}x\sin$ |
| Half-wave symmetry | $x(t+\frac{T_0}2)=-x(t)$ $\iff$ odd harmonics only |
| Quarter-wave symmetry | odd $n$ only, $\frac{8}{T_0}\int_0^{T_0/4}$ |
| Half-range series on $(0,L)$ | sine: $b_n=\frac2L\int_0^Lf\sin\frac{n\pi x}L$; cosine: $a_n=\frac2L\int_0^Lf\cos\frac{n\pi x}L$ |
| Standard series | square $\frac{4A}\pi\sum_{\text{odd}}\frac{\sin n\omega_0t}{n}$; triangle $\frac{8A}{\pi^2}\sum_{\text{odd}}\frac{\cos n\omega_0t}{n^2}$; pulse train $c_n=AD\operatorname{sinc}(nD)$ |
| Rectified sines | half-wave: $\frac{V_m}{\pi}+\frac{V_m}2\sin\omega_0t-\frac{2V_m}\pi\sum_k\frac{\cos2k\omega_0t}{4k^2-1}$; full-wave: $\frac{2V_m}\pi-\frac{4V_m}\pi\sum_k\frac{\cos2k\omega_0t}{4k^2-1}$ |
| Dirichlet theorem | $S_N(t)\to\frac12[x(t^+)+x(t^-)]$ |
| Uniform convergence | $x$ continuous and piecewise smooth, or $\sum\lvert c_n\rvert<\infty$ |
| Mean-square convergence | $\lVert x-S_N\rVert\to0$ for every finite-power $x$ |
| Gibbs phenomenon | overshoot $\to(\frac1\pi\operatorname{Si}(\pi)-\frac12)J=0.0895J$; square wave peak $\frac2\pi\operatorname{Si}(\pi)A=1.179A$ |
| Shift, modulation | $x(t-t_0)\leftrightarrow e^{-jn\omega_0t_0}c_n$; $e^{jm\omega_0t}x\leftrightarrow c_{n-m}$ |
| Differentiation, integration | $x'\leftrightarrow jn\omega_0c_n$; $\int x\leftrightarrow c_n/(jn\omega_0)$ if $c_0=0$ |
| Convolution, product | $\frac1{T_0}\int_{T_0}x(\tau)y(t-\tau)d\tau\leftrightarrow c_nd_n$; $xy\leftrightarrow\sum_mc_md_{n-m}$ |
| Jump formula | $c_n=\frac{1}{jn\omega_0T_0}\sum_iJ_ie^{-jn\omega_0t_i}+\frac{1}{jn\omega_0}c_n(x'_{\text{reg}})$ |
| Decay rate | jump in $x^{(k)}$, lower derivatives continuous $\Rightarrow\lvert c_n\rvert=O(1/n^{k+1})$ |
| DFT aliasing | $\tilde c_n=\frac1N\sum_kx_ke^{-j2\pi nk/N}=\sum_mc_{n+mN}$ |
| Parseval | $\frac1{T_0}\int_{T_0}\lvert x\rvert^2dt=\sum_n\lvert c_n\rvert^2=A_0^2+\sum_{n\ge1}\frac{A_n^2}2$ |
| Average power | $P=V_0I_0+\sum_{n\ge1}V_{n,\text{rms}}I_{n,\text{rms}}\cos\theta_n$ |
| THD | $\mathrm{THD}_F=\sqrt{\sum_{n\ge2}X_n^2}/X_1=\sqrt{X^2/X_1^2-1}$; $\mathrm{THD}_R=\mathrm{THD}_F/\sqrt{1+\mathrm{THD}_F^2}$ |
| Power factor (sinusoidal $v$) | $\mathrm{PF}=\frac{I_1}{I}\cos\theta_1=\frac{\cos\theta_1}{\sqrt{1+\mathrm{THD}_i^2}}$ |
| THD values | square 48.3%, sawtooth 80.3%, six-step 31.1%, triangle 12.1% |
| LTI steady state | $d_n=c_nH(jn\omega_0)$; $P_y=\sum_n\lvert c_n\rvert^2\lvert H(jn\omega_0)\rvert^2$ |

## Further reading

- E. Kreyszig, *Advanced Engineering Mathematics*, 10th ed., Wiley — Ch. 11: Fourier series, half-range expansions and the approximation viewpoint, with many worked examples.
- G. P. Tolstov, *Fourier Series*, Dover — a careful, readable account of convergence, the Dirichlet kernel, Gibbs phenomenon and term-by-term operations.
- A. V. Oppenheim, A. S. Willsky & S. H. Nawab, *Signals and Systems*, 2nd ed., Prentice Hall — Ch. 3: complex Fourier series, properties and the response of LTI systems to periodic inputs.
- J. Arrillaga & N. R. Watson, *Power System Harmonics*, 2nd ed., Wiley — harmonic sources, THD and power factor, and harmonic power-flow analysis.
