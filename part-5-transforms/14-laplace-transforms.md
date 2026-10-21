# Chapter 14 — Laplace Transforms

[← Previous](../part-4-differential-equations/13-partial-differential-equations.md) · [Table of Contents](../00-front/table-of-contents.md) · [Next →](15-fourier-series.md)

> **Prerequisites:** Ch. 2 (complex numbers), Ch. 5 (integral calculus, improper integrals), Ch. 11 and Ch. 12 (ordinary and linear differential equations)
> **Learning outcomes:** After this chapter you can
> - define the unilateral Laplace transform with lower limit $0^-$, decide whether it exists, and find its region of convergence;
> - use the transform properties (shifts, scaling, differentiation, integration, multiplication by $t$, periodicity) to build transforms without integrating;
> - invert rational transforms by partial fractions, including repeated and complex-conjugate poles and delayed terms $e^{-sT}$;
> - apply the convolution theorem and the initial- and final-value theorems, and know exactly when the latter is valid;
> - solve linear constant-coefficient differential, integro-differential and simultaneous equations as initial-value problems in one pass;
> - model circuits directly in the $s$-domain, form transfer functions, read poles and zeros, and judge stability from pole locations.

A linear circuit with capacitors and inductors obeys a linear differential equation. In [Ch. 12](../part-4-differential-equations/12-linear-differential-equations.md) such equations were solved by finding a complementary function, guessing a particular integral, and then fitting constants to the initial conditions — three separate steps, each with its own traps. The Laplace transform turns the whole initial-value problem into **algebra**: differentiation becomes multiplication by $s$, the initial conditions enter automatically, and the solution is recovered from a table after a partial-fraction expansion. Beyond solving equations, the transform gives the language in which engineers *think* about linear systems: impedances $sL$ and $1/(sC)$, transfer functions, poles and zeros, and stability. The Fourier transform of [Ch. 16](16-fourier-transforms.md) is the special case $s=j\omega$, and the Z-transform of [Ch. 17](17-z-transforms.md) is its discrete-time sibling. The complex-variable theory behind the inversion integral is developed in [Ch. 18](18-complex-analysis.md).

## 14.1 Definition and existence

### Motivation: from differentiation to multiplication

The exponential $e^{st}$ is the one function that differentiation leaves essentially unchanged: $\frac{d}{dt}e^{st}=s\,e^{st}$. If a signal could be written as a superposition of exponentials $e^{st}$ with weights $X(s)$, then differentiating the signal would simply multiply each weight by $s$. The Laplace transform is the recipe for computing those weights. Unlike the Fourier transform, which uses only the undamped exponentials $e^{j\omega t}$, it allows a real part $\sigma$ in $s=\sigma+j\omega$. The factor $e^{-\sigma t}$ tames signals that grow — ramps, unstable responses, $e^{3t}$ — so that many signals with no ordinary Fourier transform still have a Laplace transform.

### The unilateral transform

> [!definition] Unilateral Laplace transform
> Let $x(t)$ be defined for $t\ge 0$ (its values for $t<0$ are ignored). Its **Laplace transform** is
>
> $$
> X(s)=\mathcal{L}\{x(t)\}=\int_{0^-}^{\infty} x(t)\,e^{-st}\,dt=\lim_{\varepsilon\to0^+,\ T\to\infty}\int_{-\varepsilon}^{T}x(t)\,e^{-st}\,dt ,
> $$
>
> where $s=\sigma+j\omega$ is a complex variable. The set of $s$ for which the integral converges is the **region of convergence (ROC)**. We write $x(t)\leftrightarrow X(s)$ and call $t$ the *time domain* and $s$ the *$s$-domain* or *complex-frequency domain*.

Three points of convention deserve attention from the start.

1. **Lower limit $0^-$.** The integral starts just *before* $t=0$. An impulse $\delta(t)$ or a jump located at the origin is therefore fully inside the range of integration, so $\mathcal{L}\{\delta(t)\}=1$, and the initial conditions that appear in the derivative rule are the values $x(0^-)$ *before* any switching event. These are exactly the values an engineer knows: the capacitor voltage and inductor current just before the switch closes. With a $0^+$ lower limit, one needs instead the values just after switching, which are usually what one is trying to find. This book uses $0^-$ throughout.
2. **Causality.** Because only $t\ge0$ matters, $x(t)$ and $x(t)u(t)$ have the same transform. It is good practice to write time functions with an explicit $u(t)$, for example $e^{-2t}u(t)$, whenever the distinction could matter (shifts, convolution). Tables often omit the $u(t)$; it is always understood.
3. **Units.** Since $dt$ carries seconds, $X(s)$ has the units of $x$ times seconds, and $s$ has units $\mathrm{s^{-1}}$ (it is a complex frequency). A voltage step $V_0u(t)$ in volts has transform $V_0/s$ in V·s. Checking that every term of an $s$-domain expression has the same units is as useful here as in the time domain.

The **bilateral** transform $\int_{-\infty}^{\infty}x(t)e^{-st}dt$ is used for non-causal signals in signal theory; it requires the ROC to be stated to make the inverse unique (two different signals, $e^{-at}u(t)$ and $-e^{-at}u(-t)$, both have bilateral transform $1/(s+a)$, with ROCs $\sigma>-a$ and $\sigma<-a$). For initial-value problems and circuits the unilateral transform is the right tool and all statements below refer to it.

### First transforms from the definition

> [!example] Worked Example 14.1 — the step and the exponential
> Find $\mathcal{L}\{u(t)\}$ and $\mathcal{L}\{e^{at}\}$ for real or complex constant $a$, with their ROCs.
>
> **Step.** For $t>0$, $u(t)=1$:
>
> $$
> \int_{0^-}^{\infty}e^{-st}dt=\lim_{T\to\infty}\left[\frac{-e^{-st}}{s}\right]_0^{T}=\frac{1}{s}-\lim_{T\to\infty}\frac{e^{-\sigma T}e^{-j\omega T}}{s}.
> $$
>
> The factor $e^{-j\omega T}$ has modulus 1, so the limit is zero exactly when $\sigma>0$. Hence $\mathcal{L}\{u(t)\}=1/s$ with ROC $\operatorname{Re}\{s\}>0$.
>
> **Exponential.** $\int_{0}^{\infty}e^{at}e^{-st}dt=\int_0^\infty e^{-(s-a)t}dt$ is the same integral with $s$ replaced by $s-a$:
>
> $$
> \mathcal{L}\{e^{at}\}=\frac{1}{s-a},\qquad \operatorname{Re}\{s\}>\operatorname{Re}\{a\}.
> $$
>
> **Check.** With $a=-2$: $\int_0^\infty e^{-2t}e^{-st}dt$ at $s=1$ is $\int_0^\infty e^{-3t}dt=1/3=1/(1+2)$. ✓

> [!example] Worked Example 14.2 — the ramp, by parts
> Find $\mathcal{L}\{t\}$.
>
> Integrate by parts with $f=t$, $dg=e^{-st}dt$, valid for $\sigma>0$:
>
> $$
> \int_0^\infty t\,e^{-st}dt=\left[-\frac{t\,e^{-st}}{s}\right]_0^\infty+\frac{1}{s}\int_0^\infty e^{-st}dt=0+\frac{1}{s}\cdot\frac{1}{s}=\frac{1}{s^2}.
> $$
>
> The boundary term vanishes at infinity because $te^{-\sigma t}\to0$ for $\sigma>0$ (exponential decay beats polynomial growth). Repeating the argument $n$ times gives $\mathcal{L}\{t^n\}=n!/s^{n+1}$ (Section 14.3).
>
> **Check.** At $s=2$: $\int_0^\infty te^{-2t}dt=1/4$, and $1/s^2=1/4$. ✓

> [!example] Worked Example 14.3 — the impulse and a sifting check
> Find $\mathcal{L}\{\delta(t)\}$ and $\mathcal{L}\{\delta(t-T)\}$ for $T>0$.
>
> By the sifting property $\int g(t)\delta(t-t_0)dt=g(t_0)$ for any $t_0$ inside the range of integration,
>
> $$
> \mathcal{L}\{\delta(t)\}=\int_{0^-}^{\infty}\delta(t)e^{-st}dt=e^{0}=1,\qquad \mathcal{L}\{\delta(t-T)\}=e^{-sT},
> $$
>
> both for every $s$ (ROC: the whole plane). With a $0^+$ lower limit, the impulse at the origin would be *excluded* and the transform would be 0 — a clear sign that $0^-$ is the right choice for systems driven by impulses.

### Existence: exponential order

Not every function has a Laplace transform. The integral can fail at the lower end (a non-integrable singularity) or at infinity (growth faster than any exponential). A simple sufficient condition rules out both.

> [!definition] Piecewise continuity and exponential order
> A function $x(t)$ is **piecewise continuous** on $[0,\infty)$ if on every finite interval it has at most finitely many discontinuities, each a finite jump.
> It is of **exponential order $\alpha$** if there are constants $M>0$, $T_0\ge0$ such that
>
> $$
> \lvert x(t)\rvert\le M e^{\alpha t}\quad\text{for all } t\ge T_0 .
> $$

Polynomials, $e^{at}$, $\sin\omega t$, $\cosh t$ and all their products are of exponential order. The function $e^{t^2}$ is not: for any $\alpha$, $e^{t^2}/e^{\alpha t}=e^{t(t-\alpha)}\to\infty$.

> [!theorem] Existence theorem
> If $x(t)$ is piecewise continuous on $[0,\infty)$ and of exponential order $\alpha$, then $\mathcal{L}\{x\}$ exists for all $\operatorname{Re}\{s\}>\alpha$, the integral converging absolutely. Moreover $X(s)\to0$ as $\operatorname{Re}\{s\}\to\infty$.

*Proof sketch.* Split $\int_0^\infty=\int_0^{T_0}+\int_{T_0}^\infty$. The first part is the integral of a bounded piecewise continuous function over a finite interval, hence finite. For the second, with $\sigma=\operatorname{Re}\{s\}>\alpha$,

$$
\int_{T_0}^{\infty}\lvert x(t)e^{-st}\rvert\,dt\le M\int_{T_0}^{\infty}e^{-(\sigma-\alpha)t}dt=\frac{M e^{-(\sigma-\alpha)T_0}}{\sigma-\alpha}<\infty .
$$

Absolute convergence implies convergence. For the decay, bound $\lvert x\rvert\le M'e^{\alpha t}$ on all of $[0,\infty)$ (possible after enlarging $M$, since $x$ is bounded on $[0,T_0]$); then $\lvert X(s)\rvert\le M'/(\sigma-\alpha)\to0$. $\square$

The infimum of all $\alpha$ that work is the **abscissa of convergence** $\sigma_c$; the ROC is the open half-plane $\operatorname{Re}\{s\}>\sigma_c$ (a unilateral transform always has a right half-plane as ROC). For a rational $X(s)$ the ROC begins at the right-most pole: $\sigma_c=\max\operatorname{Re}\{p_k\}$.

The conditions are **sufficient, not necessary**. The function $t^{-1/2}$ is unbounded at $t=0$, so it is not piecewise continuous there, yet $\int_0^\infty t^{-1/2}e^{-st}dt=\sqrt{\pi/s}$ converges because the singularity is integrable. Impulses are not functions at all but are handled by the sifting property. In practice, every signal met in circuit analysis has a transform.

> [!example] Worked Example 14.4 — deciding existence and the ROC
> For each signal, decide whether the unilateral transform exists and, if so, give the ROC.
>
> **(a)** $x(t)=3e^{-t}+2e^{4t}$. Both terms are of exponential order; the larger exponent is 4, so $\lvert x\rvert\le5e^{4t}$ and the transform $\frac{3}{s+1}+\frac{2}{s-4}$ exists for $\operatorname{Re}\{s\}>4$. The ROC of a sum is (at least) the intersection of the individual ROCs — here governed by the fastest-growing term, i.e. the right-most pole.
>
> **(b)** $x(t)=e^{t^2}$. Not of exponential order. For any real $s$, the integrand $e^{t^2-st}$ grows without bound once $t>s$, so the integral diverges for every $s$: **no transform**.
>
> **(c)** $x(t)=1/t$. Of exponential order for large $t$, but near the origin $\int_0^{1}t^{-1}e^{-st}dt\ge e^{-\sigma}\int_0^1t^{-1}dt=\infty$ for real $s=\sigma$. The singularity is not integrable: **no transform**. Contrast $t^{-1/2}$, whose singularity is integrable and whose transform $\sqrt{\pi/s}$ exists.
>
> **(d)** $x(t)=\sin(e^{t^2})$. This oscillates ever faster, but $\lvert x\rvert\le1=e^{0\cdot t}$, so it is of exponential order 0 and the transform exists for $\operatorname{Re}\{s\}>0$. Rapid oscillation is harmless; only growth and non-integrable singularities matter.
>
> **(e)** $x(t)=\cosh(2t)\,e^{-3t}$. Since $\cosh2t\,e^{-3t}=\tfrac12(e^{-t}+e^{-5t})$, the transform $\frac{1}{2}\left(\frac{1}{s+1}+\frac{1}{s+5}\right)$ exists for $\operatorname{Re}\{s\}>-1$. Here the ROC includes the $j\omega$-axis, so this signal also has a Fourier transform.

> [!theorem] Uniqueness (Lerch)
> If two piecewise continuous functions of exponential order have the same Laplace transform on some half-plane $\operatorname{Re}\{s\}>\sigma_0$, they are equal at every $t>0$ where both are continuous.

Uniqueness is what makes table look-up legitimate: if we recognise $X(s)$ as the transform of $x(t)$, then $x(t)$ is *the* answer (up to values at isolated points, which no physical measurement can detect).

The behaviour of the integrand explains the ROC graphically. For $x(t)=e^{0.5t}$, the weighted integrand $x(t)e^{-\sigma t}$ decays only when $\sigma>0.5$:

```desmos-graph
left=0; right=8; bottom=-0.5; top=4
---
y=e^{0.5x}|grey|dashed
y=e^{0.5x}e^{-0.3x}|red
y=e^{0.5x}e^{-0.5x}|orange
y=e^{0.5x}e^{-1x}|blue
```
*The signal $e^{0.5t}$ (grey) multiplied by $e^{-\sigma t}$ for $\sigma=0.3$ (red, still growing), $\sigma=0.5$ (orange, constant — integral diverges) and $\sigma=1$ (blue, decaying — integral converges). The ROC is $\sigma>0.5$.*

> [!warning] Common Pitfall — evaluating outside the ROC
> The formula $1/(s-a)$ is defined for all $s\ne a$, but the *integral* converges only for $\operatorname{Re}\{s\}>\operatorname{Re}\{a\}$. For the unilateral transform this rarely causes errors, because the ROC is always a right half-plane determined by the right-most singularity. It matters when passing to the Fourier transform: $X(j\omega)$ is the Fourier transform only if the $j\omega$-axis lies in the ROC. For $e^{2t}u(t)$, $X(s)=1/(s-2)$, but $1/(j\omega-2)$ is **not** its Fourier transform — it has none ([Ch. 16](16-fourier-transforms.md)).

> [!tip] Engineering Note — why engineers accept a complex frequency
> A pole at $s=\sigma+j\omega$ corresponds to a time behaviour $e^{\sigma t}e^{j\omega t}$: $\omega$ is the oscillation frequency in rad/s and $\sigma$ is the growth rate (decay if negative) in nepers per second. Every natural response of a lumped linear circuit is a combination of such terms, so the $s$-plane is a map of every possible "mode" of a circuit, with stability determined by which half of the map the modes lie in (Section 14.9).

## 14.2 Properties

Computing every transform from the integral would be tedious. A handful of properties, each a one-line consequence of the definition, lets us build almost any transform needed in practice from a few basic pairs, and — more importantly — tells us how operations on signals (delay, differentiation, damping) appear in the $s$-domain. Throughout, $x(t)\leftrightarrow X(s)$ with ROC $\operatorname{Re}\{s\}>\sigma_x$, and similarly for $y$.

### Linearity

> [!theorem] Linearity
> For constants $a,b$:
>
> $$
> \mathcal{L}\{a\,x(t)+b\,y(t)\}=a\,X(s)+b\,Y(s),
> $$
>
> valid at least on the intersection of the two ROCs.

This is inherited directly from linearity of the integral. Combined with Euler's formula it already yields the trigonometric transforms: $\cos\omega t=\tfrac12(e^{j\omega t}+e^{-j\omega t})$ gives

$$
\mathcal{L}\{\cos\omega t\}=\frac12\left(\frac{1}{s-j\omega}+\frac{1}{s+j\omega}\right)=\frac{s}{s^2+\omega^2},\qquad
\mathcal{L}\{\sin\omega t\}=\frac{1}{2j}\left(\frac{1}{s-j\omega}-\frac{1}{s+j\omega}\right)=\frac{\omega}{s^2+\omega^2},
$$

both with ROC $\operatorname{Re}\{s\}>0$.

### Frequency shift (multiplication by an exponential)

> [!theorem] First shifting theorem ($s$-shift)
>
> $$
> \mathcal{L}\{e^{at}x(t)\}=X(s-a),\qquad \operatorname{Re}\{s\}>\sigma_x+\operatorname{Re}\{a\}.
> $$

*Proof.* $\int_{0^-}^\infty e^{at}x(t)e^{-st}dt=\int_{0^-}^\infty x(t)e^{-(s-a)t}dt=X(s-a)$. $\square$

Damping a signal by $e^{-\alpha t}$ slides its transform — and all its poles — a distance $\alpha$ to the left in the $s$-plane. Thus

$$
e^{-\alpha t}\cos\omega t\leftrightarrow\frac{s+\alpha}{(s+\alpha)^2+\omega^2},\qquad e^{-\alpha t}\sin\omega t\leftrightarrow\frac{\omega}{(s+\alpha)^2+\omega^2},\qquad t^ne^{-\alpha t}\leftrightarrow\frac{n!}{(s+\alpha)^{n+1}} .
$$

### Time shift (delay)

> [!theorem] Second shifting theorem ($t$-shift)
> For $T\ge0$,
>
> $$
> \mathcal{L}\{x(t-T)\,u(t-T)\}=e^{-sT}X(s).
> $$

*Proof.* The integrand vanishes for $t<T$; substitute $\tau=t-T$:

$$
\int_{T}^{\infty}x(t-T)e^{-st}dt=\int_{0}^{\infty}x(\tau)e^{-s(\tau+T)}d\tau=e^{-sT}X(s).\qquad\square
$$

The factor $e^{-sT}$ is the $s$-domain signature of a pure delay. The theorem applies to the *whole function* shifted, $x(t-T)u(t-T)$; it does not apply to $x(t)u(t-T)$ until that is rewritten. The safe procedure for $\mathcal{L}\{g(t)u(t-T)\}$ is

$$
\mathcal{L}\{g(t)\,u(t-T)\}=e^{-sT}\,\mathcal{L}\{g(t+T)\},
$$

which follows from the same substitution.

> [!example] Worked Example 14.5 — piecewise signals with steps
> A voltage rises linearly from 0 to 1 V during the first second and then stays at 1 V:
> $x(t)=t$ for $0\le t<1$, $x(t)=1$ for $t\ge1$. Find $X(s)$.
>
> **Write with steps.** A ramp that starts at $t=0$, minus a ramp that starts at $t=1$, levels the signal off:
>
> $$
> x(t)=t\,u(t)-(t-1)\,u(t-1).
> $$
>
> (Check: for $t\ge1$, $t-(t-1)=1$. ✓)
>
> **Transform.** Each term is a ramp, one of them delayed:
>
> $$
> X(s)=\frac{1}{s^2}-\frac{e^{-s}}{s^2}=\frac{1-e^{-s}}{s^2}.
> $$
>
> **Check by integration.** $\int_0^1te^{-st}dt+\int_1^\infty e^{-st}dt$ at $s=1$: the first integral is $1-2e^{-1}=0.26424$, the second $e^{-1}=0.36788$; the sum is $0.63212$. The formula gives $1-e^{-1}=0.63212$. ✓

> [!example] Worked Example 14.6 — a shifted sinusoid
> Find $\mathcal{L}\{\sin(t)\,u(t-\pi/2)\}$.
>
> The sine is not written as a function of $t-\pi/2$. Use $g(t)u(t-T)\leftrightarrow e^{-sT}\mathcal{L}\{g(t+T)\}$ with $g=\sin$, $T=\pi/2$: $g(t+\pi/2)=\sin(t+\pi/2)=\cos t$, so
>
> $$
> \mathcal{L}\{\sin t\;u(t-\pi/2)\}=e^{-\pi s/2}\,\frac{s}{s^2+1}.
> $$
>
> The tempting but wrong answer $e^{-\pi s/2}/(s^2+1)$ would be the transform of $\sin(t-\pi/2)u(t-\pi/2)=-\cos t\,u(t-\pi/2)$, a signal of opposite sign.

### Time scaling

> [!theorem] Scaling
> For $a>0$,
>
> $$
> \mathcal{L}\{x(at)\}=\frac{1}{a}\,X\!\left(\frac{s}{a}\right).
> $$

*Proof.* Substitute $\tau=at$: $\int_0^\infty x(at)e^{-st}dt=\frac1a\int_0^\infty x(\tau)e^{-(s/a)\tau}d\tau$. $\square$

Compressing in time ($a>1$) expands in $s$: poles move outward by the factor $a$. This is the basis of *frequency scaling* of filter designs: a prototype with cut-off 1 rad/s becomes one with cut-off $\omega_c$ by replacing $s$ with $s/\omega_c$.

### Differentiation in time

This is the property that makes the transform useful for differential equations.

> [!theorem] Derivative rule
> If $x$ is continuous on $(0^-,\infty)$ with piecewise continuous derivative, both of exponential order, then
>
> $$
> \mathcal{L}\{\dot x(t)\}=sX(s)-x(0^-),
> $$
>
> and, applying the rule repeatedly,
>
> $$
> \mathcal{L}\{x^{(n)}(t)\}=s^nX(s)-s^{n-1}x(0^-)-s^{n-2}\dot x(0^-)-\cdots-x^{(n-1)}(0^-).
> $$

*Derivation.* Integrate by parts:

$$
\int_{0^-}^{\infty}\dot x(t)e^{-st}dt=\Big[x(t)e^{-st}\Big]_{0^-}^{\infty}+s\int_{0^-}^{\infty}x(t)e^{-st}dt=-x(0^-)+sX(s),
$$

the upper boundary term vanishing in the ROC because $x$ is of exponential order. For the second derivative apply the rule to $\dot x$: $\mathcal{L}\{\ddot x\}=s\mathcal{L}\{\dot x\}-\dot x(0^-)=s^2X-sx(0^-)-\dot x(0^-)$. Induction gives the general formula. $\square$

**Jumps at the origin.** With the $0^-$ convention, the rule remains true even when $x$ jumps at $t=0$, provided $\dot x$ is understood in the generalized sense (a jump of size $J$ contributes $J\,\delta(t)$ to $\dot x$). For instance $x(t)=e^{-t}u(t)$ has $x(0^-)=0$ and $\dot x=\delta(t)-e^{-t}u(t)$, whose transform is $1-\frac{1}{s+1}=\frac{s}{s+1}=sX(s)-0$. ✓ A version of the rule using $x(0^+)=1$ reads $sX-x(0^+)=-\frac{1}{s+1}$; it describes only the derivative for $t>0$ and silently drops the impulse. Both are self-consistent, but mixing them is a classic source of error; the $0^-$ form never needs to know what happens *at* the switching instant.

A jump at a later time $t_1>0$ of size $J$ is handled the same way: $\mathcal{L}\{\dot x\}=sX-x(0^-)$ still holds if $\dot x$ includes the impulse $J\delta(t-t_1)$, while the ordinary (pointwise) derivative satisfies $\mathcal{L}\{\dot x_{\text{ord}}\}=sX-x(0^-)-Je^{-st_1}$.

### Integration in time

> [!theorem] Integral rule
>
> $$
> \mathcal{L}\left\{\int_{0^-}^{t}x(\tau)\,d\tau\right\}=\frac{X(s)}{s}.
> $$

*Proof.* Let $g(t)=\int_{0^-}^tx\,d\tau$. Then $\dot g=x$ and $g(0^-)=0$, so the derivative rule gives $X=sG-0$. $\square$

Integration divides by $s$; differentiation multiplies by $s$. The integrator $1/s$ and the differentiator $s$ are the two elementary dynamic operators, and $u(t)\leftrightarrow1/s$ is the integral of $\delta(t)\leftrightarrow1$. If the integral starts from an earlier time with a stored value, as for a capacitor voltage $v_C(t)=v_C(0^-)+\frac1C\int_{0^-}^t i\,d\tau$, the transform is $V_C=v_C(0^-)/s+I/(sC)$ — the form used for circuit models in Section 14.9.

### Multiplication by $t$ and division by $t$

> [!theorem] Differentiation in $s$
>
> $$
> \mathcal{L}\{t\,x(t)\}=-\frac{dX}{ds},\qquad \mathcal{L}\{t^n x(t)\}=(-1)^n\frac{d^nX}{ds^n}.
> $$

*Proof.* In the ROC the integral converges uniformly on closed sub-half-planes, so we may differentiate under the integral sign: $\frac{d}{ds}\int x(t)e^{-st}dt=\int(-t)x(t)e^{-st}dt$. $\square$

A by-product: inside its ROC, $X(s)$ is differentiable as a function of the complex variable $s$, i.e. **analytic** ([Ch. 18](18-complex-analysis.md)). Singularities (poles) can only occur on or to the left of the boundary of the ROC.

> [!theorem] Integration in $s$
> If $x(t)/t$ has a finite limit as $t\to0^+$,
>
> $$
> \mathcal{L}\left\{\frac{x(t)}{t}\right\}=\int_{s}^{\infty}X(p)\,dp .
> $$

*Proof sketch.* $\int_s^\infty X(p)dp=\int_s^\infty\int_0^\infty x(t)e^{-pt}dt\,dp$; exchanging the order, $\int_s^\infty e^{-pt}dp=e^{-st}/t$, which leaves $\int_0^\infty \frac{x(t)}{t}e^{-st}dt$. $\square$

> [!example] Worked Example 14.7 — combining $t$-multiplication and $s$-shift
> Find $\mathcal{L}\{t\,e^{-2t}\sin3t\}$.
>
> **Step 1 (basic pair).** $\sin3t\leftrightarrow\dfrac{3}{s^2+9}$.
>
> **Step 2 (multiply by $t$).**
>
> $$
> t\sin3t\leftrightarrow-\frac{d}{ds}\frac{3}{s^2+9}=\frac{6s}{(s^2+9)^2}.
> $$
>
> **Step 3 (damp by $e^{-2t}$: replace $s$ by $s+2$).**
>
> $$
> t\,e^{-2t}\sin3t\leftrightarrow\frac{6(s+2)}{\big((s+2)^2+9\big)^2}=\frac{6(s+2)}{(s^2+4s+13)^2}.
> $$
>
> **Check.** At $s=0$ the formula gives $12/169=0.071006$. Direct numerical integration of $\int_0^\infty te^{-2t}\sin3t\,dt$ gives $0.071006$. ✓ The order of steps 2 and 3 does not matter, but doing the $s$-shift last keeps the algebra simplest.

> [!example] Worked Example 14.8 — the sine integral by division by $t$
> Find $\mathcal{L}\{\sin t/t\}$ and deduce $\int_0^\infty\frac{\sin t}{t}dt$.
>
> Since $\sin t/t\to1$ as $t\to0$, the rule applies:
>
> $$
> \mathcal{L}\left\{\frac{\sin t}{t}\right\}=\int_s^\infty\frac{dp}{p^2+1}=\frac{\pi}{2}-\arctan s=\arctan\frac{1}{s}.
> $$
>
> Letting $s\to0^+$ (justified by Abel's theorem for Laplace integrals, since $\int_0^\infty\sin t/t\,dt$ converges) gives
>
> $$
> \int_0^\infty\frac{\sin t}{t}dt=\lim_{s\to0^+}\arctan\frac1s=\frac{\pi}{2}\approx1.5708 .
> $$
>
> This integral (the area under the unnormalized sinc) is awkward by elementary methods; here it costs two lines. It reappears in Gibbs' phenomenon in [Ch. 15](15-fourier-series.md).

### Periodic functions

> [!theorem] Transform of a periodic function
> If $x(t+T)=x(t)$ for $t\ge0$, and $x_1(t)$ denotes one period ($x$ on $[0,T)$, zero elsewhere) with transform $X_1(s)$, then
>
> $$
> X(s)=\frac{X_1(s)}{1-e^{-sT}}=\frac{1}{1-e^{-sT}}\int_0^{T}x(t)e^{-st}dt,\qquad\operatorname{Re}\{s\}>0.
> $$

*Proof.* $x(t)=\sum_{k=0}^\infty x_1(t-kT)$, so by the time-shift rule $X=X_1\sum_{k=0}^\infty e^{-ksT}$, a geometric series with ratio $\lvert e^{-sT}\rvert=e^{-\sigma T}<1$. $\square$

> [!example] Worked Example 14.9 — square wave and full-wave rectified sine
> **(a)** A unipolar square wave equals 1 on $[0,T/2)$ and 0 on $[T/2,T)$, repeating with period $T$. One period is $u(t)-u(t-T/2)$, so $X_1=(1-e^{-sT/2})/s$ and
>
> $$
> X(s)=\frac{1-e^{-sT/2}}{s(1-e^{-sT})}=\frac{1}{s\,(1+e^{-sT/2})},
> $$
>
> using $1-e^{-sT}=(1-e^{-sT/2})(1+e^{-sT/2})$.
>
> **(b)** The full-wave rectified sine $\lvert\sin t\rvert$ has period $\pi$. One period:
>
> $$
> \int_0^\pi\sin t\,e^{-st}dt=\frac{1+e^{-\pi s}}{s^2+1}
> \quad\Longrightarrow\quad
> X(s)=\frac{1+e^{-\pi s}}{(s^2+1)(1-e^{-\pi s})}=\frac{1}{s^2+1}\coth\frac{\pi s}{2}.
> $$
>
> **Check.** At $s=1$, $\coth(\pi/2)/2=0.54517$; numerical integration of $\int_0^\infty\lvert\sin t\rvert e^{-t}dt$ gives $0.54517$. ✓ As $s\to0$, $sX(s)\to2/\pi\approx0.6366$, the average value of $\lvert\sin t\rvert$ — a preview of the final-value idea in Section 14.7.

> [!warning] Common Pitfall — shifting in the wrong domain
> Multiplication by $e^{at}$ in time shifts in $s$ ($s\to s-a$); multiplication by $e^{-sT}$ in $s$ delays in time. The signs are easy to swap. A quick check: damping ($a<0$) must move poles to the **left**, and a delay must leave $X$ unchanged in magnitude on the $j\omega$-axis ($\lvert e^{-j\omega T}\rvert=1$).

> [!tip] Engineering Note — properties as building blocks
> The properties map one-to-one onto block-diagram components: $1/s$ is an integrator (an op-amp integrator stage), $s$ a differentiator, $e^{-sT}$ a transport delay (a cable, a sampling latency), and $s\to s/\omega_c$ a frequency scaling of a filter prototype. Circuit simulators that work in the $s$-domain (small-signal AC analysis) are built on exactly these rules.

## 14.3 Standard transforms

A working engineer inverts transforms by recognising patterns, so a reliable table is part of the toolkit. This section collects the standard pairs, derives the less obvious ones, and shows how the table is extended with the properties of Section 14.2. All time functions are understood to be multiplied by $u(t)$; $a,\alpha,\omega,b$ are real constants unless stated otherwise, $n$ is a non-negative integer, and $\nu>-1$ is real.

### Powers of $t$ and the Gamma function

For integer powers, repeated use of the $t$-multiplication rule on $1\leftrightarrow1/s$ gives

$$
t^n=t\cdot t^{n-1}\leftrightarrow-\frac{d}{ds}\frac{(n-1)!}{s^n}=\frac{n!}{s^{n+1}} .
$$

For non-integer powers, substitute $\tau=st$ (first for real $s>0$, then extend by analytic continuation):

$$
\int_0^\infty t^\nu e^{-st}dt=\frac{1}{s^{\nu+1}}\int_0^\infty\tau^\nu e^{-\tau}d\tau=\frac{\Gamma(\nu+1)}{s^{\nu+1}},\qquad \nu>-1,
$$

where $\Gamma(z)=\int_0^\infty\tau^{z-1}e^{-\tau}d\tau$ is the **Gamma function**. It satisfies $\Gamma(z+1)=z\Gamma(z)$ (integrate by parts), $\Gamma(n+1)=n!$, and $\Gamma(\tfrac12)=\sqrt\pi$ (from the Gaussian integral of [Ch. 6](../part-2-calculus/06-multivariable-calculus.md)).

> [!example] Worked Example 14.10 — half-integer powers
> Find $\mathcal{L}\{t^{-1/2}\}$, $\mathcal{L}\{t^{1/2}\}$ and $\mathcal{L}\{t^{3/2}\}$.
>
> $$
> t^{-1/2}\leftrightarrow\frac{\Gamma(\frac12)}{s^{1/2}}=\sqrt{\frac{\pi}{s}},\qquad
> t^{1/2}\leftrightarrow\frac{\Gamma(\frac32)}{s^{3/2}}=\frac{\tfrac12\sqrt\pi}{s^{3/2}},\qquad
> t^{3/2}\leftrightarrow\frac{\Gamma(\frac52)}{s^{5/2}}=\frac{\tfrac34\sqrt\pi}{s^{5/2}} .
> $$
>
> **Check.** At $s=1$, $\int_0^\infty t^{3/2}e^{-t}dt=\Gamma(2.5)=1.32934$, and $\tfrac34\sqrt\pi=1.32934$. ✓
>
> Transforms with $\sqrt s$ arise from distributed systems — diffusion in a semiconductor, heat flow, the $RC$ transmission line — where the "impedance" of an infinitely long line is proportional to $1/\sqrt s$, so step responses contain $\sqrt t$ and $1/\sqrt t$ rather than exponentials.

### Exponentials, sinusoids and their products

Everything in this family follows from $e^{at}\leftrightarrow1/(s-a)$ with complex $a$, plus linearity and the $s$-shift. Two further derivations illustrate the technique.

**Hyperbolic functions.** $\cosh at=\tfrac12(e^{at}+e^{-at})\leftrightarrow\tfrac12\left(\frac{1}{s-a}+\frac{1}{s+a}\right)=\frac{s}{s^2-a^2}$, and similarly $\sinh at\leftrightarrow\frac{a}{s^2-a^2}$ (ROC $\operatorname{Re}\{s\}>\lvert a\rvert$). Compare with $\cos$ and $\sin$: only the sign in the denominator changes, mirroring $\cos(jat)=\cosh at$.

**$t\sin$ and $t\cos$ (resonance terms).** By the $t$-multiplication rule,

$$
t\sin\omega t\leftrightarrow-\frac{d}{ds}\frac{\omega}{s^2+\omega^2}=\frac{2\omega s}{(s^2+\omega^2)^2},\qquad
t\cos\omega t\leftrightarrow-\frac{d}{ds}\frac{s}{s^2+\omega^2}=\frac{s^2-\omega^2}{(s^2+\omega^2)^2}.
$$

Combining, $\sin\omega t-\omega t\cos\omega t\leftrightarrow\frac{\omega(s^2+\omega^2)-\omega(s^2-\omega^2)}{(s^2+\omega^2)^2}=\frac{2\omega^3}{(s^2+\omega^2)^2}$. A *repeated* pair of imaginary poles therefore produces an oscillation whose amplitude grows linearly in $t$ — the signature of a lossless system driven at resonance.

> [!example] Worked Example 14.11 — a general damped sinusoid
> Find the transform of $x(t)=Ae^{-\alpha t}\cos(\omega t+\phi)$.
>
> Expand the cosine: $\cos(\omega t+\phi)=\cos\phi\cos\omega t-\sin\phi\sin\omega t$. Then, by the $s$-shift,
>
> $$
> X(s)=A\,\frac{(s+\alpha)\cos\phi-\omega\sin\phi}{(s+\alpha)^2+\omega^2}.
> $$
>
> **Numerical case.** $A=2$, $\alpha=1$, $\omega=2$, $\phi=\pi/3$: $\cos\phi=0.5$, $\sin\phi=0.8660$, so
>
> $$
> X(s)=\frac{(s+1)-3.4641}{s^2+2s+5}=\frac{s-2.4641}{s^2+2s+5}.
> $$
>
> **Check.** $x(0)=2\cos(\pi/3)=1$; the initial-value theorem (Section 14.7) gives $\lim_{s\to\infty}sX=1$. ✓
> This pair, read backwards, is how complex-conjugate poles are converted to a single real cosine with amplitude and phase (Section 14.5).

### Pulses, staircases and other non-smooth signals

Rectangular pulses, staircases and sawtooth waves are written with steps and handled by the time-shift rule. A rectangular pulse of height 1 from $t=a$ to $t=b$ is $u(t-a)-u(t-b)\leftrightarrow\frac{e^{-as}-e^{-bs}}{s}$. A **staircase** $\sum_{k=0}^\infty u(t-kT)$ (value $k+1$ on $[kT,(k+1)T)$) has transform $\frac{1}{s}\sum_k e^{-ksT}=\frac{1}{s(1-e^{-sT})}$. A sawtooth of slope $1/T$ resetting every $T$ is a ramp minus a staircase: $\frac{t}{T}-\sum_{k\ge1}u(t-kT)\leftrightarrow\frac{1}{Ts^2}-\frac{e^{-sT}}{s(1-e^{-sT})}$.

> [!example] Worked Example 14.12 — a finite pulse train
> A gate signal consists of three unit pulses of width 1 ms starting at $t=0$, 3 ms and 6 ms. Find its transform (work in milliseconds, so $s$ is in $\mathrm{ms^{-1}}$).
>
> Each pulse is $u(t-t_k)-u(t-t_k-1)\leftrightarrow e^{-st_k}(1-e^{-s})/s$. Summing over $t_k=0,3,6$:
>
> $$
> X(s)=\frac{1-e^{-s}}{s}\left(1+e^{-3s}+e^{-6s}\right)=\frac{(1-e^{-s})(1-e^{-9s})}{s(1-e^{-3s})},
> $$
>
> the last form using the finite geometric sum $1+r+r^2=(1-r^3)/(1-r)$ with $r=e^{-3s}$.
>
> **Check.** As $s\to0$, $X(s)\to$ total area: $\frac{(s)(9s)}{s\cdot3s}=3$ ms. Three pulses of 1 ms each. ✓

### The standard transform table

> [!theorem] Table of unilateral Laplace transforms
> (All $x(t)$ are for $t\ge0$, i.e. implicitly multiplied by $u(t)$.)

| # | $x(t)$ | $X(s)$ | ROC |
|---|---|---|---|
| 1 | $\delta(t)$ | $1$ | all $s$ |
| 2 | $\delta(t-T)$, $T\ge0$ | $e^{-sT}$ | all $s$ |
| 3 | $\delta^{(n)}(t)$ | $s^n$ | all $s$ |
| 4 | $u(t)$ | $\dfrac1s$ | $\sigma>0$ |
| 5 | $u(t-T)$ | $\dfrac{e^{-sT}}{s}$ | $\sigma>0$ |
| 6 | $u(t-a)-u(t-b)$, $0\le a<b$ | $\dfrac{e^{-as}-e^{-bs}}{s}$ | all $s$ |
| 7 | $t$ | $\dfrac{1}{s^2}$ | $\sigma>0$ |
| 8 | $t^n$ | $\dfrac{n!}{s^{n+1}}$ | $\sigma>0$ |
| 9 | $t^{\nu}$, $\nu>-1$ | $\dfrac{\Gamma(\nu+1)}{s^{\nu+1}}$ | $\sigma>0$ |
| 10 | $\dfrac{1}{\sqrt{\pi t}}$ | $\dfrac{1}{\sqrt s}$ | $\sigma>0$ |
| 11 | $2\sqrt{\dfrac{t}{\pi}}$ | $\dfrac{1}{s^{3/2}}$ | $\sigma>0$ |
| 12 | $e^{-at}$ | $\dfrac{1}{s+a}$ | $\sigma>-a$ |
| 13 | $t\,e^{-at}$ | $\dfrac{1}{(s+a)^2}$ | $\sigma>-a$ |
| 14 | $\dfrac{t^{n}}{n!}e^{-at}$ | $\dfrac{1}{(s+a)^{n+1}}$ | $\sigma>-a$ |
| 15 | $\dfrac{1}{a}\left(1-e^{-at}\right)$ | $\dfrac{1}{s(s+a)}$ | $\sigma>\max(0,-a)$ |
| 16 | $\dfrac{e^{-at}-e^{-bt}}{b-a}$, $a\ne b$ | $\dfrac{1}{(s+a)(s+b)}$ | $\sigma>\max(-a,-b)$ |
| 17 | $\dfrac{be^{-bt}-ae^{-at}}{b-a}$, $a\ne b$ | $\dfrac{s}{(s+a)(s+b)}$ | $\sigma>\max(-a,-b)$ |
| 18 | $\dfrac{1}{a^2}\left(at-1+e^{-at}\right)$ | $\dfrac{1}{s^2(s+a)}$ | $\sigma>\max(0,-a)$ |
| 19 | $\sin\omega t$ | $\dfrac{\omega}{s^2+\omega^2}$ | $\sigma>0$ |
| 20 | $\cos\omega t$ | $\dfrac{s}{s^2+\omega^2}$ | $\sigma>0$ |
| 21 | $\sin(\omega t+\phi)$ | $\dfrac{s\sin\phi+\omega\cos\phi}{s^2+\omega^2}$ | $\sigma>0$ |
| 22 | $\cos(\omega t+\phi)$ | $\dfrac{s\cos\phi-\omega\sin\phi}{s^2+\omega^2}$ | $\sigma>0$ |
| 23 | $e^{-\alpha t}\sin\omega t$ | $\dfrac{\omega}{(s+\alpha)^2+\omega^2}$ | $\sigma>-\alpha$ |
| 24 | $e^{-\alpha t}\cos\omega t$ | $\dfrac{s+\alpha}{(s+\alpha)^2+\omega^2}$ | $\sigma>-\alpha$ |
| 25 | $1-\cos\omega t$ | $\dfrac{\omega^2}{s(s^2+\omega^2)}$ | $\sigma>0$ |
| 26 | $\omega t-\sin\omega t$ | $\dfrac{\omega^3}{s^2(s^2+\omega^2)}$ | $\sigma>0$ |
| 27 | $t\sin\omega t$ | $\dfrac{2\omega s}{(s^2+\omega^2)^2}$ | $\sigma>0$ |
| 28 | $t\cos\omega t$ | $\dfrac{s^2-\omega^2}{(s^2+\omega^2)^2}$ | $\sigma>0$ |
| 29 | $\sin\omega t-\omega t\cos\omega t$ | $\dfrac{2\omega^3}{(s^2+\omega^2)^2}$ | $\sigma>0$ |
| 30 | $\sin\omega t+\omega t\cos\omega t$ | $\dfrac{2\omega s^2}{(s^2+\omega^2)^2}$ | $\sigma>0$ |
| 31 | $\sinh at$ | $\dfrac{a}{s^2-a^2}$ | $\sigma>\lvert a\rvert$ |
| 32 | $\cosh at$ | $\dfrac{s}{s^2-a^2}$ | $\sigma>\lvert a\rvert$ |
| 33 | $e^{-\alpha t}\sinh at$ | $\dfrac{a}{(s+\alpha)^2-a^2}$ | $\sigma>\lvert a\rvert-\alpha$ |
| 34 | $\sin^2\omega t$ | $\dfrac{2\omega^2}{s(s^2+4\omega^2)}$ | $\sigma>0$ |
| 35 | $\dfrac{\sin\omega t}{t}$ | $\arctan\dfrac{\omega}{s}$ | $\sigma>0$ |
| 36 | $\dfrac{e^{-bt}-e^{-at}}{t}$ | $\ln\dfrac{s+a}{s+b}$ | $\sigma>\max(-a,-b)$ |
| 37 | $\ln t$ | $-\dfrac{\gamma+\ln s}{s}$, $\gamma=0.5772\ldots$ | $\sigma>0$ |
| 38 | $J_0(\omega t)$ (Bessel) | $\dfrac{1}{\sqrt{s^2+\omega^2}}$ | $\sigma>0$ |
| 39 | $\operatorname{erfc}\dfrac{a}{2\sqrt t}$, $a>0$ | $\dfrac{e^{-a\sqrt s}}{s}$ | $\sigma>0$ |
| 40 | $\dfrac{a}{2\sqrt{\pi t^3}}e^{-a^2/4t}$, $a>0$ | $e^{-a\sqrt s}$ | $\sigma>0$ |
| 41 | square wave, 1 on $[0,\frac T2)$, 0 on $[\frac T2,T)$ | $\dfrac{1}{s(1+e^{-sT/2})}$ | $\sigma>0$ |
| 42 | $\lvert\sin\omega t\rvert$ | $\dfrac{\omega}{s^2+\omega^2}\coth\dfrac{\pi s}{2\omega}$ | $\sigma>0$ |
| 43 | staircase $\sum_{k\ge0}u(t-kT)$ | $\dfrac{1}{s(1-e^{-sT})}$ | $\sigma>0$ |
| 44 | $\sum_{k\ge0}\delta(t-kT)$ | $\dfrac{1}{1-e^{-sT}}$ | $\sigma>0$ |

**Operational properties** (companion table — these are the rules of Section 14.2):

| Operation | $x(t)$ | $X(s)$ |
|---|---|---|
| Linearity | $ax+by$ | $aX+bY$ |
| $s$-shift | $e^{at}x(t)$ | $X(s-a)$ |
| Delay | $x(t-T)u(t-T)$ | $e^{-sT}X(s)$ |
| Scaling ($a>0$) | $x(at)$ | $\frac1aX(s/a)$ |
| Derivative | $\dot x$ | $sX-x(0^-)$ |
| Second derivative | $\ddot x$ | $s^2X-sx(0^-)-\dot x(0^-)$ |
| Integral | $\int_{0^-}^tx\,d\tau$ | $X/s$ |
| Multiply by $t$ | $t\,x(t)$ | $-X'(s)$ |
| Divide by $t$ | $x(t)/t$ | $\int_s^\infty X(p)\,dp$ |
| Periodic, period $T$ | $x(t+T)=x(t)$ | $\frac{1}{1-e^{-sT}}\int_0^Tx e^{-st}dt$ |
| Convolution | $(x*h)(t)$ | $X(s)H(s)$ |

> [!example] Worked Example 14.13 — extending the table: a logarithmic pair
> Verify entry 36, $\dfrac{e^{-bt}-e^{-at}}{t}\leftrightarrow\ln\dfrac{s+a}{s+b}$.
>
> The numerator $x(t)=e^{-bt}-e^{-at}$ has $X(p)=\frac{1}{p+b}-\frac{1}{p+a}$, and $x(t)/t\to a-b$ (finite) as $t\to0$. Dividing by $t$:
>
> $$
> \int_s^\infty\left(\frac{1}{p+b}-\frac{1}{p+a}\right)dp=\Big[\ln\frac{p+b}{p+a}\Big]_s^\infty=0-\ln\frac{s+b}{s+a}=\ln\frac{s+a}{s+b}.
> $$
>
> **Corollary (Frullani integral).** Letting $s\to0$: $\int_0^\infty\frac{e^{-bt}-e^{-at}}{t}dt=\ln\frac ab$. For $a=2$, $b=1$ this is $\ln2=0.693147$; numerical quadrature gives $0.693147$. ✓

> [!example] Worked Example 14.14 — $\sin^2\omega t$ two ways
> **Identity.** $\sin^2\omega t=\tfrac12(1-\cos2\omega t)\leftrightarrow\frac12\left(\frac1s-\frac{s}{s^2+4\omega^2}\right)=\frac{2\omega^2}{s(s^2+4\omega^2)}$.
>
> **Derivative rule.** Let $x=\sin^2\omega t$: $x(0^-)=0$, $\dot x=2\omega\sin\omega t\cos\omega t=\omega\sin2\omega t\leftrightarrow\frac{2\omega^2}{s^2+4\omega^2}$. Since $\mathcal{L}\{\dot x\}=sX$, dividing by $s$ gives the same result. ✓
>
> The second route is often faster when the derivative is simpler than the function.

The table shows a pattern worth internalising: **the denominator of $X(s)$ encodes the modes, the numerator encodes how much of each mode is present.** The roots of the denominator — the poles — determine the *shape* of the signal: a real pole $-a$ gives $e^{-at}$, a pair $-\alpha\pm j\omega$ gives $e^{-\alpha t}\cos(\omega t+\phi)$, a repeated pole multiplies by powers of $t$.

```desmos-graph
left=0; right=6; bottom=-1.5; top=2.2
---
y=e^{-0.8x}|blue|0<x<6
y=e^{-0.4x}\cos(3x)|red|0<x<6
y=\cos(3x)|grey|dashed|0<x<6
y=e^{0.15x}\cos(3x)|orange|0<x<6
y=x e^{-x}|green|0<x<6
```
*Pole location and time shape: real pole at $-0.8$ (blue); poles $-0.4\pm j3$ (red, decaying oscillation); poles $\pm j3$ (grey dashed, sustained); poles $0.15\pm j3$ (orange, growing); double pole at $-1$ giving $te^{-t}$ (green).*

> [!warning] Common Pitfall — misreading $\omega$ in the denominator
> In $\dfrac{1}{s^2+9}$ the frequency is $\omega=3$, not 9, and the numerator needed for a pure sine is $\omega=3$: $\dfrac{1}{s^2+9}=\dfrac13\cdot\dfrac{3}{s^2+9}\leftrightarrow\dfrac13\sin3t$. Forgetting the $\frac13$ is one of the most frequent errors in inversion. Always check the value at $t=0^+$ or the DC gain as a sanity test.

> [!tip] Engineering Note — where the exotic entries appear
> Entries 9–11 and 39–40 (fractional powers, $\sqrt s$, $\operatorname{erfc}$) describe diffusion: carrier diffusion in a diode, heat spreading into a heatsink, and the step response of a long $RC$ line (on-chip interconnect, the original transatlantic telegraph cable). Entry 38 ($J_0$) appears in cylindrical problems such as skin effect in round wires. Entries 41–44 describe clock and switching waveforms in power electronics.

## 14.4 Inverse Laplace transforms

Solving a problem in the $s$-domain is only half the job; the answer must be returned to the time domain, where voltages and currents are measured. The inverse transform does this. In practice we almost never evaluate an inversion integral: we rearrange $X(s)$ into a sum of table entries and read off the time function, using uniqueness (Lerch's theorem, Section 14.1) to guarantee the answer is the only one.

### Definition and the inversion integral

> [!definition] Inverse Laplace transform
> If $X(s)=\mathcal{L}\{x(t)\}$, then $x(t)=\mathcal{L}^{-1}\{X(s)\}$ is the **inverse Laplace transform** of $X$. By Lerch's theorem it is unique for $t>0$ apart from values at isolated points. The inverse is linear: $\mathcal{L}^{-1}\{aX+bY\}=a\,x+b\,y$.

There is an explicit formula. Writing $s=\sigma+j\omega$ with fixed $\sigma$ in the ROC, the Laplace transform of $x(t)$ is the Fourier transform of $x(t)e^{-\sigma t}u(t)$. Applying the inverse Fourier transform of [Ch. 16](16-fourier-transforms.md) and multiplying by $e^{\sigma t}$ gives the following.

> [!theorem] Bromwich inversion integral
> If $X(s)$ is the transform of $x(t)$ with abscissa of convergence $\sigma_c$, then for any $\gamma>\sigma_c$,
>
> $$
> x(t)=\frac{1}{2\pi j}\int_{\gamma-j\infty}^{\gamma+j\infty}X(s)\,e^{st}\,ds,\qquad t>0,
> $$
>
> where the integral is along the vertical line $\operatorname{Re}\{s\}=\gamma$ (to the right of all singularities). At a jump discontinuity the integral returns the midpoint $\tfrac12[x(t^-)+x(t^+)]$; for $t<0$ it returns 0.

*Derivation sketch.* $X(\gamma+j\omega)=\int_0^\infty\big[x(t)e^{-\gamma t}\big]e^{-j\omega t}dt$ is the Fourier transform of $g(t)=x(t)e^{-\gamma t}u(t)$, which is absolutely integrable because $\gamma>\sigma_c$. Fourier inversion gives $g(t)=\frac{1}{2\pi}\int_{-\infty}^\infty X(\gamma+j\omega)e^{j\omega t}d\omega$. Multiply both sides by $e^{\gamma t}$ and change variable to $s=\gamma+j\omega$, $ds=j\,d\omega$. $\square$

This is a contour integral in the complex plane. For $t>0$ the line can be closed by a large semicircle to the left; when $X(s)$ is rational the semicircle contributes nothing, and the **residue theorem** turns the integral into a finite sum:

$$
x(t)=\sum_{\text{poles }p_k\text{ of }X}\operatorname{Res}_{s=p_k}\left[X(s)e^{st}\right],\qquad t>0.
$$

Closing the contour, Jordan's lemma, and the computation of residues are the subject of [Ch. 18](18-complex-analysis.md). The partial-fraction method of Section 14.5 is exactly this residue sum carried out with algebra, so for rational functions nothing is lost by using tables.

### Inversion by table and properties

The working method is:

1. **Separate delays.** Collect terms by their factor $e^{-sT}$; invert each coefficient and then apply the delay $t\to t-T$ with $u(t-T)$.
2. **Reduce to proper rational form.** If the degree of the numerator is not less than that of the denominator, divide (polynomial long division); the polynomial part inverts to impulses and their derivatives.
3. **Expand in partial fractions** (Section 14.5) or, for a single quadratic factor, **complete the square**.
4. **Read the table**, using the $s$-shift to absorb $s+\alpha$ and adjusting constants so each term matches a table entry exactly.
5. **Check** at least one value: $x(0^+)=\lim_{s\to\infty}sX(s)$, or the final value, or $X$ at a convenient real $s$.

> [!example] Worked Example 14.15 — completing the square
> Find $\mathcal{L}^{-1}\left\{\dfrac{2s+12}{s^2+2s+5}\right\}$.
>
> **Complete the square.** $s^2+2s+5=(s+1)^2+4$: the poles are $-1\pm j2$, so expect $e^{-t}$ times sinusoids of frequency 2 rad/s.
>
> **Rewrite the numerator in terms of $s+1$.** $2s+12=2(s+1)+10$:
>
> $$
> X(s)=2\,\frac{s+1}{(s+1)^2+2^2}+5\cdot\frac{2}{(s+1)^2+2^2}.
> $$
>
> **Table entries 23 and 24.**
>
> $$
> x(t)=e^{-t}\left(2\cos2t+5\sin2t\right)=\sqrt{29}\,e^{-t}\cos(2t-1.1903),
> $$
>
> using $A\cos\theta+B\sin\theta=\sqrt{A^2+B^2}\cos(\theta-\arctan(B/A))$ with $\arctan(5/2)=1.1903$ rad.
>
> **Check.** $x(0^+)=2$, and $\lim_{s\to\infty}sX=\lim\frac{2s^2+12s}{s^2+2s+5}=2$. ✓ Also $X(0)=12/5=2.4$; numerically $\int_0^\infty e^{-t}(2\cos2t+5\sin2t)dt=2\cdot\frac15+5\cdot\frac25=2.4$. ✓

> [!example] Worked Example 14.16 — a delayed response
> Find $\mathcal{L}^{-1}\left\{\dfrac{3e^{-2s}}{s(s+3)}\right\}$.
>
> **Invert the coefficient.** By entry 15 with $a=3$: $\dfrac{3}{s(s+3)}\leftrightarrow1-e^{-3t}$.
>
> **Apply the delay** ($T=2$): replace $t$ by $t-2$ and multiply by $u(t-2)$:
>
> $$
> x(t)=\left(1-e^{-3(t-2)}\right)u(t-2).
> $$
>
> This is the response of a first-order system with time constant $\tfrac13$ s to a step applied at $t=2$ s. A common slip is to write $(1-e^{-3t})u(t-2)$, which jumps discontinuously at $t=2$ to $1-e^{-6}$ instead of starting smoothly from 0.

> [!example] Worked Example 14.17 — repeated quadratic factor via the table
> Find $\mathcal{L}^{-1}\left\{\dfrac{s}{(s^2+4)^2}\right\}$ and $\mathcal{L}^{-1}\left\{\dfrac{1}{(s^2+4)^2}\right\}$.
>
> **First.** Entry 27 with $\omega=2$: $\dfrac{4s}{(s^2+4)^2}\leftrightarrow t\sin2t$, hence $\dfrac{s}{(s^2+4)^2}\leftrightarrow\dfrac{t\sin2t}{4}$.
>
> **Second.** Entry 29 with $\omega=2$: $\dfrac{16}{(s^2+4)^2}\leftrightarrow\sin2t-2t\cos2t$, hence
>
> $$
> \frac{1}{(s^2+4)^2}\leftrightarrow\frac{\sin2t-2t\cos2t}{16}.
> $$
>
> **Cross-check with the integral rule.** $\frac{1}{(s^2+4)^2}=\frac{1}{s}\cdot\frac{s}{(s^2+4)^2}$, so the answer must equal $\int_0^t\frac{\tau\sin2\tau}{4}d\tau=\frac14\left[\frac{\sin2\tau}{4}-\frac{\tau\cos2\tau}{2}\right]_0^t=\frac{\sin2t-2t\cos2t}{16}$. ✓

> [!example] Worked Example 14.18 — several delays: a triangular pulse
> Find and sketch $\mathcal{L}^{-1}\left\{\dfrac{1-2e^{-s}+e^{-2s}}{s^2}\right\}$, and identify the waveform.
>
> **Separate by delay.** The coefficient of each exponential is $1/s^2\leftrightarrow t$, so
>
> $$
> x(t)=t\,u(t)-2(t-1)\,u(t-1)+(t-2)\,u(t-2).
> $$
>
> **Evaluate piece by piece.**
>
> - $0\le t<1$: $x=t$ (rising ramp).
> - $1\le t<2$: $x=t-2(t-1)=2-t$ (falling ramp).
> - $t\ge2$: $x=t-2t+2+t-2=0$.
>
> The signal is a **triangular pulse** of height 1 at $t=1$ and base 2. Note also that $\frac{(1-e^{-s})^2}{s^2}=\left(\frac{1-e^{-s}}{s}\right)^2$: the triangle is the convolution of a unit rectangular pulse with itself (Section 14.6).
>
> **Check.** Area: $X(0)=\lim_{s\to0}\frac{(1-e^{-s})^2}{s^2}=1$, and the triangle's area is $\tfrac12\cdot2\cdot1=1$ ✓.

```desmos-graph
left=-0.5; right=3.5; bottom=-0.3; top=1.3
---
y=x|blue|0<x<1
y=2-x|blue|1<x<2
y=0|blue|2<x<3.5
y=x|grey|dashed|1<x<3
y=x-2(x-1)|grey|dotted|2<x<3
```
*Example 14.18: the ramp $t$ (grey dashed beyond $t=1$) is bent down by $-2(t-1)u(t-1)$ (grey dotted beyond $t=2$) and levelled by $(t-2)u(t-2)$, producing the triangular pulse (blue).*

### Non-rational transforms

Transforms involving logarithms, inverse tangents or square roots cannot be split into partial fractions. Two tricks often work.

- **Differentiate in $s$.** If $X'(s)$ is rational, invert it and use $\mathcal{L}^{-1}\{X'\}=-t\,x(t)$, so $x(t)=-\frac1t\mathcal{L}^{-1}\{X'(s)\}$.
- **Recognise a scaled or shifted table entry**, e.g. $1/\sqrt{s+a}\leftrightarrow e^{-at}/\sqrt{\pi t}$ by the $s$-shift of entry 10.

> [!example] Worked Example 14.19 — inverting an arctangent
> Find $\mathcal{L}^{-1}\left\{\arctan\dfrac{2}{s}\right\}$.
>
> **Differentiate.** $\dfrac{d}{ds}\arctan\dfrac{2}{s}=\dfrac{1}{1+4/s^2}\cdot\left(-\dfrac{2}{s^2}\right)=-\dfrac{2}{s^2+4}\leftrightarrow-\sin2t$.
>
> **Undo the $t$-multiplication.** Since $-t\,x(t)=-\sin2t$,
>
> $$
> x(t)=\frac{\sin2t}{t},
> $$
>
> agreeing with entry 35. **Check.** $x(0^+)=2$; and $\lim_{s\to\infty}s\arctan(2/s)=\lim s\cdot\frac2s=2$. ✓

> [!warning] Common Pitfall — not every function of $s$ is a transform
> By the existence theorem, a genuine transform of an ordinary function satisfies $X(s)\to0$ as $\operatorname{Re}\{s\}\to\infty$. Functions such as $s$, $s^2/(s+1)$ or $e^{+sT}$ do not; the first two contain impulses or their derivatives (handled after long division), and $e^{sT}$ would be an *advance*, which no causal signal has. If an analysis produces $e^{+sT}$, a sign error has been made.

> [!tip] Engineering Note — numerical inversion
> When $X(s)$ is known only numerically (e.g. from a measured impedance or a transmission-line model with $\sqrt s$), engineers use numerical Laplace inversion: the Talbot method deforms the Bromwich line into a contour where $e^{st}$ decays rapidly, and the Stehfest algorithm uses only real values of $X(s)$. Both are standard in circuit-simulation and geophysics codes ([Ch. 22](../part-7-numerical-computational/22-numerical-methods.md)).

## 14.5 Partial fractions

Almost every transform produced by a lumped circuit or a constant-coefficient ODE is a **rational function** — a ratio of polynomials in $s$:

$$
X(s)=\frac{N(s)}{D(s)}=\frac{b_ms^m+\cdots+b_1s+b_0}{a_ns^n+\cdots+a_1s+a_0}.
$$

Partial-fraction expansion rewrites it as a sum of simple terms, each of which is a table entry. Because each term corresponds to one *mode* (pole) of the system, the expansion is also a modal decomposition: it tells us how strongly each natural frequency is excited.

### Proper and improper rational functions

$X(s)$ is **strictly proper** if $m<n$. Only strictly proper functions are transforms of ordinary functions (Section 14.4). If $m\ge n$, long division gives

$$
X(s)=Q(s)+\frac{R(s)}{D(s)},\qquad \deg R<\deg D,
$$

and the polynomial quotient $Q(s)=c_0+c_1s+\cdots$ inverts to $c_0\delta(t)+c_1\dot\delta(t)+\cdots$. Physically, an improper transfer function is a differentiator: such terms appear for idealised elements (the current through an ideal capacitor driven by a voltage step), and they are a warning that the model has no high-frequency roll-off.

### The expansion theorem

Factor the denominator over $\mathbb{C}$: $D(s)=a_n\prod_k(s-p_k)^{r_k}$, where the $p_k$ are the distinct **poles** with multiplicities $r_k$, $\sum r_k=n$.

> [!theorem] Partial-fraction expansion
> A strictly proper rational function can be written uniquely as
>
> $$
> \frac{N(s)}{D(s)}=\sum_{k}\sum_{i=1}^{r_k}\frac{K_{k,i}}{(s-p_k)^{i}},
> $$
>
> with the coefficients of the highest power at each pole and the others given by
>
> $$
> K_{k,r_k}=\Big[(s-p_k)^{r_k}X(s)\Big]_{s=p_k},\qquad
> K_{k,r_k-j}=\frac{1}{j!}\,\frac{d^j}{ds^j}\Big[(s-p_k)^{r_k}X(s)\Big]_{s=p_k},\quad j=1,\ldots,r_k-1 .
> $$
>
> The inverse transform is then
>
> $$
> x(t)=\sum_k e^{p_kt}\sum_{i=1}^{r_k}K_{k,i}\,\frac{t^{\,i-1}}{(i-1)!},\qquad t>0 .
> $$

*Why the formulas work.* Multiply the expansion by $(s-p_k)^{r_k}$. Every term belonging to another pole acquires a factor $(s-p_k)^{r_k}$ and vanishes at $s=p_k$, while the terms belonging to $p_k$ become the polynomial $K_{k,r_k}+K_{k,r_k-1}(s-p_k)+\cdots+K_{k,1}(s-p_k)^{r_k-1}$. Evaluating it and its derivatives at $s=p_k$ isolates the coefficients one at a time — this is just Taylor's theorem. Existence and uniqueness follow from the fact that the $n$ functions $(s-p_k)^{-i}$ are linearly independent and $n$ coefficients match $n$ degrees of freedom of $N$. The inverse formula uses entry 14 with the $s$-shift. $\square$

For a **simple pole** ($r_k=1$), the coefficient $K_k$ is the **residue** of $X$ at $p_k$. Since $D(p_k)=0$, L'Hôpital's rule gives a useful alternative.

> [!theorem] Heaviside expansion formula (simple poles)
> If all poles $p_k$ of $N/D$ are simple,
>
> $$
> K_k=\left.\frac{N(s)}{D'(s)}\right|_{s=p_k},\qquad
> x(t)=\sum_k\frac{N(p_k)}{D'(p_k)}\,e^{p_kt}.
> $$

This avoids factoring $D$ symbolically: only its roots are needed, which a computer finds numerically.

### Case 1 — distinct real poles: the cover-up rule

For a simple real pole, $K_k=[(s-p_k)X(s)]_{s=p_k}$: "cover up" the factor $(s-p_k)$ in the denominator and substitute $s=p_k$ in what remains.

> [!example] Worked Example 14.20 — three real poles
> Invert $X(s)=\dfrac{3s+7}{(s+1)(s+2)(s+3)}$.
>
> **Form.** $X=\dfrac{K_1}{s+1}+\dfrac{K_2}{s+2}+\dfrac{K_3}{s+3}$.
>
> **Cover-up.**
>
> $$
> K_1=\left.\frac{3s+7}{(s+2)(s+3)}\right|_{s=-1}=\frac{4}{(1)(2)}=2,\quad
> K_2=\left.\frac{3s+7}{(s+1)(s+3)}\right|_{s=-2}=\frac{1}{(-1)(1)}=-1,\quad
> K_3=\left.\frac{3s+7}{(s+1)(s+2)}\right|_{s=-3}=\frac{-2}{(-2)(-1)}=-1.
> $$
>
> **Heaviside check.** $D=s^3+6s^2+11s+6$, $D'=3s^2+12s+11$: $D'(-1)=2$, $D'(-2)=-1$, $D'(-3)=2$, so $K_1=4/2=2$, $K_2=1/(-1)=-1$, $K_3=-2/2=-1$. ✓
>
> **Inverse.**
>
> $$
> x(t)=2e^{-t}-e^{-2t}-e^{-3t},\qquad t\ge0 .
> $$
>
> **Check.** $x(0^+)=2-1-1=0$, consistent with $\lim sX=\lim\frac{3s^2}{s^3}=0$. Also $X(0)=7/6$ and $\int_0^\infty x\,dt=2-\tfrac12-\tfrac13=\tfrac76$. ✓

### Case 2 — repeated poles

A pole of multiplicity $r$ needs $r$ terms. The top coefficient comes from the cover-up rule; the lower ones need derivatives (or, often more simply, comparison of coefficients or substitution of convenient values of $s$).

> [!example] Worked Example 14.21 — a triple pole
> Invert $X(s)=\dfrac{s^2+2s+5}{(s+1)^3(s+2)}$.
>
> **Form.** $X=\dfrac{A_3}{(s+1)^3}+\dfrac{A_2}{(s+1)^2}+\dfrac{A_1}{s+1}+\dfrac{B}{s+2}$.
>
> **Simple pole.** $B=\left.\dfrac{s^2+2s+5}{(s+1)^3}\right|_{s=-2}=\dfrac{5}{-1}=-5$.
>
> **Triple pole.** Let $g(s)=(s+1)^3X(s)=\dfrac{s^2+2s+5}{s+2}$. Division gives $g(s)=s+\dfrac{5}{s+2}$, which is easy to differentiate:
>
> $$
> g(-1)=-1+5=4,\qquad g'(s)=1-\frac{5}{(s+2)^2}\Rightarrow g'(-1)=-4,\qquad g''(s)=\frac{10}{(s+2)^3}\Rightarrow \tfrac12g''(-1)=5 .
> $$
>
> So $A_3=4$, $A_2=-4$, $A_1=5$.
>
> **Sanity check of coefficients.** For large $s$, $X\approx 1/s^2$, so the $1/s$ terms must cancel: $A_1+B=5-5=0$. ✓
>
> **Inverse** (entry 14):
>
> $$
> x(t)=\left(2t^2-4t+5\right)e^{-t}-5e^{-2t}.
> $$
>
> **Check.** $x(0^+)=0$ ✓ (matches $\lim sX=0$). At $s=0$: $X(0)=5/2$; from the expansion $4-4+5-\tfrac52=\tfrac52$. ✓

### Case 3 — complex-conjugate poles

If $D$ has real coefficients, complex poles occur in conjugate pairs $p,\bar p=-\alpha\pm j\beta$, and their coefficients are conjugates, $K$ and $\bar K$. Their combined contribution is real:

$$
Ke^{pt}+\bar K e^{\bar pt}=2\operatorname{Re}\{Ke^{pt}\}=2\lvert K\rvert e^{-\alpha t}\cos\big(\beta t+\angle K\big).
$$

> [!theorem] Complex-pole pair
> If $X(s)$ has simple poles at $-\alpha\pm j\beta$ with residue $K=[(s+\alpha-j\beta)X(s)]_{s=-\alpha+j\beta}$ at the upper pole, their contribution to $x(t)$ is
>
> $$
> 2\lvert K\rvert\,e^{-\alpha t}\cos(\beta t+\angle K).
> $$

Two equivalent working methods exist:

- **Complex residues** (above): compute one complex number per pair and write amplitude and phase directly. Efficient, and natural for numerical work.
- **Real quadratic factor:** keep $(s+\alpha)^2+\beta^2$ intact with numerator $Bs+C$, find $B,C$ by comparing coefficients, then complete the square as in Example 14.15. Avoids complex arithmetic and is often quicker by hand.

> [!example] Worked Example 14.22 — one real pole and a complex pair, both methods
> Invert $X(s)=\dfrac{s+3}{(s+1)(s^2+4s+13)}$.
>
> The quadratic has roots $-2\pm j3$.
>
> **Real pole (cover-up).** $A=\left.\dfrac{s+3}{s^2+4s+13}\right|_{s=-1}=\dfrac{2}{10}=\dfrac15$.
>
> **Method 1 — real quadratic.** Write $X=\dfrac{1/5}{s+1}+\dfrac{Bs+C}{s^2+4s+13}$ and clear denominators:
>
> $$
> s+3=\tfrac15(s^2+4s+13)+(Bs+C)(s+1).
> $$
>
> $s^2$: $0=\tfrac15+B\Rightarrow B=-\tfrac15$. Constant: $3=\tfrac{13}{5}+C\Rightarrow C=\tfrac25$. Then
>
> $$
> \frac{-\frac15s+\frac25}{(s+2)^2+9}=-\frac15\,\frac{(s+2)-4}{(s+2)^2+9}=-\frac15\,\frac{s+2}{(s+2)^2+9}+\frac{4}{15}\,\frac{3}{(s+2)^2+9},
> $$
>
> so
>
> $$
> x(t)=\tfrac15e^{-t}+e^{-2t}\left(-\tfrac15\cos3t+\tfrac{4}{15}\sin3t\right).
> $$
>
> **Method 2 — complex residue.** At $p=-2+j3$:
>
> $$
> K=\left.\frac{s+3}{(s+1)(s+2+j3)}\right|_{s=-2+j3}=\frac{1+j3}{(-1+j3)(j6)}=\frac{1+j3}{-18-j6}=-0.1-j0.1333 .
> $$
>
> $\lvert K\rvert=\tfrac16$, $\angle K=-2.2143$ rad ($-126.87^\circ$). The pair contributes $\tfrac13e^{-2t}\cos(3t-2.2143)$. Expanding, $\tfrac13[\cos(2.2143)\cos3t+\sin(2.2143)\sin3t]=\tfrac13[-0.6\cos3t+0.8\sin3t]=-\tfrac15\cos3t+\tfrac4{15}\sin3t$. ✓ Both methods agree.
>
> **Check.** $x(0^+)=\tfrac15-\tfrac15=0$ ✓; at $t=0.7$ both forms give $0.18098$.

### Case 4 — improper functions

> [!example] Worked Example 14.23 — long division first
> Invert $X(s)=\dfrac{s^3+2s^2+3s+1}{s^2+s+1}$.
>
> **Divide.** $s^3+2s^2+3s+1=(s+1)(s^2+s+1)+s$, so
>
> $$
> X(s)=s+1+\frac{s}{s^2+s+1}.
> $$
>
> **Proper part.** $s^2+s+1=(s+\tfrac12)^2+\tfrac34$, and $s=(s+\tfrac12)-\tfrac12$:
>
> $$
> \frac{s}{s^2+s+1}=\frac{s+\frac12}{(s+\frac12)^2+\frac34}-\frac{1}{\sqrt3}\cdot\frac{\frac{\sqrt3}{2}}{(s+\frac12)^2+\frac34}.
> $$
>
> **Inverse.**
>
> $$
> x(t)=\dot\delta(t)+\delta(t)+e^{-t/2}\left(\cos\frac{\sqrt3}{2}t-\frac{1}{\sqrt3}\sin\frac{\sqrt3}{2}t\right)u(t).
> $$
>
> The impulse terms are a genuine part of the answer. If $X$ were a transfer function, the system would differentiate its input — physically unrealisable at high frequencies, which is why practical transfer functions are proper.

### Choosing a method, and checking

| Situation | Fastest method |
|---|---|
| Distinct real poles | Cover-up (or $N/D'$) |
| Repeated pole, low multiplicity | Cover-up for top term, then compare coefficients or substitute a value of $s$ |
| Repeated pole, high multiplicity | Derivative formula, after simplifying $g(s)=(s-p)^rX$ |
| One complex pair, hand work | Real quadratic + completing the square |
| Several complex pairs, numerical | Complex residues, amplitude–phase form |
| $\deg N\ge\deg D$ | Long division first |

> [!warning] Common Pitfall — too few terms for a repeated pole
> Writing $\dfrac{1}{(s+1)^2(s+2)}=\dfrac{A}{(s+1)^2}+\dfrac{B}{s+2}$ is wrong: there are three unknowns' worth of freedom ($\deg D=3$) and the term $\dfrac{C}{s+1}$ is missing. The correct expansion is $\dfrac{1}{(s+1)^2}-\dfrac{1}{s+1}+\dfrac{1}{s+2}$. A quick count — the number of coefficients must equal $\deg D$ — catches this every time.

> [!tip] Engineering Note — residues as modal amplitudes
> In a circuit response, each coefficient $K_k$ is the initial amplitude of mode $e^{p_kt}$. A pole with a tiny residue is barely excited and can often be dropped (*dominant-pole approximation*); a pole and zero that nearly coincide (a "doublet") produce a small residue and a slow, low-amplitude tail that can still matter for settling time in precision amplifiers. Computer tools (`residue` in MATLAB, `scipy.signal.residue` in Python) implement the Heaviside formula with numerical root-finding.

## 14.6 Convolution theorem

Partial fractions invert *sums*; the convolution theorem inverts *products*. It answers the question "if I know the time functions for $F(s)$ and $G(s)$, what time function has transform $F(s)G(s)$?" The answer — their convolution — is also the central formula of linear system theory: the output of a linear time-invariant system is the convolution of its input with its impulse response.

### Definition

> [!definition] Convolution of causal signals
> For functions $f,g$ defined on $t\ge0$ (zero for $t<0$), the **convolution** is
>
> $$
> (f*g)(t)=\int_{0^-}^{t}f(\tau)\,g(t-\tau)\,d\tau,\qquad t\ge0 .
> $$

For general (two-sided) signals the integral runs over $(-\infty,\infty)$; for causal signals the integrand vanishes outside $0\le\tau\le t$, giving the finite limits above. The operation "flip $g$, slide it by $t$, multiply by $f$, integrate" is linear in each argument and

- commutative: $f*g=g*f$ (substitute $\tau'=t-\tau$);
- associative: $(f*g)*h=f*(g*h)$;
- distributive: $f*(g+h)=f*g+f*h$;
- has identity $\delta$: $f*\delta=f$, and $f*\delta(t-T)=f(t-T)u(t-T)$;
- $f*u$ is the running integral: $(f*u)(t)=\int_{0}^{t}f(\tau)d\tau$.

### The theorem

> [!theorem] Convolution theorem
> If $f\leftrightarrow F(s)$ and $g\leftrightarrow G(s)$, both of exponential order, then
>
> $$
> \mathcal{L}\{(f*g)(t)\}=F(s)\,G(s),
> $$
>
> on the intersection of the ROCs. Equivalently, $\mathcal{L}^{-1}\{F(s)G(s)\}=\int_0^tf(\tau)g(t-\tau)d\tau$.

*Proof.* Write the convolution with a step so that the upper limit becomes infinite:

$$
\mathcal{L}\{f*g\}=\int_{0}^{\infty}e^{-st}\int_{0}^{\infty}f(\tau)\,g(t-\tau)\,u(t-\tau)\,d\tau\,dt .
$$

For $s$ in both ROCs the double integral converges absolutely, so the order of integration can be exchanged (Fubini's theorem, [Ch. 6](../part-2-calculus/06-multivariable-calculus.md)):

$$
=\int_0^\infty f(\tau)\left[\int_{\tau}^{\infty}g(t-\tau)e^{-st}dt\right]d\tau
=\int_0^\infty f(\tau)\,e^{-s\tau}G(s)\,d\tau=F(s)G(s),
$$

where the inner integral was evaluated with $t'=t-\tau$ exactly as in the time-shift theorem. $\square$

The region of integration is the wedge $0\le\tau\le t$ in the $(\tau,t)$ plane; swapping the order means integrating first over $t$ from $\tau$ to $\infty$ — a picture worth drawing once.

### Impulse response and the system equation

Consider an LTI system initially at rest. If the input $x(t)=\delta(t)$ produces the output $h(t)$ — the **impulse response** — then by linearity and time invariance the input $x(\tau)\delta(t-\tau)\,d\tau$ produces $x(\tau)h(t-\tau)\,d\tau$, and summing over all $\tau$:

$$
y(t)=\int_{0^-}^{t}x(\tau)\,h(t-\tau)\,d\tau=(x*h)(t)\quad\Longleftrightarrow\quad Y(s)=H(s)X(s).
$$

The function $H(s)=\mathcal{L}\{h(t)\}$ is the **transfer function** (Section 14.9). Convolution in time — an integral over the whole past of the input — becomes multiplication in $s$. This is the deepest reason the Laplace transform simplifies linear-system analysis.

A related form, **Duhamel's integral**, uses the step response $a(t)=\int_0^th\,d\tau$ instead: since $H X=(sA)X=A\,(sX)$, the response is $y(t)=x(0^-)a(t)+\int_0^t\dot x(\tau)a(t-\tau)d\tau$. It is useful when the step response has been measured experimentally.

> [!example] Worked Example 14.24 — inverting a product by convolution
> Find $\mathcal{L}^{-1}\left\{\dfrac{1}{s^2(s+1)}\right\}$ by convolution and confirm by partial fractions.
>
> **Factor.** $F=1/s^2\leftrightarrow t$, $G=1/(s+1)\leftrightarrow e^{-t}$.
>
> **Convolve** (choose the order that makes the integral easiest; here put the exponential in the shifted slot):
>
> $$
> \int_0^t\tau\,e^{-(t-\tau)}d\tau=e^{-t}\int_0^t\tau e^{\tau}d\tau=e^{-t}\Big[(\tau-1)e^{\tau}\Big]_0^t=e^{-t}\big[(t-1)e^t+1\big]=t-1+e^{-t}.
> $$
>
> **Partial fractions.** $\dfrac{1}{s^2(s+1)}=\dfrac{1}{s^2}-\dfrac1s+\dfrac{1}{s+1}\leftrightarrow t-1+e^{-t}$. ✓ (Table entry 18 with $a=1$.)

> [!example] Worked Example 14.25 — a repeated quadratic by convolution
> Find $\mathcal{L}^{-1}\left\{\dfrac{1}{(s^2+1)^2}\right\}$.
>
> $\dfrac{1}{(s^2+1)^2}=\dfrac{1}{s^2+1}\cdot\dfrac{1}{s^2+1}\leftrightarrow\sin t*\sin t$. Using $\sin A\sin B=\tfrac12[\cos(A-B)-\cos(A+B)]$:
>
> $$
> \int_0^t\sin\tau\,\sin(t-\tau)\,d\tau=\frac12\int_0^t\big[\cos(2\tau-t)-\cos t\big]d\tau
> =\frac12\left[\frac{\sin(2\tau-t)}{2}\right]_0^t-\frac{t\cos t}{2}
> =\frac{\sin t}{2}-\frac{t\cos t}{2}.
> $$
>
> So $\mathcal{L}^{-1}\{(s^2+1)^{-2}\}=\tfrac12(\sin t-t\cos t)$, entry 29 with $\omega=1$. ✓ Partial fractions cannot simplify $(s^2+1)^{-2}$ further over the reals, so convolution (or the table) is the natural route.

> [!example] Worked Example 14.26 — pulse response of an RC low-pass filter
> An $RC$ low-pass filter with $RC=1$ ms has impulse response $h(t)=\frac{1}{RC}e^{-t/RC}u(t)$. Find the output for a 1 V input pulse lasting 1 ms. Work in milliseconds so $h(t)=e^{-t}$.
>
> **Time domain.** The input is $x(t)=u(t)-u(t-1)$. For $0\le t<1$ the pulse overlaps $[0,t]$:
>
> $$
> y(t)=\int_0^t1\cdot e^{-(t-\tau)}d\tau=1-e^{-t}.
> $$
>
> For $t\ge1$ only $[0,1]$ contributes:
>
> $$
> y(t)=\int_0^1e^{-(t-\tau)}d\tau=e^{-t}(e-1).
> $$
>
> **$s$-domain.** $Y=H X=\dfrac{1}{s+1}\cdot\dfrac{1-e^{-s}}{s}=(1-e^{-s})\left(\dfrac1s-\dfrac{1}{s+1}\right)$, so
>
> $$
> y(t)=(1-e^{-t})u(t)-\big(1-e^{-(t-1)}\big)u(t-1).
> $$
>
> For $t\ge1$: $1-e^{-t}-1+e^{-t+1}=e^{-t}(e-1)$. ✓ Both routes agree; the $s$-domain route needed no case analysis.
>
> **Numbers.** Peak at $t=1$ ms: $y=1-e^{-1}=0.632$ V; at $t=3$ ms: $y=(e-1)e^{-3}=0.0855$ V.

```desmos-graph
left=-0.5; right=6; bottom=-0.2; top=1.2
---
y=1|grey|dashed|0<x<1
y=1-e^{-x}|blue|0<x<1
y=(e-1)e^{-x}|blue|1<x<6
(1,0.632)|label:peak 0.632 V
```
*Response of an RC low-pass filter ($RC=1$ ms) to a 1 V, 1 ms pulse (grey): the capacitor charges toward 1 V while the pulse lasts, then discharges exponentially (Example 14.26; time in ms).*

> [!example] Worked Example 14.27 — a Volterra integral equation
> Solve $y(t)=t+\displaystyle\int_0^t\sin(t-\tau)\,y(\tau)\,d\tau$.
>
> The integral is $\sin*y$. Transforming, with $\mathcal{L}\{\sin t\}=1/(s^2+1)$:
>
> $$
> Y=\frac{1}{s^2}+\frac{Y}{s^2+1}\quad\Longrightarrow\quad Y\,\frac{s^2}{s^2+1}=\frac{1}{s^2}\quad\Longrightarrow\quad Y=\frac{s^2+1}{s^4}=\frac{1}{s^2}+\frac{1}{s^4}.
> $$
>
> $$
> y(t)=t+\frac{t^3}{6}.
> $$
>
> **Check by substitution.** $\int_0^t\sin(t-\tau)\left(\tau+\tfrac{\tau^3}{6}\right)d\tau$: using $\sin*\tau=t-\sin t$ and $\sin*\tfrac{\tau^3}{6}=\tfrac{t^3}{6}-t+\sin t$ (from $\frac{1}{s^4(s^2+1)}=\frac1{s^4}-\frac1{s^2}+\frac1{s^2+1}$), the sum is $t^3/6$, and $t+t^3/6=y$. ✓
>
> Integral equations of this convolution type arise in models with memory: hereditary (viscoelastic) materials, and circuits described by their impulse response rather than by an ODE.

> [!warning] Common Pitfall — the product of transforms is not the transform of the product
> $\mathcal{L}\{f(t)g(t)\}\ne F(s)G(s)$. Multiplication in time corresponds to a *complex convolution* in $s$ (a contour integral), which is rarely used in circuit work. For example $\mathcal{L}\{t\cdot t\}=2/s^3$, while $\frac{1}{s^2}\cdot\frac{1}{s^2}=\frac{1}{s^4}\leftrightarrow t^3/6$.

> [!tip] Engineering Note — convolution in practice
> Measuring $h(t)$ directly with an impulse is difficult, but its integral, the step response, is easy to record on an oscilloscope. Duhamel's integral then predicts the response to any input. In digital signal processing the same idea, in discrete form, is the FIR filter: the output is a weighted sum of past inputs with the impulse response as weights ([Ch. 17](17-z-transforms.md)).

## 14.7 Initial- and final-value theorems

Often an engineer does not need the whole response — only where it starts and where it ends. Does the output voltage jump when the switch closes? What is the steady-state error of a control loop? What is the DC gain? The initial- and final-value theorems read these numbers directly from $X(s)$, with no inversion. They are also the quickest way to check a partial-fraction result.

### Statements

> [!theorem] Initial-value theorem (IVT)
> If $x(t)$ and $\dot x(t)$ (excluding any impulse at the origin) are transformable and $X(s)$ is strictly proper, then
>
> $$
> x(0^+)=\lim_{s\to\infty}sX(s),
> $$
>
> the limit taken along the positive real axis (or any path with $\operatorname{Re}\{s\}\to\infty$).

> [!theorem] Final-value theorem (FVT)
> If $x(t)$ and $\dot x(t)$ are transformable and **all poles of $sX(s)$ lie in the open left half-plane** $\operatorname{Re}\{s\}<0$, then $x(\infty)=\lim_{t\to\infty}x(t)$ exists and
>
> $$
> \lim_{t\to\infty}x(t)=\lim_{s\to0}sX(s).
> $$

The condition on the poles of $sX(s)$ means: $X(s)$ may have at most a *simple* pole at the origin (cancelled by the factor $s$), and every other pole must have negative real part. Without it, the right-hand side still produces a number, but that number is meaningless.

### Proofs

Both theorems come from the derivative rule. Split the derivative integral at the origin:

$$
\int_{0^-}^{\infty}\dot x(t)e^{-st}dt=sX(s)-x(0^-).
$$

The part from $0^-$ to $0^+$ picks up the jump $x(0^+)-x(0^-)$ (the impulse in $\dot x$), so

$$
sX(s)-x(0^+)=\int_{0^+}^{\infty}\dot x(t)e^{-st}dt .
$$

**IVT.** As $s\to\infty$ along the real axis, $e^{-st}\to0$ for every $t>0$, and (dominated convergence, using exponential order) the right side tends to 0. Hence $sX(s)\to x(0^+)$. If $X$ were not strictly proper, $x$ would contain an impulse at the origin, $x(0^+)$ would not describe it, and $sX\to\infty$.

**FVT.** As $s\to0$, $e^{-st}\to1$ and the right side tends to $\int_{0^+}^{\infty}\dot x\,dt=x(\infty)-x(0^+)$ — *provided* this integral converges, which the pole condition guarantees (every mode of $\dot x$ then decays exponentially). Adding $x(0^+)$ to both sides gives $\lim_{s\to0}sX=x(\infty)$. $\square$

In terms of partial fractions the FVT is transparent: with all poles of $X$ except a simple one at $s=0$ in the left half-plane, every exponential term dies out and only the constant $K_0$ from the term $K_0/s$ survives. That constant is $\lim_{s\to0}sX(s)$ by the cover-up rule.

### Higher initial derivatives

Applying the IVT to $\dot x$ (whose transform, for $t>0$, is $sX-x(0^+)$) gives

$$
\dot x(0^+)=\lim_{s\to\infty}s\big[sX(s)-x(0^+)\big],\qquad
\ddot x(0^+)=\lim_{s\to\infty}s\big[s^2X(s)-sx(0^+)-\dot x(0^+)\big],
$$

and so on. For a rational $X=N/D$ with $\deg D-\deg N=r$ (the **relative degree**), $x$ and its first $r-2$ derivatives start at zero, and $x^{(r-1)}(0^+)=b_m/a_n$, the ratio of leading coefficients. A high relative degree means a response that leaves zero gently — a smooth start.

> [!example] Worked Example 14.28 — start, slope and end without inversion
> For $X(s)=\dfrac{2s+3}{s(s+1)(s+3)}$ find $x(0^+)$, $\dot x(0^+)$ and $x(\infty)$; then check against the inverse.
>
> **Initial value.** $sX=\dfrac{2s+3}{(s+1)(s+3)}\to\dfrac{2s}{s^2}\to0$. So $x(0^+)=0$.
>
> **Initial slope.** $\dot x(0^+)=\lim s\cdot sX=\lim\dfrac{s(2s+3)}{(s+1)(s+3)}=2$. (Relative degree 2, leading-coefficient ratio $2/1$.)
>
> **Final value.** Poles of $sX$ are $-1,-3$: both in the left half-plane, so the FVT applies: $x(\infty)=\dfrac{3}{(1)(3)}=1$.
>
> **Check.** Partial fractions give $X=\dfrac1s-\dfrac{1/2}{s+1}-\dfrac{1/2}{s+3}$, so $x(t)=1-\tfrac12e^{-t}-\tfrac12e^{-3t}$: $x(0)=0$ ✓, $\dot x(0)=\tfrac12+\tfrac32=2$ ✓, $x(\infty)=1$ ✓.

```desmos-graph
left=-0.2; right=6; bottom=-0.1; top=1.3
---
y=1-0.5e^{-x}-0.5e^{-3x}|blue|0<x<6
y=1|red|dashed
y=2x|green|dotted|0<x<0.6
```
*$x(t)=1-\tfrac12e^{-t}-\tfrac12e^{-3t}$ (blue) starts at $x(0^+)=0$ with slope 2 (green dotted tangent) and settles to $x(\infty)=1$ (red dashed) — all three read from $X(s)$ in Example 14.28.*

> [!example] Worked Example 14.29 — when the FVT lies
> Apply $\lim_{s\to0}sX(s)$ blindly to each transform, and compare with the truth.
>
> | $X(s)$ | $\lim_{s\to0}sX$ | true $x(t)$ | verdict |
> |---|---|---|---|
> | $\dfrac{1}{s-2}$ | $0$ | $e^{2t}\to\infty$ | pole at $+2$ in RHP: FVT invalid |
> | $\dfrac{\omega}{s^2+\omega^2}$ | $0$ | $\sin\omega t$, no limit | poles $\pm j\omega$ on axis: FVT invalid |
> | $\dfrac{1}{s(s^2+1)}$ | $1$ | $1-\cos t$, no limit | $sX$ has poles $\pm j$: FVT invalid |
> | $\dfrac{1}{s^2}$ | $\infty$ | $t\to\infty$ | double pole at 0: FVT invalid (but correctly signals unboundedness) |
> | $\dfrac{4}{s(s+2)}$ | $2$ | $2-2e^{-2t}\to2$ | valid ✓ |
>
> The third row is the dangerous one: the theorem returns a finite, plausible number (the *average* of $1-\cos t$) for a signal that never settles. Always check the pole condition first.

> [!example] Worked Example 14.30 — an impulse at the origin
> $X(s)=\dfrac{s+2}{s+1}$. Find $x(0^+)$.
>
> $X$ is not strictly proper; $sX\to\infty$. Divide first: $X=1+\dfrac{1}{s+1}$, so $x(t)=\delta(t)+e^{-t}u(t)$. The IVT applies to the strictly proper part: $x(0^+)=\lim_{s\to\infty}\dfrac{s}{s+1}=1$. The impulse lives at $t=0$ itself, not at $0^+$, and must be reported separately.

> [!example] Worked Example 14.31 — DC gain and initial slope of a step response
> An amplifier has transfer function $H(s)=\dfrac{10(s+5)}{s^2+3s+50}$. For a unit-step input, find the initial value, initial slope and final value of the output.
>
> The output is $Y(s)=H(s)/s$, so $sY=H(s)$.
>
> - $y(0^+)=\lim_{s\to\infty}H(s)=0$ (relative degree 1 of $H$, so the step response is continuous).
> - $\dot y(0^+)=\lim_{s\to\infty}sH(s)=10$ V/s per volt of input (ratio of leading coefficients).
> - Poles of $H$: $s^2+3s+50=0\Rightarrow s=-1.5\pm j6.91$, in the LHP, so $y(\infty)=H(0)=\dfrac{50}{50}=1$.
>
> **General rule:** for a stable system, the final value of the step response is the **DC gain** $H(0)$. The step response here overshoots substantially (damping ratio $\zeta=1.5/\sqrt{50}=0.21$), something the FVT says nothing about — it describes only where the response ends, not how it gets there.

> [!warning] Common Pitfall — using the theorems on $X$ instead of $sX$
> Both theorems involve $sX(s)$, not $X(s)$. Forgetting the factor $s$ gives $\lim_{s\to0}X(s)$, which for a step response is infinite and for an impulse response equals the *area* $\int_0^\infty x\,dt$ — a different (and sometimes useful) quantity.

> [!tip] Engineering Note — steady-state error in control loops
> In a unity-feedback loop with open-loop transfer function $G(s)$, the error for a reference $R(s)$ is $E(s)=R(s)/(1+G(s))$. The FVT gives the steady-state error $e_{ss}=\lim_{s\to0}sR/(1+G)$. For a unit step this is $1/(1+G(0))$, which is zero exactly when $G$ has a pole at $s=0$ (an integrator in the loop). This one-line calculation is the basis of the "system type" classification in control engineering (Section 14.9).

## 14.8 Solving differential equations

This is the application that justifies the whole apparatus. A linear constant-coefficient ODE with initial conditions — an **initial-value problem (IVP)** — becomes an algebraic equation for $Y(s)$. Solving it produces, in a single expression, both the natural response (fixed by the initial conditions) and the forced response (fixed by the input). Discontinuous and impulsive inputs, which are awkward for the classical method of [Ch. 12](../part-4-differential-equations/12-linear-differential-equations.md), cause no extra difficulty.

### The general procedure

Consider the $n$th-order equation

$$
a_ny^{(n)}+a_{n-1}y^{(n-1)}+\cdots+a_1\dot y+a_0y=x(t),\qquad y(0^-),\ \dot y(0^-),\ \ldots,\ y^{(n-1)}(0^-)\ \text{given}.
$$

1. **Transform** both sides, using the derivative rule for each $y^{(k)}$ and the table for $x(t)$.
2. **Collect** the terms in $Y(s)$ on the left. The coefficient of $Y$ is the **characteristic polynomial** $P(s)=a_ns^n+\cdots+a_0$ — the same polynomial whose roots gave the exponents in the classical method.
3. **Solve** algebraically for $Y(s)$, which always takes the form shown below.
4. **Invert** by partial fractions (Section 14.5) and the table.
5. **Check** the initial conditions (IVT) and, if applicable, the final value (FVT) — or substitute back.

$$
Y(s)=\underbrace{\frac{C(s)}{P(s)}}_{\text{zero-input}}+\underbrace{\frac{X(s)}{P(s)}}_{\text{zero-state}},
$$

where $C(s)$ is a polynomial of degree $<n$ built from the initial conditions.

The **zero-input response** is what the system does from its initial state with no input; the **zero-state response** is what the input does to a system initially at rest. Their sum is the complete response. The ratio $H(s)=1/P(s)$ that multiplies $X(s)$ is the transfer function of the equation.

> [!theorem] Structure of the solution
> The poles of $Y(s)$ are the roots of $P(s)$ (the **natural frequencies** or characteristic roots) together with the poles of $X(s)$ (the **forced frequencies**). The terms belonging to the roots of $P$ form the natural (transient) response; those belonging to the poles of $X$ form the forced (particular) response. If a pole of $X$ coincides with a root of $P$, the pole of $Y$ is repeated and terms $t^ke^{pt}$ appear — the resonance case that the classical method treats by "multiplying the trial solution by $t$".

> [!example] Worked Example 14.32 — first-order RC charging
> A capacitor voltage obeys $\dot y+2y=4$, $y(0^-)=1$ (time in ms, voltage in V). Find $y(t)$.
>
> **Transform.** $sY-1+2Y=\dfrac4s$.
>
> **Solve.** $Y=\dfrac{1}{s+2}+\dfrac{4}{s(s+2)}=\dfrac{s+4}{s(s+2)}$.
>
> **Partial fractions** (cover-up): at $s=0$: $4/2=2$; at $s=-2$: $2/(-2)=-1$.
>
> $$
> Y=\frac2s-\frac{1}{s+2}\quad\Longrightarrow\quad y(t)=2-e^{-2t}\ \text{V}.
> $$
>
> **Check.** $y(0^+)=1$ ✓ (no jump, as expected for a capacitor voltage); $\dot y+2y=2e^{-2t}+4-2e^{-2t}=4$ ✓; final value 2 V = forcing/coefficient ✓.

> [!example] Worked Example 14.33 — second order, zero-input and zero-state parts
> Solve $\ddot y+3\dot y+2y=e^{-3t}$ with $y(0^-)=1$, $\dot y(0^-)=0$.
>
> **Transform.**
>
> $$
> \big[s^2Y-s\cdot1-0\big]+3\big[sY-1\big]+2Y=\frac{1}{s+3}.
> $$
>
> **Collect.** $(s^2+3s+2)Y=(s+3)+\dfrac{1}{s+3}$, with $P(s)=(s+1)(s+2)$:
>
> $$
> Y=\underbrace{\frac{s+3}{(s+1)(s+2)}}_{Y_{zi}}+\underbrace{\frac{1}{(s+1)(s+2)(s+3)}}_{Y_{zs}}.
> $$
>
> **Zero-input.** Cover-up: $\dfrac{2}{s+1}-\dfrac{1}{s+2}\leftrightarrow2e^{-t}-e^{-2t}$.
>
> **Zero-state.** Cover-up: at $-1$: $\frac{1}{(1)(2)}=\frac12$; at $-2$: $\frac{1}{(-1)(1)}=-1$; at $-3$: $\frac{1}{(-2)(-1)}=\frac12$. So $\tfrac12e^{-t}-e^{-2t}+\tfrac12e^{-3t}$.
>
> **Complete response.**
>
> $$
> y(t)=\tfrac52e^{-t}-2e^{-2t}+\tfrac12e^{-3t}.
> $$
>
> **Check.** $y(0)=\tfrac52-2+\tfrac12=1$ ✓; $\dot y(0)=-\tfrac52+4-\tfrac32=0$ ✓. Substituting into the ODE: the $e^{-t}$ and $e^{-2t}$ terms are annihilated by $P$ (since $P(-1)=P(-2)=0$), and the $e^{-3t}$ term gives $\tfrac12P(-3)e^{-3t}=\tfrac12\cdot2\,e^{-3t}=e^{-3t}$ ✓.
>
> Note how the classical method's three steps (complementary function, particular integral, fitting constants) collapsed into one partial-fraction expansion.

### Discontinuous and impulsive inputs

> [!example] Worked Example 14.34 — a rectangular pulse into an undamped oscillator
> An undamped LC tank, initially at rest, obeys $\ddot y+4y=x(t)$ where $x(t)=1$ for $0\le t<\pi$ and 0 afterwards. Find $y(t)$.
>
> **Input.** $x=u(t)-u(t-\pi)\leftrightarrow\dfrac{1-e^{-\pi s}}{s}$.
>
> **Transform.** $(s^2+4)Y=\dfrac{1-e^{-\pi s}}{s}$, so $Y=(1-e^{-\pi s})\,G(s)$ with
>
> $$
> G(s)=\frac{1}{s(s^2+4)}=\frac14\left(\frac1s-\frac{s}{s^2+4}\right)\leftrightarrow g(t)=\frac{1-\cos2t}{4}.
> $$
>
> **Invert with the delay.**
>
> $$
> y(t)=\frac{1-\cos2t}{4}-\frac{1-\cos2(t-\pi)}{4}\,u(t-\pi).
> $$
>
> For $t\ge\pi$, $\cos2(t-\pi)=\cos2t$, so the two terms cancel exactly: **$y(t)=0$ for $t\ge\pi$.**
>
> **Interpretation.** The pulse lasts exactly one natural period ($2\pi/2=\pi$). The energy it pumps in during the first half-period is withdrawn in the second, leaving the oscillator at rest: $y(\pi)=0$ and $\dot y(\pi)=\tfrac12\sin2\pi=0$. A pulse of any other length would leave a residual oscillation. This "input shaping" idea is used to move resonant mechanical loads (disk-drive heads, cranes) without exciting vibration.

```desmos-graph
left=-0.3; right=7; bottom=-0.2; top=0.8
---
y=\frac{1-\cos(2x)}{4}|blue|0<x<3.14159
y=0|blue|3.14159<x<7
y=\frac{1-\cos(2x)}{4}|grey|dashed|3.14159<x<7
y=0.7|orange|0<x<3.14159
```
*Example 14.34: during the pulse (orange, scaled) the oscillator swings once (blue) and is left exactly at rest at $t=\pi$; without the pulse ending, it would continue oscillating (grey dashed).*

> [!example] Worked Example 14.35 — resonance
> Solve $\ddot y+\omega^2y=\sin\omega t$, $y(0^-)=\dot y(0^-)=0$.
>
> $$
> Y=\frac{\omega}{(s^2+\omega^2)^2}=\frac{1}{2\omega^2}\cdot\frac{2\omega^3}{(s^2+\omega^2)^2}\quad\Longrightarrow\quad y(t)=\frac{\sin\omega t-\omega t\cos\omega t}{2\omega^2}
> $$
>
> by table entry 29. The forcing frequency coincides with the natural frequency, the poles $\pm j\omega$ are doubled, and the amplitude grows linearly as $t/(2\omega)$ without bound. In a real circuit the resistance limits the growth; the ideal model warns that a lightly damped resonator driven at resonance develops very large voltages (a factor $Q$ above the drive).
>
> **Check** ($\omega=2$, $t=1$): $y=\frac{\sin2-2\cos2}{8}=\frac{0.9093+0.8323}{8}=0.2177$. Numerical integration of the ODE gives $0.2177$. ✓

### Integro-differential equations

Equations containing both derivatives and running integrals of the unknown — typical of loop equations with capacitors — are handled just as easily, because the integral rule turns $\int_0^t y\,d\tau$ into $Y/s$.

> [!example] Worked Example 14.36 — an integro-differential equation
> Solve $\dot y+3y+2\displaystyle\int_0^ty(\tau)\,d\tau=u(t)$, $y(0^-)=0$.
>
> **Transform.** $sY+3Y+\dfrac{2Y}{s}=\dfrac1s$. Multiply by $s$: $(s^2+3s+2)Y=1$.
>
> $$
> Y=\frac{1}{(s+1)(s+2)}=\frac{1}{s+1}-\frac{1}{s+2}\quad\Longrightarrow\quad y(t)=e^{-t}-e^{-2t}.
> $$
>
> This is the loop current of a series $RLC$ circuit ($L=1$, $R=3$, $1/C=2$) connected to a unit step: it rises from zero, peaks at $t=\ln2=0.693$ (value $0.25$), and decays because the charged capacitor eventually blocks DC.

### Simultaneous equations

Coupled circuits (mutual inductance, multi-loop networks) and state-space models give systems of ODEs. Each equation transforms to a linear algebraic equation, and the resulting linear system in $X(s),Y(s),\ldots$ is solved by elimination or Cramer's rule ([Ch. 8](../part-3-linear-algebra/08-matrices-and-linear-systems.md)).

In matrix form, $\dot{\mathbf{x}}=\mathbf{A}\mathbf{x}+\mathbf{b}\,v(t)$ transforms to $s\mathbf{X}-\mathbf{x}(0^-)=\mathbf{A}\mathbf{X}+\mathbf{b}V$, so

$$
\mathbf{X}(s)=(s\mathbf{I}-\mathbf{A})^{-1}\big[\mathbf{x}(0^-)+\mathbf{b}V(s)\big].
$$

Since $(s\mathbf{I}-\mathbf{A})^{-1}=\operatorname{adj}(s\mathbf{I}-\mathbf{A})/\det(s\mathbf{I}-\mathbf{A})$, the poles of every solution component are the **eigenvalues of $\mathbf{A}$** ([Ch. 10](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md)), and the matrix exponential is $e^{\mathbf{A}t}=\mathcal{L}^{-1}\{(s\mathbf{I}-\mathbf{A})^{-1}\}$.

> [!example] Worked Example 14.37 — a coupled pair
> Solve $\dot x=-2x+y$, $\dot y=x-2y+2$ for $t\ge0$, with $x(0^-)=1$, $y(0^-)=0$.
>
> **Transform.**
>
> $$
> \begin{aligned}
> (s+2)X-Y&=1,\\
> -X+(s+2)Y&=\frac2s .
> \end{aligned}
> $$
>
> **Cramer's rule.** $\Delta=(s+2)^2-1=(s+1)(s+3)$ — the eigenvalues of $\mathbf{A}=\begin{bmatrix}-2&1\\1&-2\end{bmatrix}$ are $-1$ and $-3$.
>
> $$
> X=\frac{(s+2)+\frac2s}{(s+1)(s+3)}=\frac{s^2+2s+2}{s(s+1)(s+3)},\qquad
> Y=\frac{(s+2)\frac2s+1}{(s+1)(s+3)}=\frac{3s+4}{s(s+1)(s+3)} .
> $$
>
> **Partial fractions (cover-up).**
>
> $$
> X=\frac{2/3}{s}-\frac{1/2}{s+1}+\frac{5/6}{s+3},\qquad Y=\frac{4/3}{s}-\frac{1/2}{s+1}-\frac{5/6}{s+3}.
> $$
>
> $$
> x(t)=\tfrac23-\tfrac12e^{-t}+\tfrac56e^{-3t},\qquad y(t)=\tfrac43-\tfrac12e^{-t}-\tfrac56e^{-3t}.
> $$
>
> **Check.** $x(0)=\tfrac23-\tfrac12+\tfrac56=1$ ✓, $y(0)=\tfrac43-\tfrac12-\tfrac56=0$ ✓. Steady state: setting derivatives to zero, $-2x+y=0$ and $x-2y+2=0$ give $x=\tfrac23$, $y=\tfrac43$ ✓.

### Equations with polynomial coefficients

The transform is designed for constant coefficients, but the rule $\mathcal{L}\{t\,y\}=-dY/ds$ converts a coefficient $t$ into a derivative in $s$. An ODE in $t$ with linear coefficients becomes a *first-order* ODE in $s$, which is sometimes easier.

> [!example] Worked Example 14.38 — Bessel's equation of order zero
> Solve $t\ddot y+\dot y+ty=0$ with $y(0)=1$ (bounded solution).
>
> **Transform each term.** $\mathcal{L}\{t\ddot y\}=-\frac{d}{ds}\big[s^2Y-s\,y(0)-\dot y(0)\big]=-2sY-s^2Y'+1$; $\mathcal{L}\{\dot y\}=sY-1$; $\mathcal{L}\{ty\}=-Y'$.
>
> **Sum.** $-(s^2+1)Y'-sY=0$, a separable first-order ODE in $s$:
>
> $$
> \frac{Y'}{Y}=-\frac{s}{s^2+1}\quad\Longrightarrow\quad Y=\frac{C}{\sqrt{s^2+1}} .
> $$
>
> **Constant.** IVT: $y(0^+)=\lim sY=C=1$. So $Y=1/\sqrt{s^2+1}$, i.e. $y(t)=J_0(t)$ (table entry 38), the Bessel function that describes, among other things, the current distribution in a round conductor at high frequency.

> [!warning] Common Pitfall — initial conditions at the wrong instant
> If a problem states "the switch closes at $t=0$ and $v_C(0^+)=\ldots$", the $0^-$ derivative rule cannot use that value directly. Either compute $v_C(0^-)$ from the circuit before switching (the usual case — capacitor voltages and inductor currents are continuous unless impulses act), or switch consistently to the $0^+$ convention and exclude impulses at the origin from the input. Mixing the two produces answers that are wrong by an impulse or a step.

> [!tip] Engineering Note — why simulation tools think in $s$
> SPICE's ".ac" analysis evaluates the zero-state transfer function on the $j\omega$-axis, and its ".tran" analysis integrates the ODE numerically. Laplace-domain hand analysis sits between the two: it gives the transient in closed form, exposes which poles cause ringing or slow settling, and provides the reference solution against which a numerical integrator is checked ([Ch. 22](../part-7-numerical-computational/22-numerical-methods.md)).

## 14.9 Applications to circuits and control

The previous section transformed a differential equation that had already been written down. For circuits there is a shortcut: transform the **element laws** once and for all, draw an $s$-domain circuit, and apply Kirchhoff's laws and ordinary resistive-circuit techniques (series/parallel combination, voltage dividers, nodal and mesh analysis, Thévenin equivalents) directly to it. No differential equation is ever written. This section develops that method, and then the language it leads to: transfer functions, poles and zeros, stability and feedback. As elsewhere in this book, the circuits are there to motivate the mathematics; circuit theory proper belongs to its own course.

### Kirchhoff's laws in the $s$-domain

KCL ($\sum i_k(t)=0$ at a node) and KVL ($\sum v_k(t)=0$ around a loop) are linear relations with constant coefficients, so by linearity of the transform they hold for the transforms too: $\sum I_k(s)=0$ and $\sum V_k(s)=0$. Only the element laws need translating.

### $s$-domain element models

**Resistor.** $v=Ri$ transforms to $V(s)=R\,I(s)$.

**Inductor.** $v=L\,\dfrac{di}{dt}$ transforms, by the derivative rule, to

$$
V(s)=sL\,I(s)-L\,i(0^-).
$$

This is the law of an impedance $sL$ **in series** with a voltage source $L\,i(0^-)$ (polarity opposing the stored current, as the minus sign shows). Solving instead for $I$:

$$
I(s)=\frac{V(s)}{sL}+\frac{i(0^-)}{s},
$$

an impedance $sL$ **in parallel** with a current source $i(0^-)/s$ — the transform of a constant current $i(0^-)$.

**Capacitor.** $i=C\,\dfrac{dv}{dt}$ transforms to

$$
I(s)=sC\,V(s)-C\,v(0^-)\qquad\Longleftrightarrow\qquad V(s)=\frac{I(s)}{sC}+\frac{v(0^-)}{s}.
$$

The first form is an admittance $sC$ in parallel with a current source $C\,v(0^-)$ (the transform of the impulse $Cv(0^-)\delta(t)$, i.e. the stored charge); the second is an impedance $1/(sC)$ in series with a voltage source $v(0^-)/s$ (a step of the stored voltage).

> [!definition] Impedance and admittance
> For an element or two-terminal network initially at rest, the **impedance** is $Z(s)=V(s)/I(s)$ and the **admittance** $Y(s)=1/Z(s)$:
>
> $$
> Z_R=R,\qquad Z_L=sL,\qquad Z_C=\frac{1}{sC}.
> $$
>
> Impedances combine like resistances: in series they add, in parallel their admittances add.

On the $j\omega$-axis these become the phasor impedances $j\omega L$ and $1/(j\omega C)$ of AC circuit analysis. The $s$-domain impedance is their generalisation to arbitrary (including transient and growing) signals, and it carries the initial energy through the extra sources.

**Units check.** $sL$ has units $\mathrm{s^{-1}\cdot H}=\Omega$ ✓; $1/(sC)$ has $1/(\mathrm{s^{-1}\cdot F})=\Omega$ ✓; $L\,i(0^-)$ has $\mathrm{H\cdot A}=\mathrm{V\cdot s}$, the units of a transformed voltage ✓.

| Element | Time law | $s$-domain, series form | $s$-domain, parallel form |
|---|---|---|---|
| $R$ | $v=Ri$ | $V=RI$ | $I=V/R$ |
| $L$ | $v=L\,di/dt$ | $V=sLI-Li(0^-)$ | $I=\dfrac{V}{sL}+\dfrac{i(0^-)}{s}$ |
| $C$ | $i=C\,dv/dt$ | $V=\dfrac{I}{sC}+\dfrac{v(0^-)}{s}$ | $I=sCV-Cv(0^-)$ |

> [!example] Worked Example 14.39 — an inductor with stored current
> A series circuit contains a source $v_s(t)=20u(t)$ V, a resistor $R=5\ \Omega$ and an inductor $L=0.5$ H carrying $i(0^-)=2$ A (from an earlier source configuration). Find $i(t)$.
>
> **$s$-domain circuit.** Source $20/s$; impedances $5$ and $0.5s$; the inductor's series source $L\,i(0^-)=1$ V·s, aiding the existing current direction.
>
> **KVL.** $\dfrac{20}{s}+1=(5+0.5s)\,I(s)$, so
>
> $$
> I(s)=\frac{20/s+1}{0.5(s+10)}=\frac{40+2s}{s(s+10)}=\frac{4}{s}-\frac{2}{s+10}.
> $$
>
> (Cover-up: at $0$, $40/10=4$; at $-10$, $(40-20)/(-10)=-2$.)
>
> $$
> i(t)=4-2e^{-10t}\ \text{A},\qquad t\ge0.
> $$
>
> **Check.** $i(0^+)=2$ A — continuous, as an inductor current must be ✓. $i(\infty)=20/5=4$ A, the DC value with the inductor a short circuit ✓. Time constant $L/R=0.1$ s, matching the pole at $-10\ \mathrm{s^{-1}}$ ✓.

> [!example] Worked Example 14.40 — series RLC with a charged capacitor
> A series $RLC$ circuit has $R=12\ \Omega$, $L=10$ mH, $C=100\ \mu$F. The capacitor holds $v_C(0^-)=2$ V and $i(0^-)=0$. At $t=0$ a 10 V DC source is switched in. Find $i(t)$ and $v_C(t)$.
>
> **$s$-domain circuit.** Source $10/s$; impedances $R$, $sL$, $1/(sC)$; the capacitor's series source $v_C(0^-)/s=2/s$ opposes the drive (it is charged with the same polarity as the source would charge it).
>
> **KVL.**
>
> $$
> \frac{10}{s}-\frac{2}{s}=\left(R+sL+\frac{1}{sC}\right)I(s)
> \quad\Longrightarrow\quad
> I(s)=\frac{8/s}{R+sL+\frac{1}{sC}}=\frac{8}{L\left(s^2+\frac RLs+\frac1{LC}\right)} .
> $$
>
> With $R/L=1200\ \mathrm{s^{-1}}$, $1/(LC)=10^6\ \mathrm{s^{-2}}$, $1/L=100\ \mathrm{H^{-1}}$:
>
> $$
> I(s)=\frac{800}{s^2+1200s+10^6}=\frac{800}{(s+600)^2+800^2}.
> $$
>
> **Inverse** (entry 23):
>
> $$
> i(t)=e^{-600t}\sin(800t)\ \text{A}.
> $$
>
> **Capacitor voltage.** $V_C=\dfrac{I}{sC}+\dfrac2s=\dfrac{8\times10^6}{s(s^2+1200s+10^6)}+\dfrac2s$. Partial fractions on the first term give $\dfrac8s-\dfrac{8(s+600)+6\cdot800}{(s+600)^2+800^2}$, so
>
> $$
> v_C(t)=10-e^{-600t}\left(8\cos800t+6\sin800t\right)\ \text{V}.
> $$
>
> **Checks.** $v_C(0^+)=10-8=2$ V ✓ (continuous). $v_C(\infty)=10$ V ✓ (FVT). $i(0^+)=0$ ✓. Also $C\,dv_C/dt=10^{-4}\cdot e^{-600t}\big[600(8\cos+6\sin)-800(-8\sin+6\cos)\big]=10^{-4}e^{-600t}(10000\sin800t)=i$ ✓.

### Transfer functions

> [!definition] Transfer function
> For a linear time-invariant system **initially at rest**, with input $x(t)$ and output $y(t)$, the **transfer function** is
>
> $$
> H(s)=\frac{Y(s)}{X(s)}=\mathcal{L}\{h(t)\},
> $$
>
> where $h(t)$ is the impulse response. For a lumped system, $H(s)$ is a real rational function
>
> $$
> H(s)=\frac{N(s)}{D(s)}=K\,\frac{(s-z_1)(s-z_2)\cdots(s-z_m)}{(s-p_1)(s-p_2)\cdots(s-p_n)}.
> $$
>
> The roots $z_i$ of $N$ are the **zeros**, the roots $p_k$ of $D$ the **poles**, and $K$ the gain constant. $H$ is **proper** if $m\le n$.

For the ODE $a_ny^{(n)}+\cdots+a_0y=b_mx^{(m)}+\cdots+b_0x$, transforming with zero initial conditions gives $H(s)=\dfrac{b_ms^m+\cdots+b_0}{a_ns^n+\cdots+a_0}$: the transfer function can be written down by inspection, and conversely. For a circuit it is found from the $s$-domain circuit with all initial-condition sources removed, using dividers and impedance combinations.

The **pole–zero plot** marks poles with × and zeros with ○ in the $s$-plane. Together with $K$, it determines $H$ completely. Its interpretation:

- **Poles** are the natural frequencies. Each contributes a mode $e^{p_kt}$ (or $t^je^{p_kt}$ if repeated) to the impulse response and to every transient.
- **Zeros** are frequencies at which the system blocks transmission: an input $e^{z_it}$ produces no forced output. They do not create modes; they set the *amplitudes* (residues) of the modes and shape the frequency response.
- **Frequency response.** For a stable system, the steady-state response to $\cos\omega t$ is $\lvert H(j\omega)\rvert\cos(\omega t+\angle H(j\omega))$. Geometrically, $\lvert H(j\omega)\rvert=\lvert K\rvert\frac{\prod\lvert j\omega-z_i\rvert}{\prod\lvert j\omega-p_k\rvert}$: the product of distances from the point $j\omega$ to the zeros divided by that to the poles. A pole close to the $j\omega$-axis produces a sharp resonant peak; a zero on the axis produces a notch. Bode plots and the Fourier view of $H(j\omega)$ are developed in [Ch. 16](16-fourier-transforms.md).

> [!example] Worked Example 14.41 — zeros on the $j\omega$-axis: a notch filter
> A notch filter for removing 50 Hz mains interference has
>
> $$
> H(s)=\frac{s^2+\omega_0^2}{s^2+\frac{\omega_0}{Q}s+\omega_0^2},\qquad \omega_0=2\pi\cdot50=314.16\ \text{rad/s},\quad Q=5 .
> $$
>
> Find the poles and zeros, and $\lvert H(j\omega)\rvert$ at $\omega=0$, $\omega_0$, $1.1\omega_0$ and $\infty$.
>
> **Zeros.** $s^2+\omega_0^2=0\Rightarrow s=\pm j\omega_0=\pm j314.16$: on the $j\omega$-axis.
>
> **Poles.** $s=-\frac{\omega_0}{2Q}\pm j\omega_0\sqrt{1-\frac{1}{4Q^2}}=-31.42\pm j312.58$: in the LHP, so the filter is stable, and close to the zeros.
>
> **Frequency response.** On the axis, $H(j\omega)=\dfrac{\omega_0^2-\omega^2}{\omega_0^2-\omega^2+j\omega\omega_0/Q}$.
>
> - $\omega=0$: $\lvert H\rvert=1$ (DC passes).
> - $\omega=\omega_0$: numerator zero, $\lvert H\rvert=0$ — the 50 Hz component is removed completely.
> - $\omega=1.1\omega_0$: $\omega_0^2-\omega^2=-0.21\omega_0^2$, $\omega\omega_0/Q=0.22\omega_0^2$, so $\lvert H\rvert=\frac{0.21}{\sqrt{0.21^2+0.22^2}}=0.690$.
> - $\omega\to\infty$: $\lvert H\rvert\to1$.
>
> **Geometric reading.** At $s=j\omega$ near $j\omega_0$, the distance to the zero is tiny while the distance to the nearby pole is about $\omega_0/(2Q)$, so the ratio — and $\lvert H\rvert$ — dips sharply. Raising $Q$ moves the poles closer to the axis and narrows the notch: the $-3$ dB width is $\omega_0/Q=62.8$ rad/s, i.e. 10 Hz.

### First- and second-order standard forms

Most systems are understood through their dominant first- and second-order factors.

**First order.** $H(s)=\dfrac{K}{1+s\tau}$: one real pole at $-1/\tau$. The unit-step response is $K(1-e^{-t/\tau})$: it reaches 63.2% of its final value at $t=\tau$ and 98.2% at $4\tau$. The $RC$ low-pass filter has $\tau=RC$, $K=1$.

**Second order.**

$$
H(s)=\frac{\omega_n^2}{s^2+2\zeta\omega_ns+\omega_n^2},
$$

with **natural frequency** $\omega_n>0$ and **damping ratio** $\zeta\ge0$. The poles are

$$
p_{1,2}=-\zeta\omega_n\pm\omega_n\sqrt{\zeta^2-1}.
$$

| Damping | Poles | Step response |
|---|---|---|
| $\zeta>1$ overdamped | two distinct real, negative | sum of two decaying exponentials, no overshoot |
| $\zeta=1$ critically damped | double real at $-\omega_n$ | $1-(1+\omega_nt)e^{-\omega_nt}$, fastest without overshoot |
| $0<\zeta<1$ underdamped | $-\zeta\omega_n\pm j\omega_d$, $\omega_d=\omega_n\sqrt{1-\zeta^2}$ | decaying oscillation with overshoot |
| $\zeta=0$ undamped | $\pm j\omega_n$ | $1-\cos\omega_nt$, sustained oscillation |

For $0<\zeta<1$ the poles lie on a circle of radius $\omega_n$, at angle $\theta=\arccos\zeta$ from the negative real axis. Partial fractions on $H(s)/s$ give the underdamped step response

$$
y(t)=1-e^{-\zeta\omega_nt}\left(\cos\omega_dt+\frac{\zeta}{\sqrt{1-\zeta^2}}\sin\omega_dt\right)
=1-\frac{e^{-\zeta\omega_nt}}{\sqrt{1-\zeta^2}}\sin\!\left(\omega_dt+\arccos\zeta\right).
$$

Setting $\dot y=0$ (equivalently, the impulse response $h(t)=\frac{\omega_n}{\sqrt{1-\zeta^2}}e^{-\zeta\omega_nt}\sin\omega_dt=0$) locates the first peak at $t_p=\pi/\omega_d$, where

$$
\text{overshoot}=y(t_p)-1=e^{-\pi\zeta/\sqrt{1-\zeta^2}} .
$$

The envelope $e^{-\zeta\omega_nt}$ gives the 2% **settling time** $t_s\approx4/(\zeta\omega_n)$: the real part of the poles sets the decay rate, the imaginary part the ringing frequency.

> [!example] Worked Example 14.42 — the RLC low-pass filter
> Take the series $RLC$ circuit of Example 14.40 ($R=12\ \Omega$, $L=10$ mH, $C=100\ \mu$F) with output across the capacitor. Find $H(s)$, its poles and zeros, $\omega_n$, $\zeta$, and the step-response overshoot and peak time.
>
> **Voltage divider** (zero initial conditions):
>
> $$
> H(s)=\frac{V_C}{V_s}=\frac{1/(sC)}{R+sL+1/(sC)}=\frac{1}{LCs^2+RCs+1}=\frac{1/(LC)}{s^2+\frac RLs+\frac1{LC}}=\frac{10^6}{s^2+1200s+10^6}.
> $$
>
> **Poles and zeros.** No finite zeros; poles $-600\pm j800$ (as in Example 14.40 — the same circuit has the same natural frequencies whatever the output).
>
> **Parameters.** $\omega_n=\sqrt{10^6}=1000$ rad/s; $2\zeta\omega_n=1200\Rightarrow\zeta=0.6$; $\omega_d=800$ rad/s. In circuit terms $\zeta=\frac R2\sqrt{C/L}=6\sqrt{0.01}=0.6$ ✓.
>
> **Step response.** $y(t)=1-e^{-600t}\left(\cos800t+0.75\sin800t\right)$.
>
> - Peak time $t_p=\pi/800=3.93$ ms.
> - Overshoot $e^{-\pi(0.6)/0.8}=e^{-2.356}=0.0948$, i.e. 9.5%: the peak is 1.095.
> - Settling time $t_s\approx4/600=6.7$ ms.
>
> **Check.** Evaluating $y(t)$ on a fine grid gives maximum $1.0948$ at $t=3.927$ ms ✓. DC gain $H(0)=1$ ✓ (at DC the capacitor is open, so the full input appears across it).

```desmos-graph
left=0; right=14; bottom=0; top=1.6
---
y=1-e^{-0.2x}\left(\cos(0.9798x)+0.2041\sin(0.9798x)\right)|red|0<x<14
y=1-e^{-0.6x}\left(\cos(0.8x)+0.75\sin(0.8x)\right)|blue|0<x<14
y=1-(1+x)e^{-x}|green|0<x<14
y=1|grey|dashed
```
*Unit-step responses of the standard second-order system with $\omega_n=1$ (time in units of $1/\omega_n$): $\zeta=0.2$ (red, 53% overshoot), $\zeta=0.6$ (blue, 9.5% overshoot — the RLC filter of Example 14.42 with time in ms), $\zeta=1$ (green, critically damped).*

### Stability

> [!definition] BIBO stability
> A system is **bounded-input, bounded-output (BIBO) stable** if every bounded input $\lvert x(t)\rvert\le M$ produces a bounded output. For an LTI system this holds if and only if the impulse response is absolutely integrable: $\int_0^\infty\lvert h(t)\rvert\,dt<\infty$.

*Why.* If $\int\lvert h\rvert=A<\infty$ then $\lvert y(t)\rvert\le\int_0^t\lvert h(\tau)\rvert\lvert x(t-\tau)\rvert d\tau\le MA$. Conversely, if $\int\lvert h\rvert=\infty$, the bounded input $x(t-\tau)=\operatorname{sgn}h(\tau)$ makes $y(t)=\int_0^t\lvert h\rvert d\tau$ grow without bound.

> [!theorem] Pole criterion for stability
> A system with proper rational transfer function $H(s)$ (in lowest terms) is BIBO stable **if and only if every pole has negative real part** — all poles lie in the open left half-plane.
>
> - Any pole with $\operatorname{Re}\{p\}>0$, or a repeated pole on the $j\omega$-axis, gives an unbounded impulse response: **unstable**.
> - Simple poles on the $j\omega$-axis (and none in the RHP) give a bounded but non-decaying impulse response: **marginally stable**, which is *not* BIBO stable (a bounded input at the pole frequency produces resonant growth, Example 14.35).

*Proof sketch.* By partial fractions, $h(t)=K_0\delta(t)+\sum_k\sum_iK_{k,i}\frac{t^{i-1}}{(i-1)!}e^{p_kt}$. Each term with $\operatorname{Re}\{p_k\}<0$ is absolutely integrable (exponential decay beats any polynomial). A term with $\operatorname{Re}\{p_k\}\ge0$ and $K\ne0$ is not. Since the fraction is in lowest terms, no pole has zero residue at its highest power. $\square$

Finding all roots of a high-order $D(s)$ is unnecessary to decide stability. The **Routh–Hurwitz criterion** answers it from the coefficients alone. Two consequences are easy to remember:

- **Necessary condition:** if all roots of $a_ns^n+\cdots+a_0$ (with $a_n>0$) have negative real parts, then every coefficient is **positive**. A missing or negative coefficient means instability (or a $j\omega$-axis root).
- For **second order** positivity of all coefficients is also sufficient. For **third order**, $a_3s^3+a_2s^2+a_1s+a_0$ with all $a_i>0$ is stable if and only if $a_2a_1>a_3a_0$.

> [!example] Worked Example 14.43 — a stability range
> For which gains $K$ is the system with characteristic polynomial $D(s)=s^3+2s^2+3s+K$ stable? What happens at the boundary?
>
> **Positivity.** Need $K>0$.
>
> **Third-order condition.** $a_2a_1>a_3a_0$: $2\cdot3>1\cdot K$, so $K<6$. Stable for $0<K<6$.
>
> **Boundary.** At $K=6$: $D(s)=s^3+2s^2+3s+6=(s+2)(s^2+3)$, poles $-2$ and $\pm j\sqrt3=\pm j1.732$. The system oscillates at $\sqrt3$ rad/s with constant amplitude — marginal stability. For $K>6$ the complex pair crosses into the right half-plane. This is precisely how an oscillator is designed: raise the loop gain until a pole pair reaches the $j\omega$-axis.

### Feedback and control

In a feedback loop with forward path $G(s)$ and feedback path $F(s)$, the error signal is $E=R-FY$ and $Y=GE$. Eliminating $E$:

$$
T(s)=\frac{Y(s)}{R(s)}=\frac{G(s)}{1+G(s)F(s)} .
$$

The closed-loop poles are the roots of the **characteristic equation** $1+G(s)F(s)=0$ — they are not the poles of $G$. Feedback therefore *moves poles*, which is both its power (stabilising, speeding up, desensitising) and its danger (it can destabilise a stable plant, as in Example 14.43). For **unity feedback** ($F=1$) the error is $E=R/(1+G)$, and the final-value theorem gives the steady-state errors:

| Reference | $e_{ss}$ | Error constant |
|---|---|---|
| Step $1/s$ | $\dfrac{1}{1+K_p}$ | $K_p=\lim_{s\to0}G(s)$ |
| Ramp $1/s^2$ | $\dfrac{1}{K_v}$ | $K_v=\lim_{s\to0}sG(s)$ |
| Parabola $1/s^3$ | $\dfrac{1}{K_a}$ | $K_a=\lim_{s\to0}s^2G(s)$ |

The number of pure integrators ($1/s$ factors) in $G$ is the **system type**: a type-1 loop tracks a step with zero error, a type-2 loop tracks a ramp with zero error.

> [!example] Worked Example 14.44 — a position servo
> A motor-driven positioning system has open-loop transfer function $G(s)=\dfrac{K}{s(s+4)}$ in a unity-feedback loop. Choose $K$ for a damping ratio $\zeta=0.5$; then find the overshoot, peak time, and steady-state errors for step and ramp references.
>
> **Closed loop.**
>
> $$
> T(s)=\frac{K/(s(s+4))}{1+K/(s(s+4))}=\frac{K}{s^2+4s+K}.
> $$
>
> Compare with the standard form: $\omega_n^2=K$, $2\zeta\omega_n=4$. With $\zeta=0.5$: $\omega_n=4$ rad/s, so $K=16$.
>
> **Poles.** $s^2+4s+16=0\Rightarrow s=-2\pm j2\sqrt3=-2\pm j3.464$.
>
> **Transient.** Overshoot $e^{-\pi(0.5)/\sqrt{0.75}}=e^{-1.814}=0.163$ (16.3%); peak time $t_p=\pi/\omega_d=\pi/3.464=0.907$ s; settling time $\approx4/2=2$ s. Partial fractions give
>
> $$
> y(t)=1-\frac{2}{\sqrt3}\,e^{-2t}\sin\!\left(2\sqrt3\,t+\frac{\pi}{3}\right).
> $$
>
> A fine-grid evaluation gives a maximum of $1.1630$ at $t=0.907$ s ✓.
>
> **Steady-state errors.** $G$ has one integrator: type 1. Step: $K_p=\infty$, $e_{ss}=0$. Ramp: $K_v=\lim sG=K/4=4\ \mathrm{s^{-1}}$, so $e_{ss}=1/4$ — the output lags a unit-slope ramp by 0.25 units.
>
> **Trade-off.** Raising $K$ reduces the ramp error ($e_{ss}=4/K$) but lowers $\zeta=2/\sqrt K$ and increases overshoot: the basic compromise of control design, visible entirely through pole locations.

> [!warning] Common Pitfall — pole–zero cancellation hides instability
> Mathematically $\dfrac{s-1}{(s-1)(s+2)}=\dfrac{1}{s+2}$, which looks stable. Physically, if the factor $(s-1)$ belongs to the plant and is "cancelled" by a controller zero, the unstable mode $e^{t}$ is still present inside the system; any mismatch or initial condition excites it and an internal signal grows without bound. Stability must be judged on the internal (uncancelled) characteristic polynomial, never only on the simplified input–output transfer function.

> [!tip] Engineering Note — reading a system from its pole–zero plot
> Experienced designers read pole–zero plots the way others read schematics. Poles far to the left are fast and usually negligible; the **dominant poles** nearest the $j\omega$-axis set the response time. The angle of a complex pair from the negative real axis fixes the damping ($\zeta=\cos\theta$, so poles at $\pm45^\circ$ give $\zeta=0.707$, the Butterworth value). A right-half-plane zero (common in boost converters and in amplifiers with feed-forward capacitance) produces an initial undershoot and limits the achievable bandwidth. Root-locus, Bode and Nyquist methods — the standard tools of control design — are graphical ways of tracking these pole positions as a gain varies; the Nyquist criterion itself is an application of the argument principle of [Ch. 18](18-complex-analysis.md).

## Summary

| Concept | Formula / statement |
|---|---|
| Unilateral transform | $X(s)=\int_{0^-}^\infty x(t)e^{-st}dt$, $s=\sigma+j\omega$, ROC $\operatorname{Re}\{s\}>\sigma_c$ |
| Existence | piecewise continuous and $\lvert x\rvert\le Me^{\alpha t}$ $\Rightarrow$ $X$ exists for $\operatorname{Re}\{s\}>\alpha$, $X\to0$ |
| Uniqueness (Lerch) | equal transforms $\Rightarrow$ equal functions (at points of continuity) |
| Basic pairs | $\delta\leftrightarrow1$, $u\leftrightarrow\frac1s$, $t^n\leftrightarrow\frac{n!}{s^{n+1}}$, $e^{-at}\leftrightarrow\frac{1}{s+a}$ |
| Sinusoids | $\sin\omega t\leftrightarrow\frac{\omega}{s^2+\omega^2}$, $\cos\omega t\leftrightarrow\frac{s}{s^2+\omega^2}$ |
| Damped sinusoids | $e^{-\alpha t}\sin\omega t\leftrightarrow\frac{\omega}{(s+\alpha)^2+\omega^2}$, $e^{-\alpha t}\cos\omega t\leftrightarrow\frac{s+\alpha}{(s+\alpha)^2+\omega^2}$ |
| $s$-shift | $e^{at}x(t)\leftrightarrow X(s-a)$ |
| Delay | $x(t-T)u(t-T)\leftrightarrow e^{-sT}X(s)$; $g(t)u(t-T)\leftrightarrow e^{-sT}\mathcal{L}\{g(t+T)\}$ |
| Scaling | $x(at)\leftrightarrow\frac1aX(s/a)$ |
| Derivatives | $\dot x\leftrightarrow sX-x(0^-)$; $\ddot x\leftrightarrow s^2X-sx(0^-)-\dot x(0^-)$ |
| Integral | $\int_{0^-}^tx\,d\tau\leftrightarrow X/s$ |
| $t$-multiplication / division | $tx\leftrightarrow-X'(s)$; $x/t\leftrightarrow\int_s^\infty X\,dp$ |
| Periodic | $X=\frac{1}{1-e^{-sT}}\int_0^Txe^{-st}dt$ |
| Bromwich inversion | $x(t)=\frac{1}{2\pi j}\int_{\gamma-j\infty}^{\gamma+j\infty}X(s)e^{st}ds=\sum\operatorname{Res}[Xe^{st}]$ (Ch. 18) |
| Simple-pole coefficient | $K_k=[(s-p_k)X]_{s=p_k}=N(p_k)/D'(p_k)$ |
| Repeated pole | $K_{k,r-j}=\frac1{j!}\frac{d^j}{ds^j}[(s-p_k)^rX]_{p_k}$; $\frac{1}{(s-p)^i}\leftrightarrow\frac{t^{i-1}}{(i-1)!}e^{pt}$ |
| Complex pair $-\alpha\pm j\beta$ | $2\lvert K\rvert e^{-\alpha t}\cos(\beta t+\angle K)$ |
| Convolution | $(f*g)(t)=\int_0^tf(\tau)g(t-\tau)d\tau\leftrightarrow F(s)G(s)$ |
| Initial value | $x(0^+)=\lim_{s\to\infty}sX(s)$ ($X$ strictly proper) |
| Final value | $x(\infty)=\lim_{s\to0}sX(s)$ (all poles of $sX$ in open LHP) |
| ODE solution | $Y=\frac{C(s)}{P(s)}+\frac{X(s)}{P(s)}$: zero-input + zero-state |
| State equations | $\mathbf{X}=(s\mathbf{I}-\mathbf{A})^{-1}[\mathbf{x}(0^-)+\mathbf{b}V]$, $e^{\mathbf{A}t}=\mathcal{L}^{-1}\{(s\mathbf{I}-\mathbf{A})^{-1}\}$ |
| Element models | $Z_R=R$; $V_L=sLI-Li(0^-)$; $V_C=\frac{I}{sC}+\frac{v(0^-)}{s}$ |
| Transfer function | $H(s)=Y/X$ (zero state) $=\mathcal{L}\{h\}$; $y=h*x$ |
| Second-order form | $\frac{\omega_n^2}{s^2+2\zeta\omega_ns+\omega_n^2}$; poles $-\zeta\omega_n\pm j\omega_n\sqrt{1-\zeta^2}$ |
| Step overshoot, peak, settling | $e^{-\pi\zeta/\sqrt{1-\zeta^2}}$; $t_p=\pi/\omega_d$; $t_s\approx4/(\zeta\omega_n)$ |
| BIBO stability | all poles in open LHP $\Leftrightarrow\int_0^\infty\lvert h\rvert dt<\infty$ |
| Cubic Routh condition | $a_3s^3+a_2s^2+a_1s+a_0$, all $a_i>0$, stable iff $a_2a_1>a_3a_0$ |
| Feedback | $T=\frac{G}{1+GF}$; unity-feedback step error $\frac{1}{1+G(0)}$, ramp error $\frac{1}{\lim sG}$ |

## Further reading

- E. Kreyszig, *Advanced Engineering Mathematics*, 10th ed., Wiley — Ch. 6: a clear mathematical treatment of the transform, its properties and ODE applications, with extensive tables.
- A. V. Oppenheim, A. S. Willsky & S. H. Nawab, *Signals and Systems*, 2nd ed., Prentice Hall — Ch. 9: the bilateral transform, regions of convergence, and the system-function viewpoint.
- J. W. Nilsson & S. A. Riedel, *Electric Circuits*, Pearson — the chapters on the Laplace transform in circuit analysis: $s$-domain models, transfer functions and pole–zero interpretation.
- K. Ogata, *Modern Control Engineering*, 5th ed., Pearson — transient-response specifications, steady-state errors and the Routh–Hurwitz criterion.
- R. V. Churchill, *Operational Mathematics*, 3rd ed., McGraw-Hill — the classical rigorous source for existence, uniqueness and the inversion integral.
