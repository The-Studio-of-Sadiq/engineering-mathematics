# Chapter 5 — Integral Calculus

[← Previous](../part-2-calculus/04-differential-calculus.md) · [Table of Contents](../00-front/table-of-contents.md) · [Next →](../part-2-calculus/06-multivariable-calculus.md)

> **Prerequisites:** [Ch. 1](../part-1-foundations/01-mathematical-language-and-algebra.md), [Ch. 3](../part-1-foundations/03-trigonometry-and-periodic-functions.md), [Ch. 4](../part-2-calculus/04-differential-calculus.md)
> **Learning outcomes:** After this chapter you can
> - explain the definite integral as a limit of Riemann sums and use the Fundamental Theorem of Calculus to evaluate it through antiderivatives;
> - choose and apply the standard techniques — algebraic manipulation, trigonometric identities, partial fractions, trigonometric substitution, substitution and integration by parts — to reduce an integral to a table form;
> - decide whether an improper integral converges, using direct evaluation, comparison and limit-comparison tests;
> - set up and evaluate double and triple integrals over simple regions as iterated integrals;
> - approximate integrals with the trapezoidal and Simpson's rules and estimate their error;
> - compute charge, energy, work and other accumulated quantities, and the mean and RMS values of engineering waveforms.

Differentiation answers the question *how fast is it changing right now?* Integration answers the reverse question: *given the rate, how much has accumulated?* An electrical engineer meets this question constantly. The charge on a capacitor is the accumulated current, $q(t)=q(t_0)+\int_{t_0}^{t} i(\tau)\,d\tau$. The energy delivered to a load is the accumulated power, $W=\int p(t)\,dt$. The voltage a heater "feels" is the root-mean-square value of a waveform, which is the square root of an average, and an average is an integral divided by a length. The flux through a coil, the total charge in a cloud of electrons, the expected value of a noise voltage — all are integrals.

This chapter develops single-variable integration thoroughly and introduces integrals over regions of the plane and of space. [Chapter 6](../part-2-calculus/06-multivariable-calculus.md) takes multiple integrals further with changes of variables (polar, cylindrical, spherical, general Jacobians), and [Chapter 22](../part-7-numerical-computational/22-numerical-methods.md) treats numerical integration in depth (Romberg, Gaussian quadrature, adaptive schemes). Here we lay the foundations those chapters build on.

## 5.1 Indefinite and definite integrals

### 5.1.1 Motivation: area and accumulation

Suppose a current $i(t)$ flows into a capacitor. Over a short interval $\Delta t$ during which the current is nearly constant, the charge delivered is approximately $i(t)\,\Delta t$. Over a long interval $[a,b]$ we chop time into many short pieces and add:

$$
Q \approx \sum_{k=1}^{n} i(t_k^*)\,\Delta t_k .
$$

As the pieces shrink, the approximation becomes exact. Geometrically, each term $i(t_k^*)\,\Delta t_k$ is the area of a thin rectangle, so the limit is the *area under the curve* $i(t)$ between $t=a$ and $t=b$ (counted negative where $i<0$). This limit is the definite integral. It is a *number* attached to a function and an interval.

Separately, differentiation (Chapter 4) suggests the question: which functions $F$ have a given derivative $f$? The answer is a *family of functions*, the indefinite integral. The **Fundamental Theorem of Calculus** is the remarkable statement that these two ideas — area and anti-differentiation — are the same thing.

### 5.1.2 Riemann sums and the definite integral

> [!definition] Partition and Riemann sum
> A **partition** $P$ of $[a,b]$ is a finite set of points $a=x_0<x_1<\dots<x_n=b$. Write $\Delta x_k=x_k-x_{k-1}$ and $\lVert P\rVert=\max_k \Delta x_k$ (the **mesh**). Choosing a **sample point** $x_k^*\in[x_{k-1},x_k]$ in each subinterval, the **Riemann sum** of $f$ is
> $$
> S(f,P,\{x_k^*\})=\sum_{k=1}^{n} f(x_k^*)\,\Delta x_k .
> $$

> [!definition] Definite integral
> $f$ is **(Riemann) integrable** on $[a,b]$ if there is a number $I$ such that for every $\varepsilon>0$ there is a $\delta>0$ with $\lvert S(f,P,\{x_k^*\})-I\rvert<\varepsilon$ for *every* partition with $\lVert P\rVert<\delta$ and *every* choice of sample points. Then
> $$
> \int_a^b f(x)\,dx = I = \lim_{\lVert P\rVert\to 0}\sum_{k=1}^{n} f(x_k^*)\,\Delta x_k .
> $$
> $f$ is the **integrand**, $a$ and $b$ the **limits**, and $x$ a **dummy variable**: $\int_a^b f(x)\,dx=\int_a^b f(t)\,dt$.

The notation is a fossil of the sum: $\int$ is an elongated S (*summa*), and $dx$ recalls $\Delta x$. Dimensionally, $[\int f\,dx]=[f]\,[x]$; integrating a current in amperes over time in seconds gives coulombs.

Which functions are integrable? Not every function — the function equal to 1 on rationals and 0 on irrationals is the standard counterexample — but every function an engineer is likely to write down is.

> [!theorem] Sufficient conditions for integrability
> If $f$ is continuous on $[a,b]$, or bounded on $[a,b]$ with only finitely many discontinuities, or monotonic on $[a,b]$, then $f$ is Riemann integrable on $[a,b]$.

Thus square waves, sawtooth waves, rectified sinusoids and pulse trains — all piecewise continuous — can be integrated over any finite interval.

Two common special choices of Riemann sum use uniform width $h=(b-a)/n$ and sample at the left or right ends:

$$
L_n=h\sum_{k=0}^{n-1} f(a+kh),\qquad R_n=h\sum_{k=1}^{n} f(a+kh).
$$

For an increasing function $L_n\le\int_a^b f\le R_n$, and $R_n-L_n=h\,[f(b)-f(a)]\to0$, which proves directly that monotone functions are integrable.

> [!example] Worked Example 5.1 — An integral from the definition
> Compute $\displaystyle\int_0^2 x^2\,dx$ as a limit of right-endpoint sums.
>
> **Step 1 — partition.** $h=2/n$, $x_k=2k/n$.
>
> **Step 2 — sum.** Using $\sum_{k=1}^n k^2=\tfrac{n(n+1)(2n+1)}{6}$,
> $$
> R_n=\sum_{k=1}^{n}\Big(\frac{2k}{n}\Big)^2\frac{2}{n}=\frac{8}{n^3}\cdot\frac{n(n+1)(2n+1)}{6}=\frac{4}{3}\Big(1+\frac1n\Big)\Big(2+\frac1n\Big).
> $$
>
> **Step 3 — limit.** As $n\to\infty$, $R_n\to\tfrac43\cdot 2=\tfrac83\approx2.6667$.
>
> **Check.** $n=10$: $R_{10}=\tfrac43(1.1)(2.1)=3.08$; $n=100$: $R_{100}=\tfrac43(1.01)(2.01)=2.7068$; $n=1000$: $2.6707$. The error shrinks like $1/n$, as expected for an endpoint rule (§5.7).

```desmos-graph
left=-0.3; right=2.4; bottom=-0.5; top=4.5
---
y=x^2|blue
0<y<0.16|orange|0<x<0.4
0<y<0.64|orange|0.4<x<0.8
0<y<1.44|orange|0.8<x<1.2
0<y<2.56|orange|1.2<x<1.6
0<y<4|orange|1.6<x<2
```
*Right-endpoint Riemann sum with $n=5$ for $\int_0^2x^2\,dx$: the rectangles overshoot the exact area $8/3$ because $x^2$ is increasing.*

### 5.1.3 Properties of the definite integral

All of the following follow from the corresponding properties of finite sums, by passing to the limit.

> [!theorem] Algebraic properties
> For integrable $f,g$ and constants $\alpha,\beta$:
> 1. **Linearity:** $\displaystyle\int_a^b(\alpha f+\beta g)\,dx=\alpha\int_a^b f\,dx+\beta\int_a^b g\,dx$.
> 2. **Additivity over intervals:** $\displaystyle\int_a^c f\,dx=\int_a^b f\,dx+\int_b^c f\,dx$ for any $a,b,c$.
> 3. **Orientation:** $\displaystyle\int_b^a f\,dx=-\int_a^b f\,dx$, and $\displaystyle\int_a^a f\,dx=0$.
> 4. **Monotonicity:** if $f\le g$ on $[a,b]$ ($a<b$) then $\displaystyle\int_a^b f\,dx\le\int_a^b g\,dx$.
> 5. **Triangle inequality:** $\displaystyle\Big\lvert\int_a^b f\,dx\Big\rvert\le\int_a^b \lvert f\rvert\,dx$.
> 6. **Bounds:** if $m\le f(x)\le M$ on $[a,b]$ then $m(b-a)\le\displaystyle\int_a^b f\,dx\le M(b-a)$.

Property 5 is the integral form of the fact that a phasor sum can never exceed the sum of magnitudes; it is used constantly to bound errors. Property 6 leads to the following.

> [!theorem] Mean Value Theorem for integrals
> If $f$ is continuous on $[a,b]$, there is a $c\in(a,b)$ with
> $$
> \int_a^b f(x)\,dx=f(c)\,(b-a).
> $$
> The number $\bar f=\dfrac{1}{b-a}\displaystyle\int_a^b f(x)\,dx$ is the **mean value** of $f$ on $[a,b]$.

*Proof sketch.* By property 6 with $m=\min f$, $M=\max f$, the mean $\bar f$ lies between $m$ and $M$. A continuous function takes every value between its minimum and maximum (Intermediate Value Theorem, Ch. 4), so $\bar f=f(c)$ for some $c$. $\square$

The mean value is the "DC value" of a waveform and is developed fully in §5.9.

**Symmetry.** If $f$ is **even** ($f(-x)=f(x)$) then $\int_{-a}^a f\,dx=2\int_0^a f\,dx$; if $f$ is **odd** ($f(-x)=-f(x)$) then $\int_{-a}^{a}f\,dx=0$. If $f$ is **periodic** with period $T$, then $\int_{c}^{c+T}f\,dx$ is the same for every $c$ — you may integrate over any convenient full period. These three facts eliminate a large share of the work in Fourier analysis ([Ch. 15](../part-5-transforms/15-fourier-series.md)).

### 5.1.4 Antiderivatives and the indefinite integral

> [!definition] Antiderivative and indefinite integral
> $F$ is an **antiderivative** of $f$ on an interval $I$ if $F'(x)=f(x)$ for all $x\in I$. The set of all antiderivatives is the **indefinite integral**
> $$
> \int f(x)\,dx = F(x)+C,
> $$
> where $C$ is an arbitrary **constant of integration**.

Why only a constant? If $F'=G'=f$ on an interval, then $(F-G)'=0$, and by the Mean Value Theorem of differential calculus a function with zero derivative on an interval is constant. On a domain made of several disjoint intervals, the "constant" may differ on each piece: $\int \frac{dx}{x}=\ln\lvert x\rvert+C$ really means $\ln x+C_1$ for $x>0$ and $\ln(-x)+C_2$ for $x<0$.

Every differentiation rule read backwards gives an integration rule. The table below is the working vocabulary of this chapter; every entry can be checked by differentiating the right-hand side.

| $f(x)$ | $\int f(x)\,dx$ | Condition |
|---|---|---|
| $x^n$ | $\dfrac{x^{n+1}}{n+1}+C$ | $n\ne-1$ |
| $\dfrac1x$ | $\ln\lvert x\rvert+C$ | $x\ne0$ |
| $e^{ax}$ | $\dfrac{1}{a}e^{ax}+C$ | $a\ne0$ (real or complex) |
| $a^x$ | $\dfrac{a^x}{\ln a}+C$ | $a>0,\ a\ne1$ |
| $\sin ax$ | $-\dfrac1a\cos ax+C$ | |
| $\cos ax$ | $\dfrac1a\sin ax+C$ | |
| $\sec^2 x$ | $\tan x+C$ | |
| $\tan x$ | $-\ln\lvert\cos x\rvert+C$ | |
| $\sec x$ | $\ln\lvert\sec x+\tan x\rvert+C$ | |
| $\dfrac{1}{a^2+x^2}$ | $\dfrac1a\arctan\dfrac xa+C$ | $a>0$ |
| $\dfrac{1}{\sqrt{a^2-x^2}}$ | $\arcsin\dfrac xa+C$ | $\lvert x\rvert<a$ |
| $\dfrac{1}{\sqrt{x^2\pm a^2}}$ | $\ln\big\lvert x+\sqrt{x^2\pm a^2}\big\rvert+C$ | |
| $\dfrac{1}{x^2-a^2}$ | $\dfrac{1}{2a}\ln\Big\lvert\dfrac{x-a}{x+a}\Big\rvert+C$ | $\lvert x\rvert\ne a$ |
| $\sinh ax,\ \cosh ax$ | $\dfrac1a\cosh ax,\ \dfrac1a\sinh ax$ | |

The entry $\int e^{ax}dx=e^{ax}/a$ holds for complex $a$ as well, which is why phasor methods integrate so easily: $\int e^{j\omega t}dt=\dfrac{e^{j\omega t}}{j\omega}$. Integrating a sinusoid divides its phasor by $j\omega$ — a $-90^\circ$ phase shift and a scaling by $1/\omega$. This is exactly the capacitor law $V=I/(j\omega C)$ ([Ch. 2](../part-1-foundations/02-complex-numbers.md)).

### 5.1.5 The Fundamental Theorem of Calculus

> [!theorem] Fundamental Theorem of Calculus (FTC)
> Let $f$ be continuous on $[a,b]$.
>
> **Part 1.** The accumulation function $A(x)=\displaystyle\int_a^x f(t)\,dt$ is differentiable on $(a,b)$ and
> $$
> A'(x)=f(x).
> $$
> **Part 2.** If $F$ is any antiderivative of $f$ on $[a,b]$, then
> $$
> \int_a^b f(x)\,dx = F(b)-F(a) \equiv \Big[F(x)\Big]_a^b .
> $$

*Proof of Part 1.* For $h\ne0$,

$$
\frac{A(x+h)-A(x)}{h}=\frac1h\int_x^{x+h}f(t)\,dt=f(c_h)
$$

for some $c_h$ between $x$ and $x+h$, by the Mean Value Theorem for integrals. As $h\to0$, $c_h\to x$, and continuity of $f$ gives $f(c_h)\to f(x)$. $\square$

*Proof of Part 2.* $A$ and $F$ are both antiderivatives of $f$, so $F=A+C$. Then $F(b)-F(a)=A(b)-A(a)=\int_a^b f-0$. $\square$

Part 1 says integration followed by differentiation returns the original function: integration *creates* antiderivatives, so every continuous function has one (even if, like $e^{-x^2}$, it cannot be written with elementary functions). Part 2 turns the computation of a limit of sums into a subtraction.

Combining Part 1 with the chain rule gives **Leibniz's rule** for variable limits:

$$
\frac{d}{dx}\int_{u(x)}^{v(x)} f(t)\,dt = f\big(v(x)\big)\,v'(x)-f\big(u(x)\big)\,u'(x).
$$

> [!example] Worked Example 5.2 — Evaluating definite integrals with the FTC
> (a) $\displaystyle\int_0^2x^2\,dx=\Big[\frac{x^3}{3}\Big]_0^2=\frac83$, agreeing with Example 5.1 with no limits needed.
>
> (b) $\displaystyle\int_0^{\pi}\sin x\,dx=\big[-\cos x\big]_0^{\pi}=-(-1)-(-1)=2.$
>
> (c) $\displaystyle\int_1^{e}\frac{3}{x}\,dx=3\big[\ln x\big]_1^{e}=3.$
>
> (d) $\displaystyle\int_0^{1}\frac{dx}{1+x^2}=\big[\arctan x\big]_0^1=\frac{\pi}{4}\approx0.785398.$
>
> (e) $\displaystyle\frac{d}{dx}\int_0^{x^2}\cos t\,dt=\cos(x^2)\cdot2x$ by Leibniz's rule. Check: the integral equals $\sin(x^2)$, whose derivative is $2x\cos(x^2)$. ✓

> [!example] Worked Example 5.3 — Capacitor voltage from a current waveform
> A $C=10\ \mu\text{F}$ capacitor, initially at $v(0)=2$ V, is driven by $i(t)=5e^{-200t}$ mA for $t\ge0$. Find $v(t)$ and the final voltage.
>
> **Step 1 — accumulation law.** $v(t)=v(0)+\dfrac1C\displaystyle\int_0^t i(\tau)\,d\tau$.
>
> **Step 2 — integrate.**
> $$
> \int_0^t 5\times10^{-3}e^{-200\tau}d\tau=5\times10^{-3}\Big[\frac{e^{-200\tau}}{-200}\Big]_0^t=2.5\times10^{-5}\big(1-e^{-200t}\big)\ \text{C}.
> $$
>
> **Step 3 — divide by $C$.** $v(t)=2+\dfrac{2.5\times10^{-5}}{10^{-5}}(1-e^{-200t})=2+2.5(1-e^{-200t})$ V.
>
> **Step 4 — final value.** As $t\to\infty$, $v\to4.5$ V. Units: $[\text{A}\cdot\text{s}/\text{F}]=[\text{C}/\text{F}]=[\text{V}]$ ✓. At $t=5$ ms, $v=2+2.5(1-e^{-1})=3.580$ V.

> [!example] Worked Example 5.4 — Integrating a piecewise waveform
> A current is $i(t)=2$ A for $0\le t<1$ ms, $-1$ A for $1\le t<3$ ms, and $0$ afterwards. Find the net charge $q(t)$ delivered from $t=0$.
>
> Split the integral at the breakpoints (additivity):
> - $0\le t<1$ ms: $q=2t$, reaching $2$ mC at $t=1$ ms.
> - $1\le t<3$ ms: $q=2\ \text{mC}-1\cdot(t-1\ \text{ms})$, reaching $2-2=0$ mC at $t=3$ ms.
> - $t\ge3$ ms: $q=0$.
>
> The charge is a continuous piecewise-linear function even though the current jumps. **Integration smooths**: it turns jumps into kinks and kinks into smooth bends. This is why a capacitor voltage cannot jump while its current can.

```desmos-graph
left=-0.5; right=4.5; bottom=-1.5; top=2.5
---
y=2|red|0<x<1
y=-1|red|1<x<3
y=0|red|3<x<4.5
y=2x|blue|0<x<1
y=2-(x-1)|blue|1<x<3
```
*Current $i(t)$ in A (red) and accumulated charge $q(t)$ in mC (blue), time in ms: the integral of a discontinuous function is continuous.*

> [!warning] Common Pitfall — Applying the FTC across a discontinuity or singularity
> $\int_{-1}^{1}\frac{dx}{x^2}=\big[-\frac1x\big]_{-1}^{1}=-2$ is **wrong**: a positive integrand cannot have a negative integral. The FTC requires $F'=f$ on the *whole* interval, and $1/x^2$ blows up at $x=0$. The integral is improper and in fact diverges (§5.3). Always scan the interval for points where the integrand is undefined before subtracting antiderivative values.

> [!tip] Engineering Note — Integrators in hardware
> An op-amp integrator produces $v_o(t)=-\frac{1}{RC}\int_0^t v_i(\tau)\,d\tau+v_o(0)$. The "constant of integration" is the initial capacitor voltage, which in a real circuit drifts because of input offset voltage: integrating a constant offset $V_{os}$ gives a ramp $-V_{os}t/(RC)$ that eventually saturates the amplifier. Practical integrators add a large resistor across the capacitor to bleed this off — mathematically, replacing $\int$ by a leaky integrator with a finite time constant.

## 5.2 Integration techniques

### 5.2.1 Why technique matters

Differentiation is mechanical: the product, quotient and chain rules dispose of any elementary expression. Integration is not. There is no product rule for integrals, and some innocent-looking integrands — $e^{-x^2}$, $\frac{\sin x}{x}$, $\sqrt{1-k^2\sin^2x}$ — have no elementary antiderivative at all. What we have instead is a toolkit of *transformations* that convert an integral into one already in the table. This section surveys the toolkit and the algebraic and trigonometric moves; the two most powerful tools, **integration by parts** and **integration by substitution**, each get a full section (§5.5 and §5.6), and the numerical fallback when everything fails is §5.7.

Even in an age of computer algebra, an engineer needs these techniques: to recognise the *form* of an answer (logarithm? arctangent? decaying exponential?), to sanity-check software output, to carry symbolic parameters such as $R$, $L$, $C$ through a calculation, and to see the structure that a numerical answer hides. Partial fractions, in particular, are the engine of inverse Laplace and Z-transforms ([Ch. 14](../part-5-transforms/14-laplace-transforms.md), [Ch. 17](../part-5-transforms/17-z-transforms.md)).

### 5.2.2 A strategy

Faced with $\int f(x)\,dx$, work down this list:

1. **Simplify algebraically.** Expand products, divide out improper fractions, split sums, cancel common factors, use $\sqrt{x^2}=\lvert x\rvert$ carefully.
2. **Match the table** (§5.1.4) after linearity, possibly with a linear change $x\to ax+b$.
3. **Spot a derivative.** If the integrand has the form $g(h(x))\,h'(x)$, substitute $u=h(x)$ (§5.6). Special case: $\int\frac{h'(x)}{h(x)}dx=\ln\lvert h(x)\rvert+C$.
4. **Trigonometric identities** for powers and products of sines and cosines (§5.2.4).
5. **Partial fractions** for rational functions (§5.2.5).
6. **Trigonometric or hyperbolic substitution** for $\sqrt{a^2- x^2}$, $\sqrt{a^2+x^2}$, $\sqrt{x^2-a^2}$ (§5.6.4).
7. **Integration by parts** for products of unlike functions — polynomial × exponential, polynomial × trig, logarithms, inverse trig (§5.5).
8. **Complex exponentials** for products of exponentials and sinusoids (§5.2.6).
9. **Numerical integration** when no closed form exists or only a number is needed (§5.7).

### 5.2.3 Algebraic manipulation

Many integrals yield to algebra alone.

> [!example] Worked Example 5.5 — Algebra first
> (a) $\displaystyle\int\frac{x^3+2x-\sqrt{x}}{x^2}\,dx=\int\Big(x+\frac{2}{x}-x^{-3/2}\Big)dx=\frac{x^2}{2}+2\ln\lvert x\rvert+2x^{-1/2}+C.$
>
> (b) $\displaystyle\int\frac{x^2}{x+1}\,dx$. The fraction is *improper* (degree of numerator $\ge$ degree of denominator), so divide first: $x^2=(x+1)(x-1)+1$, hence
> $$
> \int\frac{x^2}{x+1}dx=\int\Big(x-1+\frac{1}{x+1}\Big)dx=\frac{x^2}{2}-x+\ln\lvert x+1\rvert+C.
> $$
>
> (c) $\displaystyle\int\frac{dx}{x^2+4x+13}$. **Complete the square:** $x^2+4x+13=(x+2)^2+9$, so
> $$
> \int\frac{dx}{(x+2)^2+3^2}=\frac13\arctan\frac{x+2}{3}+C.
> $$
>
> (d) $\displaystyle\int_0^1\frac{2x+1}{x^2+x+1}dx$. The numerator is exactly the derivative of the denominator, so the answer is $\big[\ln(x^2+x+1)\big]_0^1=\ln3\approx1.0986$.

Completing the square deserves emphasis because it is how quadratic denominators — which in circuit theory come from second-order (RLC) characteristic polynomials $s^2+2\zeta\omega_n s+\omega_n^2$ — are put into a recognisable form:

$$
s^2+2\zeta\omega_n s+\omega_n^2=(s+\zeta\omega_n)^2+\omega_n^2(1-\zeta^2).
$$

When $\zeta<1$ the second term is positive and an arctangent (in the time domain: a damped sinusoid) appears; when $\zeta>1$ it is negative, the quadratic factors into real roots, and logarithms (decaying exponentials) appear instead.

### 5.2.4 Trigonometric integrals

The integrals of products and powers of sinusoids are the backbone of power calculations and Fourier analysis. The identities that drive them are ([Ch. 3](../part-1-foundations/03-trigonometry-and-periodic-functions.md)):

$$
\cos^2\theta=\frac{1+\cos2\theta}{2},\qquad \sin^2\theta=\frac{1-\cos2\theta}{2},\qquad \sin\theta\cos\theta=\frac{\sin2\theta}{2},
$$

$$
\begin{aligned}
\sin A\cos B&=\tfrac12\big[\sin(A+B)+\sin(A-B)\big],\\
\cos A\cos B&=\tfrac12\big[\cos(A-B)+\cos(A+B)\big],\\
\sin A\sin B&=\tfrac12\big[\cos(A-B)-\cos(A+B)\big].
\end{aligned}
$$

**Rules for $\int\sin^m x\cos^n x\,dx$:**

- If $m$ is **odd**, save one $\sin x$, convert the remaining even power of $\sin$ to $\cos$ using $\sin^2=1-\cos^2$, and substitute $u=\cos x$.
- If $n$ is **odd**, save one $\cos x$, convert the rest to $\sin$, and substitute $u=\sin x$.
- If both are **even**, reduce powers with the half-angle identities, repeatedly if necessary.

> [!example] Worked Example 5.6 — Powers of sinusoids
> (a) **Odd power.** $\displaystyle\int\sin^3x\,dx=\int(1-\cos^2x)\sin x\,dx$. With $u=\cos x$, $du=-\sin x\,dx$:
> $$
> -\int(1-u^2)\,du=-u+\frac{u^3}{3}+C=-\cos x+\frac{\cos^3x}{3}+C.
> $$
>
> (b) **Even powers.** $\displaystyle\int\sin^2x\cos^2x\,dx=\int\frac{\sin^22x}{4}dx=\int\frac{1-\cos4x}{8}dx=\frac{x}{8}-\frac{\sin4x}{32}+C.$
>
> (c) **Over a period.** For $\omega=2\pi/T$,
> $$
> \int_0^T\sin^2(\omega t)\,dt=\int_0^T\frac{1-\cos2\omega t}{2}dt=\frac T2-\Big[\frac{\sin2\omega t}{4\omega}\Big]_0^T=\frac T2 .
> $$
> The cosine term integrates to zero over a whole number of periods. This single result gives the RMS value $V_m/\sqrt2$ of a sinusoid (§5.9).
>
> **Check (b) numerically:** $\int_0^{1}\sin^2x\cos^2x\,dx=\tfrac18-\tfrac{\sin4}{32}=0.125+0.02365=0.14865$; adaptive quadrature gives $0.148650$ ✓.

> [!theorem] Orthogonality of sinusoids
> For integers $m,n\ge1$ and $T=2\pi/\omega_0$:
> $$
> \int_0^T\sin(m\omega_0t)\sin(n\omega_0t)\,dt=\int_0^T\cos(m\omega_0t)\cos(n\omega_0t)\,dt=\begin{cases}T/2,&m=n\\0,&m\ne n\end{cases}
> $$
> and $\displaystyle\int_0^T\sin(m\omega_0t)\cos(n\omega_0t)\,dt=0$ for all $m,n$.

*Proof.* Use the product-to-sum identity. For $m\ne n$, $\sin m\theta\sin n\theta=\tfrac12[\cos(m-n)\theta-\cos(m+n)\theta]$, and each cosine of a nonzero integer multiple of $\theta=\omega_0 t$ completes a whole number of cycles in $[0,T]$, integrating to zero. For $m=n$ the first cosine becomes $\cos0=1$, giving $T/2$. The other cases are identical. $\square$

Orthogonality is the reason power in different harmonics simply *adds* (Parseval's theorem, §5.9) and the reason Fourier coefficients can be computed one at a time ([Ch. 15](../part-5-transforms/15-fourier-series.md)).

```desmos-graph
left=0; right=6.2832; bottom=-1.2; top=1.2
---
y=\sin(2x)\sin(3x)|blue
y=\sin(3x)^2|red|dashed
y=0.5|red|dotted
```
*$\sin2x\sin3x$ (blue) has equal positive and negative area over $[0,2\pi]$; $\sin^23x$ (red) oscillates about its mean $1/2$.*

### 5.2.5 Partial fractions

A **rational function** $P(x)/Q(x)$ (ratio of polynomials) can always be integrated in elementary terms. The method:

1. If $\deg P\ge\deg Q$, divide to get a polynomial plus a **proper** fraction.
2. Factor $Q$ over the reals into linear factors $(x-a)$ and irreducible quadratics $(x^2+bx+c)$ with $b^2<4c$. (Over $\mathbb{C}$ every factor is linear; the real quadratics pair complex-conjugate roots.)
3. Write the proper fraction as a sum of simple terms according to this table.

| Factor of $Q$ | Terms in the expansion |
|---|---|
| $(x-a)$, simple | $\dfrac{A}{x-a}$ |
| $(x-a)^k$, repeated | $\dfrac{A_1}{x-a}+\dfrac{A_2}{(x-a)^2}+\dots+\dfrac{A_k}{(x-a)^k}$ |
| $x^2+bx+c$, irreducible | $\dfrac{Bx+C}{x^2+bx+c}$ |
| $(x^2+bx+c)^k$, repeated | $\displaystyle\sum_{i=1}^{k}\frac{B_ix+C_i}{(x^2+bx+c)^i}$ |

4. Find the constants. Fastest for simple linear factors is the **cover-up (Heaviside) rule**:
$$
A=\Big[(x-a)\frac{P(x)}{Q(x)}\Big]_{x=a}.
$$
For repeated factors, cover-up gives the highest-power coefficient $A_k$; the others follow by comparing coefficients, substituting convenient $x$ values, or differentiating: $A_{k-i}=\dfrac{1}{i!}\dfrac{d^i}{dx^i}\Big[(x-a)^kP/Q\Big]_{x=a}$.
5. Integrate each term: $\int\frac{dx}{x-a}=\ln\lvert x-a\rvert$, $\int\frac{dx}{(x-a)^k}=\frac{(x-a)^{1-k}}{1-k}$ for $k\ge2$, and quadratic terms by completing the square (split into a $\ln$ part and an $\arctan$ part).

*Why the cover-up rule works.* If $\frac{P}{Q}=\frac{A}{x-a}+R(x)$ where $R$ has no pole at $a$, multiplying by $(x-a)$ gives $(x-a)\frac PQ=A+(x-a)R(x)$; setting $x=a$ kills the second term.

> [!example] Worked Example 5.7 — Distinct linear factors
> $\displaystyle\int\frac{3x+5}{(x+1)(x+2)}dx$.
>
> Cover-up: $A=\dfrac{3(-1)+5}{-1+2}=2$, $B=\dfrac{3(-2)+5}{-2+1}=1$. So
> $$
> \int\Big(\frac{2}{x+1}+\frac{1}{x+2}\Big)dx=2\ln\lvert x+1\rvert+\ln\lvert x+2\rvert+C.
> $$
> **Check** at $x=0$: $\frac{2}{1}+\frac{1}{2}=2.5=\frac{5}{2}$ ✓.

> [!example] Worked Example 5.8 — A repeated factor
> $\displaystyle\int\frac{x^2+1}{x(x-1)^2}dx$. Write $\dfrac{x^2+1}{x(x-1)^2}=\dfrac{A}{x}+\dfrac{B}{x-1}+\dfrac{C}{(x-1)^2}$.
>
> - Cover $x$, set $x=0$: $A=\dfrac{1}{(-1)^2}=1$.
> - Cover $(x-1)^2$, set $x=1$: $C=\dfrac{1+1}{1}=2$.
> - Compare $x^2$ coefficients in $x^2+1=A(x-1)^2+Bx(x-1)+Cx$: $1=A+B$, so $B=0$.
>
> $$
> \int\Big(\frac1x+\frac{2}{(x-1)^2}\Big)dx=\ln\lvert x\rvert-\frac{2}{x-1}+C.
> $$
> **Check** at $x=2$: integrand $=5/2$; $\frac12+\frac{2}{1}=2.5$ ✓.

> [!example] Worked Example 5.9 — An irreducible quadratic
> $\displaystyle\int_0^1\frac{5x+3}{(x+1)(x^2+1)}dx$.
>
> **Step 1 — form.** $\dfrac{5x+3}{(x+1)(x^2+1)}=\dfrac{A}{x+1}+\dfrac{Bx+C}{x^2+1}$.
>
> **Step 2 — constants.** Cover-up at $x=-1$: $A=\dfrac{-5+3}{2}=-1$. Then $5x+3=-(x^2+1)+(Bx+C)(x+1)$. Comparing $x^2$: $0=-1+B\Rightarrow B=1$. Constants: $3=-1+C\Rightarrow C=4$. (Check $x$ terms: $5=B+C$ ✓.)
>
> **Step 3 — integrate.**
> $$
> \int_0^1\Big(-\frac{1}{x+1}+\frac{x}{x^2+1}+\frac{4}{x^2+1}\Big)dx=\Big[-\ln(x+1)+\tfrac12\ln(x^2+1)+4\arctan x\Big]_0^1 .
> $$
>
> **Step 4 — evaluate.** $-\ln2+\tfrac12\ln2+4\cdot\tfrac\pi4=\pi-\tfrac12\ln2\approx3.14159-0.34657=2.79502$.
>
> **Check:** numerical quadrature gives $2.795019$ ✓.

> [!example] Worked Example 5.10 — A circuit-flavoured rational integral
> The integral $I=\displaystyle\int_0^{\infty}\frac{d\omega}{(\omega^2+a^2)(\omega^2+b^2)}$ ($a,b>0$, $a\ne b$) appears in noise-bandwidth calculations for cascaded RC filters. Treat $u=\omega^2$ as the variable: cover-up on $\frac{1}{(u+a^2)(u+b^2)}$ gives
> $$
> \frac{1}{(\omega^2+a^2)(\omega^2+b^2)}=\frac{1}{b^2-a^2}\Big(\frac{1}{\omega^2+a^2}-\frac{1}{\omega^2+b^2}\Big).
> $$
> Each term integrates to an arctangent; on $[0,\infty)$, $\int_0^\infty\frac{d\omega}{\omega^2+a^2}=\frac{\pi}{2a}$ (an improper integral, §5.3). Hence
> $$
> I=\frac{1}{b^2-a^2}\Big(\frac{\pi}{2a}-\frac{\pi}{2b}\Big)=\frac{\pi}{2ab(a+b)} .
> $$
> **Check** with $a=1,b=2$: $\pi/12=0.261799$; numerical quadrature agrees.

> [!tip] Engineering Note — Partial fractions and system response
> In the Laplace domain the response of a linear circuit is a rational function $Y(s)=P(s)/Q(s)$. The roots of $Q$ are the **poles** (natural frequencies). Expanding $Y$ in partial fractions decomposes the response into **modes**: each simple real pole $p$ contributes $Ae^{pt}$, each repeated pole contributes $t^{k-1}e^{pt}$ terms, and each conjugate pair contributes a damped sinusoid. The cover-up coefficient $A$ is called the **residue**, a name that returns in [Ch. 18](../part-5-transforms/18-complex-analysis.md).

### 5.2.6 The complex-exponential method

Integrals of the form $\int e^{at}\cos bt\,dt$ and $\int e^{at}\sin bt\,dt$ can be done by parts twice (§5.5), but the fastest route uses Euler's formula. Since $e^{at}\cos bt=\operatorname{Re}\{e^{(a+jb)t}\}$ and integration commutes with taking real parts,

$$
\int e^{at}\cos bt\,dt=\operatorname{Re}\Big\{\frac{e^{(a+jb)t}}{a+jb}\Big\}=\operatorname{Re}\Big\{\frac{(a-jb)\,e^{at}(\cos bt+j\sin bt)}{a^2+b^2}\Big\},
$$

so

$$
\int e^{at}\cos bt\,dt=\frac{e^{at}\,(a\cos bt+b\sin bt)}{a^2+b^2}+C,\qquad
\int e^{at}\sin bt\,dt=\frac{e^{at}\,(a\sin bt-b\cos bt)}{a^2+b^2}+C,
$$

the second from the imaginary part.

> [!example] Worked Example 5.11 — Charge in a damped oscillation
> The current in an underdamped RLC circuit is $i(t)=e^{-2t}\sin 4t$ A (time in seconds). Find $\int_0^\infty i(t)\,dt$, the net charge transferred.
>
> With $a=-2$, $b=4$, $a^2+b^2=20$:
> $$
> \int_0^\infty e^{-2t}\sin4t\,dt=\Big[\frac{e^{-2t}(-2\sin4t-4\cos4t)}{20}\Big]_0^\infty=0-\frac{-4}{20}=0.2\ \text{C}.
> $$
> **Shortcut:** $\int_0^\infty e^{(a+jb)t}dt=\frac{-1}{a+jb}$ for $a<0$, so the answer is $\operatorname{Im}\{\frac{1}{2-4j}\}=\operatorname{Im}\{\frac{2+4j}{20}\}=0.2$ ✓. This is the Laplace transform $\frac{4}{s^2+4s+20}$ evaluated at $s=0$ — an early glimpse of [Ch. 14](../part-5-transforms/14-laplace-transforms.md).

> [!warning] Common Pitfall — "Re" and products
> $\operatorname{Re}\{\cdot\}$ commutes with integration and with multiplication by *real* quantities, but not with multiplication of two complex quantities: $\operatorname{Re}\{z_1z_2\}\ne\operatorname{Re}\{z_1\}\operatorname{Re}\{z_2\}$ in general. When computing average power $\frac1T\int v(t)i(t)\,dt$ from phasors you cannot simply multiply phasors; the correct result is $\tfrac12\operatorname{Re}\{\mathbf{V}\mathbf{I}^*\}$ (§5.9).

### 5.2.7 Integrals with no elementary antiderivative

Liouville's theory (19th century) shows that certain integrals cannot be expressed using finitely many elementary functions. Important examples define *new* special functions:

| Integral | Name | Where it appears |
|---|---|---|
| $\dfrac{2}{\sqrt\pi}\displaystyle\int_0^x e^{-t^2}dt$ | error function $\operatorname{erf}(x)$ | Gaussian noise, bit-error rate, diffusion |
| $\displaystyle\int_0^x\frac{\sin t}{t}dt$ | sine integral $\operatorname{Si}(x)$ | Gibbs phenomenon, ideal-filter step response |
| $\displaystyle\int_0^{\pi/2}\frac{d\theta}{\sqrt{1-k^2\sin^2\theta}}$ | complete elliptic integral $K(k)$ | pendulum period, coil inductance, line geometries |
| $\displaystyle\int_x^\infty\frac{e^{-t}}{t}dt$ | exponential integral $E_1(x)$ | radiation, transient heating |

These are tabulated and built into every numerical library. Recognising that an integral *is* one of them is a useful skill; trying to "find" an elementary antiderivative for them is not.

## 5.3 Improper integrals

### 5.3.1 Motivation

The Riemann integral is defined for a *bounded* function on a *bounded* interval. Engineering keeps breaking both conditions. The total energy dissipated as a capacitor discharges is an integral from $t=0$ to $t=\infty$. The Laplace and Fourier transforms are integrals over infinite ranges. The field of a line charge involves $\int\frac{dz}{(r^2+z^2)^{3/2}}$ over the whole $z$-axis. The probability that Gaussian noise exceeds a threshold is a tail integral to infinity. And the potential near a point charge involves integrands that blow up like $1/R$. **Improper integrals** extend the definition to these cases by taking limits.

### 5.3.2 Definitions

> [!definition] Improper integral of Type I (infinite interval)
> If $f$ is integrable on $[a,b]$ for every $b>a$,
> $$
> \int_a^\infty f(x)\,dx=\lim_{b\to\infty}\int_a^b f(x)\,dx ,
> $$
> provided the limit exists and is finite; the integral then **converges**, otherwise it **diverges**. Similarly $\int_{-\infty}^b f=\lim_{a\to-\infty}\int_a^b f$, and
> $$
> \int_{-\infty}^{\infty}f(x)\,dx=\int_{-\infty}^{c}f(x)\,dx+\int_{c}^{\infty}f(x)\,dx ,
> $$
> which converges only if **both** halves converge independently (the choice of $c$ does not matter).

> [!definition] Improper integral of Type II (unbounded integrand)
> If $f$ is integrable on $[a+\varepsilon,b]$ for every $\varepsilon>0$ but unbounded near $a$,
> $$
> \int_a^b f(x)\,dx=\lim_{\varepsilon\to0^+}\int_{a+\varepsilon}^{b}f(x)\,dx .
> $$
> Likewise at the right end. If the singularity is at an interior point $c$, split: $\int_a^b=\int_a^c+\int_c^b$, and require both to converge.

An integral can be improper for several reasons at once ($\int_0^\infty\frac{dx}{\sqrt x}$ is improper at both ends); split it so that each piece has only one difficulty.

### 5.3.3 The benchmark $p$-integrals

Almost every convergence question is settled by comparison with powers of $x$.

> [!theorem] $p$-integrals
> $$
> \int_1^\infty\frac{dx}{x^p}\ \begin{cases}=\dfrac{1}{p-1},&p>1\\[4pt]\text{diverges},&p\le1\end{cases}
> \qquad\qquad
> \int_0^1\frac{dx}{x^p}\ \begin{cases}=\dfrac{1}{1-p},&p<1\\[4pt]\text{diverges},&p\ge1\end{cases}
> $$

*Proof.* For $p\ne1$, $\int_1^b x^{-p}dx=\frac{b^{1-p}-1}{1-p}$. As $b\to\infty$, $b^{1-p}\to0$ if $p>1$ and $\to\infty$ if $p<1$. For $p=1$, $\int_1^b\frac{dx}{x}=\ln b\to\infty$. The $[0,1]$ case is the same computation with $\varepsilon^{1-p}$ as $\varepsilon\to0^+$. $\square$

The memory aid: **at infinity a function must decay faster than $1/x$; at a singularity it must blow up slower than $1/x$.** The function $1/x$ itself fails both tests — its integral diverges logarithmically at both ends. Exponentials decay faster than any power, so $\int_0^\infty x^n e^{-ax}dx$ converges for every $n\ge0$, $a>0$.

```desmos-graph
left=0; right=8; bottom=0; top=3
---
y=\frac{1}{x}|red
y=\frac{1}{x^2}|blue
y=\frac{1}{\sqrt{x}}|green|dashed
x=1|grey|dotted
```
*$1/x$ (red) sits on the borderline: its area diverges both near $0$ and towards $\infty$. $1/x^2$ (blue) has finite area on $[1,\infty)$; $1/\sqrt{x}$ (green) has finite area on $(0,1]$.*

> [!example] Worked Example 5.12 — Direct evaluation
> (a) $\displaystyle\int_0^\infty e^{-t/\tau}dt=\lim_{b\to\infty}\big[-\tau e^{-t/\tau}\big]_0^b=\tau$ for $\tau>0$. The area under a decaying exponential equals its time constant — the basis of the "area = height × time constant" rule of thumb.
>
> (b) $\displaystyle\int_0^1\frac{dx}{\sqrt x}=\lim_{\varepsilon\to0^+}\big[2\sqrt x\big]_\varepsilon^1=2.$ A function can be infinite at a point and still enclose finite area.
>
> (c) $\displaystyle\int_0^1\ln x\,dx=\lim_{\varepsilon\to0^+}\big[x\ln x-x\big]_\varepsilon^1=-1-\lim_{\varepsilon\to0^+}(\varepsilon\ln\varepsilon-\varepsilon)=-1$, using $\varepsilon\ln\varepsilon\to0$ (L'Hôpital, Ch. 4).
>
> (d) $\displaystyle\int_{-\infty}^\infty\frac{dx}{1+x^2}=\lim_{a\to-\infty}\big[\arctan x\big]_a^0+\lim_{b\to\infty}\big[\arctan x\big]_0^b=\frac\pi2+\frac\pi2=\pi.$

> [!example] Worked Example 5.13 — Energy dissipated by a discharging capacitor
> A capacitor $C$ charged to $V_0$ discharges through $R$: $v(t)=V_0e^{-t/RC}$. Find the total energy dissipated in $R$.
>
> **Step 1 — power.** $p(t)=v^2/R=\dfrac{V_0^2}{R}e^{-2t/RC}$.
>
> **Step 2 — improper integral.** Using Example 5.12(a) with $\tau=RC/2$,
> $$
> W=\int_0^\infty\frac{V_0^2}{R}e^{-2t/RC}dt=\frac{V_0^2}{R}\cdot\frac{RC}{2}=\frac12CV_0^2 .
> $$
>
> **Interpretation.** All of the stored energy $\tfrac12CV_0^2$ ends up as heat, independent of $R$ — a larger resistor dissipates less power but for longer. With $C=100\ \mu$F, $V_0=12$ V: $W=7.2$ mJ.

> [!example] Worked Example 5.14 — A divergent integral in disguise
> Is $\displaystyle\int_0^2\frac{dx}{(x-1)^2}$ finite?
>
> The integrand is singular at the interior point $x=1$. Split at $1$:
> $$
> \int_0^1\frac{dx}{(x-1)^2}=\lim_{\varepsilon\to0^+}\Big[-\frac{1}{x-1}\Big]_0^{1-\varepsilon}=\lim_{\varepsilon\to0^+}\Big(\frac1\varepsilon-1\Big)=\infty .
> $$
> One half diverges, so the whole integral diverges. Blindly applying the FTC gives $\big[-\frac1{x-1}\big]_0^2=-1-1=-2$, a negative "area" for a positive function — the same trap as in §5.1.

### 5.3.4 Comparison tests

Often we only need to know *whether* an improper integral converges — for example, to confirm that a Fourier or Laplace transform exists — without computing it.

> [!theorem] Direct comparison test
> Suppose $0\le f(x)\le g(x)$ for $x\ge a$.
> - If $\int_a^\infty g\,dx$ converges, so does $\int_a^\infty f\,dx$ (and $\int_a^\infty f\le\int_a^\infty g$).
> - If $\int_a^\infty f\,dx$ diverges, so does $\int_a^\infty g\,dx$.
>
> The same holds for Type II integrals near a singularity.

*Proof sketch.* $F(b)=\int_a^b f$ is non-decreasing in $b$ (since $f\ge0$) and bounded above by $\int_a^\infty g$. A bounded monotone function has a limit. $\square$

> [!theorem] Limit comparison test
> If $f,g>0$ for $x\ge a$ and $\displaystyle\lim_{x\to\infty}\frac{f(x)}{g(x)}=L$ with $0<L<\infty$, then $\int_a^\infty f$ and $\int_a^\infty g$ either both converge or both diverge.

Limit comparison formalises "only the dominant behaviour matters". To test $\int_1^\infty\frac{3x+2}{x^3+5}dx$, note the integrand behaves like $3x/x^3=3/x^2$ for large $x$; the ratio to $1/x^2$ tends to $3$, and $\int_1^\infty x^{-2}$ converges, so the original converges.

> [!theorem] Absolute convergence
> If $\int_a^\infty\lvert f\rvert\,dx$ converges then $\int_a^\infty f\,dx$ converges, and $\big\lvert\int_a^\infty f\big\rvert\le\int_a^\infty\lvert f\rvert$.

An integral that converges but not absolutely is **conditionally convergent**. The classic example is the Dirichlet integral

$$
\int_0^\infty\frac{\sin x}{x}\,dx=\frac{\pi}{2},
$$

which converges because the positive and negative lobes of $\sin x/x$ shrink and alternate, yet $\int_0^\infty\big\lvert\frac{\sin x}{x}\big\rvert dx=\infty$ (each lobe has area roughly $2/(k\pi)$, and $\sum 1/k$ diverges). This integral is the area under the $\operatorname{sinc}$ function and underlies the Gibbs phenomenon and the ideal low-pass filter ([Ch. 16](../part-5-transforms/16-fourier-transforms.md)).

> [!example] Worked Example 5.15 — Convergence decisions
> (a) $\displaystyle\int_1^\infty e^{-x^2}dx$: for $x\ge1$, $x^2\ge x$, so $0<e^{-x^2}\le e^{-x}$, and $\int_1^\infty e^{-x}dx=e^{-1}$. **Converges**, with value at most $0.3679$ (actual: $0.1394$).
>
> (b) $\displaystyle\int_1^\infty\frac{\sqrt{x}+1}{x^2}dx$: dominant term $x^{1/2}/x^2=x^{-3/2}$, $p=3/2>1$. **Converges** (in fact it equals $2+1=3$).
>
> (c) $\displaystyle\int_2^\infty\frac{dx}{\ln x}$: for $x\ge2$, $\ln x<x$, so $\frac1{\ln x}>\frac1x$, and $\int_2^\infty\frac{dx}{x}$ diverges. **Diverges.**
>
> (d) $\displaystyle\int_0^1\frac{\cos x}{x^{2/3}}dx$: near $0$, $\cos x\to1$, so the integrand behaves like $x^{-2/3}$ with $p=2/3<1$. **Converges.**
>
> (e) $\displaystyle\int_0^\infty\frac{\sin^2x}{x^2}dx$: near $0$ the integrand $\to1$ (no singularity at all); for large $x$ it is bounded by $1/x^2$. **Converges** absolutely (its value is $\pi/2$).

### 5.3.5 Cauchy principal value

Some integrals over symmetric ranges diverge in the strict sense but have a natural "balanced" value.

> [!definition] Cauchy principal value
> $$
> \operatorname{PV}\!\int_{-\infty}^{\infty}f(x)\,dx=\lim_{R\to\infty}\int_{-R}^{R}f(x)\,dx,\qquad
> \operatorname{PV}\!\int_a^b f(x)\,dx=\lim_{\varepsilon\to0^+}\Big(\int_a^{c-\varepsilon}f+\int_{c+\varepsilon}^b f\Big)
> $$
> for a singularity at $c\in(a,b)$.

For example $\int_{-\infty}^\infty x\,dx$ diverges, but its principal value is $0$; and $\operatorname{PV}\int_{-1}^{2}\frac{dx}{x}=\lim_{\varepsilon\to0^+}\big(\ln\varepsilon-0+\ln2-\ln\varepsilon\big)=\ln2$. Principal values appear in the Hilbert transform (single-sideband modulation, Kramers–Kronig relations between the real and imaginary parts of a causal response) and are evaluated systematically with contour integration in [Ch. 18](../part-5-transforms/18-complex-analysis.md).

> [!warning] Common Pitfall — Symmetric limits are not the definition
> $\int_{-\infty}^\infty f\,dx$ is defined by two *independent* limits. Computing $\lim_{R\to\infty}\int_{-R}^{R}f$ gives the principal value, which can be finite even when the integral diverges. Only when the integral is known to converge do the two agree.

### 5.3.6 Two important improper integrals

**The Gaussian integral.**

$$
\int_{-\infty}^{\infty}e^{-x^2}dx=\sqrt\pi .
$$

Convergence follows from Example 5.15(a). The value cannot be found with single-variable antiderivatives; the standard trick squares it into a double integral over the plane and switches to polar coordinates — a computation done in [Ch. 6](../part-2-calculus/06-multivariable-calculus.md). Scaling $x\to x/(\sigma\sqrt2)$ gives the normalisation of the Gaussian density, $\int_{-\infty}^\infty e^{-x^2/2\sigma^2}dx=\sigma\sqrt{2\pi}$ ([Ch. 19](../part-6-probability-statistics/19-probability.md)).

**The Gamma function.**

> [!definition] Gamma function
> For $x>0$,
> $$
> \Gamma(x)=\int_0^\infty t^{x-1}e^{-t}\,dt .
> $$

It converges at $\infty$ because the exponential wins, and at $0$ because $t^{x-1}$ has $p=1-x<1$. Integration by parts (§5.5) gives $\Gamma(x+1)=x\,\Gamma(x)$, and $\Gamma(1)=1$, so

$$
\Gamma(n+1)=n!\qquad(n=0,1,2,\dots),
$$

and $\Gamma(\tfrac12)=\sqrt\pi$ (substitute $t=u^2$ and use the Gaussian integral). The Gamma function makes "$\int_0^\infty t^ne^{-st}dt=n!/s^{n+1}$" — the Laplace transform of $t^n$ — a one-line result.

> [!tip] Engineering Note — Existence of transforms
> The Laplace transform $\int_0^\infty f(t)e^{-st}dt$ converges absolutely when $\lvert f(t)\rvert\le Me^{\alpha t}$ and $\operatorname{Re}s>\alpha$: compare with $\int_0^\infty Me^{-(\operatorname{Re}s-\alpha)t}dt$. That is the entire content of the "exponential order" existence theorem in [Ch. 14](../part-5-transforms/14-laplace-transforms.md) — a direct comparison test.

## 5.4 Multiple integrals

### 5.4.1 Motivation

Charge is rarely concentrated on a line. A semiconductor wafer carries a surface charge density $\sigma(x,y)$ in C/m²; a depletion region holds a volume charge density $\rho(x,y,z)$ in C/m³; a printed coil's inductance depends on how current is distributed across a cross-section. The total charge is obtained exactly as in §5.1 — chop, multiply by the size of each piece, add, take a limit — except that the pieces are now small rectangles of area $\Delta A$ or small boxes of volume $\Delta V$. The result is a **double** or **triple integral**.

This section introduces multiple integrals over simple regions in Cartesian coordinates and shows how to compute them as *iterated* single integrals, so that everything from §§5.1–5.3 carries over. [Chapter 6](../part-2-calculus/06-multivariable-calculus.md) develops the subject further: general changes of variables, the Jacobian, and polar, cylindrical and spherical coordinates.

### 5.4.2 Double integrals over a rectangle

> [!definition] Double integral
> Let $f$ be bounded on the rectangle $R=[a,b]\times[c,d]$. Partition $[a,b]$ and $[c,d]$, giving sub-rectangles $R_{ij}$ of area $\Delta A_{ij}=\Delta x_i\,\Delta y_j$, and pick sample points $(x_{ij}^*,y_{ij}^*)\in R_{ij}$. If the double Riemann sum
> $$
> \sum_{i}\sum_{j}f(x_{ij}^*,y_{ij}^*)\,\Delta A_{ij}
> $$
> tends to a limit $I$ as the largest sub-rectangle diagonal tends to zero, independent of the sample points, then
> $$
> \iint_R f(x,y)\,dA=I .
> $$

If $f\ge0$, the double integral is the **volume** under the surface $z=f(x,y)$ above $R$. If $f$ is a surface density (charge, mass, power per area), it is the **total amount**. With $f\equiv1$ it is the **area** of $R$. Linearity, additivity over regions, monotonicity and the triangle inequality all carry over from §5.1.3. Continuous functions — and bounded functions whose discontinuities lie on finitely many smooth curves — are integrable.

Computing a two-dimensional limit directly is hopeless, but the following theorem reduces it to two one-dimensional integrals.

> [!theorem] Fubini's theorem (rectangle)
> If $f$ is continuous on $R=[a,b]\times[c,d]$, then
> $$
> \iint_R f(x,y)\,dA=\int_a^b\Big[\int_c^d f(x,y)\,dy\Big]dx=\int_c^d\Big[\int_a^b f(x,y)\,dx\Big]dy .
> $$

*Proof idea.* Arrange the double Riemann sum in a grid. Summing first down each column (fixed $x_i$) gives approximately $\big[\int_c^d f(x_i,y)\,dy\big]\Delta x_i$; summing those columns across gives a Riemann sum for the outer $x$-integral. Summing rows first gives the other order. Both are rearrangements of the same finite sum, and continuity ensures both limits exist and agree. $\square$

In the inner integral the outer variable is held **constant**, exactly as in partial differentiation. The bracket is a function of the outer variable alone.

**Separable integrands.** If $f(x,y)=g(x)h(y)$ on a rectangle, the double integral factorises:

$$
\iint_{[a,b]\times[c,d]}g(x)h(y)\,dA=\Big(\int_a^b g(x)\,dx\Big)\Big(\int_c^d h(y)\,dy\Big).
$$

> [!example] Worked Example 5.16 — A rectangle
> Evaluate $\displaystyle\iint_R x^2y\,dA$ over $R=[0,2]\times[1,3]$.
>
> **Order $dy$ then $dx$:**
> $$
> \int_0^2\Big[\int_1^3x^2y\,dy\Big]dx=\int_0^2x^2\Big[\frac{y^2}{2}\Big]_1^3dx=\int_0^24x^2\,dx=\frac{32}{3}.
> $$
> **Check by factorising:** $\big(\int_0^2x^2dx\big)\big(\int_1^3y\,dy\big)=\frac83\cdot4=\frac{32}{3}\approx10.667$ ✓.

### 5.4.3 Double integrals over general regions

For a region $D$ that is not a rectangle, extend $f$ by zero outside $D$ and integrate over an enclosing rectangle. In practice we describe $D$ by its bounding curves.

> [!definition] Type I and Type II regions
> - **Type I (vertically simple):** $D=\{(x,y):a\le x\le b,\ g_1(x)\le y\le g_2(x)\}$. Then
> $$
> \iint_D f\,dA=\int_a^b\int_{g_1(x)}^{g_2(x)}f(x,y)\,dy\,dx .
> $$
> - **Type II (horizontally simple):** $D=\{(x,y):c\le y\le d,\ h_1(y)\le x\le h_2(y)\}$. Then
> $$
> \iint_D f\,dA=\int_c^d\int_{h_1(y)}^{h_2(y)}f(x,y)\,dx\,dy .
> $$

**Recipe for setting limits.** (1) Sketch $D$. (2) Choose the inner variable, say $y$. (3) Draw a vertical line through $D$ at a typical $x$: it enters $D$ on the lower curve $y=g_1(x)$ and leaves on the upper curve $y=g_2(x)$ — those are the inner limits, and they may depend on $x$. (4) The outer limits are the extreme $x$-values of $D$ and must be **constants**. A region that is neither type is split into pieces that are.

```desmos-graph
left=-0.2; right=1.3; bottom=-0.2; top=1.3
---
y=x|blue
y=x^2|red
x^2<y<x|grey|0<x<1
x=0.6|black|dashed
(0.6,0.36)|label:enters at y = x²
(0.6,0.6)|label:leaves at y = x
```
*The region between $y=x^2$ and $y=x$. A vertical line at $x=0.6$ enters on the parabola and leaves on the line, fixing the inner $y$-limits.*

> [!example] Worked Example 5.17 — A region between two curves, both orders
> Evaluate $\displaystyle\iint_D xy\,dA$ where $D$ is bounded by $y=x$ and $y=x^2$.
>
> **Step 1 — sketch.** The curves meet where $x=x^2$, i.e. $x=0$ and $x=1$. For $0<x<1$, $x^2<x$.
>
> **Step 2 — Type I.**
> $$
> \int_0^1\int_{x^2}^{x}xy\,dy\,dx=\int_0^1x\cdot\frac{x^2-x^4}{2}dx=\frac12\Big(\frac14-\frac16\Big)=\frac{1}{24}.
> $$
>
> **Step 3 — Type II, as a check.** A horizontal line at height $y$ enters at $x=y$ and leaves at $x=\sqrt y$:
> $$
> \int_0^1\int_{y}^{\sqrt y}xy\,dx\,dy=\int_0^1y\cdot\frac{y-y^2}{2}dy=\frac12\Big(\frac13-\frac14\Big)=\frac{1}{24}\ \checkmark
> $$
> **Area of $D$** (integrand $1$): $\int_0^1(x-x^2)dx=\frac16$. So the **mean value** of $xy$ over $D$ is $\frac{1/24}{1/6}=\frac14$.

### 5.4.4 Reversing the order of integration

Sometimes the inner integral in one order has no elementary antiderivative, but the other order is easy. Reversing the order means re-describing the *same region* the other way — the limits do not simply swap.

> [!example] Worked Example 5.18 — Reversal rescues an impossible integral
> Evaluate $\displaystyle I=\int_0^1\int_x^1e^{y^2}\,dy\,dx$.
>
> **Step 1 — the obstacle.** $\int e^{y^2}dy$ has no elementary antiderivative.
>
> **Step 2 — describe the region.** $0\le x\le1$, $x\le y\le1$: the triangle with vertices $(0,0)$, $(0,1)$, $(1,1)$, lying above the line $y=x$.
>
> **Step 3 — redescribe as Type II.** For fixed $y\in[0,1]$, $x$ runs from $0$ to $y$.
>
> **Step 4 — integrate.**
> $$
> I=\int_0^1\int_0^{y}e^{y^2}\,dx\,dy=\int_0^1y\,e^{y^2}dy=\Big[\tfrac12e^{y^2}\Big]_0^1=\frac{e-1}{2}\approx0.85914 .
> $$
> The inner integration produced the factor $y$, which is precisely the derivative needed for the substitution $u=y^2$ (§5.6).

> [!warning] Common Pitfall — Variable outer limits
> The outermost limits must be constants. Writing $\int_x^1\int_0^1 f\,dx\,dy$ is meaningless — after both integrations an $x$ would survive. When reversing order, *always* sketch the region and read the new limits off the sketch.

### 5.4.5 Triple integrals

The triple integral $\iiint_E f(x,y,z)\,dV$ is defined by Riemann sums over small boxes of volume $\Delta V=\Delta x\,\Delta y\,\Delta z$. With $f\equiv1$ it gives the volume of $E$; with $f=\rho$ (charge density) it gives the enclosed charge — the right-hand side of Gauss's law in integral form ([Ch. 7](../part-2-calculus/07-vector-calculus.md)). Fubini's theorem extends: for a solid described by

$$
E=\{(x,y,z):a\le x\le b,\ g_1(x)\le y\le g_2(x),\ u_1(x,y)\le z\le u_2(x,y)\},
$$

$$
\iiint_E f\,dV=\int_a^b\int_{g_1(x)}^{g_2(x)}\int_{u_1(x,y)}^{u_2(x,y)}f(x,y,z)\,dz\,dy\,dx .
$$

Work from the inside out: integrate in $z$ treating $x,y$ as constants, then in $y$ treating $x$ as constant, then in $x$. The limits at each stage may depend only on variables that are still to be integrated.

> [!example] Worked Example 5.19 — Charge in a tetrahedron
> A charge density $\rho=\rho_0\,z/h$ (with $\rho_0$ in C/m³) fills the tetrahedron bounded by the coordinate planes and the plane $x+y+z=h$. Find the total charge.
>
> **Step 1 — limits.** The base in the $xy$-plane is the triangle $x\ge0$, $y\ge0$, $x+y\le h$. Above a point $(x,y)$, $z$ runs from $0$ to $h-x-y$.
>
> **Step 2 — innermost.** $\displaystyle\int_0^{h-x-y}\frac{\rho_0 z}{h}dz=\frac{\rho_0}{2h}(h-x-y)^2$.
>
> **Step 3 — middle** ($y$ from $0$ to $h-x$). With $w=h-x$: $\displaystyle\int_0^{w}\frac{\rho_0}{2h}(w-y)^2dy=\frac{\rho_0}{2h}\cdot\frac{w^3}{3}$.
>
> **Step 4 — outer.** $\displaystyle\int_0^h\frac{\rho_0}{6h}(h-x)^3dx=\frac{\rho_0}{6h}\cdot\frac{h^4}{4}=\frac{\rho_0h^3}{24}.$
>
> **Check.** The volume is $h^3/6$, so the mean density is $\rho_0/4$, i.e. the mean of $z/h$ over the tetrahedron is $1/4$ — consistent with the centroid of a tetrahedron lying a quarter of the way up from the base. Units: $[\text{C/m}^3][\text{m}^3]=[\text{C}]$ ✓. With $h=1$, a Monte Carlo estimate of $\iiint z\,dV$ gives $0.0417\approx1/24$.

### 5.4.6 A first look at polar coordinates

Circular regions are awkward in Cartesian form: the disc $x^2+y^2\le a^2$ needs limits $-\sqrt{a^2-x^2}\le y\le\sqrt{a^2-x^2}$. In polar coordinates $(r,\phi)$ it is simply $0\le r\le a$, $0\le\phi\le2\pi$. The area element is **not** $dr\,d\phi$: a small polar "rectangle" has radial side $dr$ and arc side $r\,d\phi$, so

$$
dA=r\,dr\,d\phi .
$$

The extra factor $r$ is the simplest example of a **Jacobian**; [Chapter 6](../part-2-calculus/06-multivariable-calculus.md) derives it in general and extends it to cylindrical ($dV=r\,dr\,d\phi\,dz$) and spherical ($dV=R^2\sin\theta\,dR\,d\theta\,d\phi$) coordinates.

> [!example] Worked Example 5.20 — Charge on a disc
> A disc of radius $a$ carries surface charge density $\sigma(r)=\sigma_0\,r/a$. Find the total charge.
> $$
> Q=\int_0^{2\pi}\int_0^a\sigma_0\frac{r}{a}\,r\,dr\,d\phi=2\pi\cdot\frac{\sigma_0}{a}\cdot\frac{a^3}{3}=\frac{2\pi\sigma_0a^2}{3}.
> $$
> **Check.** A uniform density $\sigma_0$ would give $\pi a^2\sigma_0$; our density averages to $\tfrac23\sigma_0$ over the disc (more area sits near the rim, where $\sigma$ is largest), giving $\tfrac23\pi a^2\sigma_0$ ✓. Forgetting the factor $r$ in $dA$ would give $\pi\sigma_0a$ — dimensionally wrong (C/m instead of C), which is how the error is caught.

> [!tip] Engineering Note — Choosing coordinates
> The effort of a multiple integral is set mostly by how well the coordinates fit the boundary. Parallel-plate geometries suit Cartesian coordinates; coaxial cables, wires and solenoids suit cylindrical; point charges and antennas suit spherical. Picking coordinates in which every boundary is "coordinate = constant" turns every limit into a constant and often makes the integrand separable.

## 5.5 Integration by parts

### 5.5.1 The rule and where it comes from

Substitution (§5.6) undoes the chain rule. **Integration by parts** undoes the product rule, and it is the tool for integrands that are products of *unlike* functions: $t\,e^{-t/\tau}$, $t\sin\omega t$, $\ln x$, $x^2\cos x$. It is also the mechanism behind reduction formulas, the Gamma-function recursion, the derivative rule of the Laplace transform, and — in higher dimensions — the divergence theorem and Green's identities ([Ch. 7](../part-2-calculus/07-vector-calculus.md)).

Start from the product rule $(uv)'=u'v+uv'$ and integrate both sides:

$$
uv=\int u'v\,dx+\int uv'\,dx .
$$

Rearranging gives:

> [!theorem] Integration by parts
> If $u$ and $v$ have continuous derivatives,
> $$
> \int u\,dv=uv-\int v\,du,\qquad\text{i.e.}\qquad \int u(x)v'(x)\,dx=u(x)v(x)-\int u'(x)v(x)\,dx ,
> $$
> and for definite integrals
> $$
> \int_a^b u\,dv=\Big[uv\Big]_a^b-\int_a^b v\,du .
> $$

The rule *trades* one integral for another. It helps when $\int v\,du$ is simpler than $\int u\,dv$, which happens when differentiating $u$ simplifies it (a polynomial loses a degree; $\ln x$ becomes $1/x$) while integrating $dv$ does not make it worse (exponentials and sinusoids stay the same kind).

**Choosing $u$.** A long-standing heuristic is **LIATE**: choose $u$ as the first type present in the list

> **L**ogarithmic, **I**nverse trigonometric, **A**lgebraic (polynomials), **T**rigonometric, **E**xponential,

and let $dv$ be the rest. Functions early in the list get simpler when differentiated; those late in the list are easy to integrate.

> [!example] Worked Example 5.21 — The basic pattern
> (a) $\displaystyle\int x e^{x}dx$. Take $u=x$, $dv=e^xdx$, so $du=dx$, $v=e^x$:
> $$
> \int xe^xdx=xe^x-\int e^xdx=(x-1)e^x+C.
> $$
> **Check:** $\frac{d}{dx}(x-1)e^x=e^x+(x-1)e^x=xe^x$ ✓.
>
> (b) $\displaystyle\int x\cos3x\,dx$. $u=x$, $dv=\cos3x\,dx$, $v=\tfrac13\sin3x$:
> $$
> \int x\cos3x\,dx=\frac{x\sin3x}{3}-\frac13\int\sin3x\,dx=\frac{x\sin3x}{3}+\frac{\cos3x}{9}+C.
> $$
>
> (c) $\displaystyle\int\ln x\,dx$. There is no obvious product — take $u=\ln x$, $dv=dx$, so $du=dx/x$, $v=x$:
> $$
> \int\ln x\,dx=x\ln x-\int x\cdot\frac1x\,dx=x\ln x-x+C.
> $$
>
> (d) $\displaystyle\int\arctan x\,dx$. $u=\arctan x$, $dv=dx$:
> $$
> \int\arctan x\,dx=x\arctan x-\int\frac{x}{1+x^2}dx=x\arctan x-\tfrac12\ln(1+x^2)+C.
> $$
> Parts (c) and (d) show the trick of taking $dv=dx$: it converts the integral of a function whose *derivative* is algebraic into an algebraic integral.

> [!example] Worked Example 5.22 — Mean delay of an RC response
> The impulse response of an RC low-pass filter is $h(t)=\frac{1}{\tau}e^{-t/\tau}$, $t\ge0$, with $\int_0^\infty h\,dt=1$. Its **mean delay** (first moment, the basis of the Elmore delay used to estimate interconnect delay in integrated circuits) is $\bar t=\int_0^\infty t\,h(t)\,dt$.
>
> **Step 1 — parts.** $u=t$, $dv=\frac1\tau e^{-t/\tau}dt$, $v=-e^{-t/\tau}$:
> $$
> \int_0^\infty\frac{t}{\tau}e^{-t/\tau}dt=\Big[-te^{-t/\tau}\Big]_0^\infty+\int_0^\infty e^{-t/\tau}dt .
> $$
> **Step 2 — boundary term.** $te^{-t/\tau}\to0$ as $t\to\infty$ (exponential beats polynomial) and is $0$ at $t=0$.
>
> **Step 3 — remaining integral.** $\int_0^\infty e^{-t/\tau}dt=\tau$ (Example 5.12a). Hence $\bar t=\tau$.
>
> The mean delay equals the time constant. With $R=1$ kΩ and $C=1$ nF, $\bar t=1\ \mu$s.

### 5.5.2 Repeated parts and the tabular method

When $u$ is a polynomial of degree $n$, parts must be applied $n$ times. The bookkeeping is compressed by the **tabular method**: list $u$ and its successive derivatives in one column (until zero), $dv$ and its successive integrals in another, and combine diagonally with alternating signs $+,-,+,\dots$

> [!example] Worked Example 5.23 — Tabular integration
> $\displaystyle\int x^3e^{2x}dx$.
>
> | sign | $D$: derivatives of $x^3$ | $I$: integrals of $e^{2x}$ |
> |---|---|---|
> | $+$ | $x^3$ | $e^{2x}$ |
> | $-$ | $3x^2$ | $\tfrac12e^{2x}$ |
> | $+$ | $6x$ | $\tfrac14e^{2x}$ |
> | $-$ | $6$ | $\tfrac18e^{2x}$ |
> | | $0$ | $\tfrac1{16}e^{2x}$ |
>
> Multiply each $D$ entry by the $I$ entry one row lower, with the sign of its row:
> $$
> \int x^3e^{2x}dx=e^{2x}\Big(\frac{x^3}{2}-\frac{3x^2}{4}+\frac{6x}{8}-\frac{6}{16}\Big)+C=e^{2x}\Big(\frac{x^3}{2}-\frac{3x^2}{4}+\frac{3x}{4}-\frac38\Big)+C.
> $$
> **Check** with $\int_0^1$: formula gives $e^2(0.5-0.75+0.75-0.375)+0.375=0.125e^2+0.375=1.29863$; numerical quadrature gives $1.298632$ ✓.

### 5.5.3 Cyclic integrals

For $\int e^{ax}\sin bx\,dx$ two applications of parts return the original integral, which can then be solved for algebraically.

> [!example] Worked Example 5.24 — Solving for the integral
> $I=\displaystyle\int e^{x}\sin x\,dx$.
>
> First parts ($u=\sin x$, $dv=e^xdx$): $I=e^x\sin x-\int e^x\cos x\,dx$.
>
> Second parts ($u=\cos x$, $dv=e^xdx$): $\int e^x\cos x\,dx=e^x\cos x+\int e^x\sin x\,dx=e^x\cos x+I$.
>
> So $I=e^x\sin x-e^x\cos x-I$, giving
> $$
> I=\frac{e^x(\sin x-\cos x)}{2}+C,
> $$
> which matches the complex-exponential formula of §5.2.6 with $a=b=1$. Keep the *same type* of choice for $u$ (trig both times, or exponential both times); switching sides undoes the first step and yields the useless identity $I=I$.

### 5.5.4 Reduction formulas

Parts can lower an exponent in a whole family of integrals at once.

> [!theorem] Reduction formula for powers of sine
> For $n\ge2$,
> $$
> \int\sin^nx\,dx=-\frac{\sin^{n-1}x\cos x}{n}+\frac{n-1}{n}\int\sin^{n-2}x\,dx ,
> $$
> and on $[0,\pi/2]$, where the boundary term vanishes, $I_n=\displaystyle\int_0^{\pi/2}\sin^nx\,dx$ satisfies $I_n=\dfrac{n-1}{n}I_{n-2}$ with $I_0=\pi/2$, $I_1=1$.

*Derivation.* Write $\sin^nx=\sin^{n-1}x\cdot\sin x$, take $u=\sin^{n-1}x$, $dv=\sin x\,dx$. Then $du=(n-1)\sin^{n-2}x\cos x\,dx$, $v=-\cos x$, and

$$
\int\sin^nx\,dx=-\sin^{n-1}x\cos x+(n-1)\int\sin^{n-2}x\cos^2x\,dx .
$$

Replace $\cos^2x=1-\sin^2x$: the right side contains $(n-1)\int\sin^{n-2}x\,dx-(n-1)\int\sin^nx\,dx$. Move the last term to the left, collect $n\int\sin^nx\,dx$, and divide by $n$. $\square$

The resulting **Wallis formulas** are

$$
I_{2m}=\frac{(2m-1)(2m-3)\cdots1}{2m(2m-2)\cdots2}\cdot\frac\pi2,\qquad I_{2m+1}=\frac{2m(2m-2)\cdots2}{(2m+1)(2m-1)\cdots3}.
$$

For instance $\int_0^{\pi/2}\sin^4x\,dx=\frac34\cdot\frac12\cdot\frac\pi2=\frac{3\pi}{16}\approx0.58905$. The identical formulas hold for $\cos^nx$ on $[0,\pi/2]$.

The same idea applied to $\Gamma(x+1)=\int_0^\infty t^xe^{-t}dt$ with $u=t^x$, $dv=e^{-t}dt$ gives

$$
\Gamma(x+1)=\Big[-t^xe^{-t}\Big]_0^\infty+x\int_0^\infty t^{x-1}e^{-t}dt=x\,\Gamma(x),
$$

the recursion quoted in §5.3.6.

### 5.5.5 Fourier coefficients by parts

Fourier coefficients of polynomial-shaped waveforms (ramps, sawtooth and triangle waves) are always integrals of polynomial × sinusoid, i.e. integration by parts.

> [!example] Worked Example 5.25 — Sawtooth coefficients
> The sawtooth $f(t)=t$ on $(-\pi,\pi)$, extended with period $2\pi$, is odd, so only sine terms appear, with
> $$
> b_n=\frac1\pi\int_{-\pi}^{\pi}t\sin nt\,dt=\frac2\pi\int_0^{\pi}t\sin nt\,dt .
> $$
> **Parts** with $u=t$, $dv=\sin nt\,dt$, $v=-\frac{\cos nt}{n}$:
> $$
> \int_0^\pi t\sin nt\,dt=\Big[-\frac{t\cos nt}{n}\Big]_0^\pi+\frac1n\int_0^\pi\cos nt\,dt=-\frac{\pi\cos n\pi}{n}+0=\frac{\pi(-1)^{n+1}}{n}.
> $$
> Hence $b_n=\dfrac{2(-1)^{n+1}}{n}$: $b_1=2$, $b_2=-1$, $b_3=\tfrac23,\dots$ The $1/n$ decay reflects the jump discontinuity of the sawtooth; the triangle wave, which is continuous with a kinked derivative, has coefficients decaying like $1/n^2$ — one more integration by parts gains one more power of $n$. This link between smoothness and spectral decay is developed in [Ch. 15](../part-5-transforms/15-fourier-series.md).

> [!warning] Common Pitfall — Dropping or mishandling the boundary term
> In definite integrals the term $[uv]_a^b$ is easy to forget, and for improper integrals it must be evaluated as a *limit*. In Example 5.22, $[-te^{-t/\tau}]_0^\infty=0$ needed the fact that $te^{-t/\tau}\to0$. Had the integrand been $t\,e^{+t}$ the boundary term would diverge and the whole computation would be meaningless.

> [!tip] Engineering Note — Parts in transform theory
> The Laplace-transform rule $\mathcal{L}\{f'(t)\}=sF(s)-f(0)$ is a single integration by parts:
> $$
> \int_0^\infty f'(t)e^{-st}dt=\Big[f(t)e^{-st}\Big]_0^\infty+s\int_0^\infty f(t)e^{-st}dt=-f(0)+sF(s).
> $$
> The boundary term $-f(0)$ is how initial conditions — initial capacitor voltages and inductor currents — enter circuit analysis in the $s$-domain ([Ch. 14](../part-5-transforms/14-laplace-transforms.md)).

## 5.6 Integration by substitution

### 5.6.1 The rule

The chain rule says $\frac{d}{dx}F(g(x))=F'(g(x))\,g'(x)$. Reading it backwards: if $F'=f$, then $\int f(g(x))\,g'(x)\,dx=F(g(x))+C$. Writing $u=g(x)$ and $du=g'(x)\,dx$ makes this look like an algebraic change of symbols.

> [!theorem] Substitution rule
> If $u=g(x)$ has a continuous derivative and $f$ is continuous on the range of $g$,
> $$
> \int f(g(x))\,g'(x)\,dx=\int f(u)\,du\Big|_{u=g(x)} ,
> $$
> and for definite integrals
> $$
> \int_a^b f(g(x))\,g'(x)\,dx=\int_{g(a)}^{g(b)}f(u)\,du .
> $$

*Proof of the definite form.* Let $F'=f$. Both sides equal $F(g(b))-F(g(a))$: the left by the FTC applied to $(F\circ g)'=f(g)g'$, the right by the FTC applied to $F$. $\square$

Note that in the definite form **the limits change** and there is no need to convert back to $x$. Dimensionally, $du$ carries the units of $u$; a correct substitution keeps $f(u)\,du$ dimensionally equal to the original integrand times $dx$.

**How to spot it.** Look for a function *inside* another function whose derivative (up to a constant factor) also appears as a factor. Typical patterns:

| Pattern | Substitution | Result |
|---|---|---|
| $\int f(ax+b)\,dx$ | $u=ax+b$ | $\frac1aF(ax+b)$ |
| $\int \frac{g'(x)}{g(x)}dx$ | $u=g(x)$ | $\ln\lvert g(x)\rvert$ |
| $\int g(x)^n g'(x)\,dx$ | $u=g(x)$ | $\frac{g^{n+1}}{n+1}$, $n\ne-1$ |
| $\int e^{g(x)}g'(x)\,dx$ | $u=g(x)$ | $e^{g(x)}$ |
| $\int x\,f(x^2)\,dx$ | $u=x^2$ | $\frac12F(x^2)$ |
| $\int \frac{f(\ln x)}{x}dx$ | $u=\ln x$ | $F(\ln x)$ |

> [!example] Worked Example 5.26 — Recognising the inner function
> (a) $\displaystyle\int2x\cos(x^2)\,dx$: $u=x^2$, $du=2x\,dx$, giving $\int\cos u\,du=\sin(x^2)+C$.
>
> (b) $\displaystyle\int\frac{x}{(x^2+a^2)^{3/2}}dx$: $u=x^2+a^2$, $du=2x\,dx$:
> $$
> \frac12\int u^{-3/2}du=-u^{-1/2}+C=-\frac{1}{\sqrt{x^2+a^2}}+C .
> $$
> This integral gives the field on the axis of a charged ring or disc.
>
> (c) $\displaystyle\int\tan x\,dx=\int\frac{\sin x}{\cos x}dx$: $u=\cos x$, $du=-\sin x\,dx$, giving $-\ln\lvert\cos x\rvert+C$.
>
> (d) $\displaystyle\int\frac{e^{x}}{1+e^{2x}}dx$: $u=e^x$, $du=e^x dx$, giving $\int\frac{du}{1+u^2}=\arctan(e^x)+C$.

> [!example] Worked Example 5.27 — Definite integrals with changed limits
> (a) $\displaystyle\int_0^1x\sqrt{1-x^2}\,dx$. Let $u=1-x^2$, $du=-2x\,dx$; $x=0\mapsto u=1$, $x=1\mapsto u=0$:
> $$
> \int_1^0\sqrt u\,\Big(-\frac{du}{2}\Big)=\frac12\int_0^1u^{1/2}du=\frac12\cdot\frac23=\frac13 .
> $$
>
> (b) $\displaystyle\int_1^{e}\frac{(\ln x)^2}{x}dx$. Let $u=\ln x$, $du=dx/x$; limits $0$ to $1$: $\int_0^1u^2du=\frac13$.
>
> (c) $\displaystyle\int_0^{\infty}t\,e^{-t^2/2\sigma^2}dt$ (a Rayleigh-type integral from fading-channel statistics). Let $u=t^2/2\sigma^2$, $du=t\,dt/\sigma^2$:
> $$
> \sigma^2\int_0^\infty e^{-u}du=\sigma^2 .
> $$

> [!warning] Common Pitfall — Forgetting to change the limits (or the $dx$)
> In a definite integral either change the limits to $u$-values *or* convert the antiderivative back to $x$ before substituting the original limits — never mix the two. Equally, every $x$ must disappear: if after writing $du=g'(x)dx$ some stray $x$ remains that cannot be expressed in terms of $u$, the substitution was the wrong one.

### 5.6.2 Linear substitutions: shifting and scaling

The simplest substitutions, $u=t-t_0$ and $u=at$, are used so often in signal analysis that their effects should be automatic:

$$
\int_{-\infty}^{\infty}f(t-t_0)\,dt=\int_{-\infty}^{\infty}f(u)\,du,\qquad
\int_{-\infty}^{\infty}f(at)\,dt=\frac{1}{\lvert a\rvert}\int_{-\infty}^{\infty}f(u)\,du\quad(a\ne0).
$$

A time shift does not change the area of a pulse; compressing time by a factor $a$ divides the area by $\lvert a\rvert$. (For $a<0$ the limits swap, which produces the absolute value.) The same scaling argument shows why the Dirac impulse obeys $\delta(at)=\delta(t)/\lvert a\rvert$, and why compressing a signal in time spreads its spectrum ([Ch. 16](../part-5-transforms/16-fourier-transforms.md)).

For a $T$-periodic function, the substitution $u=t-c$ proves that the integral over any full period is the same (§5.1.3), and $u=\omega t$ converts integrals over one period to the standard interval $[0,2\pi]$:

$$
\frac1T\int_0^T g(\omega t)\,dt=\frac{1}{2\pi}\int_0^{2\pi}g(\theta)\,d\theta,\qquad \omega=\frac{2\pi}{T}.
$$

This is why averages of sinusoidal waveforms do not depend on frequency.

### 5.6.3 Symmetry by substitution

The substitution $x\to-x$ proves the even/odd rules of §5.1.3: for odd $f$, $\int_{-a}^0f(x)dx=\int_a^0f(-u)(-du)=-\int_0^af(u)du$, so the halves cancel. A less obvious use is the "king's rule" $x\to a+b-x$:

$$
\int_a^b f(x)\,dx=\int_a^b f(a+b-x)\,dx ,
$$

which evaluates, for example, $\int_0^{\pi/2}\frac{\sin x}{\sin x+\cos x}dx$: adding the integral to its transformed self gives $\int_0^{\pi/2}1\,dx=\frac\pi2$, so the integral is $\frac\pi4$.

### 5.6.4 Trigonometric and hyperbolic substitution

When the integrand contains $\sqrt{a^2-x^2}$, $\sqrt{a^2+x^2}$ or $\sqrt{x^2-a^2}$, an **inverse substitution** $x=g(\theta)$ chosen to turn the radicand into a perfect square removes the root. Here the new variable is introduced through $x$, and $dx=g'(\theta)\,d\theta$.

| Radical | Substitution | Identity used | $dx$ |
|---|---|---|---|
| $\sqrt{a^2-x^2}$ | $x=a\sin\theta$, $\lvert\theta\rvert\le\frac\pi2$ | $1-\sin^2=\cos^2$ | $a\cos\theta\,d\theta$ |
| $\sqrt{a^2+x^2}$ | $x=a\tan\theta$, $\lvert\theta\rvert<\frac\pi2$ | $1+\tan^2=\sec^2$ | $a\sec^2\theta\,d\theta$ |
| $\sqrt{a^2+x^2}$ | $x=a\sinh u$ | $1+\sinh^2=\cosh^2$ | $a\cosh u\,du$ |
| $\sqrt{x^2-a^2}$ | $x=a\sec\theta$ or $x=a\cosh u$ | $\sec^2-1=\tan^2$, $\cosh^2-1=\sinh^2$ | |

The restricted ranges of $\theta$ make the substitution invertible and fix the sign of the square root (e.g. $\cos\theta\ge0$ for $\lvert\theta\rvert\le\pi/2$).

> [!example] Worked Example 5.28 — Area of a circle
> Show that $\displaystyle\int_{-a}^{a}\sqrt{a^2-x^2}\,dx=\frac{\pi a^2}{2}$.
>
> Let $x=a\sin\theta$; $x=\pm a\mapsto\theta=\pm\frac\pi2$; $\sqrt{a^2-x^2}=a\cos\theta$; $dx=a\cos\theta\,d\theta$:
> $$
> \int_{-\pi/2}^{\pi/2}a^2\cos^2\theta\,d\theta=a^2\int_{-\pi/2}^{\pi/2}\frac{1+\cos2\theta}{2}d\theta=a^2\cdot\frac\pi2 .
> $$
> The area of a half-disc, as it must be.

> [!example] Worked Example 5.29 — Field of a finite line charge
> A uniform line charge $\lambda$ (C/m) lies on the $z$-axis from $z=-L$ to $z=L$. The radial field at distance $r$ in the mid-plane is
> $$
> E_r=\frac{\lambda r}{4\pi\varepsilon_0}\int_{-L}^{L}\frac{dz}{(r^2+z^2)^{3/2}} .
> $$
>
> **Step 1 — substitute** $z=r\tan\theta$, $dz=r\sec^2\theta\,d\theta$, $r^2+z^2=r^2\sec^2\theta$:
> $$
> \int\frac{r\sec^2\theta\,d\theta}{r^3\sec^3\theta}=\frac{1}{r^2}\int\cos\theta\,d\theta=\frac{\sin\theta}{r^2}.
> $$
>
> **Step 2 — back-substitute.** From the triangle with opposite side $z$ and adjacent side $r$: $\sin\theta=\dfrac{z}{\sqrt{r^2+z^2}}$. So
> $$
> \int_{-L}^{L}\frac{dz}{(r^2+z^2)^{3/2}}=\frac{2L}{r^2\sqrt{r^2+L^2}} ,\qquad E_r=\frac{\lambda}{4\pi\varepsilon_0}\cdot\frac{2L}{r\sqrt{r^2+L^2}} .
> $$
>
> **Step 3 — limits.** As $L\to\infty$ (an improper integral), $E_r\to\dfrac{\lambda}{2\pi\varepsilon_0r}$, the infinite-line result. For $L\ll r$, $E_r\approx\dfrac{2L\lambda}{4\pi\varepsilon_0r^2}$ — the field of a point charge $Q=2L\lambda$, as it should be.
>
> **Check:** with $r=1$, $L=2$, the integral is $4/\sqrt5=1.78885$; quadrature gives $1.788854$ ✓.

> [!example] Worked Example 5.30 — A hyperbolic substitution
> $\displaystyle\int\frac{dx}{\sqrt{x^2+a^2}}$. Let $x=a\sinh u$, $dx=a\cosh u\,du$, $\sqrt{x^2+a^2}=a\cosh u$:
> $$
> \int du=u+C=\operatorname{arsinh}\frac{x}{a}+C=\ln\big(x+\sqrt{x^2+a^2}\big)+C',
> $$
> the table entry of §5.1.4. The hyperbolic route avoids the awkward $\int\sec\theta\,d\theta$ that the tangent substitution produces. The same integral gives the capacitance per length of a two-wire line, $C'=\pi\varepsilon/\operatorname{arcosh}(D/d)$, where an inverse hyperbolic function appears for exactly this reason.

### 5.6.5 The tangent half-angle substitution

Any rational function of $\sin\theta$ and $\cos\theta$ becomes a rational function of $t$ under

$$
t=\tan\frac\theta2,\qquad \sin\theta=\frac{2t}{1+t^2},\qquad \cos\theta=\frac{1-t^2}{1+t^2},\qquad d\theta=\frac{2\,dt}{1+t^2},
$$

after which partial fractions (§5.2.5) finish the job. It is a last resort — often clumsy — but it always works.

> [!example] Worked Example 5.31 — An integral over a full period
> Evaluate $I=\displaystyle\int_0^{2\pi}\frac{d\theta}{2+\cos\theta}$.
>
> **Step 1 — avoid the singular point.** $t=\tan(\theta/2)$ blows up at $\theta=\pi$, so use symmetry ($\cos$ is even about $\pi$): $I=2\int_0^\pi\frac{d\theta}{2+\cos\theta}$, and $\theta:0\to\pi$ maps to $t:0\to\infty$.
>
> **Step 2 — transform.** $2+\cos\theta=\dfrac{2(1+t^2)+1-t^2}{1+t^2}=\dfrac{3+t^2}{1+t^2}$, so
> $$
> \frac{d\theta}{2+\cos\theta}=\frac{1+t^2}{3+t^2}\cdot\frac{2\,dt}{1+t^2}=\frac{2\,dt}{3+t^2}.
> $$
>
> **Step 3 — integrate.**
> $$
> I=2\int_0^\infty\frac{2\,dt}{t^2+3}=4\cdot\frac{1}{\sqrt3}\cdot\frac\pi2=\frac{2\pi}{\sqrt3}\approx3.6276 .
> $$
> **Check:** a 64-point trapezoidal sum over the period gives $3.627599$ — and for smooth periodic integrands the trapezoidal rule is spectacularly accurate (§5.7.5). Integrals of this form are evaluated in one line by residues in [Ch. 18](../part-5-transforms/18-complex-analysis.md).

> [!tip] Engineering Note — Substitution as change of physical variable
> Engineers substitute all the time without calling it that: converting an integral over time to one over phase ($\theta=\omega t$), over frequency to one over normalised frequency ($u=\omega/\omega_c$), or over distance to one over angle as in Example 5.29. Normalising first — making the variable dimensionless — exposes the shape of the integral and leaves the physical constants outside as a prefactor whose units are easy to check.

## 5.7 Numerical integration

### 5.7.1 Why numerical integration

Three situations force us to integrate numerically:

1. **No elementary antiderivative exists** — $\int_0^1e^{-x^2}dx$, the bit-error probability of a digital link, the inductance of a real coil.
2. **The integrand is known only at samples** — an oscilloscope or ADC record of current, voltage or temperature; the energy of a measured pulse is $\sum$ of samples, not a formula.
3. **A closed form exists but is not worth deriving** — when only a number is required, a few lines of code beat a page of algebra.

All practical methods are **quadrature rules** of the form

$$
\int_a^b f(x)\,dx\approx\sum_{k=0}^{n}w_k\,f(x_k),
$$

a weighted sum of samples at **nodes** $x_k$ with **weights** $w_k$. This section develops the three basic Newton–Cotes rules on equally spaced nodes and their error behaviour. [Chapter 22](../part-7-numerical-computational/22-numerical-methods.md) continues with Romberg integration, Gaussian quadrature, adaptive schemes and multidimensional methods.

Throughout, $h=(b-a)/n$ and $x_k=a+kh$, $f_k=f(x_k)$.

### 5.7.2 The midpoint and trapezoidal rules

The **midpoint rule** replaces $f$ on each subinterval by its value at the centre; the **trapezoidal rule** replaces it by the straight line through the end values:

$$
M_n=h\sum_{k=1}^{n}f\Big(x_{k-1}+\frac h2\Big),\qquad
T_n=h\Big[\tfrac12f_0+f_1+f_2+\dots+f_{n-1}+\tfrac12f_n\Big].
$$

Each trapezoid has area $\frac h2(f_{k-1}+f_k)$; summing them gives $T_n$, in which interior samples appear twice and therefore get weight $h$.

> [!theorem] Error of the midpoint and trapezoidal rules
> If $f''$ is continuous on $[a,b]$, there exist $\xi,\eta\in(a,b)$ with
> $$
> \int_a^b f\,dx-T_n=-\frac{(b-a)h^2}{12}f''(\xi),\qquad
> \int_a^b f\,dx-M_n=+\frac{(b-a)h^2}{24}f''(\eta).
> $$
> Both rules are **second order**: halving $h$ divides the error by about $4$. Both are exact for straight lines.

*Derivation sketch (single panel).* Expand $f$ about the panel midpoint $m$ with $s=x-m$: $f=f(m)+sf'(m)+\tfrac12s^2f''(m)+\dots$. Integrating over $s\in[-h/2,h/2]$, the odd term vanishes:

$$
\int f\,dx=hf(m)+\frac{h^3}{24}f''(m)+O(h^5).
$$

So the midpoint rule has local error $+\frac{h^3}{24}f''$. For the trapezoid, average the endpoint expansions $f(m\pm h/2)=f(m)\pm\frac h2f'+\frac{h^2}{8}f''+\dots$: $\frac h2(f_{\text{left}}+f_{\text{right}})=hf(m)+\frac{h^3}{8}f''+\dots$, so its local error is $\frac{h^3}{24}f''-\frac{h^3}{8}f''=-\frac{h^3}{12}f''$. Summing $n=(b-a)/h$ panels multiplies the local error by $n$, giving global error $\propto(b-a)h^2$. $\square$

Two practical consequences follow. First, the midpoint error is about half the trapezoidal error and of **opposite sign**, so for a convex function ($f''>0$) the trapezoid overestimates and the midpoint underestimates, bracketing the true value. Second, the weighted combination $\tfrac13(2M_n+T_n)$ cancels the leading error — that combination is Simpson's rule.

```desmos-graph
left=-0.1; right=1.1; bottom=-0.1; top=1.15
---
y=e^{-x^2}|blue
y=1+(0.77880078-1)(x/0.5)|red|0<x<0.5
y=0.77880078+(0.36787944-0.77880078)((x-0.5)/0.5)|red|0.5<x<1
x=0.5|grey|dotted
x=1|grey|dotted
```
*Trapezoidal rule with $n=2$ for $e^{-x^2}$ on $[0,1]$: the chords (red) lie below the concave part of the curve and above the convex tail, so the errors partly cancel.*

### 5.7.3 Simpson's rule

Fitting a **parabola** through three consecutive samples and integrating it exactly gives, over a pair of panels,

$$
\int_{x_0}^{x_2}f\,dx\approx\frac h3\big(f_0+4f_1+f_2\big).
$$

*Derivation.* Put $x_1=0$, so the nodes are $-h,0,h$. For $p(x)=\alpha+\beta x+\gamma x^2$, $\int_{-h}^{h}p\,dx=2\alpha h+\frac23\gamma h^3$. From the samples, $\alpha=f_1$ and $f_0+f_2=2\alpha+2\gamma h^2$, so $\gamma h^2=\frac{f_0+f_2}{2}-f_1$. Substituting: $2hf_1+\frac{2h}{3}\big(\frac{f_0+f_2}{2}-f_1\big)=\frac h3(f_0+4f_1+f_2)$. $\square$

Applying this to consecutive pairs requires **$n$ even** and gives the composite rule

$$
S_n=\frac h3\Big[f_0+4f_1+2f_2+4f_3+2f_4+\dots+2f_{n-2}+4f_{n-1}+f_n\Big].
$$

> [!theorem] Error of Simpson's rule
> If $f^{(4)}$ is continuous on $[a,b]$ and $n$ is even,
> $$
> \int_a^b f\,dx-S_n=-\frac{(b-a)h^4}{180}f^{(4)}(\xi)
> $$
> for some $\xi\in(a,b)$. Simpson's rule is **fourth order** — halving $h$ divides the error by about $16$ — and, surprisingly, exact for **cubics**, not just quadratics (the cubic term is odd about the panel centre and integrates to zero).

> [!example] Worked Example 5.32 — Comparing the rules on $\int_0^1e^{-x^2}dx$
> The exact value is $\frac{\sqrt\pi}{2}\operatorname{erf}(1)=0.746824133$. With $n=4$, $h=0.25$, the samples are
>
> | $x$ | $0$ | $0.25$ | $0.5$ | $0.75$ | $1$ |
> |---|---|---|---|---|---|
> | $f$ | $1.000000$ | $0.939413$ | $0.778801$ | $0.569783$ | $0.367879$ |
>
> **Trapezoidal:** $T_4=0.25\big[\tfrac12(1+0.367879)+0.939413+0.778801+0.569783\big]=0.742984$.
>
> **Simpson:** $S_4=\frac{0.25}{3}\big[1+4(0.939413)+2(0.778801)+4(0.569783)+0.367879\big]=0.746855$.
>
> **Midpoint** (samples at $0.125,0.375,0.625,0.875$): $M_4=0.748747$.
>
> | Rule | $n=4$ | error | $n=8$ | error | ratio |
> |---|---|---|---|---|---|
> | Trapezoidal | $0.742984$ | $-3.84\times10^{-3}$ | $0.745866$ | $-9.59\times10^{-4}$ | $4.0$ |
> | Midpoint | $0.748747$ | $+1.92\times10^{-3}$ | $0.747304$ | $+4.79\times10^{-4}$ | $4.0$ |
> | Simpson | $0.746855$ | $+3.12\times10^{-5}$ | $0.746826$ | $+1.99\times10^{-6}$ | $15.7$ |
>
> The error ratios on halving $h$ — $4$ and $16$ — confirm the orders $h^2$ and $h^4$. Simpson with five samples beats the trapezoidal rule with nine by a factor of about $30$. Also note $M$ and $T$ bracket the true value with $\lvert E_M\rvert\approx\tfrac12\lvert E_T\rvert$.

### 5.7.4 Choosing $n$ from an error bound

The error formulas give a guaranteed step size when a bound on the relevant derivative is available.

> [!example] Worked Example 5.33 — How many samples for six decimals?
> Find $n$ so that the error in $\int_0^1e^{-x^2}dx$ is below $10^{-6}$.
>
> **Trapezoidal.** $f''=(4x^2-2)e^{-x^2}$; on $[0,1]$, $\lvert f''\rvert\le2$ (at $x=0$). Need $\frac{1\cdot h^2}{12}\cdot2\le10^{-6}$, i.e. $h\le2.449\times10^{-3}$, so $n\ge409$. (Actual error at $n=410$: $3.6\times10^{-7}$.)
>
> **Simpson.** $f^{(4)}=(16x^4-48x^2+12)e^{-x^2}$; $\lvert f^{(4)}\rvert\le12$ on $[0,1]$. Need $\frac{h^4}{180}\cdot12\le10^{-6}$, i.e. $h\le0.0622$, so $n\ge16.07$; take $n=18$ (even). (Actual error: $7.8\times10^{-8}$.)
>
> The bounds are conservative — the worst-case derivative occurs only at one end — but they are *guarantees*. Twenty-odd samples for Simpson versus four hundred for the trapezoidal rule is the practical meaning of "fourth order versus second order".

### 5.7.5 Error estimation without derivatives: Richardson extrapolation

Derivative bounds are rarely available for measured data or complicated formulas. Instead, compute with $h$ and $h/2$ and use the known order. If $I=T_n+Ch^2+O(h^4)$, then $I=T_{2n}+C\frac{h^2}{4}+O(h^4)$, and eliminating $C$:

$$
I\approx T_{2n}+\frac{T_{2n}-T_n}{3},\qquad \text{error of }T_{2n}\approx\frac{T_{2n}-T_n}{3}.
$$

From Example 5.32, $\frac{T_8-T_4}{3}=\frac{0.745866-0.742984}{3}=9.61\times10^{-4}$, an excellent estimate of the actual error $9.59\times10^{-4}$; the extrapolated value $0.746826$ is identical to $S_8$ — Richardson extrapolation of the trapezoidal rule *is* Simpson's rule. Repeating the idea to eliminate $h^4$, $h^6,\dots$ gives Romberg integration ([Ch. 22](../part-7-numerical-computational/22-numerical-methods.md)).

> [!tip] Engineering Note — The trapezoidal rule and periodic signals
> For a smooth **periodic** function integrated over a whole period, the trapezoidal rule is not second order but converges *exponentially* fast: the error terms of the Euler–Maclaurin expansion are differences of derivatives at the two ends, which cancel by periodicity. Example 5.31's $64$-point sum was accurate to $15$ digits. This is why the DFT — essentially the trapezoidal rule applied to $\frac1T\int_0^Tx(t)e^{-jk\omega_0t}dt$ — computes Fourier coefficients of band-limited periodic signals exactly, and why sampled RMS and power meters work so well on steady-state waveforms.

### 5.7.6 Sampled data

> [!example] Worked Example 5.34 — Energy of a measured pulse
> An ADC samples the current through a $10\ \Omega$ resistor every $1$ ms:
>
> | $t$ (ms) | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
> |---|---|---|---|---|---|---|---|
> | $i$ (A) | 0 | 2.1 | 3.4 | 3.9 | 3.6 | 2.5 | 1.0 |
> | $p=Ri^2$ (W) | 0 | 44.1 | 115.6 | 152.1 | 129.6 | 62.5 | 10.0 |
>
> Estimate the energy delivered, $W=\int_0^{6\,\text{ms}}p\,dt$.
>
> **Trapezoidal:** $W\approx10^{-3}\big[\tfrac12(0+10)+44.1+115.6+152.1+129.6+62.5\big]=0.5089$ J.
>
> **Simpson** ($n=6$, even): $W\approx\frac{10^{-3}}{3}\big[0+4(44.1+152.1+62.5)+2(115.6+129.6)+10\big]=0.5117$ J.
>
> The $0.5\%$ difference indicates the sampling is adequate for an energy estimate to about two significant figures. Note that we integrated $p=Ri^2$, **not** $R\,(\int i\,dt)^2$ or $R\,(\bar i)^2$ — the square must be taken sample by sample before integrating (§5.9).

### 5.7.7 Difficult integrands

The error formulas assume smooth integrands. When that fails:

- **Discontinuities or kinks** (rectified waves, PWM): place nodes at the breakpoints and integrate each smooth piece separately; otherwise the error drops only like $h$.
- **Integrable singularities** such as $\int_0^1\frac{\cos x}{\sqrt x}dx$: remove the singularity by substitution ($x=u^2$ gives $\int_0^12\cos(u^2)\,du$, which is smooth), or subtract it off analytically.
- **Infinite ranges:** truncate where the tail is provably small (bound it by a comparison integral, §5.3.4), or map to a finite interval with a substitution such as $x=\frac{u}{1-u}$.
- **Rapid oscillation** ($\int f(x)\cos(\omega x)\,dx$ with large $\omega$): ensure several samples per oscillation period, or use specialised (Filon-type, FFT-based) methods.

> [!warning] Common Pitfall — Simpson's rule needs an even number of panels
> With an odd number of panels the composite Simpson pattern $1,4,2,4,\dots,4,1$ does not close. Either use Simpson's $\tfrac38$ rule, $\frac{3h}{8}(f_0+3f_1+3f_2+f_3)$, on the last three panels, or resample. Silently applying the even-$n$ formula to odd $n$ gives an answer that can be badly wrong.

## 5.8 Applications to charge, energy, and accumulated quantities

### 5.8.1 The accumulation principle

Every application in this section follows one pattern, worth stating explicitly because it is how integrals are *set up* rather than evaluated:

> [!definition] The differential-element method
> To find a total $Q$ that is built from a varying density or rate:
> 1. **Slice** the domain (time, length, area, volume) into small elements $dt$, $dx$, $dA$ or $dV$.
> 2. On one element the density is essentially constant, so write the **contribution** $dQ=(\text{density})\times(\text{size of element})$.
> 3. **Integrate** $dQ$ over the domain, choosing limits that cover each element exactly once.
> 4. **Check** units, limiting cases and signs.

The principle covers charge from current ($dq=i\,dt$), charge from charge density ($dQ=\lambda\,d\ell$, $\sigma\,dA$, $\rho\,dV$), energy from power ($dW=p\,dt$), work from force ($dW=F\,dx$), resistance of a non-uniform conductor ($dR=\rho_e\,dx/A(x)$), flux from flux density ($d\Phi=B\,dA$), and probability from a density ($dP=f_X(x)\,dx$, [Ch. 19](../part-6-probability-statistics/19-probability.md)). The table summarises the most common cases.

| Accumulated quantity | Element | Integral | Units check |
|---|---|---|---|
| Charge from current | $dq=i\,dt$ | $q=\int i\,dt$ | A·s = C |
| Charge from line density | $dQ=\lambda\,d\ell$ | $Q=\int\lambda\,d\ell$ | (C/m)·m |
| Charge from surface density | $dQ=\sigma\,dA$ | $Q=\iint\sigma\,dA$ | (C/m²)·m² |
| Charge from volume density | $dQ=\rho\,dV$ | $Q=\iiint\rho\,dV$ | (C/m³)·m³ |
| Energy from power | $dW=p\,dt$ | $W=\int p\,dt$ | W·s = J |
| Work against a force | $dW=F\,dx$ | $W=\int F\,dx$ | N·m = J |
| Flux linkage from voltage | $d\lambda=v\,dt$ | $\lambda=\int v\,dt$ | V·s = Wb |
| Resistance, non-uniform bar | $dR=\rho_e\,dx/A(x)$ | $R=\int\rho_e\,dx/A$ | (Ω·m)·m/m² = Ω |
| Battery charge | $dq=i\,dt$ | $q=\int i\,dt$ | A·h = 3600 C |

### 5.8.2 Charge from current

Since $i=dq/dt$, the FTC gives

$$
q(t)=q(t_0)+\int_{t_0}^{t}i(\tau)\,d\tau ,
$$

and the **net** charge passing a point between $t_1$ and $t_2$ is $\int_{t_1}^{t_2}i\,dt$, which may be zero even when large currents flow (pure AC). The **total** charge moved irrespective of direction is $\int\lvert i\rvert\,dt$ — the quantity that matters for electroplating, electrode wear and battery throughput.

> [!example] Worked Example 5.35 — Charge in a half-cycle
> A diode passes the positive half-cycles of $i(t)=I_m\sin\omega t$ with $I_m=10$ A at $f=50$ Hz. Find the charge passed per half-cycle and the average current.
>
> **Step 1 — integrate over one half-cycle** $0\le t\le\pi/\omega$:
> $$
> Q=\int_0^{\pi/\omega}I_m\sin\omega t\,dt=\frac{I_m}{\omega}\big[-\cos\omega t\big]_0^{\pi/\omega}=\frac{2I_m}{\omega}.
> $$
> **Step 2 — numbers.** $\omega=100\pi$ rad/s: $Q=\dfrac{20}{100\pi}=0.06366$ C $=63.66$ mC.
>
> **Step 3 — average over a full period** $T=20$ ms (the diode blocks during the negative half):
> $$
> I_{\text{avg}}=\frac{Q}{T}=\frac{0.06366}{0.02}=3.183\ \text{A}=\frac{I_m}{\pi}.
> $$
> This is the DC component of a half-wave rectified current (§5.9).

> [!example] Worked Example 5.36 — Battery capacity under a load profile
> A sensor node draws $i(t)=20$ mA for $5$ s while transmitting, then $50\ \mu$A for the remaining $595$ s of each $10$-minute cycle. How long does a $2000$ mAh battery last (ignoring self-discharge)?
>
> **Charge per cycle:** $q=20\,\text{mA}\times5\,\text{s}+0.05\,\text{mA}\times595\,\text{s}=100+29.75=129.75$ mA·s.
>
> **Average current:** $\bar i=129.75/600=0.21625$ mA.
>
> **Lifetime:** $2000\ \text{mAh}/0.21625\ \text{mA}=9249$ h $\approx385$ days.
>
> The integral of a piecewise-constant current is a sum of rectangles; the insight is that the "small" sleep current contributes almost a quarter of the total charge because it flows $119$ times longer.

### 5.8.3 Energy stored in capacitors and inductors

The power into any two-terminal element is $p=vi$, so the energy delivered between $t_1$ and $t_2$ is $W=\int_{t_1}^{t_2}v\,i\,dt$. For a capacitor, $i=C\,dv/dt$; for an inductor, $v=L\,di/dt$. Substituting and changing the variable of integration (§5.6) gives two of the most-used formulas in electrical engineering.

> [!theorem] Stored energy
> A linear capacitor charged from $0$ to $V$ and a linear inductor carrying current rising from $0$ to $I$ store
> $$
> W_C=\int_0^{t}v\,C\frac{dv}{dt}\,dt=C\int_0^{V}v\,dv=\tfrac12CV^2,\qquad
> W_L=\int_0^{t}L\frac{di}{dt}\,i\,dt=L\int_0^{I}i\,di=\tfrac12LI^2 .
> $$

The substitution $dv=\frac{dv}{dt}dt$ removes time entirely: the stored energy depends only on the final state, not on how fast or by what path it was reached. That is the mathematical signature of a **conservative** (lossless) element.

For a *nonlinear* capacitor described by $q=q(v)$, the same argument gives $W=\int_0^{V}v\,\frac{dq}{dv}\,dv=\int_0^{Q}v\,dq$, the area between the $q$–$v$ curve and the $q$-axis. Similarly, for a nonlinear inductor (iron core) $W=\int i\,d\lambda$, and over one cycle of an AC magnetisation the loop integral $\oint i\,d\lambda$ is the **hysteresis loss** per cycle — the area of the B–H loop scaled by core volume.

> [!example] Worked Example 5.37 — Charging a capacitor through a resistor: where does the energy go?
> A capacitor $C$ is charged from $0$ through $R$ by a step $V_s$. The current is $i=\frac{V_s}{R}e^{-t/RC}$.
>
> **Energy supplied by the source:**
> $$
> W_s=\int_0^\infty V_s\,i\,dt=\frac{V_s^2}{R}\int_0^\infty e^{-t/RC}dt=\frac{V_s^2}{R}\cdot RC=CV_s^2 .
> $$
> **Energy lost in the resistor:**
> $$
> W_R=\int_0^\infty i^2R\,dt=\frac{V_s^2}{R}\int_0^\infty e^{-2t/RC}dt=\tfrac12CV_s^2 .
> $$
> **Energy stored:** $W_C=\tfrac12CV_s^2$, and $W_s=W_R+W_C$ ✓.
>
> Exactly half the source energy is always lost, whatever $R$ is. This is a fundamental limit of charging a capacitor from a fixed voltage, and the reason switched-capacitor and adiabatic logic circuits charge capacitances in small voltage steps: charging in $N$ equal steps reduces the loss to $\tfrac{1}{2N}CV_s^2$.

> [!example] Worked Example 5.38 — Energy into an inductor with a non-trivial waveform
> A $2$ mH inductor carries $i(t)=3t^2$ A for $0\le t\le2$ s (an unusually slow ramp, for clarity). Compute the energy delivered by integrating power, then check against $\tfrac12LI^2$.
>
> $v=L\,di/dt=2\times10^{-3}\cdot6t=0.012t$ V, so $p=vi=0.012t\cdot3t^2=0.036t^3$ W and
> $$
> W=\int_0^2 0.036t^3dt=0.036\cdot\frac{16}{4}=0.144\ \text{J}.
> $$
> Final current $I=3\cdot4=12$ A: $\tfrac12LI^2=\tfrac12\cdot0.002\cdot144=0.144$ J ✓.

### 5.8.4 Work and potential

Moving a charge $q$ against an electric force requires work $W=-\int\mathbf{F}\cdot d\boldsymbol{\ell}$; along a radial line this is a single-variable integral. The **potential** of a point is the work per unit charge to bring a test charge from a reference point (usually infinity), so potential is an *improper integral* of field.

> [!example] Worked Example 5.39 — Potential of a point charge
> The field of a point charge $Q$ is $E(R)=\dfrac{Q}{4\pi\varepsilon_0R^2}$, radial. The potential at distance $a$ relative to infinity is
> $$
> V(a)=-\int_\infty^{a}E(R)\,dR=\frac{Q}{4\pi\varepsilon_0}\int_a^\infty\frac{dR}{R^2}=\frac{Q}{4\pi\varepsilon_0}\Big[-\frac1R\Big]_a^\infty=\frac{Q}{4\pi\varepsilon_0a}.
> $$
> The integral converges because $1/R^2$ has $p=2>1$ (§5.3.3). In two dimensions the field of a line charge falls only as $1/r$, and the analogous integral **diverges** — which is why the potential of an infinite line charge can only be defined relative to a finite reference radius $r_0$: $V(r)=\frac{\lambda}{2\pi\varepsilon_0}\ln\frac{r_0}{r}$.

### 5.8.5 Charge in non-uniform distributions

When a density depends on one coordinate only, a well-chosen element reduces a multiple integral to a single one: thin strips for a planar density varying in one direction, thin **rings** ($dA=2\pi r\,dr$) for a disc, thin **shells** ($dV=4\pi R^2dR$) for a sphere, thin **cylindrical shells** ($dV=2\pi r\ell\,dr$) for a cylinder of length $\ell$.

> [!example] Worked Example 5.40 — Charge in a sphere with radially varying density
> A sphere of radius $a$ holds $\rho(R)=\rho_0\big(1-\frac Ra\big)$. Find the total charge and the radius $R_{1/2}$ enclosing half of it.
>
> **Step 1 — shell element.** $dQ=\rho(R)\,4\pi R^2dR$.
>
> **Step 2 — integrate.**
> $$
> Q=4\pi\rho_0\int_0^a\Big(R^2-\frac{R^3}{a}\Big)dR=4\pi\rho_0\Big(\frac{a^3}{3}-\frac{a^3}{4}\Big)=\frac{\pi\rho_0a^3}{3}.
> $$
>
> **Step 3 — enclosed charge as a function of radius.** $Q(R)=4\pi\rho_0\big(\frac{R^3}{3}-\frac{R^4}{4a}\big)$. Setting $Q(R)=\frac{Q}{2}$ with $s=R/a$: $4\big(\frac{s^3}{3}-\frac{s^4}{4}\big)=\frac16$, i.e. $3s^4-4s^3+\tfrac12=0$. Newton's method (Ch. 4) gives $s\approx0.614$.
>
> **Check.** A uniform sphere with $\rho_0$ would hold $\frac43\pi\rho_0a^3$; ours holds a quarter of that, since the density averages to $\rho_0/4$ when weighted by $R^2$. Gauss's law then gives the field at any radius as $E=Q(R)/(4\pi\varepsilon_0R^2)$ ([Ch. 7](../part-2-calculus/07-vector-calculus.md)).

> [!example] Worked Example 5.41 — Energy stored in a coaxial cable
> A coaxial line with inner radius $a$, outer radius $b$ and dielectric permittivity $\varepsilon$ carries charge $\lambda$ per unit length. The field between the conductors is $E(r)=\dfrac{\lambda}{2\pi\varepsilon r}$, and the energy density is $w=\tfrac12\varepsilon E^2$ (J/m³). Find the stored energy per unit length and deduce the capacitance per unit length.
>
> **Step 1 — element.** A cylindrical shell of radius $r$, thickness $dr$ and unit length has volume $2\pi r\,dr$.
>
> **Step 2 — integrate.**
> $$
> W'=\int_a^b\frac12\varepsilon\Big(\frac{\lambda}{2\pi\varepsilon r}\Big)^2 2\pi r\,dr=\frac{\lambda^2}{4\pi\varepsilon}\int_a^b\frac{dr}{r}=\frac{\lambda^2}{4\pi\varepsilon}\ln\frac ba .
> $$
>
> **Step 3 — capacitance.** Equating $W'=\dfrac{\lambda^2}{2C'}$ (the energy of a capacitor with charge $\lambda$ per length) gives
> $$
> C'=\frac{2\pi\varepsilon}{\ln(b/a)} .
> $$
> For RG-58-like dimensions, $b/a=3.5$, polyethylene $\varepsilon=2.25\varepsilon_0$: $C'=\dfrac{2\pi(2.25)(8.854\times10^{-12})}{\ln3.5}=99.9$ pF/m, close to the datasheet value of about $100$ pF/m. The logarithm comes directly from $\int dr/r$ — the borderline $p=1$ integrand of §5.3.

### 5.8.6 Resistance of a non-uniform conductor

A conductor whose cross-section $A(x)$ varies slowly along its length behaves as a series chain of thin slices, each with resistance $dR=\rho_e\,dx/A(x)$. Series resistances add, so

$$
R=\int_0^L\frac{\rho_e\,dx}{A(x)} .
$$

> [!example] Worked Example 5.42 — A tapered wire
> A wire of resistivity $\rho_e$ and length $L$ tapers linearly in radius from $a$ to $b$. Find $R$.
>
> **Step 1 — area.** $r(x)=a+\frac{(b-a)x}{L}$, $A=\pi r^2$.
>
> **Step 2 — substitute** $u=r(x)$, $dx=\frac{L}{b-a}du$:
> $$
> R=\frac{\rho_e}{\pi}\cdot\frac{L}{b-a}\int_a^b\frac{du}{u^2}=\frac{\rho_eL}{\pi(b-a)}\Big(\frac1a-\frac1b\Big)=\frac{\rho_eL}{\pi ab}.
> $$
>
> **Check.** For $b\to a$ this reduces to $\rho_eL/(\pi a^2)$ ✓. Note that the answer uses the *geometric* mean of the end areas, $\pi ab$, not the arithmetic mean — averaging the area first and then computing $\rho_eL/\bar A$ underestimates $R$. (Copper, $\rho_e=1.72\times10^{-8}\ \Omega$·m, $L=1$ m, $a=0.5$ mm, $b=1$ mm: $R=10.95$ mΩ.)

> [!warning] Common Pitfall — Averaging the wrong quantity
> Examples 5.42 and 5.34 make the same point: the integral of a nonlinear function of a variable is not the function of the integral. $\int\frac{dx}{A}\ne\frac{L^2}{\int A\,dx}$ and $\int i^2dt\ne\frac{1}{T}\big(\int i\,dt\big)^2$. Always build $dQ$ from the local physics and integrate *that*.

> [!tip] Engineering Note — Accumulated quantities in instrumentation
> Energy meters, coulomb counters (battery "fuel gauges"), and flux meters are hardware integrators. Digital versions accumulate samples, which is a rectangular or trapezoidal rule (§5.7). Their main error sources are exactly those of numerical integration — sampling too slowly for the waveform's bandwidth — plus offset drift (§5.1 Engineering Note), which an integrator turns into an ever-growing error. Coulomb counters therefore periodically re-calibrate against a known state such as "fully charged".

## 5.9 Mean and RMS values of functions

### 5.9.1 Motivation

A mains socket is "230 V", a USB charger "5 V", a class-D amplifier output "20 V RMS". These single numbers summarise waveforms that vary continuously in time. Two summaries dominate engineering practice: the **mean** (average, DC) value, which tells how much net charge or flux a waveform delivers, and the **root-mean-square (RMS)** value, which tells how much *power* it delivers to a resistor. Both are integrals, and both appear on every multimeter, datasheet and power-supply specification.

### 5.9.2 The mean value

> [!definition] Mean (average) value
> The mean value of $f$ over $[a,b]$ is
> $$
> \bar f=\langle f\rangle=\frac{1}{b-a}\int_a^b f(t)\,dt .
> $$
> For a $T$-periodic signal the **DC value** is the mean over one period, $\displaystyle V_{\text{dc}}=\frac1T\int_{t_0}^{t_0+T}v(t)\,dt$ (any $t_0$). For a non-periodic signal on $(-\infty,\infty)$ the **time average** is $\displaystyle\lim_{T\to\infty}\frac{1}{2T}\int_{-T}^{T}v(t)\,dt$ when the limit exists.

Geometrically $\bar f$ is the height of the rectangle on $[a,b]$ with the same area as $f$ — the Mean Value Theorem for integrals (§5.1.3) guarantees that a continuous $f$ actually attains it. Physically, the mean current $\bar i$ is the steady current that would carry the same charge in the same time (Example 5.35).

Basic properties follow from linearity: $\langle\alpha f+\beta g\rangle=\alpha\langle f\rangle+\beta\langle g\rangle$; the mean of a constant is itself; and the mean of any pure sinusoid over a whole number of periods is **zero**. Hence a signal $v(t)=V_0+\sum_n c_n\cos(n\omega_0t+\theta_n)$ has DC value exactly $V_0$. A capacitor in series blocks precisely this term; an ideal inductor in steady state has zero *mean* voltage (otherwise its current would grow without bound, since $\Delta i=\frac1L\int v\,dt$), and a capacitor in steady state has zero *mean* current. These "volt-second balance" and "charge balance" conditions are the main analytic tools for switched-mode power converters.

> [!example] Worked Example 5.43 — Volt-second balance in a buck converter
> In a buck converter the inductor voltage is $V_{\text{in}}-V_o$ for a fraction $D$ of each switching period $T$ and $-V_o$ for the rest. In steady state the mean inductor voltage is zero:
> $$
> \frac1T\Big[(V_{\text{in}}-V_o)DT+(-V_o)(1-D)T\Big]=0\quad\Longrightarrow\quad V_o=DV_{\text{in}} .
> $$
> With $V_{\text{in}}=12$ V and $D=0.25$, $V_o=3$ V. One integral — a sum of two rectangles — gives the converter's fundamental conversion ratio.

### 5.9.3 The RMS value

Feed a waveform $v(t)$ into a resistor $R$. The instantaneous power is $v^2/R$ and the mean power over a period is $\frac{1}{T}\int_0^T\frac{v^2}{R}dt$. The **effective** or RMS voltage is the DC voltage that would produce the same mean power:

$$
\frac{V_{\text{rms}}^2}{R}=\frac1T\int_0^T\frac{v(t)^2}{R}\,dt .
$$

> [!definition] RMS value
> $$
> F_{\text{rms}}=\sqrt{\frac{1}{b-a}\int_a^b f(t)^2\,dt}=\sqrt{\langle f^2\rangle} ,
> $$
> literally the **root** of the **mean** of the **square**. For periodic signals take $[a,b]$ as one period. For complex-valued signals use $\lvert f\rvert^2$.

Since $f^2\ge0$, the RMS value is always non-negative, is zero only for a signal that is zero almost everywhere, and satisfies $F_{\text{rms}}\ge\lvert\bar f\rvert$ (because $\langle(f-\bar f)^2\rangle\ge0$ expands to $\langle f^2\rangle-\bar f^2\ge0$). Dimensionally, $F_{\text{rms}}$ has the same units as $f$.

**RMS of a sinusoid.** For $v=V_m\sin(\omega t+\theta)$, Example 5.6(c) gives $\frac1T\int_0^T\sin^2(\omega t+\theta)\,dt=\frac12$, so

$$
V_{\text{rms}}=\frac{V_m}{\sqrt2}\approx0.7071\,V_m ,
$$

independent of frequency and phase. A "230 V" mains supply therefore has peak $230\sqrt2=325$ V.

```desmos-graph
left=0; right=12.566; bottom=-1.3; top=1.3
---
y=\sin(x)|blue
y=\sin(x)^2|red|dashed
y=0.5|red|dotted
y=0.7071|green
```
*$\sin t$ (blue), its square (red, dashed) with mean $1/2$ (red, dotted), and the RMS level $1/\sqrt2\approx0.707$ (green).*

### 5.9.4 Mean and RMS of standard waveforms

Every entry below is a short integral over one period; amplitude $V_m$ means the peak value.

| Waveform | Mean $V_{\text{dc}}$ | Mean of $\lvert v\rvert$ | RMS | Crest factor $V_m/V_{\text{rms}}$ |
|---|---|---|---|---|
| Sine $V_m\sin\omega t$ | $0$ | $\dfrac{2V_m}{\pi}$ | $\dfrac{V_m}{\sqrt2}$ | $\sqrt2=1.414$ |
| Full-wave rectified sine $V_m\lvert\sin\omega t\rvert$ | $\dfrac{2V_m}{\pi}$ | $\dfrac{2V_m}{\pi}$ | $\dfrac{V_m}{\sqrt2}$ | $1.414$ |
| Half-wave rectified sine | $\dfrac{V_m}{\pi}$ | $\dfrac{V_m}{\pi}$ | $\dfrac{V_m}{2}$ | $2$ |
| Square wave $\pm V_m$ | $0$ | $V_m$ | $V_m$ | $1$ |
| Triangle or sawtooth, $\pm V_m$ | $0$ | $\dfrac{V_m}{2}$ | $\dfrac{V_m}{\sqrt3}$ | $\sqrt3=1.732$ |
| Pulse train, $0$ to $V_m$, duty $D$ | $DV_m$ | $DV_m$ | $\sqrt D\,V_m$ | $\dfrac{1}{\sqrt D}$ |

> [!example] Worked Example 5.44 — Deriving three table entries
> (a) **Half-wave rectified sine.** Over $T=2\pi/\omega$, $v=V_m\sin\omega t$ for $0<t<T/2$ and $0$ otherwise.
> $$
> V_{\text{dc}}=\frac1T\int_0^{T/2}V_m\sin\omega t\,dt=\frac{1}{T}\cdot\frac{2V_m}{\omega}=\frac{V_m}{\pi},
> \qquad
> V_{\text{rms}}^2=\frac1T\int_0^{T/2}V_m^2\sin^2\omega t\,dt=\frac1T\cdot\frac{V_m^2}{2}\cdot\frac T2=\frac{V_m^2}{4}.
> $$
> So $V_{\text{rms}}=V_m/2$: half the period carries the same power density as a full sine, so the mean square halves.
>
> (b) **Sawtooth** $v=V_m\,\dfrac{2t-T}{T}$ on $[0,T]$ (from $-V_m$ to $V_m$). Substitute $u=\frac{2t-T}{T}$, $dt=\frac T2du$, $u\in[-1,1]$:
> $$
> V_{\text{rms}}^2=\frac1T\int_{-1}^{1}V_m^2u^2\,\frac T2du=\frac{V_m^2}{2}\cdot\frac23=\frac{V_m^2}{3}.
> $$
> A triangle wave is made of the same linear pieces rearranged, so it has the same RMS value $V_m/\sqrt3$. RMS depends only on the *distribution of values*, not on their order in time.
>
> (c) **Pulse train** of height $V_m$ and duty cycle $D$: $v^2=V_m^2$ for a fraction $D$ of the time, so $\langle v^2\rangle=DV_m^2$ and $V_{\text{rms}}=\sqrt D\,V_m$, while $V_{\text{dc}}=DV_m$. A $12$ V PWM signal at $D=0.25$ has mean $3$ V but RMS $6$ V: a resistive heater driven by it dissipates $(6)^2/R$, four times what a $3$ V DC supply would produce.

```desmos-graph
left=0; right=12.566; bottom=-0.2; top=1.2
---
y=\left|\sin(x)\right||blue
y=\frac{2}{\pi}|red|dashed
y=0.7071|green|dotted
```
*Full-wave rectified sine (blue) with its mean $2/\pi\approx0.637$ (red) and RMS $1/\sqrt2\approx0.707$ (green). The RMS exceeds the mean, as it always does unless the waveform is constant.*

### 5.9.5 Mean, RMS and their combination

Write a periodic signal as its DC value plus an AC part with zero mean: $v=V_{\text{dc}}+v_{\text{ac}}$. Then

$$
\langle v^2\rangle=V_{\text{dc}}^2+2V_{\text{dc}}\langle v_{\text{ac}}\rangle+\langle v_{\text{ac}}^2\rangle=V_{\text{dc}}^2+V_{\text{ac,rms}}^2 ,
$$

because $\langle v_{\text{ac}}\rangle=0$. Thus

$$
V_{\text{rms}}=\sqrt{V_{\text{dc}}^2+V_{\text{ac,rms}}^2} .
$$

In statistical language ([Ch. 19](../part-6-probability-statistics/19-probability.md)), $V_{\text{dc}}$ is the mean, $V_{\text{ac,rms}}$ is the standard deviation, and this is the identity $E[X^2]=\mu^2+\sigma^2$. The RMS of a zero-mean noise voltage is its standard deviation.

Extending the argument to several harmonics uses orthogonality (§5.2.4): every cross term $\langle\cos m\omega_0t\cos n\omega_0t\rangle$ with $m\ne n$ vanishes.

> [!theorem] RMS of a sum of harmonics (Parseval)
> If $v(t)=V_0+\displaystyle\sum_{n=1}^{N}V_n\cos(n\omega_0t+\theta_n)$, then
> $$
> V_{\text{rms}}^2=V_0^2+\sum_{n=1}^{N}\frac{V_n^2}{2}=V_0^2+\sum_{n=1}^{N}V_{n,\text{rms}}^2 .
> $$
> Powers of different frequencies add; amplitudes do not.

> [!warning] Common Pitfall — Adding RMS values
> RMS values of different-frequency components add **in quadrature**, not linearly: $3$ V RMS at $1$ kHz plus $4$ V RMS at $3$ kHz is $5$ V RMS, not $7$ V. Only components at the *same* frequency combine as phasors — and then their phase matters: two $3$ V RMS sinusoids of the same frequency give anywhere from $0$ to $6$ V RMS.

> [!example] Worked Example 5.45 — A signal with DC and harmonics
> $v(t)=5+10\sin\omega t+4\sin3\omega t$ V. Find $V_{\text{dc}}$, $V_{\text{rms}}$ and the power into $8\ \Omega$.
>
> **DC:** $V_{\text{dc}}=5$ V.
>
> **RMS:** $V_{\text{rms}}^2=5^2+\frac{10^2}{2}+\frac{4^2}{2}=25+50+8=83$, so $V_{\text{rms}}=9.110$ V.
>
> **Power:** $P=V_{\text{rms}}^2/R=83/8=10.375$ W, of which $25/8=3.125$ W is due to the DC component.
>
> The **total harmonic distortion** of the AC part, $\text{THD}=\frac{\sqrt{\sum_{n\ge2}V_{n,\text{rms}}^2}}{V_{1,\text{rms}}}=\frac{4/\sqrt2}{10/\sqrt2}=40\%$, is itself an RMS ratio.

> [!example] Worked Example 5.46 — Ripple on a DC supply
> A $12$ V supply carries a triangular ripple of $1$ V peak-to-peak. Its AC part is a triangle of amplitude $0.5$ V, with RMS $0.5/\sqrt3=0.2887$ V. Then
> $$
> V_{\text{rms}}=\sqrt{12^2+0.2887^2}=12.0035\ \text{V} .
> $$
> The ripple raises the RMS by only $0.03\%$; RMS is dominated by the largest component. That is why ripple is specified separately (in mV RMS or peak-to-peak) rather than inferred from a total RMS reading.

### 5.9.6 Average power and power factor

For a two-terminal element with periodic $v(t)$ and $i(t)$, the **average (real) power** is the mean of the instantaneous power:

$$
P=\frac1T\int_0^T v(t)\,i(t)\,dt .
$$

> [!theorem] Average power for sinusoids
> If $v=V_m\cos\omega t$ and $i=I_m\cos(\omega t-\varphi)$, then
> $$
> P=\tfrac12V_mI_m\cos\varphi=V_{\text{rms}}I_{\text{rms}}\cos\varphi=\tfrac12\operatorname{Re}\{\mathbf{V}\mathbf{I}^*\},
> $$
> where $\mathbf V=V_m$ and $\mathbf I=I_me^{-j\varphi}$ are peak phasors.

*Derivation.* By the product-to-sum identity, $v\,i=\tfrac12V_mI_m\big[\cos\varphi+\cos(2\omega t-\varphi)\big]$. The second term oscillates at $2\omega$ and averages to zero over a period, leaving $\tfrac12V_mI_m\cos\varphi$. In phasor form $\mathbf V\mathbf I^*=V_mI_me^{j\varphi}$, whose real part is $V_mI_m\cos\varphi$. $\square$

The factor $\cos\varphi$ is the **power factor**. More generally the power factor is $P/(V_{\text{rms}}I_{\text{rms}})$, which also accounts for harmonic distortion: by orthogonality, only voltage and current harmonics of the *same* order contribute to $P$, while every harmonic increases $I_{\text{rms}}$.

> [!example] Worked Example 5.47 — Real power of an inductive load
> A load on $230$ V RMS mains ($V_m=325$ V) draws $I_m=10$ A lagging by $30°$.
> $$
> P=\tfrac12(325)(10)\cos30°=1407\ \text{W},\qquad S=V_{\text{rms}}I_{\text{rms}}=\tfrac12(325)(10)=1625\ \text{VA}.
> $$
> **Check by direct integration:** numerically averaging $325\cos\omega t\cdot10\cos(\omega t-\pi/6)$ over one period gives $1407.3$ W ✓.

### 5.9.7 Form factor, crest factor and meters

Two dimensionless ratios describe waveform shape:

$$
\text{form factor}=\frac{V_{\text{rms}}}{\langle\lvert v\rvert\rangle},\qquad \text{crest factor}=\frac{V_{\text{peak}}}{V_{\text{rms}}} .
$$

For a sine, the form factor is $\frac{V_m/\sqrt2}{2V_m/\pi}=\frac{\pi}{2\sqrt2}=1.1107$ and the crest factor is $\sqrt2$.

Inexpensive "average-responding" AC meters rectify the signal, measure $\langle\lvert v\rvert\rangle$ (an easy analogue integral), and multiply by $1.1107$. They read correctly only for sinusoids.

> [!example] Worked Example 5.48 — Meter error on non-sinusoidal waveforms
> (a) **Square wave $\pm V_m$.** $\langle\lvert v\rvert\rangle=V_m$, so the meter displays $1.1107V_m$; the true RMS is $V_m$. The reading is $11.1\%$ high.
>
> (b) **Triangle wave $\pm V_m$.** $\langle\lvert v\rvert\rangle=V_m/2$, display $0.5554V_m$; true RMS $V_m/\sqrt3=0.5774V_m$. The reading is $3.8\%$ low.
>
> (c) **PWM, $D=0.1$.** $\langle v\rangle=0.1V_m$, display $0.111V_m$; true RMS $\sqrt{0.1}\,V_m=0.316V_m$ — the meter reads barely a third of the true value. Signals with high crest factor (here $\sqrt{10}=3.16$) are where average-responding meters fail worst, and where a "true-RMS" meter, which computes $\sqrt{\langle v^2\rangle}$ by squaring and integrating, is essential.

### 5.9.8 Mean and RMS over a finite window

For non-periodic signals the value depends on the averaging window, which must therefore be stated.

> [!example] Worked Example 5.49 — A decaying pulse
> $v(t)=10e^{-t/\tau}$ V with $\tau=1$ ms, averaged over the window $0\le t\le T=5$ ms.
> $$
> \bar v=\frac{1}{T}\int_0^T10e^{-t/\tau}dt=\frac{10\tau}{T}\big(1-e^{-T/\tau}\big)=1.987\ \text{V},
> $$
> $$
> V_{\text{rms}}^2=\frac1T\int_0^T100e^{-2t/\tau}dt=\frac{100\tau}{2T}\big(1-e^{-2T/\tau}\big)=9.9995,\qquad V_{\text{rms}}=3.162\ \text{V}.
> $$
> Doubling the window to $10$ ms roughly halves $\bar v$ and divides $V_{\text{rms}}$ by about $\sqrt2$ — a transient's "average" is meaningful only together with its window. For such signals the window-independent quantities are the **total** integral $\int_0^\infty v\,dt=10\tau=10$ mV·s and the **energy** $\int_0^\infty v^2dt=50\tau=0.05$ V²·s (the energy delivered to $1\ \Omega$), which lead to the energy-signal concepts of [Ch. 16](../part-5-transforms/16-fourier-transforms.md).

> [!tip] Engineering Note — RMS in noise and measurement
> For random signals the time-average RMS equals the ensemble standard deviation (for ergodic processes, [Ch. 21](../part-6-probability-statistics/21-random-processes.md)). Thermal noise in a resistor has RMS voltage $\sqrt{4kTRB}$ — the square root of an *integral* of a flat power spectral density over the bandwidth $B$. Uncorrelated noise sources add in quadrature, for exactly the orthogonality reason that harmonics do: the mean of the product of two independent zero-mean signals is zero.

## Summary

| Topic | Key formula or result |
|---|---|
| Definite integral | $\displaystyle\int_a^bf\,dx=\lim_{\lVert P\rVert\to0}\sum f(x_k^*)\,\Delta x_k$ |
| FTC | $\dfrac{d}{dx}\displaystyle\int_a^xf(t)\,dt=f(x)$; $\displaystyle\int_a^bf\,dx=F(b)-F(a)$ |
| Leibniz rule | $\dfrac{d}{dx}\displaystyle\int_{u(x)}^{v(x)}f\,dt=f(v)v'-f(u)u'$ |
| Mean value | $\bar f=\dfrac{1}{b-a}\displaystyle\int_a^bf\,dx=f(c)$ for some $c$ |
| Complex exponential | $\displaystyle\int e^{(a+jb)t}dt=\frac{e^{(a+jb)t}}{a+jb}$ |
| Orthogonality | $\displaystyle\int_0^T\sin m\omega_0t\,\sin n\omega_0t\,dt=\tfrac T2\,\delta_{mn}$ |
| Cover-up rule | $A=\big[(x-a)P/Q\big]_{x=a}$ |
| $p$-integrals | $\int_1^\infty x^{-p}dx<\infty\iff p>1$; $\int_0^1x^{-p}dx<\infty\iff p<1$ |
| Gaussian, Gamma | $\int_{-\infty}^\infty e^{-x^2}dx=\sqrt\pi$; $\Gamma(n+1)=n!$ |
| Fubini | $\iint_Rf\,dA=\int_a^b\!\int_c^df\,dy\,dx=\int_c^d\!\int_a^bf\,dx\,dy$ |
| Polar area element | $dA=r\,dr\,d\phi$ |
| Integration by parts | $\int u\,dv=uv-\int v\,du$ |
| Substitution | $\int_a^bf(g(x))g'(x)\,dx=\int_{g(a)}^{g(b)}f(u)\,du$ |
| Scaling | $\int_{-\infty}^{\infty}f(at)\,dt=\frac{1}{\lvert a\rvert}\int_{-\infty}^{\infty}f(u)\,du$ |
| Trapezoidal rule | $T_n=h[\tfrac12f_0+f_1+\dots+f_{n-1}+\tfrac12f_n]$, error $-\frac{(b-a)h^2}{12}f''(\xi)$ |
| Simpson's rule | $S_n=\frac h3[f_0+4f_1+2f_2+\dots+4f_{n-1}+f_n]$, error $-\frac{(b-a)h^4}{180}f^{(4)}(\xi)$ |
| Richardson | $I\approx T_{2n}+(T_{2n}-T_n)/3$ |
| Charge, energy | $q=\int i\,dt$; $W=\int p\,dt$; $W_C=\tfrac12CV^2$; $W_L=\tfrac12LI^2$ |
| Non-uniform resistance | $R=\int_0^L\rho_e\,dx/A(x)$ |
| RMS | $F_{\text{rms}}=\sqrt{\frac1T\int_0^Tf^2dt}$; sine: $V_m/\sqrt2$ |
| DC + AC | $V_{\text{rms}}^2=V_{\text{dc}}^2+V_{\text{ac,rms}}^2=V_0^2+\sum V_n^2/2$ |
| Average power | $P=\frac1T\int_0^Tv\,i\,dt=\tfrac12V_mI_m\cos\varphi=\tfrac12\operatorname{Re}\{\mathbf V\mathbf I^*\}$ |

## Further reading

- E. Kreyszig, *Advanced Engineering Mathematics*, 10th ed., Wiley — appendix review of integration and the chapters on multiple integrals.
- G. B. Thomas, M. D. Weir and J. Hass, *Thomas' Calculus*, Pearson — thorough treatment of techniques of integration, improper integrals and multiple integrals with many worked examples.
- K. A. Stroud and D. J. Booth, *Engineering Mathematics*, Palgrave Macmillan — programmed-learning approach to integration techniques, mean and RMS values, and numerical integration.
- R. L. Burden and J. D. Faires, *Numerical Analysis*, Cengage — derivation and error analysis of Newton–Cotes rules and Richardson extrapolation (continued in [Ch. 22](../part-7-numerical-computational/22-numerical-methods.md)).
