# Chapter 2 — Complex Numbers and Engineering Representation

[← Previous](01-mathematical-language-and-algebra.md) · [Table of Contents](../00-front/table-of-contents.md) · [Next →](03-trigonometry-and-periodic-functions.md)

> **Prerequisites:** [Ch. 1](01-mathematical-language-and-algebra.md) (algebra, polynomials, exponentials and logarithms, trigonometric identities, polar coordinates)
> **Learning outcomes:** After this chapter you can
> - move freely between rectangular, polar and exponential forms of a complex number, getting the quadrant of the angle right every time;
> - add, multiply, divide, conjugate and take magnitudes of complex numbers, choosing the form that makes each operation easiest;
> - use Euler's formula to turn trigonometric identities, sums and integrals into algebra with exponentials;
> - compute integer powers and all $n$-th roots of a complex number, and use roots of unity;
> - interpret complex arithmetic geometrically and describe curves and regions of the complex plane;
> - represent sinusoids as phasors, solve sinusoidal steady-state problems algebraically, and express impedances and transfer functions as complex functions of $j\omega$ or $s$.

A sinusoidal signal is described by two numbers — an amplitude and a phase. A complex number is also two numbers, and complex multiplication combines magnitudes and angles in exactly the way that linear circuits combine amplitudes and phases. That coincidence, which turns out to be no coincidence at all, is why every alternating-current calculation in electrical engineering, every frequency response, every Fourier and Laplace transform and every digital filter is written in complex numbers.

Historically, complex numbers were not invented to solve $x^2+1=0$; mathematicians were content to say that equation had no solution. They were forced on sixteenth-century algebraists by the *cubic* formula, which for equations with three perfectly real roots (such as $x^3-15x-4=0$, whose roots include $x=4$) produces intermediate expressions like $\sqrt{-121}$. The square roots of negative numbers cancel at the end, but there is no way to avoid them on the way. The lesson carries straight into engineering: complex numbers are a *computational device* that often lets you reach a real answer by the shortest route. A circuit's voltages are real, its transient responses are real, its measured gains are real — yet the shortest path to all of them runs through $\mathbb{C}$.

This chapter is the mathematical foundation for that path. Sections 2.1–2.3 build the algebra, Sections 2.4–2.6 the exponential machinery that makes it powerful, Section 2.7 the geometry, and Sections 2.8–2.9 show how the whole apparatus becomes the phasor method and the language of impedance and transfer functions. Complex *functions* — differentiation, contour integration and residues — are deferred to [Ch. 18](../part-5-transforms/18-complex-analysis.md).

## 2.1 Rectangular and polar forms

### Why two forms?

An engineer meets complex numbers in two guises. A network analyser reports an impedance as "$40-j30\ \Omega$" — a resistance and a reactance, the **rectangular** description. An oscilloscope comparing two sinusoids reports "amplitude ratio $0.8$, phase $-36.9^\circ$" — a magnitude and an angle, the **polar** description. They are the same complex number seen in two coordinate systems, exactly as a point in the plane can be given by $(x,y)$ or by $(r,\theta)$ (Ch. 1, Section 1.7). Fluency means converting between them without thinking and knowing which one a given operation prefers.

### Definition and the complex number system

> [!definition] Complex number
> A **complex number** is an expression $z=a+jb$ with $a,b\in\mathbb{R}$, where the **imaginary unit** $j$ satisfies $j^2=-1$. The real number $a=\operatorname{Re}\{z\}$ is the **real part** and the real number $b=\operatorname{Im}\{z\}$ is the **imaginary part** (note: $b$, not $jb$). The writing $a+jb$ is the **rectangular** (Cartesian) form. The set of all complex numbers is $\mathbb{C}$.
>
> A number with $b=0$ is **real**; a number with $a=0$, $b\neq0$ is **purely imaginary**.

> [!definition] Equality
> Two complex numbers are equal if and only if their real parts are equal **and** their imaginary parts are equal:
> $$
> a+jb=c+jd\iff a=c\ \text{ and }\ b=d.
> $$

The equality rule is more powerful than it looks: **one complex equation is two real equations.** Many derivations in circuit theory exploit exactly this — a single complex KVL equation delivers both an amplitude condition and a phase condition.

A rigorous construction, for readers who distrust a symbol "whose square is $-1$": define $\mathbb{C}$ as the set of ordered pairs $(a,b)$ of reals with

$$
(a,b)+(c,d)=(a+c,\ b+d),\qquad (a,b)\cdot(c,d)=(ac-bd,\ ad+bc).
$$

One checks that these operations obey all the usual laws of arithmetic (commutative, associative, distributive, with $(0,0)$ and $(1,0)$ as additive and multiplicative identities). Pairs $(a,0)$ behave exactly like the real numbers $a$, and the pair $j=(0,1)$ satisfies $j\cdot j=(0\cdot0-1\cdot1,\ 0)=(-1,0)=-1$. Writing $(a,b)=a+jb$ recovers the familiar notation. So complex numbers are nothing more mysterious than points of the plane with a cleverly chosen multiplication rule — the rule is what makes them useful.

> [!warning] Common Pitfall — $\mathbb{C}$ is not ordered
> There is no meaningful "$z_1<z_2$" for complex numbers. Statements such as "$3+j4>2$" are meaningless. Only real quantities derived from complex numbers — magnitudes, real parts, imaginary parts, angles — can be compared. When an engineer says "this impedance is larger", they mean $|Z|$ is larger.

### The complex plane and polar coordinates

Plot $z=a+jb$ as the point $(a,b)$, with the real axis horizontal and the imaginary axis vertical. Its distance from the origin and its angle measured anticlockwise from the positive real axis give the **polar form**.

> [!definition] Magnitude, argument, polar form
> For $z=a+jb\neq0$,
> $$
> r=|z|=\sqrt{a^2+b^2},\qquad \theta=\arg z=\angle z,\quad\text{where }\cos\theta=\frac{a}{r},\ \sin\theta=\frac{b}{r}.
> $$
> Then
> $$
> z=r(\cos\theta+j\sin\theta),\qquad a=r\cos\theta,\quad b=r\sin\theta.
> $$
> Engineers abbreviate this as $z=r\angle\theta$. The quantity $r\ge0$ is the **magnitude** (modulus, absolute value) and $\theta$ the **argument** (angle, phase). The argument is defined only up to multiples of $2\pi$; the **principal argument** $\operatorname{Arg}z$ is the value in $(-\pi,\pi]$. For $z=0$ the magnitude is $0$ and the argument is undefined.

Both $\cos\theta=a/r$ **and** $\sin\theta=b/r$ are needed to fix $\theta$; either one alone leaves two candidate angles. The ratio $\tan\theta=b/a$ also leaves two candidates, $\pi$ apart, which is the source of the most common error in this whole subject.

### The quadrant problem and atan2

The calculator function $\arctan$ returns values in $(-\pi/2,\pi/2)$ — the right half-plane only. It therefore gives the correct argument only when $a>0$. The correct rule is the two-argument arctangent:

> [!definition] Two-argument arctangent
> $$
> \operatorname{atan2}(b,a)=\begin{cases}
> \arctan(b/a), & a>0,\\
> \arctan(b/a)+\pi, & a<0,\ b\ge0,\\
> \arctan(b/a)-\pi, & a<0,\ b<0,\\
> +\pi/2, & a=0,\ b>0,\\
> -\pi/2, & a=0,\ b<0.
> \end{cases}
> $$
> It returns $\operatorname{Arg}(a+jb)\in(-\pi,\pi]$. It is available as `atan2(b, a)` in C, Python, MATLAB and spreadsheets (spreadsheets often reverse the argument order — check).

The geometric content is simple: points in the second and third quadrants have arguments of magnitude greater than $90^\circ$, and $\arctan$ can never return such a value, so a correction of $\pm180^\circ$ is required there.

```desmos-graph
left=-6; right=6; bottom=-5; top=5
---
(3,4)|blue|label:3+j4
(-3,4)|red|label:-3+j4
(-3,-4)|green|label:-3-j4
(3,-4)|purple|label:3-j4
y=\frac{4}{3}x|blue|dashed|0<x<3
y=-\frac{4}{3}x|red|dashed|-3<x<0
y=\frac{4}{3}x|green|dashed|-3<x<0
y=-\frac{4}{3}x|purple|dashed|0<x<3
x^2+y^2=25|grey|dotted
```
*The four points $\pm3\pm j4$ all have magnitude 5 and $\lvert\tan\theta\rvert=4/3$, but their arguments are $53.13^\circ$, $126.87^\circ$, $-126.87^\circ$ and $-53.13^\circ$. Only the sign pattern of $(a,b)$ distinguishes them.*

### Degrees and radians

Engineers quote phase angles in degrees ($10\angle{-30^\circ}$) because instruments do; mathematics and every exponential formula require radians. Convert with $\theta_{\text{rad}}=\theta_{\text{deg}}\cdot\pi/180$. The convention in this book: when an angle carries a degree sign it is in degrees; otherwise it is in radians. Inside an exponent ($e^{j\theta}$, Section 2.4) the angle is **always** in radians.

> [!example] Worked Example 2.1 — Rectangular to polar
> Convert $3+j4$, $-3+j4$, $-1-j$, $-2$, $j5$ and $-5-j12$ to polar form.
>
> **Step 1 — magnitude.** $\lvert3+j4\rvert=\sqrt{9+16}=5$; $\lvert-3+j4\rvert=5$; $\lvert-1-j\rvert=\sqrt2\approx1.4142$; $\lvert-2\rvert=2$; $\lvert j5\rvert=5$; $\lvert-5-j12\rvert=\sqrt{25+144}=13$.
>
> **Step 2 — locate the quadrant, then the angle.**
> - $3+j4$: first quadrant, $\theta=\arctan\frac43=53.13^\circ$ ($0.9273$ rad).
> - $-3+j4$: second quadrant, $\arctan\frac{4}{-3}=-53.13^\circ$ is wrong; add $180^\circ$: $\theta=126.87^\circ$.
> - $-1-j$: third quadrant, $\arctan\frac{-1}{-1}=45^\circ$ is wrong; subtract $180^\circ$: $\theta=-135^\circ=-3\pi/4$.
> - $-2$: negative real axis, $\theta=180^\circ=\pi$ (principal value $+\pi$, not $-\pi$).
> - $j5$: positive imaginary axis, $\theta=90^\circ$; $\arctan(5/0)$ is undefined, so use the special case.
> - $-5-j12$: third quadrant, $\arctan\frac{12}{5}=67.38^\circ$, so $\theta=67.38^\circ-180^\circ=-112.62^\circ$.
>
> | $z$ | $\lvert z\rvert$ | $\operatorname{Arg} z$ |
> |---|---|---|
> | $3+j4$ | $5$ | $53.13^\circ$ |
> | $-3+j4$ | $5$ | $126.87^\circ$ |
> | $-1-j$ | $1.4142$ | $-135^\circ$ |
> | $-2$ | $2$ | $180^\circ$ |
> | $j5$ | $5$ | $90^\circ$ |
> | $-5-j12$ | $13$ | $-112.62^\circ$ |
>
> **Check** one of them backwards: $13\cos(-112.62^\circ)=13(-0.3846)=-5.00$ and $13\sin(-112.62^\circ)=13(-0.9231)=-12.00$. ✓

> [!example] Worked Example 2.2 — Polar to rectangular
> Convert $10\angle{-30^\circ}$, $25\angle143.13^\circ$, $4\angle\frac{5\pi}{6}$ and $2\angle200^\circ$.
>
> Use $a=r\cos\theta$, $b=r\sin\theta$ — no quadrant care is needed in this direction, because cosine and sine carry the signs automatically.
> - $10\angle{-30^\circ}=10(0.8660)+j10(-0.5)=8.660-j5.000$.
> - $25\angle143.13^\circ=25(-0.8000)+j25(0.6000)=-20.00+j15.00$.
> - $4\angle\frac{5\pi}{6}=4\left(-\frac{\sqrt3}{2}\right)+j4\left(\frac12\right)=-2\sqrt3+j2\approx-3.464+j2.000$.
> - $2\angle200^\circ=2(-0.9397)+j2(-0.3420)=-1.879-j0.684$. Note that $200^\circ$ is not a principal argument; the same number is $2\angle{-160^\circ}$.

> [!example] Worked Example 2.3 — Using equality of complex numbers
> Find real $x$ and $y$ such that $(x+jy)(2-j)=7+j4$.
>
> **Expand** the left side: $(x+jy)(2-j)=2x-jx+j2y-j^2y=(2x+y)+j(2y-x)$.
>
> **Equate** real and imaginary parts — one complex equation gives two real ones:
> $$
> 2x+y=7,\qquad -x+2y=4.
> $$
> From the second, $x=2y-4$; substituting, $4y-8+y=7$, so $y=3$ and $x=2$.
>
> **Check:** $(2+j3)(2-j)=4-j2+j6-j^2 3=4+3+j4=7+j4$. ✓
>
> (Section 2.3 gives a quicker route: $x+jy=(7+j4)/(2-j)$.)

> [!example] Worked Example 2.4 — An impedance measurement
> An impedance analyser reports $\lvert Z\rvert=50\ \Omega$ and $\angle Z=-36.87^\circ$ for an unknown component at $1$ kHz. Find the series resistance and reactance, and say whether the component behaves capacitively or inductively.
>
> $$
> Z=50\cos(-36.87^\circ)+j50\sin(-36.87^\circ)=40.0-j30.0\ \Omega.
> $$
> The resistance is $R=40\ \Omega$ and the reactance $X=-30\ \Omega$. A negative reactance is capacitive (Section 2.9: $Z_C=-j/(\omega C)$). If the reactance comes from a pure capacitor, $1/(\omega C)=30$, so $C=1/(2\pi\cdot1000\cdot30)=5.31\ \mu$F.
>
> **Check:** $\sqrt{40^2+30^2}=50$ ✓ and $\operatorname{atan2}(-30,40)=-36.87^\circ$ ✓.

> [!warning] Common Pitfall — the calculator arctangent
> A calculator's $\tan^{-1}$ gives $\arctan\frac{-1}{-1}=\arctan1=45^\circ$ for $-1-j$, which is the angle of $1+j$ — the diametrically opposite point. **Always sketch the point first** (or at least note the signs of $a$ and $b$), or use `atan2`. A phase error of exactly $180^\circ$ in a circuit result is almost always this mistake: it flips the sign of the entire signal.

> [!tip] Engineering Note — which form do instruments use?
> Network and impedance analysers can display either form: "$R+jX$" (rectangular), "$\lvert Z\rvert,\theta$" (polar), or on a Smith chart (a geometric map, Section 2.7). Oscilloscopes and lock-in amplifiers naturally measure magnitude and phase, or the **in-phase** ($I$) and **quadrature** ($Q$) components $a$ and $b$. In digital communications, a received symbol is literally a complex number $I+jQ$, and a constellation diagram is a plot of the complex plane. Every one of these displays is just one of the two forms of this section.


## 2.2 Complex conjugates

### Motivation

Division by a complex number, extracting real and imaginary parts, computing a magnitude without square roots, finding average power in an AC circuit, and proving that a real filter's poles come in mirror-image pairs — all of these rest on one operation: changing the sign of $j$. Because $j$ and $-j$ are algebraically indistinguishable (both square to $-1$), replacing every $j$ by $-j$ in a true equation produces another true equation. That symmetry is the complex conjugate.

> [!definition] Complex conjugate
> The **complex conjugate** of $z=a+jb=r\angle\theta$ is
> $$
> z^*=a-jb=r\angle(-\theta).
> $$
> Geometrically it is the **reflection of $z$ in the real axis**. Mathematics texts often write $\bar z$; this book uses $z^*$, which matches the conjugate-transpose notation $\mathbf{A}^\mathsf{H}$ of linear algebra.

### Algebraic properties

> [!theorem] Properties of the conjugate
> For all $z,z_1,z_2\in\mathbb{C}$:
> $$
> \begin{aligned}
> &\text{(i)}\ \ (z^*)^*=z, \qquad \text{(ii)}\ \ (z_1\pm z_2)^*=z_1^*\pm z_2^*,\qquad \text{(iii)}\ \ (z_1z_2)^*=z_1^*z_2^*,\
> &\text{(iv)}\ \ \left(\frac{z_1}{z_2}\right)^*=\frac{z_1^*}{z_2^*}\ (z_2\ne0),\qquad \text{(v)}\ \ (z^n)^*=(z^*)^n,\
> &\text{(vi)}\ \ zz^*=a^2+b^2=|z|^2,\qquad \text{(vii)}\ \ |z^*|=|z|,\
> &\text{(viii)}\ \ \operatorname{Re}\{z\}=\frac{z+z^*}{2},\qquad \text{(ix)}\ \ \operatorname{Im}\{z\}=\frac{z-z^*}{2j},\
> &\text{(x)}\ \ z\ \text{is real}\iff z^*=z,\qquad z\ \text{is purely imaginary (or 0)}\iff z^*=-z.
> \end{aligned}
> $$

*Proof sketch.* Each is a one-line computation in rectangular form. For (iii), with $z_1=a+jb$, $z_2=c+jd$:

$$
z_1z_2=(ac-bd)+j(ad+bc)\ \Rightarrow\ (z_1z_2)^*=(ac-bd)-j(ad+bc),
$$

$$
z_1^*z_2^*=(a-jb)(c-jd)=(ac-bd)-j(ad+bc).\ \checkmark
$$

Property (v) follows from (iii) by induction, (iv) from (iii) applied to $z_1=(z_1/z_2)\,z_2$, and (vi) from $(a+jb)(a-jb)=a^2-j^2b^2=a^2+b^2$. For (viii) and (ix), add and subtract $z=a+jb$ and $z^*=a-jb$. $\square$

Property (vi) deserves emphasis: **the product of a number and its conjugate is real and non-negative**. It is the identity behind division (Section 2.3), behind the energy of a signal $\int|x(t)|^2dt=\int x\,x^*\,dt$, and behind power in circuits.

> [!theorem] Conjugate of a rational expression
> Let $f$ be built from $z$ and **real** constants by addition, subtraction, multiplication, division and integer powers (for example a polynomial or rational function with real coefficients). Then
> $$
> f(z^*)=\big(f(z)\big)^*.
> $$

*Proof.* Properties (ii)–(v) let the conjugate pass through every operation, and a real constant $c$ satisfies $c^*=c$. $\square$

The most important consequence:

> [!theorem] Conjugate root theorem
> If $p(x)=c_nx^n+\dots+c_1x+c_0$ has **real** coefficients and $p(z_0)=0$, then $p(z_0^*)=0$. Non-real roots of real polynomials occur in **conjugate pairs**, and a real polynomial of odd degree has at least one real root.

*Proof.* $p(z_0^*)=\big(p(z_0)\big)^*=0^*=0$. Pairing off the non-real roots leaves an odd number of roots to be real when $n$ is odd. $\square$

Each conjugate pair $\alpha\pm j\beta$ combines into a real quadratic factor:

$$
\big(x-(\alpha+j\beta)\big)\big(x-(\alpha-j\beta)\big)=(x-\alpha)^2+\beta^2=x^2-2\alpha x+(\alpha^2+\beta^2).
$$

So every real polynomial factors into real linear and real irreducible quadratic factors — exactly the factors that appear in the partial-fraction decompositions of Ch. 1 and Ch. 14.

> [!example] Worked Example 2.5 — Basic conjugate manipulations
> Let $z=2-j3$. Find $z^*$, $zz^*$, $\operatorname{Re}\{z\}$ via (viii), and $1/z$.
>
> - $z^*=2+j3$.
> - $zz^*=(2-j3)(2+j3)=4+9=13=|z|^2$, so $|z|=\sqrt{13}\approx3.606$.
> - $\dfrac{z+z^*}{2}=\dfrac{4}{2}=2$. ✓
> - $\dfrac1z=\dfrac{z^*}{zz^*}=\dfrac{2+j3}{13}=0.1538+j0.2308$.
>
> **Check:** $(2-j3)(0.1538+j0.2308)=0.3077+0.6923+j(0.4615-0.4615)=1.000$. ✓

> [!example] Worked Example 2.6 — An equation involving $z$ and $z^*$
> Solve $z+2z^*=3-j4$.
>
> An equation containing both $z$ and $z^*$ is **not** a polynomial equation in $z$, so substitute $z=x+jy$ and use equality:
> $$
> (x+jy)+2(x-jy)=3x-jy=3-j4.
> $$
> Real parts: $3x=3$, so $x=1$. Imaginary parts: $-y=-4$, so $y=4$. Hence $z=1+j4$.
>
> **Check:** $(1+j4)+2(1-j4)=3-j4$. ✓

> [!example] Worked Example 2.7 — Completing a factorisation
> The real cubic $p(z)=z^3-5z^2+9z-5$ has a root $z=2+j$. Find the others.
>
> **Step 1.** Coefficients are real, so $2-j$ is also a root.
>
> **Step 2.** The pair gives the real quadratic factor $(z-2)^2+1=z^2-4z+5$.
>
> **Step 3.** Divide: $z^3-5z^2+9z-5=(z^2-4z+5)(z-1)$ — check by expanding: $z^3-z^2-4z^2+4z+5z-5=z^3-5z^2+9z-5$. ✓
>
> The roots are $1$, $2+j$, $2-j$. Their sum is $5$ (the negative of the $z^2$ coefficient) and their product is $1\cdot\big((2+j)(2-j)\big)=5$ (the negative of the constant term) — both consistent with Vieta's relations. ✓

> [!example] Worked Example 2.8 — Unit magnitude and the all-pass property
> **(a)** Show that if $|z|=1$ then $1/z=z^*$.
>
> From (vi), $zz^*=|z|^2=1$, so $z^*=1/z$. Points on the unit circle are inverted by conjugation — a fact used constantly in discrete-time signal processing, where $z=e^{j\Omega}$ lies on the unit circle (Ch. 17).
>
> **(b)** Let $a$ be any complex number with $|a|\neq1$. Show that
> $$
> A(z)=\frac{z-a}{1-a^*z}
> $$
> satisfies $|A(z)|=1$ whenever $|z|=1$.
>
> **Step 1.** Use $|z|=1$ to write $1=zz^*$ in the denominator: $1-a^*z=zz^*-a^*z=z\,(z^*-a^*)$.
>
> **Step 2.** Take magnitudes: $|1-a^*z|=|z|\,|z^*-a^*|=1\cdot|(z-a)^*|=|z-a|$, by property (ii) and (vii).
>
> **Step 3.** Hence $|A(z)|=|z-a|/|z-a|=1$. $\square$
>
> **Numerical check:** $a=0.5+j0.3$, $z=e^{j1}=0.5403+j0.8415$: $z-a=0.0403+j0.5415$, $|z-a|=0.5430$; $a^*z=(0.5-j0.3)(0.5403+j0.8415)=0.5226+j0.2587$, $1-a^*z=0.4774-j0.2587$, $|1-a^*z|=0.5430$. ✓
>
> This $A(z)$ is the **first-order all-pass** section of digital filtering: it changes phase but passes every frequency with unit gain.

> [!tip] Engineering Note — where conjugates appear
> - **Power:** with amplitude phasors, complex power is $S=\tfrac12\mathbf{V}\mathbf{I}^*$ (Section 2.9). The conjugate makes the angle of $S$ equal to the voltage–current phase difference.
> - **Real signals:** the Fourier transform of a real signal satisfies $X(-\omega)=X(\omega)^*$ (Ch. 16) — a direct consequence of the theorem on conjugates of expressions with real constants.
> - **Real systems:** poles and zeros of any circuit built from real $R$, $L$, $C$ and controlled sources are real or occur in conjugate pairs, so pole–zero plots are always symmetric about the real axis.
> - **Matched filtering and correlation** use $x(t)\,y^*(t)$; **maximum power transfer** requires a load equal to the conjugate of the source impedance (Section 2.9).

> [!warning] Common Pitfall — conjugating only part of an expression
> To conjugate $\dfrac{3+j2}{(1-j)\,e^{j\pi/5}}$ you must change the sign of **every** $j$: $\dfrac{3-j2}{(1+j)\,e^{-j\pi/5}}$. A frequent error is to conjugate the numerator only. Conversely, do **not** conjugate real symbols that merely look complex: if $\omega$, $R$, $L$, $C$ are real, then $(R+j\omega L)^*=R-j\omega L$, with $\omega$ untouched.

## 2.3 Complex arithmetic

### Motivation

Circuit analysis in the sinusoidal steady state is ordinary DC circuit analysis — Ohm's law, Kirchhoff's laws, series and parallel combination, nodal and mesh equations — carried out with complex numbers instead of real ones (Section 2.9). Every one of those operations reduces to the four arithmetic operations of this section. The arithmetic itself is easy; what takes practice is doing it **reliably** and choosing the form that minimises effort.

### The four operations

> [!theorem] Rectangular arithmetic
> For $z_1=a+jb$ and $z_2=c+jd$:
> $$
> \begin{aligned}
> z_1\pm z_2&=(a\pm c)+j(b\pm d),\\
> z_1z_2&=(ac-bd)+j(ad+bc),\\
> \frac{z_1}{z_2}&=\frac{z_1z_2^*}{z_2z_2^*}=\frac{(ac+bd)+j(bc-ad)}{c^2+d^2},\qquad z_2\neq0.
> \end{aligned}
> $$

**Addition** is component-wise — vector addition in the plane. **Multiplication** is ordinary expansion of brackets followed by $j^2=-1$:

$$
(a+jb)(c+jd)=ac+jad+jbc+j^2bd=(ac-bd)+j(ad+bc).
$$

**Division** uses the conjugate trick: multiplying numerator and denominator by $z_2^*$ makes the denominator the real number $|z_2|^2$, after which division is just scaling. In particular the **reciprocal** is

$$
\frac{1}{c+jd}=\frac{c-jd}{c^2+d^2},\qquad\text{e.g.}\quad \frac1j=\frac{-j}{1}=-j.
$$

> [!theorem] Polar arithmetic
> For $z_1=r_1\angle\theta_1$ and $z_2=r_2\angle\theta_2$:
> $$
> z_1z_2=r_1r_2\angle(\theta_1+\theta_2),\qquad \frac{z_1}{z_2}=\frac{r_1}{r_2}\angle(\theta_1-\theta_2),\qquad \frac{1}{z_2}=\frac{1}{r_2}\angle(-\theta_2).
> $$
> **Magnitudes multiply, angles add.**

*Proof.* Multiply the polar forms and use the compound-angle identities of Ch. 1:

$$
\begin{aligned}
z_1z_2&=r_1r_2\big[(\cos\theta_1\cos\theta_2-\sin\theta_1\sin\theta_2)+j(\sin\theta_1\cos\theta_2+\cos\theta_1\sin\theta_2)\big]\\
&=r_1r_2\big[\cos(\theta_1+\theta_2)+j\sin(\theta_1+\theta_2)\big].
\end{aligned}
$$

The quotient rule follows by checking that $\frac{r_1}{r_2}\angle(\theta_1-\theta_2)$ times $r_2\angle\theta_2$ gives $r_1\angle\theta_1$. $\square$

This proof is a little laborious; Section 2.5 gives a one-line version using exponentials. Notice that the trigonometric identities and the multiplication rule are *the same fact* in two notations — a theme Section 2.4 makes systematic.

### Powers of $j$

Because $j^2=-1$, the powers of $j$ cycle with period 4:

$$
j^0=1,\quad j^1=j,\quad j^2=-1,\quad j^3=-j,\quad j^4=1,\quad j^5=j,\ \dots
$$

so $j^n=j^{\,n\bmod 4}$. Negative powers follow from $1/j=-j$: $j^{-1}=-j$, $j^{-2}=-1$, $j^{-3}=j$. In polar form $j=1\angle90^\circ$, and $j^n=1\angle(90n)^\circ$ — each multiplication by $j$ is a quarter-turn.

### Field properties

> [!theorem] $\mathbb{C}$ is a field
> Addition and multiplication of complex numbers are commutative and associative, multiplication distributes over addition, $0$ and $1$ are identities, every $z$ has an additive inverse $-z$, and every $z\neq0$ has a multiplicative inverse $1/z$.

The practical meaning: **every algebraic manipulation you know for real numbers is valid for complex numbers** — expanding brackets, collecting terms, cross-multiplying, solving linear systems by elimination or Cramer's rule, factorising polynomials, the quadratic formula. What does *not* carry over are statements involving order ($<$, $>$) and some rules for multivalued functions such as square roots and logarithms (pitfall below, and Section 2.5).

### Magnitude rules and the triangle inequality

> [!theorem] Magnitude rules
> $$
> |z_1z_2|=|z_1|\,|z_2|,\qquad \left|\frac{z_1}{z_2}\right|=\frac{|z_1|}{|z_2|},\qquad |z^n|=|z|^n,\qquad |\operatorname{Re}\{z\}|\le|z|,\quad |\operatorname{Im}\{z\}|\le|z|.
> $$
> **Triangle inequality:**
> $$
> |z_1+z_2|\le|z_1|+|z_2|,\qquad\text{and}\qquad \big||z_1|-|z_2|\big|\le|z_1-z_2|.
> $$

*Proof of the product rule without polar form.* $|z_1z_2|^2=(z_1z_2)(z_1z_2)^*=z_1z_1^*\,z_2z_2^*=|z_1|^2|z_2|^2$.

*Proof of the triangle inequality.*

$$
\begin{aligned}
|z_1+z_2|^2&=(z_1+z_2)(z_1^*+z_2^*)=|z_1|^2+|z_2|^2+z_1z_2^*+z_1^*z_2\\
&=|z_1|^2+|z_2|^2+2\operatorname{Re}\{z_1z_2^*\}\\
&\le|z_1|^2+|z_2|^2+2|z_1z_2^*|=(|z_1|+|z_2|)^2.
\end{aligned}
$$

Take square roots. Equality holds exactly when $z_1z_2^*$ is real and non-negative, i.e. when $z_1$ and $z_2$ point in the same direction. The reverse inequality follows by applying the triangle inequality to $z_1=(z_1-z_2)+z_2$. $\square$

The middle line of that proof is worth remembering on its own:

$$
|z_1+z_2|^2=|z_1|^2+|z_2|^2+2|z_1||z_2|\cos(\theta_1-\theta_2),
$$

which is the **law of cosines**. When two sinusoids of the same frequency are added (Section 2.8), the amplitude of the sum is given by exactly this formula: maximal ($|z_1|+|z_2|$) when in phase, minimal ($\big||z_1|-|z_2|\big|$) when in antiphase.

> [!tip] Engineering Note — choose the form to suit the operation
> Rectangular for sums and differences (series impedances, KCL/KVL, adding phasors), polar for products, quotients, powers and roots (gains, voltage dividers, Ohm's law, cascaded stages). A mixed calculation such as $\dfrac{z_1z_2}{z_3+z_4}$ is done by forming $z_3+z_4$ in rectangular form, converting it to polar, then multiplying and dividing in polar. Numerical packages store complex numbers in rectangular form and convert on demand, so on a computer you simply write the expression; by hand, the choice of form is half the skill.

> [!example] Worked Example 2.9 — The four operations
> Let $z_1=3+j4$ and $z_2=1-j2$.
> - $z_1+z_2=4+j2$; $\quad z_1-z_2=2+j6$.
> - $z_1z_2=3-j6+j4-j^2 8=3+8+j(-6+4)=11-j2$.
> - $\dfrac{z_1}{z_2}=\dfrac{(3+j4)(1+j2)}{1^2+2^2}=\dfrac{3+j6+j4+j^2 8}{5}=\dfrac{-5+j10}{5}=-1+j2$.
>
> **Check the quotient** by multiplying back: $(-1+j2)(1-j2)=-1+j2+j2-j^24=3+j4$. ✓
>
> **Check with polar form:** $z_1=5\angle53.13^\circ$, $z_2=\sqrt5\angle{-63.43^\circ}$. Then $z_1z_2=5\sqrt5\angle(-10.30^\circ)=11.18\angle{-10.30^\circ}$, and $\lvert11-j2\rvert=\sqrt{125}=11.18$, $\operatorname{atan2}(-2,11)=-10.30^\circ$. ✓ Likewise $z_1/z_2=\sqrt5\angle116.57^\circ=-1+j2$. ✓

> [!example] Worked Example 2.10 — A compound expression
> Evaluate $w=\dfrac{(2+j)(3-j2)}{1+j}-j^7$ in rectangular form.
>
> **Step 1 — numerator.** $(2+j)(3-j2)=6-j4+j3-j^22=8-j$.
>
> **Step 2 — divide by $1+j$.** $\dfrac{8-j}{1+j}=\dfrac{(8-j)(1-j)}{2}=\dfrac{8-j8-j+j^2}{2}=\dfrac{7-j9}{2}=3.5-j4.5$.
>
> **Step 3 — the power of $j$.** $7\bmod4=3$, so $j^7=j^3=-j$, and $-j^7=+j$.
>
> **Step 4.** $w=3.5-j4.5+j=3.5-j3.5=3.5\sqrt2\angle{-45^\circ}\approx4.950\angle{-45^\circ}$.
>
> **Polar check of Step 2:** $\lvert8-j\rvert/\lvert1+j\rvert=\sqrt{65}/\sqrt2=5.701$ and $\lvert3.5-j4.5\rvert=\sqrt{32.5}=5.701$. ✓

> [!example] Worked Example 2.11 — A complex linear system
> Two coupled loop currents satisfy
> $$
> (3+j4)I_1-j4\,I_2=10,\qquad -j4\,I_1+(2+j2)I_2=0.
> $$
> (This is the form mesh analysis produces for two loops sharing an inductive element.) Solve for $I_1$ and $I_2$.
>
> **Step 1 — determinant.** Cramer's rule works exactly as for real systems:
> $$
> \Delta=(3+j4)(2+j2)-(-j4)(-j4)=(6+j6+j8-8)-(j^2\,16)=(-2+j14)+16=14+j14.
> $$
>
> **Step 2 — $I_1$.**
> $$
> I_1=\frac{10(2+j2)-(-j4)(0)}{\Delta}=\frac{20+j20}{14+j14}=\frac{20}{14}=1.4286\angle0^\circ.
> $$
>
> **Step 3 — $I_2$.**
> $$
> I_2=\frac{(3+j4)(0)-(-j4)(10)}{\Delta}=\frac{j40}{14+j14}=\frac{40\angle90^\circ}{19.799\angle45^\circ}=2.0203\angle45^\circ=1.4286+j1.4286.
> $$
>
> **Check** the second equation: $-j4(1.4286)+(2+j2)(1.4286+j1.4286)=-j5.714+(2.857-2.857)+j(2.857+2.857)=0$. ✓

> [!warning] Common Pitfall — $\sqrt{a}\sqrt{b}\neq\sqrt{ab}$ for negative numbers
> The rule $\sqrt a\sqrt b=\sqrt{ab}$ holds for $a,b\ge0$ but fails otherwise: $\sqrt{-4}\sqrt{-9}=(j2)(j3)=-6$, whereas $\sqrt{(-4)(-9)}=\sqrt{36}=6$. Square roots of complex numbers are two-valued (Section 2.6), and "laws" that silently pick one value can break. Safe practice: write $\sqrt{-a}=j\sqrt a$ for $a>0$ **first**, then multiply.

> [!tip] Engineering Note — numerical care
> Computing $|z|=\sqrt{a^2+b^2}$ naively can overflow or underflow when $a$ or $b$ is very large or very small, even if $|z|$ itself is representable. Library functions (`abs` for complex types, `hypot`) rescale internally. Similarly, robust complex division (Smith's algorithm) divides through by the larger of $|c|,|d|$ before forming $c^2+d^2$. You rarely need to implement these yourself, but you should know why `abs(z)` is preferred to `sqrt(z.real**2 + z.imag**2)`.

## 2.4 Euler's formula

### Motivation

Polar multiplication says angles **add** when numbers multiply. Exponentials have the same property: $e^{x}e^{y}=e^{x+y}$, exponents add when quantities multiply. Euler's formula makes this analogy exact — the angle of a complex number behaves precisely like an imaginary exponent. The payoff is enormous: sinusoids, which are awkward to differentiate, integrate, multiply and add, become exponentials, which are the easiest functions in all of calculus. This single identity underlies phasors, Fourier series and transforms, the Laplace transform, and the solution of every linear differential equation with constant coefficients.

> [!theorem] Euler's formula
> For every real $\theta$,
> $$
> e^{j\theta}=\cos\theta+j\sin\theta.
> $$

Before proving it, be clear about what is being claimed. For real $x$, $e^x$ is defined (Ch. 1, Ch. 4) and has the series $e^x=\sum_{n\ge0}x^n/n!$. For imaginary arguments $e^{j\theta}$ has no meaning until we *give* it one. There are three equivalent ways to do so, and each gives a different insight.

### Derivation 1 — power series

Define $e^{z}$ for any complex $z$ by the same series that works for real arguments,

$$
e^{z}=\sum_{n=0}^{\infty}\frac{z^n}{n!}=1+z+\frac{z^2}{2!}+\frac{z^3}{3!}+\cdots,
$$

which converges for every $z$ (Ch. 4, Ch. 18). Put $z=j\theta$ and use the cycle $j^2=-1$, $j^3=-j$, $j^4=1$:

$$
\begin{aligned}
e^{j\theta}&=1+j\theta+\frac{(j\theta)^2}{2!}+\frac{(j\theta)^3}{3!}+\frac{(j\theta)^4}{4!}+\frac{(j\theta)^5}{5!}+\cdots\\
&=1+j\theta-\frac{\theta^2}{2!}-j\frac{\theta^3}{3!}+\frac{\theta^4}{4!}+j\frac{\theta^5}{5!}-\cdots\\
&=\left(1-\frac{\theta^2}{2!}+\frac{\theta^4}{4!}-\cdots\right)+j\left(\theta-\frac{\theta^3}{3!}+\frac{\theta^5}{5!}-\cdots\right)\\
&=\cos\theta+j\sin\theta.
\end{aligned}
$$

The even powers are real and reproduce the Maclaurin series of $\cos\theta$; the odd powers are imaginary and reproduce that of $\sin\theta$. (Rearranging the terms is legitimate because the series converges absolutely.)

### Derivation 2 — a differential equation

Any sensible definition of $e^{j\theta}$ should satisfy $\frac{d}{d\theta}e^{j\theta}=je^{j\theta}$ and $e^{j0}=1$. Let $g(\theta)=\cos\theta+j\sin\theta$. Then

$$
g'(\theta)=-\sin\theta+j\cos\theta=j(\cos\theta+j\sin\theta)=j\,g(\theta),\qquad g(0)=1.
$$

So $g$ satisfies the same initial-value problem. To see that the two must coincide, consider $f(\theta)=e^{-j\theta}g(\theta)$. By the product rule $f'(\theta)=-je^{-j\theta}g+e^{-j\theta}jg=0$, so $f$ is constant, and $f(0)=1$. Hence $g(\theta)=e^{j\theta}$. (This uniqueness argument is the prototype of the uniqueness theorems for ODEs in Ch. 11.)

### Derivation 3 — geometry of motion

Think of $z(\theta)=e^{j\theta}$ as the position of a particle at "time" $\theta$. The equation $z'=jz$ says the velocity is the position rotated by $90^\circ$ — always perpendicular to the radius — and has the same magnitude as the position. Motion perpendicular to the radius keeps the distance from the origin constant, so the particle moves on the circle $|z|=1$; with unit speed, after time $\theta$ it has travelled arc length $\theta$, i.e. it sits at angle $\theta$ radians. That point is $(\cos\theta,\sin\theta)$. This picture explains *why* the angle in $e^{j\theta}$ must be in radians: radians are arc length on the unit circle.

### Immediate consequences

> [!theorem] Basic properties of $e^{j\theta}$
> For real $\theta,\phi$ and integer $k$:
> $$
> \begin{aligned}
> &|e^{j\theta}|=1,\qquad \angle e^{j\theta}=\theta,\qquad (e^{j\theta})^*=e^{-j\theta}=\frac{1}{e^{j\theta}},\\
> &e^{j\theta}e^{j\phi}=e^{j(\theta+\phi)},\qquad e^{j(\theta+2\pi k)}=e^{j\theta},\\
> &e^{j0}=1,\quad e^{j\pi/2}=j,\quad e^{j\pi}=-1,\quad e^{j3\pi/2}=e^{-j\pi/2}=-j,\quad e^{j2\pi}=1.
> \end{aligned}
> $$

The modulus follows from $\cos^2\theta+\sin^2\theta=1$. The addition rule is the polar multiplication theorem of Section 2.3 restated; periodicity follows from that of cosine and sine. The case $e^{j\pi}+1=0$ — linking $e$, $\pi$, $j$, $1$ and $0$ — is famous, but the engineer's version is more useful: **multiplying by $e^{j\pi}=-1$ is a $180^\circ$ phase shift, and multiplying by $e^{j\pi/2}=j$ is a $90^\circ$ phase advance.**

### The inverse Euler formulas

Replace $\theta$ by $-\theta$: $e^{-j\theta}=\cos\theta-j\sin\theta$. Adding and subtracting this from Euler's formula isolates the cosine and the sine.

> [!theorem] Cosine and sine as exponentials
> $$
> \cos\theta=\frac{e^{j\theta}+e^{-j\theta}}{2}=\operatorname{Re}\{e^{j\theta}\},\qquad \sin\theta=\frac{e^{j\theta}-e^{-j\theta}}{2j}=\operatorname{Im}\{e^{j\theta}\}.
> $$

These are the workhorses of signal analysis. They express the two real trigonometric functions as combinations of a **single** kind of function, the complex exponential, which differentiates, integrates and multiplies trivially:

$$
\frac{d}{d\theta}e^{ja\theta}=ja\,e^{ja\theta},\qquad \int e^{ja\theta}\,d\theta=\frac{e^{ja\theta}}{ja}+C\ (a\neq0),\qquad e^{ja\theta}e^{jb\theta}=e^{j(a+b)\theta}.
$$

The general strategy, used again and again in this book, is therefore:

1. Rewrite sines and cosines as exponentials (or as real/imaginary parts of exponentials).
2. Do the algebra or calculus with exponentials.
3. Convert back with Euler's formula, or take a real/imaginary part.

> [!tip] Engineering Note — rotating phasors and negative frequency
> As $t$ increases, $e^{j\omega t}$ traces the unit circle anticlockwise at $\omega$ rad/s, completing a revolution every $T=2\pi/\omega$ seconds. Its horizontal projection is $\cos\omega t$ and its vertical projection is $\sin\omega t$. The inverse Euler formula says a real cosine is the average of two **counter-rotating** exponentials, $e^{+j\omega t}$ and $e^{-j\omega t}$, whose imaginary parts cancel. This is the origin of the "negative frequencies" in two-sided spectra (Ch. 15–16): a real signal at $100$ Hz has spectral lines at both $+100$ Hz and $-100$ Hz, each carrying half the amplitude.

```desmos-graph
left=-1.6; right=1.6; bottom=-1.6; top=1.6
---
x^2+y^2=1|grey
(\cos(1),\sin(1))|blue|label:e^{j1}
(\cos(2.5),\sin(2.5))|red|label:e^{j2.5}
(\cos(-0.7),\sin(-0.7))|green|label:e^{-j0.7}
y=\tan(1)x|blue|dashed|0<x<0.5403
x=0.5403|blue|dotted|0<y<0.8415
```
*Euler's formula places $e^{j\theta}$ on the unit circle at angle $\theta$ radians. For $\theta=1$ the dashed radius and dotted drop line show the real part $\cos1=0.540$ and imaginary part $\sin1=0.841$.*

```desmos-graph
left=0; right=12.6; bottom=-1.4; top=1.4
---
y=\cos(x)|blue
y=\sin(x)|red
y=\cos(x-1.5708)|red|dotted
```
*The real part (blue, $\cos\theta$) and imaginary part (red, $\sin\theta$) of $e^{j\theta}$. The sine is the cosine delayed by $\pi/2$ — the dotted curve, $\cos(\theta-\pi/2)$, lies exactly on it — reflecting $\operatorname{Im}\{e^{j\theta}\}=\operatorname{Re}\{e^{j(\theta-\pi/2)}\}=\operatorname{Re}\{-je^{j\theta}\}$.*

### Trigonometry from exponentials

All the identities of Ch. 1, Section 1.6 can be derived mechanically, which means you never need to memorise more than Euler's formula.

> [!example] Worked Example 2.12 — Compound-angle formulas
> Derive $\cos(A+B)$ and $\sin(A+B)$.
>
> **Step 1.** $e^{j(A+B)}=e^{jA}e^{jB}$.
>
> **Step 2.** Expand the right side with Euler's formula:
> $$
> (\cos A+j\sin A)(\cos B+j\sin B)=(\cos A\cos B-\sin A\sin B)+j(\sin A\cos B+\cos A\sin B).
> $$
>
> **Step 3.** The left side is $\cos(A+B)+j\sin(A+B)$. Equate real and imaginary parts (one complex equation = two real ones):
> $$
> \cos(A+B)=\cos A\cos B-\sin A\sin B,\qquad \sin(A+B)=\sin A\cos B+\cos A\sin B.
> $$
>
> **Check** with $A=0.4$, $B=0.9$: $\cos1.3=0.2675$, and $\cos0.4\cos0.9-\sin0.4\sin0.9=0.9211(0.6216)-0.3894(0.7833)=0.5726-0.3050=0.2675$. ✓

> [!example] Worked Example 2.13 — Power reduction and the mixer identity
> **(a)** Express $\cos^3\theta$ in terms of cosines of multiple angles (no powers).
> $$
> \begin{aligned}
> \cos^3\theta&=\left(\frac{e^{j\theta}+e^{-j\theta}}{2}\right)^3=\frac18\left(e^{j3\theta}+3e^{j\theta}+3e^{-j\theta}+e^{-j3\theta}\right)\\
> &=\frac18\left(2\cos3\theta+6\cos\theta\right)=\frac34\cos\theta+\frac14\cos3\theta.
> \end{aligned}
> $$
> **Check** at $\theta=0.5$: $\cos^30.5=0.8776^3=0.6759$; $0.75(0.8776)+0.25\cos1.5=0.6582+0.25(0.0707)=0.6759$. ✓
>
> *Engineering reading:* passing a pure cosine through a cubic nonlinearity $y=x^3$ produces a fundamental component plus a **third harmonic** of one-third its size — harmonic distortion, computed in one line.
>
> **(b)** A mixer multiplies a received signal $\cos\omega_1t$ by a local oscillator $\cos\omega_2t$. Show what frequencies emerge.
> $$
> \cos\omega_1t\cos\omega_2t=\frac14\left(e^{j\omega_1t}+e^{-j\omega_1t}\right)\left(e^{j\omega_2t}+e^{-j\omega_2t}\right)=\frac12\cos(\omega_1-\omega_2)t+\frac12\cos(\omega_1+\omega_2)t.
> $$
> With $f_1=10.7$ MHz and $f_2=10.0$ MHz the product contains $0.7$ MHz and $20.7$ MHz components — the sum and difference frequencies on which every superheterodyne receiver relies.

> [!example] Worked Example 2.14 — Summing a trigonometric series
> Find a closed form for $C_N(\theta)=\displaystyle\sum_{k=0}^{N-1}\cos k\theta$ for $\theta$ not a multiple of $2\pi$.
>
> **Step 1 — embed in an exponential sum.** $C_N=\operatorname{Re}\{S_N\}$ where $S_N=\sum_{k=0}^{N-1}e^{jk\theta}$, a geometric series with ratio $e^{j\theta}\neq1$:
> $$
> S_N=\frac{1-e^{jN\theta}}{1-e^{j\theta}}.
> $$
>
> **Step 2 — factor out half-angles** (a standard trick that turns $1-e^{jx}$ into a sine):
> $$
> 1-e^{jx}=e^{jx/2}\left(e^{-jx/2}-e^{jx/2}\right)=-2j\,e^{jx/2}\sin\frac x2.
> $$
> Applying this to numerator and denominator:
> $$
> S_N=\frac{e^{jN\theta/2}\sin(N\theta/2)}{e^{j\theta/2}\sin(\theta/2)}=e^{j(N-1)\theta/2}\,\frac{\sin(N\theta/2)}{\sin(\theta/2)}.
> $$
>
> **Step 3 — take the real part.**
> $$
> \sum_{k=0}^{N-1}\cos k\theta=\cos\!\left(\frac{(N-1)\theta}{2}\right)\frac{\sin(N\theta/2)}{\sin(\theta/2)}.
> $$
>
> **Check** with $N=5$, $\theta=0.3$: the direct sum is $1+0.9553+0.8253+0.6216+0.3624=3.7646$; the formula gives $\cos(0.6)\,\sin(0.75)/\sin(0.15)=0.8253\times0.6816/0.1494=3.7646$. ✓
>
> The ratio $\sin(N\theta/2)/\sin(\theta/2)$ is the **Dirichlet kernel**; it describes the beam pattern of an $N$-element antenna array and the frequency response of an $N$-point moving-average filter.

> [!example] Worked Example 2.15 — An integral via exponentials
> Evaluate $I=\displaystyle\int_0^\infty e^{-2t}\cos3t\,dt$.
>
> **Step 1.** $\cos3t=\operatorname{Re}\{e^{j3t}\}$, and since $e^{-2t}$ is real and integration is linear, $I=\operatorname{Re}\left\{\int_0^\infty e^{(-2+j3)t}dt\right\}$.
>
> **Step 2.** For $a=-2+j3$, $\int_0^\infty e^{at}dt=\left[\frac{e^{at}}{a}\right]_0^\infty=0-\frac1a=\frac{1}{2-j3}$, because $|e^{at}|=e^{-2t}\to0$.
>
> **Step 3.** $\dfrac{1}{2-j3}=\dfrac{2+j3}{13}$, so $I=\dfrac{2}{13}\approx0.1538$.
>
> As a bonus the imaginary part gives $\int_0^\infty e^{-2t}\sin3t\,dt=\frac{3}{13}$ for free. Compare with the alternative — integrating by parts twice (Ch. 5). This calculation is, in fact, the Laplace transform of $\cos3t$ evaluated at $s=2$ (Ch. 14): $\mathcal{L}\{\cos\omega t\}=s/(s^2+\omega^2)=2/13$. ✓

> [!warning] Common Pitfall — degrees inside an exponent
> $e^{j30}$ means $e^{j\cdot30\ \text{rad}}$, which is $\cos30+j\sin30=0.154-j0.988$, not $0.866+j0.5$. The notation $1\angle30^\circ$ may use degrees; the exponential $e^{j\pi/6}$ must use radians. Calculators set to degree mode will happily give wrong answers to complex-exponential calculations.

> [!warning] Common Pitfall — $\cos\theta\neq e^{j\theta}$
> A cosine is the **real part** of $e^{j\theta}$, not equal to it. Replacing $\cos\omega t$ by $e^{j\omega t}$ in a calculation is legitimate only if every operation that follows is **linear with real coefficients** (addition, multiplication by real constants, differentiation, integration), so that "take the real part" can be postponed to the end. Products of two signals do not satisfy this: $\operatorname{Re}\{z_1\}\operatorname{Re}\{z_2\}\neq\operatorname{Re}\{z_1z_2\}$ in general. This is why power calculations need the conjugate formula $\tfrac12\operatorname{Re}\{\mathbf{V}\mathbf{I}^*\}$ rather than $\operatorname{Re}\{\mathbf{V}\mathbf{I}\}$ (Section 2.9).

## 2.5 Exponential representation

### Motivation

Combining the polar form of Section 2.1 with Euler's formula gives the most compact and most powerful way to write a complex number. In exponential form the rules for multiplication, division, powers and roots are simply the laws of exponents you already know, and a complex number can be fed directly into calculus. It is also the form in which signals are written: a damped oscillation, a growing instability and a pure tone are all $e^{st}$ for different complex $s$.

> [!definition] Exponential form
> Every non-zero complex number can be written
> $$
> z=re^{j\theta}=r(\cos\theta+j\sin\theta)=r\angle\theta,\qquad r=|z|>0,\ \theta=\arg z.
> $$

The notation $r\angle\theta$ is engineering shorthand for $re^{j\theta}$. Angles written with $\angle$ are often in degrees; inside $e^{j\theta}$ they must be in radians.

### Polar arithmetic from the laws of exponents

> [!theorem] Exponential-form arithmetic
> $$
> \begin{aligned}
> z_1z_2&=r_1r_2\,e^{j(\theta_1+\theta_2)}, &\frac{z_1}{z_2}&=\frac{r_1}{r_2}\,e^{j(\theta_1-\theta_2)},\\
> z^*&=r\,e^{-j\theta}, &\frac1z&=\frac1r\,e^{-j\theta},\qquad z^n=r^ne^{jn\theta}\ (n\in\mathbb{Z}).
> \end{aligned}
> $$

*Proof.* $z_1z_2=r_1e^{j\theta_1}r_2e^{j\theta_2}=r_1r_2e^{j\theta_1}e^{j\theta_2}=r_1r_2e^{j(\theta_1+\theta_2)}$, using $e^{j\theta_1}e^{j\theta_2}=e^{j(\theta_1+\theta_2)}$ from Section 2.4. The others follow similarly. $\square$

Compare this one-line proof with the trigonometric proof in Section 2.3: the exponential notation has absorbed the compound-angle identities.

An important special case: since $|e^{j\alpha}|=1$, multiplying by $e^{j\alpha}$ changes only the angle. **$e^{j\alpha}$ is a pure rotation (phase shift) by $\alpha$**, and a general multiplier $Ae^{j\alpha}$ is a scaling by $A$ combined with a rotation by $\alpha$. Every linear time-invariant system acting on a sinusoid does exactly this (Section 2.9).

> [!example] Worked Example 2.16 — Products and quotients in exponential form
> Evaluate $w=\dfrac{(2e^{j\pi/3})(6e^{-j\pi/6})}{3e^{j\pi/2}}$ and give it in rectangular form.
>
> **Magnitude:** $|w|=\dfrac{2\cdot6}{3}=4$.
>
> **Angle:** $\angle w=\dfrac\pi3-\dfrac\pi6-\dfrac\pi2=\dfrac{2\pi-\pi-3\pi}{6}=-\dfrac\pi3$.
>
> **Rectangular:** $w=4e^{-j\pi/3}=4\left(\tfrac12-j\tfrac{\sqrt3}{2}\right)=2-j2\sqrt3\approx2-j3.464$.
>
> **Check** by direct rectangular computation: $2e^{j\pi/3}=1+j1.732$, $6e^{-j\pi/6}=5.196-j3$, their product is $5.196+5.196+j(9-3)=10.392+j6$; dividing by $3e^{j\pi/2}=j3$ gives $(10.392+j6)/(j3)=2-j3.464$. ✓

> [!example] Worked Example 2.17 — A sum in exponential form
> Evaluate $w=4e^{j\pi/6}+3e^{-j\pi/3}$.
>
> Exponential form has no addition rule, so convert to rectangular:
> $$
> 4e^{j\pi/6}=3.464+j2.000,\qquad 3e^{-j\pi/3}=1.500-j2.598,\qquad w=4.964-j0.598.
> $$
> Then $|w|=\sqrt{4.964^2+0.598^2}=5.000$ and $\angle w=\operatorname{atan2}(-0.598,4.964)=-6.87^\circ$.
>
> **Geometric check:** the two terms point at $+30^\circ$ and $-60^\circ$, which are $90^\circ$ apart, so they form the legs of a right triangle: $|w|=\sqrt{4^2+3^2}=5$ ✓, and $w$ lies $\arctan\frac34=36.87^\circ$ clockwise of the $4e^{j\pi/6}$ term: $30^\circ-36.87^\circ=-6.87^\circ$. ✓

### The complex exponential function

For a general complex exponent $z=\sigma+j\omega$, the law $e^{a+b}=e^ae^b$ (which holds for the series definition) gives

> [!theorem] Exponential of a complex number
> $$
> e^{\sigma+j\omega}=e^{\sigma}e^{j\omega}=e^{\sigma}(\cos\omega+j\sin\omega),\qquad |e^{z}|=e^{\operatorname{Re}\{z\}},\quad \angle e^{z}=\operatorname{Im}\{z\}.
> $$
> Moreover $e^z\neq0$ for every $z$, and $e^{z+j2\pi k}=e^z$ for every integer $k$: the complex exponential is **periodic with period $j2\pi$**.

So the real part of the exponent controls the **size** and the imaginary part the **angle**. This splitting is the key to the Laplace variable $s=\sigma+j\omega$ ([Ch. 14](../part-5-transforms/14-laplace-transforms.md)). Replace $\sigma,\omega$ by $\sigma t,\omega t$:

$$
e^{st}=e^{\sigma t}e^{j\omega t},\qquad \operatorname{Re}\{e^{st}\}=e^{\sigma t}\cos\omega t.
$$

The single expression $e^{st}$ covers every elementary signal shape:

| $s$ | $\operatorname{Re}\{e^{st}\}$ | Behaviour |
|---|---|---|
| $0$ | $1$ | constant (DC) |
| $\sigma<0$ real | $e^{\sigma t}$ | exponential decay |
| $\sigma>0$ real | $e^{\sigma t}$ | exponential growth |
| $j\omega$ | $\cos\omega t$ | steady oscillation |
| $\sigma+j\omega$, $\sigma<0$ | $e^{\sigma t}\cos\omega t$ | damped oscillation |
| $\sigma+j\omega$, $\sigma>0$ | $e^{\sigma t}\cos\omega t$ | growing oscillation (instability) |

More generally, a damped sinusoid with arbitrary amplitude and phase is

$$
A\,e^{\sigma t}\cos(\omega t+\phi)=\operatorname{Re}\left\{Ae^{j\phi}\,e^{(\sigma+j\omega)t}\right\},
$$

a complex constant $Ae^{j\phi}$ (which will be called a phasor in Section 2.8) times a complex exponential.

```desmos-graph
left=0; right=10; bottom=-1.2; top=1.2
---
y=e^{-0.3x}\cos(3x)|blue
y=e^{-0.3x}|blue|dashed
y=-e^{-0.3x}|blue|dashed
y=0.15e^{0.2x}\cos(3x)|red
y=0.15e^{0.2x}|red|dotted
y=-0.15e^{0.2x}|red|dotted
```
*$\operatorname{Re}\{e^{st}\}$ for $s=-0.3+j3$ (blue, decaying inside its envelope $\pm e^{-0.3t}$) and $0.15\operatorname{Re}\{e^{st}\}$ for $s=0.2+j3$ (red, growing inside $\pm0.15e^{0.2t}$). The imaginary part of $s$ sets the oscillation rate; the real part sets the envelope.*

> [!example] Worked Example 2.18 — Evaluating $e^z$ and writing a signal as $\operatorname{Re}\{Ae^{st}\}$
> **(a)** Evaluate $e^{1+j2}$.
> $$
> e^{1+j2}=e^1(\cos2+j\sin2)=2.7183(-0.4161+j0.9093)=-1.1312+j2.4717.
> $$
> Its magnitude is $e^1=2.718$ and its angle $2$ rad $=114.6^\circ$ — read straight off the exponent.
>
> **(b)** Write $x(t)=3e^{-2t}\cos(5t-40^\circ)$ in the form $\operatorname{Re}\{Ae^{st}\}$.
>
> Match to $\operatorname{Re}\{Ae^{j\phi}e^{(\sigma+j\omega)t}\}$: $A=3\angle{-40^\circ}=2.298-j1.928$ and $s=-2+j5\ \mathrm{s^{-1}}$. So
> $$
> x(t)=\operatorname{Re}\left\{(2.298-j1.928)\,e^{(-2+j5)t}\right\}.
> $$
> **Check** at $t=0$: $\operatorname{Re}\{A\}=2.298$ and $3\cos(-40^\circ)=2.298$. ✓
>
> **(c)** Solve $e^z=-2$. Write $z=x+jy$; then $e^x e^{jy}=2e^{j\pi}$. Equating magnitudes, $e^x=2$, so $x=\ln2$; equating angles (up to $2\pi$), $y=\pi+2\pi k$. Hence $z=\ln2+j(2k+1)\pi$, $k\in\mathbb{Z}$ — infinitely many solutions, because $e^z$ is periodic.

### The complex logarithm

Part (c) shows that inverting the exponential is ambiguous. Solving $e^w=z$ for $z=re^{j\theta}\neq0$ in general gives $w=\ln r+j(\theta+2\pi k)$.

> [!definition] Complex logarithm
> For $z\neq0$, the **logarithm** is the multivalued expression
> $$
> \log z=\ln|z|+j(\operatorname{Arg}z+2\pi k),\qquad k\in\mathbb{Z},
> $$
> and the **principal logarithm** is the single value
> $$
> \operatorname{Log}z=\ln|z|+j\operatorname{Arg}z,\qquad -\pi<\operatorname{Arg}z\le\pi.
> $$

The real part of a logarithm is the log of the magnitude; the imaginary part is the angle. In engineering this is exactly the content of a Bode plot: $\ln H=\ln|H|+j\angle H$, so plotting $\ln|H|$ (or $20\log_{10}|H|$, in dB) and $\angle H$ is plotting the real and imaginary parts of $\ln H$, and the gains and phases of cascaded stages add because logarithms turn products into sums.

### Complex powers

For $z\neq0$ and complex $c$, define $z^c=e^{c\log z}$. This is multivalued in general (single-valued when $c$ is an integer, $n$-valued when $c=1/n$ — the $n$-th roots of Section 2.6). The **principal value** uses $\operatorname{Log}z$.

> [!example] Worked Example 2.19 — Logarithms and complex powers
> **(a)** $\operatorname{Log}(-3+j4)$. Here $|z|=5$ and $\operatorname{Arg}z=126.87^\circ=2.2143$ rad, so
> $$
> \operatorname{Log}(-3+j4)=\ln5+j2.2143=1.6094+j2.2143.
> $$
> **Check:** $e^{1.6094}e^{j2.2143}=5(\cos2.2143+j\sin2.2143)=5(-0.6+j0.8)=-3+j4$. ✓
>
> **(b)** $j^j$. Since $\log j=j(\pi/2+2\pi k)$,
> $$
> j^j=e^{j\log j}=e^{j\cdot j(\pi/2+2\pi k)}=e^{-\pi/2-2\pi k}.
> $$
> Every value is **real**; the principal value ($k=0$) is $e^{-\pi/2}=0.2079$.
>
> **(c)** Principal value of $(1+j)^j$. $\operatorname{Log}(1+j)=\ln\sqrt2+j\pi/4=0.3466+j0.7854$. Then
> $$
> (1+j)^j=e^{j(0.3466+j0.7854)}=e^{-0.7854}e^{j0.3466}=0.4559(0.9405+j0.3397)=0.4288+j0.1549.
> $$

> [!warning] Common Pitfall — logarithm laws with complex numbers
> The rule $\operatorname{Log}(z_1z_2)=\operatorname{Log}z_1+\operatorname{Log}z_2$ can fail for principal values, because the sum of two principal arguments may fall outside $(-\pi,\pi]$. Example: $z_1=z_2=-1$ gives $\operatorname{Log}(1)=0$ but $\operatorname{Log}(-1)+\operatorname{Log}(-1)=j2\pi$. The identity holds **modulo $j2\pi$**. In phase calculations this appears as the familiar "phase wrapping": a phase computed by adding contributions may differ from the instrument's principal value by a multiple of $360^\circ$, and both are correct.

### Hyperbolic functions and complex arguments

Replacing $\theta$ by $jx$ in the inverse Euler formulas links trigonometric and hyperbolic functions:

$$
\cos(jx)=\frac{e^{-x}+e^{x}}{2}=\cosh x,\qquad \sin(jx)=\frac{e^{-x}-e^{x}}{2j}=j\sinh x.
$$

Combined with the compound-angle formula this gives the real and imaginary parts of the sine and cosine of any complex number:

$$
\sin(x+jy)=\sin x\cosh y+j\cos x\sinh y,\qquad \cos(x+jy)=\cos x\cosh y-j\sin x\sinh y.
$$

For example $\sin(\pi/2+j1)=\cosh1=1.5431$ — a sine can exceed 1 when its argument is complex. These formulas appear in transmission-line theory, where the voltage along a lossy line involves $\cosh(\gamma\ell)$ and $\sinh(\gamma\ell)$ with a complex propagation constant $\gamma=\alpha+j\beta$, and in waveguides below cut-off, where a "sinusoidal" field becomes an exponentially decaying one.

> [!tip] Engineering Note — the exponential form in software
> Numerical libraries provide `exp`, `log`, `abs` and `angle`/`arg` for complex arguments, and they return **principal** values. When a quantity is built up by multiplying many complex factors (a cascade of filter stages, a long chain of transmission-line sections), it is often more robust to accumulate $\ln|H|$ and $\angle H$ separately — adding logarithms of magnitudes and adding (unwrapped) phases — than to multiply complex numbers directly, which can overflow or lose track of full turns.

## 2.6 Roots and powers of complex numbers

### Motivation

Raising a complex number to a high power by repeated rectangular multiplication is tedious and error-prone; in exponential form it is a single step. The inverse problem — finding all solutions of $w^n=z$ — is even more important in engineering, because the poles of standard filter families, the frequencies of a discrete Fourier transform, the phases of a polyphase power system and the symbols of a phase-shift-keyed modulator are all $n$-th roots of some complex number. The central fact is that a non-zero complex number has exactly $n$ distinct $n$-th roots, arranged symmetrically on a circle.

### Integer powers: De Moivre's theorem

> [!theorem] De Moivre's theorem
> For every integer $n$ (positive, zero or negative) and real $\theta$,
> $$
> (\cos\theta+j\sin\theta)^n=\cos n\theta+j\sin n\theta,\qquad\text{equivalently}\qquad (re^{j\theta})^n=r^ne^{jn\theta}.
> $$

*Proof.* For $n\ge0$, by induction: the case $n=0$ is $1=1$. If true for $n$, then

$$
(\cos\theta+j\sin\theta)^{n+1}=(\cos n\theta+j\sin n\theta)(\cos\theta+j\sin\theta)=\cos(n+1)\theta+j\sin(n+1)\theta
$$

by the polar multiplication rule (angles add). For negative $n=-m$, $(\cos\theta+j\sin\theta)^{-m}=1/(\cos m\theta+j\sin m\theta)=\cos m\theta-j\sin m\theta=\cos(-m\theta)+j\sin(-m\theta)$, using $1/e^{j\phi}=e^{-j\phi}$. With exponentials the whole theorem is just $(e^{j\theta})^n=e^{jn\theta}$. $\square$

**Recipe for $z^n$:** convert to polar, raise the magnitude to the $n$-th power, multiply the angle by $n$, reduce the angle to the principal range if desired, convert back.

> [!example] Worked Example 2.20 — High powers
> **(a)** $(1+j)^8$. Here $1+j=\sqrt2\,e^{j\pi/4}$, so $(1+j)^8=(\sqrt2)^8e^{j2\pi}=16$.
>
> *Check:* $(1+j)^2=2j$, so $(1+j)^8=(2j)^4=16j^4=16$. ✓
>
> **(b)** $(1+j\sqrt3)^{10}$. Here $1+j\sqrt3=2e^{j\pi/3}$, so
> $$
> (1+j\sqrt3)^{10}=2^{10}e^{j10\pi/3}=1024\,e^{j(10\pi/3-4\pi)}=1024\,e^{-j2\pi/3}=1024\left(-\tfrac12-j\tfrac{\sqrt3}{2}\right)=-512-j886.8.
> $$
> Expanding by the binomial theorem would have required eleven terms.
>
> **(c)** A negative power: $(1-j\sqrt3)^{-4}$. With $1-j\sqrt3=2e^{-j\pi/3}$,
> $$
> (1-j\sqrt3)^{-4}=2^{-4}e^{j4\pi/3}=\tfrac{1}{16}\left(-\tfrac12-j\tfrac{\sqrt3}{2}\right)=-0.03125-j0.05413.
> $$

### Multiple-angle formulas

Expanding the left side of De Moivre's theorem with the binomial theorem and equating real and imaginary parts produces formulas for $\cos n\theta$ and $\sin n\theta$ as polynomials in $\cos\theta$ and $\sin\theta$ — the reverse of the power-reduction of Example 2.13.

> [!example] Worked Example 2.21 — Triple-angle formulas
> Write $c=\cos\theta$, $s=\sin\theta$. Then
> $$
> (c+js)^3=c^3+3c^2(js)+3c(js)^2+(js)^3=(c^3-3cs^2)+j(3c^2s-s^3).
> $$
> By De Moivre this equals $\cos3\theta+j\sin3\theta$. Equate parts and use $s^2=1-c^2$, $c^2=1-s^2$:
> $$
> \cos3\theta=c^3-3c(1-c^2)=4\cos^3\theta-3\cos\theta,\qquad \sin3\theta=3(1-s^2)s-s^3=3\sin\theta-4\sin^3\theta.
> $$
> **Check** at $\theta=0.7$: $4\cos^30.7-3\cos0.7=-0.5048$ and $\cos2.1=-0.5048$. ✓
>
> The polynomial $T_3(x)=4x^3-3x$ with $\cos3\theta=T_3(\cos\theta)$ is a **Chebyshev polynomial**; the family $T_n$ defined by $\cos n\theta=T_n(\cos\theta)$ is the basis of Chebyshev filter design.

### $n$-th roots

> [!theorem] $n$-th roots of a complex number
> Let $z=re^{j\theta}\neq0$ and $n\ge1$ an integer. The equation $w^n=z$ has exactly $n$ distinct solutions
> $$
> w_k=r^{1/n}\exp\!\left(j\,\frac{\theta+2\pi k}{n}\right),\qquad k=0,1,\dots,n-1,
> $$
> where $r^{1/n}$ is the ordinary positive real root. They lie on the circle $|w|=r^{1/n}$, equally spaced at angular intervals of $2\pi/n$, forming the vertices of a regular $n$-gon.

*Proof.* Write the unknown as $w=\rho e^{j\phi}$. Then $w^n=\rho^ne^{jn\phi}$, and $w^n=z$ requires equal magnitudes and angles equal up to multiples of $2\pi$:

$$
\rho^n=r\ \Rightarrow\ \rho=r^{1/n};\qquad n\phi=\theta+2\pi k\ \Rightarrow\ \phi=\frac{\theta+2\pi k}{n},\quad k\in\mathbb{Z}.
$$

Values of $k$ that differ by a multiple of $n$ give angles differing by a multiple of $2\pi$ — the same point. So exactly the $n$ values $k=0,\dots,n-1$ give distinct roots. $\square$

The crucial step is keeping the "$+2\pi k$" **before** dividing by $n$. Adding $2\pi$ to $\theta$ does not change $z$, but adding $2\pi/n$ to the root's angle does change the root — that is where the extra roots come from. The root with $k=0$ and principal $\theta$ is the **principal root**.

**Recipe for the $n$-th roots of $z$:** (1) write $z=re^{j\theta}$; (2) the first root is $r^{1/n}e^{j\theta/n}$; (3) the others follow by repeatedly rotating by $2\pi/n$ (i.e. multiplying by $e^{j2\pi/n}$).

> [!example] Worked Example 2.22 — Cube roots of $8j$
> Solve $w^3=8j$.
>
> **Step 1.** $8j=8e^{j\pi/2}$, so $r^{1/3}=2$ and the angles are $\frac{\pi/2+2\pi k}{3}=\frac\pi6+\frac{2\pi k}{3}$.
>
> **Step 2.**
> $$
> w_0=2e^{j\pi/6}=\sqrt3+j,\qquad w_1=2e^{j5\pi/6}=-\sqrt3+j,\qquad w_2=2e^{j3\pi/2}=-2j.
> $$
>
> **Check** $w_2$: $(-2j)^3=-8j^3=-8(-j)=8j$. ✓ And $w_0$: $(\sqrt3+j)^2=3+2\sqrt3j-1=2+j2\sqrt3$; times $(\sqrt3+j)$ gives $2\sqrt3+j2+j6-2\sqrt3=8j$. ✓
>
> **Structural check:** the roots of $w^3-8j=0$ sum to zero (the coefficient of $w^2$ is zero): $\sqrt3+j-\sqrt3+j-2j=0$. ✓

```desmos-graph
left=-3; right=3; bottom=-2.6; top=2.6
---
x^2+y^2=4|grey|dashed
(1.7321,1)|blue|label:w_0
(-1.7321,1)|blue|label:w_1
(0,-2)|blue|label:w_2
y=\frac{1}{\sqrt{3}}x|blue|dotted|0<x<1.7321
y=-\frac{1}{\sqrt{3}}x|blue|dotted|-1.7321<x<0
x=0|blue|dotted|-2<y<0
```
*The three cube roots of $8j$ lie on the circle of radius $2$, $120^\circ$ apart, starting at $30^\circ$ (one-third of the angle $90^\circ$ of $8j$).*

### Square roots in rectangular form

For square roots it is sometimes quicker to avoid polar form. Seek $w=x+jy$ with $w^2=a+jb$: equating parts of $x^2-y^2+j2xy=a+jb$ gives $x^2-y^2=a$, $2xy=b$, and taking magnitudes gives $x^2+y^2=r=\sqrt{a^2+b^2}$. Adding and subtracting:

> [!theorem] Square root of $a+jb$
> $$
> \sqrt{a+jb}=\pm\left(\sqrt{\frac{r+a}{2}}+j\,\operatorname{sgn}(b)\sqrt{\frac{r-a}{2}}\right),\qquad r=\sqrt{a^2+b^2},
> $$
> with $\operatorname{sgn}(b)=1$ if $b\ge0$ and $-1$ if $b<0$ (the sign makes $2xy$ have the sign of $b$).

> [!example] Worked Example 2.23 — A quadratic with complex coefficients
> Solve $z^2-(3+j)z+(4+j3)=0$.
>
> **Step 1 — discriminant.** The quadratic formula is valid in $\mathbb{C}$ (it is pure field algebra):
> $$
> D=(3+j)^2-4(4+j3)=(9+j6-1)-16-j12=-8-j6.
> $$
>
> **Step 2 — square root of $D$.** $a=-8$, $b=-6$, $r=10$: $\sqrt{(10-8)/2}=1$, $\sqrt{(10+8)/2}=3$, and $b<0$, so $\sqrt D=\pm(1-j3)$.
> *Check:* $(1-j3)^2=1-j6-9=-8-j6$. ✓
>
> **Step 3 — roots.**
> $$
> z=\frac{(3+j)\pm(1-j3)}{2}:\qquad z_1=\frac{4-j2}{2}=2-j,\qquad z_2=\frac{2+j4}{2}=1+j2.
> $$
>
> **Check with Vieta:** $z_1+z_2=3+j$ ✓ and $z_1z_2=(2-j)(1+j2)=2+j4-j+2=4+j3$. ✓
>
> Note that the roots are **not** conjugates of each other: the conjugate root theorem needs real coefficients.

### Roots of unity

> [!definition] Roots of unity
> The $n$ solutions of $w^n=1$ are the **$n$-th roots of unity**
> $$
> W^k=e^{j2\pi k/n},\qquad k=0,1,\dots,n-1,\qquad W=e^{j2\pi/n}.
> $$
> $W$ is the **primitive** $n$-th root of unity; the others are its powers.

The $n$-th roots of any $z$ are the principal root multiplied by each root of unity: $w_k=w_0W^k$. Roots of unity therefore carry all the "rotational" structure of root-finding.

> [!theorem] Sums of roots of unity
> For $n\ge2$,
> $$
> \sum_{k=0}^{n-1}W^k=0,\qquad\text{and more generally}\qquad \sum_{k=0}^{n-1}W^{mk}=\begin{cases}n, & m\equiv0\pmod n,\\ 0, & \text{otherwise}.\end{cases}
> $$

*Proof.* If $m$ is a multiple of $n$, every term is $W^{mk}=e^{j2\pi mk/n}=1$ and the sum is $n$. Otherwise $W^m\ne1$ and the geometric series gives $\dfrac{1-W^{mn}}{1-W^m}=\dfrac{1-1}{1-W^m}=0$, since $W^{mn}=(W^n)^m=1$. $\square$

Geometrically: the roots are the vertices of a regular polygon centred at the origin, and arrows from the centre to the vertices of a regular polygon balance out. The second form is the **orthogonality relation** of the discrete Fourier transform ([Ch. 16](../part-5-transforms/16-fourier-transforms.md), [Ch. 24](../part-7-numerical-computational/24-computational-mathematics.md)), where $W_N=e^{-j2\pi/N}$ is called the twiddle factor; it is what allows the DFT to be inverted.

```desmos-graph
left=-1.6; right=1.6; bottom=-1.4; top=1.4
---
x^2+y^2=1|grey|dashed
(1,0)|blue|label:k=0
(\cos(1.0472),\sin(1.0472))|blue|label:k=1
(\cos(2.0944),\sin(2.0944))|blue|label:k=2
(-1,0)|blue|label:k=3
(\cos(4.1888),\sin(4.1888))|blue|label:k=4
(\cos(5.236),\sin(5.236))|blue|label:k=5
```
*The six sixth roots of unity, $e^{j2\pi k/6}$, form a regular hexagon on the unit circle; their sum is zero.*

> [!tip] Engineering Note — three-phase power and PSK
> The cube roots of unity $1$, $a=e^{j2\pi/3}$, $a^2=e^{-j2\pi/3}$ describe the $120^\circ$ spacing of a balanced three-phase supply. Because $1+a+a^2=0$, the three phase voltages (or balanced load currents) sum to zero at every instant, so the neutral conductor of a balanced system carries no current. Power engineers use the operator "$a$" throughout symmetrical-component analysis of faults. In digital communications, $M$-ary phase-shift keying transmits one of the $M$-th roots of unity (scaled) per symbol: QPSK uses $e^{j(\pi/4+k\pi/2)}$, the four fourth roots of $-1$.

### Factoring polynomials with roots

> [!example] Worked Example 2.24 — Factoring $z^4+1$ over the reals
> **Step 1 — roots.** $z^4=-1=e^{j\pi}$ gives $z_k=e^{j(\pi+2\pi k)/4}$, i.e. $e^{\pm j\pi/4}$ and $e^{\pm j3\pi/4}$, or $\frac{\pm1\pm j}{\sqrt2}$.
>
> **Step 2 — pair conjugates.** For a pair $e^{\pm j\alpha}$ the real quadratic factor is $z^2-2\cos\alpha\,z+1$:
> $$
> \alpha=\tfrac\pi4:\ z^2-\sqrt2z+1,\qquad \alpha=\tfrac{3\pi}4:\ z^2+\sqrt2z+1.
> $$
>
> **Step 3.** $z^4+1=(z^2-\sqrt2z+1)(z^2+\sqrt2z+1)$.
> *Check:* the product is $(z^2+1)^2-2z^2=z^4+2z^2+1-2z^2=z^4+1$. ✓
>
> Such factorisations are needed when integrating rational functions by partial fractions (Ch. 5) and when building filters from second-order sections.

> [!example] Worked Example 2.25 — Butterworth filter poles
> An $n$-th order Butterworth low-pass filter with cut-off $1$ rad/s has $|H(j\omega)|^2=\dfrac{1}{1+\omega^{2n}}$. Its poles are the left-half-plane roots of $1+(-s^2)^n=0$, obtained by substituting $\omega^2=-s^2$ (since $s=j\omega$). Find the third-order Butterworth polynomial.
>
> **Step 1.** For $n=3$: $1+(-s^2)^3=1-s^6=0$, so $s^6=1$: the sixth roots of unity, $s_k=e^{jk\pi/3}$.
>
> **Step 2 — keep the stable ones.** The roots with $\operatorname{Re}\{s\}<0$ are $k=2,3,4$: $e^{j2\pi/3}=-\tfrac12+j\tfrac{\sqrt3}{2}$, $-1$, and $e^{j4\pi/3}=-\tfrac12-j\tfrac{\sqrt3}{2}$. (The right-half-plane roots belong to $H(-s)$, the mirror factor of $|H|^2=H(s)H(-s)$.)
>
> **Step 3 — multiply out.** The conjugate pair gives $s^2+s+1$ (since $-2\cos(2\pi/3)=1$), and
> $$
> B_3(s)=(s+1)(s^2+s+1)=s^3+2s^2+2s+1,\qquad H(s)=\frac{1}{s^3+2s^2+2s+1}.
> $$
>
> **Check** at $\omega=1$: $B_3(j)=-j-2+2j+1=-1+j$, so $|H(j1)|^2=1/|{-1+j}|^2=1/2$, matching $1/(1+1^6)$. ✓ The poles lie on the unit circle at $120^\circ$, $180^\circ$ and $240^\circ$ — evenly spaced by $\pi/n=60^\circ$, the defining geometry of the Butterworth family.

> [!warning] Common Pitfall — keeping only the principal root
> A calculator's $\sqrt[n]{z}$, or `z**(1/n)` in software, returns only the principal root ($k=0$). Equations such as $w^4=-16$ have four roots ($\sqrt2(\pm1\pm j)$); forgetting the other three loses poles of filters, branches of a root locus, or valid solutions of a characteristic equation. Always write "$+2\pi k$" before dividing the angle by $n$.

> [!warning] Common Pitfall — $(z^{1/n})^n$ versus $(z^n)^{1/n}$
> $(z^{1/n})^n=z$ for every root, but $(z^n)^{1/n}$ need not return $z$: $\big((-1)^2\big)^{1/2}=1^{1/2}$ has principal value $1$, not $-1$. Fractional powers do not commute with integer powers once principal values are used, for the same reason the logarithm laws fail (Section 2.5).

## 2.7 Complex plane and geometric interpretation

### Motivation

Many engineering judgements are made by looking at the complex plane rather than by computing: whether a system is stable (where are its poles?), how much phase margin a feedback loop has (how close does the Nyquist curve pass to $-1$?), whether an antenna is matched (where does its reflection coefficient sit on the Smith chart?), how noisy a received constellation is (how spread are the points?). Each of these is a geometric statement about complex numbers. This section collects the dictionary between algebra and geometry so that such pictures can be read fluently.

### Complex numbers as vectors

> [!definition] Complex (Argand) plane
> The **complex plane** identifies $z=x+jy$ with the point $(x,y)$, or with the vector from the origin to that point. The horizontal axis is the **real axis**, the vertical axis the **imaginary axis**. The magnitude $|z|$ is the length of the vector and $\arg z$ its direction.

Under this identification addition, subtraction and real scaling of complex numbers are exactly vector addition, subtraction and scaling in $\mathbb{R}^2$. What the plane adds to $\mathbb{R}^2$ is **multiplication**, which has no counterpart for ordinary 2-D vectors and which is a rotation-plus-scaling.

> [!theorem] Geometric meaning of the operations
> | Operation | Geometric effect |
> |---|---|
> | $z+w$ | translation of $z$ by the vector $w$ (parallelogram rule) |
> | $z-w$ | vector from $w$ to $z$ |
> | $\lvert z-w\rvert$ | distance between the points $z$ and $w$ |
> | $\arg(z-w)$ | direction of the segment from $w$ to $z$ |
> | $kz$, $k>0$ real | scaling by $k$ about the origin |
> | $-z$ | rotation by $180^\circ$ (point reflection in the origin) |
> | $jz$ | rotation by $+90^\circ$ (anticlockwise) |
> | $e^{j\alpha}z$ | rotation by $\alpha$ about the origin |
> | $(re^{j\alpha})z$ | scaling by $r$ combined with rotation by $\alpha$ (a "spiral similarity") |
> | $z^*$ | reflection in the real axis |
> | $1/z$ | $\frac1r e^{-j\theta}$: inversion in the unit circle ($r\to1/r$) followed by reflection in the real axis |
> | $w_0+e^{j\alpha}(z-w_0)$ | rotation by $\alpha$ about the point $w_0$ |

The last entry is built from simpler ones: translate $w_0$ to the origin, rotate, translate back. Composite geometric transformations are obtained by composing complex functions.

*Why multiplication by $j$ is a quarter-turn, without polar form:* $j(x+jy)=-y+jx$, and the map $(x,y)\mapsto(-y,x)$ is the standard $90^\circ$ rotation matrix $\begin{pmatrix}0&-1\\1&0\end{pmatrix}$ of linear algebra (Ch. 9). More generally, multiplication by $a+jb$ is the linear map with matrix $\begin{pmatrix}a&-b\\b&a\end{pmatrix}$, which is $r$ times a rotation matrix.

```desmos-graph
left=-6; right=6; bottom=-2; top=6
---
(3,4)|blue|label:z=3+j4
(-4,3)|red|label:jz=-4+j3
(1.5,2)|green|label:0.5z
y=\frac{4}{3}x|blue|0<x<3
y=-\frac{3}{4}x|red|-4<x<0
x^2+y^2=25|grey|dashed
```
*Multiplying by $j$ rotates $z=3+j4$ by $90^\circ$ to $-4+j3$; both lie on the circle $\lvert z\rvert=5$. Multiplying by the real number $0.5$ only shrinks it.*

> [!example] Worked Example 2.26 — Completing an equilateral triangle
> Two vertices of an equilateral triangle are $0$ and $4+j2$. Find the third vertex.
>
> **Idea.** The third vertex is obtained by rotating the side vector $4+j2$ by $\pm60^\circ$ about $0$.
>
> **Step 1.** $e^{j\pi/3}=\tfrac12+j\tfrac{\sqrt3}{2}=0.5+j0.8660$.
>
> **Step 2.** $(4+j2)(0.5+j0.866)=2+j3.464+j1-1.732=0.268+j4.464$.
>
> **Step 3.** The other orientation: $(4+j2)(0.5-j0.866)=2-j3.464+j1+1.732=3.732-j2.464$.
>
> **Check** the first solution: its distance from $0$ is $\sqrt{0.268^2+4.464^2}=4.472=\sqrt{20}=|4+j2|$, and its distance from $4+j2$ is $|{-3.732}+j2.464|=\sqrt{13.93+6.07}=\sqrt{20}$. ✓ All three sides equal.
>
> The same idea computes the line-to-line voltages of a three-phase system (Section 2.8), whose phasors form an equilateral triangle.

### Distances, midpoints and centroids

Distances between points are magnitudes of differences, $|z_1-z_2|$. The midpoint of the segment from $z_1$ to $z_2$ is $\frac12(z_1+z_2)$, and the point dividing it in ratio $t:(1-t)$ is $(1-t)z_1+tz_2$, $0\le t\le1$ — the **parametric equation of a segment**. The centroid of points $z_1,\dots,z_n$ is their average $\frac1n\sum z_k$; for instance the centroid of the triangle $1+j2$, $5-j$, $3+j4$ is $\frac13(9+j5)=3+j1.667$. In a received constellation, the centroid of many noisy samples estimates the transmitted symbol.

### Curves and regions

Simple equations and inequalities in $z$ describe the curves and regions that recur in engineering. Here $z=x+jy$.

> [!theorem] Standard loci
> | Equation | Locus |
> |---|---|
> | $\lvert z-z_0\rvert=R$ | circle, centre $z_0$, radius $R$ |
> | $\lvert z-z_0\rvert<R$ | open disc (interior of the circle) |
> | $\lvert z-z_1\rvert=\lvert z-z_2\rvert$ | perpendicular bisector of the segment $z_1z_2$ |
> | $\lvert z-z_1\rvert=k\lvert z-z_2\rvert$, $k>0$, $k\neq1$ | a circle (circle of Apollonius) |
> | $\arg(z-z_0)=\alpha$ | ray from $z_0$ at angle $\alpha$ (excluding $z_0$) |
> | $\operatorname{Re}\{z\}=c$, $\operatorname{Im}\{z\}=c$ | vertical / horizontal line |
> | $\operatorname{Re}\{z\}<0$ | open left half-plane |
> | $\alpha<\arg z<\beta$ | open sector (wedge) |
> | $z=z_0+Re^{jt}$, $0\le t<2\pi$ | parametric form of the circle $\lvert z-z_0\rvert=R$ |

To identify an unfamiliar locus, the reliable method is to substitute $z=x+jy$, square any magnitudes (using $|w|^2=ww^*$ or simply $|a+jb|^2=a^2+b^2$), and recognise the resulting real equation in $x$ and $y$.

> [!example] Worked Example 2.27 — A half-plane
> Describe $\{z:|z-1|\le|z+1|\}$.
>
> **Geometric reading.** Points no farther from $1$ than from $-1$. The boundary $|z-1|=|z+1|$ is the set of points equidistant from $\pm1$: the perpendicular bisector, which is the imaginary axis. Testing $z=2$: $|1|\le|3|$ holds, so the region is the side containing $2$.
>
> **Algebraic confirmation.** Both sides are non-negative, so square: $(x-1)^2+y^2\le(x+1)^2+y^2\iff-2x\le2x\iff x\ge0$.
>
> The region is the closed right half-plane $\operatorname{Re}\{z\}\ge0$. ✓

> [!example] Worked Example 2.28 — A circle of Apollonius
> Identify the locus $|z-j|=2|z+1|$.
>
> **Step 1 — substitute and square.**
> $$
> x^2+(y-1)^2=4\left[(x+1)^2+y^2\right].
> $$
>
> **Step 2 — expand and collect.** $x^2+y^2-2y+1=4x^2+8x+4+4y^2$, so $3x^2+3y^2+8x+2y+3=0$, i.e.
> $$
> x^2+y^2+\tfrac83x+\tfrac23y+1=0.
> $$
>
> **Step 3 — complete the square.**
> $$
> \left(x+\tfrac43\right)^2+\left(y+\tfrac13\right)^2=\tfrac{16}{9}+\tfrac19-1=\tfrac89.
> $$
> A circle with centre $-\tfrac43-j\tfrac13$ and radius $\sqrt{8}/3\approx0.943$.
>
> **Check** a point: the rightmost point of the circle is $z=-\tfrac43+0.943-j\tfrac13=-0.391-j0.333$. Then $|z-j|=|{-0.391}-j1.333|=1.389$ and $|z+1|=|0.609-j0.333|=0.694$; ratio $2.00$. ✓
>
> The circle encloses $-1$ (the point the ratio favours) but not $j$. Circles of this kind are the constant-gain contours of a single-pole, single-zero system: $|z-j|/|z+1|$ is the magnitude of $(z-j)/(z+1)$.

```desmos-graph
left=-3.5; right=2; bottom=-2; top=1.8
---
(x+\frac{4}{3})^2+(y+\frac{1}{3})^2=\frac{8}{9}|blue
(0,1)|red|label:j
(-1,0)|green|label:-1
(-1.3333,-0.3333)|blue|label:centre
```
*The locus $\lvert z-j\rvert=2\lvert z+1\rvert$ is a circle around $-1$; every point on it is twice as far from $j$ as from $-1$.*

### Mappings: what $1/z$ does to lines

The map $w=1/z$ is a good first example of a complex function viewed as a transformation of the plane. Its key property — lines and circles map to lines and circles — is the basis of the Smith chart.

> [!example] Worked Example 2.29 — Image of a line under inversion
> Show that $w=1/z$ maps the vertical line $\operatorname{Re}\{z\}=1$ onto the circle $|w-\tfrac12|=\tfrac12$ (minus the origin).
>
> **Step 1 — parametrise the line.** $z=1+jy$, $y\in\mathbb{R}$.
>
> **Step 2 — compute $w$.**
> $$
> w=\frac{1}{1+jy}=\frac{1-jy}{1+y^2},\qquad u=\operatorname{Re}\{w\}=\frac{1}{1+y^2},\quad v=\operatorname{Im}\{w\}=\frac{-y}{1+y^2}.
> $$
>
> **Step 3 — eliminate $y$.** $u^2+v^2=\dfrac{1+y^2}{(1+y^2)^2}=\dfrac{1}{1+y^2}=u$, so $u^2-u+v^2=0$, i.e.
> $$
> \left(u-\tfrac12\right)^2+v^2=\tfrac14.
> $$
> As $y\to\pm\infty$, $w\to0$, which is approached but never reached.
>
> **Check:** $y=1$ gives $w=\frac{1-j}{2}$, and $|w-\tfrac12|=|{-j}/2|=\tfrac12$. ✓
>
> *Engineering reading:* the line $\operatorname{Re}\{Z\}=1\ \Omega$ is "all impedances with $1\ \Omega$ resistance". Its image under $Y=1/Z$ shows that the corresponding admittances fill a circle in the admittance plane — exactly the constant-resistance circles of impedance/admittance charts.

```desmos-graph
left=-1.5; right=2.5; bottom=-1.5; top=1.5
---
x=1|blue
(x-0.5)^2+y^2=0.25|red
(1,1)|blue|label:z=1+j
(0.5,-0.5)|red|label:1/z
```
*Inversion $w=1/z$ maps the vertical line $\operatorname{Re}\{z\}=1$ (blue) onto the circle of diameter 1 through the origin (red). The point $1+j$ maps to $\tfrac12-j\tfrac12$: inversion also reflects in the real axis.*

> [!tip] Engineering Note — the Smith chart and stability regions
> **Smith chart.** The reflection coefficient of a load $Z_L$ on a line of characteristic impedance $Z_0$ is $\Gamma=\dfrac{z-1}{z+1}$ with $z=Z_L/Z_0$. This map (a **bilinear** or Möbius transformation, studied in Ch. 18) sends the entire right half-plane $\operatorname{Re}\{z\}\ge0$ — every passive load — into the unit disc $|\Gamma|\le1$. A matched load $z=1$ lands at the centre $\Gamma=0$; any pure reactance $z=jx$ lands on the rim $|\Gamma|=1$ (total reflection). For example $z=1+j$ gives $\Gamma=\dfrac{j}{2+j}=0.2+j0.4$, $|\Gamma|=0.447$. The Smith chart is the unit disc with the images of the lines $\operatorname{Re}\{z\}=\text{const}$ and $\operatorname{Im}\{z\}=\text{const}$ drawn on it — all circles, by the line-to-circle property.
>
> **Stability.** A continuous-time linear time-invariant system is stable when every pole of its transfer function lies in the open left half-plane $\operatorname{Re}\{s\}<0$; a discrete-time system is stable when every pole lies inside the unit circle $|z|<1$ (Ch. 14, 17). Lines of constant damping ratio are rays $\arg s=\pm(\pi-\cos^{-1}\zeta)$; contours of constant settling time are vertical lines $\operatorname{Re}\{s\}=-\sigma$. Root locus, Nyquist and pole–zero plots are all geometry in this plane.

> [!warning] Common Pitfall — reading "$\arg$" loci
> The equation $\arg(z-z_0)=\alpha$ describes a **ray**, not a full line: the opposite ray has argument $\alpha\pm\pi$. Likewise $\arg\!\left(\frac{z-z_1}{z-z_2}\right)=\alpha$ is an **arc** of a circle through $z_1$ and $z_2$, not the whole circle (the other arc has argument $\alpha-\pi$). Squaring or taking tangents of an argument equation silently adds the missing pieces, so check candidate points against the original equation.

## 2.8 Applications to AC circuits and phasors

### Motivation

A linear circuit driven by a sinusoidal source is governed by a linear differential equation. Solving it directly — guessing $A\cos\omega t+B\sin\omega t$, substituting, and grinding through trigonometric identities — is laborious even for a single loop and impractical for a network. The **phasor method** replaces the differential equation by an algebraic equation in complex numbers. It rests on three facts already established: a sinusoid is the real part of a complex exponential (Section 2.4), complex arithmetic follows all the ordinary rules (Section 2.3), and differentiating $e^{j\omega t}$ merely multiplies it by $j\omega$.

The method applies to the **sinusoidal steady state**: the response that remains after transients have died out, when every voltage and current is a sinusoid at the source frequency $\omega$. For a stable linear time-invariant circuit, such a steady state always exists and differs from the input only in amplitude and phase. Transients are the subject of Ch. 11–12 and the Laplace transform (Ch. 14).

### Sinusoids as real parts

A sinusoid is specified by an amplitude $V_m>0$, an angular frequency $\omega$ and a phase $\phi$ ([Ch. 3](03-trigonometry-and-periodic-functions.md) treats these in detail). By Euler's formula

$$
v(t)=V_m\cos(\omega t+\phi)=\operatorname{Re}\left\{V_me^{j(\omega t+\phi)}\right\}=\operatorname{Re}\left\{\left(V_me^{j\phi}\right)e^{j\omega t}\right\}.
$$

The factor $e^{j\omega t}$ is common to every signal in the circuit. All the information that distinguishes one signal from another is in the complex constant $V_me^{j\phi}$.

> [!definition] Phasor
> The **phasor** of the sinusoid $v(t)=V_m\cos(\omega t+\phi)$ is the complex constant
> $$
> \mathbf{V}=V_me^{j\phi}=V_m\angle\phi.
> $$
> The sinusoid is recovered by
> $$
> v(t)=\operatorname{Re}\left\{\mathbf{V}e^{j\omega t}\right\}.
> $$
> The magnitude of the phasor is the amplitude; its angle is the phase. The phasor does not contain $\omega$ — the frequency is a shared, separately remembered parameter.

**Conventions.** This book uses **amplitude (peak) phasors** and the **cosine reference**. A sine is first converted to a cosine: $V_m\sin(\omega t+\phi)=V_m\cos(\omega t+\phi-90^\circ)$, so its phasor is $V_m\angle(\phi-90^\circ)=-jV_me^{j\phi}$. A negative amplitude is absorbed into the phase: $-V_m\cos(\omega t+\phi)=V_m\cos(\omega t+\phi\pm180^\circ)$. Power-systems texts commonly use **RMS phasors** with $|\mathbf{V}|=V_m/\sqrt2$; for example the UK mains is "$230$ V" RMS, a peak of $230\sqrt2=325.3$ V. Every formula in this section holds for either convention, except power formulas, which lose the factor $\tfrac12$ with RMS phasors (Section 2.9). State which convention you are using.

### The phasor transform is one-to-one and linear

> [!theorem] Uniqueness of phasors
> Let $\mathbf{A},\mathbf{B}$ be complex constants and $\omega\neq0$. If
> $$
> \operatorname{Re}\{\mathbf{A}e^{j\omega t}\}=\operatorname{Re}\{\mathbf{B}e^{j\omega t}\}\quad\text{for all }t,
> $$
> then $\mathbf{A}=\mathbf{B}$.

*Proof.* Put $\mathbf{C}=\mathbf{A}-\mathbf{B}=c_1+jc_2$; we have $\operatorname{Re}\{\mathbf{C}e^{j\omega t}\}=c_1\cos\omega t-c_2\sin\omega t=0$ for all $t$. At $t=0$ this gives $c_1=0$; at $\omega t=\pi/2$ it gives $-c_2=0$. So $\mathbf{C}=0$. $\square$

This theorem is what licenses the method: if two sinusoidal expressions at the same frequency are equal for all time, their phasors are equal, and conversely. Combined with the following rules, it converts any linear constant-coefficient relation between same-frequency sinusoids into the same relation between phasors.

> [!theorem] Phasor properties
> Let $x(t)\leftrightarrow\mathbf{X}$ and $y(t)\leftrightarrow\mathbf{Y}$ be sinusoids at a common frequency $\omega$. Then:
> - **Linearity:** $a\,x(t)+b\,y(t)\ \leftrightarrow\ a\mathbf{X}+b\mathbf{Y}$ for real constants $a,b$.
> - **Differentiation:** $\dfrac{dx}{dt}\ \leftrightarrow\ j\omega\mathbf{X}$, and $\dfrac{d^nx}{dt^n}\ \leftrightarrow\ (j\omega)^n\mathbf{X}$.
> - **Integration** (the zero-mean, steady-state antiderivative): $\displaystyle\int x\,dt\ \leftrightarrow\ \frac{\mathbf{X}}{j\omega}$.
> - **Time delay:** $x(t-\tau)\ \leftrightarrow\ e^{-j\omega\tau}\mathbf{X}$.

*Proof of the differentiation rule.* The real part commutes with differentiation with respect to the real variable $t$, so

$$
\frac{d}{dt}\operatorname{Re}\{\mathbf{X}e^{j\omega t}\}=\operatorname{Re}\left\{\mathbf{X}\frac{d}{dt}e^{j\omega t}\right\}=\operatorname{Re}\{(j\omega\mathbf{X})e^{j\omega t}\}.
$$

Linearity follows because $\operatorname{Re}$ is additive and commutes with multiplication by **real** constants. The delay rule: $\operatorname{Re}\{\mathbf{X}e^{j\omega(t-\tau)}\}=\operatorname{Re}\{(\mathbf{X}e^{-j\omega\tau})e^{j\omega t}\}$. $\square$

In words: **differentiation multiplies the amplitude by $\omega$ and advances the phase by $90^\circ$**, since $j\omega=\omega\angle90^\circ$. Check directly: $\frac{d}{dt}\cos\omega t=-\omega\sin\omega t=\omega\cos(\omega t+90^\circ)$. ✓ A delay of $\tau$ seconds is a phase lag of $\omega\tau$ radians.

### The method

1. Convert every source to a cosine and write its phasor.
2. Replace every $d/dt$ by $j\omega$ (equivalently, replace each circuit element by its impedance — Section 2.9).
3. Solve the resulting **algebraic** equations with complex arithmetic.
4. Convert each phasor answer back: $\mathbf{X}=X_m\angle\phi\ \Rightarrow\ x(t)=X_m\cos(\omega t+\phi)$.

> [!example] Worked Example 2.30 — Adding sinusoids
> Express $v(t)=10\cos(\omega t+30^\circ)+5\cos(\omega t-45^\circ)$ as a single sinusoid.
>
> **Step 1 — phasors in rectangular form** (addition needs rectangular):
> $$
> \begin{aligned}
> 10\angle30^\circ&=8.660+j5.000\\
> 5\angle{-45^\circ}&=3.536-j3.536\\
> \text{sum}&=12.196+j1.464
> \end{aligned}
> $$
>
> **Step 2 — back to polar.** $\sqrt{12.196^2+1.464^2}=12.28$; $\operatorname{atan2}(1.464,12.196)=6.85^\circ$.
>
> **Step 3.**
> $$
> v(t)=12.28\cos(\omega t+6.85^\circ).
> $$
>
> **Checks.** At $t=0$: $10\cos30^\circ+5\cos45^\circ=8.660+3.536=12.196$ and $12.28\cos6.85^\circ=12.196$. ✓ Law of cosines (Section 2.3): $\sqrt{100+25+2(10)(5)\cos75^\circ}=\sqrt{150.88}=12.28$. ✓
>
> Note $12.28<10+5$: same-frequency sinusoids add as **vectors**, not as amplitudes.

```desmos-graph
left=0; right=12.6; bottom=-15; top=15
---
y=10\cos(x+0.5236)|blue|dashed
y=5\cos(x-0.7854)|green|dashed
y=12.28\cos(x+0.1196)|red
```
*The sum (red) of $10\cos(x+30^\circ)$ (blue) and $5\cos(x-45^\circ)$ (green) is again a sinusoid of the same frequency, with amplitude $12.28$ and phase $6.85^\circ$ — exactly as the phasor sum predicts.*

> [!example] Worked Example 2.31 — Mixing sines and cosines
> Express $v(t)=4\sin\omega t-3\cos(\omega t+60^\circ)$ as a single cosine.
>
> **Step 1 — cosine reference.** $4\sin\omega t=4\cos(\omega t-90^\circ)\leftrightarrow4\angle{-90^\circ}=-j4$. The second term has phasor $3\angle60^\circ=1.500+j2.598$.
>
> **Step 2 — subtract.** $\mathbf{V}=-j4-(1.500+j2.598)=-1.500-j6.598$.
>
> **Step 3 — polar.** $|\mathbf{V}|=\sqrt{2.25+43.53}=6.766$. The point is in the **third quadrant**, so $\angle\mathbf{V}=\arctan\frac{6.598}{1.5}-180^\circ=77.19^\circ-180^\circ=-102.81^\circ$.
>
> $$
> v(t)=6.766\cos(\omega t-102.81^\circ).
> $$
>
> **Check** at $t=0$: $0-3\cos60^\circ=-1.5$ and $6.766\cos(-102.81^\circ)=-1.500$. ✓ At $\omega t=90^\circ$: $4-3\cos150^\circ=4+2.598=6.598$, and $6.766\cos(-12.81^\circ)=6.598$. ✓

> [!example] Worked Example 2.32 — A first-order differential equation
> A series RL circuit with $R=30\ \Omega$, $L=0.1$ H is driven by $v(t)=100\cos400t$ V. Its current obeys
> $$
> L\frac{di}{dt}+Ri=v(t).
> $$
> Find the steady-state current, then the complete response for $i(0)=0$.
>
> **Step 1 — phasor equation.** With $\mathbf{V}=100\angle0^\circ$ and $d/dt\to j\omega=j400$:
> $$
> (j\omega L+R)\mathbf{I}=\mathbf{V}\ \Rightarrow\ \mathbf{I}=\frac{100}{30+j40}=\frac{100\angle0^\circ}{50\angle53.13^\circ}=2\angle{-53.13^\circ}\ \text{A}.
> $$
>
> **Step 2 — back to time.** $i_{ss}(t)=2\cos(400t-53.13^\circ)$ A.
>
> **Step 3 — verify by substitution.** In rectangular form $\mathbf{I}=1.2-j1.6$, i.e. $i_{ss}=1.2\cos400t+1.6\sin400t$. Then
> $$
> L\frac{di_{ss}}{dt}=0.1(-480\sin400t+640\cos400t)=64\cos400t-48\sin400t,
> $$
> $$
> Ri_{ss}=36\cos400t+48\sin400t,
> $$
> and the sum is $100\cos400t$. ✓ The phasor method did in one division what substitution did in several lines.
>
> **Step 4 — complete response.** The homogeneous solution of $Li'+Ri=0$ is $Ke^{-(R/L)t}=Ke^{-300t}$ (Ch. 11). Then $i(0)=1.2+K=0$ gives $K=-1.2$:
> $$
> i(t)=2\cos(400t-53.13^\circ)-1.2\,e^{-300t}\ \text{A}.
> $$
> The transient decays with time constant $L/R=3.33$ ms; after about five time constants only the phasor solution remains.

> [!example] Worked Example 2.33 — A second-order differential equation
> Find the steady-state solution of
> $$
> \frac{d^2v}{dt^2}+3\frac{dv}{dt}+2v=10\cos2t.
> $$
>
> **Step 1.** With $\omega=2$, $d/dt\to j2$ and $d^2/dt^2\to(j2)^2=-4$:
> $$
> (-4+j6+2)\mathbf{V}=10\ \Rightarrow\ \mathbf{V}=\frac{10}{-2+j6}.
> $$
>
> **Step 2.** $\mathbf{V}=\dfrac{10(-2-j6)}{4+36}=-0.5-j1.5=1.581\angle{-108.43^\circ}$ (third quadrant).
>
> **Step 3.** $v_{ss}(t)=1.581\cos(2t-108.43^\circ)$.
>
> **Check** (rectangular): $v_{ss}=-0.5\cos2t+1.5\sin2t$, $v'=\sin2t+3\cos2t$, $v''=2\cos2t-6\sin2t$. Then $v''+3v'+2v=(2+9-1)\cos2t+(-6+3+3)\sin2t=10\cos2t$. ✓
>
> The complex number $(j\omega)^2+3(j\omega)+2$ is the characteristic polynomial $p(s)=s^2+3s+2$ evaluated at $s=j\omega$. In general the steady-state response of $p(D)y=x$ to $X_m\cos\omega t$ has phasor $X_m/p(j\omega)$ — the frequency response $H(j\omega)=1/p(j\omega)$ of Section 2.9.

### Phasor diagrams

Drawing phasors as arrows in the complex plane shows phase relationships at a glance. The whole diagram rotates anticlockwise at $\omega$ rad/s ($\mathbf{V}e^{j\omega t}$), and the instantaneous value of each signal is the projection of its rotating arrow onto the real axis. Relative angles never change, so the frozen diagram at $t=0$ carries all the information.

A phasor at a larger angle **leads**: its peaks occur earlier in time. If $\mathbf{V}=V_m\angle\phi_v$ and $\mathbf{I}=I_m\angle\phi_i$, the voltage leads the current by $\phi_v-\phi_i$ (a negative value means it lags). For a capacitor, $i=C\,dv/dt$ gives $\mathbf{I}=j\omega C\mathbf{V}$, so the current **leads** the voltage by $90^\circ$; for an inductor, $\mathbf{V}=j\omega L\mathbf{I}$, so the current **lags** by $90^\circ$. (The mnemonic "ELI the ICE man": voltage $E$ before current $I$ in an inductor $L$; current before voltage in a capacitor $C$.)

> [!example] Worked Example 2.34 — Three-phase line voltages
> A balanced three-phase supply has phase voltages (RMS phasors) $\mathbf{V}_a=V\angle0^\circ$, $\mathbf{V}_b=V\angle{-120^\circ}$, $\mathbf{V}_c=V\angle120^\circ$. Find the line voltage $\mathbf{V}_{ab}=\mathbf{V}_a-\mathbf{V}_b$.
>
> **Step 1.** $\mathbf{V}_{ab}=V\left(1-e^{-j2\pi/3}\right)=V\left(1+\tfrac12+j\tfrac{\sqrt3}{2}\right)=V\left(\tfrac32+j\tfrac{\sqrt3}{2}\right)$.
>
> **Step 2.** $\left|\tfrac32+j\tfrac{\sqrt3}{2}\right|=\sqrt{\tfrac94+\tfrac34}=\sqrt3$ and the angle is $\arctan\frac{1}{\sqrt3}=30^\circ$:
> $$
> \mathbf{V}_{ab}=\sqrt3\,V\angle30^\circ.
> $$
>
> With $V=230$ V the line voltage is $\sqrt3\times230=398.4$ V, quoted as the familiar "$400$ V" three-phase supply. Geometrically, the tips of $\mathbf{V}_a,\mathbf{V}_b,\mathbf{V}_c$ form an equilateral triangle and the line voltages are its sides (compare Example 2.26). Also $\mathbf{V}_a+\mathbf{V}_b+\mathbf{V}_c=V(1+a^2+a)=0$ (Section 2.6). ✓

> [!warning] Common Pitfall — different frequencies
> Phasors only combine signals of the **same** frequency. $\cos\omega_1t+\cos\omega_2t$ with $\omega_1\ne\omega_2$ cannot be written as one phasor, because the factor $e^{j\omega t}$ is no longer common. For a circuit driven at several frequencies (or by a periodic non-sinusoidal waveform), solve a separate phasor problem at each frequency and **add the time-domain results** (superposition) — the basis of Fourier-series circuit analysis (Ch. 15).

> [!warning] Common Pitfall — mixing phasor and time domains
> Expressions such as "$i(t)=2\angle{-53.13^\circ}$" or "$\mathbf{I}=2\cos(400t-53.13^\circ)$" are category errors: a phasor is a time-independent complex number, a signal is a real function of time. Keep the arrow $\leftrightarrow$ between them and convert explicitly. Likewise, never multiply two phasors and expect the phasor of the product of the signals — the product of two sinusoids at $\omega$ contains frequencies $0$ and $2\omega$ (Example 2.13b).

> [!tip] Engineering Note — what phasors buy you
> Every DC circuit-analysis technique — nodal and mesh analysis, superposition, source transformation, Thévenin and Norton equivalents — carries over unchanged to the sinusoidal steady state with resistances replaced by complex impedances and real sources by phasors. A network of fifty components at one frequency becomes a complex linear system $\mathbf{Y}\mathbf{v}=\mathbf{i}$, solved by the same Gaussian elimination as a real one (Ch. 8). This is precisely what circuit simulators do in "AC analysis" mode, sweeping $\omega$ and solving one complex linear system per frequency.

## 2.9 Complex impedance and transfer functions

### Motivation

Section 2.8 converted differential equations into algebra one equation at a time. Engineers go one step further and attach the conversion to the **components**: each resistor, inductor and capacitor is assigned a complex number, its impedance, and the circuit is then analysed exactly like a resistive network. The ratio of an output phasor to an input phasor, viewed as a function of frequency, is the **frequency response**; extended to the complex variable $s$ it becomes the **transfer function**, the central object of circuit theory, filtering and control. Both are nothing but complex-valued functions evaluated with the arithmetic of this chapter.

### Impedance and admittance

> [!definition] Impedance and admittance
> For a linear two-terminal element or network in sinusoidal steady state at angular frequency $\omega$, with terminal voltage phasor $\mathbf{V}$ and current phasor $\mathbf{I}$ (passive sign convention), the **impedance** is
> $$
> Z(j\omega)=\frac{\mathbf{V}}{\mathbf{I}}=R+jX\quad[\Omega],
> $$
> and the **admittance** is
> $$
> Y(j\omega)=\frac{1}{Z}=\frac{\mathbf{I}}{\mathbf{V}}=G+jB\quad[\mathrm{S}].
> $$
> $R=\operatorname{Re}\{Z\}$ is the **resistance**, $X=\operatorname{Im}\{Z\}$ the **reactance**, $G$ the **conductance** and $B$ the **susceptance**. $X>0$ is called inductive, $X<0$ capacitive.

In polar form $Z=|Z|\angle\theta_Z$, Ohm's law $\mathbf{V}=Z\mathbf{I}$ says two things at once: the voltage amplitude is $|Z|$ times the current amplitude, and the voltage leads the current by $\theta_Z$. Impedance is a complex number, **not** a phasor: it does not represent a sinusoid, it relates two of them.

Note that $G\ne1/R$ in general. From $Y=1/(R+jX)=(R-jX)/(R^2+X^2)$,

$$
G=\frac{R}{R^2+X^2},\qquad B=\frac{-X}{R^2+X^2}.
$$

For example $Z=30+j40\ \Omega$ gives $Y=(30-j40)/2500=0.012-j0.016$ S: the conductance is $0.012$ S, not $1/30=0.0333$ S.

### Element impedances

Apply the phasor differentiation rule ($d/dt\leftrightarrow j\omega$) to each element law:

| Element | Time-domain law | Phasor law | Impedance | $\angle Z$ |
|---|---|---|---|---|
| Resistor | $v=Ri$ | $\mathbf{V}=R\mathbf{I}$ | $Z_R=R$ | $0^\circ$ |
| Inductor | $v=L\,\dfrac{di}{dt}$ | $\mathbf{V}=j\omega L\,\mathbf{I}$ | $Z_L=j\omega L$ | $+90^\circ$ |
| Capacitor | $i=C\,\dfrac{dv}{dt}$ | $\mathbf{I}=j\omega C\,\mathbf{V}$ | $Z_C=\dfrac{1}{j\omega C}=-\dfrac{j}{\omega C}$ | $-90^\circ$ |

The inductive reactance $X_L=\omega L$ rises with frequency (an inductor is a short circuit at DC and an open circuit at very high frequency); the capacitive reactance $X_C=-1/(\omega C)$ has magnitude falling with frequency (open at DC, short at high frequency). Dimensional check: $\omega L$ has units $\mathrm{s^{-1}\cdot H}=\mathrm{s^{-1}\cdot\Omega\,s}=\Omega$, and $1/(\omega C)$ has units $1/(\mathrm{s^{-1}\cdot F})=\mathrm{s}/(\mathrm{A\,s/V})=\Omega$. ✓

### Combining impedances

Because KVL and KCL are linear relations with real ($\pm1$) coefficients, they hold for phasors (Section 2.8). Hence impedances combine exactly as resistances do:

> [!theorem] Series, parallel and divider rules
> $$
> Z_{\text{series}}=Z_1+Z_2+\cdots,\qquad \frac{1}{Z_{\text{par}}}=\frac{1}{Z_1}+\frac{1}{Z_2}+\cdots\quad\left(Z_{\text{par}}=\frac{Z_1Z_2}{Z_1+Z_2}\text{ for two}\right),
> $$
> $$
> \mathbf{V}_2=\frac{Z_2}{Z_1+Z_2}\,\mathbf{V}\quad\text{(voltage divider)},\qquad \mathbf{I}_2=\frac{Z_1}{Z_1+Z_2}\,\mathbf{I}\quad\text{(current divider)}.
> $$
> Admittances add in parallel: $Y_{\text{par}}=Y_1+Y_2+\cdots$.

> [!example] Worked Example 2.35 — Series RL: current and power
> $R=30\ \Omega$, $L=0.1$ H, source $v(t)=100\cos400t$ V (the circuit of Example 2.32). Find $i(t)$ and the power drawn.
>
> **Impedance.** $X_L=\omega L=400(0.1)=40\ \Omega$, so $Z=30+j40=50\angle53.13^\circ\ \Omega$.
>
> **Current.** $\mathbf{I}=\mathbf{V}/Z=\dfrac{100\angle0^\circ}{50\angle53.13^\circ}=2\angle{-53.13^\circ}$ A, so $i(t)=2\cos(400t-53.13^\circ)$ A — the same as Example 2.32, now with no differential equation in sight. The current lags the voltage by $53.13^\circ$ (inductive).
>
> **Power** (formulas derived below). $S=\tfrac12\mathbf{V}\mathbf{I}^*=\tfrac12(100)(2)\angle53.13^\circ=100\angle53.13^\circ=60+j80$. Average power $P=60$ W, reactive power $Q=80$ var.
>
> **Check:** $P=\tfrac12|\mathbf{I}|^2R=\tfrac12(4)(30)=60$ W ✓ and $Q=\tfrac12|\mathbf{I}|^2X_L=\tfrac12(4)(40)=80$ var. ✓ Only the resistor dissipates power.

> [!example] Worked Example 2.36 — A series–parallel network
> At $\omega=1000$ rad/s, a $10\ \Omega$ resistor is in series with the parallel combination of a $20\ \Omega$ resistor and a $50\ \mu$F capacitor. A source $\mathbf{V}=100\angle0^\circ$ V drives the network. Find the input impedance, the source current and the voltage across the parallel part.
>
> **Step 1 — capacitor.** $Z_C=\dfrac{1}{j(1000)(50\times10^{-6})}=\dfrac{1}{j0.05}=-j20\ \Omega$.
>
> **Step 2 — parallel part.**
> $$
> Z_p=\frac{20(-j20)}{20-j20}=\frac{-j400}{20-j20}=\frac{-j400(20+j20)}{800}=\frac{8000-j8000}{800}=10-j10\ \Omega.
> $$
>
> **Step 3 — total.** $Z=10+(10-j10)=20-j10=22.36\angle{-26.57^\circ}\ \Omega$ (capacitive).
>
> **Step 4 — current.** $\mathbf{I}=\dfrac{100}{22.36\angle{-26.57^\circ}}=4.472\angle26.57^\circ$ A: the current **leads** the voltage, as expected of a capacitive load.
>
> **Step 5 — divider.** $\mathbf{V}_p=Z_p\mathbf{I}=(14.14\angle{-45^\circ})(4.472\angle26.57^\circ)=63.25\angle{-18.43^\circ}=60-j20$ V.
>
> **Check** by KVL: the series resistor drops $10\mathbf{I}=44.72\angle26.57^\circ=40+j20$ V, and $(40+j20)+(60-j20)=100$. ✓ Note that $|\mathbf{V}_p|+|10\mathbf{I}|=63.25+44.72=107.97>100$: voltage **magnitudes** around a loop need not sum to the source magnitude — only the phasors do.

### Resonance

> [!example] Worked Example 2.37 — Series RLC resonance, bandwidth and $Q$
> A series RLC circuit has $R=10\ \Omega$, $L=10$ mH, $C=1\ \mu$F. Its impedance is
> $$
> Z(j\omega)=R+j\left(\omega L-\frac{1}{\omega C}\right).
> $$
>
> **Resonance.** The reactance vanishes when $\omega L=1/(\omega C)$, i.e. at
> $$
> \omega_0=\frac{1}{\sqrt{LC}}=\frac{1}{\sqrt{10^{-2}\cdot10^{-6}}}=10^4\ \text{rad/s}\quad(f_0=1.592\ \text{kHz}).
> $$
> There $Z=R=10\ \Omega$: purely resistive and minimal in magnitude, so the current for a given voltage is maximal and in phase with it.
>
> **Off resonance.** At $\omega=2\times10^4$: $\omega L=200$, $1/(\omega C)=50$, $Z=10+j150=150.3\angle86.19^\circ\ \Omega$ — strongly inductive, and the current is $15$ times smaller than at resonance.
>
> **Half-power frequencies.** The power $\tfrac12|\mathbf{I}|^2R$ falls to half its peak when $|Z|=\sqrt2R$, i.e. when the reactance equals $\pm R$: $\omega L-\dfrac{1}{\omega C}=\pm R$. Multiplying by $\omega/L$ gives $\omega^2\mp\frac RL\omega-\omega_0^2=0$, whose positive roots are
> $$
> \omega_{1,2}=\mp\frac{R}{2L}+\sqrt{\left(\frac{R}{2L}\right)^2+\omega_0^2}=\mp500+\sqrt{500^2+10^8}=9512.5,\ 10512.5\ \text{rad/s}.
> $$
> **Bandwidth** $\Delta\omega=\omega_2-\omega_1=R/L=1000$ rad/s, and **quality factor** $Q=\omega_0/\Delta\omega=\omega_0L/R=10$.
>
> **Check** at $\omega_1=9512.5$: $\omega_1L=95.125$, $1/(\omega_1C)=105.125$, reactance $-10.00\ \Omega=-R$. ✓
>
> **Voltage magnification.** At resonance, with source amplitude $V$, $|\mathbf{I}|=V/R$ and the capacitor voltage is $|\mathbf{V}_C|=\dfrac{V}{R\,\omega_0C}=QV=10V$: the capacitor sees ten times the source voltage — a design hazard and the basis of tuned circuits.

Normalising $x=\omega/\omega_0$, the current relative to its peak is

$$
\frac{|\mathbf{I}(j\omega)|}{|\mathbf{I}|_{\max}}=\frac{R}{|Z|}=\frac{1}{\sqrt{1+Q^2\left(x-\dfrac1x\right)^2}},
$$

which shows that the sharpness of the peak depends only on $Q$.

```desmos-graph
left=0; right=3; bottom=0; top=1.1
---
y=\frac{1}{\sqrt{1+100(x-\frac{1}{x})^2}}|blue|0<x<3
y=\frac{1}{\sqrt{1+4(x-\frac{1}{x})^2}}|red|0<x<3
y=0.7071|grey|dashed
```
*Normalised series-RLC current $R/\lvert Z\rvert$ against $x=\omega/\omega_0$ for $Q=10$ (blue, narrow) and $Q=2$ (red, broad). The dashed line marks the half-power level $1/\sqrt2$; its intersections give the bandwidth $\omega_0/Q$.*

### Complex power

> [!theorem] Average and complex power
> With amplitude phasors $\mathbf{V}=V_m\angle\phi_v$ and $\mathbf{I}=I_m\angle\phi_i$, the instantaneous power $p(t)=v(t)i(t)$ has average
> $$
> P=\frac{1}{T}\int_0^Tp(t)\,dt=\tfrac12V_mI_m\cos(\phi_v-\phi_i)=\operatorname{Re}\{S\},\qquad S=\tfrac12\mathbf{V}\mathbf{I}^*=P+jQ.
> $$
> $S$ is the **complex power** [VA], $P$ the **average (real) power** [W], $Q=\tfrac12V_mI_m\sin(\phi_v-\phi_i)$ the **reactive power** [var], $|S|$ the **apparent power** and $\cos(\phi_v-\phi_i)$ the **power factor**. For a load $Z$: $S=\tfrac12|\mathbf{I}|^2Z$, so $P=\tfrac12|\mathbf{I}|^2R$ and $Q=\tfrac12|\mathbf{I}|^2X$. With RMS phasors, drop the factor $\tfrac12$: $S=\mathbf{V}_{\text{rms}}\mathbf{I}_{\text{rms}}^*$.

*Derivation.* Using the product-to-sum identity (Example 2.13b),

$$
p(t)=V_mI_m\cos(\omega t+\phi_v)\cos(\omega t+\phi_i)=\tfrac12V_mI_m\left[\cos(\phi_v-\phi_i)+\cos(2\omega t+\phi_v+\phi_i)\right].
$$

The second term oscillates at $2\omega$ and averages to zero over a period, leaving $P=\tfrac12V_mI_m\cos(\phi_v-\phi_i)$. Since $\mathbf{V}\mathbf{I}^*=V_mI_m\angle(\phi_v-\phi_i)$, this is $\operatorname{Re}\{\tfrac12\mathbf{V}\mathbf{I}^*\}$. The conjugate is essential: it makes the angle of $S$ the *difference* of phases, the only thing that physically matters. For $S=\tfrac12\mathbf{V}\mathbf{I}^*$ with $\mathbf{V}=Z\mathbf{I}$: $S=\tfrac12Z\mathbf{I}\mathbf{I}^*=\tfrac12|\mathbf{I}|^2Z$. $\square$

> [!example] Worked Example 2.38 — Power factor and its correction
> A $50$ Hz, $230$ V (RMS) supply feeds a load $Z=30+j40\ \Omega$. Find $P$, $Q$ and the power factor, then the parallel capacitor that raises the power factor to 1.
>
> **Step 1 — current** (RMS phasors). $\mathbf{I}=\dfrac{230\angle0^\circ}{50\angle53.13^\circ}=4.6\angle{-53.13^\circ}$ A.
>
> **Step 2 — power.** $S=\mathbf{V}\mathbf{I}^*=230(4.6)\angle53.13^\circ=1058\angle53.13^\circ=634.8+j846.4$ VA. So $P=634.8$ W, $Q=846.4$ var (inductive), power factor $\cos53.13^\circ=0.6$ lagging.
> *Check:* $|\mathbf{I}|^2R=21.16(30)=634.8$ W ✓, $|\mathbf{I}|^2X=21.16(40)=846.4$ var. ✓
>
> **Step 3 — correction.** A parallel capacitor draws $S_C=\mathbf{V}\mathbf{I}_C^*$ with $\mathbf{I}_C=j\omega C\mathbf{V}$, so $S_C=\mathbf{V}(-j\omega C\mathbf{V}^*)=-j\omega C|\mathbf{V}|^2$: purely negative reactive power. To cancel $Q$:
> $$
> \omega C|\mathbf{V}|^2=846.4\ \Rightarrow\ C=\frac{846.4}{(2\pi\cdot50)(230^2)}=50.9\ \mu\text{F}.
> $$
> The supply then delivers the same $634.8$ W with current $634.8/230=2.76$ A instead of $4.6$ A, reducing the resistive losses in the supply cables by a factor $0.6^2=0.36$.

> [!theorem] Maximum power transfer
> A source with phasor $\mathbf{V}_s$ and internal impedance $Z_s=R_s+jX_s$ delivers the greatest average power to a load $Z_L$ when
> $$
> Z_L=Z_s^*=R_s-jX_s,\qquad P_{\max}=\frac{|\mathbf{V}_s|^2}{8R_s}\ \text{(amplitude phasor)}.
> $$

*Proof sketch.* $P=\tfrac12|\mathbf{I}|^2R_L=\dfrac{\tfrac12|\mathbf{V}_s|^2R_L}{(R_s+R_L)^2+(X_s+X_L)^2}$. For any $R_L$ the denominator is smallest when $X_L=-X_s$; then $P=\tfrac12|\mathbf{V}_s|^2R_L/(R_s+R_L)^2$, maximised over $R_L$ at $R_L=R_s$ (set the derivative to zero, Ch. 4). $\square$ For example $Z_s=50+j30\ \Omega$ and $|\mathbf{V}_s|=10$ V give $Z_L=50-j30\ \Omega$, current $10/100=0.1$ A, and $P_{\max}=\tfrac12(0.1)^2(50)=0.25$ W $=100/(8\cdot50)$. ✓ This **conjugate match** is the target of every RF matching network.

### Frequency response and transfer functions

> [!definition] Frequency response
> For a stable linear time-invariant circuit with sinusoidal input $x(t)=X_m\cos(\omega t+\phi_x)\leftrightarrow\mathbf{X}$ and steady-state output $\mathbf{Y}$, the **frequency response** is
> $$
> H(j\omega)=\frac{\mathbf{Y}}{\mathbf{X}}.
> $$
> The steady-state output is then
> $$
> y(t)=|H(j\omega)|\,X_m\cos\big(\omega t+\phi_x+\angle H(j\omega)\big).
> $$
> $|H(j\omega)|$ is the **gain** (often quoted as $20\log_{10}|H|$ dB) and $\angle H(j\omega)$ the **phase shift**.

So an LTI system acts on a sinusoid by **one complex multiplication**: scale by $|H|$, rotate by $\angle H$. $H$ is computed from the circuit by impedance methods, or from the governing differential equation by $d/dt\to j\omega$ (Example 2.33).

> [!definition] Transfer function, poles and zeros
> When $H(j\omega)$ is a ratio of polynomials in $j\omega$, replacing $j\omega$ by the complex variable $s=\sigma+j\omega$ gives the **transfer function**
> $$
> H(s)=\frac{N(s)}{D(s)}=K\,\frac{(s-z_1)(s-z_2)\cdots(s-z_m)}{(s-p_1)(s-p_2)\cdots(s-p_n)}.
> $$
> The roots $z_i$ of $N$ are the **zeros** and the roots $p_i$ of $D$ the **poles**. For real circuits $N$ and $D$ have real coefficients, so poles and zeros are real or occur in conjugate pairs (Section 2.2). In terms of impedances, each element is assigned $Z_R=R$, $Z_L=sL$, $Z_C=1/(sC)$.

The interpretation of $H(s)$ for complex $s$ — the response to $e^{st}$ — is developed with the Laplace transform in [Ch. 14](../part-5-transforms/14-laplace-transforms.md). Here what matters is that the factored form turns the evaluation of gain and phase into geometry.

> [!theorem] Gain and phase from poles and zeros
> At $s=j\omega$,
> $$
> |H(j\omega)|=|K|\,\frac{\prod_i|j\omega-z_i|}{\prod_k|j\omega-p_k|},\qquad \angle H(j\omega)=\angle K+\sum_i\angle(j\omega-z_i)-\sum_k\angle(j\omega-p_k).
> $$
> Each factor $j\omega-z_i$ is the vector from the zero $z_i$ to the point $j\omega$ on the imaginary axis.

*Proof.* Magnitudes of products and quotients multiply and divide, angles add and subtract (Section 2.3). $\square$ As $\omega$ moves up the imaginary axis, the gain rises near zeros and peaks near poles that lie close to the axis — the geometric explanation of resonance.

> [!example] Worked Example 2.39 — RC low-pass filter
> A series resistor $R$ feeds an output taken across a capacitor $C$. Find $H(j\omega)=\mathbf{V}_{\text{out}}/\mathbf{V}_{\text{in}}$, its gain and phase, and its pole.
>
> **Step 1 — divider.**
> $$
> H(j\omega)=\frac{Z_C}{R+Z_C}=\frac{1/(j\omega C)}{R+1/(j\omega C)}=\frac{1}{1+j\omega RC},\qquad H(s)=\frac{1}{1+sRC}.
> $$
>
> **Step 2 — polar form.**
> $$
> |H(j\omega)|=\frac{1}{\sqrt{1+(\omega RC)^2}},\qquad \angle H(j\omega)=-\arctan(\omega RC).
> $$
> (The denominator $1+j\omega RC$ is in the first quadrant, so $\arctan$ is safe here.)
>
> **Step 3 — cut-off.** At $\omega_c=1/(RC)$: $|H|=1/\sqrt2=0.7071$, i.e. $20\log_{10}0.7071=-3.01$ dB, and $\angle H=-45^\circ$. With $R=1\ \text{k}\Omega$, $C=1\ \mu$F: $\omega_c=1000$ rad/s, $f_c=159.2$ Hz.
>
> **Step 4 — a specific input.** For $v_{\text{in}}=10\cos2000t$ V, $\omega RC=2$: $H=1/(1+j2)=0.4472\angle{-63.43^\circ}$ ($-6.99$ dB), so $v_{\text{out}}=4.472\cos(2000t-63.43^\circ)$ V.
>
> **Step 5 — pole.** $D(s)=1+sRC=0$ at $s=-1/(RC)=-1000\ \mathrm{s^{-1}}$, on the negative real axis (left half-plane: stable). Geometrically, $|H(j\omega)|=\omega_c/|j\omega+\omega_c|$ is $\omega_c$ divided by the distance from the pole to $j\omega$; that distance is $\sqrt2\,\omega_c$ when $\omega=\omega_c$, which is where the $-3$ dB point comes from.

```desmos-graph
left=0; right=5; bottom=-1.7; top=1.2
---
y=\frac{1}{\sqrt{1+x^2}}|blue
y=-\arctan(x)|red
y=0.7071|grey|dashed
(1,0.7071)|blue|label:-3 dB
(1,-0.7854)|red|label:-45°
```
*RC low-pass response against normalised frequency $x=\omega RC$: magnitude $\lvert H\rvert$ (blue) and phase $\angle H$ in radians (red), tending to $-\pi/2$. Both shapes follow directly from the polar form of $1/(1+jx)$.*

> [!example] Worked Example 2.40 — Gain and phase from a pole–zero plot
> Evaluate the gain and phase of
> $$
> H(s)=\frac{10(s+1)}{(s+2)(s+10)}
> $$
> at $\omega=2$ rad/s.
>
> **Step 1 — poles and zeros.** Zero at $-1$; poles at $-2$ and $-10$; $K=10$.
>
> **Step 2 — vectors to $j2$.**
> $$
> \begin{aligned}
> j2-(-1)&=1+j2=2.236\angle63.43^\circ,\\
> j2-(-2)&=2+j2=2.828\angle45.00^\circ,\\
> j2-(-10)&=10+j2=10.198\angle11.31^\circ.
> \end{aligned}
> $$
>
> **Step 3 — combine.**
> $$
> |H(j2)|=\frac{10\times2.236}{2.828\times10.198}=0.7752\ (-2.21\ \text{dB}),\qquad \angle H(j2)=63.43^\circ-45.00^\circ-11.31^\circ=7.12^\circ.
> $$
>
> **Check** by direct rectangular evaluation: $(2+j2)(10+j2)=20+j4+j20-4=16+j24$, and $\dfrac{10(1+j2)}{16+j24}=\dfrac{10(1+j2)(16-j24)}{832}=\dfrac{10(64+j8)}{832}=0.7692+j0.0962$, whose magnitude is $0.7752$ and angle $7.13^\circ$. ✓ (The last-digit difference is rounding in the angles.)

> [!warning] Common Pitfall — phase of a factor in the left half of its plane
> When a factor such as $j\omega-p$ or $1-\omega^2LC+j\omega RC$ has a **negative real part**, its angle exceeds $90^\circ$, and $\arctan(\text{Im}/\text{Re})$ is off by $180^\circ$ (Section 2.1). Second-order denominators $1-(\omega/\omega_0)^2+j\cdot(\dots)$ change quadrant at $\omega=\omega_0$, which is why a naive formula makes the phase of a resonant filter jump from $-90^\circ$ to $+90^\circ$ instead of continuing smoothly to $-180^\circ$. Use `atan2` and then unwrap.

> [!tip] Engineering Note — from here to Bode, filters and control
> Bode plots are $20\log_{10}|H(j\omega)|$ and $\angle H(j\omega)$ against $\log\omega$. Because the log of a product is the sum of logs (Section 2.5), each pole and zero contributes an independent, easily sketched term, and the total is their sum. Filter design chooses pole and zero locations (Butterworth on a circle, Example 2.25; Chebyshev on an ellipse); feedback-loop stability margins are read from $|H|$ and $\angle H$ near the crossover frequency. The mathematics needed is exactly this chapter plus logarithms; the design methods belong to circuits and control texts.

## Summary

| Topic | Key formula |
|---|---|
| Rectangular form | $z=a+jb$, $a=\operatorname{Re}\{z\}$, $b=\operatorname{Im}\{z\}$, $j^2=-1$ |
| Rectangular ↔ polar | $\lvert z\rvert=\sqrt{a^2+b^2}$, $\angle z=\operatorname{atan2}(b,a)$; $a=r\cos\theta$, $b=r\sin\theta$ |
| Equality | $a+jb=c+jd\iff a=c,\ b=d$ (one complex equation = two real equations) |
| Conjugate | $z^*=a-jb=re^{-j\theta}$; $zz^*=\lvert z\rvert^2$; $\operatorname{Re}\{z\}=\frac{z+z^*}{2}$, $\operatorname{Im}\{z\}=\frac{z-z^*}{2j}$ |
| Conjugate roots | real coefficients $\Rightarrow$ non-real roots in pairs $\alpha\pm j\beta$, factor $x^2-2\alpha x+\alpha^2+\beta^2$ |
| Product | $(a+jb)(c+jd)=(ac-bd)+j(ad+bc)$; $r_1\angle\theta_1\cdot r_2\angle\theta_2=r_1r_2\angle(\theta_1+\theta_2)$ |
| Division | $\dfrac{z_1}{z_2}=\dfrac{z_1z_2^*}{\lvert z_2\rvert^2}=\dfrac{r_1}{r_2}\angle(\theta_1-\theta_2)$ |
| Powers of $j$ | $j^n=j^{\,n\bmod4}$; $1/j=-j$ |
| Triangle inequality | $\lvert z_1+z_2\rvert\le\lvert z_1\rvert+\lvert z_2\rvert$; $\lvert z_1+z_2\rvert^2=\lvert z_1\rvert^2+\lvert z_2\rvert^2+2\operatorname{Re}\{z_1z_2^*\}$ |
| Euler | $e^{j\theta}=\cos\theta+j\sin\theta$; $e^{j\pi}=-1$, $e^{j\pi/2}=j$ |
| Inverse Euler | $\cos\theta=\frac{e^{j\theta}+e^{-j\theta}}{2}$, $\sin\theta=\frac{e^{j\theta}-e^{-j\theta}}{2j}$ |
| Exponential form | $z=re^{j\theta}$; $z_1z_2=r_1r_2e^{j(\theta_1+\theta_2)}$ |
| Complex exponential | $e^{\sigma+j\omega}=e^\sigma(\cos\omega+j\sin\omega)$, $\lvert e^z\rvert=e^{\operatorname{Re}\{z\}}$, period $j2\pi$ |
| Logarithm | $\operatorname{Log}z=\ln\lvert z\rvert+j\operatorname{Arg}z$; $\log z=\operatorname{Log}z+j2\pi k$ |
| Complex power of $z$ | $z^c=e^{c\log z}$ (multivalued; principal value uses $\operatorname{Log}$) |
| De Moivre | $(re^{j\theta})^n=r^ne^{jn\theta}$, $n\in\mathbb{Z}$ |
| $n$-th roots | $w_k=r^{1/n}e^{j(\theta+2\pi k)/n}$, $k=0,\dots,n-1$ |
| Square root | $\sqrt{a+jb}=\pm\left(\sqrt{\frac{r+a}{2}}+j\operatorname{sgn}(b)\sqrt{\frac{r-a}{2}}\right)$ |
| Roots of unity | $W=e^{j2\pi/n}$; $\sum_{k=0}^{n-1}W^{mk}=n$ if $n\mid m$, else $0$ |
| Rotation | multiply by $e^{j\alpha}$; by $j$ for $+90^\circ$; about $w_0$: $w_0+e^{j\alpha}(z-w_0)$ |
| Circle | $\lvert z-z_0\rvert=R$; Apollonius $\lvert z-z_1\rvert=k\lvert z-z_2\rvert$ |
| Phasor | $V_m\cos(\omega t+\phi)\leftrightarrow V_m\angle\phi$; $v(t)=\operatorname{Re}\{\mathbf{V}e^{j\omega t}\}$ |
| Phasor calculus | $\frac{d}{dt}\leftrightarrow j\omega$, $\int dt\leftrightarrow\frac{1}{j\omega}$, delay $\tau\leftrightarrow e^{-j\omega\tau}$ |
| Impedances | $Z_R=R$, $Z_L=j\omega L$, $Z_C=\dfrac{1}{j\omega C}$; series add, admittances add in parallel |
| Series RLC | $\omega_0=1/\sqrt{LC}$, $Q=\omega_0L/R$, bandwidth $R/L$ |
| Complex power | $S=\tfrac12\mathbf{V}\mathbf{I}^*=P+jQ=\tfrac12\lvert\mathbf{I}\rvert^2Z$ (amplitude); $S=\mathbf{V}\mathbf{I}^*$ (RMS) |
| Max power transfer | $Z_L=Z_s^*$, $P_{\max}=\lvert\mathbf{V}_s\rvert^2/(8R_s)$ (amplitude) |
| Frequency response | $y(t)=\lvert H\rvert X_m\cos(\omega t+\phi_x+\angle H)$ |
| Pole–zero evaluation | $\lvert H\rvert=\lvert K\rvert\frac{\prod\lvert j\omega-z_i\rvert}{\prod\lvert j\omega-p_k\rvert}$, $\angle H=\angle K+\sum\angle(j\omega-z_i)-\sum\angle(j\omega-p_k)$ |
| RC low-pass | $H=\dfrac{1}{1+j\omega RC}$, $\lvert H\rvert=\dfrac{1}{\sqrt{1+(\omega RC)^2}}$, $\angle H=-\arctan\omega RC$ |

## Further reading

- E. Kreyszig, *Advanced Engineering Mathematics*, 10th ed., Wiley, 2011 — Chapter 13 (complex numbers and functions, the complex exponential and logarithm).
- K. A. Stroud and D. J. Booth, *Engineering Mathematics*, 8th ed., Palgrave Macmillan, 2020 — Programmes on complex numbers, with extensive step-by-step drill.
- G. James et al., *Modern Engineering Mathematics*, 6th ed., Pearson, 2020 — Chapter 3 (complex numbers), including phasors and loci.
- A. V. Oppenheim and A. S. Willsky, *Signals and Systems*, 2nd ed., Prentice Hall, 1997 — Chapter 1 on complex exponential signals and Chapter 3 on eigenfunctions of LTI systems.
- C. K. Alexander and M. N. O. Sadiku, *Fundamentals of Electric Circuits*, 7th ed., McGraw-Hill, 2021 — Chapters 9–11 for the phasor method, impedance and AC power in full circuit-analysis detail.
