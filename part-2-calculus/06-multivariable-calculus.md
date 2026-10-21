# Chapter 6 — Multivariable Calculus

[← Previous](../part-2-calculus/05-integral-calculus.md) · [Table of Contents](../00-front/table-of-contents.md) · [Next →](../part-2-calculus/07-vector-calculus.md)

> **Prerequisites:** [Ch. 4](../part-2-calculus/04-differential-calculus.md), [Ch. 5](../part-2-calculus/05-integral-calculus.md); basic vectors (dot product) and $2\times2$ determinants.
> **Learning outcomes:** After this chapter you can
> - describe functions of several variables through their domains, level curves and level surfaces, and decide whether limits exist;
> - compute partial derivatives, total differentials and gradients, and use them for linear approximation, sensitivity and tolerance analysis;
> - evaluate directional derivatives and build Jacobian matrices for linearisation, small-signal models and Newton's method;
> - set up and evaluate double and triple integrals, changing the order of integration when that helps;
> - change variables in integrals, in particular to polar, cylindrical and spherical coordinates;
> - locate and classify extrema of functions of several variables, with and without constraints, using the Hessian test and Lagrange multipliers.

Almost every quantity an engineer measures depends on more than one thing. The power in a resistor depends on voltage *and* resistance; the potential in a region of space depends on $x$, $y$ and $z$; the collector current of a transistor depends on $v_{BE}$ and $v_{CE}$; the cutoff frequency of a filter depends on every component value in it. Chapters 4 and 5 built derivatives and integrals for functions of one variable. This chapter extends both ideas to functions of several variables, and the extension is not merely notational: new phenomena appear (limits that depend on the direction of approach, derivatives in every direction at once, regions of integration with shape), and new tools (the gradient, the Jacobian, the Hessian, the Lagrange multiplier) become the everyday language of circuit sensitivity, field theory, nonlinear circuit solvers and optimisation.

The chapter has three movements. Sections 6.1–6.6 are **differential**: they describe how a function of several variables changes locally, culminating in the Jacobian matrix as *the* derivative of a vector function. Sections 6.7–6.9 are **integral**: they add up quantities over areas and volumes and show how to choose coordinates that fit the geometry. Sections 6.10–6.11 combine the two to find **optimal** points, with and without constraints. Chapter 7 then adds direction to the picture with vector fields, and almost every result there rests on something proved here.

## 6.1 Functions of several variables

### Motivation

A single-variable function $y=f(x)$ is a curve; you can sketch it and read off slopes and areas. An engineer's quantities are rarely that simple. The resonant frequency of a series RLC circuit, $\omega_0=1/\sqrt{LC}$, is a function of two variables. The electrostatic potential $V(x,y,z)$ in a device is a function of three. The bit-error rate of a link depends on signal power, noise power, bandwidth and modulation order. Before we can differentiate or integrate such functions we need a vocabulary for describing them and a way to *see* them.

### Definitions

> [!definition] Function of several variables
> A real function of $n$ variables assigns to each point $\mathbf{x}=(x_1,\dots,x_n)$ of a set $D\subseteq\mathbb{R}^n$ exactly one real number $f(\mathbf{x})$. The set $D$ is the **domain** and the set of values taken is the **range**. For $n=2$ we usually write $z=f(x,y)$; the **graph** $\{(x,y,f(x,y))\}$ is a surface in $\mathbb{R}^3$. A function of position $f(x,y,z)$ is also called a **scalar field**.

Unless stated otherwise, the domain is the largest set on which the formula makes sense: no division by zero, no square roots of negative numbers, no logarithms of non-positive numbers. Physical constraints may shrink it further (a resistance is positive, a probability lies in $[0,1]$).

A few terms about sets in $\mathbb{R}^n$ will be needed later. A **neighbourhood** (open ball) of $\mathbf{a}$ is the set of points within distance $\delta>0$ of $\mathbf{a}$: $\lVert\mathbf{x}-\mathbf{a}\rVert<\delta$. A point is an **interior point** of $D$ if some neighbourhood of it lies entirely in $D$, and a **boundary point** if every neighbourhood contains points both in and out of $D$. $D$ is **open** if all its points are interior, **closed** if it contains all its boundary points, and **bounded** if it fits inside some ball. "Closed and bounded" regions are the ones on which continuous functions are guaranteed to have a maximum and a minimum (§6.10).

> [!example] Worked Example — finding and sketching a domain
> Find the domain of $f(x,y)=\sqrt{9-x^2-y^2}+\ln(y-x)$.
>
> *Step 1 (square root):* need $9-x^2-y^2\ge0$, i.e. $x^2+y^2\le9$ — the closed disc of radius 3.
>
> *Step 2 (logarithm):* need $y-x>0$, i.e. the open half-plane above the line $y=x$.
>
> *Step 3 (intersect):* $D=\{(x,y): x^2+y^2\le9,\ y>x\}$, a half-disc. It is bounded but neither open (it contains part of the circle) nor closed (it excludes the diameter on $y=x$). Check a point: $(0,1)$ gives $\sqrt8+\ln1=2.828$, defined; $(1,0)$ gives $\ln(-1)$, undefined. ✓

```desmos-graph
left=-4; right=4; bottom=-4; top=4
---
x^2+y^2\le9|blue
y>x|red
y=x|red|dashed
```

*The domain of $\sqrt{9-x^2-y^2}+\ln(y-x)$ is where the blue disc and the red half-plane overlap; the dashed diameter itself is excluded.*

### Level curves and level surfaces

Because a surface is hard to draw, we usually draw **level curves** (contours) $f(x,y)=c$ for several constants $c$, the way a topographic map shows height. Where the contours crowd together the surface is steep; where they spread out it is flat. For $f(x,y,z)$ the analogous objects are **level surfaces** $f(x,y,z)=c$; in electrostatics these are the equipotential surfaces.

> [!definition] Level set
> The level set of $f$ at height $c$ is $\{\mathbf{x}\in D: f(\mathbf{x})=c\}$. For $n=2$ it is (generically) a curve, for $n=3$ a surface. Level sets for different values of $c$ never intersect, because $f$ has only one value at each point.

| Function | Level sets | Engineering reading |
|---|---|---|
| $f=x^2+y^2$ | circles centred on the origin | equal-distance contours |
| $f=ax+by$ | parallel straight lines | equipotentials of a uniform field |
| $f=xy$ | rectangular hyperbolas | constant-power curves $P=VI$ |
| $V=\dfrac{q}{4\pi\varepsilon_0\sqrt{x^2+y^2+z^2}}$ | spheres | equipotentials of a point charge |
| $P=V^2/R$ | curves $V^2=cR$ (parabolas in the $RV$-plane) | constant-dissipation curves |
| $\omega_0=1/\sqrt{LC}$ | hyperbolas $LC=\text{const}$ | all $(L,C)$ pairs tuned to one frequency |

The last row is worth a moment's thought: a whole curve of $(L,C)$ pairs gives the same resonant frequency, and the designer moves along that curve to trade other properties (characteristic impedance $\sqrt{L/C}$, component cost). Level sets turn such design freedom into geometry.

The equipotentials of a pair of equal and opposite point charges at $(\mp1,0)$ (in the $z=0$ plane, constants dropped) are the level curves of

$$
V(x,y)=\frac{1}{\sqrt{(x+1)^2+y^2}}-\frac{1}{\sqrt{(x-1)^2+y^2}}.
$$

```desmos-graph
left=-3; right=3; bottom=-2; top=2
---
\frac{1}{\sqrt{(x+1)^2+y^2}}-\frac{1}{\sqrt{(x-1)^2+y^2}}=0.25|blue
\frac{1}{\sqrt{(x+1)^2+y^2}}-\frac{1}{\sqrt{(x-1)^2+y^2}}=0.6|blue
\frac{1}{\sqrt{(x+1)^2+y^2}}-\frac{1}{\sqrt{(x-1)^2+y^2}}=1.5|blue
\frac{1}{\sqrt{(x+1)^2+y^2}}-\frac{1}{\sqrt{(x-1)^2+y^2}}=-0.25|red
\frac{1}{\sqrt{(x+1)^2+y^2}}-\frac{1}{\sqrt{(x-1)^2+y^2}}=-0.6|red
\frac{1}{\sqrt{(x+1)^2+y^2}}-\frac{1}{\sqrt{(x-1)^2+y^2}}=-1.5|red
x=0|grey|dashed
(-1,0)|label:+q
(1,0)|label:-q
```

*Level curves of the dipole potential: positive contours (blue) circle the positive charge, negative contours (red) the negative one, and the $V=0$ contour is the bisecting line $x=0$.*

Three-variable functions cannot be graphed at all (the graph would live in $\mathbb{R}^4$), so level surfaces are the main visual tool. The **quadric surfaces** are the level sets of quadratic functions and appear constantly as constraint surfaces and equipotentials:

| Equation | Surface |
|---|---|
| $\dfrac{x^2}{a^2}+\dfrac{y^2}{b^2}+\dfrac{z^2}{c^2}=1$ | ellipsoid (sphere if $a=b=c$) |
| $x^2+y^2=a^2$ | circular cylinder along $z$ |
| $z=x^2+y^2$ | paraboloid of revolution |
| $z^2=x^2+y^2$ | double cone |
| $z=x^2-y^2$ | hyperbolic paraboloid ("saddle") |
| $x^2+y^2-z^2=1$ | hyperboloid of one sheet |

### Limits and continuity

> [!definition] Limit and continuity
> $\displaystyle\lim_{\mathbf{x}\to\mathbf{a}}f(\mathbf{x})=L$ means: for every $\varepsilon>0$ there is a $\delta>0$ such that $\lvert f(\mathbf{x})-L\rvert<\varepsilon$ whenever $0<\lVert\mathbf{x}-\mathbf{a}\rVert<\delta$ and $\mathbf{x}\in D$. The function is **continuous** at $\mathbf{a}$ if $\lim_{\mathbf{x}\to\mathbf{a}}f(\mathbf{x})=f(\mathbf{a})$.

The definition looks exactly like the one-variable definition in Chapter 4, but $\lVert\mathbf{x}-\mathbf{a}\rVert<\delta$ is now a disc (or ball), not an interval. That means $f$ must approach $L$ **along every path** into $\mathbf{a}$ — straight lines from every direction, parabolas, spirals, anything. In one variable there are only two directions of approach (left and right); in two there are infinitely many. Two consequences follow:

- If two paths give different limiting values, the limit does **not** exist.
- To prove that a limit **does** exist, you must bound $\lvert f-L\rvert$ by something that tends to zero independently of the path. Polar coordinates centred on $\mathbf{a}$ are the standard tool: write $x-a=r\cos\theta$, $y-b=r\sin\theta$ and show $\lvert f-L\rvert\le g(r)$ with $g(r)\to0$.

> [!example] Worked Example 1 — a limit that does not exist
> Let $f(x,y)=\dfrac{xy}{x^2+y^2}$ for $(x,y)\ne(0,0)$. Along the line $y=mx$,
> $$
> f(x,mx)=\frac{mx^2}{x^2+m^2x^2}=\frac{m}{1+m^2},
> $$
> which depends on $m$: it is $0$ along the axes ($m=0$) and $\tfrac12$ along $y=x$ ($m=1$). Different paths give different values, so the limit at the origin does not exist, and no choice of $f(0,0)$ makes $f$ continuous there. In polar form $f=\tfrac12\sin2\theta$: the function depends only on the *direction* from the origin, not on the distance, so it cannot settle to a single value as $r\to0$.

> [!example] Worked Example 2 — all straight lines agree, yet no limit
> Let $f(x,y)=\dfrac{x^2y}{x^4+y^2}$. Along $y=mx$ ($m\ne0$),
> $$
> f(x,mx)=\frac{mx^3}{x^4+m^2x^2}=\frac{mx}{x^2+m^2}\to0\quad(x\to0),
> $$
> and along the axes $f=0$ identically. Every straight line gives $0$. But along the parabola $y=x^2$,
> $$
> f(x,x^2)=\frac{x^4}{x^4+x^4}=\frac12 .
> $$
> So the limit does not exist. More generally, on $y=kx^2$ the value is $k/(1+k^2)$: the level curves of $f$ are parabolas through the origin, and every one of them reaches the origin carrying its own value.

```desmos-graph
left=-2; right=2; bottom=-1.5; top=1.5
---
y=x^2|blue
y=0.5x^2|green
y=3x^2|green
y=-x^2|red
y=-0.5x^2|orange
(0,0)|label:origin
```

*Level curves of $x^2y/(x^4+y^2)$ are parabolas $y=kx^2$ with value $k/(1+k^2)$: blue ($k=1$) carries $0.5$, green ($k=0.5$ and $k=3$) carry $0.4$ and $0.3$, red ($k=-1$) and orange ($k=-0.5$) carry $-0.5$ and $-0.4$. Every one of them enters the origin, so no single limit exists — although every straight line through the origin gives 0.*

> [!example] Worked Example 3 — proving a limit exists
> Show that $\displaystyle\lim_{(x,y)\to(0,0)}\frac{x^3+2x^2y^2}{x^2+y^2}=0$.
>
> In polar coordinates $x=r\cos\theta$, $y=r\sin\theta$, the denominator is $r^2$ and
> $$
> f=\frac{r^3\cos^3\theta+2r^4\cos^2\theta\sin^2\theta}{r^2}=r\cos^3\theta+2r^2\cos^2\theta\sin^2\theta .
> $$
> Since $\lvert\cos\theta\rvert\le1$ and $\lvert\sin\theta\rvert\le1$, $\lvert f\rvert\le r+2r^2$, which tends to $0$ as $r\to0$ regardless of $\theta$. Hence the limit is $0$, and defining $f(0,0)=0$ makes $f$ continuous. Numerically, at $(0.01,0.01)$: $f=(10^{-6}+2\times10^{-8})/(2\times10^{-4})=0.0051$, small as expected.

> [!warning] Common Pitfall
> Checking straight-line paths can prove that a limit does **not** exist, but never that it does — Worked Example 2 passes every straight-line test and still fails. Likewise, a bound in polar form must be independent of $\theta$: $\lvert f\rvert\le r/\lvert\cos\theta\rvert$ is *not* good enough, because $\cos\theta$ can be arbitrarily small.

> [!theorem] Building continuous functions
> Polynomials in $x_1,\dots,x_n$ are continuous everywhere. Sums, products and compositions of continuous functions are continuous; quotients are continuous wherever the denominator is non-zero. Elementary functions ($e^u$, $\sin u$, $\ln u$, $\sqrt u$, …) of continuous expressions are continuous wherever they are defined.

So $e^{-x^2-y^2}\cos(xy)$ is continuous on all of $\mathbb{R}^2$, and $\ln(1-x^2-y^2)$ on the open unit disc. Trouble only arises at points excluded from the formula's natural domain, as at the origin above.

> [!tip] Engineering Note
> Real devices are continuous functions of their inputs, but their *models* are not always: piecewise device models (a diode switching between "off" and "on" regions, a MOSFET between triode and saturation) can be discontinuous or have discontinuous derivatives at region boundaries. Circuit simulators fail to converge exactly at such seams, which is why modern compact models are built from smooth functions that blend the regions together.

## 6.2 Partial differentiation

### Motivation

The first question an engineer asks of a multi-input function is: *if I nudge one input and hold the rest, how much does the output move?* That is a partial derivative. It answers "which resistor matters most in this filter?", "how fast does the potential change as I move east?", and "what is the transistor's transconductance at this bias?". Partial derivatives are also the building blocks of everything else in the chapter — the gradient, the Jacobian and the Hessian are just partial derivatives arranged in vectors and matrices.

### Definition and geometric meaning

> [!definition] Partial derivative
> The partial derivative of $f(x,y)$ with respect to $x$ at $(a,b)$ is the ordinary derivative obtained by holding $y$ fixed at $b$:
> $$
> \frac{\partial f}{\partial x}(a,b)=f_x(a,b)=\lim_{h\to0}\frac{f(a+h,b)-f(a,b)}{h}.
> $$
> $f_y$ is defined likewise with $x$ held fixed. For $f(x_1,\dots,x_n)$ there are $n$ first partial derivatives $\partial f/\partial x_i$, each obtained by varying one variable and freezing all the others.

Other notations in common use: $\partial_x f$, $D_1f$, $\left(\dfrac{\partial f}{\partial x}\right)_y$. The last one, with the held-constant variable written as a subscript, is used in thermodynamics and wherever there is real ambiguity about what is held fixed (see the pitfall below).

Geometrically, $f_x(a,b)$ is the slope of the surface $z=f(x,y)$ in the $x$-direction at $(a,b)$: slice the surface with the vertical plane $y=b$, and $f_x(a,b)$ is the slope of the resulting curve $z=f(x,b)$. Similarly $f_y(a,b)$ is the slope of the curve cut out by the plane $x=a$. The two partials therefore describe two tangent lines to the surface; together they span the tangent plane (§6.3).

**Computing.** Because $\partial f/\partial x$ is an ordinary derivative in disguise, every single-variable rule (power, product, quotient, chain) applies unchanged, treating the other variables as constants. For example, with $f(x,y)=x^2y^3+\sin(xy)$:

$$
f_x=2xy^3+y\cos(xy),\qquad f_y=3x^2y^2+x\cos(xy).
$$

The $y$ in front of $\cos(xy)$ in $f_x$ is the inner derivative $\partial(xy)/\partial x=y$ from the chain rule.

**Checking numerically.** A partial derivative can always be estimated by a **finite difference**, varying one input in a simulation or on the bench:

$$
f_x(a,b)\approx\frac{f(a+h,b)-f(a-h,b)}{2h}\qquad(\text{central difference, error }O(h^2)).
$$

At $(1,2)$ the formula above gives $f_x=16+2\cos2=15.1677$; the central difference with $h=10^{-3}$ gives $15.1677$ to all digits shown, while the one-sided difference $[f(a+h,b)-f(a,b)]/h$ gives $15.1739$, correct only to about three figures (error $O(h)$). This is the method circuit simulators fall back on for sensitivity analysis when analytic derivatives are not available, and the reason central differences are preferred ([Ch. 22](../part-7-numerical-computational/22-numerical-methods.md)).

> [!example] Worked Example 1 — partial derivatives of an RC low-pass response
> The magnitude response of a first-order RC low-pass filter is
> $$
> M(\omega,\tau)=\frac{1}{\sqrt{1+\omega^2\tau^2}},\qquad\tau=RC.
> $$
> *Partial with respect to frequency* (hold $\tau$): write $M=(1+\omega^2\tau^2)^{-1/2}$ and use the chain rule:
> $$
> \frac{\partial M}{\partial\omega}=-\tfrac12(1+\omega^2\tau^2)^{-3/2}\cdot2\omega\tau^2=-\frac{\omega\tau^2}{(1+\omega^2\tau^2)^{3/2}}.
> $$
> *Partial with respect to time constant* (hold $\omega$): by the same steps, $\partial M/\partial\tau=-\omega^2\tau/(1+\omega^2\tau^2)^{3/2}$.
>
> At the cutoff $\omega\tau=1$: $\partial M/\partial\omega=-\tau/2^{3/2}=-0.354\,\tau$ and $\partial M/\partial\tau=-\omega/2^{3/2}=-0.354\,\omega$. Units: $M$ is dimensionless, so $\partial M/\partial\omega$ must have units of seconds (it does: $\tau$) and $\partial M/\partial\tau$ units of $\mathrm{s^{-1}}$ (it does: $\omega$). The symmetry $\omega\,\partial M/\partial\omega=\tau\,\partial M/\partial\tau$ is no accident: $M$ depends only on the product $\omega\tau$.

> [!example] Worked Example 2 — sensitivity of a parallel combination
> For two resistors in parallel, $R=\dfrac{R_1R_2}{R_1+R_2}$. Holding $R_2$ fixed and using the quotient rule,
> $$
> \frac{\partial R}{\partial R_1}=\frac{R_2(R_1+R_2)-R_1R_2\cdot1}{(R_1+R_2)^2}=\frac{R_2^2}{(R_1+R_2)^2},
> $$
> and by symmetry $\partial R/\partial R_2=R_1^2/(R_1+R_2)^2$. With $R_1=1\,\mathrm{k\Omega}$, $R_2=4\,\mathrm{k\Omega}$: $R=0.8\,\mathrm{k\Omega}$, $\partial R/\partial R_1=16/25=0.64$ and $\partial R/\partial R_2=1/25=0.04$.
>
> The combination is 16 times more sensitive to the *smaller* resistor — a fact worth knowing when choosing which part needs the tighter tolerance. Both derivatives are dimensionless, as they must be ($\Omega/\Omega$). Check by finite difference: raising $R_1$ to $1.01\,\mathrm{k\Omega}$ gives $R=4.04/5.01=0.806387\,\mathrm{k\Omega}$, an increase of $6.39\ \Omega$ for a $10\ \Omega$ change, ratio $0.639$. ✓

### Normalised sensitivity

Absolute partial derivatives mix units, which makes it hard to compare the influence of, say, a capacitor and a resistor. Circuit designers therefore use the **normalised (relative) sensitivity**

$$
S^{f}_{x}=\frac{x}{f}\frac{\partial f}{\partial x}=\frac{\partial\ln f}{\partial\ln x}\approx\frac{\Delta f/f}{\Delta x/x},
$$

the percentage change in $f$ per percentage change in $x$. It is dimensionless, and for a power law $f=kx^ay^b$ it is simply the exponent: $S^f_x=a$, $S^f_y=b$.

> [!example] Worked Example 3 — sensitivities of a series RLC circuit
> A series RLC circuit has resonant frequency and quality factor
> $$
> \omega_0=(LC)^{-1/2},\qquad Q=\frac{1}{R}\sqrt{\frac{L}{C}}=R^{-1}L^{1/2}C^{-1/2}.
> $$
> Both are power laws, so the normalised sensitivities are read off the exponents:
> $$
> S^{\omega_0}_L=S^{\omega_0}_C=-\tfrac12,\quad S^{\omega_0}_R=0;\qquad S^{Q}_R=-1,\quad S^{Q}_L=+\tfrac12,\quad S^{Q}_C=-\tfrac12 .
> $$
> Verify one from the definition: $\partial\omega_0/\partial L=-\tfrac12L^{-3/2}C^{-1/2}$, so $S^{\omega_0}_L=\dfrac{L}{(LC)^{-1/2}}\cdot\left(-\tfrac12L^{-3/2}C^{-1/2}\right)=-\tfrac12$. ✓
>
> Reading: a 2 % high capacitor lowers $\omega_0$ by about 1 % and $Q$ by about 1 %; a 2 % high resistor leaves $\omega_0$ alone but lowers $Q$ by 2 %. A designer who needs an accurate centre frequency should spend money on $L$ and $C$, not on $R$.

> [!warning] Common Pitfall — what is being held constant?
> A partial derivative is only defined once you know which *other* variables are held fixed. For $P=VI$ with $I$ held fixed, $\partial P/\partial V=I$. For the same resistor written as $P=V^2/R$ with $R$ held fixed, $\partial P/\partial V=2V/R=2I$. Both are correct — they answer different questions (fixed current versus fixed resistance). When the context does not make this obvious, write $\left(\partial P/\partial V\right)_I$ or $\left(\partial P/\partial V\right)_R$.

### Higher-order partial derivatives

Each first partial is itself a function of $(x,y)$ and can be differentiated again, giving four second-order partials of $f(x,y)$:

$$
f_{xx}=\frac{\partial^2 f}{\partial x^2},\quad f_{yy}=\frac{\partial^2 f}{\partial y^2},\quad f_{xy}=\frac{\partial}{\partial y}\!\left(\frac{\partial f}{\partial x}\right)=\frac{\partial^2 f}{\partial y\,\partial x},\quad f_{yx}=\frac{\partial}{\partial x}\!\left(\frac{\partial f}{\partial y}\right)=\frac{\partial^2 f}{\partial x\,\partial y}.
$$

(The subscript notation is read left to right — $f_{xy}$ means "first $x$, then $y$" — while the fraction notation is read right to left. Fortunately the next theorem usually makes the order irrelevant.) $f_{xx}$ measures the curvature of the $x$-slice; $f_{xy}$ measures how the $x$-slope changes as you move in $y$, i.e. the *twist* of the surface.

> [!theorem] Clairaut (Schwarz) theorem
> If $f_{xy}$ and $f_{yx}$ both exist and are continuous on a neighbourhood of a point, they are equal there: the order of mixed differentiation does not matter. The same holds for higher derivatives: any two mixed partials involving the same number of differentiations in each variable agree if they are continuous.

*Proof sketch.* Consider the "second difference"

$$
\Delta(h,k)=f(a+h,b+k)-f(a+h,b)-f(a,b+k)+f(a,b).
$$

Grouping it as $g(a+h)-g(a)$ with $g(x)=f(x,b+k)-f(x,b)$ and applying the mean value theorem twice (first in $x$, then in $y$) gives $\Delta=hk\,f_{xy}(\xi,\eta)$ for some point $(\xi,\eta)$ in the small rectangle. Grouping it the other way round, as a difference in $y$ first, gives $\Delta=hk\,f_{yx}(\xi',\eta')$. Dividing by $hk$ and letting $h,k\to0$, continuity forces both mixed partials to the same value $\lim\Delta/(hk)$ at $(a,b)$. $\square$

> [!example] Worked Example 4 — second partials and a numerical check
> For $f=x^2y^3+\sin(xy)$, from the first partials above:
> $$
> \begin{aligned}
> f_{xx}&=2y^3-y^2\sin(xy), \qquad f_{yy}=6x^2y-x^2\sin(xy),\\
> f_{xy}&=\frac{\partial}{\partial y}\left[2xy^3+y\cos(xy)\right]=6xy^2+\cos(xy)-xy\sin(xy),\\
> f_{yx}&=\frac{\partial}{\partial x}\left[3x^2y^2+x\cos(xy)\right]=6xy^2+\cos(xy)-xy\sin(xy).
> \end{aligned}
> $$
> The mixed partials agree, as Clairaut promises. At $(1,2)$: $f_{xy}=24+\cos2-2\sin2=24-0.4161-1.8186=21.7653$. The four-point finite difference $\dfrac{f(a+h,b+h)-f(a+h,b-h)-f(a-h,b+h)+f(a-h,b-h)}{4h^2}$ with $h=10^{-3}$ gives $21.7653$. ✓

The hypotheses of Clairaut's theorem matter. The standard counterexample $f=xy(x^2-y^2)/(x^2+y^2)$, $f(0,0)=0$, has $f_{xy}(0,0)=-1$ but $f_{yx}(0,0)=+1$; its second partials exist near the origin but are not continuous there. Such functions do not arise from smooth physical models, which is why engineers can (and do) swap the order of differentiation freely.

> [!warning] Common Pitfall — partial derivatives do not guarantee continuity
> In one variable, a differentiable function is continuous. In several variables, the mere *existence* of partial derivatives guarantees nothing of the kind. The function $xy/(x^2+y^2)$ of §6.1 (with value 0 at the origin) has $f_x(0,0)=f_y(0,0)=0$, because it is identically zero along both axes — yet it is not even continuous at the origin. Partial derivatives only probe two directions. The right multivariable notion of "differentiable" is the subject of §6.3.

### Partial differential equations

Equations relating partial derivatives are **partial differential equations** (PDEs), the subject of [Ch. 13](../part-4-differential-equations/13-partial-differential-equations.md). Here we only need to *check* that a given function satisfies one, which is pure partial differentiation. The three classical equations are:

| Equation | Form (2-D, or 1-D plus time) | Typical physics |
|---|---|---|
| Laplace | $u_{xx}+u_{yy}=0$ | electrostatic potential in a charge-free region |
| Heat (diffusion) | $u_t=k\,u_{xx}$ | heat flow, carrier diffusion, RC transmission lines |
| Wave | $u_{tt}=c^2u_{xx}$ | voltage on a lossless line, sound, light |

> [!example] Worked Example 5 — checking Laplace's equation
> The potential of an infinite line charge, up to constants, is $u=\ln(x^2+y^2)$. Then
> $$
> u_x=\frac{2x}{x^2+y^2},\qquad u_{xx}=\frac{2(x^2+y^2)-2x\cdot2x}{(x^2+y^2)^2}=\frac{2y^2-2x^2}{(x^2+y^2)^2},
> $$
> and by the symmetry $x\leftrightarrow y$, $u_{yy}=\dfrac{2x^2-2y^2}{(x^2+y^2)^2}$. Hence $u_{xx}+u_{yy}=0$ away from the origin: $u$ satisfies **Laplace's equation**, as every electrostatic potential in a charge-free region must. (At the origin, where the line charge sits, $u$ is not even defined.) A function satisfying Laplace's equation is called **harmonic**.

> [!example] Worked Example 6 — travelling waves and decaying modes
> (a) Let $u(x,t)=g(x-ct)$ for any twice-differentiable $g$. With $s=x-ct$, the chain rule gives $u_x=g'(s)$, $u_{xx}=g''(s)$, $u_t=-c\,g'(s)$, $u_{tt}=c^2g''(s)$. Hence $u_{tt}=c^2u_{xx}$: **any** shape moving at speed $c$ without distortion solves the wave equation. This is why a lossless transmission line carries pulses of any shape intact.
>
> (b) Let $u(x,t)=e^{-k\beta^2t}\sin(\beta x)$. Then $u_t=-k\beta^2u$ and $u_{xx}=-\beta^2u$, so $u_t=k\,u_{xx}$: a sinusoidal temperature (or charge) profile of spatial frequency $\beta$ decays exponentially with time constant $1/(k\beta^2)$. Finer ripples (large $\beta$) die out faster, which is why diffusion smooths things. Units: $k$ is in $\mathrm{m^2/s}$ and $\beta$ in $\mathrm{m^{-1}}$, so $k\beta^2$ is in $\mathrm{s^{-1}}$. ✓

> [!tip] Engineering Note
> Partial derivatives of device equations *are* the small-signal parameters on a datasheet. A MOSFET in saturation with $i_D=\tfrac12k(v_{GS}-V_t)^2(1+\lambda v_{DS})$ has transconductance $g_m=\partial i_D/\partial v_{GS}=k(v_{GS}-V_t)(1+\lambda v_{DS})$ and output conductance $g_{ds}=\partial i_D/\partial v_{DS}=\tfrac12k(v_{GS}-V_t)^2\lambda$, each evaluated at the bias point with the other voltage held fixed. Section 6.6 assembles such derivatives into a matrix.

## 6.3 Total differentials

### Motivation

Partial derivatives answer "what if I change one input?". In practice everything changes at once: every component is off-nominal by its own tolerance, temperature shifts several parameters together, a sensor reading depends on several imperfect measurements. The **total differential** answers "what happens when all the inputs change a little, simultaneously?". It is the multivariable version of the linear approximation $\Delta y\approx f'(x)\,\Delta x$ from Chapter 4, and it leads directly to tolerance analysis, the multivariable chain rule and implicit differentiation.

### Differentiability and the tangent plane

> [!definition] Differentiability and total differential
> $f(x,y)$ is **differentiable** at $(a,b)$ if there are constants $A$, $B$ such that
> $$
> f(a+\Delta x,\,b+\Delta y)=f(a,b)+A\,\Delta x+B\,\Delta y+\varepsilon\sqrt{\Delta x^2+\Delta y^2},
> $$
> with $\varepsilon\to0$ as $(\Delta x,\Delta y)\to(0,0)$. Then necessarily $A=f_x(a,b)$ and $B=f_y(a,b)$, and the linear part
> $$
> df=\frac{\partial f}{\partial x}\,dx+\frac{\partial f}{\partial y}\,dy
> $$
> is called the **total differential** of $f$. For $n$ variables, $df=\sum_{i=1}^n \dfrac{\partial f}{\partial x_i}\,dx_i$.

In words: near a point of differentiability, $f$ is a plane plus an error that shrinks *faster than the distance* from the point. Setting $\Delta y=0$ and dividing by $\Delta x$ shows $A=f_x$; similarly $B=f_y$. So the coefficients are forced to be the partial derivatives — differentiability is the extra requirement that this plane approximates $f$ well in **all** directions, not just along the axes.

> [!theorem] A practical test for differentiability
> If $f_x$ and $f_y$ exist near $(a,b)$ and are continuous at $(a,b)$, then $f$ is differentiable at $(a,b)$. A differentiable function is continuous.

*Proof sketch.* Write $f(a+\Delta x,b+\Delta y)-f(a,b)=[f(a+\Delta x,b+\Delta y)-f(a,b+\Delta y)]+[f(a,b+\Delta y)-f(a,b)]$. Apply the mean value theorem to each bracket: the first is $f_x(\xi,b+\Delta y)\Delta x$, the second $f_y(a,\eta)\Delta y$, with $\xi$ between $a$ and $a+\Delta x$ and $\eta$ between $b$ and $b+\Delta y$. Continuity of the partials lets us replace $f_x(\xi,b+\Delta y)$ by $f_x(a,b)$ and $f_y(a,\eta)$ by $f_y(a,b)$ with errors that tend to zero, which is exactly the definition. $\square$

All the functions in this book's applications have continuous partial derivatives wherever they are defined, so in practice "differentiable" means "the partials exist and are continuous", and functions such as $xy/(x^2+y^2)$ at the origin are exotic exceptions.

Geometrically, the linear part defines the **tangent plane** to the surface $z=f(x,y)$ at $(a,b,f(a,b))$:

$$
z=f(a,b)+f_x(a,b)\,(x-a)+f_y(a,b)\,(y-b).
$$

The **linearisation** $L(x,y)$ is the right-hand side; near $(a,b)$, $f\approx L$.

> [!example] Worked Example 1 — linear approximation of a phasor magnitude
> An impedance $Z=R+jX$ has magnitude $\lvert Z\rvert=f(R,X)=\sqrt{R^2+X^2}$. At $R=3\ \Omega$, $X=4\ \Omega$, $\lvert Z\rvert=5\ \Omega$. Estimate $\lvert Z\rvert$ when $R=3.02\ \Omega$ and $X=3.97\ \Omega$.
>
> *Partials:* $f_R=R/\sqrt{R^2+X^2}=3/5=0.6$, $f_X=X/\sqrt{R^2+X^2}=4/5=0.8$ (both dimensionless).
>
> *Differential:* $df=0.6(0.02)+0.8(-0.03)=0.012-0.024=-0.012\ \Omega$.
>
> *Estimate:* $\lvert Z\rvert\approx5-0.012=4.988\ \Omega$. The exact value is $\sqrt{3.02^2+3.97^2}=4.98812\ \Omega$; the linearisation error, $1.2\times10^{-4}\ \Omega$, is of second order in the changes, as the definition promises.
>
> Note the pattern: $f_R=\cos\angle Z$ and $f_X=\sin\angle Z$. The magnitude responds to $R$ and $X$ in proportion to how much each contributes to the phasor's direction.

### Error and tolerance propagation

Interpreting $dx_i$ as small errors or tolerances $\Delta x_i$ turns the total differential into an error-propagation formula. Two estimates are standard:

$$
\Delta f_{\text{worst}}=\sum_i\left\lvert\frac{\partial f}{\partial x_i}\right\rvert\Delta x_i,
\qquad
\Delta f_{\text{RSS}}=\sqrt{\sum_i\left(\frac{\partial f}{\partial x_i}\Delta x_i\right)^2}.
$$

The **worst-case** sum assumes every error conspires in the same direction; it is guaranteed (to first order) and is what you use for safety-critical limits. The **root-sum-square** (RSS) value treats the errors as independent random variables with standard deviations $\Delta x_i$; the variance of a linear combination is the sum of the squared-weighted variances ([Ch. 19](../part-6-probability-statistics/19-probability.md), [Ch. 20](../part-6-probability-statistics/20-statistical-analysis.md)), so RSS gives the *typical* spread. RSS is never larger than worst case, and the gap widens as the number of contributing parameters grows.

For products of powers, $f=kx^ay^b z^c$, dividing $df$ by $f$ gives the handy **relative-error rule**

$$
\frac{df}{f}=a\,\frac{dx}{x}+b\,\frac{dy}{y}+c\,\frac{dz}{z},
$$

whose coefficients are the normalised sensitivities of §6.2. In words: for products and quotients, *relative* errors add, weighted by the exponents.

> [!example] Worked Example 2 — power measurement uncertainty
> A resistor $R=100\pm1\ \Omega$ carries $V=12\pm0.1\ \mathrm{V}$. Find $P=V^2/R$ and its uncertainty.
>
> *Nominal:* $P=144/100=1.44\ \mathrm{W}$.
>
> *Differential:*
> $$
> dP=\frac{2V}{R}\,dV-\frac{V^2}{R^2}\,dR=(0.24\ \mathrm{A})\,dV-(0.0144\ \mathrm{A^2})\,dR .
> $$
> *Worst case:* $\Delta P=0.24(0.1)+0.0144(1)=0.024+0.0144=0.0384\ \mathrm{W}$, i.e. $2.67\,\%$. Check with the relative rule: $2(0.1/12)+1(1/100)=0.0167+0.0100=0.0267$. ✓
>
> *RSS:* $\Delta P=\sqrt{0.024^2+0.0144^2}=0.0280\ \mathrm{W}$ ($1.94\,\%$).
>
> The voltage term dominates because $V$ enters squared: improving the voltmeter is worth more than buying a better resistor. Unit check: $(\mathrm{V}/\Omega)\cdot\mathrm{V}=\mathrm{W}$ and $(\mathrm{V^2}/\Omega^2)\cdot\Omega=\mathrm{W}$. ✓

> [!example] Worked Example 3 — tolerance of an RC cutoff frequency
> A low-pass filter uses $R=10\ \mathrm{k\Omega}\pm1\,\%$ and $C=10\ \mathrm{nF}\pm5\,\%$. Its cutoff is $f_c=1/(2\pi RC)=1591.5\ \mathrm{Hz}$.
>
> Since $f_c\propto R^{-1}C^{-1}$, the relative rule gives $\dfrac{df_c}{f_c}=-\dfrac{dR}{R}-\dfrac{dC}{C}$.
>
> *Worst case:* $1\,\%+5\,\%=6\,\%$, so $f_c=1591.5\pm95.5\ \mathrm{Hz}$.
>
> *RSS:* $\sqrt{1^2+5^2}\,\%=5.10\,\%$, so $\pm81.2\ \mathrm{Hz}$.
>
> The capacitor tolerance is almost the whole story: replacing the 1 % resistor by a 0.1 % one would change the RSS figure only from 5.10 % to 5.00 %. The exact worst case, $1/(2\pi\cdot0.99R\cdot0.95C)$, is $6.33\,\%$ high; the first-order estimate (6 %) is slightly optimistic because the true dependence is curved — differentials are accurate only for *small* changes.

> [!warning] Common Pitfall
> Two mistakes recur in tolerance work. First, adding signed terms: in Example 2 the two contributions have opposite signs, and $0.024-0.0144$ would grossly understate the uncertainty; take absolute values (worst case) or squares (RSS). Second, using RSS for errors that are *not* independent — two resistors from the same reel, or two parameters that both drift with temperature, are correlated, and RSS can then badly underestimate the spread.

### The multivariable chain rule

Suppose $z=f(x,y)$ and the inputs themselves depend on time: $x=x(t)$, $y=y(t)$. Dividing the total differential by $dt$ suggests the following, and the differentiability definition makes it rigorous.

> [!theorem] Chain rule
> If $f$ is differentiable and $x(t)$, $y(t)$ are differentiable, then
> $$
> \frac{dz}{dt}=\frac{\partial f}{\partial x}\frac{dx}{dt}+\frac{\partial f}{\partial y}\frac{dy}{dt}.
> $$
> If instead $x=x(u,v)$ and $y=y(u,v)$, then
> $$
> \frac{\partial z}{\partial u}=\frac{\partial f}{\partial x}\frac{\partial x}{\partial u}+\frac{\partial f}{\partial y}\frac{\partial y}{\partial u},\qquad
> \frac{\partial z}{\partial v}=\frac{\partial f}{\partial x}\frac{\partial x}{\partial v}+\frac{\partial f}{\partial y}\frac{\partial y}{\partial v}.
> $$

*Proof sketch.* From the definition of differentiability, $\Delta z=f_x\Delta x+f_y\Delta y+\varepsilon\sqrt{\Delta x^2+\Delta y^2}$. Divide by $\Delta t$ and let $\Delta t\to0$: $\Delta x/\Delta t\to\dot x$, $\Delta y/\Delta t\to\dot y$, the square root over $\lvert\Delta t\rvert$ stays bounded (it tends to $\sqrt{\dot x^2+\dot y^2}$), and $\varepsilon\to0$. $\square$

A **tree diagram** keeps the terms straight: draw $z$ at the top, branches to each intermediate variable ($x$, $y$), and from each of those branches down to the independent variables. Each path from $z$ to the variable of interest contributes the product of the derivatives along it; add over all paths.

> [!example] Worked Example 4 — rate of change of dissipated power
> A resistor's voltage rises at $0.5\ \mathrm{V/s}$ while self-heating raises its resistance at $2\ \Omega/\mathrm{s}$. At the instant $V=12\ \mathrm{V}$, $R=100\ \Omega$, how fast is $P=V^2/R$ changing?
> $$
> \frac{dP}{dt}=\frac{\partial P}{\partial V}\dot V+\frac{\partial P}{\partial R}\dot R=\frac{2V}{R}\dot V-\frac{V^2}{R^2}\dot R=0.24(0.5)-0.0144(2)=0.120-0.0288=0.0912\ \mathrm{W/s}.
> $$
> The rising voltage increases the power; the rising resistance reduces it; the net effect is the sum of both paths of the tree.

> [!example] Worked Example 5 — Cartesian to polar derivatives
> Let $f(x,y)=x^2+xy$ and change to polar coordinates, $x=r\cos\theta$, $y=r\sin\theta$. By the chain rule,
> $$
> \frac{\partial f}{\partial r}=f_x\cos\theta+f_y\sin\theta,\qquad
> \frac{\partial f}{\partial\theta}=-f_x\,r\sin\theta+f_y\,r\cos\theta .
> $$
> At $r=2$, $\theta=\pi/6$: $x=\sqrt3$, $y=1$, $f_x=2x+y=2\sqrt3+1=4.4641$, $f_y=x=1.7321$. Then
> $$
> \frac{\partial f}{\partial r}=4.4641(0.8660)+1.7321(0.5)=4.7321,\qquad
> \frac{\partial f}{\partial\theta}=-4.4641(1)+1.7321(1.7321)=-1.4641 .
> $$
> *Check directly:* in polar form $f=r^2(\cos^2\theta+\sin\theta\cos\theta)$, so $f_r=2r(\cos^2\theta+\sin\theta\cos\theta)=4(0.75+0.4330)=4.7321$ and $f_\theta=r^2(\cos2\theta-\sin2\theta)=4(0.5-0.8660)=-1.4641$. ✓
>
> Note that $f_\theta$ has the units of $f$ (it is a rate per radian), whereas $f_r$ has the units of $f$ per metre; the factor $r$ in the $\theta$-equation carries the difference. This observation is the seed of the scale factors in §6.9.

> [!warning] Common Pitfall — cancelling partials
> In one variable, $\frac{dz}{dx}\frac{dx}{dt}$ "cancels" to $\frac{dz}{dt}$. With partial derivatives the cancellation is wrong: $\frac{\partial z}{\partial x}\frac{\partial x}{\partial u}$ is only *one* of the terms in $\partial z/\partial u$. Likewise $\partial x/\partial r=\cos\theta$, but $\partial r/\partial x=x/r=\cos\theta$ too, not $1/\cos\theta$ — the reciprocal rule fails because different variables are held constant in the two derivatives ($\theta$ in one, $y$ in the other).

### Implicit differentiation

Often a relationship between variables is given implicitly, $F(x,y)=0$, and cannot (or need not) be solved for $y$. Along the curve $F$ is constant, so its total differential vanishes: $F_x\,dx+F_y\,dy=0$, and therefore

$$
\frac{dy}{dx}=-\frac{F_x}{F_y}\qquad (F_y\ne0).
$$

With three variables and $F(x,y,z)=0$ defining $z=z(x,y)$, holding $y$ fixed gives $F_x+F_z\,\partial z/\partial x=0$, so

$$
\frac{\partial z}{\partial x}=-\frac{F_x}{F_z},\qquad \frac{\partial z}{\partial y}=-\frac{F_y}{F_z}\qquad(F_z\ne0).
$$

> [!theorem] Implicit function theorem (two-variable form)
> If $F$ has continuous partial derivatives near $(a,b)$, $F(a,b)=0$ and $F_y(a,b)\ne0$, then near $a$ the equation $F(x,y)=0$ defines a unique differentiable function $y=g(x)$ with $g(a)=b$, and $g'(x)=-F_x/F_y$.

The condition $F_y\ne0$ is the important part: where $F_y=0$ the curve has a vertical tangent (or a self-intersection) and $y$ cannot be written as a single-valued function of $x$ there. For example, on $x^2+xy+y^2=7$ at $(1,2)$: $F_x=2x+y=4$, $F_y=x+2y=5$, so $dy/dx=-4/5$.

> [!example] Worked Example 6 — small-signal conductance of a real diode
> A diode with series resistance obeys
> $$
> F(I,V)=I-I_s\left(e^{(V-IR_s)/V_T}-1\right)=0,
> $$
> which cannot be solved for $I$ in elementary functions. Implicit differentiation still gives the small-signal conductance exactly. With $E=e^{(V-IR_s)/V_T}$:
> $$
> F_V=-\frac{I_sE}{V_T},\qquad F_I=1+\frac{I_sER_s}{V_T},\qquad
> \frac{dI}{dV}=-\frac{F_V}{F_I}=\frac{I_sE/V_T}{1+R_sI_sE/V_T}.
> $$
> From the equation itself, $I_sE=I+I_s$. Writing $r_d=V_T/(I+I_s)$ (the junction's own small-signal resistance), this simplifies to
> $$
> \frac{dI}{dV}=\frac{1}{r_d+R_s}.
> $$
> At $I=1\ \mathrm{mA}$ with $V_T=25\ \mathrm{mV}$ and $R_s=10\ \Omega$ ($I_s$ negligible): $r_d=25\ \Omega$ and $dI/dV=1/35\ \mathrm{S}=28.6\ \mathrm{mS}$. The answer is physically obvious in hindsight — the junction resistance and the series resistance add — but the calculus delivered it without ever solving for $I$.

### Exact differentials

The total differential of any $\phi(x,y)$ has the form $M\,dx+N\,dy$ with $M=\phi_x$, $N=\phi_y$. Conversely, given an expression $M\,dx+N\,dy$, is it the differential of something? If so, $M_y=\phi_{xy}=\phi_{yx}=N_x$ by Clairaut's theorem. The converse also holds on regions without holes:

> [!theorem] Test for an exact differential
> On a simply connected domain, with $M$, $N$ continuously differentiable, $M\,dx+N\,dy$ is the total differential of some function $\phi$ if and only if
> $$
> \frac{\partial M}{\partial y}=\frac{\partial N}{\partial x}.
> $$

> [!example] Worked Example 7 — recovering a potential function
> Is $(2xy+3)\,dx+(x^2+4y)\,dy$ exact? $M_y=2x$, $N_x=2x$: yes. To find $\phi$, integrate $M$ in $x$ treating $y$ as constant: $\phi=x^2y+3x+h(y)$. Differentiate in $y$ and match $N$: $x^2+h'(y)=x^2+4y$, so $h'=4y$, $h=2y^2$. Thus $\phi=x^2y+3x+2y^2$ (plus any constant). Check: $\phi_x=2xy+3$, $\phi_y=x^2+4y$. ✓

This test reappears as "curl-free implies conservative" in [Ch. 7](../part-2-calculus/07-vector-calculus.md) — the reason a static electric field can be written as the gradient of a potential — and in exact first-order ODEs in [Ch. 11](../part-4-differential-equations/11-ordinary-differential-equations.md).

### Taylor's theorem in two variables

The tangent plane is the first-order Taylor approximation. Adding second-order terms captures curvature, which we will need to classify maxima and minima in §6.10. Applying the one-variable Taylor series to $g(t)=f(a+th,\,b+tk)$ and using the chain rule twice gives

$$
f(a+h,b+k)\approx f+f_x h+f_y k+\tfrac12\left(f_{xx}h^2+2f_{xy}hk+f_{yy}k^2\right),
$$

with all derivatives evaluated at $(a,b)$. The bracket is a **quadratic form**, $\begin{bmatrix}h&k\end{bmatrix}\mathbf{H}\begin{bmatrix}h\\k\end{bmatrix}$, built from the **Hessian matrix** $\mathbf{H}$ of second partials.

> [!example] Worked Example 8 — second-order approximation
> Approximate $f=e^x\cos y$ near $(0,0)$. All needed derivatives at the origin: $f=1$, $f_x=e^x\cos y=1$, $f_y=-e^x\sin y=0$, $f_{xx}=1$, $f_{xy}=-e^x\sin y=0$, $f_{yy}=-e^x\cos y=-1$. Hence
> $$
> e^x\cos y\approx1+x+\tfrac12\left(x^2-y^2\right).
> $$
> At $(0.1,0.2)$: first order gives $1.1$, second order $1.085$, exact $1.08314$. The remaining error, $0.0019$, matches the third-order terms $x^3/6-xy^2/2=-0.0018$. Each order improves the estimate by roughly a factor of the step size, just as in one variable.

## 6.4 Gradient

### Motivation

The total differential $df=f_x\,dx+f_y\,dy+f_z\,dz$ has the shape of a dot product: one vector made of partial derivatives, the other made of the displacement $(dx,dy,dz)$. Packaging the partial derivatives into a single vector — the gradient — turns a list of numbers into a geometric object with a direction and a length. That direction is the direction of steepest ascent, it is perpendicular to the level sets, and in electrostatics it *is* (minus) the electric field. The gradient is the first of the three differential operators of vector calculus, and the only one that acts on a scalar field.

### Definition

> [!definition] Gradient
> For a differentiable scalar function $f(x,y,z)$, the gradient is the vector of its first partial derivatives
> $$
> \nabla f=\frac{\partial f}{\partial x}\hat{\mathbf{x}}+\frac{\partial f}{\partial y}\hat{\mathbf{y}}+\frac{\partial f}{\partial z}\hat{\mathbf{z}}=\left(f_x,\ f_y,\ f_z\right).
> $$
> For $f(x_1,\dots,x_n)$, $\nabla f=(\partial f/\partial x_1,\dots,\partial f/\partial x_n)$. The symbol $\nabla$ ("del" or "nabla") is the operator $\hat{\mathbf{x}}\,\partial_x+\hat{\mathbf{y}}\,\partial_y+\hat{\mathbf{z}}\,\partial_z$.

With $d\mathbf{r}=(dx,dy,dz)$, the total differential becomes

$$
df=\nabla f\cdot d\mathbf{r},
$$

and the linearisation of §6.3 becomes $f(\mathbf{r})\approx f(\mathbf{a})+\nabla f(\mathbf{a})\cdot(\mathbf{r}-\mathbf{a})$. **Units:** $\nabla$ carries $\mathrm{m^{-1}}$, so the gradient of a temperature is in $\mathrm{K/m}$ and of a potential in $\mathrm{V/m}$.

### Geometric properties

> [!theorem] Properties of the gradient
> Let $f$ be differentiable at $\mathbf{a}$ with $\nabla f(\mathbf{a})\ne\mathbf{0}$.
> 1. **Steepest ascent.** $\nabla f(\mathbf{a})$ points in the direction in which $f$ increases most rapidly, and $\lVert\nabla f(\mathbf{a})\rVert$ is that maximum rate of increase (per unit distance).
> 2. **Normal to level sets.** $\nabla f(\mathbf{a})$ is perpendicular to the level curve (or level surface) $f=f(\mathbf{a})$ through $\mathbf{a}$.
> 3. **Tangent plane.** The tangent plane to the level surface $f(x,y,z)=c$ at $\mathbf{a}$ is $\nabla f(\mathbf{a})\cdot(\mathbf{r}-\mathbf{a})=0$.

*Proof.* (1) Move a small distance $ds$ in the unit direction $\hat{\mathbf{u}}$: $df=\nabla f\cdot\hat{\mathbf{u}}\,ds=\lVert\nabla f\rVert\cos\theta\,ds$, where $\theta$ is the angle between $\hat{\mathbf{u}}$ and $\nabla f$. This is largest when $\theta=0$, giving rate $\lVert\nabla f\rVert$. (2) Let $\mathbf{r}(t)$ be any curve lying in the level set, through $\mathbf{a}$ at $t=0$. Since $f(\mathbf{r}(t))=c$ for all $t$, the chain rule gives $0=\frac{d}{dt}f(\mathbf{r}(t))=\nabla f\cdot\mathbf{r}'(0)$. So $\nabla f$ is perpendicular to the tangent vector of *every* curve in the level set, i.e. to the level set itself. (3) The tangent plane consists of points $\mathbf{r}$ with $\mathbf{r}-\mathbf{a}$ perpendicular to the normal $\nabla f(\mathbf{a})$. $\square$

Property 2 explains a picture every electrical engineer knows: field lines cross equipotentials at right angles. Property 1 explains why field strength is greatest where equipotentials crowd together — closely spaced contours mean a large change of $f$ per unit distance.

> [!example] Worked Example 1 — normal and tangent to a level curve
> Let $f(x,y)=x^2+2y^2$ and take the point $(1,1)$, which lies on the level curve $f=3$ (an ellipse). Then $\nabla f=(2x,4y)=(2,4)$ at $(1,1)$.
>
> *Tangent line:* $\nabla f\cdot(x-1,\,y-1)=0$, i.e. $2(x-1)+4(y-1)=0$, or $x+2y=3$.
>
> *Normal line* (along $\nabla f$): $(x,y)=(1,1)+t(2,4)$, i.e. $y-1=2(x-1)$.
>
> *Maximum rate of increase* at $(1,1)$: $\lVert(2,4)\rVert=\sqrt{20}\approx4.47$ units of $f$ per unit distance, in the direction $(2,4)/\sqrt{20}=(0.447,\,0.894)$.
>
> *Check by implicit differentiation (§6.3):* $dy/dx=-f_x/f_y=-2/4=-\tfrac12$, which is the slope of $x+2y=3$. ✓

```desmos-graph
left=-3; right=3; bottom=-2; top=2
---
x^2+2y^2=1|grey|dashed
x^2+2y^2=3|blue
x^2+2y^2=6|grey|dashed
x+2y=3|green
y=1+2(x-1)|red|0.6<x<1.6
(1,1)|label:(1,1)
```

*The gradient (red, along the normal) crosses each level curve at right angles; the tangent line (green) is perpendicular to it.*

> [!example] Worked Example 2 — tangent plane to a level surface
> Find the tangent plane and normal line to the ellipsoid $\dfrac{x^2}{4}+y^2+\dfrac{z^2}{9}=3$ at $(2,1,3)$.
>
> *Check the point:* $1+1+1=3$. ✓
>
> *Gradient:* with $F=\tfrac14x^2+y^2+\tfrac19z^2$, $\nabla F=\left(\tfrac{x}{2},\ 2y,\ \tfrac{2z}{9}\right)=\left(1,\ 2,\ \tfrac23\right)$.
>
> *Tangent plane:* $1(x-2)+2(y-1)+\tfrac23(z-3)=0$; multiplying by 3, $3x+6y+2z=18$. The point satisfies it: $6+6+6=18$. ✓
>
> *Normal line:* $(x,y,z)=(2,1,3)+t\left(1,2,\tfrac23\right)$.
>
> The same recipe gives the tangent plane to a graph $z=f(x,y)$: write it as the level surface $F=f(x,y)-z=0$, whose gradient $(f_x,f_y,-1)$ reproduces the tangent-plane formula of §6.3.

### Algebra of the gradient

All of the following follow at once from the corresponding single-variable rules applied component by component:

$$
\nabla(af+bg)=a\nabla f+b\nabla g,\qquad \nabla(fg)=f\nabla g+g\nabla f,\qquad \nabla\!\left(\frac fg\right)=\frac{g\nabla f-f\nabla g}{g^2},\qquad \nabla h(f)=h'(f)\,\nabla f .
$$

The last one (chain rule) is especially useful for functions of distance. With $\mathbf{R}=(x,y,z)$ and $R=\lVert\mathbf{R}\rVert=\sqrt{x^2+y^2+z^2}$, $\partial R/\partial x=x/R$ and so on, giving

$$
\nabla R=\frac{\mathbf{R}}{R}=\hat{\mathbf{R}},\qquad \nabla R^n=nR^{n-1}\hat{\mathbf{R}},\qquad \nabla\frac1R=-\frac{\hat{\mathbf{R}}}{R^2}.
$$

Any function of distance alone has a gradient pointing radially, and its level sets are spheres — the geometry of every point source.

> [!example] Worked Example 3 — Coulomb's law from a potential
> The potential of a point charge $q$ at the origin is $V=\dfrac{q}{4\pi\varepsilon_0R}$. Using $\nabla(1/R)=-\hat{\mathbf{R}}/R^2$,
> $$
> \mathbf{E}=-\nabla V=-\frac{q}{4\pi\varepsilon_0}\nabla\frac1R=\frac{q}{4\pi\varepsilon_0R^2}\hat{\mathbf{R}},
> $$
> which is Coulomb's law. Units: $V$ in volts, $\nabla$ in $\mathrm{m^{-1}}$, so $\mathbf{E}$ in $\mathrm{V/m}$. ✓ Component check at $(x,y,z)=(1,2,2)\ \mathrm{m}$, $R=3\ \mathrm{m}$: $\partial(1/R)/\partial x=-x/R^3=-1/27$, matching $-\hat{\mathbf{R}}_x/R^2=-(1/3)/9=-1/27$. ✓

> [!example] Worked Example 4 — field of the dipole on its bisector
> For the 2-D dipole potential of §6.1, $V=\dfrac1{r_1}-\dfrac1{r_2}$ with $r_1=\sqrt{(x+1)^2+y^2}$ and $r_2=\sqrt{(x-1)^2+y^2}$, find $\mathbf{E}=-\nabla V$ at $(0,1)$.
>
> Using $\partial(1/r_1)/\partial x=-(x+1)/r_1^3$, $\partial(1/r_1)/\partial y=-y/r_1^3$, and similarly for $r_2$:
> $$
> V_x=-\frac{x+1}{r_1^3}+\frac{x-1}{r_2^3},\qquad V_y=-\frac{y}{r_1^3}+\frac{y}{r_2^3}.
> $$
> At $(0,1)$, $r_1=r_2=\sqrt2$, $r^3=2^{3/2}=2.828$: $V_x=-\dfrac{1}{2.828}-\dfrac{1}{2.828}=-0.707$, $V_y=0$. Hence $\mathbf{E}=(0.707,\,0)$.
>
> The field on the bisector points in $+x$, from the positive charge towards the negative one, and parallel to the axis — perpendicular to the $V=0$ equipotential $x=0$, exactly as property 2 requires.

### The gradient in other coordinates (preview)

The definition $df=\nabla f\cdot d\mathbf{r}$ is coordinate-free, so it can be used to *derive* the gradient in any coordinate system. In polar coordinates a small displacement is $d\mathbf{r}=dr\,\hat{\mathbf{r}}+r\,d\theta\,\hat{\boldsymbol{\theta}}$ (moving $d\theta$ in angle covers a distance $r\,d\theta$), while $df=f_r\,dr+f_\theta\,d\theta$. Matching the two forms,

$$
\nabla f=\frac{\partial f}{\partial r}\hat{\mathbf{r}}+\frac1r\frac{\partial f}{\partial\theta}\hat{\boldsymbol{\theta}} .
$$

The factor $1/r$ is forced by the units: $f_\theta$ is per radian, the gradient must be per metre. The cylindrical and spherical versions follow the same way from the scale factors of §6.9 and are tabulated in [Ch. 7](../part-2-calculus/07-vector-calculus.md).

> [!warning] Common Pitfall
> Writing $\nabla f=(f_r,\ f_\theta)$ in polar coordinates — simply listing partial derivatives — is wrong, for the units reason just given. The "list of partials" recipe works only in Cartesian coordinates, where every coordinate is a length.

### Gradient descent

Since $-\nabla f$ points in the direction of steepest *decrease*, a natural minimisation algorithm is to step repeatedly downhill:

$$
\mathbf{x}_{k+1}=\mathbf{x}_k-\alpha\,\nabla f(\mathbf{x}_k),
$$

with a step size $\alpha>0$. This is the core of adaptive filtering (the LMS algorithm), neural-network training and many circuit-optimisation tools.

> [!example] Worked Example 5 — gradient descent and its step-size limit
> Minimise $f(x,y)=x^2+4y^2$ (minimum $0$ at the origin) from $(2,1)$, where $f=8$. Here $\nabla f=(2x,\,8y)$, so one step is
> $$
> x_{k+1}=(1-2\alpha)x_k,\qquad y_{k+1}=(1-8\alpha)y_k .
> $$
> *With $\alpha=0.1$:* $x$ shrinks by $0.8$ and $y$ by $0.2$ per step: $(2,1)\to(1.6,0.2)\to(1.28,0.04)\to(1.024,0.008)$, with $f=8\to2.72\to1.645\to1.049$. The steep $y$-direction converges fast; the shallow $x$-direction slowly.
>
> *With $\alpha=0.3$:* the $y$ factor is $1-2.4=-1.4$: $y$ overshoots and grows, $1\to-1.4\to1.96\to-2.744$, and $f$ *increases*, $8\to8.48\to15.47\to30.13$. The iteration diverges.
>
> Stability requires $\lvert1-8\alpha\rvert<1$, i.e. $\alpha<0.25$, set by the *steepest* direction, while convergence speed is limited by the *shallowest* one. The ratio of curvatures ($8/2=4$ here) — the condition number of the Hessian — governs how badly gradient descent zig-zags ([Ch. 23](../part-7-numerical-computational/23-optimization.md)).

```desmos-graph
left=-3; right=3; bottom=-2; top=2
---
x^2+4y^2=8|grey
x^2+4y^2=2.72|grey
x^2+4y^2=1.6448|grey
(2,1)|blue|label:start
(1.6,0.2)|blue
(1.28,0.04)|blue
(1.024,0.008)|blue
(0.8,-1.4)|red
(0.32,1.96)|red
```

*Gradient descent on $x^2+4y^2$: with $\alpha=0.1$ (blue) the iterates slide down towards the minimum, each step perpendicular to the contour it leaves; with $\alpha=0.3$ (red) they bounce across the valley onto ever higher contours.*

## 6.5 Directional derivatives

### Motivation

Partial derivatives measure change along the coordinate axes, but the axes are an arbitrary human choice. A probe moving diagonally across a circuit board, a satellite crossing a radiation pattern at an angle, an optimiser deciding which way to adjust two component values at once — each needs the rate of change in a *chosen* direction. The directional derivative supplies it, and the remarkable fact is that the gradient alone determines it for every direction at once.

### Definition and the gradient formula

> [!definition] Directional derivative
> The rate of change of $f$ at $\mathbf{a}$ in the direction of a **unit** vector $\hat{\mathbf{u}}$ is
> $$
> D_{\hat{\mathbf{u}}}f(\mathbf{a})=\lim_{h\to0}\frac{f(\mathbf{a}+h\hat{\mathbf{u}})-f(\mathbf{a})}{h}.
> $$
> It is the slope of the curve obtained by slicing the graph of $f$ with the vertical plane through $\mathbf{a}$ in the direction $\hat{\mathbf{u}}$. The partial derivatives are the special cases $\hat{\mathbf{u}}=\hat{\mathbf{x}},\hat{\mathbf{y}},\hat{\mathbf{z}}$.

> [!theorem] Directional derivative from the gradient
> If $f$ is differentiable at $\mathbf{a}$, then for every unit vector $\hat{\mathbf{u}}$
> $$
> D_{\hat{\mathbf{u}}}f(\mathbf{a})=\nabla f(\mathbf{a})\cdot\hat{\mathbf{u}}=\lVert\nabla f(\mathbf{a})\rVert\cos\theta,
> $$
> where $\theta$ is the angle between $\hat{\mathbf{u}}$ and $\nabla f(\mathbf{a})$. In two dimensions, with $\hat{\mathbf{u}}=(\cos\alpha,\sin\alpha)$, $D_{\hat{\mathbf{u}}}f=f_x\cos\alpha+f_y\sin\alpha$.

*Proof.* Let $g(h)=f(\mathbf{a}+h\hat{\mathbf{u}})$, a function of one variable. By definition $D_{\hat{\mathbf{u}}}f(\mathbf{a})=g'(0)$. The point $\mathbf{a}+h\hat{\mathbf{u}}$ has coordinates $x_i=a_i+hu_i$, so $dx_i/dh=u_i$, and the chain rule of §6.3 gives $g'(0)=\sum_i f_{x_i}(\mathbf{a})\,u_i=\nabla f(\mathbf{a})\cdot\hat{\mathbf{u}}$. $\square$

The cosine form packs the geometry of §6.4 into one line:

- the maximum value $+\lVert\nabla f\rVert$ occurs for $\hat{\mathbf{u}}$ along $\nabla f$ ($\theta=0$);
- the minimum $-\lVert\nabla f\rVert$ for $\hat{\mathbf{u}}$ opposite to $\nabla f$ ($\theta=\pi$);
- $D_{\hat{\mathbf{u}}}f=0$ for $\hat{\mathbf{u}}$ perpendicular to $\nabla f$, i.e. tangent to the level set ($\theta=\pi/2$);
- reversing the direction reverses the sign: $D_{-\hat{\mathbf{u}}}f=-D_{\hat{\mathbf{u}}}f$.

So the $n$ numbers in the gradient encode the slopes in infinitely many directions; as $\hat{\mathbf{u}}$ rotates once around the circle, $D_{\hat{\mathbf{u}}}f$ traces one period of a sinusoid of amplitude $\lVert\nabla f\rVert$.

```desmos-graph
left=0; right=6.3; bottom=-15; top=15
---
y=4\cos(x)+13\sin(x)|blue
y=\sqrt{185}|red|dashed
y=-\sqrt{185}|red|dashed
y=0|grey
(0.927,12.8)|label:direction (3,4)
(1.273,13.6)|label:gradient direction
```

*Directional derivative of $f=x^2y+y^3$ at $(1,2)$ (gradient $(4,13)$) versus the direction angle $\alpha$: a sinusoid of amplitude $\lVert\nabla f\rVert=\sqrt{185}\approx13.6$, peaking along the gradient and vanishing in the two directions tangent to the level curve.*

> [!example] Worked Example 1 — the basic calculation
> Find the rate of change of $f(x,y)=x^2y+y^3$ at $(1,2)$ in the direction of $\mathbf{v}=(3,4)$.
>
> *Gradient:* $\nabla f=(2xy,\ x^2+3y^2)=(4,\ 13)$ at $(1,2)$.
>
> *Normalise:* $\lVert\mathbf{v}\rVert=5$, so $\hat{\mathbf{u}}=(0.6,\ 0.8)$.
>
> *Dot product:*
> $$
> D_{\hat{\mathbf{u}}}f=4(0.6)+13(0.8)=2.4+10.4=12.8 .
> $$
> This is less than the maximum possible rate $\lVert\nabla f\rVert=\sqrt{185}\approx13.6$, as it must be — $(3,4)$ is only $19.8^\circ$ away from the gradient direction ($\cos\theta=12.8/13.6=0.941$).
>
> *Numerical check:* $[f(1.0006,2.0008)-f(0.9994,1.9992)]/0.002=12.8000$. ✓

> [!warning] Common Pitfall
> Forgetting to normalise the direction vector. Using $\mathbf{v}=(3,4)$ directly would give $\nabla f\cdot\mathbf{v}=64$, five times too large. The quantity $\nabla f\cdot\mathbf{v}$ for a non-unit $\mathbf{v}$ is still meaningful — it is the rate of change *per unit time* for a point moving with velocity $\mathbf{v}$ — but it is not the slope per unit distance.

> [!example] Worked Example 2 — hot spot on a circuit board
> Near a power device, the board temperature is modelled by $T(x,y)=80-x^2-2y^2$ (°C, with $x,y$ in cm, device at the origin). A probe sits at $(3,2)$, where $T=80-9-8=63\ ^\circ\mathrm{C}$.
>
> *Gradient:* $\nabla T=(-2x,\ -4y)=(-6,\ -8)\ ^\circ\mathrm{C/cm}$, with $\lVert\nabla T\rVert=10\ ^\circ\mathrm{C/cm}$.
>
> (a) *Fastest heating:* towards $\nabla T/\lVert\nabla T\rVert=(-0.6,-0.8)$, at $10\ ^\circ\mathrm{C/cm}$. Note that this is **not** straight towards the device at the origin, which would be the direction $(-3,-2)/\sqrt{13}=(-0.832,-0.555)$; the steeper $y$-dependence pulls the gradient towards the $y$-direction.
>
> (b) *Fastest cooling:* the opposite direction, $(0.6,0.8)$, at $-10\ ^\circ\mathrm{C/cm}$.
>
> (c) *No change:* perpendicular to the gradient, $\pm(0.8,-0.6)$ — along the isotherm $x^2+2y^2=17$.
>
> (d) *Heating at exactly $5\ ^\circ\mathrm{C/cm}$:* need $10\cos\theta=5$, $\theta=\pm60^\circ$ from the gradient. The gradient direction has angle $\alpha_0=\operatorname{atan2}(-8,-6)=-126.87^\circ$, so $\alpha=-66.87^\circ$ or $-186.87^\circ$, i.e. $\hat{\mathbf{u}}=(0.393,-0.920)$ or $(-0.993,0.120)$. Check the first: $-6(0.393)-8(-0.920)=-2.36+7.36=5.0$. ✓
>
> (e) *Towards the device:* $D=(-6,-8)\cdot(-0.832,-0.555)=4.99+4.44=9.43\ ^\circ\mathrm{C/cm}$, close to but below the maximum.

> [!example] Worked Example 3 — a three-dimensional potential
> A potential is $V(x,y,z)=xy^2+yz^3$ (volts, coordinates in metres). Find the rate of change of $V$ at $P=(2,-1,1)$ moving towards $Q=(4,1,2)$, and the electric field there.
>
> *Gradient:* $\nabla V=(y^2,\ 2xy+z^3,\ 3yz^2)=(1,\ -3,\ -3)\ \mathrm{V/m}$ at $P$.
>
> *Direction:* $\overrightarrow{PQ}=(2,2,1)$, length $3$, so $\hat{\mathbf{u}}=\tfrac13(2,2,1)$.
>
> *Rate:* $D_{\hat{\mathbf{u}}}V=\tfrac13\left(2-6-3\right)=-\tfrac73\approx-2.33\ \mathrm{V/m}$: the potential falls as the probe moves towards $Q$.
>
> *Field:* $\mathbf{E}=-\nabla V=(-1,\ 3,\ 3)\ \mathrm{V/m}$, magnitude $\sqrt{19}=4.36\ \mathrm{V/m}$. Since $D_{\hat{\mathbf{u}}}V=-\mathbf{E}\cdot\hat{\mathbf{u}}$, the result says the field has a component $+2.33\ \mathrm{V/m}$ along $\overrightarrow{PQ}$, so a positive test charge would be pushed towards $Q$: positive charges move "downhill" in potential.

### Rate of change along a curve

If a point moves along a curve $\mathbf{r}(t)$, the chain rule gives

$$
\frac{d}{dt}f(\mathbf{r}(t))=\nabla f\cdot\mathbf{r}'(t)=\lVert\mathbf{r}'(t)\rVert\;D_{\hat{\mathbf{T}}}f,
$$

where $\hat{\mathbf{T}}=\mathbf{r}'/\lVert\mathbf{r}'\rVert$ is the unit tangent. The rate per unit *time* is the speed times the directional derivative along the path; the rate per unit *arc length* is the directional derivative itself. Integrating $\nabla f\cdot d\mathbf{r}$ along a path gives the total change $f(\text{end})-f(\text{start})$, whatever the path — the fundamental theorem for line integrals of [Ch. 7](../part-2-calculus/07-vector-calculus.md), and the reason the work done moving a charge between two points depends only on the potential difference.

> [!example] Worked Example 4 — a sensor on a circular track
> A sensor moves on the circle $\mathbf{r}(t)=(2\cos t,\ 2\sin t)$ through the temperature field $T=80-x^2-2y^2$ of Worked Example 2. How fast does its reading change, and where is it hottest?
>
> *Directly:* $T(\mathbf{r}(t))=80-4\cos^2t-8\sin^2t=76-4\sin^2t$, so $dT/dt=-8\sin t\cos t=-4\sin2t$.
>
> *By the chain rule:* $\nabla T=(-4\cos t,\ -8\sin t)$, $\mathbf{r}'=(-2\sin t,\ 2\cos t)$, and $\nabla T\cdot\mathbf{r}'=8\sin t\cos t-16\sin t\cos t=-8\sin t\cos t$. ✓
>
> The reading is stationary where $\sin2t=0$: hottest ($76\ ^\circ\mathrm{C}$) at $t=0,\pi$ (the points $(\pm2,0)$) and coolest ($72\ ^\circ\mathrm{C}$) at $t=\pm\pi/2$. At those points $\mathbf{r}'$ is tangent to an isotherm — a first glimpse of the Lagrange condition of §6.11.

### Second directional derivative

Differentiating $g(h)=f(\mathbf{a}+h\hat{\mathbf{u}})$ twice gives the curvature of $f$ along $\hat{\mathbf{u}}$:

$$
D^2_{\hat{\mathbf{u}}}f=g''(0)=\sum_{i,j}f_{x_ix_j}u_iu_j=\hat{\mathbf{u}}^\mathsf{T}\mathbf{H}\,\hat{\mathbf{u}},
$$

with $\mathbf{H}$ the Hessian. For $T=80-x^2-2y^2$, $\mathbf{H}=\operatorname{diag}(-2,-4)$, and along $\hat{\mathbf{u}}=(\cos\alpha,\sin\alpha)$ the curvature is $-2\cos^2\alpha-4\sin^2\alpha$, between $-4$ and $-2$: the surface curves downward in every direction, which is what makes the origin a maximum (§6.10).

> [!warning] Common Pitfall — directional derivatives without differentiability
> The formula $D_{\hat{\mathbf{u}}}f=\nabla f\cdot\hat{\mathbf{u}}$ requires $f$ to be differentiable. For $f=x^2y/(x^4+y^2)$ of §6.1 (with $f(0,0)=0$), every directional derivative at the origin exists: with $\hat{\mathbf{u}}=(a,b)$, $b\ne0$,
> $$
> D_{\hat{\mathbf{u}}}f(\mathbf{0})=\lim_{h\to0}\frac{h^3a^2b}{h\,(h^4a^4+h^2b^2)}=\frac{a^2}{b},
> $$
> and $0$ when $b=0$. But $\nabla f(\mathbf{0})=(0,0)$, which would predict $0$ in every direction. The formula fails because $f$ is not even continuous at the origin. For smooth engineering functions this never bites, but it explains why "all directional derivatives exist" is not the definition of differentiability.

> [!tip] Engineering Note
> In a phased-array antenna or a MIMO receiver, the received power is a function $P(\boldsymbol{\theta})$ of several steering parameters, and tracking loops adjust them along $\nabla P$ (gradient ascent). When only measurements are available, the gradient is estimated by probing two or more directions ("dithering") and solving $D_{\hat{\mathbf{u}}_k}P=\nabla P\cdot\hat{\mathbf{u}}_k$ for $\nabla P$ — two independent directions suffice in two dimensions, precisely because the gradient formula is linear in $\hat{\mathbf{u}}$.

## 6.6 Jacobian matrices

### Motivation

So far $f$ has produced a single number. Many engineering maps produce several: a transistor maps two terminal voltages to two currents; a coordinate change maps $(r,\theta)$ to $(x,y)$; a nonlinear circuit's node equations map a vector of node voltages to a vector of current imbalances. For such vector functions the derivative is no longer a number or a vector but a **matrix** — the Jacobian. It is the object that makes "linearise about an operating point" precise, it is the engine of Newton's method in every circuit simulator, and its determinant is the area-scaling factor that §6.8 needs for changing variables in integrals.

### Definition

A vector function $\mathbf{F}:\mathbb{R}^n\to\mathbb{R}^m$ has $m$ component functions $F_1,\dots,F_m$, each depending on $x_1,\dots,x_n$.

> [!definition] Jacobian matrix
> The Jacobian of $\mathbf{F}$ is the $m\times n$ matrix of first partials,
> $$
> \mathbf{J}_{\mathbf{F}}=\frac{\partial(F_1,\dots,F_m)}{\partial(x_1,\dots,x_n)}=
> \begin{bmatrix}
> \dfrac{\partial F_1}{\partial x_1} & \cdots & \dfrac{\partial F_1}{\partial x_n}\\
> \vdots & \ddots & \vdots\\
> \dfrac{\partial F_m}{\partial x_1} & \cdots & \dfrac{\partial F_m}{\partial x_n}
> \end{bmatrix},\qquad (\mathbf{J}_{\mathbf{F}})_{ij}=\frac{\partial F_i}{\partial x_j}.
> $$
> Row $i$ is $(\nabla F_i)^\mathsf{T}$; column $j$ is the rate of change of the whole output vector with respect to input $x_j$. When $m=n$ the matrix is square and its determinant $\det\mathbf{J}=\dfrac{\partial(F_1,\dots,F_n)}{\partial(x_1,\dots,x_n)}$ is called the **Jacobian determinant** (often just "the Jacobian").

Special cases tie together everything so far. For $m=n=1$ the Jacobian is the ordinary derivative $f'(x)$. For $m=1$ it is the gradient written as a row, $\mathbf{J}_f=(\nabla f)^\mathsf{T}$. For $n=1$ (a curve $\mathbf{r}(t)$) it is the velocity column $\mathbf{r}'(t)$.

### The Jacobian as the best linear approximation

Applying the total differential to each component and stacking the results:

$$
\mathbf{F}(\mathbf{x}_0+d\mathbf{x})\approx\mathbf{F}(\mathbf{x}_0)+\mathbf{J}_{\mathbf{F}}(\mathbf{x}_0)\,d\mathbf{x},\qquad\text{i.e.}\qquad d\mathbf{F}=\mathbf{J}_{\mathbf{F}}\,d\mathbf{x}.
$$

Near $\mathbf{x}_0$, a smooth nonlinear map behaves like the linear map $d\mathbf{x}\mapsto\mathbf{J}\,d\mathbf{x}$ ([Ch. 9](../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md)). Small circles of inputs map to (approximately) small ellipses of outputs; the matrix says how they are stretched, rotated and sheared.

> [!theorem] Chain rule in matrix form
> If $\mathbf{G}:\mathbb{R}^p\to\mathbb{R}^n$ is differentiable at $\mathbf{u}$ and $\mathbf{F}:\mathbb{R}^n\to\mathbb{R}^m$ is differentiable at $\mathbf{G}(\mathbf{u})$, then
> $$
> \mathbf{J}_{\mathbf{F}\circ\mathbf{G}}(\mathbf{u})=\mathbf{J}_{\mathbf{F}}\big(\mathbf{G}(\mathbf{u})\big)\;\mathbf{J}_{\mathbf{G}}(\mathbf{u}),
> $$
> and for square Jacobians, $\det\mathbf{J}_{\mathbf{F}\circ\mathbf{G}}=\det\mathbf{J}_{\mathbf{F}}\cdot\det\mathbf{J}_{\mathbf{G}}$.

*Proof sketch.* Linearise each map: $d\mathbf{G}=\mathbf{J}_{\mathbf{G}}\,d\mathbf{u}$ and $d\mathbf{F}=\mathbf{J}_{\mathbf{F}}\,d\mathbf{G}$; substitute. Entry by entry, $(\mathbf{J}_{\mathbf{F}}\mathbf{J}_{\mathbf{G}})_{ik}=\sum_j\dfrac{\partial F_i}{\partial x_j}\dfrac{\partial x_j}{\partial u_k}$ — exactly the "sum over all paths" of the tree diagram in §6.3. The determinant statement is the product rule for determinants ([Ch. 8](../part-3-linear-algebra/08-matrices-and-linear-systems.md)). $\square$

So the multivariable chain rule, with all its terms, is just matrix multiplication.

> [!example] Worked Example 1 — the polar map
> For $x=r\cos\theta$, $y=r\sin\theta$,
> $$
> \frac{\partial(x,y)}{\partial(r,\theta)}=\begin{bmatrix}x_r & x_\theta\\ y_r & y_\theta\end{bmatrix}=\begin{bmatrix}\cos\theta & -r\sin\theta\\ \sin\theta & r\cos\theta\end{bmatrix},\qquad
> \det\mathbf{J}=r\cos^2\theta+r\sin^2\theta=r .
> $$
> Interpretation: a small change $dr$ moves the point $dr$ metres radially (first column, a unit vector); a small change $d\theta$ moves it $r\,d\theta$ metres tangentially (second column, length $r$). The two columns are perpendicular, so a small $dr\times d\theta$ rectangle maps to a small rectangle of area $r\,dr\,d\theta$. The determinant $r$ is exactly the area factor in $dA=r\,dr\,d\theta$ (§6.8). At the origin $\det\mathbf{J}=0$: the whole line $r=0$ (every $\theta$) collapses to a single point.

> [!example] Worked Example 2 — composing maps: squaring a complex number
> The complex map $w=z^2$, with $z=u+jv$ and $w=p+jq$, is $p=u^2-v^2$, $q=2uv$. Its Jacobian is
> $$
> \mathbf{J}_{\mathbf{F}}=\begin{bmatrix}2u & -2v\\ 2v & 2u\end{bmatrix},\qquad\det\mathbf{J}_{\mathbf{F}}=4(u^2+v^2)=4\lvert z\rvert^2 .
> $$
> (The pattern — equal diagonal entries, opposite off-diagonal ones — is the Cauchy–Riemann structure of an analytic function, and $\det\mathbf{J}=\lvert dw/dz\rvert^2=\lvert2z\rvert^2$; see [Ch. 18](../part-5-transforms/18-complex-analysis.md).)
>
> Now feed in polar coordinates, $u=r\cos\theta$, $v=r\sin\theta$. By the chain rule,
> $$
> \frac{\partial(p,q)}{\partial(r,\theta)}=\begin{bmatrix}2r\cos\theta & -2r\sin\theta\\ 2r\sin\theta & 2r\cos\theta\end{bmatrix}\begin{bmatrix}\cos\theta & -r\sin\theta\\ \sin\theta & r\cos\theta\end{bmatrix}=\begin{bmatrix}2r\cos2\theta & -2r^2\sin2\theta\\ 2r\sin2\theta & 2r^2\cos2\theta\end{bmatrix}.
> $$
> *Check directly:* $w=r^2e^{j2\theta}$, so $p=r^2\cos2\theta$, $q=r^2\sin2\theta$, and differentiating gives the same matrix. ✓ Determinants: $4r^2\cdot r=4r^3$ from the product, and $2r\cdot2r^2(\cos^22\theta+\sin^22\theta)=4r^3$ directly. ✓

### Small-signal models

> [!example] Worked Example 3 — the BJT as a two-port Jacobian
> A simple model of an NPN transistor in the active region is
> $$
> i_C=I_se^{v_{BE}/V_T}\left(1+\frac{v_{CE}}{V_A}\right),\qquad i_B=\frac{I_s}{\beta}e^{v_{BE}/V_T},
> $$
> a map $(v_{BE},v_{CE})\mapsto(i_B,i_C)$. Its Jacobian is
> $$
> \frac{\partial(i_B,i_C)}{\partial(v_{BE},v_{CE})}=\begin{bmatrix}\dfrac{i_B}{V_T} & 0\\[2mm] \dfrac{i_C}{V_T} & \dfrac{I_se^{v_{BE}/V_T}}{V_A}\end{bmatrix}=\begin{bmatrix}1/r_\pi & 0\\ g_m & 1/r_o\end{bmatrix}.
> $$
> *Numbers.* Bias at $I_C=1\ \mathrm{mA}$, $V_{CE}=5\ \mathrm{V}$, with $V_T=25\ \mathrm{mV}$, $V_A=100\ \mathrm{V}$, $\beta=100$. Then $g_m=I_C/V_T=40\ \mathrm{mS}$; $r_\pi\approx\beta/g_m=2.5\ \mathrm{k\Omega}$; and since $I_se^{v_{BE}/V_T}=I_C/(1+V_{CE}/V_A)$, $1/r_o=I_C/(V_A+V_{CE})$, so $r_o=105\ \mathrm{k\Omega}$.
>
> *Using it.* A $1\ \mathrm{mV}$ change in $v_{BE}$ changes $i_C$ by $g_m\cdot1\ \mathrm{mV}=40\ \mu\mathrm{A}$ (linearised). The exact change is $1\ \mathrm{mA}\,(e^{0.04}-1)=40.8\ \mu\mathrm{A}$; the 2 % discrepancy is the curvature of the exponential, and it is the origin of distortion in amplifiers driven too hard. The Jacobian entries are precisely the hybrid-$\pi$ small-signal parameters, and this matrix is the device's **$y$-parameter** matrix at the bias point.

> [!tip] Engineering Note
> Every two-port parameter set ($y$, $z$, $h$, $g$, $ABCD$) of a nonlinear device is a Jacobian at its bias point, expressed with different choices of which port variables are inputs. Converting between parameter sets is inverting or rearranging these matrices, which is why small-signal analysis is linear algebra ([Ch. 8](../part-3-linear-algebra/08-matrices-and-linear-systems.md)).

### The inverse function theorem

> [!theorem] Inverse function theorem
> Let $\mathbf{F}:\mathbb{R}^n\to\mathbb{R}^n$ have continuous partial derivatives near $\mathbf{x}_0$, with $\det\mathbf{J}_{\mathbf{F}}(\mathbf{x}_0)\ne0$. Then $\mathbf{F}$ has a differentiable inverse near $\mathbf{y}_0=\mathbf{F}(\mathbf{x}_0)$, and
> $$
> \mathbf{J}_{\mathbf{F}^{-1}}(\mathbf{y}_0)=\left[\mathbf{J}_{\mathbf{F}}(\mathbf{x}_0)\right]^{-1},\qquad \det\mathbf{J}_{\mathbf{F}^{-1}}=\frac{1}{\det\mathbf{J}_{\mathbf{F}}}.
> $$

The formula follows from the chain rule applied to $\mathbf{F}^{-1}\circ\mathbf{F}=$ identity, whose Jacobian is $\mathbf{I}$. The theorem says that a non-zero Jacobian determinant is what makes a change of variables locally reversible; where it vanishes (the origin in polar coordinates) the coordinates break down.

> [!example] Worked Example 4 — inverse polar derivatives without differentiating
> Inverting the matrix of Worked Example 1 (determinant $r$),
> $$
> \frac{\partial(r,\theta)}{\partial(x,y)}=\frac1r\begin{bmatrix}r\cos\theta & r\sin\theta\\ -\sin\theta & \cos\theta\end{bmatrix}=\begin{bmatrix}\cos\theta & \sin\theta\\ -\dfrac{\sin\theta}{r} & \dfrac{\cos\theta}{r}\end{bmatrix}.
> $$
> *Check directly:* $r=\sqrt{x^2+y^2}$ gives $r_x=x/r=\cos\theta$ and $r_y=y/r=\sin\theta$; $\theta=\operatorname{atan2}(y,x)$ gives $\theta_x=-y/(x^2+y^2)=-\sin\theta/r$ and $\theta_y=x/(x^2+y^2)=\cos\theta/r$. ✓ This resolves the puzzle of §6.3: $x_r=\cos\theta$ and $r_x=\cos\theta$ are not reciprocals; it is the *matrices* that are inverses, not their individual entries.

### Sensitivity of solutions: the implicit function theorem

Suppose a system of $n$ equations $\mathbf{F}(\mathbf{x},\mathbf{p})=\mathbf{0}$ determines $n$ unknowns $\mathbf{x}$ (node voltages, say) in terms of parameters $\mathbf{p}$ (component values). Differentiating the identity $\mathbf{F}(\mathbf{x}(\mathbf{p}),\mathbf{p})=\mathbf{0}$:

$$
\mathbf{J}_{\mathbf{x}}\,d\mathbf{x}+\mathbf{J}_{\mathbf{p}}\,d\mathbf{p}=\mathbf{0}\quad\Longrightarrow\quad \frac{\partial\mathbf{x}}{\partial\mathbf{p}}=-\mathbf{J}_{\mathbf{x}}^{-1}\mathbf{J}_{\mathbf{p}},
$$

valid wherever $\det\mathbf{J}_{\mathbf{x}}\ne0$ (the vector **implicit function theorem**). This is the matrix version of $dy/dx=-F_x/F_y$, and it is how circuit simulators compute the sensitivity of an operating point to every component in one linear solve, reusing the Jacobian they already built for Newton's method.

### Newton's method for nonlinear systems

To solve $\mathbf{F}(\mathbf{x})=\mathbf{0}$, replace $\mathbf{F}$ by its linearisation at the current guess and solve the resulting *linear* system for the step:

$$
\mathbf{J}_{\mathbf{F}}(\mathbf{x}_k)\,\Delta\mathbf{x}=-\mathbf{F}(\mathbf{x}_k),\qquad \mathbf{x}_{k+1}=\mathbf{x}_k+\Delta\mathbf{x}.
$$

Geometrically, each component $F_i$ is replaced by its tangent plane, and the step goes to where all the tangent planes cross zero simultaneously. Near a root with non-singular Jacobian the error is roughly *squared* at each step (**quadratic convergence**), so the number of correct digits about doubles per iteration. This is precisely how SPICE finds the DC operating point of a nonlinear circuit: each diode and transistor is replaced by its small-signal "companion model" (a Jacobian entry plus a source), the resulting linear circuit is solved, and the process repeats.

> [!example] Worked Example 5 — Newton iteration for two equations
> Solve $x^2+y^2=4$, $xy=1$, starting from $(x_0,y_0)=(2,\,0.5)$.
>
> *Set up:* $\mathbf{F}=(x^2+y^2-4,\ xy-1)$, $\mathbf{J}=\begin{bmatrix}2x & 2y\\ y & x\end{bmatrix}$.
>
> *Step 1:* $\mathbf{F}(2,0.5)=(0.25,\ 0)$ and
> $$
> \mathbf{J}=\begin{bmatrix}4 & 1\\ 0.5 & 2\end{bmatrix},\qquad \det\mathbf{J}=7.5,\qquad
> \Delta\mathbf{x}=\frac{1}{7.5}\begin{bmatrix}2 & -1\\ -0.5 & 4\end{bmatrix}\begin{bmatrix}-0.25\\ 0\end{bmatrix}=\begin{bmatrix}-0.0667\\ 0.0167\end{bmatrix}.
> $$
> So $(x_1,y_1)=(1.93333,\ 0.51667)$.
>
> *Further steps* (same recipe):
>
> | $k$ | $x_k$ | $y_k$ | error $\lVert\mathbf{x}_k-\mathbf{x}^*\rVert$ |
> |---|---|---|---|
> | 0 | 2.00000 | 0.50000 | $7.0\times10^{-2}$ |
> | 1 | 1.93333 | 0.51667 | $1.8\times10^{-3}$ |
> | 2 | 1.93185274 | 0.51763705 | $1.5\times10^{-6}$ |
> | 3 | 1.93185165 | 0.51763809 | $1.1\times10^{-12}$ |
>
> The exact root is $x^*=(\sqrt6+\sqrt2)/2=1.93185165$, $y^*=(\sqrt6-\sqrt2)/2=0.51763809$. Each error is roughly a constant times the square of the previous one ($1.8\times10^{-3}\approx0.36\,(7.0\times10^{-2})^2$, $1.5\times10^{-6}\approx0.47\,(1.8\times10^{-3})^2$): quadratic convergence.

```desmos-graph
left=-3; right=3; bottom=-2.5; top=2.5
---
x^2+y^2=4|blue
xy=1|red
(2,0.5)|label:start
(1.9318,0.5176)|label:root
(0.5176,1.9318)
(-1.9318,-0.5176)
(-0.5176,-1.9318)
```

*The circle and hyperbola intersect in four points; Newton's method converges to the one nearest the starting guess. Starting near the line $y=x$, where $\det\mathbf{J}=2(x^2-y^2)$ vanishes, would send the first step far away.*

> [!example] Worked Example 6 — how the solution moves with a parameter
> Replace the circle by $x^2+y^2=c$. How fast does the root found above move as $c$ increases through $4$?
>
> With $\mathbf{F}=(x^2+y^2-c,\ xy-1)$, $\mathbf{J}_{\mathbf{x}}=\begin{bmatrix}2x&2y\\y&x\end{bmatrix}$ and $\mathbf{J}_c=\begin{bmatrix}-1\\0\end{bmatrix}$. The sensitivity formula gives
> $$
> \frac{d\mathbf{x}}{dc}=-\mathbf{J}_{\mathbf{x}}^{-1}\mathbf{J}_c=\frac{1}{2(x^2-y^2)}\begin{bmatrix}x&-2y\\-y&2x\end{bmatrix}\begin{bmatrix}1\\0\end{bmatrix}=\frac{1}{2(x^2-y^2)}\begin{bmatrix}x\\-y\end{bmatrix}.
> $$
> At the root, $x^2-y^2=(x+y)(x-y)=\sqrt6\cdot\sqrt2=3.4641$, so $dx/dc=1.93185/6.9282=0.2788$ and $dy/dc=-0.0747$. Re-solving with $c=4.01$ and differencing gives $(0.2786,\ -0.0745)$. ✓ Enlarging the circle slides the intersection outward along the hyperbola, as the picture suggests.

> [!warning] Common Pitfall
> Newton's method needs a good starting point and a non-singular Jacobian. Near a point where $\det\mathbf{J}\approx0$ the step $\Delta\mathbf{x}=-\mathbf{J}^{-1}\mathbf{F}$ becomes huge and the iteration can jump to a different root or diverge. Circuit simulators protect themselves with damping (limiting the step, especially across exponential junctions), source stepping (ramping supplies up from zero) and "gmin" stepping (adding small conductances that keep $\mathbf{J}$ well conditioned). Never compute $\mathbf{J}^{-1}$ explicitly in practice: solve $\mathbf{J}\Delta\mathbf{x}=-\mathbf{F}$ by elimination ([Ch. 8](../part-3-linear-algebra/08-matrices-and-linear-systems.md)).

### The Hessian as a Jacobian

The gradient $\nabla f$ is itself a vector function $\mathbb{R}^n\to\mathbb{R}^n$, and its Jacobian is the symmetric matrix of second partials:

$$
\mathbf{H}_f=\mathbf{J}_{\nabla f}=\begin{bmatrix}f_{x_1x_1}&\cdots&f_{x_1x_n}\\ \vdots&\ddots&\vdots\\ f_{x_nx_1}&\cdots&f_{x_nx_n}\end{bmatrix},\qquad \mathbf{H}^\mathsf{T}=\mathbf{H}\ \text{(Clairaut)}.
$$

Newton's method applied to $\nabla f=\mathbf{0}$ — the condition for an extremum — reads $\mathbf{H}\,\Delta\mathbf{x}=-\nabla f$, the basis of Newton-type optimisers. The Hessian also classifies critical points (§6.10).

## 6.7 Multiple integration

### Motivation

Integration in one variable adds up a quantity distributed along a line: charge on a wire, energy over time. In two and three dimensions the same idea adds up quantities spread over areas and volumes: the total charge on a plate with non-uniform surface density, the current through a conductor cross-section with non-uniform current density, the energy stored in a field, the probability that two random noise voltages both exceed a threshold. The new ingredient is the *shape* of the region, which determines the limits of integration — and most of the skill in multiple integration lies in describing that shape correctly.

### Double integrals: definition and properties

Let $f$ be bounded on a bounded region $R$ of the plane. Cover $R$ with a grid of small rectangles; let the $k$-th rectangle inside $R$ have area $\Delta A_k$ and contain a sample point $(x_k,y_k)$.

> [!definition] Double integral
> $$
> \iint_R f(x,y)\,dA=\lim_{\max\Delta A_k\to0}\ \sum_k f(x_k,y_k)\,\Delta A_k ,
> $$
> when the limit exists and is independent of the grid and the sample points. It exists whenever $f$ is continuous on a closed, bounded region $R$ with a reasonably smooth boundary.

Interpretations depend on what $f$ is:

| Integrand $f$ | $\iint_R f\,dA$ means | Units |
|---|---|---|
| $1$ | area of $R$ | $\mathrm{m^2}$ |
| height $f\ge0$ | volume under the surface $z=f(x,y)$ | $\mathrm{m^3}$ |
| surface charge density $\rho_s$ | total charge on $R$ | $\mathrm{C}$ |
| current density $J_z$ across a cross-section | current through it | $\mathrm{A}$ |
| joint probability density $p_{XY}$ | probability that $(X,Y)\in R$ | dimensionless |

> [!theorem] Properties of double integrals
> 1. **Linearity:** $\iint_R(af+bg)\,dA=a\iint_R f\,dA+b\iint_R g\,dA$.
> 2. **Additivity:** if $R=R_1\cup R_2$ with $R_1$, $R_2$ overlapping only on boundaries, $\iint_R=\iint_{R_1}+\iint_{R_2}$.
> 3. **Monotonicity:** if $f\le g$ on $R$, then $\iint_R f\,dA\le\iint_R g\,dA$; in particular $m\,A(R)\le\iint_R f\,dA\le M\,A(R)$ when $m\le f\le M$.
> 4. **Mean value theorem:** if $f$ is continuous and $R$ connected, $\iint_R f\,dA=f(\mathbf{c})\,A(R)$ for some point $\mathbf{c}\in R$.

### Iterated integrals and Fubini's theorem

The Riemann-sum definition is not a practical way to compute. Instead, we slice: on a rectangle $[a,b]\times[c,d]$, add up the rectangles in each vertical column first (an integral in $y$ at fixed $x$), then add the columns (an integral in $x$).

> [!theorem] Fubini's theorem
> If $f$ is continuous on the rectangle $[a,b]\times[c,d]$,
> $$
> \iint_R f\,dA=\int_a^b\left[\int_c^d f(x,y)\,dy\right]dx=\int_c^d\left[\int_a^b f(x,y)\,dx\right]dy .
> $$
> More generally, for a **type I** region $R=\{a\le x\le b,\ g_1(x)\le y\le g_2(x)\}$,
> $$
> \iint_R f\,dA=\int_a^b\left[\int_{g_1(x)}^{g_2(x)}f(x,y)\,dy\right]dx ,
> $$
> and for a **type II** region $\{c\le y\le d,\ h_1(y)\le x\le h_2(y)\}$ integrate in $x$ first, between $h_1(y)$ and $h_2(y)$.

*Why it works.* The inner integral $A(x)=\int_{g_1(x)}^{g_2(x)}f(x,y)\,dy$ is the area of the cross-section of the solid under $f$ at position $x$. Integrating cross-sectional area along $x$ gives volume — the slicing method of Chapter 5 — and that volume is the double integral. Slicing the other way must give the same volume.

**Recipe for setting up limits.**

1. Sketch the region and label its boundary curves.
2. Choose the inner variable (say $y$). Draw a line across the region in the $y$-direction at a typical $x$.
3. The inner limits are the $y$-values where the line *enters* and *leaves* the region — in general functions of $x$.
4. The outer limits are the extreme values of $x$ over the whole region — always constants.
5. If the entry or exit curve changes partway across, split the region into pieces.

> [!example] Worked Example 1 — a rectangle, both orders
> Evaluate $\iint_R(x^2+xy)\,dA$ over $R=[0,2]\times[0,1]$.
>
> *$y$ first:*
> $$
> \int_0^2\left[\int_0^1(x^2+xy)\,dy\right]dx=\int_0^2\left[x^2y+\tfrac12xy^2\right]_{y=0}^{1}dx=\int_0^2\left(x^2+\tfrac12x\right)dx=\tfrac83+1=\tfrac{11}{3}.
> $$
> *$x$ first:*
> $$
> \int_0^1\left[\int_0^2(x^2+xy)\,dx\right]dy=\int_0^1\left(\tfrac83+2y\right)dy=\tfrac83+1=\tfrac{11}{3}. \checkmark
> $$
> Numerically $3.6667$. The average value of the integrand over the rectangle is $(11/3)/2=1.833$.

> [!example] Worked Example 2 — area and centroid of a region between curves
> The region between $y=x^2$ and $y=\sqrt{x}$ for $0\le x\le1$ is type I (vertical strips from $y=x^2$ up to $y=\sqrt x$).
>
> *Area:*
> $$
> A=\int_0^1\int_{x^2}^{\sqrt{x}}dy\,dx=\int_0^1\left(\sqrt{x}-x^2\right)dx=\frac23-\frac13=\frac13 .
> $$
> *Centroid:* $\bar x=\dfrac1A\iint x\,dA=3\int_0^1x\left(\sqrt x-x^2\right)dx=3\left(\tfrac25-\tfrac14\right)=\tfrac{9}{20}=0.45$. The region is symmetric about $y=x$ (swapping $x$ and $y$ swaps the two curves), so $\bar y=\bar x=0.45$ without further work; direct integration, $\bar y=3\int_0^1\tfrac12(x-x^4)\,dx=3\cdot\tfrac12\left(\tfrac12-\tfrac15\right)=0.45$, confirms it.
>
> As a type II region the same set is $0\le y\le1$, $y^2\le x\le\sqrt y$ — horizontal strips from the parabola $x=y^2$ to $x=\sqrt y$.

```desmos-graph
left=-0.2; right=1.4; bottom=-0.2; top=1.3
---
y=x^2|blue|0<x<1.2
y=\sqrt{x}|red|0<x<1.2
x^2<y<\sqrt{x}|green|0<x<1
x=0.6|grey|dashed
(0.45,0.45)|label:centroid
```

*The shaded region is type I (vertical strips such as the dashed one, entering at $y=x^2$ and leaving at $y=\sqrt{x}$) and also type II (horizontal strips from $x=y^2$ to $x=\sqrt{y}$). Its centroid lies on the symmetry line $y=x$.*

> [!example] Worked Example 3 — charge on a triangular plate
> A triangular plate occupies $0\le y\le x\le a$ and carries surface charge density $\rho_s=\rho_0\,xy/a^2$ ($\rho_0$ in $\mathrm{C/m^2}$). Find the total charge.
>
> *Region:* vertical strips at fixed $x$ run from $y=0$ to $y=x$, for $0\le x\le a$.
> $$
> Q=\frac{\rho_0}{a^2}\int_0^a\int_0^x xy\,dy\,dx=\frac{\rho_0}{a^2}\int_0^a x\cdot\frac{x^2}{2}\,dx=\frac{\rho_0}{a^2}\cdot\frac{a^4}{8}=\frac{\rho_0a^2}{8}.
> $$
> *Units:* $\mathrm{C/m^2}\times\mathrm{m^2}=\mathrm{C}$. ✓ *Sanity check:* the plate area is $a^2/2$ and the density ranges from $0$ to $\rho_0$, so $Q$ must lie between $0$ and $\rho_0a^2/2$; the answer is a quarter of the upper bound, reasonable since the density is small over most of the plate.

### Reversing the order of integration

Sometimes the inner integral cannot be done in the given order, but the other order is easy. Reversing is a matter of re-describing the *same region* the other way round — never of swapping the limits mechanically.

> [!example] Worked Example 4 — reversing to remove a non-elementary integral
> (a) $I=\displaystyle\int_0^1\int_y^1 e^{x^2}\,dx\,dy$ cannot be done as written, since $e^{x^2}$ has no elementary antiderivative.
>
> *Region:* $0\le y\le1$ and $y\le x\le1$, i.e. the triangle $0\le y\le x\le1$ with vertices $(0,0)$, $(1,0)$, $(1,1)$.
>
> *Re-describe as type I:* for each $x\in[0,1]$, $y$ runs from $0$ to $x$. Then
> $$
> I=\int_0^1\int_0^x e^{x^2}\,dy\,dx=\int_0^1 x\,e^{x^2}\,dx=\left[\tfrac12e^{x^2}\right]_0^1=\frac{e-1}{2}\approx0.859 .
> $$
> (b) The same trick handles $J=\displaystyle\int_0^\pi\int_y^\pi\frac{\sin x}{x}\,dx\,dy$, whose inner integral is the non-elementary sine integral. The region is $0\le y\le x\le\pi$, so
> $$
> J=\int_0^\pi\int_0^x\frac{\sin x}{x}\,dy\,dx=\int_0^\pi\frac{\sin x}{x}\cdot x\,dx=\int_0^\pi\sin x\,dx=2 .
> $$
> Numerical double integration in the original order gives $2.0000$. ✓ The pattern: the inner integration in the new order produces a factor ($x$) that cancels the troublesome part.

> [!warning] Common Pitfall
> The outer limits must be **constants** and the inner limits may depend only on the **outer** variable. A setup like $\int_0^1\int_0^x f\,dx\,dy$ (inner limit depending on the inner variable) is meaningless. When reversing the order, always sketch the region: swapping $\int_0^1\int_y^1$ into $\int_y^1\int_0^1$ is a common and fatal error.

### Triple integrals

For $f(x,y,z)$ on a solid $V$, the triple integral $\iiint_V f\,dV$ is defined by Riemann sums over small boxes and evaluated, by Fubini's theorem again, as three nested single integrals. A typical setup is

$$
\iiint_V f\,dV=\int_a^b\int_{g_1(x)}^{g_2(x)}\int_{h_1(x,y)}^{h_2(x,y)}f(x,y,z)\,dz\,dy\,dx ,
$$

where for each $(x,y)$ in the "shadow" of $V$ on the $xy$-plane, $z$ runs from the bottom surface $h_1$ to the top surface $h_2$, and the shadow itself is described as a 2-D region. With $f=1$ the integral gives volume; with $f=\rho_v$ (volume charge density, $\mathrm{C/m^3}$) it gives charge; with $f=w_e=\tfrac12\varepsilon\lVert\mathbf{E}\rVert^2$ it gives stored electric energy.

> [!example] Worked Example 5 — separable integrand on a cube
> A cube $0\le x,y,z\le a$ carries $\rho_v=\rho_0\,xyz/a^3$. Both the region and the integrand separate, so the triple integral is a product of single integrals:
> $$
> Q=\frac{\rho_0}{a^3}\int_0^a x\,dx\int_0^a y\,dy\int_0^a z\,dz=\frac{\rho_0}{a^3}\left(\frac{a^2}{2}\right)^3=\frac{\rho_0a^3}{8}.
> $$
> Units: $\mathrm{C/m^3}\times\mathrm{m^3}=\mathrm{C}$. ✓ Because the density is a product of three linear ramps, its average over the cube is $\rho_0/8$, and $Q=(\rho_0/8)\cdot a^3$.

> [!example] Worked Example 6 — the tetrahedron
> Let $V$ be the tetrahedron bounded by the coordinate planes and $x+y+z=1$. Find its volume and the height $\bar z$ of its centroid.
>
> *Limits:* top surface $z=1-x-y$, bottom $z=0$; the shadow is the triangle $x\ge0$, $y\ge0$, $x+y\le1$, i.e. $0\le x\le1$, $0\le y\le1-x$.
>
> *Volume:*
> $$
> \int_0^1\int_0^{1-x}(1-x-y)\,dy\,dx=\int_0^1\frac{(1-x)^2}{2}\,dx=\frac16 .
> $$
> *Moment:*
> $$
> \iiint_V z\,dV=\int_0^1\int_0^{1-x}\frac{(1-x-y)^2}{2}\,dy\,dx=\int_0^1\frac{(1-x)^3}{6}\,dx=\frac1{24}.
> $$
> Hence $\bar z=(1/24)/(1/6)=1/4$, and by symmetry the centroid is $\left(\tfrac14,\tfrac14,\tfrac14\right)$ — one quarter of the way from the base face to the opposite vertex, the standard result for any tetrahedron. Numerical triple integration gives $0.041667=1/24$. ✓

### Averages, centroids and moments

$$
\bar f=\frac{1}{A}\iint_R f\,dA,\qquad \bar x=\frac{1}{A}\iint_R x\,dA,\qquad \bar y=\frac{1}{A}\iint_R y\,dA,
$$

with the obvious extensions to volumes, and to *weighted* centres when a density $\rho$ multiplies the integrand and the total "mass" $\iint\rho\,dA$ replaces $A$. **Second moments** $\iint x^2\rho\,dA$ appear as moments of inertia in mechanics and as variances in probability ([Ch. 19](../part-6-probability-statistics/19-probability.md)); the root-mean-square value of a field over an area, $\sqrt{\frac1A\iint f^2\,dA}$, is the quantity that determines average power.

### Numerical double integration

When no antiderivative exists, the iterated form still allows numerical evaluation: apply any one-dimensional rule ([Ch. 22](../part-7-numerical-computational/22-numerical-methods.md)) to the inner integral and again to the outer one. The two-dimensional **midpoint rule** on an $n\times n$ grid of cells of size $h\times k$ is

$$
\iint_R f\,dA\approx hk\sum_{i=1}^{n}\sum_{j=1}^{n}f(x_i^\ast,y_j^\ast),
$$

with $(x_i^\ast,y_j^\ast)$ the cell centres.

> [!example] Worked Example 7 — midpoint rule for a Gaussian on the unit square
> $I=\iint_{[0,1]^2}e^{-(x^2+y^2)}\,dA=\left(\int_0^1e^{-x^2}dx\right)^2=0.746824^2=0.557746$ (separable, using the tabulated 1-D value).
>
> | grid | midpoint estimate | error |
> |---|---|---|
> | $2\times2$ | 0.56942 | $1.17\times10^{-2}$ |
> | $4\times4$ | 0.56062 | $2.88\times10^{-3}$ |
> | $8\times8$ | 0.55846 | $7.16\times10^{-4}$ |
>
> Halving the cell size divides the error by about 4: the midpoint rule is second-order accurate in 2-D, as in 1-D. But the cost — the number of function evaluations — grows as $n^2$ in 2-D and $n^d$ in $d$ dimensions. For high-dimensional integrals (bit-error rates over many noise variables, yield over many component tolerances) this "curse of dimensionality" is why Monte Carlo methods take over ([Ch. 24](../part-7-numerical-computational/24-computational-mathematics.md)).

> [!tip] Engineering Note
> Improper double integrals — over unbounded regions or with singular integrands — are routine in EE: the probability that a Gaussian noise vector lands outside a decision region, the energy of a field extending to infinity, the potential of a charged sheet. They are defined as limits of integrals over growing (or shrinking-hole) regions, and for non-negative integrands Fubini's theorem still allows any order. The most important one, $\iint_{\mathbb{R}^2}e^{-(x^2+y^2)}dA=\pi$, is evaluated in §6.8.

## 6.8 Coordinate transformations

### Motivation

Cartesian coordinates make rectangles easy and everything else awkward. A disc, a sector, a region bounded by hyperbolas, an integrand containing $x^2+y^2$ — each becomes simple in the right coordinates. Changing variables in a multiple integral is the multivariable analogue of substitution, $\int f(x)\,dx=\int f(x(u))\,x'(u)\,du$, and the role of $x'(u)$ is played by the Jacobian determinant of §6.6. The idea also has a probabilistic face: transforming random variables (from in-phase and quadrature noise to envelope and phase, say) uses exactly the same Jacobian factor.

### How a map scales area

Let $x=x(u,v)$, $y=y(u,v)$ be a one-to-one, continuously differentiable map from a region $S$ of the $uv$-plane onto a region $R$ of the $xy$-plane, and write $\mathbf{r}(u,v)=(x(u,v),y(u,v))$. A small coordinate rectangle with corners $(u,v)$, $(u+du,v)$, $(u,v+dv)$ is carried to an approximate **parallelogram** with sides

$$
\mathbf{r}(u+du,v)-\mathbf{r}(u,v)\approx\mathbf{r}_u\,du=(x_u,\,y_u)\,du,\qquad
\mathbf{r}(u,v+dv)-\mathbf{r}(u,v)\approx\mathbf{r}_v\,dv=(x_v,\,y_v)\,dv .
$$

The area of a parallelogram spanned by $(a_1,a_2)$ and $(b_1,b_2)$ is $\lvert a_1b_2-a_2b_1\rvert$, the absolute value of a $2\times2$ determinant. Hence

$$
dA=\left\lvert x_uy_v-x_vy_u\right\rvert du\,dv=\left\lvert\det\frac{\partial(x,y)}{\partial(u,v)}\right\rvert du\,dv .
$$

The Jacobian determinant is the local **area magnification factor** of the map. Summing over all the little parallelograms and taking the limit gives the theorem.

> [!theorem] Change of variables in multiple integrals
> Under the conditions above (one-to-one, continuously differentiable, Jacobian determinant non-zero except possibly on a set of zero area),
> $$
> \iint_R f(x,y)\,dx\,dy=\iint_S f\big(x(u,v),y(u,v)\big)\left\lvert\frac{\partial(x,y)}{\partial(u,v)}\right\rvert du\,dv .
> $$
> In three dimensions, with $(x,y,z)=\mathbf{r}(u,v,w)$, the small box becomes a parallelepiped of volume $\lvert\det\mathbf{J}\rvert\,du\,dv\,dw$, and
> $$
> \iiint_V f\,dx\,dy\,dz=\iiint_W f\big(\mathbf{r}(u,v,w)\big)\left\lvert\frac{\partial(x,y,z)}{\partial(u,v,w)}\right\rvert du\,dv\,dw .
> $$

Three practical facts:

- **The absolute value is essential.** Area and volume are positive; a negative Jacobian determinant only means the map reverses orientation (like a mirror).
- **Inverse Jacobians are reciprocals** (inverse function theorem, §6.6): $\dfrac{\partial(x,y)}{\partial(u,v)}=1\Big/\dfrac{\partial(u,v)}{\partial(x,y)}$. When the substitution is given as $u=u(x,y)$, $v=v(x,y)$, compute the easy Jacobian and invert it — there is no need to solve for $x$ and $y$.
- **Choose the substitution to simplify the region first.** The integrand usually follows. Good choices make the boundary curves of $R$ into coordinate lines $u=\text{const}$, $v=\text{const}$, so that $S$ is a rectangle.

> [!example] Worked Example 1 — a linear change of variables
> Evaluate $\iint_D (x+y)\,dA$ over the parallelogram $D$ bounded by $x+y=0$, $x+y=2$, $x-y=0$, $x-y=1$.
>
> *Substitution suggested by the boundaries:* $u=x+y\in[0,2]$, $v=x-y\in[0,1]$; $S$ is a rectangle.
>
> *Jacobian:* $\dfrac{\partial(u,v)}{\partial(x,y)}=\det\begin{bmatrix}1&1\\1&-1\end{bmatrix}=-2$, so $\left\lvert\dfrac{\partial(x,y)}{\partial(u,v)}\right\rvert=\dfrac12$.
>
> *Integral:*
> $$
> \iint_D(x+y)\,dA=\int_0^1\int_0^2 u\cdot\tfrac12\,du\,dv=\tfrac12\cdot2\cdot1=1 .
> $$
> *Check:* the area of $D$ is $\iint_S\tfrac12\,du\,dv=\tfrac12\cdot2\cdot1=1$, and the average of $u=x+y$ over $D$ is clearly $1$ (it runs uniformly from $0$ to $2$), so the integral is $1\times1=1$. ✓

> [!example] Worked Example 2 — the ellipse
> Find the area of the ellipse $x^2/a^2+y^2/b^2\le1$, and $\iint x^2\,dA$ over it.
>
> *Substitution:* $x=au$, $y=bv$ maps the unit disc $u^2+v^2\le1$ onto the ellipse, with $\partial(x,y)/\partial(u,v)=ab$.
>
> *Area:* $\iint_{\text{disc}}ab\,du\,dv=ab\cdot\pi=\pi ab$.
>
> *Second moment:* $\iint x^2\,dA=\iint_{\text{disc}}a^2u^2\cdot ab\,du\,dv=a^3b\iint_{\text{disc}}u^2\,du\,dv=a^3b\cdot\dfrac{\pi}{4}$, using the polar result $\iint_{\text{disc}}u^2=\int_0^{2\pi}\!\int_0^1r^2\cos^2\theta\,r\,dr\,d\theta=\tfrac14\cdot\pi$. For $a=2$, $b=1$: $2\pi=6.2832$, matching direct numerical integration. ✓ The linear stretch turned an ellipse problem into a circle problem.

> [!example] Worked Example 3 — a region bounded by hyperbolas and lines
> Find the area of the region $R$ in the first quadrant bounded by $xy=1$, $xy=4$, $y=x$ and $y=2x$.
>
> *Substitution:* $u=xy\in[1,4]$, $v=y/x\in[1,2]$ — again the boundaries become coordinate lines.
>
> *Jacobian of the easy direction:*
> $$
> \frac{\partial(u,v)}{\partial(x,y)}=\det\begin{bmatrix}y & x\\ -y/x^2 & 1/x\end{bmatrix}=\frac yx+\frac yx=2v,
> \qquad\text{so}\qquad \left\lvert\frac{\partial(x,y)}{\partial(u,v)}\right\rvert=\frac{1}{2v}.
> $$
> *Area:*
> $$
> A=\int_1^4\int_1^2\frac{1}{2v}\,dv\,du=3\cdot\tfrac12\ln2=\tfrac32\ln2\approx1.0397 .
> $$
> A direct Cartesian computation would need the region split into three pieces (the entry and exit curves change at $x=1/\sqrt2$, $1$, $\sqrt2$ and $2$); numerically it gives $1.0397$. ✓

```desmos-graph
left=0; right=3; bottom=0; top=3.5
---
xy=1|blue
xy=4|blue
y=x|red
y=2x|red
1<xy<4|green
```

*The region of Worked Example 3 is bounded by two hyperbolas (blue) and two lines (red); in the coordinates $u=xy$, $v=y/x$ it becomes the rectangle $[1,4]\times[1,2]$ — only the part between the red lines counts.*

### Polar coordinates

With $x=r\cos\theta$, $y=r\sin\theta$, the Jacobian determinant is $r$ (§6.6), so

$$
dA=r\,dr\,d\theta .
$$

The factor $r$ has a direct geometric meaning: the small polar "rectangle" between radii $r$ and $r+dr$ and angles $\theta$ and $\theta+d\theta$ has radial side $dr$ and arc side $r\,d\theta$. Use polar coordinates whenever the region is a disc, annulus, sector or is bounded by a polar curve $r=g(\theta)$, or when the integrand contains $x^2+y^2$.

> [!example] Worked Example 4 — area inside a cardioid
> The cardioid $r=1+\cos\theta$ ($0\le\theta<2\pi$) — the shape of a cardioid microphone's polar pattern — encloses the region $0\le r\le1+\cos\theta$. Its area is
> $$
> A=\int_0^{2\pi}\int_0^{1+\cos\theta}r\,dr\,d\theta=\int_0^{2\pi}\frac{(1+\cos\theta)^2}{2}\,d\theta=\frac12\int_0^{2\pi}\left(1+2\cos\theta+\cos^2\theta\right)d\theta=\frac12\left(2\pi+0+\pi\right)=\frac{3\pi}{2}.
> $$
> Numerically $4.712$, which sits sensibly between the area of the circle of radius 1 ($\pi$) and the circle of radius 2 that just encloses the cardioid ($4\pi$).

```desmos-graph
left=-1; right=2.5; bottom=-1.6; top=1.6
---
r=1+\cos(\theta)|blue
r=1|grey|dashed
```

*The cardioid $r=1+\cos\theta$ (blue) compared with the unit circle; integrating $r\,dr$ from $0$ out to the curve at each angle sweeps out its area $3\pi/2$.*

> [!example] Worked Example 5 — the Gaussian integral
> Let $I=\int_{-\infty}^{\infty}e^{-x^2}dx$. No antiderivative exists, but the *square* of $I$ is a double integral that polar coordinates make trivial:
> $$
> I^2=\int_{-\infty}^{\infty}e^{-x^2}dx\int_{-\infty}^{\infty}e^{-y^2}dy=\iint_{\mathbb{R}^2}e^{-(x^2+y^2)}\,dx\,dy=\int_0^{2\pi}\!\!\int_0^\infty e^{-r^2}\,r\,dr\,d\theta .
> $$
> The factor $r$ from the Jacobian is exactly what is needed: with $w=r^2$, $\int_0^\infty e^{-r^2}r\,dr=\tfrac12\int_0^\infty e^{-w}dw=\tfrac12$. Hence $I^2=2\pi\cdot\tfrac12=\pi$ and
> $$
> \int_{-\infty}^{\infty}e^{-x^2}dx=\sqrt{\pi},\qquad\int_{-\infty}^{\infty}e^{-x^2/(2\sigma^2)}dx=\sigma\sqrt{2\pi}.
> $$
> The second form (substitute $x\to x/(\sigma\sqrt2)$) normalises the Gaussian (normal) distribution of thermal noise ([Ch. 19](../part-6-probability-statistics/19-probability.md)) and underlies the Fourier transform of a Gaussian pulse ([Ch. 16](../part-5-transforms/16-fourier-transforms.md)).

> [!example] Worked Example 6 — the Rayleigh distribution of a noise envelope
> The in-phase and quadrature components $X$, $Y$ of narrowband noise are independent Gaussians with zero mean and variance $\sigma^2$, with joint density
> $$
> p(x,y)=\frac{1}{2\pi\sigma^2}e^{-(x^2+y^2)/(2\sigma^2)} .
> $$
> What is the probability that the envelope $R=\sqrt{X^2+Y^2}$ exceeds a threshold $\rho$?
>
> The event $R\le\rho$ is the disc of radius $\rho$, so in polar coordinates
> $$
> P(R\le\rho)=\int_0^{2\pi}\!\!\int_0^\rho\frac{1}{2\pi\sigma^2}e^{-r^2/(2\sigma^2)}\,r\,dr\,d\theta=\left[-e^{-r^2/(2\sigma^2)}\right]_0^\rho=1-e^{-\rho^2/(2\sigma^2)} .
> $$
> Differentiating in $\rho$ gives the **Rayleigh density** $p_R(\rho)=\dfrac{\rho}{\sigma^2}e^{-\rho^2/(2\sigma^2)}$ — and the leading factor $\rho$ is nothing but the polar Jacobian. For a threshold $\rho=3\sigma$, $P(R>3\sigma)=e^{-4.5}=0.0111$: about 1.1 % of the time, noise alone crosses a $3\sigma$ envelope threshold. This is the basic false-alarm calculation of a radar or envelope detector ([Ch. 21](../part-6-probability-statistics/21-random-processes.md)).

> [!warning] Common Pitfall
> Forgetting the Jacobian. "Area of the unit disc $=\int_0^{2\pi}\int_0^1dr\,d\theta=2\pi$" is wrong by a factor of 2; the correct $\int_0^{2\pi}\int_0^1r\,dr\,d\theta=\pi$. A quick check that catches this every time: integrate $f=1$ and see whether you get the right area or volume. A second trap is a map that is not one-to-one: letting $\theta$ run over $[0,4\pi]$ counts every point twice.

> [!tip] Engineering Note
> The same Jacobian factor transforms probability densities. If random variables $(X,Y)$ with joint density $p_{XY}$ are mapped one-to-one to $(U,V)$, then $p_{UV}(u,v)=p_{XY}(x,y)\left\lvert\partial(x,y)/\partial(u,v)\right\rvert$, because probability (like area times density) must be conserved. Worked Example 6 is the special case $(X,Y)\to(R,\Theta)$; it shows at once that the phase $\Theta$ is uniform on $[0,2\pi)$ and independent of the envelope.

## 6.9 Cylindrical and spherical coordinates

### Motivation

Wires, coaxial cables, solenoids and circular waveguides have an axis of symmetry; point charges, antennas seen from far away and spherical shells have a centre. Cartesian coordinates fight these geometries: the surface of a wire becomes $x^2+y^2=a^2$, and its cross-section a region with curved limits. In **cylindrical** and **spherical** coordinates the same surfaces become $r=a$ and $R=a$, the limits of integration become constants, and the integrals often separate into products of single integrals. This section sets up both systems, derives their volume elements from the Jacobian, and introduces the scale factors that Chapter 7 uses to write gradient, divergence and curl in curvilinear form.

### Definitions and conventions

To avoid clashing with charge density $\rho$, this book uses $(r,\phi,z)$ for cylindrical and $(R,\theta,\phi)$ for spherical coordinates, as in Cheng's field-theory text.

> [!definition] Cylindrical coordinates $(r,\phi,z)$
> $r\ge0$ is the perpendicular distance from the $z$-axis, $\phi\in[0,2\pi)$ the azimuth measured from the $+x$ axis towards $+y$, and $z$ the usual height:
> $$
> x=r\cos\phi,\qquad y=r\sin\phi,\qquad z=z;\qquad r=\sqrt{x^2+y^2},\quad \phi=\operatorname{atan2}(y,x).
> $$
> Coordinate surfaces: $r=$ const is a circular cylinder, $\phi=$ const a half-plane hinged on the $z$-axis, $z=$ const a horizontal plane.

> [!definition] Spherical coordinates $(R,\theta,\phi)$
> $R\ge0$ is the distance from the origin, $\theta\in[0,\pi]$ the **polar angle** measured from the $+z$ axis, and $\phi\in[0,2\pi)$ the same azimuth as in cylindrical coordinates:
> $$
> x=R\sin\theta\cos\phi,\qquad y=R\sin\theta\sin\phi,\qquad z=R\cos\theta;
> $$
> $$
> R=\sqrt{x^2+y^2+z^2},\qquad \theta=\arccos\frac{z}{R},\qquad \phi=\operatorname{atan2}(y,x).
> $$
> Coordinate surfaces: $R=$ const is a sphere, $\theta=$ const a cone about the $z$-axis (the plane $z=0$ when $\theta=\pi/2$), $\phi=$ const a half-plane.

The two systems are linked by $r=R\sin\theta$ and $z=R\cos\theta$: the cylindrical radius is the projection of the spherical radius onto the $xy$-plane. On the Earth, $\theta$ is the colatitude ($90^\circ$ minus latitude) and $\phi$ the longitude.

> [!example] Worked Example 1 — converting a point
> Express $P=(1,\,1,\,\sqrt2)$ in cylindrical and spherical coordinates.
>
> *Cylindrical:* $r=\sqrt{1+1}=\sqrt2$, $\phi=\operatorname{atan2}(1,1)=\pi/4$, $z=\sqrt2$; so $P=(\sqrt2,\ \pi/4,\ \sqrt2)$.
>
> *Spherical:* $R=\sqrt{1+1+2}=2$, $\theta=\arccos(\sqrt2/2)=\pi/4$, $\phi=\pi/4$; so $P=(2,\ \pi/4,\ \pi/4)$.
>
> *Check back:* $x=2\sin\frac\pi4\cos\frac\pi4=2\cdot\frac12=1$, $y=1$, $z=2\cos\frac\pi4=\sqrt2$. ✓ And $r=R\sin\theta=2\cdot\frac{\sqrt2}{2}=\sqrt2$. ✓
>
> For the point $(-1,-1,0)$, $\tan\phi=y/x=1$ would suggest $\phi=\pi/4$, but the point is in the third quadrant: $\operatorname{atan2}(-1,-1)=-3\pi/4$, i.e. $\phi=5\pi/4$.

> [!warning] Common Pitfall
> Mathematics texts (and some software) swap the roles of $\theta$ and $\phi$ in spherical coordinates, using $\theta$ for the azimuth and $\phi$ for the polar angle. Always check which angle is measured from the $z$-axis before using a formula from another source — the volume element $R^2\sin\theta$ contains the sine of the *polar* angle. Also, $\tan\phi=y/x$ alone does not fix the quadrant; use the signs of $x$ and $y$ (the `atan2` function).

### Jacobians and volume elements

The cylindrical Jacobian is the polar one with a trivial $z$ row and column added, so its determinant is $r$. For spherical coordinates,

$$
\frac{\partial(x,y,z)}{\partial(R,\theta,\phi)}=\begin{bmatrix}
\sin\theta\cos\phi & R\cos\theta\cos\phi & -R\sin\theta\sin\phi\\
\sin\theta\sin\phi & R\cos\theta\sin\phi & R\sin\theta\cos\phi\\
\cos\theta & -R\sin\theta & 0
\end{bmatrix}.
$$

Expanding along the bottom row,

$$
\begin{aligned}
\det\mathbf{J}&=\cos\theta\left(R^2\sin\theta\cos\theta\cos^2\phi+R^2\sin\theta\cos\theta\sin^2\phi\right)+R\sin\theta\left(R\sin^2\theta\cos^2\phi+R\sin^2\theta\sin^2\phi\right)\\
&=R^2\sin\theta\cos^2\theta+R^2\sin^3\theta=R^2\sin\theta .
\end{aligned}
$$

Since $0\le\theta\le\pi$, $\sin\theta\ge0$ and no absolute value is needed. The determinant vanishes on the $z$-axis ($\theta=0,\pi$) and at the origin, where the coordinates are degenerate (every $\phi$ describes the same point).

### Scale factors and the geometric picture

The Jacobian can be read off without any determinant, from geometry. Each column of the Jacobian is the displacement per unit change of one coordinate, and in both systems the three columns are **mutually perpendicular** (the systems are *orthogonal*). Their lengths are the **scale factors** $h_i$: the distance moved per unit change of coordinate $i$.

- Changing $r$ by $dr$ moves you $dr$; changing $z$ by $dz$ moves you $dz$; changing $\phi$ by $d\phi$ moves you along a circle of radius $r$, a distance $r\,d\phi$.
- In spherical coordinates, changing $R$ by $dR$ moves you $dR$; changing $\theta$ by $d\theta$ moves you along a great circle (a meridian) of radius $R$, a distance $R\,d\theta$; changing $\phi$ by $d\phi$ moves you along a circle of latitude, whose radius is $r=R\sin\theta$, a distance $R\sin\theta\,d\phi$.

A small coordinate box is therefore (to first order) a rectangular box with edges $h_1du_1$, $h_2du_2$, $h_3du_3$, and

$$
dV=h_1h_2h_3\,du_1\,du_2\,du_3,\qquad d\ell^2=h_1^2du_1^2+h_2^2du_2^2+h_3^2du_3^2 .
$$

| | Cylindrical $(r,\phi,z)$ | Spherical $(R,\theta,\phi)$ |
|---|---|---|
| to Cartesian | $x=r\cos\phi,\ y=r\sin\phi,\ z=z$ | $x=R\sin\theta\cos\phi,\ y=R\sin\theta\sin\phi,\ z=R\cos\theta$ |
| from Cartesian | $r=\sqrt{x^2+y^2},\ \phi=\operatorname{atan2}(y,x)$ | $R=\sqrt{x^2+y^2+z^2},\ \cos\theta=z/R$ |
| ranges | $r\ge0,\ 0\le\phi<2\pi$ | $R\ge0,\ 0\le\theta\le\pi,\ 0\le\phi<2\pi$ |
| scale factors $(h_1,h_2,h_3)$ | $(1,\ r,\ 1)$ | $(1,\ R,\ R\sin\theta)$ |
| Jacobian determinant | $r$ | $R^2\sin\theta$ |
| volume element $dV$ | $r\,dr\,d\phi\,dz$ | $R^2\sin\theta\,dR\,d\theta\,d\phi$ |
| line element $d\ell^2$ | $dr^2+r^2d\phi^2+dz^2$ | $dR^2+R^2d\theta^2+R^2\sin^2\theta\,d\phi^2$ |
| surface element on $r=a$ / $R=a$ | $a\,d\phi\,dz$ | $a^2\sin\theta\,d\theta\,d\phi$ |
| surface element on $z=$ const / $\theta=$ const | $r\,dr\,d\phi$ | $R\sin\theta\,dR\,d\phi$ |

The scale factors also generate the curvilinear forms of grad, div and curl in [Ch. 7](../part-2-calculus/07-vector-calculus.md); for instance the gradient in spherical coordinates is $\nabla f=f_R\hat{\mathbf{R}}+\dfrac1Rf_\theta\hat{\boldsymbol{\theta}}+\dfrac1{R\sin\theta}f_\phi\hat{\boldsymbol{\phi}}$, each partial divided by its scale factor, exactly as in the polar case of §6.4.

### Unit vectors depend on position

The unit vectors of curvilinear coordinates point along the directions of increasing coordinates, and those directions change from point to point:

$$
\hat{\mathbf{r}}=\cos\phi\,\hat{\mathbf{x}}+\sin\phi\,\hat{\mathbf{y}},\qquad \hat{\boldsymbol{\phi}}=-\sin\phi\,\hat{\mathbf{x}}+\cos\phi\,\hat{\mathbf{y}},\qquad \hat{\mathbf{z}}=\hat{\mathbf{z}};
$$

$$
\hat{\mathbf{R}}=\sin\theta\cos\phi\,\hat{\mathbf{x}}+\sin\theta\sin\phi\,\hat{\mathbf{y}}+\cos\theta\,\hat{\mathbf{z}},\quad
\hat{\boldsymbol{\theta}}=\cos\theta\cos\phi\,\hat{\mathbf{x}}+\cos\theta\sin\phi\,\hat{\mathbf{y}}-\sin\theta\,\hat{\mathbf{z}},\quad
\hat{\boldsymbol{\phi}}=-\sin\phi\,\hat{\mathbf{x}}+\cos\phi\,\hat{\mathbf{y}} .
$$

(They are the Jacobian columns divided by the scale factors.) Each triple is orthonormal and right-handed at every point, but $\hat{\mathbf{r}}$ at $\phi=0$ is $\hat{\mathbf{x}}$ while at $\phi=\pi$ it is $-\hat{\mathbf{x}}$.

> [!warning] Common Pitfall — integrating vectors in curvilinear coordinates
> Because $\hat{\mathbf{r}}$, $\hat{\boldsymbol{\phi}}$, $\hat{\mathbf{R}}$, $\hat{\boldsymbol{\theta}}$ vary with position, they cannot be taken outside an integral. The total field of a uniformly charged ring at its centre is $\oint\hat{\mathbf{r}}\,d\phi\cdots=\mathbf{0}$, not "$2\pi\hat{\mathbf{r}}\cdots$". To integrate a vector quantity, convert the unit vectors to the constant Cartesian $\hat{\mathbf{x}},\hat{\mathbf{y}},\hat{\mathbf{z}}$ first (using the formulas above), then integrate component by component. Scalar integrands — charge, energy, power — have no such problem.

### Worked examples in cylindrical coordinates

> [!example] Worked Example 2 — current in a wire with non-uniform density
> At high frequency, current crowds towards a conductor's surface (the skin effect). Model it in a wire of radius $a$ by $J_z=J_0\,(r/a)^2$. The cross-section is $0\le r\le a$, $0\le\phi<2\pi$, with $dA=r\,dr\,d\phi$, so
> $$
> I=\int_0^{2\pi}\!\!\int_0^a J_0\frac{r^2}{a^2}\,r\,dr\,d\phi=2\pi\frac{J_0}{a^2}\cdot\frac{a^4}{4}=\frac{\pi J_0a^2}{2},
> $$
> half what a uniform density $J_0$ would carry. Units: $\mathrm{A/m^2}\times\mathrm{m^2}=\mathrm{A}$. ✓ The factor $r$ in $dA$ weights the outer annuli more heavily — there is more area at large $r$ — which is why the result is not simply "the average of $(r/a)^2$ over $[0,1]$", which would be $1/3$.

> [!example] Worked Example 3 — energy stored in a coaxial cable
> Between the conductors of a coaxial cable (inner radius $a$, outer radius $b$, voltage $V_0$, dielectric permittivity $\varepsilon$) the field magnitude is $E=\dfrac{V_0}{r\ln(b/a)}$. Find the stored energy per unit length and hence the capacitance per unit length.
>
> *Energy density:* $w_e=\tfrac12\varepsilon E^2=\dfrac{\varepsilon V_0^2}{2r^2\ln^2(b/a)}$, a function of $r$ only.
>
> *Integrate over a length $\ell$:* $dV=r\,dr\,d\phi\,dz$, so
> $$
> W=\int_0^\ell\!\int_0^{2\pi}\!\!\int_a^b\frac{\varepsilon V_0^2}{2r^2\ln^2(b/a)}\,r\,dr\,d\phi\,dz=\ell\cdot2\pi\cdot\frac{\varepsilon V_0^2}{2\ln^2(b/a)}\int_a^b\frac{dr}{r}=\frac{\pi\varepsilon\ell V_0^2}{\ln(b/a)} .
> $$
> *Capacitance:* equating $W=\tfrac12CV_0^2$ gives $C/\ell=2\pi\varepsilon/\ln(b/a)$.
>
> *Numbers:* $a=0.5\ \mathrm{mm}$, $b=1.75\ \mathrm{mm}$, polyethylene $\varepsilon=2.25\varepsilon_0$: $C/\ell=2\pi(2.25)(8.854\times10^{-12})/\ln3.5=99.9\ \mathrm{pF/m}$, the familiar figure for such cables. Units: $(\mathrm{F/m})\cdot\mathrm{V^2}\cdot\mathrm{m}/1=\mathrm{J}$ for $W$. ✓ The Jacobian $r$ cancelled one power of $1/r^2$, leaving the logarithm.

### Worked examples in spherical coordinates

> [!example] Worked Example 4 — charge in a sphere
> A sphere of radius $a$ has $\rho_v=\rho_0(1-R/a)$. With nothing depending on $\theta$ or $\phi$, the angular integrals give
> $$
> \int_0^{2\pi}\!\!\int_0^\pi\sin\theta\,d\theta\,d\phi=2\pi\cdot2=4\pi ,
> $$
> and
> $$
> Q=4\pi\rho_0\int_0^a\left(1-\frac{R}{a}\right)R^2\,dR=4\pi\rho_0\left(\frac{a^3}{3}-\frac{a^3}{4}\right)=\frac{\pi\rho_0a^3}{3}.
> $$
> This is a quarter of what a uniform density $\rho_0$ would give ($\tfrac43\pi a^3\rho_0$). For a spherically symmetric density, $dV=4\pi R^2\,dR$: the volume of a thin shell of radius $R$ and thickness $dR$.

> [!example] Worked Example 5 — volume of an "ice-cream cone"
> Find the volume of the solid inside the sphere $R=a$ and inside the cone $\theta=\alpha$ (i.e. $0\le\theta\le\alpha$).
> $$
> V=\int_0^{2\pi}\!\!\int_0^\alpha\!\!\int_0^aR^2\sin\theta\,dR\,d\theta\,d\phi=2\pi\cdot(1-\cos\alpha)\cdot\frac{a^3}{3}=\frac{2\pi a^3}{3}(1-\cos\alpha).
> $$
> *Checks:* $\alpha=\pi$ gives $\tfrac43\pi a^3$, the whole sphere; $\alpha=\pi/2$ gives the hemisphere $\tfrac23\pi a^3$. ✓ For $\alpha=\pi/3$, $V=\pi a^3/3$. In Cartesian coordinates the same solid would need limits involving both $\sqrt{a^2-x^2-y^2}$ and $\sqrt{x^2+y^2}\cot\alpha$; in spherical coordinates all three limits are constants.

### Solid angle and antenna patterns

The angular part of the spherical surface element, $d\Omega=\sin\theta\,d\theta\,d\phi$, is the element of **solid angle**: the area of a patch on the unit sphere, measured in steradians (sr). The full sphere subtends $\int d\Omega=4\pi$ sr, and a cone of half-angle $\alpha$ subtends

$$
\Omega(\alpha)=\int_0^{2\pi}\!\!\int_0^\alpha\sin\theta\,d\theta\,d\phi=2\pi(1-\cos\alpha)\ \text{sr},
$$

which gives $2\pi$ sr for a hemisphere. Any function on the sphere — an antenna pattern, the angular distribution of scattered light — is integrated against $d\Omega$.

> [!example] Worked Example 6 — radiated power and directivity
> An antenna's radiation intensity $U(\theta,\phi)$ is the power radiated per steradian (W/sr); the total radiated power is $P=\oint U\,d\Omega=\int_0^{2\pi}\!\int_0^\pi U\sin\theta\,d\theta\,d\phi$, and the **directivity** is $D=4\pi U_{\max}/P$, the ratio of peak intensity to the average over all directions.
>
> (a) *Short (Hertzian) dipole,* $U=U_0\sin^2\theta$:
> $$
> P=2\pi U_0\int_0^\pi\sin^3\theta\,d\theta=2\pi U_0\cdot\frac43=\frac{8\pi U_0}{3},
> $$
> using $\sin^3\theta=(1-\cos^2\theta)\sin\theta$ and $w=\cos\theta$: $\int_{-1}^{1}(1-w^2)\,dw=\tfrac43$. Then $D=4\pi U_0/(8\pi U_0/3)=1.5$ (1.76 dBi), the textbook value for a short dipole.
>
> (b) *A pattern radiating only into the upper hemisphere,* $U=U_0\cos^2\theta$ for $0\le\theta\le\pi/2$ and zero below:
> $$
> P=2\pi U_0\int_0^{\pi/2}\cos^2\theta\sin\theta\,d\theta=2\pi U_0\cdot\frac13=\frac{2\pi U_0}{3},\qquad D=\frac{4\pi U_0}{2\pi U_0/3}=6\ (7.8\ \text{dBi}).
> $$
> Concentrating the same peak intensity into a narrower beam reduces the total power needed — or, for fixed power, raises the peak — and directivity measures exactly that. Numerical quadrature of both angular integrals gives $1.3333$ and $0.3333$. ✓

```desmos-graph
left=-1.2; right=1.2; bottom=-1.2; top=1.2
---
r=\sin(\theta)^2|blue
r=\cos(\theta)^2|red|dashed
```

*Polar cuts through the $z$-axis, drawn with the antenna's $z$-axis along the horizontal: $\sin^2\theta$ (blue) is the figure-eight cross-section of the short dipole's doughnut, with nulls along the axis; $\cos^2\theta$ (red) points along the axis — in Example 6(b) only its right-hand lobe ($0\le\theta\le\pi/2$) is present.*

> [!tip] Engineering Note
> The choice of coordinates is part of the modelling. A cylindrical problem (a long wire) in spherical coordinates is as awkward as a spherical one (a point charge) in Cartesian coordinates. Before setting up any multiple integral, ask: *on which coordinate surfaces are the boundaries, and on which coordinates does the integrand depend?* When both answers point to the same system, the integral usually separates into a product of three single integrals.

## 6.10 Constrained optimization

### Motivation

Design is optimisation: maximise the power delivered to a load, minimise the loss in a network, fit a model to measured data with the least error, choose dimensions that use the least material. In one variable (Chapter 4) we set $f'(x)=0$ and checked the sign of $f''$. With several variables the same ideas work, but "the sign of the second derivative" becomes a property of a matrix, and a new kind of critical point — the saddle — appears. Real problems almost always carry **constraints** as well: a fixed power budget, a fixed total current, a component that must fit in a box. This section treats unconstrained extrema first, then constraints handled by substitution and by searching boundaries; §6.11 gives the general method.

### Local and global extrema

> [!definition] Extrema and critical points
> $f$ has a **local maximum** at $\mathbf{a}$ if $f(\mathbf{x})\le f(\mathbf{a})$ for all $\mathbf{x}$ in some neighbourhood of $\mathbf{a}$, and a **global maximum** if this holds for all $\mathbf{x}$ in the domain (similarly for minima). A **critical point** (stationary point) is an interior point where $\nabla f=\mathbf{0}$, or where $\nabla f$ does not exist.

> [!theorem] First-order necessary condition
> If $f$ is differentiable at an interior point $\mathbf{a}$ and has a local extremum there, then $\nabla f(\mathbf{a})=\mathbf{0}$: every partial derivative vanishes, and the tangent plane is horizontal.

*Proof.* If $f$ has a local maximum at $\mathbf{a}$, then the one-variable function $g(t)=f(\mathbf{a}+t\hat{\mathbf{x}}_i)$ has a local maximum at $t=0$, so $g'(0)=f_{x_i}(\mathbf{a})=0$ by the one-variable result. This holds for each $i$. $\square$

The condition is necessary, not sufficient. At a **saddle point** $\nabla f=\mathbf{0}$, yet $f$ increases in some directions and decreases in others — like the middle of a mountain pass, or the origin for $f=x^2-y^2$.

### The second-derivative test

Near a critical point $\mathbf{a}$ the first-order terms of Taylor's formula (§6.3) vanish, leaving

$$
f(\mathbf{a}+\mathbf{h})-f(\mathbf{a})\approx\tfrac12\,\mathbf{h}^\mathsf{T}\mathbf{H}\,\mathbf{h}=\tfrac12\left(f_{xx}h^2+2f_{xy}hk+f_{yy}k^2\right),\qquad \mathbf{h}=(h,k).
$$

The sign of this **quadratic form** for small $\mathbf{h}$ in every direction decides the type of the point. Write $A=f_{xx}$, $B=f_{xy}$, $C=f_{yy}$ and, if $A\ne0$, complete the square:

$$
Ah^2+2Bhk+Ck^2=A\left(h+\frac BAk\right)^2+\frac{AC-B^2}{A}k^2 .
$$

If $D=AC-B^2>0$, both terms have the sign of $A$ for every $(h,k)\ne\mathbf{0}$: the form is definite, positive (minimum) if $A>0$ and negative (maximum) if $A<0$. If $D<0$, the two terms have opposite signs, so the form takes both signs: a saddle. If $D=0$, the form vanishes along some direction and higher-order terms decide.

> [!theorem] Second-derivative test (two variables)
> Let $\nabla f(\mathbf{a})=\mathbf{0}$, let $f$ have continuous second partials near $\mathbf{a}$, and let $\mathbf{H}=\begin{bmatrix}f_{xx}&f_{xy}\\f_{xy}&f_{yy}\end{bmatrix}$ with determinant $D=f_{xx}f_{yy}-f_{xy}^2$ at $\mathbf{a}$.
>
> | Condition at the critical point | Conclusion |
> |---|---|
> | $D>0,\ f_{xx}>0$ | local minimum |
> | $D>0,\ f_{xx}<0$ | local maximum |
> | $D<0$ | saddle point |
> | $D=0$ | test inconclusive |

In $n$ variables the same reasoning applies to $\mathbf{h}^\mathsf{T}\mathbf{H}\mathbf{h}$: a critical point is a local minimum if $\mathbf{H}$ is **positive definite** (all eigenvalues positive), a local maximum if it is negative definite (all eigenvalues negative), and a saddle if $\mathbf{H}$ has eigenvalues of both signs. Equivalent tests via leading principal minors (Sylvester's criterion) and the eigenvalue theory are in [Ch. 10](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md); in two variables, $D>0$ says the two eigenvalues have the same sign and $f_{xx}$ says which.

> [!example] Worked Example 1 — a minimum and a saddle
> $f(x,y)=x^3-3x+y^2$.
>
> *Critical points:* $f_x=3x^2-3=0$ and $f_y=2y=0$ give $(1,0)$ and $(-1,0)$.
>
> *Second partials:* $f_{xx}=6x$, $f_{yy}=2$, $f_{xy}=0$.
>
> At $(1,0)$: $D=6\cdot2-0=12>0$ and $f_{xx}=6>0$: **local minimum**, $f=-2$.
>
> At $(-1,0)$: $D=-6\cdot2=-12<0$: **saddle point**, $f=2$. Along the $x$-axis $f$ has a local *maximum* at $x=-1$; along the line $x=-1$ it has a *minimum* in $y$.
>
> Neither point is a global extremum: $f\to-\infty$ as $x\to-\infty$.

```desmos-graph
left=-3; right=3; bottom=-2.5; top=2.5
---
x^3-3x+y^2=2|red
x^3-3x+y^2=-1|blue
x^3-3x+y^2=-1.8|blue
x^3-3x+y^2=0.5|grey
x^3-3x+y^2=3.5|grey
(1,0)|label:min
(-1,0)|label:saddle
```

*Level curves of $x^3-3x+y^2$: small closed loops (blue) surround the minimum at $(1,0)$, while the level curve through the saddle (red, $f=2$) crosses itself at $(-1,0)$ — the signature of a saddle on a contour map.*

> [!example] Worked Example 2 — three critical points
> $f(x,y)=x^4+y^4-4xy$.
>
> *Critical points:* $f_x=4x^3-4y=0$ gives $y=x^3$; $f_y=4y^3-4x=0$ gives $x=y^3=x^9$. So $x(x^8-1)=0$: $x=0,\pm1$, giving $(0,0)$, $(1,1)$, $(-1,-1)$.
>
> *Hessian:* $\mathbf{H}=\begin{bmatrix}12x^2&-4\\-4&12y^2\end{bmatrix}$.
>
> At $(0,0)$: $D=0-16=-16<0$: saddle, $f=0$.
>
> At $(\pm1,\pm1)$: $D=144-16=128>0$, $f_{xx}=12>0$: local minima, $f=1+1-4=-2$.
>
> Since $f\to+\infty$ in every direction (the quartic terms dominate), these two minima are also the **global** minima. Eigenvalue check at $(1,1)$: $\mathbf{H}=\begin{bmatrix}12&-4\\-4&12\end{bmatrix}$ has eigenvalues $8$ and $16$, both positive. ✓

> [!example] Worked Example 3 — least-squares line fit
> A sensor is calibrated at four input levels, giving output voltages $(x_i,y_i)=(0,1),(1,3),(2,4),(3,7)$. Find the line $y=mx+b$ that minimises the sum of squared errors
> $$
> S(m,b)=\sum_{i=1}^4\left(y_i-mx_i-b\right)^2 .
> $$
> *Critical point:*
> $$
> \frac{\partial S}{\partial m}=-2\sum x_i(y_i-mx_i-b)=0,\qquad \frac{\partial S}{\partial b}=-2\sum(y_i-mx_i-b)=0,
> $$
> which are the **normal equations**
> $$
> m\sum x_i^2+b\sum x_i=\sum x_iy_i,\qquad m\sum x_i+nb=\sum y_i .
> $$
> With $\sum x_i=6$, $\sum x_i^2=14$, $\sum y_i=15$, $\sum x_iy_i=32$, $n=4$: $14m+6b=32$ and $6m+4b=15$, so $m=1.9$, $b=0.9$.
>
> *Is it a minimum?* $\mathbf{H}=2\begin{bmatrix}\sum x_i^2&\sum x_i\\\sum x_i&n\end{bmatrix}=\begin{bmatrix}28&12\\12&8\end{bmatrix}$, $D=224-144=80>0$, $S_{mm}=28>0$: yes — and since $S$ is a quadratic with positive-definite Hessian, it is the global minimum.
>
> *Residuals:* $0.1,\ 0.2,\ -0.7,\ 0.4$, which sum to zero (that is the second normal equation), with $S_{\min}=0.70\ \mathrm{V^2}$. The same calculation with many parameters is linear least squares, treated with matrices in [Ch. 9](../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md) and statistically in [Ch. 20](../part-6-probability-statistics/20-statistical-analysis.md).

> [!example] Worked Example 4 — maximum power transfer to a complex load
> A source with open-circuit RMS voltage $V_s$ and internal impedance $Z_s=R_s+jX_s$ drives a load $Z_L=R_L+jX_L$. The load power is
> $$
> P(R_L,X_L)=\frac{V_s^2R_L}{(R_s+R_L)^2+(X_s+X_L)^2},\qquad R_L>0 .
> $$
> *Reactance first:* $P_{X_L}=-\dfrac{2V_s^2R_L(X_s+X_L)}{\left[(R_s+R_L)^2+(X_s+X_L)^2\right]^2}=0$ gives $X_L=-X_s$.
>
> *Then resistance:* with $X_L=-X_s$, $P=V_s^2R_L/(R_s+R_L)^2$ and $P_{R_L}=V_s^2\dfrac{(R_s+R_L)-2R_L}{(R_s+R_L)^3}=0$ gives $R_L=R_s$.
>
> So the only critical point is $Z_L=Z_s^*$ — the **conjugate match** — with $P_{\max}=V_s^2/(4R_s)$. For $V_s=10\ \mathrm{V}$ and $Z_s=50+j30\ \Omega$: $Z_L=50-j30\ \Omega$ and $P_{\max}=100/200=0.5\ \mathrm{W}$.
>
> *Classification:* at the critical point one finds $P_{R_LR_L}=P_{X_LX_L}=-V_s^2/(8R_s^3)<0$ and $P_{R_LX_L}=0$, so $D=V_s^4/(64R_s^6)>0$ with $P_{R_LR_L}<0$: a maximum. (Units: $\mathrm{W}/\Omega^2$, as a second derivative of power with respect to resistance must be.) Physically obvious too — $P\to0$ as $R_L\to0$, as $R_L\to\infty$ or as $\lvert X_L\rvert\to\infty$, so the single critical point must be the global maximum.

> [!warning] Common Pitfall
> Checking only one variable at a time. A point can be a minimum along the $x$-direction and along the $y$-direction and still be a saddle: $f=x^2+y^2-3xy$ has $f_{xx}=f_{yy}=2>0$ at the origin, yet $D=4-9<0$, and along $y=x$, $f=-x^2$ decreases. The mixed derivative $f_{xy}$ — the twist — must be included, which is exactly what $D$ does.

### Constraints by substitution

A **constrained** problem reads

$$
\text{optimise } f(\mathbf{x})\quad\text{subject to}\quad g(\mathbf{x})=c .
$$

If the constraint can be solved for one variable, substitute and reduce the problem to an unconstrained one in fewer variables.

> [!example] Worked Example 5 — an enclosure of minimum material
> An open-top rectangular enclosure must hold $32\ \mathrm{cm^3}$. Minimise the sheet-metal area $A=xy+2xz+2yz$ subject to $xyz=32$.
>
> *Substitute* $z=32/(xy)$:
> $$
> A(x,y)=xy+\frac{64}{y}+\frac{64}{x},\qquad x,y>0 .
> $$
> *Critical point:* $A_x=y-64/x^2=0$, $A_y=x-64/y^2=0$. So $yx^2=64$ and $xy^2=64$; dividing, $x=y$, then $x^3=64$, $x=y=4\ \mathrm{cm}$, $z=2\ \mathrm{cm}$.
>
> *Classify:* $A_{xx}=128/x^3=2$, $A_{yy}=2$, $A_{xy}=1$; $D=4-1=3>0$, $A_{xx}>0$: minimum, $A=16+32+32=48\ \mathrm{cm^2}$. The open top makes the optimal box half as tall as it is wide.

### Global extrema on closed, bounded regions

> [!theorem] Extreme value theorem
> A continuous function on a **closed and bounded** region attains a global maximum and a global minimum there. They occur either at interior critical points or on the boundary.

Hence the **closed-region method**:

1. Find the critical points of $f$ inside the region and evaluate $f$ there.
2. Find the extreme values of $f$ on the boundary — a constrained problem, solved by parametrising or substituting (or by Lagrange multipliers, §6.11). Include corners and endpoints of boundary pieces.
3. The largest value from steps 1–2 is the global maximum, the smallest the global minimum.

> [!example] Worked Example 6 — extremes on a disc
> Find the extreme values of $T(x,y)=x^2+2y^2-x$ on the disc $x^2+y^2\le1$.
>
> *Interior:* $T_x=2x-1=0$, $T_y=4y=0$ give $(\tfrac12,0)$, with $T=-\tfrac14$.
>
> *Boundary:* substitute $y^2=1-x^2$: $T=x^2+2(1-x^2)-x=2-x-x^2$ for $-1\le x\le1$. Then $dT/dx=-1-2x=0$ gives $x=-\tfrac12$, $y=\pm\tfrac{\sqrt3}{2}$, $T=2+\tfrac12-\tfrac14=\tfrac94$. The endpoints of the $x$-interval give $T(1,0)=0$ and $T(-1,0)=2$.
>
> *Compare:* $\{-\tfrac14,\ \tfrac94,\ 0,\ 2\}$: global maximum $\tfrac94$ at $(-\tfrac12,\pm\tfrac{\sqrt3}{2})$; global minimum $-\tfrac14$ at $(\tfrac12,0)$.
>
> Alternatively, parametrise the circle as $(\cos t,\sin t)$: $T=\cos^2t+2\sin^2t-\cos t=1+\sin^2t-\cos t$, and $dT/dt=\sin t(2\cos t+1)=0$ gives the same four boundary candidates.

Substitution becomes awkward when the constraint cannot be solved explicitly (try solving $x^3+y^3+xy=5$ for $y$), when it breaks the symmetry between variables, or when there are several constraints. Lagrange's method, next, avoids all three problems.

> [!tip] Engineering Note
> Most practical optimisation problems are solved numerically ([Ch. 23](../part-7-numerical-computational/23-optimization.md)), but the analysis of this section remains essential: it tells you what the optimiser should find, whether the problem has one solution or several, and whether the point it returns is a minimum or merely a saddle where the gradient happens to vanish. Convex problems — those whose Hessian is positive semidefinite everywhere, like least squares — have no saddles and no spurious local minima, which is why engineers work hard to formulate problems in convex form.

## 6.11 Lagrange multipliers

### Motivation

Substitution solves a constrained problem by eliminating variables, which requires solving the constraint and usually destroys any symmetry the problem had. Lagrange's method keeps all the variables and treats them even-handedly, handles several constraints at once, and — as a bonus — produces a number, the multiplier, that measures how much the optimum would improve if the constraint were relaxed. In engineering that number is often the most useful output of the whole calculation: it is the marginal value of an extra watt of power budget, an extra square centimetre of board area, an extra volt of supply.

### The geometric idea

Consider maximising $f(x,y)$ along the curve $g(x,y)=c$. Walk along the constraint curve and watch the level curves of $f$ that you cross. If at some point you are *crossing* a level curve of $f$ transversally, then moving one way along the constraint increases $f$ and the other way decreases it — so that point is not an extremum. At a constrained extremum, therefore, the level curve of $f$ through the point must be **tangent** to the constraint curve. Tangent curves have parallel normals, and the normals are the gradients (§6.4). Hence at the extremum

$$
\nabla f=\lambda\,\nabla g
$$

for some scalar $\lambda$. The same argument in three dimensions says that $\nabla f$ must be perpendicular to the constraint surface, i.e. parallel to its normal $\nabla g$.

An equivalent analytic view: let $\mathbf{r}(t)$ be any curve lying in the constraint set and passing through the extremum $\mathbf{x}^*$ at $t=0$. Then $f(\mathbf{r}(t))$ has an extremum at $t=0$, so $\nabla f(\mathbf{x}^*)\cdot\mathbf{r}'(0)=0$. Thus $\nabla f$ is perpendicular to every tangent direction of the constraint set — and the only vectors with that property are multiples of $\nabla g$ (or, with several constraints, combinations of the $\nabla g_k$).

### The method

> [!theorem] Method of Lagrange multipliers
> Let $f$ and $g$ have continuous first partial derivatives. If $f$ has a local extremum at $\mathbf{x}^*$ subject to $g(\mathbf{x})=c$, and $\nabla g(\mathbf{x}^*)\ne\mathbf{0}$, then there is a scalar $\lambda$ such that
> $$
> \nabla f(\mathbf{x}^*)=\lambda\,\nabla g(\mathbf{x}^*),\qquad g(\mathbf{x}^*)=c .
> $$
> Equivalently, $(\mathbf{x}^*,\lambda)$ is a stationary point of the **Lagrangian**
> $$
> \mathcal{L}(\mathbf{x},\lambda)=f(\mathbf{x})-\lambda\,\big(g(\mathbf{x})-c\big):
> $$
> setting $\partial\mathcal{L}/\partial x_i=0$ reproduces $\nabla f=\lambda\nabla g$, and $\partial\mathcal{L}/\partial\lambda=0$ reproduces the constraint.
>
> With several constraints $g_1=c_1,\dots,g_m=c_m$ (and the $\nabla g_k$ linearly independent at $\mathbf{x}^*$), the condition becomes
> $$
> \nabla f=\lambda_1\nabla g_1+\dots+\lambda_m\nabla g_m,\qquad g_k=c_k .
> $$

For $n$ variables and $m$ constraints this is $n+m$ equations in $n+m$ unknowns. The procedure:

1. Write down $\nabla f=\lambda\nabla g$ component by component, together with the constraint.
2. Solve the system. There is no general algorithm — the equations are nonlinear — but useful tricks are: divide equations to eliminate $\lambda$; consider separately the cases where a factor is zero; exploit symmetry.
3. Evaluate $f$ at **every** candidate and compare. Also check points where $\nabla g=\mathbf{0}$ and, for constraint curves with ends, the endpoints.

> [!example] Worked Example 1 — a linear function on a circle
> Extremise $f=x+2y$ subject to $x^2+y^2=5$.
>
> *Conditions:* $\nabla f=(1,2)=\lambda(2x,2y)$ gives $1=2\lambda x$ and $2=2\lambda y$. Neither $x$ nor $\lambda$ can be zero, so dividing, $y/x=2$, i.e. $y=2x$.
>
> *Constraint:* $x^2+4x^2=5$, $x=\pm1$.
>
> *Candidates:* $(1,2)$ with $f=5$ (maximum) and $(-1,-2)$ with $f=-5$ (minimum); $\lambda=1/(2x)=\pm\tfrac12$.
>
> The answer is the Cauchy–Schwarz inequality in action: $\lvert(1,2)\cdot(x,y)\rvert\le\lVert(1,2)\rVert\,\lVert(x,y)\rVert=\sqrt5\cdot\sqrt5=5$, with equality when $(x,y)$ is parallel to $(1,2)$.

```desmos-graph
left=-4; right=4; bottom=-3; top=3
---
x^2+y^2=5|blue
x+2y=5|red
x+2y=-5|red
x+2y=0|grey|dashed
x+2y=2.5|grey|dashed
(1,2)|label:max
(-1,-2)|label:min
```

*Level lines of $f=x+2y$ (red and grey) sweep across the circle; the extremes occur where a level line just touches it — there the gradients of $f$ and $g$ are parallel.*

> [!example] Worked Example 2 — the disc boundary revisited
> On $x^2+y^2=1$, extremise $T=x^2+2y^2-x$ (the boundary step of §6.10, Worked Example 6).
>
> *Conditions:* $(2x-1,\ 4y)=\lambda(2x,\ 2y)$. The second component gives $4y=2\lambda y$, so **either** $y=0$ **or** $\lambda=2$.
>
> *Case $y=0$:* the constraint gives $x=\pm1$, with $T(1,0)=0$ and $T(-1,0)=2$.
>
> *Case $\lambda=2$:* the first component gives $2x-1=4x$, $x=-\tfrac12$, and the constraint $y=\pm\tfrac{\sqrt3}{2}$, with $T=\tfrac94$.
>
> These are the same four candidates found by substitution, but obtained without breaking the symmetry of the constraint — and the case split "$y=0$ or $\lambda=2$" is the standard way such systems open up. Note that the substitution method found $(\pm1,0)$ only as *endpoints* of the $x$-interval; Lagrange's method finds them as genuine critical points of the constrained problem.

> [!example] Worked Example 3 — closest point on a plane
> Find the point of the plane $x+2y+2z=9$ closest to the origin.
>
> Minimise the squared distance $f=x^2+y^2+z^2$ (same minimiser as the distance, simpler derivatives).
>
> *Conditions:* $(2x,2y,2z)=\lambda(1,2,2)$, so $(x,y,z)=\tfrac{\lambda}{2}(1,2,2)$: the closest point lies along the plane's normal, as geometry demands.
>
> *Constraint:* $\tfrac{\lambda}{2}(1+4+4)=9$, $\lambda=2$, giving $(1,2,2)$ at distance $3$.
>
> *Check:* the point-to-plane distance formula gives $9/\sqrt{1^2+2^2+2^2}=3$. ✓ No comparison of candidates is needed: there is only one, and $f$ grows without bound on the plane, so it must be the minimum.

> [!example] Worked Example 4 — the current divider minimises dissipation
> A current $I$ splits between resistors $R_1$ and $R_2$ in parallel. Ask which split $I_1+I_2=I$ **minimises** the total dissipation $P=R_1I_1^2+R_2I_2^2$.
>
> *Conditions:*
> $$
> 2R_1I_1=\lambda,\qquad 2R_2I_2=\lambda,\qquad I_1+I_2=I .
> $$
> Hence $R_1I_1=R_2I_2$ — equal voltages across both branches — and
> $$
> I_1=\frac{R_2}{R_1+R_2}I,\qquad I_2=\frac{R_1}{R_1+R_2}I,\qquad P_{\min}=\frac{R_1R_2}{R_1+R_2}I^2 .
> $$
> That is exactly the current-divider rule: in a resistive network, currents distribute themselves so as to minimise heat (Maxwell's minimum-heat theorem). Here $\lambda=2R_1I_1=2IR_1R_2/(R_1+R_2)$, twice the common branch voltage, in volts.
>
> *Numbers:* $R_1=1\ \mathrm{k\Omega}$, $R_2=4\ \mathrm{k\Omega}$, $I=10\ \mathrm{mA}$: $I_1=8\ \mathrm{mA}$, $I_2=2\ \mathrm{mA}$, $P_{\min}=800\ \Omega\times10^{-4}\ \mathrm{A^2}=80\ \mathrm{mW}$. An equal split ($5\ \mathrm{mA}$ each) would dissipate $25+100=125\ \mathrm{mW}$.

### The meaning of $\lambda$

> [!theorem] Sensitivity interpretation of the multiplier
> Let $\mathbf{x}^*(c)$ be the constrained optimum as a function of the constraint level $c$, with multiplier $\lambda(c)$, and let $f^*(c)=f(\mathbf{x}^*(c))$ be the optimal value. If $\mathbf{x}^*(c)$ depends differentiably on $c$, then
> $$
> \lambda=\frac{df^*}{dc}.
> $$

*Proof.* By the chain rule, $\dfrac{df^*}{dc}=\nabla f(\mathbf{x}^*)\cdot\dfrac{d\mathbf{x}^*}{dc}=\lambda\,\nabla g(\mathbf{x}^*)\cdot\dfrac{d\mathbf{x}^*}{dc}$. But $g(\mathbf{x}^*(c))=c$ for all $c$, so differentiating, $\nabla g\cdot\dfrac{d\mathbf{x}^*}{dc}=1$. Hence $df^*/dc=\lambda$. $\square$

So $\lambda$ is the **shadow price** of the constraint: how much the optimal value improves per unit relaxation of the constraint. Its units are (units of $f$)/(units of $g$).

*Check on Worked Example 1:* with $x^2+y^2=c$, the maximum is $f^*=\sqrt{5c}$, so $df^*/dc=\sqrt5/(2\sqrt c)=\tfrac12$ at $c=5$ — and indeed $\lambda=\tfrac12$ at $(1,2)$. ✓

*Check on Worked Example 4:* $dP_{\min}/dI=2IR_1R_2/(R_1+R_2)=\lambda$. ✓ With the numbers above, $\lambda=16\ \mathrm{V}$: each extra milliampere of total current costs about $16\ \mathrm{mW}$ of extra dissipation.

*Check on the enclosure of §6.10:* with $\mathcal{L}=xy+2xz+2yz-\lambda(xyz-V)$, the conditions $y+2z=\lambda yz$, $x+2z=\lambda xz$, $2x+2y=\lambda xy$ at $(4,4,2)$ give $\lambda=1\ \mathrm{cm^{-1}}$. Scaling says the optimal area grows as $A^*=48\,(V/32)^{2/3}$, so $dA^*/dV=\tfrac23\cdot\tfrac{48}{32}=1\ \mathrm{cm^2/cm^3}$. ✓

### Two constraints

> [!example] Worked Example 5 — extremes on the intersection of a plane and a cylinder
> Find the maximum and minimum of $f=x+2y+3z$ on the curve where the plane $x-y+z=1$ meets the cylinder $x^2+y^2=1$ (an ellipse).
>
> *Conditions:* $\nabla f=\lambda\nabla g_1+\mu\nabla g_2$ with $g_1=x-y+z$, $g_2=x^2+y^2$:
> $$
> 1=\lambda+2\mu x,\qquad 2=-\lambda+2\mu y,\qquad 3=\lambda .
> $$
> *Solve:* $\lambda=3$, so $2\mu x=-2$ and $2\mu y=5$: $x=-1/\mu$, $y=5/(2\mu)$. The cylinder gives $\dfrac{1}{\mu^2}+\dfrac{25}{4\mu^2}=1$, $\mu^2=\dfrac{29}{4}$, $\mu=\pm\dfrac{\sqrt{29}}{2}$.
>
> Then $x=\mp\dfrac{2}{\sqrt{29}}$, $y=\pm\dfrac{5}{\sqrt{29}}$, and the plane gives $z=1-x+y=1\pm\dfrac{7}{\sqrt{29}}$.
>
> *Values:* $f=x+2y+3z=\pm\dfrac{-2+10+21}{\sqrt{29}}+3=3\pm\sqrt{29}$.
>
> So the maximum is $3+\sqrt{29}\approx8.385$ at $(-0.3714,\ 0.9285,\ 2.2999)$ and the minimum $3-\sqrt{29}\approx-2.385$. A numerical constrained optimiser started from $(1,0,0)$ returns the same maximiser to six figures. ✓

### Second-order check: the bordered Hessian

Lagrange's condition is necessary, not sufficient. When comparing candidate values is impractical, a second-order test is available. For two variables and one constraint, form the **bordered Hessian**

$$
\bar{\mathbf{H}}=\begin{bmatrix}0&g_x&g_y\\ g_x&\mathcal{L}_{xx}&\mathcal{L}_{xy}\\ g_y&\mathcal{L}_{xy}&\mathcal{L}_{yy}\end{bmatrix}
$$

at the candidate. If $\det\bar{\mathbf{H}}>0$ the point is a constrained local maximum; if $\det\bar{\mathbf{H}}<0$, a constrained local minimum. (The test asks whether the Hessian of $\mathcal{L}$ is negative or positive definite *along the constraint's tangent direction*, which is all that matters.) For Worked Example 1 at $(1,2)$ with $\lambda=\tfrac12$: $\mathcal{L}_{xx}=\mathcal{L}_{yy}=-2\lambda=-1$, $\mathcal{L}_{xy}=0$, $g_x=2$, $g_y=4$, and $\det\bar{\mathbf{H}}=0-2(2\cdot(-1)-0)+4(0-4\cdot(-1))=4+16=20>0$: a maximum. ✓

> [!warning] Common Pitfall
> Three ways to lose a solution. (1) Dividing by a quantity that might be zero — in Worked Example 2, dividing $4y=2\lambda y$ by $y$ loses the candidates with $y=0$; factor instead and treat each case. (2) Ignoring points where $\nabla g=\mathbf{0}$: the theorem says nothing there (minimise $f=x$ on the cusp curve $y^2=x^3$; the minimum at the origin, where $\nabla g=\mathbf{0}$, does not satisfy $\nabla f=\lambda\nabla g$). (3) Forgetting endpoints of a constraint curve that is not closed. And always remember that a Lagrange point may be a maximum, a minimum or neither — compare values.

### Inequality constraints and water-filling

Engineering limits are usually inequalities — a power budget $\sum p_i\le P$, a current limit, non-negative powers $p_i\ge0$. The extension of Lagrange's method is the set of **Karush–Kuhn–Tucker (KKT) conditions**: each inequality gets a multiplier that must be non-negative, and that multiplier is zero unless its constraint is **active** (holds with equality) — "complementary slackness". In effect, inactive constraints are ignored and active ones are treated as equalities; the work lies in finding which are active. The full theory is in [Ch. 23](../part-7-numerical-computational/23-optimization.md); the classic EE example shows the idea.

> [!example] Worked Example 6 — water-filling power allocation
> A transmitter shares total power $P$ between two parallel channels with noise powers $N_1=1$ and $N_2=3$ (in the same units as power). The total capacity in bit/s/Hz is
> $$
> C=\log_2\!\left(1+\frac{p_1}{N_1}\right)+\log_2\!\left(1+\frac{p_2}{N_2}\right),\qquad p_1+p_2=P,\ \ p_1,p_2\ge0 .
> $$
> *Lagrange conditions (both channels active):* $\dfrac{1}{\ln2}\cdot\dfrac{1}{N_i+p_i}=\lambda$, so $N_i+p_i=\nu$ with the same "water level" $\nu=1/(\lambda\ln2)$ for both channels.
>
> *With $P=6$:* $p_1=\nu-1$, $p_2=\nu-3$, sum $2\nu-4=6$, $\nu=5$: $p_1=4$, $p_2=2$, both non-negative, so both constraints $p_i\ge0$ are inactive and the solution stands. Then $C=\log_25+\log_2\tfrac53=2.322+0.737=3.059$ bit/s/Hz, better than an equal split ($p_1=p_2=3$: $\log_24+\log_22=3$ bit/s/Hz).
>
> *With $P=1$:* the same formulas would give $\nu=2.5$ and $p_2=-0.5<0$ — impossible. So the constraint $p_2\ge0$ is active: $p_2=0$, and all power goes to the quieter channel, $p_1=1$, with water level $\nu=N_1+p_1=2$, below the noise "floor" $N_2=3$ of channel 2. Then $C=\log_22=1$ bit/s/Hz.
>
> The picture: pour total power $P$ like water into a vessel whose floor heights are the noise levels $N_i$; it settles at a common level $\nu$, and channels whose floor lies above the water get nothing. The rule "pour into the quietest channels first" is how DSL, Wi-Fi (OFDM) and MIMO systems allocate power across sub-channels.

```desmos-graph
left=0; right=8; bottom=0; top=4
---
y=\log_2(1+x)+\log_2(1+\frac{6-x}{3})|blue|0<x<6
x=4|red|dashed
(4,3.0589)|label:optimum p1=4
(3,3)|label:equal split
```

*Capacity versus the power $p_1$ given to the quieter channel (total $P=6$): the curve is concave with a single maximum at $p_1=4$, where the marginal capacity per watt — the multiplier $\lambda$ — is the same in both channels.*

> [!tip] Engineering Note
> The multiplier is often the quantity of real interest. In economic dispatch of power stations, the multiplier on the "generation equals demand" constraint is the system marginal cost (the price of one more MWh). In the water-filling example it is the marginal capacity per watt. In circuit design, it tells you which constraint — area, power, a noise specification — is binding and how much it is costing you, so you know where relaxing a specification would pay off most.

## Summary

| Concept | Key formula |
|---|---|
| Limit along all paths | polar test: $\lvert f-L\rvert\le g(r)\to0$ independent of $\theta$ |
| Partial derivative | $f_x=\lim_{h\to0}\dfrac{f(x+h,y)-f(x,y)}{h}$ (other variables fixed) |
| Normalised sensitivity | $S^f_x=\dfrac{x}{f}\dfrac{\partial f}{\partial x}=\dfrac{\partial\ln f}{\partial\ln x}$ |
| Clairaut | $f_{xy}=f_{yx}$ (continuous second partials) |
| Total differential | $df=f_x\,dx+f_y\,dy=\nabla f\cdot d\mathbf{r}$ |
| Tangent plane | $z=f(a,b)+f_x(a,b)(x-a)+f_y(a,b)(y-b)$ |
| Worst-case / RSS error | $\sum\lvert f_{x_i}\rvert\Delta x_i$ / $\sqrt{\sum(f_{x_i}\Delta x_i)^2}$ |
| Relative error, $f=kx^ay^b$ | $df/f=a\,dx/x+b\,dy/y$ |
| Chain rule | $\dfrac{dz}{dt}=f_x\dot x+f_y\dot y$; $\ \mathbf{J}_{\mathbf{F}\circ\mathbf{G}}=\mathbf{J}_{\mathbf{F}}\mathbf{J}_{\mathbf{G}}$ |
| Implicit derivative | $dy/dx=-F_x/F_y$; $\ \partial z/\partial x=-F_x/F_z$ |
| Exact differential test | $M\,dx+N\,dy$ exact $\iff M_y=N_x$ |
| Second-order Taylor | $f(\mathbf{a}+\mathbf{h})\approx f+\nabla f\cdot\mathbf{h}+\tfrac12\mathbf{h}^\mathsf{T}\mathbf{H}\mathbf{h}$ |
| Gradient | $\nabla f=(f_x,f_y,f_z)$; normal to level sets; $\lVert\nabla f\rVert$ = max rate |
| Useful gradients | $\nabla R=\hat{\mathbf{R}}$, $\ \nabla(1/R)=-\hat{\mathbf{R}}/R^2$, $\ \mathbf{E}=-\nabla V$ |
| Gradient descent | $\mathbf{x}_{k+1}=\mathbf{x}_k-\alpha\nabla f(\mathbf{x}_k)$ |
| Directional derivative | $D_{\hat{\mathbf{u}}}f=\nabla f\cdot\hat{\mathbf{u}}=\lVert\nabla f\rVert\cos\theta$ |
| Jacobian / linearisation | $(\mathbf{J})_{ij}=\partial F_i/\partial x_j$; $\ d\mathbf{F}=\mathbf{J}\,d\mathbf{x}$ |
| Inverse function | $\mathbf{J}_{\mathbf{F}^{-1}}=\mathbf{J}_{\mathbf{F}}^{-1}$ |
| Solution sensitivity | $\partial\mathbf{x}/\partial\mathbf{p}=-\mathbf{J}_{\mathbf{x}}^{-1}\mathbf{J}_{\mathbf{p}}$ |
| Newton's method | $\mathbf{J}(\mathbf{x}_k)\Delta\mathbf{x}=-\mathbf{F}(\mathbf{x}_k)$ |
| Fubini (type I) | $\iint_R f\,dA=\int_a^b\int_{g_1(x)}^{g_2(x)}f\,dy\,dx$ |
| Change of variables | $dx\,dy=\left\lvert\dfrac{\partial(x,y)}{\partial(u,v)}\right\rvert du\,dv$ |
| Polar | $dA=r\,dr\,d\theta$; $\ \int_{-\infty}^\infty e^{-x^2}dx=\sqrt\pi$ |
| Cylindrical | $dV=r\,dr\,d\phi\,dz$; scale factors $(1,r,1)$ |
| Spherical | $dV=R^2\sin\theta\,dR\,d\theta\,d\phi$; scale factors $(1,R,R\sin\theta)$ |
| Solid angle | $d\Omega=\sin\theta\,d\theta\,d\phi$; cone: $2\pi(1-\cos\alpha)$; sphere: $4\pi$ |
| Second-derivative test | $D=f_{xx}f_{yy}-f_{xy}^2$: $D>0$ extremum (sign of $f_{xx}$), $D<0$ saddle |
| Lagrange multipliers | $\nabla f=\lambda\nabla g,\ g=c$; several: $\nabla f=\sum\lambda_k\nabla g_k$ |
| Multiplier sensitivity | $\lambda=df^*/dc$ |
| Water-filling | $p_i=\max(0,\ \nu-N_i)$, $\ \sum p_i=P$ |

## Further reading

- E. Kreyszig, *Advanced Engineering Mathematics*, 10th ed., Wiley — Ch. 9–10 for multivariable and vector calculus at engineering level.
- K. A. Stroud and D. J. Booth, *Advanced Engineering Mathematics*, 6th ed., Palgrave Macmillan — programme-style treatment of partial differentiation, multiple integrals and Lagrange multipliers with many drills.
- J. Stewart, *Multivariable Calculus*, 8th ed., Cengage — clear geometric treatment of gradients, change of variables and constrained optimisation.
- D. K. Cheng, *Field and Wave Electromagnetics*, 2nd ed., Addison-Wesley — Ch. 2 for orthogonal curvilinear coordinates, scale factors and unit vectors as used in field theory.
- S. Boyd and L. Vandenberghe, *Convex Optimization*, Cambridge University Press — Lagrange duality, KKT conditions and water-filling, from the engineering viewpoint.
