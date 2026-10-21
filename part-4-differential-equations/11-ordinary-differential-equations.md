# Chapter 11 — Ordinary Differential Equations

[← Previous](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md) · [Table of Contents](../00-front/table-of-contents.md) · [Next →](12-linear-differential-equations.md)

> **Prerequisites:** Ch. 4 (Differential Calculus), Ch. 5 (Integral Calculus), Ch. 6 (Multivariable Calculus), Ch. 8–10 (Matrices, Vector Spaces, Eigenvalues)
>
> **Learning outcomes:** After this chapter you can
> - classify an ODE by order, linearity, autonomy and homogeneity, and predict the kind of solution to expect before computing anything;
> - solve first-order equations that are separable, linear, exact (possibly after an integrating factor), Bernoulli or Riccati, and read their qualitative behaviour from a direction field or phase line;
> - use the structure theory of linear $n$th-order equations — fundamental sets, the Wronskian and Abel's formula, reduction of order, Cauchy–Euler equations, $y = y_h + y_p$ and variation of parameters;
> - state and apply the existence–uniqueness theorems for initial-value problems, and recognise blow-up and non-uniqueness when the hypotheses fail;
> - decide whether a boundary-value problem has one, none or infinitely many solutions, and solve simple Sturm–Liouville eigenvalue problems;
> - write any $n$th-order equation or lumped network as a first-order system $\dot{\mathbf{x}} = \mathbf{A}\mathbf{x} + \mathbf{f}$, solve the constant-coefficient case by eigenvectors, and classify equilibria in the phase plane.

A differential equation relates an unknown function to its own derivatives. It is the natural language of anything that *stores* something — charge on a capacitor, flux in an inductor, heat in a package, carriers in a semiconductor, momentum in a motor shaft — because a storage law says how fast the stored quantity changes, not what its value is. Kirchhoff's laws plus the element laws $i_C = C\,dv/dt$ and $v_L = L\,di/dt$ turn every lumped circuit into a set of ODEs.

This chapter gives the **general theory** of ordinary differential equations (ODEs): how they are classified, which solution methods apply to which class, and when a solution exists, is unique and lasts for all time. The emphasis is on methods that work for *variable* coefficients and *nonlinear* equations, because those are what the special-purpose tools cannot handle. The special case that dominates circuit analysis — linear equations with **constant coefficients**, their characteristic equations, natural and forced responses and state-variable form — is developed fully in [Chapter 12](12-linear-differential-equations.md); the transform route to the same equations is [Chapter 14](../part-5-transforms/14-laplace-transforms.md); numerical integration is [Chapter 22](../part-7-numerical-computational/22-numerical-methods.md). Where those chapters take over, this one gives a pointer rather than a duplicate.

---

## 11.1 First-order ODEs

### Why start with first order

A first-order equation is one storage element, one rate law: an RC or RL network, a thermal mass, a single-pole amplifier, a charging battery, a population. Every higher-order equation can be rewritten as a *system* of first-order equations (§11.5, §11.9), so the ideas met here — direction fields, equilibria, existence and uniqueness — are the ideas of the whole subject in their simplest setting.

### Definitions and classification

> [!definition] Ordinary differential equation
> An **ordinary differential equation** of order $n$ is a relation
> $$
> F\big(t, y, y', y'', \dots, y^{(n)}\big) = 0
> $$
> between one independent variable $t$ (time) or $x$ (position), an unknown function $y$ and its derivatives up to order $n$. The **order** is the highest derivative present. "Ordinary" means a single independent variable; with several, the equation is a *partial* differential equation ([Ch. 13](13-partial-differential-equations.md)).

A **first-order** ODE is $F(t, y, y') = 0$. Whenever it can be solved for $y'$ we use the **normal form**

$$
y' = f(t, y).
$$

The vocabulary below is used throughout the book.

- **Solution on an interval.** A function $y(t)$ is a solution on an open interval $I$ if it is differentiable on $I$ and $y'(t) = f(t, y(t))$ for every $t \in I$. The interval is part of the answer: a formula may satisfy the equation only on a piece of the real line.
- **General solution.** A family $y = \varphi(t, C)$ containing one arbitrary constant, such that (almost) every solution is obtained for some value of $C$. Its graphs form a one-parameter family of **integral curves**.
- **Particular solution.** One member of the family, usually picked out by an **initial condition** $y(t_0) = y_0$.
- **Singular solution.** A solution not obtainable from the general solution for any finite $C$. For $y' = y^2$ the general solution is $y = 1/(C - t)$ and $y \equiv 0$ is singular.
- **Explicit vs implicit.** An explicit solution gives $y$ as a function of $t$; an implicit one is a relation $G(t, y) = 0$ that defines $y(t)$ only implicitly (by the implicit-function theorem of [Ch. 6](../part-2-calculus/06-multivariable-calculus.md), wherever $G_y \ne 0$).
- **Linear vs nonlinear.** The equation is **linear** if it can be written $a_1(t)\,y' + a_0(t)\,y = r(t)$ — the unknown and its derivative appear to the first power, are not multiplied together, and are not inside other functions. Otherwise it is **nonlinear**.
- **Autonomous.** The equation is autonomous if $f$ does not depend explicitly on $t$: $y' = f(y)$. Time-invariant circuits with constant sources are autonomous.

| Equation | Order | Linear? | Autonomous? |
|---|---|---|---|
| $RC\,\dot v + v = V_s$ | 1 | yes | yes |
| $\dot y + t\,y = \sin t$ | 1 | yes | no |
| $C\,\dot v = -I_S\big(e^{v/V_T} - 1\big)$ | 1 | no | yes |
| $y' = y^2 + t$ | 1 | no | no |
| $LC\,\ddot v + RC\,\dot v + v = v_s(t)$ | 2 | yes | if $v_s$ constant |

> [!example] Worked Example 11.1 — Checking a solution and its interval
> Show that $y = 1/(C - t)$ solves $y' = y^2$, and find the largest interval on which the solution through $y(0) = 2$ exists.
>
> **Step 1 (verify).** $y' = \dfrac{d}{dt}(C - t)^{-1} = (C - t)^{-2} = y^2$. ✓
>
> **Step 2 (fit the condition).** $y(0) = 1/C = 2$ gives $C = 1/2$, so $y = 1/(\tfrac12 - t)$.
>
> **Step 3 (interval).** The formula is undefined at $t = 1/2$. A solution must be differentiable on an *interval* containing $t_0 = 0$, so the solution is defined on $(-\infty, \tfrac12)$. The branch for $t > \tfrac12$ is a *different* solution (it is negative and passes through no point with $t = 0$).
>
> **Check.** At $t = 0.4$: $y = 1/0.1 = 10$, and the numerical slope $[y(0.401) - y(0.399)]/0.002 = (10.101 - 9.901)/0.002 \approx 100 = y^2$. ✓
>
> Notice that $f(t, y) = y^2$ is perfectly smooth everywhere; the blow-up at $t = 1/2$ is invisible in the equation. §11.7 explains why this can happen only in nonlinear equations.

### Direction fields

The equation $y' = f(t, y)$ assigns a slope to every point of the $(t, y)$ plane. Drawing a short segment of slope $f(t, y)$ at each point of a grid gives the **direction field** (slope field). Every solution curve is tangent to the field at every point, so the field is a picture of all solutions at once — available before any formula is found, and even when no formula exists.

Two auxiliary ideas make sketching fast:

- **Isoclines** are the curves $f(t, y) = k$ along which all segments have the same slope $k$. The isocline $k = 0$ (the **nullcline**) is where solutions have horizontal tangents — maxima and minima of $y(t)$.
- **Equilibrium (constant) solutions** are values $y^*$ with $f(t, y^*) = 0$ for all $t$. Their graphs are horizontal lines that other solutions can approach or leave but — when uniqueness holds — never cross.

For example $y' = 1 - y$ has the single equilibrium $y^* = 1$. Above it $y' < 0$ and below it $y' > 0$, so all solutions move towards $y = 1$. The general solution (derived in §11.3) is $y = 1 + Ce^{-t}$:

```desmos-graph
left=-0.5; right=5; bottom=-1.5; top=3.5
---
y=1|black|dashed
y=1+2e^{-x}|blue
y=1+e^{-x}|blue
y=1+0.5e^{-x}|blue
y=1-e^{-x}|red
y=1-2e^{-x}|red
y=1-0.5e^{-x}|red
```
*Integral curves of $y' = 1 - y$ for $C = \pm 0.5, \pm 1, \pm 2$: every curve approaches the stable equilibrium $y = 1$ (dashed). An RC circuit charging or discharging towards its source voltage behaves the same way.*

### Autonomous equations and the phase line

For an autonomous equation $y' = f(y)$ the slope depends only on height, so the direction field is the same on every vertical line. All the information fits on a single $y$-axis — the **phase line** — marked with the equilibria and arrows showing the sign of $f$ between them.

> [!theorem] Stability of an equilibrium of $y' = f(y)$
> Let $f$ be continuously differentiable and $f(y^*) = 0$.
> - If $f'(y^*) < 0$, the equilibrium is **asymptotically stable**: solutions starting close enough to $y^*$ tend to $y^*$ as $t \to \infty$.
> - If $f'(y^*) > 0$, it is **unstable**: nearby solutions move away.
> - If $f'(y^*) = 0$, the test is inconclusive; look at the sign of $f$ on each side.

*Proof sketch (linearisation).* Put $y = y^* + \eta$ with $\eta$ small. Taylor's theorem ([Ch. 4](../part-2-calculus/04-differential-calculus.md)) gives $\eta' = f(y^* + \eta) = f'(y^*)\,\eta + O(\eta^2)$. Dropping the quadratic term, $\eta \approx \eta_0\,e^{f'(y^*)\,t}$, which decays if $f'(y^*) < 0$ and grows if $f'(y^*) > 0$. The quantity $-1/f'(y^*)$ is the *local time constant* of the equilibrium. Making the argument rigorous uses the fact that solutions cannot cross the equilibrium line (uniqueness, §11.7). $\square$

> [!example] Worked Example 11.2 — Phase line of a cubic
> Analyse $y' = y - y^3$ without solving it.
>
> **Step 1 (equilibria).** $f(y) = y(1 - y)(1 + y) = 0$ gives $y^* = -1, 0, 1$.
>
> **Step 2 (linearise).** $f'(y) = 1 - 3y^2$, so $f'(\pm 1) = -2 < 0$ (stable) and $f'(0) = 1 > 0$ (unstable).
>
> **Step 3 (signs).** $f > 0$ on $(-\infty, -1)$ and $(0, 1)$; $f < 0$ on $(-1, 0)$ and $(1, \infty)$. So every solution with $y(0) > 0$ tends to $+1$, every solution with $y(0) < 0$ tends to $-1$, and only $y \equiv 0$ stays at zero.
>
> **Step 4 (time constant).** Near $y^* = 1$, perturbations decay like $e^{-2t}$, so the local time constant is $0.5$. Near $0$ they grow like $e^{t}$.
>
> This is the behaviour of a **bistable latch**: the cross-coupled inverter pair in an SRAM cell has two stable states and an unstable metastable point between them. A latch caught near the middle leaves it at a rate set by $f'(0)$ — the origin of the metastability time constant in synchroniser design.

```desmos-graph
left=-1.8; right=1.8; bottom=-1.2; top=1.2
---
y=x-x^{3}|blue
y=0|black
(-1,0)|label:stable
(0,0)|label:unstable
(1,0)|label:stable
```
*The function $f(y) = y - y^3$ plotted against $y$ (horizontal axis). Where the curve is above the axis solutions move right (increase); where it is below they move left. Zeros with negative slope attract, the zero with positive slope repels.*

### Exact equations

Write a first-order equation in **differential form**

$$
M(t, y)\,dt + N(t, y)\,dy = 0,
$$

which is equivalent to $y' = -M/N$. Suppose there is a function $\Phi(t, y)$ with $\Phi_t = M$ and $\Phi_y = N$. Then along any solution $y(t)$ the chain rule gives

$$
\frac{d}{dt}\Phi\big(t, y(t)\big) = \Phi_t + \Phi_y\,y' = M + N y' = 0,
$$

so $\Phi$ is constant along solutions and the solutions are the **level curves** $\Phi(t, y) = C$.

> [!theorem] Test for exactness
> Let $M, N, M_y, N_t$ be continuous on a rectangle (or any simply connected region). Then $M\,dt + N\,dy$ is the differential of some $\Phi$ — the equation is **exact** — if and only if
> $$
> \frac{\partial M}{\partial y} = \frac{\partial N}{\partial t}.
> $$

Necessity is just equality of mixed partials, $\Phi_{ty} = \Phi_{yt}$. Sufficiency is the construction below. This is exactly the test for a conservative two-dimensional vector field $(M, N)$ in [Ch. 7](../part-2-calculus/07-vector-calculus.md), and $\Phi$ is its potential.

**Construction.** Integrate $M$ in $t$ holding $y$ fixed: $\Phi = \int M\,dt + g(y)$. Differentiate in $y$, set equal to $N$, and solve for $g'(y)$; exactness guarantees that $g'$ comes out independent of $t$.

> [!example] Worked Example 11.3 — An exact equation
> Solve $(2ty + 1)\,dt + (t^2 + 3y^2)\,dy = 0$ with $y(1) = 1$.
>
> **Step 1 (test).** $M_y = 2t$ and $N_t = 2t$. Exact. ✓
>
> **Step 2 (integrate $M$).** $\Phi = \int (2ty + 1)\,dt = t^2 y + t + g(y)$.
>
> **Step 3 (match $N$).** $\Phi_y = t^2 + g'(y) = t^2 + 3y^2$, so $g' = 3y^2$ and $g = y^3$.
>
> **Step 4 (solution).** $t^2 y + t + y^3 = C$. The condition gives $1 + 1 + 1 = 3$, so the solution is the implicit curve $t^2 y + t + y^3 = 3$.
>
> **Check.** $d(t^2 y + t + y^3) = (2ty + 1)\,dt + (t^2 + 3y^2)\,dy$. ✓ Numerically, at $t = 2$ the curve gives $y^3 + 4y - 1 = 0$, whose real root is $y \approx 0.2463$; implicit differentiation predicts slope $-M/N = -(2\cdot2\cdot0.2463 + 1)/(4 + 3\cdot0.0607) \approx -0.474$, which matches a finite-difference slope of the root as $t$ varies.

### Integrating factors for non-exact equations

If $M_y \ne N_t$, multiplying by a suitable nonzero function $\mu(t, y)$ may make $\mu M\,dt + \mu N\,dy = 0$ exact. The exactness condition for $\mu M, \mu N$ is a PDE for $\mu$, usually harder than the original problem, but two special cases are routine.

> [!theorem] Integrating factors depending on one variable
> - If $\dfrac{M_y - N_t}{N}$ depends on $t$ only, call it $R(t)$; then $\mu(t) = e^{\int R(t)\,dt}$.
> - If $\dfrac{N_t - M_y}{M}$ depends on $y$ only, call it $S(y)$; then $\mu(y) = e^{\int S(y)\,dy}$.

*Derivation of the first case.* With $\mu = \mu(t)$, exactness of $\mu M\,dt + \mu N\,dy$ requires $(\mu M)_y = (\mu N)_t$, i.e. $\mu M_y = \mu' N + \mu N_t$, i.e. $\mu'/\mu = (M_y - N_t)/N$. This is consistent only if the right side is a function of $t$ alone, and then $\ln\mu = \int R\,dt$. The second case is symmetric. $\square$

> [!example] Worked Example 11.4 — Making an equation exact
> Solve $(3ty + y^2)\,dt + (t^2 + ty)\,dy = 0$ for $t > 0$.
>
> **Step 1 (test).** $M_y = 3t + 2y$, $N_t = 2t + y$. Not exact.
>
> **Step 2 (look for $\mu(t)$).** $\dfrac{M_y - N_t}{N} = \dfrac{t + y}{t(t + y)} = \dfrac1t$, a function of $t$ only. So $\mu = e^{\int dt/t} = t$.
>
> **Step 3 (exact form).** $(3t^2 y + t y^2)\,dt + (t^3 + t^2 y)\,dy = 0$. Now $M_y = 3t^2 + 2ty = N_t$. ✓
>
> **Step 4 (potential).** $\Phi = \int (3t^2y + ty^2)\,dt = t^3 y + \tfrac12 t^2 y^2 + g(y)$, and $\Phi_y = t^3 + t^2 y + g'(y)$ forces $g' = 0$. The solution is
> $$
> t^3 y + \tfrac12 t^2 y^2 = C.
> $$
>
> **Step 5 (explicit form).** This is quadratic in $y$: $y = -t \pm \sqrt{t^2 + 2C/t^2}$.
>
> **Check.** Take $C = 1/2$ and the $+$ root: at $t = 1$, $y = -1 + \sqrt{2} = 0.4142$. The ODE slope is $-M/N = -(3(0.4142) + 0.1716)/(1 + 0.4142) = -1.0000$; differentiating the explicit formula, $y' = -1 + (t - 2C/t^3)/\sqrt{t^2 + 2C/t^2} = -1 + 0/\sqrt2 = -1$. ✓

> [!warning] Common pitfall — multiplying can add or remove solutions
> An integrating factor that vanishes or blows up somewhere (here $\mu = t$ at $t = 0$) can introduce spurious solutions or hide genuine ones along that curve. State the region ($t > 0$ above) and check the curves where $\mu = 0$ or $\mu = \infty$ separately against the original equation.

> [!tip] Engineering note — conserved quantities
> An exact equation says that some function $\Phi(t, y)$ is conserved along trajectories. In circuits and mechanics this is usually an **energy**: for a lossless LC tank $\tfrac12 L i^2 + \tfrac12 C v^2$ is constant, and writing $dv/di$ from the two element laws gives an exact equation whose level curves are the energy ellipses. Spotting a conserved quantity is often the fastest route to the qualitative behaviour of a nonlinear system.

---

## 11.2 Separable equations

### Motivation

Separable equations are the first-order equations whose solution reduces to *two ordinary integrals*. They cover every autonomous first-order equation $y' = f(y)$ — a capacitor discharging into any nonlinear load, a saturating population, a cooling body — and many time-varying ones. They are also the template for the substitution methods that follow: a clever change of variable often turns an awkward equation into a separable one.

> [!definition] Separable equation
> A first-order ODE is **separable** if it can be written
> $$
> \frac{dy}{dt} = g(t)\,h(y),
> $$
> the product of a function of $t$ alone and a function of $y$ alone.

### The method and why it works

> [!theorem] Solution of a separable equation
> Let $g$ be continuous on an interval $I$ and $h$ continuous with $h \ne 0$ on an interval $J$. Let $H$ be an antiderivative of $1/h$ on $J$ and $G$ an antiderivative of $g$ on $I$. Then $y(t)$ with values in $J$ solves $y' = g(t)h(y)$ if and only if
> $$
> H\big(y(t)\big) = G(t) + C
> $$
> for some constant $C$. Informally, $\displaystyle\int \frac{dy}{h(y)} = \int g(t)\,dt + C$.

*Proof.* If $y$ solves the equation, the chain rule gives $\dfrac{d}{dt}H(y(t)) = H'(y)\,y' = \dfrac{y'}{h(y)} = g(t) = G'(t)$, so $H(y) - G$ is constant. Conversely, differentiating $H(y) = G + C$ reproduces the equation. Because $H' = 1/h$ has one sign on $J$, $H$ is strictly monotone and can in principle be inverted to give $y = H^{-1}(G(t) + C)$ explicitly. $\square$

The informal manipulation "multiply by $dt$, divide by $h(y)$" is therefore legitimate — it is shorthand for the chain-rule argument.

**Procedure.**

1. Find all constant solutions: every root $y^*$ of $h(y^*) = 0$ gives $y \equiv y^*$.
2. On each interval where $h \ne 0$, separate and integrate both sides.
3. Apply the initial condition, solve for $y$ if possible, and state the interval of validity.

> [!warning] Common pitfall — lost solutions
> Dividing by $h(y)$ discards every constant solution $y = y^*$ with $h(y^*) = 0$. They are often the equilibria — the most important solutions physically. Sometimes they are recovered by letting the constant take a special value (often $C = 0$); sometimes they are singular solutions that no value of $C$ reproduces. Always list them in step 1.

> [!example] Worked Example 11.5 — A Gaussian solution
> Solve $y' = -2ty$ with $y(0) = 3$.
>
> **Step 1.** $h(y) = y$ vanishes at $y = 0$, so $y \equiv 0$ is a solution.
>
> **Step 2.** For $y \ne 0$: $\dfrac{dy}{y} = -2t\,dt$, so $\ln\lvert y\rvert = -t^2 + C_1$, i.e. $\lvert y\rvert = e^{C_1}e^{-t^2}$. Removing the absolute value, $y = Ce^{-t^2}$ with $C = \pm e^{C_1} \ne 0$; allowing $C = 0$ restores the lost solution.
>
> **Step 3.** $y(0) = C = 3$, so $y = 3e^{-t^2}$, valid for all $t$.
>
> **Check.** $y' = -6te^{-t^2} = -2t\,(3e^{-t^2})$. ✓ At $t = 1$, $y = 3e^{-1} = 1.1036$.

> [!example] Worked Example 11.6 — Capacitor discharging through a square-law device
> A capacitor $C$ charged to $V_0$ discharges through a device whose current is $i = k v^2$ (a crude model of a space-charge-limited or square-law load). KCL at the capacitor node gives $C\,\dot v = -k v^2$, $v(0) = V_0 > 0$. Find $v(t)$ and compare with a resistor of the same initial current.
>
> **Step 1.** $v \equiv 0$ is a solution; since $V_0 > 0$ and solutions cannot cross it, $v > 0$ for all $t$.
>
> **Step 2.** $\dfrac{dv}{v^2} = -\dfrac{k}{C}\,dt$, so $-\dfrac1v = -\dfrac{k}{C}t + C_1$, and $v(0) = V_0$ gives $C_1 = -1/V_0$.
>
> **Step 3.** Solving for $v$:
> $$
> v(t) = \frac{V_0}{1 + t/T},\qquad T = \frac{C}{kV_0}.
> $$
>
> **Units.** $k$ is in $\mathrm{A/V^2}$ and $C$ in $\mathrm{A\,s/V}$, so $T = C/(kV_0)$ has units $\mathrm{(A\,s/V)/(A/V^2\cdot V)} = \mathrm{s}$. ✓
>
> **Comparison.** A resistor drawing the same initial current $kV_0^2$ has $R = 1/(kV_0)$ and gives $v = V_0e^{-t/RC} = V_0e^{-t/T}$, the same initial slope. At $t = 5T$ the resistor has left $V_0e^{-5} = 0.0067\,V_0$ while the square-law device has left $V_0/6 = 0.167\,V_0$. The nonlinear element changes the *form* of the transient — an algebraic $1/t$ tail instead of an exponential — not just its speed.

```desmos-graph
left=-0.2; right=6; bottom=-0.05; top=1.1
---
y=\frac{1}{1+x}|blue|x>0
y=e^{-x}|red|dashed|x>0
```
*Normalised discharge $v/V_0$ against $t/T$: square-law load (blue, $1/(1 + t/T)$) and a resistor with the same initial current (red dashed, $e^{-t/T}$). Same initial slope, very different tails.*

> [!example] Worked Example 11.7 — Capacitor discharging through a diode
> A capacitor $C = 1\ \mathrm{nF}$ charged to $V_0 = 0.7\ \mathrm V$ discharges through an ideal-exponential diode, $i = I_S\big(e^{v/V_T} - 1\big)$ with $I_S = 10^{-14}\ \mathrm A$ and $V_T = 25\ \mathrm{mV}$. Find $v(t)$ and the times at which $v$ reaches $0.6$, $0.5$, $0.4$ and $0.3$ V.
>
> **Step 1 (set up).** KCL: $C\,\dot v = -I_S\big(e^{v/V_T} - 1\big)$. Autonomous, hence separable; $v \equiv 0$ is the only equilibrium.
>
> **Step 2 (integrate).** Multiply numerator and denominator by $e^{-v/V_T}$:
> $$
> \int \frac{dv}{e^{v/V_T} - 1} = \int \frac{e^{-v/V_T}\,dv}{1 - e^{-v/V_T}} = V_T\ln\big(1 - e^{-v/V_T}\big).
> $$
> So $V_T\ln\big(1 - e^{-v/V_T}\big) = -\dfrac{I_S}{C}t + C_1$.
>
> **Step 3 (solve).** With $\tau = CV_T/I_S$ and the initial condition,
> $$
> 1 - e^{-v/V_T} = \big(1 - e^{-V_0/V_T}\big)e^{-t/\tau},\qquad v(t) = -V_T\ln\!\Big[1 - \big(1 - e^{-V_0/V_T}\big)e^{-t/\tau}\Big].
> $$
>
> **Step 4 (numbers).** $\tau = 10^{-9}\cdot 0.025/10^{-14} = 2500\ \mathrm s$. Inverting, $t = -\tau\ln\big[(1 - e^{-v/V_T})/(1 - e^{-V_0/V_T})\big]$:
>
> | $v$ (V) | 0.6 | 0.5 | 0.4 | 0.3 |
> |---|---|---|---|---|
> | $t$ | 92.6 ns | 5.15 µs | 0.281 ms | 15.4 ms |
>
> **Interpretation.** Each 0.1 V drop takes about $e^{0.1/V_T} = e^4 \approx 55$ times longer than the previous one. While $v \gg V_T$ the solution is approximately $v \approx -V_T\ln\big(e^{-V_0/V_T} + t/\tau\big)$ — a *logarithmic* decay. The huge $\tau$ governs only the final approach once $v \lesssim V_T$.
>
> **Check.** A stiff numerical integration (Radau, relative tolerance $10^{-10}$) of the original ODE from $0.7$ V for $15.36$ ms ends at $0.30000$ V. ✓

### Equations reducible to separable form

Two substitutions cover most textbook cases.

**(a) $y' = F(y/t)$ — "homogeneous of degree zero".** The right side is unchanged when $t$ and $y$ are scaled together. Put $v = y/t$, so $y = vt$ and $y' = v + tv'$. The equation becomes

$$
t\,\frac{dv}{dt} = F(v) - v,
$$

which is separable in $v$ and $t$.

**(b) $y' = G(at + by + c)$.** Put $u = at + by + c$. Then $u' = a + b\,G(u)$, which is autonomous and therefore separable. For example $y' = (t + y)^2$ becomes $u' = 1 + u^2$ with $u = t + y$, so $u = \tan(t + C)$ and $y = \tan(t + C) - t$.

> [!example] Worked Example 11.8 — A scale-invariant equation
> Solve $y' = \dfrac{t^2 + y^2}{ty}$ with $y(1) = 2$, for $t > 0$.
>
> **Step 1 (recognise).** Dividing top and bottom by $t^2$, $y' = \dfrac{1 + (y/t)^2}{y/t} = F(y/t)$ with $F(v) = (1 + v^2)/v$.
>
> **Step 2 (substitute).** $tv' = \dfrac{1 + v^2}{v} - v = \dfrac1v$, so $v\,dv = \dfrac{dt}{t}$.
>
> **Step 3 (integrate).** $\tfrac12 v^2 = \ln t + C$, i.e. $y^2 = t^2(2\ln t + 2C)$.
>
> **Step 4 (initial condition).** $y(1) = 2$ gives $4 = 2C$, so $y = t\sqrt{2\ln t + 4}$ (positive root, since $y(1) > 0$).
>
> **Step 5 (interval).** The square root needs $2\ln t + 4 > 0$, i.e. $t > e^{-2} \approx 0.135$. At $t = e^{-2}$ the solution reaches $y = 0$, where the equation itself is singular. The solution lives on $(e^{-2}, \infty)$.
>
> **Check.** At $t = e$: $y = e\sqrt6 = 6.658$. The ODE gives $y' = (e^2 + 6e^2)/(e\cdot e\sqrt6) = 7/\sqrt6 = 2.858$; differentiating the formula, $y' = \sqrt{2\ln t + 4} + 1/\sqrt{2\ln t + 4} = \sqrt6 + 1/\sqrt6 = 2.858$. ✓

### The logistic equation

The logistic equation $y' = ry(1 - y/K)$, with growth rate $r > 0$ and capacity $K > 0$, is the standard model of growth that saturates. Its equilibria are $y = 0$ (unstable, $f'(0) = r$) and $y = K$ (stable, $f'(K) = -r$). For $0 < y < K$, partial fractions ([Ch. 5](../part-2-calculus/05-integral-calculus.md)) give

$$
\int \frac{K\,dy}{y(K - y)} = \int\left(\frac1y + \frac{1}{K - y}\right)dy = \ln\frac{y}{K - y} = rt + C,
$$

and solving with $y(0) = y_0$,

$$
y(t) = \frac{K}{1 + \left(\dfrac{K}{y_0} - 1\right)e^{-rt}}.
$$

The curve is an S-shape with its steepest point (inflection) at $y = K/2$, where $y'' = f'(y)\,y' = 0$ because $f'(y) = r(1 - 2y/K)$ vanishes there. Logistic behaviour appears in soft-limiting amplifiers, in the build-up of oscillation amplitude in a self-starting oscillator, and in technology-adoption curves ([Ch. 25](../part-7-numerical-computational/25-mathematical-modeling.md)).

```desmos-graph
left=-1; right=10; bottom=-0.1; top=2.2
---
y=\frac{1}{1+9e^{-x}}|blue
y=\frac{1}{1-0.5e^{-x}}|red|x>0
y=1|black|dashed
y=0.5|grey|dotted
```
*Logistic solutions with $r = 1$, $K = 1$: from $y_0 = 0.1$ (blue, S-shaped, steepest at $y = K/2$) and from $y_0 = 2$ (red, decaying to $K$). Both approach the stable equilibrium $y = K$.*

> [!tip] Engineering note — when separation is the right tool
> Any single-state circuit with constant sources is autonomous: $C\,\dot v = I - g(v)$. It is always separable, so the time to reach a voltage, $t(v) = \int_{V_0}^{v} C\,dv/\big(I - g(v)\big)$, can be written down even when the integral has no closed form — evaluate it numerically. This is how the delay of a nonlinear RC timing stage, or the droop time of a hold capacitor through a leakage diode, is computed in practice.

---

## 11.3 Linear first-order equations

### Motivation

Linear first-order equations are the one large class of first-order equations that can *always* be solved in closed form (up to evaluating integrals). They describe every single-storage linear circuit — including those with time-varying parameters such as a switched-capacitor stage, a varactor, a thermistor-dependent resistance or a modulated source — and their solution formula already displays the structure of all linear systems: response = (decaying effect of the initial state) + (weighted memory of the input).

> [!definition] Linear first-order ODE
> The standard form is
> $$
> y' + p(t)\,y = r(t),
> $$
> with $p$ and $r$ continuous on an interval $I$. It is **homogeneous** if $r \equiv 0$ and **nonhomogeneous** otherwise. An equation $a_1(t)y' + a_0(t)y = b(t)$ is put in standard form by dividing by $a_1$, which is allowed only where $a_1 \ne 0$; points where $a_1 = 0$ are **singular points** of the equation.

### The homogeneous equation

$y' + p(t)y = 0$ is separable: $dy/y = -p\,dt$, so

$$
y_h(t) = C\,e^{-P(t)},\qquad P(t) = \int p(t)\,dt.
$$

Because $e^{-P}$ is never zero, the homogeneous solution never changes sign: a solution that starts positive stays positive. For constant $p = 1/\tau$ this is the familiar $Ce^{-t/\tau}$.

### Integrating factor

> [!theorem] Solution of the linear first-order equation
> Let $p, r$ be continuous on $I$ and $\mu(t) = e^{P(t)}$ with $P' = p$. Then every solution of $y' + p y = r$ on $I$ is
> $$
> y(t) = \frac{1}{\mu(t)}\left[\int \mu(t)\,r(t)\,dt + C\right].
> $$
> With the initial condition $y(t_0) = y_0$ and $P(t) = \int_{t_0}^{t} p(s)\,ds$, the unique solution is
> $$
> y(t) = y_0\,e^{-P(t)} + \int_{t_0}^{t} e^{-[P(t) - P(s)]}\,r(s)\,ds,
> $$
> and it exists on the **whole** interval $I$.

*Derivation.* The integrating factor is chosen so that the left side becomes an exact derivative. By the product rule, $(\mu y)' = \mu y' + \mu' y$. If $\mu' = p\mu$ — that is, $\mu = e^{P}$ — then $(\mu y)' = \mu(y' + py) = \mu r$. Integrating, $\mu y = \int \mu r\,dt + C$. For the definite form integrate from $t_0$ to $t$ with $\mu(t_0) = 1$: $\mu(t)y(t) - y_0 = \int_{t_0}^{t}\mu(s)r(s)\,ds$, then divide by $\mu(t) = e^{P(t)}$. Every step is reversible, so there are no other solutions, and the integrals exist wherever $p$ and $r$ are continuous. $\square$

(It is also the special case of §11.1's integrating factor: in differential form $(py - r)\,dt + dy = 0$, the test quantity $(M_y - N_t)/N = p(t)$ depends on $t$ alone.)

The definite form deserves a careful reading:

- $y_0e^{-P(t)}$ is the **zero-input response**: what the initial state does with no input.
- The integral is the **zero-state response**: each input value $r(s)$ at an earlier time $s$ contributes $r(s)\,ds$, weighted by the factor $e^{-[P(t) - P(s)]}$ by which the system has "forgotten" it by time $t$.
- For constant $p = 1/\tau$ the weight is $e^{-(t - s)/\tau}$, which depends only on the elapsed time $t - s$, and the integral is a **convolution** of $r$ with $e^{-t/\tau}$. This is the first appearance of the impulse response; [Ch. 12](12-linear-differential-equations.md) and [Ch. 14](../part-5-transforms/14-laplace-transforms.md) develop it.

Whatever the form, the general solution is $y = Ce^{-P} + y_p$: a homogeneous solution plus *one* particular solution. This is the first instance of the structure theorem of §11.6.

> [!example] Worked Example 11.9 — Exponential forcing
> Solve $y' + 2y = e^{-t}$, $y(0) = 0$.
>
> **Step 1.** $p = 2$, so $\mu = e^{2t}$.
>
> **Step 2.** $(e^{2t}y)' = e^{2t}e^{-t} = e^{t}$, hence $e^{2t}y = e^{t} + C$ and $y = e^{-t} + Ce^{-2t}$.
>
> **Step 3.** $y(0) = 1 + C = 0$, so $C = -1$:
> $$
> y(t) = e^{-t} - e^{-2t}.
> $$
>
> **Check.** $y' = -e^{-t} + 2e^{-2t}$ and $y' + 2y = -e^{-t} + 2e^{-2t} + 2e^{-t} - 2e^{-2t} = e^{-t}$. ✓
>
> **Shape.** The response rises from zero, peaks where $y' = 0$, i.e. $e^{-t} = 2e^{-2t}$, so $t = \ln 2 \approx 0.693$ with $y = \tfrac12 - \tfrac14 = \tfrac14$, then decays like the slower exponential $e^{-t}$.

```desmos-graph
left=-0.3; right=6; bottom=-0.05; top=0.35
---
y=e^{-x}-e^{-2x}|blue|x>0
(0.693,0.25)|label:peak (ln 2, 1/4)
```
*The solution $y = e^{-t} - e^{-2t}$ of Example 11.9: a two-time-constant pulse, typical of a first-order stage driven by the decaying output of another stage.*

> [!example] Worked Example 11.10 — Variable coefficient and a singular point
> Solve $t\,y' + 2y = 4t^2$ with $y(1) = 2$.
>
> **Step 1 (standard form).** Dividing by $t$ (allowed for $t \ne 0$): $y' + \dfrac2t y = 4t$. The coefficient $p = 2/t$ is continuous on $(0, \infty)$, which contains $t_0 = 1$, so the solution exists and is unique on $(0, \infty)$.
>
> **Step 2.** $\mu = e^{\int 2/t\,dt} = e^{2\ln t} = t^2$.
>
> **Step 3.** $(t^2y)' = 4t^3$, so $t^2y = t^4 + C$ and $y = t^2 + C/t^2$.
>
> **Step 4.** $y(1) = 1 + C = 2$ gives $C = 1$: $y = t^2 + t^{-2}$.
>
> **Check.** $t\,y' + 2y = t(2t - 2t^{-3}) + 2t^2 + 2t^{-2} = 4t^2$. ✓
>
> **Remark.** Every solution except the one with $C = 0$ blows up at the singular point $t = 0$. The theorem promises existence only on $(0,\infty)$; the coefficient's singularity, not the solution formula, sets the interval.

> [!example] Worked Example 11.11 — Tracking a ramp
> An RC low-pass stage with time constant $\tau$ is driven by a ramp: $\tau\dot v + v = kt$, $v(0) = 0$. Find $v(t)$ and the steady tracking error.
>
> **Step 1.** Standard form $\dot v + v/\tau = kt/\tau$, $\mu = e^{t/\tau}$.
>
> **Step 2.** $\big(e^{t/\tau}v\big)' = \dfrac{k}{\tau}te^{t/\tau}$. Integrating by parts, $\int te^{t/\tau}\,dt = \tau te^{t/\tau} - \tau^2e^{t/\tau}$, so $e^{t/\tau}v = k(t - \tau)e^{t/\tau} + C$.
>
> **Step 3.** $v(0) = -k\tau + C = 0$ gives $C = k\tau$:
> $$
> v(t) = k(t - \tau) + k\tau\,e^{-t/\tau}.
> $$
>
> **Check.** $\tau\dot v + v = \tau(k - ke^{-t/\tau}) + kt - k\tau + k\tau e^{-t/\tau} = kt$. ✓
>
> **Interpretation.** After a few time constants the output is a ramp of the same slope, delayed by exactly $\tau$: the tracking error $kt - v \to k\tau$. With $k = 1\ \mathrm{V/ms}$ and $\tau = 0.2\ \mathrm{ms}$ the error settles at $0.2\ \mathrm V$; at $t = 1\ \mathrm{ms}$ it is $0.2(1 - e^{-5}) = 0.1987\ \mathrm V$. This lag is why a first-order sample-and-hold or tracking loop always trails a ramp.

> [!example] Worked Example 11.12 — A time-varying capacitor
> A varactor whose capacitance grows linearly, $C(t) = C_0(1 + t/T)$, is charged from $0$ V through a resistor $R$ by a constant source $V_s$. Take $RC_0 = T$. Find $v(t)$.
>
> **Step 1 (model).** For a time-varying capacitor the charge is $q = C(t)v$, and the current is $i = \dot q = C\dot v + \dot C v$, not just $C\dot v$. KVL gives $Ri + v = V_s$:
> $$
> RC_0\Big(1 + \frac tT\Big)\dot v + \Big(\frac{RC_0}{T} + 1\Big)v = V_s.
> $$
>
> **Step 2 (dimensionless time).** With $s = t/T$ and $RC_0 = T$, this becomes $(1 + s)\,v' + 2v = V_s$, where $' = d/ds$.
>
> **Step 3 (integrating factor).** Standard form $v' + \dfrac{2}{1 + s}v = \dfrac{V_s}{1 + s}$, so $\mu = (1 + s)^2$ and $\big((1 + s)^2v\big)' = (1 + s)V_s$.
>
> **Step 4.** $(1 + s)^2v = \tfrac12 V_s(1 + s)^2 + K$; $v(0) = 0$ gives $K = -V_s/2$:
> $$
> v(t) = \frac{V_s}{2}\left[1 - \frac{1}{(1 + t/T)^2}\right].
> $$
>
> **Interpretation.** The voltage approaches $V_s/2$, not $V_s$, and algebraically rather than exponentially. Physically, holding $v$ constant on a growing capacitor still needs a current $\dot C v = C_0 v/T$; at $v = V_s/2$ the resistor drop $R C_0 v/T = v$ supplies exactly that. Time-varying elements break the intuition "a capacitor eventually charges to the source" — which is the principle behind parametric amplifiers and charge pumps.
>
> **Check.** $(1 + s)v' + 2v = (1 + s)\cdot V_s(1 + s)^{-3} + V_s - V_s(1 + s)^{-2} = V_s$. ✓

### Bernoulli equations

> [!definition] Bernoulli equation
> $$
> y' + p(t)\,y = q(t)\,y^n,\qquad n \ne 0, 1.
> $$
> (For $n = 0$ or $1$ it is already linear.)

Divide by $y^n$ and put $v = y^{1-n}$, so $v' = (1 - n)y^{-n}y'$:

$$
y^{-n}y' + p\,y^{1-n} = q \quad\Longrightarrow\quad v' + (1 - n)\,p(t)\,v = (1 - n)\,q(t),
$$

which is linear in $v$. If $n > 0$, $y \equiv 0$ is an extra solution lost by the division.

> [!example] Worked Example 11.13 — Amplitude of a self-starting oscillator
> Near its operating point the slowly varying amplitude $a(t)$ of many oscillators (van der Pol type) obeys $\dot a = \alpha a - \beta a^3$: linear growth from the negative-resistance element, cubic limiting from gain compression. Take $\alpha = \beta = 1$ and $a(0) = 1/2$.
>
> **Step 1 (recognise).** $\dot a - a = -a^3$ is Bernoulli with $n = 3$. Put $v = a^{-2}$: $\dot v + (1 - 3)(-1)v = (1 - 3)(-1)$, i.e. $\dot v + 2v = 2$.
>
> **Step 2 (solve linear).** $v = 1 + Ce^{-2t}$; $v(0) = a(0)^{-2} = 4$ gives $C = 3$.
>
> **Step 3 (back-substitute).**
> $$
> a(t) = \frac{1}{\sqrt{1 + 3e^{-2t}}}.
> $$
>
> **Check.** As $t \to \infty$, $a \to 1 = \sqrt{\alpha/\beta}$, the stable equilibrium of $\dot a = a - a^3$ (Example 11.2). At $t = 1$: $a = 1/\sqrt{1 + 3e^{-2}} = 0.8433$. Direct differentiation confirms $\dot a = a - a^3$ symbolically. ✓
>
> The amplitude grows like $e^{\alpha t}$ while small and saturates at $\sqrt{\alpha/\beta}$: the oscillator's steady amplitude is set by the ratio of small-signal gain to compression.

### Riccati equations

> [!definition] Riccati equation
> $$
> y' = q_0(t) + q_1(t)\,y + q_2(t)\,y^2.
> $$

Riccati equations generally have no elementary solution, but **if one particular solution $y_1$ is known**, the substitution $y = y_1 + 1/v$ turns the equation into a linear one for $v$:

$$
v' = -\big(q_1 + 2q_2y_1\big)v - q_2.
$$

(Substitute, use the fact that $y_1$ satisfies the equation, and multiply by $-v^2$.) Matrix Riccati equations are central in optimal control (the LQR gain) and Kalman filtering, and the scalar case is their prototype.

> [!example] Worked Example 11.14 — Riccati with a known solution
> Solve $y' = 1 + t^2 - 2ty + y^2$, given that $y_1 = t$ is a solution.
>
> **Step 1 (verify $y_1$).** $y_1' = 1$ and $1 + t^2 - 2t^2 + t^2 = 1$. ✓
>
> **Step 2 (substitute).** Here $q_2 = 1$, $q_1 = -2t$, so $q_1 + 2q_2y_1 = 0$ and the linear equation is $v' = -1$: $v = C - t$.
>
> **Step 3.** $y = t + \dfrac{1}{C - t}$, plus the solution $y = t$ itself (the limit $C \to \infty$).
>
> **Check.** $y' = 1 + (C - t)^{-2}$, and $1 + t^2 - 2ty + y^2 = 1 + (y - t)^2 = 1 + (C - t)^{-2}$. ✓ Each solution other than $y_1$ blows up at $t = C$ — another nonlinear finite escape.

> [!tip] Engineering note — the time constant
> Any first-order linear time-invariant system — an RC or RL network, a thermal mass, a single-pole amplifier, a first-order sensor — reduces to $\tau\dot y + y = K\,u(t)$. The homogeneous part $e^{-t/\tau}$ falls to 36.8% in one time constant, 13.5% in two, 1.8% in four and under 1% in five. The rule "settled after $5\tau$" comes from this table. For time-varying coefficients the role of $t/\tau$ is played by $P(t) = \int p\,dt$, so the system "forgets" its initial state when $P$ has grown by about 5.

> [!warning] Common pitfall — forgetting the standard form
> The integrating factor is $e^{\int p\,dt}$ where $p$ is the coefficient of $y$ **after** the coefficient of $y'$ has been made 1. For $t\,y' + 2y = 4t^2$ the factor is $t^2$, not $e^{2t}$. The same mistake in variation of parameters (§11.6) is the commonest error in second-order work.

---

## 11.4 Second-order ODEs

### Motivation

Two independent energy stores — an inductor and a capacitor, a mass and a spring, an electric and a magnetic field — produce a second-order equation. Second order is the lowest order at which **oscillation** is possible, which is why the series RLC circuit, the tuned amplifier and the mechanical resonator are all second-order systems. It is also the order of the field equations of physics in one dimension (Poisson's equation, the transmission-line equations, Schrödinger's equation), which is where boundary-value problems (§11.8) come from.

### Forms of a second-order equation

A second-order ODE is $F(t, y, y', y'') = 0$, or in normal form $y'' = f(t, y, y')$. A solution is fixed by **two** conditions — typically $y(t_0)$ and $y'(t_0)$ — because two integrations are needed. The most important class is the **linear** equation

$$
a_2(t)\,y'' + a_1(t)\,y' + a_0(t)\,y = b(t),
$$

which on any interval where $a_2 \ne 0$ can be divided through to give the **standard form**

$$
y'' + p(t)\,y' + q(t)\,y = r(t).
$$

It is convenient to name the left side as an operator: $L[y] = y'' + py' + qy$. The operator is **linear**: $L[c_1y_1 + c_2y_2] = c_1L[y_1] + c_2L[y_2]$, because differentiation and multiplication by a function are linear.

> [!theorem] Existence and uniqueness for linear second-order IVPs
> If $p, q, r$ are continuous on an open interval $I$ and $t_0 \in I$, then for any numbers $y_0, y_1$ there is exactly one solution of $y'' + py' + qy = r$ with $y(t_0) = y_0$, $y'(t_0) = y_1$, and it exists on all of $I$.

This is proved in §11.7 (by rewriting the equation as a first-order system). Everything in the structure theory below rests on it.

### Superposition and the dimension of the solution space

> [!theorem] Superposition and dimension
> If $y_1$ and $y_2$ solve the homogeneous equation $L[y] = 0$ on $I$, so does every linear combination $c_1y_1 + c_2y_2$. The set of all solutions is a vector space ([Ch. 9](../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md)) of dimension exactly 2.

*Proof.* Closure under combinations is linearity of $L$. For the dimension, consider the map that sends a solution $y$ to its initial data $\big(y(t_0), y'(t_0)\big) \in \mathbb{R}^2$. It is linear; it is onto by the existence theorem (every pair of initial data is achieved); and it is one-to-one by uniqueness (only $y \equiv 0$ has zero initial data). A linear bijection onto $\mathbb{R}^2$ means the solution space is isomorphic to $\mathbb{R}^2$. $\square$

So the general solution of $L[y] = 0$ is $y = c_1y_1 + c_2y_2$ for **any** two linearly independent solutions — a **fundamental set**, i.e. a basis.

### Linear independence and the Wronskian

Two functions are linearly dependent on $I$ if one is a constant multiple of the other there. For solutions this is tested with the **Wronskian**

$$
W(y_1, y_2)(t) = \det\begin{pmatrix} y_1 & y_2 \\ y_1' & y_2' \end{pmatrix} = y_1y_2' - y_1'y_2.
$$

The determinant is exactly the one that appears when we try to fit initial conditions: solving $c_1y_1(t_0) + c_2y_2(t_0) = y_0$, $c_1y_1'(t_0) + c_2y_2'(t_0) = y_1$ for $c_1, c_2$ is possible for *every* right side precisely when this determinant is nonzero ([Ch. 8](../part-3-linear-algebra/08-matrices-and-linear-systems.md)).

> [!theorem] Abel's formula and the Wronskian test
> Let $y_1, y_2$ solve $y'' + p(t)y' + q(t)y = 0$ on $I$, with $p, q$ continuous. Then
> $$
> W(t) = W(t_0)\,\exp\!\left(-\int_{t_0}^{t} p(s)\,ds\right).
> $$
> Consequently $W$ is either **never zero** on $I$ (and $y_1, y_2$ form a fundamental set) or **identically zero** (and they are dependent).

*Proof.* Differentiate: $W' = y_1y_2'' - y_1''y_2$ (the $y_1'y_2'$ terms cancel). Substituting $y_k'' = -py_k' - qy_k$,

$$
W' = y_1(-py_2' - qy_2) - (-py_1' - qy_1)y_2 = -p\,(y_1y_2' - y_1'y_2) = -p\,W.
$$

This is a first-order linear equation for $W$ whose solution is Abel's formula (§11.3). The exponential is never zero, so $W(t)$ vanishes everywhere or nowhere. If $W \equiv 0$, pick $t_0$; the system for $(c_1, c_2)$ with zero right side has a nontrivial solution, and $c_1y_1 + c_2y_2$ is then a solution with zero initial data, hence identically zero by uniqueness — the functions are dependent. $\square$

Abel's formula is useful even without knowing the solutions: it gives $W$ directly from the coefficient $p$, which is the key ingredient in reduction of order and variation of parameters.

> [!warning] Common pitfall — the Wronskian test needs solutions
> For arbitrary functions, $W \equiv 0$ does **not** imply dependence. On $(-1, 1)$, $f = t^2$ and $g = t\lvert t\rvert$ are independent (neither is a constant multiple of the other on the whole interval), yet $W(f, g) = 0$ for every $t$. The "never or always zero" dichotomy is a property of *solutions of one linear equation with continuous coefficients*; it tells you that $t^2$ and $t\lvert t\rvert$ cannot both solve such an equation on an interval containing $0$.

### Reduction of order

If one nonzero solution $y_1$ of $y'' + py' + qy = 0$ is known, a second one can always be found.

> [!theorem] Reduction of order
> $$
> y_2(t) = y_1(t)\int \frac{e^{-\int p\,dt}}{y_1(t)^2}\,dt
> $$
> is a second solution, independent of $y_1$, on any interval where $y_1 \ne 0$.

*Derivation.* Seek $y_2 = v(t)\,y_1(t)$. Then $y_2' = v'y_1 + vy_1'$ and $y_2'' = v''y_1 + 2v'y_1' + vy_1''$. Substituting,

$$
L[vy_1] = v''y_1 + v'(2y_1' + py_1) + v\underbrace{(y_1'' + py_1' + qy_1)}_{=0} = 0.
$$

The terms in $v$ itself cancel — this is why the method works — leaving a *first-order* linear equation for $w = v'$: $w' + (2y_1'/y_1 + p)w = 0$. Its integrating factor is $y_1^2e^{\int p\,dt}$, so $w = e^{-\int p\,dt}/y_1^2$, and $v = \int w\,dt$. Since $v$ is not constant, $y_2/y_1$ is not constant and the solutions are independent. (Alternatively: Abel's formula gives $W = y_1y_2' - y_1'y_2 = y_1^2(y_2/y_1)' = e^{-\int p}$, which is the same equation.) $\square$

> [!example] Worked Example 11.15 — Reduction of order for a Legendre equation
> $y_1 = t$ solves $(1 - t^2)y'' - 2ty' + 2y = 0$ on $-1 < t < 1$ (Legendre's equation with $n = 1$). Find a second solution.
>
> **Step 1 (verify $y_1$).** $y_1'' = 0$, and $-2t\cdot1 + 2t = 0$. ✓
>
> **Step 2 (standard form).** $p = -\dfrac{2t}{1 - t^2}$, so $-\int p\,dt = \int\dfrac{2t}{1 - t^2}\,dt = -\ln(1 - t^2)$ and $e^{-\int p\,dt} = \dfrac{1}{1 - t^2}$.
>
> **Step 3 (integrate).** By partial fractions,
> $$
> \frac{1}{t^2(1 - t^2)} = \frac{1}{t^2} + \frac{1}{1 - t^2},\qquad \int\left(\frac{1}{t^2} + \frac{1}{1 - t^2}\right)dt = -\frac1t + \frac12\ln\frac{1 + t}{1 - t}.
> $$
>
> **Step 4.**
> $$
> y_2 = t\left(-\frac1t + \frac12\ln\frac{1 + t}{1 - t}\right) = \frac t2\ln\frac{1 + t}{1 - t} - 1.
> $$
> Although the integration was done on intervals avoiding $t = 0$, the result is smooth across $0$ and solves the equation on all of $(-1, 1)$; it is the Legendre function of the second kind $Q_1$.
>
> **Check.** Symbolic substitution gives $(1 - t^2)y_2'' - 2ty_2' + 2y_2 = 0$. Numerically at $t = 0.5$: $y_2 = 0.25\ln3 - 1 = -0.72535$. Abel predicts $W = C/(1 - t^2)$, and indeed $W(t, y_2) = ty_2' - y_2 = 1/(1 - t^2)$, which at $t = 0.5$ is $1.3333$. ✓
>
> $y_2$ blows up logarithmically at $t = \pm1$, the singular points where $a_2 = 1 - t^2$ vanishes. This is why only the polynomial solution $P_1 = t$ is admissible in physical problems on a sphere, where $t = \cos\theta$ reaches $\pm 1$ at the poles ([Ch. 13](13-partial-differential-equations.md)).

### Cauchy–Euler (equidimensional) equations

The most important variable-coefficient equation with elementary solutions is

$$
t^2y'' + a\,t\,y' + b\,y = 0,\qquad t > 0,
$$

with constants $a, b$. Each term scales the same way under $t \to \lambda t$ (each derivative brings a factor $1/t$ that the power of $t$ cancels), so a power $y = t^m$ should work. Substituting, $t^2m(m - 1)t^{m-2} + a\,t\,mt^{m-1} + bt^m = t^m[m(m - 1) + am + b]$, giving the **indicial equation**

$$
m^2 + (a - 1)m + b = 0.
$$

**Why it works — the substitution $t = e^x$.** Put $x = \ln t$. Then $\dfrac{dy}{dt} = \dfrac1t\dfrac{dy}{dx}$ and $\dfrac{d^2y}{dt^2} = \dfrac{1}{t^2}\left(\dfrac{d^2y}{dx^2} - \dfrac{dy}{dx}\right)$, so $ty' = \dot y$ and $t^2y'' = \ddot y - \dot y$ (dots meaning $d/dx$). The equation becomes the constant-coefficient equation

$$
\ddot y + (a - 1)\dot y + by = 0,
$$

whose characteristic equation is exactly the indicial equation, and $e^{mx} = t^m$. The three cases of a constant-coefficient equation ([Ch. 12](12-linear-differential-equations.md)) translate directly:

| Roots of $m^2 + (a - 1)m + b = 0$ | Fundamental set on $t > 0$ |
|---|---|
| real, distinct $m_1 \ne m_2$ | $t^{m_1},\ t^{m_2}$ |
| repeated $m$ | $t^{m},\ t^{m}\ln t$ |
| complex $\alpha \pm j\beta$ | $t^{\alpha}\cos(\beta\ln t),\ t^{\alpha}\sin(\beta\ln t)$ |

For $t < 0$ replace $t$ by $\lvert t\rvert$. The point $t = 0$ is a singular point, and solutions typically blow up or lose smoothness there.

> [!example] Worked Example 11.16 — Real roots, cross-checked by reduction of order
> Solve $t^2y'' - 2ty' + 2y = 0$ for $t > 0$, and check with Abel's formula.
>
> **Step 1.** $a = -2$, $b = 2$: $m^2 - 3m + 2 = (m - 1)(m - 2) = 0$, so $y = c_1t + c_2t^2$.
>
> **Step 2 (reduction from $y_1 = t$).** Standard form has $p = -2/t$, so $e^{-\int p\,dt} = t^2$ and $y_2 = t\int t^2/t^2\,dt = t^2$. ✓
>
> **Step 3 (Wronskian).** $W = t\cdot2t - 1\cdot t^2 = t^2 \ne 0$. Abel: $W \propto e^{\int 2/t\,dt} = t^2$. ✓

> [!example] Worked Example 11.17 — Complex indicial roots
> Solve $t^2y'' + ty' + 4y = 0$ with $y(1) = 1$, $y'(1) = 0$.
>
> **Step 1.** $a = 1$, $b = 4$: $m^2 + 4 = 0$, so $m = \pm 2j$ ($\alpha = 0$, $\beta = 2$).
>
> **Step 2.** $y = c_1\cos(2\ln t) + c_2\sin(2\ln t)$, and $y' = \dfrac{2}{t}\big[-c_1\sin(2\ln t) + c_2\cos(2\ln t)\big]$.
>
> **Step 3.** At $t = 1$ ($\ln t = 0$): $c_1 = 1$ and $2c_2 = 0$. So $y = \cos(2\ln t)$.
>
> **Check.** $y' = -\dfrac2t\sin(2\ln t)$, $y'' = \dfrac{2}{t^2}\sin(2\ln t) - \dfrac{4}{t^2}\cos(2\ln t)$. Then $t^2y'' + ty' + 4y = 2\sin - 4\cos - 2\sin + 4\cos = 0$. ✓
>
> **Shape.** The solution oscillates in $\ln t$, not in $t$: its zeros are at $t = e^{\pi/4} = 2.193$, $e^{3\pi/4} = 10.55$, $e^{5\pi/4} = 50.75$, … — equally spaced on a logarithmic axis. Equations of this type arise for fields in regions bounded by coaxial cylinders or wedges, where there is no natural length scale.

```desmos-graph
left=0; right=12; bottom=-1.3; top=1.3
---
y=\cos(2\ln(x))|blue|x>0
y=0|black
(2.193,0)|label:e^{π/4}
(10.55,0)|label:e^{3π/4}
```
*$y = \cos(2\ln t)$ from Example 11.17: oscillation whose period grows in proportion to $t$, the signature of a scale-invariant (Cauchy–Euler) equation.*

> [!example] Worked Example 11.18 — A repeated indicial root
> Solve $t^2y'' - ty' + y = 0$, $t > 0$.
>
> **Step 1.** $a = -1$, $b = 1$: $m^2 - 2m + 1 = (m - 1)^2$, so $m = 1$ is a double root.
>
> **Step 2.** $y = (c_1 + c_2\ln t)\,t$.
>
> **Check by reduction of order.** $p = -1/t$, $e^{-\int p} = t$, so $y_2 = t\int t/t^2\,dt = t\ln t$. ✓ Direct check: $(t\ln t)' = \ln t + 1$, $(t\ln t)'' = 1/t$; $t^2/t - t(\ln t + 1) + t\ln t = 0$. ✓

### Nonlinear second-order equations: two reductions

Nonlinear second-order equations rarely have closed-form solutions, but two structural features reduce the order by one.

- **$y$ missing**, $F(t, y', y'') = 0$. Put $v = y'$; then $F(t, v, v') = 0$ is first order. Integrate $v$ afterwards to get $y$.
- **$t$ missing (autonomous)**, $F(y, y', y'') = 0$. Put $v = y'$ and regard $v$ as a function of $y$. By the chain rule $y'' = \dfrac{dv}{dt} = \dfrac{dv}{dy}\dfrac{dy}{dt} = v\dfrac{dv}{dy}$, so $F\!\left(y, v, v\dfrac{dv}{dy}\right) = 0$ is first order in the variables $(y, v)$ — an equation in the **phase plane** (§11.9).

For a conservative system $y'' = -U'(y)$ the second reduction gives $v\,dv = -U'(y)\,dy$, i.e. the **energy integral** $\tfrac12v^2 + U(y) = E$: kinetic plus potential energy is constant.

> [!example] Worked Example 11.19 — The pendulum energy integral
> A pendulum obeys $\ddot\theta + (g/\ell)\sin\theta = 0$. It is released from rest at $\theta_0 = 90^\circ$ with $g/\ell = 9.81\ \mathrm{s^{-2}}$. Find the angular speed at the bottom.
>
> **Step 1 ($t$ missing).** With $\omega = \dot\theta$ as a function of $\theta$: $\omega\,\dfrac{d\omega}{d\theta} = -\dfrac g\ell\sin\theta$.
>
> **Step 2 (separate).** $\tfrac12\omega^2 = \dfrac g\ell\cos\theta + E$. At release $\omega = 0$, $\theta = \theta_0$, so $E = -(g/\ell)\cos\theta_0$:
> $$
> \omega^2 = \frac{2g}{\ell}\big(\cos\theta - \cos\theta_0\big).
> $$
>
> **Step 3.** At $\theta = 0$: $\omega = \sqrt{2\cdot9.81\cdot(1 - 0)} = 4.43\ \mathrm{rad/s}$.
>
> **Check.** The small-angle linearisation $\ddot\theta + (g/\ell)\theta = 0$ would give $\omega_{\max} = \theta_0\sqrt{g/\ell} = 1.571\cdot3.132 = 4.92\ \mathrm{rad/s}$, about 11% too high — at $90^\circ$ amplitude the linear model is no longer accurate. The first-order equation in $(\theta, \omega)$ can be separated once more to give the period as an integral (an elliptic integral), which has no elementary closed form.

> [!tip] Engineering note — constant coefficients live in Chapter 12
> When $p$ and $q$ are constants, $y = e^{st}$ reduces $L[y] = 0$ to the algebraic **characteristic equation** $s^2 + ps + q = 0$, and the three root cases give the overdamped, critically damped and underdamped responses of RLC circuits. That theory, with damping ratio, natural frequency and forced response, is the subject of [Ch. 12](12-linear-differential-equations.md). The methods of this section — Wronskian, reduction of order, Cauchy–Euler — are what you use when the coefficients are **not** constant: tapered transmission lines, graded doping profiles, time-varying circuit elements, fields in curvilinear coordinates.

---

## 11.5 Higher-order ODEs

### Motivation

Each independent energy store adds one to the order of a network's equation, so an $n$th-order filter, a multi-stage amplifier with $n$ dominant poles, or a chain of $n$ coupled resonators leads to an $n$th-order ODE. Fourth-order equations govern the bending of beams and cantilevers (MEMS accelerometers and resonators), and higher-order equations appear whenever a system of first-order equations is collapsed into a single equation for one output. The theory of §11.4 carries over almost word for word; the new ingredients are bookkeeping (an $n \times n$ Wronskian) and the systematic reduction to a first-order system, which is how such equations are actually analysed and computed.

### The $n$th-order linear equation

> [!definition] Linear $n$th-order ODE
> $$
> L[y] = y^{(n)} + p_{n-1}(t)\,y^{(n-1)} + \cdots + p_1(t)\,y' + p_0(t)\,y = r(t),
> $$
> with $p_0, \dots, p_{n-1}, r$ continuous on an interval $I$. An **initial-value problem** prescribes $y(t_0), y'(t_0), \dots, y^{(n-1)}(t_0)$.

> [!theorem] Structure of the $n$th-order homogeneous equation
> 1. **Existence and uniqueness.** For every $t_0 \in I$ and every vector of initial data, the IVP has exactly one solution, defined on all of $I$.
> 2. **Dimension.** The solutions of $L[y] = 0$ form a vector space of dimension $n$. Any $n$ linearly independent solutions $y_1, \dots, y_n$ form a **fundamental set**, and the general solution is $y = c_1y_1 + \cdots + c_ny_n$.
> 3. **Wronskian.** With
> $$
> W(t) = \det\begin{pmatrix} y_1 & y_2 & \cdots & y_n \\ y_1' & y_2' & \cdots & y_n' \\ \vdots & \vdots & & \vdots \\ y_1^{(n-1)} & y_2^{(n-1)} & \cdots & y_n^{(n-1)} \end{pmatrix},
> $$
> Abel's formula reads $W(t) = W(t_0)\exp\!\left(-\int_{t_0}^{t}p_{n-1}(s)\,ds\right)$, so $W$ is never zero or always zero on $I$, and the solutions are independent exactly when $W \ne 0$.

*Proof sketch.* Parts 1 and 2 follow exactly as in §11.4, with the initial-data map now onto $\mathbb{R}^n$. For Abel's formula, differentiate the determinant row by row: $W'$ is the sum of $n$ determinants, in each of which one row is differentiated. Differentiating row $k < n$ produces a copy of row $k + 1$, so those determinants have two equal rows and vanish. Only the last row survives, becoming $\big(y_1^{(n)}, \dots, y_n^{(n)}\big)$. Substituting $y_k^{(n)} = -p_{n-1}y_k^{(n-1)} - \cdots - p_0y_k$ and subtracting multiples of the earlier rows leaves $-p_{n-1}$ times the original last row. Hence $W' = -p_{n-1}W$. $\square$

The general solution of $L[y] = r$ contains $n$ arbitrary constants, and $n$ independent conditions are needed to fix them.

> [!example] Worked Example 11.20 — A third-order constant-coefficient equation
> Solve $y''' - y' = 0$ and verify independence of the solutions.
>
> **Step 1.** With constant coefficients, $y = e^{st}$ gives $s^3 - s = s(s - 1)(s + 1) = 0$, so $y = c_1 + c_2e^{t} + c_3e^{-t}$.
>
> **Step 2 (Wronskian).**
> $$
> W = \det\begin{pmatrix} 1 & e^{t} & e^{-t} \\ 0 & e^{t} & -e^{-t} \\ 0 & e^{t} & e^{-t} \end{pmatrix} = 1\cdot\big(e^{t}e^{-t} + e^{t}e^{-t}\big) = 2 \ne 0.
> $$
>
> **Check.** The coefficient of $y''$ is $p_2 = 0$, so Abel predicts a constant $W$. ✓

> [!example] Worked Example 11.21 — A third-order Cauchy–Euler IVP
> Solve $t^3y''' + t^2y'' - 2ty' + 2y = 0$ for $t > 0$ with $y(1) = 0$, $y'(1) = 0$, $y''(1) = 6$.
>
> **Step 1 (indicial equation).** For $y = t^m$, $t^3y''' = m(m - 1)(m - 2)t^m$, $t^2y'' = m(m - 1)t^m$, $ty' = mt^m$:
> $$
> m(m - 1)(m - 2) + m(m - 1) - 2m + 2 = m^3 - 2m^2 - m + 2 = (m - 2)(m - 1)(m + 1) = 0.
> $$
> So $y = c_1t + c_2t^2 + c_3t^{-1}$.
>
> **Step 2 (Wronskian, two ways).** Directly,
> $$
> W = \det\begin{pmatrix} t & t^2 & t^{-1} \\ 1 & 2t & -t^{-2} \\ 0 & 2 & 2t^{-3}\end{pmatrix} = \frac6t.
> $$
> In standard form (divide by $t^3$) the coefficient of $y''$ is $p_2 = 1/t$, so Abel gives $W \propto e^{-\ln t} = 1/t$. ✓ It is nonzero on $t > 0$.
>
> **Step 3 (fit).** With $y' = c_1 + 2c_2t - c_3t^{-2}$ and $y'' = 2c_2 + 2c_3t^{-3}$, at $t = 1$:
> $$
> c_1 + c_2 + c_3 = 0,\qquad c_1 + 2c_2 - c_3 = 0,\qquad 2c_2 + 2c_3 = 6.
> $$
> The third gives $c_2 + c_3 = 3$, so $c_1 = -3$; then $2c_2 - c_3 = 3$, so $c_2 = 2$, $c_3 = 1$.
> $$
> y = 2t^2 - 3t + \frac1t.
> $$
>
> **Check.** $y(1) = 2 - 3 + 1 = 0$, $y'(1) = 4 - 3 - 1 = 0$, $y''(1) = 4 + 2 = 6$. ✓ The determinant of the $3 \times 3$ fitting system is $W(1) = 6 \ne 0$, which is why a unique fit was guaranteed.

### Reduction to a first-order system

> [!theorem] Companion form
> Define the **state** $x_1 = y,\ x_2 = y',\ \dots,\ x_n = y^{(n-1)}$. Then $L[y] = r$ is equivalent to the first-order system
> $$
> \dot{\mathbf{x}} = \mathbf{A}(t)\,\mathbf{x} + \mathbf{b}\,r(t),\qquad
> \mathbf{A} = \begin{pmatrix} 0 & 1 & 0 & \cdots & 0 \\ 0 & 0 & 1 & \cdots & 0 \\ \vdots & & & \ddots & \vdots \\ 0 & 0 & 0 & \cdots & 1 \\ -p_0 & -p_1 & -p_2 & \cdots & -p_{n-1}\end{pmatrix},\qquad
> \mathbf{b} = \begin{pmatrix} 0 \\ \vdots \\ 0 \\ 1\end{pmatrix}.
> $$
> For constant coefficients, $\det(s\mathbf{I} - \mathbf{A}) = s^n + p_{n-1}s^{n-1} + \cdots + p_0$, the characteristic polynomial of the scalar equation.

The first $n - 1$ rows just say "the derivative of each state is the next state"; the last row is the ODE solved for $y^{(n)}$. The initial data $y(t_0), \dots, y^{(n-1)}(t_0)$ are precisely $\mathbf{x}(t_0)$, and the Wronskian matrix is a fundamental matrix of the system (§11.9).

The conversion matters for three reasons:

- **Theory.** Existence, uniqueness and continuous dependence need only be proved once, for first-order systems (§11.7).
- **Computation.** Every numerical ODE solver ([Ch. 22](../part-7-numerical-computational/22-numerical-methods.md)) — Runge–Kutta, SPICE's transient engine, MATLAB's `ode45` — accepts only first-order systems.
- **Modelling.** The physical state of a circuit (capacitor voltages, inductor currents) is a first-order system from the start; the companion form is one particular choice of state, and [Ch. 12](12-linear-differential-equations.md) shows how to choose others.

> [!example] Worked Example 11.22 — Companion form
> Write $y''' + 2y'' - y' + 3y = \sin t$ as a first-order system and confirm the characteristic polynomial.
>
> **Step 1.** States $x_1 = y$, $x_2 = y'$, $x_3 = y''$. Then $\dot x_1 = x_2$, $\dot x_2 = x_3$, and from the equation $\dot x_3 = y''' = -3x_1 + x_2 - 2x_3 + \sin t$.
>
> **Step 2.**
> $$
> \frac{d}{dt}\begin{pmatrix}x_1\\x_2\\x_3\end{pmatrix} = \begin{pmatrix}0 & 1 & 0\\ 0 & 0 & 1\\ -3 & 1 & -2\end{pmatrix}\begin{pmatrix}x_1\\x_2\\x_3\end{pmatrix} + \begin{pmatrix}0\\0\\1\end{pmatrix}\sin t.
> $$
>
> **Check.** Expanding along the first column, $\det(s\mathbf{I} - \mathbf{A}) = s\big(s(s + 2) - 1\big) + 3\cdot1 = s^3 + 2s^2 - s + 3$. ✓

### Operator factorisation

Write $D = d/dt$. A constant-coefficient operator is a polynomial in $D$, and polynomials in $D$ **commute** and **factor** like ordinary polynomials: $D^2 - 3D + 2 = (D - 1)(D - 2)$. A factored equation can be solved as a chain of first-order linear equations, each by the integrating factor of §11.3. This is the cleanest way to see where the extra factors $t^k$ for repeated roots come from.

> [!example] Worked Example 11.23 — Solving by successive first-order steps
> Find a particular solution of $y''' - 3y'' + 3y' - y = e^{t}$.
>
> **Step 1 (factor).** The operator is $(D - 1)^3$. Put $u = (D - 1)^2y$ and $w = (D - 1)y$, so the equation is $(D - 1)u = e^t$, then $(D - 1)w = u$, then $(D - 1)y = w$.
>
> **Step 2 (three integrating factors).** Each step is $z' - z = g$, with solution $z = e^{t}\int e^{-t}g\,dt$ (constants dropped, since we need only one particular solution):
> - $u = e^t\int e^{-t}e^{t}\,dt = te^t$;
> - $w = e^t\int e^{-t}te^{t}\,dt = \tfrac12t^2e^t$;
> - $y_p = e^t\int \tfrac12t^2\,dt = \tfrac16t^3e^t$.
>
> **Check.** Symbolic substitution of $y_p = \tfrac16t^3e^t$ returns exactly $e^t$. ✓
>
> **Interpretation.** The input $e^{t}$ coincides with a triple root of the characteristic polynomial, and each coincident factor contributes one power of $t$. This is *resonance* in its purest form; [Ch. 12](12-linear-differential-equations.md) turns the observation into the modification rule of undetermined coefficients.

Variable-coefficient operators can sometimes be factored too, but the factors no longer commute: $(D + t)D\,y = y'' + ty'$ whereas $D(D + t)\,y = y'' + ty' + y$. Factor in the right order.

### Reduction of order for $n$th-order equations

If one solution $y_1$ of $L[y] = 0$ is known, substituting $y = v\,y_1$ again cancels all terms in $v$ itself (because $L[y_1] = 0$), leaving an equation of order $n - 1$ for $w = v'$. Each known solution therefore lowers the order by one. In practice this is used mostly for $n = 2$; for higher orders it is a theoretical tool showing that knowledge of solutions can be traded for reduced work.

> [!warning] Common pitfall — too few or too many conditions
> An $n$th-order equation needs exactly $n$ conditions. Fewer leaves free constants; more over-determines the problem, which usually has no solution. When a model gives you "too many" conditions, one of them is either redundant or signals that the order of the model is wrong — a classic symptom of an ignored parasitic element.

> [!tip] Engineering note — where higher order comes from
> Each independent energy-storage element (a capacitor or inductor that cannot be merged with another by series/parallel combination, and is not in a loop of capacitors or a cutset of inductors) adds one to the order. A third-order Butterworth low-pass filter obeys a third-order ODE; a fourth-order Euler–Bernoulli beam equation $EI\,w'''' = q$ describes a MEMS cantilever. In practice engineers rarely solve high-order ODEs as single equations. They factor the characteristic polynomial ([Ch. 12](12-linear-differential-equations.md)), work with transfer functions ([Ch. 14](../part-5-transforms/14-laplace-transforms.md)), or integrate the first-order system numerically.

---

## 11.6 Homogeneous and nonhomogeneous equations

### Motivation

A circuit left to itself after a switching event obeys a homogeneous equation; the same circuit driven by a source obeys a nonhomogeneous one. Linearity lets us treat the two effects separately and add the results — the mathematical content of the superposition theorem of circuit analysis. This section makes that precise and gives a method, **variation of parameters**, that produces a particular solution for *any* continuous input and *any* linear equation whose homogeneous solutions are known — including variable-coefficient equations where the guessing methods of [Ch. 12](12-linear-differential-equations.md) fail.

> [!definition] Homogeneous and nonhomogeneous linear equations
> With $L[y] = y^{(n)} + p_{n-1}y^{(n-1)} + \cdots + p_0y$, the equation $L[y] = 0$ is **homogeneous** and $L[y] = r(t)$ with $r \not\equiv 0$ is **nonhomogeneous**. The function $r(t)$ is called the *forcing*, *input* or *source term*. The general solution $y_h$ of $L[y] = 0$ is the **complementary function**; any single solution $y_p$ of $L[y] = r$ is a **particular solution**.

> [!warning] Common pitfall — two meanings of "homogeneous"
> In §11.2, "homogeneous" meant $y' = F(y/t)$, a scaling property of a first-order equation. For linear equations it means "no forcing term". The two ideas are unrelated. In this book "homogeneous" means the linear sense unless stated otherwise.

### Structure of solutions

> [!theorem] Structure theorem
> Let $y_p$ be **any** particular solution of $L[y] = r$ on $I$ and $y_h = c_1y_1 + \cdots + c_ny_n$ the general solution of $L[y] = 0$. Then every solution of $L[y] = r$ on $I$ has the form
> $$
> y = y_h + y_p = c_1y_1 + \cdots + c_ny_n + y_p.
> $$

*Proof.* If $y$ is any solution, $L[y - y_p] = L[y] - L[y_p] = r - r = 0$, so $y - y_p$ is a homogeneous solution, i.e. equals $c_1y_1 + \cdots + c_ny_n$ for some constants. Conversely every such sum is a solution. $\square$

Geometrically, the solution set of $L[y] = r$ is an **affine subspace**: the $n$-dimensional solution space of $L[y] = 0$ translated by $y_p$ — the same picture as the solution set of a consistent linear system $\mathbf{A}\mathbf{x} = \mathbf{b}$ in [Ch. 8](../part-3-linear-algebra/08-matrices-and-linear-systems.md). Different choices of $y_p$ just shift the constants $c_k$.

> [!theorem] Superposition of inputs
> If $L[y_a] = r_a$ and $L[y_b] = r_b$, then $L[\alpha y_a + \beta y_b] = \alpha r_a + \beta r_b$ for any constants $\alpha, \beta$.

A complicated input can therefore be split into simple pieces, a particular solution found for each, and the results added. This is the basis of the Fourier method ([Ch. 15](../part-5-transforms/15-fourier-series.md)): split a periodic input into sinusoids, solve for each, and sum.

> [!example] Worked Example 11.24 — Splitting the input
> Find the general solution of $y'' + 4y = 8 + 5e^{t}$.
>
> **Step 1 (homogeneous).** $s^2 + 4 = 0$ gives $y_h = c_1\cos2t + c_2\sin2t$.
>
> **Step 2 (input $8$).** A constant $y_{p1} = A$ gives $4A = 8$, so $y_{p1} = 2$.
>
> **Step 3 (input $5e^t$).** Try $y_{p2} = Be^{t}$: $Be^t + 4Be^t = 5e^t$, so $B = 1$.
>
> **Step 4.** $y = c_1\cos2t + c_2\sin2t + 2 + e^{t}$.
>
> **Check.** $y_p'' + 4y_p = e^t + 8 + 4e^t = 8 + 5e^t$. ✓

The guesses in steps 2–3 are the method of **undetermined coefficients**, which works for constant-coefficient equations with inputs built from polynomials, exponentials and sinusoids. It is developed with all its special cases in [Ch. 12](12-linear-differential-equations.md). Here we need a method with no such restrictions.

### Variation of parameters

The idea (due to Lagrange) is to look for a particular solution of the same form as $y_h$ but with the constants replaced by functions: $y_p = u_1(t)y_1(t) + u_2(t)y_2(t)$. Two unknown functions are more freedom than one equation needs, so we may impose one extra condition of our choosing — and we choose it to keep second derivatives of $u_k$ out of the calculation.

> [!theorem] Variation of parameters (second order)
> Let $y_1, y_2$ be a fundamental set for $y'' + p(t)y' + q(t)y = 0$ on $I$, with Wronskian $W = y_1y_2' - y_1'y_2$. Then a particular solution of $y'' + p(t)y' + q(t)y = r(t)$ is
> $$
> y_p(t) = -y_1(t)\int \frac{y_2(t)\,r(t)}{W(t)}\,dt + y_2(t)\int \frac{y_1(t)\,r(t)}{W(t)}\,dt.
> $$

*Derivation.* With $y_p = u_1y_1 + u_2y_2$,

$$
y_p' = u_1y_1' + u_2y_2' + \underbrace{(u_1'y_1 + u_2'y_2)}_{\text{set to } 0}.
$$

Imposing $u_1'y_1 + u_2'y_2 = 0$, differentiate again: $y_p'' = u_1y_1'' + u_2y_2'' + u_1'y_1' + u_2'y_2'$. Substituting into the equation and grouping,

$$
u_1\underbrace{(y_1'' + py_1' + qy_1)}_{0} + u_2\underbrace{(y_2'' + py_2' + qy_2)}_{0} + u_1'y_1' + u_2'y_2' = r.
$$

So $u_1', u_2'$ satisfy the linear system

$$
\begin{pmatrix} y_1 & y_2 \\ y_1' & y_2'\end{pmatrix}\begin{pmatrix}u_1'\\u_2'\end{pmatrix} = \begin{pmatrix}0\\ r\end{pmatrix},
$$

whose determinant is the Wronskian, nonzero on $I$. Cramer's rule ([Ch. 8](../part-3-linear-algebra/08-matrices-and-linear-systems.md)) gives $u_1' = -y_2r/W$ and $u_2' = y_1r/W$. Integrating gives the formula. $\square$

**General order.** For $L[y] = r$ of order $n$, impose $\sum u_k'y_k^{(j)} = 0$ for $j = 0, \dots, n - 2$; the last equation becomes $\sum u_k'y_k^{(n-1)} = r$. Cramer's rule gives $u_k' = r\,W_k/W$, where $W_k$ is $W$ with its $k$th column replaced by $(0, \dots, 0, 1)^\mathsf{T}$.

> [!warning] Common pitfall — use the standard form
> The $r(t)$ in the formula is the right side **after** dividing by the coefficient of $y''$. For $t^2y'' - 2ty' + 2y = t^3$ it is $r = t$, not $t^3$. Using the raw right side is the single most common error in variation of parameters.

> [!example] Worked Example 11.25 — An input undetermined coefficients cannot handle
> Find a particular solution of $y'' + y = \sec t$ on $\lvert t\rvert < \pi/2$.
>
> **Step 1.** $y_1 = \cos t$, $y_2 = \sin t$, $W = \cos^2 t + \sin^2 t = 1$.
>
> **Step 2.**
> $$
> \begin{aligned}
> u_1 &= -\int \sin t\sec t\,dt = -\int\tan t\,dt = \ln\lvert\cos t\rvert, \\
> u_2 &= \int\cos t\sec t\,dt = t.
> \end{aligned}
> $$
>
> **Step 3.** $y_p = \cos t\,\ln(\cos t) + t\sin t$.
>
> **Check.** $y_p' = -\sin t\ln\cos t - \sin t + \sin t + t\cos t = -\sin t\ln\cos t + t\cos t$, and $y_p'' = -\cos t\ln\cos t + \dfrac{\sin^2t}{\cos t} + \cos t - t\sin t$. Adding $y_p$: $\dfrac{\sin^2t}{\cos t} + \cos t = \dfrac{1}{\cos t} = \sec t$. ✓ The input has no finite family of derivatives, so no guess of undetermined-coefficients type can succeed.

> [!example] Worked Example 11.26 — Variable coefficients with initial conditions
> Solve $t^2y'' - 2ty' + 2y = t^3$, $t > 0$, with $y(1) = 0$, $y'(1) = 0$.
>
> **Step 1 (homogeneous).** From Example 11.16, $y_1 = t$, $y_2 = t^2$, $W = t^2$.
>
> **Step 2 (standard form).** $y'' - \dfrac2ty' + \dfrac{2}{t^2}y = t$, so $r = t$.
>
> **Step 3 (parameters).** $u_1' = -\dfrac{y_2r}{W} = -\dfrac{t^2\cdot t}{t^2} = -t$, so $u_1 = -\tfrac12t^2$; $u_2' = \dfrac{y_1r}{W} = \dfrac{t\cdot t}{t^2} = 1$, so $u_2 = t$.
>
> **Step 4.** $y_p = -\tfrac12t^2\cdot t + t\cdot t^2 = \tfrac12t^3$, and $y = c_1t + c_2t^2 + \tfrac12t^3$.
>
> **Step 5 (conditions).** $y(1) = c_1 + c_2 + \tfrac12 = 0$ and $y'(1) = c_1 + 2c_2 + \tfrac32 = 0$ give $c_2 = -1$, $c_1 = \tfrac12$:
> $$
> y = \tfrac12t - t^2 + \tfrac12t^3 = \tfrac12\,t\,(t - 1)^2.
> $$
>
> **Check.** With $y = \tfrac12(t^3 - 2t^2 + t)$: $y' = \tfrac12(3t^2 - 4t + 1)$, $y'' = 3t - 2$. Then $t^2(3t - 2) - t(3t^2 - 4t + 1) + (t^3 - 2t^2 + t) = t^3$. ✓

### The Green's-function (impulse-response) form

Variation of parameters can be written as a single integral that builds in zero initial conditions. Taking definite integrals from $t_0$,

$$
y_p(t) = \int_{t_0}^{t} G(t, s)\,r(s)\,ds,\qquad G(t, s) = \frac{y_1(s)\,y_2(t) - y_1(t)\,y_2(s)}{W(s)}.
$$

> [!theorem] Properties of the causal Green's function
> For fixed $s$, $G(t, s)$ as a function of $t$ solves the homogeneous equation with $G(s, s) = 0$ and $\partial G/\partial t\,(s, s) = 1$. Consequently $y_p$ above is the unique solution of $L[y] = r$ with $y(t_0) = y'(t_0) = 0$, and the solution with general initial data is $y_h + y_p$ with $y_h$ fitted to those data.

*Why.* $G(\cdot, s)$ is a combination of $y_1, y_2$, so it solves $L[y] = 0$; $G(s, s) = 0$ by inspection, and $\partial_tG(s, s) = [y_1(s)y_2'(s) - y_1'(s)y_2(s)]/W(s) = 1$. Differentiating $y_p$ under the integral sign (Leibniz rule, [Ch. 5](../part-2-calculus/05-integral-calculus.md)) then gives $y_p' = \int\partial_tG\,r\,ds$ and $y_p'' = r(t) + \int\partial_t^2G\,r\,ds$, so $L[y_p] = r$. $\square$

Physically, $G(t, s)$ is the response at time $t$ to a unit **impulse** of input applied at time $s$: the impulse leaves position unchanged and kicks the velocity by 1. The integral says the response to a general input is the superposition of responses to the impulses $r(s)\,ds$ that make it up. When the coefficients are constant, $G(t, s) = g(t - s)$ depends only on the elapsed time, and $y_p$ is the **convolution** $g * r$ — the impulse-response description of a linear time-invariant system used throughout [Ch. 12](12-linear-differential-equations.md), [Ch. 14](../part-5-transforms/14-laplace-transforms.md) and [Ch. 16](../part-5-transforms/16-fourier-transforms.md).

> [!example] Worked Example 11.27 — The undamped oscillator: step and resonant inputs
> For $y'' + \omega_0^2y = r(t)$ with $y(0) = y'(0) = 0$, find the response to (a) a step $r = F$, (b) a resonant input $r = \cos\omega_0t$.
>
> **Step 1 (Green's function).** $y_1 = \cos\omega_0t$, $y_2 = \sin\omega_0t$, $W = \omega_0$. Then
> $$
> G(t, s) = \frac{\cos\omega_0s\sin\omega_0t - \cos\omega_0t\sin\omega_0s}{\omega_0} = \frac{\sin\omega_0(t - s)}{\omega_0}.
> $$
>
> **Step 2 (step).** $y = \dfrac{F}{\omega_0}\displaystyle\int_0^t\sin\omega_0(t - s)\,ds = \dfrac{F}{\omega_0^2}\big(1 - \cos\omega_0t\big)$.
>
> **Step 3 (resonant input).** Using $\sin A\cos B = \tfrac12[\sin(A + B) + \sin(A - B)]$,
> $$
> y = \frac{1}{\omega_0}\int_0^t\sin\omega_0(t - s)\cos\omega_0s\,ds = \frac{1}{2\omega_0}\int_0^t\big[\sin\omega_0t + \sin\omega_0(t - 2s)\big]\,ds = \frac{t\sin\omega_0t}{2\omega_0}.
> $$
> (The second integrand integrates to zero over $[0, t]$.)
>
> **Check.** For (b), $y' = \dfrac{\sin\omega_0t}{2\omega_0} + \dfrac t2\cos\omega_0t$ and $y'' = \cos\omega_0t - \dfrac{\omega_0t}{2}\sin\omega_0t$, so $y'' + \omega_0^2y = \cos\omega_0t$; $y(0) = y'(0) = 0$. ✓
>
> **Interpretation.** The step response oscillates about the new equilibrium $F/\omega_0^2$ forever (no damping). The resonant response grows linearly without bound — the envelope $t/(2\omega_0)$ — which is why lossless resonators driven at their natural frequency are dangerous, and why real circuits always contain the resistance that [Ch. 12](12-linear-differential-equations.md) uses to limit the peak.

```desmos-graph
left=0; right=30; bottom=-16; top=16
---
y=\frac{x\sin(x)}{2}|blue
y=\frac{x}{2}|red|dashed
y=-\frac{x}{2}|red|dashed
y=1-\cos(x)|green
```
*Example 11.27 with $\omega_0 = 1$: resonant response $\tfrac12t\sin t$ (blue) inside its linearly growing envelope (red dashed), and the bounded step response $1 - \cos t$ (green).*

> [!tip] Engineering note — natural and forced response
> In circuit language $y_h$ is the **natural response** (shape set by the circuit, amplitude set by initial conditions) and $y_p$ a **forced response** (shape set by the source). For stable circuits $y_h$ decays and $y_p$ — chosen as the steady periodic solution — is what a meter reads after the transient. [Ch. 12](12-linear-differential-equations.md) develops this vocabulary for constant-coefficient circuits; the structure theorem here is what guarantees that the decomposition is always possible and complete.

---

## 11.7 Initial-value problems

### Motivation

An engineer who simulates a circuit is solving an initial-value problem: the state at $t = 0$ is known (capacitor voltages, inductor currents), and the simulator marches forward. Three questions must have good answers before the result means anything. Does a solution **exist**? Is it **unique** — or could the same initial state lead to different futures? Does it depend **continuously** on the data, so that small errors in component values or initial conditions produce small errors in the answer? A problem with all three properties is called **well-posed**. This section states when that is guaranteed, sketches why, and shows what goes wrong when the hypotheses fail.

> [!definition] Initial-value problem (IVP)
> For a first-order system, find $\mathbf{y}(t)$ with
> $$
> \mathbf{y}' = \mathbf{f}(t, \mathbf{y}),\qquad \mathbf{y}(t_0) = \mathbf{y}_0.
> $$
> For an $n$th-order scalar equation, the data are $y(t_0), y'(t_0), \dots, y^{(n-1)}(t_0)$ — all at the **same** point $t_0$ — which is the same thing after the companion-form reduction of §11.5.

We state results for scalar first-order equations $y' = f(t, y)$; with $\lvert\cdot\rvert$ replaced by a vector norm $\lVert\cdot\rVert$ everything holds verbatim for systems, and hence for $n$th-order equations.

### The Lipschitz condition

> [!definition] Lipschitz condition
> $f(t, y)$ is **Lipschitz in $y$** on a region $R$ with constant $K$ if
> $$
> \lvert f(t, y) - f(t, z)\rvert \le K\,\lvert y - z\rvert\quad\text{for all } (t, y), (t, z) \in R.
> $$

By the mean-value theorem, if $\partial f/\partial y$ is continuous on a closed rectangle then $f$ is Lipschitz there with $K = \max\lvert\partial f/\partial y\rvert$. So the practical test is: **is $\partial f/\partial y$ bounded near the initial point?** The function $f = y^{2/3}$ fails at $y = 0$ (its derivative $\tfrac23y^{-1/3}$ is unbounded); $f = \lvert y\rvert$ is Lipschitz with $K = 1$ even though it is not differentiable at $0$.

### Existence and uniqueness

> [!theorem] Picard–Lindelöf
> Let $f$ be continuous on the rectangle $R = \{\lvert t - t_0\rvert \le a,\ \lvert y - y_0\rvert \le b\}$ and Lipschitz in $y$ there. Let $M = \max_R\lvert f\rvert$ and $h = \min(a, b/M)$. Then the IVP $y' = f(t, y)$, $y(t_0) = y_0$ has **exactly one** solution on $\lvert t - t_0\rvert \le h$.

*Proof sketch.* **(1) Integral equation.** A continuous $y$ solves the IVP if and only if

$$
y(t) = y_0 + \int_{t_0}^{t} f\big(s, y(s)\big)\,ds,
$$

by the fundamental theorem of calculus. This form needs no derivatives and builds in the initial condition.

**(2) Picard iteration.** Start from $y_0(t) \equiv y_0$ and define $y_{k+1}(t) = y_0 + \int_{t_0}^{t}f(s, y_k(s))\,ds$. The choice $h \le b/M$ keeps every iterate inside $R$: $\lvert y_{k+1} - y_0\rvert \le M\lvert t - t_0\rvert \le b$.

**(3) Convergence.** By the Lipschitz condition, $\lvert y_{k+1}(t) - y_k(t)\rvert \le K\int_{t_0}^{t}\lvert y_k - y_{k-1}\rvert\,ds$. Starting from $\lvert y_1 - y_0\rvert \le M\lvert t - t_0\rvert$ and inducting,

$$
\lvert y_{k+1}(t) - y_k(t)\rvert \le \frac{MK^k\lvert t - t_0\rvert^{k+1}}{(k + 1)!}.
$$

The right side is the general term of a convergent exponential series, so $y_k$ converges uniformly to a continuous limit $y$, and passing to the limit in the iteration shows that $y$ satisfies the integral equation.

**(4) Uniqueness.** If $y$ and $z$ both solve it, $u = \lvert y - z\rvert$ satisfies $u(t) \le K\int_{t_0}^{t}u\,ds$ with $u(t_0) = 0$, and Grönwall's inequality (below) forces $u \equiv 0$. $\square$

Two weaker or complementary statements are worth knowing:

- **Peano's theorem.** Continuity of $f$ alone guarantees that *at least one* solution exists near $t_0$, but not that it is unique.
- **Continuation.** A solution can be continued as long as it stays in a region where the hypotheses hold. The maximal solution either exists for all $t \ge t_0$ or leaves every bounded region — for a scalar equation on all of $\mathbb{R}^2$, it **blows up** ($\lvert y\rvert \to \infty$) in finite time. There is no third option such as "stops for no reason".

> [!example] Worked Example 11.28 — Picard iterates for $y' = y$
> Compute the Picard iterates for $y' = y$, $y(0) = 1$.
>
> **Step 1.** $y_0 = 1$.
>
> **Step 2.** $y_1 = 1 + \int_0^t1\,ds = 1 + t$; $y_2 = 1 + \int_0^t(1 + s)\,ds = 1 + t + \tfrac{t^2}{2}$; and by induction
> $$
> y_k(t) = \sum_{m=0}^{k}\frac{t^m}{m!}.
> $$
>
> **Step 3.** The iterates are the Taylor partial sums of $e^t$ and converge to $y = e^t$ for every $t$ ([Ch. 4](../part-2-calculus/04-differential-calculus.md)).
>
> **Check against the error bound.** On $\lvert t\rvert \le 1$ with $K = 1$, the bound on $\lvert y_5 - y_4\rvert$ at $t = 1$ is $MK^4/5!$ with $M \le e$ in a suitable rectangle, i.e. about $0.023$; the actual difference is $1/5! = 0.0083$. ✓ The iteration is a constructive proof and, for simple $f$, a usable approximation method.

> [!example] Worked Example 11.29 — How long is the guaranteed interval?
> For $y' = 1 + y^2$, $y(0) = 0$, find the interval of existence promised by Picard–Lindelöf and compare with the true solution.
>
> **Step 1.** $f = 1 + y^2$ is continuous with $f_y = 2y$ bounded on any rectangle, so the theorem applies for every $a, b$. Take $a$ large; on $\lvert y\rvert \le b$, $M = 1 + b^2$, so $h = b/(1 + b^2)$.
>
> **Step 2 (best choice).** $\dfrac{d}{db}\dfrac{b}{1 + b^2} = \dfrac{1 - b^2}{(1 + b^2)^2} = 0$ at $b = 1$, giving $h = \tfrac12$.
>
> **Step 3 (true solution).** Separating, $\arctan y = t$, so $y = \tan t$, which exists on $(-\pi/2, \pi/2)$ and blows up at $\lvert t\rvert = \pi/2 \approx 1.571$.
>
> **Interpretation.** The theorem is conservative — it guarantees $\lvert t\rvert \le 0.5$ where the truth is $\lvert t\rvert < 1.571$ — but it is *local by nature*: no choice of rectangle could guarantee all $t$, because the solution genuinely escapes. Applying the theorem again from $t = 0.5$ extends the solution step by step until the blow-up.

### When the hypotheses fail

> [!example] Worked Example 11.30 — Blow-up and non-uniqueness
> **(a) Finite-time blow-up.** $y' = y^2$, $y(0) = 1$. Here $f = y^2$ is smooth everywhere, and separation gives $y = 1/(1 - t)$. The solution exists only for $t < 1$, although the equation itself is defined everywhere. The cause is the superlinear growth of $f$: $K = 2\lvert y\rvert$ grows without bound as $y$ grows.
>
> **(b) Non-uniqueness.** $y' = 3y^{2/3}$, $y(0) = 0$. Both $y \equiv 0$ and $y = t^3$ are solutions: $(t^3)' = 3t^2 = 3(t^3)^{2/3}$. In fact, for any $c \ge 0$ the function equal to $0$ for $t \le c$ and $(t - c)^3$ for $t > c$ is also a solution — infinitely many futures from the same initial state. Uniqueness fails because $\partial f/\partial y = 2y^{-1/3}$ is unbounded at $y = 0$; existence still holds (Peano).

```desmos-graph
left=-2; right=3; bottom=-2; top=8
---
y=\frac{1}{1-x}|blue|x<1
x=1|red|dashed
y=x^{3}|green|x>0
y=(x-1)^{3}|green|dashed|x>1
y=0|green|x>0
```
*Blue: $y = 1/(1 - t)$ escapes to infinity at $t = 1$ (red dashed asymptote). Green: three different solutions of $y' = 3y^{2/3}$ leave the same initial point $(0, 0)$ — $y = 0$, $y = t^3$, and $y$ that waits until $t = 1$ before departing.*

### Continuous dependence and Grönwall's inequality

> [!theorem] Grönwall's inequality
> If $u(t) \ge 0$ is continuous and $u(t) \le \alpha + K\int_{t_0}^{t}u(s)\,ds$ for $t \ge t_0$, with constants $\alpha \ge 0$, $K \ge 0$, then $u(t) \le \alpha\,e^{K(t - t_0)}$.

*Sketch.* Let $U(t) = \alpha + K\int_{t_0}^{t}u\,ds$. Then $U' = Ku \le KU$, so $(e^{-K(t - t_0)}U)' \le 0$ and $U(t) \le U(t_0)e^{K(t - t_0)} = \alpha e^{K(t - t_0)}$; finally $u \le U$. $\square$

> [!theorem] Continuous dependence on initial data
> If $f$ is Lipschitz in $y$ with constant $K$ and $y, z$ solve $y' = f(t, y)$ with initial values $y_0, z_0$, then
> $$
> \lvert y(t) - z(t)\rvert \le \lvert y_0 - z_0\rvert\,e^{K\lvert t - t_0\rvert}.
> $$

*Proof.* Subtract the two integral equations and use the Lipschitz bound: $\lvert y - z\rvert \le \lvert y_0 - z_0\rvert + K\int\lvert y - z\rvert\,ds$. Apply Grönwall. With $y_0 = z_0$ this is the uniqueness proof. $\square$

The bound is a worst case: errors *may* grow like $e^{Kt}$. For a stable system they actually shrink — for $y' = -y/\tau$ two solutions approach each other like $e^{-t/\tau}$ — but for an unstable one (an oscillator start-up, a latch near its metastable point, Example 11.2) the exponential growth of small differences is real and is exactly what limits predictability. A similar estimate shows continuous dependence on parameters in $f$: a small change in a component value changes the solution by a small amount on any finite time interval.

### Linear equations: global existence

> [!theorem] Linear IVPs
> If $p_0, \dots, p_{n-1}, r$ (or $\mathbf{A}(t), \mathbf{f}(t)$ for a system) are continuous on an open interval $I$ containing $t_0$, the linear IVP has exactly one solution, and it exists on **all** of $I$.

*Why linear equations cannot blow up.* For $\mathbf{y}' = \mathbf{A}(t)\mathbf{y} + \mathbf{f}(t)$ the Lipschitz constant on any closed subinterval $J \subset I$ is $K = \max_J\lVert\mathbf{A}\rVert$, the same for *all* $\mathbf{y}$ — no rectangle in $y$ is needed. Grönwall then bounds $\lVert\mathbf{y}(t)\rVert$ by $\big(\lVert\mathbf{y}_0\rVert + \int_J\lVert\mathbf{f}\rVert\big)e^{K\lvert t - t_0\rvert}$ on $J$, so the solution stays bounded on every closed subinterval and, by the continuation property, extends across all of $I$. The only places a linear solution can misbehave are where the *coefficients* do — the singular points — and those are visible in the equation before solving it.

> [!example] Worked Example 11.31 — Fitting circuit initial conditions
> In a series RLC circuit with $R = 4\ \Omega$, $L = 1\ \mathrm H$, $C = \tfrac13\ \mathrm F$, the capacitor holds $v_C(0^-) = 2\ \mathrm V$ and the inductor current is $i(0^-) = 0$. At $t = 0$ a $10\ \mathrm V$ source is switched into the loop. Find $v_C(t)$ for $t > 0$.
>
> **Step 1 (equation).** KVL: $Li' + Ri + v_C = 10$ with $i = C\dot v_C$ gives $LC\ddot v_C + RC\dot v_C + v_C = 10$, i.e. $\ddot v_C + 4\dot v_C + 3v_C = 30$.
>
> **Step 2 (initial data at $0^+$).** Capacitor voltage and inductor current are continuous: $v_C(0^+) = 2\ \mathrm V$ and $i(0^+) = 0$, hence $\dot v_C(0^+) = i(0^+)/C = 0$.
>
> **Step 3 (general solution).** The constant input gives $y_p = 10$; the homogeneous equation has $s^2 + 4s + 3 = (s + 1)(s + 3)$, so $v_C = 10 + Ae^{-t} + Be^{-3t}$.
>
> **Step 4 (fit).** $A + B = -8$ and $-A - 3B = 0$, so $B = 4$, $A = -12$:
> $$
> v_C(t) = 10 - 12e^{-t} + 4e^{-3t}\ \mathrm V.
> $$
>
> **Check.** $v_C(0) = 10 - 12 + 4 = 2$, $\dot v_C(0) = 12 - 12 = 0$. ✓ At $t = 1\ \mathrm s$: $v_C = 10 - 4.415 + 0.199 = 5.78\ \mathrm V$. The coefficients are constant, so this is the unique solution for all $t > 0$.

> [!tip] Engineering note — initial conditions in circuits
> A capacitor voltage cannot jump unless the current is infinite, and an inductor current cannot jump unless the voltage is infinite. The natural initial conditions for a circuit are therefore $v_C(0^-)$ and $i_L(0^-)$, the stored energies just before a switching event — exactly one number per state variable, which is why state variables make the IVP automatically well-posed. Derivatives such as $\dot v_C(0^+) = i_C(0^+)/C$ follow from the circuit at $t = 0^+$. When an idealised circuit violates the rule (a capacitor switched directly across a voltage source, an inductor current interrupted), the model needs impulses — handled systematically by the Laplace transform in [Ch. 14](../part-5-transforms/14-laplace-transforms.md).

> [!warning] Common pitfall — trusting a simulator at a non-Lipschitz point
> Models with square-root or fractional-power laws (ideal square-law MOSFET at threshold, some magnetic and friction models) can be non-Lipschitz at a point. At such points a numerical integrator may follow any of several solutions depending on step size and round-off, and the result can change with solver settings. Smoothing the model (e.g. replacing $\lvert v\rvert^{1/2}$ by $(v^2 + \varepsilon^2)^{1/4}$) restores uniqueness.

---

## 11.8 Boundary-value problems

### Motivation

When the independent variable is **position** rather than time, conditions are naturally given at two places: the potential on both electrodes, the voltage at the source and load ends of a line, the temperature at both ends of a heat sink, the clamped and free ends of a cantilever. These are **boundary-value problems** (BVPs). They behave very differently from IVPs — existence and uniqueness are no longer automatic — and their homogeneous versions, the **eigenvalue problems**, produce the modes (standing waves, resonant frequencies, orthogonal function sets) on which separation of variables ([Ch. 13](13-partial-differential-equations.md)) and Fourier series ([Ch. 15](../part-5-transforms/15-fourier-series.md)) are built.

> [!definition] Two-point boundary-value problem
> A second-order ODE on $a \le x \le b$ with one condition at **each end**. Common boundary conditions:
> - **Dirichlet:** $y(a) = \alpha$, $y(b) = \beta$ (value prescribed — e.g. electrode potential);
> - **Neumann:** $y'(a) = \alpha$, $y'(b) = \beta$ (slope prescribed — e.g. zero field or zero heat flux);
> - **Robin (mixed):** $c_1y(a) + c_2y'(a) = \alpha$ (a linear combination — e.g. convective cooling, or a line terminated in a resistance).
>
> The BVP is **homogeneous** if the equation and both boundary conditions are homogeneous (zero right sides).

### Solvability: one, none or infinitely many

For a linear equation $y'' + p(x)y' + q(x)y = r(x)$ with fundamental set $y_1, y_2$ and particular solution $y_p$, write $y = c_1y_1 + c_2y_2 + y_p$ and impose the two boundary conditions. For Dirichlet conditions this gives the $2 \times 2$ system

$$
\begin{pmatrix} y_1(a) & y_2(a) \\ y_1(b) & y_2(b)\end{pmatrix}\begin{pmatrix}c_1\\c_2\end{pmatrix} = \begin{pmatrix}\alpha - y_p(a)\\ \beta - y_p(b)\end{pmatrix}.
$$

Everything now follows from linear algebra ([Ch. 8](../part-3-linear-algebra/08-matrices-and-linear-systems.md)): if the determinant $D = y_1(a)y_2(b) - y_2(a)y_1(b)$ is nonzero there is exactly one solution for every right side; if $D = 0$ there is either none or infinitely many, depending on whether the right side lies in the column space.

> [!theorem] Fredholm alternative (two-point BVPs)
> For a linear second-order BVP, exactly one of the following holds:
> 1. The homogeneous BVP has only the trivial solution $y \equiv 0$. Then the nonhomogeneous BVP has **exactly one** solution for every $r$, $\alpha$, $\beta$.
> 2. The homogeneous BVP has a nontrivial solution. Then the nonhomogeneous BVP has **either no solution or infinitely many**, depending on the data.

*Proof.* The homogeneous BVP corresponds to a zero right side, and it has a nontrivial solution exactly when $D = 0$. Apply the corresponding statement for $2 \times 2$ linear systems. $\square$

An IVP never meets case 2, because its "determinant" is the Wronskian at $t_0$, which is never zero. A BVP meets it whenever the interval length is "tuned" to a natural oscillation of the equation — the spatial analogue of resonance.

> [!example] Worked Example 11.32 — Three possible outcomes
> The equation $y'' + y = 0$ has general solution $y = c_1\cos x + c_2\sin x$.
>
> 1. $y(0) = 0$, $y(\pi/2) = 1$: $c_1 = 0$ and $c_2 = 1$, so the **unique** solution is $y = \sin x$. Here $D = \cos0\sin\tfrac\pi2 - \sin0\cos\tfrac\pi2 = 1 \ne 0$.
> 2. $y(0) = 0$, $y(\pi) = 1$: $c_1 = 0$, then $c_2\sin\pi = 0 \ne 1$, so there is **no solution**.
> 3. $y(0) = 0$, $y(\pi) = 0$: $c_1 = 0$ and $c_2$ is arbitrary, so there are **infinitely many** solutions $y = c_2\sin x$.
>
> In cases 2 and 3, $D = \cos0\sin\pi - \sin0\cos\pi = 0$: the interval $[0, \pi]$ is exactly half a period of the natural solution $\sin x$, so the homogeneous problem has the nontrivial solution $\sin x$.

> [!example] Worked Example 11.33 — A Poisson-type BVP
> In one dimension the potential between two plates with a space-charge density obeys Poisson's equation $V'' = -\rho/\varepsilon$. As a dimensionless model, solve $y'' = 6x$ with $y(0) = 0$, $y(1) = 2$.
>
> **Step 1.** Integrate twice: $y = x^3 + c_1x + c_0$.
>
> **Step 2.** $y(0) = c_0 = 0$ and $y(1) = 1 + c_1 = 2$, so $c_1 = 1$:
> $$
> y = x^3 + x.
> $$
>
> **Check.** $y'' = 6x$, $y(0) = 0$, $y(1) = 2$. ✓ Uniqueness is guaranteed because $y'' = 0$ with $y(0) = y(1) = 0$ has only $y \equiv 0$ (a straight line through two zeros).

### Green's function for a BVP

As for IVPs (§11.6), the solution of a uniquely solvable BVP with homogeneous boundary conditions can be written as an integral against a kernel.

> [!theorem] Green's function for $-y'' = f$, $y(0) = y(1) = 0$
> $$
> y(x) = \int_0^1 G(x, s)\,f(s)\,ds,\qquad
> G(x, s) = \begin{cases} x(1 - s), & 0 \le x \le s,\\ s(1 - x), & s \le x \le 1.\end{cases}
> $$

*Construction.* For fixed $s$, $G$ must satisfy $-G'' = 0$ for $x \ne s$ (so it is linear on each side), vanish at $x = 0$ and $x = 1$, be continuous at $x = s$, and have a slope that drops by $1$ at $x = s$ (the response to a unit point source, since $-G'' = \delta(x - s)$). The linear pieces $Ax$ and $B(1 - x)$ with $As = B(1 - s)$ and $-B - A = -1$ give $A = 1 - s$, $B = s$. $\square$

$G$ is **symmetric**, $G(x, s) = G(s, x)$ — the reciprocity principle of electrostatics and of passive networks: the potential at $x$ due to a unit charge at $s$ equals the potential at $s$ due to a unit charge at $x$.

> [!example] Worked Example 11.34 — Using the Green's function
> Solve $-y'' = 1$, $y(0) = y(1) = 0$ (uniform space charge between grounded plates).
>
> **Step 1.** Split the integral at $s = x$:
> $$
> y(x) = \int_0^x s(1 - x)\,ds + \int_x^1 x(1 - s)\,ds = (1 - x)\frac{x^2}{2} + x\frac{(1 - x)^2}{2} = \frac{x(1 - x)}{2}.
> $$
>
> **Check.** $y'' = -1$ ✓, $y(0) = y(1) = 0$ ✓. The maximum $1/8$ occurs midway between the plates, as symmetry demands.

### Eigenvalue problems and Sturm–Liouville theory

A homogeneous BVP containing a parameter $\lambda$ normally has only the trivial solution; for special values of $\lambda$ — **eigenvalues** — it has nontrivial solutions, the **eigenfunctions**. This is the infinite-dimensional analogue of $\mathbf{A}\mathbf{v} = \lambda\mathbf{v}$ ([Ch. 10](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md)).

> [!theorem] The model eigenvalue problem
> $$
> y'' + \lambda y = 0,\qquad y(0) = 0,\quad y(L) = 0
> $$
> has nontrivial solutions only for
> $$
> \lambda_n = \left(\frac{n\pi}{L}\right)^2,\qquad y_n(x) = \sin\frac{n\pi x}{L},\qquad n = 1, 2, 3, \dots
> $$

*Derivation.* Consider the sign of $\lambda$.
- $\lambda = -k^2 < 0$: $y = c_1\cosh kx + c_2\sinh kx$. Then $y(0) = c_1 = 0$ and $y(L) = c_2\sinh kL = 0$ forces $c_2 = 0$, since $\sinh kL \ne 0$.
- $\lambda = 0$: $y = c_1 + c_2x$, and both conditions give $c_1 = c_2 = 0$.
- $\lambda = k^2 > 0$: $y = c_1\cos kx + c_2\sin kx$. $y(0) = 0$ gives $c_1 = 0$; $y(L) = c_2\sin kL = 0$ with $c_2 \ne 0$ requires $kL = n\pi$. $\square$

```desmos-graph
left=0; right=1; bottom=-1.3; top=1.3
---
y=\sin(\pi x)|blue
y=\sin(2\pi x)|red
y=\sin(3\pi x)|green
```
*The first three eigenfunctions $\sin(n\pi x/L)$ for $L = 1$: the $n$th has $n - 1$ interior zeros. They are the standing-wave shapes of a string fixed at both ends and the voltage patterns on a lossless line section shorted at both ends.*

The model problem is one member of a large family with the same properties.

> [!definition] Regular Sturm–Liouville problem
> $$
> \big(P(x)\,y'\big)' + Q(x)\,y + \lambda\,w(x)\,y = 0,\qquad a < x < b,
> $$
> $$
> \alpha_1y(a) + \alpha_2y'(a) = 0,\qquad \beta_1y(b) + \beta_2y'(b) = 0,
> $$
> with $P, P', Q, w$ continuous on $[a, b]$, $P > 0$ and $w > 0$ there, and each pair $(\alpha_1, \alpha_2)$, $(\beta_1, \beta_2)$ not both zero. The function $w$ is the **weight**.

Any equation $y'' + p(x)y' + q(x)y + \lambda s(x)y = 0$ can be put in this form by multiplying by $P = e^{\int p\,dx}$: then $Py'' + Py' p = (Py')'$, and $Q = Pq$, $w = Ps$.

> [!theorem] Properties of regular Sturm–Liouville problems
> 1. All eigenvalues are **real**, and they form an increasing sequence $\lambda_1 < \lambda_2 < \cdots \to \infty$.
> 2. Each eigenvalue is **simple**: its eigenfunctions are multiples of one function $y_n$, which can be taken real.
> 3. Eigenfunctions for different eigenvalues are **orthogonal with weight $w$**:
> $$
> \langle y_m, y_n\rangle_w = \int_a^b y_m(x)\,y_n(x)\,w(x)\,dx = 0,\qquad m \ne n.
> $$
> 4. The $n$th eigenfunction has exactly $n - 1$ zeros in $(a, b)$.
> 5. The eigenfunctions are **complete**: any square-integrable $f$ on $[a, b]$ has a convergent (in the mean-square sense) expansion $f = \sum c_ny_n$ with $c_n = \langle f, y_n\rangle_w/\langle y_n, y_n\rangle_w$.

*Proof of orthogonality (and reality).* Write $\mathcal{S}[y] = (Py')' + Qy$. For any two functions satisfying the boundary conditions, **Lagrange's identity**

$$
u\,\mathcal{S}[v] - v\,\mathcal{S}[u] = \big(P(uv' - vu')\big)'
$$

(verify by expanding) integrates to

$$
\int_a^b\big(u\,\mathcal{S}[v] - v\,\mathcal{S}[u]\big)\,dx = \Big[P(uv' - vu')\Big]_a^b = 0,
$$

because at each end the boundary condition makes the vectors $(u, u')$ and $(v, v')$ both orthogonal to the nonzero vector $(\alpha_1, \alpha_2)$ (resp. $(\beta_1, \beta_2)$), hence parallel, so $uv' - vu' = 0$ there. The operator is **self-adjoint** — the function-space counterpart of a symmetric matrix. Now take $u = y_m$, $v = y_n$ with $\mathcal{S}[y_k] = -\lambda_kwy_k$: the identity gives $(\lambda_m - \lambda_n)\int_a^b y_my_nw\,dx = 0$, so the integral vanishes when $\lambda_m \ne \lambda_n$. Applying the same identity to $u = y$, $v = \bar y$ for a possibly complex eigenpair gives $(\lambda - \bar\lambda)\int\lvert y\rvert^2w\,dx = 0$, so $\lambda$ is real. $\square$

Parts 1, 2, 4 and 5 need more machinery (comparison theorems and functional analysis) and are quoted. Part 3 is what makes coefficients of an eigenfunction expansion computable one at a time — exactly as with the orthogonal eigenvectors of a symmetric matrix.

> [!example] Worked Example 11.35 — Neumann conditions
> Find the eigenvalues and eigenfunctions of $y'' + \lambda y = 0$, $y'(0) = 0$, $y'(L) = 0$.
>
> **Step 1 ($\lambda < 0$).** $y = c_1\cosh kx + c_2\sinh kx$: $y'(0) = kc_2 = 0$, then $y'(L) = kc_1\sinh kL = 0$ forces $c_1 = 0$. No eigenvalues.
>
> **Step 2 ($\lambda = 0$).** $y = c_1 + c_2x$: $y' = c_2 = 0$, leaving $y = c_1$. So $\lambda_0 = 0$ **is** an eigenvalue with eigenfunction $y_0 = 1$.
>
> **Step 3 ($\lambda = k^2 > 0$).** $y'(0) = kc_2 = 0$; $y'(L) = -kc_1\sin kL = 0$ needs $kL = n\pi$.
>
> **Result.** $\lambda_n = (n\pi/L)^2$, $y_n = \cos(n\pi x/L)$, $n = 0, 1, 2, \dots$ — the cosine set of Fourier series. The constant mode $n = 0$ is the "DC term"; it exists because insulating (zero-flux) ends allow any uniform level.

> [!example] Worked Example 11.36 — A Robin condition with transcendental eigenvalues
> Find the eigenvalues of $y'' + \lambda y = 0$, $y(0) = 0$, $y'(1) + y(1) = 0$.
>
> **Step 1.** For $\lambda \le 0$ one checks as before that only $y \equiv 0$ survives (for $\lambda = -k^2$, the condition becomes $k\cosh k + \sinh k = 0$, impossible for $k > 0$; for $\lambda = 0$, $y = c_2x$ gives $2c_2 = 0$).
>
> **Step 2.** For $\lambda = k^2 > 0$, $y(0) = 0$ gives $y = \sin kx$, and the Robin condition gives $k\cos k + \sin k = 0$, i.e.
> $$
> \tan k = -k.
> $$
>
> **Step 3 (solve numerically).** The roots are the intersections of $y = \tan k$ with $y = -k$, one in each interval $\big((n - \tfrac12)\pi, n\pi\big)$. Newton's method ([Ch. 22](../part-7-numerical-computational/22-numerical-methods.md)) gives
> $$
> k_1 = 2.0288,\quad k_2 = 4.9132,\quad k_3 = 7.9787,\qquad \lambda_n = k_n^2 = 4.116,\ 24.14,\ 63.66.
> $$
>
> **Check.** $\tan(2.0288) = -2.0288$ to four figures. ✓ For large $n$, $k_n \to (n - \tfrac12)\pi$ from above: $(3 - 0.5)\pi = 7.854$ versus $7.979$. The eigenfunctions $\sin k_nx$ are orthogonal on $[0, 1]$ even though the $k_n$ are not equally spaced — guaranteed by the Sturm–Liouville theorem, not by any trigonometric identity.

```desmos-graph
left=0; right=9; bottom=-10; top=10
---
y=\tan(x)|blue
y=-x|red
(2.0288,-2.0288)|label:k₁
(4.9132,-4.9132)|label:k₂
(7.9787,-7.9787)|label:k₃
```
*Graphical solution of $\tan k = -k$ (Example 11.36): the eigenvalues of a Robin problem are the intersections of $\tan k$ (blue) with the line $-k$ (red), approaching $(n - \tfrac12)\pi$ as $n$ grows.*

> [!example] Worked Example 11.37 — A weighted Sturm–Liouville problem
> Solve $x^2y'' + xy' + \lambda y = 0$ on $1 < x < e$ with $y(1) = y(e) = 0$, and identify the weight.
>
> **Step 1 (SL form).** Dividing by $x$: $xy'' + y' + (\lambda/x)y = 0$, i.e. $(xy')' + \lambda\,\dfrac1x\,y = 0$. So $P = x$, $Q = 0$, $w = 1/x$.
>
> **Step 2 (solve).** This is Cauchy–Euler (§11.4) with indicial equation $m^2 + \lambda = 0$. For $\lambda = k^2 > 0$, $y = c_1\cos(k\ln x) + c_2\sin(k\ln x)$. $y(1) = c_1 = 0$; $y(e) = c_2\sin k = 0$ needs $k = n\pi$. (As in the model problem, $\lambda \le 0$ gives only the trivial solution.)
>
> **Step 3.** $\lambda_n = n^2\pi^2$, $y_n = \sin(n\pi\ln x)$.
>
> **Check (orthogonality with weight).** With $u = \ln x$, $du = dx/x$:
> $$
> \int_1^e\sin(m\pi\ln x)\sin(n\pi\ln x)\,\frac{dx}{x} = \int_0^1\sin(m\pi u)\sin(n\pi u)\,du = 0,\quad m \ne n. ✓
> $$
> Without the weight $1/x$ the integral is *not* zero — orthogonality is always relative to the weight of the SL form.

> [!warning] Common pitfall — forgetting $\lambda \le 0$
> Always check $\lambda < 0$, $\lambda = 0$ and $\lambda > 0$ separately. Dirichlet problems have no zero eigenvalue, but Neumann and periodic problems do, and some Robin conditions (with the "wrong" sign) produce a negative eigenvalue. Missing the $\lambda_0 = 0$ mode loses the average value of every Fourier-cosine expansion.

> [!tip] Engineering note — eigenvalues are resonances
> In a separation-of-variables solution of the wave equation on a line or in a cavity, $\lambda_n = (n\pi/L)^2$ becomes $(\omega_n/c)^2$: the eigenvalues fix the **resonant frequencies** $\omega_n = n\pi c/L$, and the eigenfunctions are the **mode shapes**. The Fredholm alternative is the precise statement of what happens if a structure is driven exactly at one of them: no steady solution exists unless the drive is orthogonal to that mode. [Ch. 13](13-partial-differential-equations.md) uses these results directly.

---

## 11.9 Systems of differential equations

### Motivation

Real networks have several state variables that interact: two coupled resonators, a multi-loop filter, a motor with electrical and mechanical dynamics, a transistor amplifier with several internal capacitances. Writing each state's rate of change in terms of all states gives a **system** of first-order equations. Systems are also the universal form: every $n$th-order equation is one (§11.5), every simulator integrates one, and the geometry of systems — trajectories in a state space, equilibria, modes — is the clearest way to understand dynamics. This section gives the general theory and the eigenvector method for constant coefficients; the matrix exponential and the full state-variable formulation of circuits are developed in [Ch. 12](12-linear-differential-equations.md).

> [!definition] Linear first-order system
> $$
> \dot{\mathbf{x}}(t) = \mathbf{A}(t)\,\mathbf{x}(t) + \mathbf{f}(t),\qquad \mathbf{x}(t) \in \mathbb{R}^n,
> $$
> with $\mathbf{A}(t)$ an $n \times n$ matrix and $\mathbf{f}(t)$ an input vector, both continuous on $I$. The system is **homogeneous** if $\mathbf{f} \equiv \mathbf{0}$. An IVP adds $\mathbf{x}(t_0) = \mathbf{x}_0$. By §11.7, the IVP has a unique solution on all of $I$.

### Structure theory

Everything in §11.4–§11.6 has a matrix counterpart.

- **Dimension.** Solutions of $\dot{\mathbf{x}} = \mathbf{A}\mathbf{x}$ form an $n$-dimensional vector space.
- **Fundamental matrix.** A matrix $\boldsymbol{\Phi}(t)$ whose $n$ columns are independent solutions satisfies $\dot{\boldsymbol{\Phi}} = \mathbf{A}\boldsymbol{\Phi}$ and $\det\boldsymbol{\Phi} \ne 0$. The general homogeneous solution is $\mathbf{x} = \boldsymbol{\Phi}(t)\mathbf{c}$, and the IVP solution is $\mathbf{x}(t) = \boldsymbol{\Phi}(t)\boldsymbol{\Phi}(t_0)^{-1}\mathbf{x}_0$.
- **Liouville's formula** (Abel's formula for systems):
$$
\det\boldsymbol{\Phi}(t) = \det\boldsymbol{\Phi}(t_0)\,\exp\!\left(\int_{t_0}^{t}\operatorname{tr}\mathbf{A}(s)\,ds\right).
$$
For the companion matrix of §11.5, $\operatorname{tr}\mathbf{A} = -p_{n-1}$ and $\boldsymbol{\Phi}$ is the Wronskian matrix, so this reproduces Abel's formula.
- **Variation of parameters.** Seeking $\mathbf{x}_p = \boldsymbol{\Phi}(t)\mathbf{u}(t)$ gives $\boldsymbol{\Phi}\dot{\mathbf{u}} = \mathbf{f}$, so
$$
\mathbf{x}(t) = \boldsymbol{\Phi}(t)\boldsymbol{\Phi}(t_0)^{-1}\mathbf{x}_0 + \boldsymbol{\Phi}(t)\int_{t_0}^{t}\boldsymbol{\Phi}(s)^{-1}\mathbf{f}(s)\,ds.
$$
The first term is the zero-input response, the second the zero-state response — the matrix version of §11.3's formula.

*Derivation of variation of parameters.* $\dot{\mathbf{x}}_p = \dot{\boldsymbol{\Phi}}\mathbf{u} + \boldsymbol{\Phi}\dot{\mathbf{u}} = \mathbf{A}\boldsymbol{\Phi}\mathbf{u} + \boldsymbol{\Phi}\dot{\mathbf{u}}$. Setting this equal to $\mathbf{A}\mathbf{x}_p + \mathbf{f} = \mathbf{A}\boldsymbol{\Phi}\mathbf{u} + \mathbf{f}$ leaves $\boldsymbol{\Phi}\dot{\mathbf{u}} = \mathbf{f}$; $\boldsymbol{\Phi}$ is invertible, so $\dot{\mathbf{u}} = \boldsymbol{\Phi}^{-1}\mathbf{f}$. $\square$

### Constant coefficients: the eigenvector method

For constant $\mathbf{A}$, try $\mathbf{x} = \mathbf{v}e^{\lambda t}$ with a constant vector $\mathbf{v} \ne \mathbf{0}$. Then $\dot{\mathbf{x}} = \lambda\mathbf{v}e^{\lambda t}$ and $\mathbf{A}\mathbf{x} = \mathbf{A}\mathbf{v}e^{\lambda t}$, so we need

$$
\mathbf{A}\mathbf{v} = \lambda\mathbf{v}:
$$

an eigenvalue problem ([Ch. 10](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md)). Each eigenpair gives a **mode** — a solution that keeps its shape $\mathbf{v}$ and only scales by $e^{\lambda t}$.

> [!theorem] Eigen-solution of $\dot{\mathbf{x}} = \mathbf{A}\mathbf{x}$
> If $\mathbf{A}$ has $n$ linearly independent eigenvectors $\mathbf{v}_1, \dots, \mathbf{v}_n$ with eigenvalues $\lambda_1, \dots, \lambda_n$, then
> $$
> \mathbf{x}(t) = c_1\mathbf{v}_1e^{\lambda_1t} + \cdots + c_n\mathbf{v}_ne^{\lambda_nt},
> $$
> and the constants are found from $\mathbf{x}(0) = \sum c_k\mathbf{v}_k$, i.e. by solving $\mathbf{V}\mathbf{c} = \mathbf{x}_0$ with $\mathbf{V} = [\mathbf{v}_1 \cdots \mathbf{v}_n]$.

Independence of the eigensolutions follows from Liouville or directly: at $t = 0$ the fundamental matrix is $\mathbf{V}$, which is invertible. Three cases arise.

**Complex eigenvalues.** If $\mathbf{A}$ is real and $\lambda = \alpha + j\beta$ with eigenvector $\mathbf{v} = \mathbf{a} + j\mathbf{b}$, then $\bar\lambda, \bar{\mathbf{v}}$ is also an eigenpair, and the real and imaginary parts of $\mathbf{v}e^{\lambda t}$ are two real independent solutions:

$$
e^{\alpha t}\big(\mathbf{a}\cos\beta t - \mathbf{b}\sin\beta t\big),\qquad e^{\alpha t}\big(\mathbf{a}\sin\beta t + \mathbf{b}\cos\beta t\big).
$$

(Real and imaginary parts of a complex solution of a real linear system are solutions, because $\mathbf{A}$ is real.)

**Repeated eigenvalues with too few eigenvectors (defective $\mathbf{A}$).** If $\lambda$ is a double eigenvalue with only one eigenvector $\mathbf{v}$, a second solution is

$$
\mathbf{x}_2 = (t\,\mathbf{v} + \mathbf{w})\,e^{\lambda t},\qquad (\mathbf{A} - \lambda\mathbf{I})\mathbf{w} = \mathbf{v},
$$

where $\mathbf{w}$ is a **generalised eigenvector**. (Substituting: $\dot{\mathbf{x}}_2 = (\mathbf{v} + \lambda t\mathbf{v} + \lambda\mathbf{w})e^{\lambda t}$, and $\mathbf{A}\mathbf{x}_2 = (\lambda t\mathbf{v} + \lambda\mathbf{w} + \mathbf{v})e^{\lambda t}$.) This is the matrix origin of the $te^{\lambda t}$ terms of repeated characteristic roots. The general treatment via $e^{\mathbf{A}t}$ and Jordan form is in [Ch. 12](12-linear-differential-equations.md).

> [!example] Worked Example 11.38 — A symmetric $2 \times 2$ system
> Solve $\dot x_1 = -3x_1 + x_2$, $\dot x_2 = x_1 - 3x_2$ with $x_1(0) = 2$, $x_2(0) = 0$.
>
> **Step 1 (eigenvalues).** $\mathbf{A} = \begin{pmatrix}-3 & 1\\ 1 & -3\end{pmatrix}$, $\det(\mathbf{A} - \lambda\mathbf{I}) = (\lambda + 3)^2 - 1 = 0$, so $\lambda_1 = -2$, $\lambda_2 = -4$.
>
> **Step 2 (eigenvectors).** For $\lambda_1 = -2$: $\begin{pmatrix}-1 & 1\\ 1 & -1\end{pmatrix}\mathbf{v} = \mathbf{0}$ gives $\mathbf{v}_1 = (1, 1)^\mathsf{T}$. For $\lambda_2 = -4$: $\mathbf{v}_2 = (1, -1)^\mathsf{T}$. (Orthogonal, as for any symmetric matrix.)
>
> **Step 3 (fit).** $c_1 + c_2 = 2$, $c_1 - c_2 = 0$, so $c_1 = c_2 = 1$:
> $$
> x_1 = e^{-2t} + e^{-4t},\qquad x_2 = e^{-2t} - e^{-4t}.
> $$
>
> **Check.** $\dot x_1 = -2e^{-2t} - 4e^{-4t}$ and $-3x_1 + x_2 = -2e^{-2t} - 4e^{-4t}$. ✓ Liouville: $\operatorname{tr}\mathbf{A} = -6$, and $\det\boldsymbol{\Phi} \propto e^{-2t}e^{-4t} = e^{-6t}$. ✓
>
> **Modes.** The symmetric mode $(1, 1)$ (both states equal) decays with time constant $0.5$; the antisymmetric mode $(1, -1)$ decays twice as fast. Any initial state is a mixture of the two.

```desmos-graph
left=-0.1; right=3; bottom=-0.1; top=2.1
---
y=e^{-2x}+e^{-4x}|blue|x>0
y=e^{-2x}-e^{-4x}|red|x>0
```
*Example 11.38: $x_1$ (blue) and $x_2$ (red). The fast mode $e^{-4t}$ dies first, after which $x_1 \approx x_2$ and both follow the slow mode $e^{-2t}$.*

> [!example] Worked Example 11.39 — Complex eigenvalues
> Solve $\dot{\mathbf{x}} = \begin{pmatrix}-1 & -4\\ 1 & -1\end{pmatrix}\mathbf{x}$, $\mathbf{x}(0) = (1, 0)^\mathsf{T}$.
>
> **Step 1.** $\det(\mathbf{A} - \lambda\mathbf{I}) = (\lambda + 1)^2 + 4 = 0$, so $\lambda = -1 \pm 2j$.
>
> **Step 2.** For $\lambda = -1 + 2j$: $(\mathbf{A} - \lambda\mathbf{I}) = \begin{pmatrix}-2j & -4\\ 1 & -2j\end{pmatrix}$; the second row gives $v_1 = 2jv_2$, so $\mathbf{v} = (2j, 1)^\mathsf{T}$, i.e. $\mathbf{a} = (0, 1)^\mathsf{T}$, $\mathbf{b} = (2, 0)^\mathsf{T}$.
>
> **Step 3 (real solutions).**
> $$
> \mathbf{x} = e^{-t}\left[c_1\begin{pmatrix}-2\sin2t\\ \cos2t\end{pmatrix} + c_2\begin{pmatrix}2\cos2t\\ \sin2t\end{pmatrix}\right].
> $$
>
> **Step 4 (fit).** At $t = 0$: $(2c_2, c_1) = (1, 0)$, so $c_1 = 0$, $c_2 = \tfrac12$:
> $$
> x_1 = e^{-t}\cos2t,\qquad x_2 = \tfrac12e^{-t}\sin2t.
> $$
>
> **Check.** $\dot x_1 = -e^{-t}\cos2t - 2e^{-t}\sin2t = -x_1 - 4x_2$ ✓; $\dot x_2 = -\tfrac12e^{-t}\sin2t + e^{-t}\cos2t = x_1 - x_2$ ✓.
>
> In the phase plane the state spirals into the origin (a stable focus), turning at $2\ \mathrm{rad/s}$ and shrinking by $e^{-1}$ per second — the behaviour of an underdamped RLC tank.

```desmos-graph
left=-0.1; right=5; bottom=-1.1; top=1.1
---
y=e^{-x}\cos(2x)|blue|x>0
y=0.5e^{-x}\sin(2x)|red|x>0
y=e^{-x}|grey|dashed|x>0
y=-e^{-x}|grey|dashed|x>0
```
*Example 11.39: $x_1 = e^{-t}\cos2t$ (blue) and $x_2 = \tfrac12e^{-t}\sin2t$ (red), a quarter-period apart inside the decaying envelope $\pm e^{-t}$ set by $\operatorname{Re}\lambda = -1$.*

> [!example] Worked Example 11.40 — A defective matrix
> Solve $\dot{\mathbf{x}} = \begin{pmatrix}1 & -1\\ 1 & 3\end{pmatrix}\mathbf{x}$, $\mathbf{x}(0) = (1, 1)^\mathsf{T}$.
>
> **Step 1.** $\operatorname{tr}\mathbf{A} = 4$, $\det\mathbf{A} = 4$, so $\lambda^2 - 4\lambda + 4 = (\lambda - 2)^2$: a double eigenvalue $\lambda = 2$.
>
> **Step 2.** $\mathbf{A} - 2\mathbf{I} = \begin{pmatrix}-1 & -1\\ 1 & 1\end{pmatrix}$ has rank 1, so there is only one eigenvector, $\mathbf{v} = (1, -1)^\mathsf{T}$.
>
> **Step 3 (generalised eigenvector).** $(\mathbf{A} - 2\mathbf{I})\mathbf{w} = \mathbf{v}$: $-w_1 - w_2 = 1$; take $\mathbf{w} = (-1, 0)^\mathsf{T}$.
>
> **Step 4.** $\mathbf{x} = c_1\mathbf{v}e^{2t} + c_2(t\mathbf{v} + \mathbf{w})e^{2t}$. At $t = 0$: $c_1\mathbf{v} + c_2\mathbf{w} = (c_1 - c_2, -c_1) = (1, 1)$, so $c_1 = -1$, $c_2 = -2$:
> $$
> \mathbf{x}(t) = e^{2t}\begin{pmatrix}1 - 2t\\ 1 + 2t\end{pmatrix}.
> $$
>
> **Check.** $\dot x_1 = e^{2t}(2 - 4t - 2) = -4te^{2t}$ and $x_1 - x_2 = -4te^{2t}$ ✓; $\dot x_2 = e^{2t}(2 + 4t + 2) = (4 + 4t)e^{2t}$ and $x_1 + 3x_2 = (4 + 4t)e^{2t}$ ✓.

### The phase plane for $2 \times 2$ systems

Plotting $x_2$ against $x_1$ — the **phase plane** — shows all trajectories at once. For $\dot{\mathbf{x}} = \mathbf{A}\mathbf{x}$ with $\det\mathbf{A} \ne 0$ the origin is the only equilibrium, and its type depends only on the eigenvalues, which in turn depend only on $T = \operatorname{tr}\mathbf{A}$ and $\Delta = \det\mathbf{A}$ through $\lambda^2 - T\lambda + \Delta = 0$:

| Condition | Eigenvalues | Equilibrium type | Stability |
|---|---|---|---|
| $\Delta < 0$ | real, opposite signs | saddle | unstable |
| $\Delta > 0$, $T^2 > 4\Delta$, $T < 0$ | real, both negative | node | asymptotically stable |
| $\Delta > 0$, $T^2 > 4\Delta$, $T > 0$ | real, both positive | node | unstable |
| $\Delta > 0$, $T^2 < 4\Delta$, $T < 0$ | complex, $\operatorname{Re}\lambda < 0$ | spiral (focus) | asymptotically stable |
| $\Delta > 0$, $T^2 < 4\Delta$, $T > 0$ | complex, $\operatorname{Re}\lambda > 0$ | spiral | unstable |
| $\Delta > 0$, $T = 0$ | purely imaginary | centre | stable, not asymptotically |

So a $2 \times 2$ linear system is asymptotically stable exactly when $\operatorname{tr}\mathbf{A} < 0$ and $\det\mathbf{A} > 0$ — the two-state case of the Routh–Hurwitz criterion. Examples 11.38, 11.39 and 11.40 are a stable node ($T = -6$, $\Delta = 8$), a stable spiral ($T = -2$, $\Delta = 5$) and an unstable degenerate node ($T = 4$, $\Delta = 4 = T^2/4$).

### Nonlinear systems and linearisation

For an autonomous nonlinear system $\dot{\mathbf{x}} = \mathbf{F}(\mathbf{x})$, an **equilibrium** $\mathbf{x}^*$ satisfies $\mathbf{F}(\mathbf{x}^*) = \mathbf{0}$. Setting $\mathbf{x} = \mathbf{x}^* + \boldsymbol{\eta}$ and expanding to first order ([Ch. 6](../part-2-calculus/06-multivariable-calculus.md)),

$$
\dot{\boldsymbol{\eta}} \approx \mathbf{J}\,\boldsymbol{\eta},\qquad \mathbf{J} = \left.\frac{\partial\mathbf{F}}{\partial\mathbf{x}}\right|_{\mathbf{x}^*},
$$

where $\mathbf{J}$ is the Jacobian matrix.

> [!theorem] Linearisation (Hartman–Grobman, informal)
> If no eigenvalue of $\mathbf{J}$ has zero real part (the equilibrium is **hyperbolic**), then near $\mathbf{x}^*$ the trajectories of the nonlinear system look like those of $\dot{\boldsymbol{\eta}} = \mathbf{J}\boldsymbol{\eta}$: same type (node, saddle, spiral) and same stability. If some eigenvalue has zero real part (for example a linear centre), the linearisation is inconclusive and the nonlinear terms decide.

> [!example] Worked Example 11.41 — Equilibria of a damped pendulum
> The damped pendulum $\ddot\theta + 0.5\dot\theta + \sin\theta = 0$ becomes, with $x_1 = \theta$, $x_2 = \dot\theta$,
> $$
> \dot x_1 = x_2,\qquad \dot x_2 = -\sin x_1 - 0.5x_2.
> $$
> Classify the equilibria $(0, 0)$ and $(\pi, 0)$.
>
> **Step 1 (Jacobian).** $\mathbf{J} = \begin{pmatrix}0 & 1\\ -\cos x_1 & -0.5\end{pmatrix}$.
>
> **Step 2 (hanging down, $(0, 0)$).** $\mathbf{J} = \begin{pmatrix}0 & 1\\ -1 & -0.5\end{pmatrix}$: $T = -0.5$, $\Delta = 1$, $T^2 < 4\Delta$. Eigenvalues $-0.25 \pm 0.968j$: a **stable spiral** — the pendulum swings with decaying amplitude.
>
> **Step 3 (inverted, $(\pi, 0)$).** $\mathbf{J} = \begin{pmatrix}0 & 1\\ 1 & -0.5\end{pmatrix}$: $\Delta = -1 < 0$, a **saddle**. Eigenvalues $\big(-0.5 \pm \sqrt{4.25}\big)/2 = 0.781,\ -1.281$. The unstable direction is the eigenvector for $0.781$: almost every small disturbance topples the pendulum.
>
> **Check.** Both equilibria are hyperbolic, so the nonlinear pendulum has the same local behaviour. The same analysis applied to a phase-locked loop (whose phase-error equation has exactly this form) shows that the loop locks at one equilibrium and is repelled from the other.

> [!tip] Engineering note — small-signal analysis is linearisation
> The small-signal model of a transistor circuit is the Jacobian of the large-signal equations at the DC operating point: $g_m = \partial i_C/\partial v_{BE}$, $r_\pi$, $r_o$ are entries of $\mathbf{J}$. Hartman–Grobman is the mathematical warrant for trusting small-signal stability analysis — and the warning that when a pole sits on the $j\omega$-axis (an oscillator at the edge of start-up), the linear model cannot predict the amplitude; nonlinear terms must be kept (compare Example 11.13).

---

## 11.10 Applications to electrical circuits

### From a circuit to a differential equation

Every lumped-circuit ODE is built from two ingredients.

1. **Element laws.** $v_R = Ri$; $v_L = L\,\dfrac{di}{dt}$; $i_C = C\,\dfrac{dv}{dt}$; plus any nonlinear laws such as a diode $i = I_S\big(e^{v/V_T} - 1\big)$, a saturating inductor $v = \dfrac{d\lambda(i)}{dt}$, or a time-varying element $i = \dfrac{d}{dt}\big(C(t)v\big)$.
2. **Kirchhoff's laws.** KVL around loops and KCL at nodes — linear algebraic constraints ([Ch. 8](../part-3-linear-algebra/08-matrices-and-linear-systems.md)).

The natural unknowns are the **state variables**: capacitor voltages and inductor currents. They carry the stored energy ($\tfrac12Cv^2$, $\tfrac12Li^2$), they are continuous in time (§11.7), and their values at one instant, together with the sources thereafter, determine everything else. The number of independent state variables is the **order** of the circuit.

**A systematic procedure.**

1. Choose the state variables: one $v_C$ per independent capacitor, one $i_L$ per independent inductor.
2. Express each $C\,\dot v_C$ as a capacitor current and each $L\,\dot i_L$ as an inductor voltage.
3. Use KVL, KCL and the resistive element laws to write those currents and voltages in terms of the states and the sources.
4. The result is $\dot{\mathbf{x}} = \mathbf{F}(\mathbf{x}, t)$ — linear, $\dot{\mathbf{x}} = \mathbf{A}\mathbf{x} + \mathbf{B}\mathbf{u}(t)$, if all elements are linear and time-invariant.
5. Eliminate states if a single higher-order equation for one output is wanted (§11.5 in reverse).

[Ch. 12](12-linear-differential-equations.md) refines steps 1–4 into the formal state-variable method and solves the resulting constant-coefficient systems completely. The table below shows which method of *this* chapter each common circuit leads to.

| Circuit | Equation | Type (this chapter) |
|---|---|---|
| Series RC, source $v_s$ | $RC\,\dot v_C + v_C = v_s(t)$ | linear first order (§11.3) |
| Series RL | $L\,\dot i + Ri = v_s(t)$ | linear first order (§11.3) |
| $C$ discharging into nonlinear load $i = g(v)$ | $C\,\dot v = -g(v)$ | separable (§11.2) |
| RC with time-varying $C(t)$ or $R(t)$ | $R\,\dfrac{d}{dt}(Cv) + v = v_s$ | linear, variable coefficient (§11.3) |
| Series RLC | $LC\,\ddot v_C + RC\,\dot v_C + v_C = v_s$ | linear second order (§11.4, Ch. 12) |
| Multi-loop or multi-node network | $\dot{\mathbf{x}} = \mathbf{A}\mathbf{x} + \mathbf{B}\mathbf{u}$ | linear system (§11.9) |
| Negative-resistance oscillator | $\ddot v - \mu(1 - v^2)\dot v + v = 0$ | nonlinear system (§11.9) |
| Distributed line or electrode at DC | $V'' = rg\,V$ on $0 < x < \ell$ | boundary-value problem (§11.8) |

> [!example] Worked Example 11.42 — RC charging
> A $10\ \mathrm{k\Omega}$ resistor and a $10\ \mu\mathrm F$ capacitor (initially uncharged) are connected to a $5\ \mathrm V$ source at $t = 0$. Find $v_C(t)$, its value at $t = 0.2\ \mathrm s$, and the time to reach $4.5\ \mathrm V$.
>
> **Step 1.** KVL: $RC\,\dot v_C + v_C = 5$, with $\tau = RC = 10^4\cdot10^{-5} = 0.1\ \mathrm s$.
>
> **Step 2.** By §11.3 with $v_C(0) = 0$:
> $$
> v_C(t) = 5\left(1 - e^{-t/0.1}\right)\ \mathrm V.
> $$
>
> **Step 3.** At $t = 0.2\ \mathrm s$: $v_C = 5(1 - e^{-2}) = 5(0.8647) = 4.32\ \mathrm V$.
>
> **Step 4.** $5(1 - e^{-t/\tau}) = 4.5$ gives $e^{-t/\tau} = 0.1$, so $t = \tau\ln10 = 0.230\ \mathrm s$.
>
> **Check.** Units: $\Omega\cdot\mathrm F = (\mathrm{V/A})(\mathrm{A\,s/V}) = \mathrm s$. ✓ Initial slope $\dot v_C(0) = 5/\tau = 50\ \mathrm{V/s}$ equals $i(0)/C = (5/10^4)/10^{-5} = 50\ \mathrm{V/s}$. ✓

```desmos-graph
left=-0.02; right=0.6; bottom=-0.3; top=5.8
---
y=5\left(1-e^{-x/0.1}\right)|blue|x>0
y=5|black|dashed
y=50x|grey|dotted|0<x<0.1
(0.1,3.161)|label:t = τ, 63%
(0.2,4.323)|label:t = 2τ, 86%
```
*RC charging with $\tau = 0.1$ s: 63% of the final value after one time constant, 86% after two and over 98% after $4\tau$. The dotted initial tangent reaches the final value exactly at $t = \tau$.*

> [!example] Worked Example 11.43 — A two-loop network as a system
> Loop 1 contains the source $v(t)$, $R_1 = 2\ \Omega$ and $L_1 = 1\ \mathrm H$. Loop 2 contains $R_3 = 2\ \Omega$ and $L_2 = 1\ \mathrm H$. The loops share $R_2 = 1\ \Omega$. Find the mesh currents when $v = 8\ \mathrm V$ is switched on at $t = 0$ with both inductors initially unenergised.
>
> **Step 1 (KVL).** With mesh currents $i_1, i_2$ (which here are also the inductor currents, hence states):
> $$
> \begin{aligned}
> L_1\dot i_1 + R_1i_1 + R_2(i_1 - i_2) &= v(t),\\
> L_2\dot i_2 + R_3i_2 + R_2(i_2 - i_1) &= 0.
> \end{aligned}
> $$
>
> **Step 2 (state form).**
> $$
> \frac{d}{dt}\begin{pmatrix}i_1\\ i_2\end{pmatrix} = \begin{pmatrix}-3 & 1\\ 1 & -3\end{pmatrix}\begin{pmatrix}i_1\\ i_2\end{pmatrix} + \begin{pmatrix}1\\ 0\end{pmatrix}v(t).
> $$
> The matrix is that of Example 11.38: eigenvalues $-2$ and $-4\ \mathrm{s^{-1}}$ (the network's **natural frequencies**), modes $(1, 1)$ and $(1, -1)$.
>
> **Step 3 (particular solution).** For constant input, try a constant: $\mathbf{0} = \mathbf{A}\mathbf{x}_p + \mathbf{b}\cdot8$, so $\mathbf{x}_p = -8\mathbf{A}^{-1}\mathbf{b}$. With $\mathbf{A}^{-1} = \tfrac18\begin{pmatrix}-3 & -1\\ -1 & -3\end{pmatrix}$, $\mathbf{x}_p = (3, 1)^\mathsf{T}$ A.
>
> **Step 4 (fit zero initial state).** $\mathbf{x} = (3, 1)^\mathsf{T} + c_1(1, 1)^\mathsf{T}e^{-2t} + c_2(1, -1)^\mathsf{T}e^{-4t}$, with $c_1 + c_2 = -3$ and $c_1 - c_2 = -1$, so $c_1 = -2$, $c_2 = -1$:
> $$
> i_1 = 3 - 2e^{-2t} - e^{-4t}\ \mathrm A,\qquad i_2 = 1 - 2e^{-2t} + e^{-4t}\ \mathrm A.
> $$
>
> **Check.** $i_1(0) = i_2(0) = 0$ ✓. DC check with inductors as shorts: loop 1, $2(3) + 1(3 - 1) = 8$ ✓; loop 2, $2(1) + 1(1 - 3) = 0$ ✓. Settling is governed by the slow mode $e^{-2t}$: within 1% of DC after about $t = \ln(200)/2 \approx 2.6\ \mathrm s$.

> [!example] Worked Example 11.44 — A leaky line at DC: a boundary-value problem
> A long cable (or a buried grounding conductor) has series resistance $r = 0.1\ \Omega/\mathrm m$ and leakage conductance $g = 10^{-5}\ \mathrm{S/m}$ to ground, length $\ell = 1\ \mathrm{km}$. A DC voltage $V_0$ is applied at $x = 0$ and the far end is open. Find $V(x)$, the far-end voltage and the input resistance.
>
> **Step 1 (model).** On a short piece $dx$: the series drop is $dV = -rI\,dx$ and the leakage is $dI = -gV\,dx$. Hence $V' = -rI$, $I' = -gV$, and eliminating $I$,
> $$
> V'' = \alpha^2V,\qquad \alpha = \sqrt{rg} = \sqrt{10^{-6}} = 10^{-3}\ \mathrm{m^{-1}}.
> $$
>
> **Step 2 (boundary conditions).** $V(0) = V_0$ (Dirichlet). The open end carries no current: $I(\ell) = -V'(\ell)/r = 0$, i.e. $V'(\ell) = 0$ (Neumann).
>
> **Step 3 (solve).** Using the basis $\cosh\alpha(\ell - x)$, $\sinh\alpha(\ell - x)$ — chosen so the Neumann condition at $x = \ell$ is easy — the condition $V'(\ell) = 0$ kills the $\sinh$ term, and $V(0) = V_0$ fixes the amplitude:
> $$
> V(x) = V_0\,\frac{\cosh\alpha(\ell - x)}{\cosh\alpha\ell}.
> $$
> The homogeneous problem ($V_0 = 0$) has only the trivial solution, since $\cosh$ never vanishes, so by the Fredholm alternative this is the unique solution.
>
> **Step 4 (numbers).** $\alpha\ell = 1$: far-end voltage $V(\ell) = V_0/\cosh1 = 0.648\,V_0$. Input current $I(0) = -V'(0)/r = V_0\alpha\tanh(\alpha\ell)/r$, so
> $$
> R_{\text{in}} = \frac{r}{\alpha\tanh\alpha\ell} = \frac{0.1}{10^{-3}\cdot0.7616} = 131\ \Omega.
> $$
>
> **Check.** Two limits confirm the formula. For a very long line ($\alpha\ell \to \infty$), $\tanh \to 1$ and $R_{\text{in}} \to \sqrt{r/g} = 100\ \Omega$, the characteristic resistance. For a very short line, $\tanh\alpha\ell \approx \alpha\ell$ and $R_{\text{in}} \approx 1/(g\ell)$, the total leakage resistance in parallel to ground, as expected when series drop is negligible. Our line is finite and open-ended, so less leakage path is available than on an infinite line and $R_{\text{in}} = 131\ \Omega$ exceeds $100\ \Omega$. ✓ Directly: $V'' = \alpha^2V$ holds since $\big(\cosh\alpha(\ell - x)\big)'' = \alpha^2\cosh\alpha(\ell - x)$, and $V'(\ell) = -\alpha V_0\sinh 0/\cosh\alpha\ell = 0$. ✓ The same equation, with $\alpha$ complex, describes AC voltage along a transmission line ([Ch. 13](13-partial-differential-equations.md)).

```desmos-graph
left=0; right=1000; bottom=0; top=1.1
---
y=\frac{\cosh(0.001(1000-x))}{\cosh(1)}|blue|0<x<1000
y=0.648|grey|dashed
```
*Normalised DC voltage $V/V_0$ along the 1 km leaky line of Example 11.44. The curve meets the open end with zero slope (zero current); the far end sits at $1/\cosh 1 = 0.648$ of the input.*

> [!example] Worked Example 11.45 — Start-up of a negative-resistance oscillator
> An LC tank in parallel with a nonlinear negative conductance (a tunnel diode, or a cross-coupled transistor pair) gives, in normalised units, the **van der Pol equation**
> $$
> \ddot v - \mu\,(1 - v^2)\,\dot v + v = 0,\qquad \mu = 0.5.
> $$
> Show that the oscillator starts by itself.
>
> **Step 1 (system).** $x_1 = v$, $x_2 = \dot v$: $\dot x_1 = x_2$, $\dot x_2 = -x_1 + \mu(1 - x_1^2)x_2$. The only equilibrium is the origin.
>
> **Step 2 (linearise).** $\mathbf{J} = \begin{pmatrix}0 & 1\\ -1 - 2\mu x_1x_2 & \mu(1 - x_1^2)\end{pmatrix}$, which at the origin is $\begin{pmatrix}0 & 1\\ -1 & \mu\end{pmatrix}$: $T = 0.5 > 0$, $\Delta = 1$, $T^2 < 4\Delta$.
>
> **Step 3 (classify).** Eigenvalues $\lambda = \tfrac12\big(\mu \pm \sqrt{\mu^2 - 4}\big) = 0.25 \pm 0.968j$: an **unstable spiral**. Any noise voltage grows as $e^{0.25t}$ while oscillating at $\approx 0.97$ (normalised) rad/s. The equilibrium is hyperbolic, so the nonlinear circuit does the same near the origin.
>
> **Step 4 (amplitude).** Growth cannot continue forever, because for $\lvert v\rvert > 1$ the term $\mu(1 - v^2)$ becomes *positive damping*. Averaging the energy balance over a cycle leads to an amplitude equation of exactly the Bernoulli form of Example 11.13, whose stable equilibrium gives the steady amplitude ($\approx 2$ for small $\mu$). Linear analysis predicts *whether* an oscillator starts; nonlinear analysis predicts *how large* it gets.

> [!warning] Common pitfall — dependent energy stores
> Two capacitors in parallel (or a loop consisting only of capacitors and voltage sources), or two inductors in series (or a cutset of only inductors and current sources), do not add an order: their states are linked by an algebraic constraint. Counting them as independent gives a singular $\mathbf{A}$ or an inconsistent set of initial conditions, and in a simulator can produce "impulsive" currents at switching instants.

> [!tip] Engineering note — what comes next
> This chapter *sets up* circuit equations and gives the general solution structure, plus the methods that survive nonlinear and time-varying elements. The systematic constant-coefficient theory — characteristic roots, damping ratio, natural and forced response, steady-state AC, state-variable form and $e^{\mathbf{A}t}$ — is [Ch. 12](12-linear-differential-equations.md). s-domain impedance methods that solve the same IVPs algebraically, impulses included, are [Ch. 14](../part-5-transforms/14-laplace-transforms.md). Nonlinear circuits beyond the reach of separation or linearisation are integrated numerically ([Ch. 22](../part-7-numerical-computational/22-numerical-methods.md)), with the existence and uniqueness theory of §11.7 as the guarantee that the numerical answer means something.

---

## Summary

| Topic | Key formula / fact |
|---|---|
| Autonomous stability | $y' = f(y)$, $f(y^*) = 0$: stable if $f'(y^*) < 0$, unstable if $f'(y^*) > 0$ |
| Exact equation | $M_y = N_t \Rightarrow \Phi(t, y) = C$ with $\Phi_t = M$, $\Phi_y = N$ |
| Integrating factor for exactness | $(M_y - N_t)/N = R(t) \Rightarrow \mu = e^{\int R\,dt}$ |
| Separable | $y' = g(t)h(y) \Rightarrow \int dy/h(y) = \int g\,dt + C$; check roots of $h$ |
| Scale-invariant substitution | $y' = F(y/t)$: $v = y/t$, $tv' = F(v) - v$ |
| Linear first order | $y' + py = r$: $y = y_0e^{-P(t)} + \int_{t_0}^{t}e^{-[P(t) - P(s)]}r(s)\,ds$, $P = \int_{t_0}^{t}p$ |
| Bernoulli | $y' + py = qy^n$: $v = y^{1-n}$, $v' + (1 - n)pv = (1 - n)q$ |
| Riccati | known $y_1$: $y = y_1 + 1/v$, $v' = -(q_1 + 2q_2y_1)v - q_2$ |
| Wronskian / Abel | $W = y_1y_2' - y_1'y_2$; $W(t) = W(t_0)e^{-\int_{t_0}^{t}p\,ds}$ (order $n$: use $p_{n-1}$) |
| Reduction of order | $y_2 = y_1\int e^{-\int p\,dt}/y_1^2\,dt$ |
| Cauchy–Euler | $t^2y'' + aty' + by = 0$: $y = t^m$, $m^2 + (a - 1)m + b = 0$; $t = e^x$ gives constant coefficients |
| Energy integral | $y'' = -U'(y)$: $\tfrac12(y')^2 + U(y) = E$ |
| Companion form | $x_k = y^{(k-1)}$; $\det(s\mathbf{I} - \mathbf{A}) = s^n + p_{n-1}s^{n-1} + \cdots + p_0$ |
| Structure theorem | $y = y_h + y_p$; $y_h$ has $n$ free constants |
| Variation of parameters | $y_p = -y_1\int y_2r/W\,dt + y_2\int y_1r/W\,dt$ ($r$ from standard form) |
| Causal Green's function | $G(t, s) = [y_1(s)y_2(t) - y_1(t)y_2(s)]/W(s)$; $y_p = \int_{t_0}^{t}G\,r\,ds$ |
| Picard–Lindelöf | $f$ continuous, Lipschitz in $y$: unique solution on $\lvert t - t_0\rvert \le \min(a, b/M)$ |
| Continuous dependence | $\lvert y - z\rvert \le \lvert y_0 - z_0\rvert e^{K\lvert t - t_0\rvert}$ |
| Linear IVP | unique solution on the whole interval where coefficients are continuous |
| BVP solvability | unique $\iff$ homogeneous BVP has only $y \equiv 0$ (Fredholm alternative) |
| Model eigenvalue problem | $y'' + \lambda y = 0$, $y(0) = y(L) = 0$: $\lambda_n = (n\pi/L)^2$, $y_n = \sin(n\pi x/L)$ |
| Sturm–Liouville | $(Py')' + Qy + \lambda wy = 0$: real simple eigenvalues, $\int y_my_nw\,dx = 0$ for $m \ne n$ |
| Systems | $\mathbf{x} = \boldsymbol{\Phi}\boldsymbol{\Phi}_0^{-1}\mathbf{x}_0 + \boldsymbol{\Phi}\int\boldsymbol{\Phi}^{-1}\mathbf{f}\,ds$; $\det\boldsymbol{\Phi} \propto e^{\int\operatorname{tr}\mathbf{A}}$ |
| Eigenvector method | $\mathbf{x} = \sum c_k\mathbf{v}_ke^{\lambda_kt}$; defective: $(t\mathbf{v} + \mathbf{w})e^{\lambda t}$, $(\mathbf{A} - \lambda\mathbf{I})\mathbf{w} = \mathbf{v}$ |
| $2 \times 2$ stability | asymptotically stable $\iff \operatorname{tr}\mathbf{A} < 0$ and $\det\mathbf{A} > 0$ |
| Linearisation | $\dot{\boldsymbol{\eta}} = \mathbf{J}\boldsymbol{\eta}$, valid locally when no eigenvalue of $\mathbf{J}$ has zero real part |
| Circuit element laws | $v_L = L\,di/dt$, $i_C = C\,dv/dt$; states $v_C$, $i_L$ are continuous; order = number of independent stores |

---

## Further reading

- E. Kreyszig, *Advanced Engineering Mathematics*, 10th ed., Wiley — Chapters 1–5 (first-order, second-order and higher-order ODEs, systems, series solutions) and Chapter 11 (Sturm–Liouville problems).
- W. E. Boyce, R. C. DiPrima and D. B. Meade, *Elementary Differential Equations and Boundary Value Problems*, 11th ed., Wiley — Chapters 2, 3, 7, 9 and 11; the standard reference for existence and uniqueness, phase-plane analysis and boundary-value problems.
- G. F. Simmons, *Differential Equations with Applications and Historical Notes*, 3rd ed., CRC Press — clear proofs of Picard's theorem, Sturm comparison and oscillation theorems.
- S. H. Strogatz, *Nonlinear Dynamics and Chaos*, 2nd ed., CRC Press — the geometric (phase-line and phase-plane) view of nonlinear systems, including oscillators and bifurcations.
