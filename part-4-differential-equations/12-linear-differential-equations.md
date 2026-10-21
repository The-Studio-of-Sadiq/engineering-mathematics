# Chapter 12 — Linear Differential Equations

[← Previous](11-ordinary-differential-equations.md) · [Table of Contents](../00-front/table-of-contents.md) · [Next →](13-partial-differential-equations.md)

> **Prerequisites:** Ch. 2 (Complex Numbers), Ch. 4 (Differential Calculus), Ch. 8 and Ch. 10 (Matrices, Eigenvalues), Ch. 11 (Ordinary Differential Equations)
>
> **Learning outcomes:** After this chapter you can
> - write a linear constant-coefficient ODE in operator form $P(D)y = r(t)$ and use linearity, superposition and time invariance to simplify it;
> - find the homogeneous solution from the roots of the characteristic polynomial in all cases (distinct, repeated, complex), and classify a second-order system as overdamped, critically damped, underdamped or undamped;
> - split a response into natural and forced parts, or into zero-input and zero-state parts, and fit initial conditions correctly;
> - find particular solutions by undetermined coefficients, the complex-exponential (phasor) method and the resonance rule;
> - separate transient from steady-state behaviour, compute frequency response, settling time and overshoot;
> - solve coupled linear equations by elimination and by modal decomposition, and write any linear lumped system as $\dot{\mathbf{x}} = \mathbf{A}\mathbf{x} + \mathbf{B}u$ with solution through the matrix exponential $e^{\mathbf{A}t}$.

[Chapter 11](11-ordinary-differential-equations.md) laid out the general theory of ordinary differential equations: the structure theorem $y = y_h + y_p$, the Wronskian, variation of parameters, existence and uniqueness, and the eigenvector method for first-order systems. This chapter specialises to the case that dominates electrical engineering: **linear equations with constant coefficients**. Every circuit built from ideal resistors, capacitors, inductors, coupled inductors, controlled sources and op-amps operating linearly is described by such equations; so are small-signal models of amplifiers, linearised motors, filters, control loops and transmission-line sections.

Constant coefficients change the mathematics fundamentally. Instead of searching for unknown functions, we search for **roots of a polynomial**. Calculus problems become algebra problems, and the algebra has a clear physical reading: the roots are the natural frequencies of the system, their real parts say how fast energy dissipates, and their imaginary parts say how fast it sloshes between storage elements. The same polynomial reappears as the denominator of the transfer function in [Chapter 14](../part-5-transforms/14-laplace-transforms.md) and as the characteristic polynomial of the state matrix in [Section 12.7](#127-state-variable-formulation).

The chapter proceeds in eight sections. Sections 12.1–12.2 build the homogeneous theory. Sections 12.3–12.5 add forcing and interpret the result physically. Sections 12.6–12.7 move to several interacting unknowns and to the state-space form that modern analysis and simulation software uses. Section 12.8 pulls everything together on RLC circuits and an electromechanical system.

---

## 12.1 Constant-coefficient equations

### Motivation

Write Kirchhoff's voltage law around a series loop containing a resistor $R$, an inductor $L$, a capacitor $C$ and a source $v_s(t)$. Using the element laws $v_R = Ri$, $v_L = L\,di/dt$, $i = C\,dv_C/dt$:

$$
LC\,\frac{d^2v_C}{dt^2} + RC\,\frac{dv_C}{dt} + v_C = v_s(t).
$$

The coefficients $LC$, $RC$ and $1$ are fixed by the component values; they do not depend on time or on $v_C$. That is the defining feature of the class studied here. The same form arises for a mass–spring–damper ($m\ddot x + b\dot x + kx = F$), a torsional shaft, a thermal network or the linearised dynamics of a phase-locked loop. Because the form is universal, a single toolkit solves all of them.

### Definitions

> [!definition] Linear constant-coefficient ODE
> A **linear ordinary differential equation of order $n$ with constant coefficients** is
> $$
> a_n y^{(n)} + a_{n-1} y^{(n-1)} + \dots + a_1 y' + a_0 y = r(t), \qquad a_n \ne 0,
> $$
> where $a_0, \dots, a_n$ are real (occasionally complex) constants and $r(t)$ is a given **forcing function**, **input** or **excitation**. The equation is **homogeneous** if $r(t) \equiv 0$ and **nonhomogeneous** (or *forced*) otherwise. Dividing by $a_n$ gives the **monic** form with $a_n = 1$.

Three words carry the content of the definition:

- *Linear*: $y$ and its derivatives appear only to the first power and are never multiplied together or placed inside a nonlinear function.
- *Constant coefficient*: the $a_k$ do not vary with $t$. Physically, the components do not age, heat up or switch during the interval considered.
- *Order $n$*: the highest derivative present. For lumped circuits $n$ equals the number of independent energy-storage elements (capacitors and inductors that cannot be combined into one equivalent element).

### The differential operator $D$

It is efficient to treat differentiation as an algebraic object. Let $D = d/dt$, so $Dy = y'$, $D^2 y = y''$ and in general $D^k y = y^{(k)}$, with $D^0 = 1$ (the identity).

> [!definition] Operator polynomial
> For the polynomial $P(s) = a_n s^n + \dots + a_1 s + a_0$, the **operator polynomial** $P(D)$ acts on a sufficiently differentiable function $y$ by
> $$
> P(D)\,y = a_n D^n y + \dots + a_1 D y + a_0 y .
> $$
> The equation of the definition above is then simply $P(D)\,y = r(t)$. The polynomial $P(s)$ is called the **characteristic polynomial** of the equation.

The key fact that makes $D$ behave like a number is that constant-coefficient operators **commute** and **factor** exactly like the polynomials that define them.

> [!theorem] Algebra of constant-coefficient operators
> If $P(s)$ and $Q(s)$ are polynomials with constant coefficients, then for any sufficiently differentiable $y$
> $$
> P(D)\big[Q(D)\,y\big] = Q(D)\big[P(D)\,y\big] = (PQ)(D)\,y .
> $$
> In particular, if $P(s) = a_n (s - s_1)(s - s_2)\cdots(s - s_n)$ then
> $P(D) = a_n (D - s_1)(D - s_2)\cdots(D - s_n)$, with the factors in any order.

*Proof sketch.* By linearity it suffices to check monomials: $D^j(D^k y) = D^{j+k} y = D^k(D^j y)$. Constants pass through $D$ because $D(cy) = cDy$ when $c$ is constant. Expanding both sides term by term gives the same coefficients as polynomial multiplication. $\square$

> [!warning] Common Pitfall
> The commuting property fails as soon as a coefficient depends on $t$. For example $D(t\,y) = y + t\,Dy \ne t\,Dy$. So $(D + t)(D - t) \ne (D - t)(D + t)$, and the factor-and-solve techniques of this chapter do **not** apply to variable-coefficient equations such as Cauchy–Euler or Bessel equations (see [Ch. 11](11-ordinary-differential-equations.md)).

### Linearity and superposition

> [!theorem] Superposition
> $P(D)$ is a linear operator: for constants $c_1, c_2$ and functions $y_1, y_2$,
> $$
> P(D)\,(c_1 y_1 + c_2 y_2) = c_1 P(D)y_1 + c_2 P(D)y_2 .
> $$
> Consequently:
> 1. any linear combination of solutions of the homogeneous equation $P(D)y = 0$ is again a solution;
> 2. if $P(D)y_1 = r_1$ and $P(D)y_2 = r_2$, then $P(D)(c_1 y_1 + c_2 y_2) = c_1 r_1 + c_2 r_2$;
> 3. the general solution of $P(D)y = r$ is $y = y_h + y_p$, where $y_h$ is the general homogeneous solution ($n$ arbitrary constants) and $y_p$ is any one particular solution.

Item 3 is the structure theorem of [Ch. 11](11-ordinary-differential-equations.md); it is restated here because everything that follows is a recipe for finding $y_h$ (Section 12.2) and $y_p$ (Section 12.4) when coefficients are constant.

### Time invariance

Constant coefficients bring a second structural property that variable-coefficient equations lack.

> [!theorem] Time invariance (shift property)
> If $y(t)$ solves $P(D)y = r(t)$, then for any fixed $t_0$ the shifted function $y(t - t_0)$ solves $P(D)y = r(t - t_0)$.

*Proof.* Let $w(t) = y(t - t_0)$. By the chain rule $w^{(k)}(t) = y^{(k)}(t - t_0)$, so $P(D)w(t) = \big[P(D)y\big](t - t_0) = r(t - t_0)$. The coefficients do not change under the shift because they are constants. $\square$

Linearity plus time invariance is the defining pair of properties of an **LTI system** (linear time-invariant). Everything in signal processing that relies on convolution, transfer functions and frequency response depends on both. A circuit with a time-varying element (a switched capacitor, a mixer, a varactor driven by a pump) is linear but not time-invariant, and its equations do not have constant coefficients.

### The exponential as eigenfunction

Why should exponentials be the building blocks? Apply $D$ to $e^{st}$, where $s$ is any complex constant:

$$
D\,e^{st} = s\,e^{st}, \qquad D^k e^{st} = s^k e^{st}, \qquad P(D)\,e^{st} = P(s)\,e^{st}.
$$

So $e^{st}$ is an **eigenfunction** of every constant-coefficient operator, with eigenvalue $P(s)$. An LTI system therefore maps an exponential input into the same exponential, merely scaled by a complex number. This single identity drives the characteristic equation (Section 12.2), the phasor method (Section 12.4), frequency response (Section 12.5) and, after integration against $e^{-st}$, the Laplace transform ([Ch. 14](../part-5-transforms/14-laplace-transforms.md)).

A useful generalisation handles exponentials multiplied by other functions.

> [!theorem] Exponential shift theorem
> For any constant $a$ and sufficiently differentiable $f$,
> $$
> P(D)\big[e^{at} f(t)\big] = e^{at}\,P(D + a)\,f(t).
> $$

*Proof.* For one factor, $D\big[e^{at}f\big] = a e^{at} f + e^{at} f' = e^{at}(D + a) f$. Applying this $k$ times gives $D^k\big[e^{at} f\big] = e^{at}(D + a)^k f$, and summing with coefficients $a_k$ gives the result. $\square$

The shift theorem lets us "pull an exponential through" an operator. It explains at once why repeated roots produce factors of $t$ (Section 12.2) and why resonant forcing produces $t\,e^{at}$ terms (Section 12.4).

### Real solutions from complex ones

When the coefficients $a_k$ are real and $y = u + jv$ is a complex-valued solution of $P(D)y = r$ with $r$ real, then taking real and imaginary parts gives $P(D)u = r$ and $P(D)v = 0$. More generally:

> [!theorem] Real and imaginary parts
> If the coefficients of $P$ are real and $P(D)z = w(t)$ with $z$, $w$ complex, then
> $$
> P(D)\,\operatorname{Re}\{z\} = \operatorname{Re}\{w\}, \qquad P(D)\,\operatorname{Im}\{z\} = \operatorname{Im}\{w\}.
> $$

This is the justification for solving with complex exponentials $e^{j\omega t}$ and then taking real parts, which is precisely what phasor analysis does.

### Worked examples

> [!example] Example 12.1 — First-order RC charging
> A capacitor $C$ charges through a resistor $R$ from a DC source $V_0$ switched on at $t = 0$; the capacitor is initially uncharged. Find $v_C(t)$.
>
> **Model.** KVL gives $Ri + v_C = V_0$ and $i = C\,dv_C/dt$, so
> $$
> RC\,\frac{dv_C}{dt} + v_C = V_0, \qquad v_C(0) = 0 .
> $$
> With $\tau = RC$ (units: $\Omega\cdot\mathrm{F} = \mathrm{s}$) the operator form is $(\tau D + 1)v_C = V_0$, with characteristic polynomial $P(s) = \tau s + 1$ and root $s_1 = -1/\tau$.
>
> **Homogeneous part.** $v_h = A e^{-t/\tau}$.
>
> **Particular part.** The input is constant, so try a constant $v_p = K$: then $0 + K = V_0$, $K = V_0$.
>
> **Initial condition.** $v_C = V_0 + A e^{-t/\tau}$; $v_C(0) = V_0 + A = 0$ gives $A = -V_0$:
> $$
> v_C(t) = V_0\left(1 - e^{-t/\tau}\right), \qquad t \ge 0 .
> $$
>
> **Check.** At $t = \tau$, $v_C = V_0(1 - e^{-1}) = 0.632\,V_0$; at $t = 5\tau$, $v_C = 0.9933\,V_0$. Substituting: $\tau \cdot (V_0/\tau)e^{-t/\tau} + V_0 - V_0e^{-t/\tau} = V_0$. ✓ With $R = 10\ \mathrm{k\Omega}$, $C = 100\ \mathrm{nF}$, $\tau = 1\ \mathrm{ms}$, the capacitor is within 1% of $V_0$ after $4.6\ \mathrm{ms}$ (since $e^{-4.6} = 0.010$).

```desmos-graph
left=-0.3; right=5; bottom=-0.1; top=1.15
---
y=1-e^{-x}|blue|0<x<5
y=1|grey|dashed
y=x|red|dashed|0<x<1.1
(1,0.632)|label:t = τ (63.2 percent)
x=1|grey|dotted
```
*Normalised RC charging curve $v_C/V_0$ against $t/\tau$; the initial tangent (red) reaches the final value exactly at $t = \tau$.*

> [!example] Example 12.2 — Deriving and factoring a second-order operator
> Show that $y'' + y' - 2y = 0$ factors as $(D - 1)(D + 2)y = 0$ and verify directly that $e^{t}$ and $e^{-2t}$ are solutions.
>
> **Factoring.** $P(s) = s^2 + s - 2 = (s - 1)(s + 2)$, so by the operator algebra theorem $P(D) = (D - 1)(D + 2) = (D + 2)(D - 1)$.
>
> **Check by expansion.** $(D - 1)(D + 2)y = (D - 1)(y' + 2y) = y'' + 2y' - y' - 2y = y'' + y' - 2y$. ✓
>
> **Solutions.** $(D - 1)e^{t} = e^{t} - e^{t} = 0$, so $(D + 2)(D - 1)e^{t} = (D+2)\,0 = 0$. Likewise $(D + 2)e^{-2t} = 0$, so $(D-1)(D+2)e^{-2t} = 0$. Alternatively, using the eigenfunction identity, $P(D)e^{st} = P(s)e^{st}$ vanishes exactly when $P(s) = 0$, i.e. $s = 1$ or $s = -2$.
>
> The general solution is $y = c_1 e^{t} + c_2 e^{-2t}$; the $e^{t}$ term grows, so this system is unstable.

> [!example] Example 12.3 — Using the exponential shift theorem
> Evaluate $P(D)\big[t e^{-3t}\big]$ for $P(D) = D^2 + 6D + 9$, and interpret.
>
> By the shift theorem with $a = -3$, $f = t$:
> $$
> P(D)\big[e^{-3t} t\big] = e^{-3t} P(D - 3)\,t = e^{-3t}\big[(D-3)^2 + 6(D-3) + 9\big]t = e^{-3t} D^2 t = 0 .
> $$
> The operator $(D - 3)^2 + 6(D - 3) + 9$ simplifies to $D^2$ because $P(s) = (s + 3)^2$, so $P(s - 3) = s^2$. Hence $t e^{-3t}$ is a solution of $y'' + 6y' + 9y = 0$, alongside $e^{-3t}$. This is the repeated-root case of Section 12.2, obtained here without any guesswork.
>
> **Direct check.** $y = te^{-3t}$, $y' = (1 - 3t)e^{-3t}$, $y'' = (-6 + 9t)e^{-3t}$. Then $y'' + 6y' + 9y = (-6 + 9t + 6 - 18t + 9t)e^{-3t} = 0$. ✓

> [!example] Example 12.4 — Series RLC in operator form
> For the series RLC loop of the motivation with $R = 100\ \Omega$, $L = 10\ \mathrm{mH}$, $C = 1\ \mu\mathrm{F}$, write the equation for $v_C$ in monic operator form and find its characteristic polynomial.
>
> Divide $LC\,v_C'' + RC\,v_C' + v_C = v_s$ by $LC = 10^{-8}\ \mathrm{s^2}$:
> $$
> \left(D^2 + \frac{R}{L}D + \frac{1}{LC}\right)v_C = \frac{1}{LC}v_s, \qquad P(s) = s^2 + 10^4 s + 10^8 .
> $$
> Dimensional check: $R/L$ has units $\Omega/\mathrm{H} = \mathrm{s^{-1}}$ and $1/LC$ has units $\mathrm{s^{-2}}$, matching $D$ and $D^2$. The roots, $s = -5000 \pm j\,8660\ \mathrm{s^{-1}}$, are found in Section 12.2 and the full step response in Section 12.8.

> [!tip] Engineering Note
> Circuit simulators (SPICE) never form the single $n$th-order equation; they assemble first-order equations for every node and storage element (modified nodal analysis), which is the state-space viewpoint of Section 12.7. Hand analysis, by contrast, almost always reduces to one equation $P(D)y = r$, because the polynomial $P(s)$ is where the design insight lives: its roots set rise time, ringing and stability.

---

## 12.2 Characteristic equations

### Motivation

Section 12.1 showed that $e^{st}$ solves $P(D)y = 0$ whenever $P(s) = 0$. Solving a homogeneous constant-coefficient ODE therefore reduces to finding the $n$ roots of a degree-$n$ polynomial and handling two complications: roots may be **repeated**, and roots of a real polynomial may be **complex**. This section treats both and then studies the second-order case in detail, because nearly every practical design specification (damping, overshoot, ringing frequency, quality factor) is phrased in terms of a dominant pair of roots.

### The characteristic equation

> [!definition] Characteristic equation and natural frequencies
> For $P(D)y = 0$ with $P(s) = a_n s^n + \dots + a_0$, the algebraic equation
> $$
> P(s) = a_n s^n + a_{n-1} s^{n-1} + \dots + a_1 s + a_0 = 0
> $$
> is the **characteristic (auxiliary) equation**. Its roots $s_1, \dots, s_n$ (counted with multiplicity) are the **characteristic roots**, also called the **natural frequencies**, **modes** or **poles** of the system. They have units of $\mathrm{s^{-1}}$ when $t$ is in seconds.

One obtains the characteristic equation mechanically by replacing $y^{(k)}$ with $s^k$. By the fundamental theorem of algebra, a degree-$n$ polynomial has exactly $n$ complex roots counted with multiplicity, and when the coefficients are real, non-real roots occur in **conjugate pairs** $\sigma \pm j\omega$.

### Case 1: distinct real roots

> [!theorem] Distinct roots
> If $P(s)$ has $n$ distinct roots $s_1, \dots, s_n$, then $e^{s_1 t}, \dots, e^{s_n t}$ form a fundamental set of solutions of $P(D)y = 0$, and the general solution is
> $$
> y_h(t) = c_1 e^{s_1 t} + c_2 e^{s_2 t} + \dots + c_n e^{s_n t}.
> $$

*Proof sketch.* Each is a solution by the eigenfunction identity. Independence follows from the Wronskian at $t = 0$, which is the Vandermonde determinant $\prod_{i<k}(s_k - s_i) \ne 0$ for distinct roots ([Ch. 8](../part-3-linear-algebra/08-matrices-and-linear-systems.md)). By the theory of [Ch. 11](11-ordinary-differential-equations.md), $n$ independent solutions of an $n$th-order linear equation span all solutions. $\square$

The theorem holds for distinct complex roots too; Case 3 converts the complex exponentials into real functions.

### Case 2: repeated roots

If $s_1$ is a root of multiplicity $m$, the list of exponentials contains $e^{s_1 t}$ only once, so we are $m - 1$ solutions short. The exponential shift theorem supplies them.

> [!theorem] Repeated roots
> If $s_1$ is a root of $P(s)$ of multiplicity $m$, then
> $$
> e^{s_1 t},\ t e^{s_1 t},\ t^2 e^{s_1 t},\ \dots,\ t^{m-1} e^{s_1 t}
> $$
> are $m$ linearly independent solutions of $P(D)y = 0$. Collecting such blocks over all distinct roots gives a fundamental set of $n$ solutions.

*Derivation.* Write $P(s) = (s - s_1)^m Q(s)$. Then $P(D) = Q(D)(D - s_1)^m$, and by the shift theorem
$$
(D - s_1)^m\big[e^{s_1 t} t^k\big] = e^{s_1 t}\,\big(D + s_1 - s_1\big)^m t^k = e^{s_1 t} D^m t^k = 0 \quad\text{for } k < m,
$$
because the $m$th derivative of a polynomial of degree $k < m$ vanishes. Applying $Q(D)$ to zero still gives zero. Independence: a relation $\sum_k c_k t^k e^{s_1 t} = 0$ for all $t$ forces the polynomial $\sum c_k t^k$ to vanish identically, so all $c_k = 0$; mixing different roots is handled by the generalised Vandermonde argument. $\square$

There is a second, intuitive derivation. For a double root, consider two nearby distinct roots $s_1$ and $s_1 + \varepsilon$. The combination $\big(e^{(s_1+\varepsilon)t} - e^{s_1 t}\big)/\varepsilon$ is a solution for every $\varepsilon \ne 0$, and as $\varepsilon \to 0$ it tends to $\partial e^{st}/\partial s\big|_{s_1} = t e^{s_1 t}$. Repeated roots are the limit of merging modes, and the factor $t$ records the merger.

### Case 3: complex-conjugate roots

For real coefficients, a complex root $s = \sigma + j\omega$ ($\omega > 0$) comes paired with $\bar s = \sigma - j\omega$. The complex solutions $e^{(\sigma \pm j\omega)t}$ are legitimate but inconvenient. Using Euler's formula ([Ch. 2](../part-1-foundations/02-complex-numbers.md)),
$$
e^{(\sigma \pm j\omega)t} = e^{\sigma t}(\cos\omega t \pm j\sin\omega t),
$$
and by the real/imaginary-parts theorem of Section 12.1, the real and imaginary parts are themselves real solutions.

> [!theorem] Complex roots
> If $\sigma \pm j\omega$ is a pair of simple complex-conjugate roots of a real polynomial $P(s)$, the corresponding contribution to $y_h$ is
> $$
> e^{\sigma t}\big(A\cos\omega t + B\sin\omega t\big) = M e^{\sigma t}\cos(\omega t - \varphi),
> $$
> with real constants $A, B$, or equivalently amplitude $M = \sqrt{A^2 + B^2}$ and phase $\varphi = \operatorname{atan2}(B, A)$. If the pair has multiplicity $m$, multiply by $1, t, \dots, t^{m-1}$ as in Case 2.

The real part $\sigma$ sets the **envelope** $e^{\sigma t}$ (decay if $\sigma < 0$, growth if $\sigma > 0$), and the imaginary part $\omega$ sets the **oscillation frequency** in rad/s. A root located in the complex $s$-plane can be read at a glance: distance to the left of the imaginary axis is decay rate, height above the real axis is ringing frequency.

> [!warning] Common Pitfall
> Writing the complex-root solution as $c_1 e^{\sigma t}\cos\omega t + c_2 e^{\sigma t}\cos\omega t$ (the same function twice) or forgetting the $e^{\sigma t}$ envelope are the two most common errors. Each conjugate *pair* contributes exactly two real functions, $e^{\sigma t}\cos\omega t$ and $e^{\sigma t}\sin\omega t$. Also note: one pair accounts for two roots, so a cubic with one complex pair has one more real root to find.

### Summary of the root-to-solution map

| Root(s) of $P(s)$ | Contribution to $y_h$ |
|---|---|
| real $s_1$, simple | $c_1 e^{s_1 t}$ |
| real $s_1$, multiplicity $m$ | $(c_1 + c_2 t + \dots + c_m t^{m-1})e^{s_1 t}$ |
| $\sigma \pm j\omega$, simple | $e^{\sigma t}(A\cos\omega t + B\sin\omega t)$ |
| $\sigma \pm j\omega$, multiplicity $m$ | $e^{\sigma t}\big[(A_1 + \dots + A_m t^{m-1})\cos\omega t + (B_1 + \dots + B_m t^{m-1})\sin\omega t\big]$ |
| $s = 0$, multiplicity $m$ | $c_1 + c_2 t + \dots + c_m t^{m-1}$ |
| $\pm j\omega$ (purely imaginary) | $A\cos\omega t + B\sin\omega t$ (sustained oscillation) |

### Stability from the roots

> [!definition] Stability of a homogeneous LTI equation
> The equation $P(D)y = 0$ is
> - **asymptotically stable** if every solution tends to $0$ as $t \to \infty$;
> - **marginally stable** if every solution stays bounded but some do not decay;
> - **unstable** if some solution is unbounded.

> [!theorem] Root criterion
> $P(D)y = 0$ is asymptotically stable if and only if every characteristic root satisfies $\operatorname{Re}\{s_k\} < 0$ (all roots in the open left half-plane). It is marginally stable if all roots satisfy $\operatorname{Re}\{s_k\} \le 0$ and those on the imaginary axis are simple. Otherwise it is unstable.

*Reason.* $\lvert t^k e^{s t}\rvert = t^k e^{\sigma t}$, which tends to zero for any $k$ when $\sigma < 0$ (exponential decay beats polynomial growth), stays bounded when $\sigma = 0$ and $k = 0$, and grows otherwise. $\square$

For a real polynomial of degree two, stability is equivalent to all coefficients having the same sign. For degree three, $s^3 + a_2 s^2 + a_1 s + a_0$ is stable if and only if $a_2, a_1, a_0 > 0$ and $a_2 a_1 > a_0$. For higher degrees, the **Routh–Hurwitz criterion** tests stability from the coefficients without computing roots; positive coefficients are necessary but not sufficient beyond degree two. Detailed Routh tables belong to control theory; the point here is that stability is a property of the characteristic polynomial alone.

### The standard second-order form

Every stable second-order equation with real coefficients can be written

$$
y'' + 2\zeta\omega_n\,y' + \omega_n^2\,y = \omega_n^2\,K\,u(t),
$$

where

- $\omega_n > 0$ is the **undamped natural frequency** (rad/s),
- $\zeta \ge 0$ is the **damping ratio** (dimensionless),
- $K$ is the DC gain (irrelevant for the homogeneous problem).

The characteristic equation $s^2 + 2\zeta\omega_n s + \omega_n^2 = 0$ has roots

$$
s_{1,2} = -\zeta\omega_n \pm \omega_n\sqrt{\zeta^2 - 1}.
$$

The product of the roots is $\omega_n^2$, so both roots lie on (or, when real, on either side of) a circle of radius $\omega_n$ in the $s$-plane; the sum is $-2\zeta\omega_n$. The damping ratio alone decides the character of the response.

> [!definition] Damping cases
> | Case | Condition | Roots | Homogeneous solution |
> |---|---|---|---|
> | Overdamped | $\zeta > 1$ | two distinct negative reals $-\alpha_1, -\alpha_2$ | $c_1 e^{-\alpha_1 t} + c_2 e^{-\alpha_2 t}$ |
> | Critically damped | $\zeta = 1$ | double root $-\omega_n$ | $(c_1 + c_2 t)e^{-\omega_n t}$ |
> | Underdamped | $0 < \zeta < 1$ | $-\alpha \pm j\omega_d$ | $e^{-\alpha t}(A\cos\omega_d t + B\sin\omega_d t)$ |
> | Undamped | $\zeta = 0$ | $\pm j\omega_n$ | $A\cos\omega_n t + B\sin\omega_n t$ |
>
> In the underdamped case, $\alpha = \zeta\omega_n$ is the **decay rate** (neper frequency) and $\omega_d = \omega_n\sqrt{1 - \zeta^2}$ is the **damped natural frequency**. Note $\alpha^2 + \omega_d^2 = \omega_n^2$.

The geometric picture is useful. As $\zeta$ increases from $0$, the root pair starts on the imaginary axis at $\pm j\omega_n$, travels along the circle of radius $\omega_n$ into the left half-plane, meets on the real axis at $-\omega_n$ when $\zeta = 1$, and then splits along the real axis, one root heading towards $0$ and the other towards $-\infty$. The angle $\theta$ between a complex root and the negative real axis satisfies $\cos\theta = \zeta$.

Consider the response from $y(0) = 1$, $y'(0) = 0$ with $\omega_n = 1$:

- $\zeta = 2$: roots $-2 \pm \sqrt{3} = -0.2679,\ -3.7321$; $y = 1.0774\,e^{-0.2679t} - 0.0774\,e^{-3.7321t}$.
- $\zeta = 1$: $y = (1 + t)e^{-t}$.
- $\zeta = 0.2$: $\omega_d = 0.9798$; $y = e^{-0.2t}\big(\cos 0.9798t + 0.2041\sin 0.9798t\big)$.
- $\zeta = 0$: $y = \cos t$.

```desmos-graph
left=-0.5; right=15; bottom=-1.1; top=1.15
---
y=1.0774e^{-0.2679x}-0.0774e^{-3.7321x}|red|0<x<15
y=(1+x)e^{-x}|green|0<x<15
y=e^{-0.2x}(\cos(0.9798x)+0.2041\sin(0.9798x))|blue|0<x<15
y=\cos(x)|grey|dashed|0<x<15
y=e^{-0.2x}|blue|dotted|0<x<15
y=-e^{-0.2x}|blue|dotted|0<x<15
```
*Free response from $y(0)=1$, $y'(0)=0$ with $\omega_n = 1$: overdamped $\zeta = 2$ (red), critically damped $\zeta = 1$ (green), underdamped $\zeta = 0.2$ (blue, with its $e^{-0.2t}$ envelope dotted) and undamped $\zeta = 0$ (grey). Critical damping returns to zero fastest without crossing it.*

Three features of the graph are worth stating as facts:

1. **Critical damping is the fastest non-oscillatory return.** The overdamped response is dominated by its slow root $-0.268$; increasing damping beyond $\zeta = 1$ actually slows the system down, because the slow root moves towards the origin ($s_1 \approx -\omega_n/(2\zeta)$ for large $\zeta$).
2. **An underdamped response crosses zero infinitely often**, at intervals of $\pi/\omega_d$, while its envelope decays with time constant $1/\alpha$.
3. **An overdamped or critically damped response crosses zero at most once.** Since $c_1 e^{-\alpha_1 t} + c_2 e^{-\alpha_2 t} = 0$ means $e^{(\alpha_2 - \alpha_1)t} = -c_2/c_1$, there is at most one solution $t$.

### Fitting initial conditions

The $n$ constants in $y_h$ are fixed by $n$ initial conditions $y(0), y'(0), \dots, y^{(n-1)}(0)$. For distinct roots this gives a Vandermonde system:

$$
\begin{bmatrix} 1 & 1 & \cdots & 1 \\ s_1 & s_2 & \cdots & s_n \\ \vdots & & & \vdots \\ s_1^{n-1} & s_2^{n-1} & \cdots & s_n^{n-1} \end{bmatrix}
\begin{bmatrix} c_1 \\ c_2 \\ \vdots \\ c_n \end{bmatrix}
=
\begin{bmatrix} y(0) \\ y'(0) \\ \vdots \\ y^{(n-1)}(0) \end{bmatrix}.
$$

For second order with distinct roots and data $y_0 = y(0)$, $v_0 = y'(0)$ the solution is

$$
c_1 = \frac{v_0 - s_2 y_0}{s_1 - s_2}, \qquad c_2 = \frac{s_1 y_0 - v_0}{s_1 - s_2}.
$$

For the underdamped case, $y = e^{-\alpha t}(A\cos\omega_d t + B\sin\omega_d t)$ gives $A = y_0$ and $-\alpha A + \omega_d B = v_0$, so

$$
A = y_0, \qquad B = \frac{v_0 + \alpha y_0}{\omega_d}.
$$

For the critically damped case $y = (c_1 + c_2 t)e^{-\omega_n t}$: $c_1 = y_0$, $c_2 = v_0 + \omega_n y_0$.

### Worked examples

> [!example] Example 12.5 — Overdamped: distinct real roots
> Solve $y'' + 5y' + 6y = 0$, $y(0) = 1$, $y'(0) = 0$.
>
> **Characteristic equation.** $s^2 + 5s + 6 = (s + 2)(s + 3) = 0$, roots $s_1 = -2$, $s_2 = -3$. (Here $\omega_n = \sqrt 6 = 2.449$, $\zeta = 5/(2\sqrt6) = 1.021$: just overdamped.)
>
> **Constants.** $c_1 = (v_0 - s_2 y_0)/(s_1 - s_2) = (0 + 3)/1 = 3$, $c_2 = (s_1 y_0 - v_0)/(s_1 - s_2) = -2$.
> $$
> y(t) = 3e^{-2t} - 2e^{-3t}.
> $$
> **Check.** $y(0) = 1$ ✓, $y'(0) = -6 + 6 = 0$ ✓. Substituting: $y'' = 12e^{-2t} - 18e^{-3t}$, $5y' = -30e^{-2t} + 30e^{-3t}$, $6y = 18e^{-2t} - 12e^{-3t}$; sum $= 0$ ✓. At $t = 1$: $y = 3(0.13534) - 2(0.04979) = 0.3064$.

> [!example] Example 12.6 — Critically damped: repeated root
> Solve $y'' + 4y' + 4y = 0$, $y(0) = 1$, $y'(0) = 0$.
>
> $s^2 + 4s + 4 = (s + 2)^2$: double root $s = -2$, so $y = (c_1 + c_2 t)e^{-2t}$. Then $c_1 = 1$ and $c_2 = v_0 + \omega_n y_0 = 0 + 2 = 2$:
> $$
> y(t) = (1 + 2t)e^{-2t}.
> $$
> **Check.** $y' = 2e^{-2t} - 2(1 + 2t)e^{-2t} = -4te^{-2t}$, so $y'(0) = 0$ ✓; $y'' = (-4 + 8t)e^{-2t}$; $y'' + 4y' + 4y = (-4 + 8t - 16t + 4 + 8t)e^{-2t} = 0$ ✓. Notice $y' < 0$ for all $t > 0$: the response decays monotonically without overshoot.

> [!example] Example 12.7 — Underdamped: complex roots
> Solve $y'' + 2y' + 10y = 0$, $y(0) = 1$, $y'(0) = 0$, and express the result in amplitude–phase form.
>
> $s^2 + 2s + 10 = 0$ gives $s = -1 \pm \sqrt{1 - 10} = -1 \pm j3$. So $\alpha = 1$, $\omega_d = 3$, $\omega_n = \sqrt{10} = 3.162$, $\zeta = 1/\sqrt{10} = 0.316$.
>
> $A = y_0 = 1$, $B = (v_0 + \alpha y_0)/\omega_d = 1/3$:
> $$
> y(t) = e^{-t}\left(\cos 3t + \tfrac13\sin 3t\right) = \frac{\sqrt{10}}{3}\,e^{-t}\cos(3t - 0.3218).
> $$
> Here $M = \sqrt{1 + 1/9} = 1.0541$ and $\varphi = \arctan(1/3) = 0.3218$ rad $= 18.43^\circ$.
>
> **Check.** $y'(0) = -1\cdot 1 + 3\cdot\tfrac13 = 0$ ✓. The oscillation period is $2\pi/3 = 2.094$ s and the envelope falls by a factor $e^{-2.094} = 0.123$ per period; this ratio, $e^{-2\pi\zeta/\sqrt{1-\zeta^2}}$, is the basis of the **logarithmic decrement** method of measuring damping from a recorded ring-down.

> [!example] Example 12.8 — Higher order with mixed roots
> Find the general solution of $y^{(4)} + 2y'' + y = 0$ and of $y''' - y'' + y' - y = 0$.
>
> **First equation.** $s^4 + 2s^2 + 1 = (s^2 + 1)^2$, so $\pm j$ are double roots. The table gives
> $$
> y = (A_1 + A_2 t)\cos t + (B_1 + B_2 t)\sin t .
> $$
> The $t\cos t$ and $t\sin t$ terms grow without bound, so the system is unstable even though all roots lie on the imaginary axis: repeated imaginary-axis roots violate the marginal-stability condition.
>
> **Second equation.** Factor by grouping: $s^3 - s^2 + s - 1 = s^2(s - 1) + (s - 1) = (s - 1)(s^2 + 1)$. Roots $1, \pm j$:
> $$
> y = c_1 e^{t} + A\cos t + B\sin t .
> $$
> **Check of factorisation.** $(s-1)(s^2+1) = s^3 + s - s^2 - 1$ ✓.

> [!example] Example 12.9 — Stability of a cubic without computing roots
> A third-order active filter has characteristic polynomial $s^3 + 2s^2 + 3s + K$, where $K > 0$ is a gain. For what $K$ is it stable? What happens at the boundary?
>
> By the cubic criterion, stability requires $2 \cdot 3 > K$, i.e. $0 < K < 6$. At $K = 6$: $s^3 + 2s^2 + 3s + 6 = s^2(s + 2) + 3(s + 2) = (s + 2)(s^2 + 3)$, with roots $-2$ and $\pm j\sqrt3$. The homogeneous solution contains $A\cos\sqrt3\,t + B\sin\sqrt3\,t$, an undamped oscillation at $\sqrt 3 = 1.732$ rad/s. This is exactly how an oscillator is designed: the loop gain is set so that a root pair sits on the imaginary axis. For $K > 6$ the pair crosses into the right half-plane and the oscillation grows until nonlinearity limits it.

> [!tip] Engineering Note
> Designers frequently specify $\zeta$ rather than component values. $\zeta = 1/\sqrt2 \approx 0.707$ gives the second-order Butterworth response (maximally flat frequency response, about 4.3% step overshoot); $\zeta = 1$ is used where overshoot is forbidden, such as in a precision positioning stage or an ADC driver that must not exceed the converter's input range; $\zeta \approx 0.5$–$0.6$ is common in control loops that trade some overshoot for speed. A probe ringing on an oscilloscope is a lightly damped pair formed by ground-lead inductance and input capacitance; the fix (a shorter ground lead or a damping resistor) raises $\zeta$.


---

## 12.3 Natural and forced responses

### Motivation

When a source drives a circuit, the response always contains two ingredients with completely different origins. One is dictated by the *circuit*: the exponentials $e^{s_k t}$ whose rates are the characteristic roots, the same no matter what drives it. The other is dictated by the *source*: a DC input produces a DC component, a sinusoid produces a sinusoid at the same frequency. Engineers call these the **natural response** and the **forced response**. Recognising which is which tells you immediately which part of a waveform you can change by redesigning the circuit and which part only by changing the input.

There is a second, independent way to split the same response, according to *cause*: the part due to energy initially stored in capacitors and inductors, and the part due to the source. Both splittings are useful, they are not the same, and confusing them is a frequent source of errors.

### Natural and forced decomposition

> [!definition] Natural and forced response
> For $P(D)y = r(t)$ with general solution $y = y_h + y_p$:
> - the **natural response** $y_n(t)$ is the homogeneous part, a combination of the characteristic modes $t^k e^{s_k t}$ with constants fixed by the initial conditions;
> - the **forced response** $y_f(t)$ is a particular solution chosen to have the same functional form as the input (constant for constant input, sinusoid for sinusoid, and so on), containing **no** natural modes unless resonance forces them in (Section 12.4).
>
> Then $y = y_n + y_f$.

The forced response is sometimes defined loosely as "any particular solution". The convention adopted here, that $y_f$ is the particular solution built from the input's own functional form, makes $y_f$ unique whenever the input's exponent does not coincide with a characteristic root. For example, for $y' + 2y = 4$ both $y_p = 2$ and $y_p = 2 + 7e^{-2t}$ are particular solutions, but only $y_f = 2$ is the forced response.

> [!theorem] Procedure for the complete response
> 1. Find the characteristic roots and write the natural response with $n$ undetermined constants.
> 2. Find the forced response $y_f$ (Section 12.4).
> 3. Write $y = y_n + y_f$ and impose the $n$ initial conditions **on the complete $y$**, not on $y_n$ alone.

> [!warning] Common Pitfall
> The constants in $y_n$ must be fitted after adding $y_f$. Fitting $y_n$ to the initial conditions first and then adding $y_f$ gives a function that satisfies the ODE but violates the initial conditions, unless $y_f$ and its derivatives happen to vanish at $t = 0$.

### Zero-input and zero-state decomposition

> [!definition] Zero-input and zero-state response
> - The **zero-input response** $y_{zi}$ solves $P(D)y = 0$ with the actual initial conditions: the response to stored energy alone, with the source switched off.
> - The **zero-state response** $y_{zs}$ solves $P(D)y = r(t)$ with **all initial conditions zero**: the response to the source alone, starting from rest.
>
> By linearity $y = y_{zi} + y_{zs}$.

> [!theorem] Relation between the two decompositions
> The zero-input response consists purely of natural modes. The zero-state response contains the whole forced response plus some natural modes:
> $$
> y_{zs} = y_f + y_{n,\text{zs}}, \qquad y_n = y_{zi} + y_{n,\text{zs}} .
> $$
> Thus the natural response is generally **not** the zero-input response; the source itself excites the natural modes when it is switched on.

The zero-state response is the one that a transfer function or a convolution integral computes: it is linear in the input, and doubling the input doubles it. The complete response is **not** linear in the input when initial conditions are nonzero; it is *affine* ($y_{zi}$ is an offset that does not scale). That is why "linear system" in signal-processing usage tacitly assumes the system starts at rest.

### The impulse response

A special natural response deserves a name. Let $h(t)$ be the zero-state response to a unit impulse $r(t) = \delta(t)$. For $t > 0$ the input is zero, so $h$ is purely natural; the impulse merely sets its initial conditions.

> [!theorem] Impulse response of $P(D)y = \delta(t)$
> With $P(s) = a_n s^n + \dots + a_0$, the impulse response is the solution of the homogeneous equation for $t > 0$ with
> $$
> h(0^+) = h'(0^+) = \dots = h^{(n-2)}(0^+) = 0, \qquad h^{(n-1)}(0^+) = \frac{1}{a_n},
> $$
> and $h(t) = 0$ for $t < 0$. The zero-state response to any input applied at $t = 0$ is then the convolution
> $$
> y_{zs}(t) = \int_0^t h(t - \tau)\,r(\tau)\,d\tau .
> $$

*Derivation sketch.* Integrate $a_n y^{(n)} + \dots + a_0 y = \delta(t)$ from $0^-$ to $0^+$. If $y$ and its first $n - 2$ derivatives are continuous and $y^{(n-1)}$ has a jump, only the highest term picks up the integral of $\delta$, giving $a_n\big[y^{(n-1)}(0^+) - 0\big] = 1$. The lower-order terms integrate to zero over the vanishing interval because they are bounded. The convolution formula follows from writing $r(t) = \int r(\tau)\delta(t - \tau)\,d\tau$ and using linearity and time invariance: each slice $r(\tau)\,d\tau\,\delta(t - \tau)$ produces $r(\tau)\,d\tau\,h(t - \tau)$, and summing the slices gives the integral. $\square$

The impulse response is the bridge to [Ch. 14](../part-5-transforms/14-laplace-transforms.md) (where $\mathcal{L}\{h\} = 1/P(s)$) and to [Ch. 16](../part-5-transforms/16-fourier-transforms.md) (where its Fourier transform is the frequency response). The convolution theorem is developed there; here it is enough that the zero-state response is fully determined by $h$.

### Worked examples

> [!example] Example 12.10 — First order: both decompositions
> Solve $y' + 2y = 4$, $y(0) = 5$, and give the natural/forced and zero-input/zero-state splits.
>
> **Natural/forced.** Root $s = -2$: $y_n = Ae^{-2t}$. Forced: constant $y_f = K$ with $2K = 4$, $K = 2$. Then $y = 2 + Ae^{-2t}$ and $y(0) = 5$ gives $A = 3$:
> $$
> y = \underbrace{2}_{y_f} + \underbrace{3e^{-2t}}_{y_n}.
> $$
> **Zero-input/zero-state.** $y_{zi}$: $y' + 2y = 0$, $y(0) = 5$, so $y_{zi} = 5e^{-2t}$. $y_{zs}$: $y' + 2y = 4$, $y(0) = 0$, so $y_{zs} = 2(1 - e^{-2t})$. Sum: $5e^{-2t} + 2 - 2e^{-2t} = 2 + 3e^{-2t}$ ✓.
>
> The natural response $3e^{-2t}$ is the zero-input response $5e^{-2t}$ plus the natural part of the zero-state response, $-2e^{-2t}$, exactly as the relation theorem says.

> [!example] Example 12.11 — Second-order step from rest
> Solve $y'' + 3y' + 2y = 10$, $y(0) = 0$, $y'(0) = 0$.
>
> **Natural.** $s^2 + 3s + 2 = (s + 1)(s + 2)$: $y_n = c_1 e^{-t} + c_2 e^{-2t}$.
>
> **Forced.** $y_f = K$, $2K = 10$, $K = 5$.
>
> **Initial conditions on the complete response.** $y = 5 + c_1 e^{-t} + c_2 e^{-2t}$:
> $$
> \begin{aligned}
> y(0) &= 5 + c_1 + c_2 = 0,\\
> y'(0) &= -c_1 - 2c_2 = 0,
> \end{aligned}
> \qquad\Rightarrow\qquad c_1 = -10,\ c_2 = 5 .
> $$
> $$
> y(t) = 5 - 10e^{-t} + 5e^{-2t}.
> $$
> **Check.** $y(0) = 0$ ✓, $y'(0) = 10 - 10 = 0$ ✓; numerically $y(1.5) = 5 - 2.2313 + 0.2489 = 3.0176$, matching a numerical ODE solver to 10 digits. Since the system starts at rest, this complete response *is* the zero-state response.

> [!example] Example 12.12 — Same system, stored energy and source together
> Solve $y'' + 3y' + 2y = 10$ with $y(0) = 2$, $y'(0) = 1$, and decompose it.
>
> **Complete response.** $y = 5 + c_1 e^{-t} + c_2 e^{-2t}$; $5 + c_1 + c_2 = 2$ and $-c_1 - 2c_2 = 1$ give $c_2 = 2$, $c_1 = -5$:
> $$
> y = 5 - 5e^{-t} + 2e^{-2t}, \qquad y_f = 5, \quad y_n = -5e^{-t} + 2e^{-2t}.
> $$
> **Zero-input.** $c_1 + c_2 = 2$, $-c_1 - 2c_2 = 1$: $c_2 = -3$, $c_1 = 5$, so $y_{zi} = 5e^{-t} - 3e^{-2t}$.
>
> **Zero-state.** From Example 12.11, $y_{zs} = 5 - 10e^{-t} + 5e^{-2t}$.
>
> **Consistency.** $y_{zi} + y_{zs} = 5 - 5e^{-t} + 2e^{-2t}$ ✓, and $y(1.5) = 3.9839$ agrees with a numerical solver. If the source were doubled to $20$, $y_{zs}$ would double but $y_{zi}$ would not, so the complete response would become $10 - 15e^{-t} + 7e^{-2t}$, not twice the original.

```desmos-graph
left=-0.3; right=6; bottom=-1.2; top=5.6
---
y=5-5e^{-x}+2e^{-2x}|blue|0<x<6
y=5e^{-x}-3e^{-2x}|green|dashed|0<x<6
y=5-10e^{-x}+5e^{-2x}|red|dashed|0<x<6
y=5|grey|dotted
```
*Example 12.12: the complete response (blue) is the sum of the zero-input response (green, decays to zero) and the zero-state response (red, rises from rest to the forced value 5).*

> [!example] Example 12.13 — Impulse response and convolution
> Find the impulse response of $y'' + 3y' + 2y = r(t)$ and use convolution to recover the unit-step response.
>
> **Impulse response.** For $t > 0$, $h = c_1 e^{-t} + c_2 e^{-2t}$ with $h(0^+) = 0$, $h'(0^+) = 1/a_2 = 1$: $c_1 + c_2 = 0$, $-c_1 - 2c_2 = 1$, so $c_1 = 1$, $c_2 = -1$:
> $$
> h(t) = \big(e^{-t} - e^{-2t}\big)\,u(t).
> $$
> **Step response by convolution.** With $r = u(t)$,
> $$
> y_{zs}(t) = \int_0^t \big(e^{-(t-\tau)} - e^{-2(t-\tau)}\big)\,d\tau = \big(1 - e^{-t}\big) - \tfrac12\big(1 - e^{-2t}\big) = \tfrac12 - e^{-t} + \tfrac12 e^{-2t}.
> $$
> Multiplying by 10 reproduces Example 12.11 ✓. Note also that the step response is the integral of the impulse response, and $h$ is the derivative of the step response: the derivative of $\tfrac12 - e^{-t} + \tfrac12 e^{-2t}$ is $e^{-t} - e^{-2t}$ ✓. Because $u' = \delta$, linearity and time invariance force the two responses to be related by the same derivative.

> [!example] Example 12.14 — Natural response excited by switching a sinusoid
> An RL branch obeys $L\,i' + R\,i = V_m\cos\omega t$ for $t \ge 0$ with $i(0) = 0$. Show that a natural component appears unless the switching instant is chosen carefully.
>
> With $\tau = L/R$, the forced response (Section 12.4) is $i_f = I_m\cos(\omega t - \theta)$ with $I_m = V_m/\sqrt{R^2 + \omega^2L^2}$ and $\theta = \arctan(\omega L/R)$. The natural response is $Ae^{-t/\tau}$, and $i(0) = 0$ gives $A = -I_m\cos\theta$:
> $$
> i(t) = I_m\big[\cos(\omega t - \theta) - \cos\theta\,e^{-t/\tau}\big].
> $$
> If instead the source is switched on at a phase $\psi$, i.e. $V_m\cos(\omega t + \psi)$, then $A = -I_m\cos(\psi - \theta)$, which vanishes when $\psi = \theta \pm \pi/2$. Switching at that point on the wave produces no transient at all; switching at $\psi = \theta$ produces the largest. This is the origin of the **DC offset** in fault currents and transformer inrush: in a highly inductive circuit ($\theta \to 90^\circ$), switching at a voltage zero gives a natural component as large as the forced amplitude, and the first current peak approaches $2I_m$.

> [!tip] Engineering Note
> On a bench, the zero-input response is what you see when you charge a circuit and then disconnect the source (a "ring-down" or discharge test); it measures the natural frequencies directly. The zero-state response is what a network analyser or a step-response test measures, because the instrument starts from rest. Datasheet specifications such as rise time, overshoot and settling time always refer to the zero-state step response.

---

## 12.4 Complementary and particular solutions

### Motivation and vocabulary

The general solution $y = y_c + y_p$ of $P(D)y = r$ has two parts that mathematicians name differently from engineers:

| Mathematics | Engineering | Found by |
|---|---|---|
| complementary function $y_c$ | natural response $y_n$ (before fitting ICs) | characteristic roots (Section 12.2) |
| particular integral $y_p$ | forced response $y_f$ | this section |

Section 12.2 settled $y_c$ completely. This section gives methods for $y_p$. For the inputs that occur most often in practice (constants, polynomials, exponentials, sinusoids and products of these) the **method of undetermined coefficients** is fast and purely algebraic. For sinusoids specifically, the **complex-exponential method** is faster still and leads straight to phasors. For arbitrary inputs, variation of parameters ([Ch. 11](11-ordinary-differential-equations.md)) and convolution (Section 12.3) always work, and the Laplace transform ([Ch. 14](../part-5-transforms/14-laplace-transforms.md)) automates the whole procedure.

### Why guessing works: the annihilator idea

Undetermined coefficients looks like guesswork, but it rests on a precise fact. The admissible inputs are exactly the functions that are themselves solutions of some constant-coefficient homogeneous equation.

> [!definition] Annihilator
> A constant-coefficient operator $Q(D)$ **annihilates** $r(t)$ if $Q(D)\,r = 0$. Examples:
>
> | $r(t)$ | Annihilator $Q(D)$ |
> |---|---|
> | $t^k$ (or any polynomial of degree $k$) | $D^{k+1}$ |
> | $e^{at}$ | $D - a$ |
> | $t^k e^{at}$ | $(D - a)^{k+1}$ |
> | $\cos\omega t$, $\sin\omega t$ | $D^2 + \omega^2$ |
> | $e^{\sigma t}\cos\omega t$, $e^{\sigma t}\sin\omega t$ | $(D - \sigma)^2 + \omega^2$ |

If $P(D)y = r$ and $Q(D)r = 0$, apply $Q(D)$ to both sides: $Q(D)P(D)\,y = 0$. So every solution $y$ of the forced equation solves a *higher-order homogeneous* equation whose characteristic polynomial is $Q(s)P(s)$. Its solutions are known from Section 12.2. Those coming from the roots of $P$ are already in $y_c$; the remaining ones, coming from the roots of $Q$, give the form of $y_p$. If a root of $Q$ coincides with a root of $P$, the multiplicities add, and the repeated-root rule automatically introduces extra factors of $t$. That is the whole theory behind the resonance (modification) rule below.

### The method of undetermined coefficients

> [!theorem] Trial forms for $y_p$
> Let $r(t) = e^{at}\big[p_k(t)\cos\omega t + q_k(t)\sin\omega t\big]$, where $p_k$, $q_k$ are polynomials of degree at most $k$ (either may be zero; $\omega = 0$ covers pure exponentials and polynomials, $a = 0$ covers pure sinusoids). Let $m$ be the multiplicity of $s_0 = a + j\omega$ as a root of $P(s)$ ($m = 0$ if it is not a root). Then there is a particular solution of the form
> $$
> y_p = t^m e^{at}\big[A_k(t)\cos\omega t + B_k(t)\sin\omega t\big],
> $$
> with $A_k$, $B_k$ polynomials of degree $k$ whose coefficients are determined by substitution.

Common cases in table form (with $m = 0$, i.e. no resonance):

| Input $r(t)$ | Trial $y_p$ |
|---|---|
| constant $K$ | $A$ |
| $t^k$ | $A_k t^k + \dots + A_1 t + A_0$ |
| $e^{at}$ | $Ae^{at}$ |
| $\cos\omega t$ or $\sin\omega t$ | $A\cos\omega t + B\sin\omega t$ |
| $e^{at}\cos\omega t$ | $e^{at}(A\cos\omega t + B\sin\omega t)$ |
| $t e^{at}$ | $(A_1 t + A_0)e^{at}$ |

> [!theorem] Modification (resonance) rule
> If any term of the trial form is a solution of the homogeneous equation, multiply the **entire** trial form by $t^m$, where $m$ is the multiplicity of the corresponding characteristic root. After the modification no term of $y_p$ duplicates a term of $y_c$.

> [!warning] Common Pitfall
> Three errors account for most failures:
> 1. **Dropping a companion term.** For $r = \sin\omega t$ the trial must contain *both* $\cos\omega t$ and $\sin\omega t$; for $r = t^2$ it must contain $t^2$, $t$ and $1$.
> 2. **Multiplying only part of the trial.** For $r = te^{-t}$ with a simple root $-1$, the trial is $t(A_1 t + A_0)e^{-t}$, not $(A_1 t^2 + A_0)e^{-t}$.
> 3. **Ignoring resonance.** Substituting $Ae^{at}$ when $P(a) = 0$ gives the impossible equation $0 = 1$; that is the signal to multiply by $t$.

### A closed-form shortcut for exponential inputs

For $r(t) = e^{at}$ the eigenfunction identity gives the answer without solving any equations:

> [!theorem] Exponential input formula
> If $P(a) \ne 0$, then $P(D)y = e^{at}$ has the particular solution
> $$
> y_p = \frac{e^{at}}{P(a)}.
> $$
> If $a$ is a root of multiplicity $m$, then
> $$
> y_p = \frac{t^m e^{at}}{P^{(m)}(a)}.
> $$

*Proof.* The first case is immediate from $P(D)e^{at} = P(a)e^{at}$. For the second, write $P(s) = (s - a)^m Q(s)$ with $Q(a) \ne 0$. By the shift theorem,
$$
P(D)\big[e^{at}t^m\big] = e^{at}P(D + a)\,t^m = e^{at}\,Q(D + a)\,D^m t^m = e^{at}\,Q(D + a)\,m! = m!\,Q(a)\,e^{at},
$$
since $D^m t^m = m!$ is a constant and $Q(D + a)$ acting on a constant multiplies it by $Q(a)$. Finally, differentiating $P(s) = (s - a)^m Q(s)$ $m$ times and setting $s = a$ gives $P^{(m)}(a) = m!\,Q(a)$. $\square$

For a constant input $K$ (the case $a = 0$), the formula gives $y_p = K/P(0) = K/a_0$: the DC gain is $1/a_0$. This is the formal version of the circuit-analysis rule "at DC, capacitors are open circuits and inductors are short circuits", which amounts to setting every derivative to zero.

### The complex-exponential (phasor) method

For a sinusoidal input $r(t) = R_m\cos(\omega t + \psi)$ write $r = \operatorname{Re}\{\mathbf{R}\,e^{j\omega t}\}$ with complex amplitude (phasor) $\mathbf{R} = R_m e^{j\psi}$. Solve the complex equation $P(D)z = \mathbf{R}e^{j\omega t}$ with the exponential formula, $z = \mathbf{R}e^{j\omega t}/P(j\omega)$, and take the real part.

> [!theorem] Sinusoidal forced response
> If $P(j\omega) \ne 0$ and $P$ has real coefficients, the forced response to $R_m\cos(\omega t + \psi)$ is
> $$
> y_f(t) = \operatorname{Re}\left\{\frac{\mathbf{R}}{P(j\omega)}e^{j\omega t}\right\} = \frac{R_m}{\lvert P(j\omega)\rvert}\cos\big(\omega t + \psi - \angle P(j\omega)\big).
> $$
> The output phasor is $\mathbf{Y} = \mathbf{R}/P(j\omega)$.

The quantity $H(j\omega) = 1/P(j\omega)$ (more generally $N(j\omega)/P(j\omega)$ when derivatives of the input appear on the right-hand side) is the **frequency response**. Section 12.5 studies it as a function of $\omega$. The method reduces sinusoidal steady-state analysis to complex arithmetic, which is why AC circuit theory uses impedances $R$, $j\omega L$ and $1/(j\omega C)$: each is the eigenvalue of an element law acting on $e^{j\omega t}$.

### Superposition of inputs

If $r = r_1 + r_2 + \dots$, find a particular solution for each term separately and add them. Each term may need its own resonance modification. This is valid only for the *particular* parts; the complementary function is written once, with $n$ constants, and fitted to the initial conditions at the end.

### When undetermined coefficients fails

The method requires $r$ to have a finite-dimensional family of derivatives. It fails for inputs such as $\tan t$, $1/t$, $\ln t$, $e^{t^2}$ or a general measured waveform. Then use:

- **variation of parameters**, which works for any continuous $r$ once $y_c$ is known (derived in [Ch. 11](11-ordinary-differential-equations.md));
- **convolution** with the impulse response (Section 12.3);
- **Laplace transforms** for piecewise inputs such as pulses, ramps and delayed steps ([Ch. 14](../part-5-transforms/14-laplace-transforms.md));
- **numerical integration** for anything else ([Ch. 22](../part-7-numerical-computational/22-numerical-methods.md)).

### Worked examples

> [!example] Example 12.15 — Polynomial input
> Find a particular solution of $y'' + 3y' + 2y = 4t^2$.
>
> No root of $s^2 + 3s + 2$ is zero, so $m = 0$ and the trial is $y_p = At^2 + Bt + C$. Then $y_p' = 2At + B$, $y_p'' = 2A$, and
> $$
> 2A + 3(2At + B) + 2(At^2 + Bt + C) = 2A\,t^2 + (6A + 2B)\,t + (2A + 3B + 2C).
> $$
> Matching coefficients with $4t^2 + 0t + 0$:
> $$
> 2A = 4,\qquad 6A + 2B = 0,\qquad 2A + 3B + 2C = 0 \quad\Rightarrow\quad A = 2,\ B = -6,\ C = 7 .
> $$
> $$
> y_p = 2t^2 - 6t + 7 .
> $$
> **Check.** $y_p'' = 4$, $3y_p' = 12t - 18$, $2y_p = 4t^2 - 12t + 14$; sum $= 4t^2$ ✓. Note that the constant and linear terms are essential even though the input contains only $t^2$.

> [!example] Example 12.16 — Sinusoid by the complex-exponential method
> Find the forced response of $y'' + 2y' + 5y = 10\cos t$.
>
> $P(s) = s^2 + 2s + 5$ and $P(j1) = -1 + 2j + 5 = 4 + 2j$, which is nonzero (the roots are $-1 \pm 2j$). The output phasor is
> $$
> \mathbf{Y} = \frac{10}{4 + 2j} = \frac{10(4 - 2j)}{20} = 2 - j .
> $$
> $$
> y_f = \operatorname{Re}\{(2 - j)(\cos t + j\sin t)\} = 2\cos t + \sin t = \sqrt5\cos(t - 0.4636).
> $$
> **Check by substitution.** $y_f' = -2\sin t + \cos t$, $y_f'' = -2\cos t - \sin t$. Then $y_f'' + 2y_f' + 5y_f = (-2 + 2 + 10)\cos t + (-1 - 4 + 5)\sin t = 10\cos t$ ✓. The amplitude $\sqrt5 = 10/\lvert 4 + 2j\rvert = 10/\sqrt{20}$ ✓, and the phase lag is $\arctan(2/4) = 26.57^\circ$.
>
> The same result by undetermined coefficients requires solving $4A + 2B = 10$, $-2A + 4B = 0$ (from the $\cos$ and $\sin$ coefficients), giving $A = 2$, $B = 1$, the same answer with more algebra.

> [!example] Example 12.17 — Pure resonance
> Solve $y'' + 4y = 3\cos 2t$, $y(0) = 0$, $y'(0) = 0$.
>
> **Complementary.** Roots $\pm 2j$: $y_c = c_1\cos 2t + c_2\sin 2t$.
>
> **Particular.** The input frequency $2$ coincides with a simple root ($m = 1$), so use the formula with $P'(s) = 2s$: the complex equation $P(D)z = 3e^{j2t}$ has
> $$
> z_p = \frac{3\,t\,e^{j2t}}{P'(2j)} = \frac{3t\,e^{j2t}}{4j} = -\frac{3j}{4}\,t\,(\cos 2t + j\sin 2t),
> $$
> whose real part is $y_p = \tfrac34\,t\sin 2t$.
>
> **Initial conditions.** $y(0) = c_1 = 0$; $y'(0) = 2c_2 + \big[\tfrac34\sin 2t + \tfrac32 t\cos 2t\big]_{t=0} = 2c_2 = 0$. So
> $$
> y(t) = \tfrac34\,t\sin 2t .
> $$
> **Check.** $y'' = 3\cos 2t - 3t\sin 2t$, so $y'' + 4y = 3\cos 2t - 3t\sin 2t + 3t\sin 2t = 3\cos 2t$ ✓. The amplitude grows linearly without bound: a lossless resonator driven at its natural frequency absorbs energy every cycle. Any real resonator has some damping, which caps the amplitude at a finite value (Section 12.5).

> [!example] Example 12.18 — Mixed input with resonance and a full IVP
> Solve $y'' - y' - 2y = e^{2t} + 3$, $y(0) = 0$, $y'(0) = 0$.
>
> **Complementary.** $s^2 - s - 2 = (s - 2)(s + 1)$: $y_c = c_1 e^{2t} + c_2 e^{-t}$.
>
> **Particular, term by term.**
> - For $e^{2t}$: $a = 2$ is a simple root, $P'(s) = 2s - 1$, $P'(2) = 3$, so $y_{p1} = \tfrac13 t e^{2t}$.
> - For $3$: $y_{p2} = 3/P(0) = 3/(-2) = -\tfrac32$.
>
> **Initial conditions.** $y = c_1 e^{2t} + c_2 e^{-t} + \tfrac13 te^{2t} - \tfrac32$:
> $$
> \begin{aligned}
> y(0) &= c_1 + c_2 - \tfrac32 = 0,\\
> y'(0) &= 2c_1 - c_2 + \tfrac13 = 0,
> \end{aligned}
> \qquad\Rightarrow\qquad c_1 = \tfrac{7}{18},\quad c_2 = \tfrac{10}{9}.
> $$
> $$
> y(t) = \tfrac{7}{18}e^{2t} + \tfrac{10}{9}e^{-t} + \tfrac13 te^{2t} - \tfrac32 .
> $$
> **Check.** At $t = 1$ the formula gives $4.24530$, agreeing with a numerical solution to nine significant figures. Verifying $y_{p1}$ directly: with $w = te^{2t}$, $w' = (1 + 2t)e^{2t}$, $w'' = (4 + 4t)e^{2t}$, and $w'' - w' - 2w = (4 + 4t - 1 - 2t - 2t)e^{2t} = 3e^{2t}$ ✓, so $w/3$ is correct.

> [!example] Example 12.19 — Exponential times polynomial via the shift theorem
> Find a particular solution of $y'' + y = t e^{t}$.
>
> Write $y_p = e^{t}f(t)$. By the shift theorem, $(D^2 + 1)\big[e^t f\big] = e^t\big[(D + 1)^2 + 1\big]f = e^t\big(D^2 + 2D + 2\big)f$. We need $(D^2 + 2D + 2)f = t$. Try $f = At + B$: $2A + 2At + 2B = t$, so $A = \tfrac12$, $B = -\tfrac12$:
> $$
> y_p = \tfrac12(t - 1)e^{t}.
> $$
> **Check.** $y_p' = \tfrac12 t e^t$, $y_p'' = \tfrac12(t + 1)e^t$; $y_p'' + y_p = \tfrac12(t + 1 + t - 1)e^t = te^t$ ✓. The shift theorem reduces an exponential-times-polynomial problem to a pure polynomial problem, which is often less error-prone than substituting $(At + B)e^t$ directly.

> [!example] Example 12.20 — Damped sinusoidal input
> Find the forced response of $y'' + 4y' + 13y = e^{-2t}\cos 3t$, and explain the result.
>
> The characteristic roots are $-2 \pm 3j$, and the input is $\operatorname{Re}\{e^{(-2 + 3j)t}\}$ with exponent $s_0 = -2 + 3j$, a simple root. So $m = 1$ and, with $P'(s) = 2s + 4$, $P'(s_0) = 6j$:
> $$
> y_p = \operatorname{Re}\left\{\frac{t\,e^{(-2+3j)t}}{6j}\right\} = \operatorname{Re}\left\{-\frac{j}{6}\,t\,e^{-2t}(\cos 3t + j\sin 3t)\right\} = \frac{1}{6}\,t\,e^{-2t}\sin 3t .
> $$
> **Check.** Using the shift theorem with $y_p = e^{-2t}g$: $P(D - 2) = (D - 2)^2 + 4(D - 2) + 13 = D^2 + 9$, so we need $g'' + 9g = \cos 3t$, which by Example 12.17's pattern has $g = \tfrac16 t\sin 3t$ ✓. The input is itself a natural mode of the system, so it resonates even though it decays.

---

## 12.5 Transient and steady-state behavior

### Motivation

Switch on a sinusoidal generator connected to a filter and watch the output on an oscilloscope. For a short time the waveform is irregular; then it settles into a clean sinusoid at the generator frequency, with a fixed amplitude and phase. The irregular part is the **transient**; the clean sinusoid is the **steady state**. Almost all AC circuit analysis, filter design and frequency-response measurement concerns the steady state, but the transient decides how long one must wait before the steady-state formulas apply, how large the switching surges are, and whether the steady state exists at all.

### Definitions

> [!definition] Transient and steady-state response
> Suppose the response of $P(D)y = r(t)$ can be written $y = y_{tr} + y_{ss}$ where
> $$
> \lim_{t \to \infty} y_{tr}(t) = 0
> $$
> and $y_{ss}$ persists (for example constant, periodic, or growing like the input). Then $y_{tr}$ is the **transient response** and $y_{ss}$ is the **steady-state response**.

> [!theorem] Transient = natural, steady state = forced (for stable systems)
> If $P(D)y = 0$ is asymptotically stable (all characteristic roots satisfy $\operatorname{Re}\{s_k\} < 0$), then for any input of the undetermined-coefficients class:
> 1. the natural response $y_n$ decays to zero for every set of initial conditions, so $y_{tr} = y_n$ and $y_{ss} = y_f$;
> 2. the steady state is **independent of the initial conditions**;
> 3. if the input is bounded, the steady state is bounded (bounded-input, bounded-output stability).

*Reason.* Every term of $y_n$ has the form $t^k e^{s_k t}$ with $\operatorname{Re}\{s_k\} < 0$ and tends to zero. Initial conditions enter only through the constants of $y_n$. If the input contains no natural modes (stable roots never coincide with a pure sinusoid or constant, which have exponents on the imaginary axis), $y_f$ has the same bounded form as the input. $\square$

When the system is not asymptotically stable, the split fails in instructive ways:

- **Unstable** ($\operatorname{Re}\{s_k\} > 0$ for some root): the natural response grows; the "steady state" of the formula $y_f$ exists on paper but is never observed, because any perturbation excites the growing mode. An unstable amplifier oscillates or latches instead of amplifying.
- **Marginally stable** (simple roots on the imaginary axis): the natural response never decays, so the initial conditions are remembered forever. A lossless LC tank rings indefinitely; an integrator ($s = 0$) holds its initial value.
- **Resonance with a marginal root**: a constant input to an integrator, or a sinusoid at an LC tank's natural frequency, produces a growing $y_f$ ($t$ or $t\sin\omega t$, Example 12.17) and there is no bounded steady state.

### How fast does the transient die?

The slowest-decaying mode dominates the late transient. With the root closest to the imaginary axis at real part $-\sigma_{\min}$, the transient is bounded by a multiple of $e^{-\sigma_{\min}t}$ for large $t$ (times a polynomial factor if that root is repeated).

> [!definition] Time constant and settling time
> - For a first-order mode $e^{-t/\tau}$, the **time constant** is $\tau = 1/\sigma$.
> - The **2% settling time** $t_s$ is the time after which the response stays within 2% of its final value. For a first-order system $t_s = \tau\ln 50 \approx 3.9\tau$, usually rounded to $4\tau$. For an underdamped second-order system the envelope $e^{-\zeta\omega_n t}/\sqrt{1 - \zeta^2}$ gives
> $$
> t_s \approx \frac{\ln\!\big(50/\sqrt{1 - \zeta^2}\big)}{\zeta\omega_n} \approx \frac{4}{\zeta\omega_n}.
> $$

For the 1% and 0.1% criteria replace $\ln 50 = 3.91$ by $\ln 100 = 4.61$ and $\ln 1000 = 6.91$. A rule worth memorising: **each time constant buys about 8.7 dB (a factor $e$)**, so reaching 16-bit accuracy ($2^{-16} \approx 1.5\times10^{-5}$) takes about $11\tau$.

### Sinusoidal steady state and frequency response

For a stable system driven by $\cos\omega t$, Section 12.4 gave the steady state directly from the complex number $P(j\omega)$. In general, when the input enters through derivatives as well, the ODE reads $P(D)y = N(D)r$, and the same argument gives:

> [!theorem] Frequency response
> For the stable LTI equation $P(D)y = N(D)r$ with input $r = R_m\cos(\omega t + \psi)$, the steady-state output is
> $$
> y_{ss}(t) = R_m\,\lvert H(j\omega)\rvert\cos\big(\omega t + \psi + \angle H(j\omega)\big), \qquad H(j\omega) = \frac{N(j\omega)}{P(j\omega)} .
> $$
> $\lvert H(j\omega)\rvert$ is the **gain** (magnitude response) and $\angle H(j\omega)$ is the **phase shift**. In phasor form $\mathbf{Y} = H(j\omega)\,\mathbf{R}$.

The frequency response is the transfer function $H(s) = N(s)/P(s)$ of [Ch. 14](../part-5-transforms/14-laplace-transforms.md) evaluated on the imaginary axis. Here it arises without any transform: it is simply the eigenvalue of the system for the eigenfunction $e^{j\omega t}$.

**First-order low-pass.** For $\tau y' + y = r$, $H(j\omega) = 1/(1 + j\omega\tau)$:
$$
\lvert H\rvert = \frac{1}{\sqrt{1 + (\omega\tau)^2}}, \qquad \angle H = -\arctan(\omega\tau).
$$
At the **corner frequency** $\omega_c = 1/\tau$ the gain is $1/\sqrt2$ ($-3$ dB) and the phase is $-45^\circ$; well above it the gain falls as $1/(\omega\tau)$, i.e. $-20$ dB per decade.

**Second-order low-pass.** For $y'' + 2\zeta\omega_n y' + \omega_n^2 y = \omega_n^2 r$, with normalised frequency $\nu = \omega/\omega_n$,
$$
H(j\omega) = \frac{1}{1 - \nu^2 + j\,2\zeta\nu}, \qquad \lvert H\rvert = \frac{1}{\sqrt{(1 - \nu^2)^2 + (2\zeta\nu)^2}}, \qquad \angle H = -\operatorname{atan2}\big(2\zeta\nu,\ 1 - \nu^2\big).
$$

> [!theorem] Resonant peak of a second-order system
> For $0 < \zeta < 1/\sqrt2$, $\lvert H(j\omega)\rvert$ has a maximum at the **resonant frequency** and with the **peak value**
> $$
> \omega_r = \omega_n\sqrt{1 - 2\zeta^2}, \qquad M_r = \frac{1}{2\zeta\sqrt{1 - \zeta^2}} .
> $$
> For $\zeta \ge 1/\sqrt2$ there is no peak; the gain decreases monotonically. At $\omega = \omega_n$ the phase is exactly $-90^\circ$ and the gain is $1/(2\zeta)$ for every $\zeta$.

*Derivation.* Maximising $\lvert H\rvert$ means minimising $g(\nu^2) = (1 - \nu^2)^2 + 4\zeta^2\nu^2$. With $x = \nu^2$, $g'(x) = -2(1 - x) + 4\zeta^2 = 0$ gives $x = 1 - 2\zeta^2$, which is positive only if $\zeta^2 < 1/2$. Substituting, $g = 4\zeta^4 + 4\zeta^2(1 - 2\zeta^2) = 4\zeta^2(1 - \zeta^2)$, so $M_r = 1/\sqrt{g} = 1/\big(2\zeta\sqrt{1 - \zeta^2}\big)$. $\square$

For light damping $M_r \approx 1/(2\zeta) = Q$, the **quality factor**, and the half-power ($-3$ dB) bandwidth of the peak is approximately $2\zeta\omega_n = \omega_n/Q$. For $\zeta = 0.1$ the exact peak is $M_r = 5.025$ at $\nu = 0.990$, and the exact half-power bandwidth is $0.202\,\omega_n$, close to the approximation $0.2\,\omega_n$.

```desmos-graph
left=0; right=2.5; bottom=0; top=5.5
---
y=\frac{1}{\sqrt{(1-x^2)^2+(0.2x)^2}}|red|0<x<2.5
y=\frac{1}{\sqrt{(1-x^2)^2+(0.6x)^2}}|orange|0<x<2.5
y=\frac{1}{\sqrt{(1-x^2)^2+(1.414x)^2}}|blue|0<x<2.5
y=\frac{1}{\sqrt{(1-x^2)^2+(2x)^2}}|green|0<x<2.5
(0.99,5.025)|label:ζ = 0.1 peak
x=1|grey|dotted
```
*Magnitude response $\lvert H\rvert$ against $\nu = \omega/\omega_n$ for $\zeta = 0.1$ (red), $0.3$ (orange), $0.707$ (blue, maximally flat) and $1$ (green). Light damping produces a tall, narrow resonant peak just below $\omega_n$.*

### Steady state for periodic inputs

A periodic input can be expanded in a Fourier series ([Ch. 15](../part-5-transforms/15-fourier-series.md)), $r(t) = c_0 + \sum_k \lvert c_k\rvert\cos(k\omega_0 t + \psi_k)$. By superposition the steady state is

$$
y_{ss}(t) = H(0)\,c_0 + \sum_k \lvert c_k\rvert\,\lvert H(jk\omega_0)\rvert\cos\big(k\omega_0 t + \psi_k + \angle H(jk\omega_0)\big).
$$

Each harmonic is scaled and shifted by the frequency response at its own frequency. A low-pass system attenuates the high harmonics, which is why a square wave through an RC filter comes out rounded. A lightly damped system can selectively amplify one harmonic that falls near its resonant peak.

### Beats: the approach to resonance

For an undamped system driven near, but not at, its natural frequency, the "transient" never decays, and the sum of natural and forced responses produces **beats**. For $y'' + \omega_0^2 y = \cos\omega t$ with $y(0) = y'(0) = 0$ and $\omega \ne \omega_0$:

$$
y = \frac{\cos\omega t - \cos\omega_0 t}{\omega_0^2 - \omega^2} = \frac{2}{\omega_0^2 - \omega^2}\,\sin\!\left(\frac{(\omega_0 - \omega)t}{2}\right)\sin\!\left(\frac{(\omega_0 + \omega)t}{2}\right),
$$

using $\cos A - \cos B = -2\sin\frac{A+B}{2}\sin\frac{A-B}{2}$. The fast oscillation at the average frequency is modulated by a slow envelope at half the difference frequency. As $\omega \to \omega_0$, the envelope period grows without bound and, near $t = 0$, $y \to \frac{t\sin\omega_0 t}{2\omega_0}$, the resonant solution of Example 12.17. With damping, the beat pattern is a transient: it decays and leaves the steady-state sinusoid.

```desmos-graph
left=0; right=130; bottom=-11.5; top=11.5
---
y=\frac{\cos(0.9x)-\cos(x)}{0.19}|blue|0<x<130
y=\frac{2}{0.19}\sin(0.05x)|red|dashed|0<x<130
y=-\frac{2}{0.19}\sin(0.05x)|red|dashed|0<x<130
```
*Beats for $\omega_0 = 1$, $\omega = 0.9$: the response (blue) oscillates near $0.95$ rad/s inside a slow envelope $\pm 10.5\sin(0.05t)$ (red), with beat period $2\pi/0.1 = 62.8$.*

### Worked examples

> [!example] Example 12.21 — RC low-pass switched onto a sinusoid
> Solve $\tau y' + y = \cos\omega t$, $y(0) = 0$, identify transient and steady state, and evaluate for $\omega\tau = 1$.
>
> **Steady state.** $H(j\omega) = 1/(1 + j\omega\tau)$:
> $$
> y_{ss} = \frac{1}{\sqrt{1 + \omega^2\tau^2}}\cos\big(\omega t - \arctan\omega\tau\big) = \frac{\cos\omega t + \omega\tau\sin\omega t}{1 + \omega^2\tau^2}.
> $$
> **Transient.** $y_{tr} = Ae^{-t/\tau}$ with $y(0) = 0$: $A = -y_{ss}(0) = -1/(1 + \omega^2\tau^2)$.
> $$
> y(t) = \frac{\cos\omega t + \omega\tau\sin\omega t - e^{-t/\tau}}{1 + \omega^2\tau^2}.
> $$
> **At $\omega\tau = 1$.** $y = \tfrac12\big(\cos\omega t + \sin\omega t - e^{-t/\tau}\big)$, with steady-state amplitude $0.7071$ ($-3$ dB) and lag $45^\circ$. Check: $y(0) = \tfrac12(1 + 0 - 1) = 0$ ✓. After $4\tau$ the transient is $0.5e^{-4} = 0.009$, under 1.3% of the steady-state amplitude.

```desmos-graph
left=0; right=16; bottom=-0.9; top=0.9
---
y=0.5\cos(x)+0.5\sin(x)-0.5e^{-x}|blue|0<x<16
y=0.5\cos(x)+0.5\sin(x)|grey|dashed|0<x<16
y=-0.5e^{-x}|red|dotted|0<x<16
y=\cos(x)|green|dotted|0<x<16
```
*Example 12.21 with $\tau = 1$, $\omega = 1$: input $\cos t$ (green), complete response (blue), steady state (grey dashed) and transient (red). The output merges with the steady state within a few time constants, attenuated by $0.707$ and lagging by $45^\circ$.*

> [!example] Example 12.22 — Second-order: transient, steady state and settling
> Solve $y'' + 2y' + 5y = 10\cos t$, $y(0) = 0$, $y'(0) = 0$, and estimate when the transient is negligible.
>
> **Steady state** (Example 12.16): $y_{ss} = 2\cos t + \sin t$.
>
> **Transient.** Roots $-1 \pm 2j$: $y_{tr} = e^{-t}(A\cos 2t + B\sin 2t)$. Conditions on the complete response:
> $$
> y(0) = 2 + A = 0 \Rightarrow A = -2, \qquad y'(0) = 1 + (-A + 2B) = 0 \Rightarrow B = -\tfrac32 .
> $$
> $$
> y(t) = 2\cos t + \sin t - e^{-t}\left(2\cos 2t + \tfrac32\sin 2t\right).
> $$
> **Check.** At $t = 3$ the formula gives $-1.91361$, matching a numerical solver. The transient amplitude is $\sqrt{4 + 2.25}\,e^{-t} = 2.5e^{-t}$; it falls below 2% of the steady amplitude $\sqrt5 = 2.236$ when $2.5e^{-t} < 0.0447$, i.e. $t > 4.0$. The transient oscillates at $2$ rad/s (the natural frequency $\omega_d$), not at the input frequency $1$ rad/s: natural modes always ring at their own frequency.

> [!example] Example 12.23 — Harmonics of a periodic input
> The input $r(t) = \cos t + \cos 3t$ (the first two terms of a crude approximation to a pulse train) drives $y' + y = r$. Find the steady state.
>
> With $\tau = 1$: at $\omega = 1$, $\lvert H\rvert = 1/\sqrt2 = 0.7071$, $\angle H = -45^\circ$; at $\omega = 3$, $\lvert H\rvert = 1/\sqrt{10} = 0.3162$, $\angle H = -71.57^\circ$.
> $$
> y_{ss} = 0.7071\cos(t - 45^\circ) + 0.3162\cos(3t - 71.57^\circ).
> $$
> The third harmonic is attenuated more (ratio $0.447$) and delayed by a different time: $45^\circ$ at 1 rad/s is a delay of $0.785$ s, whereas $71.57^\circ$ at 3 rad/s is a delay of $0.416$ s. Because different harmonics are delayed differently, the output waveform has a different *shape* from the input. This is phase distortion, which a linear-phase filter avoids.

> [!example] Example 12.24 — No steady state: unstable and marginal systems
> (a) $y'' - y = \cos t$; (b) $y'' + y' = 1$.
>
> **(a)** Roots $\pm 1$. The forced response is $y_f = \cos t/P(j) = \cos t/(-1 - 1) = -\tfrac12\cos t$, but $y_n = c_1 e^{t} + c_2 e^{-t}$. Unless the initial conditions make $c_1 = 0$ exactly (a set of measure zero, destroyed by any noise), $c_1e^t$ dominates. The formula $y_f$ is not a physically meaningful steady state.
>
> **(b)** Roots $0, -1$. The constant input matches the root $s = 0$ ($m = 1$): $y_p = t/P'(0) = t/1 = t$. The response grows as a ramp; the system integrates the constant input. In circuit terms, a capacitor charged by a constant current has no DC steady state. A DC steady state exists only if $P(0) = a_0 \ne 0$.

> [!warning] Common Pitfall
> Phasor and impedance calculations silently assume the system is stable and that enough time has passed. Applying $H(j\omega)$ to an unstable circuit, or reading a "steady-state" measurement before $4$–$5$ time constants of the slowest mode have elapsed, gives wrong answers. In very lightly damped circuits ($Q$ in the hundreds, as in crystal oscillators or high-$Q$ filters) the settling time $t_s \approx 8Q/\omega_n$ can be surprisingly long.

> [!tip] Engineering Note
> A network analyser measures $H(j\omega)$ by sweeping a sinusoid and waiting at each frequency for the transient to die. Its IF bandwidth and dwell time are chosen from exactly the settling-time estimate above; narrow-band high-$Q$ devices need slow sweeps. Similarly, a switching power supply's "load-step response" specification is a transient measurement: the regulator's dominant poles decide how quickly the output voltage recovers after a sudden change in current.

---

## 12.6 Coupled differential equations

### Motivation

Real networks have more than one loop or node, and the unknowns interact. Two mesh currents share a resistor; two coils share magnetic flux through mutual inductance $M$; a motor's current produces torque and its speed produces a back-EMF that opposes the current. Each such interaction gives a system of simultaneous linear ODEs in several unknowns. This section treats them with two tools: **elimination**, which reduces the system to a single higher-order equation of the kind already solved, and **modal decomposition**, which changes variables so that the system splits into independent scalar equations. The general first-order matrix form $\dot{\mathbf{x}} = \mathbf{A}\mathbf{x} + \mathbf{f}$ and its eigenvector solution appear in [Ch. 11](11-ordinary-differential-equations.md); here the emphasis is on equations as they arise from circuits, often second order and with derivatives on several unknowns in one equation.

### Operator-matrix form

> [!definition] Coupled linear system with constant coefficients
> A system of $k$ linear constant-coefficient ODEs in unknowns $y_1, \dots, y_k$ can be written
> $$
> \begin{bmatrix} P_{11}(D) & \cdots & P_{1k}(D) \\ \vdots & & \vdots \\ P_{k1}(D) & \cdots & P_{kk}(D) \end{bmatrix}
> \begin{bmatrix} y_1 \\ \vdots \\ y_k \end{bmatrix}
> =
> \begin{bmatrix} r_1(t) \\ \vdots \\ r_k(t) \end{bmatrix},
> \qquad\text{briefly}\qquad \mathbf{P}(D)\,\mathbf{y} = \mathbf{r},
> $$
> where each $P_{il}(D)$ is a polynomial in $D$. The polynomial $\Delta(s) = \det\mathbf{P}(s)$ is the **characteristic polynomial of the system**.

Because constant-coefficient operators commute (Section 12.1), the entries of $\mathbf{P}(D)$ can be manipulated like the entries of a matrix of polynomials. In particular **Cramer's rule** holds in operator form:

> [!theorem] Elimination by operator determinants
> If $\Delta(s) = \det\mathbf{P}(s)$ is not identically zero, every solution of $\mathbf{P}(D)\mathbf{y} = \mathbf{r}$ satisfies, for each component,
> $$
> \Delta(D)\,y_i = \det\mathbf{P}_i(D),
> $$
> where $\mathbf{P}_i$ is $\mathbf{P}$ with its $i$th column replaced by $\mathbf{r}$ (expanded along that column so that operators act on the $r_l$). Each $y_i$ therefore has natural modes given by the roots of $\Delta(s)$. The general solution of the system contains exactly $\deg\Delta$ independent arbitrary constants.

*Proof sketch.* Multiply the system on the left by the adjugate operator matrix $\operatorname{adj}\mathbf{P}(D)$. Since entries commute, $\operatorname{adj}\mathbf{P}\cdot\mathbf{P} = \Delta(D)\,\mathbf{I}$, giving $\Delta(D)\mathbf{y} = \operatorname{adj}\mathbf{P}(D)\,\mathbf{r}$, which is the column-replacement formula. The count of constants is a standard result (each component individually appears to need $\deg\Delta$ constants, but the original equations tie them together so that only $\deg\Delta$ in total are free). $\square$

> [!warning] Common Pitfall
> After eliminating, each component $y_i$ satisfies a degree-$\deg\Delta$ equation, so writing each with its own $\deg\Delta$ constants gives $k\cdot\deg\Delta$ constants, far too many. The surplus constants are **not** free: substitute back into one of the original equations to express them in terms of the others, or, more simply, solve for one unknown and *compute* the others from the original equations (as in Example 12.25). Also: if $\Delta(s) \equiv 0$ the system is degenerate (dependent equations or a redundant variable) and must be reformulated.

The degree of $\Delta$ equals the number of independent energy-storage elements, i.e. the order of the system. For two coupled inductors, $\Delta$ has degree $2$ unless the coupling is perfect ($M^2 = L_1L_2$), in which case the two inductors behave as one and the degree drops to $1$.

### Modal decomposition

The second method changes variables. Suppose the system can be written $\dot{\mathbf{x}} = \mathbf{A}\mathbf{x}$ (first order) or $\mathbf{M}\ddot{\mathbf{x}} + \mathbf{K}\mathbf{x} = \mathbf{0}$ (second order, undamped, as for LC ladders or spring–mass chains). If $\mathbf{A}$ is diagonalisable, $\mathbf{A} = \mathbf{V}\boldsymbol{\Lambda}\mathbf{V}^{-1}$ ([Ch. 10](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md)), then the **modal coordinates** $\mathbf{q} = \mathbf{V}^{-1}\mathbf{x}$ satisfy

$$
\dot{\mathbf{q}} = \boldsymbol{\Lambda}\mathbf{q}, \qquad\text{i.e.}\qquad \dot q_k = \lambda_k q_k \quad\text{independently for each } k .
$$

The system has split into scalar first-order equations, each with solution $q_k = q_k(0)e^{\lambda_k t}$, and $\mathbf{x} = \mathbf{V}\mathbf{q} = \sum_k q_k(0)\,\mathbf{v}_k e^{\lambda_k t}$. Each eigenvector $\mathbf{v}_k$ is a **mode shape**: the fixed pattern of relative amplitudes in which the unknowns move together at the rate $\lambda_k$.

For the second-order form $\ddot{\mathbf{x}} + \mathbf{K}\mathbf{x} = \mathbf{0}$ (taking $\mathbf{M} = \mathbf{I}$, or after scaling) with $\mathbf{K}$ symmetric positive definite, the eigenvalues $\kappa_k$ of $\mathbf{K}$ are positive and the modes are **normal modes** $\mathbf{v}_k\cos(\omega_k t + \varphi_k)$ with $\omega_k = \sqrt{\kappa_k}$. Symmetric matrices have orthogonal eigenvectors, which is what makes the decoupling exact and the mode shapes physically distinct.

With forcing, the same change of variables gives $\dot q_k = \lambda_k q_k + g_k(t)$ with $\mathbf{g} = \mathbf{V}^{-1}\mathbf{f}$: each mode is driven only by its own projection of the input. A mode whose projection is zero is not excited at all.

### Sinusoidal steady state of coupled systems

For a stable coupled system driven sinusoidally at frequency $\omega$, replace every $D$ by $j\omega$: the phasor amplitudes satisfy the complex linear algebraic system

$$
\mathbf{P}(j\omega)\,\mathbf{Y} = \mathbf{R},
$$

solvable by Cramer's rule or Gaussian elimination ([Ch. 8](../part-3-linear-algebra/08-matrices-and-linear-systems.md)). This is exactly mesh or nodal analysis with impedances. Its determinant $\det\mathbf{P}(j\omega) = \Delta(j\omega)$ is small near the natural frequencies (resonances), while a numerator determinant can vanish at some frequency, making one output zero: an **antiresonance** or transmission zero.

### Worked examples

> [!example] Example 12.25 — Elimination for a first-order pair
> Solve $x' = -2x + y$, $y' = x - 2y$, with $x(0) = 1$, $y(0) = 0$.
>
> **Eliminate $y$.** From the first equation $y = x' + 2x$. Substituting into the second: $x'' + 2x' = x - 2(x' + 2x)$, i.e.
> $$
> x'' + 4x' + 3x = 0, \qquad s^2 + 4s + 3 = (s + 1)(s + 3).
> $$
> (Operator check: $\Delta(s) = \det\begin{bmatrix} s + 2 & -1 \\ -1 & s + 2\end{bmatrix} = (s + 2)^2 - 1 = s^2 + 4s + 3$ ✓.)
>
> **Initial data for $x$.** $x(0) = 1$ and $x'(0) = -2x(0) + y(0) = -2$. With $x = c_1 e^{-t} + c_2 e^{-3t}$: $c_1 + c_2 = 1$, $-c_1 - 3c_2 = -2$, so $c_1 = c_2 = \tfrac12$.
>
> **Recover $y$ from the original equation** (no new constants):
> $$
> x = \tfrac12\big(e^{-t} + e^{-3t}\big), \qquad y = x' + 2x = \tfrac12\big(e^{-t} - e^{-3t}\big).
> $$
> **Modal view.** $\mathbf{A} = \begin{bmatrix} -2 & 1 \\ 1 & -2\end{bmatrix}$ has eigenpairs $\lambda = -1$, $\mathbf{v} = (1, 1)$ and $\lambda = -3$, $\mathbf{v} = (1, -1)$. The solution is $\tfrac12(1,1)e^{-t} + \tfrac12(1,-1)e^{-3t}$: the sum $x + y = e^{-t}$ decays slowly and the difference $x - y = e^{-3t}$ decays fast. Check: $y(0) = 0$ ✓, and $y' = \tfrac12(-e^{-t} + 3e^{-3t})$ equals $x - 2y = \tfrac12(e^{-t} + e^{-3t}) - (e^{-t} - e^{-3t}) = \tfrac12(-e^{-t} + 3e^{-3t})$ ✓.

> [!example] Example 12.26 — Mutually coupled coils with a step input
> Two coils with $L_1 = L_2 = 1\ \mathrm{H}$, mutual inductance $M = 0.5\ \mathrm{H}$ and winding resistances $R_1 = R_2 = 1\ \Omega$ are coupled. A $1\ \mathrm{V}$ step is applied to coil 1 at $t = 0$; coil 2 is short-circuited. Both currents start at zero. Find $i_1(t)$ and $i_2(t)$.
>
> **Equations.**
> $$
> \begin{aligned}
> L_1 i_1' + M i_2' + R_1 i_1 &= 1,\\
> M i_1' + L_2 i_2' + R_2 i_2 &= 0,
> \end{aligned}
> \qquad
> \begin{bmatrix} D + 1 & 0.5D \\ 0.5D & D + 1\end{bmatrix}
> \begin{bmatrix} i_1 \\ i_2 \end{bmatrix} = \begin{bmatrix} 1 \\ 0 \end{bmatrix}.
> $$
> **Characteristic polynomial.** $\Delta(s) = (s + 1)^2 - 0.25s^2 = 0.75s^2 + 2s + 1$, roots $s = (-2 \pm 1)/1.5$: $s_1 = -\tfrac23$, $s_2 = -2\ \mathrm{s^{-1}}$.
>
> **Cramer.** $\Delta(D)\,i_1 = (D + 1)[1] - 0.5D[0] = 1$ and $\Delta(D)\,i_2 = (D + 1)[0] - 0.5D[1] = 0$. So the forced values are $i_{1f} = 1/\Delta(0) = 1$ and $i_{2f} = 0$, as expected (at DC the coils are shorts and only $R_1$ limits $i_1$).
>
> **Initial derivatives.** At $t = 0^+$, $\begin{bmatrix}1 & 0.5\\0.5 & 1\end{bmatrix}\mathbf{i}' = \begin{bmatrix}1\\0\end{bmatrix}$, so $i_1'(0) = \tfrac43$, $i_2'(0) = -\tfrac23\ \mathrm{A/s}$.
>
> **Fit.** $i_1 = 1 + a e^{-2t/3} + b e^{-2t}$ with $1 + a + b = 0$, $-\tfrac23 a - 2b = \tfrac43$: $a = b = -\tfrac12$. For $i_2 = c e^{-2t/3} + d e^{-2t}$: $c + d = 0$, $-\tfrac23 c - 2d = -\tfrac23$: $c = -\tfrac12$, $d = \tfrac12$.
> $$
> i_1 = 1 - \tfrac12 e^{-2t/3} - \tfrac12 e^{-2t}, \qquad i_2 = -\tfrac12 e^{-2t/3} + \tfrac12 e^{-2t}\ \ \mathrm{A}.
> $$
> **Interpretation.** The sum $i_1 + i_2$ sees inductance $L + M = 1.5$ H (time constant $1.5$ s, rate $2/3$), the difference $i_1 - i_2$ sees $L - M = 0.5$ H (rate $2$): these are the modal coordinates. The secondary current is negative (Lenz's law: it opposes the change in flux), peaks, and returns to zero. A numerical eigen-solution of the matrix form gives the same eigenvalues $-0.6667$ and $-2$ and the same coefficients.

> [!example] Example 12.27 — Coupled resonators and energy exchange
> Two identical LC resonators coupled through a capacitor (or two masses joined by springs) obey, in suitable units,
> $$
> x_1'' = -2x_1 + x_2, \qquad x_2'' = x_1 - 2x_2,
> $$
> with $x_1(0) = 1$, $x_2(0) = 0$, $x_1'(0) = x_2'(0) = 0$. Find the normal modes and the motion.
>
> **Normal modes.** $\mathbf{K} = \begin{bmatrix} 2 & -1 \\ -1 & 2\end{bmatrix}$ has eigenvalues $1$ and $3$ with eigenvectors $(1, 1)$ and $(1, -1)$. So the in-phase mode oscillates at $\omega_1 = 1$ and the anti-phase mode at $\omega_2 = \sqrt3 = 1.732$.
>
> **Modal coordinates.** $q_1 = (x_1 + x_2)/2$, $q_2 = (x_1 - x_2)/2$ satisfy $q_1'' = -q_1$, $q_2'' = -3q_2$, with $q_1(0) = q_2(0) = \tfrac12$ and zero initial velocities. Hence $q_1 = \tfrac12\cos t$, $q_2 = \tfrac12\cos\sqrt3 t$, and
> $$
> x_1 = q_1 + q_2 = \tfrac12\big(\cos t + \cos\sqrt3\,t\big), \qquad x_2 = q_1 - q_2 = \tfrac12\big(\cos t - \cos\sqrt3\,t\big).
> $$
> **Check.** $x_1'' = \tfrac12(-\cos t - 3\cos\sqrt3 t)$ and $-2x_1 + x_2 = -\cos t - \cos\sqrt3 t + \tfrac12\cos t - \tfrac12\cos\sqrt3 t = \tfrac12(-\cos t - 3\cos\sqrt3 t)$ ✓.
>
> By the sum-to-product identity, $x_1 = \cos\!\big(\tfrac{\sqrt3 - 1}{2}t\big)\cos\!\big(\tfrac{\sqrt3 + 1}{2}t\big)$ and $x_2 = \sin\!\big(\tfrac{\sqrt3 - 1}{2}t\big)\sin\!\big(\tfrac{\sqrt3 + 1}{2}t\big)$. The energy moves completely from resonator 1 to resonator 2 and back, with the transfer taking $\pi/(\sqrt3 - 1) = 4.29$ time units. With weak coupling the two mode frequencies are close and the exchange is slow; this is the mechanism behind coupled-resonator bandpass filters and wireless power transfer.

```desmos-graph
left=0; right=26; bottom=-1.15; top=1.15
---
y=0.5(\cos(x)+\cos(1.7321x))|blue|0<x<26
y=0.5(\cos(x)-\cos(1.7321x))|red|0<x<26
```
*Example 12.27: displacement of resonator 1 (blue) and resonator 2 (red). Energy passes back and forth between them as the two normal modes drift in and out of phase.*

> [!example] Example 12.28 — Driven coupled resonators and an antiresonance
> Drive the system of Example 12.27 with a force on the first coordinate: $x_1'' + 2x_1 - x_2 = \cos\omega t$, $x_2'' - x_1 + 2x_2 = 0$. Find the steady-state phasor amplitudes and show that $x_1$ can be made to stand still.
>
> Replace $D^2$ by $-\omega^2$ (strictly, the undamped system has no decaying transient; in practice a small loss makes the transient decay and leaves this steady state):
> $$
> \begin{bmatrix} 2 - \omega^2 & -1 \\ -1 & 2 - \omega^2\end{bmatrix}\begin{bmatrix} X_1 \\ X_2\end{bmatrix} = \begin{bmatrix} 1 \\ 0\end{bmatrix}
> \quad\Rightarrow\quad
> X_1 = \frac{2 - \omega^2}{(2 - \omega^2)^2 - 1}, \qquad X_2 = \frac{1}{(2 - \omega^2)^2 - 1}.
> $$
> The denominator $(2 - \omega^2)^2 - 1 = (1 - \omega^2)(3 - \omega^2)$ vanishes at the natural frequencies $\omega = 1$ and $\sqrt3$ (resonances). The numerator of $X_1$ vanishes at $\omega = \sqrt2$: there $X_1 = 0$ and $X_2 = 1/(0 - 1) = -1$. The driven coordinate does not move; the second resonator, tuned to $\omega = \sqrt2$ on its own ($x_2'' + 2x_2 = x_1$), oscillates in antiphase and exactly cancels the drive. This is the principle of the mechanical vibration absorber and of the parallel-resonant **wave trap** in power-line carrier systems: a resonator tuned to an unwanted frequency blocks it.

> [!tip] Engineering Note
> In nodal analysis of an RLC network the operator matrix $\mathbf{P}(s)$ is the admittance matrix $\mathbf{Y}(s)$, and $\Delta(s) = \det\mathbf{Y}(s)$, after clearing denominators, is the characteristic polynomial of the network. Circuit simulators find natural frequencies ("pole-zero analysis") by locating the zeros of this determinant, and they find transfer-function zeros from the cofactors, exactly the numerator determinants of Cramer's rule.

---

## 12.7 State-variable formulation

### Motivation

A single $n$th-order equation is ideal for hand analysis of small systems, but it scales badly. Eliminating variables from a ten-node network is painful; systems with several inputs and outputs do not fit a scalar equation at all; and numerical ODE solvers only accept first-order systems. The **state-variable** (state-space) formulation writes any linear lumped system as a single first-order *vector* equation. It handles any number of inputs and outputs, it is what simulators and control-design software use internally, and its solution has a closed form, the **matrix exponential**, which is the exact matrix analogue of $e^{at}$ for the scalar equation $\dot x = ax$.

### State, state equation, output equation

> [!definition] State
> The **state** of a system at time $t_0$ is the smallest set of numbers $x_1(t_0), \dots, x_n(t_0)$ which, together with the input $u(t)$ for $t \ge t_0$, determines the system's behaviour for all $t \ge t_0$. The vector $\mathbf{x} = (x_1, \dots, x_n)^\mathsf{T}$ is the **state vector**; $n$ is the **order**.

> [!definition] Linear time-invariant state-space model
> $$
> \begin{aligned}
> \dot{\mathbf{x}}(t) &= \mathbf{A}\,\mathbf{x}(t) + \mathbf{B}\,\mathbf{u}(t) &&\text{(state equation)}\\
> \mathbf{y}(t) &= \mathbf{C}\,\mathbf{x}(t) + \mathbf{D}\,\mathbf{u}(t) &&\text{(output equation)}
> \end{aligned}
> $$
> with $\mathbf{x} \in \mathbb{R}^n$ (state), $\mathbf{u} \in \mathbb{R}^p$ (inputs), $\mathbf{y} \in \mathbb{R}^q$ (outputs), and constant matrices $\mathbf{A}$ ($n\times n$, the **system** or **state matrix**), $\mathbf{B}$ ($n\times p$, input matrix), $\mathbf{C}$ ($q\times n$, output matrix), $\mathbf{D}$ ($q\times p$, feedthrough matrix).

The letter $\mathbf{D}$ for the feedthrough matrix is universal in control texts; context distinguishes it from the operator $D = d/dt$.

### Choosing state variables

**From a circuit.** The natural state variables are the **capacitor voltages** and **inductor currents**: they are continuous in time (finite current cannot change a capacitor voltage instantaneously, finite voltage cannot change an inductor current), they define the stored energy $\tfrac12 Cv_C^2 + \tfrac12 Li_L^2$, and the element laws already give their derivatives:

$$
\frac{dv_C}{dt} = \frac{i_C}{C}, \qquad \frac{di_L}{dt} = \frac{v_L}{L}.
$$

The recipe: (1) choose every independent $v_C$ and $i_L$ as a state; (2) using KCL, KVL and resistor laws, express each $i_C$ and $v_L$ in terms of the states and the sources; (3) divide by $C$ or $L$. The result is $\dot{\mathbf{x}} = \mathbf{A}\mathbf{x} + \mathbf{B}\mathbf{u}$ directly, with no elimination.

**From an $n$th-order ODE.** For $y^{(n)} + a_{n-1}y^{(n-1)} + \dots + a_1 y' + a_0 y = u$ choose the **phase variables** $x_1 = y$, $x_2 = y'$, ..., $x_n = y^{(n-1)}$. Then $\dot x_k = x_{k+1}$ for $k < n$ and $\dot x_n = -a_0x_1 - \dots - a_{n-1}x_n + u$:

$$
\mathbf{A} = \begin{bmatrix} 0 & 1 & 0 & \cdots & 0 \\ 0 & 0 & 1 & \cdots & 0 \\ \vdots & & & \ddots & \vdots \\ 0 & 0 & 0 & \cdots & 1 \\ -a_0 & -a_1 & -a_2 & \cdots & -a_{n-1}\end{bmatrix}, \quad
\mathbf{B} = \begin{bmatrix} 0 \\ 0 \\ \vdots \\ 0 \\ 1 \end{bmatrix}, \quad
\mathbf{C} = \begin{bmatrix} 1 & 0 & \cdots & 0 \end{bmatrix}, \quad \mathbf{D} = 0 .
$$

This is the **companion** (controllable canonical) form. Its characteristic polynomial $\det(s\mathbf{I} - \mathbf{A})$ equals $s^n + a_{n-1}s^{n-1} + \dots + a_0 = P(s)$, so the eigenvalues of $\mathbf{A}$ are exactly the characteristic roots of the ODE. When the right-hand side contains input derivatives, $P(D)y = (b_{n-1}D^{n-1} + \dots + b_0)u$, keep the same $\mathbf{A}$ and $\mathbf{B}$, solve $P(D)w = u$ for an auxiliary variable $w$ with states $w, w', \dots$, and take $y = b_0 w + b_1 w' + \dots$, i.e. $\mathbf{C} = [\,b_0\ \ b_1\ \cdots\ b_{n-1}\,]$.

> [!theorem] State representations are not unique
> If $\mathbf{T}$ is any invertible $n\times n$ matrix and $\mathbf{z} = \mathbf{T}^{-1}\mathbf{x}$, then
> $$
> \dot{\mathbf{z}} = \mathbf{T}^{-1}\mathbf{A}\mathbf{T}\,\mathbf{z} + \mathbf{T}^{-1}\mathbf{B}\,\mathbf{u}, \qquad \mathbf{y} = \mathbf{C}\mathbf{T}\,\mathbf{z} + \mathbf{D}\mathbf{u}
> $$
> describes the same input–output behaviour. Similar matrices have the same eigenvalues, so the natural frequencies, stability and transfer function are invariant; only the internal coordinates change. Choosing $\mathbf{T} = \mathbf{V}$ (eigenvectors) gives the **modal** (diagonal) form of Section 12.6.

### The matrix exponential

For the scalar equation $\dot x = ax$ the solution is $x = e^{at}x(0)$. To mimic this for $\dot{\mathbf{x}} = \mathbf{A}\mathbf{x}$, define the exponential of a matrix by the same power series.

> [!definition] Matrix exponential
> For a square matrix $\mathbf{A}$ and real $t$,
> $$
> e^{\mathbf{A}t} = \mathbf{I} + \mathbf{A}t + \frac{\mathbf{A}^2t^2}{2!} + \frac{\mathbf{A}^3t^3}{3!} + \dots = \sum_{k=0}^{\infty}\frac{\mathbf{A}^k t^k}{k!}.
> $$
> In control theory $\boldsymbol{\Phi}(t) = e^{\mathbf{A}t}$ is called the **state-transition matrix**.

The series converges absolutely for every $\mathbf{A}$ and $t$, because $\lVert\mathbf{A}^k t^k/k!\rVert \le (\lVert\mathbf{A}\rVert\,\lvert t\rvert)^k/k!$ and the scalar exponential series converges.

> [!theorem] Properties of $e^{\mathbf{A}t}$
> 1. $e^{\mathbf{A}\cdot 0} = \mathbf{I}$.
> 2. $\dfrac{d}{dt}e^{\mathbf{A}t} = \mathbf{A}e^{\mathbf{A}t} = e^{\mathbf{A}t}\mathbf{A}$.
> 3. $e^{\mathbf{A}(t + \tau)} = e^{\mathbf{A}t}e^{\mathbf{A}\tau}$ (semigroup property).
> 4. $\big(e^{\mathbf{A}t}\big)^{-1} = e^{-\mathbf{A}t}$; in particular $e^{\mathbf{A}t}$ is always invertible, with $\det e^{\mathbf{A}t} = e^{t\operatorname{tr}\mathbf{A}}$.
> 5. If $\mathbf{A} = \mathbf{V}\boldsymbol{\Lambda}\mathbf{V}^{-1}$, then $e^{\mathbf{A}t} = \mathbf{V}e^{\boldsymbol{\Lambda}t}\mathbf{V}^{-1}$ with $e^{\boldsymbol{\Lambda}t} = \operatorname{diag}(e^{\lambda_1 t}, \dots, e^{\lambda_n t})$.
> 6. $e^{(\mathbf{A} + \mathbf{B})t} = e^{\mathbf{A}t}e^{\mathbf{B}t}$ **if** $\mathbf{AB} = \mathbf{BA}$, but not in general.

*Proof sketches.* (2): differentiate the series term by term, $\frac{d}{dt}\frac{\mathbf{A}^k t^k}{k!} = \mathbf{A}\frac{\mathbf{A}^{k-1}t^{k-1}}{(k-1)!}$; $\mathbf{A}$ commutes with its own powers. (3): both sides satisfy $\dot{\mathbf{X}} = \mathbf{A}\mathbf{X}$ in $t$ with value $e^{\mathbf{A}\tau}$ at $t = 0$, so they are equal by uniqueness. (4) follows from (3) with $\tau = -t$. (5): $\mathbf{A}^k = \mathbf{V}\boldsymbol{\Lambda}^k\mathbf{V}^{-1}$, so each term of the series factors. (6): the binomial expansion of $(\mathbf{A} + \mathbf{B})^k$ needs commuting factors. $\square$

> [!warning] Common Pitfall
> $e^{\mathbf{A}t}$ is **not** the matrix of entry-wise exponentials $[e^{a_{ik}t}]$, except when $\mathbf{A}$ is diagonal. And $e^{\mathbf{A}t}e^{\mathbf{B}t} \ne e^{(\mathbf{A} + \mathbf{B})t}$ unless $\mathbf{A}$ and $\mathbf{B}$ commute; this matters, for example, in piecewise-linear switched circuits, where each switching interval has its own state matrix and the transition matrices must be multiplied in time order.

### Solving the state equation

> [!theorem] Solution of $\dot{\mathbf{x}} = \mathbf{A}\mathbf{x} + \mathbf{B}\mathbf{u}$
> The unique solution with $\mathbf{x}(0) = \mathbf{x}_0$ is
> $$
> \mathbf{x}(t) = \underbrace{e^{\mathbf{A}t}\mathbf{x}_0}_{\text{zero-input}} + \underbrace{\int_0^t e^{\mathbf{A}(t - \tau)}\mathbf{B}\,\mathbf{u}(\tau)\,d\tau}_{\text{zero-state}},
> $$
> and the output is $\mathbf{y} = \mathbf{C}\mathbf{x} + \mathbf{D}\mathbf{u}$.

*Derivation (integrating factor).* Multiply the state equation by $e^{-\mathbf{A}t}$ and use property 2:
$$
e^{-\mathbf{A}t}\big(\dot{\mathbf{x}} - \mathbf{A}\mathbf{x}\big) = \frac{d}{dt}\Big(e^{-\mathbf{A}t}\mathbf{x}\Big) = e^{-\mathbf{A}t}\mathbf{B}\mathbf{u}.
$$
Integrate from $0$ to $t$: $e^{-\mathbf{A}t}\mathbf{x}(t) - \mathbf{x}_0 = \int_0^t e^{-\mathbf{A}\tau}\mathbf{B}\mathbf{u}(\tau)\,d\tau$. Multiply by $e^{\mathbf{A}t}$ and use property 3. Uniqueness follows from the general theory ([Ch. 11](11-ordinary-differential-equations.md)). $\square$

The formula is the vector version of the first-order integrating-factor solution, and it displays the zero-input/zero-state decomposition of Section 12.3 explicitly. The matrix $e^{\mathbf{A}t}\mathbf{B}$ is the **impulse response** of the states, and $\mathbf{C}e^{\mathbf{A}t}\mathbf{B} + \mathbf{D}\delta(t)$ is the impulse response from input to output.

> [!theorem] Constant (step) input
> If $\mathbf{u}(t) = \mathbf{u}_0$ is constant for $t \ge 0$ and $\mathbf{A}$ is invertible,
> $$
> \mathbf{x}(t) = e^{\mathbf{A}t}\mathbf{x}_0 + \mathbf{A}^{-1}\big(e^{\mathbf{A}t} - \mathbf{I}\big)\mathbf{B}\mathbf{u}_0 .
> $$
> If moreover $\mathbf{A}$ is stable (all eigenvalues in the open left half-plane), $e^{\mathbf{A}t} \to \mathbf{0}$ and the steady state is $\mathbf{x}_{ss} = -\mathbf{A}^{-1}\mathbf{B}\mathbf{u}_0$, the solution of $\mathbf{0} = \mathbf{A}\mathbf{x} + \mathbf{B}\mathbf{u}_0$.

*Proof.* $\int_0^t e^{\mathbf{A}(t - \tau)}\,d\tau = \int_0^t e^{\mathbf{A}\sigma}\,d\sigma = \mathbf{A}^{-1}\big(e^{\mathbf{A}t} - \mathbf{I}\big)$, by integrating the series term by term (or because the derivative of the right side is $e^{\mathbf{A}t}$). $\square$

**Stability in state space.** The eigenvalues of $\mathbf{A}$ are the characteristic roots: $\det(s\mathbf{I} - \mathbf{A})$ is the characteristic polynomial. The system is asymptotically stable if and only if every eigenvalue of $\mathbf{A}$ has negative real part, which is the root criterion of Section 12.2 in matrix language.

**Frequency response in state space.** For $\mathbf{u} = \mathbf{U}e^{j\omega t}$, try $\mathbf{x} = \mathbf{X}e^{j\omega t}$: $j\omega\mathbf{X} = \mathbf{A}\mathbf{X} + \mathbf{B}\mathbf{U}$, so $\mathbf{X} = (j\omega\mathbf{I} - \mathbf{A})^{-1}\mathbf{B}\mathbf{U}$ and

$$
\mathbf{H}(j\omega) = \mathbf{C}\big(j\omega\mathbf{I} - \mathbf{A}\big)^{-1}\mathbf{B} + \mathbf{D}.
$$

Replacing $j\omega$ by $s$ gives the transfer-function matrix of [Ch. 14](../part-5-transforms/14-laplace-transforms.md), where it is also shown that $\mathcal{L}\{e^{\mathbf{A}t}\} = (s\mathbf{I} - \mathbf{A})^{-1}$, the **resolvent**. That identity provides a fourth way to compute $e^{\mathbf{A}t}$: invert $s\mathbf{I} - \mathbf{A}$ and take inverse Laplace transforms entry by entry.

### Computing $e^{\mathbf{A}t}$ by hand

Four methods cover practical needs.

**1. Diagonalisation** (property 5), when $\mathbf{A}$ has $n$ independent eigenvectors.

**2. Cayley–Hamilton / Sylvester interpolation.** Every square matrix satisfies its own characteristic equation, so $\mathbf{A}^n$ and all higher powers are combinations of $\mathbf{I}, \mathbf{A}, \dots, \mathbf{A}^{n-1}$. Hence
$$
e^{\mathbf{A}t} = \alpha_0(t)\mathbf{I} + \alpha_1(t)\mathbf{A} + \dots + \alpha_{n-1}(t)\mathbf{A}^{n-1},
$$
where the scalar functions $\alpha_k(t)$ are fixed by requiring $e^{\lambda t} = \alpha_0 + \alpha_1\lambda + \dots + \alpha_{n-1}\lambda^{n-1}$ at every eigenvalue $\lambda$ (and, for a repeated eigenvalue of multiplicity $m$, also the first $m - 1$ derivatives of this identity with respect to $\lambda$). For $2\times2$ matrices this gives closed forms:
$$
e^{\mathbf{A}t} = \frac{e^{\lambda_1 t}(\mathbf{A} - \lambda_2\mathbf{I}) - e^{\lambda_2 t}(\mathbf{A} - \lambda_1\mathbf{I})}{\lambda_1 - \lambda_2} \quad (\lambda_1 \ne \lambda_2), \qquad e^{\mathbf{A}t} = e^{\lambda t}\big[\mathbf{I} + (\mathbf{A} - \lambda\mathbf{I})t\big] \quad (\lambda_1 = \lambda_2 = \lambda).
$$
The repeated-eigenvalue formula follows because $(\mathbf{A} - \lambda\mathbf{I})^2 = \mathbf{0}$ by Cayley–Hamilton, so $e^{\mathbf{A}t} = e^{\lambda t}e^{(\mathbf{A} - \lambda\mathbf{I})t}$ and the second series terminates after two terms.

**3. Jordan / nilpotent splitting.** If $\mathbf{A} = \lambda\mathbf{I} + \mathbf{N}$ with $\mathbf{N}$ nilpotent ($\mathbf{N}^m = \mathbf{0}$), then since $\lambda\mathbf{I}$ commutes with everything, $e^{\mathbf{A}t} = e^{\lambda t}\big(\mathbf{I} + \mathbf{N}t + \dots + \mathbf{N}^{m-1}t^{m-1}/(m-1)!\big)$. This is the matrix origin of the $t^k e^{\lambda t}$ terms of Section 12.2.

**4. Laplace transform:** $e^{\mathbf{A}t} = \mathcal{L}^{-1}\{(s\mathbf{I} - \mathbf{A})^{-1}\}$ ([Ch. 14](../part-5-transforms/14-laplace-transforms.md)).

For complex eigenvalues $\sigma \pm j\omega$ of a real matrix, the useful real canonical block is
$$
\mathbf{A} = \begin{bmatrix} \sigma & \omega \\ -\omega & \sigma\end{bmatrix} = \sigma\mathbf{I} + \omega\mathbf{J}, \quad \mathbf{J} = \begin{bmatrix} 0 & 1 \\ -1 & 0\end{bmatrix}, \qquad e^{\mathbf{A}t} = e^{\sigma t}\begin{bmatrix} \cos\omega t & \sin\omega t \\ -\sin\omega t & \cos\omega t\end{bmatrix},
$$
because $\mathbf{J}^2 = -\mathbf{I}$ makes the series for $e^{\omega\mathbf{J}t}$ split exactly like Euler's formula. Geometrically: a rotation at angular rate $\omega$ combined with a uniform contraction at rate $-\sigma$.

### Worked examples

> [!example] Example 12.29 — Companion form of a third-order equation
> Write $y''' + 6y'' + 11y' + 6y = u$ in state-space form and confirm the eigenvalues.
>
> With $x_1 = y$, $x_2 = y'$, $x_3 = y''$:
> $$
> \dot{\mathbf{x}} = \begin{bmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ -6 & -11 & -6\end{bmatrix}\mathbf{x} + \begin{bmatrix} 0 \\ 0 \\ 1\end{bmatrix}u, \qquad y = \begin{bmatrix} 1 & 0 & 0\end{bmatrix}\mathbf{x}.
> $$
> Expanding $\det(s\mathbf{I} - \mathbf{A})$ along the first column: $s\big[s(s + 6) + 11\big] + 6 = s^3 + 6s^2 + 11s + 6 = (s + 1)(s + 2)(s + 3)$. Eigenvalues $-1, -2, -3$: stable, overdamped, slowest time constant $1$ s. If the equation were instead $y''' + 6y'' + 11y' + 6y = 2u' + 5u$, the same $\mathbf{A}$ and $\mathbf{B}$ would be kept with $\mathbf{C} = [\,5\ \ 2\ \ 0\,]$.

> [!example] Example 12.30 — $e^{\mathbf{A}t}$ three ways
> Compute $e^{\mathbf{A}t}$ for $\mathbf{A} = \begin{bmatrix} 0 & 1 \\ -2 & -3\end{bmatrix}$ (the companion matrix of $y'' + 3y' + 2y$).
>
> **Eigenvalues.** $\det(s\mathbf{I} - \mathbf{A}) = s^2 + 3s + 2$: $\lambda_1 = -1$, $\lambda_2 = -2$.
>
> **Sylvester formula.** $\mathbf{A} - \lambda_2\mathbf{I} = \mathbf{A} + 2\mathbf{I} = \begin{bmatrix} 2 & 1 \\ -2 & -1\end{bmatrix}$, $\mathbf{A} - \lambda_1\mathbf{I} = \mathbf{A} + \mathbf{I} = \begin{bmatrix} 1 & 1 \\ -2 & -2\end{bmatrix}$, $\lambda_1 - \lambda_2 = 1$:
> $$
> e^{\mathbf{A}t} = e^{-t}\begin{bmatrix} 2 & 1 \\ -2 & -1\end{bmatrix} - e^{-2t}\begin{bmatrix} 1 & 1 \\ -2 & -2\end{bmatrix} = \begin{bmatrix} 2e^{-t} - e^{-2t} & e^{-t} - e^{-2t} \\ -2e^{-t} + 2e^{-2t} & -e^{-t} + 2e^{-2t}\end{bmatrix}.
> $$
> **Diagonalisation check.** Eigenvectors $\mathbf{v}_1 = (1, -1)$, $\mathbf{v}_2 = (1, -2)$; $\mathbf{V} = \begin{bmatrix} 1 & 1 \\ -1 & -2\end{bmatrix}$, $\mathbf{V}^{-1} = \begin{bmatrix} 2 & 1 \\ -1 & -1\end{bmatrix}$. Then $\mathbf{V}\operatorname{diag}(e^{-t}, e^{-2t})\mathbf{V}^{-1}$ reproduces the same matrix.
>
> **Sanity checks.** At $t = 0$ the matrix is $\mathbf{I}$ ✓. The derivative at $t = 0$ is $\begin{bmatrix} -2 + 2 & -1 + 2 \\ 2 - 4 & 1 - 4\end{bmatrix} = \mathbf{A}$ ✓. Entry $(1,1)$ is the zero-input $y$ for $y(0) = 1$, $y'(0) = 0$, and $(1,2)$ is the impulse response $e^{-t} - e^{-2t}$ of Example 12.13 ✓. At $t = 0.7$ the entries are $0.74657$, $0.24999$, $-0.49998$, $-0.00339$, identical to a numerical `expm` to all printed digits.

> [!example] Example 12.31 — Repeated and complex eigenvalues
> (a) $\mathbf{A} = \begin{bmatrix} -2 & 1 \\ 0 & -2\end{bmatrix}$; (b) $\mathbf{A} = \begin{bmatrix} -1 & 2 \\ -2 & -1\end{bmatrix}$.
>
> **(a)** Double eigenvalue $-2$ with only one eigenvector, so $\mathbf{A}$ is not diagonalisable. Write $\mathbf{A} = -2\mathbf{I} + \mathbf{N}$ with $\mathbf{N} = \begin{bmatrix} 0 & 1 \\ 0 & 0\end{bmatrix}$, $\mathbf{N}^2 = \mathbf{0}$:
> $$
> e^{\mathbf{A}t} = e^{-2t}(\mathbf{I} + \mathbf{N}t) = e^{-2t}\begin{bmatrix} 1 & t \\ 0 & 1\end{bmatrix}.
> $$
> The $te^{-2t}$ entry is the matrix counterpart of a critically damped response.
>
> **(b)** This is the real canonical block with $\sigma = -1$, $\omega = 2$ (eigenvalues $-1 \pm 2j$):
> $$
> e^{\mathbf{A}t} = e^{-t}\begin{bmatrix} \cos 2t & \sin 2t \\ -\sin 2t & \cos 2t\end{bmatrix}.
> $$
> At $t = 0.5$: $e^{-0.5}\cos 1 = 0.32771$, $e^{-0.5}\sin 1 = 0.51038$, matching numerical computation. Trajectories spiral into the origin, the phase-plane picture of an underdamped system ([Ch. 11](11-ordinary-differential-equations.md)).

> [!example] Example 12.32 — Forced solution by the matrix formula
> For $\dot{\mathbf{x}} = \mathbf{A}\mathbf{x} + \mathbf{B}u$ with $\mathbf{A}$ from Example 12.30, $\mathbf{B} = (0, 1)^\mathsf{T}$, $u = 1$ for $t \ge 0$ and $\mathbf{x}_0 = \mathbf{0}$, find $\mathbf{x}(t)$.
>
> $\mathbf{A}^{-1} = \tfrac12\begin{bmatrix} -3 & -1 \\ 2 & 0\end{bmatrix}$, so $\mathbf{A}^{-1}\mathbf{B} = (-\tfrac12, 0)^\mathsf{T}$. Since $\mathbf{A}^{-1}$ commutes with $e^{\mathbf{A}t}$:
> $$
> \mathbf{x}(t) = e^{\mathbf{A}t}\mathbf{A}^{-1}\mathbf{B} - \mathbf{A}^{-1}\mathbf{B} = -\tfrac12\begin{bmatrix} 2e^{-t} - e^{-2t} \\ -2e^{-t} + 2e^{-2t}\end{bmatrix} + \begin{bmatrix} \tfrac12 \\ 0\end{bmatrix} = \begin{bmatrix} \tfrac12 - e^{-t} + \tfrac12 e^{-2t} \\ e^{-t} - e^{-2t}\end{bmatrix}.
> $$
> The first state is the unit-step response of $y'' + 3y' + 2y = u$ found by scalar methods in Example 12.13 ✓, and the second state is its derivative ✓. The steady state $-\mathbf{A}^{-1}\mathbf{B} = (\tfrac12, 0)$ is the DC gain $1/P(0) = 1/2$ ✓.

> [!example] Example 12.33 — Exact discretisation for simulation
> A digital controller samples the system of Example 12.32 every $T = 0.1$ s and holds $u$ constant between samples (zero-order hold). Derive the exact discrete-time model.
>
> Over one interval $[kT, (k + 1)T]$ with $u = u[k]$ constant, the step-input formula with $\mathbf{x}_0 = \mathbf{x}[k]$ gives
> $$
> \mathbf{x}[k + 1] = \boldsymbol{\Phi}\,\mathbf{x}[k] + \boldsymbol{\Gamma}\,u[k], \qquad \boldsymbol{\Phi} = e^{\mathbf{A}T}, \quad \boldsymbol{\Gamma} = \mathbf{A}^{-1}\big(e^{\mathbf{A}T} - \mathbf{I}\big)\mathbf{B}.
> $$
> Numerically, from Example 12.30 with $t = 0.1$:
> $$
> \boldsymbol{\Phi} = \begin{bmatrix} 0.99094 & 0.08611 \\ -0.17221 & 0.73262\end{bmatrix}, \qquad \boldsymbol{\Gamma} = \begin{bmatrix} 0.004528 \\ 0.086107\end{bmatrix}.
> $$
> Iterating ten steps from rest with $u[k] = 1$ gives $x_1[10] = 0.199788$, equal to the exact continuous value $\tfrac12 - e^{-1} + \tfrac12e^{-2} = 0.199788$. The discretisation is exact at the sample instants, not an approximation like Euler's method; this is how simulators and digital-control design tools convert continuous models into difference equations ([Ch. 17](../part-5-transforms/17-z-transforms.md)). The eigenvalues of $\boldsymbol{\Phi}$ are $e^{\lambda_k T} = e^{-0.1}, e^{-0.2}$: continuous poles map to discrete poles by $z = e^{sT}$.

> [!tip] Engineering Note
> The state-space model is the native language of modern tools. SPICE assembles nodal equations that are a (descriptor) state-space system; MATLAB/Python control libraries store models as $(\mathbf{A}, \mathbf{B}, \mathbf{C}, \mathbf{D})$; Kalman filters, observers, optimal (LQR) controllers and model-order reduction all operate on these matrices. Numerically, $e^{\mathbf{A}T}$ is computed by scaling and squaring with a Padé approximant (the `expm` routine), not by summing the raw series, which suffers from cancellation when $\lVert\mathbf{A}T\rVert$ is large.

---

## 12.8 Applications to RLC circuits and dynamic systems

This section applies the whole chapter to the circuits and systems that motivated it. The methods are those of Sections 12.1–12.7; the new ingredients are the circuit-specific steps of **obtaining initial conditions** from the physics of switching and **interpreting** the parameters $\alpha$, $\omega_n$, $\zeta$ and $Q$ in terms of component values.

### Initial conditions at a switching instant

A switch operating at $t = 0$ changes the circuit, so the equation valid for $t > 0$ needs initial values at $t = 0^+$. They come from two continuity principles, both consequences of finite power:

> [!theorem] Continuity of state variables
> In a circuit with finite voltages and currents (no impulsive sources and no loops of capacitors or cut-sets of inductors forced to change state instantaneously):
> $$
> v_C(0^+) = v_C(0^-), \qquad i_L(0^+) = i_L(0^-).
> $$

*Reason.* $v_C(0^+) - v_C(0^-) = \frac1C\int_{0^-}^{0^+} i_C\,dt = 0$ when $i_C$ is bounded; similarly for $i_L$ with $v_L$. Equivalently, stored energy $\tfrac12Cv_C^2$ or $\tfrac12Li_L^2$ cannot jump without infinite power. $\square$

All other quantities (resistor currents, inductor voltages, capacitor currents) **may** jump. The derivative initial conditions needed for a second-order equation are found from these jump-allowed quantities at $t = 0^+$:

$$
\frac{dv_C}{dt}(0^+) = \frac{i_C(0^+)}{C}, \qquad \frac{di_L}{dt}(0^+) = \frac{v_L(0^+)}{L},
$$

where $i_C(0^+)$ and $v_L(0^+)$ are computed by KCL/KVL on the $t = 0^+$ circuit with the known $v_C(0^+)$ and $i_L(0^+)$. In practice, the values before switching come from the DC steady state of the old circuit, in which capacitors are open circuits and inductors are short circuits.

The complete procedure for a switched linear circuit is therefore:

1. Find $v_C(0^-)$, $i_L(0^-)$ from the pre-switching circuit (usually its DC steady state).
2. Carry them across: $v_C(0^+) = v_C(0^-)$, $i_L(0^+) = i_L(0^-)$.
3. Find the derivative initial values from the $t = 0^+$ circuit.
4. Find the final (forced) value from the post-switching steady state.
5. Find the characteristic roots of the post-switching circuit.
6. Assemble $y = y_f + y_n$ and fit the constants.

### The series RLC circuit

For the series loop (source $v_s$, all elements carrying the same current $i$), with $v_C$ as unknown:

$$
v_C'' + \frac{R}{L}v_C' + \frac{1}{LC}v_C = \frac{1}{LC}v_s .
$$

Comparing with the standard form $y'' + 2\alpha y' + \omega_n^2 y$:

$$
\alpha = \frac{R}{2L}, \qquad \omega_n = \frac{1}{\sqrt{LC}}, \qquad \zeta = \frac{\alpha}{\omega_n} = \frac{R}{2}\sqrt{\frac{C}{L}}, \qquad Q = \frac{1}{2\zeta} = \frac{1}{R}\sqrt{\frac{L}{C}} .
$$

The **critical resistance** $R_{\text{crit}} = 2\sqrt{L/C}$ separates underdamped ($R < R_{\text{crit}}$) from overdamped ($R > R_{\text{crit}}$) behaviour. The quantity $Z_0 = \sqrt{L/C}$ is the **characteristic impedance** of the LC pair; damping depends on how the series resistance compares with it. In the series circuit, *more* resistance means *more* damping.

> [!theorem] Underdamped step response of the standard second-order system
> For $y'' + 2\zeta\omega_n y' + \omega_n^2 y = \omega_n^2 V u(t)$ from rest, with $0 < \zeta < 1$, $\alpha = \zeta\omega_n$, $\omega_d = \omega_n\sqrt{1 - \zeta^2}$:
> $$
> y(t) = V\left[1 - e^{-\alpha t}\left(\cos\omega_d t + \frac{\alpha}{\omega_d}\sin\omega_d t\right)\right] = V\left[1 - \frac{e^{-\alpha t}}{\sqrt{1 - \zeta^2}}\cos(\omega_d t - \varphi)\right],\quad \sin\varphi = \zeta .
> $$
> Its derivative is $y'(t) = V\dfrac{\omega_n^2}{\omega_d}e^{-\alpha t}\sin\omega_d t$. The response has
> - **peak time** $t_p = \pi/\omega_d$,
> - **maximum (percent) overshoot** $M_p = e^{-\pi\zeta/\sqrt{1 - \zeta^2}}$,
> - **2% settling time** $t_s \approx 4/\alpha$.

*Derivation.* The forced value is $V$. With $y = V + e^{-\alpha t}(A\cos\omega_d t + B\sin\omega_d t)$, $y(0) = 0$ gives $A = -V$, and $y'(0) = -\alpha A + \omega_d B = 0$ gives $B = \alpha A/\omega_d = -V\alpha/\omega_d$. Differentiating, the cosine terms cancel and $y' = V\big(\alpha^2/\omega_d + \omega_d\big)e^{-\alpha t}\sin\omega_d t = V(\omega_n^2/\omega_d)e^{-\alpha t}\sin\omega_d t$. This vanishes first (for $t > 0$) at $\omega_d t = \pi$, the peak. There $\cos\pi = -1$, $\sin\pi = 0$, so $y(t_p) = V(1 + e^{-\alpha\pi/\omega_d})$, and $\alpha/\omega_d = \zeta/\sqrt{1 - \zeta^2}$. $\square$

The overshoot depends **only on $\zeta$**: $16.3\%$ at $\zeta = 0.5$, $4.3\%$ at $\zeta = 0.707$, $52.7\%$ at $\zeta = 0.2$, and zero for $\zeta \ge 1$. Inverting, a measured overshoot gives the damping ratio: $\zeta = -\ln M_p/\sqrt{\pi^2 + \ln^2 M_p}$.

```desmos-graph
left=-0.3; right=14; bottom=-0.05; top=1.65
---
y=1-e^{-0.2x}(\cos(0.9798x)+0.2041\sin(0.9798x))|red|0<x<14
y=1-e^{-0.5x}(\cos(0.866x)+0.5774\sin(0.866x))|orange|0<x<14
y=1-e^{-0.7071x}(\cos(0.7071x)+\sin(0.7071x))|blue|0<x<14
y=1-(1+x)e^{-x}|green|0<x<14
y=1-1.0774e^{-0.2679x}+0.0774e^{-3.7321x}|purple|0<x<14
y=1|grey|dashed
y=1.02|grey|dotted
y=0.98|grey|dotted
```
*Normalised step responses against $\omega_n t$ for $\zeta = 0.2$ (red), $0.5$ (orange), $0.707$ (blue), $1$ (green) and $2$ (purple), with the $\pm2\%$ band dotted. Lower damping rises faster but overshoots and rings; heavy damping is sluggish.*

### The parallel RLC circuit

For $R$, $L$, $C$ in parallel driven by a current source $i_s$, take the common voltage $v$ (or the inductor current) as unknown. KCL, $v/R + i_L + Cv' = i_s$ with $v = Li_L'$, gives

$$
i_L'' + \frac{1}{RC}i_L' + \frac{1}{LC}i_L = \frac{1}{LC}i_s, \qquad \alpha = \frac{1}{2RC}, \qquad \zeta = \frac{1}{2R}\sqrt{\frac{L}{C}}, \qquad Q = R\sqrt{\frac{C}{L}} .
$$

The parallel circuit is the **dual** of the series one: here a *larger* resistance means *less* damping, since a large parallel resistor draws little current and dissipates little energy. An ideal tank with $R \to \infty$ is undamped.

### Worked examples

> [!example] Example 12.34 — Series RLC step response with numbers
> $R = 100\ \Omega$, $L = 10\ \mathrm{mH}$, $C = 1\ \mu\mathrm{F}$; a $10\ \mathrm{V}$ step is applied at $t = 0$ with the circuit at rest. Find $v_C(t)$, $i(t)$, the peak voltage and its time, and the settling time.
>
> **Parameters.** $\omega_n = 1/\sqrt{10^{-8}} = 10^4\ \mathrm{rad/s}$; $\alpha = R/2L = 5000\ \mathrm{s^{-1}}$; $\zeta = 0.5$; $\omega_d = 10^4\sqrt{0.75} = 8660\ \mathrm{rad/s}$; $Q = 1$; $R_{\text{crit}} = 2\sqrt{L/C} = 200\ \Omega$, so the circuit is underdamped.
>
> **Initial conditions.** $v_C(0^+) = 0$ and $i(0^+) = 0$ (continuity), so $v_C'(0^+) = i/C = 0$.
>
> **Response.**
> $$
> v_C(t) = 10\left[1 - e^{-5000t}\big(\cos 8660t + 0.5774\sin 8660t\big)\right]\ \mathrm{V},
> $$
> $$
> i(t) = C\,v_C' = C\cdot 10\cdot\frac{\omega_n^2}{\omega_d}e^{-5000t}\sin 8660t = \frac{10}{L\omega_d}e^{-5000t}\sin 8660t = 0.1155\,e^{-5000t}\sin 8660t\ \mathrm{A}.
> $$
> **Specifications.** $t_p = \pi/\omega_d = 362.8\ \mu\mathrm{s}$; $M_p = e^{-\pi\cdot 0.5/0.866} = 0.163$, so $v_{C,\max} = 11.63\ \mathrm{V}$; $t_s \approx 4/\alpha = 0.8\ \mathrm{ms}$.
>
> **Checks.** Dimension of $10/(L\omega_d)$: $\mathrm{V}/(\mathrm{H}\cdot\mathrm{s^{-1}}) = \mathrm{V}/\Omega = \mathrm{A}$ ✓. The capacitor voltage exceeds the supply: the inductor keeps current flowing after $v_C$ reaches $10$ V. A designer who needs no overshoot would raise $R$ to at least $200\ \Omega$.

> [!example] Example 12.35 — Source-free parallel RLC
> A parallel RLC circuit with $R = 250\ \Omega$, $L = 10\ \mathrm{mH}$, $C = 1\ \mu\mathrm{F}$ has $v(0) = 10\ \mathrm{V}$ on the capacitor and $i_L(0) = 0$. Find $v(t)$.
>
> **Parameters.** $\omega_n = 10^4\ \mathrm{rad/s}$, $\alpha = 1/(2RC) = 2000\ \mathrm{s^{-1}}$, $\zeta = 0.2$, $\omega_d = \sqrt{10^8 - 4\times10^6} = 9798\ \mathrm{rad/s}$.
>
> **Derivative initial condition.** KCL at $t = 0^+$: $i_C = -v/R - i_L = -10/250 - 0 = -0.04\ \mathrm{A}$, so $v'(0^+) = i_C/C = -4\times10^4\ \mathrm{V/s}$.
>
> **Fit.** $v = e^{-2000t}(A\cos\omega_d t + B\sin\omega_d t)$ with $A = 10$ and $-2000A + 9798B = -40000$, giving $B = -2.041$:
> $$
> v(t) = e^{-2000t}\big(10\cos 9798t - 2.041\sin 9798t\big)\ \mathrm{V}.
> $$
> **Check.** At $t = 100\ \mu\mathrm{s}$ the formula gives $3.1741\ \mathrm{V}$; a numerical integration of the state equations $v' = -(v/R + i_L)/C$, $i_L' = v/L$ gives $3.1741\ \mathrm{V}$ ✓. With $Q = R\sqrt{C/L} = 2.5$, the ringing decays by a factor $e^{-2\pi\zeta/\sqrt{1 - \zeta^2}} = 0.277$ per cycle.

> [!example] Example 12.36 — Switching a source: critically damped transition
> A series RLC circuit with $R = 4\ \Omega$, $L = 1\ \mathrm{H}$, $C = 0.25\ \mathrm{F}$ has been driven by $5\ \mathrm{V}$ DC for a long time. At $t = 0$ the source steps to $15\ \mathrm{V}$. Find $v_C(t)$.
>
> **Before switching** (DC steady state): $i(0^-) = 0$, $v_C(0^-) = 5\ \mathrm{V}$.
>
> **At $0^+$.** $v_C(0^+) = 5$, $i(0^+) = 0$, so $v_C'(0^+) = i/C = 0$.
>
> **After switching.** $\omega_n = 1/\sqrt{0.25} = 2\ \mathrm{rad/s}$, $\zeta = \tfrac{R}{2}\sqrt{C/L} = 2\times0.5 = 1$: critically damped, double root $-2$. Forced value $15\ \mathrm{V}$.
>
> **Fit.** $v_C = 15 + (A + Bt)e^{-2t}$: $15 + A = 5$ gives $A = -10$; $v_C'(0) = B - 2A = 0$ gives $B = -20$:
> $$
> v_C(t) = 15 - 10(1 + 2t)e^{-2t}\ \mathrm{V}, \qquad i(t) = Cv_C' = 0.25\cdot 40\,t\,e^{-2t} = 10\,t\,e^{-2t}\ \mathrm{A}.
> $$
> **Check.** $v_C' = -20e^{-2t} + 20(1 + 2t)e^{-2t} = 40te^{-2t}$ ✓. The current peaks at $t = 0.5\ \mathrm{s}$ with $i = 5e^{-1} = 1.84\ \mathrm{A}$, and $v_C$ rises monotonically from 5 V to 15 V without overshoot, as critical damping guarantees. Only the *change* of input ($10$ V) appears in the transient: the zero-input part (initial 5 V) and the zero-state part combine into a single natural response.

> [!example] Example 12.37 — DC motor: an electromechanical coupled system
> A permanent-magnet DC motor has armature resistance $R = 2\ \Omega$, inductance $L = 0.5\ \mathrm{H}$, torque/back-EMF constant $K = 0.1\ \mathrm{N\,m/A} = 0.1\ \mathrm{V\,s/rad}$, rotor inertia $J = 0.02\ \mathrm{kg\,m^2}$ and viscous friction $b = 0.01\ \mathrm{N\,m\,s}$. A $12\ \mathrm{V}$ step is applied at rest. Find the current and speed.
>
> **Coupled equations.** Electrical (KVL with back-EMF $K\omega$) and mechanical (Newton with motor torque $Ki$):
> $$
> L\,i' + R\,i + K\omega = v, \qquad J\,\omega' + b\,\omega = K\,i .
> $$
> **State-space form.** With $\mathbf{x} = (i, \omega)^\mathsf{T}$:
> $$
> \dot{\mathbf{x}} = \begin{bmatrix} -R/L & -K/L \\ K/J & -b/J\end{bmatrix}\mathbf{x} + \begin{bmatrix} 1/L \\ 0\end{bmatrix}v = \begin{bmatrix} -4 & -0.2 \\ 5 & -0.5\end{bmatrix}\mathbf{x} + \begin{bmatrix} 2 \\ 0\end{bmatrix}v .
> $$
> **Eigenvalues.** $\operatorname{tr}\mathbf{A} = -4.5$, $\det\mathbf{A} = 2 + 1 = 3$: $\lambda^2 + 4.5\lambda + 3 = 0$, $\lambda = (-4.5 \pm\sqrt{8.25})/2 = -0.8139,\ -3.6861\ \mathrm{s^{-1}}$. Both real and negative: stable, overdamped. The slow root is mostly mechanical (inertia), the fast one mostly electrical.
>
> **Steady state.** $\mathbf{x}_{ss} = -\mathbf{A}^{-1}\mathbf{B}\,v$: from $\omega' = 0$, $i = b\omega/K = 0.1\omega$; from $i' = 0$, $2i + 0.1\omega = 12$, so $0.3\omega = 12$, $\omega_{ss} = 40\ \mathrm{rad/s}$, $i_{ss} = 4\ \mathrm{A}$.
>
> **Complete solution** (eigenvector expansion of $\mathbf{x} - \mathbf{x}_{ss}$, fitted to $\mathbf{x}(0) = \mathbf{0}$):
> $$
> \begin{aligned}
> i(t) &= 4 - 7.222\,e^{-3.686t} + 3.222\,e^{-0.814t}\ \mathrm{A},\\
> \omega(t) &= 40 + 11.334\,e^{-3.686t} - 51.334\,e^{-0.814t}\ \mathrm{rad/s}.
> \end{aligned}
> $$
> **Checks.** $i(0) = 4 - 7.222 + 3.222 = 0$ ✓, $\omega(0) = 0$ ✓; $i'(0) = 7.222(3.686) - 3.222(0.814) = 24.0 = 12/L$ ✓ (initially the full voltage appears across $L$); $\omega'(0) = -11.334(3.686) + 51.334(0.814) = 0.0$ ✓ (no current, no torque at $t = 0$). The starting current overshoots to $5.30\ \mathrm{A}$ at $t = 0.81\ \mathrm{s}$, before back-EMF builds up; the speed settles within 2% after about $5.1\ \mathrm{s}$, governed by the slow root ($4/0.814 = 4.9\ \mathrm{s}$).

> [!example] Example 12.38 — RLC in state-variable form with a sinusoidal drive
> Write the series RLC circuit of Example 12.34 in state form with states $v_C$ and $i_L$, and find the steady-state capacitor voltage for $v_s = 10\cos(10^4 t)$ V.
>
> **State equations.** $Cv_C' = i_L$ and $Li_L' = v_s - Ri_L - v_C$:
> $$
> \frac{d}{dt}\begin{bmatrix} v_C \\ i_L\end{bmatrix} = \begin{bmatrix} 0 & 1/C \\ -1/L & -R/L\end{bmatrix}\begin{bmatrix} v_C \\ i_L\end{bmatrix} + \begin{bmatrix} 0 \\ 1/L\end{bmatrix}v_s, \qquad \det(s\mathbf{I} - \mathbf{A}) = s^2 + \frac{R}{L}s + \frac{1}{LC}.
> $$
> **Frequency response at $\omega = \omega_n$.** $H(j\omega) = \mathbf{C}(j\omega\mathbf{I} - \mathbf{A})^{-1}\mathbf{B}$ with $\mathbf{C} = [1\ \ 0]$ reduces to $\omega_n^2/(\omega_n^2 - \omega^2 + j2\zeta\omega_n\omega)$. At $\omega = \omega_n$ this is $1/(j2\zeta) = -jQ = -j$. So
> $$
> v_{C,ss}(t) = 10\cos(10^4 t - 90^\circ) = 10\sin(10^4 t)\ \mathrm{V}.
> $$
> At resonance the capacitor voltage is $Q$ times the source amplitude and lags it by exactly $90^\circ$. For a high-$Q$ circuit ($Q = 50$, say) the capacitor would see $50\times$ the source voltage, which is why resonant circuits need components rated well above the drive level.

### Analogies and generality

The same mathematics governs mechanical, thermal and hydraulic systems through the standard analogies:

| Electrical (series) | Mechanical (translation) | Mechanical (rotation) |
|---|---|---|
| voltage $v$ | force $F$ | torque $T$ |
| current $i$ | velocity $\dot x$ | angular velocity $\omega$ |
| inductance $L$ | mass $m$ | inertia $J$ |
| resistance $R$ | damping $b$ | rotational damping $b$ |
| elastance $1/C$ | stiffness $k$ | torsional stiffness $k$ |
| $\zeta = \tfrac{R}{2}\sqrt{C/L}$ | $\zeta = \tfrac{b}{2\sqrt{km}}$ | $\zeta = \tfrac{b}{2\sqrt{kJ}}$ |

Any lumped linear system, of whatever physical type, reduces to $\dot{\mathbf{x}} = \mathbf{A}\mathbf{x} + \mathbf{B}\mathbf{u}$, and its behaviour is read from the eigenvalues of $\mathbf{A}$. Distributed systems such as transmission lines and heat conduction lead to partial differential equations ([Ch. 13](13-partial-differential-equations.md)); after separation of variables, each spatial mode again obeys a constant-coefficient ODE of exactly the type studied here.

> [!warning] Common Pitfall
> Three errors dominate switched-circuit problems: (1) assuming $i_C$ or $v_L$ is continuous (they are not; only $v_C$ and $i_L$ are); (2) using the *pre*-switching circuit's $\alpha$ and $\omega_n$ for the *post*-switching transient; (3) setting $v_C'(0^+) = 0$ by habit, when it equals $i_C(0^+)/C$, which is often nonzero (Example 12.35).

> [!tip] Engineering Note
> The parameters of this section are design handles. A switching converter's output LC filter is designed for a chosen $\omega_n$ (well below the switching frequency) and a damping ratio set by the load resistance and capacitor ESR; a light load (large parallel $R$) gives low damping and visible ringing on load steps. A gate-drive loop is a series RLC (trace inductance, gate capacitance, gate resistor) and the gate resistor is chosen near $R_{\text{crit}} = 2\sqrt{L/C}$ to suppress ringing without slowing the edge too much. A PLL's loop filter is chosen to give $\zeta \approx 0.7$ for the closed loop.

---

## Summary

| Topic | Key formula / fact |
|---|---|
| Operator form | $P(D)y = r$, $P(s) = a_ns^n + \dots + a_0$; $P(D)e^{st} = P(s)e^{st}$ |
| Exponential shift | $P(D)[e^{at}f] = e^{at}P(D + a)f$ |
| Distinct roots | $y_h = \sum c_k e^{s_k t}$ |
| Repeated root $s_1$ (mult. $m$) | $(c_1 + c_2t + \dots + c_mt^{m-1})e^{s_1t}$ |
| Complex pair $\sigma \pm j\omega$ | $e^{\sigma t}(A\cos\omega t + B\sin\omega t)$ |
| Stability | all roots $\operatorname{Re}\{s_k\} < 0$; cubic $s^3 + a_2s^2 + a_1s + a_0$: $a_i > 0$, $a_2a_1 > a_0$ |
| Standard 2nd order | $y'' + 2\zeta\omega_ny' + \omega_n^2y$; $s = -\zeta\omega_n \pm \omega_n\sqrt{\zeta^2 - 1}$ |
| Damping cases | $\zeta > 1$ over, $\zeta = 1$ critical, $0 < \zeta < 1$ under ($\omega_d = \omega_n\sqrt{1 - \zeta^2}$), $\zeta = 0$ undamped |
| Complete response | $y = y_n + y_f = y_{zi} + y_{zs}$; fit ICs to the complete $y$ |
| Impulse response | homogeneous with $h^{(n-1)}(0^+) = 1/a_n$; $y_{zs} = \int_0^t h(t - \tau)r(\tau)\,d\tau$ |
| Exponential input | $y_p = e^{at}/P(a)$; resonant: $y_p = t^me^{at}/P^{(m)}(a)$ |
| DC gain | $y_p = K/P(0) = K/a_0$ |
| Sinusoidal steady state | $\mathbf{Y} = H(j\omega)\mathbf{R}$, $H = N(j\omega)/P(j\omega)$ |
| Resonant peak | $\omega_r = \omega_n\sqrt{1 - 2\zeta^2}$, $M_r = 1/(2\zeta\sqrt{1 - \zeta^2})$, $Q = 1/(2\zeta)$ |
| Settling time (2%) | $t_s \approx 4\tau$ (1st order), $t_s \approx 4/(\zeta\omega_n)$ (2nd order) |
| Step overshoot | $M_p = e^{-\pi\zeta/\sqrt{1 - \zeta^2}}$, $t_p = \pi/\omega_d$ |
| Coupled systems | $\Delta(s) = \det\mathbf{P}(s)$; $\Delta(D)y_i = \det\mathbf{P}_i(D)$; $\deg\Delta$ constants |
| State space | $\dot{\mathbf{x}} = \mathbf{A}\mathbf{x} + \mathbf{B}\mathbf{u}$, $\mathbf{y} = \mathbf{C}\mathbf{x} + \mathbf{D}\mathbf{u}$ |
| Matrix exponential | $e^{\mathbf{A}t} = \sum\mathbf{A}^kt^k/k! = \mathbf{V}e^{\boldsymbol{\Lambda}t}\mathbf{V}^{-1}$ |
| State solution | $\mathbf{x} = e^{\mathbf{A}t}\mathbf{x}_0 + \int_0^te^{\mathbf{A}(t - \tau)}\mathbf{B}\mathbf{u}(\tau)\,d\tau$ |
| Step input | $\mathbf{x} = e^{\mathbf{A}t}\mathbf{x}_0 + \mathbf{A}^{-1}(e^{\mathbf{A}t} - \mathbf{I})\mathbf{B}\mathbf{u}_0$; $\mathbf{x}_{ss} = -\mathbf{A}^{-1}\mathbf{B}\mathbf{u}_0$ |
| $2\times2$ Sylvester | $e^{\mathbf{A}t} = [e^{\lambda_1t}(\mathbf{A} - \lambda_2\mathbf{I}) - e^{\lambda_2t}(\mathbf{A} - \lambda_1\mathbf{I})]/(\lambda_1 - \lambda_2)$ |
| Discretisation (ZOH) | $\boldsymbol{\Phi} = e^{\mathbf{A}T}$, $\boldsymbol{\Gamma} = \mathbf{A}^{-1}(e^{\mathbf{A}T} - \mathbf{I})\mathbf{B}$ |
| Series RLC | $\alpha = R/2L$, $\omega_n = 1/\sqrt{LC}$, $\zeta = \tfrac R2\sqrt{C/L}$, $Q = \tfrac1R\sqrt{L/C}$ |
| Parallel RLC | $\alpha = 1/(2RC)$, $\zeta = \tfrac{1}{2R}\sqrt{L/C}$, $Q = R\sqrt{C/L}$ |
| Switching | $v_C(0^+) = v_C(0^-)$, $i_L(0^+) = i_L(0^-)$; $v_C'(0^+) = i_C(0^+)/C$, $i_L'(0^+) = v_L(0^+)/L$ |

---

## Further reading

- E. Kreyszig, *Advanced Engineering Mathematics*, 10th ed., Wiley — Chapters 2–4 (constant-coefficient equations, undetermined coefficients, systems and the matrix exponential).
- W. E. Boyce, R. C. DiPrima and D. B. Meade, *Elementary Differential Equations and Boundary Value Problems*, 11th ed., Wiley — Chapters 3, 4 and 7, including fundamental matrices and the matrix exponential.
- C. K. Alexander and M. N. O. Sadiku, *Fundamentals of Electric Circuits*, 7th ed., McGraw-Hill — Chapters 7–8 on first- and second-order circuits, switching and initial conditions.
- C.-T. Chen, *Linear System Theory and Design*, 4th ed., Oxford University Press — state-space models, the matrix exponential and discretisation in depth.
