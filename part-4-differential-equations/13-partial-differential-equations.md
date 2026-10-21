# Chapter 13 — Partial Differential Equations

[← Previous](12-linear-differential-equations.md) · [Table of Contents](../00-front/table-of-contents.md) · [Next →](../part-5-transforms/14-laplace-transforms.md)

> **Prerequisites:** Ch. 6 (Multivariable Calculus), Ch. 7 (Vector Calculus), Ch. 11 and Ch. 12 (Ordinary and Linear Differential Equations); Ch. 15 (Fourier Series) is used heavily from Section 13.4 on and may be read in parallel.
>
> **Learning outcomes:** After this chapter you can
> - derive the telegrapher's equations of a transmission line and recognise the wave, diffusion and Laplace equations as its limiting cases;
> - classify a second-order linear PDE as elliptic, parabolic or hyperbolic, find its characteristics and say what kind of auxiliary data makes the problem well posed;
> - state Dirichlet, Neumann, Robin and periodic boundary conditions, homogenise non-homogeneous ones, and translate circuit terminations into boundary conditions;
> - solve the heat, wave and Laplace equations on intervals, rectangles, disks and annuli by separation of variables and Fourier-series expansion;
> - use d'Alembert's solution, characteristics and the method of images to follow pulses and reflections on lines and strings;
> - apply these solutions to electrostatic fields, semiconductor diffusion, skin effect and transient heat transfer, and estimate how many series terms a given accuracy needs.

An ordinary differential equation describes a quantity that depends on one variable — usually time. Most physical quantities in electrical engineering depend on *where* as well as *when*: the voltage along a cable, the temperature across a power transistor die, the potential between the plates of a capacitor, the concentration of boron atoms below a silicon surface, the current density inside a conductor carrying high-frequency current. Each of these obeys a **partial differential equation** (PDE), a relation between an unknown function of several variables and its partial derivatives.

This chapter concentrates on the three classical second-order linear PDEs — the **heat (diffusion) equation**, the **wave equation** and **Laplace's equation** — because between them they cover almost every field problem met in an EE curriculum, and because they are the prototypes of the three types into which all second-order linear PDEs fall. The main solution technique is **separation of variables**, which turns one PDE into several ODEs ([Ch. 11](11-ordinary-differential-equations.md), [Ch. 12](12-linear-differential-equations.md)) and then reassembles the solution as a **Fourier series** ([Ch. 15](../part-5-transforms/15-fourier-series.md)). For infinite domains we also use d'Alembert's travelling-wave solution and similarity solutions built from the error function. Transform methods for PDEs (Laplace transform in time, Fourier transform in space) are touched on only briefly here; they are developed in [Ch. 14](../part-5-transforms/14-laplace-transforms.md) and [Ch. 16](../part-5-transforms/16-fourier-transforms.md), and numerical (finite-difference) solution in [Ch. 22](../part-7-numerical-computational/22-numerical-methods.md).

### Motivation: the transmission line

We start with one physical system that contains all three classical equations as special cases. A two-conductor transmission line (coaxial cable, twisted pair, microstrip trace) is described by four constants **per unit length**: series resistance $R$ ($\Omega/\mathrm{m}$), series inductance $L$ ($\mathrm{H/m}$), shunt conductance $G$ ($\mathrm{S/m}$) and shunt capacitance $C$ ($\mathrm{F/m}$). Let $v(x,t)$ be the voltage between the conductors and $i(x,t)$ the current in the forward conductor at position $x$ and time $t$.

Cut out a short piece of length $\Delta x$. Its equivalent circuit is a series $R\Delta x$, $L\Delta x$ followed by a shunt $G\Delta x$, $C\Delta x$. Kirchhoff's voltage law around the piece gives

$$
v(x,t) - v(x+\Delta x,t) = R\,\Delta x\; i(x,t) + L\,\Delta x\,\frac{\partial i}{\partial t}(x,t),
$$

and Kirchhoff's current law at the node $x+\Delta x$ gives

$$
i(x,t) - i(x+\Delta x,t) = G\,\Delta x\; v(x+\Delta x,t) + C\,\Delta x\,\frac{\partial v}{\partial t}(x+\Delta x,t).
$$

Divide by $\Delta x$ and let $\Delta x\to 0$. The left-hand sides become $-\partial v/\partial x$ and $-\partial i/\partial x$:

> [!theorem] Telegrapher's equations
> $$
> \frac{\partial v}{\partial x} = -R\,i - L\,\frac{\partial i}{\partial t},\qquad
> \frac{\partial i}{\partial x} = -G\,v - C\,\frac{\partial v}{\partial t}.
> $$

These are two coupled *first-order* PDEs. To get a single equation for $v$, differentiate the first with respect to $x$ and substitute the second (the mixed partials $i_{xt} = i_{tx}$ commute for smooth $i$, [Ch. 6](../part-2-calculus/06-multivariable-calculus.md)):

$$
\frac{\partial^2 v}{\partial x^2} = -R\frac{\partial i}{\partial x} - L\frac{\partial}{\partial t}\frac{\partial i}{\partial x}
= R\Big(Gv + C\frac{\partial v}{\partial t}\Big) + L\Big(G\frac{\partial v}{\partial t} + C\frac{\partial^2 v}{\partial t^2}\Big).
$$

> [!theorem] Telegrapher's (transmission-line) equation
> $$
> \frac{\partial^2 v}{\partial x^2} = LC\,\frac{\partial^2 v}{\partial t^2} + (RC + GL)\,\frac{\partial v}{\partial t} + RG\,v .
> $$
> The current $i(x,t)$ satisfies the same equation.

Dimensional check: $LC$ has units $\mathrm{s^2/m^2}$, so $LC\,v_{tt}$ has units $\mathrm{V/m^2}$, as does $v_{xx}$. $RC$ and $GL$ have units $\mathrm{s/m^2}$, and $RG$ has units $\mathrm{m^{-2}}$. Every term is in $\mathrm{V/m^2}$.

Three limiting cases give the three classical equations:

1. **Lossless line** ($R = G = 0$): $v_{xx} = LC\,v_{tt}$, or
   $$
   \frac{\partial^2 v}{\partial t^2} = c^2\frac{\partial^2 v}{\partial x^2},\qquad c = \frac{1}{\sqrt{LC}} .
   $$
   This is the **wave equation**. Signals travel undistorted at speed $c$. For a coaxial cable with polyethylene dielectric, $c\approx 2\times 10^{8}\ \mathrm{m/s}$.
2. **RC line** ($L = G = 0$): the long submarine telegraph cables analysed by Kelvin in the 1850s, and on-chip interconnect today, are dominated by $R$ and $C$. Then
   $$
   \frac{\partial v}{\partial t} = \frac{1}{RC}\frac{\partial^2 v}{\partial x^2}.
   $$
   This is the **diffusion (heat) equation** with diffusivity $1/(RC)$ in $\mathrm{m^2/s}$. Signals do not travel; they *spread*, and the delay grows as the square of the length.
3. **DC steady state, lossless shunt** ($\partial/\partial t = 0$, $G = 0$): $v_{xx} = 0$, so $v$ is linear in $x$. In two or three space dimensions the analogous statement $\nabla^2 V = 0$ is **Laplace's equation**, which governs electrostatic potential in charge-free regions.

The rest of the chapter studies each of these, first in general and then with the specific techniques that solve them.

---

## 13.1 Classification of PDEs

### Basic vocabulary

> [!definition] Partial differential equation
> A **partial differential equation** for an unknown function $u(x_1,\dots,x_n)$ of $n\ge 2$ independent variables is a relation
> $$
> F\big(x_1,\dots,x_n,\ u,\ u_{x_1},\dots,u_{x_n},\ u_{x_1x_1},\ u_{x_1x_2},\dots\big) = 0
> $$
> among the variables, the function and finitely many of its partial derivatives. Subscripts denote partial derivatives: $u_x = \partial u/\partial x$, $u_{xt} = \partial^2 u/\partial t\,\partial x$.
>
> - The **order** is the order of the highest derivative present.
> - The PDE is **linear** if $F$ is linear in $u$ and all its derivatives, with coefficients that depend only on the independent variables.
> - A linear PDE is **homogeneous** if every term contains $u$ or one of its derivatives; a term depending only on the independent variables (a *source* or *forcing* term) makes it **non-homogeneous**.
> - The PDE is **quasilinear** if it is linear in the highest-order derivatives, with coefficients that may depend on $u$ and lower derivatives.

So the telegrapher's equation is second order, linear and homogeneous. The heat equation with internal Joule heating, $u_t = \alpha u_{xx} + q(x)/(\rho c_p)$, is linear but non-homogeneous. The inviscid Burgers equation $u_t + u\,u_x = 0$, which describes shock formation, is first order and quasilinear but not linear. Poisson's equation $\nabla^2 V = -\rho_v/\varepsilon$ is linear, second order and non-homogeneous whenever the charge density $\rho_v$ is non-zero.

Linearity is what makes the methods of this chapter work, through the same principle that underlies circuit analysis.

> [!theorem] Superposition principle
> Let $\mathcal{L}$ be a linear differential operator, so that $\mathcal{L}[c_1u_1 + c_2u_2] = c_1\mathcal{L}[u_1] + c_2\mathcal{L}[u_2]$ for constants $c_1, c_2$.
> 1. If $u_1,\dots,u_N$ solve the homogeneous equation $\mathcal{L}[u] = 0$, so does every linear combination $\sum c_k u_k$.
> 2. If $u_p$ solves $\mathcal{L}[u] = f$, then every solution has the form $u = u_p + u_h$ with $\mathcal{L}[u_h] = 0$.
> 3. Under convergence conditions (uniform convergence of the series and of the series of the derivatives that appear in $\mathcal{L}$), an *infinite* series $\sum_{k=1}^\infty c_k u_k$ of solutions of $\mathcal{L}[u]=0$ is again a solution.

Part 3 is the one that separation of variables relies on, and it is the reason convergence of Fourier series ([Ch. 15](../part-5-transforms/15-fourier-series.md)) matters to us.

### Second-order linear equations in two variables

The general second-order linear PDE in two independent variables $x$ and $y$ (one of which may be time) is

$$
a_{11}\,u_{xx} + 2a_{12}\,u_{xy} + a_{22}\,u_{yy} + b_1\,u_x + b_2\,u_y + c_0\,u = f,
$$

where the coefficients may depend on $(x,y)$. The factor 2 in front of $a_{12}$ is a convention that makes the formulas below clean. The second-order part, the **principal part**, controls the qualitative behaviour of solutions: how information travels, whether solutions are smooth, and what boundary data are appropriate. Associate with it the quadratic form $a_{11}\xi^2 + 2a_{12}\xi\eta + a_{22}\eta^2$, i.e. the symmetric matrix

$$
\mathbf{A} = \begin{pmatrix} a_{11} & a_{12}\\ a_{12} & a_{22}\end{pmatrix}.
$$

> [!definition] Type of a second-order linear PDE
> Let the **discriminant** be $\Delta = a_{12}^2 - a_{11}a_{22} = -\det\mathbf{A}$. At a point $(x,y)$ the equation is
> - **hyperbolic** if $\Delta > 0$ (eigenvalues of $\mathbf{A}$ of opposite sign),
> - **parabolic** if $\Delta = 0$ (one eigenvalue zero),
> - **elliptic** if $\Delta < 0$ (eigenvalues of the same sign).
>
> If the coefficients vary, the type can change from region to region.

The names come from the conic $a_{11}\xi^2 + 2a_{12}\xi\eta + a_{22}\eta^2 = 1$, which is a hyperbola, parabola or ellipse in the three cases. The lower-order terms do not affect the type.

The three prototypes:

| Equation | Standard form | $a_{11}, a_{12}, a_{22}$ | $\Delta$ | Type |
|---|---|---|---|---|
| Wave | $u_{tt} - c^2u_{xx} = 0$ (variables $x, t$) | $-c^2,\ 0,\ 1$ | $c^2>0$ | hyperbolic |
| Heat | $u_t - \alpha u_{xx} = 0$ | $-\alpha,\ 0,\ 0$ | $0$ | parabolic |
| Laplace | $u_{xx} + u_{yy} = 0$ | $1,\ 0,\ 1$ | $-1<0$ | elliptic |

The telegrapher's equation $v_{xx} - LC\,v_{tt} - (RC+GL)v_t - RGv = 0$ has $a_{11} = 1$, $a_{12} = 0$, $a_{22} = -LC$, so $\Delta = LC > 0$: it is hyperbolic whenever $L\neq 0$. The loss terms are lower order and do not change the type — they cause attenuation and dispersion, but signals still have a finite front speed $1/\sqrt{LC}$. In the RC limit $L=0$ the second-order time derivative disappears, $\Delta = 0$, and the type jumps to parabolic: information then spreads instantly (though with exponentially small amplitude) to all $x$.

### What the type tells you

| | Elliptic | Parabolic | Hyperbolic |
|---|---|---|---|
| Prototype physics | equilibrium, steady state, electrostatics | diffusion, heat conduction, RC lines | waves, lossless lines, vibrations |
| Typical data | boundary values on a closed curve | initial value + boundary values | initial value and initial velocity + boundary values |
| Speed of information | — (no time) | infinite | finite, along characteristics |
| Smoothness | solutions as smooth as the coefficients, even with rough boundary data | rough initial data smoothed instantly for $t>0$ | discontinuities propagate, never smoothed |
| Time reversal | — | ill posed backwards | reversible |

### Characteristics and canonical forms

The type can be understood through **characteristic curves**. Ask: along which curves $y = y(x)$ is it impossible to determine the second derivatives of $u$ from the values of $u$ and its first derivatives on the curve? Writing $d(u_x) = u_{xx}\,dx + u_{xy}\,dy$ and $d(u_y) = u_{xy}\,dx + u_{yy}\,dy$ together with the PDE gives a $3\times3$ linear system for $(u_{xx},u_{xy},u_{yy})$. Its determinant vanishes when

$$
a_{11}\,(dy)^2 - 2a_{12}\,dx\,dy + a_{22}\,(dx)^2 = 0,
\qquad\text{i.e.}\qquad
\frac{dy}{dx} = \frac{a_{12}\pm\sqrt{\Delta}}{a_{11}} .
$$

> [!definition] Characteristics
> The solution curves of $dy/dx = (a_{12}\pm\sqrt{\Delta})/a_{11}$ are the **characteristics** of the PDE. A hyperbolic equation has two real families, a parabolic equation one, and an elliptic equation none (the roots are complex).

Data given on a characteristic cannot be continued off it uniquely; this is why initial data for the wave equation are given on a line $t = 0$ that is *not* characteristic, and why discontinuities in hyperbolic problems travel *along* characteristics.

Changing to characteristic coordinates simplifies the principal part.

> [!theorem] Canonical forms
> With constant coefficients (or locally, with variable ones) there is a change of variables $\xi = \xi(x,y)$, $\eta = \eta(x,y)$ that reduces the principal part to
> - hyperbolic: $u_{\xi\eta} + (\text{lower order}) = 0$, with $\xi = \text{const}$, $\eta = \text{const}$ the two characteristic families;
> - parabolic: $u_{\eta\eta} + (\text{lower order}) = 0$, with $\xi = \text{const}$ the single characteristic family and $\eta$ any independent function;
> - elliptic: $u_{\xi\xi} + u_{\eta\eta} + (\text{lower order}) = 0$, with $\xi, \eta$ the real and imaginary parts of the complex characteristic.

*Proof sketch (hyperbolic, constant coefficients).* Let $m_1\ne m_2$ be the real roots of $a_{11}m^2 - 2a_{12}m + a_{22} = 0$; these are exactly the characteristic slopes. Since $m_1 + m_2 = 2a_{12}/a_{11}$ and $m_1m_2 = a_{22}/a_{11}$, the principal part factors as

$$
a_{11}\partial_x^2 + 2a_{12}\partial_x\partial_y + a_{22}\partial_y^2 = a_{11}(\partial_x + m_1\partial_y)(\partial_x + m_2\partial_y).
$$

Put $\xi = y - m_1x$ and $\eta = y - m_2x$. By the chain rule, $\partial_x + m_1\partial_y = (m_1 - m_2)\,\partial_\eta$ and $\partial_x + m_2\partial_y = (m_2 - m_1)\,\partial_\xi$ (each factor kills functions of "its own" characteristic variable). The principal part is therefore $-a_{11}(m_1-m_2)^2\,\partial_\xi\partial_\eta$, a non-zero constant times $\partial_\xi\partial_\eta$. $\square$

> [!example] Example 13.1 — Classifying and solving a constant-coefficient equation
> Classify $u_{xx} + 4u_{xy} + 3u_{yy} = 0$ and find its general solution.
>
> **Type.** $a_{11} = 1$, $2a_{12} = 4 \Rightarrow a_{12} = 2$, $a_{22} = 3$. $\Delta = 4 - 3 = 1 > 0$: hyperbolic everywhere.
>
> **Characteristics.** $dy/dx = (2\pm 1)/1$, so $dy/dx = 3$ or $1$. The families are $y - 3x = \text{const}$ and $y - x = \text{const}$.
>
> **Canonical form.** Put $\xi = y - x$, $\eta = y - 3x$. Then $u_x = -u_\xi - 3u_\eta$, $u_y = u_\xi + u_\eta$, and
> $$
> \begin{aligned}
> u_{xx} &= u_{\xi\xi} + 6u_{\xi\eta} + 9u_{\eta\eta},\\
> u_{xy} &= -u_{\xi\xi} - 4u_{\xi\eta} - 3u_{\eta\eta},\\
> u_{yy} &= u_{\xi\xi} + 2u_{\xi\eta} + u_{\eta\eta}.
> \end{aligned}
> $$
> So $u_{xx} + 4u_{xy} + 3u_{yy} = (1 - 4 + 3)u_{\xi\xi} + (6 - 16 + 6)u_{\xi\eta} + (9 - 12 + 3)u_{\eta\eta} = -4u_{\xi\eta}$. The PDE becomes $u_{\xi\eta} = 0$.
>
> **General solution.** Integrating in $\eta$ gives $u_\xi = h(\xi)$; integrating in $\xi$ gives $u = F(\xi) + G(\eta)$:
> $$
> u(x,y) = F(y - x) + G(y - 3x)
> $$
> for arbitrary twice-differentiable $F$, $G$.
>
> **Check.** For $u = G(y-3x)$: $u_{xx} = 9G''$, $u_{xy} = -3G''$, $u_{yy} = G''$, and $9 - 12 + 3 = 0$. For $u = F(y-x)$: $1 - 4 + 3 = 0$. ✓

Notice the analogy with ODEs: a second-order ODE has two arbitrary *constants*; a second-order PDE in two variables typically has two arbitrary *functions*. This is why PDE solutions are pinned down by functions (initial profiles, boundary profiles) rather than by numbers.

> [!example] Example 13.2 — A parabolic equation
> Classify $u_{xx} + 2u_{xy} + u_{yy} = 0$ and find its general solution.
>
> $a_{11} = a_{12} = a_{22} = 1$, $\Delta = 1 - 1 = 0$: parabolic. The single characteristic family is $dy/dx = 1$, i.e. $\xi = y - x$. Take $\eta = x$ as the second variable. Then $\partial_x = -\partial_\xi + \partial_\eta$, $\partial_y = \partial_\xi$, and the principal part $(\partial_x + \partial_y)^2$ becomes $\partial_\eta^2$. The PDE is $u_{\eta\eta} = 0$, so $u = F(\xi) + \eta\,G(\xi)$:
> $$
> u(x,y) = F(y - x) + x\,G(y - x).
> $$
> **Check** the second piece: with $u = xG(y-x)$, $u_{xx} = -2G' + xG''$, $u_{xy} = G' - xG''$, $u_{yy} = xG''$. Then $u_{xx} + 2u_{xy} + u_{yy} = (-2G' + 2G') + (x - 2x + x)G'' = 0$. ✓
>
> The repeated factor $x$ is the PDE analogue of the repeated-root solution $te^{rt}$ of a constant-coefficient ODE ([Ch. 12](12-linear-differential-equations.md)).

> [!example] Example 13.3 — A PDE that changes type
> The **Tricomi equation** $y\,u_{xx} + u_{yy} = 0$ models transonic flow. Here $a_{11} = y$, $a_{12} = 0$, $a_{22} = 1$, so $\Delta = -y$.
> - For $y > 0$: $\Delta < 0$, elliptic (subsonic region).
> - For $y = 0$: parabolic (sonic line).
> - For $y < 0$: hyperbolic, with characteristics $dy/dx = \pm\sqrt{-y}/y$, i.e. $\tfrac{2}{3}(-y)^{3/2} = \pm x + \text{const}$.
>
> The same physical device can therefore require different kinds of boundary data in different regions — an important warning when setting up numerical solvers.

> [!example] Example 13.4 — Classifying the RC-line and the lossy line
> (a) On-chip interconnect with $R = 2\times10^{5}\ \Omega/\mathrm{m}$ (that is, $0.2\ \Omega/\mu\mathrm{m}$) and $C = 2\times10^{-10}\ \mathrm{F/m}$ ($0.2\ \mathrm{fF}/\mu\mathrm{m}$) and negligible $L$: $v_t = v_{xx}/(RC)$ with $RC = 4\times10^{-5}\ \mathrm{s/m^2}$. Parabolic, diffusivity $1/(RC) = 2.5\times10^{4}\ \mathrm{m^2/s}$. A diffusive "front" travels a distance $\ell$ in a time of order $RC\ell^2$; for $\ell = 1\ \mathrm{mm}$ this is $4\times10^{-11}\ \mathrm{s} = 40\ \mathrm{ps}$, for $\ell = 10\ \mathrm{mm}$ it is $4\ \mathrm{ns}$ — a hundred times longer for ten times the length.
>
> (b) If the same line has $L = 4\times10^{-7}\ \mathrm{H/m}$, the equation is hyperbolic with front speed $1/\sqrt{LC} = 1/\sqrt{8\times10^{-17}} \approx 1.12\times10^{8}\ \mathrm{m/s}$. The ratio of the loss term to the inertial term at frequency $\omega$ is $RC\,\omega/(LC\,\omega^2) = R/(\omega L)$; with $\omega = 2\pi\times 1\ \mathrm{GHz}$, $\omega L \approx 2.5\times10^{3}\ \Omega/\mathrm{m}$, so $R/(\omega L)\approx 80$ and the line behaves diffusively at 1 GHz despite being formally hyperbolic. Type tells you the *ultimate* (high-frequency) behaviour; the size of the terms tells you what you will actually see.

> [!warning] Common pitfall — type is a property of the principal part only
> Do not include first-order or zeroth-order terms when computing $\Delta$. $u_t = \alpha u_{xx} - \beta u_x$ (diffusion with drift, as for minority carriers in an electric field) is just as parabolic as the plain heat equation. Also remember the factor 2: in $u_{xx} + 4u_{xy} + 3u_{yy}$, $a_{12} = 2$, not 4.

> [!tip] Engineering note — choosing a numerical method by type
> Field solvers exploit the type. Elliptic problems (electrostatic capacitance extraction) are solved as one large sparse linear system for all unknowns at once. Parabolic problems (thermal transients, process simulation of diffusion) march in time but need implicit schemes to avoid tiny time steps. Hyperbolic problems (FDTD electromagnetic simulation, signal-integrity of long lines) march in time with a step limited by the Courant condition $c\,\Delta t \le \Delta x$, which is a statement that the numerical domain of dependence must contain the characteristic one ([Ch. 22](../part-7-numerical-computational/22-numerical-methods.md)).

---

## 13.2 Boundary and initial conditions

A PDE on its own has far too many solutions — Example 13.1 produced two *arbitrary functions*. To single out the physically relevant one we must say what happens on the boundary of the region and, for evolution equations, at the initial instant. Getting these auxiliary conditions right is half the modelling work, and the wrong kind or the wrong number of conditions produces either no solution or a meaningless one.

### Initial conditions

For an equation that is $k$th order in time, we prescribe the function and its first $k-1$ time derivatives at $t = 0$, exactly as for an ODE.

- Heat equation ($u_t$ first order in $t$): one initial condition, $u(x,0) = f(x)$ — the initial temperature profile.
- Wave and telegrapher's equations ($u_{tt}$): two, $u(x,0) = f(x)$ and $u_t(x,0) = g(x)$ — initial displacement and initial velocity. On a transmission line, specifying $v(x,0)$ and $i(x,0)$ is equivalent, because the second telegrapher equation gives $v_t$ from $i_x$ and $v$: $C\,v_t = -i_x - Gv$.
- Laplace's equation has no time variable and takes no initial condition.

### Boundary conditions

Let the spatial region be $\Omega$ with boundary $\partial\Omega$ and outward unit normal $\hat{\mathbf{n}}$, and write $\partial u/\partial n = \hat{\mathbf n}\cdot\nabla u$ for the outward normal derivative. In one dimension $\Omega = (0,L)$, $\partial\Omega = \{0, L\}$, and $\partial u/\partial n$ is $-u_x$ at $x=0$ and $+u_x$ at $x=L$.

> [!definition] Classical boundary conditions
> - **Dirichlet (first kind):** the value is prescribed, $u = h$ on $\partial\Omega$. *Examples:* an electrode held at a fixed potential; a surface clamped to a heat sink at fixed temperature; a short-circuited transmission-line end, $v(\ell,t) = 0$.
> - **Neumann (second kind):** the normal derivative is prescribed, $\partial u/\partial n = h$. *Examples:* a thermally insulated surface, $\partial u/\partial n = 0$ (no heat flux); a prescribed heat flux $-k\,\partial u/\partial n = q_s$; the free end of a string; an open-circuited line end, $i(\ell,t) = 0$, which by the first telegrapher equation (lossless case) means $v_x(\ell,t) = -L\,i_t = 0$.
> - **Robin (third kind, mixed):** a linear combination, $\beta\,u + \partial u/\partial n = h$. *Example:* Newton cooling at a surface in contact with a fluid at temperature $u_\infty$, $-k\,\partial u/\partial n = h_c(u - u_\infty)$, where $h_c$ is the convective heat-transfer coefficient ($\mathrm{W\,m^{-2}K^{-1}}$).
> - **Periodic:** on an interval that represents a closed loop (angle $\phi$ around a ring or a cylinder), $u(-\pi) = u(\pi)$ and $u_\phi(-\pi) = u_\phi(\pi)$.
>
> A boundary condition is **homogeneous** if $h \equiv 0$, so that the zero function satisfies it and sums of functions satisfying it also satisfy it.

Different parts of the boundary may carry different kinds of condition. For a rectangle one might hold one edge at a potential (Dirichlet) while the others are insulating symmetry planes (Neumann).

> [!definition] Initial–boundary-value problem
> An **initial–boundary-value problem** (IBVP) consists of an evolution PDE on $\Omega\times(0,\infty)$, initial conditions at $t=0$, and boundary conditions on $\partial\Omega$ for all $t>0$. A **boundary-value problem** (BVP) consists of a steady-state PDE (typically elliptic) on $\Omega$ with conditions on $\partial\Omega$.

The model problem for the heat equation on a rod of length $L$ with ends held at zero temperature is

$$
\begin{aligned}
&u_t = \alpha\,u_{xx}, && 0<x<L,\ t>0,\\
&u(0,t) = 0,\quad u(L,t) = 0, && t>0,\\
&u(x,0) = f(x), && 0\le x\le L .
\end{aligned}
$$

### Well-posedness

> [!definition] Well-posed problem (Hadamard)
> A problem is **well posed** if
> 1. a solution **exists**,
> 2. the solution is **unique**, and
> 3. the solution **depends continuously on the data**: small changes in the initial and boundary data produce small changes in the solution.

Condition 3 is the one engineers care most about, because measured data always carry error. A problem that fails it is **ill posed**, and any computed answer is at the mercy of noise.

The right auxiliary data depend on the type:

| Type | Well-posed data | Ill-posed data |
|---|---|---|
| Elliptic (Laplace) | one condition (D, N or R) at each point of a *closed* boundary | Cauchy data ($u$ and $u_n$) on part of the boundary only |
| Parabolic (heat) | initial value + one condition at each boundary point, solved *forward* in time | solving *backward* in time from a final state |
| Hyperbolic (wave) | $u$ and $u_t$ at $t=0$ + one condition at each boundary point | boundary values on a closed curve in the $(x,t)$ plane |

> [!example] Example 13.5 — Hadamard's example: an ill-posed Laplace problem
> Consider $u_{xx} + u_{yy} = 0$ for $y > 0$ with **Cauchy data** on the line $y = 0$: $u(x,0) = 0$, $u_y(x,0) = \dfrac{\sin(nx)}{n}$.
>
> Try $u = X(x)Y(y)$ with $X = \sin(nx)$; then $Y'' = n^2Y$, $Y(0) = 0$, $Y'(0) = 1/n$, so $Y = \sinh(ny)/n^2$:
> $$
> u_n(x,y) = \frac{\sin(nx)\,\sinh(ny)}{n^2}.
> $$
> As $n\to\infty$ the data tend to zero uniformly ($\lvert u_y(x,0)\rvert\le 1/n$), and the zero data have solution $u\equiv0$. But at $y = 1$, $\max_x\lvert u_n\rvert = \sinh(n)/n^2 \approx e^n/(2n^2)$; for $n = 20$ this is about $6.1\times10^{5}$. Data of size $0.05$ produce a solution of size $6\times10^{5}$. Continuous dependence fails: Laplace's equation with Cauchy data is ill posed.
>
> *Engineering meaning:* trying to reconstruct the potential deep inside a region from measurements of $V$ and $\partial V/\partial n$ on one face (an "inverse" problem, as in impedance tomography) amplifies measurement noise exponentially and needs regularisation.

Backward heat conduction is ill posed for the same reason: a mode $\sin(n\pi x/L)$ decays like $e^{-\alpha(n\pi/L)^2 t}$ forward in time, so running time backwards multiplies high-frequency noise by $e^{+\alpha(n\pi/L)^2 t}$ (see Section 13.4).

### How many conditions?

A useful rule: a PDE that is of order $2$ in a spatial variable needs **one condition at each end** of that variable's interval; a PDE of order $k$ in time needs $k$ **initial** conditions. Heat equation on $(0,L)$: two boundary conditions + one initial condition. Wave equation on $(0,L)$: two boundary conditions + two initial conditions. Laplace on a rectangle: one condition on each of the four sides. On an infinite interval the boundary conditions are replaced by a growth or decay condition ("$u$ bounded as $\lvert x\rvert\to\infty$").

### Compatibility

The initial and boundary data must agree where they meet. For the model heat problem, continuity at the corner $(x,t) = (0,0)$ requires $f(0) = 0$. If $f(0)\ne0$ (a rod at $100\ ^\circ\mathrm{C}$ suddenly clamped to $0\ ^\circ\mathrm{C}$ ends), the solution still exists and is smooth for $t>0$, but the series solution converges non-uniformly near the corner and the boundary condition is attained only as $t\to0^+$. For the wave equation, incompatible data produce a discontinuity (a "kink") that travels along the characteristic from the corner and never smooths out.

### Homogenising non-homogeneous boundary conditions

Separation of variables (Section 13.3) needs **homogeneous** boundary conditions. When the boundary data are constant (or independent of time), subtract a steady-state solution.

> [!theorem] Steady-state subtraction
> Suppose $u_t = \alpha u_{xx}$ on $0<x<L$ with $u(0,t) = T_1$, $u(L,t) = T_2$ (constants) and $u(x,0) = f(x)$. Let $u_s(x)$ be the steady state: $u_s'' = 0$, $u_s(0) = T_1$, $u_s(L) = T_2$, so
> $$
> u_s(x) = T_1 + (T_2 - T_1)\frac{x}{L}.
> $$
> Then the **transient** $w = u - u_s$ satisfies $w_t = \alpha w_{xx}$, $w(0,t) = w(L,t) = 0$, $w(x,0) = f(x) - u_s(x)$ — a problem with homogeneous boundary conditions.

The proof is one line: $u_s$ is independent of $t$ and has $u_s''=0$, so it satisfies the PDE, and subtracting it removes the boundary values. The same idea works for Laplace's equation (split a problem with data on all four sides into four problems with data on one side each, Section 13.6), for the wave equation, and for time-dependent boundary data (subtract any function $\psi(x,t)$ that matches the boundary data; the price is a source term $-\psi_t + \alpha\psi_{xx}$ in the PDE for $w$, handled by eigenfunction expansion in Section 13.7).

> [!example] Example 13.6 — Homogenising a Neumann–Dirichlet problem
> A rod $0<x<L$ has its end $x = 0$ held at $T_0$ and receives a constant heat flux $q_s$ ($\mathrm{W/m^2}$) into the end $x = L$: $k\,u_x(L,t) = q_s$. Find a steady state that removes both boundary conditions.
>
> Steady state: $u_s'' = 0 \Rightarrow u_s = A + Bx$. $u_s(0) = T_0 \Rightarrow A = T_0$. $k u_s'(L) = q_s \Rightarrow B = q_s/k$. So
> $$
> u_s(x) = T_0 + \frac{q_s}{k}\,x .
> $$
> Then $w = u - u_s$ satisfies $w(0,t) = 0$, $w_x(L,t) = 0$, $w(x,0) = f(x) - T_0 - q_sx/k$.
>
> **Check** units: $q_s/k$ has units $(\mathrm{W/m^2})/(\mathrm{W\,m^{-1}K^{-1}}) = \mathrm{K/m}$, a temperature gradient. ✓ Numerically, a copper rod ($k = 400\ \mathrm{W\,m^{-1}K^{-1}}$, $L = 0.1\ \mathrm{m}$) receiving $q_s = 10^5\ \mathrm{W/m^2}$ ($10\ \mathrm{W/cm^2}$) settles with its hot end $q_sL/k = 25\ \mathrm{K}$ above $T_0$.
>
> *Pitfall:* with Neumann conditions at **both** ends, $k u_x(0,t) = q_0$, $k u_x(L,t) = q_L$, a linear steady state exists only if $q_0 = q_L$ (net heat flow in = out). Otherwise energy accumulates and there is no steady state. Subtract instead
> $$
> \psi(x,t) = \frac{q_L - q_0}{2kL}\,x^2 + \frac{q_0}{k}\,x + \frac{(q_L-q_0)\,\alpha}{kL}\,t,
> $$
> which matches both flux conditions and satisfies $\psi_t = \alpha\psi_{xx}$ (both sides equal $\alpha(q_L-q_0)/(kL)$), so the uniform temperature rises linearly in time.

> [!example] Example 13.7 — Circuit terminations as boundary conditions
> A lossless line of length $\ell$ ($R = G = 0$) is driven at $x = 0$ and terminated at $x = \ell$. Express the termination as a condition on $v(x,t)$ alone.
>
> The telegrapher's equation $v_x = -L\,i_t$ links $i$ and $v$.
> - **Short circuit:** $v(\ell,t) = 0$ — Dirichlet.
> - **Open circuit:** $i(\ell,t) = 0$ for all $t$, so $i_t(\ell,t) = 0$ and $v_x(\ell,t) = 0$ — Neumann.
> - **Resistor $R_L$:** $v(\ell,t) = R_L\,i(\ell,t)$. Differentiate in time: $v_t = R_L i_t = -(R_L/L)\,v_x$. So
>   $$
>   v_x(\ell,t) + \frac{L}{R_L}\,v_t(\ell,t) = 0,
>   $$
>   a boundary condition that involves a *time* derivative. It is homogeneous, but it is not of the Robin type in $x$ alone; we solve it with travelling waves in Section 13.5 and find the reflection coefficient $\Gamma = (R_L - Z_0)/(R_L + Z_0)$.
> - **Source with internal resistance $R_s$ at $x=0$:** $v(0,t) = v_s(t) - R_s\,i(0,t)$ — non-homogeneous, time-dependent.
>
> *Check limits:* $R_L\to\infty$ gives $v_x = 0$ (open). $R_L\to0$ forces $v_t(\ell,t) = 0$, i.e. $v(\ell,t)$ stays at its initial value $0$ (short). ✓ Dimensional check: $L/R_L$ has units $(\mathrm{H/m})/\Omega = \mathrm{s/m}$, so $(L/R_L)v_t$ has units $\mathrm{V/m}$, like $v_x$. ✓

> [!warning] Common pitfall — Neumann data alone do not fix the solution
> For the steady state (Laplace) with pure Neumann conditions on the whole boundary, $u$ is determined only up to an additive constant, and a solution exists only if the net flux vanishes: $\oint_{\partial\Omega}\partial u/\partial n\,ds = 0$ (by the divergence theorem, [Ch. 7](../part-2-calculus/07-vector-calculus.md)). In a field solver this shows up as a singular matrix; fix one node's potential (a "ground" reference) to remove the null space, exactly as you choose a ground node in nodal circuit analysis.

> [!tip] Engineering note — symmetry planes are Neumann boundaries
> A plane of mirror symmetry in a thermal or electrostatic problem carries zero normal flux, so it can be replaced by an insulating (Neumann, $\partial u/\partial n = 0$) boundary. A plane of *antisymmetry* (potential $+V$ on one side mirrored by $-V$ on the other) is an equipotential at $0$, a Dirichlet boundary. Modelling half or a quarter of a symmetric device this way cuts the solver's work by 2× or 4×.

### Uniqueness by the energy method

Uniqueness for linear problems reduces to showing that the *difference* of two solutions, which solves the problem with zero data, is zero.

> [!theorem] Uniqueness for the heat IBVP
> The problem $u_t = \alpha u_{xx} + q(x,t)$ on $0<x<L$, $u(x,0) = f(x)$, with Dirichlet, Neumann or Robin ($\beta\ge0$) conditions at each end, has at most one smooth solution.

*Proof.* Let $w$ be the difference of two solutions. Then $w_t = \alpha w_{xx}$ with zero initial value and homogeneous boundary conditions. Define the "energy" $E(t) = \tfrac12\int_0^L w^2\,dx\ge 0$. Then, integrating by parts,

$$
E'(t) = \int_0^L w\,w_t\,dx = \alpha\int_0^L w\,w_{xx}\,dx = \alpha\big[w\,w_x\big]_0^L - \alpha\int_0^L w_x^2\,dx .
$$

The boundary term is zero for Dirichlet and Neumann conditions and $\le 0$ for Robin conditions with $\beta \ge 0$ (at $x = L$, $w_x = -\beta w$ gives $w\,w_x = -\beta w^2\le0$; at $x=0$ the sign flip of the outward normal gives the same conclusion). Hence $E'(t)\le 0$. Since $E(0) = 0$ and $E\ge0$, $E(t) = 0$ for all $t$, so $w\equiv0$. $\square$

The same argument with $E = \tfrac12\int(w_t^2 + c^2w_x^2)\,dx$ proves uniqueness for the wave equation (there $E$ is conserved, not just non-increasing), and with $E = \int\lvert\nabla w\rvert^2\,dV$ for Laplace's equation.

---

## 13.3 Separation of variables

Separation of variables is the workhorse method for linear PDEs on simple domains (intervals, rectangles, disks, boxes, cylinders, spheres). Its idea is to look for solutions in which the dependence on each variable is *independent* of the others,

$$
u(x,t) = X(x)\,T(t),
$$

so that the PDE splits into ODEs, one per variable. Each product solution is too special to match arbitrary data, but because the problem is linear we can superpose infinitely many of them and choose the coefficients to match the data. The method is the PDE counterpart of **modal analysis**: a coupled system $\dot{\mathbf{x}} = \mathbf{A}\mathbf{x}$ is solved by expanding in the eigenvectors of $\mathbf{A}$ ([Ch. 10](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md), [Ch. 11](11-ordinary-differential-equations.md)); here the "matrix" is the differential operator $d^2/dx^2$ with its boundary conditions, and the "eigenvectors" are functions.

### The procedure

> [!definition] Separation of variables — the five steps
> 1. **Homogenise.** Make the PDE and the boundary conditions in all but one variable homogeneous (Section 13.2).
> 2. **Separate.** Substitute $u = X(x)T(t)$ (or $X(x)Y(y)$, $R(r)\Phi(\phi)$), divide by $XT$ and put all $x$-dependence on one side. A function of $x$ equal to a function of $t$ must be a constant, the **separation constant** $-\lambda$.
> 3. **Solve the eigenvalue problem.** The homogeneous boundary conditions become conditions on $X$ alone, giving a **Sturm–Liouville problem** $X'' + \lambda X = 0$ plus boundary conditions. Find all **eigenvalues** $\lambda_n$ for which a non-trivial solution exists and the corresponding **eigenfunctions** $X_n$.
> 4. **Solve the remaining ODE** for each $\lambda_n$, giving $T_n(t)$ and the **product solutions** $u_n = X_n(x)T_n(t)$.
> 5. **Superpose and match.** Write $u = \sum_n c_n X_n(x)T_n(t)$ and determine $c_n$ from the remaining (initial or non-homogeneous boundary) data using orthogonality of the $X_n$ — a Fourier-series computation (Section 13.7).

### Separating the heat equation

Take $u_t = \alpha u_{xx}$ on $0<x<L$ with $u(0,t) = u(L,t) = 0$. Substituting $u = X(x)T(t)$ gives $X\,T' = \alpha X''\,T$. Dividing by $\alpha XT$:

$$
\frac{T'(t)}{\alpha\,T(t)} = \frac{X''(x)}{X(x)} .
$$

The left side depends only on $t$, the right only on $x$. If we change $x$ while holding $t$ fixed, the left side does not change, so neither can the right; hence both equal a constant, which we call $-\lambda$:

$$
X'' + \lambda X = 0,\qquad T' + \alpha\lambda\,T = 0 .
$$

The boundary conditions $u(0,t) = X(0)T(t) = 0$ for all $t$ force $X(0) = 0$ (otherwise $T\equiv0$ and $u\equiv0$). Likewise $X(L) = 0$. This is the eigenvalue problem already met in [Ch. 11](11-ordinary-differential-equations.md):

> [!theorem] Dirichlet eigenvalue problem
> The problem $X'' + \lambda X = 0$, $X(0) = X(L) = 0$ has non-trivial solutions only for
> $$
> \lambda_n = \Big(\frac{n\pi}{L}\Big)^2,\qquad X_n(x) = \sin\frac{n\pi x}{L},\qquad n = 1, 2, 3,\dots
> $$

*Proof.* Consider the sign of $\lambda$.

- $\lambda < 0$: write $\lambda = -\mu^2$, $\mu>0$. Then $X = A\cosh\mu x + B\sinh\mu x$. $X(0) = A = 0$; $X(L) = B\sinh\mu L = 0$ with $\sinh\mu L\ne0$, so $B = 0$. Trivial.
- $\lambda = 0$: $X = A + Bx$; $X(0) = 0\Rightarrow A = 0$, $X(L) = BL = 0\Rightarrow B = 0$. Trivial.
- $\lambda > 0$: $\lambda = \mu^2$, $X = A\cos\mu x + B\sin\mu x$. $X(0) = A = 0$. $X(L) = B\sin\mu L = 0$, so a non-zero $B$ requires $\sin\mu L = 0$, i.e. $\mu L = n\pi$. $\square$

With $\lambda_n$ known, $T_n = e^{-\alpha\lambda_n t}$, and the product solutions are

$$
u_n(x,t) = e^{-\alpha(n\pi/L)^2 t}\,\sin\frac{n\pi x}{L}.
$$

Each is a spatial sine pattern that keeps its shape and decays exponentially, faster for higher $n$. Section 13.4 assembles them into the full solution.

### Common eigenvalue problems

The boundary conditions change the eigenfunctions. The following table covers nearly every case met in practice for $X'' + \lambda X = 0$ on $0<x<L$ (on $-\pi<\phi<\pi$ for the periodic case).

| Boundary conditions | Eigenvalues $\lambda_n$ | Eigenfunctions $X_n$ | $n$ |
|---|---|---|---|
| $X(0) = 0,\ X(L) = 0$ (D–D) | $(n\pi/L)^2$ | $\sin(n\pi x/L)$ | $1,2,\dots$ |
| $X'(0) = 0,\ X'(L) = 0$ (N–N) | $(n\pi/L)^2$ | $\cos(n\pi x/L)$ | $0,1,2,\dots$ |
| $X(0) = 0,\ X'(L) = 0$ (D–N) | $\big((2n-1)\pi/(2L)\big)^2$ | $\sin\big((2n-1)\pi x/(2L)\big)$ | $1,2,\dots$ |
| $X'(0) = 0,\ X(L) = 0$ (N–D) | $\big((2n-1)\pi/(2L)\big)^2$ | $\cos\big((2n-1)\pi x/(2L)\big)$ | $1,2,\dots$ |
| periodic on $(-\pi,\pi)$ | $n^2$ | $1$; $\cos n\phi$, $\sin n\phi$ (two per $n\ge1$) | $0,1,2,\dots$ |
| $X(0)=0,\ X'(L) + hX(L) = 0$ (D–R, $h>0$) | $\mu_n^2$, $\tan(\mu_nL) = -\mu_n/h$ | $\sin(\mu_n x)$ | $1,2,\dots$ |

Note that N–N admits $\lambda_0 = 0$ with the **constant** eigenfunction $X_0 = 1$. Physically, an insulated rod can hold any uniform temperature forever; this mode carries the average.

### Sturm–Liouville theory

All the problems in the table are instances of a general theory that guarantees the properties we need.

> [!definition] Regular Sturm–Liouville problem
> On $a<x<b$, the problem
> $$
> \big(p(x)\,X'\big)' + q(x)\,X + \lambda\,w(x)\,X = 0,
> $$
> $$
> \beta_1 X(a) + \beta_2 X'(a) = 0,\qquad \gamma_1 X(b) + \gamma_2 X'(b) = 0,
> $$
> with $p, p', q, w$ continuous, $p>0$ and $w>0$ on $[a,b]$, and $(\beta_1,\beta_2)\ne(0,0)$, $(\gamma_1,\gamma_2)\ne(0,0)$, is a **regular Sturm–Liouville problem**. $w$ is the **weight function**.

$X'' + \lambda X = 0$ is the case $p = w = 1$, $q = 0$. Bessel's equation, which arises for the circular disk and cylinder, is the case $p = r$, $q = -n^2/r$, $w = r$ on $0<r<a$ (singular at $r=0$, but the conclusions below still hold with a boundedness condition there).

> [!theorem] Sturm–Liouville theorem
> For a regular Sturm–Liouville problem:
> 1. All eigenvalues are **real**, and they form an infinite increasing sequence $\lambda_1<\lambda_2<\lambda_3<\dots\to\infty$.
> 2. Each eigenvalue has a single eigenfunction (up to a constant factor), and the $n$th eigenfunction has exactly $n-1$ zeros in $(a,b)$.
> 3. Eigenfunctions for different eigenvalues are **orthogonal with weight $w$**:
>    $$
>    \langle X_m, X_n\rangle_w = \int_a^b X_m(x)\,X_n(x)\,w(x)\,dx = 0,\qquad m\ne n .
>    $$
> 4. The eigenfunctions are **complete**: any square-integrable $f$ on $(a,b)$ has an expansion $f = \sum c_nX_n$ converging in the mean-square sense, with **generalised Fourier coefficients**
>    $$
>    c_n = \frac{\langle f, X_n\rangle_w}{\langle X_n, X_n\rangle_w} = \frac{\int_a^b f\,X_n\,w\,dx}{\int_a^b X_n^2\,w\,dx}.
>    $$
>    If $f$ is piecewise smooth, the series converges to $\tfrac12[f(x^+)+f(x^-)]$ at each interior point, exactly as for ordinary Fourier series ([Ch. 15](../part-5-transforms/15-fourier-series.md)).

*Proof of orthogonality and reality.* Let $\mathcal{S}[X] = (pX')' + qX$. For any two functions $X, Y$ satisfying the boundary conditions, **Lagrange's identity** gives

$$
\int_a^b \big(Y\,\mathcal{S}[X] - X\,\mathcal{S}[Y]\big)\,dx = \Big[p\,(YX' - XY')\Big]_a^b = 0,
$$

the boundary term vanishing because at each end both $X$ and $Y$ satisfy the same homogeneous condition, so the vectors $(X, X')$ and $(Y, Y')$ are both orthogonal to $(\beta_1, \beta_2)$ and are therefore parallel, making $YX' - XY' = 0$. (This "self-adjointness" is the analogue of a symmetric matrix, $\mathbf{y}^\mathsf{T}\mathbf{A}\mathbf{x} = \mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{y}$.) Now if $\mathcal{S}[X_m] = -\lambda_m wX_m$ and $\mathcal{S}[X_n] = -\lambda_n wX_n$, then

$$
0 = \int_a^b\big(X_n\mathcal{S}[X_m] - X_m\mathcal{S}[X_n]\big)dx = (\lambda_n - \lambda_m)\int_a^b X_mX_n\,w\,dx,
$$

so $\lambda_m\ne\lambda_n$ forces orthogonality. If $\lambda$ were complex with eigenfunction $X$, then $\bar\lambda$ would be an eigenvalue with eigenfunction $\bar X$, and the same identity would give $(\lambda - \bar\lambda)\int\lvert X\rvert^2w\,dx = 0$, so $\lambda = \bar\lambda$ is real. $\square$

The Rayleigh quotient tells us the sign of the eigenvalues. Multiplying $(pX')' + qX = -\lambda wX$ by $X$ and integrating by parts,

$$
\lambda = \frac{-\big[pXX'\big]_a^b + \int_a^b\big(pX'^2 - qX^2\big)dx}{\int_a^b X^2w\,dx}.
$$

For $X'' + \lambda X = 0$ with Dirichlet, Neumann or Robin ($h\ge0$) conditions, $q = 0$ and the boundary term $-[XX']_0^L$ is zero (Dirichlet, Neumann) or non-negative (Robin: $X'(L) = -hX(L)$ gives $-X(L)X'(L) = hX(L)^2\ge0$), so $\lambda\ge0$: there are no negative eigenvalues, and the case analysis "$\lambda<0$" in the proof above could have been skipped. Physically, $\lambda\ge0$ says that no mode of a passive diffusion grows.

> [!example] Example 13.8 — A Dirichlet–Neumann eigenvalue problem
> A rod $0<x<L$ is held at $0$ at $x=0$ and insulated at $x = L$. Find the eigenvalues and eigenfunctions and verify orthogonality of the first two.
>
> $X'' + \lambda X = 0$, $X(0) = 0$, $X'(L) = 0$. From the Rayleigh quotient, $\lambda\ge0$; $\lambda = 0$ gives $X = Bx$ with $X'(L) = B = 0$, trivial. For $\lambda = \mu^2 > 0$, $X(0) = 0$ gives $X = B\sin\mu x$, and $X'(L) = B\mu\cos\mu L = 0$ requires $\cos\mu L = 0$:
> $$
> \mu_n = \frac{(2n-1)\pi}{2L},\qquad X_n = \sin\frac{(2n-1)\pi x}{2L},\qquad n = 1,2,\dots
> $$
> The first eigenfunction is a quarter sine wave: zero at the clamped end, flat at the insulated end.
>
> **Orthogonality check** with $L = 1$: $\int_0^1\sin(\pi x/2)\sin(3\pi x/2)\,dx = \tfrac12\int_0^1[\cos(\pi x) - \cos(2\pi x)]dx = \tfrac12[0 - 0] = 0$. ✓ Normalisation: $\int_0^L X_n^2\,dx = L/2$.

> [!example] Example 13.9 — A Robin eigenvalue problem (convective cooling)
> A rod $0<x<1$ is held at $0$ at $x = 0$ and loses heat by convection at $x = 1$: $X'(1) + hX(1) = 0$ with $h = 1$. Find the first three eigenvalues.
>
> As before, $X = \sin\mu x$ (the $\lambda\le0$ cases are trivial since $h>0$). The Robin condition gives $\mu\cos\mu + h\sin\mu = 0$, i.e.
> $$
> \tan\mu = -\frac{\mu}{h}.
> $$
> This is **transcendental**: the eigenvalues must be found numerically, from the intersections of $y = \tan\mu$ with the line $y = -\mu/h$. One root lies in each interval $\big((n-\tfrac12)\pi,\ n\pi\big)$. Newton's method or bisection ([Ch. 22](../part-7-numerical-computational/22-numerical-methods.md)) on $g(\mu) = \sin\mu + \mu\cos\mu$ gives
> $$
> \mu_1 = 2.0288,\qquad \mu_2 = 4.9132,\qquad \mu_3 = 7.9787 .
> $$
> For large $n$, $\mu_n \approx (n-\tfrac12)\pi + \dfrac{h}{(n-\tfrac12)\pi}$; for $n=3$ this gives $7.8540 + 0.1273 = 7.9813$, close to the exact $7.9787$.
>
> **Limits.** $h\to0$ (insulated end) gives $\cos\mu = 0$, $\mu_1 = \pi/2 = 1.571$; $h\to\infty$ (end held at zero) gives $\sin\mu = 0$, $\mu_1 = \pi$. The Robin value $2.029$ lies between, as it must. The eigenfunctions $\sin\mu_nx$ are orthogonal on $(0,1)$ by the Sturm–Liouville theorem even though the $\mu_n$ are not equally spaced — a fact that is far from obvious by direct integration.

```desmos-graph
left=0; right=10; bottom=-10; top=10
---
y=\tan(x)|blue
y=-x|red
(2.0288,-2.0288)|label:mu1
(4.9132,-4.9132)|label:mu2
(7.9787,-7.9787)|label:mu3
```

*The Robin eigenvalues of Example 13.9 are the intersections of $y=\tan\mu$ (blue) with $y=-\mu/h$ (red, $h=1$); they approach $(n-\tfrac12)\pi$ from above.*

### Separating other equations

> [!example] Example 13.10 — Separating the lossy telegrapher's equation
> A line of length $\ell$ is short-circuited at both ends, $v(0,t) = v(\ell,t) = 0$. Separate $v_{xx} = LC\,v_{tt} + (RC+GL)v_t + RGv$.
>
> With $v = X(x)T(t)$:
> $$
> \frac{X''}{X} = \frac{LC\,T'' + (RC+GL)\,T' + RG\,T}{T} = -\lambda .
> $$
> The $X$ problem is the Dirichlet problem: $\lambda_n = (n\pi/\ell)^2$, $X_n = \sin(n\pi x/\ell)$. The $T$ equation becomes
> $$
> T'' + \Big(\frac{R}{L} + \frac{G}{C}\Big)T' + \frac{RG + (n\pi/\ell)^2}{LC}\,T = 0 ,
> $$
> a **damped harmonic oscillator** — the equation of a series RLC circuit ([Ch. 12](12-linear-differential-equations.md)). Each spatial mode of the line behaves like its own resonant circuit with undamped natural frequency $\omega_n = \sqrt{(RG + (n\pi/\ell)^2)/(LC)}$ and damping coefficient $\zeta_n\omega_n = \tfrac12(R/L + G/C)$.
>
> **Numbers.** $\ell = 1\ \mathrm{m}$, $L = 250\ \mathrm{nH/m}$, $C = 100\ \mathrm{pF/m}$ (so $c = 2\times10^8\ \mathrm{m/s}$, $Z_0 = 50\ \Omega$), $R = 0.5\ \Omega/\mathrm{m}$, $G = 0$. For $n = 1$: $\omega_1 = \pi/\sqrt{LC} = \pi\times 2\times10^{8} = 6.28\times10^8\ \mathrm{rad/s}$, i.e. $f_1 = 100\ \mathrm{MHz}$ (half a wavelength fits in the line). Decay rate $R/(2L) = 10^6\ \mathrm{s^{-1}}$, so the mode's quality factor is $Q = \omega_1/(2\cdot10^6) \approx 314$. This is the principle of a shorted-stub resonator.

> [!example] Example 13.11 — Separating Laplace's equation in polar coordinates
> In polar coordinates $(r,\phi)$, Laplace's equation is ([Ch. 7](../part-2-calculus/07-vector-calculus.md))
> $$
> u_{rr} + \frac{1}{r}u_r + \frac{1}{r^2}u_{\phi\phi} = 0 .
> $$
> Put $u = R(r)\Phi(\phi)$ and multiply by $r^2/(R\Phi)$:
> $$
> \frac{r^2R'' + rR'}{R} = -\frac{\Phi''}{\Phi} = \lambda .
> $$
> On a full disk or annulus $\phi$ is an angle, so $\Phi$ must be $2\pi$-periodic: $\lambda = n^2$, $\Phi = a_n\cos n\phi + b_n\sin n\phi$, $n = 0,1,2,\dots$. The radial equation $r^2R'' + rR' - n^2R = 0$ is a **Cauchy–Euler** equation ([Ch. 11](11-ordinary-differential-equations.md)): try $R = r^m$ to get $m^2 = n^2$, so
> $$
> R = \begin{cases} C_0 + D_0\ln r, & n = 0,\\ C_n r^n + D_n r^{-n}, & n\ge1. \end{cases}
> $$
> These product solutions are used in Section 13.6.

### When separation fails

Separation of variables needs three things: a separable PDE, a domain whose boundary consists of coordinate curves ($x = $ const, $r = $ const, ...), and homogeneous conditions on all but one pair of opposite boundaries (or a way to split the problem so this holds).

- **Non-separable PDE.** $u_{xx} + (x+y)u_{yy} = 0$ gives $X''Y + (x+y)XY'' = 0$, and $(x+y)$ cannot be split into a function of $x$ times or plus a function of $y$ after division. Equations with coefficients that are products, $a(x)b(y)$, often *do* separate.
- **Non-separable domain.** An L-shaped region, a microstrip cross-section with a finite-width strip, or a region bounded by a curve that is not a coordinate line requires conformal mapping ([Ch. 18](../part-5-transforms/18-complex-analysis.md)) or numerical methods.
- **Non-homogeneous conditions on every side.** Split by superposition into several problems with one non-homogeneous side each (Section 13.6).

> [!warning] Common pitfall — choosing the sign of the separation constant
> Writing the constant as $+\lambda$ rather than $-\lambda$ is not an error, but it changes which case gives the oscillatory solution, and students then reject the right case. Always examine all three cases $\lambda<0$, $\lambda=0$, $\lambda>0$ (or use the Rayleigh quotient to rule some out), and remember that N–N and periodic problems *do* have a zero eigenvalue. Losing the $\lambda_0 = 0$ constant mode loses the average value of the solution.

> [!tip] Engineering note — modes are the natural "coordinates" of a distributed system
> A distributed system (line, plate, membrane, cavity) has infinitely many degrees of freedom, but separation of variables shows that it behaves like an infinite set of *independent* lumped systems, one per eigenfunction. Engineers routinely truncate to the first few modes: a thermal model of a heat sink with one or two exponential time constants, a cavity resonator described by its lowest mode, a loudspeaker cone by its first bending modes. Sturm–Liouville orthogonality is what makes the modes decouple.

---

## 13.4 Heat equation

### Derivation

Consider a rod (or a slab, or a wire) in which heat flows only along $x$. Let $u(x,t)$ be the temperature, $\rho$ the density, $c_p$ the specific heat capacity, $k$ the thermal conductivity and $g(x,t)$ the heat generated per unit volume per unit time (Joule heating $J^2/\sigma$ in a conductor, for example). Two physical laws are needed.

1. **Conservation of energy** for the slice $(x, x+\Delta x)$ of cross-section $A$: the rate of increase of stored heat equals heat flowing in minus heat flowing out plus heat generated:
   $$
   \frac{\partial}{\partial t}\int_x^{x+\Delta x}\rho c_p\,u\,A\,ds = A\,q(x,t) - A\,q(x+\Delta x,t) + \int_x^{x+\Delta x} g\,A\,ds ,
   $$
   where $q$ is the heat flux ($\mathrm{W/m^2}$) in the $+x$ direction.
2. **Fourier's law of conduction**: heat flows down the temperature gradient, $q = -k\,u_x$.

Divide by $A\,\Delta x$, let $\Delta x\to0$ and substitute Fourier's law:

$$
\rho c_p\,u_t = (k\,u_x)_x + g .
$$

> [!definition] Heat (diffusion) equation
> For constant $k$, $\rho$, $c_p$,
> $$
> \frac{\partial u}{\partial t} = \alpha\,\frac{\partial^2u}{\partial x^2} + \frac{g}{\rho c_p},\qquad \alpha = \frac{k}{\rho c_p},
> $$
> where $\alpha$ is the **thermal diffusivity** ($\mathrm{m^2/s}$). In three dimensions, $u_t = \alpha\nabla^2u + g/(\rho c_p)$.

The same equation, with different letters, governs every process in which a conserved quantity flows down its own gradient:

| Process | Quantity $u$ | Flux law | Diffusivity |
|---|---|---|---|
| Heat conduction | temperature $T$ | Fourier, $q = -kT_x$ | $\alpha = k/(\rho c_p)$ |
| Impurity diffusion | concentration $N$ | Fick, $J = -DN_x$ | $D$ |
| Minority carriers (no field) | excess density $\Delta n$ | $J_n = qD_n\,\Delta n_x$ | $D_n$ (plus recombination term $-\Delta n/\tau_n$) |
| RC transmission line | voltage $v$ | Ohm, $i = -v_x/R$ | $1/(RC)$ |
| Magnetic field in a conductor | $H$ or $\mathbf{J}$ | Faraday + Ohm | $1/(\mu\sigma)$ |

Typical diffusivities: copper $\alpha\approx1.1\times10^{-4}\ \mathrm{m^2/s}$, silicon $8.8\times10^{-5}$, FR-4 circuit board $\approx 1.5\times10^{-7}$, water $1.4\times10^{-7}\ \mathrm{m^2/s}$.

A diffusion process has one intrinsic scale: in time $t$, the quantity spreads over a distance of order $\sqrt{\alpha t}$; equivalently, crossing a distance $L$ takes a time of order $L^2/\alpha$. Doubling the thickness of a heat spreader quadruples its thermal settling time.

### Solution on a finite interval with Dirichlet conditions

Section 13.3 gave the product solutions $e^{-\alpha\lambda_nt}\sin(n\pi x/L)$ of the model problem. Superposing them,

$$
u(x,t) = \sum_{n=1}^\infty b_n\,e^{-\alpha(n\pi/L)^2t}\sin\frac{n\pi x}{L},
$$

and imposing $u(x,0) = f(x)$ gives $f(x) = \sum b_n\sin(n\pi x/L)$: the $b_n$ are the coefficients of the **Fourier sine series** of $f$ on $(0,L)$ (the half-range expansion of [Ch. 15](../part-5-transforms/15-fourier-series.md)). Multiplying by $\sin(m\pi x/L)$ and integrating, orthogonality kills every term except $n = m$, and $\int_0^L\sin^2(m\pi x/L)\,dx = L/2$.

> [!theorem] Heat equation on $(0,L)$, zero Dirichlet ends
> If $f$ is piecewise smooth on $[0,L]$, the IBVP $u_t = \alpha u_{xx}$, $u(0,t) = u(L,t) = 0$, $u(x,0) = f(x)$ has the unique solution
> $$
> u(x,t) = \sum_{n=1}^{\infty} b_n\,e^{-\alpha n^2\pi^2t/L^2}\,\sin\frac{n\pi x}{L},\qquad b_n = \frac{2}{L}\int_0^L f(x)\sin\frac{n\pi x}{L}\,dx .
> $$
> For $t>0$ the series converges, together with all its derivatives, uniformly in $x$, and $u$ is infinitely differentiable.

The convergence statement for $t>0$ is easy: $\lvert b_n\rvert\le\frac2L\int_0^L\lvert f\rvert\,dx = M$, so the $n$th term of the series for any derivative $\partial^p_x\partial^q_t u$ is bounded by $M\,(n\pi/L)^{p+2q}\alpha^q e^{-\alpha n^2\pi^2t/L^2}$, and the exponential beats any power: the Weierstrass M-test applies. At $t = 0$ the series is just the Fourier sine series of $f$, whose convergence is governed by the theory of [Ch. 15](../part-5-transforms/15-fourier-series.md) (pointwise to the mid-value at jumps, with Gibbs overshoot).

### Properties of heat-equation solutions

1. **Exponential decay and the dominant mode.** For large $t$ the $n = 1$ term dominates:
   $$
   u(x,t)\approx b_1e^{-t/\tau_1}\sin\frac{\pi x}{L},\qquad \tau_1 = \frac{L^2}{\pi^2\alpha}.
   $$
   The ratio of the second surviving term to the first decays like $e^{-(n_2^2-1)t/\tau_1}$, so after about one time constant the profile *is* a half sine wave. $\tau_1$ is the **thermal time constant** of the rod.
2. **Smoothing.** Mode $n$ is damped by $e^{-n^2t/\tau_1}$. High spatial frequencies (sharp corners, jumps) die almost instantly; a discontinuous initial temperature becomes smooth for every $t>0$. The heat operator is a spatial *low-pass filter* whose bandwidth shrinks as $1/\sqrt{t}$.
3. **Maximum principle.**

> [!theorem] Maximum principle for the heat equation
> If $u_t = \alpha u_{xx}$ in $0<x<L$, $0<t\le T$, then the maximum (and minimum) of $u$ over the closed rectangle $[0,L]\times[0,T]$ is attained either at $t = 0$ or on the ends $x = 0$, $x = L$.

   *Physical meaning:* without sources, a point in the interior cannot become hotter than the hottest initial or boundary temperature. *Proof idea:* if $u$ had an interior maximum at $(x_0,t_0)$, then $u_{xx}\le0$ and $u_t\ge0$ there; perturbing to $v = u + \varepsilon x^2$ makes the inequality $v_t - \alpha v_{xx} = -2\alpha\varepsilon<0$ strict, contradicting $v_t\ge0\ge v_{xx}$ at its maximum. The principle gives uniqueness and continuous dependence on the data at once: two solutions whose data differ by at most $\varepsilon$ differ everywhere by at most $\varepsilon$.

4. **Irreversibility.** Running the series backwards, $e^{+\alpha n^2\pi^2\lvert t\rvert/L^2}$ multiplies the $n$th coefficient. An error of $10^{-6}$ in $b_{10}$ grows by $e^{100\,\lvert t\rvert/\tau_1}$; at $\lvert t\rvert = \tau_1$ that is $e^{100}\approx 3\times10^{43}$. You cannot reconstruct an earlier temperature profile from a later one.

> [!example] Example 13.12 — Cooling of a uniformly heated rod
> A rod of length $L = 1$ (nondimensional units, $\alpha = 1$) starts at a uniform $u = 100$ and has both ends plunged into ice at $0$. Find $u(x,t)$ and the midpoint temperature at $t = 0.05,\ 0.1,\ 0.2$.
>
> **Coefficients.**
> $$
> b_n = 2\int_0^1 100\sin(n\pi x)\,dx = \frac{200}{n\pi}\big(1 - \cos n\pi\big) = \begin{cases}\dfrac{400}{n\pi}, & n\text{ odd},\\[1ex] 0, & n\text{ even}.\end{cases}
> $$
> **Solution.**
> $$
> u(x,t) = \frac{400}{\pi}\sum_{n\ \text{odd}}\frac{1}{n}\,e^{-n^2\pi^2t}\sin(n\pi x).
> $$
> **Midpoint.** At $x = \tfrac12$, $\sin(n\pi/2) = (-1)^{(n-1)/2}$ for odd $n$. Summing 1000 terms:
>
> | $t$ | $u(\tfrac12,t)$ (series) | one-term $\tfrac{400}{\pi}e^{-\pi^2t}$ |
> |---|---|---|
> | 0.01 | 99.92 | 115.36 |
> | 0.05 | 77.23 | 77.73 |
> | 0.1 | 47.45 | 47.45 |
> | 0.2 | 17.69 | 17.69 |
> | 0.5 | 0.916 | 0.916 |
>
> For $t\ge0.1$ the one-term approximation is accurate to four figures; for very small $t$ it is useless (it even exceeds the initial temperature, violating the maximum principle) and many terms are needed. The time constant is $\tau_1 = 1/\pi^2 = 0.101$.

```desmos-graph
left=0; right=1; bottom=-5; top=120
---
y=\sum_{n=1}^{40}\frac{400}{(2n-1)\pi}e^{-(2n-1)^{2}\pi^{2}\cdot0.002}\sin((2n-1)\pi x)|red|0<x<1
y=\sum_{n=1}^{40}\frac{400}{(2n-1)\pi}e^{-(2n-1)^{2}\pi^{2}\cdot0.01}\sin((2n-1)\pi x)|orange|0<x<1
y=\sum_{n=1}^{40}\frac{400}{(2n-1)\pi}e^{-(2n-1)^{2}\pi^{2}\cdot0.05}\sin((2n-1)\pi x)|green|0<x<1
y=\sum_{n=1}^{40}\frac{400}{(2n-1)\pi}e^{-(2n-1)^{2}\pi^{2}\cdot0.1}\sin((2n-1)\pi x)|blue|0<x<1
y=\sum_{n=1}^{40}\frac{400}{(2n-1)\pi}e^{-(2n-1)^{2}\pi^{2}\cdot0.2}\sin((2n-1)\pi x)|purple|0<x<1
```

*Snapshots of Example 13.12 at $t = 0.002$ (red), $0.01$, $0.05$, $0.1$ and $0.2$ (purple): the flat-topped profile first rounds off near the cold ends, then relaxes to a decaying half sine wave.*

> [!example] Example 13.13 — Thermal time constant of a copper bus bar
> A copper bar of length $L = 0.10\ \mathrm{m}$ ($\alpha = 1.11\times10^{-4}\ \mathrm{m^2/s}$), insulated along its sides, is at $100\ ^\circ\mathrm{C}$ above ambient when both ends are suddenly clamped to ambient. Estimate its time constant and the time for the centre excess temperature to fall to $1\ \mathrm{K}$.
>
> $\tau_1 = L^2/(\pi^2\alpha) = 0.01/(9.870\times1.11\times10^{-4}) = 9.13\ \mathrm{s}$.
>
> From Example 13.12 (rescaled: $x/L$ and $t/(L^2/\alpha)$ are the nondimensional variables), the centre excess temperature for $t\gtrsim\tau_1$ is $\frac{400}{\pi}e^{-t/\tau_1}\ \mathrm{K}$. Setting it equal to $1$:
> $$
> t = \tau_1\ln\frac{400}{\pi} = 9.13\times4.847 = 44.2\ \mathrm{s}.
> $$
> For a $1\ \mathrm{m}$ bar the time constant is $100$ times longer, $15\ \mathrm{min}$ — the $L^2$ law in action.

### Insulated ends: the Neumann problem

If both ends are insulated, $u_x(0,t) = u_x(L,t) = 0$, the eigenfunctions are $\cos(n\pi x/L)$, $n = 0,1,2,\dots$ (Section 13.3), and

$$
u(x,t) = \frac{a_0}{2} + \sum_{n=1}^\infty a_n\,e^{-\alpha n^2\pi^2t/L^2}\cos\frac{n\pi x}{L},\qquad a_n = \frac{2}{L}\int_0^L f(x)\cos\frac{n\pi x}{L}\,dx,
$$

a Fourier **cosine** series. As $t\to\infty$, $u\to a_0/2 = \frac1L\int_0^Lf\,dx$, the **average** initial temperature: no heat escapes, so the total heat content $\int_0^L\rho c_pu\,dx$ is conserved and simply redistributes.

> [!example] Example 13.14 — Equalisation in an insulated rod
> An insulated rod with $L = 1$, $\alpha = 1$ starts with the linear profile $u(x,0) = x$. Find $u(x,t)$ and the end temperatures at $t = 0.1$.
>
> $a_0 = 2\int_0^1x\,dx = 1$. For $n\ge1$, integrating by parts,
> $$
> a_n = 2\int_0^1 x\cos(n\pi x)\,dx = 2\Big[\frac{x\sin n\pi x}{n\pi}\Big]_0^1 - \frac{2}{n\pi}\int_0^1\sin n\pi x\,dx = \frac{2(\cos n\pi - 1)}{n^2\pi^2},
> $$
> which is $-4/(n^2\pi^2)$ for odd $n$ and $0$ for even $n$. So
> $$
> u(x,t) = \frac12 - \frac{4}{\pi^2}\sum_{n\ \text{odd}}\frac{1}{n^2}e^{-n^2\pi^2t}\cos n\pi x .
> $$
> **Checks.** At $t = 0$, $x = 0$: $\tfrac12 - \tfrac{4}{\pi^2}\sum_{\text{odd}}1/n^2 = \tfrac12 - \tfrac4{\pi^2}\cdot\tfrac{\pi^2}{8} = 0$. ✓ At $t = 0.1$: $u(0,0.1) = 0.3489$, $u(1,0.1) = 0.6511$; they sum to $1$ by the antisymmetry of $\cos n\pi x$ ($n$ odd) about $x = \tfrac12$. Both ends are converging to the average $\tfrac12$ with time constant $1/\pi^2$.
>
> The coefficients decay like $1/n^2$ (compared with $1/n$ in Example 13.12) because the even periodic extension of $f(x) = x$ is continuous — a general rule from [Ch. 15](../part-5-transforms/15-fourier-series.md): smoother extension, faster coefficient decay.

### Non-homogeneous boundary conditions

> [!example] Example 13.15 — A rod heated at one end
> A rod with $L = 1$, $\alpha = 1$, initially at $0$, has $u(0,t) = 0$ and $u(1,t) = 100$ for $t>0$. Find $u(x,t)$ and $u(\tfrac12,t)$.
>
> **Steady state.** $u_s = 100x$.
>
> **Transient.** $w = u - 100x$ satisfies the model problem with $w(x,0) = -100x$. Its sine coefficients:
> $$
> b_n = 2\int_0^1(-100x)\sin n\pi x\,dx = -200\Big[-\frac{x\cos n\pi x}{n\pi}\Big]_0^1 - 200\int_0^1\frac{\cos n\pi x}{n\pi}dx = \frac{200(-1)^n}{n\pi}.
> $$
> **Solution.**
> $$
> u(x,t) = 100x + \frac{200}{\pi}\sum_{n=1}^\infty\frac{(-1)^n}{n}e^{-n^2\pi^2t}\sin n\pi x .
> $$
> **Midpoint values** (4000 terms): $u(\tfrac12,0) = 0.008$ (the truncated series of a function that is $0$ there; the error is the slowly converging $1/n$ tail), $u(\tfrac12,0.05) = 11.38$, $u(\tfrac12,0.1) = 26.28$, $u(\tfrac12,0.2) = 41.16$, approaching the steady value $50$. Only odd $n$ contribute at $x = \tfrac12$ and for $t\ge0.1$ the single term $50 - \frac{200}{\pi}e^{-\pi^2t}$ is accurate to better than $0.1$: at $t = 0.2$ it gives $50 - 63.66\times0.1389 = 41.16$. ✓

### Infinite and semi-infinite domains

When the rod is so long that its ends never "feel" the disturbance during the time of interest (a thick substrate under a hot spot, a semiconductor wafer seen from a surface a few micrometres deep), we replace the interval by the whole line or half-line. Fourier *series* become Fourier *integrals* ([Ch. 16](../part-5-transforms/16-fourier-transforms.md)); the result can be written in closed form.

> [!theorem] Heat kernel and solution on the whole line
> The function
> $$
> \Phi(x,t) = \frac{1}{\sqrt{4\pi\alpha t}}\,e^{-x^2/(4\alpha t)},\qquad t>0,
> $$
> solves $u_t = \alpha u_{xx}$, has $\int_{-\infty}^\infty\Phi\,dx = 1$ for every $t>0$, and tends to $\delta(x)$ as $t\to0^+$. For bounded continuous $f$, the solution of $u_t = \alpha u_{xx}$ on $-\infty<x<\infty$ with $u(x,0) = f(x)$ is the convolution
> $$
> u(x,t) = \int_{-\infty}^{\infty}\Phi(x-s,t)\,f(s)\,ds .
> $$

*Derivation sketch.* The equation is unchanged by the scaling $x\to\kappa x$, $t\to\kappa^2t$, which suggests looking for a solution of the form $t^{-1/2}F(\eta)$ with **similarity variable** $\eta = x/\sqrt{4\alpha t}$. Substituting gives the ODE $F'' + 2\eta F' + 2F = 0$, i.e. $(F' + 2\eta F)' = 0$, whose decaying solution is $F = Ce^{-\eta^2}$; $C$ is fixed by unit area. Alternatively, Fourier-transform in $x$: $\hat u_t = -\alpha\omega^2\hat u$, so $\hat u = \hat f\,e^{-\alpha\omega^2t}$, and the inverse transform of the Gaussian $e^{-\alpha\omega^2t}$ is $\Phi$ ([Ch. 16](../part-5-transforms/16-fourier-transforms.md)). In signal terms, the heat equation is a linear time-varying *Gaussian filter* whose impulse response has standard deviation $\sigma = \sqrt{2\alpha t}$.

```desmos-graph
left=-6; right=6; bottom=-0.05; top=1.2
---
y=\frac{1}{\sqrt{4\pi\cdot0.0625}}e^{-\frac{x^{2}}{4\cdot0.0625}}|red
y=\frac{1}{\sqrt{4\pi\cdot0.25}}e^{-\frac{x^{2}}{4\cdot0.25}}|orange
y=\frac{1}{\sqrt{4\pi\cdot1}}e^{-\frac{x^{2}}{4\cdot1}}|green
y=\frac{1}{\sqrt{4\pi\cdot4}}e^{-\frac{x^{2}}{4\cdot4}}|blue
```

*The heat kernel ($\alpha = 1$) at $t = 1/16,\ 1/4,\ 1,\ 4$: each quadrupling of time doubles the width and halves the peak, keeping unit area.*

For a semi-infinite body $x>0$ initially at $0$ whose surface is suddenly raised to $T_s$ at $t = 0$, the similarity solution is built from the **error function** $\operatorname{erf}z = \frac{2}{\sqrt\pi}\int_0^ze^{-s^2}ds$ and its complement $\operatorname{erfc}z = 1-\operatorname{erf}z$:

$$
u(x,t) = T_s\,\operatorname{erfc}\frac{x}{2\sqrt{\alpha t}} .
$$

Check: $\operatorname{erfc}(0) = 1$ gives $u(0,t) = T_s$; as $t\to0^+$ for fixed $x>0$ the argument $\to\infty$ and $u\to0$; and $u$ is a function of the similarity variable, so it solves the PDE by the derivation above (with $F$ replaced by its integral). The surface heat flux $-k u_x(0,t) = kT_s/\sqrt{\pi\alpha t}$ is infinite at $t = 0^+$ and decays like $t^{-1/2}$.

> [!example] Example 13.16 — Step response of a long RC interconnect
> A long on-chip line with $RC = 4\times10^{-5}\ \mathrm{s/m^2}$ (Example 13.4) is initially discharged, and at $t = 0$ its input is stepped to $V_0$. Treat it as semi-infinite. When does the voltage at $x = 1\ \mathrm{mm}$ reach $V_0/2$?
>
> The RC-line equation $v_t = v_{xx}/(RC)$ is the heat equation with $\alpha = 1/(RC)$, so
> $$
> v(x,t) = V_0\,\operatorname{erfc}\Big(\frac{x}{2}\sqrt{\frac{RC}{t}}\Big).
> $$
> Set $\operatorname{erfc}(z) = 0.5$: $z = 0.4769$. Then $t_{50} = RC\,x^2/(4z^2) = 1.099\,RC\,x^2$. With $x = 10^{-3}\ \mathrm{m}$: $t_{50} = 1.099\times4\times10^{-5}\times10^{-6} = 4.4\times10^{-11}\ \mathrm{s} = 44\ \mathrm{ps}$.
>
> The delay grows as $x^2$: repeaters (buffers) are inserted along long on-chip wires precisely because breaking a line of length $\ell$ into $k$ segments cuts the diffusive delay from $\propto\ell^2$ to $\propto k(\ell/k)^2 = \ell^2/k$. For a *finite* open-ended line of length $\ell$, the series solution (Section 13.7) gives the classic $t_{50}\approx0.38\,RC\ell^2$; the semi-infinite estimate is pessimistic because there is no reflection from the open end to speed up charging.

> [!warning] Common pitfall — one-term approximations at early times
> The dominant-mode approximation $b_1e^{-t/\tau_1}\sin(\pi x/L)$ is excellent once $t\gtrsim0.2\tau_1$ to $\tau_1$ but grossly wrong at early times (Example 13.12 shows it exceeding the initial temperature). At early times use the semi-infinite erfc solution near each boundary instead: the two approximations are complementary, short-time and long-time.

> [!tip] Engineering note — thermal RC networks
> Device datasheets give a transient thermal impedance $Z_{th}(t)$ as a sum of exponentials, the "Foster network" $Z_{th}(t) = \sum_i R_i(1 - e^{-t/\tau_i})$. That form is exactly a truncated eigenfunction expansion of the heat equation in the package: each $\tau_i$ is a modal time constant of the kind $L^2/(\mu_i^2\alpha)$, and the $R_i$ are the modal weights. The $t^{1/2}$ behaviour of $Z_{th}$ at very short times is the semi-infinite erfc regime.

---

## 13.5 Wave equation

### Derivations

**Lossless transmission line.** With $R = G = 0$ the telegrapher's equations are $v_x = -L\,i_t$ and $i_x = -C\,v_t$, and eliminating $i$ gives

$$
v_{tt} = c^2v_{xx},\qquad c = \frac{1}{\sqrt{LC}} .
$$

For a coaxial line with dielectric permittivity $\varepsilon$ and permeability $\mu$, $LC = \mu\varepsilon$, so $c = 1/\sqrt{\mu\varepsilon}$ — the speed of light in the dielectric, independent of the geometry.

**Taut string.** A string under tension $T_s$ (N) with mass per unit length $\rho_\ell$ (kg/m), displaced by a small amount $u(x,t)$, feels a net transverse force $T_s[u_x(x+\Delta x) - u_x(x)]$ on a piece of length $\Delta x$ (to first order in the slope). Newton's law $\rho_\ell\Delta x\,u_{tt} = T_s\Delta x\,u_{xx}$ gives $u_{tt} = c^2u_{xx}$ with $c = \sqrt{T_s/\rho_\ell}$.

**Electromagnetic plane waves.** In a source-free, lossless medium Maxwell's equations give $\nabla^2\mathbf{E} = \mu\varepsilon\,\mathbf{E}_{tt}$; for a field $E_y(x,t)$ that varies only with $x$, this is again the one-dimensional wave equation.

> [!definition] Wave equation
> $$
> \frac{\partial^2u}{\partial t^2} = c^2\,\frac{\partial^2u}{\partial x^2}\quad\text{(1-D)},\qquad u_{tt} = c^2\nabla^2u\quad\text{(2-D, 3-D)},
> $$
> where $c>0$ is the **wave speed**. It is hyperbolic with characteristics $x\pm ct = \text{const}$.

### d'Alembert's solution on the whole line

> [!theorem] General solution of the 1-D wave equation
> Every twice-differentiable solution of $u_{tt} = c^2u_{xx}$ on $-\infty<x<\infty$ has the form
> $$
> u(x,t) = F(x - ct) + G(x + ct)
> $$
> for some functions $F$ and $G$: the sum of a wave moving right at speed $c$ and a wave moving left at speed $c$, each keeping its shape.

*Proof.* Change to characteristic variables $\xi = x - ct$, $\eta = x + ct$. Then $\partial_x = \partial_\xi + \partial_\eta$ and $\partial_t = c(\partial_\eta - \partial_\xi)$, so

$$
u_{tt} - c^2u_{xx} = c^2\big[(\partial_\eta - \partial_\xi)^2 - (\partial_\xi + \partial_\eta)^2\big]u = -4c^2\,u_{\xi\eta}.
$$

The PDE becomes $u_{\xi\eta} = 0$, whose general solution (Example 13.1) is $F(\xi) + G(\eta)$. $\square$

A profile $F(x - ct)$ is constant along each line $x - ct = x_0$: the value that sat at $x_0$ at time $0$ sits at $x_0 + ct$ at time $t$. Now impose the initial conditions $u(x,0) = f(x)$, $u_t(x,0) = g(x)$:

$$
F(x) + G(x) = f(x),\qquad -cF'(x) + cG'(x) = g(x).
$$

Integrate the second: $-F(x) + G(x) = \frac1c\int_0^xg(s)\,ds + K$. Adding and subtracting,

$$
G(x) = \frac12f(x) + \frac1{2c}\int_0^xg\,ds + \frac K2,\qquad F(x) = \frac12f(x) - \frac1{2c}\int_0^xg\,ds - \frac K2 .
$$

Evaluate $F$ at $x - ct$ and $G$ at $x + ct$ and add; the constant $K$ cancels:

> [!theorem] d'Alembert's formula
> The solution of $u_{tt} = c^2u_{xx}$, $-\infty<x<\infty$, $u(x,0) = f(x)$, $u_t(x,0) = g(x)$ is
> $$
> u(x,t) = \frac{f(x-ct) + f(x+ct)}{2} + \frac{1}{2c}\int_{x-ct}^{x+ct}g(s)\,ds .
> $$
> If $f\in C^2$ and $g\in C^1$ this is a classical solution, and the problem is well posed: $\lvert u_1 - u_2\rvert\le\max\lvert f_1 - f_2\rvert + t\max\lvert g_1 - g_2\rvert$.

Three consequences follow from the formula.

- **Domain of dependence.** $u(x_0,t_0)$ depends only on the data on the interval $[x_0 - ct_0,\ x_0 + ct_0]$, the base of the backward characteristic triangle. Nothing outside can influence it.
- **Range of influence.** Data at a point $x_0$ influence only the forward wedge $\lvert x - x_0\rvert\le ct$. Signals travel at exactly speed $c$ — in contrast to the heat equation, whose kernel is non-zero everywhere for every $t>0$.
- **No smoothing.** A kink or jump in $f$ is carried unchanged along the characteristics. The wave equation transmits information; it does not forget it.

> [!example] Example 13.17 — A Gaussian pulse splits in two
> A string (or line) with $c = 1$ starts at rest with $f(x) = e^{-x^2}$, $g = 0$. By d'Alembert,
> $$
> u(x,t) = \tfrac12\big[e^{-(x-t)^2} + e^{-(x+t)^2}\big].
> $$
> At $t = 0$ the peak is $1$ at $x = 0$. By $t = 3$ the two halves are centred at $x = \pm3$, each of height $\tfrac12(1 + e^{-36})\approx0.5$; the overlap is negligible once $t\gtrsim2$. At $x = 0$: $u(0,t) = e^{-t^2}$, e.g. $u(0,1) = 0.368$.
>
> On a transmission line, the same statement reads: a voltage pulse with no accompanying current ($i(x,0) = 0$, so $v_t(x,0) = 0$) splits into equal forward and backward pulses of half the amplitude.

```desmos-graph
left=-6; right=6; bottom=-0.1; top=1.1
---
y=\frac{1}{2}\left(e^{-x^{2}}+e^{-x^{2}}\right)|red
y=\frac{1}{2}\left(e^{-(x-0.7)^{2}}+e^{-(x+0.7)^{2}}\right)|orange
y=\frac{1}{2}\left(e^{-(x-1.5)^{2}}+e^{-(x+1.5)^{2}}\right)|green
y=\frac{1}{2}\left(e^{-(x-3)^{2}}+e^{-(x+3)^{2}}\right)|blue
y=\frac{1}{2}\left(e^{-(x-4.5)^{2}}+e^{-(x+4.5)^{2}}\right)|purple
```

*Snapshots of Example 13.17 at $t = 0,\ 0.7,\ 1.5,\ 3,\ 4.5$: the initial pulse separates into two half-height copies travelling in opposite directions without change of shape.*

> [!example] Example 13.18 — Response to an initial velocity
> A string with $c = 2$ is initially straight, $f = 0$, and is struck so that $g(x) = \dfrac{1}{1+x^2}$. Find $u(x,t)$ and $u(0,t)$.
>
> $$
> u(x,t) = \frac{1}{2c}\int_{x-ct}^{x+ct}\frac{ds}{1+s^2} = \frac14\big[\arctan(x+2t) - \arctan(x-2t)\big].
> $$
> At $x = 0$: $u(0,t) = \tfrac14[\arctan2t + \arctan2t] = \tfrac12\arctan(2t)$. So $u(0,0.5) = \tfrac12\arctan1 = 0.3927$, $u(0,1) = \tfrac12\arctan2 = 0.5536$, and $u(0,t)\to\pi/4 = 0.785$ as $t\to\infty$.
>
> The string does **not** return to zero: an initial impulse of total momentum $\rho_\ell\int g\,dx$ produces a permanent displacement in the region the waves have passed. The total area $\int g\,dx = \pi$ is spread over a plateau of width $2ct$ and height $\frac{1}{2c}\int g = \pi/4$. On a line, the analogue is an initial current distribution leaving a residual voltage step.

### Travelling waves on a transmission line

For the lossless line, d'Alembert's form applies to $v$, and the first telegrapher equation then determines $i$. Let $v = F(x - ct) + G(x + ct)$. From $i_x = -Cv_t = -C\big[-cF'(x-ct) + cG'(x+ct)\big]$, integrating in $x$:

$$
i(x,t) = Cc\big[F(x-ct) - G(x+ct)\big] = \frac{1}{Z_0}\big[F(x-ct) - G(x+ct)\big],\qquad Z_0 = \frac{1}{Cc} = \sqrt{\frac{L}{C}} .
$$

> [!theorem] Forward and backward waves on a lossless line
> Write $v^+ = F(x-ct)$ and $v^- = G(x+ct)$. Then
> $$
> v = v^+ + v^-,\qquad i = \frac{v^+ - v^-}{Z_0},
> $$
> where $Z_0 = \sqrt{L/C}$ is the **characteristic impedance**. A single travelling wave has $v/i = +Z_0$ (forward) or $-Z_0$ (backward); a line of finite length looks resistive, of value $Z_0$, to a signal until the first reflection returns.

**Reflection at a resistive load.** At $x = \ell$ put a resistor $R_L$, so $v(\ell,t) = R_Li(\ell,t)$:

$$
v^+ + v^- = \frac{R_L}{Z_0}(v^+ - v^-)\quad\Longrightarrow\quad v^- = \Gamma\,v^+,\qquad \Gamma = \frac{R_L - Z_0}{R_L + Z_0}.
$$

This is the solution of the dynamic boundary condition found in Example 13.7. Special cases: matched load $R_L = Z_0$, $\Gamma = 0$ (no reflection — the line looks infinite); open circuit $\Gamma = +1$; short circuit $\Gamma = -1$.

> [!example] Example 13.19 — Step on a mismatched line
> A $Z_0 = 50\ \Omega$ line of length $\ell = 2\ \mathrm{m}$ with $c = 2\times10^8\ \mathrm{m/s}$ is driven by a matched source (internal resistance $50\ \Omega$) whose open-circuit voltage steps from $0$ to $10\ \mathrm{V}$ at $t = 0$, and is terminated in $R_L = 150\ \Omega$. Describe $v(x,t)$.
>
> 1. **Launch.** Initially the line looks like $Z_0$, so the source sees a $50\ \Omega$ divider: $v^+ = 10\times\frac{50}{50+50} = 5\ \mathrm{V}$, and a $5\ \mathrm{V}$ step travels toward the load, arriving at $t_1 = \ell/c = 10\ \mathrm{ns}$.
> 2. **Reflection.** $\Gamma_L = (150-50)/(150+50) = 0.5$, so $v^- = 2.5\ \mathrm{V}$ travels back. The load voltage jumps to $5 + 2.5 = 7.5\ \mathrm{V}$ at $t = 10\ \mathrm{ns}$.
> 3. **Back at the source** at $t = 20\ \mathrm{ns}$: the source is matched, $\Gamma_S = 0$, so nothing more is reflected and the line voltage is $7.5\ \mathrm{V}$ everywhere.
>
> **Check against DC.** The final state is a $10\ \mathrm{V}$ source driving $50\ \Omega + 150\ \Omega$: $v = 10\times150/200 = 7.5\ \mathrm{V}$. ✓ Current check at the load: $i = (5 - 2.5)/50 = 0.05\ \mathrm{A}$ and $7.5/150 = 0.05\ \mathrm{A}$. ✓
>
> With a mismatched source as well, the waves bounce repeatedly, multiplying by $\Gamma_L\Gamma_S$ on each round trip, and the line voltage approaches the DC value as a geometric series — the "bounce diagram" of signal-integrity analysis.

### Finite string: the method of images

On a semi-infinite string $x>0$ with a fixed end, $u(0,t) = 0$, extend the data to $x<0$ **oddly**: $f(-x) = -f(x)$, $g(-x) = -g(x)$. The d'Alembert solution of the extended problem is odd in $x$ for all $t$, so it vanishes at $x = 0$ automatically. A pulse approaching the fixed end meets its inverted mirror image coming the other way; they cancel at the end, and the image emerges as an **inverted reflected pulse** ($\Gamma = -1$, the short circuit). For a free end, $u_x(0,t) = 0$, extend **evenly**; the reflection is upright ($\Gamma = +1$, the open circuit). On $0<x<L$ with both ends fixed, extend oddly about both ends, i.e. to an odd $2L$-periodic function, and apply d'Alembert to the extension.

### Finite interval by separation of variables

For $u_{tt} = c^2u_{xx}$ on $0<x<L$ with $u(0,t) = u(L,t) = 0$, $u(x,0) = f(x)$, $u_t(x,0) = g(x)$, separation gives $X'' + \lambda X = 0$ (Dirichlet: $\lambda_n = (n\pi/L)^2$) and $T'' + c^2\lambda_nT = 0$, so $T_n = A_n\cos\omega_nt + B_n\sin\omega_nt$ with

$$
\omega_n = \frac{n\pi c}{L},\qquad f_n = \frac{\omega_n}{2\pi} = \frac{nc}{2L}.
$$

> [!theorem] Vibrating string / shorted line on $(0,L)$
> $$
> u(x,t) = \sum_{n=1}^\infty\Big(A_n\cos\frac{n\pi ct}{L} + B_n\sin\frac{n\pi ct}{L}\Big)\sin\frac{n\pi x}{L},
> $$
> $$
> A_n = \frac2L\int_0^Lf(x)\sin\frac{n\pi x}{L}\,dx,\qquad B_n = \frac{2}{n\pi c}\int_0^Lg(x)\sin\frac{n\pi x}{L}\,dx .
> $$
> The $n$th term is the $n$th **normal mode** (standing wave), with natural frequency $f_n = nc/(2L)$, an integer multiple of the **fundamental** $f_1 = c/(2L)$. The motion is periodic with period $2L/c$.

The $B_n$ formula comes from differentiating in $t$ at $t = 0$: $g = \sum B_n\omega_n\sin(n\pi x/L)$, so $B_n\omega_n$ is the sine coefficient of $g$. Using $\sin a\cos b = \tfrac12[\sin(a-b) + \sin(a+b)]$, each standing wave is the sum of two travelling waves, and the series reproduces d'Alembert's formula applied to the odd periodic extension.

> [!example] Example 13.20 — The plucked string
> A string of length $L = 1$ with $c = 1$ is pulled up at its midpoint to height $h$ and released from rest: $f(x) = 2hx$ for $0\le x\le\tfrac12$, $f(x) = 2h(1-x)$ for $\tfrac12\le x\le1$; $g = 0$.
>
> **Coefficients.** $B_n = 0$. By symmetry about $x = \tfrac12$, $A_n = 0$ for even $n$. For any $n$, integrating by parts twice (or using the formula for the sine coefficients of a triangle),
> $$
> A_n = 2\int_0^1 f(x)\sin n\pi x\,dx = \frac{8h}{n^2\pi^2}\sin\frac{n\pi}{2}.
> $$
> So $A_1 = 0.8106h$, $A_3 = -0.0901h$, $A_5 = 0.0324h$:
> $$
> u(x,t) = \frac{8h}{\pi^2}\Big[\sin\pi x\cos\pi t - \frac19\sin3\pi x\cos3\pi t + \frac1{25}\sin5\pi x\cos5\pi t - \cdots\Big].
> $$
> **Energy per mode.** The energy of mode $n$ is proportional to $\omega_n^2A_n^2\propto n^2\cdot n^{-4} = n^{-2}$. The fraction in the fundamental is $1\big/\sum_{n\ \text{odd}}n^{-2} = 8/\pi^2 = 81.1\%$. Plucking at the centre suppresses all even harmonics; plucking nearer an end strengthens the higher harmonics, which is why a guitar sounds brighter near the bridge.
>
> **Shape.** Summing 400 terms, at $t = \tfrac14$ the profile is a trapezoid with a flat top at $h/2$ between $x = 0.25$ and $0.75$; at $t = \tfrac12$ the string is momentarily straight ($u\equiv0$); at $t = 1$ it is the inverted triangle. The corners (kinks) travel at speed $c$ and reflect from the ends — d'Alembert's picture.

```desmos-graph
left=0; right=1; bottom=-1.1; top=1.1
---
y=\sum_{n=1}^{60}\frac{8}{n^{2}\pi^{2}}\sin\left(\frac{n\pi}{2}\right)\sin(n\pi x)|red|0<x<1
y=\sum_{n=1}^{60}\frac{8}{n^{2}\pi^{2}}\sin\left(\frac{n\pi}{2}\right)\sin(n\pi x)\cos(0.1n\pi)|orange|0<x<1
y=\sum_{n=1}^{60}\frac{8}{n^{2}\pi^{2}}\sin\left(\frac{n\pi}{2}\right)\sin(n\pi x)\cos(0.25n\pi)|green|0<x<1
y=\sum_{n=1}^{60}\frac{8}{n^{2}\pi^{2}}\sin\left(\frac{n\pi}{2}\right)\sin(n\pi x)\cos(0.4n\pi)|blue|0<x<1
y=\sum_{n=1}^{60}\frac{8}{n^{2}\pi^{2}}\sin\left(\frac{n\pi}{2}\right)\sin(n\pi x)\cos(0.75n\pi)|purple|0<x<1
```

*The plucked string of Example 13.20 ($h = 1$) at $t = 0,\ 0.1,\ 0.25,\ 0.4,\ 0.75$: the initial triangle splits into two kinks that run outward, leaving a flat plateau that sinks through zero at $t = 0.5$ and reappears inverted.*

> [!example] Example 13.21 — Resonances of a shorted coaxial stub
> A coaxial line of length $\ell = 0.75\ \mathrm{m}$ with $c = 2\times10^8\ \mathrm{m/s}$ is shorted at both ends (a closed resonator). Its voltage modes are $\sin(n\pi x/\ell)$ with
> $$
> f_n = \frac{nc}{2\ell} = n\times\frac{2\times10^8}{1.5}\ \mathrm{Hz} = 133.3\,n\ \mathrm{MHz}.
> $$
> If instead one end is open (D–N problem, Section 13.3), the eigenfunctions are $\sin\big((2n-1)\pi x/(2\ell)\big)$ and $f_n = (2n-1)c/(4\ell) = 66.7,\ 200,\ 333\ \mathrm{MHz},\dots$ — only odd multiples of the quarter-wave frequency. Only the odd multiples of the quarter-wave frequency resonate — the characteristic signature of a quarter-wave resonator, used in stub filters and in clock-distribution problems.

### Energy conservation

> [!theorem] Energy of the wave equation
> For $u_{tt} = c^2u_{xx}$ on $0<x<L$ with fixed or free ends, the energy
> $$
> E(t) = \frac12\int_0^L\big(u_t^2 + c^2u_x^2\big)\,dx
> $$
> is constant in time.

*Proof.* $E'(t) = \int_0^L(u_tu_{tt} + c^2u_xu_{xt})\,dx = \int_0^L u_t\,(u_{tt} - c^2u_{xx})\,dx + c^2\big[u_xu_t\big]_0^L = 0$, because the PDE kills the integral and $u_t = 0$ (fixed end) or $u_x = 0$ (free end) kills the boundary term. $\square$

For a transmission line the corresponding quantity is $\int_0^\ell\big(\tfrac12Cv^2 + \tfrac12Li^2\big)dx$, the electric plus magnetic stored energy; a resistive load makes the boundary term $-R_Li^2\le0$, the power it absorbs. Energy conservation gives uniqueness (as in Section 13.2) and, with Parseval's theorem ([Ch. 15](../part-5-transforms/15-fourier-series.md)), the energy-per-mode decomposition used in Example 13.20.

### The lossy line: distortion and the Heaviside condition

With losses, a sinusoidal wave $e^{j(\omega t - \beta x)}$ on the telegrapher's equation has a complex propagation constant $\gamma = \sqrt{(R + j\omega L)(G + j\omega C)}$, and in general both its attenuation and its phase velocity depend on frequency: a pulse spreads (disperses) as well as shrinks. There is one exception.

> [!theorem] Distortionless line (Heaviside)
> Put $v = e^{-at}w$ with $a = \tfrac12(R/L + G/C)$. Then $w$ satisfies
> $$
> LC\,w_{tt} - w_{xx} = \tfrac{LC}{4}\Big(\frac RL - \frac GC\Big)^2 w .
> $$
> If $R/L = G/C$ the right side vanishes, $w$ obeys the lossless wave equation, and
> $$
> v(x,t) = e^{-(R/L)t}\big[F(x - ct) + G(x + ct)\big],\qquad c = 1/\sqrt{LC} :
> $$
> pulses travel at speed $c$ without change of shape, merely attenuated.

*Derivation.* With $v = e^{-at}w$: $v_t = e^{-at}(w_t - aw)$, $v_{tt} = e^{-at}(w_{tt} - 2aw_t + a^2w)$. Substituting into $v_{xx} = LCv_{tt} + (RC+GL)v_t + RGv$ and dividing by $e^{-at}$:

$$
w_{xx} = LC\,w_{tt} + \big(RC + GL - 2aLC\big)w_t + \big(LCa^2 - (RC+GL)a + RG\big)w .
$$

The choice $a = (RC+GL)/(2LC) = \tfrac12(R/L + G/C)$ removes the $w_t$ term, and the $w$ coefficient becomes $RG - LCa^2 = -\tfrac{LC}{4}(R/L - G/C)^2$. $\square$

For a forward wave at position $x$ the amplitude factor $e^{-(R/L)t}$, evaluated at the arrival time $t = x/c$, is $e^{-(R/L)\sqrt{LC}\,x} = e^{-(R/Z_0)x}$: the attenuation constant is $R/Z_0 = \sqrt{RG}$ nepers per metre.

> [!example] Example 13.22 — Loading a line to make it distortionless
> A telephone pair has $L = 0.6\ \mathrm{mH/km}$, $C = 50\ \mathrm{nF/km}$, $R = 50\ \Omega/\mathrm{km}$, $G = 1\ \mu\mathrm{S/km}$. Is it distortionless, and what series inductance would make it so?
>
> $R/L = 50/(0.6\times10^{-3}) = 8.3\times10^4\ \mathrm{s^{-1}}$; $G/C = 10^{-6}/(50\times10^{-9}) = 20\ \mathrm{s^{-1}}$. Far from equal: $R/L\gg G/C$, the typical situation. Since adding $G$ wastes power, Heaviside's remedy was to *increase* $L$: $L_{\text{new}} = RC/G = 50\times50\times10^{-9}/10^{-6} = 2.5\ \mathrm{H/km}$ — impractically large, but even partial loading (Pupin coils of tens of mH every couple of km) dramatically reduced distortion on early long-distance telephone lines.
>
> With full loading, $Z_0 = \sqrt{L/C} = \sqrt{2.5/(50\times10^{-9})} = 7071\ \Omega$ and the attenuation is $\sqrt{RG} = \sqrt{50\times10^{-6}} = 7.07\times10^{-3}\ \mathrm{Np/km}$ ($0.061\ \mathrm{dB/km}$).

> [!warning] Common pitfall — applying d'Alembert directly on a finite interval
> d'Alembert's formula needs $f$ and $g$ on the whole interval $[x - ct, x + ct]$. On a bounded interval, once $ct$ exceeds the distance to an end, that interval sticks out of the domain: you must first extend $f$ and $g$ (oddly for fixed ends, evenly for free ends, periodically with period $2L$) and only then apply the formula. Forgetting this ignores the reflections.

> [!tip] Engineering note — time-domain reflectometry
> A TDR instrument launches a fast step into a cable and watches the reflected wave. By d'Alembert, a reflection seen at time $t$ comes from a discontinuity at distance $ct/2$, and its height and sign give the local $\Gamma$ and hence the impedance there. Faults (opens, shorts, crushed sections, bad connectors) on cables and PCB traces are located this way to centimetre precision. The method relies on exactly the two properties proved above: finite propagation speed and shape preservation.

---

## 13.6 Laplace equation

### Where it comes from

In electrostatics $\mathbf{E} = -\nabla V$ and Gauss's law $\nabla\cdot(\varepsilon\mathbf{E}) = \rho_v$ give, for uniform $\varepsilon$,

$$
\nabla^2V = -\frac{\rho_v}{\varepsilon}\quad\text{(Poisson)},\qquad \nabla^2V = 0\quad\text{in charge-free regions (Laplace)}.
$$

The same equation describes steady-state temperature ($u_t = 0$ in the heat equation), steady current flow in a resistive sheet or volume ($\nabla\cdot\mathbf{J} = 0$, $\mathbf{J} = -\sigma\nabla V$), magnetostatic scalar potential in current-free regions, and the steady state of any diffusion process. Laplace's equation is the **equilibrium** equation: it is what the heat and wave equations settle to once all transients have died.

> [!definition] Laplace's equation and harmonic functions
> $$
> \nabla^2u = 0,\qquad\text{in 2-D Cartesian } u_{xx} + u_{yy} = 0,\qquad\text{in polar } u_{rr} + \frac1ru_r + \frac{1}{r^2}u_{\phi\phi} = 0 .
> $$
> A twice continuously differentiable solution on a region $\Omega$ is called **harmonic** in $\Omega$.

Examples of harmonic functions: $x$, $xy$, $x^2 - y^2$, $e^{kx}\sin ky$, $\ln r$, $r^n\cos n\phi$, $r^{-n}\sin n\phi$. In two dimensions, the real and imaginary parts of any analytic function of $z = x + jy$ are harmonic ([Ch. 18](../part-5-transforms/18-complex-analysis.md)); e.g. $z^2 = (x^2 - y^2) + j\,2xy$.

### Properties of harmonic functions

> [!theorem] Mean-value property
> If $u$ is harmonic in a region containing the closed disk of radius $\rho$ centred at $\mathbf{x}_0$, then $u(\mathbf{x}_0)$ equals the average of $u$ over the circle of radius $\rho$:
> $$
> u(\mathbf{x}_0) = \frac{1}{2\pi}\int_0^{2\pi}u(\mathbf{x}_0 + \rho\,\hat{\mathbf r}(\phi))\,d\phi .
> $$
> (In 3-D, the average over a sphere.)

*Proof sketch.* Let $M(\rho)$ be the circle average. Differentiating under the integral, $M'(\rho) = \frac{1}{2\pi}\int_0^{2\pi}\partial u/\partial r\,d\phi = \frac{1}{2\pi\rho}\oint\partial u/\partial n\,ds$, which by the divergence theorem equals $\frac{1}{2\pi\rho}\iint\nabla^2u\,dA = 0$. So $M$ is constant, and $M(0^+) = u(\mathbf{x}_0)$. $\square$

> [!theorem] Maximum principle and uniqueness
> A non-constant harmonic function on a bounded connected region attains its maximum and minimum only on the boundary. Consequently the Dirichlet problem $\nabla^2u = 0$ in $\Omega$, $u = h$ on $\partial\Omega$, has at most one solution, and if $\lvert h_1 - h_2\rvert\le\varepsilon$ on the boundary then $\lvert u_1 - u_2\rvert\le\varepsilon$ throughout.

The maximum principle follows from the mean-value property: a value cannot be an average of its neighbours and strictly larger than all of them. Physically, **there is no electrostatic potential maximum in free space**, so a charge cannot be held in stable equilibrium by electrostatic forces alone (Earnshaw's theorem), and the hottest point of an unheated plate is on its edge.

### Rectangular domains

Consider the rectangle $0<x<a$, $0<y<b$ with Dirichlet data. By superposition, split a problem with data on all four sides into four problems, each with data on one side and zero on the other three. Take the case of data on the top side:

$$
u_{xx} + u_{yy} = 0,\qquad u(0,y) = u(a,y) = 0,\quad u(x,0) = 0,\quad u(x,b) = h(x).
$$

Separate, $u = X(x)Y(y)$: $X''/X = -Y''/Y = -\lambda$. The homogeneous conditions at $x = 0, a$ give the Dirichlet eigenproblem, $\lambda_n = (n\pi/a)^2$, $X_n = \sin(n\pi x/a)$. Then $Y'' - (n\pi/a)^2Y = 0$ with $Y(0) = 0$, so $Y_n = \sinh(n\pi y/a)$. (In the direction that carries the non-homogeneous data the solutions are exponential, not oscillatory — the hallmark of Laplace's equation.) Superposing,

$$
u(x,y) = \sum_{n=1}^\infty c_n\sinh\frac{n\pi y}{a}\sin\frac{n\pi x}{a},
$$

and the top condition $u(x,b) = h(x)$ makes $c_n\sinh(n\pi b/a)$ the sine coefficients of $h$.

> [!theorem] Dirichlet problem on a rectangle, data on one side
> $$
> u(x,y) = \sum_{n=1}^\infty h_n\,\frac{\sinh(n\pi y/a)}{\sinh(n\pi b/a)}\,\sin\frac{n\pi x}{a},\qquad h_n = \frac2a\int_0^ah(x)\sin\frac{n\pi x}{a}\,dx .
> $$
> The ratio $\sinh(n\pi y/a)/\sinh(n\pi b/a)\approx e^{-n\pi(b-y)/a}$ for large $n$: the $n$th harmonic of the boundary data decays exponentially with distance from the boundary, over a length $a/(n\pi)$.

That last remark is fundamental. Fine detail in the boundary potential (period $2a/n$) is invisible a short distance away; only the smoothest components penetrate. This is why the field far from a periodic electrode array (an interdigitated capacitor, a grating) is nearly uniform, and why a grounded shield with small holes still shields well.

> [!example] Example 13.23 — A single-harmonic boundary
> On the square $0<x,y<1$, $u = 0$ on three sides and $u(x,1) = V_0\sin\pi x$. Find $u$ at the centre.
>
> Only $h_1 = V_0$ is non-zero, so the series has one term:
> $$
> u(x,y) = V_0\,\frac{\sinh\pi y}{\sinh\pi}\sin\pi x .
> $$
> At the centre: $u(\tfrac12,\tfrac12) = V_0\sinh(\pi/2)/\sinh\pi = V_0\times2.3013/11.5487 = 0.1993\,V_0$.
>
> **Check** harmonicity directly: $u_{xx} = -\pi^2u$, $u_{yy} = +\pi^2u$, sum $0$. ✓ Check the maximum principle: the maximum $V_0$ is on the boundary at $(\tfrac12,1)$. ✓

> [!example] Example 13.24 — Trough with a raised lid
> A long rectangular conducting trough of width $a$ and height $b$ has its three sides grounded and is closed by a lid, insulated from the sides, at potential $V_0$. Find $V(x,y)$ and evaluate it for a square trough.
>
> With $h = V_0$: $h_n = \frac{2V_0}{a}\int_0^a\sin\frac{n\pi x}{a}dx = \frac{2V_0}{n\pi}(1 - \cos n\pi)$, i.e. $4V_0/(n\pi)$ for odd $n$, $0$ for even $n$:
> $$
> V(x,y) = \frac{4V_0}{\pi}\sum_{n\ \text{odd}}\frac1n\,\frac{\sinh(n\pi y/a)}{\sinh(n\pi b/a)}\sin\frac{n\pi x}{a}.
> $$
> **Square trough ($a = b = 1$), numerical values:**
>
> | Point $(x,y)$ | $V/V_0$ |
> |---|---|
> | $(0.5,\ 0.25)$ | $0.0954$ |
> | $(0.25,\ 0.5)$ | $0.1820$ |
> | $(0.5,\ 0.5)$ | $0.2500$ |
> | $(0.5,\ 0.75)$ | $0.5405$ |
>
> **Exact check at the centre.** Rotate the square by $90°$ four times: each rotation gives a problem with the lid on a different side. Their sum has $V_0$ on all four sides, whose solution is the constant $V_0$. By symmetry all four contribute equally at the centre, so $V(\tfrac12,\tfrac12) = V_0/4$ exactly. The series gives $0.25000$; its first term alone gives $0.2537$, already within $1.5\%$, because $\sinh(n\pi/2)/\sinh(n\pi)\approx e^{-n\pi/2}$ decays fast.
>
> Near the top corners the series converges slowly, because the boundary data jump from $V_0$ to $0$ there (the Gibbs phenomenon of [Ch. 15](../part-5-transforms/15-fourier-series.md) is present on the line $y = b$ itself, but disappears for $y<b$).

```desmos-graph
left=0; right=1; bottom=-0.05; top=1.1
---
y=\sum_{n=1}^{25}\frac{4}{(2n-1)\pi}\frac{\sinh((2n-1)\pi\cdot0.25)}{\sinh((2n-1)\pi)}\sin((2n-1)\pi x)|blue|0<x<1
y=\sum_{n=1}^{25}\frac{4}{(2n-1)\pi}\frac{\sinh((2n-1)\pi\cdot0.5)}{\sinh((2n-1)\pi)}\sin((2n-1)\pi x)|green|0<x<1
y=\sum_{n=1}^{25}\frac{4}{(2n-1)\pi}\frac{\sinh((2n-1)\pi\cdot0.75)}{\sinh((2n-1)\pi)}\sin((2n-1)\pi x)|orange|0<x<1
y=\sum_{n=1}^{25}\frac{4}{(2n-1)\pi}\frac{\sinh((2n-1)\pi\cdot0.95)}{\sinh((2n-1)\pi)}\sin((2n-1)\pi x)|red|0<x<1
y=1|black|dashed
```

*Potential profiles $V(x,y_0)/V_0$ across the square trough of Example 13.24 at heights $y_0 = 0.25,\ 0.5,\ 0.75,\ 0.95$ (blue to red). Close to the lid the profile approaches the flat boundary value with sharp edges; deeper down only the smooth fundamental survives.*

> [!example] Example 13.25 — Data on two sides; a semi-infinite strip
> (a) A square plate $0<x,y<1$ has $u(x,1) = 100$ and $u(1,y) = 100$, the other two sides at $0$. By superposition this is the sum of Example 13.24 (lid on top) and its reflection in the diagonal $y = x$ (lid on the right). At the centre, $u = 25 + 25 = 50$; at $(0.5, 0.25)$, $u = 100(0.0954) + 100\,V_{\text{lid}}(0.25, 0.5) = 9.54 + 18.20 = 27.74$, using the table (the reflection maps $(0.5,0.25)$ to $(0.25,0.5)$).
>
> (b) **Semi-infinite strip.** $0<x<a$, $y>0$, sides $x = 0, a$ grounded, bottom $u(x,0) = h(x)$, $u$ bounded as $y\to\infty$. Now $Y'' = (n\pi/a)^2Y$ must stay bounded, so $Y_n = e^{-n\pi y/a}$ and
> $$
> u(x,y) = \sum_{n=1}^\infty h_n\,e^{-n\pi y/a}\sin\frac{n\pi x}{a}.
> $$
> With $h = V_0$ and $a = 1$, at $y = 1$ the first term is $\frac{4V_0}{\pi}e^{-\pi} = 0.055\,V_0$ and the next is $\frac{4V_0}{3\pi}e^{-3\pi} = 3.4\times10^{-5}V_0$: one width above the electrode only the fundamental is left. This is the basis of the rule of thumb that a grounded slot or waveguide below cut-off attenuates fields by about $27\ \mathrm{dB}$ ($e^{-\pi}$) per width of length.

### Polar domains: disks, exteriors and annuli

From Example 13.11, the product solutions of Laplace's equation in polar coordinates are

$$
1,\quad\ln r,\quad r^n\cos n\phi,\quad r^n\sin n\phi,\quad r^{-n}\cos n\phi,\quad r^{-n}\sin n\phi\qquad(n = 1,2,\dots).
$$

Which ones are allowed depends on the region:

| Region | Allowed terms | Reason |
|---|---|---|
| Interior of disk $r<a$ | $1$, $r^n\cos n\phi$, $r^n\sin n\phi$ | $\ln r$ and $r^{-n}$ blow up at $r = 0$ |
| Exterior $r>a$, bounded at $\infty$ | $1$, $r^{-n}\cos n\phi$, $r^{-n}\sin n\phi$ | $r^n$ and $\ln r$ grow without bound |
| Annulus $a<r<b$ | all | $r = 0$ and $\infty$ excluded |

> [!theorem] Dirichlet problem for a disk
> The solution of $\nabla^2u = 0$ in $r<a$ with $u(a,\phi) = h(\phi)$ is
> $$
> u(r,\phi) = \frac{a_0}{2} + \sum_{n=1}^\infty\Big(\frac ra\Big)^n\big(a_n\cos n\phi + b_n\sin n\phi\big),
> $$
> where $a_n$, $b_n$ are the ordinary Fourier coefficients of $h$ on $(-\pi,\pi)$:
> $a_n = \frac1\pi\int_{-\pi}^{\pi}h\cos n\phi\,d\phi$, $b_n = \frac1\pi\int_{-\pi}^{\pi}h\sin n\phi\,d\phi$.
> Summing the series in closed form gives **Poisson's integral formula**
> $$
> u(r,\phi) = \frac{1}{2\pi}\int_{-\pi}^{\pi}\frac{a^2 - r^2}{a^2 - 2ar\cos(\phi - \psi) + r^2}\,h(\psi)\,d\psi .
> $$

*Derivation of Poisson's formula.* Substitute the coefficient integrals into the series and exchange sum and integral (legitimate for $r<a$ because of the geometric factor $(r/a)^n$):

$$
u = \frac{1}{2\pi}\int_{-\pi}^{\pi}h(\psi)\Big[1 + 2\sum_{n=1}^\infty\Big(\frac ra\Big)^n\cos n(\phi - \psi)\Big]d\psi .
$$

With $\rho = r/a$ and $\theta = \phi - \psi$, the bracket is $\operatorname{Re}\big\{1 + 2\sum_{n\ge1}(\rho e^{j\theta})^n\big\} = \operatorname{Re}\frac{1 + \rho e^{j\theta}}{1 - \rho e^{j\theta}} = \frac{1 - \rho^2}{1 - 2\rho\cos\theta + \rho^2}$, a geometric series. $\square$

At $r = 0$ Poisson's formula gives $u(0) = \frac{1}{2\pi}\int h\,d\psi$ — the mean-value property again. The **Poisson kernel** is positive, has unit integral and concentrates at $\psi = \phi$ as $r\to a$: it plays for the disk the role the heat kernel plays for the line.

> [!example] Example 13.26 — A split cylinder
> A long conducting cylinder of radius $a$ is split lengthwise into two halves separated by thin insulating gaps; the upper half ($0<\phi<\pi$) is held at $V_0$, the lower half at $0$. Find the potential inside.
>
> **Coefficients.** $a_0 = \frac1\pi\int_0^\pi V_0\,d\phi = V_0$; $a_n = \frac{V_0}{\pi}\int_0^\pi\cos n\phi\,d\phi = 0$; $b_n = \frac{V_0}{\pi}\int_0^\pi\sin n\phi\,d\phi = \frac{V_0}{n\pi}(1 - \cos n\pi)$, i.e. $2V_0/(n\pi)$ for odd $n$:
> $$
> V(r,\phi) = \frac{V_0}{2} + \frac{2V_0}{\pi}\sum_{n\ \text{odd}}\frac1n\Big(\frac ra\Big)^n\sin n\phi .
> $$
> **Closed form.** Using $\sum_{n\ \text{odd}}\frac{z^n}{n} = \tfrac12\ln\frac{1+z}{1-z}$ with $z = (r/a)e^{j\phi}$ and taking imaginary parts,
> $$
> V(r,\phi) = \frac{V_0}{2} + \frac{V_0}{\pi}\arctan\frac{2ar\sin\phi}{a^2 - r^2}.
> $$
> **Numbers.** At $r = a/2$, $\phi = \pi/2$: the arctan argument is $1/(3/4) = 4/3$, so $V = V_0[0.5 + \arctan(4/3)/\pi] = V_0(0.5 + 0.2952) = 0.7952\,V_0$; the series (200 terms) gives the same value. At the centre $V = V_0/2$, the average of the boundary, as the mean-value property demands. The equipotentials are arcs of circles through the two gaps.

> [!example] Example 13.27 — Coaxial line: the annulus
> Find the potential between the conductors of a coaxial line, inner radius $a$ at $V_0$ and outer radius $b$ grounded, and its capacitance per unit length.
>
> The data do not depend on $\phi$, so only the $n = 0$ terms survive: $V = A + B\ln r$. $V(b) = 0$ and $V(a) = V_0$ give
> $$
> V(r) = V_0\,\frac{\ln(b/r)}{\ln(b/a)},\qquad E_r = -\frac{dV}{dr} = \frac{V_0}{r\ln(b/a)} .
> $$
> The surface charge on the inner conductor is $\varepsilon E_r(a)$; per unit length the charge is $Q' = 2\pi a\,\varepsilon E_r(a) = 2\pi\varepsilon V_0/\ln(b/a)$, so
> $$
> C' = \frac{Q'}{V_0} = \frac{2\pi\varepsilon}{\ln(b/a)} .
> $$
> **Numbers** (RG-58-like): $a = 0.45\ \mathrm{mm}$, $b = 1.47\ \mathrm{mm}$, $\varepsilon_r = 2.25$. $\ln(b/a) = 1.1838$, $C' = 2\pi\times8.854\times10^{-12}\times2.25/1.1838 = 105.7\ \mathrm{pF/m}$. With $L' = \frac{\mu_0}{2\pi}\ln(b/a) = 236.8\ \mathrm{nH/m}$, the characteristic impedance is $Z_0 = \sqrt{L'/C'} = 47.3\ \Omega$ and $c = 1/\sqrt{L'C'} = 2.0\times10^8\ \mathrm{m/s}$, closing the loop with the transmission-line equations of Section 13.5. At $V_0 = 100\ \mathrm{V}$ the peak field (at $r = a$) is $1.88\times10^5\ \mathrm{V/m}$.

> [!example] Example 13.28 — Conducting cylinder in a uniform field (exterior problem)
> A grounded conducting cylinder of radius $a$ is placed in a uniform applied field $\mathbf{E}_0 = E_0\hat{\mathbf x}$, whose potential is $-E_0x = -E_0r\cos\phi$. Find the potential outside.
>
> We need $V(a,\phi) = 0$ and $V\to-E_0r\cos\phi$ as $r\to\infty$. The applied potential is the $n = 1$ cosine term $r\cos\phi$; add the decaying partner $r^{-1}\cos\phi$:
> $$
> V = -E_0\Big(r - \frac{a^2}{r}\Big)\cos\phi .
> $$
> It is harmonic (sum of two product solutions), vanishes at $r = a$, and approaches the applied potential far away. The surface field is $E_r(a,\phi) = -\partial V/\partial r\big|_a = 2E_0\cos\phi$: the cylinder **doubles** the field at the points facing the applied field. (For a sphere the factor is 3.) Field enhancement at curved conductors is why high-voltage hardware uses large radii and why sharp points emit corona.

> [!warning] Common pitfall — forgetting the $\ln r$ and $r^{-n}$ terms in an annulus
> In an annulus or an exterior region with non-zero net charge, the term $\ln r$ (and for the annulus also $r^{-n}$) is essential; dropping it makes the coax problem unsolvable. Conversely, in a full disk including the origin those terms must be dropped. The rule is not "always use $r^n$"; it is "keep exactly those terms that are finite everywhere in the region".

> [!tip] Engineering note — capacitance and resistance extraction
> Every capacitance or resistance calculation for an arbitrary geometry is a Laplace problem: solve $\nabla^2V = 0$ with the electrodes at fixed potentials, compute the charge (or current) by integrating the normal field over an electrode, and divide by the voltage. Closed forms exist for coax, parallel plates, the split cylinder, and a few conformally mappable shapes; everything else (PCB traces, IC interconnect, connector pins) is done by field solvers that discretise exactly this equation.

---

## 13.7 Fourier-series solutions

Sections 13.4–13.6 all ended with the same step: expand a given function (initial temperature, initial displacement, boundary potential) in the eigenfunctions of a Sturm–Liouville problem. This section collects the general framework, explains what convergence results we can rely on, and extends the method to non-homogeneous equations (sources and forcing) and to two space dimensions. The theory of Fourier series itself — coefficient formulas, even and odd extensions, Dirichlet conditions, the Gibbs phenomenon, Parseval's theorem — is developed in [Ch. 15](../part-5-transforms/15-fourier-series.md); here we use its results.

### The general recipe

For a linear evolution problem $u_t = \mathcal{A}u$ (or $u_{tt} = \mathcal{A}u$) with $\mathcal{A} = \alpha\,\partial_x^2$ (or $c^2\partial_x^2$) and homogeneous boundary conditions, let $\{X_n\}$ be the eigenfunctions of $-\partial_x^2$ with those boundary conditions, $X_n'' = -\lambda_nX_n$. Expand the unknown at every instant:

$$
u(x,t) = \sum_n T_n(t)\,X_n(x),\qquad T_n(t) = \frac{\langle u(\cdot,t), X_n\rangle}{\langle X_n, X_n\rangle} .
$$

Because $\partial_x^2$ acting on $X_n$ just multiplies by $-\lambda_n$, the PDE becomes an *uncoupled* set of ODEs, one per mode:

| PDE | Modal ODE | Modal solution |
|---|---|---|
| $u_t = \alpha u_{xx}$ | $T_n' = -\alpha\lambda_nT_n$ | $T_n(0)\,e^{-\alpha\lambda_nt}$ |
| $u_{tt} = c^2u_{xx}$ | $T_n'' = -c^2\lambda_nT_n$ | $T_n(0)\cos\omega_nt + \dfrac{T_n'(0)}{\omega_n}\sin\omega_nt$, $\omega_n = c\sqrt{\lambda_n}$ |
| $u_{xx} + u_{yy} = 0$ ($y$ plays "time") | $Y_n'' = \lambda_nY_n$ | $A_n\cosh\sqrt{\lambda_n}y + B_n\sinh\sqrt{\lambda_n}y$ |

The initial data supply $T_n(0)$, $T_n'(0)$ through the **Fourier coefficients** of $f$ and $g$. In matrix language, the eigenfunctions diagonalise the spatial operator, exactly as eigenvectors diagonalise $\mathbf{A}$ in $\dot{\mathbf{x}} = \mathbf{A}\mathbf{x}$.

### Which Fourier series?

The boundary conditions decide which expansion is natural on $(0,L)$:

| Boundary conditions | Expansion of $f$ on $(0,L)$ | Equivalent periodic extension |
|---|---|---|
| D–D | sine series $\sum b_n\sin(n\pi x/L)$ | odd, period $2L$ |
| N–N | cosine series $\frac{a_0}{2} + \sum a_n\cos(n\pi x/L)$ | even, period $2L$ |
| D–N | quarter-wave sine series $\sum b_n\sin\frac{(2n-1)\pi x}{2L}$ | odd about $0$, even about $L$, period $4L$ |
| periodic on $(-\pi,\pi)$ | full series $\frac{a_0}{2} + \sum(a_n\cos n\phi + b_n\sin n\phi)$ | $2\pi$-periodic |
| Robin | $\sum c_n X_n$ with $c_n = \int fX_n/\int X_n^2$ | none simple; use Sturm–Liouville orthogonality |

The middle column is what you compute; the right-hand column tells you, via [Ch. 15](../part-5-transforms/15-fourier-series.md), how fast the coefficients decay and where Gibbs overshoot will appear. For example, a constant initial temperature $f = 100$ with Dirichlet ends has an odd extension that is a square wave: coefficients $\propto1/n$, Gibbs ringing at the ends at $t = 0$. With Neumann ends the even extension of the same constant is constant: the expansion has a single term.

### Convergence: what we may rely on

> [!theorem] Convergence of series solutions
> Let $f$ be piecewise smooth on $[0,L]$.
> 1. **At $t = 0$** the series equals the Fourier series of the appropriate extension of $f$; it converges to $f(x)$ where the extension is continuous and to the average of the one-sided limits at jumps ([Ch. 15](../part-5-transforms/15-fourier-series.md)). Convergence is uniform if and only if the extension is continuous (no jumps, including at the ends).
> 2. **Heat equation, $t>0$:** the series and all its term-by-term derivatives converge uniformly; $u$ is $C^\infty$ and satisfies the PDE classically. Moreover $u(\cdot,t)\to f$ in mean square as $t\to0^+$.
> 3. **Wave equation:** no smoothing. The series is a classical solution if the odd extensions of $f$ and $g$ are $C^2$ and $C^1$ respectively (equivalently $\sum n^2\lvert A_n\rvert$ and $\sum n\lvert B_n\rvert$ converge); otherwise it is a **weak solution** whose kinks or jumps travel along characteristics.
> 4. **Laplace equation:** inside the region, at positive distance from the boundary that carries the data, the factors $e^{-n\pi d/a}$ (rectangle) or $(r/a)^n$ (disk) make the series converge uniformly with all derivatives, and $u$ is $C^\infty$ — even real-analytic — however rough the data.

The smoothness of $f$ controls the decay of its coefficients: jumps give $1/n$, kinks (continuous $f$, jump in $f'$) give $1/n^2$, and so on, each extra derivative of continuity buying one more power of $n$ ([Ch. 15](../part-5-transforms/15-fourier-series.md)). For the heat equation this matters only at small $t$; for the wave equation it matters for all time.

### How many terms?

Because the modal factors decay like $e^{-\alpha\lambda_nt}$, the number of terms needed for a given accuracy is set by $t$.

> [!example] Example 13.29 — Truncation error of the cooling-rod series
> In Example 13.12 ($u = \frac{400}{\pi}\sum_{n\text{ odd}}\frac1ne^{-n^2\pi^2t}\sin n\pi x$), how many terms are needed for an error below $0.01$ at $t = 0.01$ and at $t = 0.001$?
>
> **Bound.** $\lvert\sin\rvert\le1$, so the error after keeping all odd $n<N$ is at most the tail $\displaystyle\varepsilon_N(t) = \frac{400}{\pi}\sum_{n\ge N,\ \text{odd}}\frac1ne^{-n^2\pi^2t}$. Since the terms decrease faster than geometrically, the first omitted term is a good estimate.
>
> | $t$ | $N = 5$ | $N = 11$ | $N = 21$ | $N = 41$ |
> |---|---|---|---|---|
> | $0.01$ | $2.31$ | $7.6\times10^{-5}$ | $\approx10^{-18}$ | — |
> | $0.001$ | $44.5$ | $7.02$ | $0.124$ | $2.4\times10^{-7}$ |
>
> At $t = 0.01$ the odd terms $n = 1,\dots,9$ (five terms) suffice; at $t = 0.001$ about $n\le 29$ is needed (between $N = 21$ and $N = 41$). A rule of thumb: keep terms until $n^2\pi^2\alpha t/L^2\gtrsim 10$, i.e. $n\gtrsim\frac{L}{\pi}\sqrt{10/(\alpha t)}$. For very small $t$, switch to the erfc approximation of Section 13.4.

### Non-homogeneous equations: eigenfunction expansion

Sources and forcing are handled by the same expansion. Consider

$$
u_t = \alpha u_{xx} + s(x,t),\qquad u(0,t) = u(L,t) = 0,\qquad u(x,0) = f(x).
$$

Expand $u = \sum T_n(t)\sin(n\pi x/L)$ and $s = \sum s_n(t)\sin(n\pi x/L)$ with $s_n(t) = \frac2L\int_0^Ls(x,t)\sin(n\pi x/L)\,dx$. Substituting and matching coefficients (orthogonality),

$$
T_n' + \alpha\lambda_nT_n = s_n(t),\qquad T_n(0) = b_n,
$$

a first-order linear ODE for each mode ([Ch. 11](11-ordinary-differential-equations.md)) with solution

$$
T_n(t) = b_n e^{-\alpha\lambda_nt} + \int_0^t e^{-\alpha\lambda_n(t-\tau)}s_n(\tau)\,d\tau .
$$

The second term is a **convolution** of the modal source with the modal impulse response $e^{-\alpha\lambda_nt}$ — **Duhamel's principle**. Each mode is a first-order low-pass filter with cut-off $\alpha\lambda_n$, driven by the projection of the source onto that mode. The same idea for the wave equation gives a bank of undamped resonators $T_n'' + \omega_n^2T_n = s_n(t)$, and a source whose time dependence contains $\omega_n$ excites mode $n$ resonantly.

> [!example] Example 13.30 — Joule heating of a wire with cold ends
> A wire of length $L = 1$ (nondimensional, $\alpha = 1$) starts at $u = 0$, its ends are held at $0$, and a uniform current switched on at $t = 0$ produces a constant source $s = 1$. Find $u(x,t)$ and the midpoint temperature.
>
> **Modal source.** $s_n = 2\int_0^1\sin n\pi x\,dx = 4/(n\pi)$ for odd $n$, $0$ for even $n$. With $b_n = 0$ and constant $s_n$:
> $$
> T_n(t) = \frac{s_n}{\lambda_n}\big(1 - e^{-\lambda_nt}\big) = \frac{4}{n^3\pi^3}\big(1 - e^{-n^2\pi^2t}\big),\qquad n\ \text{odd}.
> $$
> $$
> u(x,t) = \frac{4}{\pi^3}\sum_{n\ \text{odd}}\frac{1 - e^{-n^2\pi^2t}}{n^3}\sin n\pi x .
> $$
> **Steady state check.** As $t\to\infty$, $u\to u_s$ with $u_s'' = -1$, $u_s(0) = u_s(1) = 0$, i.e. $u_s = x(1-x)/2$, maximum $0.125$ at the centre. The series gives $\frac{4}{\pi^3}\sum_{\text{odd}}\frac{(-1)^{(n-1)/2}}{n^3} = \frac{4}{\pi^3}\cdot\frac{\pi^3}{32} = 0.125$. ✓ (The one-term value is $4/\pi^3 = 0.1290$: the $1/n^3$ decay makes even the steady-state series converge fast.)
>
> **Transient.** $u(\tfrac12,0.05) = 0.0463$, $u(\tfrac12,0.1) = 0.0769$, $u(\tfrac12,1) = 0.12499$. The approach to steady state is governed by the slowest mode, time constant $1/\pi^2$.
>
> In dimensional form with $s = g/(\rho c_p)$, $g = J^2/\sigma$, the steady centre temperature rise is $gL^2/(8k)$: halving the length of a fuse element quarters its temperature rise, a key design relation for fuses and thin-film resistors.

> [!example] Example 13.31 — Finite open-ended RC line (the 0.38 RCℓ² rule)
> An RC line of length $\ell$, open at $x = \ell$, initially discharged, has its input stepped to $V_0$ at $t = 0$. Find $v(\ell,t)$ and the 50% delay.
>
> **Set-up.** $v_t = v_{xx}/(RC)$, $v(0,t) = V_0$ (Dirichlet), $v_x(\ell,t) = 0$ (open circuit, Section 13.2), $v(x,0) = 0$. Steady state $v_s = V_0$; transient $w = v - V_0$ with D–N eigenfunctions $\sin\mu_nx$, $\mu_n = (2n-1)\pi/(2\ell)$, and $w(x,0) = -V_0$:
> $$
> b_n = \frac{2}{\ell}\int_0^\ell(-V_0)\sin\mu_nx\,dx = -\frac{2V_0}{\ell\mu_n} = -\frac{4V_0}{(2n-1)\pi}.
> $$
> $$
> v(x,t) = V_0\Big[1 - \frac4\pi\sum_{n=1}^\infty\frac{\sin\mu_nx}{2n-1}\,e^{-\mu_n^2t/(RC)}\Big].
> $$
> **At the open end**, $\sin\mu_n\ell = (-1)^{n-1}$. In terms of $\tau = t/(RC\ell^2)$:
> $$
> \frac{v(\ell,t)}{V_0} = 1 - \frac4\pi\Big[e^{-\pi^2\tau/4} - \frac13e^{-9\pi^2\tau/4} + \frac15e^{-25\pi^2\tau/4} - \cdots\Big].
> $$
> Solving $v(\ell,t) = V_0/2$ numerically gives $\tau_{50} = 0.3787$, so $t_{50} = 0.38\,RC\ell^2$ — the standard figure used in interconnect timing. The one-term formula $\tau_{50} = \frac{4}{\pi^2}\ln\frac8\pi = 0.3788$ is already accurate to four figures. The 90% point is at $\tau = 1.03$.
>
> **Comparison.** The semi-infinite estimate of Example 13.16 gave $1.10\,RC\ell^2$; the finite open line is about three times faster because charge stops flowing out at the open end and piles up there.

```desmos-graph
left=0; right=1; bottom=-0.05; top=1.1
---
y=1-\frac{4}{\pi}\sum_{n=1}^{40}\frac{\sin\left(\frac{(2n-1)\pi x}{2}\right)}{2n-1}e^{-(2n-1)^{2}\pi^{2}\cdot\frac{0.01}{4}}|red|0<x<1
y=1-\frac{4}{\pi}\sum_{n=1}^{40}\frac{\sin\left(\frac{(2n-1)\pi x}{2}\right)}{2n-1}e^{-(2n-1)^{2}\pi^{2}\cdot\frac{0.05}{4}}|orange|0<x<1
y=1-\frac{4}{\pi}\sum_{n=1}^{40}\frac{\sin\left(\frac{(2n-1)\pi x}{2}\right)}{2n-1}e^{-(2n-1)^{2}\pi^{2}\cdot\frac{0.2}{4}}|green|0<x<1
y=1-\frac{4}{\pi}\sum_{n=1}^{40}\frac{\sin\left(\frac{(2n-1)\pi x}{2}\right)}{2n-1}e^{-(2n-1)^{2}\pi^{2}\cdot\frac{0.38}{4}}|blue|0<x<1
y=1-\frac{4}{\pi}\sum_{n=1}^{40}\frac{\sin\left(\frac{(2n-1)\pi x}{2}\right)}{2n-1}e^{-(2n-1)^{2}\pi^{2}\cdot\frac{1}{4}}|purple|0<x<1
```

*Voltage $v(x,t)/V_0$ along the open-ended RC line of Example 13.31 ($x$ in units of $\ell$) at $\tau = t/(RC\ell^2) = 0.01,\ 0.05,\ 0.2,\ 0.38,\ 1$. The profile has zero slope at the open end; at $\tau = 0.38$ the far end reaches one half.*

### Two space dimensions: double Fourier series

On a rectangle $0<x<a$, $0<y<b$ with all edges at zero, the heat equation $u_t = \alpha(u_{xx} + u_{yy})$ separates as $X(x)Y(y)T(t)$ with two separation constants. The eigenfunctions are products,

$$
\psi_{mn}(x,y) = \sin\frac{m\pi x}{a}\sin\frac{n\pi y}{b},\qquad \lambda_{mn} = \Big(\frac{m\pi}{a}\Big)^2 + \Big(\frac{n\pi}{b}\Big)^2,
$$

and the solution is a **double Fourier sine series**

$$
u = \sum_{m=1}^\infty\sum_{n=1}^\infty B_{mn}\,e^{-\alpha\lambda_{mn}t}\,\psi_{mn},\qquad B_{mn} = \frac{4}{ab}\int_0^a\!\!\int_0^b f(x,y)\,\psi_{mn}\,dy\,dx .
$$

For the wave equation (a rectangular membrane; the same eigenvalues give the cut-off frequencies of the TM modes of a rectangular waveguide) the frequencies are $\omega_{mn} = c\sqrt{\lambda_{mn}}$.

> [!example] Example 13.32 — Cooling of a square plate; degenerate modes
> A square plate of side $a = 1$ ($\alpha = 1$) with edges at zero starts at $u = 1$.
>
> $B_{mn} = 4\Big(\int_0^1\sin m\pi x\,dx\Big)\Big(\int_0^1\sin n\pi y\,dy\Big) = \dfrac{16}{mn\pi^2}$ for $m, n$ both odd, zero otherwise. The slowest mode $(1,1)$ has $\lambda_{11} = 2\pi^2$, time constant $1/(2\pi^2) = 0.0507$ — half that of a rod of the same length, because heat now escapes in two directions.
>
> **Centre temperature for large $t$**: $u(\tfrac12,\tfrac12,t)\approx\frac{16}{\pi^2}e^{-2\pi^2t}$. Note this is exactly the *square* of the 1-D one-term result $\frac4\pi e^{-\pi^2t}$: for product initial data on a product domain, the 2-D solution is the product of 1-D solutions, $u(x,y,t) = u_1(x,t)\,u_1(y,t)$.
>
> **Degeneracy.** $\lambda_{13} = \lambda_{31} = 10\pi^2$: distinct eigenfunctions share an eigenvalue. In a square drum or a square waveguide this means two modes with the same frequency; any combination $\alpha\psi_{13} + \beta\psi_{31}$ is also a mode, and small asymmetries split the pair. Circular regions lead to Bessel functions in $r$ (a singular Sturm–Liouville problem; see the Further reading), with the same structure.

> [!warning] Common pitfall — differentiating a Fourier series term by term
> Term-by-term differentiation is safe for the heat-equation series at $t>0$ and for Laplace series away from the data boundary, but **not** at $t = 0$, and not for wave-equation series with non-smooth data. Differentiating the sine series of a square wave produces a divergent series. If you need $u_x$ at a boundary (to compute a heat flux or a surface charge), either use $t>0$ or interior points and take a limit, or differentiate the closed form.

> [!tip] Engineering note — spectral solvers
> Computing the coefficients numerically with the FFT, multiplying each by its modal factor $e^{-\alpha\lambda_nt}$ or $\cos\omega_nt$, and transforming back is a **spectral method**. For simple geometries it is far more accurate per unknown than finite differences, because the error falls faster than any power of the number of modes when the solution is smooth. Spectral methods are used in thermal modelling of layered IC stacks and in optical beam propagation, which solves a heat-like (paraxial) equation.

---

## 13.8 Applications to fields, diffusion, and heat transfer

This final section applies the methods of the chapter to three problems that every electrical engineer meets. In each case the mathematics is already in hand; the work is in setting up the PDE, its domain and its auxiliary conditions, and in reading the answer. (Deeper treatment of the physics belongs to courses on electromagnetics, device physics and thermal design.)

### Fields: the skin effect as a diffusion problem

Inside a good conductor ($\sigma\gg\omega\varepsilon$) the displacement current is negligible, so Ampère's and Faraday's laws read $\nabla\times\mathbf{H} = \sigma\mathbf{E}$ and $\nabla\times\mathbf{E} = -\mu\,\partial\mathbf{H}/\partial t$. Taking the curl of the first and using $\nabla\times\nabla\times\mathbf{H} = \nabla(\nabla\cdot\mathbf{H}) - \nabla^2\mathbf{H}$ with $\nabla\cdot\mathbf{H} = 0$ ([Ch. 7](../part-2-calculus/07-vector-calculus.md)):

$$
\nabla^2\mathbf{H} = \mu\sigma\,\frac{\partial\mathbf{H}}{\partial t}.
$$

This is a **diffusion equation** with diffusivity $1/(\mu\sigma)$ — the magnetic field (and with it the current density $\mathbf{J} = \nabla\times\mathbf{H}$) *diffuses* into a conductor rather than propagating as a wave.

**Sinusoidal steady state.** For a conductor filling $x>0$ with a tangential surface field $H_y(0,t) = H_0\cos\omega t$, use phasors ([Ch. 2](../part-1-foundations/02-complex-numbers.md)): $H_y(x,t) = \operatorname{Re}\{\tilde H(x)e^{j\omega t}\}$. The PDE becomes the ODE

$$
\tilde H'' = j\omega\mu\sigma\,\tilde H,\qquad \tilde H(0) = H_0,\qquad \tilde H\ \text{bounded as } x\to\infty .
$$

Its characteristic roots are $\pm\sqrt{j\omega\mu\sigma} = \pm(1+j)/\delta$ with

$$
\delta = \sqrt{\frac{2}{\omega\mu\sigma}},
$$

since $\big((1+j)/\delta\big)^2 = 2j/\delta^2 = j\omega\mu\sigma$. Boundedness rejects the growing root:

> [!theorem] Skin-effect field in a conducting half-space
> $$
> H_y(x,t) = H_0\,e^{-x/\delta}\cos\Big(\omega t - \frac{x}{\delta}\Big),\qquad \delta = \sqrt{\frac{2}{\omega\mu\sigma}}\ \text{(skin depth)} .
> $$
> The amplitude falls by $e^{-1}$ per skin depth, and the phase lags by one radian per skin depth. The current density has the same profile, and the conductor's AC resistance per square is the **surface resistance** $R_s = 1/(\sigma\delta)$ — as if the current flowed uniformly in a layer of thickness $\delta$.

The solution is a "critically damped" wave: a travelling wave $\cos(\omega t - x/\delta)$ whose wavelength $2\pi\delta$ equals its decay length times $2\pi$, so it is attenuated by $e^{-2\pi}\approx 0.002$ within one wavelength. Exactly the same mathematics describes the daily or annual temperature wave penetrating into the ground, or into a heat sink under a pulsed load, with $\mu\sigma$ replaced by $1/\alpha$: the "thermal skin depth" is $\sqrt{2\alpha/\omega}$.

> [!example] Example 13.33 — Skin depth in copper
> Copper: $\sigma = 5.8\times10^7\ \mathrm{S/m}$, $\mu = \mu_0 = 4\pi\times10^{-7}\ \mathrm{H/m}$.
>
> | $f$ | $\delta = \sqrt{2/(2\pi f\mu_0\sigma)}$ | $R_s = 1/(\sigma\delta)$ |
> |---|---|---|
> | $60\ \mathrm{Hz}$ | $8.53\ \mathrm{mm}$ | $2.0\ \mu\Omega$ |
> | $1\ \mathrm{MHz}$ | $66.1\ \mu\mathrm{m}$ | $0.26\ \mathrm{m}\Omega$ |
> | $1\ \mathrm{GHz}$ | $2.09\ \mu\mathrm{m}$ | $8.25\ \mathrm{m}\Omega$ |
>
> Check of units: $\omega\mu\sigma$ has units $\mathrm{s^{-1}\cdot H\,m^{-1}\cdot S\,m^{-1}} = \mathrm{m^{-2}}$ (since $\mathrm{H\cdot S} = \mathrm{s}$), so $\delta$ is in metres. ✓
>
> **Consequences.** At $60\ \mathrm{Hz}$ conductors thicker than about $2\delta\approx17\ \mathrm{mm}$ waste copper, which is why large bus bars are hollow or laminated. A round wire of radius $a = 1\ \mathrm{mm}$ at $1\ \mathrm{MHz}$ ($a\gg\delta$) carries current essentially in an annulus of thickness $\delta$, giving $R_{ac}/R_{dc}\approx a/(2\delta) = 7.6$; this motivates Litz wire in switch-mode power-supply transformers. At $1\ \mathrm{GHz}$ the current flows within a few micrometres of the surface, so surface roughness of PCB copper (comparable to $\delta$) noticeably increases loss.

```desmos-graph
left=0; right=5; bottom=-1.1; top=1.1
---
y=e^{-x}\cos(0-x)|red
y=e^{-x}\cos(0.7854-x)|orange
y=e^{-x}\cos(1.5708-x)|green
y=e^{-x}\cos(3.1416-x)|blue
y=e^{-x}\cos(4.7124-x)|purple
y=e^{-x}|grey|dashed
y=-e^{-x}|grey|dashed
```

*Skin-effect field $H_y/H_0$ versus depth $x/\delta$ at phases $\omega t = 0,\ \pi/4,\ \pi/2,\ \pi,\ 3\pi/2$ (red to purple), inside the envelope $\pm e^{-x/\delta}$ (dashed). Beyond about $3\delta$ the field is negligible at every instant.*

### Diffusion: dopant profiles in semiconductors

Impurity atoms are introduced into silicon by diffusion at high temperature. With a constant diffusion coefficient $D$, the concentration $N(x,t)$ below a planar surface obeys Fick's second law,

$$
\frac{\partial N}{\partial t} = D\,\frac{\partial^2N}{\partial x^2},\qquad x>0,
$$

on what is effectively a semi-infinite domain (diffusion lengths are micrometres, wafers are hundreds of micrometres thick). Two boundary conditions model the two standard process steps.

> [!theorem] Standard diffusion profiles
> 1. **Constant-source (predeposition):** the surface is held at the solid-solubility limit, $N(0,t) = C_s$, and $N(x,0) = 0$. Then (Section 13.4)
>    $$
>    N(x,t) = C_s\,\operatorname{erfc}\frac{x}{2\sqrt{Dt}},\qquad Q(t) = \int_0^\infty N\,dx = \frac{2C_s\sqrt{Dt}}{\sqrt\pi}.
>    $$
> 2. **Limited-source (drive-in):** a fixed dose $Q$ (atoms per unit area) starts at the surface and the surface is sealed, $N_x(0,t) = 0$ (Neumann). By the method of images (even extension) the whole-line heat kernel gives
>    $$
>    N(x,t) = \frac{Q}{\sqrt{\pi Dt}}\,e^{-x^2/(4Dt)},
>    $$
>    a half-Gaussian whose surface concentration $Q/\sqrt{\pi Dt}$ falls as the dose spreads.

The **junction depth** $x_j$ is where the diffused concentration equals the opposite-type background doping $N_B$.

> [!example] Example 13.34 — Junction depths
> Boron diffuses into an $n$-type wafer with background $N_B = 10^{16}\ \mathrm{cm^{-3}}$ at a temperature where $D = 10^{-12}\ \mathrm{cm^2/s}$, for $t = 1\ \mathrm{h} = 3600\ \mathrm{s}$. Then $\sqrt{Dt} = \sqrt{3.6\times10^{-9}} = 6.0\times10^{-5}\ \mathrm{cm} = 0.60\ \mu\mathrm{m}$.
>
> **(a) Predeposition** with $C_s = 10^{20}\ \mathrm{cm^{-3}}$. Solve $\operatorname{erfc}\big(x_j/(2\sqrt{Dt})\big) = N_B/C_s = 10^{-4}$: $\operatorname{erfc}^{-1}(10^{-4}) = 2.751$, so
> $$
> x_j = 2\times0.60\ \mu\mathrm{m}\times2.751 = 3.30\ \mu\mathrm{m}.
> $$
> The dose introduced is $Q = 2C_s\sqrt{Dt/\pi} = 2\times10^{20}\times6.0\times10^{-5}/1.7725 = 6.8\times10^{15}\ \mathrm{cm^{-2}}$.
>
> **(b) Drive-in** of an implanted dose $Q = 10^{14}\ \mathrm{cm^{-2}}$ under the same $D$ and $t$. Surface concentration $N(0) = Q/\sqrt{\pi Dt} = 10^{14}/(1.7725\times6.0\times10^{-5}) = 9.40\times10^{17}\ \mathrm{cm^{-3}}$. Junction: $N(0)e^{-x_j^2/(4Dt)} = N_B$ gives
> $$
> x_j = \sqrt{4Dt\,\ln\frac{N(0)}{N_B}} = \sqrt{4\times3.6\times10^{-9}\times\ln 94.0}\ \mathrm{cm} = 2.56\ \mu\mathrm{m}.
> $$
> **Scaling check.** Both depths scale as $\sqrt{Dt}$: doubling the diffusion time deepens the junction by $\sqrt2$ (only approximately for (b), since the log factor also changes). This is the thermal-budget constraint of IC processing: every later high-temperature step adds to $Dt$ and moves all earlier junctions.

> [!tip] Engineering note — the same erfc in carrier transport
> Minority carriers injected into a long neutral region obey $\partial\Delta n/\partial t = D_n\Delta n_{xx} - \Delta n/\tau_n$, a diffusion equation with a decay term. In steady state it becomes the ODE $\Delta n'' = \Delta n/L_n^2$, $L_n = \sqrt{D_n\tau_n}$ (the diffusion length), whose decaying exponential solution underlies the ideal-diode equation. The transient (switching) behaviour of a diode is the full PDE, solved with the methods of this chapter after the substitution $\Delta n = e^{-t/\tau_n}w$, which removes the decay term exactly as in the distortionless-line trick of Section 13.5.

### Heat transfer: transient cooling with convection

A plate of thickness $2L_h$ (a ceramic substrate, a heat-spreader, a slab of potting compound) is initially at uniform excess temperature $\theta_i$ above the coolant, and at $t = 0$ both faces are exposed to a fluid with heat-transfer coefficient $h_c$. Measure $x$ from the mid-plane and let $\theta = u - u_\infty$. By symmetry we solve on $0<x<L_h$:

$$
\theta_t = \alpha\,\theta_{xx},\qquad \theta_x(0,t) = 0\ \text{(symmetry)},\qquad -k\,\theta_x(L_h,t) = h_c\,\theta(L_h,t)\ \text{(Robin)},\qquad \theta(x,0) = \theta_i .
$$

Nondimensionalise: $\xi = x/L_h$, **Fourier number** $\mathrm{Fo} = \alpha t/L_h^2$, **Biot number** $\mathrm{Bi} = h_cL_h/k$. The Biot number compares internal conduction resistance $L_h/k$ with surface convection resistance $1/h_c$.

Separation (Section 13.3) with $X'(0) = 0$ gives $X = \cos\mu\xi$, and the Robin condition $X'(1) + \mathrm{Bi}\,X(1) = 0$ gives the eigenvalue equation

$$
\mu\tan\mu = \mathrm{Bi},
$$

with one root $\mu_n$ in each interval $\big((n-1)\pi,\ (n-\tfrac12)\pi\big)$. The generalised Fourier coefficients of the constant initial value (Sturm–Liouville theorem, with $\int_0^1\cos^2\mu_n\xi\,d\xi = \tfrac12 + \tfrac{\sin2\mu_n}{4\mu_n}$) are

$$
C_n = \frac{\int_0^1\cos\mu_n\xi\,d\xi}{\int_0^1\cos^2\mu_n\xi\,d\xi} = \frac{4\sin\mu_n}{2\mu_n + \sin2\mu_n},
$$

so

$$
\frac{\theta(\xi,\mathrm{Fo})}{\theta_i} = \sum_{n=1}^\infty C_n\,e^{-\mu_n^2\mathrm{Fo}}\cos\mu_n\xi .
$$

> [!example] Example 13.35 — Biot number and the lumped-capacitance limit
> Compute $\mu_1$, $C_1$ and the centre temperature at $\mathrm{Fo} = 1$ for $\mathrm{Bi} = 0.1,\ 1,\ 10$.
>
> Solving $\mu\tan\mu = \mathrm{Bi}$ by bisection on $(0,\pi/2)$:
>
> | $\mathrm{Bi}$ | $\mu_1$ | $\mu_1^2$ | $C_1$ | $\mu_2$ | $\theta(0,1)/\theta_i$ (one term) |
> |---|---|---|---|---|---|
> | $0.1$ | $0.3111$ | $0.0968$ | $1.016$ | $3.173$ | $0.922$ |
> | $1$ | $0.8603$ | $0.7402$ | $1.119$ | $3.426$ | $0.534$ |
> | $10$ | $1.4289$ | $2.042$ | $1.262$ | $4.306$ | $0.164$ |
>
> For $\mathrm{Bi} = 1$ the full series (50 terms) gives $\theta(0,1)/\theta_i = 0.53386$ identical to five figures with the one-term value. At $\mathrm{Fo} = 0.2$ the series gives $0.9506$ and one term gives $0.9651$, still within $1.5\%$; the second mode has $\mu_2^2 = 11.7$, so it is negligible once $\mathrm{Fo}\gtrsim0.2$.
>
> **Interpretation.**
> - $\mathrm{Bi}\ll1$: $\mu_1^2\approx\mathrm{Bi}$ (here $0.0968\approx0.1$) and $C_1\approx1$. The solution reduces to $\theta\approx\theta_ie^{-\mathrm{Bi}\,\mathrm{Fo}} = \theta_ie^{-h_ct/(\rho c_pL_h)}$: the **lumped-capacitance** (single RC) model, with thermal time constant $\rho c_pL_h/h_c$ = (heat capacity per unit area)$\times$(surface resistance). The temperature is nearly uniform across the plate. The usual criterion for using a lumped model is $\mathrm{Bi}<0.1$.
> - $\mathrm{Bi}\gg1$: $\mu_1\to\pi/2$, the surface is effectively clamped at the coolant temperature (Dirichlet), and internal conduction limits the cooling.
>
> **Dimensional numbers.** An alumina substrate ($k = 30\ \mathrm{W\,m^{-1}K^{-1}}$, $\rho c_p = 3.43\times10^6\ \mathrm{J\,m^{-3}K^{-1}}$, so $\alpha = 8.74\times10^{-6}\ \mathrm{m^2/s}$) of half-thickness $L_h = 5\ \mathrm{mm}$ in forced water cooling with $h_c = 6000\ \mathrm{W\,m^{-2}K^{-1}}$ has $\mathrm{Bi} = 1$. Its dominant time constant is $L_h^2/(\mu_1^2\alpha) = 2.5\times10^{-5}/(0.7402\times8.74\times10^{-6}) = 3.86\ \mathrm{s}$; a lumped model would wrongly predict $\rho c_pL_h/h_c = 2.86\ \mathrm{s}$ and a uniform temperature.

> [!warning] Common pitfall — lumping when Bi is not small
> The single-time-constant "thermal RC" model is a one-mode truncation that is valid only when $\mathrm{Bi}\lesssim0.1$ *and* $\mathrm{Fo}$ is not small. For a power device on a thick, poorly conducting substrate, $\mathrm{Bi}$ can exceed 1, and the lumped model underestimates both the time constant and the peak junction temperature at short times. Use the distributed solution (or a multi-mode Foster/Cauer network).

### Choosing the method: a summary of the chapter's toolkit

| Situation | Method | Section |
|---|---|---|
| Finite interval/rectangle/disk, homogeneous BCs on all but one side | separation of variables + Fourier series | 13.3, 13.7 |
| Constant non-homogeneous BCs | subtract steady state, then separate | 13.2, 13.4 |
| Sources or time-dependent BCs | eigenfunction expansion + Duhamel (modal ODEs) | 13.7 |
| Infinite line, wave equation | d'Alembert, characteristics | 13.5 |
| Semi-infinite line, wave equation with fixed/free end | method of images (odd/even extension) | 13.5 |
| Infinite/semi-infinite domain, diffusion | heat kernel, erfc similarity solution | 13.4, 13.8 |
| Sinusoidal steady state | phasors turn the PDE into an ODE in space | 13.8 |
| Short-time heat transfer | erfc near each boundary; long time: one-mode series | 13.4 |
| Lossy line with $R/L = G/C$ | exponential substitution → lossless wave | 13.5 |
| Irregular geometry | numerical methods ([Ch. 22](../part-7-numerical-computational/22-numerical-methods.md)), conformal maps ([Ch. 18](../part-5-transforms/18-complex-analysis.md)) | — |

---

## Summary

| Topic | Key formula / fact |
|---|---|
| Telegrapher's equations | $v_x = -Ri - Li_t$, $\ i_x = -Gv - Cv_t$ |
| Transmission-line equation | $v_{xx} = LC\,v_{tt} + (RC+GL)\,v_t + RG\,v$ |
| Limits | lossless: $v_{tt} = c^2v_{xx}$, $c = 1/\sqrt{LC}$; RC: $v_t = v_{xx}/(RC)$ |
| Classification | $a_{11}u_{xx} + 2a_{12}u_{xy} + a_{22}u_{yy}$: $\Delta = a_{12}^2 - a_{11}a_{22}$; $>0$ hyperbolic, $=0$ parabolic, $<0$ elliptic |
| Characteristics | $dy/dx = (a_{12}\pm\sqrt\Delta)/a_{11}$ |
| Boundary conditions | Dirichlet $u = h$; Neumann $\partial u/\partial n = h$; Robin $\beta u + \partial u/\partial n = h$ |
| Steady-state subtraction | $u = u_s(x) + w$, $u_s = T_1 + (T_2-T_1)x/L$ |
| Dirichlet eigenpairs | $\lambda_n = (n\pi/L)^2$, $X_n = \sin(n\pi x/L)$ |
| D–N eigenpairs | $\lambda_n = ((2n-1)\pi/2L)^2$, $X_n = \sin((2n-1)\pi x/2L)$ |
| Robin (D–R) | $\tan\mu = -\mu/h$; slab cooling $\mu\tan\mu = \mathrm{Bi}$ |
| Sturm–Liouville coefficients | $c_n = \int fX_nw\,dx\big/\int X_n^2w\,dx$ |
| Heat equation, D–D | $u = \sum b_ne^{-\alpha n^2\pi^2t/L^2}\sin(n\pi x/L)$, $b_n = \frac2L\int_0^Lf\sin\frac{n\pi x}{L}dx$ |
| Thermal time constant | $\tau_1 = L^2/(\pi^2\alpha)$; diffusion length $\sim\sqrt{\alpha t}$ |
| Heat kernel | $\Phi = (4\pi\alpha t)^{-1/2}e^{-x^2/(4\alpha t)}$, $u = \Phi * f$ |
| Semi-infinite step | $u = T_s\operatorname{erfc}\big(x/(2\sqrt{\alpha t})\big)$ |
| d'Alembert | $u = \tfrac12[f(x-ct) + f(x+ct)] + \frac{1}{2c}\int_{x-ct}^{x+ct}g\,ds$ |
| Line waves | $v = v^+ + v^-$, $i = (v^+ - v^-)/Z_0$, $Z_0 = \sqrt{L/C}$ |
| Reflection coefficient | $\Gamma = (R_L - Z_0)/(R_L + Z_0)$ |
| String / shorted line modes | $f_n = nc/(2L)$; $B_n = \frac{2}{n\pi c}\int_0^Lg\sin\frac{n\pi x}{L}dx$ |
| Wave energy | $E = \tfrac12\int(u_t^2 + c^2u_x^2)\,dx$ conserved |
| Distortionless line | $R/L = G/C$: $v = e^{-(R/L)t}[F(x-ct) + G(x+ct)]$, attenuation $\sqrt{RG}$ Np/m |
| Laplace, rectangle | $u = \sum h_n\frac{\sinh(n\pi y/a)}{\sinh(n\pi b/a)}\sin\frac{n\pi x}{a}$ |
| Laplace, polar products | $1,\ \ln r,\ r^{\pm n}\cos n\phi,\ r^{\pm n}\sin n\phi$ |
| Disk / Poisson integral | $u = \frac{a_0}{2} + \sum(r/a)^n(a_n\cos n\phi + b_n\sin n\phi)$; kernel $\frac{a^2-r^2}{a^2 - 2ar\cos(\phi-\psi) + r^2}$ |
| Coax | $V = V_0\ln(b/r)/\ln(b/a)$, $C' = 2\pi\varepsilon/\ln(b/a)$ |
| Mean value / max principle | $u(\mathbf{x}_0)$ = circle average; extrema on the boundary |
| Duhamel (modal) | $T_n' + \alpha\lambda_nT_n = s_n(t)$, $T_n = b_ne^{-\alpha\lambda_nt} + \int_0^te^{-\alpha\lambda_n(t-\tau)}s_n\,d\tau$ |
| Open RC line delay | $t_{50} = 0.38\,RC\ell^2$ |
| Skin depth | $\delta = \sqrt{2/(\omega\mu\sigma)}$, $H = H_0e^{-x/\delta}\cos(\omega t - x/\delta)$, $R_s = 1/(\sigma\delta)$ |
| Dopant profiles | $C_s\operatorname{erfc}\big(x/2\sqrt{Dt}\big)$; $\frac{Q}{\sqrt{\pi Dt}}e^{-x^2/4Dt}$ |

---

## Further reading

- E. Kreyszig, *Advanced Engineering Mathematics*, 10th ed., Wiley — Chapter 12 (PDEs: separation of variables, heat, wave and Laplace equations, d'Alembert's solution).
- R. Haberman, *Applied Partial Differential Equations with Fourier Series and Boundary Value Problems*, 5th ed., Pearson — the clearest full-length treatment of Sturm–Liouville theory, eigenfunction expansions and Duhamel's principle, including Bessel functions for circular domains.
- W. A. Strauss, *Partial Differential Equations: An Introduction*, 2nd ed., Wiley — classification, characteristics, well-posedness, maximum principles and energy methods at a slightly more mathematical level.
- D. K. Cheng, *Field and Wave Electromagnetics*, 2nd ed., Addison-Wesley — boundary-value problems in electrostatics, transmission-line equations and skin effect, using the methods of this chapter.
