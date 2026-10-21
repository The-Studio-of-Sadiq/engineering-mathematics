# Chapter 23 — Optimization

[← Previous](22-numerical-methods.md) · [Table of Contents](../00-front/table-of-contents.md) · [Next →](24-computational-mathematics.md)

> **Prerequisites:** [Ch. 4](../part-2-calculus/04-differential-calculus.md) (single-variable extrema, §4.10), [Ch. 6](../part-2-calculus/06-multivariable-calculus.md) (gradient, Hessian, constrained extrema and Lagrange multipliers, §6.10–6.11), [Ch. 8](../part-3-linear-algebra/08-matrices-and-linear-systems.md)–[Ch. 10](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md) (linear systems, least squares, positive-definite matrices, SVD), [Ch. 22](22-numerical-methods.md) (Newton's method for roots, convergence order).
> **Learning outcomes:** After this chapter you can
> - turn an engineering design question into a precise optimization problem: decision variables, objective, constraints, and scaling;
> - solve small linear programs graphically and by the simplex method, and read the dual variables as shadow prices;
> - state and apply the first- and second-order optimality conditions, including the Karush–Kuhn–Tucker (KKT) conditions for inequality constraints;
> - run and analyse gradient descent, Newton's method and the BFGS quasi-Newton method, choosing step sizes by line search;
> - recognise convex sets, functions and problems, and explain why convexity turns "a local minimum" into "the global minimum";
> - solve linear and nonlinear least-squares problems, with and without ridge (Tikhonov) regularization, and handle trade-offs between competing objectives with Pareto fronts, weighted sums and the $\varepsilon$-constraint method.

Optimization is the mathematics of *choosing well*. Every design an engineer signs off is the answer to an implicit optimization problem: the resistor values that meet a cutoff frequency at least cost, the generator outputs that meet a load at least fuel, the filter taps that best approximate an ideal response, the parameters of a model that best explain measured data. Chapter 4 handled the one-variable case with "set the derivative to zero", and Chapter 6 extended that to several variables and to equality constraints with Lagrange multipliers. This chapter takes the subject to the point where it becomes a working tool: many variables, inequality constraints, problems with no closed-form solution, and iterative algorithms that a computer runs for you — together with the theory that tells you whether the number the computer returns can be trusted.

Throughout, the standard problem is written

$$
\begin{aligned}
\min_{\mathbf{x}\in\mathbb{R}^n}\quad & f(\mathbf{x})\\
\text{subject to}\quad & g_i(\mathbf{x})\le 0,\quad i=1,\dots,m,\\
& h_k(\mathbf{x})=0,\quad k=1,\dots,p.
\end{aligned}
$$

Maximization of $F$ is minimization of $f=-F$, and a constraint $g(\mathbf{x})\ge0$ is written $-g(\mathbf{x})\le0$, so this single template covers everything.

## 23.1 Objective functions

### Motivation

The hardest step in optimization is usually not the algorithm; it is deciding *what you are optimizing*. "Design a good low-pass filter" is not a mathematical problem. "Choose $R$ and $C$ to minimise board area subject to a cutoff of $1\ \text{kHz}\pm1\,\%$ and an output noise below $0.1\ \mu\text{V}$" is. The objective function encodes what "good" means, and a poorly chosen objective produces a design that is optimal for the wrong question — a perfect answer to a misstatement. This section sets up the vocabulary and shows how engineering goals become functions.

### Definitions

> [!definition] Decision variables, objective, optimum
> The **decision variables** $\mathbf{x}=(x_1,\dots,x_n)$ are the quantities the designer is free to choose. The **objective function** (cost, loss, or merit function) $f:D\subseteq\mathbb{R}^n\to\mathbb{R}$ assigns a single real number to each candidate design. A point $\mathbf{x}^*$ is a
> - **global minimiser** if $f(\mathbf{x}^*)\le f(\mathbf{x})$ for every feasible $\mathbf{x}$;
> - **local minimiser** if there is a $\delta>0$ with $f(\mathbf{x}^*)\le f(\mathbf{x})$ for all feasible $\mathbf{x}$ with $\lVert\mathbf{x}-\mathbf{x}^*\rVert<\delta$;
> - **strict** (local or global) minimiser if the inequality is strict for $\mathbf{x}\neq\mathbf{x}^*$.
>
> The value $f^*=f(\mathbf{x}^*)$ is the **optimal value**. If the infimum of $f$ is not attained (e.g. $f(x)=e^{-x}$ on $\mathbb{R}$) the problem has an optimal value $f^*=0$ but no minimiser.

The distinction between local and global is the central difficulty of nonlinear optimization. Every practical algorithm in this chapter is *local*: it looks at $f$ and its derivatives near the current point and moves downhill. Whether the point where it stops is the global answer depends on the global shape of $f$ — which is why §23.7 on convexity matters so much.

> [!theorem] Weierstrass extreme value theorem
> If $f$ is continuous on a nonempty, closed and bounded (compact) set $S\subset\mathbb{R}^n$, then $f$ attains both a minimum and a maximum on $S$.

A useful variant for unbounded domains: if $f$ is continuous on $\mathbb{R}^n$ and **coercive** — $f(\mathbf{x})\to\infty$ whenever $\lVert\mathbf{x}\rVert\to\infty$ — then $f$ has a global minimiser. (Proof sketch: pick any $\mathbf{x}_0$; coercivity means the sublevel set $\{\mathbf{x}: f(\mathbf{x})\le f(\mathbf{x}_0)\}$ is bounded, it is closed by continuity, so Weierstrass applies to it, and nothing outside it can do better.) Quadratic losses $\lVert\mathbf{A}\mathbf{x}-\mathbf{b}\rVert^2+\lambda\lVert\mathbf{x}\rVert^2$ with $\lambda>0$ are coercive; this is one quiet reason regularization (§23.8) is popular.

### Kinds of objective

| Objective type | Form | Typical EE instance |
|---|---|---|
| Linear | $\mathbf{c}^\mathsf{T}\mathbf{x}$ | fuel cost of generators with linear cost curves; total component cost |
| Quadratic | $\tfrac12\mathbf{x}^\mathsf{T}\mathbf{Q}\mathbf{x}+\mathbf{c}^\mathsf{T}\mathbf{x}$ | ohmic loss $\sum R_k I_k^2$; energy $\tfrac12 CV^2$; squared error |
| Sum of squares | $\sum_i r_i(\mathbf{x})^2$ | curve fitting, filter design in the $L^2$ sense |
| Minimax | $\max_i \lvert r_i(\mathbf{x})\rvert$ | equiripple (Chebyshev) filter design; worst-case tolerance |
| Ratio / log | $\log(\text{SNR})$, $\eta=P_\text{out}/P_\text{in}$ | efficiency, noise figure, capacity $\log_2(1+\text{SNR})$ |
| Expected value | $E[\ell(\mathbf{x},\boldsymbol{\xi})]$ | mean-square error of an estimator with random noise $\boldsymbol\xi$ |

The **norm** you put on an error vector $\mathbf{r}$ is a modelling choice with real consequences. The $L^2$ norm $\sum r_i^2$ is smooth and leads to linear equations (§23.8), but a few large errors dominate it. The $L^1$ norm $\sum\lvert r_i\rvert$ is robust to outliers and can be written as a linear program. The $L^\infty$ norm $\max\lvert r_i\rvert$ controls the worst case, which is what a specification mask demands; it, too, becomes a linear program after an epigraph trick (§23.2).

### Scaling and units

An objective must be a single number, so all its terms must share units. When the natural terms do not — dollars of capital and kilowatt-hours of loss, ohms and farads — you convert them to a common currency (cost per year) or nondimensionalise them by reference values. This is not cosmetic: gradient-based algorithms (§23.5) move in the direction of steepest descent *as measured in the units you chose*, and a variable measured in farads ($10^{-9}$) next to one measured in ohms ($10^{4}$) gives a problem that is numerically ill-conditioned by a factor of $10^{13}$. Rescale every variable so that its plausible range is of order one before handing it to software.

> [!example] Worked Example 23.1 — maximum power transfer as an optimization problem
> A source with open-circuit voltage $V_s$ and internal resistance $R_s$ drives a load $R_L$. Choose $R_L>0$ to maximise the power delivered.
>
> *Step 1 (variables and objective):* the decision variable is $R_L$; the load current is $I=V_s/(R_s+R_L)$, so
> $$
> P(R_L)=I^2R_L=\frac{V_s^2R_L}{(R_s+R_L)^2}.
> $$
> Units: $\mathrm{V^2\,\Omega/\Omega^2}=\mathrm{V^2/\Omega}=\mathrm{W}$. ✓
>
> *Step 2 (stationary point):* by the quotient rule
> $$
> P'(R_L)=V_s^2\frac{(R_s+R_L)^2-2R_L(R_s+R_L)}{(R_s+R_L)^4}=V_s^2\frac{R_s-R_L}{(R_s+R_L)^3},
> $$
> which vanishes only at $R_L=R_s$.
>
> *Step 3 (classify):* $P'>0$ for $R_L<R_s$ and $P'<0$ for $R_L>R_s$, so $R_L=R_s$ is the global maximiser on $(0,\infty)$, with $P_\text{max}=V_s^2/(4R_s)$.
>
> *Step 4 (check):* $V_s=10\ \text{V}$, $R_s=50\ \Omega$: $P(50)=100\cdot50/100^2=0.5\ \text{W}=V_s^2/(4R_s)$ ✓; $P(40)=4000/8100=0.494\ \text{W}$ and $P(60)=6000/12100=0.496\ \text{W}$, both smaller ✓.
>
> *Step 5 (question the objective):* at $R_L=R_s$ the efficiency is only $50\,\%$. A power utility, whose objective is efficiency, chooses $R_L\gg R_s$. Same circuit, different objective, opposite design.

```desmos-graph
left=0; right=200; bottom=0; top=0.6
---
y=\frac{100x}{(50+x)^2}|blue
y=\frac{x}{50+x}\cdot0.6|green|dashed
x=50|red|dotted
(50,0.5)|label:maximum power
```

*Load power $P(R_L)$ (blue, watts) for $V_s=10$ V, $R_s=50\ \Omega$; the dashed green curve is efficiency $R_L/(R_s+R_L)$ scaled by 0.6. The two objectives peak in different places.*

> [!example] Worked Example 23.2 — Kelvin's law for conductor sizing
> A feeder carries a steady current $I$ for $T$ hours per year. A conductor of cross-section $A$ has annualised capital cost $c_1A$ (cost per unit area per year) and resistance $\rho\ell/A$, so its yearly loss cost is $c_2 I^2\rho\ell T/A$ with $c_2$ the price of energy. Find the economic cross-section.
>
> *Step 1 (objective):* collect the constants $k=c_2I^2\rho\ell T$; the yearly cost is
> $$
> f(A)=c_1A+\frac{k}{A},\qquad A>0.
> $$
>
> *Step 2 (stationary point):* $f'(A)=c_1-k/A^2=0\Rightarrow A^*=\sqrt{k/c_1}$.
>
> *Step 3 (classify):* $f''(A)=2k/A^3>0$ for all $A>0$, so $f$ is convex on its domain and $A^*$ is the global minimiser. At the optimum $c_1A^*=k/A^*=\sqrt{c_1k}$: **the annual capital cost equals the annual loss cost** (Kelvin's law, 1881).
>
> *Step 4 (numbers):* $I=200\ \text{A}$, $\ell=1\ \text{km}$, $\rho=1.72\times10^{-8}\ \Omega\,\text{m}$, $T=8760\ \text{h}$, energy $0.10$ per kWh, capital $c_1=0.02$ per mm$^2$ per metre per year, i.e. $20$ per mm$^2$ per year for the whole kilometre. Then the loss cost is
> $$
> \frac{k}{A}=\frac{(200)^2(1.72\times10^{-8})(1000)}{A\,[\text{m}^2]}\ \text{W}\times 8760\ \text{h}\times\frac{0.10}{1000\ \text{Wh}} .
> $$
> With $A$ in mm$^2$ ($1\ \text{mm}^2=10^{-6}\ \text{m}^2$): resistive loss $=688/A\ \text{kW}$, energy $=6.03\times10^{6}/A$ kWh, cost $=6.03\times10^{5}/A$. So $k=6.03\times10^5$ and $A^*=\sqrt{6.03\times10^5/20}=174\ \text{mm}^2$. Check: capital $20\times174=3.47\times10^3$, loss $6.03\times10^5/174=3.47\times10^3$ ✓. One would then pick the next standard size (185 mm$^2$); the cost curve is flat near the optimum ($f$ rises only $0.2\,\%$), so rounding costs almost nothing.

That last remark is general and worth stating: near a smooth minimum, $f(\mathbf{x}^*+\boldsymbol\delta)\approx f^*+\tfrac12\boldsymbol\delta^\mathsf{T}\nabla^2f(\mathbf{x}^*)\boldsymbol\delta$, so a *first-order* error in the design costs only a *second-order* penalty in the objective. Optima are forgiving; this is why engineers can round to preferred values.

> [!example] Worked Example 23.3 — writing a filter-design objective
> We want a linear-phase FIR filter whose amplitude response $A(\omega;\mathbf{h})=h_0+2\sum_{k=1}^{M}h_k\cos k\omega$ approximates an ideal low-pass $D(\omega)$ (1 in the passband $[0,\omega_p]$, 0 in the stopband $[\omega_s,\pi]$). Formulate three candidate objectives.
>
> *Step 1 (variables):* $\mathbf{h}=(h_0,\dots,h_M)\in\mathbb{R}^{M+1}$. $A$ is *linear* in $\mathbf{h}$: $A(\omega;\mathbf{h})=\mathbf{a}(\omega)^\mathsf{T}\mathbf{h}$ with $\mathbf{a}(\omega)=(1,2\cos\omega,\dots,2\cos M\omega)$.
>
> *Step 2 (grid):* sample $\omega_1,\dots,\omega_N$ in the bands (not in the transition band, which is left free) and form the error $e_i(\mathbf{h})=W(\omega_i)\,[\mathbf{a}(\omega_i)^\mathsf{T}\mathbf{h}-D(\omega_i)]$ with a weight $W$ (larger in the stopband if attenuation matters more).
>
> *Step 3 (objectives):*
> - $f_2(\mathbf{h})=\sum_ie_i^2$: least squares — a quadratic, solved by one linear system (§23.8);
> - $f_\infty(\mathbf{h})=\max_i\lvert e_i\rvert$: minimax — the equiripple design; non-smooth, but an LP (§23.2–23.3);
> - $f_1(\mathbf{h})=\sum_i\lvert e_i\rvert$: least absolute error — also an LP.
>
> All three are **convex** in $\mathbf{h}$ (§23.7), so each has a single global optimum that any correct algorithm finds. The choice between them is an engineering judgement about what the specification sheet actually penalises.

> [!warning] Common Pitfall — optimizing a proxy
> Objectives are often proxies: "minimise component count" for cost, "minimise mean-square error" for audible distortion. An optimizer exploits every weakness of the proxy. If the optimal design looks strange — a resistor of $0.3\ \text{m}\Omega$, a filter with huge gain in the transition band — the objective or the constraints are incomplete, not the algorithm.

> [!tip] Engineering Note — sensitivity of the optimal value
> Once you have $\mathbf{x}^*(p)$ for a parameter $p$ (a price, a load, a tolerance), the **envelope theorem** gives the sensitivity cheaply: $\dfrac{d f^*}{dp}=\dfrac{\partial f}{\partial p}\Big|_{\mathbf{x}^*}$ — the indirect effect through $\mathbf{x}^*$ vanishes because $\nabla_\mathbf{x}f=\mathbf{0}$ there. In Example 23.2, $f^*=2\sqrt{c_1k}$ and indeed $\partial f/\partial c_1=A^*=\sqrt{k/c_1}=df^*/dc_1$. With constraints the same role is played by the Lagrange multipliers (§23.4).

## 23.2 Constraints

### Motivation

Unconstrained optimization is the exception in engineering. Resistances are positive, currents are bounded by fuse ratings, generator outputs lie between a technical minimum and a nameplate maximum, and the specification sheet is a list of inequalities: gain at least 40 dB, ripple at most 0.5 dB, power at most 2 W. Constraints change the character of the problem: the optimum is very often *on the boundary* of what is allowed, where the derivative of the objective is not zero at all. Learning to describe constraints precisely — and to see which ones actually matter — is the core of practical optimization.

### Definitions

> [!definition] Feasible set, active constraints
> For the standard problem with inequality constraints $g_i(\mathbf{x})\le0$ and equality constraints $h_k(\mathbf{x})=0$, the **feasible set** is
> $$
> S=\{\mathbf{x}\in\mathbb{R}^n:\ g_i(\mathbf{x})\le0\ \ (i=1,\dots,m),\ \ h_k(\mathbf{x})=0\ \ (k=1,\dots,p)\}.
> $$
> A point in $S$ is **feasible**. At a feasible $\mathbf{x}$ the inequality $g_i$ is **active** (binding, tight) if $g_i(\mathbf{x})=0$ and **inactive** if $g_i(\mathbf{x})<0$. The **active set** is $\mathcal{A}(\mathbf{x})=\{i: g_i(\mathbf{x})=0\}$. Equality constraints are always active. If $S=\varnothing$ the problem is **infeasible**; if $f$ is unbounded below on $S$ it is **unbounded**.

Inactive constraints have no local influence: if $g_i(\mathbf{x}^*)<0$, a small neighbourhood of $\mathbf{x}^*$ satisfies $g_i<0$ too, so $\mathbf{x}^*$ would still be a local minimiser with $g_i$ deleted. The difficulty of inequality-constrained optimization is that you do not know in advance *which* constraints will be active; the KKT conditions of §23.4 are, in effect, a systematic way to guess and verify the active set.

Common constraint types:

| Type | Form | Example |
|---|---|---|
| Bound (box) | $\ell_j\le x_j\le u_j$ | $P_{\min}\le P_g\le P_{\max}$ for a generator |
| Linear inequality | $\mathbf{a}_i^\mathsf{T}\mathbf{x}\le b_i$ | total current from a supply $\le$ rating |
| Linear equality | $\mathbf{A}\mathbf{x}=\mathbf{b}$ | Kirchhoff's current law at every node; power balance |
| Nonlinear | $g(\mathbf{x})\le0$ | $RC=\tau$; $\lvert H(j\omega_s)\rvert\le10^{-2}$; $\sqrt{kT/C}\le v_{\max}$ |
| Integer / discrete | $x_j\in\mathbb{Z}$ or $x_j\in\{\text{E12 values}\}$ | number of parallel cells; preferred resistor values |

Discrete constraints make a problem combinatorial and are outside the scope of this chapter; the standard engineering practice is to solve the continuous relaxation and then round to the nearest available values, checking the specification afterwards (the forgiving-optimum argument of §23.1 says the penalty is second order).

### Handling constraints: four techniques

**1. Elimination.** An equality constraint that can be solved for one variable removes that variable. For $\min f(x,y)$ subject to $h(x,y)=0$ with $y=\varphi(x)$, minimise $F(x)=f(x,\varphi(x))$. This is the cleanest method when it applies, and it is exactly what Chapter 6 did before introducing multipliers ([§6.10](../part-2-calculus/06-multivariable-calculus.md)).

**2. Change of variables.** Simple inequality constraints can often be removed by a substitution: $x>0$ becomes $x=e^u$ with $u\in\mathbb{R}$; $0\le x\le1$ becomes $x=\sin^2u$ or $x=1/(1+e^{-u})$. The resulting problem is unconstrained. For positive physical quantities (resistances, capacitances, frequencies) the logarithmic substitution is particularly natural because tolerances and preferred-value series are logarithmic anyway; in §23.10 it turns a component-sizing problem into a convex one.

**3. Slack variables.** An inequality $g_i(\mathbf{x})\le0$ becomes an equality $g_i(\mathbf{x})+s_i=0$ with a new variable $s_i\ge0$ (the **slack**). Linear programming uses this to put all constraints in equality form (§23.3). The slack is zero exactly when the constraint is active.

**4. Penalty and barrier functions.** Replace the constrained problem by a sequence of unconstrained ones:
- an *exterior penalty* adds $\mu\big(\sum_k h_k(\mathbf{x})^2+\sum_i\max(0,g_i(\mathbf{x}))^2\big)$, with $\mu\to\infty$; iterates approach $S$ from outside;
- an *interior (logarithmic) barrier* adds $-\frac1t\sum_i\ln(-g_i(\mathbf{x}))$, which is finite only strictly inside $S$; as $t\to\infty$ the minimiser approaches the constrained optimum from inside. Barrier methods with Newton steps are the **interior-point methods** that modern LP and convex solvers use.

### The epigraph trick

A non-smooth objective such as $\max_i r_i(\mathbf{x})$ can be made smooth by adding one variable $t$:

$$
\min_{\mathbf{x}}\ \max_i r_i(\mathbf{x})\qquad\Longleftrightarrow\qquad \min_{\mathbf{x},t}\ t\quad\text{subject to}\quad r_i(\mathbf{x})\le t,\ \ i=1,\dots,N.
$$

At the optimum of the right-hand problem $t$ is pushed down until it touches the largest $r_i$, so the two problems have the same solutions. With $r_i$ linear in $\mathbf{x}$, minimax design becomes a linear program. Similarly $\lvert r\rvert\le t$ is the pair of linear inequalities $-t\le r\le t$, which turns $L^1$ and $L^\infty$ fitting into LPs.

> [!example] Worked Example 23.4 — a voltage divider: feasibility first
> A divider $R_1$ (top), $R_2$ (bottom) from $V_\text{in}=5\ \text{V}$ must give $V_\text{out}=2\ \text{V}$ (ratio $\alpha=R_2/(R_1+R_2)=0.4$), draw at most $I_{\max}=100\ \mu\text{A}$, and have Thevenin output resistance $R_1\parallel R_2\le10\ \text{k}\Omega$. Find the feasible designs.
>
> *Step 1 (eliminate the equality):* with total resistance $R_T=R_1+R_2$, the ratio gives $R_2=\alpha R_T$, $R_1=(1-\alpha)R_T$, so there is one free variable $R_T$.
>
> *Step 2 (rewrite the inequalities):*
> - current: $V_\text{in}/R_T\le I_{\max}\Rightarrow R_T\ge 5/10^{-4}=50\ \text{k}\Omega$;
> - output resistance: $R_1R_2/R_T=\alpha(1-\alpha)R_T=0.24R_T\le10\ \text{k}\Omega\Rightarrow R_T\le41.7\ \text{k}\Omega$.
>
> *Step 3 (feasibility):* $50\le R_T\le41.7$ is impossible: **the specification is infeasible.** No algorithm can fix that; the specification must change. (This is why solvers report "infeasible" as a legitimate outcome rather than an error.)
>
> *Step 4 (relax):* allow $I_{\max}=150\ \mu\text{A}$. Then $R_T\ge33.3\ \text{k}\Omega$ and the feasible set is the interval $[33.3,\ 41.7]\ \text{k}\Omega$.
>
> *Step 5 (optimise):* to minimise quiescent power $V_\text{in}^2/R_T$, take the largest feasible $R_T=41.7\ \text{k}\Omega$: $R_1=25.0\ \text{k}\Omega$, $R_2=16.7\ \text{k}\Omega$. The output-resistance constraint is **active**, the current constraint is **inactive**. Check: $5/41.7\text{k}=120\ \mu\text{A}\le150$ ✓; $25\parallel16.7=10.0\ \text{k}\Omega$ ✓; ratio $16.7/41.7=0.40$ ✓.

> [!example] Worked Example 23.5 — quadratic penalty converging to a constrained optimum
> Minimise $f(x,y)=x^2+y^2$ subject to $x+y=1$ using the penalty function $P_\mu(x,y)=x^2+y^2+\mu(x+y-1)^2$.
>
> *Step 1 (exact answer, for reference):* by elimination ($y=1-x$), $F(x)=x^2+(1-x)^2$, $F'(x)=4x-2=0$, so $(x^*,y^*)=(\tfrac12,\tfrac12)$, $f^*=\tfrac12$.
>
> *Step 2 (penalised problem):* $\nabla P_\mu=\mathbf{0}$ gives $2x+2\mu(x+y-1)=0$ and $2y+2\mu(x+y-1)=0$. Subtracting, $x=y$; then $2x+2\mu(2x-1)=0$, so
> $$
> x_\mu=y_\mu=\frac{\mu}{1+2\mu}.
> $$
>
> *Step 3 (behaviour):* $\mu=1$: $x_\mu=0.333$; $\mu=10$: $0.476$; $\mu=100$: $0.4975$; $\mu=10^4$: $0.499975$. The constraint violation is $1-2x_\mu=1/(1+2\mu)$, which goes to zero only as $\mu\to\infty$.
>
> *Step 4 (multiplier for free):* the quantity $2\mu(x_\mu+y_\mu-1)=-2\mu/(1+2\mu)\to-1$. Comparing with the Lagrange condition $\nabla f+\lambda\nabla h=\mathbf{0}$, i.e. $2x+\lambda=0$ at $x=\tfrac12$, gives $\lambda^*=-1$: the penalty term's gradient converges to the multiplier term.
>
> *Step 5 (the catch):* the Hessian of $P_\mu$ is $\begin{pmatrix}2+2\mu&2\mu\\2\mu&2+2\mu\end{pmatrix}$ with eigenvalues $2$ and $2+4\mu$; its condition number $1+2\mu$ grows without bound. Large penalties make the unconstrained subproblems ill-conditioned — the reason practical codes use *augmented Lagrangian* or barrier methods instead of a pure penalty.

```desmos-graph
left=-0.5; right=1.5; bottom=-0.5; top=1.5
---
x+y=1|red
x^2+y^2=0.5|blue
x^2+y^2=0.25|blue|dashed
x^2+y^2=1|blue|dashed
(0.5,0.5)|label:constrained optimum
(0.333,0.333)|label:penalty mu=1
(0.476,0.476)|label:penalty mu=10
```

*Contours of $x^2+y^2$ (blue) and the constraint line (red). The constrained optimum is where a contour is tangent to the line; penalty solutions approach it from the infeasible side as $\mu$ grows.*

> [!example] Worked Example 23.6 — minimax fitting as an LP
> Fit a constant $c$ to three measurements $1.0,\ 1.4,\ 2.0$ so that the worst-case error is as small as possible.
>
> *Step 1 (epigraph form):* minimise $t$ over $(c,t)$ subject to $-t\le c-1.0\le t$, $-t\le c-1.4\le t$, $-t\le c-2.0\le t$ — six linear inequalities.
>
> *Step 2 (which constraints bind):* only the extremes matter: $c-1.0\le t$ and $2.0-c\le t$. Adding, $1.0\le 2t$, so $t\ge0.5$, with equality when $c-1.0=2.0-c$, i.e. $c=1.5$.
>
> *Step 3 (compare):* the minimax constant is the **midrange** $1.5$ (worst error $0.5$), the least-squares constant is the mean $1.467$ (worst error $0.533$), and the $L^1$ constant is the median $1.4$ (worst error $0.6$). Three norms, three answers; each is correct for its objective.

> [!warning] Common Pitfall — strict inequalities
> Constraints in optimization are always non-strict ($\le$, not $<$). A feasible set defined by $x>0$ is open, and $\min x$ over it has no minimiser. If a physical quantity must be strictly positive, impose $x\ge\varepsilon$ with a small physically meaningful $\varepsilon$ (say $1\ \Omega$) or use the substitution $x=e^u$.

> [!tip] Engineering Note — constraint qualification in plain words
> The optimality theory of §23.4 needs the active constraints to be "independent" at the optimum (their gradients linearly independent, or the feasible set having an interior). Redundant constraints — writing the same power-balance equation twice, or a bound that duplicates another — can break this and make solvers report spurious failures. Keep the constraint list minimal and non-redundant.

## 23.3 Linear programming

### Motivation

When the objective and every constraint are linear, the problem is a **linear program** (LP). Linearity sounds restrictive, but LPs model a remarkable range of engineering decisions: dispatching generators with piecewise-linear cost curves, routing power or data flows through a network, allocating bandwidth, scheduling production, and — through the epigraph trick of §23.2 — minimax and $L^1$ fitting. LPs with millions of variables are solved routinely. More importantly for us, LP is the cleanest setting in which to see three ideas that run through all of optimization: the optimum lives on the boundary, an algorithm can walk along that boundary efficiently, and every problem has a **dual** whose variables are prices.

### Standard and canonical forms

> [!definition] Linear program
> A linear program in **inequality (canonical) form** is
> $$
> \max_{\mathbf{x}}\ \mathbf{c}^\mathsf{T}\mathbf{x}\quad\text{subject to}\quad \mathbf{A}\mathbf{x}\le\mathbf{b},\ \ \mathbf{x}\ge\mathbf{0},
> $$
> and in **standard (equality) form**
> $$
> \min_{\mathbf{x}}\ \mathbf{c}^\mathsf{T}\mathbf{x}\quad\text{subject to}\quad \mathbf{A}\mathbf{x}=\mathbf{b},\ \ \mathbf{x}\ge\mathbf{0},
> $$
> with $\mathbf{A}\in\mathbb{R}^{m\times n}$ and, in standard form, $\operatorname{rank}\mathbf{A}=m<n$. Vector inequalities hold componentwise.

Any LP can be brought to either form by these moves: maximise $\mathbf{c}^\mathsf{T}\mathbf{x}$ $\leftrightarrow$ minimise $-\mathbf{c}^\mathsf{T}\mathbf{x}$; turn $\mathbf{a}^\mathsf{T}\mathbf{x}\le b$ into $\mathbf{a}^\mathsf{T}\mathbf{x}+s=b$ with slack $s\ge0$ (or $\mathbf{a}^\mathsf{T}\mathbf{x}\ge b$ into $\mathbf{a}^\mathsf{T}\mathbf{x}-s=b$ with surplus $s\ge0$); replace a free variable $x_j$ by $x_j^+-x_j^-$ with both parts nonnegative.

### Geometry: the optimum is at a vertex

The feasible set of an LP, $S=\{\mathbf{x}:\mathbf{A}\mathbf{x}\le\mathbf{b},\ \mathbf{x}\ge\mathbf{0}\}$, is an intersection of half-spaces — a **convex polyhedron** (a polygon in two dimensions). Its "corners" are **vertices** (extreme points): points of $S$ that are not the midpoint of any two other points of $S$. Algebraically, in standard form a vertex is a **basic feasible solution**: choose $m$ linearly independent columns of $\mathbf{A}$ (the *basis* $\mathbf{B}$), set the other $n-m$ variables to zero, and solve $\mathbf{B}\mathbf{x}_B=\mathbf{b}$; if $\mathbf{x}_B\ge\mathbf{0}$ the point is a vertex.

> [!theorem] Fundamental theorem of linear programming
> For an LP in standard form:
> 1. if the feasible set is nonempty, it has at least one basic feasible solution (vertex);
> 2. if the LP has an optimal solution, it has an optimal solution at a vertex.
>
> Consequently an LP is either infeasible, unbounded, or attains its optimum at one of finitely many vertices (at most $\binom{n}{m}$ of them).

*Proof sketch of part 2.* Let $\mathbf{x}^*$ be optimal but not a vertex. Then $\mathbf{x}^*=\tfrac12(\mathbf{y}+\mathbf{z})$ for distinct feasible $\mathbf{y},\mathbf{z}$, so $\mathbf{d}=\mathbf{y}-\mathbf{x}^*$ is a direction along which one can move both ways and stay feasible for a while. Since $\mathbf{c}^\mathsf{T}(\mathbf{x}^*\pm t\mathbf{d})=\mathbf{c}^\mathsf{T}\mathbf{x}^*\pm t\,\mathbf{c}^\mathsf{T}\mathbf{d}$ cannot drop below the optimum in either direction, $\mathbf{c}^\mathsf{T}\mathbf{d}=0$. Move along $\mathbf{d}$ until one more variable hits zero; the objective is unchanged and the new point has strictly more zero components. Repeating at most $n$ times ends at a vertex with the same optimal value. $\square$

The intuition: the objective $\mathbf{c}^\mathsf{T}\mathbf{x}=z$ is a family of parallel lines (planes). Slide the line in the improving direction $\mathbf{c}$ until it is about to leave $S$; the last contact is a vertex, or a whole edge if the line happens to be parallel to it.

### The graphical method

For two variables:
1. draw each constraint boundary and shade the feasible side;
2. identify the feasible polygon and its vertices (intersections of pairs of active boundaries);
3. draw one objective contour $\mathbf{c}^\mathsf{T}\mathbf{x}=z_0$ and the direction $\mathbf{c}$ of improvement;
4. slide the contour in that direction to the last vertex touched — or, equivalently, evaluate the objective at every vertex and take the best.

> [!example] Worked Example 23.7 — a production LP solved graphically
> A small plant assembles two power-supply modules. Module 1 earns 4 units of profit, module 2 earns 3. Per batch, module 1 uses 1 hour of test time, 2 hours of assembly and 1 kg of copper; module 2 uses 1 hour of test, 1 hour of assembly and 3 kg of copper. Available per week: 8 test hours, 12 assembly hours, 18 kg copper. Maximise profit.
>
> *Step 1 (formulate):* with $x_1,x_2$ batches,
> $$
> \max\ z=4x_1+3x_2\quad\text{s.t.}\quad x_1+x_2\le8,\ \ 2x_1+x_2\le12,\ \ x_1+3x_2\le18,\ \ x_1,x_2\ge0.
> $$
>
> *Step 2 (vertices):* intersect pairs of boundaries and keep the feasible points:
>
> | Vertex | Active constraints | $(x_1,x_2)$ | $z$ |
> |---|---|---|---|
> | A | $x_1=0,\ x_2=0$ | $(0,0)$ | 0 |
> | B | $x_2=0$, assembly | $(6,0)$ | 24 |
> | C | assembly, test | $(4,4)$ | **28** |
> | D | test, copper | $(3,5)$ | 27 |
> | E | copper, $x_1=0$ | $(0,6)$ | 18 |
>
> (For C: subtract $x_1+x_2=8$ from $2x_1+x_2=12$ to get $x_1=4$, then $x_2=4$. For D: subtract $x_1+x_2=8$ from $x_1+3x_2=18$ to get $x_2=5$, $x_1=3$.)
>
> *Step 3 (optimum):* $x^*=(4,4)$, $z^*=28$. At C the test and assembly constraints are active; copper is inactive with slack $18-(4+12)=2$ kg.
>
> *Step 4 (contour check):* the objective gradient $\mathbf{c}=(4,3)$ has slope $-4/3$ for its contours. At C the two active edges have slopes $-1$ (test) and $-2$ (assembly); since $-2<-4/3<-1$, the contour $4x_1+3x_2=28$ touches the polygon only at C. ✓

```desmos-graph
left=-0.5; right=10; bottom=-0.5; top=9
---
y\le8-x\left\{2x+y\le12\right\}\left\{x+3y\le18\right\}\left\{x\ge0\right\}\left\{y\ge0\right\}|blue
x+y=8|green|dashed
2x+y=12|orange|dashed
x+3y=18|purple|dashed
4x+3y=28|red
4x+3y=18|red|dotted
(4,4)|label:C optimum z=28
(3,5)|label:D
(6,0)|label:B
(0,6)|label:E
```

*Feasible region of Example 23.7 (blue) bounded by the test (green), assembly (orange) and copper (purple) constraints. The objective contour $z=18$ (dotted) slides in the direction $(4,3)$ until $z=28$ (solid red) just touches vertex C.*

> [!example] Worked Example 23.8 — a minimisation LP, an unbounded region, and alternative optima
> A utility must buy blocks from two contracts. A block from contract A supplies 1 MW of energy and 1 MW of reserve at cost 2; a block from B supplies 1 MW of energy and 3 MW of reserve at cost 3. It needs at least 4 MW of energy and 6 MW of reserve, and contract A is limited to 5 blocks.
>
> *Step 1 (formulate):* $\min\ 2x+3y$ s.t. $x+y\ge4$, $x+3y\ge6$, $x\le5$, $x,y\ge0$.
>
> *Step 2 (region):* the feasible region is *unbounded* (any large $y$ is allowed), but the objective increases in every feasible direction of escape, so a minimum exists. Its vertices are $(0,4)$, $(3,1)$ [from $x+y=4$, $x+3y=6$], and $(5,\tfrac13)$ [from $x=5$, $x+3y=6$].
>
> *Step 3 (evaluate):* $2x+3y$ gives $12$, $\mathbf{9}$, $11$. Optimum $(3,1)$ with cost 9: both the energy and reserve constraints are active.
>
> *Step 4 (alternative optima):* change B's price so that the cost is $x+3y$. Now the contours are parallel to the reserve edge $x+3y=6$; vertices give $12,\ 6,\ 6$, and **every** point on the segment from $(3,1)$ to $(5,\tfrac13)$ is optimal with cost 6. A solver will report one vertex; the existence of others is signalled by a zero reduced cost for a nonbasic variable in the final simplex tableau.

### The simplex method

The fundamental theorem suggests an algorithm: start at a vertex, move along an edge to an adjacent vertex with a better objective, and stop when no edge improves. This is Dantzig's **simplex method** (1947). Adjacent vertices differ by exchanging one basic variable for one nonbasic variable, which is a single Gauss–Jordan pivot ([Ch. 8](../part-3-linear-algebra/08-matrices-and-linear-systems.md)).

Consider the canonical maximisation problem with $\mathbf{b}\ge\mathbf{0}$. Adding slacks $\mathbf{s}$ gives $\mathbf{A}\mathbf{x}+\mathbf{s}=\mathbf{b}$, and the point $\mathbf{x}=\mathbf{0},\ \mathbf{s}=\mathbf{b}$ is an immediate basic feasible solution (the slacks form the basis). Write the objective as a row $z-\mathbf{c}^\mathsf{T}\mathbf{x}=0$. The **tableau** stores the constraint rows and this objective row.

> [!definition] Simplex iteration (maximisation, tableau form)
> 1. **Optimality test.** If every entry of the objective row (the **reduced costs**, here $-\bar c_j$) is $\ge0$, stop: the current vertex is optimal.
> 2. **Entering variable.** Otherwise choose a column $q$ with a negative entry (Dantzig's rule: the most negative). Increasing $x_q$ from zero increases $z$.
> 3. **Ratio test (leaving variable).** For each row $i$ with positive entry $a_{iq}>0$, compute $b_i/a_{iq}$. The smallest ratio identifies the row $p$ whose basic variable first reaches zero as $x_q$ grows. If no $a_{iq}>0$, the LP is **unbounded**.
> 4. **Pivot.** Divide row $p$ by $a_{pq}$, then subtract multiples of it from every other row (including the objective row) to zero column $q$. Return to step 1.

Why it works: in the current basis every basic variable is expressed in terms of the nonbasic ones; the objective row says how $z$ changes per unit increase of each nonbasic variable. The ratio test increases the entering variable as far as possible *without making any basic variable negative*, so feasibility is preserved and $z$ never decreases. With finitely many bases and strict improvement at each step (absent degeneracy), the method terminates.

> [!example] Worked Example 23.9 — simplex on Example 23.7
> *Initial tableau* (basis $s_1,s_2,s_3$; vertex A):
>
> | Basis | $x_1$ | $x_2$ | $s_1$ | $s_2$ | $s_3$ | RHS |
> |---|---|---|---|---|---|---|
> | $s_1$ | 1 | 1 | 1 | 0 | 0 | 8 |
> | $s_2$ | 2 | 1 | 0 | 1 | 0 | 12 |
> | $s_3$ | 1 | 3 | 0 | 0 | 1 | 18 |
> | $z$ | $-4$ | $-3$ | 0 | 0 | 0 | 0 |
>
> *Iteration 1.* Most negative reduced cost: $x_1$ ($-4$). Ratios: $8/1=8$, $12/2=\mathbf{6}$, $18/1=18$; $s_2$ leaves. Pivot on the entry 2:
>
> | Basis | $x_1$ | $x_2$ | $s_1$ | $s_2$ | $s_3$ | RHS |
> |---|---|---|---|---|---|---|
> | $s_1$ | 0 | 0.5 | 1 | $-0.5$ | 0 | 2 |
> | $x_1$ | 1 | 0.5 | 0 | 0.5 | 0 | 6 |
> | $s_3$ | 0 | 2.5 | 0 | $-0.5$ | 1 | 12 |
> | $z$ | 0 | $-1$ | 0 | 2 | 0 | 24 |
>
> This is vertex B $=(6,0)$, $z=24$.
>
> *Iteration 2.* $x_2$ has reduced cost $-1$. Ratios: $2/0.5=\mathbf{4}$, $6/0.5=12$, $12/2.5=4.8$; $s_1$ leaves. Pivot on 0.5:
>
> | Basis | $x_1$ | $x_2$ | $s_1$ | $s_2$ | $s_3$ | RHS |
> |---|---|---|---|---|---|---|
> | $x_2$ | 0 | 1 | 2 | $-1$ | 0 | 4 |
> | $x_1$ | 1 | 0 | $-1$ | 1 | 0 | 4 |
> | $s_3$ | 0 | 0 | $-5$ | 2 | 1 | 2 |
> | $z$ | 0 | 0 | 2 | 1 | 0 | 28 |
>
> *Optimality.* All reduced costs are $\ge0$: optimum $x_1=4$, $x_2=4$, $s_3=2$ (spare copper), $z^*=28$ — vertex C, as found graphically. The path was A → B → C: two pivots instead of examining all five vertices.
>
> *Bonus reading.* The objective-row entries under the slack columns, $(2,1,0)$, are the **shadow prices** of test time, assembly time and copper — the dual solution, explained next.

**Starting when the origin is infeasible.** If some constraints are $\ge$ or equalities, $\mathbf{x}=\mathbf{0}$ is not feasible and there is no obvious initial basis. The **two-phase method** first minimises the sum of artificial variables added to those rows (phase I); a zero optimum gives a feasible vertex for phase II, a positive optimum proves the LP infeasible. The **big-M method** does both at once by charging each artificial variable a large cost $M$.

**Degeneracy and cycling.** If a ratio test ties, some basic variable becomes zero and the next pivot may not improve $z$. In rare constructed cases the method can then cycle forever; Bland's rule (always choose the lowest-indexed eligible variable) provably prevents this.

**Complexity.** The simplex method can take exponentially many pivots on contrived examples (Klee–Minty), but on practical problems it typically needs a small multiple of $m$ pivots. **Interior-point methods** (Karmarkar 1984 and successors), which follow a smooth path through the interior using barrier functions (§23.2) and Newton steps (§23.6), have polynomial worst-case complexity and are preferred for very large sparse problems.

### Duality

Every LP has a companion LP, its **dual**, obtained by asking: *what is the best upper bound on the primal objective that can be proved by combining the constraints?*

Take Example 23.7. Multiply the three constraints by nonnegative numbers $y_1,y_2,y_3$ and add:

$$
(y_1+2y_2+y_3)\,x_1+(y_1+y_2+3y_3)\,x_2\ \le\ 8y_1+12y_2+18y_3.
$$

If the coefficients on the left are at least the profit coefficients, $y_1+2y_2+y_3\ge4$ and $y_1+y_2+3y_3\ge3$, then (because $x_1,x_2\ge0$) the left side is at least $4x_1+3x_2=z$. So every such $\mathbf{y}$ gives the bound $z\le 8y_1+12y_2+18y_3$. The tightest bound is the dual LP:

$$
\min\ 8y_1+12y_2+18y_3\quad\text{s.t.}\quad y_1+2y_2+y_3\ge4,\ \ y_1+y_2+3y_3\ge3,\ \ \mathbf{y}\ge\mathbf{0}.
$$

> [!definition] Primal–dual pair
> $$
> \text{(P)}\ \ \max\ \mathbf{c}^\mathsf{T}\mathbf{x}\ \ \text{s.t.}\ \ \mathbf{A}\mathbf{x}\le\mathbf{b},\ \mathbf{x}\ge\mathbf{0}
> \qquad\qquad
> \text{(D)}\ \ \min\ \mathbf{b}^\mathsf{T}\mathbf{y}\ \ \text{s.t.}\ \ \mathbf{A}^\mathsf{T}\mathbf{y}\ge\mathbf{c},\ \mathbf{y}\ge\mathbf{0}.
> $$
> Each primal constraint corresponds to a dual variable and each primal variable to a dual constraint. The dual of the dual is the primal.

> [!theorem] LP duality
> 1. **Weak duality.** For any primal-feasible $\mathbf{x}$ and dual-feasible $\mathbf{y}$: $\mathbf{c}^\mathsf{T}\mathbf{x}\le\mathbf{b}^\mathsf{T}\mathbf{y}$.
> 2. **Strong duality.** If either problem has a finite optimum, so does the other, and $\mathbf{c}^\mathsf{T}\mathbf{x}^*=\mathbf{b}^\mathsf{T}\mathbf{y}^*$.
> 3. **Complementary slackness.** Feasible $\mathbf{x}$ and $\mathbf{y}$ are both optimal if and only if
> $$
> y_i\,(b_i-\mathbf{a}_i^\mathsf{T}\mathbf{x})=0\ \ \text{for every } i,\qquad x_j\,\big((\mathbf{A}^\mathsf{T}\mathbf{y})_j-c_j\big)=0\ \ \text{for every } j.
> $$
> 4. If the primal is unbounded the dual is infeasible, and vice versa (both can be infeasible).

*Proof of weak duality and complementary slackness.* For feasible $\mathbf{x},\mathbf{y}$,

$$
\mathbf{c}^\mathsf{T}\mathbf{x}\ \le\ (\mathbf{A}^\mathsf{T}\mathbf{y})^\mathsf{T}\mathbf{x}\ =\ \mathbf{y}^\mathsf{T}\mathbf{A}\mathbf{x}\ \le\ \mathbf{y}^\mathsf{T}\mathbf{b},
$$

the first inequality using $\mathbf{x}\ge\mathbf{0}$ and $\mathbf{A}^\mathsf{T}\mathbf{y}\ge\mathbf{c}$, the second using $\mathbf{y}\ge\mathbf{0}$ and $\mathbf{A}\mathbf{x}\le\mathbf{b}$. Equality $\mathbf{c}^\mathsf{T}\mathbf{x}=\mathbf{b}^\mathsf{T}\mathbf{y}$ — which by weak duality certifies that both are optimal — holds exactly when both inequalities are equalities, i.e. when each product of a nonnegative variable and a nonnegative slack vanishes. Strong duality is deeper (it follows from Farkas' lemma, or from the fact that the simplex method terminates with a dual-feasible objective row); we state it without proof. $\square$

Complementary slackness says: **a constraint with spare capacity has zero price, and a resource with positive price is fully used.**

> [!example] Worked Example 23.10 — solving the dual by complementary slackness
> Use the primal optimum $\mathbf{x}^*=(4,4)$ of Example 23.7 to find the dual optimum.
>
> *Step 1 (inactive primal constraint):* copper has slack 2, so $y_3=0$.
>
> *Step 2 (positive primal variables):* $x_1,x_2>0$, so both dual constraints are tight: $y_1+2y_2=4$ and $y_1+y_2=3$. Hence $y_2=1$, $y_1=2$.
>
> *Step 3 (verify):* $\mathbf{y}^*=(2,1,0)\ge\mathbf{0}$ ✓ (dual feasible); dual objective $8(2)+12(1)+18(0)=28=z^*$ ✓. These are exactly the slack-column entries of the final simplex tableau.
>
> *Step 4 (interpret):* one more test hour is worth 2 units of profit, one more assembly hour 1 unit, more copper nothing. Check with a re-solve: with 13 assembly hours the optimum moves to $(5,3)$ and $z=29=28+1$ ✓. With 8.4 test hours it moves to $(3.6,4.8)$ and $z=28.8=28+2(0.4)$ ✓.
>
> *Step 5 (range of validity):* shadow prices are valid only while the optimal basis stays the same. Adding $\Delta$ test hours changes the basic variables by $\Delta$ times the $s_1$ column of the final tableau, $(2,-1,-5)$: copper slack becomes $2-5\Delta$, which hits zero at $\Delta=0.4$. Beyond 8.4 test hours copper becomes binding and extra test time is worth nothing: with 8.5 hours the optimum is still $28.8$.

> [!tip] Engineering Note — locational marginal prices
> In an electricity market the operator solves an LP (or QP) that dispatches generators at minimum cost subject to power balance at each bus and line flow limits. The dual variable of the balance constraint at bus $k$ is the **locational marginal price**: the cost of serving one more MW at that bus. When a transmission line is congested its limit constraint is active, its dual variable is positive, and prices on either side of the line differ — a direct, monetised reading of complementary slackness.

> [!warning] Common Pitfall — forgetting nonnegativity or sign conventions
> Most LP software assumes $\mathbf{x}\ge\mathbf{0}$ by default (SciPy's `linprog` does). A variable that can be negative — a current direction, a voltage relative to ground — must be declared free, or the solver will silently clamp it at zero and return a wrong "optimum". Equally, solvers differ in the sign they report for dual values (marginals); check against a one-unit perturbation as in Example 23.10.

A minimal Python check of Examples 23.7 and 23.10 (`linprog` minimises, so the profit is negated):

```python
from scipy.optimize import linprog
res = linprog(c=[-4, -3],
              A_ub=[[1, 1], [2, 1], [1, 3]], b_ub=[8, 12, 18],
              bounds=[(0, None), (0, None)])
print(res.x, -res.fun)          # [4. 4.] 28.0
print(-res.ineqlin.marginals)   # shadow prices [2. 1. 0.]
```

## 23.4 Nonlinear optimization

### Motivation

Most engineering objectives are nonlinear: losses go as $I^2R$, gain as a ratio of polynomials in frequency, efficiency as a ratio of powers, model error as a sum of squared exponentials. Before choosing an algorithm we need to know what a solution *looks like* — which equations it satisfies and how to tell a minimum from a saddle. Those **optimality conditions** do double duty: for small problems they are solved by hand to give the answer, and for large ones they are the stopping test that every algorithm checks. This section collects them, first without constraints, then with equality constraints (Lagrange, already met in [Ch. 6 §6.11](../part-2-calculus/06-multivariable-calculus.md)), and finally with inequality constraints (Karush–Kuhn–Tucker).

### Unconstrained optimality conditions

Let $f:\mathbb{R}^n\to\mathbb{R}$ be twice continuously differentiable, with gradient $\nabla f(\mathbf{x})$ and Hessian $\nabla^2 f(\mathbf{x})=\big[\partial^2f/\partial x_i\partial x_j\big]$ ([Ch. 6](../part-2-calculus/06-multivariable-calculus.md)). Everything rests on the second-order Taylor expansion

$$
f(\mathbf{x}^*+\mathbf{d})=f(\mathbf{x}^*)+\nabla f(\mathbf{x}^*)^\mathsf{T}\mathbf{d}+\tfrac12\,\mathbf{d}^\mathsf{T}\nabla^2f(\mathbf{x}^*)\,\mathbf{d}+o(\lVert\mathbf{d}\rVert^2).
$$

> [!theorem] Optimality conditions without constraints
> 1. **First-order necessary condition (FONC).** If $\mathbf{x}^*$ is a local minimiser, then $\nabla f(\mathbf{x}^*)=\mathbf{0}$ ($\mathbf{x}^*$ is a **stationary point**).
> 2. **Second-order necessary condition (SONC).** If $\mathbf{x}^*$ is a local minimiser, then $\nabla^2f(\mathbf{x}^*)$ is positive semidefinite.
> 3. **Second-order sufficient condition (SOSC).** If $\nabla f(\mathbf{x}^*)=\mathbf{0}$ and $\nabla^2f(\mathbf{x}^*)$ is positive definite, then $\mathbf{x}^*$ is a strict local minimiser.

*Proof.* (1) If $\mathbf{g}=\nabla f(\mathbf{x}^*)\ne\mathbf{0}$, take $\mathbf{d}=-t\mathbf{g}$: $f(\mathbf{x}^*-t\mathbf{g})=f(\mathbf{x}^*)-t\lVert\mathbf{g}\rVert^2+O(t^2)<f(\mathbf{x}^*)$ for small $t>0$, contradicting minimality. (2) With $\nabla f=\mathbf{0}$, if some $\mathbf{d}$ had $\mathbf{d}^\mathsf{T}\nabla^2f\,\mathbf{d}<0$, then $f(\mathbf{x}^*+t\mathbf{d})-f(\mathbf{x}^*)=\tfrac12t^2\mathbf{d}^\mathsf{T}\nabla^2f\,\mathbf{d}+o(t^2)<0$ for small $t$. (3) If $\nabla^2f(\mathbf{x}^*)\succ0$ with smallest eigenvalue $\lambda_{\min}>0$, then $\tfrac12\mathbf{d}^\mathsf{T}\nabla^2f\,\mathbf{d}\ge\tfrac12\lambda_{\min}\lVert\mathbf{d}\rVert^2$, which beats the $o(\lVert\mathbf{d}\rVert^2)$ remainder for $\mathbf{d}$ small enough, so $f(\mathbf{x}^*+\mathbf{d})>f(\mathbf{x}^*)$ for small $\mathbf{d}\ne\mathbf{0}$. $\square$

Definiteness is checked by eigenvalues or, for small matrices, by leading principal minors ([Ch. 10 §10.6](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md)). At a stationary point:

| Hessian | Conclusion |
|---|---|
| positive definite (all $\lambda_i>0$) | strict local minimum |
| negative definite (all $\lambda_i<0$) | strict local maximum |
| indefinite (eigenvalues of both signs) | saddle point |
| semidefinite and singular | inconclusive — look at higher-order terms |

> [!example] Worked Example 23.11 — stationary points and their type
> Classify the stationary points of $f(x,y)=x^3-3x+y^2$.
>
> *Step 1 (gradient):* $\nabla f=(3x^2-3,\ 2y)=\mathbf{0}$ gives $x=\pm1$, $y=0$.
>
> *Step 2 (Hessian):* $\nabla^2f=\begin{pmatrix}6x&0\\0&2\end{pmatrix}$.
> - At $(1,0)$: $\operatorname{diag}(6,2)\succ0$ → strict local minimum, $f=-2$.
> - At $(-1,0)$: $\operatorname{diag}(-6,2)$ indefinite → saddle, $f=2$.
>
> *Step 3 (global?):* along $y=0$, $f(x,0)=x^3-3x\to-\infty$ as $x\to-\infty$. The local minimum is **not** global; $f$ is unbounded below. A descent algorithm started at $(-3,0)$ never finds $(1,0)$ — it runs off to $-\infty$.

```desmos-graph
left=-3; right=3; bottom=-2.5; top=2.5
---
x^3-3x+y^2=-1.5|blue
x^3-3x+y^2=-1|blue
x^3-3x+y^2=0|green
x^3-3x+y^2=2|red
x^3-3x+y^2=3|orange|dashed
(1,0)|label:local min
(-1,0)|label:saddle
```

*Level curves of $f=x^3-3x+y^2$. Closed curves surround the local minimum at $(1,0)$; the red level $f=2$ crosses itself at the saddle $(-1,0)$, and to the left the function falls without bound.*

> [!warning] Common Pitfall — a semidefinite Hessian proves nothing
> $f_1=x^4+y^2$, $f_2=x^3+y^2$ and $f_3=-x^4+y^2$ all have $\nabla f=\mathbf{0}$ and Hessian $\operatorname{diag}(0,2)$ at the origin. The origin is a strict minimum of $f_1$, a degenerate saddle of $f_2$ and a saddle of $f_3$. When the Hessian is singular, the second-order test is silent; examine $f$ along the null direction (here the $x$-axis).

### Equality constraints: the Lagrange conditions

For $\min f(\mathbf{x})$ subject to $h_k(\mathbf{x})=0$ ($k=1,\dots,p$), Chapter 6 showed that at a constrained extremum the gradient of $f$ must be a linear combination of the constraint gradients, provided those gradients are linearly independent ([Ch. 6 §6.11](../part-2-calculus/06-multivariable-calculus.md)). With the **Lagrangian**

$$
\mathcal{L}(\mathbf{x},\boldsymbol\lambda)=f(\mathbf{x})+\sum_{k=1}^p\lambda_k h_k(\mathbf{x}),
$$

the conditions are $\nabla_\mathbf{x}\mathcal{L}=\mathbf{0}$ and $\nabla_{\boldsymbol\lambda}\mathcal{L}=\mathbf{h}(\mathbf{x})=\mathbf{0}$: $n+p$ equations in $n+p$ unknowns. The geometric reason: any feasible direction $\mathbf{d}$ is tangent to every constraint surface, $\nabla h_k^\mathsf{T}\mathbf{d}=0$, and at a minimum $f$ cannot decrease along such directions, so $\nabla f\perp\mathbf{d}$ for all tangent $\mathbf{d}$ — that is, $\nabla f$ lies in the span of the $\nabla h_k$. (Chapter 6 writes $\nabla f=\lambda\nabla h$; with our sign convention $\nabla f+\lambda\nabla h=\mathbf{0}$ the multiplier simply changes sign.)

The multiplier is a **sensitivity**: if the constraint is $h(\mathbf{x})=b$ rather than $0$, then $\dfrac{df^*}{db}=-\lambda^*$ in our convention. In Example 23.5, relaxing $x+y=1$ to $x+y=b$ gives $f^*=b^2/2$ and $df^*/db=1=-\lambda^*$ at $b=1$ ✓.

### Inequality constraints: the KKT conditions

Now add inequalities $g_i(\mathbf{x})\le0$. At a candidate $\mathbf{x}^*$ only the active ones matter. For an active $g_i$, a direction $\mathbf{d}$ keeps us feasible to first order if $\nabla g_i^\mathsf{T}\mathbf{d}\le0$ (we may move *into* the feasible side, not just along the boundary). At a minimum, no feasible direction may decrease $f$:

$$
\nabla g_i(\mathbf{x}^*)^\mathsf{T}\mathbf{d}\le0\ \ \text{for all active } i\quad\Longrightarrow\quad\nabla f(\mathbf{x}^*)^\mathsf{T}\mathbf{d}\ge0.
$$

By Farkas' lemma (the same theorem of alternatives that underlies LP duality), this implication holds exactly when $-\nabla f(\mathbf{x}^*)$ is a **nonnegative** combination of the active constraint gradients: $-\nabla f=\sum_{i\in\mathcal{A}}\mu_i\nabla g_i$ with $\mu_i\ge0$. Geometrically: the downhill direction $-\nabla f$ points straight out of the feasible set, into the cone spanned by the outward normals of the active walls. If it pointed anywhere else, there would be a feasible downhill direction.

Setting $\mu_i=0$ for inactive constraints and including equality constraints gives the central theorem of constrained optimization.

> [!theorem] Karush–Kuhn–Tucker (KKT) conditions
> Let $\mathbf{x}^*$ be a local minimiser of $f$ subject to $g_i(\mathbf{x})\le0$ ($i=1,\dots,m$) and $h_k(\mathbf{x})=0$ ($k=1,\dots,p$), all continuously differentiable, and suppose a **constraint qualification** holds at $\mathbf{x}^*$ — for instance LICQ: the gradients of the active inequalities and of all equalities are linearly independent. Then there exist multipliers $\mu_i$ and $\lambda_k$ such that
> $$
> \begin{aligned}
> &\text{stationarity:} && \nabla f(\mathbf{x}^*)+\sum_{i=1}^m\mu_i\nabla g_i(\mathbf{x}^*)+\sum_{k=1}^p\lambda_k\nabla h_k(\mathbf{x}^*)=\mathbf{0},\\
> &\text{primal feasibility:} && g_i(\mathbf{x}^*)\le0,\quad h_k(\mathbf{x}^*)=0,\\
> &\text{dual feasibility:} && \mu_i\ge0,\\
> &\text{complementary slackness:} && \mu_i\,g_i(\mathbf{x}^*)=0\quad\text{for every } i.
> \end{aligned}
> $$
> If the problem is **convex** ($f$ and all $g_i$ convex, all $h_k$ affine), the KKT conditions are also **sufficient**: any KKT point is a global minimiser (§23.7).

Remarks:
- **Complementary slackness** encodes the active-set logic: either a constraint is active ($g_i=0$) or its multiplier is zero. It is the nonlinear version of the LP statement "spare capacity has zero price".
- **The sign of $\mu_i$ matters**, unlike $\lambda_k$. A negative multiplier on an active inequality says the objective would improve by moving *into* the interior, so the constraint should not be active — the guess of active set was wrong.
- **Sensitivity.** If $g_i(\mathbf{x})\le b_i$, then $\partial f^*/\partial b_i=-\mu_i^*$: relaxing a constraint by one unit lowers the optimal cost by about $\mu_i^*$. The multipliers are shadow prices, exactly as in LP; for an LP, the KKT conditions *are* primal feasibility, dual feasibility and complementary slackness.
- **Second-order conditions** for constrained problems test the Hessian of the Lagrangian, $\nabla^2_{\mathbf{x}\mathbf{x}}\mathcal{L}$, only on directions tangent to the active constraints (the bordered-Hessian test of Ch. 6 is the two-variable special case).

### A recipe for small problems

1. Write all constraints as $g_i\le0$, $h_k=0$ and form $\mathcal{L}=f+\sum\mu_ig_i+\sum\lambda_kh_k$.
2. Solve the unconstrained problem; if its minimiser is feasible and the problem is convex, you are done (all $\mu_i=0$).
3. Otherwise, guess an active set $\mathcal{A}$. Solve stationarity together with $g_i=0$ ($i\in\mathcal{A}$), $h_k=0$, and $\mu_i=0$ ($i\notin\mathcal{A}$).
4. Accept the guess if all inactive constraints are satisfied and all $\mu_i\ge0$; otherwise try another active set. With $m$ inequalities there are at most $2^m$ guesses — fine by hand for $m\le3$; algorithms (active-set and interior-point methods) organise this search for large $m$.

> [!example] Worked Example 23.12 — projection onto a half-plane
> Minimise $f=(x-2)^2+(y-1)^2$ subject to $g_1=x+y-2\le0$ and $g_2=-x\le0$.
>
> *Step 1 (unconstrained):* the minimiser $(2,1)$ has $g_1=1>0$: infeasible, so $g_1$ must be active.
>
> *Step 2 (guess $\mathcal{A}=\{1\}$):* stationarity $\big(2(x-2)+\mu_1,\ 2(y-1)+\mu_1\big)=\mathbf{0}$ gives $x-2=y-1$, i.e. $x=y+1$. With $x+y=2$: $y=\tfrac12$, $x=\tfrac32$, and $\mu_1=-2(x-2)=1$.
>
> *Step 3 (verify):* $\mu_1=1\ge0$ ✓; inactive $g_2=-\tfrac32<0$ ✓ with $\mu_2=0$. KKT point $(1.5,0.5)$, $f^*=0.5$. The problem is convex, so this is the global minimum: the closest point to $(2,1)$ in the region.
>
> *Step 4 (sensitivity):* relax to $x+y\le2+\varepsilon$; the projection becomes $(1.5+\varepsilon/2,\ 0.5+\varepsilon/2)$ and $f^*=2\big(\tfrac{1-\varepsilon}{2}\big)^2=\tfrac12(1-\varepsilon)^2$, so $df^*/d\varepsilon=-1=-\mu_1$ at $\varepsilon=0$ ✓.

> [!example] Worked Example 23.13 — the minimum-norm point of a half-plane
> Minimise $f=x^2+y^2$ subject to $x+2y\ge2$.
>
> *Step 1 (standard form):* $g=2-x-2y\le0$. The unconstrained minimiser $(0,0)$ has $g=2>0$, infeasible.
>
> *Step 2 (active):* stationarity $(2x-\mu,\ 2y-2\mu)=\mathbf{0}$ gives $x=\mu/2$, $y=\mu$. Then $x+2y=\tfrac52\mu=2\Rightarrow\mu=\tfrac45$, $(x,y)=(\tfrac25,\tfrac45)$.
>
> *Step 3 (check):* $\mu=0.8\ge0$ ✓, $f^*=\tfrac4{25}+\tfrac{16}{25}=\tfrac45$. Geometric check: the distance from the origin to the line $x+2y=2$ is $2/\sqrt5$, whose square is $\tfrac45$ ✓. And $-\nabla f=(-0.8,-1.6)=\mu\nabla g=0.8(-1,-2)$: the downhill direction points out of the feasible set along the outward normal, as the theory says.

> [!example] Worked Example 23.14 — KKT with case analysis
> Minimise $f=(x-3)^2+(y-2)^2$ subject to $g_1=x^2+y^2-5\le0$ (a disc) and $g_2=x+2y-4\le0$ (a half-plane).
>
> *Step 1 (unconstrained):* $(3,2)$ violates both ($13>5$, $7>4$).
>
> *Step 2 (only $g_1$ active):* the closest point of the circle to $(3,2)$ is $(3,2)\sqrt5/\sqrt{13}=(1.861,1.240)$. But $g_2=1.861+2.481-4=0.342>0$: infeasible, reject.
>
> *Step 3 (only $g_2$ active):* the projection onto the line $x+2y=4$ is $(3,2)-\tfrac{3+4-4}{5}(1,2)=(2.4,0.8)$. But $g_1=5.76+0.64-5=1.4>0$: infeasible, reject.
>
> *Step 4 (both active):* $x=4-2y$ in $x^2+y^2=5$ gives $5y^2-16y+11=0$, $y=1$ or $y=\tfrac{11}5$. Take $(2,1)$ (the other root $(-0.4,2.2)$ is farther from $(3,2)$; one can check its multipliers are not both nonnegative). Stationarity at $(2,1)$:
> $$
> \begin{pmatrix}-2\\-2\end{pmatrix}+\mu_1\begin{pmatrix}4\\2\end{pmatrix}+\mu_2\begin{pmatrix}1\\2\end{pmatrix}=\mathbf{0}
> \ \Rightarrow\ 4\mu_1+\mu_2=2,\ \ 2\mu_1+2\mu_2=2\ \Rightarrow\ \mu_1=\tfrac13,\ \mu_2=\tfrac23.
> $$
>
> *Step 5 (accept):* both multipliers are positive ✓, so $(2,1)$ is a KKT point with $f^*=2$. The problem is convex (disc and half-plane are convex sets, $f$ is convex), so it is the global minimum. A numerical solver (`scipy.optimize.minimize` with SLSQP) returns $(2.000,1.000)$, $f=2.000$ ✓.
>
> *Step 6 (sensitivity):* enlarging the disc to radius$^2$ $5.1$ lowers $f^*$ by about $\mu_1\times0.1=0.033$; the re-solved value is $1.968$, a drop of $0.032$ ✓.

```desmos-graph
left=-1; right=4; bottom=-1; top=3.5
---
x^2+y^2\le5\left\{x+2y\le4\right\}|blue
x^2+y^2=5|blue|dashed
x+2y=4|green|dashed
(x-3)^2+(y-2)^2=2|red
(x-3)^2+(y-2)^2=0.5|red|dotted
(3,2)|label:unconstrained min
(2,1)|label:KKT point
```

*Example 23.14: the feasible set (blue) is the part of the disc below the green line. The smallest contour of $f$ (red) that touches it does so at the corner $(2,1)$, where both constraints are active.*

> [!tip] Engineering Note — water-filling
> Allocating total power $P$ across parallel channels with noise levels $N_i$ to maximise capacity, $\max\sum_i\log(1+p_i/N_i)$ s.t. $\sum p_i=P$, $p_i\ge0$, is a convex problem whose KKT conditions give $p_i=\max(0,\ \nu-N_i)$ for a common "water level" $\nu$. Channels whose noise exceeds the level get no power — their nonnegativity constraints are active. This one-line KKT result underlies power allocation in OFDM and MIMO systems.

## 23.5 Gradient-based methods

### Motivation

The optimality conditions are equations, but for realistic problems — a 64-tap filter, a 20-parameter device model, a neural-network equaliser with $10^5$ weights — they cannot be solved in closed form. Instead we **iterate**: start from a guess $\mathbf{x}_0$ and generate $\mathbf{x}_1,\mathbf{x}_2,\dots$ with $f$ decreasing at each step. The simplest and most widely used idea is to walk downhill along the negative gradient. Gradient descent is cheap (one gradient per step, no matrices), works in any dimension, and is the engine of adaptive filtering and machine learning. It is also slow on badly scaled problems, and understanding exactly why is the key to every faster method.

### Descent directions and the generic iteration

All methods in this section and the next have the form

$$
\mathbf{x}_{k+1}=\mathbf{x}_k+\alpha_k\mathbf{p}_k,
$$

with a **search direction** $\mathbf{p}_k$ and a **step size** (learning rate) $\alpha_k>0$.

> [!definition] Descent direction
> $\mathbf{p}$ is a descent direction for $f$ at $\mathbf{x}$ if $\nabla f(\mathbf{x})^\mathsf{T}\mathbf{p}<0$. Then $f(\mathbf{x}+\alpha\mathbf{p})<f(\mathbf{x})$ for all sufficiently small $\alpha>0$, because $\dfrac{d}{d\alpha}f(\mathbf{x}+\alpha\mathbf{p})\Big|_{\alpha=0}=\nabla f(\mathbf{x})^\mathsf{T}\mathbf{p}<0$.

Among unit vectors $\mathbf{p}$, the directional derivative $\nabla f^\mathsf{T}\mathbf{p}=\lVert\nabla f\rVert\cos\theta$ is most negative for $\theta=\pi$, i.e. $\mathbf{p}=-\nabla f/\lVert\nabla f\rVert$ ([Ch. 6 §6.5](../part-2-calculus/06-multivariable-calculus.md)). This gives **steepest descent** (gradient descent):

$$
\mathbf{x}_{k+1}=\mathbf{x}_k-\alpha_k\nabla f(\mathbf{x}_k).
$$

The negative gradient is perpendicular to the level curve through $\mathbf{x}_k$; each step leaves the contour at right angles. "Steepest" is relative to the Euclidean norm, i.e. to the units the variables are measured in — rescaling the variables changes the direction, which is the root of both the method's weakness and its cure (§23.6).

### Fixed step sizes: when does it converge?

Suppose $\nabla f$ is **Lipschitz continuous** with constant $L$: $\lVert\nabla f(\mathbf{x})-\nabla f(\mathbf{y})\rVert\le L\lVert\mathbf{x}-\mathbf{y}\rVert$ (for a $C^2$ function, $L$ bounds the largest Hessian eigenvalue in magnitude). Then $f$ lies below a quadratic upper bound,

$$
f(\mathbf{y})\le f(\mathbf{x})+\nabla f(\mathbf{x})^\mathsf{T}(\mathbf{y}-\mathbf{x})+\frac{L}{2}\lVert\mathbf{y}-\mathbf{x}\rVert^2 ,
$$

and substituting $\mathbf{y}=\mathbf{x}-\alpha\nabla f(\mathbf{x})$ yields the **descent lemma**

$$
f(\mathbf{x}_{k+1})\le f(\mathbf{x}_k)-\alpha\Big(1-\frac{L\alpha}{2}\Big)\lVert\nabla f(\mathbf{x}_k)\rVert^2 .
$$

> [!theorem] Convergence of gradient descent with a fixed step
> If $f$ is bounded below and $\nabla f$ is $L$-Lipschitz, then for any fixed $0<\alpha<2/L$ the iterates satisfy $f(\mathbf{x}_{k+1})<f(\mathbf{x}_k)$ whenever $\nabla f(\mathbf{x}_k)\ne\mathbf{0}$, and $\nabla f(\mathbf{x}_k)\to\mathbf{0}$. The choice $\alpha=1/L$ maximises the guaranteed decrease, $f(\mathbf{x}_{k+1})\le f(\mathbf{x}_k)-\frac1{2L}\lVert\nabla f(\mathbf{x}_k)\rVert^2$.

*Proof sketch.* The bracket $1-L\alpha/2$ is positive for $\alpha<2/L$, so $f$ decreases by a positive multiple of $\lVert\nabla f\rVert^2$. Summing the inequality from $k=0$ to $N-1$ gives $c\sum_k\lVert\nabla f(\mathbf{x}_k)\rVert^2\le f(\mathbf{x}_0)-f(\mathbf{x}_N)\le f(\mathbf{x}_0)-f_{\inf}<\infty$, so the gradient norms are square-summable and tend to zero. $\square$

Note what is *not* claimed: convergence to a global minimum. Gradient descent finds a stationary point; convexity (§23.7) is what makes that point the answer.

### The quadratic model problem

Everything about gradient descent can be seen on $f(\mathbf{x})=\tfrac12\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}-\mathbf{b}^\mathsf{T}\mathbf{x}$ with $\mathbf{A}$ symmetric positive definite, eigenvalues $0<m=\lambda_{\min}\le\dots\le\lambda_{\max}=L$. The minimiser solves $\mathbf{A}\mathbf{x}^*=\mathbf{b}$, and the gradient is $\mathbf{A}\mathbf{x}-\mathbf{b}=\mathbf{A}(\mathbf{x}-\mathbf{x}^*)$. The error $\mathbf{e}_k=\mathbf{x}_k-\mathbf{x}^*$ obeys

$$
\mathbf{e}_{k+1}=(\mathbf{I}-\alpha\mathbf{A})\,\mathbf{e}_k .
$$

Expanding $\mathbf{e}_0$ in the eigenvectors $\mathbf{v}_i$ of $\mathbf{A}$ ([Ch. 10](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md)), each component is multiplied by $1-\alpha\lambda_i$ per step. Hence:

- convergence for every start $\iff\lvert1-\alpha\lambda_i\rvert<1$ for all $i$ $\iff 0<\alpha<2/L$;
- the slowest component decays like $\max_i\lvert1-\alpha\lambda_i\rvert^k$;
- minimising $\max(\lvert1-\alpha m\rvert,\lvert1-\alpha L\rvert)$ gives the **best fixed step** $\alpha^*=\dfrac{2}{m+L}$ and the **convergence factor**
$$
\rho=\frac{L-m}{L+m}=\frac{\kappa-1}{\kappa+1},\qquad \kappa=\frac{L}{m}=\operatorname{cond}(\mathbf{A}).
$$

The **condition number** $\kappa$ governs everything. To gain one decimal digit of accuracy needs about $\ln10/\ln(1/\rho)\approx1.15\,\kappa$ iterations when $\kappa$ is large. A problem with $\kappa=10^4$ needs tens of thousands of gradient steps — while Newton's method (§23.6) would take one.

> [!example] Worked Example 23.15 — fixed step on a one-dimensional quadratic
> $f(x)=2x^2$, so $f'(x)=4x$, $L=4$, $x^*=0$, starting at $x_0=1$.
>
> *Step 1 (iteration):* $x_{k+1}=x_k-4\alpha x_k=(1-4\alpha)x_k$, so $x_k=(1-4\alpha)^k$.
>
> *Step 2 (cases):*
> - $\alpha=0.1$: factor $0.6$, monotone convergence; $x_{10}=0.006$.
> - $\alpha=0.25=1/L$: factor $0$, exact in one step (in one dimension $1/L$ is Newton's step).
> - $\alpha=0.4$: factor $-0.6$, convergence with sign flips (overshooting).
> - $\alpha=0.5=2/L$: factor $-1$, $x_k=(-1)^k$ oscillates forever.
> - $\alpha=0.6$: factor $-1.4$, divergence.
>
> *Step 3 (lesson):* the stability limit $\alpha<2/L$ is exactly the condition for a first-order discrete-time system $x_{k+1}=(1-\alpha L)x_k$ to have its pole inside the unit circle ([Ch. 17](../part-5-transforms/17-z-transforms.md)). Gradient descent *is* a discrete-time feedback loop.

> [!example] Worked Example 23.16 — step size on an anisotropic bowl
> $f(x,y)=x^2+5y^2$, i.e. $\mathbf{A}=\operatorname{diag}(2,10)$, $m=2$, $L=10$, $\kappa=5$. Start at $(5,1)$, $f_0=30$.
>
> *Step 1 (per-coordinate factors):* $x_{k+1}=(1-2\alpha)x_k$ and $y_{k+1}=(1-10\alpha)y_k$.
>
> *Step 2 (stability):* need $\lvert1-10\alpha\rvert<1$, i.e. $\alpha<0.2=2/L$. The *stiff* direction $y$ (large curvature) sets the limit, while the *soft* direction $x$ (small curvature) sets the speed.
>
> *Step 3 (compare after 20 steps):*
>
> | $\alpha$ | $x$ factor | $y$ factor | $f(\mathbf{x}_{20})$ |
> |---|---|---|---|
> | 0.05 | 0.9 | 0.5 | 0.37 |
> | 0.10 | 0.8 | 0 | $3.3\times10^{-3}$ |
> | 1/6 (optimal) | 2/3 | $-2/3$ | $2.7\times10^{-6}$ |
> | 0.18 | 0.64 | $-0.8$ | $6.7\times10^{-4}$ |
> | 0.21 | 0.58 | $-1.1$ | 226 (diverging) |
>
> For the optimal step both factors have magnitude $\rho=(5-1)/(5+1)=2/3$, so $f$ falls by $\rho^2=4/9$ per step: $30(4/9)^{20}=2.7\times10^{-6}$.
>
> *Step 4 (lesson):* the best step is a compromise dictated by the two extreme curvatures. If the problem were rescaled so that $\kappa=1$ (substitute $y=u/\sqrt5$), a single step would suffice.

### Line search

Rather than a fixed $\alpha$, choose $\alpha_k$ at each step by searching along the ray $\phi(\alpha)=f(\mathbf{x}_k+\alpha\mathbf{p}_k)$.

**Exact line search** minimises $\phi$. For the quadratic, $\phi(\alpha)=f(\mathbf{x}_k)-\alpha\mathbf{g}_k^\mathsf{T}\mathbf{g}_k+\tfrac12\alpha^2\mathbf{g}_k^\mathsf{T}\mathbf{A}\mathbf{g}_k$ along $\mathbf{p}_k=-\mathbf{g}_k$, so

$$
\alpha_k=\frac{\mathbf{g}_k^\mathsf{T}\mathbf{g}_k}{\mathbf{g}_k^\mathsf{T}\mathbf{A}\mathbf{g}_k}.
$$

At the exact minimiser along the line, $\phi'(\alpha_k)=\nabla f(\mathbf{x}_{k+1})^\mathsf{T}\mathbf{p}_k=0$: **the new gradient is orthogonal to the old search direction**. For steepest descent this means consecutive steps are at right angles — the notorious **zig-zag**. The classical bound (Kantorovich) is

$$
f(\mathbf{x}_{k+1})-f^*\le\Big(\frac{\kappa-1}{\kappa+1}\Big)^2\big(f(\mathbf{x}_k)-f^*\big),
$$

the same rate as the optimal fixed step: exact line search removes the need to know $m$ and $L$, but does not beat the condition number.

**Inexact line search** is what real codes do, because each evaluation of $f$ may be expensive. The aim is just "enough" decrease:

> [!definition] Armijo (sufficient decrease) and Wolfe conditions
> For a descent direction $\mathbf{p}_k$ and constants $0<c_1<c_2<1$ (typically $c_1=10^{-4}$, $c_2=0.9$), a step $\alpha$ satisfies
> - the **Armijo condition** if $f(\mathbf{x}_k+\alpha\mathbf{p}_k)\le f(\mathbf{x}_k)+c_1\alpha\,\nabla f(\mathbf{x}_k)^\mathsf{T}\mathbf{p}_k$;
> - the **curvature condition** if $\nabla f(\mathbf{x}_k+\alpha\mathbf{p}_k)^\mathsf{T}\mathbf{p}_k\ge c_2\,\nabla f(\mathbf{x}_k)^\mathsf{T}\mathbf{p}_k$.
>
> Both together are the **Wolfe conditions**. Armijo forbids steps that are too long (too little decrease for their length); the curvature condition forbids steps that are too short (the slope has not flattened enough).

**Backtracking** implements Armijo simply: try $\alpha=\alpha_0$ (say 1); while the Armijo condition fails, set $\alpha\leftarrow\beta\alpha$ (say $\beta=\tfrac12$). Because $\mathbf{p}_k$ is a descent direction, a small enough $\alpha$ always passes, so the loop terminates.

> [!example] Worked Example 23.17 — exact line search and the zig-zag
> Apply steepest descent with exact line search to $f=x^2+5y^2$ from $(5,1)$.
>
> *Step 1:* $\mathbf{g}_0=(2x,10y)=(10,10)$. $\mathbf{g}_0^\mathsf{T}\mathbf{g}_0=200$, $\mathbf{g}_0^\mathsf{T}\mathbf{A}\mathbf{g}_0=2(100)+10(100)=1200$, so $\alpha_0=\tfrac16$ and $\mathbf{x}_1=(5,1)-\tfrac16(10,10)=(\tfrac{10}3,-\tfrac23)$, $f_1=\tfrac{100}9+\tfrac{20}9=13.33$.
>
> *Step 2:* $\mathbf{g}_1=(\tfrac{20}3,-\tfrac{20}3)$. Check orthogonality: $\mathbf{g}_1^\mathsf{T}\mathbf{g}_0=\tfrac{200}3-\tfrac{200}3=0$ ✓. Again $\alpha_1=\tfrac16$, giving $\mathbf{x}_2=(\tfrac{20}9,\tfrac49)$, $f_2=5.93$.
>
> *Step 3 (pattern):* $\mathbf{x}_3=(\tfrac{40}{27},-\tfrac8{27})$, $\mathbf{x}_4=(\tfrac{80}{81},\tfrac{16}{81})$, and $f_k=30(\tfrac49)^k$: $30,\ 13.3,\ 5.93,\ 2.63,\ 1.17$. Each step multiplies $f$ by exactly $(\tfrac{\kappa-1}{\kappa+1})^2=\tfrac49$: this starting point is the worst case of the Kantorovich bound.
>
> *Step 4 (picture):* each step goes along $\pm45°$ and turns $90°$; the path bounces across the narrow valley.

```desmos-graph
left=-1; right=6; bottom=-1.5; top=2
---
x^2+5y^2=30|grey|dashed
x^2+5y^2=13.333|grey|dashed
x^2+5y^2=5.926|grey|dashed
x^2+5y^2=2.634|grey|dashed
x^2+5y^2=1.171|grey|dashed
y=x-4|red|3.3333<x<5
y=-x+2.6667|red|2.2222<x<3.3333
y=x-1.7778|red|1.4815<x<2.2222
y=-x+1.1852|red|0.98765<x<1.4815
(5,1)|label:x0
(3.3333,-0.6667)|label:x1
(2.2222,0.4444)|label:x2
(1.4815,-0.2963)|label:x3
(0,0)|label:minimum
```

*Steepest descent with exact line search on $f=x^2+5y^2$ (contours dashed). Every step leaves a contour at right angles and is orthogonal to the previous step, producing the characteristic zig-zag toward the minimum.*

> [!example] Worked Example 23.18 — backtracking on the Rosenbrock function
> The Rosenbrock function $f(x,y)=(1-x)^2+100(y-x^2)^2$ has its minimum $f=0$ at $(1,1)$ at the bottom of a long, curved, narrow valley. Take one backtracking step from the standard start $(-1.2,1)$ with $\alpha_0=1$, $\beta=\tfrac12$, $c_1=10^{-4}$.
>
> *Step 1 (gradient):* $\partial f/\partial x=-2(1-x)-400x(y-x^2)$ and $\partial f/\partial y=200(y-x^2)$. At $(-1.2,1)$: $y-x^2=-0.44$, so $\nabla f=(-4.4-211.2,\ -88)=(-215.6,\ -88)$, $\lVert\nabla f\rVert^2=54\,227$, $f_0=24.2$.
>
> *Step 2 (backtrack):* trial $f$ values for $\alpha=1,\tfrac12,\tfrac14,\dots$ are $2.1\times10^{11}$, $1.3\times10^{10}$, $\dots$, $543.5$ ($\alpha=2^{-6}$), $211.4$, $149.6$, $35.1$ ($\alpha=2^{-9}$) — all rejected; the Armijo bound $24.2-10^{-4}\alpha\,(54\,227)$ is about $24.2$.
>
> *Step 3 (accept):* $\alpha=2^{-10}\approx9.77\times10^{-4}$ gives $(-0.989,\ 1.086)$ with $f=5.10<24.2$ ✓.
>
> *Step 4 (the long run):* continuing with backtracking, gradient descent needs about $1.4\times10^4$ iterations to reach $\lVert\nabla f\rVert<10^{-6}$. The Hessian at the minimum is $\begin{pmatrix}802&-400\\-400&200\end{pmatrix}$, with condition number about 2500. Newton's method (§23.6) gets there in 6 iterations.

### Beyond plain gradient descent

- **Momentum (heavy ball).** $\mathbf{x}_{k+1}=\mathbf{x}_k-\alpha\nabla f(\mathbf{x}_k)+\beta(\mathbf{x}_k-\mathbf{x}_{k-1})$. The extra term is a low-pass filter on the step sequence that damps the zig-zag. On quadratics, optimally tuned $\alpha,\beta$ improve the convergence factor from $(\kappa-1)/(\kappa+1)$ to $(\sqrt\kappa-1)/(\sqrt\kappa+1)$, so the iteration count scales with $\sqrt\kappa$ instead of $\kappa$: for $\kappa=10^4$, roughly 100 times fewer iterations. Nesterov's accelerated gradient achieves this rate for general smooth convex functions.
- **Conjugate gradient (CG).** For quadratics, choosing each direction to be $\mathbf{A}$-conjugate to the previous ones ($\mathbf{p}_i^\mathsf{T}\mathbf{A}\mathbf{p}_j=0$) cures the zig-zag: CG reaches the exact minimiser in at most $n$ steps and is the standard iterative solver for large sparse SPD systems ([Ch. 24](24-computational-mathematics.md)). Nonlinear variants (Fletcher–Reeves, Polak–Ribière) extend it to general $f$.
- **Stochastic gradient.** When $f=\frac1N\sum_i f_i$ is an average over many data points, use the gradient of one (or a few) $f_i$ per step: noisy but $N$ times cheaper.

> [!tip] Engineering Note — the LMS adaptive filter is stochastic gradient descent
> An adaptive FIR filter with weights $\mathbf{w}$ and input vector $\mathbf{u}_n$ should minimise the mean-square error $E[e_n^2]$, $e_n=d_n-\mathbf{w}^\mathsf{T}\mathbf{u}_n$. The true gradient is $-2E[e_n\mathbf{u}_n]$; replacing the expectation by its instantaneous value gives Widrow's **least-mean-squares** update $\mathbf{w}_{n+1}=\mathbf{w}_n+\mu\,e_n\mathbf{u}_n$. The stability limit $\alpha<2/L$ becomes the classical LMS bound $0<\mu<2/\lambda_{\max}(\mathbf{R})$, where $\mathbf{R}=E[\mathbf{u}_n\mathbf{u}_n^\mathsf{T}]$ is the input autocorrelation matrix ([Ch. 21](../part-6-probability-statistics/21-random-processes.md)), and the convergence speed is governed by its eigenvalue spread — exactly the condition number story above.

> [!warning] Common Pitfall — stopping tests
> "Stop when $f$ stops changing" can halt in a flat region far from a minimum; "stop when $\lVert\nabla f\rVert<\varepsilon$" depends on the units of $f$ and $\mathbf{x}$. Robust codes combine a *relative* gradient test, $\lVert\nabla f(\mathbf{x}_k)\rVert\le\varepsilon\max(1,\lvert f(\mathbf{x}_k)\rvert)$, with a small-step test and an iteration limit, and report which one fired.

A minimal implementation:

```python
import numpy as np

def grad_descent(f, grad, x0, alpha0=1.0, beta=0.5, c1=1e-4, tol=1e-6, maxit=100000):
    x = np.asarray(x0, float)
    for k in range(maxit):
        g = grad(x)
        if np.linalg.norm(g) < tol:
            break
        a = alpha0
        while f(x - a*g) > f(x) - c1*a*(g @ g):   # Armijo backtracking
            a *= beta
        x = x - a*g
    return x, k
```

## 23.6 Newton and quasi-Newton methods

### Motivation

Gradient descent uses only the slope; it does not know how curved the function is, so it cannot tell a gentle valley from a steep one and must use a step that is safe for the steepest direction. Newton's method also uses the **curvature** (the Hessian): it fits a quadratic bowl to $f$ at the current point and jumps straight to the bottom of the bowl. Near a minimum this converges *quadratically* — the number of correct digits doubles every iteration — and is immune to bad scaling. The price is computing and factoring an $n\times n$ Hessian at every step. **Quasi-Newton** methods, BFGS above all, get most of the benefit while using only gradients, by learning the curvature from the history of steps. BFGS is the default unconstrained method in most scientific software, including SciPy's `minimize`.

### Newton's method

At $\mathbf{x}_k$, with $\mathbf{g}_k=\nabla f(\mathbf{x}_k)$ and $\mathbf{H}_k=\nabla^2f(\mathbf{x}_k)$, the quadratic model is

$$
m_k(\mathbf{p})=f(\mathbf{x}_k)+\mathbf{g}_k^\mathsf{T}\mathbf{p}+\tfrac12\mathbf{p}^\mathsf{T}\mathbf{H}_k\mathbf{p}.
$$

If $\mathbf{H}_k\succ0$, the model has a unique minimiser where $\nabla m_k=\mathbf{g}_k+\mathbf{H}_k\mathbf{p}=\mathbf{0}$:

> [!definition] Newton step
> $$
> \mathbf{H}_k\,\mathbf{p}_k^{\mathrm N}=-\mathbf{g}_k,\qquad \mathbf{x}_{k+1}=\mathbf{x}_k+\alpha_k\mathbf{p}_k^{\mathrm N},
> $$
> with $\alpha_k=1$ for the *pure* Newton method. The linear system is solved (never by forming $\mathbf{H}_k^{-1}$) by Cholesky factorisation, which also tests $\mathbf{H}_k\succ0$ for free: the factorisation fails if it is not ([Ch. 10 §10.6](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md)).

Two other views of the same iteration are useful:

- **Root-finding view.** Minimising $f$ means solving $\nabla f(\mathbf{x})=\mathbf{0}$. Newton's method for the nonlinear system $\mathbf{F}(\mathbf{x})=\mathbf{0}$ uses the Jacobian of $\mathbf{F}$ ([Ch. 6 §6.6](../part-2-calculus/06-multivariable-calculus.md), [Ch. 22](22-numerical-methods.md)); with $\mathbf{F}=\nabla f$ the Jacobian is the Hessian, and we recover the Newton step. In one dimension: $x_{k+1}=x_k-f'(x_k)/f''(x_k)$.
- **Rescaling view.** If $\mathbf{H}\succ0$, change variables to $\mathbf{u}=\mathbf{H}^{1/2}\mathbf{x}$. In $\mathbf{u}$ the local quadratic model is a perfectly round bowl ($\kappa=1$) and steepest descent with unit step goes straight to its centre. Translating back gives the Newton step. Newton's method is **affine invariant**: a linear change of variables $\mathbf{x}=\mathbf{T}\mathbf{z}$ does not change the iterates (mapped back), so units and scaling do not matter — in sharp contrast to gradient descent.

A **descent property**: if $\mathbf{H}_k\succ0$ then $\mathbf{g}_k^\mathsf{T}\mathbf{p}^{\mathrm N}_k=-\mathbf{g}_k^\mathsf{T}\mathbf{H}_k^{-1}\mathbf{g}_k<0$, so the Newton direction is a descent direction and a line search can be applied. If $\mathbf{H}_k$ is indefinite, the "Newton step" heads for the stationary point of an indefinite quadratic — a saddle — and may go uphill.

> [!theorem] Local quadratic convergence of Newton's method
> Let $f\in C^2$ with Lipschitz-continuous Hessian near a point $\mathbf{x}^*$ where $\nabla f(\mathbf{x}^*)=\mathbf{0}$ and $\nabla^2f(\mathbf{x}^*)\succ0$. Then for $\mathbf{x}_0$ close enough to $\mathbf{x}^*$, pure Newton iterates are well defined, converge to $\mathbf{x}^*$, and satisfy
> $$
> \lVert\mathbf{x}_{k+1}-\mathbf{x}^*\rVert\le C\,\lVert\mathbf{x}_k-\mathbf{x}^*\rVert^2
> $$
> for a constant $C$ (roughly $\tfrac12\times$ Hessian Lipschitz constant $\times\lVert\nabla^2f(\mathbf{x}^*)^{-1}\rVert$).

*Proof sketch (one dimension).* Write $e_k=x_k-x^*$. Taylor's theorem gives $0=f'(x^*)=f'(x_k)-f''(x_k)e_k+\tfrac12f'''(\xi)e_k^2$. Divide by $f''(x_k)$ and rearrange: $e_k-f'(x_k)/f''(x_k)=\tfrac{f'''(\xi)}{2f''(x_k)}e_k^2$. The left side is $e_{k+1}$. Hence $e_{k+1}\approx\dfrac{f'''(x^*)}{2f''(x^*)}e_k^2$. $\square$

Quadratic convergence is spectacular *once you are close*; far away the pure method may diverge, cycle or head for a maximum. The theorem says nothing about how close is close enough.

> [!example] Worked Example 23.19 — quadratic convergence in one dimension
> Minimise $f(x)=e^x-2x$ (minimiser $x^*=\ln2=0.693147$) from $x_0=0$.
>
> *Step 1 (iteration):* $f'(x)=e^x-2$, $f''(x)=e^x$, so $x_{k+1}=x_k-(e^{x_k}-2)/e^{x_k}=x_k-1+2e^{-x_k}$.
>
> *Step 2 (iterates):*
>
> | $k$ | $x_k$ | error $e_k=x_k-\ln2$ | $e_k/e_{k-1}^2$ |
> |---|---|---|---|
> | 0 | 0 | $-0.693$ | — |
> | 1 | 1.000000 | $3.07\times10^{-1}$ | 0.64 |
> | 2 | 0.735759 | $4.26\times10^{-2}$ | 0.45 |
> | 3 | 0.694042 | $8.95\times10^{-4}$ | 0.49 |
> | 4 | 0.6931476 | $4.00\times10^{-7}$ | 0.50 |
> | 5 | 0.69314718056 | $8.0\times10^{-14}$ | 0.50 |
>
> *Step 3 (check the constant):* the theory predicts $e_{k+1}/e_k^2\to f'''(x^*)/(2f''(x^*))=e^{x^*}/(2e^{x^*})=\tfrac12$ ✓. The number of correct digits goes 1, 1, 3, 6, 13: doubling each step.

> [!example] Worked Example 23.20 — one Newton step solves a quadratic
> For $f(\mathbf{x})=\tfrac12\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}-\mathbf{b}^\mathsf{T}\mathbf{x}$ with $\mathbf{A}\succ0$, the Hessian is $\mathbf{A}$ everywhere and $\mathbf{p}^{\mathrm N}=-\mathbf{A}^{-1}(\mathbf{A}\mathbf{x}_0-\mathbf{b})=\mathbf{x}^*-\mathbf{x}_0$. One step lands exactly on $\mathbf{x}^*$ from any start, whatever $\kappa$. On $f=x^2+5y^2$ from $(5,1)$: $\mathbf{g}=(10,10)$, $\mathbf{H}=\operatorname{diag}(2,10)$, $\mathbf{p}=-(5,1)$, $\mathbf{x}_1=(0,0)$ ✓ — compared with the endless zig-zag of Example 23.17.

> [!example] Worked Example 23.21 — pure Newton on the Rosenbrock function
> From $(-1.2,1)$, using the Hessian
> $$
> \nabla^2f=\begin{pmatrix}2-400(y-x^2)+800x^2&-400x\\-400x&200\end{pmatrix}:
> $$
>
> | $k$ | $(x_k,y_k)$ | $f(\mathbf{x}_k)$ |
> |---|---|---|
> | 0 | $(-1.2000,\ 1.0000)$ | 24.2 |
> | 1 | $(-1.1753,\ 1.3807)$ | 4.73 |
> | 2 | $(0.7631,\ -3.1750)$ | **1411.8** |
> | 3 | $(0.7634,\ 0.5828)$ | 0.0560 |
> | 4 | $(0.999995,\ 0.9440)$ | 0.313 |
> | 5 | $(0.9999957,\ 0.9999914)$ | $1.9\times10^{-11}$ |
> | 6 | $(1,\ 1)$ | $3\times10^{-20}$ |
>
> *Observations:* (i) six iterations versus about fourteen thousand for gradient descent (Example 23.18); (ii) the pure method is not a descent method — at step 2 the objective jumps from 4.7 to 1412 because the quadratic model is trusted far outside its region of validity; it happens to recover here, but on other problems such a jump diverges. (iii) The last two steps show the quadratic convergence. A line search (damped Newton) would have rejected the full step at iteration 2 and kept $f$ monotone.

### Making Newton's method robust

- **Damped Newton.** Use the Newton direction with a backtracking line search starting from $\alpha=1$ (§23.5). Far from the solution, short steps keep $f$ decreasing; near the solution $\alpha=1$ is accepted and quadratic convergence is retained.
- **Hessian modification.** If $\mathbf{H}_k$ is not positive definite, replace it by $\mathbf{H}_k+\tau\mathbf{I}$ with $\tau>0$ large enough that Cholesky succeeds. For small $\tau$ this is nearly Newton; for large $\tau$ it is a short gradient step $-\mathbf{g}/\tau$.
- **Trust region.** Minimise the model $m_k(\mathbf{p})$ only over $\lVert\mathbf{p}\rVert\le\Delta_k$, then compare the actual decrease with the predicted decrease. If they agree, accept and enlarge $\Delta$; if not, shrink $\Delta$. The KKT conditions of this subproblem give $(\mathbf{H}_k+\tau\mathbf{I})\mathbf{p}=-\mathbf{g}_k$ with $\tau\ge0$ the multiplier of the radius constraint — so Hessian modification and trust regions are the same idea. For least squares this is the **Levenberg–Marquardt** method (§23.8).

**Cost.** Each Newton step needs the Hessian ($n(n+1)/2$ second derivatives) and a Cholesky factorisation ($\approx n^3/3$ flops). For $n$ in the hundreds this is cheap; for $n=10^5$ it is prohibitive unless the Hessian is sparse — which motivates quasi-Newton and, for very large problems, limited-memory variants.

### Quasi-Newton methods and BFGS

The idea: maintain an approximation $\mathbf{B}_k\approx\nabla^2f(\mathbf{x}_k)$ (or its inverse $\mathbf{H}_k\approx\nabla^2f(\mathbf{x}_k)^{-1}$), take steps $\mathbf{p}_k=-\mathbf{H}_k\mathbf{g}_k$, and after each step update the approximation using only the gradient change.

With $\mathbf{s}_k=\mathbf{x}_{k+1}-\mathbf{x}_k$ and $\mathbf{y}_k=\mathbf{g}_{k+1}-\mathbf{g}_k$, a first-order expansion of the gradient gives $\mathbf{y}_k\approx\nabla^2f\,\mathbf{s}_k$ — exactly for a quadratic. So we require the new approximation to satisfy the

> [!definition] Secant (quasi-Newton) condition
> $$
> \mathbf{B}_{k+1}\mathbf{s}_k=\mathbf{y}_k\qquad\Longleftrightarrow\qquad \mathbf{H}_{k+1}\mathbf{y}_k=\mathbf{s}_k .
> $$

This is $n$ equations for the $n(n+1)/2$ entries of a symmetric matrix, so there is freedom. The Broyden–Fletcher–Goldfarb–Shanno (BFGS) update chooses the symmetric $\mathbf{H}_{k+1}$ satisfying the secant condition that is *closest* to $\mathbf{H}_k$ in a suitably weighted Frobenius norm. The result is a rank-two correction:

> [!theorem] BFGS update
> With $\rho_k=1/(\mathbf{y}_k^\mathsf{T}\mathbf{s}_k)$,
> $$
> \mathbf{H}_{k+1}=\big(\mathbf{I}-\rho_k\mathbf{s}_k\mathbf{y}_k^\mathsf{T}\big)\,\mathbf{H}_k\,\big(\mathbf{I}-\rho_k\mathbf{y}_k\mathbf{s}_k^\mathsf{T}\big)+\rho_k\mathbf{s}_k\mathbf{s}_k^\mathsf{T}.
> $$
> Properties:
> 1. $\mathbf{H}_{k+1}$ is symmetric and satisfies $\mathbf{H}_{k+1}\mathbf{y}_k=\mathbf{s}_k$;
> 2. if $\mathbf{H}_k\succ0$ and the **curvature condition** $\mathbf{y}_k^\mathsf{T}\mathbf{s}_k>0$ holds, then $\mathbf{H}_{k+1}\succ0$ — so every BFGS direction is a descent direction;
> 3. on a strictly convex quadratic with exact line searches, BFGS terminates in at most $n$ steps with $\mathbf{H}_n=\mathbf{A}^{-1}$;
> 4. near a minimiser with $\nabla^2f\succ0$, BFGS with a Wolfe line search converges **superlinearly**: $\lVert\mathbf{x}_{k+1}-\mathbf{x}^*\rVert/\lVert\mathbf{x}_k-\mathbf{x}^*\rVert\to0$.

*Why property 1 holds.* Multiply by $\mathbf{y}_k$: the right factor gives $(\mathbf{I}-\rho_k\mathbf{y}_k\mathbf{s}_k^\mathsf{T})\mathbf{y}_k=\mathbf{y}_k-\rho_k\mathbf{y}_k(\mathbf{s}_k^\mathsf{T}\mathbf{y}_k)=\mathbf{0}$, so only the last term survives: $\rho_k\mathbf{s}_k(\mathbf{s}_k^\mathsf{T}\mathbf{y}_k)=\mathbf{s}_k$ ✓. *Why property 2 holds.* For any $\mathbf{z}\ne\mathbf{0}$, put $\mathbf{w}=\mathbf{z}-\rho_k\mathbf{y}_k(\mathbf{s}_k^\mathsf{T}\mathbf{z})$; then $\mathbf{z}^\mathsf{T}\mathbf{H}_{k+1}\mathbf{z}=\mathbf{w}^\mathsf{T}\mathbf{H}_k\mathbf{w}+\rho_k(\mathbf{s}_k^\mathsf{T}\mathbf{z})^2$. Both terms are $\ge0$, and they cannot both vanish (if $\mathbf{s}_k^\mathsf{T}\mathbf{z}=0$ then $\mathbf{w}=\mathbf{z}\neq\mathbf{0}$). The Wolfe curvature condition of §23.5 guarantees $\mathbf{y}_k^\mathsf{T}\mathbf{s}_k>0$, which is why BFGS is always paired with a Wolfe line search.

Related methods: **DFP** is the same formula with the roles of $\mathbf{s}$ and $\mathbf{y}$ (and $\mathbf{B}$ and $\mathbf{H}$) swapped, and is less robust in practice; **SR1** is a symmetric rank-one update that need not stay positive definite (useful inside trust regions); **L-BFGS** stores only the last $m\approx5$–$20$ pairs $(\mathbf{s}_i,\mathbf{y}_i)$ and applies $\mathbf{H}_k$ implicitly, using $O(mn)$ memory — the method of choice for large smooth problems.

> [!example] Worked Example 23.22 — BFGS by hand on a quadratic
> Minimise $f=x^2+5y^2$ ($\mathbf{A}=\operatorname{diag}(2,10)$) from $\mathbf{x}_0=(5,1)$ with $\mathbf{H}_0=\mathbf{I}$ and exact line search.
>
> *Step 1 (first step = steepest descent):* $\mathbf{g}_0=(10,10)$, $\mathbf{p}_0=-(10,10)$, $\alpha_0=\tfrac16$ (as in Example 23.17), $\mathbf{x}_1=(\tfrac{10}3,-\tfrac23)$, $\mathbf{g}_1=(\tfrac{20}3,-\tfrac{20}3)$.
>
> *Step 2 (update data):* $\mathbf{s}_0=(-\tfrac53,-\tfrac53)$, $\mathbf{y}_0=\mathbf{g}_1-\mathbf{g}_0=(-\tfrac{10}3,-\tfrac{50}3)$, $\mathbf{y}_0^\mathsf{T}\mathbf{s}_0=\tfrac{50}9+\tfrac{250}9=\tfrac{100}3>0$ ✓, $\rho_0=0.03$.
>
> *Step 3 (BFGS update):* evaluating the formula gives
> $$
> \mathbf{H}_1=\frac1{36}\begin{pmatrix}53&-7\\-7&5\end{pmatrix}.
> $$
> Secant check: $\mathbf{H}_1\mathbf{y}_0=\tfrac1{36}\big(53(-\tfrac{10}3)+(-7)(-\tfrac{50}3),\ (-7)(-\tfrac{10}3)+5(-\tfrac{50}3)\big)=\tfrac1{36}(-60,-60)=(-\tfrac53,-\tfrac53)=\mathbf{s}_0$ ✓. Positive definite: trace $>0$, determinant $(265-49)/36^2>0$ ✓.
>
> *Step 4 (second step):* $\mathbf{p}_1=-\mathbf{H}_1\mathbf{g}_1=-\tfrac{20}3\cdot\tfrac1{36}(53+7,\ -7-5)=(-\tfrac{100}9,\ \tfrac{20}9)$. Exact line search: $\alpha_1=-\mathbf{g}_1^\mathsf{T}\mathbf{p}_1/(\mathbf{p}_1^\mathsf{T}\mathbf{A}\mathbf{p}_1)=0.3$, so $\mathbf{x}_2=(\tfrac{10}3-\tfrac{10}3,\ -\tfrac23+\tfrac23)=(0,0)$ — **the exact minimiser after $n=2$ steps** ✓.
>
> *Step 5 (learned curvature):* one more update gives $\mathbf{H}_2=\operatorname{diag}(0.5,\ 0.1)=\mathbf{A}^{-1}$ exactly: BFGS has learned the inverse Hessian from two gradient differences.

### Comparing the methods

On the Rosenbrock function from $(-1.2,1)$ (tolerance about $10^{-5}$ on the gradient):

| Method | Information used | Iterations | Function evaluations | Cost per iteration |
|---|---|---|---|---|
| Gradient descent + backtracking | $f,\nabla f$ | $\approx14\,000$ | many more | $O(n)$ |
| Nelder–Mead simplex (derivative-free) | $f$ | 85 | 159 | $O(n)$ |
| Nonlinear conjugate gradient | $f,\nabla f$ | 36 | 78 | $O(n)$ |
| BFGS | $f,\nabla f$ | 32 | 39 | $O(n^2)$ |
| Newton (pure) | $f,\nabla f,\nabla^2f$ | 6 | 6 | $O(n^3)$ |

(SciPy 1.18 defaults; exact counts depend on implementation and tolerances.) The rule of thumb: use Newton when Hessians are cheap and $n$ is modest; BFGS when only gradients are available; L-BFGS or CG when $n$ is large; derivative-free methods only when gradients are unavailable and $n$ is small.

> [!warning] Common Pitfall — wrong gradients
> Quasi-Newton methods trust the gradient you supply. A sign error or a missing factor of 2 produces failed line searches, "precision loss" warnings, or convergence to a non-stationary point. Always check an analytic gradient against central finite differences, $[f(\mathbf{x}+h\mathbf{e}_i)-f(\mathbf{x}-h\mathbf{e}_i)]/(2h)$ with $h\approx10^{-6}\lvert x_i\rvert$, before optimising (`scipy.optimize.check_grad` does this).

> [!tip] Engineering Note — circuit simulators run Newton's method
> SPICE's DC operating-point analysis solves the nonlinear nodal equations with Newton–Raphson and the device Jacobian; its convergence aids ("gmin stepping", "source stepping", step limiting on junction voltages) are exactly the damping and continuation ideas above. When a simulation reports "no convergence", the remedies are the optimiser's remedies: a better starting point, smaller steps, or a regularised Jacobian.

```python
import numpy as np
from scipy.optimize import minimize

def f(v):    x, y = v; return (1 - x)**2 + 100*(y - x**2)**2
def grad(v): x, y = v; return np.array([-2*(1 - x) - 400*x*(y - x**2), 200*(y - x**2)])

res = minimize(f, x0=[-1.2, 1.0], jac=grad, method="BFGS")
print(res.x, res.nit, res.nfev)     # approx [1. 1.], 32, 39
```

## 23.7 Convex optimization fundamentals

### Motivation

Sections 23.5–23.6 delivered algorithms that find *local* minima — points where no small move helps. For a general nonlinear problem that is all anyone can promise: a local method has no way to know whether a better valley exists somewhere else. There is, however, a large and very useful class of problems for which this worry disappears entirely. For **convex** problems every local minimum is global, the KKT conditions are sufficient as well as necessary, duality gives a certificate of optimality, and reliable algorithms find the solution in polynomial time. Least squares, linear programming, minimax and $L^1$ fitting, many filter-design problems, power allocation, economic dispatch with convex costs, and (after a change of variables) many component-sizing problems are all convex. The practical motto, due to Boyd and Vandenberghe, is that the great watershed in optimization is not between linear and nonlinear problems but between **convex and nonconvex** ones. Recognising convexity — or reformulating a problem until it is convex — is the most valuable modelling skill in this chapter.

### Convex sets

> [!definition] Convex set
> A set $C\subseteq\mathbb{R}^n$ is **convex** if for all $\mathbf{x},\mathbf{y}\in C$ and all $\theta\in[0,1]$,
> $$
> \theta\mathbf{x}+(1-\theta)\mathbf{y}\in C .
> $$
> That is, the line segment joining any two points of $C$ lies in $C$.

Standard examples:
- **hyperplanes** $\{\mathbf{x}:\mathbf{a}^\mathsf{T}\mathbf{x}=b\}$ and **half-spaces** $\{\mathbf{x}:\mathbf{a}^\mathsf{T}\mathbf{x}\le b\}$;
- **polyhedra** $\{\mathbf{x}:\mathbf{A}\mathbf{x}\le\mathbf{b},\ \mathbf{C}\mathbf{x}=\mathbf{d}\}$ (LP feasible sets);
- **norm balls** $\{\mathbf{x}:\lVert\mathbf{x}-\mathbf{x}_c\rVert\le r\}$ for any norm, and **ellipsoids** $\{\mathbf{x}:(\mathbf{x}-\mathbf{x}_c)^\mathsf{T}\mathbf{P}^{-1}(\mathbf{x}-\mathbf{x}_c)\le1\}$, $\mathbf{P}\succ0$;
- the **positive semidefinite cone** $\{\mathbf{X}=\mathbf{X}^\mathsf{T}:\mathbf{X}\succeq0\}$;
- sublevel sets $\{\mathbf{x}: f(\mathbf{x})\le c\}$ of a convex function $f$ (shown below).

Non-examples: an annulus $1\le x^2+y^2\le4$ (the segment between opposite points passes through the hole); the set $\{(x,y): xy\ge1\}$ including both branches; any finite set of two or more points (e.g. E12 resistor values — which is why discrete design is hard).

> [!theorem] Operations that preserve convexity of sets
> Intersections (of any number) of convex sets are convex. Images and inverse images of convex sets under affine maps $\mathbf{x}\mapsto\mathbf{T}\mathbf{x}+\mathbf{t}$ are convex.

*Proof for intersections.* If $\mathbf{x},\mathbf{y}$ lie in every $C_i$, so does each point of the segment between them, since each $C_i$ is convex. $\square$

This is why a feasible set defined by many convex constraints — $g_i(\mathbf{x})\le0$ with each $g_i$ convex, and affine equalities — is convex: it is an intersection of convex sublevel sets and hyperplanes. Note that a *nonlinear equality* $h(\mathbf{x})=0$ generally destroys convexity (a circle $x^2+y^2=1$ is not convex, though the disc $x^2+y^2\le1$ is).

### Convex functions

> [!definition] Convex function
> $f:C\to\mathbb{R}$ on a convex set $C$ is **convex** if for all $\mathbf{x},\mathbf{y}\in C$ and $\theta\in[0,1]$,
> $$
> f\big(\theta\mathbf{x}+(1-\theta)\mathbf{y}\big)\le\theta f(\mathbf{x})+(1-\theta)f(\mathbf{y}),
> $$
> i.e. the chord lies on or above the graph. It is **strictly convex** if the inequality is strict for $\mathbf{x}\ne\mathbf{y}$, $\theta\in(0,1)$, and **strongly convex** with parameter $m>0$ if $f(\mathbf{x})-\tfrac m2\lVert\mathbf{x}\rVert^2$ is convex. $f$ is **concave** if $-f$ is convex.

Equivalently, $f$ is convex iff its **epigraph** $\{(\mathbf{x},t): f(\mathbf{x})\le t\}$ is a convex set — this is how the epigraph trick of §23.2 preserves convexity. The defining inequality extends by induction to **Jensen's inequality**: $f\big(\sum\theta_i\mathbf{x}_i\big)\le\sum\theta_if(\mathbf{x}_i)$ for weights $\theta_i\ge0$ summing to 1, and in probabilistic form $f(E[X])\le E[f(X)]$ ([Ch. 19](../part-6-probability-statistics/19-probability.md)).

```desmos-graph
left=-2.5; right=2.8; bottom=-1.5; top=5
---
y=x^2|blue
y=x+2|blue|dashed|-1<x<2
(-1,1)|label:x
(2,4)|label:y
y=0.5x^4-1.5x^2+0.3x+1|red
y=0.148x+0.194|red|dotted|-1.6<x<1.4
(-1.6,-0.043)
(1.4,0.401)
```

*Blue: $f=x^2$ is convex — every chord (dashed) lies above the graph. Red: a quartic with two local minima is not convex — the dotted chord between points near its two wells passes below the hump in the middle.*

For differentiable functions there are two much more practical tests.

> [!theorem] First- and second-order characterisations
> Let $C$ be open and convex.
> 1. **First order.** A differentiable $f$ is convex on $C$ iff for all $\mathbf{x},\mathbf{y}\in C$
> $$
> f(\mathbf{y})\ge f(\mathbf{x})+\nabla f(\mathbf{x})^\mathsf{T}(\mathbf{y}-\mathbf{x}),
> $$
> i.e. every tangent plane is a global under-estimator.
> 2. **Second order.** A twice-differentiable $f$ is convex on $C$ iff $\nabla^2f(\mathbf{x})\succeq0$ for all $\mathbf{x}\in C$. If $\nabla^2f(\mathbf{x})\succ0$ everywhere, $f$ is strictly convex; if $\nabla^2f(\mathbf{x})\succeq m\mathbf{I}$, it is $m$-strongly convex.

*Proof sketch.* (1, $\Rightarrow$) From the definition, $f(\mathbf{x}+\theta(\mathbf{y}-\mathbf{x}))\le f(\mathbf{x})+\theta(f(\mathbf{y})-f(\mathbf{x}))$; rearrange to $\frac{f(\mathbf{x}+\theta(\mathbf{y}-\mathbf{x}))-f(\mathbf{x})}{\theta}\le f(\mathbf{y})-f(\mathbf{x})$ and let $\theta\to0^+$: the left side tends to $\nabla f(\mathbf{x})^\mathsf{T}(\mathbf{y}-\mathbf{x})$. (1, $\Leftarrow$) Apply the inequality at $\mathbf{z}=\theta\mathbf{x}+(1-\theta)\mathbf{y}$ toward both $\mathbf{x}$ and $\mathbf{y}$, multiply by $\theta$ and $1-\theta$, add; the gradient terms cancel. (2) Restrict $f$ to a line, $\varphi(t)=f(\mathbf{x}+t\mathbf{d})$; $f$ is convex iff every such $\varphi$ is, and a one-variable $C^2$ function is convex iff $\varphi''(t)=\mathbf{d}^\mathsf{T}\nabla^2f\,\mathbf{d}\ge0$. $\square$

The first-order condition is the heart of the matter. If $\nabla f(\mathbf{x}^*)=\mathbf{0}$, it says $f(\mathbf{y})\ge f(\mathbf{x}^*)$ for **every** $\mathbf{y}$: a stationary point of a convex function is a global minimiser. Local information certifies a global fact.

### A toolkit for recognising convex functions

Building blocks (convex on the stated domain):

| Function | Domain | Note |
|---|---|---|
| $\mathbf{a}^\mathsf{T}\mathbf{x}+b$ | $\mathbb{R}^n$ | affine: both convex and concave |
| $\tfrac12\mathbf{x}^\mathsf{T}\mathbf{Q}\mathbf{x}$ | $\mathbb{R}^n$ | iff $\mathbf{Q}\succeq0$ |
| $\lVert\mathbf{x}\rVert$ (any norm) | $\mathbb{R}^n$ | triangle inequality |
| $e^{ax}$ | $\mathbb{R}$ | |
| $x^p$ | $x>0$ | for $p\ge1$ or $p\le0$; concave for $0\le p\le1$ |
| $-\ln x$, $x\ln x$, $1/x$ | $x>0$ | |
| $\log\sum_i e^{x_i}$ | $\mathbb{R}^n$ | smooth approximation of $\max_i x_i$ |
| $\max_i\,(\mathbf{a}_i^\mathsf{T}\mathbf{x}+b_i)$ | $\mathbb{R}^n$ | piecewise-linear |

Operations that preserve convexity of functions:
1. nonnegative weighted sums $\sum w_if_i$, $w_i\ge0$;
2. composition with an affine map, $f(\mathbf{A}\mathbf{x}+\mathbf{b})$ — so $\lVert\mathbf{A}\mathbf{x}-\mathbf{b}\rVert$ and $\lVert\mathbf{A}\mathbf{x}-\mathbf{b}\rVert^2$ are convex;
3. pointwise maximum (or supremum) of convex functions — so $\max_i\lvert r_i(\mathbf{x})\rvert$ with affine $r_i$ is convex;
4. composition $h(g(\mathbf{x}))$ with $h$ convex and nondecreasing and $g$ convex (e.g. $e^{g(\mathbf{x})}$).

A product or ratio of convex functions is in general *not* convex ($xy$ is not; $x^2/y$ for $y>0$ happens to be).

> [!example] Worked Example 23.23 — testing convexity with the Hessian
> Decide which are convex on $\mathbb{R}^2$: (a) $f=x^2+xy+y^2$; (b) $f=xy$; (c) $f=x^2+4xy+y^2$; (d) $f=e^{x+y}+x^2$.
>
> *(a)* $\nabla^2f=\begin{pmatrix}2&1\\1&2\end{pmatrix}$, eigenvalues $1,3>0$: strictly (indeed 1-strongly) convex. ✓
>
> *(b)* $\nabla^2f=\begin{pmatrix}0&1\\1&0\end{pmatrix}$, eigenvalues $\pm1$: indefinite, not convex (a saddle surface). Along $y=-x$, $f=-x^2$ is concave.
>
> *(c)* $\nabla^2f=\begin{pmatrix}2&4\\4&2\end{pmatrix}$, eigenvalues $6$ and $-2$: not convex, even though both "diagonal" terms are convex. Along $y=-x$: $f=x^2-4x^2+x^2=-2x^2$.
>
> *(d)* $\nabla^2f=e^{x+y}\begin{pmatrix}1&1\\1&1\end{pmatrix}+\begin{pmatrix}2&0\\0&0\end{pmatrix}$: a sum of two PSD matrices, so PSD — convex. (Or by the toolkit: $e^{u}$ composed with the affine $u=x+y$, plus $x^2$.) Neither term alone is positive definite, but the sum has determinant $(e^{x+y}+2)e^{x+y}-e^{2(x+y)}=2e^{x+y}>0$ and positive trace, so it is positive definite: $f$ is strictly convex.

> [!example] Worked Example 23.24 — convexity in engineering variables
> (a) The conductor-cost function $f(A)=c_1A+k/A$ of Example 23.2 has $f''(A)=2k/A^3>0$ on $A>0$: strictly convex, so the stationary point found there is the unique global optimum.
>
> (b) The set of capacitances meeting a $kT/C$ noise limit, $\{C>0:\sqrt{kT/C}\le v_{\max}\}=\{C\ge kT/v_{\max}^2\}$, is a half-line: convex.
>
> (c) The power-transfer objective $P(R_L)=V_s^2R_L/(R_s+R_L)^2$ of Example 23.1 is *not* concave on $(0,\infty)$: $P''$ changes sign at $R_L=2R_s$ (check: $P''\propto 2(R_L-2R_s)/(R_s+R_L)^4$). It is, however, **unimodal**, so the stationary point is still the global maximum. Convexity is sufficient for "local = global", not necessary.
>
> (d) The set $\{(x,y): x>0,\ y>0,\ xy\ge1\}$ (the region above one branch of a hyperbola) is convex: it is $\{-\ln x-\ln y\le0\}$, a sublevel set of a convex function.

### Convex optimization problems

> [!definition] Convex optimization problem
> $$
> \min\ f(\mathbf{x})\quad\text{s.t.}\quad g_i(\mathbf{x})\le0\ (i=1,\dots,m),\qquad \mathbf{A}\mathbf{x}=\mathbf{b},
> $$
> with $f$ and every $g_i$ convex. (Equality constraints must be affine.)

> [!theorem] Local equals global
> For a convex optimization problem:
> 1. every local minimiser is a global minimiser;
> 2. the set of global minimisers is convex (a point, a segment, a face…);
> 3. if $f$ is strictly convex, there is at most one minimiser;
> 4. a feasible $\mathbf{x}^*$ is optimal iff $\nabla f(\mathbf{x}^*)^\mathsf{T}(\mathbf{y}-\mathbf{x}^*)\ge0$ for all feasible $\mathbf{y}$ (for differentiable $f$), and under Slater's condition iff it satisfies the KKT conditions.

*Proof of 1.* Suppose $\mathbf{x}^*$ is a local minimiser but some feasible $\mathbf{y}$ has $f(\mathbf{y})<f(\mathbf{x}^*)$. The feasible set is convex, so $\mathbf{z}_\theta=\theta\mathbf{y}+(1-\theta)\mathbf{x}^*$ is feasible, and by convexity $f(\mathbf{z}_\theta)\le\theta f(\mathbf{y})+(1-\theta)f(\mathbf{x}^*)<f(\mathbf{x}^*)$ for every $\theta\in(0,1]$. As $\theta\to0$, $\mathbf{z}_\theta$ comes arbitrarily close to $\mathbf{x}^*$ with a smaller value, contradicting local minimality. $\square$

*Proof of 3.* If $\mathbf{x}_1\ne\mathbf{x}_2$ were both minimisers with value $f^*$, strict convexity would give $f(\tfrac12\mathbf{x}_1+\tfrac12\mathbf{x}_2)<f^*$ at a feasible point. $\square$

Condition 4 says: at the optimum, no feasible direction goes downhill — the convex-set version of "the gradient is zero". In the interior of the feasible set it reduces to $\nabla f(\mathbf{x}^*)=\mathbf{0}$.

### Lagrangian duality

LP duality (§23.3) is a special case of a construction that works for any problem. With the Lagrangian $\mathcal{L}(\mathbf{x},\boldsymbol\mu,\boldsymbol\lambda)=f(\mathbf{x})+\sum\mu_ig_i(\mathbf{x})+\sum\lambda_kh_k(\mathbf{x})$, define the **dual function**

$$
q(\boldsymbol\mu,\boldsymbol\lambda)=\inf_{\mathbf{x}}\ \mathcal{L}(\mathbf{x},\boldsymbol\mu,\boldsymbol\lambda).
$$

For $\boldsymbol\mu\ge\mathbf{0}$ and any feasible $\mathbf{x}$, $\mathcal{L}(\mathbf{x},\boldsymbol\mu,\boldsymbol\lambda)\le f(\mathbf{x})$ (the added terms are $\le0$), so $q(\boldsymbol\mu,\boldsymbol\lambda)\le f^*$: every dual-feasible point gives a lower bound. The **dual problem** is to find the best bound, $\max_{\boldsymbol\mu\ge\mathbf{0},\boldsymbol\lambda}q$. As an infimum of affine functions of $(\boldsymbol\mu,\boldsymbol\lambda)$, $q$ is always concave — the dual is a convex problem even when the primal is not.

> [!theorem] Weak and strong duality
> **Weak duality** $q^*\le f^*$ always holds. **Strong duality** $q^*=f^*$ holds for convex problems that satisfy **Slater's condition**: there is a strictly feasible point, $g_i(\mathbf{x}_0)<0$ for all $i$ (affine inequalities need only be satisfied), with $\mathbf{A}\mathbf{x}_0=\mathbf{b}$. The gap $f(\mathbf{x})-q(\boldsymbol\mu,\boldsymbol\lambda)\ge0$ between any primal-feasible and dual-feasible pair bounds the suboptimality of $\mathbf{x}$ — a **certificate** that solvers use as their stopping test.

> [!example] Worked Example 23.25 — the dual of a small QP
> Primal (Example 23.13): $\min x^2+y^2$ s.t. $2-x-2y\le0$, with $f^*=\tfrac45$ at $(\tfrac25,\tfrac45)$.
>
> *Step 1 (Lagrangian):* $\mathcal{L}=x^2+y^2+\mu(2-x-2y)$.
>
> *Step 2 (minimise over $x,y$):* $\partial\mathcal{L}/\partial x=2x-\mu=0$, $\partial\mathcal{L}/\partial y=2y-2\mu=0$, so $x=\mu/2$, $y=\mu$. Substituting,
> $$
> q(\mu)=\frac{\mu^2}{4}+\mu^2+2\mu-\frac{\mu^2}{2}-2\mu^2=2\mu-\frac54\mu^2 .
> $$
>
> *Step 3 (dual problem):* maximise over $\mu\ge0$: $q'(\mu)=2-\tfrac52\mu=0\Rightarrow\mu^*=\tfrac45$, $q^*=\tfrac85-\tfrac45=\tfrac45$.
>
> *Step 4 (check):* $q^*=f^*=0.8$ — strong duality, as Slater predicts (the point $(2,2)$ is strictly feasible). The dual optimum $\mu^*=0.8$ is the KKT multiplier of Example 23.13. Any $\mu$ gives a bound: $\mu=0.5$ gives $q=0.6875\le0.8$, certifying that no feasible point beats $0.6875$.

### The standard convex problem classes

| Class | Objective | Constraints | Typical EE use |
|---|---|---|---|
| LP | linear | linear | dispatch, flows, minimax / $L^1$ fitting |
| QP | convex quadratic | linear | least squares with bounds, MPC, portfolio-style allocation |
| SOCP | linear | $\lVert\mathbf{A}_i\mathbf{x}+\mathbf{b}_i\rVert\le\mathbf{c}_i^\mathsf{T}\mathbf{x}+d_i$ | robust beamforming, FIR design with magnitude bounds |
| SDP | linear | $\sum x_i\mathbf{F}_i\preceq\mathbf{G}$ | Lyapunov stability, LMI control design, relaxations |
| GP (after $\log$) | posynomial | posynomial $\le1$, monomial $=1$ | component and circuit sizing (§23.10) |

Each class contains the one above it. Interior-point methods solve all of them to high accuracy in a few tens of Newton steps, nearly independent of problem size; modelling tools (CVX, CVXPY) let you write the problem in mathematical form and verify convexity automatically through the composition rules above ("disciplined convex programming").

> [!warning] Common Pitfall — "convex" in the wrong variables
> Convexity depends on the parametrisation. The gain-bandwidth product of an amplifier, a monomial like $c\,W^{0.5}L^{-1}I^{0.5}$ in transistor width, length and current, is neither convex nor concave in $(W,L,I)$, but its logarithm is *linear* in $(\ln W,\ln L,\ln I)$. Before concluding a problem is hopelessly nonconvex, try logarithmic variables, epigraph variables, or a different but equivalent objective (e.g. minimise $\lVert\cdot\rVert^2$ instead of a ratio).

> [!tip] Engineering Note — what to do with a nonconvex problem
> Many real problems (filter design with free phase, AC optimal power flow, neural-network training) are nonconvex. Practical strategies: (i) find a convex reformulation or relaxation (AC power flow → SDP or DC power flow); (ii) run a local method from many starting points and keep the best; (iii) use the convex solution of a simplified model as the starting point for the full model; (iv) use the dual bound to judge how far from optimal a local solution might be.

## 23.8 Least-squares optimization

### Motivation

Least squares is the most-used optimization problem in engineering. Fitting a calibration line to sensor readings, estimating a resistance from noisy voltage and current samples, identifying the parameters of a transfer function from a measured frequency response, designing an FIR filter, solving an overdetermined system of channel equations, training a linear equaliser — all are "minimise a sum of squared errors". It is so prevalent for good reasons: the problem is convex, its solution is given by a linear system, it is the maximum-likelihood estimator under Gaussian noise ([Ch. 20](../part-6-probability-statistics/20-statistical-analysis.md)), and it extends naturally to weighted, regularised and nonlinear versions. Chapter 9 ([§9.11](../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md)) treated least squares geometrically as orthogonal projection; here we treat it as an optimization problem, add regularization, and handle models that are nonlinear in their parameters.

### Linear least squares and the normal equations

Given $\mathbf{A}\in\mathbb{R}^{N\times n}$ ($N\ge n$ measurements, $n$ unknowns) and data $\mathbf{b}\in\mathbb{R}^N$, the system $\mathbf{A}\mathbf{x}=\mathbf{b}$ generally has no solution. We minimise the squared residual norm

$$
f(\mathbf{x})=\lVert\mathbf{A}\mathbf{x}-\mathbf{b}\rVert^2=\mathbf{x}^\mathsf{T}\mathbf{A}^\mathsf{T}\mathbf{A}\mathbf{x}-2\mathbf{b}^\mathsf{T}\mathbf{A}\mathbf{x}+\mathbf{b}^\mathsf{T}\mathbf{b}.
$$

This is a quadratic with gradient $\nabla f=2(\mathbf{A}^\mathsf{T}\mathbf{A}\mathbf{x}-\mathbf{A}^\mathsf{T}\mathbf{b})$ and Hessian $2\mathbf{A}^\mathsf{T}\mathbf{A}\succeq0$, so it is convex and every stationary point is a global minimiser (§23.7).

> [!theorem] Normal equations
> $\mathbf{x}^*$ minimises $\lVert\mathbf{A}\mathbf{x}-\mathbf{b}\rVert^2$ if and only if
> $$
> \mathbf{A}^\mathsf{T}\mathbf{A}\,\mathbf{x}^*=\mathbf{A}^\mathsf{T}\mathbf{b}.
> $$
> The minimiser is unique iff $\mathbf{A}$ has full column rank ($\operatorname{rank}\mathbf{A}=n$), in which case $\mathbf{A}^\mathsf{T}\mathbf{A}\succ0$ and $\mathbf{x}^*=(\mathbf{A}^\mathsf{T}\mathbf{A})^{-1}\mathbf{A}^\mathsf{T}\mathbf{b}$. The optimal residual $\mathbf{r}^*=\mathbf{b}-\mathbf{A}\mathbf{x}^*$ is orthogonal to every column of $\mathbf{A}$: $\mathbf{A}^\mathsf{T}\mathbf{r}^*=\mathbf{0}$.

*Proof.* Setting $\nabla f=\mathbf{0}$ gives the normal equations; convexity makes them sufficient. For uniqueness, $\mathbf{z}^\mathsf{T}\mathbf{A}^\mathsf{T}\mathbf{A}\mathbf{z}=\lVert\mathbf{A}\mathbf{z}\rVert^2$ vanishes only for $\mathbf{z}\in\operatorname{null}\mathbf{A}$, which is $\{\mathbf{0}\}$ exactly when the columns are independent. The orthogonality statement is the normal equations rewritten. $\square$

The orthogonality $\mathbf{A}^\mathsf{T}\mathbf{r}^*=\mathbf{0}$ is the geometry of Chapter 9: $\mathbf{A}\mathbf{x}^*$ is the orthogonal projection of $\mathbf{b}$ onto the column space of $\mathbf{A}$, and the residual is what is left over, perpendicular to that space. It is also a cheap correctness check on any fit.

> [!example] Worked Example 23.26 — fitting a resistor's V–I line
> A resistor with a small offset in the measurement chain gives the readings
>
> | $I$ (mA) | 1 | 2 | 3 | 4 | 5 |
> |---|---|---|---|---|---|
> | $V$ (V) | 2.3 | 4.4 | 6.8 | 8.9 | 11.1 |
>
> Fit $V=RI+V_0$.
>
> *Step 1 (matrix form):* unknowns $\mathbf{x}=(R,V_0)$; row $i$ of $\mathbf{A}$ is $(I_i,1)$, $\mathbf{b}=\mathbf{V}$.
>
> *Step 2 (normal equations):*
> $$
> \mathbf{A}^\mathsf{T}\mathbf{A}=\begin{pmatrix}\sum I_i^2&\sum I_i\\\sum I_i&N\end{pmatrix}=\begin{pmatrix}55&15\\15&5\end{pmatrix},\qquad
> \mathbf{A}^\mathsf{T}\mathbf{b}=\begin{pmatrix}\sum I_iV_i\\\sum V_i\end{pmatrix}=\begin{pmatrix}122.6\\33.5\end{pmatrix}.
> $$
>
> *Step 3 (solve):* determinant $55\cdot5-15^2=50$;
> $$
> R=\frac{5(122.6)-15(33.5)}{50}=2.21,\qquad V_0=\frac{55(33.5)-15(122.6)}{50}=0.07 .
> $$
> Units: $R$ in V/mA $=\text{k}\Omega$, so $R=2.21\ \text{k}\Omega$ and $V_0=70\ \text{mV}$.
>
> *Step 4 (check):* residuals $V_i-(2.21I_i+0.07)=(0.02,\,-0.09,\,0.10,\,-0.01,\,-0.02)$; $\sum r_i=0$ ✓ and $\sum I_ir_i=0.02-0.18+0.30-0.04-0.10=0$ ✓ (orthogonal to both columns). Sum of squares $0.019\ \text{V}^2$, so the rms residual is $\sqrt{0.019/3}=80\ \text{mV}$ (dividing by $N-n=3$ degrees of freedom).

```desmos-graph
left=0; right=6; bottom=0; top=13
---
y=2.21x+0.07|blue
(1,2.3)|red
(2,4.4)|red
(3,6.8)|red
(4,8.9)|red
(5,11.1)|red
```

*Least-squares line $V=2.21I+0.07$ through the five measurements (current in mA, voltage in V). The fit minimises the sum of squared vertical distances.*

### Solving least squares reliably

The normal equations are perfect for derivations and for small, well-conditioned problems, but they have a numerical flaw: $\operatorname{cond}(\mathbf{A}^\mathsf{T}\mathbf{A})=\operatorname{cond}(\mathbf{A})^2$. If $\mathbf{A}$ has condition number $10^6$, forming $\mathbf{A}^\mathsf{T}\mathbf{A}$ in double precision ($\approx16$ digits) leaves about 4 correct digits. Standard practice:

- **QR factorisation.** Write $\mathbf{A}=\mathbf{Q}\mathbf{R}$ with $\mathbf{Q}$ having orthonormal columns and $\mathbf{R}$ upper triangular (Gram–Schmidt or Householder, [Ch. 9 §9.10](../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md)). Since $\lVert\mathbf{A}\mathbf{x}-\mathbf{b}\rVert$ is unchanged by orthogonal transformations, the problem reduces to the triangular system $\mathbf{R}\mathbf{x}=\mathbf{Q}^\mathsf{T}\mathbf{b}$, whose sensitivity is governed by $\operatorname{cond}(\mathbf{A})$, not its square.
- **SVD.** With $\mathbf{A}=\mathbf{U}\boldsymbol\Sigma\mathbf{V}^\mathsf{T}$ ([Ch. 10 §10.7](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md)), $\mathbf{x}^*=\sum_{i}\frac{\mathbf{u}_i^\mathsf{T}\mathbf{b}}{\sigma_i}\mathbf{v}_i$. Small singular values amplify the data component $\mathbf{u}_i^\mathsf{T}\mathbf{b}$ — including its noise — by $1/\sigma_i$. The SVD also gives the minimum-norm solution when $\mathbf{A}$ is rank deficient (the pseudoinverse $\mathbf{A}^+\mathbf{b}$), and it makes regularization transparent (below).

NumPy's `np.linalg.lstsq` uses the SVD; MATLAB's backslash uses QR.

### Weighted least squares

If measurement $i$ has noise variance $\sigma_i^2$, it deserves weight $w_i=1/\sigma_i^2$:

$$
\min_\mathbf{x}\ \sum_iw_i(\mathbf{a}_i^\mathsf{T}\mathbf{x}-b_i)^2=(\mathbf{A}\mathbf{x}-\mathbf{b})^\mathsf{T}\mathbf{W}(\mathbf{A}\mathbf{x}-\mathbf{b}),\qquad \mathbf{A}^\mathsf{T}\mathbf{W}\mathbf{A}\,\mathbf{x}^*=\mathbf{A}^\mathsf{T}\mathbf{W}\mathbf{b}.
$$

More generally, with noise covariance $\boldsymbol\Sigma$, take $\mathbf{W}=\boldsymbol\Sigma^{-1}$. By the Gauss–Markov theorem this gives the minimum-variance linear unbiased estimator, with covariance $(\mathbf{A}^\mathsf{T}\mathbf{W}\mathbf{A})^{-1}$ — the formula that error bars on fitted parameters come from. Weighted LS is ordinary LS on the "whitened" data $\mathbf{W}^{1/2}\mathbf{A}$, $\mathbf{W}^{1/2}\mathbf{b}$.

### Equality-constrained least squares

To minimise $\lVert\mathbf{A}\mathbf{x}-\mathbf{b}\rVert^2$ subject to $\mathbf{C}\mathbf{x}=\mathbf{d}$ (e.g. a filter with exactly unit DC gain, a fit forced through a known point), form the Lagrangian $\lVert\mathbf{A}\mathbf{x}-\mathbf{b}\rVert^2+\boldsymbol\lambda^\mathsf{T}(\mathbf{C}\mathbf{x}-\mathbf{d})$. Its stationarity conditions are the linear **KKT system**

$$
\begin{pmatrix}2\mathbf{A}^\mathsf{T}\mathbf{A}&\mathbf{C}^\mathsf{T}\\\mathbf{C}&\mathbf{0}\end{pmatrix}\begin{pmatrix}\mathbf{x}\\\boldsymbol\lambda\end{pmatrix}=\begin{pmatrix}2\mathbf{A}^\mathsf{T}\mathbf{b}\\\mathbf{d}\end{pmatrix}.
$$

Because the problem is convex with affine constraints, any solution of this system is the global constrained optimum.

### Regularization: ridge (Tikhonov) least squares

When $\mathbf{A}$ is ill-conditioned — nearly collinear regressors, a high-order polynomial fit, deconvolution, an inverse problem in imaging or field reconstruction — the least-squares solution is extremely sensitive to noise. The cure is to add a penalty on the size of $\mathbf{x}$:

> [!definition] Ridge (Tikhonov) regression
> For $\lambda>0$,
> $$
> \mathbf{x}_\lambda=\arg\min_\mathbf{x}\ \lVert\mathbf{A}\mathbf{x}-\mathbf{b}\rVert^2+\lambda\lVert\mathbf{x}\rVert^2
> \quad\Longleftrightarrow\quad
> (\mathbf{A}^\mathsf{T}\mathbf{A}+\lambda\mathbf{I})\,\mathbf{x}_\lambda=\mathbf{A}^\mathsf{T}\mathbf{b}.
> $$
> More generally $\lambda\lVert\mathbf{L}\mathbf{x}\rVert^2$ with a matrix $\mathbf{L}$ (e.g. a difference operator, to penalise roughness rather than size).

*Derivation.* The gradient of the objective is $2(\mathbf{A}^\mathsf{T}\mathbf{A}+\lambda\mathbf{I})\mathbf{x}-2\mathbf{A}^\mathsf{T}\mathbf{b}$. Since $\mathbf{A}^\mathsf{T}\mathbf{A}+\lambda\mathbf{I}$ has eigenvalues $\sigma_i^2+\lambda\ge\lambda>0$, it is always invertible — even if $\mathbf{A}$ is rank deficient — and the objective is $2\lambda$-strongly convex, so $\mathbf{x}_\lambda$ exists and is unique. Equivalently, ridge is ordinary LS on the augmented system $\begin{pmatrix}\mathbf{A}\\\sqrt\lambda\,\mathbf{I}\end{pmatrix}\mathbf{x}\approx\begin{pmatrix}\mathbf{b}\\\mathbf{0}\end{pmatrix}$, which is how it is computed stably.

**SVD view.** Substituting $\mathbf{A}=\mathbf{U}\boldsymbol\Sigma\mathbf{V}^\mathsf{T}$,

$$
\mathbf{x}_\lambda=\sum_i\underbrace{\frac{\sigma_i^2}{\sigma_i^2+\lambda}}_{\text{filter factor }\phi_i}\ \frac{\mathbf{u}_i^\mathsf{T}\mathbf{b}}{\sigma_i}\,\mathbf{v}_i .
$$

Components with $\sigma_i^2\gg\lambda$ pass almost unchanged ($\phi_i\approx1$); components with $\sigma_i^2\ll\lambda$ — the ones that amplify noise — are suppressed ($\phi_i\approx\sigma_i^2/\lambda$). Ridge is a smooth low-pass filter in "singular-value space"; truncated SVD is the brick-wall version.

**Trade-off.** As $\lambda$ grows, the residual $\lVert\mathbf{A}\mathbf{x}_\lambda-\mathbf{b}\rVert$ increases and the norm $\lVert\mathbf{x}_\lambda\rVert$ decreases monotonically: a bias–variance trade-off, and in fact a two-objective problem whose weighted-sum solutions trace a Pareto front (§23.9). The parameter is chosen by cross-validation, by the "corner" of the log–log plot of norm versus residual (the **L-curve**), or by the discrepancy principle (choose $\lambda$ so that the residual matches the known noise level). Statistically, ridge is the MAP estimate when $\mathbf{x}$ has a Gaussian prior with variance $\sigma^2_\text{noise}/\lambda$.

> [!example] Worked Example 23.27 — ridge stabilises a nearly collinear fit
> Fit $\mathbf{x}=(x_1,x_2)$ to
> $$
> \mathbf{A}=\begin{pmatrix}1&1\\1&1.01\\1&0.99\end{pmatrix},\qquad \mathbf{b}=\begin{pmatrix}2\\2.03\\1.96\end{pmatrix}.
> $$
> The data were generated from $\mathbf{x}_\text{true}=(1,1)$ with perturbations of order $0.02$.
>
> *Step 1 (diagnose):* singular values $\sigma_1=2.449$, $\sigma_2=0.0100$; $\operatorname{cond}(\mathbf{A})=245$ and $\operatorname{cond}(\mathbf{A}^\mathsf{T}\mathbf{A})\approx6\times10^4$. The second column is almost the first.
>
> *Step 2 (plain LS):* $\mathbf{x}^*=(-1.503,\ 3.500)$ with residual norm $0.0041$. The fit is excellent and the parameters are nonsense: the small noise has been amplified by $1/\sigma_2=100$ along $\mathbf{v}_2\approx(1,-1)/\sqrt2$.
>
> *Step 3 (ridge):*
>
> | $\lambda$ | $\mathbf{x}_\lambda$ | $\lVert\mathbf{A}\mathbf{x}_\lambda-\mathbf{b}\rVert$ | $\phi_2=\sigma_2^2/(\sigma_2^2+\lambda)$ |
> |---|---|---|---|
> | 0 | $(-1.503,\ 3.500)$ | 0.0041 | 1 |
> | $10^{-3}$ | $(0.771,\ 1.226)$ | 0.0324 | 0.091 |
> | $10^{-2}$ | $(0.972,\ 1.021)$ | 0.0357 | 0.0099 |
> | $10^{-1}$ | $(0.979,\ 0.985)$ | 0.0669 | 0.001 |
> | 1 | $(0.855,\ 0.856)$ | 0.495 | $10^{-4}$ |
>
> *Step 4 (interpret):* a tiny $\lambda=0.01$ (compared with $\sigma_1^2=6$) removes the unstable component and gives parameters within 3 % of the truth at the cost of a residual comparable to the noise level ($\lVert\mathbf{b}-\mathbf{A}\mathbf{x}_\text{true}\rVert=0.028$). Too large a $\lambda$ ($=1$) starts shrinking the well-determined component as well — bias.

### Nonlinear least squares

Many models are nonlinear in their parameters: an exponential decay $v(t)=ae^{-bt}$, a resonance $\lvert H(j\omega)\rvert=K/\sqrt{(1-\omega^2/\omega_0^2)^2+(\omega/(Q\omega_0))^2}$, a diode law $I=I_s(e^{V/(nV_T)}-1)$. With residuals $r_i(\mathbf{p})=y_i-m(t_i;\mathbf{p})$ we minimise

$$
F(\mathbf{p})=\tfrac12\sum_{i=1}^Nr_i(\mathbf{p})^2=\tfrac12\lVert\mathbf{r}(\mathbf{p})\rVert^2 .
$$

With the residual Jacobian $\mathbf{J}=\partial\mathbf{r}/\partial\mathbf{p}$ ($N\times n$), the gradient and Hessian are

$$
\nabla F=\mathbf{J}^\mathsf{T}\mathbf{r},\qquad \nabla^2F=\mathbf{J}^\mathsf{T}\mathbf{J}+\sum_ir_i\nabla^2r_i .
$$

The second Hessian term is small when the residuals are small (a good fit) or the model is nearly linear. Dropping it gives an approximation that needs only first derivatives and is automatically positive semidefinite.

> [!definition] Gauss–Newton and Levenberg–Marquardt
> **Gauss–Newton:** solve the linear least-squares problem for the step,
> $$
> (\mathbf{J}_k^\mathsf{T}\mathbf{J}_k)\,\boldsymbol\delta_k=-\mathbf{J}_k^\mathsf{T}\mathbf{r}_k,\qquad \mathbf{p}_{k+1}=\mathbf{p}_k+\boldsymbol\delta_k,
> $$
> i.e. linearise the model about $\mathbf{p}_k$ and fit the linearised model exactly by least squares.
> **Levenberg–Marquardt:** $(\mathbf{J}_k^\mathsf{T}\mathbf{J}_k+\tau_k\mathbf{I})\boldsymbol\delta_k=-\mathbf{J}_k^\mathsf{T}\mathbf{r}_k$, with $\tau_k$ decreased after a successful step and increased after a failed one — a trust-region Gauss–Newton method (§23.6) that behaves like gradient descent far from the solution and like Gauss–Newton near it.

(If $\mathbf{r}=\mathbf{y}-\mathbf{m}$, then $\mathbf{J}=-\partial\mathbf{m}/\partial\mathbf{p}$; with the model Jacobian $\mathbf{J}_m=\partial\mathbf{m}/\partial\mathbf{p}$ the Gauss–Newton step reads $\mathbf{J}_m^\mathsf{T}\mathbf{J}_m\boldsymbol\delta=\mathbf{J}_m^\mathsf{T}\mathbf{r}$.) Gauss–Newton converges quadratically for zero-residual problems and linearly (often fast) for small-residual ones; it can fail on large-residual or badly initialised problems, which is when LM's damping matters. A good starting point is essential, and is often available from a linearising transformation.

> [!example] Worked Example 23.28 — fitting an RC discharge by Gauss–Newton
> A capacitor discharging through a resistor is sampled:
>
> | $t$ (ms) | 0 | 1 | 2 | 3 | 4 |
> |---|---|---|---|---|---|
> | $v$ (V) | 5.0 | 3.0 | 1.9 | 1.1 | 0.7 |
>
> Fit $v=ae^{-bt}$ (with $b=1/\tau$).
>
> *Step 1 (starting guess by linearisation):* $\ln v=\ln a-bt$ is linear. A straight-line LS fit of $\ln v_i$ against $t_i$ gives $b=0.4936\ \text{ms}^{-1}$, $\ln a=1.6048$, $a=4.977$ V. (This minimises errors in $\ln v$, which over-weights the small late samples, so it is a starting point rather than the answer.)
>
> *Step 2 (Jacobian):* the model Jacobian has rows $\big(\partial m/\partial a,\ \partial m/\partial b\big)=\big(e^{-bt_i},\ -at_ie^{-bt_i}\big)$.
>
> *Step 3 (iterate from a deliberately poor start $(a,b)=(4,\ 0.2)$):*
>
> | $k$ | $a$ (V) | $b$ (ms$^{-1}$) | $\sum r_i^2$ (V$^2$) |
> |---|---|---|---|
> | 0 | 4.0000 | 0.2000 | 4.09 |
> | 1 | 4.8385 | 0.4319 | 0.142 |
> | 2 | 4.9844 | 0.4920 | 0.00540 |
> | 3 | 4.98933 | 0.495605 | 0.00496029 |
> | 4 | 4.989372 | 0.4956304 | 0.00496027 |
>
> *Step 4 (result):* $a=4.989\ \text{V}$, $b=0.4956\ \text{ms}^{-1}$, so $\tau=1/b=2.018\ \text{ms}$; residual rms $\sqrt{0.00496/3}=41\ \text{mV}$. Starting from the log-linear guess, Gauss–Newton converges in two iterations.
>
> *Step 5 (check):* at the solution $\mathbf{J}_m^\mathsf{T}\mathbf{r}=\mathbf{0}$ to machine precision — the nonlinear normal equations hold ✓.

```desmos-graph
left=-0.3; right=5; bottom=0; top=5.5
---
y=4.989e^{-0.4956x}|blue
y=4.977e^{-0.4936x}|green|dashed
(0,5)|red
(1,3)|red
(2,1.9)|red
(3,1.1)|red
(4,0.7)|red
```

*RC discharge data (red), the nonlinear least-squares fit $4.989e^{-0.4956t}$ (blue) and the log-linear starting fit (green dashed); time in ms, voltage in V. The two curves nearly coincide here because the noise is small.*

> [!warning] Common Pitfall — fitting in transformed coordinates
> Fitting $\ln y$ instead of $y$, or $1/y$ instead of $y$, changes the error weighting: a $0.05$ V error on a $0.7$ V sample is a 7 % error in $\ln v$, but on a $5$ V sample only 1 %. Linearising transforms give excellent starting values; the final fit should minimise errors in the quantity that was actually measured with known noise (or use weights $w_i\propto y_i^2$ to compensate).

> [!tip] Engineering Note — system identification
> Fitting a transfer function $H(s)=\frac{b_0+b_1s+\dots}{1+a_1s+\dots}$ to measured frequency-response data is nonlinear in the $a_k$. Levy's method multiplies through by the denominator to obtain a *linear* least-squares problem (an equation-error formulation); Sanathanan–Koerner iteration reweights it to approximate the true output-error criterion, and the result seeds a Levenberg–Marquardt refinement. Vector fitting, used for passive macromodels of interconnects, refines this idea with pole relocation.

```python
import numpy as np
from scipy.optimize import least_squares

t = np.array([0, 1, 2, 3, 4.0]);  v = np.array([5.0, 3.0, 1.9, 1.1, 0.7])
res = least_squares(lambda p: v - p[0]*np.exp(-p[1]*t), x0=[4.0, 0.2], method="lm")
print(res.x)                                   # approx [4.989, 0.4956]

A = np.array([[1, 1], [1, 1.01], [1, 0.99]]);  b = np.array([2, 2.03, 1.96]);  lam = 1e-2
x_ridge = np.linalg.lstsq(np.vstack([A, np.sqrt(lam)*np.eye(2)]),
                          np.concatenate([b, np.zeros(2)]), rcond=None)[0]
print(x_ridge)                                 # approx [0.972, 1.021]
```

## 23.9 Multi-objective optimization

### Motivation

Real design specifications rarely reduce to one number. A low-noise amplifier should have low noise *and* low power consumption; a filter should have a sharp transition *and* low passband ripple *and* few taps; a power system should run at low cost *and* low emissions *and* high reliability; a regularised fit should have small residual *and* small parameter norm (§23.8). These goals conflict: improving one makes another worse. There is then no single "optimal" design but a family of best compromises, and the engineer's job is to *see* that family and choose from it. Multi-objective optimization provides the vocabulary (dominance, Pareto optimality) and the standard techniques (weighted sums, $\varepsilon$-constraints) that turn a vector of objectives back into ordinary single-objective problems.

### Dominance and Pareto optimality

Consider $\min_{\mathbf{x}\in S}\ \mathbf{F}(\mathbf{x})=\big(f_1(\mathbf{x}),\dots,f_q(\mathbf{x})\big)$.

> [!definition] Pareto optimality
> - A feasible $\mathbf{x}$ **dominates** $\mathbf{y}$ if $f_i(\mathbf{x})\le f_i(\mathbf{y})$ for all $i$ and $f_j(\mathbf{x})<f_j(\mathbf{y})$ for at least one $j$: $\mathbf{x}$ is at least as good in every objective and strictly better in one.
> - $\mathbf{x}^*\in S$ is **Pareto optimal** (efficient, non-dominated) if no feasible point dominates it. It is **weakly Pareto optimal** if no feasible point is strictly better in *every* objective.
> - The **Pareto set** is the set of Pareto-optimal decisions; its image $\{\mathbf{F}(\mathbf{x}^*)\}$ in objective space is the **Pareto front** (trade-off curve).
> - The **ideal (utopia) point** $\mathbf{F}^\text{id}$ has components $f_i^\text{id}=\min_{S}f_i$, each objective minimised on its own. It is usually not attainable.

At a Pareto-optimal point you cannot improve any objective without worsening another. A design that is *not* Pareto optimal is simply wasteful: some other design is better in at least one respect and worse in none. The first, uncontroversial, step of any trade-off study is therefore to discard dominated designs; only the choice *along* the front involves preferences.

The slope of a smooth two-objective front, $df_2/df_1<0$, is the **marginal rate of trade-off**: how many units of $f_2$ must be given up per unit of $f_1$ gained. Near the ends of the front this rate is often extreme — the last few percent of improvement in one objective is bought at great cost in the other — and a "knee" in the middle marks the region where a balanced compromise lives.

### The weighted-sum method

Combine the objectives with nonnegative weights:

$$
\min_{\mathbf{x}\in S}\ \sum_{i=1}^q w_if_i(\mathbf{x}),\qquad w_i\ge0,\ \ \sum w_i=1 .
$$

> [!theorem] Weighted sums and the Pareto front
> 1. If all $w_i>0$, every minimiser of the weighted sum is Pareto optimal. (With some $w_i=0$, it is at least weakly Pareto optimal.)
> 2. If the problem is **convex** (convex $S$ and convex $f_i$), then for every Pareto-optimal $\mathbf{x}^*$ there are weights $w_i\ge0$, not all zero, for which $\mathbf{x}^*$ minimises the weighted sum.
> 3. For nonconvex problems, Pareto-optimal points on non-convex parts of the front cannot be found by any choice of weights.

*Proof of 1.* If some $\mathbf{y}$ dominated a weighted-sum minimiser $\mathbf{x}^*$, then $\sum w_if_i(\mathbf{y})<\sum w_if_i(\mathbf{x}^*)$ (strict because one $f_j$ is strictly smaller and $w_j>0$), contradicting minimality. $\square$

*Why 2 and 3 hold (two objectives).* In the $(f_1,f_2)$ plane the weighted sum $w_1f_1+w_2f_2=c$ is a straight line with slope $-w_1/w_2$. Minimising $c$ slides this line toward the origin until it last touches the set of attainable objective vectors. The touching point is a point where the line **supports** the attainable set. For a convex problem the attainable set (plus everything dominated by it) is convex, and every boundary point of a convex set has a supporting line — so every front point is reachable. A part of the front that bulges *away* from the origin cannot be touched by a supporting line from below-left; the sliding line always jumps past it to the ends.

At a weighted-sum optimum on a smooth front, the line is tangent to the front, so the **front slope equals $-w_1/w_2$**: the weights *are* the trade-off rate one is willing to accept.

> [!warning] Common Pitfall — weights depend on units
> $w_1f_1+w_2f_2$ with $f_1$ in milliwatts and $f_2$ in $\text{nV}^2/\text{Hz}$ means nothing until the units are fixed; changing $f_1$ to watts changes the effective weight by 1000. Normalise each objective first, for example $\tilde f_i=(f_i-f_i^\text{id})/(f_i^\text{nad}-f_i^\text{id})$ using the ideal point and the **nadir point** (worst value of $f_i$ over the Pareto set), so that each normalised objective runs from 0 to 1 along the front. Even then, equal weights do not guarantee "equal importance"; sweep the weights and look at the front.

### The $\varepsilon$-constraint method

Keep one objective and turn the others into constraints:

$$
\min_{\mathbf{x}\in S}\ f_1(\mathbf{x})\quad\text{subject to}\quad f_i(\mathbf{x})\le\varepsilon_i,\ \ i=2,\dots,q .
$$

> [!theorem] $\varepsilon$-constraint method
> 1. Every solution of the $\varepsilon$-constraint problem is weakly Pareto optimal; if it is unique, it is Pareto optimal.
> 2. Conversely, every Pareto-optimal $\mathbf{x}^*$ solves the $\varepsilon$-constraint problem with $\varepsilon_i=f_i(\mathbf{x}^*)$ — **convex or not**.
> 3. The KKT multiplier $\mu_i$ of the constraint $f_i\le\varepsilon_i$ gives the local trade-off rate: $\partial f_1^*/\partial\varepsilon_i=-\mu_i$.

*Proof of 2.* If some feasible $\mathbf{y}$ had $f_i(\mathbf{y})\le f_i(\mathbf{x}^*)$ for $i\ge2$ and $f_1(\mathbf{y})<f_1(\mathbf{x}^*)$, it would dominate $\mathbf{x}^*$. $\square$

The $\varepsilon$-constraint form is also how specifications are usually written ("minimise power subject to noise $\le$ spec"), so it connects directly to the design process: sweeping $\varepsilon$ traces the front, and the multiplier tells the designer what each tightening of the specification costs. Its drawbacks are that it needs sensible ranges for $\varepsilon_i$ and that each point requires a constrained solve.

Other scalarisations worth knowing: **goal programming / compromise programming** minimises the distance $\lVert\mathbf{F}(\mathbf{x})-\mathbf{F}^\text{id}\rVert$ (weighted, in some norm) to the ideal point; the **Chebyshev (min–max) scalarisation** $\min\max_iw_i(f_i-f_i^\text{id})$ can reach every Pareto point even on nonconvex fronts. Population-based evolutionary algorithms (NSGA-II and relatives) approximate a whole front at once for difficult nonconvex problems.

> [!example] Worked Example 23.29 — a convex front traced two ways
> Objectives $f_1(x)=x^2$ and $f_2(x)=(x-2)^2$ over $x\in\mathbb{R}$ (two targets 0 and 2 that cannot both be met).
>
> *Step 1 (Pareto set):* for $x<0$, moving toward 0 improves both; for $x>2$, moving toward 2 improves both. For $0\le x\le2$, increasing $x$ helps $f_2$ and hurts $f_1$. Pareto set: $[0,2]$. Ideal point $(0,0)$, not attainable.
>
> *Step 2 (front):* on $[0,2]$, $x=\sqrt{f_1}$, so the front is $f_2=(\sqrt{f_1}-2)^2$, $0\le f_1\le4$ — a convex curve from $(0,4)$ to $(4,0)$.
>
> *Step 3 (weighted sum):* minimise $wx^2+(1-w)(x-2)^2$: derivative $2wx+2(1-w)(x-2)=0$ gives $x^*=2(1-w)$. As $w$ runs from 1 to 0, $x^*$ sweeps the whole Pareto set — every front point is reached, as the convex theory promises. For $w=\tfrac12$: $x^*=1$, $(f_1,f_2)=(1,1)$.
>
> *Step 4 ($\varepsilon$-constraint):* minimise $x^2$ s.t. $(x-2)^2\le\varepsilon$. For $0<\varepsilon<4$ the constraint is active: $x^*=2-\sqrt\varepsilon$, $f_1^*=(2-\sqrt\varepsilon)^2$. With $\varepsilon=1$: $x^*=1$, $f_1^*=1$ — the same point as $w=\tfrac12$.
>
> *Step 5 (trade-off rate):* $df_1^*/d\varepsilon=-(2-\sqrt\varepsilon)/\sqrt\varepsilon=-1$ at $\varepsilon=1$, so $\mu=1$. From KKT directly: $2x+\mu\cdot2(x-2)=0$ at $x=1$ gives $\mu=1$ ✓. And the weighted-sum tangent slope $-w/(1-w)=-1$ ✓: the multiplier and the weight ratio carry the same information.

```desmos-graph
left=-0.3; right=5; bottom=-0.3; top=5
---
y=(\sqrt{x}-2)^2|blue|0<x<4
y=2-x|red|dashed
(1,1)|label:w=1/2 and eps=1
(0,0)|label:ideal point
(0,4)
(4,0)
```

*Pareto front of $f_1=x^2$, $f_2=(x-2)^2$ (blue). The weighted-sum line for $w=\tfrac12$ (red dashed, slope $-1$) is tangent to the front at $(1,1)$; sweeping $w$ rolls the tangent line along the whole convex front.*

> [!example] Worked Example 23.30 — a nonconvex front defeats the weighted sum
> Objectives $f_1(x)=x$ and $f_2(x)=1-x^2$ on $S=[0,1]$.
>
> *Step 1 (front):* $f_1$ increases and $f_2$ decreases along $S$, so every $x\in[0,1]$ is Pareto optimal and the front is $f_2=1-f_1^2$, bulging *away* from the ideal point $(0,0)$.
>
> *Step 2 (weighted sum):* $\phi(x)=wx+(1-w)(1-x^2)$ has $\phi''=-2(1-w)\le0$: it is concave in $x$, so its minimum over $[0,1]$ is at an endpoint. Comparing $\phi(0)=1-w$ and $\phi(1)=w$: the minimiser is $x=1$ for $w<\tfrac12$ and $x=0$ for $w>\tfrac12$. **Only the two end points of the front are ever found**, no matter how finely $w$ is swept.
>
> *Step 3 ($\varepsilon$-constraint):* minimise $x$ s.t. $1-x^2\le\varepsilon$: $x^*=\sqrt{1-\varepsilon}$. For $\varepsilon=0.75$, $x^*=0.5$ and $(f_1,f_2)=(0.5,0.75)$, an interior front point unreachable by weights. Every front point is obtained as $\varepsilon$ sweeps $[0,1]$ ✓.

```desmos-graph
left=-0.2; right=1.4; bottom=-0.2; top=1.4
---
y=1-x^2|blue|0<x<1
y=1-x|red|dashed
y=0.75|green|dotted
(0.5,0.75)|label:eps-constraint point
(0,1)|label:weighted sum
(1,0)|label:weighted sum
```

*A nonconvex (concave-shaped) front $f_2=1-f_1^2$. The equal-weight line (red dashed) touches only the end points; the $\varepsilon$-constraint $f_2\le0.75$ (green dotted) picks out the interior point $(0.5,0.75)$.*

> [!example] Worked Example 23.31 — power versus noise in an amplifier stage
> For a bipolar input stage biased at collector current $I$ (mA) from a 3.3 V supply, model the power as $f_1=3.3I$ (mW) and the input voltage-noise density (shot-noise dominated) as $f_2=4/I$ (in $\text{nV}^2/\text{Hz}$), with $0.5\le I\le10$ mA.
>
> *Step 1 (front):* $f_1$ rises and $f_2$ falls with $I$, so every feasible $I$ is Pareto optimal. Eliminating $I$: $f_1f_2=13.2$, a hyperbola from $(1.65,8)$ to $(33,0.4)$ — a convex front.
>
> *Step 2 ($\varepsilon$-constraint, the specification form):* minimise $3.3I$ subject to $4/I\le2$ (noise spec $\sqrt2=1.41\ \text{nV}/\sqrt{\text{Hz}}$). The constraint gives $I\ge2$, so $I^*=2$ mA and $f_1^*=6.6$ mW.
>
> *Step 3 (multiplier):* stationarity $3.3+\mu(-4/I^2)=0$ at $I=2$ gives $\mu=3.3$ mW per $\text{nV}^2/\text{Hz}$. Tightening the spec to $f_2\le1.9$ should cost about $0.1\mu=0.33$ mW; the exact re-solve gives $I=2.105$ mA, $f_1=6.95$ mW, an increase of $0.35$ mW ✓.
>
> *Step 4 (weighted sum, normalised):* with weights on the raw objectives, $w(3.3I)+(1-w)(4/I)$ is minimised at $I^*=\sqrt{4(1-w)/(3.3w)}$; $w=\tfrac12$ gives $I^*=1.10$ mA, $(f_1,f_2)=(3.63,3.63)$ — a point that depends entirely on the arbitrary choice of units. Expressing the same requirement as a noise specification (Step 2) is usually more meaningful to the designer than choosing a weight.

> [!tip] Engineering Note — trade-off curves in design reviews
> Presenting a Pareto front instead of a single "optimal" design changes design conversations: it makes explicit what each specification costs, reveals knees where a small concession buys a large gain, and exposes dominated legacy designs. Typical fronts in EE: noise figure versus power in LNAs, settling time versus overshoot in control loops, filter order versus transition width, conversion efficiency versus power density in converters, and residual versus solution norm (the L-curve) in regularised estimation.

## 23.10 Applications to engineering design

### Motivation

The previous sections developed methods; this one shows how they fit together on problems that look like real design work. The pattern is always the same: (1) choose decision variables and put them in sensible units or logarithms; (2) write the objective and the specifications as functions; (3) recognise the problem class — LP, QP, convex, least squares, nonconvex; (4) solve with the matching method; (5) read the multipliers as sensitivities and check the design against the original, un-idealised requirements (discrete component values, tolerances, unmodelled effects). The three applications below — generator dispatch, component sizing, and FIR filter design — each exercise a different part of the chapter. They are motivation for the mathematics; the engineering details belong in power-systems, circuits and DSP texts.

### Application 1: economic dispatch of generators

Generators $i=1,\dots,G$ must together supply a demand $D$ (MW). Each has a convex fuel-cost curve, commonly quadratic, $C_i(P_i)=a_iP_i^2+b_iP_i+c_i$ ($\$/\text{h}$, with $a_i>0$), and operating limits $P_i^{\min}\le P_i\le P_i^{\max}$. Losses are neglected here.

$$
\min_{\mathbf{P}}\ \sum_{i=1}^G C_i(P_i)\quad\text{s.t.}\quad\sum_{i=1}^G P_i=D,\qquad P_i^{\min}\le P_i\le P_i^{\max}.
$$

This is a convex QP: a convex separable objective, one affine equality and box constraints, so the KKT conditions are necessary and sufficient (§23.4, §23.7).

**Without limits.** The Lagrangian $\mathcal{L}=\sum C_i(P_i)+\lambda\big(D-\sum P_i\big)$ gives $\partial\mathcal{L}/\partial P_i=C_i'(P_i)-\lambda=0$:

> [!theorem] Equal incremental cost criterion
> At the optimal unconstrained dispatch, every generator operates at the same **incremental cost** $dC_i/dP_i=\lambda$. For quadratic costs $P_i=(\lambda-b_i)/(2a_i)$, and the power balance gives
> $$
> \lambda=\frac{D+\sum_ib_i/(2a_i)}{\sum_i1/(2a_i)} .
> $$
> The multiplier $\lambda$ ($\$/\text{MWh}$) is the **system marginal cost**: $dC^*/dD=\lambda$.

The reason is the forgiving-optimum argument in reverse: if one unit had a lower incremental cost than another, shifting a little output from the expensive unit to the cheap one would lower total cost while keeping the balance.

**With limits.** With multipliers $\underline\mu_i,\overline\mu_i\ge0$ for the lower and upper bounds, stationarity reads $C_i'(P_i)-\lambda-\underline\mu_i+\overline\mu_i=0$. Complementary slackness then gives the textbook rule:

$$
C_i'(P_i)=\lambda\ \ \text{if}\ P_i^{\min}<P_i<P_i^{\max},\qquad C_i'(P_i)\le\lambda\ \ \text{if}\ P_i=P_i^{\max},\qquad C_i'(P_i)\ge\lambda\ \ \text{if}\ P_i=P_i^{\min}.
$$

A unit at its maximum would *like* to produce more (its incremental cost is below the system price); the gap $\overline\mu_i=\lambda-C_i'(P_i^{\max})$ is the value of one more MW of its capacity.

> [!example] Worked Example 23.32 — dispatching three units
> Units with $a=(0.005,\ 0.008,\ 0.010)\ \$/\text{MW}^2\text{h}$, $b=(6.0,\ 6.5,\ 7.0)\ \$/\text{MWh}$ (fixed costs $c_i$ do not affect the dispatch) supply $D=700$ MW.
>
> *Step 1 (no limits):* $\sum1/(2a_i)=100+62.5+50=212.5$ MW per $\$/\text{MWh}$, and $\sum b_i/(2a_i)=600+406.25+350=1356.25$ MW. So
> $$
> \lambda=\frac{700+1356.25}{212.5}=9.676\ \$/\text{MWh},\qquad P=\Big(\frac{9.676-6}{0.01},\ \frac{9.676-6.5}{0.016},\ \frac{9.676-7}{0.02}\Big)=(367.6,\ 198.5,\ 133.8)\ \text{MW}.
> $$
> Check: sum $=700.0$ ✓; incremental costs $2a_iP_i+b_i$ all equal $9.676$ ✓. Variable cost $\sum(a_iP_i^2+b_iP_i)=5603.3\ \$/\text{h}$.
>
> *Step 2 (add the limit $P_1\le300$ MW):* unit 1 violates it, so guess it is active: $P_1=300$ and units 2 and 3 share the remaining 400 MW at a common $\lambda$:
> $$
> \lambda=\frac{400+406.25+350}{62.5+50}=10.278\ \$/\text{MWh},\qquad P_2=236.1,\ \ P_3=163.9\ \text{MW}.
> $$
>
> *Step 3 (verify KKT):* unit 1's incremental cost at 300 MW is $2(0.005)(300)+6=9.0\le\lambda$ ✓, so $\overline\mu_1=10.278-9.0=1.278\ \$/\text{MWh}\ge0$ ✓. The limit raises the cost to $5646.5\ \$/\text{h}$ ($+43.2$) and raises the system marginal cost from 9.68 to 10.28 $\$/\text{MWh}$.
>
> *Step 4 (interpret):* $\overline\mu_1$ says one more MW of capacity on unit 1 would save about $\$1.28$ per hour — the information a planner needs to decide whether an uprate is worth it.

> [!example] Worked Example 23.33 — linear costs give merit-order dispatch
> If the same three units have *linear* costs $6.0,\ 6.5,\ 7.0\ \$/\text{MWh}$, each limited to 300 MW, and $D=700$ MW, the problem is an LP. By the vertex theorem (§23.3) the optimum loads units fully in order of cost: $P=(300,\ 300,\ 100)$ MW, cost $1800+1950+700=4450\ \$/\text{h}$. The dual variable of the balance constraint is $\lambda=7.0\ \$/\text{MWh}$, the cost of the **marginal unit** (the only one not at a limit), and the capacity constraints of units 1 and 2 have shadow prices $7.0-6.0=1.0$ and $7.0-6.5=0.5\ \$/\text{MWh}$, by complementary slackness.

### Application 2: sizing components — log variables and geometric programming

Component values are positive and span decades, and many circuit relationships are **monomials** $c\,x_1^{a_1}x_2^{a_2}\cdots$ (time constants $RC$, cutoff $1/(2\pi RC)$, $kT/C$ noise, gain $g_mR$). Sums of monomials with positive coefficients are **posynomials**. A **geometric program** (GP) minimises a posynomial subject to posynomial $\le1$ and monomial $=1$ constraints. It is not convex as written, but in logarithmic variables $u_j=\ln x_j$:

- a monomial becomes $e^{\mathbf{a}^\mathsf{T}\mathbf{u}+\ln c}$, and the constraint "monomial $=1$" becomes *linear*, $\mathbf{a}^\mathsf{T}\mathbf{u}+\ln c=0$;
- a posynomial becomes a sum of exponentials of affine functions, which is convex; taking its logarithm gives the convex log-sum-exp form.

So every GP is a convex problem in disguise, with all the guarantees of §23.7. GP is used for sizing transistors in op-amps, wires in clock trees and gates in digital paths.

> [!example] Worked Example 23.34 — sizing an RC low-pass filter
> Choose $R$ and $C$ for a first-order low-pass with time constant $\tau=RC=1$ ms (cutoff 159 Hz). Board area is modelled as $J=R/R_0+C/C_0$ with $R_0=1\ \text{k}\Omega$ and $C_0=10$ nF (one unit of area per kilohm and per 10 nF). The integrated output noise $\sqrt{kT/C}$ at $T=300$ K must not exceed $0.1\ \mu\text{V}$.
>
> *Step 1 (constraints in physical form):* $kT=4.14\times10^{-21}$ J, so $\sqrt{kT/C}\le10^{-7}\ \text{V}\iff C\ge kT/10^{-14}=414$ nF.
>
> *Step 2 (log variables):* let $u=\ln(R/R_0)$, $v=\ln(C/C_0)$. Then $RC=\tau$ becomes $u+v=\ln\big(\tau/(R_0C_0)\big)=\ln100=4.605$, the noise constraint becomes $v\ge\ln41.4=3.724$, and the objective is $J=e^u+e^v$ — convex. The problem is a convex program with one linear equality and one linear inequality.
>
> *Step 3 (ignore the noise constraint first):* by symmetry of $e^u+e^v$ on the line $u+v=4.605$, the optimum is $u=v=2.303$: $R=10\ \text{k}\Omega$, $C=100$ nF, $J=20$. But $C=100$ nF violates $C\ge414$ nF, so the noise constraint must be active.
>
> *Step 4 (active noise constraint):* $v=3.724$, $u=4.605-3.724=0.881$: $C=414$ nF, $R=\tau/C=2.41\ \text{k}\Omega$, $J=2.41+41.4=43.8$.
>
> *Step 5 (KKT check and sensitivity):* with $\mathcal{L}=e^u+e^v+\lambda(u+v-4.605)+\mu(3.724-v)$: $\partial_u\mathcal{L}=e^u+\lambda=0\Rightarrow\lambda=-2.41$; $\partial_v\mathcal{L}=e^v+\lambda-\mu=0\Rightarrow\mu=41.4-2.41=39.0\ge0$ ✓. In log variables $\mu$ is the cost per unit of $\ln C_{\min}$: a 1 % looser noise budget on $C_{\min}$ saves about $0.39$ area units (re-solving with $C_{\min}=410$ nF gives $J=43.44$, a saving of $0.40$ ✓). The noise specification is what drives the design.
>
> *Step 6 (round to preferred values):* the nearest E12 capacitor not below 414 nF is 470 nF; then $R=\tau/C=2.13\ \text{k}\Omega$, and the nearest E96 resistor 2.15 k$\Omega$ gives $\tau=1.011$ ms ($+1.1\,\%$). The area becomes $2.15+47=49.2$ — 12 % above the continuous optimum, mostly because the E12 capacitor series is coarse. Rounding is a discrete optimisation; the continuous solution tells you where to look.

### Application 3: FIR filter design by least squares and by linear programming

A linear-phase FIR low-pass filter of length $2M+1$ has amplitude response $A(\omega)=h_0+2\sum_{k=1}^Mh_k\cos k\omega$ (Example 23.3). Design means choosing $\mathbf{h}$ to approximate the ideal response $D(\omega)$.

**Continuous least squares.** Minimise $E(\mathbf{h})=\int_0^\pi\big(A(\omega)-D(\omega)\big)^2d\omega$. Because the functions $1,\cos\omega,\dots,\cos M\omega$ are orthogonal on $[0,\pi]$ — $\int_0^\pi\cos k\omega\cos l\omega\,d\omega=\tfrac\pi2\delta_{kl}$ for $k,l\ge1$ and $\int_0^\pi1\,d\omega=\pi$ — the normal equations are *diagonal*, and each coefficient is found independently:

$$
h_0=\frac1\pi\int_0^\pi D(\omega)\,d\omega,\qquad h_k=\frac1\pi\int_0^\pi D(\omega)\cos k\omega\,d\omega .
$$

These are the Fourier-series coefficients of $D$ ([Ch. 15](../part-5-transforms/15-fourier-series.md)): **truncating a Fourier series is the least-squares approximation** by a trigonometric polynomial of that order. For an ideal low-pass with cutoff $\omega_c$, $h_0=\omega_c/\pi$ and $h_k=\sin(k\omega_c)/(k\pi)$ — the windowed-sinc design with a rectangular window.

> [!example] Worked Example 23.35 — a 7-tap least-squares half-band filter
> Take $\omega_c=\pi/2$ and $M=3$.
>
> *Step 1 (coefficients):* $h_0=\tfrac12$, $h_1=\sin(\pi/2)/\pi=0.3183$, $h_2=\sin\pi/(2\pi)=0$, $h_3=\sin(3\pi/2)/(3\pi)=-0.1061$. The impulse response is $(h_3,h_2,h_1,h_0,h_1,h_2,h_3)=(-0.106,\ 0,\ 0.318,\ 0.5,\ 0.318,\ 0,\ -0.106)$.
>
> *Step 2 (response):* $A(\omega)=0.5+0.6366\cos\omega-0.2122\cos3\omega$. Check $A(0)=0.924$, $A(\pi/2)=0.5$ (half amplitude at the cutoff, by symmetry), $A(\pi)=0.076$.
>
> *Step 3 (error):* the minimum integrated squared error is $0.078$; the peak of $A$ is $1.100$ at $\omega=\pi/4$ — a $10\,\%$ overshoot. Raising $M$ to 10 lowers the integrated error to $0.032$ but the overshoot stays near $9\,\%$: the **Gibbs phenomenon**. Least squares minimises average error, not peak error.

```desmos-graph
left=0; right=3.1416; bottom=-0.25; top=1.25
---
y=0.5+0.6366\cos(x)-0.2122\cos(3x)|blue
y=0.5+\sum_{k=1}^{10}\frac{2\sin(k\pi/2)}{k\pi}\cos(kx)|green|dashed
y=1|red|0<x<1.5708
y=0|red|1.5708<x<3.1416
```

*Least-squares approximations of an ideal half-band low-pass (red) by 7-tap (blue) and 21-tap (green dashed) linear-phase FIR filters, against normalised frequency $\omega$ in rad/sample. More taps sharpen the transition, but the overshoot near the edge (Gibbs) does not shrink.*

**Improving the formulation.** Two changes to the *objective* do more than adding taps:
1. **Leave a transition band free.** Fit only on $[0,\omega_p]\cup[\omega_s,\pi]$ and solve the resulting discrete weighted least-squares problem on a frequency grid (§23.8). The jump disappears from the target, and with it most of the Gibbs overshoot.
2. **Change the norm.** Minimising the *peak* weighted error $\max_i\lvert A(\omega_i)-D(\omega_i)\rvert$ is an LP by the epigraph trick (§23.2); its solution is the **equiripple** (minimax) filter that the Parks–McClellan algorithm computes by a specialised exchange method.

> [!example] Worked Example 23.36 — least squares versus minimax for a 21-tap filter
> Specification: passband $[0,0.4\pi]$, stopband $[0.6\pi,\pi]$, $M=10$ (21 taps), 200 grid points in each band, equal weights.
>
> *Step 1 (least squares):* `numpy.linalg.lstsq` on the $400\times11$ system gives a peak error of $0.0268$, i.e. stopband peak $20\log_{10}0.0268=-31.4$ dB.
>
> *Step 2 (minimax LP):* minimise $t$ over $(\mathbf{h},t)$ subject to $-t\le A(\omega_i;\mathbf{h})-D(\omega_i)\le t$ for all 400 grid points — 12 variables, 800 linear inequalities. `linprog` returns $t=0.0114$: stopband peak $-38.9$ dB, about 7.5 dB better in the worst case, with the same number of taps.
>
> *Step 3 (moral):* the minimax filter has larger *average* error and smaller *peak* error than the least-squares one. Which is "optimal" depends on whether the specification is an energy criterion or a mask — the point of §23.1.

### A design-optimization checklist

1. **Variables.** What can actually be chosen? Use units of order one, or logarithms for positive quantities spanning decades.
2. **Objective.** What does "better" mean? Is it a proxy? If several goals conflict, plan a trade-off study (§23.9) rather than inventing weights.
3. **Constraints.** List every specification as $g\le0$ or $h=0$. Check feasibility first (Example 23.4). Remove redundant constraints.
4. **Classify.** Linear → LP; convex quadratic → QP; sum of squares → least squares; monomials → GP in log variables; otherwise check convexity with the toolkit of §23.7.
5. **Solve.** Use the matching solver. For nonconvex problems use several starting points, and seed them from a convex simplification.
6. **Verify.** Check KKT residuals, the sign of multipliers, and the duality gap. Perturb the data and re-solve: is the answer stable?
7. **Interpret.** Read multipliers as the price of each specification. Round to available values, re-check the specifications, and simulate the rounded design.

> [!warning] Common Pitfall — optimising the model, not the product
> An optimizer drives the design to wherever the model is most flattering — often the edge of its validity (a transistor model used outside its fitted bias range, a linear cost curve far beyond the tested load). Add constraints that keep the design inside the model's domain of validity, and verify the final design with a higher-fidelity model or measurement.

> [!tip] Engineering Note — where these problems appear in tools
> Optimal power flow and unit commitment in grid operation (LP/QP/MILP), automatic matching-network and filter synthesis in RF tools (nonlinear least squares and minimax), parameter extraction for device models (Levenberg–Marquardt), model predictive control running a QP every sample period, and analogue-circuit sizing by geometric programming are all direct applications of this chapter's methods.

## Summary

| Topic | Key result | Formula |
|---|---|---|
| Standard problem | minimise subject to inequalities and equalities | $\min f(\mathbf{x})$ s.t. $g_i(\mathbf{x})\le0,\ h_k(\mathbf{x})=0$ |
| Existence | continuous $f$ on a compact set, or coercive $f$ | Weierstrass |
| Epigraph trick | minimax / $L^1$ to LP | $\min t$ s.t. $r_i(\mathbf{x})\le t$ |
| LP primal–dual | weak/strong duality | $\max\mathbf{c}^\mathsf{T}\mathbf{x},\ \mathbf{A}\mathbf{x}\le\mathbf{b},\ \mathbf{x}\ge\mathbf{0}\ \leftrightarrow\ \min\mathbf{b}^\mathsf{T}\mathbf{y},\ \mathbf{A}^\mathsf{T}\mathbf{y}\ge\mathbf{c},\ \mathbf{y}\ge\mathbf{0}$ |
| Complementary slackness | spare capacity has zero price | $y_i(b_i-\mathbf{a}_i^\mathsf{T}\mathbf{x})=0$, $\mu_ig_i(\mathbf{x})=0$ |
| Unconstrained optimality | FONC, SOSC | $\nabla f(\mathbf{x}^*)=\mathbf{0}$, $\nabla^2f(\mathbf{x}^*)\succ0$ |
| KKT stationarity | constrained optimality | $\nabla f+\sum\mu_i\nabla g_i+\sum\lambda_k\nabla h_k=\mathbf{0}$, $\mu_i\ge0$ |
| Multiplier as sensitivity | shadow price | $\partial f^*/\partial b_i=-\mu_i$ for $g_i\le b_i$ |
| Gradient descent | stable step, quadratic rate | $\mathbf{x}_{k+1}=\mathbf{x}_k-\alpha\nabla f$, $0<\alpha<2/L$, $\rho=\dfrac{\kappa-1}{\kappa+1}$ |
| Exact line search (quadratic) | step length | $\alpha_k=\dfrac{\mathbf{g}_k^\mathsf{T}\mathbf{g}_k}{\mathbf{g}_k^\mathsf{T}\mathbf{A}\mathbf{g}_k}$ |
| Armijo condition | sufficient decrease | $f(\mathbf{x}+\alpha\mathbf{p})\le f(\mathbf{x})+c_1\alpha\nabla f^\mathsf{T}\mathbf{p}$ |
| Newton step | quadratic convergence | $\nabla^2f(\mathbf{x}_k)\,\mathbf{p}=-\nabla f(\mathbf{x}_k)$ |
| Secant condition | quasi-Newton | $\mathbf{H}_{k+1}\mathbf{y}_k=\mathbf{s}_k$ |
| BFGS update | inverse Hessian | $\mathbf{H}_{k+1}=(\mathbf{I}-\rho\mathbf{s}\mathbf{y}^\mathsf{T})\mathbf{H}_k(\mathbf{I}-\rho\mathbf{y}\mathbf{s}^\mathsf{T})+\rho\mathbf{s}\mathbf{s}^\mathsf{T}$, $\rho=1/\mathbf{y}^\mathsf{T}\mathbf{s}$ |
| Convex function | chord above graph; tangent below | $f(\mathbf{y})\ge f(\mathbf{x})+\nabla f(\mathbf{x})^\mathsf{T}(\mathbf{y}-\mathbf{x})$; $\nabla^2f\succeq0$ |
| Convex problem | local = global; KKT sufficient | Slater $\Rightarrow$ strong duality |
| Normal equations | linear least squares | $\mathbf{A}^\mathsf{T}\mathbf{A}\mathbf{x}=\mathbf{A}^\mathsf{T}\mathbf{b}$, $\mathbf{A}^\mathsf{T}\mathbf{r}^*=\mathbf{0}$ |
| Ridge regression | regularised LS | $(\mathbf{A}^\mathsf{T}\mathbf{A}+\lambda\mathbf{I})\mathbf{x}=\mathbf{A}^\mathsf{T}\mathbf{b}$, filter factors $\sigma_i^2/(\sigma_i^2+\lambda)$ |
| Gauss–Newton | nonlinear LS | $\mathbf{J}^\mathsf{T}\mathbf{J}\boldsymbol\delta=-\mathbf{J}^\mathsf{T}\mathbf{r}$ (LM: $+\tau\mathbf{I}$) |
| Weighted sum | Pareto point (convex fronts only) | $\min\sum w_if_i$, front slope $=-w_1/w_2$ |
| $\varepsilon$-constraint | any Pareto point | $\min f_1$ s.t. $f_i\le\varepsilon_i$, $\partial f_1^*/\partial\varepsilon_i=-\mu_i$ |
| Economic dispatch | equal incremental cost | $C_i'(P_i)=\lambda$ for units within limits |

## Further reading

- S. Boyd and L. Vandenberghe, *Convex Optimization*, Cambridge University Press — the standard text on convex sets, functions, duality and problem classes, with many engineering applications (freely available from the authors).
- J. Nocedal and S. J. Wright, *Numerical Optimization*, 2nd ed., Springer — line searches, trust regions, Newton, quasi-Newton (BFGS, L-BFGS), nonlinear least squares and constrained algorithms in depth.
- D. G. Luenberger and Y. Ye, *Linear and Nonlinear Programming*, Springer — simplex method, LP duality, and the convergence theory of gradient and Newton methods.
- A. J. Wood, B. F. Wollenberg and G. B. Sheblé, *Power Generation, Operation, and Control*, Wiley — economic dispatch, optimal power flow and their optimization formulations.
