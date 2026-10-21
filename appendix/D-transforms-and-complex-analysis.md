# Appendix D — Transforms, ODEs and Complex Analysis

[← Table of Contents](../00-front/table-of-contents.md)

> This appendix gives complete proofs of the results that Chapters 11 and 14–18 state with only a sketch: existence and uniqueness for initial-value problems, the analytic foundations of the Laplace, Fourier and Z-transforms, the convergence theory of Fourier series (Bessel, Parseval, Dirichlet, Gibbs), the sampling theorem, and the core of complex analysis (Cauchy–Riemann, Cauchy's theorem and integral formula, residues, Jordan's lemma).
>
> **Prerequisites.** Real analysis at the level of a first course: uniform convergence and the fact that a uniform limit of continuous functions is continuous, the Weierstrass M-test, Riemann (or Lebesgue) integration including improper integrals, and the interchange theorems of Fubini–Tonelli and dominated convergence. These are quoted, not proved. Green's theorem is taken from [Ch. 7](../part-2-calculus/07-vector-calculus.md); Laurent's theorem from [Ch. 18, §18.6](../part-5-transforms/18-complex-analysis.md#186-taylor-and-laurent-series).
>
> **Level of rigor.** Every hypothesis is stated precisely and every nontrivial step is justified. Where the cleanest general proof needs measure theory (Plancherel on all of $L^2$, the density of continuous functions in $L^2$), the proof is given in full for a large class that covers every engineering signal, and the extension is described and cited. Notation follows the [Notation and Conventions](../00-front/notation-and-conventions.md): $j^2=-1$, $\omega$ in rad/s, and the Fourier transform is $X(\omega)=\int_{-\infty}^{\infty}x(t)e^{-j\omega t}\,dt$ with inverse $x(t)=\frac{1}{2\pi}\int_{-\infty}^{\infty}X(\omega)e^{j\omega t}\,d\omega$, as in [Ch. 16](../part-5-transforms/16-fourier-transforms.md).

## D.1 Picard–Lindelöf existence and uniqueness

> [!definition] Lipschitz condition
> A function $f(t,y)$ defined on a set $R\subset\mathbb{R}^2$ is **Lipschitz in $y$** on $R$ with constant $K\ge0$ if
>
> $$
> \lvert f(t,y_1)-f(t,y_2)\rvert\le K\lvert y_1-y_2\rvert\qquad\text{for all }(t,y_1),(t,y_2)\in R .
> $$

> [!theorem] Picard–Lindelöf
> Let $R=\{(t,y):\lvert t-t_0\rvert\le a,\ \lvert y-y_0\rvert\le b\}$ with $a,b>0$. Suppose $f$ is continuous on $R$ and Lipschitz in $y$ on $R$ with constant $K$. Let $M=\max_R\lvert f\rvert$ (if $M=0$ put $h=a$) and $h=\min(a,\,b/M)$. Then the initial-value problem
>
> $$
> y'(t)=f(t,y(t)),\qquad y(t_0)=y_0
> $$
>
> has exactly one solution $y$ on $I=[t_0-h,\,t_0+h]$, and the Picard iterates $y_{k+1}(t)=y_0+\int_{t_0}^{t}f(s,y_k(s))\,ds$, $y_0(t)\equiv y_0$, converge to it uniformly on $I$.

The proof uses the contraction-mapping principle, which we prove first.

> [!theorem] Banach fixed-point theorem
> Let $(\mathcal{X},d)$ be a non-empty complete metric space and $\Phi:\mathcal{X}\to\mathcal{X}$ a **contraction**: $d(\Phi u,\Phi v)\le q\,d(u,v)$ for some $q<1$ and all $u,v$. Then $\Phi$ has exactly one fixed point $u^\star$, and for any $u_0$ the iterates $u_{k+1}=\Phi u_k$ converge to $u^\star$ with $d(u_k,u^\star)\le \dfrac{q^k}{1-q}\,d(u_1,u_0)$.

**Proof of the fixed-point theorem.** By induction $d(u_{k+1},u_k)\le q^k d(u_1,u_0)$. For $m>k$, the triangle inequality and the geometric series give

$$
d(u_m,u_k)\le\sum_{i=k}^{m-1}d(u_{i+1},u_i)\le d(u_1,u_0)\sum_{i=k}^{\infty}q^i=\frac{q^k}{1-q}\,d(u_1,u_0)\xrightarrow[k\to\infty]{}0 ,
$$

so $(u_k)$ is Cauchy and, by completeness, converges to some $u^\star$. A contraction is continuous, so $\Phi u^\star=\lim\Phi u_k=\lim u_{k+1}=u^\star$. Letting $m\to\infty$ in the display gives the error bound. If $u^\star,v^\star$ are both fixed, $d(u^\star,v^\star)=d(\Phi u^\star,\Phi v^\star)\le q\,d(u^\star,v^\star)$, which forces $d(u^\star,v^\star)=0$. ∎

**Proof of Picard–Lindelöf.**

*Step 1 — integral equation.* A continuous function $y$ on $I$ with graph in $R$ solves the IVP if and only if

$$
y(t)=y_0+\int_{t_0}^{t}f(s,y(s))\,ds\qquad(t\in I).
$$

Indeed, if $y$ solves the IVP then $y'$ is continuous (a composition of continuous functions) and integrating from $t_0$ gives the equation; conversely, the right-hand side is differentiable by the fundamental theorem of calculus, because $s\mapsto f(s,y(s))$ is continuous, and differentiating gives $y'=f(t,y)$ with $y(t_0)=y_0$.

*Step 2 — the space.* Let $\mathcal{X}$ be the set of continuous $y:I\to\mathbb{R}$ with $\lvert y(t)-y_0\rvert\le b$ on $I$. Equip it with the **weighted (Bielecki) metric**

$$
d(y,z)=\max_{t\in I}\,e^{-2K\lvert t-t_0\rvert}\,\lvert y(t)-z(t)\rvert .
$$

Since $e^{-2Kh}\le e^{-2K\lvert t-t_0\rvert}\le1$, this metric is equivalent to the sup metric: $e^{-2Kh}\lVert y-z\rVert_\infty\le d(y,z)\le\lVert y-z\rVert_\infty$. A sequence is Cauchy for $d$ iff it is uniformly Cauchy, and a uniform limit of continuous functions is continuous and still satisfies $\lvert y-y_0\rvert\le b$. Hence $(\mathcal{X},d)$ is complete. It is non-empty (it contains the constant $y_0$).

*Step 3 — $\Phi$ maps $\mathcal{X}$ into itself.* Define $(\Phi y)(t)=y_0+\int_{t_0}^{t}f(s,y(s))\,ds$. For $y\in\mathcal{X}$ the integrand is continuous, so $\Phi y$ is continuous, and

$$
\lvert(\Phi y)(t)-y_0\rvert\le M\lvert t-t_0\rvert\le Mh\le b .
$$

This is exactly why $h\le b/M$: the iterates cannot leave the rectangle where the bounds on $f$ hold.

*Step 4 — $\Phi$ is a contraction.* Take $t\ge t_0$ (the case $t<t_0$ is symmetric). By the Lipschitz condition and the definition of $d$,

$$
\begin{aligned}
\lvert(\Phi y)(t)-(\Phi z)(t)\rvert
&\le\int_{t_0}^{t}K\lvert y(s)-z(s)\rvert\,ds
\le K\,d(y,z)\int_{t_0}^{t}e^{2K(s-t_0)}\,ds\\
&=\frac{d(y,z)}{2}\big(e^{2K(t-t_0)}-1\big)
<\frac{1}{2}\,e^{2K(t-t_0)}\,d(y,z).
\end{aligned}
$$

Multiplying by $e^{-2K(t-t_0)}$ and maximising over $t$ gives $d(\Phi y,\Phi z)\le\tfrac12\,d(y,z)$. (If $K=0$ the integral is zero and the bound is trivial.)

*Step 5 — existence.* By the fixed-point theorem $\Phi$ has a fixed point $y\in\mathcal{X}$, which by Step 1 solves the IVP on $I$. The iterates $y_k=\Phi^k y_0$ converge to it in $d$, hence uniformly.

*Step 6 — uniqueness among all solutions.* The fixed-point theorem gives uniqueness only inside $\mathcal{X}$, so we must show that *every* solution $z$ on $I$ has its graph in $R$. Suppose $z$ leaves the strip $\lvert y-y_0\rvert\le b$ at some $t>t_0$ in $I$, and let $\tau=\inf\{t>t_0:\lvert z(t)-y_0\rvert>b\}$. On $[t_0,\tau]$ the graph is in $R$, so $\lvert z(\tau)-y_0\rvert\le M(\tau-t_0)$; by continuity $\lvert z(\tau)-y_0\rvert=b$. Hence $\tau-t_0\ge b/M\ge h$, so $\tau\ge t_0+h$ and $z$ does not leave the strip inside $I$ after all — a contradiction. The case $t<t_0$ is symmetric. So $z\in\mathcal{X}$ is a fixed point of $\Phi$, and $z=y$. ∎

**Remarks.**

1. *The Lipschitz condition is what buys uniqueness.* For $y'=3y^{2/3}$, $y(0)=0$, the function $f$ is continuous but not Lipschitz at $y=0$, and both $y\equiv0$ and $y=t^3$ are solutions (indeed infinitely many, $y=(t-c)^3$ for $t\ge c$ and $0$ before). Continuity alone still gives *existence* (Peano's theorem, proved by compactness rather than contraction), but not uniqueness.
2. *The result is local.* For $y'=y^2$, $y(0)=1$, $f$ is Lipschitz on every bounded rectangle, yet the solution $y=1/(1-t)$ blows up at $t=1$. The theorem's $h=\min(a,b/M)$ with $M=(1+b)^2$ is maximised at $b=1$, $a\ge 1/4$, giving $h=1/4$; repeated application extends the solution as long as it stays bounded, never past $t=1$.
3. *Why the weighted metric.* With the plain sup metric Step 4 gives the factor $Kh$, so the proof would need $h<1/K$ in addition. The weight $e^{-2K\lvert t-t_0\rvert}$ absorbs the exponential growth that Lipschitz feedback allows (the same growth that appears in Grönwall's inequality) and gives a contraction on the full interval.
4. *Systems.* The proof works verbatim for $\mathbf{y}'=\mathbf{f}(t,\mathbf{y})$ with $\lvert\cdot\rvert$ replaced by a vector norm; this covers every $n$th-order equation written as a first-order system, and hence the state equations of any lumped circuit. For a linear system $\mathbf{y}'=\mathbf{A}(t)\mathbf{y}+\mathbf{g}(t)$ with continuous coefficients, $f$ is globally Lipschitz and a solution exists on the whole interval of continuity ([Ch. 11, §11.9](../part-4-differential-equations/11-ordinary-differential-equations.md#119-systems-of-differential-equations)).
5. *Where it is used.* The theorem justifies the statement in [Ch. 11, §11.7](../part-4-differential-equations/11-ordinary-differential-equations.md#117-initial-value-problems) that an RLC circuit's response is fixed by its initial capacitor voltages and inductor currents, and underlies the error analysis of the numerical ODE solvers of [Ch. 22](../part-7-numerical-computational/22-numerical-methods.md).

```desmos-graph
left=-0.2; right=2; bottom=-0.2; top=7.5
---
y=e^{x}|black
y=1|grey|dashed
y=1+x|blue|dashed
y=1+x+\frac{x^{2}}{2}|green|dashed
y=1+x+\frac{x^{2}}{2}+\frac{x^{3}}{6}|orange|dashed
y=1+x+\frac{x^{2}}{2}+\frac{x^{3}}{6}+\frac{x^{4}}{24}|red|dashed
```

*Picard iterates $y_0,\dots,y_4$ for $y'=y$, $y(0)=1$ (dashed) approaching the fixed point $e^t$ (black); the error of $y_k$ is bounded by $\lvert t\rvert^{k+1}e^{\lvert t\rvert}/(k+1)!$, the factorial being the hallmark of the contraction argument.*

## D.2 Existence of the Laplace transform

> [!definition] Piecewise continuity and exponential order
> A function $x:[0,\infty)\to\mathbb{C}$ is **piecewise continuous** if every bounded interval $[0,\tau]$ contains only finitely many discontinuities of $x$, each a finite jump. It is of **exponential order $\alpha$** if there are constants $M>0$ and $T_0\ge0$ with $\lvert x(t)\rvert\le Me^{\alpha t}$ for all $t\ge T_0$.

> [!theorem] Existence of the Laplace transform
> If $x$ is piecewise continuous on $[0,\infty)$ and of exponential order $\alpha$, then
>
> $$
> X(s)=\int_0^\infty x(t)\,e^{-st}\,dt
> $$
>
> converges absolutely for every $s$ with $\operatorname{Re}\{s\}>\alpha$. Moreover the convergence is uniform on every half-plane $\operatorname{Re}\{s\}\ge\sigma_1>\alpha$, $X$ is analytic for $\operatorname{Re}\{s\}>\alpha$, and $X(s)\to0$ as $\operatorname{Re}\{s\}\to+\infty$.

**Proof.** Write $s=\sigma+j\omega$, so $\lvert e^{-st}\rvert=e^{-\sigma t}$.

*Absolute convergence.* On $[0,T_0]$ the function $x$ is piecewise continuous, hence bounded, say $\lvert x\rvert\le B$, and Riemann integrable; so $\int_0^{T_0}\lvert x(t)\rvert e^{-\sigma t}dt\le BT_0\max(1,e^{-\sigma T_0})<\infty$. On $[T_0,\infty)$,

$$
\int_{T_0}^{\tau}\lvert x(t)\rvert e^{-\sigma t}\,dt\le M\int_{T_0}^{\tau}e^{-(\sigma-\alpha)t}\,dt\le\frac{M\,e^{-(\sigma-\alpha)T_0}}{\sigma-\alpha}\qquad(\sigma>\alpha),
$$

a bound independent of $\tau$. A non-negative integrand with bounded partial integrals has a convergent improper integral, so $\int_0^\infty\lvert x\rvert e^{-\sigma t}dt<\infty$, and absolute convergence implies convergence.

*Uniformity.* If $\sigma\ge\sigma_1>\alpha$ then $\lvert x(t)e^{-st}\rvert\le\lvert x(t)\rvert e^{-\sigma_1 t}$, an integrable majorant independent of $s$. Hence the tails satisfy $\big\lvert\int_\tau^\infty x e^{-st}dt\big\rvert\le\int_\tau^\infty\lvert x\rvert e^{-\sigma_1t}dt\to0$ uniformly in $s$ (the Weierstrass M-test for integrals).

*Analyticity.* Each truncation $X_\tau(s)=\int_0^\tau x(t)e^{-st}dt$ is entire: differentiation under the integral sign over a bounded interval is permitted because $\partial_s(x e^{-st})=-t\,x\,e^{-st}$ is continuous in $s$ and bounded on compact $s$-sets. By the uniformity just proved, $X_\tau\to X$ uniformly on each half-plane $\operatorname{Re}\{s\}\ge\sigma_1$, and a locally uniform limit of analytic functions is analytic (each $X_\tau$ integrates to zero around every triangle by §D.14, the limit passes inside the contour integral, and Morera's theorem applies).

*Decay.* For $\sigma>\max(\alpha,0)$,

$$
\lvert X(s)\rvert\le B\int_0^{T_0}e^{-\sigma t}dt+M\int_{T_0}^\infty e^{-(\sigma-\alpha)t}dt\le\frac{B}{\sigma}+\frac{M}{\sigma-\alpha}\xrightarrow[\sigma\to\infty]{}0 .
$$

∎

**Remarks.**

1. *Sufficient, not necessary.* $x(t)=t^{-1/2}$ is not piecewise continuous at $0$ (it is unbounded), yet its transform $\sqrt{\pi/s}$ exists because the singularity is integrable. The theorem extends to any $x$ that is absolutely integrable on each $[0,\tau]$ and of exponential order at infinity, with the same proof.
2. *Growth control is essential.* $x(t)=e^{t^2}$ has no Laplace transform for any $s$: for real $s$ the integrand $e^{t(t-s)}\to\infty$. More generally, the set of $s$ where the integral converges absolutely is a right half-plane $\operatorname{Re}\{s\}>\sigma_a$ (possibly with its boundary line), because absolute convergence at $\sigma$ implies it at every $\sigma'>\sigma$ by comparison.
3. *Consequence for rational transforms.* The decay $X(s)\to0$ shows that a polynomial in $s$, or $s/(s+1)$, cannot be the transform of an ordinary function; transforms of ordinary functions with rational $X$ are strictly proper. Constants and higher powers of $s$ correspond to impulses and their derivatives at $t=0$.
4. Used in [Ch. 14, §14.1](../part-5-transforms/14-laplace-transforms.md#141-definition-and-existence).

## D.3 The Laplace derivative theorem

> [!theorem] Transform of a derivative
> Let $x$ be continuous on $[0,\infty)$ and of exponential order $\alpha$, and let $\dot x$ be piecewise continuous on $[0,\infty)$ (at finitely many points in each bounded interval $\dot x$ may fail to exist; there $x$ has a corner). Then $\mathcal{L}\{\dot x\}$ exists for $\operatorname{Re}\{s\}>\alpha$ and
>
> $$
> \mathcal{L}\{\dot x(t)\}=sX(s)-x(0).
> $$
>
> If $x,\dot x,\dots,x^{(n-1)}$ are continuous and of exponential order $\alpha$ and $x^{(n)}$ is piecewise continuous, then
>
> $$
> \mathcal{L}\{x^{(n)}\}=s^nX(s)-s^{n-1}x(0)-s^{n-2}\dot x(0)-\dots-x^{(n-1)}(0).
> $$

**Proof.** Fix $s$ with $\sigma=\operatorname{Re}\{s\}>\alpha$ and $\tau>0$. Let $0=t_0<t_1<\dots<t_m=\tau$ include all points of $[0,\tau]$ where $\dot x$ is discontinuous or undefined. On each $[t_{i-1},t_i]$ the function $x$ is continuously differentiable (with one-sided derivatives at the ends), so integration by parts is valid:

$$
\int_{t_{i-1}}^{t_i}\dot x(t)e^{-st}\,dt=\Big[x(t)e^{-st}\Big]_{t_{i-1}}^{t_i}+s\int_{t_{i-1}}^{t_i}x(t)e^{-st}\,dt .
$$

Summing over $i$, the boundary terms telescope **because $x$ is continuous** at each interior $t_i$:

$$
\int_0^\tau\dot x(t)e^{-st}\,dt=x(\tau)e^{-s\tau}-x(0)+s\int_0^\tau x(t)e^{-st}\,dt .
$$

Now let $\tau\to\infty$. For $\tau\ge T_0$, $\lvert x(\tau)e^{-s\tau}\rvert\le Me^{-(\sigma-\alpha)\tau}\to0$, and the last integral tends to $X(s)$ by §D.2. The right-hand side therefore has the limit $sX(s)-x(0)$, so the improper integral on the left converges and equals it. The $n$th-order formula follows by induction: applying the first-order result to $x^{(n-1)}$ gives $\mathcal{L}\{x^{(n)}\}=s\,\mathcal{L}\{x^{(n-1)}\}-x^{(n-1)}(0)$, and substituting the formula for $\mathcal{L}\{x^{(n-1)}\}$ yields the stated sum. ∎

**Remarks.**

1. *Continuity matters.* If $x$ has a jump of size $J$ at $t_1>0$, the boundary terms fail to telescope and an extra term $-Je^{-st_1}$ appears: the classical derivative misses the impulse $J\delta(t-t_1)$ that the generalized derivative contains. For a jump at the origin the book uses the $0^-$ convention, $\mathcal{L}\{\dot x\}=sX(s)-x(0^-)$, in which $\dot x$ includes the impulse $[x(0^+)-x(0^-)]\delta(t)$; the proof above with $x(0)$ replaced by $x(0^+)$ gives the $0^+$ form, and the two differ exactly by that impulse ([Ch. 14, §14.2](../part-5-transforms/14-laplace-transforms.md#142-properties)).
2. *$\dot x$ need not be of exponential order.* The proof never bounds $\dot x$; the existence of $\mathcal{L}\{\dot x\}$, possibly as a conditionally convergent integral, is a conclusion. Example: $x(t)=\sin(e^{t^2})$ is bounded, so $\mathcal{L}\{x\}$ exists for $\operatorname{Re}\{s\}>0$, while $\dot x(t)=2te^{t^2}\cos(e^{t^2})$ grows faster than any exponential; nevertheless $\int_0^\infty\dot x\, e^{-st}dt$ converges for $\operatorname{Re}\{s\}>0$.
3. This is the result that turns a linear constant-coefficient ODE with initial conditions into an algebraic equation in $s$ ([Ch. 14, §14.8](../part-5-transforms/14-laplace-transforms.md#148-solving-differential-equations)).

## D.4 The Laplace convolution theorem

> [!theorem] Convolution theorem (Laplace)
> Let $f$ and $g$ be piecewise continuous on $[0,\infty)$ and of exponential order $\alpha$, and define
>
> $$
> (f*g)(t)=\int_0^t f(\tau)\,g(t-\tau)\,d\tau\qquad(t\ge0).
> $$
>
> Then $f*g$ is continuous, of exponential order $\alpha+\varepsilon$ for every $\varepsilon>0$, and for $\operatorname{Re}\{s\}>\alpha$
>
> $$
> \mathcal{L}\{f*g\}(s)=F(s)\,G(s).
> $$

**Proof.** *Growth of $f*g$.* Enlarging $M$ if necessary, we may take $\lvert f(t)\rvert,\lvert g(t)\rvert\le Me^{\alpha t}$ for **all** $t\ge0$ (on $[0,T_0]$ both functions are bounded and $e^{\alpha t}$ is bounded below). Then

$$
\lvert(f*g)(t)\rvert\le M^2\int_0^te^{\alpha\tau}e^{\alpha(t-\tau)}d\tau=M^2\,t\,e^{\alpha t}\le\frac{M^2}{\varepsilon e}\,e^{(\alpha+\varepsilon)t},
$$

using $te^{-\varepsilon t}\le1/(\varepsilon e)$. For continuity, fix $\tau_0>0$ and $t,t'\in[0,\tau_0]$; extending $g$ by $0$ on $t<0$, $(f*g)(t')-(f*g)(t)=\int_0^{\tau_0}f(\tau)[g(t'-\tau)-g(t-\tau)]d\tau$. As $t'\to t$ the bracket tends to $0$ for all but finitely many $\tau$ (the points where $t-\tau$ is a jump of $g$) and is bounded, so the integral tends to $0$ by dominated convergence. Hence $\mathcal{L}\{f*g\}$ exists for $\operatorname{Re}\{s\}>\alpha$ (choose $\varepsilon<\operatorname{Re}\{s\}-\alpha$).

*The product formula.* Fix $s$ with $\sigma=\operatorname{Re}\{s\}>\alpha$. With $g=0$ on $t<0$, $(f*g)(t)=\int_0^\infty f(\tau)g(t-\tau)d\tau$, so

$$
\mathcal{L}\{f*g\}(s)=\int_0^\infty\!\!\int_0^\infty f(\tau)\,g(t-\tau)\,e^{-st}\,d\tau\,dt .
$$

For the absolute integrand, Tonelli's theorem (iterated integrals of a non-negative function agree in either order) gives

$$
\int_0^\infty\!\!\int_0^\infty\lvert f(\tau)\rvert\,\lvert g(t-\tau)\rvert\,e^{-\sigma t}\,dt\,d\tau
=\int_0^\infty\lvert f(\tau)\rvert e^{-\sigma\tau}\,d\tau\cdot\int_0^\infty\lvert g(u)\rvert e^{-\sigma u}\,du<\infty,
$$

where we substituted $u=t-\tau$ in the inner integral and used the absolute convergence of §D.2 for both factors. So Fubini's theorem allows the order of integration to be exchanged in the signed integral, and the same substitution gives

$$
\mathcal{L}\{f*g\}(s)=\int_0^\infty f(\tau)e^{-s\tau}\left(\int_0^\infty g(u)e^{-su}\,du\right)d\tau=F(s)\,G(s).
$$

∎

**Remarks.**

1. *Region of validity.* The formula holds on the intersection of the half-planes of absolute convergence of $F$ and $G$. Absolute convergence is exactly the hypothesis that licenses Fubini; with only conditional convergence the exchange of integrals can fail.
2. *Why it matters.* With $g=h$ the impulse response of a causal LTI system and $f$ its input, the theorem says the output transform is $H(s)F(s)$ — the definition of the transfer function. Geometrically, the region $0\le\tau\le t$ of the $(\tau,t)$ plane can be swept by vertical lines (convolve first) or horizontal lines (transform first). See [Ch. 14, §14.6](../part-5-transforms/14-laplace-transforms.md#146-convolution-theorem).

## D.5 Initial- and final-value theorems

> [!theorem] Initial-value theorem (IVT)
> Let $x$ be continuous on $(0,\infty)$ with a finite right limit $x(0^+)$ and of exponential order, and let $\dot x$ be piecewise continuous and of exponential order $\alpha$ on $(0,\infty)$. Then, along real $s\to+\infty$,
>
> $$
> \lim_{s\to+\infty}sX(s)=x(0^+).
> $$

> [!theorem] Final-value theorem (FVT)
> **(a) Abelian form.** If $x$ is piecewise continuous on $[0,\infty)$ and the limit $x(\infty)=\lim_{t\to\infty}x(t)$ exists and is finite, then $X(s)$ exists for $\operatorname{Re}\{s\}>0$ and
>
> $$
> \lim_{s\to0^+}sX(s)=x(\infty)\qquad(s\text{ real}).
> $$
>
> **(b) Rational form.** If $X(s)$ is a strictly proper rational function and **every pole of $sX(s)$ lies in the open left half-plane** $\operatorname{Re}\{s\}<0$, then $x(\infty)$ exists and equals $\lim_{s\to0}sX(s)$.

**Proof of the IVT.** Redefine $x(0)=x(0^+)$; then $x$ is continuous on $[0,\infty)$ and §D.3 gives, for real $s>\alpha$,

$$
sX(s)-x(0^+)=\int_0^\infty\dot x(t)\,e^{-st}\,dt=\mathcal{L}\{\dot x\}(s).
$$

Because $\dot x$ is piecewise continuous and of exponential order, §D.2 applies to $\dot x$ itself and gives $\mathcal{L}\{\dot x\}(s)\to0$ as $s\to+\infty$. ∎

**Proof of FVT (a).** Write $L=x(\infty)$. Since $x$ converges it is bounded for large $t$, and it is bounded on any $[0,T]$ by piecewise continuity; so it is of exponential order $0$ and $X(s)$ exists for $\operatorname{Re}\{s\}>0$. Using $s\int_0^\infty e^{-st}dt=1$ for $s>0$,

$$
sX(s)-L=s\int_0^\infty\big(x(t)-L\big)e^{-st}\,dt .
$$

Given $\varepsilon>0$ pick $T$ with $\lvert x(t)-L\rvert<\varepsilon$ for $t>T$, and let $B=\sup_{[0,T]}\lvert x-L\rvert<\infty$. Then for $s>0$

$$
\lvert sX(s)-L\rvert\le s\int_0^TB\,dt+s\int_T^\infty\varepsilon\, e^{-st}dt\le sBT+\varepsilon .
$$

Letting $s\to0^+$ gives $\limsup\lvert sX(s)-L\rvert\le\varepsilon$, and $\varepsilon$ was arbitrary. ∎

**Proof of FVT (b).** Let $X=N/D$ be strictly proper with distinct poles $p_1,\dots,p_r$ of multiplicities $m_1,\dots,m_r$. The partial-fraction expansion ([Ch. 14, §14.5](../part-5-transforms/14-laplace-transforms.md#145-partial-fractions)) and uniqueness of the inverse transform give

$$
x(t)=\sum_{i=1}^{r}\sum_{k=1}^{m_i}\frac{A_{ik}}{(k-1)!}\,t^{k-1}e^{p_it}\qquad(t>0).
$$

The poles of $sX(s)$ are those of $X$, except that a pole of $X$ at $0$ has its order lowered by one (and disappears if it was simple). The hypothesis therefore says: every $p_i\ne0$ has $\operatorname{Re}\{p_i\}<0$, and $0$ is at most a *simple* pole of $X$. Each term with $\operatorname{Re}\{p_i\}<0$ tends to zero, since $t^{k-1}e^{\operatorname{Re}\{p_i\}t}\to0$. If $0$ is a simple pole, the remaining term is the constant $A_0=\lim_{s\to0}sX(s)$ (the residue of $X$ at $0$); if $0$ is not a pole, there is no such term and $\lim_{s\to0}sX(s)=0\cdot X(0)=0$. In either case $x(t)\to\lim_{s\to0}sX(s)$. ∎

**Remarks — when the FVT fails.**

1. *Poles on the imaginary axis.* $x(t)=\cos t$ has $X(s)=s/(s^2+1)$ and $sX(s)=s^2/(s^2+1)\to0$ as $s\to0$, but $\cos t$ has no limit. Hypothesis (b) fails because $sX(s)$ has poles at $\pm j$. The Abelian form (a) is not contradicted: it assumes the limit exists.
2. *Poles in the right half-plane.* $x(t)=e^{t}$ has $sX(s)=s/(s-1)\to0$ as $s\to0$, a meaningless "final value" for a signal that diverges. The Laplace integral does not converge near $s=0$; the formula $s/(s-1)$ there is an analytic continuation, which is why blind substitution misleads.
3. *Double pole at the origin.* $x(t)=t$ gives $sX(s)=1/s\to\infty$, correctly signalling divergence; hypothesis (b) fails because $sX(s)$ has a pole at $0$.
4. *The converse of (a) is false.* That $sX(s)$ has a limit as $s\to0^+$ does not imply that $x(t)$ converges (Remark 1). Converse results need extra "Tauberian" hypotheses — for instance $t\,\dot x(t)$ bounded below (Hardy–Littlewood); for rational $X$, hypothesis (b) plays this role.
5. *IVT hypotheses.* If $x$ contains an impulse at $0$ then $X$ is not strictly proper (e.g. $X=1+1/(s+1)$) and $sX(s)\to\infty$; the IVT applies to the regular part only. The limit in the IVT may be taken along any path with $\operatorname{Re}\{s\}\to\infty$, by the uniform estimate of §D.2.
6. Used in [Ch. 14, §14.7](../part-5-transforms/14-laplace-transforms.md#147-initial--and-final-value-theorems) for steady-state errors of control loops and DC gains.

```desmos-graph
left=-0.5; right=12; bottom=-1.5; top=2
---
y=1-e^{-x}|blue|0<x
y=1|blue|dashed
y=\cos(x)|red|0<x
```

*Blue: $x=1-e^{-t}$, $X=1/[s(s+1)]$, whose $sX$ has its only pole at $-1$, so the FVT correctly gives $x(\infty)=1$. Red: $\cos t$, for which $sX(s)\to0$ although no final value exists — the poles $\pm j$ violate the hypothesis.*

## D.6 Fourier coefficients from orthogonality; Bessel's inequality

Throughout §§D.6–D.9, $x$ is $T_0$-periodic, $\omega_0=2\pi/T_0$, and $\int_{T_0}$ denotes integration over any interval of length $T_0$. For functions square-integrable over a period we use the normalised inner product and norm of [Ch. 15, §15.1](../part-5-transforms/15-fourier-series.md#151-periodic-function-decomposition),

$$
\langle x,y\rangle=\frac1{T_0}\int_{T_0}x(t)\,y^*(t)\,dt,\qquad\lVert x\rVert^2=\frac1{T_0}\int_{T_0}\lvert x(t)\rvert^2dt ,
$$

so $\lVert x\rVert^2$ is the mean-square value (the power) of $x$. Write $\varphi_n(t)=e^{jn\omega_0t}$.

> [!theorem] Orthonormality and the coefficient formula
> **(i)** $\langle\varphi_m,\varphi_n\rangle=\delta_{mn}$ for all integers $m,n$.
>
> **(ii)** If a series $\sum_{n=-\infty}^{\infty}d_n\varphi_n(t)$ converges **uniformly** on a period to $x(t)$, then necessarily
>
> $$
> d_n=c_n:=\langle x,\varphi_n\rangle=\frac1{T_0}\int_{T_0}x(t)\,e^{-jn\omega_0t}\,dt ,
> $$
>
> and for real $x$ the trigonometric coefficients are $a_n=2\operatorname{Re}\{c_n\}=\frac{2}{T_0}\int_{T_0}x\cos n\omega_0t\,dt$, $b_n=-2\operatorname{Im}\{c_n\}=\frac{2}{T_0}\int_{T_0}x\sin n\omega_0t\,dt$.

> [!theorem] Bessel's inequality
> If $\int_{T_0}\lvert x\rvert^2dt<\infty$ and $S_N=\sum_{n=-N}^{N}c_n\varphi_n$, then for every $N$
>
> $$
> \lVert x-S_N\rVert^2=\lVert x\rVert^2-\sum_{n=-N}^{N}\lvert c_n\rvert^2\ \ge0,\qquad\text{hence}\qquad\sum_{n=-\infty}^{\infty}\lvert c_n\rvert^2\le\frac1{T_0}\int_{T_0}\lvert x(t)\rvert^2dt .
> $$
>
> In particular $c_n\to0$ as $\lvert n\rvert\to\infty$. Moreover $S_N$ is the unique trigonometric polynomial of degree $\le N$ closest to $x$ in mean square.

**Proof.** (i) For $m=n$ the integrand is $1$. For $k=m-n\ne0$, $\int_{t_1}^{t_1+T_0}e^{jk\omega_0t}dt=e^{jk\omega_0t_1}\big(e^{j2\pi k}-1\big)/(jk\omega_0)=0$.

(ii) Multiply $x(t)=\sum_md_m\varphi_m(t)$ by the bounded function $\varphi_n^*(t)$; the product series still converges uniformly, so it may be integrated term by term (the integral of a uniform limit over a bounded interval is the limit of the integrals). By (i), $\langle x,\varphi_n\rangle=\sum_md_m\langle\varphi_m,\varphi_n\rangle=d_n$. The real formulas follow from $\cos n\omega_0t=\tfrac12(\varphi_n+\varphi_{-n})$ and $\sin n\omega_0t=\tfrac1{2j}(\varphi_n-\varphi_{-n})$.

*Bessel.* Let $p=\sum_{|n|\le N}d_n\varphi_n$ be any trigonometric polynomial. Expanding the inner product and using (i),

$$
\begin{aligned}
\lVert x-p\rVert^2&=\lVert x\rVert^2-\sum_{|n|\le N}\big(d_n^*c_n+d_nc_n^*\big)+\sum_{|n|\le N}\lvert d_n\rvert^2\\
&=\lVert x\rVert^2-\sum_{|n|\le N}\lvert c_n\rvert^2+\sum_{|n|\le N}\lvert d_n-c_n\rvert^2 ,
\end{aligned}
$$

the second line by completing the square, $\lvert d_n-c_n\rvert^2=\lvert d_n\rvert^2-d_n^*c_n-d_nc_n^*+\lvert c_n\rvert^2$. The last sum is $\ge0$ and vanishes only for $d_n=c_n$; this proves the best-approximation property. Taking $d_n=c_n$ gives the identity $\lVert x-S_N\rVert^2=\lVert x\rVert^2-\sum_{|n|\le N}\lvert c_n\rvert^2$, whose left side is non-negative. The partial sums of the non-negative series $\sum\lvert c_n\rvert^2$ are therefore bounded by $\lVert x\rVert^2$, so the series converges with sum $\le\lVert x\rVert^2$, and its terms tend to zero. ∎

**Remarks.**

1. Part (ii) explains the *form* of the coefficient formula; the hypothesis of uniform convergence is there only so that the term-by-term integration is legitimate. For an arbitrary square-integrable $x$ we *define* $c_n$ by the formula and ask afterwards whether the series represents $x$ — the subject of §§D.7–D.8.
2. Bessel's inequality is a statement about geometry, not Fourier analysis: it holds for any orthonormal set in any inner-product space, because $S_N$ is the orthogonal projection of $x$ on the span of $\varphi_{-N},\dots,\varphi_N$ and projection cannot increase length. What is special about the trigonometric system is that equality holds (completeness), proved next.
3. The decay $c_n\to0$ for square-integrable $x$ is the Fourier-series form of the Riemann–Lebesgue lemma; §D.8 proves it for merely absolutely integrable $x$.
4. Used in [Ch. 15, §15.1](../part-5-transforms/15-fourier-series.md#151-periodic-function-decomposition) and [§15.7](../part-5-transforms/15-fourier-series.md#157-fourier-coefficients).

## D.7 Parseval's theorem for Fourier series

> [!theorem] Parseval's theorem
> If $x$ is $T_0$-periodic and $\int_{T_0}\lvert x\rvert^2dt<\infty$, then $\lVert x-S_N\rVert\to0$ as $N\to\infty$ and
>
> $$
> \frac1{T_0}\int_{T_0}\lvert x(t)\rvert^2\,dt=\sum_{n=-\infty}^{\infty}\lvert c_n\rvert^2 .
> $$
>
> More generally $\langle x,y\rangle=\sum_nc_n\,d_n^*$ for two such functions with coefficients $c_n$, $d_n$.

By the identity in Bessel's inequality, Parseval's equality is *equivalent* to $\lVert x-S_N\rVert\to0$. We prove that by approximating $x$ first by a continuous function and then by Fejér means, which converge uniformly.

> [!definition] Fejér means and Fejér kernel
> Take $T_0=2\pi$ (the general case follows by the change of variable $t\mapsto\omega_0t$). The **Fejér mean** is the average of the first $N+1$ partial sums, $\sigma_N=\frac{1}{N+1}\sum_{k=0}^{N}S_k$. Since $S_k(t)=\frac1{2\pi}\int_{-\pi}^{\pi}x(t-u)D_k(u)\,du$ with the Dirichlet kernel $D_k(u)=\sum_{n=-k}^{k}e^{jnu}$ (§D.8),
>
> $$
> \sigma_N(t)=\frac1{2\pi}\int_{-\pi}^{\pi}x(t-u)\,F_N(u)\,du,\qquad F_N(u)=\frac1{N+1}\sum_{k=0}^{N}D_k(u)=\frac{1}{N+1}\left(\frac{\sin\big((N+1)u/2\big)}{\sin(u/2)}\right)^{2}.
> $$

**Lemma (closed form and properties of $F_N$).** (a) the closed form above holds; (b) $F_N\ge0$; (c) $\frac1{2\pi}\int_{-\pi}^{\pi}F_N=1$; (d) for $0<\delta\le\lvert u\rvert\le\pi$, $F_N(u)\le\dfrac{1}{(N+1)\sin^2(\delta/2)}$.

*Proof.* From §D.8, $D_k(u)=\sin\big((k+\tfrac12)u\big)/\sin(u/2)$. Using $2\sin(u/2)\sin\big((k+\tfrac12)u\big)=\cos ku-\cos(k+1)u$ and telescoping,

$$
\sum_{k=0}^{N}D_k(u)=\frac{1-\cos(N+1)u}{2\sin^2(u/2)}=\frac{\sin^2\big((N+1)u/2\big)}{\sin^2(u/2)},
$$

which is (a), and (b) is evident. Each $D_k$ has mean value $1$ over a period (only the $n=0$ term survives integration), so (c) follows. For (d), the numerator is at most $1$ and $\sin^2(u/2)\ge\sin^2(\delta/2)$ on $\delta\le\lvert u\rvert\le\pi$. $\square$

**Fejér's theorem.** *If $x$ is continuous and $2\pi$-periodic, then $\sigma_N\to x$ uniformly.*

*Proof.* By (c), $\sigma_N(t)-x(t)=\frac1{2\pi}\int_{-\pi}^{\pi}\big[x(t-u)-x(t)\big]F_N(u)\,du$. A continuous periodic function is bounded, $\lvert x\rvert\le B$, and uniformly continuous: given $\varepsilon>0$ there is $\delta>0$ with $\lvert x(t-u)-x(t)\rvert<\varepsilon$ whenever $\lvert u\rvert<\delta$, for all $t$. Splitting the integral and using (b), (c), (d),

$$
\lvert\sigma_N(t)-x(t)\rvert\le\frac{\varepsilon}{2\pi}\int_{\lvert u\rvert<\delta}F_N+\frac{2B}{2\pi}\int_{\delta\le\lvert u\rvert\le\pi}F_N\le\varepsilon+\frac{2B}{(N+1)\sin^2(\delta/2)} ,
$$

uniformly in $t$. The last term is $<\varepsilon$ for $N$ large. $\square$

**Proof of Parseval's theorem.**

*Step 1 — continuous $x$.* $\sigma_N$ is a trigonometric polynomial of degree $\le N$. By the best-approximation property of §D.6 and Fejér's theorem,

$$
\lVert x-S_N\rVert\le\lVert x-\sigma_N\rVert\le\max_t\lvert x(t)-\sigma_N(t)\rvert\xrightarrow[N\to\infty]{}0 ,
$$

the middle inequality because a mean-square value never exceeds the square of the maximum.

*Step 2 — general square-integrable $x$.* Continuous periodic functions are dense in mean square: for every $\varepsilon>0$ there is a continuous periodic $g$ with $\lVert x-g\rVert<\varepsilon$. (For Riemann-integrable $x$ this is elementary: choose a step function $\psi$ from a lower Darboux sum with $\frac1{T_0}\int\lvert x-\psi\rvert$ small and $\lvert\psi\rvert\le\sup\lvert x\rvert$, so that $\lVert x-\psi\rVert^2\le2\sup\lvert x\rvert\cdot\frac1{T_0}\int\lvert x-\psi\rvert$ is small, then replace each jump of $\psi$ by a steep linear ramp, which changes $\lVert\cdot\rVert$ by an arbitrarily small amount; unbounded $x$ are first truncated at a large level. For Lebesgue-square-integrable $x$ it is a standard theorem of measure theory.) Let $S_N^g$ be the partial sums of $g$. Since $x\mapsto S_N$ is an orthogonal projection, $\lVert S_N-S_N^g\rVert\le\lVert x-g\rVert$. Then

$$
\lVert x-S_N\rVert\le\lVert x-g\rVert+\lVert g-S_N^g\rVert+\lVert S_N^g-S_N\rVert<2\varepsilon+\lVert g-S_N^g\rVert ,
$$

and the last term tends to zero by Step 1. Hence $\limsup_N\lVert x-S_N\rVert\le2\varepsilon$ for every $\varepsilon$.

*Step 3 — the identity.* By §D.6, $\lVert x\rVert^2-\sum_{|n|\le N}\lvert c_n\rvert^2=\lVert x-S_N\rVert^2\to0$.

*Step 4 — the inner-product form.* Apply the identity to $x+y$, $x-y$, $x+jy$, $x-jy$ and combine with the polarization identity $\langle x,y\rangle=\frac14\sum_{k=0}^{3}j^k\lVert x+j^ky\rVert^2$; the coefficient sequences combine the same way. ∎

**Remarks.**

1. *Physical reading.* The left side is the average power of $x$; the right side distributes it over harmonics, $\lvert c_n\rvert^2+\lvert c_{-n}\rvert^2$ being the power at frequency $n\omega_0$. For real $x$ this is $\big(\tfrac{a_0}{2}\big)^2+\tfrac12\sum_{n\ge1}(a_n^2+b_n^2)$. Check: the unit square wave ($\pm1$) has $b_n=4/(n\pi)$ for odd $n$ and mean-square value $1$; numerically $\tfrac12\sum_{n\text{ odd}}16/(n\pi)^2=1.0000$ (summing $10^5$ terms gives $0.999998$), equivalently $\sum_{n\text{ odd}}1/n^2=\pi^2/8$.
2. *What Parseval does not say.* Mean-square convergence does not give convergence at any particular point: there are continuous functions whose Fourier series diverge at a point (du Bois-Reymond, 1873). Pointwise statements need the extra regularity of §D.8.
3. *Fejér versus Dirichlet.* $F_N$ is non-negative while $D_N$ oscillates in sign with $\frac1{2\pi}\int\lvert D_N\rvert\sim\frac{4}{\pi^2}\ln N\to\infty$; this is exactly why Cesàro averaging cures the Gibbs overshoot ([Ch. 15, §15.6](../part-5-transforms/15-fourier-series.md#156-convergence), "Taming Gibbs").
4. Used for RMS values, THD and power factor in [Ch. 15, §15.8](../part-5-transforms/15-fourier-series.md#158-power-and-energy-relationships), and in §D.11 below.

## D.8 Dirichlet's pointwise convergence theorem

> [!definition] Dirichlet kernel
> For $T_0=2\pi$ and $N\ge0$,
>
> $$
> D_N(u)=\sum_{n=-N}^{N}e^{jnu}=\frac{\sin\big((N+\tfrac12)u\big)}{\sin(u/2)}\quad(u\notin2\pi\mathbb{Z}),\qquad D_N(2\pi k)=2N+1 .
> $$
>
> The closed form follows from the geometric sum: $(e^{ju/2}-e^{-ju/2})D_N(u)=e^{j(N+\frac12)u}-e^{-j(N+\frac12)u}$. $D_N$ is even, $2\pi$-periodic, and $\frac1{2\pi}\int_{-\pi}^{\pi}D_N(u)\,du=1$, so $\frac1{2\pi}\int_0^{\pi}D_N=\frac12$ by evenness.

> [!theorem] Riemann–Lebesgue lemma
> If $g$ is absolutely integrable on an interval $(a,b)$ (bounded or not, Riemann-improper or Lebesgue), then
>
> $$
> \lim_{\lambda\to\pm\infty}\int_a^b g(u)\,e^{j\lambda u}\,du=0,
> $$
>
> and hence also $\int_a^bg(u)\sin\lambda u\,du\to0$ and $\int_a^bg(u)\cos\lambda u\,du\to0$.

**Proof of Riemann–Lebesgue.** *Step functions.* For the indicator of an interval $[\alpha,\beta]$, $\big\lvert\int_\alpha^\beta e^{j\lambda u}du\big\rvert=\lvert e^{j\lambda\beta}-e^{j\lambda\alpha}\rvert/\lvert\lambda\rvert\le2/\lvert\lambda\rvert\to0$; by linearity the result holds for every step function (finite combination of indicators).

*Approximation.* Given $\varepsilon>0$ there is a step function $\psi$ with $\int_a^b\lvert g-\psi\rvert<\varepsilon$. (For a Riemann-integrable $g$ on a bounded interval take $\psi$ from a lower Darboux sum close to the integral, treating real and imaginary parts separately; for an improper integral first cut off a neighbourhood of each singularity and of infinity where $\int\lvert g\rvert<\varepsilon/2$; for Lebesgue-integrable $g$ density of step functions is standard.) Then

$$
\left\lvert\int_a^bg\,e^{j\lambda u}du\right\rvert\le\int_a^b\lvert g-\psi\rvert\,du+\left\lvert\int_a^b\psi\,e^{j\lambda u}du\right\rvert<\varepsilon+\left\lvert\int_a^b\psi\,e^{j\lambda u}du\right\rvert ,
$$

and the last term tends to zero. So the $\limsup$ is at most $\varepsilon$ for every $\varepsilon$. The sine and cosine versions are the imaginary and real parts when $g$ is real, and linear combinations of $\lambda$ and $-\lambda$ in general. ∎

> [!theorem] Dirichlet's theorem (pointwise convergence)
> Let $x$ be $2\pi$-periodic and absolutely integrable over a period. Fix $t$ and suppose the one-sided limits $x(t^+)$, $x(t^-)$ exist and that $x$ satisfies one-sided Lipschitz (Hölder-type) bounds at $t$: there are $C,\delta>0$ with
>
> $$
> \lvert x(t+u)-x(t^+)\rvert\le Cu,\qquad\lvert x(t-u)-x(t^-)\rvert\le Cu\qquad(0<u<\delta).
> $$
>
> Then
>
> $$
> \lim_{N\to\infty}S_N(t)=\frac{x(t^+)+x(t^-)}{2}.
> $$
>
> In particular this holds at every point when $x$ is piecewise smooth (piecewise $C^1$ with one-sided limits of $x$ and $x'$ at the break points), which covers the Dirichlet conditions as used in [Ch. 15, §15.6](../part-5-transforms/15-fourier-series.md#156-convergence) for all practical waveforms.

**Proof.** *Step 1 — the partial sum as a convolution.* Inserting the coefficient formula into $S_N$ and exchanging the finite sum with the integral,

$$
S_N(t)=\sum_{n=-N}^{N}\frac1{2\pi}\int_{-\pi}^{\pi}x(\tau)e^{-jn\tau}d\tau\;e^{jnt}=\frac1{2\pi}\int_{-\pi}^{\pi}x(\tau)\,D_N(t-\tau)\,d\tau=\frac1{2\pi}\int_{-\pi}^{\pi}x(t+u)\,D_N(u)\,du,
$$

where we substituted $\tau=t+u$, used the evenness of $D_N$, and used periodicity to keep the interval $[-\pi,\pi]$.

*Step 2 — fold onto $(0,\pi)$.* Splitting at $u=0$ and replacing $u$ by $-u$ on $(-\pi,0)$,

$$
S_N(t)=\frac1{2\pi}\int_0^\pi\big[x(t+u)+x(t-u)\big]D_N(u)\,du .
$$

Since $\frac1{2\pi}\int_0^\pi D_N=\frac12$, subtracting $\frac12[x(t^+)+x(t^-)]$ gives

$$
S_N(t)-\frac{x(t^+)+x(t^-)}{2}=\frac1{2\pi}\int_0^\pi g(u)\,\sin\big((N+\tfrac12)u\big)\,du,\qquad g(u)=\frac{x(t+u)-x(t^+)+x(t-u)-x(t^-)}{\sin(u/2)} .
$$

*Step 3 — $g$ is absolutely integrable on $(0,\pi)$.* On $(0,\delta)$, the numerator is at most $2Cu$ in modulus and $\sin(u/2)\ge u/\pi$ (Jordan's inequality, §D.16), so $\lvert g(u)\rvert\le2\pi C$: bounded. On $[\delta,\pi]$, $\sin(u/2)\ge\sin(\delta/2)>0$, so $\lvert g\rvert\le\big(\lvert x(t+u)\rvert+\lvert x(t-u)\rvert+\lvert x(t^+)\rvert+\lvert x(t^-)\rvert\big)/\sin(\delta/2)$, which is integrable because $x$ is absolutely integrable over a period.

*Step 4 — conclude.* By the Riemann–Lebesgue lemma with $\lambda=N+\tfrac12\to\infty$, the integral tends to zero. ∎

**Remarks.**

1. *Localisation.* The proof shows that convergence at $t$ depends only on the behaviour of $x$ in an arbitrarily small neighbourhood of $t$ (the part of the integral over $[\delta,\pi]$ always tends to zero). Two functions that agree near $t$ have Fourier series that converge or diverge together at $t$, even though their coefficients may be completely different.
2. *Why a regularity hypothesis is needed.* Continuity alone is not enough: there are continuous periodic functions whose Fourier series diverges at a point (Remark 2 of §D.7). The Lipschitz bound can be weakened to Dini's condition $\int_0^\delta\lvert x(t\pm u)-x(t^\pm)\rvert\,du/u<\infty$, or replaced by bounded variation near $t$ (the Dirichlet–Jordan theorem, the version stated in Ch. 15 with "finitely many maxima and minima"); both proofs follow the same three steps.
3. *Midpoint value.* At a jump the series converges to the average because $D_N$ is even: it weights the two sides of the jump equally.
4. *Uniform convergence.* If $x$ is continuous and piecewise smooth, then $\sum\lvert c_n\rvert<\infty$ (since $c_n=c_n'/(jn)$ where $c_n'$ are the coefficients of $x'$, and Cauchy–Schwarz with Bessel for $x'$ gives $\sum\lvert c_n\rvert\le\lvert c_0\rvert+(\sum\lvert c'_n\rvert^2)^{1/2}(2\sum n^{-2})^{1/2}$), so by the M-test the series converges uniformly — necessarily to $x$, by Dirichlet's theorem.

```desmos-graph
left=-3.5; right=3.5; bottom=-5; top=19
---
y=\frac{\sin(8.5x)}{\sin(x/2)}|blue
y=\frac{1}{\sin(x/2)}|grey|dashed|0<x
y=-\frac{1}{\sin(x/2)}|grey|dashed|0<x
```

*The Dirichlet kernel $D_8(u)=\sin(8.5u)/\sin(u/2)$: a central peak of height $2N+1=17$ and width $\approx 2\pi/(N+\frac12)$, with side-lobes confined by the envelope $\pm1/\sin(u/2)$ (dashed). Its area is fixed but its absolute area grows like $\ln N$ — the source of both localisation and Gibbs.*

## D.9 The Gibbs overshoot

> [!definition] Sine integral
> $\displaystyle\operatorname{Si}(x)=\int_0^x\frac{\sin v}{v}\,dv$. Numerically $\operatorname{Si}(\pi)=1.851937\ldots$ and $\operatorname{Si}(\infty)=\pi/2$.

> [!theorem] Gibbs phenomenon for the square wave
> Let $x$ be the $2\pi$-periodic square wave equal to $1$ on $(0,\pi)$ and $-1$ on $(-\pi,0)$, a jump of size $J=2$ at $t=0$. Its Fourier series is $\frac4\pi\sum_{k\ge1}\frac{\sin(2k-1)t}{2k-1}$; let $S_m$ denote the sum of the first $m$ terms. Then:
>
> 1. $\displaystyle S_m(t)=\frac2\pi\int_0^t\frac{\sin 2mu}{\sin u}\,du$;
> 2. on $(0,\pi)$ the maximum of $S_m$ is attained at $t_m=\dfrac{\pi}{2m}$;
> 3. $\displaystyle\lim_{m\to\infty}S_m(t_m)=\frac2\pi\operatorname{Si}(\pi)=1.178980\ldots$
>
> So the overshoot above the value $1$ tends to $\big(\tfrac1\pi\operatorname{Si}(\pi)-\tfrac12\big)J=0.0894899\,J$, about **8.95% of the jump**, while its location $t_m\to0$.

**Proof.** *1.* Differentiating term by term (a finite sum),

$$
S_m'(t)=\frac4\pi\sum_{k=1}^{m}\cos(2k-1)t=\frac{2}{\pi}\,\frac{\sin 2mt}{\sin t}\qquad(t\notin\pi\mathbb{Z}),
$$

because $2\sin t\cos(2k-1)t=\sin 2kt-\sin(2k-2)t$, which telescopes to $\sin2mt$. The quotient extends continuously to $t=0$ (value $4m/\pi$), and $S_m(0)=0$, so integration from $0$ gives the formula.

*2.* On $(0,\pi)$ the critical points are the zeros of $\sin2mt$, $t=k\pi/(2m)$, $k=1,\dots,2m-1$. Since $\sin t>0$ there, $S_m'$ has the sign of $\sin 2mt$: positive on $(0,\pi/2m)$, then alternating. So the local maxima are at $\tau_i=(2i-1)\pi/(2m)$. Between consecutive maxima,

$$
S_m(\tau_{i+1})-S_m(\tau_i)=\frac2\pi\left(\int_{\tau_i}^{i\pi/m}+\int_{i\pi/m}^{\tau_{i+1}}\right)\frac{\sin2mu}{\sin u}\,du .
$$

On the first sub-interval $\sin2mu\le0$, on the second $\sin2mu\ge0$, and the two half-lobes of $\lvert\sin2mu\rvert$ are translates of each other by $\pi/(2m)$. If $\tau_{i+1}\le\pi/2$, the weight $1/\sin u$ is strictly decreasing over both, so the negative half-lobe outweighs the positive one and $S_m(\tau_{i+1})<S_m(\tau_i)$. Thus the maxima decrease on $(0,\pi/2]$. Finally $S_m(\pi-t)=S_m(t)$, because $\sin\big((2k-1)(\pi-t)\big)=\sin(2k-1)t$, so the maxima on $[\pi/2,\pi)$ mirror those on $(0,\pi/2]$. Hence the largest is $S_m(\tau_1)=S_m(\pi/2m)$.

*3.* Substitute $v=2mu$:

$$
S_m\!\left(\frac{\pi}{2m}\right)=\frac2\pi\int_0^{\pi}\frac{\sin v}{2m\sin\big(v/2m\big)}\,dv .
$$

For $0<v\le\pi$, $2m\sin(v/2m)\to v$ as $m\to\infty$, and by Jordan's inequality $2m\sin(v/2m)\ge 2m\cdot\frac{2}{\pi}\cdot\frac{v}{2m}=\frac{2v}{\pi}$, so the integrand is bounded by $\frac{\pi}{2}\,\frac{\sin v}{v}\le\frac\pi2$. Dominated convergence (or, elementarily, uniform convergence of the integrand on $[0,\pi]$) gives

$$
\lim_{m\to\infty}S_m\!\left(\frac{\pi}{2m}\right)=\frac2\pi\int_0^\pi\frac{\sin v}{v}\,dv=\frac2\pi\operatorname{Si}(\pi).
$$

Numerically $\frac2\pi\operatorname{Si}(\pi)=1.1789797$, an overshoot of $0.1789797=0.0894899\times J$ with $J=2$. ∎

**Numerical check.** Direct summation gives $S_{10}(\pi/20)=1.179814$, $S_{100}(\pi/200)=1.178988$, $S_{1000}(\pi/2000)=1.178980$; a dense grid search on $(0,\pi)$ for $m=10$ confirms the global maximum at $t=0.15708=\pi/20$.

**Remarks.**

1. *General jumps.* If $x$ is piecewise smooth with a jump $J=x(t_0^+)-x(t_0^-)$ at $t_0$, then $x-\frac J2\,q(t-t_0)$, with $q$ the square wave above, is continuous and piecewise smooth near $t_0$; its Fourier series converges uniformly there (§D.8, Remark 4), so it contributes no overshoot, and the overshoot of $x$ is $\frac J2\times0.1789797=0.0895\,J$ on each side.
2. *Consistency with convergence theorems.* The peak moves to the jump as $m\to\infty$, so at each fixed $t\neq t_0$ the series converges (Dirichlet), and the overshoot region has width $O(1/m)$, so its contribution to the mean-square error vanishes (Parseval). What fails is *uniform* convergence: $\sup_t\lvert S_m-x\rvert\not\to0$. No series of continuous functions can converge uniformly to a discontinuous limit.
3. *Engineering reading.* Truncating a spectrum with a brick-wall window causes 9% ringing on every edge, whatever the bandwidth — the same constant governs the step response of an ideal low-pass filter, whose peak is $\frac12+\frac1\pi\operatorname{Si}(\pi)=1.0895$ of the step. See [Ch. 15, §15.6](../part-5-transforms/15-fourier-series.md#156-convergence) and [Ch. 16, §16.5](../part-5-transforms/16-fourier-transforms.md#165-convolution-and-multiplication).

```desmos-graph
left=-0.6; right=3.8; bottom=-1.4; top=1.4
---
y=\sum_{k=1}^{10}\frac{4\sin((2k-1)x)}{\pi(2k-1)}|blue
y=1.17898|red|dashed
y=1|grey|dashed
(0.15708,1.17981)|label:first peak
```

*Partial sum $S_{10}$ of the square wave. The first peak at $t=\pi/20$ reaches $1.1798$; the red line is the limiting value $\frac2\pi\operatorname{Si}(\pi)=1.1790$, which is not reduced by taking more terms.*

## D.10 Fourier-transform convolution theorem and Plancherel–Parseval

Write $L^1$ for functions with $\int_{-\infty}^{\infty}\lvert x\rvert dt<\infty$ and $L^2$ for those with finite energy $\int\lvert x\rvert^2dt<\infty$. For $x\in L^1$ the transform $X(\omega)=\int x(t)e^{-j\omega t}dt$ exists for every $\omega$, is bounded by $\lVert x\rVert_1$ and is continuous ([Ch. 16, §16.1](../part-5-transforms/16-fourier-transforms.md#161-continuous-fourier-transform)).

> [!theorem] Convolution theorem (Fourier)
> If $x,h\in L^1$, then $y(t)=(x*h)(t)=\int_{-\infty}^{\infty}x(\tau)h(t-\tau)\,d\tau$ converges absolutely for almost every $t$ (for every $t$ if one of the two is also bounded), $y\in L^1$ with $\lVert y\rVert_1\le\lVert x\rVert_1\lVert h\rVert_1$, and
>
> $$
> Y(\omega)=X(\omega)\,H(\omega)\qquad\text{for every }\omega .
> $$

**Proof.** By Tonelli's theorem and the substitution $u=t-\tau$,

$$
\int\!\!\int\lvert x(\tau)\rvert\,\lvert h(t-\tau)\rvert\,d\tau\,dt=\int\lvert x(\tau)\rvert\left(\int\lvert h(u)\rvert du\right)d\tau=\lVert x\rVert_1\lVert h\rVert_1<\infty .
$$

Hence the inner integral $\int\lvert x(\tau)h(t-\tau)\rvert d\tau$ is finite for almost every $t$, which defines $y(t)$, and $\int\lvert y\rvert dt\le\lVert x\rVert_1\lVert h\rVert_1$. (If $\lvert h\rvert\le B$, then $\int\lvert x(\tau)h(t-\tau)\rvert d\tau\le B\lVert x\rVert_1$ for every $t$.) Since $\lvert e^{-j\omega t}\rvert=1$, the function $(\tau,t)\mapsto x(\tau)h(t-\tau)e^{-j\omega t}$ is absolutely integrable on $\mathbb{R}^2$, so Fubini's theorem allows

$$
Y(\omega)=\int\left(\int x(\tau)h(t-\tau)\,d\tau\right)e^{-j\omega t}dt=\int x(\tau)e^{-j\omega\tau}\left(\int h(t-\tau)e^{-j\omega(t-\tau)}dt\right)d\tau=H(\omega)X(\omega).
$$

∎

> [!theorem] Plancherel–Parseval
> If $x\in L^1\cap L^2$, then $X\in L^2$ and
>
> $$
> \int_{-\infty}^{\infty}\lvert x(t)\rvert^2\,dt=\frac1{2\pi}\int_{-\infty}^{\infty}\lvert X(\omega)\rvert^2\,d\omega .
> $$
>
> For $x,y\in L^1\cap L^2$, $\displaystyle\int x(t)\,y^*(t)\,dt=\frac1{2\pi}\int X(\omega)\,Y^*(\omega)\,d\omega$. The map $x\mapsto X/\sqrt{2\pi}$ therefore extends uniquely to a norm-preserving map of all of $L^2$ onto $L^2$ (Plancherel's theorem).

**Proof.** *Step 1 — the autocorrelation.* Let $\tilde x(t)=x^*(-t)$ and $g=x*\tilde x$, i.e.

$$
g(t)=\int x(\tau)\,x^*(\tau-t)\,d\tau .
$$

Because $x\in L^2$, this integral converges absolutely for **every** $t$ by Cauchy–Schwarz, with $\lvert g(t)\rvert\le\lVert x\rVert_2^2=g(0)$. Moreover $g$ is continuous: $\lvert g(t)-g(t')\rvert\le\lVert x\rVert_2\,\big(\int\lvert x(\tau-t)-x(\tau-t')\rvert^2d\tau\big)^{1/2}$, and translation is continuous in $L^2$ (true for continuous compactly supported functions by uniform continuity, and for general $x$ by density). Since $x\in L^1$, the convolution theorem gives $g\in L^1$ and, because $\tilde X(\omega)=\int x^*(-t)e^{-j\omega t}dt=X^*(\omega)$,

$$
G(\omega)=X(\omega)X^*(\omega)=\lvert X(\omega)\rvert^2\ \ge0 .
$$

So the claim is $g(0)=\frac1{2\pi}\int G(\omega)\,d\omega$: an inversion formula at the single point $t=0$ for a function whose transform is non-negative.

*Step 2 — Gaussian regularisation.* For $\varepsilon>0$ the Gaussian $e^{-\varepsilon\omega^2}$ is the transform of

$$
\phi_\varepsilon(t)=\frac1{2\pi}\int e^{-\varepsilon\omega^2}e^{j\omega t}d\omega=\frac{1}{\sqrt{4\pi\varepsilon}}\,e^{-t^2/(4\varepsilon)},
$$

the standard Gaussian integral of [Ch. 16](../part-5-transforms/16-fourier-transforms.md) (complete the square, or use $\mathcal{F}\{e^{-at^2}\}=\sqrt{\pi/a}\,e^{-\omega^2/4a}$). Note that $\phi_\varepsilon>0$, $\int\phi_\varepsilon=1$, and $\phi_\varepsilon$ concentrates at $0$ as $\varepsilon\to0$. Because $g\in L^1$ and $e^{-\varepsilon\omega^2}\in L^1$, the double integral below is absolutely convergent and Fubini gives

$$
\frac1{2\pi}\int G(\omega)e^{-\varepsilon\omega^2}d\omega=\frac1{2\pi}\int\!\!\int g(t)e^{-j\omega t}e^{-\varepsilon\omega^2}\,dt\,d\omega=\int g(t)\,\phi_\varepsilon(-t)\,dt=\int g(t)\,\phi_\varepsilon(t)\,dt .
$$

*Step 3 — the limit $\varepsilon\to0$.* On the right: $\int g\phi_\varepsilon-g(0)=\int[g(t)-g(0)]\phi_\varepsilon(t)dt$. Given $\eta>0$ choose $\delta$ with $\lvert g(t)-g(0)\rvert<\eta$ for $\lvert t\rvert<\delta$; the rest is bounded by $2g(0)\int_{\lvert t\rvert\ge\delta}\phi_\varepsilon\to0$. So the right side tends to $g(0)=\lVert x\rVert_2^2$. On the left: $G\,e^{-\varepsilon\omega^2}$ increases to $G$ as $\varepsilon\downarrow0$ because $G\ge0$, so by the monotone convergence theorem the left side tends to $\frac1{2\pi}\int G$ — which is therefore finite. Equating the limits gives $\lVert x\rVert_2^2=\frac1{2\pi}\int\lvert X\rvert^2d\omega$, and in particular $X\in L^2$.

*Step 4 — inner products and extension.* Polarization, as in §D.7 Step 4, gives the two-signal form. Since $L^1\cap L^2$ is dense in $L^2$ (truncate to $[-T,T]$), the isometry $x\mapsto X/\sqrt{2\pi}$ extends by continuity to all of $L^2$; the extended transform is the mean-square limit of the truncated transforms $X_T(\omega)=\int_{-T}^{T}x(t)e^{-j\omega t}dt$, which is how [Ch. 16, §16.1](../part-5-transforms/16-fourier-transforms.md#161-continuous-fourier-transform) defines it. That the extension is *onto* follows from applying the same argument to the inverse transform. ∎

**Remarks.**

1. *Why non-negativity is the trick.* The inversion formula can fail pointwise for general $L^1$ functions (the transform need not be integrable). Working with $G=\lvert X\rvert^2\ge0$ lets monotone convergence do the work that would otherwise require knowing $G\in L^1$ in advance.
2. *Check.* $x(t)=e^{-t}u(t)$ has energy $\int_0^\infty e^{-2t}dt=\tfrac12$ and $X(\omega)=1/(1+j\omega)$; numerically $\frac{1}{2\pi}\int\frac{d\omega}{1+\omega^2}=0.5000$.
3. *Convention.* The factor $1/(2\pi)$ belongs to the angular-frequency convention used in this book; in hertz, $\int\lvert x\rvert^2dt=\int\lvert X_f(f)\rvert^2df$.
4. Used for energy spectral density in [Ch. 16, §16.7](../part-5-transforms/16-fourier-transforms.md#167-energy-spectral-density) and for LTI filtering in [§16.5](../part-5-transforms/16-fourier-transforms.md#165-convolution-and-multiplication).

## D.11 The Nyquist–Shannon sampling theorem

Notation follows [Ch. 16, §16.8](../part-5-transforms/16-fourier-transforms.md#168-sampling-related-mathematical-foundations): sampling interval $T$, $\omega_s=2\pi/T$, and $\operatorname{sinc}(u)=\sin(\pi u)/(\pi u)$.

> [!theorem] Poisson summation (band-limited form)
> Let $X\in L^1$ vanish outside a bounded interval and be square-integrable, and let $x(t)=\frac1{2\pi}\int X(\omega)e^{j\omega t}d\omega$ (a continuous function). Then the $\omega_s$-periodic function $P(\omega)=\sum_{k=-\infty}^{\infty}X(\omega-k\omega_s)$ — a finite sum at each $\omega$ — has the mean-square convergent Fourier series
>
> $$
> \sum_{k=-\infty}^{\infty}X(\omega-k\omega_s)=T\sum_{n=-\infty}^{\infty}x(nT)\,e^{-j\omega nT}.
> $$

> [!theorem] Sampling theorem (Nyquist–Shannon)
> Let $x$ be as above with $X(\omega)=0$ for $\lvert\omega\rvert>\omega_M$, and let the sampling rate satisfy $\omega_s>2\omega_M$. Then $x$ is uniquely determined by its samples $x(nT)$, the samples are square-summable, and for every real $t$
>
> $$
> x(t)=\sum_{n=-\infty}^{\infty}x(nT)\,\operatorname{sinc}\!\left(\frac{t-nT}{T}\right),
> $$
>
> the symmetric partial sums converging uniformly in $t$ on all of $\mathbb{R}$.

**Proof of Poisson summation.** $P$ is $\omega_s$-periodic. On one period $I=[-\omega_s/2,\omega_s/2)$ only finitely many shifts $X(\omega-k\omega_s)$ are non-zero (the support of $X$ is bounded), so $P$ is square-integrable over $I$. Its Fourier series in the variable $\omega$, with period $\omega_s$ and basis $e^{-j\omega nT}$ (note $\omega_s T=2\pi$), has coefficients

$$
d_n=\frac1{\omega_s}\int_{I}P(\omega)e^{j\omega nT}d\omega=\frac1{\omega_s}\sum_k\int_I X(\omega-k\omega_s)e^{j\omega nT}d\omega=\frac1{\omega_s}\int_{-\infty}^{\infty}X(\nu)e^{j\nu nT}d\nu=\frac{2\pi}{\omega_s}\,x(nT)=T\,x(nT).
$$

In the third step we substituted $\nu=\omega-k\omega_s$, used $e^{jk\omega_snT}=e^{j2\pi kn}=1$, and observed that the translates $I+k\omega_s$ tile the line. By Parseval's theorem (§D.7, which applies to any square-integrable periodic function), the series $\sum_nd_ne^{-j\omega nT}$ converges to $P$ in mean square over a period, and $\frac{1}{\omega_s}\int_I\lvert P\rvert^2=\sum_n\lvert d_n\rvert^2$. ∎

**Proof of the sampling theorem.** *Step 1 — no aliasing.* Since $\omega_M<\omega_s/2$, for $\omega\in I=(-\omega_s/2,\omega_s/2)$ every shift with $k\neq0$ is evaluated at $\lvert\omega-k\omega_s\rvert>\omega_s/2>\omega_M$, where $X$ vanishes. So $P=X$ on $I$, and Poisson summation gives, in mean square on $I$,

$$
X(\omega)=T\sum_{n=-\infty}^{\infty}x(nT)\,e^{-j\omega nT}\qquad(\lvert\omega\rvert<\omega_s/2).
$$

This already shows that the samples determine $X$ on $I$ (and $X=0$ outside), hence determine $x$; and Parseval gives $\sum_n\lvert x(nT)\rvert^2=\frac{1}{T^2\omega_s}\int_I\lvert X\rvert^2<\infty$.

*Step 2 — inversion term by term.* Since $X=0$ outside $I$,

$$
x(t)=\frac1{2\pi}\int_I X(\omega)\,e^{j\omega t}\,d\omega=\frac{\omega_s}{2\pi}\,\big\langle X,\ e_t\big\rangle_I,\qquad e_t(\omega)=e^{-j\omega t},
$$

with $\langle F,G\rangle_I=\frac{1}{\omega_s}\int_IF\,G^*$. An inner product with a fixed square-integrable function is continuous with respect to mean-square convergence (Cauchy–Schwarz), so we may insert the series of Step 1 and integrate term by term:

$$
x(t)=\sum_nx(nT)\,\frac{T}{2\pi}\int_{-\omega_s/2}^{\omega_s/2}e^{j\omega(t-nT)}d\omega=\sum_nx(nT)\,\frac{T}{2\pi}\cdot\frac{2\sin\big(\omega_s(t-nT)/2\big)}{t-nT}.
$$

With $\omega_s/2=\pi/T$, the coefficient is $\dfrac{\sin\big(\pi(t-nT)/T\big)}{\pi(t-nT)/T}=\operatorname{sinc}\!\big((t-nT)/T\big)$ (value $1$ at $t=nT$).

*Step 3 — uniform convergence.* Let $x_N$ be the partial sum over $\lvert n\rvert\le N$. By the same computation, $x(t)-x_N(t)=\frac{\omega_s}{2\pi}\langle X-X_N,e_t\rangle_I$, where $X_N$ is the corresponding partial Fourier sum. Since $\lVert e_t\rVert_I=1$ for every $t$, Cauchy–Schwarz gives $\lvert x(t)-x_N(t)\rvert\le\frac{\omega_s}{2\pi}\lVert X-X_N\rVert_I\to0$, uniformly in $t$. ∎

**Numerical check.** $x(t)=\operatorname{sinc}^2(t/2)$ has a triangular spectrum supported on $\lvert\omega\rvert\le\pi$. With $T=0.8$ ($\omega_s=7.85>2\pi$) and $\lvert n\rvert\le 20000$, the reconstruction at $t=0.37$ gives $0.8923548471$, equal to $x(0.37)$ to ten digits.

**Remarks.**

1. *Why $\omega_s>2\omega_M$.* If $\omega_s<2\omega_M$, the shifted copies in $P$ overlap on $I$, so $P\neq X$ there: the samples determine only the folded spectrum $P$, and different signals share the same samples (aliasing). For example $\cos(\omega_1t)$ and $\cos\big((\omega_s-\omega_1)t\big)$ agree at every $t=nT$.
2. *The boundary case.* With $\omega_s=2\omega_M$ the proof goes through whenever $X\in L^2$ is supported on $[-\omega_M,\omega_M]$, because the endpoints have measure zero. But it fails for spectra containing impulses at $\pm\omega_M$: $x(t)=\sin(\omega_Mt)$ sampled at $T=\pi/\omega_M$ gives all-zero samples. This is why the theorem is stated with strict inequality.
3. *Interpretation.* The sinc functions $\operatorname{sinc}\big((t-nT)/T\big)$ are orthogonal and are the impulse responses of the ideal low-pass filter $H_r(\omega)=T\operatorname{rect}(\omega/\omega_s)$ shifted to the sampling instants. Reconstruction is "filter the impulse-sampled signal", exactly as derived with impulses in [Ch. 16, §16.8](../part-5-transforms/16-fourier-transforms.md#168-sampling-related-mathematical-foundations); the proof here avoids generalized functions by working with Fourier series in $\omega$.
4. *Practicality.* The sinc decays only like $1/\lvert t\rvert$, so truncating the series converges slowly, and real converters oversample and use smoother interpolation kernels.

```desmos-graph
left=-3.2; right=4.2; bottom=-0.5; top=1.4
---
y=\frac{\sin(\pi x)}{\pi x}|blue|dashed
y=0.8\frac{\sin(\pi(x-1))}{\pi(x-1)}|green|dashed
y=0.3\frac{\sin(\pi(x-2))}{\pi(x-2)}|orange|dashed
y=\frac{\sin(\pi x)}{\pi x}+0.8\frac{\sin(\pi(x-1))}{\pi(x-1)}+0.3\frac{\sin(\pi(x-2))}{\pi(x-2)}|black
(0,1)|label:x(0)
(1,0.8)|label:x(T)
(2,0.3)|label:x(2T)
```

*Sinc reconstruction with $T=1$ and three non-zero samples: each dashed sinc equals its sample at its own instant and vanishes at all the others, so the black sum passes exactly through every sample.*

## D.12 Z-transform: the ROC is an annulus; BIBO stability

> [!definition] Z-transform and ROC
> For a sequence $x[n]$, $n\in\mathbb{Z}$, the Z-transform is $X(z)=\sum_{n=-\infty}^{\infty}x[n]z^{-n}$. Its **region of convergence** $R_x$ is the set of $z\in\mathbb{C}\setminus\{0\}$ (with $z=\infty$ allowed when $x[n]=0$ for $n<0$, and $z=0$ allowed when $x[n]=0$ for $n>0$) at which the series converges **absolutely**. A discrete-time LTI system with impulse response $h[n]$ is **BIBO stable** if every bounded input $\lvert x[n]\rvert\le B_x$ produces a bounded output $y=h*x$.

> [!theorem] Shape of the ROC
> Let
>
> $$
> r_1=\limsup_{n\to\infty}\lvert x[n]\rvert^{1/n},\qquad\frac{1}{r_2}=\limsup_{m\to\infty}\lvert x[-m]\rvert^{1/m}
> $$
>
> (with $1/0=\infty$, $1/\infty=0$). Then the series converges absolutely for $r_1<\lvert z\rvert<r_2$ and diverges for $\lvert z\rvert<r_1$ and for $\lvert z\rvert>r_2$. Consequently $R_x$ is an annulus $\{r_1<\lvert z\rvert<r_2\}$, together possibly with some or all points of the boundary circles, and is empty if $r_1>r_2$. On the open annulus $X$ is analytic.

> [!theorem] Stability and the unit circle
> An LTI system is BIBO stable if and only if $\sum_n\lvert h[n]\rvert<\infty$, and this holds if and only if the unit circle $\lvert z\rvert=1$ lies in $R_h$. A causal system with rational $H(z)$ is therefore stable if and only if all its poles satisfy $\lvert p\rvert<1$.

**Proof of the ROC theorem.** Put $r=\lvert z\rvert$ and split $X(z)=\sum_{n\ge0}x[n]z^{-n}+\sum_{m\ge1}x[-m]z^{m}$. Absolute convergence of the whole series is equivalent to absolute convergence of both parts, and depends only on $r$ because $\lvert x[n]z^{-n}\rvert=\lvert x[n]\rvert r^{-n}$.

*Causal part.* Apply the root test to $\sum_{n\ge0}\lvert x[n]\rvert r^{-n}$: $\limsup_n\big(\lvert x[n]\rvert r^{-n}\big)^{1/n}=r_1/r$. If $r>r_1$ this is $<1$ and the series converges (choose $q$ with $r_1/r<q<1$; then $\lvert x[n]\rvert r^{-n}\le q^n$ for large $n$, and compare with the geometric series). If $r<r_1$ it is $>1$, so $\lvert x[n]\rvert r^{-n}\ge1$ for infinitely many $n$, the terms do not tend to zero, and the series diverges.

*Anticausal part.* The same test on $\sum_{m\ge1}\lvert x[-m]\rvert r^{m}$ gives $\limsup\big(\lvert x[-m]\rvert r^m\big)^{1/m}=r/r_2$: convergence for $r<r_2$, divergence for $r>r_2$.

Both parts converge for $r_1<r<r_2$ and at least one diverges outside $[r_1,r_2]$; on the circles $r=r_1$, $r=r_2$ the test is inconclusive, which is why boundary points may or may not belong to $R_x$. Analyticity: on any closed sub-annulus $r_1<\rho_1\le\lvert z\rvert\le\rho_2<r_2$ the terms are bounded by $\lvert x[n]\rvert\rho_1^{-n}$ ($n\ge0$) and $\lvert x[-m]\rvert\rho_2^{m}$ ($m\ge1$), which form convergent series by the above, so the Laurent series converges uniformly there (M-test) and its sum is analytic, each term being analytic on $\mathbb{C}\setminus\{0\}$. ∎

**Proof of the stability theorem.** *Sufficiency.* If $S=\sum\lvert h[k]\rvert<\infty$ and $\lvert x[n]\rvert\le B_x$, then

$$
\lvert y[n]\rvert=\Big\lvert\sum_kh[k]x[n-k]\Big\rvert\le B_x\sum_k\lvert h[k]\rvert=B_xS ,
$$

and the series defining $y[n]$ converges absolutely for each $n$.

*Necessity.* Suppose $\sum\lvert h[k]\rvert=\infty$. Define the bounded input $x[n]=h^*[-n]/\lvert h[-n]\rvert$ where $h[-n]\neq0$ and $x[n]=0$ otherwise, so $\lvert x[n]\rvert\le1$. Then

$$
y[0]=\sum_kh[k]x[-k]=\sum_k\frac{h[k]h^*[k]}{\lvert h[k]\rvert}=\sum_k\lvert h[k]\rvert=\infty ,
$$

so this bounded input produces an output that is not even finite at $n=0$, and the system is not BIBO stable. (To avoid a divergent defining sum, use instead the truncated inputs $x_N$, equal to $x$ for $\lvert n\rvert\le N$ and $0$ otherwise. Each satisfies $\lvert x_N\rvert\le1$ and has a well-defined output with $y_N[0]=\sum_{\lvert k\rvert\le N}\lvert h[k]\rvert$. If the system were stable, the map $x\mapsto y$, being linear and time-invariant, would satisfy a bound $\lvert y[n]\rvert\le C\,\sup\lvert x\rvert$ — for a linear system, boundedness of outputs over the unit ball of inputs is equivalent to such a bound by the uniform boundedness principle — whereas $y_N[0]\to\infty$.)

*The unit circle.* On $\lvert z\rvert=1$, $\sum_n\lvert h[n]z^{-n}\rvert=\sum_n\lvert h[n]\rvert$, so the unit circle lies in $R_h$ (with absolute convergence, as in the definition) if and only if $h$ is absolutely summable.

*Causal rational case.* If $h$ is causal, $h[-m]=0$ for $m\ge1$, so $r_2=\infty$ and $R_h\supseteq\{\lvert z\rvert>r_1\}$. For rational $H$, partial fractions give $h[n]$ as a finite sum of terms $n^kp^nu[n]$ (plus finitely many impulses), and $r_1=\max\lvert p\rvert$ over the poles; $X$ cannot be analytic at a pole, so no pole lies in the open ROC. Thus $R_h=\{\lvert z\rvert>\max\lvert p\rvert\}$ (including $\infty$), and it contains the unit circle iff $\max\lvert p\rvert<1$. A pole exactly on the circle gives a term $n^kp^n$ with $\lvert p\rvert=1$, which is not absolutely summable, so such a system is unstable even though the circle is the ROC's boundary. ∎

**Remarks.**

1. *The role of absolute convergence.* With a definition of the ROC by mere convergence, boundary points could be included for non-absolutely-summable sequences (e.g. $h[n]=u[n-1]/n$ converges on $\lvert z\rvert=1$ except at $z=1$), and the clean equivalence with stability would fail. Within the open annulus the two definitions agree.
2. *Analogy.* The continuous-time statement is that a causal LTI system is BIBO stable iff $\int\lvert h\rvert<\infty$ iff the $j\omega$-axis lies in the ROC of $H(s)$; the map $z=e^{sT}$ carries the axis onto the unit circle and the left half-plane into the disc ([Ch. 17, §17.7](../part-5-transforms/17-z-transforms.md#177-discrete-time-system-analysis)).
3. Used in [Ch. 17, §17.3](../part-5-transforms/17-z-transforms.md#173-region-of-convergence) and [§17.6](../part-5-transforms/17-z-transforms.md#176-system-function-representation).

## D.13 Cauchy–Riemann equations: necessity and sufficiency

> [!definition] Complex differentiability
> Let $f=u+jv$ be defined on an open set $D\subseteq\mathbb{C}$, with $z=x+jy$ and $u,v$ real. $f$ is **complex differentiable** at $z_0\in D$ if
>
> $$
> f'(z_0)=\lim_{\Delta z\to0}\frac{f(z_0+\Delta z)-f(z_0)}{\Delta z}
> $$
>
> exists, the limit being taken as $\Delta z\to0$ through **all** complex values. $f$ is **analytic** on $D$ if it is complex differentiable at every point of $D$.

> [!theorem] Cauchy–Riemann equations
> **(Necessity.)** If $f$ is complex differentiable at $z_0=x_0+jy_0$, then the four partial derivatives $u_x,u_y,v_x,v_y$ exist at $(x_0,y_0)$ and satisfy
>
> $$
> u_x=v_y,\qquad u_y=-v_x,\qquad\text{and}\qquad f'(z_0)=u_x+jv_x=v_y-ju_y .
> $$
>
> **(Sufficiency.)** Conversely, if $u$ and $v$ are (real-)differentiable at $(x_0,y_0)$ — in particular, if their first partials exist near $(x_0,y_0)$ and are continuous at it — and satisfy the Cauchy–Riemann equations there, then $f$ is complex differentiable at $z_0$.

**Proof of necessity.** Let $L=f'(z_0)$. Since the limit exists along every path, take $\Delta z=h$ real:

$$
L=\lim_{h\to0}\frac{u(x_0+h,y_0)-u(x_0,y_0)}{h}+j\,\frac{v(x_0+h,y_0)-v(x_0,y_0)}{h}.
$$

A complex limit exists iff the limits of its real and imaginary parts exist, so $u_x,v_x$ exist and $L=u_x+jv_x$. Now take $\Delta z=jh$, $h$ real:

$$
L=\lim_{h\to0}\frac{u(x_0,y_0+h)-u(x_0,y_0)+j\big[v(x_0,y_0+h)-v(x_0,y_0)\big]}{jh}=\frac{1}{j}\big(u_y+jv_y\big)=v_y-ju_y .
$$

Equating real and imaginary parts of the two expressions for $L$ gives $u_x=v_y$ and $v_x=-u_y$. ∎

**Proof of sufficiency.** Real differentiability of $u$ at $(x_0,y_0)$ means that, with $\Delta z=h+jk$,

$$
u(x_0+h,y_0+k)-u(x_0,y_0)=u_xh+u_yk+\varepsilon_1(h,k)\,\lvert\Delta z\rvert,\qquad\varepsilon_1\to0\text{ as }\Delta z\to0,
$$

and similarly for $v$ with a remainder $\varepsilon_2$. (If the partials are continuous near the point, this follows from the mean-value theorem applied separately in $x$ and in $y$.) Using $v_x=-u_y$ and $v_y=u_x$,

$$
\begin{aligned}
f(z_0+\Delta z)-f(z_0)&=(u_xh+u_yk)+j(v_xh+v_yk)+(\varepsilon_1+j\varepsilon_2)\lvert\Delta z\rvert\\
&=(u_xh+u_yk)+j(-u_yh+u_xk)+(\varepsilon_1+j\varepsilon_2)\lvert\Delta z\rvert\\
&=(u_x-ju_y)(h+jk)+(\varepsilon_1+j\varepsilon_2)\lvert\Delta z\rvert .
\end{aligned}
$$

(Check the last line: $(u_x-ju_y)(h+jk)=u_xh+u_yk+j(u_xk-u_yh)$.) Dividing by $\Delta z$,

$$
\frac{f(z_0+\Delta z)-f(z_0)}{\Delta z}=(u_x-ju_y)+(\varepsilon_1+j\varepsilon_2)\frac{\lvert\Delta z\rvert}{\Delta z},
$$

and the last term has modulus $\lvert\varepsilon_1+j\varepsilon_2\rvert\to0$. So $f'(z_0)=u_x-ju_y=u_x+jv_x$. ∎

**Remarks.**

1. *Geometric meaning.* Real differentiability says the map $(x,y)\mapsto(u,v)$ is approximated by a $2\times2$ Jacobian $\begin{pmatrix}u_x&u_y\\v_x&v_y\end{pmatrix}$. The Cauchy–Riemann equations say this matrix has the form $\begin{pmatrix}a&-b\\b&a\end{pmatrix}$, i.e. it is multiplication by the complex number $a+jb$: a rotation combined with a scaling. That is why analytic maps are conformal where $f'\neq0$.
2. *CR alone is not sufficient.* Let $f(z)=\sqrt{\lvert xy\rvert}$ (so $u=\sqrt{\lvert xy\rvert}$, $v=0$). On the axes $u=0$, so $u_x(0,0)=u_y(0,0)=0=v_x=v_y$ and the Cauchy–Riemann equations hold at the origin. But along $z=t(1+j)$, $t>0$, $\big(f(z)-f(0)\big)/z=t/\big(t(1+j)\big)=1/(1+j)$, whereas along the real axis the quotient is $0$; the limit does not exist. Here $u$ is not real-differentiable at $0$, which is the hypothesis that fails.
3. *Differentiability alone is not sufficient.* $f(z)=\bar z$ has $u=x$, $v=-y$, smooth, but $u_x=1\neq v_y=-1$; it is nowhere complex differentiable (the difference quotient is $\overline{\Delta z}/\Delta z$, which takes every value on the unit circle).
4. Used in [Ch. 18, §18.3](../part-5-transforms/18-complex-analysis.md#183-cauchyriemann-equations), and, through $\nabla^2u=0$, for the harmonic potentials of electrostatics.

## D.14 Cauchy's integral theorem and integral formula

> [!definition] Contours and positive orientation
> A **contour** is a continuous, piecewise continuously differentiable curve $z(\tau)$, $a\le\tau\le b$, and $\int_Cf(z)\,dz=\int_a^bf(z(\tau))\,z'(\tau)\,d\tau$. A **simple closed** contour bounds a bounded region $\Omega$ (Jordan curve theorem); it is **positively oriented** when $\Omega$ lies to its left. A domain $D$ is **simply connected** if the interior of every simple closed contour in $D$ lies in $D$.

> [!theorem] Cauchy's integral theorem (C¹ version)
> Let $f$ be analytic on a simply connected domain $D$, with $f'$ continuous on $D$. Then for every simple closed contour $C$ in $D$
>
> $$
> \oint_Cf(z)\,dz=0 .
> $$
>
> More generally, if $\Omega$ is a bounded region whose boundary $\partial\Omega$ consists of finitely many disjoint simple closed contours, oriented with $\Omega$ to the left, and $f$ is analytic with continuous $f'$ on an open set containing $\overline\Omega=\Omega\cup\partial\Omega$, then $\oint_{\partial\Omega}f\,dz=0$.

> [!theorem] Cauchy's integral formula
> Under the hypotheses of the first statement, for every $z_0$ inside $C$
>
> $$
> f(z_0)=\frac1{2\pi j}\oint_C\frac{f(z)}{z-z_0}\,dz,
> \qquad\text{and}\qquad
> f^{(n)}(z_0)=\frac{n!}{2\pi j}\oint_C\frac{f(z)}{(z-z_0)^{n+1}}\,dz\quad(n\ge1).
> $$

**Proof of Cauchy's theorem.** Write $f=u+jv$ and $dz=dx+j\,dy$. Then

$$
\oint_Cf\,dz=\oint_C(u\,dx-v\,dy)+j\oint_C(v\,dx+u\,dy).
$$

Since $f'=u_x+jv_x=v_y-ju_y$ is continuous, the partials of $u$ and $v$ are continuous on $D$, and $\overline\Omega$ (the region bounded by $C$, together with $C$) lies in $D$ by simple connectivity. Green's theorem $\oint_{\partial\Omega}(P\,dx+Q\,dy)=\iint_\Omega(Q_x-P_y)\,dx\,dy$ ([Ch. 7](../part-2-calculus/07-vector-calculus.md)), applied to each real line integral, gives

$$
\oint_Cf\,dz=\iint_\Omega(-v_x-u_y)\,dx\,dy+j\iint_\Omega(u_x-v_y)\,dx\,dy=0,
$$

since both integrands vanish identically by the Cauchy–Riemann equations (§D.13). Green's theorem holds verbatim for regions with finitely many holes (cut the region into simply connected pieces; the integrals along the cuts cancel in pairs), which proves the general statement. ∎

**Proof of the integral formula.** Fix $z_0$ inside $C$ and choose $r>0$ so small that the closed disc $\lvert z-z_0\rvert\le r$ lies inside $C$. Let $C_r$ be the circle $\lvert z-z_0\rvert=r$, positively oriented, and $\Omega_r$ the region between $C$ and $C_r$. The function $g(z)=f(z)/(z-z_0)$ is analytic with continuous derivative on a neighbourhood of $\overline{\Omega_r}$, whose boundary is $C$ (positive) together with $C_r$ (negative). By the general form of Cauchy's theorem,

$$
\oint_C\frac{f(z)}{z-z_0}\,dz=\oint_{C_r}\frac{f(z)}{z-z_0}\,dz\qquad\text{for every small }r>0 .
$$

Now split the right side:

$$
\oint_{C_r}\frac{f(z)}{z-z_0}\,dz=f(z_0)\oint_{C_r}\frac{dz}{z-z_0}+\oint_{C_r}\frac{f(z)-f(z_0)}{z-z_0}\,dz .
$$

With $z=z_0+re^{j\theta}$, $dz=jre^{j\theta}d\theta$, the first integral is $\int_0^{2\pi}j\,d\theta=2\pi j$. The second is bounded, by the ML inequality ([Ch. 18, §18.4](../part-5-transforms/18-complex-analysis.md#184-contour-integration)), by $2\pi r\cdot\max_{C_r}\lvert f(z)-f(z_0)\rvert/r=2\pi\max_{C_r}\lvert f(z)-f(z_0)\rvert$, which tends to $0$ as $r\to0$ by continuity of $f$. Since the left side does not depend on $r$, it equals $2\pi j\,f(z_0)$.

*Derivatives.* Let $F(w)=\frac1{2\pi j}\oint_C\frac{f(z)}{z-w}dz$ for $w$ inside $C$; we have shown $F=f$ there. For $w$ at distance $\ge d>0$ from $C$, the difference quotient of the integrand,

$$
\frac1h\left(\frac1{z-w-h}-\frac1{z-w}\right)=\frac{1}{(z-w-h)(z-w)}\ \longrightarrow\ \frac1{(z-w)^2}
$$

uniformly for $z\in C$ as $h\to0$ (for $\lvert h\rvert<d/2$ the denominators are bounded below by $d^2/2$, and the difference from the limit is $h/[(z-w-h)(z-w)^2]$, bounded by $2\lvert h\rvert/d^3$). Hence one may differentiate under the integral sign, and by induction $f^{(n)}(w)=\frac{n!}{2\pi j}\oint_C\frac{f(z)}{(z-w)^{n+1}}dz$. ∎

**Remarks.**

1. *Removing the C¹ hypothesis.* The Cauchy–Goursat theorem proves $\oint_Cf\,dz=0$ assuming only that $f'$ exists at every point, by subdividing triangles and estimating — no continuity of $f'$ and no Green's theorem. Combined with the derivative formula just proved, it shows that **an analytic function is automatically infinitely differentiable**, so the C¹ hypothesis used here costs no generality *a posteriori*. See [Ch. 18, §18.5](../part-5-transforms/18-complex-analysis.md#185-cauchys-integral-theorem).
2. *Simple connectivity is essential.* $f(z)=1/z$ is analytic on $\mathbb{C}\setminus\{0\}$, but $\oint_{\lvert z\rvert=1}dz/z=2\pi j\neq0$. The punctured plane is not simply connected — the hole at $0$ is exactly what the integral detects, and the residue theorem turns this into a computational tool.
3. *Mean-value property.* With $C$ a circle of radius $r$ about $z_0$, the formula reads $f(z_0)=\frac1{2\pi}\int_0^{2\pi}f(z_0+re^{j\theta})\,d\theta$: an analytic function (and hence its real part, a harmonic function) equals its average on every circle — the basis of the maximum-modulus principle and of Liouville's theorem in Ch. 18.
4. *Consequence for the Laplace transform.* The Bromwich inversion integral of [Ch. 14, §14.4](../part-5-transforms/14-laplace-transforms.md#144-inverse-laplace-transforms) is a contour integral of $X(s)e^{st}$; Cauchy's theorem is what allows the vertical line to be moved freely within the ROC.

## D.15 The residue theorem

> [!definition] Residue
> If $f$ is analytic in a punctured disc $0<\lvert z-z_k\rvert<\rho$, it has a Laurent expansion $f(z)=\sum_{n=-\infty}^{\infty}a_n(z-z_k)^n$ there, converging uniformly on every circle $\lvert z-z_k\rvert=r$, $0<r<\rho$ ([Ch. 18, §18.6](../part-5-transforms/18-complex-analysis.md#186-taylor-and-laurent-series)). The **residue** of $f$ at $z_k$ is $\operatorname{Res}_{z=z_k}f=a_{-1}$.

> [!theorem] Residue theorem
> Let $C$ be a positively oriented simple closed contour, and let $f$ be analytic on an open set containing $C$ and its interior, except at finitely many points $z_1,\dots,z_K$ inside $C$. Then
>
> $$
> \oint_Cf(z)\,dz=2\pi j\sum_{k=1}^{K}\operatorname{Res}_{z=z_k}f .
> $$

**Proof.** *Step 1 — one circle.* Let $C_k$ be the circle $\lvert z-z_k\rvert=r$, positively oriented, with $r$ small. On $C_k$ the Laurent series converges uniformly, so it may be integrated term by term. For each integer $n$, with $z=z_k+re^{j\theta}$,

$$
\oint_{C_k}(z-z_k)^n\,dz=\int_0^{2\pi}r^ne^{jn\theta}\cdot jre^{j\theta}\,d\theta=jr^{n+1}\int_0^{2\pi}e^{j(n+1)\theta}\,d\theta=\begin{cases}2\pi j,&n=-1,\\0,&n\neq-1.\end{cases}
$$

Hence $\oint_{C_k}f\,dz=2\pi j\,a_{-1}=2\pi j\operatorname{Res}_{z=z_k}f$. Only the $(z-z_k)^{-1}$ term survives, because every other power has a single-valued antiderivative on the punctured disc.

*Step 2 — deformation.* Choose $r$ so small that the closed discs $\lvert z-z_k\rvert\le r$ are disjoint and lie inside $C$. Let $\Omega$ be the region inside $C$ and outside all the discs. Its boundary is $C$ (positively oriented) together with the circles $C_k$ traversed negatively, and $f$ is analytic on a neighbourhood of $\overline\Omega$ (with continuous derivative, by Remark 1 of §D.14). By the general form of Cauchy's theorem,

$$
\oint_Cf\,dz-\sum_{k=1}^{K}\oint_{C_k}f\,dz=0 .
$$

Substituting Step 1 completes the proof. ∎

**Remarks.**

1. *Computing residues.* At a simple pole, $\operatorname{Res}=\lim_{z\to z_k}(z-z_k)f(z)$; for $f=N/D$ with $D(z_k)=0\neq D'(z_k)$ this is $N(z_k)/D'(z_k)$; at a pole of order $m$, $\operatorname{Res}=\frac{1}{(m-1)!}\lim_{z\to z_k}\frac{d^{m-1}}{dz^{m-1}}\big[(z-z_k)^mf(z)\big]$. Each follows by reading off $a_{-1}$ from the Laurent series ([Ch. 18, §18.7](../part-5-transforms/18-complex-analysis.md#187-residues)).
2. *Hypotheses.* The singularities must be isolated and finite in number inside $C$, and none may lie *on* $C$. A pole on the contour makes the integral improper; one then uses indentations and principal values, each small semicircle contributing $\pi j$ times the residue in the limit.
3. *Where it is used.* Inverse Laplace transforms by closing the Bromwich contour to the left (the Heaviside expansion is the sum of residues of $X(s)e^{st}$), inverse Z-transforms $x[n]=\frac{1}{2\pi j}\oint X(z)z^{n-1}dz$, real integrals ([Ch. 18, §18.8](../part-5-transforms/18-complex-analysis.md#188-residue-theorem)), and the argument principle behind the Nyquist criterion ([§18.9](../part-5-transforms/18-complex-analysis.md#189-applications-to-engineering-analysis)).

## D.16 Jordan's lemma

> [!theorem] Jordan's inequality and Jordan's lemma
> **(Inequality.)** For $0\le\theta\le\pi/2$, $\dfrac{2\theta}{\pi}\le\sin\theta\le\theta$.
>
> **(Lemma.)** Let $\omega>0$, let $C_R$ be the upper semicircle $z=Re^{j\theta}$, $0\le\theta\le\pi$, and let $g$ be continuous on $C_R$ for all large $R$ with $M_R=\max_{C_R}\lvert g\rvert\to0$ as $R\to\infty$. Then
>
> $$
> \left\lvert\int_{C_R}g(z)\,e^{j\omega z}\,dz\right\rvert\le\frac{\pi}{\omega}\,M_R\big(1-e^{-\omega R}\big)\le\frac{\pi M_R}{\omega}\ \xrightarrow[R\to\infty]{}\ 0 .
> $$

**Proof of the inequality.** $\sin\theta\le\theta$ because $\frac{d}{d\theta}(\theta-\sin\theta)=1-\cos\theta\ge0$ and equality holds at $0$. For the lower bound, $\sin$ is concave on $[0,\pi/2]$ (its second derivative $-\sin\theta\le0$), so its graph lies above the chord joining $(0,0)$ and $(\pi/2,1)$, which is the line $2\theta/\pi$. ∎

**Proof of the lemma.** On $C_R$, $z=R(\cos\theta+j\sin\theta)$, so $\lvert e^{j\omega z}\rvert=e^{-\omega R\sin\theta}$ and $\lvert dz\rvert=R\,d\theta$. Hence

$$
\left\lvert\int_{C_R}g\,e^{j\omega z}dz\right\rvert\le M_R\,R\int_0^{\pi}e^{-\omega R\sin\theta}\,d\theta=2M_RR\int_0^{\pi/2}e^{-\omega R\sin\theta}\,d\theta,
$$

using the symmetry $\sin(\pi-\theta)=\sin\theta$. By Jordan's inequality $e^{-\omega R\sin\theta}\le e^{-2\omega R\theta/\pi}$ on $[0,\pi/2]$, so

$$
2M_RR\int_0^{\pi/2}e^{-2\omega R\theta/\pi}\,d\theta=2M_RR\cdot\frac{\pi}{2\omega R}\big(1-e^{-\omega R}\big)=\frac{\pi M_R}{\omega}\big(1-e^{-\omega R}\big).
$$

∎

**Remarks.**

1. *Why the lemma is sharper than ML.* The ML inequality gives only $\lvert\int\rvert\le\pi RM_R$, which tends to zero only if $M_R=o(1/R)$, e.g. $g=1/z^2$. Jordan's lemma needs only $M_R\to0$, e.g. $g=1/z$: the factor $R$ from the arc length is cancelled because $e^{j\omega z}$ is exponentially small except within an angle of order $1/(\omega R)$ of the real axis. Numerically, for $\omega=1$ the key integral $\int_0^\pi e^{-R\sin\theta}d\theta$ equals $0.4208$ at $R=5$ and $0.0400$ at $R=50$, against Jordan's bounds $\frac{\pi}{R}(1-e^{-R})=0.6241$ and $0.0628$: both decay like $1/R$, as the proof predicts.
2. *The sign of $\omega$ dictates the half-plane.* For $\omega<0$ the factor $e^{j\omega z}$ is small in the *lower* half-plane, and the contour must be closed there (clockwise, introducing a minus sign). The same rule closes the Bromwich contour for $X(s)e^{st}$: to the left for $t>0$, to the right for $t<0$, which is why causal signals have all their poles to the left of the ROC.
3. *Typical use.* For $\int_{-\infty}^{\infty}\frac{\cos x}{x^2+1}dx$, integrate $e^{jz}/(z^2+1)$ over $[-R,R]$ plus $C_R$. The only enclosed pole is $z=j$, with residue $e^{-1}/(2j)$; the arc vanishes by Jordan's lemma (here even ML suffices); so the integral is $2\pi j\cdot e^{-1}/(2j)=\pi/e=1.15573$, which agrees with numerical quadrature. For $\int_{-\infty}^\infty\frac{\sin x}{x}dx=\pi$, the case $g=1/z$, ML fails and Jordan's lemma is indispensable. See [Ch. 18, §18.8](../part-5-transforms/18-complex-analysis.md#188-residue-theorem) and the Fourier inversion of rational spectra in [Ch. 16, §16.2](../part-5-transforms/16-fourier-transforms.md#162-inverse-transform).

```desmos-graph
left=-0.1; right=1.75; bottom=-0.1; top=1.15
---
y=\sin(x)|blue|0<x<1.5708
y=\frac{2x}{\pi}|red|dashed|0<x<1.5708
y=x|grey|dotted|0<x<1.2
(1.5708,1)|label:(π/2, 1)
```

*Jordan's inequality: on $[0,\pi/2]$ the concave curve $\sin\theta$ (blue) lies between its chord $2\theta/\pi$ (red, dashed) and its tangent $\theta$ (grey). The lower bound is what turns $e^{-\omega R\sin\theta}$ into an integrable exponential.*

## Further reading

- E. A. Coddington and N. Levinson, *Theory of Ordinary Differential Equations*, McGraw-Hill — existence, uniqueness and dependence on data (§D.1).
- W. Rudin, *Real and Complex Analysis*, 3rd ed., McGraw-Hill — Fourier analysis on $L^1$ and $L^2$, Plancherel's theorem, Cauchy's theorem in full generality (§§D.7, D.10, D.14).
- T. W. Körner, *Fourier Analysis*, Cambridge University Press — Fejér, Dirichlet, Gibbs, Poisson summation and sampling, with history (§§D.6–D.11).
- D. V. Widder, *The Laplace Transform*, Princeton University Press — existence, Abelian and Tauberian theorems (§§D.2–D.5).
- A. V. Oppenheim and R. W. Schafer, *Discrete-Time Signal Processing*, Pearson — the Z-transform, ROC and stability (§D.12).
- J. W. Brown and R. V. Churchill, *Complex Variables and Applications*, McGraw-Hill — Cauchy–Riemann, Cauchy–Goursat, residues and Jordan's lemma (§§D.13–D.16).
