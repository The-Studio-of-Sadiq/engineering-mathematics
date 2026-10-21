# Appendix A — Calculus and Complex Numbers

[← Table of Contents](../00-front/table-of-contents.md)

> This appendix gives complete proofs of the single-variable calculus results, the complex-exponential results and the two multivariable results that the chapters of Parts I and II state with only a sketch. The results are proved in logical order, not in the order the chapters use them: the mean value theorem comes first because almost everything else rests on it, and Euler's formula comes after Taylor's theorem because one of its proofs needs the Taylor series of $\cos$ and $\sin$.
>
> **What is assumed.** The real numbers are complete (every non-empty set bounded above has a least upper bound). From completeness we take, without proof, three standard facts: the **extreme value theorem** (a continuous function on a closed bounded interval $[a,b]$ attains a maximum and a minimum), the **intermediate value theorem** (a continuous function on $[a,b]$ takes every value between $f(a)$ and $f(b)$), and the **Riemann integrability of continuous functions** on $[a,b]$, together with the elementary properties of the integral (linearity, additivity over intervals, and monotonicity: $f\le g\Rightarrow\int f\le\int g$). We also use the differentiation rules (sum, product, quotient, chain) of [Ch. 4 §4.2–4.3](../part-2-calculus/04-differential-calculus.md), the derivatives $\sin'=\cos$, $\cos'=-\sin$, and the identity $\cos^2\theta+\sin^2\theta=1$. The implicit function theorem is stated precisely and used, but not proved.
>
> **Level of rigour.** Every hypothesis is stated; every step of every proof is justified from the assumptions above or from an earlier result of this appendix. Where a hypothesis cannot be dropped, a counterexample is given. Proofs end with ∎.
>
> **Complex-valued functions.** Several results are applied to functions $F:[a,b]\to\mathbb{C}$ of a real variable, such as $e^{j\omega t}$. For these, $F=U+jV$ with $U=\operatorname{Re}\{F\}$, $V=\operatorname{Im}\{F\}$ real, and limits, continuity, derivatives and integrals are taken component by component: $F'=U'+jV'$, $\int F=\int U+j\int V$. Because complex addition and multiplication obey the same algebraic laws as real ones and $\lvert z\rvert$ obeys the triangle inequality, the sum, product, quotient and chain rules (for a real inner function) carry over word for word.

## A.1 Rolle's theorem, the mean value theorem and Cauchy's mean value theorem

The mean value theorem (MVT) is the bridge between local information (a derivative at a point) and global information (a change over an interval). It is the tool used to prove the FTC, Taylor's theorem, L'Hôpital's rule, the uniqueness half of Euler's formula and Clairaut's theorem below.

We first need the fact that a derivative vanishes at an interior extremum.

> [!theorem] Lemma A.1 (Fermat)
> Let $f$ be defined on an open interval $I$, let $c\in I$, and suppose $f(c)\ge f(x)$ for all $x\in I$ (or $f(c)\le f(x)$ for all $x\in I$). If $f'(c)$ exists, then $f'(c)=0$.

**Proof.** Take the maximum case; the minimum case follows by applying it to $-f$. For $h>0$ small enough that $c+h\in I$, $f(c+h)-f(c)\le0$, so

$$
\frac{f(c+h)-f(c)}{h}\le 0\quad(h>0),\qquad \frac{f(c+h)-f(c)}{h}\ge 0\quad(h<0).
$$

The derivative is the limit of these quotients as $h\to0$; the limit of the right-hand quotients ($h\to0^+$) is $\le0$ and that of the left-hand ones ($h\to0^-$) is $\ge0$ (limits preserve non-strict inequalities). Both one-sided limits equal $f'(c)$, so $0\le f'(c)\le0$. ∎

> [!theorem] Theorem A.1 (Rolle, MVT, Cauchy MVT)
> Let $f$ and $g$ be real-valued, continuous on $[a,b]$ and differentiable on $(a,b)$, with $a<b$.
>
> 1. **(Rolle)** If $f(a)=f(b)$, there is $c\in(a,b)$ with $f'(c)=0$.
> 2. **(Mean value theorem)** There is $c\in(a,b)$ with
> $$
> f'(c)=\frac{f(b)-f(a)}{b-a}.
> $$
> 3. **(Cauchy's mean value theorem)** There is $c\in(a,b)$ with
> $$
> \big[f(b)-f(a)\big]\,g'(c)=\big[g(b)-g(a)\big]\,f'(c).
> $$
> If in addition $g'(x)\ne0$ on $(a,b)$, then $g(b)\ne g(a)$ and $\dfrac{f(b)-f(a)}{g(b)-g(a)}=\dfrac{f'(c)}{g'(c)}$.

**Proof.**

*Rolle.* By the extreme value theorem, $f$ attains a maximum value $M$ and a minimum value $m$ on $[a,b]$.

- If $M=m$, $f$ is constant on $[a,b]$, so $f'(c)=0$ for every $c\in(a,b)$.
- If $M>m$, at least one of $M,m$ differs from the common endpoint value $f(a)=f(b)$. Say $M\ne f(a)$ (the other case is symmetric). Then the maximum is attained at some $c$ with $c\ne a$, $c\ne b$, i.e. $c\in(a,b)$. On the open interval $(a,b)$, $f(c)\ge f(x)$ for all $x$, and $f'(c)$ exists, so $f'(c)=0$ by Lemma A.1.

*MVT.* Let $s=\dfrac{f(b)-f(a)}{b-a}$ (the slope of the chord) and define

$$
\varphi(x)=f(x)-f(a)-s\,(x-a).
$$

$\varphi$ is continuous on $[a,b]$ and differentiable on $(a,b)$ (a difference of such functions), with $\varphi(a)=0$ and $\varphi(b)=f(b)-f(a)-s(b-a)=0$. Rolle gives $c\in(a,b)$ with $0=\varphi'(c)=f'(c)-s$.

*Cauchy MVT.* Define

$$
h(x)=\big[f(b)-f(a)\big]g(x)-\big[g(b)-g(a)\big]f(x).
$$

$h$ is continuous on $[a,b]$, differentiable on $(a,b)$, and

$$
h(a)=f(b)g(a)-f(a)g(a)-g(b)f(a)+g(a)f(a)=f(b)g(a)-g(b)f(a),
$$

$$
h(b)=f(b)g(b)-f(a)g(b)-g(b)f(b)+g(a)f(b)=f(b)g(a)-g(b)f(a).
$$

So $h(a)=h(b)$, and Rolle gives $c$ with $h'(c)=[f(b)-f(a)]g'(c)-[g(b)-g(a)]f'(c)=0$, which is the claim. Finally, if $g'\ne0$ on $(a,b)$ then $g(b)\ne g(a)$, for otherwise Rolle applied to $g$ would produce a zero of $g'$; dividing by $[g(b)-g(a)]g'(c)\ne0$ gives the quotient form. ∎

The MVT is the special case $g(x)=x$ of Cauchy's theorem. Geometrically, Cauchy's theorem says that a parametric curve $t\mapsto(g(t),f(t))$ has, somewhere, a tangent parallel to the chord joining its ends.

```desmos-graph
left=-0.5; right=2.5; bottom=-2.6; top=2.6
---
y=x^{3}-3x|blue|0<x<2
y=x|red|dashed|0<x<2
y=(x-1.1547)-1.9245|green|0.3<x<2.2
(0,0)|label:a
(2,2)|label:b
(1.1547,-1.9245)|label:c
```

*Mean value theorem for $f(x)=x^3-3x$ on $[0,2]$: the chord (red, dashed) has slope $1$, and the tangent at $c=2/\sqrt3\approx1.1547$ (green) is parallel to it, since $f'(c)=3c^2-3=1$.*

> [!theorem] Corollary A.1 (consequences of the MVT)
> Let $I$ be an interval and $f$ differentiable on $I$.
> 1. If $f'(x)=0$ for all $x\in I$, then $f$ is constant on $I$. This holds also for complex-valued $f$.
> 2. If $f'=g'$ on $I$, then $f-g$ is constant on $I$.
> 3. If $f'>0$ on $I$ then $f$ is strictly increasing on $I$; if $f'\ge0$, non-decreasing. (Similarly for $<$, $\le$.)
> 4. If $\lvert f'(x)\rvert\le M$ on $I$, then $\lvert f(x_2)-f(x_1)\rvert\le M\lvert x_2-x_1\rvert$ for all $x_1,x_2\in I$.

**Proof.** For $x_1<x_2$ in $I$, $f$ is continuous on $[x_1,x_2]$ (differentiable implies continuous) and differentiable on $(x_1,x_2)$, so $f(x_2)-f(x_1)=f'(c)(x_2-x_1)$ for some $c\in(x_1,x_2)$. Items 1, 3 and 4 for real $f$ read off directly from this formula. For complex $f=U+jV$, $f'=0$ means $U'=V'=0$, so $U$ and $V$ are each constant, hence so is $f$. Item 2 is item 1 applied to $f-g$. ∎

**Remarks.**

- *Differentiability in the interior is essential.* $f(x)=\lvert x\rvert$ on $[-1,1]$ is continuous with $f(-1)=f(1)$, yet $f'(x)=\pm1$ wherever it exists: there is no $c$ with $f'(c)=0$. The failure is the single point $x=0$ where $f'$ does not exist.
- *Continuity at the endpoints is essential.* $f(x)=x$ on $[0,1)$ with $f(1)=0$ is differentiable on $(0,1)$ and has $f(0)=f(1)$, but $f'\equiv1$.
- *The MVT itself fails for complex-valued functions.* $F(t)=e^{jt}$ on $[0,2\pi]$ has $F(2\pi)-F(0)=0$, but $F'(t)=je^{jt}$ has modulus $1$ everywhere, so no $c$ gives $F'(c)\cdot2\pi=0$. (Euler's formula, used here only as an illustration, is proved in §A.6 without this example.) What survives for complex functions is Corollary A.1(1), proved through real and imaginary parts, and the inequality form (4), which can be proved from the FTC (§A.2) as $\lvert F(x_2)-F(x_1)\rvert=\lvert\int_{x_1}^{x_2}F'\rvert\le M\lvert x_2-x_1\rvert$ when $F'$ is continuous. This is why the proof of Euler's formula in §A.6 is arranged around "zero derivative implies constant".
- *Where it is used.* The MVT and its corollaries are stated in [Ch. 4 §4.2](../part-2-calculus/04-differential-calculus.md); Corollary A.1(2) is the justification of the "$+C$" of indefinite integration in [Ch. 5 §5.1.4](../part-2-calculus/05-integral-calculus.md), and item 4 is the Lipschitz estimate behind the existence–uniqueness theory of [Ch. 11 §11.7](../part-4-differential-equations/11-ordinary-differential-equations.md) and the error analysis of fixed-point and Newton iterations in [Ch. 22](../part-7-numerical-computational/22-numerical-methods.md).

## A.2 The Fundamental Theorem of Calculus

The FTC connects the two halves of calculus: the integral, defined as a limit of Riemann sums ([Ch. 5 §5.1.2](../part-2-calculus/05-integral-calculus.md)), and the derivative. We recall the definition so that the proof of Part 2 can refer to it.

> [!definition] Riemann integral
> A **partition** $P$ of $[a,b]$ is a finite set of points $a=x_0<x_1<\dots<x_n=b$, with mesh $\lVert P\rVert=\max_i(x_i-x_{i-1})$. A **Riemann sum** for $f$ is $\sum_{i=1}^n f(c_i)\,\Delta x_i$ with $\Delta x_i=x_i-x_{i-1}$ and any tags $c_i\in[x_{i-1},x_i]$. A bounded function $f$ is **Riemann integrable** on $[a,b]$, with integral $\int_a^b f$, if for every $\varepsilon>0$ there is $\delta>0$ such that every Riemann sum with $\lVert P\rVert<\delta$ (and any choice of tags) lies within $\varepsilon$ of $\int_a^b f$. We set $\int_b^a f=-\int_a^b f$ and $\int_a^a f=0$.

We first need the mean value theorem for integrals.

> [!theorem] Lemma A.2 (Mean value theorem for integrals)
> If $f$ is continuous on the closed interval with end-points $p\ne q$, there is $c$ between $p$ and $q$ (inclusive) with
> $$
> \frac{1}{q-p}\int_p^q f(t)\,dt=f(c).
> $$

**Proof.** Suppose $p<q$ (otherwise swap them; both sides of the formula are unchanged). By the extreme value theorem $f$ has a minimum $m=f(t_1)$ and a maximum $M=f(t_2)$ on $[p,q]$. Monotonicity of the integral applied to $m\le f\le M$ gives $m(q-p)\le\int_p^qf\le M(q-p)$, so the average $\mu=\frac1{q-p}\int_p^qf$ satisfies $m\le\mu\le M$. By the intermediate value theorem applied to $f$ on the closed interval between $t_1$ and $t_2$, there is $c$ in that interval, hence in $[p,q]$, with $f(c)=\mu$. ∎

> [!theorem] Theorem A.2 (Fundamental Theorem of Calculus)
> **Part 1.** Let $f$ be Riemann integrable on $[a,b]$ and define $A(x)=\displaystyle\int_a^x f(t)\,dt$ for $x\in[a,b]$. Then $A$ is continuous on $[a,b]$. If moreover $f$ is continuous at a point $x_0\in[a,b]$, then $A$ is differentiable at $x_0$ (one-sidedly if $x_0$ is an end-point) and $A'(x_0)=f(x_0)$. In particular, if $f$ is continuous on $[a,b]$, then $A'=f$ on $[a,b]$: every continuous function has an antiderivative.
>
> **Part 2.** Let $F$ be continuous on $[a,b]$ and differentiable on $(a,b)$, and suppose $F'=f$ on $(a,b)$, where $f$ is Riemann integrable on $[a,b]$. Then
> $$
> \int_a^b f(x)\,dx=F(b)-F(a).
> $$

**Proof of Part 1.** Since $f$ is integrable it is bounded, say $\lvert f\rvert\le K$. For $x,x+h\in[a,b]$, additivity of the integral gives

$$
A(x+h)-A(x)=\int_x^{x+h}f(t)\,dt,
$$

and monotonicity gives $\lvert A(x+h)-A(x)\rvert\le K\lvert h\rvert$, so $A$ is (Lipschitz) continuous.

Now let $f$ be continuous at $x_0$ and let $\varepsilon>0$. Choose $\delta>0$ with $\lvert f(t)-f(x_0)\rvert<\varepsilon$ whenever $\lvert t-x_0\rvert<\delta$, $t\in[a,b]$. For $0<\lvert h\rvert<\delta$ with $x_0+h\in[a,b]$, write $f(x_0)=\frac1h\int_{x_0}^{x_0+h}f(x_0)\,dt$; then

$$
\left\lvert\frac{A(x_0+h)-A(x_0)}{h}-f(x_0)\right\rvert=\left\lvert\frac1h\int_{x_0}^{x_0+h}\big[f(t)-f(x_0)\big]\,dt\right\rvert\le\frac1{\lvert h\rvert}\cdot\lvert h\rvert\,\varepsilon=\varepsilon,
$$

because every $t$ between $x_0$ and $x_0+h$ satisfies $\lvert t-x_0\rvert<\delta$. Hence the difference quotient tends to $f(x_0)$. (When $f$ is continuous on the whole interval one can instead quote Lemma A.2: the difference quotient equals $f(c_h)$ with $c_h$ between $x_0$ and $x_0+h$, and $f(c_h)\to f(x_0)$ — this is the version in Ch. 5.) ∎

**Proof of Part 2.** Let $\varepsilon>0$ and choose $\delta$ as in the definition of $\int_a^bf$. Take any partition $a=x_0<\dots<x_n=b$ with mesh $<\delta$. On each sub-interval $[x_{i-1},x_i]$, $F$ is continuous, and differentiable in the interior, so by the MVT (Theorem A.1) there is $c_i\in(x_{i-1},x_i)$ with $F(x_i)-F(x_{i-1})=F'(c_i)\Delta x_i=f(c_i)\Delta x_i$. Summing, the left side telescopes:

$$
F(b)-F(a)=\sum_{i=1}^n\big[F(x_i)-F(x_{i-1})\big]=\sum_{i=1}^n f(c_i)\,\Delta x_i .
$$

The right side is a Riemann sum with mesh $<\delta$, so it is within $\varepsilon$ of $\int_a^bf$. Thus $\lvert F(b)-F(a)-\int_a^bf\rvert<\varepsilon$ for every $\varepsilon>0$, i.e. the two are equal. ∎

The proof of Part 2 does not use Part 1. For continuous $f$ there is also the two-line proof of Ch. 5: $A$ and $F$ both have derivative $f$ on $(a,b)$, so $F=A+C$ by Corollary A.1(2); both are continuous on $[a,b]$, so the identity extends to the end-points, and $F(b)-F(a)=A(b)-A(a)=\int_a^bf$. The version above is stronger: it needs only integrability of $F'$, not continuity.

**Remarks.**

- *Part 1 needs continuity at the point.* With $f(t)=-1$ for $t<0$ and $f(t)=+1$ for $t\ge0$ on $[-1,1]$, $A(x)=\int_{-1}^xf=\lvert x\rvert-1$, which is continuous (as Part 1 promises for any integrable $f$) but not differentiable at $0$, the jump of $f$. In circuit language: the charge $q(t)=\int i\,dt$ on a capacitor is continuous even when the current switches abruptly, but $dq/dt$ fails to exist at the switching instant.
- *Part 2 needs $F'=f$ on the whole interval.* The formal computation $\int_{-1}^1\frac{dx}{x^2}=\big[-\frac1x\big]_{-1}^{1}=-2$ is absurd for a positive integrand. Here $F=-1/x$ is not continuous (not even defined) on $[-1,1]$, and $1/x^2$ is unbounded, hence not Riemann integrable; the integral is improper and divergent ([Ch. 5 §5.3](../part-2-calculus/05-integral-calculus.md)).
- *Part 2 needs $F'$ to be integrable.* $F(x)=x^2\sin(1/x^2)$, $F(0)=0$, is differentiable at every point of $[0,1]$, but $F'(x)=2x\sin(1/x^2)-\frac2x\cos(1/x^2)$ for $x\ne0$ is unbounded near $0$, so the Riemann integral of $F'$ over $[0,1]$ does not exist even though $F(1)-F(0)=\sin1$ is perfectly finite.
- *Complex-valued integrands.* Applying both parts to the real and imaginary parts shows that they hold verbatim for $f:[a,b]\to\mathbb{C}$; for example $\int_0^T e^{j\omega t}\,dt=\big[e^{j\omega t}/(j\omega)\big]_0^T$ for $\omega\ne0$, because $\frac{d}{dt}e^{j\omega t}=j\omega e^{j\omega t}$ (§A.6).
- *Where it is used.* The FTC is stated in [Ch. 5 §5.1.5](../part-2-calculus/05-integral-calculus.md). Part 1 underlies every "accumulated quantity" of [Ch. 5 §5.8](../part-2-calculus/05-integral-calculus.md) (charge from current, energy from power) and the integrating-factor method of [Ch. 11 §11.3](../part-4-differential-equations/11-ordinary-differential-equations.md); Part 2 is the evaluation rule used in every chapter.

## A.3 Integration by parts and integration by substitution

Each of the two main integration techniques is a differentiation rule read backwards through FTC Part 2: parts reverses the product rule, substitution reverses the chain rule. The proofs are short, but the hypotheses are where students (and computer algebra systems) go wrong.

> [!theorem] Theorem A.3 (Integration by parts)
> Let $u$ and $v$ be continuously differentiable ($C^1$) on $[a,b]$, real- or complex-valued. Then
> $$
> \int_a^b u(x)\,v'(x)\,dx=\Big[u(x)\,v(x)\Big]_a^b-\int_a^b u'(x)\,v(x)\,dx .
> $$
> On any interval where $u,v$ are $C^1$, the indefinite form $\int u\,v'\,dx=uv-\int u'v\,dx$ holds (each side an antiderivative of $uv'$, up to an additive constant).

**Proof.** By the product rule, $(uv)'=u'v+uv'$ on $[a,b]$ (one-sided at the end-points), and the right-hand side is continuous because $u,v,u',v'$ are. Continuous functions are integrable, so FTC Part 2 applied to $F=uv$ gives

$$
\Big[uv\Big]_a^b=\int_a^b (uv)'\,dx=\int_a^b u'v\,dx+\int_a^b uv'\,dx,
$$

using linearity of the integral. Rearranging gives the formula. For the indefinite form: if $G$ is an antiderivative of $u'v$ (it exists by FTC Part 1), then $(uv-G)'=u'v+uv'-u'v=uv'$, so $uv-G$ is an antiderivative of $uv'$. For complex-valued $u,v$ the product rule holds by the same algebra, and FTC applies componentwise. ∎

> [!theorem] Theorem A.4 (Integration by substitution)
> Let $g:[\alpha,\beta]\to\mathbb{R}$ be $C^1$, and let $f$ be continuous on an interval $I$ that contains $g([\alpha,\beta])$. Then
> $$
> \int_\alpha^\beta f\big(g(t)\big)\,g'(t)\,dt=\int_{g(\alpha)}^{g(\beta)}f(x)\,dx .
> $$
> $g$ need **not** be one-to-one.
>
> **Inverse form.** If in addition $g$ is strictly monotonic on $[\alpha,\beta]$, then for $a=g(\alpha)$, $b=g(\beta)$ every integral $\int_a^b f(x)\,dx$ may be computed as $\int_{g^{-1}(a)}^{g^{-1}(b)}f(g(t))\,g'(t)\,dt$, i.e. by putting $x=g(t)$, $dx=g'(t)\,dt$ and transforming the limits with $g^{-1}$.

**Proof.** Fix $x_0\in I$ and let $F(x)=\int_{x_0}^xf(s)\,ds$ for $x\in I$. Since $f$ is continuous on the interval $I$, FTC Part 1 gives $F'=f$ on $I$. The composite $F\circ g$ is defined on $[\alpha,\beta]$ because $g$ takes values in $I$, and by the chain rule

$$
\frac{d}{dt}F\big(g(t)\big)=F'\big(g(t)\big)\,g'(t)=f\big(g(t)\big)\,g'(t),
$$

which is continuous on $[\alpha,\beta]$ (composition and product of continuous functions). FTC Part 2, applied twice, gives

$$
\int_\alpha^\beta f\big(g(t)\big)g'(t)\,dt=F\big(g(\beta)\big)-F\big(g(\alpha)\big)=\int_{g(\alpha)}^{g(\beta)}f(x)\,dx .
$$

The inverse form is the same identity read from right to left; monotonicity is needed only so that $t=g^{-1}(x)$ is a well-defined function and the new limits $g^{-1}(a)=\alpha$, $g^{-1}(b)=\beta$ are determined uniquely. ∎

**Remarks.**

- *Non-injective substitutions are fine in the forward direction.* With $g(t)=t^2$ on $[-1,2]$ (which is not one-to-one: it folds $[-1,0]$ back over $[0,1]$),
$$
\int_{-1}^{2}2t\cos(t^2)\,dt=\int_{1}^{4}\cos x\,dx=\sin4-\sin1\approx-1.59827 ,
$$
which agrees with direct numerical quadrature to all printed digits. The folded part is traversed twice in opposite directions and cancels — the theorem accounts for orientation automatically, which is why no absolute value $\lvert g'\rvert$ appears (contrast the multivariable change of variables of [Ch. 6 §6.8](../part-2-calculus/06-multivariable-calculus.md), which uses $\lvert\det\mathbf{J}\rvert$ and does require injectivity).
- *The hypothesis that $g$ is $C^1$ on the whole interval is essential.* The tangent half-angle substitution $x=\tan(t/2)$ of [Ch. 5 §5.6.5](../part-2-calculus/05-integral-calculus.md) turns $\int\frac{dt}{2+\cos t}$ into $\int\frac{2\,dx}{3+x^2}$. Applied blindly to $\int_0^{2\pi}$, it maps both limits to $x=\tan0=\tan\pi=0$ and gives $0$, which is impossible for a positive integrand. The true value is $2\pi/\sqrt3\approx3.62760$ (checked by numerical quadrature). The error: $\tan(t/2)$ is not even defined at $t=\pi$, so the hypotheses fail. The correct procedure is to split at $t=\pi$; each half maps onto a half-line of $x$ and gives an improper integral equal to $\pi/\sqrt3$.
- *The hypothesis on $f$ concerns the whole range of $g$.* $f$ must be continuous on an interval containing all values $g(t)$, not merely on the interval between $g(\alpha)$ and $g(\beta)$. If $g$ wanders outside that interval into a region where $f$ is singular, the left side may diverge while the right side is finite.
- *Integration by parts with a jump.* If $u$ or $v$ has a jump discontinuity inside $[a,b]$, split the interval at the jump; the boundary terms then contribute the jump values. This is the "jump method" for Fourier coefficients in [Ch. 15 §15.7](../part-5-transforms/15-fourier-series.md).
- *Where they are used.* The rules are stated in [Ch. 5 §5.5–5.6](../part-2-calculus/05-integral-calculus.md). Parts gives the Laplace-transform derivative rule $\mathcal{L}\{x'\}=sX(s)-x(0^-)$ ([Ch. 14](../part-5-transforms/14-laplace-transforms.md)), the decay of Fourier coefficients with smoothness ([Ch. 15](../part-5-transforms/15-fourier-series.md)), the symmetry of Sturm–Liouville operators ([Ch. 11 §11.8](../part-4-differential-equations/11-ordinary-differential-equations.md)) and the integral form of Taylor's remainder (§A.4). Linear substitution $x=at+b$ gives the time-shift and time-scaling theorems of every transform in Part V.

## A.4 Taylor's theorem with Lagrange and integral remainders

Taylor's theorem quantifies how well the Taylor polynomial

$$
p_n(x)=\sum_{k=0}^{n}\frac{f^{(k)}(a)}{k!}(x-a)^k=f(a)+f'(a)(x-a)+\dots+\frac{f^{(n)}(a)}{n!}(x-a)^n
$$

approximates $f$ near $a$ ([Ch. 4 §4.7](../part-2-calculus/04-differential-calculus.md)). The **remainder** is $R_n(x)=f(x)-p_n(x)$. Two exact formulas for it are proved here: Lagrange's, which needs the weakest hypotheses and gives the familiar error bound, and the integral form, which is exact without an unknown intermediate point.

> [!theorem] Theorem A.5 (Taylor's theorem)
> Let $I$ be an open interval, $a,x\in I$, and $n\ge0$ an integer.
>
> **(a) Lagrange remainder.** If $f$ is $n+1$ times differentiable on $I$, then there is a point $c$ strictly between $a$ and $x$ (take $c=a$ if $x=a$) such that
> $$
> R_n(x)=\frac{f^{(n+1)}(c)}{(n+1)!}(x-a)^{n+1}.
> $$
> Consequently, if $\lvert f^{(n+1)}\rvert\le M$ on the interval between $a$ and $x$, then $\lvert R_n(x)\rvert\le\dfrac{M\,\lvert x-a\rvert^{n+1}}{(n+1)!}$.
>
> **(b) Integral remainder.** If $f$ is $C^{n+1}$ on $I$ (that is, $f^{(n+1)}$ exists and is continuous), then
> $$
> R_n(x)=\frac{1}{n!}\int_a^x f^{(n+1)}(t)\,(x-t)^n\,dt .
> $$

**Proof of (a).** The case $x=a$ is trivial ($R_n(a)=0$), so fix $x\ne a$. Define the constant $K$ by

$$
f(x)=p_n(x)+\frac{K}{(n+1)!}(x-a)^{n+1},
$$

which is possible because $(x-a)^{n+1}\ne0$. We must show $K=f^{(n+1)}(c)$ for some $c$ between $a$ and $x$. The trick is to treat the expansion point as the variable: for $t\in I$ define

$$
G(t)=f(x)-\sum_{k=0}^{n}\frac{f^{(k)}(t)}{k!}(x-t)^k-\frac{K}{(n+1)!}(x-t)^{n+1}.
$$

*Values at the ends.* At $t=x$ every term with $k\ge1$ and the last term vanish, and the $k=0$ term is $f(x)$, so $G(x)=0$. At $t=a$, the sum is $p_n(x)$, so $G(a)=f(x)-p_n(x)-\frac{K}{(n+1)!}(x-a)^{n+1}=0$ by the choice of $K$.

*Derivative.* Each $f^{(k)}$ with $k\le n$ is differentiable on $I$, so $G$ is differentiable on $I$ (hence continuous). By the product rule, for $k\ge1$,

$$
\frac{d}{dt}\left[\frac{f^{(k)}(t)}{k!}(x-t)^k\right]=\frac{f^{(k+1)}(t)}{k!}(x-t)^k-\frac{f^{(k)}(t)}{(k-1)!}(x-t)^{k-1},
$$

and for $k=0$ the derivative is $f'(t)$, which is the first term of the same formula with no second term. Summing over $k=0,\dots,n$, the second part of term $k$ cancels the first part of term $k-1$; the sum telescopes to its last surviving piece, $\frac{f^{(n+1)}(t)}{n!}(x-t)^n$. Therefore

$$
G'(t)=-\frac{f^{(n+1)}(t)}{n!}(x-t)^n+\frac{K}{n!}(x-t)^n=\frac{(x-t)^n}{n!}\Big[K-f^{(n+1)}(t)\Big].
$$

*Rolle.* $G$ is continuous on the closed interval between $a$ and $x$, differentiable inside it, and vanishes at both ends. By Rolle's theorem (Theorem A.1) there is $c$ strictly between $a$ and $x$ with $G'(c)=0$. Since $c\ne x$, $(x-c)^n\ne0$, so $K=f^{(n+1)}(c)$. The bound follows by taking absolute values. ∎

**Proof of (b).** Induction on $n$.

*Base case $n=0$.* The claim is $f(x)-f(a)=\int_a^xf'(t)\,dt$, which is FTC Part 2 (Theorem A.2) applied to the continuous function $f'$.

*Inductive step.* Suppose the formula holds for $n-1$ (with $n\ge1$) and that $f$ is $C^{n+1}$; then in particular $f$ is $C^n$ and

$$
R_{n-1}(x)=\frac1{(n-1)!}\int_a^xf^{(n)}(t)(x-t)^{n-1}\,dt .
$$

Integrate by parts (Theorem A.3) with $u(t)=f^{(n)}(t)$ and $v(t)=-\dfrac{(x-t)^n}{n!}$, so that $v'(t)=\dfrac{(x-t)^{n-1}}{(n-1)!}$; both $u$ and $v$ are $C^1$ because $f$ is $C^{n+1}$. Then

$$
R_{n-1}(x)=\left[-f^{(n)}(t)\frac{(x-t)^n}{n!}\right]_{t=a}^{t=x}+\int_a^x f^{(n+1)}(t)\frac{(x-t)^n}{n!}\,dt=\frac{f^{(n)}(a)}{n!}(x-a)^n+\frac1{n!}\int_a^xf^{(n+1)}(t)(x-t)^n\,dt .
$$

Since $R_n(x)=R_{n-1}(x)-\frac{f^{(n)}(a)}{n!}(x-a)^n$ (the degree-$n$ polynomial has one more term), this is exactly the formula for $R_n$. ∎

*From (b) back to (a).* Under the stronger hypothesis $f\in C^{n+1}$, (a) also follows from (b) via the weighted mean value theorem for integrals: $(x-t)^n$ does not change sign between $a$ and $x$, so $\int_a^xf^{(n+1)}(t)(x-t)^n\,dt=f^{(n+1)}(c)\int_a^x(x-t)^n\,dt=f^{(n+1)}(c)\frac{(x-a)^{n+1}}{n+1}$ for some $c$ (the proof is that of Lemma A.2, with $m\le f^{(n+1)}\le M$ multiplied by the one-signed weight). Dividing by $n!$ gives Lagrange's form.

> [!theorem] Corollary A.5 (The exponential, sine and cosine series)
> For every real $x$,
> $$
> e^x=\sum_{k=0}^\infty\frac{x^k}{k!},\qquad \cos x=\sum_{m=0}^\infty\frac{(-1)^mx^{2m}}{(2m)!},\qquad \sin x=\sum_{m=0}^\infty\frac{(-1)^mx^{2m+1}}{(2m+1)!}.
> $$

**Proof.** First, $r^n/n!\to0$ for every fixed $r\ge0$: choose an integer $N\ge2r$; for $n>N$, each factor $r/k$ with $k>N$ is at most $\tfrac12$, so $\frac{r^n}{n!}=\frac{r^N}{N!}\prod_{k=N+1}^n\frac rk\le\frac{r^N}{N!}2^{-(n-N)}\to0$.

For $\cos$ and $\sin$, every derivative is one of $\pm\sin,\pm\cos$, so $\lvert f^{(n+1)}\rvert\le1$ everywhere; by Theorem A.5(a) with $a=0$, $\lvert R_n(x)\rvert\le\lvert x\rvert^{n+1}/(n+1)!\to0$. The Taylor coefficients at $0$ follow from the cycle $\cos,\,-\sin,\,-\cos,\,\sin$ evaluated at $0$: $1,0,-1,0,\dots$ for $\cos$ and $0,1,0,-1,\dots$ for $\sin$. For $e^x$ every derivative is $e^x$, so on the interval between $0$ and $x$, $\lvert f^{(n+1)}\rvert\le e^{\lvert x\rvert}$ and $\lvert R_n(x)\rvert\le e^{\lvert x\rvert}\lvert x\rvert^{n+1}/(n+1)!\to0$. In each case $p_n(x)\to f(x)$, which is what convergence of the series to $f$ means. ∎

**Numerical check.** For $e^{0.5}$ with $n=3$: $p_3(0.5)=1+0.5+0.125+0.0208\overline{3}=1.6458\overline{3}$, and $e^{0.5}=1.6487213$, so the actual error is $0.0028879$. Lagrange's bound with $M=e^{0.5}$ is $1.6487\times0.5^4/24=0.0042935$, which is indeed larger. For $\cos\pi$ with the degree-6 polynomial (which equals $p_7$, since the $x^7$ coefficient is zero), $p_6(\pi)=-1.21135$, error $0.21135$, and the bound $\pi^8/8!=0.23533$ holds.

**Remarks.**

- *A smooth function need not equal its Taylor series.* Let $\varphi(x)=e^{-1/x^2}$ for $x\ne0$ and $\varphi(0)=0$. By induction, for $x\ne0$, $\varphi^{(k)}(x)=P_k(1/x)\,e^{-1/x^2}$ with $P_k$ a polynomial (differentiate and use the chain rule). Since $u^me^{-u^2}\to0$ as $\lvert u\rvert\to\infty$ for every $m$, the difference quotient $\varphi^{(k)}(h)/h=\frac1hP_k(\frac1h)e^{-1/h^2}\to0$, so $\varphi^{(k+1)}(0)=0$ for every $k$. The Taylor series of $\varphi$ at $0$ is identically zero, yet $\varphi(x)>0$ for $x\ne0$: here $R_n(x)=\varphi(x)$ for every $n$, and Theorem A.5 is not contradicted — the derivatives of $\varphi$ near $0$ grow so fast that the bound does not tend to $0$. Functions equal to their Taylor series near every point are called **analytic**; in complex analysis ([Ch. 18](../part-5-transforms/18-complex-analysis.md)) differentiability in the complex sense forces analyticity, which is one reason that theory is so rigid.
- *Lagrange's form does not hold for complex-valued functions* (it is a mean value statement; see the remark in §A.1). The integral form (b) does hold for complex-valued $f$, by applying it to real and imaginary parts, and gives the bound $\lvert R_n(x)\rvert\le M\lvert x-a\rvert^{n+1}/(n+1)!$ with $M=\max\lvert f^{(n+1)}\rvert$, since $\lvert\int F\rvert\le\int\lvert F\rvert$.
- *Peano form.* If $f$ is $C^{n+1}$ near $a$, part (a) gives $R_n(x)=O(\lvert x-a\rvert^{n+1})$: halving the step reduces the error of a degree-$n$ approximation by about $2^{n+1}$. This is the source of the error orders of the finite-difference formulas in [Ch. 4 §4.7](../part-2-calculus/04-differential-calculus.md), of the trapezoidal and Simpson rules in [Ch. 5 §5.7](../part-2-calculus/05-integral-calculus.md), and of Runge–Kutta methods in [Ch. 22](../part-7-numerical-computational/22-numerical-methods.md).
- *Where it is used.* Linearization and small-signal models ([Ch. 4 §4.8](../part-2-calculus/04-differential-calculus.md)), the second-derivative test and Newton's method in optimization ([Ch. 23 §23.4–23.6](../part-7-numerical-computational/23-optimization.md)) and the two-variable Taylor theorem of [Ch. 6 §6.3](../part-2-calculus/06-multivariable-calculus.md) (obtained by applying Theorem A.5 to $t\mapsto f(\mathbf{a}+t\mathbf{h})$).

```desmos-graph
left=-3; right=3; bottom=-0.2; top=1.2
---
y=e^{-1/x^{2}}|blue
y=0|red|dashed
```

*The function $e^{-1/x^2}$ (blue) is infinitely differentiable, and every derivative vanishes at $0$, so its Maclaurin series is the zero function (red, dashed) — which agrees with it only at the origin.*

## A.5 L'Hôpital's rule (the $0/0$ case)

> [!theorem] Theorem A.6 (L'Hôpital's rule, $0/0$ form)
> Let $-\infty\le a<b\le\infty$ and let $f,g$ be differentiable on $(a,b)$ with $g'(x)\ne0$ for all $x\in(a,b)$. Suppose
> $$
> \lim_{x\to a^+}f(x)=\lim_{x\to a^+}g(x)=0\qquad\text{and}\qquad\lim_{x\to a^+}\frac{f'(x)}{g'(x)}=L,
> $$
> where $L$ is a real number, $+\infty$ or $-\infty$. Then $g(x)\ne0$ on $(a,b)$ and
> $$
> \lim_{x\to a^+}\frac{f(x)}{g(x)}=L .
> $$
> The same holds for left-hand limits $x\to b^-$ (including $b=+\infty$), and hence, by combining the two sides, for two-sided limits $x\to a$ when the hypotheses hold on a punctured neighbourhood of $a$.

**Proof for finite $a$.** Extend $f$ and $g$ to $[a,b)$ by setting $f(a)=g(a)=0$. Because the limits at $a^+$ are $0$, the extended functions are continuous on $[a,x]$ for every $x\in(a,b)$, and they are differentiable on $(a,x)$.

*Step 1: $g$ does not vanish.* If $g(x)=0$ for some $x\in(a,b)$, then $g(a)=g(x)=0$ and Rolle's theorem would give a point of $(a,x)$ where $g'=0$, contrary to hypothesis.

*Step 2: Cauchy's MVT.* Fix $x\in(a,b)$. By Theorem A.1(3) on $[a,x]$ (where $g'\ne0$ on $(a,x)$), there is $\xi_x\in(a,x)$ with

$$
\frac{f(x)}{g(x)}=\frac{f(x)-f(a)}{g(x)-g(a)}=\frac{f'(\xi_x)}{g'(\xi_x)} .
$$

*Step 3: pass to the limit.* Let $V$ be any neighbourhood of $L$ (an interval $(L-\varepsilon,L+\varepsilon)$ if $L$ is finite, or $(K,\infty)$ if $L=+\infty$, or $(-\infty,-K)$ if $L=-\infty$). By hypothesis there is $\eta>0$ with $f'(s)/g'(s)\in V$ for all $s\in(a,a+\eta)$. If $x\in(a,a+\eta)$, then $\xi_x\in(a,x)\subset(a,a+\eta)$, so $f(x)/g(x)=f'(\xi_x)/g'(\xi_x)\in V$. Since $V$ was arbitrary, $f(x)/g(x)\to L$. ∎

**Proof for $a=-\infty$ (and, symmetrically, for $b=+\infty$).** Take $b=+\infty$ for definiteness, with the hypotheses holding on $(M,\infty)$, $M>0$, as $x\to+\infty$. Put $F(t)=f(1/t)$ and $G(t)=g(1/t)$ for $t\in(0,1/M)$. By the chain rule $F'(t)=-f'(1/t)/t^2$ and $G'(t)=-g'(1/t)/t^2\ne0$, and

$$
\lim_{t\to0^+}F(t)=\lim_{x\to\infty}f(x)=0,\quad\lim_{t\to0^+}G(t)=0,\quad\frac{F'(t)}{G'(t)}=\frac{f'(1/t)}{g'(1/t)}\xrightarrow[t\to0^+]{}L .
$$

The finite case at $0^+$ gives $F(t)/G(t)\to L$, i.e. $f(x)/g(x)\to L$ as $x=1/t\to\infty$. ∎

**Example.** $\displaystyle\lim_{x\to0}\frac{1-\cos x}{x^2}$: numerator and denominator tend to $0$, and $g'(x)=2x\ne0$ for $x\ne0$. The ratio of derivatives $\frac{\sin x}{2x}$ is again $0/0$ with $(2x)'=2\ne0$, and $\frac{\cos x}{2}\to\frac12$. Applying Theorem A.6 twice (inner limit first) gives $\frac12$. The Taylor series give the same answer at once: $1-\cos x=\frac{x^2}{2}-\frac{x^4}{24}+\dots$. Similarly $\lim_{x\to0}\frac{\sin x-x}{x^3}=-\frac16$.

**Remarks.**

- *The Ch. 4 sketch versus this proof.* The sketch in [Ch. 4 §4.2](../part-2-calculus/04-differential-calculus.md) assumes $f,g$ are differentiable *at* $a$ with $g'(a)\ne0$. The proof above assumes nothing at $a$ itself, which is what is needed when the rule must be applied repeatedly (after one application, $f'$ and $g'$ are typically $0/0$ at $a$ again) and for limits at infinity.
- *The rule is one-directional.* $f(x)=x^2\sin(1/x)$, $g(x)=x$ at $a=0$: $f/g=x\sin(1/x)\to0$ by the squeeze theorem, but $f'/g'=2x\sin(1/x)-\cos(1/x)$ has no limit. If $\lim f'/g'$ does not exist, the rule says nothing.
- *The hypothesis $g'\ne0$ is not a formality.* It is needed in Step 1 and in Cauchy's theorem. The classical example of Stolz, in the $\infty/\infty$ form, shows what goes wrong: $f(x)=x+\sin x\cos x$, $g(x)=e^{\sin x}f(x)$ as $x\to\infty$. Here $f'=2\cos^2x$ and $g'=e^{\sin x}\cos x\,(f+2\cos x)$, so after cancelling the common factor $\cos x$, $\frac{f'}{g'}=\frac{2e^{-\sin x}\cos x}{f+2\cos x}\to0$; but $f/g=e^{-\sin x}$ oscillates between $e^{-1}$ and $e$. The culprit: $g'$ vanishes at every zero of $\cos x$, arbitrarily far out, and the "cancellation" divided by zero.
- *The $\infty/\infty$ form.* If instead $\lvert g(x)\rvert\to\infty$ as $x\to a^+$ (no condition on $f$), the conclusion again holds. The proof applies Cauchy's MVT on $[x,y]$ with $y$ fixed close to $a$ and $x\to a^+$, writing $\frac{f(x)}{g(x)}=\frac{f(x)-f(y)}{g(x)-g(y)}\cdot\frac{g(x)-g(y)}{g(x)}+\frac{f(y)}{g(x)}$ and noting that the last two factors tend to $1$ and $0$; see Rudin's *Principles*, Thm 5.13.
- *Where it is used.* Limits in [Ch. 4 §4.1–4.2](../part-2-calculus/04-differential-calculus.md) (exponentials dominate powers), initial and final value theorems and the residue at a simple pole $\lim_{s\to p}(s-p)N(s)/D(s)=N(p)/D'(p)$ ([Ch. 14](../part-5-transforms/14-laplace-transforms.md), [Ch. 18](../part-5-transforms/18-complex-analysis.md)), and the value $\operatorname{sinc}(0)=1$ used throughout [Ch. 16](../part-5-transforms/16-fourier-transforms.md).

## A.6 Euler's formula

[Ch. 2 §2.4](../part-1-foundations/02-complex-numbers.md) derives $e^{j\theta}=\cos\theta+j\sin\theta$ in three ways. Here we make two of them rigorous. The first question is what $e^{j\theta}$ *means*: the real exponential is not defined for imaginary arguments, so we must choose a definition and then prove the formula from it.

> [!definition] The complex exponential
> For $z\in\mathbb{C}$,
> $$
> e^z=\exp(z)=\sum_{n=0}^\infty\frac{z^n}{n!}=\lim_{N\to\infty}\sum_{n=0}^N\frac{z^n}{n!}.
> $$
> The series converges for every $z$: its real and imaginary parts are real series whose terms are bounded in modulus by $\lvert z\rvert^n/n!$, and $\sum\lvert z\rvert^n/n!=e^{\lvert z\rvert}$ converges (Corollary A.5), so both are absolutely convergent and hence convergent. For real $z$ this agrees with the real exponential, by Corollary A.5.

Throughout this section write

$$
c(\theta)=\cos\theta+j\sin\theta,\qquad \bar c(\theta)=\cos\theta-j\sin\theta\qquad(\theta\in\mathbb{R}),
$$

so that $c'(\theta)=-\sin\theta+j\cos\theta=j\,c(\theta)$, $\bar c'(\theta)=-\sin\theta-j\cos\theta=-j\,\bar c(\theta)$, and $c(\theta)\,\bar c(\theta)=\cos^2\theta+\sin^2\theta=1$.

> [!theorem] Theorem A.7 (Euler's formula)
> For every real $\theta$,
> $$
> e^{j\theta}=\cos\theta+j\sin\theta .
> $$
> Moreover $\dfrac{d}{d\theta}e^{j\theta}=j\,e^{j\theta}$, and $e^{j\theta}$ is the **only** differentiable function $y:\mathbb{R}\to\mathbb{C}$ with $y'=jy$ and $y(0)=1$.

We give two independent proofs. The first uses the Taylor series of $\cos$ and $\sin$; the second uses only their derivatives.

### Proof 1 — power series

Let $S_N(\theta)=\sum_{n=0}^N\frac{(j\theta)^n}{n!}$. The powers of $j$ cycle with period 4: $j^{2m}=(j^2)^m=(-1)^m$ and $j^{2m+1}=j\,(-1)^m$. Split the partial sum with an odd number of terms, $S_{2N+1}$, into even and odd $n$:

$$
S_{2N+1}(\theta)=\sum_{m=0}^{N}\frac{(-1)^m\theta^{2m}}{(2m)!}+j\sum_{m=0}^{N}\frac{(-1)^m\theta^{2m+1}}{(2m+1)!}.
$$

No rearrangement of an infinite series is involved: this is a regrouping of finitely many terms. The real part is the $N$-th partial sum of the cosine series and the imaginary part is the $N$-th partial sum of the sine series, so by Corollary A.5,

$$
\operatorname{Re}\{S_{2N+1}(\theta)\}\to\cos\theta,\qquad\operatorname{Im}\{S_{2N+1}(\theta)\}\to\sin\theta\qquad(N\to\infty).
$$

A sequence of complex numbers converges if and only if its real and imaginary parts converge, so $S_{2N+1}(\theta)\to\cos\theta+j\sin\theta$. But the whole sequence $S_N(\theta)$ converges to $e^{j\theta}$ by definition, and every subsequence of a convergent sequence has the same limit. Hence $e^{j\theta}=\cos\theta+j\sin\theta$. ∎

### Proof 2 — the differential equation and uniqueness

This proof has two steps: the series satisfies the initial-value problem $y'=jy$, $y(0)=1$, and that problem has only one solution.

> [!theorem] Lemma A.7a (Uniqueness)
> Let $I$ be an interval containing $0$ and $y:I\to\mathbb{C}$ differentiable with $y'(\theta)=j\,y(\theta)$ on $I$. Then $y(\theta)=y(0)\,c(\theta)$ for all $\theta\in I$.

**Proof.** Let $w(\theta)=y(\theta)\,\bar c(\theta)$. By the product rule (valid for complex-valued functions),

$$
w'=y'\,\bar c+y\,\bar c'=(jy)\,\bar c+y\,(-j\bar c)=0 .
$$

By Corollary A.1(1) for complex-valued functions, $w$ is constant on $I$, equal to $w(0)=y(0)\bar c(0)=y(0)$. Multiplying $y\bar c=y(0)$ by $c$ and using $c\bar c=1$ gives $y=y(0)\,c$. ∎

> [!theorem] Lemma A.7b (Derivative of the series)
> $E(\theta)=\sum_{n=0}^\infty\frac{(j\theta)^n}{n!}$ is differentiable on $\mathbb{R}$ with $E'(\theta)=j\,E(\theta)$.

**Proof.** Fix $\theta$ and let $0<\lvert h\rvert\le1$. For each $n\ge2$, Theorem A.5(a) applied to $t\mapsto t^n$ at $a=\theta$ with $n=1$ gives

$$
(\theta+h)^n-\theta^n-n\theta^{n-1}h=\frac{n(n-1)\gamma^{n-2}}{2}h^2
$$

for some $\gamma$ between $\theta$ and $\theta+h$, so that $\lvert\gamma\rvert\le\lvert\theta\rvert+1=:r$. For $n=0,1$ the left side is $0$. The series $\sum_n\frac{j^n n\theta^{n-1}}{n!}=j\sum_{n\ge1}\frac{(j\theta)^{n-1}}{(n-1)!}=jE(\theta)$ converges, so we may combine three convergent series term by term:

$$
E(\theta+h)-E(\theta)-h\,jE(\theta)=\sum_{n=2}^\infty\frac{j^n}{n!}\Big[(\theta+h)^n-\theta^n-n\theta^{n-1}h\Big].
$$

Taking moduli and using $\lvert j^n\rvert=1$,

$$
\Big\lvert E(\theta+h)-E(\theta)-h\,jE(\theta)\Big\rvert\le\frac{h^2}{2}\sum_{n=2}^\infty\frac{n(n-1)}{n!}r^{n-2}=\frac{h^2}{2}\sum_{m=0}^\infty\frac{r^m}{m!}=\frac{h^2}{2}e^{r}.
$$

Dividing by $\lvert h\rvert$: $\left\lvert\frac{E(\theta+h)-E(\theta)}{h}-jE(\theta)\right\rvert\le\frac{\lvert h\rvert}{2}e^{\lvert\theta\rvert+1}\to0$ as $h\to0$. ∎

(The identity $\sum_{n\ge2}\frac{n(n-1)}{n!}r^{n-2}=e^r$ used in the last step was also checked numerically, e.g. both sides equal $9.97418$ at $r=2.3$.)

**Conclusion of Proof 2.** By Lemma A.7b, $y=E$ satisfies $y'=jy$ on $\mathbb{R}$ with $y(0)=1$. By Lemma A.7a, $E(\theta)=c(\theta)=\cos\theta+j\sin\theta$. The uniqueness statement of the theorem is Lemma A.7a with $y(0)=1$. ∎

```desmos-graph
left=-7; right=7; bottom=-2; top=2
---
y=\cos(x)|black
y=1-\frac{x^{2}}{2}|red|dashed
y=1-\frac{x^{2}}{2}+\frac{x^{4}}{24}|orange|dashed
y=1-\frac{x^{2}}{2}+\frac{x^{4}}{24}-\frac{x^{6}}{720}+\frac{x^{8}}{40320}|green
y=1-\frac{x^{2}}{2}+\frac{x^{4}}{24}-\frac{x^{6}}{720}+\frac{x^{8}}{40320}-\frac{x^{10}}{3628800}+\frac{x^{12}}{479001600}|blue
```

*The real part of the partial sums of $e^{j\theta}$ (degrees 2, 4, 8, 12) against $\cos\theta$ (black). Each polynomial hugs the cosine over a wider interval; by Corollary A.5 the agreement extends to any bounded interval as the degree grows.*

> [!theorem] Corollary A.7 (Consequences of Euler's formula)
> For all real $\alpha,\beta,\theta$:
> 1. $\lvert e^{j\theta}\rvert=1$ and $\overline{e^{j\theta}}=e^{-j\theta}=1/e^{j\theta}$.
> 2. $e^{j(\alpha+\beta)}=e^{j\alpha}e^{j\beta}$; taking real and imaginary parts gives $\cos(\alpha+\beta)=\cos\alpha\cos\beta-\sin\alpha\sin\beta$ and $\sin(\alpha+\beta)=\sin\alpha\cos\beta+\cos\alpha\sin\beta$.
> 3. $\cos\theta=\dfrac{e^{j\theta}+e^{-j\theta}}2$ and $\sin\theta=\dfrac{e^{j\theta}-e^{-j\theta}}{2j}$.
> 4. $e^{j\theta}=1$ if and only if $\theta=2\pi m$ for some integer $m$; hence $e^{j\theta}$ is periodic with fundamental period $2\pi$.
> 5. For real $\omega$, $\dfrac{d}{dt}e^{j\omega t}=j\omega\,e^{j\omega t}$.

**Proof.** (1) $\lvert e^{j\theta}\rvert^2=\cos^2\theta+\sin^2\theta=1$; $\overline{c(\theta)}=\cos\theta-j\sin\theta=\cos(-\theta)+j\sin(-\theta)=c(-\theta)$ because $\cos$ is even and $\sin$ odd; and $c(\theta)c(-\theta)=c(\theta)\bar c(\theta)=1$.
(2) Fix $\alpha$ and let $y(\theta)=e^{j(\alpha+\theta)}=c(\alpha+\theta)$. By the chain rule $y'(\theta)=c'(\alpha+\theta)=j\,y(\theta)$, and $y(0)=e^{j\alpha}$. Lemma A.7a gives $y(\theta)=e^{j\alpha}e^{j\theta}$. Expanding $c(\alpha)c(\beta)$ and comparing real and imaginary parts gives the compound-angle formulas — so these formulas are a *consequence* of Euler's formula, not an input to it (Proof 2 used only $\sin'=\cos$, $\cos'=-\sin$ and $\cos^2+\sin^2=1$).
(3) Add and subtract $e^{j\theta}=\cos\theta+j\sin\theta$ and $e^{-j\theta}=\cos\theta-j\sin\theta$.
(4) $e^{j\theta}=1$ means $\cos\theta=1$ and $\sin\theta=0$; and $\cos\theta=1$ exactly when $\theta\in2\pi\mathbb{Z}$, which is the statement that the point at arc length $\theta$ on the unit circle returns to $(1,0)$ only after whole turns (the unit-circle definition of $\cos$, [Ch. 3 §3.1](../part-1-foundations/03-trigonometry-and-periodic-functions.md)); conversely $\theta=2\pi m$ gives $\cos\theta=1$, $\sin\theta=0$. Then $e^{j(\theta+T)}=e^{j\theta}$ for all $\theta$ iff $e^{jT}=1$ (by (2)) iff $T\in2\pi\mathbb{Z}$.
(5) Chain rule with Theorem A.7: $\frac{d}{dt}c(\omega t)=\omega\,c'(\omega t)=j\omega\,c(\omega t)$. ∎

**Remarks.**

- *Why two proofs.* Proof 1 explains where the formula "comes from" — the even and odd parts of the exponential series are the cosine and sine series — but it leans on Taylor's theorem for $\cos$ and $\sin$. Proof 2 is the more robust argument: it shows that *any* reasonable definition of $e^{j\theta}$ — anything differentiable that satisfies $y'=jy$, $y(0)=1$ — is forced to equal $\cos\theta+j\sin\theta$. It is the simplest instance of the uniqueness theorem for linear ODEs ([Ch. 11 §11.7](../part-4-differential-equations/11-ordinary-differential-equations.md)), and the same "multiply by the inverse solution and show the derivative vanishes" trick proves that $e^{\mathbf{A}t}$ is the unique solution of $\dot{\mathbf{x}}=\mathbf{A}\mathbf{x}$ ([Ch. 11 §11.9](../part-4-differential-equations/11-ordinary-differential-equations.md)).
- *Radians are essential.* The proofs use $\sin'=\cos$, which holds only when angles are measured in radians ([Ch. 3 §3.1](../part-1-foundations/03-trigonometry-and-periodic-functions.md)). In degrees, $\frac{d}{d\theta}\sin\theta^\circ=\frac{\pi}{180}\cos\theta^\circ$, and $e^{j\theta}$ would have to be replaced by $e^{j\pi\theta/180}$.
- *The general exponential.* For $z=\sigma+j\omega$, the identity $e^z=e^\sigma e^{j\omega}=e^\sigma(\cos\omega+j\sin\omega)$ requires the addition law $e^{z+w}=e^ze^w$ for complex $z,w$, which follows from the Cauchy product of two absolutely convergent series (or from the uniqueness argument in the complex-variable setting of [Ch. 18](../part-5-transforms/18-complex-analysis.md)). It is the basis of the $s$-plane picture of $e^{st}$: $\sigma$ sets growth or decay and $\omega$ sets oscillation.
- *Where it is used.* Phasors and impedance ([Ch. 2 §2.8–2.9](../part-1-foundations/02-complex-numbers.md)); all trigonometric identities by the "master method" ([Ch. 3 §3.6](../part-1-foundations/03-trigonometry-and-periodic-functions.md)); complex Fourier series and transforms ([Ch. 15](../part-5-transforms/15-fourier-series.md), [Ch. 16](../part-5-transforms/16-fourier-transforms.md)); characteristic-root solutions of linear ODEs ([Ch. 12](../part-4-differential-equations/12-linear-differential-equations.md)).

## A.7 De Moivre's theorem and the $n$-th roots of unity

We use the polar form of [Ch. 2 §2.1](../part-1-foundations/02-complex-numbers.md): every $z\ne0$ can be written $z=r\,e^{j\varphi}$ with $r=\lvert z\rvert>0$ and $\varphi$ real (an argument of $z$), and $\varphi$ is determined up to an additive multiple of $2\pi$ (Corollary A.7(4)).

> [!theorem] Theorem A.8 (De Moivre)
> For every integer $n$ (positive, zero or negative) and every real $\theta$,
> $$
> (\cos\theta+j\sin\theta)^n=\cos n\theta+j\sin n\theta,\qquad\text{i.e.}\qquad \big(e^{j\theta}\big)^n=e^{jn\theta},
> $$
> and for $r>0$, $(re^{j\theta})^n=r^ne^{jn\theta}$.

**Proof.** *$n\ge0$, by induction.* For $n=0$ both sides equal $1$. If $(e^{j\theta})^n=e^{jn\theta}$, then by Corollary A.7(2),

$$
\big(e^{j\theta}\big)^{n+1}=e^{jn\theta}\,e^{j\theta}=e^{j(n\theta+\theta)}=e^{j(n+1)\theta}.
$$

*$n=-m<0$.* By Corollary A.7(1), $(e^{j\theta})^{-1}=e^{-j\theta}$, so $(e^{j\theta})^{-m}=\big((e^{j\theta})^{-1}\big)^m=(e^{j(-\theta)})^m=e^{-jm\theta}$ by the case already proved. The last statement follows because real and complex powers of a product factor: $(r\,e^{j\theta})^n=r^n(e^{j\theta})^n$ for integer $n$ (commutativity of multiplication). ∎

> [!definition] $n$-th roots and roots of unity
> For an integer $n\ge1$ and $z\in\mathbb{C}$, an **$n$-th root** of $z$ is any $w\in\mathbb{C}$ with $w^n=z$. The $n$-th roots of $1$ are the **$n$-th roots of unity**. Write $W=W_n=e^{j2\pi/n}$.

> [!theorem] Theorem A.9 ($n$-th roots)
> Let $n\ge1$ and $z=re^{j\varphi}\ne0$. Then $w^n=z$ has **exactly $n$** solutions,
> $$
> w_k=r^{1/n}\,e^{j(\varphi+2\pi k)/n},\qquad k=0,1,\dots,n-1,
> $$
> where $r^{1/n}$ is the positive real $n$-th root. They are $w_k=w_0W^k$, equally spaced on the circle of radius $r^{1/n}$. In particular the $n$-th roots of unity are $1,W,W^2,\dots,W^{n-1}$, and:
> 1. $w^n-1=\displaystyle\prod_{k=0}^{n-1}\big(w-W^k\big)$ for all $w\in\mathbb{C}$;
> 2. for every integer $m$,
> $$
> \sum_{k=0}^{n-1}W^{mk}=\begin{cases}n,& n\mid m,\\ 0,&\text{otherwise.}\end{cases}
> $$

**Proof.** *Each $w_k$ is a root.* By Theorem A.8, $w_k^n=r\,e^{j(\varphi+2\pi k)}=r\,e^{j\varphi}e^{j2\pi k}=z$, using Corollary A.7(2) and $e^{j2\pi k}=1$.

*Every root is some $w_k$.* Let $w^n=z$. Then $w\ne0$ (else $z=0$), so $w=\rho e^{j\psi}$ with $\rho>0$. By Theorem A.8, $\rho^ne^{jn\psi}=re^{j\varphi}$. Taking moduli (Corollary A.7(1)) gives $\rho^n=r$; since $t\mapsto t^n$ is strictly increasing on $(0,\infty)$ (its derivative $nt^{n-1}$ is positive, Corollary A.1(3)), $\rho=r^{1/n}$. Dividing out, $e^{j(n\psi-\varphi)}=1$, so by Corollary A.7(4), $n\psi-\varphi=2\pi m$ for some integer $m$, i.e. $\psi=(\varphi+2\pi m)/n$. By the division algorithm $m=qn+k$ with $q\in\mathbb{Z}$ and $0\le k\le n-1$; then $e^{j\psi}=e^{j(\varphi+2\pi k)/n}\,e^{j2\pi q}=e^{j(\varphi+2\pi k)/n}$, so $w=w_k$.

*The $w_k$ are distinct.* If $w_k=w_l$ with $0\le l<k\le n-1$, dividing gives $e^{j2\pi(k-l)/n}=1$, so $(k-l)/n$ is an integer by Corollary A.7(4). But $0<k-l<n$, a contradiction. Hence there are exactly $n$ roots. The form $w_k=w_0W^k$ is Corollary A.7(2).

*(1) Factorisation.* Let $p(w)=w^n-1$. If $p(\alpha)=0$, polynomial division gives $p(w)=(w-\alpha)q(w)$ with $q$ monic of degree $n-1$ (the factor theorem). Apply this with $\alpha=W^0=1$; each remaining $W^k$, $k\ge1$, satisfies $0=p(W^k)=(W^k-1)q(W^k)$ with $W^k-1\ne0$ (distinctness), so $W^k$ is a root of $q$. Repeating $n$ times leaves a monic polynomial of degree $0$, i.e. the constant $1$.

*(2) Sums.* If $n\mid m$, then $W^m=e^{j2\pi(m/n)}=1$ and every term is $1$: the sum is $n$. Otherwise $\rho=W^m\ne1$ (Corollary A.7(4)) and $\rho^n=(W^n)^m=1$. The finite geometric sum gives $\sum_{k=0}^{n-1}\rho^k=\frac{1-\rho^n}{1-\rho}=0$. ∎

**Example.** The cube roots of $z=8j$: $r=8$, $\varphi=\pi/2$, so $w_k=2e^{j(\pi/6+2\pi k/3)}$, i.e. $w_0=\sqrt3+j$, $w_1=-\sqrt3+j$, $w_2=-2j$. Each cubes to $8j$ (checked numerically), and their sum is $0$, as part 2 with $m=1$ predicts for $w_0(1+W+W^2)$.

```desmos-graph
left=-2.6; right=2.6; bottom=-2.4; top=2.4
---
x^{2}+y^{2}=1|grey|dashed
x^{2}+y^{2}=4|grey|dotted
(1,0)|blue|label:W^0
(0.5,0.866)|blue|label:W
(-0.5,0.866)|blue|label:W^2
(-1,0)|blue|label:W^3
(-0.5,-0.866)|blue|label:W^4
(0.5,-0.866)|blue|label:W^5
(1.732,1)|red|label:w_0
(-1.732,1)|red|label:w_1
(0,-2)|red|label:w_2
```

*The sixth roots of unity (blue) form a regular hexagon on the unit circle; the cube roots of $8j$ (red) form an equilateral triangle on the circle of radius $8^{1/3}=2$, rotated to angle $\varphi/n=30^\circ$.*

**Remarks.**

- *Why exactly $n$ and not infinitely many.* The infinitely many arguments $\varphi+2\pi m$ of $z$ give infinitely many candidate angles $(\varphi+2\pi m)/n$, but they repeat with period $n$ in $m$. This is the precise sense in which $z^{1/n}$ is "multi-valued" ([Ch. 2 §2.6](../part-1-foundations/02-complex-numbers.md), [Ch. 18](../part-5-transforms/18-complex-analysis.md)).
- *De Moivre is for integer exponents.* Applied naively to $n=\tfrac12$ it would give a contradiction: $e^{j0}=e^{j2\pi}=1$, yet "$\big(e^{j2\pi}\big)^{1/2}=e^{j\pi}=-1$" while $\big(e^{j0}\big)^{1/2}=1$. Likewise $(z^a)^b=z^{ab}$ fails for non-integer exponents: $\big((-1)^2\big)^{1/2}=1$ but $\big((-1)^{1/2}\big)^2=-1$ for either choice of $(-1)^{1/2}=\pm j$.
- *The hypothesis $z\ne0$.* For $z=0$ the only root is $w=0$, with multiplicity $n$; the polar argument is undefined.
- *Where it is used.* Part 2 is the orthogonality relation of the DFT, $\sum_{k=0}^{N-1}e^{j2\pi(m-l)k/N}=N\delta_{ml}$ for $0\le m,l<N$ ([Ch. 15 §15.7](../part-5-transforms/15-fourier-series.md), [Ch. 24](../part-7-numerical-computational/24-computational-mathematics.md)), and with $n=3$ it is the identity $1+a+a^2=0$, $a=e^{j2\pi/3}$, which makes the phasors of a balanced three-phase set sum to zero ([Ch. 3 §3.3](../part-1-foundations/03-trigonometry-and-periodic-functions.md)). Part 1 underlies root-finding and factoring in [Ch. 2 §2.6](../part-1-foundations/02-complex-numbers.md), the poles of Butterworth filters (equally spaced on a circle in the $s$-plane), and the $N$-point FFT's recursive splitting of $W_N$ into $W_{N/2}$.

## A.8 RMS of a sinusoid and of a sum of sinusoids of different frequencies

The statement "powers of components at different frequencies add" is used constantly — for harmonic distortion, for noise, for Parseval's theorem. Here it is proved in two settings: exactly over a common period when the frequencies are harmonics of one fundamental, and as a long-time average when they are not (for example $\cos t+\cos\sqrt2\,t$, which is not periodic at all).

> [!definition] Mean, mean-square and RMS
> For $x$ continuous on $[t_0,t_0+T]$, $T>0$, the **mean** and **RMS value over the window** are
> $$
> \bar x_{[t_0,T]}=\frac1T\int_{t_0}^{t_0+T}x(t)\,dt,\qquad X_{\text{rms},[t_0,T]}=\sqrt{\frac1T\int_{t_0}^{t_0+T}x(t)^2\,dt}.
> $$
> For a $T$-periodic $x$ the window is one period. For an arbitrary bounded $x$, the **long-time mean-square value** (average power in $1\ \Omega$) is $\langle x^2\rangle=\lim_{T\to\infty}\frac1T\int_{t_0}^{t_0+T}x^2\,dt$, when the limit exists.

> [!theorem] Lemma A.8 (Averages over a period do not depend on the starting point; averages of sinusoids)
> 1. If $x$ is continuous and $T$-periodic, then $\int_{t_0}^{t_0+T}x\,dt$ is the same for every $t_0$.
> 2. For $\Omega\ne0$ and any real $\psi$, $t_0$, $T>0$,
> $$
> \left\lvert\int_{t_0}^{t_0+T}\cos(\Omega t+\psi)\,dt\right\rvert\le\frac{2}{\lvert\Omega\rvert},
> $$
> and the integral is exactly $0$ when $T$ is an integer multiple of $2\pi/\lvert\Omega\rvert$.

**Proof.** (1) Let $I(t_0)=\int_{t_0}^{t_0+T}x\,dt=\int_0^{t_0+T}x\,dt-\int_0^{t_0}x\,dt$. By FTC Part 1 (Theorem A.2) and the chain rule, $I'(t_0)=x(t_0+T)-x(t_0)=0$, so $I$ is constant by Corollary A.1(1).
(2) An antiderivative of $\cos(\Omega t+\psi)$ is $\sin(\Omega t+\psi)/\Omega$, so by FTC Part 2 the integral equals $\big[\sin(\Omega(t_0+T)+\psi)-\sin(\Omega t_0+\psi)\big]/\Omega$, whose modulus is at most $2/\lvert\Omega\rvert$. If $T=2\pi q/\lvert\Omega\rvert$ with $q$ an integer, the two sines are equal because $\sin$ has period $2\pi$. ∎

> [!theorem] Theorem A.10 (RMS of sinusoids)
> **(a) Single sinusoid.** Let $x(t)=C+A\cos(\omega t+\phi)$ with $\omega>0$ and $A,C,\phi$ real, and $T=2\pi/\omega$. Over any window of length $T$,
> $$
> \bar x=C,\qquad X_{\text{rms}}=\sqrt{C^2+\tfrac12A^2};\qquad\text{for }C=0:\ X_{\text{rms}}=\frac{\lvert A\rvert}{\sqrt2}.
> $$
> **(b) Different frequencies, long-time average.** Let
> $$
> x(t)=A_0+\sum_{k=1}^{N}A_k\cos(\omega_kt+\phi_k),
> $$
> with real $A_k,\phi_k$ and **distinct** positive frequencies $\omega_1,\dots,\omega_N$ (no commensurability assumed). Then, uniformly in $t_0$,
> $$
> \left\lvert\frac1T\int_{t_0}^{t_0+T}x^2\,dt-P\right\rvert\le\frac BT,\qquad P=A_0^2+\sum_{k=1}^N\frac{A_k^2}{2},
> $$
> with the constant
> $$
> B=\sum_{k=1}^N\left(\frac{4\lvert A_0A_k\rvert}{\omega_k}+\frac{A_k^2}{2\omega_k}\right)+\sum_{1\le k<l\le N}2\lvert A_kA_l\rvert\left(\frac{1}{\lvert\omega_k-\omega_l\rvert}+\frac{1}{\omega_k+\omega_l}\right).
> $$
> Hence $\langle x^2\rangle=P$ and the long-time RMS value is $\sqrt{A_0^2+\sum_kA_k^2/2}$, independent of the phases $\phi_k$.
>
> **(c) Harmonics.** If $\omega_k=n_k\omega_0$ for distinct positive integers $n_k$, then $x$ is periodic with period $T_0=2\pi/\omega_0$, and the mean-square value over **any** window of length $T_0$ equals $P$ exactly.

**Proof.** Write $\theta_k=\omega_kt+\phi_k$. Expanding the square of a finite sum,

$$
x^2=A_0^2+2A_0\sum_kA_k\cos\theta_k+\sum_kA_k^2\cos^2\theta_k+2\sum_{k<l}A_kA_l\cos\theta_k\cos\theta_l .
$$

Use the power-reduction and product-to-sum identities (both follow from Corollary A.7(2)–(3)):

$$
\cos^2\theta_k=\tfrac12+\tfrac12\cos2\theta_k,\qquad 2\cos\theta_k\cos\theta_l=\cos(\theta_k-\theta_l)+\cos(\theta_k+\theta_l).
$$

So

$$
x^2=\underbrace{A_0^2+\sum_k\frac{A_k^2}2}_{P}+\;2A_0\sum_kA_k\cos\theta_k+\sum_k\frac{A_k^2}2\cos2\theta_k+\sum_{k<l}A_kA_l\Big[\cos(\theta_k-\theta_l)+\cos(\theta_k+\theta_l)\Big].
$$

Every term after $P$ is a constant times $\cos(\Omega t+\psi)$ with angular frequency $\Omega$ equal to one of $\omega_k$, $2\omega_k$, $\omega_k-\omega_l$ or $\omega_k+\omega_l$ ($k<l$). **None of these is zero**: the $\omega_k$ are positive, and $\omega_k-\omega_l\ne0$ because the frequencies are distinct. This is the only place the hypothesis is used.

*(b)* Integrate over $[t_0,t_0+T]$, divide by $T$ and apply Lemma A.8(2) term by term with the triangle inequality: the terms contribute at most $\frac1T\cdot2\lvert A_0A_k\rvert\frac2{\omega_k}$, $\frac1T\cdot\frac{A_k^2}{2}\cdot\frac{2}{2\omega_k}$, and $\frac1T\lvert A_kA_l\rvert\big(\frac2{\lvert\omega_k-\omega_l\rvert}+\frac2{\omega_k+\omega_l}\big)$. Their sum is $B/T$. Letting $T\to\infty$ gives $\langle x^2\rangle=P$.

*(a)* is the case $N=1$, $\omega_1=\omega$, $A_0=C$ of the expansion: the oscillating terms have frequencies $\omega$ and $2\omega$, and $T=2\pi/\omega$ is a multiple of both $2\pi/\omega$ and $2\pi/(2\omega)$, so by Lemma A.8(2) their integrals over the window vanish exactly; the mean-square is $C^2+A^2/2$. The mean is $C$ by the same lemma applied to $x$ itself.

*(c)* Each $\Omega$ above is now a non-zero integer multiple of $\omega_0$: $n_k\omega_0$, $2n_k\omega_0$, $(n_k-n_l)\omega_0$ or $(n_k+n_l)\omega_0$. So $T_0$ is an integer multiple of every $2\pi/\lvert\Omega\rvert$, every oscillating term integrates to exactly $0$ over any window of length $T_0$ (Lemma A.8(2)), and the mean-square is exactly $P$. ∎

**Numerical checks.**

- [Ch. 3 Worked Example 3.29](../part-1-foundations/03-trigonometry-and-periodic-functions.md): $v=10+8\cos\omega t+6\cos(3\omega t+40^\circ)$ gives $P=100+32+18=150\ \text{V}^2$, $V_{\text{rms}}=12.2474$ V; averaging $v^2$ over one period at 200 000 points gives $12.247449$ V.
- The non-periodic signal $x=\cos t+\cos\sqrt2\,t$ has $P=1$ and, from the formula, $B=\tfrac12+\tfrac1{2\sqrt2}+2\big(\tfrac1{\sqrt2-1}+\tfrac1{\sqrt2+1}\big)=6.51$. Numerical integration gives windowed mean-squares over $[0,T]$ of $0.7848,\ 0.98667,\ 0.99969,\ 1.00028$ for $T=10,100,1000,10\,000$; the deviations $0.215,\ 0.0133,\ 0.0003,\ 0.0003$ respect the bound $B/T=0.65,\ 0.065,\ 0.0065,\ 0.00065$.

```desmos-graph
left=0; right=60; bottom=0; top=2.2
---
y=\frac{1}{x}\int_{0}^{x}(\cos(t)+\cos(1.41421t))^{2}dt|blue|0.2<x<60
y=1|red|dashed
y=1+\frac{6.51}{x}|grey|dotted|3<x<60
y=1-\frac{6.51}{x}|grey|dotted|3<x<60
```

*Running mean-square of $\cos t+\cos\sqrt2\,t$ (blue). It never settles over a "period" (there is none) but converges to $P=\tfrac12+\tfrac12=1$ (red), staying inside the envelope $1\pm B/T$ of Theorem A.10(b) (grey).*

**Remarks.**

- *Distinct frequencies are essential.* For $x=\cos\omega t+\cos(\omega t+\pi)\equiv0$ the formula would predict $\langle x^2\rangle=\tfrac12+\tfrac12=1$, but the true value is $0$; for $\cos\omega t+\cos\omega t$ it would predict $1$ against the true $2$. The cross term $\cos(\theta_k-\theta_l)$ then has frequency $0$ and does not average out. Components at the *same* frequency must first be combined with phasors: $A\cos(\omega t+\alpha)+B\cos(\omega t+\beta)=\operatorname{Re}\{(Ae^{j\alpha}+Be^{j\beta})e^{j\omega t}\}$ has mean-square $\tfrac12\lvert Ae^{j\alpha}+Be^{j\beta}\rvert^2$ ([Ch. 3 §3.3](../part-1-foundations/03-trigonometry-and-periodic-functions.md)).
- *The $1/T$ convergence rate matters for instruments.* A true-RMS meter averaging over a window $T$ has a reading error of order $B/T$; the constant $B$ is dominated by the closest pair of frequencies, $1/\lvert\omega_k-\omega_l\rvert$. Two tones $1$ Hz apart need an averaging time of many seconds — the beat period ([Ch. 3 §3.7](../part-1-foundations/03-trigonometry-and-periodic-functions.md)) — before their powers can be read as adding.
- *Power adds, amplitude does not.* The peak value of $x$ depends strongly on the phases; only the mean-square value is phase-independent.
- *Extension.* The same proof with $N\to\infty$ (under $\sum A_k^2<\infty$ and an argument for exchanging limit and sum) gives Parseval's theorem for Fourier series ([Ch. 15 §15.8](../part-5-transforms/15-fourier-series.md)); the statistical analogue — the variances of uncorrelated random signals add — is in [Ch. 21](../part-6-probability-statistics/21-random-processes.md).
- *Where it is used.* [Ch. 3 §3.2 and §3.7](../part-1-foundations/03-trigonometry-and-periodic-functions.md), [Ch. 5 §5.9](../part-2-calculus/05-integral-calculus.md) (the DC-plus-AC rule $V_{\text{rms}}^2=V_{\text{dc}}^2+V_{\text{ac,rms}}^2$ is the case of a constant plus a zero-mean periodic part) and total harmonic distortion ([Ch. 15 §15.8](../part-5-transforms/15-fourier-series.md)).

## A.9 Equality of mixed partial derivatives (Clairaut–Schwarz)

Notation follows [Ch. 6 §6.2](../part-2-calculus/06-multivariable-calculus.md): $f_{xy}=\dfrac{\partial}{\partial y}\Big(\dfrac{\partial f}{\partial x}\Big)$ (first $x$, then $y$) and $f_{yx}=\dfrac{\partial}{\partial x}\Big(\dfrac{\partial f}{\partial y}\Big)$.

> [!theorem] Theorem A.11 (Clairaut, Schwarz)
> Let $U\subset\mathbb{R}^2$ be open, $(a,b)\in U$, and $f:U\to\mathbb{R}$.
>
> **(a) Clairaut.** If $f_x$, $f_y$, $f_{xy}$ and $f_{yx}$ exist on $U$ and $f_{xy}$, $f_{yx}$ are continuous at $(a,b)$, then $f_{xy}(a,b)=f_{yx}(a,b)$.
>
> **(b) Schwarz (sharper form).** If $f_x$, $f_y$ and $f_{xy}$ exist on $U$ and $f_{xy}$ is continuous at $(a,b)$, then $f_{yx}(a,b)$ exists and equals $f_{xy}(a,b)$.
>
> In particular, if $f\in C^2(U)$, then $f_{xy}=f_{yx}$ throughout $U$.

**Proof.** Choose $\delta>0$ such that the closed square $\lvert x-a\rvert\le\delta$, $\lvert y-b\rvert\le\delta$ lies in $U$. For $0<\lvert h\rvert,\lvert k\rvert<\delta$ form the **second difference**

$$
\Delta(h,k)=f(a+h,b+k)-f(a+h,b)-f(a,b+k)+f(a,b),
$$

the alternating sum of $f$ over the corners of a small rectangle.

*Step 1: $\Delta$ via $f_{xy}$.* Let $\varphi(x)=f(x,b+k)-f(x,b)$ for $x$ between $a$ and $a+h$. Then $\Delta=\varphi(a+h)-\varphi(a)$. Since $f_x$ exists on $U$, $\varphi$ is differentiable with $\varphi'(x)=f_x(x,b+k)-f_x(x,b)$ (and therefore continuous). By the MVT there is $\theta_1\in(0,1)$ with

$$
\Delta=h\,\varphi'(a+\theta_1h)=h\Big[f_x(a+\theta_1h,\,b+k)-f_x(a+\theta_1h,\,b)\Big].
$$

Now fix $x_1=a+\theta_1h$ and consider $y\mapsto f_x(x_1,y)$; it is differentiable with derivative $f_{xy}(x_1,y)$ because $f_{xy}$ exists on $U$. By the MVT again, there is $\theta_2\in(0,1)$ with

$$
\Delta(h,k)=hk\,f_{xy}(a+\theta_1h,\;b+\theta_2k).
$$

*Step 2 (for (a)): $\Delta$ via $f_{yx}$.* Grouping the same four terms the other way, $\Delta=\psi(b+k)-\psi(b)$ with $\psi(y)=f(a+h,y)-f(a,y)$. The same two applications of the MVT, first in $y$ (using $f_y$) and then in $x$ (using $f_{yx}$), give $\theta_3,\theta_4\in(0,1)$ with

$$
\Delta(h,k)=hk\,f_{yx}(a+\theta_4h,\;b+\theta_3k).
$$

*Conclusion of (a).* Comparing, $f_{xy}(a+\theta_1h,b+\theta_2k)=f_{yx}(a+\theta_4h,b+\theta_3k)$ for all small non-zero $h,k$. Put $k=h$ and let $h\to0$: both points tend to $(a,b)$, and continuity of $f_{xy}$ and $f_{yx}$ at $(a,b)$ gives $f_{xy}(a,b)=f_{yx}(a,b)$. ∎

*Proof of (b).* Only Step 1 is used. Let $\varepsilon>0$. By continuity of $f_{xy}$ at $(a,b)$ choose $\delta'\in(0,\delta]$ with $\lvert f_{xy}(p,q)-f_{xy}(a,b)\rvert<\varepsilon$ whenever $\lvert p-a\rvert<\delta'$, $\lvert q-b\rvert<\delta'$. The point in Step 1 satisfies this, so

$$
\left\lvert\frac{\Delta(h,k)}{hk}-f_{xy}(a,b)\right\rvert<\varepsilon\qquad\text{for all }0<\lvert h\rvert,\lvert k\rvert<\delta'.
$$

Fix $h$ and let $k\to0$. Since $f_y$ exists at $(a+h,b)$ and at $(a,b)$,

$$
\frac{\Delta(h,k)}{k}=\frac{f(a+h,b+k)-f(a+h,b)}{k}-\frac{f(a,b+k)-f(a,b)}{k}\;\longrightarrow\;f_y(a+h,b)-f_y(a,b).
$$

A non-strict inequality survives the limit, so

$$
\left\lvert\frac{f_y(a+h,b)-f_y(a,b)}{h}-f_{xy}(a,b)\right\rvert\le\varepsilon\qquad\text{for all }0<\lvert h\rvert<\delta'.
$$

This says precisely that the difference quotient defining $f_{yx}(a,b)=\frac{\partial}{\partial x}f_y\big|_{(a,b)}$ converges to $f_{xy}(a,b)$. ∎

**The counterexample.** Let

$$
f(x,y)=\begin{cases}\dfrac{xy\,(x^2-y^2)}{x^2+y^2},&(x,y)\ne(0,0),\\[2mm]0,&(x,y)=(0,0).\end{cases}
$$

In polar coordinates $f=\tfrac14r^2\sin4\phi$, so $f$ is continuous and small near the origin. For $(x,y)\ne(0,0)$ the quotient rule gives (checked with a computer algebra system)

$$
f_x=\frac{y\,(x^4+4x^2y^2-y^4)}{(x^2+y^2)^2},\qquad f_y=\frac{x\,(x^4-4x^2y^2-y^4)}{(x^2+y^2)^2}.
$$

At the origin, $f(x,0)=0$ and $f(0,y)=0$ give $f_x(0,0)=f_y(0,0)=0$. Since $\lvert x^4\pm4x^2y^2-y^4\rvert\le2(x^2+y^2)^2$, we have $\lvert f_x\rvert\le2\lvert y\rvert$ and $\lvert f_y\rvert\le2\lvert x\rvert$, so $f_x,f_y$ are continuous everywhere: $f$ is $C^1$. On the axes, $f_x(0,y)=-y$ and $f_y(x,0)=x$ (including at $0$). Therefore

$$
f_{xy}(0,0)=\lim_{k\to0}\frac{f_x(0,k)-f_x(0,0)}{k}=-1,\qquad f_{yx}(0,0)=\lim_{h\to0}\frac{f_y(h,0)-f_y(0,0)}{h}=+1 .
$$

The mixed partials exist everywhere but differ at the origin. Away from the origin both equal

$$
f_{xy}=f_{yx}=\frac{(x^2-y^2)(x^4+10x^2y^2+y^4)}{(x^2+y^2)^3}=\cos2\phi\,\big(1+2\sin^22\phi\big),
$$

which depends only on the direction $\phi$: it is $+1$ along the $x$-axis and $-1$ along the $y$-axis, so it has no limit at the origin. The continuity hypothesis of Theorem A.11 fails exactly there, and the conclusion fails with it.

**Remarks.**

- *Higher order.* If $f\in C^k$ on an open set in $\mathbb{R}^n$, any two $k$-th order partial derivatives that differentiate the same number of times in each variable are equal. Proof: any ordering can be turned into any other by swapping adjacent differentiations, and each swap is Theorem A.11 applied to a $(k-2)$-th order partial (which is $C^2$) as a function of the two variables being swapped, the others held fixed.
- *Why engineers may swap orders freely — and when not.* Potentials and fields from smooth physical models are $C^2$ away from sources. At a point charge or on a conductor edge they are not, and statements such as $\nabla\times\nabla\varphi=\mathbf{0}$ must be checked separately (or interpreted in the sense of distributions).
- *Where it is used.* Symmetry of the Hessian ([Ch. 6 §6.2, §6.10](../part-2-calculus/06-multivariable-calculus.md); [Ch. 23 §23.4](../part-7-numerical-computational/23-optimization.md)); the exactness test $M_y=N_x$ for $M\,dx+N\,dy$ ([Ch. 6 §6.3](../part-2-calculus/06-multivariable-calculus.md), [Ch. 11 §11.1](../part-4-differential-equations/11-ordinary-differential-equations.md)); $\nabla\times\nabla\varphi=\mathbf{0}$ and $\nabla\cdot(\nabla\times\mathbf{F})=0$ ([Ch. 7](../part-2-calculus/07-vector-calculus.md)); the Cauchy–Riemann equations and harmonicity of $\operatorname{Re}\{f\}$ ([Ch. 18](../part-5-transforms/18-complex-analysis.md)).

## A.10 The Lagrange multiplier condition

The method of Lagrange multipliers ([Ch. 6 §6.11](../part-2-calculus/06-multivariable-calculus.md)) rests on one geometric fact: at a constrained extremum, $\nabla f$ is perpendicular to the constraint surface. To make "surface" precise we need the implicit function theorem, which says that near a regular point the constraint set is the graph of a smooth function. Write points of $\mathbb{R}^n$ as $\mathbf{x}=(\mathbf{u},\mathbf{v})$ with $\mathbf{u}\in\mathbb{R}^{n-m}$, $\mathbf{v}\in\mathbb{R}^m$, and Jacobian blocks $D_{\mathbf{u}}$, $D_{\mathbf{v}}$; gradients are column vectors and $Df=\nabla f^{\mathsf T}$ is a row.

> [!theorem] Implicit function theorem (used without proof)
> Let $\Omega\subset\mathbb{R}^{n-m}\times\mathbb{R}^m$ be open, $\mathbf{G}:\Omega\to\mathbb{R}^m$ of class $C^1$, and $(\mathbf{u}^*,\mathbf{v}^*)\in\Omega$ with $\mathbf{G}(\mathbf{u}^*,\mathbf{v}^*)=\mathbf{0}$ and the $m\times m$ matrix $D_{\mathbf{v}}\mathbf{G}(\mathbf{u}^*,\mathbf{v}^*)$ invertible. Then there are open sets $\mathcal{U}\ni\mathbf{u}^*$ and $\mathcal{V}\ni\mathbf{v}^*$ with $\mathcal{U}\times\mathcal{V}\subset\Omega$, and a $C^1$ function $\boldsymbol{\varphi}:\mathcal{U}\to\mathcal{V}$ with $\boldsymbol{\varphi}(\mathbf{u}^*)=\mathbf{v}^*$, such that
> $$
> \{(\mathbf{u},\mathbf{v})\in\mathcal{U}\times\mathcal{V}:\ \mathbf{G}(\mathbf{u},\mathbf{v})=\mathbf{0}\}=\{(\mathbf{u},\boldsymbol{\varphi}(\mathbf{u})):\ \mathbf{u}\in\mathcal{U}\}.
> $$
> (Proof: Rudin, *Principles*, Thm 9.28, via the contraction mapping principle. The derivative formula $D\boldsymbol\varphi=-(D_{\mathbf v}\mathbf G)^{-1}D_{\mathbf u}\mathbf G$ quoted in [Ch. 6 §6.6](../part-2-calculus/06-multivariable-calculus.md) is re-derived in Step 3 below.)

> [!theorem] Theorem A.12 (Lagrange multiplier condition)
> Let $U\subset\mathbb{R}^n$ be open, $1\le m<n$, and $f,g_1,\dots,g_m:U\to\mathbb{R}$ of class $C^1$. Let $S=\{\mathbf{x}\in U:\ g_i(\mathbf{x})=c_i,\ i=1,\dots,m\}$ and suppose $\mathbf{x}^*\in S$ is a **local extremum of $f$ on $S$**: there is $\rho>0$ with $f(\mathbf{x})\le f(\mathbf{x}^*)$ for all $\mathbf{x}\in S$ with $\lVert\mathbf{x}-\mathbf{x}^*\rVert<\rho$ (or $\ge$ for all such $\mathbf{x}$). Suppose the **regularity condition** holds: $\nabla g_1(\mathbf{x}^*),\dots,\nabla g_m(\mathbf{x}^*)$ are linearly independent. Then there are unique scalars $\lambda_1,\dots,\lambda_m$ with
> $$
> \nabla f(\mathbf{x}^*)=\lambda_1\nabla g_1(\mathbf{x}^*)+\dots+\lambda_m\nabla g_m(\mathbf{x}^*).
> $$
> For a single constraint ($m=1$) the condition reads $\nabla f(\mathbf{x}^*)=\lambda\nabla g(\mathbf{x}^*)$ with $\nabla g(\mathbf{x}^*)\ne\mathbf{0}$.

**Proof.** Treat the maximum case (for a minimum replace $f$ by $-f$, which replaces $\boldsymbol\lambda$ by $-\boldsymbol\lambda$). Let $\mathbf{g}=(g_1,\dots,g_m)^{\mathsf T}$ and $\mathbf{c}=(c_1,\dots,c_m)^{\mathsf T}$.

*Step 0: choose coordinates.* The $m\times n$ Jacobian $D\mathbf{g}(\mathbf{x}^*)$ has linearly independent rows, so its rank is $m$ and it has $m$ linearly independent columns. Renumber the coordinates (which permutes the components of every gradient in the same way and does not affect the conclusion) so that these are the last $m$ columns. Writing $\mathbf{x}=(\mathbf{u},\mathbf{v})$ and $\mathbf{x}^*=(\mathbf{u}^*,\mathbf{v}^*)$,

$$
D\mathbf{g}(\mathbf{x}^*)=\big[\,\mathbf{A}\;\;\mathbf{B}\,\big],\qquad \mathbf{A}=D_{\mathbf{u}}\mathbf{g}(\mathbf{x}^*)\ (m\times(n-m)),\quad \mathbf{B}=D_{\mathbf{v}}\mathbf{g}(\mathbf{x}^*)\ (m\times m,\ \text{invertible}).
$$

*Step 1: the constraint set is a graph.* Apply the implicit function theorem to $\mathbf{G}=\mathbf{g}-\mathbf{c}$ on $\Omega=U$. We get open $\mathcal U\ni\mathbf u^*$, $\mathcal V\ni\mathbf v^*$ and a $C^1$ map $\boldsymbol\varphi:\mathcal U\to\mathcal V$ with $\boldsymbol\varphi(\mathbf u^*)=\mathbf v^*$ such that the points of $S$ in $\mathcal U\times\mathcal V$ are exactly the points $(\mathbf u,\boldsymbol\varphi(\mathbf u))$.

*Step 2: an unconstrained problem.* Define $F(\mathbf u)=f(\mathbf u,\boldsymbol\varphi(\mathbf u))$ on $\mathcal U$; it is $C^1$ by the chain rule. Because $\boldsymbol\varphi$ is continuous, there is a neighbourhood $\mathcal U'\subset\mathcal U$ of $\mathbf u^*$ on which $\lVert(\mathbf u,\boldsymbol\varphi(\mathbf u))-\mathbf x^*\rVert<\rho$. For $\mathbf u\in\mathcal U'$ the point $(\mathbf u,\boldsymbol\varphi(\mathbf u))$ lies in $S$ and within $\rho$ of $\mathbf x^*$, so $F(\mathbf u)\le f(\mathbf x^*)=F(\mathbf u^*)$: $F$ has an unconstrained local maximum at $\mathbf u^*$. For each $i=1,\dots,n-m$, the one-variable function $t\mapsto F(\mathbf u^*+t\mathbf e_i)$ is defined on an open interval around $0$, has a maximum at $t=0$ and is differentiable there, so by Fermat's lemma (Lemma A.1) $\partial F/\partial u_i(\mathbf u^*)=0$. By the chain rule,

$$
\mathbf 0^{\mathsf T}=DF(\mathbf u^*)=D_{\mathbf u}f(\mathbf x^*)+D_{\mathbf v}f(\mathbf x^*)\,D\boldsymbol\varphi(\mathbf u^*).
$$

*Step 3: the slope of the constraint set.* Differentiating the identity $\mathbf g(\mathbf u,\boldsymbol\varphi(\mathbf u))=\mathbf c$ on $\mathcal U$ at $\mathbf u^*$ (chain rule),

$$
\mathbf A+\mathbf B\,D\boldsymbol\varphi(\mathbf u^*)=\mathbf 0\qquad\Longrightarrow\qquad D\boldsymbol\varphi(\mathbf u^*)=-\mathbf B^{-1}\mathbf A .
$$

*Step 4: define the multipliers.* Let $\boldsymbol\lambda^{\mathsf T}=D_{\mathbf v}f(\mathbf x^*)\,\mathbf B^{-1}$, a $1\times m$ row vector. Then $\boldsymbol\lambda^{\mathsf T}\mathbf B=D_{\mathbf v}f(\mathbf x^*)$ by construction, and by Steps 2 and 3,

$$
\boldsymbol\lambda^{\mathsf T}\mathbf A=D_{\mathbf v}f\,\mathbf B^{-1}\mathbf A=-D_{\mathbf v}f\,D\boldsymbol\varphi(\mathbf u^*)=D_{\mathbf u}f(\mathbf x^*).
$$

Putting the two blocks side by side, $\boldsymbol\lambda^{\mathsf T}\big[\mathbf A\;\;\mathbf B\big]=\big[D_{\mathbf u}f\;\;D_{\mathbf v}f\big]$, i.e. $Df(\mathbf x^*)=\boldsymbol\lambda^{\mathsf T}D\mathbf g(\mathbf x^*)$. Transposing, $\nabla f(\mathbf x^*)=D\mathbf g(\mathbf x^*)^{\mathsf T}\boldsymbol\lambda=\sum_i\lambda_i\nabla g_i(\mathbf x^*)$.

*Uniqueness.* If also $\nabla f=\sum_i\mu_i\nabla g_i$, then $\sum_i(\lambda_i-\mu_i)\nabla g_i(\mathbf x^*)=\mathbf 0$, and linear independence forces $\lambda_i=\mu_i$. ∎

For one constraint ($m=1$), after renumbering so that $g_{x_n}(\mathbf x^*)\ne0$, the proof gives the explicit multiplier $\lambda=f_{x_n}(\mathbf x^*)/g_{x_n}(\mathbf x^*)$, and Step 3 is the familiar implicit-differentiation formula $\partial x_n/\partial x_i=-g_{x_i}/g_{x_n}$ of [Ch. 4 §4.4](../part-2-calculus/04-differential-calculus.md).

**Geometric reading.** The columns of $\begin{bmatrix}\mathbf I\\ D\boldsymbol\varphi(\mathbf u^*)\end{bmatrix}$ span the tangent space $T=\ker D\mathbf g(\mathbf x^*)$ of $S$ at $\mathbf x^*$ (an $(n-m)$-dimensional subspace, by Step 3 and the rank–nullity theorem of [Ch. 9](../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md)). Step 2 says $Df(\mathbf x^*)$ annihilates $T$: $\nabla f\perp T$. The orthogonal complement of $\ker D\mathbf g$ is the row space of $D\mathbf g$, spanned by the $\nabla g_i$ — which is Step 4.

```desmos-graph
left=-4; right=4; bottom=-3.2; top=3.2
---
x^{2}+y^{2}=5|blue
x+2y=5|red
x+2y=-5|red
x+2y=2|red|dashed
x+2y=-2|red|dashed
y=2x|green|1<x<1.6
(1,2)|label:max
(-1,-2)|label:min
```

*Extremising $f=x+2y$ on the circle $x^2+y^2=5$ ([Ch. 6](../part-2-calculus/06-multivariable-calculus.md), Worked Example 1). Level lines of $f$ (red) cross the circle where $f$ can still increase along it; at the extremes $(\pm1,\pm2)$ they are tangent, so $\nabla f=(1,2)$ (green, at the maximum) is parallel to $\nabla g=(2x,2y)$, with $\lambda=\tfrac12$ at $(1,2)$.*

**Remarks.**

- *Regularity cannot be dropped.* Minimise $f(x,y)=x$ on $S:\ y^2-x^3=0$. On $S$, $x^3=y^2\ge0$, so $x\ge0$ and the minimum is $0$, at the origin. But $\nabla g=(-3x^2,\,2y)$ vanishes there, while $\nabla f=(1,0)$; no $\lambda$ satisfies $(1,0)=\lambda(0,0)$. The cusp of $S$ at the origin is not the graph of a $C^1$ function in either variable, so Step 1 fails. (The weaker **Fritz John** condition $\lambda_0\nabla f=\lambda\nabla g$ with $(\lambda_0,\lambda)\ne(0,0)$ survives, with $\lambda_0=0$; it carries no information about $f$.)
- *The condition is necessary, not sufficient.* $f=x^3+y$ on $g=y=0$: at the origin $\nabla f=(0,1)=1\cdot\nabla g$, yet on the constraint $f=x^3$ has no extremum at $0$. Candidates must be compared, or tested with the bordered Hessian ([Ch. 6 §6.11](../part-2-calculus/06-multivariable-calculus.md)).
- *The meaning of $\lambda$.* Suppose the constrained optimum $\mathbf x^*(c)$ and its multiplier $\lambda(c)$ depend differentiably on the level $c$ of a single constraint, and let $f^*(c)=f(\mathbf x^*(c))$. Then, by the chain rule and the theorem,
$$
\frac{df^*}{dc}=\nabla f\cdot\frac{d\mathbf x^*}{dc}=\lambda\,\nabla g\cdot\frac{d\mathbf x^*}{dc}=\lambda\,\frac{d}{dc}g\big(\mathbf x^*(c)\big)=\lambda\cdot\frac{d}{dc}c=\lambda .
$$
The multiplier is the sensitivity of the optimal value to the constraint level — a "shadow price". Check: on $x^2+y^2=c$, $\max(x+2y)=\sqrt{5c}$, whose derivative $\sqrt5/(2\sqrt c)$ equals $\tfrac12=\lambda$ at $c=5$.
- *Where it is used.* [Ch. 6 §6.11](../part-2-calculus/06-multivariable-calculus.md) (the method), [Ch. 23 §23.4](../part-7-numerical-computational/23-optimization.md) (Lagrange and KKT conditions — KKT is this theorem plus sign conditions for active inequality constraints, under an analogous regularity hypothesis), equality-constrained least squares and economic dispatch (equal incremental cost $=\lambda$) in [Ch. 23 §23.8, §23.10](../part-7-numerical-computational/23-optimization.md), and the variational characterisation of eigenvalues: maximising $\mathbf x^{\mathsf T}\mathbf A\mathbf x$ for symmetric $\mathbf A$ subject to $\mathbf x^{\mathsf T}\mathbf x=1$ gives $2\mathbf A\mathbf x=\lambda\,2\mathbf x$, the eigenvalue equation, with the multiplier equal to the maximum value ([Ch. 10](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md)).

## Further reading

- T. M. Apostol, *Calculus*, Vols. I–II, 2nd ed., Wiley, 1967–1969 — rigorous proofs of the MVT, FTC, Taylor's theorem and the multivariable results at the level of this appendix.
- M. Spivak, *Calculus*, 4th ed., Publish or Perish, 2008 — the cleanest treatment of the single-variable results, including the counterexamples used here.
- W. Rudin, *Principles of Mathematical Analysis*, 3rd ed., McGraw-Hill, 1976 — L'Hôpital's rule in full generality (Thm 5.13), the inverse and implicit function theorems (Ch. 9), and the exponential function defined by its series (Ch. 8).
- T. Needham, *Visual Complex Analysis*, Oxford University Press, 1997 — the geometric meaning of $e^{j\theta}$, De Moivre's theorem and roots of unity.
