# Appendix B — Vector Calculus and Field Theorems

[← Table of Contents](../00-front/table-of-contents.md)

> **What is proved here.** The results that make vector calculus usable in field theory: the change-of-variables formula and the Jacobian as an area scale factor; the gradient as the direction of steepest ascent; the theorems of Green, Gauss (divergence) and Stokes; the equivalence "curl-free ⇔ conservative" on simply connected domains, with its failure on the punctured plane; the standard vector identities by index notation; grad, div, curl and Laplacian in cylindrical and spherical coordinates via scale factors; and the passage from Maxwell's integral laws to their differential form, the continuity equation and the vacuum wave equation.
>
> **Prerequisites.** Partial derivatives, the multivariable chain rule, Clairaut's theorem on mixed partials and Fubini's theorem ([Ch. 6](../part-2-calculus/06-multivariable-calculus.md)); determinants ([Ch. 8](../part-3-linear-algebra/08-matrices-and-linear-systems.md)); line, surface and volume integrals as defined in [Ch. 7](../part-2-calculus/07-vector-calculus.md).
>
> **Level of rigour.** Every theorem is stated with explicit hypotheses. The proofs are complete for the "elementary" regions on which the integral theorems are really proved (regions bounded by graphs). The step from elementary regions to general regions — cutting a region into finitely many elementary pieces — is carried out as an argument about cancelling boundary contributions, and the topological facts it uses (that such a decomposition exists for the piecewise-smooth regions met in engineering) are stated, not proved. In the change-of-variables theorem the linear-algebra step is proved in full and the passage to nonlinear maps is a careful sketch. Throughout, "smooth" means $C^1$ (continuous first partial derivatives) unless $C^2$ is stated; regions are bounded and their boundaries piecewise smooth.

Notation follows the book ([Notation and Conventions](../00-front/notation-and-conventions.md)): cylindrical coordinates are $(r,\phi,z)$, spherical coordinates are $(R,\theta,\phi)$ with $\theta$ the polar angle from $+z$; $\partial_i$ abbreviates $\partial/\partial x_i$ with $(x_1,x_2,x_3)=(x,y,z)$.

## B.1 Change of variables and the Jacobian as area scale factor

The change-of-variables formula is used every time the book writes $dA=r\,dr\,d\phi$ or $dV=R^2\sin\theta\,dR\,d\theta\,d\phi$ ([§6.8](../part-2-calculus/06-multivariable-calculus.md), [§6.9](../part-2-calculus/06-multivariable-calculus.md)), and every time a probability density is transformed ([Ch. 19](../part-6-probability-statistics/19-probability.md)). Its heart is a statement about linear maps, which we prove exactly; the nonlinear statement then follows because a differentiable map is, on a small enough square, almost linear.

> [!definition] Jacobian matrix and determinant
> For a map $\mathbf{T}(u,v)=\big(x(u,v),\,y(u,v)\big)$ with continuous first partials, the **Jacobian matrix** and **Jacobian determinant** at $(u,v)$ are
> $$
> D\mathbf{T}=\begin{bmatrix}x_u & x_v\\ y_u & y_v\end{bmatrix},\qquad
> \frac{\partial(x,y)}{\partial(u,v)}=\det D\mathbf{T}=x_uy_v-x_vy_u .
> $$

> [!theorem] Theorem B.1a — Linear maps scale area by $\lvert\det\rvert$
> Let $\mathbf{A}$ be a real $2\times2$ matrix and $S\subset\mathbb{R}^2$ a bounded region with an area (for instance a finite union of rectangles, or any region bounded by finitely many smooth curves). Then $\mathbf{A}(S)=\{\mathbf{A}\mathbf{s}:\mathbf{s}\in S\}$ has an area and
> $$
> \operatorname{area}\big(\mathbf{A}(S)\big)=\lvert\det\mathbf{A}\rvert\;\operatorname{area}(S).
> $$

> [!theorem] Theorem B.1b — Change of variables in a double integral
> Let $U\subset\mathbb{R}^2$ be open, $\mathbf{T}:U\to\mathbb{R}^2$ one-to-one with continuous first partial derivatives, and $\det D\mathbf{T}\ne0$ on $U$. Let $S$ be a closed bounded region in $U$ with piecewise-smooth boundary, and $R=\mathbf{T}(S)$. If $f$ is continuous on $R$, then
> $$
> \iint_R f(x,y)\,dx\,dy=\iint_S f\big(\mathbf{T}(u,v)\big)\left\lvert\frac{\partial(x,y)}{\partial(u,v)}\right\rvert du\,dv .
> $$
> The same holds in $\mathbb{R}^3$ with $\lvert\partial(x,y,z)/\partial(u,v,w)\rvert$, and the hypothesis $\det D\mathbf{T}\ne0$ may fail on a set of zero area (such as the origin in polar coordinates) without changing the conclusion.

**Proof of B.1a.**

*Step 1 — reduction to rectangles.* If the formula holds for every axis-parallel rectangle, it holds for finite unions of non-overlapping rectangles (areas add, and the images of non-overlapping rectangles under an invertible linear map are non-overlapping). A region with an area is, by definition of area (the Riemann/Jordan construction behind the double integral in [§6.7](../part-2-calculus/06-multivariable-calculus.md)), squeezed between an inner union of grid rectangles $S_-$ and an outer union $S_+$ with $\operatorname{area}(S_+)-\operatorname{area}(S_-)<\varepsilon$. Then $\mathbf{A}(S_-)\subset\mathbf{A}(S)\subset\mathbf{A}(S_+)$ and the two outer sets have areas $\lvert\det\mathbf{A}\rvert\operatorname{area}(S_\pm)$ differing by at most $\lvert\det\mathbf{A}\rvert\varepsilon$. Letting $\varepsilon\to0$ gives the formula for $S$. If $\det\mathbf{A}=0$, $\mathbf{A}$ maps the plane into a line (or a point), and any bounded subset of a line has zero area, which is what the formula says. So assume $\det\mathbf{A}\ne0$ and $S$ a rectangle $[a,b]\times[c,d]$.

*Step 2 — factor $\mathbf{A}$ into elementary matrices.* Gaussian elimination ([§8.7](../part-3-linear-algebra/08-matrices-and-linear-systems.md)) reduces any invertible matrix to the identity by row operations, each of which is left-multiplication by an elementary matrix. Hence $\mathbf{A}=\mathbf{E}_1\mathbf{E}_2\cdots\mathbf{E}_k$ with each $\mathbf{E}_i$ of one of three types:

$$
\text{scaling } \begin{bmatrix}\lambda&0\\0&1\end{bmatrix}\text{ or }\begin{bmatrix}1&0\\0&\lambda\end{bmatrix}\ (\lambda\ne0),\qquad
\text{swap } \begin{bmatrix}0&1\\1&0\end{bmatrix},\qquad
\text{shear } \begin{bmatrix}1&c\\0&1\end{bmatrix}\text{ or }\begin{bmatrix}1&0\\c&1\end{bmatrix}.
$$

Since $\det(\mathbf{A})=\prod\det\mathbf{E}_i$ and $\mathbf{A}(S)=\mathbf{E}_1(\mathbf{E}_2(\cdots\mathbf{E}_k(S)))$, it suffices to prove $\operatorname{area}(\mathbf{E}(S))=\lvert\det\mathbf{E}\rvert\operatorname{area}(S)$ for each elementary $\mathbf{E}$ and every region $S$ with an area — not just rectangles, because the intermediate images $\mathbf{E}_k(S)$, $\mathbf{E}_{k-1}\mathbf{E}_k(S),\dots$ are parallelograms. By Step 1 (applied to the single matrix $\mathbf{E}$) the general case follows from the case of rectangles, so it is enough to check rectangles for each type.

*Step 3 — the three types.*

- *Scaling* $(u,v)\mapsto(\lambda u,v)$ sends $[a,b]\times[c,d]$ to a rectangle with side lengths $\lvert\lambda\rvert(b-a)$ and $d-c$: area multiplied by $\lvert\lambda\rvert=\lvert\det\mathbf{E}\rvert$.
- *Swap* $(u,v)\mapsto(v,u)$ sends a rectangle to a congruent rectangle: area unchanged, and $\lvert\det\mathbf{E}\rvert=\lvert-1\rvert=1$.
- *Shear* $(u,v)\mapsto(u+cv,\,v)$ sends $[a,b]\times[c_1,d_1]$ to the parallelogram $P=\{(x,y):c_1\le y\le d_1,\ a+cy\le x\le b+cy\}$. By Fubini, slicing horizontally,
$$
\operatorname{area}(P)=\int_{c_1}^{d_1}\big[(b+cy)-(a+cy)\big]\,dy=(b-a)(d_1-c_1),
$$
so the area is unchanged, and $\det\mathbf{E}=1$. (This is Cavalieri's principle: every horizontal slice keeps its length.) The lower shear is the same with the roles of $x$ and $y$ exchanged.

Multiplying the factors gives $\operatorname{area}(\mathbf{A}(S))=\prod_i\lvert\det\mathbf{E}_i\rvert\operatorname{area}(S)=\lvert\det\mathbf{A}\rvert\operatorname{area}(S)$. ∎

```desmos-graph
left=-0.5; right=3.5; bottom=-0.5; top=3.5
---
y=0|blue|0<x<1
x=1|blue|0<y<1
y=1|blue|0<x<1
x=0|blue|0<y<1
y=\frac{x}{2}|red|0<x<2
y=2x|red|0<x<1
y=2x-3|red|2<x<3
y=\frac{x}{2}+1.5|red|1<x<3
```
*The unit square (blue) and its image under $\mathbf{A}=\begin{bmatrix}2&1\\1&2\end{bmatrix}$ (red): a parallelogram spanned by $(2,1)$ and $(1,2)$, of area $\lvert\det\mathbf{A}\rvert=3$.*

**Proof sketch of B.1b.**

*Step 1 — local linearity, uniformly.* Fix a compact set $K\subset U$ containing $S$ in its interior. Because the partial derivatives of $\mathbf{T}$ are continuous on $K$ they are uniformly continuous there, so for every $\varepsilon>0$ there is $h_0>0$ such that for every square $Q\subset K$ of side $h<h_0$ with centre $\mathbf{q}$ and every $\mathbf{u}\in Q$,

$$
\big\lVert\mathbf{T}(\mathbf{u})-\mathbf{T}(\mathbf{q})-D\mathbf{T}(\mathbf{q})(\mathbf{u}-\mathbf{q})\big\rVert\le\varepsilon h .
$$

(Apply the mean value theorem to each component along the segment from $\mathbf{q}$ to $\mathbf{u}$.) So $\mathbf{T}(Q)$ lies within distance $\varepsilon h$ of the parallelogram $\mathbf{T}(\mathbf{q})+D\mathbf{T}(\mathbf{q})(Q-\mathbf{q})$, whose area by Theorem B.1a is $\lvert\det D\mathbf{T}(\mathbf{q})\rvert h^2$.

*Step 2 — area of the image of a small square.* The $\varepsilon h$-neighbourhood of that parallelogram's boundary has area at most (perimeter)$\times 2\varepsilon h+\pi(\varepsilon h)^2\le C\varepsilon h^2$, where $C$ depends only on a bound for $\lVert D\mathbf{T}\rVert$ on $K$. Since $\mathbf{T}$ is one-to-one with invertible derivative, it maps the boundary of $Q$ onto the boundary of $\mathbf{T}(Q)$ (invariance of domain, here a consequence of the inverse function theorem, [§6.6](../part-2-calculus/06-multivariable-calculus.md)), and therefore $\mathbf{T}(Q)$ differs from the parallelogram only inside that thin neighbourhood. Hence

$$
\Big\lvert\operatorname{area}\big(\mathbf{T}(Q)\big)-\lvert\det D\mathbf{T}(\mathbf{q})\rvert\,h^2\Big\rvert\le C\varepsilon h^2 .
$$

*Step 3 — Riemann sums.* Cover $S$ by a grid of squares $Q_k$ of side $h$ with centres $\mathbf{q}_k$. The squares meeting $\partial S$ have total area $\to0$ as $h\to0$ (a piecewise-smooth curve has zero area), and so do their images (Step 2 bounds image areas by a constant times $h^2$). For the interior squares, the images $\mathbf{T}(Q_k)$ tile $R$ up to a set of small area, and

$$
\sum_k f\big(\mathbf{T}(\mathbf{q}_k)\big)\operatorname{area}\big(\mathbf{T}(Q_k)\big)
=\sum_k f\big(\mathbf{T}(\mathbf{q}_k)\big)\lvert\det D\mathbf{T}(\mathbf{q}_k)\rvert h^2+O\!\left(\varepsilon\max\lvert f\rvert\operatorname{area}(S)\right).
$$

The left side is a Riemann sum for $\iint_Rf\,dA$ (over a partition of $R$ into the curvilinear pieces $\mathbf{T}(Q_k)$, whose diameters $\to0$ uniformly), and the first sum on the right is a Riemann sum for $\iint_Sf(\mathbf{T})\lvert\det D\mathbf{T}\rvert\,du\,dv$. Let $h\to0$ and then $\varepsilon\to0$. If $\det D\mathbf{T}=0$ or injectivity fails only on a set $Z$ of zero area, apply the result on $S$ minus a small neighbourhood of $Z$ and let the neighbourhood shrink; both integrals change by arbitrarily little. In three dimensions the same argument runs with cubes and parallelepipeds, using the $3\times3$ elementary-matrix factorisation. ∎

**Remarks.**

- *Why the absolute value.* An integral of a positive function is positive, but $\det D\mathbf{T}$ is negative for orientation-reversing maps (e.g. the swap). The absolute value is the price of measuring unsigned area. In the integral theorems (B.3–B.5) orientation matters and the sign of the Jacobian reappears as the rule for orienting boundaries.
- *Why injectivity.* If $\mathbf{T}$ covers part of $R$ twice the right side counts that part twice. Example: $\mathbf{T}(r,\phi)=(r\cos\phi,r\sin\phi)$ on $0\le r\le1$, $0\le\phi\le4\pi$ gives $\iint r\,dr\,d\phi=2\pi$, twice the area of the unit disc.
- *Scale factors.* For an orthogonal coordinate system the columns of $D\mathbf{T}$ are perpendicular with lengths $h_1,h_2,h_3$, so $\lvert\det D\mathbf{T}\rvert=h_1h_2h_3$: this gives $r$ (cylindrical) and $R^2\sin\theta$ (spherical) without expanding a determinant, as used in [§6.9](../part-2-calculus/06-multivariable-calculus.md) and again in B.8 below.
- *Used in:* polar, cylindrical and spherical integrals ([§6.8–6.9](../part-2-calculus/06-multivariable-calculus.md)); the volume elements of [§7.8](../part-2-calculus/07-vector-calculus.md); the transformation of densities $f_Y(\mathbf{y})=f_X(\mathbf{x})/\lvert\det D\mathbf{T}\rvert$ in [Ch. 19](../part-6-probability-statistics/19-probability.md).

## B.2 The gradient as the direction of steepest ascent

[§6.4–6.5](../part-2-calculus/06-multivariable-calculus.md) define the gradient by its components and state that it points uphill. The proof rests on one fact — that for a differentiable function every directional derivative is a dot product with $\nabla f$ — and on the Cauchy–Schwarz inequality. The point of stating the hypotheses carefully is that the existence of partial derivatives alone is **not** enough.

> [!definition] Differentiability and directional derivative
> $f:U\to\mathbb{R}$ ($U\subset\mathbb{R}^n$ open) is **differentiable** at $\mathbf{a}$ if there is a vector $\mathbf{g}$ with
> $$
> f(\mathbf{a}+\mathbf{h})=f(\mathbf{a})+\mathbf{g}\cdot\mathbf{h}+o(\lVert\mathbf{h}\rVert)\quad(\mathbf{h}\to\mathbf{0}).
> $$
> Then necessarily $\mathbf{g}=\nabla f(\mathbf{a})=(\partial_1f,\dots,\partial_nf)(\mathbf{a})$. For a unit vector $\hat{\mathbf{u}}$ the **directional derivative** is $D_{\hat{\mathbf{u}}}f(\mathbf{a})=\lim_{s\to0^+}\big[f(\mathbf{a}+s\hat{\mathbf{u}})-f(\mathbf{a})\big]/s$.

> [!theorem] Theorem B.2 — Steepest ascent and normality to level sets
> Let $f$ be differentiable at $\mathbf{a}$.
> 1. For every unit vector $\hat{\mathbf{u}}$, $D_{\hat{\mathbf{u}}}f(\mathbf{a})=\nabla f(\mathbf{a})\cdot\hat{\mathbf{u}}$.
> 2. If $\nabla f(\mathbf{a})\ne\mathbf{0}$, then $\max_{\lVert\hat{\mathbf{u}}\rVert=1}D_{\hat{\mathbf{u}}}f(\mathbf{a})=\lVert\nabla f(\mathbf{a})\rVert$, attained **only** at $\hat{\mathbf{u}}=\nabla f(\mathbf{a})/\lVert\nabla f(\mathbf{a})\rVert$; the minimum $-\lVert\nabla f(\mathbf{a})\rVert$ is attained only in the opposite direction; and $D_{\hat{\mathbf{u}}}f(\mathbf{a})=0$ exactly for $\hat{\mathbf{u}}\perp\nabla f(\mathbf{a})$.
> 3. If $\mathbf{r}(t)$ is a curve, differentiable at $t=0$, lying in the level set $\{f=f(\mathbf{a})\}$ with $\mathbf{r}(0)=\mathbf{a}$, then $\nabla f(\mathbf{a})\cdot\mathbf{r}'(0)=0$.

**Proof.**

*Part 1.* Put $\mathbf{h}=s\hat{\mathbf{u}}$ in the definition: $f(\mathbf{a}+s\hat{\mathbf{u}})-f(\mathbf{a})=s\,\nabla f(\mathbf{a})\cdot\hat{\mathbf{u}}+o(s)$. Divide by $s>0$ and let $s\to0^+$; the remainder tends to $0$ because $o(\lVert\mathbf{h}\rVert)/\lVert\mathbf{h}\rVert\to0$ and $\lVert\mathbf{h}\rVert=s$.

*Part 2.* By Cauchy–Schwarz, $\lvert\nabla f\cdot\hat{\mathbf{u}}\rvert\le\lVert\nabla f\rVert\,\lVert\hat{\mathbf{u}}\rVert=\lVert\nabla f\rVert$, so $-\lVert\nabla f\rVert\le D_{\hat{\mathbf{u}}}f\le\lVert\nabla f\rVert$. Equality in Cauchy–Schwarz holds if and only if $\hat{\mathbf{u}}$ is a scalar multiple of $\nabla f$; with $\lVert\hat{\mathbf{u}}\rVert=1$ the only multiples are $\pm\nabla f/\lVert\nabla f\rVert$, which give $\pm\lVert\nabla f\rVert$ respectively. Writing $\nabla f\cdot\hat{\mathbf{u}}=\lVert\nabla f\rVert\cos\alpha$ shows the whole picture: the rate of increase is the projection of the gradient onto the direction of travel, and it is zero exactly when $\alpha=90^\circ$.

*Part 3.* The function $g(t)=f(\mathbf{r}(t))$ is constant, so $g'(0)=0$. By the chain rule (valid because $f$ is differentiable at $\mathbf{a}$ and $\mathbf{r}$ at $0$), $g'(0)=\nabla f(\mathbf{a})\cdot\mathbf{r}'(0)$. ∎

```desmos-graph
left=-2.5; right=2.5; bottom=-1.8; top=1.8
---
x^2+2y^2=3|blue
x^2+2y^2=1.5|blue|dashed
x^2+2y^2=5|blue|dashed
y-1=2(x-1)|red|0.6<x<1.6
x+2y=3|grey|0<x<2.2
(1,1)|label:a
```
*Level curves of $f=x^2+2y^2$. At $\mathbf{a}=(1,1)$ the gradient direction $(2,4)$ (red) is normal to the level curve $f=3$ and to its tangent line $x+2y=3$ (grey): it crosses the contours by the shortest route.*

**Remarks.**

- *Differentiability is essential.* Let $f(x,y)=x^2y/(x^2+y^2)$ for $(x,y)\ne\mathbf{0}$ and $f(\mathbf{0})=0$. Both partials at the origin are $0$ (the function vanishes on both axes), so the would-be gradient is $\mathbf{0}$; yet along $\hat{\mathbf{u}}=(1,1)/\sqrt2$, $f(s/\sqrt2,s/\sqrt2)/s=1/(2\sqrt2)\ne0$. The formula $D_{\hat{\mathbf{u}}}f=\nabla f\cdot\hat{\mathbf{u}}$ fails: this $f$ is continuous but not differentiable at $\mathbf{0}$. A convenient sufficient condition for differentiability is that the partials exist near $\mathbf{a}$ and are continuous at $\mathbf{a}$ ([§6.3](../part-2-calculus/06-multivariable-calculus.md)).
- *Units.* $\lVert\nabla f\rVert$ has units of $f$ per metre, which is why $\lVert\nabla V\rVert$ is a field strength in V/m. The relation $\mathbf{E}=-\nabla V$ says that $\mathbf{E}$ points in the direction of steepest **descent** of potential and is normal to the equipotentials (Part 3).
- *"Steepest" depends on the metric.* Steepest means per unit Euclidean distance. If distance is measured with a weighted norm $\lVert\mathbf{h}\rVert_{\mathbf{M}}=\sqrt{\mathbf{h}^\mathsf{T}\mathbf{M}\mathbf{h}}$ ($\mathbf{M}$ symmetric positive definite), the same Cauchy–Schwarz argument in the $\mathbf{M}$-inner product shows that the steepest direction is $\mathbf{M}^{-1}\nabla f$ — the idea behind preconditioned and Newton-type methods in [Ch. 23](../part-7-numerical-computational/23-optimization.md).
- *Used in:* gradient descent and the stationarity condition $\nabla f=\mathbf{0}$ ([§6.4](../part-2-calculus/06-multivariable-calculus.md), [Ch. 23](../part-7-numerical-computational/23-optimization.md)); Lagrange multipliers, where $\nabla f\parallel\nabla g$ follows from Part 3 applied to the constraint surface ([§6.11](../part-2-calculus/06-multivariable-calculus.md)); field lines orthogonal to equipotentials ([§7.5](../part-2-calculus/07-vector-calculus.md)).

## B.3 Green's theorem

Green's theorem is the two-dimensional prototype of all the integral theorems ([§7.9](../part-2-calculus/07-vector-calculus.md)). We prove it completely for regions that are simultaneously of "type I" and "type II", then extend it by cutting.

> [!definition] Elementary plane regions and positive orientation
> A closed bounded region $D$ is **type I** (vertically simple) if $D=\{(x,y):a\le x\le b,\ g_1(x)\le y\le g_2(x)\}$ with $g_1\le g_2$ continuous and piecewise $C^1$ on $[a,b]$; it is **type II** (horizontally simple) if $D=\{(x,y):c\le y\le d,\ h_1(y)\le x\le h_2(y)\}$ with the analogous conditions. The boundary $\partial D$ is **positively oriented** if $D$ lies to the left as one traverses it: counter-clockwise on an outer boundary, clockwise on the boundary of a hole.

> [!theorem] Theorem B.3 — Green's theorem
> Let $D$ be a closed bounded region whose boundary $\partial D$ consists of finitely many piecewise-smooth simple closed curves, positively oriented. Let $P,Q$ have continuous first partial derivatives on an open set containing $D$. Then
> $$
> \oint_{\partial D}\big(P\,dx+Q\,dy\big)=\iint_D\left(\frac{\partial Q}{\partial x}-\frac{\partial P}{\partial y}\right)dA .
> $$

**Proof.** Both sides are linear in $(P,Q)$, so it suffices to prove separately

$$
\oint_{\partial D}P\,dx=-\iint_D\frac{\partial P}{\partial y}\,dA\qquad\text{and}\qquad\oint_{\partial D}Q\,dy=\iint_D\frac{\partial Q}{\partial x}\,dA .
$$

*Step 1 — the $P$ identity on a type I region.* Since $\partial P/\partial y$ is continuous, Fubini's theorem and the fundamental theorem of calculus in $y$ give

$$
-\iint_D\frac{\partial P}{\partial y}\,dA=-\int_a^b\left[\int_{g_1(x)}^{g_2(x)}\frac{\partial P}{\partial y}(x,y)\,dy\right]dx=\int_a^b\Big[P\big(x,g_1(x)\big)-P\big(x,g_2(x)\big)\Big]dx .
$$

The positively oriented boundary of a type I region consists of (at most) four arcs:

- $C_1$, the bottom graph, $\mathbf{r}(t)=(t,g_1(t))$, $a\le t\le b$, traversed left to right. Here $dx=dt$, so $\int_{C_1}P\,dx=\int_a^bP(t,g_1(t))\,dt$.
- $C_2$, the right side $x=b$, $g_1(b)\le y\le g_2(b)$ (possibly a single point). Here $x$ is constant, so $dx=0$ and $\int_{C_2}P\,dx=0$.
- $C_3$, the top graph, traversed right to left. Parametrising by $(t,g_2(t))$ and reversing direction, $\int_{C_3}P\,dx=-\int_a^bP(t,g_2(t))\,dt$.
- $C_4$, the left side $x=a$; again $\int_{C_4}P\,dx=0$.

Adding, $\oint_{\partial D}P\,dx=\int_a^b\big[P(t,g_1(t))-P(t,g_2(t))\big]dt$, which is the expression above.

*Step 2 — the $Q$ identity on a type II region.* With $D=\{c\le y\le d,\ h_1(y)\le x\le h_2(y)\}$,

$$
\iint_D\frac{\partial Q}{\partial x}\,dA=\int_c^d\Big[Q\big(h_2(y),y\big)-Q\big(h_1(y),y\big)\Big]dy .
$$

On the positively oriented boundary the right graph $x=h_2(y)$ is traversed upward, contributing $+\int_c^dQ(h_2(y),y)\,dy$; the left graph $x=h_1(y)$ is traversed downward, contributing $-\int_c^dQ(h_1(y),y)\,dy$; the horizontal top and bottom segments contribute nothing because $dy=0$ on them. The asymmetry between Steps 1 and 2 is the origin of the minus sign in $Q_x-P_y$: going counter-clockwise, the *bottom* is traversed in the $+x$ direction but the *right* side in the $+y$ direction.

*Step 3 — regions of both types.* If $D$ is both type I and type II (a disc, a rectangle, a triangle, any convex region), Steps 1 and 2 apply to the same oriented boundary, and adding them proves the theorem for $D$.

*Step 4 — general regions by decomposition.* Suppose $D$ can be cut by finitely many piecewise-smooth arcs into regions $D_1,\dots,D_m$, each of both types, overlapping only along the cuts. (For regions bounded by finitely many piecewise-smooth curves such a decomposition exists; in practice one draws a few vertical and horizontal cuts.) Apply Step 3 to each $D_k$ with its own positive orientation and add:

$$
\sum_{k=1}^m\oint_{\partial D_k}(P\,dx+Q\,dy)=\sum_{k=1}^m\iint_{D_k}(Q_x-P_y)\,dA=\iint_D(Q_x-P_y)\,dA,
$$

the last equality because the cuts have zero area. On the left, each cut arc is a boundary arc of exactly two pieces, which lie on opposite sides of it. "Keep the region on your left" therefore traverses the arc in opposite directions for the two pieces, and since reversing a curve negates a line integral, the two contributions cancel. What remains are the arcs of $\partial D$, each traversed with $D$ on its left — that is, $\partial D$ with positive orientation. This proves the theorem for $D$.

*Step 5 — holes.* A region with holes, such as the annulus $1\le x^2+y^2\le4$, is handled by the same cutting: two radial cuts split the annulus into two pieces, each a finite union of elementary regions. After cancellation the outer circle is traversed counter-clockwise and the inner circle clockwise — exactly the positive orientation of the definition. ∎

```desmos-graph
left=-0.2; right=1.2; bottom=-0.2; top=1.2
---
y=x^2|blue|0<x<1
y=x|red|0<x<1
(0.6,0.36)|label:C1
(0.4,0.4)|label:C3
```
*The region $x^2\le y\le x$, $0\le x\le1$, which is both type I and type II. For $P=x-y^2$, $Q=xy$ both sides of Green's theorem equal $\iint_D3y\,dA=\tfrac15$ (checked symbolically): the bottom arc $C_1$ (blue, left to right) contributes $\tfrac{7}{10}$, the top arc $C_3$ (red, right to left) contributes $-\tfrac12$.*

**Remarks.**

- *The $C^1$ hypothesis cannot fail even at one point.* For $P=-y/(x^2+y^2)$, $Q=x/(x^2+y^2)$ on the closed unit disc, $Q_x-P_y=0$ wherever defined, yet $\oint_{\lvert\mathbf{r}\rvert=1}(P\,dx+Q\,dy)=2\pi$. Green's theorem is not violated: $P$ and $Q$ are not defined at the origin, which lies in $D$. Applied instead to the annulus $\varepsilon\le\lvert\mathbf{r}\rvert\le1$ it gives the correct statement — the counter-clockwise circulations around the two circles are equal. This example is the subject of B.6.
- *Area formula.* $P=-y/2$, $Q=x/2$ gives $\operatorname{area}(D)=\tfrac12\oint_{\partial D}(x\,dy-y\,dx)$, the principle of the planimeter and of the "shoelace" polygon-area formula ([§7.9](../part-2-calculus/07-vector-calculus.md)).
- *Flux form.* Replacing $(P,Q)$ by $(-F_y,F_x)$ gives $\oint_{\partial D}\mathbf{F}\cdot\hat{\mathbf{n}}\,ds=\iint_D\nabla\cdot\mathbf{F}\,dA$ with $\hat{\mathbf{n}}$ the outward normal: the 2-D divergence theorem, which B.4 lifts to three dimensions.
- *Used in:* the proof of Stokes' theorem (B.5); Cauchy's integral theorem in [Ch. 18](../part-5-transforms/18-complex-analysis.md), where Green applied to $f=u+jv$ together with the Cauchy–Riemann equations gives $\oint f\,dz=0$ (under the extra hypothesis that $f'$ is continuous, which Goursat's proof removes).

## B.4 The divergence theorem

The divergence theorem (Gauss's theorem) converts the flux through a closed surface into a volume integral of the divergence ([§7.11](../part-2-calculus/07-vector-calculus.md)). The proof has exactly the structure of B.3: one Cartesian component at a time, the fundamental theorem of calculus along one axis, and a cut-and-cancel argument for general regions.

> [!definition] Flux and outward normal
> For a piecewise-smooth surface $S$ with a chosen continuous unit normal $\hat{\mathbf{n}}$, the **flux** of $\mathbf{F}$ is $\iint_S\mathbf{F}\cdot d\mathbf{S}=\iint_S\mathbf{F}\cdot\hat{\mathbf{n}}\,dS$. If $S=\mathbf{r}(D)$ with $\mathbf{r}_u\times\mathbf{r}_v$ pointing along $\hat{\mathbf{n}}$, then $\iint_S\mathbf{F}\cdot d\mathbf{S}=\iint_D\mathbf{F}(\mathbf{r}(u,v))\cdot(\mathbf{r}_u\times\mathbf{r}_v)\,du\,dv$ ([§7.7](../part-2-calculus/07-vector-calculus.md)). For the boundary $\partial V$ of a solid $V$, $\hat{\mathbf{n}}$ is always the **outward** normal.

> [!definition] Elementary solids
> A closed bounded solid $V$ is **$z$-simple** if $V=\{(x,y,z):(x,y)\in D_{xy},\ g_1(x,y)\le z\le g_2(x,y)\}$ for a plane region $D_{xy}$ with piecewise-smooth boundary and $g_1\le g_2$ continuous on $D_{xy}$ and $C^1$ on its interior; **$x$-simple** and **$y$-simple** are defined in the same way with the roles of the coordinates permuted. $V$ is **elementary** if it is simple in all three directions (balls, boxes, tetrahedra, any convex solid with piecewise-smooth boundary).

> [!theorem] Theorem B.4 — Divergence theorem
> Let $V$ be a closed bounded solid that can be cut into finitely many elementary solids, with piecewise-smooth boundary $\partial V$ oriented by the outward normal. Let $\mathbf{F}=(F_1,F_2,F_3)$ have continuous first partial derivatives on an open set containing $V$. Then
> $$
> \oint_{\partial V}\mathbf{F}\cdot d\mathbf{S}=\iiint_V\nabla\cdot\mathbf{F}\,dV,\qquad \nabla\cdot\mathbf{F}=\partial_1F_1+\partial_2F_2+\partial_3F_3 .
> $$

**Proof.**

*Step 1 — a box.* Let $V=[a_1,b_1]\times[a_2,b_2]\times[a_3,b_3]$. Take the third term first. By Fubini and the fundamental theorem of calculus in $z$,

$$
\iiint_V\frac{\partial F_3}{\partial z}\,dV=\int_{a_1}^{b_1}\!\!\int_{a_2}^{b_2}\Big[F_3(x,y,b_3)-F_3(x,y,a_3)\Big]dy\,dx .
$$

Now compute the flux of the single component $F_3\hat{\mathbf{z}}$ out of the six faces. On the top face $z=b_3$ the outward normal is $+\hat{\mathbf{z}}$ and $dS=dx\,dy$, giving $\iint F_3(x,y,b_3)\,dx\,dy$; on the bottom face $z=a_3$ it is $-\hat{\mathbf{z}}$, giving $-\iint F_3(x,y,a_3)\,dx\,dy$; on the four side faces the normal is $\pm\hat{\mathbf{x}}$ or $\pm\hat{\mathbf{y}}$, perpendicular to $\hat{\mathbf{z}}$, so they contribute nothing. The two expressions agree. The same computation in $x$ and in $y$ gives the corresponding identities for $F_1\hat{\mathbf{x}}$ and $F_2\hat{\mathbf{y}}$; adding the three proves the theorem for the box. This case already contains the physics: the net outflow in the $z$ direction is the difference between what leaves through the top and what enters through the bottom, and $\partial_zF_3$ is that difference per unit volume.

*Step 2 — the $F_3$ term on a $z$-simple solid.* Let $V=\{(x,y)\in D_{xy},\ g_1\le z\le g_2\}$. As in Step 1,

$$
\iiint_V\frac{\partial F_3}{\partial z}\,dV=\iint_{D_{xy}}\Big[F_3\big(x,y,g_2(x,y)\big)-F_3\big(x,y,g_1(x,y)\big)\Big]dA .
$$

The boundary $\partial V$ has three parts.

- *Top* $S_2$: parametrise by $\mathbf{r}(x,y)=(x,y,g_2(x,y))$, $(x,y)\in D_{xy}$. Then $\mathbf{r}_x\times\mathbf{r}_y=(1,0,g_{2,x})\times(0,1,g_{2,y})=(-g_{2,x},\,-g_{2,y},\,1)$, whose $z$-component is positive, so it points upward — out of $V$. Hence $\iint_{S_2}F_3\hat{\mathbf{z}}\cdot d\mathbf{S}=\iint_{D_{xy}}F_3(x,y,g_2)\,(\hat{\mathbf{z}}\cdot(-g_{2,x},-g_{2,y},1))\,dA=\iint_{D_{xy}}F_3(x,y,g_2)\,dA$.
- *Bottom* $S_1$: the same parametrisation with $g_1$ gives an upward normal, but outward is now downward, so the sign flips: $\iint_{S_1}F_3\hat{\mathbf{z}}\cdot d\mathbf{S}=-\iint_{D_{xy}}F_3(x,y,g_1)\,dA$.
- *Lateral wall* $S_3$ (the vertical surface over $\partial D_{xy}$, possibly empty): it is made of vertical lines, so its normal is horizontal and $F_3\hat{\mathbf{z}}\cdot\hat{\mathbf{n}}=0$.

Adding the three parts reproduces the volume integral: $\oint_{\partial V}F_3\hat{\mathbf{z}}\cdot d\mathbf{S}=\iiint_V\partial_zF_3\,dV$. Note how cleanly the tilt of the top surface disappears: the factor $\sqrt{1+g_x^2+g_y^2}$ in $dS$ is exactly cancelled by the factor $1/\sqrt{1+g_x^2+g_y^2}$ in $\hat{\mathbf{n}}\cdot\hat{\mathbf{z}}$ — the flux through a tilted patch equals the flux through its shadow on the $xy$-plane.

*Step 3 — elementary solids.* If $V$ is $x$-, $y$- and $z$-simple, Step 2 and its two permuted versions give

$$
\oint_{\partial V}F_1\hat{\mathbf{x}}\cdot d\mathbf{S}=\iiint_V\partial_xF_1\,dV,\quad
\oint_{\partial V}F_2\hat{\mathbf{y}}\cdot d\mathbf{S}=\iiint_V\partial_yF_2\,dV,\quad
\oint_{\partial V}F_3\hat{\mathbf{z}}\cdot d\mathbf{S}=\iiint_V\partial_zF_3\,dV,
$$

and adding the three proves the theorem for $V$.

*Step 4 — cutting.* If $V=V_1\cup\dots\cup V_m$ with elementary pieces meeting only along (piecewise-smooth) faces, add the theorem for each piece. The volume integrals add to $\iiint_V\nabla\cdot\mathbf{F}\,dV$ since the faces have zero volume. Each internal face belongs to exactly two pieces, and the outward normal of one is the inward normal of the other, so its two flux contributions are equal and opposite and cancel. The remaining faces make up $\partial V$, each with the outward normal of $V$. ∎

```desmos-graph
left=-1.6; right=1.6; bottom=-1.3; top=1.3
---
x^2+y^2=1|blue
y=0|grey|dashed
x=0|grey|dashed
y=\sqrt{1-x^2}|red|-1<x<1
```
*Cross-section of the unit ball, an elementary solid: it is $z$-simple with $g_{1,2}=\mp\sqrt{1-x^2-y^2}$ (top hemisphere red). Cuts along coordinate planes (grey) split non-elementary solids, such as a hollow shell, into elementary pieces; the flux through each cut cancels in pairs. For $\mathbf{F}=(x^3,y^3,z^3)$ both sides of the theorem on the unit ball equal $12\pi/5$ (checked symbolically).*

**Remarks.**

- *Why $C^1$ on all of $V$.* For the point-charge field $\mathbf{F}=\hat{\mathbf{R}}/R^2$, $\nabla\cdot\mathbf{F}=\frac1{R^2}\partial_R(R^2\cdot R^{-2})=0$ for $R>0$, yet the flux through any sphere centred at the origin is $\frac{1}{a^2}\cdot4\pi a^2=4\pi$. The theorem does not apply to a ball containing the origin, where $\mathbf{F}$ is singular. Applied to the shell $\varepsilon\le R\le a$, whose boundary includes the inner sphere with **inward-pointing** normal (outward from the shell), it correctly gives $4\pi-4\pi=0$. The physicist's repair, $\nabla\cdot(\hat{\mathbf{R}}/R^2)=4\pi\delta^3(\mathbf{R})$, is discussed in [§7.11](../part-2-calculus/07-vector-calculus.md).
- *Coordinate-free divergence.* Applying the theorem to a small ball $B_\varepsilon(\mathbf{a})$ and using continuity of $\nabla\cdot\mathbf{F}$ (mean value theorem for integrals) gives $\nabla\cdot\mathbf{F}(\mathbf{a})=\lim_{\varepsilon\to0}\frac{1}{\operatorname{vol}(B_\varepsilon)}\oint_{\partial B_\varepsilon}\mathbf{F}\cdot d\mathbf{S}$: the divergence is outward flux per unit volume, independent of coordinates. This is the starting point of B.8.
- *Green's identities.* Applying the theorem to $\mathbf{F}=u\nabla w$ and using the product rule $\nabla\cdot(u\nabla w)=\nabla u\cdot\nabla w+u\nabla^2w$ (B.7) gives Green's first identity $\oint_{\partial V}u\,\partial_nw\,dS=\iiint_V(\nabla u\cdot\nabla w+u\nabla^2w)\,dV$, from which the uniqueness of solutions of Laplace's and Poisson's equations follows ([§7.11](../part-2-calculus/07-vector-calculus.md), [§13.6](../part-4-differential-equations/13-partial-differential-equations.md)).
- *Used in:* Gauss's law in differential form and the continuity equation (B.9); the energy method for the heat equation ([§13.4](../part-4-differential-equations/13-partial-differential-equations.md)); finite-volume discretisations ([Ch. 22](../part-7-numerical-computational/22-numerical-methods.md)).

## B.5 Stokes' theorem

Stokes' theorem lifts Green's theorem from a flat region to a curved surface in space ([§7.10](../part-2-calculus/07-vector-calculus.md)). The cleanest proof does not tile the surface at all: it pulls the whole problem back to the parameter domain, applies Green's theorem there, and recognises the integrand. All the geometry is concentrated in one algebraic lemma.

> [!definition] Parametrised surface and induced orientation
> A **smooth parametrised surface** is $S=\mathbf{r}(D)$ where $D$ is a plane region to which Green's theorem applies (B.3), $\mathbf{r}$ is one-to-one on $D$ and of class $C^2$ on an open set containing $D$, and $\mathbf{r}_u\times\mathbf{r}_v\ne\mathbf{0}$ on $D$. It is oriented by $\hat{\mathbf{n}}=\mathbf{r}_u\times\mathbf{r}_v/\lVert\mathbf{r}_u\times\mathbf{r}_v\rVert$. Its boundary curve $C=\partial S$ is $\mathbf{r}(\partial D)$, traversed in the direction inherited from the positive orientation of $\partial D$; this is the **right-hand rule**: with the thumb along $\hat{\mathbf{n}}$, the fingers curl along $C$.

> [!theorem] Lemma B.5 — The antisymmetric part of the Jacobian is a cross product
> Let $\mathbf{F}$ be $C^1$, $\mathbf{J}=D\mathbf{F}$ its Jacobian matrix ($J_{ij}=\partial_jF_i$), and $\boldsymbol{\omega}=\nabla\times\mathbf{F}$. Then for all vectors $\mathbf{a},\mathbf{b}$
> $$
> (\mathbf{J}-\mathbf{J}^\mathsf{T})\,\mathbf{a}=\boldsymbol{\omega}\times\mathbf{a},\qquad\text{hence}\qquad \mathbf{b}\cdot(\mathbf{J}\mathbf{a})-\mathbf{a}\cdot(\mathbf{J}\mathbf{b})=\boldsymbol{\omega}\cdot(\mathbf{a}\times\mathbf{b}).
> $$

*Proof of the lemma.* The first component of $(\mathbf{J}-\mathbf{J}^\mathsf{T})\mathbf{a}$ is $(\partial_2F_1-\partial_1F_2)a_2+(\partial_3F_1-\partial_1F_3)a_3=-\omega_3a_2+\omega_2a_3$, which is the first component of $\boldsymbol{\omega}\times\mathbf{a}=(\omega_2a_3-\omega_3a_2,\dots)$. The other two components follow by cyclic permutation $1\to2\to3\to1$. For the second identity, $\mathbf{a}\cdot(\mathbf{J}\mathbf{b})=\mathbf{b}\cdot(\mathbf{J}^\mathsf{T}\mathbf{a})$, so the left side is $\mathbf{b}\cdot\big((\mathbf{J}-\mathbf{J}^\mathsf{T})\mathbf{a}\big)=\mathbf{b}\cdot(\boldsymbol{\omega}\times\mathbf{a})=\boldsymbol{\omega}\cdot(\mathbf{a}\times\mathbf{b})$ by the cyclic property of the scalar triple product ([§7.2](../part-2-calculus/07-vector-calculus.md)). ∎

> [!theorem] Theorem B.5 — Stokes' theorem
> Let $S=\mathbf{r}(D)$ be a smooth parametrised surface with boundary $C$ oriented by the right-hand rule, or a finite union of such surfaces glued along boundary arcs with compatible orientations. If $\mathbf{F}$ has continuous first partial derivatives on an open set containing $S$, then
> $$
> \oint_C\mathbf{F}\cdot d\mathbf{r}=\iint_S(\nabla\times\mathbf{F})\cdot d\mathbf{S}.
> $$

**Proof.**

*Step 1 — pull the line integral back to $D$.* Let $(u(t),v(t))$ trace $\partial D$ positively. Then $\mathbf{r}(u(t),v(t))$ traces $C$, and by the chain rule $\frac{d}{dt}\mathbf{r}=\mathbf{r}_u\dot u+\mathbf{r}_v\dot v$. Hence

$$
\oint_C\mathbf{F}\cdot d\mathbf{r}=\oint_{\partial D}\big(\tilde P\,du+\tilde Q\,dv\big),\qquad \tilde P(u,v)=\mathbf{F}(\mathbf{r})\cdot\mathbf{r}_u,\quad\tilde Q(u,v)=\mathbf{F}(\mathbf{r})\cdot\mathbf{r}_v .
$$

Because $\mathbf{F}$ is $C^1$ and $\mathbf{r}$ is $C^2$, $\tilde P$ and $\tilde Q$ are $C^1$, so Green's theorem applies.

*Step 2 — Green's theorem in the $(u,v)$-plane.*

$$
\oint_{\partial D}\big(\tilde P\,du+\tilde Q\,dv\big)=\iint_D\left(\frac{\partial\tilde Q}{\partial u}-\frac{\partial\tilde P}{\partial v}\right)du\,dv .
$$

*Step 3 — differentiate.* By the chain rule $\partial_u\big[\mathbf{F}(\mathbf{r}(u,v))\big]=\mathbf{J}\mathbf{r}_u$ and $\partial_v\big[\mathbf{F}(\mathbf{r}(u,v))\big]=\mathbf{J}\mathbf{r}_v$, with $\mathbf{J}$ evaluated at $\mathbf{r}(u,v)$. The product rule gives

$$
\frac{\partial\tilde Q}{\partial u}=(\mathbf{J}\mathbf{r}_u)\cdot\mathbf{r}_v+\mathbf{F}\cdot\mathbf{r}_{vu},\qquad
\frac{\partial\tilde P}{\partial v}=(\mathbf{J}\mathbf{r}_v)\cdot\mathbf{r}_u+\mathbf{F}\cdot\mathbf{r}_{uv}.
$$

Since $\mathbf{r}$ is $C^2$, Clairaut's theorem gives $\mathbf{r}_{uv}=\mathbf{r}_{vu}$, so the second-derivative terms cancel on subtraction. This is the only place $C^2$ is used.

*Step 4 — recognise the curl.* By Lemma B.5 with $\mathbf{a}=\mathbf{r}_u$, $\mathbf{b}=\mathbf{r}_v$,

$$
\frac{\partial\tilde Q}{\partial u}-\frac{\partial\tilde P}{\partial v}=\mathbf{r}_v\cdot(\mathbf{J}\mathbf{r}_u)-\mathbf{r}_u\cdot(\mathbf{J}\mathbf{r}_v)=(\nabla\times\mathbf{F})\big(\mathbf{r}(u,v)\big)\cdot\big(\mathbf{r}_u\times\mathbf{r}_v\big).
$$

*Step 5 — identify the flux.* Combining Steps 1–4,

$$
\oint_C\mathbf{F}\cdot d\mathbf{r}=\iint_D(\nabla\times\mathbf{F})(\mathbf{r})\cdot(\mathbf{r}_u\times\mathbf{r}_v)\,du\,dv=\iint_S(\nabla\times\mathbf{F})\cdot d\mathbf{S},
$$

the last equality being the definition of flux with the orientation $\hat{\mathbf{n}}\parallel\mathbf{r}_u\times\mathbf{r}_v$.

*Step 6 — glued surfaces.* If $S=S_1\cup\dots\cup S_m$ with each $S_k$ a smooth parametrised surface and the orientations agreeing across shared edges, add the theorem for each piece. The flux integrals add. A shared edge is traversed in opposite directions by the induced orientations of its two neighbours (compatibility of orientation means precisely this), so its line integrals cancel, leaving the outer boundary $C$. This covers closed-up shapes such as a hemisphere parametrised in two patches, a cylinder with a lid, or a polyhedral surface. ∎

**Check.** For $\mathbf{F}=(z,x,y)$, $\nabla\times\mathbf{F}=(1,1,1)$. On the upper unit hemisphere with outward normal, the boundary is the unit circle in the $xy$-plane, counter-clockwise; there $\mathbf{F}=(0,\cos t,\sin t)$, $d\mathbf{r}=(-\sin t,\cos t,0)\,dt$, and $\oint\mathbf{F}\cdot d\mathbf{r}=\int_0^{2\pi}\cos^2t\,dt=\pi$. The flux of the constant vector $(1,1,1)$ through the hemisphere equals its flux through the flat unit disc (both surfaces share the boundary — see the remark below), which is $1\cdot\pi=\pi$. ✓

**Remarks.**

- *Green is the special case.* With $\mathbf{r}(u,v)=(u,v,0)$ we get $\mathbf{r}_u\times\mathbf{r}_v=\hat{\mathbf{z}}$ and $(\nabla\times\mathbf{F})\cdot\hat{\mathbf{z}}=\partial_xF_y-\partial_yF_x$: Theorem B.3. The graph case $z=g(x,y)$ sketched in [§7.10](../part-2-calculus/07-vector-calculus.md) is $\mathbf{r}(x,y)=(x,y,g(x,y))$, $\mathbf{r}_x\times\mathbf{r}_y=(-g_x,-g_y,1)$.
- *Orientation is a hypothesis.* On a non-orientable surface such as the Möbius strip no continuous $\hat{\mathbf{n}}$ exists and the theorem has no meaning; on an orientable one, reversing $\hat{\mathbf{n}}$ reverses $C$ and both sides change sign together.
- *Independence of the surface.* If $S$ and $S'$ have the same boundary $C$ with the same induced orientation, the theorem gives the same flux of $\nabla\times\mathbf{F}$ through both. This is consistent with B.4: $S\cup(-S')$ is a closed surface, and the outward flux of a curl through a closed surface is $\iiint\nabla\cdot(\nabla\times\mathbf{F})\,dV=0$ by B.7. It is exactly the consistency requirement that forced Maxwell's displacement current (B.9).
- *Coordinate-free curl.* For a small flat disc of radius $\varepsilon$ centred at $\mathbf{a}$ with normal $\hat{\mathbf{n}}$, continuity of $\nabla\times\mathbf{F}$ gives $(\nabla\times\mathbf{F})(\mathbf{a})\cdot\hat{\mathbf{n}}=\lim_{\varepsilon\to0}\frac{1}{\pi\varepsilon^2}\oint\mathbf{F}\cdot d\mathbf{r}$: curl is circulation per unit area, as in [§7.5](../part-2-calculus/07-vector-calculus.md).
- *Used in:* Faraday's and Ampère's laws in differential form (B.9); the curl test for conservative fields (B.6); EMF in moving and stationary loops ([§7.12](../part-2-calculus/07-vector-calculus.md)).

## B.6 Curl-free versus conservative fields

[§7.6](../part-2-calculus/07-vector-calculus.md) states that a field with zero curl is the gradient of a potential *if the domain is simply connected*. This section proves the equivalence of the three characterisations of a conservative field, proves the curl test on star-shaped and then simply connected domains, and shows by the standard counterexample that the topological hypothesis cannot be removed.

> [!definition] Conservative, star-shaped, simply connected
> Let $D\subset\mathbb{R}^n$ ($n=2$ or $3$) be open and connected. A continuous field $\mathbf{F}$ on $D$ is **conservative** if $\mathbf{F}=\nabla f$ for some $C^1$ function $f$ on $D$ (a **potential**; in electrostatics the potential is $V=-f$). $D$ is **star-shaped** about $\mathbf{a}\in D$ if the segment from $\mathbf{a}$ to every $\mathbf{x}\in D$ lies in $D$ (convex sets, all of $\mathbb{R}^3$, a half-space). $D$ is **simply connected** if every closed curve in $D$ can be continuously shrunk to a point without leaving $D$: there is a continuous $\mathbf{H}:[0,1]\times[0,1]\to D$ with $\mathbf{H}(s,0)=\boldsymbol{\gamma}(s)$, $\mathbf{H}(s,1)=$ a constant point, and $\mathbf{H}(0,t)=\mathbf{H}(1,t)$ for all $t$.

> [!theorem] Theorem B.6 — Conservative fields and the curl test
> Let $\mathbf{F}$ be continuous on an open connected set $D\subset\mathbb{R}^3$ (for $D\subset\mathbb{R}^2$ read $\nabla\times\mathbf{F}$ as the scalar $\partial_xF_y-\partial_yF_x$).
> 1. The following are equivalent: (a) $\mathbf{F}$ is conservative; (b) $\int_C\mathbf{F}\cdot d\mathbf{r}$ depends only on the endpoints of the piecewise-smooth curve $C\subset D$; (c) $\oint_C\mathbf{F}\cdot d\mathbf{r}=0$ for every piecewise-smooth closed curve $C\subset D$.
> 2. If $\mathbf{F}$ is $C^1$ and conservative, then $\nabla\times\mathbf{F}=\mathbf{0}$ on $D$.
> 3. (**Poincaré lemma**) If $D$ is star-shaped about $\mathbf{a}$, $\mathbf{F}$ is $C^1$ and $\nabla\times\mathbf{F}=\mathbf{0}$ on $D$, then $\mathbf{F}=\nabla f$ with
> $$
> f(\mathbf{x})=\int_0^1\mathbf{F}\big(\mathbf{a}+t(\mathbf{x}-\mathbf{a})\big)\cdot(\mathbf{x}-\mathbf{a})\,dt .
> $$
> 4. If $D$ is simply connected, $\mathbf{F}$ is $C^1$ and $\nabla\times\mathbf{F}=\mathbf{0}$ on $D$, then $\mathbf{F}$ is conservative.

**Proof.**

*Part 1, (a)⇒(b).* If $\mathbf{F}=\nabla f$ and $C$ is parametrised by $\mathbf{r}(t)$, $\alpha\le t\le\beta$, then by the chain rule $\frac{d}{dt}f(\mathbf{r}(t))=\nabla f\cdot\mathbf{r}'(t)$ on each smooth piece, and the fundamental theorem of calculus (summed over the pieces, the intermediate values telescoping) gives $\int_C\mathbf{F}\cdot d\mathbf{r}=f(\mathbf{r}(\beta))-f(\mathbf{r}(\alpha))$.

*(b)⇒(c).* A closed curve starts and ends at the same point, so it has the same endpoints as the constant path, whose integral is $0$.

*(c)⇒(b).* If $C_1,C_2$ run from $\mathbf{p}$ to $\mathbf{q}$, then $C_1$ followed by $C_2$ reversed is closed, so $0=\int_{C_1}-\int_{C_2}$.

*(b)⇒(a).* Fix $\mathbf{p}_0\in D$ and define $f(\mathbf{x})=\int_{\mathbf{p}_0}^{\mathbf{x}}\mathbf{F}\cdot d\mathbf{r}$ along any piecewise-smooth path in $D$; one exists because an open connected set is path-connected by polygonal paths, and the value is path-independent by (b). For $\mathbf{x}\in D$ choose $\delta>0$ with the ball $B_\delta(\mathbf{x})\subset D$. For $0<\lvert h\rvert<\delta$, reach $\mathbf{x}+h\hat{\mathbf{x}}$ by going first to $\mathbf{x}$ and then along the straight segment $\mathbf{x}+s\hat{\mathbf{x}}$, $0\le s\le h$:

$$
\frac{f(\mathbf{x}+h\hat{\mathbf{x}})-f(\mathbf{x})}{h}=\frac1h\int_0^hF_1(\mathbf{x}+s\hat{\mathbf{x}})\,ds\ \xrightarrow{h\to0}\ F_1(\mathbf{x})
$$

by continuity of $F_1$. So $\partial_xf=F_1$, and likewise $\partial_yf=F_2$, $\partial_zf=F_3$. The partials are continuous, so $f$ is $C^1$ and $\nabla f=\mathbf{F}$.

*Part 2.* If $\mathbf{F}=\nabla f$ with $\mathbf{F}$ of class $C^1$, then $f$ is $C^2$ and $(\nabla\times\nabla f)_3=\partial_x\partial_yf-\partial_y\partial_xf=0$ by Clairaut; similarly for the other components (see also B.7).

*Part 3.* Take $\mathbf{a}=\mathbf{0}$ (translate). Put $\mathbf{y}=t\mathbf{x}$. Since $\mathbf{F}$ is $C^1$ and $[0,1]$ is compact, we may differentiate under the integral sign:

$$
\partial_if(\mathbf{x})=\int_0^1\Big[F_i(t\mathbf{x})+\sum_jx_j\,t\,(\partial_iF_j)(t\mathbf{x})\Big]dt .
$$

The curl condition says $\partial_iF_j=\partial_jF_i$ for all $i,j$ (each such equation is one component of $\nabla\times\mathbf{F}=\mathbf{0}$). Hence $\sum_jx_j\,t\,\partial_iF_j(t\mathbf{x})=t\sum_jx_j\,\partial_jF_i(t\mathbf{x})=t\,\frac{d}{dt}\big[F_i(t\mathbf{x})\big]$, by the chain rule. Therefore

$$
\partial_if(\mathbf{x})=\int_0^1\frac{d}{dt}\Big[t\,F_i(t\mathbf{x})\Big]dt=1\cdot F_i(\mathbf{x})-0=F_i(\mathbf{x}),
$$

using the product rule $\frac{d}{dt}[tF_i]=F_i+t\frac{d}{dt}F_i$. Star-shapedness was used exactly once: to guarantee $t\mathbf{x}\in D$ for $t\in[0,1]$.

*Part 4.* By Part 1 it suffices to show $\oint_{\boldsymbol\gamma}\mathbf{F}\cdot d\mathbf{r}=0$ for every closed piecewise-smooth $\boldsymbol\gamma$ in $D$. Let $\mathbf{H}$ be a contraction of $\boldsymbol\gamma$ as in the definition. A continuous homotopy inside an open set can be replaced by a $C^2$ one with the same bottom edge, staying inside $D$ (a standard smoothing fact, which we assume; it is proved by convolving $\mathbf{H}$ with a smooth bump and using that $\mathbf{H}([0,1]^2)$ is a compact subset of the open set $D$). With $\mathbf{H}$ of class $C^2$, Steps 1–4 of the proof of B.5 apply verbatim to $\mathbf{r}=\mathbf{H}$ on the square $[0,1]^2$ — they used neither injectivity nor $\mathbf{r}_u\times\mathbf{r}_v\ne\mathbf{0}$ — and give

$$
\oint_{\partial[0,1]^2}\mathbf{F}(\mathbf{H})\cdot\big(\mathbf{H}_s\,ds+\mathbf{H}_t\,dt\big)=\iint_{[0,1]^2}(\nabla\times\mathbf{F})(\mathbf{H})\cdot(\mathbf{H}_s\times\mathbf{H}_t)\,ds\,dt=0 .
$$

Traverse the boundary of the square counter-clockwise. The bottom edge $t=0$ gives $\oint_{\boldsymbol\gamma}\mathbf{F}\cdot d\mathbf{r}$. The top edge $t=1$ is mapped to a single point, so $\mathbf{H}_s=\mathbf{0}$ there and it contributes $0$. The right edge $s=1$ (upward) and the left edge $s=0$ (downward) are the same path $t\mapsto\mathbf{H}(0,t)=\mathbf{H}(1,t)$ traversed in opposite directions, so they cancel. Hence $\oint_{\boldsymbol\gamma}\mathbf{F}\cdot d\mathbf{r}=0$. ∎

> [!example] The counterexample on the punctured plane
> On $D=\mathbb{R}^2\setminus\{\mathbf{0}\}$ let
> $$
> \mathbf{F}(x,y)=\frac{-y\,\hat{\mathbf{x}}+x\,\hat{\mathbf{y}}}{x^2+y^2}=\frac{1}{r}\hat{\boldsymbol\phi}.
> $$
> **Curl.** $\partial_x\dfrac{x}{x^2+y^2}=\dfrac{y^2-x^2}{(x^2+y^2)^2}$ and $\partial_y\dfrac{-y}{x^2+y^2}=\dfrac{y^2-x^2}{(x^2+y^2)^2}$, so $\partial_xF_y-\partial_yF_x=0$ everywhere on $D$; $\mathbf{F}$ is $C^\infty$ there.
>
> **Circulation.** On the unit circle $\mathbf{r}=(\cos t,\sin t)$, $\mathbf{F}=(-\sin t,\cos t)=\mathbf{r}'(t)$, so $\oint\mathbf{F}\cdot d\mathbf{r}=\int_0^{2\pi}1\,dt=2\pi\ne0$. By Part 1, $\mathbf{F}$ is **not** conservative on $D$.
>
> **What goes wrong.** $D$ is not simply connected: the unit circle cannot be shrunk to a point without passing through the origin. Locally $\mathbf{F}=\nabla\phi$ with $\phi$ the polar angle, and on the slit plane $D\setminus\{(x,0):x\le0\}$ — which is star-shaped about $(1,0)$ — $\phi\in(-\pi,\pi)$ is a genuine potential, as Part 3 guarantees. But $\phi$ cannot be extended continuously across the slit: it jumps by $2\pi$. For any closed curve in $D$, the circulation is $2\pi$ times its winding number about the origin; a loop that does not encircle the origin has zero circulation.

```desmos-graph
left=-2.5; right=3; bottom=-2; top=2
---
x^2+y^2=1|blue
(x-1.8)^2+y^2=0.16|green
y=0|red|dashed|x<0
(0,0)|label:singularity
```
*The punctured plane. Circulation of $\hat{\boldsymbol\phi}/r$ round the blue circle is $2\pi$; round the green circle, which does not enclose the origin, it is $0$. Removing the red slit makes the domain star-shaped and the polar angle a single-valued potential.*

**Remarks.**

- *Electromagnetic meaning.* The field of a long straight wire on the $z$-axis is $\mathbf{H}=\dfrac{I}{2\pi r}\hat{\boldsymbol\phi}$: curl-free in the region $r>0$ but with circulation $I$ around the wire (Ampère's law). The magnetic scalar potential $\psi=-I\phi/(2\pi)$ with $\mathbf{H}=-\nabla\psi$ is therefore multivalued, a fact exploited in magnetostatics by introducing a "barrier" surface — the slit above.
- *Holes of different dimension.* $\mathbb{R}^3$ minus a line is not simply connected (loops around the line), but $\mathbb{R}^3$ minus a point **is**: every loop can be slid off the point. So a $C^1$ curl-free field on $\mathbb{R}^3\setminus\{\mathbf{0}\}$, such as the Coulomb field, is always conservative, which is why the electrostatic potential $q/(4\pi\varepsilon_0R)$ exists globally. (The analogous statement for divergence-free fields fails there: $\hat{\mathbf{R}}/R^2$ has zero divergence on $\mathbb{R}^3\setminus\{\mathbf{0}\}$ but is not a curl, since its flux through a sphere is $4\pi\ne0$ while the flux of a curl through a closed surface vanishes by B.5.)
- *Kirchhoff's voltage law.* KVL is statement (c) for $\mathbf{E}$. Faraday's law (B.9) says $\nabla\times\mathbf{E}=-\partial_t\mathbf{B}$, so KVL holds only where $\partial_t\mathbf{B}$ is negligible — see the engineering note in [§7.6](../part-2-calculus/07-vector-calculus.md).
- *Used in:* finding potentials ([§7.6](../part-2-calculus/07-vector-calculus.md)); exact differential equations, where $M\,dx+N\,dy$ is exact if $M_y=N_x$ on a simply connected domain ([Ch. 11](../part-4-differential-equations/11-ordinary-differential-equations.md)); the existence of analytic logarithms on simply connected domains and the role of the winding number in the residue theorem ([Ch. 18](../part-5-transforms/18-complex-analysis.md)).

## B.7 Vector identities by index notation

The identities of [§7.5](../part-2-calculus/07-vector-calculus.md) can each be verified by writing out components, but index notation proves them all with two tools: the Levi-Civita symbol and one contraction identity.

> [!definition] Summation convention, Kronecker delta, Levi-Civita symbol
> Indices $i,j,k,\dots$ run over $1,2,3$, and a **repeated index is summed**: $a_ib_i=\sum_{i=1}^3a_ib_i$. The **Kronecker delta** is $\delta_{ij}=1$ if $i=j$ and $0$ otherwise. The **Levi-Civita symbol** $\varepsilon_{ijk}$ is $+1$ if $(i,j,k)$ is an even permutation of $(1,2,3)$, $-1$ if odd, and $0$ if any two indices coincide. Then, in Cartesian components,
> $$
> \mathbf{a}\cdot\mathbf{b}=a_ib_i,\quad(\mathbf{a}\times\mathbf{b})_i=\varepsilon_{ijk}a_jb_k,\quad\nabla\cdot\mathbf{F}=\partial_iF_i,\quad(\nabla\times\mathbf{F})_i=\varepsilon_{ijk}\partial_jF_k,\quad\nabla^2f=\partial_i\partial_if .
> $$

(For instance $(\mathbf{a}\times\mathbf{b})_1=\varepsilon_{123}a_2b_3+\varepsilon_{132}a_3b_2=a_2b_3-a_3b_2$.) Three properties of $\varepsilon_{ijk}$ are used repeatedly: it is **cyclic**, $\varepsilon_{ijk}=\varepsilon_{jki}=\varepsilon_{kij}$; it is **antisymmetric** under swapping any two indices; and it annihilates symmetric objects: if $S_{jk}=S_{kj}$ then $\varepsilon_{ijk}S_{jk}=0$ (rename $j\leftrightarrow k$ to get $\varepsilon_{ikj}S_{kj}=-\varepsilon_{ijk}S_{jk}$, so the sum equals its own negative).

> [!theorem] Lemma B.7 — The $\varepsilon$–$\delta$ identity
> $$
> \varepsilon_{kij}\,\varepsilon_{kmn}=\delta_{im}\delta_{jn}-\delta_{in}\delta_{jm}.
> $$

*Proof.* Fix $i,j,m,n$. On the left, a term $\varepsilon_{kij}\varepsilon_{kmn}$ is non-zero only if $i\ne j$, $m\ne n$, and $k$ is the third index different from both $i,j$ and from both $m,n$; so the left side vanishes unless $i\ne j$ and $\{m,n\}=\{i,j\}$, and then exactly one $k$ contributes. If $(m,n)=(i,j)$ the two symbols are equal and the product is $+1$; if $(m,n)=(j,i)$ they differ by one swap and the product is $-1$. On the right: if $i\ne j$ and $(m,n)=(i,j)$ we get $1-0=1$; if $(m,n)=(j,i)$, $0-1=-1$; if $i=j$ the two terms are $\delta_{im}\delta_{in}-\delta_{in}\delta_{im}=0$; and if $\{m,n\}\ne\{i,j\}$ with $i\ne j$, some index of $\{i,j\}$ is missing from $\{m,n\}$, which kills both terms. The two sides agree in every case. ∎

Contracting once more gives $\varepsilon_{kij}\varepsilon_{kmj}=3\delta_{im}-\delta_{im}=2\delta_{im}$; the BAC–CAB rule follows at once: $[\mathbf{a}\times(\mathbf{b}\times\mathbf{c})]_i=\varepsilon_{ijk}a_j\varepsilon_{kmn}b_mc_n=\varepsilon_{kij}\varepsilon_{kmn}a_jb_mc_n=b_i(a_jc_j)-c_i(a_jb_j)$, i.e. $\mathbf{a}\times(\mathbf{b}\times\mathbf{c})=\mathbf{b}(\mathbf{a}\cdot\mathbf{c})-\mathbf{c}(\mathbf{a}\cdot\mathbf{b})$.

> [!theorem] Theorem B.7 — Vector-calculus identities
> Let $f,g$ be scalar and $\mathbf{A},\mathbf{B}$ vector fields on an open set in $\mathbb{R}^3$. If they are $C^1$:
> 1. $\nabla(fg)=f\nabla g+g\nabla f$;
> 2. $\nabla\cdot(f\mathbf{A})=f\,\nabla\cdot\mathbf{A}+\mathbf{A}\cdot\nabla f$;
> 3. $\nabla\times(f\mathbf{A})=f\,\nabla\times\mathbf{A}+\nabla f\times\mathbf{A}$;
> 4. $\nabla\cdot(\mathbf{A}\times\mathbf{B})=\mathbf{B}\cdot(\nabla\times\mathbf{A})-\mathbf{A}\cdot(\nabla\times\mathbf{B})$;
> 5. $\nabla\times(\mathbf{A}\times\mathbf{B})=\mathbf{A}(\nabla\cdot\mathbf{B})-\mathbf{B}(\nabla\cdot\mathbf{A})+(\mathbf{B}\cdot\nabla)\mathbf{A}-(\mathbf{A}\cdot\nabla)\mathbf{B}$;
> 6. $\nabla(\mathbf{A}\cdot\mathbf{B})=(\mathbf{A}\cdot\nabla)\mathbf{B}+(\mathbf{B}\cdot\nabla)\mathbf{A}+\mathbf{A}\times(\nabla\times\mathbf{B})+\mathbf{B}\times(\nabla\times\mathbf{A})$.
>
> If they are $C^2$:
>
> 7. $\nabla\times(\nabla f)=\mathbf{0}$;
> 8. $\nabla\cdot(\nabla\times\mathbf{A})=0$;
> 9. $\nabla\times(\nabla\times\mathbf{A})=\nabla(\nabla\cdot\mathbf{A})-\nabla^2\mathbf{A}$, where $\nabla^2\mathbf{A}$ is the vector whose Cartesian components are $\nabla^2A_i$.
>
> Here $(\mathbf{B}\cdot\nabla)\mathbf{A}$ has components $B_j\partial_jA_i$.

**Proof.** Every step is the ordinary product rule $\partial_j(uv)=(\partial_ju)v+u\,\partial_jv$ for scalar components, plus the algebra of $\varepsilon$ and $\delta$.

1. $\partial_i(fg)=f\,\partial_ig+g\,\partial_if$.
2. $\partial_i(fA_i)=f\,\partial_iA_i+A_i\,\partial_if$.
3. $\varepsilon_{ijk}\partial_j(fA_k)=f\,\varepsilon_{ijk}\partial_jA_k+\varepsilon_{ijk}(\partial_jf)A_k$; the second term is $(\nabla f\times\mathbf{A})_i$.
4. $\partial_i(\varepsilon_{ijk}A_jB_k)=\varepsilon_{ijk}(\partial_iA_j)B_k+\varepsilon_{ijk}A_j(\partial_iB_k)$. By cyclicity $\varepsilon_{ijk}=\varepsilon_{kij}$, the first term is $B_k\,\varepsilon_{kij}\partial_iA_j=\mathbf{B}\cdot(\nabla\times\mathbf{A})$. By one swap $\varepsilon_{ijk}=-\varepsilon_{jik}$, the second is $-A_j\,\varepsilon_{jik}\partial_iB_k=-\mathbf{A}\cdot(\nabla\times\mathbf{B})$.
5. $[\nabla\times(\mathbf{A}\times\mathbf{B})]_i=\varepsilon_{ijk}\partial_j(\varepsilon_{kmn}A_mB_n)=\varepsilon_{kij}\varepsilon_{kmn}\partial_j(A_mB_n)=(\delta_{im}\delta_{jn}-\delta_{in}\delta_{jm})\big[(\partial_jA_m)B_n+A_m\partial_jB_n\big]$. The four resulting terms are $B_j\partial_jA_i+A_i\partial_jB_j-B_i\partial_jA_j-A_j\partial_jB_i$, which is the $i$-th component of the right side.
6. Compute $[\mathbf{A}\times(\nabla\times\mathbf{B})]_i=\varepsilon_{kij}\varepsilon_{kmn}A_j\partial_mB_n=A_j\partial_iB_j-A_j\partial_jB_i$ by Lemma B.7. Hence $A_j\partial_iB_j=[\mathbf{A}\times(\nabla\times\mathbf{B})]_i+[(\mathbf{A}\cdot\nabla)\mathbf{B}]_i$. Add the same identity with $\mathbf{A}\leftrightarrow\mathbf{B}$; the left sides add to $\partial_i(A_jB_j)$.
7. $(\nabla\times\nabla f)_i=\varepsilon_{ijk}\partial_j\partial_kf=0$, because $\partial_j\partial_kf$ is symmetric in $j,k$ (Clairaut, using $C^2$) and $\varepsilon_{ijk}$ annihilates symmetric objects.
8. $\partial_i(\varepsilon_{ijk}\partial_jA_k)=\varepsilon_{ijk}\,\partial_i\partial_jA_k=0$, for the same reason with the symmetric pair $(i,j)$.
9. $[\nabla\times(\nabla\times\mathbf{A})]_i=\varepsilon_{ijk}\partial_j(\varepsilon_{kmn}\partial_mA_n)=\varepsilon_{kij}\varepsilon_{kmn}\partial_j\partial_mA_n=(\delta_{im}\delta_{jn}-\delta_{in}\delta_{jm})\partial_j\partial_mA_n=\partial_i(\partial_jA_j)-\partial_j\partial_jA_i$. ∎

All nine were also checked symbolically for arbitrary component functions.

**Remarks.**

- *Smoothness matters for 7–9.* They rest on Clairaut's theorem. At a singularity they can fail in the distributional sense: the polar angle $\phi$ of B.6 has $\nabla\phi=\hat{\boldsymbol\phi}/r$, whose curl is zero for $r>0$ but equals $2\pi\delta(x)\delta(y)\hat{\mathbf{z}}$ as a distribution (the circulation $2\pi$ concentrated on the axis). Likewise $\nabla^2(1/R)=-4\pi\delta^3(\mathbf{R})$.
- *Identity 9 in curvilinear coordinates.* Only in Cartesian coordinates is $\nabla^2\mathbf{A}$ obtained by applying the scalar Laplacian to each component, because only there are the unit vectors constant. In cylindrical or spherical coordinates the **vector Laplacian is defined** by identity 9, $\nabla^2\mathbf{A}:=\nabla(\nabla\cdot\mathbf{A})-\nabla\times(\nabla\times\mathbf{A})$, evaluated with the formulas of B.8; e.g. $(\nabla^2\mathbf{A})_r=\nabla^2A_r-A_r/r^2-(2/r^2)\partial_\phi A_\phi$, not $\nabla^2A_r$.
- *Identity 5 is not symmetric in appearance* but is antisymmetric under $\mathbf{A}\leftrightarrow\mathbf{B}$, as it must be: a quick check catches sign errors.
- *Used in:* the wave equation (identity 9, B.9); the necessity of the displacement current (identity 8, B.9); Poynting's theorem, $\nabla\cdot(\mathbf{E}\times\mathbf{H})=\mathbf{H}\cdot(\nabla\times\mathbf{E})-\mathbf{E}\cdot(\nabla\times\mathbf{H})$ (identity 4, [§7.12](../part-2-calculus/07-vector-calculus.md)); Green's identities (identity 2, B.4); potentials $\mathbf{B}=\nabla\times\mathbf{A}$, $\mathbf{E}=-\nabla V-\partial_t\mathbf{A}$ (identities 7–8).

## B.8 Gradient, divergence, curl and Laplacian in curvilinear coordinates

The cylindrical and spherical formulas tabulated in [§7.5](../part-2-calculus/07-vector-calculus.md) and used throughout [Ch. 13](../part-4-differential-equations/13-partial-differential-equations.md) are often "derived" by drawing a small curved box and estimating fluxes. Here they are derived exactly, for every orthogonal coordinate system at once, using only the chain rule and the identities of B.7. The only geometric input is the set of three scale factors.

> [!definition] Orthogonal curvilinear coordinates and scale factors
> Let $\mathbf{x}=\mathbf{x}(u_1,u_2,u_3)$ be a one-to-one $C^2$ map from an open set of $(u_1,u_2,u_3)$-space onto an open set $\Omega\subset\mathbb{R}^3$, with non-zero Jacobian determinant. The coordinates are **orthogonal** if the three tangent vectors $\partial\mathbf{x}/\partial u_i$ are mutually perpendicular at every point. The **scale factors** and **unit vectors** are
> $$
> h_i=\left\lVert\frac{\partial\mathbf{x}}{\partial u_i}\right\rVert,\qquad\hat{\mathbf{u}}_i=\frac{1}{h_i}\frac{\partial\mathbf{x}}{\partial u_i}\qquad(i=1,2,3,\ \text{no sum}),
> $$
> and we assume the system is right-handed: $\hat{\mathbf{u}}_1\times\hat{\mathbf{u}}_2=\hat{\mathbf{u}}_3$ (and cyclically). A field is written $\mathbf{F}=F_1\hat{\mathbf{u}}_1+F_2\hat{\mathbf{u}}_2+F_3\hat{\mathbf{u}}_3$, and $\partial_i$ now means $\partial/\partial u_i$.

For cylindrical coordinates $(u_1,u_2,u_3)=(r,\phi,z)$, $\mathbf{x}=(r\cos\phi,r\sin\phi,z)$:

$$
\frac{\partial\mathbf{x}}{\partial r}=(\cos\phi,\sin\phi,0),\quad\frac{\partial\mathbf{x}}{\partial\phi}=r(-\sin\phi,\cos\phi,0),\quad\frac{\partial\mathbf{x}}{\partial z}=(0,0,1),\qquad (h_r,h_\phi,h_z)=(1,\,r,\,1).
$$

For spherical coordinates $(u_1,u_2,u_3)=(R,\theta,\phi)$, $\mathbf{x}=(R\sin\theta\cos\phi,\,R\sin\theta\sin\phi,\,R\cos\theta)$:

$$
\frac{\partial\mathbf{x}}{\partial R}=(\sin\theta\cos\phi,\sin\theta\sin\phi,\cos\theta),\quad
\frac{\partial\mathbf{x}}{\partial\theta}=R(\cos\theta\cos\phi,\cos\theta\sin\phi,-\sin\theta),\quad
\frac{\partial\mathbf{x}}{\partial\phi}=R\sin\theta(-\sin\phi,\cos\phi,0),
$$

so $(h_R,h_\theta,h_\phi)=(1,\,R,\,R\sin\theta)$. Pairwise dot products are zero by inspection, and both systems are right-handed ($\hat{\mathbf{r}}\times\hat{\boldsymbol\phi}=\hat{\mathbf{z}}$, $\hat{\mathbf{R}}\times\hat{\boldsymbol\theta}=\hat{\boldsymbol\phi}$). The maps are one-to-one with non-zero Jacobian away from the $z$-axis ($r=0$, respectively $\sin\theta=0$); the formulas below hold there.

> [!theorem] Lemma B.8 — Gradients of the coordinate functions
> Regard $u_i$ as functions of $\mathbf{x}$ (via the inverse map, which is $C^2$ by the inverse function theorem). Then
> $$
> \nabla u_i=\frac{\hat{\mathbf{u}}_i}{h_i},\qquad\text{and consequently}\qquad\hat{\mathbf{u}}_1=h_2h_3\,\nabla u_2\times\nabla u_3,\quad\hat{\mathbf{u}}_1=h_1\nabla u_1\quad(\text{and cyclically}).
> $$

*Proof.* The Jacobian matrix of $\mathbf{x}(\mathbf{u})$ has columns $h_i\hat{\mathbf{u}}_i$, so it factors as $\mathbf{Q}\operatorname{diag}(h_1,h_2,h_3)$ with $\mathbf{Q}=[\hat{\mathbf{u}}_1\ \hat{\mathbf{u}}_2\ \hat{\mathbf{u}}_3]$ orthogonal. The Jacobian of the inverse map $\mathbf{u}(\mathbf{x})$ is the inverse matrix, $\operatorname{diag}(1/h_1,1/h_2,1/h_3)\,\mathbf{Q}^\mathsf{T}$, whose $i$-th row is $\hat{\mathbf{u}}_i^\mathsf{T}/h_i$. But the $i$-th row of the Jacobian of $\mathbf{u}(\mathbf{x})$ is $(\nabla u_i)^\mathsf{T}$. The consequences follow from $\nabla u_2\times\nabla u_3=\hat{\mathbf{u}}_2\times\hat{\mathbf{u}}_3/(h_2h_3)=\hat{\mathbf{u}}_1/(h_2h_3)$. ∎

> [!theorem] Theorem B.8 — Differential operators in orthogonal coordinates
> For a $C^1$ scalar $f$ and $C^1$ field $\mathbf{F}=\sum F_i\hat{\mathbf{u}}_i$ (and $C^2$ $f$ for the Laplacian), with $h=h_1h_2h_3$:
> $$
> \begin{aligned}
> \nabla f&=\frac{1}{h_1}\frac{\partial f}{\partial u_1}\hat{\mathbf{u}}_1+\frac{1}{h_2}\frac{\partial f}{\partial u_2}\hat{\mathbf{u}}_2+\frac{1}{h_3}\frac{\partial f}{\partial u_3}\hat{\mathbf{u}}_3,\\
> \nabla\cdot\mathbf{F}&=\frac{1}{h}\left[\frac{\partial}{\partial u_1}(h_2h_3F_1)+\frac{\partial}{\partial u_2}(h_3h_1F_2)+\frac{\partial}{\partial u_3}(h_1h_2F_3)\right],\\
> \nabla\times\mathbf{F}&=\frac{1}{h}\det\begin{bmatrix}h_1\hat{\mathbf{u}}_1&h_2\hat{\mathbf{u}}_2&h_3\hat{\mathbf{u}}_3\\ \partial_1&\partial_2&\partial_3\\ h_1F_1&h_2F_2&h_3F_3\end{bmatrix},\\
> \nabla^2f&=\frac{1}{h}\left[\frac{\partial}{\partial u_1}\!\left(\frac{h_2h_3}{h_1}\frac{\partial f}{\partial u_1}\right)+\frac{\partial}{\partial u_2}\!\left(\frac{h_3h_1}{h_2}\frac{\partial f}{\partial u_2}\right)+\frac{\partial}{\partial u_3}\!\left(\frac{h_1h_2}{h_3}\frac{\partial f}{\partial u_3}\right)\right].
> \end{aligned}
> $$

**Proof.**

*Gradient.* By the chain rule, $\partial f/\partial x_j=\sum_i(\partial f/\partial u_i)(\partial u_i/\partial x_j)$, i.e. $\nabla f=\sum_i\partial_if\,\nabla u_i$. Insert Lemma B.8.

*Divergence.* It suffices (by linearity) to treat $F_1\hat{\mathbf{u}}_1$. Write it as $(h_2h_3F_1)\,(\nabla u_2\times\nabla u_3)$ using Lemma B.8. By identity 4 of B.7 and identity 7 (the $u_i$ are $C^2$),

$$
\nabla\cdot(\nabla u_2\times\nabla u_3)=\nabla u_3\cdot(\nabla\times\nabla u_2)-\nabla u_2\cdot(\nabla\times\nabla u_3)=0 .
$$

Then identity 2 of B.7 gives

$$
\nabla\cdot(F_1\hat{\mathbf{u}}_1)=\nabla(h_2h_3F_1)\cdot(\nabla u_2\times\nabla u_3)=\nabla(h_2h_3F_1)\cdot\frac{\hat{\mathbf{u}}_1}{h_2h_3}=\frac{1}{h_1h_2h_3}\frac{\partial}{\partial u_1}(h_2h_3F_1),
$$

the last step by the gradient formula (only the $\hat{\mathbf{u}}_1$-component survives the dot product). The other two terms follow cyclically.

*Curl.* Again treat $F_1\hat{\mathbf{u}}_1=(h_1F_1)\nabla u_1$. By identity 3 of B.7 and identity 7,

$$
\nabla\times\big((h_1F_1)\nabla u_1\big)=\nabla(h_1F_1)\times\nabla u_1=\left[\sum_k\frac{1}{h_k}\partial_k(h_1F_1)\hat{\mathbf{u}}_k\right]\times\frac{\hat{\mathbf{u}}_1}{h_1}
=\frac{1}{h_3h_1}\partial_3(h_1F_1)\,\hat{\mathbf{u}}_2-\frac{1}{h_1h_2}\partial_2(h_1F_1)\,\hat{\mathbf{u}}_3,
$$

using $\hat{\mathbf{u}}_3\times\hat{\mathbf{u}}_1=\hat{\mathbf{u}}_2$ and $\hat{\mathbf{u}}_2\times\hat{\mathbf{u}}_1=-\hat{\mathbf{u}}_3$. Adding the cyclic analogues for $F_2\hat{\mathbf{u}}_2$ and $F_3\hat{\mathbf{u}}_3$ and collecting the coefficient of each unit vector gives, for example,

$$
(\nabla\times\mathbf{F})_1=\frac{1}{h_2h_3}\Big[\partial_2(h_3F_3)-\partial_3(h_2F_2)\Big],
$$

which is exactly the expansion of the stated determinant along its first row.

*Laplacian.* $\nabla^2f=\nabla\cdot(\nabla f)$; apply the divergence formula to $F_i=\frac{1}{h_i}\partial_if$. ∎

Specialising (all formulas were verified symbolically against the Cartesian operators on test fields):

> [!theorem] Cylindrical coordinates $(r,\phi,z)$, $h=(1,r,1)$
> $$
> \begin{aligned}
> \nabla f&=\frac{\partial f}{\partial r}\hat{\mathbf{r}}+\frac1r\frac{\partial f}{\partial\phi}\hat{\boldsymbol\phi}+\frac{\partial f}{\partial z}\hat{\mathbf{z}},\qquad
> \nabla\cdot\mathbf{F}=\frac1r\frac{\partial}{\partial r}(rF_r)+\frac1r\frac{\partial F_\phi}{\partial\phi}+\frac{\partial F_z}{\partial z},\\
> \nabla\times\mathbf{F}&=\left(\frac1r\frac{\partial F_z}{\partial\phi}-\frac{\partial F_\phi}{\partial z}\right)\hat{\mathbf{r}}+\left(\frac{\partial F_r}{\partial z}-\frac{\partial F_z}{\partial r}\right)\hat{\boldsymbol\phi}+\frac1r\left(\frac{\partial}{\partial r}(rF_\phi)-\frac{\partial F_r}{\partial\phi}\right)\hat{\mathbf{z}},\\
> \nabla^2f&=\frac1r\frac{\partial}{\partial r}\!\left(r\frac{\partial f}{\partial r}\right)+\frac{1}{r^2}\frac{\partial^2f}{\partial\phi^2}+\frac{\partial^2f}{\partial z^2}.
> \end{aligned}
> $$

> [!theorem] Spherical coordinates $(R,\theta,\phi)$, $h=(1,R,R\sin\theta)$
> $$
> \begin{aligned}
> \nabla f&=\frac{\partial f}{\partial R}\hat{\mathbf{R}}+\frac1R\frac{\partial f}{\partial\theta}\hat{\boldsymbol\theta}+\frac{1}{R\sin\theta}\frac{\partial f}{\partial\phi}\hat{\boldsymbol\phi},\\
> \nabla\cdot\mathbf{F}&=\frac{1}{R^2}\frac{\partial}{\partial R}(R^2F_R)+\frac{1}{R\sin\theta}\frac{\partial}{\partial\theta}(\sin\theta\,F_\theta)+\frac{1}{R\sin\theta}\frac{\partial F_\phi}{\partial\phi},\\
> \nabla\times\mathbf{F}&=\frac{1}{R\sin\theta}\left(\frac{\partial}{\partial\theta}(\sin\theta\,F_\phi)-\frac{\partial F_\theta}{\partial\phi}\right)\hat{\mathbf{R}}+\frac1R\left(\frac{1}{\sin\theta}\frac{\partial F_R}{\partial\phi}-\frac{\partial}{\partial R}(RF_\phi)\right)\hat{\boldsymbol\theta}+\frac1R\left(\frac{\partial}{\partial R}(RF_\theta)-\frac{\partial F_R}{\partial\theta}\right)\hat{\boldsymbol\phi},\\
> \nabla^2f&=\frac{1}{R^2}\frac{\partial}{\partial R}\!\left(R^2\frac{\partial f}{\partial R}\right)+\frac{1}{R^2\sin\theta}\frac{\partial}{\partial\theta}\!\left(\sin\theta\frac{\partial f}{\partial\theta}\right)+\frac{1}{R^2\sin^2\theta}\frac{\partial^2f}{\partial\phi^2}.
> \end{aligned}
> $$

```desmos-graph
left=-0.5; right=2.6; bottom=-0.3; top=2.3
---
x^2+y^2=1|blue|x>0
x^2+y^2=4|blue|x>0
y=\tan(0.5)x|red|0<x<2
y=\tan(0.8)x|red|0<x<1.5
y=0|grey|0<x<2.5
```
*Cylindrical (polar) coordinate curves. The patch between $r=1,2$ and $\phi=0.5,0.8$ has radial sides of length $\Delta r$ but arcs of length $r\,\Delta\phi$ — longer on the outer circle. That growth is the scale factor $h_\phi=r$, and it is why $\partial_r(rF_r)$, not $\partial_rF_r$, appears in the divergence.*

**Remarks.**

- *Orthogonality is essential.* For non-orthogonal coordinates $\nabla u_i$ is no longer parallel to $\partial\mathbf{x}/\partial u_i$ and the formulas need the full metric tensor $g_{ij}=\partial_i\mathbf{x}\cdot\partial_j\mathbf{x}$; e.g. $\nabla\cdot\mathbf{F}=\frac{1}{\sqrt{g}}\partial_i(\sqrt{g}\,F^i)$ with contravariant components. All coordinate systems in this book are orthogonal.
- *Check by a known field.* For the Coulomb field $F_R=1/R^2$: $\nabla\cdot\mathbf{F}=\frac1{R^2}\partial_R(1)=0$ for $R>0$, consistent with B.4. For $f=1/R$: $\nabla^2f=\frac{1}{R^2}\partial_R(R^2\cdot(-R^{-2}))=0$. For the wire field $F_\phi=1/r$: $(\nabla\times\mathbf{F})_z=\frac1r\partial_r(r\cdot\frac1r)=0$, consistent with B.6.
- *A pitfall.* Do not differentiate unit vectors as if constant: $\nabla\cdot\hat{\mathbf{r}}=\frac1r\partial_r(r)=\frac1r\ne0$, even though "$F_r=1$ is constant". The scale factors encode exactly the variation of $\hat{\mathbf{u}}_i$ that a naive componentwise derivative misses.
- *Axis singularities.* On $r=0$ (or $\sin\theta=0$) the coordinates fail to be a valid chart; smooth Cartesian fields still have finite operators there, and the apparently singular terms such as $\frac{1}{r^2}\partial_\phi^2f$ combine into finite limits.
- *Used in:* separation of variables for Laplace's equation in cylinders and spheres (Bessel and Legendre equations, [§13.6](../part-4-differential-equations/13-partial-differential-equations.md)); radiation and waveguide fields ([§7.12](../part-2-calculus/07-vector-calculus.md)); the gradient in polar coordinates ([§6.4](../part-2-calculus/06-multivariable-calculus.md)).

## B.9 Maxwell's equations: from integral to differential form, continuity, and waves

[§7.12](../part-2-calculus/07-vector-calculus.md) derives the differential form of Maxwell's equations in a few lines. This section supplies the two analytical facts those lines rely on — localisation of an integral identity, and differentiation of a flux under the integral sign — and states the smoothness hypotheses under which the differential form, the continuity equation and the wave equation follow.

> [!theorem] Lemma B.9a — Localisation
> Let $\Omega\subset\mathbb{R}^3$ be open.
> 1. If $g$ is continuous on $\Omega$ and $\iiint_Bg\,dV=0$ for every closed ball $B\subset\Omega$, then $g\equiv0$ on $\Omega$.
> 2. If $\mathbf{G}$ is continuous on $\Omega$ and $\iint_S\mathbf{G}\cdot d\mathbf{S}=0$ for every flat disc $S\subset\Omega$ (of every orientation), then $\mathbf{G}\equiv\mathbf{0}$ on $\Omega$.

*Proof.* (1) Suppose $g(\mathbf{a})>0$ for some $\mathbf{a}$ (the case $<0$ is symmetric). By continuity there is $\delta>0$ with $g>g(\mathbf{a})/2$ on the ball $B_\delta(\mathbf{a})\subset\Omega$, so $\iiint_{B_\delta}g\,dV\ge\frac{g(\mathbf{a})}{2}\cdot\frac43\pi\delta^3>0$, a contradiction. (2) Suppose $\mathbf{G}(\mathbf{a})\ne\mathbf{0}$ and put $\hat{\mathbf{n}}=\mathbf{G}(\mathbf{a})/\lVert\mathbf{G}(\mathbf{a})\rVert$. The continuous function $\mathbf{G}\cdot\hat{\mathbf{n}}$ equals $\lVert\mathbf{G}(\mathbf{a})\rVert>0$ at $\mathbf{a}$, so it exceeds half that on a ball around $\mathbf{a}$; the flat disc through $\mathbf{a}$ with normal $\hat{\mathbf{n}}$ and small radius $\delta$ lies in that ball and has flux at least $\frac12\lVert\mathbf{G}(\mathbf{a})\rVert\pi\delta^2>0$, a contradiction. ∎

> [!theorem] Lemma B.9b — Differentiating a flux through a fixed surface
> Let $S=\mathbf{r}(D)$ be a fixed (time-independent) smooth parametrised surface and $\mathbf{B}(\mathbf{x},t)$ a field such that $\mathbf{B}$ and $\partial\mathbf{B}/\partial t$ are continuous on $S\times I$ for an open time interval $I$. Then
> $$
> \frac{d}{dt}\iint_S\mathbf{B}\cdot d\mathbf{S}=\iint_S\frac{\partial\mathbf{B}}{\partial t}\cdot d\mathbf{S}.
> $$
> The same holds for volume integrals $\frac{d}{dt}\iiint_V\rho\,dV=\iiint_V\partial_t\rho\,dV$ over a fixed $V$.

*Proof.* In parameters the flux is $\Phi(t)=\iint_D\mathbf{B}(\mathbf{r}(u,v),t)\cdot\mathbf{n}(u,v)\,du\,dv$ with $\mathbf{n}=\mathbf{r}_u\times\mathbf{r}_v$ independent of $t$. The integrand and its $t$-derivative $\partial_t\mathbf{B}(\mathbf{r},t)\cdot\mathbf{n}$ are continuous on the compact set $D\times[t_0-\eta,t_0+\eta]$, hence uniformly continuous. By the mean value theorem, $\big[\Phi(t_0+\tau)-\Phi(t_0)\big]/\tau=\iint_D\partial_t\mathbf{B}(\mathbf{r},t_0+\vartheta\tau)\cdot\mathbf{n}\,du\,dv$ for some $\vartheta=\vartheta(u,v)\in(0,1)$, and uniform continuity makes this converge to $\iint_D\partial_t\mathbf{B}(\mathbf{r},t_0)\cdot\mathbf{n}\,du\,dv$ as $\tau\to0$. The volume case is identical. ∎

> [!theorem] Theorem B.9 — Differential form of Maxwell's equations
> Let $\Omega\subset\mathbb{R}^3$ be open and $I$ an open time interval. Suppose $\mathbf{D},\mathbf{B},\mathbf{E},\mathbf{H}$ are $C^1$ in $(\mathbf{x},t)$ on $\Omega\times I$ and $\rho,\mathbf{J}$ are continuous. Then the integral laws
> $$
> \begin{aligned}
> \oint_{\partial V}\mathbf{D}\cdot d\mathbf{S}&=\iiint_V\rho\,dV, &\qquad \oint_{\partial V}\mathbf{B}\cdot d\mathbf{S}&=0,\\
> \oint_{\partial S}\mathbf{E}\cdot d\mathbf{l}&=-\frac{d}{dt}\iint_S\mathbf{B}\cdot d\mathbf{S}, &\qquad \oint_{\partial S}\mathbf{H}\cdot d\mathbf{l}&=\iint_S\mathbf{J}\cdot d\mathbf{S}+\frac{d}{dt}\iint_S\mathbf{D}\cdot d\mathbf{S},
> \end{aligned}
> $$
> holding for every solid $V\subset\Omega$ to which B.4 applies and every fixed oriented surface $S\subset\Omega$ to which B.5 applies, at every $t\in I$, are **equivalent** to the differential laws
> $$
> \nabla\cdot\mathbf{D}=\rho,\qquad\nabla\cdot\mathbf{B}=0,\qquad\nabla\times\mathbf{E}=-\frac{\partial\mathbf{B}}{\partial t},\qquad\nabla\times\mathbf{H}=\mathbf{J}+\frac{\partial\mathbf{D}}{\partial t}\qquad\text{on }\Omega\times I .
> $$

**Proof.**

*Integral ⇒ differential.* (Gauss.) By B.4, $\oint_{\partial V}\mathbf{D}\cdot d\mathbf{S}=\iiint_V\nabla\cdot\mathbf{D}\,dV$, so $\iiint_V(\nabla\cdot\mathbf{D}-\rho)\,dV=0$ for every ball $V\subset\Omega$. The integrand is continuous ($\mathbf{D}$ is $C^1$, $\rho$ continuous), so Lemma B.9a(1) gives $\nabla\cdot\mathbf{D}=\rho$. The same argument with $\rho=0$ gives $\nabla\cdot\mathbf{B}=0$.

(Faraday.) For a fixed surface, Lemma B.9b moves the time derivative inside, and B.5 converts the circulation:

$$
\iint_S\left(\nabla\times\mathbf{E}+\frac{\partial\mathbf{B}}{\partial t}\right)\cdot d\mathbf{S}=0\qquad\text{for every flat disc }S\subset\Omega .
$$

The bracket is continuous, so Lemma B.9a(2) gives $\nabla\times\mathbf{E}=-\partial_t\mathbf{B}$. (Ampère–Maxwell.) Identically, $\iint_S(\nabla\times\mathbf{H}-\mathbf{J}-\partial_t\mathbf{D})\cdot d\mathbf{S}=0$ for every disc, whence $\nabla\times\mathbf{H}=\mathbf{J}+\partial_t\mathbf{D}$.

*Differential ⇒ integral.* Integrate each differential law over $V$ or $S$ and apply B.4, B.5 and Lemma B.9b in the reverse direction. ∎

> [!theorem] Corollary B.9c — The continuity equation
> (a) If charge is conserved in the integral sense, $\dfrac{d}{dt}\iiint_V\rho\,dV=-\oint_{\partial V}\mathbf{J}\cdot d\mathbf{S}$ for every fixed $V\subset\Omega$, with $\rho$ of class $C^1$ and $\mathbf{J}$ of class $C^1$, then
> $$
> \frac{\partial\rho}{\partial t}+\nabla\cdot\mathbf{J}=0 .
> $$
> (b) Conversely, the continuity equation is a **consequence** of Maxwell's equations whenever $\mathbf{H}$ and $\mathbf{D}$ are $C^2$ and $\mathbf{J}$ is $C^1$.

*Proof.* (a) Lemma B.9b and B.4 give $\iiint_V(\partial_t\rho+\nabla\cdot\mathbf{J})\,dV=0$ for every ball; apply Lemma B.9a(1). (b) Take the divergence of the Ampère–Maxwell law. By identity 8 of B.7, $\nabla\cdot(\nabla\times\mathbf{H})=0$, so $0=\nabla\cdot\mathbf{J}+\nabla\cdot\partial_t\mathbf{D}$. Because $\mathbf{D}$ is $C^2$, the space and time derivatives commute (Clairaut), $\nabla\cdot\partial_t\mathbf{D}=\partial_t(\nabla\cdot\mathbf{D})=\partial_t\rho$ by Gauss's law. ∎

Part (b) is the precise form of Maxwell's argument for the displacement current ([§7.12](../part-2-calculus/07-vector-calculus.md)): without the term $\partial_t\mathbf{D}$ the same computation would force $\nabla\cdot\mathbf{J}=0$, contradicting charge accumulation on a capacitor plate. In the same way the divergence of Faraday's law gives $\partial_t(\nabla\cdot\mathbf{B})=0$, so $\nabla\cdot\mathbf{B}=0$ is preserved in time once it holds initially: the two divergence equations are constraints compatible with the two curl equations.

> [!theorem] Corollary B.9d — The wave equation in vacuum
> In a region with $\rho=0$, $\mathbf{J}=\mathbf{0}$, $\mathbf{D}=\varepsilon_0\mathbf{E}$, $\mathbf{B}=\mu_0\mathbf{H}$, and fields of class $C^2$,
> $$
> \nabla^2\mathbf{E}=\mu_0\varepsilon_0\frac{\partial^2\mathbf{E}}{\partial t^2},\qquad\nabla^2\mathbf{B}=\mu_0\varepsilon_0\frac{\partial^2\mathbf{B}}{\partial t^2},\qquad c=\frac{1}{\sqrt{\mu_0\varepsilon_0}} .
> $$

*Proof.* Take the curl of Faraday's law. Since $\mathbf{B}$ is $C^2$, $\nabla\times\partial_t\mathbf{B}=\partial_t(\nabla\times\mathbf{B})$, and by Ampère–Maxwell $\nabla\times\mathbf{B}=\mu_0\nabla\times\mathbf{H}=\mu_0\varepsilon_0\,\partial_t\mathbf{E}$. Hence

$$
\nabla\times(\nabla\times\mathbf{E})=-\partial_t(\nabla\times\mathbf{B})=-\mu_0\varepsilon_0\frac{\partial^2\mathbf{E}}{\partial t^2}.
$$

By identity 9 of B.7 (in Cartesian components), $\nabla\times(\nabla\times\mathbf{E})=\nabla(\nabla\cdot\mathbf{E})-\nabla^2\mathbf{E}=-\nabla^2\mathbf{E}$, since $\nabla\cdot\mathbf{E}=\rho/\varepsilon_0=0$. This is the first equation. For $\mathbf{B}$, take the curl of Ampère–Maxwell: $\nabla\times(\nabla\times\mathbf{B})=\mu_0\varepsilon_0\partial_t(\nabla\times\mathbf{E})=-\mu_0\varepsilon_0\partial_t^2\mathbf{B}$, and $\nabla\cdot\mathbf{B}=0$. ∎

**Numbers and a check.** With $\mu_0=4\pi\times10^{-7}$ H/m and $\varepsilon_0=8.8541878\times10^{-12}$ F/m, $1/\sqrt{\mu_0\varepsilon_0}=2.99792458\times10^8$ m/s (computed), the measured speed of light. Dimensionally, $[\mu_0\varepsilon_0]=(\mathrm{H/m})(\mathrm{F/m})=\mathrm{s^2/m^2}$ because $[LC]=\mathrm{s^2}$, so $1/\sqrt{\mu_0\varepsilon_0}$ is a speed. The plane wave $\mathbf{E}=E_0\cos(\omega t-kz)\,\hat{\mathbf{x}}$ satisfies the equation exactly when $k^2=\mu_0\varepsilon_0\omega^2$, i.e. $\omega/k=c$; it also satisfies $\nabla\cdot\mathbf{E}=\partial_xE_x=0$, and Faraday's law then gives $\mathbf{B}=(E_0/c)\cos(\omega t-kz)\,\hat{\mathbf{y}}$ — transverse and in phase, as solved in [§13.5](../part-4-differential-equations/13-partial-differential-equations.md).

**Remarks.**

- *Why the integral form is more fundamental.* The integral laws make sense for fields that are discontinuous across material interfaces, sheets of charge or surface currents, where the differential form fails (the derivatives do not exist). Applied to thin pillboxes and loops straddling an interface they give the boundary conditions $\hat{\mathbf{n}}\cdot(\mathbf{D}_1-\mathbf{D}_2)=\rho_s$, $\hat{\mathbf{n}}\times(\mathbf{E}_1-\mathbf{E}_2)=\mathbf{0}$ ([§7.12](../part-2-calculus/07-vector-calculus.md)). The equivalence in Theorem B.9 holds precisely in regions where the $C^1$ hypothesis is met.
- *Fixed versus moving surfaces.* Lemma B.9b requires $S$ to be stationary. For a moving loop, $\frac{d}{dt}\iint_{S(t)}\mathbf{B}\cdot d\mathbf{S}$ picks up an extra term $-\oint_{\partial S}(\mathbf{v}\times\mathbf{B})\cdot d\mathbf{l}$ (using $\nabla\cdot\mathbf{B}=0$); this is the motional EMF of generators. The differential law $\nabla\times\mathbf{E}=-\partial_t\mathbf{B}$ is unaffected.
- *Is the wave equation all of Maxwell?* No. Every solution of Maxwell's equations in vacuum solves the wave equation, but not conversely: a solution of $\nabla^2\mathbf{E}=c^{-2}\partial_t^2\mathbf{E}$ must also satisfy $\nabla\cdot\mathbf{E}=0$ (transversality), and $\mathbf{B}$ is then determined by Faraday's law up to a static field. Longitudinal solutions such as $\cos(\omega t-kz)\hat{\mathbf{z}}$ satisfy the wave equation but are excluded by Gauss's law.
- *Used in:* the wave equation and Helmholtz equation ([§7.12](../part-2-calculus/07-vector-calculus.md), [§13.5](../part-4-differential-equations/13-partial-differential-equations.md)); Poisson's equation $\nabla^2V=-\rho/\varepsilon$ from $\mathbf{E}=-\nabla V$, $\nabla\cdot\mathbf{D}=\rho$ ([§13.6](../part-4-differential-equations/13-partial-differential-equations.md), [§13.8](../part-4-differential-equations/13-partial-differential-equations.md)); the phasor form $\partial_t\to j\omega$ used with the Fourier transform ([Ch. 16](../part-5-transforms/16-fourier-transforms.md)).

## Further reading

- J. E. Marsden and A. Tromba, *Vector Calculus*, 6th ed., W. H. Freeman — proofs of Green, Gauss and Stokes for elementary regions, in the style used here.
- T. M. Apostol, *Calculus*, Vol. II, 2nd ed., Wiley — change of variables and the integral theorems with full rigour.
- M. Spivak, *Calculus on Manifolds*, Benjamin — the general Stokes theorem for differential forms, of which B.3–B.5 are special cases, and the Poincaré lemma.
- D. K. Cheng, *Field and Wave Electromagnetics*, 2nd ed., Addison-Wesley — curvilinear operators via scale factors and Maxwell's equations in the conventions of this book.
- H. M. Schey, *Div, Grad, Curl, and All That*, 4th ed., Norton — an informal companion that motivates the coordinate-free definitions.
