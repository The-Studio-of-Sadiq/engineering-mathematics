# Chapter 18 — Complex Analysis

[← Previous](17-z-transforms.md) · [Table of Contents](../00-front/table-of-contents.md) · [Next →](../part-6-probability-statistics/19-probability.md)

> **Prerequisites:** Ch. 2 (complex numbers), Ch. 5 (integral calculus), Ch. 6 (multivariable calculus), Ch. 7 (vector calculus: Green's theorem), Ch. 14 (Laplace transforms), Ch. 17 (Z-transforms)
> **Learning outcomes:** After this chapter you can
> - treat a function $w=f(z)$ as a mapping of the plane, compute its real and imaginary parts, and picture what it does to grid lines (conformal mapping);
> - test a function for analyticity with the Cauchy–Riemann equations and construct harmonic conjugates;
> - parametrize contours, evaluate contour integrals, and bound them with the ML inequality;
> - apply Cauchy's integral theorem and integral formula, including contour deformation;
> - expand functions in Taylor and Laurent series on discs and annuli, and classify isolated singularities;
> - compute residues and use the residue theorem to evaluate real integrals (rational, trigonometric, Fourier-type via Jordan's lemma), to invert Laplace and Z-transforms, and to count zeros and poles with the argument principle — the mathematics behind the Nyquist stability criterion.

Every transform in Part V lives in the complex plane. The Laplace variable $s=\sigma+j\omega$ of [Ch. 14](14-laplace-transforms.md), the Fourier frequency on the line $s=j\omega$ of [Ch. 16](16-fourier-transforms.md), and the $z$ of [Ch. 17](17-z-transforms.md) are complex numbers, and the functions we manipulate — transfer functions $H(s)$, spectra $X(\omega)$, system functions $H(z)$ — are complex functions of a complex variable. Until now we have treated them algebraically: partial fractions, tables, pole–zero plots. This chapter supplies the analysis underneath. It explains *why* partial-fraction coefficients are what they are (they are residues), *why* the region of convergence of a Z-transform is an annulus (it is the annulus of a Laurent series), *why* the inverse Laplace transform is an integral along a vertical line (the Bromwich integral), and *why* counting encirclements of $-1$ by a Nyquist plot reveals closed-loop poles (the argument principle).

The central surprise of complex analysis is how rigid differentiable complex functions are. A real function can be differentiable once and nowhere twice; a complex function that is differentiable once in a region is differentiable infinitely often, equals its Taylor series, and is completely determined inside a closed curve by its values on the curve. That rigidity is what makes contour integration such a powerful computational tool.

Throughout, $z = x+jy = re^{j\theta}$ denotes the independent variable and $w=f(z)=u(x,y)+jv(x,y)$ the dependent one.

## 18.1 Complex functions

### Why an engineer needs complex functions

A transfer function $H(s)=1/(s^2+2\zeta\omega_n s+\omega_n^2)$ assigns to every complex frequency $s$ a complex number: a gain and a phase. A Bode plot is a picture of $H$ restricted to the imaginary axis; a root locus is a picture of where $1+KH(s)=0$; a Nyquist plot is the *image* of the imaginary axis under $H$. All three are statements about a function from $\mathbb{C}$ to $\mathbb{C}$. To use them fluently we must understand such functions as objects in their own right — in particular as *mappings* that carry curves in one plane to curves in another.

> [!definition] Complex function
> Let $D\subseteq\mathbb{C}$. A **complex function** $f:D\to\mathbb{C}$ assigns to each $z\in D$ a complex number $w=f(z)$. Writing $z=x+jy$,
>
> $$
> f(z) = u(x,y) + j\,v(x,y),
> $$
>
> where $u=\operatorname{Re}\{f\}$ and $v=\operatorname{Im}\{f\}$ are real-valued functions of two real variables. In polar form $z=re^{j\theta}$ we may also write $f = u(r,\theta)+jv(r,\theta)$.

So a complex function is *a pair of real functions of two variables*. That viewpoint is honest but hides the structure; the payoff of complex analysis comes from treating $z$ as a single variable.

**Domains and regions.** We will need a few words about sets in the plane.

- An **open disc** $\lvert z-z_0\rvert<\rho$ is a *neighbourhood* of $z_0$; the **punctured disc** $0<\lvert z-z_0\rvert<\rho$ omits the centre.
- A set is **open** if every point has a neighbourhood inside the set. It is **connected** if any two points can be joined by a polygonal path inside it.
- A **domain** is an open connected set. A **region** is a domain possibly with some or all of its boundary points.
- A domain is **simply connected** if every closed curve in it can be shrunk continuously to a point without leaving it — informally, it has no holes. The disc $\lvert z\rvert<1$ is simply connected; the annulus $1<\lvert z\rvert<2$ and the punctured plane $z\neq0$ are not.

These distinctions are not pedantry: Cauchy's theorem (Section 18.5) holds on simply connected domains, and the Laurent series (Section 18.6) is the tool for annuli.

### Elementary functions

**Polynomials and rational functions.** $P(z)=a_0+a_1z+\dots+a_nz^n$ is defined everywhere. A rational function $R(z)=P(z)/Q(z)$ is defined except at zeros of $Q$ — the *poles* of every transfer function you have met.

**The exponential.** Defined by

$$
e^z = e^{x+jy} = e^x(\cos y + j\sin y),
$$

so $\lvert e^z\rvert = e^x$ and $\angle e^z = y$. It agrees with the real exponential on the real axis and with Euler's formula on the imaginary axis. Key properties:

- $e^{z_1+z_2}=e^{z_1}e^{z_2}$, and $e^z\neq0$ for every $z$;
- $e^z$ is **periodic with period $2\pi j$**: $e^{z+2\pi j}=e^z$. The horizontal strip $-\pi<y\le\pi$ already produces every nonzero value once.

**Trigonometric and hyperbolic functions.** Defined from the exponential,

$$
\cos z = \frac{e^{jz}+e^{-jz}}{2},\qquad \sin z=\frac{e^{jz}-e^{-jz}}{2j},\qquad \cosh z=\frac{e^z+e^{-z}}{2},\qquad \sinh z=\frac{e^z-e^{-z}}{2}.
$$

Hence $\cos(jy)=\cosh y$ and $\sin(jy)=j\sinh y$, and

$$
\cos(x+jy) = \cos x\cosh y - j\sin x\sinh y,\qquad \sin(x+jy)=\sin x\cosh y + j\cos x\sinh y .
$$

A striking consequence: $\lvert\sin z\rvert^2=\sin^2x+\sinh^2y$, which grows without bound as $\lvert y\rvert\to\infty$. **Complex sine and cosine are not bounded** — unlike their real versions. (Liouville's theorem in Section 18.5 will say no non-constant entire function can be bounded.)

**The logarithm — the first multivalued function.** We want $w=\log z$ to solve $e^w=z$. Write $z=re^{j\theta}$ and $w=u+jv$: then $e^ue^{jv}=re^{j\theta}$ gives $u=\ln r$ and $v=\theta+2\pi k$. So

$$
\log z = \ln\lvert z\rvert + j\arg z = \ln r + j(\theta+2\pi k),\qquad k\in\mathbb{Z}.
$$

Every nonzero $z$ has infinitely many logarithms, differing by integer multiples of $2\pi j$. To get a genuine (single-valued) function we must choose.

> [!definition] Principal logarithm and branch cut
> The **principal argument** $\operatorname{Arg}z\in(-\pi,\pi]$ and the **principal logarithm**
>
> $$
> \operatorname{Log}z = \ln\lvert z\rvert + j\operatorname{Arg}z,\qquad z\neq 0,
> $$
>
> define a single-valued function. It is continuous everywhere except on the negative real axis $\{x\le0,\ y=0\}$, where $\operatorname{Arg}$ jumps by $2\pi$. That excluded ray is a **branch cut**; $z=0$, around which the values wind, is a **branch point**. Any other choice of a ray from $0$ and a $2\pi$-wide interval of angles gives another **branch** of $\log$.

**Complex powers.** For complex $\alpha$, define $z^\alpha = e^{\alpha\log z}$; the **principal value** uses $\operatorname{Log}$. If $\alpha=n$ is an integer this is single-valued; if $\alpha=p/q$ is rational in lowest terms there are exactly $q$ distinct values (the $q$-th roots of Ch. 2); otherwise there are infinitely many. For example $z^{1/2}$ has two branches, $\pm\sqrt r\,e^{j\theta/2}$, which swap as $z$ circles the origin.

> [!example] Worked Example 18.1 — real and imaginary parts
> Find $u$ and $v$ for (a) $f(z)=z^2$, (b) $f(z)=1/z$, (c) $f(z)=e^{z}$, and evaluate each at $z=1+j$.
>
> **(a)** $(x+jy)^2=x^2-y^2+j\,2xy$, so $u=x^2-y^2$, $v=2xy$. At $1+j$: $u=0$, $v=2$, i.e. $f=2j$.
>
> **(b)** Multiply by the conjugate:
>
> $$
> \frac1z=\frac{\bar z}{\lvert z\rvert^2}=\frac{x-jy}{x^2+y^2},\qquad u=\frac{x}{x^2+y^2},\quad v=\frac{-y}{x^2+y^2}.
> $$
>
> At $1+j$: $u=\tfrac12$, $v=-\tfrac12$, i.e. $f=\tfrac12-\tfrac12j$. Check: $(1+j)(\tfrac12-\tfrac12j)=\tfrac12+\tfrac12=1$. ✓
>
> **(c)** $u=e^x\cos y$, $v=e^x\sin y$. At $1+j$: $e\cos1 = 1.4687$, $e\sin1=2.2874$.

> [!example] Worked Example 18.2 — logarithms and powers
> Find (a) $\operatorname{Log}(-1-j)$, (b) all values of $\log(-1)$, (c) the principal value of $j^j$.
>
> **(a)** $\lvert-1-j\rvert=\sqrt2$ and $-1-j$ lies in the third quadrant, so $\operatorname{Arg}(-1-j)=-3\pi/4$ (not $5\pi/4$, which lies outside $(-\pi,\pi]$). Hence $\operatorname{Log}(-1-j)=\tfrac12\ln2 - j\tfrac{3\pi}{4}=0.3466-2.3562j$.
>
> **(b)** $\lvert-1\rvert=1$, $\arg(-1)=\pi+2\pi k$, so $\log(-1)=j(2k+1)\pi$. The principal value is $j\pi$ — the logarithm of a negative real number is perfectly well defined in $\mathbb{C}$.
>
> **(c)** $j^j=e^{j\operatorname{Log}j}=e^{j\cdot j\pi/2}=e^{-\pi/2}=0.20788$. A complex number raised to a complex power can be real. The other values are $e^{-\pi/2-2\pi k}$.

> [!warning] Common Pitfall — "laws" that fail for multivalued functions
> $\operatorname{Log}(z_1z_2)=\operatorname{Log}z_1+\operatorname{Log}z_2$ is **false** in general. With $z_1=z_2=-1$: the left side is $\operatorname{Log}1=0$, the right side is $2\pi j$. The identity holds only modulo $2\pi j$. Likewise $(z_1z_2)^{1/2}\ne z_1^{1/2}z_2^{1/2}$ for principal branches in general, and $\sqrt{z^2}\ne z$. In software, `log`, `sqrt` and `angle` return principal values, so a phase computed as `angle(H)` jumps by $2\pi$ when the true phase crosses $\pm\pi$ — the reason phase plots must be *unwrapped*.

### Functions as mappings

A real function $y=f(x)$ is pictured by its graph in the $(x,y)$ plane. A complex function would need four real dimensions for its graph, so instead we draw **two planes**: the $z$-plane (coordinates $x,y$) and the $w$-plane (coordinates $u,v$), and we watch where $f$ sends curves. A good test pattern is the rectangular grid $x=\text{const}$, $y=\text{const}$, or the polar grid $r=\text{const}$, $\theta=\text{const}$.

**The map $w=z^2$.** In polar form $w=r^2e^{j2\theta}$: moduli are squared and angles doubled. The first quadrant $0\le\theta\le\pi/2$ opens out to the whole upper half-plane. On the rectangular grid:

- the vertical line $x=a$ ($a\neq0$) becomes $u=a^2-y^2$, $v=2ay$. Eliminating $y=v/(2a)$:
  $$
  u = a^2 - \frac{v^2}{4a^2},
  $$
  a parabola opening to the left with focus at the origin;
- the horizontal line $y=b$ ($b\neq0$) becomes $u=x^2-b^2$, $v=2bx$, i.e.
  $$
  u = \frac{v^2}{4b^2}-b^2,
  $$
  a parabola opening to the right, also with focus at the origin.

The two families cross at right angles — exactly as the original grid lines did. That is not a coincidence: it is **conformality**, developed in Section 18.3.

```desmos-graph
left=-5; right=5; bottom=-5; top=5
---
x=0.25-y^2|blue
x=1-\frac{y^2}{4}|blue
x=2.25-\frac{y^2}{9}|blue
x=4-\frac{y^2}{16}|blue
x=y^2-0.25|red
x=\frac{y^2}{4}-1|red
x=\frac{y^2}{9}-2.25|red
x=\frac{y^2}{16}-4|red
```
*Images under $w=z^2$ of the vertical lines $x=0.5,1,1.5,2$ (blue, left-opening parabolas) and horizontal lines $y=0.5,1,1.5,2$ (red, right-opening parabolas), drawn in the $w=u+jv$ plane; the two families of confocal parabolas meet at right angles.*

**The map $w=1/z$ (inversion).** In polar form $w=(1/r)e^{-j\theta}$: the modulus is inverted and the angle negated (reflection in the real axis). Points inside the unit circle go outside, and vice versa; the origin and the "point at infinity" swap. Using $u=x/(x^2+y^2)$, $v=-y/(x^2+y^2)$ and the inverse relations $x=u/(u^2+v^2)$, $y=-v/(u^2+v^2)$:

- the vertical line $x=a$ becomes $u/(u^2+v^2)=a$, i.e.
  $$
  \left(u-\frac{1}{2a}\right)^2+v^2=\frac{1}{4a^2},
  $$
  a circle through the origin with centre on the real axis;
- the horizontal line $y=b$ becomes
  $$
  u^2+\left(v+\frac{1}{2b}\right)^2=\frac{1}{4b^2},
  $$
  a circle through the origin with centre on the imaginary axis.

Lines map to circles through the origin (a line is a "circle through $\infty$", and $\infty\mapsto0$). In general, **inversion maps circles-and-lines to circles-and-lines**, a property shared by all **Möbius (bilinear) transformations** $w=(az+b)/(cz+d)$, $ad-bc\neq0$. The Smith chart is the Möbius map $\Gamma=(z_L-1)/(z_L+1)$ applied to the grid of normalized resistance and reactance; its circles are these circles.

```desmos-graph
left=-3; right=3; bottom=-3; top=3
---
(x-1)^2+y^2=1|blue
(x-0.5)^2+y^2=0.25|blue
(x+1)^2+y^2=1|blue
(x+0.5)^2+y^2=0.25|blue
x^2+(y-1)^2=1|red
x^2+(y-0.5)^2=0.25|red
x^2+(y+1)^2=1|red
x^2+(y+0.5)^2=0.25|red
```
*Images under $w=1/z$ of the lines $x=\pm0.5,\pm1$ (blue circles centred on the real axis) and $y=\pm0.5,\pm1$ (red circles centred on the imaginary axis); every circle passes through $w=0$, the image of $z=\infty$, and blue meets red at right angles.*

**The map $w=e^z$.** With $u+jv=e^xe^{jy}$: the vertical line $x=a$ maps to the circle $\lvert w\rvert=e^a$, traversed once each time $y$ advances by $2\pi$; the horizontal line $y=b$ maps to the ray $\angle w=b$. The exponential turns the Cartesian grid into the polar grid. Its inverse, $\operatorname{Log}$, does the reverse: it straightens circles and rays into vertical and horizontal lines. That is exactly why the $s$-plane and $z$-plane of sampled systems are related by $z=e^{sT}$ ([Ch. 17](17-z-transforms.md)): the imaginary axis $s=j\omega$ wraps onto the unit circle, the left half-plane maps into the unit disc, and the strip $-\pi/T<\omega\le\pi/T$ is the primary strip that maps onto the whole $z$-plane once — the complex-analysis statement of aliasing.

> [!example] Worked Example 18.3 — image of a region under $e^z$ and $z^2$
> (a) Find the image of the rectangle $0\le x\le1$, $0\le y\le\pi/2$ under $w=e^z$. (b) Find the image of the half-strip $0<x<1$, $y>0$ under $w=z^2$.
>
> **(a)** $\lvert w\rvert=e^x$ ranges over $[1,e]$ and $\angle w=y$ over $[0,\pi/2]$. The image is the quarter-annulus $1\le\lvert w\rvert\le e$, $0\le\angle w\le\pi/2$. Its area is $\tfrac14\pi(e^2-1)=5.0180$, compared with $\pi/2=1.5708$ for the rectangle. The local area magnification is $\lvert f'(z)\rvert^2=e^{2x}$, and indeed $\int_0^1\!\int_0^{\pi/2}e^{2x}\,dy\,dx=\tfrac{\pi}{2}\cdot\tfrac{e^2-1}{2}=5.0180$. ✓
>
> **(b)** The boundary $x=0$, $y>0$ maps to $u=-y^2\le0$, $v=0$: the negative real axis. The boundary $x=1$ maps to the parabola $u=1-v^2/4$ with $v=2y>0$, the upper half of the left-opening parabola. The boundary segment $y=0$, $0<x<1$ maps to $0<u<1$, $v=0$. The image is the region in the upper half-plane lying inside (to the left of) the parabola $u=1-v^2/4$ — test the interior point $z=0.5+0.5j$: $w=0.5j$, and $0<1-0.25/4$. ✓

> [!tip] Engineering Note — what a Nyquist plot really is
> A Nyquist plot draws the curve $\{L(j\omega):-\infty<\omega<\infty\}$ in the $w$-plane: the **image of the imaginary axis** under the loop transfer function $L(s)$. Similarly, a polar plot of an impedance $Z(j\omega)$ is the image of the frequency axis. Everything in this section — lines becoming circles, grids staying orthogonal — is what you are seeing in those plots. Section 18.9 uses the argument principle to read stability from such an image.

## 18.2 Analyticity

### Why differentiability is special in the plane

In one real variable, $x$ can approach $x_0$ only from the left or the right. In the complex plane, $z$ can approach $z_0$ along infinitely many directions and along spirals. Demanding that a difference quotient have the *same* limit along every one of these paths is a far stronger condition than real differentiability, and the functions that meet it — analytic functions — inherit remarkable properties. The transfer functions, impedances and transforms of circuit theory are analytic away from their poles; this is why their behaviour on the $j\omega$-axis is tied so tightly to their pole–zero pattern.

### Limits and continuity

> [!definition] Limit and continuity
> $\displaystyle\lim_{z\to z_0}f(z)=L$ means: for every $\varepsilon>0$ there is $\delta>0$ such that $\lvert f(z)-L\rvert<\varepsilon$ whenever $0<\lvert z-z_0\rvert<\delta$. The function is **continuous** at $z_0$ if $f(z_0)$ is defined and $\lim_{z\to z_0}f(z)=f(z_0)$.

Because $\lvert z-z_0\rvert<\delta$ is a disc, the limit must be the same however $z$ approaches $z_0$. Writing $f=u+jv$ and $L=a+jb$, the complex limit exists if and only if $u\to a$ and $v\to b$ as two-variable limits. Sums, products, quotients (with nonzero denominator) and compositions of continuous functions are continuous; polynomials are continuous everywhere, rational functions wherever the denominator is nonzero.

**Limits involving infinity.** We say $f(z)\to\infty$ as $z\to z_0$ if $\lvert f(z)\rvert\to\infty$, regardless of direction. Adding a single point $\infty$ to $\mathbb{C}$ gives the **extended complex plane** (the Riemann sphere). Behaviour "at infinity" is studied by substituting $z=1/\zeta$ and examining $\zeta\to0$: for instance $H(s)=1/(s+1)$ has $H(1/\zeta)=\zeta/(1+\zeta)\to0$, a zero at infinity — the "strictly proper" property of transfer functions.

> [!example] Worked Example 18.4 — a limit that does not exist
> Does $\lim_{z\to0}\dfrac{\bar z}{z}$ exist?
>
> Along the real axis, $z=x$: $\bar z/z=x/x=1$. Along the imaginary axis, $z=jy$: $\bar z/z=-jy/(jy)=-1$. Along the line $y=x$, $z=x(1+j)$: $\bar z/z=(1-j)/(1+j)=-j$. In polar form $\bar z/z=e^{-2j\theta}$, which takes every value on the unit circle depending on the direction of approach. The limit does not exist.

### The complex derivative

> [!definition] Complex derivative
> $f$ is **differentiable** at $z_0$ if
>
> $$
> f'(z_0) = \lim_{\Delta z\to0}\frac{f(z_0+\Delta z)-f(z_0)}{\Delta z}
> $$
>
> exists (as a complex limit, independent of how $\Delta z\to0$).

The algebra is identical to the real case, so all the familiar rules carry over: linearity, the product rule, the quotient rule, the chain rule, and $\dfrac{d}{dz}z^n=nz^{n-1}$ for integer $n$ (with $z\ne0$ when $n<0$). For example, from $(z+\Delta z)^2-z^2=2z\Delta z+(\Delta z)^2$ we get $(z^2)'=\lim(2z+\Delta z)=2z$, independent of direction. Likewise $(e^z)'=e^z$, $(\sin z)'=\cos z$, $(\cos z)'=-\sin z$, and on its cut plane $(\operatorname{Log}z)'=1/z$.

> [!theorem] Differentiable implies continuous
> If $f'(z_0)$ exists then $f$ is continuous at $z_0$, since $f(z_0+\Delta z)-f(z_0)=\dfrac{f(z_0+\Delta z)-f(z_0)}{\Delta z}\cdot\Delta z\to f'(z_0)\cdot0=0$.

The converse fails dramatically: there are continuous functions, like $\bar z$, that are differentiable nowhere.

> [!example] Worked Example 18.5 — $\bar z$ and $\lvert z\rvert^2$
> (a) Show $f(z)=\bar z$ is nowhere differentiable. (b) Show $g(z)=\lvert z\rvert^2$ is differentiable only at $z=0$.
>
> **(a)** $\dfrac{\overline{z+\Delta z}-\bar z}{\Delta z}=\dfrac{\overline{\Delta z}}{\Delta z}$, which by Example 18.4 has no limit as $\Delta z\to0$. So $\bar z$ — a perfectly smooth map of the plane (a reflection) — has no complex derivative anywhere.
>
> **(b)** $\lvert z+\Delta z\rvert^2-\lvert z\rvert^2=(z+\Delta z)(\bar z+\overline{\Delta z})-z\bar z=\bar z\Delta z+z\overline{\Delta z}+\Delta z\overline{\Delta z}$. Dividing,
>
> $$
> \frac{\Delta g}{\Delta z}=\bar z + z\frac{\overline{\Delta z}}{\Delta z}+\overline{\Delta z}.
> $$
>
> The last term $\to0$. The middle term has no limit unless $z=0$. So $g'(0)=0$ and $g'(z)$ does not exist for $z\ne0$. Note that $g=x^2+y^2$ is infinitely differentiable as a real function of $(x,y)$ — real smoothness is not enough.

The lesson: a function that genuinely depends on $\bar z$ is not complex-differentiable. Analytic functions are, informally, functions of $z$ alone. (One can make this precise with the Wirtinger derivative $\partial f/\partial\bar z=\tfrac12(f_x+jf_y)$; the Cauchy–Riemann equations of Section 18.3 say exactly that it vanishes.)

### Analytic functions

Differentiability at a single point, as in $\lvert z\rvert^2$ at $0$, is a curiosity. What matters is differentiability throughout a neighbourhood.

> [!definition] Analytic (holomorphic) function
> $f$ is **analytic at $z_0$** if it is differentiable at every point of some neighbourhood $\lvert z-z_0\rvert<\rho$. It is **analytic in a domain $D$** if it is analytic at every point of $D$. A function analytic in the whole finite plane is **entire**. A point where $f$ fails to be analytic, but every neighbourhood of which contains points of analyticity, is a **singular point** (singularity) of $f$.

Examples:

- Polynomials, $e^z$, $\sin z$, $\cos z$, $\sinh z$, $\cosh z$ are entire.
- A rational function $P/Q$ (in lowest terms) is analytic except at the zeros of $Q$, which are its singular points (poles).
- $\operatorname{Log}z$ is analytic in the plane cut along the negative real axis, $\mathbb{C}\setminus(-\infty,0]$. Points on the cut are singular for this branch, though another branch is analytic there; $z=0$ is singular for every branch.
- $\tan z=\sin z/\cos z$ is analytic except where $\cos z=0$; since $\lvert\cos z\rvert^2=\cos^2x+\sinh^2y$, these are exactly $z=\pi/2+k\pi$ on the real axis — the complex extension creates no new singularities.
- $\bar z$, $\operatorname{Re}z$, $\lvert z\rvert$ are analytic nowhere.

> [!theorem] Algebra of analytic functions
> Sums, products and compositions of analytic functions are analytic; a quotient $f/g$ is analytic wherever $g\neq0$. If $f$ is analytic and $f'(z)=0$ throughout a domain, $f$ is constant there.

> [!example] Worked Example 18.6 — locating singular points
> Where is $f(z)=\dfrac{z^2+1}{z^3+4z}$ analytic, and what is $f'(z)$?
>
> The denominator is $z(z^2+4)=z(z-2j)(z+2j)$. The numerator $(z-j)(z+j)$ shares no factor with it, so $f$ is analytic on $\mathbb{C}\setminus\{0,2j,-2j\}$. By the quotient rule,
>
> $$
> f'(z)=\frac{2z(z^3+4z)-(z^2+1)(3z^2+4)}{(z^3+4z)^2}=\frac{-z^4+z^2-4}{(z^3+4z)^2}.
> $$
>
> Check at $z=1$: the formula gives $-4/25=-0.16$; a central difference quotient $[f(1.0001)-f(0.9999)]/0.0002$ also gives $-0.16000$. Because $f$ is analytic at $1$, the difference quotient taken in the imaginary direction, $[f(1+0.0001j)-f(1-0.0001j)]/(0.0002j)$, gives the same value $-0.16000$. ✓

> [!example] Worked Example 18.7 — a real function explained by complex singularities
> The real function $g(x)=1/(1+x^2)$ is infinitely differentiable on all of $\mathbb{R}$, yet its Maclaurin series $1-x^2+x^4-\dots$ converges only for $\lvert x\rvert<1$. Why?
>
> Nothing on the real axis explains the radius $1$. Extend to $g(z)=1/(1+z^2)=1/[(z-j)(z+j)]$: it has singular points at $z=\pm j$, at distance $1$ from the origin. We will prove in Section 18.6 that a Taylor series about $z_0$ converges in the largest disc about $z_0$ that avoids all singularities — here $\lvert z\rvert<1$, whose intersection with the real line is $(-1,1)$. The real series "feels" poles it cannot see. Expanding instead about $x_0=1$, the nearest singularities are at distance $\lvert1-j\rvert=\sqrt2$, so that series converges for $\lvert x-1\rvert<\sqrt2=1.414$.

### What analyticity buys

The results proved later in the chapter, previewed here, show how much the single hypothesis "analytic in $D$" delivers:

1. **Harmonic parts** (Section 18.3): $u$ and $v$ satisfy Laplace's equation $\nabla^2u=\nabla^2v=0$.
2. **Conformality** (Section 18.3): wherever $f'\ne0$, $f$ preserves angles between curves.
3. **Path independence** (Section 18.5): $\oint_C f\,dz=0$ around any closed curve in a simply connected $D$.
4. **Infinite differentiability** (Section 18.5): $f'$ is analytic too, hence $f''$, $f'''$, … all exist.
5. **Power series** (Section 18.6): $f$ equals its Taylor series in the largest disc about $z_0$ contained in $D$.
6. **Identity theorem**: if two functions analytic in $D$ agree on a segment of curve (or any set with a limit point in $D$), they agree on all of $D$. This underlies *analytic continuation*: the Laplace transform $1/(s+a)$, computed as an integral only for $\operatorname{Re}s>-a$, has a unique analytic extension to $\mathbb{C}\setminus\{-a\}$, and that is the function whose "poles" we speak of.

> [!tip] Engineering Note — causality, stability and analyticity
> The Laplace transform $H(s)=\int_0^\infty h(t)e^{-st}dt$ of a causal impulse response is analytic throughout its region of convergence, a right half-plane $\operatorname{Re}s>\sigma_0$ ([Ch. 14](14-laplace-transforms.md)). A stable causal system has $\sigma_0<0$, so $H$ is analytic on and to the right of the $j\omega$-axis. Analyticity in the right half-plane is therefore the mathematical signature of a *causal, stable* system, and it is what forces relationships between gain and phase (Bode's gain–phase relation, the Kramers–Kronig relations): you cannot specify $\lvert H(j\omega)\rvert$ and $\angle H(j\omega)$ independently.

> [!warning] Common Pitfall — "analytic at a point"
> A function can be differentiable at a point without being analytic there ($\lvert z\rvert^2$ at $0$). Analyticity is a neighbourhood property. When a theorem says "analytic inside and on $C$", it means analytic on some open set containing the curve and its interior.

## 18.3 Cauchy–Riemann equations

### Motivation

The definition of the derivative is hard to check directly. We want a practical test in terms of $u$ and $v$ — something we can compute with ordinary partial derivatives. The test, the Cauchy–Riemann equations, does more than decide differentiability: it reveals that the real and imaginary parts of an analytic function are *harmonic* (they solve Laplace's equation, the equation of electrostatic potential, steady heat flow and ideal fluid flow) and that analytic maps preserve angles. Those two facts turn analytic functions into a toolkit for two-dimensional field problems.

### Derivation

Suppose $f=u+jv$ is differentiable at $z_0=x_0+jy_0$. The limit defining $f'(z_0)$ must be the same along every path, in particular along the two coordinate directions.

*Horizontal approach*, $\Delta z=\Delta x$ real:

$$
f'(z_0)=\lim_{\Delta x\to0}\frac{u(x_0+\Delta x,y_0)-u(x_0,y_0)}{\Delta x}+j\frac{v(x_0+\Delta x,y_0)-v(x_0,y_0)}{\Delta x}=u_x+jv_x .
$$

*Vertical approach*, $\Delta z=j\Delta y$:

$$
f'(z_0)=\lim_{\Delta y\to0}\frac{[u(x_0,y_0+\Delta y)-u]+j[v(x_0,y_0+\Delta y)-v]}{j\Delta y}=\frac{1}{j}(u_y+jv_y)=v_y-ju_y .
$$

Equating real and imaginary parts of the two expressions gives the equations.

> [!theorem] Cauchy–Riemann equations (necessary condition)
> If $f=u+jv$ is differentiable at $z_0$, then the first partial derivatives of $u$ and $v$ exist at $(x_0,y_0)$ and satisfy
>
> $$
> u_x = v_y,\qquad u_y = -v_x .
> $$
>
> Moreover $f'(z_0)=u_x+jv_x=v_y-ju_y$.

The equations are necessary but, on their own, not quite sufficient (there are pathological functions satisfying them at a point without being differentiable). A mild smoothness hypothesis closes the gap.

> [!theorem] Sufficient condition
> Let $u$ and $v$ have continuous first partial derivatives in a neighbourhood of $z_0$, and let them satisfy the Cauchy–Riemann equations at $z_0$. Then $f'(z_0)$ exists. If this holds at every point of a domain $D$, $f$ is analytic in $D$.

*Proof sketch.* Continuity of the partials makes $u$ and $v$ differentiable as real functions, so $\Delta u=u_x\Delta x+u_y\Delta y+\varepsilon_1\lvert\Delta z\rvert$ and $\Delta v=v_x\Delta x+v_y\Delta y+\varepsilon_2\lvert\Delta z\rvert$ with $\varepsilon_{1,2}\to0$. Substituting $u_y=-v_x$, $v_y=u_x$:

$$
\Delta f=\Delta u+j\Delta v=(u_x+jv_x)(\Delta x+j\Delta y)+(\varepsilon_1+j\varepsilon_2)\lvert\Delta z\rvert ,
$$

so $\Delta f/\Delta z\to u_x+jv_x$ because $\lvert\Delta z\rvert/\Delta z$ has modulus $1$. $\square$

**Matrix viewpoint.** As a map of $\mathbb{R}^2$, $f$ has Jacobian $\begin{pmatrix}u_x&u_y\\v_x&v_y\end{pmatrix}$. Under Cauchy–Riemann it becomes $\begin{pmatrix}a&-b\\b&a\end{pmatrix}$ with $a=u_x$, $b=v_x$ — a rotation by $\angle f'$ times a scaling by $\lvert f'\rvert$. That is precisely the matrix of multiplication by the complex number $f'=a+jb$. A complex-differentiable function is one whose local linear approximation is *multiplication by a complex number*, not an arbitrary $2\times2$ matrix. The Jacobian determinant is $a^2+b^2=\lvert f'\rvert^2$, the local area magnification used in Example 18.3.

**Polar form.** With $z=re^{j\theta}$ and $u(r,\theta)$, $v(r,\theta)$, the chain rule converts the equations to

$$
u_r=\frac{1}{r}v_\theta,\qquad v_r=-\frac{1}{r}u_\theta,\qquad f'(z)=e^{-j\theta}(u_r+jv_r),\qquad r>0 .
$$

These are convenient for functions built from $r$ and $\theta$, such as powers and logarithms.

> [!example] Worked Example 18.8 — applying the test
> Test (a) $f=z^2$, (b) $f=e^z$, (c) $f=\bar z$, (d) $f=x^2+jy^2$, (e) $f=\operatorname{Log}z$.
>
> **(a)** $u=x^2-y^2$, $v=2xy$: $u_x=2x=v_y$, $u_y=-2y=-v_x$. ✓ everywhere; $f'=u_x+jv_x=2x+2jy=2z$.
>
> **(b)** $u=e^x\cos y$, $v=e^x\sin y$: $u_x=e^x\cos y=v_y$; $u_y=-e^x\sin y=-v_x$. ✓ everywhere, $f'=e^x\cos y+je^x\sin y=e^z$. The exponential is its own derivative in $\mathbb{C}$ too.
>
> **(c)** $u=x$, $v=-y$: $u_x=1$ but $v_y=-1$. ✗ everywhere — confirming Example 18.5.
>
> **(d)** $u=x^2$, $v=y^2$: $u_x=2x$, $v_y=2y$, so the first equation needs $x=y$; the second, $0=0$, always holds. The partials are continuous, so $f$ is differentiable exactly on the line $y=x$, with $f'=2x$ there. A line has no interior, so $f$ is **analytic nowhere**.
>
> **(e)** $u=\ln r$, $v=\theta$ (principal branch, $-\pi<\theta<\pi$): $u_r=1/r=\tfrac1r v_\theta$ and $v_r=0=-\tfrac1r u_\theta$. ✓ on the cut plane, $f'=e^{-j\theta}\cdot\tfrac1r=1/z$.

### Harmonic functions

> [!theorem] Real and imaginary parts are harmonic
> If $f=u+jv$ is analytic in $D$ (and $u,v$ have continuous second partials — automatically true, as Section 18.5 shows), then
>
> $$
> \nabla^2u=u_{xx}+u_{yy}=0,\qquad \nabla^2v=v_{xx}+v_{yy}=0\quad\text{in }D .
> $$

*Proof.* Differentiate $u_x=v_y$ with respect to $x$ and $u_y=-v_x$ with respect to $y$: $u_{xx}=v_{yx}$ and $u_{yy}=-v_{xy}$. Adding, and using equality of mixed partials, $u_{xx}+u_{yy}=0$. The same argument works for $v$. $\square$

A function satisfying Laplace's equation is **harmonic**. If $u$ and $v$ are harmonic and satisfy Cauchy–Riemann together, $v$ is a **harmonic conjugate** of $u$. Every harmonic $u$ on a simply connected domain has a harmonic conjugate, unique up to an additive constant, found by integrating the Cauchy–Riemann equations.

**Orthogonality of level curves.** The gradients are $\nabla u=(u_x,u_y)$ and $\nabla v=(v_x,v_y)=(-u_y,u_x)$, so

$$
\nabla u\cdot\nabla v=-u_xu_y+u_yu_x=0 .
$$

Where $f'\neq0$, the curves $u=\text{const}$ and $v=\text{const}$ cross at right angles. In electrostatics, if $u$ is the potential then $v=\text{const}$ are the field lines (flux lines); in fluid flow they are the equipotentials and streamlines.

```desmos-graph
left=-3; right=3; bottom=-3; top=3
---
x^2-y^2=0.5|blue
x^2-y^2=2|blue
x^2-y^2=-0.5|blue
x^2-y^2=-2|blue
2xy=0.5|red
2xy=2|red
2xy=-0.5|red
2xy=-2|red
```
*Level curves of $u=\operatorname{Re}\{z^2\}=x^2-y^2$ (blue) and $v=\operatorname{Im}\{z^2\}=2xy$ (red): two families of hyperbolas that cross at right angles everywhere except at the critical point $z=0$, where $f'(0)=0$. They are the equipotentials and field lines near a right-angled conducting corner.*

> [!example] Worked Example 18.9 — constructing harmonic conjugates
> (a) Show $u=x^2-y^2+2x+1$ is harmonic, find a harmonic conjugate $v$ and the analytic $f$. (b) Repeat for $u=e^{-y}\cos x$. (c) Is $u=x^2+y^2$ the real part of an analytic function?
>
> **(a)** $u_{xx}+u_{yy}=2-2=0$. ✓ From $v_y=u_x=2x+2$, integrate in $y$: $v=2xy+2y+h(x)$. Then $v_x=2y+h'(x)$ must equal $-u_y=2y$, so $h'=0$ and $h=C$. Thus $v=2xy+2y+C$ and
>
> $$
> f=u+jv=(x^2-y^2+2x+1)+j(2xy+2y)+jC=(z+1)^2+jC .
> $$
>
> Check at $z=j$: $(1+j)^2=2j$, while $u(0,1)=-1+1=0$ and $v(0,1)=2$. ✓
>
> **(b)** $u_{xx}=-e^{-y}\cos x$, $u_{yy}=e^{-y}\cos x$; sum $0$. ✓ From $v_y=u_x=-e^{-y}\sin x$: $v=e^{-y}\sin x+h(x)$. Then $v_x=e^{-y}\cos x+h'$ must equal $-u_y=e^{-y}\cos x$, so $h=C$. Hence $f=e^{-y}(\cos x+j\sin x)=e^{-y+jx}=e^{jz}$ (plus $jC$).
>
> **(c)** $\nabla^2u=2+2=4\neq0$. Not harmonic, so it cannot be the real part of any analytic function.

**A shortcut (Milne-Thomson).** Since $f'=u_x-ju_y$, set $y=0$ and $x=z$: $f'(z)=u_x(z,0)-ju_y(z,0)$, then integrate. In (a): $u_x(z,0)=2z+2$, $u_y(z,0)=0$, so $f'=2z+2$ and $f=z^2+2z+c$ with $c=u(0,0)=1$ (up to $jC$). The method relies on the identity theorem: an analytic function is fixed by its values on the real axis.

### Conformal mapping

Here is the geometric heart of the Cauchy–Riemann equations. Let a smooth curve $z(t)$ pass through $z_0=z(t_0)$ with tangent vector $z'(t_0)\neq0$. Its image $w(t)=f(z(t))$ has tangent, by the chain rule,

$$
w'(t_0)=f'(z_0)\,z'(t_0).
$$

Multiplying by the fixed complex number $f'(z_0)=\lvert f'(z_0)\rvert e^{j\alpha}$ rotates *every* tangent direction at $z_0$ by the same angle $\alpha=\angle f'(z_0)$ and stretches it by the same factor $\lvert f'(z_0)\rvert$. Two curves crossing at $z_0$ at angle $\beta$ therefore have images crossing at angle $\beta$ too.

> [!theorem] Conformality of analytic maps
> If $f$ is analytic at $z_0$ and $f'(z_0)\neq0$, then $f$ is **conformal** at $z_0$: it preserves both the magnitude and the sense of the angle between any two smooth curves through $z_0$. Infinitesimally small figures near $z_0$ are rotated by $\angle f'(z_0)$ and scaled by $\lvert f'(z_0)\rvert$, so small squares map to (nearly) small squares, with area scaled by $\lvert f'(z_0)\rvert^2$.
>
> At a **critical point**, where $f'(z_0)=0$ and the first nonvanishing derivative is $f^{(m)}(z_0)$, angles are multiplied by $m$.

This is what the two graphs of Section 18.1 showed: the orthogonal grid stayed orthogonal under $z^2$ and $1/z$. At $z=0$, $(z^2)'=0$ with $m=2$: the right angle between the positive real and imaginary axes opens into the straight angle between the positive and negative real $w$-axes.

**Why conformal maps solve field problems.** If $\Phi(u,v)$ is harmonic in the $w$-plane and $w=f(z)$ is analytic, then $\phi(x,y)=\Phi(u(x,y),v(x,y))$ is harmonic in the $z$-plane. (Locally $\Phi=\operatorname{Re}F(w)$ for analytic $F$, and $F(f(z))$ is analytic.) Conductor boundaries map to boundaries, and the boundary values carry over. So a hard geometry in the $z$-plane can be mapped to an easy one — parallel plates, a half-plane — solved there, and mapped back. The *complex potential* $\Omega=\phi+j\psi$ packages potential $\phi$ and flux function $\psi$ into one analytic function.

> [!example] Worked Example 18.10 — coaxial-line capacitance by conformal mapping
> A coaxial line has inner conductor radius $a$ at potential $V_0$ and outer radius $b$ at $0$, filled with dielectric $\varepsilon$. Find the capacitance per unit length by mapping to a parallel-plate geometry.
>
> Use $w=\operatorname{Log}z=\ln r+j\theta$. The annulus $a<r<b$ maps to the rectangle $\ln a<u<\ln b$, $-\pi<v\le\pi$: the circle $r=a$ becomes the line $u=\ln a$ and $r=b$ the line $u=\ln b$. In the $w$-plane we have a parallel-plate capacitor of gap $d=\ln(b/a)$ and width $2\pi$ (per unit depth), whose edges $v=\pm\pi$ are glued together — so there are no fringing fields. The potential is linear in $u$:
>
> $$
> \Phi=V_0\frac{\ln b-u}{\ln(b/a)}\quad\Longrightarrow\quad \phi(r)=V_0\frac{\ln(b/r)}{\ln(b/a)} .
> $$
>
> A conformal map preserves capacitance (charge and potential are both preserved, since flux through corresponding curves is the same), so $C'=\varepsilon\cdot\dfrac{\text{width}}{\text{gap}}=\dfrac{2\pi\varepsilon}{\ln(b/a)}$.
>
> Numerically, with $a=0.5$ mm, $b=1.75$ mm and $\varepsilon_r=2.25$ (polyethylene): $C'=2\pi(8.854\times10^{-12})(2.25)/\ln3.5=99.9$ pF/m — the familiar value for a 50 $\Omega$ cable. Dimensions: $\varepsilon$ in F/m divided by a pure number gives F/m. ✓

> [!tip] Engineering Note — where conformal maps appear
> Microstrip and coplanar-waveguide impedance formulas come from Schwarz–Christoffel conformal maps (which send polygons to half-planes); the Smith chart is the Möbius map $\Gamma=(z_L-1)/(z_L+1)$; the bilinear transform $s=\tfrac{2}{T}\tfrac{z-1}{z+1}$ used to discretize analog filters is a conformal Möbius map of the left half $s$-plane onto the unit $z$-disc, which is why it preserves stability and why it warps frequency ($\omega_a=\tfrac2T\tan\tfrac{\omega_dT}{2}$).

> [!warning] Common Pitfall — checking only one equation, or only at one point
> Both Cauchy–Riemann equations must hold, and they must hold throughout an open set for analyticity. $f=x^2+jy^2$ satisfies them on a line and is analytic nowhere. Conversely, a function whose partials are discontinuous can satisfy the equations at a point and still fail to be differentiable there; quote the sufficiency theorem with its continuity hypothesis.

## 18.4 Contour integration

### Motivation

The inverse Laplace transform is an integral along a vertical line in the $s$-plane; the inverse Z-transform is an integral around a circle in the $z$-plane; the Fourier integral is an integral along the $j\omega$-axis. These are integrals of complex functions along curves — *contour integrals*. This section defines them, shows how to compute them from a parametrization, and establishes the two workhorses: the antiderivative rule and the ML bound. The next sections show that for analytic integrands the curve can be deformed almost at will, which is where the power lies.

### Curves and contours

> [!definition] Arcs, contours and orientation
> A **smooth arc** is a curve $z(t)=x(t)+jy(t)$, $a\le t\le b$, with $z'(t)$ continuous and nonzero. A **contour** (piecewise smooth path) is a finite chain of smooth arcs joined end to end. A contour is **closed** if $z(b)=z(a)$, and **simple** if it does not cross itself (apart from the endpoints of a closed one). A simple closed contour is **positively oriented** if it is traversed counter-clockwise — the interior lies to the left of the direction of travel.

Standard parametrizations:

| Curve | Parametrization |
|---|---|
| segment from $z_1$ to $z_2$ | $z(t)=z_1+(z_2-z_1)t,\ 0\le t\le1$ |
| circle $\lvert z-z_0\rvert=\rho$, counter-clockwise | $z(\theta)=z_0+\rho e^{j\theta},\ 0\le\theta\le2\pi$ |
| upper semicircle of radius $R$ | $z(\theta)=Re^{j\theta},\ 0\le\theta\le\pi$ |
| vertical line $\operatorname{Re}z=\gamma$, upward | $z(\omega)=\gamma+j\omega,\ -\infty<\omega<\infty$ |

The last is the Bromwich contour of the inverse Laplace transform (Section 18.9). On a circle, $dz=j\rho e^{j\theta}d\theta=j(z-z_0)\,d\theta$, an identity used constantly.

### Definition and evaluation

> [!definition] Contour integral
> Let $C$ be a contour $z(t)$, $a\le t\le b$, and $f$ continuous on $C$. Then
>
> $$
> \int_C f(z)\,dz=\int_a^b f\big(z(t)\big)\,z'(t)\,dt .
> $$
>
> It is the limit of Riemann sums $\sum f(\zeta_k)\,\Delta z_k$ with $\Delta z_k=z_k-z_{k-1}$, and its value does not depend on which parametrization of $C$ (with the same direction) is used. For a closed contour we write $\oint_C$.

Writing $f=u+jv$ and $dz=dx+j\,dy$ shows how the contour integral relates to the real line integrals of [Ch. 7](../part-2-calculus/07-vector-calculus.md):

$$
\int_Cf\,dz=\int_C(u\,dx-v\,dy)+j\int_C(v\,dx+u\,dy).
$$

This decomposition is what connects Cauchy's theorem to Green's theorem in the next section.

> [!theorem] Properties
> For contours $C$, $C_1$, $C_2$ and constants $\alpha,\beta$:
> 1. Linearity: $\int_C(\alpha f+\beta g)\,dz=\alpha\int_Cf\,dz+\beta\int_Cg\,dz$.
> 2. Reversal: $\int_{-C}f\,dz=-\int_Cf\,dz$, where $-C$ is $C$ traversed backwards.
> 3. Additivity: $\int_{C_1+C_2}f\,dz=\int_{C_1}f\,dz+\int_{C_2}f\,dz$ for $C_2$ starting where $C_1$ ends.
> 4. Modulus bound: $\left\lvert\int_Cf\,dz\right\rvert\le\int_C\lvert f(z)\rvert\,\lvert dz\rvert$, with $\lvert dz\rvert=\lvert z'(t)\rvert dt$ the element of arc length.

### The fundamental integral

> [!example] Worked Example 18.11 — integer powers around a circle
> Let $C$ be the circle $\lvert z-z_0\rvert=\rho$, counter-clockwise. Evaluate $\oint_C(z-z_0)^n\,dz$ for every integer $n$.
>
> Parametrize $z=z_0+\rho e^{j\theta}$, $dz=j\rho e^{j\theta}d\theta$:
>
> $$
> \oint_C(z-z_0)^n\,dz=\int_0^{2\pi}\rho^ne^{jn\theta}\,j\rho e^{j\theta}\,d\theta=j\rho^{n+1}\int_0^{2\pi}e^{j(n+1)\theta}\,d\theta .
> $$
>
> If $n\neq-1$, the integral of $e^{j(n+1)\theta}$ over a full period is $0$. If $n=-1$, the integrand is $1$ and the result is $2\pi j$:
>
> $$
> \oint_C(z-z_0)^n\,dz=\begin{cases}2\pi j,&n=-1\\0,&n\neq-1.\end{cases}
> $$
>
> The answer is independent of the radius $\rho$. This single computation is the seed of the whole residue theory: when a function is expanded in powers of $(z-z_0)$, only the $(z-z_0)^{-1}$ term survives integration around $z_0$.

### Path dependence and the antiderivative rule

> [!example] Worked Example 18.12 — $z^2$ versus $\bar z$ along two paths
> Integrate from $0$ to $1+j$ along (i) the straight segment $C_1$, and (ii) the path $C_2$ that runs along the real axis to $1$ and then straight up to $1+j$. Do this for $f=z^2$ and $f=\bar z$.
>
> **$f=z^2$, path (i):** $z=(1+j)t$, $dz=(1+j)dt$:
>
> $$
> \int_{C_1}z^2dz=\int_0^1(1+j)^3t^2\,dt=\frac{(1+j)^3}{3}=\frac{-2+2j}{3}.
> $$
>
> **$f=z^2$, path (ii):** on $[0,1]$, $\int_0^1x^2dx=\tfrac13$. On the vertical leg $z=1+jy$, $dz=j\,dy$: $\int_0^1(1+jy)^2j\,dy=j\int_0^1(1-y^2+2jy)\,dy=j\left(\tfrac23+j\right)=-1+\tfrac23j$. Total $\tfrac13-1+\tfrac23j=\tfrac{-2+2j}{3}$ — the same.
>
> **$f=\bar z$, path (i):** $\bar z=(1-j)t$: $\int_0^1(1-j)t(1+j)\,dt=2\cdot\tfrac12=1$.
>
> **$f=\bar z$, path (ii):** $\int_0^1x\,dx=\tfrac12$, then $\int_0^1(1-jy)j\,dy=\int_0^1(j+y)\,dy=j+\tfrac12$. Total $1+j$ — **different**.
>
> The analytic integrand gives a path-independent answer; the non-analytic one does not. The difference, $(1+j)-1=j$, equals $\oint\bar z\,dz$ around the closed triangle $C_2-C_1$, and in general $\oint_C\bar z\,dz=2j\times(\text{area enclosed})$ — here the triangle has area $\tfrac12$. ✓ (Around the unit circle this gives $2\pi j$, matching the direct computation $\int_0^{2\pi}e^{-j\theta}je^{j\theta}d\theta=2\pi j$.)

> [!theorem] Antiderivative (fundamental theorem for contour integrals)
> If $f$ is continuous in a domain $D$ and has an antiderivative $F$ there ($F'=f$ throughout $D$), then for every contour $C$ in $D$ from $z_1$ to $z_2$,
>
> $$
> \int_Cf(z)\,dz=F(z_2)-F(z_1),
> $$
>
> independent of the path. In particular $\oint_Cf\,dz=0$ for every closed contour in $D$.

*Proof.* $\frac{d}{dt}F(z(t))=F'(z(t))z'(t)=f(z(t))z'(t)$, so the integral is $[F(z(t))]_a^b$. $\square$

So $\int_{C}z^2dz=[z^3/3]_0^{1+j}$ immediately. Conversely, Example 18.11 shows $1/(z-z_0)$ has **no** antiderivative on any domain containing a circle around $z_0$: its "antiderivative" $\log(z-z_0)$ is multivalued, and going once around changes it by $2\pi j$.

> [!example] Worked Example 18.13 — $1/z$ and the branch of the logarithm
> Compute $\int dz/z$ from $-j$ to $j$ along (i) the right half of the unit circle and (ii) the left half.
>
> **(i)** The right half lies in the cut plane $\mathbb{C}\setminus(-\infty,0]$, where $\operatorname{Log}$ is an antiderivative: $\operatorname{Log}j-\operatorname{Log}(-j)=j\tfrac\pi2-(-j\tfrac\pi2)=j\pi$.
>
> **(ii)** The left half crosses the negative real axis, so use a branch cut along the positive real axis instead, with $\arg z\in(0,2\pi)$: then $\arg j=\pi/2$ and $\arg(-j)=3\pi/2$, giving $j\pi/2-j3\pi/2=-j\pi$.
>
> Directly: path (ii) is $z=e^{j\theta}$ with $\theta$ from $3\pi/2$ down to $\pi/2$, so $\int j\,d\theta=j(\pi/2-3\pi/2)=-j\pi$. ✓ The two paths differ by $2\pi j$: going up the right side and back down the left encircles the origin once counter-clockwise, $\oint dz/z=j\pi-(-j\pi)=2\pi j$.

### The ML inequality

Most contour arguments end by showing that some piece of the contour contributes nothing in a limit. The tool is a simple bound.

> [!theorem] ML inequality (estimation lemma)
> If $\lvert f(z)\rvert\le M$ for all $z$ on a contour $C$ of length $L$, then
>
> $$
> \left\lvert\int_Cf(z)\,dz\right\rvert\le ML .
> $$

*Proof.* $\left\lvert\int_Cf\,dz\right\rvert\le\int_C\lvert f\rvert\lvert dz\rvert\le M\int_C\lvert dz\rvert=ML$. $\square$

To bound $\lvert f\rvert$ from above we need lower bounds on denominators, usually from the reverse triangle inequality $\lvert z_1+z_2\rvert\ge\big\lvert\lvert z_1\rvert-\lvert z_2\rvert\big\rvert$.

> [!example] Worked Example 18.14 — bounding integrals
> (a) Let $C_R$ be the upper semicircle $\lvert z\rvert=R>1$. Show $\int_{C_R}\dfrac{dz}{z^2+1}\to0$ as $R\to\infty$.
> (b) Bound $\left\lvert\oint_{\lvert z\rvert=2}\dfrac{e^z}{z^2+1}dz\right\rvert$ and compare with the exact value.
>
> **(a)** On $C_R$, $\lvert z^2+1\rvert\ge\lvert z\rvert^2-1=R^2-1$, so $M=1/(R^2-1)$; $L=\pi R$. Hence $\left\lvert\int_{C_R}\right\rvert\le\dfrac{\pi R}{R^2-1}\to0$. For $R=10$ the bound is $0.0317$. The general rule: if $\lvert f\rvert\le K/R^k$ on $C_R$ with $k>1$, the integral over $C_R$ vanishes as $R\to\infty$.
>
> **(b)** On $\lvert z\rvert=2$: $\lvert e^z\rvert=e^x\le e^2$ and $\lvert z^2+1\rvert\ge4-1=3$, so $M=e^2/3$ and $L=4\pi$: the bound is $4\pi e^2/3=30.95$. The exact value (Section 18.8) is $2\pi j\sin1$, of modulus $5.287$. The ML bound is crude — its job is to prove that something tends to zero, not to estimate it accurately.

> [!warning] Common Pitfall — forgetting $z'(t)$ or the orientation
> $\int_Cf\,dz$ is not $\int f\,dt$: the factor $z'(t)$ (on circles, $j\rho e^{j\theta}$) carries both the length element and the direction. A clockwise circle gives the negative of the counter-clockwise result, so $\oint dz/z=-2\pi j$ clockwise. Also, $\int_Cf\,dz$ is not the same as the arc-length integral $\int_Cf\lvert dz\rvert$; only the latter is used in bounds.

> [!tip] Engineering Note — contour integrals you have already met
> The inverse Z-transform $x[n]=\frac{1}{2\pi j}\oint_Cz^{n-1}X(z)\,dz$ ([Ch. 17](17-z-transforms.md)) is a contour integral around a circle in the ROC; on the unit circle $z=e^{j\Omega}$, $dz=jz\,d\Omega$ and it collapses to the inverse DTFT $\frac{1}{2\pi}\int_{-\pi}^{\pi}X(e^{j\Omega})e^{j\Omega n}d\Omega$. The inverse Fourier transform $\frac1{2\pi}\int X(\omega)e^{j\omega t}d\omega$ is $\frac{1}{2\pi j}\int_{j\mathbb{R}}X(s)e^{st}ds$ along the imaginary axis. Seeing these as contour integrals is what lets us evaluate them by residues.

## 18.5 Cauchy's integral theorem

### Motivation

Example 18.12 showed an analytic integrand giving the same answer along two different paths. Cauchy's theorem says this is always so on a domain without holes: the integral of an analytic function around any closed curve is zero. Its corollary, the Cauchy integral formula, says something even stronger — the values of an analytic function *inside* a closed curve are completely determined by its values *on* the curve. These two results convert integration problems into the bookkeeping of singular points, and they deliver, almost for free, the infinite differentiability of analytic functions, Liouville's theorem, the fundamental theorem of algebra and the maximum-modulus principle.

### The theorem

> [!theorem] Cauchy's integral theorem (Cauchy–Goursat)
> Let $f$ be analytic in a simply connected domain $D$. Then for every closed contour $C$ in $D$,
>
> $$
> \oint_Cf(z)\,dz=0 .
> $$
>
> Equivalently: if $f$ is analytic inside and on a simple closed contour $C$, then $\oint_Cf\,dz=0$.

*Proof sketch (assuming $f'$ continuous).* Let $C$ be simple and positively oriented, enclosing a region $R$. Split the integral as in Section 18.4 and apply Green's theorem ([Ch. 7](../part-2-calculus/07-vector-calculus.md)), $\oint_C(P\,dx+Q\,dy)=\iint_R(Q_x-P_y)\,dA$:

$$
\oint_Cf\,dz=\oint_C(u\,dx-v\,dy)+j\oint_C(v\,dx+u\,dy)=\iint_R(-v_x-u_y)\,dA+j\iint_R(u_x-v_y)\,dA .
$$

Both integrands vanish identically by the Cauchy–Riemann equations. $\square$

Goursat's refinement removes the assumption that $f'$ is continuous: one subdivides a triangle repeatedly, uses the fact that $f$ is approximated by $f(z_0)+f'(z_0)(z-z_0)$ near a limit point $z_0$, and shows the integral is smaller than any $\varepsilon$. The distinction matters theoretically — it is why "analytic" alone suffices in every later theorem — but not computationally.

> [!theorem] Consequences on a simply connected domain
> If $f$ is analytic in a simply connected domain $D$:
> 1. $\int_{z_1}^{z_2}f\,dz$ is independent of the path in $D$;
> 2. $F(z)=\int_{z_0}^{z}f(\zeta)\,d\zeta$ is analytic in $D$ with $F'=f$ — every analytic function has an antiderivative there.

### Deformation of contours

On a domain with holes, closed curves around a hole need not give zero — $\oint dz/z=2\pi j$ around the origin. But the answer depends only on *which holes are encircled*, not on the shape of the curve.

> [!theorem] Deformation principle
> Let $C$ and $C_1$ be positively oriented simple closed contours with $C_1$ inside $C$, and let $f$ be analytic in the closed region between them (it may have singularities inside $C_1$). Then
>
> $$
> \oint_Cf\,dz=\oint_{C_1}f\,dz .
> $$
>
> More generally, if $C$ encloses disjoint contours $C_1,\dots,C_n$ and $f$ is analytic in the region inside $C$ and outside every $C_k$, then $\oint_Cf\,dz=\sum_{k=1}^n\oint_{C_k}f\,dz$.

*Proof sketch.* Join $C$ to $C_1$ by a cross-cut. The contour that runs around $C$ counter-clockwise, in along the cut, around $C_1$ clockwise and back out along the cut bounds a simply connected region where $f$ is analytic, so its integral is $0$. The two traversals of the cut cancel, leaving $\oint_C-\oint_{C_1}=0$. $\square$

The practical upshot: **shrink the contour onto the singularities.** Any contour around $z_0$ can be replaced by a tiny circle around $z_0$. Combined with Example 18.11, for any positively oriented simple closed contour $C$,

$$
\oint_C\frac{dz}{(z-a)^n}=\begin{cases}2\pi j,&n=1\text{ and }a\text{ inside }C\\0,&\text{otherwise }(n\ge1\text{ integer},\ a\notin C).\end{cases}
$$

> [!example] Worked Example 18.15 — partial fractions plus deformation
> Evaluate $\oint_C\dfrac{dz}{z^2-1}$ for (a) $C:\lvert z\rvert=2$, (b) $C:\lvert z-1\rvert=1$, (c) $C:\lvert z\rvert=\tfrac12$.
>
> Partial fractions: $\dfrac{1}{z^2-1}=\dfrac12\left(\dfrac{1}{z-1}-\dfrac{1}{z+1}\right)$.
>
> **(a)** Both $\pm1$ are inside: $\tfrac12(2\pi j-2\pi j)=0$. (The integral is zero even though the integrand is not analytic inside — Cauchy's theorem is sufficient, not necessary.)
>
> **(b)** Only $z=1$ is inside ($\lvert-1-1\rvert=2>1$): $\tfrac12\cdot2\pi j=\pi j$.
>
> **(c)** Neither is inside; $f$ is analytic on and inside $C$: $0$.

### The Cauchy integral formula

> [!theorem] Cauchy integral formula
> Let $f$ be analytic inside and on a positively oriented simple closed contour $C$, and let $z_0$ be inside $C$. Then
>
> $$
> f(z_0)=\frac{1}{2\pi j}\oint_C\frac{f(z)}{z-z_0}\,dz .
> $$

*Proof.* The integrand $f(z)/(z-z_0)$ is analytic except at $z_0$, so by deformation we may replace $C$ by a small circle $C_\rho:\lvert z-z_0\rvert=\rho$. Write

$$
\oint_{C_\rho}\frac{f(z)}{z-z_0}dz=f(z_0)\oint_{C_\rho}\frac{dz}{z-z_0}+\oint_{C_\rho}\frac{f(z)-f(z_0)}{z-z_0}dz=2\pi j\,f(z_0)+I_\rho .
$$

By continuity, $\lvert f(z)-f(z_0)\rvert<\varepsilon$ on $C_\rho$ for small $\rho$, so by ML, $\lvert I_\rho\rvert\le(\varepsilon/\rho)(2\pi\rho)=2\pi\varepsilon$. Since the left side does not depend on $\rho$, $I_\rho=0$. $\square$

This is astonishing from a real-variable viewpoint: boundary values determine the interior completely. Taking $C$ to be a circle $z=z_0+\rho e^{j\theta}$ gives the **mean-value property**

$$
f(z_0)=\frac{1}{2\pi}\int_0^{2\pi}f(z_0+\rho e^{j\theta})\,d\theta ,
$$

and its real part says that a harmonic function (an electrostatic potential in a charge-free region) equals its average over any circle centred at the point.

> [!theorem] Cauchy integral formula for derivatives
> Under the same hypotheses $f$ has derivatives of all orders at $z_0$, and
>
> $$
> f^{(n)}(z_0)=\frac{n!}{2\pi j}\oint_C\frac{f(z)}{(z-z_0)^{n+1}}\,dz,\qquad n=0,1,2,\dots
> $$

*Proof sketch.* Differentiate the integral formula with respect to $z_0$ under the integral sign: $\frac{d}{dz_0}(z-z_0)^{-1}=(z-z_0)^{-2}$, $\frac{d^2}{dz_0^2}(z-z_0)^{-1}=2(z-z_0)^{-3}$, and so on. Justifying the interchange is a uniform-convergence argument on $C$, where $\lvert z-z_0\rvert$ is bounded below. $\square$

**Corollary.** If $f$ is analytic in $D$, so are $f'$, $f''$, …. Hence $u$ and $v$ have continuous partial derivatives of all orders — the assumption made in proving harmonicity in Section 18.3 was automatically satisfied.

Read in reverse, the derivative formula is an *integration* rule, the form in which we use it:

$$
\oint_C\frac{f(z)}{(z-z_0)^{n+1}}dz=\frac{2\pi j}{n!}f^{(n)}(z_0).
$$

> [!example] Worked Example 18.16 — direct use of the integral formula
> (a) $\displaystyle\oint_{\lvert z\rvert=2}\frac{e^z}{z-1}dz$. (b) $\displaystyle\oint_{\lvert z-j\rvert=1}\frac{dz}{z^2+1}$. (c) $\displaystyle\oint_{\lvert z\rvert=3}\frac{e^z}{z^2+1}dz$.
>
> **(a)** $f=e^z$ is entire and $z_0=1$ is inside: $2\pi j\,e^1=17.079j$.
>
> **(b)** Only $z=j$ is inside the circle $\lvert z-j\rvert=1$ ($-j$ is at distance $2$). Write $\dfrac{1}{z^2+1}=\dfrac{1/(z+j)}{z-j}$, with $f(z)=1/(z+j)$ analytic inside. Result: $2\pi j\cdot\dfrac{1}{2j}=\pi$.
>
> **(c)** Both $\pm j$ are inside. Surround each by a small circle (deformation principle) and apply the formula to each:
>
> $$
> \oint=2\pi j\left[\frac{e^{j}}{2j}+\frac{e^{-j}}{-2j}\right]=2\pi j\cdot\frac{e^j-e^{-j}}{2j}=2\pi j\sin1=5.2871j .
> $$
>
> All three values were confirmed by numerical quadrature of the parametrized integrals.

> [!example] Worked Example 18.17 — higher-order formula
> (a) $\displaystyle\oint_{\lvert z\rvert=2}\frac{e^{2z}}{(z-1)^2}dz$. (b) $\displaystyle\oint_{\lvert z\rvert=2}\frac{\cos z}{(z-j)^3}dz$. (c) $\displaystyle\oint_{\lvert z\rvert=1}\frac{\sin z}{z^4}dz$.
>
> **(a)** $n=1$, $f=e^{2z}$, $f'(1)=2e^2$: result $2\pi j\cdot2e^2=4\pi e^2j=92.854j$.
>
> **(b)** $n=2$, $f=\cos z$, $f''(z)=-\cos z$, $f''(j)=-\cos j=-\cosh1$: result $\dfrac{2\pi j}{2!}(-\cosh1)=-\pi j\cosh1=-4.8477j$.
>
> **(c)** $n=3$, $f=\sin z$, $f'''(0)=-\cos0=-1$: result $\dfrac{2\pi j}{3!}(-1)=-\dfrac{\pi j}{3}=-1.0472j$. Check with the series of Section 18.6: $\sin z/z^4=z^{-3}-\tfrac16z^{-1}+\dots$, and only the $z^{-1}$ term contributes, giving $2\pi j(-\tfrac16)$. ✓

### Consequences: estimates, Liouville, and the maximum principle

**Cauchy's estimate.** If $\lvert f\rvert\le M$ on the circle $\lvert z-z_0\rvert=\rho$, the derivative formula and ML give

$$
\lvert f^{(n)}(z_0)\rvert\le\frac{n!}{2\pi}\cdot\frac{M}{\rho^{n+1}}\cdot2\pi\rho=\frac{n!\,M}{\rho^n}.
$$

> [!theorem] Liouville's theorem
> A bounded entire function is constant.

*Proof.* If $\lvert f\rvert\le M$ everywhere, Cauchy's estimate with $n=1$ gives $\lvert f'(z_0)\rvert\le M/\rho$ for every $\rho$; let $\rho\to\infty$. So $f'\equiv0$. $\square$

That is why $\sin z$ must be unbounded in $\mathbb{C}$. It also gives a two-line proof of the **fundamental theorem of algebra**: if a non-constant polynomial $P$ had no zero, $1/P$ would be entire and (since $\lvert P\rvert\to\infty$) bounded, hence constant — contradiction. So every polynomial of degree $n\ge1$ has a root, and by induction exactly $n$ roots counted with multiplicity. Every characteristic polynomial in this book therefore factors completely over $\mathbb{C}$.

> [!theorem] Maximum-modulus principle
> If $f$ is analytic and non-constant in a bounded domain $D$ and continuous on its closure, then $\lvert f\rvert$ attains its maximum on the boundary of $D$, never at an interior point.

*Proof idea.* By the mean-value property, $\lvert f(z_0)\rvert\le$ the average of $\lvert f\rvert$ over any small circle about $z_0$; if $\lvert f(z_0)\rvert$ were a maximum, $\lvert f\rvert$ would have to equal it on the whole circle, and then (by the identity theorem) everywhere. $\square$

> [!tip] Engineering Note — why the $H_\infty$ norm is measured on the $j\omega$-axis
> For a stable transfer function $H(s)$, analytic in the closed right half-plane and bounded at infinity, the maximum-modulus principle (applied on large half-discs) shows $\sup_{\operatorname{Re}s\ge0}\lvert H(s)\rvert=\sup_\omega\lvert H(j\omega)\rvert$. The peak gain over all exponentially growing inputs equals the peak of the Bode magnitude plot. This is the $H_\infty$ norm of robust control, and it is computable from frequency-response data alone.

> [!warning] Common Pitfall — hypotheses of Cauchy's theorem
> $\oint_Cf\,dz=0$ requires analyticity at **every** point inside $C$, not just on $C$. $\oint_{\lvert z\rvert=1}e^z/z^2\,dz=2\pi j\ne0$ even though $e^z/z^2$ is perfectly well behaved on the circle. Equally, a singularity *on* the contour makes the integral improper; it must be treated as a principal value with an indentation (Section 18.8), not ignored.

## 18.6 Taylor and Laurent series

### Motivation

Power series are how analytic functions are computed, approximated and classified. A Taylor series describes $f$ near a point where it is analytic. A Laurent series — a power series allowed to have negative powers — describes $f$ near, and around, its singularities, and it is valid in an **annulus** rather than a disc. Engineers already use Laurent series under another name: the Z-transform $X(z)=\sum_nx[n]z^{-n}$ *is* a Laurent series in $z$, and its region of convergence is the annulus of that series. The coefficient of $(z-z_0)^{-1}$ in a Laurent series is the *residue*, the quantity on which Sections 18.7–18.9 are built.

### Power series

> [!theorem] Convergence of power series
> For every power series $\sum_{n=0}^\infty a_n(z-z_0)^n$ there is a **radius of convergence** $R\in[0,\infty]$ such that the series converges absolutely for $\lvert z-z_0\rvert<R$ and diverges for $\lvert z-z_0\rvert>R$. It is given by
>
> $$
> \frac1R=\limsup_{n\to\infty}\lvert a_n\rvert^{1/n}\qquad\left(\text{or }R=\lim_{n\to\infty}\left\lvert\frac{a_n}{a_{n+1}}\right\rvert\text{ when this limit exists}\right).
> $$
>
> Inside the disc of convergence the sum is analytic and may be differentiated or integrated term by term, the new series having the same radius.

The region of convergence is always a disc: there is no such thing as a power series converging on an interval but not on the corresponding disc. On the circle $\lvert z-z_0\rvert=R$ itself anything can happen, point by point.

### Taylor's theorem

> [!theorem] Taylor's theorem
> Let $f$ be analytic in the disc $\lvert z-z_0\rvert<R$. Then for every $z$ in the disc
>
> $$
> f(z)=\sum_{n=0}^\infty a_n(z-z_0)^n,\qquad a_n=\frac{f^{(n)}(z_0)}{n!}=\frac{1}{2\pi j}\oint_C\frac{f(\zeta)}{(\zeta-z_0)^{n+1}}d\zeta ,
> $$
>
> where $C$ is any circle $\lvert\zeta-z_0\rvert=\rho<R$. The series converges in the **largest open disc about $z_0$ in which $f$ is analytic**: its radius is the distance from $z_0$ to the nearest singular point of $f$.

*Proof.* Fix $z$ with $\lvert z-z_0\rvert<\rho<R$ and use the Cauchy formula on the circle $C$ of radius $\rho$. For $\zeta$ on $C$, $\lvert(z-z_0)/(\zeta-z_0)\rvert<1$, so the Cauchy kernel expands as a geometric series:

$$
\frac{1}{\zeta-z}=\frac{1}{(\zeta-z_0)\left(1-\dfrac{z-z_0}{\zeta-z_0}\right)}=\sum_{n=0}^\infty\frac{(z-z_0)^n}{(\zeta-z_0)^{n+1}} .
$$

The convergence is uniform on $C$, so we may integrate term by term:

$$
f(z)=\frac{1}{2\pi j}\oint_C\frac{f(\zeta)}{\zeta-z}d\zeta=\sum_{n=0}^\infty\left[\frac{1}{2\pi j}\oint_C\frac{f(\zeta)}{(\zeta-z_0)^{n+1}}d\zeta\right](z-z_0)^n ,
$$

and the bracket is $f^{(n)}(z_0)/n!$ by the derivative formula. $\square$

This is the payoff of complex differentiability: *analytic* and *equal to a convergent power series* are the same property. Standard Maclaurin series (all valid in $\mathbb{C}$ with the stated radius):

| Function | Series | Radius |
|---|---|---|
| $\dfrac{1}{1-z}$ | $\sum_{n\ge0}z^n$ | $1$ |
| $e^z$ | $\sum_{n\ge0}\dfrac{z^n}{n!}$ | $\infty$ |
| $\sin z$ | $\sum_{n\ge0}\dfrac{(-1)^nz^{2n+1}}{(2n+1)!}$ | $\infty$ |
| $\cos z$ | $\sum_{n\ge0}\dfrac{(-1)^nz^{2n}}{(2n)!}$ | $\infty$ |
| $\operatorname{Log}(1+z)$ | $\sum_{n\ge1}\dfrac{(-1)^{n+1}z^n}{n}$ | $1$ |
| $(1+z)^\alpha$ (principal) | $\sum_{n\ge0}\dbinom{\alpha}{n}z^n$ | $1$ (unless $\alpha\in\mathbb{N}$) |

> [!example] Worked Example 18.18 — Taylor series and their radii
> (a) Expand $f(z)=1/z$ about $z_0=1+j$ and give the radius. (b) Find the first three nonzero terms of $\tan z$ about $0$ and the radius. (c) Expand $g(z)=\dfrac{1}{1+z^2}$ about $z_0=1$ to second order.
>
> **(a)** Write $z=z_0+(z-z_0)$: $\dfrac1z=\dfrac{1}{z_0}\cdot\dfrac{1}{1+(z-z_0)/z_0}=\sum_{n\ge0}\dfrac{(-1)^n(z-z_0)^n}{z_0^{n+1}}$, valid for $\lvert z-z_0\rvert<\lvert z_0\rvert=\sqrt2$. The only singularity is $z=0$, at distance $\sqrt2$, as the theorem predicts.
>
> **(b)** Divide series: $\tan z=\dfrac{z-z^3/6+z^5/120-\cdots}{1-z^2/2+z^4/24-\cdots}$. Set $\tan z=z+bz^3+cz^5$; multiplying by the denominator and matching: $z^3$: $b-\tfrac12=-\tfrac16\Rightarrow b=\tfrac13$; $z^5$: $c-\tfrac b2+\tfrac1{24}=\tfrac1{120}\Rightarrow c=\tfrac{2}{15}$. So $\tan z=z+\tfrac13z^3+\tfrac{2}{15}z^5+\cdots$. The nearest zeros of $\cos z$ are $\pm\pi/2$, so $R=\pi/2$. Check at $z=0.3$: $0.3+0.009+0.000324=0.309324$, vs $\tan0.3=0.309336$.
>
> **(c)** $g(1)=\tfrac12$; $g'=-2z/(1+z^2)^2$, $g'(1)=-\tfrac12$; $g''=(6z^2-2)/(1+z^2)^3$, $g''(1)=\tfrac12$. So $g\approx\tfrac12-\tfrac12(z-1)+\tfrac14(z-1)^2$, with radius $\lvert1\mp j\rvert=\sqrt2$ (Example 18.7).

**Zeros.** If $f$ is analytic at $z_0$, $f(z_0)=0$ and not all derivatives vanish there, the Taylor series begins $f=a_m(z-z_0)^m+\cdots$ with $a_m\neq0$: $z_0$ is a **zero of order $m$**, and $f(z)=(z-z_0)^mg(z)$ with $g$ analytic and $g(z_0)\neq0$. Since $g\ne0$ near $z_0$, **zeros of a non-constant analytic function are isolated** — the root of the identity theorem.

### Laurent series

Near a singularity a Taylor series is impossible, but an expansion in both positive and negative powers works.

> [!theorem] Laurent's theorem
> Let $f$ be analytic in the annulus $A:\ r_1<\lvert z-z_0\rvert<r_2$ (with $0\le r_1<r_2\le\infty$). Then for every $z$ in $A$
>
> $$
> f(z)=\sum_{n=-\infty}^{\infty}a_n(z-z_0)^n=\underbrace{\sum_{n=0}^\infty a_n(z-z_0)^n}_{\text{analytic (regular) part}}+\underbrace{\sum_{n=1}^\infty\frac{a_{-n}}{(z-z_0)^n}}_{\text{principal part}},
> $$
>
> with
>
> $$
> a_n=\frac{1}{2\pi j}\oint_C\frac{f(\zeta)}{(\zeta-z_0)^{n+1}}d\zeta,\qquad n\in\mathbb{Z},
> $$
>
> $C$ being any positively oriented circle $\lvert\zeta-z_0\rvert=\rho$ with $r_1<\rho<r_2$. The expansion in a given annulus is **unique**.

*Proof sketch.* For $z$ in the annulus take circles $C_2$ (radius just below $r_2$) and $C_1$ (just above $r_1$) on either side of $z$. Cauchy's formula on the region between them gives $f(z)=\frac{1}{2\pi j}\oint_{C_2}\frac{f(\zeta)}{\zeta-z}d\zeta-\frac{1}{2\pi j}\oint_{C_1}\frac{f(\zeta)}{\zeta-z}d\zeta$. On $C_2$, $\lvert z-z_0\rvert<\lvert\zeta-z_0\rvert$, and the kernel expands in nonnegative powers exactly as in Taylor's proof. On $C_1$, $\lvert\zeta-z_0\rvert<\lvert z-z_0\rvert$, so instead

$$
-\frac{1}{\zeta-z}=\frac{1}{(z-z_0)\left(1-\dfrac{\zeta-z_0}{z-z_0}\right)}=\sum_{n=1}^\infty\frac{(\zeta-z_0)^{n-1}}{(z-z_0)^n},
$$

which produces the negative powers. Both circles can be deformed to any circle $C$ in the annulus. $\square$

Note that the coefficient formula with $n=-1$ reads

$$
a_{-1}=\frac{1}{2\pi j}\oint_Cf(\zeta)\,d\zeta :
$$

the integral of $f$ around the annulus is $2\pi j$ times the coefficient of $(z-z_0)^{-1}$. That is the content of Example 18.11 seen from the other side, and it is the definition of the residue in Section 18.7.

**In practice, nobody computes Laurent coefficients from the integral formula.** Because the expansion is unique, *any* legitimate manipulation — geometric series, known Taylor series, substitution, multiplication, partial fractions — that produces a series valid in the annulus gives *the* Laurent series. The essential tool is the geometric series used in two ways:

$$
\frac{1}{1-w}=\sum_{n\ge0}w^n\ \ (\lvert w\rvert<1),\qquad \frac{1}{1-w}=-\frac1w\cdot\frac{1}{1-1/w}=-\sum_{n\ge1}w^{-n}\ \ (\lvert w\rvert>1).
$$

> [!example] Worked Example 18.19 — one function, four expansions
> Expand $f(z)=\dfrac{1}{(z-1)(z-2)}$ in (a) $\lvert z\rvert<1$, (b) $1<\lvert z\rvert<2$, (c) $\lvert z\rvert>2$, (d) $0<\lvert z-1\rvert<1$.
>
> Partial fractions: $f=\dfrac{1}{z-2}-\dfrac{1}{z-1}$. The singularities $1$ and $2$ cut the plane about $z_0=0$ into three regions; in each we expand each fraction in the convergent direction.
>
> **(a)** $\lvert z\rvert<1$ (so also $\lvert z/2\rvert<1$):
>
> $$
> \frac{1}{z-2}=-\frac12\cdot\frac{1}{1-z/2}=-\sum_{n\ge0}\frac{z^n}{2^{n+1}},\qquad -\frac{1}{z-1}=\frac{1}{1-z}=\sum_{n\ge0}z^n,
> $$
>
> $$
> f=\sum_{n\ge0}\left(1-2^{-(n+1)}\right)z^n=\tfrac12+\tfrac34z+\tfrac78z^2+\cdots\quad(\text{a Taylor series}).
> $$
>
> **(b)** $1<\lvert z\rvert<2$: keep the first expansion; for the second use $\lvert1/z\rvert<1$:
>
> $$
> -\frac{1}{z-1}=-\frac1z\cdot\frac{1}{1-1/z}=-\sum_{n\ge1}z^{-n},\qquad f=-\sum_{n\ge0}\frac{z^n}{2^{n+1}}-\sum_{n\ge1}\frac{1}{z^n}.
> $$
>
> Check at $z=1.5$: the first sum is $-\tfrac12/(1-0.75)=-2$, the second $-\tfrac{2/3}{1-2/3}=-2$; total $-4$, and indeed $f(1.5)=1/(0.5\cdot(-0.5))=-4$. ✓
>
> **(c)** $\lvert z\rvert>2$: now $\dfrac{1}{z-2}=\dfrac1z\cdot\dfrac{1}{1-2/z}=\sum_{n\ge1}2^{n-1}z^{-n}$, so
>
> $$
> f=\sum_{n\ge1}\left(2^{n-1}-1\right)z^{-n}=\frac{1}{z^2}+\frac{3}{z^3}+\frac{7}{z^4}+\cdots .
> $$
>
> **(d)** About $z_0=1$, put $w=z-1$: $f=\dfrac{1}{w(w-1)}=-\dfrac1w\cdot\dfrac{1}{1-w}=-\dfrac1w-1-w-w^2-\cdots$ for $0<\lvert w\rvert<1$. The coefficient of $(z-1)^{-1}$ is $-1$.
>
> Notice that in (c) the coefficient of $z^{-1}$ is $2^0-1=0$: by the coefficient formula $\oint_{\lvert z\rvert=3}f\,dz=0$, which agrees with Example 18.15-style reasoning ($\tfrac{1}{z-2}$ and $-\tfrac1{z-1}$ contribute $+2\pi j$ and $-2\pi j$). In (b) the coefficient of $z^{-1}$ is $-1$, so $\oint_{\lvert z\rvert=1.5}f\,dz=-2\pi j$: only the pole at $1$ is enclosed.

```desmos-graph
left=-3; right=3; bottom=-3; top=3
---
x^2+y^2=1|blue|dashed
x^2+y^2=4|red|dashed
(1,0)|label:pole
(2,0)|label:pole
(0,0)|label:centre
```
*Regions of validity for the expansions of $1/((z-1)(z-2))$ about $0$: a Taylor series inside the blue circle $\lvert z\rvert=1$, a two-sided Laurent series in the annulus between the circles, and a series in negative powers only outside the red circle $\lvert z\rvert=2$. Every boundary circle passes through a pole.*

> [!tip] Engineering Note — the Z-transform ROC is a Laurent annulus
> Put $X(z)=\sum_nx[n]z^{-n}$: the sequence values are the Laurent coefficients with $a_{-n}=x[n]$. Example 18.19 is therefore a complete catalogue of the three sequences sharing $X(z)=1/((z-1)(z-2))$ ([Ch. 17](17-z-transforms.md)): region (c) gives a **right-sided (causal)** sequence $x[n]=2^{n-1}-1$ for $n\ge1$; region (a) gives a **left-sided (anti-causal)** sequence; region (b) a **two-sided** one. The ROC is always an annulus bounded by poles, and the algebraic expression $X(z)$ alone does not determine the signal — the annulus does. A system is stable when its ROC contains the unit circle.

### Isolated singularities

> [!definition] Classification of isolated singularities
> $z_0$ is an **isolated singularity** of $f$ if $f$ is analytic in a punctured disc $0<\lvert z-z_0\rvert<\rho$. Let $\sum a_n(z-z_0)^n$ be the Laurent series there. Then $z_0$ is
> - a **removable singularity** if the principal part is zero ($a_{-n}=0$ for all $n\ge1$);
> - a **pole of order $m$** if the principal part has finitely many terms with $a_{-m}\neq0$ the last: $f=\dfrac{a_{-m}}{(z-z_0)^m}+\dots+\dfrac{a_{-1}}{z-z_0}+\sum_{n\ge0}a_n(z-z_0)^n$. A pole of order $1$ is a **simple pole**;
> - an **essential singularity** if the principal part has infinitely many nonzero terms.

Practical tests, which avoid computing series:

| Type | Behaviour as $z\to z_0$ | Test |
|---|---|---|
| removable | $f$ bounded, has a finite limit | $\lim(z-z_0)f(z)=0$ |
| pole of order $m$ | $\lvert f\rvert\to\infty$ | $(z-z_0)^mf(z)\to$ finite nonzero limit |
| essential | no limit, finite or infinite | neither of the above |

For a quotient $f=g/h$ with $g$, $h$ analytic: if $h$ has a zero of order $m$ at $z_0$ and $g(z_0)\neq0$, then $z_0$ is a pole of order $m$. More generally, if $g$ has a zero of order $k$ and $h$ of order $m$, the result is a pole of order $m-k$ if $m>k$, otherwise a removable singularity.

Near an essential singularity the behaviour is wild: by the Casorati–Weierstrass theorem, $f$ comes arbitrarily close to *every* complex value in every punctured neighbourhood (Picard's theorem says it actually takes every value, with at most one exception, infinitely often).

> [!example] Worked Example 18.20 — classifying singularities by Laurent series
> Classify the singularity at $z=0$ and find the coefficient $a_{-1}$ for (a) $\dfrac{\sin z}{z}$, (b) $\dfrac{1-\cos z}{z^4}$, (c) $\dfrac{e^z}{z^3}$, (d) $e^{1/z}$, (e) $z^2\sin(1/z)$.
>
> **(a)** $\dfrac{\sin z}{z}=1-\dfrac{z^2}{6}+\dfrac{z^4}{120}-\cdots$: removable (define the value $1$ at $0$ — the sinc function is entire). $a_{-1}=0$.
>
> **(b)** $1-\cos z=\dfrac{z^2}{2}-\dfrac{z^4}{24}+\dfrac{z^6}{720}-\cdots$, so $\dfrac{1-\cos z}{z^4}=\dfrac{1}{2z^2}-\dfrac1{24}+\dfrac{z^2}{720}-\cdots$: a pole of order $2$ (numerator zero of order 2, denominator of order 4). $a_{-1}=0$.
>
> **(c)** $\dfrac{e^z}{z^3}=\dfrac{1}{z^3}+\dfrac{1}{z^2}+\dfrac{1}{2z}+\dfrac16+\cdots$: pole of order $3$, $a_{-1}=\tfrac12$.
>
> **(d)** Substitute $1/z$ into the exponential series: $e^{1/z}=1+\dfrac1z+\dfrac{1}{2!z^2}+\dfrac{1}{3!z^3}+\cdots$, valid for $\lvert z\rvert>0$. Infinitely many negative powers: essential. $a_{-1}=1$. Along the positive real axis $e^{1/x}\to\infty$; along the negative real axis it $\to0$; along the imaginary axis $\lvert e^{1/(jy)}\rvert=1$.
>
> **(e)** $z^2\left(\dfrac1z-\dfrac{1}{6z^3}+\dfrac{1}{120z^5}-\cdots\right)=z-\dfrac{1}{6z}+\dfrac{1}{120z^3}-\cdots$: essential, $a_{-1}=-\tfrac16$.

**Non-isolated singularities.** Branch points such as $z=0$ for $\operatorname{Log}z$ or $z^{1/2}$ are not isolated singularities of any single-valued branch (the whole cut is singular), and have no Laurent expansion about them. The function $1/\sin(1/z)$ has poles at $z=1/(k\pi)$ accumulating at $0$, so $0$ is a non-isolated singularity. Rational functions — and so every lumped-element transfer function — have only poles.

> [!warning] Common Pitfall — the "right" Laurent series
> A function has a *different* Laurent series in each annulus. Asking for "the Laurent series of $f$ about $z_0$" is meaningless until the annulus is specified. Residues are read from the series in the punctured disc $0<\lvert z-z_0\rvert<\rho$ around the singularity itself — region (d) of Example 18.19, not region (b), even though both contain negative powers.

## 18.7 Residues

### Motivation

When a function is integrated around a small circle enclosing one isolated singularity, every term of its Laurent series integrates to zero except one: the $(z-z_0)^{-1}$ term, which contributes $2\pi j\,a_{-1}$. That single coefficient is all the contour "sees" of the singularity. It deserves a name and a toolkit for computing it quickly — ideally without finding the whole Laurent series. The same numbers are, not by coincidence, the coefficients of partial-fraction expansions that you have computed by the cover-up rule in [Ch. 14](14-laplace-transforms.md).

> [!definition] Residue
> Let $z_0$ be an isolated singularity of $f$, with Laurent series $\sum_na_n(z-z_0)^n$ in $0<\lvert z-z_0\rvert<\rho$. The **residue** of $f$ at $z_0$ is
>
> $$
> \operatorname{Res}_{z=z_0}f(z)=a_{-1}=\frac{1}{2\pi j}\oint_Cf(z)\,dz ,
> $$
>
> where $C$ is any positively oriented simple closed contour around $z_0$ inside the punctured disc.

At a removable singularity the residue is $0$. At an essential singularity it must be read from the Laurent series. At poles — the common case — there are formulas.

### Computing residues at poles

> [!theorem] Residue formulas
> **Simple pole.**
>
> $$
> \operatorname{Res}_{z=z_0}f=\lim_{z\to z_0}(z-z_0)f(z).
> $$
>
> **Quotient with a simple zero in the denominator.** If $f=g/h$ with $g,h$ analytic at $z_0$, $g(z_0)\neq0$, $h(z_0)=0$, $h'(z_0)\ne0$, then
>
> $$
> \operatorname{Res}_{z=z_0}\frac{g}{h}=\frac{g(z_0)}{h'(z_0)} .
> $$
>
> **Pole of order $m$.**
>
> $$
> \operatorname{Res}_{z=z_0}f=\frac{1}{(m-1)!}\lim_{z\to z_0}\frac{d^{m-1}}{dz^{m-1}}\Big[(z-z_0)^mf(z)\Big].
> $$

*Derivation.* At a pole of order $m$, $(z-z_0)^mf(z)=a_{-m}+a_{-m+1}(z-z_0)+\dots+a_{-1}(z-z_0)^{m-1}+\dots$ is analytic; $a_{-1}$ is its Taylor coefficient of order $m-1$, i.e. $\frac{1}{(m-1)!}$ times its $(m-1)$-th derivative at $z_0$. For $m=1$ this is the simple-pole formula. For the quotient rule, $h(z)=h'(z_0)(z-z_0)+\dots$, so $(z-z_0)g/h\to g(z_0)/h'(z_0)$. $\square$

The quotient rule $g/h'$ is the fastest route whenever the denominator is not already factored — for instance at the roots of $z^4+1$ or of $\sin z$.

> [!example] Worked Example 18.21 — simple poles
> (a) $f=\dfrac{ze^z}{z^2-1}$ at $z=\pm1$. (b) $f=\dfrac{1}{z^4+1}$ at $z_1=e^{j\pi/4}$. (c) $f=\cot z$ at $z=k\pi$. (d) $f=\dfrac{e^{jz}}{z^2+a^2}$ ($a>0$) at $z=ja$.
>
> **(a)** At $z=1$: $\lim(z-1)\dfrac{ze^z}{(z-1)(z+1)}=\dfrac{e}{2}=1.3591$. At $z=-1$: $\dfrac{(-1)e^{-1}}{-2}=\dfrac{1}{2e}=0.1839$.
>
> **(b)** Quotient rule: $\dfrac{1}{4z_1^3}=\dfrac{z_1}{4z_1^4}=-\dfrac{z_1}{4}=-\dfrac{e^{j\pi/4}}{4}=-0.17678-0.17678j$ (using $z_1^4=-1$).
>
> **(c)** $\cot z=\cos z/\sin z$, $\sin(k\pi)=0$, $(\sin z)'=\cos z$: residue $\cos(k\pi)/\cos(k\pi)=1$ at every integer multiple of $\pi$. (This is why $\pi\cot\pi z$ is used to sum series $\sum f(n)$: it has residue $1$ at every integer.)
>
> **(d)** $z^2+a^2=(z-ja)(z+ja)$: residue $\dfrac{e^{j(ja)}}{2ja}=\dfrac{e^{-a}}{2ja}$. We will need this in Section 18.8.

> [!example] Worked Example 18.22 — higher-order poles
> (a) $f=\dfrac{e^{2z}}{(z-1)^2}$ at $z=1$. (b) $f=\dfrac{1}{(z^2+1)^2}$ at $z=j$. (c) $f=\dfrac{1}{z^3(z-1)}$ at $z=0$ and $z=1$.
>
> **(a)** $m=2$: $\dfrac{d}{dz}e^{2z}\Big\rvert_{z=1}=2e^2=14.778$ — the same number as Example 18.17(a), as it must be.
>
> **(b)** $(z^2+1)^2=(z-j)^2(z+j)^2$, so $m=2$ and we differentiate $(z+j)^{-2}$: $-2(z+j)^{-3}$ at $z=j$ gives $-2/(2j)^3=-2/(-8j)=\dfrac{1}{4j}=-\dfrac{j}{4}$.
>
> **(c)** At $z=0$ ($m=3$): $\dfrac{1}{2!}\dfrac{d^2}{dz^2}\left(\dfrac{1}{z-1}\right)\Big\rvert_0=\dfrac12\cdot\dfrac{2}{(z-1)^3}\Big\rvert_0=-1$. Alternatively, from the series $-\dfrac{1}{z^3}\cdot\dfrac{1}{1-z}=-z^{-3}-z^{-2}-z^{-1}-\cdots$, $a_{-1}=-1$. At $z=1$ (simple): $1/1^3=1$. The two residues sum to $0$ — see the residue at infinity below.

### Residues and partial fractions

For a proper rational function $F(s)=N(s)/D(s)$ with simple poles $p_k$, the partial-fraction expansion $F=\sum_kA_k/(s-p_k)$ has

$$
A_k=\lim_{s\to p_k}(s-p_k)F(s)=\operatorname{Res}_{s=p_k}F(s),
$$

because $A_k/(s-p_k)$ is exactly the principal part of $F$ at $p_k$. In general, **a proper rational function equals the sum of its principal parts at all its poles** (the difference would be an entire function tending to zero at infinity, hence zero by Liouville). The Heaviside cover-up rule is the simple-pole residue formula. At a repeated pole of order $m$, the coefficient of $(s-p)^{-k}$ is

$$
A_{k}=\frac{1}{(m-k)!}\frac{d^{m-k}}{ds^{m-k}}\Big[(s-p)^mF(s)\Big]_{s=p},
$$

of which the residue is the case $k=1$.

> [!example] Worked Example 18.23 — partial fractions as principal parts
> Expand $F(s)=\dfrac{s+3}{(s+1)(s+2)^2}$.
>
> At $s=-1$ (simple): $A=\dfrac{s+3}{(s+2)^2}\Big\rvert_{-1}=\dfrac{2}{1}=2$.
>
> At $s=-2$ (order 2), with $G(s)=(s+2)^2F=\dfrac{s+3}{s+1}$: the $(s+2)^{-2}$ coefficient is $G(-2)=\dfrac{1}{-1}=-1$; the $(s+2)^{-1}$ coefficient (the residue) is $G'(-2)=\dfrac{(s+1)-(s+3)}{(s+1)^2}\Big\rvert_{-2}=\dfrac{-2}{1}=-2$.
>
> $$
> F(s)=\frac{2}{s+1}-\frac{2}{s+2}-\frac{1}{(s+2)^2}.
> $$
>
> Check at $s=0$: $F(0)=3/4$, and $2-1-\tfrac14=\tfrac34$. ✓ Also, the residues $2$ and $-2$ sum to zero, consistent with $F\sim1/s^2$ at infinity (no $1/s$ term), as the next result explains.

### The residue at infinity

> [!definition] Residue at infinity
> If $f$ is analytic for $\lvert z\rvert>R$, define
>
> $$
> \operatorname{Res}_{z=\infty}f=-\frac{1}{2\pi j}\oint_{\lvert z\rvert=\rho}f(z)\,dz\quad(\rho>R)=-c_{-1},
> $$
>
> where $c_{-1}$ is the coefficient of $1/z$ in the Laurent expansion of $f$ valid for $\lvert z\rvert>R$. Equivalently, $\operatorname{Res}_{z=\infty}f(z)=-\operatorname{Res}_{\zeta=0}\left[\dfrac{1}{\zeta^2}f\!\left(\dfrac1\zeta\right)\right]$.

> [!theorem] Sum of all residues
> If $f$ has finitely many singularities in the extended plane, the sum of its residues at all of them, including $\infty$, is zero.

This is useful when many poles lie inside a contour and few outside: $\oint_Cf\,dz=2\pi j\sum_{\text{inside}}\operatorname{Res}=-2\pi j\left[\sum_{\text{outside}}\operatorname{Res}+\operatorname{Res}_\infty\right]$. For a rational $f=P/Q$ with $\deg Q\ge\deg P+2$, $f=O(1/z^2)$ at infinity, the $1/z$ coefficient vanishes, and **the finite residues sum to zero** — as in Examples 18.22(c) and 18.23.

> [!example] Worked Example 18.24 — many poles inside, use infinity
> Evaluate $\displaystyle\oint_{\lvert z\rvert=2}\frac{z^{15}}{(z^2+1)^4(z^4+2)^3}dz$.
>
> The denominator has twenty roots, all with modulus $1$ or $2^{1/4}=1.189$, all inside the contour. Computing twenty residues (with poles of orders 4 and 3!) is hopeless. Instead: the degree of the denominator is $8+12=20$, so for large $\lvert z\rvert$, $f=\dfrac{z^{15}}{z^{20}(1+\cdots)}=z^{-5}(1+\cdots)$ has no $z^{-1}$ term: $c_{-1}=0$. Hence the integral is $2\pi j\,c_{-1}=0$.
>
> Contrast $\oint_{\lvert z\rvert=2}\dfrac{z^{3}}{z^4+2}dz$: here $f=z^{-1}(1+2z^{-4})^{-1}=z^{-1}-2z^{-5}+\cdots$, so $c_{-1}=1$ and the integral is $2\pi j$ — confirming that the four simple poles each have residue $z_k^3/(4z_k^3)=\tfrac14$.

> [!tip] Engineering Note — residues are the modal amplitudes
> In the inverse Laplace transform of $F(s)$ with simple poles $p_k$, the time function is $\sum_k\operatorname{Res}_{s=p_k}F\cdot e^{p_kt}$ (Section 18.9). Each residue is the amplitude with which the natural mode $e^{p_kt}$ is excited. A small residue at a pole means that mode is weakly excited or weakly observed — the near cancellation of a pole by a nearby zero shows up as a tiny residue, which is why such a pole can be dropped in model reduction.

> [!warning] Common Pitfall — using the simple-pole formula at a double pole
> $\lim(z-z_0)f(z)$ is infinite at a pole of order $2$ or more; it is not "the residue". Determine the order first (factor the denominator, or count the order of the zero). Equally, $g(z_0)/h'(z_0)$ is valid only when $h$ has a *simple* zero; at $z_0=0$ for $f=1/(z\sin z)$, $h=z\sin z$ has a double zero and the residue is $0$ (since $1/(z\sin z)=z^{-2}+\tfrac16+\cdots$), not $1/h'(0)$, which is undefined.

## 18.8 Residue theorem

### The theorem

> [!theorem] Residue theorem
> Let $C$ be a positively oriented simple closed contour, and let $f$ be analytic inside and on $C$ except at finitely many isolated singularities $z_1,\dots,z_n$ inside $C$. Then
>
> $$
> \oint_Cf(z)\,dz=2\pi j\sum_{k=1}^n\operatorname{Res}_{z=z_k}f(z).
> $$

*Proof.* Surround each $z_k$ by a small circle $C_k$ inside $C$, the circles disjoint. By the deformation principle (Section 18.5), $\oint_C=\sum_k\oint_{C_k}$. On $C_k$, integrate the Laurent series of $f$ about $z_k$ term by term: by Example 18.11 every term but $a_{-1}^{(k)}/(z-z_k)$ gives zero, and that one gives $2\pi j\,a_{-1}^{(k)}$. $\square$

Cauchy's theorem (no singularities, sum empty) and Cauchy's integral formula (one simple pole of $f(z)/(z-z_0)$ with residue $f(z_0)$) are special cases. The theorem reduces *any* contour integral of a function with isolated singularities to local algebra at the singular points.

> [!example] Worked Example 18.25 — closed contour integrals
> (a) $\displaystyle\oint_{\lvert z\rvert=2}\frac{ze^z}{z^2-1}dz$. (b) $\displaystyle\oint_{\lvert z\rvert=2}\tan z\,dz$. (c) $\displaystyle\oint_{\lvert z\rvert=1}e^{1/z}\,dz$.
>
> **(a)** Both poles $\pm1$ are inside; by Example 18.21(a), $\oint=2\pi j\left(\tfrac e2+\tfrac1{2e}\right)=2\pi j\cosh1=9.6955j$.
>
> **(b)** Inside $\lvert z\rvert=2$, $\cos z=0$ at $\pm\pi/2$ ($1.571<2$; the next, $3\pi/2$, is outside). Residues: $\sin z/(-\sin z)=-1$ at each. So $\oint=2\pi j(-2)=-4\pi j=-12.566j$.
>
> **(c)** The only singularity is the essential one at $0$, with residue $1$ (Example 18.20(d)): $\oint=2\pi j$. No residue *formula* applies; the Laurent series does.
>
> Each was confirmed by direct numerical quadrature on the parametrized circle.

The rest of this section is about the most celebrated use of the theorem: evaluating **real** definite integrals that resist real-variable methods. The strategy in every case is:

1. choose a complex function whose real or imaginary part, or whose restriction to the real axis, is the integrand;
2. choose a closed contour, part of which is the real interval of interest;
3. evaluate the closed integral by residues;
4. show that the extra parts of the contour contribute zero (or something computable) in a limit — this is where ML and Jordan's lemma come in.

### Type I — trigonometric integrals over a full period

To evaluate $\displaystyle I=\int_0^{2\pi}R(\cos\theta,\sin\theta)\,d\theta$ with $R$ rational, put $z=e^{j\theta}$. As $\theta$ runs over $[0,2\pi]$, $z$ runs once around the unit circle, and

$$
\cos\theta=\frac{z+z^{-1}}{2},\qquad\sin\theta=\frac{z-z^{-1}}{2j},\qquad d\theta=\frac{dz}{jz}.
$$

The integral becomes $\oint_{\lvert z\rvert=1}(\text{rational function of }z)\,dz$, evaluated by the residues **inside the unit circle**.

> [!example] Worked Example 18.26 — $\displaystyle\int_0^{2\pi}\frac{d\theta}{2+\cos\theta}$
> Substitute:
>
> $$
> I=\oint_{\lvert z\rvert=1}\frac{1}{2+\tfrac12(z+z^{-1})}\cdot\frac{dz}{jz}=\oint\frac{2\,dz}{j(z^2+4z+1)} .
> $$
>
> The roots of $z^2+4z+1$ are $z=-2\pm\sqrt3$: $\alpha=-2+\sqrt3=-0.268$ (inside) and $\beta=-2-\sqrt3=-3.732$ (outside). Residue at $\alpha$: $\dfrac{2}{j(\alpha-\beta)}=\dfrac{2}{j\,2\sqrt3}=\dfrac{1}{j\sqrt3}$. Hence
>
> $$
> I=2\pi j\cdot\frac{1}{j\sqrt3}=\frac{2\pi}{\sqrt3}=3.6276 .
> $$
>
> Numerical quadrature: $3.62760$. ✓ The same steps give the general result $\displaystyle\int_0^{2\pi}\frac{d\theta}{a+b\cos\theta}=\frac{2\pi}{\sqrt{a^2-b^2}}$ for $a>\lvert b\rvert$, and the same with $\sin\theta$ in place of $\cos\theta$ (e.g. $\int_0^{2\pi}d\theta/(5-4\sin\theta)=2\pi/3$).

> [!example] Worked Example 18.27 — a cosine-weighted integral
> Evaluate $\displaystyle I=\int_0^{2\pi}\frac{\cos2\theta}{5-4\cos\theta}d\theta$.
>
> It is easier to compute $J=\int_0^{2\pi}\frac{e^{j2\theta}}{5-4\cos\theta}d\theta$ and take $I=\operatorname{Re}J$ (the denominator is real). With $e^{j2\theta}=z^2$ and $5-4\cos\theta=5-2(z+z^{-1})=-\dfrac{2z^2-5z+2}{z}=-\dfrac{(2z-1)(z-2)}{z}$:
>
> $$
> J=\oint_{\lvert z\rvert=1}\frac{z^2\cdot z}{-(2z-1)(z-2)}\cdot\frac{dz}{jz}=\oint\frac{z^2\,dz}{-j(2z-1)(z-2)} .
> $$
>
> Only $z=\tfrac12$ is inside. Writing $2z-1=2(z-\tfrac12)$, the residue is $\dfrac{(1/2)^2}{-j\cdot2\cdot(\tfrac12-2)}=\dfrac{1/4}{3j}=\dfrac{1}{12j}$. So $J=2\pi j/(12j)=\pi/6$, real, and $I=\pi/6=0.5236$. ✓ (numerically $0.523599$).
>
> Generalization: $\displaystyle\int_0^{2\pi}\frac{\cos n\theta}{1-2a\cos\theta+a^2}d\theta=\frac{2\pi a^n}{1-a^2}$ for $\lvert a\rvert<1$. This is the Fourier-coefficient formula for the Poisson kernel, and in DSP terms it is the statement that $\lvert1-ae^{-j\Omega}\rvert^{-2}$ has inverse DTFT $a^{\lvert n\rvert}/(1-a^2)$.

> [!example] Worked Example 18.28 — a double pole at the origin
> Evaluate $\displaystyle I=\int_0^{2\pi}\frac{\sin^2\theta}{5+4\cos\theta}d\theta$.
>
> Now $\sin^2\theta=-\dfrac{(z-z^{-1})^2}{4}=-\dfrac{(z^2-1)^2}{4z^2}$ and $5+4\cos\theta=\dfrac{(2z+1)(z+2)}{z}$, so
>
> $$
> I=\oint_{\lvert z\rvert=1}f(z)\,dz,\qquad f(z)=-\frac{(z^2-1)^2}{4j\,z^2(2z+1)(z+2)} .
> $$
>
> Inside the unit circle: a double pole at $z=0$ and a simple pole at $z=-\tfrac12$.
>
> *At $z=-\tfrac12$:* $\operatorname{Res}=-\dfrac{(1/4-1)^2}{4j\cdot\tfrac14\cdot2\cdot\tfrac32}=-\dfrac{9/16}{3j}=-\dfrac{3}{16j}$.
>
> *At $z=0$:* let $g(z)=z^2f(z)=-\dfrac{1}{4j}\cdot\dfrac{(z^2-1)^2}{2z^2+5z+2}$. With $N=(z^2-1)^2$ ($N(0)=1$, $N'(0)=0$) and $D=2z^2+5z+2$ ($D(0)=2$, $D'(0)=5$): $(N/D)'(0)=\dfrac{0\cdot2-1\cdot5}{4}=-\dfrac54$, so $\operatorname{Res}=g'(0)=\dfrac{5}{16j}$.
>
> Sum $=\dfrac{2}{16j}=\dfrac1{8j}$, and $I=2\pi j\cdot\dfrac{1}{8j}=\dfrac\pi4=0.7854$. ✓

> [!warning] Common Pitfall — the pole on the wrong side
> Only poles **inside** $\lvert z\rvert=1$ count. For a denominator $a+b\cos\theta$ the two roots are reciprocals of each other (their product is $1$), so exactly one is inside when $a>\lvert b\rvert$. If you find both inside or both outside, there is an algebra error. If $a\le\lvert b\rvert$ a root lies *on* the circle and the real integral diverges.

### Type II — improper integrals of rational functions

To evaluate $\displaystyle\int_{-\infty}^\infty\frac{P(x)}{Q(x)}dx$, close the real segment $[-R,R]$ with the upper semicircle $C_R:\ z=Re^{j\theta}$, $0\le\theta\le\pi$.

```desmos-graph
left=-4; right=4; bottom=-1.5; top=4
---
y=\sqrt{9-x^2}|blue
y=0|blue|-3<x<3
(0,1)|label:pole inside
(0,-1)|label:pole outside
(1.2,2)|label:pole inside
```
*The standard contour for Type II and III integrals: the real segment $[-R,R]$ closed by the upper semicircle $C_R$ (here $R=3$). As $R\to\infty$ it encloses every upper-half-plane pole; lower-half-plane poles never count.*

> [!theorem] Rational integrals over the real line
> Let $P$ and $Q$ be polynomials with $\deg Q\ge\deg P+2$ and $Q(x)\neq0$ for real $x$. Then
>
> $$
> \int_{-\infty}^\infty\frac{P(x)}{Q(x)}dx=2\pi j\sum_{\operatorname{Im}z_k>0}\operatorname{Res}_{z=z_k}\frac{P(z)}{Q(z)} .
> $$

*Proof.* For $R$ large enough, the closed contour $[-R,R]+C_R$ encloses all upper-half-plane poles, so $\int_{-R}^RP/Q\,dx+\int_{C_R}P/Q\,dz=2\pi j\sum\operatorname{Res}$. The degree condition gives $\lvert P/Q\rvert\le K/R^2$ on $C_R$, so by ML $\left\lvert\int_{C_R}\right\rvert\le\pi R\cdot K/R^2\to0$. (The degree condition also guarantees the improper integral converges, so the symmetric limit $\int_{-R}^R$ is the genuine integral.) $\square$

One could equally close in the lower half-plane, traversed clockwise, giving $-2\pi j\sum_{\operatorname{Im}z_k<0}\operatorname{Res}$; for rational functions both give the same answer.

> [!example] Worked Example 18.29 — three rational integrals
> (a) $\displaystyle\int_{-\infty}^\infty\frac{dx}{(x^2+1)^2}$. (b) $\displaystyle\int_{-\infty}^\infty\frac{x^2\,dx}{(x^2+1)(x^2+4)}$. (c) $\displaystyle\int_{-\infty}^\infty\frac{dx}{x^4+1}$.
>
> **(a)** The only upper-half-plane pole is $z=j$, of order $2$, with residue $-j/4$ (Example 18.22(b)). Integral $=2\pi j(-j/4)=\pi/2=1.5708$. ✓
>
> **(b)** Degree condition: $4\ge2+2$. ✓ Upper poles $j$ and $2j$, both simple:
>
> $$
> \operatorname{Res}_{j}=\frac{j^2}{2j\,(j^2+4)}=\frac{-1}{6j},\qquad \operatorname{Res}_{2j}=\frac{(2j)^2}{((2j)^2+1)\,4j}=\frac{-4}{-12j}=\frac{1}{3j}.
> $$
>
> Sum $=\dfrac{1}{6j}$; integral $=2\pi j/(6j)=\pi/3=1.0472$. ✓
>
> **(c)** The roots of $z^4=-1$ are $e^{j\pi/4},e^{j3\pi/4},e^{j5\pi/4},e^{j7\pi/4}$; the upper two are $z_1=e^{j\pi/4}$ and $z_2=e^{j3\pi/4}$. By the quotient rule (Example 18.21(b)) the residues are $-z_k/4$:
>
> $$
> \int=2\pi j\cdot\left(-\frac{z_1+z_2}{4}\right)=-\frac{\pi j}{2}\cdot j\sqrt2=\frac{\pi}{\sqrt2}=2.2214,
> $$
>
> using $z_1+z_2=\tfrac{1+j}{\sqrt2}+\tfrac{-1+j}{\sqrt2}=j\sqrt2$. ✓ (The same method gives $\int dx/(x^6+1)=2\pi/3$.)

### Type III — Fourier-type integrals and Jordan's lemma

Integrals such as $\displaystyle\int_{-\infty}^\infty\frac{\cos\omega x}{x^2+a^2}dx$ are Fourier transforms of rational functions ([Ch. 16](16-fourier-transforms.md)). The trick is to integrate $R(z)e^{j\omega z}$ instead of $R(z)\cos\omega z$: on the upper half-plane, $\lvert e^{j\omega z}\rvert=e^{-\omega y}\le1$ for $\omega>0$, whereas $\cos\omega z$ grows exponentially there. Then take real or imaginary parts.

Here the degree condition can be relaxed by one: the oscillation of $e^{j\omega z}$ is enough to kill the semicircle when $R(z)$ merely tends to zero.

> [!theorem] Jordan's lemma
> Let $\omega>0$, and suppose $\lvert R(z)\rvert\le M_R$ on the upper semicircle $C_R$ with $M_R\to0$ as $R\to\infty$. Then
>
> $$
> \lim_{R\to\infty}\int_{C_R}R(z)\,e^{j\omega z}\,dz=0 .
> $$

*Proof.* On $C_R$, $z=Re^{j\theta}$ and $\lvert e^{j\omega z}\rvert=e^{-\omega R\sin\theta}$, so

$$
\left\lvert\int_{C_R}R\,e^{j\omega z}dz\right\rvert\le M_RR\int_0^\pi e^{-\omega R\sin\theta}d\theta=2M_RR\int_0^{\pi/2}e^{-\omega R\sin\theta}d\theta .
$$

On $[0,\pi/2]$ the concave function $\sin\theta$ lies above its chord: $\sin\theta\ge2\theta/\pi$ (**Jordan's inequality**). Therefore

$$
\int_0^{\pi/2}e^{-\omega R\sin\theta}d\theta\le\int_0^{\pi/2}e^{-2\omega R\theta/\pi}d\theta=\frac{\pi}{2\omega R}\left(1-e^{-\omega R}\right)<\frac{\pi}{2\omega R},
$$

and the whole bound is $<\pi M_R/\omega\to0$. $\square$

The ML inequality alone would give only $\pi RM_R$, which need not vanish if $M_R\sim1/R$; Jordan's lemma recovers the missing factor of $R$ by noting that $e^{-\omega R\sin\theta}$ is tiny except near $\theta=0$ and $\theta=\pi$.

```desmos-graph
left=0; right=1.7; bottom=0; top=1.1
---
y=\sin(x)|blue|0<x<1.5708
y=\frac{2x}{\pi}|red|dashed|0<x<1.5708
```
*Jordan's inequality: on $0\le\theta\le\pi/2$ the sine curve (blue) lies above its chord $2\theta/\pi$ (red, dashed), so $e^{-\omega R\sin\theta}\le e^{-2\omega R\theta/\pi}$ — the bound that makes the semicircle vanish.*

> [!theorem] Fourier-type integrals
> If $R=P/Q$ with $\deg Q\ge\deg P+1$, no real poles, and $\omega>0$, then
>
> $$
> \int_{-\infty}^\infty R(x)\,e^{j\omega x}dx=2\pi j\sum_{\operatorname{Im}z_k>0}\operatorname{Res}_{z=z_k}\Big[R(z)e^{j\omega z}\Big],
> $$
>
> and $\int R\cos\omega x\,dx$, $\int R\sin\omega x\,dx$ are its real and imaginary parts when $R$ is real on the real axis. For $\omega<0$ close in the **lower** half-plane, clockwise: $-2\pi j\sum_{\operatorname{Im}z_k<0}\operatorname{Res}$.

> [!example] Worked Example 18.30 — the Lorentzian and its Fourier transform
> Evaluate $\displaystyle\int_{-\infty}^\infty\frac{\cos\omega x}{x^2+a^2}dx$ for $a>0$ and all real $\omega$; then give the numerical value for $a=1$, $\omega=3$.
>
> For $\omega>0$: the only upper pole of $e^{j\omega z}/(z^2+a^2)$ is $ja$, with residue $e^{-\omega a}/(2ja)$ (Example 18.21(d)). So
>
> $$
> \int_{-\infty}^\infty\frac{e^{j\omega x}}{x^2+a^2}dx=2\pi j\cdot\frac{e^{-a\omega}}{2ja}=\frac{\pi}{a}e^{-a\omega}.
> $$
>
> This is real, so the cosine integral equals it and the sine integral is zero (as oddness requires). The integrand is even in $\omega$, so for all $\omega$:
>
> $$
> \int_{-\infty}^\infty\frac{\cos\omega x}{x^2+a^2}dx=\frac\pi ae^{-a\lvert\omega\rvert}.
> $$
>
> With $a=1$, $\omega=3$: $\pi e^{-3}=0.15641$; numerical quadrature gives $0.156411$. ✓ In transform language this is the pair $\dfrac{1}{t^2+a^2}\leftrightarrow\dfrac{\pi}{a}e^{-a\lvert\omega\rvert}$, the dual of $e^{-a\lvert t\rvert}\leftrightarrow\dfrac{2a}{\omega^2+a^2}$ in [Ch. 16](16-fourier-transforms.md). The kink at $\omega=0$ in the answer is the residue method's way of switching half-planes.

> [!example] Worked Example 18.31 — a case needing Jordan's lemma
> Evaluate $\displaystyle\int_{-\infty}^\infty\frac{x\sin x}{x^2+4}dx$.
>
> Here $R(z)=z/(z^2+4)$ decays only like $1/R$, so ML fails and Jordan's lemma is essential. Consider $\int\dfrac{ze^{jz}}{z^2+4}dz$. Upper pole $2j$, residue $\dfrac{2j\,e^{j(2j)}}{2\cdot2j}=\dfrac{e^{-2}}{2}$. So
>
> $$
> \int_{-\infty}^\infty\frac{xe^{jx}}{x^2+4}dx=2\pi j\cdot\frac{e^{-2}}{2}=j\pi e^{-2}.
> $$
>
> The imaginary part is the required integral: $\pi e^{-2}=0.42517$ (numerically $0.425168$ ✓). The real part, $\int x\cos x/(x^2+4)\,dx=0$, is consistent with odd symmetry. (Strictly, this integral converges only conditionally; the contour argument produces it as $\lim_{R\to\infty}\int_{-R}^R$.)

> [!example] Worked Example 18.32 — double pole, Fourier type
> Evaluate $\displaystyle\int_{-\infty}^\infty\frac{\cos2x}{(x^2+1)^2}dx$.
>
> Take $f=e^{j2z}/(z^2+1)^2$, double pole at $j$. With $g=e^{j2z}/(z+j)^2$:
>
> $$
> g'(z)=e^{j2z}\left[\frac{2j}{(z+j)^2}-\frac{2}{(z+j)^3}\right],\qquad g'(j)=e^{-2}\left[\frac{2j}{-4}-\frac{2}{-8j}\right]=e^{-2}\left[-\frac j2-\frac j4\right]=-\frac{3j}{4}e^{-2}.
> $$
>
> Integral $=2\pi j\cdot(-\tfrac{3j}{4}e^{-2})=\dfrac{3\pi}{2e^2}=0.63775$. ✓ (quadrature $0.637752$).

### Type IV — poles on the real axis: indentation and principal values

If the integrand has a simple pole on the real axis, the real integral does not exist in the ordinary sense, but its **Cauchy principal value** may:

$$
\operatorname{PV}\int_{-\infty}^\infty f(x)\,dx=\lim_{\varepsilon\to0^+}\left[\int_{-\infty}^{x_0-\varepsilon}+\int_{x_0+\varepsilon}^\infty\right]f(x)\,dx .
$$

To compute it, detour around $x_0$ on a small semicircle $\gamma_\varepsilon$ of radius $\varepsilon$ in the upper half-plane.

> [!theorem] Indentation (half-residue) lemma
> If $f$ has a simple pole at $x_0$ and $\gamma_\varepsilon$ is the upper semicircle $z=x_0+\varepsilon e^{j\theta}$ traversed clockwise ($\theta$ from $\pi$ to $0$), then
>
> $$
> \lim_{\varepsilon\to0}\int_{\gamma_\varepsilon}f(z)\,dz=-\pi j\operatorname{Res}_{z=x_0}f .
> $$
>
> More generally an arc of angle $\alpha$ contributes $\alpha j$ times the residue, with sign set by the direction.

*Proof.* Near $x_0$, $f=\dfrac{a_{-1}}{z-x_0}+g(z)$ with $g$ analytic. The $g$ part is bounded and the arc length $\pi\varepsilon\to0$. The principal-part integral is $a_{-1}\int_\pi^0\frac{j\varepsilon e^{j\theta}}{\varepsilon e^{j\theta}}d\theta=-\pi j\,a_{-1}$ exactly. $\square$

Consequently, with the upper semicircle closed as before, a simple real-axis pole counts as **half** a residue:

$$
\operatorname{PV}\int_{-\infty}^\infty f\,dx=2\pi j\sum_{\operatorname{Im}z_k>0}\operatorname{Res}f+\pi j\sum_{\text{real simple poles}}\operatorname{Res}f .
$$

> [!example] Worked Example 18.33 — the Dirichlet integral $\int_0^\infty\frac{\sin x}{x}dx$
> Consider $f(z)=e^{jz}/z$, which has a simple pole at $0$ with residue $1$ and no other singularities. Use the indented contour: $[-R,-\varepsilon]$, the small clockwise semicircle $\gamma_\varepsilon$, $[\varepsilon,R]$, and the large semicircle $C_R$. Nothing is enclosed, so the total is $0$:
>
> $$
> \left[\int_{-R}^{-\varepsilon}+\int_\varepsilon^R\right]\frac{e^{jx}}{x}dx+\int_{\gamma_\varepsilon}f\,dz+\int_{C_R}f\,dz=0 .
> $$
>
> By Jordan's lemma ($M_R=1/R$), $\int_{C_R}\to0$. By the indentation lemma, $\int_{\gamma_\varepsilon}\to-\pi j$. Hence $\operatorname{PV}\int_{-\infty}^\infty\frac{e^{jx}}{x}dx=\pi j$. The real part, $\operatorname{PV}\int\cos x/x\,dx$, is $0$ (odd integrand); the imaginary part is an ordinary convergent integral (the singularity of $\sin x/x$ is removable):
>
> $$
> \int_{-\infty}^\infty\frac{\sin x}{x}dx=\pi,\qquad\int_0^\infty\frac{\sin x}{x}dx=\frac\pi2 .
> $$
>
> This is the area under the unnormalized sinc — the integral behind the Fourier transform of the rectangular pulse, and the limiting value $\operatorname{Si}(\infty)=\pi/2$ of the sine integral that governs the Gibbs overshoot in [Ch. 15](15-fourier-series.md). The answer $\pi$ agrees with numerical quadrature to all printed digits.

```desmos-graph
left=-4; right=4; bottom=-1; top=4
---
y=\sqrt{9-x^2}|blue
y=0|blue|-3<x<-0.5
y=0|blue|0.5<x<3
y=\sqrt{0.25-x^2}|red
(0,0)|label:simple pole
```
*The indented contour for Example 18.33: the small red semicircle detours above the pole at the origin, so the pole lies outside the closed contour, and in the limit the detour contributes $-\pi j$ times the residue.*

**Branch cuts (keyhole contours), briefly.** Integrals such as $\int_0^\infty\frac{x^{a-1}}{1+x}dx$ ($0<a<1$) use a multivalued integrand $z^{a-1}$ with its cut along the positive real axis. The *keyhole* contour runs out along the top of the cut, around a large circle, back along the bottom of the cut (where $z^{a-1}$ has gained a factor $e^{j2\pi(a-1)}$) and around a small circle at the origin. The two cut integrals combine to $(1-e^{j2\pi a})I$; the residue at $z=-1=e^{j\pi}$ is $e^{j\pi(a-1)}$; solving gives

$$
\int_0^\infty\frac{x^{a-1}}{1+x}dx=\frac{\pi}{\sin\pi a}.
$$

For $a=\tfrac12$ this is $\pi$; for $a=\tfrac23$ it is $2\pi/\sqrt3=3.6276$, as numerical quadrature confirms.

> [!warning] Common Pitfall — closing on the wrong side
> For $e^{j\omega z}$ with $\omega>0$, close in the **upper** half-plane, where $\lvert e^{j\omega z}\rvert=e^{-\omega y}$ decays. Closing in the lower half-plane makes the arc integral blow up and gives a wrong (often complex) answer. In inverse-transform applications this is exactly the choice that produces causality, as the next section shows. And never integrate $\cos\omega z$ directly over a semicircle: it contains both $e^{j\omega z}$ and $e^{-j\omega z}$, one of which grows on each side.

> [!tip] Engineering Note — why these integrals matter
> Type I integrals are averages over one period of a periodic quantity (mean power in a phase-modulated signal, the response of a ring resonator over a revolution). Type II integrals give the total noise power at a filter output: for white noise into $H(s)$, the output variance is proportional to $\frac{1}{2\pi}\int\lvert H(j\omega)\rvert^2d\omega$, a rational integral solved by residues (the basis of tabulated "$I_n$" integrals in control texts). Type III integrals are inverse Fourier transforms; Type IV, with its half-residues, reappears as the $\pi\delta(\omega)$ terms in the transforms of the step and the signum.

## 18.9 Applications to engineering analysis

This section collects the places where the machinery of Sections 18.5–18.8 underpins the transform methods of Part V. The theme is the same throughout: an engineering quantity is a contour integral, the contour can be closed or deformed, and the answer is a sum of residues or a count of zeros and poles.

### The Bromwich inversion integral for the Laplace transform

[Ch. 14](14-laplace-transforms.md) inverted Laplace transforms with tables and partial fractions. The general inversion formula is a contour integral.

> [!theorem] Bromwich (complex inversion) integral
> Let $F(s)=\mathcal{L}\{f(t)\}$ converge for $\operatorname{Re}s>\sigma_0$, and let $\gamma>\sigma_0$. Then, at points where $f$ is continuous,
>
> $$
> f(t)=\frac{1}{2\pi j}\int_{\gamma-j\infty}^{\gamma+j\infty}F(s)\,e^{st}\,ds=\frac{1}{2\pi j}\lim_{W\to\infty}\int_{\gamma-jW}^{\gamma+jW}F(s)e^{st}ds ,
> $$
>
> the path being the vertical line $\operatorname{Re}s=\gamma$ (the **Bromwich contour**), which lies to the right of every singularity of $F$. At a jump of $f$ the integral gives the midpoint value.

*Derivation from Fourier inversion.* For fixed $\gamma$ in the region of convergence, $g(t)=f(t)e^{-\gamma t}u(t)$ is absolutely integrable, and its Fourier transform is $G(\omega)=\int_0^\infty f(t)e^{-(\gamma+j\omega)t}dt=F(\gamma+j\omega)$. The Fourier inversion formula of [Ch. 16](16-fourier-transforms.md) gives

$$
f(t)e^{-\gamma t}=\frac{1}{2\pi}\int_{-\infty}^\infty F(\gamma+j\omega)e^{j\omega t}d\omega .
$$

Multiply by $e^{\gamma t}$ and substitute $s=\gamma+j\omega$, $ds=j\,d\omega$: the result is the Bromwich integral. $\square$ So the Laplace inversion formula is nothing more than Fourier inversion along a line shifted into the region of convergence. The value of $\gamma$ does not matter, as long as the line stays right of all singularities — by Cauchy's theorem the integrand is analytic between any two such lines.

**Evaluation by residues.** Suppose $F$ is rational and strictly proper (or, more generally, $\lvert F(s)\rvert\to0$ uniformly as $\lvert s\rvert\to\infty$ in the left region). Close the line with a large arc $\Gamma_R$ to the left.

- For $t>0$: $\lvert e^{st}\rvert=e^{t\operatorname{Re}s}$ is small where $\operatorname{Re}s$ is very negative. A rotated form of Jordan's lemma (substitute $s=jz$; the left half-plane becomes the upper half-plane and $e^{st}=e^{jtz}$) shows $\int_{\Gamma_R}F(s)e^{st}ds\to0$. The closed contour encloses **all** the poles of $F$.
- For $t<0$: close to the **right** instead, where $e^{st}$ now decays. No singularities lie right of the Bromwich line, so $f(t)=0$ — **causality emerges from analyticity in the right half-plane**.

> [!theorem] Inversion by residues
> If $F(s)$ is rational and strictly proper, with poles $p_1,\dots,p_n$, then
>
> $$
> f(t)=\sum_{k=1}^n\operatorname{Res}_{s=p_k}\Big[F(s)e^{st}\Big],\qquad t>0,
> $$
>
> and $f(t)=0$ for $t<0$. At a simple pole $p$ the term is $\big[(s-p)F(s)\big]_{s=p}\,e^{pt}$; at a pole of order $m$ it is $\dfrac{1}{(m-1)!}\dfrac{d^{m-1}}{ds^{m-1}}\Big[(s-p)^mF(s)e^{st}\Big]_{s=p}$, which produces the $t^ke^{pt}$ terms automatically.

This one formula replaces the table look-up and the partial-fraction step together: the derivative of $e^{st}$ with respect to $s$ is $te^{st}$, which is where the polynomial-times-exponential terms of repeated poles come from.

```desmos-graph
left=-4; right=3; bottom=-4; top=4
---
x=1|blue|-3.6<y<3.6
(x-1)^2+y^2=12.96|blue|x<1
(-1,0)|label:pole
(-2,1.5)|label:pole
(-2,-1.5)|label:pole
```
*The Bromwich contour: the vertical line $\operatorname{Re}s=\gamma$ (here $\gamma=1$) to the right of all poles, closed for $t>0$ by a large arc to the left, which by Jordan's lemma contributes nothing as its radius grows; the enclosed residues of $F(s)e^{st}$ give $f(t)$.*

> [!example] Worked Example 18.34 — simple and repeated real poles
> Invert (a) $F(s)=\dfrac{1}{s(s+1)^2}$ and (b) $F(s)=\dfrac{s+3}{(s+1)(s+2)^2}$.
>
> **(a)** At $s=0$ (simple): $\left[\dfrac{e^{st}}{(s+1)^2}\right]_{s=0}=1$. At $s=-1$ (order 2):
>
> $$
> \frac{d}{ds}\left[\frac{e^{st}}{s}\right]_{s=-1}=\left[\frac{te^{st}}{s}-\frac{e^{st}}{s^2}\right]_{s=-1}=-te^{-t}-e^{-t}.
> $$
>
> So $f(t)=1-e^{-t}-te^{-t}$ for $t\ge0$. Checks: $f(0)=0$ (consistent with the initial-value theorem, $\lim sF=0$); $f(\infty)=1=\lim_{s\to0}sF$; at $t=2$, $f=1-3e^{-2}=0.59399$, matching a numerical Bromwich inversion. ✓
>
> **(b)** At $s=-1$: $\left[\dfrac{(s+3)e^{st}}{(s+2)^2}\right]_{-1}=2e^{-t}$. At $s=-2$ (order 2), with $G(s)=\dfrac{(s+3)e^{st}}{s+1}$:
>
> $$
> G'(s)=\frac{\big[e^{st}+t(s+3)e^{st}\big](s+1)-(s+3)e^{st}}{(s+1)^2},\qquad G'(-2)=\frac{(1+t)e^{-2t}(-1)-e^{-2t}}{1}=-(2+t)e^{-2t}.
> $$
>
> So $f(t)=2e^{-t}-2e^{-2t}-te^{-2t}$ — the time-domain image of the partial fractions of Example 18.23. At $t=1.3$: $0.29996$, matching numerical inversion. ✓

> [!example] Worked Example 18.35 — a complex-conjugate pair
> Invert $F(s)=\dfrac{1}{s^2+2s+5}$.
>
> Poles $p,\bar p=-1\pm2j$. At $p=-1+2j$: $\dfrac{e^{pt}}{p-\bar p}=\dfrac{e^{(-1+2j)t}}{4j}$. The residue at $\bar p$ is the conjugate (because $F$ is real on the real axis), so
>
> $$
> f(t)=2\operatorname{Re}\left\{\frac{e^{(-1+2j)t}}{4j}\right\}=\frac{e^{-t}}{2}\operatorname{Re}\{-je^{j2t}\}=\frac12e^{-t}\sin2t .
> $$
>
> In general a conjugate pair contributes $2\operatorname{Re}\{\operatorname{Res}_p\}$ — a damped sinusoid whose amplitude and phase are the modulus and angle of the residue. At $t=1.3$: $0.070245$. ✓

> [!warning] Common Pitfall — when the arc does not vanish
> The residue-sum formula needs $F(s)\to0$ on the left arc. It fails for an improper $F$ (e.g. $F=s/(s+1)=1-1/(s+1)$, whose constant part is the impulse $\delta(t)$ and must be split off first), and it must be adapted for a delay $e^{-sT}F(s)$: then $e^{s(t-T)}$ decays on the left arc only for $t>T$, so the residues give $f(t-T)$ for $t>T$ and closing right gives $0$ for $t<T$ — the shift theorem recovered. Functions with branch points, such as $F=1/\sqrt s$ (inverse $1/\sqrt{\pi t}$), need a contour that wraps around the branch cut on the negative real axis rather than a residue sum.

### Inverse Z-transform by residues

The same idea inverts Z-transforms ([Ch. 17](17-z-transforms.md)). Since $X(z)=\sum_nx[n]z^{-n}$ is a Laurent series converging in the ROC annulus, the Laurent coefficient formula gives

$$
x[n]=\frac{1}{2\pi j}\oint_CX(z)\,z^{n-1}\,dz=\sum_{\text{poles inside }C}\operatorname{Res}\Big[X(z)z^{n-1}\Big],
$$

where $C$ is a counter-clockwise circle in the ROC. For a causal sequence (ROC outside the outermost pole) $C$ encloses all poles of $X$ — but watch for an extra pole of $z^{n-1}$ at $z=0$ when $n$ is small.

> [!example] Worked Example 18.36 — causal inverse Z-transform
> Find $x[n]$ for $X(z)=\dfrac{z}{(z-0.5)(z-1)}$, $\lvert z\rvert>1$.
>
> $X(z)z^{n-1}=\dfrac{z^n}{(z-0.5)(z-1)}$, which has no pole at $z=0$ for $n\ge0$. Residues: at $z=1$, $\dfrac{1}{0.5}=2$; at $z=0.5$, $\dfrac{0.5^n}{-0.5}=-2(0.5)^n$. So
>
> $$
> x[n]=\big(2-2(0.5)^n\big)u[n]=0,\ 1,\ 1.5,\ 1.75,\ 1.875,\dots
> $$
>
> Check: $\sum_nx[n]z^{-n}$ at $z=2$ sums to $1.3333$, equal to $X(2)=2/(1.5\cdot1)$. ✓ For $n<0$ the integrand acquires a pole at $z=0$ of order $-n$ whose residue cancels the others exactly, giving $x[n]=0$ — causality again, provable more easily via the residue at infinity.

### The argument principle

The final group of applications counts zeros and poles without finding them.

> [!theorem] Argument principle
> Let $f$ be analytic inside and on a positively oriented simple closed contour $C$ except for finitely many poles inside $C$, and let $f$ have no zeros or poles on $C$. Then
>
> $$
> \frac{1}{2\pi j}\oint_C\frac{f'(z)}{f(z)}dz=Z-P ,
> $$
>
> where $Z$ and $P$ are the numbers of zeros and poles of $f$ inside $C$, each counted with multiplicity. Equivalently, $Z-P=\dfrac{1}{2\pi}\Delta_C\arg f(z)$: the net number of times the image curve $f(C)$ winds counter-clockwise around the origin.

*Proof.* The logarithmic derivative $f'/f$ is analytic except at zeros and poles of $f$. Near a zero of order $m$, $f=(z-z_0)^mg$ with $g(z_0)\ne0$, so $f'/f=\dfrac{m}{z-z_0}+\dfrac{g'}{g}$: a simple pole with residue $m$. Near a pole of order $m$, $f=(z-z_0)^{-m}g$, and the residue of $f'/f$ is $-m$. The residue theorem gives $Z-P$. For the second form, note $f'/f=\frac{d}{dz}\log f$, so $\oint_Cf'/f\,dz$ is the change of $\log f=\ln\lvert f\rvert+j\arg f$ around $C$; $\ln\lvert f\rvert$ returns to its starting value and $\arg f$ changes by $2\pi(Z-P)$. $\square$

The second form is the useful one: **draw the image of $C$ under $f$ and count how many times it encircles the origin.** That is a graphical computation — and exactly what a Nyquist plot is.

> [!theorem] Rouché's theorem
> If $f$ and $g$ are analytic inside and on $C$ and $\lvert g(z)\rvert<\lvert f(z)\rvert$ on $C$, then $f$ and $f+g$ have the same number of zeros inside $C$.

(Intuition: on $C$, $f+g=f(1+g/f)$ with $\lvert g/f\rvert<1$, so $1+g/f$ stays in the disc $\lvert w-1\rvert<1$ and cannot wind around $0$; the winding of $f+g$ equals that of $f$. This is the "dog on a leash" argument.)

> [!example] Worked Example 18.37 — locating polynomial roots
> How many roots of $p(z)=z^5+3z+1$ lie in $\lvert z\rvert<1$, and how many in $1<\lvert z\rvert<2$?
>
> *On $\lvert z\rvert=1$:* take $f=3z$, $g=z^5+1$. Then $\lvert g\rvert\le2<3=\lvert f\rvert$, so $p$ has as many zeros in the unit disc as $3z$: **one**.
>
> *On $\lvert z\rvert=2$:* take $f=z^5$, $g=3z+1$. Then $\lvert g\rvert\le7<32=\lvert f\rvert$, so $p$ has five zeros in $\lvert z\rvert<2$.
>
> Hence four roots lie in the annulus $1<\lvert z\rvert<2$. Numerically the roots have moduli $0.332$, $1.263$ (twice) and $1.374$ (twice). ✓ The integral $\frac{1}{2\pi j}\oint_{\lvert z\rvert=1}p'/p\,dz$ evaluated by quadrature gives $1.0000$. ✓ For a discrete-time characteristic polynomial, the same count on the unit circle decides stability: all roots must lie inside.

### The Nyquist stability criterion — the idea

A feedback loop with loop transfer function $L(s)$ (forward path times feedback path) has closed-loop transfer function $T(s)=\dfrac{G(s)}{1+L(s)}$. Its poles are the **zeros of $1+L(s)$**. The closed loop is stable when $1+L$ has no zeros in the closed right half-plane. Computing those zeros means factoring a polynomial; Nyquist's insight is to *count* them with the argument principle using only the frequency response $L(j\omega)$, which can be measured.

1. **The contour.** Take the *Nyquist contour* $D$: up the imaginary axis from $-jR$ to $jR$, then around the large semicircle $Re^{j\theta}$ from $\theta=\pi/2$ to $-\pi/2$. It encloses the right half-plane and is traversed **clockwise**. (If $L$ has poles on the $j\omega$-axis, indent around them with small semicircles into the left half-plane.)
2. **Apply the argument principle to $F(s)=1+L(s)$.** Inside $D$, let $Z$ = number of zeros of $1+L$ (closed-loop RHP poles, unknown) and $P$ = number of poles of $1+L$, which are the poles of $L$ (open-loop RHP poles, known from the model). Because $D$ is clockwise, the image $F(D)$ encircles the origin $N=Z-P$ times **clockwise**.
3. **Shift the origin.** $F(D)$ encircling $0$ is the same as $L(D)=F(D)-1$ encircling $-1$. For strictly proper $L$, the large semicircle maps to the single point $L(\infty)=0$, so $L(D)$ is just the Nyquist plot $L(j\omega)$, $-\infty<\omega<\infty$.

> [!theorem] Nyquist criterion
> Let $N$ be the number of clockwise encirclements of the point $-1$ by the Nyquist plot $L(j\omega)$, $-\infty<\omega<\infty$, and $P$ the number of open-loop poles of $L$ in the right half-plane. The number of closed-loop poles in the right half-plane is
>
> $$
> Z=N+P .
> $$
>
> The closed loop is stable if and only if $Z=0$, i.e. the plot encircles $-1$ exactly $P$ times **counter-clockwise**. If $L$ is open-loop stable ($P=0$), stability requires no encirclement of $-1$.

> [!example] Worked Example 18.38 — critical gain from the Nyquist plot
> Let $L(s)=\dfrac{K}{(s+1)^3}$ with $K>0$. Find the range of $K$ for closed-loop stability.
>
> All open-loop poles are at $-1$, so $P=0$ and we need no encirclement of $-1$. On the imaginary axis,
>
> $$
> L(j\omega)=\frac{K}{(1+j\omega)^3}=\frac{K}{(1+\omega^2)^{3/2}}e^{-j3\arctan\omega}.
> $$
>
> At $\omega=0$ the plot starts at $K$ on the positive real axis; as $\omega\to\infty$ the phase tends to $-3\pi/2$ and the magnitude to $0$. It crosses the negative real axis when the phase is $-\pi$: $3\arctan\omega=\pi$, $\omega_{180}=\tan(\pi/3)=\sqrt3$. There $\lvert L\rvert=K/(1+3)^{3/2}=K/8$, so the crossing point is $-K/8$.
>
> - If $K<8$, the crossing lies between $-1$ and $0$: no encirclement, $Z=0$, **stable**. The gain margin is $8/K$.
> - If $K>8$, the crossing lies left of $-1$; the plot (with its mirror image for $\omega<0$) encircles $-1$ twice clockwise, so $Z=2$: two closed-loop poles in the right half-plane.
> - At $K=8$ the plot passes through $-1$: closed-loop poles on the $j\omega$-axis at $\pm j\sqrt3$.
>
> Check by algebra: the closed-loop characteristic polynomial is $(s+1)^3+K=s^3+3s^2+3s+1+K$, and the Routh condition $3\cdot3>1+K$ gives $K<8$. With $K=8$ its roots are $-3$ and $\pm j\sqrt3=\pm1.732j$. ✓

Writing $L(j\omega)$ in polar form, $r=\lvert L\rvert$ and $\theta=\angle L=-3\arctan\omega$, we have $1/\sqrt{1+\omega^2}=\cos(\theta/3)$, so the whole Nyquist plot is the polar curve $r=K\cos^3(\theta/3)$, $\lvert\theta\rvert\le3\pi/2$.

```desmos-graph
left=-2.5; right=4.5; bottom=-3; top=3
---
r=2\left(\cos(\frac{\theta}{3})\right)^3|blue
r=8\left(\cos(\frac{\theta}{3})\right)^3|red|dashed
(-1,0)|label:critical point
```
*Nyquist plots of $L(s)=K/(s+1)^3$: for $K=2$ (blue) the curve crosses the negative real axis at $-0.25$ and does not encircle $-1$ (stable); for $K=8$ (red, dashed; only the part near the origin is visible) the crossing is exactly at $-1$, the edge of instability. Larger $K$ scales the curve outward until it encircles $-1$.*

> [!tip] Engineering Note — why Nyquist beats root-finding
> The criterion needs only $L(j\omega)$, which can be **measured** with a network analyser on real hardware, including pure delays $e^{-sT}$ (which make $1+L$ non-polynomial and defeat Routh's test). Gain and phase margins are read directly as distances of the plot from $-1$. The mathematics underneath is entirely the argument principle: closed-loop poles are counted by winding numbers.

> [!warning] Common Pitfall — orientation and the open-loop count
> Two sign conventions collide here: the Nyquist contour is traversed clockwise, and textbooks count encirclements clockwise as positive, giving $Z=N+P$. If you count counter-clockwise encirclements as positive the formula becomes $Z=P-N_{\text{ccw}}$. Also, $P$ counts **open-loop RHP poles**, not zeros; forgetting an unstable open-loop pole (e.g. an inverted pendulum or a magnetic levitator) leads to the wrong conclusion even when the plot is drawn correctly.

### Gain–phase relations (Kramers–Kronig), briefly

Let $H(\omega)$ be the frequency response of a causal, stable system, considered as a function of complex $\omega$; causality makes it analytic in the lower half $\omega$-plane for the $e^{j\omega t}$ convention ($s=j\omega$ in the right half $s$-plane), and suppose $H\to0$ as $\lvert\omega\rvert\to\infty$. Integrate $H(\zeta)/(\zeta-\omega)$ around the real axis indented at $\zeta=\omega$ and closed in the half-plane of analyticity. There are no enclosed singularities, the large arc vanishes, and the indentation contributes $\pi j$ times the residue $H(\omega)$ — exactly the Type IV computation of Example 18.33. The result is

$$
H(\omega)=\frac{j}{\pi}\operatorname{PV}\int_{-\infty}^\infty\frac{H(\zeta)}{\zeta-\omega}d\zeta\quad\Longrightarrow\quad \operatorname{Re}H(\omega)=-\frac1\pi\operatorname{PV}\int_{-\infty}^\infty\frac{\operatorname{Im}H(\zeta)}{\zeta-\omega}d\zeta .
$$

(The closed contour runs clockwise, and the small indentation passes below $\zeta=\omega$, so it turns counter-clockwise about the pole and contributes $+\pi jH(\omega)$; the PV integral plus this term is zero.) Check with $H=1/(1+j\omega)$: $\operatorname{Im}H=-\zeta/(1+\zeta^2)$, and the right side at $\omega=0$ is $\frac1\pi\int_{-\infty}^{\infty}\frac{d\zeta}{1+\zeta^2}=1=\operatorname{Re}H(0)$. ✓ The real and imaginary parts of a causal frequency response are **Hilbert transforms** of each other: specifying one determines the other. Applied to $\log H=\ln\lvert H\rvert+j\angle H$ for minimum-phase systems (no RHP zeros, so $\log H$ is analytic too), this is Bode's gain–phase relation, the reason a slope of $-20$ dB/decade goes with about $-90^\circ$ of phase.

## Summary

| Topic | Key formula or statement |
|---|---|
| Complex function | $f(z)=u(x,y)+jv(x,y)$ |
| Exponential, logarithm | $e^z=e^x(\cos y+j\sin y)$; $\operatorname{Log}z=\ln\lvert z\rvert+j\operatorname{Arg}z$, $\operatorname{Arg}\in(-\pi,\pi]$ |
| Complex power | $z^\alpha=e^{\alpha\log z}$ (principal value uses $\operatorname{Log}$) |
| Derivative | $f'(z_0)=\lim_{\Delta z\to0}[f(z_0+\Delta z)-f(z_0)]/\Delta z$, direction-independent |
| Cauchy–Riemann | $u_x=v_y,\ u_y=-v_x$; $f'=u_x+jv_x$ |
| Polar C–R | $u_r=v_\theta/r,\ v_r=-u_\theta/r$ |
| Harmonic parts | $\nabla^2u=\nabla^2v=0$; level curves of $u$, $v$ orthogonal |
| Conformality | angles preserved where $f'\ne0$; local scale $\lvert f'\rvert$, area scale $\lvert f'\rvert^2$ |
| Contour integral | $\int_Cf\,dz=\int_a^bf(z(t))z'(t)\,dt$ |
| Fundamental integral | $\oint_{\lvert z-z_0\rvert=\rho}(z-z_0)^ndz=2\pi j$ if $n=-1$, else $0$ |
| ML inequality | $\left\lvert\int_Cf\,dz\right\rvert\le ML$ |
| Cauchy's theorem | $\oint_Cf\,dz=0$ for $f$ analytic inside and on $C$ |
| Cauchy integral formula | $f^{(n)}(z_0)=\dfrac{n!}{2\pi j}\oint_C\dfrac{f(z)}{(z-z_0)^{n+1}}dz$ |
| Cauchy estimate, Liouville | $\lvert f^{(n)}(z_0)\rvert\le n!M/\rho^n$; bounded entire $\Rightarrow$ constant |
| Taylor series | $a_n=f^{(n)}(z_0)/n!$; radius = distance to nearest singularity |
| Laurent series | $f=\sum_{n=-\infty}^\infty a_n(z-z_0)^n$ in $r_1<\lvert z-z_0\rvert<r_2$; $a_n=\frac{1}{2\pi j}\oint\frac{f}{(z-z_0)^{n+1}}dz$ |
| Singularities | removable / pole of order $m$ / essential, by the principal part |
| Residue, simple pole | $\lim_{z\to z_0}(z-z_0)f(z)=g(z_0)/h'(z_0)$ |
| Residue, order $m$ | $\frac{1}{(m-1)!}\lim\frac{d^{m-1}}{dz^{m-1}}[(z-z_0)^mf]$ |
| Residue theorem | $\oint_Cf\,dz=2\pi j\sum_{\text{inside}}\operatorname{Res}f$ |
| Trig integrals | $z=e^{j\theta}$, $\cos\theta=\frac{z+z^{-1}}{2}$, $\sin\theta=\frac{z-z^{-1}}{2j}$, $d\theta=\frac{dz}{jz}$; residues in $\lvert z\rvert<1$ |
| Rational integrals | $\int_{-\infty}^\infty\frac PQdx=2\pi j\sum_{\operatorname{Im}z_k>0}\operatorname{Res}$, $\deg Q\ge\deg P+2$ |
| Jordan's lemma | $\int_{C_R}R(z)e^{j\omega z}dz\to0$ ($\omega>0$, $R\to0$ on $C_R$) |
| Fourier-type integrals | $\int_{-\infty}^\infty R(x)e^{j\omega x}dx=2\pi j\sum_{\operatorname{Im}z_k>0}\operatorname{Res}[Re^{j\omega z}]$, $\omega>0$ |
| Indentation | simple real pole contributes $\pi j\operatorname{Res}$ to the PV |
| Bromwich integral | $f(t)=\frac{1}{2\pi j}\int_{\gamma-j\infty}^{\gamma+j\infty}F(s)e^{st}ds=\sum\operatorname{Res}[F(s)e^{st}]$, $t>0$ |
| Inverse Z-transform | $x[n]=\frac{1}{2\pi j}\oint_CX(z)z^{n-1}dz=\sum_{\text{inside}}\operatorname{Res}[X(z)z^{n-1}]$ |
| Argument principle | $\frac{1}{2\pi j}\oint_C\frac{f'}{f}dz=Z-P=\frac{1}{2\pi}\Delta_C\arg f$ |
| Nyquist criterion | $Z=N+P$ (clockwise encirclements of $-1$); stable iff $Z=0$ |

## Further reading

- E. Kreyszig, *Advanced Engineering Mathematics*, 10th ed., Wiley, 2011 — Part D (complex analysis), with conformal mapping and potential theory.
- J. W. Brown and R. V. Churchill, *Complex Variables and Applications*, 9th ed., McGraw-Hill, 2014 — the standard undergraduate text; thorough on branches, residues and real integrals.
- E. B. Saff and A. D. Snider, *Fundamentals of Complex Analysis with Applications to Engineering, Science, and Mathematics*, 3rd ed., Pearson, 2003 — engineering-oriented, including the argument principle and transform inversion.
- T. Needham, *Visual Complex Analysis*, Oxford University Press, 1997 — geometric intuition for mappings, conformality and winding numbers.
