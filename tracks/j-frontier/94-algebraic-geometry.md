# Chapter 94 — Algebraic Geometry: Solving Equations Through Shapes

[← Previous](93-tensor-mathematics.md) · [Track J](README.md) · [Table of Contents](../../00-front/table-of-contents.md) · [Next →](95-representation-theory.md)

> **Level 3 · Track J — Mathematics at the Frontier**
> **Builds on:** [Ch. 1](../../part-1-foundations/01-mathematical-language-and-algebra.md) (§1.4 polynomials), [Ch. 56](../e-abstract-math/56-abstract-algebra.md) (rings and fields), [Ch. 18](../../part-5-transforms/18-complex-analysis.md) · **Leads to:** [Ch. 95](95-representation-theory.md), [Ch. 78](../h-cryptography-coding/78-elliptic-curves.md), [Ch. 98](98-mathematical-research.md) · **Related:** [Ch. 57](../e-abstract-math/57-number-theory.md), [Ch. 75](../g-graphics-simulation/75-computational-geometry.md), [Ch. 82](../h-cryptography-coding/82-zero-knowledge-proofs.md), [Ch. 93](93-tensor-mathematics.md)
> **Reading time:** ~30 min

## 94.1 The hook

A planar robot arm has its shoulder bolted to the origin, an upper arm of length 2 and a forearm of length 1. You want its hand at the point $(2,1)$. Where must the elbow go?

The elbow $(u,v)$ is at distance 2 from the shoulder and distance 1 from the hand:

$$
u^2+v^2=4,\qquad (u-2)^2+(v-1)^2=1 .
$$

Subtract the two equations and the squares cancel, leaving the line $2u+v=4$. Substitute $v=4-2u$ back into the first equation to get $5u^2-16u+12=0$, so $u=2$ or $u=1.2$. The two elbows are $(2,0)$ ("elbow down") and $(1.2,1.6)$ ("elbow up"). Every robotics engineer knows this answer.

Now look at it as a mathematician would. Each equation has degree 2: each describes a **conic**. There is a famous theorem, stated by Newton and Bézout and proved properly much later, that two curves of degrees $d$ and $e$ meet in **exactly $de$ points**. For two conics that means four. An ellipse and a hyperbola really do meet in four points. But two different circles *never* meet in more than two. You can slide them about as you like and you will never find a third intersection.

So either Bézout's theorem is false, or two intersection points of every pair of circles are hiding somewhere.

> **Where are the two missing intersection points of two circles?**

Answering it means changing what "point" means. Allow complex coordinates and add points at infinity, and polynomial equations stop behaving erratically and start obeying exact counting laws. That change of viewpoint, from equations to the *shapes* they define and back, is algebraic geometry. It now underlies robot kinematics, computer-aided design, elliptic-curve cryptography and the proof of Fermat's Last Theorem.

## 94.2 The big idea

**Idea 1: an equation is a shape.** The set of points where a polynomial vanishes is a geometric object: $x^2+y^2-1=0$ is a circle, $y^2-x^3=0$ is a curve with a sharp point. Such zero sets are called **varieties**. Their geometry, such as how many pieces they have, where they are pointed or self-crossing, and how they meet, tells you about the equations. In return, algebra about the polynomials, such as factoring, substituting and eliminating, tells you about the shapes.

![Four curves: a circle; the cusp y squared equals x cubed with a sharp point at the origin; the nodal cubic y squared equals x squared times (x plus 1), which crosses itself at the origin with tangent lines y equals x and y equals minus x; and the smooth cubic y squared equals x cubed minus x, with an oval and an unbounded branch](../../images/ch94-curve-gallery.svg)

*Four plane curves. The cusp and the node have singular points, circled in red, where the curve fails to look like a smooth line. The last curve is an elliptic curve. Its real points form two pieces, but its complex points, together with one point at infinity, form a single torus.*

**Idea 2: the shape doesn't care which equations you wrote.** The circle is described by $x^2+y^2-1=0$, but equally by $(x^2+y^2-1)(x^2+1)=0$, or by any system of equations whose consequences are the same. What actually determines the shape is the set of *all* polynomials that follow from the equations by adding and multiplying. That set is an **ideal**, and ideals are the algebraic shadows of shapes.

**Idea 3: functions on the shape form a ring, and the ring knows everything.** On the circle, $x^2$ and $1-y^2$ are the same function, because they agree at every point. Identify polynomials that agree on the shape and you get the **coordinate ring**: the algebra of polynomial functions on the variety. Two shapes are "the same" exactly when their coordinate rings are the same ring. Geometry turns into ring theory, and the ring theory of [Ch. 56](../e-abstract-math/56-abstract-algebra.md#rings-and-fields-two-operations-at-once) gets a picture.

An analogy: a variety is a building, and its coordinate ring is the complete set of measurements you could take inside it. Two buildings that give identical answers to every possible measurement are, for every practical purpose, the same building.

## 94.3 The mathematics

Fix a field $k$, usually $\mathbb{R}$ or $\mathbb{C}$, and the polynomial ring $k[x_1,\dots,x_n]$.

> [!definition] Affine variety and ideal
> For polynomials $f_1,\dots,f_m$, the **affine variety** they define is
> $$
> V(f_1,\dots,f_m)=\{\mathbf{a}\in k^n: f_1(\mathbf{a})=\dots=f_m(\mathbf{a})=0\}.
> $$
> The **ideal** they generate is $\langle f_1,\dots,f_m\rangle=\{h_1f_1+\dots+h_mf_m : h_i\in k[x_1,\dots,x_n]\}$.

Every element of the ideal vanishes wherever all the $f_i$ do. So $V(f_1,\dots,f_m)$ depends only on the ideal, and any other generators of the same ideal describe the same set. The robot arm's trick of subtracting equations was a change of generators: $\langle f_1,f_2\rangle=\langle f_1,\,f_1-f_2\rangle$. Two set operations have one-line algebraic counterparts: $V(f)\cup V(g)=V(fg)$ (a product vanishes when a factor does) and $V(f)\cap V(g)=V(f,g)$.

Going the other way, from a set $X\subseteq k^n$ to polynomials, gives the ideal $I(X)$ of all polynomials vanishing on $X$. Two classical theorems make the dictionary work.

- **Hilbert's basis theorem (1890).** Every ideal of $k[x_1,\dots,x_n]$ is generated by finitely many polynomials. So every variety, even one described by infinitely many equations, is cut out by finitely many.
- **Hilbert's Nullstellensatz (1893).** Over $\mathbb{C}$, or any algebraically closed field: $V(J)=\varnothing$ only if $1\in J$, and in general $I(V(J))$ is the **radical** $\sqrt J=\{f: f^r\in J\text{ for some }r\}$.

The Nullstellensatz fails over $\mathbb{R}$: $x^2+y^2+1$ has no real zeros, yet $1\notin\langle x^2+y^2+1\rangle$. This is the first sign that the reals are the wrong place to count solutions, just as the real numbers were the wrong place to count roots of $x^2+1$ ([Ch. 1](../../part-1-foundations/01-mathematical-language-and-algebra.md#14-polynomials-and-rational-expressions)).

### Smooth and singular points

A point $\mathbf{a}$ of the plane curve $V(f)$ is **singular** if $\nabla f(\mathbf{a})=\mathbf{0}$. At a non-singular point the implicit function theorem ([Ch. 6](../../part-2-calculus/06-multivariable-calculus.md)) makes the curve locally the graph of a smooth function. At a singular point it can fail.

- **Cusp**, $f=y^2-x^3$. Here $\nabla f=(-3x^2,2y)$, which vanishes only at the origin, and the origin lies on the curve. Near it, the lowest-order part of $f$ is $y^2$: one tangent line $y=0$, counted twice. The curve arrives along it from both sides and turns back.
- **Node**, $f=y^2-x^2(x+1)=y^2-x^2-x^3$. Here $\nabla f=(-2x-3x^2,\,2y)$, which vanishes at $(0,0)$ and at $(-\tfrac23,0)$. Only the origin lies on the curve. The lowest-order part, $y^2-x^2=(y-x)(y+x)$, gives two distinct tangent lines $y=\pm x$: two branches cross.

For the family $y^2=x^3+ax+b$, the curve is singular exactly when the cubic $x^3+ax+b$ has a repeated root, that is, when $4a^3+27b^2=0$. The cusp is $a=b=0$. For $a=-3$, $b=2$ the discriminant is $4(-3)^3+27\cdot2^2=-108+108=0$, and indeed $x^3-3x+2=(x-1)^2(x+2)$: the curve $y^2=x^3-3x+2$ has a node at $(1,0)$, just like the one in the figure. When the discriminant is nonzero the curve is smooth, and it is an **elliptic curve**, the object behind [Ch. 78](../h-cryptography-coding/78-elliptic-curves.md).

### The coordinate ring

> [!definition] Coordinate ring
> The **coordinate ring** of a variety $X\subseteq k^n$ is $k[X]=k[x_1,\dots,x_n]/I(X)$: polynomials, with two identified when they agree at every point of $X$.

A **polynomial map** $\phi:X\to Y$ between varieties is one given by polynomials in the coordinates. It induces the **pullback** $\phi^*:k[Y]\to k[X]$, $g\mapsto g\circ\phi$, a ring homomorphism. The arrows reverse direction. The central fact, which we state rather than prove, is that this correspondence is perfect: polynomial maps $X\to Y$ correspond exactly to $k$-algebra homomorphisms $k[Y]\to k[X]$. In particular $X$ and $Y$ are **isomorphic**, meaning there are polynomial maps both ways that are mutually inverse, if and only if $k[X]\cong k[Y]$.

> [!example] The parabola is a line
> The parabola $P=V(y-x^2)$ and the line $k$ are isomorphic via $t\mapsto(t,t^2)$ and its inverse $(x,y)\mapsto x$, both polynomial. Correspondingly $k[x,y]/\langle y-x^2\rangle\cong k[t]$: replace each $y$ by $x^2$.

The cusp looks just as tame at first. The map $t\mapsto(t^2,t^3)$ is a bijection from the line onto $V(y^2-x^3)$, and it is polynomial. Is the cusp a line too?

> [!theorem] The cusp is not a line
> There is no polynomial map from the cusp $C=V(y^2-x^3)$ to the line that inverts $\phi(t)=(t^2,t^3)$. In fact $C$ is not isomorphic to the line by any polynomial maps.

*Proof of the first claim.* An inverse would be a polynomial $h(x,y)$ with $h(t^2,t^3)=t$ for all $t$. But every monomial $x^ay^b$ becomes $t^{2a+3b}$, and $2a+3b$ takes the values $0,2,3,4,5,\dots$ and *never* the value 1. So $h(t^2,t^3)$ has no $t^1$ term and cannot equal $t$. $\blacksquare$

The same computation shows what the coordinate ring of the cusp *is*. Pulling back along $\phi$ identifies $k[C]$ with $k[t^2,t^3]$, the polynomials in $t$ with no linear term. That ring is not $k[t]$. The quickest invariant that tells them apart is the dimension of the tangent space $\mathfrak{m}/\mathfrak{m}^2$ at the origin: it is 2 for the cusp, spanned by $x$ and $y$, and 1 at every point of a line. So no isomorphism of any kind exists. The singularity is visible in the algebra as a "missing" function, $t=y/x$, which is defined on the cusp away from the origin but is not a polynomial.

### Bézout's theorem, and the missing points

Two lines can be parallel, a circle and a line can miss each other, and a tangent line meets a parabola only once. To count intersections exactly, three repairs are needed.

1. **Complex points.** A circle and a line that "miss" meet at complex conjugate points.
2. **Points at infinity.** Parallel lines meet at infinity. Formally, the **projective plane** $\mathbb{P}^2$ is the set of nonzero triples $(x:y:z)$ up to scaling. The ordinary point $(x,y)$ is $(x:y:1)$, and the points with $z=0$ form the line at infinity. A curve $f(x,y)=0$ of degree $d$ becomes the homogeneous equation $F(x,y,z)=z^df(x/z,y/z)=0$.
3. **Multiplicity.** A tangency counts twice, just as a double root of a polynomial counts twice.

> [!theorem] Bézout's theorem
> Two projective plane curves of degrees $d$ and $e$ over $\mathbb{C}$ with no common component meet in exactly $de$ points, counted with multiplicity.

The general proof needs the resultant or intersection theory. The case $d=1$ is a short argument that shows exactly why the three repairs are needed.

*Proof for a line and a curve.* Parametrise the projective line through two distinct points $P$ and $Q$ as $sP+tQ$, with $(s:t)\in\mathbb{P}^1$. Substituting into the homogeneous $F$ of degree $d$ gives $g(s,t)=F(sP+tQ)$, a homogeneous polynomial of degree $d$ in $(s,t)$. If $g\equiv0$, the line is a component of the curve. Otherwise $g$ factors over $\mathbb{C}$ into $d$ linear factors $\prod_{i=1}^d(\beta_is-\alpha_it)$, by the fundamental theorem of algebra ([Ch. 18](../../part-5-transforms/18-complex-analysis.md)) applied to $g(s,1)$ together with the factors of $t$ that account for roots "at $s=\infty$". Each factor is a point $(\alpha_i:\beta_i)$ of the line lying on the curve, and repeated factors are tangencies. So there are exactly $d$ intersection points counted with multiplicity. $\blacksquare$

Without complex numbers, $g$ might have no real factors. Without the projective line, the factors of $t$ would be roots "at infinity" and would be lost. Without multiplicity, repeated factors would be undercounted. The theorem is clean only after all three repairs.

**The answer to the hook.** Any circle has equation $x^2+y^2+Dx+Ey+F=0$, which homogenises to

$$
x^2+y^2+Dxz+Eyz+Fz^2=0 .
$$

On the line at infinity, $z=0$, this reduces to $x^2+y^2=0$, that is, $(x+jy)(x-jy)=0$. So *every circle*, whatever its centre and radius, passes through the two points

$$
I=(1:j:0),\qquad J=(1:-j:0),
$$

the **circular points at infinity**. Two circles share these two points, which leaves $4-2=2$ finite intersections. That is why the robot arm has at most two elbow positions, and why subtracting the two circle equations killed both squares: the difference of two circle equations has no $x^2+y^2$ term, so it is a line, the line through the two finite intersections. Conversely, a conic through $I$ and $J$ is a circle. That is how projective geometers define a circle.

![Three panels: an ellipse and a hyperbola meeting at four red points; two circles meeting at the robot arm's two elbow positions, with the note that the missing two points are (1 : plus or minus j : 0) at infinity; a parabola and a circle that are tangent at the origin, which counts twice, and cross at two more points](../../images/ch94-bezout-conics.svg)

*Bézout's count for conics. Left: four ordinary points. Middle: the robot arm's two circles meet at two real points; the other two are the circular points at infinity, shared by every circle. Right: $y=x^2$ and $x^2+(y-1)^2=1$ touch at the origin with multiplicity 2 and cross at $(\pm1,1)$, so $2+1+1=4$.*

> [!mystery] Why is the chord-and-tangent rule associative?
> Bézout says a line meets a cubic curve in exactly three points. On an elliptic curve that lets you "add" two points $P$ and $Q$: draw the line through them, find the third intersection, and reflect it. Commutativity is obvious, since the line through $P$ and $Q$ is the line through $Q$ and $P$. But why should $(P+Q)+R=P+(Q+R)$? Checked by brute algebra, it is a famously ugly computation. The conceptual proof uses Bézout again, this time for two cubics meeting in nine points: if eight of those points lie on a third cubic, so does the ninth. Why that is true, and how the group law becomes the backbone of modern cryptography, is in [Ch. 78](../h-cryptography-coding/78-elliptic-curves.md).

## 94.4 The algorithm

**Gröbner bases: Gaussian elimination for polynomials.** For linear equations, elimination ([§8.7](../../part-3-linear-algebra/08-matrices-and-linear-systems.md#87-gaussian-elimination)) turns any system into a triangular one: the last equation involves one unknown, the one before it two, and so on. A **Gröbner basis** does the same for polynomial systems. Fix a *monomial order*, a rule for which term of a polynomial is "leading". In **lex order** with $u>v$, any term containing a higher power of $u$ wins. Then:

```text
Buchberger(F):                     # F = list of polynomials generating the ideal
    G = F
    repeat:
        for each pair (f, g) in G:
            m = lcm(LT(f), LT(g))              # least common multiple of the leading terms
            S = (m / LT(f)) * f - (m / LT(g)) * g    # the S-polynomial cancels the leading terms
            r = remainder of S on division by G
            if r != 0: add r to G
    until no pair produces a new remainder
    return G
```

Buchberger's theorem (1965) says the loop always stops, by Hilbert's basis theorem, and that the result $G$ generates the same ideal. With $G$ in hand, every polynomial has a *unique* remainder on division by $G$. So $f\in\langle F\rangle$ if and only if its remainder is 0, and the ideal membership problem is solved. In lex order, the **elimination theorem** says the elements of $G$ that involve only the last variable generate all consequences in that variable alone: the system has been triangularised.

**The robot arm, by Buchberger.** With $f_1=u^2+v^2-4$ and $f_2=(u-2)^2+(v-1)^2-1=u^2-4u+v^2-2v+4$, both leading terms are $u^2$. The S-polynomial is $f_1-f_2=4u+2v-8$, which we halve to $g=2u+v-4$. Dividing $f_1$ by $g$ (substitute $u=(4-v)/2$) leaves the remainder $\tfrac14(5v^2-8v)$. The reduced Gröbner basis is

$$
\{\,2u+v-4,\ \ 5v^2-8v\,\}:
$$

a quadratic in $v$ alone, $v(5v-8)=0$, and then a linear equation for $u$. These are the two elbows, $(2,0)$ and $(1.2,1.6)$, found by the same steps a person would take, but now as an algorithm that works for any number of equations and unknowns.

**Costs and failure modes.** Gröbner bases are exact and complete, and computer algebra systems compute them routinely for kinematics and geometric constraint problems. Their worst case is brutal, though: the degrees in a Gröbner basis can be doubly exponential in the number of variables. Over floating-point numbers they are numerically fragile, because whether a coefficient is exactly zero changes the answer. Engineers who need *numbers* rather than formulas often use **resultants** instead (§94.5), or **homotopy continuation**. That method starts from a simple system whose $de\cdots$ solutions are known, by Bézout, and deforms it into the target system while tracking every solution path with an ODE solver ([§22.9](../../part-7-numerical-computational/22-numerical-methods.md#229-numerical-solution-of-odes)). The general six-joint robot arm has 16 inverse-kinematics solutions, a count found this way and then proved.

> [!mystery] Doubly exponential, yet it works
> Mayr and Meyer (1982) built ideals for which *any* algorithm deciding ideal membership needs space exponential in the number of variables, and Gröbner bases with doubly exponential degrees. Yet the systems that come from engineering are solved in seconds. Is the worst case rare, or are real problems special? How problems can be provably hard in the worst case and easy in practice is a theme of [Ch. 30](../a-algorithms/30-limits-of-computation.md).

## 94.5 The experiment

**Watch singularities form.** Desmos plots implicit equations directly. The curve $y^2=x^3-3x+b$ is shown for $b=2$, where it has a node, together with two nearby smooth members of the same family.

```desmos-graph
left=-3; right=3; bottom=-3; top=3
---
y^2=x^3-3x+2|red
y^2=x^3-3x+2.5|blue
y^2=x^3-3x+1.5|green
(1,0)|label:node
```

*At $b=2$ (red) the loop and the branch touch at the node $(1,0)$. Slightly above (blue, $b=2.5$) they merge into one smooth piece. Slightly below (green, $b=1.5$) they separate into an oval and a branch.*

Things to change:

- Replace $-3x$ by $0x$ and set all three constants to $0$, $0.3$ and $-0.3$. The red curve becomes the cusp $y^2=x^3$. Its neighbours are smooth, but they do not split into two pieces.
- The singular curves are those with $4a^3+27b^2=0$. Use this to predict, then check, which $b$ gives a node when $a=-12$.
- Add the lines `x^2+y^2=4`, `(x-2)^2+(y-1)^2=1` and `2x+y=4` to see the robot arm: the two circles and the line through their two finite intersections.

**Count every intersection.** This script finds *all* common points of two conics $Ax^2+Bxy+Cy^2+Dx+Ey+F=0$, complex ones included. It eliminates $y$ with the **resultant** of the two quadratics in $y$, a polynomial in $x$ alone that vanishes exactly where they share a root, and then recovers $y$.

```python
import numpy as np
P = np.poly1d

def intersect(c1, c2, tol=1e-4):
    """Common points of two conics A x^2 + B xy + C y^2 + D x + E y + F = 0."""
    def in_y(c):                                  # coefficients of y^2, y^1, y^0 as polys in x
        A, B, C, D, E, F = c
        return P([C]), P([B, E]), P([A, D, F])
    a2, a1, a0 = in_y(c1)
    b2, b1, b0 = in_y(c2)
    res = (a2 * b0 - a0 * b2)**2 - (a2 * b1 - a1 * b2) * (a1 * b0 - a0 * b1)   # resultant in y
    res = P(np.where(abs(res.coeffs) < 1e-12, 0, res.coeffs))
    print("  resultant has degree", res.order, "in x")
    pts = []
    for x in np.unique(np.round(res.roots, 5)):   # each x may carry several points
        for y in np.unique(np.round(np.roots([a2(x), a1(x), a0(x)]), 5)):
            if abs(b2(x) * y**2 + b1(x) * y + b0(x)) < tol * (1 + abs(x)**2 + abs(y)**2):
                pts.append((complex(x), complex(y)))
    return pts

ellipse   = (1, 0, 4, 0, 0, -4)        # x^2 + 4y^2 = 4
hyperbola = (1, 0, -1, 0, 0, -1)       # x^2 - y^2 = 1
circle1   = (1, 0, 1, 0, 0, -4)        # x^2 + y^2 = 4
circle2   = (1, 0, 1, -4, -2, 4)       # (x-2)^2 + (y-1)^2 = 1

for name, c1, c2 in [("ellipse & hyperbola", ellipse, hyperbola),
                     ("two circles", circle1, circle2)]:
    print(name)
    for x, y in intersect(c1, c2):
        print(f"   x = {x.real:+.4f}{x.imag:+.4f}j,  y = {y.real:+.4f}{y.imag:+.4f}j")
```

For the ellipse and hyperbola the resultant has degree 4 and the script prints the four points $(\pm1.2649,\pm0.7746)$, which are $(\pm\sqrt{8/5},\pm\sqrt{3/5})$. For the two circles the resultant collapses to degree 2: two solutions, the elbows $(1.2,1.6)$ and $(2,0)$. The lost degree *is* the pair of circular points at infinity. Things to change:

- Replace the hyperbola by the circle $(x-4)^2+y^2=1$, `(1, 0, 1, -8, 0, 15)`, which misses the ellipse entirely in the real plane. The script still finds four points: $(8/3,\pm0.8819j)$ and $(8,\pm3.8730j)$.
- Give `circle2` a different $x^2$ coefficient, for instance `(1.1, 0, 1, -4, -2, 4)`. It is now an ellipse, it no longer passes through $(1:\pm j:0)$, and the resultant degree jumps back to 4. The two new solutions are huge, about $(38.4\mp20.9j,\ -20.9\mp38.3j)$, with $y/x\approx\mp j$: the points that were at infinity have come in from $(1:\mp j:0)$. Try `1.01` and watch them recede.
- Make two conics tangent, for instance $y=x^2$ and $x^2+(y-1)^2=1$, `(1, 0, 0, 0, -1, 0)` and `(1, 0, 1, 0, -2, 0)`. Only three points print. Print `res.roots` to see the fourth: $x=0$ appears twice.

## 94.6 The challenge

1. **Intuition.** Three generic quadric surfaces in space meet in $2\cdot2\cdot2=8$ points. Three spheres, as in GPS trilateration, meet in at most 2. Where are the other six? (Hint: what does every sphere have in common at infinity?)
2. **Derivation.** Over $\mathbb{C}$, show that the substitution $u=x+jy$, $v=x-jy$ turns the circle $x^2+y^2=1$ into the hyperbola $uv=1$, and deduce that its coordinate ring is $\mathbb{C}[u,u^{-1}]$. Then use the lines through $(-1,0)$ to derive the rational parametrisation $x=(1-t^2)/(1+t^2)$, $y=2t/(1+t^2)$, and show that it gives every Pythagorean triple.
3. **Implementation.** Extend the resultant script to two *cubic* curves, using the $6\times6$ Sylvester determinant (evaluated at enough values of $x$ and interpolated) to get a resultant of degree 9. Verify Bézout's count of 9 on random cubics, then find a pair of cubics for which the degree drops, and identify the common points at infinity.
4. **Open investigation.** Draw the six coefficients of each of two conics independently from a standard Gaussian. Estimate by simulation the probabilities that they meet in 0, 2 or 4 real points. Then try two random curves of degree $d$ for $d=2,\dots,6$ and plot the average number of real intersection points against $d$. For one special distribution (Kostlan's) the expected count is exactly $\sqrt{de}$ (Shub and Smale, 1993). How does your distribution compare, and can you guess its growth rate?

## 94.7 Loose threads

**Why SNARKs trust random points.** The zero-knowledge proof systems of [Ch. 82](../h-cryptography-coding/82-zero-knowledge-proofs.md) encode a whole computation as an identity between polynomials and then check it at one random point. The geometric fact behind this is that a hypersurface is *thin*. The **Schwartz–Zippel lemma** says that if $f\ne0$ has total degree $d$ and $\mathbf{r}$ is drawn uniformly from $S^n$ for a finite set $S$, then $P(f(\mathbf{r})=0)\le d/\lvert S\rvert$. The proof is a short induction. For $n=1$, a nonzero polynomial of degree $d$ has at most $d$ roots. For $n>1$, write $f=\sum_{i\le k}x_1^if_i(x_2,\dots,x_n)$ with $f_k\neq0$. By induction, $f_k(r_2,\dots,r_n)=0$ with probability at most $(d-k)/\lvert S\rvert$. Otherwise $f$ is a nonzero polynomial of degree $k$ in $x_1$, which vanishes at $r_1$ with probability at most $k/\lvert S\rvert$. Adding the two gives $d/\lvert S\rvert$. Over a field with $q\approx2^{255}$ elements, a false identity of degree a million survives a random test with probability below $2^{-235}$. Seen this way, a variety is a set of measure almost zero, and a random point almost surely lies off it.

**Fermat's Last Theorem.** [Ch. 57](../e-abstract-math/57-number-theory.md) asked how geometry enters Fermat's problem. A solution of $a^n+b^n=c^n$ in positive integers is a rational point $(a/c,b/c)$ on the curve $x^n+y^n=1$. The **genus** of a smooth curve, the number of holes in its complex points viewed as a surface, controls how many rational points it can have. Genus 0 curves, like the circle, have none or infinitely many, which is why there are infinitely many Pythagorean triples. Genus 1 curves are elliptic curves, whose rational points form a finitely generated group. Curves of genus at least 2 have only *finitely* many rational points (Faltings, 1983). The Fermat curve has genus $(n-1)(n-2)/2$, which is at least 3 once $n\ge4$, so Faltings already gives finitely many solutions. Wiles' proof (1995) that there are *none* goes through an elliptic curve, Frey's curve $y^2=x(x-a^n)(x+b^n)$, which would be too strange to exist if a solution did.

> [!mystery] Can every singularity be smoothed away?
> The cusp can be repaired. Substitute $y=tx$ into $y^2=x^3$ and divide by $x^2$: the result is $t^2=x$, a smooth parabola. This "blow-up" replaces the singular point by the directions of approach to it. Hironaka proved in 1964, in a 200-page paper that earned a Fields Medal, that every variety over a field of characteristic 0, such as $\mathbb{R}$ or $\mathbb{C}$, can be made smooth by repeated blow-ups. Over fields of characteristic $p$, the setting of cryptography and coding theory, the same statement is known only in dimension at most 3. In higher dimensions it is an **open problem**, and one of the most sought-after in the field.

The other threads lead into the rest of the track. The variety of tensors of border rank at most $r$ from [Ch. 93](93-tensor-mathematics.md) is a secant variety, and finding its equations is algebraic geometry at the research frontier. Symmetry groups acting on varieties, and on the coordinate rings that describe them, bring in the representation theory of [Ch. 95](95-representation-theory.md). Grothendieck's reformulation of the whole subject, which replaces varieties by "schemes" so that any commutative ring has a geometry, shows how mathematicians build new language when the old one runs out ([Ch. 98](98-mathematical-research.md)).

## Further reading

- D. A. Cox, J. Little and D. O'Shea, *Ideals, Varieties, and Algorithms*, 4th ed., Springer, 2015. The best first book: Gröbner bases, elimination and the Nullstellensatz, with robotics examples.
- M. Reid, *Undergraduate Algebraic Geometry*, Cambridge University Press, 1988. Short, witty and geometric: conics, cubics and the dictionary.
- F. Kirwan, *Complex Algebraic Curves*, Cambridge University Press, 1992. Bézout's theorem, singularities and the topology of complex curves.
- W. Fulton, *Algebraic Curves*, freely available online. A classic, concise treatment of intersection multiplicity and Bézout.
- A. J. Sommese and C. W. Wampler, *The Numerical Solution of Systems of Polynomials Arising in Engineering and Science*, World Scientific, 2005. Homotopy continuation and robot kinematics.
