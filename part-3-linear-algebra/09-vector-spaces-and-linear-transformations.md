# Chapter 9 — Vector Spaces and Linear Transformations

[← Previous](08-matrices-and-linear-systems.md) · [Table of Contents](../00-front/table-of-contents.md) · [Next →](10-eigenvalues-and-eigenvectors.md)

> **Prerequisites:** [Ch. 2](../part-1-foundations/02-complex-numbers.md) (complex numbers), [Ch. 5](../part-2-calculus/05-integral-calculus.md) (definite integrals), [Ch. 8](08-matrices-and-linear-systems.md) (matrices, Gaussian elimination, rank)
> **Learning outcomes:** After this chapter you can
> - recognise a vector space — whether its elements are arrows, columns of numbers, polynomials, matrices or signals — and test whether a subset is a subspace;
> - decide linear independence, find a basis and the dimension of a space, and compute coordinates and change-of-basis matrices;
> - represent a linear transformation by a matrix, and compute bases for its null space, column space, row space and left null space using the rank–nullity theorem;
> - work with inner products on $\mathbb{R}^n$, $\mathbb{C}^n$ and on function spaces, and use orthogonality, norms and the Cauchy–Schwarz inequality;
> - project onto lines and subspaces, build orthonormal bases with Gram–Schmidt (including QR and Legendre polynomials);
> - solve overdetermined systems by least squares and fit lines, polynomials and exponentials to measured data.

Chapter 8 treated a matrix as a table of coefficients and a linear system as something to be solved by elimination. This chapter changes the point of view. A matrix is a *machine* that maps one vector space into another; the questions "does $\mathbf{A}\mathbf{x}=\mathbf{b}$ have a solution?" and "is it unique?" become questions about two subspaces attached to that machine. The same language then applies, almost word for word, to objects that do not look like columns of numbers at all: the set of solutions of a linear differential equation, the set of polynomials of degree at most $n$, the set of finite-energy signals. That is the real pay-off. When you later expand a periodic waveform in a Fourier series ([Ch. 15](../part-5-transforms/15-fourier-series.md)), you will be doing exactly what this chapter does with arrows in $\mathbb{R}^3$: projecting a vector onto an orthogonal basis. When you fit a calibration curve to noisy measurements, you will be projecting a data vector onto the column space of a matrix.

The chapter is organised in two halves. Sections 9.1–9.6 build the *algebraic* structure: spaces, subspaces, independence, bases, linear maps and the four fundamental subspaces. Sections 9.7–9.11 add *geometry* — lengths, angles and orthogonality — and use it for approximation: projection, Gram–Schmidt and least squares.

## 9.1 Vector spaces

### Motivation

An engineer meets many different objects that can be *added* and *scaled*: forces and field vectors, node-voltage vectors of a circuit, sampled signals $x[0],\dots,x[N-1]$, continuous waveforms $x(t)$, polynomials used in filter design, matrices of network parameters. In each case the same handful of algebraic rules hold — addition is commutative, there is a zero, scaling distributes over addition, and so on. Rather than prove every theorem separately for arrows, for sampled signals and for waveforms, we prove it once for any set obeying those rules. That set is a **vector space**.

The most important consequence is **superposition**. If $x_1(t)$ and $x_2(t)$ are responses of a linear circuit, so is $a x_1(t)+b x_2(t)$. The set of all possible responses is therefore closed under linear combinations — it is a vector space — and everything we learn about bases and dimension tells us how many "independent" responses we need to describe them all.

### Definition

> [!definition] Vector space
> A **vector space** over a field $\mathbb{F}$ (in this book $\mathbb{F}=\mathbb{R}$ or $\mathbb{C}$) is a set $V$ of objects called **vectors**, together with two operations
> - vector addition $\mathbf{u}+\mathbf{v}\in V$ for all $\mathbf{u},\mathbf{v}\in V$,
> - scalar multiplication $c\,\mathbf{u}\in V$ for all $c\in\mathbb{F}$, $\mathbf{u}\in V$,
>
> such that for all $\mathbf{u},\mathbf{v},\mathbf{w}\in V$ and $a,b\in\mathbb{F}$:
> 1. $\mathbf{u}+\mathbf{v}=\mathbf{v}+\mathbf{u}$ (commutativity)
> 2. $(\mathbf{u}+\mathbf{v})+\mathbf{w}=\mathbf{u}+(\mathbf{v}+\mathbf{w})$ (associativity)
> 3. there is a **zero vector** $\mathbf{0}\in V$ with $\mathbf{u}+\mathbf{0}=\mathbf{u}$
> 4. each $\mathbf{u}$ has an **additive inverse** $-\mathbf{u}$ with $\mathbf{u}+(-\mathbf{u})=\mathbf{0}$
> 5. $a(\mathbf{u}+\mathbf{v})=a\mathbf{u}+a\mathbf{v}$
> 6. $(a+b)\mathbf{u}=a\mathbf{u}+b\mathbf{u}$
> 7. $a(b\mathbf{u})=(ab)\mathbf{u}$
> 8. $1\,\mathbf{u}=\mathbf{u}$
>
> When $\mathbb{F}=\mathbb{R}$ we speak of a **real** vector space, when $\mathbb{F}=\mathbb{C}$ of a **complex** vector space.

The two **closure** requirements (sums and scalar multiples stay inside $V$) are the ones that fail in practice; the eight axioms are almost always inherited from ordinary arithmetic of numbers or functions. A few consequences follow directly from the axioms and are used constantly:

$$
0\,\mathbf{u}=\mathbf{0},\qquad c\,\mathbf{0}=\mathbf{0},\qquad (-1)\mathbf{u}=-\mathbf{u},\qquad c\,\mathbf{u}=\mathbf{0}\ \Rightarrow\ c=0\ \text{or}\ \mathbf{u}=\mathbf{0}.
$$

*Proof of the first.* By axiom 6, $0\mathbf{u}=(0+0)\mathbf{u}=0\mathbf{u}+0\mathbf{u}$. Add $-(0\mathbf{u})$ to both sides: $\mathbf{0}=0\mathbf{u}$. The others are similar.

### The standard examples

| Space | Vectors | Addition and scaling | Zero vector |
|---|---|---|---|
| $\mathbb{R}^n$ | columns $(x_1,\dots,x_n)^\mathsf{T}$ of reals | componentwise | $(0,\dots,0)^\mathsf{T}$ |
| $\mathbb{C}^n$ | columns of complex numbers (phasor vectors) | componentwise, complex scalars | $(0,\dots,0)^\mathsf{T}$ |
| $\mathbb{R}^{m\times n}$ | $m\times n$ real matrices | entrywise | zero matrix |
| $\mathcal{P}_n$ | polynomials $a_0+a_1t+\dots+a_nt^n$ of degree $\le n$ | add/scale coefficients | the zero polynomial |
| $C[a,b]$ | continuous real functions on $[a,b]$ | $(f+g)(t)=f(t)+g(t)$, $(cf)(t)=cf(t)$ | $f(t)\equiv 0$ |
| $\ell^2$ / $\mathbb{R}^N$ | sampled signals $x[n]$ | samplewise | the all-zero sequence |
| $L^2(a,b)$ | functions with $\int_a^b\lvert f(t)\rvert^2\,dt<\infty$ (finite-energy signals) | pointwise | $f\equiv0$ |

The table shows how broad the idea is. A polynomial is "a vector" because two polynomials can be added and scaled and the axioms hold; it does not matter that we never draw it as an arrow. The function spaces $C[a,b]$ and $L^2(a,b)$ are the setting for Fourier analysis: a periodic signal will be treated as a single point in an infinite-dimensional space.

> [!example] Worked Example 9.1 — Checking the axioms for $\mathcal{P}_2$
> Show that $\mathcal{P}_2=\{a_0+a_1t+a_2t^2 : a_k\in\mathbb{R}\}$ is a real vector space.
>
> **Closure.** If $p(t)=a_0+a_1t+a_2t^2$ and $q(t)=b_0+b_1t+b_2t^2$, then
> $$
> (p+q)(t)=(a_0+b_0)+(a_1+b_1)t+(a_2+b_2)t^2\in\mathcal{P}_2,\qquad (cp)(t)=ca_0+ca_1t+ca_2t^2\in\mathcal{P}_2 .
> $$
> **Axioms.** Each axiom reduces to the same axiom for real numbers applied coefficient by coefficient. For instance, commutativity: $(p+q)(t)$ and $(q+p)(t)$ have coefficients $a_k+b_k=b_k+a_k$. The zero vector is $0+0t+0t^2$; the inverse of $p$ is $-a_0-a_1t-a_2t^2$.
>
> **Remark.** The set of polynomials of degree *exactly* 2 is **not** a vector space: $t^2+t$ and $-t^2$ both have degree 2 but their sum $t$ does not, and the zero polynomial is missing.

### Non-examples

Most failures are failures of closure or of the zero vector:

- **The first quadrant** $\{(x,y): x\ge0,\ y\ge0\}$: closed under addition but not under scaling by $-1$.
- **A line not through the origin**, e.g. $\{(x,y): y=2x+1\}$: $(0,1)$ and $(1,3)$ are in it but their sum $(1,4)$ is not; and $(0,0)$ is not in it. Such a set is an **affine** set — a vector space shifted away from the origin.
- **Unit vectors** $\{\mathbf{x}:\lVert\mathbf{x}\rVert=1\}$: the sum of two unit vectors is generally not a unit vector.
- **Positive voltages** $\{v(t): v(t)>0\}$: not closed under multiplication by $-1$, no zero vector.
- **Solutions of a non-homogeneous equation**, e.g. all $y(t)$ with $y''+3y'+2y=1$ (worked example below).

> [!example] Worked Example 9.2 — Solution sets of differential equations
> (a) Show that the set $S_0$ of twice-differentiable functions satisfying $y''+3y'+2y=0$ is a vector space.
> (b) Show that the set $S_1$ satisfying $y''+3y'+2y=1$ is not.
>
> **(a)** $S_0$ is a subset of the vector space of all twice-differentiable functions, so the axioms are inherited; we need only closure (this is the subspace test of Section 9.2). Let $y_1,y_2\in S_0$ and $a,b\in\mathbb{R}$, and put $y=ay_1+by_2$. Because differentiation is linear,
> $$
> y''+3y'+2y=a\,(y_1''+3y_1'+2y_1)+b\,(y_2''+3y_2'+2y_2)=a\cdot0+b\cdot0=0 .
> $$
> So $y\in S_0$. The zero function is in $S_0$. Hence $S_0$ is a vector space. (In [Ch. 12](../part-4-differential-equations/12-linear-differential-equations.md) you will find $S_0=\{c_1e^{-t}+c_2e^{-2t}\}$.)
>
> **(b)** Take $y_1=y_2=\tfrac12$ (a constant: $0+0+2\cdot\tfrac12=1$, so it is in $S_1$). Then $y_1+y_2=1$ gives $0+0+2=2\neq1$. Closure fails; also $y\equiv0$ gives $0\neq1$, so the zero vector is absent. $S_1$ is an affine set: $S_1=\tfrac12+S_0$, a particular solution plus the homogeneous space.

This "particular solution plus vector space" structure is the general shape of the solution set of every consistent linear problem $L(\mathbf{x})=\mathbf{b}$, whether $L$ is a matrix or a differential operator. It reappears in Section 9.6.

### Sinusoids of a fixed frequency

The set

$$
\mathcal{S}_\omega=\{A\cos(\omega t+\phi): A\ge0,\ \phi\in\mathbb{R}\}
$$

of all sinusoids at a fixed angular frequency $\omega$ is closed under addition — the sum of two sinusoids at the same frequency is another sinusoid at that frequency — and under real scaling (a negative scale factor is a phase shift by $\pi$). With the zero signal ($A=0$) included it is a real vector space. Writing $A\cos(\omega t+\phi)=a\cos\omega t+b\sin\omega t$ with $a=A\cos\phi$, $b=-A\sin\phi$ shows that each sinusoid is fixed by the pair $(a,b)\in\mathbb{R}^2$; equivalently by the **phasor** $\mathbf{X}=Ae^{j\phi}=a-jb\in\mathbb{C}$. Adding sinusoids corresponds to adding phasors. This is the vector-space reason why phasor analysis ([Ch. 2](../part-1-foundations/02-complex-numbers.md)) works: the map "sinusoid $\mapsto$ phasor" respects addition and real scaling, so it is a linear bijection between $\mathcal{S}_\omega$ and $\mathbb{C}$ viewed as the real space $\mathbb{R}^2$ (Section 9.5 calls this an *isomorphism*).

> [!example] Worked Example 9.3 — Adding sinusoids as vectors
> Write $x(t)=3\cos(\omega t)+4\cos(\omega t-\pi/2)$ as a single sinusoid.
>
> **As vectors in the basis $\{\cos\omega t,\ \sin\omega t\}$.** $\cos(\omega t-\pi/2)=\sin\omega t$, so $x=3\cos\omega t+4\sin\omega t$, i.e. coordinates $(a,b)=(3,4)$.
>
> **Convert.** $A=\sqrt{a^2+b^2}=5$ and, from $a=A\cos\phi$, $b=-A\sin\phi$: $\cos\phi=0.6$, $\sin\phi=-0.8$, so $\phi=-0.9273$ rad $=-53.13^\circ$.
> $$
> x(t)=5\cos(\omega t-53.13^\circ).
> $$
> **Phasor check.** $3e^{j0}+4e^{-j\pi/2}=3-4j=5\angle-53.13^\circ$. Same answer, because the phasor map is linear.

> [!warning] Common Pitfall — "vector" does not mean "arrow"
> In this chapter a vector is *any element of a vector space*. The word says nothing about geometry. A polynomial, a $3\times3$ matrix and a voltage waveform are all vectors in suitable spaces. Conversely, a list of numbers is only a vector if the operations you intend to perform on it are componentwise addition and scaling: a list of (resistance, temperature, date) readings is a data record, not a meaningful element of $\mathbb{R}^3$.

> [!tip] Engineering Note — why EE lives in complex vector spaces
> Steady-state AC analysis replaces each node voltage by a phasor, so the unknown vector of nodal analysis lies in $\mathbb{C}^n$ and the admittance matrix has complex entries. Every definition in this chapter works over $\mathbb{C}$; the one change needed is in the inner product (Section 9.8), where complex conjugation is required so that lengths remain real and positive.

## 9.2 Subspaces

### Motivation

Inside a big vector space we constantly meet smaller sets that are vector spaces in their own right: the plane of all current vectors satisfying Kirchhoff's current law, the set of band-limited signals, the polynomials that vanish at a given point, the set of all outputs a system can produce. Recognising these as **subspaces** lets us use all of the vector-space machinery inside them.

### Definition and the subspace test

> [!definition] Subspace
> A subset $W$ of a vector space $V$ is a **subspace** if $W$ is itself a vector space under the operations of $V$.

Checking eight axioms again is unnecessary: they are inherited from $V$. Only closure and the presence of the zero vector can fail.

> [!theorem] Subspace test
> A subset $W\subseteq V$ is a subspace if and only if
> 1. $\mathbf{0}\in W$;
> 2. $\mathbf{u},\mathbf{v}\in W\Rightarrow \mathbf{u}+\mathbf{v}\in W$ (closed under addition);
> 3. $\mathbf{u}\in W,\ c\in\mathbb{F}\Rightarrow c\mathbf{u}\in W$ (closed under scaling).
>
> Conditions 2 and 3 can be merged: $W$ is a subspace iff $\mathbf{0}\in W$ and $a\mathbf{u}+b\mathbf{v}\in W$ for all $\mathbf{u},\mathbf{v}\in W$ and scalars $a,b$.

*Proof sketch.* Necessity is clear. For sufficiency, axioms 1, 2 and 5–8 hold in $W$ because they hold in $V$ and the operations are the same. Axiom 3 is condition 1. For axiom 4, if $\mathbf{u}\in W$ then $-\mathbf{u}=(-1)\mathbf{u}\in W$ by condition 3.

Every vector space $V$ has two **trivial subspaces**: $\{\mathbf{0}\}$ and $V$ itself. The subspaces of $\mathbb{R}^3$ are exactly: $\{\mathbf{0}\}$, lines through the origin, planes through the origin, and $\mathbb{R}^3$. The phrase "through the origin" is essential — it is condition 1.

> [!example] Worked Example 9.4 — Four subspace tests
> Decide which are subspaces.
>
> **(a)** $W_a=\{(x,y,z)\in\mathbb{R}^3: x+2y-z=0\}$. Zero: $0+0-0=0$ ✓. If $\mathbf{u}=(x_1,y_1,z_1)$ and $\mathbf{v}=(x_2,y_2,z_2)$ satisfy the equation, then $a\mathbf{u}+b\mathbf{v}$ gives $a(x_1+2y_1-z_1)+b(x_2+2y_2-z_2)=0$ ✓. **Subspace** (a plane through the origin).
>
> **(b)** $W_b=\{(x,y,z): x+2y-z=4\}$. $(0,0,0)$ gives $0\neq4$. **Not a subspace** (a parallel plane shifted off the origin).
>
> **(c)** $W_c=\{p\in\mathcal{P}_3: p(1)=0\}$. The zero polynomial vanishes at 1 ✓. $(ap+bq)(1)=a\,p(1)+b\,q(1)=0$ ✓. **Subspace.**
>
> **(d)** $W_d=\{(x,y)\in\mathbb{R}^2: xy=0\}$ (the union of the two axes). $(1,0)$ and $(0,1)$ are in $W_d$ but $(1,1)$ is not. **Not a subspace.** A union of subspaces is usually not a subspace.

A useful pattern emerges: **the solution set of a homogeneous linear condition is always a subspace; an inhomogeneous condition (right-hand side $\neq0$) or a nonlinear condition usually is not.**

### Span

> [!definition] Linear combination and span
> A **linear combination** of $\mathbf{v}_1,\dots,\mathbf{v}_k$ is any vector $c_1\mathbf{v}_1+\dots+c_k\mathbf{v}_k$. The **span** of the set is the collection of all linear combinations:
> $$
> \operatorname{span}\{\mathbf{v}_1,\dots,\mathbf{v}_k\}=\{c_1\mathbf{v}_1+\dots+c_k\mathbf{v}_k : c_i\in\mathbb{F}\}.
> $$
> If $\operatorname{span}\{\mathbf{v}_1,\dots,\mathbf{v}_k\}=W$, the vectors **span** (or **generate**) $W$.

> [!theorem] A span is a subspace
> $\operatorname{span}\{\mathbf{v}_1,\dots,\mathbf{v}_k\}$ is a subspace of $V$, and it is the *smallest* subspace containing all the $\mathbf{v}_i$.

*Proof.* The zero vector is the combination with all $c_i=0$. A linear combination of two linear combinations is again a linear combination (collect coefficients). Any subspace containing every $\mathbf{v}_i$ must, by closure, contain every combination of them.

Deciding whether $\mathbf{b}\in\operatorname{span}\{\mathbf{v}_1,\dots,\mathbf{v}_k\}$ in $\mathbb{R}^n$ is the same as deciding whether the linear system

$$
\begin{bmatrix}\mathbf{v}_1&\mathbf{v}_2&\cdots&\mathbf{v}_k\end{bmatrix}\begin{bmatrix}c_1\\ \vdots\\ c_k\end{bmatrix}=\mathbf{b}
$$

is consistent — a Chapter 8 problem solved by Gaussian elimination.

> [!example] Worked Example 9.5 — Membership in a span
> Is $\mathbf{b}=(1,4,7)^\mathsf{T}$ in $\operatorname{span}\{\mathbf{v}_1,\mathbf{v}_2\}$ with $\mathbf{v}_1=(1,2,3)^\mathsf{T}$, $\mathbf{v}_2=(0,1,2)^\mathsf{T}$? What about $\mathbf{d}=(1,4,8)^\mathsf{T}$?
>
> Row-reduce the augmented matrix with both right-hand sides:
> $$
> \left[\begin{array}{cc|cc}1&0&1&1\\2&1&4&4\\3&2&7&8\end{array}\right]
> \xrightarrow{R_2-2R_1,\ R_3-3R_1}
> \left[\begin{array}{cc|cc}1&0&1&1\\0&1&2&2\\0&2&4&5\end{array}\right]
> \xrightarrow{R_3-2R_2}
> \left[\begin{array}{cc|cc}1&0&1&1\\0&1&2&2\\0&0&0&1\end{array}\right].
> $$
> For $\mathbf{b}$ the last row reads $0=0$: consistent, $c_1=1$, $c_2=2$. Check: $1(1,2,3)+2(0,1,2)=(1,4,7)$ ✓.
>
> For $\mathbf{d}$ the last row reads $0=1$: inconsistent, so $\mathbf{d}\notin$ the span. Geometrically, the span is the plane through the origin with normal $\mathbf{v}_1\times\mathbf{v}_2=(1,-2,1)$; indeed $1-8+7=0$ for $\mathbf{b}$ but $1-8+8=1\neq0$ for $\mathbf{d}$.

### Sums and intersections of subspaces

If $U$ and $W$ are subspaces of $V$, then

- $U\cap W$ is a subspace (check the three conditions — each holds in both);
- $U+W=\{\mathbf{u}+\mathbf{w}:\mathbf{u}\in U,\ \mathbf{w}\in W\}$ is a subspace, the smallest containing both;
- $U\cup W$ is a subspace only if one contains the other (Worked Example 9.4(d)).

When $U\cap W=\{\mathbf{0}\}$ the sum is called **direct**, written $U\oplus W$, and every vector in it splits *uniquely* as $\mathbf{u}+\mathbf{w}$. Example: every function on $[-a,a]$ splits uniquely into an even and an odd part,

$$
f(t)=\underbrace{\tfrac12\big[f(t)+f(-t)\big]}_{\text{even}}+\underbrace{\tfrac12\big[f(t)-f(-t)\big]}_{\text{odd}},
$$

and the only function that is both even and odd is $0$. So the function space is the direct sum of the even subspace and the odd subspace — the fact behind the cosine-only and sine-only Fourier series of [Ch. 15](../part-5-transforms/15-fourier-series.md).

> [!example] Worked Example 9.6 — Intersection of two planes
> Let $U=\{x+2y-z=0\}$ and $W=\{x-y+z=0\}$ in $\mathbb{R}^3$. Find $U\cap W$.
>
> Solve both equations together. Adding them: $2x+y=0$, so $y=-2x$; then $z=x+2y=-3x$. Thus
> $$
> U\cap W=\operatorname{span}\{(1,-2,-3)^\mathsf{T}\},
> $$
> a line through the origin. Check: $1-4+3=0$ ✓ and $1+2-3=0$ ✓. Equivalently, the line is along $\mathbf{n}_U\times\mathbf{n}_W=(1,2,-1)\times(1,-1,1)=(1,-2,-3)$.

### Signal subspaces

Function spaces have many natural subspaces that matter in signal processing:

- **Zero-mean signals** on $[0,T]$: $\{x:\int_0^T x\,dt=0\}$ — a homogeneous linear condition, so a subspace (the "AC part" of signals).
- **$T$-periodic signals**: $x(t+T)=x(t)$ is preserved by sums and scalings.
- **Trigonometric polynomials of order $N$**: $\operatorname{span}\{1,\cos\omega_0t,\sin\omega_0t,\dots,\cos N\omega_0t,\sin N\omega_0t\}$ — the space in which a truncated Fourier series lives.
- **Band-limited signals**: those whose Fourier transform vanishes for $|\omega|>W$ ([Ch. 16](../part-5-transforms/16-fourier-transforms.md)); closed under addition because the transform is linear.
- **Causal signals** $x(t)=0$ for $t<0$.

Signals with **peak amplitude at most 1 V**, by contrast, do not form a subspace (scaling by 2 breaks the bound).

> [!warning] Common Pitfall — forgetting the origin
> Before doing any work, check whether $\mathbf{0}$ belongs to the set. If it does not, the set is not a subspace and you are done. Most "almost-subspaces" met in practice (operating points, biased signals, regions bounded by tolerances) are affine sets or cones, not subspaces.

> [!tip] Engineering Note — KCL defines a subspace
> In a circuit with $b$ branches, collect the branch currents in $\mathbf{i}\in\mathbb{R}^b$. Kirchhoff's current law at each node is a homogeneous linear equation, so the set of all KCL-compatible current vectors is a subspace of $\mathbb{R}^b$. Its dimension (Section 9.6) is the number of independent loop currents — the number of unknowns in mesh analysis.

## 9.3 Linear independence

### Motivation

A spanning set may contain redundant vectors. If a third probe signal is just a combination of the first two, measuring it tells us nothing new; if one KVL equation is the sum of two others, writing it down adds no information and makes the system singular. **Linear independence** is the precise notion of "no redundancy". It decides whether a set of node equations has a unique solution, whether two waveforms can be separated by a receiver, and whether a set of functions can serve as a basis for a general solution of a differential equation.

### Definition

> [!definition] Linear independence
> Vectors $\mathbf{v}_1,\dots,\mathbf{v}_k$ are **linearly independent** if the only solution of
> $$
> c_1\mathbf{v}_1+c_2\mathbf{v}_2+\dots+c_k\mathbf{v}_k=\mathbf{0}
> $$
> is the **trivial** one $c_1=c_2=\dots=c_k=0$. Otherwise they are **linearly dependent**, and a nontrivial solution is called a **dependence relation**.

Equivalent formulations, each useful in a different setting:

- The set is dependent iff **at least one vector is a linear combination of the others**. (If $c_j\neq0$ in a dependence relation, solve for $\mathbf{v}_j$.)
- The set is independent iff **every vector in its span has a unique representation** $c_1\mathbf{v}_1+\dots+c_k\mathbf{v}_k$. (If two representations existed, their difference would be a nontrivial dependence relation.)
- In $\mathbb{F}^n$, with $\mathbf{V}=[\mathbf{v}_1\ \cdots\ \mathbf{v}_k]$, the set is independent iff $\mathbf{V}\mathbf{c}=\mathbf{0}$ has only $\mathbf{c}=\mathbf{0}$, iff $\operatorname{rank}\mathbf{V}=k$, iff every column of the echelon form of $\mathbf{V}$ has a pivot.
- For $k=n$ (square $\mathbf{V}$): independent iff $\det\mathbf{V}\neq0$.

Some quick facts:

- Any set containing $\mathbf{0}$ is dependent ($1\cdot\mathbf{0}=\mathbf{0}$).
- Two vectors are dependent iff one is a scalar multiple of the other.
- **More than $n$ vectors in $\mathbb{F}^n$ are always dependent**: $\mathbf{V}\mathbf{c}=\mathbf{0}$ then has more unknowns than equations, so a free variable exists.
- A subset of an independent set is independent; a superset of a dependent set is dependent.

> [!example] Worked Example 9.7 — Testing vectors in $\mathbb{R}^3$
> (a) Are $\mathbf{v}_1=(1,2,0)^\mathsf{T}$, $\mathbf{v}_2=(0,1,1)^\mathsf{T}$, $\mathbf{v}_3=(1,0,-2)^\mathsf{T}$ independent?
>
> Row-reduce $\mathbf{V}=[\mathbf{v}_1\ \mathbf{v}_2\ \mathbf{v}_3]$:
> $$
> \begin{bmatrix}1&0&1\\2&1&0\\0&1&-2\end{bmatrix}\xrightarrow{R_2-2R_1}\begin{bmatrix}1&0&1\\0&1&-2\\0&1&-2\end{bmatrix}\xrightarrow{R_3-R_2}\begin{bmatrix}1&0&1\\0&1&-2\\0&0&0\end{bmatrix}.
> $$
> Column 3 has no pivot, so the set is **dependent**. Back-substitute with $c_3=1$: $c_2=2$, $c_1=-1$. Dependence relation:
> $$
> -\mathbf{v}_1+2\mathbf{v}_2+\mathbf{v}_3=\mathbf{0}\quad\Longleftrightarrow\quad\mathbf{v}_3=\mathbf{v}_1-2\mathbf{v}_2 .
> $$
> Check: $(1,2,0)-2(0,1,1)=(1,0,-2)$ ✓. (Also $\det\mathbf{V}=0$.)
>
> (b) Are $\mathbf{w}_1=(1,1,0)^\mathsf{T}$, $\mathbf{w}_2=(1,0,1)^\mathsf{T}$, $\mathbf{w}_3=(0,1,1)^\mathsf{T}$ independent?
> $$
> \det\begin{bmatrix}1&1&0\\1&0&1\\0&1&1\end{bmatrix}=1(0-1)-1(1-0)+0=-2\neq0 .
> $$
> **Independent.** Three independent vectors in $\mathbb{R}^3$ — they will turn out to be a basis.

### Independence of functions

For vectors in a function space, the coefficients $c_i$ are constants and the relation $c_1f_1(t)+\dots+c_kf_k(t)=0$ must hold **for every $t$** in the interval. Two practical tests:

**1. Evaluate at sample points.** If the functions are dependent, then for any points $t_1,\dots,t_k$ the matrix $[f_j(t_i)]$ is singular. So finding *one* set of points where $\det[f_j(t_i)]\neq0$ proves independence.

**2. The Wronskian.** For $k$ functions with $k-1$ derivatives, define

$$
W(f_1,\dots,f_k)(t)=\det\begin{bmatrix}f_1&f_2&\cdots&f_k\\ f_1'&f_2'&\cdots&f_k'\\ \vdots&&&\vdots\\ f_1^{(k-1)}&f_2^{(k-1)}&\cdots&f_k^{(k-1)}\end{bmatrix}.
$$

> [!theorem] Wronskian test
> If $W(f_1,\dots,f_k)(t_0)\neq0$ at some point $t_0$ of the interval, then $f_1,\dots,f_k$ are linearly independent there.
>
> If $f_1,\dots,f_k$ are solutions of the same $k$th-order linear homogeneous ODE with continuous coefficients, then the converse also holds: they are independent iff $W\neq0$, and in that case $W(t)\neq0$ for **every** $t$.

*Proof of the first part.* Differentiate a dependence relation $\sum c_if_i(t)=0$ repeatedly ($k-1$ times); at $t_0$ this gives the square system $\mathbf{W}(t_0)\mathbf{c}=\mathbf{0}$. If $\det\mathbf{W}(t_0)\neq0$ then $\mathbf{c}=\mathbf{0}$.

> [!warning] Common Pitfall — $W\equiv0$ does not prove dependence in general
> The functions $t^2$ and $t\lvert t\rvert$ on $[-1,1]$ have Wronskian identically zero, yet they are independent (on $t>0$ they agree, on $t<0$ they are negatives of each other, so no single constant ratio works). The converse of the Wronskian test needs the extra hypothesis that the functions solve a common linear ODE.

> [!example] Worked Example 9.8 — Exponentials and trigonometric functions
> (a) Show $e^{t},e^{2t},e^{3t}$ are independent on $\mathbb{R}$.
> $$
> W=\det\begin{bmatrix}e^t&e^{2t}&e^{3t}\\ e^t&2e^{2t}&3e^{3t}\\ e^t&4e^{2t}&9e^{3t}\end{bmatrix}=e^{6t}\det\begin{bmatrix}1&1&1\\1&2&3\\1&4&9\end{bmatrix}=2e^{6t}\neq0 .
> $$
> The $3\times3$ determinant is a Vandermonde determinant $(2-1)(3-1)(3-2)=2$. In general $e^{s_1t},\dots,e^{s_kt}$ are independent whenever the $s_i$ are distinct — the reason distinct characteristic roots give independent natural-response modes ([Ch. 12](../part-4-differential-equations/12-linear-differential-equations.md)).
>
> (b) Are $1,\cos^2t,\sin^2t$ independent? No: $1\cdot1-1\cdot\cos^2t-1\cdot\sin^2t=0$ for all $t$. Likewise $\{1,\cos^2t,\cos2t\}$ is dependent because $\cos2t=2\cos^2t-1$.
>
> (c) Are $\cos\omega t$ and $\sin\omega t$ independent? $W=\cos\omega t\cdot\omega\cos\omega t-\sin\omega t\cdot(-\omega\sin\omega t)=\omega\neq0$ for $\omega\neq0$. Yes — which is why every sinusoid at frequency $\omega$ needs exactly two numbers (amplitude and phase, or $a$ and $b$).

> [!example] Worked Example 9.9 — Independence in $\mathcal{P}_2$
> Are $p_1=1+t$, $p_2=1-t$, $p_3=t^2$ independent in $\mathcal{P}_2$?
>
> Suppose $c_1(1+t)+c_2(1-t)+c_3t^2=0$ for all $t$. Collect powers of $t$:
> $$
> (c_1+c_2)+(c_1-c_2)t+c_3t^2=0\ \Rightarrow\ c_1+c_2=0,\quad c_1-c_2=0,\quad c_3=0 .
> $$
> Hence $c_1=c_2=c_3=0$: **independent**. The method — matching coefficients — converts a question about functions into a question about the coefficient vectors $(1,1,0)$, $(1,-1,0)$, $(0,0,1)$ in $\mathbb{R}^3$. This translation is made systematic by *coordinates* in Section 9.4.

### Geometric picture

In $\mathbb{R}^2$ two vectors are independent iff they are not collinear; in $\mathbb{R}^3$ three vectors are independent iff they are not coplanar. The determinant measures the area or volume of the parallelogram/parallelepiped they span, so "nearly dependent" vectors give a small determinant relative to the product of their lengths. The graph below shows two nearly parallel vectors; any vector in the plane can still be written as a combination of them, but the coefficients become large and very sensitive to noise.

```desmos-graph
left=-1; right=6; bottom=-1; top=4
---
y=0.5x|blue|0<x<4
y=0.6x|red|0<x<4
(4,2)|label:v1=(4,2)
(4,2.4)|label:v2=(4,2.4)
(1,3)|label:b=(1,3)
```
*Nearly parallel vectors $\mathbf{v}_1=(4,2)$ and $\mathbf{v}_2=(4,2.4)$ are independent, but $\mathbf{b}=(1,3)=-6\,\mathbf{v}_1+6.25\,\mathbf{v}_2$ needs large, cancelling coefficients — the hallmark of ill-conditioning.*

> [!tip] Engineering Note — near-dependence in measurements
> In practice the question is rarely "independent or not?" but "how independent?". Two antenna-array steering vectors at nearly the same angle, two sensor channels that are almost proportional, or regressor columns $t$ and $t^2$ sampled on $[100,101]$ are all *technically* independent but numerically almost dependent. The relevant measure is the condition number of the matrix of columns, introduced in [Ch. 8](08-matrices-and-linear-systems.md) and made precise with singular values in [Ch. 10](10-eigenvalues-and-eigenvectors.md). It reappears in least squares (Section 9.11).

## 9.4 Basis and dimension

### Motivation

A basis is a minimal set of "building blocks" for a space: enough vectors to span it, with no redundancy. Once a basis is chosen, every vector is described by a unique list of numbers — its **coordinates** — and abstract vectors (polynomials, sinusoids, signals) can be handled by matrix computations. The number of vectors in a basis, the **dimension**, is the number of independent degrees of freedom: the number of initial conditions a circuit needs, the number of loop currents, the number of samples needed to pin down a band-limited signal over an interval.

### Definition

> [!definition] Basis
> A set $\mathcal{B}=\{\mathbf{b}_1,\dots,\mathbf{b}_n\}$ in $V$ is a **basis** of $V$ if
> 1. $\mathcal{B}$ is linearly independent, and
> 2. $\mathcal{B}$ spans $V$.

> [!theorem] Unique representation
> $\mathcal{B}$ is a basis of $V$ iff every $\mathbf{v}\in V$ can be written in exactly one way as
> $$
> \mathbf{v}=c_1\mathbf{b}_1+\dots+c_n\mathbf{b}_n .
> $$
> The numbers $c_1,\dots,c_n$ are the **coordinates** of $\mathbf{v}$ relative to $\mathcal{B}$, written $[\mathbf{v}]_\mathcal{B}=(c_1,\dots,c_n)^\mathsf{T}$.

Spanning gives existence; independence gives uniqueness. Note that a basis is an *ordered* list: swapping $\mathbf{b}_1$ and $\mathbf{b}_2$ swaps the first two coordinates.

**Standard bases.**

| Space | Standard basis | Dimension |
|---|---|---|
| $\mathbb{R}^n$, $\mathbb{C}^n$ | $\mathbf{e}_1=(1,0,\dots,0)^\mathsf{T},\dots,\mathbf{e}_n$ | $n$ |
| $\mathcal{P}_n$ | $1,t,t^2,\dots,t^n$ | $n+1$ |
| $\mathbb{R}^{m\times n}$ | $\mathbf{E}_{ij}$ (1 in position $(i,j)$, zeros elsewhere) | $mn$ |
| $\mathcal{S}_\omega$ (sinusoids at $\omega$) | $\cos\omega t,\ \sin\omega t$ | 2 |
| trig. polynomials of order $N$ | $1,\cos k\omega_0t,\sin k\omega_0t$ ($k=1,\dots,N$) | $2N+1$ |
| solutions of an $n$th-order linear homogeneous ODE | any $n$ independent solutions | $n$ |

### Dimension

The definition of dimension rests on a counting theorem.

> [!theorem] Exchange (replacement) lemma
> If $V$ is spanned by $m$ vectors, then any linearly independent set in $V$ has at most $m$ vectors.

*Proof sketch.* Let $\mathbf{w}_1,\dots,\mathbf{w}_m$ span $V$ and let $\mathbf{u}_1,\dots,\mathbf{u}_k$ be independent with $k>m$. Each $\mathbf{u}_j=\sum_i a_{ij}\mathbf{w}_i$, i.e. $[\mathbf{u}_1\cdots\mathbf{u}_k]=[\mathbf{w}_1\cdots\mathbf{w}_m]\mathbf{A}$ with $\mathbf{A}$ of size $m\times k$. Since $k>m$, the homogeneous system $\mathbf{A}\mathbf{c}=\mathbf{0}$ has a nonzero solution. Then $\sum_jc_j\mathbf{u}_j=[\mathbf{w}_1\cdots\mathbf{w}_m]\mathbf{A}\mathbf{c}=\mathbf{0}$ — a nontrivial dependence among the $\mathbf{u}_j$, contradiction.

> [!theorem] Invariance of dimension
> All bases of a vector space $V$ have the same number of elements. This number is the **dimension** $\dim V$. By convention $\dim\{\mathbf{0}\}=0$. If no finite set spans $V$, $V$ is **infinite-dimensional**.

*Proof.* If $\mathcal{B}_1$ has $n_1$ elements and $\mathcal{B}_2$ has $n_2$, apply the lemma both ways: $\mathcal{B}_2$ independent and $\mathcal{B}_1$ spanning gives $n_2\le n_1$; the reverse gives $n_1\le n_2$.

The function spaces $C[a,b]$ and $L^2(a,b)$ are infinite-dimensional: for every $n$ the functions $1,t,\dots,t^n$ are independent in them, so no finite set can span. This is why a Fourier series needs infinitely many terms in general, and why "convergence" becomes an issue there ([Ch. 15](../part-5-transforms/15-fourier-series.md)) that never arises in $\mathbb{R}^n$.

> [!theorem] Consequences in an $n$-dimensional space
> Let $\dim V=n$.
> 1. Any $n+1$ or more vectors are dependent; any $n-1$ or fewer cannot span.
> 2. Any $n$ independent vectors form a basis; any $n$ spanning vectors form a basis. (So for exactly $n$ vectors, check *one* of the two conditions.)
> 3. Every independent set can be **extended** to a basis; every spanning set can be **pruned** to a basis.
> 4. If $W$ is a subspace of $V$ then $\dim W\le n$, with equality iff $W=V$.

Property 2 is a great labour-saver: three independent vectors in $\mathbb{R}^3$, such as $\mathbf{w}_1,\mathbf{w}_2,\mathbf{w}_3$ of Worked Example 9.7(b), automatically form a basis.

For two subspaces $U,W$ of a finite-dimensional space there is also a counting formula,

$$
\dim(U+W)=\dim U+\dim W-\dim(U\cap W),
$$

the linear-algebra analogue of $\lvert A\cup B\rvert=\lvert A\rvert+\lvert B\rvert-\lvert A\cap B\rvert$. For the two planes of Worked Example 9.6: $\dim(U+W)=2+2-1=3$, so together they span $\mathbb{R}^3$.

> [!example] Worked Example 9.10 — Basis and dimension of a plane
> Find a basis and the dimension of $W=\{(x,y,z): x+2y-z=0\}$.
>
> The single equation has one pivot ($x$) and two free variables ($y,z$). Write $x=-2y+z$:
> $$
> \begin{bmatrix}x\\y\\z\end{bmatrix}=y\begin{bmatrix}-2\\1\\0\end{bmatrix}+z\begin{bmatrix}1\\0\\1\end{bmatrix}.
> $$
> So $W=\operatorname{span}\{(-2,1,0)^\mathsf{T},(1,0,1)^\mathsf{T}\}$. The two vectors are independent (look at the 2nd and 3rd components: each has a 1 where the other has a 0). Hence they form a basis and $\dim W=2$.
>
> **General rule.** The solution space of $\mathbf{A}\mathbf{x}=\mathbf{0}$ has dimension equal to the number of free variables, and the "special solutions" obtained by setting one free variable to 1 and the others to 0 form a basis.

> [!example] Worked Example 9.11 — Coordinates in a non-standard polynomial basis
> Show that $\mathcal{B}=\{1,\ 1+t,\ 1+t+t^2\}$ is a basis of $\mathcal{P}_2$ and find the coordinates of $p(t)=3+2t-t^2$.
>
> **Basis.** Three vectors in a 3-dimensional space: we only need independence. Their coefficient vectors relative to $\{1,t,t^2\}$ are the columns of an upper-triangular matrix with unit diagonal,
> $$
> \mathbf{P}=\begin{bmatrix}1&1&1\\0&1&1\\0&0&1\end{bmatrix},\qquad\det\mathbf{P}=1\neq0 .
> $$
> **Coordinates.** Solve $a\cdot1+b(1+t)+c(1+t+t^2)=3+2t-t^2$. Matching coefficients:
> $$
> t^2:\ c=-1;\qquad t:\ b+c=2\Rightarrow b=3;\qquad 1:\ a+b+c=3\Rightarrow a=1 .
> $$
> So $[p]_\mathcal{B}=(1,3,-1)^\mathsf{T}$. Check: $1+3(1+t)-(1+t+t^2)=3+2t-t^2$ ✓.

### Change of basis

Let $\mathcal{B}=\{\mathbf{b}_1,\dots,\mathbf{b}_n\}$ be a basis of $\mathbb{F}^n$ and $\mathbf{P}_\mathcal{B}=[\mathbf{b}_1\ \cdots\ \mathbf{b}_n]$. The equation $\mathbf{x}=c_1\mathbf{b}_1+\dots+c_n\mathbf{b}_n$ reads

$$
\mathbf{x}=\mathbf{P}_\mathcal{B}\,[\mathbf{x}]_\mathcal{B},\qquad [\mathbf{x}]_\mathcal{B}=\mathbf{P}_\mathcal{B}^{-1}\mathbf{x}.
$$

$\mathbf{P}_\mathcal{B}$ converts $\mathcal{B}$-coordinates to standard coordinates; its inverse goes the other way. Between two non-standard bases $\mathcal{B}$ and $\mathcal{C}$,

$$
[\mathbf{x}]_\mathcal{C}=\mathbf{P}_{\mathcal{C}\leftarrow\mathcal{B}}[\mathbf{x}]_\mathcal{B},\qquad \mathbf{P}_{\mathcal{C}\leftarrow\mathcal{B}}=\mathbf{P}_\mathcal{C}^{-1}\mathbf{P}_\mathcal{B},
$$

whose $j$th column is $[\mathbf{b}_j]_\mathcal{C}$.

> [!example] Worked Example 9.12 — Common-mode / differential-mode coordinates
> Two-wire signals are often described not by the individual wire voltages $(v_1,v_2)$ but by **common-mode** and **differential-mode** components. Use the basis $\mathbf{b}_1=(1,1)^\mathsf{T}$ (common mode) and $\mathbf{b}_2=(1,-1)^\mathsf{T}$ (differential mode). Find the coordinates of $\mathbf{v}=(3,1)^\mathsf{T}$ V.
>
> $$
> \mathbf{P}_\mathcal{B}=\begin{bmatrix}1&1\\1&-1\end{bmatrix},\qquad \mathbf{P}_\mathcal{B}^{-1}=\frac{1}{-2}\begin{bmatrix}-1&-1\\-1&1\end{bmatrix}=\frac12\begin{bmatrix}1&1\\1&-1\end{bmatrix}.
> $$
> $$
> [\mathbf{v}]_\mathcal{B}=\frac12\begin{bmatrix}1&1\\1&-1\end{bmatrix}\begin{bmatrix}3\\1\end{bmatrix}=\begin{bmatrix}2\\1\end{bmatrix}.
> $$
> So $\mathbf{v}=2(1,1)+1(1,-1)$: a common-mode level $v_{cm}=(v_1+v_2)/2=2$ V and a differential half-amplitude $(v_1-v_2)/2=1$ V (the full differential voltage is $v_1-v_2=2$ V). A differential amplifier with ideal common-mode rejection responds only to the second coordinate.

```desmos-graph
left=-1; right=4; bottom=-2; top=3
---
y=x|blue|dashed
y=-x|red|dashed
(3,1)|label:v=(3,1)
(2,2)|label:2b1
(1,-1)|label:b2
y=-(x-3)+1|grey|dotted|2<x<3
y=(x-3)+1|grey|dotted|1<x<3
```
*The vector $\mathbf{v}=(3,1)$ resolved along the common-mode axis (blue, $y=x$) and the differential-mode axis (red, $y=-x$): $\mathbf{v}=2\mathbf{b}_1+1\mathbf{b}_2$; the dotted lines complete the parallelogram.*

> [!example] Worked Example 9.13 — Dimension of familiar spaces
> (a) The $2\times2$ real symmetric matrices $\begin{bmatrix}a&b\\b&d\end{bmatrix}=a\begin{bmatrix}1&0\\0&0\end{bmatrix}+b\begin{bmatrix}0&1\\1&0\end{bmatrix}+d\begin{bmatrix}0&0\\0&1\end{bmatrix}$ form a space of dimension **3**.
>
> (b) $W_c=\{p\in\mathcal{P}_3: p(1)=0\}$ from Worked Example 9.4. By the factor theorem $p(t)=(t-1)q(t)$ with $q\in\mathcal{P}_2$, so a basis is $\{t-1,\ t(t-1),\ t^2(t-1)\}$ and $\dim W_c=3=\dim\mathcal{P}_3-1$: one linear constraint removes one dimension.
>
> (c) Solutions of $y''+3y'+2y=0$: basis $\{e^{-t},e^{-2t}\}$ (independent by the Wronskian $W=e^{-t}(-2e^{-2t})-e^{-2t}(-e^{-t})=-e^{-3t}\neq0$), dimension **2** — matching the two initial conditions $y(0)$, $y'(0)$ a second-order circuit needs.

> [!warning] Common Pitfall — dimension is not the number of components
> The plane $x+2y-z=0$ consists of vectors with three components, yet its dimension is 2. Dimension counts *degrees of freedom*, not entries. Likewise $\mathcal{P}_2$ has dimension 3 (three coefficients), not 2.

> [!tip] Engineering Note — counting unknowns in circuit analysis
> For a connected circuit graph with $n$ nodes and $b$ branches, the space of KCL-compatible branch-current vectors has dimension $b-n+1$ and the space of branch-voltage vectors that obey KVL has dimension $n-1$. These are exactly the number of mesh (loop) currents and of node voltages (with one reference node) that mesh and nodal analysis use. Choosing "mesh currents" or "node voltages" is choosing a **basis** for one of these subspaces. Section 9.6 proves the counts with the rank–nullity theorem.

## 9.5 Linear transformations

### Motivation

An amplifier, a filter, a differentiator, a sampler, a coordinate rotation, a Fourier transform — each takes an input from one vector space and produces an output in another, and each obeys **superposition**: the response to a weighted sum of inputs is the same weighted sum of responses. Mappings with this property are **linear transformations** (also called linear maps or linear operators). Superposition is the single most exploited property in electrical engineering: it lets us analyse a circuit one source at a time, decompose a signal into sinusoids and treat each separately, and represent the whole system by a matrix or a transfer function.

### Definition

> [!definition] Linear transformation
> Let $V$ and $W$ be vector spaces over the same field $\mathbb{F}$. A map $T:V\to W$ is **linear** if for all $\mathbf{u},\mathbf{v}\in V$ and $c\in\mathbb{F}$
> 1. $T(\mathbf{u}+\mathbf{v})=T(\mathbf{u})+T(\mathbf{v})$ (additivity),
> 2. $T(c\mathbf{u})=c\,T(\mathbf{u})$ (homogeneity).
>
> Equivalently, $T(a\mathbf{u}+b\mathbf{v})=aT(\mathbf{u})+bT(\mathbf{v})$ — the **superposition principle**. $V$ is the **domain**, $W$ the **codomain**. A linear map $V\to V$ is often called a linear **operator** on $V$.

Immediate consequences:

$$
T(\mathbf{0})=\mathbf{0},\qquad T\Big(\sum_i c_i\mathbf{v}_i\Big)=\sum_i c_iT(\mathbf{v}_i).
$$

The first gives a one-line non-linearity test: **if $T(\mathbf{0})\neq\mathbf{0}$, $T$ is not linear.** The second says a linear map is completely determined by what it does to a basis: if $\{\mathbf{b}_j\}$ is a basis of $V$ and we know $T(\mathbf{b}_j)$ for each $j$, we know $T$ everywhere.

### Examples and non-examples

| Map | Linear? | Reason |
|---|---|---|
| $T(\mathbf{x})=\mathbf{A}\mathbf{x}$, $\mathbf{A}\in\mathbb{F}^{m\times n}$ | yes | matrix algebra |
| rotation of $\mathbb{R}^2$ by angle $\theta$ | yes | the rotated parallelogram of $\mathbf{u},\mathbf{v}$ is the parallelogram of the rotated vectors |
| $D:\mathcal{P}_n\to\mathcal{P}_{n-1}$, $D(p)=p'$ | yes | $(ap+bq)'=ap'+bq'$ |
| $J:C[0,T]\to C[0,T]$, $J(f)(t)=\int_0^tf(\tau)\,d\tau$ | yes | linearity of the integral |
| time shift $(S_\tau x)(t)=x(t-\tau)$ | yes | shifts commute with sums |
| convolution with a fixed $h$: $x\mapsto h*x$ | yes | integral of a linear expression |
| evaluation $x\mapsto x(t_0)$, sampling $x\mapsto(x(0),x(T_s),\dots)$ | yes | pointwise operations |
| $T(x)=ax+b$ with $b\neq0$ ("linear" in school algebra) | **no** | $T(0)=b\neq0$; it is **affine** |
| squaring $v\mapsto v^2$, rectifier $v\mapsto\lvert v\rvert$, clipper | **no** | $(2v)^2\neq2v^2$ |
| $\mathbf{x}\mapsto\lVert\mathbf{x}\rVert$ | **no** | $\lVert-\mathbf{x}\rVert\neq-\lVert\mathbf{x}\rVert$ |
| complex conjugation $z\mapsto z^*$ on $\mathbb{C}$ as a *complex* space | **no** | $(jz)^*=-jz^*\neq jz^*$ (but it is linear on $\mathbb{C}$ viewed as the real space $\mathbb{R}^2$) |

> [!example] Worked Example 9.14 — Testing linearity
> (a) $T:\mathbb{R}^2\to\mathbb{R}^2$, $T(x,y)=(x+2y,\ 3x)$. Write $T(\mathbf{x})=\begin{bmatrix}1&2\\3&0\end{bmatrix}\mathbf{x}$: linear.
>
> (b) $T(x,y)=(x+1,\ y)$. $T(0,0)=(1,0)\neq\mathbf{0}$: **not linear**.
>
> (c) $T(x,y)=(xy,\ x)$. $T(2(1,1))=T(2,2)=(4,2)$ but $2T(1,1)=(2,2)$: **not linear**.
>
> (d) The RC low-pass filter initially at rest, $RC\,\dot y+y=x$, $y(0)=0$, defines $y=T(x)$ with
> $$
> y(t)=\frac{1}{RC}\int_0^t e^{-(t-\tau)/RC}\,x(\tau)\,d\tau .
> $$
> For inputs $x_1,x_2$ and scalars $a,b$, the integral of $a x_1+bx_2$ splits into $a$ times the first integral plus $b$ times the second: $T(ax_1+bx_2)=aT(x_1)+bT(x_2)$. **Linear.**
>
> (e) The same filter with a nonzero initial voltage $y(0)=V_0$: now $T(0)=V_0e^{-t/RC}\neq0$. **Not linear** in the strict sense — it is affine (zero-input response plus a linear zero-state response). This is why "linear system" in circuit theory always means *zero initial state* when superposition is invoked.

### The matrix of a linear map

> [!theorem] Every linear map $\mathbb{F}^n\to\mathbb{F}^m$ is a matrix
> If $T:\mathbb{F}^n\to\mathbb{F}^m$ is linear, then $T(\mathbf{x})=\mathbf{A}\mathbf{x}$ for the unique $m\times n$ matrix
> $$
> \mathbf{A}=\begin{bmatrix}T(\mathbf{e}_1)&T(\mathbf{e}_2)&\cdots&T(\mathbf{e}_n)\end{bmatrix},
> $$
> called the **standard matrix** of $T$.

*Proof.* $\mathbf{x}=\sum_jx_j\mathbf{e}_j$, so $T(\mathbf{x})=\sum_jx_jT(\mathbf{e}_j)$, which is the matrix–vector product of $[T(\mathbf{e}_1)\cdots T(\mathbf{e}_n)]$ with $\mathbf{x}$.

**Recipe:** *the $j$th column of the matrix is the image of the $j$th basis vector.*

More generally, for $T:V\to W$ with bases $\mathcal{B}=\{\mathbf{b}_1,\dots,\mathbf{b}_n\}$ of $V$ and $\mathcal{C}$ of $W$, the **matrix of $T$ relative to $\mathcal{B}$ and $\mathcal{C}$** is

$$
[T]_{\mathcal{C}\leftarrow\mathcal{B}}=\begin{bmatrix}[T(\mathbf{b}_1)]_\mathcal{C}&\cdots&[T(\mathbf{b}_n)]_\mathcal{C}\end{bmatrix},\qquad [T(\mathbf{v})]_\mathcal{C}=[T]_{\mathcal{C}\leftarrow\mathcal{B}}\,[\mathbf{v}]_\mathcal{B}.
$$

So *any* linear map between finite-dimensional spaces — differentiation of polynomials, a filter acting on a finite set of sinusoids — becomes a matrix once bases are chosen.

> [!example] Worked Example 9.15 — Rotation and projection in the plane
> (a) **Rotation by $\theta$.** $\mathbf{e}_1=(1,0)$ rotates to $(\cos\theta,\sin\theta)$ and $\mathbf{e}_2=(0,1)$ to $(-\sin\theta,\cos\theta)$, so
> $$
> \mathbf{R}_\theta=\begin{bmatrix}\cos\theta&-\sin\theta\\ \sin\theta&\cos\theta\end{bmatrix}.
> $$
> This is the matrix form of multiplying a phasor by $e^{j\theta}$: $(x+jy)e^{j\theta}=(x\cos\theta-y\sin\theta)+j(x\sin\theta+y\cos\theta)$.
>
> (b) **Orthogonal projection onto the line at angle $\theta$** through the origin, with unit direction $\hat{\mathbf{u}}=(\cos\theta,\sin\theta)$. Section 9.9 shows $P(\mathbf{x})=(\hat{\mathbf{u}}\cdot\mathbf{x})\hat{\mathbf{u}}$. Images of the basis vectors: $P(\mathbf{e}_1)=\cos\theta\,\hat{\mathbf{u}}$, $P(\mathbf{e}_2)=\sin\theta\,\hat{\mathbf{u}}$, so
> $$
> \mathbf{P}_\theta=\begin{bmatrix}\cos^2\theta&\cos\theta\sin\theta\\ \cos\theta\sin\theta&\sin^2\theta\end{bmatrix}.
> $$
> For $\theta=30^\circ$: $\mathbf{P}=\begin{bmatrix}3/4&\sqrt3/4\\ \sqrt3/4&1/4\end{bmatrix}$ and $\mathbf{P}(2,0)^\mathsf{T}=(1.5,\ 0.866)^\mathsf{T}$. Check: this point lies on $y=x\tan30^\circ=0.5774x$ ($1.5\times0.5774=0.866$ ✓), and the residual $(0.5,-0.866)$ is perpendicular to $\hat{\mathbf{u}}$: $0.5\cdot0.866-0.866\cdot0.5=0$ ✓.

```desmos-graph
left=-0.5; right=3; bottom=-1; top=2
---
y=0.57735x|blue
(2,0)|label:x=(2,0)
(1.5,0.866)|label:Px
y=-1.732(x-2)|red|dashed|1.5<x<2
x^2+y^2=4|grey|dotted
(1.732,1)|label:R30 x
```
*The vector $\mathbf{x}=(2,0)$, its rotation by $30^\circ$ (on the dotted circle of radius 2) and its orthogonal projection onto the $30^\circ$ line; the dashed segment is perpendicular to the line.*

> [!example] Worked Example 9.16 — The differentiation matrix
> Find the matrix of $D:\mathcal{P}_3\to\mathcal{P}_2$, $D(p)=p'$, relative to the monomial bases $\{1,t,t^2,t^3\}$ and $\{1,t,t^2\}$, and use it to differentiate $p(t)=2-t+4t^2+t^3$.
>
> Images of the basis vectors: $D(1)=0$, $D(t)=1$, $D(t^2)=2t$, $D(t^3)=3t^2$. In coordinates these are $(0,0,0)$, $(1,0,0)$, $(0,2,0)$, $(0,0,3)$, so
> $$
> [D]=\begin{bmatrix}0&1&0&0\\0&0&2&0\\0&0&0&3\end{bmatrix}.
> $$
> With $[p]=(2,-1,4,1)^\mathsf{T}$:
> $$
> [D][p]=\begin{bmatrix}0\cdot2+1\cdot(-1)\\ 2\cdot4\\ 3\cdot1\end{bmatrix}=\begin{bmatrix}-1\\8\\3\end{bmatrix}\ \Rightarrow\ p'(t)=-1+8t+3t^2 .
> $$
> Direct check: $\frac{d}{dt}(2-t+4t^2+t^3)=-1+8t+3t^2$ ✓. The matrix has a zero first column: constants are sent to zero, which is the null space of $D$ (Section 9.6).

### Composition, inverses and isomorphism

If $S:U\to V$ and $T:V\to W$ are linear, so is the composition $T\circ S$, and in coordinates

$$
[T\circ S]=[T]\,[S].
$$

Matrix multiplication is *defined* so that this holds — it is the reason the row-times-column rule looks the way it does. A cascade of linear stages (pre-amplifier, filter, ADC) is represented by the product of their matrices, applied **right to left**.

A linear map $T:V\to W$ is

- **one-to-one** (injective) if $T(\mathbf{u})=T(\mathbf{v})\Rightarrow\mathbf{u}=\mathbf{v}$;
- **onto** (surjective) if every $\mathbf{w}\in W$ equals $T(\mathbf{v})$ for some $\mathbf{v}$;
- an **isomorphism** if it is both; then $T^{-1}$ exists and is linear.

> [!theorem] Coordinates are an isomorphism
> If $\dim V=n$ and $\mathcal{B}$ is a basis, the coordinate map $\mathbf{v}\mapsto[\mathbf{v}]_\mathcal{B}$ is an isomorphism $V\to\mathbb{F}^n$. Consequently any two vector spaces of the same finite dimension over the same field are isomorphic.

This theorem is what licenses replacing $\mathcal{P}_2$ by $\mathbb{R}^3$, sinusoids at $\omega$ by $\mathbb{R}^2$ (or by phasors in $\mathbb{C}$), and so on: every finite-dimensional calculation can be done with columns of numbers.

### Change of basis for a linear operator

Let $T:\mathbb{F}^n\to\mathbb{F}^n$ have standard matrix $\mathbf{A}$, and let $\mathcal{B}$ be another basis with $\mathbf{P}=\mathbf{P}_\mathcal{B}$. Since $\mathbf{x}=\mathbf{P}[\mathbf{x}]_\mathcal{B}$,

$$
[T(\mathbf{x})]_\mathcal{B}=\mathbf{P}^{-1}\mathbf{A}\mathbf{x}=\mathbf{P}^{-1}\mathbf{A}\mathbf{P}\,[\mathbf{x}]_\mathcal{B}\quad\Longrightarrow\quad [T]_\mathcal{B}=\mathbf{P}^{-1}\mathbf{A}\mathbf{P}.
$$

Matrices related by $\mathbf{B}=\mathbf{P}^{-1}\mathbf{A}\mathbf{P}$ are **similar**: they represent the *same* operator in different coordinates, so they share every basis-independent property — determinant, trace, rank and (Chapter 10) eigenvalues.

> [!example] Worked Example 9.17 — A basis that makes an operator simple
> Let $\mathbf{A}=\begin{bmatrix}1&1\\0&2\end{bmatrix}$ and $\mathcal{B}=\{(1,0)^\mathsf{T},(1,1)^\mathsf{T}\}$. Find $[T]_\mathcal{B}$.
>
> $$
> \mathbf{P}=\begin{bmatrix}1&1\\0&1\end{bmatrix},\quad\mathbf{P}^{-1}=\begin{bmatrix}1&-1\\0&1\end{bmatrix},\quad\mathbf{A}\mathbf{P}=\begin{bmatrix}1&2\\0&2\end{bmatrix},
> $$
> $$
> \mathbf{P}^{-1}\mathbf{A}\mathbf{P}=\begin{bmatrix}1&-1\\0&1\end{bmatrix}\begin{bmatrix}1&2\\0&2\end{bmatrix}=\begin{bmatrix}1&0\\0&2\end{bmatrix}.
> $$
> In the new basis the operator simply scales the first coordinate by 1 and the second by 2: $\mathbf{A}(1,0)=(1,0)$ and $\mathbf{A}(1,1)=(2,2)$. These special basis vectors are **eigenvectors**; finding them systematically is the subject of [Ch. 10](10-eigenvalues-and-eigenvectors.md). Check invariants: $\det\mathbf{A}=2=\det\operatorname{diag}(1,2)$ and $\operatorname{tr}\mathbf{A}=3=1+2$ ✓.

> [!warning] Common Pitfall — "linear" versus "affine"
> The straight-line function $y=mx+c$ is called linear in school algebra and in device data sheets, but with $c\neq0$ it is **affine**, not linear: doubling the input does not double the output. A biased transistor stage is affine about its operating point; small-signal analysis subtracts the operating point so that the remaining map *is* linear. Always remove offsets (DC bias, initial conditions) before invoking superposition.

> [!tip] Engineering Note — LTI systems are linear operators on signal spaces
> A linear time-invariant (LTI) system is a linear operator $T$ on a signal space that commutes with every time shift: $T\circ S_\tau=S_\tau\circ T$. One can show that every such operator is a convolution $y=h*x$ with its impulse response $h$. Complex exponentials $e^{j\omega t}$ satisfy $T(e^{j\omega t})=H(j\omega)e^{j\omega t}$ — they are the "eigenvectors" of every LTI operator, the infinite-dimensional analogue of Worked Example 9.17. Expanding inputs in that basis is exactly Fourier analysis ([Ch. 15](../part-5-transforms/15-fourier-series.md), [Ch. 16](../part-5-transforms/16-fourier-transforms.md)).

## 9.6 Null space and column space

### Motivation

For a linear map $T(\mathbf{x})=\mathbf{A}\mathbf{x}$ two questions dominate:

1. **Existence:** for which outputs $\mathbf{b}$ does $\mathbf{A}\mathbf{x}=\mathbf{b}$ have a solution?
2. **Uniqueness:** if a solution exists, how much freedom is left in it?

The answer to the first is a subspace of the output space, the **column space**; the answer to the second is a subspace of the input space, the **null space**. Their dimensions are tied together by the **rank–nullity theorem**, which in circuit theory counts independent node voltages and loop currents, in control theory counts uncontrollable modes, and in signal processing counts the components of a signal a measurement system cannot see.

### Definitions

> [!definition] The four fundamental subspaces
> For $\mathbf{A}\in\mathbb{R}^{m\times n}$:
> - **Null space (kernel):** $\operatorname{null}(\mathbf{A})=\{\mathbf{x}\in\mathbb{R}^n:\mathbf{A}\mathbf{x}=\mathbf{0}\}\subseteq\mathbb{R}^n$.
> - **Column space (range):** $\operatorname{col}(\mathbf{A})=\{\mathbf{A}\mathbf{x}:\mathbf{x}\in\mathbb{R}^n\}=\operatorname{span}\{\text{columns of }\mathbf{A}\}\subseteq\mathbb{R}^m$.
> - **Row space:** $\operatorname{row}(\mathbf{A})=\operatorname{col}(\mathbf{A}^\mathsf{T})\subseteq\mathbb{R}^n$.
> - **Left null space:** $\operatorname{null}(\mathbf{A}^\mathsf{T})=\{\mathbf{y}\in\mathbb{R}^m:\mathbf{A}^\mathsf{T}\mathbf{y}=\mathbf{0}\}\subseteq\mathbb{R}^m$.
>
> For a general linear map $T:V\to W$, $\ker T=\{\mathbf{v}:T(\mathbf{v})=\mathbf{0}\}$ and $\operatorname{range}T=\{T(\mathbf{v})\}$. The dimension of the column space is the **rank** $r$; the dimension of the null space is the **nullity**. (Over $\mathbb{C}$ replace $\mathsf{T}$ by $\mathsf{H}$.)

Both are subspaces: $\mathbf{A}\mathbf{0}=\mathbf{0}$; if $\mathbf{A}\mathbf{x}=\mathbf{A}\mathbf{y}=\mathbf{0}$ then $\mathbf{A}(a\mathbf{x}+b\mathbf{y})=\mathbf{0}$; and a span is always a subspace.

The column space viewpoint comes from reading $\mathbf{A}\mathbf{x}$ **column by column**:

$$
\mathbf{A}\mathbf{x}=x_1\mathbf{a}_1+x_2\mathbf{a}_2+\dots+x_n\mathbf{a}_n .
$$

So $\mathbf{A}\mathbf{x}=\mathbf{b}$ asks: *which combination of the columns produces $\mathbf{b}$?*

> [!theorem] Existence and structure of solutions
> 1. $\mathbf{A}\mathbf{x}=\mathbf{b}$ is solvable iff $\mathbf{b}\in\operatorname{col}(\mathbf{A})$.
> 2. If $\mathbf{x}_p$ is one solution, the complete solution set is
> $$
> \{\mathbf{x}_p+\mathbf{x}_n:\ \mathbf{x}_n\in\operatorname{null}(\mathbf{A})\}.
> $$
> 3. The solution is unique iff $\operatorname{null}(\mathbf{A})=\{\mathbf{0}\}$, i.e. iff the columns of $\mathbf{A}$ are independent. More generally, $T$ is one-to-one iff $\ker T=\{\mathbf{0}\}$.

*Proof of 2.* If $\mathbf{A}\mathbf{x}=\mathbf{b}$ then $\mathbf{A}(\mathbf{x}-\mathbf{x}_p)=\mathbf{b}-\mathbf{b}=\mathbf{0}$, so $\mathbf{x}-\mathbf{x}_p\in\operatorname{null}(\mathbf{A})$; conversely $\mathbf{A}(\mathbf{x}_p+\mathbf{x}_n)=\mathbf{b}+\mathbf{0}$.

This is the same "particular plus homogeneous" structure as the solution of a linear ODE (Worked Example 9.2): the null space of the operator $L=D^2+3D+2$ is the space of natural responses, and any forced response is one particular solution plus an element of that null space.

### Computing bases

All four subspaces are read off from the reduced row echelon form $\mathbf{R}=\operatorname{rref}(\mathbf{A})$ ([Ch. 8](08-matrices-and-linear-systems.md)).

| Subspace | Basis | Dimension |
|---|---|---|
| $\operatorname{col}(\mathbf{A})$ | the **pivot columns of the original** $\mathbf{A}$ | $r$ |
| $\operatorname{row}(\mathbf{A})$ | the nonzero rows of $\mathbf{R}$ | $r$ |
| $\operatorname{null}(\mathbf{A})$ | special solutions: one per free variable | $n-r$ |
| $\operatorname{null}(\mathbf{A}^\mathsf{T})$ | null space of $\mathbf{A}^\mathsf{T}$ (or last $m-r$ rows of $\mathbf{E}$ where $\mathbf{E}\mathbf{A}=\mathbf{R}$) | $m-r$ |

Why the pivot columns of $\mathbf{A}$ (not of $\mathbf{R}$)? Row operations change the column space, but they do **not** change linear relations among the columns, because $\mathbf{A}\mathbf{x}=\mathbf{0}\iff\mathbf{R}\mathbf{x}=\mathbf{0}$. In $\mathbf{R}$ the pivot columns are obviously independent and every other column is a combination of them; the same relations hold among the columns of $\mathbf{A}$.

Row operations replace rows by combinations of rows, so the row space *is* unchanged, and the nonzero rows of $\mathbf{R}$ are independent (each has a leading 1 where the others have 0).

### Rank–nullity

> [!theorem] Rank–nullity theorem
> For $\mathbf{A}\in\mathbb{F}^{m\times n}$ (or a linear map $T:V\to W$ with $\dim V=n$),
> $$
> \operatorname{rank}(\mathbf{A})+\dim\operatorname{null}(\mathbf{A})=n\quad(\text{number of columns}).
> $$
> Moreover $\dim\operatorname{row}(\mathbf{A})=\dim\operatorname{col}(\mathbf{A})=r$ ("row rank equals column rank"), and $\dim\operatorname{null}(\mathbf{A}^\mathsf{T})=m-r$.

*Proof (matrix form).* Every column of $\mathbf{R}$ is either a pivot column ($r$ of them) or a free column. The null space has one special solution per free column, and these are independent (each has a 1 in its own free position and 0 in the other free positions), so the nullity is $n-r$. The column space has the $r$ pivot columns as basis. The row space has the $r$ nonzero rows as basis. Applying the result to $\mathbf{A}^\mathsf{T}$ ($n$ rows, $m$ columns, rank $r$) gives $\dim\operatorname{null}(\mathbf{A}^\mathsf{T})=m-r$.

*Proof (abstract form).* Choose a basis $\mathbf{k}_1,\dots,\mathbf{k}_q$ of $\ker T$ and extend it to a basis $\mathbf{k}_1,\dots,\mathbf{k}_q,\mathbf{v}_1,\dots,\mathbf{v}_s$ of $V$, so $q+s=n$. Then $T(\mathbf{v}_1),\dots,T(\mathbf{v}_s)$ span the range (the $T(\mathbf{k}_i)$ are zero) and are independent (a relation $\sum c_iT(\mathbf{v}_i)=\mathbf{0}$ puts $\sum c_i\mathbf{v}_i$ in the kernel, which forces all $c_i=0$ by independence of the full basis). So $\dim\operatorname{range}T=s=n-q$.

In words: **each input dimension either survives to the output or is lost in the null space.** The map cannot create dimensions.

> [!example] Worked Example 9.18 — All four subspaces of a $3\times4$ matrix
> $$
> \mathbf{A}=\begin{bmatrix}1&2&1&0\\2&4&0&2\\3&6&1&2\end{bmatrix}.
> $$
> **Row-reduce.** $R_2-2R_1$, $R_3-3R_1$ give rows $(0,0,-2,2)$ and $(0,0,-2,2)$; then $R_3-R_2$ kills row 3; scale row 2 by $-\tfrac12$ and clear above:
> $$
> \mathbf{R}=\begin{bmatrix}1&2&0&1\\0&0&1&-1\\0&0&0&0\end{bmatrix}.
> $$
> Pivots in columns 1 and 3, so $r=2$; free variables $x_2,x_4$.
>
> **Column space** (pivot columns of $\mathbf{A}$): $\{(1,2,3)^\mathsf{T},(1,0,1)^\mathsf{T}\}$, dimension 2 — a plane in $\mathbb{R}^3$.
>
> **Row space:** $\{(1,2,0,1),(0,0,1,-1)\}$, dimension 2.
>
> **Null space.** From $\mathbf{R}$: $x_1=-2x_2-x_4$, $x_3=x_4$. Special solutions:
> $$
> x_2=1,x_4=0:\ \mathbf{n}_1=(-2,1,0,0)^\mathsf{T};\qquad x_2=0,x_4=1:\ \mathbf{n}_2=(-1,0,1,1)^\mathsf{T}.
> $$
> Dimension $2=n-r=4-2$ ✓. Check: $\mathbf{A}\mathbf{n}_2=(-1+1,\ -2+2,\ -3+1+2)=\mathbf{0}$ ✓.
>
> **Left null space.** Solve $\mathbf{A}^\mathsf{T}\mathbf{y}=\mathbf{0}$; one finds $\mathbf{y}=(-1,-1,1)^\mathsf{T}$ (row 3 $=$ row 1 $+$ row 2). Dimension $1=m-r=3-2$ ✓.
>
> **Solvability.** $\mathbf{A}\mathbf{x}=\mathbf{b}$ is solvable iff $\mathbf{b}\in\operatorname{col}(\mathbf{A})$, which here is equivalent to $\mathbf{y}^\mathsf{T}\mathbf{b}=0$, i.e. $b_3=b_1+b_2$. For $\mathbf{b}=(2,2,4)^\mathsf{T}$: $4=2+2$ ✓. Reducing $[\mathbf{A}\,|\,\mathbf{b}]$ gives $x_1+2x_2+x_4=1$, $x_3-x_4=1$; with the free variables zero, $\mathbf{x}_p=(1,0,1,0)^\mathsf{T}$ (check: $\mathbf{A}\mathbf{x}_p=(2,2,4)^\mathsf{T}$ ✓). Complete solution:
> $$
> \mathbf{x}=\begin{bmatrix}1\\0\\1\\0\end{bmatrix}+s\begin{bmatrix}-2\\1\\0\\0\end{bmatrix}+t\begin{bmatrix}-1\\0\\1\\1\end{bmatrix},\qquad s,t\in\mathbb{R}.
> $$
> For $\mathbf{b}=(1,1,1)^\mathsf{T}$: $1\neq1+1$, no solution — $\mathbf{b}$ lies off the column-space plane. Section 9.11 finds the best approximate solution.

The left null space acts as a set of **consistency checks** on the right-hand side: $\mathbf{A}\mathbf{x}=\mathbf{b}$ is solvable iff $\mathbf{y}^\mathsf{T}\mathbf{b}=0$ for every $\mathbf{y}\in\operatorname{null}(\mathbf{A}^\mathsf{T})$. (Proof: $\mathbf{y}^\mathsf{T}\mathbf{A}\mathbf{x}=(\mathbf{A}^\mathsf{T}\mathbf{y})^\mathsf{T}\mathbf{x}=0$ for every $\mathbf{x}$, so every $\mathbf{b}$ in the column space is orthogonal to the left null space; Section 9.7 shows the converse.)

> [!example] Worked Example 9.19 — The incidence matrix of a circuit graph
> A circuit has 4 nodes and 5 branches: $1\!:\,1\to2$, $2\!:\,2\to3$, $3\!:\,3\to1$, $4\!:\,3\to4$, $5\!:\,4\to1$. The **node–branch incidence matrix** $\mathbf{A}$ ($4\times5$) has $a_{kb}=+1$ if branch $b$ leaves node $k$, $-1$ if it enters, 0 otherwise:
> $$
> \mathbf{A}=\begin{bmatrix}1&0&-1&0&-1\\-1&1&0&0&0\\0&-1&1&1&0\\0&0&0&-1&1\end{bmatrix}.
> $$
> **KCL** says $\mathbf{A}\mathbf{i}=\mathbf{0}$ for the branch-current vector $\mathbf{i}$: allowed currents form $\operatorname{null}(\mathbf{A})$.
>
> **Rank.** Every column has one $+1$ and one $-1$, so the four rows sum to zero: $\mathbf{y}=(1,1,1,1)^\mathsf{T}$ is in the left null space and $r\le3$. Deleting the reference node (row 4) leaves three independent rows, so $r=3=n-1$.
>
> **Null space.** Dimension $b-r=5-3=2$: two independent loop currents. A basis:
> $$
> \boldsymbol{\ell}_1=(1,1,1,0,0)^\mathsf{T}\ (\text{loop }1\to2\to3\to1),\qquad \boldsymbol{\ell}_2=(1,1,0,1,1)^\mathsf{T}\ (\text{loop }1\to2\to3\to4\to1).
> $$
> Check $\mathbf{A}\boldsymbol{\ell}_2$: row 1: $1-0-1=0$; row 2: $-1+1=0$; row 3: $-1+1=0$; row 4: $-1+1=0$ ✓.
>
> **Interpretation.** Every KCL-compatible current distribution is $\mathbf{i}=I_1\boldsymbol{\ell}_1+I_2\boldsymbol{\ell}_2$ — mesh analysis needs exactly $b-n+1=2$ unknowns. Dually, branch voltages derived from node potentials, $\mathbf{v}=\mathbf{A}^\mathsf{T}\mathbf{e}$, fill $\operatorname{row}(\mathbf{A})$ of dimension $n-1=3$, and the node potentials are determined only up to the left-null-space vector $(1,1,1,1)$ — a common offset, removed by grounding one node.

> [!example] Worked Example 9.20 — Kernel and range of differentiation
> For $D:\mathcal{P}_3\to\mathcal{P}_2$ (Worked Example 9.16): $\ker D=\{\text{constants}\}=\operatorname{span}\{1\}$, nullity 1. Every $q\in\mathcal{P}_2$ is the derivative of some cubic (integrate it), so $\operatorname{range}D=\mathcal{P}_2$, rank 3. Rank–nullity: $3+1=4=\dim\mathcal{P}_3$ ✓. The matrix $[D]$ of Worked Example 9.16 shows the same thing: three pivots, one free (first) column. The "constant of integration" of calculus is precisely the freedom given by the null space.

> [!theorem] The invertible matrix theorem (vector-space form)
> For a square $\mathbf{A}\in\mathbb{F}^{n\times n}$ the following are equivalent:
> $\mathbf{A}$ is invertible; $\operatorname{rank}\mathbf{A}=n$; $\operatorname{null}(\mathbf{A})=\{\mathbf{0}\}$; $\operatorname{col}(\mathbf{A})=\mathbb{F}^n$; the columns form a basis of $\mathbb{F}^n$; the rows form a basis; $\det\mathbf{A}\neq0$; $\mathbf{A}\mathbf{x}=\mathbf{b}$ has exactly one solution for every $\mathbf{b}$; $T(\mathbf{x})=\mathbf{A}\mathbf{x}$ is one-to-one; $T$ is onto.

For a square matrix, one-to-one and onto are *equivalent* — a direct consequence of rank–nullity ($r=n\iff n-r=0$). For rectangular matrices they are not: a tall matrix ($m>n$) can be one-to-one but never onto; a wide matrix ($m<n$) can be onto but never one-to-one.

> [!warning] Common Pitfall — reading the column space from the wrong matrix
> The pivot columns of $\mathbf{R}=\operatorname{rref}(\mathbf{A})$ are standard unit vectors such as $(1,0,0)^\mathsf{T}$ and $(0,1,0)^\mathsf{T}$; they span a *different* plane from $\operatorname{col}(\mathbf{A})$ in general. Use the pivot **positions** from $\mathbf{R}$ but take the **columns** from $\mathbf{A}$. (For the row space it is the other way round: the rows of $\mathbf{R}$ are fine.)

> [!tip] Engineering Note — what a measurement cannot see
> If a sensor array produces readings $\mathbf{y}=\mathbf{H}\mathbf{x}$ of an unknown source vector $\mathbf{x}$, any component of $\mathbf{x}$ in $\operatorname{null}(\mathbf{H})$ is invisible: adding it changes nothing in $\mathbf{y}$. In tomography, phased arrays and state estimation (observability in control), the first design question is therefore "what is the null space of the measurement matrix?" A wide $\mathbf{H}$ (fewer measurements than unknowns) always has a nontrivial null space, so extra assumptions — minimum norm, sparsity, smoothness — are needed to choose among the infinitely many consistent $\mathbf{x}$.

## 9.7 Orthogonality

### Motivation

So far nothing has been said about lengths or angles; a vector space is pure algebra. Yet the most useful bases in engineering are **orthogonal** ones: the $\hat{\mathbf{x}},\hat{\mathbf{y}},\hat{\mathbf{z}}$ axes of field theory, the in-phase and quadrature carriers $\cos\omega_ct$ and $\sin\omega_ct$ of a modulator, the harmonics of a Fourier series, the subcarriers of OFDM. With an orthogonal basis, coordinates are found by a single dot product each instead of by solving a linear system, energy splits into a sum of independent parts (Pythagoras), and signals sent along orthogonal directions do not interfere. This section develops orthogonality in $\mathbb{R}^n$ with the ordinary dot product; Section 9.8 generalises it to complex vectors and to functions.

### Dot product, length and angle

> [!definition] Dot product, norm, distance
> For $\mathbf{x},\mathbf{y}\in\mathbb{R}^n$,
> $$
> \mathbf{x}\cdot\mathbf{y}=\mathbf{x}^\mathsf{T}\mathbf{y}=\sum_{k=1}^nx_ky_k,\qquad \lVert\mathbf{x}\rVert=\sqrt{\mathbf{x}\cdot\mathbf{x}},\qquad d(\mathbf{x},\mathbf{y})=\lVert\mathbf{x}-\mathbf{y}\rVert .
> $$
> A **unit vector** has $\lVert\mathbf{x}\rVert=1$; **normalising** $\mathbf{x}\neq\mathbf{0}$ means forming $\hat{\mathbf{x}}=\mathbf{x}/\lVert\mathbf{x}\rVert$.

> [!theorem] Cauchy–Schwarz inequality
> $$
> \lvert\mathbf{x}\cdot\mathbf{y}\rvert\le\lVert\mathbf{x}\rVert\,\lVert\mathbf{y}\rVert,
> $$
> with equality iff $\mathbf{x}$ and $\mathbf{y}$ are parallel (one is a scalar multiple of the other).

*Proof.* If $\mathbf{y}=\mathbf{0}$ both sides vanish. Otherwise, for every real $s$,
$$
0\le\lVert\mathbf{x}-s\mathbf{y}\rVert^2=\lVert\mathbf{x}\rVert^2-2s\,\mathbf{x}\cdot\mathbf{y}+s^2\lVert\mathbf{y}\rVert^2 .
$$
Choose the minimising value $s=\mathbf{x}\cdot\mathbf{y}/\lVert\mathbf{y}\rVert^2$: $0\le\lVert\mathbf{x}\rVert^2-(\mathbf{x}\cdot\mathbf{y})^2/\lVert\mathbf{y}\rVert^2$, which rearranges to the inequality. Equality means $\mathbf{x}-s\mathbf{y}=\mathbf{0}$.

The proof is worth remembering, because the same "minimise $\lVert\mathbf{x}-s\mathbf{y}\rVert^2$" step is a projection (Section 9.9) and it works verbatim in any inner-product space.

Cauchy–Schwarz guarantees $-1\le\dfrac{\mathbf{x}\cdot\mathbf{y}}{\lVert\mathbf{x}\rVert\lVert\mathbf{y}\rVert}\le1$, so the **angle** between nonzero vectors can be *defined* in any dimension by

$$
\cos\theta=\frac{\mathbf{x}\cdot\mathbf{y}}{\lVert\mathbf{x}\rVert\,\lVert\mathbf{y}\rVert},\qquad 0\le\theta\le\pi .
$$

It also gives the **triangle inequality** $\lVert\mathbf{x}+\mathbf{y}\rVert\le\lVert\mathbf{x}\rVert+\lVert\mathbf{y}\rVert$: expand $\lVert\mathbf{x}+\mathbf{y}\rVert^2=\lVert\mathbf{x}\rVert^2+2\mathbf{x}\cdot\mathbf{y}+\lVert\mathbf{y}\rVert^2\le(\lVert\mathbf{x}\rVert+\lVert\mathbf{y}\rVert)^2$.

### Orthogonal vectors and sets

> [!definition] Orthogonality
> $\mathbf{x}$ and $\mathbf{y}$ are **orthogonal**, $\mathbf{x}\perp\mathbf{y}$, if $\mathbf{x}\cdot\mathbf{y}=0$. A set $\{\mathbf{u}_1,\dots,\mathbf{u}_k\}$ is **orthogonal** if $\mathbf{u}_i\cdot\mathbf{u}_j=0$ for all $i\neq j$, and **orthonormal** if in addition every $\lVert\mathbf{u}_i\rVert=1$, i.e. $\mathbf{u}_i\cdot\mathbf{u}_j=\delta_{ij}$.

> [!theorem] Pythagoras
> $\mathbf{x}\perp\mathbf{y}\iff\lVert\mathbf{x}+\mathbf{y}\rVert^2=\lVert\mathbf{x}\rVert^2+\lVert\mathbf{y}\rVert^2$. For an orthogonal set, $\big\lVert\sum_i\mathbf{u}_i\big\rVert^2=\sum_i\lVert\mathbf{u}_i\rVert^2$.

> [!theorem] Orthogonal sets are independent
> An orthogonal set of nonzero vectors is linearly independent.

*Proof.* Take the dot product of $\sum_ic_i\mathbf{u}_i=\mathbf{0}$ with $\mathbf{u}_j$: all cross terms vanish, leaving $c_j\lVert\mathbf{u}_j\rVert^2=0$, so $c_j=0$.

Hence $n$ mutually orthogonal nonzero vectors in $\mathbb{R}^n$ are automatically a basis. Their real advantage is how easy coordinates become:

> [!theorem] Coordinates in an orthogonal basis
> If $\{\mathbf{u}_1,\dots,\mathbf{u}_n\}$ is an orthogonal basis then
> $$
> \mathbf{x}=\sum_{j=1}^n\frac{\mathbf{x}\cdot\mathbf{u}_j}{\mathbf{u}_j\cdot\mathbf{u}_j}\,\mathbf{u}_j ,
> $$
> and if it is orthonormal, $\mathbf{x}=\sum_j(\mathbf{x}\cdot\mathbf{u}_j)\,\mathbf{u}_j$ and $\lVert\mathbf{x}\rVert^2=\sum_j(\mathbf{x}\cdot\mathbf{u}_j)^2$ (Parseval).

*Proof.* Write $\mathbf{x}=\sum c_i\mathbf{u}_i$ and dot with $\mathbf{u}_j$: $\mathbf{x}\cdot\mathbf{u}_j=c_j\,\mathbf{u}_j\cdot\mathbf{u}_j$.

Compare: in a general basis we must solve an $n\times n$ system; in an orthogonal basis each coordinate is computed independently. This is *exactly* the structure of the Fourier coefficient formulas of [Ch. 15](../part-5-transforms/15-fourier-series.md), where the dot product becomes an integral.

> [!example] Worked Example 9.21 — Lengths and angles
> Find the angle between $\mathbf{x}=(1,2,2)^\mathsf{T}$ and $\mathbf{y}=(2,0,1)^\mathsf{T}$ and verify Cauchy–Schwarz.
>
> $\mathbf{x}\cdot\mathbf{y}=2+0+2=4$, $\lVert\mathbf{x}\rVert=\sqrt{1+4+4}=3$, $\lVert\mathbf{y}\rVert=\sqrt{5}=2.236$.
> $$
> \cos\theta=\frac{4}{3\sqrt5}=0.5963,\qquad\theta=53.4^\circ .
> $$
> Cauchy–Schwarz: $4\le3\times2.236=6.708$ ✓. Triangle inequality: $\lVert\mathbf{x}+\mathbf{y}\rVert=\lVert(3,2,3)\rVert=\sqrt{22}=4.690\le3+2.236$ ✓.

> [!example] Worked Example 9.22 — Expanding in an orthonormal basis
> Verify that
> $$
> \mathbf{u}_1=\tfrac{1}{\sqrt3}(1,1,1),\qquad\mathbf{u}_2=\tfrac{1}{\sqrt2}(1,-1,0),\qquad\mathbf{u}_3=\tfrac{1}{\sqrt6}(1,1,-2)
> $$
> is an orthonormal basis of $\mathbb{R}^3$, and expand $\mathbf{x}=(3,1,2)$.
>
> **Orthogonality.** $(1,1,1)\cdot(1,-1,0)=0$, $(1,1,1)\cdot(1,1,-2)=0$, $(1,-1,0)\cdot(1,1,-2)=0$ ✓. Each vector has been divided by its length ($\sqrt3$, $\sqrt2$, $\sqrt6$), so all are unit vectors.
>
> **Coordinates** by dot products:
> $$
> \mathbf{x}\cdot\mathbf{u}_1=\tfrac{6}{\sqrt3}=2\sqrt3,\qquad \mathbf{x}\cdot\mathbf{u}_2=\tfrac{2}{\sqrt2}=\sqrt2,\qquad \mathbf{x}\cdot\mathbf{u}_3=\tfrac{3+1-4}{\sqrt6}=0 .
> $$
> So $\mathbf{x}=2\sqrt3\,\mathbf{u}_1+\sqrt2\,\mathbf{u}_2=2(1,1,1)+(1,-1,0)=(3,1,2)$ ✓.
>
> **Parseval.** $(2\sqrt3)^2+(\sqrt2)^2+0^2=12+2=14=\lVert\mathbf{x}\rVert^2=9+1+4$ ✓.
>
> **Interpretation.** $\mathbf{u}_1$ is the "average" (DC) direction: $\mathbf{x}\cdot\mathbf{u}_1\,\mathbf{u}_1=(2,2,2)$ is the mean value of the entries repeated. The rest, $(1,-1,0)$, is the zero-mean (AC) part. Energy splits as $12$ (DC) $+\ 2$ (AC).

### Orthogonal matrices

> [!definition] Orthogonal matrix
> A square real matrix $\mathbf{Q}$ is **orthogonal** if its columns are orthonormal, equivalently
> $$
> \mathbf{Q}^\mathsf{T}\mathbf{Q}=\mathbf{I},\qquad\text{i.e.}\qquad\mathbf{Q}^{-1}=\mathbf{Q}^\mathsf{T}.
> $$
> (The complex analogue, $\mathbf{U}^\mathsf{H}\mathbf{U}=\mathbf{I}$, is a **unitary** matrix.)

Orthogonal matrices preserve dot products, lengths and angles:

$$
(\mathbf{Q}\mathbf{x})\cdot(\mathbf{Q}\mathbf{y})=\mathbf{x}^\mathsf{T}\mathbf{Q}^\mathsf{T}\mathbf{Q}\mathbf{y}=\mathbf{x}\cdot\mathbf{y},\qquad\lVert\mathbf{Q}\mathbf{x}\rVert=\lVert\mathbf{x}\rVert .
$$

They are the rotations and reflections, and $\det\mathbf{Q}=\pm1$. The rotation $\mathbf{R}_\theta$ of Worked Example 9.15 is orthogonal; so is the matrix $[\mathbf{u}_1\ \mathbf{u}_2\ \mathbf{u}_3]$ of Worked Example 9.22. Numerically they are ideal: $\kappa(\mathbf{Q})=1$, so solving $\mathbf{Q}\mathbf{x}=\mathbf{b}$ by $\mathbf{x}=\mathbf{Q}^\mathsf{T}\mathbf{b}$ never amplifies errors. This is why orthogonal transforms (the DFT scaled by $1/\sqrt N$, the discrete cosine transform, Householder reflections in QR) dominate numerical signal processing.

### Orthogonal complements

> [!definition] Orthogonal complement
> For a subspace $W$ of $\mathbb{R}^n$, its **orthogonal complement** is
> $$
> W^\perp=\{\mathbf{x}\in\mathbb{R}^n:\ \mathbf{x}\cdot\mathbf{w}=0\ \text{for all}\ \mathbf{w}\in W\}.
> $$

$W^\perp$ is a subspace, and it suffices to check orthogonality against a *spanning set* of $W$. Key properties:

$$
W\cap W^\perp=\{\mathbf{0}\},\qquad \dim W+\dim W^\perp=n,\qquad (W^\perp)^\perp=W,\qquad \mathbb{R}^n=W\oplus W^\perp .
$$

> [!theorem] Fundamental theorem of linear algebra (orthogonality part)
> For $\mathbf{A}\in\mathbb{R}^{m\times n}$:
> $$
> \operatorname{null}(\mathbf{A})=\operatorname{row}(\mathbf{A})^\perp\ \text{in}\ \mathbb{R}^n,\qquad \operatorname{null}(\mathbf{A}^\mathsf{T})=\operatorname{col}(\mathbf{A})^\perp\ \text{in}\ \mathbb{R}^m .
> $$

*Proof.* $\mathbf{A}\mathbf{x}=\mathbf{0}$ says that the dot product of $\mathbf{x}$ with every row of $\mathbf{A}$ is zero, i.e. $\mathbf{x}$ is orthogonal to the row space. Apply this to $\mathbf{A}^\mathsf{T}$ for the second statement.

Together with rank–nullity this gives the complete picture: $\mathbb{R}^n$ splits orthogonally into the row space ($\dim r$) and the null space ($\dim n-r$); $\mathbb{R}^m$ splits into the column space ($\dim r$) and the left null space ($\dim m-r$). $\mathbf{A}$ maps the row space one-to-one onto the column space and sends the null space to zero. The second statement also completes the solvability test of Section 9.6: $\mathbf{b}\in\operatorname{col}(\mathbf{A})$ **iff** $\mathbf{b}\perp\operatorname{null}(\mathbf{A}^\mathsf{T})$.

> [!example] Worked Example 9.23 — Complements and the four subspaces
> (a) Find $W^\perp$ for $W=\operatorname{span}\{(1,1,1)\}$.
>
> $\mathbf{x}\perp(1,1,1)\iff x+y+z=0$. So $W^\perp$ is the plane $x+y+z=0$, with basis $\{(-1,1,0),(-1,0,1)\}$; $\dim W+\dim W^\perp=1+2=3$ ✓. In the language of Worked Example 9.22: the orthogonal complement of the DC direction is the subspace of zero-mean vectors.
>
> (b) For $\mathbf{A}$ of Worked Example 9.18, check $\operatorname{null}(\mathbf{A})\perp\operatorname{row}(\mathbf{A})$:
> $$
> (-2,1,0,0)\cdot(1,2,0,1)=0,\ \ (-2,1,0,0)\cdot(0,0,1,-1)=0,\ \ (-1,0,1,1)\cdot(1,2,0,1)=0,\ \ (-1,0,1,1)\cdot(0,0,1,-1)=0\ ✓.
> $$
> And $\operatorname{null}(\mathbf{A}^\mathsf{T})=\operatorname{span}\{(-1,-1,1)\}\perp\operatorname{col}(\mathbf{A})$: $(-1,-1,1)\cdot(1,2,3)=0$ and $(-1,-1,1)\cdot(1,0,1)=0$ ✓. Dimensions $2+2=4=n$ and $2+1=3=m$ ✓.

> [!warning] Common Pitfall — orthogonal is not the same as independent
> Orthogonal nonzero vectors are independent, but independent vectors need not be orthogonal: $(1,0)$ and $(1,1)$ are independent and meet at $45^\circ$. Also, "orthogonal matrix" means ortho*normal* columns; a matrix with merely orthogonal columns satisfies $\mathbf{Q}^\mathsf{T}\mathbf{Q}=\text{diagonal}$, not $\mathbf{I}$, and its inverse is *not* its transpose.

> [!tip] Engineering Note — I/Q modulation and orthogonal channels
> Over one symbol period $T$ containing an integer number of carrier cycles, $\cos\omega_ct$ and $\sin\omega_ct$ are orthogonal (Section 9.8). A QAM transmitter sends $s(t)=a\cos\omega_ct-b\sin\omega_ct$; the receiver recovers $a$ and $b$ separately by taking inner products with each carrier — exactly the orthonormal-coordinate formula of this section. OFDM extends this to hundreds of mutually orthogonal subcarriers, and CDMA to orthogonal spreading codes (rows of a Hadamard matrix, which satisfies $\mathbf{H}^\mathsf{T}\mathbf{H}=N\mathbf{I}$).

## 9.8 Inner products

### Motivation

The dot product gave $\mathbb{R}^n$ lengths, angles and orthogonality. To do the same for complex phasor vectors, for vectors whose components have different importance (weights), and above all for **signals**, we abstract its essential properties. The result, an **inner product**, turns the function space of finite-energy signals into a geometric space in which a signal has a length (its RMS value or the square root of its energy), two signals have an angle (their correlation), and the harmonics $\cos k\omega_0t$, $\sin k\omega_0t$ are mutually perpendicular axes. This is the bridge from linear algebra to Fourier analysis.

### Definition

> [!definition] Inner product
> An **inner product** on a vector space $V$ over $\mathbb{F}$ assigns to each pair $\mathbf{u},\mathbf{v}$ a scalar $\langle\mathbf{u},\mathbf{v}\rangle\in\mathbb{F}$ such that for all $\mathbf{u},\mathbf{v},\mathbf{w}$ and scalars $a,b$:
> 1. **Linearity in the first argument:** $\langle a\mathbf{u}+b\mathbf{w},\mathbf{v}\rangle=a\langle\mathbf{u},\mathbf{v}\rangle+b\langle\mathbf{w},\mathbf{v}\rangle$;
> 2. **Conjugate symmetry:** $\langle\mathbf{u},\mathbf{v}\rangle=\langle\mathbf{v},\mathbf{u}\rangle^*$ (plain symmetry when $\mathbb{F}=\mathbb{R}$);
> 3. **Positivity:** $\langle\mathbf{u},\mathbf{u}\rangle>0$ for $\mathbf{u}\neq\mathbf{0}$ (and $\langle\mathbf{0},\mathbf{0}\rangle=0$).
>
> A vector space with an inner product is an **inner-product space**. The induced **norm** is $\lVert\mathbf{u}\rVert=\sqrt{\langle\mathbf{u},\mathbf{u}\rangle}$, and $\mathbf{u}\perp\mathbf{v}$ means $\langle\mathbf{u},\mathbf{v}\rangle=0$.

Properties 1 and 2 together give **conjugate linearity in the second argument**: $\langle\mathbf{u},a\mathbf{v}\rangle=a^*\langle\mathbf{u},\mathbf{v}\rangle$. This book follows the engineering and signal-processing convention (linear in the first slot, conjugated in the second); physics texts often use the opposite convention. Both give the same norms and the same notion of orthogonality.

> [!theorem] Results that hold in every inner-product space
> - **Cauchy–Schwarz:** $\lvert\langle\mathbf{u},\mathbf{v}\rangle\rvert\le\lVert\mathbf{u}\rVert\,\lVert\mathbf{v}\rVert$, equality iff $\mathbf{u},\mathbf{v}$ are parallel.
> - **Triangle inequality:** $\lVert\mathbf{u}+\mathbf{v}\rVert\le\lVert\mathbf{u}\rVert+\lVert\mathbf{v}\rVert$.
> - **Pythagoras:** if $\mathbf{u}\perp\mathbf{v}$ then $\lVert\mathbf{u}+\mathbf{v}\rVert^2=\lVert\mathbf{u}\rVert^2+\lVert\mathbf{v}\rVert^2$.
> - **Parallelogram law:** $\lVert\mathbf{u}+\mathbf{v}\rVert^2+\lVert\mathbf{u}-\mathbf{v}\rVert^2=2\lVert\mathbf{u}\rVert^2+2\lVert\mathbf{v}\rVert^2$.
> - Orthogonal nonzero vectors are independent, and coordinates in an orthogonal basis are $c_j=\langle\mathbf{x},\mathbf{u}_j\rangle/\langle\mathbf{u}_j,\mathbf{u}_j\rangle$.

The proofs of Section 9.7 go through unchanged, with $\mathbf{x}\cdot\mathbf{y}$ replaced by $\langle\mathbf{x},\mathbf{y}\rangle$. (In the complex case of Cauchy–Schwarz, minimise $\lVert\mathbf{u}-s\mathbf{v}\rVert^2$ over complex $s$; the minimiser is $s=\langle\mathbf{u},\mathbf{v}\rangle/\lVert\mathbf{v}\rVert^2$.)

### Inner products on $\mathbb{R}^n$ and $\mathbb{C}^n$

**Standard complex inner product.**

$$
\langle\mathbf{u},\mathbf{v}\rangle=\sum_{k=1}^nu_k\,v_k^*=\mathbf{v}^\mathsf{H}\mathbf{u},\qquad \lVert\mathbf{u}\rVert^2=\sum_k\lvert u_k\rvert^2 .
$$

The conjugate is essential. Without it, $(1,j)$ would have "squared length" $1^2+j^2=0$, violating positivity.

**Weighted inner product.** For a symmetric positive-definite matrix $\mathbf{W}$ ([Ch. 10](10-eigenvalues-and-eigenvectors.md)), $\langle\mathbf{u},\mathbf{v}\rangle_\mathbf{W}=\mathbf{v}^\mathsf{T}\mathbf{W}\mathbf{u}$. The diagonal case $\mathbf{W}=\operatorname{diag}(w_1,\dots,w_n)$ with $w_k>0$ weights each component by its reliability — the basis of weighted least squares (Section 9.11).

> [!example] Worked Example 9.24 — A complex inner product
> Let $\mathbf{u}=(1+j,\ 2)^\mathsf{T}$ and $\mathbf{v}=(j,\ 1-j)^\mathsf{T}$ in $\mathbb{C}^2$. Compute $\langle\mathbf{u},\mathbf{v}\rangle$, both norms, and check Cauchy–Schwarz.
>
> $$
> \langle\mathbf{u},\mathbf{v}\rangle=(1+j)(j)^*+2(1-j)^*=(1+j)(-j)+2(1+j)=(-j+1)+(2+2j)=3+j .
> $$
> $$
> \lVert\mathbf{u}\rVert^2=\lvert1+j\rvert^2+\lvert2\rvert^2=2+4=6,\qquad\lVert\mathbf{v}\rVert^2=\lvert j\rvert^2+\lvert1-j\rvert^2=1+2=3 .
> $$
> Cauchy–Schwarz: $\lvert3+j\rvert=\sqrt{10}=3.162\le\sqrt6\sqrt3=\sqrt{18}=4.243$ ✓. Conjugate symmetry: $\langle\mathbf{v},\mathbf{u}\rangle=j(1-j)+(1-j)\cdot2=(j+1)+(2-2j)=3-j=\langle\mathbf{u},\mathbf{v}\rangle^*$ ✓.

### Inner products of functions and signals

For real continuous functions on $[a,b]$, the natural inner product replaces the sum over components by an integral over $t$:

$$
\langle f,g\rangle=\int_a^bf(t)\,g(t)\,dt,\qquad\lVert f\rVert^2=\int_a^bf(t)^2\,dt .
$$

For complex-valued signals, $\langle f,g\rangle=\int_a^bf(t)\,g^*(t)\,dt$. Positivity holds because a continuous function with $\int f^2=0$ is identically zero. (On $L^2$, functions that differ only on a set of zero length — say at isolated points — are regarded as the same vector, so that positivity still holds.)

The analogy with $\mathbb{R}^n$ is exact if you think of a function as a vector with a continuum of components indexed by $t$: sampling $f$ and $g$ at spacing $\Delta t$ gives vectors whose dot product times $\Delta t$ is a Riemann sum for $\int fg\,dt$.

Two inner products are standard in signal analysis:

| Signal class | Inner product | Squared norm | Physical meaning |
|---|---|---|---|
| energy signals on $\mathbb{R}$ | $\langle x,y\rangle=\displaystyle\int_{-\infty}^{\infty}x(t)y^*(t)\,dt$ | $E_x=\displaystyle\int\lvert x\rvert^2dt$ | energy delivered to a $1\ \Omega$ load |
| $T$-periodic (power) signals | $\langle x,y\rangle_T=\dfrac1T\displaystyle\int_{T}x(t)y^*(t)\,dt$ | $P_x=\dfrac1T\displaystyle\int_T\lvert x\rvert^2dt$ | average power in $1\ \Omega$; $\lVert x\rVert_T=$ RMS value |

So **the RMS value of a periodic waveform is its norm** in the power inner product, and the energy of a pulse is its squared norm.

### Orthogonality of harmonics

> [!theorem] Orthogonality of the trigonometric system
> Let $\omega_0=2\pi/T$. For integers $m,n\ge0$, over any interval of length $T$:
> $$
> \int_T\cos m\omega_0t\,\cos n\omega_0t\,dt=\begin{cases}0,&m\neq n\\ T/2,&m=n\ge1\\ T,&m=n=0\end{cases}\qquad
> \int_T\sin m\omega_0t\,\sin n\omega_0t\,dt=\begin{cases}0,&m\neq n\\ T/2,&m=n\ge1\end{cases}
> $$
> $$
> \int_T\cos m\omega_0t\,\sin n\omega_0t\,dt=0\quad\text{for all }m,n .
> $$
> In complex form, with the power inner product,
> $$
> \big\langle e^{jm\omega_0t},e^{jn\omega_0t}\big\rangle_T=\frac1T\int_Te^{j(m-n)\omega_0t}\,dt=\delta_{mn}\qquad(m,n\in\mathbb{Z}),
> $$
> so $\{e^{jn\omega_0t}\}_{n\in\mathbb{Z}}$ is an **orthonormal** set.

*Derivation.* The complex form is the cleanest: if $m=n$ the integrand is 1 and the integral is $T$. If $m\neq n$, put $k=m-n\neq0$:
$$
\int_{t_0}^{t_0+T}e^{jk\omega_0t}\,dt=\frac{e^{jk\omega_0t_0}\big(e^{jk\omega_0T}-1\big)}{jk\omega_0}=0,
$$
because $k\omega_0T=2\pi k$. The real formulas follow from product-to-sum identities, e.g. $\cos m\omega_0t\cos n\omega_0t=\tfrac12[\cos(m-n)\omega_0t+\cos(m+n)\omega_0t]$: each cosine integrates to zero over a period unless its frequency is zero.

```desmos-graph
left=0; right=6.2832; bottom=-1.2; top=1.2
---
y=\cos(x)\cos(2x)|blue
y=\cos(x)|grey|dashed
y=\cos(2x)|grey|dotted
y=\cos(x)^{2}|red
```
*Over one period $[0,2\pi]$, the product $\cos t\cos2t$ (blue) has equal positive and negative areas, so the inner product is $0$; the product $\cos^2t$ (red) is never negative and integrates to $\pi=T/2$.*

Consequently, if a signal is a finite trigonometric sum $x(t)=\sum_{n=-N}^Nc_ne^{jn\omega_0t}$, its coefficients are its coordinates in an orthonormal basis:

$$
c_n=\big\langle x,e^{jn\omega_0t}\big\rangle_T=\frac1T\int_Tx(t)\,e^{-jn\omega_0t}\,dt,\qquad P_x=\sum_n\lvert c_n\rvert^2 .
$$

These are the **Fourier coefficient formula** and **Parseval's relation** of [Ch. 15](../part-5-transforms/15-fourier-series.md) — here obtained purely from the orthonormal-coordinate theorem. What Chapter 15 adds is the analysis of what happens when $N\to\infty$.

> [!example] Worked Example 9.25 — Angles between polynomials
> (a) On $[-1,1]$ with $\langle f,g\rangle=\int_{-1}^1fg\,dt$: $\langle1,t\rangle=\int_{-1}^1t\,dt=0$ and $\langle t,t^2\rangle=\int_{-1}^1t^3\,dt=0$. Any even function is orthogonal to any odd function on a symmetric interval — the subspaces of even and odd functions are orthogonal complements. Norms: $\lVert1\rVert^2=2$, $\lVert t\rVert^2=2/3$, $\lVert t^2\rVert^2=2/5$. But $\langle1,t^2\rangle=2/3\neq0$: the monomials are *not* an orthogonal set (Section 9.10 fixes this).
>
> (b) On $[0,1]$ the same functions behave very differently:
> $$
> \langle t,t^2\rangle=\int_0^1t^3\,dt=\tfrac14,\quad\lVert t\rVert=\sqrt{\int_0^1t^2dt}=\tfrac{1}{\sqrt3},\quad\lVert t^2\rVert=\tfrac1{\sqrt5},
> $$
> $$
> \cos\theta=\frac{1/4}{(1/\sqrt3)(1/\sqrt5)}=\frac{\sqrt{15}}{4}=0.9682,\qquad\theta=14.5^\circ .
> $$
> On $[0,1]$, $t$ and $t^2$ are nearly parallel. This geometric fact is the root of the ill-conditioning of polynomial fitting with monomials (Section 9.11). The inner product — and hence orthogonality — depends on the interval.

> [!example] Worked Example 9.26 — Energy, correlation and distance of two pulses
> Let $x_1(t)=e^{-t}u(t)$ and $x_2(t)=e^{-2t}u(t)$ (volts, $t$ in seconds, $1\ \Omega$ reference).
>
> $$
> E_1=\lVert x_1\rVert^2=\int_0^\infty e^{-2t}dt=\tfrac12\ \text{J},\qquad E_2=\int_0^\infty e^{-4t}dt=\tfrac14\ \text{J},\qquad \langle x_1,x_2\rangle=\int_0^\infty e^{-3t}dt=\tfrac13 .
> $$
> **Correlation coefficient** (cosine of the angle):
> $$
> \rho=\frac{\langle x_1,x_2\rangle}{\lVert x_1\rVert\lVert x_2\rVert}=\frac{1/3}{\sqrt{1/8}}=\frac{2\sqrt2}{3}=0.943 .
> $$
> **Distance** (energy of the difference): $\lVert x_1-x_2\rVert^2=E_1+E_2-2\langle x_1,x_2\rangle=\tfrac12+\tfrac14-\tfrac23=\tfrac1{12}$ J. A receiver deciding between these two pulses in noise has a hard job: the pulses are nearly parallel and close together. In digital communications, error probability is governed by exactly this distance between signal vectors.

> [!example] Worked Example 9.27 — RMS value by Pythagoras
> Find the RMS value of $v(t)=3+4\cos\omega_0t-2\sin3\omega_0t$ V.
>
> The three terms are mutually orthogonal in the power inner product, so their powers add:
> $$
> P=\lVert3\rVert_T^2+\lVert4\cos\omega_0t\rVert_T^2+\lVert2\sin3\omega_0t\rVert_T^2=9+\frac{16}{2}+\frac{4}{2}=19\ \text{V}^2,
> $$
> using $\frac1T\int_T\cos^2=\frac1T\int_T\sin^2=\frac12$. Hence $V_{\text{rms}}=\sqrt{19}=4.359$ V. Without orthogonality we would have to expand $v^2$ and integrate nine cross terms; orthogonality says the six cross terms all vanish.

### Cauchy–Schwarz and the matched filter

Cauchy–Schwarz has a celebrated engineering reading. A known pulse $s(t)$ arrives in white noise and is passed through a filter whose output at the sampling instant is an inner product $\langle s,h\rangle$ with some template $h$. The signal-to-noise ratio turns out to be proportional to

$$
\frac{\lvert\langle s,h\rangle\rvert^2}{\lVert h\rVert^2}\le\frac{\lVert s\rVert^2\lVert h\rVert^2}{\lVert h\rVert^2}=\lVert s\rVert^2=E_s,
$$

with equality iff $h$ is parallel to $s$. The optimal template is the signal itself (the **matched filter**), and the best achievable SNR depends only on the pulse energy, not its shape. The derivation of the noise side belongs to [Ch. 21](../part-6-probability-statistics/21-random-processes.md); the optimisation step is pure Cauchy–Schwarz.

> [!warning] Common Pitfall — forgetting the conjugate, or the interval
> With complex vectors or signals, $\sum u_kv_k$ (no conjugate) is *not* an inner product, and "orthogonality" computed with it is meaningless. Similarly, $\cos t$ and $\cos2t$ are orthogonal over $[0,2\pi]$ but not over $[0,1]$: orthogonality of sinusoids requires an interval containing a whole number of periods of each (more precisely, of their difference frequency). Truncating a measurement to a non-integer number of periods destroys orthogonality and causes spectral leakage in the DFT.

> [!tip] Engineering Note — signal space
> Digital modulation is analysed in **signal space**: choose an orthonormal basis $\phi_1(t),\dots,\phi_N(t)$ for the set of transmitted waveforms (for QPSK, $\sqrt{2/T}\cos\omega_ct$ and $\sqrt{2/T}\sin\omega_ct$ over one symbol). Each waveform becomes a point in $\mathbb{R}^N$ — the constellation diagram — energies become squared lengths, and the receiver's correlators compute the coordinates $\langle r,\phi_k\rangle$. Every geometric statement of this chapter has a direct meaning on that diagram.

## 9.9 Projection

### Motivation

Most approximation problems in engineering have the same shape: we want to represent a complicated object (a measured data vector, a distorted waveform, a field) by the closest member of a simpler family (a straight line, a few harmonics, a low-order model). When the simpler family is a subspace and "closest" is measured by an inner-product norm — energy, RMS error, sum of squared residuals — the answer is always the **orthogonal projection**, and it is characterised by one geometric fact: *the error is perpendicular to the subspace*. Least squares (Section 9.11), truncated Fourier series ([Ch. 15](../part-5-transforms/15-fourier-series.md)) and Wiener filtering are all instances.

### Projection onto a line

Let $\mathbf{a}\neq\mathbf{0}$ span a line through the origin. The point $\mathbf{p}=\hat c\,\mathbf{a}$ closest to $\mathbf{b}$ is determined by requiring the error $\mathbf{e}=\mathbf{b}-\hat c\mathbf{a}$ to be orthogonal to $\mathbf{a}$:

$$
\langle\mathbf{b}-\hat c\,\mathbf{a},\mathbf{a}\rangle=0\quad\Longrightarrow\quad\hat c=\frac{\langle\mathbf{b},\mathbf{a}\rangle}{\langle\mathbf{a},\mathbf{a}\rangle},\qquad \mathbf{p}=\operatorname{proj}_\mathbf{a}\mathbf{b}=\frac{\langle\mathbf{b},\mathbf{a}\rangle}{\langle\mathbf{a},\mathbf{a}\rangle}\,\mathbf{a}.
$$

In $\mathbb{R}^n$ the projection is a matrix acting on $\mathbf{b}$:

$$
\mathbf{p}=\frac{\mathbf{a}\mathbf{a}^\mathsf{T}}{\mathbf{a}^\mathsf{T}\mathbf{a}}\,\mathbf{b}=\mathbf{P}\mathbf{b},\qquad\mathbf{P}=\frac{\mathbf{a}\mathbf{a}^\mathsf{T}}{\mathbf{a}^\mathsf{T}\mathbf{a}}\ \ (\text{rank }1).
$$

The scalar $\langle\mathbf{b},\mathbf{a}\rangle/\lVert\mathbf{a}\rVert$ is the **scalar component** of $\mathbf{b}$ along $\mathbf{a}$, equal to $\lVert\mathbf{b}\rVert\cos\theta$.

> [!example] Worked Example 9.28 — Projecting onto the DC direction
> Project $\mathbf{b}=(1,2,3)^\mathsf{T}$ onto $\mathbf{a}=(1,1,1)^\mathsf{T}$.
>
> $\hat c=\dfrac{\mathbf{b}\cdot\mathbf{a}}{\mathbf{a}\cdot\mathbf{a}}=\dfrac{6}{3}=2$, so $\mathbf{p}=(2,2,2)^\mathsf{T}$ and $\mathbf{e}=\mathbf{b}-\mathbf{p}=(-1,0,1)^\mathsf{T}$.
>
> **Check:** $\mathbf{e}\cdot\mathbf{a}=-1+0+1=0$ ✓. **Pythagoras:** $\lVert\mathbf{b}\rVert^2=14=\lVert\mathbf{p}\rVert^2+\lVert\mathbf{e}\rVert^2=12+2$ ✓.
>
> $\hat c=2$ is the **mean** of the entries of $\mathbf{b}$: projecting onto the all-ones vector extracts the average, and the error is the zero-mean residual. The distance from $\mathbf{b}$ to the line is $\lVert\mathbf{e}\rVert=\sqrt2$ — $\sqrt3$ times the (population) standard deviation of the entries $1,2,3$, which is $\sqrt{2/3}$.

### Projection onto a subspace

> [!theorem] Best approximation (orthogonal projection) theorem
> Let $W$ be a finite-dimensional subspace of an inner-product space $V$ and $\mathbf{b}\in V$. There is a unique $\mathbf{p}\in W$ closest to $\mathbf{b}$:
> $$
> \lVert\mathbf{b}-\mathbf{p}\rVert<\lVert\mathbf{b}-\mathbf{w}\rVert\quad\text{for every }\mathbf{w}\in W,\ \mathbf{w}\neq\mathbf{p}.
> $$
> It is characterised by the **orthogonality principle**: $\mathbf{b}-\mathbf{p}\perp W$. We write $\mathbf{p}=\operatorname{proj}_W\mathbf{b}$.

*Proof.* Suppose $\mathbf{p}\in W$ with $\mathbf{e}=\mathbf{b}-\mathbf{p}\perp W$. For any $\mathbf{w}\in W$, $\mathbf{p}-\mathbf{w}\in W$, so by Pythagoras
$$
\lVert\mathbf{b}-\mathbf{w}\rVert^2=\lVert(\mathbf{b}-\mathbf{p})+(\mathbf{p}-\mathbf{w})\rVert^2=\lVert\mathbf{e}\rVert^2+\lVert\mathbf{p}-\mathbf{w}\rVert^2\ge\lVert\mathbf{e}\rVert^2,
$$
with equality only if $\mathbf{w}=\mathbf{p}$. Existence: construct $\mathbf{p}$ with an orthogonal basis (next formula) or by the normal equations.

**With an orthogonal basis** $\{\mathbf{u}_1,\dots,\mathbf{u}_k\}$ of $W$ the projection is the sum of the projections onto each basis vector:

$$
\operatorname{proj}_W\mathbf{b}=\sum_{j=1}^k\frac{\langle\mathbf{b},\mathbf{u}_j\rangle}{\langle\mathbf{u}_j,\mathbf{u}_j\rangle}\,\mathbf{u}_j\qquad\Big(=\sum_j\langle\mathbf{b},\mathbf{u}_j\rangle\mathbf{u}_j\ \text{if orthonormal}\Big).
$$

*Check of orthogonality:* $\langle\mathbf{b}-\mathbf{p},\mathbf{u}_i\rangle=\langle\mathbf{b},\mathbf{u}_i\rangle-\frac{\langle\mathbf{b},\mathbf{u}_i\rangle}{\langle\mathbf{u}_i,\mathbf{u}_i\rangle}\langle\mathbf{u}_i,\mathbf{u}_i\rangle=0$, the cross terms having vanished by orthogonality.

> [!warning] Common Pitfall — the sum formula needs an orthogonal basis
> Adding up the projections onto each basis vector gives the projection onto their span **only** if the basis is orthogonal. For a non-orthogonal basis the individual projections "double count" the overlap. Either orthogonalise first (Gram–Schmidt, Section 9.10) or use the normal-equation formula below.

**With an arbitrary basis** — the columns of $\mathbf{A}\in\mathbb{R}^{m\times k}$, assumed independent — write $\mathbf{p}=\mathbf{A}\hat{\mathbf{x}}$. Orthogonality of the error to every column means $\mathbf{A}^\mathsf{T}(\mathbf{b}-\mathbf{A}\hat{\mathbf{x}})=\mathbf{0}$:

$$
\mathbf{A}^\mathsf{T}\mathbf{A}\,\hat{\mathbf{x}}=\mathbf{A}^\mathsf{T}\mathbf{b}\qquad(\text{normal equations}),
$$

$$
\hat{\mathbf{x}}=(\mathbf{A}^\mathsf{T}\mathbf{A})^{-1}\mathbf{A}^\mathsf{T}\mathbf{b},\qquad\mathbf{p}=\mathbf{P}\mathbf{b},\qquad\mathbf{P}=\mathbf{A}(\mathbf{A}^\mathsf{T}\mathbf{A})^{-1}\mathbf{A}^\mathsf{T}.
$$

$\mathbf{A}^\mathsf{T}\mathbf{A}$ is invertible precisely when the columns of $\mathbf{A}$ are independent: if $\mathbf{A}^\mathsf{T}\mathbf{A}\mathbf{x}=\mathbf{0}$ then $\mathbf{x}^\mathsf{T}\mathbf{A}^\mathsf{T}\mathbf{A}\mathbf{x}=\lVert\mathbf{A}\mathbf{x}\rVert^2=0$, so $\mathbf{A}\mathbf{x}=\mathbf{0}$, so $\mathbf{x}=\mathbf{0}$. For complex data replace $\mathsf{T}$ by $\mathsf{H}$.

> [!theorem] Properties of projection matrices
> The orthogonal projection matrix $\mathbf{P}$ onto a subspace $W$ satisfies
> 1. $\mathbf{P}^2=\mathbf{P}$ (**idempotent**: projecting twice changes nothing);
> 2. $\mathbf{P}^\mathsf{T}=\mathbf{P}$ (**symmetric**);
> 3. $\operatorname{col}(\mathbf{P})=W$, $\operatorname{null}(\mathbf{P})=W^\perp$, $\operatorname{rank}\mathbf{P}=\dim W$;
> 4. $\mathbf{I}-\mathbf{P}$ is the projection onto $W^\perp$, and $\mathbf{b}=\mathbf{P}\mathbf{b}+(\mathbf{I}-\mathbf{P})\mathbf{b}$ is the orthogonal splitting $\mathbb{R}^m=W\oplus W^\perp$.
>
> Conversely, any matrix with $\mathbf{P}^2=\mathbf{P}=\mathbf{P}^\mathsf{T}$ is the orthogonal projection onto its column space. If $W$ has an orthonormal basis forming the columns of $\mathbf{Q}$, then simply $\mathbf{P}=\mathbf{Q}\mathbf{Q}^\mathsf{T}$.

*Proof of 1.* $\mathbf{P}^2=\mathbf{A}(\mathbf{A}^\mathsf{T}\mathbf{A})^{-1}(\mathbf{A}^\mathsf{T}\mathbf{A})(\mathbf{A}^\mathsf{T}\mathbf{A})^{-1}\mathbf{A}^\mathsf{T}=\mathbf{P}$. Property 2 follows because $(\mathbf{A}^\mathsf{T}\mathbf{A})^{-1}$ is symmetric.

A matrix with $\mathbf{P}^2=\mathbf{P}$ but $\mathbf{P}^\mathsf{T}\neq\mathbf{P}$ is an **oblique** projection: it maps onto a subspace along a direction that is not perpendicular to it.

> [!example] Worked Example 9.29 — Projection onto a plane in $\mathbb{R}^3$
> Project $\mathbf{b}=(6,0,0)^\mathsf{T}$ onto the column space of
> $$
> \mathbf{A}=\begin{bmatrix}1&0\\1&1\\1&2\end{bmatrix}.
> $$
> **Normal equations.**
> $$
> \mathbf{A}^\mathsf{T}\mathbf{A}=\begin{bmatrix}3&3\\3&5\end{bmatrix},\qquad\mathbf{A}^\mathsf{T}\mathbf{b}=\begin{bmatrix}6\\0\end{bmatrix},\qquad \begin{bmatrix}3&3\\3&5\end{bmatrix}\begin{bmatrix}\hat x_1\\ \hat x_2\end{bmatrix}=\begin{bmatrix}6\\0\end{bmatrix}.
> $$
> Subtracting the first equation from the second: $2\hat x_2=-6$, $\hat x_2=-3$; then $3\hat x_1-9=6$, $\hat x_1=5$.
>
> **Projection and error.** $\mathbf{p}=\mathbf{A}\hat{\mathbf{x}}=5(1,1,1)-3(0,1,2)=(5,2,-1)^\mathsf{T}$, $\mathbf{e}=\mathbf{b}-\mathbf{p}=(1,-2,1)^\mathsf{T}$.
>
> **Checks.** $\mathbf{A}^\mathsf{T}\mathbf{e}=(1-2+1,\ 0-2+2)=(0,0)$ ✓ — the error is orthogonal to both columns. $\lVert\mathbf{b}\rVert^2=36=\lVert\mathbf{p}\rVert^2+\lVert\mathbf{e}\rVert^2=30+6$ ✓.
>
> **Projection matrix.** $(\mathbf{A}^\mathsf{T}\mathbf{A})^{-1}=\frac16\begin{bmatrix}5&-3\\-3&3\end{bmatrix}$ and
> $$
> \mathbf{P}=\mathbf{A}(\mathbf{A}^\mathsf{T}\mathbf{A})^{-1}\mathbf{A}^\mathsf{T}=\frac16\begin{bmatrix}5&2&-1\\2&2&2\\-1&2&5\end{bmatrix}.
> $$
> Its first column is $\mathbf{P}\mathbf{e}_1=\frac16(5,2,-1)$, consistent with $\mathbf{P}\mathbf{b}=6\,\mathbf{P}\mathbf{e}_1=(5,2,-1)$ ✓. Its trace is $12/6=2=\operatorname{rank}\mathbf{P}$ (the trace of a projection equals the dimension of the subspace). This same matrix reappears in Section 9.11 as fitting a straight line to three data points.

> [!example] Worked Example 9.30 — Projection with an orthonormal basis
> Project $\mathbf{x}=(3,1,2)^\mathsf{T}$ onto the plane $W:\ x+y+z=0$.
>
> From Worked Example 9.22, $\mathbf{u}_2=\frac1{\sqrt2}(1,-1,0)$ and $\mathbf{u}_3=\frac1{\sqrt6}(1,1,-2)$ form an orthonormal basis of $W$ (both are orthogonal to $(1,1,1)$). Then
> $$
> \operatorname{proj}_W\mathbf{x}=(\mathbf{x}\cdot\mathbf{u}_2)\mathbf{u}_2+(\mathbf{x}\cdot\mathbf{u}_3)\mathbf{u}_3=\sqrt2\cdot\tfrac1{\sqrt2}(1,-1,0)+0=(1,-1,0).
> $$
> **Alternative** via the complement: $W^\perp=\operatorname{span}\{(1,1,1)\}$, so $\operatorname{proj}_W\mathbf{x}=\mathbf{x}-\operatorname{proj}_{(1,1,1)}\mathbf{x}=(3,1,2)-(2,2,2)=(1,-1,0)$ ✓. When the complement is smaller than the subspace, projecting onto the complement and subtracting is quicker.

### Projection in function spaces and Bessel's inequality

Everything above works in a function space with the integral inner product. If $\{\phi_1,\dots,\phi_N\}$ is an orthogonal set of functions, the best approximation of $f$ in their span (minimum $\int\lvert f-\hat f\rvert^2$) is

$$
\hat f=\sum_{n=1}^Nc_n\phi_n,\qquad c_n=\frac{\langle f,\phi_n\rangle}{\langle\phi_n,\phi_n\rangle},
$$

and Pythagoras gives the error energy and **Bessel's inequality**:

$$
\lVert f-\hat f\rVert^2=\lVert f\rVert^2-\sum_{n=1}^N\lvert c_n\rvert^2\lVert\phi_n\rVert^2\ \ \ge0\quad\Longrightarrow\quad\sum_{n=1}^N\lvert c_n\rvert^2\lVert\phi_n\rVert^2\le\lVert f\rVert^2 .
$$

Two important consequences:

1. **The coefficients do not change when $N$ grows.** Adding $\phi_{N+1}$ to the basis adds a term but leaves $c_1,\dots,c_N$ alone — a property unique to orthogonal bases, and the reason Fourier coefficients can be computed one at a time.
2. **The error energy decreases monotonically with $N$.** If it tends to zero for every $f$, the orthogonal set is **complete** and Bessel's inequality becomes **Parseval's identity** $\sum\lvert c_n\rvert^2\lVert\phi_n\rVert^2=\lVert f\rVert^2$ — the central convergence statement of [Ch. 15](../part-5-transforms/15-fourier-series.md).

A simple but important special case: the projection of $f$ onto the constant function $1$ on $[a,b]$ is $\langle f,1\rangle/\langle1,1\rangle=\frac1{b-a}\int_a^bf\,dt$, the **average value**. The DC component of a signal is its projection onto the constants.

> [!example] Worked Example 9.31 — Best sinusoidal approximation of a square wave
> Let $q(t)=\operatorname{sgn}(\sin t)$ on $[-\pi,\pi]$ (a $\pm1$ square wave, period $2\pi$). Find its best approximation (least error energy) in $\operatorname{span}\{\sin t\}$ and then in $\operatorname{span}\{\sin t,\sin3t\}$, and the fraction of energy left in the error.
>
> **Inner products.** $\lVert q\rVert^2=\int_{-\pi}^{\pi}1\,dt=2\pi$, $\lVert\sin nt\rVert^2=\pi$, and since $q(t)\sin nt$ is even,
> $$
> \langle q,\sin nt\rangle=2\int_0^\pi\sin nt\,dt=\frac{2(1-\cos n\pi)}{n}=\begin{cases}4/n,&n\text{ odd}\\0,&n\text{ even.}\end{cases}
> $$
> **One term.** $c_1=\dfrac{4}{\pi}=1.2732$, so $\hat q_1(t)=\dfrac4\pi\sin t$. Note that the best amplitude is *larger* than the square wave's amplitude 1.
>
> Error energy: $\lVert q\rVert^2-c_1^2\lVert\sin t\rVert^2=2\pi-\dfrac{16}{\pi^2}\pi=2\pi-\dfrac{16}{\pi}=6.2832-5.0930=1.1902$. Fraction of energy in the error: $1-\dfrac{8}{\pi^2}=0.1894$ (18.9 %).
>
> **Two terms.** $\sin t\perp\sin3t$, so $c_1$ is unchanged and $c_3=\dfrac{4/3}{\pi}=\dfrac{4}{3\pi}=0.4244$. Error fraction:
> $$
> 1-\frac{8}{\pi^2}\Big(1+\frac19\Big)=0.0994\quad(9.9\ \%).
> $$
> Adding the third harmonic halves the error energy, and $c_1$ did not have to be recomputed.

```desmos-graph
left=-3.1416; right=3.1416; bottom=-1.6; top=1.6
---
y=\operatorname{sign}(\sin(x))|black
y=\frac{4}{\pi}\sin(x)|blue|dashed
y=\frac{4}{\pi}\left(\sin(x)+\frac{\sin(3x)}{3}\right)|red
```
*A $\pm1$ square wave (black), its orthogonal projection onto $\operatorname{span}\{\sin t\}$ (blue dashed, amplitude $4/\pi$) and onto $\operatorname{span}\{\sin t,\sin3t\}$ (red). Each added orthogonal direction reduces the error energy without altering the earlier coefficients.*

> [!tip] Engineering Note — the orthogonality principle in estimation
> The condition "error $\perp$ every allowed direction" is the **orthogonality principle** of estimation theory. In a Wiener filter or a linear minimum-mean-square-error estimator, the vector space consists of random variables with inner product $\langle X,Y\rangle=E[XY]$, and the optimal linear estimate makes the estimation error uncorrelated with every observation used. Only the inner product has changed; the geometry is that of this section ([Ch. 21](../part-6-probability-statistics/21-random-processes.md)).

## 9.10 Gram–Schmidt orthogonalization

### Motivation

Orthogonal bases make coordinates, projections and approximations trivial — but the bases that arise naturally (the columns of a data matrix, the monomials $1,t,t^2,\dots$, a set of measured waveforms) are rarely orthogonal. The **Gram–Schmidt process** converts any independent set into an orthonormal set with the same span, one vector at a time. In matrix form it is the **QR factorisation** $\mathbf{A}=\mathbf{Q}\mathbf{R}$, the standard numerical route to least squares; applied to polynomials it produces the **Legendre polynomials**, and applied to other weights it produces the other classical orthogonal polynomials used in filter design and spectral methods.

### The algorithm

Given independent $\mathbf{a}_1,\dots,\mathbf{a}_k$ in an inner-product space, subtract from each new vector its projections onto the orthogonal vectors already built:

$$
\begin{aligned}
\mathbf{v}_1&=\mathbf{a}_1,\\
\mathbf{v}_2&=\mathbf{a}_2-\frac{\langle\mathbf{a}_2,\mathbf{v}_1\rangle}{\langle\mathbf{v}_1,\mathbf{v}_1\rangle}\mathbf{v}_1,\\
\mathbf{v}_3&=\mathbf{a}_3-\frac{\langle\mathbf{a}_3,\mathbf{v}_1\rangle}{\langle\mathbf{v}_1,\mathbf{v}_1\rangle}\mathbf{v}_1-\frac{\langle\mathbf{a}_3,\mathbf{v}_2\rangle}{\langle\mathbf{v}_2,\mathbf{v}_2\rangle}\mathbf{v}_2,\\
&\ \,\vdots\\
\mathbf{v}_j&=\mathbf{a}_j-\sum_{i=1}^{j-1}\frac{\langle\mathbf{a}_j,\mathbf{v}_i\rangle}{\langle\mathbf{v}_i,\mathbf{v}_i\rangle}\mathbf{v}_i=\mathbf{a}_j-\operatorname{proj}_{W_{j-1}}\mathbf{a}_j,
\end{aligned}
$$

where $W_{j-1}=\operatorname{span}\{\mathbf{a}_1,\dots,\mathbf{a}_{j-1}\}$. Finally normalise: $\mathbf{q}_j=\mathbf{v}_j/\lVert\mathbf{v}_j\rVert$.

> [!theorem] Gram–Schmidt
> If $\mathbf{a}_1,\dots,\mathbf{a}_k$ are linearly independent, the vectors $\mathbf{q}_1,\dots,\mathbf{q}_k$ produced above are orthonormal and, for every $j$,
> $$
> \operatorname{span}\{\mathbf{q}_1,\dots,\mathbf{q}_j\}=\operatorname{span}\{\mathbf{a}_1,\dots,\mathbf{a}_j\}.
> $$
> Consequently every finite-dimensional inner-product space has an orthonormal basis.

*Proof by induction.* Suppose $\mathbf{v}_1,\dots,\mathbf{v}_{j-1}$ are orthogonal, nonzero and span $W_{j-1}$. Then $\mathbf{v}_j=\mathbf{a}_j-\operatorname{proj}_{W_{j-1}}\mathbf{a}_j$ is orthogonal to $W_{j-1}$ by the best-approximation theorem (the error of a projection is perpendicular to the subspace). It is nonzero, because otherwise $\mathbf{a}_j\in W_{j-1}$, contradicting independence. And $\mathbf{v}_j$ is $\mathbf{a}_j$ minus an element of $W_{j-1}$, so $\operatorname{span}\{\mathbf{v}_1,\dots,\mathbf{v}_j\}=\operatorname{span}\{\mathbf{a}_1,\dots,\mathbf{a}_j\}$.

If the input vectors are dependent, some $\mathbf{v}_j$ comes out as $\mathbf{0}$; dropping it and continuing yields an orthonormal basis of the span. In numerical work a "nearly zero" $\mathbf{v}_j$ is a warning of near-dependence.

> [!example] Worked Example 9.32 — Gram–Schmidt in $\mathbb{R}^3$
> Orthonormalise $\mathbf{a}_1=(1,1,0)$, $\mathbf{a}_2=(1,0,1)$, $\mathbf{a}_3=(0,1,1)$ (independent by Worked Example 9.7(b)).
>
> **Step 1.** $\mathbf{v}_1=(1,1,0)$, $\lVert\mathbf{v}_1\rVert^2=2$.
>
> **Step 2.** $\langle\mathbf{a}_2,\mathbf{v}_1\rangle=1$:
> $$
> \mathbf{v}_2=(1,0,1)-\tfrac12(1,1,0)=\big(\tfrac12,-\tfrac12,1\big),\qquad\lVert\mathbf{v}_2\rVert^2=\tfrac14+\tfrac14+1=\tfrac32 .
> $$
> **Step 3.** $\langle\mathbf{a}_3,\mathbf{v}_1\rangle=1$, $\langle\mathbf{a}_3,\mathbf{v}_2\rangle=-\tfrac12+1=\tfrac12$:
> $$
> \mathbf{v}_3=(0,1,1)-\tfrac12(1,1,0)-\frac{1/2}{3/2}\big(\tfrac12,-\tfrac12,1\big)=\big(-\tfrac12-\tfrac16,\ \tfrac12+\tfrac16,\ 1-\tfrac13\big)=\big(-\tfrac23,\tfrac23,\tfrac23\big),
> $$
> $\lVert\mathbf{v}_3\rVert^2=\tfrac43$.
>
> **Checks.** $\mathbf{v}_1\cdot\mathbf{v}_2=\tfrac12-\tfrac12+0=0$; $\mathbf{v}_1\cdot\mathbf{v}_3=-\tfrac23+\tfrac23=0$; $\mathbf{v}_2\cdot\mathbf{v}_3=-\tfrac13-\tfrac13+\tfrac23=0$ ✓.
>
> **Normalise.**
> $$
> \mathbf{q}_1=\tfrac1{\sqrt2}(1,1,0),\qquad\mathbf{q}_2=\tfrac1{\sqrt6}(1,-1,2),\qquad\mathbf{q}_3=\tfrac1{\sqrt3}(-1,1,1).
> $$
> (Scaling $\mathbf{v}_2$ by 2 and $\mathbf{v}_3$ by $\tfrac32$ before normalising avoids fractions — any positive rescaling of a $\mathbf{v}_j$ is harmless.)

### QR factorisation

Gram–Schmidt can be recorded as a matrix factorisation. Since each $\mathbf{a}_j$ is a combination of $\mathbf{q}_1,\dots,\mathbf{q}_j$ only,

$$
\mathbf{a}_j=r_{1j}\mathbf{q}_1+r_{2j}\mathbf{q}_2+\dots+r_{jj}\mathbf{q}_j,\qquad r_{ij}=\langle\mathbf{a}_j,\mathbf{q}_i\rangle\ (i<j),\quad r_{jj}=\lVert\mathbf{v}_j\rVert>0 .
$$

> [!theorem] QR factorisation
> Every $\mathbf{A}\in\mathbb{R}^{m\times n}$ with independent columns can be written
> $$
> \mathbf{A}=\mathbf{Q}\mathbf{R},
> $$
> where $\mathbf{Q}\in\mathbb{R}^{m\times n}$ has orthonormal columns ($\mathbf{Q}^\mathsf{T}\mathbf{Q}=\mathbf{I}_n$) spanning $\operatorname{col}(\mathbf{A})$, and $\mathbf{R}\in\mathbb{R}^{n\times n}$ is upper triangular with positive diagonal. With these sign conventions the factorisation is unique.

Since $\mathbf{Q}^\mathsf{T}\mathbf{Q}=\mathbf{I}$, the triangular factor is obtained as $\mathbf{R}=\mathbf{Q}^\mathsf{T}\mathbf{A}$, and the projection onto the column space is $\mathbf{P}=\mathbf{Q}\mathbf{Q}^\mathsf{T}$.

> [!example] Worked Example 9.33 — QR of a $3\times2$ matrix
> Factor $\mathbf{A}=\begin{bmatrix}1&0\\1&1\\1&2\end{bmatrix}$ (the matrix of Worked Example 9.29).
>
> $\mathbf{a}_1=(1,1,1)$: $r_{11}=\sqrt3$, $\mathbf{q}_1=\frac1{\sqrt3}(1,1,1)$.
>
> $\mathbf{a}_2=(0,1,2)$: $r_{12}=\langle\mathbf{a}_2,\mathbf{q}_1\rangle=\frac{3}{\sqrt3}=\sqrt3$. Then $\mathbf{v}_2=\mathbf{a}_2-r_{12}\mathbf{q}_1=(0,1,2)-(1,1,1)=(-1,0,1)$, $r_{22}=\sqrt2$, $\mathbf{q}_2=\frac1{\sqrt2}(-1,0,1)$.
> $$
> \mathbf{A}=\underbrace{\begin{bmatrix}1/\sqrt3&-1/\sqrt2\\1/\sqrt3&0\\1/\sqrt3&1/\sqrt2\end{bmatrix}}_{\mathbf{Q}}\underbrace{\begin{bmatrix}\sqrt3&\sqrt3\\0&\sqrt2\end{bmatrix}}_{\mathbf{R}}.
> $$
> **Check** the second column: $\sqrt3\,\mathbf{q}_1+\sqrt2\,\mathbf{q}_2=(1,1,1)+(-1,0,1)=(0,1,2)$ ✓. Notice that $\mathbf{v}_2=(-1,0,1)$ is the second column "centred" by removing its mean — Gram–Schmidt against the all-ones column is mean removal.

### Orthogonal polynomials: Legendre

Apply Gram–Schmidt to $1,t,t^2,t^3$ in $C[-1,1]$ with $\langle f,g\rangle=\int_{-1}^1fg\,dt$. Useful integrals: $\int_{-1}^1t^k\,dt=\frac{2}{k+1}$ for even $k$, $0$ for odd $k$.

- $v_0=1$, $\lVert v_0\rVert^2=2$.
- $v_1=t-\frac{\langle t,1\rangle}{2}\cdot1=t$ (since $\langle t,1\rangle=0$), $\lVert v_1\rVert^2=\frac23$.
- $v_2=t^2-\frac{\langle t^2,1\rangle}{2}\cdot1-\frac{\langle t^2,t\rangle}{2/3}t=t^2-\frac{2/3}{2}=t^2-\frac13$.
- $v_3=t^3-\frac{\langle t^3,1\rangle}{2}-\frac{\langle t^3,t\rangle}{2/3}t-\frac{\langle t^3,v_2\rangle}{\lVert v_2\rVert^2}v_2=t^3-\frac{2/5}{2/3}t=t^3-\frac35t$ (the first and last inner products vanish by odd symmetry).

Rescaled so that $P_n(1)=1$, these are the **Legendre polynomials**:

$$
P_0=1,\qquad P_1=t,\qquad P_2=\tfrac12(3t^2-1),\qquad P_3=\tfrac12(5t^3-3t),\qquad \int_{-1}^1P_mP_n\,dt=\frac{2}{2n+1}\delta_{mn}.
$$

They satisfy the three-term recurrence $(n+1)P_{n+1}=(2n+1)tP_n-nP_{n-1}$, which generates them far more cheaply than Gram–Schmidt. Changing the interval or inserting a weight function $w(t)$ in the inner product produces the other classical families: Chebyshev ($w=1/\sqrt{1-t^2}$ — used in Chebyshev filters and minimax approximation), Laguerre ($w=e^{-t}$ on $[0,\infty)$), Hermite ($w=e^{-t^2}$ on $\mathbb{R}$).

```desmos-graph
left=-1; right=1; bottom=-1.2; top=1.2
---
y=1|grey
y=x|blue
y=\frac{3x^{2}-1}{2}|red
y=\frac{5x^{3}-3x}{2}|green
```
*Legendre polynomials $P_0$ (grey), $P_1$ (blue), $P_2$ (red) and $P_3$ (green) on $[-1,1]$: each has exactly $n$ zeros in the interval and they are mutually orthogonal under $\int_{-1}^1fg\,dt$.*

> [!example] Worked Example 9.34 — Checking Legendre orthogonality and norms
> Verify $\langle P_0,P_2\rangle=0$, $\langle P_1,P_3\rangle=0$ and $\lVert P_2\rVert^2=\frac25$.
>
> $$
> \langle P_0,P_2\rangle=\tfrac12\int_{-1}^1(3t^2-1)\,dt=\tfrac12(2-2)=0\ ✓
> $$
> $$
> \langle P_1,P_3\rangle=\tfrac12\int_{-1}^1(5t^4-3t^2)\,dt=\tfrac12\big(5\cdot\tfrac25-3\cdot\tfrac23\big)=\tfrac12(2-2)=0\ ✓
> $$
> $$
> \lVert P_2\rVert^2=\tfrac14\int_{-1}^1(9t^4-6t^2+1)\,dt=\tfrac14\big(\tfrac{18}{5}-4+2\big)=\tfrac14\cdot\tfrac85=\tfrac25=\tfrac{2}{2\cdot2+1}\ ✓
> $$
> (Pairs of opposite parity, such as $P_0,P_1$ or $P_2,P_3$, are orthogonal automatically by symmetry.)

### Numerical behaviour: classical versus modified Gram–Schmidt

In floating-point arithmetic the formula above (**classical** Gram–Schmidt, CGS) computes all the coefficients $\langle\mathbf{a}_j,\mathbf{v}_i\rangle$ from the *original* vector $\mathbf{a}_j$. When the columns are nearly dependent, rounding errors make the computed $\mathbf{q}_j$ lose orthogonality badly. **Modified** Gram–Schmidt (MGS) subtracts the projections one at a time, using the partially updated vector:

$$
\mathbf{w}\leftarrow\mathbf{a}_j;\quad\text{for }i=1,\dots,j-1:\ \ r_{ij}=\langle\mathbf{w},\mathbf{q}_i\rangle,\ \ \mathbf{w}\leftarrow\mathbf{w}-r_{ij}\mathbf{q}_i;\quad r_{jj}=\lVert\mathbf{w}\rVert,\ \ \mathbf{q}_j=\mathbf{w}/r_{jj}.
$$

In exact arithmetic CGS and MGS are identical; in floating point MGS loses orthogonality roughly in proportion to $\kappa(\mathbf{A})\,\varepsilon_{\text{mach}}$ rather than $\kappa(\mathbf{A})^2\varepsilon_{\text{mach}}$. Library QR routines use Householder reflections, which are more stable still ([Ch. 22](../part-7-numerical-computational/22-numerical-methods.md)).

> [!warning] Common Pitfall — projecting onto the original vectors
> A frequent hand-calculation error is to subtract projections onto the *original* $\mathbf{a}_i$ instead of onto the already-orthogonalised $\mathbf{v}_i$. The $\mathbf{a}_i$ are not orthogonal, so the sum-of-projections formula does not apply to them and the result is not orthogonal. Always project onto the $\mathbf{v}_i$ (or $\mathbf{q}_i$), and check each new vector against all previous ones.

> [!tip] Engineering Note — orthogonalising measured signals
> Gram–Schmidt is used directly in adaptive signal processing: a **Gram–Schmidt (lattice) canceller** in a radar or antenna array orthogonalises the auxiliary channels one stage at a time, so that interference can be removed from the main channel stage by stage. The same idea appears as "prewhitening" and "decorrelation" in multichannel receivers, and in the QR-based recursive least-squares (QR-RLS) adaptive filter.

## 9.11 Least-squares approximation

### Motivation

Measurements are noisy and models are idealised, so an engineer usually has **more equations than unknowns** and no exact solution. Calibrating a sensor from twenty readings with a two-parameter straight-line model, estimating a resistor's value and an amplifier's offset from many voltage–current pairs, identifying the time constant of an RC network from a sampled decay, fitting a polynomial to a thermistor curve — each is an overdetermined linear system $\mathbf{A}\mathbf{x}=\mathbf{b}$ with $\mathbf{b}\notin\operatorname{col}(\mathbf{A})$. The **least-squares** solution chooses $\mathbf{x}$ to make the residual as small as possible in the Euclidean norm. Geometrically it is a projection (Section 9.9); computationally it is solved by the normal equations or, better, by QR (Section 9.10).

### The least-squares problem

> [!definition] Least-squares solution
> For $\mathbf{A}\in\mathbb{R}^{m\times n}$ and $\mathbf{b}\in\mathbb{R}^m$, a **least-squares solution** of $\mathbf{A}\mathbf{x}=\mathbf{b}$ is any $\hat{\mathbf{x}}$ minimising
> $$
> J(\mathbf{x})=\lVert\mathbf{b}-\mathbf{A}\mathbf{x}\rVert^2=\sum_{i=1}^m\big(b_i-(\mathbf{A}\mathbf{x})_i\big)^2 .
> $$
> The vector $\mathbf{r}=\mathbf{b}-\mathbf{A}\hat{\mathbf{x}}$ is the **residual**, and $\lVert\mathbf{r}\rVert^2$ is the **sum of squared errors** (SSE).

> [!theorem] Normal equations
> $\hat{\mathbf{x}}$ is a least-squares solution of $\mathbf{A}\mathbf{x}=\mathbf{b}$ if and only if
> $$
> \mathbf{A}^\mathsf{T}\mathbf{A}\,\hat{\mathbf{x}}=\mathbf{A}^\mathsf{T}\mathbf{b}.
> $$
> The normal equations are always consistent. The least-squares solution is unique iff the columns of $\mathbf{A}$ are independent ($\operatorname{rank}\mathbf{A}=n$), and then
> $$
> \hat{\mathbf{x}}=(\mathbf{A}^\mathsf{T}\mathbf{A})^{-1}\mathbf{A}^\mathsf{T}\mathbf{b},\qquad \mathbf{A}\hat{\mathbf{x}}=\operatorname{proj}_{\operatorname{col}(\mathbf{A})}\mathbf{b}.
> $$

**Geometric derivation.** As $\mathbf{x}$ ranges over $\mathbb{R}^n$, $\mathbf{A}\mathbf{x}$ ranges over $\operatorname{col}(\mathbf{A})$. The closest point of that subspace to $\mathbf{b}$ is the orthogonal projection $\mathbf{p}$, characterised by $\mathbf{b}-\mathbf{p}\perp\operatorname{col}(\mathbf{A})$, i.e. $\mathbf{b}-\mathbf{A}\hat{\mathbf{x}}\in\operatorname{col}(\mathbf{A})^\perp=\operatorname{null}(\mathbf{A}^\mathsf{T})$ (Section 9.7). That is exactly $\mathbf{A}^\mathsf{T}(\mathbf{b}-\mathbf{A}\hat{\mathbf{x}})=\mathbf{0}$. Since $\mathbf{p}$ always exists, the normal equations are always consistent; if $\operatorname{null}(\mathbf{A})\neq\{\mathbf{0}\}$, any null-space vector can be added to $\hat{\mathbf{x}}$ without changing $\mathbf{A}\hat{\mathbf{x}}$.

**Calculus derivation.** Expand
$$
J(\mathbf{x})=\mathbf{b}^\mathsf{T}\mathbf{b}-2\mathbf{x}^\mathsf{T}\mathbf{A}^\mathsf{T}\mathbf{b}+\mathbf{x}^\mathsf{T}\mathbf{A}^\mathsf{T}\mathbf{A}\mathbf{x}.
$$
Setting the gradient ([Ch. 6](../part-2-calculus/06-multivariable-calculus.md)) to zero gives $\nabla J=-2\mathbf{A}^\mathsf{T}\mathbf{b}+2\mathbf{A}^\mathsf{T}\mathbf{A}\mathbf{x}=\mathbf{0}$, the normal equations. The Hessian $2\mathbf{A}^\mathsf{T}\mathbf{A}$ is positive semidefinite ($\mathbf{x}^\mathsf{T}\mathbf{A}^\mathsf{T}\mathbf{A}\mathbf{x}=\lVert\mathbf{A}\mathbf{x}\rVert^2\ge0$), so the stationary point is a minimum — strictly so when the columns are independent. The two derivations agree, as they must: minimising a squared distance *is* projecting.

The minimum value follows from Pythagoras, $\lVert\mathbf{b}\rVert^2=\lVert\mathbf{A}\hat{\mathbf{x}}\rVert^2+\lVert\mathbf{r}\rVert^2$, and can be computed as $\lVert\mathbf{r}\rVert^2=\mathbf{b}^\mathsf{T}\mathbf{b}-\hat{\mathbf{x}}^\mathsf{T}\mathbf{A}^\mathsf{T}\mathbf{b}$.

### Fitting a straight line

Given data $(t_i,y_i)$, $i=1,\dots,m$, the model $y\approx c_0+c_1t$ produces one equation per data point:

$$
\underbrace{\begin{bmatrix}1&t_1\\1&t_2\\ \vdots&\vdots\\1&t_m\end{bmatrix}}_{\mathbf{A}}\begin{bmatrix}c_0\\c_1\end{bmatrix}\approx\underbrace{\begin{bmatrix}y_1\\y_2\\ \vdots\\y_m\end{bmatrix}}_{\mathbf{b}},\qquad
\mathbf{A}^\mathsf{T}\mathbf{A}=\begin{bmatrix}m&\sum t_i\\ \sum t_i&\sum t_i^2\end{bmatrix},\quad\mathbf{A}^\mathsf{T}\mathbf{b}=\begin{bmatrix}\sum y_i\\ \sum t_iy_i\end{bmatrix}.
$$

Solving the $2\times2$ system gives the familiar regression formulas

$$
c_1=\frac{\sum(t_i-\bar t)(y_i-\bar y)}{\sum(t_i-\bar t)^2},\qquad c_0=\bar y-c_1\bar t,
$$

where $\bar t,\bar y$ are the sample means. (They follow most quickly from the QR remark of Worked Example 9.33: orthogonalising the $t$ column against the ones column replaces $t_i$ by $t_i-\bar t$, which decouples the two unknowns.) In particular **the least-squares line passes through the centroid $(\bar t,\bar y)$**, and the residuals sum to zero (the residual is orthogonal to the column of ones).

> [!example] Worked Example 9.35 — The smallest line-fitting problem
> Fit $y=c_0+c_1t$ to the points $(0,6)$, $(1,0)$, $(2,0)$.
>
> This is exactly Worked Example 9.29 with $\mathbf{A}=\begin{bmatrix}1&0\\1&1\\1&2\end{bmatrix}$ and $\mathbf{b}=(6,0,0)^\mathsf{T}$. The normal equations gave $c_0=5$, $c_1=-3$:
> $$
> \hat y=5-3t,\qquad\text{fitted values }(5,2,-1),\qquad\text{residuals }(1,-2,1),\qquad\text{SSE}=6 .
> $$
> **Checks.** Residuals sum to zero ✓. The line passes through the centroid $(\bar t,\bar y)=(1,2)$: $5-3=2$ ✓. The regression formula: $\bar t=1$, $\sum(t_i-\bar t)(y_i-\bar y)=(-1)(4)+0+(1)(-2)=-6$, $\sum(t_i-\bar t)^2=2$, so $c_1=-3$ ✓.
>
> **Via QR** (Worked Example 9.33): $\mathbf{Q}^\mathsf{T}\mathbf{b}=\big(\tfrac{6}{\sqrt3},-\tfrac{6}{\sqrt2}\big)=(2\sqrt3,-3\sqrt2)$. Back-substitute in $\mathbf{R}\hat{\mathbf{x}}=\mathbf{Q}^\mathsf{T}\mathbf{b}$: $\sqrt2\,c_1=-3\sqrt2\Rightarrow c_1=-3$; $\sqrt3\,c_0+\sqrt3(-3)=2\sqrt3\Rightarrow c_0=5$ ✓.

> [!example] Worked Example 9.36 — Resistance and offset from V–I measurements
> A resistor is measured with a voltmeter that may have an offset. The data are
>
> | $I$ (mA) | 0 | 1 | 2 | 3 | 4 |
> |---|---|---|---|---|---|
> | $V$ (V) | 0.1 | 1.9 | 4.2 | 5.8 | 8.1 |
>
> Fit $V=V_{\text{off}}+RI$.
>
> **Sums.** $m=5$, $\sum I=10$, $\sum I^2=30$, $\sum V=20.1$, $\sum IV=0+1.9+8.4+17.4+32.4=60.1$.
>
> **Normal equations.**
> $$
> \begin{bmatrix}5&10\\10&30\end{bmatrix}\begin{bmatrix}V_{\text{off}}\\R\end{bmatrix}=\begin{bmatrix}20.1\\60.1\end{bmatrix}.
> $$
> Determinant $5\cdot30-10\cdot10=50$. By Cramer's rule,
> $$
> V_{\text{off}}=\frac{30(20.1)-10(60.1)}{50}=\frac{603-601}{50}=0.04\ \text{V},\qquad R=\frac{5(60.1)-10(20.1)}{50}=\frac{300.5-201}{50}=1.99\ \text{V/mA}=1.99\ \text{k}\Omega .
> $$
> **Residuals.** Fitted values $0.04,2.03,4.02,6.01,8.00$; residuals $0.06,-0.13,0.18,-0.21,0.10$ (sum $0$ ✓); SSE $=0.107\ \text{V}^2$.
>
> **Quality of fit.** The total variation about the mean is $\sum(V_i-\bar V)^2=\sum V_i^2-m\bar V^2=120.51-5(4.02)^2=39.71$, so the coefficient of determination is $R^2=1-0.107/39.71=0.9973$. **Units check:** $V_{\text{off}}$ in volts, $R$ in V/mA — each column of $\mathbf{A}$ carries its own units, and each unknown takes the units that make $\mathbf{A}\mathbf{x}$ volts.

```desmos-graph
left=-0.5; right=4.5; bottom=-0.5; top=9
---
y=0.04+1.99x|blue
(0,0.1)|black
(1,1.9)|black
(2,4.2)|black
(3,5.8)|black
(4,8.1)|black
(2,4.02)|label:centroid
```
*Least-squares line $V=0.04+1.99I$ through the five measured points; it passes through the centroid $(\bar I,\bar V)=(2,4.02)$, and the vertical residuals are what is minimised.*

### Polynomial and general linear-in-the-parameters models

The same machinery fits *any* model that is **linear in its unknown parameters**, however nonlinear it is in $t$:

$$
y(t)\approx c_1\phi_1(t)+c_2\phi_2(t)+\dots+c_n\phi_n(t),\qquad \mathbf{A}=\big[\phi_j(t_i)\big]_{m\times n}.
$$

Examples: polynomials ($\phi_j=t^{j-1}$, giving a **Vandermonde** matrix), sinusoids of known frequency ($\phi=\cos\omega t,\sin\omega t$, which estimates amplitude and phase), exponentials with known rates, splines. The column $j$ of $\mathbf{A}$ is simply the $j$th basis function sampled at the measurement times.

> [!example] Worked Example 9.37 — A quadratic fit
> Fit $y=c_0+c_1t+c_2t^2$ to
>
> | $t$ | $-2$ | $-1$ | $0$ | $1$ | $2$ |
> |---|---|---|---|---|---|
> | $y$ | 8.2 | 2.4 | 0.9 | 3.6 | 10.1 |
>
> **Matrix.** Rows $(1,t_i,t_i^2)$. Because the $t_i$ are symmetric about 0, all odd power sums vanish: $\sum t_i=\sum t_i^3=0$, $\sum t_i^2=10$, $\sum t_i^4=34$.
> $$
> \mathbf{A}^\mathsf{T}\mathbf{A}=\begin{bmatrix}5&0&10\\0&10&0\\10&0&34\end{bmatrix},\qquad\mathbf{A}^\mathsf{T}\mathbf{b}=\begin{bmatrix}\sum y_i\\ \sum t_iy_i\\ \sum t_i^2y_i\end{bmatrix}=\begin{bmatrix}25.2\\5.0\\79.2\end{bmatrix}.
> $$
> ($\sum t_iy_i=-16.4-2.4+0+3.6+20.2=5.0$; $\sum t_i^2y_i=32.8+2.4+0+3.6+40.4=79.2$.)
>
> **Solve.** The middle equation decouples: $10c_1=5\Rightarrow c_1=0.5$. From the first, $c_0=5.04-2c_2$; substitute in the third: $10(5.04-2c_2)+34c_2=79.2\Rightarrow14c_2=28.8\Rightarrow c_2=2.0571$, and $c_0=0.9257$.
> $$
> \hat y=0.9257+0.5t+2.0571t^2 .
> $$
> **Residuals** $(0.046,-0.083,-0.026,0.117,-0.054)$, SSE $=0.0263$. Symmetric sample points made $t$ orthogonal to both $1$ and $t^2$ — a small taste of the simplification orthogonal bases bring.

> [!example] Worked Example 9.38 — Time constant of an RC discharge (linearised fit)
> A capacitor discharging through a resistor is sampled:
>
> | $t$ (ms) | 0 | 1 | 2 | 3 | 4 |
> |---|---|---|---|---|---|
> | $v$ (V) | 5.00 | 3.02 | 1.85 | 1.10 | 0.68 |
>
> Estimate $V_0$ and $\tau$ in $v=V_0e^{-t/\tau}$.
>
> **Linearise.** The model is not linear in $\tau$, but $\ln v=\ln V_0-t/\tau$ is linear in the parameters $a=\ln V_0$ and $b=-1/\tau$. Take logs: $\ln v=(1.6094,\ 1.1053,\ 0.6152,\ 0.0953,\ -0.3857)$.
>
> **Normal equations** (same $\mathbf{A}^\mathsf{T}\mathbf{A}$ as Worked Example 9.36 since the $t_i$ are again $0,\dots,4$):
> $$
> \begin{bmatrix}5&10\\10&30\end{bmatrix}\begin{bmatrix}a\\b\end{bmatrix}=\begin{bmatrix}\sum\ln v_i\\ \sum t_i\ln v_i\end{bmatrix}=\begin{bmatrix}3.0395\\1.0789\end{bmatrix}.
> $$
> $$
> a=\frac{30(3.0395)-10(1.0789)}{50}=1.6079,\qquad b=\frac{5(1.0789)-10(3.0395)}{50}=-0.5000\ \text{ms}^{-1}.
> $$
> **Result.** $V_0=e^{1.6079}=4.99$ V and $\tau=-1/b=2.00$ ms. Fitted values $4.99,3.03,1.84,1.11,0.68$ V agree with the data to within $0.02$ V.
>
> **Caveat.** Minimising squared errors in $\ln v$ is *not* the same as minimising squared errors in $v$: the log transform gives small voltages (late samples, where relative noise is largest) as much weight as large ones. When that matters, use weighted least squares (below) with weights $w_i\propto v_i^2$, or a genuinely nonlinear least-squares fit ([Ch. 23](../part-7-numerical-computational/23-optimization.md)).

### Weighted least squares

If some measurements are more reliable than others — say measurement $i$ has noise standard deviation $\sigma_i$ — minimise the weighted sum

$$
J_\mathbf{W}(\mathbf{x})=\sum_iw_i\big(b_i-(\mathbf{A}\mathbf{x})_i\big)^2=(\mathbf{b}-\mathbf{A}\mathbf{x})^\mathsf{T}\mathbf{W}(\mathbf{b}-\mathbf{A}\mathbf{x}),\qquad\mathbf{W}=\operatorname{diag}(w_i),\ w_i=1/\sigma_i^2 .
$$

This is ordinary least squares in the weighted inner product $\langle\mathbf{u},\mathbf{v}\rangle_\mathbf{W}=\mathbf{v}^\mathsf{T}\mathbf{W}\mathbf{u}$ of Section 9.8, so the orthogonality principle gives the **weighted normal equations**

$$
\mathbf{A}^\mathsf{T}\mathbf{W}\mathbf{A}\,\hat{\mathbf{x}}=\mathbf{A}^\mathsf{T}\mathbf{W}\mathbf{b}.
$$

Equivalently, scale row $i$ of $\mathbf{A}$ and $\mathbf{b}$ by $\sqrt{w_i}$ and solve an ordinary problem. The statistical justification (best linear unbiased estimation) and confidence intervals for the fitted parameters belong to [Ch. 20](../part-6-probability-statistics/20-statistical-analysis.md).

> [!example] Worked Example 9.39 — Averaging measurements of unequal quality
> A voltage is measured three times: $5.02$ V and $4.98$ V with a meter of $\sigma=0.02$ V, and $5.20$ V with a meter of $\sigma=0.10$ V. Find the weighted least-squares estimate.
>
> The model is $v\approx x$ for each reading: $\mathbf{A}=(1,1,1)^\mathsf{T}$, $\mathbf{b}=(5.02,4.98,5.20)^\mathsf{T}$, weights $w=(2500,2500,100)$ V$^{-2}$.
> $$
> \hat x=\frac{\mathbf{A}^\mathsf{T}\mathbf{W}\mathbf{b}}{\mathbf{A}^\mathsf{T}\mathbf{W}\mathbf{A}}=\frac{2500(5.02)+2500(4.98)+100(5.20)}{2500+2500+100}=\frac{25\,520}{5100}=5.0039\ \text{V}.
> $$
> The unweighted mean would be $5.0667$ V, pulled up by the poor reading. Weighted least squares with a single constant column is simply the **inverse-variance-weighted mean**.

### Numerical aspects: normal equations versus QR

Forming $\mathbf{A}^\mathsf{T}\mathbf{A}$ squares the condition number:

$$
\kappa(\mathbf{A}^\mathsf{T}\mathbf{A})=\kappa(\mathbf{A})^2 .
$$

If $\kappa(\mathbf{A})=10^6$, the normal equations have $\kappa=10^{12}$ and, in double precision, only about four correct digits can be expected. The QR route avoids the squaring. With $\mathbf{A}=\mathbf{Q}\mathbf{R}$,

$$
\mathbf{A}^\mathsf{T}\mathbf{A}\hat{\mathbf{x}}=\mathbf{A}^\mathsf{T}\mathbf{b}\ \Longrightarrow\ \mathbf{R}^\mathsf{T}\mathbf{Q}^\mathsf{T}\mathbf{Q}\mathbf{R}\hat{\mathbf{x}}=\mathbf{R}^\mathsf{T}\mathbf{Q}^\mathsf{T}\mathbf{b}\ \Longrightarrow\ \mathbf{R}\hat{\mathbf{x}}=\mathbf{Q}^\mathsf{T}\mathbf{b},
$$

a triangular system solved by back substitution (Worked Example 9.35). This is what `numpy.linalg.lstsq`, MATLAB's backslash and similar routines do (some use the SVD of [Ch. 10](10-eigenvalues-and-eigenvectors.md), which also handles rank deficiency gracefully).

The classic ill-conditioned case is polynomial fitting with the monomial basis. In the continuous version — approximating $f$ on $[0,1]$ by $\sum_{j=0}^{n-1}c_jt^j$ — the normal-equation matrix has entries $\langle t^i,t^j\rangle=\int_0^1t^{i+j}dt=\frac{1}{i+j+1}$: the **Hilbert matrix**. Its condition number explodes:

| $n$ | 3 | 6 | 10 |
|---|---|---|---|
| $\kappa(\mathbf{H}_n)$ | $5.2\times10^{2}$ | $1.5\times10^{7}$ | $1.6\times10^{13}$ |

The cause is the near-parallelism of $t^i$ and $t^j$ seen in Worked Example 9.25(b). The cure is to use an orthogonal basis, which makes the normal-equation matrix diagonal.

### Continuous least squares and orthogonal polynomials

To approximate a function $f$ on $[a,b]$ by an element of a finite-dimensional subspace $W$ in the sense of least integrated squared error, project in the function-space inner product (Section 9.9). With an orthogonal basis $\{\phi_k\}$ of $W$ there is no system to solve at all:

$$
\hat f=\sum_kc_k\phi_k,\qquad c_k=\frac{\langle f,\phi_k\rangle}{\langle\phi_k,\phi_k\rangle},\qquad\lVert f-\hat f\rVert^2=\lVert f\rVert^2-\sum_k\lvert c_k\rvert^2\lVert\phi_k\rVert^2 .
$$

> [!example] Worked Example 9.40 — Best quadratic approximation of $e^t$ on $[-1,1]$
> Find the quadratic $q(t)$ minimising $\int_{-1}^1(e^t-q(t))^2\,dt$ and compare it with the Taylor polynomial $1+t+t^2/2$.
>
> **Use Legendre polynomials** (Section 9.10), $\lVert P_k\rVert^2=\frac{2}{2k+1}$, so $c_k=\frac{2k+1}{2}\int_{-1}^1e^tP_k\,dt$. The integrals ($\int te^tdt=(t-1)e^t$, $\int t^2e^tdt=(t^2-2t+2)e^t$):
> $$
> \int_{-1}^1e^t\,dt=e-e^{-1},\qquad\int_{-1}^1te^t\,dt=\frac{2}{e},\qquad\int_{-1}^1t^2e^t\,dt=e-\frac5e .
> $$
> Hence
> $$
> c_0=\tfrac12(e-e^{-1})=\sinh1=1.1752,\qquad c_1=\tfrac32\cdot\tfrac2e=\tfrac3e=1.1036,
> $$
> $$
> c_2=\tfrac52\int_{-1}^1\tfrac{3t^2-1}{2}e^t\,dt=\tfrac52\Big[\tfrac32\big(e-\tfrac5e\big)-\tfrac12\big(e-\tfrac1e\big)\Big]=\tfrac52\Big(e-\tfrac7e\Big)=0.3578 .
> $$
> **Assemble.** $q=c_0+c_1t+c_2\cdot\frac{3t^2-1}{2}$:
> $$
> q(t)=\big(1.1752-0.1789\big)+1.1036\,t+0.5367\,t^2=0.9963+1.1036\,t+0.5367\,t^2 .
> $$
> **Error.** Numerically, $\int_{-1}^1(e^t-q)^2dt=1.44\times10^{-3}$ (RMS error $0.027$ over the interval), whereas the Taylor polynomial gives $8.99\times10^{-3}$ — more than six times larger. Taylor is perfect at $t=0$ and poor at the ends; least squares spreads the error over the whole interval (errors $-0.061$, $0.004$, $0.082$ at $t=-1,0,1$).
>
> **Monomial check.** Solving the $3\times3$ normal equations in the basis $1,t,t^2$ (Gram matrix entries $\int_{-1}^1t^{i+j}dt$) gives the same $q$, as it must — the projection does not depend on the basis — but requires solving a coupled system.

```desmos-graph
left=-1; right=1; bottom=0.2; top=2.9
---
y=e^{x}|black
y=0.9963+1.1036x+0.5367x^{2}|blue|dashed
y=1+x+\frac{x^{2}}{2}|red|dotted
```
*$e^t$ (black), its least-squares quadratic on $[-1,1]$ (blue dashed, nearly indistinguishable) and the Taylor quadratic about $0$ (red dotted), which drifts away towards $t=\pm1$.*

### Fourier series as least squares

Take $W=\operatorname{span}\{1,\cos\omega_0t,\sin\omega_0t,\dots,\cos N\omega_0t,\sin N\omega_0t\}$ on one period. This basis is orthogonal (Section 9.8), so the best approximation of a periodic $f$ in $W$ — the trigonometric polynomial of order $N$ with least mean-square error — has coefficients

$$
a_0=\frac1T\int_Tf\,dt,\qquad a_k=\frac2T\int_Tf\cos k\omega_0t\,dt,\qquad b_k=\frac2T\int_Tf\sin k\omega_0t\,dt .
$$

These are precisely the Fourier coefficients. In other words, **the truncated Fourier series is the least-squares approximation of $f$ by harmonics up to order $N$**, the error is orthogonal to every retained harmonic, and adding more harmonics never changes the earlier coefficients (Worked Example 9.31). In discrete form, fitting $N$ samples with the DFT basis is a least-squares problem whose matrix has orthogonal columns, which is why the DFT can be inverted by a conjugate transpose. [Ch. 15](../part-5-transforms/15-fourier-series.md) takes up convergence as $N\to\infty$, Gibbs' phenomenon and power relations; [Ch. 16](../part-5-transforms/16-fourier-transforms.md) extends the idea to aperiodic signals.

> [!warning] Common Pitfall — extrapolating and over-fitting
> A least-squares polynomial is only as good as the range of the data. A high-degree polynomial fitted to $m$ points can pass close to every one of them while oscillating wildly between and beyond them; with degree $m-1$ it interpolates exactly, the residual is zero, and noise is fitted along with signal. Choose the model order from physics or from residual behaviour on held-out data, not by driving the SSE to zero, and never trust a polynomial fit outside the measured interval.

> [!tip] Engineering Note — where least squares runs in hardware
> Least squares is computed continuously inside modern electronics: channel estimation in a Wi-Fi or LTE receiver solves a small least-squares problem against known pilot symbols on every frame; GPS receivers solve linearised least squares for position and clock bias from four or more pseudoranges; battery-management ICs fit equivalent-circuit parameters; adaptive equalisers run recursive least squares (RLS), which updates $\hat{\mathbf{x}}$ as each new row of $\mathbf{A}$ arrives without re-solving from scratch. The mathematics is always this section's: project the measurement vector onto the column space of the model.

## Summary

| Concept | Key formula / statement |
|---|---|
| Vector space | closed under $\mathbf{u}+\mathbf{v}$ and $c\mathbf{u}$; eight axioms; contains $\mathbf{0}$ |
| Subspace test | $\mathbf{0}\in W$ and $a\mathbf{u}+b\mathbf{v}\in W$ for all $\mathbf{u},\mathbf{v}\in W$ |
| Span | $\operatorname{span}\{\mathbf{v}_i\}=\{\sum c_i\mathbf{v}_i\}$, the smallest subspace containing the $\mathbf{v}_i$ |
| Linear independence | $\sum c_i\mathbf{v}_i=\mathbf{0}\Rightarrow$ all $c_i=0$; in $\mathbb{F}^n$: $\operatorname{rank}[\mathbf{v}_1\cdots\mathbf{v}_k]=k$ |
| Wronskian | $W(f_1,\dots,f_k)(t_0)\neq0\Rightarrow$ independent |
| Basis, coordinates | independent + spanning; $\mathbf{x}=\mathbf{P}_\mathcal{B}[\mathbf{x}]_\mathcal{B}$ |
| Dimension | number of vectors in any basis; $\dim\mathcal{P}_n=n+1$, $\dim\mathbb{R}^{m\times n}=mn$ |
| Linear map | $T(a\mathbf{u}+b\mathbf{v})=aT(\mathbf{u})+bT(\mathbf{v})$; matrix columns $=T(\mathbf{e}_j)$ |
| Change of basis | $[T]_\mathcal{B}=\mathbf{P}^{-1}\mathbf{A}\mathbf{P}$ (similar matrices) |
| Rank–nullity | $\operatorname{rank}\mathbf{A}+\dim\operatorname{null}\mathbf{A}=n$ |
| Four subspaces | $\dim\operatorname{col}=\dim\operatorname{row}=r$, $\dim\operatorname{null}=n-r$, $\dim\operatorname{null}(\mathbf{A}^\mathsf{T})=m-r$ |
| Orthogonal complements | $\operatorname{null}(\mathbf{A})=\operatorname{row}(\mathbf{A})^\perp$, $\operatorname{null}(\mathbf{A}^\mathsf{T})=\operatorname{col}(\mathbf{A})^\perp$ |
| Solution structure | $\mathbf{x}=\mathbf{x}_p+\mathbf{x}_n$, $\mathbf{x}_n\in\operatorname{null}(\mathbf{A})$; solvable iff $\mathbf{b}\in\operatorname{col}(\mathbf{A})$ |
| Cauchy–Schwarz | $\lvert\langle\mathbf{u},\mathbf{v}\rangle\rvert\le\lVert\mathbf{u}\rVert\lVert\mathbf{v}\rVert$; $\cos\theta=\langle\mathbf{u},\mathbf{v}\rangle/(\lVert\mathbf{u}\rVert\lVert\mathbf{v}\rVert)$ |
| Complex inner product | $\langle\mathbf{u},\mathbf{v}\rangle=\mathbf{v}^\mathsf{H}\mathbf{u}=\sum u_kv_k^*$ |
| Signal inner products | $\langle x,y\rangle=\int x\,y^*\,dt$ (energy); $\langle x,y\rangle_T=\frac1T\int_Tx\,y^*\,dt$ (power, norm $=$ RMS) |
| Orthogonal harmonics | $\frac1T\int_Te^{jm\omega_0t}e^{-jn\omega_0t}dt=\delta_{mn}$ |
| Orthogonal coordinates | $c_j=\langle\mathbf{x},\mathbf{u}_j\rangle/\langle\mathbf{u}_j,\mathbf{u}_j\rangle$; Parseval $\lVert\mathbf{x}\rVert^2=\sum\lvert c_j\rvert^2\lVert\mathbf{u}_j\rVert^2$ |
| Orthogonal matrix | $\mathbf{Q}^\mathsf{T}\mathbf{Q}=\mathbf{I}$, $\mathbf{Q}^{-1}=\mathbf{Q}^\mathsf{T}$, $\lVert\mathbf{Q}\mathbf{x}\rVert=\lVert\mathbf{x}\rVert$ |
| Projection onto a line | $\operatorname{proj}_\mathbf{a}\mathbf{b}=\dfrac{\langle\mathbf{b},\mathbf{a}\rangle}{\langle\mathbf{a},\mathbf{a}\rangle}\mathbf{a}$ |
| Projection onto $\operatorname{col}(\mathbf{A})$ | $\mathbf{P}=\mathbf{A}(\mathbf{A}^\mathsf{T}\mathbf{A})^{-1}\mathbf{A}^\mathsf{T}=\mathbf{Q}\mathbf{Q}^\mathsf{T}$; $\mathbf{P}^2=\mathbf{P}=\mathbf{P}^\mathsf{T}$ |
| Best approximation | $\mathbf{b}-\operatorname{proj}_W\mathbf{b}\perp W$; Bessel $\sum\lvert c_k\rvert^2\lVert\phi_k\rVert^2\le\lVert f\rVert^2$ |
| Gram–Schmidt | $\mathbf{v}_j=\mathbf{a}_j-\sum_{i<j}\frac{\langle\mathbf{a}_j,\mathbf{v}_i\rangle}{\langle\mathbf{v}_i,\mathbf{v}_i\rangle}\mathbf{v}_i$; $\mathbf{A}=\mathbf{Q}\mathbf{R}$ |
| Legendre polynomials | $1,\ t,\ \frac12(3t^2-1),\ \frac12(5t^3-3t)$; $\int_{-1}^1P_mP_n\,dt=\frac{2}{2n+1}\delta_{mn}$ |
| Normal equations | $\mathbf{A}^\mathsf{T}\mathbf{A}\hat{\mathbf{x}}=\mathbf{A}^\mathsf{T}\mathbf{b}$; via QR: $\mathbf{R}\hat{\mathbf{x}}=\mathbf{Q}^\mathsf{T}\mathbf{b}$ |
| Line fit | $c_1=\dfrac{\sum(t_i-\bar t)(y_i-\bar y)}{\sum(t_i-\bar t)^2}$, $c_0=\bar y-c_1\bar t$ |
| Weighted LS | $\mathbf{A}^\mathsf{T}\mathbf{W}\mathbf{A}\hat{\mathbf{x}}=\mathbf{A}^\mathsf{T}\mathbf{W}\mathbf{b}$, $w_i=1/\sigma_i^2$ |
| Conditioning | $\kappa(\mathbf{A}^\mathsf{T}\mathbf{A})=\kappa(\mathbf{A})^2$ — prefer QR or SVD |

## Further reading

- G. Strang, *Introduction to Linear Algebra*, 6th ed., Wellesley-Cambridge Press, 2023 — the four fundamental subspaces, projections and least squares, with exceptionally clear geometry.
- D. C. Lay, S. R. Lay and J. J. McDonald, *Linear Algebra and Its Applications*, 6th ed., Pearson, 2021 — careful treatment of vector spaces, bases, coordinates and inner-product spaces.
- E. Kreyszig, *Advanced Engineering Mathematics*, 10th ed., Wiley, 2011 — linear algebra and orthogonal functions in the broader engineering-mathematics setting.
- T. K. Moon and W. C. Stirling, *Mathematical Methods and Algorithms for Signal Processing*, Prentice Hall, 2000 — vector spaces, Hilbert spaces, projection and least squares applied to signals and estimation.
- L. N. Trefethen and D. Bau, *Numerical Linear Algebra*, SIAM, 1997 — QR, Gram–Schmidt stability and least-squares conditioning.
