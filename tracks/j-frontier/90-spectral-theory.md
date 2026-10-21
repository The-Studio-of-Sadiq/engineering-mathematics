# Chapter 90 — Spectral Theory: Understanding Operators Through Their Modes

[← Previous](89-pde-beyond-the-basics.md) · [Track J](README.md) · [Table of Contents](../../00-front/table-of-contents.md) · [Next →](91-inverse-problems.md)

> **Level 3 · Track J — Mathematics at the Frontier**
> **Builds on:** [Ch. 10](../../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md), [Ch. 13](../../part-4-differential-equations/13-partial-differential-equations.md), [Ch. 61](../e-abstract-math/61-functional-analysis.md), [Ch. 89](89-pde-beyond-the-basics.md) · **Leads to:** [Ch. 91](91-inverse-problems.md), [Ch. 92](92-mathematical-physics.md), [Ch. 95](95-representation-theory.md) · **Related:** [Ch. 27](../a-algorithms/27-graph-theory.md), [Ch. 34](../b-probability-information/34-random-walks-and-markov-chains.md), [Ch. 38](../b-probability-information/38-monte-carlo.md), [Ch. 99](../e-abstract-math/99-abstract-linear-algebra.md)
> **Reading time:** ~30 min

## 90.1 The hook

In 1910 the physicist Hendrik Lorentz gave a series of lectures in Göttingen. He posed a problem that came out of the theory of black-body radiation: take a cavity of any shape, count how many standing electromagnetic waves it supports below a given frequency, and show that for high frequencies the count depends only on the cavity's *volume*, not its shape. David Hilbert, in the audience, is said to have predicted that he would not live to see it proved. His student Hermann Weyl proved it in 1911.

Fifty-five years later, Mark Kac turned Weyl's theorem into a question that every mathematician now knows by its title: *Can one hear the shape of a drum?* A drum is a membrane clamped along its rim. When struck, it vibrates in a superposition of standing modes, each with its own frequency. An ideal ear hears the complete list of frequencies, the **spectrum**, and nothing else. Weyl showed that this list reveals the drum's area. Later work showed it also reveals the perimeter, and the number of holes. So how much more does it reveal? Is the shape itself encoded in the sound?

This is not just a riddle about drums. The same question, "what does the list of eigenvalues tell me about the object?", sits underneath a remarkable range of things:

- A **graph's** Laplacian eigenvalues tell you whether the network has a bottleneck ([Ch. 27](../a-algorithms/27-graph-theory.md)), and Google-style ranking and spectral clustering are built on them.
- A **Markov chain's** second eigenvalue tells you how fast it forgets its starting point, which is how long a Monte Carlo simulation must run ([Ch. 34](../b-probability-information/34-random-walks-and-markov-chains.md), [Ch. 38](../b-probability-information/38-monte-carlo.md)).
- An **atom's** spectrum, the frequencies of light it emits, was decoded into the energy levels of a Schrödinger operator. Spectroscopy reads the structure of matter from a list of numbers.
- A **bridge's** or **aircraft wing's** lowest natural frequencies decide whether wind or engines can shake it apart, which is the modal analysis of [§10.8](../../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md#108-modal-analysis).

> **How much of an operator, and of the shape it lives on, can you recover from its spectrum alone? And when an operator has no eigenvectors at all, what is its spectrum made of?**

The second half of that question is where matrices stop being a good guide. The first half was answered for drums in 1992, and the answer is a surprise.

## 90.2 The big idea

**An operator is an orchestra of modes.** For a symmetric matrix the spectral theorem says: there is an orthonormal basis of eigenvectors, and in that basis the matrix just multiplies each coordinate by its eigenvalue ([§10.5](../../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md#105-symmetric-matrices)). Spectral theory is the same statement for operators on functions. For the Laplacian on a drum the "eigenvectors" are mode shapes, the **eigenfunctions** $\phi_k(x,y)$, and the eigenvalues are squared frequencies:

$$
-\nabla^2\phi_k=\lambda_k\phi_k\ \text{ inside the drum},\qquad \phi_k=0\ \text{ on the rim}.
$$

Once you know the modes, any question about the operator becomes a question about numbers. Heat flow multiplies mode $k$ by $e^{-\lambda_k t}$. A struck membrane oscillates mode $k$ as $\cos(\sqrt{\lambda_k}\,t)$. Solving $-\nabla^2u=f$ divides mode $k$ of $f$ by $\lambda_k$. Every one of these is "apply a function to the eigenvalues". That is exactly why separation of variables works ([§13.3](../../part-4-differential-equations/13-partial-differential-equations.md#133-separation-of-variables)).

![Eight square heat maps showing vibration modes of a square drum, coloured red and blue with black nodal lines: the (1,1) fundamental, the (2,1) and (1,2) modes with one straight nodal line, the (2,2) mode, two combinations of (2,1) and (1,2) whose nodal lines are the diagonals, a (3,1)+(1,3) combination with a closed nodal curve, and the (3,2) mode](../../images/ch90-drum-modes.svg)

*Modes of the unit square, $\phi_{mn}=\sin(m\pi x)\sin(n\pi y)$ with $\lambda=\pi^2(m^2+n^2)$. Because the square is symmetric, $(2,1)$ and $(1,2)$ share an eigenvalue, and every combination of them is also a mode: the bottom row shows diagonal nodal lines and even a closed nodal curve, none of which appears in any single product of sines.*

**The low notes know the big picture; the high notes know the details.** The lowest eigenvalue $\lambda_1$ is a single number that summarizes the drum's overall size and shape: how hard it is to bend. High eigenvalues come in a dense crowd, and their *density*, how many there are below $\lambda$, is controlled by the area. Between these extremes lies the subtle information: the perimeter, the corners, the holes, and possibly the shape.

**When there are no modes.** On a bounded drum the modes are a discrete list. On an infinite string there are no standing waves at all: a wave $e^{j\omega x}$ has the right shape but infinite energy, so it is not an eigenfunction. The spectrum does not disappear; it smears into a continuum. Spectral theory's deepest result says that even then, the operator is still "multiplication by a number" in a suitable sense, and the Fourier transform of [Ch. 16](../../part-5-transforms/16-fourier-transforms.md) is the prototype.

## 90.3 The mathematics

Throughout, $\Omega$ is a bounded region in the plane, and $-\nabla^2$ acts on functions that vanish on its boundary $\partial\Omega$ (a **Dirichlet** drum). The inner product is $\langle u,v\rangle=\iint_\Omega u\,v\,dA$.

### Self-adjointness: why the modes are real and orthogonal

> [!theorem] The Dirichlet Laplacian is symmetric and positive
> For smooth $u,v$ vanishing on $\partial\Omega$: $\langle -\nabla^2u,v\rangle=\langle u,-\nabla^2v\rangle=\iint_\Omega\nabla u\cdot\nabla v\,dA$. Consequently its eigenvalues are real and positive, and eigenfunctions with different eigenvalues are orthogonal.

*Proof.* Green's first identity ([§7.11](../../part-2-calculus/07-vector-calculus.md#711-divergence-theorem)) gives $\iint(-\nabla^2u)v\,dA=\iint\nabla u\cdot\nabla v\,dA-\oint v\,\partial_nu\,ds$, and the boundary term vanishes because $v=0$ there. The result is symmetric in $u$ and $v$. If $-\nabla^2\phi=\lambda\phi$ with $\phi\ne0$, take $u=v=\phi$: $\lambda\lVert\phi\rVert^2=\iint\lvert\nabla\phi\rvert^2>0$, since a function with zero gradient that vanishes on the boundary is zero. If $-\nabla^2\phi=\lambda\phi$ and $-\nabla^2\psi=\mu\psi$, then $\lambda\langle\phi,\psi\rangle=\langle-\nabla^2\phi,\psi\rangle=\langle\phi,-\nabla^2\psi\rangle=\mu\langle\phi,\psi\rangle$, so $\lambda\ne\mu$ forces $\langle\phi,\psi\rangle=0$. ∎

This is the matrix proof of [§10.5](../../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md#105-symmetric-matrices) with the transpose replaced by integration by parts. What it does *not* give is existence: a matrix always has eigenvalues, an operator need not. For the drum, existence comes from [Ch. 61](../e-abstract-math/61-functional-analysis.md)'s spectral theorem for compact self-adjoint operators, applied to the inverse $(-\nabla^2)^{-1}$, which solves the weak problem of [Ch. 89](89-pde-beyond-the-basics.md). That inverse is compact because finite-energy functions on a bounded domain cannot oscillate arbitrarily fast without paying for it (Rellich's theorem). The conclusion:

> [!theorem] Eigenfunction expansion (stated)
> On a bounded domain, the Dirichlet Laplacian has eigenvalues $0<\lambda_1\le\lambda_2\le\cdots\to\infty$, each of finite multiplicity, with an orthonormal basis of eigenfunctions $\phi_k$ of $L^2(\Omega)$. Every $f\in L^2(\Omega)$ satisfies $f=\sum_k\langle f,\phi_k\rangle\phi_k$, with convergence in the $L^2$ norm.

**Applied to the heat equation.** If $u_t=\nabla^2u$ with $u=0$ on the boundary and $u(\cdot,0)=f$, then

$$
u(\cdot,t)=\sum_k e^{-\lambda_kt}\langle f,\phi_k\rangle\,\phi_k ,\qquad \lVert u(\cdot,t)\rVert\le e^{-\lambda_1t}\lVert f\rVert .
$$

The inequality follows from Parseval: every coefficient shrinks by at least $e^{-\lambda_1t}$. So $\lambda_1$ is the cooling rate, and in the long run $u\approx e^{-\lambda_1t}\langle f,\phi_1\rangle\phi_1$: everything except the fundamental mode is forgotten. For an insulated (Neumann) body, $\lambda_1=0$ is the constant mode, the final uniform temperature, and the approach to it is governed by the **spectral gap** $\lambda_2$. Replace "temperature" by "probability" and this is the convergence of a Markov chain to equilibrium, which is why mixing-time bounds in [Ch. 38](../b-probability-information/38-monte-carlo.md) are spectral-gap bounds.

### The variational principle

> [!theorem] Rayleigh quotient and min–max
> Let $R(v)=\dfrac{\iint\lvert\nabla v\rvert^2\,dA}{\iint v^2\,dA}$ for nonzero $v$ in $H^1_0(\Omega)$. Then
> $$
> \lambda_1=\min_{v\ne0}R(v),\qquad \lambda_k=\min_{\dim S=k}\ \max_{v\in S,\,v\ne0}R(v),
> $$
> where $S$ ranges over $k$-dimensional subspaces of $H^1_0(\Omega)$ (Courant–Fischer).

*Proof of the first statement.* Expand $v=\sum_kc_k\phi_k$. By orthogonality, $\iint v^2=\sum c_k^2$. Because $\iint\nabla\phi_i\cdot\nabla\phi_k=\langle-\nabla^2\phi_i,\phi_k\rangle=\lambda_i\delta_{ik}$, also $\iint\lvert\nabla v\rvert^2=\sum\lambda_kc_k^2$. Hence

$$
R(v)=\frac{\sum_k\lambda_kc_k^2}{\sum_kc_k^2}\ \ge\ \frac{\sum_k\lambda_1c_k^2}{\sum_kc_k^2}=\lambda_1 ,
$$

with equality for $v=\phi_1$. ∎ The min–max formula is proved the same way: any $k$-dimensional $S$ contains a nonzero $v$ orthogonal to $\phi_1,\dots,\phi_{k-1}$ (a dimension count), and for that $v$ the same estimate gives $R(v)\ge\lambda_k$; the subspace spanned by $\phi_1,\dots,\phi_k$ attains it.

This one principle settles several questions at once.

- **Poincaré's constant.** On $(0,1)$, $\lambda_1=\pi^2$, so $\int v^2\le\pi^{-2}\int(v')^2$ for every $v\in H^1_0$: the sharp constant $1/\pi$ that [Ch. 89](89-pde-beyond-the-basics.md) promised.
- **Cheap estimates.** Any trial shape gives an upper bound. For $v=x(1-x)$ on $(0,1)$, $R(v)=\frac{1/3}{1/30}=10$, within $1.3\%$ of $\pi^2=9.8696$, with no differential equation solved.
- **Galerkin eigenvalues are always too high.** A finite-element method restricts the subspaces $S$ in the min–max formula to subspaces of $V_h$. Minimizing over fewer candidates cannot give a smaller minimum, so $\lambda_k^h\ge\lambda_k$ for every $k$ and every mesh. That is the clue left open in Ch. 89, where ten elements gave $9.95>\pi^2$.
- **Domain monotonicity.** If $\Omega_1\subset\Omega_2$, every trial function for $\Omega_1$ (extended by zero) is a trial function for $\Omega_2$, so $\lambda_k(\Omega_2)\le\lambda_k(\Omega_1)$. Bigger drums sound lower.
- **Bottlenecks.** For a dumbbell, two rooms joined by a thin corridor, a Neumann trial function that is $+1$ in one room and $-1$ in the other has tiny gradient energy, so $\lambda_2$ is tiny. Cheeger's inequality (1970) makes this two-sided: $\lambda_2\ge h^2/4$, where $h$ is the smallest ratio of "cut length" to "area cut off". A small gap *forces* a bottleneck. The graph version is the clue left open in [Ch. 27](../a-algorithms/27-graph-theory.md).

### Weyl's law: hearing the area

Let $N(\lambda)$ be the number of eigenvalues $\le\lambda$, counted with multiplicity.

> [!theorem] Weyl's law
> For a bounded plane domain of area $A$, $\displaystyle N(\lambda)\sim\frac{A}{4\pi}\,\lambda$ as $\lambda\to\infty$. In $d$ dimensions, $N(\lambda)\sim\dfrac{\omega_d\,\lvert\Omega\rvert}{(2\pi)^d}\lambda^{d/2}$, where $\omega_d$ is the volume of the unit ball.

*Proof for a rectangle* $[0,a]\times[0,b]$. The eigenfunctions are $\sin(m\pi x/a)\sin(n\pi y/b)$ with $\lambda=\pi^2(m^2/a^2+n^2/b^2)$, $m,n\ge1$. So $N(\lambda)$ counts lattice points $(m,n)$ with $m,n\ge1$ inside the quarter-ellipse $E=\{m^2/a^2+n^2/b^2\le\lambda/\pi^2,\ m,n\ge0\}$, whose area is $\tfrac{\pi}{4}\cdot\tfrac{a\sqrt\lambda}{\pi}\cdot\tfrac{b\sqrt\lambda}{\pi}=\tfrac{ab}{4\pi}\lambda$. Attach to each counted point the unit square $[m-1,m]\times[n-1,n]$ that has it as its top-right corner. These squares do not overlap, and each lies inside $E$ (moving down or left only decreases $m^2/a^2+n^2/b^2$). So $N(\lambda)\le\operatorname{area}(E)=\tfrac{ab}{4\pi}\lambda$. Conversely, a point of $E$ with positive coordinates lies in the square whose top-right corner is the lattice point just above and to the right of it. If the point is more than $\sqrt2$ from the curved edge, that corner is in $E$ and is counted. So the squares cover all of $E$ except a strip of width $\sqrt2$ along the curved edge, whose area is $O(\sqrt\lambda)$. Hence $N(\lambda)=\tfrac{ab}{4\pi}\lambda+O(\sqrt\lambda)$, and the leading term is $A\lambda/(4\pi)$. ∎

The general case glues rectangles together, using domain monotonicity to squeeze the count from above and below. The $O(\sqrt\lambda)$ correction is not random. If lattice points on the two axes were allowed, the count would overshoot the area by about half of the points lying on those axes, and there are about $\tfrac{(a+b)\sqrt\lambda}{\pi}$ of them. Excluding them ($m,n\ge1$) leaves a net deficit of about $\tfrac{(a+b)\sqrt\lambda}{2\pi}$, which is the second term of Weyl's conjecture,

$$
N(\lambda)\approx\frac{A}{4\pi}\lambda-\frac{P}{4\pi}\sqrt{\lambda},
$$

where $P$ is the perimeter. For the unit square at $\lambda=10{,}000$ the true count is $764$; one term predicts $795.8$, two terms predict $763.9$.

![Left: the staircase counting function of the unit square up to lambda = 1500, hugging the two-term Weyl curve, with a discrete 10 by 10 grid's staircase drifting above and saturating at 100. Right: scaled counting functions for 10, 20 and 40 grids collapsing onto one curve that leaves the straight Weyl line and reaches 1 at lambda h squared = 8](../../images/ch90-weyl-law.svg)

*Left: the eigenvalue count of a square drum follows Weyl's law, and the two-term version almost exactly. A $10\times10$ finite-difference grid agrees at first, then overcounts and stops at $100$ modes. Right: every grid, rescaled, follows the same curve, which bends away from Weyl long before its last mode at $\lambda=8/h^2$. A discretization only hears the low notes correctly.*

A cleaner way to hear the geometry is the **heat trace** $Z(t)=\sum_ke^{-\lambda_kt}$, which is determined by the spectrum. For a domain with a smooth boundary and $g$ holes, as $t\to0^+$,

$$
Z(t)=\frac{A}{4\pi t}-\frac{P}{8\sqrt{\pi t}}+\frac{1-g}{6}+O(\sqrt t).
$$

Area, perimeter and the number of holes can all be heard. (Corners add their own constants: for the square the constant is $\tfrac14$, which you can check from $Z=\big(\sum_me^{-\pi^2m^2t}\big)^2$.)

**And the shape?** In 1992 Carolyn Gordon, David Webb and Scott Wolpert built two different polygons, each assembled from seven congruent half-squares, with exactly the same spectrum. One cannot hear the shape of a drum. The proof is a transplantation: cut an eigenfunction on one drum into its seven pieces and reassemble them, with sign changes, on the other.

> [!mystery] Can one hear the shape of a *smooth* drum?
> The Gordon–Webb–Wolpert drums have corners. Nobody knows whether two different smooth convex drums can share a spectrum. Zelditch (2009) proved that among analytic domains with one mirror symmetry, the spectrum does determine the shape, and Hezari and Zelditch later extended the result to nearly circular ellipses. The general question, which Kac really asked, is an open problem. The tools that might settle it come from billiards: the heat trace and its cousin, the wave trace, carry the lengths of periodic billiard paths inside the drum.

> [!mystery] The second term of Weyl's law
> The two-term formula above is called Weyl's *conjecture*. Ivrii proved it in 1980, but only under a condition on billiards: the periodic trajectories bouncing around inside the drum must form a set of measure zero. Whether every domain satisfies this condition is unknown. It is open even for convex domains with smooth boundary in the plane. A question about counting musical notes is waiting on a question about pool tables.

### When there are no eigenvectors: continuous spectrum

On the whole line, $-d^2/dx^2$ acting on $L^2(\mathbb{R})$ has no eigenfunctions: the solutions of $-u''=\lambda u$ are $e^{\pm j\sqrt\lambda x}$, or exponentials, none of them square-integrable. But the Fourier transform ([§16.3](../../part-5-transforms/16-fourier-transforms.md#163-transform-properties)) turns differentiation into multiplication, $\mathcal{F}\{-u''\}=\omega^2\,\mathcal{F}\{u\}$. In the frequency domain the operator is just "multiply by $\omega^2$": diagonal, with a continuum of diagonal entries $[0,\infty)$ instead of a list. This is the general picture, and it answers the clue of [Ch. 61](../e-abstract-math/61-functional-analysis.md).

> [!theorem] Spectral theorem for self-adjoint operators, multiplication form (stated)
> For every self-adjoint operator $T$ on a Hilbert space $\mathcal{H}$, possibly unbounded, there is a measure space $(M,\mu)$, a real function $g$ on $M$ and a unitary map $U:\mathcal{H}\to L^2(M,\mu)$ with $UTU^{-1}=$ multiplication by $g$. Equivalently, $T=\int\lambda\,dE(\lambda)$ for a projection-valued measure $E$.

Eigenvalues are the atoms of the measure, where a single value of $g$ carries positive weight; continuous spectrum is the rest. The physicist's "eigenstates" $\lvert x_0\rangle$ and plane waves $e^{j\omega x}$ are not vectors in $\mathcal{H}$; they are the kernel of $U$, a coordinate system for the multiplication picture. Real operators often mix both kinds:

![Three spectra drawn on number lines: the string on [0,1] with isolated eigenvalues running off to infinity, the free particle on the whole line with a continuous band from 0 to infinity, and the hydrogen atom with eigenvalues accumulating at 0 from below followed by a continuous band](../../images/ch90-spectrum-types.svg)

*Point, continuous and mixed spectra. Hydrogen's bound states $-1/(2n^2)$ (in atomic units) are the spectral lines of the Balmer and Lyman series; the continuum above $0$ is an electron that has escaped.*

**Beyond self-adjoint.** For non-symmetric operators, such as convection-diffusion or the linearized equations of fluid flow, eigenvalues can be dangerously misleading. A flow whose eigenvalues all predict decay can still amplify a disturbance by a factor of thousands before decaying, which is one route to turbulence. The robust object is the **$\varepsilon$-pseudospectrum**, the set of $z$ where $\lVert(T-zI)^{-1}\rVert>1/\varepsilon$, introduced in [Ch. 99](../e-abstract-math/99-abstract-linear-algebra.md)'s loose threads. For a normal operator it is just the $\varepsilon$-neighbourhood of the spectrum. For a non-normal one it can be enormously larger, and it, not the spectrum, predicts the transient growth.

## 90.4 The algorithm

Discretize a drum on an $n\times n$ grid and the Laplacian becomes a sparse symmetric matrix of size $N=n^2$. For a square it is a **Kronecker sum**, $\mathbf{A}=\mathbf{I}\otimes\mathbf{T}+\mathbf{T}\otimes\mathbf{I}$, where $\mathbf{T}=h^{-2}\operatorname{tridiag}(-1,2,-1)$ is the 1-D second difference. Its eigenvalues are known exactly:

$$
\lambda^h_{mn}=\frac{4}{h^2}\Big(\sin^2\frac{m\pi h}{2}+\sin^2\frac{n\pi h}{2}\Big),\qquad 1\le m,n\le n_{\text{grid}} .
$$

For small $mh$ this is $\pi^2(m^2+n^2)$ minus a correction of order $h^2$, so finite differences *underestimate* the low eigenvalues, the opposite of finite elements. For large $m$ the formula saturates at $8/h^2$: the grid has no notes above that. This is why the right panel of the Weyl figure bends away.

You never want *all* $N$ eigenvalues of a large drum. A dense symmetric eigensolver costs $O(N^3)$ time and $O(N^2)$ memory, so $N=40{,}000$ already needs about 13 GB just to store the matrix. Engineers want the lowest few modes. The workhorse is **inverse iteration with a shift**: factor $\mathbf{A}-\sigma\mathbf{I}$ once, then repeatedly solve with it. Each solve multiplies the component along mode $k$ by $1/(\lambda_k-\sigma)$, so the mode nearest $\sigma$ wins.

```python
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu, eigsh

n = 199; h = 1.0 / (n + 1)                       # interior grid of the unit square
T = sp.diags([-1.0, 2.0, -1.0], [-1, 0, 1], shape=(n, n)) / h**2
I = sp.identity(n)
A = (sp.kron(I, T) + sp.kron(T, I)).tocsc()     # 2-D Laplacian: a Kronecker sum

def inverse_iteration(A, sigma=0.0, iters=50, tol=1e-12):
    """Eigenpair of A nearest the shift sigma."""
    lu = splu(A - sigma * sp.identity(A.shape[0], format="csc"))   # factor once
    v = np.random.default_rng(1).standard_normal(A.shape[0])
    lam = 0.0
    for _ in range(iters):
        w = lu.solve(v)                       # w = (A - sigma I)^{-1} v
        v = w / np.linalg.norm(w)
        lam_new = v @ (A @ v)                 # Rayleigh quotient
        if abs(lam_new - lam) < tol * abs(lam_new): break
        lam = lam_new
    return lam_new, v

lam1, v1 = inverse_iteration(A)
print("inverse iteration:", lam1, " exact 2*pi^2 =", 2 * np.pi**2)
print("Lanczos, 6 smallest:", np.round(eigsh(A, k=6, sigma=0, which="LM")[0], 3))
```

On this $39{,}601$-unknown problem, inverse iteration returns $\lambda_1^h=19.7388$ against the exact $2\pi^2=19.7392$. The last line uses the **Lanczos** method in shift-invert mode, which builds the best approximations from the whole Krylov space $\{\mathbf{v},\mathbf{A}^{-1}\mathbf{v},\mathbf{A}^{-2}\mathbf{v},\dots\}$ instead of keeping only the latest vector. It finds $19.739,\ 49.345,\ 49.345,\ 78.95,\ 98.679,\ 98.679$: the double eigenvalues $5\pi^2$ and $10\pi^2$ appear twice, as they should.

**Cost and failure modes.**

- Inverse iteration converges with ratio $\lvert\lambda_1-\sigma\rvert/\lvert\lambda_2-\sigma\rvert$ per step. A good shift makes it fast; a shift halfway between two eigenvalues makes it stall. Rayleigh-quotient iteration, which updates $\sigma$ each step, converges cubically, at the price of refactoring every step.
- **Degenerate eigenvalues.** For $5\pi^2$ the solver returns *some* orthonormal pair from a two-dimensional eigenspace, and which pair depends on rounding and on the random start. Comparing "mode 2" between two runs, or between two meshes, can fail for this reason alone. Only the eigenspace is well defined.
- **Non-symmetric matrices.** For convection-dominated problems, eigenvalues of the matrix can sit far from those of the operator, and tiny perturbations move them a lot (the pseudospectra of §90.3). Check them, don't trust them.
- **Spectral pollution.** For operators with continuous spectrum, or for Maxwell's equations discretized with the wrong elements, Galerkin methods can produce eigenvalues that converge to nothing at all, sitting in a gap of the true spectrum. The min–max principle protects the bottom of the spectrum, never the interior. Detecting and avoiding pollution is current research.

## 90.5 The experiment

**Hear the area and the perimeter.** For an $a\times b$ rectangle the counting function is a finite sum of floor functions. With $a=b=1$, the blue staircase is the exact $N(\lambda)$ and the red curve is the two-term Weyl law. Then make the drum long and thin with the same area: set $a=4$, $b=0.25$. In the sum, $b$ is the factor in front of the square root and $a$ is the base of `1^{2}`; in the Weyl curves, $ab$ is `1\cdot1` and the perimeter $2(a+b)$ is `2(1+1)`. The one-term line does not move, because the area is unchanged; the staircase drops (at $\lambda=2000$, from 146 modes to 131), and the perimeter term explains by how much.

```desmos-graph
left=0; right=2000; bottom=0; top=170
---
y=\sum_{m=1}^{80}\operatorname{floor}\left(1\cdot\sqrt{\max\left(0,\frac{x}{\pi^{2}}-\frac{m^{2}}{1^{2}}\right)}\right)|blue
y=\frac{1\cdot1}{4\pi}x|grey|dashed
y=\frac{1\cdot1}{4\pi}x-\frac{2(1+1)}{4\pi}\sqrt{x}|red
```

*Exact eigenvalue count of a rectangular drum (blue) against Weyl's law with one term (dashed) and two terms (red). The same area gives the same slope; the perimeter sets the offset.*

**Your own drum.** This script builds the five-point Laplacian on any shape you can describe by a mask, computes all its eigenvalues with a dense solver, and compares the count with Weyl's law. The shapes live in $[-1,1]^2$: a square of side 2, the L-shape of Ch. 89 (area 3), and the unit disk.

```python
import numpy as np

n = 30                                   # knob: grid points per unit length
shape = "L"                              # knob: "square", "L" or "disk"
h = 1.0 / n
g = (np.arange(-n, n) + 0.5) * h         # cell centres covering [-1, 1]
X, Y = np.meshgrid(g, g)
masks = {"square": (np.abs(X) < 1) & (np.abs(Y) < 1),
         "L":      ~((X > 0) & (Y < 0)),
         "disk":   X**2 + Y**2 < 1}
area_perim = {"square": (4, 8), "L": (3, 8), "disk": (np.pi, 2 * np.pi)}
inside = masks[shape]
idx = -np.ones(X.shape, int); idx[inside] = np.arange(inside.sum())
N = inside.sum()

# five-point Dirichlet Laplacian on the cells inside the drum
A = np.zeros((N, N))
for i, j in zip(*np.nonzero(inside)):
    p = idx[i, j]; A[p, p] = 4 / h**2
    for di, dj in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
        a, b = i + di, j + dj
        if 0 <= a < 2 * n and 0 <= b < 2 * n and inside[a, b]:
            A[p, idx[a, b]] = -1 / h**2
lam = np.linalg.eigvalsh(A)              # N x N dense: keep N below ~4000

print(f"{shape}: {N} unknowns; first eigenvalues:", np.round(lam[:6], 2))
Area, P = area_perim[shape]
print("   lambda   count   Weyl 1-term   Weyl 2-term")
for L in [50, 100, 200, 400, 800]:
    w1 = Area * L / (4 * np.pi); w2 = w1 - P * np.sqrt(L) / (4 * np.pi)
    print(f"{L:8d} {np.sum(lam <= L):7d} {w1:12.1f} {w2:12.1f}")
```

For the L-shape with $n=30$ the first eigenvalue comes out as $9.08$. The true value is $9.6397$: this is the famous L-shaped membrane whose first mode is the MATLAB logo, computed to high accuracy by Fox, Henrici and Moler in 1967. The grid is crude, and because the grid's drum is effectively half a cell larger than the true one, its eigenvalues sit below the truth. At $\lambda=400$ the count is $90$, between the one-term ($95.5$) and two-term ($82.8$) predictions. Knobs to turn:

- Raise `n` to 35 (about 3,700 unknowns) and watch $\lambda_1$ creep toward $9.64$. How fast? Compare with the square, whose exact $\lambda_1$ is $\pi^2/2\approx4.935$, and explain the difference using the corner singularity of Ch. 89.
- Look at `lam[1]` and `lam[2]` for the square and the disk. Which are degenerate, and which symmetry is responsible?
- Add your own mask: an annulus, a dumbbell, two disks joined by a thin channel. For the dumbbell, rerun with Neumann conditions (set each diagonal entry to the number of inside neighbours divided by $h^2$) and watch $\lambda_2$ collapse as the channel narrows.

## 90.6 The challenge

1. **Intuition.** Without computing, rank these drums of equal area by their fundamental frequency $\sqrt{\lambda_1}$: a disk, a square, an equilateral triangle, a $4\times\tfrac14$ rectangle. Which principle of §90.3 helps, and which pairs does it fail to decide?
2. **Derivation.** For the equilateral triangle and the square of the same area, compute the two-term Weyl predictions and say which drum has more modes below a large $\lambda$. Then verify the heat-trace constant $\tfrac14$ for the unit square from $Z(t)=\big(\sum_{m\ge1}e^{-\pi^2m^2t}\big)^2$, using $\sum_{m\ge1}e^{-\pi^2m^2t}\approx\tfrac{1}{2\sqrt{\pi t}}-\tfrac12$ for small $t$.
3. **Implementation.** Build the Gordon–Webb–Wolpert isospectral pair from seven half-squares each (look up their coordinates), discretize both with the script of §90.5 on triangular-friendly grids, and compare their first 20 eigenvalues. How close do two genuinely isospectral drums look after discretization, and can you tell them apart from a *non*-isospectral pair with the same area and perimeter?
4. **Open investigation.** Compute the first 2,000 eigenvalues of a rectangle with irrational aspect ratio and of a "stadium" (a rectangle with semicircular ends). Unfold each spectrum by subtracting the two-term Weyl law, form the gaps between consecutive eigenvalues, and plot their histograms. One histogram should look like $e^{-s}$; the other should vanish at $s=0$. Which is which, and why might the shape of a billiard table decide whether eigenvalues repel each other?

## 90.7 Loose threads

Spectral theory turned operators into lists of numbers and the lists into geometry. The bottom of the spectrum measured stiffness, cooling rates and bottlenecks through the min–max principle. The density of the spectrum measured area, and its first correction, the perimeter. And when eigenvectors disappeared, the spectral theorem survived as "multiplication by a function". The challenge's last problem points at the most mysterious clue of all.

> [!mystery] Why do some spectra repel?
> Unfold the eigenvalues of a rectangle and their gaps look like independent random numbers: a Poisson process, with many near-collisions. Do the same for a stadium, whose billiard paths are chaotic, and the gaps follow the statistics of eigenvalues of large random symmetric matrices, the Wigner surmise $\tfrac{\pi s}{2}e^{-\pi s^2/4}$, which vanishes at $s=0$. Bohigas, Giannoni and Schmit conjectured in 1984 that this holds for every chaotic system, and the same statistics appear in heavy nuclei and in the zeros of the Riemann zeta function. Nobody has proved it for a single drum. The physics of it is in [Ch. 92](92-mathematical-physics.md); the proof is an open problem.

Further pointers:

- **Inverse problems.** Singular values are the spectral theory of $\mathbf{K}^\mathsf{T}\mathbf{K}$, and the speed at which they decay decides whether an inverse problem is mildly or catastrophically unstable: [Ch. 91](91-inverse-problems.md).
- **Symmetry and degeneracy.** The double eigenvalue $5\pi^2$ is the square's rotation symmetry in disguise. In general, the multiplicities of a symmetric operator are the dimensions of the irreducible representations of its symmetry group: [Ch. 95](95-representation-theory.md).
- **Graphs.** Everything here has a discrete twin. The graph Laplacian has a min–max principle, a Cheeger inequality and even isospectral pairs, and spectral clustering uses its second eigenvector to cut networks ([Ch. 27](../a-algorithms/27-graph-theory.md), [Ch. 87](../i-dynamical-systems/87-network-science.md)).

## Further reading

- M. Kac, "Can one hear the shape of a drum?", *American Mathematical Monthly* 73 (1966), 1–23. Still the best invitation to the subject.
- C. Gordon, D. Webb and S. Wolpert, "One cannot hear the shape of a drum", *Bulletin of the AMS* 27 (1992), 134–138. Four pages, with pictures of the two drums.
- D. S. Grebenkov and B.-T. Nguyen, "Geometrical structure of Laplacian eigenfunctions", *SIAM Review* 55 (2013), 601–667. A survey of what is known about modes, nodal lines and localization, with many figures.
- E. B. Davies, *Spectral Theory and Differential Operators*, Cambridge University Press, 1995. A rigorous course on the spectral theorem and its use for PDEs.
- L. N. Trefethen and M. Embree, *Spectra and Pseudospectra*, Princeton University Press, 2005. What eigenvalues do not tell you about non-normal matrices and operators.
- V. Ivrii, "100 years of Weyl's law", *Bulletin of Mathematical Sciences* 6 (2016), 379–452. The history and current state of eigenvalue asymptotics, by the person who proved the two-term formula.
