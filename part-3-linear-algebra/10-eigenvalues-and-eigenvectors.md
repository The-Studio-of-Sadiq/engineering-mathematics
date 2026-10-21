# Chapter 10 — Eigenvalues and Eigenvectors

[← Previous](09-vector-spaces-and-linear-transformations.md) · [Table of Contents](../00-front/table-of-contents.md) · [Next →](../part-4-differential-equations/11-ordinary-differential-equations.md)

> **Prerequisites:** [Ch. 2](../part-1-foundations/02-complex-numbers.md) (complex numbers, polar form), [Ch. 8](08-matrices-and-linear-systems.md) (determinants, inverses, Gaussian elimination, rank), [Ch. 9](09-vector-spaces-and-linear-transformations.md) (bases, change of basis, null space, inner products, orthogonality, Gram–Schmidt)
> **Learning outcomes:** After this chapter you can
> - form and solve the characteristic equation of a square matrix, and use trace, determinant and the Cayley–Hamilton theorem as checks and shortcuts;
> - compute eigenvalues, eigenvectors and eigenspaces (real and complex), and decide whether a matrix is diagonalizable;
> - diagonalize a matrix and use the factorization to compute powers $\mathbf{A}^k$, the matrix exponential $e^{\mathbf{A}t}$ and the long-run behaviour of difference equations;
> - apply the spectral theorem to symmetric and Hermitian matrices, find principal axes of quadratic forms, and test a matrix for positive definiteness;
> - compute and interpret the singular value decomposition, including matrix norm, condition number, rank, low-rank approximation and the pseudoinverse;
> - decouple vibrating and oscillating systems into normal modes, and judge the stability of linear (and linearized) continuous- and discrete-time systems from their eigenvalues.

Chapters 8 and 9 asked what a matrix *does* to a vector: it maps $\mathbf{x}$ to $\mathbf{A}\mathbf{x}$, generally rotating it, stretching it and mixing its components. This chapter asks a sharper question: **are there special directions that the matrix does not rotate at all?** Along such a direction the matrix acts like a plain number, $\mathbf{A}\mathbf{v}=\lambda\mathbf{v}$. The number $\lambda$ is an *eigenvalue* (German *eigen*, "own" or "characteristic") and $\mathbf{v}$ an *eigenvector*.

The idea sounds modest, but it is one of the most productive in applied mathematics. If a matrix has enough eigenvectors to form a basis, then in that basis the matrix is *diagonal*: a coupled problem in $n$ unknowns falls apart into $n$ independent scalar problems. A system of coupled differential equations $\dot{\mathbf{x}}=\mathbf{A}\mathbf{x}$ becomes $n$ separate first-order equations $\dot y_i=\lambda_iy_i$, each solved by an exponential $e^{\lambda_it}$. Coupled masses or coupled resonant circuits separate into *normal modes*, each vibrating at its own frequency. Repeated multiplication $\mathbf{A}^k$ — the heart of every discrete-time recursion, Markov chain and iterative algorithm — becomes repeated multiplication of numbers, $\lambda_i^k$. And the question "is this system stable?" reduces to "do all the eigenvalues lie in the left half-plane (or inside the unit circle)?"

For an electrical engineer, eigenvalues are hiding in plain sight. The **natural frequencies** of a circuit, the **poles** of its transfer function, the **resonant modes** of a cavity or a coupled-line filter, the **principal axes** of a polarisation ellipse or of an anisotropic dielectric tensor, the **principal components** of a set of sensor measurements, the **singular values** that measure how close a MIMO channel is to losing a spatial stream — all are eigenvalue problems. Chapters 11–14 use the results of this chapter constantly.

The plan is as follows. Sections 10.1–10.4 develop the general theory: the characteristic equation, eigenvectors, diagonalization and matrix powers (including the matrix exponential). Sections 10.5–10.7 specialise to the symmetric and Hermitian matrices that arise from energy, power and covariance, culminating in the singular value decomposition, which works for *every* matrix. Sections 10.8–10.9 turn the theory onto dynamics: normal-mode analysis and stability.

## 10.1 Characteristic equations

### Motivation

Suppose we want a nonzero vector $\mathbf{v}$ and a scalar $\lambda$ with $\mathbf{A}\mathbf{v}=\lambda\mathbf{v}$. Unlike $\mathbf{A}\mathbf{x}=\mathbf{b}$, this is not a linear system for $\mathbf{v}$ in the usual sense, because $\lambda$ is unknown too and multiplies $\mathbf{v}$. The trick is to separate the two unknowns. We first ask *for which* $\lambda$ the equation has a nonzero solution at all; that condition turns out to be a single polynomial equation in $\lambda$ — the **characteristic equation**. Only after the eigenvalues are known do we go back and find the eigenvectors, which is then an ordinary homogeneous linear system.

The characteristic equation is not new to you. When you solved $a\ddot y+b\dot y+cy=0$ by trying $y=e^{st}$ and obtained $as^2+bs+c=0$, that quadratic was the characteristic equation of a $2\times2$ matrix in disguise (Example 10.3 below makes the connection exact). The roots of the denominator of a transfer function — the poles — are likewise eigenvalues of the system's state matrix (Section 10.9).

### Derivation

Rewrite $\mathbf{A}\mathbf{v}=\lambda\mathbf{v}$ as

$$
\mathbf{A}\mathbf{v}-\lambda\mathbf{I}\mathbf{v}=\mathbf{0}
\quad\Longleftrightarrow\quad
(\mathbf{A}-\lambda\mathbf{I})\,\mathbf{v}=\mathbf{0}.
$$

The identity $\mathbf{I}$ is needed: $\mathbf{A}-\lambda$ (matrix minus number) is meaningless. For a fixed number $\lambda$ this *is* a homogeneous linear system. It always has the trivial solution $\mathbf{v}=\mathbf{0}$, and by the invertible matrix theorem of [Ch. 8](08-matrices-and-linear-systems.md) it has a nonzero solution if and only if the coefficient matrix $\mathbf{A}-\lambda\mathbf{I}$ is singular, i.e. if and only if its determinant vanishes.

> [!definition] Characteristic polynomial and characteristic equation
> Let $\mathbf{A}$ be an $n\times n$ matrix (real or complex). Its **characteristic polynomial** is
> $$
> p(\lambda)=\det(\mathbf{A}-\lambda\mathbf{I}),
> $$
> a polynomial of degree exactly $n$ in $\lambda$. The equation $p(\lambda)=0$ is the **characteristic equation**, and its roots are the **eigenvalues** of $\mathbf{A}$. The multiset of eigenvalues is the **spectrum** of $\mathbf{A}$, and $\rho(\mathbf{A})=\max_i\lvert\lambda_i\rvert$ is the **spectral radius**.

Some texts use $\det(\lambda\mathbf{I}-\mathbf{A})$ instead; the two differ by the factor $(-1)^n$ and have the same roots. The form $\det(\lambda\mathbf{I}-\mathbf{A})$ has leading coefficient $+1$ (it is *monic*), which is convenient for control-system work where it matches the denominator of a transfer function. We shall use whichever is tidier and say which.

Why is the degree exactly $n$? Expanding the determinant as a sum of signed products, one per permutation ([Ch. 8](08-matrices-and-linear-systems.md)), the only product that contains $n$ factors of $\lambda$ is the product of the diagonal entries $(a_{11}-\lambda)(a_{22}-\lambda)\cdots(a_{nn}-\lambda)$, which contributes $(-\lambda)^n$. Every other permutation leaves out at least two diagonal entries and therefore has degree at most $n-2$.

### Coefficients: trace, determinant and principal minors

The same expansion gives two coefficients immediately. The coefficient of $\lambda^{n-1}$ also comes only from the diagonal product, and equals $(-1)^{n-1}(a_{11}+\dots+a_{nn})$. The constant term is $p(0)=\det\mathbf{A}$. On the other hand, by the fundamental theorem of algebra $p$ factors over $\mathbb{C}$ as

$$
p(\lambda)=(\lambda_1-\lambda)(\lambda_2-\lambda)\cdots(\lambda_n-\lambda),
$$

whose $\lambda^{n-1}$ coefficient is $(-1)^{n-1}(\lambda_1+\dots+\lambda_n)$ and whose constant term is $\lambda_1\lambda_2\cdots\lambda_n$. Comparing:

> [!theorem] Trace and determinant
> If $\lambda_1,\dots,\lambda_n$ are the eigenvalues of $\mathbf{A}$, repeated according to algebraic multiplicity, then
> $$
> \lambda_1+\lambda_2+\dots+\lambda_n=\operatorname{tr}\mathbf{A}=\sum_{i=1}^n a_{ii},
> \qquad
> \lambda_1\lambda_2\cdots\lambda_n=\det\mathbf{A}.
> $$
> In particular $\mathbf{A}$ is singular if and only if $0$ is an eigenvalue.

More generally, the coefficient of $\lambda^{n-k}$ in $\det(\lambda\mathbf{I}-\mathbf{A})$ is $(-1)^kE_k$, where $E_k$ is the sum of all $k\times k$ **principal minors** of $\mathbf{A}$ (determinants of submatrices obtained by keeping the same $k$ rows and columns). This gives the very useful closed forms

$$
\begin{aligned}
n=2:&\quad \det(\lambda\mathbf{I}-\mathbf{A})=\lambda^2-(\operatorname{tr}\mathbf{A})\,\lambda+\det\mathbf{A},\\[2pt]
n=3:&\quad \det(\lambda\mathbf{I}-\mathbf{A})=\lambda^3-(\operatorname{tr}\mathbf{A})\,\lambda^2+E_2\,\lambda-\det\mathbf{A},
\end{aligned}
$$

with $E_2=(a_{11}a_{22}-a_{12}a_{21})+(a_{11}a_{33}-a_{13}a_{31})+(a_{22}a_{33}-a_{23}a_{32})$.

For a $2\times2$ matrix the quadratic formula then gives the eigenvalues in one line:

$$
\lambda_{1,2}=\frac{\operatorname{tr}\mathbf{A}}{2}\pm\sqrt{\left(\frac{\operatorname{tr}\mathbf{A}}{2}\right)^2-\det\mathbf{A}} .
$$

The sign of the discriminant $(\operatorname{tr}\mathbf{A})^2-4\det\mathbf{A}$ decides whether the eigenvalues are real and distinct, real and repeated, or a complex-conjugate pair — the same three cases (overdamped, critically damped, underdamped) you meet for second-order circuits in [Ch. 12](../part-4-differential-equations/12-linear-differential-equations.md). Section 10.9 turns this observation into a complete classification of two-dimensional dynamics.

### Algebraic multiplicity and complex roots

> [!definition] Algebraic multiplicity
> If $p(\lambda)=(\lambda_0-\lambda)^m q(\lambda)$ with $q(\lambda_0)\ne0$, the eigenvalue $\lambda_0$ has **algebraic multiplicity** $m$. The multiplicities of the distinct eigenvalues add up to $n$.

Two facts follow from elementary polynomial algebra ([Ch. 1](../part-1-foundations/01-mathematical-language-and-algebra.md), [Ch. 2](../part-1-foundations/02-complex-numbers.md)):

- **Every $n\times n$ matrix has exactly $n$ eigenvalues in $\mathbb{C}$**, counted with multiplicity. A real matrix may have *no* real eigenvalues (a rotation of the plane is the standard example), so eigenvalue theory naturally lives over the complex numbers.
- **A real matrix has real characteristic-polynomial coefficients, so its non-real eigenvalues occur in conjugate pairs** $\lambda,\lambda^*$. A real $3\times3$ matrix therefore always has at least one real eigenvalue.

### Special structures

Several structures let you read eigenvalues off without solving anything:

- **Triangular and diagonal matrices.** If $\mathbf{A}$ is upper or lower triangular, so is $\mathbf{A}-\lambda\mathbf{I}$, and its determinant is the product of its diagonal entries. Hence *the eigenvalues of a triangular matrix are its diagonal entries.*
- **Block-triangular matrices.** If $\mathbf{A}=\begin{bmatrix}\mathbf{A}_{11}&\mathbf{A}_{12}\\\mathbf{0}&\mathbf{A}_{22}\end{bmatrix}$ with square diagonal blocks, then $\det(\mathbf{A}-\lambda\mathbf{I})=\det(\mathbf{A}_{11}-\lambda\mathbf{I})\det(\mathbf{A}_{22}-\lambda\mathbf{I})$: the spectrum is the union of the spectra of the diagonal blocks. This is how a system made of a cascade of subsystems inherits the natural frequencies of each stage.
- **Transpose.** $\det(\mathbf{A}^\mathsf{T}-\lambda\mathbf{I})=\det\big((\mathbf{A}-\lambda\mathbf{I})^\mathsf{T}\big)=\det(\mathbf{A}-\lambda\mathbf{I})$, so $\mathbf{A}$ and $\mathbf{A}^\mathsf{T}$ have the same eigenvalues (but usually different eigenvectors).
- **Similar matrices.** If $\mathbf{B}=\mathbf{S}^{-1}\mathbf{A}\mathbf{S}$ for some invertible $\mathbf{S}$, then
$$
\det(\mathbf{B}-\lambda\mathbf{I})=\det\big(\mathbf{S}^{-1}(\mathbf{A}-\lambda\mathbf{I})\mathbf{S}\big)=\det(\mathbf{S}^{-1})\det(\mathbf{A}-\lambda\mathbf{I})\det(\mathbf{S})=\det(\mathbf{A}-\lambda\mathbf{I}).
$$
Similar matrices have the same characteristic polynomial, hence the same eigenvalues, trace and determinant. Since $\mathbf{S}$ is a change of basis ([Ch. 9](09-vector-spaces-and-linear-transformations.md)), this says that **eigenvalues belong to the linear transformation, not to the particular matrix that represents it.** Choosing different state variables for the same circuit changes the state matrix but never its eigenvalues — the natural frequencies are physical.

### Companion matrices

Every monic polynomial is the characteristic polynomial of some matrix. For $q(\lambda)=\lambda^n+c_{n-1}\lambda^{n-1}+\dots+c_1\lambda+c_0$, the **companion matrix**

$$
\mathbf{C}=\begin{bmatrix}
0&1&0&\cdots&0\\
0&0&1&\cdots&0\\
\vdots&&&\ddots&\vdots\\
0&0&0&\cdots&1\\
-c_0&-c_1&-c_2&\cdots&-c_{n-1}
\end{bmatrix}
$$

satisfies $\det(\lambda\mathbf{I}-\mathbf{C})=q(\lambda)$ (expand along the last row, or verify for $n=2$ and induct). The companion matrix is exactly what you obtain when you convert the scalar ODE $y^{(n)}+c_{n-1}y^{(n-1)}+\dots+c_0y=0$ to first-order form with the state $\mathbf{x}=(y,\dot y,\dots,y^{(n-1)})$. It is also how numerical software finds polynomial roots: `numpy.roots` builds the companion matrix and calls an eigenvalue routine, because eigenvalue algorithms are more robust than direct root-finding.

### The Cayley–Hamilton theorem

> [!theorem] Cayley–Hamilton
> Every square matrix satisfies its own characteristic equation: if $p(\lambda)=\det(\lambda\mathbf{I}-\mathbf{A})=\lambda^n+c_{n-1}\lambda^{n-1}+\dots+c_0$, then
> $$
> p(\mathbf{A})=\mathbf{A}^n+c_{n-1}\mathbf{A}^{n-1}+\dots+c_1\mathbf{A}+c_0\mathbf{I}=\mathbf{0}.
> $$

*Proof sketch.* The adjugate identity of [Ch. 8](08-matrices-and-linear-systems.md) gives $(\lambda\mathbf{I}-\mathbf{A})\operatorname{adj}(\lambda\mathbf{I}-\mathbf{A})=p(\lambda)\mathbf{I}$. The entries of the adjugate are polynomials of degree at most $n-1$ in $\lambda$, so we can write $\operatorname{adj}(\lambda\mathbf{I}-\mathbf{A})=\mathbf{B}_{n-1}\lambda^{n-1}+\dots+\mathbf{B}_1\lambda+\mathbf{B}_0$ with constant matrices $\mathbf{B}_k$. Multiplying out and equating coefficients of each power of $\lambda$ gives $n+1$ matrix equations: $\mathbf{B}_{n-1}=\mathbf{I}$, $\mathbf{B}_{k-1}-\mathbf{A}\mathbf{B}_k=c_k\mathbf{I}$, and $-\mathbf{A}\mathbf{B}_0=c_0\mathbf{I}$. Multiply the equation for $\lambda^k$ on the left by $\mathbf{A}^k$ and add them all: the left side telescopes to $\mathbf{0}$ and the right side is $p(\mathbf{A})$. (Simply "substituting $\lambda=\mathbf{A}$ into $\det(\lambda\mathbf{I}-\mathbf{A})$" is *not* a valid proof, because the determinant of a matrix with matrix entries is not defined that way.) $\square$

Cayley–Hamilton has three practical uses. First, $\mathbf{A}^n$ — and therefore every higher power — is a linear combination of $\mathbf{I},\mathbf{A},\dots,\mathbf{A}^{n-1}$, so any polynomial or power series in $\mathbf{A}$ (including $e^{\mathbf{A}t}$, Section 10.4) can be reduced to a polynomial of degree at most $n-1$. Second, if $c_0=\pm\det\mathbf{A}\ne0$ we can solve for the inverse:

$$
\mathbf{A}^{-1}=-\frac{1}{c_0}\left(\mathbf{A}^{n-1}+c_{n-1}\mathbf{A}^{n-2}+\dots+c_1\mathbf{I}\right).
$$

Third, it underlies the controllability and observability tests of state-space control theory, where powers of $\mathbf{A}$ beyond $n-1$ never add new directions.

### Worked examples

> [!example] Worked Example 10.1 — A $2\times2$ characteristic equation
> Find the eigenvalues of $\mathbf{A}=\begin{bmatrix}4&1\\2&3\end{bmatrix}$.
>
> **Step 1 — form $\mathbf{A}-\lambda\mathbf{I}$.** Subtract $\lambda$ from the diagonal only:
> $$
> \mathbf{A}-\lambda\mathbf{I}=\begin{bmatrix}4-\lambda&1\\2&3-\lambda\end{bmatrix}.
> $$
> **Step 2 — take the determinant.**
> $$
> p(\lambda)=(4-\lambda)(3-\lambda)-(1)(2)=\lambda^2-7\lambda+12-2=\lambda^2-7\lambda+10 .
> $$
> **Step 3 — check with the shortcut.** $\operatorname{tr}\mathbf{A}=7$ and $\det\mathbf{A}=12-2=10$, so $\lambda^2-7\lambda+10$ is confirmed.
>
> **Step 4 — solve.** $\lambda^2-7\lambda+10=(\lambda-5)(\lambda-2)$, so $\lambda_1=5$, $\lambda_2=2$.
>
> **Check:** $\lambda_1+\lambda_2=7=\operatorname{tr}\mathbf{A}$ and $\lambda_1\lambda_2=10=\det\mathbf{A}$. $\checkmark$

> [!example] Worked Example 10.2 — A $3\times3$ tridiagonal matrix
> Find the eigenvalues of
> $$
> \mathbf{B}=\begin{bmatrix}2&1&0\\1&2&1\\0&1&2\end{bmatrix}.
> $$
> **Step 1 — expand along the first row.**
> $$
> \det(\mathbf{B}-\lambda\mathbf{I})=(2-\lambda)\begin{vmatrix}2-\lambda&1\\1&2-\lambda\end{vmatrix}-1\cdot\begin{vmatrix}1&1\\0&2-\lambda\end{vmatrix}
> =(2-\lambda)\big[(2-\lambda)^2-1\big]-(2-\lambda).
> $$
> **Step 2 — factor before expanding.** The common factor $(2-\lambda)$ is the key labour-saver:
> $$
> \det(\mathbf{B}-\lambda\mathbf{I})=(2-\lambda)\big[(2-\lambda)^2-2\big].
> $$
> **Step 3 — solve.** Either $\lambda=2$, or $(2-\lambda)^2=2$, i.e. $\lambda=2\pm\sqrt2$. The eigenvalues are
> $$
> \lambda_1=2+\sqrt2\approx3.4142,\qquad \lambda_2=2,\qquad\lambda_3=2-\sqrt2\approx0.5858 .
> $$
> **Step 4 — check with principal minors.** $\operatorname{tr}\mathbf{B}=6$; $E_2=(4-1)+(4-0)+(4-1)=10$; $\det\mathbf{B}=2(4-1)-1(2-0)=4$. So $\det(\lambda\mathbf{I}-\mathbf{B})=\lambda^3-6\lambda^2+10\lambda-4$. Sum of roots $=6$ $\checkmark$; product $=2\,(4-2)=4$ $\checkmark$; and substituting $\lambda=2$: $8-24+20-4=0$ $\checkmark$.
>
> Matrices of this "$2$ on the diagonal, $1$ off it" type (with sign changes) arise whenever identical elements are chained — an $LC$ ladder, a discretised transmission line, a finite-difference Laplacian ([Ch. 13](../part-4-differential-equations/13-partial-differential-equations.md)). Their eigenvalues are known in closed form: for the $n\times n$ version, $\lambda_k=2+2\cos\frac{k\pi}{n+1}$, $k=1,\dots,n$, which for $n=3$ gives $2+\sqrt2$, $2$, $2-\sqrt2$.

> [!example] Worked Example 10.3 — The companion matrix of a third-order ODE
> The differential equation $\dddot y+6\ddot y+11\dot y+6y=0$ is put into state form with $x_1=y$, $x_2=\dot y$, $x_3=\ddot y$:
> $$
> \dot{\mathbf{x}}=\begin{bmatrix}0&1&0\\0&0&1\\-6&-11&-6\end{bmatrix}\mathbf{x}=\mathbf{C}\mathbf{x}.
> $$
> **Step 1 — characteristic polynomial.** Expand $\det(\lambda\mathbf{I}-\mathbf{C})$ along the first column:
> $$
> \det\begin{bmatrix}\lambda&-1&0\\0&\lambda&-1\\6&11&\lambda+6\end{bmatrix}
> =\lambda\big[\lambda(\lambda+6)+11\big]+6\big[(-1)(-1)-0\cdot\lambda\big]
> =\lambda^3+6\lambda^2+11\lambda+6 .
> $$
> This is precisely the characteristic polynomial you would obtain by substituting $y=e^{\lambda t}$ into the ODE.
>
> **Step 2 — roots.** Trying $\lambda=-1$: $-1+6-11+6=0$. Dividing, $\lambda^3+6\lambda^2+11\lambda+6=(\lambda+1)(\lambda^2+5\lambda+6)=(\lambda+1)(\lambda+2)(\lambda+3)$.
>
> The eigenvalues $-1,-2,-3$ are the natural frequencies (in s$^{-1}$, if $t$ is in seconds) of the system; all are negative, so every solution decays (Section 10.9). **Check:** $\operatorname{tr}\mathbf{C}=-6=-1-2-3$ $\checkmark$, $\det\mathbf{C}=-6=(-1)(-2)(-3)$ $\checkmark$.

> [!example] Worked Example 10.4 — Inverse and powers by Cayley–Hamilton
> For $\mathbf{A}=\begin{bmatrix}4&1\\2&3\end{bmatrix}$ of Example 10.1, Cayley–Hamilton gives $\mathbf{A}^2-7\mathbf{A}+10\mathbf{I}=\mathbf{0}$.
>
> **Inverse.** Rearranging, $\mathbf{A}(7\mathbf{I}-\mathbf{A})=10\mathbf{I}$, so
> $$
> \mathbf{A}^{-1}=\frac{1}{10}(7\mathbf{I}-\mathbf{A})=\frac1{10}\begin{bmatrix}3&-1\\-2&4\end{bmatrix}.
> $$
> Check: $\mathbf{A}\mathbf{A}^{-1}=\frac1{10}\begin{bmatrix}12-2&-4+4\\6-6&-2+12\end{bmatrix}=\mathbf{I}$ $\checkmark$.
>
> **Cube.** $\mathbf{A}^2=7\mathbf{A}-10\mathbf{I}$, so $\mathbf{A}^3=7\mathbf{A}^2-10\mathbf{A}=7(7\mathbf{A}-10\mathbf{I})-10\mathbf{A}=39\mathbf{A}-70\mathbf{I}$:
> $$
> \mathbf{A}^3=\begin{bmatrix}156-70&39\\78&117-70\end{bmatrix}=\begin{bmatrix}86&39\\78&47\end{bmatrix},
> $$
> which agrees with direct multiplication. Every power $\mathbf{A}^k$ is of the form $\alpha_k\mathbf{A}+\beta_k\mathbf{I}$; Section 10.4 finds $\alpha_k,\beta_k$ in closed form.

The characteristic polynomial of Example 10.2 is worth seeing as a graph. Its three real roots are where the cubic crosses the axis; for a symmetric matrix (Section 10.5) such crossings are guaranteed to be real.

```desmos-graph
left=-0.5; right=4.5; bottom=-3; top=3
---
y=-x^{3}+6x^{2}-10x+4|blue
(0.5858,0)|label:2−√2
(2,0)|label:2
(3.4142,0)|label:2+√2
(0,4)|label:p(0)=det B=4
```
*The characteristic polynomial $p(\lambda)=\det(\mathbf{B}-\lambda\mathbf{I})=-\lambda^3+6\lambda^2-10\lambda+4$ of Example 10.2 (horizontal axis $\lambda$); its zeros are the eigenvalues and its intercept is $\det\mathbf{B}$.*

> [!warning] Common Pitfall — subtracting $\lambda$ everywhere, or expanding too early
> Two errors account for most wrong characteristic polynomials. (i) $\lambda$ is subtracted from the **diagonal entries only**; writing $a_{12}-\lambda$ is wrong. (ii) Multiplying everything out and then hunting for roots of a cubic is laborious and error-prone; look for a common factor first (Example 10.2) or row-reduce cleverly. Always finish with the trace and determinant checks — they cost seconds and catch almost every arithmetic slip.

> [!warning] Common Pitfall — row reduction changes eigenvalues
> Gaussian elimination preserves the *solution set* of $\mathbf{A}\mathbf{x}=\mathbf{b}$ but not the eigenvalues of $\mathbf{A}$. The echelon form of $\begin{bmatrix}4&1\\2&3\end{bmatrix}$ is $\begin{bmatrix}4&1\\0&2.5\end{bmatrix}$, with diagonal $4,2.5$ — not the eigenvalues $5,2$. You may row-reduce $\mathbf{A}-\lambda\mathbf{I}$ to find *eigenvectors* once $\lambda$ is fixed, but never row-reduce $\mathbf{A}$ to read off eigenvalues. Only *similarity* transformations $\mathbf{S}^{-1}\mathbf{A}\mathbf{S}$ preserve them.

> [!tip] Engineering Note — never compute large characteristic polynomials numerically
> For hand work with $n\le3$ the characteristic polynomial is the right tool. For $n$ in the tens or thousands — a power-system network, a finite-element model, a circuit simulator's Jacobian — the roots of a polynomial are extraordinarily sensitive to rounding of its coefficients (Wilkinson's famous example: perturbing one coefficient of $\prod_{k=1}^{20}(\lambda-k)$ by $2^{-23}$ makes several roots complex). Production software (LAPACK, MATLAB `eig`, NumPy) never forms $p(\lambda)$; it applies a sequence of orthogonal similarity transformations — the QR algorithm — that drive $\mathbf{A}$ to triangular form, whose diagonal holds the eigenvalues. The method is treated in [Ch. 22](../part-7-numerical-computational/22-numerical-methods.md); this chapter is about what the answers mean.


## 10.2 Eigenvalues and eigenvectors

### Motivation

Once the eigenvalues are known, the eigenvectors tell us *where* the matrix acts as pure scaling. Geometrically: a $2\times2$ matrix turns the unit circle into an ellipse and generally rotates every direction; an eigenvector is a direction that comes out parallel to itself — stretched by $\lambda$, shrunk, or reversed if $\lambda<0$. Physically: in a vibrating structure an eigenvector is a *mode shape*, the pattern of relative amplitudes with which all parts move together at one frequency; in a Markov chain the eigenvector for $\lambda=1$ is the long-run probability distribution; in an antenna array or a MIMO channel the dominant eigenvector of the correlation matrix is the best beam-forming direction.

### Definitions

> [!definition] Eigenvalue, eigenvector, eigenspace
> Let $\mathbf{A}$ be $n\times n$. A scalar $\lambda\in\mathbb{C}$ is an **eigenvalue** of $\mathbf{A}$ if there is a **nonzero** vector $\mathbf{v}\in\mathbb{C}^n$ with
> $$
> \mathbf{A}\mathbf{v}=\lambda\mathbf{v}.
> $$
> Such a $\mathbf{v}$ is an **eigenvector** belonging to $\lambda$, and $(\lambda,\mathbf{v})$ is an **eigenpair**. The set
> $$
> E_\lambda=\operatorname{null}(\mathbf{A}-\lambda\mathbf{I})=\{\mathbf{v}:\mathbf{A}\mathbf{v}=\lambda\mathbf{v}\}
> $$
> is the **eigenspace** of $\lambda$; it consists of all eigenvectors for $\lambda$ together with $\mathbf{0}$. Its dimension $g=\dim E_\lambda=n-\operatorname{rank}(\mathbf{A}-\lambda\mathbf{I})$ is the **geometric multiplicity** of $\lambda$.

Three remarks clarify the definition.

1. **The zero vector is excluded as an eigenvector** (otherwise every $\lambda$ would be an eigenvalue), but $\lambda=0$ is a perfectly good eigenvalue: it means $\mathbf{A}\mathbf{v}=\mathbf{0}$ for some $\mathbf{v}\ne\mathbf{0}$, i.e. $\mathbf{A}$ is singular and $\mathbf{v}$ lies in its null space.
2. **Eigenvectors are never unique.** If $\mathbf{v}$ is an eigenvector, so is $c\mathbf{v}$ for every $c\ne0$ (including complex $c$). An eigenvector really specifies a *direction* (a line through the origin), and when the eigenspace has dimension $g>1$, any nonzero vector in that $g$-dimensional subspace qualifies. Different software normalises differently — unit length, largest entry $1$, first entry $1$ — so do not be alarmed when your hand answer differs from MATLAB's by a scale factor or a sign.
3. **The eigenspace is a subspace** ([Ch. 9](09-vector-spaces-and-linear-transformations.md)): it is the null space of a matrix, so it is closed under addition and scaling.

### Finding eigenvectors

The procedure is a direct application of null-space computation:

1. Solve $\det(\mathbf{A}-\lambda\mathbf{I})=0$ for the eigenvalues.
2. For each distinct eigenvalue $\lambda_i$, row-reduce $\mathbf{A}-\lambda_i\mathbf{I}$. Because $\lambda_i$ makes this matrix singular, at least one row must reduce to zero — **if it doesn't, your eigenvalue is wrong.**
3. Express the basic variables in terms of the free variables; one free variable per dimension of $E_{\lambda_i}$. Setting each free variable to $1$ in turn (others $0$) gives a basis of the eigenspace.
4. Check $\mathbf{A}\mathbf{v}=\lambda\mathbf{v}$ by direct multiplication.

For a $2\times2$ matrix step 2 is trivial. If $\mathbf{A}-\lambda\mathbf{I}=\begin{bmatrix}a&b\\c&d\end{bmatrix}$ is singular, its rows are proportional, and $\mathbf{v}=\begin{bmatrix}b\\-a\end{bmatrix}$ (or $\begin{bmatrix}d\\-c\end{bmatrix}$) is an eigenvector whenever it is nonzero — it is perpendicular to the first row.

### Algebraic versus geometric multiplicity

> [!theorem] Multiplicity inequality
> For every eigenvalue $\lambda_0$ of $\mathbf{A}$,
> $$
> 1\le g(\lambda_0)\le m(\lambda_0),
> $$
> where $g$ is the geometric and $m$ the algebraic multiplicity.

*Proof sketch.* $g\ge1$ because $\mathbf{A}-\lambda_0\mathbf{I}$ is singular. For the upper bound, take a basis $\mathbf{v}_1,\dots,\mathbf{v}_g$ of $E_{\lambda_0}$ and extend it to a basis of $\mathbb{C}^n$; put these vectors as the columns of $\mathbf{S}$. In the new basis the first $g$ columns of $\mathbf{S}^{-1}\mathbf{A}\mathbf{S}$ are $\lambda_0\mathbf{e}_1,\dots,\lambda_0\mathbf{e}_g$, so $\mathbf{S}^{-1}\mathbf{A}\mathbf{S}=\begin{bmatrix}\lambda_0\mathbf{I}_g&\ast\\\mathbf{0}&\mathbf{C}\end{bmatrix}$. This block-triangular matrix has characteristic polynomial $(\lambda_0-\lambda)^g\det(\mathbf{C}-\lambda\mathbf{I})$, and similar matrices share characteristic polynomials, so $(\lambda_0-\lambda)^g$ divides $p(\lambda)$: $g\le m$. $\square$

An eigenvalue with $g<m$ is called **defective**, and so is a matrix having one. The standard example is the **shear**

$$
\mathbf{N}=\begin{bmatrix}1&1\\0&1\end{bmatrix}:\qquad p(\lambda)=(1-\lambda)^2,\quad m=2,\qquad \mathbf{N}-\mathbf{I}=\begin{bmatrix}0&1\\0&0\end{bmatrix},\quad g=1 .
$$

The only eigenvector direction is $\mathbf{e}_1=(1,0)$: the shear slides every horizontal line along itself, so only horizontal vectors keep their direction. Defective matrices are exactly the ones that cannot be diagonalized (Section 10.3); in dynamics they produce the $te^{\lambda t}$ terms familiar from critically damped circuits.

### Independence of eigenvectors

> [!theorem] Eigenvectors for distinct eigenvalues are independent
> If $\mathbf{v}_1,\dots,\mathbf{v}_k$ are eigenvectors of $\mathbf{A}$ belonging to **distinct** eigenvalues $\lambda_1,\dots,\lambda_k$, then $\{\mathbf{v}_1,\dots,\mathbf{v}_k\}$ is linearly independent.

*Proof.* Suppose not, and let $r$ be the smallest number such that $\mathbf{v}_{r+1}$ is a combination of the independent set $\mathbf{v}_1,\dots,\mathbf{v}_r$: $\mathbf{v}_{r+1}=c_1\mathbf{v}_1+\dots+c_r\mathbf{v}_r$. Apply $\mathbf{A}$: $\lambda_{r+1}\mathbf{v}_{r+1}=c_1\lambda_1\mathbf{v}_1+\dots+c_r\lambda_r\mathbf{v}_r$. Subtract $\lambda_{r+1}$ times the first relation:

$$
\mathbf{0}=c_1(\lambda_1-\lambda_{r+1})\mathbf{v}_1+\dots+c_r(\lambda_r-\lambda_{r+1})\mathbf{v}_r .
$$

Since $\mathbf{v}_1,\dots,\mathbf{v}_r$ are independent and every $\lambda_i-\lambda_{r+1}\ne0$, all $c_i=0$, so $\mathbf{v}_{r+1}=\mathbf{0}$ — impossible for an eigenvector. $\square$

The same argument shows that if we take bases of several eigenspaces for distinct eigenvalues and pool them, the pooled set is still independent. This is the engine of the diagonalization theorem in Section 10.3.

### Eigenvalues of related matrices

If $\mathbf{A}\mathbf{v}=\lambda\mathbf{v}$, then repeated application gives $\mathbf{A}^k\mathbf{v}=\lambda^k\mathbf{v}$, and linear combinations give the following table. The eigenvector $\mathbf{v}$ is the **same** in every row except the last.

| Matrix | Eigenvalue | Condition |
|---|---|---|
| $\mathbf{A}^k$, $k=1,2,\dots$ | $\lambda^k$ | — |
| $\mathbf{A}^{-1}$ | $1/\lambda$ | $\mathbf{A}$ invertible ($\lambda\ne0$) |
| $\mathbf{A}+c\mathbf{I}$ | $\lambda+c$ | — (spectral shift) |
| $c\mathbf{A}$ | $c\lambda$ | — |
| $q(\mathbf{A})=\sum_k q_k\mathbf{A}^k$ | $q(\lambda)$ | any polynomial $q$ |
| $(\mathbf{A}-\mu\mathbf{I})^{-1}$ | $1/(\lambda-\mu)$ | $\mu$ not an eigenvalue |
| $e^{\mathbf{A}t}$ | $e^{\lambda t}$ | (Section 10.4) |
| $\mathbf{A}^\mathsf{T}$ | $\lambda$ | eigenvector generally different |

The shifted-inverse row is the basis of the *inverse iteration* algorithm: choosing $\mu$ close to a wanted eigenvalue makes $1/(\lambda-\mu)$ huge, so that eigenvalue dominates. The last row has an important companion: vectors with $\mathbf{w}^\mathsf{T}\mathbf{A}=\lambda\mathbf{w}^\mathsf{T}$, i.e. eigenvectors of $\mathbf{A}^\mathsf{T}$, are called **left eigenvectors**. A left eigenvector for $\lambda_i$ is orthogonal to every right eigenvector for a different eigenvalue $\lambda_j$, since $\lambda_i\mathbf{w}_i^\mathsf{T}\mathbf{v}_j=\mathbf{w}_i^\mathsf{T}\mathbf{A}\mathbf{v}_j=\lambda_j\mathbf{w}_i^\mathsf{T}\mathbf{v}_j$ forces $\mathbf{w}_i^\mathsf{T}\mathbf{v}_j=0$. Left eigenvectors reappear in modal analysis (Section 10.8) as the tools that extract modal amplitudes.

### Complex eigenvalues of real matrices

When a real matrix has a complex eigenvalue $\lambda=\sigma+j\omega$ ($\omega\ne0$), its eigenvector is necessarily complex, and conjugating $\mathbf{A}\mathbf{v}=\lambda\mathbf{v}$ (with $\mathbf{A}$ real) gives $\mathbf{A}\mathbf{v}^*=\lambda^*\mathbf{v}^*$: **the conjugate eigenvalue has the conjugate eigenvector**, so only one of each pair needs to be computed. Geometrically, a complex pair signals *rotation*: no real direction is preserved. The prototype is

$$
\mathbf{A}=\begin{bmatrix}\sigma&-\omega\\\omega&\sigma\end{bmatrix}=r\begin{bmatrix}\cos\theta&-\sin\theta\\\sin\theta&\cos\theta\end{bmatrix},
\qquad r=\sqrt{\sigma^2+\omega^2},\ \theta=\angle(\sigma+j\omega),
$$

with eigenvalues $\sigma\pm j\omega=re^{\pm j\theta}$ and eigenvectors $(1,\mp j)$. This matrix is multiplication by the complex number $re^{j\theta}$ written in real coordinates ([Ch. 2](../part-1-foundations/02-complex-numbers.md)): rotate by $\theta$, scale by $r$. Section 10.3 shows that *every* real $2\times2$ matrix with complex eigenvalues is similar to this rotation-scaling form.

### Worked examples

> [!example] Worked Example 10.5 — Eigenvectors of a $2\times2$ matrix
> Find the eigenvectors of $\mathbf{A}=\begin{bmatrix}4&1\\2&3\end{bmatrix}$ (eigenvalues $5$ and $2$, Example 10.1).
>
> **$\lambda_1=5$.** $\mathbf{A}-5\mathbf{I}=\begin{bmatrix}-1&1\\2&-2\end{bmatrix}$. The second row is $-2$ times the first, as it must be. The first row says $-v_1+v_2=0$, so $v_2=v_1$:
> $$
> \mathbf{v}_1=\begin{bmatrix}1\\1\end{bmatrix},\qquad \mathbf{A}\mathbf{v}_1=\begin{bmatrix}5\\5\end{bmatrix}=5\mathbf{v}_1\ \checkmark
> $$
> **$\lambda_2=2$.** $\mathbf{A}-2\mathbf{I}=\begin{bmatrix}2&1\\2&1\end{bmatrix}$, giving $2v_1+v_2=0$:
> $$
> \mathbf{v}_2=\begin{bmatrix}1\\-2\end{bmatrix},\qquad \mathbf{A}\mathbf{v}_2=\begin{bmatrix}2\\-4\end{bmatrix}=2\mathbf{v}_2\ \checkmark
> $$
> The eigenvectors are independent (distinct eigenvalues) but **not** orthogonal: $\mathbf{v}_1\cdot\mathbf{v}_2=1-2=-1$. Orthogonality is a privilege of symmetric matrices (Section 10.5).

The graph below shows the geometry. The matrix maps the unit circle onto an ellipse: the image of $\lVert\mathbf{x}\rVert=1$ is the set of $\mathbf{y}$ with $\lVert\mathbf{A}^{-1}\mathbf{y}\rVert=1$, i.e. $\mathbf{y}^\mathsf{T}(\mathbf{A}\mathbf{A}^\mathsf{T})^{-1}\mathbf{y}=1$, which here is $13x^2-22xy+17y^2=100$. Points of the circle on the eigenvector lines land on the same lines, scaled by $5$ and $2$; the point $(1,0)$ is rotated off its line to $(4,2)$.

```desmos-graph
left=-5; right=5; bottom=-4; top=4
---
x^2+y^2=1|grey|dashed
13x^2-22xy+17y^2=100|blue
y=x|red|dotted
y=-2x|green|dotted
(0.7071,0.7071)|label:v1 (unit)
(3.5355,3.5355)|label:A v1 = 5 v1
(0.4472,-0.8944)|label:v2 (unit)
(0.8944,-1.7889)|label:A v2 = 2 v2
(1,0)|label:e1
(4,2)|label:A e1
```
*The unit circle (dashed) and its image under $\mathbf{A}=\begin{bmatrix}4&1\\2&3\end{bmatrix}$ (blue ellipse). Vectors on the eigenvector lines $y=x$ (red, $\lambda=5$) and $y=-2x$ (green, $\lambda=2$) are only stretched; every other direction, such as $\mathbf{e}_1\mapsto(4,2)$, is also turned. The eigenvector lines are not the ellipse's axes — that happens only for symmetric matrices.*

> [!example] Worked Example 10.6 — A repeated eigenvalue with a two-dimensional eigenspace
> Find all eigenvalues and eigenspaces of
> $$
> \mathbf{L}=\begin{bmatrix}4&-1&6\\2&1&6\\2&-1&8\end{bmatrix}.
> $$
> **Step 1 — eigenvalues.** $\operatorname{tr}\mathbf{L}=13$, $E_2=(4+2)+(32-12)+(8+6)=40$, $\det\mathbf{L}=4(8+6)+1(16-12)+6(-2-2)=56+4-24=36$. Hence
> $$
> \det(\lambda\mathbf{I}-\mathbf{L})=\lambda^3-13\lambda^2+40\lambda-36 .
> $$
> Trying $\lambda=2$: $8-52+80-36=0$. Dividing out, $\lambda^3-13\lambda^2+40\lambda-36=(\lambda-2)(\lambda^2-11\lambda+18)=(\lambda-2)^2(\lambda-9)$. So $\lambda=2$ with $m=2$ and $\lambda=9$ with $m=1$.
>
> **Step 2 — eigenspace for $\lambda=2$.**
> $$
> \mathbf{L}-2\mathbf{I}=\begin{bmatrix}2&-1&6\\2&-1&6\\2&-1&6\end{bmatrix}\ \longrightarrow\ \begin{bmatrix}2&-1&6\\0&0&0\\0&0&0\end{bmatrix}.
> $$
> Rank $1$, so $g=3-1=2$. The single equation $2v_1-v_2+6v_3=0$ gives $v_2=2v_1+6v_3$ with $v_1,v_3$ free. Choosing $(v_1,v_3)=(1,0)$ and $(0,1)$:
> $$
> E_2=\operatorname{span}\left\{\begin{bmatrix}1\\2\\0\end{bmatrix},\begin{bmatrix}0\\6\\1\end{bmatrix}\right\}.
> $$
> (Equally valid: $(1,2,0)$ and $(-3,0,1)$ — any two independent solutions of the plane equation.)
>
> **Step 3 — eigenspace for $\lambda=9$.**
> $$
> \mathbf{L}-9\mathbf{I}=\begin{bmatrix}-5&-1&6\\2&-8&6\\2&-1&-1\end{bmatrix}.
> $$
> The third row gives $v_2=2v_1-v_3$; substituting into the first, $-5v_1-2v_1+v_3+6v_3=0$, so $v_3=v_1$ and then $v_2=v_1$. Thus $E_9=\operatorname{span}\{(1,1,1)\}$. Check with the second row: $2-8+6=0$ $\checkmark$, and $\mathbf{L}(1,1,1)^\mathsf{T}=(9,9,9)^\mathsf{T}$ $\checkmark$.
>
> **Conclusion.** Here $g(2)=m(2)=2$, so the repeated eigenvalue is *not* defective, and the three eigenvectors $(1,2,0)$, $(0,6,1)$, $(1,1,1)$ form a basis of $\mathbb{R}^3$. Section 10.3 uses this to diagonalize $\mathbf{L}$. Contrast the shear $\mathbf{N}$ above, where the double eigenvalue had only a one-dimensional eigenspace.

> [!example] Worked Example 10.7 — Complex eigenvalues and eigenvectors
> Find the eigenpairs of $\mathbf{R}=\begin{bmatrix}1&-2\\1&3\end{bmatrix}$.
>
> **Step 1 — eigenvalues.** $\operatorname{tr}\mathbf{R}=4$, $\det\mathbf{R}=3+2=5$, so $\lambda^2-4\lambda+5=0$ and
> $$
> \lambda=2\pm\sqrt{4-5}=2\pm j .
> $$
> **Step 2 — eigenvector for $\lambda=2+j$.**
> $$
> \mathbf{R}-(2+j)\mathbf{I}=\begin{bmatrix}-1-j&-2\\1&1-j\end{bmatrix}.
> $$
> Use the simpler second row: $v_1+(1-j)v_2=0$. With $v_2=1$, $v_1=-1+j$, so
> $$
> \mathbf{v}=\begin{bmatrix}-1+j\\1\end{bmatrix}.
> $$
> Check the first row, which must also vanish: $(-1-j)(-1+j)-2=(1+1)-2=0$ $\checkmark$.
>
> **Step 3 — the conjugate pair.** For $\lambda^*=2-j$ the eigenvector is $\mathbf{v}^*=(-1-j,\,1)$, with no further work.
>
> **Interpretation.** $\lvert\lambda\rvert=\sqrt5$ and $\angle\lambda=\arctan(1/2)=26.57^\circ$: in a suitable (non-orthogonal) coordinate system, $\mathbf{R}$ rotates by $26.57^\circ$ and scales by $\sqrt5\approx2.236$. In continuous-time dynamics $\dot{\mathbf{x}}=\mathbf{R}\mathbf{x}$, the pair $2\pm j$ means oscillation at $1$ rad/s inside a growing envelope $e^{2t}$ — an unstable spiral (Section 10.9).

> [!example] Worked Example 10.8 — Eigenvalues of functions of a matrix
> With $\mathbf{A}=\begin{bmatrix}4&1\\2&3\end{bmatrix}$ (eigenvalues $5,2$), find the eigenvalues of $\mathbf{B}=\mathbf{A}^2-3\mathbf{A}+\mathbf{I}$, of $\mathbf{A}^{-1}$ and of $\mathbf{A}-6\mathbf{I}$, and decide whether $\mathbf{B}$ is invertible — without computing $\mathbf{B}$.
>
> - $q(\lambda)=\lambda^2-3\lambda+1$: $q(5)=11$, $q(2)=-1$. Eigenvalues of $\mathbf{B}$: $11,-1$. Since neither is zero, $\mathbf{B}$ is invertible, with $\det\mathbf{B}=-11$.
> - $\mathbf{A}^{-1}$: $1/5,\ 1/2$ (compare Example 10.4: $\operatorname{tr}\mathbf{A}^{-1}=(3+4)/10=0.7=0.2+0.5$ $\checkmark$).
> - $\mathbf{A}-6\mathbf{I}$: $-1,\,-4$. Shifting moves the whole spectrum without changing any eigenvector.
>
> **Direct check for $\mathbf{B}$.** $\mathbf{A}^2=\begin{bmatrix}18&7\\14&11\end{bmatrix}$, so $\mathbf{B}=\begin{bmatrix}18-12+1&7-3\\14-6&11-9+1\end{bmatrix}=\begin{bmatrix}7&4\\8&3\end{bmatrix}$; $\operatorname{tr}\mathbf{B}=10=11+(-1)$ $\checkmark$, $\det\mathbf{B}=21-32=-11$ $\checkmark$. And $\mathbf{B}(1,1)^\mathsf{T}=(11,11)^\mathsf{T}$: same eigenvector as for $\mathbf{A}$ $\checkmark$.

> [!warning] Common Pitfall — eigenvalues of sums and products
> The table above concerns functions of a **single** matrix. In general the eigenvalues of $\mathbf{A}+\mathbf{B}$ are **not** sums of eigenvalues of $\mathbf{A}$ and $\mathbf{B}$, and those of $\mathbf{A}\mathbf{B}$ are not products — unless $\mathbf{A}$ and $\mathbf{B}$ share a full set of eigenvectors (which happens, for diagonalizable matrices, exactly when $\mathbf{A}\mathbf{B}=\mathbf{B}\mathbf{A}$). Example: $\begin{bmatrix}0&1\\0&0\end{bmatrix}$ and $\begin{bmatrix}0&0\\1&0\end{bmatrix}$ both have only the eigenvalue $0$, but their sum $\begin{bmatrix}0&1\\1&0\end{bmatrix}$ has eigenvalues $\pm1$. What *is* true: $\operatorname{tr}(\mathbf{A}+\mathbf{B})=\operatorname{tr}\mathbf{A}+\operatorname{tr}\mathbf{B}$, and $\mathbf{A}\mathbf{B}$ and $\mathbf{B}\mathbf{A}$ have the same eigenvalues.

> [!tip] Engineering Note — what eigenvalues mean physically
> When $\mathbf{A}$ maps a physical state to its rate of change ($\dot{\mathbf{x}}=\mathbf{A}\mathbf{x}$), the eigenvalues have units of s$^{-1}$: a real eigenvalue $-1/\tau$ is a decay rate with time constant $\tau$, and a complex pair $-\alpha\pm j\omega_d$ is a damped oscillation at $\omega_d$ rad/s. When $\mathbf{A}$ maps a state at one sample to the next ($\mathbf{x}[k+1]=\mathbf{A}\mathbf{x}[k]$), the eigenvalues are dimensionless multipliers per step. When $\mathbf{A}$ is a covariance matrix of voltages, its eigenvalues are in V$^2$ — the variance (power) carried along each eigenvector direction. Always ask what $\mathbf{A}$ maps from and to; the units of the eigenvalues follow.

## 10.3 Diagonalization

### Motivation

A diagonal matrix is the simplest matrix there is: it scales the $i$-th coordinate by $d_i$ and leaves the coordinates uncoupled. Its powers, inverse, exponential and every other function are computed entry by entry. A general matrix couples the coordinates — the rate of change of the capacitor voltage depends on the inductor current and vice versa — and that coupling is what makes systems hard to analyse. **Diagonalization** removes the coupling by changing to a basis of eigenvectors. In that basis the matrix is diagonal with the eigenvalues on the diagonal, and every computation becomes $n$ scalar computations. Afterwards we change back.

This is the matrix version of a strategy used all through this book: transform a hard problem into a domain where it is easy (phasors, Laplace, Fourier), solve it there, and transform back. Here the "transform" is multiplication by $\mathbf{P}^{-1}$, the "easy domain" is modal coordinates, and the "inverse transform" is multiplication by $\mathbf{P}$.

### The factorization

Suppose $\mathbf{A}$ has $n$ linearly independent eigenvectors $\mathbf{v}_1,\dots,\mathbf{v}_n$ with eigenvalues $\lambda_1,\dots,\lambda_n$ (not necessarily distinct). Place the eigenvectors as columns of $\mathbf{P}$ and the eigenvalues, **in the same order**, on the diagonal of $\mathbf{D}$. Then, column by column,

$$
\mathbf{A}\mathbf{P}=\mathbf{A}\begin{bmatrix}\mathbf{v}_1&\cdots&\mathbf{v}_n\end{bmatrix}
=\begin{bmatrix}\lambda_1\mathbf{v}_1&\cdots&\lambda_n\mathbf{v}_n\end{bmatrix}
=\begin{bmatrix}\mathbf{v}_1&\cdots&\mathbf{v}_n\end{bmatrix}\begin{bmatrix}\lambda_1&&\\&\ddots&\\&&\lambda_n\end{bmatrix}=\mathbf{P}\mathbf{D}.
$$

Since the columns of $\mathbf{P}$ are independent, $\mathbf{P}$ is invertible, and

$$
\mathbf{A}=\mathbf{P}\mathbf{D}\mathbf{P}^{-1},\qquad \mathbf{D}=\mathbf{P}^{-1}\mathbf{A}\mathbf{P}.
$$

Conversely, if $\mathbf{A}=\mathbf{P}\mathbf{D}\mathbf{P}^{-1}$ with $\mathbf{D}$ diagonal, then $\mathbf{A}\mathbf{P}=\mathbf{P}\mathbf{D}$ says exactly that each column of $\mathbf{P}$ is an eigenvector with eigenvalue the corresponding diagonal entry of $\mathbf{D}$.

> [!definition] Diagonalizable matrix
> A square matrix $\mathbf{A}$ is **diagonalizable** if it is similar to a diagonal matrix: there exist an invertible $\mathbf{P}$ and a diagonal $\mathbf{D}$ with $\mathbf{A}=\mathbf{P}\mathbf{D}\mathbf{P}^{-1}$. The columns of $\mathbf{P}$ form an **eigenbasis**.

> [!theorem] Diagonalization theorem
> For an $n\times n$ matrix $\mathbf{A}$ the following are equivalent:
> 1. $\mathbf{A}$ is diagonalizable;
> 2. $\mathbf{A}$ has $n$ linearly independent eigenvectors;
> 3. for every eigenvalue, the geometric multiplicity equals the algebraic multiplicity, $g(\lambda_i)=m(\lambda_i)$;
> 4. the dimensions of the eigenspaces add up to $n$.
>
> In particular, **an $n\times n$ matrix with $n$ distinct eigenvalues is diagonalizable.**

*Proof sketch.* (1)$\Leftrightarrow$(2) is the computation above. For (2)$\Leftrightarrow$(3)$\Leftrightarrow$(4): the sum of algebraic multiplicities is always $n$ and $g_i\le m_i$, so $\sum g_i=n$ if and only if $g_i=m_i$ for all $i$. Pooling bases of the eigenspaces gives an independent set (Section 10.2) with $\sum g_i$ vectors, and no independent set of eigenvectors can be larger. The final statement follows because each distinct eigenvalue has $g=m=1$. $\square$

The "distinct eigenvalues" criterion is sufficient but not necessary: the identity matrix has a single eigenvalue repeated $n$ times yet is already diagonal, and $\mathbf{L}$ of Example 10.6 has a double eigenvalue but is diagonalizable.

### Interpretation as a change of basis

Read $\mathbf{A}\mathbf{x}=\mathbf{P}\mathbf{D}\mathbf{P}^{-1}\mathbf{x}$ from right to left:

1. $\mathbf{y}=\mathbf{P}^{-1}\mathbf{x}$ — find the coordinates of $\mathbf{x}$ in the eigenbasis, i.e. write $\mathbf{x}=y_1\mathbf{v}_1+\dots+y_n\mathbf{v}_n$ ([Ch. 9](09-vector-spaces-and-linear-transformations.md));
2. $\mathbf{D}\mathbf{y}$ — scale each coordinate by its eigenvalue: $y_i\mapsto\lambda_iy_i$;
3. $\mathbf{P}(\mathbf{D}\mathbf{y})$ — reassemble: $\mathbf{A}\mathbf{x}=\lambda_1y_1\mathbf{v}_1+\dots+\lambda_ny_n\mathbf{v}_n$.

The same reading applied to the linear system $\dot{\mathbf{x}}=\mathbf{A}\mathbf{x}$ is the whole of modal analysis. Substituting $\mathbf{x}=\mathbf{P}\mathbf{y}$ gives $\mathbf{P}\dot{\mathbf{y}}=\mathbf{A}\mathbf{P}\mathbf{y}$, i.e.

$$
\dot{\mathbf{y}}=\mathbf{P}^{-1}\mathbf{A}\mathbf{P}\,\mathbf{y}=\mathbf{D}\mathbf{y}
\quad\Longleftrightarrow\quad
\dot y_i=\lambda_iy_i,\quad i=1,\dots,n,
$$

$n$ uncoupled first-order equations with solutions $y_i(t)=y_i(0)e^{\lambda_it}$. The coordinates $y_i$ are the **modal coordinates** (Section 10.8).

### Diagonalizing with complex eigenvalues; the real block form

A real matrix with complex eigenvalues is diagonalizable over $\mathbb{C}$ (if its eigenvalues are distinct), with complex $\mathbf{P}$ and $\mathbf{D}$. Often we prefer to stay real. For a $2\times2$ real matrix with eigenvalues $\sigma\pm j\omega$ ($\omega\ne0$) and eigenvector $\mathbf{v}=\mathbf{a}+j\mathbf{b}$ for $\sigma+j\omega$ (with $\mathbf{a},\mathbf{b}$ real), separate $\mathbf{A}(\mathbf{a}+j\mathbf{b})=(\sigma+j\omega)(\mathbf{a}+j\mathbf{b})$ into real and imaginary parts:

$$
\mathbf{A}\mathbf{a}=\sigma\mathbf{a}-\omega\mathbf{b},\qquad \mathbf{A}\mathbf{b}=\omega\mathbf{a}+\sigma\mathbf{b}.
$$

With $\mathbf{P}=\begin{bmatrix}\mathbf{a}&\mathbf{b}\end{bmatrix}$ (real and invertible, since $\mathbf{a},\mathbf{b}$ are independent when $\omega\ne0$), these two lines say

$$
\mathbf{P}^{-1}\mathbf{A}\mathbf{P}=\begin{bmatrix}\sigma&\omega\\-\omega&\sigma\end{bmatrix}.
$$

This **real block form** is a rotation-scaling matrix. For an $n\times n$ real matrix, each complex pair contributes such a $2\times2$ block and each real eigenvalue a $1\times1$ block; the result is block-diagonal and entirely real. In state-space terms the two real coordinates of a block are the in-phase and quadrature components of an oscillating mode.

### When diagonalization fails: the Jordan form

A defective matrix has too few eigenvectors to form a basis. The best that similarity can achieve is the **Jordan canonical form**: $\mathbf{A}=\mathbf{P}\mathbf{J}\mathbf{P}^{-1}$ with $\mathbf{J}$ block-diagonal, each block of the form

$$
\mathbf{J}_k(\lambda)=\begin{bmatrix}\lambda&1&&\\&\lambda&\ddots&\\&&\ddots&1\\&&&\lambda\end{bmatrix}\quad(k\times k).
$$

The columns of $\mathbf{P}$ for a $2\times2$ Jordan block are an eigenvector $\mathbf{v}$ and a **generalized eigenvector** $\mathbf{w}$ satisfying

$$
(\mathbf{A}-\lambda\mathbf{I})\mathbf{v}=\mathbf{0},\qquad(\mathbf{A}-\lambda\mathbf{I})\mathbf{w}=\mathbf{v}.
$$

Then $\mathbf{A}\mathbf{v}=\lambda\mathbf{v}$ and $\mathbf{A}\mathbf{w}=\mathbf{v}+\lambda\mathbf{w}$, which is the column statement of $\mathbf{A}\mathbf{P}=\mathbf{P}\mathbf{J}_2(\lambda)$. Jordan forms are numerically fragile — an arbitrarily small perturbation makes a defective matrix diagonalizable with nearly parallel eigenvectors — so they are a theoretical tool rather than a computational one. Their practical message is the appearance of polynomial-times-exponential terms, $te^{\lambda t}$ and $k\lambda^{k-1}$ (Section 10.4).

### Simultaneous diagonalization

Two diagonalizable matrices $\mathbf{A}$ and $\mathbf{B}$ can be diagonalized by the **same** $\mathbf{P}$ if and only if they commute, $\mathbf{A}\mathbf{B}=\mathbf{B}\mathbf{A}$. One direction is immediate: diagonal matrices commute, so $\mathbf{P}\mathbf{D}_1\mathbf{P}^{-1}\mathbf{P}\mathbf{D}_2\mathbf{P}^{-1}=\mathbf{P}\mathbf{D}_2\mathbf{D}_1\mathbf{P}^{-1}$. This is why, for instance, all circulant matrices (which commute with the cyclic shift) are diagonalized by the same matrix — the DFT matrix — a fact that underlies fast convolution ([Ch. 16](../part-5-transforms/16-fourier-transforms.md), [Ch. 24](../part-7-numerical-computational/24-computational-mathematics.md)).

### Worked examples

> [!example] Worked Example 10.9 — Diagonalizing a $2\times2$ matrix
> Diagonalize $\mathbf{A}=\begin{bmatrix}4&1\\2&3\end{bmatrix}$.
>
> **Step 1 — eigenpairs** (Examples 10.1, 10.5): $(5,(1,1))$ and $(2,(1,-2))$.
>
> **Step 2 — assemble.**
> $$
> \mathbf{P}=\begin{bmatrix}1&1\\1&-2\end{bmatrix},\qquad \mathbf{D}=\begin{bmatrix}5&0\\0&2\end{bmatrix}.
> $$
> **Step 3 — invert $\mathbf{P}$.** $\det\mathbf{P}=-2-1=-3$, so
> $$
> \mathbf{P}^{-1}=\frac{1}{-3}\begin{bmatrix}-2&-1\\-1&1\end{bmatrix}=\frac13\begin{bmatrix}2&1\\1&-1\end{bmatrix}.
> $$
> **Step 4 — verify.**
> $$
> \mathbf{P}\mathbf{D}\mathbf{P}^{-1}=\begin{bmatrix}5&2\\5&-4\end{bmatrix}\cdot\frac13\begin{bmatrix}2&1\\1&-1\end{bmatrix}=\frac13\begin{bmatrix}12&3\\6&9\end{bmatrix}=\begin{bmatrix}4&1\\2&3\end{bmatrix}\ \checkmark
> $$
> **Remarks.** Reordering the columns of $\mathbf{P}$ (and correspondingly the diagonal of $\mathbf{D}$) or rescaling any column of $\mathbf{P}$ gives an equally valid factorization; software might return $\mathbf{P}=\begin{bmatrix}-1&1\\2&1\end{bmatrix}$, $\mathbf{D}=\operatorname{diag}(2,5)$. The rows of $\mathbf{P}^{-1}$, $(2,1)/3$ and $(1,-1)/3$, are left eigenvectors of $\mathbf{A}$: indeed $(2,1)\mathbf{A}=(10,5)=5(2,1)$.

> [!example] Worked Example 10.10 — A $3\times3$ matrix with a repeated eigenvalue
> Diagonalize $\mathbf{L}=\begin{bmatrix}4&-1&6\\2&1&6\\2&-1&8\end{bmatrix}$ from Example 10.6.
>
> **Step 1 — check diagonalizability.** $m(2)=2=g(2)$ and $m(9)=1=g(9)$, so the eigenspace dimensions sum to $3$: diagonalizable.
>
> **Step 2 — assemble** (eigenvectors in the order $2,2,9$):
> $$
> \mathbf{P}=\begin{bmatrix}1&0&1\\2&6&1\\0&1&1\end{bmatrix},\qquad\mathbf{D}=\begin{bmatrix}2&0&0\\0&2&0\\0&0&9\end{bmatrix}.
> $$
> **Step 3 — invert.** $\det\mathbf{P}=1(6-1)-0+1(2-0)=7$, and the adjugate gives
> $$
> \mathbf{P}^{-1}=\frac17\begin{bmatrix}5&1&-6\\-2&1&1\\2&-1&6\end{bmatrix}.
> $$
> **Step 4 — verify without inverting** (cheaper and just as conclusive): check $\mathbf{L}\mathbf{P}=\mathbf{P}\mathbf{D}$.
> $$
> \mathbf{L}\mathbf{P}=\begin{bmatrix}4-2&0-6+6&4-1+6\\2+2&0+6+6&2+1+6\\2-2&0-6+8&2-1+8\end{bmatrix}=\begin{bmatrix}2&0&9\\4&12&9\\0&2&9\end{bmatrix}=\mathbf{P}\mathbf{D}\ \checkmark
> $$
> (The columns of $\mathbf{P}\mathbf{D}$ are $2\times$, $2\times$ and $9\times$ the columns of $\mathbf{P}$.)

> [!example] Worked Example 10.11 — A defective matrix and its Jordan form
> Show that $\mathbf{A}=\begin{bmatrix}3&1\\-1&1\end{bmatrix}$ is not diagonalizable and find its Jordan form.
>
> **Step 1 — eigenvalues.** $\operatorname{tr}\mathbf{A}=4$, $\det\mathbf{A}=3+1=4$: $\lambda^2-4\lambda+4=(\lambda-2)^2$. Double eigenvalue $\lambda=2$, $m=2$.
>
> **Step 2 — eigenspace.** $\mathbf{A}-2\mathbf{I}=\begin{bmatrix}1&1\\-1&-1\end{bmatrix}$ has rank $1$, so $g=1<m$. Defective: not diagonalizable. Eigenvector $\mathbf{v}=(1,-1)$.
>
> **Step 3 — generalized eigenvector.** Solve $(\mathbf{A}-2\mathbf{I})\mathbf{w}=\mathbf{v}$: $w_1+w_2=1$. Take $\mathbf{w}=(1,0)$.
>
> **Step 4 — Jordan form.**
> $$
> \mathbf{P}=\begin{bmatrix}1&1\\-1&0\end{bmatrix},\qquad \mathbf{J}=\begin{bmatrix}2&1\\0&2\end{bmatrix},\qquad \mathbf{A}=\mathbf{P}\mathbf{J}\mathbf{P}^{-1}.
> $$
> Check $\mathbf{A}\mathbf{P}=\begin{bmatrix}2&3\\-2&-1\end{bmatrix}$ and $\mathbf{P}\mathbf{J}=\begin{bmatrix}2&1+2\\-2&-1+0\end{bmatrix}=\begin{bmatrix}2&3\\-2&-1\end{bmatrix}$ $\checkmark$.
>
> **Consequence.** The system $\dot{\mathbf{x}}=\mathbf{A}\mathbf{x}$ has solutions containing $te^{2t}$; indeed (Section 10.4) $e^{\mathbf{A}t}=e^{2t}\begin{bmatrix}1+t&t\\-t&1-t\end{bmatrix}$. In circuit language: a critically damped response.

> [!example] Worked Example 10.12 — Real block form of a rotation-scaling matrix
> Bring $\mathbf{R}=\begin{bmatrix}1&-2\\1&3\end{bmatrix}$ (Example 10.7, eigenvalues $2\pm j$) to real block form.
>
> **Step 1 — split the eigenvector.** For $\lambda=2+j$, $\mathbf{v}=(-1+j,\,1)=\underbrace{(-1,1)}_{\mathbf{a}}+j\underbrace{(1,0)}_{\mathbf{b}}$.
>
> **Step 2 — form $\mathbf{P}=[\mathbf{a}\ \mathbf{b}]$.** $\mathbf{P}=\begin{bmatrix}-1&1\\1&0\end{bmatrix}$, $\det\mathbf{P}=-1$, $\mathbf{P}^{-1}=\begin{bmatrix}0&1\\1&1\end{bmatrix}$.
>
> **Step 3 — compute.**
> $$
> \mathbf{R}\mathbf{P}=\begin{bmatrix}-3&1\\2&1\end{bmatrix},\qquad \mathbf{P}^{-1}\mathbf{R}\mathbf{P}=\begin{bmatrix}0&1\\1&1\end{bmatrix}\begin{bmatrix}-3&1\\2&1\end{bmatrix}=\begin{bmatrix}2&1\\-1&2\end{bmatrix}=\begin{bmatrix}\sigma&\omega\\-\omega&\sigma\end{bmatrix}
> $$
> with $\sigma=2$, $\omega=1$, as predicted. In the coordinates defined by $\mathbf{a}$ and $\mathbf{b}$, $\mathbf{R}$ is exactly "scale by $\sqrt5$, rotate by $26.57^\circ$" (clockwise, given this sign arrangement).

> [!warning] Common Pitfall — mismatched order, and "distinct" misread
> The $i$-th column of $\mathbf{P}$ must be the eigenvector for the $i$-th diagonal entry of $\mathbf{D}$; swapping one without the other produces a matrix that is not $\mathbf{A}$. Also, "distinct eigenvalues $\Rightarrow$ diagonalizable" is one-way: a repeated eigenvalue does **not** imply non-diagonalizable. Always compute $\operatorname{rank}(\mathbf{A}-\lambda\mathbf{I})$ for each repeated eigenvalue before deciding.

> [!tip] Engineering Note — nearly defective means sensitive
> Matrices that are close to defective have eigenvector matrices $\mathbf{P}$ that are close to singular (nearly parallel eigenvectors). The condition number $\kappa(\mathbf{P})$ then bounds how much eigenvalues can move under a perturbation of $\mathbf{A}$ (the Bauer–Fike theorem: each eigenvalue of $\mathbf{A}+\mathbf{E}$ lies within $\kappa(\mathbf{P})\lVert\mathbf{E}\rVert$ of an eigenvalue of $\mathbf{A}$). In control design, closed-loop poles placed with a nearly defective $\mathbf{A}$ are very sensitive to component tolerances, and large $\kappa(\mathbf{P})$ also predicts large transient overshoot before an eventually stable response decays (Section 10.9). Symmetric matrices, whose $\mathbf{P}$ can be orthogonal with $\kappa=1$, are the opposite extreme.

## 10.4 Matrix powers

### Motivation

Many engineering processes advance in steps: a digital filter updates its internal state once per sample, an iterative solver refines its estimate once per sweep, a Markov model of a communication channel moves between "good" and "bad" states once per packet, a power-flow solver repeats a Newton step. If each step is linear, $\mathbf{x}[k+1]=\mathbf{A}\mathbf{x}[k]$, then after $k$ steps

$$
\mathbf{x}[k]=\mathbf{A}^k\mathbf{x}[0].
$$

Everything about the long-run behaviour — growth, decay, oscillation, convergence to a steady state — is encoded in $\mathbf{A}^k$. Multiplying out $\mathbf{A}^k$ directly costs $k-1$ matrix products and reveals nothing about the trend. Diagonalization turns it into $n$ scalar powers and makes the trend obvious. The same machinery then gives the **matrix exponential** $e^{\mathbf{A}t}$, which plays for continuous-time systems $\dot{\mathbf{x}}=\mathbf{A}\mathbf{x}$ exactly the role that $\mathbf{A}^k$ plays for discrete-time ones.

### Powers by diagonalization

If $\mathbf{A}=\mathbf{P}\mathbf{D}\mathbf{P}^{-1}$, the inner factors cancel in pairs:

$$
\mathbf{A}^k=(\mathbf{P}\mathbf{D}\mathbf{P}^{-1})(\mathbf{P}\mathbf{D}\mathbf{P}^{-1})\cdots(\mathbf{P}\mathbf{D}\mathbf{P}^{-1})=\mathbf{P}\mathbf{D}^k\mathbf{P}^{-1},
\qquad
\mathbf{D}^k=\begin{bmatrix}\lambda_1^k&&\\&\ddots&\\&&\lambda_n^k\end{bmatrix}.
$$

The formula also holds for negative $k$ when $\mathbf{A}$ is invertible, and suggests a definition of $\mathbf{A}^{1/2}$ when all $\lambda_i\ge0$ (Section 10.6).

An equivalent and often more illuminating form expands the initial state in the eigenbasis. If $\mathbf{x}[0]=c_1\mathbf{v}_1+\dots+c_n\mathbf{v}_n$, i.e. $\mathbf{c}=\mathbf{P}^{-1}\mathbf{x}[0]$, then

$$
\mathbf{x}[k]=\mathbf{A}^k\mathbf{x}[0]=c_1\lambda_1^k\mathbf{v}_1+c_2\lambda_2^k\mathbf{v}_2+\dots+c_n\lambda_n^k\mathbf{v}_n .
$$

Each eigen-direction evolves independently, multiplied by $\lambda_i$ at every step. This **modal expansion** answers the qualitative questions at a glance:

- terms with $\lvert\lambda_i\rvert<1$ die out geometrically, those with $\lvert\lambda_i\rvert>1$ grow, those with $\lvert\lambda_i\rvert=1$ persist;
- a negative real $\lambda_i$ makes its term alternate in sign every step; a complex pair $re^{\pm j\theta}$ makes the pair rotate by $\theta$ per step while scaling by $r$;
- for large $k$, the term with the largest $\lvert\lambda_i\rvert$ (if unique, and if $c_i\ne0$) **dominates**: $\mathbf{x}[k]\approx c_1\lambda_1^k\mathbf{v}_1$, so the state lines up with the dominant eigenvector and grows or decays like $\rho(\mathbf{A})^k$.

The last observation is the **power method** for computing a dominant eigenpair: start from any vector, multiply repeatedly by $\mathbf{A}$ and normalise. It converges at the rate $\lvert\lambda_2/\lambda_1\rvert^k$ and is the basis of Google's PageRank and of many large-scale eigen-solvers ([Ch. 22](../part-7-numerical-computational/22-numerical-methods.md)).

### Powers of a Jordan block and the spectral-radius theorem

For a defective matrix we use the Jordan form instead. Write $\mathbf{J}_2(\lambda)=\lambda\mathbf{I}+\mathbf{N}$ with $\mathbf{N}=\begin{bmatrix}0&1\\0&0\end{bmatrix}$, $\mathbf{N}^2=\mathbf{0}$. Because $\lambda\mathbf{I}$ commutes with $\mathbf{N}$, the binomial theorem applies and stops after two terms:

$$
\mathbf{J}_2(\lambda)^k=\lambda^k\mathbf{I}+k\lambda^{k-1}\mathbf{N}=\begin{bmatrix}\lambda^k&k\lambda^{k-1}\\0&\lambda^k\end{bmatrix}.
$$

For an $m\times m$ block, entries $\binom{k}{i}\lambda^{k-i}$ appear for $i<m$. Polynomial factors like $k$ are eventually overwhelmed by $\lvert\lambda\rvert^k$ when $\lvert\lambda\rvert<1$, but they can make $\mathbf{A}^k$ grow substantially before it decays. Combining the diagonalizable and Jordan cases:

> [!theorem] Convergence of matrix powers
> For any square matrix $\mathbf{A}$:
> 1. $\mathbf{A}^k\to\mathbf{0}$ as $k\to\infty$ **if and only if** $\rho(\mathbf{A})<1$ (all eigenvalues strictly inside the unit circle).
> 2. If $\rho(\mathbf{A})>1$, $\lVert\mathbf{A}^k\mathbf{x}\rVert\to\infty$ for almost every starting vector $\mathbf{x}$.
> 3. If $\rho(\mathbf{A})=1$, the powers stay bounded if and only if every eigenvalue with $\lvert\lambda\rvert=1$ is non-defective.
> 4. (Gelfand) For any matrix norm, $\lVert\mathbf{A}^k\rVert^{1/k}\to\rho(\mathbf{A})$.

Part 1 is the stability condition for discrete-time linear systems (Section 10.9) and the convergence condition for stationary iterative solvers such as Jacobi and Gauss–Seidel, whose error obeys $\mathbf{e}[k+1]=\mathbf{G}\mathbf{e}[k]$ for an iteration matrix $\mathbf{G}$.

### The matrix exponential

For a scalar, $\dot x=ax$ has solution $x(t)=e^{at}x(0)$. For the system $\dot{\mathbf{x}}=\mathbf{A}\mathbf{x}$ we want a matrix that does the same job. Define, by analogy with the Taylor series of [Ch. 4](../part-2-calculus/04-differential-calculus.md),

$$
e^{\mathbf{A}t}=\mathbf{I}+\mathbf{A}t+\frac{(\mathbf{A}t)^2}{2!}+\frac{(\mathbf{A}t)^3}{3!}+\dots=\sum_{k=0}^\infty\frac{\mathbf{A}^kt^k}{k!}.
$$

The series converges for every square $\mathbf{A}$ and every $t$ (each entry is bounded by the scalar series for $e^{\lVert\mathbf{A}\rVert t}$). Differentiating term by term,

$$
\frac{d}{dt}e^{\mathbf{A}t}=\mathbf{A}+\mathbf{A}^2t+\frac{\mathbf{A}^3t^2}{2!}+\dots=\mathbf{A}e^{\mathbf{A}t}=e^{\mathbf{A}t}\mathbf{A},
$$

so $\mathbf{x}(t)=e^{\mathbf{A}t}\mathbf{x}(0)$ satisfies $\dot{\mathbf{x}}=\mathbf{A}\mathbf{x}$ with the right initial value; by uniqueness ([Ch. 11](../part-4-differential-equations/11-ordinary-differential-equations.md)) it is *the* solution. In control texts $e^{\mathbf{A}t}$ is called the **state-transition matrix** $\boldsymbol{\Phi}(t)$.

> [!theorem] Properties of the matrix exponential
> For square $\mathbf{A},\mathbf{B}$ and scalars $t,s$:
> 1. $e^{\mathbf{0}}=\mathbf{I}$ and $e^{\mathbf{A}(t+s)}=e^{\mathbf{A}t}e^{\mathbf{A}s}$;
> 2. $e^{\mathbf{A}t}$ is always invertible, with $(e^{\mathbf{A}t})^{-1}=e^{-\mathbf{A}t}$;
> 3. if $\mathbf{A}=\mathbf{P}\mathbf{D}\mathbf{P}^{-1}$ then $e^{\mathbf{A}t}=\mathbf{P}e^{\mathbf{D}t}\mathbf{P}^{-1}$ with $e^{\mathbf{D}t}=\operatorname{diag}(e^{\lambda_1t},\dots,e^{\lambda_nt})$;
> 4. the eigenvalues of $e^{\mathbf{A}t}$ are $e^{\lambda_it}$, so $\det e^{\mathbf{A}t}=e^{(\operatorname{tr}\mathbf{A})t}$;
> 5. $e^{(\mathbf{A}+\mathbf{B})t}=e^{\mathbf{A}t}e^{\mathbf{B}t}$ **if** $\mathbf{A}\mathbf{B}=\mathbf{B}\mathbf{A}$ — but not in general.

Property 3 follows from $\mathbf{A}^k=\mathbf{P}\mathbf{D}^k\mathbf{P}^{-1}$ applied term by term. With the modal expansion $\mathbf{x}(0)=\sum c_i\mathbf{v}_i$ it gives the continuous-time analogue of the discrete result:

$$
\mathbf{x}(t)=c_1e^{\lambda_1t}\mathbf{v}_1+c_2e^{\lambda_2t}\mathbf{v}_2+\dots+c_ne^{\lambda_nt}\mathbf{v}_n .
$$

Now the dividing line is $\operatorname{Re}\lambda_i=0$ instead of $\lvert\lambda_i\rvert=1$: $e^{\lambda t}$ decays if $\operatorname{Re}\lambda<0$. The two pictures are linked by sampling: if $\mathbf{x}(t)$ is sampled every $T$ seconds, $\mathbf{x}[k+1]=e^{\mathbf{A}T}\mathbf{x}[k]$, and the eigenvalue $\lambda$ maps to $z=e^{\lambda T}$ — the left half of the $s$-plane maps to the inside of the unit circle of the $z$-plane, exactly the $s\to z$ correspondence of [Ch. 17](../part-5-transforms/17-z-transforms.md).

For a Jordan block, $e^{\mathbf{J}_2(\lambda)t}=e^{\lambda t}e^{\mathbf{N}t}=e^{\lambda t}(\mathbf{I}+\mathbf{N}t)=e^{\lambda t}\begin{bmatrix}1&t\\0&1\end{bmatrix}$ — the origin of the $te^{\lambda t}$ terms. A rotation generator gives the familiar oscillation:

$$
\exp\!\left(\begin{bmatrix}0&-\omega\\\omega&0\end{bmatrix}t\right)=\begin{bmatrix}\cos\omega t&-\sin\omega t\\\sin\omega t&\cos\omega t\end{bmatrix},
$$

the real-matrix form of $e^{j\omega t}$ — a rotating phasor.

### Functions of a matrix via Cayley–Hamilton (Sylvester interpolation)

By Cayley–Hamilton, every power $\mathbf{A}^k$ with $k\ge n$ reduces to a combination of $\mathbf{I},\mathbf{A},\dots,\mathbf{A}^{n-1}$, and so does any convergent power series. Hence for an analytic function $f$ (such as $\lambda^k$ or $e^{\lambda t}$),

$$
f(\mathbf{A})=\alpha_0\mathbf{I}+\alpha_1\mathbf{A}+\dots+\alpha_{n-1}\mathbf{A}^{n-1}
$$

for some scalars $\alpha_i$. Applying both sides to an eigenvector shows that the polynomial $r(\lambda)=\alpha_0+\alpha_1\lambda+\dots+\alpha_{n-1}\lambda^{n-1}$ must agree with $f$ at each eigenvalue:

$$
f(\lambda_i)=r(\lambda_i),\qquad i=1,\dots,n.
$$

With distinct eigenvalues this is $n$ linear equations for the $n$ unknowns $\alpha_i$. For an eigenvalue of multiplicity $m$, also match the first $m-1$ derivatives, $f^{(q)}(\lambda_i)=r^{(q)}(\lambda_i)$. This method avoids eigenvectors altogether and is convenient for $2\times2$ and $3\times3$ hand calculations.

### Worked examples

> [!example] Worked Example 10.13 — A tenth power
> Compute $\mathbf{A}^{10}$ for $\mathbf{A}=\begin{bmatrix}4&1\\2&3\end{bmatrix}$.
>
> **Step 1 — use Example 10.9.** $\mathbf{P}=\begin{bmatrix}1&1\\1&-2\end{bmatrix}$, $\mathbf{D}=\operatorname{diag}(5,2)$, $\mathbf{P}^{-1}=\frac13\begin{bmatrix}2&1\\1&-1\end{bmatrix}$.
>
> **Step 2 — multiply.**
> $$
> \mathbf{A}^{10}=\begin{bmatrix}1&1\\1&-2\end{bmatrix}\begin{bmatrix}5^{10}&0\\0&2^{10}\end{bmatrix}\frac13\begin{bmatrix}2&1\\1&-1\end{bmatrix}
> =\frac13\begin{bmatrix}2\cdot5^{10}+2^{10}&5^{10}-2^{10}\\2\cdot5^{10}-2\cdot2^{10}&5^{10}+2\cdot2^{10}\end{bmatrix}.
> $$
> **Step 3 — numbers.** $5^{10}=9\,765\,625$, $2^{10}=1024$:
> $$
> \mathbf{A}^{10}=\begin{bmatrix}6\,510\,758&3\,254\,867\\6\,509\,734&3\,255\,891\end{bmatrix}.
> $$
> **Checks.** $\operatorname{tr}\mathbf{A}^{10}=9\,766\,649=5^{10}+2^{10}$ $\checkmark$. The general formula at $k=1$ returns $\frac13\begin{bmatrix}12&3\\6&9\end{bmatrix}=\mathbf{A}$ $\checkmark$. Notice how the $5^{10}$ term dominates: the columns of $\mathbf{A}^{10}$ are both nearly parallel to $\mathbf{v}_1=(1,1)$ — the power method at work.
>
> **Sylvester alternative.** Write $\mathbf{A}^{10}=\alpha_0\mathbf{I}+\alpha_1\mathbf{A}$ with $5^{10}=\alpha_0+5\alpha_1$, $2^{10}=\alpha_0+2\alpha_1$. Then $\alpha_1=(5^{10}-2^{10})/3=3\,254\,867$ and $\alpha_0=2^{10}-2\alpha_1=-6\,508\,710$. The $(1,1)$ entry is $\alpha_0+4\alpha_1=6\,510\,758$ $\checkmark$.

> [!example] Worked Example 10.14 — Fibonacci numbers in closed form
> The recursion $F_{k+1}=F_k+F_{k-1}$, $F_0=0$, $F_1=1$, is the prototype of a second-order linear difference equation (the same structure as an IIR filter with no input). With $\mathbf{x}[k]=(F_{k+1},F_k)$,
> $$
> \mathbf{x}[k+1]=\begin{bmatrix}1&1\\1&0\end{bmatrix}\mathbf{x}[k]=\mathbf{F}\mathbf{x}[k],\qquad \mathbf{x}[0]=\begin{bmatrix}1\\0\end{bmatrix}.
> $$
> **Step 1 — eigenvalues.** $\lambda^2-\lambda-1=0$: $\varphi=\frac{1+\sqrt5}2\approx1.618$ and $\psi=\frac{1-\sqrt5}2\approx-0.618$.
>
> **Step 2 — eigenvectors.** $\mathbf{F}-\lambda\mathbf{I}=\begin{bmatrix}1-\lambda&1\\1&-\lambda\end{bmatrix}$; the second row gives $\mathbf{v}=(\lambda,1)$. So $\mathbf{v}_\varphi=(\varphi,1)$, $\mathbf{v}_\psi=(\psi,1)$.
>
> **Step 3 — expand the initial state.** $(1,0)=c_1(\varphi,1)+c_2(\psi,1)$ gives $c_2=-c_1$ and $c_1(\varphi-\psi)=1$, so $c_1=1/\sqrt5$.
>
> **Step 4 — read off the second component.**
> $$
> F_k=c_1\varphi^k+c_2\psi^k=\frac{\varphi^k-\psi^k}{\sqrt5}\qquad\text{(Binet's formula).}
> $$
> **Check.** $F_{10}=(1.618034^{10}-(-0.618034)^{10})/\sqrt5=(122.9919-0.0081)/2.23607=55.000$ $\checkmark$, and directly $\mathbf{F}^{10}=\begin{bmatrix}89&55\\55&34\end{bmatrix}$. Because $\lvert\psi\rvert<1$, $F_k$ is the integer nearest to $\varphi^k/\sqrt5$, and $F_{k+1}/F_k\to\varphi$: the dominant eigenvalue sets the growth rate.

> [!example] Worked Example 10.15 — A two-state Markov channel and its steady state
> A burst-error channel (Gilbert–Elliott model) is either Good (G) or Bad (B) during each packet. From G it stays G with probability $0.9$; from B it returns to G with probability $0.4$. With $\mathbf{x}[k]=(P(\text{G}),P(\text{B}))$ at packet $k$,
> $$
> \mathbf{x}[k+1]=\mathbf{M}\mathbf{x}[k],\qquad \mathbf{M}=\begin{bmatrix}0.9&0.4\\0.1&0.6\end{bmatrix}
> $$
> (each column sums to $1$). Find $\mathbf{x}[k]$ if the channel starts Bad, and the long-run fraction of Good packets.
>
> **Step 1 — eigenvalues.** $\operatorname{tr}\mathbf{M}=1.5$, $\det\mathbf{M}=0.54-0.04=0.5$; $\lambda^2-1.5\lambda+0.5=(\lambda-1)(\lambda-0.5)$. A column-stochastic matrix always has $\lambda=1$, because $(1,1)\mathbf{M}=(1,1)$: the all-ones row is a left eigenvector.
>
> **Step 2 — eigenvectors.** $\lambda=1$: $\mathbf{M}-\mathbf{I}=\begin{bmatrix}-0.1&0.4\\0.1&-0.4\end{bmatrix}$ gives $\mathbf{v}_1\propto(4,1)$, normalised as a probability vector $\boldsymbol{\pi}=(0.8,0.2)$. $\lambda=0.5$: $\mathbf{M}-0.5\mathbf{I}=\begin{bmatrix}0.4&0.4\\0.1&0.1\end{bmatrix}$ gives $\mathbf{v}_2=(1,-1)$ (entries sum to zero, as they must for a mode that does not change total probability).
>
> **Step 3 — expand $\mathbf{x}[0]=(0,1)$.** $(0,1)=\boldsymbol{\pi}+c_2(1,-1)$ forces $c_2=-0.8$.
>
> **Step 4 — solution.**
> $$
> \mathbf{x}[k]=\begin{bmatrix}0.8\\0.2\end{bmatrix}-0.8\,(0.5)^k\begin{bmatrix}1\\-1\end{bmatrix},\qquad P_k(\text{G})=0.8\big(1-0.5^k\big).
> $$
> **Check.** $k=1$: $P(\text{G})=0.4$, matching one step from B $\checkmark$. $k=10$: $0.8(1-1/1024)=0.79922$ (numerical $\mathbf{M}^{10}(0,1)^\mathsf{T}$ agrees). Long run: **80 % of packets are Good**, whatever the initial state. The second eigenvalue $0.5$ sets the **mixing rate**: the memory of the initial state halves every packet. A channel with $\lambda_2$ close to $1$ has long error bursts — which is why interleaver depth is chosen from the second eigenvalue.

> [!example] Worked Example 10.16 — The state-transition matrix $e^{\mathbf{A}t}$
> Compute $e^{\mathbf{A}t}$ for $\mathbf{A}=\begin{bmatrix}0&1\\-2&-3\end{bmatrix}$ (the companion matrix of $\ddot y+3\dot y+2y=0$) and find $\mathbf{x}(t)$ for $\mathbf{x}(0)=(1,0)$.
>
> **Method 1 — diagonalization.** $\lambda^2+3\lambda+2=(\lambda+1)(\lambda+2)$. Eigenvectors of a companion matrix are $(1,\lambda)$: $\mathbf{v}_1=(1,-1)$, $\mathbf{v}_2=(1,-2)$. With $\mathbf{P}=\begin{bmatrix}1&1\\-1&-2\end{bmatrix}$, $\det\mathbf{P}=-1$, $\mathbf{P}^{-1}=\begin{bmatrix}2&1\\-1&-1\end{bmatrix}$:
> $$
> e^{\mathbf{A}t}=\mathbf{P}\begin{bmatrix}e^{-t}&0\\0&e^{-2t}\end{bmatrix}\mathbf{P}^{-1}
> =\begin{bmatrix}2e^{-t}-e^{-2t}&e^{-t}-e^{-2t}\\-2e^{-t}+2e^{-2t}&-e^{-t}+2e^{-2t}\end{bmatrix}.
> $$
> **Method 2 — Sylvester.** $e^{\mathbf{A}t}=\alpha_0\mathbf{I}+\alpha_1\mathbf{A}$ with $e^{-t}=\alpha_0-\alpha_1$ and $e^{-2t}=\alpha_0-2\alpha_1$, so $\alpha_1=e^{-t}-e^{-2t}$, $\alpha_0=2e^{-t}-e^{-2t}$. Then $e^{\mathbf{A}t}=\begin{bmatrix}\alpha_0&\alpha_1\\-2\alpha_1&\alpha_0-3\alpha_1\end{bmatrix}$, identical to Method 1 (e.g. $\alpha_0-3\alpha_1=-e^{-t}+2e^{-2t}$ $\checkmark$).
>
> **Response.** $\mathbf{x}(t)=e^{\mathbf{A}t}(1,0)^\mathsf{T}$ is the first column: $y(t)=2e^{-t}-e^{-2t}$, $\dot y(t)=-2e^{-t}+2e^{-2t}$.
>
> **Checks.** $t=0$ gives $\mathbf{I}$ $\checkmark$. $\det e^{\mathbf{A}t}=e^{-3t}=e^{(\operatorname{tr}\mathbf{A})t}$ (expand: $(2e^{-t}-e^{-2t})(-e^{-t}+2e^{-2t})+(e^{-t}-e^{-2t})(2e^{-t}-2e^{-2t})=e^{-3t}$) $\checkmark$. At $t=1$: $2e^{-1}-e^{-2}=0.6004$ and $e^{-1}-e^{-2}=0.2325$, matching a numerical `expm` $\checkmark$.

> [!example] Worked Example 10.17 — Exponential of a defective matrix
> For the defective $\mathbf{A}=\begin{bmatrix}3&1\\-1&1\end{bmatrix}$ of Example 10.11 (double eigenvalue $2$), Sylvester's method needs the derivative condition. Write $e^{\mathbf{A}t}=\alpha_0\mathbf{I}+\alpha_1\mathbf{A}$ and match $f(\lambda)=e^{\lambda t}$ and $f'(\lambda)=te^{\lambda t}$ (derivative with respect to $\lambda$) at $\lambda=2$:
> $$
> e^{2t}=\alpha_0+2\alpha_1,\qquad te^{2t}=\alpha_1 .
> $$
> So $\alpha_1=te^{2t}$, $\alpha_0=(1-2t)e^{2t}$, and
> $$
> e^{\mathbf{A}t}=(1-2t)e^{2t}\mathbf{I}+te^{2t}\begin{bmatrix}3&1\\-1&1\end{bmatrix}=e^{2t}\begin{bmatrix}1+t&t\\-t&1-t\end{bmatrix}.
> $$
> **Check** via Jordan form: $\mathbf{P}e^{\mathbf{J}t}\mathbf{P}^{-1}$ with $e^{\mathbf{J}t}=e^{2t}\begin{bmatrix}1&t\\0&1\end{bmatrix}$ gives the same matrix, and $\frac{d}{dt}$ at $t=0$ returns $\begin{bmatrix}2+1&1\\-1&2-1\end{bmatrix}=\mathbf{A}$ $\checkmark$.

The graph shows the discrete-time counterpart of a complex eigenvalue pair. For $\mathbf{A}=0.9\begin{bmatrix}\cos30^\circ&-\sin30^\circ\\\sin30^\circ&\cos30^\circ\end{bmatrix}$ (eigenvalues $0.9e^{\pm j\pi/6}$), the iterates $\mathbf{x}[k]=\mathbf{A}^k(1,0)$ rotate $30^\circ$ per step and shrink by $10\%$; they lie on the logarithmic spiral $r=0.9^{6\theta/\pi}=e^{-0.2012\theta}$.

```desmos-graph
left=-1; right=1.2; bottom=-0.6; top=1
---
r=e^{-0.2012\theta}|blue|dashed|0<\theta<12.57
(1,0)|label:x[0]
(0.7794,0.45)|label:x[1]
(0.405,0.7015)|label:x[2]
(0,0.729)|label:x[3]
(-0.328,0.5682)|label:x[4]
(-0.5114,0.2952)|label:x[5]
(-0.5314,0)|label:x[6]
(-0.4142,-0.2391)|label:x[7]
(-0.2152,-0.3728)|label:x[8]
(0,-0.3874)|label:x[9]
(0.1743,-0.302)|label:x[10]
(0.2824,0)|label:x[12]
x^2+y^2=1|grey|dotted
```
*Iterates of a discrete-time system with eigenvalues $0.9e^{\pm j\pi/6}$: each step rotates by $\angle\lambda=30^\circ$ and scales by $\lvert\lambda\rvert=0.9$, so the state spirals into the origin along $r=e^{-0.2012\theta}$; the dotted unit circle is the stability boundary $\lvert\lambda\rvert=1$.*

> [!warning] Common Pitfall — exponentiating entry by entry
> $e^{\mathbf{A}t}$ is **not** the matrix of $e^{a_{ij}t}$ (except for diagonal $\mathbf{A}$), and $\mathbf{A}^k$ is not the matrix of $a_{ij}^k$. Equally, $e^{(\mathbf{A}+\mathbf{B})t}\ne e^{\mathbf{A}t}e^{\mathbf{B}t}$ unless $\mathbf{A}$ and $\mathbf{B}$ commute — which matters when a state matrix is split into "nominal" and "perturbation" parts, or when switching between circuit topologies in a switched-mode converter.

> [!tip] Engineering Note — computing $e^{\mathbf{A}T}$ in practice
> Discretising a continuous-time model for a digital controller or a circuit simulator needs $\boldsymbol{\Phi}=e^{\mathbf{A}T}$ for the sample period $T$. Diagonalization is fine for small, well-conditioned problems, but robust software (MATLAB `expm`, SciPy `scipy.linalg.expm`) uses **scaling and squaring** with a Padé approximant: compute $e^{\mathbf{A}T/2^s}$ for an $s$ that makes the argument small, then square $s$ times. This avoids both the ill-conditioning of nearly defective $\mathbf{P}$ and the catastrophic cancellation of the raw Taylor series for large $\lVert\mathbf{A}T\rVert$ (Moler and Van Loan's classic survey lists nineteen methods and why most of them are dubious).

## 10.5 Symmetric matrices

### Motivation

A remarkable fraction of the matrices an engineer meets are symmetric, $\mathbf{A}^\mathsf{T}=\mathbf{A}$, or their complex cousins, Hermitian, $\mathbf{A}^\mathsf{H}=\mathbf{A}$. The nodal conductance matrix of a resistive network is symmetric (reciprocity: the current at node $j$ due to a voltage at node $i$ equals the current at $i$ due to the same voltage at $j$). Inductance and capacitance matrices of coupled conductors are symmetric because they come from stored energy. Covariance and correlation matrices of noise are symmetric (real signals) or Hermitian (complex baseband). Hessian matrices of second derivatives are symmetric ([Ch. 6](../part-2-calculus/06-multivariable-calculus.md)). The impedance matrix of a lossless reciprocal network at a given frequency is $j$ times a real symmetric matrix.

For all of these, eigen-theory is as good as it can possibly be: the eigenvalues are real, the eigenvectors are orthogonal, and diagonalization never fails. This section proves these facts and draws out their geometry.

### Real eigenvalues and orthogonal eigenvectors

> [!theorem] Symmetric and Hermitian matrices: real eigenvalues
> If $\mathbf{A}$ is Hermitian ($\mathbf{A}^\mathsf{H}=\mathbf{A}$; in particular if $\mathbf{A}$ is real symmetric), every eigenvalue of $\mathbf{A}$ is real.

*Proof.* Let $\mathbf{A}\mathbf{v}=\lambda\mathbf{v}$ with $\mathbf{v}\ne\mathbf{0}$, possibly complex. Multiply on the left by $\mathbf{v}^\mathsf{H}$: $\mathbf{v}^\mathsf{H}\mathbf{A}\mathbf{v}=\lambda\,\mathbf{v}^\mathsf{H}\mathbf{v}=\lambda\lVert\mathbf{v}\rVert^2$. The scalar $\mathbf{v}^\mathsf{H}\mathbf{A}\mathbf{v}$ is real, because it equals its own conjugate: $(\mathbf{v}^\mathsf{H}\mathbf{A}\mathbf{v})^*=(\mathbf{v}^\mathsf{H}\mathbf{A}\mathbf{v})^\mathsf{H}=\mathbf{v}^\mathsf{H}\mathbf{A}^\mathsf{H}\mathbf{v}=\mathbf{v}^\mathsf{H}\mathbf{A}\mathbf{v}$. Since $\lVert\mathbf{v}\rVert^2>0$ is real too, $\lambda=\mathbf{v}^\mathsf{H}\mathbf{A}\mathbf{v}/\lVert\mathbf{v}\rVert^2$ is real. $\square$

For a real symmetric matrix, real $\lambda$ means $\mathbf{A}-\lambda\mathbf{I}$ is real, so the eigenvectors can be chosen real.

> [!theorem] Orthogonality of eigenvectors
> Eigenvectors of a Hermitian matrix belonging to distinct eigenvalues are orthogonal.

*Proof.* Let $\mathbf{A}\mathbf{v}_1=\lambda_1\mathbf{v}_1$, $\mathbf{A}\mathbf{v}_2=\lambda_2\mathbf{v}_2$, $\lambda_1\ne\lambda_2$ (both real). Then
$$
\lambda_1\,\mathbf{v}_2^\mathsf{H}\mathbf{v}_1=\mathbf{v}_2^\mathsf{H}\mathbf{A}\mathbf{v}_1=(\mathbf{A}^\mathsf{H}\mathbf{v}_2)^\mathsf{H}\mathbf{v}_1=(\mathbf{A}\mathbf{v}_2)^\mathsf{H}\mathbf{v}_1=\lambda_2\,\mathbf{v}_2^\mathsf{H}\mathbf{v}_1 ,
$$
so $(\lambda_1-\lambda_2)\langle\mathbf{v}_1,\mathbf{v}_2\rangle=0$ and $\langle\mathbf{v}_1,\mathbf{v}_2\rangle=0$. $\square$

Within a repeated eigenvalue's eigenspace the eigenvectors are not automatically orthogonal, but the eigenspace is a subspace, so Gram–Schmidt ([Ch. 9](09-vector-spaces-and-linear-transformations.md)) produces an orthonormal basis of it — and the orthogonalised vectors are still eigenvectors.

### The spectral theorem

> [!theorem] Spectral theorem
> A real $n\times n$ matrix $\mathbf{A}$ is symmetric **if and only if** it can be factored as
> $$
> \mathbf{A}=\mathbf{Q}\boldsymbol{\Lambda}\mathbf{Q}^\mathsf{T},
> $$
> where $\boldsymbol{\Lambda}=\operatorname{diag}(\lambda_1,\dots,\lambda_n)$ is real and $\mathbf{Q}$ is **orthogonal** ($\mathbf{Q}^\mathsf{T}\mathbf{Q}=\mathbf{I}$, so $\mathbf{Q}^{-1}=\mathbf{Q}^\mathsf{T}$) with orthonormal eigenvectors $\mathbf{q}_1,\dots,\mathbf{q}_n$ as columns. Likewise a complex matrix is Hermitian if and only if $\mathbf{A}=\mathbf{U}\boldsymbol{\Lambda}\mathbf{U}^\mathsf{H}$ with $\mathbf{U}$ unitary ($\mathbf{U}^\mathsf{H}\mathbf{U}=\mathbf{I}$) and $\boldsymbol{\Lambda}$ real.

*Proof sketch.* "If": $(\mathbf{Q}\boldsymbol{\Lambda}\mathbf{Q}^\mathsf{T})^\mathsf{T}=\mathbf{Q}\boldsymbol{\Lambda}^\mathsf{T}\mathbf{Q}^\mathsf{T}=\mathbf{Q}\boldsymbol{\Lambda}\mathbf{Q}^\mathsf{T}$. "Only if", by induction on $n$: $\mathbf{A}$ has a real eigenvalue $\lambda_1$ with a unit eigenvector $\mathbf{q}_1$. Extend $\mathbf{q}_1$ to an orthonormal basis and let $\mathbf{Q}_1$ have these as columns. The first column of $\mathbf{Q}_1^\mathsf{T}\mathbf{A}\mathbf{Q}_1$ is $\lambda_1\mathbf{e}_1$; because $\mathbf{Q}_1^\mathsf{T}\mathbf{A}\mathbf{Q}_1$ is again symmetric, its first row is $\lambda_1\mathbf{e}_1^\mathsf{T}$ too:
$$
\mathbf{Q}_1^\mathsf{T}\mathbf{A}\mathbf{Q}_1=\begin{bmatrix}\lambda_1&\mathbf{0}^\mathsf{T}\\\mathbf{0}&\mathbf{A}_1\end{bmatrix},
$$
with $\mathbf{A}_1$ symmetric of size $n-1$. Apply the induction hypothesis to $\mathbf{A}_1$ and assemble. Note what symmetry bought: the zero block in the *top right*, which is what prevents Jordan blocks. $\square$

So a symmetric matrix is never defective, and its eigenvector matrix is as well-conditioned as a matrix can be: $\kappa(\mathbf{Q})=1$.

### Spectral decomposition

Multiplying out $\mathbf{Q}\boldsymbol{\Lambda}\mathbf{Q}^\mathsf{T}$ column-by-row gives the **spectral decomposition**

$$
\mathbf{A}=\lambda_1\mathbf{q}_1\mathbf{q}_1^\mathsf{T}+\lambda_2\mathbf{q}_2\mathbf{q}_2^\mathsf{T}+\dots+\lambda_n\mathbf{q}_n\mathbf{q}_n^\mathsf{T}.
$$

Each $\mathbf{P}_i=\mathbf{q}_i\mathbf{q}_i^\mathsf{T}$ is the orthogonal projector onto the line through $\mathbf{q}_i$ ([Ch. 9](09-vector-spaces-and-linear-transformations.md)); these projectors satisfy $\mathbf{P}_i^2=\mathbf{P}_i$, $\mathbf{P}_i\mathbf{P}_k=\mathbf{0}$ for $i\ne k$ and $\sum_i\mathbf{P}_i=\mathbf{I}$. A symmetric matrix is therefore a weighted sum of mutually perpendicular projections: to apply $\mathbf{A}$, resolve $\mathbf{x}$ into its components along the eigenvectors, scale each by its eigenvalue, and add. Functions of $\mathbf{A}$ follow at once: $f(\mathbf{A})=\sum_if(\lambda_i)\mathbf{q}_i\mathbf{q}_i^\mathsf{T}$, e.g. $\mathbf{A}^{-1}=\sum\lambda_i^{-1}\mathbf{q}_i\mathbf{q}_i^\mathsf{T}$.

### Quadratic forms and principal axes

A **quadratic form** is a function $Q(\mathbf{x})=\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}=\sum_{i,k}a_{ik}x_ix_k$. Only the symmetric part of $\mathbf{A}$ contributes (since $\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}=\mathbf{x}^\mathsf{T}\mathbf{A}^\mathsf{T}\mathbf{x}$), so we always take $\mathbf{A}$ symmetric: for $Q=ax^2+2bxy+cy^2$, $\mathbf{A}=\begin{bmatrix}a&b\\b&c\end{bmatrix}$ — note the **half** of the cross coefficient off the diagonal. Energy expressions are quadratic forms: $\tfrac12\mathbf{i}^\mathsf{T}\mathbf{L}\mathbf{i}$ for coupled inductors, $\tfrac12\mathbf{v}^\mathsf{T}\mathbf{C}\mathbf{v}$ for a capacitance matrix, $\mathbf{v}^\mathsf{T}\mathbf{G}\mathbf{v}$ for power dissipated in a conductance network.

Substituting $\mathbf{x}=\mathbf{Q}\mathbf{y}$ (a rotation, possibly with a reflection) and using $\mathbf{Q}^\mathsf{T}\mathbf{A}\mathbf{Q}=\boldsymbol{\Lambda}$:

$$
\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}=\mathbf{y}^\mathsf{T}\boldsymbol{\Lambda}\mathbf{y}=\lambda_1y_1^2+\lambda_2y_2^2+\dots+\lambda_ny_n^2 .
$$

> [!theorem] Principal axes theorem
> Every real quadratic form can be reduced by an orthogonal change of variables $\mathbf{x}=\mathbf{Q}\mathbf{y}$ to a sum of squares $\sum\lambda_iy_i^2$ with no cross terms. The new axes — the eigenvectors $\mathbf{q}_i$ — are the **principal axes** of the form. For $n=2$, the level curve $\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}=1$ is an ellipse with semi-axes $1/\sqrt{\lambda_1}$, $1/\sqrt{\lambda_2}$ along $\mathbf{q}_1,\mathbf{q}_2$ when both eigenvalues are positive, and a hyperbola when they have opposite signs.

### The Rayleigh quotient and extremal eigenvalues

> [!definition] Rayleigh quotient
> For symmetric $\mathbf{A}$ and $\mathbf{x}\ne\mathbf{0}$, the **Rayleigh quotient** is
> $$
> R(\mathbf{x})=\frac{\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}}{\mathbf{x}^\mathsf{T}\mathbf{x}} .
> $$

Expand $\mathbf{x}=\sum y_i\mathbf{q}_i$ with eigenvalues ordered $\lambda_1\ge\lambda_2\ge\dots\ge\lambda_n$:

$$
R(\mathbf{x})=\frac{\lambda_1y_1^2+\dots+\lambda_ny_n^2}{y_1^2+\dots+y_n^2},
$$

a weighted average of the eigenvalues with weights $y_i^2\ge0$. Hence:

> [!theorem] Rayleigh's principle
> $$
> \lambda_{\min}=\lambda_n\le R(\mathbf{x})\le\lambda_1=\lambda_{\max}\quad\text{for all }\mathbf{x}\ne\mathbf{0},
> $$
> with $\max R=\lambda_1$ attained at $\mathbf{x}=\mathbf{q}_1$ and $\min R=\lambda_n$ at $\mathbf{x}=\mathbf{q}_n$. More generally (Courant–Fischer), $\lambda_k$ is the maximum of $R$ over vectors orthogonal to $\mathbf{q}_1,\dots,\mathbf{q}_{k-1}$.

Two consequences are used constantly. First, **$R$ is stationary at every eigenvector**, so an approximate eigenvector with error $\varepsilon$ gives an eigenvalue estimate with error of order $\varepsilon^2$ — the reason Rayleigh's method gives excellent natural-frequency estimates from rough guesses of mode shape. Second, maximising a quadratic form on the unit sphere is an eigenvalue problem: the direction of maximum output power of a beamformer, the direction of largest variance in data (principal component analysis), the polarisation that maximises received power — each is the top eigenvector of a suitable symmetric or Hermitian matrix.

### Worked examples

> [!example] Worked Example 10.18 — Orthogonal diagonalization of a $2\times2$ symmetric matrix
> Diagonalize $\mathbf{S}=\begin{bmatrix}3&1\\1&3\end{bmatrix}$ orthogonally and write its spectral decomposition.
>
> **Step 1 — eigenvalues.** $\lambda^2-6\lambda+8=0$: $\lambda_1=4$, $\lambda_2=2$.
>
> **Step 2 — unit eigenvectors.** $\lambda=4$: $\mathbf{S}-4\mathbf{I}=\begin{bmatrix}-1&1\\1&-1\end{bmatrix}$, $\mathbf{q}_1=\frac1{\sqrt2}(1,1)$. $\lambda=2$: $\mathbf{S}-2\mathbf{I}=\begin{bmatrix}1&1\\1&1\end{bmatrix}$, $\mathbf{q}_2=\frac1{\sqrt2}(1,-1)$. They are orthogonal, as guaranteed.
>
> **Step 3 — factor.**
> $$
> \mathbf{Q}=\frac1{\sqrt2}\begin{bmatrix}1&1\\1&-1\end{bmatrix},\quad\boldsymbol{\Lambda}=\begin{bmatrix}4&0\\0&2\end{bmatrix},\quad\mathbf{S}=\mathbf{Q}\boldsymbol{\Lambda}\mathbf{Q}^\mathsf{T}.
> $$
> No matrix inversion was needed: $\mathbf{Q}^{-1}=\mathbf{Q}^\mathsf{T}$.
>
> **Step 4 — spectral decomposition.**
> $$
> \mathbf{S}=4\cdot\frac12\begin{bmatrix}1&1\\1&1\end{bmatrix}+2\cdot\frac12\begin{bmatrix}1&-1\\-1&1\end{bmatrix}=\begin{bmatrix}2&2\\2&2\end{bmatrix}+\begin{bmatrix}1&-1\\-1&1\end{bmatrix}=\begin{bmatrix}3&1\\1&3\end{bmatrix}\ \checkmark
> $$
> **Interpretation.** This is the matrix of a symmetric two-port or a pair of identical coupled lines. The eigenvectors are the **common mode** $(1,1)$ and **differential mode** $(1,-1)$ of [Ch. 9](09-vector-spaces-and-linear-transformations.md); any symmetric structure with two identical halves has them as eigenvectors, which is why even/odd-mode analysis works.

> [!example] Worked Example 10.19 — A repeated eigenvalue and Gram–Schmidt
> Orthogonally diagonalize $\mathbf{A}=\begin{bmatrix}2&1&1\\1&2&1\\1&1&2\end{bmatrix}$.
>
> **Step 1 — eigenvalues.** Write $\mathbf{A}=\mathbf{I}+\mathbf{J}$ with $\mathbf{J}$ the all-ones matrix. $\mathbf{J}$ has rank $1$, so eigenvalue $0$ twice, and $\operatorname{tr}\mathbf{J}=3$ is its third eigenvalue (eigenvector $(1,1,1)$). Adding $1$: the eigenvalues of $\mathbf{A}$ are $4,1,1$. (Check: $\operatorname{tr}\mathbf{A}=6$ $\checkmark$; $\det\mathbf{A}=2(3)-1(1)+1(-1)=4$ $\checkmark$.)
>
> **Step 2 — $\lambda=4$.** $\mathbf{q}_1=\frac1{\sqrt3}(1,1,1)$.
>
> **Step 3 — $\lambda=1$.** $\mathbf{A}-\mathbf{I}=\mathbf{J}$, so the eigenspace is the plane $x_1+x_2+x_3=0$. A natural basis is $\mathbf{w}_1=(1,-1,0)$, $\mathbf{w}_2=(1,0,-1)$ — valid eigenvectors but **not orthogonal** ($\mathbf{w}_1\cdot\mathbf{w}_2=1$). Apply Gram–Schmidt inside the eigenspace:
> $$
> \mathbf{q}_2=\frac{\mathbf{w}_1}{\lVert\mathbf{w}_1\rVert}=\frac1{\sqrt2}(1,-1,0),\qquad
> \mathbf{w}_2-(\mathbf{w}_2\cdot\mathbf{q}_2)\mathbf{q}_2=(1,0,-1)-\tfrac12(1,-1,0)=\left(\tfrac12,\tfrac12,-1\right),
> $$
> normalised: $\mathbf{q}_3=\frac1{\sqrt6}(1,1,-2)$. It still satisfies $x_1+x_2+x_3=0$, so it is still an eigenvector.
>
> **Step 4 — result.**
> $$
> \mathbf{Q}=\begin{bmatrix}1/\sqrt3&1/\sqrt2&1/\sqrt6\\1/\sqrt3&-1/\sqrt2&1/\sqrt6\\1/\sqrt3&0&-2/\sqrt6\end{bmatrix},\qquad\boldsymbol{\Lambda}=\operatorname{diag}(4,1,1).
> $$
> **Check via projectors.** $\mathbf{q}_1\mathbf{q}_1^\mathsf{T}=\frac13\mathbf{J}$ and, since the three projectors sum to $\mathbf{I}$, $\mathbf{q}_2\mathbf{q}_2^\mathsf{T}+\mathbf{q}_3\mathbf{q}_3^\mathsf{T}=\mathbf{I}-\frac13\mathbf{J}$. So $\mathbf{Q}\boldsymbol{\Lambda}\mathbf{Q}^\mathsf{T}=4\cdot\frac13\mathbf{J}+1\cdot(\mathbf{I}-\frac13\mathbf{J})=\mathbf{I}+\mathbf{J}=\mathbf{A}$ $\checkmark$. In a balanced three-phase system this same matrix structure gives the **zero-sequence** vector $(1,1,1)$ and a two-dimensional space of positive/negative-sequence vectors — the symmetrical-components decomposition is a spectral decomposition.
>
> **Rayleigh check.** $\mathbf{x}=(1,1,0)$: $\mathbf{A}\mathbf{x}=(3,3,2)$, $R=6/2=3$, which lies in $[\lambda_{\min},\lambda_{\max}]=[1,4]$ $\checkmark$.

> [!example] Worked Example 10.20 — Principal axes of a conic
> Identify and sketch the curve $5x^2+4xy+2y^2=6$.
>
> **Step 1 — matrix.** $Q=\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}$ with $\mathbf{A}=\begin{bmatrix}5&2\\2&2\end{bmatrix}$ (half of $4$ off the diagonal).
>
> **Step 2 — eigen-decomposition.** $\lambda^2-7\lambda+6=0$: $\lambda_1=6$, $\lambda_2=1$. $\lambda=6$: $\mathbf{A}-6\mathbf{I}=\begin{bmatrix}-1&2\\2&-4\end{bmatrix}$, $\mathbf{q}_1=\frac1{\sqrt5}(2,1)$. $\lambda=1$: $\mathbf{q}_2=\frac1{\sqrt5}(1,-2)$.
>
> **Step 3 — rotate.** With $\mathbf{x}=u\mathbf{q}_1+v\mathbf{q}_2$, the equation becomes $6u^2+v^2=6$, i.e.
> $$
> \frac{u^2}{1}+\frac{v^2}{6}=1 :
> $$
> an ellipse with semi-major axis $\sqrt6\approx2.449$ along $\mathbf{q}_2$ (the direction of the **smaller** eigenvalue) and semi-minor axis $1$ along $\mathbf{q}_1$.
>
> **Step 4 — check a vertex.** $\sqrt6\,\mathbf{q}_2=(1.0954,-2.1909)$: $5(1.2)+4(1.0954)(-2.1909)+2(4.8)=6-9.6+9.6=6$ $\checkmark$. The tilt of the major axis is $\arctan(-2)=-63.43^\circ$.

```desmos-graph
left=-3; right=3; bottom=-3; top=3
---
5x^2+4xy+2y^2=6|blue
y=0.5x|red|dotted
y=-2x|green|dotted
(0.8944,0.4472)|label:q1, semi-axis 1
(1.0954,-2.1909)|label:√6 q2
6x^2+y^2=6|grey|dashed
```
*The conic $5x^2+4xy+2y^2=6$ (blue) and its principal axes along the eigenvectors $\mathbf{q}_1\propto(2,1)$ ($\lambda=6$, red) and $\mathbf{q}_2\propto(1,-2)$ ($\lambda=1$, green). The dashed ellipse $6u^2+v^2=6$ is the same curve drawn in principal-axis coordinates before rotation: the eigenvector change of variables simply rotates it into place.*

> [!example] Worked Example 10.21 — A Hermitian matrix
> Find the eigenvalues and orthonormal eigenvectors of $\mathbf{H}=\begin{bmatrix}2&1+j\\1-j&3\end{bmatrix}$.
>
> **Step 1 — eigenvalues.** $\operatorname{tr}\mathbf{H}=5$, $\det\mathbf{H}=6-(1+j)(1-j)=6-2=4$. So $\lambda^2-5\lambda+4=0$: $\lambda=4,1$ — real, as promised, even though $\mathbf{H}$ is complex.
>
> **Step 2 — eigenvectors.** $\lambda=4$: $\mathbf{H}-4\mathbf{I}=\begin{bmatrix}-2&1+j\\1-j&-1\end{bmatrix}$. From the first row, $\mathbf{v}_1=(1+j,\,2)$; check the second row: $(1-j)(1+j)-2=0$ $\checkmark$. $\lambda=1$: first row $v_1+(1+j)v_2=0$, so $\mathbf{v}_2=(1+j,\,-1)$.
>
> **Step 3 — orthogonality with the Hermitian inner product.** $\mathbf{v}_1^\mathsf{H}\mathbf{v}_2=(1-j)(1+j)+2(-1)=2-2=0$ $\checkmark$. Norms: $\lVert\mathbf{v}_1\rVert^2=2+4=6$, $\lVert\mathbf{v}_2\rVert^2=2+1=3$.
> $$
> \mathbf{U}=\begin{bmatrix}(1+j)/\sqrt6&(1+j)/\sqrt3\\2/\sqrt6&-1/\sqrt3\end{bmatrix},\qquad \mathbf{H}=\mathbf{U}\begin{bmatrix}4&0\\0&1\end{bmatrix}\mathbf{U}^\mathsf{H}.
> $$
> **Check.** $\mathbf{H}\mathbf{v}_1=\big(2(1+j)+2(1+j),\ (1-j)(1+j)+6\big)=(4+4j,\,8)=4\mathbf{v}_1$ $\checkmark$. Hermitian matrices of this kind are the spatial covariance matrices of two-antenna receivers; the eigenvalue $4$ is the power captured by the best beamformer $\mathbf{u}_1$, and $1$ the power in the orthogonal direction.

> [!example] Worked Example 10.22 — Rayleigh's quotient as an eigenvalue estimator
> Estimate the largest eigenvalue of the tridiagonal $\mathbf{B}$ of Example 10.2 ($\lambda_{\max}=2+\sqrt2=3.414214$) from trial vectors.
>
> - $\mathbf{x}=(1,1,1)$: $\mathbf{B}\mathbf{x}=(3,4,3)$, $R=10/3=3.3333$ — error $0.081$.
> - $\mathbf{x}=(1,1.4,1)$ (a better mode-shape guess, centre bulging): $\mathbf{B}\mathbf{x}=(3.4,4.8,3.4)$, $\mathbf{x}^\mathsf{T}\mathbf{B}\mathbf{x}=3.4+6.72+3.4=13.52$, $\mathbf{x}^\mathsf{T}\mathbf{x}=3.96$, $R=3.414141$ — error $7\times10^{-5}$.
>
> The exact eigenvector is $(1,\sqrt2,1)$; the second trial vector has its middle entry about $1\%$ too small, yet the eigenvalue is accurate to five significant figures. Both estimates are **below** $\lambda_{\max}$, as Rayleigh's principle guarantees. This quadratic accuracy is why a rough guess of a mode shape gives a good resonant frequency.

> [!warning] Common Pitfall — "symmetric" complex matrices
> A complex matrix with $\mathbf{A}^\mathsf{T}=\mathbf{A}$ (without conjugation) is **complex symmetric**, not Hermitian, and none of the theorems of this section apply: $\begin{bmatrix}1&j\\j&-1\end{bmatrix}$ is complex symmetric with the double eigenvalue $0$ and only one eigenvector — it is defective. Impedance and admittance matrices of *lossy* reciprocal networks are complex symmetric. Their eigenvalues are generally complex and the eigenvectors are orthogonal only in the unconjugated sense $\mathbf{v}_i^\mathsf{T}\mathbf{v}_k=0$. Check which transpose the theorem needs before using it.

> [!tip] Engineering Note — principal component analysis of sensor data
> Given $N$ zero-mean measurement vectors $\mathbf{x}_k\in\mathbb{R}^n$ (for example, the outputs of $n$ microphones), the sample covariance $\mathbf{C}=\frac1N\sum_k\mathbf{x}_k\mathbf{x}_k^\mathsf{T}$ is symmetric with non-negative eigenvalues. Its spectral decomposition sorts the signal space into orthogonal directions by power: $\lambda_i$ is the average power along $\mathbf{q}_i$. Projecting onto the top few eigenvectors keeps most of the energy and discards mostly noise — principal component analysis, also known as the Karhunen–Loève transform. The statistics side is in [Ch. 20](../part-6-probability-statistics/20-statistical-analysis.md) and [Ch. 21](../part-6-probability-statistics/21-random-processes.md); the mathematics is the spectral theorem.

## 10.6 Positive-definite matrices

### Motivation

Energy is never negative. The energy stored in a set of coupled inductors, $W=\tfrac12\mathbf{i}^\mathsf{T}\mathbf{L}\mathbf{i}$, is positive for every nonzero current vector; so is the electrostatic energy $\tfrac12\mathbf{v}^\mathsf{T}\mathbf{C}\mathbf{v}$ of a system of conductors, and the power $\mathbf{v}^\mathsf{T}\mathbf{G}\mathbf{v}$ dissipated in a passive resistive network. The variance of any linear combination $\mathbf{w}^\mathsf{T}\mathbf{x}$ of random voltages is $\mathbf{w}^\mathsf{T}\mathbf{C}_{xx}\mathbf{w}\ge0$. A smooth function has a strict local minimum where its gradient vanishes and its Hessian quadratic form is positive for every displacement ([Ch. 6](../part-2-calculus/06-multivariable-calculus.md), [Ch. 23](../part-7-numerical-computational/23-optimization.md)). In all these cases a symmetric matrix $\mathbf{A}$ has the property $\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}>0$ — it is **positive definite**. Such matrices are the "positive numbers" of matrix algebra: they have square roots, they can be factored as $\mathbf{R}^\mathsf{T}\mathbf{R}$, linear systems with them can be solved stably without pivoting, and they define inner products and distances.

### Definitions

> [!definition] Definiteness
> A real symmetric (or complex Hermitian) $n\times n$ matrix $\mathbf{A}$ is
> - **positive definite** (PD, $\mathbf{A}\succ0$) if $\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}>0$ for every $\mathbf{x}\ne\mathbf{0}$;
> - **positive semidefinite** (PSD, $\mathbf{A}\succeq0$) if $\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}\ge0$ for every $\mathbf{x}$;
> - **negative definite / semidefinite** if $-\mathbf{A}$ is positive definite / semidefinite;
> - **indefinite** if $\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}$ takes both signs.
>
> For complex Hermitian matrices replace $\mathbf{x}^\mathsf{T}$ by $\mathbf{x}^\mathsf{H}$; the form $\mathbf{x}^\mathsf{H}\mathbf{A}\mathbf{x}$ is then real.

Definiteness is normally defined only for symmetric matrices, because the quadratic form only sees the symmetric part $\tfrac12(\mathbf{A}+\mathbf{A}^\mathsf{T})$.

### Equivalent tests

> [!theorem] Tests for positive definiteness
> For a real symmetric $\mathbf{A}$, the following are equivalent:
> 1. $\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}>0$ for all $\mathbf{x}\ne\mathbf{0}$;
> 2. **all eigenvalues of $\mathbf{A}$ are positive**;
> 3. **all leading principal minors are positive** (Sylvester's criterion): $\Delta_1=a_{11}>0$, $\Delta_2=\det\begin{bmatrix}a_{11}&a_{12}\\a_{21}&a_{22}\end{bmatrix}>0$, …, $\Delta_n=\det\mathbf{A}>0$;
> 4. **all pivots** of Gaussian elimination without row exchanges are positive;
> 5. $\mathbf{A}=\mathbf{R}^\mathsf{T}\mathbf{R}$ for some matrix $\mathbf{R}$ with independent columns — in particular, the **Cholesky factorization** $\mathbf{A}=\mathbf{L}\mathbf{L}^\mathsf{T}$ exists with $\mathbf{L}$ lower triangular with positive diagonal.

*Proof sketch.*
(1)$\Leftrightarrow$(2): with $\mathbf{x}=\mathbf{Q}\mathbf{y}$, $\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}=\sum\lambda_iy_i^2$; this is positive for all $\mathbf{y}\ne\mathbf{0}$ exactly when every $\lambda_i>0$ (take $\mathbf{y}=\mathbf{e}_i$ for necessity).
(2)$\Rightarrow$(3): the upper-left $k\times k$ block $\mathbf{A}_k$ is itself PD (restrict (1) to vectors whose last $n-k$ entries are zero), so its eigenvalues are positive and $\Delta_k=\det\mathbf{A}_k$ — their product — is positive.
(3)$\Leftrightarrow$(4): elimination without row exchanges gives pivots $d_k=\Delta_k/\Delta_{k-1}$ (with $\Delta_0=1$), because the determinant of each leading block is the product of its pivots.
(4)$\Rightarrow$(5): elimination gives $\mathbf{A}=\mathbf{L}_1\mathbf{D}\mathbf{L}_1^\mathsf{T}$ with unit lower-triangular $\mathbf{L}_1$ and $\mathbf{D}=\operatorname{diag}(d_k)$; set $\mathbf{L}=\mathbf{L}_1\mathbf{D}^{1/2}$.
(5)$\Rightarrow$(1): $\mathbf{x}^\mathsf{T}\mathbf{R}^\mathsf{T}\mathbf{R}\mathbf{x}=\lVert\mathbf{R}\mathbf{x}\rVert^2\ge0$, with equality only if $\mathbf{R}\mathbf{x}=\mathbf{0}$, i.e. $\mathbf{x}=\mathbf{0}$ because $\mathbf{R}$ has independent columns. $\square$

For semidefiniteness the tests weaken to: all eigenvalues $\ge0$; $\mathbf{A}=\mathbf{R}^\mathsf{T}\mathbf{R}$ for *some* $\mathbf{R}$; or **all** principal minors (not just the leading ones) $\ge0$.

In practice: for a $2\times2$ matrix, check $a_{11}>0$ and $\det\mathbf{A}>0$. For $3\times3$ by hand, use leading minors or pivots. For large matrices, attempt a Cholesky factorization — it succeeds if and only if the matrix is PD, and it costs about half as much as an LU factorization.

### Gram matrices and the energy form

The test (5) has a very useful reading. For **any** $m\times n$ matrix $\mathbf{R}$, the **Gram matrix** $\mathbf{R}^\mathsf{T}\mathbf{R}$ is symmetric positive semidefinite, and it is positive definite iff the columns of $\mathbf{R}$ are linearly independent. The normal-equation matrix $\mathbf{A}^\mathsf{T}\mathbf{A}$ of least squares ([Ch. 9](09-vector-spaces-and-linear-transformations.md)), sample covariance matrices and the matrices $\mathbf{A}^\mathsf{T}\mathbf{A}$, $\mathbf{A}\mathbf{A}^\mathsf{T}$ of the SVD (Section 10.7) are all Gram matrices.

Physical matrices often arrive in exactly this form. The conductance matrix of a resistor network is $\mathbf{G}=\mathbf{B}^\mathsf{T}\operatorname{diag}(g_e)\mathbf{B}$, where $\mathbf{B}$ is the branch–node incidence matrix (each branch voltage is a difference of node voltages) and $g_e>0$ are the branch conductances; the dissipated power $\mathbf{v}^\mathsf{T}\mathbf{G}\mathbf{v}=\sum_eg_e(\text{branch voltage})^2$ is a sum of squares. Grounding at least one node makes $\mathbf{B}$'s columns independent and $\mathbf{G}$ positive definite.

### Geometry: level curves

The level sets $\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}=c$ show definiteness at a glance. In principal-axis coordinates they read $\sum\lambda_iy_i^2=c$:

- **PD**: ellipses (ellipsoids) centred at the origin; $Q$ is a bowl with a unique minimum.
- **PSD but singular**: a zero eigenvalue makes the ellipse degenerate into a pair of parallel lines (a cylinder in 3-D); $Q$ is a trough, flat along the null-space direction.
- **Indefinite**: hyperbolas; $Q$ is a saddle.

```desmos-graph
left=-3; right=3; bottom=-2.5; top=2.5
---
x^2-xy+y^2=1|blue
(x+y)^2=1|green
x^2+4xy+y^2=1|red
x^2+4xy+y^2=-1|red|dashed
```
*Level curves of three quadratic forms: positive definite $x^2-xy+y^2=1$ (blue ellipse, eigenvalues $1.5$ and $0.5$); positive semidefinite $(x+y)^2=1$ (green parallel lines, eigenvalues $2$ and $0$); indefinite $x^2+4xy+y^2=\pm1$ (red hyperbolas, eigenvalues $3$ and $-1$).*

### Square roots and the energy norm

A PSD matrix has a unique PSD **square root** $\mathbf{A}^{1/2}=\mathbf{Q}\boldsymbol{\Lambda}^{1/2}\mathbf{Q}^\mathsf{T}$, with $\boldsymbol{\Lambda}^{1/2}=\operatorname{diag}(\sqrt{\lambda_i})$; indeed $\mathbf{A}^{1/2}\mathbf{A}^{1/2}=\mathbf{Q}\boldsymbol{\Lambda}\mathbf{Q}^\mathsf{T}=\mathbf{A}$. Both $\mathbf{A}^{1/2}$ and the Cholesky factor $\mathbf{L}^\mathsf{T}$ are "square roots" in the sense $\mathbf{A}=\mathbf{R}^\mathsf{T}\mathbf{R}$; they differ by an orthogonal factor. A PD matrix also defines an inner product $\langle\mathbf{x},\mathbf{y}\rangle_\mathbf{A}=\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{y}$ and an **energy norm** $\lVert\mathbf{x}\rVert_\mathbf{A}=\sqrt{\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}}$; the Mahalanobis distance of statistics is the energy norm of the inverse covariance, and "whitening" noise means multiplying by $\mathbf{C}^{-1/2}$ or $\mathbf{L}^{-1}$ so that the noise covariance becomes $\mathbf{I}$.

### Worked examples

> [!example] Worked Example 10.23 — Three tests on one matrix
> Show that $\mathbf{K}=\begin{bmatrix}2&-1&0\\-1&2&-1\\0&-1&2\end{bmatrix}$ is positive definite, in three ways.
>
> **Leading minors.** $\Delta_1=2$, $\Delta_2=4-1=3$, $\Delta_3=2(4-1)-(-1)(-2-0)=6-2=4$. All positive $\Rightarrow$ PD.
>
> **Pivots.** $d_k=\Delta_k/\Delta_{k-1}$: $2,\ 3/2,\ 4/3$ — all positive. (Elimination directly: $R_2+\tfrac12R_1$ gives pivot $2-\tfrac12=\tfrac32$; $R_3+\tfrac23R_2$ gives $2-\tfrac23=\tfrac43$ $\checkmark$.)
>
> **Eigenvalues.** $\mathbf{K}=\mathbf{S}\mathbf{B}\mathbf{S}$ with $\mathbf{S}=\operatorname{diag}(1,-1,1)=\mathbf{S}^{-1}$ and $\mathbf{B}$ the matrix of Example 10.2, so $\mathbf{K}$ is similar to $\mathbf{B}$ and has eigenvalues $2-\sqrt2,\ 2,\ 2+\sqrt2$, all positive. (Check: product $=4=\Delta_3$ $\checkmark$.)
>
> **Energy.** Better than any test is to *see* the sum of squares:
> $$
> \mathbf{x}^\mathsf{T}\mathbf{K}\mathbf{x}=2x_1^2+2x_2^2+2x_3^2-2x_1x_2-2x_2x_3=x_1^2+(x_1-x_2)^2+(x_2-x_3)^2+x_3^2 ,
> $$
> which vanishes only if $x_1=0$, $x_1=x_2$, $x_2=x_3$, $x_3=0$. This is the stiffness (or conductance) matrix of a chain of four identical springs (resistors) between two fixed walls (grounded ends): the four squares are the four element energies.

> [!example] Worked Example 10.24 — Cholesky factorization
> Factor $\mathbf{A}=\begin{bmatrix}4&2&-2\\2&10&2\\-2&2&6\end{bmatrix}$ as $\mathbf{L}\mathbf{L}^\mathsf{T}$.
>
> Write $\mathbf{L}=\begin{bmatrix}l_{11}&0&0\\l_{21}&l_{22}&0\\l_{31}&l_{32}&l_{33}\end{bmatrix}$ and match entries of $\mathbf{L}\mathbf{L}^\mathsf{T}$ column by column:
> - $l_{11}^2=4\Rightarrow l_{11}=2$; $\ l_{21}l_{11}=2\Rightarrow l_{21}=1$; $\ l_{31}l_{11}=-2\Rightarrow l_{31}=-1$.
> - $l_{21}^2+l_{22}^2=10\Rightarrow l_{22}=\sqrt{9}=3$; $\ l_{31}l_{21}+l_{32}l_{22}=2\Rightarrow l_{32}=(2+1)/3=1$.
> - $l_{31}^2+l_{32}^2+l_{33}^2=6\Rightarrow l_{33}=\sqrt{6-1-1}=2$.
> $$
> \mathbf{L}=\begin{bmatrix}2&0&0\\1&3&0\\-1&1&2\end{bmatrix}.
> $$
> Every square root had a positive argument, which certifies that $\mathbf{A}$ is PD. **By-products:** $\det\mathbf{A}=(2\cdot3\cdot2)^2=144$, and the pivots are $l_{kk}^2=4,9,4$. To solve $\mathbf{A}\mathbf{x}=\mathbf{b}$, solve $\mathbf{L}\mathbf{y}=\mathbf{b}$ forward and $\mathbf{L}^\mathsf{T}\mathbf{x}=\mathbf{y}$ backward. **Check:** row 3 of $\mathbf{L}$ dotted with itself, $1+1+4=6=a_{33}$ $\checkmark$; with row 2, $-1+3+0=2=a_{32}$ $\checkmark$.

> [!example] Worked Example 10.25 — For which parameter values is a matrix PD?
> For which real $k$ is $\mathbf{A}(k)=\begin{bmatrix}1&k&0\\k&1&k\\0&k&1\end{bmatrix}$ positive definite?
>
> **Leading minors.** $\Delta_1=1>0$; $\Delta_2=1-k^2>0\iff\lvert k\rvert<1$; $\Delta_3=1(1-k^2)-k(k-0)=1-2k^2>0\iff\lvert k\rvert<1/\sqrt2$. The binding condition is the last:
> $$
> \mathbf{A}(k)\succ0\iff\lvert k\rvert<\tfrac{1}{\sqrt2}\approx0.7071 .
> $$
> **Eigenvalue check.** Expanding $\det(\mathbf{A}-\lambda\mathbf{I})=(1-\lambda)[(1-\lambda)^2-2k^2]$ gives $\lambda=1,\ 1\pm\sqrt2\,k$. The smallest, $1-\sqrt2\lvert k\rvert$, is positive exactly when $\lvert k\rvert<1/\sqrt2$ $\checkmark$. At $k=1/\sqrt2$ the matrix is PSD and singular, with null vector $(1,-\sqrt2,1)$.
>
> **Moral.** Positivity of all entries, or diagonal entries larger than off-diagonal ones (here $1>k$ for $k=0.8$), does **not** imply definiteness; at $k=0.8$, $\Delta_3=-0.28<0$.

> [!example] Worked Example 10.26 — Coupled inductors: the energy condition $k<1$
> Two coupled coils have inductance matrix $\mathbf{L}=\begin{bmatrix}L_1&M\\M&L_2\end{bmatrix}$ and stored energy $W=\tfrac12\mathbf{i}^\mathsf{T}\mathbf{L}\mathbf{i}=\tfrac12(L_1i_1^2+2Mi_1i_2+L_2i_2^2)$.
>
> **Condition.** Passivity requires $W>0$ for all nonzero currents: $\mathbf{L}\succ0\iff L_1>0$ and $L_1L_2-M^2>0$, i.e. the **coupling coefficient** $k=\lvert M\rvert/\sqrt{L_1L_2}<1$.
>
> **Numbers.** $L_1=4$ mH, $L_2=9$ mH, $M=5$ mH: $k=5/6=0.833$, $\det\mathbf{L}=36-25=11$ mH$^2>0$. Eigenvalues: $\lambda=\tfrac12\big(13\pm\sqrt{169-44}\big)=12.09$ mH and $0.910$ mH. For opposing currents $\mathbf{i}=(1,-1)$ A, $W=\tfrac12(4-10+9)=1.5$ mJ $>0$, and indeed no current vector gives less than $\tfrac12(0.910\text{ mH})\lVert\mathbf{i}\rVert^2$.
>
> **The limit $k=1$.** With $M=6$ mH, $\det\mathbf{L}=0$: PSD, not PD. The null vector $\mathbf{i}=(3,-2)$ A gives $W=\tfrac12(4\cdot9-2\cdot6\cdot6+9\cdot4)=0$: the two coils' fluxes cancel exactly. An ideal transformer is this singular limit; a value $k>1$ in a simulator netlist describes a device that could deliver energy from nothing, and the simulation will typically blow up.

> [!warning] Common Pitfall — leading minors and semidefiniteness
> Sylvester's criterion with $\ge$ instead of $>$ does **not** test semidefiniteness. $\begin{bmatrix}0&0\\0&-1\end{bmatrix}$ has leading minors $0$ and $0$, both $\ge0$, yet it is negative semidefinite. For PSD you must check *all* principal minors (here the $1\times1$ minor $-1$ fails) or the eigenvalues. Also, for **negative** definiteness the leading minors must *alternate* in sign, $\Delta_1<0,\ \Delta_2>0,\ \Delta_3<0,\dots$ — not all be negative.

> [!tip] Engineering Note — positive definiteness as a modelling check
> Numerically extracted parasitic matrices (field-solver $\mathbf{L}$, $\mathbf{C}$ and $\mathbf{G}$ for packages and interconnect), measured covariance matrices with missing samples, and estimated Hessians in optimisation all lose positive definiteness through rounding or estimation error. A non-PD inductance or capacitance matrix makes a passive model non-passive, and transient simulators then produce growing oscillations. Standard practice is to test with a Cholesky attempt and, if it fails, repair the matrix by clipping negative eigenvalues to a small positive floor: $\mathbf{A}\leftarrow\mathbf{Q}\max(\boldsymbol{\Lambda},\epsilon\mathbf{I})\mathbf{Q}^\mathsf{T}$, the nearest PD matrix in the Frobenius norm.

## 10.7 Singular value decomposition

### Motivation

Eigen-decomposition has three limitations: it applies only to square matrices; it may fail (defective matrices); and, for non-symmetric matrices, the eigenvectors are not orthogonal, so $\mathbf{P}$ can be badly conditioned. Many engineering matrices are not square — a $4\times2$ MIMO channel matrix from 2 transmit to 4 receive antennas, a $1000\times3$ regression matrix of measurements against model parameters, an $m\times n$ image. The **singular value decomposition** (SVD) repairs all three defects at once. Every matrix, of any shape, factors as

$$
\mathbf{A}=\mathbf{U}\boldsymbol{\Sigma}\mathbf{V}^\mathsf{T}
$$

with **orthogonal** $\mathbf{U}$ and $\mathbf{V}$ and a non-negative diagonal $\boldsymbol{\Sigma}$. The price is that we use two different orthonormal bases — one for the input space and one for the output space — instead of one eigenbasis. The SVD is arguably the most useful single factorization in applied linear algebra: it gives the rank, the norm and the condition number of a matrix, the best low-rank approximation, the pseudoinverse and the least-squares solution, and the capacity-achieving decomposition of a MIMO channel into parallel independent sub-channels.

### Statement

> [!theorem] Singular value decomposition
> Every real $m\times n$ matrix $\mathbf{A}$ of rank $r$ can be factored as
> $$
> \mathbf{A}=\mathbf{U}\boldsymbol{\Sigma}\mathbf{V}^\mathsf{T},
> $$
> where $\mathbf{U}$ ($m\times m$) and $\mathbf{V}$ ($n\times n$) are orthogonal and $\boldsymbol{\Sigma}$ ($m\times n$) is zero except for its leading diagonal entries
> $$
> \sigma_1\ge\sigma_2\ge\dots\ge\sigma_r>0=\sigma_{r+1}=\dots .
> $$
> The $\sigma_i$ are the **singular values**, the columns $\mathbf{u}_i$ of $\mathbf{U}$ the **left singular vectors** and the columns $\mathbf{v}_i$ of $\mathbf{V}$ the **right singular vectors**. They satisfy
> $$
> \mathbf{A}\mathbf{v}_i=\sigma_i\mathbf{u}_i,\qquad \mathbf{A}^\mathsf{T}\mathbf{u}_i=\sigma_i\mathbf{v}_i,
> $$
> and $\sigma_i^2$ are the nonzero eigenvalues of both $\mathbf{A}^\mathsf{T}\mathbf{A}$ and $\mathbf{A}\mathbf{A}^\mathsf{T}$. For complex $\mathbf{A}$, replace transposes by conjugate transposes and "orthogonal" by "unitary"; the $\sigma_i$ are still real and non-negative.

The singular values are unique; the singular vectors are unique up to sign (or a phase factor, in the complex case) when the $\sigma_i$ are distinct.

### Derivation

The construction uses only the spectral theorem. The $n\times n$ matrix $\mathbf{A}^\mathsf{T}\mathbf{A}$ is symmetric and positive semidefinite (a Gram matrix, Section 10.6), so it has orthonormal eigenvectors $\mathbf{v}_1,\dots,\mathbf{v}_n$ with eigenvalues $\mu_1\ge\dots\ge\mu_n\ge0$. Exactly $r=\operatorname{rank}\mathbf{A}$ of them are positive, because $\operatorname{null}(\mathbf{A}^\mathsf{T}\mathbf{A})=\operatorname{null}(\mathbf{A})$. Define

$$
\sigma_i=\sqrt{\mu_i},\qquad \mathbf{u}_i=\frac{\mathbf{A}\mathbf{v}_i}{\sigma_i}\quad(i=1,\dots,r).
$$

These $\mathbf{u}_i$ are orthonormal:

$$
\mathbf{u}_i^\mathsf{T}\mathbf{u}_k=\frac{\mathbf{v}_i^\mathsf{T}\mathbf{A}^\mathsf{T}\mathbf{A}\mathbf{v}_k}{\sigma_i\sigma_k}=\frac{\mu_k\,\mathbf{v}_i^\mathsf{T}\mathbf{v}_k}{\sigma_i\sigma_k}=\begin{cases}1&i=k\\0&i\ne k.\end{cases}
$$

Extend $\mathbf{u}_1,\dots,\mathbf{u}_r$ to an orthonormal basis of $\mathbb{R}^m$ (Gram–Schmidt). Then $\mathbf{A}\mathbf{v}_i=\sigma_i\mathbf{u}_i$ for $i\le r$ and $\mathbf{A}\mathbf{v}_i=\mathbf{0}$ for $i>r$; in matrix form $\mathbf{A}\mathbf{V}=\mathbf{U}\boldsymbol{\Sigma}$, and multiplying on the right by $\mathbf{V}^\mathsf{T}=\mathbf{V}^{-1}$ gives the SVD. $\square$

The derivation is also a hand-computation recipe: eigen-decompose the smaller of $\mathbf{A}^\mathsf{T}\mathbf{A}$ and $\mathbf{A}\mathbf{A}^\mathsf{T}$, then get the other set of singular vectors from $\mathbf{u}_i=\mathbf{A}\mathbf{v}_i/\sigma_i$ (or $\mathbf{v}_i=\mathbf{A}^\mathsf{T}\mathbf{u}_i/\sigma_i$). Do **not** compute $\mathbf{U}$ and $\mathbf{V}$ independently from the two Gram matrices: each eigenvector is only determined up to sign, and independently chosen signs may not match.

### Geometry

Reading $\mathbf{A}\mathbf{x}=\mathbf{U}\boldsymbol{\Sigma}\mathbf{V}^\mathsf{T}\mathbf{x}$ right to left, **every linear map is a rotation (or reflection), followed by axis-aligned stretching, followed by another rotation**. In particular, $\mathbf{A}$ maps the unit sphere in $\mathbb{R}^n$ onto an ellipsoid (possibly flattened) in $\mathbb{R}^m$ whose semi-axes are $\sigma_i\mathbf{u}_i$; the right singular vectors $\mathbf{v}_i$ are the orthogonal input directions that map onto those axes. Compare the graph in Section 10.2: eigenvectors are directions that are not rotated; singular vectors are orthogonal directions that *stay* orthogonal.

### The SVD and the four fundamental subspaces

The SVD supplies orthonormal bases for all four fundamental subspaces of [Ch. 9](09-vector-spaces-and-linear-transformations.md) at once:

| Subspace | Orthonormal basis | Dimension |
|---|---|---|
| column space $\mathcal{C}(\mathbf{A})$ | $\mathbf{u}_1,\dots,\mathbf{u}_r$ | $r$ |
| left null space $\mathcal{N}(\mathbf{A}^\mathsf{T})$ | $\mathbf{u}_{r+1},\dots,\mathbf{u}_m$ | $m-r$ |
| row space $\mathcal{C}(\mathbf{A}^\mathsf{T})$ | $\mathbf{v}_1,\dots,\mathbf{v}_r$ | $r$ |
| null space $\mathcal{N}(\mathbf{A})$ | $\mathbf{v}_{r+1},\dots,\mathbf{v}_n$ | $n-r$ |

Keeping only the first $r$ columns gives the **compact** (or "thin") SVD and the **outer-product expansion**

$$
\mathbf{A}=\mathbf{U}_r\boldsymbol{\Sigma}_r\mathbf{V}_r^\mathsf{T}=\sigma_1\mathbf{u}_1\mathbf{v}_1^\mathsf{T}+\sigma_2\mathbf{u}_2\mathbf{v}_2^\mathsf{T}+\dots+\sigma_r\mathbf{u}_r\mathbf{v}_r^\mathsf{T},
$$

a sum of $r$ rank-one matrices in decreasing order of importance.

### Norms, condition number and rank

The SVD measures matrices. Because orthogonal matrices preserve lengths,

$$
\lVert\mathbf{A}\rVert_2=\max_{\mathbf{x}\ne\mathbf{0}}\frac{\lVert\mathbf{A}\mathbf{x}\rVert}{\lVert\mathbf{x}\rVert}=\sigma_1,
\qquad
\lVert\mathbf{A}\rVert_F=\sqrt{\sum_{i,k}a_{ik}^2}=\sqrt{\sigma_1^2+\dots+\sigma_r^2},
$$

the spectral norm (maximum gain) and the Frobenius norm. For a square invertible matrix, $\lvert\det\mathbf{A}\rvert=\sigma_1\sigma_2\cdots\sigma_n$, $\lVert\mathbf{A}^{-1}\rVert_2=1/\sigma_n$, and the **condition number** introduced in [Ch. 8](08-matrices-and-linear-systems.md) is

$$
\kappa_2(\mathbf{A})=\lVert\mathbf{A}\rVert_2\lVert\mathbf{A}^{-1}\rVert_2=\frac{\sigma_{\max}}{\sigma_{\min}} .
$$

It bounds error amplification in solving $\mathbf{A}\mathbf{x}=\mathbf{b}$: a relative perturbation $\epsilon$ in $\mathbf{b}$ can cause a relative change of up to $\kappa\epsilon$ in $\mathbf{x}$. Rule of thumb: you lose about $\log_{10}\kappa$ significant digits.

The **numerical rank** is the number of singular values above a tolerance (typically $\max(m,n)\,\sigma_1\,\epsilon_{\text{mach}}$). Unlike the exact rank, which jumps discontinuously under perturbation, the singular values move continuously: $\lvert\sigma_i(\mathbf{A}+\mathbf{E})-\sigma_i(\mathbf{A})\rvert\le\lVert\mathbf{E}\rVert_2$. The smallest singular value is the **distance to the nearest singular matrix**: $\min\{\lVert\mathbf{E}\rVert_2:\mathbf{A}+\mathbf{E}\text{ singular}\}=\sigma_n$.

> [!theorem] Eckart–Young (best low-rank approximation)
> For $k<r$, the truncated sum $\mathbf{A}_k=\sum_{i=1}^k\sigma_i\mathbf{u}_i\mathbf{v}_i^\mathsf{T}$ is the closest matrix of rank at most $k$ to $\mathbf{A}$, in both the spectral and Frobenius norms:
> $$
> \min_{\operatorname{rank}\mathbf{B}\le k}\lVert\mathbf{A}-\mathbf{B}\rVert_2=\sigma_{k+1},\qquad
> \min_{\operatorname{rank}\mathbf{B}\le k}\lVert\mathbf{A}-\mathbf{B}\rVert_F=\sqrt{\sigma_{k+1}^2+\dots+\sigma_r^2}.
> $$

This underlies image and data compression, noise reduction by truncating small singular values, and model-order reduction of large circuits (balanced truncation discards states associated with small *Hankel* singular values).

### The pseudoinverse and least squares

> [!definition] Moore–Penrose pseudoinverse
> If $\mathbf{A}=\mathbf{U}_r\boldsymbol{\Sigma}_r\mathbf{V}_r^\mathsf{T}$ (compact SVD), its **pseudoinverse** is
> $$
> \mathbf{A}^+=\mathbf{V}_r\boldsymbol{\Sigma}_r^{-1}\mathbf{U}_r^\mathsf{T}=\sum_{i=1}^r\frac{1}{\sigma_i}\mathbf{v}_i\mathbf{u}_i^\mathsf{T}.
> $$

For any $\mathbf{b}$, $\mathbf{x}^+=\mathbf{A}^+\mathbf{b}$ is the **minimum-norm least-squares solution**: among all $\mathbf{x}$ that minimise $\lVert\mathbf{A}\mathbf{x}-\mathbf{b}\rVert$, it is the one with smallest $\lVert\mathbf{x}\rVert$. *Why:* in the rotated coordinates $\mathbf{y}=\mathbf{V}^\mathsf{T}\mathbf{x}$, $\mathbf{c}=\mathbf{U}^\mathsf{T}\mathbf{b}$, the residual is $\lVert\boldsymbol{\Sigma}\mathbf{y}-\mathbf{c}\rVert^2=\sum_{i\le r}(\sigma_iy_i-c_i)^2+\sum_{i>r}c_i^2$. The first sum is zeroed by $y_i=c_i/\sigma_i$; the second cannot be changed; the remaining $y_i$ ($i>r$) do not affect the residual and are set to zero to minimise the norm. Special cases: if $\mathbf{A}$ has independent columns, $\mathbf{A}^+=(\mathbf{A}^\mathsf{T}\mathbf{A})^{-1}\mathbf{A}^\mathsf{T}$ (the normal-equation solution of [Ch. 9](09-vector-spaces-and-linear-transformations.md)); if it has independent rows, $\mathbf{A}^+=\mathbf{A}^\mathsf{T}(\mathbf{A}\mathbf{A}^\mathsf{T})^{-1}$; if it is square and invertible, $\mathbf{A}^+=\mathbf{A}^{-1}$.

The formula also exposes the danger of tiny singular values: noise in $c_i$ is amplified by $1/\sigma_i$. **Truncated SVD** (drop terms with $\sigma_i$ below a threshold) and **Tikhonov regularisation** (replace $1/\sigma_i$ by $\sigma_i/(\sigma_i^2+\alpha^2)$) tame this amplification, at the cost of a small bias.

### Singular values versus eigenvalues

For a symmetric matrix the SVD and the spectral decomposition coincide up to signs: $\sigma_i=\lvert\lambda_i\rvert$, and a negative $\lambda_i$ is absorbed by flipping $\mathbf{u}_i=-\mathbf{q}_i$. For a symmetric PSD matrix they are identical. For non-symmetric matrices the two sets of numbers can be very different. For $\mathbf{A}=\begin{bmatrix}4&1\\2&3\end{bmatrix}$ the eigenvalues are $5,2$ but the singular values are $5.117,1.954$ (square roots of the eigenvalues $15\pm\sqrt{125}$ of $\mathbf{A}\mathbf{A}^\mathsf{T}$); their products agree ($\lvert\det\mathbf{A}\rvert=10$) but nothing else need. Always $\sigma_{\min}\le\lvert\lambda_i\rvert\le\sigma_{\max}$. Eigenvalues govern what happens under *repeated* application of $\mathbf{A}$ (powers, dynamics); singular values govern what happens under a *single* application (gain, sensitivity, rank).

### Worked examples

> [!example] Worked Example 10.27 — SVD of a $2\times2$ matrix
> Find the SVD of $\mathbf{A}=\begin{bmatrix}3&0\\4&5\end{bmatrix}$.
>
> **Step 1 — $\mathbf{A}^\mathsf{T}\mathbf{A}$.**
> $$
> \mathbf{A}^\mathsf{T}\mathbf{A}=\begin{bmatrix}3&4\\0&5\end{bmatrix}\begin{bmatrix}3&0\\4&5\end{bmatrix}=\begin{bmatrix}25&20\\20&25\end{bmatrix}.
> $$
> **Step 2 — its eigen-decomposition.** Eigenvalues $25\pm20$: $\mu_1=45$, $\mu_2=5$, with $\mathbf{v}_1=\frac1{\sqrt2}(1,1)$, $\mathbf{v}_2=\frac1{\sqrt2}(-1,1)$. Singular values: $\sigma_1=\sqrt{45}=3\sqrt5\approx6.708$, $\sigma_2=\sqrt5\approx2.236$.
>
> **Step 3 — left singular vectors from $\mathbf{u}_i=\mathbf{A}\mathbf{v}_i/\sigma_i$.**
> $$
> \mathbf{u}_1=\frac{1}{3\sqrt5}\cdot\frac{1}{\sqrt2}\begin{bmatrix}3\\9\end{bmatrix}=\frac{1}{\sqrt{10}}\begin{bmatrix}1\\3\end{bmatrix},\qquad
> \mathbf{u}_2=\frac{1}{\sqrt5}\cdot\frac1{\sqrt2}\begin{bmatrix}-3\\1\end{bmatrix}=\frac1{\sqrt{10}}\begin{bmatrix}-3\\1\end{bmatrix}.
> $$
> They are orthonormal: $\mathbf{u}_1\cdot\mathbf{u}_2=(-3+3)/10=0$ $\checkmark$.
>
> **Step 4 — assemble and verify by outer products.**
> $$
> \sigma_1\mathbf{u}_1\mathbf{v}_1^\mathsf{T}+\sigma_2\mathbf{u}_2\mathbf{v}_2^\mathsf{T}
> =\frac{3\sqrt5}{\sqrt{20}}\begin{bmatrix}1&1\\3&3\end{bmatrix}+\frac{\sqrt5}{\sqrt{20}}\begin{bmatrix}3&-3\\-1&1\end{bmatrix}
> =\begin{bmatrix}1.5&1.5\\4.5&4.5\end{bmatrix}+\begin{bmatrix}1.5&-1.5\\-0.5&0.5\end{bmatrix}=\begin{bmatrix}3&0\\4&5\end{bmatrix}\ \checkmark
> $$
> **Derived quantities.** $\lVert\mathbf{A}\rVert_2=6.708$; $\lVert\mathbf{A}\rVert_F=\sqrt{45+5}=\sqrt{50}$ (directly: $\sqrt{9+16+25}$ $\checkmark$); $\lvert\det\mathbf{A}\rvert=15=\sigma_1\sigma_2$ $\checkmark$; $\kappa_2=3$. Note that the eigenvalues of $\mathbf{A}$ are $3$ and $5$ (triangular), different from the singular values.
>
> **Best rank-1 approximation.** $\mathbf{A}_1=\begin{bmatrix}1.5&1.5\\4.5&4.5\end{bmatrix}$, with $\lVert\mathbf{A}-\mathbf{A}_1\rVert_2=\sigma_2=\sqrt5$: no rank-1 matrix is closer.

```desmos-graph
left=-8; right=8; bottom=-7.5; top=7.5
---
x^2+y^2=1|grey|dashed
41x^2-24xy+9y^2=225|blue
y=3x|red|dotted
y=-\frac{x}{3}|green|dotted
(0.7071,0.7071)|label:v1
(-0.7071,0.7071)|label:v2
(2.1213,6.364)|label:σ1 u1
(-2.1213,0.7071)|label:σ2 u2
```
*The matrix $\mathbf{A}=\begin{bmatrix}3&0\\4&5\end{bmatrix}$ maps the unit circle (dashed) onto the blue ellipse $41x^2-24xy+9y^2=225$. The orthogonal input directions $\mathbf{v}_1,\mathbf{v}_2$ land on the ellipse's semi-axes $\sigma_1\mathbf{u}_1$ (length $6.708$, red line) and $\sigma_2\mathbf{u}_2$ (length $2.236$, green line), which are again orthogonal.*

> [!example] Worked Example 10.28 — SVD of a rectangular matrix
> Find the SVD of $\mathbf{A}=\begin{bmatrix}1&1&0\\0&1&1\end{bmatrix}$ ($2\times3$).
>
> **Step 1 — use the smaller Gram matrix.** $\mathbf{A}\mathbf{A}^\mathsf{T}=\begin{bmatrix}2&1\\1&2\end{bmatrix}$, eigenvalues $3$ and $1$, eigenvectors $\mathbf{u}_1=\frac1{\sqrt2}(1,1)$, $\mathbf{u}_2=\frac1{\sqrt2}(1,-1)$. So $\sigma_1=\sqrt3$, $\sigma_2=1$, rank $2$.
>
> **Step 2 — right singular vectors from $\mathbf{v}_i=\mathbf{A}^\mathsf{T}\mathbf{u}_i/\sigma_i$.**
> $$
> \mathbf{v}_1=\frac{1}{\sqrt3}\cdot\frac{1}{\sqrt2}\begin{bmatrix}1\\2\\1\end{bmatrix}=\frac1{\sqrt6}\begin{bmatrix}1\\2\\1\end{bmatrix},\qquad
> \mathbf{v}_2=\frac{1}{1}\cdot\frac1{\sqrt2}\begin{bmatrix}1\\0\\-1\end{bmatrix}.
> $$
> **Step 3 — complete $\mathbf{V}$.** The third right singular vector spans $\mathcal{N}(\mathbf{A})$: $x_1+x_2=0$, $x_2+x_3=0$ gives $\mathbf{v}_3=\frac1{\sqrt3}(1,-1,1)$. Check orthogonality: $\mathbf{v}_1\cdot\mathbf{v}_3\propto1-2+1=0$, $\mathbf{v}_2\cdot\mathbf{v}_3\propto1-1=0$ $\checkmark$.
>
> **Step 4 — assemble.**
> $$
> \mathbf{A}=\underbrace{\frac1{\sqrt2}\begin{bmatrix}1&1\\1&-1\end{bmatrix}}_{\mathbf{U}}
> \underbrace{\begin{bmatrix}\sqrt3&0&0\\0&1&0\end{bmatrix}}_{\boldsymbol{\Sigma}}
> \underbrace{\begin{bmatrix}1/\sqrt6&2/\sqrt6&1/\sqrt6\\1/\sqrt2&0&-1/\sqrt2\\1/\sqrt3&-1/\sqrt3&1/\sqrt3\end{bmatrix}}_{\mathbf{V}^\mathsf{T}}.
> $$
> **Check** an entry: $(\mathbf{A})_{12}=\sigma_1u_{11}v_{21}+\sigma_2u_{12}v_{22}=\sqrt3\cdot\frac1{\sqrt2}\cdot\frac2{\sqrt6}+1\cdot\frac1{\sqrt2}\cdot0=\frac{2\sqrt3}{\sqrt{12}}=1$ $\checkmark$.
>
> **Minimum-norm solution of $\mathbf{A}\mathbf{x}=(1,2)$.** Independent rows, so $\mathbf{A}^+=\mathbf{A}^\mathsf{T}(\mathbf{A}\mathbf{A}^\mathsf{T})^{-1}$ with $(\mathbf{A}\mathbf{A}^\mathsf{T})^{-1}=\frac13\begin{bmatrix}2&-1\\-1&2\end{bmatrix}$:
> $$
> (\mathbf{A}\mathbf{A}^\mathsf{T})^{-1}\begin{bmatrix}1\\2\end{bmatrix}=\frac13\begin{bmatrix}0\\3\end{bmatrix}=\begin{bmatrix}0\\1\end{bmatrix},\qquad \mathbf{x}^+=\mathbf{A}^\mathsf{T}\begin{bmatrix}0\\1\end{bmatrix}=\begin{bmatrix}0\\1\\1\end{bmatrix}.
> $$
> $\mathbf{A}\mathbf{x}^+=(1,2)$ $\checkmark$. Every other solution is $\mathbf{x}^++t\mathbf{v}_3$ with $\lVert\mathbf{x}\rVert^2=2+t^2\cdot1$ (since $\mathbf{x}^+\perp\mathbf{v}_3$), so $\mathbf{x}^+$ with norm $\sqrt2$ is the shortest. In an under-determined estimation problem (more unknowns than measurements), this is the "least-energy" estimate.

> [!example] Worked Example 10.29 — Condition number and sensitivity
> Solve $\mathbf{A}\mathbf{x}=\mathbf{b}$ with $\mathbf{A}=\begin{bmatrix}1&1\\1&1.0001\end{bmatrix}$ for $\mathbf{b}=(2,\,2.0001)$ and for the slightly perturbed $\mathbf{b}'=(2,\,2.0002)$.
>
> **Solutions.** Subtracting the equations: $0.0001x_2=b_2-b_1$. For $\mathbf{b}$: $x_2=1$, $x_1=1$. For $\mathbf{b}'$: $x_2=2$, $x_1=0$.
>
> **Sizes.** $\lVert\delta\mathbf{b}\rVert/\lVert\mathbf{b}\rVert=10^{-4}/2.83=3.5\times10^{-5}$, while $\lVert\delta\mathbf{x}\rVert/\lVert\mathbf{x}\rVert=\sqrt2/\sqrt2=1$: an amplification of about $28\,000$.
>
> **Explanation by SVD.** $\mathbf{A}$ is symmetric with eigenvalues $\lambda=\tfrac12\big(2.0001\pm\sqrt{2.0001^2-4\times0.0001}\big)$, i.e. $\sigma_1=2.00005$, $\sigma_2=0.00005$, so $\kappa_2=40\,002$. The observed amplification ($28\,000$) is below the bound $\kappa$, as it must be. The perturbation $\delta\mathbf{b}=(0,10^{-4})$ has a large component along the weak direction $\mathbf{u}_2\approx\frac1{\sqrt2}(1,-1)$, which $\mathbf{A}^{-1}$ amplifies by $1/\sigma_2=20\,000$. Two sensors measuring almost the same combination of two unknowns produce exactly this situation (compare the near-dependence note in [Ch. 9](09-vector-spaces-and-linear-transformations.md)).

> [!example] Worked Example 10.30 — MIMO channel as parallel sub-channels
> A $2\times2$ narrowband MIMO link has $\mathbf{y}=\mathbf{H}\mathbf{s}+\mathbf{n}$, and suppose $\mathbf{H}=\begin{bmatrix}3&0\\4&5\end{bmatrix}$ (in a real channel the entries are complex gains; the same steps apply with $\mathsf{H}$ in place of $\mathsf{T}$). With $\mathbf{H}=\mathbf{U}\boldsymbol{\Sigma}\mathbf{V}^\mathsf{T}$ from Example 10.27, let the transmitter send $\mathbf{s}=\mathbf{V}\tilde{\mathbf{s}}$ (precoding) and the receiver compute $\tilde{\mathbf{y}}=\mathbf{U}^\mathsf{T}\mathbf{y}$:
> $$
> \tilde{\mathbf{y}}=\mathbf{U}^\mathsf{T}\mathbf{U}\boldsymbol{\Sigma}\mathbf{V}^\mathsf{T}\mathbf{V}\tilde{\mathbf{s}}+\mathbf{U}^\mathsf{T}\mathbf{n}=\boldsymbol{\Sigma}\tilde{\mathbf{s}}+\tilde{\mathbf{n}},
> \qquad\text{i.e.}\qquad
> \tilde y_1=6.708\,\tilde s_1+\tilde n_1,\quad \tilde y_2=2.236\,\tilde s_2+\tilde n_2 .
> $$
> Because $\mathbf{U}$ is orthogonal, $\tilde{\mathbf{n}}$ has the same white statistics as $\mathbf{n}$. The coupled $2\times2$ channel has become **two independent scalar channels** with power gains $\sigma_1^2=45$ and $\sigma_2^2=5$ (16.5 dB and 7.0 dB). Allocating transmit power across them ("water-filling") is then a one-variable problem per sub-channel. The number of nonzero singular values is the number of spatial streams the channel can carry; a rank-deficient (keyhole) channel has only one.

> [!warning] Common Pitfall — confusing singular values with eigenvalues, and sign mismatches
> Singular values are always real and non-negative; eigenvalues need not be either. "The matrix is nearly singular because an eigenvalue is small" is correct, but "the matrix is well conditioned because its eigenvalues are similar in size" is **not** — $\begin{bmatrix}1&1000\\0&1\end{bmatrix}$ has both eigenvalues equal to $1$ but $\kappa_2\approx10^6$. And when computing an SVD by hand, find one set of singular vectors from a Gram matrix and derive the other from $\mathbf{A}\mathbf{v}_i=\sigma_i\mathbf{u}_i$; computing both sets independently risks inconsistent signs and a factorization that does not reproduce $\mathbf{A}$.

> [!tip] Engineering Note — the SVD in practice
> SVD algorithms (Golub–Kahan bidiagonalization followed by an implicit QR iteration) never form $\mathbf{A}^\mathsf{T}\mathbf{A}$, which would square the condition number and lose half the significant digits. Typical engineering uses: computing rank and null spaces robustly (e.g. finding the number of independent loops or the controllable subspace of a large model); least-squares fitting when the normal equations are ill-conditioned; direction-of-arrival estimation (MUSIC and ESPRIT split the SVD of the array data matrix into signal and noise subspaces); principal component analysis; and reduced-order modelling of interconnect and power-delivery networks.

## 10.8 Modal analysis

### Motivation

Couple two identical pendulums with a weak spring, start one swinging and leave the other at rest: the motion slowly transfers to the second pendulum and back again. Couple two identical LC resonators and the single resonance peak splits into two. Excite one conductor of a multi-conductor transmission line and signals appear on its neighbours with two different propagation speeds. Each of these puzzling behaviours becomes simple once the system is viewed in its **normal modes** — the eigenvectors of the system matrices. In a normal mode every part of the system moves at the same frequency with fixed relative amplitudes, as if the system were a single oscillator. Any motion is a superposition of modes, each evolving independently. **Modal analysis** is the systematic use of this fact: find the modes (an eigenvalue problem), project the initial conditions and forcing onto them, solve $n$ scalar problems, and superpose.

### Undamped second-order systems: the generalized eigenproblem

Lossless mechanical structures and lossless LC networks obey equations of the form

$$
\mathbf{M}\ddot{\mathbf{x}}+\mathbf{K}\mathbf{x}=\mathbf{0},
$$

where $\mathbf{M}$ is symmetric positive definite (mass, or a capacitance matrix) and $\mathbf{K}$ is symmetric positive semidefinite (stiffness, or a reciprocal-inductance matrix $\boldsymbol{\Gamma}=\mathbf{L}^{-1}$). Both symmetries come from energy: kinetic energy $\tfrac12\dot{\mathbf{x}}^\mathsf{T}\mathbf{M}\dot{\mathbf{x}}$ and potential energy $\tfrac12\mathbf{x}^\mathsf{T}\mathbf{K}\mathbf{x}$.

Seek synchronous solutions $\mathbf{x}(t)=\boldsymbol{\phi}\cos(\omega t+\theta)$, in which all coordinates oscillate in phase with a fixed shape $\boldsymbol{\phi}$. Substituting, $(-\omega^2\mathbf{M}\boldsymbol{\phi}+\mathbf{K}\boldsymbol{\phi})\cos(\omega t+\theta)=\mathbf{0}$ for all $t$, so

$$
\mathbf{K}\boldsymbol{\phi}=\omega^2\mathbf{M}\boldsymbol{\phi},
\qquad\text{with}\qquad
\det(\mathbf{K}-\omega^2\mathbf{M})=0 .
$$

This is a **generalized eigenvalue problem** with eigenvalue $\omega^2$; the determinant condition is the **frequency equation**. Multiplying by $\mathbf{M}^{-1}$ turns it into the standard problem $\mathbf{M}^{-1}\mathbf{K}\boldsymbol{\phi}=\omega^2\boldsymbol{\phi}$, but $\mathbf{M}^{-1}\mathbf{K}$ is generally not symmetric, so it is better to keep the symmetric pair.

> [!theorem] Normal modes of a conservative system
> Let $\mathbf{M}\succ0$ and $\mathbf{K}=\mathbf{K}^\mathsf{T}\succeq0$ be $n\times n$. Then the generalized problem $\mathbf{K}\boldsymbol{\phi}=\omega^2\mathbf{M}\boldsymbol{\phi}$ has $n$ real non-negative eigenvalues $\omega_1^2\le\dots\le\omega_n^2$ and real eigenvectors (mode shapes) that can be chosen **$\mathbf{M}$-orthonormal**:
> $$
> \boldsymbol{\phi}_i^\mathsf{T}\mathbf{M}\boldsymbol{\phi}_k=\delta_{ik},\qquad \boldsymbol{\phi}_i^\mathsf{T}\mathbf{K}\boldsymbol{\phi}_k=\omega_i^2\,\delta_{ik}.
> $$
> With the modal matrix $\boldsymbol{\Phi}=[\boldsymbol{\phi}_1\ \cdots\ \boldsymbol{\phi}_n]$: $\boldsymbol{\Phi}^\mathsf{T}\mathbf{M}\boldsymbol{\Phi}=\mathbf{I}$ and $\boldsymbol{\Phi}^\mathsf{T}\mathbf{K}\boldsymbol{\Phi}=\boldsymbol{\Omega}^2=\operatorname{diag}(\omega_i^2)$.

*Proof sketch.* Factor $\mathbf{M}=\mathbf{L}\mathbf{L}^\mathsf{T}$ (Cholesky) and put $\boldsymbol{\psi}=\mathbf{L}^\mathsf{T}\boldsymbol{\phi}$. Then $\mathbf{K}\boldsymbol{\phi}=\omega^2\mathbf{M}\boldsymbol{\phi}$ becomes $\tilde{\mathbf{K}}\boldsymbol{\psi}=\omega^2\boldsymbol{\psi}$ with $\tilde{\mathbf{K}}=\mathbf{L}^{-1}\mathbf{K}\mathbf{L}^{-\mathsf{T}}$, which is symmetric PSD. The spectral theorem (Section 10.5) gives real $\omega_i^2\ge0$ and orthonormal $\boldsymbol{\psi}_i$; transforming back, $\boldsymbol{\psi}_i^\mathsf{T}\boldsymbol{\psi}_k=\boldsymbol{\phi}_i^\mathsf{T}\mathbf{L}\mathbf{L}^\mathsf{T}\boldsymbol{\phi}_k=\boldsymbol{\phi}_i^\mathsf{T}\mathbf{M}\boldsymbol{\phi}_k$. $\square$

$\mathbf{M}$-orthogonality is orthogonality in the kinetic-energy inner product (Section 10.6): modes do not share energy. A zero eigenvalue $\omega^2=0$ is a **rigid-body mode** (or, in a circuit, a DC loop current or floating node voltage that the network does not restore).

### Modal coordinates and the general solution

Substitute $\mathbf{x}=\boldsymbol{\Phi}\mathbf{q}$ into $\mathbf{M}\ddot{\mathbf{x}}+\mathbf{K}\mathbf{x}=\mathbf{f}(t)$ (now with a forcing term) and multiply on the left by $\boldsymbol{\Phi}^\mathsf{T}$:

$$
\underbrace{\boldsymbol{\Phi}^\mathsf{T}\mathbf{M}\boldsymbol{\Phi}}_{\mathbf{I}}\ddot{\mathbf{q}}+\underbrace{\boldsymbol{\Phi}^\mathsf{T}\mathbf{K}\boldsymbol{\Phi}}_{\boldsymbol{\Omega}^2}\mathbf{q}=\boldsymbol{\Phi}^\mathsf{T}\mathbf{f}(t)
\quad\Longleftrightarrow\quad
\ddot q_i+\omega_i^2q_i=\boldsymbol{\phi}_i^\mathsf{T}\mathbf{f}(t),\quad i=1,\dots,n.
$$

The **modal coordinates** $q_i$ obey uncoupled single-degree-of-freedom equations. For free vibration,

$$
\mathbf{x}(t)=\sum_{i=1}^n\boldsymbol{\phi}_i\big(a_i\cos\omega_it+b_i\sin\omega_it\big),\qquad
a_i=\boldsymbol{\phi}_i^\mathsf{T}\mathbf{M}\mathbf{x}(0),\quad b_i=\frac{\boldsymbol{\phi}_i^\mathsf{T}\mathbf{M}\dot{\mathbf{x}}(0)}{\omega_i},
$$

where the formulas for $a_i,b_i$ use $\mathbf{M}$-orthonormality to invert $\mathbf{x}(0)=\boldsymbol{\Phi}\mathbf{a}$ without solving a linear system: $\mathbf{a}=\boldsymbol{\Phi}^{-1}\mathbf{x}(0)=\boldsymbol{\Phi}^\mathsf{T}\mathbf{M}\mathbf{x}(0)$.

The forced equations reveal two key facts. **Resonance**: a sinusoidal force at frequency $\omega$ drives mode $i$ with amplitude proportional to $\boldsymbol{\phi}_i^\mathsf{T}\mathbf{f}/(\omega_i^2-\omega^2)$, which blows up (in the undamped model) as $\omega\to\omega_i$. **Selective excitation**: if $\boldsymbol{\phi}_i^\mathsf{T}\mathbf{f}=0$ — for instance a symmetric force on a symmetric structure and an antisymmetric mode — that mode is not excited at all, regardless of frequency. Likewise a sensor placed at a node of a mode shape (where $\phi_{ik}=0$) cannot see that mode.

### Damping

With a damping matrix, $\mathbf{M}\ddot{\mathbf{x}}+\mathbf{C}\dot{\mathbf{x}}+\mathbf{K}\mathbf{x}=\mathbf{f}$, the transformation $\boldsymbol{\Phi}^\mathsf{T}\mathbf{C}\boldsymbol{\Phi}$ is in general not diagonal and the modes couple. Two cases are standard:

- **Proportional (Rayleigh) damping**, $\mathbf{C}=\alpha\mathbf{M}+\beta\mathbf{K}$: then $\boldsymbol{\Phi}^\mathsf{T}\mathbf{C}\boldsymbol{\Phi}=\alpha\mathbf{I}+\beta\boldsymbol{\Omega}^2$ is diagonal, and each mode obeys $\ddot q_i+2\zeta_i\omega_i\dot q_i+\omega_i^2q_i=\boldsymbol{\phi}_i^\mathsf{T}\mathbf{f}$ with modal damping ratio $\zeta_i=\frac{\alpha}{2\omega_i}+\frac{\beta\omega_i}{2}$. The same mode shapes survive; each decays at its own rate.
- **General damping**: convert to first-order state form with $\mathbf{z}=(\mathbf{x},\dot{\mathbf{x}})$, $\dot{\mathbf{z}}=\begin{bmatrix}\mathbf{0}&\mathbf{I}\\-\mathbf{M}^{-1}\mathbf{K}&-\mathbf{M}^{-1}\mathbf{C}\end{bmatrix}\mathbf{z}$, and apply the first-order modal analysis below. The eigenvalues come in pairs $-\zeta_i\omega_i\pm j\omega_i\sqrt{1-\zeta_i^2}$, and the mode shapes become complex (different parts of the structure have different phases).

### First-order systems: modal decomposition with left and right eigenvectors

For $\dot{\mathbf{x}}=\mathbf{A}\mathbf{x}+\mathbf{B}\mathbf{u}$ with diagonalizable $\mathbf{A}=\mathbf{P}\boldsymbol{\Lambda}\mathbf{P}^{-1}$, let $\mathbf{v}_i$ be the columns of $\mathbf{P}$ (right eigenvectors) and $\mathbf{w}_i^\mathsf{T}$ the rows of $\mathbf{W}=\mathbf{P}^{-1}$ (left eigenvectors, normalised so that $\mathbf{w}_i^\mathsf{T}\mathbf{v}_k=\delta_{ik}$). Then the zero-input response is

$$
\mathbf{x}(t)=e^{\mathbf{A}t}\mathbf{x}(0)=\sum_{i=1}^n\mathbf{v}_i\,\big(\mathbf{w}_i^\mathsf{T}\mathbf{x}(0)\big)\,e^{\lambda_it}.
$$

Each term has a clear reading: $\mathbf{w}_i^\mathsf{T}\mathbf{x}(0)$ is **how much the initial state excites** mode $i$, $e^{\lambda_it}$ is **how the mode evolves**, and $\mathbf{v}_i$ is **how the mode appears** in the state variables. With input and output $\mathbf{y}=\mathbf{C}\mathbf{x}$, mode $i$ contributes $\dfrac{(\mathbf{C}\mathbf{v}_i)(\mathbf{w}_i^\mathsf{T}\mathbf{B})}{s-\lambda_i}$ to the transfer matrix — the partial-fraction expansion of $\mathbf{C}(s\mathbf{I}-\mathbf{A})^{-1}\mathbf{B}$ ([Ch. 14](../part-5-transforms/14-laplace-transforms.md)), with residues formed from eigenvectors. If $\mathbf{w}_i^\mathsf{T}\mathbf{B}=\mathbf{0}$ the mode is **uncontrollable** from the input; if $\mathbf{C}\mathbf{v}_i=\mathbf{0}$ it is **unobservable** at the output; either way it cancels from the transfer function. In power-system small-signal stability, the products $p_{ki}=v_{ki}w_{ik}$ are the **participation factors** that tell which generator state participates in which electromechanical oscillation mode.

### Worked examples

> [!example] Worked Example 10.31 — Two masses, three springs
> Two equal masses $m=1$ kg are connected wall–spring–mass–spring–mass–spring–wall, all springs $k=100$ N/m. Find the natural frequencies and mode shapes.
>
> **Step 1 — equations.** Mass 1 feels $-kx_1$ from the left spring and $-k(x_1-x_2)$ from the middle one; similarly for mass 2:
> $$
> m\ddot{\mathbf{x}}+\mathbf{K}\mathbf{x}=\mathbf{0},\qquad \mathbf{K}=k\begin{bmatrix}2&-1\\-1&2\end{bmatrix}.
> $$
> **Step 2 — frequency equation.** With $\mathbf{M}=m\mathbf{I}$, $\det(\mathbf{K}-\omega^2m\mathbf{I})=0$ gives $(2k-m\omega^2)^2=k^2$, so $2k-m\omega^2=\pm k$:
> $$
> \omega_1^2=\frac{k}{m},\qquad \omega_2^2=\frac{3k}{m}.
> $$
> Numerically $\omega_1=10$ rad/s ($1.59$ Hz) and $\omega_2=17.32$ rad/s ($2.76$ Hz).
>
> **Step 3 — mode shapes.** $\omega_1$: $(\mathbf{K}-k\mathbf{I})\boldsymbol{\phi}=k\begin{bmatrix}1&-1\\-1&1\end{bmatrix}\boldsymbol{\phi}=\mathbf{0}$, $\boldsymbol{\phi}_1\propto(1,1)$: **in phase** — the middle spring is never stretched, so each mass sees only one spring, $\omega^2=k/m$. $\omega_2$: $\boldsymbol{\phi}_2\propto(1,-1)$: **antiphase** — the middle spring's midpoint is stationary, so each mass sees $k$ plus a half-length spring of stiffness $2k$, $\omega^2=3k/m$.
>
> **Step 4 — $\mathbf{M}$-normalise.** $\boldsymbol{\phi}_i^\mathsf{T}m\boldsymbol{\phi}_i=1$ gives $\boldsymbol{\phi}_1=\frac{1}{\sqrt{2m}}(1,1)$, $\boldsymbol{\phi}_2=\frac1{\sqrt{2m}}(1,-1)$. **Check:** $\boldsymbol{\phi}_2^\mathsf{T}\mathbf{K}\boldsymbol{\phi}_2=\frac{k}{2m}(2+2+2)=3k/m=\omega_2^2$ $\checkmark$. The physical reasoning in Step 3 — recognising which springs carry force in each mode — is a valuable check on any computed mode.

> [!example] Worked Example 10.32 — Two capacitively coupled LC resonators
> Two identical parallel tanks ($L=10\ \mu$H, $C=1$ nF) have their top nodes joined by a coupling capacitor $C_c=0.1$ nF. Find the resonant frequencies and the mode shapes.
>
> **Step 1 — node equations.** KCL at node 1: $C\dot v_1+C_c(\dot v_1-\dot v_2)+i_{L1}=0$ with $\dot i_{L1}=v_1/L$. Differentiating once to eliminate $i_{L1}$, and similarly at node 2:
> $$
> \mathbf{C}_n\ddot{\mathbf{v}}+\frac1L\mathbf{v}=\mathbf{0},\qquad \mathbf{C}_n=\begin{bmatrix}C+C_c&-C_c\\-C_c&C+C_c\end{bmatrix}.
> $$
> This is $\mathbf{M}\ddot{\mathbf{x}}+\mathbf{K}\mathbf{x}=\mathbf{0}$ with the **capacitance matrix as mass** and the **reciprocal inductance as stiffness**.
>
> **Step 2 — modes by symmetry.** Both matrices commute with the swap of nodes 1 and 2, so the even and odd vectors are eigenvectors of both (Section 10.3):
> - **Even mode** $(1,1)$: $\mathbf{C}_n(1,1)=C(1,1)$; no current in $C_c$. $\omega_e^2=\dfrac{1}{LC}$.
> - **Odd mode** $(1,-1)$: $\mathbf{C}_n(1,-1)=(C+2C_c)(1,-1)$; the midpoint of $C_c$ is a virtual ground, so each tank sees an extra $2C_c$. $\omega_o^2=\dfrac{1}{L(C+2C_c)}$.
>
> **Step 3 — numbers.** $f_e=\dfrac{1}{2\pi\sqrt{10^{-5}\cdot10^{-9}}}=1.5915$ MHz; $f_o=\dfrac{1}{2\pi\sqrt{10^{-5}\cdot1.2\times10^{-9}}}=1.4529$ MHz. The single resonance of an isolated tank (1.5915 MHz) has **split** into two, separated by $139$ kHz.
>
> **Step 4 — coupling coefficient.** The standard filter-design measure is $k=\dfrac{f_e^2-f_o^2}{f_e^2+f_o^2}=\dfrac{1-1/1.2}{1+1/1.2}=0.0909=\dfrac{C_c}{C+C_c}$ $\checkmark$. Coupled-resonator bandpass filters are designed by choosing these mode splittings; the even/odd decomposition is exactly the eigen-decomposition of the symmetric pair $(\mathbf{K},\mathbf{M})$.

> [!example] Worked Example 10.33 — Beating and energy exchange between weakly coupled oscillators
> Two unit masses on springs $k=1$ N/m, coupled by a weak spring $k_c=0.105$ N/m, have $\mathbf{K}=\begin{bmatrix}1+k_c&-k_c\\-k_c&1+k_c\end{bmatrix}$. Mass 1 is displaced by $1$ m and released; mass 2 starts at rest. Find $\mathbf{x}(t)$.
>
> **Step 1 — modes.** As in Example 10.31: in-phase $(1,1)$ with $\omega_1^2=1$, so $\omega_1=1$ rad/s; antiphase $(1,-1)$ with $\omega_2^2=1+2k_c=1.21$, so $\omega_2=1.1$ rad/s.
>
> **Step 2 — modal amplitudes.** $\mathbf{x}(0)=(1,0)=\tfrac12(1,1)+\tfrac12(1,-1)$ and $\dot{\mathbf{x}}(0)=\mathbf{0}$, so only cosine terms appear, each with amplitude $\tfrac12$.
>
> **Step 3 — superpose and simplify** with sum-to-product identities ([Ch. 3](../part-1-foundations/03-trigonometry-and-periodic-functions.md)):
> $$
> \begin{aligned}
> x_1(t)&=\tfrac12(\cos t+\cos1.1t)=\cos(0.05t)\cos(1.05t),\\
> x_2(t)&=\tfrac12(\cos t-\cos1.1t)=\sin(0.05t)\sin(1.05t).
> \end{aligned}
> $$
> **Interpretation.** Both masses oscillate at the mean frequency $1.05$ rad/s under slowly varying envelopes $\lvert\cos0.05t\rvert$ and $\lvert\sin0.05t\rvert$. At $t=\pi/0.1=31.4$ s mass 1 is momentarily still and mass 2 carries all the motion; at $62.8$ s it has all come back. The exchange time $\pi/(\omega_2-\omega_1)$ is set by the **mode splitting**: the weaker the coupling, the slower the exchange. The same mathematics describes energy transfer in a directional coupler and between coupled optical waveguides, with distance in place of time.

```desmos-graph
left=0; right=130; bottom=-2.5; top=2.5
---
y=1.2+0.5(\cos(x)+\cos(1.1x))|blue
y=1.2+\cos(0.05x)|blue|dashed
y=1.2-\cos(0.05x)|blue|dashed
y=-1.2+0.5(\cos(x)-\cos(1.1x))|red
y=-1.2+\sin(0.05x)|red|dashed
y=-1.2-\sin(0.05x)|red|dashed
```
*Free motion of two weakly coupled oscillators (Example 10.33), plotted against $t$ in seconds: $x_1(t)$ (blue, offset $+1.2$) and $x_2(t)$ (red, offset $-1.2$). The dashed envelopes $\pm\cos0.05t$ and $\pm\sin0.05t$ show energy passing completely from one oscillator to the other every $31.4$ s — the beat between the normal-mode frequencies $1$ and $1.1$ rad/s.*

> [!example] Worked Example 10.34 — Modes of an RC ladder (first-order modal analysis)
> A two-section RC ladder ($R=1$ k$\Omega$, $C=1\ \mu$F per section, so $RC=1$ ms) has its input shorted and both capacitors initially at $1$ V. Find $v_1(t)$ and $v_2(t)$.
>
> **Step 1 — state equations.** KCL at node 1: $C\dot v_1=-v_1/R-(v_1-v_2)/R$; at node 2: $C\dot v_2=(v_1-v_2)/R$. With $t$ in ms,
> $$
> \dot{\mathbf{v}}=\mathbf{A}\mathbf{v},\qquad \mathbf{A}=\begin{bmatrix}-2&1\\1&-1\end{bmatrix}.
> $$
> **Step 2 — eigenvalues.** $\lambda^2+3\lambda+1=0$: $\lambda=\dfrac{-3\pm\sqrt5}{2}$, i.e. $\lambda_s=-0.382$ and $\lambda_f=-2.618$ ms$^{-1}$. Time constants: $\tau_s=2.618$ ms (slow) and $\tau_f=0.382$ ms (fast).
>
> **Step 3 — eigenvectors.** $\mathbf{A}$ is symmetric, so they are orthogonal. From the second row, $v_1=(1+\lambda)v_2$: $\mathbf{v}_s\propto(0.618,1)\propto(1,1.618)$ (both capacitors discharging together, the far one higher) and $\mathbf{v}_f\propto(1,-0.618)$ (charge sloshing between the capacitors). Normalised: $\mathbf{q}_s=(0.5257,0.8507)$, $\mathbf{q}_f=(0.8507,-0.5257)$.
>
> **Step 4 — project the initial state.** Orthonormal eigenvectors make the left eigenvectors equal to the right ones: $c_s=\mathbf{q}_s\cdot(1,1)=1.3764$, $c_f=\mathbf{q}_f\cdot(1,1)=0.3249$.
>
> **Step 5 — response.**
> $$
> \begin{aligned}
> v_1(t)&=0.7236\,e^{-t/2.618}+0.2764\,e^{-t/0.382},\\
> v_2(t)&=1.1708\,e^{-t/2.618}-0.1708\,e^{-t/0.382}\qquad(t\text{ in ms, }v\text{ in V}).
> \end{aligned}
> $$
> **Checks.** $t=0$: $v_1=1.0000$, $v_2=1.0000$ $\checkmark$. At $t=1$ ms: $v_1=0.7236(0.6825)+0.2764(0.0729)=0.5140$ V and $v_2=0.7866$ V, matching a numerical $e^{\mathbf{A}t}\mathbf{v}(0)$ $\checkmark$. After about $1$ ms the fast mode has essentially vanished and the circuit behaves as a single-time-constant system with $\tau=2.618$ ms — the basis of **dominant-pole approximation**. Note that neither time constant equals $RC$: the sections interact, which is why "$n$ cascaded RC sections" is not $n$ identical poles.

> [!warning] Common Pitfall — normalising modes and projecting initial conditions
> When the mass (or capacitance) matrix is not a multiple of $\mathbf{I}$, mode shapes are orthogonal **with respect to $\mathbf{M}$**, not in the ordinary dot product. The modal amplitude is $a_i=\boldsymbol{\phi}_i^\mathsf{T}\mathbf{M}\mathbf{x}(0)$ for $\mathbf{M}$-normalised modes; using $\boldsymbol{\phi}_i^\mathsf{T}\mathbf{x}(0)$ gives wrong answers. For non-symmetric first-order systems, use the **left** eigenvectors (rows of $\mathbf{P}^{-1}$) to project; the right eigenvectors are not orthogonal.

> [!tip] Engineering Note — modal analysis at scale
> Finite-element models of PCBs, packages, enclosures and antennas have $10^4$–$10^7$ degrees of freedom, but only the few modes in the frequency band of interest matter. Solvers (Lanczos and Arnoldi methods, shift-invert about a target frequency) extract just those eigenpairs from the sparse pair $(\mathbf{K},\mathbf{M})$. Cavity resonances of a metal enclosure, which can turn a shielded box into an efficient radiator at specific frequencies, are found this way; so are the even and odd propagation modes of coupled microstrip lines (eigenvectors of the product of per-unit-length $\mathbf{L}$ and $\mathbf{C}$ matrices), whose different velocities cause far-end crosstalk.

## 10.9 Applications to system dynamics and stability

### Motivation

Every circuit, control loop and oscillator must answer one question before any other: **if it is disturbed, does it return to rest, settle into a steady oscillation, or run away?** An amplifier that oscillates, a phase-locked loop that loses lock, a switched-mode converter whose output voltage rings ever larger, a power grid whose generators swing apart — all are stability failures, and for linear (or linearized) models all are diagnosed by eigenvalues. This section collects the results: stability criteria for continuous- and discrete-time systems, the complete geometric classification of two-dimensional systems, linearization of nonlinear systems about an equilibrium, and Lyapunov's energy-like test, which brings positive-definite matrices back into the story. Chapters [11](../part-4-differential-equations/11-ordinary-differential-equations.md) and [12](../part-4-differential-equations/12-linear-differential-equations.md) develop the differential-equation side; here we concentrate on what the matrix tells us.

### Stability of linear systems

Consider $\dot{\mathbf{x}}=\mathbf{A}\mathbf{x}$, whose equilibrium is $\mathbf{x}=\mathbf{0}$. From Section 10.4 every solution is a combination of terms $t^ke^{\lambda_it}\mathbf{v}$, with polynomial factors $t^k$ only for defective eigenvalues. Since $\lvert e^{\lambda t}\rvert=e^{(\operatorname{Re}\lambda)t}$, the real parts decide everything.

> [!definition] Stability of an equilibrium
> The equilibrium $\mathbf{x}=\mathbf{0}$ of $\dot{\mathbf{x}}=\mathbf{A}\mathbf{x}$ is
> - **asymptotically stable** if every solution tends to $\mathbf{0}$ as $t\to\infty$;
> - **marginally stable** (stable but not asymptotically) if every solution stays bounded but some do not tend to $\mathbf{0}$;
> - **unstable** if some solution grows without bound.

> [!theorem] Eigenvalue criteria for stability
> **Continuous time**, $\dot{\mathbf{x}}=\mathbf{A}\mathbf{x}$:
> - asymptotically stable $\iff$ $\operatorname{Re}\lambda_i<0$ for every eigenvalue (all eigenvalues in the open left half-plane; $\mathbf{A}$ is then called **Hurwitz**);
> - marginally stable $\iff$ $\operatorname{Re}\lambda_i\le0$ for all $i$, and every eigenvalue on the imaginary axis is non-defective;
> - unstable otherwise.
>
> **Discrete time**, $\mathbf{x}[k+1]=\mathbf{A}\mathbf{x}[k]$: the same statements with "$\operatorname{Re}\lambda<0$" replaced by "$\lvert\lambda\rvert<1$" and "the imaginary axis" by "the unit circle" ($\mathbf{A}$ is then called **Schur**).

For a linear time-invariant system with input and output, the eigenvalues of $\mathbf{A}$ are the **poles** of the transfer function $\mathbf{C}(s\mathbf{I}-\mathbf{A})^{-1}\mathbf{B}+\mathbf{D}$ — except for any modes that are uncontrollable or unobservable and therefore cancel (Section 10.8). This is why "all poles in the left half-plane" is the same condition as "all eigenvalues in the left half-plane" for a minimal realisation, and why a hidden unstable mode is dangerous: it does not appear in the transfer function, yet the internal state still grows.

### Classifying two-dimensional systems: the trace–determinant plane

For a real $2\times2$ matrix, write $\tau=\operatorname{tr}\mathbf{A}$ and $\Delta=\det\mathbf{A}$. The eigenvalues are

$$
\lambda_{1,2}=\frac{\tau\pm\sqrt{\tau^2-4\Delta}}{2},\qquad \lambda_1+\lambda_2=\tau,\quad\lambda_1\lambda_2=\Delta,
$$

so the two numbers $(\tau,\Delta)$ determine the qualitative behaviour completely:

| Region of the $(\tau,\Delta)$ plane | Eigenvalues | Equilibrium type |
|---|---|---|
| $\Delta<0$ | real, opposite signs | **saddle** (unstable) |
| $\Delta>0$, $\tau^2>4\Delta$, $\tau<0$ | real, both negative | **stable node** |
| $\Delta>0$, $\tau^2>4\Delta$, $\tau>0$ | real, both positive | **unstable node** |
| $\Delta>0$, $\tau^2<4\Delta$, $\tau<0$ | complex, $\operatorname{Re}<0$ | **stable spiral** (focus) |
| $\Delta>0$, $\tau^2<4\Delta$, $\tau>0$ | complex, $\operatorname{Re}>0$ | **unstable spiral** |
| $\Delta>0$, $\tau=0$ | $\pm j\sqrt\Delta$ | **centre** (marginally stable) |
| $\tau^2=4\Delta$ | real, repeated | degenerate (star) node — boundary case |
| $\Delta=0$ | one zero eigenvalue | line of equilibria — boundary case |

The compact stability statement for $n=2$ is: **asymptotically stable $\iff\tau<0$ and $\Delta>0$**. For a second-order circuit with characteristic polynomial $\lambda^2+2\zeta\omega_0\lambda+\omega_0^2$, $\tau=-2\zeta\omega_0$ and $\Delta=\omega_0^2$, and the parabola $\tau^2=4\Delta$ is the critical damping $\zeta=1$: overdamped circuits are stable nodes, underdamped ones stable spirals, and a lossless LC tank is a centre.

```desmos-graph
left=-4; right=4; bottom=-2; top=5
---
y=\frac{x^2}{4}|red
y=0|black
x=0|black|0<y<5
(-3,1.2)|label:stable node
(3,1.2)|label:unstable node
(-1.2,3.5)|label:stable spiral
(1.2,3.5)|label:unstable spiral
(0,-1)|label:saddle
(0,2.5)|label:centre (τ=0)
```
*The trace–determinant plane (horizontal axis $\tau=\operatorname{tr}\mathbf{A}$, vertical axis $\Delta=\det\mathbf{A}$). Below the $\tau$-axis lie saddles; above it, the parabola $\Delta=\tau^2/4$ (red, critical damping) separates nodes (outside) from spirals (inside); the positive $\Delta$-axis carries centres. The stable region is the open upper-left quadrant $\tau<0,\ \Delta>0$.*

### Phase portraits from eigenvectors

The **phase portrait** draws solution trajectories in the $(x_1,x_2)$ plane. Eigenvectors organise it: a trajectory that starts on a real eigenvector line stays on it (it is the solution $c\,e^{\lambda t}\mathbf{v}$), moving inward if $\lambda<0$ and outward if $\lambda>0$. In eigenvector coordinates $(y_1,y_2)$ a system with real eigenvalues is $\dot y_1=\lambda_1y_1$, $\dot y_2=\lambda_2y_2$, so $y_i=c_ie^{\lambda_it}$, and eliminating $t$ gives the **explicit trajectory curves**

$$
y_2=C\,\lvert y_1\rvert^{\lambda_2/\lambda_1}.
$$

For a stable node with $\lambda_1=-1$, $\lambda_2=-2$: $y_2=Cy_1^2$, a family of parabolas all tangent at the origin to the **slow** eigenvector ($y_1$-axis). For a saddle with $\lambda_1=1$, $\lambda_2=-1$: $y_1y_2=C$, hyperbolas whose asymptotes are the eigenvector lines. For a complex pair $\sigma\pm j\omega$, the real block form $\begin{bmatrix}\sigma&-\omega\\\omega&\sigma\end{bmatrix}$ in polar coordinates gives $\dot r=\sigma r$, $\dot\theta=\omega$, hence $r=r_0e^{(\sigma/\omega)(\theta-\theta_0)}$: **logarithmic spirals**, shrinking if $\sigma<0$. In the original coordinates these pictures are mapped through $\mathbf{P}$, which skews them but preserves their type.

```desmos-graph
left=-3; right=3; bottom=-3; top=3
---
y=0.5x^2|blue
y=-0.5x^2|blue
y=2x^2|blue
y=-2x^2|blue
y=0|blue
x=0|green|dashed
xy=1|red
xy=-1|red
xy=0.25|red
xy=-0.25|red
```
*Explicit trajectory families in eigenvector coordinates. Blue: stable node with $\lambda_1=-1$, $\lambda_2=-2$, trajectories $y_2=Cy_1^2$, all entering the origin tangent to the slow eigen-direction (the $y_1$-axis); the fast direction is the $y_2$-axis (green dashed). Red: saddle with $\lambda_1=1$, $\lambda_2=-1$, trajectories $y_1y_2=C$, whose asymptotes are the two eigenvector axes.*

```desmos-graph
left=-1.5; right=1.5; bottom=-1.5; top=1.5
---
r=e^{-0.2\theta}|blue|0<\theta<30
r=0.6e^{-0.2\theta}|purple|0<\theta<30
r=1.3|grey|dashed
```
*Trajectories of a stable spiral with eigenvalues $-0.2\pm j$ (real block form): $r=r_0e^{-0.2\theta}$, two starting radii shown. Each turn shrinks the radius by $e^{-0.2\cdot2\pi}=0.285$; with $\operatorname{Re}\lambda=0$ the curves would be circles (a centre, dashed reference), and with $\operatorname{Re}\lambda>0$ they would spiral outward.*

### Routh–Hurwitz: stability without eigenvalues

Often only the characteristic polynomial $\det(s\mathbf{I}-\mathbf{A})=s^n+a_{n-1}s^{n-1}+\dots+a_0$ is needed, and only the *signs* of the real parts of its roots. The Routh–Hurwitz criterion decides this from the coefficients. Its lowest-order cases are worth memorising:

- $n=1,2$: all roots in the left half-plane $\iff$ all coefficients positive.
- $n=3$: $s^3+a_2s^2+a_1s+a_0$ is Hurwitz $\iff a_2,a_1,a_0>0$ **and** $a_2a_1>a_0$.

Positivity of all coefficients is **necessary** for every $n$ (a product of factors $(s+p)$ and $(s^2+bs+c)$ with $p,b,c>0$ has only positive coefficients) but not sufficient for $n\ge3$: $s^3+s^2+s+6$ has positive coefficients but $a_2a_1=1<6$, and it has the unstable roots $0.5\pm j1.66$. The full Routh array for higher orders is part of control theory and is treated with the Laplace transform in [Ch. 14](../part-5-transforms/14-laplace-transforms.md).

### Linearization of nonlinear systems

Most real systems are nonlinear: $\dot{\mathbf{x}}=\mathbf{f}(\mathbf{x})$. An **equilibrium** $\mathbf{x}^*$ satisfies $\mathbf{f}(\mathbf{x}^*)=\mathbf{0}$ — the DC operating point of a circuit. For small deviations $\boldsymbol{\xi}=\mathbf{x}-\mathbf{x}^*$, the multivariable Taylor expansion of [Ch. 6](../part-2-calculus/06-multivariable-calculus.md) gives

$$
\dot{\boldsymbol{\xi}}=\mathbf{f}(\mathbf{x}^*+\boldsymbol{\xi})\approx\mathbf{J}\boldsymbol{\xi},\qquad
\mathbf{J}=\left.\frac{\partial\mathbf{f}}{\partial\mathbf{x}}\right|_{\mathbf{x}^*}=\left[\frac{\partial f_i}{\partial x_k}\right]_{\mathbf{x}^*},
$$

the **Jacobian matrix**. This is exactly the small-signal model of circuit analysis: transistors and diodes are replaced by their incremental conductances at the operating point.

> [!theorem] Linearized stability (Lyapunov's indirect method; Hartman–Grobman)
> Let $\mathbf{J}$ be the Jacobian at an equilibrium $\mathbf{x}^*$ of $\dot{\mathbf{x}}=\mathbf{f}(\mathbf{x})$, with $\mathbf{f}$ continuously differentiable.
> - If every eigenvalue of $\mathbf{J}$ has $\operatorname{Re}\lambda<0$, $\mathbf{x}^*$ is locally asymptotically stable.
> - If some eigenvalue has $\operatorname{Re}\lambda>0$, $\mathbf{x}^*$ is unstable.
> - If no eigenvalue has zero real part (a *hyperbolic* equilibrium), the nonlinear phase portrait near $\mathbf{x}^*$ has the same type (node, saddle, spiral) as the linearization.
>
> If some eigenvalue has $\operatorname{Re}\lambda=0$ and none is positive, the linearization is **inconclusive**: the neglected nonlinear terms decide.

### Lyapunov's equation: stability through energy

A physical way to prove stability is to exhibit an "energy" that always decreases. For a linear system try a quadratic form $V(\mathbf{x})=\mathbf{x}^\mathsf{T}\mathbf{P}\mathbf{x}$ with $\mathbf{P}\succ0$. Along solutions,

$$
\dot V=\dot{\mathbf{x}}^\mathsf{T}\mathbf{P}\mathbf{x}+\mathbf{x}^\mathsf{T}\mathbf{P}\dot{\mathbf{x}}=\mathbf{x}^\mathsf{T}(\mathbf{A}^\mathsf{T}\mathbf{P}+\mathbf{P}\mathbf{A})\mathbf{x}.
$$

> [!theorem] Lyapunov stability theorem (linear case)
> $\mathbf{A}$ is Hurwitz (all $\operatorname{Re}\lambda<0$) **if and only if** for some (equivalently, every) symmetric $\mathbf{Q}\succ0$ the **Lyapunov equation**
> $$
> \mathbf{A}^\mathsf{T}\mathbf{P}+\mathbf{P}\mathbf{A}=-\mathbf{Q}
> $$
> has a symmetric positive-definite solution $\mathbf{P}$. The solution is then unique and equals $\mathbf{P}=\int_0^\infty e^{\mathbf{A}^\mathsf{T}t}\mathbf{Q}e^{\mathbf{A}t}\,dt$.

*Why it works:* with such a $\mathbf{P}$, $V>0$ away from the origin while $\dot V=-\mathbf{x}^\mathsf{T}\mathbf{Q}\mathbf{x}<0$, so $V$ — and with it $\mathbf{x}$ — decays to zero; the level curves of $V$ are nested ellipses that trajectories cross only inwards. The discrete-time counterpart is $\mathbf{A}^\mathsf{T}\mathbf{P}\mathbf{A}-\mathbf{P}=-\mathbf{Q}$. The Lyapunov approach matters because it extends to nonlinear and time-varying systems, where eigenvalues alone say little, and because $\mathbf{P}$ itself has meanings: for a stable system driven by white noise, the solution of a closely related Lyapunov equation is the steady-state **covariance** of the state ([Ch. 21](../part-6-probability-statistics/21-random-processes.md)).

### Worked examples

> [!example] Worked Example 10.35 — Natural frequencies of a series RLC circuit
> A series RLC circuit with $R=100\ \Omega$, $L=10$ mH, $C=1\ \mu$F is left to discharge. With state $\mathbf{x}=(i,v_C)$, find the eigenvalues and characterise the response.
>
> **Step 1 — state equations.** KVL: $L\,di/dt=-Ri-v_C$; capacitor: $C\,dv_C/dt=i$.
> $$
> \dot{\mathbf{x}}=\mathbf{A}\mathbf{x},\qquad \mathbf{A}=\begin{bmatrix}-R/L&-1/L\\1/C&0\end{bmatrix}=\begin{bmatrix}-10^4&-100\\10^6&0\end{bmatrix}.
> $$
> (Units: $R/L$ in s$^{-1}$; $1/L$ in A/(V·s); $1/C$ in V/(A·s). The off-diagonal units differ, but their product $1/(LC)$ is in s$^{-2}$, as an entry of $\det\mathbf{A}$ must be.)
>
> **Step 2 — trace and determinant.** $\tau=-R/L=-10^4$ s$^{-1}$, $\Delta=1/(LC)=10^8$ s$^{-2}$. Since $\tau<0$ and $\Delta>0$: **asymptotically stable**. Since $\tau^2=10^8<4\Delta=4\times10^8$: **spiral** (underdamped).
>
> **Step 3 — eigenvalues.**
> $$
> \lambda=-5000\pm j\sqrt{10^8-2.5\times10^7}=-5000\pm j\,8660\ \text{s}^{-1}.
> $$
> **Interpretation.** $\omega_0=\sqrt\Delta=10^4$ rad/s, $\zeta=-\tau/(2\omega_0)=0.5$; the current rings at $\omega_d=8660$ rad/s ($1.378$ kHz) inside an envelope $e^{-t/0.2\,\text{ms}}$. Doubling $R$ to $200\ \Omega$ gives $\tau^2=4\Delta$: critical damping, a defective double eigenvalue at $-10^4$ s$^{-1}$, and $te^{-10^4t}$ terms (compare Example 10.11).

> [!example] Worked Example 10.36 — Reading the type off trace and determinant
> Classify the equilibrium of $\dot{\mathbf{x}}=\mathbf{A}\mathbf{x}$ for each matrix.
>
> | $\mathbf{A}$ | $\tau$ | $\Delta$ | $\tau^2-4\Delta$ | Eigenvalues | Type |
> |---|---|---|---|---|---|
> | $\begin{bmatrix}-3&1\\1&-3\end{bmatrix}$ | $-6$ | $8$ | $4$ | $-2,\,-4$ | stable node |
> | $\begin{bmatrix}1&2\\3&-4\end{bmatrix}$ | $-3$ | $-10$ | $49$ | $2,\,-5$ | saddle (unstable) |
> | $\begin{bmatrix}-1&-4\\1&-1\end{bmatrix}$ | $-2$ | $5$ | $-16$ | $-1\pm2j$ | stable spiral |
> | $\begin{bmatrix}0&1\\-4&0\end{bmatrix}$ | $0$ | $4$ | $-16$ | $\pm2j$ | centre |
> | $\begin{bmatrix}1&-2\\1&3\end{bmatrix}$ | $4$ | $5$ | $-4$ | $2\pm j$ | unstable spiral |
>
> **Notes.** The saddle has negative trace yet is unstable: a negative trace alone proves nothing. Its eigenvector for $\lambda=2$ is $(2,1)$ (from $\mathbf{A}-2\mathbf{I}=\begin{bmatrix}-1&2\\3&-6\end{bmatrix}$) — the single direction along which almost every trajectory eventually escapes. For the spiral, the sense of rotation follows from any one velocity: at $\mathbf{x}=(1,0)$, $\dot{\mathbf{x}}=(-1,1)$ has positive $\dot x_2$, so the rotation is counter-clockwise.

> [!example] Worked Example 10.37 — Linearizing the damped pendulum
> A pendulum (or, equivalently, the swing equation of a synchronous generator against an infinite bus) obeys, in normalised units, $\ddot\theta+0.5\dot\theta+\sin\theta=0$. With $\mathbf{x}=(\theta,\omega)$, $\omega=\dot\theta$:
> $$
> \mathbf{f}(\mathbf{x})=\begin{bmatrix}\omega\\-\sin\theta-0.5\,\omega\end{bmatrix},\qquad
> \mathbf{J}(\theta,\omega)=\begin{bmatrix}0&1\\-\cos\theta&-0.5\end{bmatrix}.
> $$
> **Equilibria.** $\omega=0$, $\sin\theta=0$: $\theta=0$ (hanging down) and $\theta=\pi$ (inverted).
>
> **At $(0,0)$.** $\mathbf{J}=\begin{bmatrix}0&1\\-1&-0.5\end{bmatrix}$, $\tau=-0.5$, $\Delta=1$, $\tau^2-4\Delta=-3.75<0$. Eigenvalues $-0.25\pm j\,0.968$: a **stable spiral** (decaying swings). Hyperbolic, so the nonlinear pendulum behaves the same way near the bottom.
>
> **At $(\pi,0)$.** $\mathbf{J}=\begin{bmatrix}0&1\\1&-0.5\end{bmatrix}$, $\Delta=-1<0$: a **saddle**. Eigenvalues from $\lambda^2+0.5\lambda-1=0$: $\lambda=\frac{-0.5\pm\sqrt{4.25}}{2}=0.781,\ -1.281$. Eigenvectors of this companion-type matrix are $(1,\lambda)$: the unstable direction is $(1,0.781)$, the stable direction $(1,-1.281)$. The stable eigen-line is the linear approximation of the **separatrix** — the knife-edge set of initial conditions that just reach the inverted position. In power systems the analogous saddle (the unstable equilibrium of the swing equation) defines the boundary of the region from which a generator can recover after a fault.

```desmos-graph
left=-2; right=2; bottom=-2; top=2
---
y=0.781x|red
y=-1.281x|blue
(1,0.781)|label:unstable, λ=0.781
(-1,1.281)|label:stable, λ=−1.281
```
*Eigenvector lines of the linearized pendulum at the inverted equilibrium, plotted in deviation coordinates $(\theta-\pi,\ \omega)$: trajectories approach along the blue stable line and leave along the red unstable line, the local skeleton of the saddle.*

> [!example] Worked Example 10.38 — Start-up condition of a negative-resistance oscillator
> A parallel LC tank loaded by a nonlinear negative-resistance device (tunnel diode, cross-coupled transistor pair) obeys the van der Pol equation $\ddot v-\varepsilon(1-v^2)\dot v+\omega_0^2v=0$, with $\varepsilon>0$ set by how much negative conductance exceeds the tank losses. Will the oscillation start from noise?
>
> **Linearize at $v=\dot v=0$.** With $\mathbf{x}=(v,\dot v)$, the term $\varepsilon v^2\dot v$ is cubic and drops out:
> $$
> \mathbf{J}=\begin{bmatrix}0&1\\-\omega_0^2&\varepsilon\end{bmatrix},\qquad \tau=\varepsilon>0,\ \Delta=\omega_0^2>0 .
> $$
> **Classify.** $\tau>0$, $\Delta>0$: unstable. For $\varepsilon<2\omega_0$ the eigenvalues $\frac{\varepsilon}{2}\pm j\sqrt{\omega_0^2-\varepsilon^2/4}$ are complex: an **unstable spiral**. Any tiny noise voltage grows as $e^{\varepsilon t/2}$ while oscillating at nearly $\omega_0$ — the oscillator starts.
>
> **What eigenvalues cannot tell.** The growth does not continue for ever: as $v$ grows, the nonlinear term reduces the effective negative conductance and the trajectory settles onto a **limit cycle** of amplitude about $2$. The linearization predicts *start-up* (Barkhausen's loop-gain condition is the same statement) but not the final amplitude, which is a nonlinear question (Ch. 11).

> [!example] Worked Example 10.39 — Stability of a digital resonator
> The recursion $y[n]=1.8\,y[n-1]-0.9\,y[n-2]$ is the homogeneous part of a second-order IIR filter. Is it stable, and at what frequency does it ring?
>
> **Step 1 — state form.** With $\mathbf{x}[n]=(y[n],y[n-1])$:
> $$
> \mathbf{x}[n+1]=\begin{bmatrix}1.8&-0.9\\1&0\end{bmatrix}\mathbf{x}[n],\qquad \det(\lambda\mathbf{I}-\mathbf{A})=\lambda^2-1.8\lambda+0.9 .
> $$
> **Step 2 — eigenvalues.** $\lambda=0.9\pm j\sqrt{0.9-0.81}=0.9\pm j\,0.3$, so $\lvert\lambda\rvert=\sqrt{0.9}=0.9487$ and $\angle\lambda=\arctan(0.3/0.9)=0.3218$ rad $=18.43^\circ$.
>
> **Step 3 — conclusions.** $\lvert\lambda\rvert<1$: **stable**; the ringing decays by a factor $0.9487$ per sample (to $1\%$ in $\ln0.01/\ln0.9487\approx87$ samples) and oscillates at $0.3218$ rad/sample, i.e. $0.3218/(2\pi)=0.0512$ cycles/sample — $2.46$ kHz at a $48$ kHz sample rate. Note $\lvert\lambda\rvert^2=\det\mathbf{A}=0.9$: for any second-order section $y[n]=a_1y[n-1]-a_2y[n-2]$ with complex poles, the pole radius is $\sqrt{a_2}$, so $a_2<1$ is required.
>
> **Step 4 — the boundary.** Changing $0.9$ to $1.0$ gives $\lambda=e^{\pm j\,0.451}$ on the unit circle: a marginally stable digital oscillator. In fixed-point arithmetic, rounding of the coefficient can push the poles slightly outside the unit circle — a practical reason for using structures (such as coupled-form or lattice filters) whose eigenvalues are less sensitive to coefficient quantisation.

> [!example] Worked Example 10.40 — Solving a Lyapunov equation
> For $\mathbf{A}=\begin{bmatrix}0&1\\-2&-3\end{bmatrix}$ (eigenvalues $-1,-2$) solve $\mathbf{A}^\mathsf{T}\mathbf{P}+\mathbf{P}\mathbf{A}=-\mathbf{I}$ and confirm that $\mathbf{P}\succ0$.
>
> **Step 1 — unknowns.** $\mathbf{P}=\begin{bmatrix}p_{1}&p_{2}\\p_{2}&p_{3}\end{bmatrix}$. Compute $\mathbf{A}^\mathsf{T}\mathbf{P}=\begin{bmatrix}-2p_2&-2p_3\\p_1-3p_2&p_2-3p_3\end{bmatrix}$; $\mathbf{P}\mathbf{A}$ is its transpose.
>
> **Step 2 — equate entries of the symmetric sum to $-\mathbf{I}$.**
> $$
> (1,1):\ -4p_2=-1,\qquad (2,2):\ 2p_2-6p_3=-1,\qquad (1,2):\ p_1-3p_2-2p_3=0 .
> $$
> So $p_2=0.25$, $p_3=(1+0.5)/6=0.25$, $p_1=3(0.25)+2(0.25)=1.25$:
> $$
> \mathbf{P}=\begin{bmatrix}1.25&0.25\\0.25&0.25\end{bmatrix}.
> $$
> **Step 3 — definiteness.** $p_1=1.25>0$, $\det\mathbf{P}=0.3125-0.0625=0.25>0$: $\mathbf{P}\succ0$, confirming that $\mathbf{A}$ is Hurwitz without computing an eigenvalue. (A numerical `solve_continuous_lyapunov` returns the same $\mathbf{P}$.)
>
> **Step 4 — interpretation.** $V(\mathbf{x})=1.25x_1^2+0.5x_1x_2+0.25x_2^2$ is an energy-like function with $\dot V=-(x_1^2+x_2^2)$ along every trajectory. Its level sets are ellipses that every trajectory crosses inward.

> [!warning] Common Pitfall — stable eigenvalues, large transients
> Eigenvalues describe the behaviour as $t\to\infty$, not the size of the response on the way. For the stable but strongly **non-normal** matrix $\mathbf{A}=\begin{bmatrix}-1&50\\0&-2\end{bmatrix}$ (eigenvalues $-1,-2$),
> $$
> e^{\mathbf{A}t}=\begin{bmatrix}e^{-t}&50(e^{-t}-e^{-2t})\\0&e^{-2t}\end{bmatrix},
> $$
> and the $(1,2)$ entry reaches $50(0.5-0.25)=12.5$ at $t=\ln2$: a unit initial condition on $x_2$ produces a twelve-fold excursion in $x_1$ before decaying. If the model came from linearization, such an excursion can drive the real system out of the region where the linear model is valid — or saturate an amplifier. The eigenvector matrix of such an $\mathbf{A}$ is ill-conditioned (Section 10.3); singular values of $e^{\mathbf{A}t}$, or a Lyapunov function, quantify the transient.

> [!tip] Engineering Note — stability margins and eigenvalue location
> Asymptotic stability is a yes/no property, but designs need margins. Requiring all eigenvalues to satisfy $\operatorname{Re}\lambda\le-\alpha$ guarantees decay at least as fast as $e^{-\alpha t}$ (settling time about $4/\alpha$); requiring them inside a sector $\lvert\angle(-\lambda)\rvert\le\arccos\zeta_{\min}$ guarantees a minimum damping ratio and hence bounded overshoot. In power electronics and power systems, the eigenvalue plot of the linearized model versus a parameter (load, gain, line impedance) — a **root locus** — shows how close each mode is to crossing the imaginary axis; the crossing point, where a complex pair reaches $\operatorname{Re}\lambda=0$, is a **Hopf bifurcation**, the birth of a sustained oscillation as in Example 10.38.

## Summary

| Concept | Key formula or statement |
|---|---|
| Eigenpair | $\mathbf{A}\mathbf{v}=\lambda\mathbf{v}$, $\mathbf{v}\ne\mathbf{0}$; eigenspace $E_\lambda=\operatorname{null}(\mathbf{A}-\lambda\mathbf{I})$ |
| Characteristic equation | $\det(\mathbf{A}-\lambda\mathbf{I})=0$, degree $n$, $n$ roots in $\mathbb{C}$ |
| $2\times2$ shortcut | $\lambda^2-(\operatorname{tr}\mathbf{A})\lambda+\det\mathbf{A}=0$ |
| $3\times3$ shortcut | $\lambda^3-(\operatorname{tr}\mathbf{A})\lambda^2+E_2\lambda-\det\mathbf{A}=0$, $E_2$ = sum of $2\times2$ principal minors |
| Trace and determinant | $\sum\lambda_i=\operatorname{tr}\mathbf{A}$, $\prod\lambda_i=\det\mathbf{A}$ |
| Cayley–Hamilton | $p(\mathbf{A})=\mathbf{0}$; e.g. $\mathbf{A}^2-(\operatorname{tr}\mathbf{A})\mathbf{A}+(\det\mathbf{A})\mathbf{I}=\mathbf{0}$ |
| Multiplicities | $1\le g(\lambda)\le m(\lambda)$; defective if $g<m$ |
| Functions of $\mathbf{A}$ | eigenvalues of $q(\mathbf{A})$, $\mathbf{A}^{-1}$, $e^{\mathbf{A}t}$ are $q(\lambda)$, $1/\lambda$, $e^{\lambda t}$ (same eigenvectors) |
| Diagonalization | $\mathbf{A}=\mathbf{P}\mathbf{D}\mathbf{P}^{-1}$ iff $n$ independent eigenvectors; distinct eigenvalues suffice |
| Real block form | $\lambda=\sigma\pm j\omega$, $\mathbf{v}=\mathbf{a}+j\mathbf{b}$: $[\mathbf{a}\ \mathbf{b}]^{-1}\mathbf{A}[\mathbf{a}\ \mathbf{b}]=\begin{bmatrix}\sigma&\omega\\-\omega&\sigma\end{bmatrix}$ |
| Jordan block powers | $\mathbf{J}_2(\lambda)^k=\begin{bmatrix}\lambda^k&k\lambda^{k-1}\\0&\lambda^k\end{bmatrix}$, $e^{\mathbf{J}_2t}=e^{\lambda t}\begin{bmatrix}1&t\\0&1\end{bmatrix}$ |
| Matrix powers | $\mathbf{A}^k=\mathbf{P}\mathbf{D}^k\mathbf{P}^{-1}$; $\mathbf{x}[k]=\sum c_i\lambda_i^k\mathbf{v}_i$; $\mathbf{A}^k\to\mathbf{0}\iff\rho(\mathbf{A})<1$ |
| Matrix exponential | $e^{\mathbf{A}t}=\sum\mathbf{A}^kt^k/k!=\mathbf{P}e^{\mathbf{D}t}\mathbf{P}^{-1}$; $\mathbf{x}(t)=e^{\mathbf{A}t}\mathbf{x}(0)=\sum c_ie^{\lambda_it}\mathbf{v}_i$ |
| Sylvester interpolation | $f(\mathbf{A})=\sum_{k<n}\alpha_k\mathbf{A}^k$ with $f(\lambda_i)=\sum\alpha_k\lambda_i^k$ (plus derivatives for repeated $\lambda_i$) |
| Spectral theorem | $\mathbf{A}=\mathbf{A}^\mathsf{T}\iff\mathbf{A}=\mathbf{Q}\boldsymbol{\Lambda}\mathbf{Q}^\mathsf{T}=\sum\lambda_i\mathbf{q}_i\mathbf{q}_i^\mathsf{T}$, $\lambda_i$ real, $\mathbf{Q}^\mathsf{T}\mathbf{Q}=\mathbf{I}$ |
| Rayleigh quotient | $\lambda_{\min}\le\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}/\mathbf{x}^\mathsf{T}\mathbf{x}\le\lambda_{\max}$ |
| Principal axes | $\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}=\sum\lambda_iy_i^2$ with $\mathbf{y}=\mathbf{Q}^\mathsf{T}\mathbf{x}$; ellipse semi-axes $1/\sqrt{\lambda_i}$ for $\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}=1$ |
| Positive definite | $\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}>0$ $\iff$ all $\lambda_i>0$ $\iff$ all leading minors $>0$ $\iff$ all pivots $>0$ $\iff$ $\mathbf{A}=\mathbf{L}\mathbf{L}^\mathsf{T}$ |
| $2\times2$ PD test | $a_{11}>0$ and $a_{11}a_{22}-a_{12}^2>0$ (coupled inductors: $M^2<L_1L_2$) |
| SVD | $\mathbf{A}=\mathbf{U}\boldsymbol{\Sigma}\mathbf{V}^\mathsf{T}=\sum\sigma_i\mathbf{u}_i\mathbf{v}_i^\mathsf{T}$, $\sigma_i^2=\lambda_i(\mathbf{A}^\mathsf{T}\mathbf{A})$, $\mathbf{A}\mathbf{v}_i=\sigma_i\mathbf{u}_i$ |
| Norms and conditioning | $\lVert\mathbf{A}\rVert_2=\sigma_1$, $\lVert\mathbf{A}\rVert_F^2=\sum\sigma_i^2$, $\kappa_2=\sigma_{\max}/\sigma_{\min}$ |
| Eckart–Young | $\min_{\operatorname{rank}\mathbf{B}\le k}\lVert\mathbf{A}-\mathbf{B}\rVert_2=\sigma_{k+1}$, attained by truncated SVD |
| Pseudoinverse | $\mathbf{A}^+=\sum_{i\le r}\sigma_i^{-1}\mathbf{v}_i\mathbf{u}_i^\mathsf{T}$; $\mathbf{A}^+\mathbf{b}$ = minimum-norm least-squares solution |
| Normal modes | $\mathbf{K}\boldsymbol{\phi}=\omega^2\mathbf{M}\boldsymbol{\phi}$; $\boldsymbol{\Phi}^\mathsf{T}\mathbf{M}\boldsymbol{\Phi}=\mathbf{I}$, $\boldsymbol{\Phi}^\mathsf{T}\mathbf{K}\boldsymbol{\Phi}=\boldsymbol{\Omega}^2$; $\ddot q_i+\omega_i^2q_i=\boldsymbol{\phi}_i^\mathsf{T}\mathbf{f}$ |
| Modal response (1st order) | $\mathbf{x}(t)=\sum\mathbf{v}_i(\mathbf{w}_i^\mathsf{T}\mathbf{x}(0))e^{\lambda_it}$, $\mathbf{w}_i^\mathsf{T}$ = rows of $\mathbf{P}^{-1}$ |
| Continuous-time stability | asymptotically stable $\iff\operatorname{Re}\lambda_i<0$ for all $i$ |
| Discrete-time stability | asymptotically stable $\iff\lvert\lambda_i\rvert<1$ for all $i$ |
| $2\times2$ classification | stable $\iff\tau<0,\ \Delta>0$; saddle if $\Delta<0$; spiral if $\tau^2<4\Delta$; node if $\tau^2>4\Delta$ |
| Third-order Routh–Hurwitz | $s^3+a_2s^2+a_1s+a_0$ stable $\iff a_i>0$ and $a_2a_1>a_0$ |
| Linearization | $\dot{\boldsymbol{\xi}}=\mathbf{J}\boldsymbol{\xi}$, $\mathbf{J}=\partial\mathbf{f}/\partial\mathbf{x}$ at equilibrium; conclusive if no $\operatorname{Re}\lambda=0$ |
| Lyapunov equation | $\mathbf{A}$ Hurwitz $\iff\mathbf{A}^\mathsf{T}\mathbf{P}+\mathbf{P}\mathbf{A}=-\mathbf{Q}$ has $\mathbf{P}\succ0$ for $\mathbf{Q}\succ0$ |

## Further reading

- G. Strang, *Introduction to Linear Algebra*, 6th ed., Wellesley–Cambridge Press, 2023 — chapters on eigenvalues, symmetric and positive-definite matrices and the SVD, with exceptional geometric intuition.
- E. Kreyszig, *Advanced Engineering Mathematics*, 10th ed., Wiley, 2011 — matrix eigenvalue problems, quadratic forms and systems of ODEs from the engineering side.
- C. D. Meyer, *Matrix Analysis and Applied Linear Algebra*, 2nd ed., SIAM, 2023 — rigorous treatment of Jordan forms, matrix functions, the SVD and Perron–Frobenius theory for Markov chains.
- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — the reference for how eigenvalues, the SVD and matrix exponentials are actually computed.
- C.-T. Chen, *Linear System Theory and Design*, 4th ed., Oxford University Press, 2013 — state-space modal decomposition, Lyapunov equations and stability of linear systems for electrical engineers.

---

[← Previous](09-vector-spaces-and-linear-transformations.md) · [Table of Contents](../00-front/table-of-contents.md) · [Next →](../part-4-differential-equations/11-ordinary-differential-equations.md)
