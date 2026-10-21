# Chapter 8 — Matrices and Linear Systems

[← Previous](../part-2-calculus/07-vector-calculus.md) · [Table of Contents](../00-front/table-of-contents.md) · [Next →](../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md)

> **Prerequisites:** [Ch. 1](../part-1-foundations/01-mathematical-language-and-algebra.md), [Ch. 2](../part-1-foundations/02-complex-numbers.md); the Jacobian matrix of [Ch. 6](../part-2-calculus/06-multivariable-calculus.md) is helpful but not essential.
> **Learning outcomes:** After this chapter you can
> - manipulate matrices fluently (sums, products, transposes, conjugate transposes, block forms) and read circuit equations as matrix equations;
> - evaluate determinants by cofactor expansion and by row reduction, and use them to test invertibility and to apply Cramer's rule;
> - compute inverses, determine rank, and classify a linear system as having a unique solution, no solution or infinitely many solutions;
> - solve systems by Gaussian elimination with partial pivoting, by Gauss–Jordan elimination, and by $\mathbf{LU}$ and Cholesky factorisation, and count the work each needs;
> - measure the sensitivity of a solution with norms and the condition number, and recognise ill-conditioned problems;
> - exploit sparsity (banded and tridiagonal solvers, compressed storage, iterative methods) in the large systems produced by circuit simulators and field solvers.

A linear circuit with $n$ unknown node voltages produces $n$ linear equations in $n$ unknowns. A finite-difference model of the potential in a capacitor produces ten thousand. A power-flow study of a national grid produces a hundred thousand, solved again at every iteration of Newton's method. Whatever the size, the mathematical object is the same: a **linear system** $\mathbf{A}\mathbf{x}=\mathbf{b}$. This chapter develops the algebra of matrices needed to write such systems compactly, the theory that tells us when they can be solved, and the algorithms — elimination and factorisation — that every circuit simulator, field solver and estimation routine actually runs. Chapter 9 then looks at the same objects geometrically, as vector spaces and linear maps, and Chapter 10 studies their eigenvalues.

Throughout the chapter one small circuit serves as a running example. It has three non-reference nodes and six resistive branches:

| Branch | Connects | Conductance |
|---|---|---|
| $a$ | node 1 – ground | $1\ \mathrm{mS}$ |
| $b$ | node 1 – node 2 | $2\ \mathrm{mS}$ |
| $c$ | node 2 – node 3 | $1\ \mathrm{mS}$ |
| $d$ | node 2 – ground | $1\ \mathrm{mS}$ |
| $e$ | node 3 – ground | $2\ \mathrm{mS}$ |
| $f$ | node 1 – node 3 | $1\ \mathrm{mS}$ |

Independent current sources inject $8\ \mathrm{mA}$, $2\ \mathrm{mA}$ and $1\ \mathrm{mA}$ into nodes 1, 2 and 3. Kirchhoff's current law (KCL) at each node, with currents leaving the node counted positive, gives

$$
\begin{aligned}
(1+2+1)v_1 - 2v_2 - 1v_3 &= 8,\\
-2v_1 + (2+1+1)v_2 - 1v_3 &= 2,\\
-1v_1 - 1v_2 + (1+2+1)v_3 &= 1,
\end{aligned}
$$

with conductances in mS, currents in mA and therefore voltages in V (since $\mathrm{mA}/\mathrm{mS}=\mathrm{V}$). In matrix form,

$$
\underbrace{\begin{bmatrix} 4 & -2 & -1\\ -2 & 4 & -1\\ -1 & -1 & 4\end{bmatrix}}_{\mathbf{G}}
\underbrace{\begin{bmatrix} v_1\\ v_2\\ v_3\end{bmatrix}}_{\mathbf{v}}
=
\underbrace{\begin{bmatrix} 8\\ 2\\ 1\end{bmatrix}}_{\mathbf{i}_s}.
$$

We will find that $\mathbf{v}=(4,3,2)^\mathsf{T}\ \mathrm{V}$ in at least six different ways before the chapter is over, and each way teaches something different about $\mathbf{G}$.

## 8.1 Matrix operations

### Why matrices

Writing the nodal equations above as three separate lines is fine for three nodes. For three thousand it is hopeless — and, more importantly, it hides structure. The matrix $\mathbf{G}$ is symmetric; its diagonal entries dominate its rows; its off-diagonal entries are non-positive; most of its entries in a large circuit are zero. Each of these properties has a name, a theorem attached to it, and an algorithmic payoff. Matrix notation is what lets us see them.

### Definitions

> [!definition] Matrix
> An $m\times n$ **matrix** $\mathbf{A}$ is a rectangular array of $mn$ numbers (real or complex) arranged in $m$ rows and $n$ columns. The entry in row $i$, column $j$ is written $a_{ij}$ or $[\mathbf{A}]_{ij}$:
> $$
> \mathbf{A}=\begin{bmatrix} a_{11} & a_{12} & \cdots & a_{1n}\\ a_{21} & a_{22} & \cdots & a_{2n}\\ \vdots & \vdots & \ddots & \vdots\\ a_{m1} & a_{m2} & \cdots & a_{mn}\end{bmatrix}=[a_{ij}].
> $$
> The pair $(m,n)$ is the **size** (or shape). A matrix with $m=n$ is **square** of **order** $n$. An $n\times 1$ matrix is a **column vector**, a $1\times n$ matrix a **row vector**. The set of real $m\times n$ matrices is written $\mathbb{R}^{m\times n}$, the complex ones $\mathbb{C}^{m\times n}$.

Vectors in this book are columns unless stated otherwise: $\mathbf{x}=(x_1,\dots,x_n)^\mathsf{T}$ is a convenient way to print a column on one line.

Some square matrices occur so often they have names.

| Name | Definition | Example ($3\times3$ pattern) |
|---|---|---|
| Diagonal | $a_{ij}=0$ for $i\ne j$ | $\operatorname{diag}(d_1,d_2,d_3)$ |
| Identity $\mathbf{I}$ (or $\mathbf{I}_n$) | diagonal with all $d_i=1$ | $\delta_{ij}$ (Kronecker delta) |
| Zero $\mathbf{0}$ | all entries zero | |
| Upper triangular | $a_{ij}=0$ for $i>j$ | nonzeros on and above diagonal |
| Lower triangular | $a_{ij}=0$ for $i<j$ | nonzeros on and below diagonal |
| Tridiagonal | $a_{ij}=0$ for $\lvert i-j\rvert>1$ | ladder networks, 1-D finite differences |
| Symmetric | $a_{ij}=a_{ji}$ | nodal matrix of a reciprocal circuit |
| Skew-symmetric | $a_{ij}=-a_{ji}$ (so $a_{ii}=0$) | gyrator, cross-product matrix |
| Hermitian | $a_{ij}=a_{ji}^*$ | power/covariance matrices with complex entries |

The **main diagonal** of a square matrix is $a_{11},a_{22},\dots,a_{nn}$. The **trace** is the sum of the diagonal entries, $\operatorname{tr}\mathbf{A}=\sum_i a_{ii}$; for $\mathbf{G}$ above, $\operatorname{tr}\mathbf{G}=12\ \mathrm{mS}$, which is twice the total conductance of branches $b,c,f$ plus once that of $a,d,e$ (each inter-node branch is counted at both its ends).

### Equality, addition and scalar multiplication

> [!definition] Linear operations on matrices
> Two matrices are **equal** if they have the same size and all corresponding entries are equal. For $\mathbf{A},\mathbf{B}\in\mathbb{C}^{m\times n}$ and a scalar $c$:
> $$
> [\mathbf{A}+\mathbf{B}]_{ij}=a_{ij}+b_{ij},\qquad [c\mathbf{A}]_{ij}=c\,a_{ij}.
> $$
> Addition of matrices of different sizes is undefined.

These operations inherit every rule of ordinary arithmetic, entry by entry:

$$
\begin{aligned}
\mathbf{A}+\mathbf{B}&=\mathbf{B}+\mathbf{A}, & (\mathbf{A}+\mathbf{B})+\mathbf{C}&=\mathbf{A}+(\mathbf{B}+\mathbf{C}),\\
\mathbf{A}+\mathbf{0}&=\mathbf{A}, & \mathbf{A}+(-\mathbf{A})&=\mathbf{0},\\
c(\mathbf{A}+\mathbf{B})&=c\mathbf{A}+c\mathbf{B}, & (c+d)\mathbf{A}&=c\mathbf{A}+d\mathbf{A},\quad (cd)\mathbf{A}=c(d\mathbf{A}).
\end{aligned}
$$

In the language of Chapter 9, $\mathbb{R}^{m\times n}$ is a vector space of dimension $mn$.

Addition has a direct circuit meaning. Each branch of the running example contributes a small "stamp" to $\mathbf{G}$, and $\mathbf{G}$ is the sum of the stamps. Branch $b$ (2 mS between nodes 1 and 2) contributes

$$
\mathbf{G}_b = 2\begin{bmatrix} 1 & -1 & 0\\ -1 & 1 & 0\\ 0&0&0\end{bmatrix},
$$

branch $a$ (1 mS from node 1 to ground) contributes $1$ at position $(1,1)$ only, and so on. Summing the six stamps reproduces $\mathbf{G}$. This **stamping** procedure is exactly how SPICE-type simulators assemble their matrices: loop over the elements, add each element's stamp into the global matrix.

### Transpose and conjugate transpose

> [!definition] Transpose, conjugate transpose
> The **transpose** of $\mathbf{A}\in\mathbb{C}^{m\times n}$ is the $n\times m$ matrix $\mathbf{A}^\mathsf{T}$ with $[\mathbf{A}^\mathsf{T}]_{ij}=a_{ji}$ (rows become columns). The **conjugate transpose** (Hermitian transpose) is $\mathbf{A}^\mathsf{H}=(\mathbf{A}^*)^\mathsf{T}$, i.e. $[\mathbf{A}^\mathsf{H}]_{ij}=a_{ji}^*$. For real matrices $\mathbf{A}^\mathsf{H}=\mathbf{A}^\mathsf{T}$.

Properties, each checked entry by entry:

$$
(\mathbf{A}^\mathsf{T})^\mathsf{T}=\mathbf{A},\qquad (\mathbf{A}+\mathbf{B})^\mathsf{T}=\mathbf{A}^\mathsf{T}+\mathbf{B}^\mathsf{T},\qquad (c\mathbf{A})^\mathsf{T}=c\mathbf{A}^\mathsf{T},\qquad (c\mathbf{A})^\mathsf{H}=c^*\mathbf{A}^\mathsf{H}.
$$

The rule for products, $(\mathbf{A}\mathbf{B})^\mathsf{T}=\mathbf{B}^\mathsf{T}\mathbf{A}^\mathsf{T}$, is proved in §8.2.

A square matrix is symmetric exactly when $\mathbf{A}^\mathsf{T}=\mathbf{A}$, skew-symmetric when $\mathbf{A}^\mathsf{T}=-\mathbf{A}$, and Hermitian when $\mathbf{A}^\mathsf{H}=\mathbf{A}$. A Hermitian matrix has a real diagonal, since $a_{ii}=a_{ii}^*$.

> [!theorem] Symmetric/skew decomposition
> Every square matrix splits uniquely into a symmetric and a skew-symmetric part:
> $$
> \mathbf{A}=\underbrace{\tfrac12(\mathbf{A}+\mathbf{A}^\mathsf{T})}_{\text{symmetric}}+\underbrace{\tfrac12(\mathbf{A}-\mathbf{A}^\mathsf{T})}_{\text{skew-symmetric}}.
> $$
> Likewise every complex square matrix is a Hermitian part plus $j$ times another Hermitian part: $\mathbf{A}=\mathbf{H}_1+j\mathbf{H}_2$ with $\mathbf{H}_1=\tfrac12(\mathbf{A}+\mathbf{A}^\mathsf{H})$, $\mathbf{H}_2=\tfrac1{2j}(\mathbf{A}-\mathbf{A}^\mathsf{H})$.

*Proof sketch.* Transposing the first bracket returns it unchanged; transposing the second negates it; their sum is $\mathbf{A}$. For uniqueness, if $\mathbf{A}=\mathbf{S}+\mathbf{K}$ with $\mathbf{S}$ symmetric and $\mathbf{K}$ skew, then $\mathbf{A}^\mathsf{T}=\mathbf{S}-\mathbf{K}$; adding and subtracting recovers the two brackets. $\square$

> [!example] Worked Example 8.1 — Basic operations
> Let
> $$
> \mathbf{A}=\begin{bmatrix} 1 & 2 & 0\\ 3 & -1 & 4\end{bmatrix},\qquad \mathbf{B}=\begin{bmatrix} 0 & 1 & -2\\ 5 & 2 & 1\end{bmatrix}.
> $$
> Compute $2\mathbf{A}-3\mathbf{B}$ and $(\mathbf{A}+\mathbf{B})^\mathsf{T}$.
>
> **Step 1.** Both are $2\times3$, so linear combinations are defined. Entry by entry,
> $$
> 2\mathbf{A}-3\mathbf{B}=\begin{bmatrix} 2-0 & 4-3 & 0+6\\ 6-15 & -2-6 & 8-3\end{bmatrix}=\begin{bmatrix} 2 & 1 & 6\\ -9 & -8 & 5\end{bmatrix}.
> $$
> **Step 2.** $\mathbf{A}+\mathbf{B}=\begin{bmatrix} 1 & 3 & -2\\ 8 & 1 & 5\end{bmatrix}$, so
> $$
> (\mathbf{A}+\mathbf{B})^\mathsf{T}=\begin{bmatrix} 1 & 8\\ 3 & 1\\ -2 & 5\end{bmatrix},
> $$
> a $3\times 2$ matrix. Check: $\mathbf{A}^\mathsf{T}+\mathbf{B}^\mathsf{T}$ gives the same thing.

> [!example] Worked Example 8.2 — A nonreciprocal circuit and its symmetric part
> A transistor small-signal model adds a voltage-controlled current source $g_m v_1$ flowing *out* of node 2 to the running example with $g_m=3\ \mathrm{mS}$ (the source is controlled by $v_1$). KCL at node 2 gains a term $+3v_1$, so the nodal matrix becomes
> $$
> \mathbf{Y}=\begin{bmatrix} 4 & -2 & -1\\ 1 & 4 & -1\\ -1 & -1 & 4\end{bmatrix}.
> $$
> Split it into symmetric and skew parts.
>
> **Step 1.** $\mathbf{Y}^\mathsf{T}=\begin{bmatrix} 4 & 1 & -1\\ -2 & 4 & -1\\ -1 & -1 & 4\end{bmatrix}$.
>
> **Step 2.**
> $$
> \mathbf{S}=\tfrac12(\mathbf{Y}+\mathbf{Y}^\mathsf{T})=\begin{bmatrix} 4 & -\tfrac12 & -1\\ -\tfrac12 & 4 & -1\\ -1 & -1 & 4\end{bmatrix},\qquad
> \mathbf{K}=\tfrac12(\mathbf{Y}-\mathbf{Y}^\mathsf{T})=\begin{bmatrix} 0 & -\tfrac32 & 0\\ \tfrac32 & 0 & 0\\ 0&0&0\end{bmatrix}.
> $$
> **Step 3 (check).** $\mathbf{S}+\mathbf{K}=\mathbf{Y}$, $\mathbf{S}^\mathsf{T}=\mathbf{S}$, $\mathbf{K}^\mathsf{T}=-\mathbf{K}$.
>
> **Reading.** The skew part is entirely due to the controlled source: a passive reciprocal network has a symmetric nodal matrix, and the skew part measures how far a network is from reciprocity. In §8.2 we will see that the power absorbed by the network, $\mathbf{v}^\mathsf{T}\mathbf{Y}\mathbf{v}$, depends only on $\mathbf{S}$.

> [!example] Worked Example 8.3 — Hermitian part of an AC admittance matrix
> At angular frequency $\omega$ a two-node network has phasor nodal equations $\mathbf{Y}\mathbf{V}=\mathbf{I}$ with
> $$
> \mathbf{Y}=\begin{bmatrix} 2+j3 & -j3\\ -j3 & 1+j4\end{bmatrix}\ \mathrm{mS}.
> $$
> Find $\mathbf{Y}^\mathsf{H}$ and the Hermitian part $\mathbf{H}_1=\tfrac12(\mathbf{Y}+\mathbf{Y}^\mathsf{H})$.
>
> **Step 1.** Conjugate every entry, then transpose: $\mathbf{Y}^*=\begin{bmatrix} 2-j3 & j3\\ j3 & 1-j4\end{bmatrix}$, which is already symmetric, so $\mathbf{Y}^\mathsf{H}=\mathbf{Y}^*$.
>
> **Step 2.** $\mathbf{H}_1=\tfrac12(\mathbf{Y}+\mathbf{Y}^*)=\operatorname{Re}\{\mathbf{Y}\}=\begin{bmatrix} 2&0\\0&1\end{bmatrix}\ \mathrm{mS}$.
>
> **Reading.** For a reciprocal network $\mathbf{Y}$ is complex *symmetric* ($\mathbf{Y}^\mathsf{T}=\mathbf{Y}$) but not Hermitian. Its Hermitian part is the conductance matrix: average power is $P=\tfrac12\operatorname{Re}\{\mathbf{V}^\mathsf{H}\mathbf{I}\}=\tfrac12\mathbf{V}^\mathsf{H}\mathbf{H}_1\mathbf{V}$, which only involves the resistive branches. The $j3$ susceptance (a capacitor between the nodes) and the extra $j1$ at node 2 (a capacitor to ground) store energy but dissipate none.

> [!warning] Common Pitfall
> For complex matrices, $\mathbf{A}^\mathsf{T}$ and $\mathbf{A}^\mathsf{H}$ are different, and MATLAB's `'` operator is the *conjugate* transpose; `.'` is the plain transpose. Using `'` on a complex symmetric admittance matrix silently conjugates it — a classic source of sign errors in reactive power. NumPy's `.T` is the plain transpose; use `.conj().T` for $\mathbf{A}^\mathsf{H}$.

> [!tip] Engineering Note
> Matrix symmetry is not a mathematical nicety; it is reciprocity. Tellegen's theorem implies that any network made of resistors, capacitors, inductors and coupled inductors has a symmetric nodal admittance matrix. Simulators check symmetry to choose faster solvers (Cholesky rather than general $\mathbf{LU}$, §8.9), and a non-symmetric matrix in a model that should be passive is a quick bug detector.

## 8.2 Matrix multiplication

### Motivation: composition

Multiplication of matrices is defined the way it is because it represents *composition* of linear maps. If $\mathbf{y}=\mathbf{B}\mathbf{x}$ and $\mathbf{z}=\mathbf{A}\mathbf{y}$, we want a single matrix $\mathbf{C}$ with $\mathbf{z}=\mathbf{C}\mathbf{x}$. Writing out the sums,

$$
z_i=\sum_k a_{ik}y_k=\sum_k a_{ik}\sum_j b_{kj}x_j=\sum_j\Big(\sum_k a_{ik}b_{kj}\Big)x_j,
$$

so the right definition is $c_{ij}=\sum_k a_{ik}b_{kj}$. Two-port networks in cascade, coordinate rotations followed by scalings, and a filter followed by a sampler are all compositions, and all become products.

> [!definition] Matrix product
> If $\mathbf{A}$ is $m\times p$ and $\mathbf{B}$ is $p\times n$, the product $\mathbf{C}=\mathbf{A}\mathbf{B}$ is the $m\times n$ matrix
> $$
> c_{ij}=\sum_{k=1}^{p}a_{ik}b_{kj},\qquad i=1,\dots,m,\ j=1,\dots,n.
> $$
> The product is defined only when the inner dimensions agree: $(m\times p)(p\times n)\to(m\times n)$.

The entry $c_{ij}$ is the **dot product of row $i$ of $\mathbf{A}$ with column $j$ of $\mathbf{B}$**. Three other readings are equally important:

1. **Column view.** Column $j$ of $\mathbf{A}\mathbf{B}$ is $\mathbf{A}$ times column $j$ of $\mathbf{B}$. In particular $\mathbf{A}\mathbf{x}=x_1\mathbf{a}_1+x_2\mathbf{a}_2+\dots+x_n\mathbf{a}_n$ is a **linear combination of the columns** $\mathbf{a}_k$ of $\mathbf{A}$. Solving $\mathbf{A}\mathbf{x}=\mathbf{b}$ means finding the weights that combine the columns into $\mathbf{b}$.
2. **Row view.** Row $i$ of $\mathbf{A}\mathbf{B}$ is row $i$ of $\mathbf{A}$ times $\mathbf{B}$: a linear combination of the rows of $\mathbf{B}$. Row operations in elimination (§8.7) are left-multiplications.
3. **Outer-product view.** $\mathbf{A}\mathbf{B}=\sum_k \mathbf{a}_k\mathbf{b}_k^\mathsf{T}$, a sum of rank-one matrices, column $k$ of $\mathbf{A}$ times row $k$ of $\mathbf{B}$. This is the view behind $\mathbf{LU}$ factorisation and the Sherman–Morrison formula.

For the running example the column view reads

$$
\mathbf{G}\mathbf{v}=4\begin{bmatrix}4\\-2\\-1\end{bmatrix}+3\begin{bmatrix}-2\\4\\-1\end{bmatrix}+2\begin{bmatrix}-1\\-1\\4\end{bmatrix}=\begin{bmatrix}16-6-2\\-8+12-2\\-4-3+8\end{bmatrix}=\begin{bmatrix}8\\2\\1\end{bmatrix},
$$

which confirms the claimed solution $\mathbf{v}=(4,3,2)^\mathsf{T}$. Column $k$ of $\mathbf{G}$ is the vector of currents that must be injected to raise node $k$ to 1 V with every other node held at 0 V — the circuit meaning of the column view.

### Algebraic rules

> [!theorem] Properties of the matrix product
> Whenever the sizes allow the products:
> 1. Associative: $(\mathbf{A}\mathbf{B})\mathbf{C}=\mathbf{A}(\mathbf{B}\mathbf{C})$.
> 2. Distributive: $\mathbf{A}(\mathbf{B}+\mathbf{C})=\mathbf{A}\mathbf{B}+\mathbf{A}\mathbf{C}$, $(\mathbf{A}+\mathbf{B})\mathbf{C}=\mathbf{A}\mathbf{C}+\mathbf{B}\mathbf{C}$.
> 3. Scalars commute: $c(\mathbf{A}\mathbf{B})=(c\mathbf{A})\mathbf{B}=\mathbf{A}(c\mathbf{B})$.
> 4. Identity: $\mathbf{I}_m\mathbf{A}=\mathbf{A}=\mathbf{A}\mathbf{I}_n$ for $\mathbf{A}\in\mathbb{C}^{m\times n}$.
> 5. Transpose reverses order: $(\mathbf{A}\mathbf{B})^\mathsf{T}=\mathbf{B}^\mathsf{T}\mathbf{A}^\mathsf{T}$, $(\mathbf{A}\mathbf{B})^\mathsf{H}=\mathbf{B}^\mathsf{H}\mathbf{A}^\mathsf{H}$.
> 6. Trace is cyclic: $\operatorname{tr}(\mathbf{A}\mathbf{B})=\operatorname{tr}(\mathbf{B}\mathbf{A})$ (when both products are square).
>
> But in general $\mathbf{A}\mathbf{B}\ne\mathbf{B}\mathbf{A}$ — matrix multiplication is **not commutative**.

*Proofs.* (1) The $(i,j)$ entry of both sides is $\sum_k\sum_l a_{ik}b_{kl}c_{lj}$; only the order of summation differs. (5) $[(\mathbf{A}\mathbf{B})^\mathsf{T}]_{ij}=[\mathbf{A}\mathbf{B}]_{ji}=\sum_k a_{jk}b_{ki}=\sum_k [\mathbf{B}^\mathsf{T}]_{ik}[\mathbf{A}^\mathsf{T}]_{kj}=[\mathbf{B}^\mathsf{T}\mathbf{A}^\mathsf{T}]_{ij}$. (6) $\operatorname{tr}(\mathbf{A}\mathbf{B})=\sum_i\sum_k a_{ik}b_{ki}=\sum_k\sum_i b_{ki}a_{ik}=\operatorname{tr}(\mathbf{B}\mathbf{A})$. $\square$

Other familiar rules of arithmetic also fail:

- $\mathbf{A}\mathbf{B}=\mathbf{0}$ does **not** imply $\mathbf{A}=\mathbf{0}$ or $\mathbf{B}=\mathbf{0}$. Example: $\begin{bmatrix}1&1\\1&1\end{bmatrix}\begin{bmatrix}1&-1\\-1&1\end{bmatrix}=\mathbf{0}$.
- $\mathbf{A}\mathbf{B}=\mathbf{A}\mathbf{C}$ does **not** imply $\mathbf{B}=\mathbf{C}$ unless $\mathbf{A}$ is invertible (§8.4).
- $(\mathbf{A}+\mathbf{B})^2=\mathbf{A}^2+\mathbf{A}\mathbf{B}+\mathbf{B}\mathbf{A}+\mathbf{B}^2$, which equals $\mathbf{A}^2+2\mathbf{A}\mathbf{B}+\mathbf{B}^2$ only if $\mathbf{A}$ and $\mathbf{B}$ commute.

### Quadratic forms and power

For a square $\mathbf{A}$ and vector $\mathbf{x}$ the scalar $\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}=\sum_i\sum_j a_{ij}x_ix_j$ is a **quadratic form**. Because $\mathbf{x}^\mathsf{T}\mathbf{K}\mathbf{x}$ is a scalar equal to its own transpose $\mathbf{x}^\mathsf{T}\mathbf{K}^\mathsf{T}\mathbf{x}=-\mathbf{x}^\mathsf{T}\mathbf{K}\mathbf{x}$ for skew $\mathbf{K}$, it vanishes, and only the symmetric part of $\mathbf{A}$ contributes to a quadratic form. For a resistive network, the total power absorbed is

$$
P=\mathbf{v}^\mathsf{T}\mathbf{i}_s=\mathbf{v}^\mathsf{T}\mathbf{G}\mathbf{v}.
$$

For the running example $P=(4,3,2)\cdot(8,2,1)=32+6+2=40\ \mathrm{mW}$. Since every resistor absorbs non-negative power, $\mathbf{v}^\mathsf{T}\mathbf{G}\mathbf{v}>0$ for every nonzero $\mathbf{v}$ when the circuit is connected to ground: $\mathbf{G}$ is **positive definite**. That property is studied in Chapter 10 and exploited by the Cholesky factorisation in §8.9.

### Incidence matrices: building $\mathbf{G}$ by multiplication

The cleanest derivation of nodal analysis is a matrix product. Number the nodes and branches and give each branch a reference direction (from its first to its second node in the table above).

> [!definition] Reduced incidence matrix
> For a network with $n$ non-reference nodes and $b$ branches, the **reduced incidence matrix** $\mathbf{A}$ is $n\times b$ with
> $$
> a_{ik}=\begin{cases}+1 & \text{branch } k \text{ leaves node } i,\\ -1 & \text{branch } k \text{ enters node } i,\\ 0 & \text{otherwise.}\end{cases}
> $$
> The row for the ground node is omitted.

With branch currents $\mathbf{i}_b$, branch voltages $\mathbf{v}_b$ and branch conductances collected in $\mathbf{G}_b=\operatorname{diag}(g_1,\dots,g_b)$, the three laws of circuit theory are three matrix equations:

$$
\begin{aligned}
\text{KCL:}&\quad \mathbf{A}\mathbf{i}_b=\mathbf{i}_s,\\
\text{KVL:}&\quad \mathbf{v}_b=\mathbf{A}^\mathsf{T}\mathbf{v},\\
\text{Ohm:}&\quad \mathbf{i}_b=\mathbf{G}_b\mathbf{v}_b.
\end{aligned}
$$

Substituting, $\mathbf{A}\mathbf{G}_b\mathbf{A}^\mathsf{T}\mathbf{v}=\mathbf{i}_s$, so

$$
\mathbf{G}=\mathbf{A}\mathbf{G}_b\mathbf{A}^\mathsf{T}.
$$

Symmetry is now automatic: $\mathbf{G}^\mathsf{T}=(\mathbf{A}^\mathsf{T})^\mathsf{T}\mathbf{G}_b^\mathsf{T}\mathbf{A}^\mathsf{T}=\mathbf{A}\mathbf{G}_b\mathbf{A}^\mathsf{T}=\mathbf{G}$. And positivity is automatic too: $\mathbf{v}^\mathsf{T}\mathbf{G}\mathbf{v}=\mathbf{v}_b^\mathsf{T}\mathbf{G}_b\mathbf{v}_b=\sum_k g_kv_{b,k}^2\ge0$.

> [!example] Worked Example 8.4 — Assembling the running example from its incidence matrix
> Branches $a,b,c,d,e,f$ in that order, directed from the first-named node to the second:
> $$
> \mathbf{A}=\begin{bmatrix}
> 1 & 1 & 0 & 0 & 0 & 1\\
> 0 & -1 & 1 & 1 & 0 & 0\\
> 0 & 0 & -1 & 0 & 1 & -1
> \end{bmatrix},\qquad \mathbf{G}_b=\operatorname{diag}(1,2,1,1,2,1).
> $$
> **Step 1.** $\mathbf{A}\mathbf{G}_b$ scales column $k$ of $\mathbf{A}$ by $g_k$ (right-multiplying by a diagonal matrix scales columns):
> $$
> \mathbf{A}\mathbf{G}_b=\begin{bmatrix}
> 1 & 2 & 0 & 0 & 0 & 1\\
> 0 & -2 & 1 & 1 & 0 & 0\\
> 0 & 0 & -1 & 0 & 2 & -1
> \end{bmatrix}.
> $$
> **Step 2.** Multiply by $\mathbf{A}^\mathsf{T}$ ($6\times 3$). Entry $(1,1)$ is row 1 of $\mathbf{A}\mathbf{G}_b$ dotted with row 1 of $\mathbf{A}$: $1+2+1=4$. Entry $(1,2)$: row 1 of $\mathbf{A}\mathbf{G}_b$ with row 2 of $\mathbf{A}$: $2\cdot(-1)=-2$. Entry $(1,3)$: $1\cdot(-1)=-1$. Entry $(2,2)$: $(-2)(-1)+1+1=4$. Entry $(2,3)$: $1\cdot(-1)=-1$. Entry $(3,3)$: $1+2+1=4$. By symmetry the rest follow:
> $$
> \mathbf{A}\mathbf{G}_b\mathbf{A}^\mathsf{T}=\begin{bmatrix} 4&-2&-1\\-2&4&-1\\-1&-1&4\end{bmatrix}=\mathbf{G}.\ \checkmark
> $$
> **Step 3.** Branch voltages: $\mathbf{v}_b=\mathbf{A}^\mathsf{T}\mathbf{v}=(4,\,1,\,1,\,3,\,2,\,2)^\mathsf{T}\ \mathrm{V}$. Branch currents $\mathbf{i}_b=\mathbf{G}_b\mathbf{v}_b=(4,2,1,3,4,2)^\mathsf{T}\ \mathrm{mA}$. KCL check: $\mathbf{A}\mathbf{i}_b=(4+2+2,\,-2+1+3,\,-1+4-2)^\mathsf{T}=(8,2,1)^\mathsf{T}$. ✓ Power: $\sum g_kv_{b,k}^2=16+2+1+9+8+4=40\ \mathrm{mW}$, matching $\mathbf{v}^\mathsf{T}\mathbf{i}_s$.

### Non-commutativity in practice: cascaded two-ports

A two-port is described by its **ABCD (transmission) matrix**, relating input voltage and current to output voltage and current:

$$
\begin{bmatrix}V_1\\ I_1\end{bmatrix}=\begin{bmatrix}A & B\\ C & D\end{bmatrix}\begin{bmatrix}V_2\\ I_2\end{bmatrix}.
$$

Because the output of one stage is the input of the next, a cascade is a product, taken in the physical order of the stages. A series impedance $Z$ has $\begin{bmatrix}1&Z\\0&1\end{bmatrix}$, a shunt admittance $Y$ has $\begin{bmatrix}1&0\\Y&1\end{bmatrix}$.

> [!example] Worked Example 8.5 — Order matters (an L-section both ways round)
> A series resistor $R=2\ \Omega$ is followed by a shunt conductance $G=0.5\ \mathrm{S}$. Find the cascade ABCD matrix, then repeat with the order reversed.
>
> **Series then shunt:**
> $$
> \begin{bmatrix}1&2\\0&1\end{bmatrix}\begin{bmatrix}1&0\\0.5&1\end{bmatrix}=\begin{bmatrix}1+1 & 0+2\\ 0+0.5 & 0+1\end{bmatrix}=\begin{bmatrix}2&2\\0.5&1\end{bmatrix}.
> $$
> **Shunt then series:**
> $$
> \begin{bmatrix}1&0\\0.5&1\end{bmatrix}\begin{bmatrix}1&2\\0&1\end{bmatrix}=\begin{bmatrix}1&2\\0.5&2\end{bmatrix}.
> $$
> **Reading.** The open-circuit voltage gain is $V_2/V_1=1/A$ (set $I_2=0$): $0.5$ in the first case (a divider) and $1$ in the second (the shunt element sits across the source and does nothing to $V_2$). Different physical circuits, different products. Both matrices have determinant $AD-BC=1$, the hallmark of a reciprocal two-port (§8.3).

### Block matrices

Partitioning matrices into blocks and multiplying the blocks as if they were scalars is legitimate provided the partitions are conformable and the order of block factors is kept:

$$
\begin{bmatrix}\mathbf{A}_{11}&\mathbf{A}_{12}\\ \mathbf{A}_{21}&\mathbf{A}_{22}\end{bmatrix}\begin{bmatrix}\mathbf{B}_{11}&\mathbf{B}_{12}\\ \mathbf{B}_{21}&\mathbf{B}_{22}\end{bmatrix}=\begin{bmatrix}\mathbf{A}_{11}\mathbf{B}_{11}+\mathbf{A}_{12}\mathbf{B}_{21} & \mathbf{A}_{11}\mathbf{B}_{12}+\mathbf{A}_{12}\mathbf{B}_{22}\\ \mathbf{A}_{21}\mathbf{B}_{11}+\mathbf{A}_{22}\mathbf{B}_{21} & \mathbf{A}_{21}\mathbf{B}_{12}+\mathbf{A}_{22}\mathbf{B}_{22}\end{bmatrix}.
$$

Block structure appears naturally in **modified nodal analysis** (MNA), where unknown node voltages $\mathbf{v}$ and unknown voltage-source currents $\mathbf{j}$ are stacked:

$$
\begin{bmatrix}\mathbf{G} & \mathbf{B}\\ \mathbf{B}^\mathsf{T} & \mathbf{0}\end{bmatrix}\begin{bmatrix}\mathbf{v}\\ \mathbf{j}\end{bmatrix}=\begin{bmatrix}\mathbf{i}_s\\ \mathbf{e}\end{bmatrix},
$$

where $\mathbf{B}$ says which nodes each voltage source connects and $\mathbf{e}$ holds the source voltages. Eliminating one block in terms of the other gives the **Schur complement**, a block version of Gaussian elimination (§8.7).

> [!example] Worked Example 8.6 — Matrix powers and a recursive network
> A ladder of identical L-sections from Worked Example 8.5 (series $2\ \Omega$ then shunt $0.5\ \mathrm{S}$) has $\mathbf{T}=\begin{bmatrix}2&2\\0.5&1\end{bmatrix}$. Find the ABCD matrix of three sections and the open-circuit voltage gain.
>
> **Step 1.** $\mathbf{T}^2=\begin{bmatrix}2&2\\0.5&1\end{bmatrix}\begin{bmatrix}2&2\\0.5&1\end{bmatrix}=\begin{bmatrix}4+1 & 4+2\\ 1+0.5 & 1+1\end{bmatrix}=\begin{bmatrix}5&6\\1.5&2\end{bmatrix}$.
>
> **Step 2.** $\mathbf{T}^3=\mathbf{T}^2\mathbf{T}=\begin{bmatrix}10+3 & 10+6\\ 3+1 & 3+2\end{bmatrix}=\begin{bmatrix}13&16\\4&5\end{bmatrix}$.
>
> **Step 3.** Gain $V_2/V_1=1/A=1/13\approx0.0769$. Check reciprocity: $\det\mathbf{T}^3=65-64=1$. ✓
>
> Powers of a matrix are computed far more cleverly via eigenvalues (Chapter 10); here the attenuation per section tends to the reciprocal of the larger eigenvalue of $\mathbf{T}$, which is $(3+\sqrt5)/2\approx2.618$.

### Cost of multiplication

Each entry of an $m\times p$ by $p\times n$ product costs $p$ multiplications and $p-1$ additions, so the whole product costs about $2mpn$ floating-point operations (**flops**). For two $n\times n$ matrices this is $2n^3$; a matrix–vector product is $2n^2$. These counts set the scale for everything in §§8.7–8.11: a dense $n=10^4$ matrix–matrix product is $2\times10^{12}$ flops, a couple of seconds on a modern CPU core with vector units, while a dense $n=10^6$ problem is out of reach — which is why sparsity (§8.11) matters.

Associativity lets us choose the cheapest bracketing. For $\mathbf{A}\in\mathbb{R}^{n\times n}$, $\mathbf{B}\in\mathbb{R}^{n\times n}$ and $\mathbf{x}\in\mathbb{R}^n$, $(\mathbf{A}\mathbf{B})\mathbf{x}$ costs $2n^3+2n^2$ but $\mathbf{A}(\mathbf{B}\mathbf{x})$ costs only $4n^2$.

> [!warning] Common Pitfall
> "Multiplying matrices entry by entry" (the **Hadamard product** $\mathbf{A}\odot\mathbf{B}$, MATLAB `.*`, NumPy `*`) is a different operation from the matrix product (MATLAB `*`, NumPy `@`). In NumPy, `A*x` with a 2-D array `A` and 1-D `x` silently broadcasts and returns a matrix, not $\mathbf{A}\mathbf{x}$. Always check the shape of the result.

> [!tip] Engineering Note
> Chain matrices (ABCD), scattering-transfer matrices (T-parameters) in RF design, Jones matrices in optics and the ray-transfer matrices of lens systems are all multiplied in the physical order in which the signal meets the elements. Write the product from left to right in the order the wave travels from input to output, and check the first and last factors against a sketch.

## 8.3 Determinants

### Motivation

The determinant is a single number computed from a square matrix that answers the most important yes/no question about it: *is the matrix invertible?* — equivalently, does $\mathbf{A}\mathbf{x}=\mathbf{b}$ have exactly one solution for every $\mathbf{b}$? It also measures how the linear map $\mathbf{x}\mapsto\mathbf{A}\mathbf{x}$ scales areas and volumes (the Jacobian determinant of Chapter 6 is exactly this, applied to the local linearisation of a nonlinear map), it gives closed-form solutions of small systems (Cramer's rule), and its vanishing defines the eigenvalues of Chapter 10 through $\det(\mathbf{A}-\lambda\mathbf{I})=0$. In circuit theory the determinant of the nodal or mesh matrix is the denominator of every network function, so its zeros as a function of $s$ are the natural frequencies (poles) of the network.

### Definition

For small matrices:

$$
\det\begin{bmatrix}a&b\\c&d\end{bmatrix}=ad-bc,
$$

$$
\det\begin{bmatrix}a_{11}&a_{12}&a_{13}\\a_{21}&a_{22}&a_{23}\\a_{31}&a_{32}&a_{33}\end{bmatrix}=a_{11}(a_{22}a_{33}-a_{23}a_{32})-a_{12}(a_{21}a_{33}-a_{23}a_{31})+a_{13}(a_{21}a_{32}-a_{22}a_{31}).
$$

The $3\times3$ formula is a pattern that generalises.

> [!definition] Minor, cofactor, determinant
> Let $\mathbf{A}$ be $n\times n$. The **minor** $M_{ij}$ is the determinant of the $(n-1)\times(n-1)$ matrix obtained by deleting row $i$ and column $j$. The **cofactor** is $C_{ij}=(-1)^{i+j}M_{ij}$. The **determinant** is defined recursively by $\det[a]=a$ for $1\times1$ and, for $n\ge2$, by **cofactor (Laplace) expansion along the first row**:
> $$
> \det\mathbf{A}=\sum_{j=1}^{n}a_{1j}C_{1j}.
> $$
> We also write $\lvert\mathbf{A}\rvert$ or $\Delta$.

The sign factor $(-1)^{i+j}$ follows a chessboard pattern beginning with $+$ in the top-left corner:

$$
\begin{bmatrix}+&-&+&\cdots\\-&+&-&\cdots\\+&-&+&\cdots\\\vdots&\vdots&\vdots&\ddots\end{bmatrix}.
$$

> [!theorem] Laplace expansion along any row or column
> For any fixed row $i$ or any fixed column $j$,
> $$
> \det\mathbf{A}=\sum_{k=1}^{n}a_{ik}C_{ik}=\sum_{k=1}^{n}a_{kj}C_{kj}.
> $$
> Moreover, expanding a row against the cofactors of a *different* row gives zero: $\sum_k a_{ik}C_{lk}=0$ for $i\ne l$.

The equivalent closed form, which makes the theorem transparent, is the **Leibniz formula**

$$
\det\mathbf{A}=\sum_{\sigma}\operatorname{sgn}(\sigma)\,a_{1\sigma(1)}a_{2\sigma(2)}\cdots a_{n\sigma(n)},
$$

summed over all $n!$ permutations $\sigma$ of $\{1,\dots,n\}$, with $\operatorname{sgn}\sigma=\pm1$ according to whether $\sigma$ is an even or odd number of swaps. Each term takes exactly one entry from every row and every column. Expanding along any row or column just groups these $n!$ terms differently. The "different row gives zero" statement holds because $\sum_k a_{ik}C_{lk}$ is the determinant of a matrix whose row $l$ has been replaced by a copy of row $i$ — a matrix with two equal rows, whose determinant is zero (property 3 below).

The $n!$ terms make the definition useless for computation beyond $n\approx 5$: for $n=25$, $25!\approx1.6\times10^{25}$ products would take longer than the age of the universe at $10^{9}$ products per second. Row reduction (below) needs only about $\tfrac23n^3$ flops.

### Properties

> [!theorem] Properties of determinants
> For $n\times n$ matrices $\mathbf{A},\mathbf{B}$ and scalar $c$:
> 1. $\det\mathbf{A}^\mathsf{T}=\det\mathbf{A}$; $\det\mathbf{A}^\mathsf{H}=(\det\mathbf{A})^*$.
> 2. **Row swap** changes the sign of the determinant.
> 3. If two rows (or columns) are equal, or a row is all zeros, $\det\mathbf{A}=0$.
> 4. **Row scaling**: multiplying one row by $c$ multiplies the determinant by $c$. Hence $\det(c\mathbf{A})=c^n\det\mathbf{A}$.
> 5. **Row replacement**: adding a multiple of one row to another leaves the determinant unchanged.
> 6. The determinant of a triangular (in particular diagonal) matrix is the product of its diagonal entries.
> 7. **Product rule**: $\det(\mathbf{A}\mathbf{B})=\det\mathbf{A}\,\det\mathbf{B}$.
> 8. $\mathbf{A}$ is invertible if and only if $\det\mathbf{A}\neq0$, and then $\det(\mathbf{A}^{-1})=1/\det\mathbf{A}$.
> 9. The determinant is linear in each row separately (multilinearity).

*Proof sketches.* (1) Transposing maps each Leibniz term to a term of the same sign. (2) Swapping two rows composes every permutation with a transposition, flipping every sign. (3) Swapping two equal rows leaves $\mathbf{A}$ unchanged but negates $\det\mathbf{A}$, so $\det\mathbf{A}=-\det\mathbf{A}=0$; a zero row puts a zero factor in every Leibniz term. (4) and (9) Each Leibniz term contains exactly one factor from each row. (5) By multilinearity, $\det(\dots,\mathbf{r}_i+c\mathbf{r}_l,\dots)=\det(\dots,\mathbf{r}_i,\dots)+c\det(\dots,\mathbf{r}_l,\dots)$, and the second determinant has two equal rows. (6) Expand repeatedly down the first column: only the diagonal entry survives each time. (7) Row-reduce $\mathbf{A}$ with elementary matrices whose determinants are known (§8.4), then use (6). (8) follows from (7), $\det\mathbf{I}=1$, and the fact that row reduction of a singular matrix produces a zero row. $\square$

Properties 2, 4 and 5 are exactly the three row operations of Gaussian elimination. That gives the practical algorithm:

> **Determinant by elimination.** Reduce $\mathbf{A}$ to upper triangular $\mathbf{U}$ using row replacements and $s$ row swaps (no scaling). Then $\det\mathbf{A}=(-1)^s\,u_{11}u_{22}\cdots u_{nn}$, the signed product of the pivots.

### Geometric meaning

For a $2\times2$ matrix with columns $\mathbf{a}_1=(a,c)^\mathsf{T}$ and $\mathbf{a}_2=(b,d)^\mathsf{T}$, $\lvert ad-bc\rvert$ is the area of the parallelogram spanned by the columns; the sign records orientation (positive if $\mathbf{a}_1$ turns anticlockwise to reach $\mathbf{a}_2$). For $3\times3$, $\lvert\det\mathbf{A}\rvert$ is the volume of the parallelepiped spanned by the columns — the scalar triple product $\mathbf{a}_1\cdot(\mathbf{a}_2\times\mathbf{a}_3)$ of Chapter 7. In general the linear map $\mathbf{x}\mapsto\mathbf{A}\mathbf{x}$ multiplies every $n$-dimensional volume by $\lvert\det\mathbf{A}\rvert$. A zero determinant means the columns are squashed into a lower-dimensional set — they are linearly dependent — and information is lost; no inverse map can recover it. This is the geometric content of property 8, and it also explains the product rule: applying $\mathbf{B}$ then $\mathbf{A}$ scales volumes by $\lvert\det\mathbf{B}\rvert$ and then by $\lvert\det\mathbf{A}\rvert$.

### Cramer's rule

> [!theorem] Cramer's rule
> If $\det\mathbf{A}\ne0$, the unique solution of $\mathbf{A}\mathbf{x}=\mathbf{b}$ is
> $$
> x_k=\frac{\det\mathbf{A}_k}{\det\mathbf{A}},\qquad k=1,\dots,n,
> $$
> where $\mathbf{A}_k$ is $\mathbf{A}$ with column $k$ replaced by $\mathbf{b}$.

*Proof.* Write $\mathbf{b}=\sum_j x_j\mathbf{a}_j$ (column view). Replace column $k$ of $\mathbf{A}$ by this sum and use multilinearity in columns: $\det\mathbf{A}_k=\sum_j x_j\det(\mathbf{a}_1,\dots,\mathbf{a}_j,\dots,\mathbf{a}_n)$ with $\mathbf{a}_j$ in slot $k$. Every term with $j\ne k$ has a repeated column and vanishes, leaving $\det\mathbf{A}_k=x_k\det\mathbf{A}$. $\square$

Cramer's rule is the right tool for $2\times2$ and $3\times3$ systems with symbolic entries — for example, to obtain a transfer function in terms of $s$ and component values — and the wrong tool for numerical systems larger than that, where it costs $n+1$ determinants and is numerically inferior to elimination.

> [!example] Worked Example 8.7 — Determinant of the running example, two ways
> Find $\det\mathbf{G}$ for $\mathbf{G}=\begin{bmatrix}4&-2&-1\\-2&4&-1\\-1&-1&4\end{bmatrix}$.
>
> **Method 1: cofactor expansion along row 1.**
> $$
> \begin{aligned}
> \det\mathbf{G}&=4\begin{vmatrix}4&-1\\-1&4\end{vmatrix}-(-2)\begin{vmatrix}-2&-1\\-1&4\end{vmatrix}+(-1)\begin{vmatrix}-2&4\\-1&-1\end{vmatrix}\\
> &=4(16-1)+2(-8-1)-(2+4)=60-18-6=36.
> \end{aligned}
> $$
> **Method 2: elimination.** $R_2\leftarrow R_2+\tfrac12R_1$ and $R_3\leftarrow R_3+\tfrac14R_1$ give rows $(0,3,-\tfrac32)$ and $(0,-\tfrac32,\tfrac{15}4)$. Then $R_3\leftarrow R_3+\tfrac12R_2$ gives $(0,0,3)$. No swaps, so
> $$
> \det\mathbf{G}=4\cdot3\cdot3=36.
> $$
> The units are $\mathrm{mS}^3$. Method 2 is the one that scales.

> [!example] Worked Example 8.8 — A determinant needing a row swap
> The modified nodal formulation of §8.7 (Worked Example 8.22) leads to
> $$
> \mathbf{M}=\begin{bmatrix}0&1&0&0\\1&4&-2&-1\\0&-2&4&-1\\0&-1&-1&4\end{bmatrix}.
> $$
> Find $\det\mathbf{M}$.
>
> **Step 1.** $m_{11}=0$, so swap $R_1\leftrightarrow R_2$ ($s=1$):
> $$
> \begin{bmatrix}1&4&-2&-1\\0&1&0&0\\0&-2&4&-1\\0&-1&-1&4\end{bmatrix}.
> $$
> **Step 2.** $R_3\leftarrow R_3+2R_2$, $R_4\leftarrow R_4+R_2$: rows 3 and 4 become $(0,0,4,-1)$ and $(0,0,-1,4)$.
>
> **Step 3.** $R_4\leftarrow R_4+\tfrac14R_3$: row 4 becomes $(0,0,0,\tfrac{15}4)$.
>
> **Step 4.** $\det\mathbf{M}=(-1)^1\cdot1\cdot1\cdot4\cdot\tfrac{15}{4}=-15$.
>
> **Check by expansion along row 1:** only $m_{12}=1$ is nonzero, so $\det\mathbf{M}=(-1)^{1+2}M_{12}$, and $M_{12}=\det\begin{bmatrix}1&-2&-1\\0&4&-1\\0&-1&4\end{bmatrix}=1\cdot(16-1)=15$. Hence $-15$. ✓ (Choosing the sparsest row for the expansion is always worth a glance.)

### Two-mesh example and its geometry

Consider a two-mesh resistive circuit: mesh 1 contains a $2\ \mathrm{k\Omega}$ resistor and a $7\ \mathrm{V}$ source, mesh 2 contains a $4\ \mathrm{k\Omega}$ resistor and a $1\ \mathrm{V}$ source, and the two meshes share a $3\ \mathrm{k\Omega}$ resistor. With clockwise mesh currents $i_1,i_2$ in mA, KVL gives the **mesh equations** $\mathbf{R}\mathbf{i}=\mathbf{v}_s$:

$$
\begin{bmatrix}5&-3\\-3&7\end{bmatrix}\begin{bmatrix}i_1\\i_2\end{bmatrix}=\begin{bmatrix}7\\1\end{bmatrix}.
$$

The diagonal entries are the total resistance around each mesh, the off-diagonal entries minus the shared resistance — the mesh-analysis mirror image of the nodal rule. Units: $\mathrm{k\Omega}\times\mathrm{mA}=\mathrm{V}$.

> [!example] Worked Example 8.9 — Cramer's rule for the two-mesh circuit
> **Step 1.** $\det\mathbf{R}=5\cdot7-(-3)(-3)=35-9=26\ \mathrm{k\Omega}^2$.
>
> **Step 2.**
> $$
> i_1=\frac{\begin{vmatrix}7&-3\\1&7\end{vmatrix}}{26}=\frac{49+3}{26}=2\ \mathrm{mA},\qquad
> i_2=\frac{\begin{vmatrix}5&7\\-3&1\end{vmatrix}}{26}=\frac{5+21}{26}=1\ \mathrm{mA}.
> $$
> **Step 3 (check).** $5(2)-3(1)=7$ ✓, $-3(2)+7(1)=1$ ✓. The shared resistor carries $i_1-i_2=1\ \mathrm{mA}$.
>
> **Geometry.** Each equation is a straight line in the $(i_1,i_2)$ plane; the solution is their intersection. Because $\det\mathbf{R}\ne0$ the lines are not parallel and meet in exactly one point.

```desmos-graph
left=-1; right=4; bottom=-1; top=3
---
5x-3y=7|blue
-3x+7y=1|red
(2,1)|label:(2, 1)
```

*The two mesh equations as lines in the $(i_1,i_2)$ plane (horizontal axis $i_1$, vertical axis $i_2$, both in mA). They cross at the unique solution $(2,1)$; a nonzero determinant guarantees the slopes $5/3$ and $3/7$ differ.*

> [!example] Worked Example 8.10 — Cramer's rule for one node voltage, symbolically
> In the running example, suppose only node 1 is driven, by a current $I$, and we want the transfer resistance $v_3/I$.
>
> **Step 1.** $\mathbf{i}_s=(I,0,0)^\mathsf{T}$. Replace column 3 of $\mathbf{G}$:
> $$
> \det\mathbf{G}_3=\begin{vmatrix}4&-2&I\\-2&4&0\\-1&-1&0\end{vmatrix}=I\begin{vmatrix}-2&4\\-1&-1\end{vmatrix}=I(2+4)=6I,
> $$
> expanding along column 3, where only the top entry is nonzero and its sign factor is $(-1)^{1+3}=+1$.
>
> **Step 2.** $v_3=6I/36=I/6$, so the transfer resistance is $v_3/I=1/6\ \mathrm{k\Omega}\approx167\ \Omega$ (since $1/\mathrm{mS}=\mathrm{k\Omega}$).
>
> **Step 3 (full solution check).** With the original sources $(8,2,1)$, Cramer gives $\det\mathbf{G}_1=144$, $\det\mathbf{G}_2=108$, $\det\mathbf{G}_3=72$, so $\mathbf{v}=(144,108,72)/36=(4,3,2)\ \mathrm{V}$. ✓

### Determinant as a function of a parameter

When a component value varies, the determinant becomes a polynomial in that value, and its zeros mark the values at which the circuit equations become singular. Consider $\mathbf{A}-\lambda\mathbf{I}$ with $\mathbf{A}=\begin{bmatrix}4&-2\\-2&4\end{bmatrix}$ (the nodal matrix of nodes 1 and 2 of a circuit like the running example with node 3 tied to ground):

$$
p(\lambda)=\det\begin{bmatrix}4-\lambda&-2\\-2&4-\lambda\end{bmatrix}=(4-\lambda)^2-4=\lambda^2-8\lambda+12=(\lambda-2)(\lambda-6).
$$

```desmos-graph
left=-1; right=8; bottom=-6; top=14
---
y=(4-x)^2-4|blue
(2,0)|label:λ = 2
(6,0)|label:λ = 6
(0,12)|label:det A = 12
```

*$\det(\mathbf{A}-\lambda\mathbf{I})$ as a function of $\lambda$ (horizontal axis). At $\lambda=0$ it equals $\det\mathbf{A}=12$; it vanishes at $\lambda=2$ and $\lambda=6$, where $\mathbf{A}-\lambda\mathbf{I}$ is singular. These are the eigenvalues of $\mathbf{A}$ (Chapter 10).*

> [!warning] Common Pitfall
> A small determinant does **not** mean a matrix is nearly singular. $\det(10^{-3}\mathbf{I}_{10})=10^{-30}$, yet $10^{-3}\mathbf{I}$ is perfectly conditioned; and the determinant changes with units ($\mathbf{G}$ in $\mathrm{S}$ instead of $\mathrm{mS}$ has a determinant $10^{9}$ times smaller, with no change in the difficulty of solving the system). Use the condition number (§8.10), not the determinant, to judge numerical near-singularity.

> [!tip] Engineering Note
> For a network with impedances depending on $s$, $\Delta(s)=\det\mathbf{Y}(s)$ is, up to a factor, the characteristic polynomial of the network; its roots are the natural frequencies. Every node voltage is $\Delta_k(s)/\Delta(s)$ by Cramer's rule, which is why all transfer functions of a network share the same denominator. Mason's gain formula for signal-flow graphs is a combinatorial way of evaluating these same determinants.

## 8.4 Matrix inverse

### Motivation

For a scalar equation $ax=b$ we write $x=a^{-1}b$ provided $a\ne0$. The matrix analogue $\mathbf{x}=\mathbf{A}^{-1}\mathbf{b}$ is the natural way to *express* the solution of a linear system and to reason about it, even though (as we shall stress) it is rarely the way to *compute* it. The inverse also has direct physical meaning: the inverse of a nodal conductance matrix is the **nodal resistance (impedance) matrix**, whose entries are driving-point and transfer resistances; the inverse of a mesh resistance matrix is a conductance matrix; the inverse of a two-port's $\mathbf{Z}$-parameter matrix is its $\mathbf{Y}$-parameter matrix.

### Definition and uniqueness

> [!definition] Inverse
> A square matrix $\mathbf{A}$ is **invertible** (nonsingular) if there is a matrix $\mathbf{A}^{-1}$ with
> $$
> \mathbf{A}\mathbf{A}^{-1}=\mathbf{A}^{-1}\mathbf{A}=\mathbf{I}.
> $$
> A square matrix with no inverse is **singular**.

The inverse, when it exists, is unique: if $\mathbf{B}$ and $\mathbf{C}$ are both inverses then $\mathbf{B}=\mathbf{B}\mathbf{I}=\mathbf{B}(\mathbf{A}\mathbf{C})=(\mathbf{B}\mathbf{A})\mathbf{C}=\mathbf{I}\mathbf{C}=\mathbf{C}$. For square matrices a one-sided inverse is automatically two-sided: if $\mathbf{B}\mathbf{A}=\mathbf{I}$ then $\mathbf{A}\mathbf{B}=\mathbf{I}$ too (this follows from the rank arguments of §8.5). Non-square matrices never have two-sided inverses; their one-sided and generalised inverses (the pseudo-inverse) belong to least-squares theory in Chapter 9.

If $\mathbf{A}$ is invertible, $\mathbf{A}\mathbf{x}=\mathbf{b}$ has the unique solution $\mathbf{x}=\mathbf{A}^{-1}\mathbf{b}$: it is a solution because $\mathbf{A}(\mathbf{A}^{-1}\mathbf{b})=\mathbf{b}$, and it is the only one because any solution satisfies $\mathbf{x}=\mathbf{A}^{-1}\mathbf{A}\mathbf{x}=\mathbf{A}^{-1}\mathbf{b}$.

### Formulas

For a $2\times2$ matrix the inverse is: *swap the diagonal, negate the off-diagonal, divide by the determinant*,

$$
\begin{bmatrix}a&b\\c&d\end{bmatrix}^{-1}=\frac{1}{ad-bc}\begin{bmatrix}d&-b\\-c&a\end{bmatrix},\qquad ad-bc\ne0.
$$

In general:

> [!theorem] Adjugate formula
> Let $\mathbf{C}=[C_{ij}]$ be the matrix of cofactors of $\mathbf{A}$ and define the **adjugate** $\operatorname{adj}\mathbf{A}=\mathbf{C}^\mathsf{T}$. Then
> $$
> \mathbf{A}\,\operatorname{adj}\mathbf{A}=\operatorname{adj}\mathbf{A}\,\mathbf{A}=(\det\mathbf{A})\,\mathbf{I},
> $$
> so if $\det\mathbf{A}\ne0$,
> $$
> \mathbf{A}^{-1}=\frac{1}{\det\mathbf{A}}\operatorname{adj}\mathbf{A}.
> $$

*Proof.* The $(i,l)$ entry of $\mathbf{A}\operatorname{adj}\mathbf{A}$ is $\sum_k a_{ik}C_{lk}$. By the Laplace expansion theorem of §8.3 this equals $\det\mathbf{A}$ when $i=l$ and $0$ when $i\ne l$. $\square$

The adjugate formula is the matrix version of Cramer's rule — indeed $(\operatorname{adj}\mathbf{A}\,\mathbf{b})_k=\det\mathbf{A}_k$ — and has the same uses (symbolic $2\times2$ and $3\times3$ work, theory) and the same drawback (cost) for anything larger.

> [!theorem] Algebra of inverses
> If $\mathbf{A}$ and $\mathbf{B}$ are invertible $n\times n$ matrices and $c\ne0$:
> 1. $(\mathbf{A}^{-1})^{-1}=\mathbf{A}$.
> 2. $(\mathbf{A}\mathbf{B})^{-1}=\mathbf{B}^{-1}\mathbf{A}^{-1}$ (order reverses — "socks and shoes").
> 3. $(\mathbf{A}^\mathsf{T})^{-1}=(\mathbf{A}^{-1})^\mathsf{T}$ and $(\mathbf{A}^\mathsf{H})^{-1}=(\mathbf{A}^{-1})^\mathsf{H}$; hence the inverse of a symmetric matrix is symmetric.
> 4. $(c\mathbf{A})^{-1}=c^{-1}\mathbf{A}^{-1}$.
> 5. $\mathbf{A}^{-1}$ of a diagonal matrix is the diagonal matrix of reciprocals; the inverse of an upper (lower) triangular matrix is upper (lower) triangular.

*Proof of 2 and 3.* $(\mathbf{A}\mathbf{B})(\mathbf{B}^{-1}\mathbf{A}^{-1})=\mathbf{A}\mathbf{I}\mathbf{A}^{-1}=\mathbf{I}$. Transposing $\mathbf{A}\mathbf{A}^{-1}=\mathbf{I}$ gives $(\mathbf{A}^{-1})^\mathsf{T}\mathbf{A}^\mathsf{T}=\mathbf{I}$. $\square$

### Elementary matrices

Each row operation of elimination is a left-multiplication by an **elementary matrix**, obtained by doing that operation to $\mathbf{I}$:

| Row operation | Elementary matrix ($3\times3$ example) | Inverse | Determinant |
|---|---|---|---|
| $R_2\leftarrow R_2+mR_1$ | $\mathbf{E}=\begin{bmatrix}1&0&0\\m&1&0\\0&0&1\end{bmatrix}$ | same with $-m$ | $1$ |
| $R_2\leftarrow cR_2$ | $\operatorname{diag}(1,c,1)$ | $\operatorname{diag}(1,1/c,1)$ | $c$ |
| $R_1\leftrightarrow R_2$ | $\mathbf{P}=\begin{bmatrix}0&1&0\\1&0&0\\0&0&1\end{bmatrix}$ | $\mathbf{P}$ itself | $-1$ |

Every elementary matrix is invertible, and its inverse undoes the operation. If a sequence of row operations reduces $\mathbf{A}$ to $\mathbf{I}$, say $\mathbf{E}_k\cdots\mathbf{E}_2\mathbf{E}_1\mathbf{A}=\mathbf{I}$, then $\mathbf{A}^{-1}=\mathbf{E}_k\cdots\mathbf{E}_1$ — the same operations applied to $\mathbf{I}$. This is the Gauss–Jordan method for the inverse (§8.8), and it also proves the product rule for determinants claimed in §8.3.

### The invertible matrix theorem

> [!theorem] Invertible matrix theorem
> For an $n\times n$ matrix $\mathbf{A}$ the following are equivalent:
> 1. $\mathbf{A}$ is invertible.
> 2. $\det\mathbf{A}\ne0$.
> 3. $\mathbf{A}\mathbf{x}=\mathbf{b}$ has a unique solution for every $\mathbf{b}$.
> 4. $\mathbf{A}\mathbf{x}=\mathbf{0}$ has only the trivial solution $\mathbf{x}=\mathbf{0}$.
> 5. $\operatorname{rank}\mathbf{A}=n$ (§8.5); the columns (and rows) are linearly independent.
> 6. Gaussian elimination produces $n$ nonzero pivots (with row swaps if needed).
> 7. $\mathbf{A}$ is a product of elementary matrices.
> 8. $0$ is not an eigenvalue of $\mathbf{A}$ (Chapter 10).

For circuits, statement 4 has a direct reading: the nodal equations with all independent sources switched off have only the zero solution. A network in which some node "floats" (no DC path to ground) violates it — the floating part can sit at any voltage — and its nodal matrix is singular (§8.5, Worked Example 8.17).

> [!example] Worked Example 8.11 — Inverse of the mesh matrix
> For $\mathbf{R}=\begin{bmatrix}5&-3\\-3&7\end{bmatrix}\ \mathrm{k\Omega}$ (§8.3), find $\mathbf{R}^{-1}$ and use it to solve $\mathbf{R}\mathbf{i}=(7,1)^\mathsf{T}$.
>
> **Step 1.** $\det\mathbf{R}=26$. Swap diagonal, negate off-diagonal:
> $$
> \mathbf{R}^{-1}=\frac{1}{26}\begin{bmatrix}7&3\\3&5\end{bmatrix}\ \mathrm{mS}.
> $$
> **Step 2.** $\mathbf{i}=\mathbf{R}^{-1}\mathbf{v}_s=\tfrac1{26}(49+3,\ 21+5)^\mathsf{T}=(2,1)^\mathsf{T}\ \mathrm{mA}$. ✓
>
> **Step 3 (check).** $\mathbf{R}\mathbf{R}^{-1}=\tfrac1{26}\begin{bmatrix}35-9&15-15\\-21+21&-9+35\end{bmatrix}=\mathbf{I}$. ✓
>
> **Reading.** $[\mathbf{R}^{-1}]_{12}=3/26\ \mathrm{mS}$ is the current in mesh 1 per volt of source in mesh 2. Its equality with $[\mathbf{R}^{-1}]_{21}$ is reciprocity, and it follows from property 3 of the inverse.

> [!example] Worked Example 8.12 — The nodal resistance matrix of the running example
> Find $\mathbf{Z}=\mathbf{G}^{-1}$ by the adjugate formula and interpret its entries.
>
> **Step 1 (cofactors).** Using $\mathbf{G}=\begin{bmatrix}4&-2&-1\\-2&4&-1\\-1&-1&4\end{bmatrix}$:
> $$
> \begin{aligned}
> C_{11}&=+\begin{vmatrix}4&-1\\-1&4\end{vmatrix}=15, & C_{12}&=-\begin{vmatrix}-2&-1\\-1&4\end{vmatrix}=9, & C_{13}&=+\begin{vmatrix}-2&4\\-1&-1\end{vmatrix}=6,\\
> C_{22}&=+\begin{vmatrix}4&-1\\-1&4\end{vmatrix}=15, & C_{23}&=-\begin{vmatrix}4&-2\\-1&-1\end{vmatrix}=6, & C_{33}&=+\begin{vmatrix}4&-2\\-2&4\end{vmatrix}=12,
> \end{aligned}
> $$
> and by symmetry $C_{21}=9$, $C_{31}=6$, $C_{32}=6$.
>
> **Step 2.** $\det\mathbf{G}=36$ (Worked Example 8.7), so
> $$
> \mathbf{Z}=\mathbf{G}^{-1}=\frac{1}{36}\begin{bmatrix}15&9&6\\9&15&6\\6&6&12\end{bmatrix}=\begin{bmatrix}\tfrac{5}{12}&\tfrac14&\tfrac16\\\tfrac14&\tfrac{5}{12}&\tfrac16\\\tfrac16&\tfrac16&\tfrac13\end{bmatrix}\ \mathrm{k\Omega}.
> $$
> **Step 3 (check).** Row 1 of $\mathbf{G}$ times column 1 of $\operatorname{adj}\mathbf{G}$: $4(15)-2(9)-1(6)=36$ ✓; row 1 times column 2: $4(9)-2(15)-6=0$ ✓.
>
> **Step 4 (solve).** $\mathbf{v}=\mathbf{Z}\mathbf{i}_s=\tfrac1{36}(120+18+6,\ 72+30+6,\ 48+12+12)^\mathsf{T}=(4,3,2)^\mathsf{T}\ \mathrm{V}$. ✓
>
> **Reading.** $Z_{11}=5/12\ \mathrm{k\Omega}\approx417\ \Omega$ is the resistance seen looking into node 1 from ground (the Thévenin resistance at node 1); $Z_{31}=1/6\ \mathrm{k\Omega}$ is the transfer resistance found by Cramer's rule in Worked Example 8.10. All entries are positive and each diagonal entry is the largest in its column — typical of resistive networks: injecting current at a node raises that node most.

### Updating an inverse: Sherman–Morrison

Design work often changes one component and re-solves. Changing a conductance $g$ between nodes $p$ and $q$ by $\delta$ changes $\mathbf{G}$ by a **rank-one** matrix $\delta\,\mathbf{u}\mathbf{u}^\mathsf{T}$ with $\mathbf{u}=\mathbf{e}_p-\mathbf{e}_q$ (or $\mathbf{u}=\mathbf{e}_p$ for a branch to ground), where $\mathbf{e}_k$ is column $k$ of $\mathbf{I}$.

> [!theorem] Sherman–Morrison formula
> If $\mathbf{A}$ is invertible and $1+\mathbf{v}^\mathsf{T}\mathbf{A}^{-1}\mathbf{u}\ne0$, then
> $$
> (\mathbf{A}+\mathbf{u}\mathbf{v}^\mathsf{T})^{-1}=\mathbf{A}^{-1}-\frac{\mathbf{A}^{-1}\mathbf{u}\,\mathbf{v}^\mathsf{T}\mathbf{A}^{-1}}{1+\mathbf{v}^\mathsf{T}\mathbf{A}^{-1}\mathbf{u}}.
> $$
> Consequently the new solution of $(\mathbf{A}+\mathbf{u}\mathbf{v}^\mathsf{T})\mathbf{x}'=\mathbf{b}$ is
> $$
> \mathbf{x}'=\mathbf{x}-\frac{\mathbf{v}^\mathsf{T}\mathbf{x}}{1+\mathbf{v}^\mathsf{T}\mathbf{w}}\,\mathbf{w},\qquad \mathbf{x}=\mathbf{A}^{-1}\mathbf{b},\ \ \mathbf{w}=\mathbf{A}^{-1}\mathbf{u}.
> $$

*Proof.* Multiply $\mathbf{A}+\mathbf{u}\mathbf{v}^\mathsf{T}$ by the claimed inverse. With $\beta=\mathbf{v}^\mathsf{T}\mathbf{A}^{-1}\mathbf{u}$ (a scalar),
$$
(\mathbf{A}+\mathbf{u}\mathbf{v}^\mathsf{T})\Big(\mathbf{A}^{-1}-\frac{\mathbf{A}^{-1}\mathbf{u}\mathbf{v}^\mathsf{T}\mathbf{A}^{-1}}{1+\beta}\Big)=\mathbf{I}+\mathbf{u}\mathbf{v}^\mathsf{T}\mathbf{A}^{-1}-\frac{\mathbf{u}\mathbf{v}^\mathsf{T}\mathbf{A}^{-1}+\beta\,\mathbf{u}\mathbf{v}^\mathsf{T}\mathbf{A}^{-1}}{1+\beta}=\mathbf{I}.\ \square
$$

The update costs $O(n^2)$ given $\mathbf{w}$ (one extra solve with the existing factorisation), instead of $O(n^3)$ for a fresh solve. Its generalisation to rank-$k$ changes is the **Woodbury identity** $(\mathbf{A}+\mathbf{U}\mathbf{C}\mathbf{V}^\mathsf{T})^{-1}=\mathbf{A}^{-1}-\mathbf{A}^{-1}\mathbf{U}(\mathbf{C}^{-1}+\mathbf{V}^\mathsf{T}\mathbf{A}^{-1}\mathbf{U})^{-1}\mathbf{V}^\mathsf{T}\mathbf{A}^{-1}$, used in Kalman filtering and in contingency analysis of power networks (removing a line is a rank-one change with $\delta=-g$).

> [!example] Worked Example 8.13 — Doubling one resistor's conductance
> In the running example, branch $f$ (between nodes 1 and 3) is changed from $1\ \mathrm{mS}$ to $2\ \mathrm{mS}$. Find the new node voltages without re-solving.
>
> **Step 1.** $\delta=1$, $\mathbf{u}=\mathbf{v}=\mathbf{e}_1-\mathbf{e}_3=(1,0,-1)^\mathsf{T}$, so $\mathbf{G}'=\mathbf{G}+\mathbf{u}\mathbf{u}^\mathsf{T}$.
>
> **Step 2.** $\mathbf{w}=\mathbf{Z}\mathbf{u}=$ column 1 minus column 3 of $\mathbf{Z}$ $=(\tfrac5{12}-\tfrac16,\ \tfrac14-\tfrac16,\ \tfrac16-\tfrac13)^\mathsf{T}=(\tfrac14,\ \tfrac1{12},\ -\tfrac16)^\mathsf{T}$.
>
> **Step 3.** $\mathbf{u}^\mathsf{T}\mathbf{w}=\tfrac14+\tfrac16=\tfrac5{12}$ — this is the resistance (in $\mathrm{k\Omega}$) seen by the branch between nodes 1 and 3 — and $\mathbf{u}^\mathsf{T}\mathbf{v}_{\text{old}}=4-2=2\ \mathrm{V}$.
>
> **Step 4.**
> $$
> \mathbf{v}'=\begin{bmatrix}4\\3\\2\end{bmatrix}-\frac{2}{1+\tfrac5{12}}\begin{bmatrix}\tfrac14\\\tfrac1{12}\\-\tfrac16\end{bmatrix}=\begin{bmatrix}4\\3\\2\end{bmatrix}-\frac{24}{17}\begin{bmatrix}\tfrac14\\\tfrac1{12}\\-\tfrac16\end{bmatrix}=\frac1{17}\begin{bmatrix}62\\49\\38\end{bmatrix}\approx\begin{bmatrix}3.647\\2.882\\2.235\end{bmatrix}\ \mathrm{V}.
> $$
> **Step 5 (check).** $\mathbf{G}'=\begin{bmatrix}5&-2&-2\\-2&4&-1\\-2&-1&5\end{bmatrix}$; row 1 times $(62,49,38)/17$ is $(310-98-76)/17=136/17=8$ ✓; row 3: $(-124-49+190)/17=17/17=1$ ✓.
>
> **Reading.** Stiffening the link between nodes 1 and 3 pulls their voltages together: the difference falls from $2\ \mathrm{V}$ to $24/17\approx1.41\ \mathrm{V}$.

### Block inverses, Schur complements and Kron reduction

Partition a system into the unknowns we care about ($\mathbf{x}_1$) and those we would like to eliminate ($\mathbf{x}_2$):

$$
\begin{bmatrix}\mathbf{A}_{11}&\mathbf{A}_{12}\\\mathbf{A}_{21}&\mathbf{A}_{22}\end{bmatrix}\begin{bmatrix}\mathbf{x}_1\\\mathbf{x}_2\end{bmatrix}=\begin{bmatrix}\mathbf{b}_1\\\mathbf{b}_2\end{bmatrix}.
$$

If $\mathbf{A}_{22}$ is invertible, the second block row gives $\mathbf{x}_2=\mathbf{A}_{22}^{-1}(\mathbf{b}_2-\mathbf{A}_{21}\mathbf{x}_1)$; substituting into the first,

$$
\underbrace{(\mathbf{A}_{11}-\mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{A}_{21})}_{\mathbf{S},\ \text{Schur complement}}\mathbf{x}_1=\mathbf{b}_1-\mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{b}_2.
$$

Moreover $\det\mathbf{A}=\det\mathbf{A}_{22}\,\det\mathbf{S}$. In network theory this is **Kron reduction**: eliminating internal nodes leaves an equivalent network seen from the retained nodes. The star–delta transformation is the special case of eliminating a single node with three neighbours.

> [!example] Worked Example 8.14 — Kron reduction of node 3
> Eliminate $v_3$ from the running example and solve the reduced $2\times2$ system.
>
> **Step 1.** $\mathbf{A}_{11}=\begin{bmatrix}4&-2\\-2&4\end{bmatrix}$, $\mathbf{A}_{12}=\begin{bmatrix}-1\\-1\end{bmatrix}=\mathbf{A}_{21}^\mathsf{T}$, $\mathbf{A}_{22}=[4]$, $\mathbf{b}_1=(8,2)^\mathsf{T}$, $b_2=1$.
>
> **Step 2.**
> $$
> \mathbf{S}=\begin{bmatrix}4&-2\\-2&4\end{bmatrix}-\frac14\begin{bmatrix}-1\\-1\end{bmatrix}\begin{bmatrix}-1&-1\end{bmatrix}=\begin{bmatrix}3.75&-2.25\\-2.25&3.75\end{bmatrix},\qquad
> \mathbf{b}_1'=\begin{bmatrix}8\\2\end{bmatrix}-\frac14\begin{bmatrix}-1\\-1\end{bmatrix}(1)=\begin{bmatrix}8.25\\2.25\end{bmatrix}.
> $$
> **Step 3.** $\det\mathbf{S}=3.75^2-2.25^2=9$, and $\det\mathbf{A}_{22}\det\mathbf{S}=4\cdot9=36=\det\mathbf{G}$ ✓. Then
> $$
> v_1=\frac{8.25(3.75)+2.25(2.25)}{9}=\frac{36}{9}=4,\qquad v_2=\frac{3.75(2.25)+2.25(8.25)}{9}=\frac{27}{9}=3\ \mathrm{V}.
> $$
> **Step 4.** Back-substitute: $v_3=\tfrac14(1+v_1+v_2)=\tfrac14(8)=2\ \mathrm{V}$. ✓
>
> **Reading.** $\mathbf{S}$ is the nodal matrix of an equivalent two-node network: the off-diagonal $-2.25$ says the direct $2\ \mathrm{mS}$ link between nodes 1 and 2 is supplemented by an effective $0.25\ \mathrm{mS}$ through the eliminated node (the series combination of the two $1\ \mathrm{mS}$ branches meeting there, weighted by node 3's $4\ \mathrm{mS}$ total — exactly the star–delta rule $g_1g_2/\sum g$). The modified source vector shows that node 3's $1\ \mathrm{mA}$ source is redistributed as $0.25\ \mathrm{mA}$ into each of nodes 1 and 2.

> [!warning] Common Pitfall
> Do not compute $\mathbf{A}^{-1}$ in order to solve $\mathbf{A}\mathbf{x}=\mathbf{b}$. Forming the inverse costs about $2n^3$ flops versus $\tfrac23n^3$ for an $\mathbf{LU}$ solve, is less accurate in floating point, and for sparse matrices is catastrophic: the inverse of a sparse nodal matrix is typically completely full (every entry of $\mathbf{Z}$ above is nonzero although $\mathbf{G}$ could be sparse). Write $\mathbf{A}^{-1}\mathbf{b}$ in formulas; compute it as `A\b` (MATLAB) or `numpy.linalg.solve(A, b)`.

> [!tip] Engineering Note
> The nodal resistance matrix $\mathbf{Z}=\mathbf{G}^{-1}$ is the natural object in **sensitivity analysis**. By differentiating $\mathbf{G}\mathbf{v}=\mathbf{i}_s$, $\partial\mathbf{v}/\partial g=-\mathbf{Z}(\partial\mathbf{G}/\partial g)\mathbf{v}$; with the stamp $\partial\mathbf{G}/\partial g=\mathbf{u}\mathbf{u}^\mathsf{T}$ this is the derivative form of Sherman–Morrison. Circuit simulators compute such sensitivities with one extra *adjoint* solve using $\mathbf{G}^\mathsf{T}$, reusing the factorisation — the same trick that underlies back-propagation in neural networks.

## 8.5 Rank

### Motivation

A square matrix is either invertible or not; the determinant is a yes/no test. But engineers constantly meet matrices that are rectangular (more measurements than unknowns, more branches than nodes) or singular for structural reasons (a floating subcircuit, a redundant sensor, a dependent equation written twice). The **rank** is the single number that says how much independent information a matrix carries. It decides whether a linear system is solvable, how many free parameters its solution has, how many independent KCL or KVL equations a network supplies, and — in Chapter 9 — the dimensions of the column space and null space.

### Echelon form

> [!definition] Row echelon form, pivots
> A matrix is in **row echelon form (REF)** if
> 1. all zero rows are at the bottom, and
> 2. the first nonzero entry (the **pivot** or leading entry) of each nonzero row lies strictly to the right of the pivot of the row above.
>
> The columns containing pivots are **pivot columns**; the others are **free columns**.

A typical REF pattern, with $\blacksquare$ a nonzero pivot and $*$ anything:

$$
\begin{bmatrix}\blacksquare&*&*&*&*\\0&0&\blacksquare&*&*\\0&0&0&\blacksquare&*\\0&0&0&0&0\end{bmatrix}.
$$

Gaussian elimination (§8.7) reduces every matrix to an REF using row operations. The REF itself is not unique (it depends on the operations chosen), but the *number* of pivots and the *positions* of the pivot columns are.

> [!definition] Rank
> The **rank** of $\mathbf{A}$, written $\operatorname{rank}\mathbf{A}$, is the number of pivots in any row echelon form of $\mathbf{A}$. Equivalently, it is
> - the maximum number of linearly independent rows of $\mathbf{A}$ (the **row rank**);
> - the maximum number of linearly independent columns of $\mathbf{A}$ (the **column rank**);
> - the order of the largest square submatrix with nonzero determinant (the largest nonzero **minor**).

A set of vectors $\mathbf{a}_1,\dots,\mathbf{a}_k$ is **linearly independent** if $c_1\mathbf{a}_1+\dots+c_k\mathbf{a}_k=\mathbf{0}$ forces $c_1=\dots=c_k=0$; otherwise one of them is a linear combination of the others. Chapter 9 develops this notion fully; here we need only the computational test supplied by elimination.

*Why the characterisations agree (sketch).* Row operations replace rows by combinations of rows and are reversible, so they preserve the set of all combinations of rows and hence the row rank; in an REF the nonzero rows are visibly independent, so row rank equals the number of pivots. Row operations also preserve every linear relation among the columns (if $\mathbf{A}\mathbf{c}=\mathbf{0}$ then $\mathbf{E}\mathbf{A}\mathbf{c}=\mathbf{0}$, and conversely since $\mathbf{E}$ is invertible), so the pivot columns of $\mathbf{A}$ are independent and every free column is a combination of them: column rank equals the number of pivots too. Finally, a nonzero $r\times r$ minor exists in the REF (the pivot rows and columns form a triangular block with nonzero diagonal), and no larger one can, since any larger set of rows is dependent. $\square$

The equality of row rank and column rank is one of the genuinely surprising facts of linear algebra: a $3\times100$ matrix cannot have more than three independent columns, because it has only three rows.

### Properties

> [!theorem] Properties of rank
> For $\mathbf{A}\in\mathbb{C}^{m\times n}$:
> 1. $0\le\operatorname{rank}\mathbf{A}\le\min(m,n)$; $\operatorname{rank}\mathbf{A}=0$ only for $\mathbf{A}=\mathbf{0}$. If $\operatorname{rank}\mathbf{A}=\min(m,n)$, $\mathbf{A}$ has **full rank**.
> 2. $\operatorname{rank}\mathbf{A}^\mathsf{T}=\operatorname{rank}\mathbf{A}^\mathsf{H}=\operatorname{rank}\mathbf{A}$.
> 3. $\operatorname{rank}(\mathbf{A}\mathbf{B})\le\min(\operatorname{rank}\mathbf{A},\operatorname{rank}\mathbf{B})$.
> 4. If $\mathbf{P}$ and $\mathbf{Q}$ are invertible, $\operatorname{rank}(\mathbf{P}\mathbf{A}\mathbf{Q})=\operatorname{rank}\mathbf{A}$.
> 5. A nonzero outer product $\mathbf{u}\mathbf{v}^\mathsf{T}$ has rank 1; conversely every rank-1 matrix is an outer product.
> 6. $\operatorname{rank}(\mathbf{A}+\mathbf{B})\le\operatorname{rank}\mathbf{A}+\operatorname{rank}\mathbf{B}$.
> 7. **Rank–nullity:** the solutions of $\mathbf{A}\mathbf{x}=\mathbf{0}$ form a set described by exactly $n-\operatorname{rank}\mathbf{A}$ free parameters (the *nullity*; see Chapter 9).
> 8. A square $n\times n$ matrix is invertible if and only if its rank is $n$.

Property 3 follows from the column view: every column of $\mathbf{A}\mathbf{B}$ is a combination of the columns of $\mathbf{A}$, and (transposing) every row of $\mathbf{A}\mathbf{B}$ is a combination of the rows of $\mathbf{B}$. It explains why $\mathbf{G}=\mathbf{A}\mathbf{G}_b\mathbf{A}^\mathsf{T}$ can be no better than the incidence matrix $\mathbf{A}$ it is built from.

> [!example] Worked Example 8.15 — Rank of a rectangular matrix
> Find the rank of
> $$
> \mathbf{M}=\begin{bmatrix}1&2&1&3\\2&4&0&2\\3&6&1&5\end{bmatrix}
> $$
> and express the free columns in terms of the pivot columns.
>
> **Step 1.** $R_2\leftarrow R_2-2R_1$, $R_3\leftarrow R_3-3R_1$:
> $$
> \begin{bmatrix}1&2&1&3\\0&0&-2&-4\\0&0&-2&-4\end{bmatrix}.
> $$
> **Step 2.** $R_3\leftarrow R_3-R_2$ gives a zero third row. The REF has pivots in columns 1 and 3, so $\operatorname{rank}\mathbf{M}=2$.
>
> **Step 3.** Row 3 of $\mathbf{M}$ is row 1 plus row 2 — the dependency the elimination discovered. Columns: $\mathbf{m}_2=2\mathbf{m}_1$ by inspection; for $\mathbf{m}_4=\alpha\mathbf{m}_1+\beta\mathbf{m}_3$, the second component gives $2\alpha=2$, so $\alpha=1$, and the first gives $1+\beta=3$, so $\beta=2$. Check the third: $3+2(1)=5$ ✓.
>
> **Reading.** Although $\mathbf{M}$ has four columns in $\mathbb{R}^3$, they span only a plane. A system $\mathbf{M}\mathbf{x}=\mathbf{b}$ can be solved only for $\mathbf{b}$ in that plane (those with $b_3=b_1+b_2$), and then with $4-2=2$ free parameters.

> [!example] Worked Example 8.16 — Rank depending on a parameter
> For which $k$ does $\mathbf{B}=\begin{bmatrix}1&1&1\\1&2&3\\1&3&k\end{bmatrix}$ have rank 3? What is the rank otherwise?
>
> **Step 1.** $R_2\leftarrow R_2-R_1$, $R_3\leftarrow R_3-R_1$: rows $(0,1,2)$ and $(0,2,k-1)$.
>
> **Step 2.** $R_3\leftarrow R_3-2R_2$: row $(0,0,k-5)$.
>
> **Step 3.** The pivots are $1,1,k-5$, so $\det\mathbf{B}=k-5$. For $k\ne5$, $\operatorname{rank}\mathbf{B}=3$; for $k=5$ the third row vanishes and $\operatorname{rank}\mathbf{B}=2$ (the $2\times2$ minor $\begin{vmatrix}1&1\\1&2\end{vmatrix}=1\ne0$ confirms it is not lower). At $k=5$, row 3 $=2\,$row 2 $-$ row 1. We use this matrix again in §8.6.

### Rank in network topology

For a connected network with $N$ nodes (including ground) and $b$ branches, the full $N\times b$ incidence matrix has every column containing one $+1$ and one $-1$, so its rows sum to zero: the rows are dependent and the rank is at most $N-1$. In fact it is exactly $N-1$ for a connected graph — which is why we delete the ground row and keep the $(N-1)\times b$ **reduced** incidence matrix, of full row rank. The number of independent KCL equations is therefore $N-1$, and the number of independent KVL (loop/mesh) equations is $b-(N-1)$, the nullity of the reduced incidence matrix. For the running example $N=4$, $b=6$: three nodal equations, three independent loops.

> [!example] Worked Example 8.17 — A floating network has a singular nodal matrix
> Remove the three branches to ground ($a,d,e$) from the running example, leaving branches $b$ (2 mS, 1–2), $c$ (1 mS, 2–3), $f$ (1 mS, 1–3), with no connection to ground. Write the nodal matrix and find its rank.
>
> **Step 1.** Stamping the three branches:
> $$
> \mathbf{L}=\begin{bmatrix}3&-2&-1\\-2&3&-1\\-1&-1&2\end{bmatrix}.
> $$
> **Step 2.** Every row sums to zero: $\mathbf{L}(1,1,1)^\mathsf{T}=\mathbf{0}$. A nonzero solution of $\mathbf{L}\mathbf{x}=\mathbf{0}$ exists, so $\mathbf{L}$ is singular and $\operatorname{rank}\mathbf{L}\le2$.
>
> **Step 3.** Elimination: $R_2\leftarrow R_2+\tfrac23R_1$ gives $(0,\tfrac53,-\tfrac53)$; $R_3\leftarrow R_3+\tfrac13R_1$ gives $(0,-\tfrac53,\tfrac53)$; then $R_3\leftarrow R_3+R_2$ gives $\mathbf{0}$. Two pivots: $\operatorname{rank}\mathbf{L}=2$.
>
> **Reading.** Physically, adding the same voltage to every node of a floating network changes no branch voltage and therefore no current; the solution can only ever be determined up to that common offset. Choosing any node as the reference (deleting its row and column) restores full rank. Matrices with zero row sums and non-positive off-diagonals like $\mathbf{L}$ are **graph Laplacians**; their rank is $N-c$ where $c$ is the number of connected pieces.

### Numerical rank

In exact arithmetic rank is an integer and elimination finds it. In floating point, or with measured data, "zero" pivots come out as $10^{-14}$ or as noise, and the question becomes *how close is the matrix to one of lower rank?* The matrix $\begin{bmatrix}1&1\\1&1.0001\end{bmatrix}$ has rank 2, but a perturbation of size $10^{-4}$ in one entry makes it rank 1; for most engineering purposes it behaves like a rank-1 matrix. The robust tool is the singular value decomposition (Chapter 10): the singular values of this matrix are about $2.00005$ and $5.0\times10^{-5}$, and the **numerical rank** is the number of singular values above a tolerance tied to the data accuracy. The ratio of the largest to the smallest singular value is the 2-norm condition number of §8.10.

> [!warning] Common Pitfall
> Rank is not the number of nonzero rows of the *original* matrix, nor the number of nonzero entries on its diagonal. $\begin{bmatrix}1&2\\2&4\end{bmatrix}$ has no zero rows and a nonzero diagonal but rank 1. Always reduce to echelon form (or use singular values) before counting.

> [!tip] Engineering Note
> Rank deficiency in a circuit simulator usually means a modelling error, and SPICE reports it as a "singular matrix" error. The usual culprits: a node with no DC path to ground (e.g. the junction of two series capacitors), a loop consisting only of ideal voltage sources and inductors (KVL makes one of the source equations redundant or contradictory), or a cut-set of only current sources and capacitors. Each is a rank defect of a structural submatrix and can be detected from the network topology before any number is computed.

## 8.6 Systems of linear equations

### Formulation

> [!definition] Linear system
> A system of $m$ linear equations in $n$ unknowns,
> $$
> \begin{aligned}
> a_{11}x_1+a_{12}x_2+\dots+a_{1n}x_n&=b_1,\\
> &\ \ \vdots\\
> a_{m1}x_1+a_{m2}x_2+\dots+a_{mn}x_n&=b_m,
> \end{aligned}
> $$
> is written $\mathbf{A}\mathbf{x}=\mathbf{b}$ with **coefficient matrix** $\mathbf{A}\in\mathbb{C}^{m\times n}$ and **right-hand side** $\mathbf{b}\in\mathbb{C}^m$. The **augmented matrix** is $[\mathbf{A}\mid\mathbf{b}]$. The system is **homogeneous** if $\mathbf{b}=\mathbf{0}$. It is **consistent** if it has at least one solution, **inconsistent** otherwise.

Systems are classified by shape: **square** ($m=n$, the usual case for well-posed circuit equations), **overdetermined** ($m>n$, typical of measurements; usually inconsistent and solved in the least-squares sense of Chapter 9) and **underdetermined** ($m<n$; consistent systems have infinitely many solutions, and a particular one is chosen by an extra criterion such as minimum norm).

### Existence and uniqueness

> [!theorem] Rouché–Capelli (Kronecker–Capelli) theorem
> Let $\mathbf{A}$ be $m\times n$ and $r=\operatorname{rank}\mathbf{A}$.
> 1. $\mathbf{A}\mathbf{x}=\mathbf{b}$ is consistent if and only if $\operatorname{rank}[\mathbf{A}\mid\mathbf{b}]=\operatorname{rank}\mathbf{A}$.
> 2. If it is consistent and $r=n$, the solution is unique.
> 3. If it is consistent and $r<n$, there are infinitely many solutions, described by $n-r$ free parameters.

*Proof sketch.* Reduce $[\mathbf{A}\mid\mathbf{b}]$ to echelon form. If the augmented column contains a pivot, some row reads $0=c$ with $c\ne0$: inconsistent, and the augmented rank exceeds $r$. Otherwise the $r$ nonzero rows can be solved by back substitution for the $r$ pivot variables in terms of the $n-r$ free variables, which may be chosen arbitrarily. $\square$

In terms of the column view: $\mathbf{A}\mathbf{x}=\mathbf{b}$ is solvable exactly when $\mathbf{b}$ is a combination of the columns of $\mathbf{A}$, and appending such a $\mathbf{b}$ cannot raise the rank.

> [!theorem] Structure of the solution set
> If $\mathbf{x}_p$ is any one solution of $\mathbf{A}\mathbf{x}=\mathbf{b}$, every solution has the form
> $$
> \mathbf{x}=\mathbf{x}_p+\mathbf{x}_h,
> $$
> where $\mathbf{x}_h$ ranges over all solutions of the homogeneous system $\mathbf{A}\mathbf{x}_h=\mathbf{0}$.

*Proof.* If $\mathbf{A}\mathbf{x}=\mathbf{b}$ then $\mathbf{A}(\mathbf{x}-\mathbf{x}_p)=\mathbf{b}-\mathbf{b}=\mathbf{0}$; conversely $\mathbf{A}(\mathbf{x}_p+\mathbf{x}_h)=\mathbf{b}+\mathbf{0}$. $\square$

This "particular plus homogeneous" structure reappears verbatim for linear differential equations (Chapters 11–12), where it is the split into forced and natural response. Its companion is **superposition**: if $\mathbf{A}\mathbf{x}_1=\mathbf{b}_1$ and $\mathbf{A}\mathbf{x}_2=\mathbf{b}_2$ then $\mathbf{A}(\alpha\mathbf{x}_1+\beta\mathbf{x}_2)=\alpha\mathbf{b}_1+\beta\mathbf{b}_2$. A circuit's response to several sources is the sum of its responses to each — which is just the linearity of $\mathbf{x}=\mathbf{A}^{-1}\mathbf{b}$ in $\mathbf{b}$.

### Homogeneous systems

A homogeneous system always has the **trivial solution** $\mathbf{x}=\mathbf{0}$. By Rouché–Capelli it has nontrivial solutions if and only if $\operatorname{rank}\mathbf{A}<n$. Two consequences:

- A homogeneous system with **more unknowns than equations** ($n>m$) always has nontrivial solutions, since $r\le m<n$.
- A **square** homogeneous system has nontrivial solutions if and only if $\det\mathbf{A}=0$. This is how eigenvalue problems $(\mathbf{A}-\lambda\mathbf{I})\mathbf{x}=\mathbf{0}$ lead to the characteristic equation $\det(\mathbf{A}-\lambda\mathbf{I})=0$ (Chapter 10), and how a circuit's natural frequencies arise from $\det\mathbf{Y}(s)=0$: a nonzero response with no sources needs a singular system matrix.

### Geometry

With two unknowns, each equation $a_{i1}x+a_{i2}y=b_i$ is a line in the plane. Two lines either cross at one point (unique solution, $\det\mathbf{A}\ne0$, the case drawn in §8.3), are parallel and distinct (no solution), or coincide (infinitely many). With three unknowns each equation is a plane in space; three planes in general position meet in a single point, but they may also share a common line (rank 2, consistent), be arranged like the faces of a triangular prism with no common point (rank 2, inconsistent), contain two parallel planes, or all coincide (rank 1).

```desmos-graph
left=-2; right=6; bottom=-1; top=4
---
x+2y=4|blue
x+2y=1|red
2x+4y=8|green|dashed
```

*Singular $2\times2$ systems. The blue line $x+2y=4$ and red line $x+2y=1$ are parallel: the system has no solution ($\operatorname{rank}\mathbf{A}=1$, $\operatorname{rank}[\mathbf{A}\mid\mathbf{b}]=2$). The dashed green line $2x+4y=8$ lies exactly on the blue one: together they give infinitely many solutions (both ranks 1, one free parameter).*

> [!example] Worked Example 8.18 — The three cases for two unknowns
> Classify (a) $x+2y=4,\ 3x-y=5$; (b) $x+2y=4,\ x+2y=1$; (c) $x+2y=4,\ 2x+4y=8$.
>
> **(a)** $\det\begin{bmatrix}1&2\\3&-1\end{bmatrix}=-7\ne0$: unique solution. Cramer: $x=\dfrac{-4-10}{-7}=2$, $y=\dfrac{5-12}{-7}=1$.
>
> **(b)** $R_2\leftarrow R_2-R_1$ gives $[\,0\ \ 0\mid-3\,]$, i.e. $0=-3$. $\operatorname{rank}\mathbf{A}=1<\operatorname{rank}[\mathbf{A}\mid\mathbf{b}]=2$: inconsistent.
>
> **(c)** $R_2\leftarrow R_2-2R_1$ gives $[\,0\ \ 0\mid0\,]$. Both ranks are 1, $n-r=1$ free parameter. Let $y=t$; then $x=4-2t$, so $(x,y)=(4,0)+t(-2,1)$ — a particular solution plus any multiple of the homogeneous solution $(-2,1)$, which satisfies $x+2y=0$.

> [!example] Worked Example 8.19 — A system with parameters
> Discuss, in terms of $k$ and $m$, the solutions of
> $$
> \begin{aligned}
> x+y+z&=3,\\
> x+2y+3z&=6,\\
> x+3y+kz&=m.
> \end{aligned}
> $$
> **Step 1.** The coefficient matrix is $\mathbf{B}$ of Worked Example 8.16. Apply the same row operations to the augmented matrix:
> $$
> \left[\begin{array}{ccc|c}1&1&1&3\\1&2&3&6\\1&3&k&m\end{array}\right]\to
> \left[\begin{array}{ccc|c}1&1&1&3\\0&1&2&3\\0&2&k-1&m-3\end{array}\right]\to
> \left[\begin{array}{ccc|c}1&1&1&3\\0&1&2&3\\0&0&k-5&m-9\end{array}\right].
> $$
> **Step 2 (classify).**
> - $k\ne5$: three pivots, rank 3 $=n$: unique solution $z=\dfrac{m-9}{k-5}$, $y=3-2z$, $x=3-y-z=z$.
> - $k=5$, $m\ne9$: last row reads $0=m-9\ne0$; $\operatorname{rank}\mathbf{A}=2<\operatorname{rank}[\mathbf{A}\mid\mathbf{b}]=3$: no solution.
> - $k=5$, $m=9$: both ranks 2, one free parameter. With $z=t$: $y=3-2t$, $x=3-(3-2t)-t=t$, so $(x,y,z)=(0,3,0)+t(1,-2,1)$.
>
> **Step 3 (check).** For $k=6$, $m=10$: $z=1$, $y=1$, $x=1$; indeed $1+1+1=3$, $1+2+3=6$, $1+3+6=10$ ✓. For the infinite family, $(t,3-2t,t)$ in the third equation with $k=5$: $t+9-6t+5t=9$ ✓.
>
> **Geometry.** For $k=5$ the three planes share a common line when $m=9$ and form a "prism" (pairwise intersecting lines, no common point) when $m\ne9$.

> [!example] Worked Example 8.20 — Solving the floating network: offsets and superposition
> For the floating network $\mathbf{L}$ of Worked Example 8.17, inject $1\ \mathrm{mA}$ into node 1 and extract it from node 3: $\mathbf{L}\mathbf{v}=(1,0,-1)^\mathsf{T}$. Find all solutions.
>
> **Step 1 (consistency).** Since every column of $\mathbf{L}$ sums to zero, so does every vector $\mathbf{L}\mathbf{v}$: the system is solvable only if the injected currents sum to zero. Here $1+0-1=0$ — physically, current into a floating network must also leave it. Reduction of $[\mathbf{L}\mid\mathbf{b}]$ confirms both ranks are 2.
>
> **Step 2 (particular solution).** Fix the free variable by grounding node 3, $v_3=0$; the first two equations become $3v_1-2v_2=1$, $-2v_1+3v_2=0$, with solution $v_1=\tfrac35$, $v_2=\tfrac25\ \mathrm{V}$. Check row 3: $-\tfrac35-\tfrac25+0=-1$ ✓.
>
> **Step 3 (general solution).** The homogeneous solutions are $t(1,1,1)^\mathsf{T}$, so
> $$
> \mathbf{v}=\begin{bmatrix}0.6\\0.4\\0\end{bmatrix}+t\begin{bmatrix}1\\1\\1\end{bmatrix},\qquad t\in\mathbb{R}.
> $$
> Every physically meaningful quantity — the branch voltages $v_1-v_2=0.2$, $v_2-v_3=0.4$, $v_1-v_3=0.6\ \mathrm{V}$ — is independent of $t$. The resistance between nodes 1 and 3 is $0.6\ \mathrm{V}/1\ \mathrm{mA}=0.6\ \mathrm{k\Omega}$, which matches the series–parallel calculation: $(\tfrac12+1)\ \mathrm{k\Omega}$ in parallel with $1\ \mathrm{k\Omega}$ is $\tfrac{1.5}{2.5}=0.6\ \mathrm{k\Omega}$ ✓.
>
> **Step 4 (superposition).** Adding a second injection pattern $(0,1,-1)^\mathsf{T}$ (with solution $(0.4,0.6,0)^\mathsf{T}+t(1,1,1)^\mathsf{T}$ by the same method — check: $3(0.4)-2(0.6)=0$, $-2(0.4)+3(0.6)=1$ ✓) gives the response to $(1,1,-2)^\mathsf{T}$ as the sum $(1,1,0)^\mathsf{T}+t(1,1,1)^\mathsf{T}$.

> [!warning] Common Pitfall
> Counting equations and unknowns does not settle solvability. Three equations in three unknowns may have no solution or infinitely many (Worked Example 8.19 with $k=5$); five equations in three unknowns may have exactly one (if two are redundant). Only the ranks of $\mathbf{A}$ and $[\mathbf{A}\mid\mathbf{b}]$ decide.

> [!tip] Engineering Note
> Well-posed circuit formulations (nodal, mesh, modified nodal) are constructed to give square nonsingular systems when the circuit itself is well-posed. Overdetermined systems arise on purpose in **state estimation**: a power-system control centre measures more line flows and voltages than there are unknown bus voltages, precisely so that bad measurements can be detected from the inconsistency, and solves the system in the weighted least-squares sense. Underdetermined systems arise in **compressed sensing** and in redundant actuator allocation, where the extra freedom is spent on minimising energy or sparsity.

## 8.7 Gaussian elimination

### Motivation

Gaussian elimination is the workhorse algorithm of computational linear algebra. Every dense linear solver — MATLAB's backslash, LAPACK's `dgesv`, the DC operating-point and transient solvers inside SPICE — is a carefully engineered version of it. The idea is the one taught in school for two equations: use one equation to eliminate a variable from the others, repeat until one equation in one unknown remains, then substitute back. Doing this systematically turns $\mathbf{A}\mathbf{x}=\mathbf{b}$ into an equivalent **upper triangular** system $\mathbf{U}\mathbf{x}=\mathbf{c}$, which is trivial to solve.

### Elementary row operations and equivalence

> [!definition] Elementary row operations
> 1. **Replacement:** add a multiple of one row to another, $R_i\leftarrow R_i+cR_k$.
> 2. **Interchange:** swap two rows, $R_i\leftrightarrow R_k$.
> 3. **Scaling:** multiply a row by a nonzero constant, $R_i\leftarrow cR_i$.
>
> Applied to the augmented matrix $[\mathbf{A}\mid\mathbf{b}]$, they produce an **equivalent** system: one with exactly the same solution set.

Equivalence holds because each operation is a left-multiplication by an invertible elementary matrix $\mathbf{E}$ (§8.4): $\mathbf{A}\mathbf{x}=\mathbf{b}\iff\mathbf{E}\mathbf{A}\mathbf{x}=\mathbf{E}\mathbf{b}$. In the plane, a replacement swaps one line for a different line through the *same* intersection point — the solution never moves.

```desmos-graph
left=-1; right=4; bottom=-1; top=3
---
5x-3y=7|blue
-3x+7y=1|red|dashed
y=1|green
(2,1)|label:(2, 1)
```

*Elimination in pictures for the mesh system of §8.3. The replacement $R_2\leftarrow R_2+\tfrac35R_1$ turns the red dashed line $-3x+7y=1$ into the green line $\tfrac{26}{5}y=\tfrac{26}{5}$, i.e. $y=1$, which passes through the same point $(2,1)$ but no longer involves $x$. Back substitution then reads off $x$ from the blue line.*

### The algorithm

**Forward elimination.** For $k=1,2,\dots,n-1$ (the current pivot column):

1. The **pivot** is $a_{kk}^{(k)}$, the current $(k,k)$ entry. (If it is zero, or small, interchange rows first — see pivoting below.)
2. For each row $i=k+1,\dots,n$ compute the **multiplier** $m_{ik}=a_{ik}^{(k)}/a_{kk}^{(k)}$.
3. Replace $R_i\leftarrow R_i-m_{ik}R_k$, which makes $a_{ik}^{(k+1)}=0$:
$$
a_{ij}^{(k+1)}=a_{ij}^{(k)}-m_{ik}a_{kj}^{(k)},\quad j=k+1,\dots,n,\qquad b_i^{(k+1)}=b_i^{(k)}-m_{ik}b_k^{(k)}.
$$

After $n-1$ stages the matrix is upper triangular, $\mathbf{U}=[u_{ij}]$ with $u_{ij}=a_{ij}^{(i)}$, and the right-hand side has become $\mathbf{c}$.

**Back substitution.** Solve $\mathbf{U}\mathbf{x}=\mathbf{c}$ from the bottom up:

$$
x_n=\frac{c_n}{u_{nn}},\qquad x_i=\frac{1}{u_{ii}}\Big(c_i-\sum_{j=i+1}^{n}u_{ij}x_j\Big),\quad i=n-1,\dots,1.
$$

The process succeeds — in exact arithmetic, with row interchanges allowed — exactly when $\mathbf{A}$ is nonsingular, and the pivots $u_{11},\dots,u_{nn}$ multiply to $\pm\det\mathbf{A}$.

### Operation count

At stage $k$ there are $n-k$ rows to update, each with $n-k$ matrix entries (plus one right-hand-side entry), each costing one multiply and one subtract. Ignoring lower-order terms,

$$
\text{forward elimination}\approx\sum_{k=1}^{n-1}2(n-k)^2=2\sum_{p=1}^{n-1}p^2\approx\frac{2n^3}{3}\ \text{flops},
$$

while back substitution costs $\sum_i 2(n-i)\approx n^2$ flops, as does applying the elimination to one extra right-hand side. Doubling $n$ multiplies the elimination time by eight. For $n=1000$, $\tfrac23n^3\approx6.7\times10^8$ flops: well under a second. For $n=10^5$ dense it is $6.7\times10^{14}$ — minutes to hours, and $80\ \mathrm{GB}$ just to store the matrix in double precision. Large problems are only tractable because they are sparse (§8.11).

> [!example] Worked Example 8.21 — The running example by Gaussian elimination
> Solve $\mathbf{G}\mathbf{v}=\mathbf{i}_s$ with the augmented matrix
> $$
> \left[\begin{array}{ccc|c}4&-2&-1&8\\-2&4&-1&2\\-1&-1&4&1\end{array}\right].
> $$
> **Stage 1** (pivot $4$). Multipliers $m_{21}=-2/4=-\tfrac12$, $m_{31}=-\tfrac14$.
> $$
> \begin{aligned}
> R_2&\leftarrow R_2+\tfrac12R_1: & (0,\ 4-1,\ -1-\tfrac12\mid 2+4)&=(0,\ 3,\ -\tfrac32\mid 6),\\
> R_3&\leftarrow R_3+\tfrac14R_1: & (0,\ -1-\tfrac12,\ 4-\tfrac14\mid 1+2)&=(0,\ -\tfrac32,\ \tfrac{15}4\mid 3).
> \end{aligned}
> $$
> **Stage 2** (pivot $3$). Multiplier $m_{32}=(-\tfrac32)/3=-\tfrac12$.
> $$
> R_3\leftarrow R_3+\tfrac12R_2:\quad(0,\ 0,\ \tfrac{15}4-\tfrac34\mid 3+3)=(0,\ 0,\ 3\mid 6).
> $$
> **Triangular system.**
> $$
> \left[\begin{array}{ccc|c}4&-2&-1&8\\0&3&-\tfrac32&6\\0&0&3&6\end{array}\right].
> $$
> **Back substitution.** $v_3=6/3=2$; $v_2=(6+\tfrac32\cdot2)/3=9/3=3$; $v_1=(8+2\cdot3+1\cdot2)/4=16/4=4$. So $\mathbf{v}=(4,3,2)^\mathsf{T}\ \mathrm{V}$. ✓
>
> **By-products.** $\det\mathbf{G}=4\cdot3\cdot3=36$. The multipliers $-\tfrac12,-\tfrac14,-\tfrac12$ will reappear as the entries of $\mathbf{L}$ in §8.9.
>
> **Circuit reading.** Stage 1 is Kron reduction of node 1: the $2\times2$ block $\begin{bmatrix}3&-\tfrac32\\-\tfrac32&\tfrac{15}4\end{bmatrix}$ is the nodal matrix, seen from nodes 2 and 3, of the network with node 1 eliminated, and $(6,3)$ are the correspondingly redistributed source currents. Gaussian elimination *is* repeated node elimination.

### Pivoting

If a pivot is zero the algorithm cannot proceed as written, even though the matrix may be perfectly nonsingular; if a pivot is merely small, the multipliers are large and rounding errors are amplified.

> [!definition] Partial pivoting
> At stage $k$, choose the row $p\ge k$ for which $\lvert a_{pk}^{(k)}\rvert$ is largest, and interchange rows $p$ and $k$ before eliminating. Then every multiplier satisfies $\lvert m_{ik}\rvert\le1$.

**Complete pivoting** searches the whole remaining submatrix and swaps columns as well as rows; it is more stable in theory but rarely needed. The quantity that controls stability is the **growth factor** $\rho=\max_{i,j,k}\lvert a_{ij}^{(k)}\rvert/\max_{i,j}\lvert a_{ij}\rvert$: backward error analysis (Wilkinson) shows that the computed solution is the exact solution of a nearby system $(\mathbf{A}+\delta\mathbf{A})\hat{\mathbf{x}}=\mathbf{b}$ with $\lVert\delta\mathbf{A}\rVert\lesssim c\,n\,\rho\,u\lVert\mathbf{A}\rVert$, where $u\approx1.1\times10^{-16}$ is the double-precision unit roundoff. With partial pivoting $\rho\le2^{n-1}$ in the worst case, but in practice $\rho$ is almost always small, which is why partial pivoting is the universal default.

> [!theorem] When pivoting is unnecessary
> If $\mathbf{A}$ is **strictly diagonally dominant by columns** ($\lvert a_{jj}\rvert>\sum_{i\ne j}\lvert a_{ij}\rvert$ for every $j$), or **symmetric positive definite**, then Gaussian elimination without interchanges never meets a zero pivot, the growth factor is bounded ($\rho\le2$ and $\rho\le1$ respectively), and the matrix is nonsingular. Diagonal dominance by rows also guarantees nonsingularity and nonzero pivots.

*Sketch for diagonal dominance.* After one elimination stage, the reduced $(n-1)\times(n-1)$ block is again diagonally dominant (a short calculation with the triangle inequality), so induction shows every pivot is nonzero. $\square$

Nodal matrices of resistive networks with at least one branch to ground at every node are strictly diagonally dominant (diagonal = sum of all conductances at the node, off-diagonal magnitudes = the conductances to other nodes) and symmetric positive definite, so they can be eliminated in any order without pivoting for stability. Simulators exploit this freedom to choose an elimination order that preserves sparsity (§8.11).

> [!example] Worked Example 8.22 — A zero pivot from a voltage source
> Replace the $8\ \mathrm{mA}$ source at node 1 of the running example by an ideal $4\ \mathrm{V}$ voltage source from ground to node 1, and let $j$ be the current (mA) flowing from node 1 into the source. With unknowns ordered $(j,v_1,v_2,v_3)$ and the source equation written first, modified nodal analysis gives
> $$
> \left[\begin{array}{cccc|c}0&1&0&0&4\\1&4&-2&-1&0\\0&-2&4&-1&2\\0&-1&-1&4&1\end{array}\right].
> $$
> **Step 1.** The $(1,1)$ pivot is $0$. Interchange $R_1\leftrightarrow R_2$:
> $$
> \left[\begin{array}{cccc|c}1&4&-2&-1&0\\0&1&0&0&4\\0&-2&4&-1&2\\0&-1&-1&4&1\end{array}\right].
> $$
> Column 1 is already clear below the pivot.
>
> **Step 2** (pivot $1$ in column 2; strict partial pivoting would pick the $-2$ in row 3, which gives the same exact answer). $R_3\leftarrow R_3+2R_2$, $R_4\leftarrow R_4+R_2$:
> rows 3, 4 become $(0,0,4,-1\mid10)$ and $(0,0,-1,4\mid5)$.
>
> **Step 3.** $R_4\leftarrow R_4+\tfrac14R_3$: $(0,0,0,\tfrac{15}4\mid 7.5)$.
>
> **Step 4 (back substitution).** $v_3=7.5/3.75=2$; row 3 reads $4v_2-v_3=10$, so $v_2=3$; row 2 gives $v_1=4$; and $j+4(4)-2(3)-1(2)=0\Rightarrow j=-8\ \mathrm{mA}$.
>
> **Reading.** The node voltages are unchanged, as they must be: the voltage source holds node 1 at exactly the voltage the current source produced. The source current $j=-8\ \mathrm{mA}$ means $8\ \mathrm{mA}$ flows *out* of the source into node 1 — the same $8\ \mathrm{mA}$ the current source injected. The determinant is $(-1)^1\cdot1\cdot1\cdot4\cdot\tfrac{15}4=-15$, agreeing with Worked Example 8.8. MNA matrices are symmetric but **indefinite** (zero diagonal blocks), so they need pivoting where pure nodal matrices do not.

> [!example] Worked Example 8.23 — Why small pivots are dangerous
> Solve
> $$
> \begin{aligned}
> 0.0001\,x_1+x_2&=1,\\
> x_1+x_2&=2
> \end{aligned}
> $$
> in **three-significant-figure** arithmetic (each result rounded to 3 significant figures), first without and then with pivoting. The exact solution is $x_1=1/0.9999=1.00010\ldots$, $x_2=0.9998/0.9999=0.99990\ldots$.
>
> **Without pivoting.** $m_{21}=1/0.0001=10000$. Then $a_{22}=1-10000=-9999\to-1.00\times10^4$ and $b_2=2-10000=-9998\to-1.00\times10^4$. So $x_2=1.00$, and $x_1=(1-x_2)/0.0001=0/0.0001=0$. The answer $x_1=0$ is completely wrong: the information in the second equation ("$x_1+x_2=2$") was swamped by the huge multiple of the first.
>
> **With partial pivoting.** $\lvert1\rvert>\lvert0.0001\rvert$, so swap the equations. Now $m_{21}=0.0001$, $a_{22}=1-0.0001=0.9999\to1.00$, $b_2=1-0.0002=0.9998\to1.00$. So $x_2=1.00$ and $x_1=2-1.00=1.00$. Both components are correct to the working precision.
>
> **Moral.** The matrix is well conditioned (condition number about 4, §8.10); the failure without pivoting is entirely the algorithm's fault. In double precision the same thing happens with a pivot of $10^{-20}$.

> [!example] Worked Example 8.24 — Elimination with complex (phasor) entries
> A two-node AC network has a $1\ \mathrm{S}$ conductance from each node to ground and a capacitive susceptance $j1\ \mathrm{S}$ between the nodes. A phasor current $1\angle0^\circ\ \mathrm{A}$ is injected at node 1. Find $V_1$ and $V_2$.
>
> **Step 1.** $\mathbf{Y}=\begin{bmatrix}1+j&-j\\-j&1+j\end{bmatrix}$, $\mathbf{I}=(1,0)^\mathsf{T}$.
>
> **Step 2.** $m_{21}=\dfrac{-j}{1+j}=\dfrac{-j(1-j)}{2}=\dfrac{-1-j}{2}$. Then
> $$
> u_{22}=(1+j)-m_{21}(-j)=(1+j)+\frac{(-1-j)j}{2}=(1+j)+\frac{1-j}{2}=\frac{3+j}{2},\qquad c_2=0-m_{21}\cdot1=\frac{1+j}{2}.
> $$
> **Step 3.**
> $$
> V_2=\frac{(1+j)/2}{(3+j)/2}=\frac{(1+j)(3-j)}{10}=\frac{4+2j}{10}=0.4+j0.2\ \mathrm{V},
> $$
> $$
> V_1=\frac{1+jV_2}{1+j}=\frac{1+j0.4-0.2}{1+j}=\frac{(0.8+j0.4)(1-j)}{2}=0.6-j0.2\ \mathrm{V}.
> $$
> **Step 4 (check).** $(1+j)(0.6-j0.2)-j(0.4+j0.2)=(0.8+j0.4)+(0.2-j0.4)=1$ ✓; $-j(0.6-j0.2)+(1+j)(0.4+j0.2)=(-0.2-j0.6)+(0.2+j0.6)=0$ ✓.
>
> **Reading.** Complex arithmetic changes nothing in the algorithm; each complex flop costs about four real ones. $\lvert V_2\rvert\approx0.447\ \mathrm{V}$ with phase $+26.6^\circ$: the capacitor couples a leading current into node 2.

> [!warning] Common Pitfall
> When doing elimination by hand, apply every row operation to the right-hand side as well (or carry the augmented matrix), and do **not** scale rows by convenience factors in the middle of a determinant calculation without tracking the factor. Scaling is harmless for solving but changes the determinant.

> [!tip] Engineering Note
> Elimination order has a physical meaning in circuits: eliminating a node is Kron reduction, and the "fill-in" it creates is the set of new branches the star–mesh transformation adds between the node's neighbours. A node with many neighbours (a ground plane, a power rail) creates a dense clique when eliminated, which is why simulators eliminate such nodes last.

## 8.8 Gauss–Jordan elimination

### Motivation

Gaussian elimination stops at a triangular matrix and finishes with back substitution. **Gauss–Jordan elimination** keeps going: it scales every pivot to 1 and clears the entries *above* each pivot as well as below, until the coefficient part is as simple as possible. For a nonsingular square system that simplest form is the identity, and the solution can be read off directly. The method costs more than Gaussian elimination for a single solve, but it has three uses that make it indispensable: it computes inverses (by reducing $[\mathbf{A}\mid\mathbf{I}]$), it produces the unique **reduced row echelon form** that settles questions of rank and dependence unambiguously, and it writes the complete solution of an underdetermined system in terms of free parameters.

### Reduced row echelon form

> [!definition] Reduced row echelon form (RREF)
> A matrix is in **reduced row echelon form** if it is in row echelon form and, in addition,
> 1. every pivot equals 1, and
> 2. every pivot is the only nonzero entry in its column.

For example,

$$
\begin{bmatrix}1&2&0&1&0\\0&0&1&2&0\\0&0&0&0&1\\0&0&0&0&0\end{bmatrix}
$$

is in RREF, with pivot columns 1, 3, 5.

> [!theorem] Uniqueness of the RREF
> Every matrix is row-equivalent to exactly one matrix in reduced row echelon form, written $\operatorname{rref}(\mathbf{A})$.

*Proof sketch.* Row operations preserve all linear relations among columns (§8.5). In an RREF those relations are visible: each pivot column is a standard basis vector $\mathbf{e}_1,\mathbf{e}_2,\dots$ (in order) and each free column lists the coefficients expressing it in terms of the earlier pivot columns. Since the relations are determined by $\mathbf{A}$, so is the RREF. $\square$

Uniqueness is what makes the RREF a reliable bookkeeping device: two people reducing the same matrix in different orders reach the same answer (whereas their REFs may differ).

### The algorithm

For each column from left to right:

1. Find a row at or below the current row with a nonzero entry in this column (with partial pivoting, the one of largest magnitude); swap it into the current row. If there is none, the column is free: move to the next column, staying on the same row.
2. **Normalise:** divide the current row by its pivot so the pivot becomes 1.
3. **Eliminate:** subtract multiples of the current row from *every other row*, above and below, to make the rest of the column zero.
4. Move to the next row and column.

Applied to $[\mathbf{A}\mid\mathbf{b}]$ with $\mathbf{A}$ nonsingular, the result is $[\mathbf{I}\mid\mathbf{x}]$.

**Cost.** Clearing above as well as below the pivots roughly doubles the work in the elimination phase: about $n^3$ flops for one right-hand side, versus $\tfrac23n^3+n^2$ for Gaussian elimination with back substitution. For solving a single system Gauss–Jordan is therefore 50% more expensive and is not used in production solvers; its value is conceptual and for the special tasks below.

> [!example] Worked Example 8.25 — The running example by Gauss–Jordan
> Reduce $\left[\begin{array}{ccc|c}4&-2&-1&8\\-2&4&-1&2\\-1&-1&4&1\end{array}\right]$ to RREF.
>
> **Column 1.** $R_1\leftarrow\tfrac14R_1$: $(1,-\tfrac12,-\tfrac14\mid2)$. Then $R_2\leftarrow R_2+2R_1$: $(0,3,-\tfrac32\mid6)$; $R_3\leftarrow R_3+R_1$: $(0,-\tfrac32,\tfrac{15}4\mid3)$.
>
> **Column 2.** $R_2\leftarrow\tfrac13R_2$: $(0,1,-\tfrac12\mid2)$. Then $R_1\leftarrow R_1+\tfrac12R_2$: $(1,0,-\tfrac12\mid3)$; $R_3\leftarrow R_3+\tfrac32R_2$: $(0,0,3\mid6)$.
>
> **Column 3.** $R_3\leftarrow\tfrac13R_3$: $(0,0,1\mid2)$. Then $R_1\leftarrow R_1+\tfrac12R_3$: $(1,0,0\mid4)$; $R_2\leftarrow R_2+\tfrac12R_3$: $(0,1,0\mid3)$.
>
> **Result.** $\left[\begin{array}{ccc|c}1&0&0&4\\0&1&0&3\\0&0&1&2\end{array}\right]$, so $\mathbf{v}=(4,3,2)^\mathsf{T}\ \mathrm{V}$ with no back substitution. ✓

### Computing the inverse

Since $\mathbf{A}\mathbf{X}=\mathbf{I}$ is $n$ systems $\mathbf{A}\mathbf{x}_j=\mathbf{e}_j$ sharing one matrix, reduce them all at once:

$$
[\mathbf{A}\mid\mathbf{I}]\ \xrightarrow{\ \text{Gauss–Jordan}\ }\ [\mathbf{I}\mid\mathbf{A}^{-1}].
$$

If the left block cannot be reduced to $\mathbf{I}$ (a zero row appears), $\mathbf{A}$ is singular. In terms of elementary matrices: the operations that turn $\mathbf{A}$ into $\mathbf{I}$, accumulated, are $\mathbf{E}_k\cdots\mathbf{E}_1=\mathbf{A}^{-1}$, and the right block records exactly that product applied to $\mathbf{I}$. The cost is about $2n^3$ flops (fewer if one exploits the zeros of $\mathbf{I}$).

> [!example] Worked Example 8.26 — Inverting $\mathbf{G}$ by Gauss–Jordan
> Reduce $[\mathbf{G}\mid\mathbf{I}]$.
>
> **Start.**
> $$
> \left[\begin{array}{ccc|ccc}4&-2&-1&1&0&0\\-2&4&-1&0&1&0\\-1&-1&4&0&0&1\end{array}\right]
> $$
> **Column 1.** $R_1\leftarrow\tfrac14R_1$, $R_2\leftarrow R_2+2R_1$, $R_3\leftarrow R_3+R_1$:
> $$
> \left[\begin{array}{ccc|ccc}1&-\tfrac12&-\tfrac14&\tfrac14&0&0\\0&3&-\tfrac32&\tfrac12&1&0\\0&-\tfrac32&\tfrac{15}4&\tfrac14&0&1\end{array}\right]
> $$
> **Column 2.** $R_2\leftarrow\tfrac13R_2$, $R_1\leftarrow R_1+\tfrac12R_2$, $R_3\leftarrow R_3+\tfrac32R_2$:
> $$
> \left[\begin{array}{ccc|ccc}1&0&-\tfrac12&\tfrac13&\tfrac16&0\\0&1&-\tfrac12&\tfrac16&\tfrac13&0\\0&0&3&\tfrac12&\tfrac12&1\end{array}\right]
> $$
> **Column 3.** $R_3\leftarrow\tfrac13R_3$, $R_1\leftarrow R_1+\tfrac12R_3$, $R_2\leftarrow R_2+\tfrac12R_3$:
> $$
> \left[\begin{array}{ccc|ccc}1&0&0&\tfrac5{12}&\tfrac14&\tfrac16\\0&1&0&\tfrac14&\tfrac5{12}&\tfrac16\\0&0&1&\tfrac16&\tfrac16&\tfrac13\end{array}\right]
> $$
> The right block is $\mathbf{Z}=\mathbf{G}^{-1}$, agreeing with the adjugate calculation of Worked Example 8.12. ✓ Each column of $\mathbf{Z}$ is the vector of node voltages produced by a $1\ \mathrm{mA}$ source at one node.

### Several right-hand sides

Gauss–Jordan (and Gaussian elimination) can process any number of right-hand sides in one pass: reduce $[\mathbf{A}\mid\mathbf{b}_1\ \mathbf{b}_2\ \cdots\ \mathbf{b}_k]$. Each extra column costs only $O(n^2)$. In circuit terms this computes the response to several source patterns at once. In production codes the same economy is obtained more flexibly by storing an $\mathbf{LU}$ factorisation (§8.9), since the right-hand sides are often not all known in advance.

### General solutions of consistent singular systems

When $\operatorname{rank}\mathbf{A}=r<n$, the RREF of $[\mathbf{A}\mid\mathbf{b}]$ splits the unknowns into $r$ **basic** (pivot) variables and $n-r$ **free** variables. Each nonzero row of the RREF contains exactly one basic variable, so it expresses that variable directly in terms of the free ones. Setting the free variables to parameters $s,t,\dots$ gives the solution in the form

$$
\mathbf{x}=\mathbf{x}_p+s\,\mathbf{n}_1+t\,\mathbf{n}_2+\cdots,
$$

where $\mathbf{x}_p$ (all free variables zero) is a particular solution and $\mathbf{n}_1,\mathbf{n}_2,\dots$ are homogeneous solutions forming a basis of the **null space** of $\mathbf{A}$ (Chapter 9). This is the constructive content of the Rouché–Capelli theorem.

> [!example] Worked Example 8.27 — An underdetermined system
> Solve $\mathbf{M}\mathbf{x}=\mathbf{b}$ with $\mathbf{M}$ from Worked Example 8.15 and $\mathbf{b}=(2,2,4)^\mathsf{T}$.
>
> **Step 1 (consistency).** $\mathbf{b}$ satisfies $b_3=b_1+b_2$, the compatibility condition found in Worked Example 8.15, so a solution exists.
>
> **Step 2 (reduce).**
> $$
> \left[\begin{array}{cccc|c}1&2&1&3&2\\2&4&0&2&2\\3&6&1&5&4\end{array}\right]\to
> \left[\begin{array}{cccc|c}1&2&1&3&2\\0&0&-2&-4&-2\\0&0&0&0&0\end{array}\right]\to
> \left[\begin{array}{cccc|c}1&2&0&1&1\\0&0&1&2&1\\0&0&0&0&0\end{array}\right].
> $$
> (The last step is $R_2\leftarrow-\tfrac12R_2$, then $R_1\leftarrow R_1-R_2$.)
>
> **Step 3 (read off).** Basic variables $x_1,x_3$; free variables $x_2=s$, $x_4=t$:
> $$
> x_1=1-2s-t,\qquad x_3=1-2t,
> $$
> $$
> \mathbf{x}=\begin{bmatrix}1\\0\\1\\0\end{bmatrix}+s\begin{bmatrix}-2\\1\\0\\0\end{bmatrix}+t\begin{bmatrix}-1\\0\\-2\\1\end{bmatrix}.
> $$
> **Step 4 (check).** $\mathbf{M}(1,0,1,0)^\mathsf{T}=(2,2,4)^\mathsf{T}$ ✓; $\mathbf{M}(-2,1,0,0)^\mathsf{T}=(0,0,0)^\mathsf{T}$ ✓; $\mathbf{M}(-1,0,-2,1)^\mathsf{T}=(-1-2+3,\ -2+0+2,\ -3-2+5)^\mathsf{T}=\mathbf{0}$ ✓. The two free parameters match $n-r=4-2$.
>
> **Reading.** The RREF also displays the column relations of Worked Example 8.15: the free column 2 reads $(2,0)$, i.e. $\mathbf{m}_2=2\mathbf{m}_1$; the free column 4 reads $(1,2)$, i.e. $\mathbf{m}_4=\mathbf{m}_1+2\mathbf{m}_3$.

> [!warning] Common Pitfall
> Inverting by Gauss–Jordan in floating point still needs pivoting. And having $\mathbf{A}^{-1}$ in hand is rarely useful by itself: if you need $\mathbf{A}^{-1}\mathbf{B}$, solve $\mathbf{A}\mathbf{X}=\mathbf{B}$; if you need one entry of the inverse (say a transfer impedance), solve for one column. Compute the full inverse only when its entries are themselves the deliverable — a nodal impedance matrix for a short-circuit study, or a covariance matrix.

> [!tip] Engineering Note
> In power-system fault analysis the **bus impedance matrix** $\mathbf{Z}_{\text{bus}}=\mathbf{Y}_{\text{bus}}^{-1}$ is genuinely wanted: its diagonal entries are the Thévenin impedances at each bus, which set the short-circuit currents. For large networks it is built not by Gauss–Jordan but by adding one branch at a time with rank-one updates — Sherman–Morrison again — or only the needed columns are computed from an $\mathbf{LU}$ factorisation.

## 8.9 LU decomposition

### Motivation

In practice the same coefficient matrix meets many right-hand sides. A transient simulation with fixed time step solves $(\mathbf{G}+\mathbf{C}/h)\mathbf{v}_{k+1}=\mathbf{i}_{k+1}+(\mathbf{C}/h)\mathbf{v}_k$ at every step with the *same* matrix; a sweep over source positions, a Newton iteration with a frozen Jacobian, or an adjoint sensitivity computation all have the same pattern. Gaussian elimination spends $\tfrac23n^3$ flops on the matrix and only $O(n^2)$ on each right-hand side. **$\mathbf{LU}$ decomposition** packages the expensive part — the record of the elimination — so it can be stored and reused.

### Elimination as factorisation

Each stage of Gaussian elimination without interchanges is a left-multiplication by a unit lower triangular matrix. For $n=3$:

$$
\mathbf{L}_1=\begin{bmatrix}1&0&0\\-m_{21}&1&0\\-m_{31}&0&1\end{bmatrix},\qquad \mathbf{L}_2=\begin{bmatrix}1&0&0\\0&1&0\\0&-m_{32}&1\end{bmatrix},\qquad \mathbf{L}_2\mathbf{L}_1\mathbf{A}=\mathbf{U}.
$$

Inverting, $\mathbf{A}=\mathbf{L}_1^{-1}\mathbf{L}_2^{-1}\mathbf{U}$. The inverses simply flip the signs of the multipliers, and — the pleasant surprise — their product just collects the multipliers in place:

$$
\mathbf{L}=\mathbf{L}_1^{-1}\mathbf{L}_2^{-1}=\begin{bmatrix}1&0&0\\m_{21}&1&0\\m_{31}&m_{32}&1\end{bmatrix}.
$$

> [!theorem] LU factorisation
> If Gaussian elimination on $\mathbf{A}$ ($n\times n$) succeeds without row interchanges, then
> $$
> \mathbf{A}=\mathbf{L}\mathbf{U},
> $$
> where $\mathbf{L}$ is **unit lower triangular** (ones on the diagonal) with $\ell_{ik}=m_{ik}$, the elimination multipliers, and $\mathbf{U}$ is the upper triangular result of elimination. This happens if and only if the leading principal minors $\det\mathbf{A}_{1:k,1:k}$, $k=1,\dots,n-1$, are nonzero, and then (for nonsingular $\mathbf{A}$) the factorisation is unique.
>
> For **any** nonsingular $\mathbf{A}$ there is a permutation matrix $\mathbf{P}$ (recording the row interchanges) such that
> $$
> \mathbf{P}\mathbf{A}=\mathbf{L}\mathbf{U}.
> $$
> With partial pivoting, all $\lvert\ell_{ik}\rvert\le1$.

*Why the leading minors.* The pivots satisfy $u_{11}u_{22}\cdots u_{kk}=\det\mathbf{A}_{1:k,1:k}$, because row replacements using rows $1,\dots,k$ do not change the leading $k\times k$ determinant. So the $k$-th pivot is nonzero exactly when the $k$-th leading minor is. Uniqueness: if $\mathbf{L}_1\mathbf{U}_1=\mathbf{L}_2\mathbf{U}_2$ then $\mathbf{L}_2^{-1}\mathbf{L}_1=\mathbf{U}_2\mathbf{U}_1^{-1}$ is both unit lower and upper triangular, hence $\mathbf{I}$. $\square$

### Solving with the factors

With $\mathbf{P}\mathbf{A}=\mathbf{L}\mathbf{U}$, the system $\mathbf{A}\mathbf{x}=\mathbf{b}$ becomes $\mathbf{L}\mathbf{U}\mathbf{x}=\mathbf{P}\mathbf{b}$, solved in two triangular sweeps:

$$
\begin{aligned}
&\text{forward substitution:} & \mathbf{L}\mathbf{y}&=\mathbf{P}\mathbf{b}, & y_i&=(\mathbf{P}\mathbf{b})_i-\sum_{j<i}\ell_{ij}y_j,\\
&\text{back substitution:} & \mathbf{U}\mathbf{x}&=\mathbf{y}, & x_i&=\frac1{u_{ii}}\Big(y_i-\sum_{j>i}u_{ij}x_j\Big).
\end{aligned}
$$

Each sweep costs about $n^2$ flops. So: **factor once for $\tfrac23n^3$, then solve each right-hand side for $2n^2$.** For $n=1000$ that is $6.7\times10^8$ flops for the factorisation and $2\times10^6$ per solve — a ratio of over 300. The determinant also comes free: $\det\mathbf{A}=\det\mathbf{P}^{-1}\det\mathbf{L}\det\mathbf{U}=\pm\prod_i u_{ii}$.

### Direct (Doolittle) formulas

Equating entries of $\mathbf{A}=\mathbf{L}\mathbf{U}$ row by row gives the factors without writing out elimination steps. For $k=1,\dots,n$:

$$
u_{kj}=a_{kj}-\sum_{p=1}^{k-1}\ell_{kp}u_{pj}\quad(j\ge k),\qquad
\ell_{ik}=\frac1{u_{kk}}\Big(a_{ik}-\sum_{p=1}^{k-1}\ell_{ip}u_{pk}\Big)\quad(i>k).
$$

This is the **Doolittle** algorithm ($\mathbf{L}$ unit diagonal). **Crout**'s variant puts the ones on the diagonal of $\mathbf{U}$ instead. Both perform the same arithmetic as elimination in a different order, which matters for memory access patterns on real hardware but not for the result.

> [!example] Worked Example 8.28 — LU of the running example and reuse
> **Step 1 (factor).** From Worked Example 8.21 the multipliers are $m_{21}=-\tfrac12$, $m_{31}=-\tfrac14$, $m_{32}=-\tfrac12$ and the final triangle is $\mathbf{U}$:
> $$
> \mathbf{G}=\underbrace{\begin{bmatrix}1&0&0\\-\tfrac12&1&0\\-\tfrac14&-\tfrac12&1\end{bmatrix}}_{\mathbf{L}}\underbrace{\begin{bmatrix}4&-2&-1\\0&3&-\tfrac32\\0&0&3\end{bmatrix}}_{\mathbf{U}}.
> $$
> Check row 3: $-\tfrac14(4,-2,-1)-\tfrac12(0,3,-\tfrac32)+(0,0,3)=(-1,\ \tfrac12-\tfrac32,\ \tfrac14+\tfrac34+3)=(-1,-1,4)$ ✓.
>
> **Step 2 (solve for the original sources).** $\mathbf{L}\mathbf{y}=(8,2,1)^\mathsf{T}$: $y_1=8$, $y_2=2+\tfrac12(8)=6$, $y_3=1+\tfrac14(8)+\tfrac12(6)=6$. Then $\mathbf{U}\mathbf{v}=\mathbf{y}$ gives $v_3=2$, $v_2=(6+3)/3=3$, $v_1=(8+6+2)/4=4$ ✓ — the same arithmetic as before, now split into two clean sweeps.
>
> **Step 3 (reuse for a new source pattern).** A $1\ \mathrm{mA}$ source at node 1 only: $\mathbf{L}\mathbf{y}=(1,0,0)^\mathsf{T}$ gives $y=(1,\tfrac12,\tfrac14+\tfrac14)=(1,\tfrac12,\tfrac12)$. Then $v_3=\tfrac12/3=\tfrac16$, $v_2=(\tfrac12+\tfrac32\cdot\tfrac16)/3=(\tfrac12+\tfrac14)/3=\tfrac14$, $v_1=(1+2\cdot\tfrac14+\tfrac16)/4=\tfrac{5/3}{4}=\tfrac5{12}$. This is the first column of $\mathbf{Z}$ from Worked Example 8.12 ✓, obtained for about $2n^2$ operations.
>
> **Step 4.** $\det\mathbf{G}=4\cdot3\cdot3=36$ ✓.

> [!example] Worked Example 8.29 — A factorisation that needs a permutation
> Factor $\mathbf{A}=\begin{bmatrix}1&1&1\\2&2&5\\4&6&8\end{bmatrix}$.
>
> **Step 1.** Eliminate column 1 with $m_{21}=2$, $m_{31}=4$: rows 2 and 3 become $(0,0,3)$ and $(0,2,4)$.
>
> **Step 2.** The $(2,2)$ pivot is $0$, so no $\mathbf{L}\mathbf{U}$ factorisation of $\mathbf{A}$ exists (indeed the leading minor $\begin{vmatrix}1&1\\2&2\end{vmatrix}=0$). Interchange rows 2 and 3 — together with their stored multipliers — to get pivot $2$; column 2 is then already clear below it.
>
> **Step 3.** With $\mathbf{P}$ swapping rows 2 and 3,
> $$
> \mathbf{P}\mathbf{A}=\begin{bmatrix}1&1&1\\4&6&8\\2&2&5\end{bmatrix}=\underbrace{\begin{bmatrix}1&0&0\\4&1&0\\2&0&1\end{bmatrix}}_{\mathbf{L}}\underbrace{\begin{bmatrix}1&1&1\\0&2&4\\0&0&3\end{bmatrix}}_{\mathbf{U}}.
> $$
> Check row 2: $4(1,1,1)+(0,2,4)=(4,6,8)$ ✓; row 3: $2(1,1,1)+(0,0,3)=(2,2,5)$ ✓.
>
> **Step 4.** $\det\mathbf{A}=(-1)^1\cdot1\cdot2\cdot3=-6$. Direct check by cofactors: $1(16-30)-1(16-20)+1(12-8)=-14+4+4=-6$ ✓.
>
> (A library routine with partial pivoting would instead start by swapping row 3, whose entry $4$ is the largest in column 1, to the top; it produces different but equally valid $\mathbf{P},\mathbf{L},\mathbf{U}$ with all $\lvert\ell_{ik}\rvert\le1$.)

### Symmetric matrices: $\mathbf{L}\mathbf{D}\mathbf{L}^\mathsf{T}$ and Cholesky

If $\mathbf{A}$ is symmetric and has an $\mathbf{L}\mathbf{U}$ factorisation, write $\mathbf{U}=\mathbf{D}\tilde{\mathbf{U}}$ with $\mathbf{D}=\operatorname{diag}(u_{11},\dots,u_{nn})$ and $\tilde{\mathbf{U}}$ unit upper triangular. Then $\mathbf{A}=\mathbf{L}\mathbf{D}\tilde{\mathbf{U}}$, and symmetry plus uniqueness force $\tilde{\mathbf{U}}=\mathbf{L}^\mathsf{T}$:

$$
\mathbf{A}=\mathbf{L}\mathbf{D}\mathbf{L}^\mathsf{T}.
$$

If moreover $\mathbf{A}$ is **symmetric positive definite** (SPD: $\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}>0$ for all $\mathbf{x}\ne\mathbf{0}$), all pivots are positive and we can split $\mathbf{D}=\mathbf{D}^{1/2}\mathbf{D}^{1/2}$.

> [!theorem] Cholesky factorisation
> A real symmetric matrix $\mathbf{A}$ is positive definite if and only if it can be written
> $$
> \mathbf{A}=\mathbf{R}\mathbf{R}^\mathsf{T}
> $$
> with $\mathbf{R}$ lower triangular with positive diagonal entries; $\mathbf{R}$ is then unique. (For complex Hermitian positive definite $\mathbf{A}$, $\mathbf{A}=\mathbf{R}\mathbf{R}^\mathsf{H}$.) The entries are computed column by column:
> $$
> r_{kk}=\sqrt{a_{kk}-\sum_{p<k}r_{kp}^2},\qquad r_{ik}=\frac1{r_{kk}}\Big(a_{ik}-\sum_{p<k}r_{ip}r_{kp}\Big)\quad(i>k).
> $$

Cholesky needs about $\tfrac13n^3$ flops — half of $\mathbf{L}\mathbf{U}$ — and half the storage, requires no pivoting for stability, and doubles as a **test for positive definiteness**: if a square root of a non-positive number is ever requested, $\mathbf{A}$ is not SPD. Nodal matrices of passive resistive networks with a ground connection are SPD (§8.2), so Cholesky is the natural direct solver for them.

> [!example] Worked Example 8.30 — $\mathbf{L}\mathbf{D}\mathbf{L}^\mathsf{T}$ and Cholesky for $\mathbf{G}$
> **Step 1 ($\mathbf{L}\mathbf{D}\mathbf{L}^\mathsf{T}$).** From Worked Example 8.28, the pivots are $4,3,3$, so
> $$
> \mathbf{G}=\begin{bmatrix}1&0&0\\-\tfrac12&1&0\\-\tfrac14&-\tfrac12&1\end{bmatrix}\begin{bmatrix}4&0&0\\0&3&0\\0&0&3\end{bmatrix}\begin{bmatrix}1&-\tfrac12&-\tfrac14\\0&1&-\tfrac12\\0&0&1\end{bmatrix}.
> $$
> Indeed $\mathbf{D}\mathbf{L}^\mathsf{T}=\mathbf{U}$: row 2 of $\mathbf{D}\mathbf{L}^\mathsf{T}$ is $3(0,1,-\tfrac12)=(0,3,-\tfrac32)$ ✓.
>
> **Step 2 (Cholesky by the formulas).**
> $$
> \begin{aligned}
> r_{11}&=\sqrt4=2, & r_{21}&=-2/2=-1, & r_{31}&=-1/2,\\
> r_{22}&=\sqrt{4-(-1)^2}=\sqrt3, & r_{32}&=\frac{-1-(-\tfrac12)(-1)}{\sqrt3}=-\frac{3/2}{\sqrt3}=-\frac{\sqrt3}{2}, & &\\
> r_{33}&=\sqrt{4-\tfrac14-\tfrac34}=\sqrt3. & & & &
> \end{aligned}
> $$
> $$
> \mathbf{R}=\begin{bmatrix}2&0&0\\-1&\sqrt3&0\\-\tfrac12&-\tfrac{\sqrt3}2&\sqrt3\end{bmatrix}=\mathbf{L}\mathbf{D}^{1/2}.
> $$
> **Step 3 (check).** Row 3 of $\mathbf{R}$ dotted with itself: $\tfrac14+\tfrac34+3=4=g_{33}$ ✓; row 3 with row 2: $\tfrac12-\tfrac32+0=-1=g_{32}$ ✓.
>
> **Reading.** All square-root arguments were positive, confirming $\mathbf{G}$ is positive definite — the algebraic face of "every resistor absorbs power". Units: $\mathbf{R}$ has units $\mathrm{mS}^{1/2}$.

> [!example] Worked Example 8.31 — Factor once, step many times
> A node with conductance $G=1\ \mathrm{mS}$ to ground and capacitance $C=1\ \mu\mathrm{F}$ to ground, coupled through $2\ \mathrm{mS}$ to a second identical node, is simulated with backward Euler at step $h=1\ \mathrm{ms}$. Find the system matrix and explain the cost per step.
>
> **Step 1.** Backward Euler replaces $\mathbf{C}\dot{\mathbf{v}}$ by $\mathbf{C}(\mathbf{v}_{k+1}-\mathbf{v}_k)/h$, giving $(\mathbf{G}+\mathbf{C}/h)\mathbf{v}_{k+1}=\mathbf{i}_{k+1}+(\mathbf{C}/h)\mathbf{v}_k$. Here $C/h=10^{-6}/10^{-3}=1\ \mathrm{mS}$, so
> $$
> \mathbf{G}+\mathbf{C}/h=\begin{bmatrix}1+2+1&-2\\-2&1+2+1\end{bmatrix}=\begin{bmatrix}4&-2\\-2&4\end{bmatrix}\ \mathrm{mS}.
> $$
> **Step 2.** Factor once: $\ell_{21}=-\tfrac12$, $u_{22}=4-1=3$, so $\mathbf{L}=\begin{bmatrix}1&0\\-\tfrac12&1\end{bmatrix}$, $\mathbf{U}=\begin{bmatrix}4&-2\\0&3\end{bmatrix}$.
>
> **Step 3.** Each time step is now one forward and one back substitution. With a $3\ \mathrm{mA}$ step input at node 1 and $\mathbf{v}_0=\mathbf{0}$: the right-hand side is $(3,0)$; $\mathbf{y}=(3,1.5)$; $v_2=0.5$, $v_1=(3+1)/4=1\ \mathrm{V}$. Next step, right-hand side $(3,0)+1\cdot(1,0.5)=(4,0.5)$, $\mathbf{y}=(4,2.5)$, $v_2=\tfrac56\approx0.833$, $v_1=(4+\tfrac53)/4=\tfrac{17}{12}\approx1.417\ \mathrm{V}$. The voltages approach the DC solution of $\begin{bmatrix}3&-2\\-2&3\end{bmatrix}\mathbf{v}=(3,0)^\mathsf{T}$, namely $(1.8,1.2)\ \mathrm{V}$.
>
> **Cost.** For an $n$-node circuit simulated over $N$ steps the work is $\tfrac23n^3+2Nn^2$ instead of $N\cdot\tfrac23n^3$. Changing $h$ changes the matrix and forces a refactorisation — one reason adaptive-step simulators change $h$ reluctantly.

> [!warning] Common Pitfall
> The $\mathbf{L}$ and $\mathbf{U}$ returned by a library are tied to its permutation: `scipy.linalg.lu` returns $\mathbf{P},\mathbf{L},\mathbf{U}$ with $\mathbf{A}=\mathbf{P}\mathbf{L}\mathbf{U}$, while MATLAB's `[L,U,P]=lu(A)` satisfies $\mathbf{P}\mathbf{A}=\mathbf{L}\mathbf{U}$. Forgetting to apply the permutation to $\mathbf{b}$ gives a wrong answer with no error message. Use the library's paired factor/solve routines (`lu_factor`/`lu_solve`, `decomposition`) rather than handling $\mathbf{P}$ yourself.

> [!tip] Engineering Note
> SPICE's inner loop is: stamp the linearised companion models into the MNA matrix, $\mathbf{LU}$-factor it (sparse, §8.11), forward/back substitute, update the nonlinear devices, repeat until Newton converges, then advance time. Because the sparsity *pattern* never changes, the symbolic part of the factorisation (which entries will be nonzero, which pivot order to use) is computed once at the start of the run, and only the numerical values are refactored at each iteration.

## 8.10 Numerical conditioning

### Motivation

Every number that enters a real linear system carries an error. Component values have tolerances of $1\%$ or $5\%$; measured data carry noise; and even exact inputs are rounded to about 16 significant digits when stored in double precision. The question is how much those small input errors can change the solution. For some systems the answer is "about as much as the inputs changed"; for others a change in the fifth digit of the data changes the first digit of the answer. This sensitivity is a property of the *problem*, not of the algorithm used to solve it, and it is measured by the **condition number**. Distinguishing an ill-conditioned problem from an unstable algorithm (Worked Example 8.23 was the latter) is the central skill of numerical linear algebra.

### Vector and matrix norms

To talk about "small" errors in vectors and matrices we need sizes.

> [!definition] Vector norms
> For $\mathbf{x}\in\mathbb{C}^n$:
> $$
> \lVert\mathbf{x}\rVert_1=\sum_i\lvert x_i\rvert,\qquad \lVert\mathbf{x}\rVert_2=\Big(\sum_i\lvert x_i\rvert^2\Big)^{1/2},\qquad \lVert\mathbf{x}\rVert_\infty=\max_i\lvert x_i\rvert.
> $$
> Each satisfies $\lVert\mathbf{x}\rVert\ge0$ with equality only for $\mathbf{x}=\mathbf{0}$, $\lVert c\mathbf{x}\rVert=\lvert c\rvert\lVert\mathbf{x}\rVert$, and the triangle inequality $\lVert\mathbf{x}+\mathbf{y}\rVert\le\lVert\mathbf{x}\rVert+\lVert\mathbf{y}\rVert$.

> [!definition] Induced matrix norm
> Given a vector norm, the **induced** (operator) norm of $\mathbf{A}$ is the largest factor by which $\mathbf{A}$ can stretch a vector:
> $$
> \lVert\mathbf{A}\rVert=\max_{\mathbf{x}\ne\mathbf{0}}\frac{\lVert\mathbf{A}\mathbf{x}\rVert}{\lVert\mathbf{x}\rVert}.
> $$
> For the three common vector norms this gives
> $$
> \lVert\mathbf{A}\rVert_1=\max_j\sum_i\lvert a_{ij}\rvert\ \ (\text{max column sum}),\qquad
> \lVert\mathbf{A}\rVert_\infty=\max_i\sum_j\lvert a_{ij}\rvert\ \ (\text{max row sum}),
> $$
> $$
> \lVert\mathbf{A}\rVert_2=\sigma_{\max}(\mathbf{A})=\sqrt{\lambda_{\max}(\mathbf{A}^\mathsf{H}\mathbf{A})}\ \ (\text{spectral norm}).
> $$
> The **Frobenius norm** $\lVert\mathbf{A}\rVert_F=\big(\sum_{i,j}\lvert a_{ij}\rvert^2\big)^{1/2}$ is not induced but is easy to compute and satisfies $\lVert\mathbf{A}\rVert_2\le\lVert\mathbf{A}\rVert_F$.

*Why the max-row-sum formula.* For $\lVert\mathbf{x}\rVert_\infty=1$, $\lvert(\mathbf{A}\mathbf{x})_i\rvert\le\sum_j\lvert a_{ij}\rvert\lvert x_j\rvert\le\sum_j\lvert a_{ij}\rvert$, so $\lVert\mathbf{A}\rVert_\infty\le\max_i\sum_j\lvert a_{ij}\rvert$; and choosing $x_j=\operatorname{sgn}a_{pj}$ for the row $p$ attaining the maximum achieves it. The 1-norm case is the same argument on columns. $\square$

Induced norms satisfy $\lVert\mathbf{A}\mathbf{x}\rVert\le\lVert\mathbf{A}\rVert\lVert\mathbf{x}\rVert$ and are **submultiplicative**, $\lVert\mathbf{A}\mathbf{B}\rVert\le\lVert\mathbf{A}\rVert\lVert\mathbf{B}\rVert$, with $\lVert\mathbf{I}\rVert=1$. For the running example $\lVert\mathbf{G}\rVert_1=\lVert\mathbf{G}\rVert_\infty=4+2+1=7\ \mathrm{mS}$ and, $\mathbf{G}$ being symmetric, $\lVert\mathbf{G}\rVert_2$ equals its largest eigenvalue, $6\ \mathrm{mS}$.

### The condition number

> [!definition] Condition number
> For a nonsingular square matrix,
> $$
> \kappa(\mathbf{A})=\lVert\mathbf{A}\rVert\,\lVert\mathbf{A}^{-1}\rVert,
> $$
> with a subscript naming the norm used ($\kappa_1,\kappa_2,\kappa_\infty$). For singular $\mathbf{A}$ one sets $\kappa=\infty$.

Basic facts: $\kappa(\mathbf{A})\ge\lVert\mathbf{A}\mathbf{A}^{-1}\rVert=1$; $\kappa(c\mathbf{A})=\kappa(\mathbf{A})$, so unlike the determinant it is **independent of overall scale and units**; $\kappa_2(\mathbf{A})=\sigma_{\max}/\sigma_{\min}$, the ratio of the largest to smallest singular values; for symmetric positive definite $\mathbf{A}$, $\kappa_2=\lambda_{\max}/\lambda_{\min}$. Geometrically, $\mathbf{A}$ maps the unit sphere to an ellipsoid, and $\kappa_2$ is the ratio of its longest to shortest axis: a very flat ellipsoid means some directions are almost annihilated, and recovering them (inverting) amplifies errors enormously. Moreover $1/\kappa(\mathbf{A})$ is the relative distance from $\mathbf{A}$ to the nearest singular matrix: $\min\{\lVert\delta\mathbf{A}\rVert/\lVert\mathbf{A}\rVert:\ \mathbf{A}+\delta\mathbf{A}\text{ singular}\}=1/\kappa(\mathbf{A})$ in the 2-norm.

### Perturbation theorems

> [!theorem] Sensitivity to the right-hand side
> If $\mathbf{A}\mathbf{x}=\mathbf{b}$ with $\mathbf{b}\ne\mathbf{0}$ and $\mathbf{A}(\mathbf{x}+\delta\mathbf{x})=\mathbf{b}+\delta\mathbf{b}$, then
> $$
> \frac{\lVert\delta\mathbf{x}\rVert}{\lVert\mathbf{x}\rVert}\le\kappa(\mathbf{A})\,\frac{\lVert\delta\mathbf{b}\rVert}{\lVert\mathbf{b}\rVert},
> $$
> and for every $\mathbf{A}$ there are $\mathbf{b},\delta\mathbf{b}$ for which equality holds.

*Proof.* Subtracting, $\mathbf{A}\,\delta\mathbf{x}=\delta\mathbf{b}$, so $\lVert\delta\mathbf{x}\rVert=\lVert\mathbf{A}^{-1}\delta\mathbf{b}\rVert\le\lVert\mathbf{A}^{-1}\rVert\lVert\delta\mathbf{b}\rVert$. Also $\lVert\mathbf{b}\rVert=\lVert\mathbf{A}\mathbf{x}\rVert\le\lVert\mathbf{A}\rVert\lVert\mathbf{x}\rVert$, i.e. $1/\lVert\mathbf{x}\rVert\le\lVert\mathbf{A}\rVert/\lVert\mathbf{b}\rVert$. Multiply the two inequalities. Equality needs $\mathbf{x}$ to be a direction that $\mathbf{A}$ stretches most and $\delta\mathbf{b}$ a direction that $\mathbf{A}^{-1}$ stretches most. $\square$

> [!theorem] Sensitivity to the matrix
> If $(\mathbf{A}+\delta\mathbf{A})(\mathbf{x}+\delta\mathbf{x})=\mathbf{b}$ and $\kappa(\mathbf{A})\lVert\delta\mathbf{A}\rVert/\lVert\mathbf{A}\rVert<1$, then
> $$
> \frac{\lVert\delta\mathbf{x}\rVert}{\lVert\mathbf{x}\rVert}\le\frac{\kappa(\mathbf{A})\,\lVert\delta\mathbf{A}\rVert/\lVert\mathbf{A}\rVert}{1-\kappa(\mathbf{A})\,\lVert\delta\mathbf{A}\rVert/\lVert\mathbf{A}\rVert}\approx\kappa(\mathbf{A})\frac{\lVert\delta\mathbf{A}\rVert}{\lVert\mathbf{A}\rVert}.
> $$

*Sketch.* Expand: $\mathbf{A}\,\delta\mathbf{x}=-\delta\mathbf{A}(\mathbf{x}+\delta\mathbf{x})$, so $\lVert\delta\mathbf{x}\rVert\le\lVert\mathbf{A}^{-1}\rVert\lVert\delta\mathbf{A}\rVert(\lVert\mathbf{x}\rVert+\lVert\delta\mathbf{x}\rVert)$; collect the $\lVert\delta\mathbf{x}\rVert$ terms and divide by $\lVert\mathbf{x}\rVert$. $\square$

**Rule of thumb.** If the data are accurate to $d$ significant digits and $\kappa\approx10^k$, expect the solution to be accurate to about $d-k$ digits. With double precision ($d\approx16$) a system with $\kappa=10^{10}$ still yields about 6 good digits; with component values known to $1\%$ ($d=2$) even $\kappa=100$ leaves nothing trustworthy in the worst case.

### Residual versus error

Having computed $\hat{\mathbf{x}}$, we can always evaluate the **residual** $\mathbf{r}=\mathbf{b}-\mathbf{A}\hat{\mathbf{x}}$, but not the **error** $\mathbf{e}=\mathbf{x}-\hat{\mathbf{x}}$ (we do not know $\mathbf{x}$). Since $\mathbf{A}\mathbf{e}=\mathbf{r}$, the first perturbation theorem applies with $\delta\mathbf{b}=\mathbf{r}$:

$$
\frac{\lVert\mathbf{e}\rVert}{\lVert\mathbf{x}\rVert}\le\kappa(\mathbf{A})\frac{\lVert\mathbf{r}\rVert}{\lVert\mathbf{b}\rVert}.
$$

A small residual guarantees a small error only when $\kappa$ is modest. Gaussian elimination with partial pivoting is **backward stable**: it always produces a small residual, $\lVert\mathbf{r}\rVert\lesssim n\rho u\lVert\mathbf{A}\rVert\lVert\hat{\mathbf{x}}\rVert$. So for a well-conditioned system its answer is accurate, and for an ill-conditioned one its answer is as accurate as the problem allows — no algorithm can do substantially better with the same data.

```desmos-graph
left=-1; right=3; bottom=-0.5; top=2.5
---
x+y=2|blue
x+1.05y=2.05|red
x+1.05y=2.1|red|dashed
(1,1)|label:original solution
(0,2)|label:perturbed solution
```

*An ill-conditioned $2\times2$ system: nearly parallel lines $x+y=2$ (blue) and $x+1.05y=2.05$ (red) meet at $(1,1)$. Raising the second right-hand side by only $0.05$ (dashed red) slides the intersection along the blue line to $(0,2)$. The narrow angle between the lines is what a large condition number looks like.*

> [!example] Worked Example 8.32 — Condition number of the running example
> Compute $\kappa_1(\mathbf{G})$ and $\kappa_2(\mathbf{G})$.
>
> **Step 1.** $\lVert\mathbf{G}\rVert_1=\max(4+2+1,\ 2+4+1,\ 1+1+4)=7$.
>
> **Step 2.** From Worked Example 8.12, the column sums of $\lvert\mathbf{Z}\rvert$ are $\tfrac{30}{36},\tfrac{30}{36},\tfrac{24}{36}$, so $\lVert\mathbf{Z}\rVert_1=\tfrac56$.
>
> **Step 3.** $\kappa_1(\mathbf{G})=7\cdot\tfrac56=\tfrac{35}6\approx5.83$.
>
> **Step 4 (2-norm).** $\mathbf{G}$ is symmetric positive definite with eigenvalues $3-\sqrt3\approx1.268$, $3+\sqrt3\approx4.732$ and $6$ (Chapter 10). So $\kappa_2=6/(3-\sqrt3)=3+\sqrt3\approx4.73$.
>
> **Reading.** $\kappa\approx5$: a $1\%$ uncertainty in the source currents produces at most about $5\%$ relative uncertainty in the node voltages. Diagonally dominant nodal matrices of modest size are typically this benign.

> [!example] Worked Example 8.33 — Nearly parallel lines
> For $\mathbf{A}=\begin{bmatrix}1&1\\1&1.05\end{bmatrix}$, $\mathbf{b}=(2,\ 2.05)^\mathsf{T}$, verify the perturbation bound in the $\infty$-norm for $\delta\mathbf{b}=(0,\ 0.05)^\mathsf{T}$, and show that a small residual can hide a large error.
>
> **Step 1.** $\det\mathbf{A}=0.05$, $\mathbf{A}^{-1}=\dfrac1{0.05}\begin{bmatrix}1.05&-1\\-1&1\end{bmatrix}=\begin{bmatrix}21&-20\\-20&20\end{bmatrix}$.
>
> **Step 2.** $\lVert\mathbf{A}\rVert_\infty=2.05$, $\lVert\mathbf{A}^{-1}\rVert_\infty=41$, so $\kappa_\infty=84.05$.
>
> **Step 3.** $\mathbf{x}=\mathbf{A}^{-1}\mathbf{b}=(42-41,\ -40+41)^\mathsf{T}=(1,1)^\mathsf{T}$. $\delta\mathbf{x}=\mathbf{A}^{-1}\delta\mathbf{b}=(-1,\ 1)^\mathsf{T}$, so $\mathbf{x}+\delta\mathbf{x}=(0,2)^\mathsf{T}$.
>
> **Step 4 (bound).** Relative input change $\lVert\delta\mathbf{b}\rVert_\infty/\lVert\mathbf{b}\rVert_\infty=0.05/2.05\approx0.0244$; relative output change $\lVert\delta\mathbf{x}\rVert_\infty/\lVert\mathbf{x}\rVert_\infty=1$, i.e. $100\%$. The amplification is $1/0.0244=41$, within the bound $84.05$ ✓.
>
> **Step 5 (residual).** Suppose some method returned $\hat{\mathbf{x}}=(0,2)^\mathsf{T}$ for the original system. Its residual $\mathbf{b}-\mathbf{A}\hat{\mathbf{x}}=(2-2,\ 2.05-2.1)^\mathsf{T}=(0,-0.05)^\mathsf{T}$ is about $2\%$ of $\lVert\mathbf{b}\rVert$, yet the error is $100\%$. Checking that "the answer satisfies the equations" is not enough for an ill-conditioned system.

> [!example] Worked Example 8.34 — Wide spread of conductances
> Two nodes are joined by a low-resistance strap of $1000\ \mathrm{mS}$ ($1\ \Omega$), and each has $1\ \mathrm{mS}$ to ground. Analyse
> $$
> \mathbf{G}_s=\begin{bmatrix}1001&-1000\\-1000&1001\end{bmatrix}\ \mathrm{mS}.
> $$
> **Step 1.** The eigenvalues are $1001\mp1000$, i.e. $1$ and $2001$, with eigenvectors $(1,1)$ (both nodes moving together, opposed only by the ground conductances) and $(1,-1)$ (the nodes moving apart, opposed by the strap). Hence $\kappa_2=2001$.
>
> **Step 2.** With $1\ \mathrm{mA}$ into node 1, $\mathbf{v}=\tfrac1{2001}(1001,1000)^\mathsf{T}\approx(0.50025,\ 0.49975)\ \mathrm{V}$.
>
> **Step 3.** The *sum* $v_1+v_2$ is well determined, but the *difference* $v_1-v_2=1/2001\approx0.5\ \mathrm{mV}$ — the voltage across the strap — is the small result of subtracting two nearly equal numbers. A $1\%$ change in the strap conductance changes $v_1-v_2$ by about $1\%$, while a $10^{-6}$ relative rounding error in $v_1$ and $v_2$ changes it by about $0.1\%$: the relative accuracy lost is roughly the factor $\kappa$.
>
> **Reading.** Circuits with conductances spanning many decades (milliohm interconnect next to megaohm bias resistors, or the $\mathrm{G}\Omega$ "gmin" conductances SPICE adds to ground) have large condition numbers. Engineers often replace an extremely stiff element by an ideal short (merging the nodes), which removes the ill-conditioning analytically.

> [!example] Worked Example 8.35 — Scaling can remove apparent ill-conditioning
> The system of Worked Example 8.23, with its first equation multiplied by $10^4$, reads
> $$
> \begin{bmatrix}1&10^4\\1&1\end{bmatrix}\mathbf{x}=\begin{bmatrix}10^4\\2\end{bmatrix}.
> $$
> **Step 1.** $\lVert\mathbf{A}\rVert_\infty=10001$, $\mathbf{A}^{-1}=\dfrac{1}{1-10^4}\begin{bmatrix}1&-10^4\\-1&1\end{bmatrix}$, $\lVert\mathbf{A}^{-1}\rVert_\infty=10001/9999\approx1.0002$, so $\kappa_\infty\approx1.0\times10^4$.
>
> **Step 2.** Dividing row 1 by $10^4$ (**row equilibration**, making the largest entry in each row equal to 1) gives back $\begin{bmatrix}10^{-4}&1\\1&1\end{bmatrix}$ with $\kappa_\infty\approx4.0$.
>
> **Reading.** The solution $(1.0001,\ 0.9999)$ is the same either way; only the bookkeeping changed. Part of a large condition number may therefore come from poor scaling — for example, MNA equations mixing currents in amperes with voltage constraints in volts, or conductances in siemens next to capacitances in farads — and is cured by rescaling rows or columns (equilibration), which good solvers do automatically. The part that survives every scaling, as in Worked Example 8.33, is intrinsic to the problem.

### A notoriously ill-conditioned family

The **Hilbert matrix** $\mathbf{H}_n$, with $h_{ij}=1/(i+j-1)$, arises in least-squares fitting of polynomials $1,t,t^2,\dots$ on $[0,1]$ — a reason never to fit high-degree polynomials in the monomial basis. Its condition numbers grow explosively:

| $n$ | 3 | 5 | 8 | 10 |
|---|---|---|---|---|
| $\kappa_2(\mathbf{H}_n)$ | $5.2\times10^{2}$ | $4.8\times10^{5}$ | $1.5\times10^{10}$ | $1.6\times10^{13}$ |

At $n=10$ a double-precision solve keeps only about three correct digits. The cure is a better-conditioned formulation (orthogonal polynomials, QR factorisation — Chapter 9), not a cleverer elimination.

> [!warning] Common Pitfall
> Do not judge conditioning by the determinant (§8.3) or by the size of the entries. And do not confuse an **ill-conditioned problem** (large $\kappa$: every method struggles) with an **unstable algorithm** (elimination without pivoting on Worked Example 8.23: a good method would have succeeded). The remedies are different: reformulate or measure better for the first; pivot or use a stable factorisation for the second.

> [!tip] Engineering Note
> Libraries estimate $\kappa_1$ cheaply from an existing $\mathbf{LU}$ factorisation (LAPACK `dgecon`, MATLAB `rcond`, which returns $1/\kappa_1$); MATLAB's warning "Matrix is close to singular or badly scaled. RCOND = 1e-17" is this estimate. In tolerance analysis the condition number gives a quick worst-case bound before running a Monte Carlo study: if $\kappa\approx10^3$ and components are $1\%$ parts, some outputs may be completely unreliable, and the design should be changed (for instance, a bridge measurement should be read as a difference directly rather than computed from two large node voltages).

## 8.11 Sparse matrices

### Motivation

A circuit with $n$ nodes has an $n\times n$ nodal matrix, but each node connects to only a few others — typically three to five. The matrix therefore has about $5n$ nonzero entries out of $n^2$. For a post-layout simulation with $n=10^6$ nodes, the dense matrix would have $10^{12}$ entries (8 terabytes) and dense elimination would need $\tfrac23\times10^{18}$ flops; the sparse matrix has perhaps $5\times10^6$ nonzeros (tens of megabytes) and can be factored in seconds. Finite-difference and finite-element field solvers, power-grid models and the normal equations of large estimation problems are all sparse in the same way. **Sparse matrix methods** store and operate on only the nonzeros, and choose algorithms that keep it that way.

> [!definition] Sparse matrix
> A matrix is **sparse** when the number of nonzero entries, $\operatorname{nnz}(\mathbf{A})$, is small enough compared with $mn$ that exploiting the zeros pays off — in practice $\operatorname{nnz}=O(n)$ for an $n\times n$ matrix. The **sparsity pattern** is the set of positions of the nonzeros. A matrix is **banded** with bandwidth $w$ if $a_{ij}=0$ for $\lvert i-j\rvert>w$; **tridiagonal** is $w=1$.

The sparsity pattern of a nodal matrix *is* the circuit graph: $g_{ij}\ne0$ ($i\ne j$) exactly when nodes $i$ and $j$ are joined by a branch.

### Storage formats

- **Coordinate (COO, "triplet")**: three arrays of length $\operatorname{nnz}$ holding (row, column, value). Easy to build by stamping — duplicates are simply summed afterwards — and the usual assembly format.
- **Compressed sparse row (CSR)**: the values and column indices row by row, plus a row-pointer array of length $n+1$ marking where each row starts. Storage $2\operatorname{nnz}+n+1$; matrix–vector products touch each nonzero once ($2\operatorname{nnz}$ flops).
- **Compressed sparse column (CSC)**: the same by columns; the native format of MATLAB and most direct solvers.
- **Banded storage**: the $2w+1$ diagonals as dense vectors — ideal for ladder networks and 1-D discretisations.

> [!example] Worked Example 8.36 — CSR storage of a ladder network
> A four-node resistive ladder has $1\ \mathrm{S}$ series links between consecutive nodes and $1\ \mathrm{S}$ to ground at every node, with the left end driven through a further $1\ \mathrm{S}$ from a $10\ \mathrm{V}$ source (Norton equivalent: $10\ \mathrm{A}$ into node 1). Its nodal matrix is
> $$
> \mathbf{T}=\begin{bmatrix}3&-1&0&0\\-1&3&-1&0\\0&-1&3&-1\\0&0&-1&2\end{bmatrix},\qquad \mathbf{b}=\begin{bmatrix}10\\0\\0\\0\end{bmatrix}.
> $$
> Give its CSR representation (0-based indices).
>
> **Solution.** Reading row by row,
>
> | Array | Contents |
> |---|---|
> | values | $3,-1,\ -1,3,-1,\ -1,3,-1,\ -1,2$ |
> | column index | $0,1,\ 0,1,2,\ 1,2,3,\ 2,3$ |
> | row pointer | $0,2,5,8,10$ |
>
> Row $i$ occupies positions $\text{ptr}[i]$ to $\text{ptr}[i+1]-1$; e.g. row 2 is positions 5–7: values $-1,3,-1$ in columns $1,2,3$. Ten nonzeros out of sixteen entries here; for a 1000-node ladder it would be $2998$ out of $10^6$.

### Fill-in and ordering

Eliminating a variable creates nonzeros wherever two of its neighbours were not already connected: the update $a_{ij}\leftarrow a_{ij}-a_{ik}a_{kj}/a_{kk}$ makes $a_{ij}$ nonzero whenever $a_{ik}$ and $a_{kj}$ are. In graph terms, eliminating node $k$ joins all its neighbours into a clique — exactly the star–mesh transformation of circuit theory. These new nonzeros are **fill-in**. The factors $\mathbf{L}$ and $\mathbf{U}$ can be far denser than $\mathbf{A}$, and the amount of fill depends dramatically on the **elimination order** (the symmetric permutation $\mathbf{P}\mathbf{A}\mathbf{P}^\mathsf{T}$).

> [!example] Worked Example 8.37 — A star network: order is everything
> A hub node connects through $1\ \mathrm{S}$ to each of four leaf nodes, and each leaf has $1\ \mathrm{S}$ to ground. Compare elimination with the hub numbered first and last.
>
> **Hub first.**
> $$
> \mathbf{A}=\begin{bmatrix}4&-1&-1&-1&-1\\-1&2&0&0&0\\-1&0&2&0&0\\-1&0&0&2&0\\-1&0&0&0&2\end{bmatrix}.
> $$
> Eliminating the hub replaces the leaf block by $2\mathbf{I}-\tfrac14\mathbf{1}\mathbf{1}^\mathsf{T}$ (where $\mathbf{1}$ is the all-ones vector): every one of its 12 off-diagonal zeros becomes $-\tfrac14$. The factor $\mathbf{L}$ is completely full — the eliminated hub has "meshed" all four leaves together.
>
> **Hub last.** The matrix is the "arrowhead" $\begin{bmatrix}2\mathbf{I}_4&-\mathbf{1}\\-\mathbf{1}^\mathsf{T}&4\end{bmatrix}$. Eliminating leaf $k$ affects only the hub's diagonal: $4\leftarrow4-(-1)(-1)/2$. After all four leaves the hub pivot is $4-4\cdot\tfrac12=2$, and **no fill-in** occurs: $\mathbf{L}$ has exactly the nonzeros of the lower triangle of $\mathbf{A}$.
>
> **Check.** Both orders must give the same determinant: hub last, $2^4\cdot2=32$; hub first, $4\cdot\det(2\mathbf{I}-\tfrac14\mathbf{1}\mathbf{1}^\mathsf{T})=4\cdot2^3(2-\tfrac44)=4\cdot8\cdot1=32$ ✓ (using $\det(c\mathbf{I}-\alpha\mathbf{1}\mathbf{1}^\mathsf{T})=c^{n-1}(c-n\alpha)$).
>
> **Reading.** Eliminate low-degree nodes first and high-degree nodes (hubs, supply rails) last. This greedy rule is the **minimum-degree** ordering; its refinements (approximate minimum degree, AMD) and **nested dissection** (recursively cutting the graph with small separators, which are numbered last) are what production sparse solvers use.

For a nodal matrix that is SPD and diagonally dominant, any symmetric ordering is numerically safe, so the ordering can be chosen purely to minimise fill. For MNA and other indefinite matrices the solver must balance sparsity against pivoting for stability (threshold pivoting: accept a pivot if it is not much smaller than the largest candidate).

### Banded and tridiagonal systems

If $\mathbf{A}$ has bandwidth $w$, elimination without pivoting creates no fill outside the band, and the factorisation costs about $2nw^2$ flops instead of $\tfrac23n^3$. The tridiagonal case ($w=1$) is the **Thomas algorithm**: for the system with subdiagonal $a_i$, diagonal $d_i$, superdiagonal $c_i$, right-hand side $b_i$,

$$
\begin{aligned}
&\text{forward:} & c_1'&=\frac{c_1}{d_1},\quad b_1'=\frac{b_1}{d_1},\quad c_i'=\frac{c_i}{d_i-a_ic_{i-1}'},\quad b_i'=\frac{b_i-a_ib_{i-1}'}{d_i-a_ic_{i-1}'},\\
&\text{backward:} & x_n&=b_n',\quad x_i=b_i'-c_i'x_{i+1}.
\end{aligned}
$$

About $8n$ flops in total — linear in $n$. It is stable without pivoting when the matrix is diagonally dominant.

> [!example] Worked Example 8.38 — Solving the ladder by the Thomas algorithm
> Solve $\mathbf{T}\mathbf{v}=\mathbf{b}$ from Worked Example 8.36. Here $d=(3,3,3,2)$, $a_i=c_i=-1$.
>
> **Forward sweep.**
> $$
> \begin{aligned}
> c_1'&=-\tfrac13, & b_1'&=\tfrac{10}3,\\
> d_2-a_2c_1'&=3-\tfrac13=\tfrac83, & c_2'&=-\tfrac38, & b_2'&=\frac{0+\tfrac{10}3}{8/3}=\tfrac54,\\
> d_3-a_3c_2'&=3-\tfrac38=\tfrac{21}8, & c_3'&=-\tfrac8{21}, & b_3'&=\frac{0+\tfrac54}{21/8}=\tfrac{10}{21},\\
> d_4-a_4c_3'&=2-\tfrac8{21}=\tfrac{34}{21}, & & & b_4'&=\frac{0+\tfrac{10}{21}}{34/21}=\tfrac5{17}.
> \end{aligned}
> $$
> **Backward sweep.** $v_4=\tfrac5{17}$; $v_3=\tfrac{10}{21}+\tfrac8{21}\cdot\tfrac5{17}=\tfrac{210}{357}=\tfrac{10}{17}$; $v_2=\tfrac54+\tfrac38\cdot\tfrac{10}{17}=\tfrac{25}{17}$; $v_1=\tfrac{10}3+\tfrac13\cdot\tfrac{25}{17}=\tfrac{65}{17}$.
>
> So $\mathbf{v}=\tfrac1{17}(65,25,10,5)^\mathsf{T}\approx(3.824,\ 1.471,\ 0.588,\ 0.294)\ \mathrm{V}$.
>
> **Check.** Row 1: $3\cdot65-25=170=17\cdot10$ ✓. Row 4: $-10+2\cdot5=0$ ✓. The pivots $3,\tfrac83,\tfrac{21}8,\tfrac{34}{21}$ multiply to $\det\mathbf{T}=34$. Successive node voltages fall by factors $2.6$, $2.5$ and $2$ (the last is an end effect of the open termination); far from the ends of a long ladder the ratio approaches $r$ with $r+1/r=3$, i.e. $(3+\sqrt5)/2\approx2.618$ of an infinite ladder, the same number met in Worked Example 8.6.

### Iterative methods

For very large sparse systems — especially 3-D field problems, where even the best ordering produces heavy fill — **iterative methods** replace factorisation by repeated cheap matrix–vector products. Split $\mathbf{A}=\mathbf{D}+\mathbf{L}_A+\mathbf{U}_A$ into diagonal, strictly lower and strictly upper parts.

> [!definition] Jacobi and Gauss–Seidel iterations
> Starting from a guess $\mathbf{x}^{(0)}$,
> $$
> \text{Jacobi:}\quad x_i^{(k+1)}=\frac1{a_{ii}}\Big(b_i-\sum_{j\ne i}a_{ij}x_j^{(k)}\Big),
> $$
> $$
> \text{Gauss–Seidel:}\quad x_i^{(k+1)}=\frac1{a_{ii}}\Big(b_i-\sum_{j<i}a_{ij}x_j^{(k+1)}-\sum_{j>i}a_{ij}x_j^{(k)}\Big).
> $$
> In matrix form $\mathbf{x}^{(k+1)}=\mathbf{M}\mathbf{x}^{(k)}+\mathbf{c}$ with $\mathbf{M}_J=-\mathbf{D}^{-1}(\mathbf{L}_A+\mathbf{U}_A)$ or $\mathbf{M}_{GS}=-(\mathbf{D}+\mathbf{L}_A)^{-1}\mathbf{U}_A$.

Gauss–Seidel uses each new value as soon as it is available. In circuit terms, Jacobi sets every node voltage to the conductance-weighted average of its neighbours' old voltages plus its source contribution — a relaxation process.

> [!theorem] Convergence of stationary iterations
> The iteration $\mathbf{x}^{(k+1)}=\mathbf{M}\mathbf{x}^{(k)}+\mathbf{c}$ converges for every starting vector if and only if the **spectral radius** $\rho(\mathbf{M})=\max\lvert\lambda_i(\mathbf{M})\rvert<1$; the error then shrinks asymptotically by the factor $\rho(\mathbf{M})$ per iteration. Sufficient conditions: if $\mathbf{A}$ is strictly diagonally dominant, both Jacobi and Gauss–Seidel converge; if $\mathbf{A}$ is symmetric positive definite, Gauss–Seidel converges.

*Sketch.* The error $\mathbf{e}^{(k)}=\mathbf{x}^{(k)}-\mathbf{x}$ obeys $\mathbf{e}^{(k+1)}=\mathbf{M}\mathbf{e}^{(k)}$, so $\mathbf{e}^{(k)}=\mathbf{M}^k\mathbf{e}^{(0)}$, and $\mathbf{M}^k\to\mathbf{0}$ exactly when $\rho(\mathbf{M})<1$ (Chapter 10). For diagonal dominance, $\lVert\mathbf{M}_J\rVert_\infty=\max_i\sum_{j\ne i}\lvert a_{ij}\rvert/\lvert a_{ii}\rvert<1$ directly. $\square$

> [!example] Worked Example 8.39 — Jacobi and Gauss–Seidel on the running example
> Solve $\mathbf{G}\mathbf{v}=(8,2,1)^\mathsf{T}$ from $\mathbf{v}^{(0)}=\mathbf{0}$.
>
> **Jacobi, iteration 1.** $v_1=\tfrac14(8+2\cdot0+0)=2$, $v_2=\tfrac14(2+0+0)=0.5$, $v_3=\tfrac14(1+0+0)=0.25$.
>
> **Jacobi, iteration 2.** $v_1=\tfrac14(8+2(0.5)+0.25)=2.3125$, $v_2=\tfrac14(2+2(2)+0.25)=1.5625$, $v_3=\tfrac14(1+2+0.5)=0.875$.
>
> **Gauss–Seidel, iteration 1.** $v_1=\tfrac14(8)=2$; then, using the new $v_1$, $v_2=\tfrac14(2+2\cdot2+0)=1.5$; then $v_3=\tfrac14(1+2+1.5)=1.125$.
>
> **Gauss–Seidel, iteration 2.** $v_1=\tfrac14(8+3+1.125)=3.03125$, $v_2=\tfrac14(2+6.0625+1.125)=2.296875$, $v_3=\tfrac14(1+3.03125+2.296875)=1.58203$.
>
> **Progress.** Maximum error $\lVert\mathbf{v}^{(k)}-(4,3,2)^\mathsf{T}\rVert_\infty$:
>
> | $k$ | 1 | 2 | 3 | 4 | 5 | 6 |
> |---|---|---|---|---|---|---|
> | Jacobi | 2.500 | 1.688 | 1.125 | 0.758 | 0.512 | 0.347 |
> | Gauss–Seidel | 2.000 | 0.969 | 0.456 | 0.216 | 0.102 | 0.048 |
>
> **Rates.** $\rho(\mathbf{M}_J)=(1+\sqrt3)/4\approx0.683$ and $\rho(\mathbf{M}_{GS})\approx0.472$: Gauss–Seidel removes roughly twice as many digits per iteration. Reaching $10^{-6}$ accuracy takes about $\ln(10^{-6}/4)/\ln0.683\approx40$ Jacobi or $21$ Gauss–Seidel iterations — far more work than a $3\times3$ elimination, which is why iterative methods are only worthwhile for large sparse systems.

```desmos-graph
left=0; right=9; bottom=0; top=4.2
---
y=4(0.683)^{x}|blue|dashed
y=4(0.4724)^{x}|red|dashed
(1,2.5)|blue
(2,1.6875)|blue
(3,1.125)|blue
(4,0.758)|blue
(5,0.512)|blue
(6,0.347)|blue
(7,0.235)|blue
(8,0.160)|blue
(1,2.0)|red
(2,0.969)|red
(3,0.456)|red
(4,0.216)|red
(5,0.102)|red
(6,0.048)|red
```

*Maximum error against iteration number $k$ for Jacobi (blue) and Gauss–Seidel (red) on the running example. The dashed curves $4\rho^k$ with $\rho=0.683$ and $0.472$ show the geometric decay set by the spectral radius of each iteration matrix.*

Modern practice goes well beyond these classical methods. **Successive over-relaxation** (SOR) extrapolates the Gauss–Seidel step by a factor $1<\omega<2$; the **conjugate gradient** method (for SPD matrices) and **GMRES** (general matrices) are Krylov-subspace methods that converge in a number of iterations governed by $\sqrt{\kappa}$ or by the eigenvalue distribution; and **preconditioning** — solving $\mathbf{P}^{-1}\mathbf{A}\mathbf{x}=\mathbf{P}^{-1}\mathbf{b}$ with $\mathbf{P}\approx\mathbf{A}$ cheap to invert, such as an incomplete $\mathbf{LU}$ factorisation that discards fill — is what makes them fast. These are treated in Chapter 22.

> [!warning] Common Pitfall
> Never form the inverse, or any dense intermediate, of a sparse matrix: $\mathbf{A}^{-1}$ is generally full (Worked Example 8.12 shows a dense $\mathbf{Z}$). Likewise, building a sparse matrix by inserting entries one at a time into a CSR/CSC structure is quadratic in cost; assemble in COO (triplet) form and convert once. And beware of accidental densification: in NumPy/SciPy, `A + 0` or `A * dense_array` may silently return a dense result.

> [!tip] Engineering Note
> SPICE-class simulators use sparse $\mathbf{LU}$ with a Markowitz-type ordering chosen once per circuit (KLU, designed specifically for circuit matrices, exploits their near-block-triangular structure). Field solvers using finite elements produce SPD sparse matrices solved by preconditioned conjugate gradients or multigrid. Power-system load flow uses sparse $\mathbf{LU}$ of the Jacobian with ordering by minimum degree — introduced for exactly this purpose by Tinney and Walker in 1967. In every case the decisive step is the same: recognise the graph behind the matrix and order the unknowns to keep the factors sparse.

## Summary

| Concept | Key formula or fact |
|---|---|
| Product | $c_{ij}=\sum_k a_{ik}b_{kj}$; $\mathbf{A}\mathbf{x}=\sum_k x_k\mathbf{a}_k$; $\mathbf{A}\mathbf{B}\ne\mathbf{B}\mathbf{A}$ in general |
| Transpose rules | $(\mathbf{A}\mathbf{B})^\mathsf{T}=\mathbf{B}^\mathsf{T}\mathbf{A}^\mathsf{T}$, $(\mathbf{A}\mathbf{B})^\mathsf{H}=\mathbf{B}^\mathsf{H}\mathbf{A}^\mathsf{H}$ |
| Nodal matrix | $\mathbf{G}=\mathbf{A}\mathbf{G}_b\mathbf{A}^\mathsf{T}$ (symmetric, positive definite if grounded) |
| $2\times2$ determinant / inverse | $ad-bc$; $\dfrac{1}{ad-bc}\begin{bmatrix}d&-b\\-c&a\end{bmatrix}$ |
| Cofactor expansion | $\det\mathbf{A}=\sum_k a_{ik}C_{ik}$, $C_{ik}=(-1)^{i+k}M_{ik}$ |
| Determinant rules | $\det(\mathbf{A}\mathbf{B})=\det\mathbf{A}\det\mathbf{B}$; swap $\to-$; replacement $\to$ unchanged; $\det=\pm\prod u_{ii}$ |
| Cramer's rule | $x_k=\det\mathbf{A}_k/\det\mathbf{A}$ |
| Adjugate inverse | $\mathbf{A}^{-1}=\operatorname{adj}\mathbf{A}/\det\mathbf{A}$ |
| Product inverse | $(\mathbf{A}\mathbf{B})^{-1}=\mathbf{B}^{-1}\mathbf{A}^{-1}$ |
| Sherman–Morrison | $(\mathbf{A}+\mathbf{u}\mathbf{v}^\mathsf{T})^{-1}=\mathbf{A}^{-1}-\dfrac{\mathbf{A}^{-1}\mathbf{u}\mathbf{v}^\mathsf{T}\mathbf{A}^{-1}}{1+\mathbf{v}^\mathsf{T}\mathbf{A}^{-1}\mathbf{u}}$ |
| Schur complement (Kron) | $\mathbf{S}=\mathbf{A}_{11}-\mathbf{A}_{12}\mathbf{A}_{22}^{-1}\mathbf{A}_{21}$, $\det\mathbf{A}=\det\mathbf{A}_{22}\det\mathbf{S}$ |
| Rank | number of pivots = independent rows = independent columns |
| Rouché–Capelli | consistent $\iff\operatorname{rank}\mathbf{A}=\operatorname{rank}[\mathbf{A}\mid\mathbf{b}]$; $n-r$ free parameters |
| Solution structure | $\mathbf{x}=\mathbf{x}_p+\mathbf{x}_h$, $\mathbf{A}\mathbf{x}_h=\mathbf{0}$ |
| Gaussian elimination | $m_{ik}=a_{ik}/a_{kk}$; cost $\tfrac23n^3$; back substitution $n^2$ |
| Partial pivoting | largest $\lvert a_{pk}\rvert$ as pivot; $\lvert m_{ik}\rvert\le1$ |
| Gauss–Jordan | $[\mathbf{A}\mid\mathbf{I}]\to[\mathbf{I}\mid\mathbf{A}^{-1}]$; cost $\approx n^3$ per solve |
| LU | $\mathbf{P}\mathbf{A}=\mathbf{L}\mathbf{U}$; $\mathbf{L}\mathbf{y}=\mathbf{P}\mathbf{b}$, $\mathbf{U}\mathbf{x}=\mathbf{y}$; $2n^2$ per extra right-hand side |
| Cholesky | $\mathbf{A}=\mathbf{R}\mathbf{R}^\mathsf{T}$ for SPD; $r_{kk}=\sqrt{a_{kk}-\sum_{p<k}r_{kp}^2}$; cost $\tfrac13n^3$ |
| Induced norms | $\lVert\mathbf{A}\rVert_1$ = max column sum, $\lVert\mathbf{A}\rVert_\infty$ = max row sum, $\lVert\mathbf{A}\rVert_2=\sigma_{\max}$ |
| Condition number | $\kappa=\lVert\mathbf{A}\rVert\lVert\mathbf{A}^{-1}\rVert$; $\kappa_2=\sigma_{\max}/\sigma_{\min}$ |
| Perturbation bound | $\lVert\delta\mathbf{x}\rVert/\lVert\mathbf{x}\rVert\le\kappa\,\lVert\delta\mathbf{b}\rVert/\lVert\mathbf{b}\rVert$; lose $\approx\log_{10}\kappa$ digits |
| Thomas algorithm | tridiagonal solve in $\approx8n$ flops |
| Stationary iteration | converges $\iff\rho(\mathbf{M})<1$; diagonal dominance suffices |

## Further reading

- E. Kreyszig, *Advanced Engineering Mathematics*, 10th ed., Wiley — Ch. 7 (matrices, determinants, linear systems) and Ch. 20 (numerical linear algebra).
- G. Strang, *Introduction to Linear Algebra*, 5th ed., Wellesley-Cambridge Press — the column picture, elimination as $\mathbf{LU}$, and an engineering-minded treatment of rank and solvability.
- L. N. Trefethen and D. Bau, *Numerical Linear Algebra*, SIAM — conditioning, stability, pivoting and the backward-error viewpoint, in short self-contained lectures.
- T. A. Davis, *Direct Methods for Sparse Linear Systems*, SIAM — sparse storage, fill-reducing orderings and sparse $\mathbf{LU}$/Cholesky; with C. W. Ho and A. E. Ruehli's 1975 paper on modified nodal analysis as the classic circuit-side reference.
