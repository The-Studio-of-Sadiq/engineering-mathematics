# Appendix C — Linear Algebra

[← Table of Contents](../00-front/table-of-contents.md)

> This appendix gives complete proofs of the central theorems of Part III ([Ch. 8](../part-3-linear-algebra/08-matrices-and-linear-systems.md), [Ch. 9](../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md), [Ch. 10](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md)), which the chapters state with only sketches, together with the two results from numerical analysis and dynamics that lean on them most heavily: the condition-number error bound ([§8.10](../part-3-linear-algebra/08-matrices-and-linear-systems.md), [Ch. 22](../part-7-numerical-computational/22-numerical-methods.md)) and the eigenvalue test for stability of $\dot{\mathbf{x}}=\mathbf{A}\mathbf{x}$ ([§10.9](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md), [§12.7](../part-4-differential-equations/12-linear-differential-equations.md)).
>
> **Prerequisites.** Matrix algebra, row reduction and echelon form (Ch. 8); vector spaces, span, independence, basis and dimension (Ch. 9), in particular the facts that every independent set in a finite-dimensional space can be extended to a basis and that any $n$ independent vectors in an $n$-dimensional space form a basis. The fundamental theorem of algebra (every non-constant complex polynomial has a complex root) is used without proof.
>
> **Rigor.** Hypotheses are stated exactly and every step is justified. Scalars are $\mathbb{F}=\mathbb{R}$ or $\mathbb{C}$. Where a result holds for complex matrices with $\mathsf{T}$ replaced by $\mathsf{H}$, we say so; proofs are written in whichever setting is the more natural. The sections are ordered so that each proof uses only earlier ones.

---

## C.1 Determinants: the product rule and invertibility

Chapter 8 defines $\det\mathbf{A}$ by cofactor expansion and notes that it equals the Leibniz sum. Here we take the Leibniz sum as the working definition, because the proofs are cleanest from it.

> [!definition] Permutations and their sign
> A **permutation** of $\{1,\dots,n\}$ is a bijection $\sigma$ of that set to itself; the set of all $n!$ of them is $S_n$. A **transposition** swaps two elements and fixes the rest. Every permutation is a product of transpositions, and although the factorisation is not unique, the **parity** of the number of transpositions is; $\operatorname{sgn}\sigma=+1$ for an even number, $-1$ for an odd number. Hence $\operatorname{sgn}(\sigma\tau)=\operatorname{sgn}\sigma\,\operatorname{sgn}\tau$ and $\operatorname{sgn}\sigma^{-1}=\operatorname{sgn}\sigma$.

> [!definition] Determinant (Leibniz form)
> For $\mathbf{A}=[a_{ij}]\in\mathbb{F}^{n\times n}$,
> $$
> \det\mathbf{A}=\sum_{\sigma\in S_n}\operatorname{sgn}(\sigma)\,a_{1\sigma(1)}a_{2\sigma(2)}\cdots a_{n\sigma(n)} .
> $$

We regard $\det$ as a function $D(\mathbf{r}_1,\dots,\mathbf{r}_n)$ of the rows $\mathbf{r}_i$ of $\mathbf{A}$.

> [!definition] Multilinear, alternating
> A function $f$ of $n$ row vectors in $\mathbb{F}^n$ is **multilinear** if it is linear in each row when the others are held fixed, and **alternating** if $f=0$ whenever two rows are equal.

> [!theorem] Theorem C.1 — Product rule and invertibility
> Let $\mathbf{A},\mathbf{B}\in\mathbb{F}^{n\times n}$. Then
> 1. $\det(\mathbf{A}\mathbf{B})=\det\mathbf{A}\,\det\mathbf{B}$;
> 2. $\det\mathbf{A}^\mathsf{T}=\det\mathbf{A}$;
> 3. $\mathbf{A}$ is invertible (there is $\mathbf{X}$ with $\mathbf{A}\mathbf{X}=\mathbf{X}\mathbf{A}=\mathbf{I}$) if and only if $\det\mathbf{A}\neq0$, and then $\det\mathbf{A}^{-1}=1/\det\mathbf{A}$.

### Proof

**Step 1: $\det$ is multilinear, alternating, and $\det\mathbf{I}=1$.** Each Leibniz term contains exactly one factor $a_{i\sigma(i)}$ from row $i$, so each term — and hence the sum — is linear in row $i$. For $\mathbf{I}$, the only non-zero term is $\sigma=\mathrm{id}$, giving $1$. For the alternating property, suppose rows $p\neq q$ are equal and let $\tau$ be the transposition swapping $p$ and $q$. The map $\sigma\mapsto\sigma\tau$ pairs the permutations into pairs $\{\sigma,\sigma\tau\}$ (it is its own inverse and has no fixed points). In the pair, the terms are

$$
\operatorname{sgn}(\sigma)\prod_i a_{i\sigma(i)}\quad\text{and}\quad\operatorname{sgn}(\sigma\tau)\prod_i a_{i\sigma\tau(i)}=-\operatorname{sgn}(\sigma)\prod_i a_{i\sigma\tau(i)} .
$$

The two products differ only in the factors for $i=p,q$: the first has $a_{p\sigma(p)}a_{q\sigma(q)}$, the second $a_{p\sigma(q)}a_{q\sigma(p)}$. Since rows $p$ and $q$ are equal, $a_{p k}=a_{q k}$ for all $k$, so the products are equal and the pair cancels. Hence $\det\mathbf{A}=0$.

**Step 2: alternating implies antisymmetric.** If $f$ is multilinear and alternating, swapping two rows negates $f$. Indeed, put $\mathbf{u}+\mathbf{v}$ in both positions $p$ and $q$ (other rows fixed) and expand by multilinearity:

$$
0=f(\dots,\mathbf{u}+\mathbf{v},\dots,\mathbf{u}+\mathbf{v},\dots)=f(\dots,\mathbf{u},\dots,\mathbf{u},\dots)+f(\dots,\mathbf{u},\dots,\mathbf{v},\dots)+f(\dots,\mathbf{v},\dots,\mathbf{u},\dots)+f(\dots,\mathbf{v},\dots,\mathbf{v},\dots).
$$

The first and last terms vanish (equal rows), so $f(\dots,\mathbf{u},\dots,\mathbf{v},\dots)=-f(\dots,\mathbf{v},\dots,\mathbf{u},\dots)$.

**Step 3: uniqueness lemma.** *If $f$ is multilinear and alternating in the rows of an $n\times n$ matrix, then $f(\mathbf{A})=f(\mathbf{I})\det\mathbf{A}$ for every $\mathbf{A}$.* Write row $i$ as $\mathbf{r}_i=\sum_{k}a_{ik}\mathbf{e}_k^\mathsf{T}$ and expand $f$ by multilinearity in all $n$ rows:

$$
f(\mathbf{A})=\sum_{k_1=1}^n\cdots\sum_{k_n=1}^n a_{1k_1}\cdots a_{nk_n}\,f(\mathbf{e}_{k_1}^\mathsf{T},\dots,\mathbf{e}_{k_n}^\mathsf{T}).
$$

If two indices $k_p=k_q$ coincide, the corresponding $f$ value has two equal rows and vanishes. The surviving index tuples are exactly those with $k_i=\sigma(i)$ for a permutation $\sigma$. For such a tuple, the rows $\mathbf{e}_{\sigma(1)}^\mathsf{T},\dots,\mathbf{e}_{\sigma(n)}^\mathsf{T}$ can be sorted into $\mathbf{I}$ by a sequence of row swaps whose number has the parity of $\sigma$; by Step 2 each swap flips the sign, so $f(\mathbf{e}_{\sigma(1)}^\mathsf{T},\dots,\mathbf{e}_{\sigma(n)}^\mathsf{T})=\operatorname{sgn}(\sigma)f(\mathbf{I})$. Therefore

$$
f(\mathbf{A})=f(\mathbf{I})\sum_{\sigma\in S_n}\operatorname{sgn}(\sigma)a_{1\sigma(1)}\cdots a_{n\sigma(n)}=f(\mathbf{I})\det\mathbf{A}.
$$

**Step 4: product rule.** Fix $\mathbf{B}$ and define $f(\mathbf{A})=\det(\mathbf{A}\mathbf{B})$. Row $i$ of $\mathbf{A}\mathbf{B}$ is $\mathbf{r}_i\mathbf{B}$, a linear function of $\mathbf{r}_i$ alone. Hence $f$ is multilinear in the rows of $\mathbf{A}$ (composition of the linear map $\mathbf{r}_i\mapsto\mathbf{r}_i\mathbf{B}$ with the multilinear $\det$), and alternating (equal rows of $\mathbf{A}$ give equal rows of $\mathbf{A}\mathbf{B}$). By Step 3, $\det(\mathbf{A}\mathbf{B})=f(\mathbf{I})\det\mathbf{A}=\det\mathbf{B}\,\det\mathbf{A}$. This proves (1).

**Step 5: transpose.** $\det\mathbf{A}^\mathsf{T}=\sum_\sigma\operatorname{sgn}\sigma\prod_i a_{\sigma(i)i}$. Re-index each product by $k=\sigma(i)$, i.e. $i=\sigma^{-1}(k)$: $\prod_i a_{\sigma(i)i}=\prod_k a_{k\sigma^{-1}(k)}$. As $\sigma$ runs over $S_n$ so does $\pi=\sigma^{-1}$, with $\operatorname{sgn}\pi=\operatorname{sgn}\sigma$. So $\det\mathbf{A}^\mathsf{T}=\sum_\pi\operatorname{sgn}\pi\prod_k a_{k\pi(k)}=\det\mathbf{A}$. This proves (2). In particular $\det$ is also multilinear and alternating in the **columns**.

**Step 6: invertible $\Rightarrow$ $\det\neq0$.** If $\mathbf{A}\mathbf{X}=\mathbf{I}$ then by (1), $\det\mathbf{A}\det\mathbf{X}=\det\mathbf{I}=1$, so $\det\mathbf{A}\neq0$ and $\det\mathbf{A}^{-1}=1/\det\mathbf{A}$.

**Step 7: $\det\neq0$ $\Rightarrow$ invertible.** First, *if the columns $\mathbf{a}_1,\dots,\mathbf{a}_n$ of $\mathbf{A}$ are dependent, then $\det\mathbf{A}=0$.* Some column is a combination of the others, say $\mathbf{a}_k=\sum_{j\neq k}c_j\mathbf{a}_j$. By column multilinearity,

$$
\det\mathbf{A}=\sum_{j\neq k}c_j\det[\mathbf{a}_1,\dots,\underbrace{\mathbf{a}_j}_{\text{position }k},\dots,\mathbf{a}_n],
$$

and each determinant on the right has column $\mathbf{a}_j$ in two positions, so vanishes. Contrapositively, $\det\mathbf{A}\neq0$ implies the $n$ columns are independent, hence a basis of $\mathbb{F}^n$ (Ch. 9). So each $\mathbf{e}_i$ is a combination $\mathbf{A}\mathbf{x}_i$, and $\mathbf{X}=[\mathbf{x}_1,\dots,\mathbf{x}_n]$ satisfies $\mathbf{A}\mathbf{X}=\mathbf{I}$. By Step 6, $\det\mathbf{X}=1/\det\mathbf{A}\neq0$, so by the same argument $\mathbf{X}$ has a right inverse $\mathbf{Y}$, $\mathbf{X}\mathbf{Y}=\mathbf{I}$. Then

$$
\mathbf{A}=\mathbf{A}(\mathbf{X}\mathbf{Y})=(\mathbf{A}\mathbf{X})\mathbf{Y}=\mathbf{Y},
$$

so $\mathbf{X}\mathbf{A}=\mathbf{X}\mathbf{Y}=\mathbf{I}$ as well: $\mathbf{X}$ is a two-sided inverse. This proves (3). $\blacksquare$

### Remarks

- **What the proof really shows.** Step 3 says the determinant is the *only* multilinear alternating function with $f(\mathbf{I})=1$. That characterisation, not the $n!$-term formula, is what makes the determinant useful: it is why row replacement leaves $\det$ unchanged, row swaps negate it, and $\det\mathbf{A}=(-1)^s u_{11}\cdots u_{nn}$ from elimination ([§8.3](../part-3-linear-algebra/08-matrices-and-linear-systems.md)). Cofactor expansion along row $1$ is also multilinear, alternating and equals $1$ at $\mathbf{I}$, which is the cleanest proof that the two definitions in Ch. 8 agree.
- **Square matrices only.** For rectangular $\mathbf{A}$ ($m\times n$) and $\mathbf{B}$ ($n\times m$) with $m<n$, $\mathbf{A}\mathbf{B}$ can be invertible though neither factor has a determinant; the replacement is the Cauchy–Binet formula. For $m>n$, $\det(\mathbf{A}\mathbf{B})=0$ always, since $\operatorname{rank}(\mathbf{A}\mathbf{B})\le n<m$.
- **Sum rule fails.** $\det(\mathbf{A}+\mathbf{B})\neq\det\mathbf{A}+\det\mathbf{B}$ in general: $\mathbf{A}=\mathbf{B}=\mathbf{I}_2$ gives $4\neq2$.
- **Step 7 is the theoretical half of the invertible matrix theorem** ([§8.4](../part-3-linear-algebra/08-matrices-and-linear-systems.md)): it also shows that a one-sided inverse of a square matrix is automatically two-sided. Numerically, though, $\det\mathbf{A}$ is a poor measure of "nearly singular" — $\det(10^{-1}\mathbf{I}_{30})=10^{-30}$ for a perfectly conditioned matrix — which is why §C.11 uses the condition number instead.
- **Uses.** The characteristic polynomial $\det(\lambda\mathbf{I}-\mathbf{A})$ and the fact that $\lambda$ is an eigenvalue iff $\mathbf{A}-\lambda\mathbf{I}$ is singular ([§10.1](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md)) rest on part 3; similarity invariance $\det(\mathbf{S}^{-1}\mathbf{A}\mathbf{S})=\det\mathbf{A}$ and $\det\mathbf{A}=\prod\lambda_i$ rest on part 1; Cramer's rule and the Jacobian change of variables ([Ch. 6](../part-2-calculus/06-multivariable-calculus.md)) use both.


## C.2 The rank–nullity theorem

> [!definition] Kernel, range, rank, nullity
> For a linear map $T:V\to W$, the **kernel** (null space) is $\ker T=\{\mathbf{v}\in V:T\mathbf{v}=\mathbf{0}\}$ and the **range** is $\operatorname{range}T=\{T\mathbf{v}:\mathbf{v}\in V\}$. Both are subspaces. Their dimensions are the **nullity** and the **rank**. For a matrix $\mathbf{A}\in\mathbb{F}^{m\times n}$ acting by $\mathbf{x}\mapsto\mathbf{A}\mathbf{x}$, the kernel is $\operatorname{null}(\mathbf{A})$ and the range is the column space $\operatorname{col}(\mathbf{A})$.

> [!theorem] Theorem C.2 — Rank–nullity
> Let $V$ be finite-dimensional and $T:V\to W$ linear. Then
> $$
> \dim\operatorname{range}T+\dim\ker T=\dim V .
> $$
> For $\mathbf{A}\in\mathbb{F}^{m\times n}$: $\operatorname{rank}\mathbf{A}+\dim\operatorname{null}\mathbf{A}=n$. Moreover the row rank equals the column rank, $\dim\operatorname{row}\mathbf{A}=\dim\operatorname{col}\mathbf{A}$, and consequently $\dim\operatorname{null}\mathbf{A}^\mathsf{T}=m-\operatorname{rank}\mathbf{A}$.

### Proof

**Rank–nullity.** Let $n=\dim V$. The kernel is a subspace of $V$, hence finite-dimensional; choose a basis $\mathbf{k}_1,\dots,\mathbf{k}_q$ of it ($q=0$ if $\ker T=\{\mathbf{0}\}$). An independent set extends to a basis of $V$, so there are $\mathbf{v}_1,\dots,\mathbf{v}_s$ with $\mathbf{k}_1,\dots,\mathbf{k}_q,\mathbf{v}_1,\dots,\mathbf{v}_s$ a basis of $V$; thus $q+s=n$. We claim that $T\mathbf{v}_1,\dots,T\mathbf{v}_s$ is a basis of $\operatorname{range}T$, which gives $\dim\operatorname{range}T=s=n-q$.

*Spanning.* Any $\mathbf{w}\in\operatorname{range}T$ is $T\mathbf{v}$ for some $\mathbf{v}=\sum_i a_i\mathbf{k}_i+\sum_j b_j\mathbf{v}_j$. By linearity and $T\mathbf{k}_i=\mathbf{0}$, $\mathbf{w}=\sum_j b_jT\mathbf{v}_j$.

*Independence.* Suppose $\sum_j c_jT\mathbf{v}_j=\mathbf{0}$. Then $T\big(\sum_j c_j\mathbf{v}_j\big)=\mathbf{0}$, so $\sum_j c_j\mathbf{v}_j\in\ker T$ and equals some $\sum_i d_i\mathbf{k}_i$. The relation $\sum_j c_j\mathbf{v}_j-\sum_i d_i\mathbf{k}_i=\mathbf{0}$ among basis vectors of $V$ forces all $c_j=0$ (and all $d_i=0$).

**Row rank equals column rank.** Let $r=\dim\operatorname{col}\mathbf{A}$ and let $\mathbf{c}_1,\dots,\mathbf{c}_r\in\mathbb{F}^m$ be a basis of $\operatorname{col}\mathbf{A}$; put $\mathbf{C}=[\mathbf{c}_1,\dots,\mathbf{c}_r]$ ($m\times r$). Each column $\mathbf{a}_j$ of $\mathbf{A}$ is a combination $\mathbf{a}_j=\mathbf{C}\mathbf{f}_j$; collecting the coefficient vectors into the $r\times n$ matrix $\mathbf{F}=[\mathbf{f}_1,\dots,\mathbf{f}_n]$ gives the **rank factorisation**

$$
\mathbf{A}=\mathbf{C}\mathbf{F}.
$$

Read row-wise, row $i$ of $\mathbf{A}$ equals $\sum_{k=1}^r c_{ik}\,(\text{row }k\text{ of }\mathbf{F})$: every row of $\mathbf{A}$ lies in the span of the $r$ rows of $\mathbf{F}$. Hence $\dim\operatorname{row}\mathbf{A}\le r=\dim\operatorname{col}\mathbf{A}$. Applying this inequality to $\mathbf{A}^\mathsf{T}$, whose rows are the columns of $\mathbf{A}$, gives the reverse inequality. So the two are equal.

**Left null space.** Apply the matrix form of rank–nullity to $\mathbf{A}^\mathsf{T}\in\mathbb{F}^{n\times m}$, which has $m$ columns and rank $\dim\operatorname{col}\mathbf{A}^\mathsf{T}=\dim\operatorname{row}\mathbf{A}=r$: $\dim\operatorname{null}\mathbf{A}^\mathsf{T}=m-r$. $\blacksquare$

### Remarks

- **Finite dimension of $V$ is essential; of $W$ it is not.** On the infinite-dimensional space of polynomials, differentiation $D$ has a one-dimensional kernel (the constants) and is onto, so the identity would read $\infty=\infty+1$ and carries no information. Worse, a map with zero kernel need not be onto there: multiplication by $t$ is one-to-one but misses the constants. For square matrices, by contrast, rank–nullity gives "one-to-one $\iff$ onto", the backbone of the invertible matrix theorem ([§8.4](../part-3-linear-algebra/08-matrices-and-linear-systems.md)).
- **Nothing depends on row reduction.** The chapter proof ([§9.6](../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md)) counts pivot and free columns of the reduced echelon form; the proof above shows the result is a statement about the map, independent of any algorithm or basis.
- **Uses.** Existence and uniqueness for $\mathbf{A}\mathbf{x}=\mathbf{b}$ ([§8.6](../part-3-linear-algebra/08-matrices-and-linear-systems.md)): the solution set is $\mathbf{x}_p+\operatorname{null}\mathbf{A}$, an affine set of dimension $n-r$. In circuit theory, for a connected graph with $N$ nodes and $B$ branches the reduced incidence matrix has rank $N-1$, so the loop-current space has dimension $B-N+1$, the number of independent meshes ([§8.5](../part-3-linear-algebra/08-matrices-and-linear-systems.md)). Sections C.6 and C.10 need $\operatorname{rank}(\mathbf{A}^\mathsf{T}\mathbf{A})=\operatorname{rank}\mathbf{A}$, which follows from rank–nullity once the null spaces are shown equal.

## C.3 Existence and uniqueness of the LU factorisation

> [!definition] Leading principal submatrices and minors
> For $\mathbf{A}\in\mathbb{F}^{n\times n}$, $\mathbf{A}_k=\mathbf{A}_{1:k,1:k}$ is the **leading principal submatrix** of order $k$, and $\Delta_k=\det\mathbf{A}_k$ the **leading principal minor**; set $\Delta_0=1$. A **unit lower triangular** matrix is lower triangular with ones on the diagonal. A **permutation matrix** is $\mathbf{I}$ with its rows reordered.

Three facts about triangular matrices are used throughout. Products and inverses of (unit) lower triangular matrices are (unit) lower triangular, and likewise for upper. A triangular matrix is invertible iff its diagonal entries are non-zero, its determinant being their product ([§8.3](../part-3-linear-algebra/08-matrices-and-linear-systems.md)). And the **leading-block rule**: if $\mathbf{L}$ is lower and $\mathbf{U}$ upper triangular then

$$
(\mathbf{L}\mathbf{U})_k=\mathbf{L}_k\mathbf{U}_k ,
$$

because $(\mathbf{L}\mathbf{U})_{ij}=\sum_{p}\ell_{ip}u_{pj}$ and for $i,j\le k$ every non-zero term has $p\le\min(i,j)\le k$.

> [!theorem] Theorem C.3 — LU without pivoting
> Let $\mathbf{A}\in\mathbb{F}^{n\times n}$.
> 1. If $\Delta_1,\dots,\Delta_{n-1}$ are all non-zero, then $\mathbf{A}=\mathbf{L}\mathbf{U}$ with $\mathbf{L}$ unit lower triangular and $\mathbf{U}$ upper triangular, and this factorisation is **unique**. The pivots are $u_{kk}=\Delta_k/\Delta_{k-1}$.
> 2. Conversely, if $\mathbf{A}$ is nonsingular and has such a factorisation, then all $\Delta_k\neq0$, $k=1,\dots,n$.
>
> Thus a nonsingular $\mathbf{A}$ has an LU factorisation (Gaussian elimination runs without row interchanges) if and only if every leading principal minor is non-zero.

> [!theorem] Theorem C.4 — LU with pivoting
> For every nonsingular $\mathbf{A}\in\mathbb{F}^{n\times n}$ there is a permutation matrix $\mathbf{P}$ such that $\mathbf{P}\mathbf{A}=\mathbf{L}\mathbf{U}$, with $\mathbf{L}$ unit lower triangular and $\mathbf{U}$ upper triangular and nonsingular. $\mathbf{P}$ can be chosen (partial pivoting) so that $\lvert\ell_{ik}\rvert\le1$ for all $i,k$.

### Proof of Theorem C.3

**Existence, by induction on $n$.** For $n=1$, $[a]=[1][a]$. Suppose the result holds for order $n-1$ and let $\Delta_1,\dots,\Delta_{n-1}\neq0$. Partition

$$
\mathbf{A}=\begin{bmatrix}\mathbf{A}_{n-1}&\mathbf{b}\\ \mathbf{c}^\mathsf{T}&d\end{bmatrix}.
$$

The matrix $\mathbf{A}_{n-1}$ has leading minors $\Delta_1,\dots,\Delta_{n-2}$, all non-zero, so by hypothesis $\mathbf{A}_{n-1}=\mathbf{L}'\mathbf{U}'$. Since $\det\mathbf{U}'=\det\mathbf{L}'\det\mathbf{U}'=\Delta_{n-1}\neq0$, $\mathbf{U}'$ is invertible, and $\mathbf{L}'$ is invertible as a unit triangular matrix. Seek

$$
\mathbf{A}=\begin{bmatrix}\mathbf{L}'&\mathbf{0}\\ \boldsymbol{\ell}^\mathsf{T}&1\end{bmatrix}\begin{bmatrix}\mathbf{U}'&\mathbf{u}\\ \mathbf{0}^\mathsf{T}&u_{nn}\end{bmatrix}=\begin{bmatrix}\mathbf{L}'\mathbf{U}'&\mathbf{L}'\mathbf{u}\\ \boldsymbol{\ell}^\mathsf{T}\mathbf{U}'&\boldsymbol{\ell}^\mathsf{T}\mathbf{u}+u_{nn}\end{bmatrix}.
$$

Matching blocks requires

$$
\mathbf{L}'\mathbf{u}=\mathbf{b},\qquad \boldsymbol{\ell}^\mathsf{T}\mathbf{U}'=\mathbf{c}^\mathsf{T},\qquad u_{nn}=d-\boldsymbol{\ell}^\mathsf{T}\mathbf{u},
$$

which have the solutions $\mathbf{u}=\mathbf{L}'^{-1}\mathbf{b}$, $\boldsymbol{\ell}^\mathsf{T}=\mathbf{c}^\mathsf{T}\mathbf{U}'^{-1}$, $u_{nn}=d-\mathbf{c}^\mathsf{T}\mathbf{A}_{n-1}^{-1}\mathbf{b}$. The two factors are unit lower and upper triangular, completing the induction. Note that $\Delta_n$ is not needed: $u_{nn}$ may be zero.

**Uniqueness, by the same induction.** Suppose $\mathbf{A}=\mathbf{L}\mathbf{U}$ is any such factorisation. By the leading-block rule $\mathbf{A}_{n-1}=\mathbf{L}_{n-1}\mathbf{U}_{n-1}$, so by uniqueness at order $n-1$, $\mathbf{L}_{n-1}=\mathbf{L}'$ and $\mathbf{U}_{n-1}=\mathbf{U}'$. The remaining blocks must then satisfy the three displayed equations, each of which has exactly one solution because $\mathbf{L}'$ and $\mathbf{U}'$ are invertible. So $\mathbf{L},\mathbf{U}$ are determined.

**Pivots.** By the leading-block rule, $\Delta_k=\det\mathbf{L}_k\det\mathbf{U}_k=u_{11}u_{22}\cdots u_{kk}$ for every $k$. Dividing consecutive minors gives $u_{kk}=\Delta_k/\Delta_{k-1}$.

**Converse.** If $\mathbf{A}=\mathbf{L}\mathbf{U}$ is nonsingular then $0\neq\det\mathbf{A}=\prod_{i=1}^nu_{ii}$, so every $u_{ii}\neq0$ and $\Delta_k=\prod_{i\le k}u_{ii}\neq0$ for every $k$. $\blacksquare$

### Proof of Theorem C.4

Induction on $n$; $n=1$ is trivial with $\mathbf{P}=[1]$. Let $\mathbf{A}$ be nonsingular of order $n$. Its first column is non-zero (otherwise $\det\mathbf{A}=0$ by column linearity), so choose a row $p$ with $\lvert a_{p1}\rvert=\max_i\lvert a_{i1}\rvert>0$ and let $\mathbf{P}_1$ swap rows $1$ and $p$. Write

$$
\mathbf{P}_1\mathbf{A}=\begin{bmatrix}\alpha&\mathbf{w}^\mathsf{T}\\ \mathbf{v}&\mathbf{A}'\end{bmatrix}=\begin{bmatrix}1&\mathbf{0}^\mathsf{T}\\ \mathbf{v}/\alpha&\mathbf{I}\end{bmatrix}\begin{bmatrix}\alpha&\mathbf{w}^\mathsf{T}\\ \mathbf{0}&\mathbf{S}\end{bmatrix},\qquad \mathbf{S}=\mathbf{A}'-\frac{\mathbf{v}\mathbf{w}^\mathsf{T}}{\alpha},
$$

which is verified by multiplying out. The entries of $\mathbf{v}/\alpha$ have modulus at most $1$ by the choice of pivot. Taking determinants, $\pm\det\mathbf{A}=\alpha\det\mathbf{S}$, so the **Schur complement** $\mathbf{S}$, of order $n-1$, is nonsingular. By the induction hypothesis $\mathbf{P}_2\mathbf{S}=\mathbf{L}_2\mathbf{U}_2$ with multipliers of modulus at most $1$. Put

$$
\mathbf{P}=\begin{bmatrix}1&\mathbf{0}^\mathsf{T}\\ \mathbf{0}&\mathbf{P}_2\end{bmatrix}\mathbf{P}_1,\qquad
\mathbf{L}=\begin{bmatrix}1&\mathbf{0}^\mathsf{T}\\ \mathbf{P}_2\mathbf{v}/\alpha&\mathbf{L}_2\end{bmatrix},\qquad
\mathbf{U}=\begin{bmatrix}\alpha&\mathbf{w}^\mathsf{T}\\ \mathbf{0}&\mathbf{U}_2\end{bmatrix}.
$$

Then

$$
\mathbf{P}\mathbf{A}=\begin{bmatrix}1&\mathbf{0}^\mathsf{T}\\ \mathbf{0}&\mathbf{P}_2\end{bmatrix}\begin{bmatrix}\alpha&\mathbf{w}^\mathsf{T}\\ \mathbf{v}&\mathbf{A}'\end{bmatrix}=\begin{bmatrix}\alpha&\mathbf{w}^\mathsf{T}\\ \mathbf{P}_2\mathbf{v}&\mathbf{P}_2\mathbf{A}'\end{bmatrix},
\qquad
\mathbf{L}\mathbf{U}=\begin{bmatrix}\alpha&\mathbf{w}^\mathsf{T}\\ \mathbf{P}_2\mathbf{v}&\mathbf{P}_2\mathbf{v}\mathbf{w}^\mathsf{T}/\alpha+\mathbf{L}_2\mathbf{U}_2\end{bmatrix},
$$

and the bottom-right blocks agree because $\mathbf{L}_2\mathbf{U}_2=\mathbf{P}_2\mathbf{S}=\mathbf{P}_2\mathbf{A}'-\mathbf{P}_2\mathbf{v}\mathbf{w}^\mathsf{T}/\alpha$. $\mathbf{P}$ is a product of permutation matrices, hence a permutation matrix; $\mathbf{L}$ is unit lower triangular with all multipliers of modulus at most $1$ (the permutation $\mathbf{P}_2$ only reorders the entries of $\mathbf{v}/\alpha$); and $\mathbf{U}$ is upper triangular with $\det\mathbf{U}=\alpha\det\mathbf{U}_2\neq0$. $\blacksquare$

### Remarks

- **The hypothesis cannot be dropped.** $\mathbf{A}=\begin{bmatrix}0&1\\1&0\end{bmatrix}$ is nonsingular but $\Delta_1=0$, and it has no LU factorisation: $u_{11}=a_{11}=0$ and $\ell_{21}u_{11}=a_{21}=1$ are incompatible. One row swap fixes it, as Theorem C.4 guarantees.
- **Uniqueness needs $\Delta_1,\dots,\Delta_{n-1}\neq0$.** $\begin{bmatrix}0&0\\0&0\end{bmatrix}=\begin{bmatrix}1&0\\ \ell&1\end{bmatrix}\begin{bmatrix}0&0\\0&0\end{bmatrix}$ for every $\ell$. Without the unit-diagonal normalisation of $\mathbf{L}$ uniqueness also fails, since $(\mathbf{L}\mathbf{D}^{-1})(\mathbf{D}\mathbf{U})$ is another factorisation for any nonsingular diagonal $\mathbf{D}$.
- **Existence is not stability.** $\mathbf{A}=\begin{bmatrix}\varepsilon&1\\1&1\end{bmatrix}$ with $0<\lvert\varepsilon\rvert\ll1$ has $\Delta_1=\varepsilon\neq0$, so LU without pivoting exists, but $\ell_{21}=1/\varepsilon$ and $u_{22}=1-1/\varepsilon$ are huge and rounding destroys the answer (see the pivoting discussion in [§8.7](../part-3-linear-algebra/08-matrices-and-linear-systems.md)). The bound $\lvert\ell_{ik}\rvert\le1$ in Theorem C.4 keeps the computed factors representative of $\mathbf{A}$; growth in $\mathbf{U}$ can still occur in contrived cases but is rare in practice ([§24.3](../part-7-numerical-computational/24-computational-mathematics.md)).
- **A check.** For $\mathbf{A}=\begin{bmatrix}2&1&1\\4&3&3\\8&7&9\end{bmatrix}$: $\Delta_1=2$, $\Delta_2=2$, $\Delta_3=4$, so the pivots without pivoting are $2,\ 1,\ 2$ and $\det\mathbf{A}=4$. Partial pivoting instead brings row 3 to the top and produces pivots $8,\ -0.75,\ -0.\overline{6}$ (product $4$, since the row permutation used is a cyclic shift with determinant $+1$) and multipliers $0.25,\ 0.5,\ 0.\overline{6}$, all at most $1$.
- **Uses.** Solving nodal equations ([§8.9](../part-3-linear-algebra/08-matrices-and-linear-systems.md)); the pivot formula $u_{kk}=\Delta_k/\Delta_{k-1}$ is the bridge to Sylvester's criterion in §C.9; the Schur complement $\mathbf{S}$ is exactly Kron reduction of a network ([§8.4](../part-3-linear-algebra/08-matrices-and-linear-systems.md)).

## C.4 Cauchy–Schwarz and triangle inequalities

We use the book's convention ([§9.8](../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md)): an inner product on a vector space $V$ over $\mathbb{F}$ is linear in its first argument, conjugate-symmetric, $\langle\mathbf{u},\mathbf{v}\rangle=\langle\mathbf{v},\mathbf{u}\rangle^*$, and positive, $\langle\mathbf{u},\mathbf{u}\rangle>0$ for $\mathbf{u}\neq\mathbf{0}$. It follows that $\langle\mathbf{u},a\mathbf{v}\rangle=a^*\langle\mathbf{u},\mathbf{v}\rangle$. The induced norm is $\lVert\mathbf{u}\rVert=\sqrt{\langle\mathbf{u},\mathbf{u}\rangle}$. No finite-dimensionality is assumed in this section: the results hold verbatim for signals with $\langle x,y\rangle=\int x(t)y^*(t)\,dt$.

> [!theorem] Theorem C.5 — Cauchy–Schwarz inequality
> For all $\mathbf{u},\mathbf{v}$ in an inner-product space,
> $$
> \lvert\langle\mathbf{u},\mathbf{v}\rangle\rvert\le\lVert\mathbf{u}\rVert\,\lVert\mathbf{v}\rVert ,
> $$
> with equality if and only if $\mathbf{u}$ and $\mathbf{v}$ are linearly dependent (one is a scalar multiple of the other).

> [!theorem] Theorem C.6 — Triangle inequality
> For all $\mathbf{u},\mathbf{v}$, $\ \lVert\mathbf{u}+\mathbf{v}\rVert\le\lVert\mathbf{u}\rVert+\lVert\mathbf{v}\rVert$, with equality if and only if one vector is a non-negative real multiple of the other. Consequently the induced norm is a norm, and $d(\mathbf{u},\mathbf{v})=\lVert\mathbf{u}-\mathbf{v}\rVert$ is a distance.

### Proof of Theorem C.5

If $\mathbf{v}=\mathbf{0}$ both sides are $0$ (since $\langle\mathbf{u},\mathbf{0}\rangle=\langle\mathbf{u},0\cdot\mathbf{0}\rangle=0$) and the vectors are dependent, so the claim holds. Assume $\mathbf{v}\neq\mathbf{0}$. For any scalar $s$, positivity gives

$$
0\le\lVert\mathbf{u}-s\mathbf{v}\rVert^2=\langle\mathbf{u}-s\mathbf{v},\mathbf{u}-s\mathbf{v}\rangle=\lVert\mathbf{u}\rVert^2-s^*\langle\mathbf{u},\mathbf{v}\rangle-s\langle\mathbf{v},\mathbf{u}\rangle+\lvert s\rvert^2\lVert\mathbf{v}\rVert^2 ,
$$

where we expanded using linearity in the first slot and conjugate-linearity in the second. Choose $s=\langle\mathbf{u},\mathbf{v}\rangle/\lVert\mathbf{v}\rVert^2$, the coefficient of the projection of $\mathbf{u}$ onto $\mathbf{v}$. Writing $c=\langle\mathbf{u},\mathbf{v}\rangle$, so that $\langle\mathbf{v},\mathbf{u}\rangle=c^*$, the three last terms are

$$
s^*c=\frac{c^*c}{\lVert\mathbf{v}\rVert^2},\qquad s\,c^*=\frac{c\,c^*}{\lVert\mathbf{v}\rVert^2},\qquad \lvert s\rvert^2\lVert\mathbf{v}\rVert^2=\frac{\lvert c\rvert^2}{\lVert\mathbf{v}\rVert^2},
$$

so

$$
0\le\lVert\mathbf{u}-s\mathbf{v}\rVert^2=\lVert\mathbf{u}\rVert^2-\frac{\lvert\langle\mathbf{u},\mathbf{v}\rangle\rvert^2}{\lVert\mathbf{v}\rVert^2}.
$$

Multiplying by $\lVert\mathbf{v}\rVert^2>0$ and taking square roots gives the inequality.

*Equality.* With $\mathbf{v}\neq\mathbf{0}$, equality holds iff $\lVert\mathbf{u}-s\mathbf{v}\rVert=0$, i.e. (by positivity) $\mathbf{u}=s\mathbf{v}$. Conversely, if $\mathbf{u}=a\mathbf{v}$ then $\lvert\langle\mathbf{u},\mathbf{v}\rangle\rvert=\lvert a\rvert\lVert\mathbf{v}\rVert^2=\lVert\mathbf{u}\rVert\lVert\mathbf{v}\rVert$. $\blacksquare$

### Proof of Theorem C.6

Expand:

$$
\lVert\mathbf{u}+\mathbf{v}\rVert^2=\lVert\mathbf{u}\rVert^2+\langle\mathbf{u},\mathbf{v}\rangle+\langle\mathbf{v},\mathbf{u}\rangle+\lVert\mathbf{v}\rVert^2=\lVert\mathbf{u}\rVert^2+2\operatorname{Re}\{\langle\mathbf{u},\mathbf{v}\rangle\}+\lVert\mathbf{v}\rVert^2 ,
$$

since $z+z^*=2\operatorname{Re}\{z\}$. Now $\operatorname{Re}\{z\}\le\lvert z\rvert$ and Cauchy–Schwarz give

$$
\lVert\mathbf{u}+\mathbf{v}\rVert^2\le\lVert\mathbf{u}\rVert^2+2\lVert\mathbf{u}\rVert\lVert\mathbf{v}\rVert+\lVert\mathbf{v}\rVert^2=\big(\lVert\mathbf{u}\rVert+\lVert\mathbf{v}\rVert\big)^2 .
$$

Take square roots (both sides are non-negative).

*Equality* requires equality in both steps: $\operatorname{Re}\{\langle\mathbf{u},\mathbf{v}\rangle\}=\lvert\langle\mathbf{u},\mathbf{v}\rangle\rvert$, i.e. $\langle\mathbf{u},\mathbf{v}\rangle$ is real and $\ge0$, and dependence. If $\mathbf{v}\neq\mathbf{0}$ then $\mathbf{u}=a\mathbf{v}$ and $\langle\mathbf{u},\mathbf{v}\rangle=a\lVert\mathbf{v}\rVert^2\ge0$ forces $a\ge0$; the case $\mathbf{v}=\mathbf{0}$ is trivial. Conversely $\mathbf{u}=a\mathbf{v}$, $a\ge0$, gives $\lVert(1+a)\mathbf{v}\rVert=(1+a)\lVert\mathbf{v}\rVert$.

*Norm axioms.* $\lVert\mathbf{u}\rVert\ge0$ with equality only at $\mathbf{0}$ is positivity; $\lVert a\mathbf{u}\rVert=\sqrt{aa^*\langle\mathbf{u},\mathbf{u}\rangle}=\lvert a\rvert\lVert\mathbf{u}\rVert$; the triangle inequality was just proved. For the distance, symmetry is $\lVert-\mathbf{w}\rVert=\lVert\mathbf{w}\rVert$ and the triangle inequality for $d$ is Theorem C.6 applied to $(\mathbf{u}-\mathbf{w})+(\mathbf{w}-\mathbf{v})$. $\blacksquare$

The proof is easiest to see in the real case, where $q(t)=\lVert\mathbf{u}-t\mathbf{v}\rVert^2=\lVert\mathbf{u}\rVert^2-2t\langle\mathbf{u},\mathbf{v}\rangle+t^2\lVert\mathbf{v}\rVert^2$ is a parabola that never goes below zero, so its discriminant $4\langle\mathbf{u},\mathbf{v}\rangle^2-4\lVert\mathbf{u}\rVert^2\lVert\mathbf{v}\rVert^2$ is $\le0$. For $\mathbf{u}=(3,1)$, $\mathbf{v}=(1,2)$: $\lVert\mathbf{u}\rVert^2=10$, $\langle\mathbf{u},\mathbf{v}\rangle=5$, $\lVert\mathbf{v}\rVert^2=5$, so $q(t)=10-10t+5t^2$ with minimum $q(1)=5=10-5^2/5$.

```desmos-graph
left=-1; right=3; bottom=-1; top=12
---
y=10-10x+5x^2|blue
y=5|red|dashed
(1,5)|label:minimum q(1) = 5
y=0|black
```
*The real Cauchy–Schwarz parabola $q(t)=\lVert\mathbf{u}-t\mathbf{v}\rVert^2$ for $\mathbf{u}=(3,1)$, $\mathbf{v}=(1,2)$: it stays above the axis, its minimum value is the squared distance from $\mathbf{u}$ to the line through $\mathbf{v}$, and it touches zero only when $\mathbf{u}$ is a multiple of $\mathbf{v}$.*

### Remarks

- **Positivity is the whole proof.** For an indefinite "inner product" such as the Minkowski form $\langle\mathbf{u},\mathbf{v}\rangle=u_1v_1-u_2v_2$ on $\mathbb{R}^2$, Cauchy–Schwarz fails: $\mathbf{u}=(1,0)$, $\mathbf{v}=(1,1)$ give $\lvert\langle\mathbf{u},\mathbf{v}\rangle\rvert=1$ but $\langle\mathbf{v},\mathbf{v}\rangle=0$. Positive *semi*-definiteness is enough for the inequality (the proof goes through after a limiting argument) but not for the equality case.
- **The conjugate is essential in the complex case.** Without it, $\langle(1,j),(1,j)\rangle=1+j^2=0$ and positivity fails ([§9.8](../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md)).
- **Uses.** Cauchy–Schwarz defines angles, $\cos\theta=\langle\mathbf{u},\mathbf{v}\rangle/(\lVert\mathbf{u}\rVert\lVert\mathbf{v}\rVert)\in[-1,1]$ ([§9.7](../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md)); it proves the **matched filter** is optimal, because the output SNR is proportional to $\lvert\langle s,h\rangle\rvert^2/\lVert h\rVert^2\le\lVert s\rVert^2$ with equality iff $h\propto s$ ([§9.8](../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md), [Ch. 21](../part-6-probability-statistics/21-random-processes.md)); and it gives $\lvert\rho_{XY}\rvert\le1$ for the correlation coefficient, with $\langle X,Y\rangle=E[XY]$ on zero-mean random variables ([Ch. 19](../part-6-probability-statistics/19-probability.md)). The triangle inequality is what makes $\lVert\cdot\rVert$ usable as an error measure everywhere in Part VII.

## C.5 Gram–Schmidt orthonormalisation and the QR factorisation

> [!theorem] Theorem C.7 — Gram–Schmidt
> Let $\mathbf{a}_1,\dots,\mathbf{a}_n$ be linearly independent vectors in an inner-product space. Define recursively, for $k=1,\dots,n$,
> $$
> \mathbf{w}_k=\mathbf{a}_k-\sum_{i=1}^{k-1}\langle\mathbf{a}_k,\mathbf{q}_i\rangle\,\mathbf{q}_i,\qquad \mathbf{q}_k=\frac{\mathbf{w}_k}{\lVert\mathbf{w}_k\rVert}.
> $$
> Then every $\mathbf{w}_k\neq\mathbf{0}$ (so the division is legitimate), the vectors $\mathbf{q}_1,\dots,\mathbf{q}_n$ are orthonormal, and for every $k$
> $$
> \operatorname{span}\{\mathbf{q}_1,\dots,\mathbf{q}_k\}=\operatorname{span}\{\mathbf{a}_1,\dots,\mathbf{a}_k\}.
> $$
> In particular every finite-dimensional inner-product space has an orthonormal basis, and every orthonormal set in it extends to an orthonormal basis.

> [!theorem] Theorem C.8 — QR factorisation
> Let $\mathbf{A}\in\mathbb{R}^{m\times n}$ have linearly independent columns ($m\ge n$, $\operatorname{rank}\mathbf{A}=n$). Then $\mathbf{A}=\mathbf{Q}\mathbf{R}$ where $\mathbf{Q}\in\mathbb{R}^{m\times n}$ has orthonormal columns ($\mathbf{Q}^\mathsf{T}\mathbf{Q}=\mathbf{I}_n$) and $\mathbf{R}\in\mathbb{R}^{n\times n}$ is upper triangular with positive diagonal. This factorisation is unique. There is also a **full** QR factorisation $\mathbf{A}=\widetilde{\mathbf{Q}}\begin{bmatrix}\mathbf{R}\\ \mathbf{0}\end{bmatrix}$ with $\widetilde{\mathbf{Q}}$ $m\times m$ orthogonal. (Complex case: $\mathsf{T}\to\mathsf{H}$, orthogonal $\to$ unitary.)

### Proof of Theorem C.7

Induction on $k$, with the statement $S_k$: "$\mathbf{w}_1,\dots,\mathbf{w}_k\neq\mathbf{0}$, $\mathbf{q}_1,\dots,\mathbf{q}_k$ are orthonormal, and $\operatorname{span}\{\mathbf{q}_1,\dots,\mathbf{q}_j\}=\operatorname{span}\{\mathbf{a}_1,\dots,\mathbf{a}_j\}$ for $j\le k$".

*Base.* $\mathbf{w}_1=\mathbf{a}_1\neq\mathbf{0}$, since an independent set contains no zero vector. $\mathbf{q}_1$ is a unit vector and spans the same line as $\mathbf{a}_1$.

*Step.* Assume $S_{k-1}$.

1. **$\mathbf{w}_k\neq\mathbf{0}$.** If $\mathbf{w}_k=\mathbf{0}$, then $\mathbf{a}_k=\sum_{i<k}\langle\mathbf{a}_k,\mathbf{q}_i\rangle\mathbf{q}_i\in\operatorname{span}\{\mathbf{q}_1,\dots,\mathbf{q}_{k-1}\}=\operatorname{span}\{\mathbf{a}_1,\dots,\mathbf{a}_{k-1}\}$, contradicting independence of the $\mathbf{a}$'s.
2. **Orthogonality.** For $l<k$, using linearity in the first slot and $\langle\mathbf{q}_i,\mathbf{q}_l\rangle=\delta_{il}$ for $i,l<k$,

$$
\langle\mathbf{w}_k,\mathbf{q}_l\rangle=\langle\mathbf{a}_k,\mathbf{q}_l\rangle-\sum_{i<k}\langle\mathbf{a}_k,\mathbf{q}_i\rangle\langle\mathbf{q}_i,\mathbf{q}_l\rangle=\langle\mathbf{a}_k,\mathbf{q}_l\rangle-\langle\mathbf{a}_k,\mathbf{q}_l\rangle=0 .
$$

Dividing by the positive real $\lVert\mathbf{w}_k\rVert$ keeps this, and makes $\lVert\mathbf{q}_k\rVert=1$.
3. **Spans.** By the defining formula $\mathbf{q}_k$ is a combination of $\mathbf{a}_k$ and $\mathbf{q}_1,\dots,\mathbf{q}_{k-1}$, which lie in $\operatorname{span}\{\mathbf{a}_1,\dots,\mathbf{a}_{k-1}\}$; so $\operatorname{span}\{\mathbf{q}_1,\dots,\mathbf{q}_k\}\subseteq\operatorname{span}\{\mathbf{a}_1,\dots,\mathbf{a}_k\}$. Conversely $\mathbf{a}_k=\lVert\mathbf{w}_k\rVert\mathbf{q}_k+\sum_{i<k}\langle\mathbf{a}_k,\mathbf{q}_i\rangle\mathbf{q}_i$ lies in the left-hand span. Both spans are therefore equal.

This proves $S_n$. For the final assertions: apply the process to any basis of the space; and to extend an orthonormal set $\mathbf{q}_1,\dots,\mathbf{q}_p$, extend it to a basis $\mathbf{q}_1,\dots,\mathbf{q}_p,\mathbf{a}_{p+1},\dots,\mathbf{a}_N$ and run the process — it leaves the first $p$ vectors unchanged, because for orthonormal input $\mathbf{w}_k=\mathbf{q}_k$. $\blacksquare$

### Proof of Theorem C.8

*Existence.* Run Gram–Schmidt on the columns $\mathbf{a}_1,\dots,\mathbf{a}_n$ with the standard inner product. Rearranging the definition of $\mathbf{w}_k$,

$$
\mathbf{a}_k=\sum_{i=1}^{k-1}r_{ik}\mathbf{q}_i+r_{kk}\mathbf{q}_k,\qquad r_{ik}=\langle\mathbf{a}_k,\mathbf{q}_i\rangle=\mathbf{q}_i^\mathsf{T}\mathbf{a}_k\ (i<k),\qquad r_{kk}=\lVert\mathbf{w}_k\rVert>0 .
$$

Setting $r_{ik}=0$ for $i>k$, this says column $k$ of $\mathbf{A}$ is column $k$ of $\mathbf{Q}\mathbf{R}$, with $\mathbf{Q}=[\mathbf{q}_1,\dots,\mathbf{q}_n]$ and $\mathbf{R}=[r_{ik}]$ upper triangular with positive diagonal. Orthonormality of the $\mathbf{q}_i$ is $\mathbf{Q}^\mathsf{T}\mathbf{Q}=\mathbf{I}$.

*Uniqueness.* Suppose $\mathbf{Q}_1\mathbf{R}_1=\mathbf{Q}_2\mathbf{R}_2$ are two such factorisations. Both $\mathbf{R}_i$ are invertible (positive diagonal). Then $\mathbf{Q}_1=\mathbf{Q}_2\mathbf{M}$ with $\mathbf{M}=\mathbf{R}_2\mathbf{R}_1^{-1}$ upper triangular with positive diagonal. Also

$$
\mathbf{I}=\mathbf{Q}_1^\mathsf{T}\mathbf{Q}_1=\mathbf{M}^\mathsf{T}\mathbf{Q}_2^\mathsf{T}\mathbf{Q}_2\mathbf{M}=\mathbf{M}^\mathsf{T}\mathbf{M},
$$

so $\mathbf{M}^{-1}=\mathbf{M}^\mathsf{T}$. The left side is upper triangular and the right side lower triangular, so $\mathbf{M}$ is diagonal; then $\mathbf{M}^2=\mathbf{I}$ with positive diagonal forces $\mathbf{M}=\mathbf{I}$. Hence $\mathbf{R}_1=\mathbf{R}_2$ and $\mathbf{Q}_1=\mathbf{Q}_2$.

*Full QR.* Extend $\mathbf{q}_1,\dots,\mathbf{q}_n$ to an orthonormal basis $\mathbf{q}_1,\dots,\mathbf{q}_m$ of $\mathbb{R}^m$ (Theorem C.7) and put $\widetilde{\mathbf{Q}}=[\mathbf{q}_1,\dots,\mathbf{q}_m]$; the extra columns multiply the zero block. $\blacksquare$

**Check.** For $\mathbf{A}=\begin{bmatrix}1&1\\1&0\\0&1\end{bmatrix}$: $\mathbf{q}_1=(1,1,0)/\sqrt2$, $r_{11}=\sqrt2$; $r_{12}=\langle\mathbf{a}_2,\mathbf{q}_1\rangle=1/\sqrt2$; $\mathbf{w}_2=(1,0,1)-\tfrac12(1,1,0)=(\tfrac12,-\tfrac12,1)$, $r_{22}=\sqrt{3/2}=1.2247$. A library QR returns the same $\lvert r_{ik}\rvert$, with signs flipped where it does not enforce a positive diagonal.

### Remarks

- **Independence is necessary.** If $\mathbf{a}_k\in\operatorname{span}\{\mathbf{a}_1,\dots,\mathbf{a}_{k-1}\}$ the process produces $\mathbf{w}_k=\mathbf{0}$ and stops. A rank-deficient $\mathbf{A}$ still has a factorisation $\mathbf{A}=\mathbf{Q}\mathbf{R}$ (skip such columns, or use column pivoting, $\mathbf{A}\boldsymbol{\Pi}=\mathbf{Q}\mathbf{R}$), but $\mathbf{R}$ then has zero diagonal entries and uniqueness is lost.
- **Exact versus floating-point.** The theorem is about exact arithmetic. In floating point, classical Gram–Schmidt loses orthogonality roughly in proportion to $\kappa(\mathbf{A})^2$ times the unit round-off; the modified version (subtract each projection from the running vector) and Householder reflections behave far better ([§9.10](../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md), [§24.3](../part-7-numerical-computational/24-computational-mathematics.md)). All three produce the same $\mathbf{Q}$ and $\mathbf{R}$ in exact arithmetic, by the uniqueness just proved.
- **Uses.** Orthogonal polynomials (Legendre, Chebyshev) are Gram–Schmidt applied to $1,t,t^2,\dots$ under an integral inner product ([§9.10](../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md), [§22.6](../part-7-numerical-computational/22-numerical-methods.md)). QR is the stable route to least squares (§C.6) and the engine of the QR eigenvalue algorithm. The extension statement is used in the induction step of the spectral theorem (§C.8) and in building $\mathbf{U}$ in the SVD (§C.10).

## C.6 The orthogonal projection theorem and least squares

> [!definition] Orthogonal complement and projection
> For a subspace $W$ of an inner-product space $V$, $W^\perp=\{\mathbf{z}\in V:\langle\mathbf{z},\mathbf{w}\rangle=0\ \text{for all }\mathbf{w}\in W\}$. A vector $\mathbf{p}\in W$ with $\mathbf{b}-\mathbf{p}\in W^\perp$ is called an **orthogonal projection** of $\mathbf{b}$ onto $W$.

> [!theorem] Theorem C.9 — Orthogonal projection theorem
> Let $W$ be a **finite-dimensional** subspace of an inner-product space $V$, and $\mathbf{b}\in V$. Then:
> 1. $\mathbf{b}$ has exactly one orthogonal projection $\mathbf{p}$ onto $W$; if $\mathbf{q}_1,\dots,\mathbf{q}_k$ is an orthonormal basis of $W$,
> $$
> \mathbf{p}=\sum_{i=1}^k\langle\mathbf{b},\mathbf{q}_i\rangle\,\mathbf{q}_i ;
> $$
> 2. $\mathbf{p}$ is the unique closest point of $W$ to $\mathbf{b}$: for every $\mathbf{w}\in W$,
> $$
> \lVert\mathbf{b}-\mathbf{w}\rVert^2=\lVert\mathbf{b}-\mathbf{p}\rVert^2+\lVert\mathbf{p}-\mathbf{w}\rVert^2 ;
> $$
> 3. (Bessel) $\lVert\mathbf{p}\rVert^2=\sum_i\lvert\langle\mathbf{b},\mathbf{q}_i\rangle\rvert^2\le\lVert\mathbf{b}\rVert^2$.

> [!theorem] Theorem C.10 — Normal equations
> Let $\mathbf{A}\in\mathbb{R}^{m\times n}$ and $\mathbf{b}\in\mathbb{R}^m$. A vector $\hat{\mathbf{x}}$ minimises $J(\mathbf{x})=\lVert\mathbf{b}-\mathbf{A}\mathbf{x}\rVert^2$ if and only if
> $$
> \mathbf{A}^\mathsf{T}\mathbf{A}\,\hat{\mathbf{x}}=\mathbf{A}^\mathsf{T}\mathbf{b}.
> $$
> The normal equations are always consistent, $\mathbf{A}\hat{\mathbf{x}}$ is the same vector (the projection of $\mathbf{b}$ onto $\operatorname{col}\mathbf{A}$) for every minimiser, and the minimiser is **unique if and only if $\operatorname{rank}\mathbf{A}=n$** (independent columns), in which case $\hat{\mathbf{x}}=(\mathbf{A}^\mathsf{T}\mathbf{A})^{-1}\mathbf{A}^\mathsf{T}\mathbf{b}=\mathbf{R}^{-1}\mathbf{Q}^\mathsf{T}\mathbf{b}$ with $\mathbf{A}=\mathbf{Q}\mathbf{R}$. (Complex case: $\mathsf{T}\to\mathsf{H}$.)

### Proof of Theorem C.9

*Existence.* $W$ has an orthonormal basis $\mathbf{q}_1,\dots,\mathbf{q}_k$ by Theorem C.7. Define $\mathbf{p}$ by the formula. It lies in $W$, and for each $l$

$$
\langle\mathbf{b}-\mathbf{p},\mathbf{q}_l\rangle=\langle\mathbf{b},\mathbf{q}_l\rangle-\sum_i\langle\mathbf{b},\mathbf{q}_i\rangle\langle\mathbf{q}_i,\mathbf{q}_l\rangle=0 .
$$

An arbitrary $\mathbf{w}\in W$ is $\sum_lc_l\mathbf{q}_l$, and $\langle\mathbf{b}-\mathbf{p},\mathbf{w}\rangle=\sum_lc_l^*\langle\mathbf{b}-\mathbf{p},\mathbf{q}_l\rangle=0$. So $\mathbf{b}-\mathbf{p}\in W^\perp$.

*Uniqueness.* If $\mathbf{p},\mathbf{p}'$ are both orthogonal projections, then $\mathbf{p}-\mathbf{p}'=(\mathbf{b}-\mathbf{p}')-(\mathbf{b}-\mathbf{p})$ lies in $W$ (difference of two elements of $W$) and in $W^\perp$ (difference of two elements of the subspace $W^\perp$). A vector orthogonal to itself has $\lVert\mathbf{p}-\mathbf{p}'\rVert^2=0$, so $\mathbf{p}=\mathbf{p}'$.

*Best approximation.* For $\mathbf{w}\in W$, write $\mathbf{b}-\mathbf{w}=(\mathbf{b}-\mathbf{p})+(\mathbf{p}-\mathbf{w})$. The first summand is in $W^\perp$ and the second in $W$, so their inner products vanish and the expansion of $\lVert\mathbf{b}-\mathbf{w}\rVert^2$ (as in the proof of C.6) reduces to Pythagoras:

$$
\lVert\mathbf{b}-\mathbf{w}\rVert^2=\lVert\mathbf{b}-\mathbf{p}\rVert^2+\lVert\mathbf{p}-\mathbf{w}\rVert^2\ \ge\ \lVert\mathbf{b}-\mathbf{p}\rVert^2 ,
$$

with equality iff $\lVert\mathbf{p}-\mathbf{w}\rVert=0$, i.e. $\mathbf{w}=\mathbf{p}$.

*Bessel.* Pythagoras for orthonormal vectors gives $\lVert\mathbf{p}\rVert^2=\sum_i\lvert\langle\mathbf{b},\mathbf{q}_i\rangle\rvert^2$, and the identity above with $\mathbf{w}=\mathbf{0}$ gives $\lVert\mathbf{b}\rVert^2=\lVert\mathbf{b}-\mathbf{p}\rVert^2+\lVert\mathbf{p}\rVert^2\ge\lVert\mathbf{p}\rVert^2$. $\blacksquare$

```desmos-graph
left=-1.5; right=4; bottom=-0.5; top=4
---
y=0.5x|blue
(1,3)|label:b
(2,1)|label:p = proj b
y=3-2(x-1)|red|dashed|1<x<2
(x-1)^2+(y-3)^2=5|grey|dotted
(0.4,0.2)|label:w (any other point of W)
y=3+\frac{-2.8}{-0.6}(x-1)|green|dashed|0.4<x<1
```
*Projection of $\mathbf{b}=(1,3)$ onto the line $W$ spanned by $(2,1)$: $\mathbf{p}=\tfrac{5}{5}(2,1)=(2,1)$ and the residual $(-1,2)$ is perpendicular to $W$. The circle of radius $\lVert\mathbf{b}-\mathbf{p}\rVert=\sqrt5$ about $\mathbf{b}$ touches $W$ only at $\mathbf{p}$, so any other point $\mathbf{w}$ of $W$ (green segment) is farther away — Pythagoras in picture form.*

### Proof of Theorem C.10

Let $W=\operatorname{col}\mathbf{A}=\{\mathbf{A}\mathbf{x}:\mathbf{x}\in\mathbb{R}^n\}$, a subspace of $\mathbb{R}^m$ of dimension $\le n$. By Theorem C.9 there is a unique $\mathbf{p}\in W$ minimising $\lVert\mathbf{b}-\mathbf{w}\rVert$ over $\mathbf{w}\in W$. Since $J(\mathbf{x})=\lVert\mathbf{b}-\mathbf{A}\mathbf{x}\rVert^2$ and $\mathbf{A}\mathbf{x}$ ranges over exactly $W$, $\hat{\mathbf{x}}$ minimises $J$ iff $\mathbf{A}\hat{\mathbf{x}}=\mathbf{p}$. By part 1 of C.9 (uniqueness of the orthogonal projection), this holds iff $\mathbf{b}-\mathbf{A}\hat{\mathbf{x}}\perp W$, i.e. iff $\langle\mathbf{b}-\mathbf{A}\hat{\mathbf{x}},\mathbf{a}_j\rangle=0$ for each column $\mathbf{a}_j$ (orthogonality to a spanning set suffices, by linearity). Stacking these $n$ conditions,

$$
\mathbf{A}^\mathsf{T}(\mathbf{b}-\mathbf{A}\hat{\mathbf{x}})=\mathbf{0}\iff\mathbf{A}^\mathsf{T}\mathbf{A}\hat{\mathbf{x}}=\mathbf{A}^\mathsf{T}\mathbf{b}.
$$

*Consistency.* $\mathbf{p}\in W$, so $\mathbf{p}=\mathbf{A}\mathbf{x}_0$ for some $\mathbf{x}_0$, which then solves the normal equations. All minimisers share $\mathbf{A}\hat{\mathbf{x}}=\mathbf{p}$.

*Uniqueness criterion.* First, $\operatorname{null}(\mathbf{A}^\mathsf{T}\mathbf{A})=\operatorname{null}\mathbf{A}$: clearly $\mathbf{A}\mathbf{x}=\mathbf{0}\Rightarrow\mathbf{A}^\mathsf{T}\mathbf{A}\mathbf{x}=\mathbf{0}$; conversely $\mathbf{A}^\mathsf{T}\mathbf{A}\mathbf{x}=\mathbf{0}$ implies $0=\mathbf{x}^\mathsf{T}\mathbf{A}^\mathsf{T}\mathbf{A}\mathbf{x}=\lVert\mathbf{A}\mathbf{x}\rVert^2$, so $\mathbf{A}\mathbf{x}=\mathbf{0}$. The solution set of the consistent system $\mathbf{A}^\mathsf{T}\mathbf{A}\mathbf{x}=\mathbf{A}^\mathsf{T}\mathbf{b}$ is $\mathbf{x}_0+\operatorname{null}(\mathbf{A}^\mathsf{T}\mathbf{A})=\mathbf{x}_0+\operatorname{null}\mathbf{A}$. It is a single point iff $\operatorname{null}\mathbf{A}=\{\mathbf{0}\}$, i.e. by rank–nullity (C.2) iff $\operatorname{rank}\mathbf{A}=n$. In that case the $n\times n$ matrix $\mathbf{A}^\mathsf{T}\mathbf{A}$ has trivial null space, hence is invertible, and $\hat{\mathbf{x}}=(\mathbf{A}^\mathsf{T}\mathbf{A})^{-1}\mathbf{A}^\mathsf{T}\mathbf{b}$.

*QR form.* With $\mathbf{A}=\mathbf{Q}\mathbf{R}$ (Theorem C.8), $\mathbf{A}^\mathsf{T}\mathbf{A}=\mathbf{R}^\mathsf{T}\mathbf{R}$ and $\mathbf{A}^\mathsf{T}\mathbf{b}=\mathbf{R}^\mathsf{T}\mathbf{Q}^\mathsf{T}\mathbf{b}$; cancelling the invertible $\mathbf{R}^\mathsf{T}$ gives $\mathbf{R}\hat{\mathbf{x}}=\mathbf{Q}^\mathsf{T}\mathbf{b}$. $\blacksquare$

**Check.** Fitting $y=c_0+c_1x$ to $(0,1),(1,2),(2,2),(3,4)$: $\mathbf{A}^\mathsf{T}\mathbf{A}=\begin{bmatrix}4&6\\6&14\end{bmatrix}$, $\mathbf{A}^\mathsf{T}\mathbf{b}=(9,18)^\mathsf{T}$, so $c_0=c_1=0.9$. The residual $(0.1,\,0.2,\,-0.7,\,0.4)$ sums to zero and has zero inner product with $(0,1,2,3)$: it is orthogonal to both columns, exactly as the proof requires.

### Remarks

- **Finite dimension of $W$ (or completeness) is needed.** In the space of continuous functions on $[-1,1]$ with $\langle f,g\rangle=\int fg$, the subspace $W$ of polynomials is infinite-dimensional and dense, so a non-polynomial $b$ such as $\lvert t\rvert$ has no closest polynomial: every $\varepsilon>0$ is beaten by some polynomial, but $0$ is not attained. In Hilbert spaces the theorem holds for every *closed* subspace; that is the setting of Fourier series ([Ch. 15](../part-5-transforms/15-fourier-series.md)), where Bessel's inequality becomes Parseval's identity for a complete orthonormal set.
- **Rank deficiency.** If $\operatorname{rank}\mathbf{A}<n$ the fitted values $\mathbf{A}\hat{\mathbf{x}}$ are still unique but the parameters are not; the conventional choice is the minimum-norm solution $\mathbf{A}^+\mathbf{b}$ given by the SVD ([§10.7](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md)). In parameter estimation this is non-identifiability ([§24.7](../part-7-numerical-computational/24-computational-mathematics.md)).
- **Normal equations square the condition number.** $\kappa_2(\mathbf{A}^\mathsf{T}\mathbf{A})=\kappa_2(\mathbf{A})^2$ (from the SVD, §C.10), which is why QR or SVD is preferred numerically ([§23.8](../part-7-numerical-computational/23-optimization.md), [§24.3](../part-7-numerical-computational/24-computational-mathematics.md)).
- **The projector.** With full column rank, $\mathbf{P}=\mathbf{A}(\mathbf{A}^\mathsf{T}\mathbf{A})^{-1}\mathbf{A}^\mathsf{T}=\mathbf{Q}\mathbf{Q}^\mathsf{T}$ satisfies $\mathbf{P}^2=\mathbf{P}=\mathbf{P}^\mathsf{T}$; conversely any symmetric idempotent matrix is the orthogonal projector onto its column space.
- **Uses.** Least-squares curve fitting and calibration ([§9.11](../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md), [§22.6](../part-7-numerical-computational/22-numerical-methods.md)), linear regression ([Ch. 20](../part-6-probability-statistics/20-statistical-analysis.md)), the orthogonality principle of Wiener filtering ([Ch. 21](../part-6-probability-statistics/21-random-processes.md)), and Gauss–Newton steps ([§24.7](../part-7-numerical-computational/24-computational-mathematics.md)).

## C.7 Independence of eigenvectors and diagonalizability

> [!definition] Eigenpairs, multiplicities, diagonalizability
> $\lambda\in\mathbb{C}$ is an **eigenvalue** of $\mathbf{A}\in\mathbb{C}^{n\times n}$ with **eigenvector** $\mathbf{v}\neq\mathbf{0}$ if $\mathbf{A}\mathbf{v}=\lambda\mathbf{v}$; equivalently (Theorem C.1) $\det(\lambda\mathbf{I}-\mathbf{A})=0$. The **algebraic multiplicity** $a(\lambda)$ is the multiplicity of $\lambda$ as a root of the characteristic polynomial $p(s)=\det(s\mathbf{I}-\mathbf{A})$; the **geometric multiplicity** $g(\lambda)=\dim\operatorname{null}(\mathbf{A}-\lambda\mathbf{I})$. $\mathbf{A}$ is **diagonalizable** if $\mathbf{A}=\mathbf{V}\boldsymbol{\Lambda}\mathbf{V}^{-1}$ for some invertible $\mathbf{V}$ and diagonal $\boldsymbol{\Lambda}$.

> [!theorem] Theorem C.11 — Eigenvectors for distinct eigenvalues are independent
> If $\mathbf{v}_1,\dots,\mathbf{v}_k$ are eigenvectors of $\mathbf{A}$ for **pairwise distinct** eigenvalues $\lambda_1,\dots,\lambda_k$, then $\mathbf{v}_1,\dots,\mathbf{v}_k$ are linearly independent. More generally, if for each distinct eigenvalue $\lambda_i$ we take a basis of its eigenspace, the union of these bases is independent.

> [!theorem] Theorem C.12 — Diagonalizability
> For $\mathbf{A}\in\mathbb{C}^{n\times n}$:
> 1. $\mathbf{A}$ is diagonalizable if and only if it has $n$ linearly independent eigenvectors; then the columns of $\mathbf{V}$ are such eigenvectors and $\boldsymbol{\Lambda}$ holds the corresponding eigenvalues.
> 2. $1\le g(\lambda)\le a(\lambda)$ for every eigenvalue, and $\mathbf{A}$ is diagonalizable iff $g(\lambda)=a(\lambda)$ for every eigenvalue.
> 3. In particular, if $\mathbf{A}$ has $n$ distinct eigenvalues, it is diagonalizable.

### Proof of Theorem C.11

Suppose not, and among all dependent sets of eigenvectors with distinct eigenvalues choose one with the fewest vectors, $\mathbf{v}_1,\dots,\mathbf{v}_k$; then $k\ge2$ (a single eigenvector is non-zero, hence independent), and there is a relation

$$
c_1\mathbf{v}_1+\dots+c_k\mathbf{v}_k=\mathbf{0}
$$

in which every $c_i\neq0$ — if some $c_i$ were zero, the remaining $k-1$ vectors would be a smaller dependent set. Apply $\mathbf{A}$, and separately multiply by $\lambda_k$:

$$
c_1\lambda_1\mathbf{v}_1+\dots+c_k\lambda_k\mathbf{v}_k=\mathbf{0},\qquad c_1\lambda_k\mathbf{v}_1+\dots+c_k\lambda_k\mathbf{v}_k=\mathbf{0}.
$$

Subtracting eliminates $\mathbf{v}_k$:

$$
c_1(\lambda_1-\lambda_k)\mathbf{v}_1+\dots+c_{k-1}(\lambda_{k-1}-\lambda_k)\mathbf{v}_{k-1}=\mathbf{0}.
$$

The coefficients $c_i(\lambda_i-\lambda_k)$ are non-zero because $c_i\neq0$ and the eigenvalues are distinct. So $\mathbf{v}_1,\dots,\mathbf{v}_{k-1}$ is a smaller dependent set, a contradiction.

*General form.* Let $\mathcal{B}_i$ be a basis of the eigenspace $E_i=\operatorname{null}(\mathbf{A}-\lambda_i\mathbf{I})$ and suppose a combination of all the vectors in $\bigcup_i\mathcal{B}_i$ vanishes. Grouping terms by eigenvalue, it reads $\mathbf{z}_1+\dots+\mathbf{z}_k=\mathbf{0}$ with $\mathbf{z}_i\in E_i$. Those $\mathbf{z}_i$ that are non-zero are eigenvectors for distinct eigenvalues summing to zero, which the first part forbids; so every $\mathbf{z}_i=\mathbf{0}$. Each $\mathbf{z}_i$ is a combination of the independent set $\mathcal{B}_i$, so all its coefficients vanish. $\blacksquare$

### Proof of Theorem C.12

*Part 1.* Write $\mathbf{V}=[\mathbf{v}_1,\dots,\mathbf{v}_n]$ and $\boldsymbol{\Lambda}=\operatorname{diag}(\lambda_1,\dots,\lambda_n)$. Column $j$ of $\mathbf{A}\mathbf{V}$ is $\mathbf{A}\mathbf{v}_j$ and column $j$ of $\mathbf{V}\boldsymbol{\Lambda}$ is $\lambda_j\mathbf{v}_j$, so

$$
\mathbf{A}\mathbf{V}=\mathbf{V}\boldsymbol{\Lambda}\iff\mathbf{A}\mathbf{v}_j=\lambda_j\mathbf{v}_j\ \text{ for all }j .
$$

If the $\mathbf{v}_j$ are $n$ independent eigenvectors, $\mathbf{V}$ is invertible (Theorem C.1, Step 7 and rank–nullity) and $\mathbf{A}=\mathbf{V}\boldsymbol{\Lambda}\mathbf{V}^{-1}$. Conversely, if $\mathbf{A}=\mathbf{V}\boldsymbol{\Lambda}\mathbf{V}^{-1}$ then $\mathbf{A}\mathbf{V}=\mathbf{V}\boldsymbol{\Lambda}$, so the columns of $\mathbf{V}$ — non-zero and independent because $\mathbf{V}$ is invertible — are eigenvectors.

*Part 2.* $g(\lambda)\ge1$ because $\mathbf{A}-\lambda\mathbf{I}$ is singular. Let $g=g(\lambda)$, choose a basis $\mathbf{v}_1,\dots,\mathbf{v}_g$ of the eigenspace, extend it to a basis of $\mathbb{C}^n$ and let $\mathbf{S}$ be the invertible matrix with these columns. Since $\mathbf{A}\mathbf{v}_i=\lambda\mathbf{v}_i$ for $i\le g$, the first $g$ columns of $\mathbf{S}^{-1}\mathbf{A}\mathbf{S}$ are $\lambda\mathbf{e}_1,\dots,\lambda\mathbf{e}_g$:

$$
\mathbf{S}^{-1}\mathbf{A}\mathbf{S}=\begin{bmatrix}\lambda\mathbf{I}_g&\mathbf{B}\\ \mathbf{0}&\mathbf{C}\end{bmatrix}.
$$

Similar matrices have the same characteristic polynomial, by the product rule: $\det(s\mathbf{I}-\mathbf{S}^{-1}\mathbf{A}\mathbf{S})=\det\big(\mathbf{S}^{-1}(s\mathbf{I}-\mathbf{A})\mathbf{S}\big)=\det(s\mathbf{I}-\mathbf{A})$. The determinant of a block upper-triangular matrix is the product of the determinants of its diagonal blocks (expand down the first $g$ columns, or apply the uniqueness lemma of C.1 to the lower block), so

$$
p(s)=(s-\lambda)^g\det(s\mathbf{I}-\mathbf{C}),
$$

and $(s-\lambda)^g$ divides $p$: $a(\lambda)\ge g$.

Now let $\mu_1,\dots,\mu_k$ be the distinct eigenvalues. By the fundamental theorem of algebra $p(s)=\prod_i(s-\mu_i)^{a(\mu_i)}$, so $\sum_ia(\mu_i)=n$. By Theorem C.11 the union of eigenspace bases is an independent set of $\sum_ig(\mu_i)$ eigenvectors, and every eigenvector lies in some eigenspace, so the largest number of independent eigenvectors is exactly $\sum_ig(\mu_i)$. Hence, by Part 1,

$$
\mathbf{A}\ \text{diagonalizable}\iff\sum_ig(\mu_i)=n=\sum_ia(\mu_i)\iff g(\mu_i)=a(\mu_i)\ \text{for all }i,
$$

the last step because $g(\mu_i)\le a(\mu_i)$ term by term.

*Part 3.* With $n$ distinct eigenvalues each has $a=1$, so $1\le g\le a=1$ forces $g=a$. (Directly: one eigenvector per eigenvalue gives $n$ independent eigenvectors by C.11.) $\blacksquare$

### Remarks

- **Distinctness is sufficient, not necessary.** $\mathbf{I}$ has a single eigenvalue of multiplicity $n$ and is diagonal. Symmetric matrices are always diagonalizable however their eigenvalues repeat (§C.8).
- **The standard counterexample.** $\mathbf{N}=\begin{bmatrix}1&1\\0&1\end{bmatrix}$ has $a(1)=2$ but $\operatorname{null}(\mathbf{N}-\mathbf{I})=\operatorname{span}\{\mathbf{e}_1\}$, so $g(1)=1<2$: it is **defective** and not diagonalizable. Its $\dot{\mathbf{x}}=(\mathbf{N}-2\mathbf{I})\mathbf{x}$ solutions contain $te^{-t}$ — the critically damped RLC response ([§12.2](../part-4-differential-equations/12-linear-differential-equations.md)).
- **Generic behaviour.** Matrices with a repeated eigenvalue form a set of measure zero (the discriminant of $p$ vanishes), so almost every matrix met in practice is diagonalizable — but *nearly* defective matrices have ill-conditioned eigenvector matrices $\mathbf{V}$, with consequences in §C.12.
- **Field matters.** The rotation $\begin{bmatrix}0&-1\\1&0\end{bmatrix}$ has eigenvalues $\pm j$; it is diagonalizable over $\mathbb{C}$ but not over $\mathbb{R}$. The real block form of [§10.3](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md) is the real substitute.
- **Uses.** Diagonalization ([§10.3](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md)), matrix powers and $e^{\mathbf{A}t}$ ([§10.4](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md)), modal decomposition of coupled circuits and vibrating systems ([§10.8](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md), [§12.6](../part-4-differential-equations/12-linear-differential-equations.md)), and the stability theorem of §C.12.

## C.8 The spectral theorem for real symmetric matrices

> [!theorem] Theorem C.13 — Spectral theorem
> Let $\mathbf{A}\in\mathbb{R}^{n\times n}$ be symmetric, $\mathbf{A}^\mathsf{T}=\mathbf{A}$. Then
> 1. every eigenvalue of $\mathbf{A}$ is real, and to each there corresponds a real eigenvector;
> 2. eigenvectors belonging to distinct eigenvalues are orthogonal;
> 3. there is an orthogonal matrix $\mathbf{Q}$ ($\mathbf{Q}^\mathsf{T}\mathbf{Q}=\mathbf{I}$) and a real diagonal $\boldsymbol{\Lambda}$ with
> $$
> \mathbf{A}=\mathbf{Q}\boldsymbol{\Lambda}\mathbf{Q}^\mathsf{T}=\sum_{i=1}^n\lambda_i\,\mathbf{q}_i\mathbf{q}_i^\mathsf{T}.
> $$
>
> Equivalently, $\mathbb{R}^n$ has an orthonormal basis of eigenvectors of $\mathbf{A}$. For complex Hermitian $\mathbf{A}=\mathbf{A}^\mathsf{H}$ the same holds with a unitary $\mathbf{U}$: $\mathbf{A}=\mathbf{U}\boldsymbol{\Lambda}\mathbf{U}^\mathsf{H}$, $\boldsymbol{\Lambda}$ real.

### Proof

**Part 1: real eigenvalues.** By the fundamental theorem of algebra the characteristic polynomial has a root $\lambda\in\mathbb{C}$, with an eigenvector $\mathbf{x}\in\mathbb{C}^n$, $\mathbf{x}\neq\mathbf{0}$. Consider the scalar $\mathbf{x}^\mathsf{H}\mathbf{A}\mathbf{x}$. Since $\mathbf{A}$ is real and symmetric, $\mathbf{A}^\mathsf{H}=\mathbf{A}$, and the conjugate of this $1\times1$ matrix is its conjugate transpose:

$$
(\mathbf{x}^\mathsf{H}\mathbf{A}\mathbf{x})^*=(\mathbf{x}^\mathsf{H}\mathbf{A}\mathbf{x})^\mathsf{H}=\mathbf{x}^\mathsf{H}\mathbf{A}^\mathsf{H}\mathbf{x}=\mathbf{x}^\mathsf{H}\mathbf{A}\mathbf{x},
$$

so $\mathbf{x}^\mathsf{H}\mathbf{A}\mathbf{x}$ is real. But $\mathbf{x}^\mathsf{H}\mathbf{A}\mathbf{x}=\lambda\,\mathbf{x}^\mathsf{H}\mathbf{x}=\lambda\lVert\mathbf{x}\rVert^2$ with $\lVert\mathbf{x}\rVert^2>0$ real. Hence $\lambda=\mathbf{x}^\mathsf{H}\mathbf{A}\mathbf{x}/\lVert\mathbf{x}\rVert^2$ is real.

*Real eigenvector.* $\mathbf{A}-\lambda\mathbf{I}$ is now a real matrix with determinant zero. By Theorem C.1 applied over $\mathbb{R}$, it is not invertible as a real matrix, so (Step 7 of that proof, contrapositive) its columns are dependent over $\mathbb{R}$: there is a real $\mathbf{q}\neq\mathbf{0}$ with $(\mathbf{A}-\lambda\mathbf{I})\mathbf{q}=\mathbf{0}$. (Alternatively, write $\mathbf{x}=\mathbf{y}+j\mathbf{z}$ with $\mathbf{y},\mathbf{z}$ real; since $\mathbf{A}$ and $\lambda$ are real, $\mathbf{A}\mathbf{y}=\lambda\mathbf{y}$ and $\mathbf{A}\mathbf{z}=\lambda\mathbf{z}$, and $\mathbf{y},\mathbf{z}$ are not both zero.)

**Part 2: orthogonality.** Let $\mathbf{A}\mathbf{u}=\lambda\mathbf{u}$, $\mathbf{A}\mathbf{v}=\mu\mathbf{v}$ with real $\lambda\neq\mu$ and real vectors. Using symmetry,

$$
\lambda\,\mathbf{v}^\mathsf{T}\mathbf{u}=\mathbf{v}^\mathsf{T}\mathbf{A}\mathbf{u}=(\mathbf{A}^\mathsf{T}\mathbf{v})^\mathsf{T}\mathbf{u}=(\mathbf{A}\mathbf{v})^\mathsf{T}\mathbf{u}=\mu\,\mathbf{v}^\mathsf{T}\mathbf{u},
$$

so $(\lambda-\mu)\mathbf{v}^\mathsf{T}\mathbf{u}=0$ and $\mathbf{v}^\mathsf{T}\mathbf{u}=0$.

**Part 3: orthogonal diagonalization, by induction on $n$.** For $n=1$ take $\mathbf{Q}=[1]$. Assume the statement for symmetric matrices of order $n-1$, and let $\mathbf{A}$ be symmetric of order $n$.

*Step (a): split off one eigenvector.* By Part 1 there is a real eigenvalue $\lambda_1$ with a real eigenvector, which we normalise to $\lVert\mathbf{q}_1\rVert=1$. By Theorem C.7 extend $\mathbf{q}_1$ to an orthonormal basis $\mathbf{q}_1,\mathbf{p}_2,\dots,\mathbf{p}_n$ of $\mathbb{R}^n$, and set $\mathbf{P}=[\mathbf{p}_2,\dots,\mathbf{p}_n]$ ($n\times(n-1)$) and $\mathbf{Q}_1=[\mathbf{q}_1\ \mathbf{P}]$, an orthogonal matrix.

*Step (b): the similarity is block diagonal.* Compute

$$
\mathbf{Q}_1^\mathsf{T}\mathbf{A}\mathbf{Q}_1=\begin{bmatrix}\mathbf{q}_1^\mathsf{T}\mathbf{A}\mathbf{q}_1&\mathbf{q}_1^\mathsf{T}\mathbf{A}\mathbf{P}\\ \mathbf{P}^\mathsf{T}\mathbf{A}\mathbf{q}_1&\mathbf{P}^\mathsf{T}\mathbf{A}\mathbf{P}\end{bmatrix}.
$$

Here $\mathbf{q}_1^\mathsf{T}\mathbf{A}\mathbf{q}_1=\lambda_1\mathbf{q}_1^\mathsf{T}\mathbf{q}_1=\lambda_1$; $\mathbf{P}^\mathsf{T}\mathbf{A}\mathbf{q}_1=\lambda_1\mathbf{P}^\mathsf{T}\mathbf{q}_1=\mathbf{0}$ because each $\mathbf{p}_i\perp\mathbf{q}_1$; and $\mathbf{q}_1^\mathsf{T}\mathbf{A}\mathbf{P}=(\mathbf{P}^\mathsf{T}\mathbf{A}^\mathsf{T}\mathbf{q}_1)^\mathsf{T}=(\mathbf{P}^\mathsf{T}\mathbf{A}\mathbf{q}_1)^\mathsf{T}=\mathbf{0}^\mathsf{T}$ — **this is where symmetry is used a second time.** So

$$
\mathbf{Q}_1^\mathsf{T}\mathbf{A}\mathbf{Q}_1=\begin{bmatrix}\lambda_1&\mathbf{0}^\mathsf{T}\\ \mathbf{0}&\mathbf{B}\end{bmatrix},\qquad \mathbf{B}=\mathbf{P}^\mathsf{T}\mathbf{A}\mathbf{P},
$$

and $\mathbf{B}^\mathsf{T}=\mathbf{P}^\mathsf{T}\mathbf{A}^\mathsf{T}\mathbf{P}=\mathbf{B}$ is symmetric of order $n-1$.

*Step (c): induct.* By hypothesis $\mathbf{B}=\mathbf{Q}_2\boldsymbol{\Lambda}_2\mathbf{Q}_2^\mathsf{T}$ with $\mathbf{Q}_2$ orthogonal and $\boldsymbol{\Lambda}_2$ real diagonal. Put

$$
\mathbf{Q}=\mathbf{Q}_1\begin{bmatrix}1&\mathbf{0}^\mathsf{T}\\ \mathbf{0}&\mathbf{Q}_2\end{bmatrix},\qquad \boldsymbol{\Lambda}=\begin{bmatrix}\lambda_1&\mathbf{0}^\mathsf{T}\\ \mathbf{0}&\boldsymbol{\Lambda}_2\end{bmatrix}.
$$

$\mathbf{Q}$ is a product of orthogonal matrices, hence orthogonal ($(\mathbf{Q}_a\mathbf{Q}_b)^\mathsf{T}\mathbf{Q}_a\mathbf{Q}_b=\mathbf{Q}_b^\mathsf{T}\mathbf{Q}_b=\mathbf{I}$), and

$$
\mathbf{Q}^\mathsf{T}\mathbf{A}\mathbf{Q}=\begin{bmatrix}1&\mathbf{0}^\mathsf{T}\\ \mathbf{0}&\mathbf{Q}_2^\mathsf{T}\end{bmatrix}\begin{bmatrix}\lambda_1&\mathbf{0}^\mathsf{T}\\ \mathbf{0}&\mathbf{B}\end{bmatrix}\begin{bmatrix}1&\mathbf{0}^\mathsf{T}\\ \mathbf{0}&\mathbf{Q}_2\end{bmatrix}=\begin{bmatrix}\lambda_1&\mathbf{0}^\mathsf{T}\\ \mathbf{0}&\mathbf{Q}_2^\mathsf{T}\mathbf{B}\mathbf{Q}_2\end{bmatrix}=\boldsymbol{\Lambda}.
$$

Hence $\mathbf{A}=\mathbf{Q}\boldsymbol{\Lambda}\mathbf{Q}^\mathsf{T}$. Writing the product column by column gives $\mathbf{A}=\sum_i\lambda_i\mathbf{q}_i\mathbf{q}_i^\mathsf{T}$, and $\mathbf{A}\mathbf{Q}=\mathbf{Q}\boldsymbol{\Lambda}$ shows the columns $\mathbf{q}_i$ are orthonormal eigenvectors.

*Hermitian case.* Replace $\mathsf{T}$ by $\mathsf{H}$ and "orthogonal" by "unitary" throughout; Part 1 was already proved for $\mathbf{A}^\mathsf{H}=\mathbf{A}$, and Gram–Schmidt works in $\mathbb{C}^n$. $\blacksquare$

```desmos-graph
left=-1; right=1; bottom=-1; top=1
---
3x^2+2xy+3y^2=1|blue
y=x|red|dashed
y=-x|green|dashed
(0.35355,0.35355)|label:q1 axis, length 1/2
(0.5,-0.5)|label:q2 axis, length 1/sqrt 2
```
*Level curve $\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}=1$ for $\mathbf{A}=\begin{bmatrix}3&1\\1&3\end{bmatrix}$. The eigenvectors $(1,1)/\sqrt2$ (eigenvalue 4) and $(1,-1)/\sqrt2$ (eigenvalue 2) are orthogonal and lie along the principal axes, with semi-axes $1/\sqrt{\lambda_i}=0.5$ and $0.707$: in the rotated coordinates $\mathbf{y}=\mathbf{Q}^\mathsf{T}\mathbf{x}$ the form is $4y_1^2+2y_2^2$, with no cross term.*

### Remarks

- **Repeated eigenvalues are harmless here.** The induction never assumes distinct eigenvalues: it produces $n$ orthonormal eigenvectors regardless. Contrast Theorem C.12, where distinctness was the easy sufficient condition. Within an eigenspace of dimension $>1$ the orthonormal eigenvectors are not unique (any orthonormal basis of the eigenspace will do), but the spectral projectors $\sum_{\lambda_i=\lambda}\mathbf{q}_i\mathbf{q}_i^\mathsf{T}$ are.
- **Symmetry is essential, twice.** Without it eigenvalues can be complex ($\begin{bmatrix}0&-1\\1&0\end{bmatrix}$), eigenvectors can be non-orthogonal ($\begin{bmatrix}1&1\\0&2\end{bmatrix}$ has eigenvectors $(1,0)$ and $(1,1)$), and the matrix can be defective ($\begin{bmatrix}1&1\\0&1\end{bmatrix}$). In the induction, symmetry is what kills the top-right block in Step (b); for a general matrix the same construction gives only a block upper-triangular form, which iterated yields the **Schur decomposition** $\mathbf{A}=\mathbf{U}\mathbf{T}\mathbf{U}^\mathsf{H}$ with $\mathbf{T}$ upper triangular. A matrix is unitarily diagonalizable iff it is **normal**, $\mathbf{A}\mathbf{A}^\mathsf{H}=\mathbf{A}^\mathsf{H}\mathbf{A}$.
- **Rayleigh quotient.** From $\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}=\sum\lambda_iy_i^2$ with $\mathbf{y}=\mathbf{Q}^\mathsf{T}\mathbf{x}$, $\lVert\mathbf{y}\rVert=\lVert\mathbf{x}\rVert$, one gets $\lambda_{\min}\le\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}/\mathbf{x}^\mathsf{T}\mathbf{x}\le\lambda_{\max}$, with the extremes attained at the corresponding eigenvectors ([§10.5](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md)).
- **Uses.** Principal axes of quadratic forms and energy functions ([§10.5](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md)); normal modes of LC networks and mass–spring systems, where $\mathbf{M}^{-1/2}\mathbf{K}\mathbf{M}^{-1/2}$ is symmetric ([§10.8](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md)); principal component analysis and Karhunen–Loève expansions of covariance matrices ([Ch. 20](../part-6-probability-statistics/20-statistical-analysis.md), [Ch. 21](../part-6-probability-statistics/21-random-processes.md)); the convergence analysis of gradient descent on quadratics ([§23.5](../part-7-numerical-computational/23-optimization.md)). It is the foundation of §C.9 and §C.10.

## C.9 Characterisations of positive-definite matrices

> [!definition] Positive definite
> A real symmetric $\mathbf{A}\in\mathbb{R}^{n\times n}$ is **positive definite** ($\mathbf{A}\succ0$) if $\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}>0$ for every $\mathbf{x}\neq\mathbf{0}$, and **positive semidefinite** ($\mathbf{A}\succeq0$) if $\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}\ge0$ for every $\mathbf{x}$.

> [!theorem] Theorem C.14 — Tests for positive definiteness
> For real symmetric $\mathbf{A}\in\mathbb{R}^{n\times n}$ the following are equivalent:
> 1. (energy) $\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}>0$ for all $\mathbf{x}\neq\mathbf{0}$;
> 2. (eigenvalues) every eigenvalue of $\mathbf{A}$ is positive;
> 3. (Sylvester's criterion) every leading principal minor is positive: $\Delta_k=\det\mathbf{A}_k>0$, $k=1,\dots,n$;
> 4. (pivots) Gaussian elimination without row interchanges runs to completion with all pivots positive;
> 5. (Cholesky) $\mathbf{A}=\mathbf{L}\mathbf{L}^\mathsf{T}$ with $\mathbf{L}$ lower triangular with positive diagonal entries; this $\mathbf{L}$ is unique;
> 6. (Gram) $\mathbf{A}=\mathbf{R}^\mathsf{T}\mathbf{R}$ for some (possibly rectangular) $\mathbf{R}$ with linearly independent columns.

### Proof

We prove $1\Leftrightarrow2$, then the cycle $1\Rightarrow3\Rightarrow4\Rightarrow5\Rightarrow6\Rightarrow1$.

**$1\Leftrightarrow2$.** By Theorem C.13, $\mathbf{A}=\mathbf{Q}\boldsymbol{\Lambda}\mathbf{Q}^\mathsf{T}$. Substituting $\mathbf{y}=\mathbf{Q}^\mathsf{T}\mathbf{x}$, which is a bijection of $\mathbb{R}^n$ mapping non-zero vectors to non-zero vectors,

$$
\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}=\mathbf{y}^\mathsf{T}\boldsymbol{\Lambda}\mathbf{y}=\sum_{i=1}^n\lambda_iy_i^2 .
$$

If all $\lambda_i>0$, this is positive whenever $\mathbf{y}\neq\mathbf{0}$. Conversely, if $1$ holds, take $\mathbf{x}=\mathbf{q}_i$ (i.e. $\mathbf{y}=\mathbf{e}_i$): $\lambda_i=\mathbf{q}_i^\mathsf{T}\mathbf{A}\mathbf{q}_i>0$.

**$1\Rightarrow3$.** Fix $k$. For $\mathbf{z}\in\mathbb{R}^k$, $\mathbf{z}\neq\mathbf{0}$, pad it with $n-k$ zeros to $\mathbf{x}=(\mathbf{z},\mathbf{0})$; then $\mathbf{z}^\mathsf{T}\mathbf{A}_k\mathbf{z}=\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}>0$. So the symmetric matrix $\mathbf{A}_k$ satisfies $1$, hence (by $1\Rightarrow2$ applied to $\mathbf{A}_k$) has positive eigenvalues $\mu_1,\dots,\mu_k$. By the spectral theorem and the product rule, $\Delta_k=\det(\mathbf{Q}_k\mathbf{M}\mathbf{Q}_k^\mathsf{T})=\det\mathbf{M}=\mu_1\cdots\mu_k>0$, where $\det\mathbf{Q}_k\det\mathbf{Q}_k^\mathsf{T}=\det(\mathbf{Q}_k\mathbf{Q}_k^\mathsf{T})=1$.

**$3\Rightarrow4$.** All $\Delta_k\neq0$, so Theorem C.3 gives $\mathbf{A}=\mathbf{L}_1\mathbf{U}$ (elimination without interchanges succeeds), with pivots $u_{kk}=\Delta_k/\Delta_{k-1}>0$.

**$4\Rightarrow5$.** Let $\mathbf{A}=\mathbf{L}_1\mathbf{U}$ with $\mathbf{L}_1$ unit lower triangular and pivots $d_k=u_{kk}>0$. Put $\mathbf{D}=\operatorname{diag}(d_1,\dots,d_n)$ and $\mathbf{U}_1=\mathbf{D}^{-1}\mathbf{U}$, which is unit upper triangular. Then $\mathbf{A}=\mathbf{L}_1\mathbf{D}\mathbf{U}_1$ and, by symmetry,

$$
\mathbf{A}=\mathbf{A}^\mathsf{T}=\mathbf{U}_1^\mathsf{T}\,(\mathbf{D}\mathbf{L}_1^\mathsf{T}).
$$

This is a second factorisation of $\mathbf{A}$ into a unit lower triangular matrix $\mathbf{U}_1^\mathsf{T}$ and an upper triangular matrix $\mathbf{D}\mathbf{L}_1^\mathsf{T}$. Since the pivots are non-zero, $\Delta_k=d_1\cdots d_k\neq0$ and Theorem C.3 says the LU factorisation is unique; hence $\mathbf{U}_1^\mathsf{T}=\mathbf{L}_1$, i.e.

$$
\mathbf{A}=\mathbf{L}_1\mathbf{D}\mathbf{L}_1^\mathsf{T}.
$$

With $\mathbf{D}^{1/2}=\operatorname{diag}(\sqrt{d_k})$ — real because $d_k>0$ — set $\mathbf{L}=\mathbf{L}_1\mathbf{D}^{1/2}$: lower triangular with diagonal $\sqrt{d_k}>0$, and $\mathbf{L}\mathbf{L}^\mathsf{T}=\mathbf{L}_1\mathbf{D}\mathbf{L}_1^\mathsf{T}=\mathbf{A}$.

*Uniqueness of $\mathbf{L}$.* If $\mathbf{L}\mathbf{L}^\mathsf{T}=\mathbf{G}\mathbf{G}^\mathsf{T}$ with both factors lower triangular with positive diagonal, then $\mathbf{M}=\mathbf{G}^{-1}\mathbf{L}=\mathbf{G}^\mathsf{T}\mathbf{L}^{-\mathsf{T}}$. The left expression is lower triangular, the right upper triangular, so $\mathbf{M}$ is diagonal with positive diagonal; and $\mathbf{M}\mathbf{M}^\mathsf{T}=\mathbf{G}^{-1}\mathbf{L}\mathbf{L}^\mathsf{T}\mathbf{G}^{-\mathsf{T}}=\mathbf{I}$ gives $m_{ii}^2=1$, so $\mathbf{M}=\mathbf{I}$ and $\mathbf{L}=\mathbf{G}$.

**$5\Rightarrow6$.** Take $\mathbf{R}=\mathbf{L}^\mathsf{T}$; it is triangular with non-zero diagonal, hence invertible, so its columns are independent.

**$6\Rightarrow1$.** $\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}=\mathbf{x}^\mathsf{T}\mathbf{R}^\mathsf{T}\mathbf{R}\mathbf{x}=\lVert\mathbf{R}\mathbf{x}\rVert^2\ge0$, and it vanishes only if $\mathbf{R}\mathbf{x}=\mathbf{0}$, which by independence of the columns of $\mathbf{R}$ means $\mathbf{x}=\mathbf{0}$. $\blacksquare$

**Check.** $\mathbf{A}=\begin{bmatrix}4&2&0\\2&5&3\\0&3&10\end{bmatrix}$ has $\Delta_1=4$, $\Delta_2=16$, $\Delta_3=124$, so it is positive definite with pivots $4,\ 4,\ 7.75$. The Cholesky factor is $\mathbf{L}=\begin{bmatrix}2&0&0\\1&2&0\\0&1.5&\sqrt{7.75}\end{bmatrix}$ with $\sqrt{7.75}=2.7839$ (its diagonal is the square roots of the pivots, as the proof predicts), and the eigenvalues $1.941,\ 5.553,\ 11.507$ are indeed all positive with product $124$.

### Remarks

- **Symmetry is part of the hypothesis.** For non-symmetric $\mathbf{A}$, only the symmetric part $\tfrac12(\mathbf{A}+\mathbf{A}^\mathsf{T})$ enters $\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}$, and the tests come apart: $\begin{bmatrix}1&-4\\0&1\end{bmatrix}$ has positive eigenvalues and positive leading minors, yet $\mathbf{x}=(1,1)$ gives $\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}=-2$.
- **Semidefinite: leading minors are not enough.** $\mathbf{A}=\begin{bmatrix}0&0\\0&-1\end{bmatrix}$ has $\Delta_1=\Delta_2=0\ge0$ but is negative semidefinite. The correct semidefinite version of Sylvester's criterion requires **all** principal minors (every subset of rows with the same columns) to be $\ge0$; the eigenvalue and Gram tests carry over with $>$ replaced by $\ge$ and "independent columns" dropped. The proof of $3\Rightarrow1$ fails exactly because a zero pivot halts elimination.
- **Negative definite** is not "all $\Delta_k<0$" but alternating signs, $(-1)^k\Delta_k>0$, from applying the theorem to $-\mathbf{A}$ ($\det(-\mathbf{A}_k)=(-1)^k\Delta_k$).
- **Engineering reading.** For a passive RLC network, the inductance and capacitance matrices are positive definite because stored magnetic and electric energy $\tfrac12\mathbf{i}^\mathsf{T}\mathbf{L}\mathbf{i}$, $\tfrac12\mathbf{v}^\mathsf{T}\mathbf{C}\mathbf{v}$ is positive for any non-zero state. Nodal conductance matrices of connected resistive networks with a ground path are positive definite too, which is why Cholesky (half the work of LU, and stable without pivoting) can solve them ([§8.9](../part-3-linear-algebra/08-matrices-and-linear-systems.md)).
- **Uses.** Second-derivative test and convexity: a $C^2$ function is strictly convex where its Hessian is positive definite ([§23.4](../part-7-numerical-computational/23-optimization.md), [§23.7](../part-7-numerical-computational/23-optimization.md)); covariance matrices are positive semidefinite ([Ch. 19](../part-6-probability-statistics/19-probability.md)); the normal-equation matrix $\mathbf{A}^\mathsf{T}\mathbf{A}$ is positive definite exactly when $\mathbf{A}$ has full column rank, by test 6 (§C.6); weighted inner products need $\mathbf{W}\succ0$ ([§9.8](../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md)); Lyapunov stability, $\mathbf{A}^\mathsf{T}\mathbf{P}+\mathbf{P}\mathbf{A}\prec0$ with $\mathbf{P}\succ0$ ([§10.9](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md)).

## C.10 Existence of the singular value decomposition

> [!theorem] Theorem C.15 — Singular value decomposition
> Let $\mathbf{A}\in\mathbb{R}^{m\times n}$ have rank $r$. There are orthogonal $\mathbf{U}\in\mathbb{R}^{m\times m}$, $\mathbf{V}\in\mathbb{R}^{n\times n}$ and numbers $\sigma_1\ge\sigma_2\ge\dots\ge\sigma_r>0$ such that
> $$
> \mathbf{A}=\mathbf{U}\boldsymbol{\Sigma}\mathbf{V}^\mathsf{T}=\sum_{i=1}^r\sigma_i\,\mathbf{u}_i\mathbf{v}_i^\mathsf{T},\qquad \boldsymbol{\Sigma}=\begin{bmatrix}\operatorname{diag}(\sigma_1,\dots,\sigma_r)&\mathbf{0}\\ \mathbf{0}&\mathbf{0}\end{bmatrix}\in\mathbb{R}^{m\times n}.
> $$
> The $\sigma_i^2$ are the non-zero eigenvalues of $\mathbf{A}^\mathsf{T}\mathbf{A}$ (and of $\mathbf{A}\mathbf{A}^\mathsf{T}$), so the singular values are uniquely determined by $\mathbf{A}$. Moreover
> $$
> \lVert\mathbf{A}\rVert_2=\max_{\mathbf{x}\neq\mathbf{0}}\frac{\lVert\mathbf{A}\mathbf{x}\rVert_2}{\lVert\mathbf{x}\rVert_2}=\sigma_1 ,
> $$
> and if $\mathbf{A}$ is square and invertible ($r=n=m$), $\lVert\mathbf{A}^{-1}\rVert_2=1/\sigma_n$. (Complex case: $\mathsf{T}\to\mathsf{H}$, orthogonal $\to$ unitary; the $\sigma_i$ remain real and positive.)

### Proof

**Step 1: eigen-decompose the Gram matrix.** $\mathbf{A}^\mathsf{T}\mathbf{A}$ is $n\times n$ and symmetric, so by Theorem C.13 it has orthonormal eigenvectors $\mathbf{v}_1,\dots,\mathbf{v}_n$ with real eigenvalues $\mu_1,\dots,\mu_n$, which we order so that $\mu_1\ge\mu_2\ge\dots\ge\mu_n$. Each is non-negative:

$$
\mu_i=\mu_i\mathbf{v}_i^\mathsf{T}\mathbf{v}_i=\mathbf{v}_i^\mathsf{T}\mathbf{A}^\mathsf{T}\mathbf{A}\mathbf{v}_i=\lVert\mathbf{A}\mathbf{v}_i\rVert^2\ge0 .
$$

**Step 2: exactly $r$ of them are positive.** In §C.6 we showed $\operatorname{null}(\mathbf{A}^\mathsf{T}\mathbf{A})=\operatorname{null}\mathbf{A}$, so by rank–nullity $\operatorname{rank}(\mathbf{A}^\mathsf{T}\mathbf{A})=n-\dim\operatorname{null}\mathbf{A}=r$. The rank of the diagonalizable matrix $\mathbf{A}^\mathsf{T}\mathbf{A}=\mathbf{V}\mathbf{M}\mathbf{V}^\mathsf{T}$ equals the rank of $\mathbf{M}=\operatorname{diag}(\mu_i)$ (multiplying by invertible matrices does not change rank), which is the number of non-zero $\mu_i$. So $\mu_1\ge\dots\ge\mu_r>0=\mu_{r+1}=\dots=\mu_n$. Define $\sigma_i=\sqrt{\mu_i}$ for $i\le r$.

**Step 3: build orthonormal left vectors.** For $i\le r$ put $\mathbf{u}_i=\mathbf{A}\mathbf{v}_i/\sigma_i\in\mathbb{R}^m$. Then for $i,k\le r$

$$
\mathbf{u}_i^\mathsf{T}\mathbf{u}_k=\frac{\mathbf{v}_i^\mathsf{T}\mathbf{A}^\mathsf{T}\mathbf{A}\mathbf{v}_k}{\sigma_i\sigma_k}=\frac{\mu_k\,\mathbf{v}_i^\mathsf{T}\mathbf{v}_k}{\sigma_i\sigma_k}=\frac{\sigma_k^2\,\delta_{ik}}{\sigma_i\sigma_k}=\delta_{ik},
$$

so $\mathbf{u}_1,\dots,\mathbf{u}_r$ are orthonormal. By Theorem C.7 extend them to an orthonormal basis $\mathbf{u}_1,\dots,\mathbf{u}_m$ of $\mathbb{R}^m$ and let $\mathbf{U}=[\mathbf{u}_1,\dots,\mathbf{u}_m]$, $\mathbf{V}=[\mathbf{v}_1,\dots,\mathbf{v}_n]$; both are orthogonal.

**Step 4: verify $\mathbf{A}\mathbf{V}=\mathbf{U}\boldsymbol{\Sigma}$.** Column $i$ of $\mathbf{U}\boldsymbol{\Sigma}$ is $\sigma_i\mathbf{u}_i$ for $i\le r$ and $\mathbf{0}$ for $i>r$. Column $i$ of $\mathbf{A}\mathbf{V}$ is $\mathbf{A}\mathbf{v}_i$, which equals $\sigma_i\mathbf{u}_i$ for $i\le r$ by definition, and for $i>r$ satisfies $\lVert\mathbf{A}\mathbf{v}_i\rVert^2=\mu_i=0$, so $\mathbf{A}\mathbf{v}_i=\mathbf{0}$. Hence $\mathbf{A}\mathbf{V}=\mathbf{U}\boldsymbol{\Sigma}$; multiplying on the right by $\mathbf{V}^\mathsf{T}=\mathbf{V}^{-1}$ gives $\mathbf{A}=\mathbf{U}\boldsymbol{\Sigma}\mathbf{V}^\mathsf{T}$. Expanding the product by columns gives the rank-one sum.

**Step 5: uniqueness of the $\sigma_i$.** From any SVD, $\mathbf{A}^\mathsf{T}\mathbf{A}=\mathbf{V}\boldsymbol{\Sigma}^\mathsf{T}\boldsymbol{\Sigma}\mathbf{V}^\mathsf{T}$, an orthogonal diagonalization with diagonal $\sigma_1^2,\dots,\sigma_r^2,0,\dots$; similarly $\mathbf{A}\mathbf{A}^\mathsf{T}=\mathbf{U}\boldsymbol{\Sigma}\boldsymbol{\Sigma}^\mathsf{T}\mathbf{U}^\mathsf{T}$. So the $\sigma_i^2$ are the non-zero eigenvalues of either product, which depend only on $\mathbf{A}$.

**Step 6: the 2-norm.** For any $\mathbf{x}$ let $\mathbf{y}=\mathbf{V}^\mathsf{T}\mathbf{x}$; orthogonal matrices preserve length ($\lVert\mathbf{Q}\mathbf{z}\rVert^2=\mathbf{z}^\mathsf{T}\mathbf{Q}^\mathsf{T}\mathbf{Q}\mathbf{z}=\lVert\mathbf{z}\rVert^2$), so $\lVert\mathbf{y}\rVert=\lVert\mathbf{x}\rVert$ and

$$
\lVert\mathbf{A}\mathbf{x}\rVert^2=\lVert\mathbf{U}\boldsymbol{\Sigma}\mathbf{y}\rVert^2=\lVert\boldsymbol{\Sigma}\mathbf{y}\rVert^2=\sum_{i=1}^r\sigma_i^2y_i^2\le\sigma_1^2\lVert\mathbf{y}\rVert^2=\sigma_1^2\lVert\mathbf{x}\rVert^2 ,
$$

with equality at $\mathbf{x}=\mathbf{v}_1$. Hence $\lVert\mathbf{A}\rVert_2=\sigma_1$. If $\mathbf{A}$ is square and invertible, $\mathbf{A}^{-1}=\mathbf{V}\boldsymbol{\Sigma}^{-1}\mathbf{U}^\mathsf{T}$ is an SVD (after reversing the order of the columns) with singular values $1/\sigma_n\ge\dots\ge1/\sigma_1$, so $\lVert\mathbf{A}^{-1}\rVert_2=1/\sigma_n$. $\blacksquare$

**Check.** $\mathbf{A}=\begin{bmatrix}3&0\\4&5\end{bmatrix}$: $\mathbf{A}^\mathsf{T}\mathbf{A}=\begin{bmatrix}25&20\\20&25\end{bmatrix}$ has eigenvalues $45$ and $5$ with eigenvectors $(1,1)/\sqrt2$ and $(1,-1)/\sqrt2$, so $\sigma_1=3\sqrt5=6.708$, $\sigma_2=\sqrt5=2.236$, and $\mathbf{u}_1=\mathbf{A}\mathbf{v}_1/\sigma_1=(1,3)/\sqrt{10}$, $\mathbf{u}_2=(3,-1)/\sqrt{10}$, which are orthogonal as Step 3 promises. Note $\sigma_1\sigma_2=15=\lvert\det\mathbf{A}\rvert$.

```desmos-graph
left=-8; right=8; bottom=-8; top=8
---
41x^2-24xy+9y^2=225|blue
x^2+y^2=1|grey
y=3x|red|dashed
y=-\frac{x}{3}|green|dashed
(2.1213,6.364)|label:sigma1 u1
(2.1213,-0.7071)|label:sigma2 u2
(0.7071,0.7071)|label:v1
```
*The SVD as geometry: $\mathbf{A}=\begin{bmatrix}3&0\\4&5\end{bmatrix}$ maps the unit circle onto the ellipse $\mathbf{y}^\mathsf{T}(\mathbf{A}\mathbf{A}^\mathsf{T})^{-1}\mathbf{y}=1$, i.e. $41x^2-24xy+9y^2=225$. The right singular vector $\mathbf{v}_1$ on the circle lands on the long semi-axis $\sigma_1\mathbf{u}_1$ (length $6.708$); $\mathbf{v}_2$ lands on the short one, $\sigma_2\mathbf{u}_2$ (length $2.236$); the axes are perpendicular.*

### Remarks

- **No hypothesis at all.** Unlike diagonalization, the SVD exists for every matrix — rectangular, singular, defective — because it uses two different orthonormal bases and always works through the symmetric matrix $\mathbf{A}^\mathsf{T}\mathbf{A}$. The defective $\begin{bmatrix}1&1\\0&1\end{bmatrix}$, with no eigenvector basis, has singular values $(1\pm\sqrt5)/2$ in absolute value, i.e. $1.618$ and $0.618$.
- **Singular values are not eigenvalue magnitudes** except for normal matrices. For $\begin{bmatrix}0&100\\0&0\end{bmatrix}$ both eigenvalues are $0$ but $\sigma_1=100$. This gap between $\lvert\lambda\rvert$ and $\sigma$ is exactly the transient-growth phenomenon of §C.12.
- **Do not compute the SVD this way.** Forming $\mathbf{A}^\mathsf{T}\mathbf{A}$ squares the condition number and loses small singular values below $\sqrt{u}\,\sigma_1$ (with $u$ the unit round-off). Practical algorithms bidiagonalise $\mathbf{A}$ with orthogonal transformations ([§24.3](../part-7-numerical-computational/24-computational-mathematics.md)); the proof above is an existence proof.
- **Uses.** Rank and the four fundamental subspaces, the pseudoinverse $\mathbf{A}^+=\mathbf{V}\boldsymbol{\Sigma}^+\mathbf{U}^\mathsf{T}$ and minimum-norm least squares ([§10.7](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md)); best rank-$k$ approximation (Eckart–Young) behind PCA and data compression ([Ch. 20](../part-6-probability-statistics/20-statistical-analysis.md)); numerical rank ([§8.5](../part-3-linear-algebra/08-matrices-and-linear-systems.md)); and $\kappa_2=\sigma_1/\sigma_n$ in §C.11.

## C.11 The condition-number error bound

> [!definition] Induced norm and condition number
> Given a vector norm $\lVert\cdot\rVert$ on $\mathbb{F}^n$, the induced matrix norm is $\lVert\mathbf{A}\rVert=\max_{\mathbf{x}\neq\mathbf{0}}\lVert\mathbf{A}\mathbf{x}\rVert/\lVert\mathbf{x}\rVert$, and for nonsingular $\mathbf{A}$ the **condition number** is $\kappa(\mathbf{A})=\lVert\mathbf{A}\rVert\,\lVert\mathbf{A}^{-1}\rVert$ ([§8.10](../part-3-linear-algebra/08-matrices-and-linear-systems.md)).

The maximum exists because $\mathbf{x}\mapsto\lVert\mathbf{A}\mathbf{x}\rVert$ is continuous and the unit sphere of a finite-dimensional space is compact; by homogeneity the ratio takes the same values on the sphere as on all of $\mathbb{F}^n\setminus\{\mathbf{0}\}$. Two consequences are immediate from the definition: $\lVert\mathbf{A}\mathbf{x}\rVert\le\lVert\mathbf{A}\rVert\lVert\mathbf{x}\rVert$ for every $\mathbf{x}$, and $\lVert\mathbf{A}\mathbf{B}\rVert\le\lVert\mathbf{A}\rVert\lVert\mathbf{B}\rVert$ (apply the first twice). With $\mathbf{B}=\mathbf{A}^{-1}$ this gives $\kappa(\mathbf{A})\ge\lVert\mathbf{I}\rVert=1$.

> [!theorem] Theorem C.16 — Perturbation bounds for $\mathbf{A}\mathbf{x}=\mathbf{b}$
> Let $\mathbf{A}$ be nonsingular, $\mathbf{b}\neq\mathbf{0}$, and $\mathbf{A}\mathbf{x}=\mathbf{b}$.
> 1. (right-hand side) If $\mathbf{A}(\mathbf{x}+\delta\mathbf{x})=\mathbf{b}+\delta\mathbf{b}$, then
> $$
> \frac{\lVert\delta\mathbf{x}\rVert}{\lVert\mathbf{x}\rVert}\le\kappa(\mathbf{A})\,\frac{\lVert\delta\mathbf{b}\rVert}{\lVert\mathbf{b}\rVert}.
> $$
> In the 2-norm the bound is sharp: there exist $\mathbf{b}$, $\delta\mathbf{b}$ attaining equality.
> 2. (residual) For any approximate solution $\hat{\mathbf{x}}$ with residual $\mathbf{r}=\mathbf{b}-\mathbf{A}\hat{\mathbf{x}}$: $\ \lVert\mathbf{x}-\hat{\mathbf{x}}\rVert/\lVert\mathbf{x}\rVert\le\kappa(\mathbf{A})\lVert\mathbf{r}\rVert/\lVert\mathbf{b}\rVert$.
> 3. (matrix) If $(\mathbf{A}+\delta\mathbf{A})(\mathbf{x}+\delta\mathbf{x})=\mathbf{b}$, then $\ \dfrac{\lVert\delta\mathbf{x}\rVert}{\lVert\mathbf{x}+\delta\mathbf{x}\rVert}\le\kappa(\mathbf{A})\dfrac{\lVert\delta\mathbf{A}\rVert}{\lVert\mathbf{A}\rVert}$.

### Proof

**Part 1.** Subtracting $\mathbf{A}\mathbf{x}=\mathbf{b}$ from the perturbed equation gives $\mathbf{A}\,\delta\mathbf{x}=\delta\mathbf{b}$, so $\delta\mathbf{x}=\mathbf{A}^{-1}\delta\mathbf{b}$ and

$$
\lVert\delta\mathbf{x}\rVert\le\lVert\mathbf{A}^{-1}\rVert\,\lVert\delta\mathbf{b}\rVert .
$$

Also $\lVert\mathbf{b}\rVert=\lVert\mathbf{A}\mathbf{x}\rVert\le\lVert\mathbf{A}\rVert\lVert\mathbf{x}\rVert$; since $\mathbf{b}\neq\mathbf{0}$ we have $\mathbf{x}\neq\mathbf{0}$, and

$$
\frac{1}{\lVert\mathbf{x}\rVert}\le\frac{\lVert\mathbf{A}\rVert}{\lVert\mathbf{b}\rVert}.
$$

Multiplying the two inequalities (both sides non-negative) gives the bound.

*Sharpness in the 2-norm.* Take the SVD of $\mathbf{A}$ (Theorem C.15), square of order $n$, with $\sigma_n>0$. Choose $\mathbf{x}=\mathbf{v}_1$, so $\mathbf{b}=\sigma_1\mathbf{u}_1$ and $\lVert\mathbf{A}\mathbf{x}\rVert=\lVert\mathbf{A}\rVert_2\lVert\mathbf{x}\rVert$; and $\delta\mathbf{b}=\varepsilon\mathbf{u}_n$, so $\delta\mathbf{x}=\varepsilon\mathbf{v}_n/\sigma_n$ and $\lVert\delta\mathbf{x}\rVert=\lVert\mathbf{A}^{-1}\rVert_2\lVert\delta\mathbf{b}\rVert$. Both inequalities are then equalities:

$$
\frac{\lVert\delta\mathbf{x}\rVert}{\lVert\mathbf{x}\rVert}=\frac{\varepsilon}{\sigma_n},\qquad \kappa_2(\mathbf{A})\frac{\lVert\delta\mathbf{b}\rVert}{\lVert\mathbf{b}\rVert}=\frac{\sigma_1}{\sigma_n}\cdot\frac{\varepsilon}{\sigma_1}=\frac{\varepsilon}{\sigma_n}.
$$

**Part 2.** $\hat{\mathbf{x}}$ solves $\mathbf{A}\hat{\mathbf{x}}=\mathbf{b}-\mathbf{r}$ exactly; apply Part 1 with $\delta\mathbf{b}=-\mathbf{r}$ and $\delta\mathbf{x}=\hat{\mathbf{x}}-\mathbf{x}$.

**Part 3.** Expanding, $\mathbf{A}\mathbf{x}+\mathbf{A}\,\delta\mathbf{x}+\delta\mathbf{A}(\mathbf{x}+\delta\mathbf{x})=\mathbf{b}$, and $\mathbf{A}\mathbf{x}=\mathbf{b}$, so $\delta\mathbf{x}=-\mathbf{A}^{-1}\delta\mathbf{A}(\mathbf{x}+\delta\mathbf{x})$. Taking norms, $\lVert\delta\mathbf{x}\rVert\le\lVert\mathbf{A}^{-1}\rVert\lVert\delta\mathbf{A}\rVert\lVert\mathbf{x}+\delta\mathbf{x}\rVert$; divide by $\lVert\mathbf{x}+\delta\mathbf{x}\rVert$ (if it is zero, then $\mathbf{b}=\mathbf{0}$, excluded) and write $\lVert\mathbf{A}^{-1}\rVert=\kappa(\mathbf{A})/\lVert\mathbf{A}\rVert$. $\blacksquare$

**Check.** $\mathbf{A}=\begin{bmatrix}1&1\\1&1.1\end{bmatrix}$, $\mathbf{b}=(2,2.1)$, $\mathbf{x}=(1,1)$. Changing $b_2$ to $2.2$ ($\lVert\delta\mathbf{b}\rVert_2/\lVert\mathbf{b}\rVert_2=0.1/2.900=0.0345$) moves the solution to $(0,2)$, a relative change $\lVert(-1,1)\rVert/\lVert(1,1)\rVert=1$: an amplification of $29$. The theorem allows up to $\kappa_2(\mathbf{A})=42.08$, i.e. a relative change of $1.45$, so the bound holds and is of the right order. Making the $(2,2)$ entry $1.0001$ instead raises $\kappa_2$ to about $4.0\times10^4$.

```desmos-graph
left=-1; right=3; bottom=-0.5; top=3
---
y=2-x|blue
y=\frac{2.1-x}{1.1}|red
y=\frac{2.2-x}{1.1}|red|dashed
(1,1)|label:x
(0,2)|label:x + dx
```
*Ill-conditioning seen geometrically: the two equations of $\mathbf{A}\mathbf{x}=\mathbf{b}$ are nearly parallel lines (blue and red). Shifting the red line by $0.1$ in $b_2$ (dashed) slides the intersection from $(1,1)$ to $(0,2)$ — a $3.4\%$ data change producing a $100\%$ solution change.*

### Remarks

- **Sharpness in other norms.** For $\lVert\cdot\rVert_1$ and $\lVert\cdot\rVert_\infty$ the bound is also attained for suitable $\mathbf{b},\delta\mathbf{b}$ (pick a vector achieving $\lVert\mathbf{A}\rVert$ and one achieving $\lVert\mathbf{A}^{-1}\rVert$; any induced norm attains its maximum). For a particular $\mathbf{b}$ the actual amplification can be much smaller — the bound is worst case over the data.
- **Part 3 is "first order".** Bounding $\lVert\delta\mathbf{x}\rVert/\lVert\mathbf{x}\rVert$ rather than $/\lVert\mathbf{x}+\delta\mathbf{x}\rVert$ requires $\kappa(\mathbf{A})\lVert\delta\mathbf{A}\rVert/\lVert\mathbf{A}\rVert<1$, which also guarantees $\mathbf{A}+\delta\mathbf{A}$ is nonsingular; then the bound is $\kappa\epsilon_A/(1-\kappa\epsilon_A)$ with $\epsilon_A=\lVert\delta\mathbf{A}\rVert/\lVert\mathbf{A}\rVert$ ([§8.10](../part-3-linear-algebra/08-matrices-and-linear-systems.md)). Without that condition there is no bound at all: $\mathbf{A}+\delta\mathbf{A}$ may be singular.
- **Small residual does not mean small error.** Part 2 shows a residual at round-off level only guarantees about $\kappa(\mathbf{A})\cdot u$ relative accuracy. This is the rule of thumb "you lose $\log_{10}\kappa$ digits" ([§22.2](../part-7-numerical-computational/22-numerical-methods.md)), and combined with backward stability of LU with pivoting it is the basis of all error estimates for linear solvers ([§24.3](../part-7-numerical-computational/24-computational-mathematics.md)).
- **Determinant is not the measure.** $\det\mathbf{A}=0.1$ above, but scaling $\mathbf{A}$ by $10$ changes the determinant by $100$ and leaves $\kappa$ unchanged; $\kappa$ is scale- and unit-free, the determinant is neither.
- **Uses.** Sensitivity of nodal analysis to component tolerances, error analysis of least squares ($\kappa(\mathbf{A}^\mathsf{T}\mathbf{A})=\kappa_2(\mathbf{A})^2$), convergence rates of iterative solvers and gradient descent, which scale with $\kappa$ or $\sqrt\kappa$ ([§23.5](../part-7-numerical-computational/23-optimization.md), [§24.4](../part-7-numerical-computational/24-computational-mathematics.md)).

## C.12 Stability of $\dot{\mathbf{x}}=\mathbf{A}\mathbf{x}$

For $\mathbf{A}\in\mathbb{C}^{n\times n}$, the initial-value problem $\dot{\mathbf{x}}=\mathbf{A}\mathbf{x}$, $\mathbf{x}(0)=\mathbf{x}_0$ has the unique solution $\mathbf{x}(t)=e^{\mathbf{A}t}\mathbf{x}_0$, where $e^{\mathbf{A}t}=\sum_{k\ge0}\mathbf{A}^kt^k/k!$ ([§10.4](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md), [§12.7](../part-4-differential-equations/12-linear-differential-equations.md)). We use two properties of the series: $e^{\mathbf{S}\mathbf{B}\mathbf{S}^{-1}t}=\mathbf{S}e^{\mathbf{B}t}\mathbf{S}^{-1}$ (because $(\mathbf{S}\mathbf{B}\mathbf{S}^{-1})^k=\mathbf{S}\mathbf{B}^k\mathbf{S}^{-1}$ term by term), and $e^{(\lambda\mathbf{I}+\mathbf{N})t}=e^{\lambda t}e^{\mathbf{N}t}$ (because $\lambda\mathbf{I}$ commutes with $\mathbf{N}$). The **spectral abscissa** is $\alpha(\mathbf{A})=\max_i\operatorname{Re}\{\lambda_i\}$. Stability is defined as in [§10.9](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md): asymptotically stable if every solution tends to $\mathbf{0}$, marginally stable if every solution is bounded but not all tend to $\mathbf{0}$, unstable if some solution is unbounded.

> [!theorem] Theorem C.17 — Eigenvalue test, diagonalizable case
> Let $\mathbf{A}=\mathbf{V}\boldsymbol{\Lambda}\mathbf{V}^{-1}$ be diagonalizable. Then for all $t\ge0$
> $$
> \lVert e^{\mathbf{A}t}\mathbf{x}_0\rVert_2\le\kappa_2(\mathbf{V})\,e^{\alpha(\mathbf{A})t}\,\lVert\mathbf{x}_0\rVert_2 ,
> $$
> and
> 1. $\dot{\mathbf{x}}=\mathbf{A}\mathbf{x}$ is asymptotically stable $\iff$ $\operatorname{Re}\{\lambda_i\}<0$ for all $i$;
> 2. every solution is bounded on $[0,\infty)$ $\iff$ $\operatorname{Re}\{\lambda_i\}\le0$ for all $i$.
>
> If $\mathbf{A}$ is real, the same equivalences hold when only real initial conditions are allowed.

> [!theorem] Theorem C.18 — General case (Jordan form)
> For arbitrary $\mathbf{A}\in\mathbb{C}^{n\times n}$:
> 1. $\dot{\mathbf{x}}=\mathbf{A}\mathbf{x}$ is asymptotically stable $\iff$ $\operatorname{Re}\{\lambda_i\}<0$ for all $i$; in that case, for every $\beta$ with $\alpha(\mathbf{A})<\beta<0$ there is $C$ with $\lVert e^{\mathbf{A}t}\rVert\le Ce^{\beta t}$ for $t\ge0$;
> 2. all solutions are bounded $\iff$ $\operatorname{Re}\{\lambda_i\}\le0$ for all $i$ and every eigenvalue with $\operatorname{Re}\{\lambda_i\}=0$ is **semisimple** ($g=a$, i.e. all its Jordan blocks are $1\times1$).

### Proof of Theorem C.17

**The bound.** By the similarity property, $e^{\mathbf{A}t}=\mathbf{V}e^{\boldsymbol{\Lambda}t}\mathbf{V}^{-1}$, and $e^{\boldsymbol{\Lambda}t}=\operatorname{diag}(e^{\lambda_1t},\dots,e^{\lambda_nt})$ since powers of a diagonal matrix are diagonal. Equivalently, in modal coordinates $\mathbf{y}=\mathbf{V}^{-1}\mathbf{x}$ the system decouples into $\dot y_i=\lambda_iy_i$, with $y_i(t)=e^{\lambda_it}y_i(0)$. The 2-norm of a diagonal matrix is its largest entry modulus (Theorem C.15, or directly), and $\lvert e^{\lambda_it}\rvert=e^{\operatorname{Re}\{\lambda_i\}t}\le e^{\alpha t}$ for $t\ge0$. Submultiplicativity gives

$$
\lVert e^{\mathbf{A}t}\mathbf{x}_0\rVert\le\lVert\mathbf{V}\rVert\,\lVert e^{\boldsymbol{\Lambda}t}\rVert\,\lVert\mathbf{V}^{-1}\rVert\,\lVert\mathbf{x}_0\rVert\le\kappa_2(\mathbf{V})\,e^{\alpha t}\lVert\mathbf{x}_0\rVert .
$$

**Sufficiency.** If all $\operatorname{Re}\{\lambda_i\}<0$ then $\alpha<0$, so the right-hand side tends to $0$ for every $\mathbf{x}_0$. If all $\operatorname{Re}\{\lambda_i\}\le0$ then $\alpha\le0$ and every solution is bounded by $\kappa_2(\mathbf{V})\lVert\mathbf{x}_0\rVert$.

**Necessity.** Suppose some eigenvalue has $\operatorname{Re}\{\lambda_k\}\ge0$, with eigenvector $\mathbf{v}_k$. The solution from $\mathbf{x}_0=\mathbf{v}_k$ is $e^{\lambda_kt}\mathbf{v}_k$ (it satisfies the equation and the initial condition, and solutions are unique), whose norm $e^{\operatorname{Re}\{\lambda_k\}t}\lVert\mathbf{v}_k\rVert\ge\lVert\mathbf{v}_k\rVert>0$ does not tend to $0$. So asymptotic stability forces all $\operatorname{Re}\{\lambda_i\}<0$. Likewise if $\operatorname{Re}\{\lambda_k\}>0$ this solution is unbounded, so boundedness forces all $\operatorname{Re}\{\lambda_i\}\le0$.

**Real initial conditions.** Let $\mathbf{A}$ be real and $\mathbf{v}_k=\mathbf{a}+j\mathbf{b}$ with $\mathbf{a},\mathbf{b}$ real. Because $\mathbf{A}$ is real, $e^{\mathbf{A}t}$ is real, so $e^{\lambda_kt}\mathbf{v}_k=e^{\mathbf{A}t}\mathbf{a}+j\,e^{\mathbf{A}t}\mathbf{b}$ splits into the two real solutions from $\mathbf{a}$ and from $\mathbf{b}$. If both tended to $\mathbf{0}$ (or both were bounded), so would their complex combination, contradicting the previous paragraph. Hence one of the real initial conditions $\mathbf{a}$, $\mathbf{b}$ already exhibits the failure. $\blacksquare$

### Proof of Theorem C.18

We use the Jordan form ([§10.3](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md)): $\mathbf{A}=\mathbf{S}\mathbf{J}\mathbf{S}^{-1}$ with $\mathbf{J}$ block diagonal, each block $\mathbf{J}_m=\lambda\mathbf{I}_m+\mathbf{N}_m$ of size $m$, where $\mathbf{N}_m$ has ones on the superdiagonal and $\mathbf{N}_m^m=\mathbf{0}$. Since $e^{\mathbf{A}t}=\mathbf{S}e^{\mathbf{J}t}\mathbf{S}^{-1}$ and $e^{\mathbf{J}t}$ is block diagonal with blocks

$$
e^{\mathbf{J}_mt}=e^{\lambda t}\sum_{k=0}^{m-1}\frac{t^k}{k!}\mathbf{N}_m^k=e^{\lambda t}\begin{bmatrix}1&t&\tfrac{t^2}{2}&\cdots&\tfrac{t^{m-1}}{(m-1)!}\\ &1&t&\cdots&\vdots\\ &&\ddots&\ddots&\tfrac{t^2}{2}\\ &&&1&t\\ &&&&1\end{bmatrix}
$$

(the series for $e^{\mathbf{N}_mt}$ terminates), every entry of $e^{\mathbf{A}t}$ is a finite combination of functions $t^ke^{\lambda t}$, where $k<m$ for an eigenvalue whose largest Jordan block has size $m$.

**Part 1, sufficiency and the rate.** Let $\alpha<\beta<0$. For each term, $\lvert t^ke^{\lambda t}\rvert=t^ke^{\operatorname{Re}\{\lambda\}t}=e^{\beta t}\cdot t^ke^{-(\beta-\operatorname{Re}\{\lambda\})t}$, and since $\beta-\operatorname{Re}\{\lambda\}>0$ the factor $t^ke^{-(\beta-\operatorname{Re}\{\lambda\})t}$ is continuous on $[0,\infty)$ and tends to $0$ (exponentials beat polynomials), hence is bounded. A finite sum of such bounds, combined with $\lVert\mathbf{S}\rVert\lVert\mathbf{S}^{-1}\rVert$, gives $\lVert e^{\mathbf{A}t}\rVert\le Ce^{\beta t}\to0$.

**Part 1, necessity.** As in C.17: an eigenvector $\mathbf{v}$ with $\operatorname{Re}\{\lambda\}\ge0$ gives the non-decaying solution $e^{\lambda t}\mathbf{v}$.

**Part 2, sufficiency.** If $\operatorname{Re}\{\lambda\}<0$, the block terms are bounded by Part 1's argument. If $\operatorname{Re}\{\lambda\}=0$ and $\lambda$ is semisimple, its blocks are $1\times1$ and contribute $e^{\lambda t}$, of modulus $1$. So all entries of $e^{\mathbf{J}t}$, and hence of $e^{\mathbf{A}t}$, are bounded.

**Part 2, necessity.** If some $\operatorname{Re}\{\lambda\}>0$, the eigenvector solution is unbounded. If $\operatorname{Re}\{\lambda\}=0$ and $\lambda$ has a Jordan block of size $m\ge2$, let $\mathbf{s}_1,\mathbf{s}_2$ be the columns of $\mathbf{S}$ corresponding to the first two positions of that block, so $(\mathbf{A}-\lambda\mathbf{I})\mathbf{s}_1=\mathbf{0}$ and $(\mathbf{A}-\lambda\mathbf{I})\mathbf{s}_2=\mathbf{s}_1$. From the block formula, the solution with $\mathbf{x}_0=\mathbf{s}_2$ is

$$
\mathbf{x}(t)=e^{\lambda t}(\mathbf{s}_2+t\,\mathbf{s}_1),
$$

which can be verified directly: $\dot{\mathbf{x}}=\lambda\mathbf{x}+e^{\lambda t}\mathbf{s}_1$ and $\mathbf{A}\mathbf{x}=e^{\lambda t}(\lambda\mathbf{s}_2+\mathbf{s}_1+t\lambda\mathbf{s}_1)=\lambda\mathbf{x}+e^{\lambda t}\mathbf{s}_1$. Its norm is $\lVert\mathbf{s}_2+t\mathbf{s}_1\rVert\ge t\lVert\mathbf{s}_1\rVert-\lVert\mathbf{s}_2\rVert\to\infty$ (triangle inequality, C.6). For real $\mathbf{A}$, real initial conditions suffice by the splitting argument of C.17. $\blacksquare$

```desmos-graph
left=0; right=6; bottom=-0.2; top=3
---
y=10(e^{-x}-e^{-2x})|blue
y=e^{-x}|grey|dashed
y=xe^{-x}|red
(0.6931,2.5)|label:peak 2.5 at t = ln 2
```
*Asymptotic stability says nothing about the size of the transient. For $\mathbf{A}=\begin{bmatrix}-1&10\\0&-2\end{bmatrix}$ (eigenvalues $-1,-2$, eigenvector matrix with $\kappa_2(\mathbf{V})\approx102$), the first state component from $\mathbf{x}_0=(0,1)$ is $10(e^{-t}-e^{-2t})$ (blue): it rises to $2.5$ before decaying at the eigenvalue rate $e^{-t}$ (grey). The red curve $te^{-t}$ is the extra polynomial factor carried by a $2\times2$ Jordan block at $-1$: still decaying, but more slowly than $e^{-t}$.*

### Remarks

- **Where the hypotheses bite.** The single criterion "$\operatorname{Re}\{\lambda\}<0$" decides asymptotic stability whether or not $\mathbf{A}$ is diagonalizable; Jordan structure matters only on the boundary $\operatorname{Re}\{\lambda\}=0$. The double integrator $\mathbf{A}=\begin{bmatrix}0&1\\0&0\end{bmatrix}$ has both eigenvalues at $0$ but $x_1(t)=x_1(0)+t\,x_2(0)$ is unbounded: a defective imaginary-axis eigenvalue, exactly the counterexample of Theorem C.18(2). An undamped LC tank, $\begin{bmatrix}0&1\\-\omega_0^2&0\end{bmatrix}$, has simple eigenvalues $\pm j\omega_0$ and is marginally stable. The same dichotomy appears for scalar ODEs as "repeated imaginary-axis roots are unstable" ([§12.2](../part-4-differential-equations/12-linear-differential-equations.md)).
- **Transient growth and non-normality.** The constant $\kappa_2(\mathbf{V})$ in Theorem C.17 is not decoration. For the matrix in the graph, $\max_{t\ge0}\lVert e^{\mathbf{A}t}\rVert_2\approx2.56$ (near $t=0.66$) even though every eigenvalue is negative; for highly non-normal matrices the peak can be orders of magnitude, which matters for saturating amplifiers and in fluid and power-system models. For normal (e.g. symmetric) $\mathbf{A}$, $\mathbf{V}$ can be taken unitary, $\kappa_2(\mathbf{V})=1$, and $\lVert e^{\mathbf{A}t}\rVert_2=e^{\alpha t}$ exactly.
- **Discrete time.** The same proofs, with $e^{\lambda t}$ replaced by $\lambda^k$ and $\lvert e^{\lambda t}\rvert=e^{\operatorname{Re}\{\lambda\}t}$ by $\lvert\lambda^k\rvert=\lvert\lambda\rvert^k$, show that $\mathbf{x}[k+1]=\mathbf{A}\mathbf{x}[k]$ is asymptotically stable iff all $\lvert\lambda_i\rvert<1$ ([§10.4](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md), [Ch. 17](../part-5-transforms/17-z-transforms.md)); this is also the convergence condition $\rho(\mathbf{G})<1$ for stationary iterative solvers ([§24.4](../part-7-numerical-computational/24-computational-mathematics.md)).
- **Nonlinear systems.** By Lyapunov's indirect method, an equilibrium of $\dot{\mathbf{x}}=\mathbf{f}(\mathbf{x})$ is locally asymptotically stable if the Jacobian there has all $\operatorname{Re}\{\lambda\}<0$, and unstable if some $\operatorname{Re}\{\lambda\}>0$; the boundary case is inconclusive ([§10.9](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md), [Ch. 11](../part-4-differential-equations/11-ordinary-differential-equations.md)).
- **Uses.** Every stability statement in the book for LTI state models: natural response of RLC circuits ([§12.8](../part-4-differential-equations/12-linear-differential-equations.md)), pole locations and BIBO stability ([Ch. 14](../part-5-transforms/14-laplace-transforms.md)), and modal analysis ([§10.8](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md)).

## Further reading

- G. Strang, *Linear Algebra and Its Applications*, 4th ed., Cengage, 2006 — the four subspaces, LU, least squares and the SVD from the same geometric viewpoint as Part III.
- S. Axler, *Linear Algebra Done Right*, 4th ed., Springer, 2024 — determinant-free proofs of the spectral theorem and of rank–nullity in abstract form.
- R. A. Horn and C. R. Johnson, *Matrix Analysis*, 2nd ed., Cambridge University Press, 2013 — the definitive reference for positive-definite matrices, Sylvester's criterion, Schur and Jordan forms, and norms.
- L. N. Trefethen and D. Bau, *Numerical Linear Algebra*, SIAM, 1997 — QR, SVD, conditioning and the stability of Gaussian elimination, with the floating-point view of §C.3, §C.5 and §C.11.
- C. D. Meyer, *Matrix Analysis and Applied Linear Algebra*, 2nd ed., SIAM, 2023 — readable proofs throughout, including the matrix exponential and stability.
