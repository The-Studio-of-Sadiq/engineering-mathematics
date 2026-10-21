# Chapter 24 — Computational Mathematics

[← Previous](23-optimization.md) · [Table of Contents](../00-front/table-of-contents.md) · [Next →](25-mathematical-modeling.md)

> **Prerequisites:** [Ch. 8](../part-3-linear-algebra/08-matrices-and-linear-systems.md) (elimination, LU, norms, condition number, sparse storage, Jacobi and Gauss–Seidel), [Ch. 9](../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md) (least squares, QR), [Ch. 10](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md) (eigenvalues, SVD), [Ch. 19](../part-6-probability-statistics/19-probability.md)–[Ch. 20](../part-6-probability-statistics/20-statistical-analysis.md) (random variables, variance, estimation), [Ch. 22](22-numerical-methods.md) (floating point, numerical error, Newton's method).
> **Learning outcomes:** After this chapter you can
> - estimate the cost of an algorithm in flops and memory traffic, read and use big-O notation, and restructure a computation (reordering, vectorization) to change its cost;
> - decide when a symbolic computation is worth doing, recognise expression swell, and combine symbolic derivation with numerical evaluation;
> - compute QR and SVD factorisations, count their cost, and use backward error and the condition number to judge whether a computed answer can be trusted;
> - choose and tune iterative solvers (SOR, conjugate gradients, GMRES) and preconditioners, and predict their convergence from spectral information;
> - store and factor sparse matrices efficiently, controlling fill-in by reordering;
> - design Monte Carlo experiments with quantified error and variance reduction, fit nonlinear models by Gauss–Newton and Levenberg–Marquardt, propagate uncertainty by linearization (GUM), Monte Carlo and adjoint sensitivities, and make every such experiment reproducible.

Chapters 8 and 22 taught the *building blocks* of numerical computation: how a floating-point number is stored, how Gaussian elimination works, how Newton's method converges, how a system is ill-conditioned. This chapter is about *putting them together at scale*. A post-layout transistor-level simulation solves a sparse linear system with a million unknowns thousands of times; a finite-element model of a motor or an antenna produces systems with tens of millions; a yield study runs the same simulation ten thousand times with randomly perturbed component values; a device model is calibrated by fitting twenty parameters to a few thousand measured points. In every case the questions are the same. How much work will this take, and how does the work grow with problem size? Can I trust the answer, and how much of the error comes from the data, the model, the algorithm and the arithmetic? Could someone else — or I, in six months — get the same answer again?

The chapter therefore mixes three kinds of material: **cost** (Sections 24.1, 24.5), **algorithms with guaranteed accuracy** (24.2–24.4, 24.7) and **statistical computation** (24.6, 24.8), finishing with the discipline that ties them together, **reproducibility** (24.9). Short Python fragments using NumPy, SciPy and SymPy illustrate each idea; they are written to be read, not just run, and every number quoted in a worked example has been checked by such a script.

## 24.1 Algorithmic thinking

### Motivation

Two programs that compute the same mathematical quantity can differ in running time by a factor of a million. The difference rarely comes from clever micro-optimisation; it comes from the *algorithm* — the sequence of operations chosen — and from how that sequence interacts with the hardware. A discrete Fourier transform of $2^{20}$ samples computed from the definition takes about $10^{12}$ complex multiply–adds; the FFT takes about $10^{8}$ flops. Multiplying three matrices in the wrong order can cost a thousand times more than in the right order. Algorithmic thinking is the habit of asking, *before* writing code, "how does the work grow with $n$?" and "is there an equivalent formulation that does less?"

### Measuring cost

> [!definition] Flop, cost model
> A **flop** is one floating-point addition, subtraction, multiplication or division. The **arithmetic cost** of an algorithm is the number of flops it performs as a function of the problem size $n$ (number of unknowns, samples, grid points). The **memory cost** is the number of words it stores; the **communication cost** is the number of words moved between levels of the memory hierarchy (main memory, cache, registers, or between processors).

Flop counts are the classical measure and remain the right first estimate. On modern hardware, however, moving a number from main memory costs far more time than an arithmetic operation on it — typically 10–100 flop-times. An algorithm whose ratio of flops to memory words (its **arithmetic intensity**) is high, like dense matrix–matrix multiplication ($2n^3$ flops on $3n^2$ words), runs near peak speed; one with intensity $O(1)$, like a sparse matrix–vector product ($2\operatorname{nnz}$ flops on about $\operatorname{nnz}$ words), is limited by memory bandwidth no matter how fast the processor is.

### Asymptotic notation

We usually care about how cost grows for large $n$, not about exact counts.

> [!definition] Big-O, big-Omega, big-Theta
> For positive functions $f(n)$, $g(n)$:
> - $f=O(g)$ if there are constants $c>0$, $n_0$ with $f(n)\le c\,g(n)$ for all $n\ge n_0$ (grows *no faster* than $g$);
> - $f=\Omega(g)$ if $f(n)\ge c\,g(n)$ for $n\ge n_0$ (grows *at least as fast*);
> - $f=\Theta(g)$ if both hold (grows *at the same rate*).
>
> We also write $f\sim g$ if $f(n)/g(n)\to1$; this keeps the leading constant, as in "LU costs $\sim\tfrac23n^3$ flops".

Big-O hides constants, so it answers *scaling* questions ("if I double $n$, how much longer?") rather than *absolute* ones ("which of two $O(n^3)$ algorithms is faster?"). For the latter we keep the leading term with its constant, $\sim c\,n^p$. Typical growth classes in numerical work are:

| Class | Example |
|---|---|
| $O(1)$ | one evaluation of a closed-form formula |
| $O(\log n)$ | bisection to fixed absolute accuracy on an interval of length $n$; binary search |
| $O(n)$ | vector addition, inner product, tridiagonal (Thomas) solve, one sparse matrix–vector product |
| $O(n\log n)$ | FFT, sorting, nested-dissection fill for 2-D grids |
| $O(n^2)$ | dense matrix–vector product, triangular solve, direct DFT |
| $O(n^3)$ | dense LU, QR, SVD, eigenvalues, matrix–matrix product |
| $O(2^n)$, $O(n!)$ | exhaustive search over subsets or orderings; cofactor-expansion determinant |

```desmos-graph
left=0; right=20; bottom=0; top=400
---
y=x|blue
y=x\log_2(x)|green
y=x^2|orange
y=x^3/10|red
y=2^x/1000|purple
```

*Growth of $n$ (blue), $n\log_2n$ (green), $n^2$ (orange), $n^3/10$ (red) and $2^n/1000$ (purple). The constants are deliberately chosen so that the "worse" algorithms start out faster: for small $n$ any of them may win, but beyond a crossover the exponent (or the exponential) always dominates.*

> [!theorem] Empirical scaling exponent
> If a cost behaves as $T(n)\approx c\,n^p$, then measurements at two sizes give
> $$
> p\approx\frac{\ln\big(T(n_2)/T(n_1)\big)}{\ln(n_2/n_1)},
> $$
> the slope of $T$ against $n$ on log–log axes. Predicting the time at a new size is then $T(n_3)\approx T(n_2)\,(n_3/n_2)^p$.

This is the single most useful back-of-envelope tool in computational work: time a code at two modest sizes, compute $p$, and you know whether the production run will take minutes or months — and whether the code really has the complexity you think it has.

### Recurrences and divide-and-conquer

Many fast algorithms split a problem of size $n$ into $a$ subproblems of size $n/b$ and combine the results with $O(n^d)$ extra work, giving the recurrence $T(n)=a\,T(n/b)+O(n^d)$.

> [!theorem] Master theorem (simplified)
> For $T(n)=a\,T(n/b)+\Theta(n^d)$ with $a\ge1$, $b>1$:
> $$
> T(n)=\begin{cases}\Theta(n^d) & \text{if } d>\log_b a,\\ \Theta(n^d\log n) & \text{if } d=\log_b a,\\ \Theta(n^{\log_b a}) & \text{if } d<\log_b a.\end{cases}
> $$

*Sketch.* Unrolling the recursion gives $\log_b n$ levels; level $\ell$ has $a^\ell$ subproblems each costing $(n/b^\ell)^d$, i.e. total $n^d(a/b^d)^\ell$ at that level. This is a geometric series in $\ell$ with ratio $a/b^d$: if the ratio is below 1 the top level dominates, if it equals 1 all $\log_b n$ levels contribute equally, and if it exceeds 1 the $a^{\log_b n}=n^{\log_b a}$ leaves dominate. $\square$

The radix-2 FFT has $a=b=2$, $d=1$, so $T=\Theta(n\log n)$. Strassen's matrix multiplication replaces eight half-size products by seven, $a=7$, $b=2$, $d=2$, giving $\Theta(n^{\log_27})=\Theta(n^{2.807})$.

> [!example] Worked Example 24.1 — Polynomial evaluation: naive versus Horner
> Count the flops needed to evaluate $p(x)=a_0+a_1x+\dots+a_nx^n$ (a) by computing each power $x^k$ from scratch, (b) by accumulating powers, (c) by Horner's rule $p=a_0+x(a_1+x(a_2+\dots+x\,a_n))$.
>
> **(a)** The term $a_kx^k$ needs $k$ multiplications ($k-1$ for the power and one for the coefficient), so the total is $\sum_{k=1}^n k=\tfrac12n(n+1)$ multiplications plus $n$ additions: $\Theta(n^2)$.
>
> **(b)** Keeping $x^k=x\cdot x^{k-1}$ costs one multiplication per power and one per coefficient: $2n$ multiplications and $n$ additions, $3n$ flops.
>
> **(c)** Horner performs one multiplication and one addition per coefficient: $2n$ flops — optimal, since any method must at least touch every coefficient.
>
> **Check, $n=3$, $x=2$, $p=1-3x+0x^2+2x^3$.** Horner: $b_3=2$, $b_2=0+2\cdot2=4$, $b_1=-3+2\cdot4=5$, $b_0=1+2\cdot5=11$; directly $1-6+0+16=11$ ✓. For a degree-20 filter polynomial the counts are $230$, $60$ and $40$ flops. Horner is also the most accurate of the three in floating point, and its intermediate values $b_k$ are the coefficients of the deflated polynomial $p(x)/(x-x_0)$ — two benefits for the price of one loop.

> [!example] Worked Example 24.2 — Predicting a run time from two measurements
> A dense solver takes $0.21\ \mathrm{s}$ for $n=2000$ and $1.62\ \mathrm{s}$ for $n=4000$. Estimate the scaling exponent and the time for $n=20\,000$. Would $n=10^6$ (a modest post-layout netlist) be feasible?
>
> **Exponent.** $p=\ln(1.62/0.21)/\ln2=\ln7.714/0.6931=2.043/0.6931\approx2.95$ — consistent with an $O(n^3)$ LU factorisation.
>
> **Prediction.** $T(20\,000)\approx1.62\times(20\,000/4000)^{2.95}=1.62\times5^{2.95}=1.62\times115.5\approx187\ \mathrm{s}$, about three minutes.
>
> **Extrapolation to $10^6$.** $T\approx1.62\times(250)^{3}\approx2.5\times10^{7}\ \mathrm{s}$, about ten months — and the matrix alone would need $8n^2=8\times10^{12}$ bytes. Dense methods are out; the problem must be solved by exploiting sparsity (Section 24.5), where the cost for circuit matrices is typically close to $O(n^{1.1\text{–}1.5})$.
>
> **Reading.** The measured $p=2.95$ rather than exactly 3 is normal: at moderate $n$ lower-order terms and cache effects still contribute. Always measure at sizes large enough that the leading term dominates, and treat the prediction as an order of magnitude.

> [!example] Worked Example 24.3 — Associativity changes the cost
> Let $\mathbf{A},\mathbf{B}\in\mathbb{R}^{n\times n}$ and $\mathbf{x}\in\mathbb{R}^n$. Compare the cost of $(\mathbf{A}\mathbf{B})\mathbf{x}$ and $\mathbf{A}(\mathbf{B}\mathbf{x})$ for $n=1000$. Then consider $\mathbf{y}=(\mathbf{I}+\mathbf{u}\mathbf{v}^\mathsf{T})\mathbf{x}$.
>
> **Products.** A matrix–matrix product costs $2n^3$ flops ($n^2$ entries, each an inner product of length $n$: $n$ multiplications and $n-1$ additions); a matrix–vector product costs $2n^2$.
> - $(\mathbf{A}\mathbf{B})\mathbf{x}$: $2n^3+2n^2\approx2.002\times10^9$ flops.
> - $\mathbf{A}(\mathbf{B}\mathbf{x})$: $2n^2+2n^2=4\times10^6$ flops.
>
> The results are mathematically identical; the second is 500 times cheaper, and the ratio grows as $n/2$.
>
> **Rank-one update.** Forming $\mathbf{I}+\mathbf{u}\mathbf{v}^\mathsf{T}$ explicitly and multiplying costs $n^2$ (outer product) $+\,2n^2$ (product) flops and $n^2$ words of memory. Writing instead $\mathbf{y}=\mathbf{x}+\mathbf{u}(\mathbf{v}^\mathsf{T}\mathbf{x})$ costs $2n$ (inner product) $+\,2n$ (scaled add) $=4n$ flops and no extra storage. For $n=10^6$ the explicit matrix would not even fit in memory.
>
> **Rule.** Never form a matrix you only need to multiply by. The same rule appears as "never form $\mathbf{A}^{-1}$" ([Ch. 8](../part-3-linear-algebra/08-matrices-and-linear-systems.md)), "never form $\mathbf{A}^\mathsf{T}\mathbf{A}$" (Section 24.3) and "never form the dense Jacobian" (Section 24.8).

### Vectorization and the memory hierarchy

Interpreted languages such as Python or MATLAB execute an explicit loop body at perhaps $10^7$ iterations per second; the same operation expressed as a whole-array operation runs in compiled, cache-friendly library code at $10^9$–$10^{11}$ flops per second. **Vectorization** means expressing an algorithm in terms of whole-array operations — elementwise arithmetic, reductions, matrix products — so that the inner loops run inside optimised libraries (BLAS, LAPACK, FFTW).

```python
import numpy as np
rng = np.random.default_rng(0)
t = np.linspace(0.0, 1e-3, 1_000_000)          # 1 ms, 1 M samples
# loop version: ~1 s in CPython
v = np.empty_like(t)
for k in range(t.size):
    v[k] = 5.0 * (1.0 - np.exp(-t[k] / 1e-4))
# vectorized version: ~10 ms
v = 5.0 * (1.0 - np.exp(-t / 1e-4))
```

The libraries are themselves organised in three levels. **Level-1 BLAS** (vector–vector, $O(n)$ work on $O(n)$ data), **level-2** (matrix–vector, $O(n^2)$ on $O(n^2)$) and **level-3** (matrix–matrix, $O(n^3)$ on $O(n^2)$). Only level-3 operations have enough arithmetic intensity to reuse data in cache, so high-performance LU, QR and Cholesky are written as **blocked algorithms**: the matrix is partitioned into $b\times b$ tiles (with $b$ of order 100) and most of the work is recast as tile–tile products. This is why `numpy.linalg.solve` on a $4000\times4000$ matrix runs at a large fraction of the processor's peak, while a hand-written triple loop runs a thousand times slower *with the same flop count*.

> [!example] Worked Example 24.4 — Direct DFT versus FFT, and a vectorized DFT
> A spectrum is needed for a record of $N=2^{20}\approx1.05\times10^6$ samples. Compare the costs of the direct DFT $X_k=\sum_{n=0}^{N-1}x_ne^{-j2\pi kn/N}$ and the radix-2 FFT.
>
> **Direct DFT.** $N$ outputs, each a sum of $N$ complex products: $N^2\approx1.1\times10^{12}$ complex multiply–adds, about $8\times10^{12}$ real flops (a complex multiply–add is 4 real multiplications and 4 additions).
>
> **FFT.** About $5N\log_2N=5\times2^{20}\times20\approx1.05\times10^8$ real flops.
>
> **Ratio.** $8N^2/(5N\log_2N)=8N/(5\log_2 N)\approx8.4\times10^4$. At $10^{10}$ flops per second the DFT takes about 13 minutes, the FFT about 10 ms.
>
> **Vectorization is not a substitute for a better algorithm.** The DFT can be vectorized as a matrix–vector product, `X = np.exp(-2j*np.pi*np.outer(k, n)/N) @ x`, which is fast for $N=1024$ — but the matrix has $N^2$ entries ($16\ \mathrm{TB}$ in complex double precision for $N=2^{20}$). Vectorization reduces the constant; only the FFT reduces the exponent.

> [!warning] Common Pitfall
> Big-O comparisons hold only *asymptotically*. Insertion sort ($O(n^2)$) beats merge sort for $n\lesssim20$, which is why library sorts switch to it for small subarrays; Strassen's algorithm only wins for $n$ in the thousands and is less accurate; a sparse solver can be slower than a dense one for a $50\times50$ matrix. Conversely, a "fast" algorithm whose constant is small but whose exponent is large will lose eventually. Measure at realistic sizes, and report both the exponent and the constant.

> [!tip] Engineering Note
> SPICE-class simulators spend most of their time in two places: evaluating device models (an $O(n)$ loop over transistors, vectorized and parallelized in modern simulators) and solving the sparse linear system at each Newton iteration. For large digital blocks the solve dominates; for analog blocks with complex transistor models the evaluation dominates. Profiling to find out *which* — before optimising anything — is itself part of algorithmic thinking. The same holds for an FPGA or DSP implementation: a filter whose coefficients are evaluated by Horner's rule needs $n$ multipliers rather than $n(n+1)/2$, and an FFT-based convolution beats direct convolution once the filter length exceeds roughly 64 taps.

### Accuracy is part of the algorithm

Cost is not the only property of an algorithm. Two algorithms with the same flop count can deliver very different accuracy: the textbook variance formula $\tfrac1N\sum x_i^2-\bar x^2$ and the two-pass formula $\tfrac1N\sum(x_i-\bar x)^2$ both cost $O(N)$, but the first suffers catastrophic cancellation when the mean is large compared with the spread ([Ch. 22](22-numerical-methods.md)). An algorithm is judged by the triple **(cost, accuracy, robustness)** — and the next sections develop the tools (backward error, condition number, convergence rate) for the second and third.

> [!example] Worked Example 24.5 — One-pass variance: cost versus accuracy
> Ten resistor readings near $10\ \mathrm{k\Omega}$ are recorded as $x_i=10\,000+\delta_i\ \Omega$ with $\delta_i\in\{0.1,0.2,\dots,1.0\}$. Compute the variance in single precision ($\varepsilon\approx6\times10^{-8}$) by the one-pass and two-pass formulas.
>
> **Exact.** The $\delta_i$ are $0.1k$, $k=1,\dots,10$, with mean $0.55$ and population variance $0.01\times\operatorname{Var}(k)=0.01\times8.25=0.0825\ \Omega^2$.
>
> **One-pass.** $\tfrac1N\sum x_i^2\approx1.0000110\times10^8$; each $x_i^2$ is about $10^8$, and single precision represents numbers near $10^8$ only to the nearest multiple of $8$. The difference $\tfrac1N\sum x_i^2-\bar x^2$ is therefore a difference of two numbers each uncertain by several units, whose true value is $0.0825$: the result is pure rounding noise. Run in NumPy `float32`, it returns $-8.0$ — a *negative* variance.
>
> **Two-pass.** Subtracting $\bar x\approx10\,000.55$ first leaves residuals between $-0.45$ and $0.45$ with small relative error; the result is $0.08254$, which agrees to five digits with the variance of the stored `float32` values computed in double precision. (It differs from $0.0825$ in the fourth digit only because $10\,000.1$ etc. are themselves stored with errors up to $5\times10^{-4}$.)
>
> **Welford's algorithm** obtains the two-pass accuracy in one pass by updating the mean and the sum of squared deviations recursively, $\bar x_k=\bar x_{k-1}+(x_k-\bar x_{k-1})/k$, $M_k=M_{k-1}+(x_k-\bar x_{k-1})(x_k-\bar x_k)$ — the method of choice for streaming data from an ADC.

## 24.2 Symbolic versus numerical computation

### Motivation

A computer can manipulate mathematics in two quite different ways. A **computer algebra system** (CAS) such as SymPy, Mathematica or Maple manipulates *expressions*: it can differentiate $\ln(1+e^{x/V_T})$, factor $s^3+2s^2+2s+1$, or solve a nodal system with symbolic $R$ and $C$ and return a transfer function. A **numerical** library such as NumPy, SciPy or LAPACK manipulates *floating-point numbers*: it evaluates, factorises and integrates, quickly and approximately. Engineers need both, and the skill is knowing which to use for which part of a problem. Symbolic results give insight (which component controls the pole?), exactness (no rounding) and formulas that can be reused; numerical results give speed and scale. Each fails in characteristic ways, and the best workflows combine them.

### Two models of computation

> [!definition] Symbolic and numerical computation
> **Symbolic computation** represents mathematical objects exactly — integers of arbitrary length, rationals, algebraic numbers such as $\sqrt2$, symbols, and expression trees built from them — and transforms them by rewriting rules (expansion, simplification, differentiation, Gröbner bases, symbolic integration). **Numerical computation** represents real numbers by fixed-precision floating-point approximations ([Ch. 22](22-numerical-methods.md)) and evaluates algorithms whose every step is a rounded arithmetic operation. **Exact numeric** (arbitrary-precision rational or integer) arithmetic and **arbitrary-precision floating point** (e.g. `mpmath` with 50 digits) are intermediate options.

The trade-off is between three resources:

| | Symbolic / exact | Floating point |
|---|---|---|
| Cost per operation | grows with the *size* of the operands (digits, terms) | constant (one machine instruction) |
| Error | none (or explicit) | rounding, $\lesssim\varepsilon\approx1.1\times10^{-16}$ per operation in double precision |
| Output | formula valid for all parameter values | number valid for one parameter value |
| Scale | tens of unknowns | millions of unknowns |

### Expression swell

The central failure mode of symbolic computation is that intermediate and final expressions can grow without bound even when the inputs are small.

> [!definition] Expression swell
> **Expression swell** is the growth in size (number of terms, nesting depth, or number of digits) of the expressions produced during a symbolic or exact computation, so that time and memory grow much faster than the operation count of the corresponding floating-point algorithm.

Three sources recur.

1. **Combinatorial growth.** The determinant of a general symbolic $n\times n$ matrix has $n!$ terms in expanded form: 24 for $n=4$, $3\,628\,800$ for $n=10$. Floating-point LU computes the numerical value in $\tfrac23n^3\approx670$ flops for $n=10$.
2. **Digit growth.** Exact rational arithmetic roughly doubles the number of digits at each multiplication. Gaussian elimination on an integer matrix, done naively over the rationals, can produce entries whose length grows exponentially with $n$ (fraction-free Bareiss elimination keeps the growth polynomial, but still far from $O(1)$ per operation).
3. **Closed-form blow-up.** Formulas that exist in principle can be useless in practice. A quartic has closed-form roots, but SymPy's expression for a single root of $x^4+x+1$ already contains 72 operations of nested radicals; by the Abel–Ruffini theorem degree five and above have no general radical formula at all.

> [!example] Worked Example 24.6 — Digit growth in exact Newton iteration
> Apply Newton's method $x_{k+1}=\tfrac12(x_k+2/x_k)$ for $\sqrt2$ from $x_0=1$ in exact rational arithmetic, and compare with floating point.
>
> **Iterates.** $x_1=\tfrac32$, $x_2=\tfrac{17}{12}$, $x_3=\tfrac{577}{408}$, $x_4=\tfrac{665\,857}{470\,832}$, …
>
> **Errors and sizes.**
>
> | $k$ | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
> |---|---|---|---|---|---|---|---|
> | digits in numerator | 1 | 2 | 3 | 6 | 12 | 25 | 49 |
> | error $x_k-\sqrt2$ | $8.6\times10^{-2}$ | $2.5\times10^{-3}$ | $2.1\times10^{-6}$ | $1.6\times10^{-12}$ | $\approx10^{-24}$ | $\approx10^{-48}$ | $\approx10^{-97}$ |
>
> Quadratic convergence doubles the number of correct digits each step; in exact arithmetic it *also* doubles the number of digits stored, so the cost of each step doubles too. In double precision the iteration reaches full accuracy at $k=5$ and every step costs the same three flops.
>
> **Reading.** Exactness was not buying anything: after step 5 the extra digits exceed any physical meaning. When the answer is going to be used as a number, use floating point (or arbitrary precision with a *fixed* number of digits) and control the error analytically.

> [!example] Worked Example 24.7 — Exact versus floating-point inversion of Hilbert matrices
> The Hilbert matrix $\mathbf{H}_n$, $h_{ij}=1/(i+j-1)$, is the classic ill-conditioned matrix of [Ch. 8](../part-3-linear-algebra/08-matrices-and-linear-systems.md). Compare exact rational and floating-point inversion.
>
> **Exact.** $\mathbf{H}_n^{-1}$ has *integer* entries, computable exactly in SymPy. Their size grows rapidly: the largest entry is $192$ for $n=3$, $179\,200$ for $n=5$, $4.25\times10^{9}$ for $n=8$ and $3.48\times10^{12}$ for $n=10$.
>
> **Floating point.** The 2-norm condition numbers are $\kappa_2(\mathbf{H}_5)=4.8\times10^5$, $\kappa_2(\mathbf{H}_8)=1.5\times10^{10}$, $\kappa_2(\mathbf{H}_{10})=1.6\times10^{13}$. The worst-case rule of thumb "relative error $\lesssim\kappa\varepsilon$" allows about $2\times10^{-6}$ for $n=8$, $2\times10^{-3}$ for $n=10$ and order 1 for $n=12$ ($\kappa_2\approx2\times10^{16}$); the measured relative errors of `numpy.linalg.inv` against the exact inverse are $2.1\times10^{-9}$, $8.8\times10^{-5}$ and $3.1\times10^{-2}$ — better than the bound, as is usual, but growing at the same exponential rate.
>
> **Which wins?** If the matrix entries are *exact data* (as in a proof, or a combinatorial identity), exact arithmetic is the right tool for small $n$. If the entries come from measurements with three significant figures, exact inversion is meaningless: the problem itself amplifies data errors by $\kappa$, and no arithmetic can recover what the data do not contain. The cure is reformulation (an orthogonal polynomial basis instead of monomials, Section 24.3), not more digits.

### What symbolic computation does well

- **Deriving formulas once, evaluating many times.** A symbolic transfer function, a Jacobian for Newton's method, or a sensitivity $\partial y/\partial p$ derived symbolically and then compiled to numerical code (SymPy `lambdify`, code generation to C) is exact in form and fast in evaluation.
- **Insight into parameter dependence.** A symbolic pole location shows which component moves it; a series expansion shows the leading behaviour.
- **Removing cancellation before it happens.** Rewriting an expression analytically (e.g. $1-\cos x=2\sin^2(x/2)$) eliminates catastrophic cancellation that no amount of numerical care can fix afterwards.
- **Verification.** Symbolic results provide exact references against which numerical code is tested (Section 24.9).

### What numerical computation does well

- Problems with **no closed form**: roots of degree-5 polynomials, nonlinear ODEs, most integrals of measured or tabulated data.
- **Large problems**: $10^6$ unknowns, where any symbolic representation is hopeless.
- Problems whose **data are inexact anyway**: measurements, tolerances, fitted models.
- Problems whose answer is wanted **only as a number** at specific parameter values.

> [!example] Worked Example 24.8 — Cancellation removed symbolically
> Evaluate $f(x)=\dfrac{1-\cos x}{x^2}$ at $x=10^{-4}$ and $x=10^{-8}$ in double precision. Fix it with a symbolic rewrite.
>
> **Direct evaluation.** At $x=10^{-8}$, $\cos x=1-5\times10^{-17}$ rounds to exactly $1$ in double precision (the spacing of doubles just below 1 is $1.1\times10^{-16}$), so the numerator is $0$ and $f=0$. The true value is $\tfrac12$: 100% error. At $x=10^{-4}$ the direct formula gives $0.499\,999\,996\,96$, relative error $6\times10^{-9}$ — eight digits lost to cancellation.
>
> **Symbolic rewrite.** $1-\cos x=2\sin^2(x/2)$, so $f(x)=\dfrac{2\sin^2(x/2)}{x^2}=\dfrac12\Big(\dfrac{\sin(x/2)}{x/2}\Big)^2$. There is no subtraction; the evaluation gives $0.5$ at $x=10^{-8}$ and $0.499\,999\,999\,583$ at $x=10^{-4}$, correct to all digits shown.
>
> **Series check.** SymPy's `series((1-cos(x))/x**2, x, 0, 6)` returns $\tfrac12-\tfrac{x^2}{24}+\tfrac{x^4}{720}+O(x^6)$; at $x=10^{-4}$, $\tfrac12-4.17\times10^{-10}=0.499\,999\,999\,583$ ✓.
>
> **Lesson.** The symbolic step costs nothing at run time and removes the error at its source. This is the typical successful division of labour: algebra *before* the numbers.

> [!example] Worked Example 24.9 — A hybrid workflow: symbolic transfer function, numerical evaluation
> A two-section RC ladder (series $R$, shunt $C$, series $R$, shunt $C$) is driven by $v_\text{in}$. Derive $H(s)=V_2/V_\text{in}$ symbolically, then evaluate $\lvert H(j\omega)\rvert$ numerically over 1000 frequencies.
>
> **Symbolic nodal analysis.** With node voltages $V_1$ (first capacitor) and $V_2$ (output) and $V_\text{in}=1$,
> $$
> \begin{bmatrix}\tfrac2R+sC & -\tfrac1R\\ -\tfrac1R & \tfrac1R+sC\end{bmatrix}\begin{bmatrix}V_1\\V_2\end{bmatrix}=\begin{bmatrix}\tfrac1R\\0\end{bmatrix}.
> $$
>
> ```python
> import sympy as sp, numpy as np
> s, R, C = sp.symbols('s R C', positive=True)
> G = sp.Matrix([[2/R + s*C, -1/R], [-1/R, 1/R + s*C]])
> V = G.LUsolve(sp.Matrix([1/R, 0]))
> H = sp.factor(sp.simplify(V[1]))      # 1/(C**2*R**2*s**2 + 3*C*R*s + 1)
> Hnum = sp.lambdify((s, R, C), H, 'numpy')
> w = np.logspace(1, 6, 1000)
> mag = np.abs(Hnum(1j*w, 1e3, 100e-9))  # vectorized over w
> ```
>
> **Result.** $H(s)=\dfrac{1}{R^2C^2s^2+3RCs+1}$. The symbolic form tells us immediately that the poles are at $s=-\dfrac{3\pm\sqrt5}{2RC}$, i.e. $-0.382/RC$ and $-2.618/RC$ (the same golden-ratio numbers that appear in infinite ladders): the loading of the first section by the second splits the double pole of two buffered sections apart by a factor $2.618^2\approx6.85$.
>
> **Numerical check.** With $R=1\ \mathrm{k\Omega}$, $C=100\ \mathrm{nF}$ ($RC=10^{-4}\ \mathrm{s}$), at $\omega=10^4\ \mathrm{rad/s}$: $sRC=j$, so $H=1/(-1+3j+1)=1/(3j)$, $\lvert H\rvert=1/3$ and the phase is $-90^\circ$ — what `Hnum(1j*1e4, 1e3, 100e-9)` returns.
>
> **Why hybrid?** The symbolic solve of a $2\times2$ system is instant and gives insight; repeating it for a 40-node filter would produce an enormous rational function that is slower to evaluate, and less accurate, than a numerical sparse solve at each frequency. The rule: go symbolic for small subsystems whose *structure* you want to understand; go numerical for the full system.

### Choosing between them

> [!tip] Engineering Note — a decision checklist
> 1. **Is a formula the deliverable** (a design equation, a sensitivity expression, a model to hand to someone else)? Use symbolic computation, on the smallest subsystem that captures the effect.
> 2. **Are there more than about ten coupled unknowns, or nonlinear equations without closed-form solution?** Use numerical methods.
> 3. **Will the expression be evaluated many times?** Derive symbolically, simplify, then generate numerical code (`lambdify`, `codegen`) — or use automatic differentiation (Section 24.8), which gives exact derivatives without symbolic expression swell.
> 4. **Is the computation losing digits to cancellation?** Rewrite symbolically first.
> 5. **Are the data exact integers or rationals and the size small** (combinatorics, coding theory, exact filter coefficients)? Exact arithmetic is appropriate.
> 6. **Is the problem ill-conditioned?** Neither exact arithmetic nor more digits fixes ill-conditioned *data*; reformulate.

> [!warning] Common Pitfall — floats leaking into symbolic code
> In SymPy, `sp.Rational(1, 10)` is exactly $\tfrac1{10}$, but `0.1` is the binary float $0.1000000000000000055\ldots$ and turns the whole expression into floating-point arithmetic, silently losing exactness (and `sp.nsimplify` may not recover the intended rational). Similarly, `sp.sqrt(2)` is exact but `np.sqrt(2)` is not. Conversely, never call `simplify` inside a loop over thousands of parameter values: simplification is heuristic and can take seconds; derive once, `lambdify`, then loop numerically.

> [!warning] Common Pitfall — symbolic solutions with hidden assumptions
> A CAS returns *generic* answers. Solving $ax=b$ symbolically gives $x=b/a$, valid only for $a\ne0$; the symbolic inverse of a matrix ignores the parameter values at which the determinant vanishes. When a symbolic result will be evaluated across a parameter range, check the denominators (SymPy `fraction`, `solve(..., dict=True)` with assumptions) or verify numerically at the extremes of the range.

## 24.3 Numerical linear algebra

### Motivation

Almost every large computation in engineering ends in linear algebra: Newton's method for a nonlinear circuit solves a linear system per iteration; a least-squares fit factorises a tall matrix; a modal analysis computes eigenvalues; a model-order reduction or a beamformer computes an SVD. [Chapter 8](../part-3-linear-algebra/08-matrices-and-linear-systems.md) developed Gaussian elimination, LU and Cholesky factorisation, norms and the condition number; [Chapter 9](../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md) introduced Gram–Schmidt and QR for least squares; [Chapter 10](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md) introduced the SVD. This section asks the computational questions those chapters postponed: *how* are QR and the SVD computed in floating point, what do they cost, and in what precise sense are the answers "correct"?

### Forward error, backward error and stability

Let $f$ be the exact map from data to answer (for a linear system, $\mathbf{b}\mapsto\mathbf{A}^{-1}\mathbf{b}$) and $\hat f$ the computed version.

> [!definition] Forward error, backward error
> The **forward error** of a computed answer $\hat{\mathbf{y}}=\hat f(\mathbf{x})$ is $\lVert\hat{\mathbf{y}}-f(\mathbf{x})\rVert$ (absolute) or $\lVert\hat{\mathbf{y}}-f(\mathbf{x})\rVert/\lVert f(\mathbf{x})\rVert$ (relative). The **backward error** is the size of the smallest perturbation of the *data* that makes the computed answer exact:
> $$
> \eta(\hat{\mathbf{y}})=\min\Big\{\frac{\lVert\Delta\mathbf{x}\rVert}{\lVert\mathbf{x}\rVert}\ :\ \hat{\mathbf{y}}=f(\mathbf{x}+\Delta\mathbf{x})\Big\}.
> $$
> An algorithm is **backward stable** if for every input its backward error is $O(\varepsilon)$ — that is, it computes the exact answer to a slightly perturbed problem, with the perturbation comparable to rounding the data.

The two are linked by the condition number $\kappa$ of the *problem*:

> [!theorem] Rule of thumb: forward error ≤ condition × backward error
> To first order,
> $$
> \frac{\lVert\hat{\mathbf{y}}-\mathbf{y}\rVert}{\lVert\mathbf{y}\rVert}\ \lesssim\ \kappa\cdot\eta.
> $$
> In particular, a backward-stable algorithm applied to a problem with condition number $\kappa$ delivers a relative forward error of about $\kappa\varepsilon$: roughly $\log_{10}\kappa$ significant digits are lost, and no algorithm can do systematically better with data known only to precision $\varepsilon$.

This separates responsibilities cleanly. **Conditioning** is a property of the problem — the circuit, the measurement geometry, the choice of basis — and is the engineer's business. **Stability** is a property of the algorithm and is the library writer's business. A backward-stable algorithm is "as good as the data deserve".

For linear systems the backward error has a computable formula (Rigal–Gaches): for an approximate solution $\hat{\mathbf{x}}$ with residual $\mathbf{r}=\mathbf{b}-\mathbf{A}\hat{\mathbf{x}}$, the smallest normwise perturbation $\Delta\mathbf{A}$ with $(\mathbf{A}+\Delta\mathbf{A})\hat{\mathbf{x}}=\mathbf{b}$ has
$$
\eta_A(\hat{\mathbf{x}})=\frac{\lVert\mathbf{r}\rVert}{\lVert\mathbf{A}\rVert\,\lVert\hat{\mathbf{x}}\rVert}.
$$
A *small residual* therefore certifies a *small backward error* — but, as [Ch. 8](../part-3-linear-algebra/08-matrices-and-linear-systems.md) showed, not a small forward error.

> [!example] Worked Example 24.10 — Small residual, large error
> $\mathbf{A}=\begin{bmatrix}1&1\\1&1.0001\end{bmatrix}$, $\mathbf{b}=\begin{bmatrix}2\\2.0001\end{bmatrix}$, exact solution $\mathbf{x}=(1,1)^\mathsf{T}$. A candidate answer is $\hat{\mathbf{x}}=(2,0)^\mathsf{T}$. Compute its residual, backward error and forward error in the $\infty$-norm, and check the rule of thumb.
>
> **Residual.** $\mathbf{A}\hat{\mathbf{x}}=(2,2)^\mathsf{T}$, so $\mathbf{r}=(0,\ 0.0001)^\mathsf{T}$, $\lVert\mathbf{r}\rVert_\infty=10^{-4}$.
>
> **Backward error.** $\lVert\mathbf{A}\rVert_\infty=2.0001$, $\lVert\hat{\mathbf{x}}\rVert_\infty=2$: $\eta=10^{-4}/(2.0001\times2)=2.5\times10^{-5}$. The wrong answer is the *exact* answer to a problem whose matrix differs from $\mathbf{A}$ in the fifth significant digit.
>
> **Forward error.** $\hat{\mathbf{x}}-\mathbf{x}=(1,-1)^\mathsf{T}$, relative error $1/1=100\%$.
>
> **Condition number.** $\mathbf{A}^{-1}=10^4\begin{bmatrix}1.0001&-1\\-1&1\end{bmatrix}$, $\lVert\mathbf{A}^{-1}\rVert_\infty=2.0001\times10^4$, so $\kappa_\infty=4.0004\times10^4$. Then $\kappa\eta=4.0\times10^4\times2.5\times10^{-5}=1.0$ ✓ — the bound is attained.
>
> **Reading.** If the matrix entries are resistor values known to 0.01%, both $(1,1)$ and $(2,0)$ are equally valid answers. The singular values $\sigma_1=2.00005$, $\sigma_2=5.0\times10^{-5}$ show why: the direction $\mathbf{v}_2\approx(1,-1)/\sqrt2$ is almost in the null space.

### Householder QR

Gram–Schmidt ([Ch. 9](../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md)) builds $\mathbf{Q}$ column by column by subtracting projections. In floating point, *classical* Gram–Schmidt loses orthogonality in proportion to $\kappa(\mathbf{A})^2\varepsilon$ and *modified* Gram–Schmidt in proportion to $\kappa(\mathbf{A})\varepsilon$. Production codes instead apply a sequence of orthogonal transformations to $\mathbf{A}$ that introduce zeros below the diagonal — orthogonal *triangularisation* rather than triangular *orthogonalisation*.

> [!definition] Householder reflector
> For a nonzero vector $\mathbf{v}\in\mathbb{R}^m$, the **Householder reflector**
> $$
> \mathbf{H}=\mathbf{I}-2\frac{\mathbf{v}\mathbf{v}^\mathsf{T}}{\mathbf{v}^\mathsf{T}\mathbf{v}}
> $$
> is symmetric and orthogonal ($\mathbf{H}^\mathsf{T}\mathbf{H}=\mathbf{H}^2=\mathbf{I}$); it reflects vectors in the hyperplane orthogonal to $\mathbf{v}$. Given $\mathbf{a}\ne\mathbf{0}$, the choice
> $$
> \mathbf{v}=\mathbf{a}+\operatorname{sign}(a_1)\lVert\mathbf{a}\rVert_2\,\mathbf{e}_1
> $$
> gives $\mathbf{H}\mathbf{a}=-\operatorname{sign}(a_1)\lVert\mathbf{a}\rVert_2\,\mathbf{e}_1$: all entries but the first are zeroed.

*Derivation.* $\mathbf{H}\mathbf{a}=\mathbf{a}-2\mathbf{v}(\mathbf{v}^\mathsf{T}\mathbf{a})/(\mathbf{v}^\mathsf{T}\mathbf{v})$. With $\mathbf{v}=\mathbf{a}+\alpha\mathbf{e}_1$, $\alpha=\operatorname{sign}(a_1)\lVert\mathbf{a}\rVert$, we have $\mathbf{v}^\mathsf{T}\mathbf{a}=\lVert\mathbf{a}\rVert^2+\alpha a_1$ and $\mathbf{v}^\mathsf{T}\mathbf{v}=\lVert\mathbf{a}\rVert^2+2\alpha a_1+\alpha^2=2(\lVert\mathbf{a}\rVert^2+\alpha a_1)$, so $2\mathbf{v}^\mathsf{T}\mathbf{a}/\mathbf{v}^\mathsf{T}\mathbf{v}=1$ and $\mathbf{H}\mathbf{a}=\mathbf{a}-\mathbf{v}=-\alpha\mathbf{e}_1$. The sign choice makes $\lvert v_1\rvert=\lvert a_1\rvert+\lVert\mathbf{a}\rVert$ a sum of same-signed quantities, avoiding cancellation when $\mathbf{a}$ is nearly parallel to $\mathbf{e}_1$. $\square$

Applying reflectors $\mathbf{H}_1,\dots,\mathbf{H}_n$ to zero columns $1,\dots,n$ below the diagonal in turn gives $\mathbf{H}_n\cdots\mathbf{H}_1\mathbf{A}=\mathbf{R}$, i.e. $\mathbf{A}=\mathbf{Q}\mathbf{R}$ with $\mathbf{Q}=\mathbf{H}_1\cdots\mathbf{H}_n$. The reflectors are never formed as matrices: $\mathbf{H}\mathbf{x}=\mathbf{x}-\beta\mathbf{v}(\mathbf{v}^\mathsf{T}\mathbf{x})$ with $\beta=2/\mathbf{v}^\mathsf{T}\mathbf{v}$ costs $4m$ flops (Worked Example 24.3 again), and $\mathbf{Q}$ is stored implicitly as the list of vectors $\mathbf{v}_k$.

> [!theorem] Cost and stability of Householder QR
> For $\mathbf{A}\in\mathbb{R}^{m\times n}$, $m\ge n$, Householder QR costs $\sim2mn^2-\tfrac23n^3$ flops ($\tfrac43n^3$ for square $\mathbf{A}$, twice LU). It is backward stable: the computed $\hat{\mathbf{R}}$ is the exact triangular factor of $\mathbf{A}+\Delta\mathbf{A}$ with $\lVert\Delta\mathbf{A}\rVert\le c\,mn\,\varepsilon\lVert\mathbf{A}\rVert$, and the computed $\hat{\mathbf{Q}}$ is orthogonal to working precision, $\lVert\hat{\mathbf{Q}}^\mathsf{T}\hat{\mathbf{Q}}-\mathbf{I}\rVert=O(\varepsilon)$, *independently of* $\kappa(\mathbf{A})$.

*Cost sketch.* Step $k$ applies a reflector of length $m-k+1$ to the remaining $n-k+1$ columns, at $4(m-k+1)(n-k+1)$ flops. Summing, $\sum_{k=1}^n4(m-k)(n-k)\approx4\int_0^n(m-t)(n-t)\,dt=4(mn^2-\tfrac12(m+n)n^2+\tfrac13n^3)=2mn^2-\tfrac23n^3$. $\square$

**Givens rotations** $\begin{bmatrix}c&s\\-s&c\end{bmatrix}$ with $c^2+s^2=1$ zero one entry at a time. They cost more for dense matrices (about $3mn^2-n^3$) but are ideal when only a few entries need zeroing — sparse or banded matrices, or *updating* a QR factorisation when one new measurement row arrives (recursive least squares).

> [!example] Worked Example 24.11 — Householder QR by hand
> Compute the QR factorisation of $\mathbf{A}=\begin{bmatrix}3&1\\4&2\end{bmatrix}$ with one Householder reflector.
>
> **Reflector.** $\mathbf{a}=(3,4)^\mathsf{T}$, $\lVert\mathbf{a}\rVert=5$, $\operatorname{sign}(a_1)=+1$, so $\mathbf{v}=(3+5,\ 4)^\mathsf{T}=(8,4)^\mathsf{T}$, $\mathbf{v}^\mathsf{T}\mathbf{v}=80$, $\beta=2/80=0.025$.
>
> **First column.** $\mathbf{H}\mathbf{a}=-5\mathbf{e}_1=(-5,0)^\mathsf{T}$ by construction. (Check: $\mathbf{v}^\mathsf{T}\mathbf{a}=24+16=40$, $\mathbf{a}-0.025\times40\,\mathbf{v}=(3,4)-(8,4)=(-5,0)$ ✓.)
>
> **Second column.** $\mathbf{x}=(1,2)^\mathsf{T}$: $\mathbf{v}^\mathsf{T}\mathbf{x}=8+8=16$, $\mathbf{H}\mathbf{x}=(1,2)-0.025\times16\,(8,4)=(1-3.2,\ 2-1.6)=(-2.2,\ 0.4)$.
>
> **Result.**
> $$
> \mathbf{R}=\begin{bmatrix}-5&-2.2\\0&0.4\end{bmatrix},\qquad \mathbf{Q}=\mathbf{H}=\mathbf{I}-0.025\,\mathbf{v}\mathbf{v}^\mathsf{T}=\begin{bmatrix}-0.6&-0.8\\-0.8&0.6\end{bmatrix}.
> $$
>
> **Checks.** $\mathbf{Q}\mathbf{R}$: first row $(-0.6)(-5)=3$, $(-0.6)(-2.2)+(-0.8)(0.4)=1.32-0.32=1$ ✓; second row $4$, $1.76+0.24=2$ ✓. $\det\mathbf{A}=2$, $\det\mathbf{Q}\det\mathbf{R}=(-1)(-2)=2$ ✓ (a reflector has determinant $-1$). `numpy.linalg.qr` returns exactly these factors; the negative diagonal of $\mathbf{R}$ is the Householder sign convention, not an error.

> [!example] Worked Example 24.12 — Loss of orthogonality: Gram–Schmidt versus Householder
> The Läuchli matrix with $\epsilon=10^{-8}$,
> $$
> \mathbf{L}=\begin{bmatrix}1&1&1\\\epsilon&0&0\\0&\epsilon&0\\0&0&\epsilon\end{bmatrix},
> $$
> has nearly parallel columns ($\kappa_2\approx1.7\times10^8$). Orthogonalise its columns by classical Gram–Schmidt (CGS), modified Gram–Schmidt (MGS) and Householder QR, and measure the loss of orthogonality $\lVert\hat{\mathbf{Q}}^\mathsf{T}\hat{\mathbf{Q}}-\mathbf{I}\rVert_F$ (Frobenius norm).
>
> **CGS by hand.** $1+\epsilon^2$ rounds to 1, so $\mathbf{q}_1=(1,\epsilon,0,0)$. Then $\mathbf{a}_2-(\mathbf{q}_1^\mathsf{T}\mathbf{a}_2)\mathbf{q}_1=(0,-\epsilon,\epsilon,0)$, $\mathbf{q}_2=(0,-1,1,0)/\sqrt2$. For $\mathbf{a}_3$ CGS projects the *original* column on each $\mathbf{q}_i$: $\mathbf{q}_1^\mathsf{T}\mathbf{a}_3=1$, $\mathbf{q}_2^\mathsf{T}\mathbf{a}_3=0$, giving $(0,-\epsilon,0,\epsilon)$, $\mathbf{q}_3=(0,-1,0,1)/\sqrt2$. But $\mathbf{q}_2^\mathsf{T}\mathbf{q}_3=\tfrac12$ — the "orthonormal" vectors are at $60^\circ$!
>
> **MGS.** Projects the *updated* vector: after removing $\mathbf{q}_1$, $\mathbf{a}_3$ becomes $(0,-\epsilon,0,\epsilon)$, and removing its component along $\mathbf{q}_2$ ($=\epsilon/\sqrt2$) gives $(0,-\epsilon/2,-\epsilon/2,\epsilon)$, orthogonal to $\mathbf{q}_2$.
>
> **Measured** (NumPy, double precision):
>
> | Method | $\lVert\hat{\mathbf{Q}}^\mathsf{T}\hat{\mathbf{Q}}-\mathbf{I}\rVert_F$ |
> |---|---|
> | classical Gram–Schmidt | $0.71$ |
> | modified Gram–Schmidt | $1.2\times10^{-8}$ |
> | Householder | $3.5\times10^{-16}$ |
>
> MGS loses orthogonality like $\kappa\varepsilon$ ($1.7\times10^8\times1.1\times10^{-16}\approx2\times10^{-8}$ ✓); Householder is orthogonal to machine precision regardless of $\kappa$.

### Least squares: normal equations versus QR

For the overdetermined problem $\min\lVert\mathbf{A}\mathbf{x}-\mathbf{b}\rVert_2$ ([Ch. 9](../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md)), the normal equations $\mathbf{A}^\mathsf{T}\mathbf{A}\mathbf{x}=\mathbf{A}^\mathsf{T}\mathbf{b}$ solved by Cholesky cost $\sim mn^2+\tfrac13n^3$ flops, about half of QR. But
$$
\kappa_2(\mathbf{A}^\mathsf{T}\mathbf{A})=\kappa_2(\mathbf{A})^2,
$$
because the singular values of $\mathbf{A}^\mathsf{T}\mathbf{A}$ are the squares of those of $\mathbf{A}$. Forming $\mathbf{A}^\mathsf{T}\mathbf{A}$ in floating point squares the condition number *before* any solving begins: if $\kappa(\mathbf{A})=10^8$, the normal equations have $\kappa=10^{16}$ and no digits survive.

The QR route avoids this: with $\mathbf{A}=\mathbf{Q}\mathbf{R}$ (reduced, $\mathbf{Q}\in\mathbb{R}^{m\times n}$),
$$
\lVert\mathbf{A}\mathbf{x}-\mathbf{b}\rVert_2^2=\lVert\mathbf{R}\mathbf{x}-\mathbf{Q}^\mathsf{T}\mathbf{b}\rVert_2^2+\lVert(\mathbf{I}-\mathbf{Q}\mathbf{Q}^\mathsf{T})\mathbf{b}\rVert_2^2,
$$
so the solution is $\mathbf{R}\hat{\mathbf{x}}=\mathbf{Q}^\mathsf{T}\mathbf{b}$, a triangular solve whose sensitivity is governed by $\kappa(\mathbf{A})$ (plus a term $\kappa^2\lVert\mathbf{r}\rVert$ that only matters when the residual is large, and is intrinsic to the problem).

> [!example] Worked Example 24.13 — A line fit by QR
> Fit $V=x_1+x_2I$ to $(I,V)=(1,2.1),(2,3.9),(3,6.2),(4,7.8)$ by QR.
>
> **Factorisation.** For $\mathbf{A}=\begin{bmatrix}1&1&1&1\\1&2&3&4\end{bmatrix}^\mathsf{T}$, Householder QR gives
> $$
> \mathbf{R}=\begin{bmatrix}-2&-5\\0&-2.2361\end{bmatrix},\qquad \mathbf{Q}^\mathsf{T}\mathbf{b}=\begin{bmatrix}-10\\-4.3380\end{bmatrix}.
> $$
> (The first column of $\mathbf{Q}$ is $-\tfrac12(1,1,1,1)$, so $(\mathbf{Q}^\mathsf{T}\mathbf{b})_1=-\tfrac12\sum V_i=-10$; the second is $\propto(-3,-1,1,3)$, the centred $I$ values.)
>
> **Back substitution.** $x_2=-4.3380/-2.2361=1.940$; $x_1=(-10+5\times1.940)/(-2)=0.150$.
>
> **Check against the closed form.** $\bar I=2.5$, $\bar V=5$, $\sum(I-\bar I)(V-\bar V)=(-1.5)(-2.9)+(-0.5)(-1.1)+(0.5)(1.2)+(1.5)(2.8)=9.7$, $\sum(I-\bar I)^2=5$, slope $1.94$ ✓, intercept $5-1.94\times2.5=0.15$ ✓.
>
> **Conditioning.** $\kappa_2(\mathbf{A})=7.47$, $\kappa_2(\mathbf{A}^\mathsf{T}\mathbf{A})=55.8=7.47^2$ ✓. Harmless here; the next example is not.

> [!example] Worked Example 24.14 — Polynomial fitting: when the normal equations fail
> Fit a degree-$d$ polynomial in the monomial basis to $m=50$ equally spaced points on $[0,1]$, with data generated from the polynomial whose coefficients are all 1 (so the exact answer is known). Compare normal equations and QR.
>
> | degree $d$ | $\kappa_2(\mathbf{V})$ | rel. error, normal equations | rel. error, QR |
> |---|---|---|---|
> | 6 | $2.0\times10^4$ | $2.7\times10^{-9}$ | $5.3\times10^{-13}$ |
> | 10 | $2.0\times10^7$ | $1.0\times10^{-3}$ | $7.1\times10^{-11}$ |
> | 14 | $2.3\times10^{10}$ | $6.0$ (no correct digits) | $1.4\times10^{-6}$ |
>
> **Reading.** Each column follows the $\kappa\varepsilon$ rule: QR errors track $\kappa(\mathbf{V})\varepsilon$, normal-equation errors track $\kappa(\mathbf{V})^2\varepsilon$ (for $d=14$, $\kappa^2\varepsilon\approx6\times10^4$: hopeless). Even QR is losing ten digits at $d=14$, which is a property of the *problem*: the monomials $t^{13}$ and $t^{14}$ are nearly parallel on $[0,1]$. The real cure is a better basis — Chebyshev or Legendre polynomials on $[-1,1]$, for which $\kappa$ grows only slowly with $d$ — exactly the reformulation advice of Worked Example 24.7.

### The SVD computationally

The SVD $\mathbf{A}=\mathbf{U}\boldsymbol{\Sigma}\mathbf{V}^\mathsf{T}$ ([Ch. 10](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md)) is computed in two phases. **Bidiagonalisation**: Householder reflectors applied alternately from the left and the right reduce $\mathbf{A}$ to upper bidiagonal form $\mathbf{B}=\mathbf{U}_1^\mathsf{T}\mathbf{A}\mathbf{V}_1$ (Golub–Kahan), at $\sim4mn^2-\tfrac43n^3$ flops. **Diagonalisation**: an implicitly shifted QR iteration (or a divide-and-conquer method) drives the superdiagonal of $\mathbf{B}$ to zero, at $O(n^2)$ flops per sweep for the values. Forming $\mathbf{U}$ and $\mathbf{V}$ explicitly multiplies the cost several times.

Never compute the SVD from the eigenvalues of $\mathbf{A}^\mathsf{T}\mathbf{A}$: as with the normal equations this squares the condition number, and small singular values below $\sqrt{\varepsilon}\,\sigma_1$ are lost completely.

> [!theorem] Accuracy of computed singular values
> The computed singular values $\hat\sigma_i$ of a backward-stable SVD satisfy $\lvert\hat\sigma_i-\sigma_i\rvert\le c\,\varepsilon\,\sigma_1$. Singular values are perfectly conditioned in the *absolute* sense (Weyl's inequality: $\lvert\sigma_i(\mathbf{A}+\mathbf{E})-\sigma_i(\mathbf{A})\rvert\le\lVert\mathbf{E}\rVert_2$), so large singular values are computed to high relative accuracy and singular values below $\varepsilon\sigma_1$ are indistinguishable from zero.

This makes the SVD the reliable tool for **numerical rank**: the number of singular values above a tolerance $\tau$, typically $\tau=\max(m,n)\,\varepsilon\,\sigma_1$ for exact data, or the noise level $\lVert\mathbf{E}\rVert$ for measured data. The **truncated SVD** solution of a least-squares problem keeps only the $k$ components with $\sigma_i>\tau$:
$$
\hat{\mathbf{x}}_k=\sum_{i=1}^{k}\frac{\mathbf{u}_i^\mathsf{T}\mathbf{b}}{\sigma_i}\mathbf{v}_i .
$$
Discarding the small $\sigma_i$ removes the components that would divide noise by a tiny number — a form of *regularisation*.

> [!example] Worked Example 24.15 — Truncated SVD stabilises a near-singular solve
> Return to $\mathbf{A}=\begin{bmatrix}1&1\\1&1.0001\end{bmatrix}$. The right-hand side is measured as $\mathbf{b}=(2,\ 2.0002)^\mathsf{T}$ — the true $(2,2.0001)$ plus an error of $10^{-4}$ in the second entry. Compare the exact solve with the rank-1 truncated SVD.
>
> **Exact solve.** $x_1+x_2=2$, $x_1+1.0001x_2=2.0002$ give $0.0001x_2=0.0002$: $\mathbf{x}=(0,2)^\mathsf{T}$, far from the true $(1,1)$.
>
> **SVD.** $\sigma_1=2.00005$, $\sigma_2=4.99988\times10^{-5}$; $\mathbf{v}_1\approx(0.70709,\ 0.70712)$, $\mathbf{v}_2\approx(-0.70712,\ 0.70709)$ (up to sign). The coefficients are $\mathbf{u}_1^\mathsf{T}\mathbf{b}=2.82857$ and $\mathbf{u}_2^\mathsf{T}\mathbf{b}=7.07\times10^{-5}$ (signs as returned by NumPy).
>
> - Component 1 contributes $(2.82857/2.00005)\mathbf{v}_1=(1.0000,\ 1.0000)$.
> - Component 2 contributes $(7.07\times10^{-5}/5.0\times10^{-5})\mathbf{v}_2\approx1.414\,\mathbf{v}_2=(-1,\ 1)$ — almost entirely the amplified measurement error.
>
> **Truncated.** Keeping $k=1$: $\hat{\mathbf{x}}_1=(1.0000,\ 1.00005)$, within $5\times10^{-5}$ of the truth.
>
> **Reading.** With data error of order $10^{-4}$, any singular value below about $10^{-4}$ cannot be trusted; here $\sigma_2=5\times10^{-5}$ is below that level and is discarded. Choosing the truncation level is a *modelling* decision based on the noise, not a numerical one.

### Cost summary

| Task | Algorithm | Flops (leading term) | Stable? |
|---|---|---|---|
| Solve $\mathbf{A}\mathbf{x}=\mathbf{b}$, general | LU, partial pivoting | $\tfrac23n^3$ | in practice (growth factor) |
| Solve, SPD | Cholesky | $\tfrac13n^3$ | yes |
| Solve, extra right-hand side | triangular solves with stored factors | $2n^2$ | yes |
| Least squares $m\times n$ | normal equations + Cholesky | $mn^2+\tfrac13n^3$ | no: $\kappa^2$ |
| Least squares $m\times n$ | Householder QR | $2mn^2-\tfrac23n^3$ | yes |
| Singular values only | Golub–Kahan + QR iteration | $4mn^2-\tfrac43n^3$ | yes |
| Full SVD, square | with $\mathbf{U},\mathbf{V}$ | $\approx21n^3$ | yes |
| Symmetric eigenvalues | tridiagonalisation + QR/divide & conquer | $\tfrac43n^3$ ($\approx9n^3$ with vectors) | yes |
| Nonsymmetric eigenvalues | Hessenberg + shifted QR | $\approx10n^3$ ($\approx25n^3$ with vectors) | yes (but eigenvalues may be ill-conditioned) |
| Determinant | from LU | $\tfrac23n^3$ | yes (but beware overflow; use $\log\lvert\det\rvert$) |

*(Approximate counts from Golub and Van Loan; constants for the iterative phases depend on the matrix.)*

> [!theorem] Gaussian elimination with partial pivoting is backward stable in practice
> GEPP computes $\hat{\mathbf{x}}$ with $(\mathbf{A}+\Delta\mathbf{A})\hat{\mathbf{x}}=\mathbf{b}$, $\lVert\Delta\mathbf{A}\rVert_\infty\le c\,n^3\rho_n\varepsilon\lVert\mathbf{A}\rVert_\infty$, where the **growth factor** $\rho_n=\max_{i,j,k}\lvert a_{ij}^{(k)}\rvert/\max_{i,j}\lvert a_{ij}\rvert$ measures how large the entries become during elimination. The worst case is $\rho_n=2^{n-1}$ — attained by Wilkinson's matrix (1 on the diagonal and in the last column, $-1$ below the diagonal), for which $n=20$ gives $\rho=2^{19}=524\,288$ — but matrices arising in practice almost never show growth beyond $\rho\approx10$.

### Eigenvalues: a computational note

Eigenvalues are roots of the characteristic polynomial, but **no one computes them that way**: polynomial roots are violently ill-conditioned functions of the coefficients (Wilkinson's polynomial $\prod_{k=1}^{20}(x-k)$: changing one coefficient by $2^{-23}$ moves roots into the complex plane). Instead the **QR algorithm** iterates $\mathbf{A}_k=\mathbf{Q}_k\mathbf{R}_k$, $\mathbf{A}_{k+1}=\mathbf{R}_k\mathbf{Q}_k$ (an orthogonal similarity) on a Hessenberg or tridiagonal form with shifts, converging to triangular (Schur) form. For one dominant eigenvalue of a huge sparse matrix the **power method** $\mathbf{v}_{k+1}=\mathbf{A}\mathbf{v}_k/\lVert\mathbf{A}\mathbf{v}_k\rVert$ suffices, converging at the rate $\lvert\lambda_2/\lambda_1\rvert$ per step; Lanczos and Arnoldi methods (the eigenvalue cousins of the CG and GMRES solvers of Section 24.4) find a few extreme eigenvalues of matrices with millions of rows.

> [!example] Worked Example 24.16 — Power iteration
> Apply the power method with Rayleigh quotient $\lambda^{(k)}=\mathbf{v}_k^\mathsf{T}\mathbf{A}\mathbf{v}_k$ to $\mathbf{A}=\begin{bmatrix}2&1\\1&3\end{bmatrix}$ from $\mathbf{v}_0=(1,0)$.
>
> **Iterates.** The Rayleigh quotients (using the current normalised vector) are $2,\ 3,\ 3.5,\ 3.6,\ 3.6154,\dots$, and the vectors $(0.894,0.447)$, $(0.707,0.707)$, $(0.6,0.8)$, $(0.555,0.832)$, $(0.537,0.844)$.
>
> **Exact.** $\lambda=(5\pm\sqrt5)/2=3.6180,\ 1.3820$; dominant eigenvector $\propto(0.526,0.851)$.
>
> **Rate.** $\lvert\lambda_2/\lambda_1\rvert=0.382$: the vector error shrinks by $0.382$ per step and, for a symmetric matrix, the Rayleigh-quotient error by $0.382^2=0.146$. The errors $0.118$, $0.018$, $0.0027$ in the last three quotients fall by about $0.15$ per step ✓.

> [!warning] Common Pitfall
> "The residual is tiny, so the answer is right" is false for ill-conditioned problems (Worked Example 24.10), and "the answer changed in the fifth digit when I switched libraries, so one of them is buggy" is usually false too — both may be backward stable answers to a problem with $\kappa\approx10^{11}$. Before blaming an algorithm, estimate $\kappa$ (`numpy.linalg.cond`, or the cheap LAPACK estimate from an existing factorisation) and compare the observed discrepancy with $\kappa\varepsilon$.

> [!tip] Engineering Note
> Scaling matters. A modified nodal analysis matrix mixing conductances in siemens ($10^{-3}$) with unit entries for voltage sources, or a state-space model mixing millivolts and kiloamps, can have an enormous condition number that is purely an artefact of units. **Equilibration** — scaling rows and columns by diagonal matrices so that entries are of comparable size, $\mathbf{D}_r\mathbf{A}\mathbf{D}_c$ — often reduces $\kappa$ by many orders of magnitude at negligible cost, and LAPACK's expert drivers (`dgesvx`) do it automatically. It is the matrix version of choosing sensible units, the topic of nondimensionalisation in [Ch. 25](25-mathematical-modeling.md).

## 24.4 Iterative solvers

### Motivation

A direct solver (LU, Cholesky) delivers the solution of $\mathbf{A}\mathbf{x}=\mathbf{b}$ in a fixed number of operations. For a 3-D finite-element model of a transformer with $n=10^7$ unknowns, however, even the best fill-reducing ordering (Section 24.5) leaves factors with $10^{10}$ or more nonzeros — too much memory and too many flops. **Iterative methods** never factor $\mathbf{A}$. They need only the ability to compute $\mathbf{A}\mathbf{v}$ for a given vector $\mathbf{v}$ — $O(\operatorname{nnz})$ work and no extra storage — and improve an approximate solution step by step. Their cost is (number of iterations) $\times$ (cost per iteration), and the whole art lies in making the first factor small. [Chapter 8](../part-3-linear-algebra/08-matrices-and-linear-systems.md) introduced the Jacobi and Gauss–Seidel iterations and their convergence criterion; here we develop the general theory, over-relaxation, the Krylov methods (conjugate gradients and GMRES) that dominate modern practice, and preconditioning.

### Stationary iterations: the general picture

Split $\mathbf{A}=\mathbf{M}-\mathbf{N}$ with $\mathbf{M}$ easy to invert. Then $\mathbf{A}\mathbf{x}=\mathbf{b}$ is equivalent to $\mathbf{M}\mathbf{x}=\mathbf{N}\mathbf{x}+\mathbf{b}$, suggesting
$$
\mathbf{x}^{(k+1)}=\mathbf{M}^{-1}\big(\mathbf{N}\mathbf{x}^{(k)}+\mathbf{b}\big)=\mathbf{x}^{(k)}+\mathbf{M}^{-1}\mathbf{r}^{(k)},\qquad \mathbf{r}^{(k)}=\mathbf{b}-\mathbf{A}\mathbf{x}^{(k)}.
$$
The second form shows the structure: *correct the current iterate by an approximate solve of the residual equation*. With $\mathbf{A}=\mathbf{D}+\mathbf{L}_A+\mathbf{U}_A$ (diagonal, strictly lower, strictly upper), Jacobi is $\mathbf{M}=\mathbf{D}$ and Gauss–Seidel is $\mathbf{M}=\mathbf{D}+\mathbf{L}_A$.

> [!theorem] Convergence of stationary iterations
> The error $\mathbf{e}^{(k)}=\mathbf{x}^{(k)}-\mathbf{x}$ satisfies $\mathbf{e}^{(k+1)}=\mathbf{T}\mathbf{e}^{(k)}$ with **iteration matrix** $\mathbf{T}=\mathbf{M}^{-1}\mathbf{N}=\mathbf{I}-\mathbf{M}^{-1}\mathbf{A}$. The iteration converges for every $\mathbf{x}^{(0)}$ if and only if the spectral radius $\rho(\mathbf{T})<1$, and asymptotically
> $$
> \lVert\mathbf{e}^{(k)}\rVert\approx C\,\rho(\mathbf{T})^k,\qquad k_\epsilon\approx\frac{\ln\epsilon}{\ln\rho(\mathbf{T})}
> $$
> iterations reduce the error by a factor $\epsilon$. The quantity $R=-\ln\rho(\mathbf{T})$ is the **asymptotic rate of convergence**.

(Proof in [Ch. 8](../part-3-linear-algebra/08-matrices-and-linear-systems.md): $\mathbf{T}^k\to\mathbf{0}$ iff $\rho(\mathbf{T})<1$, which follows from the Jordan form of [Ch. 10](../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md).) When $\rho$ is close to 1, $R\approx1-\rho$ and the iteration count is $\approx\ln(1/\epsilon)/(1-\rho)$: halving $1-\rho$ doubles the work.

Two caveats are worth stating. First, $\rho(\mathbf{T})<1$ controls the *asymptotic* rate; for non-normal $\mathbf{T}$ the error can grow for many steps before decaying. Second, the iteration stops on the residual, not the error: since $\mathbf{e}=-\mathbf{A}^{-1}\mathbf{r}$,
$$
\frac{\lVert\mathbf{e}^{(k)}\rVert}{\lVert\mathbf{x}\rVert}\le\kappa(\mathbf{A})\,\frac{\lVert\mathbf{r}^{(k)}\rVert}{\lVert\mathbf{b}\rVert},
$$
so a relative-residual tolerance of $10^{-8}$ guarantees only $10^{-8}\kappa$ in the error.

### The model problem

The standard test case is the discrete Laplacian. In one dimension, $-u''=f$ on $(0,1)$ with $n$ interior points and spacing $h=1/(n+1)$ gives the tridiagonal matrix $\mathbf{A}=\operatorname{tridiag}(-1,2,-1)$ — also the nodal matrix of a uniform resistive ladder with grounded ends. Its eigenvalues are
$$
\lambda_k=2-2\cos(k\pi h)=4\sin^2\!\Big(\frac{k\pi h}{2}\Big),\qquad k=1,\dots,n,
$$
with eigenvectors $\sin(k\pi ih)$. The Jacobi iteration matrix $\mathbf{I}-\tfrac12\mathbf{A}$ has eigenvalues $\cos(k\pi h)$, so
$$
\rho_J=\cos(\pi h)\approx1-\tfrac12\pi^2h^2 .
$$
The convergence factor approaches 1 as the grid is refined: the number of Jacobi iterations grows like $h^{-2}\propto n^2$. The slowest modes are the *smooth* ones ($k$ small) — the iteration damps oscillatory errors quickly and smooth errors very slowly, an observation that is the starting point of **multigrid** methods.

### Successive over-relaxation (SOR)

Gauss–Seidel moves each unknown to the value that satisfies its own equation. Over-relaxation moves it *further* in the same direction by a factor $\omega$:
$$
x_i^{(k+1)}=(1-\omega)x_i^{(k)}+\frac{\omega}{a_{ii}}\Big(b_i-\sum_{j<i}a_{ij}x_j^{(k+1)}-\sum_{j>i}a_{ij}x_j^{(k)}\Big).
$$
In matrix form $\mathbf{M}_\omega=\tfrac1\omega\mathbf{D}+\mathbf{L}_A$ and
$$
\mathbf{T}_\omega=(\mathbf{D}+\omega\mathbf{L}_A)^{-1}\big((1-\omega)\mathbf{D}-\omega\mathbf{U}_A\big).
$$
$\omega=1$ is Gauss–Seidel; $\omega>1$ is over-relaxation, $\omega<1$ under-relaxation.

> [!theorem] SOR convergence (Kahan, Ostrowski–Reich, Young)
> 1. (**Kahan**) $\rho(\mathbf{T}_\omega)\ge\lvert\omega-1\rvert$, so SOR can converge only for $0<\omega<2$.
> 2. (**Ostrowski–Reich**) If $\mathbf{A}$ is symmetric positive definite, SOR converges for every $0<\omega<2$.
> 3. (**Young**) If $\mathbf{A}$ is *consistently ordered* (true for tridiagonal matrices and for the 5-point Laplacian in natural or red–black ordering) and the Jacobi eigenvalues are real with $\rho_J<1$, then the optimal relaxation factor and the resulting spectral radius are
> $$
> \omega_\text{opt}=\frac{2}{1+\sqrt{1-\rho_J^2}},\qquad \rho(\mathbf{T}_{\omega_\text{opt}})=\omega_\text{opt}-1 .
> $$
> Moreover $\rho(\mathbf{T}_1)=\rho_J^2$: Gauss–Seidel converges exactly twice as fast as Jacobi.

*Sketch of 1.* $\det\mathbf{T}_\omega=\det(\mathbf{D}+\omega\mathbf{L}_A)^{-1}\det((1-\omega)\mathbf{D}-\omega\mathbf{U}_A)=\det(\mathbf{D})^{-1}(1-\omega)^n\det\mathbf{D}=(1-\omega)^n$ because both factors are triangular. The determinant is the product of the eigenvalues, so at least one has modulus $\ge\lvert1-\omega\rvert$. $\square$

*Sketch of 3.* For consistently ordered matrices, each Jacobi eigenvalue $\mu$ and SOR eigenvalue $\lambda$ are related by Young's equation $(\lambda+\omega-1)^2=\lambda\omega^2\mu^2$. Solving this quadratic in $\sqrt\lambda$ gives $\sqrt\lambda=\tfrac12\big(\omega\mu\pm\sqrt{\omega^2\mu^2-4(\omega-1)}\big)$. For $\omega$ below the value making the discriminant zero for $\mu=\rho_J$, the largest $\lvert\lambda\rvert$ decreases with $\omega$; above it the roots are complex with $\lvert\lambda\rvert=\omega-1$, increasing with $\omega$. The minimum is where the discriminant vanishes, $\omega^2\rho_J^2=4(\omega-1)$, whose root in $(1,2)$ is $\omega_\text{opt}$. $\square$

```desmos-graph
left=1; right=2; bottom=0; top=1
---
y=\left(\frac{0.95x+\sqrt{0.9025x^{2}-4(x-1)}}{2}\right)^{2}|blue|1<x<1.5241
y=x-1|blue|1.5241<x<2
y=0.9025|grey|dashed
(1.5241,0.5241)|label:optimum
```

*Spectral radius of SOR against $\omega$ for a consistently ordered matrix with $\rho_J=0.95$. At $\omega=1$ (Gauss–Seidel) $\rho=0.9025$ (dashed); the minimum $\rho=0.524$ at $\omega_\text{opt}=1.524$ is a sharp corner. Overestimating $\omega$ costs little (slope 1 on the right); underestimating it costs a lot (infinite slope on the left).*

> [!example] Worked Example 24.17 — SOR on the three-node network
> Apply SOR with $\omega=1.2$ to the nodal system of [Ch. 8](../part-3-linear-algebra/08-matrices-and-linear-systems.md), Worked Example 8.39:
> $$
> \begin{bmatrix}4&-2&-1\\-2&4&-1\\-1&-1&4\end{bmatrix}\mathbf{v}=\begin{bmatrix}8\\2\\1\end{bmatrix},\qquad \mathbf{v}=(4,3,2)^\mathsf{T},
> $$
> from $\mathbf{v}^{(0)}=\mathbf{0}$, and compare spectral radii.
>
> **Iteration 1.**
> - $v_1$: Gauss–Seidel value $\tfrac14(8)=2$; SOR $v_1=-0.2\cdot0+1.2\cdot2=2.4$.
> - $v_2$: GS value $\tfrac14(2+2\cdot2.4+0)=1.7$; $v_2=1.2\times1.7=2.04$.
> - $v_3$: GS value $\tfrac14(1+2.4+2.04)=1.36$; $v_3=1.2\times1.36=1.632$.
>
> **Iterations 2 and 3.** $\mathbf{v}^{(2)}=(3.6336,\ 2.8618,\ 1.9222)$, $\mathbf{v}^{(3)}=(3.9670,\ 2.9845,\ 2.0010)$. Maximum errors: $1.6$, $0.366$, $0.033$ — against $2.0$, $0.969$, $0.456$ for Gauss–Seidel.
>
> **Spectral radii** (computed from $\mathbf{T}_\omega$):
>
> | $\omega$ | 1.0 | 1.1 | 1.16 | 1.2 | 1.3 | 1.5 |
> |---|---|---|---|---|---|---|
> | $\rho(\mathbf{T}_\omega)$ | 0.472 | 0.350 | 0.230 | 0.271 | 0.367 | 0.555 |
>
> The optimum is near $\omega=1.16$, where the iterations needed for $10^{-6}$ fall from 18 (GS) to 9.5. Young's formula with $\rho_J=0.683$ would predict $\omega_\text{opt}=1.155$; it applies only approximately because this matrix is not consistently ordered, but it is a good starting estimate.

> [!example] Worked Example 24.18 — The model problem: Jacobi, Gauss–Seidel and optimal SOR
> For the 1-D model problem with $n=50$ interior points, how many iterations does each method need to reduce the error by $10^{-6}$? Repeat for $n=100$.
>
> **$n=50$.** $h=1/51$, $\rho_J=\cos(\pi/51)=0.998103$. Then
> - Jacobi: $k=\ln10^{-6}/\ln0.998103\approx7277$;
> - Gauss–Seidel: $\rho=\rho_J^2=0.996210$, $k\approx3639$;
> - SOR: $\omega_\text{opt}=2/(1+\sqrt{1-\rho_J^2})=1.884$, $\rho=0.884$, $k\approx112$.
>
> **$n=100$.** $\rho_J=0.999516$: Jacobi $28\,554$, Gauss–Seidel $14\,277$, SOR ($\omega_\text{opt}=1.940$) $222$.
>
> **Scaling.** Doubling $n$ multiplies Jacobi/GS work per unit accuracy by 4 (since $1-\rho\propto h^2$) but optimal SOR only by 2 (since $1-\rho_\text{SOR}\approx2\pi h$). Optimal SOR turns an $O(n^2)$ iteration count into $O(n)$ — the same improvement, we shall see, that conjugate gradients achieve *without* needing to know $\rho_J$.

### Krylov subspace methods

Stationary methods use only the latest residual. **Krylov methods** use the whole history: after $k$ steps from $\mathbf{x}_0$ (residual $\mathbf{r}_0$) they choose the best approximation from the affine space $\mathbf{x}_0+\mathcal{K}_k$, where
$$
\mathcal{K}_k(\mathbf{A},\mathbf{r}_0)=\operatorname{span}\{\mathbf{r}_0,\ \mathbf{A}\mathbf{r}_0,\ \mathbf{A}^2\mathbf{r}_0,\dots,\mathbf{A}^{k-1}\mathbf{r}_0\}
$$
is the **Krylov subspace**. Every vector in it is $p(\mathbf{A})\mathbf{r}_0$ for a polynomial $p$ of degree $<k$, so the error after $k$ steps has the form $\mathbf{e}_k=q_k(\mathbf{A})\mathbf{e}_0$ with $q_k$ a polynomial of degree $k$ and $q_k(0)=1$. Krylov methods differ in *which* $q_k$ they choose — that is, in what "best" means.

#### Conjugate gradients

For SPD $\mathbf{A}$, solving $\mathbf{A}\mathbf{x}=\mathbf{b}$ is equivalent to minimising the quadratic "energy"
$$
\phi(\mathbf{x})=\tfrac12\mathbf{x}^\mathsf{T}\mathbf{A}\mathbf{x}-\mathbf{b}^\mathsf{T}\mathbf{x},\qquad \nabla\phi=\mathbf{A}\mathbf{x}-\mathbf{b}=-\mathbf{r}.
$$
(For a resistive network, $\phi$ is half the dissipated power minus the power delivered by the sources, and the solution minimises it — Maxwell's minimum-heat theorem.) Steepest descent searches along $\mathbf{r}$ each time, which zigzags when the level sets are elongated. CG instead searches along directions $\mathbf{p}_0,\mathbf{p}_1,\dots$ that are **$\mathbf{A}$-conjugate**, $\mathbf{p}_i^\mathsf{T}\mathbf{A}\mathbf{p}_j=0$ for $i\ne j$, so that minimising along a new direction never spoils the minimisation along the old ones.

> [!definition] Conjugate gradient algorithm
> Given SPD $\mathbf{A}$, $\mathbf{b}$, $\mathbf{x}_0$: set $\mathbf{r}_0=\mathbf{b}-\mathbf{A}\mathbf{x}_0$, $\mathbf{p}_0=\mathbf{r}_0$. For $k=0,1,2,\dots$:
> $$
> \begin{aligned}
> \alpha_k&=\frac{\mathbf{r}_k^\mathsf{T}\mathbf{r}_k}{\mathbf{p}_k^\mathsf{T}\mathbf{A}\mathbf{p}_k}, &
> \mathbf{x}_{k+1}&=\mathbf{x}_k+\alpha_k\mathbf{p}_k, &
> \mathbf{r}_{k+1}&=\mathbf{r}_k-\alpha_k\mathbf{A}\mathbf{p}_k,\\
> \beta_k&=\frac{\mathbf{r}_{k+1}^\mathsf{T}\mathbf{r}_{k+1}}{\mathbf{r}_k^\mathsf{T}\mathbf{r}_k}, &
> \mathbf{p}_{k+1}&=\mathbf{r}_{k+1}+\beta_k\mathbf{p}_k .
> \end{aligned}
> $$
> Each iteration costs one matrix–vector product, two inner products and three vector updates — about $2\operatorname{nnz}+10n$ flops — and stores four vectors.

> [!theorem] Properties of CG
> In exact arithmetic:
> 1. The residuals are mutually orthogonal and the directions mutually $\mathbf{A}$-conjugate.
> 2. $\mathbf{x}_k$ minimises the energy norm of the error $\lVert\mathbf{e}\rVert_\mathbf{A}=\sqrt{\mathbf{e}^\mathsf{T}\mathbf{A}\mathbf{e}}$ over $\mathbf{x}_0+\mathcal{K}_k$; equivalently
> $$
> \lVert\mathbf{e}_k\rVert_\mathbf{A}=\min_{q\in\mathcal{P}_k,\ q(0)=1}\lVert q(\mathbf{A})\mathbf{e}_0\rVert_\mathbf{A}\le\min_{q}\max_{\lambda\in\sigma(\mathbf{A})}\lvert q(\lambda)\rvert\ \lVert\mathbf{e}_0\rVert_\mathbf{A}.
> $$
> 3. CG terminates with the exact solution in at most $n$ steps — in fact in at most $d$ steps, where $d$ is the number of distinct eigenvalues of $\mathbf{A}$.
> 4. With $\kappa=\lambda_\text{max}/\lambda_\text{min}$,
> $$
> \frac{\lVert\mathbf{e}_k\rVert_\mathbf{A}}{\lVert\mathbf{e}_0\rVert_\mathbf{A}}\le2\left(\frac{\sqrt\kappa-1}{\sqrt\kappa+1}\right)^k .
> $$

*Sketch.* Property 2 follows because $\phi(\mathbf{x})-\phi(\mathbf{x}^*)=\tfrac12\lVert\mathbf{x}-\mathbf{x}^*\rVert_\mathbf{A}^2$ and conjugacy makes the minimisation over the growing Krylov space separable, one direction at a time. Expanding $\mathbf{e}_0$ in eigenvectors, $\lVert q(\mathbf{A})\mathbf{e}_0\rVert_\mathbf{A}^2=\sum_i\lambda_iq(\lambda_i)^2c_i^2\le\max_i q(\lambda_i)^2\,\lVert\mathbf{e}_0\rVert_\mathbf{A}^2$. Property 3: the polynomial $q(\lambda)=\prod_{i=1}^d(1-\lambda/\lambda_i)$ vanishes on the spectrum. Property 4: choose $q$ to be the Chebyshev polynomial on $[\lambda_\text{min},\lambda_\text{max}]$ scaled so $q(0)=1$; its maximum on the interval is $1/T_k\big(\tfrac{\kappa+1}{\kappa-1}\big)\le2\big(\tfrac{\sqrt\kappa-1}{\sqrt\kappa+1}\big)^k$. $\square$

The bound gives an iteration count of about $k\approx\tfrac12\sqrt\kappa\ln(2/\epsilon)$, compared with $\tfrac12\kappa\ln(1/\epsilon)$ for steepest descent. For the 2-D model problem $\kappa\propto h^{-2}$, so CG needs $O(h^{-1})=O(\sqrt n)$ iterations — the same as optimal SOR, with no parameter to tune.

> [!example] Worked Example 24.19 — Conjugate gradients by hand
> Solve $\begin{bmatrix}4&1\\1&3\end{bmatrix}\mathbf{x}=\begin{bmatrix}1\\2\end{bmatrix}$ by CG from $\mathbf{x}_0=\mathbf{0}$.
>
> **Step 0.** $\mathbf{r}_0=\mathbf{p}_0=(1,2)$, $\mathbf{r}_0^\mathsf{T}\mathbf{r}_0=5$, $\mathbf{A}\mathbf{p}_0=(6,7)$, $\mathbf{p}_0^\mathsf{T}\mathbf{A}\mathbf{p}_0=6+14=20$, $\alpha_0=0.25$.
> $\mathbf{x}_1=(0.25,\ 0.5)$, $\mathbf{r}_1=(1,2)-0.25(6,7)=(-0.5,\ 0.25)$, $\mathbf{r}_1^\mathsf{T}\mathbf{r}_1=0.3125$, $\beta_0=0.3125/5=0.0625$, $\mathbf{p}_1=(-0.5,0.25)+0.0625(1,2)=(-0.4375,\ 0.375)$.
>
> **Step 1.** $\mathbf{A}\mathbf{p}_1=(-1.75+0.375,\ -0.4375+1.125)=(-1.375,\ 0.6875)$, $\mathbf{p}_1^\mathsf{T}\mathbf{A}\mathbf{p}_1=0.6016+0.2578=0.859375$, $\alpha_1=0.3125/0.859375=0.363636=\tfrac4{11}$.
> $\mathbf{x}_2=(0.25,0.5)+\tfrac4{11}(-0.4375,0.375)=(0.090909,\ 0.636364)=(\tfrac1{11},\tfrac7{11})$, and $\mathbf{r}_2=\mathbf{0}$.
>
> **Checks.** $4/11+7/11=1$ ✓, $1/11+21/11=2$ ✓. Termination in $n=2$ steps ✓. Orthogonality: $\mathbf{r}_0^\mathsf{T}\mathbf{r}_1=-0.5+0.5=0$ ✓; conjugacy: $\mathbf{p}_0^\mathsf{T}\mathbf{A}\mathbf{p}_1=(1,2)\cdot(-1.375,0.6875)=0$ ✓.

> [!example] Worked Example 24.20 — Iteration counts from the condition number
> An SPD system from a capacitance extraction has $\kappa=10^4$. Estimate the iterations for a $10^{-6}$ reduction of the energy-norm error with steepest descent and with CG.
>
> **Steepest descent.** Error factor $\le\big(\tfrac{\kappa-1}{\kappa+1}\big)^k\approx e^{-2k/\kappa}$: $k\approx\tfrac{\kappa}2\ln10^6=5000\times13.8\approx69\,000$.
>
> **CG.** $k\approx\tfrac{\sqrt\kappa}{2}\ln(2\times10^6)=50\times14.5\approx725$.
>
> **Reading.** CG's dependence on $\sqrt\kappa$ rather than $\kappa$ is worth a factor of 100 here. And the bound is pessimistic: if the eigenvalues are clustered (a few outliers and the rest close together), CG "removes" each outlier in about one iteration and then converges at the rate set by the cluster. That is exactly the effect a good preconditioner aims for.

#### GMRES

For nonsymmetric or indefinite matrices — modified nodal analysis with voltage sources, convection–diffusion models, integral-equation formulations of electromagnetic scattering — there is no energy to minimise. **GMRES** (generalised minimal residual) chooses $\mathbf{x}_k\in\mathbf{x}_0+\mathcal{K}_k$ to minimise the residual norm $\lVert\mathbf{b}-\mathbf{A}\mathbf{x}_k\rVert_2$.

> [!definition] Arnoldi process and GMRES
> The **Arnoldi process** builds an orthonormal basis $\mathbf{q}_1,\dots,\mathbf{q}_k$ of $\mathcal{K}_k$ by modified Gram–Schmidt: $\mathbf{q}_1=\mathbf{r}_0/\beta$ with $\beta=\lVert\mathbf{r}_0\rVert$, and for each $j$, $\mathbf{w}=\mathbf{A}\mathbf{q}_j$ is orthogonalised against $\mathbf{q}_1,\dots,\mathbf{q}_j$ (coefficients $h_{ij}$) and normalised to give $\mathbf{q}_{j+1}$ ($h_{j+1,j}=\lVert\mathbf{w}\rVert$). In matrix form $\mathbf{A}\mathbf{Q}_k=\mathbf{Q}_{k+1}\tilde{\mathbf{H}}_k$ with $\tilde{\mathbf{H}}_k$ a $(k+1)\times k$ upper Hessenberg matrix. Writing $\mathbf{x}_k=\mathbf{x}_0+\mathbf{Q}_k\mathbf{y}$,
> $$
> \lVert\mathbf{b}-\mathbf{A}\mathbf{x}_k\rVert_2=\lVert\beta\mathbf{e}_1-\tilde{\mathbf{H}}_k\mathbf{y}\rVert_2,
> $$
> a small $(k+1)\times k$ least-squares problem solved by Givens rotations. **GMRES** returns the minimising $\mathbf{y}$.

> [!theorem] GMRES convergence
> GMRES residuals are non-increasing and satisfy
> $$
> \frac{\lVert\mathbf{r}_k\rVert}{\lVert\mathbf{r}_0\rVert}=\min_{q\in\mathcal{P}_k,\ q(0)=1}\frac{\lVert q(\mathbf{A})\mathbf{r}_0\rVert}{\lVert\mathbf{r}_0\rVert}\le\kappa(\mathbf{V})\min_{q}\max_{\lambda\in\sigma(\mathbf{A})}\lvert q(\lambda)\rvert
> $$
> for diagonalisable $\mathbf{A}=\mathbf{V}\boldsymbol\Lambda\mathbf{V}^{-1}$. GMRES terminates in at most $d$ steps if $\mathbf{A}$ has a minimal polynomial of degree $d$. Convergence is fast when the eigenvalues are clustered away from the origin and $\mathbf{V}$ is well conditioned.

The price of GMRES is that every step must orthogonalise against *all* previous basis vectors: step $k$ costs $O(kn)$ flops and the method stores $k$ vectors. In practice it is **restarted** every $m$ steps (GMRES($m$), with $m$ typically 20–100), which bounds the memory but can slow or even stall convergence. Short-recurrence alternatives for nonsymmetric systems — BiCGSTAB, QMR, IDR($s$) — have fixed cost per step but less predictable convergence.

> [!example] Worked Example 24.21 — GMRES and the minimal polynomial
> A $50\times50$ nonsymmetric matrix $\mathbf{A}=\mathbf{S}\boldsymbol\Lambda\mathbf{S}^{-1}$ is constructed with random $\mathbf{S}$ and eigenvalues $1$ (25 times) and $3$ (25 times). Predict, and then observe, the GMRES iteration count.
>
> **Prediction.** $(\mathbf{A}-\mathbf{I})(\mathbf{A}-3\mathbf{I})=\mathbf{S}(\boldsymbol\Lambda-\mathbf{I})(\boldsymbol\Lambda-3\mathbf{I})\mathbf{S}^{-1}=\mathbf{0}$, so the polynomial $q(\lambda)=(1-\lambda)(1-\lambda/3)$, with $q(0)=1$, annihilates $\mathbf{A}$. GMRES must converge in 2 steps regardless of $n$.
>
> **Observation.** SciPy's `gmres` reports relative residuals $0.988$ after one step and $1.6\times10^{-13}$ after two ✓.
>
> **Contrast.** For a $1024$-unknown convection–diffusion matrix (32×32 grid) with scattered eigenvalues, full GMRES needed 74 iterations for a $10^{-8}$ residual reduction while GMRES(20) needed 229 — restarting throws away the subspace that was making progress. This is why preconditioning, which clusters the spectrum, matters even more for GMRES than for CG.

### Preconditioning

The convergence of every Krylov method depends on the spectrum of $\mathbf{A}$. **Preconditioning** replaces the system by an equivalent one with a better spectrum:
$$
\text{left: }\ \mathbf{P}^{-1}\mathbf{A}\mathbf{x}=\mathbf{P}^{-1}\mathbf{b},\qquad\text{right: }\ \mathbf{A}\mathbf{P}^{-1}\mathbf{y}=\mathbf{b},\ \ \mathbf{x}=\mathbf{P}^{-1}\mathbf{y},
$$
where the preconditioner $\mathbf{P}$ satisfies two competing requirements: $\mathbf{P}\approx\mathbf{A}$ (so $\mathbf{P}^{-1}\mathbf{A}\approx\mathbf{I}$, eigenvalues clustered near 1) and $\mathbf{P}\mathbf{z}=\mathbf{r}$ cheap to solve. $\mathbf{P}=\mathbf{I}$ is cheap but useless; $\mathbf{P}=\mathbf{A}$ is perfect but as hard as the original problem. For CG, $\mathbf{P}$ must be SPD, and the preconditioned algorithm (PCG) is CG in the $\mathbf{P}$-inner product: the only change is that each iteration solves $\mathbf{P}\mathbf{z}_k=\mathbf{r}_k$ and uses $\mathbf{z}_k$ in place of $\mathbf{r}_k$ in the direction update and in $\alpha_k,\beta_k$ (via $\mathbf{r}_k^\mathsf{T}\mathbf{z}_k$).

Common preconditioners, in rough order of cost and power:

| Preconditioner | $\mathbf{P}$ | Comment |
|---|---|---|
| Jacobi (diagonal) | $\mathbf{D}$ | removes bad *scaling*; free; weak otherwise |
| SSOR | $(\mathbf{D}+\omega\mathbf{L}_A)\mathbf{D}^{-1}(\mathbf{D}+\omega\mathbf{U}_A)$ | symmetric version of SOR; no setup |
| Incomplete Cholesky IC(0) / incomplete LU ILU(0) | $\tilde{\mathbf{L}}\tilde{\mathbf{U}}$ with the sparsity of $\mathbf{A}$ | elimination that discards all fill; workhorse of circuit and FEM codes |
| ILUT($\tau$, $p$) | incomplete LU keeping fill above threshold $\tau$, at most $p$ per row | tunable between ILU(0) and full LU |
| Multigrid / algebraic multigrid (AMG) | one V-cycle | $O(n)$ work and iteration count independent of $h$ for elliptic problems |
| Domain decomposition | local solves on subdomains plus a coarse problem | parallel |

> [!example] Worked Example 24.22 — Diagonal scaling as a preconditioner
> $\mathbf{A}=\begin{bmatrix}10^4&1\\1&1\end{bmatrix}$ (a stiff spring coupled to a soft one, or a $10\ \mathrm{k\Omega}$ and a $1\ \Omega$ branch sharing a node). Compare $\kappa_2(\mathbf{A})$ with that of the symmetrically Jacobi-preconditioned matrix $\mathbf{D}^{-1/2}\mathbf{A}\mathbf{D}^{-1/2}$.
>
> **Unpreconditioned.** Eigenvalues $\approx10^4$ and $\approx0.9999$: $\kappa_2\approx1.0\times10^4$.
>
> **Scaled.** $\mathbf{D}^{-1/2}=\operatorname{diag}(10^{-2},1)$ gives $\begin{bmatrix}1&0.01\\0.01&1\end{bmatrix}$, eigenvalues $1\pm0.01$, $\kappa_2=1.01/0.99=1.020$.
>
> **Reading.** CG on the scaled system converges in 2 steps anyway ($n=2$), but the iteration-count bound drops from about $\tfrac12\sqrt{10^4}\ln(2/\epsilon)$ to essentially 1 step per digit. For large systems diagonal scaling never hurts and fixes purely unit-induced ill-conditioning; it cannot fix the *geometric* ill-conditioning of the model problem, whose diagonal is constant.

> [!example] Worked Example 24.23 — Incomplete Cholesky on the 2-D Poisson problem
> The 5-point Laplacian on a $k\times k$ grid ($n=k^2$ unknowns) is solved by CG to a relative residual of $10^{-8}$, with and without an IC(0) preconditioner. Observed iteration counts:
>
> | grid | $n$ | $\kappa(\mathbf{A})$ | CG | PCG with IC(0) |
> |---|---|---|---|---|
> | $16\times16$ | 256 | 116 | 28 | 17 |
> | $32\times32$ | 1024 | 441 | 59 | 29 |
> | $64\times64$ | 4096 | 1712 | 119 | — |
>
> **Reading.** Without preconditioning the count doubles each time $k$ doubles, tracking $\sqrt\kappa\propto k$ ($\sqrt{441}=21$, $\sqrt{1712}=41$); the bound $\tfrac12\sqrt\kappa\ln(2/\epsilon)\approx200$ for $k=32$ is pessimistic by a factor of 3, as usual. IC(0) halves the count but does not change its growth rate (its preconditioned condition number still grows like $h^{-2}$, with a smaller constant). Each PCG iteration costs about twice a CG iteration (two triangular solves), so the net gain is modest here; it is much larger for problems with variable coefficients (materials with different permittivities or conductivities), where IC also corrects the scaling. Only multigrid makes the count independent of $k$.

> [!warning] Common Pitfall — wrong method for the matrix
> CG requires $\mathbf{A}$ (and the preconditioner) to be symmetric positive definite. Applied to a nonsymmetric MNA matrix, or with a nonsymmetric ILU preconditioner, it may stall or diverge without any error message — in a test with an ILU preconditioner on the Poisson problem SciPy's `cg` ran to its iteration limit and returned garbage, while the same preconditioner inside GMRES worked. Check symmetry ($\lVert\mathbf{A}-\mathbf{A}^\mathsf{T}\rVert$) before choosing CG, use symmetric preconditioners (Jacobi, SSOR, IC) with it, and always check the returned convergence flag (`info` in SciPy).

> [!warning] Common Pitfall — stopping on the wrong quantity
> A relative residual of $10^{-6}$ does not mean six correct digits: the error bound is $\kappa$ times larger. Conversely, demanding $10^{-14}$ in a problem whose data are known to $10^{-3}$ wastes iterations. Choose the tolerance from the accuracy actually needed in the *output quantity* — often a single node voltage or a capacitance — and, if possible, monitor that quantity directly.

> [!tip] Engineering Note
> In IC design, power-grid analysis (IR drop) solves SPD systems with $10^8$ unknowns by PCG with algebraic multigrid or incomplete Cholesky preconditioners; full-wave EM solvers using integral equations (method of moments) solve dense complex systems by GMRES with matrix–vector products accelerated by the fast multipole method; transient circuit simulation, where each Newton step needs only a modest improvement on the previous solution, uses the previous LU factors as a preconditioner for a few GMRES steps. The common thread: the matrix–vector product and the preconditioner are designed together with the physics.

```python
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla
k = 64
T = sp.diags([-1.0, 2.0, -1.0], [-1, 0, 1], shape=(k, k))
A = (sp.kron(sp.eye(k), T) + sp.kron(T, sp.eye(k))).tocsr()   # 2-D Laplacian
b = np.ones(k * k)
its = []
x, info = spla.cg(A, b, rtol=1e-8, callback=lambda xk: its.append(1))
print(info, len(its))          # 0 (converged), 119 iterations
```

## 24.5 Sparse computation

### Motivation

The matrices of engineering are overwhelmingly **sparse**: a node in a circuit connects to a handful of others, a finite-element basis function overlaps only its neighbours, a bus in a power grid connects to a few lines. [Chapter 8](../part-3-linear-algebra/08-matrices-and-linear-systems.md) introduced sparse storage (COO and CSR), fill-in and the effect of elimination order through the star-network example, and the Thomas algorithm for tridiagonal systems. This section makes those ideas quantitative: how the storage formats support the kernels that iterative and direct solvers need, how fill-in is predicted from the graph before any arithmetic is done, how reordering algorithms work, and how the cost of sparse factorisation scales with the geometry of the problem.

### Storage formats and their kernels

> [!definition] CSR and CSC
> For an $m\times n$ matrix with $\operatorname{nnz}$ nonzeros, **compressed sparse row** (CSR) stores three arrays: `data` (the nonzero values, row by row), `indices` (the column index of each value) and `indptr` (length $m+1$; row $i$ occupies positions `indptr[i]` to `indptr[i+1]-1`). **Compressed sparse column** (CSC) stores the same information column by column: `indices` then holds row indices and `indptr` has length $n+1$. Both need $2\operatorname{nnz}+(\text{dimension}+1)$ words.

The choice between them is dictated by the access pattern of the kernel:

| Kernel | Natural format | Why |
|---|---|---|
| $\mathbf{y}=\mathbf{A}\mathbf{x}$ | CSR | each $y_i$ is an inner product of row $i$ with $\mathbf{x}$ |
| $\mathbf{y}=\mathbf{A}^\mathsf{T}\mathbf{x}$ | CSC (or CSR of $\mathbf{A}$ read as CSC of $\mathbf{A}^\mathsf{T}$) | column access |
| left-looking LU, Cholesky, column pivoting | CSC | algorithms proceed column by column |
| row slicing, row-wise assembly | CSR | |
| assembly by stamping | COO | append (row, col, value) triplets; duplicates summed on conversion |
| ladder/1-D problems | banded / diagonal (DIA) | contiguous diagonals |

The CSR matrix–vector product is the most executed loop in scientific computing:

```python
def csr_matvec(data, indices, indptr, x):
    y = [0.0] * (len(indptr) - 1)
    for i in range(len(y)):
        s = 0.0
        for p in range(indptr[i], indptr[i + 1]):
            s += data[p] * x[indices[p]]      # indirect access to x
        y[i] = s
    return y
```

It performs $2\operatorname{nnz}$ flops and reads $\operatorname{nnz}$ values, $\operatorname{nnz}$ indices and (irregularly) entries of $\mathbf{x}$: an arithmetic intensity well below one flop per word, so its speed is set by memory bandwidth (Section 24.1). Orderings that cluster nonzeros near the diagonal improve cache reuse of $\mathbf{x}$ — a second, purely performance-related reason to reorder.

> [!example] Worked Example 24.24 — CSR and CSC of a nonsymmetric matrix, and a product
> Store
> $$
> \mathbf{M}=\begin{bmatrix}4&0&-1&0\\1&3&0&-2\\0&0&5&0\\0&-2&0&6\end{bmatrix}
> $$
> in CSR and CSC (0-based), and compute $\mathbf{M}\mathbf{x}$ for $\mathbf{x}=(1,2,3,4)^\mathsf{T}$ from the CSR arrays.
>
> **CSR.** `data` $=[4,-1,\ 1,3,-2,\ 5,\ -2,6]$, `indices` $=[0,2,\ 0,1,3,\ 2,\ 1,3]$, `indptr` $=[0,2,5,6,8]$.
>
> **CSC.** Reading down the columns: column 0 has $4$ (row 0), $1$ (row 1); column 1 has $3$ (row 1), $-2$ (row 3); column 2 has $-1$ (row 0), $5$ (row 2); column 3 has $-2$ (row 1), $6$ (row 3). So `data` $=[4,1,\ 3,-2,\ -1,5,\ -2,6]$, `indices` $=[0,1,\ 1,3,\ 0,2,\ 1,3]$, `indptr` $=[0,2,4,6,8]$.
>
> **Product.** Row 0: positions 0–1, $4\cdot x_0+(-1)\cdot x_2=4-3=1$. Row 1: positions 2–4, $1\cdot1+3\cdot2+(-2)\cdot4=-1$. Row 2: position 5, $5\cdot3=15$. Row 3: positions 6–7, $-2\cdot2+6\cdot4=20$. So $\mathbf{M}\mathbf{x}=(1,-1,15,20)^\mathsf{T}$ — 16 flops for 8 nonzeros.
>
> **Check.** A dense product gives the same result with 28 flops (16 multiplications, 12 additions), half of them on zeros. For a $10^6$-node circuit with 5 nonzeros per row the ratio is $10^6/5=200\,000$.

### The graph of a sparse matrix and the elimination game

> [!definition] Adjacency graph and fill
> The **graph** $G(\mathbf{A})$ of a structurally symmetric $n\times n$ matrix has vertices $1,\dots,n$ and an edge $\{i,j\}$ whenever $a_{ij}\ne0$, $i\ne j$. Eliminating vertex $k$ (one step of Gaussian elimination) removes $k$ and adds an edge between every pair of its remaining neighbours that were not already adjacent; these new edges are the **fill-in**. The **filled graph** $G^+(\mathbf{A})$ is $G(\mathbf{A})$ plus all fill edges; its edges are exactly the nonzeros of $\mathbf{L}+\mathbf{L}^\mathsf{T}$ (barring numerical cancellation).

This purely combinatorial description has three consequences that sparse direct solvers exploit.

1. **Symbolic factorisation.** The pattern of $\mathbf{L}$ can be computed from the graph alone, before any floating-point work, in time roughly proportional to $\operatorname{nnz}(\mathbf{L})$. Memory is allocated once.
2. **Reuse.** In Newton iterations, frequency sweeps, transient time steps and Monte Carlo runs, the *values* of the matrix change but its *pattern* does not. Ordering and symbolic factorisation are done once; only the numeric factorisation is repeated. SPICE has worked this way since the 1970s.
3. **The elimination tree.** Define the parent of column $j$ as the row index of the first off-diagonal nonzero in column $j$ of $\mathbf{L}$. The resulting tree encodes all dependencies: columns in disjoint subtrees can be factored independently, in parallel. Groups of consecutive columns with identical pattern (**supernodes**) are processed as dense blocks with level-3 BLAS, recovering the cache efficiency of dense code.

### Reordering

The amount of fill depends on the elimination order, i.e. on the symmetric permutation $\mathbf{P}\mathbf{A}\mathbf{P}^\mathsf{T}$. Finding the order that minimises fill is NP-hard, so heuristics are used.

- **Bandwidth and profile reduction: (reverse) Cuthill–McKee.** Number the vertices by breadth-first search from a peripheral vertex, visiting neighbours in order of increasing degree; reversing the order (RCM) usually reduces the profile further. Elimination on a banded matrix creates fill only inside the band, so the fill is at most $n\times$(bandwidth).
- **Minimum degree** (and its fast approximation AMD): at each step eliminate the vertex with fewest neighbours in the current graph, since eliminating a degree-$d$ vertex creates at most $d(d-1)/2$ fill edges. This is the greedy rule of the star example in [Ch. 8](../part-3-linear-algebra/08-matrices-and-linear-systems.md) — leaves before hubs.
- **Nested dissection.** Find a small vertex **separator** $S$ that splits the graph into two halves $A$ and $B$ with no edges between them; number $A$, then $B$, then $S$ last; recurse on $A$ and $B$. Since $A$ and $B$ are not connected, no fill ever links them; fill is confined to the separator blocks.
- **Column orderings for unsymmetric LU** (COLAMD) order the columns to limit fill in the factors of $\mathbf{A}^\mathsf{T}\mathbf{A}$'s pattern, leaving row pivoting free for numerical stability.

> [!theorem] Fill and work for model grids
> For the 5-point (2-D) and 7-point (3-D) Laplacians on grids with $n$ unknowns, Cholesky factorisation costs, to leading order in $n$:
>
> | Ordering | 2-D: $\operatorname{nnz}(\mathbf{L})$ | 2-D: flops | 3-D: $\operatorname{nnz}(\mathbf{L})$ | 3-D: flops |
> |---|---|---|---|---|
> | natural (banded) | $n^{3/2}$ | $n^2$ | $n^{5/3}$ | $n^{7/3}$ |
> | nested dissection | $n\log n$ | $n^{3/2}$ | $n^{4/3}$ | $n^2$ |
>
> Nested dissection is asymptotically optimal for these grids (George, 1973; Lipton–Rose–Tarjan).

*Sketch, 2-D.* A $k\times k$ grid ($n=k^2$) has a separator of $k$ vertices (a grid line). After the two halves are eliminated the separator block is dense, $k\times k$, costing $O(k^3)=O(n^{3/2})$ flops to factor; the recursion $W(n)=4W(n/4)+O(n^{3/2})$ (splitting in both directions) gives $W=O(n^{3/2})$ by the master theorem (Section 24.1, $a=4$, $b=4$, $d=3/2>1$). The separators at level $\ell$ contribute $4^\ell$ dense blocks of size $(k/2^\ell)^2$, i.e. $O(n)$ fill per level over $\log n$ levels: $O(n\log n)$. With natural ordering the bandwidth is $k=\sqrt n$, the band fills completely ($n\cdot k=n^{3/2}$ entries) and the banded Cholesky costs $n k^2=n^2$. $\square$

In 3-D the separators are planes of $n^{2/3}$ vertices, and their dense blocks cost $(n^{2/3})^3=n^2$ — this is why 3-D field problems are the domain of iterative methods while 2-D problems and circuits (whose graphs have small separators) are routinely solved directly.

> [!example] Worked Example 24.25 — Reverse Cuthill–McKee on a badly numbered chain
> Six components form a chain whose nodes, in physical order along the chain, were numbered $0,3,5,1,4,2$ by a netlister. The nodal matrix has 2 on the diagonal and $-1$ for each link. Find its bandwidth, then renumber by RCM.
>
> **Original.** The links are $\{0,3\},\{3,5\},\{5,1\},\{1,4\},\{4,2\}$, so the nonzero off-diagonals are in positions with $\lvert i-j\rvert=3,2,4,3,2$. The bandwidth is $4$; banded elimination would treat the matrix as having up to 9 diagonals.
>
> **RCM.** Start from a peripheral vertex (degree 1): vertex 2 or vertex 0. BFS from 0 visits $0,3,5,1,4,2$ — the physical order. Reversed: $2,4,1,5,3,0$. (SciPy's `reverse_cuthill_mckee` returns exactly this permutation.)
>
> **Result.** In the new numbering consecutive chain nodes have consecutive numbers: the matrix is tridiagonal, bandwidth 1, and can be solved by the Thomas algorithm in $8n$ flops with zero fill.
>
> **Reading.** For a chain this is obvious by eye; for a $10^5$-node mesh it is not, and RCM does it in $O(\operatorname{nnz})$ time.

> [!example] Worked Example 24.26 — Measured fill on a 2-D grid
> The 2-D Laplacian on $k\times k$ grids is factored by SciPy's SuperLU with three column orderings. Nonzeros in $\mathbf{L}+\mathbf{U}$:
>
> | $k$ | $n$ | $\operatorname{nnz}(\mathbf{A})$ | natural | COLAMD | minimum degree on $\mathbf{A}+\mathbf{A}^\mathsf{T}$ |
> |---|---|---|---|---|---|
> | 30 | 900 | 4 380 | 54 058 | 30 338 | 20 196 |
> | 60 | 3 600 | 17 760 | 432 118 | 179 720 | 108 430 |
> | 100 | 10 000 | 49 600 | 2 000 198 | 645 750 | 371 346 |
>
> **Check against theory.** Natural ordering fills the band: $\operatorname{nnz}(\mathbf{L}+\mathbf{U})\approx2nk=2k^3$, i.e. $54\,000$, $432\,000$, $2\,000\,000$ ✓ — growing as $n^{3/2}$ (factor 8 per doubling of $k$). Minimum degree grows by a factor of $5.4$ from $k=30$ to $k=60$, close to the $n\log n$ prediction $4\times\log(3600)/\log(900)\approx4.8$.
>
> **Consequence.** At $k=100$ the minimum-degree factors need $5.4$ times less memory than natural ordering, and the flop savings are larger still (the flop count scales as the sum of squared column counts). At $k=1000$ ($n=10^6$) natural ordering would need $2\times10^9$ nonzeros (16 GB) and $\sim10^{12}$ flops, nested dissection about $\tfrac{31}{4}n\log_2n\approx1.5\times10^{8}$ nonzeros (1.2 GB) and $\tfrac{829}{84}n^{3/2}\approx10^{10}$ flops (George's constants for the model grid).

> [!tip] Engineering Note — the circuit matrix is special
> Circuit matrices differ from FEM matrices in ways that matter for solver choice: they are very sparse (often under 5 nonzeros per row), have a few dense rows and columns (supply rails, clock nets), are unsymmetric when controlled sources and voltage sources are present, and have small separators (circuits are "almost planar" locally). The KLU solver, designed for them, first permutes to **block triangular form** (each strongly connected component of the graph becomes a diagonal block, so only those blocks need factoring), then orders each block with AMD and factors it with partial pivoting. Power-grid analysis matrices, by contrast, are SPD and mesh-like, and are solved with sparse Cholesky (for $n\lesssim10^7$) or PCG with multigrid.

> [!warning] Common Pitfall — accidental densification
> Several innocent operations silently produce dense matrices: `A.toarray()` or `A.todense()` "for debugging"; `scipy.sparse.linalg.inv(A)` (the inverse of a sparse matrix is almost always dense — a node voltage depends on every source); `A + c` with a scalar $c\ne0$; `numpy.linalg.solve(A.toarray(), b)` when a sparse solve was intended; or building a sparse matrix by element insertion into CSR (`A[i, j] = v` in a loop, quadratic cost — assemble in COO and convert once). Also beware of explicit zeros: a stamped matrix may store entries that cancel to exactly zero; call `eliminate_zeros()` only if the pattern will not change later (otherwise the stored zero keeps the symbolic factorisation valid).

> [!warning] Common Pitfall — ordering versus pivoting
> For SPD matrices any symmetric ordering is numerically stable, so the ordering can be chosen purely for sparsity. For unsymmetric or indefinite matrices (MNA with voltage sources has zeros on the diagonal), the pivot sequence must also guard against small pivots. Sparse LU uses **threshold pivoting**: accept the sparsity-preferred pivot $a_{kk}$ if $\lvert a_{kk}\rvert\ge u\max_i\lvert a_{ik}\rvert$ with $u$ around $0.1$ (or $0.001$ in circuit simulators), otherwise pivot. Setting $u=0$ for speed can give catastrophic growth.

### Direct or iterative?

Putting Sections 24.4 and 24.5 together gives a practical rule:

| Situation | Usual choice |
|---|---|
| circuits, 1-D and 2-D problems up to $\sim10^6$–$10^7$ unknowns; many right-hand sides; need robustness | sparse direct (LU, Cholesky) with fill-reducing ordering |
| 3-D field problems, $>10^6$ unknowns | preconditioned Krylov (PCG for SPD, GMRES/BiCGSTAB otherwise), multigrid |
| matrix not available explicitly (only $\mathbf{A}\mathbf{v}$, e.g. fast multipole) | Krylov |
| sequence of slowly varying systems | direct factors of one reused as preconditioner for the others |
| low accuracy acceptable, cheap iterations | Krylov with loose tolerance |

## 24.6 Monte Carlo methods

### Motivation

Some quantities are expectations that are hard to compute any other way: the probability that a filter built from 5% components meets its specification, the bit-error rate of a receiver with a nonlinear equaliser, the average power absorbed by a randomly oriented antenna, an integral over a ten-dimensional parameter space. **Monte Carlo methods** estimate such quantities by *simulation*: draw random samples of the inputs, evaluate the quantity for each, and average. They are slow to converge — the error falls only as $1/\sqrt N$ — but that rate does not depend on the dimension of the problem or the smoothness of the function, and the method needs nothing more than the ability to evaluate the model. For high-dimensional problems, discontinuous pass/fail criteria, and black-box simulators, Monte Carlo is often the only practical method. Random variables, expectation, variance and the generation of samples from a given distribution (inverse transform) were developed in [Ch. 19](../part-6-probability-statistics/19-probability.md); sampling distributions and confidence intervals in [Ch. 20](../part-6-probability-statistics/20-statistical-analysis.md).

### The basic estimator

> [!definition] Monte Carlo estimator
> To estimate $\mu=E[f(\mathbf{X})]$, where $\mathbf{X}$ is a random vector with known distribution, draw independent samples $\mathbf{X}_1,\dots,\mathbf{X}_N$ and form
> $$
> \hat\mu_N=\frac1N\sum_{i=1}^Nf(\mathbf{X}_i),\qquad s_N^2=\frac1{N-1}\sum_{i=1}^N\big(f(\mathbf{X}_i)-\hat\mu_N\big)^2 .
> $$
> An integral $\int_\Omega g(\mathbf{x})\,d\mathbf{x}$ over a region of volume $V$ is the special case $\mu=V\,E[g(\mathbf{X})]$ with $\mathbf{X}$ uniform on $\Omega$. A probability $P(A)$ is the special case $f=\mathbf{1}_A$, the indicator of the event.

> [!theorem] Convergence of Monte Carlo
> If $\sigma^2=\operatorname{Var}(f(\mathbf{X}))<\infty$, then
> 1. $\hat\mu_N$ is unbiased, $E[\hat\mu_N]=\mu$, with $\operatorname{Var}(\hat\mu_N)=\sigma^2/N$;
> 2. $\hat\mu_N\to\mu$ almost surely (strong law of large numbers);
> 3. $\sqrt N(\hat\mu_N-\mu)/\sigma\to\mathcal{N}(0,1)$ in distribution (central limit theorem), so for large $N$
> $$
> \mu\in\Big[\hat\mu_N-z\frac{s_N}{\sqrt N},\ \hat\mu_N+z\frac{s_N}{\sqrt N}\Big]
> $$
> with confidence about 95% for $z=1.96$ and 99.7% for $z=3$.
>
> The quantity $s_N/\sqrt N$ is the **standard error** of the estimate.

*Proof of 1.* Linearity of expectation gives $E[\hat\mu_N]=\tfrac1N\sum E[f(\mathbf{X}_i)]=\mu$; independence gives $\operatorname{Var}(\hat\mu_N)=\tfrac1{N^2}\sum\operatorname{Var}(f(\mathbf{X}_i))=\sigma^2/N$. $\square$

Three features stand out. The error decreases as $N^{-1/2}$: **each additional correct digit costs 100 times more samples**. The constant is $\sigma$, the standard deviation of the *integrand* — so reducing $\sigma$ (variance reduction) is as good as increasing $N$. And neither the dimension $d$ of $\mathbf{X}$ nor the smoothness of $f$ appears. Compare a product trapezoidal rule in $d$ dimensions, whose error is $O(h^2)=O(N^{-2/d})$ for $N$ points: it beats Monte Carlo for $d<4$, ties at $d=4$, and loses for $d>4$.

```desmos-graph
left=0; right=8; bottom=-6; top=0.5
---
y=-0.5x|blue
y=-2x|green|dashed
y=-\frac{2}{4}x-0.02|orange|dotted
y=-\frac{2}{10}x|red|dashed
```

*Log–log error versus number of function evaluations ($x=\log_{10}N$, $y=\log_{10}$ error, constants set to 1): Monte Carlo (blue, slope $-\tfrac12$ in any dimension) against the trapezoidal rule in 1 dimension (green, slope $-2$), 4 dimensions (orange, slope $-\tfrac12$, coinciding with Monte Carlo) and 10 dimensions (red, slope $-0.2$). In high dimension Monte Carlo wins.*

> [!example] Worked Example 24.27 — A one-dimensional integral and its error bar
> Estimate $I=\int_0^1e^x\,dx=e-1=1.718\,282$ by Monte Carlo with $N=10^4$ uniform samples, and predict the standard error.
>
> **Variance of the integrand.** With $U\sim\mathcal{U}(0,1)$, $E[e^{2U}]=\tfrac12(e^2-1)=3.19453$, so $\sigma^2=3.19453-(e-1)^2=0.242036$, $\sigma=0.49197$.
>
> **Predicted standard error.** $\sigma/\sqrt N=0.49197/100=0.00492$.
>
> **One run** (NumPy `default_rng(2024)`): $\hat I=1.71993$, $s_N/\sqrt N=0.00493$. The actual error $0.00165$ is $0.33$ standard errors — consistent with the CLT.
>
> **Reading.** Three correct digits from $10^4$ evaluations; the trapezoidal rule gets six digits with 40 points (error $9\times10^{-5}$ with $n=40$, Section 24.9). In one dimension Monte Carlo is a poor integrator; its value lies elsewhere.

### Variance reduction

Since the error is $\sigma/\sqrt N$, any transformation that keeps the mean but reduces the variance of the averaged quantity by a factor $r$ is equivalent to multiplying $N$ by $r$. Four standard techniques:

> [!definition] Variance-reduction techniques
> - **Antithetic variates.** Average $f$ over pairs of negatively correlated samples, e.g. $\tfrac12\big(f(U)+f(1-U)\big)$. The variance of a pair mean is $\tfrac12\sigma^2(1+\rho)$, where $\rho=\operatorname{Corr}(f(U),f(1-U))$; it helps when $f$ is monotone.
> - **Control variates.** If $h(\mathbf{X})$ has *known* mean $\nu$ and is correlated with $f$, estimate $E[f-c(h-\nu)]$. The optimal coefficient and the resulting variance are
> $$
> c^*=\frac{\operatorname{Cov}(f,h)}{\operatorname{Var}(h)},\qquad \operatorname{Var}\big(f-c^*(h-\nu)\big)=\sigma^2(1-\rho_{fh}^2).
> $$
> - **Stratified sampling.** Partition the domain into strata, sample each stratum separately (proportionally to its probability) and combine. The between-strata component of the variance is eliminated.
> - **Importance sampling.** Sample from a different density $q$ that puts more mass where $\lvert f\rvert p$ is large, and reweight: $\mu=E_p[f(\mathbf{X})]=E_q\big[f(\mathbf{X})\,w(\mathbf{X})\big]$ with likelihood ratio $w=p/q$.

*Derivation of the control-variate optimum.* $\operatorname{Var}(f-ch)=\sigma_f^2-2c\operatorname{Cov}(f,h)+c^2\sigma_h^2$, a quadratic in $c$ minimised at $c^*=\operatorname{Cov}(f,h)/\sigma_h^2$, where it equals $\sigma_f^2-\operatorname{Cov}(f,h)^2/\sigma_h^2=\sigma_f^2(1-\rho^2)$. A correlation of $0.99$ reduces the variance 50-fold. $\square$

> [!example] Worked Example 24.28 — Three variance reductions on the same integral
> Apply antithetic variates, a control variate and stratification to $I=\int_0^1e^x\,dx$, at a cost of $10^4$ evaluations of $e^x$ each.
>
> **Antithetic.** $\operatorname{Cov}(e^U,e^{1-U})=E[e]-(e-1)^2=2.71828-2.95249=-0.23421$, so a pair mean has variance $\tfrac12(0.24204-0.23421)=0.003912$. Two independent samples averaged would have variance $0.24204/2=0.12102$: a reduction by a factor $30.9$ for the same work. With 5000 pairs: $\hat I=1.71789\pm0.00088$.
>
> **Control variate $h(U)=U$** (mean $\tfrac12$, variance $\tfrac1{12}$). $\operatorname{Cov}(e^U,U)=E[Ue^U]-\tfrac12(e-1)=1-0.85914=0.14086$, so $c^*=0.14086\times12=1.6903$ and the residual variance is $0.24204-0.14086^2\times12=0.003940$: a factor $61$. With $10^4$ samples: $\hat I=1.718462\pm0.00063$.
>
> **Stratified, one sample in each of $10^4$ equal strata.** $\hat I=1.718\,281\,80$, error $3\times10^{-8}$ — the estimator is now essentially a randomised midpoint rule, whose error for smooth $f$ falls as $N^{-3/2}$.
>
> **Reading.** Variance reduction exploits knowledge of the integrand (monotonicity, a correlated function with known mean, smoothness). The cost is a little analysis; the gain can be orders of magnitude.

### Rare events and importance sampling

Estimating a small probability $p$ with plain Monte Carlo is expensive: the indicator has variance $p(1-p)$, so the *relative* standard error is
$$
\frac{\sqrt{p(1-p)/N}}{p}\approx\frac1{\sqrt{Np}},
$$
and a 10% relative error needs $N\approx100/p$ samples — about 100 *observed* events.

> [!example] Worked Example 24.29 — A tail probability and a bit-error rate
> (a) Estimate $p=P(X>4)$ for $X\sim\mathcal{N}(0,1)$ (exact $p=Q(4)=3.167\times10^{-5}$) with $N=10^4$ samples, plainly and by importance sampling from $\mathcal{N}(4,1)$. (b) How many simulated bits does a plain Monte Carlo BER measurement at $P_e=10^{-9}$ need?
>
> **(a) Plain.** Relative standard error $\approx1/\sqrt{Np}=1/\sqrt{0.317}=1.78$ — with $10^4$ samples we expect to see $0.3$ events; the estimate is usually exactly zero. A 10% error would need $3.2\times10^6$ samples.
>
> **Importance sampling.** Draw $Y_i\sim\mathcal{N}(4,1)$; the likelihood ratio is $w(y)=\dfrac{e^{-y^2/2}}{e^{-(y-4)^2/2}}=e^{-4y+8}$. The estimator $\tfrac1N\sum\mathbf{1}\{Y_i>4\}w(Y_i)$ is unbiased. Its second moment is $E_q[\mathbf{1}w^2]=e^{16}Q(8)=4.53\times10^{-9}$, giving a relative standard error $\sqrt{(4.53\times10^{-9}-p^2)/N}/p=2.1\%$ at $N=10^4$. One run: $\hat p=3.12\times10^{-5}\pm0.07\times10^{-5}$ ✓ — about $7000$ times fewer samples than plain Monte Carlo for the same accuracy.
>
> **(b) BER.** $N\approx100/P_e=10^{11}$ bits for a 10% error. At $10^6$ simulated bits per second that is 28 hours per point on the BER curve. Importance sampling (biasing the noise variance or mean towards the decision boundary) or semi-analytic methods (simulate the noiseless waveform, compute the error probability analytically for each symbol) cut this to seconds.

> [!warning] Common Pitfall — importance sampling can make things worse
> The variance of an importance-sampling estimator is $E_q[f^2w^2]-\mu^2$, which can be *infinite* if $q$ has lighter tails than $fp$ (e.g. sampling a heavy-tailed $p$ from a Gaussian $q$). The estimator then converges erratically, with occasional enormous weights, and the sample standard error is meaningless. Always choose $q$ with tails at least as heavy as $p$ where $f\ne0$, and inspect the distribution of the weights; the **effective sample size** $(\sum w_i)^2/\sum w_i^2$ should not be tiny compared with $N$.

### Quasi-Monte Carlo

Pseudo-random points cluster and leave gaps by chance. **Low-discrepancy sequences** (Halton, Sobol', lattice rules) fill the cube more evenly by design. For integrands of bounded variation the Koksma–Hlawka inequality bounds the error by the variation times the discrepancy of the point set, $O\big((\log N)^d/N\big)$ — asymptotically almost $1/N$ instead of $1/\sqrt N$. The points are deterministic, so error bars are obtained by **randomised QMC** (random shifts or scrambling, repeated a few times). For smooth integrands in moderate dimension (up to tens), QMC routinely gains one to two orders of magnitude; for discontinuous indicators (pass/fail) the gain is smaller. SciPy provides `scipy.stats.qmc.Sobol`.

### Monte Carlo tolerance analysis of a circuit

The most common engineering use of Monte Carlo is **statistical tolerance analysis**: components are manufactured with random deviations from their nominal values; what is the distribution of a circuit's performance, and what fraction of units (the **yield**) meets specification?

> [!definition] Tolerance analysis
> Let a performance measure $y=g(p_1,\dots,p_m)$ depend on component values $p_k$ with nominal values $p_k^0$ and tolerances $\pm t_kp_k^0$. **Worst-case analysis** finds the extreme values of $y$ over the tolerance box. **Statistical analysis** models each $p_k$ as a random variable (commonly uniform on the tolerance band for precision parts, or Gaussian with $3\sigma$ equal to the tolerance), and estimates the distribution of $y$; the **yield** is $Y=P(y\in[y_\text{min},y_\text{max}])$.

The Monte Carlo procedure is:

1. Choose distributions (and correlations, e.g. resistors in the same thin-film network track each other) for all component values.
2. For $i=1,\dots,N$: draw a component set, evaluate $y_i=g(\mathbf{p}_i)$ (a formula or a full circuit simulation).
3. Report $\hat\mu_y$, $s_y$, percentiles and the yield estimate $\hat Y=\tfrac1N\#\{y_i\text{ in spec}\}$ with standard error $\sqrt{\hat Y(1-\hat Y)/N}$.

> [!example] Worked Example 24.30 — Gain spread of a non-inverting amplifier
> A non-inverting op-amp stage has gain $G=1+R_2/R_1$ with $R_1=1\ \mathrm{k\Omega}$ and $R_2=9\ \mathrm{k\Omega}$, both $\pm1\%$ resistors modelled as independent and uniform on their tolerance bands. The specification is $G=10\pm1\%$. Find the worst case, the standard deviation (by linearisation and by Monte Carlo) and the yield.
>
> **Worst case.** $G_\text{max}=1+9\times1.01/0.99=10.1818$, $G_\text{min}=1+9\times0.99/1.01=9.8218$: $\pm1.8\%$. Worst-case design says the spec cannot be guaranteed.
>
> **Linearisation.** The relative sensitivities are $S^G_{R_2}=\dfrac{R_2}{G}\dfrac{\partial G}{\partial R_2}=\dfrac{R_2/R_1}{G}=0.9$ and $S^G_{R_1}=-0.9$. A uniform $\pm1\%$ deviation has standard deviation $\sigma_r=0.01/\sqrt3=0.5774\%$, so
> $$
> \frac{\sigma_G}{G}\approx\sqrt{0.9^2+0.9^2}\;\sigma_r=1.273\times0.5774\%=0.735\%,\qquad\sigma_G\approx0.0735 .
> $$
>
> **Monte Carlo**, $N=10^5$:
>
> ```python
> rng = np.random.default_rng(12345)
> N = 100_000
> R1 = 1e3 * (1 + 0.01 * rng.uniform(-1, 1, N))
> R2 = 9e3 * (1 + 0.01 * rng.uniform(-1, 1, N))
> G = 1 + R2 / R1
> yield_ = np.mean(np.abs(G - 10) <= 0.1)
> ```
>
> Results: mean $10.0001$, standard deviation $0.0738$ (linearisation ✓), range $[9.8225,\ 10.1809]$ (inside the worst case, as it must be), yield within $\pm1\%$: $\hat Y=80.0\%$ (standard error $0.13\%$); within $\pm1.5\%$: $97.2\%$.
>
> **Why not just use the Gaussian approximation?** Treating $G$ as Gaussian with $\sigma=0.0735$ predicts a $\pm1\%$ yield of $2\Phi(0.1/0.0735)-1=82.6\%$. The true distribution of the sum of two uniform deviations is *triangular*, with thinner tails and more mass in the shoulders; the Gaussian overestimates the yield by 2.6 points. Monte Carlo makes no such assumption.
>
> **How many samples?** With only $N=1000$ samples, three seeds gave $78.6\%$, $79.3\%$ and $80.2\%$, consistent with the standard error $\sqrt{0.8\times0.2/1000}=1.3\%$. To quote yield to $\pm0.5\%$ (one standard error) needs $N=0.8\times0.2/0.005^2=6400$ runs.

> [!example] Worked Example 24.31 — A nonlinear performance: RC cutoff frequency
> An RC low-pass filter with $R=10\ \mathrm{k\Omega}\pm5\%$ and $C=10\ \mathrm{nF}\pm10\%$ (both uniform) has nominal $f_c=1/(2\pi RC)=1591.5\ \mathrm{Hz}$. The spec is $f_c$ within $\pm10\%$ of nominal. Compare linearised and Monte Carlo results.
>
> **Linearised.** $S^{f_c}_R=S^{f_c}_C=-1$, so $\sigma_{f_c}/f_c=\sqrt{(0.05/\sqrt3)^2+(0.10/\sqrt3)^2}=6.45\%$, $\sigma_{f_c}=102.7\ \mathrm{Hz}$, mean predicted $1591.5\ \mathrm{Hz}$.
>
> **Monte Carlo** ($N=10^5$). Mean $1598.4\ \mathrm{Hz}$, standard deviation $103.6\ \mathrm{Hz}$, $0.5\%$ and $99.5\%$ percentiles $1396.7$ and $1832.8\ \mathrm{Hz}$, i.e. $-12.2\%$ and $+15.2\%$; yield within $\pm10\%$: $87.4\%$.
>
> **Reading.** Because $f_c\propto1/(RC)$ is convex, the mean is *shifted upward* by $0.43\%$ (second-order effect: $E[1/X]\approx\tfrac1{\mu}(1+\sigma^2/\mu^2)$, and $(0.0645)^2\approx0.42\%$ ✓), and the distribution is skewed towards high frequencies. The linearised standard deviation is still accurate to 1%, but the asymmetric tails — which are what the yield depends on — are visible only in the simulation. Section 24.8 compares the two approaches systematically.

> [!tip] Engineering Note
> SPICE simulators implement this directly (`.MC` analysis in many dialects; Monte Carlo with `mismatch` and `process` statistical models in foundry PDKs). Foundry models separate **global (process)** variation, common to all devices on a die, from **local (mismatch)** variation, independent per device and scaling as $1/\sqrt{\text{area}}$ (Pelgrom's law); a correct Monte Carlo draws one global sample per run and independent local samples per device. Because each run is a full circuit simulation, the techniques of this section — fewer runs through variance reduction, QMC, or importance sampling towards the failure region for high-sigma SRAM yield (failure probabilities of $10^{-6}$–$10^{-9}$ per cell) — translate directly into days of compute saved.

> [!warning] Common Pitfall — reporting Monte Carlo results without error bars
> A Monte Carlo yield of "93.4%" from 200 runs has a standard error of $1.8\%$; the next digit is noise, and comparing it with "92.1%" from another design proves nothing. Always report the standard error (or a confidence interval), and when comparing two designs use **common random numbers** — the same random component deviations for both — so that the comparison is not swamped by sampling noise.

## 24.7 Parameter estimation

### Motivation

Models contain parameters that must be learned from measurements: the time constant and final value of a charging capacitor, the saturation current and ideality factor of a diode, the ESR and inductance of a capacitor from an impedance sweep, the poles and zeros of a loudspeaker from its frequency response, the twenty parameters of a transistor compact model. When the model is linear in the parameters, the problem is linear least squares ([Ch. 9](../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md)), solved in one step by QR (Section 24.3). Most physical models are *nonlinear* in at least some parameters — $\tau$ in $e^{-t/\tau}$, $n$ in $e^{V/(nV_T)}$, a resonance frequency — and need an iterative method. This section develops the Gauss–Newton and Levenberg–Marquardt algorithms, which exploit the least-squares structure far better than a general optimiser ([Ch. 23](23-optimization.md)), and shows how to attach uncertainties to the fitted parameters. The statistical justification — least squares is maximum likelihood under Gaussian noise — is in [Ch. 20](../part-6-probability-statistics/20-statistical-analysis.md).

### Problem statement

> [!definition] Nonlinear least squares
> Given data $(t_i,y_i)$, $i=1,\dots,m$, and a model $y=\phi(t;\boldsymbol\theta)$ with parameters $\boldsymbol\theta\in\mathbb{R}^n$ ($m>n$), define the **residuals** $r_i(\boldsymbol\theta)=y_i-\phi(t_i;\boldsymbol\theta)$ and the cost
> $$
> S(\boldsymbol\theta)=\tfrac12\sum_{i=1}^m w_i\,r_i(\boldsymbol\theta)^2=\tfrac12\,\mathbf{r}^\mathsf{T}\mathbf{W}\mathbf{r},
> $$
> with weights $w_i=1/\sigma_i^2$ (inverse noise variances; $\mathbf{W}=\mathbf{I}$ if the noise is uniform). The **Jacobian** $\mathbf{J}\in\mathbb{R}^{m\times n}$ of the model has entries $J_{ik}=\partial\phi(t_i;\boldsymbol\theta)/\partial\theta_k$, so that $\partial\mathbf{r}/\partial\boldsymbol\theta=-\mathbf{J}$.

The gradient and Hessian of $S$ (with $\mathbf{W}=\mathbf{I}$) are
$$
\nabla S=-\mathbf{J}^\mathsf{T}\mathbf{r},\qquad \nabla^2S=\mathbf{J}^\mathsf{T}\mathbf{J}-\sum_{i=1}^mr_i\nabla^2\phi(t_i;\boldsymbol\theta).
$$
The first Hessian term needs only first derivatives; the second is small when the residuals are small or the model is nearly linear. This structure is what the specialised methods exploit.

### Gauss–Newton

Linearise the model about the current estimate $\boldsymbol\theta_k$: $\phi(t_i;\boldsymbol\theta_k+\boldsymbol\delta)\approx\phi(t_i;\boldsymbol\theta_k)+\mathbf{J}_i\boldsymbol\delta$. The residual becomes $\mathbf{r}-\mathbf{J}\boldsymbol\delta$, and minimising $\lVert\mathbf{r}-\mathbf{J}\boldsymbol\delta\rVert^2$ is a *linear* least-squares problem.

> [!definition] Gauss–Newton iteration
> $$
> \boldsymbol\theta_{k+1}=\boldsymbol\theta_k+\boldsymbol\delta_k,\qquad \boldsymbol\delta_k=\arg\min_{\boldsymbol\delta}\lVert\mathbf{r}_k-\mathbf{J}_k\boldsymbol\delta\rVert_2,\quad\text{i.e.}\quad (\mathbf{J}_k^\mathsf{T}\mathbf{J}_k)\boldsymbol\delta_k=\mathbf{J}_k^\mathsf{T}\mathbf{r}_k .
> $$
> The step should be computed by QR of $\mathbf{J}_k$, not by forming $\mathbf{J}_k^\mathsf{T}\mathbf{J}_k$ (Section 24.3), whenever $\mathbf{J}$ is ill-conditioned.

Gauss–Newton is Newton's method for $\nabla S=\mathbf{0}$ with the second Hessian term dropped. Consequently:

> [!theorem] Local convergence of Gauss–Newton
> Near a minimiser $\boldsymbol\theta^*$ with $\mathbf{J}(\boldsymbol\theta^*)$ of full rank, Gauss–Newton converges
> - **quadratically** if the residuals at the solution are zero (the model fits exactly);
> - **linearly** otherwise, with rate approximately $\rho\big((\mathbf{J}^\mathsf{T}\mathbf{J})^{-1}\sum_ir_i\nabla^2\phi_i\big)$ evaluated at $\boldsymbol\theta^*$ — fast when the residuals are small relative to the model's curvature, and possibly divergent when they are large.
>
> Far from the solution it may diverge: the linearised model can be a poor approximation, and the full step can overshoot.

> [!example] Worked Example 24.32 — One Gauss–Newton step by hand
> Fit $v(t)=a\,e^{-bt}$ to $(t,v)=(0,5.0),\ (1,3.1),\ (2,1.8)$ starting from $(a,b)=(5,\ 0.5)$.
>
> **Model, residuals, Jacobian.** $\phi=(5,\ 5e^{-0.5},\ 5e^{-1})=(5,\ 3.03265,\ 1.83940)$; $\mathbf{r}=(0,\ 0.06735,\ -0.03940)$; $S_0=\sum r_i^2=0.006088$ (we quote $\sum r_i^2$, i.e. $2S$). With $\partial\phi/\partial a=e^{-bt}$ and $\partial\phi/\partial b=-a\,t\,e^{-bt}$,
> $$
> \mathbf{J}=\begin{bmatrix}1&0\\0.60653&-3.03265\\0.36788&-3.67879\end{bmatrix}.
> $$
>
> **Normal equations of the step.**
> $$
> \mathbf{J}^\mathsf{T}\mathbf{J}=\begin{bmatrix}1.50321&-3.19275\\-3.19275&22.73051\end{bmatrix},\qquad \mathbf{J}^\mathsf{T}\mathbf{r}=\begin{bmatrix}0.02635\\-0.05930\end{bmatrix}.
> $$
> The determinant is $1.50321\times22.73051-3.19275^2=23.975$, so
> $$
> \boldsymbol\delta=\frac1{23.975}\begin{bmatrix}22.73051&3.19275\\3.19275&1.50321\end{bmatrix}\begin{bmatrix}0.02635\\-0.05930\end{bmatrix}=\begin{bmatrix}0.01709\\-0.00021\end{bmatrix}.
> $$
>
> **New estimate.** $(a,b)=(5.01709,\ 0.49979)$, $\sum r_i^2=0.005625$. A library solver (`scipy.optimize.least_squares`) returns the same values to five digits and the same cost: one step has reached the minimum, because the starting point was close and the residuals small.
>
> **Contrast with the log-linear fit.** Fitting $\ln v=\ln a-bt$ by linear least squares gives $(a,b)=(5.055,\ 0.511)$ — different, because the log transform reweights the data (small values gain influence). It is an excellent *starting point*, not the answer to the original problem ([Ch. 9](../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md)).

### Levenberg–Marquardt

To make Gauss–Newton robust far from the solution, Levenberg (1944) and Marquardt (1963) added a damping term:

> [!definition] Levenberg–Marquardt step
> $$
> \big(\mathbf{J}_k^\mathsf{T}\mathbf{J}_k+\lambda_k\mathbf{D}_k\big)\boldsymbol\delta_k=\mathbf{J}_k^\mathsf{T}\mathbf{r}_k,
> $$
> with damping parameter $\lambda_k\ge0$ and a positive diagonal scaling $\mathbf{D}_k$ (Levenberg: $\mathbf{D}=\mathbf{I}$; Marquardt: $\mathbf{D}=\operatorname{diag}(\mathbf{J}^\mathsf{T}\mathbf{J})$, which makes the method invariant to the units of the parameters). A typical update: if the step reduces $S$, accept it and decrease $\lambda$ (e.g. $\lambda\leftarrow\lambda/10$); otherwise reject it, increase $\lambda$ ($\lambda\leftarrow10\lambda$) and recompute.

The damping interpolates between two methods:

- $\lambda\to0$: the Gauss–Newton step (fast near the solution);
- $\lambda\to\infty$: $\boldsymbol\delta\approx\lambda^{-1}\mathbf{D}^{-1}\mathbf{J}^\mathsf{T}\mathbf{r}=-\lambda^{-1}\mathbf{D}^{-1}\nabla S$, a short (scaled) steepest-descent step, which always reduces $S$ for $\lambda$ large enough.

Equivalently, the LM step minimises the linearised cost $\lVert\mathbf{r}-\mathbf{J}\boldsymbol\delta\rVert^2$ subject to $\lVert\mathbf{D}^{1/2}\boldsymbol\delta\rVert\le\Delta$ for some trust radius $\Delta$ that decreases as $\lambda$ increases: LM is a **trust-region** method, and modern implementations (MINPACK, `scipy.optimize.least_squares` with `method='lm'` or `'trf'`) adjust $\Delta$ directly by comparing the actual decrease of $S$ with the decrease predicted by the linear model. The LM system is again best solved as an augmented least-squares problem $\min\big\lVert\big[\begin{smallmatrix}\mathbf{J}\\ \sqrt\lambda\mathbf{D}^{1/2}\end{smallmatrix}\big]\boldsymbol\delta-\big[\begin{smallmatrix}\mathbf{r}\\ \mathbf{0}\end{smallmatrix}\big]\big\rVert$ by QR.

> [!example] Worked Example 24.33 — Fitting an RC charging curve: Gauss–Newton diverges, Levenberg–Marquardt converges
> A capacitor charging through a resistor is sampled every $0.25\ \mathrm{ms}$:
>
> | $t$ (ms) | 0.25 | 0.50 | 0.75 | 1.00 | 1.25 | 1.50 | 1.75 | 2.00 |
> |---|---|---|---|---|---|---|---|---|
> | $v$ (V) | 1.342 | 2.330 | 3.036 | 3.550 | 3.943 | 4.213 | 4.440 | 4.616 |
>
> (The data were synthesised from $V=5\ \mathrm{V}$, $\tau=0.8\ \mathrm{ms}$ with Gaussian noise of $\sigma=20\ \mathrm{mV}$.) Fit $v=V(1-e^{-t/\tau})$ from the poor initial guess $V=2\ \mathrm{V}$, $\tau=5\ \mathrm{ms}$.
>
> **Jacobian.** $\partial\phi/\partial V=1-e^{-t/\tau}$, $\partial\phi/\partial\tau=-V\,\dfrac{t}{\tau^2}e^{-t/\tau}$.
>
> **Gauss–Newton.** The initial cost is $\sum r_i^2=80.3\ \mathrm{V^2}$. The first GN step jumps to $V=-42.4\ \mathrm{V}$, $\tau=-167\ \mathrm{ms}$ — the linearisation about a nearly straight initial curve (with $\tau=5\ \mathrm{ms}$ the model is almost linear in $t$ over the data) wildly extrapolates — and the iteration then diverges to $\lvert V\rvert\sim10^{15}$.
>
> **Levenberg–Marquardt** (Marquardt scaling, $\lambda_0=10^{-2}$, factor 10):
>
> | iter | $V$ (V) | $\tau$ (ms) | $\sum r_i^2$ (V²) | $\lambda$ after step |
> |---|---|---|---|---|
> | 0 | 2.000 | 5.000 | 80.27 | $10^{-2}$ |
> | 1 | 3.197 | 1.466 | 28.67 | $1$ |
> | 2 | 4.358 | 0.5445 | 0.4349 | $10^{-1}$ |
> | 3 | 4.644 | 0.6725 | 0.09904 | $10^{-2}$ |
> | 4 | 4.939 | 0.7815 | 0.004695 | $10^{-3}$ |
> | 5 | 5.0127 | 0.80667 | 0.0013139 | $10^{-4}$ |
> | 9 | 5.015855 | 0.807757 | 0.0013084 | $10^{-7}$ |
>
> In the first iteration three trial steps were rejected ($\lambda$ raised to $10$) before a short, safe step was accepted; thereafter $\lambda$ falls by ten each iteration and the method becomes Gauss–Newton, converging rapidly. From the better start $(4\ \mathrm{V},\ 0.3\ \mathrm{ms})$ plain Gauss–Newton also converges, in about 6 iterations.
>
> **Reading.** The final residual $\sqrt{0.0013084/6}=14.8\ \mathrm{mV}$ (dividing by $m-n=6$ degrees of freedom) is consistent with the 20 mV noise. Fitted values: $V=5.016\ \mathrm{V}$, $\tau=0.808\ \mathrm{ms}$. How good are they? That needs the next subsection.

### Uncertainty of the fitted parameters

Near the solution the model is approximately linear, $\mathbf{y}\approx\boldsymbol\phi(\boldsymbol\theta^*)+\mathbf{J}(\boldsymbol\theta-\boldsymbol\theta^*)+\text{noise}$, so the linear least-squares theory applies with $\mathbf{J}$ in place of the design matrix.

> [!theorem] Asymptotic covariance of least-squares estimates
> If the measurement errors are independent with variance $\sigma^2$ and the model is correct, the estimate $\hat{\boldsymbol\theta}$ has approximately
> $$
> \operatorname{Cov}(\hat{\boldsymbol\theta})\approx\sigma^2(\mathbf{J}^\mathsf{T}\mathbf{J})^{-1},\qquad \hat\sigma^2=\frac{\sum_ir_i(\hat{\boldsymbol\theta})^2}{m-n},
> $$
> with $\mathbf{J}$ evaluated at $\hat{\boldsymbol\theta}$. With weights $\mathbf{W}=\operatorname{diag}(1/\sigma_i^2)$ the covariance is $(\mathbf{J}^\mathsf{T}\mathbf{W}\mathbf{J})^{-1}$. Standard errors are the square roots of the diagonal; the **correlation** between $\hat\theta_j$ and $\hat\theta_k$ is $C_{jk}/\sqrt{C_{jj}C_{kk}}$. For small samples, approximate 95% intervals use Student's $t_{m-n}$ quantile ([Ch. 20](../part-6-probability-statistics/20-statistical-analysis.md)).

*Sketch.* Linearising, $\hat{\boldsymbol\theta}-\boldsymbol\theta^*\approx(\mathbf{J}^\mathsf{T}\mathbf{J})^{-1}\mathbf{J}^\mathsf{T}\boldsymbol\epsilon$, whose covariance is $(\mathbf{J}^\mathsf{T}\mathbf{J})^{-1}\mathbf{J}^\mathsf{T}(\sigma^2\mathbf{I})\mathbf{J}(\mathbf{J}^\mathsf{T}\mathbf{J})^{-1}=\sigma^2(\mathbf{J}^\mathsf{T}\mathbf{J})^{-1}$. This is the linear error propagation of Section 24.8 applied to the map from data to estimate. $\square$

> [!example] Worked Example 24.34 — Error bars and correlation for the RC fit
> Continue Worked Example 24.33.
>
> **Covariance.** With $\hat\sigma=14.8\ \mathrm{mV}$ and $\mathbf{J}$ at the solution, $\hat\sigma^2(\mathbf{J}^\mathsf{T}\mathbf{J})^{-1}$ gives standard errors $u(V)=0.020\ \mathrm{V}$ and $u(\tau)=0.0077\ \mathrm{ms}$, and a correlation coefficient $+0.93$.
>
> **Intervals.** With $t_{6,0.975}=2.447$: $V=5.016\pm0.049\ \mathrm{V}$, $\tau=0.808\pm0.019\ \mathrm{ms}$ (95%). Both contain the true values $5\ \mathrm{V}$ and $0.8\ \mathrm{ms}$ ✓.
>
> **The correlation matters.** The data stop at $t=2\ \mathrm{ms}=2.5\tau$, before the curve has settled; a larger final value with a slower time constant fits almost as well as a smaller one with a faster time constant. The individual error bars understate how well the *combination* is determined, and any quantity derived from both (e.g. the initial slope $V/\tau$) must be propagated with the full covariance (Section 24.8). Extending the record to $5\tau$ would decorrelate them.
>
> **Scaling.** In SI units the columns of $\mathbf{J}$ differ by a factor of order $10^3$–$10^4$ (V per V versus V per second), and $\kappa(\mathbf{J}^\mathsf{T}\mathbf{J})=5.1\times10^7$. Measuring $\tau$ in milliseconds gives $\kappa=65$. The answer is the same in exact arithmetic, but a well-scaled problem is kinder to every algorithm — Marquardt's diagonal scaling achieves the same effect automatically.

### Identifiability and good formulations

> [!definition] Identifiability
> A parameter vector is **(locally) identifiable** from an experiment if distinct parameter values near $\boldsymbol\theta^*$ produce distinct model outputs, i.e. $\mathbf{J}(\boldsymbol\theta^*)$ has full column rank. It is **practically identifiable** if, in addition, $\mathbf{J}^\mathsf{T}\mathbf{J}$ is well enough conditioned that the noise does not produce huge parameter uncertainties.

Practical non-identifiability is common and shows up as near-singular $\mathbf{J}^\mathsf{T}\mathbf{J}$, correlations near $\pm1$ and enormous standard errors. The remedies are experimental (measure over a wider range, excite the relevant dynamics) or mathematical (reparametrise, fix a parameter, add prior information as regularisation).

> [!example] Worked Example 24.35 — Diode parameters: correlation by design
> The forward diode law $I=I_s e^{V/(nV_T)}$ becomes linear after taking logarithms: $\ln I=\ln I_s+\dfrac{1}{nV_T}V$, with unknowns $\theta_1=\ln I_s$ and $\theta_2=1/(nV_T)$. Measurements are taken at 11 equally spaced voltages from $0.5\ \mathrm{V}$ to $0.7\ \mathrm{V}$. Find the correlation between the estimates, and fix it.
>
> **Correlation.** The design matrix has columns $\mathbf{1}$ and $\mathbf{V}$. With equal noise, $\operatorname{Cov}\propto(\mathbf{X}^\mathsf{T}\mathbf{X})^{-1}$, and the correlation is $-\bar V/\sqrt{\overline{V^2}}=-0.6/\sqrt{0.36+0.004}=-0.994$. The intercept $\ln I_s$ is an extrapolation to $V=0$, far outside the data; any error in the slope swings it strongly.
>
> **Reparametrisation.** Write $\ln I=\ln I_0+\theta_2(V-0.6)$, where $I_0$ is the current at the centre of the data. The columns $\mathbf{1}$ and $\mathbf{V}-\bar V$ are orthogonal, the correlation is exactly 0, and $\ln I_0$ is determined to the noise level divided by $\sqrt{11}$. $I_s=I_0e^{-0.6\theta_2}$ can be computed afterwards — with its large but now honestly propagated uncertainty.
>
> **Reading.** Parametrise models by quantities the data actually measure: a value at the centre of the data rather than an extrapolated intercept, $\log$ of scale parameters that span decades (like $I_s\sim10^{-14}\ \mathrm{A}$), time constants rather than poles when the data are in the time domain.

> [!warning] Common Pitfall — local minima and bad starting values
> Nonlinear least squares can have several local minima — fitting a sum of two exponentials, a resonance frequency, or a phase is notorious — and every method here finds only a local one. Get starting values from physics (the final value of a step response, the $-3\ \mathrm{dB}$ point of a response), from a linearising transformation (Worked Example 24.32), or by a coarse grid search over the nonlinear parameters; and for models linear in some parameters, eliminate those analytically (**variable projection**: for each trial $\tau$, the best $V$ is a one-line linear least-squares solve), leaving a lower-dimensional search.

> [!tip] Engineering Note
> Device-model extraction (BSIM transistor models, IBIS buffer models, battery equivalent circuits) is large-scale weighted nonlinear least squares, usually with LM and carefully staged parameter subsets so that each stage is identifiable. Impedance spectroscopy fits equivalent circuits to complex data by stacking real and imaginary residuals, $\mathbf{r}=[\operatorname{Re}\{\mathbf{Z}_\text{meas}-\mathbf{Z}_\text{model}\},\ \operatorname{Im}\{\cdots\}]$, often weighted by $1/\lvert Z\rvert$ so that every decade of frequency counts equally. Vector fitting, used to build rational macromodels of interconnects from S-parameter data, is an iterated *linear* least-squares method that relocates poles — an example of choosing a formulation that avoids nonlinear search altogether.

## 24.8 Sensitivity and uncertainty propagation

### Motivation

Every computed result inherits uncertainty from its inputs: component tolerances, measured parameters, fitted model coefficients, environmental conditions. Two related questions follow. **Sensitivity analysis** asks how much the output changes per unit change of each input — which resistor sets the gain error, which permittivity dominates the resonance shift — and therefore where to spend money on tighter tolerances. **Uncertainty propagation** asks what the resulting uncertainty of the output is, so that a measurement or a simulation can be reported as "$y=\hat y\pm U$". Both depend on derivatives of the output with respect to the inputs, and computing those derivatives accurately and cheaply for large models is itself a computational problem. This section treats derivative computation (finite differences, complex step, automatic differentiation and the adjoint method), then linearised propagation in the framework of the *Guide to the Expression of Uncertainty in Measurement* (GUM), then Monte Carlo propagation and the conditions under which the two disagree.

### Sensitivities

> [!definition] Sensitivity
> For an output $y=f(p_1,\dots,p_m)$, the **absolute sensitivity** to $p_k$ is $\partial y/\partial p_k$ and the **relative (normalised, Bode) sensitivity** is
> $$
> S^y_{p_k}=\frac{p_k}{y}\frac{\partial y}{\partial p_k}=\frac{\partial\ln\lvert y\rvert}{\partial\ln\lvert p_k\rvert},
> $$
> the percentage change in $y$ per percent change in $p_k$. The vector of absolute sensitivities is the gradient $\nabla_\mathbf{p}y$; for several outputs it is the **sensitivity matrix** (Jacobian) $\mathbf{S}=\partial\mathbf{y}/\partial\mathbf{p}$.

Relative sensitivities are dimensionless and obey simple rules: $S^{y_1y_2}_p=S^{y_1}_p+S^{y_2}_p$, $S^{1/y}_p=-S^y_p$, $S^{y^n}_p=nS^y_p$. For an RC cutoff $f_c=1/(2\pi RC)$, $S^{f_c}_R=S^{f_c}_C=-1$; for an inverting amplifier $G=-R_2/R_1$, $S^G_{R_2}=1$, $S^G_{R_1}=-1$.

### Computing derivatives numerically

**Finite differences.** The forward difference $\big(f(x+h)-f(x)\big)/h$ has truncation error $\tfrac12hf''$ and rounding error about $2\varepsilon\lvert f\rvert/h$ ([Ch. 22](22-numerical-methods.md)). Balancing them gives $h_\text{opt}\approx\sqrt{\varepsilon}\lvert x\rvert\approx10^{-8}\lvert x\rvert$ and a best relative accuracy of about $\sqrt\varepsilon\approx10^{-8}$ — half the digits. The central difference $\big(f(x+h)-f(x-h)\big)/(2h)$ has truncation error $\tfrac16h^2f'''$, giving $h_\text{opt}\approx\varepsilon^{1/3}$ and accuracy about $\varepsilon^{2/3}\approx10^{-11}$. Each gradient of $m$ parameters costs $m$ (forward) or $2m$ (central) extra evaluations of $f$.

**Complex step.** If $f$ is real-analytic and implemented with operations that extend to complex arguments, Taylor expansion gives $f(x+jh)=f(x)+jhf'(x)-\tfrac12h^2f''(x)-\tfrac16jh^3f'''(x)+\cdots$, so
$$
f'(x)=\frac{\operatorname{Im}\{f(x+jh)\}}{h}+O(h^2).
$$
There is **no subtraction**, hence no cancellation: $h$ can be made as small as $10^{-20}$ or $10^{-100}$ and the derivative is accurate to machine precision. The cost is one complex evaluation per parameter. Circuit and EM codes that already use complex arithmetic for phasors can often be differentiated this way with almost no change — but functions such as `abs`, `max` or comparisons must be redefined for complex arguments, and the code must not take the conjugate of the perturbed variable.

> [!example] Worked Example 24.36 — Finite differences versus complex step
> Differentiate $f(x)=\dfrac{e^x}{\sqrt{\sin^3x+\cos^3x}}$ at $x=1.5$ (exact $f'(1.5)=4.053\,427\,894$). Relative errors:
>
> | $h$ | forward difference | central difference | complex step |
> |---|---|---|---|
> | $10^{-2}$ | $1.2\times10^{-2}$ | $1.3\times10^{-4}$ | $1.3\times10^{-4}$ |
> | $10^{-4}$ | $1.2\times10^{-4}$ | $1.3\times10^{-8}$ | $1.3\times10^{-8}$ |
> | $10^{-6}$ | $1.2\times10^{-6}$ | $1.4\times10^{-10}$ | $1.3\times10^{-12}$ |
> | $10^{-8}$ | $2.5\times10^{-8}$ | $8.4\times10^{-9}$ | $<10^{-16}$ |
> | $10^{-10}$ | $3.5\times10^{-6}$ | $2.0\times10^{-7}$ | $<10^{-16}$ |
> | $10^{-12}$ | $2.7\times10^{-4}$ | $5.4\times10^{-5}$ | $<10^{-16}$ |
> | $10^{-20}$ | (meaningless) | (meaningless) | $2.2\times10^{-16}$ |
>
> **Reading.** Finite differences reach their best accuracy at an intermediate $h$ ($10^{-8}$ forward, about $10^{-6}$ central) and get *worse* as $h$ decreases further — the V-shaped error curve of [Ch. 22](22-numerical-methods.md). The complex step improves monotonically and saturates at machine precision. For both large and small $h$, the central difference and complex step share the same $O(h^2)$ truncation error (identical first two rows).

**Automatic differentiation (AD).** A program is a composition of elementary operations whose derivatives are known. AD applies the chain rule to the *program*, operation by operation, producing derivatives exact to rounding error with no step size at all. **Forward mode** propagates, alongside each intermediate value $v$, its derivative $\dot v=\partial v/\partial p$ with respect to *one* input (dual numbers $v+\dot v\epsilon$ with $\epsilon^2=0$; the complex step is a close relative); a full gradient needs $m$ passes. **Reverse mode** first evaluates the program and records it, then propagates *adjoints* $\bar v=\partial y/\partial v$ backwards from *one* output; a full gradient with respect to all $m$ inputs costs a small constant multiple (typically 2–5) of one function evaluation, independent of $m$. Reverse-mode AD is "backpropagation" in neural networks and is what JAX, PyTorch and similar tools implement.

### The adjoint method

The reverse-mode idea, applied to a model defined by a linear system, is the **adjoint method**. It is the key to sensitivity analysis of large circuits and field models, where there may be thousands of parameters but only a few outputs.

Let the state $\mathbf{x}\in\mathbb{R}^n$ be defined by $\mathbf{A}(\mathbf{p})\mathbf{x}=\mathbf{b}(\mathbf{p})$ and the output be $y=\mathbf{c}^\mathsf{T}\mathbf{x}$ (a node voltage, a port current). Differentiating the system with respect to $p_k$:
$$
\mathbf{A}\frac{\partial\mathbf{x}}{\partial p_k}=\frac{\partial\mathbf{b}}{\partial p_k}-\frac{\partial\mathbf{A}}{\partial p_k}\mathbf{x}.
$$
The **direct** method solves this for $\partial\mathbf{x}/\partial p_k$ — one solve per parameter — and forms $\partial y/\partial p_k=\mathbf{c}^\mathsf{T}\partial\mathbf{x}/\partial p_k$. The **adjoint** method notices that only the combination $\mathbf{c}^\mathsf{T}\mathbf{A}^{-1}$ is needed.

> [!theorem] Adjoint sensitivity formula
> Let $\boldsymbol\lambda$ solve the **adjoint system** $\mathbf{A}^\mathsf{T}\boldsymbol\lambda=\mathbf{c}$. Then for every parameter $p_k$
> $$
> \frac{\partial y}{\partial p_k}=\boldsymbol\lambda^\mathsf{T}\Big(\frac{\partial\mathbf{b}}{\partial p_k}-\frac{\partial\mathbf{A}}{\partial p_k}\mathbf{x}\Big).
> $$
> All $m$ sensitivities of one output cost **one** extra solve (with the transposed matrix, whose LU factors are those of $\mathbf{A}$ transposed) plus $m$ cheap inner products — instead of $m$ solves.

*Proof.* $\dfrac{\partial y}{\partial p_k}=\mathbf{c}^\mathsf{T}\dfrac{\partial\mathbf{x}}{\partial p_k}=\mathbf{c}^\mathsf{T}\mathbf{A}^{-1}\Big(\dfrac{\partial\mathbf{b}}{\partial p_k}-\dfrac{\partial\mathbf{A}}{\partial p_k}\mathbf{x}\Big)$, and $\mathbf{c}^\mathsf{T}\mathbf{A}^{-1}=(\mathbf{A}^{-\mathsf{T}}\mathbf{c})^\mathsf{T}=\boldsymbol\lambda^\mathsf{T}$. $\square$

The rule of thumb: **direct** sensitivities when there are few parameters and many outputs; **adjoint** sensitivities when there are many parameters and few outputs (the common case in design optimisation and tolerance analysis).

For a conductance $g$ connected between nodes $a$ and $b$ of a nodal matrix, $\partial\mathbf{G}/\partial g=(\mathbf{e}_a-\mathbf{e}_b)(\mathbf{e}_a-\mathbf{e}_b)^\mathsf{T}$ (the stamp of a unit conductance), so
$$
\frac{\partial y}{\partial g}=-(\lambda_a-\lambda_b)(v_a-v_b):
$$
the product of the branch voltage in the original circuit and the branch voltage in the adjoint circuit. This is the classical **adjoint-network** result derived from Tellegen's theorem (Director and Rohrer, 1969). For a reciprocal (symmetric) network the adjoint network is the original network, driven by a unit current source at the output node.

> [!example] Worked Example 24.37 — Sensitivities of a node voltage to every conductance with one extra solve
> The three-node network of Worked Example 24.17 has conductances (in siemens) $g_{12}=2$, $g_{13}=1$, $g_{23}=1$ between nodes and $g_1=1$, $g_2=1$, $g_3=2$ from each node to ground, and source currents $\mathbf{b}=(8,2,1)\ \mathrm{A}$, giving $\mathbf{G}=\begin{bmatrix}4&-2&-1\\-2&4&-1\\-1&-1&4\end{bmatrix}$ and $\mathbf{v}=(4,3,2)\ \mathrm{V}$. Find $\partial v_3/\partial g$ for all six conductances and $\partial v_3/\partial b_k$ for all three sources.
>
> **Adjoint solve.** $\mathbf{G}^\mathsf{T}\boldsymbol\lambda=\mathbf{e}_3$ ($\mathbf{G}$ is symmetric). Check that $\boldsymbol\lambda=(\tfrac16,\tfrac16,\tfrac13)$: row 1, $\tfrac46-\tfrac26-\tfrac13=0$ ✓; row 2 similarly $0$ ✓; row 3, $-\tfrac16-\tfrac16+\tfrac43=1$ ✓.
>
> **Source sensitivities.** $\partial\mathbf{b}/\partial b_k=\mathbf{e}_k$, so $\partial v_3/\partial b_k=\lambda_k$: $(\tfrac16,\tfrac16,\tfrac13)\ \Omega$ — these are transfer resistances, by reciprocity.
>
> **Conductance sensitivities** $-(\lambda_a-\lambda_b)(v_a-v_b)$ (for a grounded element, $\lambda_b=v_b=0$):
>
> | element | $\lambda_a-\lambda_b$ | $v_a-v_b$ (V) | $\partial v_3/\partial g$ (V/S) | relative $S^{v_3}_g$ |
> |---|---|---|---|---|
> | $g_{12}$ | $0$ | $1$ | $0$ | $0$ |
> | $g_{13}$ | $-\tfrac16$ | $2$ | $+0.3333$ | $+0.167$ |
> | $g_{23}$ | $-\tfrac16$ | $1$ | $+0.1667$ | $+0.083$ |
> | $g_1$ | $\tfrac16$ | $4$ | $-0.6667$ | $-0.333$ |
> | $g_2$ | $\tfrac16$ | $3$ | $-0.5000$ | $-0.250$ |
> | $g_3$ | $\tfrac13$ | $2$ | $-0.6667$ | $-0.667$ |
>
> **Check by finite difference.** Increasing $g_{13}$ by $10^{-6}$ and re-solving changes $v_3$ by $3.3333\times10^{-7}$ ✓; every row agrees to six digits.
>
> **Reading.** One factorisation, one extra (transposed) solve, and six inner products give every sensitivity; the direct method would need six solves, and for a circuit with $10^4$ elements the adjoint saves a factor of $10^4$. The table also explains the physics: $g_{12}$ does not affect $v_3$ at all because nodes 1 and 2 are at the same potential in the adjoint circuit, and the shunt $g_3$ at the output node has the largest relative influence.

### Linearised uncertainty propagation (GUM)

> [!definition] Standard uncertainty and the GUM law of propagation
> The **standard uncertainty** $u(x)$ of an input estimate $x$ is its estimated standard deviation — evaluated statistically from repeated observations (**Type A**) or from other information such as a datasheet tolerance or calibration certificate (**Type B**; a uniform interval $\pm a$ gives $u=a/\sqrt3$). For a measurement model $Y=f(X_1,\dots,X_n)$ with estimates $x_i$, the **combined standard uncertainty** of $y=f(x_1,\dots,x_n)$ is, to first order,
> $$
> u_c^2(y)=\sum_{i=1}^nc_i^2u^2(x_i)+2\sum_{i<j}c_ic_j\,u(x_i,x_j),\qquad c_i=\frac{\partial f}{\partial x_i}\Big|_{\mathbf{x}},
> $$
> where $c_i$ are the **sensitivity coefficients** and $u(x_i,x_j)$ the covariances. In matrix form, for a vector output, $\mathbf{U}_y=\mathbf{J}\mathbf{U}_x\mathbf{J}^\mathsf{T}$. The **expanded uncertainty** is $U=k\,u_c(y)$ with **coverage factor** $k$, usually $k=2$ (about 95% coverage for a Gaussian output).

*Derivation.* Expand $f$ to first order, $Y-y\approx\sum_ic_i(X_i-x_i)$; the variance of a linear combination is $\sum_{i,j}c_ic_j\operatorname{Cov}(X_i,X_j)$. $\square$

For products and quotients of powers, $y=\prod x_i^{a_i}$ with uncorrelated inputs, the rule takes the convenient relative form $\big(u_c(y)/y\big)^2=\sum_ia_i^2\big(u(x_i)/x_i\big)^2$ — the relative sensitivities are just the exponents.

> [!example] Worked Example 24.38 — Power from voltage and resistance
> The power dissipated in a resistor is computed as $P=V^2/R$ from $V=5.000\ \mathrm{V}$ with $u(V)=0.020\ \mathrm{V}$ and $R=100.0\ \Omega$ with $u(R)=0.50\ \Omega$ (independent). Find $u_c(P)$ and the expanded uncertainty.
>
> **Estimate.** $P=25/100=0.2500\ \mathrm{W}$.
>
> **Relative form.** Exponents $2$ and $-1$: $\big(u_c/P\big)^2=(2\times0.004)^2+(0.005)^2=6.4\times10^{-5}+2.5\times10^{-5}=8.9\times10^{-5}$, so $u_c/P=0.943\%$ and $u_c(P)=2.36\ \mathrm{mW}$.
>
> **Budget.** The voltage contributes $64/89=72\%$ of the variance despite its smaller relative uncertainty, because of the square. A better voltmeter is the cheaper improvement.
>
> **Result.** $P=(250.0\pm4.7)\ \mathrm{mW}$ with $k=2$.
>
> **Monte Carlo check** ($10^5$ Gaussian draws): mean $250.00\ \mathrm{mW}$, standard deviation $2.35\ \mathrm{mW}$ ✓. The model is nearly linear over $\pm1\%$, so linearisation is excellent.

> [!example] Worked Example 24.39 — Correlated inputs: a quantity derived from fitted parameters
> From the RC fit of Worked Example 24.34, compute the initial slope $s=V/\tau$ and its uncertainty, with and without the correlation $\rho=0.931$.
>
> **Estimate.** $s=5.0159/0.80776\times10^{-3}=6210\ \mathrm{V/s}$.
>
> **Sensitivities (relative).** $s=V\tau^{-1}$: exponents $+1$ and $-1$. Relative uncertainties $u(V)/V=0.398\%$, $u(\tau)/\tau=0.949\%$.
>
> **Ignoring correlation.** $\sqrt{0.398^2+0.949^2}\%=1.03\%$, $u(s)=64\ \mathrm{V/s}$.
>
> **With correlation.** The cross term is $2\times(+1)(-1)\times\rho\times0.398\%\times0.949\%$:
> $$
> \Big(\frac{u(s)}{s}\Big)^2=0.398^2+0.949^2-2(0.931)(0.398)(0.949)=0.1584+0.9006-0.7033=0.3557\ (\%)^2,
> $$
> so $u(s)/s=0.596\%$, $u(s)=37\ \mathrm{V/s}$.
>
> **Reading.** $V$ and $\tau$ err in the *same* direction, so their ratio — which the early data points measure directly — is much better determined than either. Dropping the covariance overstates the uncertainty by 73% here; with $\rho<0$ it would understate it. When propagating from fitted parameters, always use the full covariance matrix $\mathbf{J}_s\mathbf{C}\mathbf{J}_s^\mathsf{T}$.

### Monte Carlo propagation and the limits of linearisation

The **GUM Supplement 1** approach replaces the first-order formula by simulation: assign a probability distribution to each input (Gaussian, uniform, triangular, or a joint distribution for correlated inputs), draw $N$ samples (typically $10^5$–$10^6$), evaluate the model for each, and summarise the output distribution by its mean, standard deviation and a *coverage interval* taken from the empirical quantiles (e.g. the 2.5% and 97.5% percentiles). It is exactly the tolerance analysis of Section 24.6, with measurement uncertainties in place of tolerances.

Linearisation and Monte Carlo agree when the model is approximately linear over the region where the input distributions have most of their mass, i.e. when second-order terms $\tfrac12f''u^2$ are small compared with $f'u$. They disagree when:

1. **The model is strongly nonlinear over the uncertainty range** (large relative uncertainties; exponentials, ratios with uncertain denominators near zero).
2. **The first derivative vanishes** at the estimate (an extremum: maximum power transfer, the peak of a resonance, a balanced bridge), so the linear formula predicts zero uncertainty.
3. **The output distribution is far from Gaussian** (bounded, skewed or bimodal), so that $y\pm2u_c$ is not a 95% interval.
4. **The quantity of interest is a probability or a quantile** (yield, coverage) rather than a standard deviation.

For case 2 the GUM provides a second-order correction: for a single Gaussian input, $u_c^2=c_1^2u^2+\tfrac12\big(f''\big)^2u^4$.

> [!example] Worked Example 24.40 — Matched load: where linearisation predicts zero
> A source with $R_s=50\ \Omega$ drives a load $R_L$; the delivered power relative to the maximum is $\eta=\dfrac{4R_sR_L}{(R_s+R_L)^2}$. The load is nominally $50\ \Omega$ with a 5% Gaussian standard uncertainty. Find the uncertainty of $\eta$ by GUM (first and second order) and by Monte Carlo.
>
> **First order.** With $x=R_L/R_s$, $\eta=4x/(1+x)^2$, $d\eta/dx=4(1-x)/(1+x)^3=0$ at $x=1$. The linear law gives $u_c(\eta)=0$: "the power is exactly 100%".
>
> **Second order.** Expanding about $x=1+e$: $\eta=\dfrac{1+e}{(1+e/2)^2}\approx(1+e)(1-e+\tfrac34e^2)=1-\tfrac14e^2$. So $\eta''=-\tfrac12$, and with $e\sim\mathcal{N}(0,0.05^2)$:
> $$
> E[\eta]\approx1-\tfrac14(0.05)^2=0.999375,\qquad u(\eta)=\tfrac14\sqrt{\operatorname{Var}(e^2)}=\tfrac14\sqrt2(0.05)^2=0.000884 .
> $$
>
> **Monte Carlo** ($10^5$ samples). Mean $0.999373$, standard deviation $0.000896$, median $0.99972$, 95% coverage interval $[0.99685,\ 0.99999]$.
>
> **Reading.** The second-order GUM matches the mean and standard deviation, but the distribution is one-sided (it can never exceed 1) and strongly skewed: $\eta\pm2u$ would include impossible values above 1. The Monte Carlo coverage interval is the honest answer. Mismatch is always a *loss*, never a gain — a qualitative fact the linear method cannot express.

```desmos-graph
left=0.7; right=1.3; bottom=0.97; top=1.005
---
y=\frac{4x}{(1+x)^{2}}|blue
y=1-\frac{(x-1)^{2}}{4}|red|dashed
y=1|grey|dotted
x=0.95|green|dashed
x=1.05|green|dashed
```

*Relative delivered power $\eta=4x/(1+x)^2$ against $x=R_L/R_s$ (blue) with its quadratic approximation $1-(x-1)^2/4$ (red, dashed). The green lines mark $\pm1$ standard uncertainty of the load. At the maximum the tangent is horizontal, so a first-order propagation sees no uncertainty at all; the curvature maps a symmetric input spread into a one-sided output spread below 1.*

| | Linearised (GUM) | Monte Carlo (GUM-S1) |
|---|---|---|
| Cost | one evaluation plus $m$ derivatives (one adjoint solve for many inputs) | $10^4$–$10^6$ evaluations |
| Nonlinearity | first order (second order with extra effort) | exact for the assumed input distributions |
| Output distribution | assumed Gaussian (or Student $t$) | estimated |
| Correlated inputs | covariance terms | joint sampling |
| Insight | uncertainty budget: contribution of each input | needs extra runs (one-at-a-time or variance-based indices) |
| Best for | small uncertainties, many inputs, design iterations | large uncertainties, extrema, yield and coverage questions, validation of the linear result |

> [!tip] Engineering Note — global sensitivity
> Derivatives describe *local* sensitivity at the nominal point. When inputs vary over wide ranges, **variance-based (Sobol') indices** quantify how much of $\operatorname{Var}(y)$ each input explains on its own ($S_i=\operatorname{Var}(E[y\mid X_i])/\operatorname{Var}(y)$) and through interactions. They are estimated by structured Monte Carlo designs and are the standard way to decide which of fifty process parameters actually matter for a circuit's yield. For a linear model with independent inputs they reduce to the normalised terms $c_i^2u_i^2/u_c^2$ of the GUM uncertainty budget.

> [!warning] Common Pitfall — tolerance is not standard uncertainty
> A "$\pm1\%$" resistor tolerance is a bound, not a standard deviation. Treating it as $u=1\%$ overstates the uncertainty by $\sqrt3$ for a uniform distribution, and the true distribution of a purchased batch is often neither uniform nor Gaussian (parts near nominal may have been sorted into a tighter tolerance grade, leaving a bimodal distribution). State the assumed distribution explicitly, and where it matters, measure a sample.

## 24.9 Reproducible numerical experiments

### Motivation

A numerical result is only useful if it can be trusted, and it can only be trusted if it can be checked. Six months after a tolerance study, someone asks why the yield was quoted as 87% when a re-run gives 85%. A colleague's filter-design script produces different coefficients on a new laptop. A regression test fails after a library upgrade, and nobody knows whether the old or the new answer is right. Each of these is a failure of **reproducibility**, and each has a mathematical component (random-number streams, floating-point non-associativity, ill-conditioning) as well as an organisational one (versions, environments, records). This section treats computations as *experiments* that must be designed, recorded and tested like laboratory measurements.

> [!definition] Levels of reproducibility
> - **Repeatability (bitwise reproducibility):** the same code, data, parameters and environment produce bit-for-bit identical output.
> - **Reproducibility:** a re-run by someone else, possibly on different hardware or library versions, produces results that agree *within a stated tolerance* justified by the numerical and statistical error analysis.
> - **Replicability:** an independent implementation of the same method, or a new experiment, leads to the same scientific or engineering conclusion.
>
> Bitwise repeatability is a debugging tool; reproducibility within tolerance is the engineering requirement.

### Sources of irreproducibility

| Source | Mechanism | Control |
|---|---|---|
| random numbers | different seeds, shared global generator, parallel streams overlapping | explicit seeded generators, independent spawned streams, seeds recorded |
| floating-point order | $(a+b)+c\ne a+(b+c)$; parallel reductions, multithreaded BLAS, vectorised loops, fused multiply–add, GPU atomics change the order | accept and bound the variation; deterministic modes where needed; compensated summation |
| library versions | changed algorithms, defaults (tolerances, orderings, pivoting thresholds), bug fixes | pinned versions, environment files, containers |
| ill-conditioning | tiny rounding differences amplified by $\kappa$ | estimate $\kappa$; set test tolerances accordingly; reformulate |
| data and parameters | edited input files, manual steps, undocumented options | version control for code *and* inputs; configuration files; scripted pipelines |
| hardware and OS | different CPU instruction sets, math-library implementations of $\exp$, $\sin$ | record platform; tolerance-based comparison |

### Floating-point non-associativity

Floating-point addition is commutative but not associative: each operation rounds, and the rounding depends on the operands' magnitudes ([Ch. 22](22-numerical-methods.md)). Any change in the *order* of a reduction can therefore change the last few bits of the result — and multithreaded libraries choose the order dynamically, depending on the number of threads and on how the data are split.

> [!example] Worked Example 24.41 — The order of a sum
> (a) Compare $(0.1+0.2)+0.3$ with $0.1+(0.2+0.3)$, and $(10^{16}+1)-10^{16}$ with $(10^{16}-10^{16})+1$. (b) Sum $10^6$ copies of $0.1$ by a naive loop, by NumPy's `sum`, and by `math.fsum`. (c) Sum $10^6$ standard-normal samples in forward, reversed and sorted order.
>
> **(a)** $(0.1+0.2)+0.3=0.6000000000000001$ but $0.1+(0.2+0.3)=0.6$. $(10^{16}+1)-10^{16}=0$ (the spacing of doubles near $10^{16}$ is 2, so $10^{16}+1$ rounds back to $10^{16}$), but $(10^{16}-10^{16})+1=1$.
>
> **(b)** Naive loop: $100\,000.000\,001\,332\,88$ (relative error $1.3\times10^{-11}$: each addition rounds, and errors accumulate roughly as $n\varepsilon$ in the worst case). NumPy `sum`: $100\,000.000\,000\,000\,03$ — NumPy uses **pairwise summation**, whose error grows like $\varepsilon\log_2n$. `math.fsum`: exactly $100\,000.0$, the correctly rounded sum of the stored values.
>
> **(c)** Forward: $998.570\,649\,438\,621\,3$; reversed: $998.570\,649\,438\,621\,4$; sorted ascending: $998.570\,649\,438\,654\,1$; `fsum`: $998.570\,649\,438\,621\,3$. The results differ in the 12th–14th significant digit.
>
> **Reading.** None of these is a bug. A test that compares such results with `==` will fail on another machine, or with a different thread count. Results that must agree should be compared with a tolerance derived from the error analysis — here about $n\varepsilon\max\lvert x_i\rvert/\lvert\sum x_i\rvert$, far larger than $\varepsilon$ when there is cancellation. When bitwise repeatability matters (regression debugging, audit), fix the thread count, use deterministic library modes, or use compensated (Kahan) or exact summation in the few places that need it.

### Random-number streams

Pseudo-random number generators (PRNGs) are deterministic: the same seed gives the same sequence. Reproducible Monte Carlo therefore requires that every random draw comes from a generator whose seed is chosen deliberately and recorded.

> [!definition] Seeded generators and independent streams
> A **seed** is the initial state from which a PRNG's sequence is generated. A **stream** is the sequence produced by one generator. For parallel or multi-component simulations, each worker or component needs its own **statistically independent stream**, obtained not by seeds $1,2,3,\dots$ chosen ad hoc but by a mechanism designed to avoid overlap and correlation — in NumPy, `SeedSequence.spawn`, which derives child seeds by hashing; in other libraries, jump-ahead or counter-based generators (Philox, Threefry).

```python
import numpy as np
root = np.random.SeedSequence(20241021)          # record this number
children = root.spawn(8)                          # one per worker / per component group
rngs = [np.random.default_rng(s) for s in children]
# worker k draws only from rngs[k]; results are independent of scheduling
R1 = 1e3 * (1 + 0.01 * rngs[0].uniform(-1, 1, 10_000))
```

Three rules follow. Never use the legacy global generator (`np.random.seed`, `np.random.rand`): any library call that draws from it silently shifts every later draw. Never seed from the clock in production code (do it, if at all, once — and *log* the seed). And structure the draws so that adding a new random component does not change existing ones (give each component its own stream), so that a design change can be compared with the old design using **common random numbers** (Section 24.6).

### Versions and environments

The result of a computation is a function of the code, its inputs and the environment — interpreter, libraries, compiler flags, hardware. All three must be pinned and recorded.

- **Version control** (git) for code, configuration and small input data. Every reported result should carry the commit hash of the code that produced it; a "dirty" working tree (uncommitted changes) should be flagged in the output.
- **Environment specification.** A `requirements.txt` with exact versions (`numpy==2.1.3`), a conda `environment.yml` or lock file, or a container image. Library defaults change: tolerances of iterative solvers, the default ordering of a sparse LU, the algorithm behind a random-number distribution. A pinned environment can be rebuilt; "latest" cannot.
- **Metadata with every output.** Write alongside each result file the parameters, the seeds, the commit hash, the library versions and the platform. This costs a few lines:

```python
import json, platform, subprocess, numpy, scipy
meta = {
    "commit": subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip(),
    "numpy": numpy.__version__, "scipy": scipy.__version__,
    "python": platform.python_version(), "machine": platform.machine(),
    "seed": 20241021, "N": 100_000, "tolerance_model": "uniform, 1%",
}
json.dump(meta, open("results/mc_gain_meta.json", "w"), indent=2)
```

- **Scripted pipelines.** Figures and tables should be generated by scripts from raw data with a single command, not by interactive sessions or edits in a spreadsheet. A notebook is acceptable if it is run top to bottom from a clean kernel before results are recorded.

### Testing numerical code

Numerical code fails silently: a wrong sign in a Jacobian slows Newton's method rather than crashing it; an off-by-one in a stencil changes the convergence order from 2 to 1 while still giving "reasonable" answers. Tests must therefore check *mathematical* properties, with tolerances justified by error analysis.

1. **Exact or analytic references.** Compare against cases with known solutions — a symbolic result (Section 24.2), a closed-form circuit response, a problem constructed so that the answer is known (data generated from known parameters, as in Worked Example 24.33).
2. **Convergence-order tests.** For a method of order $p$, the error behaves as $Ch^p$; halving $h$ should divide it by $2^p$. The **observed order** $p_\text{obs}=\log_2\big(e(h)/e(h/2)\big)$ is a far more sensitive check than the error itself.
3. **Method of manufactured solutions** (for ODE/PDE solvers): choose an arbitrary smooth solution $u^*$, substitute it into the equation to compute the source term that makes it exact, run the solver with that source, and verify both the error and its convergence order.
4. **Invariants and properties.** Kirchhoff's current law at every node (residual of $\mathbf{G}\mathbf{v}-\mathbf{i}$), symmetry and positive definiteness of matrices that must have them, reciprocity ($\mathbf{Z}=\mathbf{Z}^\mathsf{T}$ for passive networks), conservation of energy or charge, orthogonality of computed $\mathbf{Q}$. Property-based testing frameworks generate many random inputs and check such properties automatically.
5. **Regression tests.** Store reference outputs and compare future runs against them with a tolerance — catching unintended changes, not establishing correctness.
6. **Statistical tests for stochastic code.** A Monte Carlo estimate is a random variable; test that it lies within $z$ standard errors of the known answer.

> [!example] Worked Example 24.42 — An observed-order test for the trapezoidal rule
> Test a trapezoidal-rule implementation on $\int_0^1e^x\,dx=e-1$ with $n=10$, $20$, $40$ subintervals.
>
> **Errors.** $T_{10}-I=1.4317\times10^{-3}$, $T_{20}-I=3.5796\times10^{-4}$, $T_{40}-I=8.9493\times10^{-5}$.
>
> **Observed order.** $\log_2(1.4317/0.35796)=\log_2(3.9996)=1.9998$, and $\log_2(3.5796/0.89493)=2.0000$ ✓ — the expected $p=2$.
>
> **Leading constant.** The Euler–Maclaurin formula predicts $T_n-I\approx\tfrac{h^2}{12}\big(f'(1)-f'(0)\big)=\tfrac{h^2}{12}(e-1)$, i.e. $1.4319\times10^{-3}$ for $h=0.1$ ✓.
>
> **As a test.** `assert abs(p_obs - 2) < 0.05`. A typical bug — weighting both end points by 1 instead of $\tfrac12$ — gives an error of order $h$ and $p_\text{obs}\approx1$, which this test catches immediately even though $T_{40}$ would still be correct to two digits.

> [!example] Worked Example 24.43 — Setting tolerances for deterministic and stochastic tests
> (a) A regression test checks the solution of a linear system with $\kappa(\mathbf{A})\approx10^6$ computed by LU. What relative tolerance is appropriate? (b) A Monte Carlo yield estimator is tested on a case with known yield $Y=0.80$, using $N=10^4$ samples. What acceptance band gives a false-failure probability below $10^{-4}$?
>
> **(a)** A backward-stable solver guarantees a relative error of about $\kappa\varepsilon\approx10^6\times1.1\times10^{-16}\approx10^{-10}$; differences between platforms or thread counts can be of this size. A tolerance of `rtol=1e-8` (a safety factor of 100) is appropriate; `rtol=1e-15` would fail spuriously, and `rtol=1e-3` would hide real bugs.
>
> **(b)** The standard error is $\sqrt{0.8\times0.2/10^4}=0.004$. Accepting $\lvert\hat Y-0.80\rvert\le4\times0.004=0.016$ has false-failure probability $2Q(4)=6.3\times10^{-5}$ (a 3-standard-error band would fail spuriously $0.27\%$ of the time — about once in every 370 runs of the test suite). With a fixed seed the test is also bitwise repeatable, but the band should still be set statistically, so that the test survives a legitimate change of seed or generator.
>
> **Rule.** Tolerances come from error analysis: $\kappa\varepsilon$ for deterministic linear algebra, the truncation error $Ch^p$ for discretisations, $z\cdot\text{SE}$ for Monte Carlo. A tolerance chosen by trial and error ("1e-6 seemed to work") is a test that has not been designed.

```python
import numpy as np
from numpy.testing import assert_allclose

def test_mc_gain_yield():
    rng = np.random.default_rng(7)
    N = 10_000
    R1 = 1e3 * (1 + 0.01 * rng.uniform(-1, 1, N))
    R2 = 9e3 * (1 + 0.01 * rng.uniform(-1, 1, N))
    G = 1 + R2 / R1
    Y = np.mean(np.abs(G - 10) <= 0.1)
    se = np.sqrt(Y * (1 - Y) / N)
    assert abs(Y - 0.80) < 4 * se          # reference 0.80 from a 10^6-sample run

def test_kcl_residual(G, v, i):
    assert_allclose(G @ v, i, rtol=1e-10, atol=1e-12)
```

> [!warning] Common Pitfall — the reference answer is also a computation
> A regression test that compares a new result with an old one only shows that nothing *changed*. If the reference was produced by the same buggy code, the test enshrines the bug. Independent references — analytic results, a different algorithm (sparse versus dense solve, Monte Carlo versus linearisation), a published benchmark, or a measurement — are what establish correctness. Keep at least one independent check for every major result.

> [!warning] Common Pitfall — "it worked on my machine"
> An analysis run in an interactive session with variables left over from earlier cells, a file edited by hand, a parameter changed and changed back, or an unrecorded library upgrade cannot be reproduced even by its author. The cure is mundane: scripts under version control, a pinned environment, a recorded seed, and a single command that regenerates every number in the report.

> [!tip] Engineering Note
> In circuit design, reproducibility has an extra layer: the simulator version, the model files (PDK release), the simulator options (`reltol`, `abstol`, integration method, Gmin stepping) and the netlist all affect the numbers. Sign-off flows record all of them, and the same Monte Carlo seed is used for "before" and "after" runs of an engineering change so that differences reflect the design, not the sampling. Measurement laboratories accredited to ISO/IEC 17025 must document their uncertainty calculations (Section 24.8) and the software that performs them, which is validated like any other instrument.

## Summary

| Topic | Key formula or fact |
|---|---|
| Asymptotic cost | $f=O(g)$: $f\le cg$ for large $n$; scaling exponent $p\approx\ln(T_2/T_1)/\ln(n_2/n_1)$ |
| Master theorem | $T=aT(n/b)+n^d$: $n^d$, $n^d\log n$ or $n^{\log_ba}$ as $d\gtrless\log_ba$ |
| Typical costs | matrix–vector $2n^2$; matrix–matrix $2n^3$; LU $\tfrac23n^3$; Cholesky $\tfrac13n^3$; FFT $5N\log_2N$; sparse matvec $2\operatorname{nnz}$ |
| Expression swell | symbolic $\det$ has $n!$ terms; exact Newton doubles digits per step |
| Backward error (linear system) | $\eta=\lVert\mathbf{r}\rVert/(\lVert\mathbf{A}\rVert\lVert\hat{\mathbf{x}}\rVert)$; forward error $\lesssim\kappa\eta$ |
| Householder reflector | $\mathbf{H}=\mathbf{I}-2\mathbf{v}\mathbf{v}^\mathsf{T}/\mathbf{v}^\mathsf{T}\mathbf{v}$, $\mathbf{v}=\mathbf{a}+\operatorname{sign}(a_1)\lVert\mathbf{a}\rVert\mathbf{e}_1$ |
| QR cost | $2mn^2-\tfrac23n^3$; backward stable; $\lVert\hat{\mathbf{Q}}^\mathsf{T}\hat{\mathbf{Q}}-\mathbf{I}\rVert=O(\varepsilon)$ |
| Normal equations | $\kappa(\mathbf{A}^\mathsf{T}\mathbf{A})=\kappa(\mathbf{A})^2$ — use QR or SVD |
| Truncated SVD | $\hat{\mathbf{x}}_k=\sum_{i\le k}(\mathbf{u}_i^\mathsf{T}\mathbf{b}/\sigma_i)\mathbf{v}_i$ |
| Stationary iteration | $\mathbf{e}^{(k+1)}=\mathbf{T}\mathbf{e}^{(k)}$; converges iff $\rho(\mathbf{T})<1$; $k_\epsilon\approx\ln\epsilon/\ln\rho$ |
| Model problem | $\rho_J=\cos\pi h$; $\rho_{GS}=\rho_J^2$ |
| Optimal SOR | $\omega_\text{opt}=2/(1+\sqrt{1-\rho_J^2})$, $\rho=\omega_\text{opt}-1$ |
| CG | $\lVert\mathbf{e}_k\rVert_\mathbf{A}\le2\big(\tfrac{\sqrt\kappa-1}{\sqrt\kappa+1}\big)^k\lVert\mathbf{e}_0\rVert_\mathbf{A}$; $k\approx\tfrac12\sqrt\kappa\ln(2/\epsilon)$ |
| GMRES | minimises $\lVert\mathbf{r}_k\rVert_2$ over $\mathcal{K}_k$; exact in $\deg(\text{minimal polynomial})$ steps |
| Preconditioning | solve $\mathbf{P}^{-1}\mathbf{A}\mathbf{x}=\mathbf{P}^{-1}\mathbf{b}$ with $\mathbf{P}\approx\mathbf{A}$, $\mathbf{P}\mathbf{z}=\mathbf{r}$ cheap |
| Sparse storage | CSR/CSC: $2\operatorname{nnz}+n+1$ words |
| Fill, 2-D grid | natural: $n^{3/2}$ fill, $n^2$ flops; nested dissection: $n\log n$ fill, $n^{3/2}$ flops |
| Monte Carlo | $\hat\mu_N=\tfrac1N\sum f(\mathbf{X}_i)$, $\text{SE}=s_N/\sqrt N$ |
| Control variate | $c^*=\operatorname{Cov}(f,h)/\operatorname{Var}(h)$; variance $\times(1-\rho^2)$ |
| Importance sampling | $E_p[f]=E_q[f\,p/q]$ |
| Rare event | relative SE $\approx1/\sqrt{Np}$; $N\approx100/p$ for 10% |
| Yield | $\hat Y=\#\text{pass}/N$, $\text{SE}=\sqrt{\hat Y(1-\hat Y)/N}$ |
| Gauss–Newton | $(\mathbf{J}^\mathsf{T}\mathbf{J})\boldsymbol\delta=\mathbf{J}^\mathsf{T}\mathbf{r}$ |
| Levenberg–Marquardt | $(\mathbf{J}^\mathsf{T}\mathbf{J}+\lambda\operatorname{diag}(\mathbf{J}^\mathsf{T}\mathbf{J}))\boldsymbol\delta=\mathbf{J}^\mathsf{T}\mathbf{r}$ |
| Parameter covariance | $\hat\sigma^2(\mathbf{J}^\mathsf{T}\mathbf{J})^{-1}$, $\hat\sigma^2=\sum r_i^2/(m-n)$ |
| Relative sensitivity | $S^y_p=(p/y)\,\partial y/\partial p$ |
| Complex step | $f'(x)\approx\operatorname{Im}\{f(x+jh)\}/h$, no cancellation |
| Adjoint sensitivity | $\mathbf{A}^\mathsf{T}\boldsymbol\lambda=\mathbf{c}$; $\partial y/\partial p=\boldsymbol\lambda^\mathsf{T}(\partial_p\mathbf{b}-\partial_p\mathbf{A}\,\mathbf{x})$; branch: $-(\lambda_a-\lambda_b)(v_a-v_b)$ |
| GUM propagation | $u_c^2=\sum c_i^2u_i^2+2\sum_{i<j}c_ic_ju_{ij}$; $\mathbf{U}_y=\mathbf{J}\mathbf{U}_x\mathbf{J}^\mathsf{T}$; $U=ku_c$ |
| Uniform tolerance | $\pm a\Rightarrow u=a/\sqrt3$ |
| Observed order | $p_\text{obs}=\log_2\big(e(h)/e(h/2)\big)$ |
| Test tolerances | $\kappa\varepsilon$ (linear algebra), $Ch^p$ (discretisation), $z\cdot\text{SE}$ (Monte Carlo) |

## Further reading

- L. N. Trefethen and D. Bau III, *Numerical Linear Algebra*, 25th anniversary ed., SIAM, 2022 — QR, SVD, conditioning, backward stability, Arnoldi, GMRES and CG, in short, beautifully argued lectures.
- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins University Press, 2013 — the reference for algorithms, flop counts and implementation details.
- Y. Saad, *Iterative Methods for Sparse Linear Systems*, 2nd ed., SIAM, 2003 — Krylov methods, preconditioning and sparse storage; T. A. Davis, *Direct Methods for Sparse Linear Systems*, SIAM, 2006, is the companion for sparse factorisation and ordering.
- N. J. Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd ed., SIAM, 2002 — backward error analysis, summation, and testing numerical software.
- A. B. Owen, *Monte Carlo Theory, Methods and Examples*, online text, 2013– — variance reduction, importance sampling and quasi-Monte Carlo; with JCGM 100:2008 (*GUM*) and JCGM 101:2008 (*GUM Supplement 1*) for uncertainty propagation, freely available from the BIPM.

---

[← Previous](23-optimization.md) · [Table of Contents](../00-front/table-of-contents.md) · [Next →](25-mathematical-modeling.md)
