# Chapter 22 — Numerical Methods

[← Previous](../part-6-probability-statistics/21-random-processes.md) · [Table of Contents](../00-front/table-of-contents.md) · [Next →](23-optimization.md)

> **Prerequisites:** Ch. 4 (Differential Calculus — Taylor's theorem), Ch. 5 (Integral Calculus), Ch. 8 (Matrices and Linear Systems), Ch. 10 (Eigenvalues and Eigenvectors), Ch. 11–13 (ODEs and PDEs)
> **Learning outcomes:** After this chapter you can …
> - explain how IEEE-754 floating-point numbers are stored, bound the rounding error of an arithmetic operation, and recognise catastrophic cancellation;
> - separate truncation error from rounding error, and distinguish the *conditioning* of a problem from the *stability* of an algorithm;
> - solve nonlinear equations by bisection, Newton–Raphson, secant and fixed-point iteration, and predict their order of convergence and failure modes;
> - build interpolating polynomials (Lagrange, Newton) and cubic splines, explain the Runge phenomenon and cure it with Chebyshev nodes, and compute least-squares polynomial approximations;
> - derive finite-difference and quadrature rules with their error orders, and accelerate them by Richardson extrapolation;
> - integrate ODEs with Euler and Runge–Kutta methods, recognise stiffness, and discretise the heat, wave and Laplace equations with step sizes that respect stability (CFL) limits;
> - state the consistency–stability–convergence relationship (Lax equivalence) and carry out a von Neumann stability analysis.

Almost every quantitative result an engineer produces today passes through a computer. SPICE solves the nonlinear equations of a circuit by Newton's method; an electromagnetic field solver discretises Maxwell's equations on a grid; a control engineer simulates a plant by integrating ODEs with an adaptive Runge–Kutta method; a DSP engineer worries about whether a 16-bit accumulator will overflow. In each case the answer on the screen is *not* the exact mathematical answer — it is an approximation, and its quality depends on choices that the user often leaves at default values without understanding them.

This chapter is about those choices. The guiding questions are always the same three:

1. **What is computed?** (the algorithm, i.e. the discrete approximation to a continuous problem);
2. **How wrong is it?** (truncation error from the discretisation, rounding error from finite-precision arithmetic, and how both propagate);
3. **How does the error behave as we refine?** (order of convergence, stability, cost).

Throughout, small Python snippets show the idea in executable form; every number quoted in a worked example was checked by running such code in IEEE-754 double precision.

## 22.1 Floating-point arithmetic

### Why an engineer must care

Real numbers have infinitely many digits; computer words have 16, 32 or 64 bits. Every real number that enters a computation is therefore replaced by a nearby *representable* number, and every arithmetic operation on representable numbers is followed by a rounding back to the representable set. Usually the resulting errors are around $10^{-16}$ relative and harmless. Occasionally they are amplified into total nonsense — a subtraction of nearly equal numbers, a recursion that multiplies errors by $n$ at every step, a matrix that is nearly singular. To tell the two situations apart you need a precise model of floating-point arithmetic.

### Floating-point number systems

> [!definition] Floating-point number system
> A (normalised, binary) floating-point number system $\mathbb{F}(p, e_{\min}, e_{\max})$ consists of zero and all numbers of the form
> $$
> x = \pm\,(1.d_1 d_2 \ldots d_{p-1})_2 \times 2^{e},\qquad d_i\in\{0,1\},\quad e_{\min}\le e\le e_{\max}.
> $$
> The integer $p$ is the **precision** (number of significant bits, including the leading 1), $(1.d_1\ldots d_{p-1})_2$ is the **significand** (or mantissa) and $e$ the **exponent**.

Two features follow at once.

- **The numbers are not evenly spaced.** Between $2^e$ and $2^{e+1}$ there are exactly $2^{p-1}$ equally spaced numbers with spacing $2^{e-p+1}$. Each time you cross a power of two the spacing doubles. Relative spacing, however, is roughly constant: this is the whole point of floating point.
- **There is a largest and a smallest positive normalised number**, $x_{\max} = (2-2^{1-p})\,2^{e_{\max}}$ and $x_{\min} = 2^{e_{\min}}$. Results larger than $x_{\max}$ **overflow**; results smaller than $x_{\min}$ **underflow**.

> [!definition] Machine epsilon, unit roundoff, ulp
> - **Machine epsilon** $\varepsilon_M = 2^{1-p}$ is the distance from 1 to the next larger floating-point number.
> - **Unit roundoff** $u = \tfrac12\varepsilon_M = 2^{-p}$ is the maximum relative error committed when a real number in the normal range is rounded to nearest.
> - The **ulp** (unit in the last place) of $x$ is the spacing of floating-point numbers near $x$: $\operatorname{ulp}(x) = 2^{e-p+1}$ when $2^e\le\lvert x\rvert<2^{e+1}$.

Be aware that some texts (and MATLAB's `eps`) call $\varepsilon_M$ "machine epsilon" while others reserve that name for $u$. The factor of two rarely matters, but always check which one a bound uses.

### The IEEE-754 standard

Since 1985 essentially all hardware has implemented the IEEE-754 standard (revised 2008 and 2019). The formats an EE meets most often are:

| Format | Total bits | Sign / exponent / fraction bits | Precision $p$ | $\varepsilon_M=2^{1-p}$ | Largest normal | Smallest normal |
|---|---|---|---|---|---|---|
| binary16 (half) | 16 | 1 / 5 / 10 | 11 | $9.77\times10^{-4}$ | $65504$ | $6.10\times10^{-5}$ |
| binary32 (single, `float`) | 32 | 1 / 8 / 23 | 24 | $1.19\times10^{-7}$ | $3.40\times10^{38}$ | $1.18\times10^{-38}$ |
| binary64 (double) | 64 | 1 / 11 / 52 | 53 | $2.22\times10^{-16}$ | $1.80\times10^{308}$ | $2.23\times10^{-308}$ |

So single precision carries about 7 significant decimal digits and double precision about 16.

**Encoding.** A binary64 word is laid out as

$$
\underbrace{s}_{1\ \text{bit}}\ \underbrace{E}_{11\ \text{bits}}\ \underbrace{f}_{52\ \text{bits}},\qquad
x = (-1)^s\,(1.f)_2\times 2^{E-1023}\quad (1\le E\le 2046).
$$

The exponent is stored with a **bias** of $1023$ ($127$ for single), so that the stored field $E$ is an unsigned integer. The leading 1 of the significand is *not* stored ("hidden bit"), which is why 52 stored fraction bits give $p=53$. The two extreme exponent fields are reserved:

| Stored $E$ | $f$ | Meaning |
|---|---|---|
| $0$ | $0$ | $\pm 0$ (signed zero) |
| $0$ | $\neq 0$ | **subnormal**: $(-1)^s (0.f)_2 \times 2^{-1022}$ |
| $1\ldots2046$ | any | normal number |
| $2047$ | $0$ | $\pm\infty$ |
| $2047$ | $\neq 0$ | NaN (Not a Number) |

Subnormals fill the gap between $0$ and $x_{\min}$ evenly ("gradual underflow"); the smallest positive subnormal double is $2^{-1074}\approx 4.94\times10^{-324}$. Infinities and NaNs make arithmetic *closed*: $1/0=+\infty$, $0/0=\text{NaN}$, $\sqrt{-1}=\text{NaN}$, and any comparison involving NaN is false (so `x != x` is the classic NaN test).

> [!example] Worked Example 22.1 — Encoding a number in binary64
> Encode $x = 5.75$.
>
> **Step 1 — binary expansion.** $5 = 101_2$ and $0.75 = 0.11_2$, so $x = 101.11_2$.
>
> **Step 2 — normalise.** $101.11_2 = 1.0111_2\times 2^{2}$, so $s=0$, $e=2$, fraction $f = 0111\,000\ldots0$ (52 bits).
>
> **Step 3 — bias the exponent.** $E = 2 + 1023 = 1025 = 100\,0000\,0001_2$.
>
> **Step 4 — assemble.** $0\;10000000001\;0111\,0000\ldots$ Grouping in fours: $0100\ 0000\ 0001\ 0111\ 0000\ldots$, i.e. the hexadecimal word `0x4017000000000000`.
>
> **Check.** In Python, `struct.pack('>d', 5.75).hex()` returns `'4017000000000000'`. Because $5.75$ has a finite binary expansion, it is stored *exactly*.

> [!example] Worked Example 22.2 — Why $0.1+0.2\neq0.3$
> In binary, $0.1 = 0.0\overline{0011}_2$ repeats forever, so it must be rounded. The stored double nearest to $0.1$ is
> $$
> \operatorname{fl}(0.1) = 0.1000000000000000055511151231257827\ldots
> $$
> (hex `0x3FB999999999999A` — note the final `A`: the repeating `9` pattern rounded up). The absolute error is $5.55\times10^{-18}$, so the relative error is $5.55\times10^{-18}/0.1\approx 5.6\times10^{-17}<u = 1.11\times10^{-16}$, as promised.
>
> Adding $\operatorname{fl}(0.1)+\operatorname{fl}(0.2)$ and rounding gives `0.30000000000000004`, which is one ulp above $\operatorname{fl}(0.3)$. Likewise, adding $0.1$ ten times in a loop gives `0.9999999999999999`. Neither result is a bug; each is the correctly rounded result of operations on numbers that were never exactly $0.1$, $0.2$, $0.3$.
>
> **Lesson.** Never test floating-point results for exact equality. Use a tolerance, e.g. `abs(a-b) <= rtol*max(abs(a),abs(b)) + atol`.

### Rounding and the standard model

IEEE-754 requires the five basic operations $+,-,\times,\div,\sqrt{\ }$ to be **correctly rounded**: the result is the exact real result, rounded to the representable set according to the current rounding mode (default: *round to nearest, ties to even*). This single requirement gives the most useful model in numerical analysis.

> [!theorem] Standard model of floating-point arithmetic
> Let $\operatorname{fl}(\cdot)$ denote rounding to nearest, and let $u=2^{-p}$. If $x$ lies in the normal range, then
> $$
> \operatorname{fl}(x) = x(1+\delta),\qquad \lvert\delta\rvert\le u.
> $$
> Consequently, for floating-point numbers $a,b$ and $\circ\in\{+,-,\times,\div\}$ (no overflow or underflow),
> $$
> \operatorname{fl}(a\circ b) = (a\circ b)(1+\delta),\qquad \lvert\delta\rvert\le u.
> $$

*Proof sketch.* Suppose $2^e\le x<2^{e+1}$. Representable numbers in this interval are spaced $\Delta = 2^{e-p+1}$ apart, so the nearest one is within $\Delta/2 = 2^{e-p}$ of $x$. Dividing by $x\ge 2^e$ gives a relative error at most $2^{-p}=u$. Correct rounding means the computed $a\circ b$ is just $\operatorname{fl}$ applied to the exact result. $\square$

The model says that *each individual operation* is accurate to about 16 digits in double precision. Large errors can only arise from the way errors from many operations accumulate or are amplified, which is the subject of Section 22.2.

### Consequences: non-associativity and summation error

Floating-point addition is commutative but **not associative**. In double precision, $(10^{16} + 1) - 10^{16} = 0$ because $10^{16}+1$ rounds to $10^{16}$ (the ulp at $10^{16}$ is 2), whereas $10^{16} + (1 - 10^{16}) = 1$. The order of operations matters.

For a sum $s_n = x_1+\cdots+x_n$ computed left to right, repeated application of the standard model gives

$$
\lvert \hat s_n - s_n\rvert \;\le\; (n-1)\,u\sum_{i=1}^{n}\lvert x_i\rvert + O(u^2).
$$

The bound grows linearly in $n$; in practice, with randomly signed rounding errors, the error typically grows like $\sqrt{n}\,u$. For a million-term sum in single precision ($u\approx 6\times10^{-8}$) the worst-case bound is already of order $6\times 10^{-2}$ relative — a serious problem for, say, accumulating energy in a long FFT or a power meter. Two cures:

- **Pairwise summation** (add in a binary tree): error bound $O(u\log_2 n)$. NumPy's `sum` does this.
- **Compensated (Kahan) summation**: carry a correction term that recaptures the low-order bits lost in each addition; the error bound becomes $2u\sum\lvert x_i\rvert + O(nu^2)$, essentially independent of $n$.

```python
def kahan_sum(xs):
    s, c = 0.0, 0.0           # running sum and compensation
    for x in xs:
        y = x - c             # subtract the error from last time
        t = s + y             # low bits of y may be lost here ...
        c = (t - s) - y       # ... and are recovered algebraically
        s = t
    return s
```

### Catastrophic cancellation

> [!definition] Cancellation
> When two nearly equal numbers $a\approx b$ that already carry relative errors are subtracted, the leading digits cancel and the relative error of $a-b$ is magnified by the factor $\lvert a\rvert/\lvert a-b\rvert$. The subtraction itself is exact (Sterbenz's lemma: if $b/2\le a\le 2b$ then $a-b$ is exactly representable); it merely *exposes* the errors that were already present.

> [!example] Worked Example 22.3 — Evaluating $1-\cos x$ for small $x$
> Compute $g(x) = (1-\cos x)/x^2$ at $x = 1.2\times10^{-5}$. The exact value is very close to $\tfrac12$ (Taylor: $\tfrac12 - x^2/24$).
>
> **Naive formula.** $\cos x = 1 - 7.2\times 10^{-11}+\ldots$ is computed with absolute error up to about $u\approx 1.1\times10^{-16}$. The subtraction $1-\cos x\approx 7.2\times10^{-11}$ then has relative error up to $1.1\times10^{-16}/7.2\times10^{-11}\approx 1.5\times10^{-6}$. In double precision the naive formula returns $0.4999997330$ — only six correct digits.
>
> **Rewritten formula.** Use the identity $1-\cos x = 2\sin^2(x/2)$, which involves no subtraction. This returns $0.49999999999400$, correct to full precision.
>
> At $x = 10^{-8}$ the naive formula gives $1-\cos x = 0$ exactly (because $\cos x$ rounds to 1), a 100 % error, while $2\sin^2(x/2) = 5.0\times10^{-17}$ is correct.

> [!example] Worked Example 22.4 — A stable quadratic formula
> Solve $x^2 - 10^8 x + 1 = 0$. The roots are $x_1\approx 10^8$ and $x_2\approx 10^{-8}$.
>
> **Textbook formula.** $x_{2} = \big(10^8 - \sqrt{10^{16}-4}\big)/2$. Here $\sqrt{10^{16}-4}$ rounds to a number within one ulp ($\approx 1.5\times10^{-8}$) of $10^8$, and subtracting it from $10^8$ cancels all the good digits. Double precision returns $x_2 = 7.45\times10^{-9}$, a 25 % error.
>
> **Stable approach.** Compute the large-magnitude root without cancellation, $x_1 = \big(-b + \operatorname{sign}(-b)\sqrt{b^2-4ac}\big)/(2a) = 10^8$, then use Vieta's relation $x_1x_2 = c/a$: $x_2 = c/(a x_1) = 10^{-8}$, correct to full precision.
>
> The general rule: **rearrange formulas so that nearly equal quantities that carry error are never subtracted.**

### Overflow, underflow and scaling

Intermediate results can overflow even when the final answer is modest. Computing the magnitude of a phasor $\lvert a+jb\rvert=\sqrt{a^2+b^2}$ with $a = b = 10^{200}$ overflows at $a^2$, although the answer $1.41\times10^{200}$ is representable. Library functions such as `hypot(a,b)` scale first: with $m=\max(\lvert a\rvert,\lvert b\rvert)$, $\sqrt{a^2+b^2} = m\sqrt{(a/m)^2+(b/m)^2}$. Similarly, probabilities of long sequences (Ch. 19) and likelihoods underflow, which is why they are accumulated as logarithms.

> [!warning] Common Pitfall — Comparing floats and counting with floats
> Loops such as `t = 0.0; while t != 1.0: t += 0.1` never terminate, because $t$ steps over $1.0$. Loop over an integer counter and compute $t = n\,h$ from it. The same applies to time stepping in simulations: computing $t_n = n h$ avoids the drift that $t \mathrel{+}= h$ accumulates.

> [!tip] Engineering Note — Fixed point versus floating point in DSP
> Many low-cost DSPs and FPGAs use **fixed-point** arithmetic, e.g. the Q15 format, a 16-bit two's-complement integer $k$ representing $k/2^{15}\in[-1,1)$. The spacing is uniform ($2^{-15}\approx 3.05\times10^{-5}$), so the *absolute* error is constant but the *relative* error of small signals is large, and sums can overflow. Designers therefore scale filter coefficients and signals carefully and analyse quantisation noise statistically (each rounding is modelled as additive uniform noise of variance $\Delta^2/12$ — see Ch. 19–21). Floating point removes most scaling headaches at the cost of silicon area and power; half precision (binary16) and bfloat16 are now common in machine-learning accelerators, where $\varepsilon_M\approx10^{-3}$ is acceptable.

## 22.2 Numerical error

### Measuring error

> [!definition] Absolute and relative error
> If $\hat x$ approximates $x$, the **absolute error** is $\lvert \hat x - x\rvert$ and, for $x\neq0$, the **relative error** is
> $$
> \rho(\hat x) = \frac{\lvert\hat x - x\rvert}{\lvert x\rvert}.
> $$
> For vectors replace $\lvert\cdot\rvert$ by a norm. $\hat x$ has about $k$ **correct significant digits** when $\rho\lesssim 0.5\times10^{-k}$.

Relative error is dimensionless and scale-invariant, which is why it is the natural measure for floating-point results ($\rho\le u$ after a single rounding). Absolute error is the right measure near zero, where relative error is meaningless: a computed voltage of $3\times10^{-17}$ V where the exact answer is $0$ has an infinite relative error but is perfectly acceptable.

### Sources of error

Any numerical answer is contaminated by four distinct kinds of error, and it pays to keep them apart:

1. **Modelling error** — the equations idealise reality (an ideal diode, a lumped-element line, a linearised amplifier). Not a numerical issue, but often the largest error (Ch. 25).
2. **Data error** — component tolerances, measurement noise, constants known to finite accuracy.
3. **Truncation (discretisation) error** — the error of replacing a limiting process by a finite one: truncating a Taylor series, replacing a derivative by a difference quotient, an integral by a finite sum, an infinite iteration by a finite number of steps. It exists even in exact arithmetic and typically behaves like $C h^p$ for a step $h$.
4. **Rounding error** — from finite-precision arithmetic (Section 22.1). It exists even for an exact algorithm and typically behaves like $C' u/h^q$ when it is amplified by small steps.

> [!example] Worked Example 22.5 — Truncation error of a Taylor series
> How many terms of $e^x = \sum_k x^k/k!$ are needed to compute $e$ with error below $10^{-6}$?
>
> **Step 1 — remainder.** By Taylor's theorem with Lagrange remainder (Ch. 4), truncating after the $x^n$ term leaves
> $$
> R_n(1) = \frac{e^{\xi}}{(n+1)!},\qquad 0<\xi<1,\qquad\text{so}\quad 0< R_n(1) < \frac{e}{(n+1)!}.
> $$
>
> **Step 2 — solve the inequality.** We need $(n+1)! > e\times10^{6}\approx 2.72\times10^6$. Since $9! = 362\,880$ and $10! = 3\,628\,800$, $n+1 = 10$, i.e. $n = 9$ (ten terms, $k=0,\ldots,9$).
>
> **Check.** $\sum_{k=0}^{9}1/k! = 2.718281526$, and $e - 2.718281526 = 3.03\times10^{-7}<10^{-6}$. The truncation error is a property of the *formula*; rounding contributes only $\sim10^{-16}$ here.

### Order notation

We write $E(h) = O(h^p)$ as $h\to0$ if $\lvert E(h)\rvert\le C h^p$ for some constant $C$ and all sufficiently small $h$. A method with error $O(h^p)$ is said to be **of order $p$**. The order tells you how the error responds to refinement: halving $h$ divides the error by about $2^p$. Conversely, from errors $E_1,E_2$ at steps $h_1, h_2$ you can *measure* the order:

$$
p \approx \frac{\ln(E_1/E_2)}{\ln(h_1/h_2)}.
$$

This "observed order" test is the single most useful debugging tool for numerical code: if a method that should be second order shows $p\approx1$, there is a bug (often in a boundary condition).

For iterations $x_k\to x^\ast$ the analogous notion is the order of convergence $q$: $\lvert x_{k+1}-x^\ast\rvert\approx C\lvert x_k - x^\ast\rvert^{q}$ (Section 22.3).

### Propagation of errors

If $y = f(x_1,\ldots,x_n)$ and each input carries a small error $\Delta x_i$, a first-order Taylor expansion (Ch. 6) gives

$$
\Delta y \approx \sum_{i=1}^{n}\frac{\partial f}{\partial x_i}\,\Delta x_i .
$$

Bounding term by term gives the **worst-case** estimate $\lvert\Delta y\rvert\le\sum_i\lvert\partial f/\partial x_i\rvert\,\lvert\Delta x_i\rvert$. If the errors are independent random variables with standard deviations $\sigma_i$, the variance formula (Ch. 19–20) gives the **root-sum-square** estimate $\sigma_y^2\approx\sum_i(\partial f/\partial x_i)^2\sigma_i^2$. For products and quotients it is convenient to work with relative errors: if $y = x_1^{a_1}x_2^{a_2}\cdots$, then $\Delta y/y \approx \sum_i a_i\,\Delta x_i/x_i$.

> [!example] Worked Example 22.6 — Error in a computed power
> A resistor $R = 50\ \Omega\pm2\%$ carries a measured voltage $V = 10\ \mathrm{V}\pm1\%$. Estimate the error in $P = V^2/R$.
>
> **Nominal value.** $P = 100/50 = 2\ \mathrm{W}$.
>
> **Relative sensitivities.** $P = V^{2}R^{-1}$, so $\Delta P/P\approx 2\,\Delta V/V - \Delta R/R$.
>
> **Worst case.** $2(1\%)+1(2\%) = 4\%$, i.e. $P = 2.00\pm0.08\ \mathrm{W}$.
>
> **Statistical (independent tolerances treated as standard deviations).** $\sqrt{(2\cdot1\%)^2+(2\%)^2} = 2.83\%$, i.e. $\pm0.057\ \mathrm{W}$.
>
> The squared voltage doubles its relative error — the exponent in a power law is exactly the condition number defined next.

### Conditioning of a problem

Error propagation describes how the *problem* responds to perturbations of its data, independently of any algorithm. This intrinsic sensitivity is captured by a condition number.

> [!definition] Condition number of a scalar function
> The (relative) condition number of evaluating a differentiable $f$ at $x$ is
> $$
> \kappa_f(x) = \lim_{\varepsilon\to0}\ \sup_{\lvert\Delta x\rvert\le\varepsilon\lvert x\rvert}\frac{\lvert f(x+\Delta x)-f(x)\rvert/\lvert f(x)\rvert}{\lvert\Delta x\rvert/\lvert x\rvert} = \left\lvert\frac{x\,f'(x)}{f(x)}\right\rvert.
> $$
> The problem is **well-conditioned** if $\kappa$ is modest and **ill-conditioned** if $\kappa$ is large. A relative data error $\rho$ produces a relative result error of about $\kappa\rho$.

> [!example] Worked Example 22.7 — Condition numbers of elementary operations
> | Problem | $\kappa_f(x)$ | Comment |
> |---|---|---|
> | $f(x)=\sqrt{x}$ | $\tfrac12$ | always well-conditioned |
> | $f(x)=x^n$ | $\lvert n\rvert$ | squaring doubles relative error (Example 22.6) |
> | $f(x)=e^x$ | $\lvert x\rvert$ | ill-conditioned for large $\lvert x\rvert$: $e^{700}$ is sensitive |
> | $f(x)=\ln x$ | $1/\lvert\ln x\rvert$ | ill-conditioned near $x=1$ |
> | $f(x)=x-1$ | $\lvert x/(x-1)\rvert$ | at $x = 1.0001$: $\kappa = 10\,001$ |
> | $f(x)=\tan x$ | $\lvert x/(\sin x\cos x)\rvert$ | ill-conditioned near $x = \pi/2$ |
>
> The entry for $x-1$ is cancellation (Section 22.1) seen as a *conditioning* statement: subtracting nearly equal numbers is an ill-conditioned problem, so no algorithm can save you if $x$ itself is inexact.

For a linear system $\mathbf{A}\mathbf{x} = \mathbf{b}$ the corresponding quantity is the **matrix condition number** (Ch. 8)

$$
\kappa(\mathbf{A}) = \lVert\mathbf{A}\rVert\,\lVert\mathbf{A}^{-1}\rVert = \frac{\sigma_{\max}}{\sigma_{\min}}\quad(\text{2-norm}),
$$

and perturbation theory gives, for $\mathbf{A}(\mathbf{x}+\Delta\mathbf{x}) = \mathbf{b}+\Delta\mathbf{b}$,

$$
\frac{\lVert\Delta\mathbf{x}\rVert}{\lVert\mathbf{x}\rVert}\le\kappa(\mathbf{A})\,\frac{\lVert\Delta\mathbf{b}\rVert}{\lVert\mathbf{b}\rVert}.
$$

*Derivation.* $\Delta\mathbf{x} = \mathbf{A}^{-1}\Delta\mathbf{b}$ gives $\lVert\Delta\mathbf{x}\rVert\le\lVert\mathbf{A}^{-1}\rVert\lVert\Delta\mathbf{b}\rVert$, and $\mathbf{b} = \mathbf{A}\mathbf{x}$ gives $\lVert\mathbf{b}\rVert\le\lVert\mathbf{A}\rVert\lVert\mathbf{x}\rVert$, i.e. $1/\lVert\mathbf{x}\rVert\le\lVert\mathbf{A}\rVert/\lVert\mathbf{b}\rVert$. Multiply the two inequalities. $\square$

A rule of thumb follows: **solving a linear system in floating point loses about $\log_{10}\kappa(\mathbf{A})$ decimal digits.**

> [!example] Worked Example 22.8 — An ill-conditioned circuit matrix
> Consider
> $$
> \mathbf{A}=\begin{bmatrix}1&1\\1&1.0001\end{bmatrix},\qquad \mathbf{b}=\begin{bmatrix}2\\2.0001\end{bmatrix}\ \Rightarrow\ \mathbf{x}=\begin{bmatrix}1\\1\end{bmatrix}.
> $$
> Such matrices arise when two nodes are joined by a very small resistance compared with everything else, so their equations are nearly identical.
>
> **Condition number.** $\kappa_2(\mathbf{A})\approx 4.0\times10^{4}$.
>
> **Perturb the data.** Change $b_2$ by $10^{-4}$ to $2.0002$ — a relative change $\lVert\Delta\mathbf b\rVert/\lVert\mathbf b\rVert = 10^{-4}/2.83 = 3.5\times10^{-5}$. The solution becomes $\mathbf{x} = (0,2)^\mathsf{T}$: a relative change of $\lVert(-1,1)\rVert/\lVert(1,1)\rVert = 100\%$.
>
> **Compare with the bound.** $\kappa\cdot3.5\times10^{-5}\approx1.4$, so the observed amplification ($\approx 2.8\times10^4$) is close to the worst case. Nothing is wrong with the solver — the *problem* is ill-conditioned.
>
> **Hilbert matrices** $H_{ij} = 1/(i+j-1)$, which appear in least-squares polynomial fitting (Section 22.6), are a notorious family: $\kappa_2(\mathbf H_3)\approx 524$, $\kappa_2(\mathbf H_5)\approx4.8\times10^5$, $\kappa_2(\mathbf H_{8})\approx1.5\times10^{10}$, $\kappa_2(\mathbf H_{12})\approx 1.6\times10^{16}$ — at $n=12$ double precision has no correct digits left.

### Stability of an algorithm

Conditioning belongs to the problem; **stability** belongs to the algorithm. Let $f$ be the exact map from data to result and $\hat f$ the algorithm executed in floating point.

> [!definition] Forward error, backward error, backward stability
> - The **forward error** is $\lvert\hat f(x) - f(x)\rvert$ (relative: divided by $\lvert f(x)\rvert$).
> - The **backward error** is the smallest $\lvert\Delta x\rvert/\lvert x\rvert$ such that $\hat f(x) = f(x+\Delta x)$: how much would the data have to change for the computed answer to be *exact*?
> - An algorithm is **backward stable** if its backward error is $O(u)$ for all inputs. It is (forward) **unstable** if it can produce forward errors much larger than $\kappa_f(x)\,u$.

These ideas combine in the most important inequality in numerical analysis:

$$
\text{forward error}\ \lesssim\ \text{condition number}\times\text{backward error}.
$$

A backward-stable algorithm therefore delivers the best accuracy the problem allows: about $\kappa\,u$. Gaussian elimination with partial pivoting, Householder QR and the standard algorithms for $\sqrt{\ }$, $\sin$, the FFT, etc. are backward stable (or nearly so). When a computed answer is bad, ask first: *is the problem ill-conditioned, or is my algorithm unstable?* The cures are completely different — reformulate the problem (or get better data) in the first case, change the algorithm in the second.

> [!example] Worked Example 22.9 — An unstable recurrence and its stable reversal
> Let $I_n = \int_0^1 x^n e^{x-1}\,dx$. Integration by parts gives the exact recurrence
> $$
> I_n = 1 - n\,I_{n-1},\qquad I_0 = 1 - e^{-1}\approx 0.632121.
> $$
> Since $0<x^ne^{x-1}\le x^n$, we know $0<I_n<1/(n+1)$, so $I_n\to0$ slowly.
>
> **Forward recurrence in double precision.**
>
> | $n$ | 0 | 1 | 5 | 10 | 15 | 17 | 18 | 19 | 20 |
> |---|---|---|---|---|---|---|---|---|---|
> | $\hat I_n$ | 0.632121 | 0.367879 | 0.145533 | 0.083877 | 0.059034 | 0.057192 | $-0.029454$ | 1.55962 | $-30.1924$ |
>
> By $n = 17$ the value has started to *increase* and at $n = 18$ it is negative, both impossible.
>
> **Explanation.** An error $\epsilon_0$ in $I_0$ (about $u\cdot I_0\approx 7\times10^{-17}$ from rounding) propagates as $\epsilon_n = -n\,\epsilon_{n-1}$, so $\lvert\epsilon_n\rvert = n!\,\lvert\epsilon_0\rvert$. With $20!\approx 2.4\times10^{18}$, $\lvert\epsilon_{20}\rvert\sim 10^{2}$ — consistent in order of magnitude with the observed $-30$. The *problem* (computing $I_{20}$) is perfectly well conditioned; the *algorithm* is unstable.
>
> **Stable reversal.** Run the recurrence backwards, $I_{n-1} = (1-I_n)/n$. Now errors are *divided* by $n$ at every step. Starting from the crude guess $I_{30}=0$ (error $<1/31$) and recurring down, the error at $n=20$ is reduced by $30!/20!\approx 1.1\times10^{14}$ and is below $10^{-15}$. This gives $I_{20}=0.0455449$, $I_{10}=0.0838771$, and $I_0 = 0.6321205588285577$, which agrees with $1-e^{-1}$ to all 16 digits.

> [!warning] Common Pitfall — Blaming the computer for an ill-conditioned problem
> Using higher precision (quad, arbitrary precision) cures *instability* only by brute force and does nothing for a problem whose *data* are uncertain. If the entries of $\mathbf b$ in Example 22.8 are known only to $10^{-4}$, the solution is unknown to 100 % no matter how the arithmetic is done. Conversely, a large residual from a backward-stable solver is a sign of a bug, not of ill-conditioning: backward-stable methods always produce small residuals $\lVert\mathbf b - \mathbf A\hat{\mathbf x}\rVert\approx u\lVert\mathbf A\rVert\lVert\hat{\mathbf x}\rVert$ even when $\hat{\mathbf x}$ is inaccurate.

> [!tip] Engineering Note — Ill-conditioning in practice
> Ill-conditioned problems are common in EE: polynomial root-finding from coefficients (Wilkinson showed that changing the $x^{19}$ coefficient of $\prod_{k=1}^{20}(x-k)$ by $2^{-23}$ moves some roots far into the complex plane), high-order IIR filters implemented in direct form (pole positions are ill-conditioned functions of the polynomial coefficients — hence the universal use of cascaded second-order sections), deconvolution and inverse problems (which need regularisation), and nodal analysis with resistances spanning many decades. Recognising ill-conditioning — and reformulating to avoid it — is a design skill.

## 22.3 Root-finding methods

### The problem

Given a continuous $f:\mathbb{R}\to\mathbb{R}$, find $x^\ast$ with $f(x^\ast) = 0$. Nonlinear equations without closed-form solutions appear constantly in engineering: the operating point of a diode or transistor circuit, the resonant frequency of a distributed structure ($\tan\beta\ell = \ldots$), the propagation constant of a waveguide mode, the time at which a transient crosses a threshold, Kepler-type equations in orbital mechanics, the break frequency at which a gain equals $-3$ dB. All practical methods are **iterative**: they generate a sequence $x_0, x_1, x_2,\ldots$ that (one hopes) converges to $x^\ast$.

Two questions organise the subject: **Does the iteration converge** (from which starting points)? And **how fast**?

### Existence and bracketing

> [!theorem] Intermediate value theorem (root form)
> If $f$ is continuous on $[a,b]$ and $f(a)f(b)<0$, then $f$ has at least one root in $(a,b)$.

An interval with a sign change is called a **bracket**. Bracketing methods (bisection, regula falsi, Brent) keep a bracket at every step and therefore *cannot fail* to converge to a root; open methods (Newton, secant, fixed-point) use only the current iterates, converge faster near the root, but can diverge from a poor start. Robust library routines (`scipy.optimize.brentq`, MATLAB `fzero`) combine both: open steps when they behave, bisection steps when they do not.

### Order of convergence

> [!definition] Order of convergence
> Let $x_k\to x^\ast$ and $e_k = x_k - x^\ast$. The sequence converges with **order** $q\ge1$ and **asymptotic constant** $C$ if
> $$
> \lim_{k\to\infty}\frac{\lvert e_{k+1}\rvert}{\lvert e_k\rvert^{q}} = C,\qquad 0<C<\infty,
> $$
> with $C<1$ required when $q = 1$. The cases $q=1$, $1<q<2$ and $q=2$ are called **linear**, **superlinear** and **quadratic** convergence.

The practical meaning: with linear convergence the number of correct digits grows by a *fixed amount* $-\log_{10}C$ per iteration; with order $q$ the number of correct digits is *multiplied by* $q$ per iteration. Quadratic convergence roughly doubles the digits each step, so from 2 correct digits one reaches full double precision in about three steps.

### Fixed-point iteration

Rewrite $f(x)=0$ in the equivalent form $x = g(x)$ and iterate

$$
x_{k+1} = g(x_k).
$$

A **fixed point** of $g$ is a point with $g(x^\ast) = x^\ast$. There are many ways to rewrite one equation: for $x^2 - 2 = 0$ one could use $g(x) = 2/x$, $g(x) = x - (x^2-2)/3$, or $g(x) = \tfrac12(x + 2/x)$ — and they behave very differently.

> [!theorem] Contraction mapping (Banach fixed-point) theorem, scalar form
> Suppose $g$ maps the closed interval $I=[a,b]$ into itself and is a contraction there:
> $$
> \lvert g(x) - g(y)\rvert\le L\lvert x-y\rvert\quad\text{for all } x,y\in I,\qquad 0\le L<1
> $$
> (sufficient: $g$ differentiable with $\lvert g'(x)\rvert\le L<1$ on $I$). Then
> 1. $g$ has exactly one fixed point $x^\ast\in I$;
> 2. for every $x_0\in I$ the iteration $x_{k+1} = g(x_k)$ converges to $x^\ast$;
> 3. the error satisfies $\lvert x_k - x^\ast\rvert\le L^k\lvert x_0 - x^\ast\rvert$ and the computable bound $\lvert x_k - x^\ast\rvert\le\dfrac{L}{1-L}\lvert x_k - x_{k-1}\rvert$.

*Proof sketch.* $\lvert x_{k+1}-x_k\rvert = \lvert g(x_k) - g(x_{k-1})\rvert\le L\lvert x_k - x_{k-1}\rvert$, so the step lengths decrease geometrically and $\sum\lvert x_{k+1}-x_k\rvert$ converges; hence $(x_k)$ is Cauchy and has a limit $x^\ast\in I$ ($I$ is closed). Continuity of $g$ gives $g(x^\ast) = x^\ast$. If $y^\ast$ were a second fixed point, $\lvert x^\ast - y^\ast\rvert = \lvert g(x^\ast)-g(y^\ast)\rvert\le L\lvert x^\ast-y^\ast\rvert$ forces $x^\ast=y^\ast$. Finally $\lvert x_k - x^\ast\rvert = \lvert g(x_{k-1})-g(x^\ast)\rvert\le L\lvert x_{k-1}-x^\ast\rvert$, and the a posteriori bound follows from $\lvert x_{k-1}-x^\ast\rvert\le\lvert x_{k-1}-x_k\rvert + \lvert x_k - x^\ast\rvert$. $\square$

**Local behaviour.** Near the fixed point, Taylor expansion gives

$$
e_{k+1} = g(x^\ast + e_k) - g(x^\ast) = g'(x^\ast)\,e_k + \tfrac12 g''(\xi_k)e_k^2 .
$$

Hence:

- if $0<\lvert g'(x^\ast)\rvert<1$: **linear** convergence with $C = \lvert g'(x^\ast)\rvert$; if $g'(x^\ast)<0$ the iterates alternate sides of $x^\ast$, if $g'(x^\ast)>0$ they approach monotonically;
- if $\lvert g'(x^\ast)\rvert>1$: the fixed point **repels** and the iteration diverges from it;
- if $g'(x^\ast)=0$ and $g''(x^\ast)\neq0$: **quadratic** convergence. More generally, if $g'(x^\ast)=\cdots=g^{(q-1)}(x^\ast)=0\neq g^{(q)}(x^\ast)$, the order is $q$. Newton's method (Section 22.4) is precisely the choice of $g$ that makes $g'(x^\ast)=0$.

> [!example] Worked Example 22.10 — Solving $x = \cos x$
> **Step 1 — check the theorem.** Take $I = [0,1]$. Then $\cos$ maps $I$ into $[\cos1,1] = [0.540,1]\subset I$, and $\lvert g'(x)\rvert = \lvert\sin x\rvert\le\sin 1 = 0.841<1$ on $I$. So there is a unique fixed point and the iteration converges from any $x_0\in I$.
>
> **Step 2 — iterate from $x_0 = 1$.**
>
> | $k$ | 1 | 2 | 3 | 4 | 5 | 6 | 10 |
> |---|---|---|---|---|---|---|---|
> | $x_k$ | 0.540302 | 0.857553 | 0.654290 | 0.793480 | 0.701369 | 0.763960 | 0.744237 |
>
> **Step 3 — rate.** The fixed point is $x^\ast = 0.7390851332$ and $g'(x^\ast) = -\sin x^\ast = -0.6736$. The iterates alternate about $x^\ast$ (negative $g'$) and each error is about $0.67$ times the previous one: only $-\log_{10}0.674\approx 0.17$ digits per iteration, so about 90 iterations are needed for full double precision. Linear convergence with a constant this close to 1 is slow — the motivation for Newton's method.

The **cobweb diagram** makes the geometry visible: from $(x_k, x_k)$ on the line $y = x$, move vertically to the curve $y = g(x)$, then horizontally back to $y=x$; the abscissa is $x_{k+1}$. A slope $g'$ between $-1$ and $0$ produces an inward spiral.

```desmos-graph
left=0.4; right=1.1; bottom=0.4; top=1.1
---
y=\cos(x)|blue
y=x|black|dashed
x=1|red|0.5403<y<1
y=0.5403|red|0.5403<x<1
x=0.5403|red|0.5403<y<0.8576
y=0.8576|red|0.5403<x<0.8576
x=0.8576|red|0.6543<y<0.8576
y=0.6543|red|0.6543<x<0.8576
x=0.6543|red|0.6543<y<0.7935
y=0.7935|red|0.6543<x<0.7935
(0.7391,0.7391)|label:x*
```

*Cobweb diagram for $x_{k+1}=\cos x_k$ from $x_0=1$: because $-1<g'(x^\ast)<0$ the path spirals inward to the fixed point $x^\ast\approx0.739$.*

> [!example] Worked Example 22.11 — Three rearrangements of $x^2 = 2$
> | $g(x)$ | $g'(x^\ast)$ at $x^\ast=\sqrt2$ | Behaviour from $x_0 = 1$ |
> |---|---|---|
> | $2/x$ | $-1$ | $1, 2, 1, 2,\ldots$ — oscillates forever |
> | $x - (x^2-2)/3$ | $1 - 2\sqrt2/3 = 0.0572$ | $1.3333, 1.4074, 1.41381, 1.414190, 1.4142122$ — linear, gains $\approx1.24$ digits per step |
> | $\tfrac12(x+2/x)$ | $0$ | $1.5, 1.41667, 1.4142157, 1.41421356237469, 1.414213562373095$ — quadratic |
>
> The third rearrangement is Newton's method (Heron's ancient square-root algorithm): its errors are $8.6\times10^{-2}, 2.5\times10^{-3}, 2.1\times10^{-6}, 1.6\times10^{-12}$ and then machine precision — digits double each step. The choice of $g$ is everything.

### The secant method

Newton's method (next section) needs $f'$. If the derivative is unavailable or expensive — $f$ may be the output of a circuit simulation — replace it by the slope of the chord through the last two iterates:

$$
x_{k+1} = x_k - f(x_k)\,\frac{x_k - x_{k-1}}{f(x_k) - f(x_{k-1})}.
$$

Geometrically, $x_{k+1}$ is where the secant line through $(x_{k-1},f(x_{k-1}))$ and $(x_k, f(x_k))$ crosses the axis. Only one *new* function evaluation is needed per step.

> [!theorem] Convergence of the secant method
> If $f\in C^2$ near a simple root $x^\ast$ ($f'(x^\ast)\neq0$) and $x_0,x_1$ are sufficiently close to $x^\ast$, the secant iterates converge with order
> $$
> q = \varphi = \frac{1+\sqrt5}{2}\approx1.618,\qquad\text{and}\quad e_{k+1}\approx\frac{f''(x^\ast)}{2f'(x^\ast)}\,e_k\,e_{k-1}.
> $$

*Derivation sketch.* Expanding $f$ about $x^\ast$ to second order in the secant formula gives $e_{k+1}\approx M e_k e_{k-1}$ with $M = f''(x^\ast)/(2f'(x^\ast))$. Try $\lvert e_{k+1}\rvert\approx A\lvert e_k\rvert^q$; then $\lvert e_{k-1}\rvert\approx(\lvert e_k\rvert/A)^{1/q}$ and substitution gives $A\lvert e_k\rvert^q\approx M A^{-1/q}\lvert e_k\rvert^{1+1/q}$. Matching exponents, $q = 1 + 1/q$, i.e. $q^2 - q - 1 = 0$, whose positive root is the golden ratio. $\square$

Per function evaluation the secant method is actually more efficient than Newton when $f'$ costs as much as $f$: two secant steps (two evaluations) multiply the digits by $\varphi^2 = 2.618 > 2$, whereas one Newton step (two evaluations: $f$ and $f'$) multiplies them by $2$.

> [!example] Worked Example 22.12 — Secant method for $x^3 - x - 2 = 0$
> $f(1) = -2$, $f(2) = 4$, so $[1,2]$ is a bracket. Start with $x_0 = 1$, $x_1 = 2$.
>
> **Step 1.** $x_2 = 2 - 4\cdot\dfrac{2-1}{4-(-2)} = 2 - \dfrac{4}{6} = 1.333333$.
>
> **Subsequent steps** (root $x^\ast = 1.5213797068$):
>
> | $k$ | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
> |---|---|---|---|---|---|---|---|
> | $x_k$ | 1.333333 | 1.462687 | 1.531169 | 1.520926 | 1.521376 | 1.5213797080 | 1.5213797068 |
> | $\lvert e_k\rvert$ | $1.9\times10^{-1}$ | $5.9\times10^{-2}$ | $9.8\times10^{-3}$ | $4.5\times10^{-4}$ | $3.4\times10^{-6}$ | $1.2\times10^{-9}$ | $3\times10^{-15}$ |
>
> **Order check.** $\ln(1.2\times10^{-9}/3.4\times10^{-6})/\ln(3.4\times10^{-6}/4.5\times10^{-4}) = (-7.95)/(-4.89) = 1.63\approx\varphi$. Note that the secant method does not maintain a bracket: $x_2$ and $x_3$ both lie below the root, and $x_4$ overshoots it.

**Regula falsi** (false position) uses the same chord formula but keeps the two points that bracket the root. It always converges, but typically only linearly, because one endpoint tends to get "stuck". The Illinois modification (halve the stored function value at a stuck endpoint) and Brent's method repair this.

### Stopping criteria

An iteration needs a rule for stopping. Common tests, usually combined:

- **step size:** $\lvert x_{k+1}-x_k\rvert\le \text{atol} + \text{rtol}\,\lvert x_{k+1}\rvert$;
- **residual:** $\lvert f(x_k)\rvert\le\tau$;
- **bracket width** (bracketing methods): $b_k - a_k\le\text{tol}$;
- **iteration cap** — always, to guard against non-convergence.

Neither the step test nor the residual test alone is reliable. For a function with a very flat zero ($f'(x^\ast)$ small), the residual can be tiny far from the root; for a steep function, the residual can be large even at the floating-point number closest to the root. The **attainable accuracy** of any root-finder is limited by conditioning: if $f$ can only be evaluated to absolute accuracy $\delta$, the root is determined only to about $\delta/\lvert f'(x^\ast)\rvert$ (the "condition number of a root" is $1/\lvert f'(x^\ast)\rvert$), and much worse for a multiple root.

### Polynomial roots

For polynomials, special methods exist. The most robust general approach — used by `numpy.roots` and MATLAB `roots` — is to form the **companion matrix** whose characteristic polynomial is $p(x)$ and compute its eigenvalues by the QR algorithm (Ch. 10). Once one root $r$ is found, **deflation** divides $p(x)$ by $(x-r)$ (synthetic division) to obtain a lower-degree polynomial for the remaining roots; deflation should proceed from the smallest-magnitude roots to limit error growth. Remember Wilkinson's warning (Section 22.2): roots are often ill-conditioned functions of the coefficients, so if you have the polynomial in factored or state-space form, do not multiply it out.

> [!warning] Common Pitfall — Multiple and clustered roots
> At a root of multiplicity $m>1$, $f$ touches the axis without (for even $m$) changing sign, so bracketing methods may not even see it, and $f$ is so flat that rounding errors of size $\delta$ smear the root over an interval of width about $(\delta\,m!/\lvert f^{(m)}(x^\ast)\rvert)^{1/m}$ — in double precision a double root is typically determined only to about 8 digits. Newton and secant methods slow to linear convergence there (Section 22.4).

> [!tip] Engineering Note — Root-finding inside simulators
> SPICE's DC operating-point analysis is a multidimensional root-finding problem $\mathbf{F}(\mathbf{v}) = \mathbf{0}$ (Kirchhoff's current law at every node). It uses Newton's method with damping, "source stepping" (ramping supplies from zero) and "gmin stepping" (adding small conductances to ground) — all techniques to bring the starting point inside the region of convergence. The message "no convergence in DC analysis" is a root-finder giving up.

## 22.4 Bisection and Newton–Raphson methods

These two methods are the extremes of the root-finding spectrum: bisection is slow but guaranteed, Newton–Raphson is fast but needs a good start. Every serious root-finder is built from one or both.

### The bisection method

Given a bracket $[a_0,b_0]$ with $f(a_0)f(b_0)<0$, repeat for $k = 0,1,2,\ldots$:

1. $m_k = \tfrac12(a_k + b_k)$;
2. if $f(a_k)f(m_k)<0$ set $[a_{k+1},b_{k+1}] = [a_k,m_k]$, otherwise $[a_{k+1},b_{k+1}] = [m_k,b_k]$.

The bracket halves every step and always contains a root.

> [!theorem] Bisection error bound
> If $f$ is continuous on $[a_0,b_0]$ with $f(a_0)f(b_0)<0$, the midpoints satisfy
> $$
> \lvert m_k - x^\ast\rvert\le\frac{b_0 - a_0}{2^{k+1}}
> $$
> for some root $x^\ast\in[a_0,b_0]$. Hence an absolute accuracy $\tau$ is guaranteed after
> $$
> k\ \ge\ \log_2\frac{b_0-a_0}{\tau} - 1\quad\text{iterations}.
> $$

*Proof.* By induction $b_k - a_k = (b_0-a_0)/2^k$ and $[a_k,b_k]$ contains a root by the intermediate value theorem; the midpoint is within half the width of every point of the bracket. $\square$

The error bound halves each step, so bisection behaves like a linearly convergent method with $C = \tfrac12$: one binary digit, or $\log_{10}2\approx0.3$ decimal digits, per function evaluation. Its virtues: it needs only continuity and sign evaluations, the iteration count is known *in advance*, and it is immune to round-off misbehaviour (as long as the sign of $f$ is computed correctly). In floating point compute the midpoint as `a + (b - a)/2` and test signs with `copysign` or `sign(fa) != sign(fm)` rather than `fa*fm < 0`, which can underflow to zero.

```python
def bisect(f, a, b, tol=1e-12, maxit=200):
    fa = f(a)
    for _ in range(maxit):
        m = a + (b - a) / 2
        fm = f(m)
        if fm == 0 or (b - a) / 2 < tol:
            return m
        if (fa < 0) != (fm < 0):   # sign change in [a, m]
            b = m
        else:
            a, fa = m, fm
    return a + (b - a) / 2
```

> [!example] Worked Example 22.13 — Bisection for $x^3 - x - 2 = 0$
> **Bracket.** $f(1) = -2<0$, $f(2) = 4>0$.
>
> | $k$ | $a_k$ | $b_k$ | $m_k$ | $f(m_k)$ |
> |---|---|---|---|---|
> | 0 | 1 | 2 | 1.5 | $-0.125$ |
> | 1 | 1.5 | 2 | 1.75 | $1.609$ |
> | 2 | 1.5 | 1.75 | 1.625 | $0.666$ |
> | 3 | 1.5 | 1.625 | 1.5625 | $0.252$ |
> | 4 | 1.5 | 1.5625 | 1.53125 | $0.0591$ |
> | 5 | 1.5 | 1.53125 | 1.515625 | $-0.0341$ |
> | 6 | 1.515625 | 1.53125 | 1.5234375 | $0.0123$ |
> | 7 | 1.515625 | 1.5234375 | 1.51953125 | $-0.0110$ |
>
> **Error check.** After 8 midpoints, $\lvert m_7 - x^\ast\rvert = \lvert1.51953 - 1.52138\rvert = 1.8\times10^{-3}\le 1/2^{8} = 3.9\times10^{-3}$. ✓
>
> **Iterations for $\tau = 10^{-6}$.** $\log_2(1/10^{-6}) - 1 = 19.93 - 1$, so $k = 19$, i.e. 20 midpoint evaluations; for $\tau = 10^{-12}$ about 40. Compare Newton below.

### Newton–Raphson: derivation

Replace $f$ near the current iterate by its tangent line (first-order Taylor polynomial, Ch. 4):

$$
f(x)\approx f(x_k) + f'(x_k)(x - x_k).
$$

Setting the linear model to zero gives the next iterate

$$
x_{k+1} = x_k - \frac{f(x_k)}{f'(x_k)}.
$$

Equivalently, Newton is fixed-point iteration with $g(x) = x - f(x)/f'(x)$; differentiating, $g'(x) = f(x)f''(x)/f'(x)^2$, which **vanishes at a simple root**. By Section 22.3 this already suggests quadratic convergence.

```desmos-graph
left=1.3; right=2.15; bottom=-1; top=5
---
y=x^{3}-x-2|blue
y=0|black
y=4+11(x-2)|red|1.55<x<2.1
y=0.7453+7.0331(x-1.63636)|orange|1.45<x<1.75
y=0.05394+6.0263(x-1.53039)|purple|1.4<x<1.65
x=2|grey|dashed|0<y<4
x=1.63636|grey|dashed|0<y<0.7453
x=1.53039|grey|dashed|0<y<0.05394
(2,4)|label:x0
(1.63636,0)|label:x1
(1.53039,0)|label:x2
(1.52138,0)|label:root
```

*Newton–Raphson for $f(x)=x^3-x-2$ from $x_0=2$: each tangent line (red, orange, purple) meets the axis at the next iterate, $x_1=1.6364$, $x_2=1.5304$, $x_3=1.52144$; the tangents become ever better models as the iterates close in on the root at $1.52138$.*

### Quadratic convergence

> [!theorem] Local quadratic convergence of Newton's method
> Let $f\in C^2$ on an interval containing a simple root $x^\ast$ ($f(x^\ast)=0$, $f'(x^\ast)\neq0$). Then there is $\delta>0$ such that for every $x_0$ with $\lvert x_0 - x^\ast\rvert<\delta$ the Newton iterates converge to $x^\ast$ and
> $$
> \lim_{k\to\infty}\frac{e_{k+1}}{e_k^2} = \frac{f''(x^\ast)}{2f'(x^\ast)}.
> $$

*Proof.* Taylor's theorem about $x_k$ evaluated at $x^\ast$:

$$
0 = f(x^\ast) = f(x_k) + f'(x_k)(x^\ast - x_k) + \tfrac12 f''(\xi_k)(x^\ast - x_k)^2 .
$$

Divide by $f'(x_k)$ and rearrange using the Newton formula $x_{k+1} = x_k - f(x_k)/f'(x_k)$:

$$
x_{k+1} - x^\ast = \frac{f''(\xi_k)}{2f'(x_k)}\,(x_k - x^\ast)^2,\qquad\text{i.e.}\quad e_{k+1} = \frac{f''(\xi_k)}{2f'(x_k)}\,e_k^2 .
$$

On a small interval around $x^\ast$, $\lvert f'\rvert\ge m>0$ and $\lvert f''\rvert\le M$, so $\lvert e_{k+1}\rvert\le\frac{M}{2m}\lvert e_k\rvert^2$. If $\frac{M}{2m}\lvert e_0\rvert<1$, then $\lvert e_{k+1}\rvert\le\left(\frac{M}{2m}\lvert e_0\rvert\right)\lvert e_k\rvert$ and the errors decrease to zero; as $x_k\to x^\ast$ also $\xi_k\to x^\ast$, giving the limit. $\square$

The constant $M/(2m)$ explains *when* Newton is happy: small curvature $f''$ and steep slope $f'$. The radius $\delta\approx 2m/M$ of guaranteed convergence can be very small for strongly curved functions — the exponential of a diode being the prime example. (The Newton–Kantorovich theorem gives a checkable sufficient condition based on $f(x_0)$, $f'(x_0)$ and a Lipschitz bound on $f'$.)

> [!example] Worked Example 22.14 — Newton for $x^3 - x - 2 = 0$
> $f'(x) = 3x^2 - 1$. Start from the bisection midpoint $x_0 = 1.5$.
>
> **Step 1.** $f(1.5) = 3.375 - 1.5 - 2 = -0.125$, $f'(1.5) = 5.75$, $x_1 = 1.5 + 0.125/5.75 = 1.5217391304$.
>
> **Further steps.**
>
> | $k$ | $x_k$ | $\lvert e_k\rvert$ | $\lvert e_k\rvert/\lvert e_{k-1}\rvert^2$ |
> |---|---|---|---|
> | 0 | 1.5 | $2.14\times10^{-2}$ | — |
> | 1 | 1.521739130434783 | $3.59\times10^{-4}$ | 0.79 |
> | 2 | 1.521379805964786 | $9.9\times10^{-8}$ | 0.77 |
> | 3 | 1.521379706804575 | $7.5\times10^{-15}$ | 0.77 |
> | 4 | 1.521379706804568 | $<10^{-16}$ | — |
>
> **Check against theory.** $f''(x^\ast)/(2f'(x^\ast)) = 6(1.52138)/(2\cdot5.9438) = 0.768$, matching the last column. The number of correct digits goes $2\to3\to7\to14\to16$: four Newton steps versus about fifty bisection steps for the same accuracy.

### Newton's method in circuit analysis

> [!example] Worked Example 22.15 — Diode operating point
> A $5\ \mathrm{V}$ source drives a $1\ \mathrm{k\Omega}$ resistor in series with a silicon diode, $I_D = I_S(e^{V/V_T}-1)$ with $I_S = 10^{-14}\ \mathrm{A}$ and $V_T = 25.852\ \mathrm{mV}$ (300 K). Find the diode voltage $V$.
>
> **Step 1 — KCL.** Current through the resistor equals diode current:
> $$
> F(V) = \frac{5 - V}{1000} - I_S\left(e^{V/V_T}-1\right) = 0,\qquad F'(V) = -\frac{1}{1000} - \frac{I_S}{V_T}e^{V/V_T}.
> $$
>
> **Step 2 — good start.** Engineering judgement says $V\approx0.7\ \mathrm{V}$.
>
> | $k$ | $V_k$ (V) | $F(V_k)$ (A) |
> |---|---|---|
> | 1 | 0.693518214 | $-1.7\times10^{-4}$ |
> | 2 | 0.692561670 | $-3.0\times10^{-6}$ |
> | 3 | 0.692543639 | $-1.0\times10^{-9}$ |
> | 4 | 0.692543633 | $\approx10^{-16}$ |
>
> So $V = 0.69254\ \mathrm{V}$ and $I_D = (5-0.69254)/1000 = 4.307\ \mathrm{mA}$. Quadratic convergence is evident in the residuals.
>
> **Step 3 — bad start.** From $V_0 = 0$ the diode term is negligible, the tangent is essentially the resistor's line, and $V_1 = 4.99999999807\ \mathrm{V}$. At $5\ \mathrm{V}$ the exponential is astronomically steep ($e^{5/0.0259}\approx 10^{84}$), so each subsequent Newton step reduces $V$ by only about $V_T$: $4.9741, 4.9483, 4.9224,\ldots$ It would take about 165 iterations to crawl down to $0.7\ \mathrm{V}$ — and with a slightly larger overshoot $e^{V/V_T}$ would overflow to $\infty$.
>
> **Cure.** SPICE-type simulators **limit** the change in junction voltage per iteration (e.g. to a few $V_T$ on a logarithmic scale, the "pnjlim" algorithm) or **damp** the step, $V_{k+1} = V_k - \lambda F/F'$ with $\lambda\le1$ chosen so that $\lvert F\rvert$ decreases.

### When Newton fails

Newton's method is only *locally* convergent. The classic failure modes are:

1. **Zero or tiny derivative.** If $f'(x_k)\approx0$ the step is huge (the tangent is nearly horizontal).
2. **Cycling.** For $f(x) = x^3 - 2x + 2$ and $x_0 = 0$: $x_1 = 0 - 2/(-2) = 1$, $x_2 = 1 - 1/1 = 0$, and the iterates alternate $0,1,0,1,\ldots$ forever, although a real root exists at $x^\ast\approx -1.7693$.
3. **Divergence.** For $f(x) = \arctan x$ (root $0$), Newton converges only for $\lvert x_0\rvert<1.3917$; from $x_0 = 1.5$ the iterates oscillate with growing magnitude ($-1.69, 2.32, -5.11,\ldots$).
4. **Convergence to an unwanted root** — the basins of attraction of different roots can be intricately interleaved (in the complex plane they form fractals).
5. **Overflow** of $f$ itself in the first steps (Example 22.15).

**Safeguarding.** Robust implementations combine Newton with a bracket: if the Newton step would leave the current bracket or fails to reduce $\lvert f\rvert$ sufficiently, take a bisection step instead. This "Newton–bisection hybrid" (and Brent's method, which does the same with inverse quadratic interpolation) inherits bisection's guarantee and Newton's speed.

### Multiple roots

At a root of multiplicity $m$, $f(x) = (x-x^\ast)^m h(x)$ with $h(x^\ast)\neq0$. Then $g(x) = x - f/f'$ satisfies $g'(x^\ast) = 1 - 1/m$, so Newton converges only **linearly** with constant $(m-1)/m$ — for a double root, the error just halves each step. Two remedies:

- **Modified Newton**, if $m$ is known: $x_{k+1} = x_k - m\,f(x_k)/f'(x_k)$ restores quadratic convergence (now $g'(x^\ast) = 0$).
- Apply Newton to $\mu(x) = f(x)/f'(x)$, which has a *simple* root at $x^\ast$ whatever $m$ is (requires $f''$).

> [!example] Worked Example 22.16 — A double root
> $f(x) = (x-1)^2(x+2)$ has a double root at $1$. From $x_0 = 2$:
>
> | $k$ | 1 | 2 | 3 | 4 | 5 |
> |---|---|---|---|---|---|
> | Newton $e_k$ | 0.556 | 0.298 | 0.155 | 0.0796 | 0.0403 |
> | Modified ($m=2$) $e_k$ | 0.111 | $1.95\times10^{-3}$ | $6.3\times10^{-7}$ | $6.7\times10^{-14}$ | — |
>
> Plain Newton's errors approach the ratio $\tfrac12$ predicted by $(m-1)/m$; the modified iteration is quadratic again. A practical symptom of a multiple root: Newton suddenly converges linearly with ratio $\approx 1/2$ or $2/3$.

### Newton's method for systems

For $\mathbf F:\mathbb R^n\to\mathbb R^n$, the tangent line becomes the tangent *plane* given by the **Jacobian** $\mathbf J(\mathbf x) = [\partial F_i/\partial x_j]$ (Ch. 6):

$$
\mathbf F(\mathbf x)\approx\mathbf F(\mathbf x_k) + \mathbf J(\mathbf x_k)(\mathbf x - \mathbf x_k)\quad\Longrightarrow\quad
\mathbf J(\mathbf x_k)\,\Delta\mathbf x_k = -\mathbf F(\mathbf x_k),\qquad \mathbf x_{k+1} = \mathbf x_k + \Delta\mathbf x_k .
$$

Each step solves a *linear* system (never form $\mathbf J^{-1}$). Convergence is again quadratic near a root at which $\mathbf J$ is nonsingular. This is exactly the inner loop of nonlinear circuit simulation, power-flow analysis in power systems (the Newton–Raphson load flow), and implicit ODE solvers (Section 22.9). When $\mathbf J$ is expensive, **quasi-Newton** methods (Broyden) update an approximation to it — a theme of Ch. 23.

```python
import numpy as np
def newton_system(F, J, x, tol=1e-12, maxit=50):
    for _ in range(maxit):
        dx = np.linalg.solve(J(x), -F(x))
        x = x + dx
        if np.linalg.norm(dx) <= tol * (1 + np.linalg.norm(x)):
            break
    return x
```

> [!example] Worked Example 22.17 — A 2×2 nonlinear system
> Solve $x^2 + y^2 = 4$, $xy = 1$ (intersection of a circle and a hyperbola) starting from $(x_0,y_0) = (2, 0.5)$.
>
> **Set-up.** $\mathbf F = \begin{bmatrix}x^2+y^2-4\\ xy-1\end{bmatrix}$, $\mathbf J = \begin{bmatrix}2x&2y\\ y&x\end{bmatrix}$.
>
> **Step 1.** $\mathbf F(2,0.5) = (0.25, 0)^\mathsf T$, $\mathbf J = \begin{bmatrix}4&1\\0.5&2\end{bmatrix}$. Solving $\mathbf J\Delta = -\mathbf F$: $\det\mathbf J = 7.5$, $\Delta = \frac{1}{7.5}\begin{bmatrix}2&-1\\-0.5&4\end{bmatrix}\begin{bmatrix}-0.25\\0\end{bmatrix} = (-0.06667, 0.01667)^\mathsf T$, so $(x_1,y_1) = (1.93333, 0.51667)$.
>
> **Further steps.** $\lVert\mathbf F\rVert$: $0.25\to4.9\times10^{-3}\to3.4\times10^{-6}\to2.5\times10^{-12}\to 4\times10^{-16}$. The solution is $(1.9318517, 0.5176381)$.
>
> **Check.** Exactly, $x^2+y^2=4$ and $xy=1$ give $(x+y)^2 = 6$, $(x-y)^2 = 2$, so $x = (\sqrt6+\sqrt2)/2 = 1.9318517$, $y = (\sqrt6-\sqrt2)/2 = 0.5176381$. ✓

> [!warning] Common Pitfall — Trusting a converged Newton iteration blindly
> Newton converging tells you that you found *a* root — not necessarily the physical one. A transistor circuit with positive feedback (a latch, a Schmitt trigger) has several DC operating points; which one Newton finds depends on the starting guess. Always sanity-check the answer (signs of currents, device regions) and, for multi-stable circuits, perform a transient or continuation analysis.

> [!tip] Engineering Note — Newton in hardware
> Reciprocal and square-root units in many processors and DSPs start from a small lookup-table estimate and refine it with one or two Newton steps, using only multiplications: for $1/a$, apply Newton to $f(x) = 1/x - a$ to get $x_{k+1} = x_k(2 - a x_k)$; for $1/\sqrt a$, $x_{k+1} = \tfrac12 x_k(3 - a x_k^2)$ (the iteration behind the famous "fast inverse square root"). An 8-bit table estimate becomes 16 bits after one step and 32 bits after two — quadratic convergence in silicon.

## 22.5 Interpolation

### Motivation

Interpolation constructs a function that passes *exactly* through given data points. Engineers use it to read between the entries of a datasheet table (a transistor's $h_{FE}$ versus $I_C$, a thermistor's resistance versus temperature), to implement function generators by table lookup, to resample a signal at a new rate or a fractional delay, and to draw smooth curves through simulation output. Equally important, interpolation is the *theoretical engine* behind most of this chapter: numerical differentiation (Section 22.7) differentiates an interpolant, and quadrature (Section 22.8) integrates one. The error of those methods is inherited from the interpolation error formula derived here.

### Existence and uniqueness

> [!theorem] Polynomial interpolation theorem
> Given $n+1$ distinct nodes $x_0,\ldots,x_n$ and values $y_0,\ldots,y_n$, there is exactly one polynomial $p_n$ of degree at most $n$ with $p_n(x_i) = y_i$ for $i = 0,\ldots,n$.

*Proof.* Writing $p_n(x) = \sum_{k=0}^{n} a_kx^k$, the interpolation conditions form the linear system $\mathbf V\mathbf a = \mathbf y$ with the **Vandermonde matrix** $V_{ik} = x_i^k$, whose determinant is $\prod_{i<j}(x_j - x_i)\neq0$ for distinct nodes. (Alternatively: if $p$ and $q$ both interpolate, $p-q$ has degree $\le n$ and $n+1$ roots, hence is zero; existence follows from the Lagrange construction below.) $\square$

Solving the Vandermonde system is a poor way to *compute* the interpolant: $\mathbf V$ is notoriously ill-conditioned (for 20 equally spaced nodes on $[0,1]$, $\kappa(\mathbf V)\approx10^{16}$). The Lagrange and Newton forms below represent the same unique polynomial in better-conditioned bases.

### Lagrange form

> [!definition] Lagrange basis polynomials
> For distinct nodes $x_0,\ldots,x_n$,
> $$
> \ell_i(x) = \prod_{\substack{j=0\\ j\neq i}}^{n}\frac{x - x_j}{x_i - x_j},\qquad\text{so that}\quad\ell_i(x_k) = \delta_{ik},
> $$
> and the interpolating polynomial is
> $$
> p_n(x) = \sum_{i=0}^{n} y_i\,\ell_i(x).
> $$

Each $\ell_i$ is a "cardinal" function: 1 at its own node, 0 at all others. The Lagrange form is ideal for theory (quadrature weights are $\int\ell_i$, differentiation weights are $\ell_i'$). For computation, the algebraically equivalent **barycentric formula**

$$
p_n(x) = \frac{\displaystyle\sum_{i=0}^{n}\frac{w_i}{x - x_i}\,y_i}{\displaystyle\sum_{i=0}^{n}\frac{w_i}{x-x_i}},\qquad w_i = \frac{1}{\prod_{j\neq i}(x_i - x_j)},
$$

costs $O(n)$ per evaluation once the weights are known and is numerically stable; it is what modern software uses.

### Newton's divided-difference form

> [!definition] Divided differences
> $$
> f[x_i] = f(x_i),\qquad f[x_i,\ldots,x_{i+k}] = \frac{f[x_{i+1},\ldots,x_{i+k}] - f[x_i,\ldots,x_{i+k-1}]}{x_{i+k}-x_i}.
> $$

> [!theorem] Newton form of the interpolant
> $$
> p_n(x) = f[x_0] + f[x_0,x_1](x-x_0) + f[x_0,x_1,x_2](x-x_0)(x-x_1) + \cdots + f[x_0,\ldots,x_n]\prod_{j=0}^{n-1}(x - x_j).
> $$

*Proof sketch.* Let $p_{k}$ interpolate at $x_0,\ldots,x_k$. Then $p_{k+1} - p_k$ has degree $\le k+1$ and vanishes at $x_0,\ldots,x_k$, so $p_{k+1} = p_k + c_{k+1}\prod_{j=0}^{k}(x - x_j)$, and $c_{k+1}$ is the leading coefficient of $p_{k+1}$. Comparing the leading coefficients of two interpolants on overlapping node sets yields the recursive definition of the divided difference. $\square$

Two practical advantages: adding a new data point adds one term without recomputing the others, and the form is evaluated cheaply by nested multiplication (Horner's scheme). Divided differences also approximate derivatives: if $f\in C^k$, then $f[x_0,\ldots,x_k] = f^{(k)}(\xi)/k!$ for some $\xi$ in the span of the nodes.

> [!example] Worked Example 22.18 — Cubic through four points, two ways
> Data: $(0,1), (1,3), (2,2), (3,5)$. Find $p_3$ and evaluate it at $x = 1.5$.
>
> **Divided-difference table.**
>
> | $x_i$ | $f[x_i]$ | first | second | third |
> |---|---|---|---|---|
> | 0 | 1 | | | |
> | | | $(3-1)/1 = 2$ | | |
> | 1 | 3 | | $(-1-2)/2 = -1.5$ | |
> | | | $(2-3)/1 = -1$ | | $(2-(-1.5))/3 = 7/6$ |
> | 2 | 2 | | $(3-(-1))/2 = 2$ | |
> | | | $(5-2)/1 = 3$ | | |
> | 3 | 5 | | | |
>
> **Newton form.** The top diagonal gives
> $$
> p_3(x) = 1 + 2x - 1.5\,x(x-1) + \tfrac76\,x(x-1)(x-2).
> $$
>
> **Evaluate at 1.5.** $1 + 3 - 1.5(0.75) + \tfrac76(1.5)(0.5)(-0.5) = 4 - 1.125 - 0.4375 = 2.4375$.
>
> **Lagrange check.** With $\ell_0(1.5) = \frac{(0.5)(-0.5)(-1.5)}{(-1)(-2)(-3)} = -0.0625$, $\ell_1(1.5) = \frac{(1.5)(-0.5)(-1.5)}{(1)(-1)(-2)} = 0.5625$, $\ell_2(1.5) = \frac{(1.5)(0.5)(-1.5)}{(2)(1)(-1)} = 0.5625$, $\ell_3(1.5) = \frac{(1.5)(0.5)(-0.5)}{(3)(2)(1)} = -0.0625$ (they sum to 1, as they must):
> $p_3(1.5) = 1(-0.0625) + 3(0.5625) + 2(0.5625) + 5(-0.0625) = 2.4375$. ✓
>
> Expanded, $p_3(x) = \tfrac76x^3 - 5x^2 + \tfrac{35}{6}x + 1$, which matches `numpy.polyfit` with degree 3.

### The interpolation error

> [!theorem] Interpolation error formula
> Let $f\in C^{n+1}[a,b]$ and let $p_n$ interpolate $f$ at distinct nodes $x_0,\ldots,x_n\in[a,b]$. For every $x\in[a,b]$ there is $\xi = \xi(x)\in(a,b)$ with
> $$
> f(x) - p_n(x) = \frac{f^{(n+1)}(\xi)}{(n+1)!}\,\omega_{n+1}(x),\qquad \omega_{n+1}(x) = \prod_{i=0}^{n}(x - x_i).
> $$

*Proof.* Fix $x$ not a node and define $\phi(t) = f(t) - p_n(t) - \dfrac{f(x)-p_n(x)}{\omega_{n+1}(x)}\,\omega_{n+1}(t)$. Then $\phi$ vanishes at the $n+2$ points $x_0,\ldots,x_n,x$. By Rolle's theorem applied repeatedly, $\phi^{(n+1)}$ has a zero $\xi$. Since $p_n^{(n+1)} = 0$ and $\omega_{n+1}^{(n+1)} = (n+1)!$, $0 = f^{(n+1)}(\xi) - (n+1)!\,(f(x)-p_n(x))/\omega_{n+1}(x)$. $\square$

The formula separates two influences: the **smoothness of $f$** (through $f^{(n+1)}$) and the **placement of the nodes** (through $\omega_{n+1}$). The second factor is ours to choose.

> [!example] Worked Example 22.19 — Quadratic interpolation of $\sin x$
> Interpolate $\sin x$ at $0,\ \pi/4,\ \pi/2$ and estimate the error at $x = \pi/8$.
>
> **Divided differences.** $f[x_0,x_1] = 0.70711/0.78540 = 0.90032$, $f[x_1,x_2] = (1-0.70711)/0.78540 = 0.37292$, $f[x_0,x_1,x_2] = (0.37292-0.90032)/1.5708 = -0.33575$.
>
> **Interpolant.** $p_2(x) = 0.90032x - 0.33575\,x(x - \pi/4)$; at $\pi/8 = 0.392699$: $p_2 = 0.35355 + 0.05178 = 0.40533$.
>
> **Actual error.** $\sin(\pi/8) = 0.38268$, so $\lvert f - p_2\rvert = 0.02265$.
>
> **Bound.** $\lvert f'''\rvert = \lvert\cos\rvert\le1$, so $\lvert f-p_2\rvert\le\frac{1}{3!}\lvert(\pi/8)(\pi/8-\pi/4)(\pi/8-\pi/2)\rvert = \frac{1}{6}(0.3927)(0.3927)(1.1781) = 0.0303$. The bound holds and is not wildly pessimistic.

### The Runge phenomenon

It is tempting to think that more nodes always give a better fit. For **equally spaced** nodes this is false, even for perfectly smooth functions. Runge's classic example is

$$
f(x) = \frac{1}{1+25x^2},\qquad x\in[-1,1].
$$

| $n$ (degree) | 4 | 8 | 12 | 16 | 20 |
|---|---|---|---|---|---|
| $\max\lvert f-p_n\rvert$, equispaced nodes | 0.438 | 1.045 | 3.66 | 14.4 | 59.8 |
| $\max\lvert f-p_n\rvert$, Chebyshev nodes | 0.402 | 0.171 | 0.069 | 0.033 | 0.015 |

With equispaced nodes the error near the ends of the interval *grows exponentially* with $n$; with Chebyshev nodes (defined below) it decreases steadily.

```desmos-graph
left=-1.05; right=1.05; bottom=-0.5; top=2
---
y=\frac{1}{1+25x^{2}}|black
y=-220.941742x^{10}+494.909502x^{8}-381.433824x^{6}+123.359729x^{4}-16.8552036x^{2}+1|red|-1<x<1
y=-46.6329171x^{10}+130.105839x^{8}-133.444756x^{6}+61.4430186x^{4}-12.4765115x^{2}+1|blue|dashed|-1<x<1
(-1,0.03846)|red
(-0.8,0.05882)|red
(-0.6,0.1)|red
(-0.4,0.2)|red
(-0.2,0.5)|red
(0,1)|red
(0.2,0.5)|red
(0.4,0.2)|red
(0.6,0.1)|red
(0.8,0.05882)|red
(1,0.03846)|red
```

*Runge phenomenon: $f(x)=1/(1+25x^2)$ (black), its degree-10 interpolant at 11 equally spaced nodes (red, nodes marked) overshooting to $\approx1.96$ near $x=\pm0.94$ (an error of $1.92$), and the degree-10 interpolant at Chebyshev nodes (blue dashed) staying within $0.11$ of $f$ everywhere.*

**Why it happens.** Two factors in the error formula conspire. First, $\lvert\omega_{n+1}(x)\rvert$ for equispaced nodes is much larger near the ends of the interval than in the middle (by a factor that grows like $2^n$). Second, $f$ has complex poles at $x = \pm j/5$, close to the real interval, so its derivatives grow like $(n+1)!\,5^{n+1}$ — fast enough to beat the $(n+1)!$ in the denominator. Complex analysis (Ch. 18) makes this precise: polynomial interpolants at equispaced nodes converge only inside a certain lens-shaped region of the complex plane, and if a singularity of $f$ lies inside that region, interpolation diverges near the endpoints.

### Chebyshev nodes

Since $f^{(n+1)}$ is beyond our control, minimise the other factor: choose nodes that make $\max_{[-1,1]}\lvert\omega_{n+1}\rvert$ as small as possible.

> [!definition] Chebyshev polynomials and nodes
> The Chebyshev polynomials are $T_k(x) = \cos(k\arccos x)$ on $[-1,1]$, equivalently $T_0 = 1$, $T_1 = x$, $T_{k+1} = 2xT_k - T_{k-1}$. The **Chebyshev nodes** (of the first kind) for degree $n$ are the $n+1$ zeros of $T_{n+1}$:
> $$
> x_i = \cos\!\left(\frac{(2i+1)\pi}{2(n+1)}\right),\qquad i = 0,\ldots,n.
> $$
> On a general interval $[a,b]$ use $\tfrac{a+b}{2}+\tfrac{b-a}{2}x_i$.

The nodes are projections onto the axis of equally spaced points on a semicircle, so they *cluster near the endpoints* — exactly where equispaced interpolation goes wrong.

> [!theorem] Minimax property of Chebyshev nodes
> Among all monic polynomials of degree $n+1$, $2^{-n}T_{n+1}(x)$ has the smallest maximum magnitude on $[-1,1]$, namely $2^{-n}$. Consequently, for Chebyshev nodes $\omega_{n+1} = 2^{-n}T_{n+1}$ and
> $$
> \max_{[-1,1]}\lvert f - p_n\rvert\le\frac{\max\lvert f^{(n+1)}\rvert}{2^n(n+1)!}.
> $$

*Proof sketch.* $2^{-n}T_{n+1}$ is monic and attains $\pm2^{-n}$ alternately at the $n+2$ points $\cos(k\pi/(n+1))$. If a monic $q$ of the same degree had smaller maximum, $2^{-n}T_{n+1}-q$ (degree $\le n$) would alternate in sign at those $n+2$ points and so have $n+1$ zeros — hence be identically zero, a contradiction. $\square$

A finer measure of node quality is the **Lebesgue constant** $\Lambda_n = \max_x\sum_i\lvert\ell_i(x)\rvert$, which bounds how much interpolation can amplify errors in the data and relates the interpolation error to the best possible: $\lVert f - p_n\rVert_\infty\le(1+\Lambda_n)\,\lVert f - p_n^{\ast}\rVert_\infty$, where $p_n^\ast$ is the best (minimax) approximation. For equispaced nodes $\Lambda_n\sim 2^{n+1}/(e\,n\ln n)$ (exponential growth); for Chebyshev nodes $\Lambda_n\le\tfrac{2}{\pi}\ln(n+1)+1$ (logarithmic growth, essentially optimal). With Chebyshev nodes, interpolation of any Lipschitz-continuous function converges, and for analytic functions it converges geometrically. This is the basis of the "Chebfun" philosophy of computing with functions.

### Piecewise interpolation and splines

The other cure for Runge's phenomenon is to **keep the degree low and subdivide**: interpolate by a low-degree polynomial on each subinterval between nodes.

**Piecewise linear interpolation** (what a plotting program does when it "connects the dots") on nodes with spacing $h$ has error

$$
\max\lvert f - s\rvert\le\frac{h^2}{8}\max\lvert f''\rvert,
$$

from the error formula with $n = 1$ and $\max_{[x_i,x_{i+1}]}\lvert(x-x_i)(x-x_{i+1})\rvert = h^2/4$. It converges for any $C^2$ function, but the interpolant has corners.

> [!definition] Cubic spline
> A **cubic spline** on nodes $a = x_0<x_1<\cdots<x_n = b$ is a function $S$ that is a cubic polynomial on each $[x_i,x_{i+1}]$ and has continuous first and second derivatives on $[a,b]$. An **interpolating** cubic spline also satisfies $S(x_i) = y_i$. Two extra conditions are required, typically:
> - **natural**: $S''(a) = S''(b) = 0$;
> - **clamped**: $S'(a)$, $S'(b)$ prescribed;
> - **not-a-knot**: $S'''$ continuous at $x_1$ and $x_{n-1}$ (MATLAB/SciPy default).

**Counting.** $n$ cubics have $4n$ coefficients. Interpolation at both ends of each piece gives $2n$ conditions, continuity of $S'$ and $S''$ at the $n-1$ interior nodes gives $2n-2$; two end conditions complete the count.

**Derivation (moment form).** Let $M_i = S''(x_i)$ and $h_i = x_{i+1}-x_i$. Since $S''$ is linear on each piece, integrating twice and imposing $S(x_i) = y_i$, $S(x_{i+1}) = y_{i+1}$ gives on $[x_i,x_{i+1}]$:

$$
S(x) = M_i\frac{(x_{i+1}-x)^3}{6h_i} + M_{i+1}\frac{(x-x_i)^3}{6h_i} + \left(\frac{y_i}{h_i}-\frac{M_ih_i}{6}\right)(x_{i+1}-x) + \left(\frac{y_{i+1}}{h_i}-\frac{M_{i+1}h_i}{6}\right)(x-x_i).
$$

Requiring $S'$ to be continuous at each interior node $x_i$ yields the **tridiagonal** system

$$
h_{i-1}M_{i-1} + 2(h_{i-1}+h_i)M_i + h_iM_{i+1} = 6\left(\frac{y_{i+1}-y_i}{h_i} - \frac{y_i-y_{i-1}}{h_{i-1}}\right),\qquad i = 1,\ldots,n-1,
$$

which is strictly diagonally dominant, hence uniquely solvable in $O(n)$ operations by the Thomas algorithm (Gaussian elimination without pivoting on a tridiagonal matrix).

> [!theorem] Properties of cubic splines
> 1. **Accuracy.** For $f\in C^4$ and the clamped (or not-a-knot) spline with maximum spacing $h$: $\lVert f - S\rVert_\infty\le\tfrac{5}{384}h^4\lVert f^{(4)}\rVert_\infty$, with $O(h^3)$ for $S'$ and $O(h^2)$ for $S''$.
> 2. **Minimum curvature.** Among all $C^2$ functions interpolating the data (with the same end conditions, natural case), the natural cubic spline minimises $\int_a^b (g'')^2\,dx$ — the bending energy of a thin elastic beam, which is where the name "spline" (a draughtsman's flexible strip) comes from.
> 3. **Locality.** Changing one data value changes $S$ everywhere, but the influence decays geometrically (by a factor of about $2-\sqrt3\approx0.27$ per interval for uniform nodes).

> [!example] Worked Example 22.20 — A natural cubic spline
> Data: $(0,0), (1,1), (2,0), (3,1)$, so $h = 1$ and $M_0 = M_3 = 0$.
>
> **Step 1 — system for $M_1, M_2$.** With $h_i = 1$ the equations read $M_{i-1} + 4M_i + M_{i+1} = 6(y_{i+1} - 2y_i + y_{i-1})$:
> $$
> \begin{aligned}
> 4M_1 + M_2 &= 6(0 - 2 + 0) = -12,\\
> M_1 + 4M_2 &= 6(1 - 0 + 1) = 12,
> \end{aligned}
> $$
> giving $M_1 = -4$, $M_2 = 4$.
>
> **Step 2 — pieces.** On $[0,1]$: $S_0(x) = -4\frac{x^3}{6} + \left(1+\frac46\right)x = -\tfrac23x^3 + \tfrac53x$. On $[1,2]$: $S_1(x) = -\tfrac23(2-x)^3 + \tfrac23(x-1)^3 + \tfrac53(2-x) - \tfrac23(x-1)$. On $[2,3]$, by the symmetry $y_i\mapsto 1-y_{3-i}$, $S_2(x) = 1 - S_0(3-x)$.
>
> **Step 3 — check smoothness at $x=1$.** $S_0'(1) = -2 + \tfrac53 = -\tfrac13$; $S_1'(1) = 2(2-1)^2 + 0 - \tfrac53 - \tfrac23 = -\tfrac13$. ✓ Second derivatives: $S_0''(1) = -4 = M_1 = S_1''(1)$. ✓ And $S_1(1.5) = 0.5$, as symmetry requires.

### Hermite interpolation

If derivatives are known as well as values (for example position *and* velocity, or a circuit waveform and its slope), **Hermite interpolation** matches both. Matching $f$ and $f'$ at $n+1$ nodes determines a unique polynomial of degree $\le 2n+1$, with error $f^{(2n+2)}(\xi)\,\omega_{n+1}(x)^2/(2n+2)!$. The piecewise cubic Hermite interpolant (matching $f, f'$ at both ends of every interval) underlies "PCHIP", a shape-preserving interpolant that avoids the overshoots a cubic spline can produce for monotone data.

> [!warning] Common Pitfall — Extrapolation and high-degree fits
> The interpolation error formula contains $\omega_{n+1}(x)$, which grows like $\lvert x\rvert^{n+1}$ outside the node interval. **Extrapolating** an interpolating polynomial beyond the data is therefore dangerous, and the higher the degree the worse. Likewise, never fit a single polynomial of degree 15 through 16 equally spaced measurements; use a spline, or least squares with a low degree (Section 22.6).

> [!tip] Engineering Note — Interpolation in signal processing
> A DAC with a zero-order hold performs piecewise-constant interpolation; a first-order hold is piecewise linear. The ideal band-limited reconstruction of the sampling theorem (Ch. 16) is interpolation with the sinc kernel, which is the limit of Lagrange interpolation on infinitely many equally spaced nodes. Practical **fractional-delay** filters and sample-rate converters use Lagrange interpolation over 4–8 neighbouring samples, often in the efficient **Farrow** structure, whose polynomial coefficients are exactly the Newton/Lagrange coefficients of this section.

## 22.6 Polynomial approximation

### Interpolation versus approximation

Interpolation forces the approximant through every data point. That is the wrong goal when the data are **noisy** (fitting a polynomial through every noisy sample fits the noise), when there are **more data than parameters** (a thousand measurements of a sensor characteristic that is nearly quadratic), or when we want the **best possible accuracy** for a given degree (a function-evaluation routine in a math library or an FPGA). Approximation instead minimises a measure of the error in a chosen norm:

- the **uniform (minimax, $L^\infty$) norm** $\lVert f-p\rVert_\infty = \max_{x\in[a,b]}\lvert f(x)-p(x)\rvert$ — the worst-case error, used for library functions and lookup tables;
- the **least-squares ($L^2$) norm** $\lVert f-p\rVert_2 = \left(\int_a^b w(x)\,(f-p)^2dx\right)^{1/2}$ or, for data, $\left(\sum_i w_i(y_i - p(x_i))^2\right)^{1/2}$ — the energy of the error, used for data fitting because it is optimal for Gaussian noise (Ch. 20) and leads to *linear* equations.

That polynomials are a sensible choice at all is guaranteed by a classical theorem.

> [!theorem] Weierstrass approximation theorem
> If $f$ is continuous on $[a,b]$, then for every $\varepsilon>0$ there is a polynomial $p$ with $\lVert f - p\rVert_\infty<\varepsilon$.

A constructive proof uses the **Bernstein polynomials** $B_n f(x) = \sum_{k=0}^{n}f(k/n)\binom nk x^k(1-x)^{n-k}$ on $[0,1]$, which converge uniformly to $f$ — but only slowly ($O(1/n)$ for smooth $f$). The theorem guarantees existence; the methods below achieve good accuracy at low degree.

### Taylor polynomials: good locally, poor globally

A Taylor polynomial (Ch. 4) is exact at one point with all derivatives up to order $n$ and its error, $f^{(n+1)}(\xi)(x-x_0)^{n+1}/(n+1)!$, grows rapidly away from $x_0$. For $e^x$ on $[0,1]$, the degree-1 Taylor polynomial $1+x$ has maximum error $e - 2 = 0.718$ at $x = 1$; we shall see that other degree-1 polynomials do seven times better. Taylor polynomials are the right tool for analysis near a point (linearisation, small-signal models) and the wrong tool for approximation over an interval.

### Best uniform (minimax) approximation

> [!theorem] Chebyshev equioscillation theorem
> Let $f\in C[a,b]$. A polynomial $p^\ast$ of degree $\le n$ is the best uniform approximation to $f$ if and only if there are $n+2$ points $a\le t_0<t_1<\cdots<t_{n+1}\le b$ at which the error $e = f-p^\ast$ attains its maximum magnitude with alternating signs:
> $$
> e(t_i) = \pm(-1)^i\lVert f-p^\ast\rVert_\infty .
> $$
> The best approximation exists and is unique.

The theorem is the same alternation argument that proved the minimax property of $T_{n+1}$ in Section 22.5 (which is the special case "best approximation of $x^{n+1}$ by lower degree"). In practice the minimax polynomial is computed by the **Remez exchange algorithm**, which iteratively adjusts a set of $n+2$ reference points until the error equioscillates. The same equioscillation principle defines **equiripple (Parks–McClellan) FIR filters** and **elliptic (Cauer) analog filters** — optimal filter design is minimax approximation.

> [!example] Worked Example 22.21 — Best straight line for $e^x$ on $[0,1]$
> Seek $p^\ast(x) = mx + c$ with $n+2 = 3$ alternation points. Because $e^x$ is convex, the error $e^x - mx - c$ has its extreme values at the endpoints $0, 1$ and at one interior point $\xi$ where $e' = 0$, i.e. $e^{\xi} = m$.
>
> **Equioscillation conditions.** $e(0) = +E$, $e(\xi) = -E$, $e(1) = +E$:
> $$
> 1 - c = E,\qquad e^{\xi} - m\xi - c = -E,\qquad e - m - c = E .
> $$
>
> **Solve.** Subtract the first from the third: $m = e - 1 = 1.718282$. Then $\xi = \ln m = 0.541325$. Adding the first two equations: $1 + m - m\xi - 2c = 0$, so $c = \tfrac12(1 + m - m\xi) = 0.894067$, and $E = 1 - c = 0.105933$.
>
> **Result.** $p^\ast(x) = 1.718282x + 0.894067$ with $\lVert e^x - p^\ast\rVert_\infty = 0.1059$ — versus $0.718$ for the Taylor line $1+x$. A numerical scan of $\lvert e^x - p^\ast(x)\rvert$ on $[0,1]$ confirms the maximum $0.10593$ at $x = 0,\ 0.5413,\ 1$.

### Discrete least squares

Given data $(x_i,y_i)$, $i = 1,\ldots,N$, and a model $p(x) = \sum_{k=0}^{n}a_k\phi_k(x)$ that is *linear in the parameters* (for polynomials $\phi_k = x^k$, but sinusoids, exponentials with known rates, or any fixed basis work identically), least squares chooses $\mathbf a$ to minimise

$$
S(\mathbf a) = \sum_{i=1}^{N}\big(y_i - p(x_i)\big)^2 = \lVert\mathbf y - \mathbf A\mathbf a\rVert_2^2,\qquad A_{ik} = \phi_k(x_i).
$$

Setting $\nabla S = \mathbf 0$ (Ch. 6, Ch. 23) gives the **normal equations**

$$
\mathbf A^\mathsf T\mathbf A\,\mathbf a = \mathbf A^\mathsf T\mathbf y,
$$

which have a unique solution when the columns of $\mathbf A$ are linearly independent (for polynomials: at least $n+1$ distinct $x_i$). Geometrically, $\mathbf A\mathbf a$ is the orthogonal projection of $\mathbf y$ onto the column space of $\mathbf A$: the residual $\mathbf r = \mathbf y - \mathbf A\mathbf a$ is orthogonal to every column, $\mathbf A^\mathsf T\mathbf r = \mathbf 0$ (Ch. 9).

> [!example] Worked Example 22.22 — Fitting a resistor's I–V line
> Measured voltage–current pairs ($V$ in volts, $I$ in mA): $(0,0.1), (1,1.9), (2,4.2), (3,5.8), (4,8.1)$. Fit $I = a_0 + a_1V$.
>
> **Normal equations.** $\sum1 = 5$, $\sum V_i = 10$, $\sum V_i^2 = 30$, $\sum I_i = 20.1$, $\sum V_iI_i = 0 + 1.9 + 8.4 + 17.4 + 32.4 = 60.1$:
> $$
> \begin{bmatrix}5&10\\10&30\end{bmatrix}\begin{bmatrix}a_0\\a_1\end{bmatrix} = \begin{bmatrix}20.1\\60.1\end{bmatrix}.
> $$
>
> **Solve.** Determinant $150 - 100 = 50$; $a_1 = (5\cdot60.1 - 10\cdot20.1)/50 = (300.5 - 201)/50 = 1.99$ mA/V, $a_0 = (20.1 - 10a_1)/5 = 0.04$ mA.
>
> **Interpretation.** Conductance $1.99$ mS, i.e. $R = 1/1.99\ \mathrm{mS} = 503\ \Omega$; offset $0.04$ mA (instrument offset). Residuals $0.06, -0.13, 0.18, -0.21, 0.10$ sum to zero (as orthogonality to the constant column requires) and their sum of squares is $0.107$, giving an estimated noise standard deviation $\sqrt{0.107/(5-2)} = 0.19$ mA (Ch. 20).

**Numerical caution.** Forming $\mathbf A^\mathsf T\mathbf A$ *squares* the condition number: $\kappa(\mathbf A^\mathsf T\mathbf A) = \kappa(\mathbf A)^2$. For the quadratic fit to the data above, $\kappa(\mathbf A) = 27$ but $\kappa(\mathbf A^\mathsf T\mathbf A) = 735$; for high-degree monomial fits the squaring is catastrophic. Stable practice is to solve the least-squares problem by a **QR factorisation** $\mathbf A = \mathbf Q\mathbf R$ (then $\mathbf R\mathbf a = \mathbf Q^\mathsf T\mathbf y$) or the SVD, which is what `numpy.linalg.lstsq` and MATLAB's backslash do; and to **centre and scale** the abscissae (use $(x - \bar x)/s$) before fitting polynomials.

### Continuous least squares and orthogonal polynomials

To approximate a *function* $f$ on $[a,b]$ in the $L^2$ sense with $p(x) = \sum_k a_kx^k$, minimise $\int_a^b(f-p)^2dx$. The normal equations become

$$
\sum_{k=0}^{n}\left(\int_a^b x^{j+k}dx\right)a_k = \int_a^b x^jf(x)\,dx,\qquad j = 0,\ldots,n.
$$

On $[0,1]$ the matrix is the **Hilbert matrix** $H_{jk} = 1/(j+k+1)$ — which we met in Section 22.2 as the archetype of ill-conditioning. The monomials are a terrible basis: $x^{10}$ and $x^{11}$ look almost alike on $[0,1]$.

The cure is an **orthogonal basis**. If $\langle\phi_j,\phi_k\rangle = \int_a^bw\,\phi_j\phi_k\,dx = 0$ for $j\neq k$, the normal-equation matrix is diagonal and the coefficients decouple:

$$
a_k = \frac{\langle f,\phi_k\rangle}{\langle\phi_k,\phi_k\rangle}.
$$

This is exactly the Fourier-coefficient formula of Ch. 15 with polynomials in place of sinusoids — and, as there, adding a term does not change the previous coefficients, and the error decreases monotonically as terms are added (Bessel's inequality, Parseval).

| Family | Interval | Weight $w(x)$ | Recurrence | Norm $\langle\phi_k,\phi_k\rangle$ |
|---|---|---|---|---|
| Legendre $P_k$ | $[-1,1]$ | $1$ | $(k+1)P_{k+1} = (2k+1)xP_k - kP_{k-1}$ | $2/(2k+1)$ |
| Chebyshev $T_k$ | $[-1,1]$ | $(1-x^2)^{-1/2}$ | $T_{k+1} = 2xT_k - T_{k-1}$ | $\pi$ ($k=0$), $\pi/2$ ($k\ge1$) |
| Laguerre $L_k$ | $[0,\infty)$ | $e^{-x}$ | $(k+1)L_{k+1} = (2k+1-x)L_k - kL_{k-1}$ | $1$ |
| Hermite $H_k$ | $(-\infty,\infty)$ | $e^{-x^2}$ | $H_{k+1} = 2xH_k - 2kH_{k-1}$ | $\sqrt\pi\,2^kk!$ |

The first Legendre polynomials are $P_0 = 1$, $P_1 = x$, $P_2 = \tfrac12(3x^2-1)$, $P_3 = \tfrac12(5x^3-3x)$. Any family of orthogonal polynomials can be generated from the monomials by Gram–Schmidt (Ch. 9), and every such family satisfies a three-term recurrence.

> [!example] Worked Example 22.23 — Least-squares straight line for $e^x$ on $[0,1]$
> **Normal equations (monomial basis).**
> $$
> \begin{bmatrix}1&\tfrac12\\ \tfrac12&\tfrac13\end{bmatrix}\begin{bmatrix}a_0\\a_1\end{bmatrix} = \begin{bmatrix}\int_0^1e^x dx\\ \int_0^1xe^xdx\end{bmatrix} = \begin{bmatrix}e-1\\1\end{bmatrix}.
> $$
>
> **Solve.** Determinant $\tfrac13-\tfrac14 = \tfrac1{12}$; $a_0 = 12\left(\tfrac{e-1}{3} - \tfrac12\right) = 4e - 10 = 0.873127$, $a_1 = 12\left(1 - \tfrac{e-1}{2}\right) = 18 - 6e = 1.690309$.
>
> **Compare with the minimax line of Example 22.21.**
>
> | Approximation | $\lVert e^x - p\rVert_\infty$ | $\lVert e^x - p\rVert_2$ |
> |---|---|---|
> | Taylor $1 + x$ | 0.718 | — |
> | Least squares $0.8731 + 1.6903x$ | 0.155 | **0.0628** |
> | Minimax $0.8941 + 1.7183x$ | **0.106** | 0.0723 |
>
> Each is optimal in its own norm, and they are close to each other — a general phenomenon for smooth functions.

> [!example] Worked Example 22.24 — Legendre expansion of $e^x$ on $[-1,1]$
> **Coefficients.** $a_k = \frac{2k+1}{2}\int_{-1}^{1}e^xP_k(x)\,dx$:
> $$
> a_0 = \tfrac12(e - e^{-1}) = 1.175201,\quad a_1 = \tfrac32\cdot2e^{-1} = 1.103638,\quad a_2 = \tfrac52\,(e - 7e^{-1}) = 0.357814.
> $$
> (The integrals follow by repeated integration by parts: $\int_{-1}^1xe^x dx = 2/e$, $\int_{-1}^1x^2e^xdx = e - 5/e$.)
>
> **Quadratic approximation.** $p_2 = a_0 + a_1x + a_2\cdot\tfrac12(3x^2-1) = 0.996294 + 1.103638x + 0.536722x^2$.
>
> **Accuracy.** $\max_{[-1,1]}\lvert e^x - p_2\rvert = 0.0816$, compared with $0.218$ for the Taylor quadratic $1 + x + x^2/2$. Adding $a_3P_3$ would not change $a_0,a_1,a_2$.

### Chebyshev series and economisation

Expanding in Chebyshev polynomials (weight $(1-x^2)^{-1/2}$) gives a series $f = \sum c_kT_k$ whose truncation is **near-minimax**: since $\lvert T_k\rvert\le1$, the truncation error is at most $\sum_{k>n}\lvert c_k\rvert$, and because the first omitted term $c_{n+1}T_{n+1}$ usually dominates and equioscillates $n+2$ times, the truncated series is very close to the true minimax polynomial. Interpolating at Chebyshev nodes (Section 22.5) gives almost the same polynomial and is easier to compute (by FFT).

**Economisation** turns a Taylor polynomial into a better approximation of lower degree: express the highest power through Chebyshev polynomials and drop the $T_{n}$ term, whose magnitude is known exactly.

> [!example] Worked Example 22.25 — Economising the cubic Taylor polynomial of $e^x$
> **Step 1.** $e^x\approx1 + x + \tfrac12x^2 + \tfrac16x^3$ on $[-1,1]$ (max error $0.0516$).
>
> **Step 2.** From $T_3 = 4x^3 - 3x$: $x^3 = \tfrac14(T_3 + 3x)$, so $\tfrac16x^3 = \tfrac1{24}T_3 + \tfrac18x$.
>
> **Step 3.** Drop $\tfrac{1}{24}T_3$ (error at most $1/24 = 0.0417$): $e^x\approx1 + 1.125x + 0.5x^2$.
>
> **Accuracy.** A numerical scan gives $\max\lvert e^x - (1 + 1.125x + 0.5x^2)\rvert = 0.093$, close to the predicted $0.0516 + 0.0417$ and far better than the Taylor quadratic ($0.218$). For comparison, interpolating at the three Chebyshev nodes gives $1 + 1.12977x + 0.53204x^2$ with error $0.056$.

### Choosing the degree

In data fitting, the residual always decreases as the degree increases, until at degree $N-1$ the fit interpolates and the residual is zero — while the polynomial oscillates wildly between data points. This is **overfitting**: the extra parameters model the noise. Choose the degree by physics (a resistor is linear; a thermistor's $\ln R$ is nearly linear in $1/T$), by checking that residuals look like white noise, or by statistical criteria (adjusted $R^2$, cross-validation, the $F$-test of Ch. 20). Weighted least squares (weights $w_i = 1/\sigma_i^2$) is appropriate when measurements have unequal uncertainties.

> [!warning] Common Pitfall — Fitting in the wrong variables
> Least squares minimises errors in the variable you fit. Linearising an exponential model $y = Ae^{bx}$ by fitting $\ln y$ versus $x$ is convenient but changes the error weighting — small $y$ values (often the noisiest) get magnified influence. If the noise is additive in $y$, either weight the log-fit by $y_i^2$ or fit the nonlinear model directly (Gauss–Newton, Ch. 23).

> [!tip] Engineering Note — Polynomial approximation in hardware and firmware
> The `sin`, `exp` and `log` routines in a math library (or an FPGA's CORDIC alternative) typically reduce the argument to a small interval and then evaluate a **minimax** polynomial or rational function of modest degree, coefficients precomputed by the Remez algorithm. Sensor linearisation in firmware — thermocouple voltage to temperature (the NIST ITS-90 polynomials are degree 8–10), strain-gauge calibration, RTD Callendar–Van Dusen curves — uses **least-squares** polynomial fits, evaluated by Horner's scheme to minimise multiplications and rounding error.

## 22.7 Numerical differentiation

### Motivation

Derivatives are needed numerically whenever a function is known only through samples (a digitised waveform, a measured C–V curve whose slope gives doping profile), or through a black-box program (the Jacobian for Newton's method, sensitivities $\partial V_{\text{out}}/\partial R_k$ of a simulated circuit), and — most importantly — as the building block for solving differential equations on grids (Sections 22.9–22.10). Numerical differentiation is easy to do and easy to do badly: it is the one basic numerical operation in which **truncation and rounding errors pull in opposite directions**.

### Finite-difference formulas from Taylor series

Let $f$ be smooth and $h>0$ a step. Taylor expansions about $x$ give

$$
\begin{aligned}
f(x+h) &= f + hf' + \tfrac{h^2}{2}f'' + \tfrac{h^3}{6}f''' + \tfrac{h^4}{24}f^{(4)} + \cdots\\
f(x-h) &= f - hf' + \tfrac{h^2}{2}f'' - \tfrac{h^3}{6}f''' + \tfrac{h^4}{24}f^{(4)} - \cdots
\end{aligned}
$$

(all derivatives evaluated at $x$). Rearranging the first line:

$$
\textbf{forward difference:}\quad D_+f(x) = \frac{f(x+h)-f(x)}{h} = f'(x) + \frac h2f''(\xi),
$$

so the truncation error is $O(h)$ — first order. The **backward difference** $D_-f = (f(x)-f(x-h))/h$ is likewise first order. Subtracting the two expansions cancels the even terms:

$$
\textbf{central difference:}\quad D_0f(x) = \frac{f(x+h)-f(x-h)}{2h} = f'(x) + \frac{h^2}{6}f'''(\xi),
$$

which is **second order**. Adding them cancels the odd terms and gives the standard second-derivative formula:

$$
D_+D_-f(x) = \frac{f(x+h)-2f(x)+f(x-h)}{h^2} = f''(x) + \frac{h^2}{12}f^{(4)}(\xi).
$$

The general pattern: a symmetric stencil kills alternate error terms, gaining one order for free.

> [!definition] Order of a difference approximation
> A difference formula $D_hf$ approximates $f^{(m)}(x)$ with **order $p$** if $D_hf - f^{(m)}(x) = Ch^pf^{(m+p)}(\xi)$ with $C\neq0$ for all smooth $f$. The quantity $D_hf - f^{(m)}(x)$ is the **truncation error**; $C$ is its error constant.

| Approximation | Formula | Truncation error |
|---|---|---|
| $f'$ forward | $\dfrac{f_1 - f_0}{h}$ | $\tfrac{h}{2}f''$ |
| $f'$ backward | $\dfrac{f_0 - f_{-1}}{h}$ | $-\tfrac{h}{2}f''$ |
| $f'$ central | $\dfrac{f_1 - f_{-1}}{2h}$ | $\tfrac{h^2}{6}f'''$ |
| $f'$ one-sided, 2nd order | $\dfrac{-3f_0 + 4f_1 - f_2}{2h}$ | $-\tfrac{h^2}{3}f'''$ |
| $f'$ central, 4th order | $\dfrac{-f_2 + 8f_1 - 8f_{-1} + f_{-2}}{12h}$ | $-\tfrac{h^4}{30}f^{(5)}$ |
| $f''$ central | $\dfrac{f_1 - 2f_0 + f_{-1}}{h^2}$ | $\tfrac{h^2}{12}f^{(4)}$ |
| $f''$ central, 4th order | $\dfrac{-f_2 + 16f_1 - 30f_0 + 16f_{-1} - f_{-2}}{12h^2}$ | $-\tfrac{h^4}{90}f^{(6)}$ |
| $f'''$ central | $\dfrac{f_2 - 2f_1 + 2f_{-1} - f_{-2}}{2h^3}$ | $\tfrac{h^2}{4}f^{(5)}$ |

Here $f_k = f(x + kh)$. The one-sided second-order formula is what one uses at a boundary, where no point to the left exists (Neumann boundary conditions in Section 22.10).

**Method of undetermined coefficients.** Any stencil can be derived systematically. To approximate $f'(0)$ from $f_0,f_1,f_2$, write $f'(0)\approx(af_0 + bf_1 + cf_2)/h$, expand each $f_k$ in Taylor series, and demand that the coefficients of $f$, $f'$, $f''$ be $0,1,0$:

$$
a+b+c = 0,\qquad b + 2c = 1,\qquad \tfrac12b + 2c = 0\quad\Longrightarrow\quad (a,b,c) = \left(-\tfrac32,\ 2,\ -\tfrac12\right),
$$

which is the one-sided formula in the table. Equivalently, one can differentiate the interpolating polynomial through the stencil points (Section 22.5): the weights are $\ell_i'(x)$. That viewpoint generalises to unequally spaced points.

> [!example] Worked Example 22.26 — Derivatives of $e^x$ at $x = 0$
> Exact values: $f'(0) = f''(0) = 1$.
>
> | $h$ | forward | error | central | error | $D_+D_-$ for $f''$ | error |
> |---|---|---|---|---|---|---|
> | 0.1 | 1.0517092 | $5.2\times10^{-2}$ | 1.0016675 | $1.67\times10^{-3}$ | 1.0008336 | $8.3\times10^{-4}$ |
> | 0.05 | 1.0254219 | $2.5\times10^{-2}$ | 1.0004167 | $4.17\times10^{-4}$ | 1.0002084 | $2.1\times10^{-4}$ |
>
> **Order check.** Halving $h$ halves the forward error (ratio 2.0, order 1) and quarters the central and second-difference errors (ratio 4.0, order 2).
>
> **Error constants.** Predicted forward error $\tfrac h2f'' = 0.05$ ✓; central $\tfrac{h^2}6f''' = 1.67\times10^{-3}$ ✓; second difference $\tfrac{h^2}{12}f^{(4)} = 8.3\times10^{-4}$ ✓.
>
> **Higher order.** At $h = 0.1$ the one-sided second-order formula gives an error of $-3.6\times10^{-3}$ (predicted $-\tfrac{h^2}{3} = -3.3\times10^{-3}$) — much better than the forward difference's $5.2\times10^{-2}$ with the same points available on one side — and the five-point fourth-order formula gives $-3.3\times10^{-6}$ (predicted $-h^4/30 = -3.3\times10^{-6}$).

### Richardson extrapolation

If an approximation has an error expansion in known powers of $h$,

$$
D(h) = L + c_1h^p + c_2h^{q} + \cdots\qquad(q>p),
$$

then combining two step sizes eliminates the leading error term:

$$
D(h) - L\approx c_1h^p,\quad D(h/2) - L\approx c_1(h/2)^p\quad\Longrightarrow\quad L\approx R(h) = \frac{2^pD(h/2) - D(h)}{2^p - 1},
$$

with error $O(h^{q})$. For the central difference, $p = 2$ and the expansion contains only even powers, so $R(h) = \frac{4D_0(h/2) - D_0(h)}{3}$ is **fourth order**. Repeating the process builds a triangular table (the same idea gives Romberg integration in Section 22.8). Richardson extrapolation also provides an **error estimate** for free: $D(h/2) - L\approx\frac{D(h/2)-D(h)}{2^p-1}$.

> [!example] Worked Example 22.27 — Richardson for the derivative of $e^x$
> From Example 22.26, $D_0(0.1) = 1.00166750$ and $D_0(0.05) = 1.00041672$.
>
> **Extrapolate.** $R = \dfrac{4(1.00041672) - 1.00166750}{3} = 0.99999979$, error $-2.1\times10^{-7}$ — an improvement by a factor of 2000 over $D_0(0.05)$, using no new function values.
>
> **Error estimate.** $(D_0(0.05) - D_0(0.1))/3 = -4.17\times10^{-4}$, so $D_0(0.05)\approx L + 4.17\times10^{-4}$: the estimate matches the actual error of $D_0(0.05)$ to three digits.
>
> (Algebraically $R(h)$ equals the five-point formula with step $h/2$, whose predicted error $-(0.05)^4/30 = -2.1\times10^{-7}$ agrees.)

### Rounding error: the trade-off

The formulas above assume exact function values. In floating point, $f(x\pm h)$ carries an absolute error of about $u\lvert f\rvert$ (more if $f$ is computed by a long program). For the forward difference these errors are divided by $h$, giving a total error

$$
E(h)\approx\underbrace{\frac{h}{2}\lvert f''\rvert}_{\text{truncation}} + \underbrace{\frac{2u\lvert f\rvert}{h}}_{\text{rounding}} .
$$

The first term falls with $h$, the second rises. Minimising over $h$ ($E'(h) = 0$):

$$
h_{\text{opt}} = 2\sqrt{\frac{u\lvert f\rvert}{\lvert f''\rvert}}\sim\sqrt{u}\approx10^{-8},\qquad E_{\min} = 2\sqrt{u\lvert f\rvert\lvert f''\rvert}\sim\sqrt u\approx10^{-8}.
$$

So in double precision **a forward difference can never deliver more than about 8 correct digits.** For the central difference, $E(h)\approx\tfrac{h^2}{6}\lvert f'''\rvert + u\lvert f\rvert/h$, giving $h_{\text{opt}} = (3u\lvert f\rvert/\lvert f'''\rvert)^{1/3}\sim u^{1/3}\approx 5\times10^{-6}$ and $E_{\min}\sim u^{2/3}\approx10^{-11}$. Higher-order formulas push the optimum to larger $h$ and smaller minimum error, but every formula for an $m$-th derivative suffers rounding amplification like $u/h^m$.

> [!example] Worked Example 22.28 — The error "V" for $\frac{d}{dx}\sin x$ at $x = 1$
> Exact: $\cos1 = 0.540302306$.
>
> | $h$ | $10^{-1}$ | $10^{-2}$ | $10^{-4}$ | $10^{-5}$ | $10^{-6}$ | $10^{-8}$ | $10^{-10}$ | $10^{-12}$ |
> |---|---|---|---|---|---|---|---|---|
> | forward error | $4.3\times10^{-2}$ | $4.2\times10^{-3}$ | $4.2\times10^{-5}$ | $4.2\times10^{-6}$ | $4.2\times10^{-7}$ | $3.0\times10^{-9}$ | $5.8\times10^{-8}$ | $4.3\times10^{-5}$ |
> | central error | $9.0\times10^{-4}$ | $9.0\times10^{-6}$ | $9.0\times10^{-10}$ | $1.1\times10^{-11}$ | $2.8\times10^{-11}$ | $2.6\times10^{-9}$ | $5.8\times10^{-8}$ | $1.2\times10^{-5}$ |
>
> **Reading the table.** Both errors first fall at their theoretical rates (one decade per decade of $h$ for forward, two for central). The forward error bottoms out near $h = 10^{-8}$ at about $3\times10^{-9}$ (theory: $h_{\text{opt}} = 2\sqrt u = 2.1\times10^{-8}$, $E_{\min}\approx\sqrt{u}\cdot2\sin1\approx 2\times10^{-8}$; the observed value is smaller because rounding errors partially cancel). The central error bottoms out near $h = 10^{-5}$ at about $10^{-11}$ (theory: $h_{\text{opt}} = (3u\tan 1)^{1/3} = 8\times10^{-6}$). For smaller $h$ both curves *rise* along the same line, slope $-1$, because rounding dominates.

```desmos-graph
left=-13; right=0; bottom=-12; top=0
---
y=\log(0.4207\cdot10^{x}+1.87\cdot10^{-16}\cdot10^{-x})|red
y=\log(0.09005\cdot10^{2x}+0.935\cdot10^{-16}\cdot10^{-x})|blue
(-1,-1.367)|red
(-2,-2.375)|red
(-4,-4.376)|red
(-6,-6.376)|red
(-8,-8.527)|red
(-10,-7.233)|red
(-12,-4.364)|red
(-1,-3.046)|blue
(-2,-5.046)|blue
(-4,-9.046)|blue
(-5,-10.953)|blue
(-6,-10.557)|blue
(-8,-8.588)|blue
(-12,-4.911)|blue
```

*Total error of the forward (red) and central (blue) difference for $(\sin x)'$ at $x=1$, on log–log axes ($x=\log_{10}h$, $y=\log_{10}$ error): model curves $\tfrac h2\lvert f''\rvert+2u\lvert f\rvert/h$ and $\tfrac{h^2}{6}\lvert f'''\rvert+u\lvert f\rvert/h$ with measured errors as dots. Truncation dominates on the right (slopes 1 and 2), rounding on the left (slope $-1$).*

### The complex-step derivative

Engineers comfortable with $j$ have a neat escape from the rounding trade-off when $f$ is analytic and implemented in complex arithmetic. Expand with an *imaginary* step:

$$
f(x + jh) = f(x) + jhf'(x) - \tfrac{h^2}{2}f''(x) - j\tfrac{h^3}{6}f'''(x) + \cdots
\quad\Longrightarrow\quad
f'(x) = \frac{\operatorname{Im}\{f(x+jh)\}}{h} + \frac{h^2}{6}f'''(x) + \cdots
$$

There is **no subtraction**, hence no cancellation, so $h$ can be taken absurdly small: with $h = 10^{-20}$ the truncation error $\sim10^{-40}$ is invisible and the result for $(\sin x)'$ at $x = 1$ is correct to the last bit (computed error exactly $0$). The method is widely used for accurate sensitivities in optimisation codes; it requires that the code for $f$ be "complex-safe" (no `abs`, no comparisons on real parts that break analyticity).

### Differentiating noisy data

Measured data contain noise of size $\sigma$ that does not decrease with $h$. The rounding analysis then applies with $u\lvert f\rvert$ replaced by $\sigma$: a central difference of samples spaced $h$ apart has noise standard deviation $\sigma/(\sqrt2\,h)$, and a second difference $\sqrt6\,\sigma/h^2$. Differentiation is a **high-pass** operation — its frequency response is $j\omega$ (Ch. 16) — so it amplifies high-frequency noise. Practical remedies smooth before (or while) differentiating: fit a local low-degree least-squares polynomial and differentiate it (the **Savitzky–Golay** filter), differentiate a smoothing spline, or use a band-limited differentiator FIR filter.

> [!warning] Common Pitfall — "Smaller $h$ is always more accurate"
> It is not. Below the optimum step, every halving of $h$ *doubles* the rounding error (for a first derivative) or quadruples it (for a second). A finite-difference Jacobian in a Newton solver with $h = 10^{-14}\lvert x\rvert$ can be pure noise and stall the iteration. A sound default is $h = \sqrt{u}\,\max(\lvert x\rvert, x_{\text{typ}})$ for forward differences and $h = u^{1/3}\max(\lvert x\rvert,x_{\text{typ}})$ for central differences, with $h$ rounded so that $x+h$ is exactly representable (`h = (x + h) - x`).

> [!tip] Engineering Note — Differentiators in practice
> An oscilloscope's "derivative" math channel, a PID controller's D term and a digital tachometer all implement numerical differentiation of noisy samples. That is why PID implementations filter the derivative term (a first-order low-pass on $D$, "derivative kick" avoidance) and why rotational speed is often estimated by counting encoder pulses over a longer window (which is a long-baseline difference, trading bandwidth for noise) rather than by differentiating successive samples.

## 22.8 Numerical integration

### Motivation

Definite integrals without elementary antiderivatives are everywhere: the error function in bit-error-rate calculations (Ch. 19), the RMS value of a measured waveform, the energy of a pulse, Fourier coefficients of sampled signals (Ch. 15), the inductance of a loop from Neumann's formula, radiation patterns of antennas, and integrals of tabulated data. **Quadrature** (the classical name for numerical integration) approximates

$$
I(f) = \int_a^b f(x)\,dx\qquad\text{by}\qquad Q(f) = \sum_{i=0}^{n}w_if(x_i),
$$

a weighted sum of function values at **nodes** $x_i$ with **weights** $w_i$. All the methods of this section are of this form; they differ in how nodes and weights are chosen.

> [!definition] Degree of exactness
> A quadrature rule has **degree of exactness** (or degree of precision) $d$ if $Q(p) = I(p)$ for every polynomial $p$ of degree $\le d$, but not for some polynomial of degree $d+1$.

Since $Q$ and $I$ are both linear, it suffices to check the monomials $1,x,\ldots,x^d$. A rule obtained by integrating the interpolating polynomial through $n+1$ nodes — an **interpolatory rule**, with weights $w_i = \int_a^b\ell_i(x)\,dx$ — has degree at least $n$.

### Newton–Cotes rules

Interpolatory rules on **equally spaced** nodes are called Newton–Cotes rules.

**Trapezoidal rule** (linear interpolation through the endpoints):

$$
\int_a^b f\,dx\approx\frac{b-a}{2}\big(f(a)+f(b)\big),\qquad E = -\frac{(b-a)^3}{12}f''(\xi).
$$

*Derivation of the error.* By the interpolation error formula (Section 22.5), $f - p_1 = \tfrac12f''(\xi_x)(x-a)(x-b)$. Since $(x-a)(x-b)\le0$ on $[a,b]$ does not change sign, the integral mean-value theorem gives $\int_a^b(f - p_1)\,dx = \tfrac12f''(\xi)\int_a^b(x-a)(x-b)\,dx = \tfrac12f''(\xi)\left(-\tfrac{(b-a)^3}{6}\right)$. $\square$

**Midpoint rule** $\int_a^bf\approx(b-a)f\!\left(\tfrac{a+b}2\right)$ has error $+\tfrac{(b-a)^3}{24}f''(\xi)$ — half the trapezoidal error and of opposite sign, a fact exploited by Simpson's rule.

**Simpson's rule** (quadratic interpolation through $a$, $m = \tfrac{a+b}2$, $b$; $h = \tfrac{b-a}{2}$):

$$
\int_a^b f\,dx\approx\frac h3\big(f(a) + 4f(m) + f(b)\big),\qquad E = -\frac{h^5}{90}f^{(4)}(\xi).
$$

Although it interpolates with a quadratic, Simpson's rule integrates **cubics exactly** (degree of exactness 3): by symmetry about $m$, the rule integrates $(x-m)^3$ exactly (both sides are zero). This "bonus degree" occurs for every symmetric rule with an odd number of points. Simpson's rule can also be seen as the combination $S = \tfrac23M + \tfrac13T$ of midpoint and trapezoid, which cancels their $f''$ errors.

Higher-order Newton–Cotes rules exist (Simpson's 3/8 rule, Boole's rule), but beyond about 8 points their weights become negative and large — the integration counterpart of the Runge phenomenon — so in practice one uses low-order rules **composed** over many subintervals.

### Composite rules

Divide $[a,b]$ into $n$ subintervals of width $h = (b-a)/n$, nodes $x_i = a + ih$, and apply the basic rule on each:

$$
\begin{aligned}
T_n &= h\left[\tfrac12f_0 + f_1 + f_2 + \cdots + f_{n-1} + \tfrac12f_n\right], & E_T &= -\frac{(b-a)h^2}{12}f''(\xi),\\
S_n &= \frac h3\left[f_0 + 4f_1 + 2f_2 + 4f_3 + \cdots + 2f_{n-2} + 4f_{n-1} + f_n\right]\ (n\text{ even}), & E_S &= -\frac{(b-a)h^4}{180}f^{(4)}(\xi).
\end{aligned}
$$

The composite error follows by summing the $n$ local errors $-\tfrac{h^3}{12}f''(\xi_i)$ and replacing the average $\frac1n\sum f''(\xi_i)$ by $f''(\xi)$ (intermediate value theorem). Each local error is one power of $h$ higher than the global error: $n\propto1/h$ local errors add up. So the composite trapezoidal rule is **second order** and composite Simpson **fourth order**.

```python
import numpy as np
def trapezoid(f, a, b, n):
    x = np.linspace(a, b, n + 1); y = f(x)
    return (b - a) / n * (y.sum() - 0.5 * (y[0] + y[-1]))
def simpson(f, a, b, n):                     # n must be even
    x = np.linspace(a, b, n + 1); y = f(x); h = (b - a) / n
    return h / 3 * (y[0] + y[-1] + 4 * y[1:-1:2].sum() + 2 * y[2:-1:2].sum())
```

> [!example] Worked Example 22.29 — $\int_0^1e^{-x^2}dx$ by trapezoid and Simpson
> The exact value is $\tfrac{\sqrt\pi}{2}\operatorname{erf}(1) = 0.746824132812$.
>
> **Trapezoid, $n = 4$** ($h = 0.25$). $f$ values: $1,\ 0.939413,\ 0.778801,\ 0.569783,\ 0.367879$.
> $T_4 = 0.25\,[0.5 + 0.939413 + 0.778801 + 0.569783 + 0.183940] = 0.742984$.
>
> **Simpson, $n = 4$.** $S_4 = \tfrac{0.25}{3}\,[1 + 4(0.939413) + 2(0.778801) + 4(0.569783) + 0.367879] = 0.746855$.
>
> | $n$ | $T_n$ | error | $S_n$ | error |
> |---|---|---|---|---|
> | 2 | 0.731370 | $-1.55\times10^{-2}$ | 0.747180 | $3.6\times10^{-4}$ |
> | 4 | 0.742984 | $-3.84\times10^{-3}$ | 0.746855 | $3.1\times10^{-5}$ |
> | 8 | 0.745866 | $-9.59\times10^{-4}$ | 0.746826 | $2.0\times10^{-6}$ |
>
> **Orders.** Trapezoid errors fall by 4.0 per halving ($h^2$); Simpson by $11.4$ and $15.7$, approaching 16 ($h^4$).
>
> **Error prediction.** The Euler–Maclaurin formula below predicts $T_n - I\approx\tfrac{h^2}{12}\big(f'(1) - f'(0)\big) = \tfrac{h^2}{12}(-2e^{-1})$, i.e. $-3.832\times10^{-3}$ for $n=4$, matching the observed $-3.840\times10^{-3}$.

### The Euler–Maclaurin formula and Romberg integration

> [!theorem] Euler–Maclaurin formula
> If $f\in C^{2m+2}[a,b]$, the composite trapezoidal rule satisfies
> $$
> T_n - \int_a^bf\,dx = \sum_{k=1}^{m}\frac{B_{2k}h^{2k}}{(2k)!}\Big(f^{(2k-1)}(b) - f^{(2k-1)}(a)\Big) + O(h^{2m+2}),
> $$
> where $B_2 = \tfrac16$, $B_4 = -\tfrac1{30}$, $B_6 = \tfrac{1}{42},\ldots$ are Bernoulli numbers. In particular the error has an expansion in **even** powers of $h$: $T_n = I + c_1h^2 + c_2h^4 + c_3h^6 + \cdots$.

Two important consequences:

1. **Romberg integration.** Since the error expansion is in powers of $h^2$, Richardson extrapolation (Section 22.7) applies repeatedly. With $R_{k,0} = T_{2^k}$,
$$
R_{k,m} = R_{k,m-1} + \frac{R_{k,m-1} - R_{k-1,m-1}}{4^m - 1}.
$$
The column $R_{k,1}$ is exactly composite Simpson; $R_{k,2}$ is composite Boole (order 6); and so on.

2. **Periodic functions.** If $f$ is smooth and periodic with period $b - a$, all the derivative differences $f^{(2k-1)}(b) - f^{(2k-1)}(a)$ vanish, so the trapezoidal error is smaller than *every* power of $h$. For analytic periodic functions the error decays **geometrically** in $n$. The humble trapezoidal rule is the best possible method for integrating smooth periodic functions over a full period.

> [!example] Worked Example 22.30 — Romberg table for $\int_0^1e^{-x^2}dx$
> | $k$ | $n = 2^k$ | $R_{k,0}$ (trapezoid) | $R_{k,1}$ (Simpson) | $R_{k,2}$ | $R_{k,3}$ |
> |---|---|---|---|---|---|
> | 0 | 1 | 0.6839397 | | | |
> | 1 | 2 | 0.7313703 | 0.7471804 | | |
> | 2 | 4 | 0.7429841 | 0.7468554 | 0.7468337 | |
> | 3 | 8 | 0.7458656 | 0.7468261 | 0.74682417 | 0.74682402 |
>
> For example $R_{1,1} = 0.7313703 + (0.7313703 - 0.6839397)/3 = 0.7471804$. With 9 function evaluations, $R_{3,2}$ is correct to $4\times10^{-8}$, versus $10^{-3}$ for $T_8$ alone. The last entry $R_{3,3}$ ($-1.1\times10^{-7}$) is slightly worse than $R_{3,2}$: extrapolation assumes the asymptotic regime, and the crude $R_{0,0}$ (one interval!) contaminates the far corner of the table. Practical Romberg codes stop when successive diagonal entries agree to tolerance.

> [!example] Worked Example 22.31 — Spectral accuracy for a periodic integrand
> $I = \int_0^{2\pi}\frac{dx}{2 + \cos x} = \frac{2\pi}{\sqrt3} = 3.6275987285$ (residue calculus, Ch. 18). For a periodic function the trapezoidal rule with $n$ points is simply $\tfrac{2\pi}{n}\sum_{i=0}^{n-1}f(2\pi i/n)$:
>
> | $n$ | 2 | 4 | 8 | 16 |
> |---|---|---|---|---|
> | $T_n$ | 4.18879 | 3.66519 | 3.62779 | 3.6275987336 |
> | error | $5.6\times10^{-1}$ | $3.8\times10^{-2}$ | $1.9\times10^{-4}$ | $5.1\times10^{-9}$ |
>
> Each doubling of $n$ *squares* the error, roughly — geometric convergence (the error behaves like $(2-\sqrt3)^n\approx 0.268^n$, set by the location of the integrand's poles at $\cos x = -2$). A fourth-order method would only gain a factor of 16 per doubling.

### Gaussian quadrature

Newton–Cotes rules fix the nodes and choose $n+1$ weights, achieving degree $n$ (or $n+1$). If the **nodes are free too**, there are $2n+2$ parameters and one can hope for degree $2n+1$. Gauss showed this is achievable, and optimal.

> [!theorem] Gauss–Legendre quadrature
> For each $n\ge1$ there is a unique $n$-point rule
> $$
> \int_{-1}^{1}f(x)\,dx\approx\sum_{i=1}^{n}w_if(x_i)
> $$
> with degree of exactness $2n - 1$; no $n$-point rule has higher degree. Its nodes are the zeros of the Legendre polynomial $P_n$, all real, distinct and in $(-1,1)$, and its weights $w_i = \int_{-1}^1\ell_i(x)\,dx = \dfrac{2}{(1-x_i^2)[P_n'(x_i)]^2}$ are all positive. The error is
> $$
> E_n = \frac{2^{2n+1}(n!)^4}{(2n+1)[(2n)!]^3}\,f^{(2n)}(\xi).
> $$

*Proof sketch of exactness.* Let $p$ have degree $\le 2n-1$. Divide by $P_n$: $p = qP_n + r$ with $\deg q,\deg r\le n-1$. Then $\int_{-1}^1 p = \int qP_n + \int r = 0 + \int r$, because $P_n$ is orthogonal to all polynomials of lower degree (Section 22.6). At the nodes $P_n(x_i) = 0$, so $\sum w_ip(x_i) = \sum w_ir(x_i) = \int r$, the last equality because the rule is interpolatory with $n$ nodes and $\deg r\le n-1$. *Optimality:* $p = \prod(x-x_i)^2$ has degree $2n$, is positive, and is integrated to $0$ by the rule. $\square$

| $n$ | nodes $x_i$ | weights $w_i$ | exact for degree |
|---|---|---|---|
| 1 | $0$ | $2$ | 1 |
| 2 | $\pm1/\sqrt3 = \pm0.5773503$ | $1,\ 1$ | 3 |
| 3 | $0,\ \pm\sqrt{3/5} = \pm0.7745967$ | $8/9,\ 5/9,\ 5/9$ | 5 |
| 4 | $\pm0.3399810,\ \pm0.8611363$ | $0.6521452,\ 0.3478548$ | 7 |

**Deriving the 2-point rule directly.** Require $w_1f(x_1) + w_2f(x_2)$ to integrate $1, x, x^2, x^3$ exactly on $[-1,1]$: $w_1 + w_2 = 2$, $w_1x_1 + w_2x_2 = 0$, $w_1x_1^2 + w_2x_2^2 = \tfrac23$, $w_1x_1^3 + w_2x_2^3 = 0$. Symmetry suggests $x_2 = -x_1$, $w_1 = w_2 = 1$; then $2x_1^2 = \tfrac23$ gives $x_{1,2} = \pm1/\sqrt3$.

**General interval.** Map $x\in[-1,1]$ to $t = \tfrac{a+b}{2} + \tfrac{b-a}{2}x$:

$$
\int_a^bg(t)\,dt\approx\frac{b-a}{2}\sum_{i=1}^{n}w_i\,g\!\left(\frac{a+b}{2}+\frac{b-a}{2}x_i\right).
$$

> [!example] Worked Example 22.32 — Gauss–Legendre for $\int_0^1e^{-x^2}dx$ and $\int_0^\pi\sin x\,dx$
> **2-point rule on $[0,1]$.** Nodes $t = \tfrac12 \pm\tfrac{1}{2\sqrt3} = 0.211325,\ 0.788675$; weights $\tfrac12$ each.
> $Q_2 = \tfrac12\big(e^{-0.044658} + e^{-0.622008}\big) = \tfrac12(0.956324 + 0.536865) = 0.746595$, error $-2.3\times10^{-4}$.
>
> **3- and 4-point rules.** $Q_3 = 0.7468146$ (error $-9.5\times10^{-6}$), $Q_4 = 0.74682447$ (error $3.4\times10^{-7}$).
>
> **Comparison at equal cost.** $Q_4$ with 4 evaluations beats Simpson $S_8$ with 9 evaluations ($2.0\times10^{-6}$).
>
> **Second integral (exact value 2).** Trapezoid $T_4 = 1.896119$, Simpson $S_4 = 2.004560$, Gauss $Q_2 = 1.935820$, $Q_3 = 2.001389$. Three function values in the right places do better than five equally spaced ones.

**Other Gaussian rules.** The same construction with a weight function $w(x)$ uses the corresponding orthogonal polynomials (Section 22.6): **Gauss–Chebyshev** for $\int_{-1}^1f(x)(1-x^2)^{-1/2}dx$ (nodes are Chebyshev nodes, equal weights $\pi/n$); **Gauss–Laguerre** for $\int_0^\infty f(x)e^{-x}dx$; **Gauss–Hermite** for $\int_{-\infty}^\infty f(x)e^{-x^2}dx$ — the natural rule for expectations over Gaussian noise (Ch. 19), e.g. average bit-error rate over a Gaussian-distributed parameter. Variants that include endpoints as nodes (Gauss–Lobatto, Gauss–Radau) underlie spectral-element field solvers and some implicit Runge–Kutta methods.

### Adaptive quadrature and difficult integrands

A fixed composite rule wastes evaluations where $f$ is smooth and under-resolves where it is not. **Adaptive** quadrature estimates the error on each subinterval — e.g. by comparing Simpson on $[a,b]$ with Simpson on its two halves, $\lvert S_{\text{halves}} - S_{\text{whole}}\rvert/15$ being an error estimate by Richardson's argument — and recursively bisects only the subintervals that fail the tolerance. Library routines (QUADPACK's `qags`, behind `scipy.integrate.quad`; MATLAB `integral`) use adaptive Gauss–Kronrod pairs, which reuse the Gauss nodes in a higher-order Kronrod extension so the error estimate costs few extra evaluations.

> [!example] Worked Example 22.33 — An endpoint singularity
> $\int_0^1\sqrt{x}\,dx = \tfrac23$. The integrand is continuous but $f'$ is infinite at $0$, so the error constants of the composite rules are infinite and convergence degrades: composite Simpson gives errors $-1.3\times10^{-3}$, $-1.6\times10^{-4}$, $-2.0\times10^{-5}$ for $n = 16, 64, 256$ — order $h^{1.5}$ instead of $h^4$. An adaptive Simpson code with tolerance $10^{-8}$ concentrates its subintervals near $x = 0$ and reaches an error of $2.6\times10^{-11}$ with 317 evaluations. Better still is to remove the singularity analytically — here the substitution $x = t^2$ gives $\int_0^12t^2dt$, which Simpson integrates *exactly*.

**Infinite intervals and singularities.** Transform them away when possible: substitutions ($x = t^2$, $x = \tan t$, $x = t/(1-t)$), subtracting off the singular part and integrating it analytically, or Gaussian rules with the singular weight built in. For oscillatory integrands such as Fourier integrals, specialised (Filon, Levin) methods or the FFT are preferable to general-purpose quadrature.

> [!warning] Common Pitfall — Quadrature on data you cannot refine
> When $f$ is a table of measurements, you cannot choose the nodes, so Gaussian quadrature is unavailable, and the "error formulas" involve derivatives of an unknown function. Use composite trapezoid or Simpson (or integrate a spline fit), estimate the error by comparing the rule on all points with the rule on every second point, and remember that noise in the data integrates to a random error of order $\sigma\sqrt{b-a}\sqrt{h}$ — integration, unlike differentiation, *averages* noise down.

> [!tip] Engineering Note — The FFT is a trapezoidal rule
> The DFT coefficient $X_k = \sum_{n=0}^{N-1}x[n]e^{-j2\pi kn/N}$ is, up to the factor $T/N$, the trapezoidal rule for the Fourier-series coefficient $\frac1T\int_0^Tx(t)e^{-jk\omega_0t}dt$ of a periodic signal. Example 22.31 explains why FFT-based spectral analysis of smooth periodic signals is so accurate — and why a *non*-periodic record (a window that does not contain an integer number of periods) loses this accuracy and produces leakage: the endpoint terms of Euler–Maclaurin no longer cancel. RMS meters, power analysers and energy meters integrate sampled $v(t)i(t)$ with exactly these rules.

## 22.9 Numerical solution of ODEs

### The initial-value problem

Most ODEs that describe real systems — nonlinear circuits, motors with saturation, power converters, oscillators, chemical kinetics in batteries — have no closed-form solution. We consider the **initial-value problem (IVP)**

$$
\mathbf y'(t) = \mathbf f(t,\mathbf y),\qquad\mathbf y(t_0) = \mathbf y_0,
$$

where $\mathbf y\in\mathbb R^m$ is the state vector (Ch. 11–12). If $\mathbf f$ is continuous and **Lipschitz** in $\mathbf y$, $\lVert\mathbf f(t,\mathbf y) - \mathbf f(t,\mathbf z)\rVert\le L\lVert\mathbf y - \mathbf z\rVert$, the solution exists and is unique (Picard–Lindelöf). Higher-order equations are reduced to first-order systems: the series RLC equation $L\ddot q + R\dot q + q/C = v(t)$ becomes, with $\mathbf y = (q, i)$, $i = \dot q$,

$$
\frac{d}{dt}\begin{bmatrix}q\\ i\end{bmatrix} = \begin{bmatrix}i\\ \big(v(t) - Ri - q/C\big)/L\end{bmatrix}.
$$

A numerical method produces approximations $\mathbf y_n\approx\mathbf y(t_n)$ on a grid $t_n = t_0 + nh$ (or with variable steps $h_n$).

### Euler's method

Replace the derivative at $t_n$ by a forward difference (Section 22.7):

$$
\mathbf y_{n+1} = \mathbf y_n + h\,\mathbf f(t_n,\mathbf y_n).
$$

Geometrically, follow the tangent of the solution curve through $(t_n,\mathbf y_n)$ for a time $h$. Each step is cheap — one evaluation of $\mathbf f$ — but, as we shall see, inaccurate.

> [!definition] Local and global truncation error
> - The **local truncation error** $\tau_{n+1}$ is the error made in one step starting from the exact solution: for Euler, $\tau_{n+1} = y(t_{n+1}) - y(t_n) - hf(t_n, y(t_n)) = \tfrac{h^2}{2}y''(\xi_n) = O(h^2)$.
> - The **global error** is $e_n = y(t_n) - y_n$, the accumulated error at a fixed time $t_n$.
> - A one-step method is **of order $p$** if its local truncation error is $O(h^{p+1})$; then (under a Lipschitz condition) its global error is $O(h^p)$.

> [!theorem] Convergence of Euler's method
> Let $f$ be Lipschitz in $y$ with constant $L$ and let $\lvert y''\rvert\le M$ on $[t_0,T]$. Then the Euler iterates satisfy
> $$
> \lvert e_n\rvert\le\frac{Mh}{2L}\left(e^{L(t_n - t_0)} - 1\right),\qquad t_0\le t_n\le T,
> $$
> so the global error is $O(h)$: Euler's method is first order.

*Proof.* Subtract the Euler step from the exact relation $y(t_{n+1}) = y(t_n) + hf(t_n,y(t_n)) + \tau_{n+1}$:

$$
e_{n+1} = e_n + h\big[f(t_n,y(t_n)) - f(t_n,y_n)\big] + \tau_{n+1}\ \Rightarrow\ \lvert e_{n+1}\rvert\le(1 + hL)\lvert e_n\rvert + \tfrac{M}{2}h^2.
$$

With $e_0 = 0$, induction gives $\lvert e_n\rvert\le\tfrac M2h^2\sum_{k=0}^{n-1}(1+hL)^k = \tfrac{Mh}{2L}\big((1+hL)^n - 1\big)$, and $(1+hL)^n\le e^{nhL}$. $\square$

The proof displays the general mechanism, which we shall meet again in Section 22.11: **global error = (local errors, $O(h^{p+1})$ each) × (number of steps, $O(1/h)$) × (an amplification factor that stays bounded — stability)**.

> [!example] Worked Example 22.34 — Charging an RC circuit with Euler's method
> A capacitor charges through $R$ from a $5\ \mathrm V$ step, $\tau = RC = 1\ \mathrm{ms}$: $\dot v = (5 - v)/\tau$, $v(0) = 0$, exact $v = 5(1 - e^{-t/\tau})$. Use $h = 0.2\ \mathrm{ms}$.
>
> **Step 1.** $v_1 = 0 + 0.2\cdot(5 - 0)/1 = 1.0\ \mathrm V$ (exact $0.906$).
>
> **Step 2.** $v_2 = 1 + 0.2(5 - 1) = 1.8\ \mathrm V$ (exact $1.648$).
>
> | $t$ (ms) | 0.2 | 0.4 | 0.6 | 0.8 | 1.0 |
> |---|---|---|---|---|---|
> | Euler $v_n$ (V) | 1.000 | 1.800 | 2.440 | 2.952 | 3.362 |
> | exact (V) | 0.906 | 1.648 | 2.256 | 2.753 | 3.161 |
>
> **Interpretation.** Here $v_n = 5\big(1 - (1 - h/\tau)^n\big)$, so Euler replaces $e^{-t/\tau}$ by $(1 - h/\tau)^{t/h} = 0.8^{5} = 0.328$ at $t = \tau$ instead of $0.368$ — a 6 % error in the voltage at one time constant. Halving $h$ halves the error (first order). Notice also that $0<1-h/\tau<1$ is needed for the numerical solution to decay monotonically like the true one; this is a *stability* requirement, discussed below.

### Runge–Kutta methods

Euler uses the slope only at the start of the step. Runge–Kutta (RK) methods sample the slope at several points within the step and combine the samples so that the Taylor expansion of the numerical solution matches the true one to higher order.

**Heun's method (improved Euler, explicit trapezoid), order 2:**

$$
\mathbf k_1 = \mathbf f(t_n,\mathbf y_n),\quad\mathbf k_2 = \mathbf f(t_n + h,\mathbf y_n + h\mathbf k_1),\quad\mathbf y_{n+1} = \mathbf y_n + \tfrac h2(\mathbf k_1 + \mathbf k_2).
$$

**Explicit midpoint method**, also order 2: $\mathbf k_2 = \mathbf f(t_n + \tfrac h2,\mathbf y_n + \tfrac h2\mathbf k_1)$, $\mathbf y_{n+1} = \mathbf y_n + h\mathbf k_2$.

*Why order 2.* For scalar $y' = f(t,y)$, the exact solution satisfies $y(t+h) = y + hf + \tfrac{h^2}{2}(f_t + f_yf) + O(h^3)$. Expanding Heun's $k_2 = f + h(f_t + f_yf) + O(h^2)$ gives $y_{n+1} = y + hf + \tfrac{h^2}{2}(f_t + f_yf) + O(h^3)$ — the two expansions agree through $h^2$, so the local error is $O(h^3)$.

**The classical fourth-order Runge–Kutta method (RK4):**

$$
\begin{aligned}
\mathbf k_1 &= \mathbf f(t_n,\ \mathbf y_n), & \mathbf k_2 &= \mathbf f\!\left(t_n + \tfrac h2,\ \mathbf y_n + \tfrac h2\mathbf k_1\right),\\
\mathbf k_3 &= \mathbf f\!\left(t_n + \tfrac h2,\ \mathbf y_n + \tfrac h2\mathbf k_2\right), & \mathbf k_4 &= \mathbf f\!\left(t_n + h,\ \mathbf y_n + h\mathbf k_3\right),\\
\mathbf y_{n+1} &= \mathbf y_n + \tfrac h6\left(\mathbf k_1 + 2\mathbf k_2 + 2\mathbf k_3 + \mathbf k_4\right).
\end{aligned}
$$

When $\mathbf f$ depends on $t$ only, RK4 reduces to Simpson's rule (weights $\tfrac16,\tfrac46,\tfrac16$, with the two midpoint samples sharing the middle weight). Four function evaluations per step buy fourth-order accuracy: halving $h$ reduces the error by 16 for only twice the work.

```python
def rk4_step(f, t, y, h):
    k1 = f(t, y)
    k2 = f(t + h/2, y + h/2 * k1)
    k3 = f(t + h/2, y + h/2 * k2)
    k4 = f(t + h,   y + h   * k3)
    return y + h/6 * (k1 + 2*k2 + 2*k3 + k4)
```

> [!example] Worked Example 22.35 — One RK4 step, and observed orders
> Solve $y' = y - t^2 + 1$, $y(0) = 0.5$, whose exact solution is $y = (t+1)^2 - \tfrac12e^t$.
>
> **One step with $h = 0.2$.**
> $k_1 = f(0, 0.5) = 1.5$;
> $k_2 = f(0.1,\ 0.5 + 0.1\cdot1.5) = f(0.1, 0.65) = 0.65 - 0.01 + 1 = 1.64$;
> $k_3 = f(0.1,\ 0.5 + 0.1\cdot1.64) = f(0.1, 0.664) = 1.654$;
> $k_4 = f(0.2,\ 0.5 + 0.2\cdot1.654) = f(0.2, 0.8308) = 0.8308 - 0.04 + 1 = 1.7908$;
> $y_1 = 0.5 + \tfrac{0.2}{6}(1.5 + 3.28 + 3.308 + 1.7908) = 0.5 + 0.329293 = 0.829293$.
> Exact: $y(0.2) = 1.44 - 0.5e^{0.2} = 0.829299$. The one-step error is $5.3\times10^{-6}$.
>
> **Global errors at $t = 2$** (exact $y(2) = 5.305472$):
>
> | $h$ | Euler | Heun | RK4 |
> |---|---|---|---|
> | 0.2 | $-0.440$ | $-7.24\times10^{-2}$ | $-1.09\times10^{-4}$ |
> | 0.1 | $-0.242$ | $-1.89\times10^{-2}$ | $-6.99\times10^{-6}$ |
> | 0.05 | $-0.127$ | $-4.82\times10^{-3}$ | $-4.42\times10^{-7}$ |
> | ratio per halving | $\to2$ (order 1) | $3.9$ (order 2) | $15.8$ (order 4) |
>
> To reach an error of $10^{-6}$, Euler would need $h\approx4\times10^{-7}$ (about five million steps and evaluations); RK4 needs $h\approx0.06$ — about 33 steps, 130 evaluations.

### Step-size control

A fixed step is wasteful: transients need small steps, quiescent intervals large ones. Production integrators (MATLAB `ode45`, SciPy `solve_ivp` with `RK45`) use **embedded Runge–Kutta pairs**, such as Dormand–Prince 5(4): two RK formulas of orders $p$ and $p-1$ share the same stage evaluations, and their difference estimates the local error $\mathrm{err}$ of the lower-order one. The step is accepted if $\mathrm{err}\le\mathrm{tol}$ and the next step is chosen as

$$
h_{\text{new}} = h\cdot\min\!\left(f_{\max},\ \max\!\left(f_{\min},\ 0.9\left(\frac{\mathrm{tol}}{\mathrm{err}}\right)^{1/p}\right)\right),
$$

where the exponent follows from $\mathrm{err}\propto h^{p}$ and the safety factor $0.9$ and limits $f_{\min}, f_{\max}$ (e.g. $0.2$, $5$) prevent oscillation. The tolerance is usually mixed, $\mathrm{tol} = \text{atol} + \text{rtol}\cdot\lvert y\rvert$. Note that the controller bounds the *local* error per step; the global error is typically of the same order as the tolerance for well-behaved problems, but is not guaranteed by it.

### Multistep methods

Instead of sampling within the step, **linear multistep methods** reuse past values. The second-order **Adams–Bashforth** method,

$$
\mathbf y_{n+1} = \mathbf y_n + \frac h2\big(3\mathbf f_n - \mathbf f_{n-1}\big),
$$

integrates the linear interpolant of the last two slopes and costs only *one* new evaluation per step. Implicit **Adams–Moulton** methods and the **backward differentiation formulas (BDF)**, e.g. BDF2 $\tfrac32\mathbf y_{n+1} - 2\mathbf y_n + \tfrac12\mathbf y_{n-1} = h\mathbf f_{n+1}$, are the workhorses for stiff problems (below). Multistep methods need a starting procedure (a few RK steps) and are more awkward to use with variable steps.

### Absolute stability: the test equation

Convergence is a statement about $h\to0$. In practice $h$ is finite, and a method can produce wildly wrong, exponentially growing answers for a decaying problem if $h$ is too large. The analysis uses the scalar **test equation**

$$
y' = \lambda y,\qquad\lambda\in\mathbb C,\ \operatorname{Re}\lambda<0,
$$

whose solution decays. (It models each eigen-mode of a linear system $\mathbf y' = \mathbf A\mathbf y$, with $\lambda$ an eigenvalue of $\mathbf A$ — Ch. 10.) Applied to it, every RK method gives $y_{n+1} = R(z)\,y_n$ with $z = h\lambda$, where $R$ is the **stability function**:

| Method | $R(z)$ | Real stability interval |
|---|---|---|
| Euler (explicit) | $1 + z$ | $-2\le z\le0$ |
| Heun / midpoint (RK2) | $1 + z + \tfrac{z^2}{2}$ | $-2\le z\le0$ |
| RK4 | $1 + z + \tfrac{z^2}2 + \tfrac{z^3}6 + \tfrac{z^4}{24}$ | $-2.785\le z\le0$ |
| backward (implicit) Euler | $\dfrac{1}{1 - z}$ | all $z\le0$ |
| trapezoidal (Crank–Nicolson) | $\dfrac{1 + z/2}{1 - z/2}$ | all $z\le0$ |

> [!definition] Region of absolute stability; A-stability
> The **region of absolute stability** of a method is $\{z\in\mathbb C:\lvert R(z)\rvert\le1\}$. A method is **A-stable** if this region contains the entire left half-plane $\operatorname{Re}z\le0$, so that the numerical solution decays whenever the true one does, for every step size.

For explicit Euler the region is the disk $\lvert 1 + z\rvert\le1$ of radius 1 centred at $z = -1$. For real $\lambda<0$ the condition is $h\le2/\lvert\lambda\rvert$; for monotone (non-oscillating) decay, $h\le1/\lvert\lambda\rvert$. Explicit methods always have bounded stability regions (their $R$ is a polynomial, which is unbounded); the implicit Euler and trapezoidal rules are A-stable.

```desmos-graph
left=-0.1; right=3.2; bottom=-1.7; top=2.4
---
y=e^{-2x}|black
y=0.5^{4x}|blue|dashed
(0.25,0.5)|blue
(0.5,0.25)|blue
(0.75,0.125)|blue
(1,0.0625)|blue
y=1-2x|orange|0<x<0.75
y=-0.5+(x-0.75)|orange|0.75<x<1.5
y=0.25-0.5(x-1.5)|orange|1.5<x<2.25
y=-0.125+0.25(x-2.25)|orange|2.25<x<3
y=1-2x|red|0<x<1.25
y=-1.5+3(x-1.25)|red|1.25<x<2.5
y=2.25-4.5(x-2.5)|red|2.5<x<3.2
```

*Euler's method on $y'=-2y$, $y(0)=1$ (exact $e^{-2t}$, black). With $h=0.25$ ($z=-0.5$, blue) the numerical solution $0.5^{n}$ decays monotonically but too fast; with $h=0.75$ ($z=-1.5$, orange) it oscillates in sign while decaying ($\lvert R\rvert=0.5$); with $h=1.25$ ($z=-2.5$, red) $\lvert R\rvert=1.5>1$ and it grows without bound — outside the stability region.*

### Stiffness

> [!definition] Stiff problem (working definition)
> An ODE system is **stiff** when it contains modes with very different time scales — eigenvalues of the Jacobian $\partial\mathbf f/\partial\mathbf y$ with $\operatorname{Re}\lambda<0$ and $\max\lvert\operatorname{Re}\lambda\rvert\gg\min\lvert\operatorname{Re}\lambda\rvert$ — such that the step size of an explicit method is limited by **stability** for the fast, already-decayed modes rather than by **accuracy** for the slow mode of interest. The ratio $\max\lvert\operatorname{Re}\lambda\rvert/\min\lvert\operatorname{Re}\lambda\rvert$ is the stiffness ratio.

Stiffness is endemic in electronics: a circuit with a 1 ns parasitic time constant and a 1 ms thermal or bias time constant has stiffness ratio $10^6$; an explicit method would need a million steps per millisecond merely to remain stable, long after the fast transient has died.

**Implicit methods.** Backward Euler,

$$
\mathbf y_{n+1} = \mathbf y_n + h\,\mathbf f(t_{n+1},\mathbf y_{n+1}),
$$

requires solving a (generally nonlinear) equation for $\mathbf y_{n+1}$ at every step — by Newton's method (Section 22.4) with Jacobian $\mathbf I - h\,\partial\mathbf f/\partial\mathbf y$ — but is A-stable, so $h$ can be chosen for accuracy alone. Higher-order stiff solvers (BDF methods in `ode15s` and `solve_ivp(method='BDF')`, implicit RK methods such as Radau IIA) follow the same pattern. The **Dahlquist barrier** limits what is possible: no explicit linear multistep method is A-stable, and an A-stable linear multistep method has order at most 2 (the trapezoidal rule has the smallest error constant among them).

> [!example] Worked Example 22.36 — A stiff scalar problem
> $y' = -1000\,(y - \cos t)$, $y(0) = 0$. After a fast initial transient ($e^{-1000t}$, time constant 1 ms) the solution follows $y\approx\cos t$ closely. Exact value at $t = 0.1$: $y(0.1) = 0.9951030$. Here $\lambda = -1000$, so explicit Euler is stable only for $h\le2/1000 = 0.002$.
>
> | $h$ | steps to $t = 0.1$ | explicit Euler $y(0.1)$ | implicit Euler $y(0.1)$ |
> |---|---|---|---|
> | 0.001 | 100 | 0.9951035 | 0.9951025 |
> | 0.0021 | 48 | $-96.0$ | 0.9950226 |
> | 0.01 | 10 | $-3.5\times10^{9}$ | 0.9950980 |
>
> **Implicit Euler step.** Because $f$ is linear in $y$, the implicit equation $y_{n+1} = y_n - 1000h\,(y_{n+1} - \cos t_{n+1})$ is solved in closed form: $y_{n+1} = \dfrac{y_n + 1000h\cos t_{n+1}}{1 + 1000h}$.
>
> **Lesson.** Just $5\%$ beyond the stability limit ($h = 0.0021$, $z = -2.1$, $\lvert R\rvert = 1.1$) explicit Euler explodes; implicit Euler gives a 0.0005 % accurate answer with ten steps. The interesting solution ($\cos t$) varies on a time scale of seconds; only the stability of the method forced explicit Euler to millisecond steps.

> [!warning] Common Pitfall — Energy drift in oscillators
> For an undamped LC oscillator, $\ddot x + x = 0$, the eigenvalues are $\pm j$ — on the imaginary axis, which lies *outside* the stability region of explicit Euler except at $z = 0$. Euler multiplies the "energy" $x^2 + v^2$ by exactly $1 + h^2$ every step: with $h = 0.1$, after one period (62 steps) the amplitude has grown by a factor $\sqrt{1.853} = 1.36$. Backward Euler errs the other way, damping a lossless circuit. The trapezoidal rule ($\lvert R(j\omega h)\rvert = 1$ exactly) conserves amplitude, and *symplectic* methods (symplectic Euler, Störmer–Verlet/leapfrog) conserve a nearby energy for long times — the methods of choice for lossless systems such as orbits and particle accelerators.

> [!tip] Engineering Note — How SPICE integrates
> SPICE's transient analysis uses implicit methods throughout — the **trapezoidal rule** by default and **Gear (BDF) methods** of order 1–2 as alternatives — with Newton iteration at every time point and step control based on a local-truncation-error estimate. The trapezoidal rule's tendency to produce a spurious point-to-point ringing on stiff modes (since $R(z)\to-1$ as $z\to-\infty$: the fast mode is not damped, just sign-flipped) is the reason SPICE offers `method=gear`, whose $R(z)\to0$ ("L-stability") damps fast modes strongly. When a simulation shows ringing at the time-step frequency, switch the integration method before blaming the circuit.

## 22.10 Numerical solution of PDEs

### Grids and notation

The three model PDEs of Ch. 13 — the **heat (diffusion)** equation, the **wave** equation and the **Laplace/Poisson** equation — are prototypes of the parabolic, hyperbolic and elliptic classes, and each class demands a different numerical strategy. EE instances abound: diffusion of dopants and heat in a chip, the telegrapher's equations on a transmission line, Maxwell's equations in an FDTD simulator, the potential in an electrostatic lens or around a microstrip.

The **finite-difference method** replaces the continuous domain by a grid $x_i = i\,\Delta x$, $t_n = n\,\Delta t$, the unknown function by grid values $u_i^n\approx u(x_i,t_n)$, and every derivative by a difference quotient from Section 22.7. Throughout we use

$$
\delta_x^2u_i^n = u_{i+1}^n - 2u_i^n + u_{i-1}^n,\qquad\frac{\delta_x^2u_i^n}{\Delta x^2} = u_{xx}(x_i,t_n) + O(\Delta x^2).
$$

### The heat equation: explicit scheme (FTCS)

For $u_t = \alpha u_{xx}$ on $0<x<L$ with boundary values $u(0,t)$, $u(L,t)$ and initial data $u(x,0)$, use a **f**orward difference in **t**ime and a **c**entral difference in **s**pace:

$$
\frac{u_i^{n+1} - u_i^n}{\Delta t} = \alpha\,\frac{u_{i+1}^n - 2u_i^n + u_{i-1}^n}{\Delta x^2}
\quad\Longleftrightarrow\quad
u_i^{n+1} = u_i^n + r\left(u_{i+1}^n - 2u_i^n + u_{i-1}^n\right),\qquad r = \frac{\alpha\,\Delta t}{\Delta x^2}.
$$

The dimensionless **mesh ratio** $r$ controls everything. The scheme is explicit — each new value comes directly from three old ones — and its truncation error is $O(\Delta t) + O(\Delta x^2)$.

**Stability by Fourier modes.** Substitute a single grid Fourier mode $u_i^n = G^n e^{j\theta i}$ ($\theta = k\Delta x$) into the scheme:

$$
G = 1 + r\left(e^{j\theta} - 2 + e^{-j\theta}\right) = 1 - 2r(1-\cos\theta) = 1 - 4r\sin^2\!\frac\theta2 .
$$

The true solution mode decays, so we need $\lvert G\rvert\le1$ for all $\theta$. The worst case is the sawtooth mode $\theta = \pi$ ($u_i\propto(-1)^i$), where $G = 1 - 4r$, and $\lvert1 - 4r\rvert\le1$ requires

$$
\boxed{\ r = \frac{\alpha\,\Delta t}{\Delta x^2}\le\frac12\ }.
$$

This is a severe restriction: **halving $\Delta x$ forces $\Delta t$ to be quartered**, so refining the grid by 2 costs 8 times the work in one dimension. (Section 22.11 turns this calculation into the general von Neumann method.)

> [!example] Worked Example 22.37 — One FTCS step
> $u_t = u_{xx}$ on $[0,1]$, $u = 0$ at both ends, $u(x,0) = \sin\pi x$; exact solution $u = e^{-\pi^2t}\sin\pi x$. Take $\Delta x = 0.25$ and $r = 0.25$, so $\Delta t = r\Delta x^2 = 0.015625$.
>
> **Initial values.** $u^0 = (0,\ 0.707107,\ 1,\ 0.707107,\ 0)$.
>
> **Step.** $u_1^1 = 0.707107 + 0.25(1 - 2\cdot0.707107 + 0) = 0.603553$; $u_2^1 = 1 + 0.25(0.707107 - 2 + 0.707107) = 0.853553$; $u_3^1 = u_1^1$ by symmetry.
>
> **Compare.** Exact: $e^{-\pi^2\cdot0.015625} = 0.857090$, so $u(0.25,\Delta t) = 0.606054$, $u(0.5,\Delta t) = 0.857090$. The numerical mode is multiplied by $G = 1 - 4(0.25)\sin^2(\pi/8) = 0.853553$ per step versus the exact factor $0.857090$ — an error of 0.4 % per step from the coarse grid.

> [!example] Worked Example 22.38 — Watching FTCS go unstable
> Same problem with $\Delta x = 0.05$ (21 points), run to $t\approx0.1$, with a tiny perturbation $10^{-6}(-1)^i$ (think rounding noise) added to the initial data.
>
> | $r$ | steps | $\max_i\lvert u_i - u_{\text{exact}}\rvert$ |
> |---|---|---|
> | 0.50 | 80 | $1.5\times10^{-3}$ (truncation error only) |
> | 0.60 | 67 | $3.9\times10^{3}$ |
>
> At $r = 0.6$ the sawtooth mode has $G = 1 - 2.4 = -1.4$; after 67 steps it has grown by $1.4^{67}\approx 6\times10^{9}$, turning the $10^{-6}$ perturbation into thousands. Instability amplifies *whatever* high-frequency content is present — rounding error suffices — so an unstable scheme is useless however small the perturbation.

### Implicit schemes: BTCS and Crank–Nicolson

Evaluating the spatial difference at the *new* time level gives the **backward-time (BTCS) scheme**

$$
-r\,u_{i-1}^{n+1} + (1 + 2r)\,u_i^{n+1} - r\,u_{i+1}^{n+1} = u_i^n,
$$

a **tridiagonal** linear system per step (solved in $O(N)$ by the Thomas algorithm, as for splines). Its amplification factor $G = 1/(1 + 4r\sin^2(\theta/2))$ satisfies $0<G\le1$ for *every* $r$: the scheme is **unconditionally stable**, though still only first order in time.

Averaging the explicit and implicit schemes gives **Crank–Nicolson**:

$$
\frac{u_i^{n+1} - u_i^n}{\Delta t} = \frac{\alpha}{2}\left(\frac{\delta_x^2u_i^{n+1}}{\Delta x^2} + \frac{\delta_x^2u_i^n}{\Delta x^2}\right),\qquad G = \frac{1 - 2r\sin^2(\theta/2)}{1 + 2r\sin^2(\theta/2)},
$$

second order in both $\Delta t$ and $\Delta x$ and unconditionally stable ($\lvert G\rvert\le1$ always). It is the trapezoidal rule of Section 22.9 in disguise, and it inherits the trapezoidal rule's weakness: for large $r$ the high-frequency modes have $G\to-1$ (at $r = 5$, $\theta = \pi$: $G = -9/11 = -0.82$), so sharp initial data (a step in temperature, a switched voltage) produce slowly decaying point-to-point oscillations. A few initial BTCS steps ("Rannacher smoothing") suppress them.

### The method-of-lines view

Discretising *only in space* turns the PDE into a system of ODEs, $\dot{\mathbf u} = \frac{\alpha}{\Delta x^2}\mathbf T\mathbf u$, where $\mathbf T = \operatorname{tridiag}(1,-2,1)$. The eigenvalues of $\mathbf T$ (Ch. 10) are $-4\sin^2\!\big(\tfrac{k\pi}{2(N+1)}\big)$, so the ODE eigenvalues range from about $-\alpha\pi^2/L^2$ (the slow fundamental mode) to about $-4\alpha/\Delta x^2$ (the fast grid-scale mode). The system is therefore **stiff**, with stiffness ratio $\sim(2L/(\pi\Delta x))^2$, growing as the grid is refined. FTCS is explicit Euler on this system and its stability limit $\Delta t\cdot4\alpha/\Delta x^2\le2$ is exactly $r\le\tfrac12$; BTCS and Crank–Nicolson are implicit Euler and the trapezoidal rule. The method of lines thus lets any stiff ODE solver from Section 22.9 serve as a PDE solver — a common and powerful approach.

### The wave equation and the CFL condition

For $u_{tt} = c^2u_{xx}$, central differences in both variables give the explicit **leapfrog** scheme

$$
u_i^{n+1} = 2u_i^n - u_i^{n-1} + C^2\left(u_{i+1}^n - 2u_i^n + u_{i-1}^n\right),\qquad C = \frac{c\,\Delta t}{\Delta x},
$$

second order in both variables. $C$ is the **Courant number**. The first step needs $u^{-1}$; with initial velocity $u_t(x,0) = g(x)$, the central-difference condition $u_i^{1} - u_i^{-1} = 2\Delta t\,g_i$ gives $u_i^1 = u_i^0 + \Delta t\,g_i + \tfrac{C^2}{2}\delta_x^2u_i^0$.

> [!theorem] Courant–Friedrichs–Lewy (CFL) condition
> A necessary condition for the convergence of an explicit difference scheme for a hyperbolic equation is that the **numerical domain of dependence contain the physical domain of dependence**. For the leapfrog scheme for $u_{tt} = c^2u_{xx}$ (and for upwind schemes for $u_t + cu_x = 0$) this is
> $$
> C = \frac{c\,\Delta t}{\Delta x}\le1 ,
> $$
> and for leapfrog it is also sufficient for stability.

*Reasoning.* The exact solution at $(x_i,t_n)$ depends on initial data in $[x_i - ct_n,\ x_i + ct_n]$ (d'Alembert, Ch. 13). The scheme's value $u_i^n$ depends only on initial grid values in $[x_i - n\Delta x,\ x_i + n\Delta x]$ since information moves one cell per step. If $c\Delta t>\Delta x$ the numerical interval is smaller than the physical one, and changing the initial data in the gap changes the true solution but not the numerical one — so the scheme cannot converge. Von Neumann analysis confirms it: the amplification factors satisfy $G^2 - 2\big(1 - 2C^2\sin^2\tfrac\theta2\big)G + 1 = 0$, whose roots have $\lvert G\rvert = 1$ if $C\le1$, while for $C>1$ one root has $\lvert G\rvert>1$ at $\theta = \pi$. $\square$

> [!example] Worked Example 22.39 — The CFL limit in action
> A plucked string: $c = 1$, $u = \sin\pi x$ at $t = 0$, $u_t = 0$, fixed ends, $\Delta x = 0.05$, run to $t = 2$ (one full period); exact $u = \cos(\pi t)\sin(\pi x)$. A $10^{-6}$ sawtooth perturbation is added to the initial data.
>
> | Courant number $C$ | steps | max error at $t\approx2$ |
> |---|---|---|
> | 1.00 | 40 | $1.0\times10^{-6}$ (only the perturbation) |
> | 0.50 | 80 | $1.2\times10^{-5}$ |
> | 1.05 | 38 | $7.2\times10^{3}$ |
>
> **Why $C = 1$ is exact.** With $C = 1$ the scheme becomes $u_i^{n+1} = u_{i+1}^n + u_{i-1}^n - u_i^{n-1}$, which is satisfied exactly by any $F(x - ct) + G(x + ct)$ on the grid: the "magic time step" reproduces d'Alembert's solution. For $C<1$ the scheme is stable but **dispersive** (different wavelengths travel at slightly wrong, different speeds). For $C = 1.05$, only 5 % over the limit, the sawtooth mode grows explosively.

**First-order hyperbolic equations.** For the advection equation $u_t + a\,u_x = 0$ ($a>0$, a signal travelling along a lossless line), the natural explicit scheme takes the spatial difference **upwind** (from the side the information comes from):

$$
u_i^{n+1} = u_i^n - C\left(u_i^n - u_{i-1}^n\right),\qquad C = \frac{a\Delta t}{\Delta x}\le1 .
$$

It is stable for $0\le C\le1$, first order, and numerically diffusive (it smears sharp fronts). The tempting FTCS version with a central difference $u_{i+1}^n - u_{i-1}^n$ is **unconditionally unstable** ($\lvert G\rvert^2 = 1 + C^2\sin^2\theta>1$) — a scheme that is perfectly consistent but useless, illustrating that consistency is not enough (Section 22.11). Second-order schemes (Lax–Wendroff, leapfrog) reduce diffusion at the cost of dispersive wiggles near discontinuities.

> [!tip] Engineering Note — FDTD and the Courant limit
> The finite-difference time-domain (FDTD, Yee) method for Maxwell's equations is a staggered leapfrog scheme for $\mathbf E$ and $\mathbf H$. Its stability limit in 3-D on a cubic grid is $c\,\Delta t\le\Delta x/\sqrt3$ — the CFL condition with the 3-D domain of dependence. One consequence: halving the cell size to resolve a fine feature (a thin trace, a narrow slot) halves the allowable time step for the *whole* simulation, so run time scales like $(\Delta x)^{-4}$. This is why FDTD codes offer sub-gridding and why frequency-domain (FEM, MoM) solvers compete for electrically small fine structure.

### Laplace's and Poisson's equations

Elliptic problems have no time variable: the potential in a region is determined all at once by its boundary values. For Poisson's equation $\nabla^2V = -\rho/\varepsilon$ on a square grid with spacing $h$, the **five-point stencil**

$$
\frac{V_{i+1,j} + V_{i-1,j} + V_{i,j+1} + V_{i,j-1} - 4V_{i,j}}{h^2} = -\frac{\rho_{i,j}}{\varepsilon}
$$

has truncation error $O(h^2)$. For Laplace's equation ($\rho = 0$) it says each interior value is the **average of its four neighbours** — the discrete mean-value property, which also yields a discrete maximum principle (no interior extrema), mirroring the continuous theory of Ch. 13. Writing one equation per interior node gives a large, sparse, symmetric positive-definite linear system $\mathbf A\mathbf V = \mathbf b$, with boundary values moved into $\mathbf b$.

> [!example] Worked Example 22.40 — Potential in a square trough
> A square region has its top side at $100\ \mathrm V$ and the other three sides grounded. Use $h = \tfrac14$ of the side, giving a $3\times3$ array of interior nodes.
>
> **Equations.** For each interior node, $4V_{i,j} - (\text{sum of neighbours}) = (\text{sum of neighbouring boundary values})$. For example the top-left node: $4V_{1,1} - V_{1,2} - V_{2,1} = 100$; the centre: $4V_{2,2} - V_{1,2} - V_{3,2} - V_{2,1} - V_{2,3} = 0$.
>
> **Solution** (by symmetry 6 unknowns, solved directly):
>
> | row | left | centre | right |
> |---|---|---|---|
> | top ($y = \tfrac34$) | 42.857 | 52.679 | 42.857 |
> | middle ($y = \tfrac12$) | 18.750 | 25.000 | 18.750 |
> | bottom ($y = \tfrac14$) | 7.143 | 9.821 | 7.143 |
>
> **Check.** The centre value $25\ \mathrm V$ is exact even on this coarse grid: by superposition, the four problems with one side at $100\ \mathrm V$ sum to a box entirely at $100\ \mathrm V$, so each contributes $25\ \mathrm V$ at the centre — and the discrete problem has the same symmetry. The Fourier-series solution (Ch. 13) gives $25.000$ as well. Each value is also the average of its neighbours, e.g. $(100 + 18.75 + 52.679 + 0)/4 = 42.857$ for the top-left node.

**Iterative solution.** Real grids have $10^4$–$10^8$ nodes, and direct elimination fills in the sparse matrix. The classical iterations sweep the grid repeatedly:

- **Jacobi:** replace every $V_{i,j}$ by the average of its neighbours' *old* values.
- **Gauss–Seidel:** use new values as soon as they are available (sweeping in order); converges twice as fast as Jacobi for this problem.
- **Successive over-relaxation (SOR):** $V_{i,j}\leftarrow(1-\omega)V_{i,j} + \omega\cdot(\text{Gauss–Seidel value})$ with $1<\omega<2$; for the square with $N$ intervals per side the optimal factor is $\omega_{\text{opt}} = 2/(1 + \sin(\pi/N))$.

For the model problem the Jacobi iteration matrix has spectral radius $\cos(\pi h)\approx1 - \tfrac{\pi^2h^2}{2}$, so the number of iterations to reduce the error by a fixed factor grows like $h^{-2}$ (Gauss–Seidel: same order, half the count); optimal SOR's radius $\approx1 - 2\pi h$ makes it grow only like $h^{-1}$. In Example 22.40, reducing the error below $10^{-6}$ from a zero start takes 51 Jacobi, 25 Gauss–Seidel and 12 SOR sweeps ($\omega_{\text{opt}} = 2/(1 + \sin(\pi/4)) = 1.172$). Modern solvers use preconditioned conjugate gradients and multigrid, whose cost is nearly proportional to the number of unknowns (Ch. 24).

```python
import numpy as np
def laplace_sor(V, fixed, omega=1.8, tol=1e-8, maxit=10000):
    """V: 2-D array with boundary values set; fixed: boolean mask of Dirichlet nodes."""
    for it in range(maxit):
        delta = 0.0
        for i in range(1, V.shape[0] - 1):
            for j in range(1, V.shape[1] - 1):
                if fixed[i, j]:
                    continue
                gs = 0.25 * (V[i+1, j] + V[i-1, j] + V[i, j+1] + V[i, j-1])
                new = (1 - omega) * V[i, j] + omega * gs
                delta = max(delta, abs(new - V[i, j]))
                V[i, j] = new
        if delta < tol:
            return V, it + 1
    return V, maxit
```

**Boundary conditions.** Dirichlet values (fixed potentials, electrode voltages) enter $\mathbf b$ directly. A Neumann condition $\partial V/\partial n = g$ (a symmetry plane, an insulating boundary) is imposed with a **ghost point** outside the domain: at $x = 0$, the central difference $(V_{1,j} - V_{-1,j})/(2h) = g$ eliminates $V_{-1,j}$ from the stencil while keeping second-order accuracy. Curved boundaries and material interfaces (dielectric steps, where $\varepsilon\,\partial V/\partial n$ must be continuous) need modified stencils; this is where finite-element methods, which handle unstructured meshes naturally, usually take over.

> [!warning] Common Pitfall — Accuracy near corners and edges
> The $O(h^2)$ error estimate assumes the solution is smooth. Near re-entrant corners and electrode edges, fields are singular (e.g. $E\propto r^{-1/3}$ near a $90^\circ$ re-entrant corner, $r^{-1/2}$ near a thin conductor's edge) and convergence degrades sharply. Quantities computed from derivatives — fields, capacitance per unit length from surface charge — suffer most. Refine locally near singularities, or extract capacitance from a variational (energy) expression, which is less sensitive.

## 22.11 Stability and convergence

### Three properties of a discretisation

Every section of this chapter has met the same pattern: a discretisation that looks accurate (small truncation error) may still produce garbage (Examples 22.9, 22.36, 22.38, 22.39). This final section isolates the concepts that explain the pattern and states the theorems that make them precise. We phrase them for a linear time-dependent problem $u_t = \mathcal Lu$ (with $\mathcal L$ a spatial differential operator) discretised by a one-step scheme

$$
\mathbf u^{n+1} = \mathbf C(\Delta t)\,\mathbf u^{n},
$$

where $\mathbf u^n$ is the vector of grid values and $\mathbf C$ the **solution (amplification) operator** of the scheme, with the grid spacing $\Delta x$ tied to $\Delta t$ by a fixed rule (e.g. $r$ or $C$ held constant). The same ideas apply to ODE solvers ($\mathbf C = R(h\mathbf A)$ for $\mathbf y' = \mathbf A\mathbf y$) and to elliptic problems.

> [!definition] Consistency, stability, convergence
> - **Consistency.** The scheme is consistent if its **local truncation error** — the residual $\boldsymbol\tau^n = \big(\mathbf u_{\text{exact}}^{n+1} - \mathbf C\,\mathbf u_{\text{exact}}^n\big)/\Delta t$ obtained by inserting the exact solution — tends to zero as $\Delta t, \Delta x\to0$. It is consistent **of order $(p,q)$** if $\lVert\boldsymbol\tau\rVert = O(\Delta t^p + \Delta x^q)$.
> - **Stability.** The scheme is stable (in the sense of Lax–Richtmyer) if for every $T>0$ there is a constant $K_T$, independent of $\Delta t$ and $\Delta x$, such that
> $$
> \lVert\mathbf C(\Delta t)^n\rVert\le K_T\qquad\text{for all } n\text{ with } n\Delta t\le T.
> $$
> - **Convergence.** The scheme is convergent if, for every initial condition (in a suitable class), $\max_{n\Delta t\le T}\lVert\mathbf u^n - \mathbf u_{\text{exact}}(t_n)\rVert\to0$ as $\Delta t,\Delta x\to0$.

Consistency is a local, *one-step* property and is checked by Taylor expansion. Stability is a global property of *many steps*: errors introduced at any time (truncation, rounding, data) must not be amplified without bound. Convergence is what we actually want, but it is hard to verify directly, because it involves the unknown exact solution.

### The Lax equivalence theorem

> [!theorem] Lax equivalence theorem (Lax–Richtmyer, 1956)
> For a consistent finite-difference approximation to a **well-posed linear** initial-value problem,
> $$
> \text{stability}\iff\text{convergence}.
> $$
> Moreover, a stable scheme that is consistent of order $(p,q)$ converges with the same order: $\lVert\mathbf e^n\rVert = O(\Delta t^p + \Delta x^q)$.

*Proof of the important direction (stability + consistency ⇒ convergence).* Let $\mathbf e^n = \mathbf u_{\text{exact}}^n - \mathbf u^n$. By the definition of $\boldsymbol\tau$,

$$
\mathbf u_{\text{exact}}^{n+1} = \mathbf C\,\mathbf u_{\text{exact}}^{n} + \Delta t\,\boldsymbol\tau^n,\qquad\mathbf u^{n+1} = \mathbf C\,\mathbf u^{n}\quad\Longrightarrow\quad\mathbf e^{n+1} = \mathbf C\,\mathbf e^n + \Delta t\,\boldsymbol\tau^n .
$$

Unrolling with $\mathbf e^0 = \mathbf 0$:

$$
\mathbf e^{n} = \Delta t\sum_{k=0}^{n-1}\mathbf C^{\,n-1-k}\boldsymbol\tau^k\quad\Longrightarrow\quad\lVert\mathbf e^n\rVert\le\Delta t\sum_{k=0}^{n-1}K_T\lVert\boldsymbol\tau^k\rVert\le K_T\,T\max_k\lVert\boldsymbol\tau^k\rVert .
$$

Consistency makes the right side tend to zero, at the rate of the truncation error. $\square$

This is the same "local error × number of steps × bounded amplification" argument as for Euler's method in Section 22.9 — now in its general form. The converse (convergence ⇒ stability) uses the uniform boundedness principle of functional analysis.

The practical value of the theorem is enormous: **to prove that a linear scheme converges, check consistency (Taylor series) and stability (usually von Neumann analysis); you never need the exact solution.** It also explains the failures we saw — FTCS for the heat equation with $r>\tfrac12$, leapfrog with $C>1$, central-space FTCS for advection — all consistent, all unstable, therefore all non-convergent.

> [!warning] Common Pitfall — The scope of the Lax theorem
> The theorem is for **linear** problems. For nonlinear PDEs (shock-forming conservation laws, nonlinear diffusion) stability of the linearised scheme is necessary but not sufficient, and convergence may be to a wrong (non-physical) weak solution unless the scheme is also *conservative* and satisfies an *entropy* condition (Lax–Wendroff theorem). Also, "stable" is relative to a norm: a scheme can be stable in $\ell^2$ but produce overshoots in max-norm.

### Von Neumann (Fourier) stability analysis

For linear constant-coefficient schemes on uniform grids (with periodic boundaries, or away from boundaries), Fourier modes $u_i^n = G^n e^{j\theta i}$ diagonalise the solution operator, and $\lVert\mathbf C^n\rVert_2 = \max_\theta\lvert G(\theta)\rvert^n$ by Parseval's theorem (Ch. 15–16).

> [!theorem] Von Neumann stability condition
> A one-step, constant-coefficient scheme with amplification factor $G(\theta,\Delta t)$ is stable in $\ell^2$ if and only if
> $$
> \lvert G(\theta,\Delta t)\rvert\le1 + K\Delta t\qquad\text{for all }\theta\in[-\pi,\pi]
> $$
> and some constant $K$ independent of $\Delta t$, $\Delta x$. For problems whose exact solutions do not grow, the practical form is $\lvert G\rvert\le1$.

The allowance $1 + K\Delta t$ permits physically growing solutions (e.g. $u_t = u_{xx} + u$): then $\lvert G\rvert^n\le e^{KT}$, still bounded. For multi-level schemes (leapfrog), one analyses the matrix amplification factor or the characteristic equation for $G$, requiring its roots to satisfy $\lvert G\rvert\le1$ (and simple roots on the unit circle).

**Recipe.**

1. Substitute $u_i^n = G^ne^{j\theta i}$ into the homogeneous scheme.
2. Divide by $G^ne^{j\theta i}$ and solve for $G(\theta)$.
3. Find the conditions on the mesh ratios under which $\max_\theta\lvert G\rvert\le1$.

> [!example] Worked Example 22.41 — Von Neumann analysis of the upwind scheme
> Scheme: $u_i^{n+1} = u_i^n - C(u_i^n - u_{i-1}^n)$ for $u_t + au_x = 0$, $C = a\Delta t/\Delta x$.
>
> **Step 1–2.** $G = 1 - C\left(1 - e^{-j\theta}\right) = (1 - C) + Ce^{-j\theta}$.
>
> **Step 3.** $\lvert G\rvert^2 = (1 - C + C\cos\theta)^2 + C^2\sin^2\theta = 1 - 2C(1 - C)(1 - \cos\theta)$.
> For $0\le C\le1$, $C(1 - C)\ge0$ and $\lvert G\rvert\le1$ for all $\theta$: **stable**. For $C>1$ (or $C<0$, i.e. "downwind"), $C(1-C)<0$ and $\lvert G\rvert>1$ for $\theta\neq0$: **unstable**. The von Neumann condition reproduces the CFL condition.
>
> **Numerical check** ($C = 0.5$): $\lvert G\rvert^2 = 1,\ 0.854,\ 0.5,\ 0.146,\ 0$ at $\theta = 0,\ \pi/4,\ \pi/2,\ 3\pi/4,\ \pi$. The short-wavelength modes are strongly damped — the numerical diffusion of the upwind scheme.

> [!example] Worked Example 22.42 — Stability of the heat-equation schemes compared
> From Section 22.10, with $s = \sin^2(\theta/2)\in[0,1]$:
>
> | Scheme | $G(\theta)$ | Condition | Order |
> |---|---|---|---|
> | FTCS (explicit) | $1 - 4rs$ | $r\le\tfrac12$ | $(1,2)$ |
> | BTCS (implicit) | $\dfrac{1}{1 + 4rs}$ | none (unconditional) | $(1,2)$ |
> | Crank–Nicolson | $\dfrac{1 - 2rs}{1 + 2rs}$ | none (unconditional) | $(2,2)$ |
>
> The exact amplification of mode $\theta$ over one step is $e^{-\alpha k^2\Delta t} = e^{-r\theta^2}$; all three schemes match it for small $\theta$ (well-resolved modes), and differ only in how they treat grid-scale modes. Only FTCS can *amplify* them.

```desmos-graph
left=0; right=3.3; bottom=-1.6; top=1.3
---
y=1-4\cdot0.25\sin(x/2)^{2}|blue|0<x<3.1416
y=1-4\cdot0.5\sin(x/2)^{2}|green|0<x<3.1416
y=1-4\cdot0.6\sin(x/2)^{2}|red|0<x<3.1416
y=\frac{1-2\cdot5\sin(x/2)^{2}}{1+2\cdot5\sin(x/2)^{2}}|purple|dashed|0<x<3.1416
y=1|black|dotted
y=-1|black|dotted
```

*FTCS amplification factor $G(\theta)=1-4r\sin^2(\theta/2)$ versus $\theta\in[0,\pi]$ for $r=0.25$ (blue), $0.5$ (green, just touching $-1$ at $\theta=\pi$) and $0.6$ (red, dropping below $-1$: unstable), with Crank–Nicolson at $r=5$ (purple dashed), which stays inside $[-1,1]$ but approaches $-1$ for grid-scale modes.*

### Stability for ODE methods: zero-stability and absolute stability

For ODE solvers two different stability notions are needed, and it is important not to confuse them.

**Zero-stability** concerns $h\to0$ on a *fixed* interval and is the ODE counterpart of Lax–Richtmyer stability. For a linear multistep method $\sum_{k=0}^{s}\alpha_ky_{n+k} = h\sum_{k=0}^{s}\beta_kf_{n+k}$, it is governed by the **first characteristic polynomial** $\rho(\zeta) = \sum\alpha_k\zeta^k$.

> [!theorem] Dahlquist equivalence theorem
> A linear multistep method is convergent if and only if it is **consistent** ($\rho(1) = 0$, $\rho'(1) = \sum\beta_k$) and **zero-stable**: all roots of $\rho(\zeta)$ satisfy $\lvert\zeta\rvert\le1$, and roots on the unit circle are simple (the **root condition**). All one-step methods (Euler, RK) are automatically zero-stable.

> [!example] Worked Example 22.43 — A high-order method that does not converge
> The explicit two-step method
> $$
> y_{n+2} + 4y_{n+1} - 5y_n = h\,(4f_{n+1} + 2f_n)
> $$
> has the highest possible order (3) for its form — its truncation error is tiny. But $\rho(\zeta) = \zeta^2 + 4\zeta - 5 = (\zeta - 1)(\zeta + 5)$ has the root $-5$: the root condition fails.
>
> **Experiment** on $y' = -y$, $y(0) = 1$, with exact starting values, integrating to $t = 1$ (exact $e^{-1} = 0.368$):
>
> | $h$ | 0.1 | 0.05 | 0.025 |
> |---|---|---|---|
> | $y_N$ | $-6.68$ | $-4.65\times10^{6}$ | $-2.88\times10^{19}$ |
>
> The results get *worse* as $h$ decreases: the parasitic solution $(-5)^n$ is excited by rounding and truncation errors and grows faster when there are more steps. Consistency with high order, yet no convergence — the Dahlquist (and Lax) theorems in action.

**Absolute stability** (Section 22.9) concerns a *fixed* step $h>0$ applied to a decaying problem, and requires $h\lambda$ to lie in the stability region for every eigenvalue $\lambda$ of the Jacobian. A method can be zero-stable and convergent but useless at practical step sizes on stiff problems because its absolute stability region is small.

```desmos-graph
left=-3.5; right=2.5; bottom=-2; top=2
---
(x+1)^{2}+y^{2}\le1|blue
(x-1)^{2}+y^{2}\ge1|orange|dashed
x=0|black|dotted
```

*Absolute stability regions in the complex $z=h\lambda$ plane: explicit Euler, $\lvert1+z\rvert\le1$ (blue disk centred at $-1$); implicit Euler, $\lvert1-z\rvert\ge1$ (orange, everything outside the disk centred at $+1$), which contains the whole left half-plane — A-stability. The imaginary axis (dotted) is where undamped oscillators live, outside the explicit Euler region.*

### Dissipation, dispersion and the modified equation

Stability guarantees that errors stay bounded; it says nothing about their *character*. Writing $G = \lvert G\rvert e^{-j\omega_{\text{num}}\Delta t}$ for a wave-like problem and comparing with the exact factor $e^{-j\omega\Delta t}$:

- $\lvert G\rvert<1$ where the exact mode is undamped is **numerical dissipation** (amplitude error): sharp features are smeared;
- $\omega_{\text{num}}\neq\omega$ is **numerical dispersion** (phase error): different wavelengths travel at wrong and different speeds, producing trailing or leading wiggles.

The **modified equation** makes these effects quantitative: it is the PDE that the scheme solves *more accurately* than the original one, obtained by Taylor-expanding the scheme and eliminating time derivatives. For the upwind scheme,

$$
u_t + au_x = \frac{a\Delta x}{2}(1 - C)\,u_{xx} + O(\Delta x^2),
$$

so upwinding is, to leading order, the advection equation plus an artificial diffusion with coefficient $\tfrac{a\Delta x}{2}(1 - C)$ — positive (hence stable) exactly when $0\le C\le1$, and zero at $C = 1$, the exact "magic" step. For leapfrog and Lax–Wendroff the leading error term is $u_{xxx}$, a dispersive term, which explains their oscillations near fronts. In FDTD simulations numerical dispersion makes waves travel slightly slower than $c$, with an error depending on direction and on the number of cells per wavelength — the reason for the rule of thumb of at least 10–20 cells per wavelength.

### Verifying convergence in practice

The theorems tell us what *should* happen; verification checks that our code does it.

1. **Grid-refinement study.** Run on grids $h, h/2, h/4$; with a known exact solution, compute the observed order $p\approx\log_2(E_h/E_{h/2})$ (Section 22.2). Without one, use three grids: $p\approx\log_2\dfrac{u_h - u_{h/2}}{u_{h/2} - u_{h/4}}$.
2. **Method of manufactured solutions.** Pick any smooth $u_{\text{m}}(x,t)$, compute the source term $f = (u_{\text{m}})_t - \mathcal Lu_{\text{m}}$ symbolically, add $f$ to the code, and check that the numerical solution converges to $u_{\text{m}}$ at the design order. This tests every term, including boundary conditions.
3. **Richardson extrapolation** of the finest results to estimate the remaining discretisation error (as in Sections 22.7–22.8).
4. **Conservation checks.** Total charge, energy or mass that the PDE conserves should be conserved by the scheme to the expected accuracy.

> [!tip] Engineering Note — Reading a simulator's warnings
> "Time step too small", "timestep reduced below minimum", "Courant number exceeded", "solution diverging" — these messages are the stability theory of this section surfacing in commercial tools. An explicit electromagnetic or thermal solver chooses its time step from the CFL or $r\le\tfrac12$ condition applied to the *smallest* cell; an implicit solver chooses it for accuracy but may fail when its Newton iteration does not converge. Before trusting a simulation, perform at least one refinement study: halve the mesh size (and time step) and check that the quantity you care about changes by roughly the amount the method's order predicts.

## Summary

| Topic | Key formula or result | Notes |
|---|---|---|
| Floating point | $x = \pm(1.f)_2\times2^e$; $\operatorname{fl}(a\circ b) = (a\circ b)(1+\delta)$, $\lvert\delta\rvert\le u$ | double: $u = 2^{-53}\approx1.1\times10^{-16}$, $\varepsilon_M = 2^{-52}$ |
| Cancellation | relative error of $a-b$ amplified by $\lvert a\rvert/\lvert a-b\rvert$ | rewrite formulas (e.g. $1-\cos x = 2\sin^2\tfrac x2$) |
| Condition number | $\kappa_f = \lvert xf'(x)/f(x)\rvert$; $\kappa(\mathbf A) = \lVert\mathbf A\rVert\lVert\mathbf A^{-1}\rVert$ | lose $\approx\log_{10}\kappa$ digits |
| Stability (algorithm) | forward error $\lesssim\kappa\times$ backward error | backward stable ⇒ error $\approx\kappa u$ |
| Observed order | $p\approx\ln(E_1/E_2)/\ln(h_1/h_2)$ | the key debugging test |
| Bisection | $\lvert m_k - x^\ast\rvert\le(b-a)/2^{k+1}$ | guaranteed, 1 bit per step |
| Fixed point | $x_{k+1} = g(x_k)$; converges if $\lvert g'\rvert\le L<1$; $e_{k+1}\approx g'(x^\ast)e_k$ | linear unless $g'(x^\ast) = 0$ |
| Newton–Raphson | $x_{k+1} = x_k - f(x_k)/f'(x_k)$; $e_{k+1}\approx\frac{f''}{2f'}e_k^2$ | quadratic; local; linear at multiple roots |
| Newton (systems) | $\mathbf J(\mathbf x_k)\Delta\mathbf x = -\mathbf F(\mathbf x_k)$ | solve, never invert |
| Secant | $x_{k+1} = x_k - f_k\frac{x_k - x_{k-1}}{f_k - f_{k-1}}$ | order $\varphi = 1.618$ |
| Lagrange interpolation | $p_n = \sum y_i\ell_i$, $\ell_i = \prod_{j\neq i}\frac{x-x_j}{x_i-x_j}$ | unique, degree $\le n$ |
| Interpolation error | $f - p_n = \frac{f^{(n+1)}(\xi)}{(n+1)!}\prod(x - x_i)$ | Runge phenomenon for equispaced nodes |
| Chebyshev nodes | $x_i = \cos\frac{(2i+1)\pi}{2(n+1)}$; $\max\lvert\omega_{n+1}\rvert = 2^{-n}$ | near-optimal; $\Lambda_n\sim\frac2\pi\ln n$ |
| Cubic spline | tridiagonal system for $M_i = S''(x_i)$; error $\tfrac{5}{384}h^4\lVert f^{(4)}\rVert$ | minimum curvature |
| Least squares | $\mathbf A^\mathsf T\mathbf A\mathbf a = \mathbf A^\mathsf T\mathbf y$ (solve by QR) | orthogonal polynomials decouple coefficients |
| Minimax | error equioscillates at $n+2$ points | Remez; equiripple filters |
| Forward difference | $\frac{f(x+h)-f(x)}{h} = f' + \frac h2f''$ | $h_{\text{opt}}\sim\sqrt u$ |
| Central difference | $\frac{f(x+h)-f(x-h)}{2h} = f' + \frac{h^2}{6}f'''$ | $h_{\text{opt}}\sim u^{1/3}$ |
| Second difference | $\frac{f_1 - 2f_0 + f_{-1}}{h^2} = f'' + \frac{h^2}{12}f^{(4)}$ | |
| Richardson | $R = \frac{2^pD(h/2) - D(h)}{2^p - 1}$ | also gives error estimate |
| Composite trapezoid | $T_n = h[\tfrac12f_0 + f_1 + \cdots + \tfrac12f_n]$, $E = -\frac{(b-a)h^2}{12}f''$ | spectral for periodic $f$ |
| Composite Simpson | $S_n = \frac h3[f_0 + 4f_1 + 2f_2 + \cdots + f_n]$, $E = -\frac{(b-a)h^4}{180}f^{(4)}$ | exact for cubics |
| Gauss–Legendre | $n$ nodes = zeros of $P_n$; exact for degree $2n-1$ | 2-pt: $\pm1/\sqrt3$, weights 1 |
| Euler | $\mathbf y_{n+1} = \mathbf y_n + h\mathbf f(t_n,\mathbf y_n)$ | order 1; stable for $\lvert1 + h\lambda\rvert\le1$ |
| RK4 | $\mathbf y_{n+1} = \mathbf y_n + \frac h6(\mathbf k_1 + 2\mathbf k_2 + 2\mathbf k_3 + \mathbf k_4)$ | order 4; real stability limit $h\lambda\ge-2.785$ |
| Implicit Euler / trapezoid | $R(z) = \frac{1}{1-z}$, $\frac{1+z/2}{1-z/2}$ | A-stable; for stiff problems |
| Heat equation FTCS | $u_i^{n+1} = u_i^n + r\,\delta_x^2u_i^n$, $r = \frac{\alpha\Delta t}{\Delta x^2}\le\frac12$ | BTCS, Crank–Nicolson unconditionally stable |
| Wave equation leapfrog | $u_i^{n+1} = 2u_i^n - u_i^{n-1} + C^2\delta_x^2u_i^n$, $C = \frac{c\Delta t}{\Delta x}\le1$ | CFL condition |
| Laplace 5-point | $V_{i,j} = \tfrac14(V_{i+1,j} + V_{i-1,j} + V_{i,j+1} + V_{i,j-1})$ | Jacobi, Gauss–Seidel, SOR $\omega_{\text{opt}} = \frac{2}{1+\sin(\pi/N)}$ |
| Lax equivalence | consistent + stable $\iff$ convergent (linear, well-posed) | order of convergence = order of consistency |
| Von Neumann | $u_i^n = G^ne^{j\theta i}$; need $\lvert G(\theta)\rvert\le1 + K\Delta t$ | Fourier modes diagonalise the scheme |
| Dahlquist | consistent + root condition $\iff$ convergent (multistep) | A-stable LMM has order $\le2$ |

## Further reading

- R. L. Burden, J. D. Faires and A. M. Burden, *Numerical Analysis*, 10th ed., Cengage, 2016 — a thorough, example-driven introduction to every topic of this chapter.
- N. J. Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd ed., SIAM, 2002 — the definitive reference on floating-point arithmetic, conditioning and backward error analysis.
- L. N. Trefethen, *Approximation Theory and Approximation Practice*, extended ed., SIAM, 2019 — Chebyshev interpolation, the Runge phenomenon and Gaussian quadrature from a modern viewpoint.
- R. J. LeVeque, *Finite Difference Methods for Ordinary and Partial Differential Equations*, SIAM, 2007 — stability regions, stiffness, the CFL condition, von Neumann analysis and the Lax equivalence theorem.
- E. Hairer, S. P. Nørsett and G. Wanner, *Solving Ordinary Differential Equations I: Nonstiff Problems*, and E. Hairer and G. Wanner, *II: Stiff and Differential-Algebraic Problems*, Springer, 1993/1996 — the standard references on Runge–Kutta and multistep methods.
