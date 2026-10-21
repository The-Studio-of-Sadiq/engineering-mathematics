# Chapter 1 — Mathematical Language and Algebra

[← Table of Contents](../00-front/table-of-contents.md) · [Next →](02-complex-numbers.md)

> **Prerequisites:** Secondary-school arithmetic and algebra; read [Notation & Conventions](../00-front/notation-and-conventions.md) first.
> **Learning outcomes:** After this chapter you can
> - classify numbers, reason about precision and rounding, and use interval, set and absolute-value notation correctly;
> - rearrange, factor and simplify algebraic and rational expressions, including partial-fraction decomposition;
> - solve linear, quadratic, radical, exponential, logarithmic and absolute-value equations and inequalities, and recognise extraneous or lost roots;
> - work fluently with exponential, logarithmic, hyperbolic and trigonometric functions and their identities, including decibels and harmonic addition;
> - convert points and equations between Cartesian, polar, cylindrical $(r,\phi,z)$ and spherical $(R,\theta,\phi)$ coordinates;
> - state domains and ranges, build inverse and composite functions, read mathematical notation precisely, and check any formula for dimensional consistency.

Every later chapter in this book is written in the language set out here. Very little of it is new mathematics, but experience in design reviews, laboratories and examination halls is unanimous: engineering errors are far more often algebra slips, sign errors, misplaced factors of ten and unit mistakes than misunderstandings of calculus. A transfer function with one wrong sign predicts a stable amplifier that oscillates; a decibel figure computed with $10\log$ instead of $20\log$ is off by a factor of ten in power. This chapter is therefore deliberately thorough. It revisits the familiar material with the precision that later chapters assume — exact definitions, the reasons the rules hold, and the traps that catch even experienced engineers — and it points forward to where each tool is used. Treat it as the toolkit you will reach for on every page that follows.

## 1.1 Real and complex number systems

### Why the number system matters

An engineer measures a voltage as $3.30\ \text{V}$, a computer stores it as a 64-bit floating-point pattern, a theorem about limits assumes it is an exact real number, and an AC analysis represents its amplitude and phase as a single complex number. These are four different number systems, and moving between them without noticing is the source of many subtle mistakes. Knowing what each system can and cannot represent — and what operations are legal in each — is the first piece of mathematical language.

### The standard number sets

> [!definition] The standard number sets
> $$
> \mathbb{N}=\{1,2,3,\dots\}\subset\mathbb{Z}=\{\dots,-2,-1,0,1,2,\dots\}\subset\mathbb{Q}=\left\{\tfrac{p}{q}: p,q\in\mathbb{Z},\ q\neq 0\right\}\subset\mathbb{R}\subset\mathbb{C}
> $$
> - $\mathbb{N}$, the **natural numbers**, are the counting numbers.
> - $\mathbb{Z}$, the **integers**, add zero and the negatives, so subtraction is always possible.
> - $\mathbb{Q}$, the **rationals**, are ratios of integers, so division by any non-zero number is possible.
> - $\mathbb{R}$, the **reals**, fill the gaps in $\mathbb{Q}$; real numbers that are not rational are **irrational** (e.g. $\sqrt2$, $\pi$, $e$).
> - $\mathbb{C}$, the **complex numbers**, have the form $z=a+jb$ with $a,b\in\mathbb{R}$ and $j^2=-1$.

Some authors include $0$ in $\mathbb{N}$; this book does not, and writes $\mathbb{N}_0=\{0,1,2,\dots\}$ when zero is needed (for example the index $n$ of a discrete-time signal that starts at $n=0$). Each enlargement of the number system is motivated by an equation that could not be solved in the smaller system: $x+3=1$ has no solution in $\mathbb{N}$, $3x=1$ none in $\mathbb{Z}$, $x^2=2$ none in $\mathbb{Q}$, and $x^2=-1$ none in $\mathbb{R}$.

### Decimal expansions of rationals

> [!theorem] Decimal characterisation of $\mathbb{Q}$
> A real number is rational if and only if its decimal expansion either terminates or eventually repeats.

*Proof sketch.* Long division of $p$ by $q$ produces at each step a remainder in $\{0,1,\dots,q-1\}$. After at most $q$ steps a remainder must recur (there are only $q$ possibilities), and from then on the digits repeat. Conversely, a repeating decimal can be turned into a fraction by the shift-and-subtract trick in the next example. $\square$

> [!example] Worked Example 1.1 — Repeating decimals to fractions
> **(a)** $x=0.454545\ldots$ The repeating block has length 2, so multiply by $10^2$:
> $$
> 100x=45.4545\ldots,\qquad 100x-x=45\quad\Longrightarrow\quad x=\frac{45}{99}=\frac{5}{11}.
> $$
> **(b)** $y=1.2363636\ldots$ Here one digit ($2$) precedes the repeating block $36$. Shift so that the repeat starts at the decimal point, then shift by one block:
> $$
> 10y=12.363636\ldots,\qquad 1000y=1236.363636\ldots
> $$
> Subtracting, $990y=1224$, so $y=\dfrac{1224}{990}=\dfrac{68}{55}$.
>
> **Check:** $68\div55=1.236363\ldots$ ✓ and $5\div11=0.454545\ldots$ ✓

### Irrational numbers

**Why $\sqrt2$ is irrational.** Suppose $\sqrt2=p/q$ with $p,q$ integers having no common factor. Then $p^2=2q^2$, so $p^2$ is even, hence $p$ is even (the square of an odd number is odd), say $p=2k$. Then $4k^2=2q^2$, so $q^2=2k^2$ and $q$ is even too — contradicting "no common factor". $\square$

The same argument shows $\sqrt{n}$ is irrational whenever the integer $n$ is not a perfect square. The numbers $\pi$ and $e$ are also irrational (indeed **transcendental**: they are not roots of any polynomial with integer coefficients), though the proofs are harder. In practice this means that $\pi$, $e$, $\sqrt2$ and $\sqrt3$ — which appear constantly in RMS values, three-phase power and time constants — can only ever be *approximated* in a calculation. Keep them symbolic as long as possible and round only at the end.

### The real numbers as a complete ordered field

What exactly *is* $\mathbb{R}$? For engineering purposes it is enough to know the three groups of properties that every manipulation in this book relies on.

> [!definition] Axioms of the real numbers (summary)
> **Field axioms.** Addition and multiplication are commutative and associative; multiplication distributes over addition, $a(b+c)=ab+ac$; there are identities $0$ and $1$; every $a$ has a negative $-a$, and every $a\ne0$ has a reciprocal $a^{-1}$.
>
> **Order axioms.** For any $a,b$ exactly one of $a<b$, $a=b$, $a>b$ holds; if $a<b$ then $a+c<b+c$; if $a<b$ and $c>0$ then $ac<bc$.
>
> **Completeness.** Every non-empty set of reals that is bounded above has a **least upper bound** (supremum) in $\mathbb{R}$.

$\mathbb{Q}$ satisfies the field and order axioms but not completeness: the set $\{x\in\mathbb{Q}:x^2<2\}$ is bounded above but has no least upper bound *in* $\mathbb{Q}$ — the "missing" bound is $\sqrt2$. Completeness is equivalent to the statement that every bounded, increasing sequence converges to a real limit, and it is what makes limits, continuity, the intermediate value theorem and hence all of calculus work ([Ch. 4](../part-2-calculus/04-differential-calculus.md)). It also guarantees that numerical algorithms such as bisection converge to an actual number.

Every rule of algebra in Section 1.2 is a consequence of the field axioms. For example, $a\cdot0=0$ follows from $a\cdot0=a(0+0)=a\cdot0+a\cdot0$ and subtracting $a\cdot0$ from both sides; and $(-1)(-1)=1$ follows from $0=(-1)\cdot0=(-1)(1+(-1))=-1+(-1)(-1)$. Nothing is arbitrary.

### Intervals

Subsets of $\mathbb{R}$ that appear as domains, tolerance bands and solution sets are usually intervals.

| Notation | Set | Notation | Set |
|---|---|---|---|
| $[a,b]$ | $a\le x\le b$ (closed) | $(a,b)$ | $a<x<b$ (open) |
| $[a,b)$ | $a\le x<b$ (half-open) | $(a,b]$ | $a<x\le b$ |
| $(a,\infty)$ | $x>a$ | $[a,\infty)$ | $x\ge a$ |
| $(-\infty,b)$ | $x<b$ | $\mathbb{R}$ | $(-\infty,\infty)$ |

The symbol $\infty$ is never a real number and never sits next to a square bracket. Unions such as $(-\infty,-2)\cup[1,\infty)$ describe solution sets of inequalities (Section 1.3).

### Absolute value

> [!definition] Absolute value
> $$
> |x|=\begin{cases}x, & x\ge 0\\ -x, & x<0\end{cases}
> $$
> Equivalently $|x|=\sqrt{x^2}$. Geometrically, $|x-a|$ is the distance between $x$ and $a$ on the number line.

> [!theorem] Properties of absolute value
> For all real $x,y$:
> $$
> |x|\ge0,\qquad |xy|=|x|\,|y|,\qquad \left|\frac{x}{y}\right|=\frac{|x|}{|y|}\ (y\ne0),
> $$
> $$
> |x+y|\le |x|+|y|\ \ \text{(triangle inequality)},\qquad \big||x|-|y|\big|\le |x-y|\ \ \text{(reverse triangle inequality)}.
> $$

*Proof of the triangle inequality.* Since $-|x|\le x\le|x|$ and $-|y|\le y\le|y|$, adding gives $-(|x|+|y|)\le x+y\le |x|+|y|$, which is exactly $|x+y|\le|x|+|y|$. For the reverse form, write $x=(x-y)+y$ so $|x|\le|x-y|+|y|$, i.e. $|x|-|y|\le|x-y|$; swapping $x$ and $y$ gives $|y|-|x|\le|x-y|$, and together these are the claim. $\square$

The triangle inequality is the source of every **worst-case tolerance** bound in circuit design: if two errors are at most $\varepsilon_1$ and $\varepsilon_2$ in magnitude, their sum is at most $\varepsilon_1+\varepsilon_2$ in magnitude — and no better guarantee is possible without further information.

> [!example] Worked Example 1.2 — Worst-case tolerance
> A $1.00\ \text{k}\Omega$ resistor with $\pm5\%$ tolerance is placed in series with a $2.20\ \text{k}\Omega$ resistor with $\pm1\%$ tolerance. Bound the error in the total resistance.
>
> **Step 1 — individual error bounds.** $|e_1|\le0.05\times1000=50\ \Omega$ and $|e_2|\le0.01\times2200=22\ \Omega$.
>
> **Step 2 — triangle inequality.** The total is $R=3200+e_1+e_2$, and
> $$
> |e_1+e_2|\le|e_1|+|e_2|\le 72\ \Omega.
> $$
> **Step 3 — relative bound.** $72/3200=2.25\%$. So $R\in[3128,\,3272]\ \Omega$.
>
> **Check:** the extreme values $1050+2222=3272$ and $950+2178=3128$ are attained, so the bound is tight.

### Precision, significant figures and engineering notation

Measured quantities carry finite precision. The number of **significant figures** is the number of digits known reliably: $3.30\ \text{V}$ has three, $0.0047\ \mu\text{F}$ has two. **Scientific notation** writes a number as $m\times10^k$ with $1\le|m|<10$; **engineering notation** restricts $k$ to multiples of 3 so that the SI prefixes apply directly: $4.7\times10^{-9}\ \text{F}=4.7\ \text{nF}$, $2.2\times10^{4}\ \Omega=22\ \text{k}\Omega$. A result should not be quoted to more significant figures than the least precise input justifies, but intermediate calculations should carry two or three extra digits to avoid accumulating rounding error.

### Real numbers in a computer

Computers do not store real numbers; they store a finite set of **floating-point** numbers of the form $\pm m\times2^{e}$ with a fixed number of binary digits in $m$. IEEE 754 double precision uses a 53-bit significand, giving about 15–16 significant decimal digits; the gap between $1$ and the next representable number is the **machine epsilon** $\varepsilon_{\text{mach}}=2^{-52}\approx2.2\times10^{-16}$.

Most decimal fractions are not representable exactly in binary, just as $\tfrac13$ is not representable exactly in decimal. In binary, $0.1=0.000110011001100\ldots_2$ repeats forever, so the stored value is slightly off.

> [!example] Worked Example 1.3 — Floating-point surprises
> **(a)** In double precision, `0.1 + 0.2` evaluates to $0.30000000000000004$, so the test `0.1 + 0.2 == 0.3` is **false**. Every term was rounded to the nearest binary fraction and the rounding errors did not cancel. Comparisons of floating-point results must use a tolerance: $|a-b|\le\varepsilon\,\max(|a|,|b|)$.
>
> **(b)** Subtracting nearly equal numbers destroys relative accuracy. If $a=1.0000001$ and $b=1.0000000$ are each known to 8 significant figures, then $a-b=1\times10^{-7}$ is known to only **one** significant figure. This is **catastrophic cancellation**; Section 1.3 shows how it ruins the textbook quadratic formula, and [Ch. 22](../part-7-numerical-computational/22-numerical-methods.md) treats it systematically.

> [!tip] Engineering Note — fixed point and quantisation
> Microcontrollers and DSP chips often use **fixed-point** arithmetic: an $N$-bit integer interpreted as a fraction with a fixed binary point. An analogue-to-digital converter maps a continuous voltage range onto $2^N$ integer codes, so it replaces a real number by the nearest of a finite set — a quantisation error of at most half a step. An ideal 12-bit ADC spanning $3.3\ \text{V}$ has step $3.3/4096=0.806\ \text{mV}$ and therefore a worst-case quantisation error of $0.403\ \text{mV}$. Rational and real numbers are idealisations; the hardware works with integers.

### Complex numbers — a preview

The equation $x^2+1=0$ has no real solution, because $x^2\ge0$ for every real $x$. Adjoining a new number $j$ with $j^2=-1$ and insisting that the field axioms still hold produces $\mathbb{C}=\{a+jb:a,b\in\mathbb{R}\}$, with

$$
(a+jb)+(c+jd)=(a+c)+j(b+d),\qquad (a+jb)(c+jd)=(ac-bd)+j(ad+bc).
$$

The payoff is remarkable: in $\mathbb{C}$, **every** polynomial equation has a solution (the fundamental theorem of algebra, Section 1.4). Complex numbers are the natural language of AC circuits, filters and signals and are developed fully in [Ch. 2](02-complex-numbers.md). For now note only that $z=a+jb$ has real part $\operatorname{Re}\{z\}=a$ and imaginary part $\operatorname{Im}\{z\}=b$ (both real), that $z^*=a-jb$ is its **conjugate**, and that $|z|=\sqrt{a^2+b^2}$ is its **magnitude**, so that $zz^*=|z|^2$ is real and non-negative.

> [!warning] Common Pitfall — $\mathbb{C}$ is not ordered
> There is no way to define "$<$" on $\mathbb{C}$ that respects the order axioms: if $j>0$ then $j\cdot j=-1>0$, absurd; if $j<0$ then $-j>0$ and $(-j)(-j)=-1>0$, equally absurd. An inequality between complex numbers, such as $Z_1<Z_2$ for two impedances, is meaningless. Compare **magnitudes** $|Z_1|<|Z_2|$ or **real parts** instead, and say which.

> [!warning] Common Pitfall — the letter $j$
> Mathematics texts and many programming languages write $i$ for $\sqrt{-1}$; electrical engineering uses $j$ because $i$ denotes current. Python writes `1j`, MATLAB accepts both `1i` and `1j`. Never use $j$ simultaneously as a summation index and the imaginary unit.

## 1.2 Algebraic manipulation

### Why fluency matters

Engineering formulas are rarely used in the form in which they are derived. Ohm's law is written $V=IR$ but used as $R=V/I$; a resonance formula gives frequency in terms of $L$ and $C$ but the designer needs $C$ in terms of frequency; a transfer function must be put over a common denominator before its poles can be read off. Algebraic manipulation is the skill of transforming an expression or equation into an *equivalent* one that is more useful. Every transformation must be justified by the field axioms of Section 1.1 — and most mistakes come from "rules" that are not among them.

### Order of operations and the distributive law

Expressions are evaluated in the conventional order: brackets, then powers and roots, then multiplication and division (left to right), then addition and subtraction (left to right). Two conventions trip people up:

- $-x^2$ means $-(x^2)$, not $(-x)^2$. So $-3^2=-9$.
- $a/bc$ is ambiguous in print. In this book a fraction bar or explicit brackets are always used: $\dfrac{a}{bc}$ or $(a/b)c$.

The **distributive law** $a(b+c)=ab+ac$ is the only axiom that links addition and multiplication, so it is behind every expansion and every factorisation. Its two directions are **expanding** (left to right) and **factoring** (right to left). Repeated use gives $(a+b)(c+d)=ac+ad+bc+bd$ and, more generally, the product of two sums is the sum of all pairwise products of their terms.

### Laws of exponents

For a positive integer $n$, $a^n$ means $a\cdot a\cdots a$ ($n$ factors). The laws below are obvious for positive integer exponents by counting factors; the definitions $a^0=1$, $a^{-n}=1/a^n$ and $a^{p/q}=\sqrt[q]{a^p}$ are then *chosen* precisely so that the laws continue to hold.

> [!theorem] Laws of exponents
> For $a,b>0$ and real $m,n$:
> $$
> a^m a^n=a^{m+n},\qquad \frac{a^m}{a^n}=a^{m-n},\qquad (a^m)^n=a^{mn},\qquad (ab)^n=a^nb^n,\qquad \left(\frac{a}{b}\right)^n=\frac{a^n}{b^n},
> $$
> $$
> a^0=1,\qquad a^{-n}=\frac{1}{a^n},\qquad a^{1/n}=\sqrt[n]{a},\qquad a^{m/n}=\left(\sqrt[n]{a}\right)^m=\sqrt[n]{a^m}.
> $$

For example, requiring $a^{1/2}a^{1/2}=a^{1}$ forces $a^{1/2}=\sqrt a$, and requiring $a^0a^n=a^n$ forces $a^0=1$. Irrational exponents such as $2^{\sqrt2}$ are defined by continuity — as the limit of $2^{q}$ over rationals $q\to\sqrt2$ — or, more cleanly, as $e^{\sqrt2\ln2}$ (Section 1.5).

> [!warning] Common Pitfall — negative bases
> The laws above require $a,b>0$. With negative bases they fail: $\sqrt{(-2)^2}=2\ne-2$, and $\big((-8)^{2}\big)^{1/6}=2$ while $\big((-8)^{1/6}\big)^{2}$ is not even real. In general $\sqrt{x^2}=|x|$, not $x$. Calculators and programming languages also differ on whether $(-8)^{1/3}$ returns the real cube root $-2$ or a complex number.

### Expansion and factoring identities

The following identities should be recognisable instantly in both directions:

$$
\begin{aligned}
(a\pm b)^2 &= a^2\pm 2ab+b^2\\
(a+b)^3 &= a^3+3a^2b+3ab^2+b^3\\
a^2-b^2 &= (a-b)(a+b)\\
a^3\pm b^3 &= (a\pm b)(a^2\mp ab+b^2)\\
a^n-b^n &= (a-b)\left(a^{n-1}+a^{n-2}b+\dots+ab^{n-2}+b^{n-1}\right)
\end{aligned}
$$

The last identity, verified by multiplying out (all middle terms cancel in pairs), contains the geometric-series formula of Section 1.10 as the special case $a=1$, $b=r$.

### The binomial theorem

> [!theorem] Binomial theorem
> For any positive integer $n$,
> $$
> (a+b)^n=\sum_{k=0}^{n}\binom{n}{k}a^{n-k}b^k,\qquad \binom{n}{k}=\frac{n!}{k!\,(n-k)!},
> $$
> where $n!=1\cdot2\cdots n$ and $0!=1$.

*Proof sketch.* Expanding $(a+b)(a+b)\cdots(a+b)$ ($n$ factors) means choosing $a$ or $b$ from each factor and multiplying. A term $a^{n-k}b^k$ arises once for each way of choosing which $k$ of the $n$ factors supply $b$, and there are $\binom{n}{k}$ such choices. $\square$

The coefficients form **Pascal's triangle**, in which each entry is the sum of the two above it, $\binom{n}{k}=\binom{n-1}{k-1}+\binom{n-1}{k}$:

| $n$ | coefficients $\binom{n}{0},\dots,\binom{n}{n}$ |
|---|---|
| 0 | 1 |
| 1 | 1, 1 |
| 2 | 1, 2, 1 |
| 3 | 1, 3, 3, 1 |
| 4 | 1, 4, 6, 4, 1 |
| 5 | 1, 5, 10, 10, 5, 1 |

The binomial coefficients reappear as counts of outcomes in the binomial probability distribution ([Ch. 19](../part-6-probability-statistics/19-probability.md)).

> [!example] Worked Example 1.4 — A specific binomial term
> Find the term independent of $x$ in $\left(2x-\dfrac1x\right)^6$.
>
> **Step 1 — general term.** With $a=2x$, $b=-1/x$:
> $$
> T_k=\binom{6}{k}(2x)^{6-k}\left(-\frac1x\right)^k=\binom{6}{k}2^{6-k}(-1)^k\,x^{6-2k}.
> $$
> **Step 2 — choose $k$.** The power of $x$ is zero when $6-2k=0$, i.e. $k=3$.
>
> **Step 3 — evaluate.** $T_3=\binom63\,2^3(-1)^3=20\cdot8\cdot(-1)=-160$.
>
> **Check:** substituting $x=1$ in the full expansion gives $(2-1)^6=1$; the seven coefficients $64,-192,240,-160,60,-12,1$ do sum to $1$. ✓

### The binomial approximation

For small $x$ the binomial expansion is dominated by its first two terms:

$$
(1+x)^n\approx 1+nx,\qquad |x|\ll1.
$$

Remarkably, this holds for **any real** $n$ — negative and fractional included — as the Taylor series of [Ch. 4](../part-2-calculus/04-differential-calculus.md) shows. The next term is $\tfrac{n(n-1)}{2}x^2$, which estimates the error. Special cases used constantly:

$$
\frac{1}{1+x}\approx1-x,\qquad \sqrt{1+x}\approx1+\frac x2,\qquad \frac{1}{\sqrt{1+x}}\approx1-\frac x2.
$$

> [!example] Worked Example 1.5 — Loading of a voltage divider
> A divider $R_1=R_2=10\ \text{k}\Omega$ feeds a load $R_L=100\ \text{k}\Omega$. Its output with the load attached can be written
> $$
> V_{\text{out}}=V_{\text{in}}\frac{R_2}{R_1+R_2}\cdot\frac{1}{1+R_{\text{th}}/R_L},\qquad R_{\text{th}}=\frac{R_1R_2}{R_1+R_2}=5\ \text{k}\Omega.
> $$
> **Step 1.** $x=R_{\text{th}}/R_L=0.05$ is small, so $\dfrac{1}{1+x}\approx1-x=0.95$.
>
> **Step 2.** The loaded output is about $5\%$ below the unloaded value $V_{\text{in}}/2$.
>
> **Check:** the exact factor is $1/1.05=0.95238$; the approximation error $0.0024$ matches the size of the neglected term $x^2=0.0025$. ✓ The approximation turns a messy expression into a rule of thumb: *loading error $\approx R_{\text{th}}/R_L$*.

### Completing the square

> [!theorem] Completing the square
> For $a\ne0$,
> $$
> ax^2+bx+c=a\left(x+\frac{b}{2a}\right)^2+c-\frac{b^2}{4a}.
> $$

*Derivation.* Factor out $a$: $a\left(x^2+\frac bax\right)+c$. Since $\left(x+\frac{b}{2a}\right)^2=x^2+\frac bax+\frac{b^2}{4a^2}$, adding and subtracting $\frac{b^2}{4a^2}$ inside the bracket gives the result. $\square$

This single manipulation yields the quadratic formula (Section 1.3), locates the vertex $\left(-\frac{b}{2a},\,c-\frac{b^2}{4a}\right)$ of a parabola, converts the equation of a circle to centre–radius form (Section 1.7), and in [Ch. 14](../part-5-transforms/14-laplace-transforms.md) turns a second-order denominator $s^2+2\zeta\omega_ns+\omega_n^2$ into $(s+\zeta\omega_n)^2+\omega_d^2$, the form that identifies a damped sinusoid.

> [!example] Worked Example 1.6 — Completing the square
> **(a)** $2x^2-8x+3=2(x^2-4x)+3=2\big[(x-2)^2-4\big]+3=2(x-2)^2-5$. Minimum value $-5$ at $x=2$.
>
> **(b)** $s^2+6s+13=(s+3)^2-9+13=(s+3)^2+2^2$. As a sum of squares it never vanishes for real $s$; its complex zeros are $s=-3\pm j2$.
>
> **Check (a):** expanding, $2(x^2-4x+4)-5=2x^2-8x+3$. ✓ **Check (b):** $(-3+j2)^2+6(-3+j2)+13=(5-j12)+(-18+j12)+13=0$. ✓

### Algebraic fractions

The rules for fractions follow from the field axioms:

$$
\frac ab+\frac cd=\frac{ad+bc}{bd},\qquad \frac ab\cdot\frac cd=\frac{ac}{bd},\qquad \frac{a/b}{c/d}=\frac{ad}{bc},\qquad \frac{ac}{bc}=\frac ab\ (c\ne0).
$$

Three working principles:

1. **Combine over the lowest common denominator**, factoring denominators first so you do not create needlessly large products.
2. **Cancel only common factors**, never common terms. In $\dfrac{x(x+1)}{x(x+2)}$ the factor $x$ cancels (for $x\ne0$); in $\dfrac{x+1}{x+2}$ nothing cancels.
3. **Rationalise** surd denominators by multiplying by the conjugate, using $(\sqrt a-\sqrt b)(\sqrt a+\sqrt b)=a-b$:

$$
\frac{1}{\sqrt{a}+\sqrt{b}}=\frac{\sqrt a-\sqrt b}{a-b}.
$$

The same conjugate trick, applied to $a+jb$, is how complex division is done in [Ch. 2](02-complex-numbers.md).

> [!example] Worked Example 1.7 — Simplifying fractions
> **(a) Compound fraction.** Simplify $\dfrac{\frac1x-\frac1y}{\frac1{x^2}-\frac1{y^2}}$.
> Numerator $=\dfrac{y-x}{xy}$; denominator $=\dfrac{y^2-x^2}{x^2y^2}=\dfrac{(y-x)(y+x)}{x^2y^2}$. Dividing,
> $$
> \frac{y-x}{xy}\cdot\frac{x^2y^2}{(y-x)(y+x)}=\frac{xy}{x+y}\qquad(x\ne y).
> $$
> **Check** with $x=2$, $y=3$: original $=\dfrac{1/2-1/3}{1/4-1/9}=\dfrac{1/6}{5/36}=1.2$; result $=6/5=1.2$. ✓ (The result is the "parallel combination" of $x$ and $y$.)
>
> **(b) Rationalising.** $\dfrac{1}{\sqrt5-\sqrt3}=\dfrac{\sqrt5+\sqrt3}{5-3}=\dfrac{\sqrt5+\sqrt3}{2}\approx1.98406$. **Check:** $1/(2.23607-1.73205)=1/0.50402=1.98406$. ✓

> [!example] Worked Example 1.8 — Parallel resistors
> Two resistors in parallel satisfy $\dfrac{1}{R}=\dfrac{1}{R_1}+\dfrac{1}{R_2}$. (a) Express $R$ explicitly. (b) What $R_2$ must be paralleled with $R_1=10\ \text{k}\Omega$ to give $R=4\ \text{k}\Omega$?
>
> **(a)** Common denominator: $\dfrac1R=\dfrac{R_2+R_1}{R_1R_2}$, so the familiar "product over sum"
> $$
> R=\frac{R_1R_2}{R_1+R_2}.
> $$
> **(b)** Solve $\dfrac1{R_2}=\dfrac1R-\dfrac1{R_1}=\dfrac{R_1-R}{RR_1}$, giving
> $$
> R_2=\frac{RR_1}{R_1-R}=\frac{4\times10}{10-4}\ \text{k}\Omega=6.667\ \text{k}\Omega.
> $$
> **Check:** $\dfrac{10\times6.667}{16.667}=4.000\ \text{k}\Omega$. ✓ Note that the formula requires $R<R_1$ — a parallel combination is always smaller than either branch.

### Changing the subject of a formula

To make a variable the **subject**, apply the same operation to both sides repeatedly, undoing the operations that act on the variable in reverse order. Keep track of conditions (division by zero, square roots of negatives, sign of quantities) as you go.

> [!example] Worked Example 1.9 — Two rearrangements
> **(a) Resonance.** An $LC$ tank resonates at $f_0=\dfrac{1}{2\pi\sqrt{LC}}$. Find $C$ for $f_0=1\ \text{MHz}$ with $L=10\ \mu\text{H}$.
> Square both sides: $f_0^2=\dfrac{1}{4\pi^2LC}$, then solve: $C=\dfrac{1}{4\pi^2f_0^2L}$. Numerically
> $$
> C=\frac{1}{4\pi^2\,(10^6)^2\,(10\times10^{-6})}=2.533\times10^{-9}\ \text{F}=2.533\ \text{nF}.
> $$
> **(b) Feedback.** A negative-feedback amplifier has closed-loop gain $A_f=\dfrac{A}{1+A\beta}$. Find $\beta$ for $A=10^5$ and $A_f=10$.
> Multiply out: $A_f+A_fA\beta=A$, so $\beta=\dfrac{A-A_f}{AA_f}=\dfrac1{A_f}-\dfrac1A=0.1-0.00001=0.09999$.
>
> **Checks:** (a) $1/(2\pi\sqrt{10^{-5}\times2.533\times10^{-9}})=1.000\times10^6\ \text{Hz}$ ✓; (b) $10^5/(1+10^5\times0.09999)=10^5/10000=10$ ✓. The form $\beta=1/A_f-1/A$ also shows why $\beta\approx1/A_f$ when $A$ is huge.

> [!warning] Common Pitfall — illegal "simplifications"
> $\dfrac{a+b}{a+c}\neq\dfrac{b}{c}$, $\ \sqrt{a^2+b^2}\neq a+b$, $\ (a+b)^{-1}\neq a^{-1}+b^{-1}$, $\ (a+b)^2\ne a^2+b^2$. Each of these cancels or distributes an operation across a **sum**, which is never allowed. A two-second numerical spot check (put $a=1$, $b=2$, $c=3$) exposes all of them.

> [!tip] Engineering Note — keep it symbolic
> Rearrange symbolically first and substitute numbers last. The symbolic answer can be dimension-checked (Section 1.10), tested at limiting cases ($R_2\to\infty$ in a parallel combination should give $R\to R_1$), and reused for other component values. Plugging numbers in early hides all of these checks.

## 1.3 Equations and inequalities

### Equivalent equations

Solving an equation means finding every value of the unknown that makes it true — the **solution set**. The method is to replace the equation by a chain of simpler **equivalent** equations (equations with the same solution set). The safe operations are:

- adding the same quantity to both sides;
- multiplying or dividing both sides by the same **non-zero** quantity;
- applying a **one-to-one** function to both sides (e.g. $\ln$ to positive quantities, $\exp$, cubing).

Two operations are *not* equivalence-preserving and must be handled with care:

- **Squaring** both sides (or applying any function that is not one-to-one) can introduce **extraneous roots**: $x=2$ has one solution, but $x^2=4$ has two.
- **Dividing by an expression that may be zero** can **lose roots**: dividing $x^2=3x$ by $x$ gives $x=3$ and loses $x=0$. Factor instead: $x(x-3)=0$.

Whenever a non-equivalent step is used, every candidate must be checked in the *original* equation.

### Linear equations

A linear equation $ax+b=0$ with $a\ne0$ has the single root $x=-b/a$. If $a=0$ there is no solution ($b\ne0$) or every $x$ is a solution ($b=0$). This trivial trichotomy — one, none, or infinitely many — reappears unchanged for systems of linear equations in [Ch. 8](../part-3-linear-algebra/08-matrices-and-linear-systems.md).

### Quadratic equations

> [!theorem] Quadratic formula
> The roots of $ax^2+bx+c=0$, $a\ne0$, are
> $$
> x=\frac{-b\pm\sqrt{\Delta}}{2a},\qquad \Delta=b^2-4ac.
> $$
> - $\Delta>0$: two distinct real roots;
> - $\Delta=0$: one repeated real root $x=-b/(2a)$;
> - $\Delta<0$: a complex-conjugate pair $x=\dfrac{-b}{2a}\pm j\dfrac{\sqrt{-\Delta}}{2a}$.

*Derivation.* By completing the square (Section 1.2), $ax^2+bx+c=0$ becomes $a\left(x+\frac b{2a}\right)^2=\frac{b^2-4ac}{4a}$, i.e. $\left(x+\frac b{2a}\right)^2=\frac{\Delta}{4a^2}$. Taking square roots, $x+\frac{b}{2a}=\pm\frac{\sqrt\Delta}{2a}$. $\square$

> [!theorem] Vieta's relations for a quadratic
> If $x_1,x_2$ are the roots of $ax^2+bx+c=0$, then
> $$
> x_1+x_2=-\frac ba,\qquad x_1x_2=\frac ca.
> $$

These follow from comparing coefficients in $a(x-x_1)(x-x_2)=ax^2-a(x_1+x_2)x+ax_1x_2$. They are the cheapest check on any computed pair of roots.

The sign of the discriminant decides whether a second-order circuit is over-, critically or under-damped ([Ch. 12](../part-4-differential-equations/12-linear-differential-equations.md)) — the same three cases, because the characteristic equation $LCs^2+RCs+1=0$ of a series RLC circuit is a quadratic.

> [!example] Worked Example 1.10 — RLC characteristic roots
> A series RLC circuit with $R=200\ \Omega$, $L=0.1\ \text{H}$, $C=2\ \mu\text{F}$ has characteristic equation $s^2+\frac RLs+\frac1{LC}=0$, i.e. $s^2+2000s+5\times10^{6}=0$. Solve it.
>
> **Step 1.** $\Delta=2000^2-4\cdot5\times10^6=4\times10^6-2\times10^7=-1.6\times10^{7}<0$: complex roots, an under-damped response.
>
> **Step 2.** $\sqrt{-\Delta}=4000$, so
> $$
> s=\frac{-2000\pm j4000}{2}=-1000\pm j2000\ \text{s}^{-1}.
> $$
> **Check (Vieta):** sum $=-2000=-b/a$ ✓; product $=1000^2+2000^2=5\times10^6=c/a$ ✓. The real part $-1000$ is the decay rate and the imaginary part $2000\ \text{rad/s}$ the ringing frequency.

### Numerically stable roots

The textbook formula is exact in exact arithmetic but can be disastrous in floating point. When $b^2\gg|4ac|$, $\sqrt\Delta\approx|b|$, and one of $-b\pm\sqrt\Delta$ subtracts two nearly equal numbers — catastrophic cancellation (Section 1.1).

> [!theorem] Stable quadratic formula
> For real coefficients with $\Delta\ge0$, compute
> $$
> q=-\tfrac12\left(b+\operatorname{sgn}(b)\sqrt{\Delta}\right),\qquad x_1=\frac{q}{a},\qquad x_2=\frac{c}{q}.
> $$
> The sum inside $q$ never cancels, and $x_2$ follows from Vieta's product $x_1x_2=c/a$.

> [!example] Worked Example 1.11 — Cancellation in the quadratic formula
> Solve $x^2-10^8x+1=0$ in double precision.
>
> **Naive formula.** $\sqrt{\Delta}=\sqrt{10^{16}-4}$ rounds to exactly $10^8$ in double precision. The small root $\big(10^8-\sqrt\Delta\big)/2$ is computed as $7.45\times10^{-9}$ — wrong by $25\%$ (the true value is $1.0000\times10^{-8}$).
>
> **Stable formula.** $q=-\tfrac12\left(-10^8-10^8\right)=10^8$, so $x_1=10^8$ and $x_2=1/10^8=10^{-8}$, correct to full precision.
>
> **Check:** product $x_1x_2=1=c/a$ ✓; sum $\approx10^8=-b/a$ ✓.

### Equations reducible to quadratics; radical equations

Many equations become quadratic after a substitution: $x^4-5x^2+4=0$ with $u=x^2$ gives $u^2-5u+4=0$, so $u=1$ or $u=4$ and $x\in\{-2,-1,1,2\}$. Equations in $e^x$ (Section 1.5) and $\sin x$ (Section 1.6) yield to the same trick.

Equations containing square roots are solved by isolating the root and squaring — a non-equivalent step, so every candidate must be checked.

> [!example] Worked Example 1.12 — An extraneous root
> Solve $\sqrt{x+5}=x-1$.
>
> **Step 1 — square.** $x+5=x^2-2x+1$, i.e. $x^2-3x-4=0$, so $(x-4)(x+1)=0$ and the candidates are $x=4$ and $x=-1$.
>
> **Step 2 — check in the original.** $x=4$: $\sqrt9=3=4-1$ ✓. $x=-1$: $\sqrt4=2$ but $-1-1=-2$ ✗.
>
> **Solution:** $x=4$ only. The root $x=-1$ solves $\sqrt{x+5}=-(x-1)$, which squares to the same quadratic. Note the original equation also requires $x-1\ge0$, which rules out $x=-1$ in advance.

### Simultaneous equations

Two linear equations in two unknowns are solved by **substitution** (solve one for a variable, substitute in the other) or **elimination** (add multiples of the equations to remove a variable). Eliminating in general gives **Cramer's rule** for the $2\times2$ case:

$$
\begin{aligned}
a_1x+b_1y&=c_1\\
a_2x+b_2y&=c_2
\end{aligned}
\qquad\Longrightarrow\qquad
x=\frac{c_1b_2-c_2b_1}{a_1b_2-a_2b_1},\quad y=\frac{a_1c_2-a_2c_1}{a_1b_2-a_2b_1},
$$

provided the **determinant** $D=a_1b_2-a_2b_1\ne0$. Geometrically each equation is a line; $D\ne0$ means the lines cross at one point, while $D=0$ means they are parallel (no solution) or coincident (infinitely many). Larger systems, and the matrix language that organises them, are the subject of [Ch. 8](../part-3-linear-algebra/08-matrices-and-linear-systems.md).

> [!example] Worked Example 1.13 — Two-mesh circuit
> Mesh analysis of a two-loop resistive circuit ($R_1=2\ \Omega$, $R_2=4\ \Omega$ shared, $R_3=6\ \Omega$, source $10\ \text{V}$ in loop 1) gives
> $$
> \begin{aligned}
> 6I_1-4I_2&=10\\
> -4I_1+10I_2&=0
> \end{aligned}
> $$
> **Step 1 — determinant.** $D=6\cdot10-(-4)(-4)=60-16=44$.
>
> **Step 2 — Cramer.** $I_1=\dfrac{10\cdot10-0\cdot(-4)}{44}=\dfrac{100}{44}=2.273\ \text{A}$, $\quad I_2=\dfrac{6\cdot0-(-4)\cdot10}{44}=\dfrac{40}{44}=0.909\ \text{A}$.
>
> **Check:** $6(2.273)-4(0.909)=13.636-3.636=10.000$ ✓ and $-4(2.273)+10(0.909)=0.000$ ✓.

Nonlinear simultaneous equations are usually reduced by substitution to a single equation in one unknown.

> [!example] Worked Example 1.14 — Line meets circle
> Find the intersections of $y=x+1$ and $x^2+y^2=25$.
>
> Substitute: $x^2+(x+1)^2=25\Rightarrow 2x^2+2x-24=0\Rightarrow x^2+x-12=0\Rightarrow(x+4)(x-3)=0$.
> So $x=3,\ y=4$ or $x=-4,\ y=-3$: the points $(3,4)$ and $(-4,-3)$.
>
> **Check:** $9+16=25$ ✓, $16+9=25$ ✓. A negative discriminant here would have meant the line misses the circle.

### Inequalities

Inequalities obey the same rules as equations **except**:

1. multiplying or dividing by a **negative** number reverses the inequality ($-2<3$ but $4>-6$);
2. taking reciprocals of two quantities of the same sign reverses it ($0<a<b\Rightarrow 1/a>1/b$);
3. squaring preserves an inequality only between **non-negative** quantities ($0\le a<b\Rightarrow a^2<b^2$, but $-3<2$ while $9>4$).

Each rule is a consequence of the order axioms in Section 1.1. For instance, if $a<b$ and $c<0$, then $-c>0$ so $-ac<-bc$, i.e. $ac>bc$.

| Inequality | Solution set |
|---|---|
| $\lvert x-a\rvert<\delta$ | $a-\delta<x<a+\delta$ |
| $\lvert x-a\rvert>\delta$ | $x<a-\delta$ or $x>a+\delta$ |
| $(x-r_1)(x-r_2)<0,\ r_1<r_2$ | $r_1<x<r_2$ (between the roots) |
| $(x-r_1)(x-r_2)>0,\ r_1<r_2$ | $x<r_1$ or $x>r_2$ (outside the roots) |

For polynomial or rational inequalities use a **sign chart**: bring everything to one side as a single fraction compared with $0$, find all zeros of numerator and denominator, split the real line at them, and determine the sign on each piece (each linear factor changes sign at its own zero; factors with even multiplicity do not).

> [!example] Worked Example 1.15 — Three inequalities
> **(a)** $|x-2|<3 \iff -3<x-2<3 \iff -1<x<5$, i.e. $x\in(-1,5)$.
>
> **(b)** $x^2-5x+6\le 0 \iff (x-2)(x-3)\le0 \iff x\in[2,3]$.
>
> **(c)** $\dfrac{x-1}{x+2}\ge0$. Critical points $x=1$ (zero of numerator, included) and $x=-2$ (zero of denominator, always excluded).
>
> | Interval | $x-1$ | $x+2$ | quotient |
> |---|---|---|---|
> | $x<-2$ | $-$ | $-$ | $+$ |
> | $-2<x<1$ | $-$ | $+$ | $-$ |
> | $x>1$ | $+$ | $+$ | $+$ |
>
> Solution: $x\in(-\infty,-2)\cup[1,\infty)$. We may **not** multiply through by $x+2$ without knowing its sign — doing so blindly gives $x\ge1$ and loses the whole interval $x<-2$.

```desmos-graph
left=-8; right=6; bottom=-6; top=6
---
y=\frac{x-1}{x+2}|blue
x=-2|red|dashed
y=0|grey
y=1|grey|dashed
(1,0)|label:zero at x=1
```
*$y=(x-1)/(x+2)$: the curve is on or above the axis for $x<-2$ and $x\ge1$, matching the sign chart. The red dashed line is the vertical asymptote where the denominator vanishes.*

> [!example] Worked Example 1.16 — Two absolute values
> Solve $|x-1|+|x+2|\le5$.
>
> The expressions inside change sign at $x=1$ and $x=-2$, so treat three cases.
>
> - $x\ge1$: $(x-1)+(x+2)=2x+1\le5\Rightarrow x\le2$. Combined: $1\le x\le2$.
> - $-2\le x<1$: $-(x-1)+(x+2)=3\le5$, always true. Combined: $-2\le x<1$.
> - $x<-2$: $-(x-1)-(x+2)=-2x-1\le5\Rightarrow x\ge-3$. Combined: $-3\le x<-2$.
>
> **Solution:** $x\in[-3,2]$. **Check:** at $x=2$, $1+4=5$ ✓; at $x=-3$, $4+1=5$ ✓; at $x=0$, $1+2=3\le5$ ✓. Geometrically, the sum of distances from $x$ to $1$ and to $-2$ is at most $5$.

### The AM–GM inequality and a design bound

> [!theorem] Arithmetic–geometric mean inequality
> For $a,b\ge0$,
> $$
> \frac{a+b}{2}\ge\sqrt{ab},
> $$
> with equality if and only if $a=b$.

*Proof.* $(\sqrt a-\sqrt b)^2\ge0$ expands to $a-2\sqrt{ab}+b\ge0$. Equality holds only when $\sqrt a=\sqrt b$. $\square$

Inequalities like this one prove that a design is optimal *without calculus*.

> [!example] Worked Example 1.17 — Maximum power transfer
> A source of EMF $V$ and internal resistance $R_s$ drives a load $R_L$, delivering $P=\dfrac{V^2R_L}{(R_s+R_L)^2}$. Show $P\le\dfrac{V^2}{4R_s}$ and find when equality holds.
>
> **Step 1.** By AM–GM with $a=R_s$, $b=R_L$: $R_s+R_L\ge2\sqrt{R_sR_L}$, so $(R_s+R_L)^2\ge4R_sR_L$.
>
> **Step 2.** Therefore $P=\dfrac{V^2R_L}{(R_s+R_L)^2}\le\dfrac{V^2R_L}{4R_sR_L}=\dfrac{V^2}{4R_s}$, with equality iff $R_L=R_s$.
>
> **Numerical check:** $V=10\ \text{V}$, $R_s=50\ \Omega$: the bound is $100/200=0.5\ \text{W}$, attained at $R_L=50\ \Omega$ ($100\cdot50/100^2=0.5$ ✓). At $R_L=100\ \Omega$: $100\cdot100/150^2=0.444\ \text{W}<0.5$ ✓.

### Inequalities as specifications

> [!example] Worked Example 1.18 — Choosing a component range
> An RC timer needs $\tau=RC\in[0.85,\,1.15]\ \text{ms}$ for **every** capacitor in a batch of nominal $100\ \text{nF}\pm10\%$. Which resistances guarantee this?
>
> **Step 1.** $C\in[90,110]\ \text{nF}$ and $R>0$, so $\tau$ ranges over $[90R,\,110R]$ (nF times $R$).
>
> **Step 2.** Both ends must satisfy the spec: $90\times10^{-9}R\ge0.85\times10^{-3}$ and $110\times10^{-9}R\le1.15\times10^{-3}$.
>
> **Step 3.** $R\ge9.444\ \text{k}\Omega$ and $R\le10.455\ \text{k}\Omega$. So $R\in[9.44,\,10.45]\ \text{k}\Omega$, and a $10.0\ \text{k}\Omega\ \pm1\%$ resistor (range $9.9$–$10.1\ \text{k}\Omega$) qualifies.
>
> **Check:** worst cases $9.9\ \text{k}\Omega\times90\ \text{nF}=0.891\ \text{ms}\ge0.85$ ✓ and $10.1\ \text{k}\Omega\times110\ \text{nF}=1.111\ \text{ms}\le1.15$ ✓.

> [!tip] Engineering Note
> Design specifications are inequalities: "gain within $\pm0.5$ dB", "ripple below $50\ \text{mV}$", "phase margin at least $45^\circ$". Solving them gives the admissible range of component values, and the worst-case analysis above (check every combination of extremes) is the systematic way to guarantee them. When there are many parameters, the same problem becomes constrained optimisation ([Ch. 23](../part-7-numerical-computational/23-optimization.md)).

> [!warning] Common Pitfall
> Never multiply both sides of an inequality by an expression whose sign you do not know, and never square both sides unless both are known to be non-negative. When in doubt, move everything to one side and use a sign chart.

## 1.4 Polynomials and rational expressions

### Why polynomials and rational functions

Polynomials are the simplest functions that can be computed with nothing but addition and multiplication, which is why they are used to approximate everything else (Taylor series, [Ch. 4](../part-2-calculus/04-differential-calculus.md); curve fitting, [Ch. 22](../part-7-numerical-computational/22-numerical-methods.md)). Ratios of polynomials — **rational functions** — are even more central to electrical engineering: the impedance of any network of resistors, capacitors and inductors, and the transfer function of any linear time-invariant circuit, is a rational function of the complex frequency $s$. Its **zeros** and **poles** (roots of numerator and denominator) determine the circuit's frequency response and stability. Factoring polynomials and decomposing rational functions are therefore everyday skills.

### Definitions and arithmetic

> [!definition] Polynomial
> A **polynomial** of **degree** $n$ in $x$ is
> $$
> p(x)=a_nx^n+a_{n-1}x^{n-1}+\dots+a_1x+a_0,\qquad a_n\ne0.
> $$
> The $a_k$ are the **coefficients**, $a_n$ the **leading coefficient** (the polynomial is **monic** if $a_n=1$), and $a_0$ the **constant term**. A **root** or **zero** is a number $r$ with $p(r)=0$.

Sums and products of polynomials are polynomials, with $\deg(pq)=\deg p+\deg q$. Division is like integer division: it leaves a quotient and a remainder.

> [!theorem] Division algorithm
> For polynomials $p$ and $d$ with $d\ne0$ there are unique polynomials $q$ (quotient) and $r$ (remainder) with
> $$
> p(x)=d(x)\,q(x)+r(x),\qquad \deg r<\deg d\ \ (\text{or } r=0).
> $$

The quotient and remainder are found by **long division**: repeatedly divide the leading term of what remains by the leading term of $d$, multiply back, and subtract. When $d(x)=x-r$ the bookkeeping shrinks to **synthetic division**, which keeps only the coefficients.

> [!example] Worked Example 1.19 — Long and synthetic division
> **(a)** Divide $p(x)=x^3+2$ by $d(x)=x^2-1$.
> Leading terms: $x^3/x^2=x$. Subtract $x(x^2-1)=x^3-x$ from $x^3+0x^2+0x+2$ to leave $x+2$, whose degree is below 2. So
> $$
> x^3+2=(x^2-1)\,x+(x+2),\qquad q(x)=x,\quad r(x)=x+2.
> $$
> **(b)** Divide $2x^3-3x^2-11x+6$ by $x-3$ synthetically. Write the coefficients, bring down the leading one, and repeatedly multiply by $3$ and add:
>
> | | $2$ | $-3$ | $-11$ | $6$ |
> |---|---|---|---|---|
> | $\times3$ | | $6$ | $9$ | $-6$ |
> | sum | $2$ | $3$ | $-2$ | $\mathbf{0}$ |
>
> Quotient $2x^2+3x-2$, remainder $0$.
>
> **Check:** $(x-3)(2x^2+3x-2)=2x^3+3x^2-2x-6x^2-9x+6=2x^3-3x^2-11x+6$ ✓.

### Remainder and factor theorems

> [!theorem] Remainder and factor theorems
> - **Remainder theorem:** dividing $p(x)$ by $(x-r)$ leaves the constant remainder $p(r)$.
> - **Factor theorem:** $(x-r)$ is a factor of $p(x)$ if and only if $p(r)=0$.

*Proof.* By the division algorithm $p(x)=(x-r)q(x)+c$ with $c$ constant (degree below 1). Setting $x=r$ gives $p(r)=c$. The factor theorem is the case $c=0$. $\square$

So the remainder of $x^5-3x+2$ on division by $x+1$ is $p(-1)=-1+3+2=4$, found without dividing.

### The fundamental theorem of algebra

> [!theorem] Fundamental theorem of algebra
> Every polynomial of degree $n\ge1$ with complex (in particular, real) coefficients has exactly $n$ complex roots counted with multiplicity, and factors completely as
> $$
> p(x)=a_n(x-r_1)(x-r_2)\cdots(x-r_n).
> $$
> If all coefficients are **real**, non-real roots occur in **conjugate pairs** $\alpha\pm j\beta$.

The existence of at least one root is the deep part (it needs completeness of $\mathbb{R}$ and is proved with complex analysis in [Ch. 18](../part-5-transforms/18-complex-analysis.md)); the full factorisation then follows by applying the factor theorem repeatedly. The conjugate-pair property has a one-line proof: conjugating $\sum a_kr^k=0$ with real $a_k$ gives $\sum a_k(r^*)^k=0$, so $r^*$ is a root too.

Consequently every **real** polynomial factors into real linear factors and real **irreducible quadratics** $x^2+bx+c$ with $b^2<4c$, since $(x-\alpha-j\beta)(x-\alpha+j\beta)=(x-\alpha)^2+\beta^2$. A real polynomial of odd degree therefore always has at least one real root.

A root $r$ has **multiplicity** $m$ if $(x-r)^m$ divides $p$ but $(x-r)^{m+1}$ does not. On a graph, the curve crosses the axis at a root of odd multiplicity and touches it (turning back) at a root of even multiplicity; higher multiplicity makes the curve flatter at the root.

> [!theorem] Vieta's relations
> For $p(x)=a_nx^n+\dots+a_0$ with roots $r_1,\dots,r_n$:
> $$
> \sum_i r_i=-\frac{a_{n-1}}{a_n},\qquad \sum_{i<k}r_ir_k=\frac{a_{n-2}}{a_n},\qquad \prod_i r_i=(-1)^n\frac{a_0}{a_n}.
> $$

These follow from expanding $a_n\prod(x-r_i)$ and comparing coefficients. In a control system, the sum of the closed-loop poles is fixed by one coefficient, which is why moving one pole to the left forces another to the right.

### Finding roots: the rational-root test

> [!theorem] Rational-root test
> If $p$ has integer coefficients and a rational root $u/v$ in lowest terms, then $u$ divides $a_0$ and $v$ divides $a_n$.

*Proof.* Multiply $p(u/v)=0$ by $v^n$: $a_nu^n+a_{n-1}u^{n-1}v+\dots+a_0v^n=0$. Every term except $a_0v^n$ contains $u$, so $u$ divides $a_0v^n$; as $u,v$ share no factor, $u$ divides $a_0$. Similarly $v$ divides $a_nu^n$, hence $a_n$. $\square$

> [!example] Worked Example 1.20 — Factoring a cubic
> Factor $p(x)=2x^3-3x^2-11x+6$ completely.
>
> **Step 1 — candidates.** $u\mid6$ and $v\mid2$: $\pm1,\pm2,\pm3,\pm6,\pm\tfrac12,\pm\tfrac32$.
>
> **Step 2 — test.** $p(1)=-6$, $p(-1)=12$, $p(3)=54-27-33+6=0$. ✓
>
> **Step 3 — deflate.** By the synthetic division of Example 1.19(b), $p(x)=(x-3)(2x^2+3x-2)$, and $2x^2+3x-2=(2x-1)(x+2)$.
>
> **Result:** $p(x)=(x-3)(2x-1)(x+2)$, roots $3,\ \tfrac12,\ -2$.
>
> **Check (Vieta):** sum $3+\tfrac12-2=\tfrac32=-(-3)/2$ ✓; product $3\cdot\tfrac12\cdot(-2)=-3=(-1)^3\cdot6/2$ ✓.

> [!example] Worked Example 1.21 — Characteristic polynomial with complex roots
> A third-order system has characteristic polynomial $D(s)=s^3+4s^2+6s+4$. Find its roots.
>
> **Step 1.** Rational candidates $\pm1,\pm2,\pm4$. All coefficients are positive, so no positive real root is possible; try negatives: $D(-1)=-1+4-6+4=1$, $D(-2)=-8+16-12+4=0$. ✓
>
> **Step 2.** Deflate: $D(s)=(s+2)(s^2+2s+2)$.
>
> **Step 3.** Complete the square: $s^2+2s+2=(s+1)^2+1$, so $s=-1\pm j$.
>
> **Roots:** $-2,\ -1+j,\ -1-j$ — all with negative real part, so the system is stable ([Ch. 14](../part-5-transforms/14-laplace-transforms.md)). **Check:** sum $=-4=-a_2/a_3$ ✓; product $(-2)\big((-1)^2+1^2\big)=-4=(-1)^3\cdot4$ ✓.

For degree 5 and above there is no general formula in radicals (Abel–Ruffini theorem), and even the cubic and quartic formulas are rarely worth using. In practice, high-degree roots are found numerically — Newton's method or the eigenvalues of a companion matrix ([Ch. 22](../part-7-numerical-computational/22-numerical-methods.md)).

```desmos-graph
left=-3.5; right=4; bottom=-25; top=20
---
y=2x^{3}-3x^{2}-11x+6|blue
y=(x-1)^{2}(x+2)|red|dashed
(3,0)|label:x=3
(0.5,0)|label:x=1/2
(-2,0)|label:x=-2
```
*The cubic of Example 1.20 (blue) crosses the axis at its three simple roots. The red dashed cubic $(x-1)^2(x+2)$ touches the axis at its double root $x=1$ without crossing.*

### Rational functions

> [!definition] Rational function
> A **rational function** is $R(x)=\dfrac{N(x)}{D(x)}$ with polynomials $N,D$ and $D\not\equiv0$. It is **proper** if $\deg N<\deg D$ and **improper** otherwise. After cancelling common factors, zeros of $N$ are the **zeros** of $R$ and zeros of $D$ are its **poles**.

An improper rational function is reduced by long division to a polynomial plus a proper part: $R(x)=q(x)+\dfrac{r(x)}{D(x)}$. This also identifies the **asymptotes** of its graph:

- **vertical asymptotes** at the real poles, where $|R(x)|\to\infty$;
- a **horizontal asymptote** $y=a_n/b_n$ if $\deg N=\deg D$, or $y=0$ if $\deg N<\deg D$;
- an **oblique** (slant) asymptote $y=q(x)$ if $\deg N=\deg D+1$.

> [!example] Worked Example 1.22 — Asymptotes
> $R(x)=\dfrac{x^2+1}{x-1}$. Division gives $x^2+1=(x-1)(x+1)+2$, so
> $$
> R(x)=x+1+\frac{2}{x-1}.
> $$
> Vertical asymptote $x=1$ (pole); oblique asymptote $y=x+1$, approached from above as $x\to+\infty$ (since $\frac2{x-1}>0$) and from below as $x\to-\infty$. No real zeros, as $x^2+1>0$.
>
> **Check:** at $x=3$, $R=10/2=5$ and $3+1+2/2=5$ ✓.

```desmos-graph
left=-6; right=8; bottom=-8; top=12
---
y=\frac{x^{2}+1}{x-1}|blue
y=x+1|red|dashed
x=1|grey|dashed
```
*$R(x)=(x^2+1)/(x-1)$ (blue) hugs its oblique asymptote $y=x+1$ (red) far from the pole at $x=1$.*

> [!warning] Common Pitfall — cancelled factors
> $\dfrac{x^2-1}{x-1}=x+1$ only for $x\ne1$. The left side is undefined at $x=1$ (a "hole", not a pole). In a transfer function, cancelling a pole against a zero is mathematically legal but can hide an internal mode that still exists in the hardware — an unstable cancelled pole is still unstable.

### Partial-fraction decomposition

Adding fractions over a common denominator is easy; partial fractions run the process **backwards**, splitting a complicated proper rational function into a sum of simple terms, each of which can be integrated, inverse-transformed or interpreted on its own.

> [!theorem] Partial-fraction decomposition
> Every proper real rational function $N(x)/D(x)$ can be written uniquely as a sum of terms, one group for each factor of $D$:
>
> | Factor in $D(x)$ | Terms in the decomposition |
> |---|---|
> | $(x-a)$ | $\dfrac{A}{x-a}$ |
> | $(x-a)^m$ | $\dfrac{A_1}{x-a}+\dfrac{A_2}{(x-a)^2}+\dots+\dfrac{A_m}{(x-a)^m}$ |
> | $x^2+bx+c$ (irreducible) | $\dfrac{Bx+C}{x^2+bx+c}$ |
> | $(x^2+bx+c)^m$ | $\displaystyle\sum_{k=1}^{m}\frac{B_kx+C_k}{(x^2+bx+c)^k}$ |

The number of unknown constants always equals $\deg D$, which is why the decomposition is unique: multiplying out gives $\deg D$ linear equations in $\deg D$ unknowns. There are three ways to find them, often combined:

1. **Cover-up rule** (distinct linear factor). Multiply both sides by $(x-a)$ and set $x=a$; every other term vanishes:
$$
A=\Big[(x-a)R(x)\Big]_{x=a}=\left.\frac{N(x)}{D(x)/(x-a)}\right|_{x=a}.
$$
2. **Substituting convenient values** of $x$ into the cleared-of-fractions identity.
3. **Equating coefficients** of like powers of $x$.

For a repeated factor $(x-a)^m$, the cover-up rule gives the top coefficient $A_m=\big[(x-a)^mR(x)\big]_{x=a}$, and the others follow by differentiation:

$$
A_{m-k}=\frac{1}{k!}\,\frac{d^k}{dx^k}\Big[(x-a)^mR(x)\Big]_{x=a},\qquad k=1,\dots,m-1.
$$

(Differentiation is introduced in [Ch. 4](../part-2-calculus/04-differential-calculus.md); until then, method 3 always works.)

> [!example] Worked Example 1.23 — Distinct and repeated linear factors
> **(a)** $\dfrac{s+3}{(s+1)(s+2)}=\dfrac{A}{s+1}+\dfrac{B}{s+2}$.
> Cover-up: $A=\dfrac{-1+3}{-1+2}=2$, $B=\dfrac{-2+3}{-2+1}=-1$.
> **Check:** $2(s+2)-(s+1)=s+3$ ✓.
>
> **(b)** $\dfrac{1}{s(s+1)^2}=\dfrac{A}{s}+\dfrac{B}{s+1}+\dfrac{C}{(s+1)^2}$.
> Clear fractions: $1=A(s+1)^2+Bs(s+1)+Cs$.
> $s=0$: $A=1$. $\ s=-1$: $C=-1$. Coefficient of $s^2$: $0=A+B\Rightarrow B=-1$.
> (By the derivative formula: $B=\frac{d}{ds}\big[1/s\big]_{s=-1}=-1/s^2\big|_{-1}=-1$ ✓.)
> **Check** at $s=1$: LHS $=\tfrac14$; RHS $=1-\tfrac12-\tfrac14=\tfrac14$ ✓.

> [!example] Worked Example 1.24 — An irreducible quadratic factor
> Decompose $\dfrac{s+2}{(s+1)(s^2+2s+5)}$.
>
> **Step 1.** $s^2+2s+5=(s+1)^2+4$ has $\Delta=4-20<0$, so it is irreducible. Write
> $$
> \frac{s+2}{(s+1)(s^2+2s+5)}=\frac{A}{s+1}+\frac{Bs+C}{s^2+2s+5}.
> $$
> **Step 2 — cover-up for $A$.** $A=\dfrac{-1+2}{1-2+5}=\dfrac14$.
>
> **Step 3 — equate coefficients** in $s+2=\tfrac14(s^2+2s+5)+(Bs+C)(s+1)$:
> $s^2$: $0=\tfrac14+B\Rightarrow B=-\tfrac14$. Constant: $2=\tfrac54+C\Rightarrow C=\tfrac34$.
>
> **Step 4 — check the unused coefficient** ($s^1$): $\tfrac24+B+C=\tfrac12-\tfrac14+\tfrac34=1$ ✓.
>
> **Result**, with the quadratic term arranged for later use in Laplace inversion:
> $$
> \frac{1/4}{s+1}+\frac{-\tfrac14 s+\tfrac34}{(s+1)^2+4}=\frac{1/4}{s+1}-\frac14\cdot\frac{s+1}{(s+1)^2+2^2}+\frac12\cdot\frac{2}{(s+1)^2+2^2}.
> $$

> [!example] Worked Example 1.25 — An improper fraction
> Decompose $\dfrac{x^3+2}{x^2-1}$.
>
> **Step 1 — divide first** (Example 1.19a): $\dfrac{x^3+2}{x^2-1}=x+\dfrac{x+2}{(x-1)(x+1)}$.
>
> **Step 2 — cover-up** on the proper part: at $x=1$, $\dfrac{1+2}{1+1}=\dfrac32$; at $x=-1$, $\dfrac{-1+2}{-1-1}=-\dfrac12$.
>
> **Result:** $\dfrac{x^3+2}{x^2-1}=x+\dfrac{3/2}{x-1}-\dfrac{1/2}{x+1}$.
>
> **Check** at $x=2$: LHS $=10/3=3.333$; RHS $=2+1.5-\tfrac16=3.333$ ✓.

Partial fractions are the workhorse of Laplace and Z-transform inversion ([Ch. 14](../part-5-transforms/14-laplace-transforms.md), [Ch. 17](../part-5-transforms/17-z-transforms.md)) and of integrating rational functions ([Ch. 5](../part-2-calculus/05-integral-calculus.md)); in transfer functions the variable is conventionally $s$, and the coefficients $A$ are called **residues**.

> [!tip] Engineering Note — partial fractions as parallel decomposition
> Each term $\frac{A}{s-p}$ corresponds to a first-order mode $Ae^{pt}$ in the time response. A partial-fraction expansion therefore splits a high-order system into a **parallel bank of simple subsystems**, one per pole. Software such as MATLAB's `residue` or SciPy's `scipy.signal.residue` performs the decomposition numerically; the hand method remains essential for checking it and for symbolic work.

> [!warning] Common Pitfall
> Partial fractions apply only to **proper** fractions: $\dfrac{x^2}{x^2-1}$ must first become $1+\dfrac{1}{x^2-1}$. Also, an irreducible quadratic needs a **linear** numerator $Bx+C$, not a constant — using only $C$ leaves too few unknowns and the equations become inconsistent.

## 1.5 Exponential and logarithmic functions

### Why exponentials and logarithms

Whenever the rate of change of a quantity is proportional to the quantity itself, the result is exponential: a capacitor discharging through a resistor, current building up in an inductor, a radioactive source decaying, a population of charge carriers being swept out of a junction. The diode current grows exponentially with voltage, and the open-loop growth of an oscillator is exponential until it saturates. Logarithms, the inverse operation, turn multiplication into addition, which is why engineers measure gains in decibels, plot frequency responses on logarithmic axes, and quote signal-to-noise ratios as logarithms. These two functions are as fundamental to engineering as the four arithmetic operations.

### Exponential functions

> [!definition] Exponential function
> For a fixed base $a>0$, $a\ne1$, the **exponential function** $f(x)=a^x$ is defined for all real $x$ by the exponent laws of Section 1.2 (extended to irrational $x$ by continuity). It satisfies
> $$
> a^{x+y}=a^xa^y,\qquad a^0=1,\qquad a^x>0\ \text{for all }x.
> $$
> For $a>1$ it is strictly increasing, with $a^x\to\infty$ as $x\to\infty$ and $a^x\to0$ as $x\to-\infty$; for $0<a<1$ it is strictly decreasing.

The functional equation $f(x+y)=f(x)f(y)$ is the essential property: equal steps in $x$ multiply $f$ by equal factors. A quantity that halves every $5\ \text{ms}$ is $f(t)=f(0)\,2^{-t/5\,\text{ms}}$.

### The natural base $e$

Among all bases one is singled out by calculus. Consider interest at a nominal rate of $100\%$ per year compounded $n$ times a year: one unit grows to $\left(1+\frac1n\right)^n$.

| $n$ | 1 | 10 | 100 | 1000 | $10^6$ |
|---|---|---|---|---|---|
| $\left(1+\frac1n\right)^n$ | 2 | 2.59374 | 2.70481 | 2.71692 | 2.71828 |

> [!definition] Euler's number
> $$
> e=\lim_{n\to\infty}\left(1+\frac1n\right)^n=\sum_{k=0}^{\infty}\frac{1}{k!}=1+1+\frac12+\frac16+\frac1{24}+\dots=2.718281828\ldots
> $$
> The **natural exponential** is $\exp(x)=e^x$, and more generally $e^x=\lim_{n\to\infty}\left(1+\frac xn\right)^n=\sum_{k=0}^\infty\frac{x^k}{k!}$.

The series converges very fast: the first eight terms already give $2.71825$. What makes $e$ natural is that $e^x$ is the unique function equal to its own derivative with value $1$ at $x=0$ ([Ch. 4](../part-2-calculus/04-differential-calculus.md)). Consequently the solution of $\dot y=-y/\tau$ is $y=y(0)e^{-t/\tau}$, which is why $e$ — not 2 or 10 — appears in every first-order circuit ([Ch. 11](../part-4-differential-equations/11-ordinary-differential-equations.md)). Any other exponential can be rewritten in base $e$: $a^x=e^{x\ln a}$.

### Logarithms

> [!definition] Logarithm
> For $a>0$, $a\ne1$, the **logarithm to base $a$** is the inverse of $a^x$:
> $$
> y=a^x\iff x=\log_a y,\qquad y>0.
> $$
> The **natural logarithm** is $\ln x=\log_e x$; the **common logarithm** is $\log_{10}x$; $\log_2x$ is used in information theory and digital design.

So $\log_a y$ answers the question "to what power must $a$ be raised to give $y$?" Immediately, $\log_a1=0$, $\log_aa=1$, $a^{\log_a y}=y$ and $\log_a(a^x)=x$. The domain of $\log_a$ is $(0,\infty)$ and its range is all of $\mathbb{R}$.

> [!theorem] Laws of logarithms ($x,y>0$, any real $p$)
> $$
> \log_a(xy)=\log_a x+\log_a y,\qquad \log_a\frac{x}{y}=\log_a x-\log_a y,\qquad \log_a x^p=p\log_a x,
> $$
> $$
> \log_a x=\frac{\log_b x}{\log_b a}\quad\text{(change of base)},\qquad\text{in particular}\quad \log_a x=\frac{\ln x}{\ln a}.
> $$

*Proof.* Let $u=\log_ax$ and $v=\log_ay$, so $x=a^u$, $y=a^v$. Then $xy=a^{u+v}$, so $\log_a(xy)=u+v$; $x^p=a^{pu}$, so $\log_ax^p=pu$. For change of base, take $\log_b$ of $x=a^u$: $\log_bx=u\log_ba$. $\square$

Each law is just an exponent law read backwards. Two useful consequences: $\log_{10}2=0.30103$, so $2^{10}=1024\approx10^3$ ("ten bits per factor of a thousand"), and the number of bits needed to index $N$ items is $\lceil\log_2N\rceil$; a million items need $\lceil19.93\rceil=20$ bits.

```desmos-graph
left=-3; right=5; bottom=-3; top=5
---
y=e^{x}|blue
y=\ln(x)|red
y=x|grey|dashed
y=2^{x}|blue|dotted
y=\frac{\ln(x)}{\ln(2)}|red|dotted
```
*$e^x$ (blue) and $\ln x$ (red) are mirror images in $y=x$, as every function and its inverse must be (Section 1.9). Base-2 versions (dotted) have the same shapes. Exponentials all pass through $(0,1)$ and logarithms through $(1,0)$.*

### Solving exponential and logarithmic equations

The strategy is to isolate the exponential (or logarithm), then apply the inverse function. Because $\exp$ and $\ln$ are one-to-one, these steps are equivalences — provided arguments of logarithms are positive.

> [!example] Worked Example 1.26 — Exponential equations
> **(a)** Solve $2^{x+1}=3^x$. Take natural logs: $(x+1)\ln2=x\ln3$, so $x(\ln3-\ln2)=\ln2$ and
> $$
> x=\frac{\ln2}{\ln3-\ln2}=\frac{0.69315}{0.40547}=1.7095.
> $$
> **Check:** $2^{2.7095}=6.541$ and $3^{1.7095}=6.541$ ✓.
>
> **(b)** Solve $e^{2x}-3e^x+2=0$. With $u=e^x>0$: $u^2-3u+2=(u-1)(u-2)=0$, so $e^x=1$ or $e^x=2$, i.e. $x=0$ or $x=\ln2=0.6931$. (Had a root $u\le0$ appeared, it would be rejected, since $e^x>0$.)
>
> **(c)** Solve $\log_{10}(x)+\log_{10}(x-3)=1$. Combine: $\log_{10}\big(x(x-3)\big)=1$, so $x^2-3x-10=0$, $(x-5)(x+2)=0$. The candidate $x=-2$ is rejected (logarithm of a negative number), leaving $x=5$. **Check:** $\log_{10}5+\log_{10}2=\log_{10}10=1$ ✓.

### First-order transients

The step and natural responses of every first-order circuit are built from $e^{-t/\tau}$, where the **time constant** is $\tau=RC$ or $\tau=L/R$. Useful landmarks follow directly from the logarithm:

| Elapsed time | $e^{-t/\tau}$ | $1-e^{-t/\tau}$ |
|---|---|---|
| $\tau$ | $0.368$ | $0.632$ |
| $2\tau$ | $0.135$ | $0.865$ |
| $3\tau$ | $0.050$ | $0.950$ |
| $5\tau$ | $0.0067$ | $0.9933$ |
| $\tau\ln2=0.693\tau$ | $0.5$ | $0.5$ |

> [!example] Worked Example 1.27 — Capacitor discharge and rise time
> **(a)** A capacitor discharges as $v(t)=V_0e^{-t/\tau}$ with $R=10\ \text{k}\Omega$, $C=100\ \text{nF}$. How long until $v$ falls to $10\%$ of $V_0$?
> $\tau=10^4\times10^{-7}=10^{-3}\ \text{s}$. Set $e^{-t/\tau}=0.1$ and take $\ln$: $-t/\tau=\ln0.1=-\ln10$, so
> $$
> t=\tau\ln10=2.303\ \text{ms}.
> $$
> **(b)** The **10–90 % rise time** of a charging response $v=V_\infty(1-e^{-t/\tau})$: solve $1-e^{-t_1/\tau}=0.1$ and $1-e^{-t_2/\tau}=0.9$, giving $t_1=\tau\ln\frac{1}{0.9}$ and $t_2=\tau\ln10$. Hence
> $$
> t_r=t_2-t_1=\tau\left(\ln10-\ln\tfrac{10}{9}\right)=\tau\ln9=2.197\,\tau.
> $$
> **Check:** $e^{-2.303}=0.100$ ✓ and $\ln 9=2.1972$ ✓. Combined with the bandwidth of the same circuit, $f_{3\text{dB}}=1/(2\pi\tau)$, this gives the well-known rule $t_rf_{3\text{dB}}\approx0.35$.

```desmos-graph
left=-0.5; right=6; bottom=-0.1; top=1.15
---
y=1-e^{-x}|blue|0<x<6
y=e^{-x}|red|0<x<6
y=1|grey|dashed
y=1-x|red|dotted|0<x<1
x=1|grey|dotted
```
*Charging $1-e^{-t/\tau}$ (blue) and discharging $e^{-t/\tau}$ (red) against $t/\tau$. The dotted tangent at $t=0$ reaches zero at exactly one time constant; after $5\tau$ the transient is within $1\%$ of its final value.*

> [!example] Worked Example 1.28 — Inverting the diode equation
> The Shockley diode equation is $I=I_S\left(e^{V/(nV_T)}-1\right)$. Find $V$ for $I=1\ \text{mA}$ when $I_S=10^{-14}\ \text{A}$, $n=1$, $V_T=25.85\ \text{mV}$, and find how much $V$ changes per decade of current.
>
> **Step 1 — isolate the exponential.** $e^{V/(nV_T)}=1+I/I_S$, so
> $$
> V=nV_T\ln\left(1+\frac{I}{I_S}\right).
> $$
> **Step 2 — evaluate.** $I/I_S=10^{11}$, $\ln(10^{11}+1)=25.328$, so $V=0.02585\times25.328=0.6547\ \text{V}$.
>
> **Step 3 — per decade.** For $I\gg I_S$, $V\approx nV_T\ln(I/I_S)$, and multiplying $I$ by 10 adds $nV_T\ln10=59.5\ \text{mV}$.
>
> **Check:** $10^{-14}\left(e^{0.6547/0.02585}-1\right)\approx1.00\times10^{-3}\ \text{A}$ ✓. The "60 mV per decade" rule of analogue design is just $\ln10$ in disguise.

### Decibels and logarithmic scales

Power ratios in electronics span many orders of magnitude — a radio receiver may handle signals from $10^{-15}\ \text{W}$ to $10^{-3}\ \text{W}$ — so engineers use the **decibel**:

> [!definition] Decibel
> $$
> G_{\text{dB}}=10\log_{10}\frac{P_2}{P_1}=20\log_{10}\frac{V_2}{V_1}\quad(\text{voltages across equal resistances}).
> $$
> Absolute levels use a stated reference: dBm is relative to $1\ \text{mW}$, dBW to $1\ \text{W}$, dBV to $1\ \text{V}$, dB$\mu$V to $1\ \mu\text{V}$.

The factor 20 for voltage comes from $P\propto V^2$ and the law $\log x^2=2\log x$. Useful anchors: $\times2$ in power $=+3.01$ dB; $\times2$ in voltage $=+6.02$ dB; $\times10$ in voltage $=+20$ dB; $1/\sqrt2$ in voltage $=-3.01$ dB (the "half-power" or "$-3$ dB" point that defines bandwidth). Since $0\ \text{dBm}=1\ \text{mW}$, $+30\ \text{dBm}=1\ \text{W}$.

Gains in cascade **multiply**, so their dB values **add** — the whole point of the logarithm.

> [!example] Worked Example 1.29 — A receiver chain
> A signal of $-90\ \text{dBm}$ passes through a low-noise amplifier ($+15$ dB), a mixer ($-7$ dB conversion loss) and an IF amplifier ($+25$ dB). Find the output power in dBm and in watts, and the overall voltage gain factor.
>
> **Step 1 — add.** $-90+15-7+25=-57\ \text{dBm}$.
>
> **Step 2 — convert.** $P=10^{-57/10}\ \text{mW}=1.995\times10^{-6}\ \text{mW}=1.995\ \text{nW}$.
>
> **Step 3 — overall gain.** $+33$ dB; as a power ratio $10^{3.3}=1995$; as a voltage ratio (equal impedances) $10^{33/20}=44.7$.
>
> **Check:** multiplying the linear power gains $31.62\times0.1995\times316.2=1995$ ✓.

**Logarithmic axes.** Plotting $\log y$ against $\log x$ turns a **power law** $y=ax^b$ into a straight line, $\log y=\log a+b\log x$, of slope $b$. Plotting $\log y$ against $x$ (semi-log) turns an exponential $y=ae^{kx}$ into a straight line of slope $k\log e$. Bode plots ([Ch. 16](../part-5-transforms/16-fourier-transforms.md)) use dB against $\log\omega$, so that a pole contributes a straight-line slope of $-20$ dB per decade.

> [!example] Worked Example 1.30 — Fitting a power law from two points
> A measured quantity gives $y=5.0$ at $x=2$ and $y=158.1$ at $x=20$. Assuming $y=ax^b$, find $a$ and $b$.
>
> **Step 1 — slope on log–log axes.** $b=\dfrac{\log_{10}158.1-\log_{10}5.0}{\log_{10}20-\log_{10}2}=\dfrac{\log_{10}31.62}{1}=1.500$.
>
> **Step 2 — intercept.** $a=5.0/2^{1.5}=1.768$.
>
> **Check:** $1.768\times20^{1.5}=1.768\times89.44=158.1$ ✓.

### Hyperbolic functions

Certain combinations of exponentials occur so often — in transmission-line theory, catenary cables, and the solution of $y''=k^2y$ — that they have names.

> [!definition] Hyperbolic functions
> $$
> \cosh x=\frac{e^x+e^{-x}}{2},\qquad \sinh x=\frac{e^x-e^{-x}}{2},\qquad \tanh x=\frac{\sinh x}{\cosh x}=\frac{e^{x}-e^{-x}}{e^{x}+e^{-x}}.
> $$

> [!theorem] Hyperbolic identities
> $$
> \cosh^2x-\sinh^2x=1,\qquad \cosh x+\sinh x=e^x,\qquad \sinh2x=2\sinh x\cosh x,\qquad \cosh2x=\cosh^2x+\sinh^2x.
> $$

*Proof of the first.* $\cosh^2x-\sinh^2x=\frac14\left[(e^{2x}+2+e^{-2x})-(e^{2x}-2+e^{-2x})\right]=1$. $\square$

They mirror the trigonometric identities with some sign changes — no coincidence: with complex arguments $\cosh(jx)=\cos x$ and $\sinh(jx)=j\sin x$ ([Ch. 2](02-complex-numbers.md)). The parametrisation $(\cosh t,\sinh t)$ traces the hyperbola $x^2-y^2=1$, just as $(\cos t,\sin t)$ traces the circle. Note $\cosh$ is even, $\sinh$ and $\tanh$ are odd, $\cosh x\ge1$, and $\tanh x\in(-1,1)$ — the soft-saturation shape used to model a differential pair (Section 1.8). Numerically, $\cosh1=1.5431$, $\sinh1=1.1752$, $\tanh1=0.7616$.

```desmos-graph
left=-3; right=3; bottom=-4; top=5
---
y=\cosh(x)|blue
y=\sinh(x)|red
y=\tanh(x)|green
y=\frac{e^{x}}{2}|grey|dashed
```
*$\cosh x$ (blue), $\sinh x$ (red) and $\tanh x$ (green). For large $x$ both $\cosh x$ and $\sinh x$ approach $e^x/2$ (grey dashed), while $\tanh x$ saturates at $\pm1$.*

> [!warning] Common Pitfall
> $\ln(x+y)\neq\ln x+\ln y$, $\ \ln(xy)\ne\ln x\cdot\ln y$, and $\ln$ is undefined for $x\le0$. Before taking logs of both sides of an equation, make sure both sides are positive. In decibels, never use $20\log$ for a power ratio or $10\log$ for a voltage ratio, and remember that "$\text{dB}$" alone is a ratio while "dBm" is an absolute power.

> [!tip] Engineering Note — logarithms of dimensional quantities
> Strictly, only dimensionless numbers have logarithms (Section 1.10). An expression like $\ln I$ is shorthand for $\ln(I/I_0)$ with some reference $I_0$, and $\ln I_1-\ln I_2=\ln(I_1/I_2)$ is then reference-independent. This is exactly why decibel levels always state a reference (dBm, dBV): without it the number is meaningless.

## 1.6 Trigonometric functions and identities

### Why trigonometry

Alternating current, radio carriers, audio tones, rotating machines and the components of a Fourier series are all sinusoids. Angles also enter through geometry: the phase difference between two phasors, the direction of a field vector, the bearing of an antenna beam. Trigonometric identities are the algebra of these quantities — they let a product of two signals be rewritten as a sum of tones (the principle of every mixer), a sum of a sine and a cosine be rewritten as a single shifted sinusoid, and a squared signal be rewritten in terms of its average power. This section gathers the functions and identities used throughout the book; [Ch. 3](03-trigonometry-and-periodic-functions.md) develops sinusoids and periodic signals in depth, and [Ch. 2](02-complex-numbers.md) shows how Euler's formula derives every identity below from one exponential law.

### Angles and radians

> [!definition] Radian
> The **radian** measure of an angle is the length of the arc it subtends on a circle divided by the radius: $\theta=s/r$. A full turn is $2\pi$ rad, so
> $$
> \pi\ \text{rad}=180^\circ,\qquad 1\ \text{rad}=57.2958^\circ,\qquad 1^\circ=0.0174533\ \text{rad}.
> $$

Being a ratio of two lengths, the radian is dimensionless, which is what allows $\sin\theta$ to have a power-series expansion (Section 1.10). In calculus, angles are **always** in radians: the derivative $\frac{d}{d\theta}\sin\theta=\cos\theta$ holds only in radians (in degrees it acquires a factor $\pi/180$). The arc length $s=r\theta$ and the sector area $A=\tfrac12r^2\theta$ also assume radians.

### Definitions on the unit circle

> [!definition] Sine, cosine and tangent
> Let the ray at angle $\theta$, measured anticlockwise from the positive $x$-axis, meet the unit circle $x^2+y^2=1$ at the point $P$. Then
> $$
> P=(\cos\theta,\ \sin\theta),\qquad \tan\theta=\frac{\sin\theta}{\cos\theta}\quad(\cos\theta\ne0).
> $$
> The reciprocal functions are $\sec\theta=1/\cos\theta$, $\csc\theta=1/\sin\theta$, $\cot\theta=1/\tan\theta$.

For an acute angle in a right triangle this agrees with "opposite over hypotenuse" etc., but the unit-circle definition works for every real $\theta$, including negative angles and angles beyond $2\pi$. Immediately:

- $\sin$ and $\cos$ are **periodic** with period $2\pi$; $\tan$ has period $\pi$.
- $-1\le\sin\theta\le1$ and $-1\le\cos\theta\le1$; $\tan$ takes every real value.
- $\cos$ is **even** and $\sin$, $\tan$ are **odd**: $\cos(-\theta)=\cos\theta$, $\sin(-\theta)=-\sin\theta$.
- Shifts: $\sin\left(\theta+\frac\pi2\right)=\cos\theta$, $\cos\left(\theta+\frac\pi2\right)=-\sin\theta$, $\sin(\theta+\pi)=-\sin\theta$, $\cos(\pi-\theta)=-\cos\theta$, $\sin(\pi-\theta)=\sin\theta$.
- Signs by quadrant: all positive in the first; only $\sin$ in the second; only $\tan$ in the third; only $\cos$ in the fourth.

The exact values below follow from the $45^\circ$ isosceles right triangle and the half-equilateral $30^\circ$–$60^\circ$ triangle.

| $\theta$ | $0$ | $\pi/6$ | $\pi/4$ | $\pi/3$ | $\pi/2$ | $2\pi/3$ | $3\pi/4$ | $\pi$ |
|---|---|---|---|---|---|---|---|---|
| $\sin\theta$ | $0$ | $1/2$ | $\sqrt2/2$ | $\sqrt3/2$ | $1$ | $\sqrt3/2$ | $\sqrt2/2$ | $0$ |
| $\cos\theta$ | $1$ | $\sqrt3/2$ | $\sqrt2/2$ | $1/2$ | $0$ | $-1/2$ | $-\sqrt2/2$ | $-1$ |
| $\tan\theta$ | $0$ | $1/\sqrt3$ | $1$ | $\sqrt3$ | undef. | $-\sqrt3$ | $-1$ | $0$ |

```desmos-graph
left=-7; right=7; bottom=-2.5; top=2.5
---
y=\sin(x)|blue
y=\cos(x)|red
y=\tan(x)|green|dashed
x=1.5708|grey|dotted
x=-1.5708|grey|dotted
```
*$\sin x$ (blue), $\cos x$ (red, the sine shifted left by $\pi/2$) and $\tan x$ (green, period $\pi$, vertical asymptotes where $\cos x=0$, marked at $\pm\pi/2$).*

### The fundamental identities

> [!theorem] Pythagorean identities
> $$
> \sin^2\theta+\cos^2\theta=1,\qquad 1+\tan^2\theta=\sec^2\theta,\qquad 1+\cot^2\theta=\csc^2\theta.
> $$

*Proof.* The point $(\cos\theta,\sin\theta)$ lies on $x^2+y^2=1$. Divide by $\cos^2\theta$ or $\sin^2\theta$ for the other two. $\square$

> [!theorem] Compound-angle formulas
> $$
> \begin{aligned}
> \cos(A\pm B)&=\cos A\cos B\mp\sin A\sin B\\
> \sin(A\pm B)&=\sin A\cos B\pm\cos A\sin B\\
> \tan(A\pm B)&=\frac{\tan A\pm\tan B}{1\mp\tan A\tan B}
> \end{aligned}
> $$

*Proof of $\cos(A-B)$.* The points $P=(\cos A,\sin A)$ and $Q=(\cos B,\sin B)$ on the unit circle are separated by the angle $A-B$. By the distance formula,
$$
|PQ|^2=(\cos A-\cos B)^2+(\sin A-\sin B)^2=2-2(\cos A\cos B+\sin A\sin B).
$$
Rotating the picture so that $Q$ lies at $(1,0)$ moves $P$ to $(\cos(A-B),\sin(A-B))$ without changing distances, so also $|PQ|^2=(\cos(A-B)-1)^2+\sin^2(A-B)=2-2\cos(A-B)$. Comparing gives $\cos(A-B)=\cos A\cos B+\sin A\sin B$. Replacing $B$ by $-B$ gives $\cos(A+B)$; replacing $A$ by $\frac\pi2-A$ converts cosines to sines and gives the sine formulas; dividing gives the tangent formula. $\square$

Every other identity in this section follows from these by algebra.

> [!theorem] Double-angle and power-reduction formulas
> $$
> \begin{aligned}
> \sin2\theta&=2\sin\theta\cos\theta\\
> \cos2\theta&=\cos^2\theta-\sin^2\theta=1-2\sin^2\theta=2\cos^2\theta-1\\
> \sin^2\theta&=\tfrac12(1-\cos2\theta),\qquad \cos^2\theta=\tfrac12(1+\cos2\theta)
> \end{aligned}
> $$

(Set $A=B=\theta$ in the compound formulas, then use $\sin^2+\cos^2=1$.)

> [!theorem] Product-to-sum and sum-to-product
> $$
> \begin{aligned}
> \cos A\cos B&=\tfrac12[\cos(A-B)+\cos(A+B)]\\
> \sin A\sin B&=\tfrac12[\cos(A-B)-\cos(A+B)]\\
> \sin A\cos B&=\tfrac12[\sin(A-B)+\sin(A+B)]\\
> \sin P+\sin Q&=2\sin\tfrac{P+Q}2\cos\tfrac{P-Q}2,\qquad \cos P+\cos Q=2\cos\tfrac{P+Q}2\cos\tfrac{P-Q}2
> \end{aligned}
> $$

Each product formula follows by adding or subtracting a pair of compound-angle formulas (e.g. $\cos(A-B)+\cos(A+B)=2\cos A\cos B$). The sum-to-product forms are the same identities read backwards with $P=A+B$, $Q=A-B$.

> [!example] Worked Example 1.31 — Exact values from identities
> **(a)** $\cos75^\circ=\cos(45^\circ+30^\circ)=\dfrac{\sqrt2}{2}\cdot\dfrac{\sqrt3}{2}-\dfrac{\sqrt2}{2}\cdot\dfrac12=\dfrac{\sqrt6-\sqrt2}{4}$.
> **Check:** $(2.44949-1.41421)/4=0.25882=\cos75^\circ$ ✓.
>
> **(b)** If $\sin\theta=\tfrac35$ with $\theta$ in the second quadrant, find $\cos\theta$ and $\sin2\theta$. From $\cos^2\theta=1-\tfrac9{25}=\tfrac{16}{25}$ and $\cos\theta<0$ in quadrant II, $\cos\theta=-\tfrac45$. Then $\sin2\theta=2\cdot\tfrac35\cdot\left(-\tfrac45\right)=-\tfrac{24}{25}$.
> **Check:** $\theta=\pi-\arcsin0.6=2.4981$ rad, $\sin(4.9962)=-0.96=-24/25$ ✓.

### Harmonic addition: $a\cos\theta+b\sin\theta$

A sum of a cosine and a sine at the **same** frequency is a single sinusoid with a phase shift. This identity is the bridge between the rectangular and polar descriptions of a sinusoid.

> [!theorem] Harmonic addition
> For real $a,b$ not both zero,
> $$
> a\cos\theta+b\sin\theta=R\cos(\theta-\alpha),\qquad R=\sqrt{a^2+b^2},\quad \alpha=\operatorname{atan2}(b,a).
> $$

*Derivation.* Expand the right side: $R\cos\theta\cos\alpha+R\sin\theta\sin\alpha$. Matching coefficients requires $R\cos\alpha=a$ and $R\sin\alpha=b$, i.e. $(a,b)$ has polar coordinates $(R,\alpha)$ (Section 1.7). Squaring and adding gives $R^2=a^2+b^2$; the quadrant of $\alpha$ is fixed by the signs of $a$ and $b$, which is exactly what $\operatorname{atan2}$ provides. $\square$

> [!example] Worked Example 1.32 — Combining two components
> Write $v(t)=3\cos\omega t+4\sin\omega t$ as a single cosine.
>
> **Step 1.** $R=\sqrt{9+16}=5$.
>
> **Step 2.** $\alpha=\operatorname{atan2}(4,3)=0.9273$ rad $=53.13^\circ$ (first quadrant, since $a,b>0$).
>
> **Result:** $v(t)=5\cos(\omega t-0.9273)$: amplitude $5$, lagging the reference cosine by $53.13^\circ$.
>
> **Check** at $\omega t=0.5$: $3\cos0.5+4\sin0.5=4.5504$; $5\cos(0.5-0.9273)=4.5504$ ✓.
>
> Had the signal been $-3\cos\omega t+4\sin\omega t$, $\alpha=\operatorname{atan2}(4,-3)=126.87^\circ$ — **not** $\arctan(4/(-3))=-53.13^\circ$, which would put the phase in the wrong quadrant.

```desmos-graph
left=-1; right=13; bottom=-6; top=6
---
y=3\cos(x)|blue|dashed
y=4\sin(x)|red|dashed
y=3\cos(x)+4\sin(x)|purple
y=5|grey|dotted
(0.9273,5)|label:peak at 0.927 rad
```
*$3\cos x$ (blue) plus $4\sin x$ (red) is a single sinusoid (purple) of amplitude 5, peaking at $x=\alpha=0.927$ rad: harmonic addition made visible.*

### Products of sinusoids: mixing, beats and power

> [!example] Worked Example 1.33 — A mixer and beats
> **(a) Mixing.** Multiplying two tones,
> $$
> \cos(\omega_1t)\cos(\omega_2t)=\tfrac12\cos\big((\omega_1-\omega_2)t\big)+\tfrac12\cos\big((\omega_1+\omega_2)t\big),
> $$
> produces components at the **difference** and **sum** frequencies. With $f_1=10.7$ MHz and $f_2=10.0$ MHz the outputs are at $0.7$ MHz and $20.7$ MHz. This identity is the entire mathematical content of frequency conversion in a superheterodyne receiver.
>
> **(b) Beats.** Adding two nearby tones, sum-to-product gives
> $$
> \cos(2\pi\cdot440\,t)+\cos(2\pi\cdot444\,t)=2\cos(2\pi\cdot2\,t)\cos(2\pi\cdot442\,t):
> $$
> a $442$ Hz tone whose envelope $|2\cos(2\pi\cdot2t)|$ swells and fades **4 times per second** (the envelope magnitude has period $0.25$ s).
>
> **Check (b)** at $t=0.1$ s: LHS $=\cos(88\pi)+\cos(88.8\pi)=1+\cos(0.8\pi)=1-0.80902=0.19098$; RHS $=2\cos(0.4\pi)\cos(88.4\pi)=2(0.30902)(0.30902)=0.19098$ ✓.

> [!example] Worked Example 1.34 — Average power of a sinusoid
> A voltage $v(t)=V_m\cos\omega t$ across a resistor $R$ delivers instantaneous power
> $$
> p(t)=\frac{v^2}{R}=\frac{V_m^2}{R}\cos^2\omega t=\frac{V_m^2}{2R}\big(1+\cos2\omega t\big)
> $$
> by power reduction. The $\cos2\omega t$ term averages to zero over a period, so the mean power is $\bar P=\dfrac{V_m^2}{2R}$. Writing $\bar P=V_{\text{rms}}^2/R$ defines the **RMS value** $V_{\text{rms}}=V_m/\sqrt2$. For the $230\ \text{V}$ (RMS) mains, $V_m=230\sqrt2=325.3\ \text{V}$. Note also that the power fluctuates at **twice** the supply frequency — the origin of the $100$/$120$ Hz hum in poorly filtered supplies. (Averages of functions are treated formally in [Ch. 5](../part-2-calculus/05-integral-calculus.md).)

### Solving trigonometric equations

Because trigonometric functions are periodic, a solvable trigonometric equation has **infinitely many** solutions. Find all solutions in one period, then add multiples of the period.

> [!theorem] General solutions
> For an integer $k$:
> $$
> \sin x=\sin\alpha\iff x=\alpha+2k\pi\ \text{ or }\ x=\pi-\alpha+2k\pi;\qquad
> \cos x=\cos\alpha\iff x=\pm\alpha+2k\pi;\qquad
> \tan x=\tan\alpha\iff x=\alpha+k\pi.
> $$

> [!example] Worked Example 1.35 — A quadratic in $\sin x$
> Solve $2\sin^2x-\sin x-1=0$.
>
> **Step 1.** With $u=\sin x$: $2u^2-u-1=(2u+1)(u-1)=0$, so $\sin x=1$ or $\sin x=-\tfrac12$.
>
> **Step 2.** $\sin x=1\Rightarrow x=\tfrac\pi2+2k\pi$. $\ \sin x=-\tfrac12=\sin\left(-\tfrac\pi6\right)\Rightarrow x=-\tfrac\pi6+2k\pi$ or $x=\pi+\tfrac\pi6+2k\pi=\tfrac{7\pi}6+2k\pi$.
>
> **In $[0,2\pi)$:** $x\in\left\{\tfrac\pi2,\ \tfrac{7\pi}6,\ \tfrac{11\pi}6\right\}$.
>
> **Check** $x=7\pi/6$: $\sin x=-0.5$, $2(0.25)+0.5-1=0$ ✓.

### Small-angle approximations

For $|\theta|\ll1$ (in radians),

$$
\sin\theta\approx\theta,\qquad \tan\theta\approx\theta,\qquad \cos\theta\approx1-\frac{\theta^2}{2}.
$$

These are the first terms of the Taylor series ([Ch. 4](../part-2-calculus/04-differential-calculus.md)). At $\theta=0.1$ rad ($5.7^\circ$): $\sin0.1=0.099833$ (error $0.17\%$), $\tan0.1=0.100335$ (error $0.33\%$), $\cos0.1=0.995004$ versus $0.995$. They linearise pendulums, phase-locked loops near lock, and small phase errors in modulators.

### Triangles: the laws of sines and cosines

For a triangle with sides $a,b,c$ opposite angles $A,B,C$:

$$
\frac{a}{\sin A}=\frac{b}{\sin B}=\frac{c}{\sin C},\qquad c^2=a^2+b^2-2ab\cos C.
$$

The law of cosines generalises Pythagoras ($C=90^\circ$) and is the quickest way to add two phasors graphically.

> [!example] Worked Example 1.36 — Adding two phasors geometrically
> Two voltages of magnitude $10$ V and $6$ V differ in phase by $60^\circ$. Find the magnitude and phase of their sum relative to the $10$ V phasor.
>
> **Step 1 — triangle.** Placing the phasors head to tail, the angle *inside* the triangle between them is $180^\circ-60^\circ=120^\circ$. The law of cosines gives
> $$
> |V|^2=10^2+6^2-2(10)(6)\cos120^\circ=100+36+60=196,\qquad |V|=14\ \text{V}.
> $$
> **Step 2 — angle.** By the law of sines, $\sin\beta=\dfrac{6\sin120^\circ}{14}=0.3712$, so $\beta=21.79^\circ$.
>
> **Check by components:** $10+6\cos60^\circ=13$, $6\sin60^\circ=5.196$; $\sqrt{13^2+5.196^2}=14.00$ ✓ and $\operatorname{atan2}(5.196,13)=21.79^\circ$ ✓. [Ch. 2](02-complex-numbers.md) does such sums with complex numbers, which avoids the triangle entirely.

> [!warning] Common Pitfall
> Calculators default to degrees or radians unpredictably — $\sin(30)$ is $0.5$ in degree mode and $-0.988$ in radian mode. Programming languages always use radians. Also, $\sin(A+B)\ne\sin A+\sin B$ and $\sin2\theta\ne2\sin\theta$; and an equation like $\sin x=0.5$ has infinitely many solutions, not just the one the calculator returns.

> [!tip] Engineering Note — identities as circuit operations
> Product-to-sum is what a nonlinear device (a diode or transistor mixer) does to two inputs; power reduction explains why squaring a sinusoid produces DC plus a second harmonic, the basis of RMS detectors and frequency doublers; harmonic addition is how two quadrature signals ($I\cos\omega t+Q\sin\omega t$) carry amplitude and phase in a digital radio.

## 1.7 Coordinate systems

### Why several coordinate systems

A coordinate system is a rule that labels each point of the plane or of space by numbers. The Cartesian system is the default, but it is not always the most convenient: the field of a long straight wire depends only on the distance from the wire, the radiation pattern of an antenna depends on direction, and a complex number is often best described by a length and an angle. Choosing coordinates that match the **symmetry** of a problem can turn an intractable calculation into a short one. This section fixes the conventions used throughout the book — cylindrical $(r,\phi,z)$ and spherical $(R,\theta,\phi)$, as in [Notation & Conventions](../00-front/notation-and-conventions.md) — and gives the conversion formulas; the calculus in these coordinates (Jacobians, gradient, divergence) comes in [Ch. 6](../part-2-calculus/06-multivariable-calculus.md) and [Ch. 7](../part-2-calculus/07-vector-calculus.md).

### The Cartesian plane

A point $P$ is labelled $(x,y)$ by its signed distances from two perpendicular axes. The basic formulas follow from Pythagoras:

$$
d=\sqrt{(x_2-x_1)^2+(y_2-y_1)^2},\qquad M=\left(\frac{x_1+x_2}{2},\ \frac{y_1+y_2}{2}\right),
$$

for the distance between two points and their midpoint.

> [!definition] Straight lines
> - **Slope** through $(x_1,y_1)$ and $(x_2,y_2)$: $m=\dfrac{y_2-y_1}{x_2-x_1}$ (undefined for a vertical line).
> - **Slope–intercept form:** $y=mx+c$. **Point–slope form:** $y-y_1=m(x-x_1)$. **General form:** $ax+by+c=0$.
> - Lines are **parallel** iff $m_1=m_2$ and **perpendicular** iff $m_1m_2=-1$.
> - **Distance from a point** $(x_0,y_0)$ to the line $ax+by+c=0$:
> $$
> d=\frac{|ax_0+by_0+c|}{\sqrt{a^2+b^2}}.
> $$

> [!definition] Circle
> The circle with centre $(h,k)$ and radius $\rho_0>0$ is the set of points at distance $\rho_0$ from the centre:
> $$
> (x-h)^2+(y-k)^2=\rho_0^2.
> $$

> [!example] Worked Example 1.37 — Lines and circles
> **(a) Load line.** A transistor's collector circuit satisfies $V_{CE}=V_{CC}-I_CR_C$. With $V_{CC}=12\ \text{V}$, $R_C=2\ \text{k}\Omega$, plotted with $I_C$ (mA) on the vertical axis and $V_{CE}$ (V) on the horizontal, this is the line $I_C=6-0.5V_{CE}$: slope $-1/R_C=-0.5\ \text{mA/V}$, intercepts $(0,6)$ and $(12,0)$. The operating point must lie on it.
>
> **(b) Circle from a general equation.** Find the centre and radius of $x^2+y^2-4x+6y-12=0$. Complete the square in each variable:
> $$
> (x^2-4x+4)+(y^2+6y+9)=12+4+9\quad\Longrightarrow\quad(x-2)^2+(y+3)^2=25.
> $$
> Centre $(2,-3)$, radius $5$.
>
> **(c) Point-to-line distance.** The distance from $(1,2)$ to $3x+4y-1=0$ is $\dfrac{|3+8-1|}{\sqrt{9+16}}=\dfrac{10}{5}=2$.
>
> **Check (b):** the point $(7,-3)$ is $5$ units right of the centre; $49+9-28-18-12=0$ ✓.

### Plane polar coordinates

> [!definition] Polar coordinates
> A point is labelled $(r,\theta)$, where $r\ge0$ is its distance from the origin (the **pole**) and $\theta$ the angle from the positive $x$-axis, anticlockwise positive.
> $$
> x=r\cos\theta,\quad y=r\sin\theta;\qquad r=\sqrt{x^2+y^2},\quad \theta=\operatorname{atan2}(y,x).
> $$

The angle is not unique — $(r,\theta)$ and $(r,\theta+2k\pi)$ are the same point — so a **principal value** is chosen, usually $\theta\in(-\pi,\pi]$. The origin has $r=0$ and any $\theta$.

> [!definition] Two-argument arctangent
> $$
> \operatorname{atan2}(y,x)=\begin{cases}
> \arctan(y/x), & x>0\\
> \arctan(y/x)+\pi, & x<0,\ y\ge0\\
> \arctan(y/x)-\pi, & x<0,\ y<0\\
> +\pi/2, & x=0,\ y>0\\
> -\pi/2, & x=0,\ y<0
> \end{cases}
> $$
> It returns the angle in $(-\pi,\pi]$ in the **correct quadrant** (undefined at the origin).

The single-argument $\arctan(y/x)$ only sees the ratio $y/x$, so it cannot distinguish $(1,1)$ from $(-1,-1)$; it always returns an angle in $(-\pi/2,\pi/2)$. Every programming language provides `atan2(y, x)` — note the order of arguments.

> [!example] Worked Example 1.38 — Polar conversion
> **(a)** $(3,4)$: $r=5$, $\theta=\operatorname{atan2}(4,3)=\arctan(4/3)=0.9273$ rad $=53.13^\circ$.
>
> **(b)** $(-3,-4)$: $r=5$, but $\theta=\operatorname{atan2}(-4,-3)=\arctan(4/3)-\pi=-2.2143$ rad $=-126.87^\circ$, not $53.13^\circ$.
>
> **(c)** $(r,\theta)=(2,\ 5\pi/6)$: $x=2\cos150^\circ=-\sqrt3=-1.732$, $y=2\sin150^\circ=1$.
>
> **Check (b):** $5\cos(-2.2143)=-3.000$, $5\sin(-2.2143)=-4.000$ ✓.

This is exactly the rectangular-to-polar conversion of complex numbers in [Ch. 2](02-complex-numbers.md), where $(x,y)\leftrightarrow x+jy$ and $(r,\theta)\leftrightarrow re^{j\theta}$.

**Curves in polar form.** Equations convert using $x=r\cos\theta$, $y=r\sin\theta$, $x^2+y^2=r^2$. For instance the circle $x^2+y^2=4x$ (centre $(2,0)$, radius 2) becomes $r^2=4r\cos\theta$, i.e. $r=4\cos\theta$. Polar plots are the natural way to show **radiation patterns**: the cardioid $r=1+\cos\theta$ is the pattern of a simple directional (cardioid) microphone or of two suitably phased antenna elements, with full response forward ($\theta=0$) and a null behind ($\theta=\pi$).

```desmos-graph
left=-2.5; right=2.5; bottom=-2; top=2
---
r=1+\cos(\theta)|blue
r=\operatorname{abs}\left(\cos(\theta)\right)|red|dashed
r=1|grey|dotted
```
*Polar patterns: the cardioid $r=1+\cos\theta$ (blue) has a single null at $\theta=\pi$; the figure-eight $r=\lvert\cos\theta\rvert$ (red) of a dipole-like element has nulls at $\pm\pi/2$; the grey circle is an omnidirectional pattern.*

### Three-dimensional Cartesian coordinates

A point in space is $(x,y,z)$ with respect to three mutually perpendicular axes forming a **right-handed** system: curling the fingers of the right hand from $+x$ towards $+y$, the thumb points along $+z$. The distance formula extends directly: $d=\sqrt{\Delta x^2+\Delta y^2+\Delta z^2}$.

### Cylindrical coordinates $(r,\phi,z)$

> [!definition] Cylindrical coordinates
> $r\ge0$ is the perpendicular distance from the $z$-axis, $\phi\in[0,2\pi)$ (or $(-\pi,\pi]$) the azimuthal angle measured from the $+x$ axis in the $xy$-plane, and $z$ the usual height.
> $$
> x=r\cos\phi,\quad y=r\sin\phi,\quad z=z;\qquad r=\sqrt{x^2+y^2},\quad \phi=\operatorname{atan2}(y,x).
> $$

Cylindrical coordinates are plane polar coordinates with a $z$ axis added. The surfaces $r=\text{const}$ are coaxial cylinders, $\phi=\text{const}$ half-planes hinged on the $z$-axis, and $z=\text{const}$ horizontal planes. This book uses $r$ (not $\rho$) for the cylindrical radius so that $\rho$ remains free for charge density.

### Spherical coordinates $(R,\theta,\phi)$

> [!definition] Spherical coordinates
> $R\ge0$ is the distance from the origin, $\theta\in[0,\pi]$ the **polar angle** measured from the $+z$ axis, and $\phi$ the **azimuth**, the same angle as in cylindrical coordinates.
> $$
> x=R\sin\theta\cos\phi,\quad y=R\sin\theta\sin\phi,\quad z=R\cos\theta;
> $$
> $$
> R=\sqrt{x^2+y^2+z^2},\qquad \theta=\arccos\frac{z}{R},\qquad \phi=\operatorname{atan2}(y,x).
> $$

*Derivation.* Drop a perpendicular from $P$ to the $xy$-plane. The right triangle formed by the origin, $P$ and the $z$-axis has hypotenuse $R$ and angle $\theta$ at the origin, so $z=R\cos\theta$ and the distance from the $z$-axis is $r=R\sin\theta$. Within the $xy$-plane, that distance $r$ and the azimuth $\phi$ give $x=r\cos\phi$, $y=r\sin\phi$. $\square$

The bridge between the two systems is therefore

$$
r=R\sin\theta,\qquad z=R\cos\theta,\qquad \phi\ \text{unchanged}.
$$

The surfaces $R=\text{const}$ are spheres, $\theta=\text{const}$ cones about the $z$-axis (the plane $z=0$ when $\theta=\pi/2$), and $\phi=\text{const}$ half-planes. On the Earth, $\theta$ is the **colatitude** ($90^\circ$ minus latitude) and $\phi$ the longitude.

| System | Coordinates | To Cartesian | Line element $ds^2$ | Volume element |
|---|---|---|---|---|
| Cartesian | $(x,y,z)$ | — | $dx^2+dy^2+dz^2$ | $dx\,dy\,dz$ |
| Cylindrical | $(r,\phi,z)$ | $x=r\cos\phi,\ y=r\sin\phi,\ z=z$ | $dr^2+r^2d\phi^2+dz^2$ | $r\,dr\,d\phi\,dz$ |
| Spherical | $(R,\theta,\phi)$ | $x=R\sin\theta\cos\phi,\ y=R\sin\theta\sin\phi,\ z=R\cos\theta$ | $dR^2+R^2d\theta^2+R^2\sin^2\theta\,d\phi^2$ | $R^2\sin\theta\,dR\,d\theta\,d\phi$ |

The line elements follow geometrically: in cylindrical coordinates a small change $d\phi$ moves a point along an arc of radius $r$, a distance $r\,d\phi$; in spherical coordinates $d\theta$ moves it along a circle of radius $R$ (distance $R\,d\theta$), and $d\phi$ along a circle of radius $R\sin\theta$ (distance $R\sin\theta\,d\phi$). Multiplying the three perpendicular displacements gives the volume elements. For example, integrating $R^2\sin\theta$ over $\theta\in[0,\pi]$ and $\phi\in[0,2\pi)$ gives the sphere's surface area $4\pi R^2$ ([Ch. 6](../part-2-calculus/06-multivariable-calculus.md)).

> [!example] Worked Example 1.39 — Converting points
> **(a) Cartesian to cylindrical.** $(x,y,z)=(-2,\ 2\sqrt3,\ 5)$: $r=\sqrt{4+12}=4$, $\phi=\operatorname{atan2}(2\sqrt3,-2)=2\pi/3$ (second quadrant), $z=5$. So $(r,\phi,z)=(4,\ 2\pi/3,\ 5)$.
>
> **(b) Cartesian to spherical.** $(1,1,\sqrt2)$: $R=\sqrt{1+1+2}=2$; $\theta=\arccos(\sqrt2/2)=\pi/4$; $\phi=\operatorname{atan2}(1,1)=\pi/4$.
>
> **(c) Spherical to cylindrical.** $(R,\theta,\phi)=(6,\ \pi/3,\ \pi/2)$: $r=6\sin\tfrac\pi3=5.196$, $z=6\cos\tfrac\pi3=3$, $\phi=\pi/2$.
>
> **Checks:** (a) $4\cos(2\pi/3)=-2$, $4\sin(2\pi/3)=3.464=2\sqrt3$ ✓. (b) $x=2\sin\tfrac\pi4\cos\tfrac\pi4=2\cdot\tfrac12=1$ ✓, $z=2\cos\tfrac\pi4=\sqrt2$ ✓. (c) $r^2+z^2=27+9=36=R^2$ ✓.

> [!example] Worked Example 1.40 — Distance between two points given in spherical form
> Find the distance between $P_1=(R,\theta,\phi)=(2,\ \pi/2,\ 0)$ and $P_2=(3,\ \pi/3,\ \pi/2)$.
>
> **Step 1 — convert.** $P_1=(2\cdot1\cdot1,\ 0,\ 0)=(2,0,0)$. $P_2=(3\sin\tfrac\pi3\cos\tfrac\pi2,\ 3\sin\tfrac\pi3\sin\tfrac\pi2,\ 3\cos\tfrac\pi3)=(0,\ 2.598,\ 1.5)$.
>
> **Step 2 — distance.** $d=\sqrt{2^2+2.598^2+1.5^2}=\sqrt{4+6.75+2.25}=\sqrt{13}=3.606$.
>
> **Check by the spherical law of cosines.** The angle $\gamma$ between the two position vectors satisfies $\cos\gamma=\sin\theta_1\sin\theta_2\cos(\phi_1-\phi_2)+\cos\theta_1\cos\theta_2=1\cdot0.866\cdot0+0=0$, so $d^2=R_1^2+R_2^2-2R_1R_2\cos\gamma=4+9=13$ ✓.
>
> Note that one may **not** subtract spherical coordinates component by component: $(3-2,\ \pi/3-\pi/2,\ \pi/2-0)$ is meaningless as a displacement.

> [!example] Worked Example 1.41 — Converting surfaces
> | Surface | Cartesian | Cylindrical | Spherical |
> |---|---|---|---|
> | Sphere of radius $a$ | $x^2+y^2+z^2=a^2$ | $r^2+z^2=a^2$ | $R=a$ |
> | Cylinder of radius 3 about $z$ | $x^2+y^2=9$ | $r=3$ | $R\sin\theta=3$ |
> | Cone, half-angle $45^\circ$ | $z=\sqrt{x^2+y^2}$ | $z=r$ | $\theta=\pi/4$ |
> | Half-plane $y=x,\ x>0$ | $y=x,\ x>0$ | $\phi=\pi/4$ | $\phi=\pi/4$ |
>
> For the cone: $z=r$ means $R\cos\theta=R\sin\theta$, so $\tan\theta=1$ and $\theta=\pi/4$. Each surface has the simplest equation in the system matched to its symmetry.

> [!warning] Common Pitfall — conventions differ
> Many mathematics texts swap the roles of $\theta$ and $\phi$ in spherical coordinates (using $\theta$ for azimuth), and some use $\rho$ for the spherical radius. This book follows the physics/electromagnetics convention: $\theta$ is the polar angle from $+z$, $\phi$ the azimuth. Also, the cylindrical $r$ (distance from the $z$-axis) is **not** the spherical $R$ (distance from the origin): $r=R\sin\theta$. Always check a source's convention before using its formulas.

> [!tip] Engineering Note — let symmetry choose
> Choose coordinates to match the geometry: cylindrical for coaxial cables, long wires, solenoids and optical fibres; spherical for point charges, antennas and far-field radiation; Cartesian for parallel-plate capacitors and rectangular waveguides. The field of a line charge, a messy expression in $(x,y,z)$, is simply proportional to $1/r$ in cylindrical coordinates; the far field of an antenna is a function of $(\theta,\phi)$ times $e^{-jkR}/R$.

## 1.8 Functions, domains, and ranges

### Why the function concept

A circuit, a sensor, an amplifier or a filter is, to a mathematician, a **function**: a rule that takes an input (a voltage, a temperature, a signal) and produces exactly one output. Signals themselves are functions of time. Speaking precisely about functions — what inputs are allowed, what outputs can occur, how graphs are transformed by delays and gains — is the vocabulary of signals and systems. Questions like "for what input range is this amplifier linear?" or "can this sensor reading be inverted to give a unique temperature?" are questions about domains, ranges and one-to-one functions.

### Definition

> [!definition] Function
> A **function** $f:A\to B$ assigns to each element $x$ of the set $A$ (the **domain**) exactly one element $f(x)$ of the set $B$ (the **codomain**). The **range** (or image) is the set of values actually taken,
> $$
> f(A)=\{f(x):x\in A\}\subseteq B.
> $$
> $x$ is the **independent variable** (argument, input) and $y=f(x)$ the **dependent variable** (value, output).

The essential word is *exactly one*. The relation $y^2=x$ does not define $y$ as a function of $x$, because $x=4$ gives $y=\pm2$; it defines two functions, $y=\sqrt x$ and $y=-\sqrt x$. Graphically, a curve in the $xy$-plane is the graph of a function of $x$ if and only if every vertical line meets it at most once — the **vertical-line test**.

A function is not the same thing as a formula. The same formula on different domains gives different functions (with different ranges and different inverses), and one function may need several formulas (a piecewise definition). Notation such as $f:\mathbb{R}\to\mathbb{R},\ x\mapsto x^2$ ("$x$ maps to $x^2$") states domain, codomain and rule together.

### Natural domain

When a formula is given without a stated domain, the **natural domain** is the largest set of reals for which it makes sense. Exclude:

- zeros of denominators;
- negative arguments of even roots ($\sqrt{\ },\ \sqrt[4]{\ },\dots$);
- non-positive arguments of logarithms;
- points where $\tan$ or $\sec$ is undefined ($\cos=0$), or $\cot$, $\csc$ ($\sin=0$);
- arguments outside $[-1,1]$ for $\arcsin$ and $\arccos$.

In a physical model the domain is further restricted by physics: a resistance is positive, a probability lies in $[0,1]$, a model of a transistor is valid only in its active region.

> [!example] Worked Example 1.42 — Natural domains
> **(a)** $f(x)=\dfrac{\ln(4-x^2)}{\sqrt{x+1}}$ needs $4-x^2>0$, i.e. $-2<x<2$, and $x+1>0$ (strictly, since the root is in a denominator), i.e. $x>-1$. Intersection: domain $(-1,2)$.
>
> **(b)** $g(x)=\sqrt{\dfrac{x-1}{x+2}}$ needs $\dfrac{x-1}{x+2}\ge0$, which by the sign chart of Example 1.15(c) gives $(-\infty,-2)\cup[1,\infty)$.
>
> **(c)** $h(x)=\arcsin(2x-1)$ needs $-1\le2x-1\le1$, i.e. $0\le x\le1$.

### Finding the range

The range is usually harder than the domain. Three systematic methods:

1. **Solve for $x$.** Set $y=f(x)$ and ask for which $y$ the equation has a solution $x$ in the domain.
2. **Complete the square** or otherwise rewrite $f$ in a form whose extreme values are obvious.
3. **Use monotonicity**: a continuous increasing function on $[a,b]$ has range $[f(a),f(b)]$. (Calculus gives a general method in [Ch. 4](../part-2-calculus/04-differential-calculus.md).)

> [!example] Worked Example 1.43 — Three ranges
> **(a) By completing the square.** $f(x)=x^2-4x+7$ on $[0,5]$. Since $f(x)=(x-2)^2+3$, the minimum is $3$ at $x=2$; the maximum over $[0,5]$ is at the end farther from $2$, $f(5)=12$ (while $f(0)=7$). Range $[3,12]$.
>
> **(b) By solving for $x$.** $f(x)=\dfrac{x}{x^2+1}$ on $\mathbb{R}$. Setting $y=\dfrac{x}{x^2+1}$ gives $yx^2-x+y=0$. If $y=0$, then $x=0$ works. If $y\ne0$, this quadratic in $x$ has a real root iff its discriminant $1-4y^2\ge0$, i.e. $|y|\le\tfrac12$. Range $\left[-\tfrac12,\tfrac12\right]$, with $y=\tfrac12$ at $x=1$ and $y=-\tfrac12$ at $x=-1$.
>
> **(c) From bounded functions.** $g(x)=2+3\sin x$: since $-1\le\sin x\le1$, the range is $[-1,5]$.
>
> **Check (b):** $f(1)=1/2$ ✓; $f(3)=0.3<0.5$ ✓.

### Properties worth naming

> [!definition] Even, odd, periodic, monotonic, bounded
> - **Even:** $f(-x)=f(x)$ for all $x$ (graph symmetric about the $y$-axis; $x^2$, $\cos x$, $\lvert x\rvert$, $\cosh x$).
> - **Odd:** $f(-x)=-f(x)$ (graph symmetric under a half-turn about the origin; $x^3$, $\sin x$, $\tan x$, $\sinh x$).
> - **Periodic** with period $T>0$: $f(x+T)=f(x)$ for all $x$; the smallest such $T$ is the **fundamental period**.
> - **Increasing** on an interval: $x_1<x_2\Rightarrow f(x_1)\le f(x_2)$ (**strictly** if $<$); **decreasing** similarly. **Monotonic** means one or the other.
> - **Bounded:** $|f(x)|\le M$ for some $M$ and all $x$ in the domain.

> [!theorem] Even–odd decomposition
> Every function $f$ defined on a symmetric domain splits uniquely as $f=f_e+f_o$ with $f_e$ even and $f_o$ odd:
> $$
> f_e(x)=\tfrac12\big[f(x)+f(-x)\big],\qquad f_o(x)=\tfrac12\big[f(x)-f(-x)\big].
> $$

*Proof.* The formulas clearly give an even and an odd function summing to $f$. For uniqueness, if $f=g_e+g_o$, then $f(-x)=g_e(x)-g_o(x)$; adding and subtracting recovers the formulas. $\square$

Products follow sign rules: even$\times$even and odd$\times$odd are even, even$\times$odd is odd. This is used constantly in Fourier analysis ([Ch. 15](../part-5-transforms/15-fourier-series.md)), where an even signal has only cosine terms and an odd one only sine terms.

> [!example] Worked Example 1.44 — Even and odd parts
> **(a)** For $f(x)=e^x$: $f_e=\tfrac12(e^x+e^{-x})=\cosh x$ and $f_o=\tfrac12(e^x-e^{-x})=\sinh x$, recovering $e^x=\cosh x+\sinh x$ (Section 1.5).
>
> **(b)** For the decaying exponential pulse $x(t)=e^{-t}u(t)$ (zero for $t<0$):
> $$
> x_e(t)=\tfrac12e^{-|t|},\qquad x_o(t)=\tfrac12\operatorname{sgn}(t)\,e^{-|t|},
> $$
> where $\operatorname{sgn}(t)$ is $+1$ for $t>0$ and $-1$ for $t<0$.
>
> **Check** at $t=1$: $x_e+x_o=\tfrac12e^{-1}+\tfrac12e^{-1}=e^{-1}=x(1)$ ✓; at $t=-1$: $\tfrac12e^{-1}-\tfrac12e^{-1}=0=x(-1)$ ✓.

### Transformations of graphs

Given the graph of $y=f(x)$, the graph of

$$
y=A\,f\big(b(x-c)\big)+d
$$

is obtained by: horizontal scaling by a factor $1/b$ (compression if $|b|>1$, with reflection in the $y$-axis if $b<0$); horizontal **shift** right by $c$; vertical scaling by $A$ (reflection in the $x$-axis if $A<0$); and vertical shift up by $d$. Changes **inside** the function act on $x$ and work "backwards" ($f(x-3)$ moves the graph *right*); changes **outside** act on $y$ and work as expected.

In signal language with $t$ for time: $x(t-t_0)$ is a **delay** by $t_0$ (for $t_0>0$), $x(at)$ is time compression ($a>1$) or expansion ($0<a<1$), $x(-t)$ is time reversal, and $Ax(t)$ is amplification.

> [!example] Worked Example 1.45 — Order of shift and scale
> Describe $y(t)=x(2t-4)$ in terms of $x(t)$.
>
> **Step 1 — factor the inner coefficient.** $x(2t-4)=x\big(2(t-2)\big)$.
>
> **Step 2 — interpret.** Compress $x$ by a factor 2 in time, then delay by $2$ — **not** by $4$. A feature of $x$ at time $t=\tau$ appears in $y$ at the time where $2t-4=\tau$, i.e. $t=(\tau+4)/2$.
>
> **Check:** if $x$ has a peak at $t=0$, then $y$ peaks where $2t-4=0$, i.e. $t=2$ ✓.

### Piecewise-defined functions

Many engineering functions are defined by different formulas on different sub-domains. Two that matter everywhere in EE are the unit step and the saturation (clipping) function:

$$
u(t)=\begin{cases}0,&t<0\\1,&t\ge0\end{cases}\qquad
\operatorname{sat}(x)=\begin{cases}-1,&x<-1\\ x,&|x|\le1\\ 1,&x>1\end{cases}
$$

(The value of $u(0)$ is a convention; some texts use $\tfrac12$. It never affects integrals.) Steps are the building blocks of piecewise signals: the rectangular pulse of width $T$ starting at $0$ is $u(t)-u(t-T)$, and the ramp is $r(t)=t\,u(t)$. Absolute value, $\operatorname{sgn}$, and $\max(0,x)$ (the ideal rectifier, and the "ReLU" of machine learning) are other common piecewise functions: $\max(0,x)=\tfrac12(x+|x|)$.

> [!example] Worked Example 1.46 — Building a signal from steps
> Express the signal that is $0$ for $t<0$, rises linearly to $2$ at $t=1$, stays at $2$ until $t=3$, then drops to $0$.
>
> **Step 1.** The ramp $2t\,u(t)$ produces the rise but keeps going.
>
> **Step 2.** At $t=1$ cancel the slope with $-2(t-1)u(t-1)$; for $t>1$ the sum is $2t-2(t-1)=2$ ✓.
>
> **Step 3.** At $t=3$ subtract the level: $-2u(t-3)$.
>
> **Result:** $x(t)=2t\,u(t)-2(t-1)\,u(t-1)-2\,u(t-3)$.
>
> **Check** at $t=0.5$: $1-0-0=1$ ✓; at $t=2$: $4-2-0=2$ ✓; at $t=4$: $8-6-2=0$ ✓.

```desmos-graph
left=-1; right=5; bottom=-0.5; top=2.6
---
y=2x|blue|0<x<1
y=2|blue|1<x<3
y=0|blue|x<0
y=0|blue|x>3
x=3|blue|dashed|0<y<2
```
*The trapezoidal-ramp signal of Example 1.46, built from a ramp, a delayed ramp and a delayed step. The dashed segment marks the jump at $t=3$.*

```desmos-graph
left=-3; right=3; bottom=-1.6; top=1.6
---
y=-1|blue|x<-1
y=x|blue|-1<x<1
y=1|blue|x>1
y=\tanh(2x)|red|dashed
```
*The ideal saturation function (blue) and the smooth $\tanh(2x)$ (red), often used to model an amplifier driven into its supply rails. Domain $\mathbb{R}$, range $[-1,1]$; the linear region is the domain sub-interval $\lvert x\rvert\le1$.*

### Functions of several variables

A function may take several inputs: the power $P(V,R)=V^2/R$ has domain $\{(V,R):R>0\}$ and the potential of a point charge $\Phi(x,y,z)$ is defined everywhere except at the charge. The ideas of domain and range carry over unchanged; differentiation and integration of such functions are the subject of [Ch. 6](../part-2-calculus/06-multivariable-calculus.md).

> [!warning] Common Pitfall — the formula is not the domain
> Simplifying a formula can silently enlarge its domain. $\dfrac{x^2-1}{x-1}$ and $x+1$ agree except at $x=1$, where the first is undefined; $\left(\sqrt x\right)^2$ and $x$ agree only for $x\ge0$; $e^{\ln x}$ and $x$ agree only for $x>0$. State the domain before simplifying, and keep it afterwards.

> [!tip] Engineering Note — domains as operating regions
> A device model is a function with a restricted domain. The small-signal model of a transistor is valid only near its bias point; the Shockley equation fails at reverse breakdown; an op-amp's linear model $v_o=A(v_+-v_-)$ holds only while $|v_o|$ is below the rails. Many design failures are evaluations of a function outside its domain.

## 1.9 Inverse and composite functions

### Why composition and inversion

Engineering systems are built by connecting blocks in series: a sensor converts temperature to resistance, a bridge converts resistance to voltage, an amplifier scales the voltage, an ADC converts it to a number. The overall input–output relation is a **composite function**. Measurement is the reverse problem: from the number read out, recover the temperature. That requires the **inverse** of each block — and it is possible only if each block is one-to-one on the range of interest. Composition also underlies the chain rule of differentiation ([Ch. 4](../part-2-calculus/04-differential-calculus.md)) and substitution in integration ([Ch. 5](../part-2-calculus/05-integral-calculus.md)), both of which require recognising a function as "an outer function of an inner function".

### Composite functions

> [!definition] Composite function
> The **composition** of $f$ with $g$ is
> $$
> (f\circ g)(x)=f\big(g(x)\big),
> $$
> read "$f$ after $g$". Its domain is the set of $x$ in the domain of $g$ for which $g(x)$ lies in the domain of $f$.

Composition is **associative**, $(f\circ g)\circ h=f\circ(g\circ h)$, but in general **not commutative**: $f\circ g\neq g\circ f$. A signal chain $x\to g\to f\to y$ is the composite $y=(f\circ g)(x)$ — note that the block applied first is written on the right.

> [!example] Worked Example 1.47 — Composition and its domain
> Let $f(x)=\sqrt{x}$ and $g(x)=1-x^2$.
>
> **(a)** $(f\circ g)(x)=\sqrt{1-x^2}$. Need $g(x)\ge0$: $1-x^2\ge0$, so the domain is $[-1,1]$; the range is $[0,1]$ (the upper unit semicircle).
>
> **(b)** $(g\circ f)(x)=1-(\sqrt x)^2=1-x$, but only for $x$ in the domain of $f$: the domain is $[0,\infty)$, even though the simplified formula $1-x$ looks defined everywhere.
>
> **(c)** $f\circ g\ne g\circ f$: at $x=0.5$, $(f\circ g)(0.5)=\sqrt{0.75}=0.866$ while $(g\circ f)(0.5)=0.5$.

> [!example] Worked Example 1.48 — Decomposing a function
> Write $h(t)=e^{-(t-2)^2}$ and $k(x)=\ln\big(1+\cos^2x\big)$ as compositions of simpler functions.
>
> $h=f_1\circ f_2\circ f_3$ with $f_3(t)=t-2$ (delay), $f_2(u)=-u^2$, $f_1(w)=e^w$.
>
> $k=F_1\circ F_2\circ F_3$ with $F_3(x)=\cos x$, $F_2(u)=1+u^2$, $F_1(w)=\ln w$. The domain is all of $\mathbb{R}$, since $1+\cos^2x\ge1>0$.
>
> Recognising this nesting — innermost to outermost — is exactly the skill needed to apply the chain rule.

### One-to-one functions

> [!definition] One-to-one (injective), onto (surjective)
> $f:A\to B$ is **one-to-one** if $f(x_1)=f(x_2)\Rightarrow x_1=x_2$ — distinct inputs give distinct outputs. Graphically, every horizontal line meets the graph at most once (the **horizontal-line test**). $f$ is **onto** $B$ if its range is all of $B$. A function that is both is a **bijection**.

> [!theorem] Strict monotonicity implies one-to-one
> A function that is strictly increasing (or strictly decreasing) on an interval is one-to-one there.

*Proof.* If $x_1\ne x_2$, say $x_1<x_2$, then strict increase gives $f(x_1)<f(x_2)$, so $f(x_1)\ne f(x_2)$. $\square$

So $x^3$, $e^x$, $\ln x$ and $\tanh x$ are one-to-one on their whole domains, while $x^2$ and $\sin x$ are not — but become so when restricted to $[0,\infty)$ and $[-\pi/2,\pi/2]$ respectively.

### Inverse functions

> [!definition] Inverse function
> If $f$ is one-to-one with domain $A$ and range $B$, its **inverse** $f^{-1}:B\to A$ is defined by
> $$
> f^{-1}(y)=x\iff f(x)=y.
> $$
> Equivalently, $f^{-1}\big(f(x)\big)=x$ for all $x\in A$ and $f\big(f^{-1}(y)\big)=y$ for all $y\in B$. The domain of $f^{-1}$ is the range of $f$, and vice versa.

The graph of $f^{-1}$ is the reflection of the graph of $f$ in the line $y=x$, because the point $(a,b)$ lies on the graph of $f$ exactly when $(b,a)$ lies on the graph of $f^{-1}$. If $f$ is increasing, so is $f^{-1}$.

**Procedure.** To find $f^{-1}$: write $y=f(x)$; solve for $x$ in terms of $y$ (if the solution is not unique, $f$ is not one-to-one — restrict its domain); state the domain of $f^{-1}$ as the range of $f$; and, if desired, rename the variable.

> [!example] Worked Example 1.49 — A fractional-linear inverse
> $f(x)=\dfrac{2x+1}{x-3}$, $x\ne3$.
>
> **Step 1.** $y(x-3)=2x+1\Rightarrow xy-3y=2x+1\Rightarrow x(y-2)=3y+1$.
>
> **Step 2.** For $y\ne2$,
> $$
> f^{-1}(y)=\frac{3y+1}{y-2}.
> $$
> **Step 3 — domain/range.** The step fails only at $y=2$, so the range of $f$ is $\mathbb{R}\setminus\{2\}$ — the horizontal asymptote of $f$, which it never reaches.
>
> **Check:** $f(0)=-\tfrac13$ and $f^{-1}\!\left(-\tfrac13\right)=\dfrac{-1+1}{-\tfrac73}=0$ ✓.

> [!example] Worked Example 1.50 — Inverting the logistic and decibel functions
> **(a)** The logistic (sigmoid) function $\sigma(x)=\dfrac{1}{1+e^{-x}}$ maps $\mathbb{R}$ onto $(0,1)$ and is strictly increasing. Solve $y=\dfrac1{1+e^{-x}}$: $1+e^{-x}=\dfrac1y$, $e^{-x}=\dfrac{1-y}{y}$, so
> $$
> \sigma^{-1}(y)=\ln\frac{y}{1-y},\qquad 0<y<1.
> $$
> **Check:** $\sigma^{-1}(0.8)=\ln4=1.3863$ and $\sigma(1.3863)=1/(1+0.25)=0.8$ ✓.
>
> **(b)** Voltage gain in dB is $G=20\log_{10}A$ for $A>0$. Its inverse is $A=10^{G/20}$; a gain of $26$ dB is a voltage ratio $10^{1.3}=19.95$.

> [!example] Worked Example 1.51 — Inverting a measurement chain
> A platinum RTD has $R(T)=R_0(1+\alpha T)$ with $R_0=100\ \Omega$, $\alpha=0.00385\ {}^\circ\text{C}^{-1}$. It forms a divider with a fixed $100\ \Omega$ resistor across $5\ \text{V}$, giving $V=\dfrac{5R}{R+100}$, read by a 10-bit ADC as the integer $N=\left\lfloor 1024\,V/5\right\rfloor$.
>
> **(a) Forward (composite).** At $T=50.39\ ^\circ\text{C}$: $R=119.4\ \Omega$, $V=5\times119.4/219.4=2.7211\ \text{V}$, $N=\lfloor557.27\rfloor=557$.
>
> **(b) Inverse.** Each block is strictly increasing on its operating range, so invert them in reverse order:
> $$
> V\approx\frac{5N}{1024},\qquad R=\frac{100\,V}{5-V},\qquad T=\frac{R/R_0-1}{\alpha}.
> $$
> From $N=557$: $V=2.7197\ \text{V}$, $R=119.27\ \Omega$, $T=50.06\ ^\circ\text{C}$.
>
> **Interpretation.** The recovered $50.06\ ^\circ\text{C}$ differs from the true $50.39\ ^\circ\text{C}$ because the floor function in the ADC is **not** one-to-one: every voltage in a $4.9\ \text{mV}$ band maps to the same code. Its "inverse" can only return a representative value — here a $0.33\ ^\circ\text{C}$ quantisation error. The general rule: $(f\circ g)^{-1}=g^{-1}\circ f^{-1}$ — undo the last operation first.

> [!theorem] Inverse of a composite
> If $f$ and $g$ are one-to-one and $f\circ g$ is defined, then
> $$
> (f\circ g)^{-1}=g^{-1}\circ f^{-1}.
> $$

(To undo "put on socks, then shoes", take off the shoes first.)

### Inverse trigonometric functions

$\sin$, $\cos$ and $\tan$ are periodic, hence far from one-to-one. Their inverses are defined by restricting each to a **principal** interval on which it is strictly monotonic and takes every value once.

| Function | Domain | Range (principal values) | Restriction inverted |
|---|---|---|---|
| $\arcsin x$ | $[-1,1]$ | $[-\pi/2,\pi/2]$ | $\sin$ on $[-\pi/2,\pi/2]$ |
| $\arccos x$ | $[-1,1]$ | $[0,\pi]$ | $\cos$ on $[0,\pi]$ |
| $\arctan x$ | $\mathbb{R}$ | $(-\pi/2,\pi/2)$ | $\tan$ on $(-\pi/2,\pi/2)$ |

Useful identities, valid on the stated domains:

$$
\arcsin x+\arccos x=\frac\pi2,\qquad \arctan x+\arctan\frac1x=\frac\pi2\ (x>0),\qquad \arcsin(-x)=-\arcsin x,\qquad \arccos(-x)=\pi-\arccos x.
$$

> [!example] Worked Example 1.52 — Compositions with inverse trig functions
> **(a)** Simplify $\cos(\arctan x)$. Let $\theta=\arctan x\in(-\pi/2,\pi/2)$, so $\tan\theta=x$ and $\cos\theta>0$. From $1+\tan^2\theta=\sec^2\theta$,
> $$
> \cos(\arctan x)=\frac{1}{\sqrt{1+x^2}}.
> $$
> (A right triangle with opposite side $x$ and adjacent side $1$ gives the same result.) **Check** at $x=2$: $\cos(1.1071)=0.4472=1/\sqrt5$ ✓.
>
> **(b)** Similarly $\sin(\arccos x)=\sqrt{1-x^2}$, positive because $\arccos x\in[0,\pi]$.
>
> **(c)** $\arcsin(\sin2.5)$ is **not** $2.5$, since $2.5\notin[-\pi/2,\pi/2]$. Using $\sin(\pi-\theta)=\sin\theta$, $\sin2.5=\sin(\pi-2.5)$ with $\pi-2.5=0.6416$ in the principal range, so $\arcsin(\sin2.5)=0.6416$.

```desmos-graph
left=-2; right=2; bottom=-2; top=2
---
y=\sin(x)|blue|-1.5708<x<1.5708
y=\arcsin(x)|red
y=x|grey|dashed
```
*$\sin x$ restricted to $[-\pi/2,\pi/2]$ (blue) and its inverse $\arcsin x$ (red) are reflections of each other in $y=x$.*

### Inverse hyperbolic functions

$\sinh$ is strictly increasing on $\mathbb{R}$ and $\tanh$ on $\mathbb{R}$ (onto $(-1,1)$); $\cosh$ is one-to-one on $[0,\infty)$. Since they are built from exponentials, their inverses are logarithms.

> [!theorem] Inverse hyperbolic functions
> $$
> \operatorname{arsinh}x=\ln\left(x+\sqrt{x^2+1}\right),\qquad
> \operatorname{arcosh}x=\ln\left(x+\sqrt{x^2-1}\right)\ (x\ge1),\qquad
> \operatorname{artanh}x=\frac12\ln\frac{1+x}{1-x}\ (|x|<1).
> $$

*Derivation of $\operatorname{arsinh}$.* Set $x=\sinh y=\frac12(e^y-e^{-y})$ and multiply by $2e^y$: $e^{2y}-2xe^y-1=0$, a quadratic in $u=e^y$. Its roots are $u=x\pm\sqrt{x^2+1}$; only the $+$ root is positive, so $e^y=x+\sqrt{x^2+1}$ and $y=\ln\left(x+\sqrt{x^2+1}\right)$. $\square$

**Check:** $\operatorname{arsinh}2=\ln(2+\sqrt5)=1.4436$ and $\sinh(1.4436)=2.000$ ✓. Inverse hyperbolic functions appear in the capacitance of two parallel wires, $C'=\pi\varepsilon/\operatorname{arcosh}(d/2a)$, and in the design of Chebyshev filters.

> [!warning] Common Pitfall — the $-1$ superscript
> $\sin^{-1}x$ means $\arcsin x$, **not** $1/\sin x=\csc x$; yet $\sin^2x$ means $(\sin x)^2$. Likewise $f^{-1}(x)$ is the inverse function, which is in general **not** $1/f(x)$. Prefer the unambiguous names $\arcsin$, $\arccos$, $\arctan$. Also, $f^{-1}(f(x))=x$ holds only on the domain on which $f$ was inverted.

> [!tip] Engineering Note — calibration and linearisation
> Inverting a sensor characteristic is called **calibration** or **linearisation**. When the forward model has no closed-form inverse (a thermocouple polynomial, the full diode equation with series resistance), the inverse is computed numerically — by a lookup table with interpolation, or by Newton's method ([Ch. 22](../part-7-numerical-computational/22-numerical-methods.md)). Either way the forward function must be one-to-one on the measuring range, or the reading is ambiguous.

## 1.10 Mathematical notation and dimensional consistency

### Why notation and dimensions

Mathematics is a language with a compact, precise grammar. Datasheets, standards, textbooks and papers state results as "for all $\varepsilon>0$ there exists $\delta>0$ such that…", or "$P\Rightarrow Q$", or with nested sums; misreading a quantifier or the direction of an implication changes the meaning entirely. At the same time, every equation an engineer writes describes physical quantities, and physical quantities carry **dimensions**. Requiring both sides of every equation to have the same dimensions is the cheapest, fastest error check there is — it catches most algebra slips in seconds, before any numbers are substituted.

### Logic and set notation

| Symbol | Meaning | Symbol | Meaning |
|---|---|---|---|
| $\forall$ | for all | $\exists$ | there exists |
| $\exists!$ | there exists exactly one | $\neg P$ | not $P$ |
| $P\Rightarrow Q$ | $P$ implies $Q$ | $P\iff Q$ | $P$ if and only if $Q$ |
| $\land,\ \lor$ | and, or (inclusive) | $:=$ or $\equiv$ | is defined as |
| $\in,\ \notin$ | is / is not an element of | $\subseteq,\ \subset$ | subset of, proper subset of |
| $\cup,\ \cap$ | union, intersection | $A\setminus B$ | elements of $A$ not in $B$ |
| $\varnothing$ | empty set | $A\times B$ | Cartesian product (ordered pairs) |
| $\approx$ | approximately equal | $\propto$ | proportional to |
| $\sim$ | asymptotically equal / of the order of | $\ll,\ \gg$ | much less, much greater than |

**Set-builder notation** describes a set by a property: $\{x\in\mathbb{R}: x^2<4\}=(-2,2)$, read "the set of real $x$ such that $x^2<4$". **Cartesian products** build multi-dimensional spaces: $\mathbb{R}^2=\mathbb{R}\times\mathbb{R}$ is the plane, $\mathbb{R}^3$ is space.

**Implication.** $P\Rightarrow Q$ says that whenever $P$ is true, $Q$ is true. Then $P$ is **sufficient** for $Q$ and $Q$ is **necessary** for $P$. The implication is equivalent to its **contrapositive** $\neg Q\Rightarrow\neg P$, but **not** to its **converse** $Q\Rightarrow P$. For example, "if a system is linear and time-invariant, its response to a sinusoid is a sinusoid at the same frequency" is true; the converse is false. "$P\iff Q$" asserts both directions.

**Quantifiers.** The order matters. "$\forall x\ \exists y: y>x$" (every number has a larger one) is true; "$\exists y\ \forall x: y>x$" (some number is larger than every number) is false. To negate a quantified statement, swap $\forall$ and $\exists$ and negate the inner statement: $\neg(\forall x: P(x))\iff\exists x:\neg P(x)$. So disproving "all transistors in this batch meet spec" needs only one counter-example.

### Proof techniques

Four patterns cover nearly every argument in this book:

- **Direct proof:** assume $P$, deduce $Q$ (the proof of the triangle inequality in Section 1.1).
- **Contrapositive:** assume $\neg Q$, deduce $\neg P$.
- **Contradiction:** assume the statement false and derive an absurdity (the irrationality of $\sqrt2$).
- **Induction:** to prove $P(n)$ for all $n\in\mathbb{N}$, prove $P(1)$ (base case) and $P(n)\Rightarrow P(n+1)$ (inductive step).

> [!example] Worked Example 1.53 — Proof by induction
> Prove $\displaystyle\sum_{k=1}^{n}k=\frac{n(n+1)}{2}$ for all $n\in\mathbb{N}$.
>
> **Base case.** $n=1$: LHS $=1$, RHS $=\tfrac{1\cdot2}{2}=1$ ✓.
>
> **Inductive step.** Assume the formula for $n$. Then
> $$
> \sum_{k=1}^{n+1}k=\frac{n(n+1)}{2}+(n+1)=\frac{(n+1)(n+2)}{2},
> $$
> which is the formula for $n+1$. By induction it holds for all $n$. $\square$
>
> **Check** at $n=100$: $5050$ — the sum Gauss is said to have found as a schoolboy by pairing $1+100$, $2+99$, … .

### Sums and products

$$
\sum_{k=1}^{n}a_k=a_1+a_2+\dots+a_n,\qquad \prod_{k=1}^{n}a_k=a_1a_2\cdots a_n.
$$

The index $k$ is a **dummy variable**: $\sum_{k=1}^na_k=\sum_{i=1}^na_i$. Rules that follow from the field axioms:

$$
\sum_{k}(ca_k+b_k)=c\sum_ka_k+\sum_kb_k,\qquad \sum_{k=1}^{n}a_k=\sum_{m=0}^{n-1}a_{m+1}\ \ (\text{shift } m=k-1),\qquad \sum_{k=1}^n c=nc.
$$

But $\sum a_kb_k\ne\left(\sum a_k\right)\left(\sum b_k\right)$ in general. Useful closed forms:

$$
\sum_{k=1}^{n}k=\frac{n(n+1)}{2},\qquad \sum_{k=1}^{n}k^2=\frac{n(n+1)(2n+1)}{6},\qquad \sum_{k=1}^{n}k^3=\left(\frac{n(n+1)}{2}\right)^2.
$$

(Check: $\sum_{k=1}^{10}k^2=385=\tfrac{10\cdot11\cdot21}{6}$ ✓.)

> [!theorem] Geometric series
> For $r\ne1$,
> $$
> \sum_{k=0}^{n-1}r^k=1+r+r^2+\dots+r^{n-1}=\frac{1-r^n}{1-r},
> $$
> and for $|r|<1$ the infinite series converges:
> $$
> \sum_{k=0}^{\infty}r^k=\frac{1}{1-r}.
> $$

*Proof.* Let $S=1+r+\dots+r^{n-1}$. Then $rS=r+r^2+\dots+r^n$, and subtracting, $S-rS=1-r^n$ — all other terms cancel ("telescope"). Divide by $1-r$. If $|r|<1$ then $r^n\to0$ as $n\to\infty$, giving the infinite sum. $\square$

> [!example] Worked Example 1.54 — Geometric sums in digital systems
> **(a) Binary full scale.** An 8-bit unsigned register holds at most $\sum_{k=0}^{7}2^k=\dfrac{1-2^8}{1-2}=255$.
>
> **(b) Echoes.** A signal is reflected repeatedly so that each echo has half the amplitude of the previous one. The total of the original plus all echoes is $\sum_{k=0}^{\infty}0.5^k=\dfrac{1}{1-0.5}=2$ times the original amplitude; the first ten terms already give $1.998$.
>
> **(c) Feedback.** In a loop with forward gain $A$ and feedback $\beta$, a signal circulating repeatedly sums to $A\sum_{k=0}^{\infty}(-A\beta)^k=\dfrac{A}{1+A\beta}$ when $|A\beta|<1$ — the closed-loop gain of Example 1.9(b) seen as a geometric series.
>
> **Check (b):** $(1-0.5^{10})/0.5=1.998047$ ✓.

The geometric series underlies discrete-time systems, the Z-transform ([Ch. 17](../part-5-transforms/17-z-transforms.md)) and the analysis of convergence of iterative algorithms.

### Order-of-magnitude notation

**Big-O** notation states how fast a quantity grows or shrinks: $f(x)=O(g(x))$ as $x\to0$ means $|f(x)|\le M|g(x)|$ for some constant $M$ and all small enough $x$. Thus $\sin x=x+O(x^3)$ and $(1+x)^n=1+nx+O(x^2)$, which makes precise the size of the error in the approximations of Sections 1.2 and 1.6. As $n\to\infty$, Big-O describes algorithm cost: Gaussian elimination needs $O(n^3)$ operations, the FFT $O(n\log n)$ ([Ch. 24](../part-7-numerical-computational/24-computational-mathematics.md)).

### Dimensions and units

Every physical quantity has a **dimension** built from the seven SI base dimensions: M (mass), L (length), T (time), I (electric current), $\Theta$ (temperature), N (amount of substance), J (luminous intensity). The dimension of a quantity $Q$ is written $[Q]$. Derived dimensions follow from defining equations:

| Quantity | Defining relation | Dimension | SI unit |
|---|---|---|---|
| Charge $Q$ | $Q=It$ | $\mathrm{T\,I}$ | C |
| Force $F$ | $F=ma$ | $\mathrm{M\,L\,T^{-2}}$ | N |
| Energy $W$ | force $\times$ distance | $\mathrm{M\,L^2\,T^{-2}}$ | J |
| Power $P$ | $W/t$ | $\mathrm{M\,L^2\,T^{-3}}$ | W |
| Voltage $V$ | $P/I$ | $\mathrm{M\,L^2\,T^{-3}\,I^{-1}}$ | V |
| Resistance $R$ | $V/I$ | $\mathrm{M\,L^2\,T^{-3}\,I^{-2}}$ | $\Omega$ |
| Conductance $G$ | $1/R$ | $\mathrm{M^{-1}\,L^{-2}\,T^{3}\,I^{2}}$ | S |
| Capacitance $C$ | $Q/V$ | $\mathrm{M^{-1}\,L^{-2}\,T^{4}\,I^{2}}$ | F |
| Inductance $L$ | $V\,t/I$ | $\mathrm{M\,L^2\,T^{-2}\,I^{-2}}$ | H |
| Magnetic flux $\Phi$ | $V\,t$ | $\mathrm{M\,L^2\,T^{-2}\,I^{-1}}$ | Wb |
| Electric field $E$ | $V/\text{length}$ | $\mathrm{M\,L\,T^{-3}\,I^{-1}}$ | V/m |
| Permittivity $\varepsilon$ | $C\cdot$length/area | $\mathrm{M^{-1}\,L^{-3}\,T^{4}\,I^{2}}$ | F/m |
| Permeability $\mu$ | $L\cdot$length/area | $\mathrm{M\,L\,T^{-2}\,I^{-2}}$ | H/m |
| Conductivity $\sigma$ | $G\cdot$length/area | $\mathrm{M^{-1}\,L^{-3}\,T^{3}\,I^{2}}$ | S/m |
| Frequency $f$, $\omega$ | cycles or radians per second | $\mathrm{T^{-1}}$ | Hz, rad/s |

> [!theorem] Rules of dimensional consistency
> 1. Only quantities of the **same dimension** may be added, subtracted, equated or compared.
> 2. Dimensions **multiply and divide** like algebraic symbols, and powers apply to them: $[x^n]=[x]^n$.
> 3. The arguments of $\exp$, $\ln$, $\sin$, $\cos$, $\tanh$ and of any other function defined by a power series must be **dimensionless**, and so are the values of these functions.
> 4. Pure numbers ($2$, $\pi$, $\tfrac12$) and angles in radians are dimensionless.

*Why rule 3 holds.* $e^x=1+x+\tfrac{x^2}{2}+\dots$ adds $1$, $x$ and $x^2$; by rule 1 this is possible only if $x$ is dimensionless. The same argument applies to every function with a power series. $\square$

> [!example] Worked Example 1.55 — Checking time constants and energy
> **(a)** $[RC]=\mathrm{M\,L^2\,T^{-3}\,I^{-2}}\cdot\mathrm{M^{-1}L^{-2}T^{4}I^{2}}=\mathrm{T}$, so $e^{-t/RC}$ has a dimensionless exponent ✓.
>
> **(b)** $[L/R]=\dfrac{\mathrm{M\,L^2\,T^{-2}\,I^{-2}}}{\mathrm{M\,L^2\,T^{-3}\,I^{-2}}}=\mathrm{T}$ ✓.
>
> **(c)** $[LC]=\mathrm{T^2}$, so $\omega_0=1/\sqrt{LC}$ has dimension $\mathrm{T^{-1}}$ — a frequency ✓.
>
> **(d)** Capacitor energy $W=\tfrac12CV^2$: $[C][V]^2=\mathrm{M^{-1}L^{-2}T^4I^2}\cdot\mathrm{M^2L^4T^{-6}I^{-2}}=\mathrm{M\,L^2\,T^{-2}}$, an energy ✓. Numerically $\tfrac12(100\ \mu\text{F})(20\ \text{V})^2=0.02\ \text{J}$.
>
> **(e)** A formula claiming $\tau=R/C$ fails immediately: $[R/C]=\mathrm{M^2L^4T^{-7}I^{-4}}\ne\mathrm{T}$.

> [!example] Worked Example 1.56 — Checking an electromagnetic formula
> **(a) Speed of light.** $[\mu\varepsilon]=\mathrm{M\,L\,T^{-2}I^{-2}}\cdot\mathrm{M^{-1}L^{-3}T^4I^2}=\mathrm{L^{-2}T^{2}}$, so $1/\sqrt{\mu_0\varepsilon_0}$ is a speed. Numerically $1/\sqrt{(8.854\times10^{-12})(4\pi\times10^{-7})}=2.998\times10^8\ \text{m/s}$ ✓.
>
> **(b) Skin depth.** The claimed formula $\delta=\sqrt{2/(\omega\mu\sigma)}$ should be a length. $[\omega\mu\sigma]=\mathrm{T^{-1}}\cdot\mathrm{M\,L\,T^{-2}I^{-2}}\cdot\mathrm{M^{-1}L^{-3}T^3I^2}=\mathrm{L^{-2}}$, so $[\delta]=\mathrm{L}$ ✓. For copper ($\sigma=5.8\times10^7\ \text{S/m}$, $\mu=\mu_0$) at $60$ Hz, $\delta=8.53\ \text{mm}$; at $1$ MHz, $\delta=66\ \mu\text{m}$.

### Dimensional analysis as a design tool

Dimensional consistency can do more than check formulas: it can predict their form. If a quantity depends only on a known list of others, its formula must be a product of powers with matching dimensions.

> [!example] Worked Example 1.57 — Predicting the resonance formula
> Suppose the resonant angular frequency of an ideal $LC$ circuit depends only on $L$ and $C$: $\omega_0=k\,L^aC^b$ with $k$ dimensionless. Find $a$ and $b$.
>
> **Step 1 — dimensions.** $\mathrm{T^{-1}}=\left(\mathrm{M\,L^2T^{-2}I^{-2}}\right)^a\left(\mathrm{M^{-1}L^{-2}T^{4}I^{2}}\right)^b$.
>
> **Step 2 — match exponents.** M: $a-b=0$. L: $2a-2b=0$. I: $-2a+2b=0$. T: $-2a+4b=-1$.
>
> **Step 3 — solve.** $a=b$ and $2a=-1$, so $a=b=-\tfrac12$, giving $\omega_0=k/\sqrt{LC}$.
>
> The analysis cannot find $k$ (circuit analysis gives $k=1$), but it pins down the structure. Its systematic form, the Buckingham $\pi$ theorem, is developed in [Ch. 25](../part-7-numerical-computational/25-mathematical-modeling.md).

### Units, prefixes and conversions

Dimensions say *what kind* of quantity; **units** say *how much*. SI prefixes scale units by powers of $10^3$:

| Prefix | p | n | $\mu$ | m | k | M | G |
|---|---|---|---|---|---|---|---|
| Factor | $10^{-12}$ | $10^{-9}$ | $10^{-6}$ | $10^{-3}$ | $10^{3}$ | $10^{6}$ | $10^{9}$ |

Convert by multiplying by factors equal to one, such as $\frac{10^{-9}\ \text{F}}{1\ \text{nF}}$, so that units cancel algebraically. Prefixes bind to the unit before powers apply: $1\ \text{cm}^2=(10^{-2}\ \text{m})^2=10^{-4}\ \text{m}^2$, not $10^{-2}\ \text{m}^2$. Treat units as algebraic symbols throughout a calculation and the final unit appears automatically — if it is not the expected one, there is an error.

> [!warning] Common Pitfall — dimensionless is not unitless in practice
> Angles are dimensionless, but radians and degrees are different units: $\omega t$ must be in radians for $\sin\omega t$ to be correct with $\omega$ in rad/s. Hertz and rad/s have the same dimension $\mathrm{T^{-1}}$ but differ by $2\pi$ — a dimensional check cannot catch $f$ written for $\omega$. Likewise decibels are dimensionless but need a stated reference (Section 1.5).

> [!tip] Engineering Note
> Dimensional checking cannot detect missing dimensionless factors ($2\pi$, $\tfrac12$, $\sqrt2$), but it catches the majority of algebra slips in seconds. Make it a habit: check the dimensions of every derived formula before substituting numbers, check limiting cases (what happens as $R\to0$ or $\omega\to\infty$?), and estimate the order of magnitude of the answer before trusting a calculator.

## Summary

| Topic | Key result |
|---|---|
| Number sets | $\mathbb{N}\subset\mathbb{Z}\subset\mathbb{Q}\subset\mathbb{R}\subset\mathbb{C}$; $\mathbb{C}$ is not ordered |
| Triangle inequality | $\lvert x+y\rvert\le\lvert x\rvert+\lvert y\rvert$, $\ \big\lvert\lvert x\rvert-\lvert y\rvert\big\rvert\le\lvert x-y\rvert$ |
| Machine epsilon (double) | $\varepsilon_{\text{mach}}=2^{-52}\approx2.2\times10^{-16}$ |
| Binomial theorem | $(a+b)^n=\sum_{k=0}^n\binom nk a^{n-k}b^k$, $\ \binom nk=\dfrac{n!}{k!(n-k)!}$ |
| Binomial approximation | $(1+x)^n\approx1+nx,\ \lvert x\rvert\ll1$, any real $n$ |
| Completing the square | $ax^2+bx+c=a\left(x+\frac{b}{2a}\right)^2+c-\frac{b^2}{4a}$ |
| Quadratic formula | $x=\dfrac{-b\pm\sqrt{b^2-4ac}}{2a}$; stable form $q=-\tfrac12\big(b+\operatorname{sgn}(b)\sqrt\Delta\big)$, $x_1=q/a$, $x_2=c/q$ |
| Vieta (quadratic) | $x_1+x_2=-b/a$, $\ x_1x_2=c/a$ |
| Cramer ($2\times2$) | $x=\dfrac{c_1b_2-c_2b_1}{a_1b_2-a_2b_1}$, $\ y=\dfrac{a_1c_2-a_2c_1}{a_1b_2-a_2b_1}$ |
| AM–GM | $\tfrac12(a+b)\ge\sqrt{ab}$, equality iff $a=b$ |
| Remainder / factor theorem | $p(x)\bmod(x-r)=p(r)$; $\ p(r)=0\iff(x-r)\mid p(x)$ |
| Fundamental theorem of algebra | $p(x)=a_n\prod_{i=1}^n(x-r_i)$; real coefficients $\Rightarrow$ conjugate pairs |
| Cover-up rule | $A=\Big[(x-a)R(x)\Big]_{x=a}$ |
| Repeated pole | $A_{m-k}=\dfrac1{k!}\dfrac{d^k}{dx^k}\Big[(x-a)^mR(x)\Big]_{x=a}$ |
| Euler's number | $e=\lim_{n\to\infty}(1+1/n)^n=\sum 1/k!=2.71828\ldots$ |
| Logarithms | $\log(xy)=\log x+\log y$, $\ \log x^p=p\log x$, $\ \log_a x=\ln x/\ln a$ |
| First-order transient | $e^{-t/\tau}$; $t_{10\%}=\tau\ln10$; 10–90 % rise time $=\tau\ln9=2.197\tau$ |
| Decibels | $10\log_{10}(P_2/P_1)=20\log_{10}(V_2/V_1)$; $0\ \text{dBm}=1\ \text{mW}$ |
| Hyperbolic | $\cosh x=\tfrac12(e^x+e^{-x})$, $\sinh x=\tfrac12(e^x-e^{-x})$, $\cosh^2x-\sinh^2x=1$ |
| Pythagorean identity | $\sin^2\theta+\cos^2\theta=1$, $\ 1+\tan^2\theta=\sec^2\theta$ |
| Compound angle | $\cos(A\pm B)=\cos A\cos B\mp\sin A\sin B$, $\ \sin(A\pm B)=\sin A\cos B\pm\cos A\sin B$ |
| Double angle / power reduction | $\sin2\theta=2\sin\theta\cos\theta$, $\ \cos^2\theta=\tfrac12(1+\cos2\theta)$ |
| Product-to-sum | $\cos A\cos B=\tfrac12[\cos(A-B)+\cos(A+B)]$ |
| Harmonic addition | $a\cos\theta+b\sin\theta=\sqrt{a^2+b^2}\cos\big(\theta-\operatorname{atan2}(b,a)\big)$ |
| Law of cosines | $c^2=a^2+b^2-2ab\cos C$ |
| Polar | $r=\sqrt{x^2+y^2}$, $\ \theta=\operatorname{atan2}(y,x)$ |
| Cylindrical $(r,\phi,z)$ | $x=r\cos\phi,\ y=r\sin\phi$; $\ dV=r\,dr\,d\phi\,dz$ |
| Spherical $(R,\theta,\phi)$ | $x=R\sin\theta\cos\phi,\ y=R\sin\theta\sin\phi,\ z=R\cos\theta$; $\ dV=R^2\sin\theta\,dR\,d\theta\,d\phi$ |
| Even–odd split | $f_e=\tfrac12[f(x)+f(-x)]$, $\ f_o=\tfrac12[f(x)-f(-x)]$ |
| Graph transformation | $A\,f\big(b(x-c)\big)+d$: shift right $c$, scale by $1/b$ horizontally, $A$ vertically, up $d$ |
| Inverse of composite | $(f\circ g)^{-1}=g^{-1}\circ f^{-1}$ |
| Inverse hyperbolic | $\operatorname{arsinh}x=\ln\big(x+\sqrt{x^2+1}\big)$, $\ \operatorname{artanh}x=\tfrac12\ln\frac{1+x}{1-x}$ |
| Geometric series | $\sum_{k=0}^{n-1}r^k=\dfrac{1-r^n}{1-r}$; $\ \sum_{k=0}^\infty r^k=\dfrac1{1-r}$, $\lvert r\rvert<1$ |
| Dimensions | $[V]=\mathrm{M\,L^2T^{-3}I^{-1}}$; $[RC]=[L/R]=\mathrm{T}$; $[LC]=\mathrm{T^2}$; arguments of $\exp,\ln,\sin$ dimensionless |

## Further reading

- E. Kreyszig, *Advanced Engineering Mathematics*, 10th ed., Wiley, 2011 — Appendix 3 (formulas and identities) for reference.
- K. A. Stroud and D. J. Booth, *Engineering Mathematics*, 8th ed., Palgrave Macmillan, 2020 — Foundation programmes F.1–F.11 and Programmes 1–3 give extensive practice in algebra, logarithms, trigonometry and functions.
- G. James et al., *Modern Engineering Mathematics*, 6th ed., Pearson, 2020 — Chapters 1–2 (number, algebra, functions, coordinate geometry).
- G. Strang, *Calculus*, 3rd ed., Wellesley-Cambridge Press, 2017 — Chapter 1 on functions and their graphs (freely available from MIT OpenCourseWare).
